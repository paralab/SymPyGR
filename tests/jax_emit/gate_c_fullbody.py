#!/usr/bin/env python
"""Gate C' -- full-body differential.

    python tests/jax_emit/gate_c_fullbody.py [config.py ...]

Checks every emitted statement of a real solver's RHS against the sympy
expression it was printed from. Gate B does this on a 14-expression battery;
this runs the actual population (BSSN: 537 temps + 24 outputs).

The two sides share only leaf VALUES -- re-parsing the emitted text with sympy
would just check it against itself. CPU-only; numpy stands in for jnp.

No byte-diff against DendroJAX's hand-written RHS: it descends from
bssneqs_SSL_HD_dxsq.cpp, which predates the current generator (537 vs 826
temps, different naming and addressing). It fixed the transliteration rules
run_all.py encodes, but is not a regeneration target.
"""

import argparse
import ast
import importlib.util
import math
import os
import random
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import numpy as np  # noqa: E402
import sympy as sym  # noqa: E402

import dendrosym  # noqa: E402
import dendrosym.derivs  # noqa: E402
from dendrosym.codegen_jax import atomize_derivs  # noqa: E402
from dendrosym.jax_printer import jax_symbol_name  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCRATCH = os.environ.get("JAX_GATE_SCRATCH", "/tmp/dendro_jax_gate")


def load_config(path):
    """Import a config, tolerating its scratch writes.

    bssn_eqns.py has no dendrosym.run entry point and writes fragments to cwd,
    so import from a scratch dir.
    """
    path = os.path.abspath(path)          # cwd moves below; resolve first
    os.makedirs(SCRATCH, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(SCRATCH)
    try:
        spec = importlib.util.spec_from_file_location("gate_cfg", path)
        mod = importlib.util.module_from_spec(spec)
        saved = sys.argv
        sys.argv = [os.path.basename(path)]
        try:
            spec.loader.exec_module(mod)
        except SystemExit:
            pass
        finally:
            sys.argv = saved
    finally:
        os.chdir(cwd)
    return mod


def find_config_object(mod):
    for name in dir(mod):
        obj = getattr(mod, name)
        if isinstance(obj, dendrosym.general_configs.DendroConfiguration):
            return obj
    raise RuntimeError("no DendroConfigs instance found in the config module")


def source_names(src):
    """Bare identifiers the emitted source reads, jnp excluded."""
    names = set()
    for node in ast.walk(ast.parse(src, mode="eval")):
        if isinstance(node, ast.Name):
            names.add(node.id)
    names.discard("jnp")
    return names


# BSSN's packed 6-component metric names its diagonal 0/3/5; the two-index
# forms (gt00, gt_mat00, gammat11) are detected structurally below.
PACKED_DIAGONAL = ("gt0", "gt3", "gt5", "At0", "At3", "At5")
UNIT_FIELDS = ("alpha", "chi", "psi", "phi", "W")
_TWO_INDEX = re.compile(r"^([A-Za-z_]+?)(\d)(\d)$")


def leaf_value(name, rng):
    """Physically-shaped random leaf: perturb around flat space.

    Uniform values make the metric determinant negative, the first sqrt NaNs,
    and the NaN eats every downstream statement -- gutting coverage rather than
    failing. Diagonals detected structurally, not by name: BSSN packs gt0..gt5,
    CCZ4 writes gt00..gt22.
    """
    if name.startswith(("grad", "agrad", "kograd", "d2", "mixed")):
        return rng.uniform(-0.02, 0.02)                 # small gradients
    base = name.split("[")[0]
    m = _TWO_INDEX.match(base)
    if (m and m.group(2) == m.group(3)) or base in PACKED_DIAGONAL:
        return 1.0 + rng.uniform(-0.05, 0.05)           # tensor diagonal
    if base in UNIT_FIELDS:
        return 1.0 + rng.uniform(-0.05, 0.05)           # lapse / conformal factor
    v = rng.uniform(0.02, 0.05)                         # bounded away from 0
    return v if rng.random() < 0.5 else -v


def check_var_type(cfg, vt, trials, tol, seed):
    failures = []
    cfg.find_derivatives(vt)        # else this checks a body the emitter never emits
    body = cfg.generate_rhs_code(vt, arc_type="jax")
    print(f"  {vt:12s} {len(body.statements)} statements "
          f"({body.n_temps} temps, {len(body.outputs)} outputs)")

    # ---- structural: valid, parseable Python; no numpy leak ----
    for lhs, rhs in body.statements:
        if not lhs.isidentifier():
            failures.append(f"{vt}:lhs:{lhs}")
            print(f"      LHS not an identifier: {lhs!r}")
            break
        try:
            ast.parse(f"{lhs} = {rhs}")
        except SyntaxError as exc:
            failures.append(f"{vt}:syntax:{lhs}")
            print(f"      SYNTAX {lhs}: {exc}")
            break
    leaked = [lhs for lhs, rhs in body.statements if "numpy." in rhs]
    if leaked:
        failures.append(f"{vt}:numpy")
        print(f"      NUMPY LEAK in {len(leaked)}, e.g. {leaked[:3]}")

    # ---- atomize the source expressions to the emitter's own leaf names ----
    atom = [atomize_derivs(e) for e in body.exprs]

    # Leaves come from the atomized sympy side (authoritative: the source was
    # printed from it). `lf[1]` is one symbol there but an array subscript in
    # the source; both must see the same value.
    defined, leaves = set(), set()
    for (lhs, _rhs), a in zip(body.statements, atom):
        for s in a.free_symbols:
            n = jax_symbol_name(s.name)
            if n not in defined:
                leaves.add(n)
        defined.add(lhs)

    scalars, arrays = set(), {}
    for n in leaves:
        m = re.fullmatch(r"(\w+)\[(\d+)\]", n)
        if m:
            base, idx = m.group(1), int(m.group(2))
            arrays[base] = max(arrays.get(base, -1), idx)
        else:
            scalars.add(n)

    rng = random.Random(seed)
    worst, worst_at, checked = 0.0, None, 0
    nonfinite = set()
    for _t in range(trials):
        scope, ref_env = {"jnp": np}, {}
        for n in sorted(scalars):
            v = leaf_value(n, rng)
            scope[n] = v
            ref_env[sym.Symbol(n)] = sym.Float(v, 40)
        for base, hi_idx in sorted(arrays.items()):
            vec = [leaf_value(base, rng) for _ in range(hi_idx + 1)]
            scope[base] = np.array(vec, dtype=float)
            for i, v in enumerate(vec):
                ref_env[sym.Symbol(f"{base}[{i}]")] = sym.Float(v, 40)

        for (lhs, rhs), a in zip(body.statements, atom):
            try:
                got = float(eval(rhs, {"__builtins__": {}}, scope))  # noqa: S307
            except Exception as exc:
                failures.append(f"{vt}:eval:{lhs}")
                print(f"      EVAL {lhs}: {type(exc).__name__}: {exc}")
                return failures
            scope[lhs] = got

            # the independent side: original sympy tree, high precision
            sub = {s: ref_env[sym.Symbol(jax_symbol_name(s.name))]
                   for s in a.free_symbols
                   if sym.Symbol(jax_symbol_name(s.name)) in ref_env}
            try:
                ref = float(a.subs(sub).evalf(40))
            except Exception:
                ref = None
            if ref is not None and math.isfinite(ref) and math.isfinite(got):
                checked += 1
                rel = abs(got - ref) / max(abs(ref), 1.0)
                if rel > worst:
                    worst, worst_at = rel, lhs
            ref_env[sym.Symbol(lhs)] = sym.Float(got, 40)

        nonfinite.update(k for k in body.outputs
                         if not math.isfinite(float(scope.get(k, float("nan")))))

    expected = len(body.statements) * trials
    cover = checked / expected if expected else 0.0
    print(f"      differential: {checked}/{expected} statement-evals "
          f"({cover:.1%}), worst rel {worst:.3e}"
          + (f" at {worst_at}" if worst_at else ""))

    # Non-finite is a property of the random point (psi4 carries 1/r), not the
    # emitter. Report it; let coverage be what fails.
    if nonfinite:
        print(f"      note: {len(nonfinite)} output(s) non-finite at the random "
              f"point, e.g. {sorted(nonfinite)[:4]}")
    if checked == 0:
        failures.append(f"{vt}:nothing-checked")
        print("      NOTHING CHECKED -- the gate proved nothing")
    elif cover < 0.95:
        failures.append(f"{vt}:coverage")
        print(f"      COVERAGE TOO LOW ({cover:.1%}) -- NaNs are eating the body")
    if worst > tol:
        failures.append(f"{vt}:tolerance")
        print(f"      EXCEEDS TOLERANCE {tol:.1e}")
    return failures


def run_config(path, trials, tol, seed):
    print(f"\n=== {os.path.basename(path)} ===")
    cfg = find_config_object(load_config(path))
    failures = []
    for vt in [v for v in cfg.all_var_names if v != "general"]:
        try:
            failures += check_var_type(cfg, vt, trials, tol, seed)
        except ValueError as exc:
            # "parameter" is a declared var type with no equations attached; the
            # extractor rejects it. Not an emitter failure.
            if "Cannot extract expressions" in str(exc):
                print(f"  {vt:12s} skipped (no equations)")
                continue
            import traceback
            traceback.print_exc()
            failures.append(f"{os.path.basename(path)}:{vt}:{type(exc).__name__}")
        except Exception as exc:
            import traceback
            traceback.print_exc()
            failures.append(f"{os.path.basename(path)}:{vt}:{type(exc).__name__}")
    return failures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("configs", nargs="*",
                    default=[os.path.join(REPO, "bssn_eqns.py")])
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--seed", type=int, default=20260903)
    args = ap.parse_args()

    failures = []
    for path in args.configs:
        failures += run_config(path, args.trials, args.tol, args.seed)

    print()
    if failures:
        print(f"GATE C' FAILED ({len(failures)}): {failures}")
        return 1
    print("GATE C' PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
