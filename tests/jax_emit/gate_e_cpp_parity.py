#!/usr/bin/env python
"""Gate E -- the emitted JAX body must be the C++ body, statement for statement.

    python tests/jax_emit/gate_e_cpp_parity.py <config.py> [config.py ...]

Both backends share one cse_list, so the bodies must correspond exactly once the
struct/naming rewrites are undone. Catches the two paths seeing different
*expressions*, not just different text.
"""

import math
import os
import random
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import dendrosym  # noqa: E402
import dendrosym.project_generator  # noqa: E402
from dendrosym.codegen_jax import change_deriv_names_jax  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from gate_c_fullbody import find_config_object, load_config  # noqa: E402


def cpp_statements(src):
    """`double X = expr;` -> (lhs, rhs)."""
    out = []
    for line in src.splitlines():
        line = line.strip()
        if not line or line.startswith("//"):
            continue
        if not line.endswith(";") or " = " not in line:
            raise AssertionError(f"unparsed C++ statement: {line[:120]}")
        lhs, rhs = line[:-1].split(" = ", 1)
        out.append((lhs.strip(), rhs.strip()))
    return out


def norm_cpp(lhs, rhs, idx_str):
    """C++ statement -> the form the jax body prints."""
    def strip(s):
        s = s.replace(idx_str, "")
        s = re.sub(r"\b(?:in|out)\.", "", s)
        s = re.sub(r"^double\s+", "", s)
        s = s.replace("lambda[", "lambda_param[")
        return s

    lhs = strip(lhs)
    if not lhs.startswith("DENDRO_"):       # out.alpha is the jax alpha_rhs
        lhs += "_rhs"
    rhs = change_deriv_names_jax(strip(rhs))
    # accepted rewrite: x**(-1/2) is 1/sqrt(x) in jnp, 1 ULP apart at worst
    rhs = re.sub(r"pow\(([^,]+), -1\.0/2\.0\)", r"1/sqrt(\1)", rhs)
    return lhs, rhs


def norm_jax(lhs, rhs):
    return lhs, rhs.replace("jnp.", "")


def compare(cfg, var_type):
    cfg.find_derivatives(var_type)
    cpu = cfg.generate_rhs_code(var_type, arc_type="cpu")
    jax = cfg.generate_rhs_code(var_type, arc_type="jax")

    C = cpp_statements(cpu)
    J = list(jax.statements)
    idx_str = dendrosym.derivs.idx_str

    diffs = []
    if len(C) != len(J):
        diffs.append(f"statement count: cpp {len(C)} vs jax {len(J)}")
    for i, (c, j) in enumerate(zip(C, J)):
        a, b = norm_cpp(*c, idx_str), norm_jax(*j)
        if a != b:
            diffs.append(f"stmt {i}\n    cpp: {a[0]} = {a[1][:220]}\n    jax: {b[0]} = {b[1][:220]}")

    # a derivative of a CSE temp is never a buffer the driver can fill
    orphan = sorted({
        m for _lhs, rhs in J
        for m in re.findall(r"\b(?:a?grad|grad2|kograd)(?:_\d)+_(DENDRO_\w+)", rhs)
    })
    if orphan:
        diffs.append(f"deriv buffers named after CSE temps (no driver can fill these): {orphan}")
    return len(C), diffs


def _ic_maps(cfg):
    """C++ name -> jax name, for the two things initial data renames."""
    from dendrosym.jax_printer import safe_name

    sub = {}
    for plist in cfg.all_vars.get("parameter", {}).values():
        for pvar in plist:
            sub[f"{cfg.project_upper}_{pvar.var_name.upper()}"] = safe_name(pvar.var_name)
    for sy, cpp in (getattr(cfg, "runtime_symbol_map", {}) or {}).items():
        sub[str(cpp)] = str(sy)
    return sub


def _diff_numerically(cpp_src, jax_src, trials=6, tol=1e-13):
    """Worst relative difference over random leaf values, or "" if they agree."""
    # a name followed by "(" is a call, never a leaf -- keying on that beats a
    # hand-kept function list, which missed C99's cbrt and bound it to a float
    both = cpp_src + " " + jax_src
    names = set(re.findall(r"\b[A-Za-z_]\w*\b(?!\s*\()", both))
    unknown = set(re.findall(r"\b([A-Za-z_]\w*)\s*\(", both)) - set(_FUNCS)
    if unknown:
        return f"evaluator has no {sorted(unknown)}"
    rng = random.Random(5)
    worst = 0.0
    for _ in range(trials):
        env = {n: rng.uniform(0.3, 1.7) for n in names}
        try:
            a, b = _numeric(cpp_src, env), _numeric(jax_src, env)
        except Exception as exc:
            return f"could not evaluate ({type(exc).__name__}: {exc})"
        if not (math.isfinite(a) and math.isfinite(b)):
            continue
        worst = max(worst, abs(a - b) / max(abs(a), 1.0))
    return "" if worst <= tol else f"differ by {worst:.3e}"


#: what both printers may emit as a call. Extend when a gate says it cannot.
_FUNCS = {
    "sqrt": math.sqrt, "cbrt": lambda v: math.copysign(abs(v) ** (1 / 3), v),
    "exp": math.exp, "log": math.log, "log10": math.log10,
    "sin": math.sin, "cos": math.cos, "tan": math.tan,
    "asin": math.asin, "acos": math.acos, "atan": math.atan, "atan2": math.atan2,
    "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
    "fabs": abs, "abs": abs, "pow": math.pow, "hypot": math.hypot,
    "erf": math.erf, "fmax": max, "fmin": min, "maximum": max, "minimum": min,
}


def _numeric(src, env):
    """Evaluate a printed expression with `env` bound. C99 and jnp both parse."""
    ns = dict(_FUNCS, M_PI=math.pi)
    ns.update(env)
    return float(eval(src, {"__builtins__": {}}, ns))  # noqa: S307


def compare_initial_data(cfg):
    """The emitted IC bodies must be the C++ IC bodies, expression for expression.

    Compared **numerically**, unlike the RHS. The RHS renames at print time so
    both sides share one expression tree; the C++ IC path must `subs()` params to
    reach `PROJECT_PARAM`, and that reorders sympy's Mul args, so the two print
    the same value in a different order. Values are what matter here.
    """
    from dendrosym.project_generator import _jax_initial_data

    jax_ic = _jax_initial_data(cfg)
    if not jax_ic:
        return 0, []

    # the C++ side of the same entries, rendered through the real generator path
    gen = dendrosym.project_generator.DendroProjectGenerator(cfg)
    cpp_ctx = gen._build_context()
    cpp_by_id = {e.get("id"): e for e in cpp_ctx.get("initial_data_types", [])}

    sub = _ic_maps(cfg)
    rename = re.compile("|".join(re.escape(k) for k in sorted(sub, key=len, reverse=True))) \
        if sub else None

    def norm_cpp_line(line):
        line = line.strip().rstrip(";")
        lhs, rhs = line.split(" = ", 1)
        name = lhs[lhs.index("U_") + 2: lhs.rindex("]")].lower() if "U_" in lhs else lhs
        if rename:
            rhs = rename.sub(lambda m: sub[m.group(0)], rhs)
        rhs = re.sub(r"pow\(([^,]+), -1\.0/2\.0\)", r"1/sqrt(\1)", rhs)
        return name, rhs

    diffs, n = [], 0
    for entry in jax_ic["entries"]:
        cpp = cpp_by_id.get(entry["id"])
        if cpp is None or not cpp.get("code"):
            diffs.append(f"id {entry['id']}: no C++ counterpart to compare")
            continue
        cpp_lines = [norm_cpp_line(l) for l in cpp["code"].splitlines() if l.strip()]
        jax_lines = [(lhs.lower(), rhs.replace("jnp.", ""))
                     for lhs, rhs in entry["lines"]]
        if len(cpp_lines) != len(jax_lines):
            diffs.append(f"id {entry['id']}: cpp {len(cpp_lines)} vs jax {len(jax_lines)} fields")
            continue
        for (cn, cr), (jn, jr) in zip(cpp_lines, jax_lines):
            n += 1
            if cn != jn:
                diffs.append(f"id {entry['id']}: field {cn} vs {jn}")
                continue
            bad = _diff_numerically(cr, jr)
            if bad:
                diffs.append(f"id {entry['id']} field {jn}: {bad}")
    # the same check for the single-dict surfaces (symbolic_initial_data /
    # symbolic_analytical_solution), which a config may use instead of ids
    for key, cpp_key in (("symbolic", "symbolic_init_code"),
                         ("analytical", "symbolic_analytical_code")):
        entry = jax_ic[key]
        cpp_code = cpp_ctx.get(cpp_key, "")
        if not entry:
            continue
        if not cpp_code:
            diffs.append(f"{key}: emitted for jax but not for C++")
            continue
        cpp_lines = [norm_cpp_line(l) for l in cpp_code.splitlines() if l.strip()]
        jax_lines = [(lhs.lower(), rhs.replace("jnp.", "")) for lhs, rhs in entry["lines"]]
        if len(cpp_lines) != len(jax_lines):
            diffs.append(f"{key}: cpp {len(cpp_lines)} vs jax {len(jax_lines)} fields")
            continue
        for (cn, cr), (jn, jr) in zip(cpp_lines, jax_lines):
            n += 1
            if cn != jn:
                diffs.append(f"{key}: field {cn} vs {jn}")
                continue
            bad = _diff_numerically(cr, jr)
            if bad:
                diffs.append(f"{key} field {jn}: {bad}")

    for e in jax_ic["entries"] + [x for x in (jax_ic["symbolic"], jax_ic["analytical"]) if x]:
        if e.get("unknown"):
            diffs.append(f"{e.get('func', 'analytical')}: undeclared symbol(s) {e['unknown']}")
    return n, diffs


def main():
    paths = sys.argv[1:]
    if not paths:
        print("usage: gate_e_cpp_parity.py <config.py> [...]")
        return 2

    failures = 0
    for path in paths:
        cfg = find_config_object(load_config(path))
        name = os.path.basename(path)
        for vt in cfg.all_var_names:
            if cfg.all_rhs_functions.get(vt) is None:
                continue
            n, diffs = compare(cfg, vt)
            if diffs:
                failures += 1
                print(f"FAIL  {name}:{vt}  ({n} statements)")
                for d in diffs[:10]:
                    print("      " + d.replace("\n", "\n      "))
                if len(diffs) > 10:
                    print(f"      ... {len(diffs) - 10} more")
            else:
                print(f"PASS  {name}:{vt}  {n} statements identical")

        n, diffs = compare_initial_data(cfg)
        if diffs:
            failures += 1
            print(f"FAIL  {name}:initial_data")
            for d in diffs[:8]:
                print("      " + d.replace("\n", "\n      "))
        elif n:
            print(f"PASS  {name}:initial_data  {n} expressions identical")

    print("\nall gate E checks passed" if not failures
          else f"\n{failures} var_type(s) diverge from the C++ backend")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
