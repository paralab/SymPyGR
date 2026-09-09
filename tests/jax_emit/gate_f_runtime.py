#!/usr/bin/env python
"""Gate F -- the emitted package must import, jit and run under real JAX.

    python tests/jax_emit/gate_f_runtime.py <emitted_pkg_dir>

Gates B-E all reason about the source. This one executes it: real jax, float64,
random near-Minkowski state, every declared kernel jitted and evaluated. Needs
jax + dendrojax installed; skips (exit 0) if they are not.
"""

import importlib
import pathlib
import sys

import numpy as np


def main(pkg_dir):
    pkg_dir = pathlib.Path(pkg_dir).resolve()
    if not (pkg_dir / "__init__.py").exists():
        cands = [p for p in pkg_dir.iterdir()
                 if p.is_dir() and (p / "__init__.py").exists()]
        if len(cands) != 1:
            print(f"no single package under {pkg_dir}: {cands}")
            return 2
        pkg_dir = cands[0]

    sys.path.insert(0, str(pkg_dir.parent))
    try:
        import dendrojax  # noqa: F401  must precede any jax op (x64)
        import jax
        import jax.numpy as jnp
    except ImportError as exc:
        print(f"SKIP  jax/dendrojax not importable ({exc})")
        return 0

    m = importlib.import_module(pkg_dir.name)
    print(f"imported {m.__name__}, dtype {jnp.zeros(1).dtype}")
    assert jnp.zeros(1).dtype == np.float64, "x64 is off"

    rng = np.random.default_rng(7)
    shape = (8, 8, 8)
    flat = {"alpha": 1.0, "chi": 1.0, "gt00": 1.0, "gt11": 1.0, "gt22": 1.0,
            "u": 1.0}

    ok = True
    # keyed on _OUTPUTS: `<name>_rhs` also matches the rhs submodule
    vts = [n[: -len("_OUTPUTS")].lower() for n in dir(m) if n.endswith("_OUTPUTS")]
    assert vts, "the package declares no <VAR_TYPE>_OUTPUTS"

    params = next(getattr(m, n) for n in dir(m) if n.endswith("Params"))()
    for vt in sorted(vts):
        up = vt.upper()
        FIELDS = getattr(m, f"{up}_FIELDS")
        DERIVS = getattr(m, f"{up}_DERIVS")
        EXTRA = getattr(m, f"{up}_EXTRA_INPUTS")
        OUTS = getattr(m, f"{up}_OUTPUTS")

        u = {n: jnp.asarray(flat.get(n, 0.0) + 0.01 * rng.standard_normal(shape))
             for n in FIELDS}
        u.update({n: jnp.asarray(1.0 + rng.random(shape)) for n in EXTRA})
        d = {n: jnp.asarray(0.001 * rng.standard_normal(shape)) for n in DERIVS}

        fn = jax.jit(getattr(m, f"{vt}_rhs"))
        out = [np.asarray(o) for o in fn(u, d, params)]
        bad = [n for n, o in zip(OUTS, out) if not np.isfinite(o).all()]
        shapes = all(o.shape == shape for o in out)

        n0 = fn._cache_size()
        fn(u, d, params)                      # identical shapes must not retrace
        retraced = fn._cache_size() != n0

        # the BC table must be fillable from what the module told the driver
        BCS = getattr(m, f"{up}_BCS", ())
        bc_ok = ({g for r in BCS for g in r[2:5]} <= set(DERIVS)
                 and {r[0] for r in BCS} <= set(OUTS)
                 and {r[1] for r in BCS} <= set(FIELDS))

        good = (not bad and shapes and len(out) == len(OUTS)
                and not retraced and bc_ok)
        print(f"{'PASS' if good else 'FAIL'}  {vt:11s} {len(out)}/{len(OUTS)} outputs, "
              f"shapes {'ok' if shapes else 'BAD'}, "
              f"{'retraced' if retraced else 'no retrace'}, "
              f"{len(BCS)} bc rows {'ok' if bc_ok else 'UNFILLABLE'}"
              + (f", nonfinite {bad}" if bad else ""))
        ok &= good

    # initial data: every entry must produce exactly the evolved fields, shaped
    if hasattr(m, "INITIAL_DATA"):
        main = "evolution" if "evolution" in vts else sorted(vts)[0]
        want = set(getattr(m, f"{main.upper()}_FIELDS"))
        ax = jnp.linspace(-2.0, 2.0, shape[0])
        Z, Y, X = jnp.meshgrid(ax, ax, ax, indexing="ij")
        rt = {n: 0.5 for grp in m.INITIAL_DATA.values() for n in grp[2]}
        for id_type, (nm, _fn, needs) in sorted(m.INITIAL_DATA.items()):
            out = jax.jit(m.initial_data, static_argnums=0)(
                id_type, X, Y, Z, params, **{n: rt[n] for n in needs})
            fields_ok = set(out) == want
            shapes_ok = all(v.shape == shape for v in out.values())
            finite = all(np.isfinite(np.asarray(v)).all() for v in out.values())
            good = fields_ok and shapes_ok and finite
            print(f"{'PASS' if good else 'FAIL'}  id {id_type:<3} {nm[:38]:<38} "
                  f"{len(out)} fields"
                  + ("" if fields_ok else f", MISMATCH {sorted(set(out) ^ want)}")
                  + ("" if shapes_ok else ", BAD SHAPE")
                  + ("" if finite else ", NONFINITE"))
            ok &= good
        for id_type, why in sorted(getattr(m, "INITIAL_DATA_UNAVAILABLE", {}).items()):
            try:
                m.initial_data(id_type, X, Y, Z, params)
            except NotImplementedError:
                print(f"PASS  id {id_type:<3} refused with a reason, as it should")
            else:
                print(f"FAIL  id {id_type} silently returned something")
                ok = False

    # constraint enforcement: names only, but they must name real fields
    if hasattr(m, "METRIC_VARS"):
        every = {n for vt in vts for n in getattr(m, f"{vt.upper()}_FIELDS")}
        named = (set(m.METRIC_VARS)
                 | {n for g in m.TRACE_FREE_VARS for n in g}
                 | {n for n, _ in m.POS_FLOOR_VARS})
        enf_ok = named <= every and len(m.METRIC_VARS) == 6
        print(f"{'PASS' if enf_ok else 'FAIL'}  enforcement  "
              f"{len(m.METRIC_VARS)} metric, {len(m.TRACE_FREE_VARS)} trace-free, "
              f"{len(m.POS_FLOOR_VARS)} floored"
              + ("" if enf_ok else f", not fields: {sorted(named - every)}"))
        ok &= enf_ok

    print("\nall gate F checks passed" if ok else "\ngate F FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
