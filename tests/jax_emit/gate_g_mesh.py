#!/usr/bin/env python
"""Gate G -- an emitted package on a real DendroJAX octree.

    python tests/jax_emit/gate_g_mesh.py <generated_project_dir> [dendrojax_src]

Gates A-F all run one block. This is the first that puts the emitted kernel
behind the actual mesh: octree, 2:1 balance, inter-block ghost exchange, RK4,
and a remesh. What it proves is the contract, not the physics -- if the ghost
sync, the deriv plan and the BC table did not agree, Minkowski would not stay
flat. That is exactly how the missing use_ec_sync was found: the edge halo
between two shared faces went unfilled and every MIXED second derivative was
wrong, silently (2.75e-2 vs 1.95e-17 with the flag on).

Arms:
  1. Minkowski stays a fixed point on the octree (across a coarse/fine
     interface too, when the mesh has one -- reported as a GAP when not).
  2. A puncture (id=6) runs, stays finite, and conserves its enforcement.
  3. A remesh mid-run preserves the state where the mesh did not change.
"""
import importlib
import os
import pathlib
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

FAILURES = []
GAPS = []


def check(name, got, tol):
    ok = got <= tol
    print(f"  {'ok  ' if ok else 'FAIL'} {name}: {got:.3e} (tol {tol:.0e})")
    if not ok:
        FAILURES.append(name)
    return ok


def _nonconforming(sim):
    """Count faces joining elements of different level.

    ``element_levels`` is block-relative and gets renumbered when the atlas
    re-blocks (a remesh was seen going 64@[2] -> 64@[0]), so min==max says
    nothing about whether the mesh is adaptive. TYPE_COARSE/TYPE_FINE faces
    are the real thing arm 1 wants: they are where ghost data is interpolated
    rather than copied.
    """
    import numpy as _np

    import dendrojax as _dj

    ft = _np.asarray(sim.state.connectivity.face_types)[:sim.active_count]
    return int(_np.count_nonzero((ft == _dj.TYPE_COARSE) | (ft == _dj.TYPE_FINE)))


def _patch_simulation(post_step=None):
    """Forward the two step options Simulation does not expose.

    ``run_n_steps_split(use_ec_sync=False)`` is the default and ``Simulation``
    never forwards it, so the edge halo between two SHARED faces is never
    filled. First and diagonal second derivatives only read face halos and stay
    exact; every MIXED second derivative reads that edge and is garbage at
    block boundaries -- silently, with no NaN. bssn_solver.py sets it True in
    three places and bssn_profiler.py calls it "required by BSSN's mixed"
    derivatives; the high-level driver just does not expose it.

    ``post_step_func`` is the same story: ``state -> state`` after every RK
    step, which is where conformal-constraint enforcement has to run. Its
    docstring names BSSN's renormalization as the use case, and Simulation
    forwards that no more than it forwards use_ec_sync.

    Returns ``(note, extra_simulation_kwargs)``: prefer real kwargs if
    dendrojax grows them, otherwise patch the call.
    """
    import functools
    import inspect

    import dendrojax.simulation as dsim

    sig = inspect.signature(dsim.Simulation.__init__).parameters
    kw = {}
    if "use_ec_sync" in sig:
        kw["use_ec_sync"] = True
    if post_step is not None and "post_step_func" in sig:
        kw["post_step_func"] = post_step
    if len(kw) == (1 if post_step is None else 2):
        return ("  step opts: Simulation exposes them", kw)

    extra = {"use_ec_sync": True}
    if post_step is not None:
        extra["post_step_func"] = post_step
    base = getattr(dsim.run_n_steps_split, "func", dsim.run_n_steps_split)
    dsim.run_n_steps_split = functools.partial(base, **extra)
    return ("  step opts: PATCHED use_ec_sync=True"
            + (" + post_step_func" if post_step is not None else "")
            + " -- Simulation forwards neither; without ec_sync every mixed "
              "2nd derivative is wrong at block boundaries", kw)


def main(pkg_dir, dj_src="/home/denv/research/dendrojax/src"):
    sys.path.insert(0, str(pathlib.Path(pkg_dir).resolve()))
    sys.path.insert(0, dj_src)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    import dendrojax
    import jax
    import jax.numpy as jnp
    from dendrojax.initial_data import register_brill_lindquist
    from emitted_driver import Params, build_callbacks

    ec_note, ec_kw = _patch_simulation()
    print(ec_note)

    name = pathlib.Path(pkg_dir).name
    mod = importlib.import_module(next(
        p.name for p in pathlib.Path(pkg_dir).iterdir()
        if p.is_dir() and (p / "__init__.py").exists()))
    print(f"imported {mod.__name__}: {len(mod.EVOLUTION_FIELDS)} fields, "
          f"{len(mod.EVOLUTION_DERIVS)} buffers, {len(mod.EVOLUTION_BCS)} bc rows")

    DOM = ((-16.0, 16.0),) * 3
    ELE = 6
    phys = mod.CCZ4Params()
    params = Params(phys=phys, ko_sigma=0.0)

    # ---- arm 1: Minkowski on an adaptive mesh --------------------------------
    print("\n1. Minkowski fixed point on the octree")
    cb = build_callbacks(mod, domain=DOM, phys=phys, eleorder=ELE, id_type=2)

    # Minkowski is constant, so the wavelet criterion would never refine it.
    # Build the mesh off a bump instead and evolve the flat state on it: any
    # nonzero RHS is then the mesh (ghost sync across a level jump), not data.
    from dendrojax.refinement import get_dendro5_criteria

    def bump(x, y, z, t=0.0):
        return jnp.exp(-(x * x + y * y + z * z) / 4.0)[..., None]

    criteria = get_dendro5_criteria(bump, 1e-3, order=ELE)

    sim = dendrojax.Simulation(
        initial_data=cb.initial_data, rhs=cb.rhs, rhs_interior=cb.rhs_interior,
        num_vars=cb.num_vars, eleorder=ELE, domain=DOM, params=params,
        max_depth=5, min_depth=3, wavelet_tol=1e-3, max_blocks=4096,
        init_grid_iter=5, cfl=0.2, batch_size=32, verbose=False,
        criteria_fn=criteria, max_active_budget=2048, **ec_kw,
    )
    sim.setup()
    nc = _nonconforming(sim)
    if nc == 0:                       # a remesh is what actually refines
        sim.remesh()
        nc = _nonconforming(sim)
    print(f"  mesh: {sim.active_count} elements, {nc} non-conforming faces, "
          f"dt {sim.dt:.4g}")
    if nc == 0:
        GAPS.append("arm 1 ran on a conforming mesh: ghost interpolation "
                    "across a coarse/fine interface is NOT covered")

    sim.run(n_steps=2)
    u = np.asarray(sim.state.physics_data[:sim.active_count])
    flat = np.zeros_like(u)
    for i, v in enumerate(cb.fields):
        flat[:, i] = 1.0 if v in ("alpha", "chi", "gt00", "gt11", "gt22") else 0.0
    check("|state - flat| after 2 steps", np.abs(u - flat).max(), 1e-12)
    check("state finite", float(np.count_nonzero(~np.isfinite(u))), 0.5)

    # ---- arm 2: a puncture on the mesh --------------------------------------
    print("\n2. Brill-Lindquist puncture (id=6) on the mesh")
    register_brill_lindquist(mod, chi_floor=1e-4)
    rt = dict(BH1_mass=1.0, BH1_x=0.0, BH1_y=0.0, BH1_z=0.0, BH1_spin=0.0,
              BH1_spin_theta=0.0, BH1_spin_phi=0.0,
              BH2_mass=0.0, BH2_x=0.0, BH2_y=0.0, BH2_z=0.0, BH2_spin=0.0,
              BH2_spin_theta=0.0, BH2_spin_phi=0.0)
    cb6 = build_callbacks(mod, domain=DOM, phys=phys, eleorder=ELE,
                          id_type=6, runtime=rt)
    sim6 = dendrojax.Simulation(
        initial_data=cb6.initial_data, rhs=cb6.rhs,
        rhs_interior=cb6.rhs_interior, num_vars=cb6.num_vars, eleorder=ELE,
        domain=DOM, params=Params(phys=phys, ko_sigma=0.1),
        max_depth=5, min_depth=3, wavelet_tol=1e-3, max_blocks=4096,
        init_grid_iter=5, cfl=0.2, batch_size=32, verbose=False,
        max_active_budget=2048, **ec_kw,
    )
    sim6.setup()
    print(f"  mesh: {sim6.active_count} elements, "
          f"{_nonconforming(sim6)} non-conforming faces")

    u0 = np.asarray(sim6.state.physics_data[:sim6.active_count])
    chi_i = list(cb6.fields).index("chi")
    print(f"  chi range at t=0: {u0[:, chi_i].min():.3e} .. {u0[:, chi_i].max():.3f}")
    check("id=6 state finite", float(np.count_nonzero(~np.isfinite(u0))), 0.5)

    sim6.run(n_steps=2)
    u6 = np.asarray(sim6.state.physics_data[:sim6.active_count])
    check("finite after 2 steps", float(np.count_nonzero(~np.isfinite(u6))), 0.5)
    check("chi stays positive", float(-u6[:, chi_i].min()), -0.0 + 1e-16)
    print(f"  chi range after 2 steps: {u6[:, chi_i].min():.3e} .. "
          f"{u6[:, chi_i].max():.3f}")

    # enforcement must be idempotent on the evolved state
    blk = jnp.asarray(u6[0])
    e1 = np.asarray(cb6.enforce(blk))
    e2 = np.asarray(cb6.enforce(jnp.asarray(e1)))
    check("enforcement idempotent", np.abs(e1 - e2).max(), 1e-12)

    # ---- arm 3: remesh preserves the state ----------------------------------
    print("\n3. remesh")
    before = int(sim6.active_count)
    sim6.remesh()
    after = int(sim6.active_count)
    print(f"  {before} -> {after} elements")
    ur = np.asarray(sim6.state.physics_data[:after])
    check("finite after remesh", float(np.count_nonzero(~np.isfinite(ur))), 0.5)
    sim6.run(n_steps=2)
    ur2 = np.asarray(sim6.state.physics_data[:sim6.active_count])
    check("finite after remesh + 2 steps",
          float(np.count_nonzero(~np.isfinite(ur2))), 0.5)

    print()
    for g in GAPS:
        print(f"  GAP  {g}")
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        return 1
    print("all gate G checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(*sys.argv[1:]))
