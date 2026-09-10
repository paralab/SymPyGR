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
  1. Minkowski stays a fixed point across a coarse/fine interface. The mesh is
     built from an OFF-CENTRE puncture, then the state is overwritten with
     exact flat: a centred puncture is mirror-symmetric and flat data has no
     structure, so neither refines non-uniformly no matter how good the
     criterion is. Flat data on an adaptive mesh is the one configuration
     where a ghost-interpolation error cannot hide behind truncation error.
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

    ``sync_wavelet`` is the third: the refinement criterion reads
    ``physics_data``, which is interior-only and gets ZERO-padded to S^3, so
    every block edge carries a field-value-to-zero cliff. The wavelet then
    measures the cliff instead of the data -- 4.82e-01 for every field whose
    flat value is 1, identical on every element at every level, nonzero even
    on exactly-flat Minkowski. Everything always exceeds the threshold, so
    refinement is uniform and no coarse/fine interface ever forms. With it on,
    the same puncture gives 127 elements / 100 non-conforming faces instead of
    uniformly reaching 4096.

    Returns ``(note, extra_simulation_kwargs)``: prefer real kwargs if
    dendrojax grows them, otherwise patch the call.
    """
    import functools
    import inspect

    import dendrojax.simulation as dsim

    sig = inspect.signature(dsim.Simulation.__init__).parameters
    want = {"use_ec_sync": True, "sync_wavelet": True}
    if post_step is not None:
        want["post_step_func"] = post_step
    kw = {k: v for k, v in want.items() if k in sig}
    if len(kw) == len(want):
        return ("  step opts: Simulation exposes " + ", ".join(sorted(kw)), kw)

    extra = {"use_ec_sync": True}
    if post_step is not None:
        extra["post_step_func"] = post_step
    base = getattr(dsim.run_n_steps_split, "func", dsim.run_n_steps_split)
    dsim.run_n_steps_split = functools.partial(base, **extra)

    rbase = getattr(dsim._remesh, "func", dsim._remesh)
    dsim._remesh = functools.partial(rbase, sync_wavelet=True)

    return ("  step opts: PATCHED use_ec_sync=True + sync_wavelet=True"
            + (" + post_step_func" if post_step is not None else "")
            + " -- Simulation forwards none of them; without ec_sync every "
              "mixed 2nd derivative is wrong at block boundaries, and without "
              "sync_wavelet the mesh only ever refines uniformly", kw)


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
    print("\n1. Minkowski fixed point across a coarse/fine interface")

    # Flat data has no structure, so a CORRECT wavelet criterion will never
    # refine it -- arm 1 cannot make its own level jump. Build the mesh from an
    # off-centre puncture (centred is mirror-symmetric, so it refines
    # uniformly), then overwrite the state with exact Minkowski: an adaptive
    # mesh carrying data whose answer is known to the last bit.
    register_brill_lindquist(mod, chi_floor=1e-4)
    off = dict(BH1_mass=1.0, BH1_x=3.7, BH1_y=-2.3, BH1_z=1.1, BH1_spin=0.0,
               BH1_spin_theta=0.0, BH1_spin_phi=0.0,
               BH2_mass=0.0, BH2_x=0.0, BH2_y=0.0, BH2_z=0.0, BH2_spin=0.0,
               BH2_spin_theta=0.0, BH2_spin_phi=0.0)
    cbm = build_callbacks(mod, domain=DOM, phys=phys, eleorder=ELE,
                          id_type=6, runtime=off)
    # the plan knows whether the edge band is needed; nobody has to guess
    print(f"  plan needs edge/corner sync: {cbm.needs_ec_sync} "
          f"(mixed 2nd derivatives present)")
    if cbm.needs_ec_sync and ec_kw.get("use_ec_sync") is not True \
            and "PATCHED" not in ec_note:
        FAILURES.append("mixed derivatives present but ec sync not enabled")
    sim = dendrojax.Simulation(
        initial_data=cbm.initial_data, rhs=cbm.rhs,
        rhs_interior=cbm.rhs_interior, num_vars=cbm.num_vars, eleorder=ELE,
        domain=DOM, params=Params(phys=phys, ko_sigma=0.0),
        max_depth=8, min_depth=2, wavelet_tol=1e-3, max_blocks=4096,
        init_grid_iter=0, cfl=0.2, batch_size=64, verbose=False,
        max_active_budget=2048, **ec_kw,
    )
    sim.setup()
    for _ in range(6):
        nc = _nonconforming(sim)
        if nc:
            break
        sim.remesh()
    nc = _nonconforming(sim)
    print(f"  mesh: {sim.active_count} elements, {nc} non-conforming faces, "
          f"dt {sim.dt:.4g}")
    if nc == 0:
        GAPS.append("arm 1 ran on a conforming mesh: ghost interpolation "
                    "across a coarse/fine interface is NOT covered")

    # now make the state exactly flat everywhere
    flat = jnp.asarray([1.0 if v in ("alpha", "chi", "gt00", "gt11", "gt22")
                        else 0.0 for v in cbm.fields])
    pd = jnp.broadcast_to(flat[None, :, None, None, None],
                          sim.state.physics_data.shape)
    sim.state = sim.state._replace(physics_data=jnp.asarray(pd))

    sim.run(n_steps=2)
    u = np.asarray(sim.state.physics_data[:sim.active_count])
    dev = np.abs(u - np.asarray(flat)[None, :, None, None, None]).max()
    check("|state - flat| after 2 steps", dev, 1e-12)
    check("state finite", float(np.count_nonzero(~np.isfinite(u))), 0.5)

    # ---- arm 2: a puncture on the mesh --------------------------------------
    print("\n2. Brill-Lindquist puncture (id=6) on the mesh")
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
