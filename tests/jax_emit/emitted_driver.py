"""Turn a dendrosym-emitted JAX package into DendroJAX ``Simulation`` callbacks.

    cbs = build_callbacks(ccz4, domain=..., eleorder=6, ko_sigma=..., ...)
    sim = dendrojax.Simulation(initial_data=cbs.initial_data, rhs=cbs.rhs, ...)

This is the ``build_rhs`` of HANDOFF_emitted_solver_api.md section 3 and belongs
in ``dendrojax.emitted`` -- it is here only because dendrojax is frozen for
optimization work. Nothing in it is CCZ4-specific: every name comes off the
emitted tables.

The contract it bridges, learned from ``Simulation._wrap_rhs`` and
``bssn_rhs_optimized``:

  * ``rhs(f_block, *, hx, hy, hz, cx, cy, cz, bflag, t, dt, params)``
  * ``f_block`` is the PADDED (V, S, S, S) block, S = 2*eleorder+1
  * the return is the TIGHT interior (V, E, E, E), E = eleorder+1
  * ``hx`` is world spacing but ``cx`` is an octree [0,1] corner
  * KO and the Sommerfeld rows are the driver's job -- the emitted body has
    neither (``kograd`` appears nowhere in it)
"""

from typing import NamedTuple

import jax
import jax.numpy as jnp

import dendrojax
from dendrojax import emitted as em


class Params(NamedTuple):
    """What ``Simulation`` passes through to the RHS.

    The emitted params carry no KO strength (KO is a C++-side runtime knob, not
    an equation parameter), so it rides alongside rather than being invented as
    a field of the emitted NamedTuple.
    """

    phys: object
    ko_sigma: float = 0.0


class Callbacks(NamedTuple):
    initial_data: object
    rhs: object
    rhs_interior: object
    num_vars: int
    fields: tuple
    enforce: object
    enforce_state: object


def build_callbacks(mod, *, domain, phys, eleorder=6, id_type=None,
                    runtime=None, var_type="evolution", chi_floor=1e-4,
                    alpha_floor=1e-4, ko_order=4):
    """Resolve the emitted name tables to indices once, outside ``jit``."""
    vt = var_type.upper()
    fields = list(getattr(mod, f"{vt}_FIELDS"))
    outputs = list(getattr(mod, f"{vt}_OUTPUTS"))
    derivs = list(getattr(mod, f"{vt}_DERIVS"))
    bcs = list(getattr(mod, f"{vt}_BCS"))
    extra = list(getattr(mod, f"{vt}_EXTRA_INPUTS"))
    rhs_fn = getattr(mod, f"{var_type}_rhs")

    if extra:
        raise NotImplementedError(f"{vt}_EXTRA_INPUTS not wired: {extra}")

    plan = em.plan_derivs(derivs, fields)          # remesh-stable, hashable
    out_idx = tuple(outputs.index(f"{v}_rhs") for v in fields)
    idx = {v: i for i, v in enumerate(fields)}
    bc_rows = tuple((idx[f], gx, gy, gz, fo, asy) for _, f, gx, gy, gz, fo, asy in bcs)

    metric_idx, tf_idx, floor_names = em.conformal_indices(
        fields, mod.METRIC_VARS, mod.TRACE_FREE_VARS, mod.POS_FLOOR_VARS)
    floors = tuple((fields.index(v), {"chi_floor": chi_floor}.get(p, alpha_floor))
                   for v, p in mod.POS_FLOOR_VARS)

    (xmin, xmax), (ymin, ymax), (zmin, zmax) = domain
    grad_order = 2 * (eleorder // 2)
    # KO carries its own order: dendrojax ships one stencil (KO_WEIGHTS[4]),
    # and passing the derivative grad_order here is a KeyError, not a fallback
    assert ko_order == 4, "dendrojax has only the 4th-order KO stencil"
    pw = eleorder // 2
    S = 2 * eleorder + 1
    E = eleorder + 1
    interior = slice(pw, pw + E)

    def initial_data(X, Y, Z):
        """-> (E, E, E, V); Simulation transposes to (V, E, E, E)."""
        st = mod.initial_data(id_type, X, Y, Z, phys, **(runtime or {}))
        return jnp.stack([jnp.broadcast_to(jnp.asarray(st[v]), X.shape)
                          for v in fields], axis=-1)

    def _body(f_block, hx, hy, hz, bflag, params):
        d = em.compute_derivs(f_block, hx, hy, hz, bflag, plan,
                              eleorder=eleorder, grad_order=grad_order)
        u = {v: f_block[i] for i, v in enumerate(fields)}
        out = rhs_fn(u, d, params.phys)
        return jnp.stack([out[i] for i in out_idx], axis=0), d

    def rhs(f_block, *, hx, hy, hz, cx, cy, cz, bflag, t, dt, params):
        stack, d = _body(f_block, hx, hy, hz, bflag, params)

        # world coords of the PADDED block: the ghost layer sits at negative
        # offsets, so the -pw is what puts the outer face at the domain edge
        xs = (xmin + cx * (xmax - xmin)) + (jnp.arange(S) - pw) * hx
        ys = (ymin + cy * (ymax - ymin)) + (jnp.arange(S) - pw) * hy
        zs = (zmin + cz * (zmax - zmin)) + (jnp.arange(S) - pw) * hz

        def apply_bcs(cur):
            for i, gx, gy, gz, fo, asy in bc_rows:
                cur = cur.at[i].set(em.apply_bc_row(
                    cur[i], f_block[i], d[gx], d[gy], d[gz],
                    xs, ys, zs, bflag, fo, asy, pw))
            return cur

        stack = jax.lax.cond(bflag > 0, apply_bcs, lambda s: s, operand=stack)
        ko = dendrojax.derivs.compute_ko_diss_xyz(
            f_block, hx, hy, hz, bflag, eleorder, ko_order)
        stack = stack + ko * params.ko_sigma
        return stack[:, interior, interior, interior]

    def rhs_interior(f_block, *, hx, hy, hz, cx, cy, cz, t, dt, params):
        """bflag == 0 batches: no BC branch, no one-sided KO closures."""
        stack, _ = _body(f_block, hx, hy, hz, 0, params)
        ko = dendrojax.derivs.compute_ko_diss_xyz_interior(
            f_block, hx, hy, hz, eleorder, ko_order)
        stack = stack + ko * params.ko_sigma
        return stack[:, interior, interior, interior]

    def enforce(block):
        """det(g~) -> 1, trace-free A~, positive chi/alpha. One (V,...) block."""
        return em.enforce_conformal_constraints(block, metric_idx, tf_idx, floors)

    def enforce_state(state):
        """``state -> state`` for run_n_steps_split's ``post_step_func``.

        Built once and reused: post_step_func is a JIT cache key, so handing in
        a fresh closure per call re-traces the whole step loop.
        """
        return state._replace(
            physics_data=jax.vmap(enforce)(state.physics_data))

    return Callbacks(initial_data, rhs, rhs_interior, len(fields),
                     tuple(fields), enforce, enforce_state)
