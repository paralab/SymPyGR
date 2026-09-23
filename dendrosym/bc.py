"""Outer boundary conditions for generated solvers.

A config declares one or more named *boundary modes*; each maps evolved
variables (or "all") to a kind, and the par file picks the mode at run time
(`[boundary] mode = "<name>"`, default: the first). The kinds are implemented in
dendrolib (`include/boundary_conditions.h`, namespace dendro_bc); the generator
only wires them.

    cfg.boundary_modes = {
        "outflow": {"all": outflow(wave_outflow_velocity(K))},
        "sommerfeld": {u: sommerfeld(1.0, 0.0), Pi: sommerfeld(2.0, 0.0)},
        "walls": {"all": neumann()},
    }

Kinds (all act on the physical-face nodes of each boundary block, through the
rhs, so the time integrator carries them):

  static()                   hold the face at its current value (f_t = 0)
  dirichlet(g)               time-dependent Dirichlet: g is a sympy expression in
                             x, y, z, t (and config parameters); the generator
                             takes d/dt itself. The face follows g if the initial
                             data already equals g(t = 0) there -- otherwise the
                             offset is kept, not corrected.
  dirichlet_value(g)         Dirichlet by overwriting: the face is set to g on
                             every stage state (at the stage time) and after each
                             step. g is a sympy expression or a C++ expression
                             string in x, y, z, t -- for data with no closed-form
                             d/dt (tables, other codes). Exact on the face; with
                             time-dependent g it can cost RK some time order at the
                             boundary, which the rhs form (dirichlet) does not.
  exact()                    dirichlet() from the id type's exact solution (tests)
  outflow(vel)               f_t = -v . grad f, v the outgoing characteristic
                             velocity per face (6 x 3, sympy allowed); exact for
                             waves meeting the face head-on, partial reflection
                             for oblique ones
  neumann()                  normal derivative held at its initial value (zero
                             flux if the initial data has zero flux)
  robin(a, b)                a f + b d_n f held at its initial value, d_n the
                             outward normal derivative; a, b scalars or 6 per-face
                             values in OCT_DIR order. b = 0 holds f; a = 0 is neumann
  sommerfeld(falloff, f_inf) radiative condition for fields falling off as
                             f_inf + O(1/r^falloff) from a central source

Data-free: static, outflow, neumann, robin, sommerfeld. Needs data: dirichlet,
dirichlet_value, exact.
"""

import sympy as sym

FACE_NORMALS_BY_OCT_DIR = ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1))

KINDS = ("static", "dirichlet", "dirichlet_value", "exact", "outflow", "neumann",
         "robin", "sommerfeld")


def static():
    return ("static",)


def dirichlet(g):
    return ("dirichlet", sym.sympify(g))


def dirichlet_value(g):
    return ("dirichlet_value", g if isinstance(g, str) else sym.sympify(g))


def exact():
    return ("exact",)


def outflow(vel):
    vel = [[sym.sympify(c) for c in row] for row in vel]
    if len(vel) != 6 or any(len(r) != 3 for r in vel):
        raise ValueError("outflow velocity must be 6 faces x 3 components")
    return ("outflow", vel)


def neumann():
    return ("neumann",)


def robin(a, b):
    def per_face(c):
        c = list(c) if isinstance(c, (list, tuple)) else [c] * 6
        if len(c) != 6:
            raise ValueError("robin coefficients: a scalar or 6 per-face values")
        return [sym.sympify(x) for x in c]
    return ("robin", per_face(a), per_face(b))


def sommerfeld(falloff, f_inf):
    return ("sommerfeld", sym.sympify(falloff), sym.sympify(f_inf))


def wave_outflow_velocity(K=None, c=1):
    """Per-face outgoing velocity for u_tt = K_ij d_i d_j u (or speed c if K is
    None): v = K n / sqrt(n.K.n) with n the outward normal, faces in dendrolib's
    OCT_DIR order (-x, +x, -y, +y, -z, +z). For K = c^2 I this is c n."""
    out = []
    for n in FACE_NORMALS_BY_OCT_DIR:
        nv = sym.Matrix(n)
        if K is None:
            v = c * nv
        else:
            K = sym.Matrix(K)
            v = K * nv / sym.sqrt((nv.T * K * nv)[0])
        out.append([v[0], v[1], v[2]])
    return out
