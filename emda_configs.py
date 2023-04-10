"""EMDA-BSSN Core Variable generation

Uses the Sympy framework to generate the variables
"""

import dendrosym
import sympy as sym

ONLY_BSSN_VARIABLES = True

# shortcuts for using these specific equations
pi = sym.pi
exp = sym.exp
pow = sym.Pow

# DEFINE: dendro config class to use for generating code
dendroConfigs = dendrosym.NRConfig("emda")

# the indexing used for the dendro configs
idx_str = "[pp]"

# save the index string
dendroConfigs.set_idx_str(idx_str)

# ==========
# PARAMETER SPECIFICATION
#
# these are parameters that can be freely chosen for the supplemental gauge
# conditions
# ===
# sets the lambda params which are involved with the gauge
lambda_param = dendrosym.dtypes.ParameterVariable(
    "lambda", dtype="unsigned int", num_params=4
)
lambda1, lambda2, lambda3, lambda4 = lambda_param.get_symbolic_repr()

# sets the lf parameters which are also involved with the gauge.
lf_param = dendrosym.dtypes.ParameterVariable(
    "lf", dtype="double", num_params=2
)
lf0, lf1 = lf_param.get_symbolic_repr()

# sets the eta parameters which come in with the constraint damping
eta_param = dendrosym.dtypes.ParameterVariable(
    "eta", dtype="double", num_params=2
)
eta1, eta2 = eta_param.get_symbolic_repr()

# sets eta damping which is involved with gaugeB
etadamp_param = dendrosym.dtypes.ParameterVariable("etadamp", dtype="double")
etadamp = etadamp_param.get_symbolic_repr()

# These are the theory parameters for EMDA. For plain Kerr-Newman, these will
# both be 0.  For the Kerr-Sen black hole (low energy string theory),
# both of these will be 1.
alpha_param = dendrosym.dtypes.ParameterVariable(
    "alpha_theory", dtype="double", num_params=2
)
alpha0, alpha1 = alpha_param.get_symbolic_repr()

# == END PARAMETERS ==

param_info = []  # TODO: store defaults?
dendroConfigs.add_parameter_variables(
    [lambda_param, eta_param, etadamp_param, alpha_param, lf_param],
    "evolution",
)
dendroConfigs.add_parameter_variables([alpha_param], "constraint")

# ==========
# CONSTRAINT VARIABLES
#
# These are the variables used when formulating the constraint equations. This
# includes the Hamiltonian constraint and the 3 components of the momentum
# constraint. We will also lump in to this the variables associated with Psi_4
# (even though its calculation is not a constraint per se).
# ===
ham = dendrosym.dtypes.scalar("ham" + idx_str)
mom = dendrosym.dtypes.vec3("mom" + idx_str)
psi4_real = dendrosym.dtypes.scalar("psi4_real" + idx_str)
psi4_imag = dendrosym.dtypes.scalar("psi4_imag" + idx_str)
dendroConfigs.add_constraint_variables([ham, mom, psi4_real, psi4_imag])

# [EWH]  We are going to move entirely to using chi.  We will not use psi^p
#   As a result, there will be no need for p_expo.  I am commenting these out.
# p_expo_param = dendrosym.dtypes.ParameterVariable("p_expo", dtype="double")
# p_expo = p_expo_param.get_symbolic_repr()
# dendroConfigs.add_parameter_variables([p_expo_param], "constraint")

# ==========
# PHYSICAL CONSTANTS
#
# Various fractions that are assigned later to the different weights and
# equations
# ===
one_half = sym.Rational(1, 2)
one_third = sym.Rational(1, 3)
three_fourths = sym.Rational(3, 4)
two_thirds = sym.Rational(2, 3)

weight = -two_thirds
weight_Gh = two_thirds

# the lie derivative weights
# [EWH]  These are redundant with the above.
# weight_lie = -two_thirds
# weight_Gammah = two_thirds


# == END CONSTANTS ==

# ==========
# DERIVATIVE SPECIFICATION
#
# These are the derivative functions that will be used throughout the program
# setting them stores them inside the dendro package
# ===
# first derivative is the gradient, and first argument is direction
d_ = dendrosym.nr.set_first_derivative("grad")
# second derivative is second order gradient, and first 2
# arguments are direction
# d2s_ = dendrosym.nr.set_second_derivative('grad2')
d2_ = dendrosym.nr.set_second_derivative("grad2")
# advective derivate, first argument is direction
ad_ = dendrosym.nr.set_advective_derivative("agrad")
# and then we set the kreiss oliger dissipation
kod_ = dendrosym.nr.set_kreiss_oliger_dissipation("kograd")
# == END DERIVATIVES ==

# ==========
# VARIABLE INITIALIZATION
#
# These are the variables that will be used in the dendro C++
# code. They need to be initialized and include their indexing string
# ===

# lapse
alpha = dendrosym.dtypes.scalar("alpha" + idx_str)

# shift
beta = dendrosym.dtypes.vec3("beta" + idx_str)

# extrinsic curvature trace
trK = dendrosym.dtypes.scalar("trK" + idx_str)

# conformal factor
chi = dendrosym.dtypes.scalar("chi" + idx_str)

# shift
Gt = dendrosym.dtypes.vec3("CAP_Gt" + idx_str)

dendroConfigs.add_evolution_variables([alpha, beta, trK, chi, Gt])

dendroConfigs.set_advective_derivative_var(beta)

# gaugeB (capital b!)
gaugeB = dendrosym.dtypes.vec3("gaugeB" + idx_str)
dendroConfigs.add_evolution_variables(gaugeB)

# A tilde
At = dendrosym.dtypes.sym_3x3("At" + idx_str)
dendroConfigs.add_evolution_variables(At)

# gamma tilde (the conformally rescaled metric); renamed gt
gt = dendrosym.dtypes.sym_3x3("gt" + idx_str)
dendroConfigs.add_evolution_variables(gt)  # gammat -> gt

# The dilaton field; an additional scalar field in EMDA theory
dilatonPhi = dendrosym.dtypes.scalar("dilatonPhi" + idx_str)
dendroConfigs.add_evolution_variables(dilatonPhi)

# The axion field; an additional (pseudo-)scalar field in EMDA theory;
# (Note that in some previous codes, we used the name "kappa" for something
# very different associated with constraint damping.)
kappa = dendrosym.dtypes.scalar("kappa" + idx_str)
dendroConfigs.add_evolution_variables(kappa)

# The (negative) directional (or Lie) derivative of the dilaton field in the
# direction normal to the spatial hypersurfaces.  Used as an auxiliary
# variable to turn the 2nd order pde for dilatonPhi into two first order eqns.
capitalPi = dendrosym.dtypes.scalar("capitalPi" + idx_str)
dendroConfigs.add_evolution_variables(capitalPi)

# The (negative) directional (or Lie) derivative of the axion field in the
# direction normal to the spatial hypersurfaces.  Used as an auxiliary
# variable to turn the 2nd order pde for kappa into two first order eqns.
capitalXi = dendrosym.dtypes.scalar("capitalXi" + idx_str)
dendroConfigs.add_evolution_variables(capitalXi)

# The 3D (spacelike) electric field
perpE = dendrosym.dtypes.vec3("perpE" + idx_str)
dendroConfigs.add_evolution_variables(perpE)

# The 3D (spacelike) magnetic field
perpB = dendrosym.dtypes.vec3("perpB" + idx_str)
dendroConfigs.add_evolution_variables(perpB)

# An extra scalar field to effect damping of electric constraint violations
dampingPsi = dendrosym.dtypes.scalar("dampingPsi" + idx_str)
dendroConfigs.add_evolution_variables(dampingPsi)

# An extra scalar field to effect damping of magnetic constraint violations
dampingPhi = dendrosym.dtypes.scalar("dampingPhi" + idx_str)
dendroConfigs.add_evolution_variables(dampingPhi)


# We here set the (conformal) metric for dendro
dendroConfigs.set_metric(gt)
# and then we get the inverse (conformal) metric
igt = dendrosym.nr.get_inverse_metric()
# as well as the two Christoffel symbols built from the conformal metric
C1 = dendrosym.nr.get_first_christoffel()
C2 = dendrosym.nr.get_second_christoffel()
# compute phi from the conformal factor which we will be evolving
# [EWH]: C3 is not used in emda-gr and its definition in dendro.py
#        requires the use of chi and not phi.  In addition, phi, while
#        correctly defined here in terms of chi, actually conflicts with
#        the angular coordinate defined as part of the calculation of
#        Psi_4 below.  As a result, I am commenting these out.
# phi = -sym.Rational(1, 4) * sym.log(chi)
# NOTE: [DFV] - C3 needs to be computed to calculate DiDj
C3 = dendrosym.nr.get_complete_christoffel(chi)


# standard Levi-Civita symbol
leviCivita = [
    [[0, 0, 0], [0, 0, 1], [0, -1, 0]],
    [[0, 0, -1], [0, 0, 0], [1, 0, 0]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 0]],
]

# set up all 10 parts of the stress energy tensor
perpTpart1 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            2 * d_(i, dilatonPhi) * d_(j, dilatonPhi)
            for i, j in dendrosym.nr.e_ij
        ]
    )
)

perpTpart2 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            -gt[i, j]
            * sum(
                [
                    d_(l, dilatonPhi) * d_(m, dilatonPhi) * igt[l, m]
                    for l in dendrosym.nr.e_i
                    for m in dendrosym.nr.e_i
                ]
            )
            for i, j in dendrosym.nr.e_ij
        ]
    )
)  # increase efficiency by moving sum?!

perpTpart3 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [gt[i, j] * pow(capitalPi, 2.0) / chi for i, j in dendrosym.nr.e_ij]
    )
)

perpTpart4 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            2
            * exp(-2 * alpha0 * dilatonPhi)
            / (pow(chi, 2.0))
            * sum(
                [
                    perpB[k]
                    * perpB[m]
                    * leviCivita[i][j][k]
                    * leviCivita[l][f][m]
                    * igt[j, f]
                    for k in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                    for m in dendrosym.nr.e_i
                    for f in dendrosym.nr.e_i
                ]
            )
            for i, l in dendrosym.nr.e_ij
        ]
    )
)

perpTpart5 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            -2
            * exp(-2 * alpha0 * dilatonPhi)
            / (chi**2)
            * sum(
                [
                    gt[i, k] * gt[j, l] * perpE[k] * perpE[l]
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            for i, j in dendrosym.nr.e_ij
        ]
    )
)
perpTpart6 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            -exp(-2 * alpha0 * dilatonPhi)
            / (chi**2)
            * sum(
                [
                    perpB[l] * perpB[m] * gt[l, m]
                    for l in dendrosym.nr.e_i
                    for m in dendrosym.nr.e_i
                ]
            )
            * gt[i, j]
            for i, j in dendrosym.nr.e_ij
        ]
    )
)
perpTpart7 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            exp(-2 * alpha0 * dilatonPhi)
            / (chi**2)
            * sum(
                [
                    perpE[l] * perpE[m] * gt[l, m]
                    for l in dendrosym.nr.e_i
                    for m in dendrosym.nr.e_i
                ]
            )
            * gt[i, j]
            for i, j in dendrosym.nr.e_ij
        ]
    )
)
perpTpart8 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            sym.Rational(1, 2)
            * exp(4 * alpha1 * dilatonPhi)
            * d_(i, kappa)
            * d_(j, kappa)
            for i, j in dendrosym.nr.e_ij
        ]
    )
)
perpTpart9 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            -sym.Rational(1, 4)
            * exp(4 * alpha1 * dilatonPhi)
            * gt[i, j]
            * sum(
                [
                    d_(l, kappa) * d_(m, kappa) * igt[l, m]
                    for l in dendrosym.nr.e_i
                    for m in dendrosym.nr.e_i
                ]
            )
            for i, j in dendrosym.nr.e_ij
        ]
    )
)
perpTpart10 = (
    1
    / (8 * pi)
    * sym.Matrix(
        [
            sym.Rational(1, 4)
            * exp(4 * alpha1 * dilatonPhi)
            * gt[i, j]
            * capitalXi**2
            / chi
            for i, j in dendrosym.nr.e_ij
        ]
    )
)

perpTpart1 = perpTpart1.reshape(3, 3)
perpTpart2 = perpTpart2.reshape(3, 3)
perpTpart3 = perpTpart3.reshape(3, 3)
perpTpart4 = perpTpart4.reshape(3, 3)
perpTpart5 = perpTpart5.reshape(3, 3)
perpTpart6 = perpTpart6.reshape(3, 3)
perpTpart7 = perpTpart7.reshape(3, 3)
perpTpart8 = perpTpart8.reshape(3, 3)
perpTpart9 = perpTpart9.reshape(3, 3)
perpTpart10 = perpTpart10.reshape(3, 3)

# assemble stress energy tensor, shows up in At_rhs
perpT = sym.Matrix(
    [
        perpTpart1[i, j]
        + perpTpart2[i, j]
        + perpTpart3[i, j]
        + perpTpart4[i, j]
        + perpTpart5[i, j]
        + perpTpart6[i, j]
        + perpTpart7[i, j]
        + perpTpart8[i, j]
        + perpTpart9[i, j]
        + perpTpart10[i, j]
        for i, j in dendrosym.nr.e_ij
    ]
)
perpT = perpT.reshape(
    3, 3
)  # very key. Before, it would return a 9x1 matrix I think
# stress energy density, shows up in trK_rhs
rho = (
    1
    / (8 * pi)
    * (
        pow(capitalPi, 2.0)
        + chi
        * sum(
            [
                d_(i, dilatonPhi) * igt[i, j] * d_(j, dilatonPhi)
                for i in dendrosym.nr.e_i
                for j in dendrosym.nr.e_i
            ]
        )
        + exp(-2 * alpha0 * dilatonPhi)
        * sum(
            [
                perpE[i] * perpE[j] * gt[i, j]
                for i in dendrosym.nr.e_i
                for j in dendrosym.nr.e_i
            ]
        )
        / chi
        + exp(-2 * alpha0 * dilatonPhi)
        * sum(
            [
                perpB[i] * perpB[j] * gt[i, j]
                for i in dendrosym.nr.e_i
                for j in dendrosym.nr.e_i
            ]
        )
        / chi
        + sym.Rational(1, 4)
        * exp(4 * alpha1 * dilatonPhi)
        * pow(capitalXi, 2.0)
        + sym.Rational(1, 4)
        * exp(4 * alpha1 * dilatonPhi)
        * chi
        * sum(
            [
                d_(i, kappa) * d_(j, kappa) * igt[i, j]
                for i in dendrosym.nr.e_i
                for j in dendrosym.nr.e_i
            ]
        )
    )
)
# stress energy current, shows up in Gt_rhs
stressCurrent = [
    1
    / (8 * pi)
    * 2
    * capitalPi
    * sum([igt[i, j] * d_(j, dilatonPhi) * chi for j in dendrosym.nr.e_i])
    + 1
    / (8 * pi)
    * sym.Rational(1, 2)
    * exp(4 * alpha1 * dilatonPhi)
    * capitalXi
    * chi
    * sum([igt[i, j] * d_(j, kappa) for j in dendrosym.nr.e_i])
    + 1
    / (8 * pi)
    * 2
    * (chi) ** (-sym.Rational(1, 2))
    * exp(-2 * alpha0 * dilatonPhi)
    * sum(
        [
            leviCivita[i][j][k] * gt[j, l] * perpE[l] * gt[k, m] * perpB[m]
            for j in dendrosym.nr.e_i
            for k in dendrosym.nr.e_i
            for l in dendrosym.nr.e_i
            for m in dendrosym.nr.e_i
        ]
    )
    for i in dendrosym.nr.e_i
]

# [EWH]:  I don't think these two quantities (the covariant version of the
#         normal vector to the spacelike hypersurfaces and the mean extrinsic
#         curvature used in the gauge conditions) get used, so I will comment
#         them out.
# nc = (-a, 0, 0, 0)
# trK0 = 0


# NOTE: f is actually a function of alpha, according to Dr. Hirschmann
# but in the original paper they apparently set it to 1, in regards to
# the equation for beta_rhs below (they techically set it to 3/4, and don't
# have the 3/4 constant out front)
# [EWH]:  f(alpha) shows up as a function in the gauge conditions.  In
#         particular, it shows up in the Gamma driver shift condition as
#         a term in the evolution equation for the shift.  In most
#         implementations, it is just a constant (usually 1) but in some
#         early work on BBH simulations, it was set proportional to the
#         lapse, alpha. In an effort to capture both cases, we have defined
#         it as f(alpha) = lf0 + lf1*alpha, i.e. linear in alpha with the
#         constants lf0 and lf1 (specified above) either 0 or 1 to be able
#         to toggle between these two cases.  As we hard code this below,
#         the definition here is superfluous.  I will comment it out here.
# def f(a):
#    return 1


# The conformal connection functions.
# We will need a storage Gt variable for use with dendro, even if it's
# calculated and used later.
# Gt_var = dendrosym.dtypes.vec3("CAP_Gt" + idx_str)
# dendroConfigs.add_evolution_variables(Gt_var)

# this is the full calculation of it, for use in the equations below
# Gt = tuple([
#    -sum([d_(j, igt[i, j]) for j in dendrosym.nr.e_i])
#    for i in dendrosym.nr.e_i
# ])

# calculate R
# [DFVK]: I changed phi to chi here for computing the ricci, EWH please confirm
R, Rt, Rphi, CalGt = dendrosym.nr.compute_ricci(Gt, chi)
R_trace = dendrosym.nr.trace(
    R
)  # this should be (3)R or just R with no indices

# compute the trace of the stress energy tensor, this shows up in trK_rhs
perpTTrace = dendrosym.nr.trace(perpT)

# we also need to calculate psi on its own
# NOTE: before, psi was initialized to be a "scalar" + idxstr
# this doesn't work because Psi is actually det(gt)**(1/(3*P_EXPO))
# psi = dendrosym.dtypes.scalar("psi" + idx_str)
# Sebastian - I actually made a mistake, sorry about that.
# since chi * gd = gt, and gd = psi**(p_expo) gt, we need chi**(-3) = psi**(3 p_expo),
# so psi = chi**(- 1/p_expo).
# psi = chi**(-1/p_expo)

# == END VARIABLE INITIALIZATION ==

# ====================================================================
# ====================================================================
# =================== EVOLUTION EQUATIONS ============================
# ====================================================================

# throw these inside a function so that it doesn't get generated right away


def evolution_rhs_eqns():

    # ==========
    # GAMMA (tilde) RHS
    # ===
    # gamma tilde's right hand side, (14)

    gt_rhs = (
        dendrosym.nr.lie(beta, gt, weight) - 2 * alpha * At
    )  # can I add this

    # == END gAMMA (tilde) RHS

    # ==========
    # At RHS
    # ===
    AikAkj = sym.Matrix(
        [
            sum(
                [
                    At[ii, kk]
                    * sum(
                        [
                            dendrosym.nr.inv_metric[kk, ll] * At[ll, jj]
                            for ll in dendrosym.nr.e_i
                        ]
                    )
                    for kk in dendrosym.nr.e_i
                ]
            )
            for ii, jj in dendrosym.nr.e_ij
        ]
    ).reshape(3, 3)

    # At_rhs = dendrosym.nr.lie(beta, At, weight) + chi * dendrosym.nr.trace_free(
    #    alpha * R - dendrosym.nr.DiDj(alpha) -
    #    8 * pi * alpha * perpT) + alpha * (trK * At - 2 * AikAkj)
    # At_rhs = ( dendrosym.nr.lie(beta, At, weight) ) + chi * dendrosym.nr.trace_free( alpha*At - dendrosym.nr.DiDj(alpha) - 0.0 * 8 * pi * alpha * perpT     ) + alpha * (trK * At - 2 * AikAkj)
    At_rhs = (
        (dendrosym.nr.lie(beta, At, weight))
        + chi
        * dendrosym.nr.trace_free(
            alpha * R - dendrosym.nr.DiDj(alpha) - 0.0 * 8 * pi * alpha * perpT
        )
        + alpha * (trK * At - 2 * AikAkj)
    )

    # == END At RHS ==

    # ==========
    # CHI RHS
    # ===
    # chi's right hand side
    # Here we operate in terms of chi.
    chi_rhs = dendrosym.nr.lie(beta, chi, weight) + sym.Rational(2, 3) * (
        chi * alpha * trK
    )

    # == END CHI RHS ==

    # ==========
    # trK RHS
    # ===
    # trK rhs
    # NOTE: "covariant divergence" might be incorrect!
    trK_rhs = (
        dendrosym.nr.lie(beta, trK)
        - dendrosym.nr.laplacian(alpha, chi)
        + alpha * (trK * trK / 3 + dendrosym.nr.sqr(At))
        + 4 * pi * alpha * (rho + perpTTrace)
    )
    # == END K RHS ==

    # ==========
    # Gamma Hat RHS
    # ===
    # get the up up version of At
    At_UU = dendrosym.nr.up_up(At)
    # compute Gt, the contracted conformal connection coefficient
    Gt_rhs = (
        sym.Matrix(
            (
                [
                    sum(beta[j] * d_(j, Gt[i]) for j in dendrosym.nr.e_i)
                    for i in dendrosym.nr.e_i
                ]
            )
        )
        - sym.Matrix(
            (
                [
                    sum(CalGt[j] * d_(j, beta[i]) for j in dendrosym.nr.e_i)
                    for i in dendrosym.nr.e_i
                ]
            )
        )
        + sym.Rational(2, 3)
        * sym.Matrix(
            (
                [
                    CalGt[i] * sum(d_(j, beta[j]) for j in dendrosym.nr.e_i)
                    for i in dendrosym.nr.e_i
                ]
            )
        )
        + sym.Matrix(
            (
                [
                    sum(
                        [
                            igt[j, k] * d2_(j, k, beta[i])
                            + igt[i, j] * d2_(j, k, beta[k]) / 3
                            for j, k in dendrosym.nr.e_ij
                        ]
                    )
                    for i in dendrosym.nr.e_i
                ]
            )
        )
        - sym.Matrix(
            (
                [
                    sum(
                        [
                            2 * At_UU[i, j] * d_(j, alpha)
                            for j in dendrosym.nr.e_i
                        ]
                    )
                    for i in dendrosym.nr.e_i
                ]
            )
        )
        + sym.Matrix(
            (
                [
                    sum(
                        [
                            2 * alpha * C2[i, j, k] * At_UU[j, k]
                            for j, k in dendrosym.nr.e_ij
                        ]
                    )
                    for i in dendrosym.nr.e_i
                ]
            )
        )
        - sym.Matrix(
            (
                [
                    sum(
                        [
                            alpha
                            * (
                                3 / chi * At_UU[i, j] * d_(j, chi)
                                + sym.Rational(4, 3)
                                * dendrosym.nr.inv_metric[i, j]
                                * d_(j, trK)
                            )
                            for j in dendrosym.nr.e_i
                        ]
                    )
                    for i in dendrosym.nr.e_i
                ]
            )
        )
        + sym.Matrix(
            (
                [
                    -16 * pi * alpha / chi * stressCurrent[i]
                    for i in dendrosym.nr.e_i
                ]
            )
        )
    )  # is dendrosym.nr.inv_metric ok?
    # == END GAMMA TILDE RHS ==

    # dilatonPhi RHS
    # this is the dilaton field, a scalar field
    # dilatonPhi_rhs = -alpha * capitalPi + dendrosym.nr.vec_j_del_j(beta, dilatonPhi)
    dilatonPhi_rhs = -alpha * capitalPi + dendrosym.nr.vec_j_del_j(
        beta, dilatonPhi
    )
    # END dilatonPhi RHS

    # Axion field RHS, we call it kappa
    # kappa_rhs = -alpha * capitalXi + dendrosym.nr.vec_j_del_j(beta, kappa)
    kappa_rhs = -alpha * capitalXi + dendrosym.nr.vec_j_del_j(beta, kappa)
    # end axion field RHS

    # capitalPi_rhs
    if not ONLY_BSSN_VARIABLES:
        # this is the negative derivative of the dilaton field
        # in the normal direction
        capitalPi_rhs = (
            dendrosym.nr.vec_j_del_j(beta, capitalPi)
            + trK * alpha * capitalPi
            - alpha
            * chi
            * sum(
                [
                    igt[i, j] * d2_(i, j, dilatonPhi)
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            - chi
            * sum(
                [
                    igt[i, j] * d_(i, alpha) * d_(j, dilatonPhi)
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            + alpha
            * chi
            * sum([CalGt[i] * d_(i, dilatonPhi) for i in dendrosym.nr.e_i])
            + 1
            / 2
            * alpha
            * sum(
                [
                    igt[i, j] * d_(i, chi) * d_(j, dilatonPhi)
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            - alpha0
            * alpha
            * exp(-2 * alpha0 * dilatonPhi)
            * sum(
                [
                    gt[i, j] * perpB[i] * perpB[j]
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            / chi
            + alpha0
            * alpha
            * exp(-2 * alpha0 * dilatonPhi)
            * sum(
                [
                    gt[i, j] * perpE[i] * perpE[j]
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            / chi
            + 1
            / 2
            * alpha1
            * alpha
            * exp(4 * alpha1 * dilatonPhi)
            * sum(
                [
                    igt[i, j] * d_(i, kappa) * d_(j, kappa)
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            * chi
            - 1
            / 2
            * alpha1
            * alpha
            * exp(4 * alpha1 * dilatonPhi)
            * capitalXi**2
        )
    else:
        capitalPi_rhs = 0.0 * alpha
    # end capitalPi_rhs

    # capitalXi_rhs
    if not ONLY_BSSN_VARIABLES:
        # this is the negative normal derivative of the axion field
        capitalXi_rhs = (
            dendrosym.nr.vec_j_del_j(beta, capitalXi)
            + trK * alpha * capitalXi
            - alpha
            * chi
            * sum(
                [
                    igt[i, j] * d2_(i, j, kappa)
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            - chi
            * sum(
                [
                    igt[i, j] * d_(i, alpha) * d_(j, kappa)
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            + alpha
            * chi
            * sum([CalGt[i] * d_(i, kappa) for i in dendrosym.nr.e_i])
            + 1
            / 2
            * alpha
            * sum(
                [
                    igt[i, j] * d_(i, kappa) * d_(j, chi)
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            - 4
            * alpha
            * alpha1
            * sum(
                [
                    igt[i, j] * d_(i, kappa) * d_(j, dilatonPhi)
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            * chi
            + 4 * alpha * alpha1 * capitalPi * capitalXi
            + 4
            * alpha
            * exp(-4 * alpha1 * dilatonPhi)
            * sum(
                [
                    perpB[i] * perpE[j] * gt[i, j]
                    for i in dendrosym.nr.e_i
                    for j in dendrosym.nr.e_i
                ]
            )
            / chi
        )
    else:
        capitalXi_rhs = 0.0 * alpha
    # end of capitalXi_rhs

    # dampingPsi_rhs
    if not ONLY_BSSN_VARIABLES:
        # this is related to Gauss's equation
        dampingPsi_rhs = (
            dendrosym.nr.vec_j_del_j(beta, dampingPsi)
            - alpha * sum([d_(i, perpE[i]) for i in dendrosym.nr.e_i])
            + 3
            * alpha
            / (2 * chi)
            * sum([perpE[i] * d_(i, chi) for i in dendrosym.nr.e_i])
            + 2
            * alpha0
            * alpha
            * sum([perpE[i] * d_(i, dilatonPhi) for i in dendrosym.nr.e_i])
            + alpha
            * exp(2 * alpha0 * dilatonPhi)
            * sum([perpB[i] * d_(i, kappa) for i in dendrosym.nr.e_i])
            - alpha * eta1 * dampingPsi
        )
    else:
        dampingPsi_rhs = 0.0 * alpha
    # end dampingPsi_rhs

    # dampingPhi_rhs
    if not ONLY_BSSN_VARIABLES:
        # this is related to the "no magnetic monopoles law"
        dampingPhi_rhs = (
            dendrosym.nr.vec_j_del_j(beta, dampingPhi)
            + alpha * sum([d_(i, perpB[i]) for i in dendrosym.nr.e_i])
            - 3
            * alpha
            / (2 * chi)
            * sum([perpB[i] * d_(i, chi) for i in dendrosym.nr.e_i])
            - alpha * eta2 * dampingPhi
        )
    else:
        dampingPhi_rhs = 0.0 * alpha
    # end dampingPhi_rhs

    # perpE_rhs
    if not ONLY_BSSN_VARIABLES:
        perpE_rhs = [
            dendrosym.nr.lie(beta, perpE)[i]
            + alpha * trK * perpE[i]
            + 2
            * alpha
            * alpha0
            * pow(chi, sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k]
                    * gt[k, l]
                    * perpB[l]
                    * d_(j, dilatonPhi)
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            - 2 * alpha * alpha0 * perpE[i] * capitalPi
            + alpha
            * exp(2 * alpha0 * dilatonPhi)
            * pow(chi, sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k] * gt[k, l] * perpE[l] * d_(j, kappa)
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            - alpha * exp(2 * alpha0 * dilatonPhi) * capitalXi * perpB[i]
            - alpha
            * chi
            * sum([igt[i, j] * d_(j, dampingPsi) for j in dendrosym.nr.e_i])
            - pow(chi, sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k] * gt[k, l] * perpB[l] * d_(j, alpha)
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            - alpha
            * pow(chi, sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k] * perpB[l] * d_(j, gt[k, l])
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            - alpha
            * pow(chi, sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k] * gt[k, l] * d_(j, perpB[l])
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            + alpha
            * pow(chi, -sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k] * gt[k, l] * perpB[l] * d_(j, chi)
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            for i in dendrosym.nr.e_i
        ]
    else:
        perpE_rhs = [0.0 * alpha, 0.0 * alpha, 0.0 * alpha]
    # end perpE_rhs

    # perpB_rhs
    if not ONLY_BSSN_VARIABLES:
        perpB_rhs = [
            alpha * trK * perpB[i]
            + alpha
            * chi
            * sum([igt[i, j] * d_(j, dampingPhi) for j in dendrosym.nr.e_i])
            - pow(chi, sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k] * gt[k, l] * perpE[l] * d_(j, alpha)
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            + pow(chi, -sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k]
                    * alpha
                    * gt[k, l]
                    * perpE[l]
                    * d_(j, chi)
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            - alpha
            * pow(chi, sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k] * perpE[l] * d_(j, gt[k, l])
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            - alpha
            * pow(chi, sym.Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k] * gt[k, l] * d_(j, perpE[l])
                    for j in dendrosym.nr.e_i
                    for k in dendrosym.nr.e_i
                    for l in dendrosym.nr.e_i
                ]
            )
            + dendrosym.nr.lie(beta, perpB)[i]
            for i in dendrosym.nr.e_i
        ]
    else:
        perpB_rhs = [0.0 * alpha, 0.0 * alpha, 0.0 * alpha]
    # end perpB_rhs

    # ====================================================================
    # ====================================================================
    # ======================= GAUGE EQUATIONS ============================
    # ====================================================================

    # ==========
    # Alpha RHS
    # ===
    # LAPSE EQUATION

    alpha_rhs = lambda1 * dendrosym.nr.lie(beta, alpha) - 2 * alpha * trK

    # == END ALPHA RHS ==

    # ==========
    # beta RHS
    # ===
    # SHIFT EQUATION

    # beta_rhs = [(sym.Rational(3, 4) * (lf0 + lf1 * alpha) * gaugeB[i] +
    #          lambda2 * dendrosym.nr.vec_j_del_j(beta, beta[i]))
    #         for i in dendrosym.nr.e_i]

    beta_rhs = [
        (
            sym.Rational(3, 4) * (lf0 + lf1 * alpha) * gaugeB[i]
            + lambda2 * dendrosym.nr.vec_j_del_j(beta, beta[i])
        )
        for i in dendrosym.nr.e_i
    ]

    # == END beta RHS ==

    # ==========
    # gaugeB RHS
    # ===
    # NOTE: gaugeB_rhs is a constraint/gauge condition but relies on another
    # formula so it is placed after the evolution eqns
    # gaugeB right hand side equation

    # gaugeB_rhs = [(Gt_rhs[i] - etadamp * gaugeB[i] +
    #          lambda3 * dendrosym.nr.vec_j_del_j(beta, gaugeB[i]) -
    #          lambda4 * dendrosym.nr.vec_j_del_j(beta, Gt[i]))
    #         for i in dendrosym.nr.e_i]

    gaugeB_rhs = [
        (
            Gt_rhs[i]
            - etadamp * gaugeB[i]
            + lambda3 * dendrosym.nr.vec_j_del_j(beta, gaugeB[i])
            - lambda4 * dendrosym.nr.vec_j_del_j(beta, Gt[i])
        )
        for i in dendrosym.nr.e_i
    ]
    # == END B RHS ==

    # ====================================================================
    # ====================================================================
    # SET UP THE RETURNS
    rhs_list = [
        alpha_rhs,
        beta_rhs,
        gt_rhs,
        chi_rhs,
        At_rhs,
        trK_rhs,
        Gt_rhs,
        gaugeB_rhs,
        dilatonPhi_rhs,
        kappa_rhs,
        capitalPi_rhs,
        capitalXi_rhs,
        perpE_rhs,
        perpB_rhs,
        dampingPsi_rhs,
        dampingPhi_rhs,
    ]
    var_list = [
        # alpha, beta, gt, chi, At, trK, Gt_var, gaugeB, dilatonPhi, kappa, capitalPi, capitalXi,
        alpha,
        beta,
        gt,
        chi,
        At,
        trK,
        Gt,
        gaugeB,
        dilatonPhi,
        kappa,
        capitalPi,
        capitalXi,
        perpE,
        perpB,
        dampingPsi,
        dampingPhi,
    ]

    return rhs_list, var_list


def physical_constraint_eqns():
    # define the coordinates
    x, y, z = sym.symbols("x, y, z")

    invsqrt2 = 0.7071067811865475244

    r_vec = sym.Matrix([[x, y, z]])
    theta = sym.Matrix([[x * z, y * z, -(x * x + y * y)]])
    phi = sym.Matrix([[-y, x, 0.0]])

    # Gram-Schmidt for orthonormal basis. This if the non-conformally scaled
    # metric we renamed gamma tilde to gt for short.
    # gd = gt * psi**(p_expo)
    # [EWH]:  alternative using chi instead of psi^p:
    gd = gt / chi

    # for r_vec
    inner_product = sum(
        [
            sum(
                [gd[ii, jj] * r_vec[ii] * r_vec[jj] for ii in dendrosym.nr.e_i]
            )
            for jj in dendrosym.nr.e_i
        ]
    )

    r_vec /= sym.sqrt(inner_product)

    # now for theta
    inner_product_1 = sum(
        [
            sum(
                [gd[ii, jj] * theta[ii] * theta[jj] for ii in dendrosym.nr.e_i]
            )
            for jj in dendrosym.nr.e_i
        ]
    )
    inner_product_2 = sum(
        [
            sum(
                [gd[ii, jj] * r_vec[ii] * theta[jj] for ii in dendrosym.nr.e_i]
            )
            for jj in dendrosym.nr.e_i
        ]
    )

    theta -= inner_product_2 * r_vec
    theta /= sym.sqrt(inner_product_1 - inner_product_2 * inner_product_2)

    # now for phi
    inner_product_1 = sum(
        [
            sum([gd[ii, jj] * phi[ii] * phi[jj] for ii in dendrosym.nr.e_i])
            for jj in dendrosym.nr.e_i
        ]
    )
    inner_product_2 = sum(
        [
            sum([gd[ii, jj] * r_vec[ii] * phi[jj] for ii in dendrosym.nr.e_i])
            for jj in dendrosym.nr.e_i
        ]
    )
    inner_product_3 = sum(
        [
            sum([gd[ii, jj] * theta[ii] * phi[jj] for ii in dendrosym.nr.e_i])
            for jj in dendrosym.nr.e_i
        ]
    )

    phi -= inner_product_2 * r_vec + inner_product_3 * theta
    phi /= sym.sqrt(
        inner_product_1
        - inner_product_2 * inner_product_2
        - inner_product_3 * inner_product_3
    )

    # END TETRAD CONSTRUCTION

    # now we calculate the Weyl scalar, which is for Psi4 wave extraction

    # rename tetrad quantities for Psi_4
    r_np = sym.Matrix([[r_vec[0], r_vec[1], r_vec[2]]])
    m_np_real = sym.Matrix([[theta[0], theta[1], theta[2]]]) * invsqrt2
    m_np_img = sym.Matrix([[phi[0], phi[1], phi[2]]]) * invsqrt2

    At_UD = dendrosym.nr.up_down(At)

    # Auxilary variables, MM and NN, MR and NR
    # MM and NN are 3x3 symmetric matrices
    # MR and NR are 3x3 antisymmetric matrices
    MM = sym.Matrix(
        [
            m_np_real[ii] * m_np_real[jj] - m_np_img[ii] * m_np_img[jj]
            for ii, jj in dendrosym.nr.e_ij
        ]
    )
    MM = MM.reshape(3, 3)

    NN = sym.Matrix(
        [
            m_np_real[ii] * m_np_img[jj] - m_np_real[jj] * m_np_img[ii]
            for ii, jj in dendrosym.nr.e_ij
        ]
    )
    NN = NN.reshape(3, 3)

    # additional intermediate variables
    Atr_vec = [
        sum([At[ii, jj] * r_np[jj] for jj in dendrosym.nr.e_i])
        for ii in dendrosym.nr.e_i
    ]
    Atrr = sum([Atr_vec[ii] * r_np[ii] for ii in dendrosym.nr.e_i])

    # further intermediate variables
    Uu = sym.Matrix(
        [
            sum(
                [
                    r_np[kk]
                    * (
                        d_(kk, At[ii, jj])
                        - sum(
                            [
                                C2[mm, kk, ii] * At[mm, jj]
                                for mm in dendrosym.nr.e_i
                            ]
                        )
                    )
                    for kk in dendrosym.nr.e_i
                ]
            )
            for ii, jj in dendrosym.nr.e_ij
        ]
    )
    Uu = Uu.reshape(3, 3)

    Vv = sym.Matrix(
        [
            sum(
                [
                    r_np[kk]
                    * (
                        d_(jj, At[ii, kk])
                        - sum(
                            [
                                C2[mm, jj, ii] * At[mm, kk]
                                for mm in dendrosym.nr.e_i
                            ]
                        )
                    )
                    for kk in dendrosym.nr.e_i
                ]
            )
            for ii, jj in dendrosym.nr.e_ij
        ]
    )
    Vv = Vv.reshape(3, 3)

    # r_d_psi
    # r_d_psi = sum([r_np[ii] * d_(ii, psi) for ii in dendrosym.nr.e_i])
    # r_d_chi
    r_d_chi = sum([r_np[ii] * d_(ii, chi) for ii in dendrosym.nr.e_i])

    # A temporary
    # A_temp = psi**(p_expo) * (one_third * K + psi**(p_expo) * Atrr -
    #                          0.5 * p_expo * r_d_psi / psi)
    # [EWH]:  alternative using chi instead of psi^p:
    A_temp = (1 / chi) * (one_third * trK + Atrr / chi + 0.5 * r_d_chi / chi)

    # now actually calculate Psi4
    # Psi4_temp = sym.Matrix([
    #    R[ii, jj] + At[ii, jj] * A_temp - Atr_vec[ii] *
    #    (psi**(2 * p_expo) * Atr_vec[jj] - 0.5 * p_expo * psi**
    #     (p_expo - 1) * d_(jj, psi)) - psi**(p_expo) *
    #    (Uu[ii, jj] - Vv[ii, jj]) for ii, jj in dendrosym.nr.e_ij
    # ])
    # [EWH]:  alternative using chi instead of psi^p:
    Psi4_temp = sym.Matrix(
        [
            R[ii, jj]
            + At[ii, jj] * A_temp
            - Atr_vec[ii]
            * ((1 / (chi * chi)) * Atr_vec[jj] + 0.5 * chi * chi * d_(jj, chi))
            - (1 / chi) * (Uu[ii, jj] - Vv[ii, jj])
            for ii, jj in dendrosym.nr.e_ij
        ]
    )

    Psi4_temp = Psi4_temp.reshape(3, 3)

    psi4_real_rhs = sum(
        [Psi4_temp[ii, jj] * MM[ii, jj] for ii, jj in dendrosym.nr.e_ij]
    )
    psi4_imag_rhs = sum(
        [Psi4_temp[ii, jj] * NN[ii, jj] for ii, jj in dendrosym.nr.e_ij]
    )

    # ACTUAL CONSTRAINT EQUATIONS
    # ham_rhs = sum(psi**(-p_expo) * igt[jj, kk] * R[jj, kk]
    #              for jj, kk in dendrosym.nr.e_ij) - dendrosym.nr.sqr(
    #                  At) + two_thirds * K**2 - 16 * pi * rho
    # [EWH]: alternative using chi instead of psi^p
    ham_rhs = (
        sum(chi * igt[jj, kk] * R[jj, kk] for jj, kk in dendrosym.nr.e_ij)
        - dendrosym.nr.sqr(At)
        + two_thirds * trK**2
        - 16 * pi * rho
    )

    mom_rhs = (
        sym.Matrix(
            [
                sum(
                    [
                        igt[jj, kk]
                        * (
                            d_(jj, At[kk, ii])
                            - sum(
                                C2[mm, jj, ii] * At[kk, mm]
                                for mm in dendrosym.nr.e_i
                            )
                        )
                        for jj, kk in dendrosym.nr.e_ij
                    ]
                )
                for ii in dendrosym.nr.e_i
            ]
        )
        - sym.Matrix(
            [
                sum([CalGt[jj] * At[ii, jj] for jj in dendrosym.nr.e_i])
                for ii in dendrosym.nr.e_i
                # ]) + (1.5) * p_expo * sym.Matrix([
                # [EWH]: alternative using chi instead of psi^p
            ]
        )
        - (1.5)
        * sym.Matrix(
            [
                sum(
                    [
                        # igt[jj, kk] * At[kk, ii] * d_(jj, psi) / psi
                        igt[jj, kk] * At[kk, ii] * d_(jj, chi) / chi
                        for jj, kk in dendrosym.nr.e_ij
                    ]
                )
                for ii in dendrosym.nr.e_i
            ]
        )
        - two_thirds * sym.Matrix([d_(ii, trK) for ii in dendrosym.nr.e_i])
        - 8
        * pi
        * sym.Matrix(
            [
                sum([stressCurrent[i] * gt[i, j] for i in dendrosym.nr.e_i])
                for j in dendrosym.nr.e_i
            ]
        )
    )
    mom_rhs = [item for sublist in mom_rhs.tolist() for item in sublist]

    # then we build up the outs

    rhs_list = [psi4_real_rhs, psi4_imag_rhs, ham_rhs, mom_rhs]
    var_list = [psi4_real, psi4_imag, ham, mom]

    return rhs_list, var_list


dendroConfigs.set_rhs_equation_function("evolution", evolution_rhs_eqns)

dendroConfigs.set_rhs_equation_function("constraint", physical_constraint_eqns)

# ==========
# BCS INFORMATION FOR FALLOFF AND ASYMPTOTIC VALUES
# ==========
gt_f_and_a = [
    [[1.0, 1.0], [1.0, 0.0], [1.0, 0.0]],
    [[], [1.0, 1.0], [1.0, 0.0]],
    [[], [], [1.0, 1.0]],
]
At_f_and_a = [2.0, 0.0]
# NOTE: this was originally named phi_f_and_a
chi_f_and_a = [1.0, 1.0]
trK_f_and_a = [1.0, 0.0]
Theta_f_and_a = [1.0, 0.0]
# NOTE: this was originally Gamma_hat_f_and_a
Gt_f_and_a = [2.0, 0.0]
alpha_f_and_a = [1.0, 1.0]
beta_f_and_a = [1.0, 0.0]
gaugeB_f_and_a = [1.0, 0.0]
dilatonPhi_f_and_a = [1.0, 0.0]
capitalPi_f_and_a = [1.0, 0.0]
kappa_f_and_a = [1.0, 0.0]
capitalXi_f_and_a = [1.0, 0.0]
perpB_f_and_a = [1.0, 0.0]
perpE_f_and_a = [1.0, 0.0]
dampingPsi_f_and_a = [1.0, 0.0]
dampingPhi_f_and_a = [1.0, 0.0]
# standard stuff, all of the new gadgets and gismos go to zero
# like 1/r at infinity
#
# we may eventually take the electric and magnetic fields and force them to go
# to zero at like 1/r^2,
# but that assumes that we are not throwing off charged matter

dendroConfigs.set_bhs_falloff_and_asymptotic(
    "evolution",
    [
        # gt, At, chi, trK, Gt_var, alpha, beta, gaugeB, dilatonPhi, kappa, capitalPi, capitalXi,
        gt,
        At,
        chi,
        trK,
        Gt,
        alpha,
        beta,
        gaugeB,
        dilatonPhi,
        kappa,
        capitalPi,
        capitalXi,
        perpE,
        perpB,
        dampingPsi,
        dampingPhi,
    ],
    [
        gt_f_and_a,
        At_f_and_a,
        chi_f_and_a,
        trK_f_and_a,
        Gt_f_and_a,
        alpha_f_and_a,
        beta_f_and_a,
        gaugeB_f_and_a,
        dilatonPhi_f_and_a,
        kappa_f_and_a,
        capitalPi_f_and_a,
        capitalXi_f_and_a,
        perpE_f_and_a,
        perpB_f_and_a,
        dampingPsi_f_and_a,
        dampingPhi_f_and_a,
    ],
)

# add a few evolution constraints
dendroConfigs.add_evolution_constraint(At, "trace_zero")
dendroConfigs.add_evolution_constraint(chi, "pos_floor")
dendroConfigs.add_evolution_constraint(alpha, "pos_floor")
