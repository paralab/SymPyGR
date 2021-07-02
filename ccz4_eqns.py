"""CCz4 Core Variable generation

Uses the Sympy framework to generate the variables

The formulation originally used is detailed in the paper
Conformal and covariant formulation of the Z4 system with
constraint-violation damping
available at https://arxiv.org/pdf/1106.2254.pdf
The evolution equations begin with equation 14, with the
constraing equations being 20-22

There is a "note" included in this repository on ccz4 as derived by
Dr. Hirschmann (not fully added)
"""

import dendrosym
import sympy as sym

# DEFINE: dendro config class to use for generating code
dendroConfigs = dendrosym.NRConfig("ccz4")

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
# lambda1, lambda2, lambda3, lambda4, eta = sym.symbols(
#     'lambda[0] lambda[1] lambda[2] lambda[3] eta')
# kappa1, kappa2, kappa3, tau = sym.symbols('kappa[0], kappa[1], kappa[2], tau')
# == END PARAMETERS ==
lambda_param = dendrosym.dtypes.ParameterVariable("lambda",
                                                  dtype="unsigned int",
                                                  num_params=4)
lambda1, lambda2, lambda3, lambda4 = lambda_param.get_symbolic_repr()
kappa_param = dendrosym.dtypes.ParameterVariable("kappa",
                                                 dtype="unsigned int",
                                                 num_params=3)
kappa1, kappa2, kappa3 = kappa_param.get_symbolic_repr()
eta_param = dendrosym.dtypes.ParameterVariable("eta", dtype="double")
eta = eta_param.get_symbolic_repr()
tau_param = dendrosym.dtypes.ParameterVariable("tau", dtype="unsigned int")
tau = tau_param.get_symbolic_repr()

# set some information about these parameters
param_info = []  # TODO: store defaults?
dendroConfigs.add_parameter_variables(
    [lambda_param, kappa_param, tau_param, eta_param], "evolution")

# ==========
# CONSTRAINT VARIABLES
#
# These are the variables used in the formulations as constraints. This
# includes the Hamiltonian constraint, the 3-component momentum constraint,
# as well as the psi4 variable constraint.
# ===
ham = dendrosym.dtypes.scalar("ham" + idx_str)
mom = dendrosym.dtypes.vec3("mom" + idx_str)
psi4_real = dendrosym.dtypes.scalar("psi4_real" + idx_str)
psi4_imag = dendrosym.dtypes.scalar("psi4_imag" + idx_str)
dendroConfigs.add_constraint_variables([ham, mom, psi4_real, psi4_imag])
psi = dendrosym.dtypes.scalar("psi" + idx_str)

p_expo_param = dendrosym.dtypes.ParameterVariable("p_expo", dtype="double")
p_expo = p_expo_param.get_symbolic_repr()
dendroConfigs.add_parameter_variables([p_expo_param], "constraint")

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

# the lie derivative weight
weight_lie = -two_thirds
weight_Gammah = two_thirds
# == END CONSTANTS ==

# ==========
# DERIVATIVE SPECIFICATION
#
# These are the derivative functions that will be used throughout the program
# setting them stores them inside the dendro package
# ===
# first derivative is the gradient, and first argument is direction
d_ = dendrosym.nr.set_first_derivative('grad')
# second derivative is second order gradient, and first 2
# arguments are direction
d2s_ = dendrosym.nr.set_second_derivative('grad2')
# advective derivate, first argument is direction
ad_ = dendrosym.nr.set_advective_derivative('agrad')
# and then we set the dreiss oliger dissipation
kod_ = dendrosym.nr.set_kreiss_oliger_dissipation('kograd')
# == END DERIVATIVES ==

# ==========
# VARIABLE INITIALIZATION
#
# These are the variables that will be used in the dendro C++
# code. They need to be initialized and include their indexing string
# ===

alpha = dendrosym.dtypes.scalar("alpha" + idx_str)
beta = dendrosym.dtypes.vec3("beta" + idx_str)
K = dendrosym.dtypes.scalar("K" + idx_str)
phi = dendrosym.dtypes.scalar("phi" + idx_str)
dendroConfigs.add_evolution_variables([alpha, beta, K, phi])
dendroConfigs.set_advective_derivative_var(beta)

# theta is an expression
Theta = dendrosym.dtypes.scalar("Theta" + idx_str)
dendroConfigs.add_evolution_variables(Theta)

# B (capital b!)
B = dendrosym.dtypes.vec3("B" + idx_str)
dendroConfigs.add_evolution_variables(B)

# A tilde
At = dendrosym.dtypes.sym_3x3("At" + idx_str)
dendroConfigs.add_evolution_variables(At)
# gamma tilde
gammat = dendrosym.dtypes.sym_3x3("gammat" + idx_str)
dendroConfigs.add_evolution_variables(gammat)

# Z
Z = dendrosym.dtypes.vec3("Z" + idx_str)

# S is zero: NOTE: for a vacuum S is just all 0
S = sym.Matrix([[0 for ii in dendrosym.nr.e_i] for jj in dendrosym.nr.e_i])

nc = (alpha, 0, 0, 0)

K0 = 0


# NOTE: f is actually a function of alpha, according to Dr. Hirschmann
# but in the original paper they apparentlyl set it to 1, in regards to
# the equation for beta_rhs below (they techically set it to 3/4, and don't
# have the 3/4 constant out front)
def f(alpha):
    return 1


# set the metric for dendro
# dendrosym.nr.set_metric(gammat)
dendroConfigs.set_metric(gammat)
# and then get the inverse metric
igammat = dendrosym.nr.get_inverse_metric()
C1 = dendrosym.nr.get_first_christoffel()
C2 = dendrosym.nr.get_second_christoffel()
C3 = dendrosym.nr.get_complete_christoffel(phi)

# calculate Gamma tilde from gamma hat
# Gamma hat
Gamma_hat = dendrosym.dtypes.vec3("Gammahat" + idx_str)
# so then gamma tilde is given by:
Gammat = tuple([
    Gamma_hat[ii] - 2 * sum(igammat[ii, jj] * Z[jj] for jj in dendrosym.nr.e_i)
    for ii in dendrosym.nr.e_i
])
dendroConfigs.add_evolution_variables(Gamma_hat)

# calculate R
R, Rt, Rphi, CalGt = dendrosym.nr.compute_ricci(Gammat, phi)
R_trace = dendrosym.nr.trace(
    R)  # this should be (3)R or just R with no indices

# NOTE: this will need to be updated if we want to not use a vacuum
S_trace = dendrosym.nr.trace(S)

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
    M_ij = sym.Matrix([
        sum([gammat[ii, kk] * d_(jj, beta[kk]) for kk in dendrosym.nr.e_i])
        for ii, jj in dendrosym.nr.e_ij
    ]).reshape(3, 3)
    M_ij = 2 * dendrosym.nr.calc_symmetric_part_rank2(M_ij)
    # NOTE: Hirschmann's original python implementation equation doesn't include
    #       the M_ij piece I have here
    gammat_rhs = - 2 * alpha * dendrosym.nr.trace_free(At) + \
        M_ij - \
        two_thirds * gammat * sum(
            [d_(kk, beta[kk]) for kk in dendrosym.nr.e_i]) + \
        sym.Matrix([sum(
            [beta[kk] * d_(kk, gammat[ii, jj]) for kk in dendrosym.nr.e_i]
                    ) for ii, jj in dendrosym.nr.e_ij]).reshape(3, 3)
    # dendrosym.nr.lie(beta, gammat)
    # == END GAMMA (tilde) RHS

    # ==========
    # At RHS
    # ===
    # stack the gradients for dizj and djzi
    DiZj = sym.Matrix([[
        d_(ii, Z[jj]) - sum(C2[kk, ii, jj] * Z[kk] for kk in dendrosym.nr.e_i)
        for jj in dendrosym.nr.e_i
    ] for ii in dendrosym.nr.e_i
                       ])  # TODO: please look up the right cristoffel symbol
    # then the summation part
    AilALj = sym.Matrix([
        sum([
            At[ii, kk] * sum([
                dendrosym.nr.inv_metric[kk, ll] * At[ll, jj]
                for ll in dendrosym.nr.e_i
            ]) for kk in dendrosym.nr.e_i
        ]) for ii, jj in dendrosym.nr.e_ij
    ]).reshape(3, 3)
    # same idea as M_ij above, but this time with Atilde instead of gammatilde
    M_ij = sym.Matrix([
        sum([At[ii, kk] * d_(jj, beta[kk]) for kk in dendrosym.nr.e_i])
        for ii, jj in dendrosym.nr.e_ij
    ]).reshape(3, 3)
    M_ij = 2 * dendrosym.nr.calc_symmetric_part_rank2(M_ij)
    # A tilde i j's right hand side
    At_rhs = phi * phi * dendrosym.nr.trace_free(
        -dendrosym.nr.DiDj(alpha) + alpha * (R + 2 * DiZj - 8 * sym.pi * S)
    ) + alpha * At * (K - 2 * Theta) \
        - 2 * alpha * AilALj \
        + M_ij \
        - two_thirds * At * sum([d_(kk, beta[kk]) for kk in dendrosym.nr.e_i]) \
        + sym.Matrix([
            sum(
                [beta[kk] * d_(kk, At[ii, jj]) for kk in dendrosym.nr.e_i]
                ) for ii, jj in dendrosym.nr.e_ij]).reshape(3, 3)
    # == END At RHS ==

    # ==========
    # PHI RHS
    # ===
    # phi's right hand side
    # NOTE: 2 - the dendro note for CCZ4 (potentially) uses psi for this
    # equation and it's slightly different
    phi_rhs = dendrosym.nr.lie(beta, phi) + \
        one_third * phi * (alpha * K -
                        sum([d_(kk, beta[kk]) for kk in dendrosym.nr.e_i]))
    # NOTE: the old implemenation includes p_expo in the denominator
    # == END PHI RHS ==

    # ==========
    # K RHS
    # ===
    # K rhs
    # NOTE: "covariant divergence" might be incorrect!
    K_rhs = -dendrosym.nr.laplacian(alpha, phi) + alpha * (
        R_trace + 2 * dendrosym.nr.covariant_divergence(Z) + K * K -
        2 * Theta * K) + dendrosym.nr.lie(beta, K) - 3 * alpha * kappa1 * (
            1 + kappa2) * Theta + 4 * sym.pi * alpha * (S_trace - 3 * tau)
    # == END K RHS ==

    # ==========
    # Theta RHS
    # ===
    # Theta rhs
    # NOTE: "covariant divergence" might be incorrect
    Aij_AIJ = dendrosym.nr.sqr(At)
    Theta_rhs = alpha / 2 * (
        R_trace + 2 * dendrosym.nr.covariant_divergence(Z) -
        Aij_AIJ + two_thirds * K * K - 2 * Theta * K) - dendrosym.nr.lie(
            Z, alpha) + dendrosym.nr.lie(beta, Theta) - alpha * kappa1 * (
                2 + kappa2) * Theta - 8 * sym.pi * alpha * tau
    # == END Theta RHS ==

    # ==========
    # Gamma Hat RHS
    # ===
    # now for Gamma hat (oh boy)
    # get the up up version of At
    At_UU = dendrosym.nr.up_up(At)
    # NOTE: Ghat = Gt + 2 gammat^{ij} Z_j
    Gamma_hat_rhs = (
        # the 2 * alpha terms, I've gone ahead and absorbed
        # the alphas into each term
        sym.Matrix([
            sum([2 * alpha * C2[ii, jj, kk] * \
                At_UU[jj, kk] for jj, kk in dendrosym.nr.e_ij]
                ) for ii in dendrosym.nr.e_i]) -
        sym.Matrix([
            sum((6 * alpha / phi) * \
                At_UU[ii, jj] * d_(jj, phi) for jj in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) -
        sym.Matrix([
            sum(2 * two_thirds * \
                alpha * igammat[ii, jj] * d_(jj, K) for jj in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) +
        # the 2 gammat inverse terms
        sym.Matrix([
            sum(2 * alpha * \
                igammat[kk, ii] * d_(kk, Theta) for kk in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) -
        sym.Matrix([
            sum(2 * Theta * igammat[kk, ii] * \
                d_(kk, alpha) for kk in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) -
        sym.Matrix([
            sum(2 * two_thirds * alpha * \
                K * igammat[kk, ii] * Z[kk] for kk in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) -
        # now the 2 A^(ij) dj alpha
        sym.Matrix([
            sum(2 * At_UU[ii, jj] * \
                d_(jj, alpha) for jj in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) +
        # gamma ^kl times the double derivative
        sym.Matrix([
            sum([igammat[jj, kk] * d2s_(jj, kk, beta[ii]) + igammat[ii, jj] * \
                d2s_(jj, kk, beta[kk])/3 for jj, kk in dendrosym.nr.e_ij]
                ) for ii in dendrosym.nr.e_i]) +
        # now for the two thirds * Gammai * derivative of beta in k
        two_thirds * sym.Matrix([
            CalGt[ii] * sum(
                d_(jj, beta[jj]) for jj in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) -
        sym.Matrix([
            sum(CalGt[kk] * d_(kk, beta[ii]) for kk in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) +
        # now for the 2 * kappa term
        2 * kappa3 * (
            two_thirds * sym.Matrix([
                sum(igammat[ii, jj] * Z[jj] * sum(
                        d_(kk, beta[kk]) for kk in dendrosym.nr.e_i
                        ) for jj in dendrosym.nr.e_i
                    ) for ii in dendrosym.nr.e_i]) -
            sym.Matrix([
                sum(igammat[jj, kk] * Z[jj] * d_(kk, beta[ii])
                    for jj, kk in dendrosym.nr.e_ij) for ii in dendrosym.nr.e_i])
        ) +
        # this is essentially the lie derivative and is the second to last term
        sym.Matrix([
            sum(beta[jj]*ad_(jj, Gamma_hat[ii]) for jj in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) -
        sym.Matrix([
            2 * alpha * kappa1 * sum(
                    igammat[ii, jj] * Z[jj] for jj in dendrosym.nr.e_i
                ) for ii in dendrosym.nr.e_i]) -
        16 * sym.pi * alpha * \
        sym.Matrix([sum(igammat[ii, jj] * S[jj] for jj in dendrosym.nr.e_i)
                    for ii in dendrosym.nr.e_i])
    )
    # then force this "matrix" to be a list since it's just a vector, we use the
    # matrix so that it can easily be modified/added/subtracted
    Gamma_hat_rhs = [
        item for sublist in Gamma_hat_rhs.tolist() for item in sublist
    ]
    # == END GAMMA HAT RHS ==

    # bkdkGt = sym.Matrix([
    #               sum(beta[j] * ad_(j, Gammat[i])
    #                   for j in dendrosym.nr.e_i) for i in dendrosym.nr.e_i])

    # ====================================================================
    # ====================================================================
    # ======================= GAUGE EQUATIONS ============================
    # ====================================================================

    # ==========
    # Alpha RHS
    # ===
    # alpha right hand side equation (20 in the original paper)
    # LAPSE EQUATION
    alpha_rhs = lambda1 * dendrosym.nr.lie(beta, alpha) - \
        2 * alpha * (K - 2 * Theta - K0)
    # print(alpha_rhs)
    # == END ALPHA RHS ==

    # ==========
    # beta RHS
    # ===
    # beta right hand side equation (as given by Dr. Hirschmann, but also eqn 21)
    # SHIFT EQUATION
    beta_rhs = [
        lambda2 * dendrosym.nr.vec_j_ad_j(beta, beta[ii]) +
        three_fourths * f(alpha) * B[ii] for ii in dendrosym.nr.e_i
    ]
    # NOTE: the original formualtion of beta RHS via the original code used
    # f(alpha) = (lf0 + lf1 * alpha) -> which is a linear combination of
    # parameters lambda_f1 and lambda_f2 with alpha
    # == END beta RHS ==

    # ==========
    # B RHS
    # ===
    # NOTE: B_rhs is a constraint/gauge condition but relies on another
    # formula so it is placed after the evolution eqns
    # B right hand side equation
    # This one is tricky because beta includes Gamma hat
    # (as given by Dr. Hirschmann, but also eqn (22))
    B_rhs = [(Gamma_hat_rhs[ii] - eta * B[ii] +
              lambda3 * dendrosym.nr.vec_j_ad_j(beta, B[ii]) -
              lambda4 * dendrosym.nr.vec_j_ad_j(beta, Gamma_hat[ii]))
             for ii in dendrosym.nr.e_i]
    # == END B RHS ==

    # ====================================================================
    # ====================================================================
    # SET UP THE RETURNS

    rhs_list = [
        alpha_rhs, beta_rhs, B_rhs, gammat_rhs, At_rhs, phi_rhs, K_rhs,
        Theta_rhs, Gamma_hat_rhs
    ]
    var_list = [alpha, beta, B, gammat, At, phi, K, Theta, Gamma_hat]

    return rhs_list, var_list


def physical_constraint_eqns():
    # TODO: physical constraint equations based on the psi4, ham, and momentum

    # define the coordinates
    x, y, z = sym.symbols("x, y, z")

    invsqrt2 = 0.7071067811865475244

    r_vec = sym.Matrix([[x, y, z]])
    theta = sym.Matrix([[x * z, y * z, -(x * x + y * y)]])
    phi = sym.Matrix([[-y, x, 0.0]])

    # Gram-Schmidt for orthonormal basis. THis if the non-conformally scaled metric
    gd = gammat * psi**(p_expo)

    # for r_vec
    inner_product = sum([
        sum([gd[ii, jj] * r_vec[ii] * r_vec[jj] for ii in dendrosym.nr.e_i])
        for jj in dendrosym.nr.e_i
    ])

    r_vec /= sym.sqrt(inner_product)

    # now for theta
    inner_product_1 = sum([
        sum([gd[ii, jj] * theta[ii] * theta[jj] for ii in dendrosym.nr.e_i])
        for jj in dendrosym.nr.e_i
    ])
    inner_product_2 = sum([
        sum([gd[ii, jj] * r_vec[ii] * theta[jj] for ii in dendrosym.nr.e_i])
        for jj in dendrosym.nr.e_i
    ])

    theta -= inner_product_2 * r_vec
    theta /= sym.sqrt(inner_product_1 - inner_product_2 * inner_product_2)

    # now for phi
    inner_product_1 = sum([
        sum([gd[ii, jj] * phi[ii] * phi[jj] for ii in dendrosym.nr.e_i])
        for jj in dendrosym.nr.e_i
    ])
    inner_product_2 = sum([
        sum([gd[ii, jj] * r_vec[ii] * phi[jj] for ii in dendrosym.nr.e_i])
        for jj in dendrosym.nr.e_i
    ])
    inner_product_3 = sum([
        sum([gd[ii, jj] * theta[ii] * phi[jj] for ii in dendrosym.nr.e_i])
        for jj in dendrosym.nr.e_i
    ])

    phi -= inner_product_2 * r_vec + inner_product_3 * theta
    phi /= sym.sqrt(inner_product_1 - inner_product_2 * inner_product_2 -
                    inner_product_3 * inner_product_3)

    # END TETRAD CONSTRUCTION

    # now we calculate the Weyl scalar, which is for Psi4 wave extraction

    # rename tetrad quantities for psi-4
    r_np = sym.Matrix([[r_vec[0], r_vec[1], r_vec[2]]])
    m_np_real = sym.Matrix([[theta[0], theta[1], theta[2]]]) * invsqrt2
    m_np_img = sym.Matrix([[phi[0], phi[1], phi[2]]]) * invsqrt2

    At_UD = dendrosym.nr.up_down(At)

    # Auxilary variables, MM and NN, MR and NR
    # MM and NN are 3x3 symmetric matrices
    # MR and NR are 3x3 antisymmetric matrices
    MM = sym.Matrix([
        m_np_real[ii] * m_np_real[jj] - m_np_img[ii] * m_np_img[jj]
        for ii, jj in dendrosym.nr.e_ij
    ])
    MM = MM.reshape(3, 3)

    NN = sym.Matrix([
        m_np_real[ii] * m_np_img[jj] - m_np_real[jj] * m_np_img[ii]
        for ii, jj in dendrosym.nr.e_ij
    ])
    NN = NN.reshape(3, 3)

    # additional intermediate variables
    Atr_vec = [
        sum([At[ii, jj] * r_np[jj] for jj in dendrosym.nr.e_i])
        for ii in dendrosym.nr.e_i
    ]
    Atrr = sum([Atr_vec[ii] * r_np[ii] for ii in dendrosym.nr.e_i])

    # further intermediate variables
    Uu = sym.Matrix([
        sum([
            r_np[kk] *
            (d_(kk, At[ii, jj]) -
             sum([C2[mm, kk, ii] * At[mm, jj] for mm in dendrosym.nr.e_i]))
            for kk in dendrosym.nr.e_i
        ]) for ii, jj in dendrosym.nr.e_ij
    ])
    Uu = Uu.reshape(3, 3)

    Vv = sym.Matrix([
        sum([
            r_np[kk] *
            (d_(jj, At[ii, kk]) -
             sum([C2[mm, jj, ii] * At[mm, kk] for mm in dendrosym.nr.e_i]))
            for kk in dendrosym.nr.e_i
        ]) for ii, jj in dendrosym.nr.e_ij
    ])
    Vv = Vv.reshape(3, 3)

    # r_d_psi
    r_d_psi = sum([r_np[ii] * d_(ii, psi) for ii in dendrosym.nr.e_i])

    # A temporary
    A_temp = psi**(p_expo) * (one_third * K + psi**(p_expo) * Atrr -
                              0.5 * p_expo * r_d_psi / psi)

    # now actually calculate Psi4
    Psi4_temp = sym.Matrix([
        R[ii, jj] + At[ii, jj] * A_temp - Atr_vec[ii] *
        (psi**(2 * p_expo) * Atr_vec[jj] - 0.5 * p_expo * psi**
         (p_expo - 1) * d_(jj, psi)) - psi**(p_expo) *
        (Uu[ii, jj] - Vv[ii, jj]) for ii, jj in dendrosym.nr.e_ij
    ])
    Psi4_temp = Psi4_temp.reshape(3, 3)

    psi4_real_rhs = sum(
        [Psi4_temp[ii, jj] * MM[ii, jj] for ii, jj in dendrosym.nr.e_ij])
    psi4_imag_rhs = sum(
        [Psi4_temp[ii, jj] * NN[ii, jj] for ii, jj in dendrosym.nr.e_ij])

    # ACTUAL CONSTRAINT EQUATIONS
    ham_rhs = sum(psi**(-p_expo) * igammat[jj, kk] * R[jj, kk] for jj, kk in
                  dendrosym.nr.e_ij) - dendrosym.nr.sqr(At) + two_thirds * K**2

    # momentum constraint
    mom_rhs = (sym.Matrix([
        sum([
            igammat[jj, kk] *
            (d_(jj, At[kk, ii]) - sum(C2[mm, jj, ii] * At[kk, mm]
                                      for mm in dendrosym.nr.e_i))
            for jj, kk in dendrosym.nr.e_ij
        ]) for ii in dendrosym.nr.e_i
    ]) - sym.Matrix([
        sum([CalGt[jj] * At[ii, jj] for jj in dendrosym.nr.e_i])
        for ii in dendrosym.nr.e_i
    ]) + (1.5) * p_expo * sym.Matrix([
        sum([
            igammat[jj, kk] * At[kk, ii] * d_(jj, psi) / psi
            for jj, kk in dendrosym.nr.e_ij
        ]) for ii in dendrosym.nr.e_i
    ]) - two_thirds * sym.Matrix([d_(ii, K) for ii in dendrosym.nr.e_i]))
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
gammat_f_and_a = [[[1.0, 1.0], [1.0, 0.0], [1.0, 0.0]],
                  [[], [1.0, 1.0], [1.0, 0.0]], [[], [], [1.0, 1.0]]]
At_f_and_a = [2.0, 0.0
              ]  # At is a matrix, but you can also have them all be the same
phi_f_and_a = [1.0, 1.0]
K_f_and_a = [1.0, 0.0]
Theta_f_and_a = [1.0, 0.0]
Gamma_hat_f_and_a = [[2.0, 0.0], [2.0, 0.0], [2.0, 0.0]]
alpha_f_and_a = [1.0, 1.0]
beta_f_and_a = [1.0, 1.0]  # this is a 3 vec, but trying with them all the same
B_f_and_a = [1.0,
             0.0]  # this is also a 3 vec, but trying with them all the same

dendroConfigs.set_bhs_falloff_and_asymptotic(
    "evolution", [gammat, At, phi, K, Theta, Gamma_hat, alpha, beta, B], [
        gammat_f_and_a, At_f_and_a, phi_f_and_a, K_f_and_a, Theta_f_and_a,
        Gamma_hat_f_and_a, alpha_f_and_a, beta_f_and_a, B_f_and_a
    ])

# add a few evolution constraints
dendroConfigs.add_evolution_constraint(At, "trace_zero")
dendroConfigs.add_evolution_constraint(phi, "pos_floor")
dendroConfigs.add_evolution_constraint(alpha, "pos_floor")
