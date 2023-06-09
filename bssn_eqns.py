import dendrosym
import sympy as sym

# GLOBAL PARAMETERS
# which eta to use, either eta or eta_func
ETA_SELECTION = "eta"
IS_STAGED = False

# shortcuts for using these specific equations
pi = sym.pi
exp = sym.exp
pow = sym.Pow
sqrt = sym.sqrt
Rational = sym.Rational

# DEFINE: dendro config class to use for generating code
dendroConfigs = dendrosym.NRConfig("bssn")

# the indexing used for the dendro configs
idx_str = "[pp]"

# save the index string
dendroConfigs.set_idx_str(idx_str)

# BSSN parameter variables
lambda_param = dendrosym.dtypes.ParameterVariable(
    "lambda", dtype="unsigned int", num_params=4
)
lambda1, lambda2, lambda3, lambda4 = lambda_param.get_symbolic_repr()

l1 = lambda1
l2 = lambda2
l3 = lambda3
l4 = lambda4

# sets the lambda f parameters which are also involved with the gauge.
lf_param = dendrosym.dtypes.ParameterVariable(
    "lf", dtype="double", num_params=2
)
lf0, lf1 = lf_param.get_symbolic_repr()

# sets the eta parameters which come in with the constraint damping
eta_param = dendrosym.dtypes.ParameterVariable(
    "eta", dtype="double", num_params=1
)
eta = eta_param.get_symbolic_repr()

# damping parameters
R0_param = dendrosym.dtypes.ParameterVariable(
    "damping_eta_r0", dtype="double", num_params=1
)
R0 = R0_param.get_symbolic_repr()

eta_power_param = dendrosym.dtypes.ParameterVariable(
    "damping_eta_power", dtype="double", num_params=2
)
ep1, ep2 = eta_power_param.get_symbolic_repr()

xi_param = dendrosym.dtypes.ParameterVariable(
    "damping_xi", dtype="unsigned int", num_params=3
)
xi1, xi2, xi3 = xi_param.get_symbolic_repr()

# add the parameters to evolution
dendroConfigs.add_parameter_variables(
    [lambda_param, lf_param, eta_param, R0_param, eta_power_param, xi_param],
    "evolution",
)
# TODO: add parameters to constraints
#

# ==================
# CONSTRAINT VARIABLES
ham_constraint = dendrosym.dtypes.scalar("ham" + idx_str)
mom_constraint = dendrosym.dtypes.vec3("mom" + idx_str)
psi4_real_constraint = dendrosym.dtypes.scalar("psi4_real" + idx_str)
psi4_imag_constraint = dendrosym.dtypes.scalar("psi4_imag" + idx_str)
dendroConfigs.add_constraint_variables(
    [
        ham_constraint,
        mom_constraint,
        psi4_real_constraint,
        psi4_imag_constraint,
    ]
)

# =================
# EVOLUTION VARIABLES

a = dendrosym.dtypes.scalar("alpha" + idx_str)
chi = dendrosym.dtypes.scalar("chi" + idx_str)
K = dendrosym.dtypes.scalar("K" + idx_str)

Gt = dendrosym.dtypes.vec3("Gt" + idx_str)
b = dendrosym.dtypes.vec3("beta" + idx_str)
B = dendrosym.dtypes.vec3("B" + idx_str)

gt = dendrosym.dtypes.sym_3x3("gt" + idx_str)
At = dendrosym.dtypes.sym_3x3("At" + idx_str)

# Gt_rhs  = dendrosym.dtypes.vec3("Gt_rhs" + idx_str)

dendroConfigs.add_evolution_variables([a, chi, K, Gt, b, B, gt, At])

dendroConfigs.set_advective_derivative_var(b)

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
ad_ = dendrosym.nr.set_advective_derivative("grad")
# and then we set the kreiss oliger dissipation
kod_ = dendrosym.nr.set_kreiss_oliger_dissipation("kograd")
# == END DERIVATIVES ==

# Lie derivative weight
weight = -sym.Rational(2, 3)
weight_Gt = sym.Rational(2, 3)

dendroConfigs.set_metric(gt)
# and then we get the inverse (conformal) metric
igt = dendrosym.nr.get_inverse_metric()
inv_metric = igt
# as well as the two Christoffel symbols built from the conformal metric
C1 = dendrosym.nr.get_first_christoffel()
C2 = dendrosym.nr.get_second_christoffel()
# and the third Christoffel symbol
C3 = dendrosym.nr.get_complete_christoffel(chi)

C2_spatial = C3  # the spatial C2 is the complete christoffel
R, Rt, Rphi, CalGt = dendrosym.nr.compute_ricci(Gt, chi)

# define eta_func
eta_func = (
    R0
    * sqrt(
        sum(
            [igt[i, j] * d_(i, chi) * d_(j, chi) for i, j in dendrosym.nr.e_ij]
        )
    )
    / ((1 - chi**ep1) ** ep2)
)

"""
BSSN puncture gauge (HAD/ traditional BSSN puncture gaugue) with const eta damping 
"""

if ETA_SELECTION == "eta":
    eta_damp = eta
elif ETA_SELECTION == "eta_func":
    eta_damp = eta_func


def bssn_puncture_gauge():

    if not IS_STAGED:

        a_rhs = l1 * dendrosym.nr.lie(b, a) - 2 * a * K

        b_rhs = [
            (
                Rational(3, 4) * (lf0 + lf1 * a) * B[i]
                + l2 * dendrosym.nr.vec_j_ad_j(b, b[i])
            )
            for i in dendrosym.nr.e_i
        ]

        gt_rhs = dendrosym.nr.lie(b, gt, weight) - 2 * a * At

        chi_rhs = dendrosym.nr.lie(b, chi, weight) + Rational(2, 3) * (
            chi * a * K
        )

        AikAkj = sym.Matrix(
            [
                sum(
                    [
                        At[i, k]
                        * sum(
                            [
                                dendrosym.nr.inv_metric[k, l] * At[l, j]
                                for l in dendrosym.nr.e_i
                            ]
                        )
                        for k in dendrosym.nr.e_i
                    ]
                )
                for i, j in dendrosym.nr.e_ij
            ]
        )

        At_rhs = (
            dendrosym.nr.lie(b, At, weight)
            + chi * dendrosym.nr.trace_free(a * R - dendrosym.nr.DiDj(a))
            + a * (K * At - 2 * AikAkj.reshape(3, 3))
        )

        K_rhs = (
            dendrosym.nr.lie(b, K)
            - dendrosym.nr.laplacian(a, chi)
            + a * (K * K / 3 + dendrosym.nr.sqr(At))
        )

        At_UU = dendrosym.nr.up_up(At)

        Gt_rhs = (
            sym.Matrix(
                [
                    sum(b[j] * ad_(j, Gt[i]) for j in dendrosym.nr.e_i)
                    for i in dendrosym.nr.e_i
                ]
            )
            - sym.Matrix(
                [
                    sum(CalGt[j] * d_(j, b[i]) for j in dendrosym.nr.e_i)
                    for i in dendrosym.nr.e_i
                ]
            )
            + Rational(2, 3)
            * sym.Matrix(
                [
                    CalGt[i] * sum(d_(j, b[j]) for j in dendrosym.nr.e_i)
                    for i in dendrosym.nr.e_i
                ]
            )
            + sym.Matrix(
                [
                    sum(
                        [
                            igt[j, k] * d2_(j, k, b[i])
                            + igt[i, j] * d2_(j, k, b[k]) / 3
                            for j, k in dendrosym.nr.e_ij
                        ]
                    )
                    for i in dendrosym.nr.e_i
                ]
            )
            - sym.Matrix(
                [
                    sum([2 * At_UU[i, j] * d_(j, a) for j in dendrosym.nr.e_i])
                    for i in dendrosym.nr.e_i
                ]
            )
            + sym.Matrix(
                [
                    sum(
                        [
                            2 * a * dendrosym.nr.C2[i, j, k] * At_UU[j, k]
                            for j, k in dendrosym.nr.e_ij
                        ]
                    )
                    for i in dendrosym.nr.e_i
                ]
            )
            - sym.Matrix(
                [
                    sum(
                        [
                            a
                            * (
                                3 / chi * At_UU[i, j] * d_(j, chi)
                                + Rational(4, 3) * inv_metric[i, j] * d_(j, K)
                            )
                            for j in dendrosym.nr.e_i
                        ]
                    )
                    for i in dendrosym.nr.e_i
                ]
            )
        )

        Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        B_rhs = [
            (
                Gt_rhs[i]
                - eta_damp * B[i]
                + l3 * dendrosym.nr.vec_j_ad_j(b, B[i])
                - l4 * dendrosym.nr.vec_j_ad_j(b, Gt[i])
            )
            for i in dendrosym.nr.e_i
        ]

        ###################################################################
        # generate code
        ###################################################################

        rhs_list = [
            a_rhs,
            b_rhs,
            gt_rhs,
            chi_rhs,
            At_rhs,
            K_rhs,
            Gt_rhs,
            B_rhs,
        ]
        vnames = [
            "a_rhs",
            "b_rhs",
            "gt_rhs",
            "chi_rhs",
            "At_rhs",
            "K_rhs",
            "Gt_rhs",
            "B_rhs",
        ]
        var_list = [a, b, gt, chi, At, K, Gt, B]

        return rhs_list, var_list

    else:
        pass
        # TODO: DendroSym does not yet support staged generation when using all of the internal pieces
        # note: these are just the symbolic vars that is being used to generate the
        # Gt_rhs by satges

        # _Gt_rhs_s1 = dendro.vec3("Gt_rhs_s1_", "[pp]")
        # _Gt_rhs_s2 = dendro.vec3("Gt_rhs_s2_", "[pp]")
        # _Gt_rhs_s3 = dendro.vec3("Gt_rhs_s3_", "[pp]")
        # _Gt_rhs_s4 = dendro.vec3("Gt_rhs_s4_", "[pp]")
        # _Gt_rhs_s5 = dendro.vec3("Gt_rhs_s5_", "[pp]")
        # _Gt_rhs_s6 = dendro.vec3("Gt_rhs_s6_", "[pp]")
        # _Gt_rhs_s7 = dendro.vec3("Gt_rhs_s7_", "[pp]")
        # _CalGt = dendro.vec3("CalGt", "[pp]")
        # _Gt_rhs = dendro.vec3("Gt_rhs", "[pp]")

        # # Gt_rhs staged vars that is being used to generate the code.
        # At_UU = dendro.sym_3x3("At_UU", "[pp]")
        # CalGt = dendro.vec3("CalGt", "[pp]")
        # Gt_rhs_s1 = dendro.vec3("Gt_rhs_s1_", "[pp]")
        # Gt_rhs_s2 = dendro.vec3("Gt_rhs_s2_", "[pp]")
        # Gt_rhs_s3 = dendro.vec3("Gt_rhs_s3_", "[pp]")
        # Gt_rhs_s4 = dendro.vec3("Gt_rhs_s4_", "[pp]")
        # Gt_rhs_s5 = dendro.vec3("Gt_rhs_s5_", "[pp]")
        # Gt_rhs_s6 = dendro.vec3("Gt_rhs_s6_", "[pp]")
        # Gt_rhs_s7 = dendro.vec3("Gt_rhs_s7_", "[pp]")

        # C1 = dendro.get_first_christoffel()
        # C2 = dendro.get_second_christoffel()
        # C2_spatial = dendro.get_complete_christoffel(chi)
        # [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

        # a_rhs = l1 * dendro.lie(b, a) - 2 * a * K

        # b_rhs = [
        #     (Rational(3, 4) * (lf0 + lf1 * a) * B[i] + l2 * dendro.vec_j_ad_j(b, b[i]))
        #     for i in dendro.e_i
        # ]

        # gt_rhs = dendro.lie(b, gt, weight) - 2 * a * At

        # chi_rhs = dendro.lie(b, chi, weight) + Rational(2, 3) * (chi * a * K)

        # AikAkj = Matrix(
        #     [
        #         sum(
        #             [
        #                 At[i, k]
        #                 * sum([inv_metric[k, l] * At[l, j] for l in dendro.e_i])
        #                 for k in dendro.e_i
        #             ]
        #         )
        #         for i, j in dendro.e_ij
        #     ]
        # )

        # At_rhs = (
        #     dendro.lie(b, At, weight)
        #     + chi * dendro.trace_free(a * R - dendro.DiDj(a))
        #     + a * (K * At - 2 * AikAkj.reshape(3, 3))
        # )

        # K_rhs = (
        #     dendro.lie(b, K)
        #     - dendro.laplacian(a, chi)
        #     + a * (K * K / 3 + dendro.sqr(At))
        # )

        # At_UU = dendro.up_up(At)

        # Gt_rhs_s1 = [sum(b[j] * ad(j, Gt[i]) for j in dendro.e_i) for i in dendro.e_i]
        # Gt_rhs_s2 = [
        #     sum(_CalGt[j] * d(j, b[i]) for j in dendro.e_i) for i in dendro.e_i
        # ]
        # Gt_rhs_s3 = [
        #     _CalGt[i] * sum(d(j, b[j]) for j in dendro.e_i) for i in dendro.e_i
        # ]
        # Gt_rhs_s4 = [
        #     sum(
        #         [
        #             igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k]) / 3
        #             for j, k in dendro.e_ij
        #         ]
        #     )
        #     for i in dendro.e_i
        # ]
        # Gt_rhs_s5 = [
        #     sum([2 * At_UU[i, j] * d(j, a) for j in dendro.e_i]) for i in dendro.e_i
        # ]
        # Gt_rhs_s6 = [
        #     sum([2 * a * dendro.C2[i, j, k] * At_UU[j, k] for j, k in dendro.e_ij])
        #     for i in dendro.e_i
        # ]
        # Gt_rhs_s7 = [
        #     sum(
        #         [
        #             a
        #             * (
        #                 3 / chi * At_UU[i, j] * d(j, chi)
        #                 + Rational(4, 3) * inv_metric[i, j] * d(j, K)
        #             )
        #             for j in dendro.e_i
        #         ]
        #     )
        #     for i in dendro.e_i
        # ]

        # Gt_rhs = (
        #     Matrix(_Gt_rhs_s1)
        #     - Matrix(_Gt_rhs_s2)
        #     + Rational(2, 3) * Matrix(_Gt_rhs_s3)
        #     + Matrix(_Gt_rhs_s4)
        #     - Matrix(_Gt_rhs_s5)
        #     + Matrix(_Gt_rhs_s6)
        #     - Matrix(_Gt_rhs_s7)
        # )

        # Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        # B_rhs = [
        #     (
        #         Gt_rhs[i]
        #         - eta_damp * B[i]
        #         + l3 * dendro.vec_j_ad_j(b, B[i])
        #         - l4 * dendro.vec_j_ad_j(b, Gt[i])
        #     )
        #     for i in dendro.e_i
        # ]

        # outs = [
        #     a_rhs,
        #     b_rhs,
        #     gt_rhs,
        #     chi_rhs,
        #     At_rhs,
        #     K_rhs,
        #     CalGt,
        #     Gt_rhs_s1,
        #     Gt_rhs_s2,
        #     Gt_rhs_s3,
        #     Gt_rhs_s4,
        #     Gt_rhs_s5,
        #     Gt_rhs_s6,
        #     Gt_rhs_s7,
        #     Gt_rhs,
        #     B_rhs,
        # ]
        # vnames = [
        #     "a_rhs",
        #     "b_rhs",
        #     "gt_rhs",
        #     "chi_rhs",
        #     "At_rhs",
        #     "K_rhs",
        #     "CalGt",
        #     "Gt_rhs_s1_",
        #     "Gt_rhs_s2_",
        #     "Gt_rhs_s3_",
        #     "Gt_rhs_s4_",
        #     "Gt_rhs_s5_",
        #     "Gt_rhs_s6_",
        #     "Gt_rhs_s7_",
        #     "Gt_rhs",
        #     "B_rhs",
        # ]


dendroConfigs.set_rhs_equation_function("evolution", bssn_puncture_gauge)


def bssn_constraint_func():

    # Define coordinates
    x, y, z = sym.symbols("x, y, z")

    # Some other values
    invsqrt2 = 1 / sqrt(2)  # 0.7071067811865475244
    inv_chi = 1 / chi

    # Define the original spatial vectors in our tetrad
    r_vec = sym.Matrix([[x, y, z]])
    theta = sym.Matrix([[x * z, y * z, -(x * x + y * y)]])
    phi = sym.Matrix([[-y, x, 0.0]])

    # We use Gram-Schmidt to make the basis orthonormal.
    # Note that we use the original (not conformally rescaled) metric to define
    # the tetrad and correspondingly Psi4.
    gd = gt * inv_chi

    # For r_vec
    inner_product = 0.0
    inner_product = sum(
        [
            sum([gd[i, j] * r_vec[i] * r_vec[j] for i in dendrosym.nr.e_i])
            for j in dendrosym.nr.e_i
        ]
    )

    r_vec /= sqrt(inner_product)

    # For theta
    inner_product_1 = 0.0
    inner_product_2 = 0.0

    inner_product_1 = sum(
        [
            sum([gd[i, j] * theta[i] * theta[j] for i in dendrosym.nr.e_i])
            for j in dendrosym.nr.e_i
        ]
    )
    inner_product_2 = sum(
        [
            sum([gd[i, j] * r_vec[i] * theta[j] for i in dendrosym.nr.e_i])
            for j in dendrosym.nr.e_i
        ]
    )

    theta -= inner_product_2 * r_vec
    theta /= sqrt(inner_product_1 - inner_product_2 * inner_product_2)

    # For phi
    inner_product_1 = 0.0
    inner_product_2 = 0.0
    inner_product_3 = 0.0

    inner_product_1 = sum(
        [
            sum([gd[i, j] * phi[i] * phi[j] for i in dendrosym.nr.e_i])
            for j in dendrosym.nr.e_i
        ]
    )
    inner_product_2 = sum(
        [
            sum([gd[i, j] * r_vec[i] * phi[j] for i in dendrosym.nr.e_i])
            for j in dendrosym.nr.e_i
        ]
    )
    inner_product_3 = sum(
        [
            sum([gd[i, j] * theta[i] * phi[j] for i in dendrosym.nr.e_i])
            for j in dendrosym.nr.e_i
        ]
    )

    phi -= inner_product_2 * r_vec + inner_product_3 * theta
    phi /= sqrt(
        inner_product_1
        - inner_product_2 * inner_product_2
        - inner_product_3 * inner_product_3
    )

    # This completes the tetrad construction.
    ###################################################################
    # Calculate the Weyl scalar, Psi4, for graviational wave extraction
    ###################################################################

    # Rename the tetrad quantities for calculating Psi4
    r_np = sym.Matrix([[r_vec[0], r_vec[1], r_vec[2]]])
    m_np_real = sym.Matrix([[theta[0], theta[1], theta[2]]]) * invsqrt2
    m_np_img = sym.Matrix([[phi[0], phi[1], phi[2]]]) * invsqrt2

    # Some auxilary variables
    # MM and NN are symmetric 2nd rank objects and
    # MR and NR are anti-symmetric 2nd rank objects

    MM = sym.Matrix(
        [
            m_np_real[i] * m_np_real[j] - m_np_img[i] * m_np_img[j]
            for i, j in dendrosym.nr.e_ij
        ]
    )
    MM = MM.reshape(3, 3)
    NN = sym.Matrix(
        [
            m_np_real[i] * m_np_img[j] + m_np_real[j] * m_np_img[i]
            for i, j in dendrosym.nr.e_ij
        ]
    )
    NN = NN.reshape(3, 3)
    MR = sym.Matrix(
        [
            m_np_real[i] * r_np[j] - m_np_real[j] * r_np[i]
            for i, j in dendrosym.nr.e_ij
        ]
    )
    MR = MR.reshape(3, 3)
    NR = sym.Matrix(
        [
            m_np_img[i] * r_np[j] - m_np_img[j] * r_np[i]
            for i, j in dendrosym.nr.e_ij
        ]
    )
    NR = NR.reshape(3, 3)

    # Additional intermediate variables
    # A_vec = Matrix([[sum([At[j,0]*r_np[j] for j in dendro.e_i]), sum([At[j,1]*r_np[j] for j in dendro.e_i]),sum([At[j,2]*r_np[j] for j in dendro.e_i])]])
    A_vec = [
        sum([At[i, j] * r_np[j] for j in dendrosym.nr.e_i])
        for i in dendrosym.nr.e_i
    ]

    Uu = sym.Matrix(
        [
            sum(
                [
                    m_np_real[k]
                    * (
                        d_(j, At[k, i])
                        + sum(
                            [C2[m, k, i] * At[m, j] for m in dendrosym.nr.e_i]
                        )
                    )
                    for k in dendrosym.nr.e_i
                ]
            )
            for i, j in dendrosym.nr.e_ij
        ]
    )
    Uu = Uu.reshape(3, 3)
    Vv = sym.Matrix(
        [
            sum(
                [
                    m_np_img[k]
                    * (
                        d_(j, At[k, i])
                        + sum(
                            [C2[m, k, i] * At[m, j] for m in dendrosym.nr.e_i]
                        )
                    )
                    for k in dendrosym.nr.e_i
                ]
            )
            for i, j in dendrosym.nr.e_ij
        ]
    )
    Vv = Vv.reshape(3, 3)

    r_d_chi = sum([r_np[i] * d_(i, chi) for i in dendrosym.nr.e_i])

    A_temp = (
        inv_chi
        * inv_chi
        * (
            sum([A_vec[i] * r_np[i] for i in dendrosym.nr.e_i])
            + K * chi * Rational(1, 3)
            + Rational(1, 2) * r_d_chi
        )
    )

    m_real_d_chi = sum([m_np_real[i] * d_(i, chi) for i in dendrosym.nr.e_i])
    m_img_d_chi = sum([m_np_img[i] * d_(i, chi) for i in dendrosym.nr.e_i])

    m_real_A_vec = sum([m_np_real[i] * A_vec[i] for i in dendrosym.nr.e_i])
    m_img_A_vec = sum([m_np_img[i] * A_vec[i] for i in dendrosym.nr.e_i])

    # Calculate Psi4

    psi4_1_real = sum(
        [R[i, i] * MM[i, i] for i in dendrosym.nr.e_i]
    ) + 2 * sum(
        [
            sum([R[i, j] * MM[i, j] for j in range(i + 1, 3)])
            for i in range(0, 2)
        ]
    )
    psi4_1_img = sum([R[i, i] * NN[i, i] for i in dendrosym.nr.e_i]) + 2 * sum(
        [
            sum([R[i, j] * NN[i, j] for j in range(i + 1, 3)])
            for i in range(0, 2)
        ]
    )

    psi4_2_real = A_temp * (
        sum([At[i, i] * MM[i, i] for i in dendrosym.nr.e_i])
        + 2
        * sum(
            [
                sum([At[i, j] * MM[i, j] for j in range(i + 1, 3)])
                for i in range(0, 2)
            ]
        )
    )
    psi4_2_img = A_temp * (
        sum([At[i, i] * NN[i, i] for i in dendrosym.nr.e_i])
        + 2
        * sum(
            [
                sum([At[i, j] * NN[i, j] for j in range(i + 1, 3)])
                for i in range(0, 2)
            ]
        )
    )

    psi4_3_real = inv_chi * sum(
        [
            sum(
                [
                    MR[i, j] * Uu[i, j] - NR[i, j] * Vv[i, j]
                    for i in dendrosym.nr.e_i
                ]
            )
            for j in dendrosym.nr.e_i
        ]
    )
    psi4_3_img = inv_chi * sum(
        [
            sum(
                [
                    NR[i, j] * Uu[i, j] + MR[i, j] * Vv[i, j]
                    for i in dendrosym.nr.e_i
                ]
            )
            for j in dendrosym.nr.e_i
        ]
    )

    # 12/31/2020: 0.5 is replaced with rational.
    psi4_4_real = (
        inv_chi
        * inv_chi
        * (
            m_real_A_vec * (m_real_A_vec + Rational(1, 2) * m_real_d_chi)
            - m_img_A_vec * (m_img_A_vec + Rational(1, 2) * m_img_d_chi)
        )
    )
    psi4_4_img = (
        inv_chi
        * inv_chi
        * (
            m_real_A_vec * (m_img_A_vec - Rational(1, 2) * m_img_d_chi)
            + m_img_A_vec * (m_real_A_vec - Rational(1, 2) * m_real_d_chi)
        )
    )

    # Adding previous auxilary Psi4 calculations
    # 12/31/2020 : There is a - sign convention issue to match the sign with the LazEv Code.
    # psi4_real =     psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real
    # psi4_img  = - ( psi4_1_img  + psi4_2_img  - psi4_3_img  - psi4_4_img  )
    psi4_real = -(psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real)
    psi4_img = psi4_1_img + psi4_2_img - psi4_3_img - psi4_4_img
    ###################################################################
    # Constraint Equations
    ###################################################################

    # The Hamiltonian constraint
    ham = (
        sum(chi * igt[j, k] * R[j, k] for j, k in dendrosym.nr.e_ij)
        - dendrosym.nr.sqr(At)
        + Rational(2, 3) * K**2
    )

    # The momentum  constraints
    mom = (
        sym.Matrix(
            [
                sum(
                    [
                        igt[j, k]
                        * (
                            d_(k, At[i, j])
                            - sum(
                                C2[m, k, i] * At[j, m]
                                for m in dendrosym.nr.e_i
                            )
                        )
                        for j, k in dendrosym.nr.e_ij
                    ]
                )
                for i in dendrosym.nr.e_i
            ]
        )
        - sym.Matrix(
            [
                sum([Gt[j] * At[i, j] for j in dendrosym.nr.e_i])
                for i in dendrosym.nr.e_i
            ]
        )
        - Rational(3, 2)
        * sym.Matrix(
            [
                sum(
                    [
                        igt[j, k] * At[k, i] * d_(j, chi) / chi
                        for j, k in dendrosym.nr.e_ij
                    ]
                )
                for i in dendrosym.nr.e_i
            ]
        )
        - Rational(2, 3) * sym.Matrix([d_(i, K) for i in dendrosym.nr.e_i])
    )
    mom = [item for sublist in mom.tolist() for item in sublist]

    rhs_list = [psi4_real, psi4_img, ham, mom]
    var_list = [
        psi4_real_constraint,
        psi4_imag_constraint,
        ham_constraint,
        mom_constraint,
    ]

    return rhs_list, var_list


dendroConfigs.set_rhs_equation_function("constraint", bssn_constraint_func)


# BOUNDARY INFORMATION

a_f_and_a = [1.0, 1.0]

chi_f_and_a = [1.0, 1.0]

K_f_and_a = [1.0, 1.0]

b_f_and_a = [[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]]

Gt_f_and_a = [[2.0, 0.0], [2.0, 0.0], [2.0, 0.0]]

B_f_and_a = [[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]]

At_f_and_a = [
    [[2.0, 0.0], [2.0, 0.0], [2.0, 0.0]],
    [[2.0, 0.0], [2.0, 0.0], [2.0, 0.0]],
    [[2.0, 0.0], [2.0, 0.0], [2.0, 0.0]],
]

gt_f_and_a = [
    [[1.0, 1.0], [1.0, 0.0], [1.0, 0.0]],
    [[1.0, 0.0], [1.0, 1.0], [1.0, 0.0]],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
]

dendroConfigs.set_bhs_falloff_and_asymptotic(
    "evolution",
    [a, chi, K, b, Gt, B, At, gt],
    [
        a_f_and_a,
        chi_f_and_a,
        K_f_and_a,
        b_f_and_a,
        Gt_f_and_a,
        B_f_and_a,
        At_f_and_a,
        gt_f_and_a,
    ],
)

# add a few constraints, like At needs to have a trace of zero, chi needs a positive floor, and a needs a positive floor as well
dendroConfigs.add_evolution_constraint(At, "trace_zero")
dendroConfigs.add_evolution_constraint(chi, "pos_floor")
dendroConfigs.add_evolution_constraint(a, "pos_floor")

dendroConfigs.replace_derivatives_with_stencil("evolution", 6)

dendroConfigs.replace_derivatives_with_stencil("constraint", 6)

with open("output.cpp", "w") as f:
    generated_code = dendroConfigs.generate_rhs_code("evolution")
    f.write(generated_code)

with open("output_constraint.cpp", "w") as f:
    generated_code = dendroConfigs.generate_rhs_code("constraint")
    f.write(generated_code)

# # TODO: REMOVE THIS BEFORE GENERATING THE REAL C++ CODE
evolution_var_extraction = dendroConfigs.generate_variable_extraction(
    "evolution", use_const=True
)
evolution_var_rhs_extraction = dendroConfigs.generate_rhs_var_extraction(
    "evolution", zip_var_name="unzipVarsRHS"
)

print("// EVOLUTION VARIABLE EXTRACTION CODE")
print(evolution_var_extraction)
print(evolution_var_rhs_extraction)

evolution_parameters = dendroConfigs.gen_parameter_code("evolution")
print()
print("// PARAMETER EXTRACTION")
print(evolution_parameters)

# we know that BSSN doesn't have "nested" derivatives
dendroConfigs.override_derivative_expansion()


(
    intermediate_grad_str,
    deallocate_intermediate_grad_str,
) = dendroConfigs.generate_pre_necessary_derivatives(
    "evolution", dtype="double", include_byte_declaration=False
)

print()
# print(intermediate_grad_str)

(
    deriv_alloc,
    deriv_calc,
    deriv_dealloc,
) = dendroConfigs.generate_deriv_allocation_and_calc(
    "evolution", include_byte_declaration=False
)

evolution_rhs_code = dendroConfigs.generate_rhs_code("evolution")

with open("temporary_rhs_output.cpp", "w") as f:
    f.write(evolution_rhs_code)
