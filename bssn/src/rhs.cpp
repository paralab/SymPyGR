#include "rhs.h"


using namespace std;
using namespace bssn;

/*----------------------------------------------------------------------;
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void bssnrhs(double **unzipVarsRHS, const double **uZipVars,
             const unsigned int& offset,
             const double *pmin, const double *pmax, const unsigned int *sz,
             const unsigned int& bflag)
{



    const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
    const double *chi = &uZipVars[VAR::U_CHI][offset];
    const double *K = &uZipVars[VAR::U_K][offset];
    const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
    const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
    const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
    const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
    const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
    const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
    const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
    const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
    const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
    const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
    const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
    const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
    const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
    const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
    const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
    const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
    const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
    const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
    const double *B0 = &uZipVars[VAR::U_B0][offset];
    const double *B1 = &uZipVars[VAR::U_B1][offset];
    const double *B2 = &uZipVars[VAR::U_B2][offset];

    double *a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
    double *chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
    double *K_rhs = &unzipVarsRHS[VAR::U_K][offset];
    double *gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
    double *gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
    double *gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
    double *gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
    double *gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
    double *gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
    double *b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
    double *b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
    double *b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
    double *At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
    double *At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
    double *At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
    double *At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
    double *At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
    double *At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
    double *Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
    double *Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
    double *Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
    double *B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
    double *B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
    double *B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx = (pmax[0] - pmin[0]) / (nx - 1);
    double hy = (pmax[1] - pmin[1]) / (ny - 1);
    double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int lambda[4] = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                    BSSN_LAMBDA[2], BSSN_LAMBDA[3]
                                   };
    const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};



    int idx[3];

    unsigned int n = sz[0]*sz[1]*sz[2];

    bssn::timer::t_deriv.start();


    /*[[[cog
    import cog
    import bssnDerivs as bssnDerivs

    for deriv in bssnDerivs.FUNC_D_I:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")

    for deriv in bssnDerivs.FUNC_D_IJ:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")

    for deriv in bssnDerivs.FUNC_AD_I:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")


    for var in bssnDerivs.D:
        cog.outl("\t deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[0] + var ,var))
        cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[1] + var ,var))
        cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[2] + var ,var))

        if var in bssnDerivs.DD:
            cog.outl("\t deriv_xx(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[0] + var ,var))
            cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[1] + var , bssnDerivs.PREFIX_D[0] + var ))
            cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[2] + var , bssnDerivs.PREFIX_D[0] + var ))

            cog.outl("\t deriv_yy(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[3] + var ,var))
            cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[4] + var , bssnDerivs.PREFIX_D[1] + var))

            cog.outl("\t deriv_zz(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[5] + var ,var))

        if var in bssnDerivs.AD:
            cog.outl("\t adv_deriv_x(%s, %s, hx, sz, beta0, bflag);" %(bssnDerivs.PREFIX_AD[0] + var ,var))
            cog.outl("\t adv_deriv_y(%s, %s, hx, sz, beta1, bflag);" %(bssnDerivs.PREFIX_AD[1] + var ,var))
            cog.outl("\t adv_deriv_z(%s, %s, hx, sz, beta2, bflag);" %(bssnDerivs.PREFIX_AD[2] + var ,var))



    ]]]*/
    double* grad_0_alpha = (double*)malloc(sizeof(double)*n);
    double* grad_1_alpha = (double*)malloc(sizeof(double)*n);
    double* grad_2_alpha = (double*)malloc(sizeof(double)*n);
    double* grad_0_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_B0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_B0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_B0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_B1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_B1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_B1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_B2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_B2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_B2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_chi = (double*)malloc(sizeof(double)*n);
    double* grad_1_chi = (double*)malloc(sizeof(double)*n);
    double* grad_2_chi = (double*)malloc(sizeof(double)*n);
    double* grad_0_Gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_Gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_Gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_Gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_Gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_Gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_Gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_Gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_Gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_K = (double*)malloc(sizeof(double)*n);
    double* grad_1_K = (double*)malloc(sizeof(double)*n);
    double* grad_2_K = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At3 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At3 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At3 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At4 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At4 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At4 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At5 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At5 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At5 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt3 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt3 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt4 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt4 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt5 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt5 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At3 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At3 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At3 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At4 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At4 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At4 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At5 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At5 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At5 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_alpha = (double*)malloc(sizeof(double)*n);
    double* agrad_1_alpha = (double*)malloc(sizeof(double)*n);
    double* agrad_2_alpha = (double*)malloc(sizeof(double)*n);
    double* agrad_0_beta0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_beta0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_beta1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_beta1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_beta2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_beta2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_chi = (double*)malloc(sizeof(double)*n);
    double* agrad_1_chi = (double*)malloc(sizeof(double)*n);
    double* agrad_2_chi = (double*)malloc(sizeof(double)*n);
    double* agrad_0_Gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_Gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_Gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_Gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_Gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_Gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_Gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_Gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_Gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_K = (double*)malloc(sizeof(double)*n);
    double* agrad_1_K = (double*)malloc(sizeof(double)*n);
    double* agrad_2_K = (double*)malloc(sizeof(double)*n);
    double* agrad_0_B0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_B0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_B0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_B1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_B1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_B1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_B2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_B2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_B2 = (double*)malloc(sizeof(double)*n);
    deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
    deriv_y(grad_1_alpha, alpha, hx, sz, bflag);
    deriv_z(grad_2_alpha, alpha, hx, sz, bflag);
    deriv_xx(grad2_0_0_alpha, alpha, hx, sz, bflag);
    deriv_y(grad2_0_1_alpha, grad_0_alpha, hx, sz, bflag);
    deriv_z(grad2_0_2_alpha, grad_0_alpha, hx, sz, bflag);
    deriv_yy(grad2_1_1_alpha, alpha, hx, sz, bflag);
    deriv_z(grad2_1_2_alpha, grad_1_alpha, hx, sz, bflag);
    deriv_zz(grad2_2_2_alpha, alpha, hx, sz, bflag);
    adv_deriv_x(agrad_0_alpha, alpha, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_alpha, alpha, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_alpha, alpha, hx, sz, beta2, bflag);
    deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
    deriv_y(grad_1_beta0, beta0, hx, sz, bflag);
    deriv_z(grad_2_beta0, beta0, hx, sz, bflag);
    deriv_xx(grad2_0_0_beta0, beta0, hx, sz, bflag);
    deriv_y(grad2_0_1_beta0, grad_0_beta0, hx, sz, bflag);
    deriv_z(grad2_0_2_beta0, grad_0_beta0, hx, sz, bflag);
    deriv_yy(grad2_1_1_beta0, beta0, hx, sz, bflag);
    deriv_z(grad2_1_2_beta0, grad_1_beta0, hx, sz, bflag);
    deriv_zz(grad2_2_2_beta0, beta0, hx, sz, bflag);
    adv_deriv_x(agrad_0_beta0, beta0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_beta0, beta0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_beta0, beta0, hx, sz, beta2, bflag);
    deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
    deriv_y(grad_1_beta1, beta1, hx, sz, bflag);
    deriv_z(grad_2_beta1, beta1, hx, sz, bflag);
    deriv_xx(grad2_0_0_beta1, beta1, hx, sz, bflag);
    deriv_y(grad2_0_1_beta1, grad_0_beta1, hx, sz, bflag);
    deriv_z(grad2_0_2_beta1, grad_0_beta1, hx, sz, bflag);
    deriv_yy(grad2_1_1_beta1, beta1, hx, sz, bflag);
    deriv_z(grad2_1_2_beta1, grad_1_beta1, hx, sz, bflag);
    deriv_zz(grad2_2_2_beta1, beta1, hx, sz, bflag);
    adv_deriv_x(agrad_0_beta1, beta1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_beta1, beta1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_beta1, beta1, hx, sz, beta2, bflag);
    deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
    deriv_y(grad_1_beta2, beta2, hx, sz, bflag);
    deriv_z(grad_2_beta2, beta2, hx, sz, bflag);
    deriv_xx(grad2_0_0_beta2, beta2, hx, sz, bflag);
    deriv_y(grad2_0_1_beta2, grad_0_beta2, hx, sz, bflag);
    deriv_z(grad2_0_2_beta2, grad_0_beta2, hx, sz, bflag);
    deriv_yy(grad2_1_1_beta2, beta2, hx, sz, bflag);
    deriv_z(grad2_1_2_beta2, grad_1_beta2, hx, sz, bflag);
    deriv_zz(grad2_2_2_beta2, beta2, hx, sz, bflag);
    adv_deriv_x(agrad_0_beta2, beta2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_beta2, beta2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_beta2, beta2, hx, sz, beta2, bflag);
    deriv_x(grad_0_B0, B0, hx, sz, bflag);
    deriv_y(grad_1_B0, B0, hx, sz, bflag);
    deriv_z(grad_2_B0, B0, hx, sz, bflag);
    adv_deriv_x(agrad_0_B0, B0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_B0, B0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_B0, B0, hx, sz, beta2, bflag);
    deriv_x(grad_0_B1, B1, hx, sz, bflag);
    deriv_y(grad_1_B1, B1, hx, sz, bflag);
    deriv_z(grad_2_B1, B1, hx, sz, bflag);
    adv_deriv_x(agrad_0_B1, B1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_B1, B1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_B1, B1, hx, sz, beta2, bflag);
    deriv_x(grad_0_B2, B2, hx, sz, bflag);
    deriv_y(grad_1_B2, B2, hx, sz, bflag);
    deriv_z(grad_2_B2, B2, hx, sz, bflag);
    adv_deriv_x(agrad_0_B2, B2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_B2, B2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_B2, B2, hx, sz, beta2, bflag);
    deriv_x(grad_0_chi, chi, hx, sz, bflag);
    deriv_y(grad_1_chi, chi, hx, sz, bflag);
    deriv_z(grad_2_chi, chi, hx, sz, bflag);
    deriv_xx(grad2_0_0_chi, chi, hx, sz, bflag);
    deriv_y(grad2_0_1_chi, grad_0_chi, hx, sz, bflag);
    deriv_z(grad2_0_2_chi, grad_0_chi, hx, sz, bflag);
    deriv_yy(grad2_1_1_chi, chi, hx, sz, bflag);
    deriv_z(grad2_1_2_chi, grad_1_chi, hx, sz, bflag);
    deriv_zz(grad2_2_2_chi, chi, hx, sz, bflag);
    adv_deriv_x(agrad_0_chi, chi, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_chi, chi, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_chi, chi, hx, sz, beta2, bflag);
    deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
    deriv_y(grad_1_Gt0, Gt0, hx, sz, bflag);
    deriv_z(grad_2_Gt0, Gt0, hx, sz, bflag);
    adv_deriv_x(agrad_0_Gt0, Gt0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_Gt0, Gt0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_Gt0, Gt0, hx, sz, beta2, bflag);
    deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
    deriv_y(grad_1_Gt1, Gt1, hx, sz, bflag);
    deriv_z(grad_2_Gt1, Gt1, hx, sz, bflag);
    adv_deriv_x(agrad_0_Gt1, Gt1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_Gt1, Gt1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_Gt1, Gt1, hx, sz, beta2, bflag);
    deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
    deriv_y(grad_1_Gt2, Gt2, hx, sz, bflag);
    deriv_z(grad_2_Gt2, Gt2, hx, sz, bflag);
    adv_deriv_x(agrad_0_Gt2, Gt2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_Gt2, Gt2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_Gt2, Gt2, hx, sz, beta2, bflag);
    deriv_x(grad_0_K, K, hx, sz, bflag);
    deriv_y(grad_1_K, K, hx, sz, bflag);
    deriv_z(grad_2_K, K, hx, sz, bflag);
    adv_deriv_x(agrad_0_K, K, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_K, K, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_K, K, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
    deriv_y(grad_1_gt0, gt0, hx, sz, bflag);
    deriv_z(grad_2_gt0, gt0, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt0, gt0, hx, sz, bflag);
    deriv_y(grad2_0_1_gt0, grad_0_gt0, hx, sz, bflag);
    deriv_z(grad2_0_2_gt0, grad_0_gt0, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt0, gt0, hx, sz, bflag);
    deriv_z(grad2_1_2_gt0, grad_1_gt0, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt0, gt0, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt0, gt0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt0, gt0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt0, gt0, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
    deriv_y(grad_1_gt1, gt1, hx, sz, bflag);
    deriv_z(grad_2_gt1, gt1, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt1, gt1, hx, sz, bflag);
    deriv_y(grad2_0_1_gt1, grad_0_gt1, hx, sz, bflag);
    deriv_z(grad2_0_2_gt1, grad_0_gt1, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt1, gt1, hx, sz, bflag);
    deriv_z(grad2_1_2_gt1, grad_1_gt1, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt1, gt1, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt1, gt1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt1, gt1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt1, gt1, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
    deriv_y(grad_1_gt2, gt2, hx, sz, bflag);
    deriv_z(grad_2_gt2, gt2, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt2, gt2, hx, sz, bflag);
    deriv_y(grad2_0_1_gt2, grad_0_gt2, hx, sz, bflag);
    deriv_z(grad2_0_2_gt2, grad_0_gt2, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt2, gt2, hx, sz, bflag);
    deriv_z(grad2_1_2_gt2, grad_1_gt2, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt2, gt2, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt2, gt2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt2, gt2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt2, gt2, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
    deriv_y(grad_1_gt3, gt3, hx, sz, bflag);
    deriv_z(grad_2_gt3, gt3, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt3, gt3, hx, sz, bflag);
    deriv_y(grad2_0_1_gt3, grad_0_gt3, hx, sz, bflag);
    deriv_z(grad2_0_2_gt3, grad_0_gt3, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt3, gt3, hx, sz, bflag);
    deriv_z(grad2_1_2_gt3, grad_1_gt3, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt3, gt3, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt3, gt3, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt3, gt3, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt3, gt3, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
    deriv_y(grad_1_gt4, gt4, hx, sz, bflag);
    deriv_z(grad_2_gt4, gt4, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt4, gt4, hx, sz, bflag);
    deriv_y(grad2_0_1_gt4, grad_0_gt4, hx, sz, bflag);
    deriv_z(grad2_0_2_gt4, grad_0_gt4, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt4, gt4, hx, sz, bflag);
    deriv_z(grad2_1_2_gt4, grad_1_gt4, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt4, gt4, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt4, gt4, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt4, gt4, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt4, gt4, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
    deriv_y(grad_1_gt5, gt5, hx, sz, bflag);
    deriv_z(grad_2_gt5, gt5, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt5, gt5, hx, sz, bflag);
    deriv_y(grad2_0_1_gt5, grad_0_gt5, hx, sz, bflag);
    deriv_z(grad2_0_2_gt5, grad_0_gt5, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt5, gt5, hx, sz, bflag);
    deriv_z(grad2_1_2_gt5, grad_1_gt5, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt5, gt5, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt5, gt5, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt5, gt5, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt5, gt5, hx, sz, beta2, bflag);
    deriv_x(grad_0_At0, At0, hx, sz, bflag);
    deriv_y(grad_1_At0, At0, hx, sz, bflag);
    deriv_z(grad_2_At0, At0, hx, sz, bflag);
    adv_deriv_x(agrad_0_At0, At0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At0, At0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At0, At0, hx, sz, beta2, bflag);
    deriv_x(grad_0_At1, At1, hx, sz, bflag);
    deriv_y(grad_1_At1, At1, hx, sz, bflag);
    deriv_z(grad_2_At1, At1, hx, sz, bflag);
    adv_deriv_x(agrad_0_At1, At1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At1, At1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At1, At1, hx, sz, beta2, bflag);
    deriv_x(grad_0_At2, At2, hx, sz, bflag);
    deriv_y(grad_1_At2, At2, hx, sz, bflag);
    deriv_z(grad_2_At2, At2, hx, sz, bflag);
    adv_deriv_x(agrad_0_At2, At2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At2, At2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At2, At2, hx, sz, beta2, bflag);
    deriv_x(grad_0_At3, At3, hx, sz, bflag);
    deriv_y(grad_1_At3, At3, hx, sz, bflag);
    deriv_z(grad_2_At3, At3, hx, sz, bflag);
    adv_deriv_x(agrad_0_At3, At3, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At3, At3, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At3, At3, hx, sz, beta2, bflag);
    deriv_x(grad_0_At4, At4, hx, sz, bflag);
    deriv_y(grad_1_At4, At4, hx, sz, bflag);
    deriv_z(grad_2_At4, At4, hx, sz, bflag);
    adv_deriv_x(agrad_0_At4, At4, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At4, At4, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At4, At4, hx, sz, beta2, bflag);
    deriv_x(grad_0_At5, At5, hx, sz, bflag);
    deriv_y(grad_1_At5, At5, hx, sz, bflag);
    deriv_z(grad_2_At5, At5, hx, sz, bflag);
    adv_deriv_x(agrad_0_At5, At5, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At5, At5, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At5, At5, hx, sz, beta2, bflag);
    //[[[end]]]

    bssn::timer::t_deriv.stop();

    
    double eta;

    //cout << "begin loop" << endl;
    for (unsigned int k = 3; k < nz-3; k++) {
        const double z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            const double y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                const double x = pmin[0] + i*hx;
                const unsigned int pp = i + nx*(j + ny*k);
                const double r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }


                bssn::timer::t_rhs.start();

                /*[[[cog

                import dendro
                import bssn

                # modify to use the eta const value.
                bssn.eta_func=bssn.eta
                bssn.B_rhs = [bssn.Gt_rhs[i] - bssn.eta_func * bssn.B[i] +
                         bssn.l3 * dendro.vec_j_ad_j(bssn.b, bssn.B[i]) -
                         bssn.l4 * dendro.vec_j_ad_j(bssn.b, bssn.Gt[i])
                         for i in dendro.e_i]

                outs = [bssn.a_rhs, bssn.b_rhs, bssn.gt_rhs, bssn.chi_rhs, bssn.At_rhs, bssn.K_rhs, bssn.Gt_rhs, bssn.B_rhs]
                vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']
                #outs=outs[0:3]
                #vnames=vnames[0:3]
                #dendro.generate_cpu(outs, vnames, '[pp]')


                ]]]*/
                //[[[end]]]

                //#include "bssn.cpp"

                bssn::timer::t_rhs.stop();

            }
        }
    }


    if (bflag != 0) {

        bssn::timer::t_bdyc.start();

        bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, pmin, pmax,
                 1.0, 1.0, sz, bflag);


        bssn::timer::t_bdyc.stop();
    }


    bssn::timer::t_deriv.start();
    // ko derives.
    /*[[[cog
    import cog
    import bssnDerivs as bssnDerivs

    for var in bssnDerivs.KO:
        cog.outl("\t ko_deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[0] + var ,var))
        cog.outl("\t ko_deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[1] + var ,var))
        cog.outl("\t ko_deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[2] + var ,var))

    ]]]*/
    ko_deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
    ko_deriv_y(grad_1_gt0, gt0, hx, sz, bflag);
    ko_deriv_z(grad_2_gt0, gt0, hx, sz, bflag);
    ko_deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
    ko_deriv_y(grad_1_gt1, gt1, hx, sz, bflag);
    ko_deriv_z(grad_2_gt1, gt1, hx, sz, bflag);
    ko_deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
    ko_deriv_y(grad_1_gt2, gt2, hx, sz, bflag);
    ko_deriv_z(grad_2_gt2, gt2, hx, sz, bflag);
    ko_deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
    ko_deriv_y(grad_1_gt3, gt3, hx, sz, bflag);
    ko_deriv_z(grad_2_gt3, gt3, hx, sz, bflag);
    ko_deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
    ko_deriv_y(grad_1_gt4, gt4, hx, sz, bflag);
    ko_deriv_z(grad_2_gt4, gt4, hx, sz, bflag);
    ko_deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
    ko_deriv_y(grad_1_gt5, gt5, hx, sz, bflag);
    ko_deriv_z(grad_2_gt5, gt5, hx, sz, bflag);
    ko_deriv_x(grad_0_At0, At0, hx, sz, bflag);
    ko_deriv_y(grad_1_At0, At0, hx, sz, bflag);
    ko_deriv_z(grad_2_At0, At0, hx, sz, bflag);
    ko_deriv_x(grad_0_At1, At1, hx, sz, bflag);
    ko_deriv_y(grad_1_At1, At1, hx, sz, bflag);
    ko_deriv_z(grad_2_At1, At1, hx, sz, bflag);
    ko_deriv_x(grad_0_At2, At2, hx, sz, bflag);
    ko_deriv_y(grad_1_At2, At2, hx, sz, bflag);
    ko_deriv_z(grad_2_At2, At2, hx, sz, bflag);
    ko_deriv_x(grad_0_At3, At3, hx, sz, bflag);
    ko_deriv_y(grad_1_At3, At3, hx, sz, bflag);
    ko_deriv_z(grad_2_At3, At3, hx, sz, bflag);
    ko_deriv_x(grad_0_At4, At4, hx, sz, bflag);
    ko_deriv_y(grad_1_At4, At4, hx, sz, bflag);
    ko_deriv_z(grad_2_At4, At4, hx, sz, bflag);
    ko_deriv_x(grad_0_At5, At5, hx, sz, bflag);
    ko_deriv_y(grad_1_At5, At5, hx, sz, bflag);
    ko_deriv_z(grad_2_At5, At5, hx, sz, bflag);
    ko_deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
    ko_deriv_y(grad_1_alpha, alpha, hx, sz, bflag);
    ko_deriv_z(grad_2_alpha, alpha, hx, sz, bflag);
    ko_deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
    ko_deriv_y(grad_1_beta0, beta0, hx, sz, bflag);
    ko_deriv_z(grad_2_beta0, beta0, hx, sz, bflag);
    ko_deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
    ko_deriv_y(grad_1_beta1, beta1, hx, sz, bflag);
    ko_deriv_z(grad_2_beta1, beta1, hx, sz, bflag);
    ko_deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
    ko_deriv_y(grad_1_beta2, beta2, hx, sz, bflag);
    ko_deriv_z(grad_2_beta2, beta2, hx, sz, bflag);
    ko_deriv_x(grad_0_chi, chi, hx, sz, bflag);
    ko_deriv_y(grad_1_chi, chi, hx, sz, bflag);
    ko_deriv_z(grad_2_chi, chi, hx, sz, bflag);
    ko_deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
    ko_deriv_y(grad_1_Gt0, Gt0, hx, sz, bflag);
    ko_deriv_z(grad_2_Gt0, Gt0, hx, sz, bflag);
    ko_deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
    ko_deriv_y(grad_1_Gt1, Gt1, hx, sz, bflag);
    ko_deriv_z(grad_2_Gt1, Gt1, hx, sz, bflag);
    ko_deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
    ko_deriv_y(grad_1_Gt2, Gt2, hx, sz, bflag);
    ko_deriv_z(grad_2_Gt2, Gt2, hx, sz, bflag);
    ko_deriv_x(grad_0_K, K, hx, sz, bflag);
    ko_deriv_y(grad_1_K, K, hx, sz, bflag);
    ko_deriv_z(grad_2_K, K, hx, sz, bflag);
    ko_deriv_x(grad_0_B0, B0, hx, sz, bflag);
    ko_deriv_y(grad_1_B0, B0, hx, sz, bflag);
    ko_deriv_z(grad_2_B0, B0, hx, sz, bflag);
    ko_deriv_x(grad_0_B1, B1, hx, sz, bflag);
    ko_deriv_y(grad_1_B1, B1, hx, sz, bflag);
    ko_deriv_z(grad_2_B1, B1, hx, sz, bflag);
    ko_deriv_x(grad_0_B2, B2, hx, sz, bflag);
    ko_deriv_y(grad_1_B2, B2, hx, sz, bflag);
    ko_deriv_z(grad_2_B2, B2, hx, sz, bflag);
    //[[[end]]]

    bssn::timer::t_deriv.stop();

    bssn::timer::t_rhs.start();

    const  double sigma = KO_DISS_SIGMA;


    for (unsigned int k = 3; k < nz-3; k++) {
        for (unsigned int j = 3; j < ny-3; j++) {
            for (unsigned int i = 3; i < nx-3; i++) {
                const unsigned int pp = i + nx*(j + ny*k); 

                a_rhs[pp]  += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
                b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
                b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
                b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

                gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

                chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

                At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

                K_rhs[pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                Gt_rhs0[pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                Gt_rhs1[pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                Gt_rhs2[pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

                B_rhs0[pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                B_rhs1[pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                B_rhs2[pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
            }
        }
    }

    bssn::timer::t_rhs.stop();



    bssn::timer::t_deriv.start();
    /*[[[cog
    import cog
    import bssnDerivs as bssnDerivs

    for deriv in bssnDerivs.FUNC_D_I:
        cog.outl("\t free(%s);" %(deriv))

    for deriv in bssnDerivs.FUNC_D_IJ:
        cog.outl("\t free(%s);" %(deriv))

    for deriv in bssnDerivs.FUNC_AD_I:
        cog.outl("\t free(%s);" %(deriv))

    ]]]*/
    free(grad_0_alpha);
    free(grad_1_alpha);
    free(grad_2_alpha);
    free(grad_0_beta0);
    free(grad_1_beta0);
    free(grad_2_beta0);
    free(grad_0_beta1);
    free(grad_1_beta1);
    free(grad_2_beta1);
    free(grad_0_beta2);
    free(grad_1_beta2);
    free(grad_2_beta2);
    free(grad_0_B0);
    free(grad_1_B0);
    free(grad_2_B0);
    free(grad_0_B1);
    free(grad_1_B1);
    free(grad_2_B1);
    free(grad_0_B2);
    free(grad_1_B2);
    free(grad_2_B2);
    free(grad_0_chi);
    free(grad_1_chi);
    free(grad_2_chi);
    free(grad_0_Gt0);
    free(grad_1_Gt0);
    free(grad_2_Gt0);
    free(grad_0_Gt1);
    free(grad_1_Gt1);
    free(grad_2_Gt1);
    free(grad_0_Gt2);
    free(grad_1_Gt2);
    free(grad_2_Gt2);
    free(grad_0_K);
    free(grad_1_K);
    free(grad_2_K);
    free(grad_0_gt0);
    free(grad_1_gt0);
    free(grad_2_gt0);
    free(grad_0_gt1);
    free(grad_1_gt1);
    free(grad_2_gt1);
    free(grad_0_gt2);
    free(grad_1_gt2);
    free(grad_2_gt2);
    free(grad_0_gt3);
    free(grad_1_gt3);
    free(grad_2_gt3);
    free(grad_0_gt4);
    free(grad_1_gt4);
    free(grad_2_gt4);
    free(grad_0_gt5);
    free(grad_1_gt5);
    free(grad_2_gt5);
    free(grad_0_At0);
    free(grad_1_At0);
    free(grad_2_At0);
    free(grad_0_At1);
    free(grad_1_At1);
    free(grad_2_At1);
    free(grad_0_At2);
    free(grad_1_At2);
    free(grad_2_At2);
    free(grad_0_At3);
    free(grad_1_At3);
    free(grad_2_At3);
    free(grad_0_At4);
    free(grad_1_At4);
    free(grad_2_At4);
    free(grad_0_At5);
    free(grad_1_At5);
    free(grad_2_At5);
    free(grad2_0_0_gt0);
    free(grad2_0_1_gt0);
    free(grad2_0_2_gt0);
    free(grad2_1_1_gt0);
    free(grad2_1_2_gt0);
    free(grad2_2_2_gt0);
    free(grad2_0_0_gt1);
    free(grad2_0_1_gt1);
    free(grad2_0_2_gt1);
    free(grad2_1_1_gt1);
    free(grad2_1_2_gt1);
    free(grad2_2_2_gt1);
    free(grad2_0_0_gt2);
    free(grad2_0_1_gt2);
    free(grad2_0_2_gt2);
    free(grad2_1_1_gt2);
    free(grad2_1_2_gt2);
    free(grad2_2_2_gt2);
    free(grad2_0_0_gt3);
    free(grad2_0_1_gt3);
    free(grad2_0_2_gt3);
    free(grad2_1_1_gt3);
    free(grad2_1_2_gt3);
    free(grad2_2_2_gt3);
    free(grad2_0_0_gt4);
    free(grad2_0_1_gt4);
    free(grad2_0_2_gt4);
    free(grad2_1_1_gt4);
    free(grad2_1_2_gt4);
    free(grad2_2_2_gt4);
    free(grad2_0_0_gt5);
    free(grad2_0_1_gt5);
    free(grad2_0_2_gt5);
    free(grad2_1_1_gt5);
    free(grad2_1_2_gt5);
    free(grad2_2_2_gt5);
    free(grad2_0_0_chi);
    free(grad2_0_1_chi);
    free(grad2_0_2_chi);
    free(grad2_1_1_chi);
    free(grad2_1_2_chi);
    free(grad2_2_2_chi);
    free(grad2_0_0_alpha);
    free(grad2_0_1_alpha);
    free(grad2_0_2_alpha);
    free(grad2_1_1_alpha);
    free(grad2_1_2_alpha);
    free(grad2_2_2_alpha);
    free(grad2_0_0_beta0);
    free(grad2_0_1_beta0);
    free(grad2_0_2_beta0);
    free(grad2_1_1_beta0);
    free(grad2_1_2_beta0);
    free(grad2_2_2_beta0);
    free(grad2_0_0_beta1);
    free(grad2_0_1_beta1);
    free(grad2_0_2_beta1);
    free(grad2_1_1_beta1);
    free(grad2_1_2_beta1);
    free(grad2_2_2_beta1);
    free(grad2_0_0_beta2);
    free(grad2_0_1_beta2);
    free(grad2_0_2_beta2);
    free(grad2_1_1_beta2);
    free(grad2_1_2_beta2);
    free(grad2_2_2_beta2);
    free(agrad_0_gt0);
    free(agrad_1_gt0);
    free(agrad_2_gt0);
    free(agrad_0_gt1);
    free(agrad_1_gt1);
    free(agrad_2_gt1);
    free(agrad_0_gt2);
    free(agrad_1_gt2);
    free(agrad_2_gt2);
    free(agrad_0_gt3);
    free(agrad_1_gt3);
    free(agrad_2_gt3);
    free(agrad_0_gt4);
    free(agrad_1_gt4);
    free(agrad_2_gt4);
    free(agrad_0_gt5);
    free(agrad_1_gt5);
    free(agrad_2_gt5);
    free(agrad_0_At0);
    free(agrad_1_At0);
    free(agrad_2_At0);
    free(agrad_0_At1);
    free(agrad_1_At1);
    free(agrad_2_At1);
    free(agrad_0_At2);
    free(agrad_1_At2);
    free(agrad_2_At2);
    free(agrad_0_At3);
    free(agrad_1_At3);
    free(agrad_2_At3);
    free(agrad_0_At4);
    free(agrad_1_At4);
    free(agrad_2_At4);
    free(agrad_0_At5);
    free(agrad_1_At5);
    free(agrad_2_At5);
    free(agrad_0_alpha);
    free(agrad_1_alpha);
    free(agrad_2_alpha);
    free(agrad_0_beta0);
    free(agrad_1_beta0);
    free(agrad_2_beta0);
    free(agrad_0_beta1);
    free(agrad_1_beta1);
    free(agrad_2_beta1);
    free(agrad_0_beta2);
    free(agrad_1_beta2);
    free(agrad_2_beta2);
    free(agrad_0_chi);
    free(agrad_1_chi);
    free(agrad_2_chi);
    free(agrad_0_Gt0);
    free(agrad_1_Gt0);
    free(agrad_2_Gt0);
    free(agrad_0_Gt1);
    free(agrad_1_Gt1);
    free(agrad_2_Gt1);
    free(agrad_0_Gt2);
    free(agrad_1_Gt2);
    free(agrad_2_Gt2);
    free(agrad_0_K);
    free(agrad_1_K);
    free(agrad_2_K);
    free(agrad_0_B0);
    free(agrad_1_B0);
    free(agrad_2_B0);
    free(agrad_0_B1);
    free(agrad_1_B1);
    free(agrad_2_B1);
    free(agrad_0_B2);
    free(agrad_1_B2);
    free(agrad_2_B2);
    //[[[end]]]
    bssn::timer::t_deriv.stop();


}



void bssnrhs_auto_sep(double **unzipVarsRHS, const double **uZipVars,
             const unsigned int& offset,
             const double *pmin, const double *pmax, const unsigned int *sz,
             const unsigned int& bflag)
{



    const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
    const double *chi = &uZipVars[VAR::U_CHI][offset];
    const double *K = &uZipVars[VAR::U_K][offset];
    const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
    const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
    const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
    const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
    const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
    const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
    const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
    const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
    const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
    const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
    const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
    const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
    const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
    const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
    const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
    const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
    const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
    const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
    const double *B0 = &uZipVars[VAR::U_B0][offset];
    const double *B1 = &uZipVars[VAR::U_B1][offset];
    const double *B2 = &uZipVars[VAR::U_B2][offset];

    double *a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
    double *chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
    double *K_rhs = &unzipVarsRHS[VAR::U_K][offset];
    double *gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
    double *gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
    double *gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
    double *gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
    double *gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
    double *gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
    double *b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
    double *b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
    double *b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
    double *At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
    double *At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
    double *At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
    double *At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
    double *At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
    double *At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
    double *Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
    double *Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
    double *Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
    double *B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
    double *B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
    double *B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx = (pmax[0] - pmin[0]) / (nx - 1);
    double hy = (pmax[1] - pmin[1]) / (ny - 1);
    double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int lambda[4] = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                    BSSN_LAMBDA[2], BSSN_LAMBDA[3]
                                   };
    const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};



    int idx[3];

    unsigned int n = sz[0]*sz[1]*sz[2];

    bssn::timer::t_deriv.start();


    /*[[[cog
    import cog
    import bssnDerivs as bssnDerivs

    for deriv in bssnDerivs.FUNC_D_I:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")

    for deriv in bssnDerivs.FUNC_D_IJ:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")

    for deriv in bssnDerivs.FUNC_AD_I:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")


    for var in bssnDerivs.D:
        cog.outl("\t deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[0] + var ,var))
        cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[1] + var ,var))
        cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[2] + var ,var))

        if var in bssnDerivs.DD:
            cog.outl("\t deriv_xx(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[0] + var ,var))
            cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[1] + var , bssnDerivs.PREFIX_D[0] + var ))
            cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[2] + var , bssnDerivs.PREFIX_D[0] + var ))

            cog.outl("\t deriv_yy(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[3] + var ,var))
            cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[4] + var , bssnDerivs.PREFIX_D[1] + var))

            cog.outl("\t deriv_zz(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[5] + var ,var))

        if var in bssnDerivs.AD:
            cog.outl("\t adv_deriv_x(%s, %s, hx, sz, beta0, bflag);" %(bssnDerivs.PREFIX_AD[0] + var ,var))
            cog.outl("\t adv_deriv_y(%s, %s, hx, sz, beta1, bflag);" %(bssnDerivs.PREFIX_AD[1] + var ,var))
            cog.outl("\t adv_deriv_z(%s, %s, hx, sz, beta2, bflag);" %(bssnDerivs.PREFIX_AD[2] + var ,var))



    ]]]*/
    double* grad_0_alpha = (double*)malloc(sizeof(double)*n);
    double* grad_1_alpha = (double*)malloc(sizeof(double)*n);
    double* grad_2_alpha = (double*)malloc(sizeof(double)*n);
    double* grad_0_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_B0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_B0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_B0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_B1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_B1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_B1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_B2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_B2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_B2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_chi = (double*)malloc(sizeof(double)*n);
    double* grad_1_chi = (double*)malloc(sizeof(double)*n);
    double* grad_2_chi = (double*)malloc(sizeof(double)*n);
    double* grad_0_Gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_Gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_Gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_Gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_Gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_Gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_Gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_Gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_Gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_K = (double*)malloc(sizeof(double)*n);
    double* grad_1_K = (double*)malloc(sizeof(double)*n);
    double* grad_2_K = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad_0_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad_1_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At0 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At0 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At0 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At1 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At1 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At1 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At2 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At2 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At2 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At3 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At3 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At3 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At4 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At4 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At4 = (double*)malloc(sizeof(double)*n);
    double* grad_0_At5 = (double*)malloc(sizeof(double)*n);
    double* grad_1_At5 = (double*)malloc(sizeof(double)*n);
    double* grad_2_At5 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_chi = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_alpha = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_0_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_1_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_0_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_1_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_1_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* grad2_2_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt3 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt3 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt3 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt4 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt4 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt4 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_gt5 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_gt5 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_gt5 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At3 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At3 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At3 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At4 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At4 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At4 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_At5 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_At5 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_At5 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_alpha = (double*)malloc(sizeof(double)*n);
    double* agrad_1_alpha = (double*)malloc(sizeof(double)*n);
    double* agrad_2_alpha = (double*)malloc(sizeof(double)*n);
    double* agrad_0_beta0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_beta0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_beta0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_beta1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_beta1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_beta1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_beta2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_beta2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_beta2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_chi = (double*)malloc(sizeof(double)*n);
    double* agrad_1_chi = (double*)malloc(sizeof(double)*n);
    double* agrad_2_chi = (double*)malloc(sizeof(double)*n);
    double* agrad_0_Gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_Gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_Gt0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_Gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_Gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_Gt1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_Gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_Gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_Gt2 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_K = (double*)malloc(sizeof(double)*n);
    double* agrad_1_K = (double*)malloc(sizeof(double)*n);
    double* agrad_2_K = (double*)malloc(sizeof(double)*n);
    double* agrad_0_B0 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_B0 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_B0 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_B1 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_B1 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_B1 = (double*)malloc(sizeof(double)*n);
    double* agrad_0_B2 = (double*)malloc(sizeof(double)*n);
    double* agrad_1_B2 = (double*)malloc(sizeof(double)*n);
    double* agrad_2_B2 = (double*)malloc(sizeof(double)*n);
    deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
    deriv_y(grad_1_alpha, alpha, hx, sz, bflag);
    deriv_z(grad_2_alpha, alpha, hx, sz, bflag);
    deriv_xx(grad2_0_0_alpha, alpha, hx, sz, bflag);
    deriv_y(grad2_0_1_alpha, grad_0_alpha, hx, sz, bflag);
    deriv_z(grad2_0_2_alpha, grad_0_alpha, hx, sz, bflag);
    deriv_yy(grad2_1_1_alpha, alpha, hx, sz, bflag);
    deriv_z(grad2_1_2_alpha, grad_1_alpha, hx, sz, bflag);
    deriv_zz(grad2_2_2_alpha, alpha, hx, sz, bflag);
    adv_deriv_x(agrad_0_alpha, alpha, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_alpha, alpha, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_alpha, alpha, hx, sz, beta2, bflag);
    deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
    deriv_y(grad_1_beta0, beta0, hx, sz, bflag);
    deriv_z(grad_2_beta0, beta0, hx, sz, bflag);
    deriv_xx(grad2_0_0_beta0, beta0, hx, sz, bflag);
    deriv_y(grad2_0_1_beta0, grad_0_beta0, hx, sz, bflag);
    deriv_z(grad2_0_2_beta0, grad_0_beta0, hx, sz, bflag);
    deriv_yy(grad2_1_1_beta0, beta0, hx, sz, bflag);
    deriv_z(grad2_1_2_beta0, grad_1_beta0, hx, sz, bflag);
    deriv_zz(grad2_2_2_beta0, beta0, hx, sz, bflag);
    adv_deriv_x(agrad_0_beta0, beta0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_beta0, beta0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_beta0, beta0, hx, sz, beta2, bflag);
    deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
    deriv_y(grad_1_beta1, beta1, hx, sz, bflag);
    deriv_z(grad_2_beta1, beta1, hx, sz, bflag);
    deriv_xx(grad2_0_0_beta1, beta1, hx, sz, bflag);
    deriv_y(grad2_0_1_beta1, grad_0_beta1, hx, sz, bflag);
    deriv_z(grad2_0_2_beta1, grad_0_beta1, hx, sz, bflag);
    deriv_yy(grad2_1_1_beta1, beta1, hx, sz, bflag);
    deriv_z(grad2_1_2_beta1, grad_1_beta1, hx, sz, bflag);
    deriv_zz(grad2_2_2_beta1, beta1, hx, sz, bflag);
    adv_deriv_x(agrad_0_beta1, beta1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_beta1, beta1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_beta1, beta1, hx, sz, beta2, bflag);
    deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
    deriv_y(grad_1_beta2, beta2, hx, sz, bflag);
    deriv_z(grad_2_beta2, beta2, hx, sz, bflag);
    deriv_xx(grad2_0_0_beta2, beta2, hx, sz, bflag);
    deriv_y(grad2_0_1_beta2, grad_0_beta2, hx, sz, bflag);
    deriv_z(grad2_0_2_beta2, grad_0_beta2, hx, sz, bflag);
    deriv_yy(grad2_1_1_beta2, beta2, hx, sz, bflag);
    deriv_z(grad2_1_2_beta2, grad_1_beta2, hx, sz, bflag);
    deriv_zz(grad2_2_2_beta2, beta2, hx, sz, bflag);
    adv_deriv_x(agrad_0_beta2, beta2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_beta2, beta2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_beta2, beta2, hx, sz, beta2, bflag);
    deriv_x(grad_0_B0, B0, hx, sz, bflag);
    deriv_y(grad_1_B0, B0, hx, sz, bflag);
    deriv_z(grad_2_B0, B0, hx, sz, bflag);
    adv_deriv_x(agrad_0_B0, B0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_B0, B0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_B0, B0, hx, sz, beta2, bflag);
    deriv_x(grad_0_B1, B1, hx, sz, bflag);
    deriv_y(grad_1_B1, B1, hx, sz, bflag);
    deriv_z(grad_2_B1, B1, hx, sz, bflag);
    adv_deriv_x(agrad_0_B1, B1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_B1, B1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_B1, B1, hx, sz, beta2, bflag);
    deriv_x(grad_0_B2, B2, hx, sz, bflag);
    deriv_y(grad_1_B2, B2, hx, sz, bflag);
    deriv_z(grad_2_B2, B2, hx, sz, bflag);
    adv_deriv_x(agrad_0_B2, B2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_B2, B2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_B2, B2, hx, sz, beta2, bflag);
    deriv_x(grad_0_chi, chi, hx, sz, bflag);
    deriv_y(grad_1_chi, chi, hx, sz, bflag);
    deriv_z(grad_2_chi, chi, hx, sz, bflag);
    deriv_xx(grad2_0_0_chi, chi, hx, sz, bflag);
    deriv_y(grad2_0_1_chi, grad_0_chi, hx, sz, bflag);
    deriv_z(grad2_0_2_chi, grad_0_chi, hx, sz, bflag);
    deriv_yy(grad2_1_1_chi, chi, hx, sz, bflag);
    deriv_z(grad2_1_2_chi, grad_1_chi, hx, sz, bflag);
    deriv_zz(grad2_2_2_chi, chi, hx, sz, bflag);
    adv_deriv_x(agrad_0_chi, chi, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_chi, chi, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_chi, chi, hx, sz, beta2, bflag);
    deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
    deriv_y(grad_1_Gt0, Gt0, hx, sz, bflag);
    deriv_z(grad_2_Gt0, Gt0, hx, sz, bflag);
    adv_deriv_x(agrad_0_Gt0, Gt0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_Gt0, Gt0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_Gt0, Gt0, hx, sz, beta2, bflag);
    deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
    deriv_y(grad_1_Gt1, Gt1, hx, sz, bflag);
    deriv_z(grad_2_Gt1, Gt1, hx, sz, bflag);
    adv_deriv_x(agrad_0_Gt1, Gt1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_Gt1, Gt1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_Gt1, Gt1, hx, sz, beta2, bflag);
    deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
    deriv_y(grad_1_Gt2, Gt2, hx, sz, bflag);
    deriv_z(grad_2_Gt2, Gt2, hx, sz, bflag);
    adv_deriv_x(agrad_0_Gt2, Gt2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_Gt2, Gt2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_Gt2, Gt2, hx, sz, beta2, bflag);
    deriv_x(grad_0_K, K, hx, sz, bflag);
    deriv_y(grad_1_K, K, hx, sz, bflag);
    deriv_z(grad_2_K, K, hx, sz, bflag);
    adv_deriv_x(agrad_0_K, K, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_K, K, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_K, K, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
    deriv_y(grad_1_gt0, gt0, hx, sz, bflag);
    deriv_z(grad_2_gt0, gt0, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt0, gt0, hx, sz, bflag);
    deriv_y(grad2_0_1_gt0, grad_0_gt0, hx, sz, bflag);
    deriv_z(grad2_0_2_gt0, grad_0_gt0, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt0, gt0, hx, sz, bflag);
    deriv_z(grad2_1_2_gt0, grad_1_gt0, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt0, gt0, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt0, gt0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt0, gt0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt0, gt0, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
    deriv_y(grad_1_gt1, gt1, hx, sz, bflag);
    deriv_z(grad_2_gt1, gt1, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt1, gt1, hx, sz, bflag);
    deriv_y(grad2_0_1_gt1, grad_0_gt1, hx, sz, bflag);
    deriv_z(grad2_0_2_gt1, grad_0_gt1, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt1, gt1, hx, sz, bflag);
    deriv_z(grad2_1_2_gt1, grad_1_gt1, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt1, gt1, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt1, gt1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt1, gt1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt1, gt1, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
    deriv_y(grad_1_gt2, gt2, hx, sz, bflag);
    deriv_z(grad_2_gt2, gt2, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt2, gt2, hx, sz, bflag);
    deriv_y(grad2_0_1_gt2, grad_0_gt2, hx, sz, bflag);
    deriv_z(grad2_0_2_gt2, grad_0_gt2, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt2, gt2, hx, sz, bflag);
    deriv_z(grad2_1_2_gt2, grad_1_gt2, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt2, gt2, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt2, gt2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt2, gt2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt2, gt2, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
    deriv_y(grad_1_gt3, gt3, hx, sz, bflag);
    deriv_z(grad_2_gt3, gt3, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt3, gt3, hx, sz, bflag);
    deriv_y(grad2_0_1_gt3, grad_0_gt3, hx, sz, bflag);
    deriv_z(grad2_0_2_gt3, grad_0_gt3, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt3, gt3, hx, sz, bflag);
    deriv_z(grad2_1_2_gt3, grad_1_gt3, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt3, gt3, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt3, gt3, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt3, gt3, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt3, gt3, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
    deriv_y(grad_1_gt4, gt4, hx, sz, bflag);
    deriv_z(grad_2_gt4, gt4, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt4, gt4, hx, sz, bflag);
    deriv_y(grad2_0_1_gt4, grad_0_gt4, hx, sz, bflag);
    deriv_z(grad2_0_2_gt4, grad_0_gt4, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt4, gt4, hx, sz, bflag);
    deriv_z(grad2_1_2_gt4, grad_1_gt4, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt4, gt4, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt4, gt4, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt4, gt4, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt4, gt4, hx, sz, beta2, bflag);
    deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
    deriv_y(grad_1_gt5, gt5, hx, sz, bflag);
    deriv_z(grad_2_gt5, gt5, hx, sz, bflag);
    deriv_xx(grad2_0_0_gt5, gt5, hx, sz, bflag);
    deriv_y(grad2_0_1_gt5, grad_0_gt5, hx, sz, bflag);
    deriv_z(grad2_0_2_gt5, grad_0_gt5, hx, sz, bflag);
    deriv_yy(grad2_1_1_gt5, gt5, hx, sz, bflag);
    deriv_z(grad2_1_2_gt5, grad_1_gt5, hx, sz, bflag);
    deriv_zz(grad2_2_2_gt5, gt5, hx, sz, bflag);
    adv_deriv_x(agrad_0_gt5, gt5, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_gt5, gt5, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_gt5, gt5, hx, sz, beta2, bflag);
    deriv_x(grad_0_At0, At0, hx, sz, bflag);
    deriv_y(grad_1_At0, At0, hx, sz, bflag);
    deriv_z(grad_2_At0, At0, hx, sz, bflag);
    adv_deriv_x(agrad_0_At0, At0, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At0, At0, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At0, At0, hx, sz, beta2, bflag);
    deriv_x(grad_0_At1, At1, hx, sz, bflag);
    deriv_y(grad_1_At1, At1, hx, sz, bflag);
    deriv_z(grad_2_At1, At1, hx, sz, bflag);
    adv_deriv_x(agrad_0_At1, At1, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At1, At1, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At1, At1, hx, sz, beta2, bflag);
    deriv_x(grad_0_At2, At2, hx, sz, bflag);
    deriv_y(grad_1_At2, At2, hx, sz, bflag);
    deriv_z(grad_2_At2, At2, hx, sz, bflag);
    adv_deriv_x(agrad_0_At2, At2, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At2, At2, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At2, At2, hx, sz, beta2, bflag);
    deriv_x(grad_0_At3, At3, hx, sz, bflag);
    deriv_y(grad_1_At3, At3, hx, sz, bflag);
    deriv_z(grad_2_At3, At3, hx, sz, bflag);
    adv_deriv_x(agrad_0_At3, At3, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At3, At3, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At3, At3, hx, sz, beta2, bflag);
    deriv_x(grad_0_At4, At4, hx, sz, bflag);
    deriv_y(grad_1_At4, At4, hx, sz, bflag);
    deriv_z(grad_2_At4, At4, hx, sz, bflag);
    adv_deriv_x(agrad_0_At4, At4, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At4, At4, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At4, At4, hx, sz, beta2, bflag);
    deriv_x(grad_0_At5, At5, hx, sz, bflag);
    deriv_y(grad_1_At5, At5, hx, sz, bflag);
    deriv_z(grad_2_At5, At5, hx, sz, bflag);
    adv_deriv_x(agrad_0_At5, At5, hx, sz, beta0, bflag);
    adv_deriv_y(agrad_1_At5, At5, hx, sz, beta1, bflag);
    adv_deriv_z(agrad_2_At5, At5, hx, sz, beta2, bflag);
    //[[[end]]]

    bssn::timer::t_deriv.stop();

    
    double eta;

    //cout << "begin loop" << endl;
    #include "bssn_auto_sep.cpp"


    if (bflag != 0) {

        bssn::timer::t_bdyc.start();

        bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, pmin, pmax,
                 1.0, 1.0, sz, bflag);


        bssn::timer::t_bdyc.stop();
    }


    bssn::timer::t_deriv.start();
    // ko derives.
    /*[[[cog
    import cog
    import bssnDerivs as bssnDerivs

    for var in bssnDerivs.KO:
        cog.outl("\t ko_deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[0] + var ,var))
        cog.outl("\t ko_deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[1] + var ,var))
        cog.outl("\t ko_deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[2] + var ,var))

    ]]]*/
    ko_deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
    ko_deriv_y(grad_1_gt0, gt0, hx, sz, bflag);
    ko_deriv_z(grad_2_gt0, gt0, hx, sz, bflag);
    ko_deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
    ko_deriv_y(grad_1_gt1, gt1, hx, sz, bflag);
    ko_deriv_z(grad_2_gt1, gt1, hx, sz, bflag);
    ko_deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
    ko_deriv_y(grad_1_gt2, gt2, hx, sz, bflag);
    ko_deriv_z(grad_2_gt2, gt2, hx, sz, bflag);
    ko_deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
    ko_deriv_y(grad_1_gt3, gt3, hx, sz, bflag);
    ko_deriv_z(grad_2_gt3, gt3, hx, sz, bflag);
    ko_deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
    ko_deriv_y(grad_1_gt4, gt4, hx, sz, bflag);
    ko_deriv_z(grad_2_gt4, gt4, hx, sz, bflag);
    ko_deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
    ko_deriv_y(grad_1_gt5, gt5, hx, sz, bflag);
    ko_deriv_z(grad_2_gt5, gt5, hx, sz, bflag);
    ko_deriv_x(grad_0_At0, At0, hx, sz, bflag);
    ko_deriv_y(grad_1_At0, At0, hx, sz, bflag);
    ko_deriv_z(grad_2_At0, At0, hx, sz, bflag);
    ko_deriv_x(grad_0_At1, At1, hx, sz, bflag);
    ko_deriv_y(grad_1_At1, At1, hx, sz, bflag);
    ko_deriv_z(grad_2_At1, At1, hx, sz, bflag);
    ko_deriv_x(grad_0_At2, At2, hx, sz, bflag);
    ko_deriv_y(grad_1_At2, At2, hx, sz, bflag);
    ko_deriv_z(grad_2_At2, At2, hx, sz, bflag);
    ko_deriv_x(grad_0_At3, At3, hx, sz, bflag);
    ko_deriv_y(grad_1_At3, At3, hx, sz, bflag);
    ko_deriv_z(grad_2_At3, At3, hx, sz, bflag);
    ko_deriv_x(grad_0_At4, At4, hx, sz, bflag);
    ko_deriv_y(grad_1_At4, At4, hx, sz, bflag);
    ko_deriv_z(grad_2_At4, At4, hx, sz, bflag);
    ko_deriv_x(grad_0_At5, At5, hx, sz, bflag);
    ko_deriv_y(grad_1_At5, At5, hx, sz, bflag);
    ko_deriv_z(grad_2_At5, At5, hx, sz, bflag);
    ko_deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
    ko_deriv_y(grad_1_alpha, alpha, hx, sz, bflag);
    ko_deriv_z(grad_2_alpha, alpha, hx, sz, bflag);
    ko_deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
    ko_deriv_y(grad_1_beta0, beta0, hx, sz, bflag);
    ko_deriv_z(grad_2_beta0, beta0, hx, sz, bflag);
    ko_deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
    ko_deriv_y(grad_1_beta1, beta1, hx, sz, bflag);
    ko_deriv_z(grad_2_beta1, beta1, hx, sz, bflag);
    ko_deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
    ko_deriv_y(grad_1_beta2, beta2, hx, sz, bflag);
    ko_deriv_z(grad_2_beta2, beta2, hx, sz, bflag);
    ko_deriv_x(grad_0_chi, chi, hx, sz, bflag);
    ko_deriv_y(grad_1_chi, chi, hx, sz, bflag);
    ko_deriv_z(grad_2_chi, chi, hx, sz, bflag);
    ko_deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
    ko_deriv_y(grad_1_Gt0, Gt0, hx, sz, bflag);
    ko_deriv_z(grad_2_Gt0, Gt0, hx, sz, bflag);
    ko_deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
    ko_deriv_y(grad_1_Gt1, Gt1, hx, sz, bflag);
    ko_deriv_z(grad_2_Gt1, Gt1, hx, sz, bflag);
    ko_deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
    ko_deriv_y(grad_1_Gt2, Gt2, hx, sz, bflag);
    ko_deriv_z(grad_2_Gt2, Gt2, hx, sz, bflag);
    ko_deriv_x(grad_0_K, K, hx, sz, bflag);
    ko_deriv_y(grad_1_K, K, hx, sz, bflag);
    ko_deriv_z(grad_2_K, K, hx, sz, bflag);
    ko_deriv_x(grad_0_B0, B0, hx, sz, bflag);
    ko_deriv_y(grad_1_B0, B0, hx, sz, bflag);
    ko_deriv_z(grad_2_B0, B0, hx, sz, bflag);
    ko_deriv_x(grad_0_B1, B1, hx, sz, bflag);
    ko_deriv_y(grad_1_B1, B1, hx, sz, bflag);
    ko_deriv_z(grad_2_B1, B1, hx, sz, bflag);
    ko_deriv_x(grad_0_B2, B2, hx, sz, bflag);
    ko_deriv_y(grad_1_B2, B2, hx, sz, bflag);
    ko_deriv_z(grad_2_B2, B2, hx, sz, bflag);
    //[[[end]]]

    bssn::timer::t_deriv.stop();

    bssn::timer::t_rhs.start();

    const  double sigma = KO_DISS_SIGMA;


    for (unsigned int k = 3; k < nz-3; k++) {
        for (unsigned int j = 3; j < ny-3; j++) {
            for (unsigned int i = 3; i < nx-3; i++) {
                const unsigned int pp = i + nx*(j + ny*k);

                a_rhs[pp]  += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
                b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
                b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
                b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

                gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

                chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

                At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

                K_rhs[pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                Gt_rhs0[pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                Gt_rhs1[pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                Gt_rhs2[pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

                B_rhs0[pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                B_rhs1[pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                B_rhs2[pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
            }
        }
    }

    bssn::timer::t_rhs.stop();



    bssn::timer::t_deriv.start();
    /*[[[cog
    import cog
    import bssnDerivs as bssnDerivs

    for deriv in bssnDerivs.FUNC_D_I:
        cog.outl("\t free(%s);" %(deriv))

    for deriv in bssnDerivs.FUNC_D_IJ:
        cog.outl("\t free(%s);" %(deriv))

    for deriv in bssnDerivs.FUNC_AD_I:
        cog.outl("\t free(%s);" %(deriv))

    ]]]*/
    free(grad_0_alpha);
    free(grad_1_alpha);
    free(grad_2_alpha);
    free(grad_0_beta0);
    free(grad_1_beta0);
    free(grad_2_beta0);
    free(grad_0_beta1);
    free(grad_1_beta1);
    free(grad_2_beta1);
    free(grad_0_beta2);
    free(grad_1_beta2);
    free(grad_2_beta2);
    free(grad_0_B0);
    free(grad_1_B0);
    free(grad_2_B0);
    free(grad_0_B1);
    free(grad_1_B1);
    free(grad_2_B1);
    free(grad_0_B2);
    free(grad_1_B2);
    free(grad_2_B2);
    free(grad_0_chi);
    free(grad_1_chi);
    free(grad_2_chi);
    free(grad_0_Gt0);
    free(grad_1_Gt0);
    free(grad_2_Gt0);
    free(grad_0_Gt1);
    free(grad_1_Gt1);
    free(grad_2_Gt1);
    free(grad_0_Gt2);
    free(grad_1_Gt2);
    free(grad_2_Gt2);
    free(grad_0_K);
    free(grad_1_K);
    free(grad_2_K);
    free(grad_0_gt0);
    free(grad_1_gt0);
    free(grad_2_gt0);
    free(grad_0_gt1);
    free(grad_1_gt1);
    free(grad_2_gt1);
    free(grad_0_gt2);
    free(grad_1_gt2);
    free(grad_2_gt2);
    free(grad_0_gt3);
    free(grad_1_gt3);
    free(grad_2_gt3);
    free(grad_0_gt4);
    free(grad_1_gt4);
    free(grad_2_gt4);
    free(grad_0_gt5);
    free(grad_1_gt5);
    free(grad_2_gt5);
    free(grad_0_At0);
    free(grad_1_At0);
    free(grad_2_At0);
    free(grad_0_At1);
    free(grad_1_At1);
    free(grad_2_At1);
    free(grad_0_At2);
    free(grad_1_At2);
    free(grad_2_At2);
    free(grad_0_At3);
    free(grad_1_At3);
    free(grad_2_At3);
    free(grad_0_At4);
    free(grad_1_At4);
    free(grad_2_At4);
    free(grad_0_At5);
    free(grad_1_At5);
    free(grad_2_At5);
    free(grad2_0_0_gt0);
    free(grad2_0_1_gt0);
    free(grad2_0_2_gt0);
    free(grad2_1_1_gt0);
    free(grad2_1_2_gt0);
    free(grad2_2_2_gt0);
    free(grad2_0_0_gt1);
    free(grad2_0_1_gt1);
    free(grad2_0_2_gt1);
    free(grad2_1_1_gt1);
    free(grad2_1_2_gt1);
    free(grad2_2_2_gt1);
    free(grad2_0_0_gt2);
    free(grad2_0_1_gt2);
    free(grad2_0_2_gt2);
    free(grad2_1_1_gt2);
    free(grad2_1_2_gt2);
    free(grad2_2_2_gt2);
    free(grad2_0_0_gt3);
    free(grad2_0_1_gt3);
    free(grad2_0_2_gt3);
    free(grad2_1_1_gt3);
    free(grad2_1_2_gt3);
    free(grad2_2_2_gt3);
    free(grad2_0_0_gt4);
    free(grad2_0_1_gt4);
    free(grad2_0_2_gt4);
    free(grad2_1_1_gt4);
    free(grad2_1_2_gt4);
    free(grad2_2_2_gt4);
    free(grad2_0_0_gt5);
    free(grad2_0_1_gt5);
    free(grad2_0_2_gt5);
    free(grad2_1_1_gt5);
    free(grad2_1_2_gt5);
    free(grad2_2_2_gt5);
    free(grad2_0_0_chi);
    free(grad2_0_1_chi);
    free(grad2_0_2_chi);
    free(grad2_1_1_chi);
    free(grad2_1_2_chi);
    free(grad2_2_2_chi);
    free(grad2_0_0_alpha);
    free(grad2_0_1_alpha);
    free(grad2_0_2_alpha);
    free(grad2_1_1_alpha);
    free(grad2_1_2_alpha);
    free(grad2_2_2_alpha);
    free(grad2_0_0_beta0);
    free(grad2_0_1_beta0);
    free(grad2_0_2_beta0);
    free(grad2_1_1_beta0);
    free(grad2_1_2_beta0);
    free(grad2_2_2_beta0);
    free(grad2_0_0_beta1);
    free(grad2_0_1_beta1);
    free(grad2_0_2_beta1);
    free(grad2_1_1_beta1);
    free(grad2_1_2_beta1);
    free(grad2_2_2_beta1);
    free(grad2_0_0_beta2);
    free(grad2_0_1_beta2);
    free(grad2_0_2_beta2);
    free(grad2_1_1_beta2);
    free(grad2_1_2_beta2);
    free(grad2_2_2_beta2);
    free(agrad_0_gt0);
    free(agrad_1_gt0);
    free(agrad_2_gt0);
    free(agrad_0_gt1);
    free(agrad_1_gt1);
    free(agrad_2_gt1);
    free(agrad_0_gt2);
    free(agrad_1_gt2);
    free(agrad_2_gt2);
    free(agrad_0_gt3);
    free(agrad_1_gt3);
    free(agrad_2_gt3);
    free(agrad_0_gt4);
    free(agrad_1_gt4);
    free(agrad_2_gt4);
    free(agrad_0_gt5);
    free(agrad_1_gt5);
    free(agrad_2_gt5);
    free(agrad_0_At0);
    free(agrad_1_At0);
    free(agrad_2_At0);
    free(agrad_0_At1);
    free(agrad_1_At1);
    free(agrad_2_At1);
    free(agrad_0_At2);
    free(agrad_1_At2);
    free(agrad_2_At2);
    free(agrad_0_At3);
    free(agrad_1_At3);
    free(agrad_2_At3);
    free(agrad_0_At4);
    free(agrad_1_At4);
    free(agrad_2_At4);
    free(agrad_0_At5);
    free(agrad_1_At5);
    free(agrad_2_At5);
    free(agrad_0_alpha);
    free(agrad_1_alpha);
    free(agrad_2_alpha);
    free(agrad_0_beta0);
    free(agrad_1_beta0);
    free(agrad_2_beta0);
    free(agrad_0_beta1);
    free(agrad_1_beta1);
    free(agrad_2_beta1);
    free(agrad_0_beta2);
    free(agrad_1_beta2);
    free(agrad_2_beta2);
    free(agrad_0_chi);
    free(agrad_1_chi);
    free(agrad_2_chi);
    free(agrad_0_Gt0);
    free(agrad_1_Gt0);
    free(agrad_2_Gt0);
    free(agrad_0_Gt1);
    free(agrad_1_Gt1);
    free(agrad_2_Gt1);
    free(agrad_0_Gt2);
    free(agrad_1_Gt2);
    free(agrad_2_Gt2);
    free(agrad_0_K);
    free(agrad_1_K);
    free(agrad_2_K);
    free(agrad_0_B0);
    free(agrad_1_B0);
    free(agrad_2_B0);
    free(agrad_0_B1);
    free(agrad_1_B1);
    free(agrad_2_B1);
    free(agrad_0_B2);
    free(agrad_1_B2);
    free(agrad_2_B2);
    //[[[end]]]
    bssn::timer::t_deriv.stop();


}





/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void bssn_bcs(double *f_rhs, const double *f,
              const double *dxf, const double *dyf, const double *dzf,
              const double *pmin, const double *pmax,
              const double f_falloff, const double f_asymptotic,
              const unsigned int *sz, const unsigned int &bflag)
{

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx = (pmax[0] - pmin[0]) / (nx - 1);
    double hy = (pmax[1] - pmin[1]) / (ny - 1);
    double hz = (pmax[2] - pmin[2]) / (nz - 1);

    unsigned int ib = 3;
    unsigned int jb = 3;
    unsigned int kb = 3;
    unsigned int ie = sz[0]-3;
    unsigned int je = sz[1]-3;
    unsigned int ke = sz[2]-3;

    double x,y,z;
    unsigned int pp;
    double inv_r;

    //std::cout<<"boundary bssnrhs: size [ "<<nx<<", "<<ny<<", "<<nz<<" ]"<<std::endl;
    //std::cout<<"boundary bssnrhs: pmin [ "<<pmin[0]<<", "<<pmin[1]<<", "<<pmin[2]<<" ]"<<std::endl;
    //std::cout<<"boundary bssnrhs: pmax [ "<<pmax[0]<<", "<<pmax[1]<<", "<<pmax[2]<<" ]"<<std::endl;

    if (bflag & (1u<<OCT_DIR_LEFT)) {
        double x = pmin[0] + ib*hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k*hz;
            for (unsigned int j = jb; j < je; j++) {
                y = pmin[1] + j*hy;
                pp = IDX(ib,j,k);
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);

                f_rhs[pp] = -  inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_RIGHT)) {
        x = pmin[0] + (ie-1)*hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k*hz;
            for (unsigned int j = jb; j < je; j++) {
                y = pmin[1] + j*hy;
                pp = IDX((ie-1),j,k);
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);

                f_rhs[pp] = -  inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_DOWN)) {
        y = pmin[1] + jb*hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k*hz;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i*hx;
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);
                pp = IDX(i,jb,k);

                f_rhs[pp] = -  inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_UP)) {
        y = pmin[1] + (je-1)*hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k*hz;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i*hx;
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);
                pp = IDX(i,(je-1),k);

                f_rhs[pp] = -  inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_BACK)) {
        z = pmin[2] + kb*hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j*hy;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i*hx;
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);
                pp = IDX(i,j,kb);

                f_rhs[pp] = - inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_FRONT)) {
        z = pmin[2] + (ke-1)*hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j*hy;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i*hx;
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);
                pp = IDX(i,j,(ke-1));

                f_rhs[pp] = - inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void freeze_bcs(double *f_rhs, const unsigned int *sz, const unsigned int &bflag)
{

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    unsigned int ib = 3;
    unsigned int jb = 3;
    unsigned int kb = 3;
    unsigned int ie = sz[0]-3;
    unsigned int je = sz[1]-3;
    unsigned int ke = sz[2]-3;

    unsigned int pp;

    if (bflag & (1u<<OCT_DIR_LEFT)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int j = jb; j < je; j++) {
                pp = IDX(ib,j,k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_RIGHT)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int j = jb; j < je; j++) {
                pp = IDX((ie-1),j,k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_DOWN)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i,jb,k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_UP)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i,(je-1),k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_BACK)) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i,j,kb);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_FRONT)) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i,j,(ke-1));
                f_rhs[pp] = 0.0;
            }
        }
    }

}

#if 0
/*--------------------------------------------------------------
 * Kerr-Schild data
 *--------------------------------------------------------------*/

void ks_initial_data(double x, double y, double z, double *u)
{

    u[VAR::U_ALPHA] = 0.0;
    u[VAR::U_BETA0] = 0.0;
    u[VAR::U_BETA1] = pi/5.0*cos(y)*sin(z+x);
    u[VAR::U_BETA2] = 4.0/17.0*sin(f2*x)*sin(z);

    u[VAR::U_B0] = 0.0;
    u[VAR::U_B1] = 0.0;
    u[VAR::U_B2] = 0.0;

    u[VAR::U_GT0] = Gamt_1;
    u[VAR::U_GT1] = Gamt_2;
    u[VAR::U_GT2] = Gamt_3;

    u[VAR::U_CHI] = 1.0 + exp(-4.0*cos(x)*sin(y));

    u[VAR::U_SYMGT0] = 1.00+0.2*sin(x+z)*cos(y);
    u[VAR::U_SYMGT3] = 1.00+0.2*cos(y)*cos(z+ x);
    u[VAR::U_SYMGT5] = 1.00 / ( u[VAR::U_SYMGT0] + u[VAR::U_SYMGT3]);
    u[VAR::U_SYMGT1] = 0.7*cos(x*x + y*y);
    u[VAR::U_SYMGT2] = 0.3*sin(z)*cos(x);
    u[VAR::U_SYMGT4] = -0.5*sin(x*x)*cos(y)*cos(z);

    u[VAR::U_K] = 5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)
                  +5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)
                  +0.4*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))
                  *exp(-4.0*cos(x)*sin(y))*cos(z);

    u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x) -0.3333333333*exp(4.0*cos(x)*sin(y)) *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y) +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
    u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
    u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);
    u[VAR::U_SYMAT3] = exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
    u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);
    u[VAR::U_SYMAT5] = exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));

}
#endif
/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
//Time to generate ../../Bobby_research/rhs.cpp.in : 5.29 sec
