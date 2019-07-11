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

#if 0
    double vars[24];
    pmin[0] = 0.0;
    pmin[1] = 0.0;
    pmin[2] = 0.0;
    pmax[0] = 2.0;
    pmax[1] = 2.0;
    pmax[2] = 2.0;

    hx = (pmax[0] - pmin[0]) / (nx - 1);
    hy = (pmax[1] - pmin[1]) / (ny - 1);
    hz = (pmax[2] - pmin[2]) / (nz - 1);

    for (unsigned int k = 0; k < nz; k++) {
        double z = pmin[2] + k*hz;
        for (unsigned int j = 0; j < ny; j++) {
            double y = pmin[1] + j*hy;
            for (unsigned int i = 0; i < nx; i++) {
                double x = pmin[0] + i*hx;
                int pp = i + nx*(j + k*ny);
                fake_initial_data(x, y, z, vars);
                for (unsigned int m = 0; m < 24; m++) {
                    uZipVars[m][offset+pp] = vars[m];
                }
            }
        }
    }
#endif

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

    register double x;
    register double y;
    register double z;
    register unsigned int pp;

    double r_coord;
    double eta;

    //cout << "begin loop" << endl;
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }


                bssn::timer::t_rhs.start();

#ifdef USE_ETA_FUNC
                /*[[[cog

                import dendro
                import bssn

                outs = [bssn.a_rhs, bssn.b_rhs, bssn.gt_rhs, bssn.chi_rhs, bssn.At_rhs, bssn.K_rhs, bssn.Gt_rhs, bssn.B_rhs]
                vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']
                #outs=outs[0:3]
                #vnames=vnames[0:3]
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 668203 
                // Dendro: printing temp variables
                double DENDRO_0 = 2*alpha[pp];
                double DENDRO_1 = (3.0/4.0)*alpha[pp]*lambda_f[1] + (3.0/4.0)*lambda_f[0];
                double DENDRO_2 = grad_0_beta0[pp];
                double DENDRO_3 = (4.0/3.0)*DENDRO_2;
                double DENDRO_4 = grad_1_beta1[pp];
                double DENDRO_5 = (2.0/3.0)*gt0[pp];
                double DENDRO_6 = grad_2_beta2[pp];
                double DENDRO_7 = grad_0_beta1[pp];
                double DENDRO_8 = 2*gt1[pp];
                double DENDRO_9 = grad_0_beta2[pp];
                double DENDRO_10 = 2*gt2[pp];
                double DENDRO_11 = grad_1_beta0[pp];
                double DENDRO_12 = grad_1_beta2[pp];
                double DENDRO_13 = (1.0/3.0)*gt1[pp];
                double DENDRO_14 = (2.0/3.0)*DENDRO_6;
                double DENDRO_15 = grad_2_beta0[pp];
                double DENDRO_16 = grad_2_beta1[pp];
                double DENDRO_17 = (1.0/3.0)*gt2[pp];
                double DENDRO_18 = (2.0/3.0)*DENDRO_4;
                double DENDRO_19 = (2.0/3.0)*DENDRO_2;
                double DENDRO_20 = (4.0/3.0)*DENDRO_4;
                double DENDRO_21 = 2*gt4[pp];
                double DENDRO_22 = (1.0/3.0)*gt4[pp];
                double DENDRO_23 = (4.0/3.0)*DENDRO_6;
                double DENDRO_24 = (2.0/3.0)*chi[pp];
                double DENDRO_25 = DENDRO_2 + DENDRO_4 + DENDRO_6;
                double DENDRO_26 = 2*At1[pp];
                double DENDRO_27 = 2*At2[pp];
                double DENDRO_28 = gt2[pp]*gt4[pp];
                double DENDRO_29 = -DENDRO_28 + gt1[pp]*gt5[pp];
                double DENDRO_30 = At0[pp]*DENDRO_29;
                double DENDRO_31 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
                double DENDRO_32 = At2[pp]*DENDRO_31;
                double DENDRO_33 = gt0[pp]*gt5[pp];
                double DENDRO_34 = pow(gt2[pp], 2);
                double DENDRO_35 = DENDRO_33 - DENDRO_34;
                double DENDRO_36 = At1[pp]*DENDRO_35;
                double DENDRO_37 = pow(gt4[pp], 2);
                double DENDRO_38 = pow(gt1[pp], 2);
                double DENDRO_39 = gt3[pp]*gt5[pp];
                double DENDRO_40 = -DENDRO_28*DENDRO_8 + DENDRO_34*gt3[pp] + DENDRO_37*gt0[pp] + DENDRO_38*gt5[pp] - DENDRO_39*gt0[pp];
                double DENDRO_41 = 1.0/DENDRO_40;
                double DENDRO_42 = DENDRO_26*DENDRO_41;
                double DENDRO_43 = At1[pp]*DENDRO_29;
                double DENDRO_44 = gt1[pp]*gt4[pp];
                double DENDRO_45 = gt2[pp]*gt3[pp];
                double DENDRO_46 = DENDRO_44 - DENDRO_45;
                double DENDRO_47 = At2[pp]*DENDRO_46;
                double DENDRO_48 = -DENDRO_47;
                double DENDRO_49 = -DENDRO_37 + DENDRO_39;
                double DENDRO_50 = At0[pp]*DENDRO_49;
                double DENDRO_51 = 2*DENDRO_41;
                double DENDRO_52 = At0[pp]*DENDRO_51;
                double DENDRO_53 = At1[pp]*DENDRO_31;
                double DENDRO_54 = gt0[pp]*gt3[pp];
                double DENDRO_55 = -DENDRO_38 + DENDRO_54;
                double DENDRO_56 = At2[pp]*DENDRO_55;
                double DENDRO_57 = DENDRO_27*DENDRO_41;
                double DENDRO_58 = grad2_0_0_alpha[pp];
                double DENDRO_59 = grad_1_chi[pp];
                double DENDRO_60 = grad_2_chi[pp];
                double DENDRO_61 = grad_0_chi[pp];
                double DENDRO_62 = DENDRO_29*DENDRO_61 + DENDRO_31*DENDRO_60;
                double DENDRO_63 = -DENDRO_35*DENDRO_59 + DENDRO_62;
                double DENDRO_64 = 1.0/chi[pp];
                double DENDRO_65 = 0.5*DENDRO_64;
                double DENDRO_66 = DENDRO_65*gt0[pp];
                double DENDRO_67 = grad_0_gt1[pp];
                double DENDRO_68 = 1.0*DENDRO_67;
                double DENDRO_69 = grad_1_gt0[pp];
                double DENDRO_70 = 0.5*DENDRO_69;
                double DENDRO_71 = DENDRO_68 - DENDRO_70;
                double DENDRO_72 = grad_0_gt2[pp];
                double DENDRO_73 = 1.0*DENDRO_72;
                double DENDRO_74 = grad_2_gt0[pp];
                double DENDRO_75 = 0.5*DENDRO_74;
                double DENDRO_76 = DENDRO_73 - DENDRO_75;
                double DENDRO_77 = grad_0_gt0[pp];
                double DENDRO_78 = 0.5*DENDRO_77;
                double DENDRO_79 = DENDRO_29*DENDRO_78 + DENDRO_31*DENDRO_76;
                double DENDRO_80 = -DENDRO_35*DENDRO_71 + DENDRO_79;
                double DENDRO_81 = DENDRO_63*DENDRO_66 + DENDRO_80;
                double DENDRO_82 = grad_1_alpha[pp];
                double DENDRO_83 = 12*DENDRO_82;
                double DENDRO_84 = DENDRO_41*DENDRO_83;
                double DENDRO_85 = DENDRO_31*DENDRO_59;
                double DENDRO_86 = DENDRO_46*DENDRO_61;
                double DENDRO_87 = DENDRO_55*DENDRO_60;
                double DENDRO_88 = DENDRO_85 - DENDRO_86 - DENDRO_87;
                double DENDRO_89 = DENDRO_66*DENDRO_88;
                double DENDRO_90 = DENDRO_46*DENDRO_78;
                double DENDRO_91 = DENDRO_55*DENDRO_76;
                double DENDRO_92 = DENDRO_31*DENDRO_71;
                double DENDRO_93 = DENDRO_90 + DENDRO_91 - DENDRO_92;
                double DENDRO_94 = grad_2_alpha[pp];
                double DENDRO_95 = 12*DENDRO_94;
                double DENDRO_96 = DENDRO_41*DENDRO_95;
                double DENDRO_97 = DENDRO_49*DENDRO_78;
                double DENDRO_98 = DENDRO_41*DENDRO_97;
                double DENDRO_99 = DENDRO_46*DENDRO_76;
                double DENDRO_100 = DENDRO_41*DENDRO_99;
                double DENDRO_101 = DENDRO_29*DENDRO_71;
                double DENDRO_102 = DENDRO_101*DENDRO_41;
                double DENDRO_103 = 1.0*DENDRO_61;
                double DENDRO_104 = DENDRO_41*gt0[pp];
                double DENDRO_105 = -DENDRO_44 + DENDRO_45;
                double DENDRO_106 = DENDRO_29*DENDRO_59;
                double DENDRO_107 = DENDRO_37 - DENDRO_39;
                double DENDRO_108 = DENDRO_105*DENDRO_60 + DENDRO_106 + DENDRO_107*DENDRO_61;
                double DENDRO_109 = 0.5*DENDRO_108;
                double DENDRO_110 = DENDRO_104*DENDRO_109;
                double DENDRO_111 = DENDRO_64*(-DENDRO_103 + DENDRO_110);
                double DENDRO_112 = grad_0_alpha[pp];
                double DENDRO_113 = 12*DENDRO_112;
                double DENDRO_114 = grad2_0_0_chi[pp];
                double DENDRO_115 = -DENDRO_114;
                double DENDRO_116 = -DENDRO_90 - DENDRO_91 + DENDRO_92;
                double DENDRO_117 = DENDRO_41*DENDRO_60;
                double DENDRO_118 = DENDRO_41*DENDRO_59;
                double DENDRO_119 = DENDRO_101 - DENDRO_97 - DENDRO_99;
                double DENDRO_120 = DENDRO_41*DENDRO_61;
                double DENDRO_121 = 2*DENDRO_64;
                double DENDRO_122 = grad_0_gt3[pp];
                double DENDRO_123 = 0.5*DENDRO_122;
                double DENDRO_124 = grad_1_gt1[pp];
                double DENDRO_125 = 1.0*DENDRO_124;
                double DENDRO_126 = -DENDRO_125;
                double DENDRO_127 = DENDRO_123 + DENDRO_126;
                double DENDRO_128 = grad_0_gt4[pp];
                double DENDRO_129 = grad_2_gt1[pp];
                double DENDRO_130 = grad_1_gt2[pp];
                double DENDRO_131 = DENDRO_128 + DENDRO_129 - DENDRO_130;
                double DENDRO_132 = grad_0_gt5[pp];
                double DENDRO_133 = DENDRO_132*DENDRO_31 + DENDRO_29*DENDRO_74;
                double DENDRO_134 = -DENDRO_131*DENDRO_35 + DENDRO_133;
                double DENDRO_135 = DENDRO_127*DENDRO_134;
                double DENDRO_136 = DENDRO_122*DENDRO_35;
                double DENDRO_137 = DENDRO_128 - DENDRO_129 + DENDRO_130;
                double DENDRO_138 = DENDRO_137*DENDRO_31;
                double DENDRO_139 = DENDRO_29*DENDRO_69;
                double DENDRO_140 = DENDRO_138 + DENDRO_139;
                double DENDRO_141 = -DENDRO_136 + DENDRO_140;
                double DENDRO_142 = DENDRO_131*DENDRO_141;
                double DENDRO_143 = 0.25*DENDRO_142;
                double DENDRO_144 = pow(DENDRO_40, -2);
                double DENDRO_145 = DENDRO_144*DENDRO_31;
                double DENDRO_146 = 4*DENDRO_145;
                double DENDRO_147 = 0.5*DENDRO_132;
                double DENDRO_148 = grad_2_gt2[pp];
                double DENDRO_149 = 1.0*DENDRO_148;
                double DENDRO_150 = -DENDRO_149;
                double DENDRO_151 = DENDRO_147 + DENDRO_150;
                double DENDRO_152 = DENDRO_122*DENDRO_31;
                double DENDRO_153 = DENDRO_46*DENDRO_69;
                double DENDRO_154 = DENDRO_137*DENDRO_55;
                double DENDRO_155 = DENDRO_152 - DENDRO_153 - DENDRO_154;
                double DENDRO_156 = DENDRO_151*DENDRO_155;
                double DENDRO_157 = DENDRO_46*DENDRO_74;
                double DENDRO_158 = DENDRO_132*DENDRO_55;
                double DENDRO_159 = DENDRO_131*DENDRO_31;
                double DENDRO_160 = -DENDRO_157 - DENDRO_158 + DENDRO_159;
                double DENDRO_161 = DENDRO_137*DENDRO_160;
                double DENDRO_162 = 0.25*DENDRO_161;
                double DENDRO_163 = -DENDRO_162;
                double DENDRO_164 = -DENDRO_152 + DENDRO_153 + DENDRO_154;
                double DENDRO_165 = 0.25*DENDRO_132;
                double DENDRO_166 = DENDRO_164*DENDRO_165;
                double DENDRO_167 = -DENDRO_128 + DENDRO_129 + DENDRO_130;
                double DENDRO_168 = DENDRO_157 + DENDRO_158 - DENDRO_159;
                double DENDRO_169 = 0.5*DENDRO_168;
                double DENDRO_170 = DENDRO_122*DENDRO_29;
                double DENDRO_171 = DENDRO_49*DENDRO_69;
                double DENDRO_172 = DENDRO_137*DENDRO_46;
                double DENDRO_173 = -DENDRO_170 + DENDRO_171 + DENDRO_172;
                double DENDRO_174 = 0.25*DENDRO_74;
                double DENDRO_175 = DENDRO_173*DENDRO_174;
                double DENDRO_176 = DENDRO_132*DENDRO_46;
                double DENDRO_177 = DENDRO_49*DENDRO_74;
                double DENDRO_178 = DENDRO_131*DENDRO_29;
                double DENDRO_179 = DENDRO_176 + DENDRO_177 - DENDRO_178;
                double DENDRO_180 = 0.25*DENDRO_69;
                double DENDRO_181 = DENDRO_179*DENDRO_180;
                double DENDRO_182 = 2*DENDRO_151;
                double DENDRO_183 = DENDRO_144*DENDRO_46;
                double DENDRO_184 = 4*DENDRO_183;
                double DENDRO_185 = DENDRO_173*DENDRO_77;
                double DENDRO_186 = 0.25*DENDRO_185;
                double DENDRO_187 = -DENDRO_101 + DENDRO_97 + DENDRO_99;
                double DENDRO_188 = DENDRO_187*DENDRO_69;
                double DENDRO_189 = DENDRO_144*DENDRO_29;
                double DENDRO_190 = 4*DENDRO_189;
                double DENDRO_191 = 0.5*DENDRO_76;
                double DENDRO_192 = DENDRO_179*DENDRO_77;
                double DENDRO_193 = 0.25*DENDRO_192;
                double DENDRO_194 = DENDRO_187*DENDRO_74;
                double DENDRO_195 = DENDRO_168*DENDRO_74;
                double DENDRO_196 = DENDRO_132*DENDRO_93;
                double DENDRO_197 = 2.0*DENDRO_183;
                double DENDRO_198 = DENDRO_164*DENDRO_74;
                double DENDRO_199 = DENDRO_137*DENDRO_93;
                double DENDRO_200 = 2.0*DENDRO_189;
                double DENDRO_201 = DENDRO_134*DENDRO_46;
                double DENDRO_202 = grad_2_gt3[pp];
                double DENDRO_203 = DENDRO_202*DENDRO_35;
                double DENDRO_204 = grad_1_gt5[pp];
                double DENDRO_205 = DENDRO_204*DENDRO_31;
                double DENDRO_206 = DENDRO_167*DENDRO_29;
                double DENDRO_207 = DENDRO_205 + DENDRO_206;
                double DENDRO_208 = -DENDRO_203 + DENDRO_207;
                double DENDRO_209 = DENDRO_208*DENDRO_31;
                double DENDRO_210 = DENDRO_141*DENDRO_29;
                double DENDRO_211 = grad_2_gt5[pp];
                double DENDRO_212 = 0.5*DENDRO_211;
                double DENDRO_213 = DENDRO_212*DENDRO_31;
                double DENDRO_214 = 0.5*DENDRO_204;
                double DENDRO_215 = grad_2_gt4[pp];
                double DENDRO_216 = 1.0*DENDRO_215;
                double DENDRO_217 = -DENDRO_216;
                double DENDRO_218 = DENDRO_214 + DENDRO_217;
                double DENDRO_219 = -DENDRO_151*DENDRO_29 + DENDRO_213 + DENDRO_218*DENDRO_35;
                double DENDRO_220 = DENDRO_219*DENDRO_55;
                double DENDRO_221 = grad_1_gt3[pp];
                double DENDRO_222 = 0.5*DENDRO_221;
                double DENDRO_223 = DENDRO_222*DENDRO_35;
                double DENDRO_224 = grad_1_gt4[pp];
                double DENDRO_225 = 1.0*DENDRO_224;
                double DENDRO_226 = 0.5*DENDRO_202;
                double DENDRO_227 = DENDRO_225 - DENDRO_226;
                double DENDRO_228 = DENDRO_227*DENDRO_31;
                double DENDRO_229 = DENDRO_127*DENDRO_29;
                double DENDRO_230 = -DENDRO_223 + DENDRO_228 - DENDRO_229;
                double DENDRO_231 = DENDRO_230*DENDRO_35;
                double DENDRO_232 = DENDRO_49*DENDRO_80;
                double DENDRO_233 = DENDRO_144*(DENDRO_201 - 1.0*DENDRO_209 - 1.0*DENDRO_210 + DENDRO_220 + DENDRO_231 + DENDRO_232);
                double DENDRO_234 = 2.0*DENDRO_69;
                double DENDRO_235 = -DENDRO_176 - DENDRO_177 + DENDRO_178;
                double DENDRO_236 = DENDRO_235*DENDRO_46;
                double DENDRO_237 = DENDRO_202*DENDRO_29;
                double DENDRO_238 = DENDRO_204*DENDRO_46;
                double DENDRO_239 = DENDRO_167*DENDRO_49;
                double DENDRO_240 = DENDRO_237 - DENDRO_238 - DENDRO_239;
                double DENDRO_241 = DENDRO_240*DENDRO_31;
                double DENDRO_242 = DENDRO_170 - DENDRO_171 - DENDRO_172;
                double DENDRO_243 = DENDRO_242*DENDRO_29;
                double DENDRO_244 = DENDRO_212*DENDRO_46;
                double DENDRO_245 = DENDRO_218*DENDRO_29;
                double DENDRO_246 = DENDRO_151*DENDRO_49;
                double DENDRO_247 = -DENDRO_244 - DENDRO_245 + DENDRO_246;
                double DENDRO_248 = DENDRO_247*DENDRO_55;
                double DENDRO_249 = DENDRO_222*DENDRO_29;
                double DENDRO_250 = DENDRO_127*DENDRO_49 - DENDRO_227*DENDRO_46 + DENDRO_249;
                double DENDRO_251 = DENDRO_250*DENDRO_35;
                double DENDRO_252 = DENDRO_119*DENDRO_49;
                double DENDRO_253 = DENDRO_144*(DENDRO_236 - 1.0*DENDRO_241 - 1.0*DENDRO_243 + DENDRO_248 + DENDRO_251 + DENDRO_252);
                double DENDRO_254 = 2.0*DENDRO_77;
                double DENDRO_255 = DENDRO_160*DENDRO_46;
                double DENDRO_256 = DENDRO_202*DENDRO_31;
                double DENDRO_257 = DENDRO_204*DENDRO_55;
                double DENDRO_258 = DENDRO_167*DENDRO_46;
                double DENDRO_259 = DENDRO_256 - DENDRO_257 - DENDRO_258;
                double DENDRO_260 = DENDRO_259*DENDRO_31;
                double DENDRO_261 = DENDRO_155*DENDRO_29;
                double DENDRO_262 = DENDRO_212*DENDRO_55;
                double DENDRO_263 = DENDRO_151*DENDRO_46;
                double DENDRO_264 = DENDRO_218*DENDRO_31;
                double DENDRO_265 = -DENDRO_262 + DENDRO_263 - DENDRO_264;
                double DENDRO_266 = DENDRO_265*DENDRO_55;
                double DENDRO_267 = DENDRO_222*DENDRO_31;
                double DENDRO_268 = DENDRO_127*DENDRO_46 - DENDRO_227*DENDRO_55 + DENDRO_267;
                double DENDRO_269 = DENDRO_268*DENDRO_35;
                double DENDRO_270 = DENDRO_116*DENDRO_49;
                double DENDRO_271 = DENDRO_144*(DENDRO_255 - 1.0*DENDRO_260 - 1.0*DENDRO_261 + DENDRO_266 + DENDRO_269 + DENDRO_270);
                double DENDRO_272 = 2.0*DENDRO_74;
                double DENDRO_273 = grad2_2_2_chi[pp];
                double DENDRO_274 = pow(DENDRO_60, 2);
                double DENDRO_275 = 3*DENDRO_64;
                double DENDRO_276 = DENDRO_55*(2*DENDRO_273 - DENDRO_274*DENDRO_275);
                double DENDRO_277 = grad2_1_1_chi[pp];
                double DENDRO_278 = pow(DENDRO_59, 2);
                double DENDRO_279 = DENDRO_35*(-DENDRO_275*DENDRO_278 + 2*DENDRO_277);
                double DENDRO_280 = pow(DENDRO_61, 2);
                double DENDRO_281 = DENDRO_49*(2*DENDRO_114 - DENDRO_275*DENDRO_280);
                double DENDRO_282 = grad2_1_2_chi[pp];
                double DENDRO_283 = DENDRO_275*DENDRO_60;
                double DENDRO_284 = 2*DENDRO_31;
                double DENDRO_285 = DENDRO_284*(2*DENDRO_282 - DENDRO_283*DENDRO_59);
                double DENDRO_286 = grad2_0_2_chi[pp];
                double DENDRO_287 = 2*DENDRO_46;
                double DENDRO_288 = DENDRO_287*(-DENDRO_283*DENDRO_61 + 2*DENDRO_286);
                double DENDRO_289 = grad2_0_1_chi[pp];
                double DENDRO_290 = DENDRO_59*DENDRO_61;
                double DENDRO_291 = 2*DENDRO_29;
                double DENDRO_292 = DENDRO_291*(-DENDRO_275*DENDRO_290 + 2*DENDRO_289);
                double DENDRO_293 = -1.0*DENDRO_201 + DENDRO_209 + DENDRO_210 - DENDRO_220 - DENDRO_231 - DENDRO_232;
                double DENDRO_294 = 2*DENDRO_118*DENDRO_293;
                double DENDRO_295 = -1.0*DENDRO_255 + DENDRO_260 + DENDRO_261 - DENDRO_266 - DENDRO_269 - DENDRO_270;
                double DENDRO_296 = 2*DENDRO_117*DENDRO_295;
                double DENDRO_297 = -1.0*DENDRO_236 + DENDRO_241 + DENDRO_243 - DENDRO_248 - DENDRO_251 - DENDRO_252;
                double DENDRO_298 = 2*DENDRO_120*DENDRO_297;
                double DENDRO_299 = DENDRO_276 + DENDRO_279 + DENDRO_281 - DENDRO_285 + DENDRO_288 - DENDRO_292 + DENDRO_294 + DENDRO_296 + DENDRO_298;
                double DENDRO_300 = DENDRO_104*DENDRO_64;
                double DENDRO_301 = DENDRO_144*DENDRO_55;
                double DENDRO_302 = 4*DENDRO_301;
                double DENDRO_303 = DENDRO_160*DENDRO_302;
                double DENDRO_304 = 0.25*DENDRO_122;
                double DENDRO_305 = DENDRO_144*DENDRO_35;
                double DENDRO_306 = 4*DENDRO_305;
                double DENDRO_307 = DENDRO_141*DENDRO_306;
                double DENDRO_308 = 0.25*DENDRO_128;
                double DENDRO_309 = -DENDRO_308;
                double DENDRO_310 = 0.75*DENDRO_130;
                double DENDRO_311 = 0.25*DENDRO_129;
                double DENDRO_312 = DENDRO_306*(DENDRO_309 + DENDRO_310 + DENDRO_311);
                double DENDRO_313 = DENDRO_73 + DENDRO_75;
                double DENDRO_314 = DENDRO_144*DENDRO_49;
                double DENDRO_315 = 4*DENDRO_314;
                double DENDRO_316 = DENDRO_179*DENDRO_74;
                double DENDRO_317 = 3.0*DENDRO_301;
                double DENDRO_318 = DENDRO_173*DENDRO_69;
                double DENDRO_319 = 3.0*DENDRO_305;
                double DENDRO_320 = 6.0*DENDRO_77;
                double DENDRO_321 = pow(chi[pp], -2);
                double DENDRO_322 = grad_0_Gt0[pp];
                double DENDRO_323 = grad_0_Gt1[pp];
                double DENDRO_324 = 4*gt1[pp];
                double DENDRO_325 = grad_0_Gt2[pp];
                double DENDRO_326 = 4*gt2[pp];
                double DENDRO_327 = 0.5*DENDRO_71;
                double DENDRO_328 = DENDRO_134*DENDRO_327;
                double DENDRO_329 = DENDRO_41*DENDRO_46;
                double DENDRO_330 = 4*DENDRO_329;
                double DENDRO_331 = DENDRO_134*DENDRO_304;
                double DENDRO_332 = 0.5*DENDRO_167;
                double DENDRO_333 = DENDRO_141*DENDRO_327;
                double DENDRO_334 = DENDRO_41*DENDRO_55;
                double DENDRO_335 = 2.0*DENDRO_334;
                double DENDRO_336 = DENDRO_35*DENDRO_41;
                double DENDRO_337 = 2.0*DENDRO_336;
                double DENDRO_338 = DENDRO_41*DENDRO_49;
                double DENDRO_339 = 2.0*DENDRO_338;
                double DENDRO_340 = DENDRO_141*DENDRO_69;
                double DENDRO_341 = DENDRO_122*DENDRO_80;
                double DENDRO_342 = DENDRO_134*DENDRO_69;
                double DENDRO_343 = DENDRO_131*DENDRO_80;
                double DENDRO_344 = 4.0*DENDRO_41;
                double DENDRO_345 = DENDRO_31*DENDRO_344;
                double DENDRO_346 = DENDRO_29*DENDRO_344;
                double DENDRO_347 = 0.25*DENDRO_130;
                double DENDRO_348 = 0.75*DENDRO_129;
                double DENDRO_349 = 4*DENDRO_144;
                double DENDRO_350 = -DENDRO_134*DENDRO_302*(DENDRO_309 + DENDRO_347 + DENDRO_348) + DENDRO_146*(DENDRO_141*DENDRO_332 + DENDRO_331) - DENDRO_184*(DENDRO_167*DENDRO_80 + DENDRO_328) + DENDRO_190*(-2*DENDRO_127*DENDRO_80 + DENDRO_333) - DENDRO_197*(DENDRO_342 + DENDRO_343) + DENDRO_200*(DENDRO_340 + DENDRO_341) - DENDRO_232*DENDRO_349*(DENDRO_68 + DENDRO_70) - DENDRO_280*DENDRO_321 + 4*DENDRO_322*gt0[pp] + DENDRO_323*DENDRO_324 + DENDRO_325*DENDRO_326 + DENDRO_330*grad2_0_2_gt0[pp] + DENDRO_335*grad2_2_2_gt0[pp] + DENDRO_337*grad2_1_1_gt0[pp] + DENDRO_339*grad2_0_0_gt0[pp] - DENDRO_345*grad2_1_2_gt0[pp] - DENDRO_346*grad2_0_1_gt0[pp];
                double DENDRO_351 = 3*alpha[pp];
                double DENDRO_352 = grad2_2_2_alpha[pp];
                double DENDRO_353 = DENDRO_65*gt5[pp];
                double DENDRO_354 = DENDRO_219 + DENDRO_353*DENDRO_63;
                double DENDRO_355 = 4*DENDRO_82;
                double DENDRO_356 = DENDRO_355*DENDRO_41;
                double DENDRO_357 = DENDRO_46*DENDRO_60;
                double DENDRO_358 = DENDRO_49*DENDRO_61;
                double DENDRO_359 = DENDRO_106 - DENDRO_357 - DENDRO_358;
                double DENDRO_360 = DENDRO_353*DENDRO_359;
                double DENDRO_361 = 4*DENDRO_112;
                double DENDRO_362 = DENDRO_361*DENDRO_41;
                double DENDRO_363 = DENDRO_262*DENDRO_41;
                double DENDRO_364 = DENDRO_263*DENDRO_41;
                double DENDRO_365 = DENDRO_264*DENDRO_41;
                double DENDRO_366 = 1.0*DENDRO_60;
                double DENDRO_367 = DENDRO_41*gt5[pp];
                double DENDRO_368 = DENDRO_38 - DENDRO_54;
                double DENDRO_369 = DENDRO_105*DENDRO_61 + DENDRO_368*DENDRO_60 + DENDRO_85;
                double DENDRO_370 = 0.5*DENDRO_369;
                double DENDRO_371 = DENDRO_367*DENDRO_370;
                double DENDRO_372 = DENDRO_64*(-DENDRO_366 + DENDRO_371);
                double DENDRO_373 = 4*DENDRO_94;
                double DENDRO_374 = -DENDRO_273;
                double DENDRO_375 = DENDRO_212*DENDRO_368;
                double DENDRO_376 = -DENDRO_214;
                double DENDRO_377 = DENDRO_216 + DENDRO_376;
                double DENDRO_378 = DENDRO_31*DENDRO_377;
                double DENDRO_379 = -DENDRO_147;
                double DENDRO_380 = DENDRO_149 + DENDRO_379;
                double DENDRO_381 = DENDRO_105*DENDRO_380;
                double DENDRO_382 = -DENDRO_33 + DENDRO_34;
                double DENDRO_383 = DENDRO_213 + DENDRO_29*DENDRO_380 + DENDRO_377*DENDRO_382;
                double DENDRO_384 = DENDRO_105*DENDRO_212 + DENDRO_107*DENDRO_380 + DENDRO_29*DENDRO_377;
                double DENDRO_385 = 0.5*DENDRO_218;
                double DENDRO_386 = DENDRO_134*DENDRO_385;
                double DENDRO_387 = -DENDRO_386;
                double DENDRO_388 = DENDRO_137*DENDRO_219;
                double DENDRO_389 = 0.5*DENDRO_151;
                double DENDRO_390 = -DENDRO_235*DENDRO_389;
                double DENDRO_391 = 2*DENDRO_76;
                double DENDRO_392 = DENDRO_160*DENDRO_211;
                double DENDRO_393 = 0.25*DENDRO_392;
                double DENDRO_394 = DENDRO_132*DENDRO_265;
                double DENDRO_395 = DENDRO_240*DENDRO_389;
                double DENDRO_396 = -DENDRO_395;
                double DENDRO_397 = DENDRO_137*DENDRO_247;
                double DENDRO_398 = DENDRO_211*DENDRO_259;
                double DENDRO_399 = 0.25*DENDRO_398;
                double DENDRO_400 = DENDRO_204*DENDRO_265;
                double DENDRO_401 = DENDRO_240*DENDRO_76;
                double DENDRO_402 = DENDRO_167*DENDRO_235;
                double DENDRO_403 = 0.25*DENDRO_402;
                double DENDRO_404 = 0.25*DENDRO_204;
                double DENDRO_405 = DENDRO_160*DENDRO_404;
                double DENDRO_406 = DENDRO_165*DENDRO_259;
                double DENDRO_407 = DENDRO_174*DENDRO_240;
                double DENDRO_408 = 0.5*DENDRO_137;
                double DENDRO_409 = DENDRO_167*DENDRO_247;
                double DENDRO_410 = 2.0*DENDRO_145;
                double DENDRO_411 = DENDRO_144*DENDRO_293;
                double DENDRO_412 = 2.0*DENDRO_411;
                double DENDRO_413 = DENDRO_144*DENDRO_295;
                double DENDRO_414 = 2.0*DENDRO_413;
                double DENDRO_415 = DENDRO_144*DENDRO_297;
                double DENDRO_416 = 2.0*DENDRO_415;
                double DENDRO_417 = DENDRO_247*DENDRO_74;
                double DENDRO_418 = -DENDRO_276 - DENDRO_279 - DENDRO_281 + DENDRO_285 - DENDRO_288 + DENDRO_292 - DENDRO_294 - DENDRO_296 - DENDRO_298;
                double DENDRO_419 = DENDRO_418*DENDRO_64;
                double DENDRO_420 = DENDRO_248*DENDRO_349;
                double DENDRO_421 = DENDRO_220*DENDRO_349;
                double DENDRO_422 = -DENDRO_311;
                double DENDRO_423 = DENDRO_306*(DENDRO_308 + DENDRO_310 + DENDRO_422);
                double DENDRO_424 = DENDRO_315*(-DENDRO_174 + DENDRO_73);
                double DENDRO_425 = DENDRO_204*DENDRO_259;
                double DENDRO_426 = DENDRO_132*DENDRO_160;
                double DENDRO_427 = 3.0*DENDRO_314;
                double DENDRO_428 = 6.0*DENDRO_144;
                double DENDRO_429 = grad_2_Gt0[pp];
                double DENDRO_430 = grad_2_Gt1[pp];
                double DENDRO_431 = 4*gt4[pp];
                double DENDRO_432 = grad_2_Gt2[pp];
                double DENDRO_433 = -DENDRO_208*DENDRO_385;
                double DENDRO_434 = DENDRO_134*DENDRO_227;
                double DENDRO_435 = DENDRO_131*DENDRO_208;
                double DENDRO_436 = 0.25*DENDRO_435;
                double DENDRO_437 = 0.25*DENDRO_202;
                double DENDRO_438 = DENDRO_134*DENDRO_437;
                double DENDRO_439 = DENDRO_204*DENDRO_208;
                double DENDRO_440 = DENDRO_202*DENDRO_219;
                double DENDRO_441 = DENDRO_134*DENDRO_204;
                double DENDRO_442 = DENDRO_131*DENDRO_219;
                double DENDRO_443 = 0.75*DENDRO_128;
                double DENDRO_444 = -DENDRO_134*DENDRO_315*(DENDRO_347 + DENDRO_422 + DENDRO_443) + DENDRO_146*(2*DENDRO_219*DENDRO_227 + DENDRO_433) + DENDRO_190*(DENDRO_434 + DENDRO_436) + DENDRO_190*(DENDRO_208*DENDRO_408 + DENDRO_438) - DENDRO_197*(DENDRO_441 + DENDRO_442) - DENDRO_208*DENDRO_306*(DENDRO_225 - DENDRO_437) - DENDRO_274*DENDRO_321 + DENDRO_326*DENDRO_429 + DENDRO_330*grad2_0_2_gt5[pp] + DENDRO_335*grad2_2_2_gt5[pp] + DENDRO_337*grad2_1_1_gt5[pp] + DENDRO_339*grad2_0_0_gt5[pp] - DENDRO_345*grad2_1_2_gt5[pp] - DENDRO_346*grad2_0_1_gt5[pp] + DENDRO_410*(DENDRO_439 + DENDRO_440) + DENDRO_430*DENDRO_431 + 4*DENDRO_432*gt5[pp];
                double DENDRO_445 = grad2_1_1_alpha[pp];
                double DENDRO_446 = DENDRO_65*gt3[pp];
                double DENDRO_447 = DENDRO_373*DENDRO_41;
                double DENDRO_448 = 1.0*DENDRO_59;
                double DENDRO_449 = -DENDRO_448;
                double DENDRO_450 = DENDRO_382*DENDRO_59 + DENDRO_62;
                double DENDRO_451 = DENDRO_41*gt3[pp];
                double DENDRO_452 = 0.5*DENDRO_451;
                double DENDRO_453 = DENDRO_450*DENDRO_452;
                double DENDRO_454 = DENDRO_228*DENDRO_41;
                double DENDRO_455 = -DENDRO_223*DENDRO_41 - DENDRO_229*DENDRO_41 + DENDRO_454;
                double DENDRO_456 = -DENDRO_277;
                double DENDRO_457 = -DENDRO_123;
                double DENDRO_458 = DENDRO_125 + DENDRO_457;
                double DENDRO_459 = DENDRO_105*DENDRO_458 + DENDRO_227*DENDRO_368 + DENDRO_267;
                double DENDRO_460 = DENDRO_222*DENDRO_382;
                double DENDRO_461 = DENDRO_29*DENDRO_458;
                double DENDRO_462 = DENDRO_105*DENDRO_227 + DENDRO_107*DENDRO_458 + DENDRO_249;
                double DENDRO_463 = DENDRO_240*DENDRO_71;
                double DENDRO_464 = DENDRO_167*DENDRO_242;
                double DENDRO_465 = 0.25*DENDRO_464;
                double DENDRO_466 = DENDRO_155*DENDRO_218;
                double DENDRO_467 = DENDRO_137*DENDRO_259;
                double DENDRO_468 = 0.25*DENDRO_467;
                double DENDRO_469 = DENDRO_155*DENDRO_404;
                double DENDRO_470 = 0.5*DENDRO_131;
                double DENDRO_471 = DENDRO_180*DENDRO_240;
                double DENDRO_472 = 0.5*DENDRO_127;
                double DENDRO_473 = DENDRO_240*DENDRO_472;
                double DENDRO_474 = -DENDRO_473;
                double DENDRO_475 = DENDRO_131*DENDRO_250;
                double DENDRO_476 = 0.5*DENDRO_227;
                double DENDRO_477 = DENDRO_259*DENDRO_476;
                double DENDRO_478 = 2*DENDRO_218*DENDRO_268;
                double DENDRO_479 = DENDRO_208*DENDRO_221;
                double DENDRO_480 = 0.25*DENDRO_479;
                double DENDRO_481 = DENDRO_202*DENDRO_230;
                double DENDRO_482 = DENDRO_155*DENDRO_476;
                double DENDRO_483 = DENDRO_131*DENDRO_268;
                double DENDRO_484 = -DENDRO_242*DENDRO_472;
                double DENDRO_485 = 2*DENDRO_250*DENDRO_71;
                double DENDRO_486 = DENDRO_141*DENDRO_221;
                double DENDRO_487 = 0.25*DENDRO_486;
                double DENDRO_488 = DENDRO_122*DENDRO_230;
                double DENDRO_489 = DENDRO_167*DENDRO_250;
                double DENDRO_490 = DENDRO_204*DENDRO_268;
                double DENDRO_491 = DENDRO_137*DENDRO_268;
                double DENDRO_492 = DENDRO_250*DENDRO_69;
                double DENDRO_493 = DENDRO_259*DENDRO_302;
                double DENDRO_494 = -DENDRO_347;
                double DENDRO_495 = DENDRO_302*(DENDRO_308 + DENDRO_348 + DENDRO_494);
                double DENDRO_496 = DENDRO_251*DENDRO_349;
                double DENDRO_497 = DENDRO_315*(-DENDRO_180 + DENDRO_68);
                double DENDRO_498 = DENDRO_315*(DENDRO_311 + DENDRO_443 + DENDRO_494);
                double DENDRO_499 = grad_1_Gt0[pp];
                double DENDRO_500 = grad_1_Gt1[pp];
                double DENDRO_501 = grad_1_Gt2[pp];
                double DENDRO_502 = DENDRO_141*DENDRO_437;
                double DENDRO_503 = DENDRO_208*DENDRO_304;
                double DENDRO_504 = DENDRO_202*DENDRO_208;
                double DENDRO_505 = DENDRO_122*DENDRO_141;
                double DENDRO_506 = -DENDRO_184*(DENDRO_123*DENDRO_208 + DENDRO_502) - DENDRO_184*(DENDRO_141*DENDRO_226 + DENDRO_503) - DENDRO_269*DENDRO_349*(DENDRO_225 + DENDRO_226) - DENDRO_278*DENDRO_321 - DENDRO_317*DENDRO_504 + DENDRO_324*DENDRO_499 + DENDRO_330*grad2_0_2_gt3[pp] + DENDRO_335*grad2_2_2_gt3[pp] + DENDRO_337*grad2_1_1_gt3[pp] + DENDRO_339*grad2_0_0_gt3[pp] - DENDRO_345*grad2_1_2_gt3[pp] - DENDRO_346*grad2_0_1_gt3[pp] - DENDRO_427*DENDRO_505 + DENDRO_431*DENDRO_501 + 4*DENDRO_500*gt3[pp];
                double DENDRO_507 = DENDRO_105*DENDRO_78 + DENDRO_368*DENDRO_76 + DENDRO_92;
                double DENDRO_508 = DENDRO_382*DENDRO_71 + DENDRO_79;
                double DENDRO_509 = DENDRO_107*DENDRO_78;
                double DENDRO_510 = DENDRO_105*DENDRO_76;
                double DENDRO_511 = DENDRO_160*DENDRO_191;
                double DENDRO_512 = DENDRO_235*DENDRO_77;
                double DENDRO_513 = 0.25*DENDRO_512;
                double DENDRO_514 = DENDRO_119*DENDRO_74;
                double DENDRO_515 = DENDRO_155*DENDRO_165;
                double DENDRO_516 = DENDRO_180*DENDRO_235;
                double DENDRO_517 = DENDRO_174*DENDRO_242;
                double DENDRO_518 = DENDRO_155*DENDRO_191;
                double DENDRO_519 = DENDRO_242*DENDRO_77;
                double DENDRO_520 = 0.25*DENDRO_519;
                double DENDRO_521 = DENDRO_119*DENDRO_69;
                double DENDRO_522 = DENDRO_116*DENDRO_137;
                double DENDRO_523 = DENDRO_116*DENDRO_132;
                double DENDRO_524 = DENDRO_235*DENDRO_74;
                double DENDRO_525 = DENDRO_242*DENDRO_69;
                double DENDRO_526 = DENDRO_259*DENDRO_74;
                double DENDRO_527 = DENDRO_132*DENDRO_155;
                double DENDRO_528 = DENDRO_161 + DENDRO_527;
                double DENDRO_529 = DENDRO_208*DENDRO_69;
                double DENDRO_530 = DENDRO_137*DENDRO_141;
                double DENDRO_531 = DENDRO_122*DENDRO_134;
                double DENDRO_532 = DENDRO_530 + DENDRO_531;
                double DENDRO_533 = DENDRO_165*DENDRO_235;
                double DENDRO_534 = DENDRO_151*DENDRO_265;
                double DENDRO_535 = -DENDRO_534;
                double DENDRO_536 = DENDRO_219*DENDRO_332 + 0.25*DENDRO_441;
                double DENDRO_537 = DENDRO_137*DENDRO_242;
                double DENDRO_538 = 0.25*DENDRO_537;
                double DENDRO_539 = DENDRO_141*DENDRO_476;
                double DENDRO_540 = -DENDRO_208*DENDRO_472;
                double DENDRO_541 = DENDRO_539 + DENDRO_540;
                double DENDRO_542 = DENDRO_119*DENDRO_76;
                double DENDRO_543 = DENDRO_160*DENDRO_174;
                double DENDRO_544 = 0.25*DENDRO_342 + DENDRO_408*DENDRO_80;
                double DENDRO_545 = 0.25*DENDRO_134;
                double DENDRO_546 = DENDRO_167*DENDRO_545 + DENDRO_214*DENDRO_80;
                double DENDRO_547 = DENDRO_116*DENDRO_212;
                double DENDRO_548 = -DENDRO_160*DENDRO_389;
                double DENDRO_549 = DENDRO_131*DENDRO_545;
                double DENDRO_550 = DENDRO_137*DENDRO_545 + DENDRO_219*DENDRO_70;
                double DENDRO_551 = DENDRO_174*DENDRO_235;
                double DENDRO_552 = DENDRO_247*DENDRO_78;
                double DENDRO_553 = DENDRO_191*DENDRO_235 + DENDRO_552;
                double DENDRO_554 = DENDRO_127*DENDRO_219;
                double DENDRO_555 = 0.5*DENDRO_434;
                double DENDRO_556 = -DENDRO_554 + DENDRO_555;
                double DENDRO_557 = DENDRO_167*DENDRO_208;
                double DENDRO_558 = 0.25*DENDRO_557;
                double DENDRO_559 = DENDRO_123*DENDRO_219;
                double DENDRO_560 = DENDRO_141*DENDRO_404 + DENDRO_559;
                double DENDRO_561 = 0.25*DENDRO_137;
                double DENDRO_562 = DENDRO_247*DENDRO_70;
                double DENDRO_563 = DENDRO_235*DENDRO_561 + DENDRO_562;
                double DENDRO_564 = DENDRO_265*DENDRO_332;
                double DENDRO_565 = DENDRO_405 + DENDRO_406;
                double DENDRO_566 = DENDRO_259*DENDRO_389;
                double DENDRO_567 = -DENDRO_566;
                double DENDRO_568 = 0.25*DENDRO_155;
                double DENDRO_569 = DENDRO_211*DENDRO_568;
                double DENDRO_570 = DENDRO_265*DENDRO_408 + DENDRO_569;
                double DENDRO_571 = DENDRO_165*DENDRO_242 + DENDRO_562;
                double DENDRO_572 = DENDRO_208*DENDRO_327;
                double DENDRO_573 = -0.5*DENDRO_135 + DENDRO_227*DENDRO_80;
                double DENDRO_574 = DENDRO_116*DENDRO_214;
                double DENDRO_575 = DENDRO_191*DENDRO_259;
                double DENDRO_576 = 0.25*DENDRO_167;
                double DENDRO_577 = DENDRO_160*DENDRO_576;
                double DENDRO_578 = 0.25*DENDRO_240;
                double DENDRO_579 = DENDRO_578*DENDRO_77;
                double DENDRO_580 = DENDRO_191*DENDRO_242;
                double DENDRO_581 = DENDRO_579 + DENDRO_580;
                double DENDRO_582 = DENDRO_119*DENDRO_408 + DENDRO_516;
                double DENDRO_583 = DENDRO_167*DENDRO_259;
                double DENDRO_584 = DENDRO_155*DENDRO_204;
                double DENDRO_585 = DENDRO_467 + DENDRO_584;
                double DENDRO_586 = 1.0*DENDRO_305;
                double DENDRO_587 = -DENDRO_286;
                double DENDRO_588 = DENDRO_132*DENDRO_368;
                double DENDRO_589 = DENDRO_105*DENDRO_74;
                double DENDRO_590 = 0.5*DENDRO_117;
                double DENDRO_591 = DENDRO_131*DENDRO_382 + DENDRO_133;
                double DENDRO_592 = 0.5*DENDRO_118;
                double DENDRO_593 = DENDRO_105*DENDRO_132;
                double DENDRO_594 = DENDRO_107*DENDRO_74;
                double DENDRO_595 = 0.5*DENDRO_120;
                double DENDRO_596 = 2.0*DENDRO_429;
                double DENDRO_597 = DENDRO_596*gt0[pp];
                double DENDRO_598 = 2.0*DENDRO_430;
                double DENDRO_599 = DENDRO_598*gt1[pp];
                double DENDRO_600 = 2.0*gt2[pp];
                double DENDRO_601 = DENDRO_322*DENDRO_600;
                double DENDRO_602 = DENDRO_432*DENDRO_600;
                double DENDRO_603 = 2.0*gt4[pp];
                double DENDRO_604 = DENDRO_323*DENDRO_603;
                double DENDRO_605 = 2.0*gt5[pp];
                double DENDRO_606 = DENDRO_325*DENDRO_605;
                double DENDRO_607 = DENDRO_321*DENDRO_60;
                double DENDRO_608 = -DENDRO_607*DENDRO_61;
                double DENDRO_609 = DENDRO_330*grad2_0_2_gt2[pp];
                double DENDRO_610 = DENDRO_335*grad2_2_2_gt2[pp];
                double DENDRO_611 = DENDRO_337*grad2_1_1_gt2[pp];
                double DENDRO_612 = DENDRO_339*grad2_0_0_gt2[pp];
                double DENDRO_613 = -DENDRO_345*grad2_1_2_gt2[pp];
                double DENDRO_614 = -DENDRO_346*grad2_0_1_gt2[pp];
                double DENDRO_615 = DENDRO_41*gt2[pp];
                double DENDRO_616 = -DENDRO_121*(DENDRO_587 + DENDRO_590*(DENDRO_159 + DENDRO_588 + DENDRO_589) + DENDRO_591*DENDRO_592 + DENDRO_595*(DENDRO_178 + DENDRO_593 + DENDRO_594)) + DENDRO_130*DENDRO_412 + DENDRO_148*DENDRO_414 + DENDRO_416*DENDRO_72 + DENDRO_419*DENDRO_615 + DENDRO_597 + DENDRO_599 + DENDRO_601 + DENDRO_602 + DENDRO_604 + DENDRO_606 + DENDRO_608 + DENDRO_609 + DENDRO_610 + DENDRO_611 + DENDRO_612 + DENDRO_613 + DENDRO_614;
                double DENDRO_617 = grad2_0_2_alpha[pp];
                double DENDRO_618 = DENDRO_157*DENDRO_41;
                double DENDRO_619 = DENDRO_158*DENDRO_41;
                double DENDRO_620 = DENDRO_159*DENDRO_41;
                double DENDRO_621 = -DENDRO_61;
                double DENDRO_622 = DENDRO_369*DENDRO_615;
                double DENDRO_623 = DENDRO_64*(DENDRO_621 + DENDRO_622);
                double DENDRO_624 = 2.0*DENDRO_94;
                double DENDRO_625 = DENDRO_176*DENDRO_41;
                double DENDRO_626 = DENDRO_177*DENDRO_41;
                double DENDRO_627 = DENDRO_178*DENDRO_41;
                double DENDRO_628 = -DENDRO_60;
                double DENDRO_629 = DENDRO_108*DENDRO_615;
                double DENDRO_630 = DENDRO_64*(DENDRO_628 + DENDRO_629);
                double DENDRO_631 = 2.0*DENDRO_112;
                double DENDRO_632 = 2.0*DENDRO_82;
                double DENDRO_633 = DENDRO_64*gt2[pp];
                double DENDRO_634 = DENDRO_41*(DENDRO_134 + DENDRO_63*DENDRO_633);
                double DENDRO_635 = -4*DENDRO_617 + DENDRO_624*(-DENDRO_618 - DENDRO_619 + DENDRO_620 + DENDRO_623) + DENDRO_631*(-DENDRO_625 - DENDRO_626 + DENDRO_627 + DENDRO_630) + DENDRO_632*DENDRO_634;
                double DENDRO_636 = DENDRO_240*DENDRO_74;
                double DENDRO_637 = DENDRO_132*DENDRO_242 + DENDRO_636;
                double DENDRO_638 = 0.5*DENDRO_392;
                double DENDRO_639 = 0.25*DENDRO_583;
                double DENDRO_640 = DENDRO_160*DENDRO_165 + DENDRO_547;
                double DENDRO_641 = DENDRO_141*DENDRO_385;
                double DENDRO_642 = -DENDRO_641;
                double DENDRO_643 = -DENDRO_242*DENDRO_389;
                double DENDRO_644 = DENDRO_405 + DENDRO_569;
                double DENDRO_645 = DENDRO_515 + DENDRO_574;
                double DENDRO_646 = DENDRO_174*DENDRO_259;
                double DENDRO_647 = DENDRO_119*DENDRO_332;
                double DENDRO_648 = DENDRO_240*DENDRO_69;
                double DENDRO_649 = DENDRO_537 + DENDRO_648;
                double DENDRO_650 = DENDRO_134*DENDRO_202;
                double DENDRO_651 = DENDRO_141*DENDRO_204 + DENDRO_650;
                double DENDRO_652 = 0.25*DENDRO_530;
                double DENDRO_653 = DENDRO_226*DENDRO_80;
                double DENDRO_654 = DENDRO_180*DENDRO_208 + DENDRO_653;
                double DENDRO_655 = DENDRO_145*(DENDRO_557 + DENDRO_651) - DENDRO_184*(DENDRO_546 + DENDRO_549) - DENDRO_184*(-DENDRO_218*DENDRO_80 + DENDRO_550) + DENDRO_190*(DENDRO_143 + DENDRO_573) + DENDRO_190*(DENDRO_652 + DENDRO_654) - DENDRO_302*(DENDRO_387 + DENDRO_536) - DENDRO_306*(DENDRO_502 + DENDRO_541) - DENDRO_315*(0.5*DENDRO_343 + DENDRO_544);
                double DENDRO_656 = DENDRO_160*DENDRO_202;
                double DENDRO_657 = DENDRO_122*DENDRO_235;
                double DENDRO_658 = 0.25*DENDRO_439;
                double DENDRO_659 = DENDRO_218*DENDRO_265;
                double DENDRO_660 = -DENDRO_659;
                double DENDRO_661 = DENDRO_165*DENDRO_240 + DENDRO_247*DENDRO_470;
                double DENDRO_662 = DENDRO_227*DENDRO_230;
                double DENDRO_663 = DENDRO_259*DENDRO_437;
                double DENDRO_664 = DENDRO_250*DENDRO_408;
                double DENDRO_665 = DENDRO_240*DENDRO_304 + DENDRO_664;
                double DENDRO_666 = DENDRO_235*DENDRO_327;
                double DENDRO_667 = DENDRO_580 + DENDRO_666;
                double DENDRO_668 = DENDRO_247*DENDRO_71 + 0.5*DENDRO_401;
                double DENDRO_669 = DENDRO_131*DENDRO_235;
                double DENDRO_670 = 0.25*DENDRO_669;
                double DENDRO_671 = DENDRO_208*DENDRO_561 + DENDRO_559;
                double DENDRO_672 = DENDRO_265*DENDRO_470;
                double DENDRO_673 = DENDRO_160*DENDRO_385;
                double DENDRO_674 = -DENDRO_673;
                double DENDRO_675 = DENDRO_212*DENDRO_268;
                double DENDRO_676 = -DENDRO_259*DENDRO_385 + DENDRO_675;
                double DENDRO_677 = DENDRO_147*DENDRO_250;
                double DENDRO_678 = DENDRO_131*DENDRO_578 + DENDRO_677;
                double DENDRO_679 = DENDRO_240*DENDRO_576;
                double DENDRO_680 = DENDRO_123*DENDRO_247 + DENDRO_240*DENDRO_561;
                double DENDRO_681 = DENDRO_208*DENDRO_476;
                double DENDRO_682 = DENDRO_219*DENDRO_222;
                double DENDRO_683 = DENDRO_208*DENDRO_437 + DENDRO_682;
                double DENDRO_684 = -DENDRO_235*DENDRO_472;
                double DENDRO_685 = DENDRO_250*DENDRO_76;
                double DENDRO_686 = 0.5*DENDRO_463 + DENDRO_685;
                double DENDRO_687 = DENDRO_160*DENDRO_476;
                double DENDRO_688 = 0.25*DENDRO_131;
                double DENDRO_689 = DENDRO_147*DENDRO_268;
                double DENDRO_690 = DENDRO_259*DENDRO_688 + DENDRO_689;
                double DENDRO_691 = DENDRO_221*DENDRO_545;
                double DENDRO_692 = DENDRO_503 + DENDRO_691;
                double DENDRO_693 = DENDRO_230*DENDRO_408;
                double DENDRO_694 = DENDRO_131*DENDRO_160;
                double DENDRO_695 = 1.0*DENDRO_314;
                double DENDRO_696 = -DENDRO_282;
                double DENDRO_697 = DENDRO_204*DENDRO_368;
                double DENDRO_698 = DENDRO_105*DENDRO_167;
                double DENDRO_699 = DENDRO_202*DENDRO_382;
                double DENDRO_700 = DENDRO_105*DENDRO_204 + DENDRO_107*DENDRO_167 + DENDRO_237;
                double DENDRO_701 = DENDRO_41*gt4[pp];
                double DENDRO_702 = DENDRO_330*grad2_0_2_gt4[pp] + DENDRO_335*grad2_2_2_gt4[pp] + DENDRO_337*grad2_1_1_gt4[pp] + DENDRO_339*grad2_0_0_gt4[pp] - DENDRO_345*grad2_1_2_gt4[pp] - DENDRO_346*grad2_0_1_gt4[pp] + DENDRO_432*DENDRO_603 + DENDRO_499*DENDRO_600 + DENDRO_500*DENDRO_603 + DENDRO_501*DENDRO_605 - DENDRO_59*DENDRO_607 + DENDRO_596*gt1[pp] + DENDRO_598*gt3[pp];
                double DENDRO_703 = -DENDRO_121*(DENDRO_590*(DENDRO_256 + DENDRO_697 + DENDRO_698) + DENDRO_592*(DENDRO_207 + DENDRO_699) + DENDRO_595*DENDRO_700 + DENDRO_696) + DENDRO_128*DENDRO_416 + DENDRO_215*DENDRO_414 + DENDRO_224*DENDRO_412 + DENDRO_419*DENDRO_701 + DENDRO_702;
                double DENDRO_704 = grad2_1_2_alpha[pp];
                double DENDRO_705 = DENDRO_256*DENDRO_41;
                double DENDRO_706 = DENDRO_257*DENDRO_41;
                double DENDRO_707 = DENDRO_258*DENDRO_41;
                double DENDRO_708 = -DENDRO_59;
                double DENDRO_709 = DENDRO_369*DENDRO_701;
                double DENDRO_710 = DENDRO_64*(DENDRO_708 + DENDRO_709);
                double DENDRO_711 = DENDRO_450*DENDRO_701;
                double DENDRO_712 = DENDRO_205*DENDRO_41 + DENDRO_206*DENDRO_41;
                double DENDRO_713 = -DENDRO_203*DENDRO_41 + DENDRO_712;
                double DENDRO_714 = DENDRO_64*gt4[pp];
                double DENDRO_715 = DENDRO_359*DENDRO_714;
                double DENDRO_716 = DENDRO_41*DENDRO_631*(DENDRO_240 + DENDRO_715) + DENDRO_624*(DENDRO_705 - DENDRO_706 - DENDRO_707 + DENDRO_710) + DENDRO_632*(DENDRO_64*(DENDRO_628 + DENDRO_711) + DENDRO_713) - 4*DENDRO_704;
                double DENDRO_717 = 0.5*DENDRO_398;
                double DENDRO_718 = 1.0*DENDRO_490;
                double DENDRO_719 = 0.5*DENDRO_489;
                double DENDRO_720 = 0.25*DENDRO_694;
                double DENDRO_721 = DENDRO_406 + DENDRO_569;
                double DENDRO_722 = DENDRO_151*DENDRO_250;
                double DENDRO_723 = DENDRO_681 + DENDRO_682;
                double DENDRO_724 = DENDRO_259*DENDRO_404;
                double DENDRO_725 = DENDRO_250*DENDRO_75;
                double DENDRO_726 = DENDRO_235*DENDRO_304 + DENDRO_725;
                double DENDRO_727 = DENDRO_160*DENDRO_437 + DENDRO_689;
                double DENDRO_728 = DENDRO_502 + DENDRO_503;
                double DENDRO_729 = DENDRO_230*DENDRO_470 + DENDRO_691;
                double DENDRO_730 = 1.0*DENDRO_183;
                double DENDRO_731 = -DENDRO_184*(DENDRO_642 + DENDRO_671) - DENDRO_302*(DENDRO_219*DENDRO_226 + DENDRO_433 + DENDRO_658) - DENDRO_695*(DENDRO_142 + DENDRO_532) - DENDRO_730*(DENDRO_435 + DENDRO_651);
                double DENDRO_732 = DENDRO_567 + DENDRO_674;
                double DENDRO_733 = 0.5*DENDRO_486;
                double DENDRO_734 = DENDRO_127*DENDRO_230;
                double DENDRO_735 = -DENDRO_734;
                double DENDRO_736 = DENDRO_268*DENDRO_332;
                double DENDRO_737 = DENDRO_155*DENDRO_437 + DENDRO_736;
                double DENDRO_738 = DENDRO_242*DENDRO_304;
                double DENDRO_739 = DENDRO_250*DENDRO_70;
                double DENDRO_740 = DENDRO_119*DENDRO_71;
                double DENDRO_741 = DENDRO_116*DENDRO_470 + DENDRO_155*DENDRO_174;
                double DENDRO_742 = DENDRO_116*DENDRO_218;
                double DENDRO_743 = 0.5*DENDRO_156;
                double DENDRO_744 = -DENDRO_742 - DENDRO_743;
                double DENDRO_745 = DENDRO_119*DENDRO_470 + DENDRO_517;
                double DENDRO_746 = DENDRO_579 + DENDRO_666;
                double DENDRO_747 = DENDRO_151*DENDRO_268;
                double DENDRO_748 = 0.5*DENDRO_466;
                double DENDRO_749 = -DENDRO_747 - DENDRO_748;
                double DENDRO_750 = DENDRO_242*DENDRO_688 + DENDRO_725;
                double DENDRO_751 = DENDRO_230*DENDRO_332;
                double DENDRO_752 = DENDRO_502 + DENDRO_691;
                double DENDRO_753 = DENDRO_268*DENDRO_75;
                double DENDRO_754 = DENDRO_131*DENDRO_568 + DENDRO_753;
                double DENDRO_755 = DENDRO_250*DENDRO_78;
                double DENDRO_756 = DENDRO_242*DENDRO_327 + DENDRO_755;
                double DENDRO_757 = DENDRO_155*DENDRO_561;
                double DENDRO_758 = DENDRO_116*DENDRO_226 + DENDRO_167*DENDRO_568;
                double DENDRO_759 = DENDRO_222*DENDRO_80;
                double DENDRO_760 = DENDRO_141*DENDRO_304 + DENDRO_759;
                double DENDRO_761 = 1.0*DENDRO_301;
                double DENDRO_762 = -DENDRO_289;
                double DENDRO_763 = DENDRO_105*DENDRO_69 + DENDRO_137*DENDRO_368 + DENDRO_152;
                double DENDRO_764 = DENDRO_122*DENDRO_382;
                double DENDRO_765 = DENDRO_107*DENDRO_69;
                double DENDRO_766 = DENDRO_105*DENDRO_137;
                double DENDRO_767 = 2.0*DENDRO_499*gt0[pp];
                double DENDRO_768 = 2.0*gt1[pp];
                double DENDRO_769 = DENDRO_322*DENDRO_768;
                double DENDRO_770 = DENDRO_500*DENDRO_768;
                double DENDRO_771 = DENDRO_501*DENDRO_600;
                double DENDRO_772 = 2.0*DENDRO_323*gt3[pp];
                double DENDRO_773 = DENDRO_325*DENDRO_603;
                double DENDRO_774 = -DENDRO_290*DENDRO_321;
                double DENDRO_775 = DENDRO_330*grad2_0_2_gt1[pp];
                double DENDRO_776 = DENDRO_335*grad2_2_2_gt1[pp];
                double DENDRO_777 = DENDRO_337*grad2_1_1_gt1[pp];
                double DENDRO_778 = DENDRO_339*grad2_0_0_gt1[pp];
                double DENDRO_779 = -DENDRO_345*grad2_1_2_gt1[pp];
                double DENDRO_780 = -DENDRO_346*grad2_0_1_gt1[pp];
                double DENDRO_781 = DENDRO_41*gt1[pp];
                double DENDRO_782 = -DENDRO_121*(DENDRO_590*DENDRO_763 + DENDRO_592*(DENDRO_140 + DENDRO_764) + DENDRO_595*(DENDRO_170 + DENDRO_765 + DENDRO_766) + DENDRO_762) + DENDRO_124*DENDRO_412 + DENDRO_129*DENDRO_414 + DENDRO_416*DENDRO_67 + DENDRO_419*DENDRO_781 + DENDRO_767 + DENDRO_769 + DENDRO_770 + DENDRO_771 + DENDRO_772 + DENDRO_773 + DENDRO_774 + DENDRO_775 + DENDRO_776 + DENDRO_777 + DENDRO_778 + DENDRO_779 + DENDRO_780;
                double DENDRO_783 = 0.25*DENDRO_340;
                double DENDRO_784 = DENDRO_141*DENDRO_576 + DENDRO_653;
                double DENDRO_785 = -DENDRO_141*DENDRO_472;
                double DENDRO_786 = DENDRO_146*(DENDRO_540 + DENDRO_752) - DENDRO_184*(DENDRO_331 + DENDRO_654) - DENDRO_184*(DENDRO_331 + DENDRO_784) + DENDRO_190*(DENDRO_760 + DENDRO_785) - DENDRO_302*(DENDRO_134*DENDRO_226 + DENDRO_558) - DENDRO_315*(1.0*DENDRO_341 + DENDRO_783);
                double DENDRO_787 = grad2_0_1_alpha[pp];
                double DENDRO_788 = DENDRO_450*DENDRO_781;
                double DENDRO_789 = DENDRO_138*DENDRO_41 + DENDRO_139*DENDRO_41;
                double DENDRO_790 = -DENDRO_136*DENDRO_41 + DENDRO_789;
                double DENDRO_791 = DENDRO_170*DENDRO_41;
                double DENDRO_792 = DENDRO_171*DENDRO_41;
                double DENDRO_793 = DENDRO_172*DENDRO_41;
                double DENDRO_794 = DENDRO_108*DENDRO_781;
                double DENDRO_795 = DENDRO_64*(DENDRO_708 + DENDRO_794);
                double DENDRO_796 = DENDRO_64*gt1[pp];
                double DENDRO_797 = DENDRO_796*DENDRO_88;
                double DENDRO_798 = DENDRO_41*DENDRO_624*(DENDRO_155 + DENDRO_797) + DENDRO_631*(DENDRO_791 - DENDRO_792 - DENDRO_793 + DENDRO_795) + DENDRO_632*(DENDRO_64*(DENDRO_621 + DENDRO_788) + DENDRO_790) - 4*DENDRO_787;
                double DENDRO_799 = DENDRO_180*DENDRO_242;
                double DENDRO_800 = -DENDRO_29*(DENDRO_798 + alpha[pp]*(DENDRO_145*(DENDRO_464 + DENDRO_648 + DENDRO_657) + DENDRO_145*(DENDRO_583 + DENDRO_584 + DENDRO_656) + DENDRO_146*(DENDRO_684 + DENDRO_750) + DENDRO_146*(DENDRO_687 + DENDRO_749) + DENDRO_146*(DENDRO_751 + DENDRO_752) - DENDRO_184*(DENDRO_162 + DENDRO_744) - DENDRO_184*(DENDRO_516 + DENDRO_745) - DENDRO_184*(DENDRO_647 + DENDRO_746) - DENDRO_184*(DENDRO_574 + DENDRO_646 + DENDRO_720) + DENDRO_190*(DENDRO_757 + DENDRO_758) + DENDRO_190*(DENDRO_116*DENDRO_227 + DENDRO_754) + DENDRO_190*(-DENDRO_119*DENDRO_127 + DENDRO_756) + DENDRO_190*(DENDRO_230*DENDRO_70 + DENDRO_760) + DENDRO_200*(DENDRO_119*DENDRO_122 + DENDRO_525) - DENDRO_302*(DENDRO_405 + DENDRO_732) - DENDRO_306*(DENDRO_482 + DENDRO_737) - DENDRO_306*(DENDRO_733 + DENDRO_735) - DENDRO_306*(DENDRO_484 + DENDRO_738 + DENDRO_739) - DENDRO_315*(0.5*DENDRO_522 + DENDRO_741) - DENDRO_315*(DENDRO_119*DENDRO_70 + DENDRO_520 + DENDRO_740) - DENDRO_761*(DENDRO_402 + DENDRO_636 + DENDRO_669) + DENDRO_782 + DENDRO_786)) - DENDRO_29*(DENDRO_798 + alpha[pp]*(DENDRO_146*(DENDRO_468 + DENDRO_749) + DENDRO_146*(DENDRO_471 + DENDRO_726) + DENDRO_146*(DENDRO_471 + DENDRO_750) + DENDRO_146*(DENDRO_540 + DENDRO_729) + DENDRO_146*(DENDRO_639 + DENDRO_727) + DENDRO_146*(DENDRO_728 + DENDRO_751) - DENDRO_184*(DENDRO_517 + DENDRO_746) - DENDRO_184*(DENDRO_572 + DENDRO_784) - DENDRO_184*(DENDRO_575 + DENDRO_744) - DENDRO_184*(DENDRO_579 + DENDRO_745) + DENDRO_190*(DENDRO_754 + DENDRO_757) + DENDRO_190*(DENDRO_756 + DENDRO_799) + DENDRO_190*(DENDRO_268*DENDRO_76 + DENDRO_758) + DENDRO_190*(DENDRO_119*DENDRO_123 + DENDRO_755 + DENDRO_799) + DENDRO_190*(DENDRO_230*DENDRO_71 + DENDRO_759 + DENDRO_785) + DENDRO_200*(DENDRO_230*DENDRO_69 + DENDRO_505) - DENDRO_302*(DENDRO_406 + DENDRO_732) - DENDRO_302*(DENDRO_240*DENDRO_75 + DENDRO_670) - DENDRO_306*(0.5*DENDRO_491 + DENDRO_737) - DENDRO_306*(1.0*DENDRO_492 + DENDRO_738) - DENDRO_306*(DENDRO_123*DENDRO_230 + DENDRO_487 + DENDRO_735) - DENDRO_315*(DENDRO_518 + DENDRO_741) - DENDRO_315*(0.5*DENDRO_519 + DENDRO_740) - DENDRO_315*(DENDRO_123*DENDRO_80 + DENDRO_333 + DENDRO_783) - DENDRO_730*(DENDRO_142 + DENDRO_529 + DENDRO_531) - DENDRO_730*(DENDRO_526 + DENDRO_527 + DENDRO_694) - DENDRO_761*(DENDRO_435 + DENDRO_557 + DENDRO_650) + DENDRO_782)) - DENDRO_31*(DENDRO_716 + alpha[pp]*(DENDRO_146*(DENDRO_676 + DENDRO_724) + DENDRO_146*(DENDRO_678 + DENDRO_679) + DENDRO_146*(DENDRO_680 - DENDRO_722) + DENDRO_146*(-DENDRO_218*DENDRO_230 + DENDRO_723) + DENDRO_146*(DENDRO_226*DENDRO_265 + DENDRO_675 + DENDRO_724) - DENDRO_184*(DENDRO_643 + DENDRO_668) - DENDRO_184*(DENDRO_672 + DENDRO_721) - DENDRO_184*(DENDRO_674 + DENDRO_721) + DENDRO_190*(DENDRO_465 + DENDRO_686) + DENDRO_190*(DENDRO_469 + DENDRO_690) + DENDRO_190*(DENDRO_469 + DENDRO_727) + DENDRO_190*(DENDRO_538 + DENDRO_726) + DENDRO_190*(DENDRO_539 + DENDRO_729) + DENDRO_190*(DENDRO_693 + DENDRO_728) - DENDRO_302*(DENDRO_396 + DENDRO_661) - DENDRO_302*(DENDRO_660 + DENDRO_717) - DENDRO_306*(DENDRO_663 + DENDRO_718) - DENDRO_306*(DENDRO_665 + DENDRO_719) - DENDRO_306*(DENDRO_226*DENDRO_230 + DENDRO_480 + DENDRO_662) - DENDRO_315*(DENDRO_517 + DENDRO_667) - DENDRO_315*(DENDRO_147*DENDRO_155 + DENDRO_720) + DENDRO_410*(DENDRO_204*DENDRO_230 + DENDRO_504) + DENDRO_703 - DENDRO_730*(DENDRO_637 + DENDRO_669) + DENDRO_731)) - DENDRO_31*(DENDRO_716 + alpha[pp]*(DENDRO_146*(DENDRO_679 + DENDRO_680) + DENDRO_146*(DENDRO_681 + DENDRO_683) + DENDRO_146*(-DENDRO_127*DENDRO_247 + DENDRO_678) + DENDRO_146*(DENDRO_214*DENDRO_230 + DENDRO_683) + DENDRO_146*(DENDRO_227*DENDRO_265 + DENDRO_676) - DENDRO_184*(DENDRO_403 + DENDRO_668) - DENDRO_184*(DENDRO_438 + DENDRO_560) - DENDRO_184*(DENDRO_438 + DENDRO_671) - DENDRO_184*(DENDRO_565 + DENDRO_672) - DENDRO_184*(DENDRO_570 + DENDRO_674) - DENDRO_184*(DENDRO_571 + DENDRO_670) + DENDRO_189*(DENDRO_585 + DENDRO_656) + DENDRO_189*(DENDRO_649 + DENDRO_657) + DENDRO_190*(DENDRO_539 + DENDRO_692) + DENDRO_190*(DENDRO_684 + DENDRO_686) + DENDRO_190*(DENDRO_687 + DENDRO_690) + DENDRO_190*(DENDRO_692 + DENDRO_693) - DENDRO_302*(0.5*DENDRO_409 + DENDRO_661) - DENDRO_302*(1.0*DENDRO_440 + DENDRO_658) - DENDRO_302*(DENDRO_214*DENDRO_265 + DENDRO_399 + DENDRO_660) - DENDRO_306*(DENDRO_474 + DENDRO_665) - DENDRO_306*(0.5*DENDRO_479 + DENDRO_662) - DENDRO_306*(DENDRO_214*DENDRO_268 + DENDRO_477 + DENDRO_663) - DENDRO_315*(DENDRO_516 + DENDRO_667) - DENDRO_315*(DENDRO_123*DENDRO_134 + DENDRO_652) + DENDRO_410*(DENDRO_202*DENDRO_265 + DENDRO_425) - DENDRO_695*(DENDRO_528 + DENDRO_694) + DENDRO_703)) + DENDRO_35*(DENDRO_355*(DENDRO_455 + DENDRO_64*(DENDRO_449 + DENDRO_453)) + DENDRO_362*(DENDRO_250 + DENDRO_359*DENDRO_446) - 4*DENDRO_445 + DENDRO_447*(DENDRO_268 + DENDRO_446*DENDRO_88) + alpha[pp]*(-DENDRO_121*(DENDRO_117*DENDRO_459 + DENDRO_118*(DENDRO_228 + DENDRO_460 + DENDRO_461) + DENDRO_120*DENDRO_462 + DENDRO_456) + DENDRO_122*DENDRO_416 + DENDRO_146*(DENDRO_474 + DENDRO_475) + DENDRO_146*(DENDRO_477 - DENDRO_478) + DENDRO_146*(DENDRO_480 + 1.0*DENDRO_481) - DENDRO_155*DENDRO_498 - DENDRO_184*(DENDRO_463 + DENDRO_465) - DENDRO_184*(-1.0*DENDRO_466 + DENDRO_468) - DENDRO_184*(DENDRO_242*DENDRO_470 + DENDRO_471) - DENDRO_184*(DENDRO_259*DENDRO_470 + DENDRO_469) + DENDRO_190*(DENDRO_482 + DENDRO_483) + DENDRO_190*(DENDRO_484 + DENDRO_485) + DENDRO_190*(DENDRO_487 + 1.0*DENDRO_488) + DENDRO_200*(DENDRO_486 + DENDRO_488) + DENDRO_200*(DENDRO_122*DENDRO_242 + DENDRO_492) + DENDRO_200*(DENDRO_155*DENDRO_202 + DENDRO_491) + DENDRO_202*DENDRO_414 - DENDRO_221*DENDRO_231*DENDRO_428 + DENDRO_221*DENDRO_412 - DENDRO_240*DENDRO_495 - DENDRO_242*DENDRO_497 + DENDRO_410*(DENDRO_479 + DENDRO_481) + DENDRO_410*(DENDRO_122*DENDRO_240 + DENDRO_489) + DENDRO_410*(DENDRO_202*DENDRO_259 + DENDRO_490) + DENDRO_419*DENDRO_451 - DENDRO_493*(DENDRO_216 - DENDRO_404) - DENDRO_496*(DENDRO_123 + DENDRO_125) + DENDRO_506)) + DENDRO_46*(DENDRO_635 + alpha[pp]*(DENDRO_145*(DENDRO_402 + DENDRO_637) + DENDRO_146*(DENDRO_556 + DENDRO_642) + DENDRO_146*(DENDRO_563 + DENDRO_643) + DENDRO_146*(DENDRO_564 + DENDRO_644) + DENDRO_146*(DENDRO_567 + DENDRO_644) - DENDRO_184*(DENDRO_548 + DENDRO_640) - DENDRO_184*(-DENDRO_119*DENDRO_151 + DENDRO_553) - DENDRO_184*(DENDRO_265*DENDRO_75 + DENDRO_640) + DENDRO_190*(DENDRO_517 + DENDRO_582) + DENDRO_190*(DENDRO_577 + DENDRO_645) + DENDRO_190*(DENDRO_581 + DENDRO_647) + DENDRO_190*(DENDRO_645 + DENDRO_646) - DENDRO_197*(DENDRO_119*DENDRO_132 + DENDRO_524) - DENDRO_302*(DENDRO_535 + DENDRO_638) - DENDRO_302*(DENDRO_247*DENDRO_75 + DENDRO_390 + DENDRO_533) - DENDRO_306*(DENDRO_155*DENDRO_214 + DENDRO_639) - DENDRO_315*(1.0*DENDRO_523 + DENDRO_543) - DENDRO_315*(DENDRO_119*DENDRO_75 + DENDRO_513 + DENDRO_542) - DENDRO_586*(DENDRO_464 + DENDRO_649) + DENDRO_616 + DENDRO_655)) + DENDRO_46*(DENDRO_635 + alpha[pp]*(DENDRO_146*(DENDRO_407 + DENDRO_563) + DENDRO_146*(DENDRO_407 + DENDRO_571) + DENDRO_146*(DENDRO_436 + DENDRO_556) + DENDRO_146*(DENDRO_558 + DENDRO_560) + DENDRO_146*(DENDRO_564 + DENDRO_565) + DENDRO_146*(DENDRO_567 + DENDRO_570) - DENDRO_184*(DENDRO_549 + DENDRO_550) - DENDRO_184*(DENDRO_551 + DENDRO_553) - DENDRO_184*(DENDRO_219*DENDRO_71 + DENDRO_546) - DENDRO_184*(DENDRO_119*DENDRO_147 + DENDRO_551 + DENDRO_552) - DENDRO_184*(DENDRO_265*DENDRO_76 + DENDRO_547 + DENDRO_548) + DENDRO_189*(DENDRO_526 + DENDRO_528) + DENDRO_189*(DENDRO_529 + DENDRO_532) + DENDRO_190*(DENDRO_516 + DENDRO_581) + DENDRO_190*(DENDRO_572 + DENDRO_573) + DENDRO_190*(DENDRO_579 + DENDRO_582) + DENDRO_190*(DENDRO_574 + DENDRO_575 + DENDRO_577) - DENDRO_197*(DENDRO_265*DENDRO_74 + DENDRO_426) - DENDRO_302*(1.0*DENDRO_417 + DENDRO_533) - DENDRO_302*(0.5*DENDRO_442 + DENDRO_536) - DENDRO_302*(DENDRO_147*DENDRO_265 + DENDRO_393 + DENDRO_535) - DENDRO_306*(DENDRO_503 + DENDRO_541) - DENDRO_306*(DENDRO_240*DENDRO_70 + DENDRO_538) - DENDRO_315*(DENDRO_328 + DENDRO_544) - DENDRO_315*(0.5*DENDRO_512 + DENDRO_542) - DENDRO_315*(DENDRO_116*DENDRO_147 + DENDRO_511 + DENDRO_543) - DENDRO_586*(DENDRO_583 + DENDRO_585) + DENDRO_616)) + DENDRO_49*(DENDRO_356*DENDRO_81 + DENDRO_361*(-DENDRO_100 + DENDRO_102 + DENDRO_111 - DENDRO_98) + DENDRO_447*(DENDRO_116 + DENDRO_89) - 4*DENDRO_58 + alpha[pp]*(-DENDRO_121*(DENDRO_115 + DENDRO_117*DENDRO_507 + DENDRO_118*DENDRO_508 + DENDRO_120*(DENDRO_101 + DENDRO_509 + DENDRO_510)) - DENDRO_144*DENDRO_252*DENDRO_320 + DENDRO_146*(-1.0*DENDRO_135 + DENDRO_143) + DENDRO_146*(-1.0*DENDRO_156 + DENDRO_162) + DENDRO_146*(DENDRO_160*DENDRO_332 + DENDRO_515) + DENDRO_146*(DENDRO_235*DENDRO_70 + DENDRO_517) + DENDRO_146*(DENDRO_242*DENDRO_75 + DENDRO_516) - DENDRO_155*DENDRO_312 - DENDRO_184*(DENDRO_513 + 1.0*DENDRO_514) - DENDRO_184*(-DENDRO_116*DENDRO_182 + DENDRO_511) + DENDRO_190*(DENDRO_520 + 1.0*DENDRO_521) + DENDRO_190*(DENDRO_116*DENDRO_167 + DENDRO_518) - DENDRO_197*(DENDRO_512 + DENDRO_514) - DENDRO_197*(DENDRO_160*DENDRO_74 + DENDRO_523) + DENDRO_200*(DENDRO_519 + DENDRO_521) + DENDRO_200*(DENDRO_155*DENDRO_74 + DENDRO_522) + DENDRO_234*DENDRO_411 + DENDRO_254*DENDRO_415 - DENDRO_270*DENDRO_313*DENDRO_349 + DENDRO_272*DENDRO_413 + DENDRO_300*DENDRO_418 - DENDRO_303*(DENDRO_149 - DENDRO_165) - DENDRO_307*(DENDRO_125 - DENDRO_304) - DENDRO_317*DENDRO_524 - DENDRO_319*DENDRO_525 + DENDRO_350)) + DENDRO_55*(-4*DENDRO_352 + DENDRO_354*DENDRO_356 + DENDRO_362*(DENDRO_247 + DENDRO_360) + DENDRO_373*(-DENDRO_363 + DENDRO_364 - DENDRO_365 + DENDRO_372) + alpha[pp]*(-DENDRO_121*(DENDRO_117*(DENDRO_375 + DENDRO_378 + DENDRO_381) + DENDRO_118*DENDRO_383 + DENDRO_120*DENDRO_384 + DENDRO_374) + DENDRO_132*DENDRO_416 + DENDRO_146*(DENDRO_396 + DENDRO_397) + DENDRO_146*(DENDRO_399 + 1.0*DENDRO_400) - DENDRO_184*(DENDRO_387 + DENDRO_388) - DENDRO_184*(DENDRO_393 + 1.0*DENDRO_394) - DENDRO_184*(DENDRO_247*DENDRO_391 + DENDRO_390) + DENDRO_190*(DENDRO_401 + DENDRO_403) + DENDRO_190*(DENDRO_147*DENDRO_259 + DENDRO_405) + DENDRO_190*(DENDRO_160*DENDRO_214 + DENDRO_406) + DENDRO_190*(DENDRO_235*DENDRO_408 + DENDRO_407) - DENDRO_197*(DENDRO_392 + DENDRO_394) - DENDRO_197*(DENDRO_132*DENDRO_235 + DENDRO_417) + DENDRO_204*DENDRO_412 - DENDRO_211*DENDRO_266*DENDRO_428 + DENDRO_211*DENDRO_414 - DENDRO_235*DENDRO_424 - DENDRO_240*DENDRO_423 - DENDRO_319*DENDRO_425 + DENDRO_367*DENDRO_419 + DENDRO_410*(DENDRO_398 + DENDRO_400) + DENDRO_410*(DENDRO_132*DENDRO_240 + DENDRO_409) - DENDRO_420*(DENDRO_147 + DENDRO_149) - DENDRO_421*(DENDRO_214 + DENDRO_216) - DENDRO_426*DENDRO_427 + DENDRO_444));
                double DENDRO_801 = DENDRO_41*DENDRO_800;
                double DENDRO_802 = (1.0/12.0)*chi[pp];
                double DENDRO_803 = (1.0/3.0)*At1[pp];
                double DENDRO_804 = At4[pp]*DENDRO_31;
                double DENDRO_805 = At3[pp]*DENDRO_35;
                double DENDRO_806 = DENDRO_43 + DENDRO_804 - DENDRO_805;
                double DENDRO_807 = At3[pp]*DENDRO_29;
                double DENDRO_808 = At4[pp]*DENDRO_46;
                double DENDRO_809 = -At1[pp]*DENDRO_49 + DENDRO_807 - DENDRO_808;
                double DENDRO_810 = At4[pp]*DENDRO_55;
                double DENDRO_811 = -At1[pp]*DENDRO_46 + At3[pp]*DENDRO_31 - DENDRO_810;
                double DENDRO_812 = 6.0*DENDRO_94;
                double DENDRO_813 = 6.0*DENDRO_112;
                double DENDRO_814 = 6.0*DENDRO_82;
                double DENDRO_815 = DENDRO_131*DENDRO_179;
                double DENDRO_816 = DENDRO_167*DENDRO_179;
                double DENDRO_817 = -DENDRO_237 + DENDRO_238 + DENDRO_239;
                double DENDRO_818 = DENDRO_74*DENDRO_817;
                double DENDRO_819 = DENDRO_816 + DENDRO_818;
                double DENDRO_820 = DENDRO_173*DENDRO_472;
                double DENDRO_821 = -DENDRO_164*DENDRO_476;
                double DENDRO_822 = DENDRO_223 - DENDRO_228 + DENDRO_229;
                double DENDRO_823 = DENDRO_168*DENDRO_688;
                double DENDRO_824 = -DENDRO_256 + DENDRO_257 + DENDRO_258;
                double DENDRO_825 = DENDRO_214*DENDRO_93;
                double DENDRO_826 = DENDRO_174*DENDRO_824 + DENDRO_825;
                double DENDRO_827 = DENDRO_179*DENDRO_327;
                double DENDRO_828 = DENDRO_187*DENDRO_332 + 0.25*DENDRO_77*DENDRO_817;
                double DENDRO_829 = DENDRO_175 + DENDRO_181;
                double DENDRO_830 = DENDRO_167*DENDRO_824;
                double DENDRO_831 = 1.0*DENDRO_145;
                double DENDRO_832 = DENDRO_167*DENDRO_173;
                double DENDRO_833 = DENDRO_69*DENDRO_817 + DENDRO_832;
                double DENDRO_834 = 2.0*DENDRO_233;
                double DENDRO_835 = 2.0*DENDRO_253;
                double DENDRO_836 = 2.0*DENDRO_271;
                double DENDRO_837 = DENDRO_299*DENDRO_64;
                double DENDRO_838 = (1.0/3.0)*At2[pp];
                double DENDRO_839 = At5[pp]*DENDRO_31;
                double DENDRO_840 = At2[pp]*DENDRO_29 - At4[pp]*DENDRO_35 + DENDRO_839;
                double DENDRO_841 = At5[pp]*DENDRO_46;
                double DENDRO_842 = -At2[pp]*DENDRO_49 + At4[pp]*DENDRO_29 - DENDRO_841;
                double DENDRO_843 = -At5[pp]*DENDRO_55 + DENDRO_48 + DENDRO_804;
                double DENDRO_844 = DENDRO_137*DENDRO_173;
                double DENDRO_845 = DENDRO_179*DENDRO_389;
                double DENDRO_846 = DENDRO_244 + DENDRO_245 - DENDRO_246;
                double DENDRO_847 = 0.25*DENDRO_164*DENDRO_211;
                double DENDRO_848 = DENDRO_168*DENDRO_404;
                double DENDRO_849 = DENDRO_262 - DENDRO_263 + DENDRO_264;
                double DENDRO_850 = DENDRO_165*DENDRO_168;
                double DENDRO_851 = DENDRO_212*DENDRO_93;
                double DENDRO_852 = DENDRO_173*DENDRO_191;
                double DENDRO_853 = -DENDRO_847;
                double DENDRO_854 = DENDRO_173*DENDRO_389;
                double DENDRO_855 = DENDRO_132*DENDRO_173;
                double DENDRO_856 = 2*At4[pp];
                double DENDRO_857 = At3[pp]*DENDRO_51;
                double DENDRO_858 = DENDRO_41*DENDRO_856;
                double DENDRO_859 = DENDRO_113*DENDRO_41;
                double DENDRO_860 = DENDRO_202*DENDRO_822;
                double DENDRO_861 = 1.0*DENDRO_817;
                double DENDRO_862 = 0.25*DENDRO_832;
                double DENDRO_863 = DENDRO_164*DENDRO_404;
                double DENDRO_864 = DENDRO_122*DENDRO_822;
                double DENDRO_865 = (1.0/3.0)*At4[pp];
                double DENDRO_866 = DENDRO_165*DENDRO_824;
                double DENDRO_867 = -DENDRO_404*DENDRO_824 + DENDRO_675;
                double DENDRO_868 = DENDRO_689 - DENDRO_863;
                double DENDRO_869 = DENDRO_211*DENDRO_824;
                double DENDRO_870 = DENDRO_204*DENDRO_849;
                double DENDRO_871 = DENDRO_168*DENDRO_211;
                double DENDRO_872 = DENDRO_132*DENDRO_849;
                double DENDRO_873 = DENDRO_41*DENDRO_82;
                double DENDRO_874 = DENDRO_109*DENDRO_64;
                double DENDRO_875 = DENDRO_112*DENDRO_41;
                double DENDRO_876 = DENDRO_370*DENDRO_64;
                double DENDRO_877 = DENDRO_41*DENDRO_94;
                double DENDRO_878 = 0.5*DENDRO_94;
                double DENDRO_879 = 0.5*DENDRO_82;
                double DENDRO_880 = DENDRO_284*DENDRO_41;
                double DENDRO_881 = 0.5*DENDRO_112;
                double DENDRO_882 = DENDRO_287*DENDRO_41;
                double DENDRO_883 = DENDRO_291*DENDRO_41;
                double DENDRO_884 = pow(DENDRO_29, 2);
                double DENDRO_885 = pow(DENDRO_46, 2);
                double DENDRO_886 = 2*DENDRO_49;
                double DENDRO_887 = At0[pp]*pow(DENDRO_49, 2) + At3[pp]*DENDRO_884 + At5[pp]*DENDRO_885 - DENDRO_291*DENDRO_808 - DENDRO_43*DENDRO_886 + DENDRO_47*DENDRO_886;
                double DENDRO_888 = 3*DENDRO_144;
                double DENDRO_889 = pow(DENDRO_31, 2);
                double DENDRO_890 = 2*DENDRO_35;
                double DENDRO_891 = At0[pp]*DENDRO_884 + At3[pp]*pow(DENDRO_35, 2) + At5[pp]*DENDRO_889 + DENDRO_291*DENDRO_32 - DENDRO_43*DENDRO_890 - DENDRO_804*DENDRO_890;
                double DENDRO_892 = 2*DENDRO_55;
                double DENDRO_893 = At0[pp]*DENDRO_885 + At3[pp]*DENDRO_889 + At5[pp]*pow(DENDRO_55, 2) - DENDRO_287*DENDRO_53 + DENDRO_47*DENDRO_892 - DENDRO_804*DENDRO_892;
                double DENDRO_894 = At2[pp]*DENDRO_885 - DENDRO_29*DENDRO_810 + DENDRO_31*DENDRO_807 - DENDRO_31*DENDRO_808 - DENDRO_43*DENDRO_46 + DENDRO_46*DENDRO_50 - DENDRO_49*DENDRO_53 + DENDRO_49*DENDRO_56 + DENDRO_55*DENDRO_841;
                double DENDRO_895 = 6*DENDRO_144;
                double DENDRO_896 = At1[pp]*DENDRO_884;
                double DENDRO_897 = DENDRO_29*DENDRO_804;
                double DENDRO_898 = DENDRO_29*DENDRO_47;
                double DENDRO_899 = DENDRO_35*DENDRO_808;
                double DENDRO_900 = DENDRO_31*DENDRO_841;
                double DENDRO_901 = DENDRO_29*DENDRO_50;
                double DENDRO_902 = DENDRO_36*DENDRO_49;
                double DENDRO_903 = DENDRO_32*DENDRO_49;
                double DENDRO_904 = DENDRO_35*DENDRO_807;
                double DENDRO_905 = DENDRO_896 + DENDRO_897 - DENDRO_898 + DENDRO_899 - DENDRO_900 - DENDRO_901 + DENDRO_902 - DENDRO_903 - DENDRO_904;
                double DENDRO_906 = At4[pp]*DENDRO_889;
                double DENDRO_907 = DENDRO_31*DENDRO_43;
                double DENDRO_908 = DENDRO_30*DENDRO_46;
                double DENDRO_909 = DENDRO_36*DENDRO_46;
                double DENDRO_910 = DENDRO_31*DENDRO_47;
                double DENDRO_911 = DENDRO_29*DENDRO_56;
                double DENDRO_912 = DENDRO_31*DENDRO_805;
                double DENDRO_913 = DENDRO_35*DENDRO_810;
                double DENDRO_914 = DENDRO_55*DENDRO_839;
                double DENDRO_915 = DENDRO_906 + DENDRO_907 - DENDRO_908 + DENDRO_909 - DENDRO_910 - DENDRO_911 - DENDRO_912 + DENDRO_913 - DENDRO_914;
                double DENDRO_916 = (1.0/3.0)*alpha[pp];
                double DENDRO_917 = 2*DENDRO_144;
                double DENDRO_918 = DENDRO_112*DENDRO_917;
                double DENDRO_919 = DENDRO_82*DENDRO_917;
                double DENDRO_920 = DENDRO_917*DENDRO_94;
                double DENDRO_921 = grad2_0_2_beta0[pp];
                double DENDRO_922 = (7.0/3.0)*DENDRO_329;
                double DENDRO_923 = grad2_0_0_beta0[pp];
                double DENDRO_924 = 4*grad_0_K[pp];
                double DENDRO_925 = 9*DENDRO_64;
                double DENDRO_926 = DENDRO_120*DENDRO_925;
                double DENDRO_927 = DENDRO_41*DENDRO_916;
                double DENDRO_928 = 4*grad_2_K[pp];
                double DENDRO_929 = DENDRO_117*DENDRO_925;
                double DENDRO_930 = 4*grad_1_K[pp];
                double DENDRO_931 = DENDRO_118*DENDRO_925;
                double DENDRO_932 = grad2_0_1_beta1[pp];
                double DENDRO_933 = (1.0/3.0)*DENDRO_338;
                double DENDRO_934 = grad2_0_2_beta2[pp];
                double DENDRO_935 = grad2_1_2_beta1[pp];
                double DENDRO_936 = (1.0/3.0)*DENDRO_329;
                double DENDRO_937 = grad2_2_2_beta2[pp];
                double DENDRO_938 = grad2_1_1_beta1[pp];
                double DENDRO_939 = DENDRO_29*DENDRO_41;
                double DENDRO_940 = (1.0/3.0)*DENDRO_939;
                double DENDRO_941 = grad2_1_2_beta2[pp];
                double DENDRO_942 = (2.0/3.0)*DENDRO_25;
                double DENDRO_943 = grad2_0_1_beta0[pp];
                double DENDRO_944 = (7.0/3.0)*DENDRO_939;
                double DENDRO_945 = pow(DENDRO_40, -3);
                double DENDRO_946 = DENDRO_0*DENDRO_945;
                double DENDRO_947 = DENDRO_891*DENDRO_946;
                double DENDRO_948 = DENDRO_893*DENDRO_946;
                double DENDRO_949 = DENDRO_887*DENDRO_946;
                double DENDRO_950 = 2.0*DENDRO_945*alpha[pp];
                double DENDRO_951 = DENDRO_894*DENDRO_950;
                double DENDRO_952 = DENDRO_905*DENDRO_950;
                double DENDRO_953 = DENDRO_915*DENDRO_950;
                double DENDRO_954 = beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp];
                double DENDRO_955 = -DENDRO_11*DENDRO_411 + DENDRO_119*DENDRO_949 - DENDRO_15*DENDRO_413 - DENDRO_2*DENDRO_415 + DENDRO_235*DENDRO_951 + DENDRO_240*DENDRO_953 + DENDRO_242*DENDRO_952 + DENDRO_247*DENDRO_948 + DENDRO_250*DENDRO_947 - DENDRO_334*grad2_2_2_beta0[pp] - DENDRO_336*grad2_1_1_beta0[pp] - 4.0/3.0*DENDRO_338*DENDRO_923 + DENDRO_415*DENDRO_942 + DENDRO_880*grad2_1_2_beta0[pp] - DENDRO_887*DENDRO_918 - DENDRO_894*DENDRO_920 - DENDRO_905*DENDRO_919 - DENDRO_921*DENDRO_922 - DENDRO_927*(DENDRO_29*DENDRO_930 + DENDRO_905*DENDRO_931) - DENDRO_927*(-DENDRO_46*DENDRO_928 + DENDRO_894*DENDRO_929) - DENDRO_927*(-DENDRO_49*DENDRO_924 + DENDRO_887*DENDRO_926) - DENDRO_932*DENDRO_933 - DENDRO_933*DENDRO_934 - DENDRO_935*DENDRO_936 - DENDRO_936*DENDRO_937 + DENDRO_938*DENDRO_940 + DENDRO_940*DENDRO_941 + DENDRO_943*DENDRO_944 + DENDRO_954;
                double DENDRO_956 = -DENDRO_896 - DENDRO_897 + DENDRO_898 - DENDRO_899 + DENDRO_900 + DENDRO_901 - DENDRO_902 + DENDRO_903 + DENDRO_904;
                double DENDRO_957 = -DENDRO_906 - DENDRO_907 + DENDRO_908 - DENDRO_909 + DENDRO_910 + DENDRO_911 + DENDRO_912 - DENDRO_913 + DENDRO_914;
                double DENDRO_958 = DENDRO_31*DENDRO_928;
                double DENDRO_959 = DENDRO_29*DENDRO_924;
                double DENDRO_960 = DENDRO_35*DENDRO_930;
                double DENDRO_961 = DENDRO_891*DENDRO_931;
                double DENDRO_962 = (1.0/3.0)*DENDRO_336;
                double DENDRO_963 = DENDRO_31*DENDRO_41;
                double DENDRO_964 = (1.0/3.0)*DENDRO_921;
                double DENDRO_965 = (1.0/3.0)*DENDRO_963;
                double DENDRO_966 = (7.0/3.0)*DENDRO_963;
                double DENDRO_967 = beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp];
                double DENDRO_968 = DENDRO_134*DENDRO_951 + DENDRO_219*DENDRO_948 - DENDRO_334*grad2_2_2_beta1[pp] - 4.0/3.0*DENDRO_336*DENDRO_938 - DENDRO_338*grad2_0_0_beta1[pp] + DENDRO_80*DENDRO_949 - DENDRO_882*grad2_0_2_beta1[pp] - DENDRO_891*DENDRO_919 + DENDRO_923*DENDRO_940 + DENDRO_932*DENDRO_944 + DENDRO_934*DENDRO_940 + DENDRO_935*DENDRO_966 + DENDRO_937*DENDRO_965 - DENDRO_941*DENDRO_962 - DENDRO_943*DENDRO_962 + DENDRO_963*DENDRO_964 + DENDRO_967;
                double DENDRO_969 = beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp];
                double DENDRO_970 = DENDRO_116*DENDRO_949 - DENDRO_12*DENDRO_411 + DENDRO_155*DENDRO_952 + DENDRO_160*DENDRO_951 + DENDRO_259*DENDRO_953 + DENDRO_265*DENDRO_948 + DENDRO_268*DENDRO_947 - 1.0/3.0*DENDRO_334*DENDRO_935 - 4.0/3.0*DENDRO_334*DENDRO_937 - DENDRO_334*DENDRO_964 - DENDRO_336*grad2_1_1_beta2[pp] - DENDRO_338*grad2_0_0_beta2[pp] - DENDRO_413*DENDRO_6 + DENDRO_413*DENDRO_942 - DENDRO_415*DENDRO_9 + DENDRO_883*grad2_0_1_beta2[pp] - DENDRO_893*DENDRO_920 - DENDRO_894*DENDRO_918 - DENDRO_915*DENDRO_919 - DENDRO_922*DENDRO_934 - DENDRO_923*DENDRO_936 - DENDRO_927*(DENDRO_31*DENDRO_930 + DENDRO_915*DENDRO_931) - DENDRO_927*(-DENDRO_46*DENDRO_924 + DENDRO_894*DENDRO_926) - DENDRO_927*(-DENDRO_55*DENDRO_928 + DENDRO_893*DENDRO_929) - DENDRO_932*DENDRO_936 + DENDRO_938*DENDRO_965 + DENDRO_941*DENDRO_966 + DENDRO_943*DENDRO_965 + DENDRO_969;
                double DENDRO_971 = 2*DENDRO_61;
                double DENDRO_972 = BSSN_ETA_R0*sqrt(DENDRO_41*(DENDRO_106*DENDRO_971 - DENDRO_274*DENDRO_55 - DENDRO_278*DENDRO_35 - DENDRO_280*DENDRO_49 - DENDRO_357*DENDRO_971 + 2*DENDRO_60*DENDRO_85))*pow(-pow(chi[pp], BSSN_ETA_POWER[0]) + 1, -BSSN_ETA_POWER[1]);

                // Dendro: printing variables
                //--
                a_rhs[pp] = -DENDRO_0*K[pp] + lambda[0]*(beta0[pp]*agrad_0_alpha[pp] + beta1[pp]*agrad_1_alpha[pp] + beta2[pp]*agrad_2_alpha[pp]);
                //--
                b_rhs0[pp] = B0[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta0[pp] + beta1[pp]*agrad_1_beta0[pp] + beta2[pp]*agrad_2_beta0[pp]);
                //--
                b_rhs1[pp] = B1[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta1[pp] + beta1[pp]*agrad_1_beta1[pp] + beta2[pp]*agrad_2_beta1[pp]);
                //--
                b_rhs2[pp] = B2[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta2[pp] + beta1[pp]*agrad_1_beta2[pp] + beta2[pp]*agrad_2_beta2[pp]);
                //--
                gt_rhs00[pp] = -At0[pp]*DENDRO_0 + DENDRO_10*DENDRO_9 + DENDRO_3*gt0[pp] - DENDRO_4*DENDRO_5 - DENDRO_5*DENDRO_6 + DENDRO_7*DENDRO_8 + beta0[pp]*agrad_0_gt0[pp] + beta1[pp]*agrad_1_gt0[pp] + beta2[pp]*agrad_2_gt0[pp];
                //--
                gt_rhs01[pp] = -At1[pp]*DENDRO_0 + DENDRO_11*gt0[pp] + DENDRO_12*gt2[pp] + DENDRO_13*DENDRO_2 + DENDRO_13*DENDRO_4 - DENDRO_14*gt1[pp] + DENDRO_7*gt3[pp] + DENDRO_9*gt4[pp] + beta0[pp]*agrad_0_gt1[pp] + beta1[pp]*agrad_1_gt1[pp] + beta2[pp]*agrad_2_gt1[pp];
                //--
                gt_rhs02[pp] = -At2[pp]*DENDRO_0 + DENDRO_15*gt0[pp] + DENDRO_16*gt1[pp] + DENDRO_17*DENDRO_2 + DENDRO_17*DENDRO_6 - DENDRO_18*gt2[pp] + DENDRO_7*gt4[pp] + DENDRO_9*gt5[pp] + beta0[pp]*agrad_0_gt2[pp] + beta1[pp]*agrad_1_gt2[pp] + beta2[pp]*agrad_2_gt2[pp];
                //--
                gt_rhs11[pp] = -At3[pp]*DENDRO_0 + DENDRO_11*DENDRO_8 + DENDRO_12*DENDRO_21 - DENDRO_14*gt3[pp] - DENDRO_19*gt3[pp] + DENDRO_20*gt3[pp] + beta0[pp]*agrad_0_gt3[pp] + beta1[pp]*agrad_1_gt3[pp] + beta2[pp]*agrad_2_gt3[pp];
                //--
                gt_rhs12[pp] = -At4[pp]*DENDRO_0 + DENDRO_11*gt2[pp] + DENDRO_12*gt5[pp] + DENDRO_15*gt1[pp] + DENDRO_16*gt3[pp] - DENDRO_19*gt4[pp] + DENDRO_22*DENDRO_4 + DENDRO_22*DENDRO_6 + beta0[pp]*agrad_0_gt4[pp] + beta1[pp]*agrad_1_gt4[pp] + beta2[pp]*agrad_2_gt4[pp];
                //--
                gt_rhs22[pp] = -At5[pp]*DENDRO_0 + DENDRO_10*DENDRO_15 + DENDRO_16*DENDRO_21 - DENDRO_18*gt5[pp] - DENDRO_19*gt5[pp] + DENDRO_23*gt5[pp] + beta0[pp]*agrad_0_gt5[pp] + beta1[pp]*agrad_1_gt5[pp] + beta2[pp]*agrad_2_gt5[pp];
                //--
                chi_rhs[pp] = -DENDRO_24*DENDRO_25 + DENDRO_24*K[pp]*alpha[pp] + beta0[pp]*agrad_0_chi[pp] + beta1[pp]*agrad_1_chi[pp] + beta2[pp]*agrad_2_chi[pp];
                //--
                At_rhs00[pp] = -At0[pp]*DENDRO_14 - At0[pp]*DENDRO_18 + At0[pp]*DENDRO_3 + DENDRO_26*DENDRO_7 + DENDRO_27*DENDRO_9 + DENDRO_802*(-DENDRO_113*(DENDRO_100 - DENDRO_102 - DENDRO_111 + DENDRO_98) + DENDRO_351*(-DENDRO_121*(DENDRO_115 + DENDRO_116*DENDRO_117 + DENDRO_118*DENDRO_80 + DENDRO_119*DENDRO_120) - DENDRO_146*(DENDRO_135 - DENDRO_143) - DENDRO_146*(DENDRO_156 + DENDRO_163) - DENDRO_146*(DENDRO_166 + DENDRO_167*DENDRO_169) - DENDRO_146*(DENDRO_175 + DENDRO_179*DENDRO_70) - DENDRO_146*(DENDRO_173*DENDRO_75 + DENDRO_181) + DENDRO_164*DENDRO_312 + DENDRO_184*(DENDRO_193 + 1.0*DENDRO_194) - DENDRO_184*(-DENDRO_169*DENDRO_76 + DENDRO_182*DENDRO_93) + DENDRO_187*DENDRO_314*DENDRO_320 - DENDRO_190*(DENDRO_186 + 1.0*DENDRO_188) - DENDRO_190*(DENDRO_164*DENDRO_191 + 1.0*DENDRO_167*DENDRO_93) + DENDRO_197*(DENDRO_192 + DENDRO_194) + DENDRO_197*(DENDRO_195 + DENDRO_196) - DENDRO_200*(DENDRO_185 + DENDRO_188) - DENDRO_200*(DENDRO_198 + DENDRO_199) - DENDRO_233*DENDRO_234 - DENDRO_253*DENDRO_254 - DENDRO_271*DENDRO_272 - DENDRO_299*DENDRO_300 + DENDRO_303*(DENDRO_150 + DENDRO_165) + DENDRO_307*(DENDRO_126 + DENDRO_304) + DENDRO_313*DENDRO_315*DENDRO_93 + DENDRO_316*DENDRO_317 + DENDRO_318*DENDRO_319 + DENDRO_350) - 12*DENDRO_58 + DENDRO_801*gt0[pp] + DENDRO_81*DENDRO_84 - DENDRO_96*(-DENDRO_89 + DENDRO_93)) - alpha[pp]*(-At0[pp]*K[pp] + DENDRO_42*(DENDRO_30 + DENDRO_32 - DENDRO_36) + DENDRO_52*(DENDRO_43 + DENDRO_48 - DENDRO_50) + DENDRO_57*(-At0[pp]*DENDRO_46 + DENDRO_53 - DENDRO_56)) + beta0[pp]*agrad_0_At0[pp] + beta1[pp]*agrad_1_At0[pp] + beta2[pp]*agrad_2_At0[pp];
                //--
                At_rhs01[pp] = At0[pp]*DENDRO_11 - At1[pp]*DENDRO_14 + At2[pp]*DENDRO_12 + At3[pp]*DENDRO_7 + At4[pp]*DENDRO_9 + DENDRO_2*DENDRO_803 + DENDRO_4*DENDRO_803 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_141*DENDRO_592 + DENDRO_155*DENDRO_590 + DENDRO_242*DENDRO_595 + DENDRO_762) - DENDRO_124*DENDRO_834 - DENDRO_129*DENDRO_836 + DENDRO_146*(-DENDRO_332*DENDRO_822 + DENDRO_752) - DENDRO_146*(-DENDRO_687 + DENDRO_747 + DENDRO_748) + DENDRO_146*(-DENDRO_173*DENDRO_688 + DENDRO_179*DENDRO_472 + DENDRO_725) + DENDRO_184*(DENDRO_823 + DENDRO_826) + DENDRO_184*(DENDRO_827 + DENDRO_828) + DENDRO_184*(DENDRO_187*DENDRO_470 + DENDRO_829) + DENDRO_184*(DENDRO_163 + DENDRO_742 + DENDRO_743) + DENDRO_190*(-DENDRO_70*DENDRO_822 + DENDRO_760) + DENDRO_190*(DENDRO_127*DENDRO_187 - DENDRO_173*DENDRO_327 + DENDRO_755) - DENDRO_190*(DENDRO_164*DENDRO_561 + DENDRO_164*DENDRO_576 + DENDRO_226*DENDRO_93) + DENDRO_190*(-DENDRO_164*DENDRO_688 - DENDRO_227*DENDRO_93 + DENDRO_753) - DENDRO_200*(DENDRO_122*DENDRO_187 + DENDRO_318) + DENDRO_301*(DENDRO_815 + DENDRO_819) + DENDRO_302*(-DENDRO_405 + DENDRO_566 + DENDRO_673) + DENDRO_306*(-DENDRO_733 + DENDRO_734) - DENDRO_306*(-DENDRO_164*DENDRO_437 + DENDRO_736 + DENDRO_821) - DENDRO_306*(-DENDRO_173*DENDRO_304 + DENDRO_739 + DENDRO_820) + DENDRO_315*(DENDRO_186 + DENDRO_187*DENDRO_70 + DENDRO_187*DENDRO_71) + DENDRO_315*(0.25*DENDRO_198 + 0.5*DENDRO_199 + DENDRO_470*DENDRO_93) - DENDRO_67*DENDRO_835 + DENDRO_767 + DENDRO_769 + DENDRO_770 + DENDRO_771 + DENDRO_772 + DENDRO_773 + DENDRO_774 + DENDRO_775 + DENDRO_776 + DENDRO_777 + DENDRO_778 + DENDRO_779 + DENDRO_780 - DENDRO_781*DENDRO_837 + DENDRO_786 - DENDRO_831*(DENDRO_122*DENDRO_179 + DENDRO_833) - DENDRO_831*(DENDRO_164*DENDRO_204 + DENDRO_168*DENDRO_202 + DENDRO_830)) - DENDRO_41*DENDRO_812*(DENDRO_164 - DENDRO_797) + DENDRO_781*DENDRO_800 - 12*DENDRO_787 - DENDRO_813*(-DENDRO_791 + DENDRO_792 + DENDRO_793 - DENDRO_795) + DENDRO_814*(DENDRO_64*(DENDRO_621 + DENDRO_63*DENDRO_781) + DENDRO_790)) - alpha[pp]*(-At1[pp]*K[pp] + DENDRO_42*DENDRO_806 + DENDRO_52*DENDRO_809 + DENDRO_57*DENDRO_811) + beta0[pp]*agrad_0_At1[pp] + beta1[pp]*agrad_1_At1[pp] + beta2[pp]*agrad_2_At1[pp];
                //--
                At_rhs02[pp] = At0[pp]*DENDRO_15 + At1[pp]*DENDRO_16 - At2[pp]*DENDRO_18 + At4[pp]*DENDRO_7 + At5[pp]*DENDRO_9 + DENDRO_2*DENDRO_838 + DENDRO_6*DENDRO_838 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_134*DENDRO_592 + DENDRO_160*DENDRO_590 + DENDRO_235*DENDRO_595 + DENDRO_587) - DENDRO_130*DENDRO_834 - DENDRO_146*(DENDRO_554 - DENDRO_555 + DENDRO_641) + DENDRO_146*(-DENDRO_179*DENDRO_561 - DENDRO_70*DENDRO_846 + DENDRO_854) - DENDRO_146*(DENDRO_332*DENDRO_849 + DENDRO_847 + DENDRO_848) + DENDRO_146*(DENDRO_389*DENDRO_824 - DENDRO_848 + DENDRO_853) - DENDRO_148*DENDRO_836 - DENDRO_184*(DENDRO_151*DENDRO_169 - DENDRO_850 - DENDRO_851) - DENDRO_184*(DENDRO_151*DENDRO_187 - DENDRO_179*DENDRO_191 - DENDRO_78*DENDRO_846) + DENDRO_184*(DENDRO_75*DENDRO_849 + DENDRO_850 + DENDRO_851) - DENDRO_190*(DENDRO_166 + DENDRO_826) - DENDRO_190*(DENDRO_828 + DENDRO_852) - DENDRO_190*(DENDRO_187*DENDRO_408 + DENDRO_829) - DENDRO_190*(DENDRO_166 + DENDRO_168*DENDRO_576 + DENDRO_825) + DENDRO_197*(DENDRO_132*DENDRO_187 + DENDRO_316) + DENDRO_302*(DENDRO_534 - DENDRO_638) - DENDRO_302*(-DENDRO_165*DENDRO_179 - DENDRO_75*DENDRO_846 + DENDRO_845) + DENDRO_305*(DENDRO_833 + DENDRO_844) + DENDRO_306*(DENDRO_164*DENDRO_214 + 0.25*DENDRO_830) + DENDRO_315*(0.25*DENDRO_195 + 1.0*DENDRO_196) + DENDRO_315*(DENDRO_187*DENDRO_75 + DENDRO_187*DENDRO_76 + DENDRO_193) + DENDRO_597 + DENDRO_599 + DENDRO_601 + DENDRO_602 + DENDRO_604 + DENDRO_606 + DENDRO_608 + DENDRO_609 + DENDRO_610 + DENDRO_611 + DENDRO_612 + DENDRO_613 + DENDRO_614 - DENDRO_615*DENDRO_837 + DENDRO_655 - DENDRO_72*DENDRO_835 - DENDRO_831*(DENDRO_819 + DENDRO_855)) + DENDRO_615*DENDRO_800 - 12*DENDRO_617 + DENDRO_634*DENDRO_814 - DENDRO_812*(DENDRO_618 + DENDRO_619 - DENDRO_620 - DENDRO_623) - DENDRO_813*(DENDRO_625 + DENDRO_626 - DENDRO_627 - DENDRO_630)) - alpha[pp]*(-At2[pp]*K[pp] + DENDRO_42*DENDRO_840 + DENDRO_52*DENDRO_842 + DENDRO_57*DENDRO_843) + beta0[pp]*agrad_0_At2[pp] + beta1[pp]*agrad_1_At2[pp] + beta2[pp]*agrad_2_At2[pp];
                //--
                At_rhs11[pp] = -At3[pp]*DENDRO_14 - At3[pp]*DENDRO_19 + At3[pp]*DENDRO_20 + DENDRO_11*DENDRO_26 + DENDRO_12*DENDRO_856 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_117*DENDRO_268 + DENDRO_118*DENDRO_230 + DENDRO_120*DENDRO_250 + DENDRO_456) - DENDRO_122*DENDRO_835 - DENDRO_146*(DENDRO_473 - 1.0*DENDRO_475) + DENDRO_146*(DENDRO_480 - 1.0*DENDRO_860) - DENDRO_146*(DENDRO_476*DENDRO_824 + DENDRO_478) + DENDRO_164*DENDRO_498 + DENDRO_173*DENDRO_497 + DENDRO_184*(DENDRO_466 - DENDRO_468) + DENDRO_184*(DENDRO_173*DENDRO_470 + DENDRO_180*DENDRO_817) + DENDRO_184*(DENDRO_470*DENDRO_824 + DENDRO_863) + DENDRO_184*(DENDRO_71*DENDRO_861 + DENDRO_862) + DENDRO_190*(DENDRO_483 + DENDRO_821) + DENDRO_190*(DENDRO_485 + DENDRO_820) + DENDRO_190*(DENDRO_487 - 1.0*DENDRO_864) + DENDRO_200*(DENDRO_486 - DENDRO_864) + DENDRO_200*(-DENDRO_122*DENDRO_173 + DENDRO_492) + DENDRO_200*(-DENDRO_164*DENDRO_202 + DENDRO_491) - DENDRO_202*DENDRO_836 + 6.0*DENDRO_221*DENDRO_305*DENDRO_822 - DENDRO_221*DENDRO_834 + DENDRO_410*(DENDRO_479 - DENDRO_860) + DENDRO_410*(-DENDRO_122*DENDRO_817 + DENDRO_489) + DENDRO_410*(-DENDRO_202*DENDRO_824 + DENDRO_490) - DENDRO_451*DENDRO_837 + DENDRO_493*(DENDRO_217 + DENDRO_404) + DENDRO_495*DENDRO_817 + DENDRO_496*(DENDRO_126 + DENDRO_457) + DENDRO_506) - 12*DENDRO_445 + DENDRO_801*gt3[pp] + DENDRO_83*(DENDRO_455 + DENDRO_64*(DENDRO_449 + DENDRO_452*DENDRO_63)) + DENDRO_859*(DENDRO_250 - DENDRO_446*(-DENDRO_106 + DENDRO_357 + DENDRO_358)) + DENDRO_96*(DENDRO_268 - DENDRO_446*(-DENDRO_85 + DENDRO_86 + DENDRO_87))) - alpha[pp]*(-At3[pp]*K[pp] + DENDRO_42*DENDRO_809 + DENDRO_806*DENDRO_857 + DENDRO_811*DENDRO_858) + beta0[pp]*agrad_0_At3[pp] + beta1[pp]*agrad_1_At3[pp] + beta2[pp]*agrad_2_At3[pp];
                //--
                At_rhs12[pp] = At1[pp]*DENDRO_15 + At2[pp]*DENDRO_11 + At3[pp]*DENDRO_16 - At4[pp]*DENDRO_19 + At5[pp]*DENDRO_12 + DENDRO_4*DENDRO_865 + DENDRO_6*DENDRO_865 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_208*DENDRO_592 + DENDRO_240*DENDRO_595 + DENDRO_259*DENDRO_590 + DENDRO_696) - DENDRO_128*DENDRO_835 + DENDRO_146*(DENDRO_218*DENDRO_822 + DENDRO_723) + DENDRO_146*(-DENDRO_226*DENDRO_849 + DENDRO_867) + DENDRO_146*(DENDRO_385*DENDRO_824 + DENDRO_867) - DENDRO_146*(DENDRO_123*DENDRO_846 + DENDRO_561*DENDRO_817 + DENDRO_722) + DENDRO_146*(-DENDRO_576*DENDRO_817 + DENDRO_677 - DENDRO_688*DENDRO_817) + DENDRO_183*(DENDRO_815 + DENDRO_818 + DENDRO_855) - DENDRO_184*(DENDRO_169*DENDRO_218 + DENDRO_853 - DENDRO_866) - DENDRO_184*(-DENDRO_191*DENDRO_817 - DENDRO_71*DENDRO_846 + DENDRO_854) + DENDRO_184*(DENDRO_470*DENDRO_849 + DENDRO_847 + DENDRO_866) + DENDRO_190*(-DENDRO_168*DENDRO_437 + DENDRO_868) + DENDRO_190*(-DENDRO_408*DENDRO_822 + DENDRO_728) + DENDRO_190*(-DENDRO_688*DENDRO_824 + DENDRO_868) + DENDRO_190*(-DENDRO_179*DENDRO_304 + DENDRO_725 - 0.25*DENDRO_844) + DENDRO_190*(-DENDRO_327*DENDRO_817 + DENDRO_685 - DENDRO_862) + DENDRO_190*(-DENDRO_470*DENDRO_822 + DENDRO_539 + DENDRO_691) - DENDRO_215*DENDRO_836 - DENDRO_224*DENDRO_834 + DENDRO_302*(DENDRO_659 - DENDRO_717) - DENDRO_302*(-DENDRO_165*DENDRO_817 + DENDRO_389*DENDRO_817 - DENDRO_470*DENDRO_846) - DENDRO_306*(-DENDRO_437*DENDRO_824 + DENDRO_718) - DENDRO_306*(-DENDRO_226*DENDRO_822 - DENDRO_227*DENDRO_822 + DENDRO_480) - DENDRO_306*(-DENDRO_304*DENDRO_817 + DENDRO_664 + DENDRO_719) + DENDRO_315*(DENDRO_147*DENDRO_164 + DENDRO_823) + DENDRO_315*(DENDRO_175 + DENDRO_827 + DENDRO_852) + DENDRO_410*(-DENDRO_204*DENDRO_822 + DENDRO_504) - DENDRO_701*DENDRO_837 + DENDRO_702 + DENDRO_731) - DENDRO_41*DENDRO_813*(-DENDRO_715 + DENDRO_817) + DENDRO_701*DENDRO_800 - 12*DENDRO_704 - DENDRO_812*(-DENDRO_705 + DENDRO_706 + DENDRO_707 - DENDRO_710) + DENDRO_814*(DENDRO_64*(DENDRO_628 + DENDRO_63*DENDRO_701) + DENDRO_713)) - alpha[pp]*(-At4[pp]*K[pp] + DENDRO_42*DENDRO_842 + DENDRO_840*DENDRO_857 + DENDRO_843*DENDRO_858) + beta0[pp]*agrad_0_At4[pp] + beta1[pp]*agrad_1_At4[pp] + beta2[pp]*agrad_2_At4[pp];
                //--
                At_rhs22[pp] = -At5[pp]*DENDRO_18 - At5[pp]*DENDRO_19 + At5[pp]*DENDRO_23 + DENDRO_15*DENDRO_27 + DENDRO_16*DENDRO_856 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_117*DENDRO_265 + DENDRO_118*DENDRO_219 + DENDRO_120*DENDRO_247 + DENDRO_374) + DENDRO_132*DENDRO_168*DENDRO_427 - DENDRO_132*DENDRO_835 - DENDRO_146*(DENDRO_395 - 1.0*DENDRO_397) - DENDRO_146*(0.25*DENDRO_869 + 1.0*DENDRO_870) + DENDRO_179*DENDRO_424 + DENDRO_184*(DENDRO_386 - 1.0*DENDRO_388) + DENDRO_184*(0.25*DENDRO_871 + 1.0*DENDRO_872) - DENDRO_184*(-DENDRO_391*DENDRO_846 + DENDRO_845) - DENDRO_190*(DENDRO_147*DENDRO_824 + DENDRO_848) - DENDRO_190*(DENDRO_168*DENDRO_214 + DENDRO_866) - DENDRO_190*(DENDRO_174*DENDRO_817 + DENDRO_179*DENDRO_408) - DENDRO_190*(DENDRO_76*DENDRO_861 + 0.25*DENDRO_816) + DENDRO_197*(DENDRO_871 + DENDRO_872) + DENDRO_197*(DENDRO_132*DENDRO_179 + DENDRO_74*DENDRO_846) + DENDRO_204*DENDRO_319*DENDRO_824 - DENDRO_204*DENDRO_834 + 6.0*DENDRO_211*DENDRO_301*DENDRO_849 - DENDRO_211*DENDRO_836 - DENDRO_367*DENDRO_837 - DENDRO_410*(DENDRO_869 + DENDRO_870) - DENDRO_410*(DENDRO_132*DENDRO_817 + DENDRO_167*DENDRO_846) + DENDRO_420*(DENDRO_150 + DENDRO_379) + DENDRO_421*(DENDRO_217 + DENDRO_376) + DENDRO_423*DENDRO_817 + DENDRO_444) - 12*DENDRO_352 + DENDRO_354*DENDRO_84 + DENDRO_801*gt5[pp] - DENDRO_859*(-DENDRO_360 + DENDRO_846) - DENDRO_95*(DENDRO_363 - DENDRO_364 + DENDRO_365 - DENDRO_372)) - alpha[pp]*(At5[pp]*DENDRO_51*DENDRO_843 - At5[pp]*K[pp] + DENDRO_57*DENDRO_842 + DENDRO_840*DENDRO_858) + beta0[pp]*agrad_0_At5[pp] + beta1[pp]*agrad_1_At5[pp] + beta2[pp]*agrad_2_At5[pp];
                //--
                K_rhs[pp] = -DENDRO_334*chi[pp]*(-DENDRO_352 + DENDRO_873*(DENDRO_353*DENDRO_450 + DENDRO_383) + DENDRO_875*(DENDRO_384 + DENDRO_874*gt5[pp]) + DENDRO_94*(DENDRO_375*DENDRO_41 + DENDRO_378*DENDRO_41 + DENDRO_381*DENDRO_41 - DENDRO_64*(DENDRO_366 - DENDRO_371))) - DENDRO_336*chi[pp]*(-DENDRO_445 + DENDRO_82*(DENDRO_41*DENDRO_460 + DENDRO_41*DENDRO_461 + DENDRO_454 - DENDRO_64*(DENDRO_448 - DENDRO_453)) + DENDRO_875*(DENDRO_462 + DENDRO_874*gt3[pp]) + DENDRO_877*(DENDRO_459 + DENDRO_876*gt3[pp])) - DENDRO_338*chi[pp]*(DENDRO_112*(DENDRO_102 + DENDRO_41*DENDRO_509 + DENDRO_41*DENDRO_510 - DENDRO_64*(DENDRO_103 - DENDRO_110)) - DENDRO_58 + DENDRO_873*(DENDRO_450*DENDRO_66 + DENDRO_508) + DENDRO_877*(DENDRO_507 + DENDRO_876*gt0[pp])) + DENDRO_880*chi[pp]*(-DENDRO_704 + 0.5*DENDRO_875*(DENDRO_108*DENDRO_714 + DENDRO_700) + DENDRO_878*(DENDRO_41*DENDRO_697 + DENDRO_41*DENDRO_698 - DENDRO_64*(DENDRO_59 - DENDRO_709) + DENDRO_705) + DENDRO_879*(DENDRO_41*DENDRO_699 - DENDRO_64*(DENDRO_60 - DENDRO_711) + DENDRO_712)) - DENDRO_882*chi[pp]*(-DENDRO_617 + 0.5*DENDRO_873*(DENDRO_450*DENDRO_633 + DENDRO_591) + DENDRO_878*(DENDRO_41*DENDRO_588 + DENDRO_41*DENDRO_589 + DENDRO_620 - DENDRO_64*(DENDRO_61 - DENDRO_622)) + DENDRO_881*(DENDRO_41*DENDRO_593 + DENDRO_41*DENDRO_594 + DENDRO_627 - DENDRO_64*(DENDRO_60 - DENDRO_629))) + DENDRO_883*chi[pp]*(-DENDRO_787 + 0.5*DENDRO_877*(DENDRO_369*DENDRO_796 + DENDRO_763) + DENDRO_879*(DENDRO_41*DENDRO_764 - DENDRO_64*(DENDRO_61 - DENDRO_788) + DENDRO_789) + DENDRO_881*(DENDRO_41*DENDRO_765 + DENDRO_41*DENDRO_766 - DENDRO_64*(DENDRO_59 - DENDRO_794) + DENDRO_791)) + DENDRO_916*(At0[pp]*DENDRO_887*DENDRO_888 + At1[pp]*DENDRO_895*DENDRO_905 + At2[pp]*DENDRO_894*DENDRO_895 + At3[pp]*DENDRO_888*DENDRO_891 + At4[pp]*DENDRO_895*DENDRO_915 + At5[pp]*DENDRO_888*DENDRO_893 + pow(K[pp], 2)) + beta0[pp]*agrad_0_K[pp] + beta1[pp]*agrad_1_K[pp] + beta2[pp]*agrad_2_K[pp];
                //--
                Gt_rhs0[pp] = DENDRO_955;
                //--
                Gt_rhs1[pp] = -DENDRO_141*DENDRO_950*DENDRO_956 + DENDRO_16*DENDRO_271 - DENDRO_208*DENDRO_950*DENDRO_957 + DENDRO_233*DENDRO_4 - DENDRO_233*DENDRO_942 + DENDRO_253*DENDRO_7 - DENDRO_822*DENDRO_947 + DENDRO_918*DENDRO_956 + DENDRO_920*DENDRO_957 + DENDRO_927*(DENDRO_960 - DENDRO_961) - DENDRO_927*(-DENDRO_926*DENDRO_956 + DENDRO_959) - DENDRO_927*(-DENDRO_929*DENDRO_957 + DENDRO_958) + DENDRO_968;
                //--
                Gt_rhs2[pp] = DENDRO_970;
                //--
                B_rhs0[pp] = -B0[pp]*DENDRO_972 - DENDRO_954*lambda[3] + DENDRO_955 + lambda[2]*(beta0[pp]*agrad_0_B0[pp] + beta1[pp]*agrad_1_B0[pp] + beta2[pp]*agrad_2_B0[pp]);
                //--
                B_rhs1[pp] = -B1[pp]*DENDRO_972 + DENDRO_141*DENDRO_952 - DENDRO_16*DENDRO_413 + DENDRO_208*DENDRO_953 + DENDRO_230*DENDRO_947 - DENDRO_4*DENDRO_411 + DENDRO_411*DENDRO_942 - DENDRO_415*DENDRO_7 - DENDRO_905*DENDRO_918 - DENDRO_915*DENDRO_920 - DENDRO_927*(-DENDRO_960 + DENDRO_961) - DENDRO_927*(DENDRO_905*DENDRO_926 + DENDRO_959) - DENDRO_927*(DENDRO_915*DENDRO_929 + DENDRO_958) - DENDRO_967*lambda[3] + DENDRO_968 + lambda[2]*(beta0[pp]*agrad_0_B1[pp] + beta1[pp]*agrad_1_B1[pp] + beta2[pp]*agrad_2_B1[pp]);
                //--
                B_rhs2[pp] = -B2[pp]*DENDRO_972 - DENDRO_969*lambda[3] + DENDRO_970 + lambda[2]*(beta0[pp]*agrad_0_B2[pp] + beta1[pp]*agrad_1_B2[pp] + beta2[pp]*agrad_2_B2[pp]);
                // Dendro: reduced ops: 4158
                // Dendro: }}} 
                //[[[end]]]
#else
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
                dendro.generate_cpu(outs, vnames, '[pp]')


                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 667753 
                // Dendro: printing temp variables
		  //std::cout<<"printer"<<std::endl;
             double DENDRO_240 = grad_2_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-grad_1_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double DENDRO_235 = -(grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double DENDRO_208 = -(grad_2_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+grad_1_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double DENDRO_173 = -(grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double DENDRO_824 = -(grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double DENDRO_349 = 4*pow(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp],-(2));
double DENDRO_595 = 0.5*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_0_chi[pp];
double DENDRO_592 = 0.5*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_1_chi[pp];
double DENDRO_590 = 0.5*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_2_chi[pp];
double DENDRO_42 = 2*At1[pp]*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_545 = 0.25*(-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp]);
double DENDRO_486 = (-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp])*grad_1_gt3[pp];
double DENDRO_57 = 2*At2[pp]*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_392 = (-((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp])-grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt5[pp];
double DENDRO_858 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*2*At4[pp];
double DENDRO_888 = 3*pow(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp],-(2));
double DENDRO_519 = (grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*grad_0_gt0[pp];
double DENDRO_848 = ((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.25*grad_1_gt5[pp];
double DENDRO_895 = 6*pow(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp],-(2));
double DENDRO_398 = grad_2_gt5[pp]*(grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double DENDRO_917 = 2*pow(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp],-(2));
double DENDRO_569 = grad_2_gt5[pp]*0.25*(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double DENDRO_469 = (grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*0.25*grad_1_gt5[pp];
double DENDRO_52 = At0[pp]*2*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_426 = grad_0_gt5[pp]*(-((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp])-grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double DENDRO_526 = (grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*grad_2_gt0[pp];
double DENDRO_84 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*12*grad_1_alpha[pp];
double DENDRO_650 = (-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp])*grad_2_gt3[pp];
double DENDRO_859 = 12*grad_0_alpha[pp]*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_525 = (grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*grad_1_gt0[pp];
double DENDRO_531 = grad_0_gt3[pp]*(-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp]);
double DENDRO_447 = 4*grad_2_alpha[pp]*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_515 = (grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*0.25*grad_0_gt5[pp];
double DENDRO_89 = 0.5*(1.0/chi[pp])*gt0[pp]*((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp]-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_0_chi[pp]-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad_2_chi[pp]);
double DENDRO_428 = 6.0*pow(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp],-(2));
double DENDRO_192 = (grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*grad_0_gt0[pp];
double DENDRO_452 = 0.5*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*gt3[pp];
double DENDRO_874 = 0.5*((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_0_chi[pp])*(1.0/chi[pp]);
double DENDRO_356 = 4*grad_1_alpha[pp]*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_850 = 0.25*grad_0_gt5[pp]*((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double DENDRO_656 = (-((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp])-grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt3[pp];
double DENDRO_316 = (grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*grad_2_gt0[pp];
double DENDRO_857 = At3[pp]*2*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_505 = grad_0_gt3[pp]*(-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp]);
double DENDRO_441 = (-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp])*grad_1_gt5[pp];
double DENDRO_362 = 4*grad_0_alpha[pp]*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_340 = (-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp])*grad_1_gt0[pp];
double DENDRO_96 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*12*grad_2_alpha[pp];
double DENDRO_700 = (-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_1_gt5[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+grad_2_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double DENDRO_425 = grad_1_gt5[pp]*(grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double DENDRO_342 = (-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp])*grad_1_gt0[pp];
double DENDRO_715 = ((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_chi[pp]-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_0_chi[pp])*(1.0/chi[pp])*gt4[pp];
double DENDRO_198 = (-(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*grad_2_gt0[pp];
double _915 = (grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*0.25*grad_2_gt0[pp];
double _2963 = (-(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*grad_1_gt5[pp];
double _1785 = (-((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp])-grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt0[pp];
double _3053 = -(grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3964 = -(-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp]);
double _3873 = 0.5*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_1_alpha[pp];
double _2483 = -(At1[pp])*2*alpha[pp]+grad_1_beta0[pp]*gt0[pp]+grad_1_beta2[pp]*gt2[pp]+(1.0/3.0)*gt1[pp]*grad_0_beta0[pp]+(1.0/3.0)*gt1[pp]*grad_1_beta1[pp]-(2.0/3.0)*grad_2_beta2[pp]*gt1[pp];
double _356 = -(-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp]);
double _700 = DENDRO_595*((grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_gt5[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_2_gt0[pp]);
double _3406 = (grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*DENDRO_590;
double _696 = ((grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp])*DENDRO_592;
double _958 = DENDRO_592*((grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp]+grad_0_gt3[pp]*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2)));
double _2965 = _2963+((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt3[pp];
double _730 = (-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*(1.0/chi[pp])*gt2[pp];
double _2959 = grad_0_gt3[pp]*(grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _2514 = -(At2[pp])*2*alpha[pp]+grad_2_beta0[pp]*gt0[pp]+grad_2_beta1[pp]*gt1[pp]+(1.0/3.0)*gt2[pp]*grad_0_beta0[pp]+(1.0/3.0)*gt2[pp]*grad_2_beta2[pp]-(2.0/3.0)*grad_1_beta1[pp]*gt2[pp];
double _3902 = 0.5*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_2_alpha[pp];
double _3186 = -(At2[pp])*K[pp]+DENDRO_42*(At2[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-At4[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _1420 = 0.5*(1.0/chi[pp])*gt3[pp]*((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp]-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_0_chi[pp]-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad_2_chi[pp]);
double _750 = grad_0_gt5[pp]*(grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _3188 = _3186+DENDRO_52*(-(At2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])+At4[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _1014 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*2.0*grad_2_alpha[pp];
double _565 = -(grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _3287 = -(-(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _2816 = At0[pp]*grad_1_beta0[pp]-At1[pp]*(2.0/3.0)*grad_2_beta2[pp]+At2[pp]*grad_1_beta2[pp]+At3[pp]*grad_0_beta1[pp]+At4[pp]*grad_0_beta2[pp]+grad_0_beta0[pp]*(1.0/3.0)*At1[pp]+grad_1_beta1[pp]*(1.0/3.0)*At1[pp];
double _3674 = ((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*grad_1_gt5[pp];
double _482 = -(-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp]);
double _3129 = 0.25*((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt0[pp];
double _874 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*2.0*grad_0_alpha[pp];
double _302 = 2*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_1_chi[pp];
double _2887 = -(-(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _2913 = -(-(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _3692 = grad_0_gt5[pp]*(grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _2825 = (grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*DENDRO_595;
double _1578 = (grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*0.5*grad_1_gt5[pp];
double _3557 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*6.0*grad_0_alpha[pp];
double _3875 = ((-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*(1.0/chi[pp])*gt2[pp];
double _756 = -(grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _3481 = -(grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _1998 = 0.25*(-(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _1958 = -(-(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _3397 = At1[pp]*grad_2_beta0[pp]+At2[pp]*grad_1_beta0[pp]+At3[pp]*grad_2_beta1[pp]-At4[pp]*(2.0/3.0)*grad_0_beta0[pp]+At5[pp]*grad_1_beta2[pp]+grad_1_beta1[pp]*(1.0/3.0)*At4[pp]+grad_2_beta2[pp]*(1.0/3.0)*At4[pp];
double _852 = DENDRO_590*(grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+grad_1_gt5[pp]*(pow(gt1[pp],2)-gt0[pp]*gt3[pp])+(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]));
double _2998 = -(At1[pp])*K[pp]+DENDRO_42*(At1[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-At3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _3000 = _2998+DENDRO_52*(-(At1[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])+At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-At4[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _1289 = 0.5*grad_0_gt5[pp]*(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _959 = DENDRO_590*((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(pow(gt1[pp],2)-gt0[pp]*gt3[pp])+grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+_958;
double _3039 = (-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp])*DENDRO_592;
double _3874 = _3875+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp];
double _799 = -(grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _2968 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*6.0*grad_2_alpha[pp];
double _3658 = 0.25*((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt5[pp];
double _320 = 2*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_0_chi[pp];
double _2822 = (-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp])*DENDRO_592;
double _2824 = _2822+(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*DENDRO_590;
double _3846 = 0.5*grad_0_alpha[pp]*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double _2125 = 2.0*pow(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp],-(3));
double _1483 = grad_0_gt3[pp]*(grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _3033 = At0[pp]*grad_2_beta0[pp]+At1[pp]*grad_2_beta1[pp]-At2[pp]*(2.0/3.0)*grad_1_beta1[pp]+At4[pp]*grad_0_beta1[pp]+At5[pp]*grad_0_beta2[pp]+grad_0_beta0[pp]*(1.0/3.0)*At2[pp]+grad_2_beta2[pp]*(1.0/3.0)*At2[pp];
double _2576 = -(At4[pp])*2*alpha[pp]+grad_1_beta0[pp]*gt2[pp]+grad_1_beta2[pp]*gt5[pp]+grad_2_beta0[pp]*gt1[pp]+grad_2_beta1[pp]*gt3[pp]-(2.0/3.0)*grad_0_beta0[pp]*gt4[pp]+(1.0/3.0)*gt4[pp]*grad_1_beta1[pp];
double _3040 = (-((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp])-grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_590;
double _1209 = DENDRO_526+grad_0_gt5[pp]*(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _3179 = _3188+DENDRO_57*(-(At5[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+-(At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))+At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _311 = 2*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_2_chi[pp];
double _2695 = -(0.5*((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _3805 = 1.0*grad_1_chi[pp]-((-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*DENDRO_452;
double _3748 = At5[pp]*2*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double _641 = -(-((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp])-grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _729 = -(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp]+_730;
double _2773 = -(At0[pp])*K[pp]+DENDRO_42*(At0[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+At2[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-At1[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _2777 = _2773+DENDRO_52*(At1[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+-(At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))-At0[pp]*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3848 = ((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_0_chi[pp])*(1.0/chi[pp])*gt4[pp];
double _913 = (grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*0.25*grad_2_gt3[pp];
double _2786 = alpha[pp]*(_2777+DENDRO_57*(-(At0[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+At1[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-At2[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])));
double _757 = (-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp])*grad_1_gt5[pp];
double _962 = DENDRO_595*(grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_1_gt0[pp]+(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]));
double _3003 = alpha[pp]*(_3000+DENDRO_57*(-(At1[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+At3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-At4[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])));
double _2023 = At0[pp]*pow(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp],2)+At3[pp]*pow(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp],2)+At5[pp]*pow(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp],2);
double _988 = -(-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp]);
double _3837 = 0.5*((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp])*(1.0/chi[pp])*gt0[pp];
double _1509 = grad_2_gt3[pp]*(grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _853 = grad_1_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+grad_2_gt3[pp]*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2));
double _3816 = 0.5*((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp])*(1.0/chi[pp])*gt3[pp];
double DENDRO_259 = grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double DENDRO_160 = -((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp])-grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double DENDRO_242 = grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double DENDRO_141 = -(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp];
double DENDRO_144 = pow(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp],-(2));
double DENDRO_179 = grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double DENDRO_164 = -(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _2969 = DENDRO_164-(1.0/chi[pp])*gt1[pp]*((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp]-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_0_chi[pp]-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad_2_chi[pp]);
double DENDRO_134 = -(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+grad_0_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_2_gt0[pp];
double DENDRO_118 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_1_chi[pp];
double DENDRO_117 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_2_chi[pp];
double DENDRO_120 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_0_chi[pp];
double DENDRO_781 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*gt1[pp];
double _2985 = -(grad_0_chi[pp])+(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*DENDRO_781;
double _1013 = -(grad_1_chi[pp])+((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_0_chi[pp])*DENDRO_781;
double _1023 = -(grad_0_chi[pp])+((-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*DENDRO_781;
double DENDRO_701 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*gt4[pp];
double _3856 = (1.0/chi[pp])*(grad_1_chi[pp]-((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp])*DENDRO_701);
double _883 = -(grad_2_chi[pp])+((-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*DENDRO_701;
double _3576 = -(grad_2_chi[pp])+(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*DENDRO_701;
double DENDRO_615 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*gt2[pp];
double _726 = -(grad_2_chi[pp])+((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_0_chi[pp])*DENDRO_615;
double _3884 = (1.0/chi[pp])*(grad_0_chi[pp]-((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp])*DENDRO_615);
double DENDRO_367 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*gt5[pp];
double _3792 = 1.0*grad_2_chi[pp]-DENDRO_367*0.5*((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp]);
double DENDRO_169 = 0.5*((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double DENDRO_877 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_2_alpha[pp];
double DENDRO_873 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*grad_1_alpha[pp];
double DENDRO_568 = 0.25*(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double DENDRO_875 = grad_0_alpha[pp]*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_451 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*gt3[pp];
double DENDRO_344 = 4.0*(1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double DENDRO_104 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*gt0[pp];
double DENDRO_110 = DENDRO_104*0.5*((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_0_chi[pp]);
double DENDRO_584 = (grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*grad_1_gt5[pp];
double DENDRO_634 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]))*_729;
double DENDRO_871 = ((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt5[pp];
double _723 = -(grad_0_chi[pp])+((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp])*DENDRO_615;
double _3588 = DENDRO_42*(-(At2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])+At4[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _990 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(DENDRO_208)*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))+DENDRO_141*0.25*grad_2_gt3[pp]+grad_1_gt3[pp]*DENDRO_545);
double _1487 = (grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*grad_2_gt3[pp];
double _2455 = -(At0[pp])*2*alpha[pp]+2*gt2[pp]*grad_0_beta2[pp]+(4.0/3.0)*grad_0_beta0[pp]*gt0[pp]-grad_1_beta1[pp]*(2.0/3.0)*gt0[pp]-(2.0/3.0)*gt0[pp]*grad_2_beta2[pp];
double _3632 = grad_0_gt5[pp]*((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double DENDRO_797 = (1.0/chi[pp])*gt1[pp]*((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp]-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_0_chi[pp]-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad_2_chi[pp]);
double _1016 = _1014*(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+DENDRO_797);
double _1169 = 0.25*grad_0_gt5[pp]*DENDRO_259+-(DENDRO_259*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])))+-(DENDRO_160*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])));
double _1645 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(0.25*grad_1_gt0[pp]*DENDRO_235+0.25*DENDRO_240*grad_0_gt0[pp]+0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])*DENDRO_242);
double _3754 = DENDRO_57*(-(At2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])+At4[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _3909 = grad_0_chi[pp]-((-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*DENDRO_781;
double DENDRO_945 = pow(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp],-(3));
double _2650 = -(At0[pp])*(2.0/3.0)*grad_2_beta2[pp]-At0[pp]*(2.0/3.0)*grad_1_beta1[pp]+At0[pp]*(4.0/3.0)*grad_0_beta0[pp]+2*At1[pp]*grad_0_beta1[pp]+2*At2[pp]*grad_0_beta2[pp];
double _3862 = grad_2_chi[pp]-((-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*DENDRO_701;
double _2052 = At0[pp]*pow(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp],2)+At3[pp]*pow(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp],2)+At5[pp]*pow(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp],2);
double _2038 = At0[pp]*pow(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp],2)+At3[pp]*pow(gt0[pp]*gt5[pp]-pow(gt2[pp],2),2)+At5[pp]*pow(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp],2);
double _1738 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(1.0)*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*DENDRO_134+0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_141);
double _3217 = -(At3[pp])*(2.0/3.0)*grad_2_beta2[pp]-At3[pp]*(2.0/3.0)*grad_0_beta0[pp]+At3[pp]*(4.0/3.0)*grad_1_beta1[pp]+grad_1_beta0[pp]*2*At1[pp]+grad_1_beta2[pp]*2*At4[pp];
double _1083 = DENDRO_160*0.25*grad_1_gt5[pp]+-(DENDRO_259*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])))+-(DENDRO_160*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])));
double _3592 = (-(At5[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+-(At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))+At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_858;
double _3892 = grad_2_chi[pp]-((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_0_chi[pp])*DENDRO_615;
double _3620 = -(At5[pp])*(2.0/3.0)*grad_1_beta1[pp]-At5[pp]*(2.0/3.0)*grad_0_beta0[pp]+At5[pp]*(4.0/3.0)*grad_2_beta2[pp]+grad_2_beta0[pp]*2*At2[pp]+grad_2_beta1[pp]*2*At4[pp];
double _3361 = DENDRO_42*(-(At1[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])+At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-At4[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]));
double _694 = DENDRO_590*((grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+grad_0_gt5[pp]*(pow(gt1[pp],2)-gt0[pp]*gt3[pp])+(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_gt0[pp]);
double DENDRO_527 = grad_0_gt5[pp]*(grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _3749 = _3748*(-(At5[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+-(At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))+At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3918 = grad_1_chi[pp]-((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_0_chi[pp])*DENDRO_781;
double _3463 = -((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3904 = ((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp])*(1.0/chi[pp])*gt1[pp];
double _3905 = _3902*(_3904+(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(pow(gt1[pp],2)-gt0[pp]*gt3[pp])+grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _1792 = (grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*grad_2_gt0[pp];
double _868 = -(grad_1_chi[pp])+((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp])*DENDRO_701;
double DENDRO_195 = ((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt0[pp];
double DENDRO_41 = (1.0/(-(gt2[pp]*gt4[pp])*2*gt1[pp]+pow(gt2[pp],2)*gt3[pp]+pow(gt4[pp],2)*gt0[pp]+pow(gt1[pp],2)*gt5[pp]-gt3[pp]*gt5[pp]*gt0[pp]));
double _1708 = -(DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])*DENDRO_41;
double _547 = -(0.5*grad_1_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*DENDRO_41-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41;
double _3803 = DENDRO_41*0.5*grad_1_gt3[pp]*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))+DENDRO_41*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_1_gt1[pp]+-(0.5*grad_0_gt3[pp]));
double _1833 = -(0.5*grad_2_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_41)+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41;
double _2655 = DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])*DENDRO_41;
double DENDRO_155 = grad_0_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_1_gt0[pp]-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _1458 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(-(1.0)*DENDRO_155*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))+0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259);
double _1745 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(1.0)*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*DENDRO_155+0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_160);
double DENDRO_168 = (gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double _1134 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.25*grad_2_gt0[pp]*DENDRO_242+0.25*DENDRO_240*grad_0_gt0[pp]+DENDRO_235*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _2216 = -(At1[pp]*pow(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp],2))-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double _2218 = _2216+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At4[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double DENDRO_453 = ((-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp])*DENDRO_452;
double _3339 = -((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp])+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_chi[pp]+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_0_chi[pp];
double _2544 = -(At3[pp])*2*alpha[pp]+grad_1_beta0[pp]*2*gt1[pp]+grad_1_beta2[pp]*2*gt4[pp]-(2.0/3.0)*grad_2_beta2[pp]*gt3[pp]-(2.0/3.0)*grad_0_beta0[pp]*gt3[pp];
double _248 = DENDRO_235*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-1.0*DENDRO_240*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-1.0*DENDRO_242*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _1411 = 4*grad_1_alpha[pp]*(_547+(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41+(1.0/chi[pp])*(-(1.0*grad_1_chi[pp])+DENDRO_453));
double _3334 = DENDRO_452*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp]);
double _265 = DENDRO_160*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-1.0*DENDRO_259*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-1.0*DENDRO_155*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _3920 = DENDRO_41*(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_1_gt0[pp]+DENDRO_41*(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])-(1.0/chi[pp])*_3918;
double _3446 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_169*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))+-(_1998*grad_2_gt5[pp])-0.25*grad_0_gt5[pp]*DENDRO_824);
double _2657 = -(12*grad_0_alpha[pp])*(_2655-(1.0/chi[pp])*(-(1.0*grad_0_chi[pp])+DENDRO_110)+DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*0.5*grad_0_gt0[pp]);
double _3336 = 12*grad_1_alpha[pp]*(_547+(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41+(1.0/chi[pp])*(-(1.0*grad_1_chi[pp])+_3334));
double _1350 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_141*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+DENDRO_208*0.25*grad_0_gt3[pp]+grad_1_gt3[pp]*DENDRO_545);
double _317 = -(1.0)*DENDRO_235*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+DENDRO_240*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+DENDRO_242*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _227 = DENDRO_134*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-1.0*DENDRO_208*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-1.0*DENDRO_141*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _299 = -(1.0)*DENDRO_134*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+DENDRO_208*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+DENDRO_141*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _555 = (-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*(1.0*grad_1_gt1[pp]+-(0.5*grad_0_gt3[pp]));
double _1538 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(DENDRO_259*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])))+DENDRO_160*0.25*grad_1_gt5[pp]+DENDRO_569);
double _552 = (-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(1.0*grad_1_gt1[pp]+-(0.5*grad_0_gt3[pp]))+(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(pow(gt1[pp],2)-gt0[pp]*gt3[pp]);
double DENDRO_371 = DENDRO_367*0.5*((-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_chi[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*grad_2_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp]);
double _1835 = 4*grad_2_alpha[pp]*(_1833-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41+(1.0/chi[pp])*(-(1.0*grad_2_chi[pp])+DENDRO_371));
double _258 = (0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double DENDRO_749 = -((0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])))-0.5*DENDRO_155*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]));
double _1355 = DENDRO_160*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+DENDRO_259*0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+0.5*grad_0_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3421 = 0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*DENDRO_824+-(0.25*grad_1_gt5[pp])*DENDRO_824+0.5*grad_2_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _1257 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_469+DENDRO_160*0.25*grad_2_gt3[pp]+0.5*grad_0_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _2841 = -(DENDRO_160*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp]))+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3289 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_3287*grad_2_gt3[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _3522 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(0.25*grad_2_gt3[pp])*DENDRO_824+1.0*grad_1_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _1222 = _799*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))+0.5*grad_2_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_259*0.25*grad_1_gt5[pp];
double _1086 = DENDRO_155*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+_913+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]);
double _1440 = DENDRO_259*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])-2*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3241 = 0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*DENDRO_824+2*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _1488 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_1487+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _1127 = 0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259+DENDRO_160*0.25*grad_2_gt3[pp]+0.5*grad_0_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _308 = -(1.0)*DENDRO_160*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+DENDRO_259*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+DENDRO_155*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _2672 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*((0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*DENDRO_155+-(0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_160));
double _1108 = 1.0*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*((-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_235+DENDRO_240*grad_2_gt0[pp]+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_235);
double _3741 = 0.5*grad_2_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_41-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41;
double _3743 = 12*grad_2_alpha[pp]*(_3741+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41-(1.0/chi[pp])*(-(1.0*grad_2_chi[pp])+DENDRO_371));
double _1145 = (grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_568+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*grad_2_gt0[pp]+DENDRO_155*0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]);
double _2603 = -(At5[pp])*2*alpha[pp]+2*gt2[pp]*grad_2_beta0[pp]+grad_2_beta1[pp]*2*gt4[pp]-(2.0/3.0)*grad_1_beta1[pp]*gt5[pp]-(2.0/3.0)*grad_0_beta0[pp]*gt5[pp];
double _3365 = (-(At1[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+At3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-At4[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*DENDRO_858;
double _2900 = -(DENDRO_160*0.25*grad_1_gt5[pp])+DENDRO_259*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))+DENDRO_160*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]));
double _241 = (0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])-(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _1506 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(grad_0_gt3[pp]*DENDRO_240+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double _3268 = 2*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_173*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]));
double _1092 = _565*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))+DENDRO_242*0.25*grad_0_gt3[pp]+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_1_gt0[pp];
double _1471 = _565*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))+2*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _3231 = DENDRO_240*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))-1.0*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3323 = (_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_349*(-(1.0*grad_1_gt1[pp])+-(0.5*grad_0_gt3[pp]));
double _1437 = -(DENDRO_240*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp])))+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _1148 = DENDRO_242*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_0_gt0[pp]+0.25*grad_1_gt0[pp]*DENDRO_242;
double _3283 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(-(grad_0_gt3[pp])*DENDRO_173+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*grad_1_gt0[pp]);
double _1184 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*grad_1_gt0[pp]+DENDRO_242*0.25*grad_0_gt3[pp]);
double DENDRO_686 = 0.5*DENDRO_240*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _1259 = 0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_242+DENDRO_235*0.25*grad_0_gt3[pp]+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_2_gt0[pp];
double _1121 = 0.25*grad_1_gt0[pp]*DENDRO_240+DENDRO_242*0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_2_gt0[pp];
double _3479 = _3481*0.25*grad_0_gt3[pp]+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_2_gt0[pp]-0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_173;
double _1248 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(-(DENDRO_160*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])))+0.25*grad_0_gt5[pp]*DENDRO_259+DENDRO_569);
double _1872 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_240*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])+0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_235);
double _1275 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(DENDRO_259*0.25*grad_2_gt3[pp]+1.0*grad_1_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _1384 = 0.5*grad_1_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_259*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+DENDRO_259*0.25*grad_2_gt3[pp];
double _3473 = -(0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp]))*DENDRO_824+0.5*grad_0_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))-DENDRO_164*0.25*grad_1_gt5[pp];
chi_rhs[pp] = -((2.0/3.0)*chi[pp])*(grad_0_beta0[pp]+grad_1_beta1[pp]+grad_2_beta2[pp])+(2.0/3.0)*chi[pp]*K[pp]*alpha[pp]+beta0[pp]*agrad_0_chi[pp]+beta1[pp]*agrad_1_chi[pp]+beta2[pp]*agrad_2_chi[pp];
double DENDRO_63 = -(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp];
double DENDRO_450 = (-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*grad_1_chi[pp]+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_0_chi[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_2_chi[pp];
double DENDRO_359 = (-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_chi[pp]-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_chi[pp]-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_0_chi[pp];
double _3347 = -((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*grad_1_chi[pp])+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_0_chi[pp]+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad_2_chi[pp];
double _3067 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*DENDRO_824-DENDRO_848+-(_1998*grad_2_gt5[pp]));
double _3590 = (At2[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-At4[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_857;
double _2225 = -(At4[pp]*pow(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp],2))-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At1[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _2227 = _2225+At0[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-At1[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _238 = -(0.5*grad_2_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _1891 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(grad_0_gt5[pp]*DENDRO_235+(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*grad_2_gt0[pp]);
double DENDRO_553 = 0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])*DENDRO_235+(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*0.5*grad_0_gt0[pp];
double _1603 = 0.25*grad_2_gt0[pp]*DENDRO_240+0.25*grad_0_gt5[pp]*DENDRO_242+(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*0.5*grad_1_gt0[pp];
double DENDRO_661 = 0.25*grad_0_gt5[pp]*DENDRO_240+(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp]);
double DENDRO_680 = 0.5*grad_0_gt3[pp]*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+DENDRO_240*0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]);
double _1228 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_680-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double DENDRO_563 = DENDRO_235*0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*0.5*grad_1_gt0[pp];
double DENDRO_420 = (_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_349;
double _1911 = grad_0_gt5[pp]*DENDRO_240+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _1828 = DENDRO_362*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])+0.5*(1.0/chi[pp])*gt5[pp]*DENDRO_359);
double _1710 = 4*grad_0_alpha[pp]*(_1708+(1.0/chi[pp])*(-(1.0*grad_0_chi[pp])+DENDRO_110)-DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*0.5*grad_0_gt0[pp]);
double _318 = _317-(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _1605 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(0.25*grad_2_gt0[pp]*DENDRO_240+DENDRO_563)+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1603;
double _3857 = DENDRO_41*grad_1_gt5[pp]*(pow(gt1[pp],2)-gt0[pp]*gt3[pp])+DENDRO_41*(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])-_3856;
double _1418 = _1411+DENDRO_362*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+DENDRO_359*0.5*(1.0/chi[pp])*gt3[pp])-4*grad2_1_1_alpha[pp];
double _1484 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_1483+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*grad_1_gt0[pp]);
double _255 = -(0.5*grad_2_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _1865 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.25*DENDRO_392+1.0*grad_0_gt5[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _1857 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(0.25*DENDRO_398+1.0*grad_1_gt5[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _1334 = DENDRO_160*0.25*grad_1_gt5[pp]+0.25*grad_0_gt5[pp]*DENDRO_259+(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp]);
double _1534 = (_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_160*0.25*grad_1_gt5[pp]+DENDRO_569;
double _1887 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_392+grad_0_gt5[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _1244 = (_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+0.25*grad_0_gt5[pp]*DENDRO_259+DENDRO_569;
double _1271 = -((0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])))+0.5*DENDRO_398;
double _1659 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*((_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_2_gt0[pp]+DENDRO_426);
double _1569 = -((0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])))+0.5*DENDRO_392;
double _1909 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_398+grad_1_gt5[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _1398 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(grad_2_gt3[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_425);
double _3506 = (0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))-0.5*DENDRO_398;
double _3824 = (-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])*DENDRO_41+DENDRO_41*(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*0.5*grad_0_gt0[pp];
double _3830 = grad_0_alpha[pp]*(_3824+DENDRO_41*(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])-(1.0/chi[pp])*(1.0*grad_0_chi[pp]-DENDRO_110))-grad2_0_0_alpha[pp];
double _249 = _248+(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _3363 = (At1[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-At3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*DENDRO_857;
double _1896 = grad_2_gt5[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_428;
double _1314 = -(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3107 = (0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))-0.5*DENDRO_392;
double _1214 = 1.0*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*((grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_208+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_208+DENDRO_650);
double _3313 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(grad_2_gt3[pp])*DENDRO_824+grad_1_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _2061 = At2[pp]*pow(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp],2)-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At4[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _2065 = _2061+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At4[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _780 = DENDRO_141*0.25*grad_2_gt3[pp]+DENDRO_141*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+-(DENDRO_208)*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]));
double _3789 = 0.5*grad_2_gt5[pp]*(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*DENDRO_41+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_2_gt4[pp]+-(0.5*grad_1_gt5[pp]))*DENDRO_41;
double _1451 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_240*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_242);
double _198 = -((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*0.5*grad_0_gt0[pp];
double _2875 = (0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))-DENDRO_173*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _3079 = (0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))-DENDRO_179*0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _2691 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.25*DENDRO_192+1.0*(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_2_gt0[pp]);
double _3094 = (_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_173*0.25*grad_2_gt0[pp]+DENDRO_179*0.25*grad_1_gt0[pp];
double _2722 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_173*grad_0_gt0[pp]+(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_1_gt0[pp]);
double _2858 = (_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+DENDRO_173*0.25*grad_2_gt0[pp]+DENDRO_179*0.25*grad_1_gt0[pp];
double _2895 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(grad_0_gt3[pp]*(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))+DENDRO_173*grad_1_gt0[pp]);
double _2876 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_2875+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_0_gt0[pp]);
double _2716 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_192+(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_2_gt0[pp]);
double _2701 = (_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*6.0*grad_0_gt0[pp];
double _1592 = 1.0*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*((-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_242+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_242+DENDRO_240*grad_1_gt0[pp]);
double _3250 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_155*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))-0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259);
double _1518 = (_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_349*(0.5*grad_0_gt3[pp]+1.0*grad_1_gt1[pp]);
double _2077 = At1[pp]*pow(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp],2)+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double _2079 = _2077-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At4[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _3756 = (At2[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-At4[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_858;
double _1361 = 0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+DENDRO_661;
double _493 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_134*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_208);
double _2669 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*((0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*DENDRO_134-0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_141);
double _3105 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(grad_0_gt5[pp]*(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))+DENDRO_316);
double _213 = -(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+0.5*grad_2_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double _500 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_441+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))));
double _1829 = -(4)*grad2_2_2_alpha[pp]+(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+0.5*(1.0/chi[pp])*gt5[pp]*DENDRO_63)*DENDRO_356+_1828;
double _1919 = (_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_349*(0.5*grad_1_gt5[pp]+1.0*grad_2_gt4[pp]);
double _1331 = DENDRO_134*0.25*grad_2_gt3[pp]+DENDRO_208*0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+0.5*grad_0_gt3[pp]*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _903 = (_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*0.5*grad_2_gt3[pp]+-(DENDRO_208)*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]));
double _1618 = (grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_545+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_545+(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*0.5*grad_1_gt0[pp];
double _529 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(grad_1_gt5[pp]*DENDRO_208+grad_2_gt3[pp]*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))));
double DENDRO_723 = DENDRO_208*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*0.5*grad_1_gt3[pp];
double _1328 = DENDRO_134*0.25*grad_2_gt3[pp]+DENDRO_141*0.25*grad_1_gt5[pp]+0.5*grad_0_gt3[pp]*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double DENDRO_536 = (_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+0.25*DENDRO_441;
double DENDRO_671 = DENDRO_208*0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+0.5*grad_0_gt3[pp]*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _1676 = DENDRO_208*0.25*grad_0_gt3[pp]+DENDRO_141*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+-(DENDRO_208)*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]));
double DENDRO_712 = grad_1_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41;
double _900 = -(4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*(-(DENDRO_141*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])))+DENDRO_671);
double _778 = 4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(-(DENDRO_134*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])))+DENDRO_536);
double _1661 = 1.0*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*grad_2_gt0[pp]+0.25*grad_0_gt5[pp]*DENDRO_235;
double _1697 = 1.0*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*((-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259+DENDRO_584);
double _437 = 0.5*grad_2_gt5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_2_gt2[pp]+-(0.5*grad_0_gt5[pp]));
double _3780 = -(grad2_2_2_alpha[pp])+DENDRO_873*(0.5*(1.0/chi[pp])*gt5[pp]*DENDRO_450+_437+(1.0*grad_2_gt4[pp]+-(0.5*grad_1_gt5[pp]))*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2)));
double DENDRO_571 = 0.25*grad_0_gt5[pp]*DENDRO_242+(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*0.5*grad_1_gt0[pp];
double DENDRO_756 = DENDRO_242*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_0_gt0[pp];
double _266 = _265+(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _907 = 1.0*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*((grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_141+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_141+DENDRO_531);
double _1996 = 0.5*grad_2_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _3081 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_3079-0.5*grad_0_gt0[pp]*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])));
double _3694 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_3692+grad_2_gt0[pp]*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])));
double _3738 = DENDRO_859*(-(0.5*(1.0/chi[pp])*gt5[pp]*DENDRO_359)+_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3056 = _3053*0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])-0.5*grad_1_gt0[pp]*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3666 = -(2*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3115 = -(0.25*grad_0_gt5[pp])*DENDRO_179-0.5*grad_2_gt0[pp]*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _153 = (-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*0.5*grad_0_gt0[pp];
double _1764 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.25*DENDRO_235*grad_0_gt0[pp]+1.0*(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_2_gt0[pp]);
double _1775 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(0.25*DENDRO_519+1.0*(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_1_gt0[pp]);
double _1052 = 0.25*grad_1_gt0[pp]*DENDRO_235+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+0.25*grad_2_gt0[pp]*DENDRO_242;
double _1650 = 0.25*DENDRO_240*grad_0_gt0[pp]+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+0.25*grad_1_gt0[pp]*DENDRO_235;
double _1142 = 0.25*DENDRO_240*grad_0_gt0[pp]+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+0.25*grad_2_gt0[pp]*DENDRO_242;
double _1553 = 0.25*grad_2_gt0[pp]*DENDRO_242+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+0.25*grad_1_gt0[pp]*DENDRO_235;
double _1789 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_519+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_1_gt0[pp]);
double _1567 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*((_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_0_gt5[pp]+DENDRO_235*grad_2_gt0[pp]);
double _1081 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*((_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_0_gt3[pp]+DENDRO_525);
double _1782 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_235*grad_0_gt0[pp]+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_2_gt0[pp]);
double _1731 = DENDRO_144*(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*6.0*grad_0_gt0[pp];
double _1686 = 0.5*DENDRO_235*grad_0_gt0[pp]+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double DENDRO_421 = (_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_349;
double _1353 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(-(DENDRO_235)*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))+DENDRO_686);
double _3265 = (grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+_1958*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp]);
double _587 = (_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_349*(1.0*grad_1_gt4[pp]+0.5*grad_2_gt3[pp]);
double DENDRO_550 = (grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_545+(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*0.5*grad_1_gt0[pp];
double _647 = -((0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))));
double _1608 = _1605+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_208+_647+0.5*DENDRO_134*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp]));
double _1401 = 1.0*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*((grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_160+DENDRO_527+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_160);
double _1311 = _1314+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*0.25*DENDRO_240+0.5*grad_0_gt5[pp]*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _300 = _299-(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double DENDRO_660 = -((0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double DENDRO_535 = -((0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _228 = _227+(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _1868 = (_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*2*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _1002 = 4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(DENDRO_134*0.5*grad_2_gt3[pp]+0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_208);
double _2085 = At4[pp]*pow(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp],2)+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At1[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]);
double _2087 = _2085-At0[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+At1[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _309 = _308-(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _1069 = -(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))+DENDRO_756;
double _3794 = grad_2_alpha[pp]*(_3789+(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(1.0*grad_2_gt2[pp]+-(0.5*grad_0_gt5[pp]))*DENDRO_41-(1.0/chi[pp])*_3792);
double _441 = (-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*0.5*grad_2_gt5[pp]+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*(1.0*grad_2_gt2[pp]+-(0.5*grad_0_gt5[pp]));
double _1543 = -(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))+DENDRO_553;
double _1999 = 0.5*grad_2_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _3661 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_3658+1.0*grad_0_gt5[pp]*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _3456 = 0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+_1998*grad_2_gt5[pp]+0.25*grad_0_gt5[pp]*DENDRO_824;
double _3711 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(grad_2_gt5[pp]*DENDRO_824+grad_1_gt5[pp]*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _3689 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_871+grad_0_gt5[pp]*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _3704 = 6.0*grad_2_gt5[pp]*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3059 = 0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+_1998*grad_2_gt5[pp]+DENDRO_848;
double _3730 = (_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+0.5*(1.0/chi[pp])*gt5[pp]*DENDRO_63)*DENDRO_84;
double _1269 = 4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(-(DENDRO_240*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])))+DENDRO_661);
double _3808 = -(grad2_1_1_alpha[pp])+grad_1_alpha[pp]*(_3803+(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41-(1.0/chi[pp])*_3805);
double _3774 = _3780+DENDRO_875*(_441+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_2_gt4[pp]+-(0.5*grad_1_gt5[pp]))+DENDRO_874*gt5[pp])+_3794;
double _152 = -((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*0.5*grad_0_gt0[pp])-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _1191 = DENDRO_155*0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+_915;
double _2863 = -(0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_160)+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]));
double _1691 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*grad_0_gt5[pp]+DENDRO_160*0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])+DENDRO_160*0.25*grad_2_gt0[pp];
double _1793 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_1792+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]));
double _1058 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*grad_1_gt5[pp]+0.25*grad_2_gt0[pp]*DENDRO_259+0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_160;
double _1777 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_155*0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _1800 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(1.0*grad_0_gt2[pp]+0.5*grad_2_gt0[pp])*DENDRO_349;
double _1556 = DENDRO_160*0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_515+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*grad_1_gt5[pp];
double _1654 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*grad_1_gt5[pp]+0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])*DENDRO_259;
double _1786 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_1785+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*grad_0_gt5[pp]);
double DENDRO_758 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*grad_2_gt3[pp]+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_568;
double _1153 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*((_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])+DENDRO_758);
double _1769 = -(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*2*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]));
double DENDRO_676 = _799*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))+0.5*grad_2_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _1468 = DENDRO_155*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _1118 = 0.25*grad_1_gt0[pp]*DENDRO_240+DENDRO_235*0.25*grad_0_gt3[pp]+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_2_gt0[pp];
double _1120 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259+DENDRO_749)+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1118;
double DENDRO_665 = DENDRO_240*0.25*grad_0_gt3[pp]+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]);
double _1278 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(DENDRO_665+0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double DENDRO_570 = (_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_569;
double DENDRO_683 = DENDRO_208*0.25*grad_2_gt3[pp]+(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*0.5*grad_1_gt3[pp];
double DENDRO_79 = (-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*0.5*grad_0_gt0[pp]+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _384 = (-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)*DENDRO_349*(1.0*grad_0_gt1[pp]+0.5*grad_1_gt0[pp]);
double _365 = (-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)+DENDRO_134*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _1199 = 0.5*grad_0_gt3[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)+DENDRO_141*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+0.25*DENDRO_340;
double _774 = 0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_141+0.25*grad_1_gt0[pp]*DENDRO_208+0.5*grad_2_gt3[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _377 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_342+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79));
double _380 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_340+grad_0_gt3[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79));
double DENDRO_784 = DENDRO_141*0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+0.5*grad_2_gt3[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _991 = DENDRO_134*0.25*grad_0_gt3[pp]+0.25*grad_1_gt0[pp]*DENDRO_208+0.5*grad_2_gt3[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _996 = _990-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_991-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_134*0.25*grad_0_gt3[pp]+DENDRO_784);
double DENDRO_546 = (-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_545+0.5*grad_1_gt5[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _3834 = _3830+DENDRO_873*(DENDRO_450*0.5*(1.0/chi[pp])*gt0[pp]+(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _1711 = DENDRO_356*(DENDRO_63*0.5*(1.0/chi[pp])*gt0[pp]+-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)+_1710;
double _765 = -(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)+DENDRO_550;
double _373 = -(2)*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double DENDRO_544 = 0.25*DENDRO_342+0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _2756 = (DENDRO_63*0.5*(1.0/chi[pp])*gt0[pp]+-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)*DENDRO_84;
double DENDRO_760 = DENDRO_141*0.25*grad_0_gt3[pp]+0.5*grad_1_gt3[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _1003 = _996+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_760+_988*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp])))-_1002;
double DENDRO_232 = (-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double DENDRO_789 = (grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*grad_1_gt0[pp]*DENDRO_41;
double _3913 = -(grad2_0_1_alpha[pp])+_3905+0.5*grad_1_alpha[pp]*(DENDRO_41*grad_0_gt3[pp]*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))-(1.0/chi[pp])*_3909+DENDRO_789);
double _3922 = 2*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41*chi[pp]*(_3913+0.5*grad_0_alpha[pp]*(_3920+grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41));
double _1371 = 0.5*grad_1_gt5[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+0.25*DENDRO_398+DENDRO_660;
double _486 = _482*4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(0.25*grad_1_gt2[pp]+-(0.25*grad_2_gt1[pp])+0.75*grad_0_gt4[pp]);
double DENDRO_564 = (_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]);
double DENDRO_750 = DENDRO_242*0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_2_gt0[pp];
double _1497 = DENDRO_240*4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(0.25*grad_0_gt4[pp]+0.75*grad_2_gt1[pp]+-(0.25*grad_1_gt2[pp]));
double _783 = 0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)+DENDRO_544;
double _497 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_208*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_134*0.25*grad_2_gt3[pp]);
double _220 = -(0.5*grad_1_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double _702 = grad_1_gt2[pp]*2.0*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _2727 = DENDRO_144*(_228+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_232)*2.0*grad_1_gt0[pp];
double _863 = grad_1_gt4[pp]*2.0*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _3504 = grad_1_gt4[pp]*2.0*DENDRO_144*(_228+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_232);
double _1495 = grad_1_gt3[pp]*2.0*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _2828 = grad_1_gt1[pp]*2.0*DENDRO_144*(_228+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_232);
double _1319 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(0.5*grad_1_gt5[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+DENDRO_683);
double _1502 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_208*grad_1_gt3[pp]+grad_2_gt3[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double _2139 = -(grad_1_beta0[pp])*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _1167 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*((_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*grad_1_gt0[pp]+DENDRO_505);
double _3045 = grad_1_gt2[pp]*2.0*DENDRO_144*(_228+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_232);
double _3225 = DENDRO_117*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_118*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _1480 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_486+grad_0_gt3[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double _1795 = 2.0*grad_1_gt0[pp]*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _1893 = grad_1_gt5[pp]*2.0*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _3298 = grad_1_gt3[pp]*2.0*DENDRO_144*(_228+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_232);
double _965 = grad_1_gt1[pp]*2.0*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _1294 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(grad_1_gt5[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+grad_2_gt3[pp]*DENDRO_208);
double _1358 = DENDRO_208*0.25*grad_0_gt3[pp]+grad_1_gt3[pp]*DENDRO_545+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]);
double _3699 = grad_1_gt5[pp]*2.0*DENDRO_144*(_228+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_232);
double _4035 = grad_1_beta1[pp]*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _1089 = 0.5*DENDRO_486+-((0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double _2303 = grad_1_beta2[pp]*DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _3972 = DENDRO_144*(_228+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_232)*grad_1_beta1[pp];
double _2906 = -(0.5*DENDRO_486)+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double DENDRO_294 = _302*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double _1230 = -(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+DENDRO_723;
double DENDRO_411 = DENDRO_144*(_300-(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_232);
double DENDRO_233 = DENDRO_144*(_228+(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_232);
double _1379 = 0.5*DENDRO_208*grad_1_gt3[pp]+(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double DENDRO_729 = (_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+grad_1_gt3[pp]*DENDRO_545;
double _1493 = grad_1_gt3[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_428;
double DENDRO_735 = -((0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double DENDRO_751 = (_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]);
double _1443 = 0.25*DENDRO_208*grad_1_gt3[pp]+1.0*grad_2_gt3[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3228 = -(2*(1.0/chi[pp]))*(_3225+DENDRO_120*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+-(grad2_1_1_chi[pp]));
double _1186 = 0.5*grad_0_gt3[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+0.25*DENDRO_486+DENDRO_735;
double _1282 = 0.5*grad_2_gt3[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+0.25*DENDRO_208*grad_1_gt3[pp];
double _1280 = _1282+(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _2915 = _2913*0.25*grad_2_gt3[pp]+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]);
double _3715 = (-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _1194 = 0.5*DENDRO_519+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _1338 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_570+-(DENDRO_160*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))));
double _1528 = _647+0.5*DENDRO_134*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+-(DENDRO_141*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])));
double _1530 = DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*((-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_235+_750+DENDRO_240*grad_2_gt0[pp])+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1528;
double _1429 = (1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+0.5*grad_1_gt3[pp]*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2));
double _1431 = DENDRO_117*(_552+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_118*(_1429+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_1_gt1[pp]+-(0.5*grad_0_gt3[pp])));
double _1434 = -(2*(1.0/chi[pp]))*(_1431+DENDRO_120*(_555+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+-(grad2_1_1_chi[pp]));
double _1727 = (-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*0.5*grad_0_gt0[pp];
double _1842 = 0.5*grad_2_gt5[pp]*(pow(gt1[pp],2)-gt0[pp]*gt3[pp])+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_2_gt4[pp]+-(0.5*grad_1_gt5[pp]));
double _3543 = DENDRO_173*0.25*grad_2_gt0[pp]+DENDRO_179*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_173*0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _3452 = (1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _1869 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_1868+-(DENDRO_235)*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])));
double _1616 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(DENDRO_259*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])))+DENDRO_570);
double _3541 = 4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(0.5*grad_0_gt5[pp]*DENDRO_164+DENDRO_168*0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp]));
double _1366 = 1.0*grad_2_gt3[pp]*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+0.25*grad_1_gt5[pp]*DENDRO_208;
double _363 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_141*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_134*0.25*grad_0_gt3[pp]);
double _364 = _356*4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(-(0.25*grad_0_gt4[pp])+0.25*grad_1_gt2[pp]+0.75*grad_2_gt1[pp])+_363;
double DENDRO_554 = (0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _3049 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_554-0.5*DENDRO_134*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+DENDRO_141*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])));
double _1161 = (_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _1160 = _1161+0.5*grad_1_gt3[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)+_988*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]));
double _1510 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_1509+grad_1_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])));
double _135 = (gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*0.5*grad_0_gt0[pp]+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _3098 = DENDRO_164*0.25*grad_0_gt5[pp]+DENDRO_168*0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+0.5*grad_1_gt5[pp]*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _2693 = _2695*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])+2*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _2708 = DENDRO_164*0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])+1.0*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _2725 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_198+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])));
double _2742 = (1.0*grad_0_gt2[pp]+0.5*grad_2_gt0[pp])*4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _2852 = DENDRO_168*0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+0.25*grad_2_gt0[pp]*DENDRO_824+0.5*grad_1_gt5[pp]*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _3132 = 4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(_3129+1.0*grad_0_gt5[pp]*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])));
double _3071 = (0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*DENDRO_169-DENDRO_850-0.5*grad_2_gt5[pp]*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _2890 = _2887*0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])-(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _2891 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_2890+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*grad_2_gt0[pp]);
double _3088 = DENDRO_164*0.25*grad_0_gt5[pp]+0.25*grad_2_gt0[pp]*DENDRO_824+0.5*grad_1_gt5[pp]*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _2719 = 2.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_195+grad_0_gt5[pp]*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])));
double _2936 = 0.25*DENDRO_198+0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _1321 = (1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_676;
double _2980 = -(grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41)+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]*DENDRO_41;
double _3116 = 4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(_3115+DENDRO_179*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])));
double _623 = (-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*0.5*grad_0_gt0[pp]+(pow(gt1[pp],2)-gt0[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _3839 = DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*chi[pp]*(_3834+DENDRO_877*(_623+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+_3837));
double _1389 = 0.25*grad_1_gt0[pp]*DENDRO_235+0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])*DENDRO_242+DENDRO_235*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _792 = (_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _1242 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_756*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))+_792+0.5*DENDRO_240*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]));
double _1326 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_235+_792+0.5*DENDRO_240*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]));
double _918 = -((_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])));
double _1140 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])*DENDRO_259+_918-0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*DENDRO_155);
double _1050 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_160+_918-0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*DENDRO_155);
double _1884 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_235*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+0.25*grad_2_gt0[pp]*DENDRO_240);
double _3245 = DENDRO_164*4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(0.25*grad_2_gt1[pp]+0.75*grad_0_gt4[pp]+-(0.25*grad_1_gt2[pp]));
double DENDRO_678 = (grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*0.25*DENDRO_240+0.5*grad_0_gt5[pp]*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _1226 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1222+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_678+DENDRO_240*0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]));
double _1377 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(DENDRO_240*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp])))+DENDRO_665);
double _1624 = (_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_546;
double DENDRO_560 = DENDRO_141*0.25*grad_1_gt5[pp]+0.5*grad_0_gt3[pp]*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _3461 = _3463*0.25*grad_2_gt3[pp]+0.5*grad_0_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))-DENDRO_164*0.25*grad_1_gt5[pp];
double _3639 = 1.0*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3515 = 0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _1125 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(DENDRO_208)*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))+DENDRO_729);
double _1129 = _1120+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1121+_1125+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1127;
double _1135 = _1129+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_141*0.25*grad_2_gt3[pp]+DENDRO_208*0.25*grad_0_gt3[pp]+DENDRO_751)-_1134;
double _1141 = _1135-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_208*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_784)-_1140;
double _1147 = _1141-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1142+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1145;
double _1041 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(DENDRO_235)*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))+DENDRO_750);
double _489 = 2*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp]);
double _501 = _486+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_489+-(DENDRO_208)*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])))+_493+_497-_500;
double _509 = _501-DENDRO_208*4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_1_gt4[pp]-0.25*grad_2_gt3[pp])-pow(grad_2_chi[pp],2)*pow(chi[pp],-(2))+4*gt2[pp]*grad_2_Gt0[pp];
double _515 = _509+4*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_gt5[pp]+2.0*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_gt5[pp];
double _521 = _515+2.0*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_gt5[pp]+2.0*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_gt5[pp];
double _532 = _521-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344*grad2_1_2_gt5[pp]-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344*grad2_0_1_gt5[pp]+_529+grad_2_Gt1[pp]*4*gt4[pp];
double _1474 = 0.25*DENDRO_486+1.0*grad_0_gt3[pp]*(_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _1671 = 0.5*grad_0_gt5[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+0.25*DENDRO_392+DENDRO_535;
double _820 = 4*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_gt4[pp]+2.0*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_gt4[pp];
double _826 = _820+2.0*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_gt4[pp]+2.0*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_gt4[pp];
double _836 = _826-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344*grad2_1_2_gt4[pp]-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344*grad2_0_1_gt4[pp]+grad_2_Gt2[pp]*2.0*gt4[pp]+grad_1_Gt0[pp]*2.0*gt2[pp];
double DENDRO_702 = _836+grad_1_Gt1[pp]*2.0*gt4[pp]+grad_1_Gt2[pp]*2.0*gt5[pp]-grad_1_chi[pp]*pow(chi[pp],-(2))*grad_2_chi[pp]+2.0*grad_2_Gt0[pp]*gt1[pp]+2.0*grad_2_Gt1[pp]*gt3[pp];
double _3642 = 0.25*grad_2_gt5[pp]*DENDRO_824+1.0*grad_1_gt5[pp]*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _1562 = DENDRO_515+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*grad_1_gt5[pp]+0.25*grad_2_gt0[pp]*DENDRO_259;
double _1448 = DENDRO_155*4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(0.25*grad_2_gt1[pp]+0.75*grad_0_gt4[pp]+-(0.25*grad_1_gt2[pp]));
double _1174 = 4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(DENDRO_240*0.5*grad_2_gt0[pp]+0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_235);
double _1253 = DENDRO_469+DENDRO_259*0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+0.5*grad_0_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double DENDRO_693 = (_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]);
double _1394 = 4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(0.5*grad_0_gt3[pp]*DENDRO_134+0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_141);
double _1285 = 0.25*grad_2_gt0[pp]*DENDRO_242+0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])*DENDRO_242+DENDRO_235*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _1666 = 0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+DENDRO_536;
double _1533 = _1530+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_563+_756*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])));
double _3125 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(DENDRO_164*0.5*grad_1_gt5[pp]+0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_824);
b_rhs2[pp] = B2[pp]*((3.0/4.0)*alpha[pp]*lambda_f[1]+(3.0/4.0)*lambda_f[0])+lambda[1]*(beta0[pp]*agrad_0_beta2[pp]+beta1[pp]*agrad_1_beta2[pp]+beta2[pp]*agrad_2_beta2[pp]);
b_rhs0[pp] = B0[pp]*((3.0/4.0)*alpha[pp]*lambda_f[1]+(3.0/4.0)*lambda_f[0])+lambda[1]*(beta0[pp]*agrad_0_beta0[pp]+beta1[pp]*agrad_1_beta0[pp]+beta2[pp]*agrad_2_beta0[pp]);
b_rhs1[pp] = B1[pp]*((3.0/4.0)*alpha[pp]*lambda_f[1]+(3.0/4.0)*lambda_f[0])+lambda[1]*(beta0[pp]*agrad_0_beta1[pp]+beta1[pp]*agrad_1_beta1[pp]+beta2[pp]*agrad_2_beta1[pp]);
double DENDRO_252 = (_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double _2168 = DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252)*(2.0/3.0)*(grad_0_beta0[pp]+grad_1_beta1[pp]+grad_2_beta2[pp]);
double _3230 = _3228-grad_0_gt3[pp]*2.0*DENDRO_144*(_249+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_252);
double _1436 = _1434+grad_0_gt3[pp]*2.0*DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252);
double _3635 = grad_0_gt5[pp]*2.0*DENDRO_144*(_249+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_252);
double _969 = 2.0*DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252)*grad_0_gt1[pp];
double _3156 = grad_0_gt2[pp]*2.0*DENDRO_144*(_249+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_252);
double _1849 = grad_0_gt5[pp]*2.0*DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252);
double _1797 = 2.0*grad_0_gt0[pp]*DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252);
double _2940 = grad_0_gt1[pp]*2.0*DENDRO_144*(_249+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_252);
double _2729 = DENDRO_144*(_249+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_252)*2.0*grad_0_gt0[pp];
double _3409 = grad_0_gt4[pp]*2.0*DENDRO_144*(_249+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_252);
double _859 = grad_0_gt4[pp]*2.0*DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252);
double _706 = 2.0*DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252)*grad_0_gt2[pp];
double _3976 = DENDRO_144*(_249+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_252)*grad_0_beta1[pp];
double _2144 = grad_0_beta0[pp]*DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252);
double _2339 = DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252)*grad_0_beta2[pp];
double _4039 = DENDRO_144*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252)*grad_0_beta1[pp];
double DENDRO_298 = _320*(_318-(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_252);
double _1442 = _1436+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1437+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1440;
double DENDRO_292 = 2*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(-(3*(1.0/chi[pp]))*grad_1_chi[pp]*grad_0_chi[pp]+2*grad2_0_1_chi[pp]);
double DENDRO_640 = DENDRO_160*0.25*grad_0_gt5[pp]+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*grad_2_gt5[pp];
double _1551 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*((_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*grad_2_gt0[pp]+DENDRO_640);
double DENDRO_647 = (_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]);
double _1560 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(0.25*DENDRO_240*grad_0_gt0[pp]+0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])*DENDRO_242+DENDRO_647);
double _1056 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_647+0.25*DENDRO_240*grad_0_gt0[pp]+DENDRO_235*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double DENDRO_270 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double _2337 = DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270)*(2.0/3.0)*(grad_0_beta0[pp]+grad_1_beta1[pp]+grad_2_beta2[pp]);
double _1802 = 2.0*grad_2_gt0[pp]*DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _861 = grad_2_gt4[pp]*2.0*DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _3291 = grad_2_gt3[pp]*2.0*DENDRO_144*(_266+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_270);
double _967 = grad_2_gt1[pp]*2.0*DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _3706 = grad_2_gt5[pp]*2.0*DENDRO_144*(_266+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_270);
double _1490 = grad_2_gt3[pp]*2.0*DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _3069 = grad_2_gt2[pp]*2.0*DENDRO_144*(_266+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_270);
double _2830 = grad_2_gt1[pp]*2.0*DENDRO_144*(_266+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_270);
double _3502 = grad_2_gt4[pp]*2.0*DENDRO_144*(_266+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_270);
double _1898 = grad_2_gt5[pp]*2.0*DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _2731 = DENDRO_144*(_266+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_270)*2.0*grad_2_gt0[pp];
double _704 = grad_2_gt2[pp]*2.0*DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _2142 = grad_2_beta0[pp]*DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _3967 = grad_2_beta1[pp]*DENDRO_144*(_266+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))+DENDRO_270);
double _4029 = grad_2_beta1[pp]*DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _2335 = DENDRO_144*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270)*grad_2_beta2[pp];
double DENDRO_296 = _311*(_309-(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-DENDRO_270);
double _1347 = DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*((grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_242+DENDRO_240*grad_1_gt0[pp]+grad_0_gt3[pp]*DENDRO_235);
double _1630 = (_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*grad_0_gt5[pp]+0.25*grad_2_gt0[pp]*DENDRO_235;
double _1631 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_1630+(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*0.5*grad_0_gt0[pp]);
double _2982 = 6.0*grad_0_alpha[pp]*(_2980+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41-(1.0/chi[pp])*_1013);
double _1266 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_693+DENDRO_141*0.25*grad_2_gt3[pp]+DENDRO_208*0.25*grad_0_gt3[pp]);
double _2927 = 0.25*DENDRO_173*grad_0_gt0[pp]+(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*grad_1_gt0[pp];
double _1588 = (_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*grad_2_gt0[pp]+0.25*DENDRO_235*grad_0_gt0[pp];
double _1759 = DENDRO_155*4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(0.25*grad_0_gt4[pp])+0.75*grad_1_gt2[pp]+0.25*grad_2_gt1[pp]);
double _1581 = 1.0*(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*grad_0_gt5[pp]+DENDRO_160*0.25*grad_2_gt0[pp];
double _1586 = _1588+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _745 = -(grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41)-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]*DENDRO_41;
double _747 = 2.0*grad_0_alpha[pp]*(_745+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41+(1.0/chi[pp])*_726);
double _1655 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_1654+DENDRO_160*0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]));
double _3667 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_3666+DENDRO_179*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])));
double DENDRO_495 = 4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(0.25*grad_0_gt4[pp]+0.75*grad_2_gt1[pp]+-(0.25*grad_1_gt2[pp]));
double _1959 = 0.5*grad_1_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))-(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double _3301 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_208*grad_1_gt3[pp]-grad_2_gt3[pp]*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double _3277 = 2.0*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_486-grad_0_gt3[pp]*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])));
double _3296 = 6.0*grad_1_gt3[pp]*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3411 = (0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+DENDRO_723;
double _3414 = -(2*(1.0/chi[pp]))*(DENDRO_208*DENDRO_592+DENDRO_240*DENDRO_595+_3406+-(grad2_1_2_chi[pp]))-_3409+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3411;
double _2835 = -(0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]))*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3470 = -(0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]))*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3498 = -(0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp]))*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3500 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_3498+DENDRO_141*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+grad_1_gt3[pp]*DENDRO_545);
double _3236 = 0.25*DENDRO_208*grad_1_gt3[pp]-1.0*grad_2_gt3[pp]*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _2866 = -(0.5*grad_1_gt0[pp])*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+DENDRO_760;
double _3240 = _3230-4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3231+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3236;
double _3271 = 0.25*DENDRO_486-1.0*grad_0_gt3[pp]*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3528 = (1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3524 = -(0.5*grad_2_gt3[pp])*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))-_3528+0.25*DENDRO_208*grad_1_gt3[pp];
double _3547 = -(grad_1_gt5[pp])*(_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+grad_2_gt3[pp]*DENDRO_208;
double _1717 = _1711+DENDRO_447*(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_89)-4*grad2_0_0_alpha[pp];
double _3571 = -(grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41)+grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_41;
double _3573 = 6.0*grad_2_alpha[pp]*(_3571+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41-(1.0/chi[pp])*_868);
double _2703 = 0.25*DENDRO_173*grad_0_gt0[pp]+1.0*(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*grad_1_gt0[pp];
double _1075 = (_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_1_gt0[pp]+DENDRO_760;
double DENDRO_397 = (grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _1852 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(DENDRO_240*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])))+DENDRO_397);
double _2686 = DENDRO_164*4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(0.25*grad_0_gt4[pp])+0.75*grad_1_gt2[pp]+0.25*grad_2_gt1[pp]);
double _2922 = -(DENDRO_173)*0.25*grad_0_gt3[pp]+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_1_gt0[pp];
double _1018 = grad_0_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_1_gt0[pp]*DENDRO_41;
double _1021 = _1016+2.0*grad_0_alpha[pp]*(_1018-(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41+(1.0/chi[pp])*_1013);
double DENDRO_798 = _1021+2.0*grad_1_alpha[pp]*((1.0/chi[pp])*_1023+-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*DENDRO_41+DENDRO_789)-4*grad2_0_1_alpha[pp];
double _1251 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_242+DENDRO_686);
double _2676 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_164*0.25*grad_0_gt5[pp]+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_169);
double _2925 = _2927+(_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _3057 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_3056+DENDRO_173*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])));
double _3063 = -(2*(1.0/chi[pp]))*(_3039+_3040+DENDRO_235*DENDRO_595+-(grad2_0_2_chi[pp]))-_3045-_3049+_3057-4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3059;
double _1611 = _1608+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_208+DENDRO_560);
double _1617 = _1611+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_564+DENDRO_160*0.25*grad_1_gt5[pp]+0.25*grad_0_gt5[pp]*DENDRO_259)+_1616;
double _1623 = _1617-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1618-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.25*grad_2_gt0[pp]*DENDRO_235+DENDRO_553);
double _3640 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_240*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))-_3639);
double _2923 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(_2922+DENDRO_173*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp])));
double _3471 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_3470+DENDRO_141*0.25*grad_2_gt3[pp]+DENDRO_208*0.25*grad_0_gt3[pp]);
double _2842 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_2841+0.5*DENDRO_155*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])));
double _1634 = (_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _1633 = _1634+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*grad_2_gt5[pp]+_641*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]));
double _1637 = _1623-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1624-_1631-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1633;
double _3653 = 1.0*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _3654 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_134*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))-_3653);
double _1035 = DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*((-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_242+DENDRO_240*grad_1_gt0[pp]+grad_0_gt3[pp]*DENDRO_235);
double _1042 = _1035+DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*((-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259+DENDRO_584+DENDRO_656)+_1041;
double _1045 = _1042+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_160*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+DENDRO_749);
double _3859 = -(grad2_1_2_alpha[pp])+_3846*(_3848+DENDRO_700)+0.5*grad_2_alpha[pp]*(_3857+grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41);
double _3842 = _3859+0.5*grad_1_alpha[pp]*(DENDRO_41*grad_2_gt3[pp]*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))-(1.0/chi[pp])*_3862+DENDRO_712);
double _1541 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_641*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))+DENDRO_640);
double _1552 = _1533+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1534+_1538-_1541-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1543-_1551;
double _1561 = _1552+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1553+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1556+_1560;
double _1571 = _1561+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1562-_1567-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1569;
double _374 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_373+DENDRO_141*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _392 = _364-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_365+_374-_377+_380-_384-pow(grad_0_chi[pp],2)*pow(chi[pp],-(2))+4*grad_0_Gt0[pp]*gt0[pp];
double _3259 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_824+DENDRO_164*0.25*grad_1_gt5[pp]);
double _1681 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(DENDRO_240*0.5*grad_1_gt0[pp]+0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_242);
double _1004 = 1.0*grad_0_gt3[pp]*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)+0.25*DENDRO_340;
double _1051 = _1045+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_751+DENDRO_141*0.25*grad_2_gt3[pp]+grad_1_gt3[pp]*DENDRO_545)-_1050;
double _1061 = _1051-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1052-_1056-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1058;
double _1064 = _1061+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_155*0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_758);
double _1902 = DENDRO_240*4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(0.25*grad_0_gt4[pp]+0.75*grad_1_gt2[pp]+-(0.25*grad_2_gt1[pp]));
double _738 = -((gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]*DENDRO_41)-grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_41;
double _748 = -(4)*grad2_0_2_alpha[pp]+2.0*grad_2_alpha[pp]*(_738+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41+(1.0/chi[pp])*_723)+_747;
double _1684 = 4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(DENDRO_134*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_544);
double _1263 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_141*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+DENDRO_729);
double _1462 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_242*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+0.25*grad_1_gt0[pp]*DENDRO_240);
double _2864 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_2863+0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*DENDRO_155);
double DENDRO_741 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+_915;
double _1096 = 0.5*(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_741;
double _864 = -(2*(1.0/chi[pp]))*(_852+DENDRO_592*_853+DENDRO_595*DENDRO_700+-(grad2_1_2_chi[pp]))+_859+_861+_863;
double _660 = (1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _1647 = DENDRO_208*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+-(0.5)*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*DENDRO_134+_660;
double _771 = 0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_141+-(0.5)*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*DENDRO_134+_660;
double _3344 = DENDRO_859*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])-0.5*(1.0/chi[pp])*gt3[pp]*_3339);
double _399 = _392+grad_0_Gt1[pp]*4*gt1[pp]+grad_0_Gt2[pp]*4*gt2[pp]+4*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_gt0[pp];
double _405 = _399+2.0*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_gt0[pp]+2.0*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_gt0[pp];
double _411 = _405+2.0*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_gt0[pp]-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344*grad2_1_2_gt0[pp];
double _2931 = _2936+0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _1964 = (_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]);
double DENDRO_288 = 2*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(-(3*(1.0/chi[pp])*grad_2_chi[pp])*grad_0_chi[pp]+2*grad2_0_2_chi[pp]);
double _3318 = DENDRO_259*4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(-(1.0*grad_2_gt4[pp])+0.25*grad_1_gt5[pp]);
double _2916 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(_2915+_1958*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp]));
double _1499 = DENDRO_242*4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(-(0.25*grad_1_gt0[pp])+1.0*grad_0_gt1[pp]);
double _1290 = 4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(_1289+0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_160);
double DENDRO_719 = 0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _1307 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_240*0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_680);
double _1310 = _1307+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_208*0.5*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp])+DENDRO_683);
double _1327 = _1310+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1311+_1319+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1321-_1326;
double _1333 = _1327-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1328-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1331;
double _1066 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(1.0*grad_1_gt4[pp]-0.5*grad_2_gt3[pp]);
double _1065 = _1066+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_568+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*grad_2_gt0[pp];
double _1074 = _1064+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1065+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1069;
double _1085 = _1074+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1075+_1081-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1083;
double _1091 = _1085-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_1086-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_1089;
double _1100 = _1091-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_1092-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1096;
double DENDRO_388 = (grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _1860 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(-(DENDRO_134*0.5*(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp])))+DENDRO_388);
double _1206 = 1.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*((grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_141+DENDRO_208*grad_1_gt0[pp]+DENDRO_531);
double DENDRO_740 = (_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]);
double _1101 = (_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*grad_1_gt0[pp]+0.25*DENDRO_519+DENDRO_740;
double _3649 = DENDRO_179*4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(-(0.25*grad_2_gt0[pp])+1.0*grad_0_gt2[pp]);
double _3890 = DENDRO_41*(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_0_gt5[pp]+DENDRO_41*(pow(gt4[pp],2)-gt3[pp]*gt5[pp])*grad_2_gt0[pp];
double _3894 = 0.5*grad_0_alpha[pp]*(_3890+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41-(1.0/chi[pp])*_3892);
double _2849 = -(DENDRO_173)*0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+DENDRO_179*0.5*(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]));
double _2850 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_2849+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_2_gt0[pp]);
double _1579 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(_1578+0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259);
double _1640 = _1637+DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_526+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_160+DENDRO_527);
double _1646 = _1640+DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_208*grad_1_gt0[pp]+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_141+DENDRO_531)+_1645;
double _1660 = _1646+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1647+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1650+_1655-_1659;
double _1670 = _1660-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1661-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1666;
double _1685 = _1670-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1671-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_1676-_1681-_1684;
double _1698 = _1685-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1686-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1691-_1697;
double _3352 = DENDRO_96*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])-0.5*(1.0/chi[pp])*gt3[pp]*_3347);
double _885 = 2.0*grad_1_alpha[pp]*((1.0/chi[pp])*_883+-(grad_2_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*DENDRO_41+DENDRO_712);
double _3169 = (gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad_2_gt0[pp]*DENDRO_41+grad_0_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_41;
double _3171 = 6.0*grad_2_alpha[pp]*(_3169-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41-(1.0/chi[pp])*_723);
double _3426 = 0.5*grad_0_gt3[pp]*(_1996-(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _1749 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_160*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_515);
double _580 = -(4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*(0.5*grad_0_gt3[pp]*DENDRO_208+DENDRO_141*0.25*grad_2_gt3[pp]);
double _590 = _580-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_141*0.5*grad_2_gt3[pp]+DENDRO_208*0.25*grad_0_gt3[pp])-_587-pow(grad_1_chi[pp],2)*pow(chi[pp],-(2));
double _597 = _590-3.0*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad_2_gt3[pp]*DENDRO_208+4*gt1[pp]*grad_1_Gt0[pp]+4*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_gt3[pp];
double _603 = _597+2.0*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_gt3[pp]+2.0*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_gt3[pp];
double _609 = _603+2.0*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_gt3[pp]-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344*grad2_1_2_gt3[pp];
double _616 = _609-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344*grad2_0_1_gt3[pp]-3.0*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*DENDRO_505+4*gt4[pp]*grad_1_Gt2[pp];
double _3085 = 0.5*grad_2_gt0[pp]*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_850;
double _3086 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_3085+0.5*grad_2_gt5[pp]*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])));
double _3090 = _3063+_3067-_3069-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_3071-_3081+_3086-4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3088;
double _908 = _900-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(_903+0.25*grad_1_gt5[pp]*DENDRO_208)-_907;
double DENDRO_731 = _908-1.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*((grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_208+_757+DENDRO_650);
double _1880 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_160*0.5*grad_1_gt5[pp]+0.25*grad_0_gt5[pp]*DENDRO_259);
double _2736 = DENDRO_160*4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(-(1.0*grad_2_gt2[pp])+0.25*grad_0_gt5[pp]);
double _459 = -((-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(2*grad2_2_2_chi[pp]-pow(grad_2_chi[pp],2)*3*(1.0/chi[pp])));
double _460 = _459-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(3*(1.0/chi[pp]))*pow(grad_1_chi[pp],2)+2*grad2_1_1_chi[pp]);
double DENDRO_867 = -(0.25*grad_1_gt5[pp])*DENDRO_824+0.5*grad_2_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3419 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(0.5*grad_2_gt3[pp])*(_1999+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_867);
double _1770 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_1769+DENDRO_160*0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]));
double _3247 = DENDRO_173*4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(-(0.25*grad_1_gt0[pp])+1.0*grad_0_gt1[pp]);
double _3174 = grad_0_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41+(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad_2_gt0[pp]*DENDRO_41;
double _3176 = 6.0*grad_0_alpha[pp]*(_3174-(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41-(1.0/chi[pp])*_726);
double _3136 = (_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _3134 = (_198+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*grad_2_gt0[pp]+_3136+0.25*DENDRO_192;
double DENDRO_737 = _913+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*0.5*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]);
double _1179 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_737);
double DENDRO_722 = (0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3579 = 6.0*grad_1_alpha[pp]*((1.0/chi[pp])*_3576+-(grad_2_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*DENDRO_41+DENDRO_712);
double _1900 = DENDRO_235*4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(-(0.25*grad_2_gt0[pp])+1.0*grad_0_gt2[pp]);
double DENDRO_423 = 4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(0.25*grad_0_gt4[pp]+0.75*grad_1_gt2[pp]+-(0.25*grad_2_gt1[pp]));
double _2988 = 6.0*grad_1_alpha[pp]*((1.0/chi[pp])*_2985+-(grad_0_gt3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*DENDRO_41+DENDRO_789);
double _1298 = 1.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_750+DENDRO_240*grad_2_gt0[pp]+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_235);
double _461 = _460-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(2*grad2_0_0_chi[pp]-3*(1.0/chi[pp])*pow(grad_0_chi[pp],2));
double DENDRO_419 = (_461+2*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(2*grad2_1_2_chi[pp]-3*(1.0/chi[pp])*grad_2_chi[pp]*grad_1_chi[pp])-DENDRO_288+DENDRO_292-DENDRO_294-DENDRO_296-DENDRO_298)*(1.0/chi[pp]);
double _711 = -(2*(1.0/chi[pp]))*(-(grad2_0_2_chi[pp])+_694+_696+_700)+_702+_704+_706+DENDRO_419*DENDRO_615+2.0*grad_2_Gt0[pp]*gt0[pp]+2.0*grad_2_Gt1[pp]*gt1[pp];
double _716 = _711+grad_0_Gt0[pp]*2.0*gt2[pp]+grad_2_Gt2[pp]*2.0*gt2[pp]+grad_0_Gt1[pp]*2.0*gt4[pp]+grad_0_Gt2[pp]*2.0*gt5[pp]+-(pow(chi[pp],-(2))*grad_2_chi[pp])*grad_0_chi[pp];
double _974 = -(2*(1.0/chi[pp]))*(_959+_962+-(grad2_0_1_chi[pp]))+_965+_967+_969+DENDRO_419*DENDRO_781+2.0*grad_1_Gt0[pp]*gt0[pp]+grad_0_Gt0[pp]*2.0*gt1[pp];
double _979 = _974+grad_1_Gt1[pp]*2.0*gt1[pp]+grad_1_Gt2[pp]*2.0*gt2[pp]+2.0*grad_0_Gt1[pp]*gt3[pp]+grad_0_Gt2[pp]*2.0*gt4[pp]+-(grad_1_chi[pp]*grad_0_chi[pp])*pow(chi[pp],-(2));
double DENDRO_418 = _461+2*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(2*grad2_1_2_chi[pp]-3*(1.0/chi[pp])*grad_2_chi[pp]*grad_1_chi[pp])-DENDRO_288+DENDRO_292-DENDRO_294-DENDRO_296-DENDRO_298;
double _981 = _979+4*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_gt1[pp]+2.0*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_gt1[pp];
double _983 = _981+2.0*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_gt1[pp]+2.0*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_gt1[pp];
double DENDRO_782 = _983+-((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344)*grad2_1_2_gt1[pp]+-((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344)*grad2_0_1_gt1[pp];
double _1032 = _1100-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1101-_1108+DENDRO_782+_1003-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1004;
double _718 = _716+4*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_gt2[pp]+2.0*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_gt2[pp];
double _720 = _718+2.0*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_gt2[pp]+2.0*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_gt2[pp];
double _1700 = alpha[pp]*(_1698+_720+-((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344)*grad2_1_2_gt2[pp]+-((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344)*grad2_0_1_gt2[pp]);
double DENDRO_616 = _720+-((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344)*grad2_1_2_gt2[pp]+-((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344)*grad2_0_1_gt2[pp];
double _1876 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(0.5*grad_0_gt5[pp]*DENDRO_259+DENDRO_160*0.25*grad_1_gt5[pp]);
double _878 = grad_2_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41-grad_1_gt5[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_41;
double _886 = _874*(DENDRO_240+DENDRO_715)+2.0*grad_2_alpha[pp]*(_878-(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41+(1.0/chi[pp])*_868)+_885;
double _1573 = (_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*0.5*grad_2_gt0[pp];
double _1576 = _1571-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(_1573+-(DENDRO_235)*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))+0.25*grad_0_gt5[pp]*DENDRO_235);
double _1594 = _1576-_1579-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1581-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1586-_1592+DENDRO_616;
double _3881 = DENDRO_41*grad_0_gt5[pp]*(pow(gt1[pp],2)-gt0[pp]*gt3[pp])+DENDRO_41*(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*grad_2_gt0[pp];
double _3869 = -(grad2_0_2_alpha[pp])+_3873*_3874+0.5*grad_2_alpha[pp]*(_3881+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41-_3884)+_3894;
double _3813 = _3808+DENDRO_875*(_555+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])+DENDRO_874*gt3[pp]);
double _3818 = (gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*chi[pp]*(_3813+DENDRO_877*(_552+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+_3816));
double _3867 = -(DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*chi[pp]*_3774-_3818-_3839+2*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41*chi[pp]*_3842;
double DENDRO_281 = (-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(2*grad2_0_0_chi[pp]-3*(1.0/chi[pp])*pow(grad_0_chi[pp],2));
double _1844 = DENDRO_118*(_437+(1.0*grad_2_gt4[pp]+-(0.5*grad_1_gt5[pp]))*(-(gt0[pp]*gt5[pp])+pow(gt2[pp],2)));
double _1845 = DENDRO_117*(_1842+(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(1.0*grad_2_gt2[pp]+-(0.5*grad_0_gt5[pp])))+_1844;
double _1850 = -(2*(1.0/chi[pp]))*(_1845+DENDRO_120*(_441+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(1.0*grad_2_gt4[pp]+-(0.5*grad_1_gt5[pp])))+-(grad2_2_2_chi[pp]))+_1849;
double _1903 = _1850+_1852+_1857-_1860-_1865-_1869+_1872+_1876+_1880+_1884-_1887-_1891+_1893-_1896+_1898-_1900-_1902;
double _1914 = _1903-3.0*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_425+DENDRO_367*DENDRO_419+_1909+2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1911;
double _1924 = alpha[pp]*(_1914-DENDRO_420*(0.5*grad_0_gt5[pp]+1.0*grad_2_gt2[pp])-_1919-DENDRO_426*3.0*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])+_532+4*grad_2_Gt2[pp]*gt5[pp]);
double _1753 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_235*0.5*grad_1_gt0[pp]+0.25*grad_2_gt0[pp]*DENDRO_242);
double DENDRO_489 = (-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _1466 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_259*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+DENDRO_469);
double _1467 = _1442+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1443-_1448-_1451-_1458-_1462-_1466;
double _1473 = _1467+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1468+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1471;
double _1511 = _1473+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1474+_1480+_1484+_1488+_1490-_1493+_1495-_1497-_1499+_1502+_1506+_1510;
double _1423 = _1511+DENDRO_419*DENDRO_451-DENDRO_259*4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(1.0*grad_2_gt4[pp]-0.25*grad_1_gt5[pp])-_1518+_616+4*grad_1_Gt1[pp]*gt3[pp];
double _1522 = (gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(_1418+DENDRO_447*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+_1420)+alpha[pp]*_1423);
double _1344 = DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*((grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_259+DENDRO_584+DENDRO_656);
double _1807 = DENDRO_160*4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(1.0*grad_2_gt2[pp]-0.25*grad_0_gt5[pp]);
double DENDRO_279 = (gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(3*(1.0/chi[pp]))*pow(grad_1_chi[pp],2)+2*grad2_1_1_chi[pp]);
double _322 = (-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(2*grad2_2_2_chi[pp]-pow(grad_2_chi[pp],2)*3*(1.0/chi[pp]))+DENDRO_279+DENDRO_281;
double DENDRO_837 = (_322-2*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(2*grad2_1_2_chi[pp]-3*(1.0/chi[pp])*grad_2_chi[pp]*grad_1_chi[pp])+DENDRO_288-DENDRO_292+DENDRO_294+DENDRO_296+DENDRO_298)*(1.0/chi[pp]);
double DENDRO_299 = _322-2*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(2*grad2_1_2_chi[pp]-3*(1.0/chi[pp])*grad_2_chi[pp]*grad_1_chi[pp])+DENDRO_288-DENDRO_292+DENDRO_294+DENDRO_296+DENDRO_298;
double _2761 = DENDRO_96*(-(DENDRO_89)+_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double _1723 = -(grad2_0_0_chi[pp])+DENDRO_117*(_623+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]));
double DENDRO_664 = (_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]);
double _1757 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_242*0.5*grad_2_gt0[pp]+0.25*grad_1_gt0[pp]*DENDRO_235);
double DENDRO_685 = (_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _2026 = At1[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*2*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double _2027 = _2023-2*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At4[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-_2026;
double _2197 = (_2027+At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*2*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*DENDRO_120*9*(1.0/chi[pp]);
double DENDRO_949 = (_2027+At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*2*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*2*alpha[pp]*DENDRO_945;
double _2173 = (_2027+At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*2*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*grad_0_alpha[pp]*DENDRO_917;
double _3927 = At0[pp]*(_2027+At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*2*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*DENDRO_888;
double _2198 = DENDRO_41*(1.0/3.0)*alpha[pp]*(-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*4*grad_0_K[pp]+_2197);
double _2145 = _2139+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*DENDRO_949-_2142-_2144;
double _1237 = 0.5*grad_2_gt3[pp]*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _1239 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_1237+0.5*grad_2_gt5[pp]*(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+DENDRO_259*0.25*grad_1_gt5[pp]);
double _1252 = _1226+_1228+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_1230+_1239-_1242-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1244-_1248+_1251;
double _1267 = _1252+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1253+_1257+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1259+_1263+_1266;
double _1284 = _1267-_1269-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1271-_1275-_1278-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_1280;
double _1301 = alpha[pp]*(_1284-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1285-_1290+_1294+_864+DENDRO_419*DENDRO_701+DENDRO_702-_1298+DENDRO_731);
double _2663 = -(grad2_0_0_chi[pp])+(_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*DENDRO_117;
double _2665 = _2663+DENDRO_118*(-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _2677 = -(2*(1.0/chi[pp]))*(_2665+(_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*DENDRO_120)-_2669-_2672-_2676;
double _2681 = _2677-4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_173*0.25*grad_2_gt0[pp]+DENDRO_179*0.5*grad_1_gt0[pp]);
double _2692 = _2681-4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(DENDRO_173*0.5*grad_2_gt0[pp]+DENDRO_179*0.25*grad_1_gt0[pp])+_2686+_2691;
double _2707 = _2692-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_2693+_2701-4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_2703;
double _2734 = _2707-4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_2708+_2716+_2719-_2722-_2725-_2727-_2729-_2731-DENDRO_299*DENDRO_104*(1.0/chi[pp]);
double _2743 = _2734+_2736+DENDRO_141*4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(1.0*grad_1_gt1[pp])+0.25*grad_0_gt3[pp])+_2742;
double _2747 = _2743+DENDRO_316*3.0*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])+DENDRO_173*grad_1_gt0[pp]*3.0*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2));
double _2752 = _2657+3*alpha[pp]*(_2747+_411-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344*grad2_0_1_gt0[pp])-12*grad2_0_0_alpha[pp];
double _2966 = 1.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_2965+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_824);
double _1725 = _1723+DENDRO_118*((-(gt0[pp]*gt5[pp])+pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79);
double _1750 = -(2*(1.0/chi[pp]))*(_1725+DENDRO_120*(_1727+(-(gt1[pp]*gt4[pp])+gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])))-_1731+_1738+_1745+_1749;
double _1796 = _1750+_1753+_1757-_1759-_1764-_1770+_1775+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1777-_1782-_1786+_1789+_1793+_1795;
double _1811 = _1796+_1797-_1800+_1802+DENDRO_104*(1.0/chi[pp])*DENDRO_418-_1807-DENDRO_141*4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_1_gt1[pp]-0.25*grad_0_gt3[pp]);
double _1815 = _1811-3.0*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*DENDRO_235*grad_2_gt0[pp]-3.0*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_525;
double _1818 = (-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*(_1717+alpha[pp]*(_1815+_411-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344*grad2_0_1_gt0[pp]));
double _1956 = -(grad_2_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+grad_1_gt5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _2897 = (grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_179+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_179+grad_2_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3118 = grad_1_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_173+(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_173;
double _3091 = _1964+0.25*grad_0_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+DENDRO_173*0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]);
double _2855 = DENDRO_179*0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+_1964+0.25*grad_0_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3307 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(grad_0_gt3[pp])*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+DENDRO_489);
double _3158 = (-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_179+grad_2_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+grad_0_gt5[pp]*DENDRO_173;
double _3716 = 2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(grad_0_gt5[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+_3715);
double _3677 = 0.25*grad_2_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+DENDRO_179*0.5*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp]);
double _3453 = -(0.5*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))-_3452;
double _3439 = (grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_179+grad_2_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+grad_0_gt5[pp]*DENDRO_173;
double _2958 = _2959+grad_1_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_173;
double _3097 = _3090-4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3091-4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3094;
double _3117 = _3097-4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3098+_3105+4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_3107-_3116;
double _3140 = _3117+DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_3118+_3125+_3132+4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_3134+2.0*grad_2_Gt0[pp]*gt0[pp];
double _3252 = DENDRO_173*0.5*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])+0.25*grad_1_gt0[pp]*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3260 = _3240-4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3241+_3245+_3247+_3250+4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_3252+_3259;
double _3492 = -(0.5*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+DENDRO_685;
double _3425 = _3426+0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+DENDRO_722;
double _3145 = _3140+2.0*grad_2_Gt1[pp]*gt1[pp]+grad_0_Gt0[pp]*2.0*gt2[pp]+grad_2_Gt2[pp]*2.0*gt2[pp]+grad_0_Gt1[pp]*2.0*gt4[pp]+grad_0_Gt2[pp]*2.0*gt5[pp];
double _3430 = _3414+_3419+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3421-4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3425;
double _3434 = -(0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]))*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3513 = 0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3683 = (1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp])*1.0*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3262 = (1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])*1.0*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3509 = -(0.25*grad_0_gt5[pp])*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+_3513-_3515;
double _3532 = -(0.25*grad_0_gt3[pp])*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+DENDRO_664+DENDRO_719;
double _3686 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_3683+0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_179);
double _3264 = _3260+4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_3262+0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_173);
double _3270 = _3264+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3265+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3268;
double _3316 = _3270+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3271+_3277+_3283+_3289-_3291+_3296-_3298+_3301+_3307+_3313-DENDRO_451*DENDRO_837;
double _3325 = 3*alpha[pp]*(_3316+_3318+DENDRO_495*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))+_3323+_616+4*grad_1_Gt1[pp]*gt3[pp]);
double _3454 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_3453+DENDRO_173*0.5*(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp])));
double _3147 = _3145+-(pow(chi[pp],-(2))*grad_2_chi[pp])*grad_0_chi[pp]+4*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_gt2[pp];
double _3149 = _3147+2.0*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_gt2[pp]+2.0*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_gt2[pp];
double _3151 = _3149+2.0*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_gt2[pp]+-((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344)*grad2_1_2_gt2[pp];
double _3436 = 0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3438 = _3430+4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_3434+0.5*grad_0_gt5[pp]*(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))-_3436);
double _3460 = _3438+DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_3439-_3446-_3454+4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_3456;
double _3478 = _3460+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3461+_3471+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3473;
double _3493 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_3492-0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_173);
double _3561 = _3557*(-(DENDRO_715)+_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3505 = _3478+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3479+_3493+_3500-_3502-_3504;
double _3523 = _3505+4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_3506-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_3509-_3522;
double _3542 = _3523-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_3524-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_3532+_3541;
double _3399 = _3542+4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_3543+2.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3547-DENDRO_701*DENDRO_837+DENDRO_702+DENDRO_731;
double _1341 = 4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_571+0.25*(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_235);
double _1354 = _1333-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_1334-_1338-_1341+_1344+_1347+_1350+_1353;
double _1360 = _1354+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1355+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1358;
double _1370 = _1360-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1361-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1366;
double _1383 = _1370-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1371-_1377-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_1379;
double _1402 = _1383-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_1384-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1389-_1394+_1398-_1401;
double _1405 = (gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_886-4*grad2_1_2_alpha[pp]+alpha[pp]*(_1402+_864+DENDRO_419*DENDRO_701+DENDRO_702));
double _2836 = 4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_2835+DENDRO_141*0.25*grad_2_gt3[pp]+grad_1_gt3[pp]*DENDRO_545);
double _2854 = -(2*(1.0/chi[pp]))*(_2824+_2825+-(grad2_0_1_chi[pp]))-_2828-_2830+_2836-_2842+_2850+4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_2852;
double _2865 = _2854+4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_2855+4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_2858+_2864;
double _2040 = _2038+2*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At2[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double _3629 = DENDRO_120*(_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _2042 = _2040-At1[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*2*(gt0[pp]*gt5[pp]-pow(gt2[pp],2));
double DENDRO_961 = (_2042-At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*2*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*DENDRO_118*9*(1.0/chi[pp]);
double DENDRO_947 = (_2042-At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*2*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*2*alpha[pp]*DENDRO_945;
double _2274 = (_2042-At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*2*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*grad_1_alpha[pp]*DENDRO_917;
double _3938 = At3[pp]*DENDRO_888*(_2042-At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*2*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)));
double _4048 = DENDRO_41*(1.0/3.0)*alpha[pp]*(-((gt0[pp]*gt5[pp]-pow(gt2[pp],2))*4*grad_1_K[pp])+DENDRO_961);
double _3978 = (_1959+(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*DENDRO_947;
double _4033 = (_220-(0.5*grad_0_gt3[pp]+-(1.0*grad_1_gt1[pp]))*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*DENDRO_947;
double _761 = DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*((-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*DENDRO_208+_757+DENDRO_650);
double _770 = _761-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(DENDRO_546+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_545)-4*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*_765;
double _779 = _770+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_771+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_774-_778;
double _1596 = alpha[pp]*(_1594+_779-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_780-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_783);
double DENDRO_655 = _779-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_780-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_783;
double _3035 = _3151+-((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344)*grad2_0_1_gt2[pp]-DENDRO_615*DENDRO_837+DENDRO_655-_3156-1.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3158;
double DENDRO_901 = (-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At0[pp]*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double DENDRO_911 = (-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*At2[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _2230 = _2227+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+DENDRO_911+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2));
double _2090 = _2087-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-DENDRO_911-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At3[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2));
double _2269 = (-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])+DENDRO_79)*DENDRO_949;
double _1156 = (_153-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(1.0*grad_0_gt2[pp]-0.5*grad_2_gt0[pp]))*0.5*grad_0_gt3[pp];
double _1155 = _1156+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*0.5*grad_0_gt0[pp]+0.25*grad_1_gt0[pp]*DENDRO_242;
double _1159 = _1147+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1148+_1153+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1155;
double _1180 = _1159+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1160+_1167-4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_1169-_1174-_1179;
double _1193 = _1180-_1184-4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_1186-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1191;
double _1207 = _1193-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1194-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1199-_1206;
double _1114 = DENDRO_798+alpha[pp]*(_1207-1.0*DENDRO_144*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_1209+(grad_0_gt4[pp]+grad_2_gt1[pp]-grad_1_gt2[pp])*DENDRO_160)-_1214+DENDRO_782);
double _1219 = -(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(DENDRO_798+alpha[pp]*_1032)-(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_1114;
double _1598 = _1219-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(_886-4*grad2_1_2_alpha[pp]+_1301)-_1405+_1522+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_748+2.0*grad_1_alpha[pp]*DENDRO_634+_1596);
double DENDRO_801 = DENDRO_41*(_1598+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_748+2.0*grad_1_alpha[pp]*DENDRO_634+_1700)+_1818+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(_1829+_1835+_1924));
double _3161 = DENDRO_615*(_1598+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_748+2.0*grad_1_alpha[pp]*DENDRO_634+_1700)+_1818+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(_1829+_1835+_1924));
double _2972 = DENDRO_781*(_1598+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_748+2.0*grad_1_alpha[pp]*DENDRO_634+_1700)+_1818+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(_1829+_1835+_1924));
double _3563 = DENDRO_701*(_1598+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*(_748+2.0*grad_1_alpha[pp]*DENDRO_634+_1700)+_1818+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(_1829+_1835+_1924));
double _3595 = _3397+(1.0/12.0)*chi[pp]*(3*alpha[pp]*_3399-_3561+_3563-12*grad2_1_2_alpha[pp]-_3573+_3579)-alpha[pp]*(-(At4[pp])*K[pp]+_3588+_3590+_3592);
double _3196 = _3033+(1.0/12.0)*chi[pp]*(3*alpha[pp]*_3035+_3161-12*grad2_0_2_alpha[pp]+DENDRO_634*6.0*grad_1_alpha[pp]-_3171-_3176)-alpha[pp]*_3179+beta0[pp]*agrad_0_At2[pp];
double _3368 = _3217+(1.0/12.0)*chi[pp]*(_3325-12*grad2_1_1_alpha[pp]+DENDRO_801*gt3[pp]+_3336+_3344+_3352)-alpha[pp]*(-(At3[pp])*K[pp]+_3361+_3363+_3365);
At_rhs00[pp] = _2650+(1.0/12.0)*chi[pp]*(_2752+DENDRO_801*gt0[pp]+_2756-_2761)-_2786+beta0[pp]*agrad_0_At0[pp]+beta1[pp]*agrad_1_At0[pp]+beta2[pp]*agrad_2_At0[pp];
double _2220 = _2218+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+DENDRO_901;
double _2072 = (-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*At2[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _3675 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_3674+0.25*grad_0_gt5[pp]*DENDRO_824);
double _3626 = DENDRO_117*(_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3631 = -(2*(1.0/chi[pp]))*(_3626+DENDRO_118*(_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))+_3629+-(grad2_2_2_chi[pp]));
double _3655 = _3631+_3632*3.0*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])-_3635-_3640-4*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_3642+_3649+_3654;
double _3676 = _3655+_3661-_3667-4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(0.5*grad_0_gt5[pp]*DENDRO_824+DENDRO_848)-_3675;
double _3698 = _3676-4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_3677-_3686+_3689+_3694+grad_1_gt5[pp]*3.0*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_824;
double _3720 = _3698-_3699+_3704-_3706-DENDRO_367*DENDRO_837-_3711-_3716+DENDRO_420*(-(1.0*grad_2_gt2[pp])+-(0.5*grad_0_gt5[pp]));
double _3725 = _3720+DENDRO_421*(-(1.0*grad_2_gt4[pp])+-(0.5*grad_1_gt5[pp]))+DENDRO_423*(_1956+(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]));
double _3745 = _3620+(1.0/12.0)*chi[pp]*(3*alpha[pp]*(_3725+_532+4*grad_2_Gt2[pp]*gt5[pp])-12*grad2_2_2_alpha[pp]+_3730+DENDRO_801*gt5[pp]-_3738-_3743);
At_rhs22[pp] = _3745-alpha[pp]*(_3749-At5[pp]*K[pp]+_3754+_3756)+beta0[pp]*agrad_0_At5[pp]+beta1[pp]*agrad_1_At5[pp]+beta2[pp]*agrad_2_At5[pp];
double _2054 = _2052-2*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*At1[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]);
double _2056 = _2054+At2[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*2*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _2370 = (_2056-At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*2*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*DENDRO_117*9*(1.0/chi[pp]);
double DENDRO_948 = (_2056-At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*2*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*2*alpha[pp]*DENDRO_945;
double _3946 = At5[pp]*DENDRO_888*(_2056-At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*2*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]));
double _2344 = (_2056-At4[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*2*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]))*DENDRO_917*grad_2_alpha[pp];
double _2371 = DENDRO_41*(1.0/3.0)*alpha[pp]*(-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*4*grad_2_K[pp]+_2370);
double _2152 = (_238+(0.5*grad_0_gt5[pp]+-(1.0*grad_2_gt2[pp]))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]))*DENDRO_948;
double _2311 = (_255-(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_948;
double _3985 = DENDRO_41*(1.0/3.0)*alpha[pp]*((gt0[pp]*gt5[pp]-pow(gt2[pp],2))*4*grad_1_K[pp]-DENDRO_961);
double _2221 = _2220-At1[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double _2304 = (_152+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp]))*DENDRO_949-_2303;
double _2881 = DENDRO_164*0.25*(grad_0_gt4[pp]-grad_2_gt1[pp]+grad_1_gt2[pp])+DENDRO_164*0.25*(-(grad_0_gt4[pp])+grad_2_gt1[pp]+grad_1_gt2[pp]);
double _2883 = 4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(_2881+0.5*grad_2_gt3[pp]*(_135-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(1.0*grad_0_gt1[pp]-0.5*grad_1_gt0[pp])));
double _2899 = _2865+4*DENDRO_144*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*_2866+_2876-_2883+_2891-_2895+DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_2897;
double _2924 = _2899+4*DENDRO_144*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*_2900+4*DENDRO_144*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*_2906-_2916-_2923;
double _2941 = _2924+4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_2925+4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_2931-_2940;
double _2947 = _2941+2.0*grad_1_Gt0[pp]*gt0[pp]+grad_0_Gt0[pp]*2.0*gt1[pp]+grad_1_Gt1[pp]*2.0*gt1[pp]+grad_1_Gt2[pp]*2.0*gt2[pp]+2.0*grad_0_Gt1[pp]*gt3[pp]+grad_0_Gt2[pp]*2.0*gt4[pp];
double _2949 = _2947+-(grad_1_chi[pp]*grad_0_chi[pp])*pow(chi[pp],-(2))+4*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_gt1[pp];
double _2951 = _2949+2.0*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_gt1[pp]+2.0*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_gt1[pp];
double _2953 = _2951+2.0*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_gt1[pp]+-((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_344)*grad2_1_2_gt1[pp];
double _2957 = _2953+-((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_344)*grad2_0_1_gt1[pp]-DENDRO_781*DENDRO_837+_1003-4*DENDRO_144*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*_1004;
double _2817 = 3*alpha[pp]*(_2957-1.0*DENDRO_144*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*_2958-_2966)-_2968*_2969+_2972-12*grad2_0_1_alpha[pp]-_2982+_2988;
At_rhs01[pp] = _2816+(1.0/12.0)*chi[pp]*_2817-_3003+beta0[pp]*agrad_0_At1[pp]+beta1[pp]*agrad_1_At1[pp]+beta2[pp]*agrad_2_At1[pp];
double _2091 = _2090+(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At4[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _2356 = (_2091-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_118*9*(1.0/chi[pp]);
double _4055 = (_2091-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_117*9*(1.0/chi[pp]);
double _2309 = DENDRO_259*(_2091-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*_2125*alpha[pp];
double _4031 = DENDRO_208*(_2091-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*_2125*alpha[pp];
double _2148 = DENDRO_240*(_2091-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*_2125*alpha[pp];
double _2348 = (_2091-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*grad_1_alpha[pp]*DENDRO_917;
double _4043 = (_2091-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_917*grad_2_alpha[pp];
double _3942 = At4[pp]*DENDRO_895*(_2091-(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _4056 = DENDRO_41*(1.0/3.0)*alpha[pp]*(_4055+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*4*grad_2_K[pp]);
double _2222 = _2221+At2[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double _3990 = -(DENDRO_120*9*(1.0/chi[pp]))*(_2222+(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3968 = _3964*_2125*alpha[pp]*(_2222+(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))+_3967;
double _3980 = grad_0_alpha[pp]*DENDRO_917*(_2222+(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _3991 = DENDRO_41*(1.0/3.0)*alpha[pp]*(_3990+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*4*grad_0_K[pp]);
double _2081 = _2079-(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-DENDRO_901;
double _2082 = _2081+At1[pp]*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double _2083 = _2082-At2[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double _4034 = -(B1[pp])*eta+DENDRO_141*(_2083-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*_2125*alpha[pp]-_4029+_4031+_4033;
double _2306 = _2304+DENDRO_155*(_2083-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*_2125*alpha[pp];
double _2183 = (_2083-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*DENDRO_118*9*(1.0/chi[pp]);
double _4051 = (_2083-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*DENDRO_120*9*(1.0/chi[pp]);
double _3931 = _3927+At1[pp]*DENDRO_895*(_2083-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]));
double _2150 = DENDRO_242*(_2083-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*_2125*alpha[pp];
double _2177 = (_2083-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*grad_1_alpha[pp]*DENDRO_917;
double _4041 = (_2083-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*grad_0_alpha[pp]*DENDRO_917;
double _4052 = DENDRO_41*(1.0/3.0)*alpha[pp]*(_4051+(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*4*grad_0_K[pp]);
double _4057 = _4034-_4035+DENDRO_411*(2.0/3.0)*(grad_0_beta0[pp]+grad_1_beta1[pp]+grad_2_beta2[pp])-_4039-_4041-_4043-_4048-_4052-_4056;
double _2184 = DENDRO_41*(1.0/3.0)*alpha[pp]*((-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*4*grad_1_K[pp]+_2183);
double _2231 = _2230-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*At4[pp]*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp]);
double _3973 = _3968-DENDRO_208*_2125*alpha[pp]*(_2231+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))+_3972;
double _3996 = -(DENDRO_117*9*(1.0/chi[pp]))*(_2231+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3982 = DENDRO_917*grad_2_alpha[pp]*(_2231+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]));
double _3992 = _3973-DENDRO_233*(2.0/3.0)*(grad_0_beta0[pp]+grad_1_beta1[pp]+grad_2_beta2[pp])+_3976-_3978+_3980+_3982+_3985-_3991;
double _3998 = _3992-DENDRO_41*(1.0/3.0)*alpha[pp]*(_3996+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*4*grad_2_K[pp]);
double _2255 = (_213+(0.5*grad_1_gt5[pp]+-(1.0*grad_2_gt4[pp]))*(gt0[pp]*gt5[pp]-pow(gt2[pp],2)))*DENDRO_948;
double _2067 = _2065-At1[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]);
double _2069 = _2067+(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*At0[pp]*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp]);
double _2073 = _2069-(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*At1[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])+_2072;
double _3924 = _3931+At2[pp]*(_2073+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*DENDRO_895+_3938+_3942+_3946+pow(K[pp],2);
K_rhs[pp] = _3867-2*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41*chi[pp]*_3869+_3922+(1.0/3.0)*alpha[pp]*_3924+beta0[pp]*agrad_0_K[pp]+beta1[pp]*agrad_1_K[pp]+beta2[pp]*agrad_2_K[pp];
double _2153 = _2145+DENDRO_235*(_2073+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*_2125*alpha[pp]+_2148+_2150+_2152;
double _2158 = _2153+(_241+0.5*grad_1_gt3[pp]*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp]))*DENDRO_947-DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_beta0[pp];
double _2169 = _2158-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_beta0[pp]-(4.0/3.0)*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_beta0[pp]+_2168;
double _2312 = _2306+DENDRO_160*(_2073+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*_2125*alpha[pp]+_2309+_2311;
double _2320 = _2312+(_258+0.5*grad_1_gt3[pp]*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp]))*DENDRO_947-(1.0/3.0)*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_1_2_beta1[pp];
double _2328 = _2320-(4.0/3.0)*DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_beta2[pp]-DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*(1.0/3.0)*grad2_0_2_beta0[pp];
double _2340 = _2328-(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_beta2[pp]-DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_beta2[pp]-_2335+_2337-_2339;
double _2363 = (_2073+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*DENDRO_120*9*(1.0/chi[pp]);
double _2256 = DENDRO_134*(_2073+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*_2125*alpha[pp]+_2255;
double _2265 = _2256-DENDRO_41*(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*grad2_2_2_beta1[pp]-(4.0/3.0)*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41*grad2_1_1_beta1[pp];
double _2275 = _2265-DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_0_beta1[pp]+_2269-2*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*DENDRO_41*grad2_0_2_beta1[pp]-_2274;
double _2279 = _2275+grad2_0_0_beta0[pp]*(1.0/3.0)*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41+grad2_0_1_beta1[pp]*(7.0/3.0)*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41;
double _2283 = _2279+grad2_0_2_beta2[pp]*(1.0/3.0)*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41+grad2_1_2_beta1[pp]*(7.0/3.0)*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41;
double _2287 = _2283+grad2_2_2_beta2[pp]*(1.0/3.0)*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41-grad2_1_2_beta2[pp]*(1.0/3.0)*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41;
double _2291 = _2287-grad2_0_1_beta0[pp]*(1.0/3.0)*(gt0[pp]*gt5[pp]-pow(gt2[pp],2))*DENDRO_41+(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41*(1.0/3.0)*grad2_0_2_beta0[pp];
double _2190 = (_2073+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*DENDRO_117*9*(1.0/chi[pp]);
double _4061 = _4057-(beta0[pp]*agrad_0_Gt1[pp]+beta1[pp]*agrad_1_Gt1[pp]+beta2[pp]*agrad_2_Gt1[pp])*lambda[3]+_2291+beta0[pp]*agrad_0_Gt1[pp]+beta1[pp]*agrad_1_Gt1[pp]+beta2[pp]*agrad_2_Gt1[pp];
double _2175 = (_2073+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*DENDRO_917*grad_2_alpha[pp];
double _2185 = _2169+2*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41*grad2_1_2_beta0[pp]-_2173-_2175-_2177-grad2_0_2_beta0[pp]*(7.0/3.0)*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-_2184;
double _2199 = _2185-DENDRO_41*(1.0/3.0)*alpha[pp]*(-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*4*grad_2_K[pp]+_2190)-_2198;
double _2203 = _2199-grad2_0_1_beta1[pp]*(1.0/3.0)*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])-(1.0/3.0)*DENDRO_41*(-(pow(gt4[pp],2))+gt3[pp]*gt5[pp])*grad2_0_2_beta2[pp];
double _2207 = _2203-grad2_1_2_beta1[pp]*(1.0/3.0)*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-(1.0/3.0)*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_2_2_beta2[pp];
double _2211 = _2207+grad2_1_1_beta1[pp]*(1.0/3.0)*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41+(1.0/3.0)*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41*grad2_1_2_beta2[pp];
Gt_rhs0[pp] = 1.0*(_2211+grad2_0_1_beta0[pp]*(7.0/3.0)*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41+beta0[pp]*agrad_0_Gt0[pp]+beta1[pp]*agrad_1_Gt0[pp]+beta2[pp]*agrad_2_Gt0[pp]);
double DENDRO_955 = _2211+grad2_0_1_beta0[pp]*(7.0/3.0)*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41+beta0[pp]*agrad_0_Gt0[pp]+beta1[pp]*agrad_1_Gt0[pp]+beta2[pp]*agrad_2_Gt0[pp];
B_rhs0[pp] = -(B0[pp])*eta-(beta0[pp]*agrad_0_Gt0[pp]+beta1[pp]*agrad_1_Gt0[pp]+beta2[pp]*agrad_2_Gt0[pp])*lambda[3]+DENDRO_955+lambda[2]*(beta0[pp]*agrad_0_B0[pp]+beta1[pp]*agrad_1_B0[pp]+beta2[pp]*agrad_2_B0[pp]);
double _2346 = (_2073+(-(pow(gt1[pp],2))+gt0[pp]*gt3[pp])*At5[pp]*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp]))*grad_0_alpha[pp]*DENDRO_917;
double _2351 = _2340+2*(-(gt2[pp]*gt4[pp])+gt1[pp]*gt5[pp])*DENDRO_41*grad2_0_1_beta2[pp]-_2344-_2346-_2348-(7.0/3.0)*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*grad2_0_2_beta2[pp];
double _2358 = _2351-grad2_0_0_beta0[pp]*(1.0/3.0)*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])-DENDRO_41*(1.0/3.0)*alpha[pp]*((gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*4*grad_1_K[pp]+_2356);
double _2372 = _2358-DENDRO_41*(1.0/3.0)*alpha[pp]*(-(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])*4*grad_0_K[pp]+_2363)-_2371;
double _2376 = _2372-grad2_0_1_beta1[pp]*(1.0/3.0)*DENDRO_41*(gt1[pp]*gt4[pp]-gt2[pp]*gt3[pp])+grad2_1_1_beta1[pp]*(1.0/3.0)*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41;
double _2380 = _2376+grad2_1_2_beta2[pp]*(7.0/3.0)*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41+grad2_0_1_beta0[pp]*(1.0/3.0)*(gt0[pp]*gt4[pp]-gt1[pp]*gt2[pp])*DENDRO_41;
At_rhs11[pp] = _3368+beta0[pp]*agrad_0_At3[pp]+beta1[pp]*agrad_1_At3[pp]+beta2[pp]*agrad_2_At3[pp];
Gt_rhs1[pp] = _3998+_2291+beta0[pp]*agrad_0_Gt1[pp]+beta1[pp]*agrad_1_Gt1[pp]+beta2[pp]*agrad_2_Gt1[pp];
B_rhs2[pp] = -(B2[pp])*eta-(beta0[pp]*agrad_0_Gt2[pp]+beta1[pp]*agrad_1_Gt2[pp]+beta2[pp]*agrad_2_Gt2[pp])*lambda[3]+_2380+beta0[pp]*agrad_0_Gt2[pp]+beta1[pp]*agrad_1_Gt2[pp]+beta2[pp]*agrad_2_Gt2[pp]+lambda[2]*(beta0[pp]*agrad_0_B2[pp]+beta1[pp]*agrad_1_B2[pp]+beta2[pp]*agrad_2_B2[pp]);
gt_rhs11[pp] = _2544+(4.0/3.0)*grad_1_beta1[pp]*gt3[pp]+beta0[pp]*agrad_0_gt3[pp]+beta1[pp]*agrad_1_gt3[pp]+beta2[pp]*agrad_2_gt3[pp];
a_rhs[pp] = -(2*alpha[pp])*K[pp]+lambda[0]*(beta0[pp]*agrad_0_alpha[pp]+beta1[pp]*agrad_1_alpha[pp]+beta2[pp]*agrad_2_alpha[pp]);
gt_rhs01[pp] = _2483+grad_0_beta1[pp]*gt3[pp]+grad_0_beta2[pp]*gt4[pp]+beta0[pp]*agrad_0_gt1[pp]+beta1[pp]*agrad_1_gt1[pp]+beta2[pp]*agrad_2_gt1[pp];
At_rhs02[pp] = _3196+beta1[pp]*agrad_1_At2[pp]+beta2[pp]*agrad_2_At2[pp];
gt_rhs12[pp] = _2576+(1.0/3.0)*gt4[pp]*grad_2_beta2[pp]+beta0[pp]*agrad_0_gt4[pp]+beta1[pp]*agrad_1_gt4[pp]+beta2[pp]*agrad_2_gt4[pp];
At_rhs12[pp] = _3595+beta0[pp]*agrad_0_At4[pp]+beta1[pp]*agrad_1_At4[pp]+beta2[pp]*agrad_2_At4[pp];
B_rhs1[pp] = _4061+lambda[2]*(beta0[pp]*agrad_0_B1[pp]+beta1[pp]*agrad_1_B1[pp]+beta2[pp]*agrad_2_B1[pp]);
Gt_rhs2[pp] = 1.0*(_2380+beta0[pp]*agrad_0_Gt2[pp]+beta1[pp]*agrad_1_Gt2[pp]+beta2[pp]*agrad_2_Gt2[pp]);
gt_rhs22[pp] = _2603+(4.0/3.0)*grad_2_beta2[pp]*gt5[pp]+beta0[pp]*agrad_0_gt5[pp]+beta1[pp]*agrad_1_gt5[pp]+beta2[pp]*agrad_2_gt5[pp];
gt_rhs02[pp] = _2514+grad_0_beta1[pp]*gt4[pp]+grad_0_beta2[pp]*gt5[pp]+beta0[pp]*agrad_0_gt2[pp]+beta1[pp]*agrad_1_gt2[pp]+beta2[pp]*agrad_2_gt2[pp];
gt_rhs00[pp] = _2455+grad_0_beta1[pp]*2*gt1[pp]+beta0[pp]*agrad_0_gt0[pp]+beta1[pp]*agrad_1_gt0[pp]+beta2[pp]*agrad_2_gt0[pp];

//[[[end]]]
#endif



                bssn::timer::t_rhs.stop();

                /* debugging */
                unsigned int qi = 46 - 1;
                unsigned int qj = 10 - 1;
                unsigned int qk = 60 - 1;
                unsigned int qidx = qi + nx*(qj + ny*qk);
                if (0 && qidx == pp) {
                    std::cout << ".... end OPTIMIZED debug stuff..." << std::endl;
                }

            }
        }
    }


    if (bflag != 0) {

        bssn::timer::t_bdyc.start();


#ifdef BSSN_KERR_SCHILD_TEST

        freeze_bcs(a_rhs, sz, bflag);
        freeze_bcs(chi_rhs, sz, bflag);
        freeze_bcs(K_rhs, sz, bflag);
        freeze_bcs(b_rhs0, sz, bflag);
        freeze_bcs(b_rhs1, sz, bflag);
        freeze_bcs(b_rhs2, sz, bflag);
        freeze_bcs(Gt_rhs0, sz, bflag);
        freeze_bcs(Gt_rhs1, sz, bflag);
        freeze_bcs(Gt_rhs2, sz, bflag);
        freeze_bcs(B_rhs0, sz, bflag);
        freeze_bcs(B_rhs1, sz, bflag);
        freeze_bcs(B_rhs2, sz, bflag);
        freeze_bcs(At_rhs00, sz, bflag);
        freeze_bcs(At_rhs01, sz, bflag);
        freeze_bcs(At_rhs02, sz, bflag);
        freeze_bcs(At_rhs11, sz, bflag);
        freeze_bcs(At_rhs12, sz, bflag);
        freeze_bcs(At_rhs22, sz, bflag);
        freeze_bcs(gt_rhs00, sz, bflag);
        freeze_bcs(gt_rhs01, sz, bflag);
        freeze_bcs(gt_rhs02, sz, bflag);
        freeze_bcs(gt_rhs11, sz, bflag);
        freeze_bcs(gt_rhs12, sz, bflag);
        freeze_bcs(gt_rhs22, sz, bflag);

#else

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

#endif



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
                pp = i + nx*(j + ny*k);

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

#if 0
    for (unsigned int m = 0; m < 24; m++) {
        std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
    }
#endif



}



void bssnrhs_sep(double **unzipVarsRHS, const double **uZipVars,
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
    const double eta_power[2] = {BSSN_ETA_POWER[0], BSSN_ETA_POWER[1]};



    int idx[3];

    unsigned int n = sz[0]*sz[1]*sz[2];

    double * CalGt0 =new double[n];
    double * CalGt1 =new double[n];
    double * CalGt2 =new double[n];

    double *Gt_rhs_s1_0 =new double [n];
    double *Gt_rhs_s1_1 =new double [n];
    double *Gt_rhs_s1_2 =new double [n];

    double *Gt_rhs_s2_0 =new double [n];
    double *Gt_rhs_s2_1 =new double [n];
    double *Gt_rhs_s2_2 =new double [n];

    double *Gt_rhs_s3_0 =new double [n];
    double *Gt_rhs_s3_1 =new double [n];
    double *Gt_rhs_s3_2 =new double [n];

    double *Gt_rhs_s4_0 =new double [n];
    double *Gt_rhs_s4_1 =new double [n];
    double *Gt_rhs_s4_2 =new double [n];

    double *Gt_rhs_s5_0 =new double [n];
    double *Gt_rhs_s5_1 =new double [n];
    double *Gt_rhs_s5_2 =new double [n];

    double *Gt_rhs_s6_0 =new double [n];
    double *Gt_rhs_s6_1 =new double [n];
    double *Gt_rhs_s6_2 =new double [n];

    double *Gt_rhs_s7_0 =new double [n];
    double *Gt_rhs_s7_1 =new double [n];
    double *Gt_rhs_s7_2 =new double [n];



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

    register double x;
    register double y;
    register double z;
    register unsigned int pp;

    double r_coord;
    double eta;

    //cout << "begin loop" << endl;
    bssn::timer::t_rhs.start();

    bssn::timer::t_rhs_a.start();

    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }

                // a_rhs

                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.a_rhs]
                vnames = ['a_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 12 
                // Dendro: printing temp variables

                // Dendro: printing variables
                //--
                a_rhs[pp] = -2*K[pp]*alpha[pp] + lambda[0]*(beta0[pp]*agrad_0_alpha[pp] + beta1[pp]*agrad_1_alpha[pp] + beta2[pp]*agrad_2_alpha[pp]);
                // Dendro: reduced ops: 12
                // Dendro: }}} 
                //[[[end]]]

            }
        }
    }


    bssn::timer::t_rhs_a.stop();

    bssn::timer::t_rhs_b.start();
    // b_rhs


    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }

                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.b_rhs]
                vnames = ['b_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 51 
                // Dendro: printing temp variables
                double DENDRO_0 = (3.0/4.0)*alpha[pp]*lambda_f[1] + (3.0/4.0)*lambda_f[0];

                // Dendro: printing variables
                //--
                b_rhs0[pp] = B0[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta0[pp] + beta1[pp]*agrad_1_beta0[pp] + beta2[pp]*agrad_2_beta0[pp]);
                //--
                b_rhs1[pp] = B1[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta1[pp] + beta1[pp]*agrad_1_beta1[pp] + beta2[pp]*agrad_2_beta1[pp]);
                //--
                b_rhs2[pp] = B2[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta2[pp] + beta1[pp]*agrad_1_beta2[pp] + beta2[pp]*agrad_2_beta2[pp]);
                // Dendro: reduced ops: 39
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }
    bssn::timer::t_rhs_b.stop();

    bssn::timer::t_rhs_gt.start();

//gt_rhs

    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.gt_rhs]
                vnames = ['gt_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 210 
                // Dendro: printing temp variables
                double DENDRO_0 = 2*alpha[pp];
                double DENDRO_1 = grad_0_beta0[pp];
                double DENDRO_2 = grad_1_beta1[pp];
                double DENDRO_3 = (2.0/3.0)*gt0[pp];
                double DENDRO_4 = grad_2_beta2[pp];
                double DENDRO_5 = grad_0_beta1[pp];
                double DENDRO_6 = 2*gt1[pp];
                double DENDRO_7 = grad_0_beta2[pp];
                double DENDRO_8 = 2*gt2[pp];
                double DENDRO_9 = grad_1_beta0[pp];
                double DENDRO_10 = grad_1_beta2[pp];
                double DENDRO_11 = (1.0/3.0)*gt1[pp];
                double DENDRO_12 = (2.0/3.0)*DENDRO_4;
                double DENDRO_13 = grad_2_beta0[pp];
                double DENDRO_14 = grad_2_beta1[pp];
                double DENDRO_15 = (1.0/3.0)*gt2[pp];
                double DENDRO_16 = (2.0/3.0)*DENDRO_2;
                double DENDRO_17 = (2.0/3.0)*DENDRO_1;
                double DENDRO_18 = 2*gt4[pp];
                double DENDRO_19 = (1.0/3.0)*gt4[pp];

                // Dendro: printing variables
                //--
                gt_rhs00[pp] = -At0[pp]*DENDRO_0 + (4.0/3.0)*DENDRO_1*gt0[pp] - DENDRO_2*DENDRO_3 - DENDRO_3*DENDRO_4 + DENDRO_5*DENDRO_6 + DENDRO_7*DENDRO_8 + beta0[pp]*agrad_0_gt0[pp] + beta1[pp]*agrad_1_gt0[pp] + beta2[pp]*agrad_2_gt0[pp];
                //--
                gt_rhs01[pp] = -At1[pp]*DENDRO_0 + DENDRO_1*DENDRO_11 + DENDRO_10*gt2[pp] + DENDRO_11*DENDRO_2 - DENDRO_12*gt1[pp] + DENDRO_5*gt3[pp] + DENDRO_7*gt4[pp] + DENDRO_9*gt0[pp] + beta0[pp]*agrad_0_gt1[pp] + beta1[pp]*agrad_1_gt1[pp] + beta2[pp]*agrad_2_gt1[pp];
                //--
                gt_rhs02[pp] = -At2[pp]*DENDRO_0 + DENDRO_1*DENDRO_15 + DENDRO_13*gt0[pp] + DENDRO_14*gt1[pp] + DENDRO_15*DENDRO_4 - DENDRO_16*gt2[pp] + DENDRO_5*gt4[pp] + DENDRO_7*gt5[pp] + beta0[pp]*agrad_0_gt2[pp] + beta1[pp]*agrad_1_gt2[pp] + beta2[pp]*agrad_2_gt2[pp];
                //--
                gt_rhs11[pp] = -At3[pp]*DENDRO_0 + DENDRO_10*DENDRO_18 - DENDRO_12*gt3[pp] - DENDRO_17*gt3[pp] + (4.0/3.0)*DENDRO_2*gt3[pp] + DENDRO_6*DENDRO_9 + beta0[pp]*agrad_0_gt3[pp] + beta1[pp]*agrad_1_gt3[pp] + beta2[pp]*agrad_2_gt3[pp];
                //--
                gt_rhs12[pp] = -At4[pp]*DENDRO_0 + DENDRO_10*gt5[pp] + DENDRO_13*gt1[pp] + DENDRO_14*gt3[pp] - DENDRO_17*gt4[pp] + DENDRO_19*DENDRO_2 + DENDRO_19*DENDRO_4 + DENDRO_9*gt2[pp] + beta0[pp]*agrad_0_gt4[pp] + beta1[pp]*agrad_1_gt4[pp] + beta2[pp]*agrad_2_gt4[pp];
                //--
                gt_rhs22[pp] = -At5[pp]*DENDRO_0 + DENDRO_13*DENDRO_8 + DENDRO_14*DENDRO_18 - DENDRO_16*gt5[pp] - DENDRO_17*gt5[pp] + (4.0/3.0)*DENDRO_4*gt5[pp] + beta0[pp]*agrad_0_gt5[pp] + beta1[pp]*agrad_1_gt5[pp] + beta2[pp]*agrad_2_gt5[pp];
                // Dendro: reduced ops: 162
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }

    bssn::timer::t_rhs_gt.stop();

    bssn::timer::t_rhs_chi.start();
    // chi_rhs
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.chi_rhs]
                vnames = ['chi_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 22 
                // Dendro: printing temp variables
                double DENDRO_0 = (2.0/3.0)*chi[pp];

                // Dendro: printing variables
                //--
                chi_rhs[pp] = DENDRO_0*K[pp]*alpha[pp] - DENDRO_0*(grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp]) + beta0[pp]*agrad_0_chi[pp] + beta1[pp]*agrad_1_chi[pp] + beta2[pp]*agrad_2_chi[pp];
                // Dendro: reduced ops: 20
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }

    bssn::timer::t_rhs_chi.stop();

    bssn::timer::t_rhs_At.start();
    // At_rhs
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.At_rhs]
                vnames = ['At_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 630012 
                // Dendro: printing temp variables
                double DENDRO_0 = grad_0_beta0[pp];
                double DENDRO_1 = grad_1_beta1[pp];
                double DENDRO_2 = (2.0/3.0)*At0[pp];
                double DENDRO_3 = grad_2_beta2[pp];
                double DENDRO_4 = grad_0_beta1[pp];
                double DENDRO_5 = 2*At1[pp];
                double DENDRO_6 = grad_0_beta2[pp];
                double DENDRO_7 = 2*At2[pp];
                double DENDRO_8 = gt2[pp]*gt4[pp];
                double DENDRO_9 = -DENDRO_8 + gt1[pp]*gt5[pp];
                double DENDRO_10 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
                double DENDRO_11 = gt0[pp]*gt5[pp];
                double DENDRO_12 = pow(gt2[pp], 2);
                double DENDRO_13 = DENDRO_11 - DENDRO_12;
                double DENDRO_14 = pow(gt4[pp], 2);
                double DENDRO_15 = pow(gt1[pp], 2);
                double DENDRO_16 = gt3[pp]*gt5[pp];
                double DENDRO_17 = DENDRO_12*gt3[pp] + DENDRO_14*gt0[pp] + DENDRO_15*gt5[pp] - DENDRO_16*gt0[pp] - 2*DENDRO_8*gt1[pp];
                double DENDRO_18 = 1.0/DENDRO_17;
                double DENDRO_19 = DENDRO_18*DENDRO_5;
                double DENDRO_20 = At1[pp]*DENDRO_9;
                double DENDRO_21 = gt1[pp]*gt4[pp];
                double DENDRO_22 = gt2[pp]*gt3[pp];
                double DENDRO_23 = DENDRO_21 - DENDRO_22;
                double DENDRO_24 = -At2[pp]*DENDRO_23;
                double DENDRO_25 = -DENDRO_14 + DENDRO_16;
                double DENDRO_26 = 2*DENDRO_18;
                double DENDRO_27 = At0[pp]*DENDRO_26;
                double DENDRO_28 = gt0[pp]*gt3[pp];
                double DENDRO_29 = -DENDRO_15 + DENDRO_28;
                double DENDRO_30 = DENDRO_18*DENDRO_7;
                double DENDRO_31 = grad2_0_0_alpha[pp];
                double DENDRO_32 = 1.0/chi[pp];
                double DENDRO_33 = grad_1_chi[pp];
                double DENDRO_34 = grad_2_chi[pp];
                double DENDRO_35 = grad_0_chi[pp];
                double DENDRO_36 = DENDRO_10*DENDRO_34 + DENDRO_35*DENDRO_9;
                double DENDRO_37 = -DENDRO_13*DENDRO_33 + DENDRO_36;
                double DENDRO_38 = DENDRO_32*DENDRO_37;
                double DENDRO_39 = 0.5*gt0[pp];
                double DENDRO_40 = grad_0_gt1[pp];
                double DENDRO_41 = 1.0*DENDRO_40;
                double DENDRO_42 = grad_1_gt0[pp];
                double DENDRO_43 = 0.5*DENDRO_42;
                double DENDRO_44 = DENDRO_41 - DENDRO_43;
                double DENDRO_45 = grad_0_gt2[pp];
                double DENDRO_46 = 1.0*DENDRO_45;
                double DENDRO_47 = grad_2_gt0[pp];
                double DENDRO_48 = 0.5*DENDRO_47;
                double DENDRO_49 = DENDRO_46 - DENDRO_48;
                double DENDRO_50 = grad_0_gt0[pp];
                double DENDRO_51 = 0.5*DENDRO_50;
                double DENDRO_52 = DENDRO_10*DENDRO_49 + DENDRO_51*DENDRO_9;
                double DENDRO_53 = -DENDRO_13*DENDRO_44 + DENDRO_52;
                double DENDRO_54 = DENDRO_38*DENDRO_39 + DENDRO_53;
                double DENDRO_55 = grad_1_alpha[pp];
                double DENDRO_56 = 12*DENDRO_55;
                double DENDRO_57 = DENDRO_18*DENDRO_56;
                double DENDRO_58 = DENDRO_10*DENDRO_33;
                double DENDRO_59 = DENDRO_23*DENDRO_35;
                double DENDRO_60 = DENDRO_29*DENDRO_34;
                double DENDRO_61 = DENDRO_32*(DENDRO_58 - DENDRO_59 - DENDRO_60);
                double DENDRO_62 = DENDRO_39*DENDRO_61;
                double DENDRO_63 = DENDRO_23*DENDRO_51;
                double DENDRO_64 = DENDRO_29*DENDRO_49;
                double DENDRO_65 = DENDRO_10*DENDRO_44;
                double DENDRO_66 = DENDRO_63 + DENDRO_64 - DENDRO_65;
                double DENDRO_67 = grad_2_alpha[pp];
                double DENDRO_68 = 12*DENDRO_67;
                double DENDRO_69 = DENDRO_18*DENDRO_68;
                double DENDRO_70 = DENDRO_25*DENDRO_51;
                double DENDRO_71 = DENDRO_18*DENDRO_70;
                double DENDRO_72 = DENDRO_23*DENDRO_49;
                double DENDRO_73 = DENDRO_18*DENDRO_72;
                double DENDRO_74 = DENDRO_44*DENDRO_9;
                double DENDRO_75 = DENDRO_18*DENDRO_74;
                double DENDRO_76 = -DENDRO_21 + DENDRO_22;
                double DENDRO_77 = DENDRO_33*DENDRO_9;
                double DENDRO_78 = DENDRO_14 - DENDRO_16;
                double DENDRO_79 = DENDRO_34*DENDRO_76 + DENDRO_35*DENDRO_78 + DENDRO_77;
                double DENDRO_80 = DENDRO_18*gt0[pp];
                double DENDRO_81 = DENDRO_32*(-1.0*DENDRO_35 + 0.5*DENDRO_79*DENDRO_80);
                double DENDRO_82 = grad_0_alpha[pp];
                double DENDRO_83 = 12*DENDRO_82;
                double DENDRO_84 = grad2_0_0_chi[pp];
                double DENDRO_85 = -DENDRO_84;
                double DENDRO_86 = -DENDRO_63 - DENDRO_64 + DENDRO_65;
                double DENDRO_87 = DENDRO_18*DENDRO_34;
                double DENDRO_88 = DENDRO_18*DENDRO_33;
                double DENDRO_89 = -DENDRO_70 - DENDRO_72 + DENDRO_74;
                double DENDRO_90 = DENDRO_18*DENDRO_35;
                double DENDRO_91 = 2*DENDRO_32;
                double DENDRO_92 = grad_0_gt3[pp];
                double DENDRO_93 = 0.5*DENDRO_92;
                double DENDRO_94 = grad_1_gt1[pp];
                double DENDRO_95 = 1.0*DENDRO_94;
                double DENDRO_96 = -DENDRO_95;
                double DENDRO_97 = DENDRO_93 + DENDRO_96;
                double DENDRO_98 = grad_0_gt4[pp];
                double DENDRO_99 = grad_2_gt1[pp];
                double DENDRO_100 = grad_1_gt2[pp];
                double DENDRO_101 = -DENDRO_100 + DENDRO_98 + DENDRO_99;
                double DENDRO_102 = grad_0_gt5[pp];
                double DENDRO_103 = DENDRO_10*DENDRO_102 + DENDRO_47*DENDRO_9;
                double DENDRO_104 = -DENDRO_101*DENDRO_13 + DENDRO_103;
                double DENDRO_105 = DENDRO_104*DENDRO_97;
                double DENDRO_106 = DENDRO_13*DENDRO_92;
                double DENDRO_107 = DENDRO_100 + DENDRO_98 - DENDRO_99;
                double DENDRO_108 = DENDRO_10*DENDRO_107;
                double DENDRO_109 = DENDRO_42*DENDRO_9;
                double DENDRO_110 = DENDRO_108 + DENDRO_109;
                double DENDRO_111 = -DENDRO_106 + DENDRO_110;
                double DENDRO_112 = DENDRO_101*DENDRO_111;
                double DENDRO_113 = 0.25*DENDRO_112;
                double DENDRO_114 = pow(DENDRO_17, -2);
                double DENDRO_115 = DENDRO_10*DENDRO_114;
                double DENDRO_116 = 4*DENDRO_115;
                double DENDRO_117 = 0.5*DENDRO_102;
                double DENDRO_118 = grad_2_gt2[pp];
                double DENDRO_119 = 1.0*DENDRO_118;
                double DENDRO_120 = -DENDRO_119;
                double DENDRO_121 = DENDRO_117 + DENDRO_120;
                double DENDRO_122 = DENDRO_10*DENDRO_92;
                double DENDRO_123 = DENDRO_23*DENDRO_42;
                double DENDRO_124 = DENDRO_107*DENDRO_29;
                double DENDRO_125 = DENDRO_122 - DENDRO_123 - DENDRO_124;
                double DENDRO_126 = DENDRO_121*DENDRO_125;
                double DENDRO_127 = DENDRO_23*DENDRO_47;
                double DENDRO_128 = DENDRO_102*DENDRO_29;
                double DENDRO_129 = DENDRO_10*DENDRO_101;
                double DENDRO_130 = -DENDRO_127 - DENDRO_128 + DENDRO_129;
                double DENDRO_131 = DENDRO_107*DENDRO_130;
                double DENDRO_132 = 0.25*DENDRO_131;
                double DENDRO_133 = -DENDRO_132;
                double DENDRO_134 = -DENDRO_122 + DENDRO_123 + DENDRO_124;
                double DENDRO_135 = 0.25*DENDRO_102;
                double DENDRO_136 = DENDRO_134*DENDRO_135;
                double DENDRO_137 = DENDRO_100 - DENDRO_98 + DENDRO_99;
                double DENDRO_138 = DENDRO_127 + DENDRO_128 - DENDRO_129;
                double DENDRO_139 = 0.5*DENDRO_138;
                double DENDRO_140 = DENDRO_9*DENDRO_92;
                double DENDRO_141 = DENDRO_25*DENDRO_42;
                double DENDRO_142 = DENDRO_107*DENDRO_23;
                double DENDRO_143 = -DENDRO_140 + DENDRO_141 + DENDRO_142;
                double DENDRO_144 = 0.25*DENDRO_47;
                double DENDRO_145 = DENDRO_143*DENDRO_144;
                double DENDRO_146 = DENDRO_102*DENDRO_23;
                double DENDRO_147 = DENDRO_25*DENDRO_47;
                double DENDRO_148 = DENDRO_101*DENDRO_9;
                double DENDRO_149 = DENDRO_146 + DENDRO_147 - DENDRO_148;
                double DENDRO_150 = 0.25*DENDRO_42;
                double DENDRO_151 = DENDRO_149*DENDRO_150;
                double DENDRO_152 = 2*DENDRO_121;
                double DENDRO_153 = DENDRO_114*DENDRO_23;
                double DENDRO_154 = 4*DENDRO_153;
                double DENDRO_155 = DENDRO_143*DENDRO_50;
                double DENDRO_156 = 0.25*DENDRO_155;
                double DENDRO_157 = DENDRO_70 + DENDRO_72 - DENDRO_74;
                double DENDRO_158 = DENDRO_157*DENDRO_42;
                double DENDRO_159 = DENDRO_114*DENDRO_9;
                double DENDRO_160 = 4*DENDRO_159;
                double DENDRO_161 = 0.5*DENDRO_49;
                double DENDRO_162 = DENDRO_149*DENDRO_50;
                double DENDRO_163 = 0.25*DENDRO_162;
                double DENDRO_164 = DENDRO_157*DENDRO_47;
                double DENDRO_165 = DENDRO_138*DENDRO_47;
                double DENDRO_166 = DENDRO_102*DENDRO_66;
                double DENDRO_167 = 2.0*DENDRO_153;
                double DENDRO_168 = DENDRO_134*DENDRO_47;
                double DENDRO_169 = DENDRO_107*DENDRO_66;
                double DENDRO_170 = 2.0*DENDRO_159;
                double DENDRO_171 = DENDRO_104*DENDRO_23;
                double DENDRO_172 = grad_2_gt3[pp];
                double DENDRO_173 = DENDRO_13*DENDRO_172;
                double DENDRO_174 = grad_1_gt5[pp];
                double DENDRO_175 = DENDRO_10*DENDRO_174;
                double DENDRO_176 = DENDRO_137*DENDRO_9;
                double DENDRO_177 = DENDRO_175 + DENDRO_176;
                double DENDRO_178 = -DENDRO_173 + DENDRO_177;
                double DENDRO_179 = DENDRO_10*DENDRO_178;
                double DENDRO_180 = DENDRO_111*DENDRO_9;
                double DENDRO_181 = grad_2_gt5[pp];
                double DENDRO_182 = 0.5*DENDRO_181;
                double DENDRO_183 = DENDRO_10*DENDRO_182;
                double DENDRO_184 = 0.5*DENDRO_174;
                double DENDRO_185 = grad_2_gt4[pp];
                double DENDRO_186 = 1.0*DENDRO_185;
                double DENDRO_187 = -DENDRO_186;
                double DENDRO_188 = DENDRO_184 + DENDRO_187;
                double DENDRO_189 = -DENDRO_121*DENDRO_9 + DENDRO_13*DENDRO_188 + DENDRO_183;
                double DENDRO_190 = DENDRO_189*DENDRO_29;
                double DENDRO_191 = grad_1_gt3[pp];
                double DENDRO_192 = 0.5*DENDRO_191;
                double DENDRO_193 = DENDRO_13*DENDRO_192;
                double DENDRO_194 = grad_1_gt4[pp];
                double DENDRO_195 = 1.0*DENDRO_194;
                double DENDRO_196 = 0.5*DENDRO_172;
                double DENDRO_197 = DENDRO_195 - DENDRO_196;
                double DENDRO_198 = DENDRO_10*DENDRO_197;
                double DENDRO_199 = DENDRO_9*DENDRO_97;
                double DENDRO_200 = -DENDRO_193 + DENDRO_198 - DENDRO_199;
                double DENDRO_201 = DENDRO_13*DENDRO_200;
                double DENDRO_202 = DENDRO_25*DENDRO_53;
                double DENDRO_203 = DENDRO_171 - 1.0*DENDRO_179 - 1.0*DENDRO_180 + DENDRO_190 + DENDRO_201 + DENDRO_202;
                double DENDRO_204 = 2.0*DENDRO_114;
                double DENDRO_205 = DENDRO_204*DENDRO_42;
                double DENDRO_206 = -DENDRO_146 - DENDRO_147 + DENDRO_148;
                double DENDRO_207 = DENDRO_206*DENDRO_23;
                double DENDRO_208 = DENDRO_172*DENDRO_9;
                double DENDRO_209 = DENDRO_174*DENDRO_23;
                double DENDRO_210 = DENDRO_137*DENDRO_25;
                double DENDRO_211 = DENDRO_208 - DENDRO_209 - DENDRO_210;
                double DENDRO_212 = DENDRO_10*DENDRO_211;
                double DENDRO_213 = DENDRO_140 - DENDRO_141 - DENDRO_142;
                double DENDRO_214 = DENDRO_213*DENDRO_9;
                double DENDRO_215 = DENDRO_182*DENDRO_23;
                double DENDRO_216 = DENDRO_188*DENDRO_9;
                double DENDRO_217 = DENDRO_121*DENDRO_25;
                double DENDRO_218 = -DENDRO_215 - DENDRO_216 + DENDRO_217;
                double DENDRO_219 = DENDRO_218*DENDRO_29;
                double DENDRO_220 = DENDRO_192*DENDRO_9;
                double DENDRO_221 = -DENDRO_197*DENDRO_23 + DENDRO_220 + DENDRO_25*DENDRO_97;
                double DENDRO_222 = DENDRO_13*DENDRO_221;
                double DENDRO_223 = DENDRO_25*DENDRO_89;
                double DENDRO_224 = DENDRO_207 - 1.0*DENDRO_212 - 1.0*DENDRO_214 + DENDRO_219 + DENDRO_222 + DENDRO_223;
                double DENDRO_225 = DENDRO_204*DENDRO_50;
                double DENDRO_226 = DENDRO_130*DENDRO_23;
                double DENDRO_227 = DENDRO_10*DENDRO_172;
                double DENDRO_228 = DENDRO_174*DENDRO_29;
                double DENDRO_229 = DENDRO_137*DENDRO_23;
                double DENDRO_230 = DENDRO_227 - DENDRO_228 - DENDRO_229;
                double DENDRO_231 = DENDRO_10*DENDRO_230;
                double DENDRO_232 = DENDRO_125*DENDRO_9;
                double DENDRO_233 = DENDRO_182*DENDRO_29;
                double DENDRO_234 = DENDRO_121*DENDRO_23;
                double DENDRO_235 = DENDRO_10*DENDRO_188;
                double DENDRO_236 = -DENDRO_233 + DENDRO_234 - DENDRO_235;
                double DENDRO_237 = DENDRO_236*DENDRO_29;
                double DENDRO_238 = DENDRO_10*DENDRO_192;
                double DENDRO_239 = -DENDRO_197*DENDRO_29 + DENDRO_23*DENDRO_97 + DENDRO_238;
                double DENDRO_240 = DENDRO_13*DENDRO_239;
                double DENDRO_241 = DENDRO_25*DENDRO_86;
                double DENDRO_242 = DENDRO_226 - 1.0*DENDRO_231 - 1.0*DENDRO_232 + DENDRO_237 + DENDRO_240 + DENDRO_241;
                double DENDRO_243 = DENDRO_204*DENDRO_47;
                double DENDRO_244 = grad2_2_2_chi[pp];
                double DENDRO_245 = pow(DENDRO_34, 2);
                double DENDRO_246 = 3*DENDRO_32;
                double DENDRO_247 = DENDRO_29*(2*DENDRO_244 - DENDRO_245*DENDRO_246);
                double DENDRO_248 = grad2_1_1_chi[pp];
                double DENDRO_249 = pow(DENDRO_33, 2);
                double DENDRO_250 = DENDRO_13*(-DENDRO_246*DENDRO_249 + 2*DENDRO_248);
                double DENDRO_251 = pow(DENDRO_35, 2);
                double DENDRO_252 = DENDRO_25*(-DENDRO_246*DENDRO_251 + 2*DENDRO_84);
                double DENDRO_253 = grad2_1_2_chi[pp];
                double DENDRO_254 = DENDRO_246*DENDRO_34;
                double DENDRO_255 = 2*DENDRO_10*(2*DENDRO_253 - DENDRO_254*DENDRO_33);
                double DENDRO_256 = grad2_0_2_chi[pp];
                double DENDRO_257 = 2*DENDRO_23*(-DENDRO_254*DENDRO_35 + 2*DENDRO_256);
                double DENDRO_258 = grad2_0_1_chi[pp];
                double DENDRO_259 = DENDRO_33*DENDRO_35;
                double DENDRO_260 = 2*DENDRO_9*(-DENDRO_246*DENDRO_259 + 2*DENDRO_258);
                double DENDRO_261 = -1.0*DENDRO_171 + DENDRO_179 + DENDRO_180 - DENDRO_190 - DENDRO_201 - DENDRO_202;
                double DENDRO_262 = 2*DENDRO_261*DENDRO_88;
                double DENDRO_263 = -1.0*DENDRO_226 + DENDRO_231 + DENDRO_232 - DENDRO_237 - DENDRO_240 - DENDRO_241;
                double DENDRO_264 = 2*DENDRO_263*DENDRO_87;
                double DENDRO_265 = -1.0*DENDRO_207 + DENDRO_212 + DENDRO_214 - DENDRO_219 - DENDRO_222 - DENDRO_223;
                double DENDRO_266 = 2*DENDRO_265*DENDRO_90;
                double DENDRO_267 = DENDRO_247 + DENDRO_250 + DENDRO_252 - DENDRO_255 + DENDRO_257 - DENDRO_260 + DENDRO_262 + DENDRO_264 + DENDRO_266;
                double DENDRO_268 = DENDRO_32*DENDRO_80;
                double DENDRO_269 = DENDRO_114*DENDRO_29;
                double DENDRO_270 = 4*DENDRO_269;
                double DENDRO_271 = DENDRO_130*DENDRO_270;
                double DENDRO_272 = 0.25*DENDRO_92;
                double DENDRO_273 = DENDRO_114*DENDRO_13;
                double DENDRO_274 = 4*DENDRO_273;
                double DENDRO_275 = DENDRO_111*DENDRO_274;
                double DENDRO_276 = 0.25*DENDRO_98;
                double DENDRO_277 = -DENDRO_276;
                double DENDRO_278 = 0.75*DENDRO_100;
                double DENDRO_279 = 0.25*DENDRO_99;
                double DENDRO_280 = DENDRO_274*(DENDRO_277 + DENDRO_278 + DENDRO_279);
                double DENDRO_281 = DENDRO_46 + DENDRO_48;
                double DENDRO_282 = DENDRO_114*DENDRO_25;
                double DENDRO_283 = 4*DENDRO_282;
                double DENDRO_284 = DENDRO_149*DENDRO_47;
                double DENDRO_285 = 3.0*DENDRO_269;
                double DENDRO_286 = DENDRO_143*DENDRO_42;
                double DENDRO_287 = 3.0*DENDRO_273;
                double DENDRO_288 = 6.0*DENDRO_50;
                double DENDRO_289 = pow(chi[pp], -2);
                double DENDRO_290 = grad_0_Gt0[pp];
                double DENDRO_291 = grad_0_Gt1[pp];
                double DENDRO_292 = 4*gt1[pp];
                double DENDRO_293 = grad_0_Gt2[pp];
                double DENDRO_294 = 4*gt2[pp];
                double DENDRO_295 = 0.5*DENDRO_44;
                double DENDRO_296 = DENDRO_104*DENDRO_295;
                double DENDRO_297 = 4*DENDRO_18*DENDRO_23;
                double DENDRO_298 = DENDRO_104*DENDRO_272;
                double DENDRO_299 = 0.5*DENDRO_137;
                double DENDRO_300 = DENDRO_111*DENDRO_295;
                double DENDRO_301 = 2.0*DENDRO_18;
                double DENDRO_302 = DENDRO_29*DENDRO_301;
                double DENDRO_303 = DENDRO_13*DENDRO_301;
                double DENDRO_304 = DENDRO_25*DENDRO_301;
                double DENDRO_305 = DENDRO_111*DENDRO_42;
                double DENDRO_306 = DENDRO_53*DENDRO_92;
                double DENDRO_307 = DENDRO_104*DENDRO_42;
                double DENDRO_308 = DENDRO_101*DENDRO_53;
                double DENDRO_309 = 4.0*DENDRO_18;
                double DENDRO_310 = DENDRO_10*DENDRO_309;
                double DENDRO_311 = DENDRO_309*DENDRO_9;
                double DENDRO_312 = 0.25*DENDRO_100;
                double DENDRO_313 = 0.75*DENDRO_99;
                double DENDRO_314 = 4*DENDRO_114;
                double DENDRO_315 = -DENDRO_104*DENDRO_270*(DENDRO_277 + DENDRO_312 + DENDRO_313) + DENDRO_116*(DENDRO_111*DENDRO_299 + DENDRO_298) - DENDRO_154*(DENDRO_137*DENDRO_53 + DENDRO_296) + DENDRO_160*(DENDRO_300 - 2*DENDRO_53*DENDRO_97) - DENDRO_167*(DENDRO_307 + DENDRO_308) + DENDRO_170*(DENDRO_305 + DENDRO_306) - DENDRO_202*DENDRO_314*(DENDRO_41 + DENDRO_43) - DENDRO_251*DENDRO_289 + 4*DENDRO_290*gt0[pp] + DENDRO_291*DENDRO_292 + DENDRO_293*DENDRO_294 + DENDRO_297*grad2_0_2_gt0[pp] + DENDRO_302*grad2_2_2_gt0[pp] + DENDRO_303*grad2_1_1_gt0[pp] + DENDRO_304*grad2_0_0_gt0[pp] - DENDRO_310*grad2_1_2_gt0[pp] - DENDRO_311*grad2_0_1_gt0[pp];
                double DENDRO_316 = 3*alpha[pp];
                double DENDRO_317 = grad2_2_2_alpha[pp];
                double DENDRO_318 = 0.5*gt5[pp];
                double DENDRO_319 = DENDRO_189 + DENDRO_318*DENDRO_38;
                double DENDRO_320 = 4*DENDRO_55;
                double DENDRO_321 = DENDRO_18*DENDRO_320;
                double DENDRO_322 = DENDRO_23*DENDRO_34;
                double DENDRO_323 = DENDRO_25*DENDRO_35;
                double DENDRO_324 = DENDRO_32*(-DENDRO_322 - DENDRO_323 + DENDRO_77);
                double DENDRO_325 = DENDRO_318*DENDRO_324;
                double DENDRO_326 = 4*DENDRO_82;
                double DENDRO_327 = DENDRO_18*DENDRO_326;
                double DENDRO_328 = DENDRO_18*DENDRO_233;
                double DENDRO_329 = DENDRO_18*DENDRO_234;
                double DENDRO_330 = DENDRO_18*DENDRO_235;
                double DENDRO_331 = DENDRO_15 - DENDRO_28;
                double DENDRO_332 = DENDRO_331*DENDRO_34 + DENDRO_35*DENDRO_76 + DENDRO_58;
                double DENDRO_333 = DENDRO_18*gt5[pp];
                double DENDRO_334 = DENDRO_32*(0.5*DENDRO_332*DENDRO_333 - 1.0*DENDRO_34);
                double DENDRO_335 = 4*DENDRO_67;
                double DENDRO_336 = -DENDRO_244;
                double DENDRO_337 = -DENDRO_184;
                double DENDRO_338 = DENDRO_186 + DENDRO_337;
                double DENDRO_339 = -DENDRO_117;
                double DENDRO_340 = DENDRO_119 + DENDRO_339;
                double DENDRO_341 = -DENDRO_11 + DENDRO_12;
                double DENDRO_342 = 0.5*DENDRO_188;
                double DENDRO_343 = DENDRO_104*DENDRO_342;
                double DENDRO_344 = -DENDRO_343;
                double DENDRO_345 = DENDRO_107*DENDRO_189;
                double DENDRO_346 = 0.5*DENDRO_121;
                double DENDRO_347 = -DENDRO_206*DENDRO_346;
                double DENDRO_348 = 2*DENDRO_49;
                double DENDRO_349 = DENDRO_130*DENDRO_181;
                double DENDRO_350 = 0.25*DENDRO_349;
                double DENDRO_351 = DENDRO_102*DENDRO_236;
                double DENDRO_352 = DENDRO_211*DENDRO_346;
                double DENDRO_353 = -DENDRO_352;
                double DENDRO_354 = DENDRO_107*DENDRO_218;
                double DENDRO_355 = DENDRO_181*DENDRO_230;
                double DENDRO_356 = 0.25*DENDRO_355;
                double DENDRO_357 = DENDRO_174*DENDRO_236;
                double DENDRO_358 = DENDRO_211*DENDRO_49;
                double DENDRO_359 = DENDRO_137*DENDRO_206;
                double DENDRO_360 = 0.25*DENDRO_359;
                double DENDRO_361 = 0.25*DENDRO_174;
                double DENDRO_362 = DENDRO_130*DENDRO_361;
                double DENDRO_363 = DENDRO_135*DENDRO_230;
                double DENDRO_364 = DENDRO_144*DENDRO_211;
                double DENDRO_365 = 0.5*DENDRO_107;
                double DENDRO_366 = DENDRO_137*DENDRO_218;
                double DENDRO_367 = 2.0*DENDRO_115;
                double DENDRO_368 = DENDRO_204*DENDRO_261;
                double DENDRO_369 = DENDRO_204*DENDRO_263;
                double DENDRO_370 = DENDRO_204*DENDRO_265;
                double DENDRO_371 = DENDRO_218*DENDRO_47;
                double DENDRO_372 = -DENDRO_247 - DENDRO_250 - DENDRO_252 + DENDRO_255 - DENDRO_257 + DENDRO_260 - DENDRO_262 - DENDRO_264 - DENDRO_266;
                double DENDRO_373 = DENDRO_32*DENDRO_372;
                double DENDRO_374 = DENDRO_219*DENDRO_314;
                double DENDRO_375 = DENDRO_190*DENDRO_314;
                double DENDRO_376 = -DENDRO_279;
                double DENDRO_377 = DENDRO_274*(DENDRO_276 + DENDRO_278 + DENDRO_376);
                double DENDRO_378 = DENDRO_283*(-DENDRO_144 + DENDRO_46);
                double DENDRO_379 = DENDRO_174*DENDRO_230;
                double DENDRO_380 = DENDRO_102*DENDRO_130;
                double DENDRO_381 = 3.0*DENDRO_282;
                double DENDRO_382 = 6.0*DENDRO_114;
                double DENDRO_383 = grad_2_Gt0[pp];
                double DENDRO_384 = grad_2_Gt1[pp];
                double DENDRO_385 = 4*gt4[pp];
                double DENDRO_386 = grad_2_Gt2[pp];
                double DENDRO_387 = -DENDRO_178*DENDRO_342;
                double DENDRO_388 = DENDRO_104*DENDRO_197;
                double DENDRO_389 = DENDRO_101*DENDRO_178;
                double DENDRO_390 = 0.25*DENDRO_389;
                double DENDRO_391 = 0.25*DENDRO_172;
                double DENDRO_392 = DENDRO_104*DENDRO_391;
                double DENDRO_393 = DENDRO_174*DENDRO_178;
                double DENDRO_394 = DENDRO_172*DENDRO_189;
                double DENDRO_395 = DENDRO_104*DENDRO_174;
                double DENDRO_396 = DENDRO_101*DENDRO_189;
                double DENDRO_397 = 0.75*DENDRO_98;
                double DENDRO_398 = -DENDRO_104*DENDRO_283*(DENDRO_312 + DENDRO_376 + DENDRO_397) + DENDRO_116*(2*DENDRO_189*DENDRO_197 + DENDRO_387) + DENDRO_160*(DENDRO_388 + DENDRO_390) + DENDRO_160*(DENDRO_178*DENDRO_365 + DENDRO_392) - DENDRO_167*(DENDRO_395 + DENDRO_396) - DENDRO_178*DENDRO_274*(DENDRO_195 - DENDRO_391) - DENDRO_245*DENDRO_289 + DENDRO_294*DENDRO_383 + DENDRO_297*grad2_0_2_gt5[pp] + DENDRO_302*grad2_2_2_gt5[pp] + DENDRO_303*grad2_1_1_gt5[pp] + DENDRO_304*grad2_0_0_gt5[pp] - DENDRO_310*grad2_1_2_gt5[pp] - DENDRO_311*grad2_0_1_gt5[pp] + DENDRO_367*(DENDRO_393 + DENDRO_394) + DENDRO_384*DENDRO_385 + 4*DENDRO_386*gt5[pp];
                double DENDRO_399 = grad2_1_1_alpha[pp];
                double DENDRO_400 = 0.5*gt3[pp];
                double DENDRO_401 = DENDRO_18*DENDRO_335;
                double DENDRO_402 = -1.0*DENDRO_33;
                double DENDRO_403 = DENDRO_33*DENDRO_341 + DENDRO_36;
                double DENDRO_404 = DENDRO_18*gt3[pp];
                double DENDRO_405 = 0.5*DENDRO_404;
                double DENDRO_406 = -DENDRO_18*DENDRO_193 + DENDRO_18*DENDRO_198 - DENDRO_18*DENDRO_199;
                double DENDRO_407 = -DENDRO_248;
                double DENDRO_408 = -DENDRO_93;
                double DENDRO_409 = DENDRO_408 + DENDRO_95;
                double DENDRO_410 = DENDRO_211*DENDRO_44;
                double DENDRO_411 = DENDRO_137*DENDRO_213;
                double DENDRO_412 = 0.25*DENDRO_411;
                double DENDRO_413 = DENDRO_125*DENDRO_188;
                double DENDRO_414 = DENDRO_107*DENDRO_230;
                double DENDRO_415 = 0.25*DENDRO_414;
                double DENDRO_416 = DENDRO_125*DENDRO_361;
                double DENDRO_417 = 0.5*DENDRO_101;
                double DENDRO_418 = DENDRO_150*DENDRO_211;
                double DENDRO_419 = 0.5*DENDRO_97;
                double DENDRO_420 = DENDRO_211*DENDRO_419;
                double DENDRO_421 = -DENDRO_420;
                double DENDRO_422 = DENDRO_101*DENDRO_221;
                double DENDRO_423 = 0.5*DENDRO_197;
                double DENDRO_424 = DENDRO_230*DENDRO_423;
                double DENDRO_425 = 2*DENDRO_188*DENDRO_239;
                double DENDRO_426 = DENDRO_178*DENDRO_191;
                double DENDRO_427 = 0.25*DENDRO_426;
                double DENDRO_428 = DENDRO_172*DENDRO_200;
                double DENDRO_429 = DENDRO_125*DENDRO_423;
                double DENDRO_430 = DENDRO_101*DENDRO_239;
                double DENDRO_431 = -DENDRO_213*DENDRO_419;
                double DENDRO_432 = 2*DENDRO_221*DENDRO_44;
                double DENDRO_433 = DENDRO_111*DENDRO_191;
                double DENDRO_434 = 0.25*DENDRO_433;
                double DENDRO_435 = DENDRO_200*DENDRO_92;
                double DENDRO_436 = DENDRO_137*DENDRO_221;
                double DENDRO_437 = DENDRO_174*DENDRO_239;
                double DENDRO_438 = DENDRO_107*DENDRO_239;
                double DENDRO_439 = DENDRO_221*DENDRO_42;
                double DENDRO_440 = DENDRO_230*DENDRO_270;
                double DENDRO_441 = -DENDRO_312;
                double DENDRO_442 = DENDRO_270*(DENDRO_276 + DENDRO_313 + DENDRO_441);
                double DENDRO_443 = DENDRO_222*DENDRO_314;
                double DENDRO_444 = DENDRO_283*(-DENDRO_150 + DENDRO_41);
                double DENDRO_445 = DENDRO_283*(DENDRO_279 + DENDRO_397 + DENDRO_441);
                double DENDRO_446 = grad_1_Gt0[pp];
                double DENDRO_447 = grad_1_Gt1[pp];
                double DENDRO_448 = grad_1_Gt2[pp];
                double DENDRO_449 = DENDRO_111*DENDRO_391;
                double DENDRO_450 = DENDRO_178*DENDRO_272;
                double DENDRO_451 = DENDRO_172*DENDRO_178;
                double DENDRO_452 = DENDRO_111*DENDRO_92;
                double DENDRO_453 = -DENDRO_154*(DENDRO_111*DENDRO_196 + DENDRO_450) - DENDRO_154*(DENDRO_178*DENDRO_93 + DENDRO_449) - DENDRO_240*DENDRO_314*(DENDRO_195 + DENDRO_196) - DENDRO_249*DENDRO_289 - DENDRO_285*DENDRO_451 + DENDRO_292*DENDRO_446 + DENDRO_297*grad2_0_2_gt3[pp] + DENDRO_302*grad2_2_2_gt3[pp] + DENDRO_303*grad2_1_1_gt3[pp] + DENDRO_304*grad2_0_0_gt3[pp] - DENDRO_310*grad2_1_2_gt3[pp] - DENDRO_311*grad2_0_1_gt3[pp] - DENDRO_381*DENDRO_452 + DENDRO_385*DENDRO_448 + 4*DENDRO_447*gt3[pp];
                double DENDRO_454 = DENDRO_130*DENDRO_161;
                double DENDRO_455 = DENDRO_206*DENDRO_50;
                double DENDRO_456 = 0.25*DENDRO_455;
                double DENDRO_457 = DENDRO_47*DENDRO_89;
                double DENDRO_458 = DENDRO_125*DENDRO_135;
                double DENDRO_459 = DENDRO_150*DENDRO_206;
                double DENDRO_460 = DENDRO_144*DENDRO_213;
                double DENDRO_461 = DENDRO_125*DENDRO_161;
                double DENDRO_462 = DENDRO_213*DENDRO_50;
                double DENDRO_463 = 0.25*DENDRO_462;
                double DENDRO_464 = DENDRO_42*DENDRO_89;
                double DENDRO_465 = DENDRO_107*DENDRO_86;
                double DENDRO_466 = DENDRO_102*DENDRO_86;
                double DENDRO_467 = DENDRO_206*DENDRO_47;
                double DENDRO_468 = DENDRO_213*DENDRO_42;
                double DENDRO_469 = DENDRO_230*DENDRO_47;
                double DENDRO_470 = DENDRO_102*DENDRO_125;
                double DENDRO_471 = DENDRO_131 + DENDRO_470;
                double DENDRO_472 = DENDRO_178*DENDRO_42;
                double DENDRO_473 = DENDRO_107*DENDRO_111;
                double DENDRO_474 = DENDRO_104*DENDRO_92;
                double DENDRO_475 = DENDRO_473 + DENDRO_474;
                double DENDRO_476 = DENDRO_135*DENDRO_206;
                double DENDRO_477 = DENDRO_121*DENDRO_236;
                double DENDRO_478 = -DENDRO_477;
                double DENDRO_479 = DENDRO_189*DENDRO_299 + 0.25*DENDRO_395;
                double DENDRO_480 = DENDRO_107*DENDRO_213;
                double DENDRO_481 = 0.25*DENDRO_480;
                double DENDRO_482 = DENDRO_111*DENDRO_423;
                double DENDRO_483 = -DENDRO_178*DENDRO_419;
                double DENDRO_484 = DENDRO_482 + DENDRO_483;
                double DENDRO_485 = DENDRO_49*DENDRO_89;
                double DENDRO_486 = DENDRO_130*DENDRO_144;
                double DENDRO_487 = 0.25*DENDRO_307 + DENDRO_365*DENDRO_53;
                double DENDRO_488 = 0.25*DENDRO_104;
                double DENDRO_489 = DENDRO_137*DENDRO_488 + DENDRO_184*DENDRO_53;
                double DENDRO_490 = DENDRO_182*DENDRO_86;
                double DENDRO_491 = -DENDRO_130*DENDRO_346;
                double DENDRO_492 = DENDRO_101*DENDRO_488;
                double DENDRO_493 = DENDRO_107*DENDRO_488 + DENDRO_189*DENDRO_43;
                double DENDRO_494 = DENDRO_144*DENDRO_206;
                double DENDRO_495 = DENDRO_218*DENDRO_51;
                double DENDRO_496 = DENDRO_161*DENDRO_206 + DENDRO_495;
                double DENDRO_497 = DENDRO_189*DENDRO_97;
                double DENDRO_498 = 0.5*DENDRO_388;
                double DENDRO_499 = -DENDRO_497 + DENDRO_498;
                double DENDRO_500 = DENDRO_137*DENDRO_178;
                double DENDRO_501 = 0.25*DENDRO_500;
                double DENDRO_502 = DENDRO_189*DENDRO_93;
                double DENDRO_503 = DENDRO_111*DENDRO_361 + DENDRO_502;
                double DENDRO_504 = 0.25*DENDRO_107;
                double DENDRO_505 = DENDRO_218*DENDRO_43;
                double DENDRO_506 = DENDRO_206*DENDRO_504 + DENDRO_505;
                double DENDRO_507 = DENDRO_236*DENDRO_299;
                double DENDRO_508 = DENDRO_362 + DENDRO_363;
                double DENDRO_509 = DENDRO_230*DENDRO_346;
                double DENDRO_510 = -DENDRO_509;
                double DENDRO_511 = 0.25*DENDRO_125;
                double DENDRO_512 = DENDRO_181*DENDRO_511;
                double DENDRO_513 = DENDRO_236*DENDRO_365 + DENDRO_512;
                double DENDRO_514 = DENDRO_135*DENDRO_213 + DENDRO_505;
                double DENDRO_515 = DENDRO_178*DENDRO_295;
                double DENDRO_516 = -0.5*DENDRO_105 + DENDRO_197*DENDRO_53;
                double DENDRO_517 = DENDRO_184*DENDRO_86;
                double DENDRO_518 = DENDRO_161*DENDRO_230;
                double DENDRO_519 = 0.25*DENDRO_137;
                double DENDRO_520 = DENDRO_130*DENDRO_519;
                double DENDRO_521 = 0.25*DENDRO_211;
                double DENDRO_522 = DENDRO_50*DENDRO_521;
                double DENDRO_523 = DENDRO_161*DENDRO_213;
                double DENDRO_524 = DENDRO_522 + DENDRO_523;
                double DENDRO_525 = DENDRO_365*DENDRO_89 + DENDRO_459;
                double DENDRO_526 = DENDRO_137*DENDRO_230;
                double DENDRO_527 = DENDRO_125*DENDRO_174;
                double DENDRO_528 = DENDRO_414 + DENDRO_527;
                double DENDRO_529 = 1.0*DENDRO_273;
                double DENDRO_530 = -DENDRO_256;
                double DENDRO_531 = 0.5*DENDRO_87;
                double DENDRO_532 = 0.5*DENDRO_88;
                double DENDRO_533 = 0.5*DENDRO_90;
                double DENDRO_534 = 2.0*DENDRO_383;
                double DENDRO_535 = DENDRO_534*gt0[pp];
                double DENDRO_536 = 2.0*DENDRO_384;
                double DENDRO_537 = DENDRO_536*gt1[pp];
                double DENDRO_538 = 2.0*gt2[pp];
                double DENDRO_539 = DENDRO_290*DENDRO_538;
                double DENDRO_540 = DENDRO_386*DENDRO_538;
                double DENDRO_541 = 2.0*gt4[pp];
                double DENDRO_542 = DENDRO_291*DENDRO_541;
                double DENDRO_543 = 2.0*gt5[pp];
                double DENDRO_544 = DENDRO_293*DENDRO_543;
                double DENDRO_545 = DENDRO_289*DENDRO_34;
                double DENDRO_546 = -DENDRO_35*DENDRO_545;
                double DENDRO_547 = DENDRO_297*grad2_0_2_gt2[pp];
                double DENDRO_548 = DENDRO_302*grad2_2_2_gt2[pp];
                double DENDRO_549 = DENDRO_303*grad2_1_1_gt2[pp];
                double DENDRO_550 = DENDRO_304*grad2_0_0_gt2[pp];
                double DENDRO_551 = -DENDRO_310*grad2_1_2_gt2[pp];
                double DENDRO_552 = -DENDRO_311*grad2_0_1_gt2[pp];
                double DENDRO_553 = DENDRO_18*gt2[pp];
                double DENDRO_554 = DENDRO_100*DENDRO_368 + DENDRO_118*DENDRO_369 + DENDRO_370*DENDRO_45 + DENDRO_373*DENDRO_553 + DENDRO_535 + DENDRO_537 + DENDRO_539 + DENDRO_540 + DENDRO_542 + DENDRO_544 + DENDRO_546 + DENDRO_547 + DENDRO_548 + DENDRO_549 + DENDRO_550 + DENDRO_551 + DENDRO_552 - DENDRO_91*(DENDRO_530 + DENDRO_531*(DENDRO_102*DENDRO_331 + DENDRO_129 + DENDRO_47*DENDRO_76) + DENDRO_532*(DENDRO_101*DENDRO_341 + DENDRO_103) + DENDRO_533*(DENDRO_102*DENDRO_76 + DENDRO_148 + DENDRO_47*DENDRO_78));
                double DENDRO_555 = grad2_0_2_alpha[pp];
                double DENDRO_556 = DENDRO_127*DENDRO_18;
                double DENDRO_557 = DENDRO_128*DENDRO_18;
                double DENDRO_558 = DENDRO_129*DENDRO_18;
                double DENDRO_559 = -DENDRO_35;
                double DENDRO_560 = DENDRO_32*(DENDRO_332*DENDRO_553 + DENDRO_559);
                double DENDRO_561 = 2.0*DENDRO_67;
                double DENDRO_562 = DENDRO_146*DENDRO_18;
                double DENDRO_563 = DENDRO_147*DENDRO_18;
                double DENDRO_564 = DENDRO_148*DENDRO_18;
                double DENDRO_565 = -DENDRO_34;
                double DENDRO_566 = DENDRO_32*(DENDRO_553*DENDRO_79 + DENDRO_565);
                double DENDRO_567 = 2.0*DENDRO_82;
                double DENDRO_568 = 2.0*DENDRO_55;
                double DENDRO_569 = DENDRO_18*(DENDRO_104 + DENDRO_38*gt2[pp]);
                double DENDRO_570 = -4*DENDRO_555 + DENDRO_561*(-DENDRO_556 - DENDRO_557 + DENDRO_558 + DENDRO_560) + DENDRO_567*(-DENDRO_562 - DENDRO_563 + DENDRO_564 + DENDRO_566) + DENDRO_568*DENDRO_569;
                double DENDRO_571 = DENDRO_211*DENDRO_47;
                double DENDRO_572 = DENDRO_102*DENDRO_213 + DENDRO_571;
                double DENDRO_573 = 0.5*DENDRO_349;
                double DENDRO_574 = 0.25*DENDRO_526;
                double DENDRO_575 = DENDRO_130*DENDRO_135 + DENDRO_490;
                double DENDRO_576 = DENDRO_111*DENDRO_342;
                double DENDRO_577 = -DENDRO_576;
                double DENDRO_578 = -DENDRO_213*DENDRO_346;
                double DENDRO_579 = DENDRO_362 + DENDRO_512;
                double DENDRO_580 = DENDRO_458 + DENDRO_517;
                double DENDRO_581 = DENDRO_144*DENDRO_230;
                double DENDRO_582 = DENDRO_299*DENDRO_89;
                double DENDRO_583 = DENDRO_211*DENDRO_42;
                double DENDRO_584 = DENDRO_480 + DENDRO_583;
                double DENDRO_585 = DENDRO_104*DENDRO_172;
                double DENDRO_586 = DENDRO_111*DENDRO_174 + DENDRO_585;
                double DENDRO_587 = 0.25*DENDRO_473;
                double DENDRO_588 = DENDRO_196*DENDRO_53;
                double DENDRO_589 = DENDRO_150*DENDRO_178 + DENDRO_588;
                double DENDRO_590 = DENDRO_115*(DENDRO_500 + DENDRO_586) - DENDRO_154*(DENDRO_489 + DENDRO_492) - DENDRO_154*(-DENDRO_188*DENDRO_53 + DENDRO_493) + DENDRO_160*(DENDRO_113 + DENDRO_516) + DENDRO_160*(DENDRO_587 + DENDRO_589) - DENDRO_270*(DENDRO_344 + DENDRO_479) - DENDRO_274*(DENDRO_449 + DENDRO_484) - DENDRO_283*(0.5*DENDRO_308 + DENDRO_487);
                double DENDRO_591 = DENDRO_130*DENDRO_172;
                double DENDRO_592 = DENDRO_206*DENDRO_92;
                double DENDRO_593 = 0.25*DENDRO_393;
                double DENDRO_594 = DENDRO_188*DENDRO_236;
                double DENDRO_595 = -DENDRO_594;
                double DENDRO_596 = DENDRO_135*DENDRO_211 + DENDRO_218*DENDRO_417;
                double DENDRO_597 = DENDRO_197*DENDRO_200;
                double DENDRO_598 = DENDRO_230*DENDRO_391;
                double DENDRO_599 = DENDRO_221*DENDRO_365;
                double DENDRO_600 = DENDRO_211*DENDRO_272 + DENDRO_599;
                double DENDRO_601 = DENDRO_206*DENDRO_295;
                double DENDRO_602 = DENDRO_523 + DENDRO_601;
                double DENDRO_603 = DENDRO_218*DENDRO_44 + 0.5*DENDRO_358;
                double DENDRO_604 = DENDRO_101*DENDRO_206;
                double DENDRO_605 = 0.25*DENDRO_604;
                double DENDRO_606 = DENDRO_178*DENDRO_504 + DENDRO_502;
                double DENDRO_607 = DENDRO_236*DENDRO_417;
                double DENDRO_608 = DENDRO_130*DENDRO_342;
                double DENDRO_609 = -DENDRO_608;
                double DENDRO_610 = DENDRO_182*DENDRO_239;
                double DENDRO_611 = -DENDRO_230*DENDRO_342 + DENDRO_610;
                double DENDRO_612 = DENDRO_117*DENDRO_221;
                double DENDRO_613 = DENDRO_101*DENDRO_521 + DENDRO_612;
                double DENDRO_614 = DENDRO_211*DENDRO_519;
                double DENDRO_615 = DENDRO_211*DENDRO_504 + DENDRO_218*DENDRO_93;
                double DENDRO_616 = DENDRO_178*DENDRO_423;
                double DENDRO_617 = DENDRO_189*DENDRO_192;
                double DENDRO_618 = DENDRO_178*DENDRO_391 + DENDRO_617;
                double DENDRO_619 = -DENDRO_206*DENDRO_419;
                double DENDRO_620 = DENDRO_221*DENDRO_49;
                double DENDRO_621 = 0.5*DENDRO_410 + DENDRO_620;
                double DENDRO_622 = DENDRO_130*DENDRO_423;
                double DENDRO_623 = 0.25*DENDRO_101;
                double DENDRO_624 = DENDRO_117*DENDRO_239;
                double DENDRO_625 = DENDRO_230*DENDRO_623 + DENDRO_624;
                double DENDRO_626 = DENDRO_191*DENDRO_488;
                double DENDRO_627 = DENDRO_450 + DENDRO_626;
                double DENDRO_628 = DENDRO_200*DENDRO_365;
                double DENDRO_629 = DENDRO_101*DENDRO_130;
                double DENDRO_630 = 1.0*DENDRO_282;
                double DENDRO_631 = -DENDRO_253;
                double DENDRO_632 = DENDRO_18*gt4[pp];
                double DENDRO_633 = DENDRO_297*grad2_0_2_gt4[pp] + DENDRO_302*grad2_2_2_gt4[pp] + DENDRO_303*grad2_1_1_gt4[pp] + DENDRO_304*grad2_0_0_gt4[pp] - DENDRO_310*grad2_1_2_gt4[pp] - DENDRO_311*grad2_0_1_gt4[pp] - DENDRO_33*DENDRO_545 + DENDRO_386*DENDRO_541 + DENDRO_446*DENDRO_538 + DENDRO_447*DENDRO_541 + DENDRO_448*DENDRO_543 + DENDRO_534*gt1[pp] + DENDRO_536*gt3[pp];
                double DENDRO_634 = DENDRO_185*DENDRO_369 + DENDRO_194*DENDRO_368 + DENDRO_370*DENDRO_98 + DENDRO_373*DENDRO_632 + DENDRO_633 - DENDRO_91*(DENDRO_531*(DENDRO_137*DENDRO_76 + DENDRO_174*DENDRO_331 + DENDRO_227) + DENDRO_532*(DENDRO_172*DENDRO_341 + DENDRO_177) + DENDRO_533*(DENDRO_137*DENDRO_78 + DENDRO_174*DENDRO_76 + DENDRO_208) + DENDRO_631);
                double DENDRO_635 = grad2_1_2_alpha[pp];
                double DENDRO_636 = DENDRO_18*DENDRO_227;
                double DENDRO_637 = DENDRO_18*DENDRO_228;
                double DENDRO_638 = DENDRO_18*DENDRO_229;
                double DENDRO_639 = -DENDRO_33;
                double DENDRO_640 = DENDRO_32*(DENDRO_332*DENDRO_632 + DENDRO_639);
                double DENDRO_641 = -DENDRO_173*DENDRO_18 + DENDRO_175*DENDRO_18 + DENDRO_176*DENDRO_18;
                double DENDRO_642 = DENDRO_324*gt4[pp];
                double DENDRO_643 = DENDRO_18*DENDRO_567*(DENDRO_211 + DENDRO_642) + DENDRO_561*(DENDRO_636 - DENDRO_637 - DENDRO_638 + DENDRO_640) + DENDRO_568*(DENDRO_32*(DENDRO_403*DENDRO_632 + DENDRO_565) + DENDRO_641) - 4*DENDRO_635;
                double DENDRO_644 = 0.5*DENDRO_355;
                double DENDRO_645 = 1.0*DENDRO_437;
                double DENDRO_646 = 0.5*DENDRO_436;
                double DENDRO_647 = 0.25*DENDRO_629;
                double DENDRO_648 = DENDRO_363 + DENDRO_512;
                double DENDRO_649 = DENDRO_121*DENDRO_221;
                double DENDRO_650 = DENDRO_616 + DENDRO_617;
                double DENDRO_651 = DENDRO_230*DENDRO_361;
                double DENDRO_652 = DENDRO_221*DENDRO_48;
                double DENDRO_653 = DENDRO_206*DENDRO_272 + DENDRO_652;
                double DENDRO_654 = DENDRO_130*DENDRO_391 + DENDRO_624;
                double DENDRO_655 = DENDRO_449 + DENDRO_450;
                double DENDRO_656 = DENDRO_200*DENDRO_417 + DENDRO_626;
                double DENDRO_657 = 1.0*DENDRO_153;
                double DENDRO_658 = -DENDRO_154*(DENDRO_577 + DENDRO_606) - DENDRO_270*(DENDRO_189*DENDRO_196 + DENDRO_387 + DENDRO_593) - DENDRO_630*(DENDRO_112 + DENDRO_475) - DENDRO_657*(DENDRO_389 + DENDRO_586);
                double DENDRO_659 = DENDRO_510 + DENDRO_609;
                double DENDRO_660 = 0.5*DENDRO_433;
                double DENDRO_661 = DENDRO_200*DENDRO_97;
                double DENDRO_662 = -DENDRO_661;
                double DENDRO_663 = DENDRO_239*DENDRO_299;
                double DENDRO_664 = DENDRO_125*DENDRO_391 + DENDRO_663;
                double DENDRO_665 = DENDRO_213*DENDRO_272;
                double DENDRO_666 = DENDRO_221*DENDRO_43;
                double DENDRO_667 = DENDRO_44*DENDRO_89;
                double DENDRO_668 = DENDRO_125*DENDRO_144 + DENDRO_417*DENDRO_86;
                double DENDRO_669 = DENDRO_188*DENDRO_86;
                double DENDRO_670 = 0.5*DENDRO_126;
                double DENDRO_671 = -DENDRO_669 - DENDRO_670;
                double DENDRO_672 = DENDRO_417*DENDRO_89 + DENDRO_460;
                double DENDRO_673 = DENDRO_522 + DENDRO_601;
                double DENDRO_674 = DENDRO_121*DENDRO_239;
                double DENDRO_675 = 0.5*DENDRO_413;
                double DENDRO_676 = -DENDRO_674 - DENDRO_675;
                double DENDRO_677 = DENDRO_213*DENDRO_623 + DENDRO_652;
                double DENDRO_678 = DENDRO_200*DENDRO_299;
                double DENDRO_679 = DENDRO_449 + DENDRO_626;
                double DENDRO_680 = DENDRO_239*DENDRO_48;
                double DENDRO_681 = DENDRO_101*DENDRO_511 + DENDRO_680;
                double DENDRO_682 = DENDRO_221*DENDRO_51;
                double DENDRO_683 = DENDRO_213*DENDRO_295 + DENDRO_682;
                double DENDRO_684 = DENDRO_125*DENDRO_504;
                double DENDRO_685 = DENDRO_137*DENDRO_511 + DENDRO_196*DENDRO_86;
                double DENDRO_686 = DENDRO_192*DENDRO_53;
                double DENDRO_687 = DENDRO_111*DENDRO_272 + DENDRO_686;
                double DENDRO_688 = 1.0*DENDRO_269;
                double DENDRO_689 = -DENDRO_258;
                double DENDRO_690 = 2.0*DENDRO_446*gt0[pp];
                double DENDRO_691 = 2.0*gt1[pp];
                double DENDRO_692 = DENDRO_290*DENDRO_691;
                double DENDRO_693 = DENDRO_447*DENDRO_691;
                double DENDRO_694 = DENDRO_448*DENDRO_538;
                double DENDRO_695 = 2.0*DENDRO_291*gt3[pp];
                double DENDRO_696 = DENDRO_293*DENDRO_541;
                double DENDRO_697 = -DENDRO_259*DENDRO_289;
                double DENDRO_698 = DENDRO_297*grad2_0_2_gt1[pp];
                double DENDRO_699 = DENDRO_302*grad2_2_2_gt1[pp];
                double DENDRO_700 = DENDRO_303*grad2_1_1_gt1[pp];
                double DENDRO_701 = DENDRO_304*grad2_0_0_gt1[pp];
                double DENDRO_702 = -DENDRO_310*grad2_1_2_gt1[pp];
                double DENDRO_703 = -DENDRO_311*grad2_0_1_gt1[pp];
                double DENDRO_704 = DENDRO_18*gt1[pp];
                double DENDRO_705 = DENDRO_368*DENDRO_94 + DENDRO_369*DENDRO_99 + DENDRO_370*DENDRO_40 + DENDRO_373*DENDRO_704 + DENDRO_690 + DENDRO_692 + DENDRO_693 + DENDRO_694 + DENDRO_695 + DENDRO_696 + DENDRO_697 + DENDRO_698 + DENDRO_699 + DENDRO_700 + DENDRO_701 + DENDRO_702 + DENDRO_703 - DENDRO_91*(DENDRO_531*(DENDRO_107*DENDRO_331 + DENDRO_122 + DENDRO_42*DENDRO_76) + DENDRO_532*(DENDRO_110 + DENDRO_341*DENDRO_92) + DENDRO_533*(DENDRO_107*DENDRO_76 + DENDRO_140 + DENDRO_42*DENDRO_78) + DENDRO_689);
                double DENDRO_706 = 0.25*DENDRO_305;
                double DENDRO_707 = DENDRO_111*DENDRO_519 + DENDRO_588;
                double DENDRO_708 = -DENDRO_111*DENDRO_419;
                double DENDRO_709 = DENDRO_116*(DENDRO_483 + DENDRO_679) - DENDRO_154*(DENDRO_298 + DENDRO_589) - DENDRO_154*(DENDRO_298 + DENDRO_707) + DENDRO_160*(DENDRO_687 + DENDRO_708) - DENDRO_270*(DENDRO_104*DENDRO_196 + DENDRO_501) - DENDRO_283*(1.0*DENDRO_306 + DENDRO_706);
                double DENDRO_710 = grad2_0_1_alpha[pp];
                double DENDRO_711 = -DENDRO_106*DENDRO_18 + DENDRO_108*DENDRO_18 + DENDRO_109*DENDRO_18;
                double DENDRO_712 = DENDRO_140*DENDRO_18;
                double DENDRO_713 = DENDRO_141*DENDRO_18;
                double DENDRO_714 = DENDRO_142*DENDRO_18;
                double DENDRO_715 = DENDRO_32*(DENDRO_639 + DENDRO_704*DENDRO_79);
                double DENDRO_716 = DENDRO_61*gt1[pp];
                double DENDRO_717 = DENDRO_18*DENDRO_561*(DENDRO_125 + DENDRO_716) + DENDRO_567*(DENDRO_712 - DENDRO_713 - DENDRO_714 + DENDRO_715) + DENDRO_568*(DENDRO_32*(DENDRO_403*DENDRO_704 + DENDRO_559) + DENDRO_711) - 4*DENDRO_710;
                double DENDRO_718 = DENDRO_150*DENDRO_213;
                double DENDRO_719 = -DENDRO_10*(DENDRO_643 + alpha[pp]*(DENDRO_116*(DENDRO_611 + DENDRO_651) + DENDRO_116*(DENDRO_613 + DENDRO_614) + DENDRO_116*(DENDRO_615 - DENDRO_649) + DENDRO_116*(-DENDRO_188*DENDRO_200 + DENDRO_650) + DENDRO_116*(DENDRO_196*DENDRO_236 + DENDRO_610 + DENDRO_651) - DENDRO_154*(DENDRO_578 + DENDRO_603) - DENDRO_154*(DENDRO_607 + DENDRO_648) - DENDRO_154*(DENDRO_609 + DENDRO_648) + DENDRO_160*(DENDRO_412 + DENDRO_621) + DENDRO_160*(DENDRO_416 + DENDRO_625) + DENDRO_160*(DENDRO_416 + DENDRO_654) + DENDRO_160*(DENDRO_481 + DENDRO_653) + DENDRO_160*(DENDRO_482 + DENDRO_656) + DENDRO_160*(DENDRO_628 + DENDRO_655) - DENDRO_270*(DENDRO_353 + DENDRO_596) - DENDRO_270*(DENDRO_595 + DENDRO_644) - DENDRO_274*(DENDRO_598 + DENDRO_645) - DENDRO_274*(DENDRO_600 + DENDRO_646) - DENDRO_274*(DENDRO_196*DENDRO_200 + DENDRO_427 + DENDRO_597) - DENDRO_283*(DENDRO_460 + DENDRO_602) - DENDRO_283*(DENDRO_117*DENDRO_125 + DENDRO_647) + DENDRO_367*(DENDRO_174*DENDRO_200 + DENDRO_451) + DENDRO_634 - DENDRO_657*(DENDRO_572 + DENDRO_604) + DENDRO_658)) - DENDRO_10*(DENDRO_643 + alpha[pp]*(DENDRO_116*(DENDRO_614 + DENDRO_615) + DENDRO_116*(DENDRO_616 + DENDRO_618) + DENDRO_116*(DENDRO_184*DENDRO_200 + DENDRO_618) + DENDRO_116*(DENDRO_197*DENDRO_236 + DENDRO_611) + DENDRO_116*(-DENDRO_218*DENDRO_97 + DENDRO_613) - DENDRO_154*(DENDRO_360 + DENDRO_603) - DENDRO_154*(DENDRO_392 + DENDRO_503) - DENDRO_154*(DENDRO_392 + DENDRO_606) - DENDRO_154*(DENDRO_508 + DENDRO_607) - DENDRO_154*(DENDRO_513 + DENDRO_609) - DENDRO_154*(DENDRO_514 + DENDRO_605) + DENDRO_159*(DENDRO_528 + DENDRO_591) + DENDRO_159*(DENDRO_584 + DENDRO_592) + DENDRO_160*(DENDRO_482 + DENDRO_627) + DENDRO_160*(DENDRO_619 + DENDRO_621) + DENDRO_160*(DENDRO_622 + DENDRO_625) + DENDRO_160*(DENDRO_627 + DENDRO_628) - DENDRO_270*(0.5*DENDRO_366 + DENDRO_596) - DENDRO_270*(1.0*DENDRO_394 + DENDRO_593) - DENDRO_270*(DENDRO_184*DENDRO_236 + DENDRO_356 + DENDRO_595) - DENDRO_274*(DENDRO_421 + DENDRO_600) - DENDRO_274*(0.5*DENDRO_426 + DENDRO_597) - DENDRO_274*(DENDRO_184*DENDRO_239 + DENDRO_424 + DENDRO_598) - DENDRO_283*(DENDRO_459 + DENDRO_602) - DENDRO_283*(DENDRO_104*DENDRO_93 + DENDRO_587) + DENDRO_367*(DENDRO_172*DENDRO_236 + DENDRO_379) - DENDRO_630*(DENDRO_471 + DENDRO_629) + DENDRO_634)) + DENDRO_13*(DENDRO_320*(DENDRO_32*(DENDRO_402 + DENDRO_403*DENDRO_405) + DENDRO_406) + DENDRO_327*(DENDRO_221 + DENDRO_324*DENDRO_400) - 4*DENDRO_399 + DENDRO_401*(DENDRO_239 + DENDRO_400*DENDRO_61) + alpha[pp]*(DENDRO_116*(DENDRO_421 + DENDRO_422) + DENDRO_116*(DENDRO_424 - DENDRO_425) + DENDRO_116*(DENDRO_427 + 1.0*DENDRO_428) - DENDRO_125*DENDRO_445 - DENDRO_154*(DENDRO_410 + DENDRO_412) - DENDRO_154*(-1.0*DENDRO_413 + DENDRO_415) - DENDRO_154*(DENDRO_213*DENDRO_417 + DENDRO_418) - DENDRO_154*(DENDRO_230*DENDRO_417 + DENDRO_416) + DENDRO_160*(DENDRO_429 + DENDRO_430) + DENDRO_160*(DENDRO_431 + DENDRO_432) + DENDRO_160*(DENDRO_434 + 1.0*DENDRO_435) + DENDRO_170*(DENDRO_433 + DENDRO_435) + DENDRO_170*(DENDRO_125*DENDRO_172 + DENDRO_438) + DENDRO_170*(DENDRO_213*DENDRO_92 + DENDRO_439) + DENDRO_172*DENDRO_369 - DENDRO_191*DENDRO_201*DENDRO_382 + DENDRO_191*DENDRO_368 - DENDRO_211*DENDRO_442 - DENDRO_213*DENDRO_444 + DENDRO_367*(DENDRO_426 + DENDRO_428) + DENDRO_367*(DENDRO_172*DENDRO_230 + DENDRO_437) + DENDRO_367*(DENDRO_211*DENDRO_92 + DENDRO_436) + DENDRO_370*DENDRO_92 + DENDRO_373*DENDRO_404 - DENDRO_440*(DENDRO_186 - DENDRO_361) - DENDRO_443*(DENDRO_93 + DENDRO_95) + DENDRO_453 - DENDRO_91*(DENDRO_407 + DENDRO_87*(DENDRO_197*DENDRO_331 + DENDRO_238 + DENDRO_409*DENDRO_76) + DENDRO_88*(DENDRO_192*DENDRO_341 + DENDRO_198 + DENDRO_409*DENDRO_9) + DENDRO_90*(DENDRO_197*DENDRO_76 + DENDRO_220 + DENDRO_409*DENDRO_78)))) + DENDRO_23*(DENDRO_570 + alpha[pp]*(DENDRO_115*(DENDRO_359 + DENDRO_572) + DENDRO_116*(DENDRO_499 + DENDRO_577) + DENDRO_116*(DENDRO_506 + DENDRO_578) + DENDRO_116*(DENDRO_507 + DENDRO_579) + DENDRO_116*(DENDRO_510 + DENDRO_579) - DENDRO_154*(DENDRO_491 + DENDRO_575) - DENDRO_154*(-DENDRO_121*DENDRO_89 + DENDRO_496) - DENDRO_154*(DENDRO_236*DENDRO_48 + DENDRO_575) + DENDRO_160*(DENDRO_460 + DENDRO_525) + DENDRO_160*(DENDRO_520 + DENDRO_580) + DENDRO_160*(DENDRO_524 + DENDRO_582) + DENDRO_160*(DENDRO_580 + DENDRO_581) - DENDRO_167*(DENDRO_102*DENDRO_89 + DENDRO_467) - DENDRO_270*(DENDRO_478 + DENDRO_573) - DENDRO_270*(DENDRO_218*DENDRO_48 + DENDRO_347 + DENDRO_476) - DENDRO_274*(DENDRO_125*DENDRO_184 + DENDRO_574) - DENDRO_283*(1.0*DENDRO_466 + DENDRO_486) - DENDRO_283*(DENDRO_456 + DENDRO_48*DENDRO_89 + DENDRO_485) - DENDRO_529*(DENDRO_411 + DENDRO_584) + DENDRO_554 + DENDRO_590)) + DENDRO_23*(DENDRO_570 + alpha[pp]*(DENDRO_116*(DENDRO_364 + DENDRO_506) + DENDRO_116*(DENDRO_364 + DENDRO_514) + DENDRO_116*(DENDRO_390 + DENDRO_499) + DENDRO_116*(DENDRO_501 + DENDRO_503) + DENDRO_116*(DENDRO_507 + DENDRO_508) + DENDRO_116*(DENDRO_510 + DENDRO_513) - DENDRO_154*(DENDRO_492 + DENDRO_493) - DENDRO_154*(DENDRO_494 + DENDRO_496) - DENDRO_154*(DENDRO_189*DENDRO_44 + DENDRO_489) - DENDRO_154*(DENDRO_117*DENDRO_89 + DENDRO_494 + DENDRO_495) - DENDRO_154*(DENDRO_236*DENDRO_49 + DENDRO_490 + DENDRO_491) + DENDRO_159*(DENDRO_469 + DENDRO_471) + DENDRO_159*(DENDRO_472 + DENDRO_475) + DENDRO_160*(DENDRO_459 + DENDRO_524) + DENDRO_160*(DENDRO_515 + DENDRO_516) + DENDRO_160*(DENDRO_522 + DENDRO_525) + DENDRO_160*(DENDRO_517 + DENDRO_518 + DENDRO_520) - DENDRO_167*(DENDRO_236*DENDRO_47 + DENDRO_380) - DENDRO_270*(1.0*DENDRO_371 + DENDRO_476) - DENDRO_270*(0.5*DENDRO_396 + DENDRO_479) - DENDRO_270*(DENDRO_117*DENDRO_236 + DENDRO_350 + DENDRO_478) - DENDRO_274*(DENDRO_450 + DENDRO_484) - DENDRO_274*(DENDRO_211*DENDRO_43 + DENDRO_481) - DENDRO_283*(DENDRO_296 + DENDRO_487) - DENDRO_283*(0.5*DENDRO_455 + DENDRO_485) - DENDRO_283*(DENDRO_117*DENDRO_86 + DENDRO_454 + DENDRO_486) - DENDRO_529*(DENDRO_526 + DENDRO_528) + DENDRO_554)) + DENDRO_25*(-4*DENDRO_31 + DENDRO_321*DENDRO_54 + DENDRO_326*(-DENDRO_71 - DENDRO_73 + DENDRO_75 + DENDRO_81) + DENDRO_401*(DENDRO_62 + DENDRO_86) + alpha[pp]*(-DENDRO_114*DENDRO_223*DENDRO_288 + DENDRO_116*(-1.0*DENDRO_105 + DENDRO_113) + DENDRO_116*(-1.0*DENDRO_126 + DENDRO_132) + DENDRO_116*(DENDRO_130*DENDRO_299 + DENDRO_458) + DENDRO_116*(DENDRO_206*DENDRO_43 + DENDRO_460) + DENDRO_116*(DENDRO_213*DENDRO_48 + DENDRO_459) - DENDRO_125*DENDRO_280 - DENDRO_154*(DENDRO_456 + 1.0*DENDRO_457) - DENDRO_154*(-DENDRO_152*DENDRO_86 + DENDRO_454) + DENDRO_160*(DENDRO_463 + 1.0*DENDRO_464) + DENDRO_160*(DENDRO_137*DENDRO_86 + DENDRO_461) - DENDRO_167*(DENDRO_455 + DENDRO_457) - DENDRO_167*(DENDRO_130*DENDRO_47 + DENDRO_466) + DENDRO_170*(DENDRO_462 + DENDRO_464) + DENDRO_170*(DENDRO_125*DENDRO_47 + DENDRO_465) + DENDRO_205*DENDRO_261 + DENDRO_225*DENDRO_265 - DENDRO_241*DENDRO_281*DENDRO_314 + DENDRO_243*DENDRO_263 + DENDRO_268*DENDRO_372 - DENDRO_271*(DENDRO_119 - DENDRO_135) - DENDRO_275*(-DENDRO_272 + DENDRO_95) - DENDRO_285*DENDRO_467 - DENDRO_287*DENDRO_468 + DENDRO_315 - DENDRO_91*(DENDRO_85 + DENDRO_87*(DENDRO_331*DENDRO_49 + DENDRO_51*DENDRO_76 + DENDRO_65) + DENDRO_88*(DENDRO_341*DENDRO_44 + DENDRO_52) + DENDRO_90*(DENDRO_49*DENDRO_76 + DENDRO_51*DENDRO_78 + DENDRO_74)))) + DENDRO_29*(-4*DENDRO_317 + DENDRO_319*DENDRO_321 + DENDRO_327*(DENDRO_218 + DENDRO_325) + DENDRO_335*(-DENDRO_328 + DENDRO_329 - DENDRO_330 + DENDRO_334) + alpha[pp]*(DENDRO_102*DENDRO_370 + DENDRO_116*(DENDRO_353 + DENDRO_354) + DENDRO_116*(DENDRO_356 + 1.0*DENDRO_357) - DENDRO_154*(DENDRO_344 + DENDRO_345) - DENDRO_154*(DENDRO_350 + 1.0*DENDRO_351) - DENDRO_154*(DENDRO_218*DENDRO_348 + DENDRO_347) + DENDRO_160*(DENDRO_358 + DENDRO_360) + DENDRO_160*(DENDRO_117*DENDRO_230 + DENDRO_362) + DENDRO_160*(DENDRO_130*DENDRO_184 + DENDRO_363) + DENDRO_160*(DENDRO_206*DENDRO_365 + DENDRO_364) - DENDRO_167*(DENDRO_349 + DENDRO_351) - DENDRO_167*(DENDRO_102*DENDRO_206 + DENDRO_371) + DENDRO_174*DENDRO_368 - DENDRO_181*DENDRO_237*DENDRO_382 + DENDRO_181*DENDRO_369 - DENDRO_206*DENDRO_378 - DENDRO_211*DENDRO_377 - DENDRO_287*DENDRO_379 + DENDRO_333*DENDRO_373 + DENDRO_367*(DENDRO_355 + DENDRO_357) + DENDRO_367*(DENDRO_102*DENDRO_211 + DENDRO_366) - DENDRO_374*(DENDRO_117 + DENDRO_119) - DENDRO_375*(DENDRO_184 + DENDRO_186) - DENDRO_380*DENDRO_381 + DENDRO_398 - DENDRO_91*(DENDRO_336 + DENDRO_87*(DENDRO_10*DENDRO_338 + DENDRO_182*DENDRO_331 + DENDRO_340*DENDRO_76) + DENDRO_88*(DENDRO_183 + DENDRO_338*DENDRO_341 + DENDRO_340*DENDRO_9) + DENDRO_90*(DENDRO_182*DENDRO_76 + DENDRO_338*DENDRO_9 + DENDRO_340*DENDRO_78)))) - DENDRO_9*(DENDRO_717 + alpha[pp]*(DENDRO_115*(DENDRO_411 + DENDRO_583 + DENDRO_592) + DENDRO_115*(DENDRO_526 + DENDRO_527 + DENDRO_591) + DENDRO_116*(DENDRO_619 + DENDRO_677) + DENDRO_116*(DENDRO_622 + DENDRO_676) + DENDRO_116*(DENDRO_678 + DENDRO_679) - DENDRO_154*(DENDRO_132 + DENDRO_671) - DENDRO_154*(DENDRO_459 + DENDRO_672) - DENDRO_154*(DENDRO_582 + DENDRO_673) - DENDRO_154*(DENDRO_517 + DENDRO_581 + DENDRO_647) + DENDRO_160*(DENDRO_683 - DENDRO_89*DENDRO_97) + DENDRO_160*(DENDRO_684 + DENDRO_685) + DENDRO_160*(DENDRO_197*DENDRO_86 + DENDRO_681) + DENDRO_160*(DENDRO_200*DENDRO_43 + DENDRO_687) + DENDRO_170*(DENDRO_468 + DENDRO_89*DENDRO_92) - DENDRO_270*(DENDRO_362 + DENDRO_659) - DENDRO_274*(DENDRO_429 + DENDRO_664) - DENDRO_274*(DENDRO_660 + DENDRO_662) - DENDRO_274*(DENDRO_431 + DENDRO_665 + DENDRO_666) - DENDRO_283*(0.5*DENDRO_465 + DENDRO_668) - DENDRO_283*(DENDRO_43*DENDRO_89 + DENDRO_463 + DENDRO_667) - DENDRO_688*(DENDRO_359 + DENDRO_571 + DENDRO_604) + DENDRO_705 + DENDRO_709)) - DENDRO_9*(DENDRO_717 + alpha[pp]*(DENDRO_116*(DENDRO_415 + DENDRO_676) + DENDRO_116*(DENDRO_418 + DENDRO_653) + DENDRO_116*(DENDRO_418 + DENDRO_677) + DENDRO_116*(DENDRO_483 + DENDRO_656) + DENDRO_116*(DENDRO_574 + DENDRO_654) + DENDRO_116*(DENDRO_655 + DENDRO_678) - DENDRO_154*(DENDRO_460 + DENDRO_673) - DENDRO_154*(DENDRO_515 + DENDRO_707) - DENDRO_154*(DENDRO_518 + DENDRO_671) - DENDRO_154*(DENDRO_522 + DENDRO_672) + DENDRO_160*(DENDRO_681 + DENDRO_684) + DENDRO_160*(DENDRO_683 + DENDRO_718) + DENDRO_160*(DENDRO_239*DENDRO_49 + DENDRO_685) + DENDRO_160*(DENDRO_682 + DENDRO_718 + DENDRO_89*DENDRO_93) + DENDRO_160*(DENDRO_200*DENDRO_44 + DENDRO_686 + DENDRO_708) + DENDRO_170*(DENDRO_200*DENDRO_42 + DENDRO_452) - DENDRO_270*(DENDRO_363 + DENDRO_659) - DENDRO_270*(DENDRO_211*DENDRO_48 + DENDRO_605) - DENDRO_274*(0.5*DENDRO_438 + DENDRO_664) - DENDRO_274*(1.0*DENDRO_439 + DENDRO_665) - DENDRO_274*(DENDRO_200*DENDRO_93 + DENDRO_434 + DENDRO_662) - DENDRO_283*(DENDRO_461 + DENDRO_668) - DENDRO_283*(0.5*DENDRO_462 + DENDRO_667) - DENDRO_283*(DENDRO_300 + DENDRO_53*DENDRO_93 + DENDRO_706) - DENDRO_657*(DENDRO_112 + DENDRO_472 + DENDRO_474) - DENDRO_657*(DENDRO_469 + DENDRO_470 + DENDRO_629) - DENDRO_688*(DENDRO_389 + DENDRO_500 + DENDRO_585) + DENDRO_705));
                double DENDRO_720 = DENDRO_18*DENDRO_719;
                double DENDRO_721 = (1.0/12.0)*chi[pp];
                double DENDRO_722 = grad_1_beta0[pp];
                double DENDRO_723 = grad_1_beta2[pp];
                double DENDRO_724 = (1.0/3.0)*At1[pp];
                double DENDRO_725 = (2.0/3.0)*DENDRO_3;
                double DENDRO_726 = At4[pp]*DENDRO_10;
                double DENDRO_727 = -At3[pp]*DENDRO_13 + DENDRO_20 + DENDRO_726;
                double DENDRO_728 = -At1[pp]*DENDRO_25 + At3[pp]*DENDRO_9 - At4[pp]*DENDRO_23;
                double DENDRO_729 = -At1[pp]*DENDRO_23 + At3[pp]*DENDRO_10 - At4[pp]*DENDRO_29;
                double DENDRO_730 = 6.0*DENDRO_67;
                double DENDRO_731 = 6.0*DENDRO_82;
                double DENDRO_732 = 6.0*DENDRO_55;
                double DENDRO_733 = DENDRO_101*DENDRO_149;
                double DENDRO_734 = DENDRO_137*DENDRO_149;
                double DENDRO_735 = -DENDRO_208 + DENDRO_209 + DENDRO_210;
                double DENDRO_736 = DENDRO_47*DENDRO_735;
                double DENDRO_737 = DENDRO_734 + DENDRO_736;
                double DENDRO_738 = DENDRO_143*DENDRO_419;
                double DENDRO_739 = -DENDRO_134*DENDRO_423;
                double DENDRO_740 = DENDRO_193 - DENDRO_198 + DENDRO_199;
                double DENDRO_741 = DENDRO_138*DENDRO_623;
                double DENDRO_742 = -DENDRO_227 + DENDRO_228 + DENDRO_229;
                double DENDRO_743 = DENDRO_184*DENDRO_66;
                double DENDRO_744 = DENDRO_144*DENDRO_742 + DENDRO_743;
                double DENDRO_745 = DENDRO_149*DENDRO_295;
                double DENDRO_746 = DENDRO_157*DENDRO_299 + 0.25*DENDRO_50*DENDRO_735;
                double DENDRO_747 = DENDRO_145 + DENDRO_151;
                double DENDRO_748 = DENDRO_137*DENDRO_742;
                double DENDRO_749 = 1.0*DENDRO_115;
                double DENDRO_750 = DENDRO_137*DENDRO_143;
                double DENDRO_751 = DENDRO_42*DENDRO_735 + DENDRO_750;
                double DENDRO_752 = DENDRO_203*DENDRO_204;
                double DENDRO_753 = DENDRO_204*DENDRO_224;
                double DENDRO_754 = DENDRO_204*DENDRO_242;
                double DENDRO_755 = DENDRO_267*DENDRO_32;
                double DENDRO_756 = grad_2_beta0[pp];
                double DENDRO_757 = grad_2_beta1[pp];
                double DENDRO_758 = (1.0/3.0)*At2[pp];
                double DENDRO_759 = (2.0/3.0)*DENDRO_1;
                double DENDRO_760 = At2[pp]*DENDRO_9 - At4[pp]*DENDRO_13 + At5[pp]*DENDRO_10;
                double DENDRO_761 = -At2[pp]*DENDRO_25 + At4[pp]*DENDRO_9 - At5[pp]*DENDRO_23;
                double DENDRO_762 = -At5[pp]*DENDRO_29 + DENDRO_24 + DENDRO_726;
                double DENDRO_763 = DENDRO_107*DENDRO_143;
                double DENDRO_764 = DENDRO_149*DENDRO_346;
                double DENDRO_765 = DENDRO_215 + DENDRO_216 - DENDRO_217;
                double DENDRO_766 = 0.25*DENDRO_134*DENDRO_181;
                double DENDRO_767 = DENDRO_138*DENDRO_361;
                double DENDRO_768 = DENDRO_233 - DENDRO_234 + DENDRO_235;
                double DENDRO_769 = DENDRO_135*DENDRO_138;
                double DENDRO_770 = DENDRO_182*DENDRO_66;
                double DENDRO_771 = DENDRO_143*DENDRO_161;
                double DENDRO_772 = -DENDRO_766;
                double DENDRO_773 = DENDRO_143*DENDRO_346;
                double DENDRO_774 = DENDRO_102*DENDRO_143;
                double DENDRO_775 = (2.0/3.0)*DENDRO_0;
                double DENDRO_776 = 2*At4[pp];
                double DENDRO_777 = At3[pp]*DENDRO_26;
                double DENDRO_778 = DENDRO_18*DENDRO_776;
                double DENDRO_779 = DENDRO_32*DENDRO_400;
                double DENDRO_780 = DENDRO_18*DENDRO_83;
                double DENDRO_781 = DENDRO_172*DENDRO_740;
                double DENDRO_782 = 1.0*DENDRO_735;
                double DENDRO_783 = 0.25*DENDRO_750;
                double DENDRO_784 = DENDRO_134*DENDRO_361;
                double DENDRO_785 = DENDRO_740*DENDRO_92;
                double DENDRO_786 = (1.0/3.0)*At4[pp];
                double DENDRO_787 = DENDRO_135*DENDRO_742;
                double DENDRO_788 = -DENDRO_361*DENDRO_742 + DENDRO_610;
                double DENDRO_789 = DENDRO_624 - DENDRO_784;
                double DENDRO_790 = DENDRO_181*DENDRO_742;
                double DENDRO_791 = DENDRO_174*DENDRO_768;
                double DENDRO_792 = DENDRO_138*DENDRO_181;
                double DENDRO_793 = DENDRO_102*DENDRO_768;

                // Dendro: printing variables
                //--
                At_rhs00[pp] = (4.0/3.0)*At0[pp]*DENDRO_0 - DENDRO_1*DENDRO_2 - DENDRO_2*DENDRO_3 + DENDRO_4*DENDRO_5 + DENDRO_6*DENDRO_7 + DENDRO_721*(-12*DENDRO_31 + DENDRO_316*(-DENDRO_116*(DENDRO_105 - DENDRO_113) - DENDRO_116*(DENDRO_126 + DENDRO_133) - DENDRO_116*(DENDRO_136 + DENDRO_137*DENDRO_139) - DENDRO_116*(DENDRO_145 + DENDRO_149*DENDRO_43) - DENDRO_116*(DENDRO_143*DENDRO_48 + DENDRO_151) + DENDRO_134*DENDRO_280 + DENDRO_154*(DENDRO_163 + 1.0*DENDRO_164) - DENDRO_154*(-DENDRO_139*DENDRO_49 + DENDRO_152*DENDRO_66) + DENDRO_157*DENDRO_282*DENDRO_288 - DENDRO_160*(DENDRO_156 + 1.0*DENDRO_158) - DENDRO_160*(DENDRO_134*DENDRO_161 + 1.0*DENDRO_137*DENDRO_66) + DENDRO_167*(DENDRO_162 + DENDRO_164) + DENDRO_167*(DENDRO_165 + DENDRO_166) - DENDRO_170*(DENDRO_155 + DENDRO_158) - DENDRO_170*(DENDRO_168 + DENDRO_169) - DENDRO_203*DENDRO_205 - DENDRO_224*DENDRO_225 - DENDRO_242*DENDRO_243 - DENDRO_267*DENDRO_268 + DENDRO_271*(DENDRO_120 + DENDRO_135) + DENDRO_275*(DENDRO_272 + DENDRO_96) + DENDRO_281*DENDRO_283*DENDRO_66 + DENDRO_284*DENDRO_285 + DENDRO_286*DENDRO_287 + DENDRO_315 - DENDRO_91*(DENDRO_53*DENDRO_88 + DENDRO_85 + DENDRO_86*DENDRO_87 + DENDRO_89*DENDRO_90)) + DENDRO_54*DENDRO_57 - DENDRO_69*(-DENDRO_62 + DENDRO_66) + DENDRO_720*gt0[pp] - DENDRO_83*(DENDRO_71 + DENDRO_73 - DENDRO_75 - DENDRO_81)) - alpha[pp]*(-At0[pp]*K[pp] + DENDRO_19*(At0[pp]*DENDRO_9 - At1[pp]*DENDRO_13 + At2[pp]*DENDRO_10) + DENDRO_27*(-At0[pp]*DENDRO_25 + DENDRO_20 + DENDRO_24) + DENDRO_30*(-At0[pp]*DENDRO_23 + At1[pp]*DENDRO_10 - At2[pp]*DENDRO_29)) + beta0[pp]*agrad_0_At0[pp] + beta1[pp]*agrad_1_At0[pp] + beta2[pp]*agrad_2_At0[pp];
                //--
                At_rhs01[pp] = At0[pp]*DENDRO_722 - At1[pp]*DENDRO_725 + At2[pp]*DENDRO_723 + At3[pp]*DENDRO_4 + At4[pp]*DENDRO_6 + DENDRO_0*DENDRO_724 + DENDRO_1*DENDRO_724 + DENDRO_721*(-DENDRO_18*DENDRO_730*(DENDRO_134 - DENDRO_716) + DENDRO_316*(DENDRO_116*(-DENDRO_299*DENDRO_740 + DENDRO_679) - DENDRO_116*(-DENDRO_622 + DENDRO_674 + DENDRO_675) + DENDRO_116*(-DENDRO_143*DENDRO_623 + DENDRO_149*DENDRO_419 + DENDRO_652) + DENDRO_154*(DENDRO_741 + DENDRO_744) + DENDRO_154*(DENDRO_745 + DENDRO_746) + DENDRO_154*(DENDRO_157*DENDRO_417 + DENDRO_747) + DENDRO_154*(DENDRO_133 + DENDRO_669 + DENDRO_670) + DENDRO_160*(-DENDRO_43*DENDRO_740 + DENDRO_687) - DENDRO_160*(DENDRO_134*DENDRO_504 + DENDRO_134*DENDRO_519 + DENDRO_196*DENDRO_66) + DENDRO_160*(-DENDRO_134*DENDRO_623 - DENDRO_197*DENDRO_66 + DENDRO_680) + DENDRO_160*(-DENDRO_143*DENDRO_295 + DENDRO_157*DENDRO_97 + DENDRO_682) - DENDRO_170*(DENDRO_157*DENDRO_92 + DENDRO_286) + DENDRO_269*(DENDRO_733 + DENDRO_737) + DENDRO_270*(-DENDRO_362 + DENDRO_509 + DENDRO_608) + DENDRO_274*(-DENDRO_660 + DENDRO_661) - DENDRO_274*(-DENDRO_134*DENDRO_391 + DENDRO_663 + DENDRO_739) - DENDRO_274*(-DENDRO_143*DENDRO_272 + DENDRO_666 + DENDRO_738) + DENDRO_283*(DENDRO_156 + DENDRO_157*DENDRO_43 + DENDRO_157*DENDRO_44) + DENDRO_283*(0.25*DENDRO_168 + 0.5*DENDRO_169 + DENDRO_417*DENDRO_66) - DENDRO_40*DENDRO_753 + DENDRO_690 + DENDRO_692 + DENDRO_693 + DENDRO_694 + DENDRO_695 + DENDRO_696 + DENDRO_697 + DENDRO_698 + DENDRO_699 + DENDRO_700 + DENDRO_701 + DENDRO_702 + DENDRO_703 - DENDRO_704*DENDRO_755 + DENDRO_709 - DENDRO_749*(DENDRO_149*DENDRO_92 + DENDRO_751) - DENDRO_749*(DENDRO_134*DENDRO_174 + DENDRO_138*DENDRO_172 + DENDRO_748) - DENDRO_752*DENDRO_94 - DENDRO_754*DENDRO_99 - DENDRO_91*(DENDRO_111*DENDRO_532 + DENDRO_125*DENDRO_531 + DENDRO_213*DENDRO_533 + DENDRO_689)) + DENDRO_704*DENDRO_719 - 12*DENDRO_710 - DENDRO_731*(-DENDRO_712 + DENDRO_713 + DENDRO_714 - DENDRO_715) + DENDRO_732*(DENDRO_32*(DENDRO_37*DENDRO_704 + DENDRO_559) + DENDRO_711)) - alpha[pp]*(-At1[pp]*K[pp] + DENDRO_19*DENDRO_727 + DENDRO_27*DENDRO_728 + DENDRO_30*DENDRO_729) + beta0[pp]*agrad_0_At1[pp] + beta1[pp]*agrad_1_At1[pp] + beta2[pp]*agrad_2_At1[pp];
                //--
                At_rhs02[pp] = At0[pp]*DENDRO_756 + At1[pp]*DENDRO_757 - At2[pp]*DENDRO_759 + At4[pp]*DENDRO_4 + At5[pp]*DENDRO_6 + DENDRO_0*DENDRO_758 + DENDRO_3*DENDRO_758 + DENDRO_721*(DENDRO_316*(-DENDRO_100*DENDRO_752 - DENDRO_116*(DENDRO_497 - DENDRO_498 + DENDRO_576) + DENDRO_116*(-DENDRO_149*DENDRO_504 - DENDRO_43*DENDRO_765 + DENDRO_773) - DENDRO_116*(DENDRO_299*DENDRO_768 + DENDRO_766 + DENDRO_767) + DENDRO_116*(DENDRO_346*DENDRO_742 - DENDRO_767 + DENDRO_772) - DENDRO_118*DENDRO_754 - DENDRO_154*(DENDRO_121*DENDRO_139 - DENDRO_769 - DENDRO_770) - DENDRO_154*(DENDRO_121*DENDRO_157 - DENDRO_149*DENDRO_161 - DENDRO_51*DENDRO_765) + DENDRO_154*(DENDRO_48*DENDRO_768 + DENDRO_769 + DENDRO_770) - DENDRO_160*(DENDRO_136 + DENDRO_744) - DENDRO_160*(DENDRO_746 + DENDRO_771) - DENDRO_160*(DENDRO_157*DENDRO_365 + DENDRO_747) - DENDRO_160*(DENDRO_136 + DENDRO_138*DENDRO_519 + DENDRO_743) + DENDRO_167*(DENDRO_102*DENDRO_157 + DENDRO_284) + DENDRO_270*(DENDRO_477 - DENDRO_573) - DENDRO_270*(-DENDRO_135*DENDRO_149 - DENDRO_48*DENDRO_765 + DENDRO_764) + DENDRO_273*(DENDRO_751 + DENDRO_763) + DENDRO_274*(DENDRO_134*DENDRO_184 + 0.25*DENDRO_748) + DENDRO_283*(0.25*DENDRO_165 + 1.0*DENDRO_166) + DENDRO_283*(DENDRO_157*DENDRO_48 + DENDRO_157*DENDRO_49 + DENDRO_163) - DENDRO_45*DENDRO_753 + DENDRO_535 + DENDRO_537 + DENDRO_539 + DENDRO_540 + DENDRO_542 + DENDRO_544 + DENDRO_546 + DENDRO_547 + DENDRO_548 + DENDRO_549 + DENDRO_550 + DENDRO_551 + DENDRO_552 - DENDRO_553*DENDRO_755 + DENDRO_590 - DENDRO_749*(DENDRO_737 + DENDRO_774) - DENDRO_91*(DENDRO_104*DENDRO_532 + DENDRO_130*DENDRO_531 + DENDRO_206*DENDRO_533 + DENDRO_530)) + DENDRO_553*DENDRO_719 - 12*DENDRO_555 + DENDRO_569*DENDRO_732 - DENDRO_730*(DENDRO_556 + DENDRO_557 - DENDRO_558 - DENDRO_560) - DENDRO_731*(DENDRO_562 + DENDRO_563 - DENDRO_564 - DENDRO_566)) - alpha[pp]*(-At2[pp]*K[pp] + DENDRO_19*DENDRO_760 + DENDRO_27*DENDRO_761 + DENDRO_30*DENDRO_762) + beta0[pp]*agrad_0_At2[pp] + beta1[pp]*agrad_1_At2[pp] + beta2[pp]*agrad_2_At2[pp];
                //--
                At_rhs11[pp] = (4.0/3.0)*At3[pp]*DENDRO_1 - At3[pp]*DENDRO_725 - At3[pp]*DENDRO_775 + DENDRO_5*DENDRO_722 + DENDRO_721*(DENDRO_316*(-DENDRO_116*(DENDRO_420 - 1.0*DENDRO_422) + DENDRO_116*(DENDRO_427 - 1.0*DENDRO_781) - DENDRO_116*(DENDRO_423*DENDRO_742 + DENDRO_425) + DENDRO_134*DENDRO_445 + DENDRO_143*DENDRO_444 + DENDRO_154*(DENDRO_413 - DENDRO_415) + DENDRO_154*(DENDRO_143*DENDRO_417 + DENDRO_150*DENDRO_735) + DENDRO_154*(DENDRO_417*DENDRO_742 + DENDRO_784) + DENDRO_154*(DENDRO_44*DENDRO_782 + DENDRO_783) + DENDRO_160*(DENDRO_430 + DENDRO_739) + DENDRO_160*(DENDRO_432 + DENDRO_738) + DENDRO_160*(DENDRO_434 - 1.0*DENDRO_785) + DENDRO_170*(DENDRO_433 - DENDRO_785) + DENDRO_170*(-DENDRO_134*DENDRO_172 + DENDRO_438) + DENDRO_170*(-DENDRO_143*DENDRO_92 + DENDRO_439) - DENDRO_172*DENDRO_754 + 6.0*DENDRO_191*DENDRO_273*DENDRO_740 - DENDRO_191*DENDRO_752 + DENDRO_367*(DENDRO_426 - DENDRO_781) + DENDRO_367*(DENDRO_436 - DENDRO_735*DENDRO_92) + DENDRO_367*(-DENDRO_172*DENDRO_742 + DENDRO_437) - DENDRO_404*DENDRO_755 + DENDRO_440*(DENDRO_187 + DENDRO_361) + DENDRO_442*DENDRO_735 + DENDRO_443*(DENDRO_408 + DENDRO_96) + DENDRO_453 - DENDRO_753*DENDRO_92 - DENDRO_91*(DENDRO_200*DENDRO_88 + DENDRO_221*DENDRO_90 + DENDRO_239*DENDRO_87 + DENDRO_407)) - 12*DENDRO_399 + DENDRO_56*(DENDRO_32*(DENDRO_37*DENDRO_405 + DENDRO_402) + DENDRO_406) + DENDRO_69*(DENDRO_239 - DENDRO_779*(-DENDRO_58 + DENDRO_59 + DENDRO_60)) + DENDRO_720*gt3[pp] + DENDRO_780*(DENDRO_221 - DENDRO_779*(DENDRO_322 + DENDRO_323 - DENDRO_77))) + DENDRO_723*DENDRO_776 - alpha[pp]*(-At3[pp]*K[pp] + DENDRO_19*DENDRO_728 + DENDRO_727*DENDRO_777 + DENDRO_729*DENDRO_778) + beta0[pp]*agrad_0_At3[pp] + beta1[pp]*agrad_1_At3[pp] + beta2[pp]*agrad_2_At3[pp];
                //--
                At_rhs12[pp] = At1[pp]*DENDRO_756 + At2[pp]*DENDRO_722 + At3[pp]*DENDRO_757 - At4[pp]*DENDRO_775 + At5[pp]*DENDRO_723 + DENDRO_1*DENDRO_786 + DENDRO_3*DENDRO_786 + DENDRO_721*(-DENDRO_18*DENDRO_731*(-DENDRO_642 + DENDRO_735) + DENDRO_316*(DENDRO_116*(DENDRO_188*DENDRO_740 + DENDRO_650) + DENDRO_116*(-DENDRO_196*DENDRO_768 + DENDRO_788) + DENDRO_116*(DENDRO_342*DENDRO_742 + DENDRO_788) - DENDRO_116*(DENDRO_504*DENDRO_735 + DENDRO_649 + DENDRO_765*DENDRO_93) + DENDRO_116*(-DENDRO_519*DENDRO_735 + DENDRO_612 - DENDRO_623*DENDRO_735) + DENDRO_153*(DENDRO_733 + DENDRO_736 + DENDRO_774) - DENDRO_154*(DENDRO_139*DENDRO_188 + DENDRO_772 - DENDRO_787) - DENDRO_154*(-DENDRO_161*DENDRO_735 - DENDRO_44*DENDRO_765 + DENDRO_773) + DENDRO_154*(DENDRO_417*DENDRO_768 + DENDRO_766 + DENDRO_787) + DENDRO_160*(-DENDRO_138*DENDRO_391 + DENDRO_789) + DENDRO_160*(-DENDRO_365*DENDRO_740 + DENDRO_655) + DENDRO_160*(-DENDRO_623*DENDRO_742 + DENDRO_789) + DENDRO_160*(-DENDRO_149*DENDRO_272 + DENDRO_652 - 0.25*DENDRO_763) + DENDRO_160*(-DENDRO_295*DENDRO_735 + DENDRO_620 - DENDRO_783) + DENDRO_160*(-DENDRO_417*DENDRO_740 + DENDRO_482 + DENDRO_626) - DENDRO_185*DENDRO_754 - DENDRO_194*DENDRO_752 + DENDRO_270*(DENDRO_594 - DENDRO_644) - DENDRO_270*(-DENDRO_135*DENDRO_735 + DENDRO_346*DENDRO_735 - DENDRO_417*DENDRO_765) - DENDRO_274*(-DENDRO_391*DENDRO_742 + DENDRO_645) - DENDRO_274*(-DENDRO_196*DENDRO_740 - DENDRO_197*DENDRO_740 + DENDRO_427) - DENDRO_274*(-DENDRO_272*DENDRO_735 + DENDRO_599 + DENDRO_646) + DENDRO_283*(DENDRO_117*DENDRO_134 + DENDRO_741) + DENDRO_283*(DENDRO_145 + DENDRO_745 + DENDRO_771) + DENDRO_367*(-DENDRO_174*DENDRO_740 + DENDRO_451) - DENDRO_632*DENDRO_755 + DENDRO_633 + DENDRO_658 - DENDRO_753*DENDRO_98 - DENDRO_91*(DENDRO_178*DENDRO_532 + DENDRO_211*DENDRO_533 + DENDRO_230*DENDRO_531 + DENDRO_631)) + DENDRO_632*DENDRO_719 - 12*DENDRO_635 - DENDRO_730*(-DENDRO_636 + DENDRO_637 + DENDRO_638 - DENDRO_640) + DENDRO_732*(DENDRO_32*(DENDRO_37*DENDRO_632 + DENDRO_565) + DENDRO_641)) - alpha[pp]*(-At4[pp]*K[pp] + DENDRO_19*DENDRO_761 + DENDRO_760*DENDRO_777 + DENDRO_762*DENDRO_778) + beta0[pp]*agrad_0_At4[pp] + beta1[pp]*agrad_1_At4[pp] + beta2[pp]*agrad_2_At4[pp];
                //--
                At_rhs22[pp] = (4.0/3.0)*At5[pp]*DENDRO_3 - At5[pp]*DENDRO_759 - At5[pp]*DENDRO_775 + DENDRO_7*DENDRO_756 + DENDRO_721*(DENDRO_316*(DENDRO_102*DENDRO_138*DENDRO_381 - DENDRO_102*DENDRO_753 - DENDRO_116*(DENDRO_352 - 1.0*DENDRO_354) - DENDRO_116*(0.25*DENDRO_790 + 1.0*DENDRO_791) + DENDRO_149*DENDRO_378 + DENDRO_154*(DENDRO_343 - 1.0*DENDRO_345) + DENDRO_154*(0.25*DENDRO_792 + 1.0*DENDRO_793) - DENDRO_154*(-DENDRO_348*DENDRO_765 + DENDRO_764) - DENDRO_160*(DENDRO_117*DENDRO_742 + DENDRO_767) - DENDRO_160*(DENDRO_138*DENDRO_184 + DENDRO_787) - DENDRO_160*(DENDRO_144*DENDRO_735 + DENDRO_149*DENDRO_365) - DENDRO_160*(DENDRO_49*DENDRO_782 + 0.25*DENDRO_734) + DENDRO_167*(DENDRO_792 + DENDRO_793) + DENDRO_167*(DENDRO_102*DENDRO_149 + DENDRO_47*DENDRO_765) + DENDRO_174*DENDRO_287*DENDRO_742 - DENDRO_174*DENDRO_752 + 6.0*DENDRO_181*DENDRO_269*DENDRO_768 - DENDRO_181*DENDRO_754 - DENDRO_333*DENDRO_755 - DENDRO_367*(DENDRO_790 + DENDRO_791) - DENDRO_367*(DENDRO_102*DENDRO_735 + DENDRO_137*DENDRO_765) + DENDRO_374*(DENDRO_120 + DENDRO_339) + DENDRO_375*(DENDRO_187 + DENDRO_337) + DENDRO_377*DENDRO_735 + DENDRO_398 - DENDRO_91*(DENDRO_189*DENDRO_88 + DENDRO_218*DENDRO_90 + DENDRO_236*DENDRO_87 + DENDRO_336)) - 12*DENDRO_317 + DENDRO_319*DENDRO_57 - DENDRO_68*(DENDRO_328 - DENDRO_329 + DENDRO_330 - DENDRO_334) + DENDRO_720*gt5[pp] - DENDRO_780*(-DENDRO_325 + DENDRO_765)) + DENDRO_757*DENDRO_776 - alpha[pp]*(At5[pp]*DENDRO_26*DENDRO_762 - At5[pp]*K[pp] + DENDRO_30*DENDRO_761 + DENDRO_760*DENDRO_778) + beta0[pp]*agrad_0_At5[pp] + beta1[pp]*agrad_1_At5[pp] + beta2[pp]*agrad_2_At5[pp];
                // Dendro: reduced ops: 3279
                // Dendro: }}} 
//[[[end]]]
            }
        }
    }
    bssn::timer::t_rhs_At.stop();

    bssn::timer::t_rhs_K.start();
    //K_rhs
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.K_rhs]
                vnames = ['K_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 3960 
                // Dendro: printing temp variables
                double DENDRO_0 = gt0[pp]*gt3[pp];
                double DENDRO_1 = pow(gt1[pp], 2);
                double DENDRO_2 = DENDRO_0 - DENDRO_1;
                double DENDRO_3 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
                double DENDRO_4 = 0.5*grad_2_gt5[pp];
                double DENDRO_5 = pow(gt2[pp], 2);
                double DENDRO_6 = gt0[pp]*gt5[pp];
                double DENDRO_7 = DENDRO_5 - DENDRO_6;
                double DENDRO_8 = grad_1_gt5[pp];
                double DENDRO_9 = -0.5*DENDRO_8 + 1.0*grad_2_gt4[pp];
                double DENDRO_10 = gt2[pp]*gt4[pp];
                double DENDRO_11 = -DENDRO_10 + gt1[pp]*gt5[pp];
                double DENDRO_12 = grad_0_gt5[pp];
                double DENDRO_13 = -0.5*DENDRO_12 + 1.0*grad_2_gt2[pp];
                double DENDRO_14 = 1.0/chi[pp];
                double DENDRO_15 = grad_2_chi[pp];
                double DENDRO_16 = grad_1_chi[pp];
                double DENDRO_17 = grad_0_chi[pp];
                double DENDRO_18 = DENDRO_11*DENDRO_17 + DENDRO_15*DENDRO_3 + DENDRO_16*DENDRO_7;
                double DENDRO_19 = DENDRO_14*DENDRO_18;
                double DENDRO_20 = 0.5*gt5[pp];
                double DENDRO_21 = grad_1_alpha[pp];
                double DENDRO_22 = pow(gt4[pp], 2);
                double DENDRO_23 = gt3[pp]*gt5[pp];
                double DENDRO_24 = DENDRO_1*gt5[pp] - 2*DENDRO_10*gt1[pp] + DENDRO_22*gt0[pp] - DENDRO_23*gt0[pp] + DENDRO_5*gt3[pp];
                double DENDRO_25 = 1.0/DENDRO_24;
                double DENDRO_26 = DENDRO_21*DENDRO_25;
                double DENDRO_27 = gt2[pp]*gt3[pp];
                double DENDRO_28 = gt1[pp]*gt4[pp];
                double DENDRO_29 = DENDRO_27 - DENDRO_28;
                double DENDRO_30 = DENDRO_22 - DENDRO_23;
                double DENDRO_31 = DENDRO_11*DENDRO_16 + DENDRO_15*DENDRO_29 + DENDRO_17*DENDRO_30;
                double DENDRO_32 = DENDRO_14*DENDRO_31;
                double DENDRO_33 = grad_0_alpha[pp];
                double DENDRO_34 = DENDRO_25*DENDRO_33;
                double DENDRO_35 = grad_2_alpha[pp];
                double DENDRO_36 = -DENDRO_0 + DENDRO_1;
                double DENDRO_37 = DENDRO_25*DENDRO_36;
                double DENDRO_38 = DENDRO_25*DENDRO_3;
                double DENDRO_39 = DENDRO_25*DENDRO_29;
                double DENDRO_40 = DENDRO_15*DENDRO_36 + DENDRO_16*DENDRO_3 + DENDRO_17*DENDRO_29;
                double DENDRO_41 = DENDRO_25*chi[pp];
                double DENDRO_42 = -DENDRO_5 + DENDRO_6;
                double DENDRO_43 = 0.5*grad_1_gt3[pp];
                double DENDRO_44 = grad_2_gt3[pp];
                double DENDRO_45 = -0.5*DENDRO_44 + 1.0*grad_1_gt4[pp];
                double DENDRO_46 = grad_0_gt3[pp];
                double DENDRO_47 = -0.5*DENDRO_46 + 1.0*grad_1_gt1[pp];
                double DENDRO_48 = DENDRO_14*DENDRO_40;
                double DENDRO_49 = 0.5*gt3[pp];
                double DENDRO_50 = DENDRO_25*DENDRO_35;
                double DENDRO_51 = DENDRO_25*DENDRO_7;
                double DENDRO_52 = DENDRO_11*DENDRO_25;
                double DENDRO_53 = -DENDRO_22 + DENDRO_23;
                double DENDRO_54 = 0.5*grad_0_gt0[pp];
                double DENDRO_55 = grad_2_gt0[pp];
                double DENDRO_56 = -0.5*DENDRO_55 + 1.0*grad_0_gt2[pp];
                double DENDRO_57 = grad_1_gt0[pp];
                double DENDRO_58 = -0.5*DENDRO_57 + 1.0*grad_0_gt1[pp];
                double DENDRO_59 = 0.5*gt0[pp];
                double DENDRO_60 = DENDRO_25*DENDRO_30;
                double DENDRO_61 = grad_1_gt2[pp];
                double DENDRO_62 = grad_2_gt1[pp];
                double DENDRO_63 = grad_0_gt4[pp];
                double DENDRO_64 = DENDRO_61 + DENDRO_62 - DENDRO_63;
                double DENDRO_65 = DENDRO_25*gt4[pp];
                double DENDRO_66 = 0.5*DENDRO_35;
                double DENDRO_67 = 0.5*DENDRO_21;
                double DENDRO_68 = 2*chi[pp];
                double DENDRO_69 = -DENDRO_61 + DENDRO_62 + DENDRO_63;
                double DENDRO_70 = DENDRO_25*gt2[pp];
                double DENDRO_71 = 0.5*DENDRO_33;
                double DENDRO_72 = -DENDRO_27 + DENDRO_28;
                double DENDRO_73 = 2*DENDRO_72;
                double DENDRO_74 = DENDRO_61 - DENDRO_62 + DENDRO_63;
                double DENDRO_75 = DENDRO_25*gt1[pp];
                double DENDRO_76 = pow(DENDRO_11, 2);
                double DENDRO_77 = pow(DENDRO_72, 2);
                double DENDRO_78 = At4[pp]*DENDRO_72;
                double DENDRO_79 = 2*DENDRO_11;
                double DENDRO_80 = At1[pp]*DENDRO_53;
                double DENDRO_81 = At2[pp]*DENDRO_53;
                double DENDRO_82 = pow(DENDRO_24, -2);
                double DENDRO_83 = 3*DENDRO_82;
                double DENDRO_84 = pow(DENDRO_3, 2);
                double DENDRO_85 = DENDRO_11*DENDRO_3;
                double DENDRO_86 = At1[pp]*DENDRO_42;
                double DENDRO_87 = 2*DENDRO_3;
                double DENDRO_88 = DENDRO_3*DENDRO_72;
                double DENDRO_89 = At2[pp]*DENDRO_2;
                double DENDRO_90 = At4[pp]*DENDRO_2;
                double DENDRO_91 = At0[pp]*DENDRO_53;
                double DENDRO_92 = DENDRO_11*DENDRO_72;
                double DENDRO_93 = At5[pp]*DENDRO_72;
                double DENDRO_94 = 6*DENDRO_82;
                double DENDRO_95 = At3[pp]*DENDRO_42;

                // Dendro: printing variables
                //--
                K_rhs[pp] = -DENDRO_2*DENDRO_41*(DENDRO_26*(DENDRO_11*DENDRO_13 + DENDRO_19*DENDRO_20 + DENDRO_3*DENDRO_4 + DENDRO_7*DENDRO_9) + DENDRO_34*(DENDRO_11*DENDRO_9 + DENDRO_13*DENDRO_30 + DENDRO_20*DENDRO_32 + DENDRO_29*DENDRO_4) + DENDRO_35*(DENDRO_13*DENDRO_39 - DENDRO_14*(1.0*DENDRO_15 - DENDRO_20*DENDRO_25*DENDRO_40) + DENDRO_37*DENDRO_4 + DENDRO_38*DENDRO_9) - grad2_2_2_alpha[pp]) + DENDRO_38*DENDRO_68*(0.5*DENDRO_34*(DENDRO_11*DENDRO_44 + DENDRO_29*DENDRO_8 + DENDRO_30*DENDRO_64 + DENDRO_32*gt4[pp]) + DENDRO_66*(-DENDRO_14*(DENDRO_16 - DENDRO_40*DENDRO_65) + DENDRO_37*DENDRO_8 + DENDRO_38*DENDRO_44 + DENDRO_39*DENDRO_64) + DENDRO_67*(-DENDRO_14*(DENDRO_15 - DENDRO_18*DENDRO_65) + DENDRO_38*DENDRO_8 + DENDRO_44*DENDRO_51 + DENDRO_52*DENDRO_64) - grad2_1_2_alpha[pp]) - DENDRO_41*DENDRO_42*(DENDRO_21*(-DENDRO_14*(1.0*DENDRO_16 - DENDRO_18*DENDRO_25*DENDRO_49) + DENDRO_38*DENDRO_45 + DENDRO_43*DENDRO_51 + DENDRO_47*DENDRO_52) + DENDRO_34*(DENDRO_11*DENDRO_43 + DENDRO_29*DENDRO_45 + DENDRO_30*DENDRO_47 + DENDRO_32*DENDRO_49) + DENDRO_50*(DENDRO_29*DENDRO_47 + DENDRO_3*DENDRO_43 + DENDRO_36*DENDRO_45 + DENDRO_48*DENDRO_49) - grad2_1_1_alpha[pp]) - DENDRO_41*DENDRO_53*(DENDRO_26*(DENDRO_11*DENDRO_54 + DENDRO_19*DENDRO_59 + DENDRO_3*DENDRO_56 + DENDRO_58*DENDRO_7) + DENDRO_33*(-DENDRO_14*(1.0*DENDRO_17 - DENDRO_25*DENDRO_31*DENDRO_59) + DENDRO_39*DENDRO_56 + DENDRO_52*DENDRO_58 + DENDRO_54*DENDRO_60) + DENDRO_50*(DENDRO_29*DENDRO_54 + DENDRO_3*DENDRO_58 + DENDRO_36*DENDRO_56 + DENDRO_48*DENDRO_59) - grad2_0_0_alpha[pp]) - DENDRO_41*DENDRO_73*(0.5*DENDRO_26*(DENDRO_11*DENDRO_55 + DENDRO_12*DENDRO_3 + DENDRO_19*gt2[pp] + DENDRO_69*DENDRO_7) + DENDRO_66*(DENDRO_12*DENDRO_37 - DENDRO_14*(DENDRO_17 - DENDRO_40*DENDRO_70) + DENDRO_38*DENDRO_69 + DENDRO_39*DENDRO_55) + DENDRO_71*(DENDRO_12*DENDRO_39 - DENDRO_14*(DENDRO_15 - DENDRO_31*DENDRO_70) + DENDRO_52*DENDRO_69 + DENDRO_55*DENDRO_60) - grad2_0_2_alpha[pp]) + DENDRO_52*DENDRO_68*(0.5*DENDRO_50*(DENDRO_29*DENDRO_57 + DENDRO_3*DENDRO_46 + DENDRO_36*DENDRO_74 + DENDRO_48*gt1[pp]) + DENDRO_67*(-DENDRO_14*(DENDRO_17 - DENDRO_18*DENDRO_75) + DENDRO_38*DENDRO_74 + DENDRO_46*DENDRO_51 + DENDRO_52*DENDRO_57) + DENDRO_71*(-DENDRO_14*(DENDRO_16 - DENDRO_31*DENDRO_75) + DENDRO_39*DENDRO_74 + DENDRO_46*DENDRO_52 + DENDRO_57*DENDRO_60) - grad2_0_1_alpha[pp]) + (1.0/3.0)*alpha[pp]*(At0[pp]*DENDRO_83*(At0[pp]*pow(DENDRO_53, 2) + At3[pp]*DENDRO_76 + At5[pp]*DENDRO_77 + DENDRO_73*DENDRO_81 - DENDRO_78*DENDRO_79 - DENDRO_79*DENDRO_80) + At1[pp]*DENDRO_94*(At1[pp]*DENDRO_76 - At2[pp]*DENDRO_92 + At4[pp]*DENDRO_85 - DENDRO_11*DENDRO_91 - DENDRO_11*DENDRO_95 - DENDRO_3*DENDRO_81 - DENDRO_3*DENDRO_93 + DENDRO_42*DENDRO_78 + DENDRO_42*DENDRO_80) + At2[pp]*DENDRO_94*(-At1[pp]*DENDRO_92 + At2[pp]*DENDRO_77 + At3[pp]*DENDRO_85 - DENDRO_11*DENDRO_90 + DENDRO_2*DENDRO_81 + DENDRO_2*DENDRO_93 - DENDRO_3*DENDRO_78 - DENDRO_3*DENDRO_80 + DENDRO_72*DENDRO_91) + At3[pp]*DENDRO_83*(At0[pp]*DENDRO_76 + 2*At2[pp]*DENDRO_85 + At3[pp]*pow(DENDRO_42, 2) - At4[pp]*DENDRO_42*DENDRO_87 + At5[pp]*DENDRO_84 - DENDRO_79*DENDRO_86) + At4[pp]*DENDRO_94*(-At0[pp]*DENDRO_92 + At1[pp]*DENDRO_85 - At2[pp]*DENDRO_88 + At4[pp]*DENDRO_84 - At5[pp]*DENDRO_2*DENDRO_3 - DENDRO_11*DENDRO_89 - DENDRO_3*DENDRO_95 + DENDRO_42*DENDRO_90 + DENDRO_72*DENDRO_86) + At5[pp]*DENDRO_83*(At0[pp]*DENDRO_77 - 2*At1[pp]*DENDRO_88 + At3[pp]*DENDRO_84 + At5[pp]*pow(DENDRO_2, 2) + DENDRO_73*DENDRO_89 - DENDRO_87*DENDRO_90) + pow(K[pp], 2)) + beta0[pp]*agrad_0_K[pp] + beta1[pp]*agrad_1_K[pp] + beta2[pp]*agrad_2_K[pp];
                // Dendro: reduced ops: 479
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }
    bssn::timer::t_rhs_K.stop();

    bssn::timer::t_rhs_Gt.start();

    //#include "CalGt.cpp"
    //CalGt
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.CalGt]
                vnames = ['CalGt']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 1716 
                // Dendro: printing temp variables
                double DENDRO_0 = pow(gt4[pp], 2);
                double DENDRO_1 = pow(gt1[pp], 2);
                double DENDRO_2 = pow(gt2[pp], 2);
                double DENDRO_3 = gt3[pp]*gt5[pp];
                double DENDRO_4 = gt1[pp]*gt2[pp];
                double DENDRO_5 = pow(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt0[pp] - 2*DENDRO_4*gt4[pp], -2);
                double DENDRO_6 = -DENDRO_4 + gt0[pp]*gt4[pp];
                double DENDRO_7 = grad_2_gt3[pp];
                double DENDRO_8 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
                double DENDRO_9 = grad_1_gt5[pp];
                double DENDRO_10 = gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp];
                double DENDRO_11 = -DENDRO_0 + DENDRO_3;
                double DENDRO_12 = grad_1_gt2[pp];
                double DENDRO_13 = grad_2_gt1[pp];
                double DENDRO_14 = grad_0_gt4[pp];
                double DENDRO_15 = DENDRO_12 + DENDRO_13 - DENDRO_14;
                double DENDRO_16 = grad_0_gt3[pp];
                double DENDRO_17 = grad_1_gt0[pp];
                double DENDRO_18 = DENDRO_12 - DENDRO_13 + DENDRO_14;
                double DENDRO_19 = grad_0_gt5[pp];
                double DENDRO_20 = grad_2_gt0[pp];
                double DENDRO_21 = -DENDRO_12 + DENDRO_13 + DENDRO_14;
                double DENDRO_22 = 1.0*DENDRO_10;
                double DENDRO_23 = -DENDRO_1 + gt0[pp]*gt3[pp];
                double DENDRO_24 = 0.5*grad_2_gt5[pp];
                double DENDRO_25 = 0.5*DENDRO_9 - 1.0*grad_2_gt4[pp];
                double DENDRO_26 = 0.5*DENDRO_19 - 1.0*grad_2_gt2[pp];
                double DENDRO_27 = -DENDRO_2 + gt0[pp]*gt5[pp];
                double DENDRO_28 = 0.5*grad_1_gt3[pp];
                double DENDRO_29 = -0.5*DENDRO_7 + 1.0*grad_1_gt4[pp];
                double DENDRO_30 = 0.5*DENDRO_16 - 1.0*grad_1_gt1[pp];
                double DENDRO_31 = 0.5*grad_0_gt0[pp];
                double DENDRO_32 = -0.5*DENDRO_17 + 1.0*grad_0_gt1[pp];
                double DENDRO_33 = -0.5*DENDRO_20 + 1.0*grad_0_gt2[pp];

                // Dendro: printing variables
                //--
                CalGt0[pp] = DENDRO_5*(-DENDRO_11*(-DENDRO_10*DENDRO_33 - DENDRO_11*DENDRO_31 + DENDRO_32*DENDRO_8) - DENDRO_22*(-DENDRO_10*DENDRO_19 - DENDRO_11*DENDRO_20 + DENDRO_21*DENDRO_8) - DENDRO_23*(-DENDRO_10*DENDRO_24 + DENDRO_11*DENDRO_26 - DENDRO_25*DENDRO_8) - DENDRO_27*(-DENDRO_10*DENDRO_29 + DENDRO_11*DENDRO_30 + DENDRO_28*DENDRO_8) + DENDRO_6*(-DENDRO_10*DENDRO_9 - DENDRO_11*DENDRO_15 + DENDRO_7*DENDRO_8) + DENDRO_8*(-DENDRO_10*DENDRO_18 - DENDRO_11*DENDRO_17 + DENDRO_16*DENDRO_8));
                //--
                CalGt1[pp] = DENDRO_5*(-DENDRO_11*(-DENDRO_27*DENDRO_32 + DENDRO_31*DENDRO_8 + DENDRO_33*DENDRO_6) - DENDRO_22*(DENDRO_19*DENDRO_6 + DENDRO_20*DENDRO_8 - DENDRO_21*DENDRO_27) - DENDRO_23*(DENDRO_24*DENDRO_6 + DENDRO_25*DENDRO_27 - DENDRO_26*DENDRO_8) - DENDRO_27*(-DENDRO_27*DENDRO_28 + DENDRO_29*DENDRO_6 - DENDRO_30*DENDRO_8) + DENDRO_6*(DENDRO_15*DENDRO_8 - DENDRO_27*DENDRO_7 + DENDRO_6*DENDRO_9) + DENDRO_8*(-DENDRO_16*DENDRO_27 + DENDRO_17*DENDRO_8 + DENDRO_18*DENDRO_6));
                //--
                CalGt2[pp] = DENDRO_5*(-DENDRO_11*(-DENDRO_10*DENDRO_31 - DENDRO_23*DENDRO_33 + DENDRO_32*DENDRO_6) - DENDRO_22*(-DENDRO_10*DENDRO_20 - DENDRO_19*DENDRO_23 + DENDRO_21*DENDRO_6) - DENDRO_23*(DENDRO_10*DENDRO_26 - DENDRO_23*DENDRO_24 - DENDRO_25*DENDRO_6) - DENDRO_27*(DENDRO_10*DENDRO_30 - DENDRO_23*DENDRO_29 + DENDRO_28*DENDRO_6) + DENDRO_6*(-DENDRO_10*DENDRO_15 - DENDRO_23*DENDRO_9 + DENDRO_6*DENDRO_7) + DENDRO_8*(-DENDRO_10*DENDRO_17 + DENDRO_16*DENDRO_6 - DENDRO_18*DENDRO_23));
                // Dendro: reduced ops: 202
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }


    //#include "Gt_rhs_s1_.cpp"
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.Gt_rhs_s1]
                vnames = ['Gt_rhs_s1']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 24 
                // Dendro: printing temp variables

                // Dendro: printing variables
                //--
                Gt_rhs_s1_0[pp] = beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp];
                //--
                Gt_rhs_s1_1[pp] = beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp];
                //--
                Gt_rhs_s1_2[pp] = beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp];
                // Dendro: reduced ops: 24
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }

    //#include "Gt_rhs_s2_.cpp"

    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.Gt_rhs_s2]
                vnames = ['Gt_rhs_s2']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 24 
                // Dendro: printing temp variables

                // Dendro: printing variables
                //--
                Gt_rhs_s2_0[pp] = CalGt0[pp]*grad_0_beta0[pp] + CalGt1[pp]*grad_1_beta0[pp] + CalGt2[pp]*grad_2_beta0[pp];
                //--
                Gt_rhs_s2_1[pp] = CalGt0[pp]*grad_0_beta1[pp] + CalGt1[pp]*grad_1_beta1[pp] + CalGt2[pp]*grad_2_beta1[pp];
                //--
                Gt_rhs_s2_2[pp] = CalGt0[pp]*grad_0_beta2[pp] + CalGt1[pp]*grad_1_beta2[pp] + CalGt2[pp]*grad_2_beta2[pp];
                // Dendro: reduced ops: 24
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }

    //#include "Gt_rhs_s3_.cpp"

    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.Gt_rhs_s3]
                vnames = ['Gt_rhs_s3']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 18 
                // Dendro: printing temp variables
                double DENDRO_0 = grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp];

                // Dendro: printing variables
                //--
                Gt_rhs_s3_0[pp] = CalGt0[pp]*DENDRO_0;
                //--
                Gt_rhs_s3_1[pp] = CalGt1[pp]*DENDRO_0;
                //--
                Gt_rhs_s3_2[pp] = CalGt2[pp]*DENDRO_0;
                // Dendro: reduced ops: 8
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }

    //#include "Gt_rhs_s4_.cpp"
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.Gt_rhs_s4]
                vnames = ['Gt_rhs_s4']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 909 
                // Dendro: printing temp variables
                double DENDRO_0 = pow(gt4[pp], 2);
                double DENDRO_1 = pow(gt1[pp], 2);
                double DENDRO_2 = pow(gt2[pp], 2);
                double DENDRO_3 = gt0[pp]*gt3[pp];
                double DENDRO_4 = gt1[pp]*gt2[pp];
                double DENDRO_5 = 1.0/(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt5[pp] - 2*DENDRO_4*gt4[pp]);
                double DENDRO_6 = -DENDRO_4 + gt0[pp]*gt4[pp];
                double DENDRO_7 = grad2_0_2_beta0[pp];
                double DENDRO_8 = gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp];
                double DENDRO_9 = (7.0/3.0)*DENDRO_8;
                double DENDRO_10 = grad2_1_2_beta1[pp];
                double DENDRO_11 = (1.0/3.0)*DENDRO_8;
                double DENDRO_12 = grad2_2_2_beta2[pp];
                double DENDRO_13 = grad2_0_1_beta0[pp];
                double DENDRO_14 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
                double DENDRO_15 = (7.0/3.0)*DENDRO_14;
                double DENDRO_16 = grad2_1_1_beta1[pp];
                double DENDRO_17 = (1.0/3.0)*DENDRO_14;
                double DENDRO_18 = grad2_1_2_beta2[pp];
                double DENDRO_19 = -DENDRO_1 + DENDRO_3;
                double DENDRO_20 = -DENDRO_2 + gt0[pp]*gt5[pp];
                double DENDRO_21 = grad2_0_0_beta0[pp];
                double DENDRO_22 = -DENDRO_0 + gt3[pp]*gt5[pp];
                double DENDRO_23 = grad2_0_1_beta1[pp];
                double DENDRO_24 = (1.0/3.0)*DENDRO_22;
                double DENDRO_25 = grad2_0_2_beta2[pp];
                double DENDRO_26 = (1.0/3.0)*DENDRO_6;
                double DENDRO_27 = (7.0/3.0)*DENDRO_6;
                double DENDRO_28 = (1.0/3.0)*DENDRO_20;
                double DENDRO_29 = (1.0/3.0)*DENDRO_19;

                // Dendro: printing variables
                //--
                Gt_rhs_s4_0[pp] = DENDRO_5*(-DENDRO_10*DENDRO_11 - DENDRO_11*DENDRO_12 + DENDRO_13*DENDRO_15 + DENDRO_16*DENDRO_17 + DENDRO_17*DENDRO_18 - DENDRO_19*grad2_2_2_beta0[pp] - DENDRO_20*grad2_1_1_beta0[pp] - 4.0/3.0*DENDRO_21*DENDRO_22 - DENDRO_23*DENDRO_24 - DENDRO_24*DENDRO_25 + 2*DENDRO_6*grad2_1_2_beta0[pp] - DENDRO_7*DENDRO_9);
                //--
                Gt_rhs_s4_1[pp] = DENDRO_5*(DENDRO_10*DENDRO_27 + DENDRO_12*DENDRO_26 - DENDRO_13*DENDRO_28 + DENDRO_15*DENDRO_23 - 4.0/3.0*DENDRO_16*DENDRO_20 + DENDRO_17*DENDRO_21 + DENDRO_17*DENDRO_25 - DENDRO_18*DENDRO_28 - DENDRO_19*grad2_2_2_beta1[pp] - DENDRO_22*grad2_0_0_beta1[pp] + DENDRO_26*DENDRO_7 - 2*DENDRO_8*grad2_0_2_beta1[pp]);
                //--
                Gt_rhs_s4_2[pp] = DENDRO_5*(-DENDRO_10*DENDRO_29 - DENDRO_11*DENDRO_21 - DENDRO_11*DENDRO_23 - 4.0/3.0*DENDRO_12*DENDRO_19 + DENDRO_13*DENDRO_26 + 2*DENDRO_14*grad2_0_1_beta2[pp] + DENDRO_16*DENDRO_26 + DENDRO_18*DENDRO_27 - DENDRO_20*grad2_1_1_beta2[pp] - DENDRO_22*grad2_0_0_beta2[pp] - DENDRO_25*DENDRO_9 - DENDRO_29*DENDRO_7);
                // Dendro: reduced ops: 140
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }

    //#include "Gt_rhs_s5_.cpp"
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.Gt_rhs_s5]
                vnames = ['Gt_rhs_s5']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 1914 
                // Dendro: printing temp variables
                double DENDRO_0 = grad_0_alpha[pp];
                double DENDRO_1 = gt2[pp]*gt4[pp];
                double DENDRO_2 = -DENDRO_1 + gt1[pp]*gt5[pp];
                double DENDRO_3 = pow(DENDRO_2, 2);
                double DENDRO_4 = gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp];
                double DENDRO_5 = pow(DENDRO_4, 2);
                double DENDRO_6 = gt3[pp]*gt5[pp];
                double DENDRO_7 = pow(gt4[pp], 2);
                double DENDRO_8 = DENDRO_6 - DENDRO_7;
                double DENDRO_9 = DENDRO_2*DENDRO_4;
                double DENDRO_10 = 2*At4[pp];
                double DENDRO_11 = At1[pp]*DENDRO_8;
                double DENDRO_12 = 2*DENDRO_2;
                double DENDRO_13 = At2[pp]*DENDRO_8;
                double DENDRO_14 = 2*DENDRO_4;
                double DENDRO_15 = grad_2_alpha[pp];
                double DENDRO_16 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
                double DENDRO_17 = At3[pp]*DENDRO_2;
                double DENDRO_18 = At0[pp]*DENDRO_8;
                double DENDRO_19 = At4[pp]*DENDRO_4;
                double DENDRO_20 = pow(gt1[pp], 2);
                double DENDRO_21 = -DENDRO_20 + gt0[pp]*gt3[pp];
                double DENDRO_22 = At5[pp]*DENDRO_4;
                double DENDRO_23 = At4[pp]*DENDRO_2;
                double DENDRO_24 = -At1[pp]*DENDRO_9 + At2[pp]*DENDRO_5 - DENDRO_11*DENDRO_16 + DENDRO_13*DENDRO_21 + DENDRO_16*DENDRO_17 - DENDRO_16*DENDRO_19 + DENDRO_18*DENDRO_4 + DENDRO_21*DENDRO_22 - DENDRO_21*DENDRO_23;
                double DENDRO_25 = grad_1_alpha[pp];
                double DENDRO_26 = pow(gt2[pp], 2);
                double DENDRO_27 = -DENDRO_26 + gt0[pp]*gt5[pp];
                double DENDRO_28 = -At1[pp]*DENDRO_3 + At2[pp]*DENDRO_9 - DENDRO_11*DENDRO_27 + DENDRO_13*DENDRO_16 + DENDRO_16*DENDRO_22 - DENDRO_16*DENDRO_23 + DENDRO_17*DENDRO_27 + DENDRO_18*DENDRO_2 - DENDRO_19*DENDRO_27;
                double DENDRO_29 = 2/pow(-2*DENDRO_1*gt1[pp] + DENDRO_20*gt5[pp] + DENDRO_26*gt3[pp] - DENDRO_6*gt0[pp] + DENDRO_7*gt0[pp], 2);
                double DENDRO_30 = pow(DENDRO_16, 2);
                double DENDRO_31 = At2[pp]*DENDRO_2;
                double DENDRO_32 = At1[pp]*DENDRO_27;
                double DENDRO_33 = DENDRO_16*DENDRO_27;
                double DENDRO_34 = DENDRO_16*DENDRO_4;
                double DENDRO_35 = DENDRO_16*DENDRO_21;
                double DENDRO_36 = At0[pp]*DENDRO_9 - At1[pp]*DENDRO_16*DENDRO_2 + At2[pp]*DENDRO_34 + At3[pp]*DENDRO_33 - At4[pp]*DENDRO_21*DENDRO_27 - At4[pp]*DENDRO_30 + At5[pp]*DENDRO_35 + DENDRO_21*DENDRO_31 - DENDRO_32*DENDRO_4;

                // Dendro: printing variables
                //--
                Gt_rhs_s5_0[pp] = DENDRO_29*(DENDRO_0*(At0[pp]*pow(DENDRO_8, 2) + At3[pp]*DENDRO_3 + At5[pp]*DENDRO_5 - DENDRO_10*DENDRO_9 - DENDRO_11*DENDRO_12 + DENDRO_13*DENDRO_14) + DENDRO_15*DENDRO_24 - DENDRO_25*DENDRO_28);
                //--
                Gt_rhs_s5_1[pp] = DENDRO_29*(-DENDRO_0*DENDRO_28 - DENDRO_15*DENDRO_36 + DENDRO_25*(At0[pp]*DENDRO_3 + At3[pp]*pow(DENDRO_27, 2) + At5[pp]*DENDRO_30 - DENDRO_10*DENDRO_33 - DENDRO_12*DENDRO_32 + 2*DENDRO_16*DENDRO_31));
                //--
                Gt_rhs_s5_2[pp] = DENDRO_29*(DENDRO_0*DENDRO_24 + DENDRO_15*(At0[pp]*DENDRO_5 - 2*At1[pp]*DENDRO_34 + At2[pp]*DENDRO_14*DENDRO_21 + At3[pp]*DENDRO_30 + At5[pp]*pow(DENDRO_21, 2) - DENDRO_10*DENDRO_35) - DENDRO_25*DENDRO_36);
                // Dendro: reduced ops: 162
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }

    //#include "Gt_rhs_s6_.cpp"

    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.Gt_rhs_s6]
                vnames = ['Gt_rhs_s6']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 4812 
                // Dendro: printing temp variables
                double DENDRO_0 = grad_1_gt3[pp];
                double DENDRO_1 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
                double DENDRO_2 = 0.5*DENDRO_1;
                double DENDRO_3 = gt1[pp]*gt4[pp];
                double DENDRO_4 = DENDRO_3 - gt2[pp]*gt3[pp];
                double DENDRO_5 = grad_2_gt3[pp];
                double DENDRO_6 = -0.5*DENDRO_5 + 1.0*grad_1_gt4[pp];
                double DENDRO_7 = gt3[pp]*gt5[pp];
                double DENDRO_8 = pow(gt4[pp], 2);
                double DENDRO_9 = DENDRO_7 - DENDRO_8;
                double DENDRO_10 = grad_0_gt3[pp];
                double DENDRO_11 = 0.5*DENDRO_10 - 1.0*grad_1_gt1[pp];
                double DENDRO_12 = pow(DENDRO_1, 2);
                double DENDRO_13 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
                double DENDRO_14 = pow(DENDRO_13, 2);
                double DENDRO_15 = pow(gt2[pp], 2);
                double DENDRO_16 = -DENDRO_15 + gt0[pp]*gt5[pp];
                double DENDRO_17 = At2[pp]*DENDRO_1;
                double DENDRO_18 = 2*DENDRO_13;
                double DENDRO_19 = At1[pp]*DENDRO_16;
                double DENDRO_20 = 2*DENDRO_1;
                double DENDRO_21 = At4[pp]*DENDRO_16;
                double DENDRO_22 = At0[pp]*DENDRO_12 + At3[pp]*pow(DENDRO_16, 2) + At5[pp]*DENDRO_14 + DENDRO_17*DENDRO_18 - DENDRO_18*DENDRO_21 - DENDRO_19*DENDRO_20;
                double DENDRO_23 = grad_0_gt0[pp];
                double DENDRO_24 = 0.5*DENDRO_23;
                double DENDRO_25 = grad_2_gt0[pp];
                double DENDRO_26 = -0.5*DENDRO_25 + 1.0*grad_0_gt2[pp];
                double DENDRO_27 = grad_1_gt0[pp];
                double DENDRO_28 = -0.5*DENDRO_27 + 1.0*grad_0_gt1[pp];
                double DENDRO_29 = pow(DENDRO_4, 2);
                double DENDRO_30 = At4[pp]*DENDRO_4;
                double DENDRO_31 = At1[pp]*DENDRO_1;
                double DENDRO_32 = DENDRO_4*DENDRO_9;
                double DENDRO_33 = At0[pp]*pow(DENDRO_9, 2) + 2*At2[pp]*DENDRO_32 + At3[pp]*DENDRO_12 + At5[pp]*DENDRO_29 - DENDRO_20*DENDRO_30 - 2*DENDRO_31*DENDRO_9;
                double DENDRO_34 = 0.5*grad_2_gt5[pp];
                double DENDRO_35 = grad_1_gt5[pp];
                double DENDRO_36 = 0.5*DENDRO_35 - 1.0*grad_2_gt4[pp];
                double DENDRO_37 = grad_0_gt5[pp];
                double DENDRO_38 = 0.5*DENDRO_37 - 1.0*grad_2_gt2[pp];
                double DENDRO_39 = pow(gt1[pp], 2);
                double DENDRO_40 = -DENDRO_39 + gt0[pp]*gt3[pp];
                double DENDRO_41 = At1[pp]*DENDRO_13;
                double DENDRO_42 = 2*DENDRO_4;
                double DENDRO_43 = At2[pp]*DENDRO_40;
                double DENDRO_44 = At4[pp]*DENDRO_40;
                double DENDRO_45 = At0[pp]*DENDRO_29 + At3[pp]*DENDRO_14 + At5[pp]*pow(DENDRO_40, 2) - DENDRO_18*DENDRO_44 - DENDRO_41*DENDRO_42 + DENDRO_42*DENDRO_43;
                double DENDRO_46 = grad_0_gt4[pp];
                double DENDRO_47 = grad_1_gt2[pp];
                double DENDRO_48 = grad_2_gt1[pp];
                double DENDRO_49 = DENDRO_46 + DENDRO_47 - DENDRO_48;
                double DENDRO_50 = At5[pp]*DENDRO_13;
                double DENDRO_51 = At0[pp]*DENDRO_1;
                double DENDRO_52 = At2[pp]*DENDRO_13;
                double DENDRO_53 = At3[pp]*DENDRO_16;
                double DENDRO_54 = DENDRO_1*DENDRO_13;
                double DENDRO_55 = -At1[pp]*DENDRO_12 - At4[pp]*DENDRO_54 + DENDRO_1*DENDRO_53 + DENDRO_17*DENDRO_4 - DENDRO_19*DENDRO_9 - DENDRO_21*DENDRO_4 + DENDRO_4*DENDRO_50 + DENDRO_51*DENDRO_9 + DENDRO_52*DENDRO_9;
                double DENDRO_56 = -DENDRO_46 + DENDRO_47 + DENDRO_48;
                double DENDRO_57 = -At1[pp]*DENDRO_54 - At4[pp]*DENDRO_14 + DENDRO_13*DENDRO_53 + DENDRO_17*DENDRO_40 - DENDRO_19*DENDRO_4 - DENDRO_21*DENDRO_40 + DENDRO_4*DENDRO_51 + DENDRO_4*DENDRO_52 + DENDRO_40*DENDRO_50;
                double DENDRO_58 = DENDRO_46 - DENDRO_47 + DENDRO_48;
                double DENDRO_59 = At0[pp]*DENDRO_32 + At2[pp]*DENDRO_29 + At3[pp]*DENDRO_54 + At5[pp]*DENDRO_4*DENDRO_40 - DENDRO_1*DENDRO_44 - DENDRO_13*DENDRO_30 - DENDRO_31*DENDRO_4 - DENDRO_41*DENDRO_9 + DENDRO_43*DENDRO_9;
                double DENDRO_60 = 1.0*DENDRO_59;
                double DENDRO_61 = 2*alpha[pp]/pow(DENDRO_15*gt3[pp] - 2*DENDRO_3*gt2[pp] + DENDRO_39*gt5[pp] - DENDRO_7*gt0[pp] + DENDRO_8*gt0[pp], 3);
                double DENDRO_62 = 0.5*DENDRO_0;

                // Dendro: printing variables
                //--
                Gt_rhs_s6_0[pp] = DENDRO_61*(DENDRO_22*(DENDRO_0*DENDRO_2 + DENDRO_11*DENDRO_9 - DENDRO_4*DENDRO_6) - DENDRO_33*(-DENDRO_1*DENDRO_28 + DENDRO_24*DENDRO_9 + DENDRO_26*DENDRO_4) - DENDRO_45*(DENDRO_1*DENDRO_36 + DENDRO_34*DENDRO_4 - DENDRO_38*DENDRO_9) + DENDRO_55*(-DENDRO_1*DENDRO_10 + DENDRO_27*DENDRO_9 + DENDRO_4*DENDRO_49) + DENDRO_57*(-DENDRO_1*DENDRO_5 + DENDRO_35*DENDRO_4 + DENDRO_56*DENDRO_9) - DENDRO_60*(-DENDRO_1*DENDRO_58 + DENDRO_25*DENDRO_9 + DENDRO_37*DENDRO_4));
                //--
                Gt_rhs_s6_1[pp] = DENDRO_61*(-DENDRO_22*(DENDRO_1*DENDRO_11 - DENDRO_13*DENDRO_6 + DENDRO_16*DENDRO_62) + DENDRO_33*(DENDRO_13*DENDRO_26 - DENDRO_16*DENDRO_28 + DENDRO_2*DENDRO_23) + DENDRO_45*(-DENDRO_1*DENDRO_38 + DENDRO_13*DENDRO_34 + DENDRO_16*DENDRO_36) - 1.0*DENDRO_55*(DENDRO_1*DENDRO_27 - DENDRO_10*DENDRO_16 + DENDRO_13*DENDRO_49) - 1.0*DENDRO_57*(DENDRO_1*DENDRO_56 + DENDRO_13*DENDRO_35 - DENDRO_16*DENDRO_5) + DENDRO_59*(DENDRO_1*DENDRO_25 + DENDRO_13*DENDRO_37 - DENDRO_16*DENDRO_58));
                //--
                Gt_rhs_s6_2[pp] = DENDRO_61*(DENDRO_22*(DENDRO_11*DENDRO_4 + DENDRO_13*DENDRO_62 - DENDRO_40*DENDRO_6) - DENDRO_33*(-DENDRO_13*DENDRO_28 + DENDRO_24*DENDRO_4 + DENDRO_26*DENDRO_40) - DENDRO_45*(DENDRO_13*DENDRO_36 + DENDRO_34*DENDRO_40 - DENDRO_38*DENDRO_4) + DENDRO_55*(-DENDRO_10*DENDRO_13 + DENDRO_27*DENDRO_4 + DENDRO_40*DENDRO_49) + DENDRO_57*(-DENDRO_13*DENDRO_5 + DENDRO_35*DENDRO_40 + DENDRO_4*DENDRO_56) - DENDRO_60*(-DENDRO_13*DENDRO_58 + DENDRO_25*DENDRO_4 + DENDRO_37*DENDRO_40));
                // Dendro: reduced ops: 316
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }

    //#include "Gt_rhs_s7_.cpp"

    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.Gt_rhs_s7]
                vnames = ['Gt_rhs_s7']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 2121 
                // Dendro: printing temp variables
                double DENDRO_0 = gt1[pp]*gt4[pp];
                double DENDRO_1 = DENDRO_0 - gt2[pp]*gt3[pp];
                double DENDRO_2 = (4.0/3.0)*grad_2_K[pp];
                double DENDRO_3 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
                double DENDRO_4 = (4.0/3.0)*grad_1_K[pp];
                double DENDRO_5 = gt3[pp]*gt5[pp];
                double DENDRO_6 = pow(gt4[pp], 2);
                double DENDRO_7 = DENDRO_5 - DENDRO_6;
                double DENDRO_8 = (4.0/3.0)*grad_0_K[pp];
                double DENDRO_9 = pow(DENDRO_3, 2);
                double DENDRO_10 = pow(DENDRO_1, 2);
                double DENDRO_11 = DENDRO_1*DENDRO_3;
                double DENDRO_12 = 2*At4[pp];
                double DENDRO_13 = At1[pp]*DENDRO_7;
                double DENDRO_14 = 2*DENDRO_3;
                double DENDRO_15 = At2[pp]*DENDRO_7;
                double DENDRO_16 = 2*DENDRO_1;
                double DENDRO_17 = grad_0_chi[pp];
                double DENDRO_18 = pow(gt1[pp], 2);
                double DENDRO_19 = pow(gt2[pp], 2);
                double DENDRO_20 = 1.0/(-2*DENDRO_0*gt2[pp] + DENDRO_18*gt5[pp] + DENDRO_19*gt3[pp] - DENDRO_5*gt0[pp] + DENDRO_6*gt0[pp]);
                double DENDRO_21 = 3*DENDRO_20/chi[pp];
                double DENDRO_22 = DENDRO_17*DENDRO_21;
                double DENDRO_23 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
                double DENDRO_24 = At3[pp]*DENDRO_3;
                double DENDRO_25 = At0[pp]*DENDRO_7;
                double DENDRO_26 = At4[pp]*DENDRO_1;
                double DENDRO_27 = -DENDRO_18 + gt0[pp]*gt3[pp];
                double DENDRO_28 = At5[pp]*DENDRO_1;
                double DENDRO_29 = At4[pp]*DENDRO_3;
                double DENDRO_30 = -At1[pp]*DENDRO_11 + At2[pp]*DENDRO_10 + DENDRO_1*DENDRO_25 - DENDRO_13*DENDRO_23 + DENDRO_15*DENDRO_27 + DENDRO_23*DENDRO_24 - DENDRO_23*DENDRO_26 + DENDRO_27*DENDRO_28 - DENDRO_27*DENDRO_29;
                double DENDRO_31 = DENDRO_21*grad_2_chi[pp];
                double DENDRO_32 = grad_1_chi[pp];
                double DENDRO_33 = -DENDRO_19 + gt0[pp]*gt5[pp];
                double DENDRO_34 = DENDRO_21*(-At1[pp]*DENDRO_9 + At2[pp]*DENDRO_11 - DENDRO_13*DENDRO_33 + DENDRO_15*DENDRO_23 + DENDRO_23*DENDRO_28 - DENDRO_23*DENDRO_29 + DENDRO_24*DENDRO_33 + DENDRO_25*DENDRO_3 - DENDRO_26*DENDRO_33);
                double DENDRO_35 = DENDRO_20*alpha[pp];
                double DENDRO_36 = pow(DENDRO_23, 2);
                double DENDRO_37 = At2[pp]*DENDRO_3;
                double DENDRO_38 = At1[pp]*DENDRO_33;
                double DENDRO_39 = DENDRO_23*DENDRO_33;
                double DENDRO_40 = DENDRO_21*DENDRO_32;
                double DENDRO_41 = DENDRO_1*DENDRO_23;
                double DENDRO_42 = DENDRO_23*DENDRO_27;
                double DENDRO_43 = At0[pp]*DENDRO_11 - At1[pp]*DENDRO_23*DENDRO_3 + At2[pp]*DENDRO_41 + At3[pp]*DENDRO_39 - At4[pp]*DENDRO_27*DENDRO_33 - At4[pp]*DENDRO_36 + At5[pp]*DENDRO_42 - DENDRO_1*DENDRO_38 + DENDRO_27*DENDRO_37;

                // Dendro: printing variables
                //--
                Gt_rhs_s7_0[pp] = DENDRO_35*(-DENDRO_1*DENDRO_2 + DENDRO_22*(At0[pp]*pow(DENDRO_7, 2) + At3[pp]*DENDRO_9 + At5[pp]*DENDRO_10 - DENDRO_11*DENDRO_12 - DENDRO_13*DENDRO_14 + DENDRO_15*DENDRO_16) + DENDRO_3*DENDRO_4 + DENDRO_30*DENDRO_31 - DENDRO_32*DENDRO_34 - DENDRO_7*DENDRO_8);
                //--
                Gt_rhs_s7_1[pp] = DENDRO_35*(-DENDRO_17*DENDRO_34 + DENDRO_2*DENDRO_23 + DENDRO_3*DENDRO_8 - DENDRO_31*DENDRO_43 - DENDRO_33*DENDRO_4 + DENDRO_40*(At0[pp]*DENDRO_9 + At3[pp]*pow(DENDRO_33, 2) + At5[pp]*DENDRO_36 - DENDRO_12*DENDRO_39 - DENDRO_14*DENDRO_38 + 2*DENDRO_23*DENDRO_37));
                //--
                Gt_rhs_s7_2[pp] = DENDRO_35*(-DENDRO_1*DENDRO_8 - DENDRO_2*DENDRO_27 + DENDRO_22*DENDRO_30 + DENDRO_23*DENDRO_4 + DENDRO_31*(At0[pp]*DENDRO_10 - 2*At1[pp]*DENDRO_41 + At2[pp]*DENDRO_16*DENDRO_27 + At3[pp]*DENDRO_36 + At5[pp]*pow(DENDRO_27, 2) - DENDRO_12*DENDRO_42) - DENDRO_40*DENDRO_43);
                // Dendro: reduced ops: 195
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }
    //#include "Gt_rhs.cpp"

    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.Gt_rhs]
                vnames = ['Gt_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 24 
                // Dendro: printing temp variables

                // Dendro: printing variables
                //--
                Gt_rhs0[pp] = Gt_rhs_s1_0[pp] - Gt_rhs_s2_0[pp] + (2.0/3.0)*Gt_rhs_s3_0[pp] + Gt_rhs_s4_0[pp] - Gt_rhs_s5_0[pp] + Gt_rhs_s6_0[pp] - Gt_rhs_s7_0[pp];
                //--
                Gt_rhs1[pp] = Gt_rhs_s1_1[pp] - Gt_rhs_s2_1[pp] + (2.0/3.0)*Gt_rhs_s3_1[pp] + Gt_rhs_s4_1[pp] - Gt_rhs_s5_1[pp] + Gt_rhs_s6_1[pp] - Gt_rhs_s7_1[pp];
                //--
                Gt_rhs2[pp] = Gt_rhs_s1_2[pp] - Gt_rhs_s2_2[pp] + (2.0/3.0)*Gt_rhs_s3_2[pp] + Gt_rhs_s4_2[pp] - Gt_rhs_s5_2[pp] + Gt_rhs_s6_2[pp] - Gt_rhs_s7_2[pp];
                // Dendro: reduced ops: 24
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }


    bssn::timer::t_rhs_Gt.stop();

    bssn::timer::t_rhs_B.start();
#ifdef USE_ETA_FUNC
//#include "B_rhs_eta_func.cpp"
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                outs = [bssn_stages.B_rhs]
                vnames = ['B_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 516 
                // Dendro: printing temp variables
                double DENDRO_0 = pow(gt4[pp], 2);
                double DENDRO_1 = pow(gt1[pp], 2);
                double DENDRO_2 = pow(gt2[pp], 2);
                double DENDRO_3 = gt3[pp]*gt5[pp];
                double DENDRO_4 = gt1[pp]*gt4[pp];
                double DENDRO_5 = grad_2_chi[pp];
                double DENDRO_6 = grad_1_chi[pp];
                double DENDRO_7 = grad_0_chi[pp];
                double DENDRO_8 = 2*DENDRO_5;
                double DENDRO_9 = BSSN_ETA_R0*sqrt((-pow(DENDRO_5, 2)*(-DENDRO_1 + gt0[pp]*gt3[pp]) - pow(DENDRO_6, 2)*(-DENDRO_2 + gt0[pp]*gt5[pp]) + 2*DENDRO_6*DENDRO_7*(gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp]) + DENDRO_6*DENDRO_8*(gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp]) - pow(DENDRO_7, 2)*(-DENDRO_0 + DENDRO_3) - DENDRO_7*DENDRO_8*(DENDRO_4 - gt2[pp]*gt3[pp]))/(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt0[pp] - 2*DENDRO_4*gt2[pp]))*pow(-pow(chi[pp], BSSN_ETA_POWER[0]) + 1, -BSSN_ETA_POWER[1]);

                // Dendro: printing variables
                //--
                B_rhs0[pp] = -B0[pp]*DENDRO_9 + Gt_rhs0[pp] + lambda[2]*(beta0[pp]*agrad_0_B0[pp] + beta1[pp]*agrad_1_B0[pp] + beta2[pp]*agrad_2_B0[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp]);
                //--
                B_rhs1[pp] = -B1[pp]*DENDRO_9 + Gt_rhs1[pp] + lambda[2]*(beta0[pp]*agrad_0_B1[pp] + beta1[pp]*agrad_1_B1[pp] + beta2[pp]*agrad_2_B1[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp]);
                //--
                B_rhs2[pp] = -B2[pp]*DENDRO_9 + Gt_rhs2[pp] + lambda[2]*(beta0[pp]*agrad_0_B2[pp] + beta1[pp]*agrad_1_B2[pp] + beta2[pp]*agrad_2_B2[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp]);
                // Dendro: reduced ops: 124
                // Dendro: }}} 
//[[[end]]]
            }
        }
    }
#else

//#include "B_rhs_eta_const.cpp"
    for (unsigned int k = 3; k < nz-3; k++) {
        z = pmin[2] + k*hz;

        for (unsigned int j = 3; j < ny-3; j++) {
            y = pmin[1] + j*hy;

            for (unsigned int i = 3; i < nx-3; i++) {
                x = pmin[0] + i*hx;
                pp = i + nx*(j + ny*k);
                r_coord = sqrt(x*x + y*y + z*z);
                eta=ETA_CONST;
                if (r_coord >= ETA_R0) {
                    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
                }
                /*[[[cog

                import dendro
                import bssn_stages

                bssn_stages.eta_func=bssn_stages.eta
                bssn_stages.B_rhs = [bssn_stages.Gt_rhs[i] - bssn_stages.eta_func * bssn_stages.B[i] +
                         bssn_stages.l3 * dendro.vec_j_ad_j(bssn_stages.b, bssn_stages.B[i]) -
                         bssn_stages.l4 * dendro.vec_j_ad_j(bssn_stages.b, bssn_stages.Gt[i])
                         for i in dendro.e_i]

                outs = [bssn_stages.B_rhs]
                vnames = ['B_rhs']
                dendro.generate_cpu(outs, vnames, '[pp]')

                ]]]*/
                // Dendro: {{{ 
                // Dendro: original ops: 90 
                // Dendro: printing temp variables

                // Dendro: printing variables
                //--
                B_rhs0[pp] = -B0[pp]*eta + Gt_rhs_s1_0[pp] - Gt_rhs_s2_0[pp] + (2.0/3.0)*Gt_rhs_s3_0[pp] + Gt_rhs_s4_0[pp] - Gt_rhs_s5_0[pp] + Gt_rhs_s6_0[pp] - Gt_rhs_s7_0[pp] + lambda[2]*(beta0[pp]*agrad_0_B0[pp] + beta1[pp]*agrad_1_B0[pp] + beta2[pp]*agrad_2_B0[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp]);
                //--
                B_rhs1[pp] = -B1[pp]*eta + Gt_rhs_s1_1[pp] - Gt_rhs_s2_1[pp] + (2.0/3.0)*Gt_rhs_s3_1[pp] + Gt_rhs_s4_1[pp] - Gt_rhs_s5_1[pp] + Gt_rhs_s6_1[pp] - Gt_rhs_s7_1[pp] + lambda[2]*(beta0[pp]*agrad_0_B1[pp] + beta1[pp]*agrad_1_B1[pp] + beta2[pp]*agrad_2_B1[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp]);
                //--
                B_rhs2[pp] = -B2[pp]*eta + Gt_rhs_s1_2[pp] - Gt_rhs_s2_2[pp] + (2.0/3.0)*Gt_rhs_s3_2[pp] + Gt_rhs_s4_2[pp] - Gt_rhs_s5_2[pp] + Gt_rhs_s6_2[pp] - Gt_rhs_s7_2[pp] + lambda[2]*(beta0[pp]*agrad_0_B2[pp] + beta1[pp]*agrad_1_B2[pp] + beta2[pp]*agrad_2_B2[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp]);
                // Dendro: reduced ops: 90
                // Dendro: }}} 
                //[[[end]]]
            }
        }
    }
#endif

    bssn::timer::t_rhs_B.stop();

    bssn::timer::t_rhs.stop();


    delete [] CalGt0;
    delete [] CalGt1;
    delete [] CalGt2;

    delete [] Gt_rhs_s1_0;
    delete [] Gt_rhs_s1_1;
    delete [] Gt_rhs_s1_2;

    delete [] Gt_rhs_s2_0;
    delete [] Gt_rhs_s2_1;
    delete [] Gt_rhs_s2_2;

    delete [] Gt_rhs_s3_0;
    delete [] Gt_rhs_s3_1;
    delete [] Gt_rhs_s3_2;

    delete [] Gt_rhs_s4_0;
    delete [] Gt_rhs_s4_1;
    delete [] Gt_rhs_s4_2;

    delete [] Gt_rhs_s5_0;
    delete [] Gt_rhs_s5_1;
    delete [] Gt_rhs_s5_2;

    delete [] Gt_rhs_s6_0;
    delete [] Gt_rhs_s6_1;
    delete [] Gt_rhs_s6_2;

    delete [] Gt_rhs_s7_0;
    delete [] Gt_rhs_s7_1;
    delete [] Gt_rhs_s7_2;



    if (bflag != 0) {

        bssn::timer::t_bdyc.start();

#ifdef BSSN_KERR_SCHILD_TEST

        freeze_bcs(a_rhs, sz, bflag);
        freeze_bcs(chi_rhs, sz, bflag);
        freeze_bcs(K_rhs, sz, bflag);
        freeze_bcs(b_rhs0, sz, bflag);
        freeze_bcs(b_rhs1, sz, bflag);
        freeze_bcs(b_rhs2, sz, bflag);
        freeze_bcs(Gt_rhs0, sz, bflag);
        freeze_bcs(Gt_rhs1, sz, bflag);
        freeze_bcs(Gt_rhs2, sz, bflag);
        freeze_bcs(B_rhs0, sz, bflag);
        freeze_bcs(B_rhs1, sz, bflag);
        freeze_bcs(B_rhs2, sz, bflag);
        freeze_bcs(At_rhs00, sz, bflag);
        freeze_bcs(At_rhs01, sz, bflag);
        freeze_bcs(At_rhs02, sz, bflag);
        freeze_bcs(At_rhs11, sz, bflag);
        freeze_bcs(At_rhs12, sz, bflag);
        freeze_bcs(At_rhs22, sz, bflag);
        freeze_bcs(gt_rhs00, sz, bflag);
        freeze_bcs(gt_rhs01, sz, bflag);
        freeze_bcs(gt_rhs02, sz, bflag);
        freeze_bcs(gt_rhs11, sz, bflag);
        freeze_bcs(gt_rhs12, sz, bflag);
        freeze_bcs(gt_rhs22, sz, bflag);

#else
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
#endif
        bssn::timer::t_bdyc.stop();
    }


    bssn::timer::t_deriv.start();

// ko derivs

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
                pp = i + nx*(j + ny*k);

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

    //bssn::writeBLockToBinary((const double**)unzipVarsRHS,offset,pmin,pmax,bxMin,bxMax,sz,15,1.0,"af_ko");

    bssn::timer::t_deriv.start();
    //dealloc derivative variables.
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

#if 0
    for (unsigned int m = 0; m < 24; m++) {
        std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
    }
#endif



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
