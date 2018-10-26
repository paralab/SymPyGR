#include "rhs.h"

using namespace std;
enum VAR {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};
enum VAR_CONSTRAINT {C_HAM=0, C_MOM0, C_MOM1,C_MOM2};

double ETA_CONST=0.1;
double ETA_R0=0.1;
double ETA_DAMPING_EXP=0.1;
double KO_DISS_SIGMA=1e-4;
unsigned int BSSN_LAMBDA[4]={1,2,3,4};
double BSSN_LAMBDA_F[2]={0.8,0.9};

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
                                  BSSN_LAMBDA[2], BSSN_LAMBDA[3]};
  const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};



  int idx[3];

  unsigned int n = sz[0]*sz[1]*sz[2];

#if 0
  double vars[24];
  pmin[0] = 0.0; pmin[1] = 0.0; pmin[2] = 0.0;
  pmax[0] = 2.0; pmax[1] = 2.0; pmax[2] = 2.0;

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

#include "bssnrhs_memalloc.h"
#include "bssnrhs_memalloc_adv.h"
#include "bssnrhs_derivs.h"
#include "bssnrhs_derivs_adv.h"

bssn::timer::t_deriv.stop();

  register double x;
  register double y;
  register double z;
  register unsigned int pp;

  double r_coord;
  double eta;

  bssn::timer::t_rhs.start();
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


        #include "bssneqs.cpp"


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
  bssn::timer::t_rhs.stop();


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
#include "bssnrhs_ko_derivs.h"
bssn::timer::t_deriv.stop();


bssn::timer::t_rhs.start();

  const  double sigma = KO_DISS_SIGMA;


  // for (unsigned int k = 3; k < nz-3; k++) {
  //   for (unsigned int j = 3; j < ny-3; j++) {
  //     for (unsigned int i = 3; i < nx-3; i++) {
  //       pp = i + nx*(j + ny*k);

  //       a_rhs[pp]  += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
  //       b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
  //       b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
  //       b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

  //       gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
  //       gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
  //       gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
  //       gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
  //       gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
  //       gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

  //       chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

  //       At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
  //       At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
  //       At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
  //       At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
  //       At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
  //       At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

  //       K_rhs[pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

  //       Gt_rhs0[pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
  //       Gt_rhs1[pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
  //       Gt_rhs2[pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

  //       B_rhs0[pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
  //       B_rhs1[pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
  //       B_rhs2[pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
  //     }
  //   }
  // }

bssn::timer::t_rhs.stop();



bssn::timer::t_deriv.start();
#include "bssnrhs_dealloc.h"
#include "bssnrhs_dealloc_adv.h"
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
     x = pmin[0] + ie*hx;
    for (unsigned int k = kb; k < ke; k++) {
       z = pmin[2] + k*hz;
      for (unsigned int j = jb; j < je; j++) {
         y = pmin[1] + j*hy;
         pp = IDX(ie,j,k);
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
     y = pmin[1] + je*hy;
    for (unsigned int k = kb; k < ke; k++) {
       z = pmin[2] + k*hz;
      for (unsigned int i = ib; i < ie; i++) {
         x = pmin[0] + i*hx;
         inv_r = 1.0 / sqrt(x*x + y*y + z*z);
         pp = IDX(i,je,k);

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
    z = pmin[2] + ke*hz;
    for (unsigned int j = jb; j < je; j++) {
      y = pmin[1] + j*hy;
      for (unsigned int i = ib; i < ie; i++) {
        x = pmin[0] + i*hx;
        inv_r = 1.0 / sqrt(x*x + y*y + z*z);
        pp = IDX(i,j,ke);

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
void fake_initial_data(double x, double y, double z, double *u)
{

  const double pi = acos(-1.0);
  const double f1 = 31.0/17.0;
  const double f2 = 37.0/11.0;

  u[VAR::U_ALPHA] = 1.0 - 0.25*sin(f1*x);
  //u[F_ALPHA][pp] = 1.0;
  u[VAR::U_BETA0] = 4.0/17.0*sin(x)*cos(z);
  u[VAR::U_BETA1] = pi/5.0*cos(y)*sin(z+x);
  u[VAR::U_BETA2] = 4.0/17.0*sin(f2*x)*sin(z);

  u[VAR::U_B0] = 31.0*x*cos(f1*z+y);
  u[VAR::U_B1] = 7.0*y*sin(f1*x+y) + 3.0*cos(z);
  u[VAR::U_B2] = 5.0*z*cos(f1*x+y) + 7.0*sin(z+y+x) + 1.0;

  u[VAR::U_GT0] = 5.0*cos(x)/(10.0*sin(x+z)+26.0-1.0*cos(x*z)*cos(x));
  u[VAR::U_GT1] = -5.0*sin(y)/(25.0+10.0*cos(y+z)+cos(y)*cos(y*z));
  u[VAR::U_GT2] = -5.0*sin(z)/(25.0+10.0*cos(y+x)+cos(y*x)*cos(z));

  u[VAR::U_CHI] = 1.0 + exp(-4.0*cos(x)*sin(y));
  //u[F_CHI][pp] = 2.0;

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

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
