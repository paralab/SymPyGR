#include <cmath>
#include <iostream>
#include <dollar.h>
#include <time.h>
#include "eqtest.h"

#define grad(i, s) d ## i ## _ ## s
#define grad2(i, j, s) d ## i ## j ## _ ## s
#define agrad(i, s) ad ## i ## _ ## s

using namespace std;

/*----------------------------------------------------------------------;
 *
 * optimized RHS
 *
 * @author Hari Sundar  hari@cs.utah.edu
 * @date   12 May 2017
 *
 *----------------------------------------------------------------------*/
void rhs_opt(double *a_rhs_3D, double *b_rhs0_3D, double *b_rhs1_3D,
         double *b_rhs2_3D, double *chi_rhs_3D, double *K_rhs_3D,
         double *gt_rhs00_3D, double *gt_rhs01_3D, double *gt_rhs02_3D,
         double *gt_rhs11_3D, double *gt_rhs12_3D, double *gt_rhs22_3D,
         double *At_rhs00_3D, double *At_rhs01_3D, double *At_rhs02_3D,
         double *At_rhs11_3D, double *At_rhs12_3D, double *At_rhs22_3D,
         double *Gt_rhs0_3D, double *Gt_rhs1_3D, double *Gt_rhs2_3D,
         double *B_rhs0_3D, double *B_rhs1_3D, double *B_rhs2_3D,
         double *alpha_3D, double *beta0_3D, double *beta1_3D,
         double *beta2_3D, double *chi_3D, double *K_3D,
         double *gt0_3D, double *gt1_3D, double *gt2_3D,
         double *gt3_3D, double *gt4_3D, double *gt5_3D,
         double *At0_3D, double *At1_3D, double *At2_3D,
         double *At3_3D, double *At4_3D, double *At5_3D,
         double *Gt0_3D, double *Gt1_3D, double *Gt2_3D,
         double *B0_3D, double *B1_3D, double *B2_3D,
         int *lambda, int *lambda_f,
         double eta_const, int eta_damping, int eta_damping_exp, double eta_R0,
         double trK0, double chi_floor,
         double *xi[], int sz[])
{ $

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];

  double *x1d = xi[0];
  double *y1d = xi[1];
  double *z1d = xi[2];
  double h = x1d[1] - x1d[0];

  int idx[3];

  /* Hari - Notes
	 *
	 * 1. precompute gradients
	 *   - grad
	 *   		- alpha
	 *   		- beta [0,1,2]
	 *   		- chi
	 *   		- gt [0,1,2,3,4,5]
	 *   		- Gt [0,1,2]
	 *   		- K
	 *   - grad2
	 *   		- alpha
	 *   		- chi
	 *   		- beta [0,1,2]
	 *   		- gt [0, 1, 2, 3, 4, 5]
	 *   - agrad
	 *   		- Gt [0,1,2]
	 *   		- alpha
	 *   		- beta [0,1,2]
	 *   		- gt [0,1,2,3,4,5]
	 *   		- chi
	 *   		- At [0,1,2,3,4,5]
	 *   		- K
	 *   		- B [0,1,2]
	 *
	 * 2. all other variables
	 *    -- can be vectorized.
	 *
	 **/
  int n = sz[0]*sz[1]*sz[2];
  size_t align = 32;

#ifdef __APPLE__
  auto t1 = std::chrono::steady_clock::now();
#else
  struct timespec t1, t2, t3;
  GetRdtscTime(&t1);
#endif

#include "opt_memalloc.h"
#include "opt_memalloc_adv.h"
#include "opt_derivs.h"
#include "opt_derivs_adv.h"

  cout << "done taking derivatives" << endl;
  avx_apply_stencil_x(gt0_3D, dx_gt0, sz, h);
  apply_stencil_x(gt0_3D, dx_gt1, sz, h);

  cout << "check x-derivative" << endl;
  for (int k = 3; k < nz-3; k++) {
  for (int j = 3; j < ny-3; j++) {
  for (int i = 3; i < nx-3; i++) {
    int pp = i + nx*(j + ny*k);
    double df = dx_gt0[pp] - dx_gt1[pp];
    if (fabs(df) > 1.0e-10) {
      cout << "std=" << dx_gt1[pp] << ",  opt=" << dx_gt0[pp] << ",  diff=" << fabs(df) << endl;
    }
  }
  }
  }

#ifdef __APPLE__
  auto t2 = std::chrono::steady_clock::now();
#else
  GetRdtscTime(&t2);
#endif

  cout << "begin loop" << endl;
  for (int k = 3; k < nz-3; k++) {
    double z = xi[2][k];
    idx[2] = k;

    for (int j = 3; j < ny-3; j++) {
      double y = xi[1][j];
      idx[1] = j;

      for (int i = 3; i < nx-3; i++) {
        double x = xi[0][i];
        idx[0] = i;

        int pp = i + nx*(j + ny*k);

        double r_coord = sqrt(x*x + y*y + z*z);
        double eta = eta_const;
        if (r_coord >= eta_R0) {
          eta *= pow( (eta_R0/r_coord), eta_damping_exp);
        }

        double alpha = alpha_3D[pp];
        double beta0 = beta0_3D[pp];
        double beta1 = beta1_3D[pp];
        double beta2 = beta2_3D[pp];

        double B0 = B0_3D[pp];
        double B1 = B1_3D[pp];
        double B2 = B2_3D[pp];

        double chi = fmax( chi_3D[pp] , chi_floor );
        double gt0 = gt0_3D[pp];
        double gt1 = gt1_3D[pp];
        double gt2 = gt2_3D[pp];
        double gt3 = gt3_3D[pp];
        double gt4 = gt4_3D[pp];
        double gt5 = gt5_3D[pp];

        double K = K_3D[pp];
        double At0 = At0_3D[pp];
        double At1 = At1_3D[pp];
        double At2 = At2_3D[pp];
        double At3 = At3_3D[pp];
        double At4 = At4_3D[pp];
        double At5 = At5_3D[pp];

        double Gt0 = Gt0_3D[pp];
        double Gt1 = Gt1_3D[pp];
        double Gt2 = Gt2_3D[pp];

#include "opt_set_derivs.h"

        double a_rhs, b_rhs0, b_rhs1, b_rhs2, B_rhs0, B_rhs1, B_rhs2;
        double chi_rhs, K_rhs;
        double gt_rhs00, gt_rhs01, gt_rhs02, gt_rhs11, gt_rhs12, gt_rhs22;
        double At_rhs00, At_rhs01, At_rhs02, At_rhs11, At_rhs12, At_rhs22;
        double Gt_rhs0, Gt_rhs1, Gt_rhs2;

        double K2_rhs, K3_rhs;
        double At_UU00, At_UU01, At_UU02, At_UU11, At_UU12, At_UU22;
        double gtu00, gtu01, gtu02, gtu11, gtu12, gtu22;
        double gt00, gt01, gt02, gt11, gt12, gt22;
        double igt[9];
        double R00, R01, R02, R11, R12, R22;
        double Rt00, Rt01, Rt02, Rt11, Rt12, Rt22;
        double Rphi00, Rphi01, Rphi02, Rphi11, Rphi12, Rphi22;
        double CalGt0, CalGt1, CalGt2;
        double c1x00, c1x01, c1x02, c1x11, c1x12, c1x22;
        double c1y00, c1y01, c1y02, c1y11, c1y12, c1y22;
        double c1z00, c1z01, c1z02, c1z11, c1z12, c1z22;
        double c2x00, c2x01, c2x02, c2x11, c2x12, c2x22;
        double c2y00, c2y01, c2y02, c2y11, c2y12, c2y22;
        double c2z00, c2z01, c2z02, c2z11, c2z12, c2z22;
        double At2_rhs00, At2_rhs01, At2_rhs02, At2_rhs11, At2_rhs12, At2_rhs22;

//#include "bssn.cpp"

        a_rhs_3D[pp] = a_rhs;

        b_rhs0_3D[pp] = b_rhs0;
				b_rhs1_3D[pp] = b_rhs1;
				b_rhs2_3D[pp] = b_rhs2;

        B_rhs0_3D[pp] = B_rhs0;
        B_rhs1_3D[pp] = B_rhs1;
        B_rhs2_3D[pp] = B_rhs2;

        chi_rhs_3D[pp] = chi_rhs;
        gt_rhs00_3D[pp] = gt_rhs00;
        gt_rhs01_3D[pp] = gt_rhs01;
        gt_rhs02_3D[pp] = gt_rhs02;
        gt_rhs11_3D[pp] = gt_rhs11;
        gt_rhs12_3D[pp] = gt_rhs12;
        gt_rhs22_3D[pp] = gt_rhs22;

        K_rhs_3D[pp] = K_rhs;
        At_rhs00_3D[pp] = At_rhs00;
        At_rhs01_3D[pp] = At_rhs01;
        At_rhs02_3D[pp] = At_rhs02;
        At_rhs11_3D[pp] = At_rhs11;
        At_rhs12_3D[pp] = At_rhs12;
        At_rhs22_3D[pp] = At_rhs22;

        Gt_rhs0_3D[pp] = Gt_rhs0;
        Gt_rhs1_3D[pp] = Gt_rhs1;
        Gt_rhs2_3D[pp] = Gt_rhs2;

        cout << "A. dxgtxx=" << d0_gt0 << ",   dygtxy=" << d0_gt1
               << ",  dzgtxz=" << d0_gt2 << endl;
        cout << "A. dygtxx=" << d1_gt0 << ",   dygtxy=" << d1_gt1
               << ",  dygtxz=" << d1_gt2 << endl;
        cout << "A. dzgtxx=" << d2_gt0 << ",   dzgtxy=" << d2_gt1
               << ",  dzgtxz=" << d2_gt2 << endl;

       /* debugging */
        int qi = 46 - 1;
        int qj = 10 - 1;
        int qk = 60 - 1;
        int qpp = qi + nx*(qj + ny*qk);
        if (qpp == pp) {
          double Rt[3][3], Rp[3][3], CalGt[3];
          cout << ".... OPTIMIZED debug stuff..." << endl;
          cout << "dxgtxx=" << d0_gt0 << ",   dygtxx=" << d0_gt1
               << ",  dzgtxx=" << d0_gt2 << endl;
          CalGt[0] = CalGt0;
          CalGt[1] = CalGt1;
          CalGt[2] = CalGt2;

          Rt[0][0] = Rt00;
          Rt[0][1] = Rt01;
          Rt[0][2] = Rt02;
          Rt[1][0] = Rt01;
          Rt[1][1] = Rt11;
          Rt[1][2] = Rt12;
          Rt[2][0] = Rt02;
          Rt[2][1] = Rt12;
          Rt[2][2] = Rt22;

          Rp[0][0] = Rphi00;
          Rp[0][1] = Rphi01;
          Rp[0][2] = Rphi02;
          Rp[1][0] = Rphi01;
          Rp[1][1] = Rphi11;
          Rp[1][2] = Rphi12;
          Rp[2][0] = Rphi02;
          Rp[2][1] = Rphi12;
          Rp[2][2] = Rphi22;

          std::cout << ".... end OPTIMIZED debug stuff..." << std::endl;

        }

      }
    }
  }

#ifdef __APPLE__
  auto t3 = std::chrono::steady_clock::now();
#else
  GetRdtscTime(&t3);
#endif

#include "opt_dealloc.h"
#include "opt_dealloc_adv.h"

#ifdef __APPLE__
auto dt1 = t2 - t1;
auto dt2 = t3 - t2;
  std::cout << "Dendro Opt RHS: derivs        " << std::chrono::duration <double, std::milli> (dt1).count() << " ms" << std::endl;
  std::cout << "Dendro Opt RHS: rhs           " << std::chrono::duration <double, std::milli> (dt2).count() << " ms" << std::endl;
#else
  std::cout << "Dendro Opt RHS: derivs     " << TimeSpecDiff(&t2,&t1)->tv_sec << "." << TimeSpecDiff(&t2,&t1)->tv_nsec << std::endl;
  std::cout << "Dendro Opt RHS: rhs        " << TimeSpecDiff(&t3,&t2)->tv_sec << "." << TimeSpecDiff(&t3,&t2)->tv_nsec << std::endl;
#endif

}
