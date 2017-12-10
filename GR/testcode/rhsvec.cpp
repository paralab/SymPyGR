#include <cmath>
#include <iostream>
#include <dollar.h>
#include <time.h>
#include "eqtest.h"

#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )
#define deriv_x deriv42_x
#define deriv_y deriv42_y
#define deriv_z deriv42_z

#define deriv_xx deriv42_xx
#define deriv_yy deriv42_yy
#define deriv_zz deriv42_zz

#define adv_deriv_x deriv42adv_x
#define adv_deriv_y deriv42adv_y
#define adv_deriv_z deriv42adv_z


using namespace std;

/*----------------------------------------------------------------------;
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void rhs_vec(double *a_rhs, double *b_rhs0, double *b_rhs1,
         double *b_rhs2, double *chi_rhs, double *K_rhs,
         double *gt_rhs00, double *gt_rhs01, double *gt_rhs02,
         double *gt_rhs11, double *gt_rhs12, double *gt_rhs22,
         double *At_rhs00, double *At_rhs01, double *At_rhs02,
         double *At_rhs11, double *At_rhs12, double *At_rhs22,
         double *Gt_rhs0, double *Gt_rhs1, double *Gt_rhs2,
         double *B_rhs0, double *B_rhs1, double *B_rhs2,
         double *alpha, double *beta0, double *beta1,
         double *beta2, double *chi, double *K,
         double *gt0, double *gt1, double *gt2,
         double *gt3, double *gt4, double *gt5,
         double *At0, double *At1, double *At2,
         double *At3, double *At4, double *At5,
         double *Gt0, double *Gt1, double *Gt2,
         double *B0, double *B1, double *B2,
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

  auto dadd = std::plus<double>();
  auto dmul = std::multiplies<double>();

  const double negone = -1.0;

#ifdef __APPLE__
  auto t1 = std::chrono::steady_clock::now();
#else
  struct timespec t1, t2, t3;
  GetRdtscTime(&t1);
#endif

#include "vec_memalloc.h"
#include "vec_memalloc_adv.h"
#include "vec_derivs.h"
#include "vec_derivs_adv.h"

#ifdef __APPLE__
  auto t2 = std::chrono::steady_clock::now();
#else
  GetRdtscTime(&t2);
#endif

  cout << "begin loop" << endl;
  for (int k = 3; k < nz-3; k++) {
    double z = xi[2][k];

    for (int j = 3; j < ny-3; j++) {
      double y = xi[1][j];

      for (int i = 3; i < nx-3; i++) {
        double x = xi[0][i];

        int pp = i + nx*(j + ny*k);

        double r_coord = sqrt(x*x + y*y + z*z);
        double eta = eta_const;
        if (r_coord >= eta_R0) {
          eta *= pow( (eta_R0/r_coord), eta_damping_exp);
        }

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

#include "bssn_vec.cpp"

       /* debugging */
        int qi = 46 - 1;
        int qj = 10 - 1;
        int qk = 60 - 1;
        int qidx = qi + nx*(qj + ny*qk);
        if (qidx == pp) {
          double Rt[3][3], Rp[3][3], CalGt[3];

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

#include "vec_dealloc.h"
#include "vec_dealloc_adv.h"

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
