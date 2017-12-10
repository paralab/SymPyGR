#include <cmath>
#include <iostream>
#include <dollar.h>

double grad(int dir, double *u, double h, int idx[], int shp[]);
double grad2(int dir1, int dir2, double *u, double h, int idx[], int shp[]);
double agrad(int dir1, double *beta, double *u, double h, int idx[], int shp[]);

#define BETA(x) beta##x

#define grad(i, s) grad(i, s, h, idx, shp)
#define grad2(i, j, s) grad2(i, j, s, h, idx, shp)
#define agrad(i, s) agrad(i, BETA(i), s, h, idx, shp)

using namespace std;

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void rhs(double *a_rhs, double *b_rhs0, double *b_rhs1,
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
         double *xi[], int shp[])
{ $

  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];

  double *x1d = xi[0];
  double *y1d = xi[1];
  double *z1d = xi[2];
  double h = x1d[1] - x1d[0];

  int idx[3];

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

        chi[pp] = fmax( chi[pp] , chi_floor );

#if 0
       /* debugging vars */
        double dxgtxx = grad(0, gt0);
        double dygtxx = grad(1, gt0);
        double dzgtxx = grad(2, gt0);

        double adxgtxx = agrad(0, gt0);
        double adygtxx = agrad(1, gt0);
        double adzgtxx = agrad(2, gt0);

        double dxbetax = grad(0, beta0);
        double dxbetay = grad(0, beta1);
        double dxbetaz = grad(0, beta2);
        double divbeta = grad(0, beta0) + grad(1,beta1) + grad(2,beta2);

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
#endif


#include "bssn.cpp"



#if 1
       /* debugging */
        int qi = 46 - 1;
        int qj = 10 - 1;
        int qk = 60 - 1;
        int qpp = qi + nx*(qj + ny*qk);
        if (qpp == pp) {
          cout << ".... DENDRO debug stuff..." << endl;

          std::cout << ".... end DENDRO debug stuff..." << std::endl;

        }
#endif

      }
    }
  }


}
