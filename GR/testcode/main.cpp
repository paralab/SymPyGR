#include <iostream>
#include <cmath>
#include <assert.h>
#include <dollar.h>

#include <stdint.h> /* for uint64_t */
#include <time.h>  /* for struct timespec */

#include "vars.h"
#include "eqtest.h"
// #include "sdf.h"

#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

using namespace std;


#ifdef __APPLE__
  /* OS X doesn't have RDTSC, so use the chrono library in C++-11 */
#include <chrono>
#else

void InitRdtsc();

/* assembly code to read the TSC */
static inline uint64_t RDTSC()
{
  unsigned int hi, lo;
  __asm__ volatile("rdtsc" : "=a" (lo), "=d" (hi));
  return ((uint64_t)hi << 32) | lo;
}

const int NANO_SECONDS_IN_SEC = 1000000000;
void GetRdtscTime(struct timespec *ts);
struct timespec *TimeSpecDiff(struct timespec *ts1, struct timespec *ts2);
#endif

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

#ifndef __APPLE__
  InitRdtsc();
#endif

  int nx = 128;
  int ny = 128;
  int nz = 128;
  double xmin = 0.0;
  double xmax = 2.0;
  double ymin = 0.0;
  double ymax = 2.0;
  double zmin = 0.0;
  double zmax = 2.0;

 /*---------------------------------------------------------------
  *  allocate memory & define coordinates
  *---------------------------------------------------------------*/
  double *u[NU], *dtu0[NU], *dtu1[NU], *dtu2[NU], *dtu3[NU], *dtu4[NU];
  double *drhs[NU], *xi[3];
  int npts = nx * ny * nz;
  for (int i = 0; i < NU; i++) {
    u[i] = new double [npts];
    dtu0[i] = new double [npts];
    dtu1[i] = new double [npts];
    dtu2[i] = new double [npts];
    dtu3[i] = new double [npts];
    dtu4[i] = new double [npts];
    drhs[i] = new double [npts];
  }
  double *coords = new double [nx+ny+nz];
  xi[0] = coords;
  xi[1] = coords+nx;
  xi[2] = coords+nx+ny;

  double dx, dy, dz;
  define_coordinates(xi[0], xmin, xmax, dx, nx);
  define_coordinates(xi[1], ymin, ymax, dy, ny);
  define_coordinates(xi[2], zmin, zmax, dz, nz);

  int shp[3] = {nx, ny, nz};

 /*---------------------------------------------------------------
  *  calculate initial data
  *---------------------------------------------------------------*/
  initial_data(u, xi, shp);
  std::cout << std::endl << std::endl;
  std::cout << "L2 Norms of the Initial Data" << std::endl;
  std::cout << "----------------------------" << std::endl;
  for (int m = 0; m < NU; m++) {
    std::cout << "|| u[" << m << "] || = " << L2norm(u[m],shp) << std::endl;
  }

  for (int i = 0; i < NU; i++) {
    set_zero(dtu0[i], shp);
    set_zero(dtu1[i], shp);
    set_zero(dtu2[i], shp);
  }

  int lambda[4] = {1, 1, 1, 1};
  int lambda_f[2] = {1, 1};
  double eta = exp(1.5);
  double eta_R0 = 0.75;
  int eta_damping = 1;
  int eta_damping_exp = 1;
  double trK0 = 0.0;
  double chi_floor = 0.1;

#ifdef __APPLE__
  auto t1 = std::chrono::steady_clock::now();
#else
  struct timespec t1, t2, t3, t4, t5, t6;
  GetRdtscTime(&t1);
#endif
 /* H A D   R H S   ( P O I N T - W I S E ) */
  cout << "...calling point-wise had rhs" << endl;
  cal_bssn_rhs_( dtu0[F_ALPHA], dtu0[F_SHIFT1], dtu0[F_SHIFT2], dtu0[F_SHIFT3],
                 dtu0[F_CHI], dtu0[F_TRK],
                 dtu0[F_GT11], dtu0[F_GT12], dtu0[F_GT13],
                 dtu0[F_GT22], dtu0[F_GT23], dtu0[F_GT33],
                 dtu0[F_A11], dtu0[F_A12], dtu0[F_A13],
                 dtu0[F_A22], dtu0[F_A23], dtu0[F_A33],
                 dtu0[F_GAM1], dtu0[F_GAM2], dtu0[F_GAM3],
                 dtu0[F_GB1], dtu0[F_GB2], dtu0[F_GB3],
                 u[F_ALPHA], u[F_SHIFT1], u[F_SHIFT2], u[F_SHIFT3],
                 u[F_CHI], u[F_TRK],
                 u[F_GT11], u[F_GT12], u[F_GT13],
                 u[F_GT22], u[F_GT23], u[F_GT33],
                 u[F_A11], u[F_A12], u[F_A13],
                 u[F_A22], u[F_A23], u[F_A33],
                 u[F_GAM1], u[F_GAM2], u[F_GAM3],
                 u[F_GB1], u[F_GB2], u[F_GB3],
                 lambda, lambda_f,
                 &eta, &eta_damping, &eta_damping_exp, &eta_R0,
                 &trK0, &chi_floor,
                 xi[0], xi[1], xi[2], shp, (shp+1), (shp+2));
#ifdef __APPLE__
  auto t2 = std::chrono::steady_clock::now();
#else
  GetRdtscTime(&t2);
#endif
 /* D E N D R O   R H S */
  cout << "...calling DENDRO rhs" << endl;
  rhs(           dtu1[F_ALPHA], dtu1[F_SHIFT1], dtu1[F_SHIFT2], dtu1[F_SHIFT3],
                 dtu1[F_CHI], dtu1[F_TRK],
                 dtu1[F_GT11], dtu1[F_GT12], dtu1[F_GT13],
                 dtu1[F_GT22], dtu1[F_GT23], dtu1[F_GT33],
                 dtu1[F_A11], dtu1[F_A12], dtu1[F_A13],
                 dtu1[F_A22], dtu1[F_A23], dtu1[F_A33],
                 dtu1[F_GAM1], dtu1[F_GAM2], dtu1[F_GAM3],
                 dtu1[F_GB1], dtu1[F_GB2], dtu1[F_GB3],
                 u[F_ALPHA], u[F_SHIFT1], u[F_SHIFT2], u[F_SHIFT3],
                 u[F_CHI], u[F_TRK],
                 u[F_GT11], u[F_GT12], u[F_GT13],
                 u[F_GT22], u[F_GT23], u[F_GT33],
                 u[F_A11], u[F_A12], u[F_A13],
                 u[F_A22], u[F_A23], u[F_A33],
                 u[F_GAM1], u[F_GAM2], u[F_GAM3],
                 u[F_GB1], u[F_GB2], u[F_GB3],
                 lambda, lambda_f,
                 eta, eta_damping, eta_damping_exp, eta_R0,
                 trK0, chi_floor,
                 xi, shp);
#ifdef __APPLE__
  auto t3 = std::chrono::steady_clock::now();
# else
  GetRdtscTime(&t3);
#endif
  bool do_transpose = false;
  if (do_transpose) {
    cout << "...calling had_transpose rhs" << endl;
    hadrhs_transpose(dtu2[F_ALPHA], dtu2[F_SHIFT1], dtu2[F_SHIFT2], dtu2[F_SHIFT3],
                     dtu2[F_CHI], dtu2[F_TRK],
                     dtu2[F_GT11], dtu2[F_GT12], dtu2[F_GT13],
                     dtu2[F_GT22], dtu2[F_GT23], dtu2[F_GT33],
                     dtu2[F_A11], dtu2[F_A12], dtu2[F_A13],
                     dtu2[F_A22], dtu2[F_A23], dtu2[F_A33],
                     dtu2[F_GAM1], dtu2[F_GAM2], dtu2[F_GAM3],
                     dtu2[F_GB1], dtu2[F_GB2], dtu2[F_GB3],
                     u[F_ALPHA], u[F_SHIFT1], u[F_SHIFT2], u[F_SHIFT3],
                     u[F_CHI], u[F_TRK],
                     u[F_GT11], u[F_GT12], u[F_GT13],
                     u[F_GT22], u[F_GT23], u[F_GT33],
                     u[F_A11], u[F_A12], u[F_A13],
                     u[F_A22], u[F_A23], u[F_A33],
                     u[F_GAM1], u[F_GAM2], u[F_GAM3],
                     u[F_GB1], u[F_GB2], u[F_GB3],
                     lambda, lambda_f,
                     eta, eta_damping, eta_damping_exp, eta_R0,
                     trK0, chi_floor,
                     xi, shp);
  }

#ifdef __APPLE__
  auto t4 = std::chrono::steady_clock::now();
# else
  GetRdtscTime(&t4);
#endif
  cout << "...calling hadrhs (derivatives stored in arrays)" << endl;
  hadrhs(dtu3[F_ALPHA], dtu3[F_SHIFT1], dtu3[F_SHIFT2], dtu3[F_SHIFT3],
                 dtu3[F_CHI], dtu3[F_TRK],
                 dtu3[F_GT11], dtu3[F_GT12], dtu3[F_GT13],
                 dtu3[F_GT22], dtu3[F_GT23], dtu3[F_GT33],
                 dtu3[F_A11], dtu3[F_A12], dtu3[F_A13],
                 dtu3[F_A22], dtu3[F_A23], dtu3[F_A33],
                 dtu3[F_GAM1], dtu3[F_GAM2], dtu3[F_GAM3],
                 dtu3[F_GB1], dtu3[F_GB2], dtu3[F_GB3],
                 u[F_ALPHA], u[F_SHIFT1], u[F_SHIFT2], u[F_SHIFT3],
                 u[F_CHI], u[F_TRK],
                 u[F_GT11], u[F_GT12], u[F_GT13],
                 u[F_GT22], u[F_GT23], u[F_GT33],
                 u[F_A11], u[F_A12], u[F_A13],
                 u[F_A22], u[F_A23], u[F_A33],
                 u[F_GAM1], u[F_GAM2], u[F_GAM3],
                 u[F_GB1], u[F_GB2], u[F_GB3],
                 lambda, lambda_f,
                 eta, eta_damping, eta_damping_exp, eta_R0,
                 trK0, chi_floor,
                 xi, shp);

#ifdef __APPLE__
  auto t5 = std::chrono::steady_clock::now();
# else
  GetRdtscTime(&t5);
#endif
 /*   D E N D R O   V E C T O R I Z E D   R H S   */
  cout << "...calling vectorized DENDRO vec" << endl;
  rhs_vec(       dtu2[F_ALPHA], dtu2[F_SHIFT1], dtu2[F_SHIFT2], dtu2[F_SHIFT3],
                 dtu2[F_CHI], dtu2[F_TRK],
                 dtu2[F_GT11], dtu2[F_GT12], dtu2[F_GT13],
                 dtu2[F_GT22], dtu2[F_GT23], dtu2[F_GT33],
                 dtu2[F_A11], dtu2[F_A12], dtu2[F_A13],
                 dtu2[F_A22], dtu2[F_A23], dtu2[F_A33],
                 dtu2[F_GAM1], dtu2[F_GAM2], dtu2[F_GAM3],
                 dtu2[F_GB1], dtu2[F_GB2], dtu2[F_GB3],
                 u[F_ALPHA], u[F_SHIFT1], u[F_SHIFT2], u[F_SHIFT3],
                 u[F_CHI], u[F_TRK],
                 u[F_GT11], u[F_GT12], u[F_GT13],
                 u[F_GT22], u[F_GT23], u[F_GT33],
                 u[F_A11], u[F_A12], u[F_A13],
                 u[F_A22], u[F_A23], u[F_A33],
                 u[F_GAM1], u[F_GAM2], u[F_GAM3],
                 u[F_GB1], u[F_GB2], u[F_GB3],
                 lambda, lambda_f,
                 eta, eta_damping, eta_damping_exp, eta_R0,
                 trK0, chi_floor,
                 xi, shp);

#ifdef __APPLE__
  auto t6 = std::chrono::steady_clock::now();
# else
  GetRdtscTime(&t6);
#endif

 /* bounds = 1  : include boundaries in analysis of differences
  * bounds = 0  : no boundaries
  */
  int bounds = 0;
  double dn01[NU], dn02[NU], dn12[NU], dn03[NU], dn04[NU];
  cal_diffs(drhs, dtu0, dtu1, shp, bounds);
  cal_norms(dn01, drhs, shp, bounds);

  cal_diffs(drhs, dtu0, dtu2, shp, bounds);
  cal_norms(dn02, drhs, shp, bounds);

  cal_diffs(drhs, dtu1, dtu2, shp, bounds);
  cal_norms(dn12, drhs, shp, bounds);

  cal_diffs(drhs, dtu0, dtu3, shp, bounds);
  cal_norms(dn03, drhs, shp, bounds);

 /* compare had and Dendro RHS */
  std::cout << std::endl << std::endl;
  std::cout << "L2 Norms of had (pt) and Dendro RHS differences" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  for (int m = 0; m < NU; m++) {
    std::cout << "  || diff " << m << " || = " << dn01[m]
              << ", || rhs had || = " << L2norm(dtu0[m],shp)
              << ", || rhs Dendro || = " << L2norm(dtu1[m],shp)
              << std::endl;
  }

 /* compare Dendro and optimized Dendro RHS */
  std::cout << std::endl << std::endl;
  std::cout << "L2 Norms of Dendro and vectorized Dendro RHS differences" << std::endl;
  std::cout << "--------------------------------------------------------" << std::endl;
  for (int m = 0; m < NU; m++) {
    std::cout << "  || diff " << m << " || = " << dn12[m]
              << ", || rhs || = " << L2norm(dtu1[m],shp)
              << ", || rhs opt || = " << L2norm(dtu2[m],shp)
              << std::endl;

  }

 /* compare Dendro and optimized Dendro RHS */
  std::cout << std::endl << std::endl;
  std::cout << "L2 Norms of  had (pt) and V E C T O R I Z E D  Dendro RHS differences" << std::endl;
  std::cout << "---------------------------------------------------------------------" << std::endl;
  for (int m = 0; m < NU; m++) {
    std::cout << "  || diff " << m << " || = " << dn02[m]
              << ", || rhs || = " << L2norm(dtu0[m],shp)
              << ", || rhs opt || = " << L2norm(dtu2[m],shp)
              << std::endl;

  }

 /* compare Dendro and optimized Dendro RHS */
  std::cout << std::endl << std::endl;
  std::cout << "L2 Norms of  A R R A Y  H A D  RHS differences" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for (int m = 0; m < NU; m++) {
    std::cout << "  || diff " << m << " || = " << dn03[m]
              << ", || rhs || = " << L2norm(dtu0[m],shp)
              << ", || rhs opt || = " << L2norm(dtu3[m],shp)
              << std::endl;

  }

 /*---------------------------------------------------------------
  *  Timing information
  *---------------------------------------------------------------*/
  std::cout << std::endl << std::endl;
  std::cout << "Timing Information" << std::endl;
  std::cout << "------------------" << std::endl;

#ifdef __APPLE__
  auto dt1 = t2 - t1;
  auto dt2 = t3 - t2;
  auto dt3 = t4 - t3;
  auto dt4 = t5 - t4;
  auto dt5 = t6 - t5;

  std::cout << "Fortran (pt derivs): " << std::chrono::duration <double, std::milli> (dt1).count() << " ms" << std::endl;
  std::cout << "Fortran (3D derivs/transpose): " << std::chrono::duration <double, std::milli> (dt3).count() << " ms" << std::endl;
  std::cout << "Fortran (3D derivs): " << std::chrono::duration <double, std::milli> (dt4).count() << " ms" << std::endl;
  std::cout << "Dendro:  " << std::chrono::duration <double, std::milli> (dt2).count() << " ms" << std::endl;
  std::cout << "Dendro (opt, 3D derivs): " << std::chrono::duration <double, std::milli> (dt5).count() << " ms" << std::endl;
#else
  std::cout << "Fortran (pt derivs): " << TimeSpecDiff(&t2,&t1)->tv_sec << "." << TimeSpecDiff(&t2,&t1)->tv_nsec << std::endl;
  std::cout << "Dendro: " << TimeSpecDiff(&t3,&t2)->tv_sec << "." << TimeSpecDiff(&t3,&t2)->tv_nsec << std::endl;
  std::cout << "Dendro (opt): " << TimeSpecDiff(&t6,&t5)->tv_sec << "." << TimeSpecDiff(&t6,&t5)->tv_nsec << std::endl;
  std::cout << "Fortran (3D derivs/transpose): " << TimeSpecDiff(&t4,&t3)->tv_sec << "." << TimeSpecDiff(&t4,&t3)->tv_nsec << std::endl;
  std::cout << "Fortran (3D derivs): " << TimeSpecDiff(&t5,&t4)->tv_sec << "." << TimeSpecDiff(&t5,&t4)->tv_nsec << std::endl;
#endif

  dollar::text(std::cout);
  dollar::clear();

  //printdiff(dtu, dth, shp);

  //sdfoutput("gt11rhs", shp, xi[0], dtu1[F_GT11]);
  //sdfoutput("gt11rhs_had", shp, xi[0], dtu0[F_GT11]);

 /* free memory */
  for (int i = 0; i < NU; i++) {
    delete [] u[i];
    delete [] dtu0[i];
    delete [] dtu1[i];
    delete [] dtu2[i];
    delete [] drhs[i];
  }

  std::cout << std::endl << "Finis." << std::endl;

  return 0;
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void define_coordinates(double *x, double xmin, double xmax, double &h, int n)
{
  h = (xmax - xmin) / (n - 1);
  for (int i = 0; i < n; i++) {
    x[i] = xmin + i*h;
  }
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void initial_data(double *u[], double *xi[], int shp[])
{

  const double pi = acos(-1.0);
  const double f1 = 31.0/17.0;
  const double f2 = 37.0/11.0;

  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];

  for (int k = 0; k < nz; k++) {
    double z = xi[2][k];
    for (int j = 0; j < ny; j++) {
      double y = xi[1][j];
      for (int i = 0; i < nx; i++) {
        double x = xi[0][i];
        int pp = i + nx*(j + ny*k);

        u[F_ALPHA][pp] = 1.0 - 0.25*sin(f1*x);
        //u[F_ALPHA][pp] = 1.0;
        u[F_SHIFT1][pp] = 4.0/17.0*sin(x)*cos(z);
        u[F_SHIFT2][pp] = pi/5.0*cos(y)*sin(z+x);
        u[F_SHIFT3][pp] = 4.0/17.0*sin(f2*x)*sin(z);

        u[F_GB1][pp] = 31.0*x*cos(f1*z+y);
        u[F_GB2][pp] = 7.0*y*sin(f1*x+y) + 3.0*cos(z);
        u[F_GB3][pp] = 5.0*z*cos(f1*x+y) + 7.0*sin(z+y+x) + 1.0;

        u[F_GAM1][pp] = 5.0*cos(x)/(10.0*sin(x+z)+26.0-1.0*cos(x*z)*cos(x));
        u[F_GAM2][pp] = -5.0*sin(y)/(25.0+10.0*cos(y+z)+cos(y)*cos(y*z));
        u[F_GAM3][pp] = -5.0*sin(z)/(25.0+10.0*cos(y+x)+cos(y*x)*cos(z));

        u[F_CHI][pp] = 1.0 + exp(-4.0*cos(x)*sin(y));
        //u[F_CHI][pp] = 2.0;

        u[F_GT11][pp] = 1.00+0.2*sin(x+z)*cos(y);
        u[F_GT22][pp] = 1.00+0.2*cos(y)*cos(z+ x);
        u[F_GT33][pp] = 1.00 / ( u[F_GT11][pp] +  u[F_GT22][pp] );
        u[F_GT12][pp] = 0.7*cos(x*x + y*y);
        u[F_GT13][pp] = 0.3*sin(z)*cos(x);
        u[F_GT23][pp] = -0.5*sin(x*x)*cos(y)*cos(z);

        u[F_TRK][pp] = 5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)
               +5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)
               +0.4*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))
                 *exp(-4.0*cos(x)*sin(y))*cos(z);

        u[F_A11][pp] = exp(-4.0*cos(x)*sin(y))*(cos(x) -0.3333333333*exp(4.0*cos(x)*sin(y)) *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y) +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
        u[F_A12][pp] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
        u[F_A13][pp] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);

        u[F_A22][pp] = exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
        u[F_A23][pp] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);

        u[F_A33][pp] = exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));


      }
    }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
double grad(int dir, double *u, double h, int idx[], int shp[])
{

  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];

  const int i = idx[0];
  const int j = idx[1];
  const int k = idx[2];

  assert(i > 2);
  assert(i < nx-2);
  assert(j > 2);
  assert(j < ny-2);
  assert(k > 2);
  assert(k < nz-2);

  int pp = i + nx*(j + ny*k);
  double du;

  switch(dir) {
    case 0:
      du = (u[pp-2] - 8.0 * u[pp-1] + 8.0 * u[pp+1] - u[pp+2]) / (12.0*h);
      break;
    case 1:
      du = (u[pp-2*nx] - 8.0 * u[pp-1*nx] + 8.0 * u[pp+1*nx]
                       - u[pp+2*nx]) / (12.0*h);
      break;
    case 2:
      du = (u[pp-2*nx*ny] - 8.0 * u[pp-1*nx*ny] + 8.0 * u[pp+1*nx*ny]
                       - u[pp+2*nx*ny]) / (12.0*h);
      break;
  }

  return du;
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
double grad2(int dir1, int dir2, double *u, double h, int idx[], int shp[])
{
  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];

  const int i = idx[0];
  const int j = idx[1];
  const int k = idx[2];

  assert(i > 2);
  assert(i < nx-2);
  assert(j > 2);
  assert(j < ny-2);
  assert(k > 2);
  assert(k < nz-2);

  int pp = i + nx * (j + ny*k);
  double ddu;


  if (dir1 > dir2) {
    int t = dir1;
    dir1 = dir2;
    dir2 = t;
  }

  if (dir1 == dir2) {
    switch (dir1) {
    case 0:
      ddu = ( -u[pp-2] + 16.0 * u[pp-1] - 30.0 * u[pp]
                       + 16.0 * u[pp+1] - u[pp+2] ) / (12.0*h*h);
      break;
    case 1:
      ddu = ( -u[pp-2*nx] + 16.0 * u[pp-1*nx] - 30.0 * u[pp]
                       + 16.0 * u[pp+1*nx] - u[pp+2*nx] ) / (12.0*h*h);
      break;
    case 2:
      ddu = ( -u[pp-2*nx*ny] + 16.0 * u[pp-1*nx*ny] - 30.0 * u[pp]
                       + 16.0 * u[pp+1*nx*ny] - u[pp+2*nx*ny] ) / (12.0*h*h);
      break;
    }
  }
  else {
    if (dir1 == 0 && dir2 == 1) {
      ddu = (
         u[IDX(i-2,j-2,k)] * (0.00694444444444444)
       + u[IDX(i-2,j-1,k)] * (-0.0555555555555556)
       + u[IDX(i-2,j+1,k)] * (0.0555555555555556)
       + u[IDX(i-2,j+2,k)] * (-0.00694444444444444)
       + u[IDX(i-1,j-2,k)] * (-0.0555555555555556)
       + u[IDX(i-1,j-1,k)] * (0.444444444444444)
       + u[IDX(i-1,j+1,k)] * (-0.444444444444444)
       + u[IDX(i-1,j+2,k)] * (0.0555555555555556)
       + u[IDX(i+1,j-2,k)] * (0.0555555555555556)
       + u[IDX(i+1,j-1,k)] * (-0.444444444444444)
       + u[IDX(i+1,j+1,k)] * (0.444444444444444)
       + u[IDX(i+1,j+2,k)] * (-0.0555555555555556)
       + u[IDX(i+2,j-2,k)] * (-0.00694444444444444)
       + u[IDX(i+2,j-1,k)] * (0.0555555555555556)
       + u[IDX(i+2,j+1,k)] * (-0.0555555555555556)
       + u[IDX(i+2,j+2,k)] * (0.00694444444444444)
       ) / (h*h);
    }
    else if (dir1 == 0 && dir2 == 2) {
      ddu = (
         u[IDX(i-2,j,k-2)] * (0.00694444444444444)
       + u[IDX(i-2,j,k-1)] * (-0.0555555555555556)
       + u[IDX(i-2,j,k+1)] * (0.0555555555555556)
       + u[IDX(i-2,j,k+2)] * (-0.00694444444444444)
       + u[IDX(i-1,j,k-2)] * (-0.0555555555555556)
       + u[IDX(i-1,j,k-1)] * (0.444444444444444)
       + u[IDX(i-1,j,k+1)] * (-0.444444444444444)
       + u[IDX(i-1,j,k+2)] * (0.0555555555555556)
       + u[IDX(i+1,j,k-2)] * (0.0555555555555556)
       + u[IDX(i+1,j,k-1)] * (-0.444444444444444)
       + u[IDX(i+1,j,k+1)] * (0.444444444444444)
       + u[IDX(i+1,j,k+2)] * (-0.0555555555555556)
       + u[IDX(i+2,j,k-2)] * (-0.00694444444444444)
       + u[IDX(i+2,j,k-1)] * (0.0555555555555556)
       + u[IDX(i+2,j,k+1)] * (-0.0555555555555556)
       + u[IDX(i+2,j,k+2)] * (0.00694444444444444)
       ) / (h*h);
    }
    else if (dir1 == 1 && dir2 == 2) {
      ddu = (
         u[IDX(i,j-2,k-2)] * (0.00694444444444444)
       + u[IDX(i,j-2,k-1)] * (-0.0555555555555556)
       + u[IDX(i,j-2,k+1)] * (0.0555555555555556)
       + u[IDX(i,j-2,k+2)] * (-0.00694444444444444)
       + u[IDX(i,j-1,k-2)] * (-0.0555555555555556)
       + u[IDX(i,j-1,k-1)] * (0.444444444444444)
       + u[IDX(i,j-1,k+1)] * (-0.444444444444444)
       + u[IDX(i,j-1,k+2)] * (0.0555555555555556)
       + u[IDX(i,j+1,k-2)] * (0.0555555555555556)
       + u[IDX(i,j+1,k-1)] * (-0.444444444444444)
       + u[IDX(i,j+1,k+1)] * (0.444444444444444)
       + u[IDX(i,j+1,k+2)] * (-0.0555555555555556)
       + u[IDX(i,j+2,k-2)] * (-0.00694444444444444)
       + u[IDX(i,j+2,k-1)] * (0.0555555555555556)
       + u[IDX(i,j+2,k+1)] * (-0.0555555555555556)
       + u[IDX(i,j+2,k+2)] * (0.00694444444444444)
       ) / (h*h);
    }
  }

  return ddu;
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
double agrad(int dir, double *beta, double *u, double h, int idx[], int shp[])
{

  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];

  const int i = idx[0];
  const int j = idx[1];
  const int k = idx[2];

  assert(i > 2);
  assert(i < nx-3);
  assert(j > 2);
  assert(j < ny-3);
  assert(k > 2);
  assert(k < nz-3);

  int pp = i + nx*(j + ny*k);
  double du;

  switch(dir) {
    case 0:
      if (beta[pp] >= 0.0) {
        du = (
                  u[pp-1]*(-0.25)
                + u[pp  ]*(-0.833333333333333)
                + u[pp+1]*(1.5)
                + u[pp+2]*(-0.5)
                + u[pp+3]*(0.0833333333333333)
             ) / h;
      }
      else {
        du = (
                  u[pp-3]*(-0.0833333333333333)
                + u[pp-2]*(0.5)
                + u[pp-1]*(-1.5)
                + u[pp  ]*(0.833333333333333)
                + u[pp+1]*(0.25)
             ) / h;
      }

      break;
    case 1:
      if (beta[pp] >= 0.0) {
        du = (
                  u[pp-1*nx]*(-0.25)
                + u[pp     ]*(-0.833333333333333)
                + u[pp+1*nx]*(1.5)
                + u[pp+2*nx]*(-0.5)
                + u[pp+3*nx]*(0.0833333333333333)
             ) / h;
      }
      else {
        du = (
                  u[pp-3*nx]*(-0.0833333333333333)
                + u[pp-2*nx]*(0.5)
                + u[pp-1*nx]*(-1.5)
                + u[pp     ]*(0.833333333333333)
                + u[pp+1*nx]*(0.25)
             ) / h;
      }
      break;
    case 2:
      if (beta[pp] >= 0.0) {
        du = (
                  u[pp-1*nx*ny]*(-0.25)
                + u[pp        ]*(-0.833333333333333)
                + u[pp+1*nx*ny]*(1.5)
                + u[pp+2*nx*ny]*(-0.5)
                + u[pp+3*nx*ny]*(0.0833333333333333)
             ) / h;
      }
      else {
        du = (
                  u[pp-3*nx*ny]*(-0.0833333333333333)
                + u[pp-2*nx*ny]*(0.5)
                + u[pp-1*nx*ny]*(-1.5)
                + u[pp        ]*(0.833333333333333)
                + u[pp+1*nx*ny]*(0.25)
             ) / h;
      }


      break;
  }

  return du;
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void cal_norms(double df[], double *f[], int shp[], int bounds)
{

  for (int m = 0; m < NU; m++) {
    df[m] = L2norm(f[m], shp);
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void cal_diffs(double *df[], double *f1[], double *f2[], int shp[], int bounds)
{

  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];
  const int w = 3;
  const int npts = nx * ny * nz;

  for (int m = 0; m < NU; m++) {
    for (int i = 0; i < npts; i++) {
      df[m][i] = fabs(f1[m][i] - f2[m][i]);
    }
  }

  if (bounds != 1) {

   /* zero xmin */
    for (int m = 0; m < NU; m++) {

      for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
          for (int i = 0; i < w; i++) {
            int pp = i + nx * (j + ny * k);
            df[m][pp] = 0.0;
          }
        }
      }

     /* zero xmax */
      for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
          for (int i = nx - w; i < nx; i++) {
            int pp = i + nx * (j + ny * k);
            df[m][pp] = 0.0;
          }
        }
      }

     /* zero ymin */
      for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
          for (int j = 0; j < w; j++) {
            int pp = i + nx * (j + ny * k);
            df[m][pp] = 0.0;
          }
        }
      }

     /* zero ymax */
      for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
          for (int j = ny - w; j < ny; j++) {
            int pp = i + nx * (j + ny * k);
            df[m][pp] = 0.0;
          }
        }
      }

     /* zero zmin */
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
          for (int k = 0; k < w; k++) {
            int pp = i + nx * (j + ny * k);
            df[m][pp] = 0.0;
          }
        }
      }

     /* zero zmax */
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
          for (int k = nz - w; k < nz; k++) {
            int pp = i + nx * (j + ny * k);
            df[m][pp] = 0.0;
          }
        }
      }

    }  // m loop
  }    // if (!bounds)

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void set_zero(double *u, int shp[])
{
  int n = shp[0] * shp[1] * shp[2];
  for (int i = 0; i < n; i++) {
    u[i] = 0.0;
  }
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
double L2norm(double * const u, const int shp[])
{

  const int n = shp[0] * shp[1] * shp[2];
  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s += u[i]*u[i];
  }
  return sqrt(s/n);

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void sdfoutput(char *fname, int *shp, double *coords, double *data)
{

  char *cname = "x|y|z";
  int rank = 3;
  double time = 0.0;

  // int rc = gft_out_full(fname, time, shp, cname, rank, coords, data);

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void printdiff(double **f1, double **f2, int *shp)
{
  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];

  for (int k = 0; k < nz; k++) {
  for (int j = 0; j < ny; j++) {
  for (int i = 0; i < nx; i++) {
    int pp = i + nx * (j + ny * k);
    double diff = fabs(f1[0][pp] - f2[0][pp]);
    if (diff > 1.0e-4) {
      std::cout << "... gtxx.  Diff = " << diff <<
       ",   i=" << i << ", j=" << j << ", k =" << k << std::endl;
    }
  }
  }
  }

}


#ifndef __APPLE__
/* returns a static buffer of struct timespec with the time difference of ts1 and ts2
   ts1 is assumed to be greater than ts2 */
struct timespec *TimeSpecDiff(struct timespec *ts1, struct timespec *ts2)
{
  static struct timespec ts;
  ts.tv_sec = ts1->tv_sec - ts2->tv_sec;
  ts.tv_nsec = ts1->tv_nsec - ts2->tv_nsec;
  if (ts.tv_nsec < 0) {
    ts.tv_sec--;
    ts.tv_nsec += NANO_SECONDS_IN_SEC;
  }
  return &ts;
}

double g_TicksPerNanoSec;
static void CalibrateTicks()
{
  struct timespec begints, endts;
  uint64_t begin = 0, end = 0;
  clock_gettime(CLOCK_MONOTONIC, &begints);
  begin = RDTSC();
  uint64_t i;
  for (i = 0; i < 1000000; i++); /* must be CPU intensive */
  end = RDTSC();
  clock_gettime(CLOCK_MONOTONIC, &endts);
  struct timespec *tmpts = TimeSpecDiff(&endts, &begints);
  uint64_t nsecElapsed = tmpts->tv_sec * 1000000000LL + tmpts->tv_nsec;
  g_TicksPerNanoSec = (double)(end - begin)/(double)nsecElapsed;
}

/* Call once before using RDTSC, has side effect of binding process to CPU1 */
void InitRdtsc()
{
  cpu_set_t cpuMask;
 // cpuMask = 2; // bind to cpu 1
  CPU_ZERO(&cpuMask);
  CPU_SET(1, &cpuMask);
  sched_setaffinity(0, sizeof(cpuMask), &cpuMask);
  CalibrateTicks();
}

void GetTimeSpec(struct timespec *ts, uint64_t nsecs)
{
  ts->tv_sec = nsecs / NANO_SECONDS_IN_SEC;
  ts->tv_nsec = nsecs % NANO_SECONDS_IN_SEC;
}

/* ts will be filled with time converted from TSC reading */
void GetRdtscTime(struct timespec *ts)
{
  GetTimeSpec(ts, RDTSC() / g_TicksPerNanoSec);
}

#endif
