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

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void hadrhs( double *  dt_alpha, double *  dt_shiftx,
             double *  dt_shifty, double *  dt_shiftz,
             double *  dt_chi, double *  dt_trK,
             double *  dt_gtxx, double *  dt_gtxy,
             double *  dt_gtxz, double *  dt_gtyy,
             double *  dt_gtyz, double *  dt_gtzz,
             double *  dt_Atxx, double *  dt_Atxy,
             double *  dt_Atxz, double *  dt_Atyy,
             double *  dt_Atyz, double *  dt_Atzz,
             double *  dt_Gamtx, double *  dt_Gamty,
             double *  dt_Gamtz, double *  dt_gbx,
             double *  dt_gby, double *  dt_gbz,
             double *  alpha,  double *  shiftx,
             double *  shifty,  double *  shiftz,
             double *  chi,  double *  trK,
             double *  gtxx,  double *  gtxy,
             double *  gtxz,  double *  gtyy,
             double *  gtyz,  double *  gtzz,
             double *  Atxx,  double *  Atxy,
             double *  Atxz,  double *  Atyy,
             double *  Atyz,  double *  Atzz,
             double *  Gamtx,  double *  Gamty,
             double *  Gamtz,  double *  gbx,
             double *  gby,  double *  gbz,
             int *  lambda,  int *  lambda_f,
             double eta_const,  int eta_damping,
             int eta_damping_exp,  double eta_R0,
             double trK0,  double chi_floor,
             double *  xi[],  int shp[] )
{ $

  const bool ltrace = false;

  int nx = shp[0];
  int ny = shp[1];
  int nz = shp[2];
  int nd = nx * ny * nz;

  double * x1d = xi[0];
  double * y1d = xi[1];
  double * z1d = xi[2];

  const double dx = x1d[1] - x1d[0];
  const double dy = y1d[1] - y1d[0];
  const double dz = z1d[1] - z1d[0];


#ifdef __APPLE__
  auto t1 = std::chrono::steady_clock::now();
#else
  struct timespec t1, t2, t3, t4, t5;
  GetRdtscTime(&t1);
#endif

#include "NewMR_alloc.h"

#include "MR_derivs.h"

#ifdef __APPLE__
  auto t2 = std::chrono::steady_clock::now();
#else
  GetRdtscTime(&t2);
#endif


  if (ltrace) std::cout << "calling cal_had_bssn_rhs" << std::endl;

  cal_had_bssn_rhs_(
       dt_alpha, dt_shiftx, dt_shifty, dt_shiftz, dt_chi,
       dt_trK, dt_gtxx, dt_gtxy, dt_gtxz, dt_gtyy, dt_gtyz,
       dt_gtzz, dt_Atxx, dt_Atxy, dt_Atxz, dt_Atyy, dt_Atyz,
       dt_Atzz, dt_Gamtx, dt_Gamty, dt_Gamtz, dt_gbx, dt_gby, dt_gbz,
       alpha, shiftx, shifty, shiftz, chi, trK,
       gtxx, gtxy, gtxz, gtyy, gtyz, gtzz,
       Atxx, Atxy, Atxz, Atyy, Atyz,
       Atzz, Gamtx, Gamty, Gamtz, gbx, gby, gbz,
       dx_alpha,dy_alpha,dz_alpha,dx_shiftx,
       dy_shiftx,dz_shiftx,dx_shifty,dy_shifty,
       dz_shifty,dx_shiftz,dy_shiftz,dz_shiftz,
       dx_gbx,dy_gbx,dz_gbx,dx_gby,
       dy_gby,dz_gby,dx_gbz,dy_gbz,
       dz_gbz,dx_chi,dy_chi,dz_chi,
       dx_Gamtx,dy_Gamtx,dz_Gamtx,dx_Gamty,
       dy_Gamty,dz_Gamty,dx_Gamtz,dy_Gamtz,
       dz_Gamtz,dx_trK,dy_trK,dz_trK,
       dx_gtxx,dy_gtxx,dz_gtxx,dx_gtxy,
       dy_gtxy,dz_gtxy,dx_gtxz,dy_gtxz,
       dz_gtxz,dx_gtyy,dy_gtyy,dz_gtyy,
       dx_gtyz,dy_gtyz,dz_gtyz,dx_gtzz,
       dy_gtzz,dz_gtzz,dx_Atxx,dy_Atxx,
       dz_Atxx,dx_Atxy,dy_Atxy,dz_Atxy,
       dx_Atxz,dy_Atxz,dz_Atxz,dx_Atyy,
       dy_Atyy,dz_Atyy,dx_Atyz,dy_Atyz,
       dz_Atyz,dx_Atzz,dy_Atzz,dz_Atzz,
       dxx_gtxx,dxy_gtxx,dxz_gtxx,dyy_gtxx,
       dyz_gtxx,dzz_gtxx,dxx_gtxy,dxy_gtxy,
       dxz_gtxy,dyy_gtxy,dyz_gtxy,dzz_gtxy,
       dxx_gtxz,dxy_gtxz,dxz_gtxz,dyy_gtxz,
       dyz_gtxz,dzz_gtxz,dxx_gtyy,dxy_gtyy,
       dxz_gtyy,dyy_gtyy,dyz_gtyy,dzz_gtyy,
       dxx_gtyz,dxy_gtyz,dxz_gtyz,dyy_gtyz,
       dyz_gtyz,dzz_gtyz,dxx_gtzz,dxy_gtzz,
       dxz_gtzz,dyy_gtzz,dyz_gtzz,dzz_gtzz,
       dxx_chi,dxy_chi,dxz_chi,dyy_chi,
       dyz_chi,dzz_chi,dxx_alpha,dxy_alpha,
       dxz_alpha,dyy_alpha,dyz_alpha,dzz_alpha,
       dxx_shiftx,dxy_shiftx,dxz_shiftx,dyy_shiftx,
       dyz_shiftx,dzz_shiftx,dxx_shifty,dxy_shifty,
       dxz_shifty,dyy_shifty,dyz_shifty,dzz_shifty,
       dxx_shiftz,dxy_shiftz,dxz_shiftz,dyy_shiftz,
       dyz_shiftz,dzz_shiftz,
       lambda, lambda_f,
       &eta_const, &eta_damping, &eta_damping_exp, &eta_R0,
       &trK0, &chi_floor,
       x1d, y1d, z1d, &nx, &ny, &nz);

#ifdef __APPLE__
  auto t3 = std::chrono::steady_clock::now();
#else
  GetRdtscTime(&t3);
#endif


#include "NewMR_dealloc.h"

#include "NewMR_alloc_adv.h"

#include "MR_derivs_adv.h"
#ifdef __APPLE__
  auto t4 = std::chrono::steady_clock::now();
#else
  GetRdtscTime(&t4);
#endif



  cal_had_bssn_rhs_lie_(
        dt_alpha, dt_shiftx, dt_shifty, dt_shiftz, dt_chi,
        dt_trK, dt_gtxx, dt_gtxy, dt_gtxz, dt_gtyy, dt_gtyz,
        dt_gtzz, dt_Atxx, dt_Atxy, dt_Atxz, dt_Atyy, dt_Atyz,
        dt_Atzz, dt_Gamtx, dt_Gamty, dt_Gamtz, dt_gbx, dt_gby, dt_gbz,
        shiftx, shifty, shiftz,
        adx_gtxx,ady_gtxx,adz_gtxx,adx_gtxy,
        ady_gtxy,adz_gtxy,adx_gtxz,ady_gtxz,
        adz_gtxz,adx_gtyy,ady_gtyy,adz_gtyy,
        adx_gtyz,ady_gtyz,adz_gtyz,adx_gtzz,
        ady_gtzz,adz_gtzz,adx_Atxx,ady_Atxx,
        adz_Atxx,adx_Atxy,ady_Atxy,adz_Atxy,
        adx_Atxz,ady_Atxz,adz_Atxz,adx_Atyy,
        ady_Atyy,adz_Atyy,adx_Atyz,ady_Atyz,
        adz_Atyz,adx_Atzz,ady_Atzz,adz_Atzz,
        adx_alpha,ady_alpha,adz_alpha,adx_shiftx,
        ady_shiftx,adz_shiftx,adx_shifty,ady_shifty,
        adz_shifty,adx_shiftz,ady_shiftz,adz_shiftz,
        adx_chi,ady_chi,adz_chi,adx_Gamtx,
        ady_Gamtx,adz_Gamtx,adx_Gamty,ady_Gamty,
        adz_Gamty,adx_Gamtz,ady_Gamtz,adz_Gamtz,
        adx_trK,ady_trK,adz_trK,adx_gbx,
        ady_gbx,adz_gbx,adx_gby,ady_gby,
        adz_gby,adx_gbz,ady_gbz,adz_gbz,
        lambda, lambda_f,
        x1d, y1d, z1d, &nx, &ny, &nz);

#ifdef __APPLE__
  auto t5 = std::chrono::steady_clock::now();
#else
  GetRdtscTime(&t5);
#endif


#include "NewMR_dealloc_adv.h"

#ifdef __APPLE__
auto dt1 = t2 - t1;
auto dt2 = t3 - t2;
auto dt3 = t4 - t3;
auto dt4 = t5 - t4;

auto dtderivs = dt1 + dt3;
auto dtrhs = dt2 + dt4;

  std::cout << "hadrhs: derivs (1)       " << std::chrono::duration <double, std::milli> (dt1).count() << " ms" << std::endl;
  std::cout << "hadrhs: rhs (1) " << std::chrono::duration <double, std::milli> (dt2).count() << " ms" << std::endl;
  std::cout << "hadrhs: adv derivs (2)   " << std::chrono::duration <double, std::milli> (dt3).count() << " ms" << std::endl;
  std::cout << "hadrhs: adv rhs (2)      " << std::chrono::duration <double, std::milli> (dt4).count() << " ms" << std::endl;
  std::cout << "-------------------------------------------------------------" << std::endl;
  std::cout << "hadrhs: total derivs     " << std::chrono::duration <double, std::milli> (dtderivs).count() << " ms" << std::endl;
  std::cout << "hadrhs: total rhs        " << std::chrono::duration <double, std::milli> (dtrhs).count() << " ms" << std::endl;
#else
  std::cout << "hadrhs: derivs (1) " << TimeSpecDiff(&t2,&t1)->tv_sec << "." << TimeSpecDiff(&t2,&t1)->tv_nsec << std::endl;
  std::cout << "hadrhs: rhs (1) " << TimeSpecDiff(&t3,&t2)->tv_sec << "." << TimeSpecDiff(&t3,&t2)->tv_nsec << std::endl;
  std::cout << "hadrhs: adv derivs (2) " << TimeSpecDiff(&t4,&t3)->tv_sec << "." << TimeSpecDiff(&t4,&t3)->tv_nsec << std::endl;
  std::cout << "hadrhs: adv rhs (2) " << TimeSpecDiff(&t5,&t4)->tv_sec << "." << TimeSpecDiff(&t5,&t4)->tv_nsec << std::endl;
#endif






}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void hadrhs_transpose( double *  dt_alpha, double *  dt_shiftx,
             double *  dt_shifty, double *  dt_shiftz,
             double *  dt_chi, double *  dt_trK,
             double *  dt_gtxx, double *  dt_gtxy,
             double *  dt_gtxz, double *  dt_gtyy,
             double *  dt_gtyz, double *  dt_gtzz,
             double *  dt_Atxx, double *  dt_Atxy,
             double *  dt_Atxz, double *  dt_Atyy,
             double *  dt_Atyz, double *  dt_Atzz,
             double *  dt_Gamtx, double *  dt_Gamty,
             double *  dt_Gamtz, double *  dt_gbx,
             double *  dt_gby, double *  dt_gbz,
             double *  alpha,  double *  shiftx,
             double *  shifty,  double *  shiftz,
             double *  chi,  double *  trK,
             double *  gtxx,  double *  gtxy,
             double *  gtxz,  double *  gtyy,
             double *  gtyz,  double *  gtzz,
             double *  Atxx,  double *  Atxy,
             double *  Atxz,  double *  Atyy,
             double *  Atyz,  double *  Atzz,
             double *  Gamtx,  double *  Gamty,
             double *  Gamtz,  double *  gbx,
             double *  gby,  double *  gbz,
             int *  lambda,  int *  lambda_f,
             double eta_const,  int eta_damping,
             int eta_damping_exp,  double eta_R0,
             double trK0,  double chi_floor,
             double *  xi[],  int shp[] )
{ $

  const bool ltrace = false;

  int nx = shp[0];
  int ny = shp[1];
  int nz = shp[2];
  int nd = nx * ny * nz;

  double * x1d = xi[0];
  double * y1d = xi[1];
  double * z1d = xi[2];

  const double dx = x1d[1] - x1d[0];
  const double dy = y1d[1] - y1d[0];
  const double dz = z1d[1] - z1d[0];



#include "NewMR_alloc.h"

#include "NewMR_derivs.h"

  if (ltrace) std::cout << "calling cal_had_bssn_rhs" << std::endl;

  cal_had_bssn_rhs_(
       dt_alpha, dt_shiftx, dt_shifty, dt_shiftz, dt_chi,
       dt_trK, dt_gtxx, dt_gtxy, dt_gtxz, dt_gtyy, dt_gtyz,
       dt_gtzz, dt_Atxx, dt_Atxy, dt_Atxz, dt_Atyy, dt_Atyz,
       dt_Atzz, dt_Gamtx, dt_Gamty, dt_Gamtz, dt_gbx, dt_gby, dt_gbz,
       alpha, shiftx, shifty, shiftz, chi, trK,
       gtxx, gtxy, gtxz, gtyy, gtyz, gtzz,
       Atxx, Atxy, Atxz, Atyy, Atyz,
       Atzz, Gamtx, Gamty, Gamtz, gbx, gby, gbz,
       dx_alpha,dy_alpha,dz_alpha,dx_shiftx,
       dy_shiftx,dz_shiftx,dx_shifty,dy_shifty,
       dz_shifty,dx_shiftz,dy_shiftz,dz_shiftz,
       dx_gbx,dy_gbx,dz_gbx,dx_gby,
       dy_gby,dz_gby,dx_gbz,dy_gbz,
       dz_gbz,dx_chi,dy_chi,dz_chi,
       dx_Gamtx,dy_Gamtx,dz_Gamtx,dx_Gamty,
       dy_Gamty,dz_Gamty,dx_Gamtz,dy_Gamtz,
       dz_Gamtz,dx_trK,dy_trK,dz_trK,
       dx_gtxx,dy_gtxx,dz_gtxx,dx_gtxy,
       dy_gtxy,dz_gtxy,dx_gtxz,dy_gtxz,
       dz_gtxz,dx_gtyy,dy_gtyy,dz_gtyy,
       dx_gtyz,dy_gtyz,dz_gtyz,dx_gtzz,
       dy_gtzz,dz_gtzz,dx_Atxx,dy_Atxx,
       dz_Atxx,dx_Atxy,dy_Atxy,dz_Atxy,
       dx_Atxz,dy_Atxz,dz_Atxz,dx_Atyy,
       dy_Atyy,dz_Atyy,dx_Atyz,dy_Atyz,
       dz_Atyz,dx_Atzz,dy_Atzz,dz_Atzz,
       dxx_gtxx,dxy_gtxx,dxz_gtxx,dyy_gtxx,
       dyz_gtxx,dzz_gtxx,dxx_gtxy,dxy_gtxy,
       dxz_gtxy,dyy_gtxy,dyz_gtxy,dzz_gtxy,
       dxx_gtxz,dxy_gtxz,dxz_gtxz,dyy_gtxz,
       dyz_gtxz,dzz_gtxz,dxx_gtyy,dxy_gtyy,
       dxz_gtyy,dyy_gtyy,dyz_gtyy,dzz_gtyy,
       dxx_gtyz,dxy_gtyz,dxz_gtyz,dyy_gtyz,
       dyz_gtyz,dzz_gtyz,dxx_gtzz,dxy_gtzz,
       dxz_gtzz,dyy_gtzz,dyz_gtzz,dzz_gtzz,
       dxx_chi,dxy_chi,dxz_chi,dyy_chi,
       dyz_chi,dzz_chi,dxx_alpha,dxy_alpha,
       dxz_alpha,dyy_alpha,dyz_alpha,dzz_alpha,
       dxx_shiftx,dxy_shiftx,dxz_shiftx,dyy_shiftx,
       dyz_shiftx,dzz_shiftx,dxx_shifty,dxy_shifty,
       dxz_shifty,dyy_shifty,dyz_shifty,dzz_shifty,
       dxx_shiftz,dxy_shiftz,dxz_shiftz,dyy_shiftz,
       dyz_shiftz,dzz_shiftz,
       lambda, lambda_f,
       &eta_const, &eta_damping, &eta_damping_exp, &eta_R0,
       &trK0, &chi_floor,
       x1d, y1d, z1d, &nx, &ny, &nz);

#include "NewMR_dealloc.h"

#include "NewMR_alloc_adv.h"

#include "NewMR_derivs_adv.h"


  cal_had_bssn_rhs_lie_(
        dt_alpha, dt_shiftx, dt_shifty, dt_shiftz, dt_chi,
        dt_trK, dt_gtxx, dt_gtxy, dt_gtxz, dt_gtyy, dt_gtyz,
        dt_gtzz, dt_Atxx, dt_Atxy, dt_Atxz, dt_Atyy, dt_Atyz,
        dt_Atzz, dt_Gamtx, dt_Gamty, dt_Gamtz, dt_gbx, dt_gby, dt_gbz,
        shiftx, shifty, shiftz,
        adx_gtxx,ady_gtxx,adz_gtxx,adx_gtxy,
        ady_gtxy,adz_gtxy,adx_gtxz,ady_gtxz,
        adz_gtxz,adx_gtyy,ady_gtyy,adz_gtyy,
        adx_gtyz,ady_gtyz,adz_gtyz,adx_gtzz,
        ady_gtzz,adz_gtzz,adx_Atxx,ady_Atxx,
        adz_Atxx,adx_Atxy,ady_Atxy,adz_Atxy,
        adx_Atxz,ady_Atxz,adz_Atxz,adx_Atyy,
        ady_Atyy,adz_Atyy,adx_Atyz,ady_Atyz,
        adz_Atyz,adx_Atzz,ady_Atzz,adz_Atzz,
        adx_alpha,ady_alpha,adz_alpha,adx_shiftx,
        ady_shiftx,adz_shiftx,adx_shifty,ady_shifty,
        adz_shifty,adx_shiftz,ady_shiftz,adz_shiftz,
        adx_chi,ady_chi,adz_chi,adx_Gamtx,
        ady_Gamtx,adz_Gamtx,adx_Gamty,ady_Gamty,
        adz_Gamty,adx_Gamtz,ady_Gamtz,adz_Gamtz,
        adx_trK,ady_trK,adz_trK,adx_gbx,
        ady_gbx,adz_gbx,adx_gby,ady_gby,
        adz_gby,adx_gbz,ady_gbz,adz_gbz,
        lambda, lambda_f,
        x1d, y1d, z1d, &nx, &ny, &nz);

#include "NewMR_dealloc_adv.h"

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void transpose_xy(double * const restrict fout,
                  const double * const restrict fin,
                  const int nx, const int ny, const int nz)
{ $

  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < nx; j++) {
      for (int i = 0; i < ny; i++) {
        fout[i + j*ny + k*nx*ny] = fin[j + i*nx + k*nx*ny];
      }
    }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void transpose_xz(double * const restrict fout,
                  const double * const restrict fin,
                  const int nx, const int ny, const int nz)
{ $

  for (int k = 0; k < nx; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nz; i++) {
        fout[i + j*nz + k*nz*ny] = fin[k + j*nx + i*nx*ny];
      }
    }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_x(double * const restrict Dxu,
               const double * const restrict u,
               const double dx,
               const int nx, const int ny, const int nz)
{ $

  const double idx = 1.0/dx;
  const double idx_by_2 = 0.5 * idx;
  const double idx_by_12 = idx / 12.0;

  for (int k = 0; k < nz; k++) {
  for (int j = 0; j < ny; j++) {

    Dxu[IDX(0,j,k)] = ( -  3.0 * u[IDX(0,j,k)]
                        +  4.0 * u[IDX(1,j,k)]
                        -        u[IDX(2,j,k)]
                      ) * idx_by_2;

    Dxu[IDX(1,j,k)] = ( - u[IDX(0,j,k)]
                        + u[IDX(2,j,k)]
                      ) * idx_by_2;

    for (int i = 2; i < nx-2; i++) {
      Dxu[IDX(i,j,k)] = (          u[IDX(i-2,j,k)]
                           - 8.0 * u[IDX(i-1,j,k)]
                           + 8.0 * u[IDX(i+1,j,k)]
                           -       u[IDX(i+2,j,k)]
                        ) * idx_by_12;
    }

    Dxu[IDX(nx-2,j,k)] = ( - u[IDX(nx-3,j,k)]
                           + u[IDX(nx-1,j,k)]
                         ) * idx_by_2;

    Dxu[IDX(nx-1,j,k)] = (        u[IDX(nx-3,j,k)]
                          - 4.0 * u[IDX(nx-2,j,k)]
                          + 3.0 * u[IDX(nx-1,j,k)]
                      ) * idx_by_2;
  }
  }

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_y(double * const restrict Dyu,
               const double * const restrict u,
               const double dy,
               const int nx, const int ny, const int nz)
{ $
  const double idy = 1.0/dy;
  const double idy_by_2 = 0.50 * idy;
  const double idy_by_12 = idy / 12.0;

  for (int k = 0; k < nz; k++) {
  for (int i = 0; i < nx; i++) {
    Dyu[IDX(i, 0,k)] = ( - 3.0 * u[IDX(i,0,k)]
                        +  4.0 * u[IDX(i,1,k)]
                        -        u[IDX(i,2,k)]
                      ) * idy_by_2;

    Dyu[IDX(i,1,k)] = ( - u[IDX(i,0,k)]
                        + u[IDX(i,2,k)]
                      ) * idy_by_2;

    for (int j = 2; j < ny-2; j++) {
      Dyu[IDX(i,j,k)] = (          u[IDX(i,j-2,k)]
                           - 8.0 * u[IDX(i,j-1,k)]
                           + 8.0 * u[IDX(i,j+1,k)]
                           -       u[IDX(i,j+2,k)]
                        ) * idy_by_12;
    }

    Dyu[IDX(i,ny-2,k)] = ( - u[IDX(i,ny-3,k)]
                           + u[IDX(i,ny-1,k)]
                         ) * idy_by_2;

    Dyu[IDX(i,ny-1,k)] = (        u[IDX(i,ny-3,k)]
                          - 4.0 * u[IDX(i,ny-2,k)]
                          + 3.0 * u[IDX(i,ny-1,k)]
                      ) * idy_by_2;
  }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_z(double * const restrict Dzu,
               const double * const restrict u,
               const double dz,
               const int nx, const int ny, const int nz)
{ $
  const double idz = 1.0/dz;
  const double idz_by_2 = 0.50 * idz;
  const double idz_by_12 = idz / 12.0;

  for (int j = 0; j < ny; j++) {
  for (int i = 0; i < nx; i++) {
    Dzu[IDX(i, j, 0)] = ( - 3.0 *  u[IDX(i,j,0)]
                          +  4.0 * u[IDX(i,j,1)]
                          -        u[IDX(i,j,2)]
                        ) * idz_by_2;

    Dzu[IDX(i,j,1)] = ( - u[IDX(i,j,0)]
                        + u[IDX(i,j,2)]
                      ) * idz_by_2;

    for (int k = 2; k < nz-2; k++) {
      Dzu[IDX(i,j,k)] = (          u[IDX(i,j,k-2)]
                           - 8.0 * u[IDX(i,j,k-1)]
                           + 8.0 * u[IDX(i,j,k+1)]
                           -       u[IDX(i,j,k+2)]
                        ) * idz_by_12;
    }

    Dzu[IDX(i,j,nz-2)] = ( - u[IDX(i,j,nz-3)]
                           + u[IDX(i,j,nz-1)]
                         ) * idz_by_2;

    Dzu[IDX(i,j,nz-1)] = (        u[IDX(i,j,nz-3)]
                          - 4.0 * u[IDX(i,j,nz-2)]
                          + 3.0 * u[IDX(i,j,nz-1)]
                      ) * idz_by_2;
  }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_xx(double * const restrict DxDxu,
                const double * const restrict u,
                const double dx, const int nx, const int ny, const int nz)
{ $

  const double idx_sqrd = 1.0/(dx*dx);
  const double idx_sqrd_by_12 = idx_sqrd / 12.0;

  for (int k = 0; k < nz; k++) {
  for (int j = 0; j < ny; j++) {

    DxDxu[IDX(0,j,k)] = (   2.0 * u[IDX(0,j,k)]
                           - 5.0 * u[IDX(1,j,k)]
                           + 4.0 * u[IDX(2,j,k)]
                           -        u[IDX(3,j,k)]
                        ) * idx_sqrd;

    DxDxu[IDX(1,j,k)] = (          u[IDX(0,j,k)]
                           - 2.0 * u[IDX(1  ,j,k)]
                           +       u[IDX(2,j,k)]
                        ) * idx_sqrd;
    for (int i = 2; i < nx-2; i++) {
      DxDxu[IDX(i,j,k)] = ( -         u[IDX(i-2,j,k)]
                             + 16.0 * u[IDX(i-1,j,k)]
                             - 30.0 * u[IDX(i  ,j,k)]
                             + 16.0 * u[IDX(i+1,j,k)]
                             -         u[IDX(i+2,j,k)]
                          ) * idx_sqrd_by_12;
    }

    DxDxu[IDX(nx-2,j,k)] = (       u[IDX(nx-3,j,k)]
                           - 2.0 * u[IDX(nx-2,j,k)]
                           +       u[IDX(nx-1,j,k)]
                          ) * idx_sqrd;

    DxDxu[IDX(nx-1,j,k)] = ( -    u[IDX(nx-4,j,k)]
                          + 4.0 * u[IDX(nx-3,j,k)]
                          - 5.0 * u[IDX(nx-2,j,k)]
                          + 2.0 * u[IDX(nx-1,j,k)]
                        ) * idx_sqrd;
  }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_yy(double * const restrict Du,
                const double * const restrict u,
                const double dy, const int nx, const int ny, const int nz)
{ $

  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_12 = idy_sqrd / 12.0;

  for (int k = 0; k < nz; k++) {
  for (int i = 0; i < nx; i++) {

    Du[IDX(i,0,k)] = (       2.0 * u[IDX(i,0,k)]
                           - 5.0 * u[IDX(i,1,k)]
                           + 4.0 * u[IDX(i,2,k)]
                           -       u[IDX(i,3,k)]
                        ) * idy_sqrd;

    Du[IDX(i,1,k)] = (             u[IDX(i,0,k)]
                           - 2.0 * u[IDX(i,1,k)]
                           +       u[IDX(i,2,k)]
                        ) * idy_sqrd;
    for (int j = 2; j < ny-2; j++) {
      Du[IDX(i,j,k)] = ( -            u[IDX(i,j-2,k)]
                             + 16.0 * u[IDX(i,j-1,k)]
                             - 30.0 * u[IDX(i,j  ,k)]
                             + 16.0 * u[IDX(i,j+1,k)]
                             -        u[IDX(i,j+2,k)]
                          ) * idy_sqrd_by_12;
    }

    Du[IDX(i,ny-2,k)] = (          u[IDX(i,ny-3,k)]
                           - 2.0 * u[IDX(i,ny-2,k)]
                           +       u[IDX(i,ny-1,k)]
                          ) * idy_sqrd;

    Du[IDX(i,ny-1,k)] = ( -       u[IDX(i,ny-4,k)]
                          + 4.0 * u[IDX(i,ny-3,k)]
                          - 5.0 * u[IDX(i,ny-2,k)]
                          + 2.0 * u[IDX(i,ny-1,k)]
                        ) * idy_sqrd;
  }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_zz(double * const restrict Du,
                const double * const restrict u,
                const double dz, const int nx, const int ny, const int nz)
{ $

  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_12 = idz_sqrd / 12.0;

  for (int j = 0; j < ny; j++) {
  for (int i = 0; i < nx; i++) {

    Du[IDX(i,j,0)] = (       2.0 * u[IDX(i,j,0)]
                           - 5.0 * u[IDX(i,j,1)]
                           + 4.0 * u[IDX(i,j,2)]
                           -       u[IDX(i,j,3)]
                        ) * idz_sqrd;

    Du[IDX(i,j,1)] = (             u[IDX(i,j,0)]
                           - 2.0 * u[IDX(i,j,1)]
                           +       u[IDX(i,j,2)]
                        ) * idz_sqrd;
    for (int k = 2; k < nz-2; k++) {
      Du[IDX(i,j,k)] = ( -            u[IDX(i,j,k-2)]
                             + 16.0 * u[IDX(i,j,k-1)]
                             - 30.0 * u[IDX(i,j,k  )]
                             + 16.0 * u[IDX(i,j,k+1)]
                             -        u[IDX(i,j,k+2)]
                          ) * idz_sqrd_by_12;
    }

    Du[IDX(i,j,nz-2)] = (          u[IDX(i,j,nz-3)]
                           - 2.0 * u[IDX(i,j,nz-2)]
                           +       u[IDX(i,j,nz-1)]
                          ) * idz_sqrd;

    Du[IDX(i,j,nz-1)] = ( -       u[IDX(i,j,nz-4)]
                          + 4.0 * u[IDX(i,j,nz-3)]
                          - 5.0 * u[IDX(i,j,nz-2)]
                          + 2.0 * u[IDX(i,j,nz-1)]
                        ) * idz_sqrd;
  }
  }

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_x(double * const restrict Dxu,
                  const double * const restrict u,
                  const double dx,
                  const int nx, const int ny, const int nz,
                  const double * const betax)
{ $

  const double idx = 1.0/dx;
  const double idx_by_2 = 0.50 * idx;
  const double idx_by_12 = idx / 12.0;


  for (int k = 0; k < nz; k++) {
  for (int j = 0; j < ny; j++) {

    Dxu[IDX(0,j,k)] = ( -  3.0 * u[IDX(0  ,j,k)]
                        +  4.0 * u[IDX(1,j,k)]
                        -        u[IDX(2,j,k)]
                      ) * idx_by_2;

    if (betax[IDX(1,j,k)] > 0.0) {
      Dxu[IDX(1,j,k)] = ( -  3.0 * u[IDX(1,j,k)]
                          +  4.0 * u[IDX(2,j,k)]
                          -        u[IDX(3,j,k)]
                        ) * idx_by_2;
    }
    else {
      Dxu[IDX(1,j,k)] = ( -         u[IDX(0,j,k)]
                           +        u[IDX(2,j,k)]
                        ) * idx_by_2;
    }

    if (betax[IDX(2,j,k)] > 0.0 ) {
      Dxu[IDX(2,j,k)] = ( -  3.0 * u[IDX(1,j,k)]
                         - 10.0 * u[IDX(2  ,j,k)]
                         + 18.0 * u[IDX(3,j,k)]
                         -  6.0 * u[IDX(4,j,k)]
                         +         u[IDX(5,j,k)]
                       ) * idx_by_12;
    }
    else {
      Dxu[IDX(2,j,k)] = (           u[IDX(0,j,k)]
                           -  4.0 * u[IDX(1,j,k)]
                           +  3.0 * u[IDX(2,j,k)]
                        ) * idx_by_2;
    }

    for (int i = 3; i < nx-3; i++) {
      if (betax[IDX(i,j,k)] > 0.0 ) {
        Dxu[IDX(i,j,k)] = ( -  3.0 * u[IDX(i-1,j,k)]
                            - 10.0 * u[IDX(i  ,j,k)]
                            + 18.0 * u[IDX(i+1,j,k)]
                            -  6.0 * u[IDX(i+2,j,k)]
                            +        u[IDX(i+3,j,k)]
                          ) * idx_by_12;
      }
      else {
        Dxu[IDX(i,j,k)] = ( -        u[IDX(i-3,j,k)]
                            +  6.0 * u[IDX(i-2,j,k)]
                            - 18.0 * u[IDX(i-1,j,k)]
                            + 10.0 * u[IDX(i  ,j,k)]
                            +  3.0 * u[IDX(i+1,j,k)]
                          ) * idx_by_12;
      }
    }

    if ( betax[IDX(nx-3,j,k)] < 0.0 ) {
      Dxu[IDX(nx-3,j,k)] = (  - 3.0 * u[IDX(nx-3,j,k)]
                              + 4.0 * u[IDX(nx-2,j,k)]
                              -       u[IDX(nx-1,j,k)]
                           ) * idx_by_2;
    }
    else {
      Dxu[IDX(nx-3,j,k)] = ( -   u[IDX(nx-6,j,k)]
                        +  6.0 * u[IDX(nx-5,j,k)]
                        - 18.0 * u[IDX(nx-4,j,k)]
                        + 10.0 * u[IDX(nx-3  ,j,k)]
                        +  3.0 * u[IDX(nx-2,j,k)]
                      ) * idx_by_12;
    }

    if (betax[IDX(nx-2,j,k)] > 0.0 ) {
      Dxu[IDX(nx-2,j,k)] = (  -  u[IDX(nx-3,j,k)]
                              +  u[IDX(nx-1,j,k)]
                           ) * idx_by_2;
    }
    else {
      Dxu[IDX(nx-2,j,k)] = (     u[IDX(nx-4,j,k)]
                         - 4.0 * u[IDX(nx-3,j,k)]
                         + 3.0 * u[IDX(nx-2,j,k)]
                           ) * idx_by_2;
    }

    Dxu[IDX(nx-1,j,k)] = (          u[IDX(nx-3,j,k)]
                            - 4.0 * u[IDX(nx-2,j,k)]
                            + 3.0 * u[IDX(nx-1,j,k)]
                         ) * idx_by_2;
  }
  }

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_y(double * const restrict Dyu,
                  const double * const restrict u,
                  const double dy,
                  const int nx, const int ny, const int nz,
                  const double * const betay)
{ $

  const double idy = 1.0/dy;
  const double idy_by_2 = 0.50 * idy;
  const double idy_by_12 = idy / 12.0;


  for (int k = 0; k < nz; k++) {
  for (int i = 0; i < nx; i++) {

    Dyu[IDX(i,0,k)] = ( -  3.0 * u[IDX(i,0,k)]
                        +  4.0 * u[IDX(i,1,k)]
                        -        u[IDX(i,2,k)]
                      ) * idy_by_2;

    if (betay[IDX(i,1,k)] > 0.0) {
      Dyu[IDX(i,1,k)] = ( -  3.0 * u[IDX(i,1,k)]
                          +  4.0 * u[IDX(i,2,k)]
                          -        u[IDX(i,3,k)]
                        ) * idy_by_2;
    }
    else {
      Dyu[IDX(i,1,k)] = ( -         u[IDX(i,0,k)]
                           +        u[IDX(i,2,k)]
                        ) * idy_by_2;
    }

    if (betay[IDX(i,2,k)] > 0.0 ) {
      Dyu[IDX(i,2,k)] = ( -  3.0 * u[IDX(i,1,k)]
                          - 10.0 * u[IDX(i,2,k)]
                          + 18.0 * u[IDX(i,3,k)]
                          -  6.0 * u[IDX(i,4,k)]
                         +         u[IDX(i,5,k)]
                       ) * idy_by_12;
    }
    else {
      Dyu[IDX(i,2,k)] = (           u[IDX(i,0,k)]
                           -  4.0 * u[IDX(i,1,k)]
                           +  3.0 * u[IDX(i,2,k)]
                        ) * idy_by_2;
    }

    for (int j = 3; j < ny-3; j++) {
      if (betay[IDX(i,j,k)] > 0.0 ) {
        Dyu[IDX(i,j,k)] = ( -  3.0 * u[IDX(i,j-1,k)]
                            - 10.0 * u[IDX(i,j  ,k)]
                            + 18.0 * u[IDX(i,j+1,k)]
                            -  6.0 * u[IDX(i,j+2,k)]
                            +        u[IDX(i,j+3,k)]
                          ) * idy_by_12;
      }
      else {
        Dyu[IDX(i,j,k)] = ( -        u[IDX(i,j-3,k)]
                            +  6.0 * u[IDX(i,j-2,k)]
                            - 18.0 * u[IDX(i,j-1,k)]
                            + 10.0 * u[IDX(i,j  ,k)]
                            +  3.0 * u[IDX(i,j+1,k)]
                          ) * idy_by_12;
      }
    }

    if ( betay[IDX(i,ny-3,k)] < 0.0 ) {
      Dyu[IDX(i,ny-3,k)] = (  - 3.0 * u[IDX(i,ny-3,k)]
                              + 4.0 * u[IDX(i,ny-2,k)]
                              -       u[IDX(i,ny-1,k)]
                           ) * idy_by_2;
    }
    else {
      Dyu[IDX(i,ny-3,k)] = ( -   u[IDX(i,ny-6,k)]
                        +  6.0 * u[IDX(i,ny-5,k)]
                        - 18.0 * u[IDX(i,ny-4,k)]
                        + 10.0 * u[IDX(i,ny-3,k)]
                        +  3.0 * u[IDX(i,ny-2,k)]
                      ) * idy_by_12;
    }

    if (betay[IDX(i,ny-2,k)] > 0.0 ) {
      Dyu[IDX(i,ny-2,k)] = (  -  u[IDX(i,ny-3,k)]
                              +  u[IDX(i,ny-1,k)]
                           ) * idy_by_2;
    }
    else {
      Dyu[IDX(i,ny-2,k)] = (     u[IDX(i,ny-4,k)]
                         - 4.0 * u[IDX(i,ny-3,k)]
                         + 3.0 * u[IDX(i,ny-2,k)]
                           ) * idy_by_2;
    }

    Dyu[IDX(i,ny-1,k)] = (          u[IDX(i,ny-3,k)]
                            - 4.0 * u[IDX(i,ny-2,k)]
                            + 3.0 * u[IDX(i,ny-1,k)]
                         ) * idy_by_2;
  }
  }

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_z(double * const restrict Dzu,
                  const double * const restrict u,
                  const double dz,
                  const int nx, const int ny, const int nz,
                  const double * const betaz)
{ $

  const double idz = 1.0/dz;
  const double idz_by_2 = 0.50 * idz;
  const double idz_by_12 = idz / 12.0;


  for (int j = 0; j < ny; j++) {
  for (int i = 0; i < nx; i++) {

    Dzu[IDX(i,j,0)] = ( -  3.0 * u[IDX(i,j,0)]
                        +  4.0 * u[IDX(i,j,1)]
                        -        u[IDX(i,j,2)]
                      ) * idz_by_2;

    if (betaz[IDX(i,j,1)] > 0.0) {
      Dzu[IDX(i,j,1)] = ( -  3.0 * u[IDX(i,j,1)]
                          +  4.0 * u[IDX(i,j,2)]
                          -        u[IDX(i,j,3)]
                        ) * idz_by_2;
    }
    else {
      Dzu[IDX(i,j,1)] = ( -         u[IDX(i,j,0)]
                           +        u[IDX(i,j,2)]
                        ) * idz_by_2;
    }

    if (betaz[IDX(i,j,2)] > 0.0 ) {
      Dzu[IDX(i,j,2)] = ( -  3.0 * u[IDX(i,j,1)]
                          - 10.0 * u[IDX(i,j,2)]
                          + 18.0 * u[IDX(i,j,3)]
                          -  6.0 * u[IDX(i,j,4)]
                         +         u[IDX(i,j,5)]
                       ) * idz_by_12;
    }
    else {
      Dzu[IDX(i,j,2)] = (           u[IDX(i,j,0)]
                           -  4.0 * u[IDX(i,j,1)]
                           +  3.0 * u[IDX(i,j,2)]
                        ) * idz_by_2;
    }

    for (int k = 3; k < nz-3; k++) {
      if (betaz[IDX(i,j,k)] > 0.0 ) {
        Dzu[IDX(i,j,k)] = ( -  3.0 * u[IDX(i,j,k-1)]
                            - 10.0 * u[IDX(i,j,k  )]
                            + 18.0 * u[IDX(i,j,k+1)]
                            -  6.0 * u[IDX(i,j,k+2)]
                            +        u[IDX(i,j,k+3)]
                          ) * idz_by_12;
      }
      else {
        Dzu[IDX(i,j,k)] = ( -        u[IDX(i,j,k-3)]
                            +  6.0 * u[IDX(i,j,k-2)]
                            - 18.0 * u[IDX(i,j,k-1)]
                            + 10.0 * u[IDX(i,j,k  )]
                            +  3.0 * u[IDX(i,j,k+1)]
                          ) * idz_by_12;
      }
    }

    if ( betaz[IDX(i,j,nz-3)] < 0.0 ) {
      Dzu[IDX(i,j,nz-3)] = (  - 3.0 * u[IDX(i,j,nz-3)]
                              + 4.0 * u[IDX(i,j,nz-2)]
                              -       u[IDX(i,j,nz-1)]
                           ) * idz_by_2;
    }
    else {
      Dzu[IDX(i,j,nz-3)] = ( -   u[IDX(i,j,nz-6)]
                        +  6.0 * u[IDX(i,j,nz-5)]
                        - 18.0 * u[IDX(i,j,nz-4)]
                        + 10.0 * u[IDX(i,j,nz-3)]
                        +  3.0 * u[IDX(i,j,nz-2)]
                      ) * idz_by_12;
    }

    if (betaz[IDX(i,j,nz-2)] > 0.0 ) {
      Dzu[IDX(i,j,nz-2)] = (  -  u[IDX(i,j,nz-3)]
                              +  u[IDX(i,j,nz-1)]
                           ) * idz_by_2;
    }
    else {
      Dzu[IDX(i,j,nz-2)] = (     u[IDX(i,j,nz-4)]
                         - 4.0 * u[IDX(i,j,nz-3)]
                         + 3.0 * u[IDX(i,j,nz-2)]
                           ) * idz_by_2;
    }

    Dzu[IDX(i,j,nz-1)] = (          u[IDX(i,j,nz-3)]
                            - 4.0 * u[IDX(i,j,nz-2)]
                            + 3.0 * u[IDX(i,j,nz-1)]
                         ) * idz_by_2;
  }
  }

}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv42_x(double * const restrict Du,
                const double * const restrict u,
                const double dx, 
                const int nx, const int ny, const int nz,
                int mode)
{ $

  double pre_factor_6_dx = -1.0 / 64.0 / dx;

  double smr3=59.0/48.0*64*dx;
  double smr2=43.0/48.0*64*dx;
  double smr1=49.0/48.0*64*dx;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;


  const int lef = 3;
  const int rig = nx-3;

  if (mode == 0) {
    int n = nx * ny * nz;
    for (int i = 0; i < n; i++) {
      Du[i] = 0.0;
    }
    return;
  }

  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {

       Du[IDX(0,j,k)] =  (      u[IDX(3,j,k)] 
                          - 3.0*u[IDX(2,j,k)] 
                          + 3.0*u[IDX(1,j,k)]
                          -     u[IDX(0,j,k)]
                         )/smr3;
       Du[IDX(1,j,k)] =  (
                                u[IDX(4,j,k)]
                         -  6.0*u[IDX(3,j,k)]
                         + 12.0*u[IDX(2,j,k)]
                         - 10.0*u[IDX(1,j,k)]
                         +  3.0*u[IDX(0,j,k)]
                         )/smr2;
       Du[IDX(3,j,k)] =  (
                                u[IDX(5,j,k)]
                         -  6.0*u[IDX(4,j,k)]
                         + 15.0*u[IDX(3,j,k)]
                         - 19.0*u[IDX(2,j,k)]
                         + 12.0*u[IDX(1,j,k)]
                         -  3.0*u[IDX(0,j,k)]
                         )/smr1;

       for (int i = lef; i < rig; i++) {
          Du[IDX(i,j,k)] = pre_factor_6_dx * 
                         ( 
                         -      u[IDX(i-3,j,k)]
                         +  6.0*u[IDX(i-2,j,k)]
                         - 15.0*u[IDX(i-1,j,k)]
                         + 20.0*u[IDX(i  ,j,k)]
                         - 15.0*u[IDX(i+1,j,k)]
                         +  6.0*u[IDX(i+2,j,k)]
                         -      u[IDX(i+3,j,k)]
                         );
       }

       Du[IDX(nx-3,j,k)] = (
                                u[IDX(rig-3,j,k)]
                          - 6.0*u[IDX(rig-2,j,k)]
                          + 15.0*u[IDX(rig-1,j,k)]
                          - 19.0*u[IDX(rig,j,k)]
                          + 12.0*u[IDX(rig+1,j,k)]
                          -  3.0*u[IDX(rig+2,j,k)]
                           )/spr1;

       Du[IDX(nx-2,j,k)] = (
                                 u[IDX(rig-2,j,k)]
                          -  6.0*u[IDX(rig-1,j,k)]
                          + 12.0*u[IDX(rig,j,k)]
                          - 10.0*u[IDX(rig+1,j,k)]
                          +  3.0*u[IDX(rig+2,j,k)]
                           )/spr2;

       Du[IDX(nx-1,j,k)] = (
                                 u[IDX(rig-1,j,k)]
                          -  3.0*u[IDX(rig,j,k)]
                          +  3.0*u[IDX(rig+1,j,k)]
                          -      u[IDX(rig+2,j,k)]
                           )/spr3;
    }
  }

}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv42_y(double * const restrict Du,
                const double * const restrict u,
                const double dy, 
                const int nx, const int ny, const int nz,
                int mode)
{ $

  double pre_factor_6_dy = -1.0 / 64.0 / dy;

  double smr3=59.0/48.0*64*dy;
  double smr2=43.0/48.0*64*dy;
  double smr1=49.0/48.0*64*dy;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;


  const int lef = 3;
  const int rig = ny-3;

  if (mode == 0) {
    int n = nx * ny * nz;
    for (int i = 0; i < n; i++) {
      Du[i] = 0.0;
    }
    return;
  }

  for (int k = 0; k < nz; k++) {
    for (int i = 0; i < nx; i++) {

       Du[IDX(i,0,k)] =  (      u[IDX(i,3,k)] 
                          - 3.0*u[IDX(i,2,k)] 
                          + 3.0*u[IDX(i,1,k)]
                          -     u[IDX(i,0,k)]
                         )/smr3;
       Du[IDX(i,1,k)] =  (
                                u[IDX(i,4,k)]
                         -  6.0*u[IDX(i,3,k)]
                         + 12.0*u[IDX(i,2,k)]
                         - 10.0*u[IDX(i,1,k)]
                         +  3.0*u[IDX(i,0,k)]
                         )/smr2;
       Du[IDX(i,3,k)] =  (
                                u[IDX(i,5,k)]
                         -  6.0*u[IDX(i,4,k)]
                         + 15.0*u[IDX(i,3,k)]
                         - 19.0*u[IDX(i,2,k)]
                         + 12.0*u[IDX(i,1,k)]
                         -  3.0*u[IDX(i,0,k)]
                         )/smr1;

       for (int j = lef; j < rig; j++) {
          Du[IDX(i,j,k)] = pre_factor_6_dy * 
                         ( 
                         -      u[IDX(i,j-3,k)]
                         +  6.0*u[IDX(i,j-2,k)]
                         - 15.0*u[IDX(i,j-1,k)]
                         + 20.0*u[IDX(i,j  ,k)]
                         - 15.0*u[IDX(i,j+1,k)]
                         +  6.0*u[IDX(i,j+2,k)]
                         -      u[IDX(i,j+3,k)]
                         );
       }

       Du[IDX(i,ny-3,k)] = (
                                 u[IDX(i,rig-3,k)]
                          -  6.0*u[IDX(i,rig-2,k)]
                          + 15.0*u[IDX(i,rig-1,k)]
                          - 19.0*u[IDX(i,rig  ,k)]
                          + 12.0*u[IDX(i,rig+1,k)]
                          -  3.0*u[IDX(i,rig+2,k)]
                           )/spr1;

       Du[IDX(i,ny-2,k)] = (
                                 u[IDX(i,rig-2,k)]
                          -  6.0*u[IDX(i,rig-1,k)]
                          + 12.0*u[IDX(i,rig  ,k)]
                          - 10.0*u[IDX(i,rig+1,k)]
                          +  3.0*u[IDX(i,rig+2,k)]
                           )/spr2;

       Du[IDX(i,ny-1,k)] = (
                                 u[IDX(i,rig-1,k)]
                          -  3.0*u[IDX(i,rig,k)]
                          +  3.0*u[IDX(i,rig+1,k)]
                          -      u[IDX(i,rig+2,k)]
                           )/spr3;
    }
  }


}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv42_z(double * const restrict Du,
                const double * const restrict u,
                const double dz, 
                const int nx, const int ny, const int nz,
                int mode)
{ $

  double pre_factor_6_dz = -1.0 / 64.0 / dz;

  double smr3=59.0/48.0*64*dz;
  double smr2=43.0/48.0*64*dz;
  double smr1=49.0/48.0*64*dz;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;


  const int lef = 3;
  const int rig = nz-3;

  if (mode == 0) {
    int n = nx * ny * nz;
    for (int i = 0; i < n; i++) {
      Du[i] = 0.0;
    }
    return;
  }

  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {

       Du[IDX(i,j,0)] =  (      u[IDX(i,j,3)] 
                          - 3.0*u[IDX(i,j,2)] 
                          + 3.0*u[IDX(i,j,1)]
                          -     u[IDX(i,j,0)]
                         )/smr3;
       Du[IDX(i,j,1)] =  (
                                u[IDX(i,j,4)]
                         -  6.0*u[IDX(i,j,3)]
                         + 12.0*u[IDX(i,j,2)]
                         - 10.0*u[IDX(i,j,1)]
                         +  3.0*u[IDX(i,j,0)]
                         )/smr2;
       Du[IDX(i,j,3)] =  (
                                u[IDX(i,j,5)]
                         -  6.0*u[IDX(i,j,4)]
                         + 15.0*u[IDX(i,j,3)]
                         - 19.0*u[IDX(i,j,2)]
                         + 12.0*u[IDX(i,j,1)]
                         -  3.0*u[IDX(i,j,0)]
                         )/smr1;

       for (int k = lef; k < rig; k++) {
          Du[IDX(i,j,j)] = pre_factor_6_dz * 
                         ( 
                         -      u[IDX(i,j,k-3)]
                         +  6.0*u[IDX(i,j,k-2)]
                         - 15.0*u[IDX(i,j,k-1)]
                         + 20.0*u[IDX(i,j,k  )]
                         - 15.0*u[IDX(i,j,k+1)]
                         +  6.0*u[IDX(i,j,k+2)]
                         -      u[IDX(i,j,k+3)]
                         );
       }

       Du[IDX(i,j,nz-3)] = (
                                 u[IDX(i,j,rig-3)]
                          -  6.0*u[IDX(i,j,rig-2)]
                          + 15.0*u[IDX(i,j,rig-1)]
                          - 19.0*u[IDX(i,j,rig  )]
                          + 12.0*u[IDX(i,j,rig+1)]
                          -  3.0*u[IDX(i,j,rig+2)]
                           )/spr1;

       Du[IDX(i,j,nz-2)] = (
                                 u[IDX(i,j,rig-2)]
                          -  6.0*u[IDX(i,j,rig-1)]
                          + 12.0*u[IDX(i,j,rig  )]
                          - 10.0*u[IDX(i,j,rig+1)]
                          +  3.0*u[IDX(i,j,rig+2)]
                           )/spr2;

       Du[IDX(i,j,nz-1)] = (
                                 u[IDX(i,j,rig-1)]
                          -  3.0*u[IDX(i,j,rig)]
                          +  3.0*u[IDX(i,j,rig+1)]
                          -      u[IDX(i,j,rig+2)]
                           )/spr3;

    }
  }

}





