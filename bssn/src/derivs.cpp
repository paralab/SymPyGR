#include <cmath>
#include <iostream>
#include "derivs.h"

unsigned int DERIV_PW=3;

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_x(double * const  Dxu, const double * const  u,
               const double dx, const unsigned int *sz, unsigned bflag)
{

  const double idx = 1.0/dx;
  const double idx_by_2 = 0.5 * idx;
  const double idx_by_12 = idx / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 1;
  const int kb = 1;
  const int ie = sz[0]-3;
  const int je = sz[1]-1;
  const int ke = sz[2]-1;
  const int n=1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        Dxu[pp] = (u[pp-2] -8.0*u[pp-1] + 8.0*u[pp+1] - u[pp+2] ) * idx_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Dxu[IDX(3,j,k)] = ( -  3.0 * u[IDX(3,j,k)]
                            +  4.0 * u[IDX(4,j,k)]
                            -        u[IDX(5,j,k)]
                          ) * idx_by_2;
        Dxu[IDX(4,j,k)] = ( - u[IDX(3,j,k)]
                            + u[IDX(5,j,k)]
                          ) * idx_by_2;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Dxu[IDX(ie-2,j,k)] = ( - u[IDX(ie-3,j,k)]
                               + u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

        Dxu[IDX(ie-1,j,k)] = (        u[IDX(ie-3,j,k)]
                              - 4.0 * u[IDX(ie-2,j,k)]
                              + 3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
#pragma message("DEBUG_DERIVS_COMP: ON")
  for (int k = 3; k < sz[2]-3; k++) {
    for (int j = 3; j < sz[1]-3; j++) {
      for (int i = 3; i < sz[0]-3; i++) {
        const int pp = IDX(i,j,k);
         if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif


}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_y(double * const  Dyu, const double * const  u,
               const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy = 1.0/dy;
  const double idy_by_2 = 0.50 * idy;
  const double idy_by_12 = idy / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 1;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-1;

  const int n=nx;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        Dyu[pp] = (u[pp-2*nx] - 8.0*u[pp-nx] + 8.0*u[pp+nx] - u[pp+2*nx])*idy_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        Dyu[IDX(i, 3,k)] = ( - 3.0 * u[IDX(i,3,k)]
                            +  4.0 * u[IDX(i,4,k)]
                            -        u[IDX(i,5,k)]
                          ) * idy_by_2;

        Dyu[IDX(i,4,k)] = ( - u[IDX(i,3,k)]
                            + u[IDX(i,5,k)]
                          ) * idy_by_2;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        Dyu[IDX(i,je-2,k)] = ( - u[IDX(i,je-3,k)]
                               + u[IDX(i,je-1,k)]
                             ) * idy_by_2;

        Dyu[IDX(i,je-1,k)] = (        u[IDX(i,je-3,k)]
                              - 4.0 * u[IDX(i,je-2,k)]
                              + 3.0 * u[IDX(i,je-1,k)]
                          ) * idy_by_2;
      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = 3; k < sz[2]-3; k++) {
    for (int j = 3; j < sz[1]-3; j++) {
      for (int i = 3; i < sz[0]-3; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_z(double * const  Dzu, const double * const  u,
               const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz = 1.0/dz;
  const double idz_by_2 = 0.50 * idz;
  const double idz_by_12 = idz / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-3;

  const int n = nx*ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        Dzu[pp] = (u[pp-2*n] - 8.0*u[pp-n] + 8.0*u[pp+n] - u[pp+2*n]) * idz_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i, j, 3)] = ( - 3.0 *  u[IDX(i,j,3)]
                              +  4.0 * u[IDX(i,j,4)]
                              -        u[IDX(i,j,5)]
                            ) * idz_by_2;

        Dzu[IDX(i,j,4)] = ( - u[IDX(i,j,3)]
                            + u[IDX(i,j,5)]
                          ) * idz_by_2;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i,j,ke-2)] = ( - u[IDX(i,j,ke-3)]
                               + u[IDX(i,j,ke-1)]
                             ) * idz_by_2;

        Dzu[IDX(i,j,ke-1)] = (        u[IDX(i,j,ke-3)]
                              - 4.0 * u[IDX(i,j,ke-2)]
                              + 3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_2;

      }
    }

  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_xx(double * const  DxDxu, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  const double idx_sqrd = 1.0/(dx*dx);
  const double idx_sqrd_by_12 = idx_sqrd / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        DxDxu[pp] = (   -        u[pp-2]
                        + 16.0 * u[pp-1]
                        - 30.0 * u[pp]
                        + 16.0 * u[pp+1]
                        -        u[pp+2]
                    ) * idx_sqrd_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        DxDxu[IDX(3,j,k)] = (   2.0 * u[IDX(3,j,k)]
                              - 5.0 * u[IDX(4,j,k)]
                              + 4.0 * u[IDX(5,j,k)]
                              -       u[IDX(6,j,k)]
                            ) * idx_sqrd;

        DxDxu[IDX(4,j,k)] = (         u[IDX(3,j,k)]
                              - 2.0 * u[IDX(4,j,k)]
                              +       u[IDX(5,j,k)]
                            ) * idx_sqrd;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        DxDxu[IDX(ie-2,j,k)] = (       u[IDX(ie-3,j,k)]
                               - 2.0 * u[IDX(ie-2,j,k)]
                               +       u[IDX(ie-1,j,k)]
                              ) * idx_sqrd;

        DxDxu[IDX(ie-1,j,k)] = ( -    u[IDX(ie-4,j,k)]
                          + 4.0 * u[IDX(ie-3,j,k)]
                          - 5.0 * u[IDX(ie-2,j,k)]
                          + 2.0 * u[IDX(ie-1,j,k)]
                        ) * idx_sqrd;

      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_yy(double * const  DyDyu, const double * const  u,
                const double dy, const unsigned int *sz, unsigned bflag)
{

  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_12 = idy_sqrd / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        DyDyu[pp] = ( -u[pp-2*nx] + 16.0 * u[pp-nx] - 30.0 * u[pp]
                               + 16.0 * u[pp+nx] - u[pp+2*nx]
                 ) * idy_sqrd_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        DyDyu[IDX(i,3,k)] = (   2.0 * u[IDX(i,3,k)]
                           - 5.0 * u[IDX(i,4,k)]
                           + 4.0 * u[IDX(i,5,k)]
                           -       u[IDX(i,6,k)]
                        ) * idy_sqrd;

        DyDyu[IDX(i,4,k)] = (         u[IDX(i,3,k)]
                           - 2.0 * u[IDX(i,4,k)]
                           +       u[IDX(i,5,k)]
                        ) * idy_sqrd;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        DyDyu[IDX(i,je-2,k)] = (      u[IDX(i,je-3,k)]
                           - 2.0 * u[IDX(i,je-2,k)]
                           +       u[IDX(i,je-1,k)]
                          ) * idy_sqrd;

        DyDyu[IDX(i,je-1,k)] = ( -   u[IDX(i,je-4,k)]
                          + 4.0 * u[IDX(i,je-3,k)]
                          - 5.0 * u[IDX(i,je-2,k)]
                          + 2.0 * u[IDX(i,je-1,k)]
                        ) * idy_sqrd;

      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_zz(double * const  DzDzu, const double * const  u,
                const double dz, const unsigned int *sz, unsigned bflag)
{

  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_12 = idz_sqrd / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        DzDzu[pp] = ( - u[pp-2*n] + 16.0 * u[pp-n] - 30.0 * u[pp]
                   + 16.0 * u[pp+n] - u[pp+2*n] ) * idz_sqrd_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        DzDzu[IDX(i,j,3)] = (   2.0 * u[IDX(i,j,3)]
                           - 5.0 * u[IDX(i,j,4)]
                           + 4.0 * u[IDX(i,j,5)]
                           -       u[IDX(i,j,6)]
                        ) * idz_sqrd;

        DzDzu[IDX(i,j,4)] = (         u[IDX(i,j,3)]
                           - 2.0 * u[IDX(i,j,4)]
                           +       u[IDX(i,j,5)]
                        ) * idz_sqrd;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        DzDzu[IDX(i,j,ke-2)] = (      u[IDX(i,j,ke-3)]
                           - 2.0 * u[IDX(i,j,ke-2)]
                           +       u[IDX(i,j,ke-1)]
                          ) * idz_sqrd;

        DzDzu[IDX(i,j,ke-1)] = ( -   u[IDX(i,j,ke-4)]
                          + 4.0 * u[IDX(i,j,ke-3)]
                          - 5.0 * u[IDX(i,j,ke-2)]
                          + 2.0 * u[IDX(i,j,ke-1)]
                        ) * idz_sqrd;

      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_x(double * const  Dxu, const double * const  u,
                  const double dx, const unsigned int *sz,
                  const double * const betax, unsigned bflag)
{

  const double idx = 1.0/dx;
  const double idx_by_2 = 0.50 * idx;
  const double idx_by_12 = idx / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if (betax[pp] > 0.0 ) {
          Dxu[pp] = ( -  3.0 * u[pp-1]
                      - 10.0 * u[pp]
                      + 18.0 * u[pp+1]
                      -  6.0 * u[pp+2]
                      +        u[pp+3]
                    ) * idx_by_12;
        }
        else {
          Dxu[pp] = ( -        u[pp-3]
                      +  6.0 * u[pp-2]
                      - 18.0 * u[pp-1]
                      + 10.0 * u[pp]
                      +  3.0 * u[pp+1]
                    ) * idx_by_12;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Dxu[IDX(3,j,k)] = ( -  3.0 * u[IDX(3,j,k)]
                            +  4.0 * u[IDX(4,j,k)]
                            -        u[IDX(5,j,k)]
                          ) * idx_by_2;

        if (betax[IDX(4,j,k)] > 0.0) {
          Dxu[IDX(4,j,k)] = ( -  3.0 * u[IDX(4,j,k)]
                              +  4.0 * u[IDX(5,j,k)]
                              -        u[IDX(6,j,k)]
                            ) * idx_by_2;
        }
        else {
          Dxu[IDX(4,j,k)] = ( -         u[IDX(3,j,k)]
                               +        u[IDX(5,j,k)]
                            ) * idx_by_2;
        }

        if (betax[IDX(5,j,k)] > 0.0 ) {
          Dxu[IDX(5,j,k)] = (-  3.0 * u[IDX(4,j,k)]
                             - 10.0 * u[IDX(5,j,k)]
                             + 18.0 * u[IDX(6,j,k)]
                             -  6.0 * u[IDX(7,j,k)]
                             +        u[IDX(8,j,k)]
                           ) * idx_by_12;
        }
        else {
          Dxu[IDX(5,j,k)] = (           u[IDX(3,j,k)]
                               -  4.0 * u[IDX(4,j,k)]
                               +  3.0 * u[IDX(5,j,k)]
                            ) * idx_by_2;
        }

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        if ( betax[IDX(ie-3,j,k)] > 0.0 ) {
          Dxu[IDX(ie-3,j,k)] = (  - 3.0 * u[IDX(ie-3,j,k)]
                                  + 4.0 * u[IDX(ie-2,j,k)]
                                  -       u[IDX(ie-1,j,k)]
                               ) * idx_by_2;
        }
        else {
          Dxu[IDX(ie-3,j,k)] = ( -   u[IDX(ie-6,j,k)]
                            +  6.0 * u[IDX(ie-5,j,k)]
                            - 18.0 * u[IDX(ie-4,j,k)]
                            + 10.0 * u[IDX(ie-3  ,j,k)]
                            +  3.0 * u[IDX(ie-2,j,k)]
                          ) * idx_by_12;
        }

        if (betax[IDX(ie-2,j,k)] > 0.0 ) {
          Dxu[IDX(ie-2,j,k)] = (  -  u[IDX(ie-3,j,k)]
                                  +  u[IDX(ie-1,j,k)]
                               ) * idx_by_2;
        }
        else {
          Dxu[IDX(ie-2,j,k)] = (     u[IDX(ie-4,j,k)]
                             - 4.0 * u[IDX(ie-3,j,k)]
                             + 3.0 * u[IDX(ie-2,j,k)]
                               ) * idx_by_2;
        }

        Dxu[IDX(ie-1,j,k)] = (          u[IDX(ie-3,j,k)]
                                - 4.0 * u[IDX(ie-2,j,k)]
                                + 3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

      }
    }
  }


#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_y(double * const  Dyu, const double * const  u,
                  const double dy, const unsigned int *sz,
                  const double * const betay, unsigned bflag)
{

  const double idy = 1.0/dy;
  const double idy_by_2 = 0.50 * idy;
  const double idy_by_12 = idy / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        if (betay[pp] > 0.0 ) {
          Dyu[pp] = ( -  3.0 * u[pp-nx]
                      - 10.0 * u[pp]
                      + 18.0 * u[pp+nx]
                      -  6.0 * u[pp+2*nx]
                      +        u[pp+3*nx]
                    ) * idy_by_12;
        }
        else {
          Dyu[pp] = ( - u[pp-3*nx]
                      +  6.0 * u[pp-2*nx]
                      - 18.0 * u[pp-nx]
                      + 10.0 * u[pp]
                      +  3.0 * u[pp+nx]
                    ) * idy_by_12;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        Dyu[IDX(i,3,k)] = ( -  3.0 * u[IDX(i,3,k)]
                            +  4.0 * u[IDX(i,4,k)]
                            -        u[IDX(i,5,k)]
                          ) * idy_by_2;

        if (betay[IDX(i,4,k)] > 0.0) {
          Dyu[IDX(i,4,k)] = ( -  3.0 * u[IDX(i,4,k)]
                              +  4.0 * u[IDX(i,5,k)]
                              -        u[IDX(i,6,k)]
                            ) * idy_by_2;
        }
        else {
          Dyu[IDX(i,4,k)] = ( -         u[IDX(i,3,k)]
                               +        u[IDX(i,5,k)]
                            ) * idy_by_2;
        }

        if (betay[IDX(i,5,k)] > 0.0 ) {
          Dyu[IDX(i,5,k)] = ( -  3.0 * u[IDX(i,4,k)]
                              - 10.0 * u[IDX(i,5,k)]
                              + 18.0 * u[IDX(i,6,k)]
                              -  6.0 * u[IDX(i,7,k)]
                             +         u[IDX(i,8,k)]
                           ) * idy_by_12;
        }
        else {
          Dyu[IDX(i,5,k)] = (           u[IDX(i,3,k)]
                               -  4.0 * u[IDX(i,4,k)]
                               +  3.0 * u[IDX(i,5,k)]
                            ) * idy_by_2;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        if ( betay[IDX(i,je-3,k)] > 0.0 ) {
          Dyu[IDX(i,je-3,k)] = (  - 3.0 * u[IDX(i,je-3,k)]
                                + 4.0 * u[IDX(i,je-2,k)]
                                -       u[IDX(i,je-1,k)]
                             ) * idy_by_2;
        }
        else {
          Dyu[IDX(i,je-3,k)] = ( -   u[IDX(i,je-6,k)]
                            +  6.0 * u[IDX(i,je-5,k)]
                            - 18.0 * u[IDX(i,je-4,k)]
                            + 10.0 * u[IDX(i,je-3,k)]
                            +  3.0 * u[IDX(i,je-2,k)]
                          ) * idy_by_12;
        }

        if (betay[IDX(i,je-2,k)] > 0.0 ) {
          Dyu[IDX(i,je-2,k)] = (  -  u[IDX(i,je-3,k)]
                                  +  u[IDX(i,je-1,k)]
                               ) * idy_by_2;
        }
        else {
          Dyu[IDX(i,je-2,k)] = (     u[IDX(i,je-4,k)]
                             - 4.0 * u[IDX(i,je-3,k)]
                             + 3.0 * u[IDX(i,je-2,k)]
                               ) * idy_by_2;
        }

        Dyu[IDX(i,je-1,k)] = (          u[IDX(i,je-3,k)]
                                - 4.0 * u[IDX(i,je-2,k)]
                                + 3.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_2;

      }
    }
  }


#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_z(double * const  Dzu, const double * const  u,
                  const double dz, const unsigned int *sz,
                  const double * const betaz, unsigned bflag)
{

  const double idz = 1.0/dz;
  const double idz_by_2 = 0.50 * idz;
  const double idz_by_12 = idz / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        if (betaz[pp] > 0.0 ) {
          Dzu[pp] = ( -  3.0 * u[pp-n]
                      - 10.0 * u[pp]
                      + 18.0 * u[pp+n]
                      -  6.0 * u[pp+2*n]
                      +        u[pp+3*n]
                    ) * idz_by_12;
        }
        else {
          Dzu[pp] = ( -        u[pp-3*n]
                      +  6.0 * u[pp-2*n]
                      - 18.0 * u[pp-n]
                      + 10.0 * u[pp]
                      +  3.0 * u[pp+n]
                    ) * idz_by_12;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i,j,3)] = ( -  3.0 * u[IDX(i,j,3)]
                            +  4.0 * u[IDX(i,j,4)]
                            -        u[IDX(i,j,5)]
                          ) * idz_by_2;

        if (betaz[IDX(i,j,4)] > 0.0) {
          Dzu[IDX(i,j,4)] = ( -  3.0 * u[IDX(i,j,4)]
                              +  4.0 * u[IDX(i,j,5)]
                              -        u[IDX(i,j,6)]
                            ) * idz_by_2;
        }
        else {
          Dzu[IDX(i,j,4)] = ( -         u[IDX(i,j,3)]
                               +        u[IDX(i,j,5)]
                            ) * idz_by_2;
        }

        if (betaz[IDX(i,j,5)] > 0.0 ) {
          Dzu[IDX(i,j,5)] = ( -  3.0 * u[IDX(i,j,4)]
                              - 10.0 * u[IDX(i,j,5)]
                              + 18.0 * u[IDX(i,j,6)]
                              -  6.0 * u[IDX(i,j,7)]
                             +         u[IDX(i,j,8)]
                           ) * idz_by_12;
        }
        else {
          Dzu[IDX(i,j,5)] = (           u[IDX(i,j,3)]
                               -  4.0 * u[IDX(i,j,4)]
                               +  3.0 * u[IDX(i,j,5)]
                            ) * idz_by_2;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        if ( betaz[IDX(i,j,ke-3)] > 0.0 ) {
          Dzu[IDX(i,j,ke-3)] = (  - 3.0 * u[IDX(i,j,ke-3)]
                                  + 4.0 * u[IDX(i,j,ke-2)]
                                  -       u[IDX(i,j,ke-1)]
                               ) * idz_by_2;
        }
        else {
          Dzu[IDX(i,j,ke-3)] = ( -   u[IDX(i,j,ke-6)]
                            +  6.0 * u[IDX(i,j,ke-5)]
                            - 18.0 * u[IDX(i,j,ke-4)]
                            + 10.0 * u[IDX(i,j,ke-3)]
                            +  3.0 * u[IDX(i,j,ke-2)]
                          ) * idz_by_12;
        }

        if (betaz[IDX(i,j,ke-2)] > 0.0 ) {
          Dzu[IDX(i,j,ke-2)] = (  -  u[IDX(i,j,ke-3)]
                                  +  u[IDX(i,j,ke-1)]
                               ) * idz_by_2;
        }
        else {
          Dzu[IDX(i,j,ke-2)] = (     u[IDX(i,j,ke-4)]
                             - 4.0 * u[IDX(i,j,ke-3)]
                             + 3.0 * u[IDX(i,j,ke-2)]
                               ) * idz_by_2;
        }

        Dzu[IDX(i,j,ke-1)] = (          u[IDX(i,j,ke-3)]
                                - 4.0 * u[IDX(i,j,ke-2)]
                                + 3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_2;
      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv42_x(double * const  Du, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_6_dx = -1.0 / 64.0 / dx;

  double smr3=59.0/48.0*64*dx;
  double smr2=43.0/48.0*64*dx;
  double smr1=49.0/48.0*64*dx;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
       for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_6_dx *
                         (
                         -      u[pp-3]
                         +  6.0*u[pp-2]
                         - 15.0*u[pp-1]
                         + 20.0*u[pp]
                         - 15.0*u[pp+1]
                         +  6.0*u[pp+2]
                         -      u[pp+3]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Du[IDX(3,j,k)] =  (      u[IDX(6,j,k)]
                           - 3.0*u[IDX(5,j,k)]
                           + 3.0*u[IDX(4,j,k)]
                           -     u[IDX(3,j,k)]
                          )/smr3;
        Du[IDX(4,j,k)] =  (
                                 u[IDX(7,j,k)]
                          -  6.0*u[IDX(6,j,k)]
                          + 12.0*u[IDX(5,j,k)]
                          - 10.0*u[IDX(4,j,k)]
                          +  3.0*u[IDX(3,j,k)]
                          )/smr2;
        Du[IDX(5,j,k)] =  (
                                 u[IDX(8,j,k)]
                          -  6.0*u[IDX(7,j,k)]
                          + 15.0*u[IDX(6,j,k)]
                          - 19.0*u[IDX(5,j,k)]
                          + 12.0*u[IDX(4,j,k)]
                          -  3.0*u[IDX(3,j,k)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
         Du[IDX(ie-3,j,k)] = (
                                 u[IDX(ie-6,j,k)]
                          -  6.0*u[IDX(ie-5,j,k)]
                          + 15.0*u[IDX(ie-4,j,k)]
                          - 19.0*u[IDX(ie-3,j,k)]
                          + 12.0*u[IDX(ie-2,j,k)]
                          -  3.0*u[IDX(ie-1,j,k)]
                           )/spr1;

         Du[IDX(ie-2,j,k)] = (
                                 u[IDX(ie-5,j,k)]
                          -  6.0*u[IDX(ie-4,j,k)]
                          + 12.0*u[IDX(ie-3,j,k)]
                          - 10.0*u[IDX(ie-2,j,k)]
                          +  3.0*u[IDX(ie-1,j,k)]
                           )/spr2;

         Du[IDX(ie-1,j,k)] = (
                                 u[IDX(ie-4,j,k)]
                          -  3.0*u[IDX(ie-3,j,k)]
                          +  3.0*u[IDX(ie-2,j,k)]
                          -      u[IDX(ie-1,j,k)]
                           )/spr3;
      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}



/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/

void ko_deriv42_y(double * const  Du, const double * const  u,
                const double dy, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_6_dy = -1.0 / 64.0 / dy;

  double smr3=59.0/48.0*64*dy;
  double smr2=43.0/48.0*64*dy;
  double smr1=49.0/48.0*64*dy;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
       for (int j = jb; j < je; j++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_6_dy *
                         (
                         -      u[pp-3*nx]
                         +  6.0*u[pp-2*nx]
                         - 15.0*u[pp-nx]
                         + 20.0*u[pp]
                         - 15.0*u[pp+nx]
                         +  6.0*u[pp+2*nx]
                         -      u[pp+3*nx]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,3,k)] =  (      u[IDX(i,6,k)]
                           - 3.0*u[IDX(i,5,k)]
                           + 3.0*u[IDX(i,4,k)]
                           -     u[IDX(i,3,k)]
                          )/smr3;
        Du[IDX(i,4,k)] =  (
                                 u[IDX(i,7,k)]
                          -  6.0*u[IDX(i,6,k)]
                          + 12.0*u[IDX(i,5,k)]
                          - 10.0*u[IDX(i,4,k)]
                          +  3.0*u[IDX(i,3,k)]
                          )/smr2;
        Du[IDX(i,5,k)] =  (
                                 u[IDX(i,8,k)]
                          -  6.0*u[IDX(i,7,k)]
                          + 15.0*u[IDX(i,6,k)]
                          - 19.0*u[IDX(i,5,k)]
                          + 12.0*u[IDX(i,4,k)]
                          -  3.0*u[IDX(i,3,k)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,je-3,k)] = (
                                 u[IDX(i,je-6,k)]
                          -  6.0*u[IDX(i,je-5,k)]
                          + 15.0*u[IDX(i,je-4,k)]
                          - 19.0*u[IDX(i,je-3,k)]
                          + 12.0*u[IDX(i,je-2,k)]
                          -  3.0*u[IDX(i,je-1,k)]
                           )/spr1;

       Du[IDX(i,je-2,k)] = (
                                 u[IDX(i,je-5,k)]
                          -  6.0*u[IDX(i,je-4,k)]
                          + 12.0*u[IDX(i,je-3,k)]
                          - 10.0*u[IDX(i,je-2,k)]
                          +  3.0*u[IDX(i,je-1,k)]
                           )/spr2;

       Du[IDX(i,je-1,k)] = (
                                 u[IDX(i,je-4,k)]
                          -  3.0*u[IDX(i,je-3,k)]
                          +  3.0*u[IDX(i,je-2,k)]
                          -      u[IDX(i,je-1,k)]
                           )/spr3;

      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}



/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/

void ko_deriv42_z(double * const  Du, const double * const  u,
                const double dz, const unsigned *sz, unsigned bflag)
{

  double pre_factor_6_dz = -1.0 / 64.0 / dz;

  double smr3=59.0/48.0*64*dz;
  double smr2=43.0/48.0*64*dz;
  double smr1=49.0/48.0*64*dz;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
       for (int k = kb; k < ke; k++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_6_dz *
                         (
                         -      u[pp-3*n]
                         +  6.0*u[pp-2*n]
                         - 15.0*u[pp-n]
                         + 20.0*u[pp]
                         - 15.0*u[pp+n]
                         +  6.0*u[pp+2*n]
                         -      u[pp+3*n]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,j,3)] =  (      u[IDX(i,j,6)]
                           - 3.0*u[IDX(i,j,5)]
                           + 3.0*u[IDX(i,j,4)]
                           -     u[IDX(i,j,3)]
                          )/smr3;
        Du[IDX(i,j,4)] =  (
                                 u[IDX(i,j,7)]
                          -  6.0*u[IDX(i,j,6)]
                          + 12.0*u[IDX(i,j,5)]
                          - 10.0*u[IDX(i,j,4)]
                          +  3.0*u[IDX(i,j,3)]
                          )/smr2;
        Du[IDX(i,j,5)] =  (
                                 u[IDX(i,j,8)]
                          -  6.0*u[IDX(i,j,7)]
                          + 15.0*u[IDX(i,j,6)]
                          - 19.0*u[IDX(i,j,5)]
                          + 12.0*u[IDX(i,j,4)]
                          -  3.0*u[IDX(i,j,3)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,j,ke-3)] = (
                                 u[IDX(i,j,ke-6)]
                          -  6.0*u[IDX(i,j,ke-5)]
                          + 15.0*u[IDX(i,j,ke-4)]
                          - 19.0*u[IDX(i,j,ke-3)]
                          + 12.0*u[IDX(i,j,ke-2)]
                          -  3.0*u[IDX(i,j,ke-1)]
                           )/spr1;

         Du[IDX(i,j,ke-2)] = (
                                 u[IDX(i,j,ke-5)]
                          -  6.0*u[IDX(i,j,ke-4)]
                          + 12.0*u[IDX(i,j,ke-3)]
                          - 10.0*u[IDX(i,j,ke-2)]
                          +  3.0*u[IDX(i,j,ke-1)]
                           )/spr2;

         Du[IDX(i,j,ke-1)] = (
                                 u[IDX(i,j,ke-4)]
                          -  3.0*u[IDX(i,j,ke-3)]
                          +  3.0*u[IDX(i,j,ke-2)]
                          -      u[IDX(i,j,ke-1)]
                           )/spr3;
      }
    }
  }

#ifdef DEBUG_DERIVS_COMP
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}



/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void cpy_unzip_padd(double * const  Du, const double * const  u,const unsigned int *sz, unsigned bflag)
{
    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];

    for(unsigned int k=0;k<sz[2];k++)
      for(unsigned int j=0;j<sz[1];j++)
        for(unsigned int i=0;i<sz[0];i++)
            if((i<3||i>=sz[0]-3) || (j<3||j>=sz[0]-3)|| (k<3||k>=sz[2]-3))
                Du[IDX(i,j,k)]=u[IDX(i,j,k)];


}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void disstvb3_x(double * const  Du, const double * const  u,
                const double * const lam,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  const double pre_factor = -1.0 / 12.0 / dx;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const double smr3=59.0/48.0*64*dx;
  const double smr2=43.0/48.0*64*dx;
  const double smr1=49.0/48.0*64*dx;
  const double spr3=smr3;
  const double spr2=smr2;
  const double spr1=smr1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
       for (int i = ib+1; i < ie-1; i++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor *
                         (
                               lam[pp-2]*u[pp-2]
                         - 4.0*lam[pp-1]*u[pp-1]
                         + 6.0*lam[pp]*u[pp]
                         - 4.0*lam[pp+1]*u[pp+1]
                         +     lam[pp+2]*u[pp+2]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Du[IDX(3,j,k)] =  (      u[IDX(6,j,k)]
                           - 3.0*u[IDX(5,j,k)]
                           + 3.0*u[IDX(4,j,k)]
                           -     u[IDX(3,j,k)]
                          )/smr3;
        Du[IDX(4,j,k)] =  (
                                 u[IDX(7,j,k)]
                          -  6.0*u[IDX(6,j,k)]
                          + 12.0*u[IDX(5,j,k)]
                          - 10.0*u[IDX(4,j,k)]
                          +  3.0*u[IDX(3,j,k)]
                          )/smr2;
        Du[IDX(5,j,k)] =  (
                                 u[IDX(8,j,k)]
                          -  6.0*u[IDX(7,j,k)]
                          + 15.0*u[IDX(6,j,k)]
                          - 19.0*u[IDX(5,j,k)]
                          + 12.0*u[IDX(4,j,k)]
                          -  3.0*u[IDX(3,j,k)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
         Du[IDX(ie-3,j,k)] = (
                                 u[IDX(ie-6,j,k)]
                          -  6.0*u[IDX(ie-5,j,k)]
                          + 15.0*u[IDX(ie-4,j,k)]
                          - 19.0*u[IDX(ie-3,j,k)]
                          + 12.0*u[IDX(ie-2,j,k)]
                          -  3.0*u[IDX(ie-1,j,k)]
                           )/spr1;

         Du[IDX(ie-2,j,k)] = (
                                 u[IDX(ie-5,j,k)]
                          -  6.0*u[IDX(ie-4,j,k)]
                          + 12.0*u[IDX(ie-3,j,k)]
                          - 10.0*u[IDX(ie-2,j,k)]
                          +  3.0*u[IDX(ie-1,j,k)]
                           )/spr2;

         Du[IDX(ie-1,j,k)] = (
                                 u[IDX(ie-4,j,k)]
                          -  3.0*u[IDX(ie-3,j,k)]
                          +  3.0*u[IDX(ie-2,j,k)]
                          -      u[IDX(ie-1,j,k)]
                           )/spr3;
      }
    }
  }

}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void disstvb3_y(double * const  Du, const double * const  u,
                const double * const lam,
                const double dy, const unsigned int *sz, unsigned bflag)
{

  const double pre_factor = -1.0 / 12.0 / dy;

  double smr3=59.0/48.0*64*dy;
  double smr2=43.0/48.0*64*dy;
  double smr1=49.0/48.0*64*dy;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
       for (int j = jb+1; j < je-1; j++) {
         const int pp = IDX(i,j,k);
         Du[pp] = pre_factor *
                         (
                               lam[pp-2*nx]*u[pp-2*nx]
                         - 4.0*lam[pp-nx]*u[pp-nx]
                         + 6.0*lam[pp]*u[pp]
                         - 4.0*lam[pp+nx]*u[pp+nx]
                         +     lam[pp+2*nx]*u[pp+2*nx]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,3,k)] =  (      u[IDX(i,6,k)]
                           - 3.0*u[IDX(i,5,k)]
                           + 3.0*u[IDX(i,4,k)]
                           -     u[IDX(i,3,k)]
                          )/smr3;
        Du[IDX(i,4,k)] =  (
                                 u[IDX(i,7,k)]
                          -  6.0*u[IDX(i,6,k)]
                          + 12.0*u[IDX(i,5,k)]
                          - 10.0*u[IDX(i,4,k)]
                          +  3.0*u[IDX(i,3,k)]
                          )/smr2;
        Du[IDX(i,5,k)] =  (
                                 u[IDX(i,8,k)]
                          -  6.0*u[IDX(i,7,k)]
                          + 15.0*u[IDX(i,6,k)]
                          - 19.0*u[IDX(i,5,k)]
                          + 12.0*u[IDX(i,4,k)]
                          -  3.0*u[IDX(i,3,k)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,je-3,k)] = (
                                 u[IDX(i,je-6,k)]
                          -  6.0*u[IDX(i,je-5,k)]
                          + 15.0*u[IDX(i,je-4,k)]
                          - 19.0*u[IDX(i,je-3,k)]
                          + 12.0*u[IDX(i,je-2,k)]
                          -  3.0*u[IDX(i,je-1,k)]
                           )/spr1;

       Du[IDX(i,je-2,k)] = (
                                 u[IDX(i,je-5,k)]
                          -  6.0*u[IDX(i,je-4,k)]
                          + 12.0*u[IDX(i,je-3,k)]
                          - 10.0*u[IDX(i,je-2,k)]
                          +  3.0*u[IDX(i,je-1,k)]
                           )/spr2;

       Du[IDX(i,je-1,k)] = (
                                 u[IDX(i,je-4,k)]
                          -  3.0*u[IDX(i,je-3,k)]
                          +  3.0*u[IDX(i,je-2,k)]
                          -      u[IDX(i,je-1,k)]
                           )/spr3;

      }
    }
  }

}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void disstvb3_z(double * const  Du, const double * const  u,
                const double * const lam,
                const double dz, const unsigned *sz, unsigned bflag)
{

  const double pre_factor = -1.0 / 12.0 / dz;

  double smr3=59.0/48.0*64*dz;
  double smr2=43.0/48.0*64*dz;
  double smr1=49.0/48.0*64*dz;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
       for (int k = kb+1; k < ke-1; k++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor *
                         (
                               lam[pp-2*n]*u[pp-2*n]
                         - 4.0*lam[pp-n]*u[pp-n]
                         + 6.0*lam[pp]*u[pp]
                         - 4.0*lam[pp+n]*u[pp+n]
                         +     lam[pp+2*n]*u[pp+2*n]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,j,3)] =  (      u[IDX(i,j,6)]
                           - 3.0*u[IDX(i,j,5)]
                           + 3.0*u[IDX(i,j,4)]
                           -     u[IDX(i,j,3)]
                          )/smr3;
        Du[IDX(i,j,4)] =  (
                                 u[IDX(i,j,7)]
                          -  6.0*u[IDX(i,j,6)]
                          + 12.0*u[IDX(i,j,5)]
                          - 10.0*u[IDX(i,j,4)]
                          +  3.0*u[IDX(i,j,3)]
                          )/smr2;
        Du[IDX(i,j,5)] =  (
                                 u[IDX(i,j,8)]
                          -  6.0*u[IDX(i,j,7)]
                          + 15.0*u[IDX(i,j,6)]
                          - 19.0*u[IDX(i,j,5)]
                          + 12.0*u[IDX(i,j,4)]
                          -  3.0*u[IDX(i,j,3)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,j,ke-3)] = (
                                 u[IDX(i,j,ke-6)]
                          -  6.0*u[IDX(i,j,ke-5)]
                          + 15.0*u[IDX(i,j,ke-4)]
                          - 19.0*u[IDX(i,j,ke-3)]
                          + 12.0*u[IDX(i,j,ke-2)]
                          -  3.0*u[IDX(i,j,ke-1)]
                           )/spr1;

         Du[IDX(i,j,ke-2)] = (
                                 u[IDX(i,j,ke-5)]
                          -  6.0*u[IDX(i,j,ke-4)]
                          + 12.0*u[IDX(i,j,ke-3)]
                          - 10.0*u[IDX(i,j,ke-2)]
                          +  3.0*u[IDX(i,j,ke-1)]
                           )/spr2;

         Du[IDX(i,j,ke-1)] = (
                                 u[IDX(i,j,ke-4)]
                          -  3.0*u[IDX(i,j,ke-3)]
                          +  3.0*u[IDX(i,j,ke-2)]
                          -      u[IDX(i,j,ke-1)]
                           )/spr3;
      }
    }
  }

}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void disstvb5_x(double * const  Du, const double * const  u,
                const double * const lam,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_6_dx = 2.0 / 75.0 / dx;

  double smr3=59.0/48.0*64*dx;
  double smr2=43.0/48.0*64*dx;
  double smr1=49.0/48.0*64*dx;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      int pp = IDX(ib,j,k);
      Du[pp] =  (
                                 lam[pp+3]*u[pp+3]
                          -  6.0*lam[pp+2]*u[pp+2]
                          + 15.0*lam[pp+1]*u[pp+1]
                          - 19.0*lam[pp]*u[pp]
                          + 12.0*lam[pp-1]*u[pp-1]
                          -  3.0*lam[pp-2]*u[pp-2]
                          ) / smr1;
        
      for (int i = ib+1; i < ie-1; i++) {
        pp = IDX(i,j,k);
        Du[pp] = pre_factor_6_dx *
                 ( 
                          std::max(lam[pp+3],lam[pp+2])*(u[pp+3] - u[pp+2])
                   -  5.0*std::max(lam[pp+2],lam[pp+1])*(u[pp+2] - u[pp+1])
                   + 10.0*std::max(lam[pp+1],lam[pp])*(u[pp+1] - u[pp])
                   - 10.0*std::max(lam[pp],lam[pp-1])*(u[pp] - u[pp-1])
                   +  5.0*std::max(lam[pp-1],lam[pp-2])*(u[pp-1] - u[pp-2])
                   -      std::max(lam[pp-2],lam[pp-3])*(u[pp-2] - u[pp-3])
                 );
      }
       
      pp = IDX(ie-1,j,k);
      Du[pp] = (
                                 lam[pp-4]*u[pp-4]
                          -  6.0*lam[pp-3]*u[pp-3]
                          + 15.0*lam[pp-2]*u[pp-2]
                          - 19.0*lam[pp-1]*u[pp-1]
                          + 12.0*lam[pp]*u[pp]
                          -  3.0*lam[pp+1]*u[pp+1]
                           ) / spr1;
       
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Du[IDX(3,j,k)] =  (      u[IDX(6,j,k)]
                           - 3.0*u[IDX(5,j,k)]
                           + 3.0*u[IDX(4,j,k)]
                           -     u[IDX(3,j,k)]
                          )/smr3;
        Du[IDX(4,j,k)] =  (
                                 u[IDX(7,j,k)]
                          -  6.0*u[IDX(6,j,k)]
                          + 12.0*u[IDX(5,j,k)]
                          - 10.0*u[IDX(4,j,k)]
                          +  3.0*u[IDX(3,j,k)]
                          )/smr2;
        Du[IDX(5,j,k)] =  (
                                 u[IDX(8,j,k)]
                          -  6.0*u[IDX(7,j,k)]
                          + 15.0*u[IDX(6,j,k)]
                          - 19.0*u[IDX(5,j,k)]
                          + 12.0*u[IDX(4,j,k)]
                          -  3.0*u[IDX(3,j,k)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
         Du[IDX(ie-3,j,k)] = (
                                 u[IDX(ie-6,j,k)]
                          -  6.0*u[IDX(ie-5,j,k)]
                          + 15.0*u[IDX(ie-4,j,k)]
                          - 19.0*u[IDX(ie-3,j,k)]
                          + 12.0*u[IDX(ie-2,j,k)]
                          -  3.0*u[IDX(ie-1,j,k)]
                           )/spr1;

         Du[IDX(ie-2,j,k)] = (
                                 u[IDX(ie-5,j,k)]
                          -  6.0*u[IDX(ie-4,j,k)]
                          + 12.0*u[IDX(ie-3,j,k)]
                          - 10.0*u[IDX(ie-2,j,k)]
                          +  3.0*u[IDX(ie-1,j,k)]
                           )/spr2;

         Du[IDX(ie-1,j,k)] = (
                                 u[IDX(ie-4,j,k)]
                          -  3.0*u[IDX(ie-3,j,k)]
                          +  3.0*u[IDX(ie-2,j,k)]
                          -      u[IDX(ie-1,j,k)]
                           )/spr3;
      }
    }
  }

}



/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void disstvb5_y(double * const  Du, const double * const  u,
                const double * const lam,
                const double dy, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_6_dy = 2.0 / 75.0 / dy;

  double smr3=59.0/48.0*64*dy;
  double smr2=43.0/48.0*64*dy;
  double smr1=49.0/48.0*64*dy;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      int pp = IDX(i,jb,k);

      Du[pp] =  (
                                 lam[pp+3*nx]*u[pp+3*nx]
                          -  6.0*lam[pp+2*nx]*u[pp+2*nx]
                          + 15.0*lam[pp+nx]*u[pp+nx]
                          - 19.0*lam[pp]*u[pp]
                          + 12.0*lam[pp-nx]*u[pp-nx]
                          -  3.0*lam[pp-2*nx]*u[pp-2*nx]
                          )/smr1;

      for (int j = jb+1; j < je-1; j++) {
        pp = IDX(i,j,k);
        Du[pp] = pre_factor_6_dy *
                 (
                          std::max(lam[pp+3*nx],lam[pp+2*nx])*(u[pp+3*nx] - u[pp+2*nx])
                   -  5.0*std::max(lam[pp+2*nx],lam[pp+nx])*(u[pp+2*nx] - u[pp+nx])
                   + 10.0*std::max(lam[pp+nx],lam[pp])*(u[pp+nx] - u[pp])
                   - 10.0*std::max(lam[pp],lam[pp-nx])*(u[pp] - u[pp-nx])
                   +  5.0*std::max(lam[pp-nx],lam[pp-2*nx])*(u[pp-nx] - u[pp-2*nx])
                   -      std::max(lam[pp-2*nx],lam[pp-3*nx])*(u[pp-2*nx] - u[pp-3*nx])
                 );
      }

      pp = IDX(i,je-1,k);
      Du[pp] = (
                                 lam[pp-4*nx]*u[pp-4*nx]
                          -  6.0*lam[pp-3*nx]*u[pp-3*nx]
                          + 15.0*lam[pp-2*nx]*u[pp-2*nx]
                          - 19.0*lam[pp-nx]*u[pp-nx]
                          + 12.0*lam[pp]*u[pp]
                          -  3.0*lam[pp+nx]*u[pp+nx]
                           ) / spr1;
       
       
       
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,3,k)] =  (      u[IDX(i,6,k)]
                           - 3.0*u[IDX(i,5,k)]
                           + 3.0*u[IDX(i,4,k)]
                           -     u[IDX(i,3,k)]
                          )/smr3;
        Du[IDX(i,4,k)] =  (
                                 u[IDX(i,7,k)]
                          -  6.0*u[IDX(i,6,k)]
                          + 12.0*u[IDX(i,5,k)]
                          - 10.0*u[IDX(i,4,k)]
                          +  3.0*u[IDX(i,3,k)]
                          )/smr2;
        Du[IDX(i,5,k)] =  (
                                 u[IDX(i,8,k)]
                          -  6.0*u[IDX(i,7,k)]
                          + 15.0*u[IDX(i,6,k)]
                          - 19.0*u[IDX(i,5,k)]
                          + 12.0*u[IDX(i,4,k)]
                          -  3.0*u[IDX(i,3,k)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,je-3,k)] = (
                                 u[IDX(i,je-6,k)]
                          -  6.0*u[IDX(i,je-5,k)]
                          + 15.0*u[IDX(i,je-4,k)]
                          - 19.0*u[IDX(i,je-3,k)]
                          + 12.0*u[IDX(i,je-2,k)]
                          -  3.0*u[IDX(i,je-1,k)]
                           )/spr1;

       Du[IDX(i,je-2,k)] = (
                                 u[IDX(i,je-5,k)]
                          -  6.0*u[IDX(i,je-4,k)]
                          + 12.0*u[IDX(i,je-3,k)]
                          - 10.0*u[IDX(i,je-2,k)]
                          +  3.0*u[IDX(i,je-1,k)]
                           )/spr2;

       Du[IDX(i,je-1,k)] = (
                                 u[IDX(i,je-4,k)]
                          -  3.0*u[IDX(i,je-3,k)]
                          +  3.0*u[IDX(i,je-2,k)]
                          -      u[IDX(i,je-1,k)]
                           )/spr3;

      }
    }
  }

}



/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/

void disstvb5_z(double * const  Du, const double * const  u,
                const double * const lam,
                const double dz, const unsigned *sz, unsigned bflag)
{

  double pre_factor_6_dz = -1.0 / 64.0 / dz;

  double smr3=59.0/48.0*64*dz;
  double smr2=43.0/48.0*64*dz;
  double smr1=49.0/48.0*64*dz;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {

      int pp = IDX(i,j,kb);  
      Du[pp] =  (
                                 lam[pp+3*n]*u[pp+3*n]
                          -  6.0*lam[pp+2*n]*u[pp+2*n]
                          + 15.0*lam[pp+n]*u[pp+n]
                          - 19.0*lam[pp]*u[pp]
                          + 12.0*lam[pp-n]*u[pp-n]
                          -  3.0*lam[pp-2*n]*u[pp-2*n]
                          )/smr1;

        
      for (int k = kb+1; k < ke-1; k++) {
        pp = IDX(i,j,k);
        Du[pp] = pre_factor_6_dz *
                 (
                          std::max(lam[pp+3*n],lam[pp+2*n])*(u[pp+3*n] - u[pp+2*n])
                   -  5.0*std::max(lam[pp+2*n],lam[pp+n])*(u[pp+2*n] - u[pp+n])
                   + 10.0*std::max(lam[pp+n],lam[pp])*(u[pp+n] - u[pp])
                   - 10.0*std::max(lam[pp],lam[pp-n])*(u[pp] - u[pp-n])
                   +  5.0*std::max(lam[pp-n],lam[pp-2*n])*(u[pp-n] - u[pp-2*n])
                   -      std::max(lam[pp-2*n],lam[pp-3*n])*(u[pp-2*n] - u[pp-3*n])
                 );
      }

      pp = IDX(i,j,ke-1);
      Du[pp] = (
                                 lam[pp-4*n]*u[pp-4*n]
                          -  6.0*lam[pp-3*n]*u[pp-3*n]
                          + 15.0*lam[pp-2*n]*u[pp-2*n]
                          - 19.0*lam[pp-n]*u[pp-n]
                          + 12.0*lam[pp]*u[pp]
                          -  3.0*lam[pp+n]*u[pp+n]
                           ) / spr1;
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,j,3)] =  (      u[IDX(i,j,6)]
                           - 3.0*u[IDX(i,j,5)]
                           + 3.0*u[IDX(i,j,4)]
                           -     u[IDX(i,j,3)]
                          )/smr3;
        Du[IDX(i,j,4)] =  (
                                 u[IDX(i,j,7)]
                          -  6.0*u[IDX(i,j,6)]
                          + 12.0*u[IDX(i,j,5)]
                          - 10.0*u[IDX(i,j,4)]
                          +  3.0*u[IDX(i,j,3)]
                          )/smr2;
        Du[IDX(i,j,5)] =  (
                                 u[IDX(i,j,8)]
                          -  6.0*u[IDX(i,j,7)]
                          + 15.0*u[IDX(i,j,6)]
                          - 19.0*u[IDX(i,j,5)]
                          + 12.0*u[IDX(i,j,4)]
                          -  3.0*u[IDX(i,j,3)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,j,ke-3)] = (
                                 u[IDX(i,j,ke-6)]
                          -  6.0*u[IDX(i,j,ke-5)]
                          + 15.0*u[IDX(i,j,ke-4)]
                          - 19.0*u[IDX(i,j,ke-3)]
                          + 12.0*u[IDX(i,j,ke-2)]
                          -  3.0*u[IDX(i,j,ke-1)]
                           )/spr1;

         Du[IDX(i,j,ke-2)] = (
                                 u[IDX(i,j,ke-5)]
                          -  6.0*u[IDX(i,j,ke-4)]
                          + 12.0*u[IDX(i,j,ke-3)]
                          - 10.0*u[IDX(i,j,ke-2)]
                          +  3.0*u[IDX(i,j,ke-1)]
                           )/spr2;

         Du[IDX(i,j,ke-1)] = (
                                 u[IDX(i,j,ke-4)]
                          -  3.0*u[IDX(i,j,ke-3)]
                          +  3.0*u[IDX(i,j,ke-2)]
                          -      u[IDX(i,j,ke-1)]
                           )/spr3;
      }
    }
  }

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64_x(double * const  Dxu, const double * const  u,
               const double dx, const unsigned int *sz, unsigned bflag)
{

  const double idx = 1.0/dx;
  const double idx_by_12 = idx/12.0;
  const double idx_by_60 = idx/60.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 0;
  const int kb = 0;
  const int ie = sz[0]-3;
  const int je = sz[1]-0;
  const int ke = sz[2]-0;
  const int n=1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        Dxu[pp] = ( - 1.0  * u[pp-3] 
                    + 9.0  * u[pp-2]
                    -45.0  * u[pp-1]
                    +45.0  * u[pp+1]
                    - 9.0  * u[pp+2]
                    + 1.0  * u[pp+3] ) * idx_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Dxu[IDX(3,j,k)] = ( -  25.0 * u[IDX(3,j,k)]
                            +  48.0 * u[IDX(4,j,k)]
                            -  36.0 * u[IDX(5,j,k)]
                            +  16.0 * u[IDX(6,j,k)]
                            -   3.0 * u[IDX(7,j,k)]
                          ) * idx_by_12;
        
        Dxu[IDX(4,j,k)] = ( - 3.0 * u[IDX(3,j,k)]
                            -10.0 * u[IDX(4,j,k)]
                            +18.0 * u[IDX(5,j,k)]
                            - 6.0 * u[IDX(6,j,k)]
                            + 1.0 * u[IDX(7,j,k)]
                          ) * idx_by_12;

        Dxu[IDX(5,j,k)] = ( + 1.0 * u[IDX(3,j,k)]
                            - 8.0 * u[IDX(4,j,k)]
                            + 8.0 * u[IDX(6,j,k)]
                            - 1.0 * u[IDX(7,j,k)]
                          ) * idx_by_12;
        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {


        Dxu[IDX(ie-3,j,k)] = ( +  1.0 * u[IDX(ie-5,j,k)]
                               -  8.0 * u[IDX(ie-4,j,k)]
                               +  8.0 * u[IDX(ie-2,j,k)]
                               -  1.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

        Dxu[IDX(ie-2,j,k)] = ( -  1.0 * u[IDX(ie-5,j,k)]
                               +  6.0 * u[IDX(ie-4,j,k)]
                               - 18.0 * u[IDX(ie-3,j,k)]
                               + 10.0 * u[IDX(ie-2,j,k)]
                               +  3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_12;
                             

        Dxu[IDX(ie-1,j,k)] = ( + 3.0 * u[IDX(ie-5,j,k)]
                              - 16.0 * u[IDX(ie-4,j,k)]
                              + 36.0 * u[IDX(ie-3,j,k)]
                              - 48.0 * u[IDX(ie-2,j,k)]
                              + 25.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

                             

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = 3; k < sz[2]-3; k++) {
      for (int j = 3; j < sz[1]-3; j++) {
        for (int i = 3; i < sz[0]-3; i++) {
          const int pp = IDX(i,j,k);
          if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif


}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64_y(double * const  Dyu, const double * const  u,
               const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy = 1.0/dy;
  const double idy_by_12 = idy / 12.0;
  const double idy_by_60 = idy / 60.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 0;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-0;
  const int n=nx;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        Dyu[pp] = ( - 1.0  * u[pp-3*nx] 
                    + 9.0  * u[pp-2*nx]
                    -45.0  * u[pp-1*nx]
                    +45.0  * u[pp+1*nx]
                    - 9.0  * u[pp+2*nx]
                    + 1.0  * u[pp+3*nx]) * idy_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        Dyu[IDX(i,3,k)] =  ( - 25.0 * u[IDX(i,3,k)]
                            +  48.0 * u[IDX(i,4,k)]
                            -  36.0 * u[IDX(i,5,k)]
                            +  16.0 * u[IDX(i,6,k)]
                            -   3.0 * u[IDX(i,7,k)]
                          ) * idy_by_12;

        Dyu[IDX(i,4,k)] = ( - 3.0 * u[IDX(i,3,k)]
                            -10.0 * u[IDX(i,4,k)]
                            +18.0 * u[IDX(i,5,k)]
                            - 6.0 * u[IDX(i,6,k)]
                            + 1.0 * u[IDX(i,7,k)]
                          ) * idy_by_12;

        Dyu[IDX(i,5,k)] = ( + 1.0 * u[IDX(i,3,k)]
                            - 8.0 * u[IDX(i,4,k)]
                            + 8.0 * u[IDX(i,6,k)]
                            - 1.0 * u[IDX(i,7,k)]
                          ) * idy_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        

        Dyu[IDX(i,je-3,k)] = ( +  1.0 * u[IDX(i,je-5,k)]
                               -  8.0 * u[IDX(i,je-4,k)]
                               +  8.0 * u[IDX(i,je-2,k)]
                               -  1.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_12;


        Dyu[IDX(i,je-2,k)] = ( -  1.0 * u[IDX(i,je-5,k)]
                               +  6.0 * u[IDX(i,je-4,k)]
                               - 18.0 * u[IDX(i,je-3,k)]
                               + 10.0 * u[IDX(i,je-2,k)]
                               +  3.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_12;
                             

        Dyu[IDX(i,je-1,k)] = ( + 3.0 * u[IDX(i,je-5,k)]
                              - 16.0 * u[IDX(i,je-4,k)]
                              + 36.0 * u[IDX(i,je-3,k)]
                              - 48.0 * u[IDX(i,je-2,k)]
                              + 25.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_12;
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
    for (int k = 3; k < sz[2]-3; k++) {
      for (int j = 3; j < sz[1]-3; j++) {
        for (int i = 3; i < sz[0]-3; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64_z(double * const  Dzu, const double * const  u,
               const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz = 1.0/dz;
  const double idz_by_12 = idz / 12.0;
  const double idz_by_60 = idz / 60.0;
  

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-3;

  const int n = nx*ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        Dzu[pp] = ( - 1.0  * u[pp-3*n] 
                    + 9.0  * u[pp-2*n]
                    -45.0  * u[pp-1*n]
                    +45.0  * u[pp+1*n]
                    - 9.0  * u[pp+2*n]
                    + 1.0  * u[pp+3*n]) * idz_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i, j,3)] =  (- 25.0 * u[IDX(i,j,3)]
                            +  48.0 * u[IDX(i,j,4)]
                            -  36.0 * u[IDX(i,j,5)]
                            +  16.0 * u[IDX(i,j,6)]
                            -   3.0 * u[IDX(i,j,7)]
                          ) * idz_by_12;

        Dzu[IDX(i,j,4)] = ( - 3.0 * u[IDX(i,j,3)]
                            -10.0 * u[IDX(i,j,4)]
                            +18.0 * u[IDX(i,j,5)]
                            - 6.0 * u[IDX(i,j,6)]
                            + 1.0 * u[IDX(i,j,7)]
                          ) * idz_by_12;

        Dzu[IDX(i,j,5)] = ( + 1.0 * u[IDX(i,j,3)]
                            - 8.0 * u[IDX(i,j,4)]
                            + 8.0 * u[IDX(i,j,6)]
                            - 1.0 * u[IDX(i,j,7)]
                          ) * idz_by_12;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        

        Dzu[IDX(i,j,ke-3)] = ( +  1.0 * u[IDX(i,j,ke-5)]
                               -  8.0 * u[IDX(i,j,ke-4)]
                               +  8.0 * u[IDX(i,j,ke-2)]
                               -  1.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_12;

        Dzu[IDX(i,j,ke-2)] = ( -  1.0 * u[IDX(i,j,ke-5)]
                               +  6.0 * u[IDX(i,j,ke-4)]
                               - 18.0 * u[IDX(i,j,ke-3)]
                               + 10.0 * u[IDX(i,j,ke-2)]
                               +  3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_12;
                             

        Dzu[IDX(i,j,ke-1)] = ( + 3.0 * u[IDX(i,j,ke-5)]
                              - 16.0 * u[IDX(i,j,ke-4)]
                              + 36.0 * u[IDX(i,j,ke-3)]
                              - 48.0 * u[IDX(i,j,ke-2)]
                              + 25.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_12;

      }
    }

  }

  #ifdef DEBUG_DERIVS_COMP
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64_xx(double * const  DxDxu, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  const double idx_sqrd = 1.0/(dx*dx);
  const double idx_sqrd_by_180 = idx_sqrd / 180.0;
  const double idx_sqrd_by_12  = idx_sqrd / 12.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);

        DxDxu[pp] = (      2.0  * u[pp-3]
                        - 27.0  * u[pp-2]
                        +270.0  * u[pp-1]
                        -490.0  * u[pp]
                        +270.0  * u[pp+1]
                        - 27.0  * u[pp+2]
                        +  2.0  * u[pp+3]
                    ) * idx_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {


        DxDxu[IDX(3,j,k)] = (    35.0 * u[IDX(3,j,k)]
                              - 104.0 * u[IDX(4,j,k)]
                              + 114.0 * u[IDX(5,j,k)]
                              -  56.0 * u[IDX(6,j,k)]
                              +  11.0 * u[IDX(7,j,k)]
                            ) * idx_sqrd_by_12;
        
        DxDxu[IDX(4,j,k)] = (    11.0 * u[IDX(3,j,k)]
                              -  20.0 * u[IDX(4,j,k)]
                              +   6.0 * u[IDX(5,j,k)]
                              +   4.0 * u[IDX(6,j,k)]
                              -   1.0 * u[IDX(7,j,k)]
                            ) * idx_sqrd_by_12;
        
        DxDxu[IDX(5,j,k)] = ( -   1.0 * u[IDX(3,j,k)]
                              +  16.0 * u[IDX(4,j,k)]
                              -  30.0 * u[IDX(5,j,k)]
                              +  16.0 * u[IDX(6,j,k)]
                              -   1.0 * u[IDX(7,j,k)]
                            ) * idx_sqrd_by_12;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {


        DxDxu[IDX(ie-3,j,k)] = ( -   1.0 * u[IDX(ie-5,j,k)]
                                 +  16.0 * u[IDX(ie-4,j,k)]
                                 -  30.0 * u[IDX(ie-3,j,k)]
                                 +  16.0 * u[IDX(ie-2,j,k)]
                                 -   1.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;
        
        DxDxu[IDX(ie-2,j,k)] = ( -   1.0 * u[IDX(ie-5,j,k)]
                                 +   4.0 * u[IDX(ie-4,j,k)]
                                 +   6.0 * u[IDX(ie-3,j,k)]
                                 -  20.0 * u[IDX(ie-2,j,k)]
                                 +  11.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;

        DxDxu[IDX(ie-1,j,k)] = ( +  11.0 * u[IDX(ie-5,j,k)]
                                 -  56.0 * u[IDX(ie-4,j,k)]
                                 + 114.0 * u[IDX(ie-3,j,k)]
                                 - 104.0 * u[IDX(ie-2,j,k)]
                                 +  35.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;

      }
    }
  }

  #ifdef DEBUG_DERIVS_COMP
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64_yy(double * const  DyDyu, const double * const  u,
                const double dy, const unsigned int *sz, unsigned bflag)
{

  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_180 = idy_sqrd / 180.0;
  const double idy_sqrd_by_12  = idy_sqrd /  12.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        DyDyu[pp] = (      2.0  * u[pp-3*nx]
                        - 27.0  * u[pp-2*nx]
                        +270.0  * u[pp-1*nx]
                        -490.0  * u[pp]
                        +270.0  * u[pp+1*nx]
                        - 27.0  * u[pp+2*nx]
                        +  2.0  * u[pp+3*nx]
                    ) * idy_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        
        DyDyu[IDX(i,3,k)] = (    35.0 * u[IDX(i,3,k)]
                              - 104.0 * u[IDX(i,4,k)]
                              + 114.0 * u[IDX(i,5,k)]
                              -  56.0 * u[IDX(i,6,k)]
                              +  11.0 * u[IDX(i,7,k)]
                            ) * idy_sqrd_by_12;
        
        DyDyu[IDX(i,4,k)] = (    11.0 * u[IDX(i,3,k)]
                              -  20.0 * u[IDX(i,4,k)]
                              +   6.0 * u[IDX(i,5,k)]
                              +   4.0 * u[IDX(i,6,k)]
                              -   1.0 * u[IDX(i,7,k)]
                            ) * idy_sqrd_by_12;
        
        DyDyu[IDX(i,5,k)] = ( -   1.0 * u[IDX(i,3,k)]
                              +  16.0 * u[IDX(i,4,k)]
                              -  30.0 * u[IDX(i,5,k)]
                              +  16.0 * u[IDX(i,6,k)]
                              -   1.0 * u[IDX(i,7,k)]
                            ) * idy_sqrd_by_12;


        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        

        DyDyu[IDX(i,je-3,k)] = ( -   1.0 * u[IDX(i,je-5,k)]
                                 +  16.0 * u[IDX(i,je-4,k)]
                                 -  30.0 * u[IDX(i,je-3,k)]
                                 +  16.0 * u[IDX(i,je-2,k)]
                                 -   1.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;
        
        DyDyu[IDX(i,je-2,k)] = ( -   1.0 * u[IDX(i,je-5,k)]
                                 +   4.0 * u[IDX(i,je-4,k)]
                                 +   6.0 * u[IDX(i,je-3,k)]
                                 -  20.0 * u[IDX(i,je-2,k)]
                                 +  11.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;

        DyDyu[IDX(i,je-1,k)] = ( +  11.0 * u[IDX(i,je-5,k)]
                                 -  56.0 * u[IDX(i,je-4,k)]
                                 + 114.0 * u[IDX(i,je-3,k)]
                                 - 104.0 * u[IDX(i,je-2,k)]
                                 +  35.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;

      }
    }
  }

  #ifdef DEBUG_DERIVS_COMP
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64_zz(double * const  DzDzu, const double * const  u,
                const double dz, const unsigned int *sz, unsigned bflag)
{

  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_180 = idz_sqrd / 180.0;
  const double idz_sqrd_by_12  = idz_sqrd /  12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        DzDzu[pp] = (      2.0  * u[pp-3*n]
                        - 27.0  * u[pp-2*n]
                        +270.0  * u[pp-1*n]
                        -490.0  * u[pp]
                        +270.0  * u[pp+1*n]
                        - 27.0  * u[pp+2*n]
                        +  2.0  * u[pp+3*n]
                    ) * idz_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        
        DzDzu[IDX(i,j,3)] = (    35.0 * u[IDX(i,j,3)]
                              - 104.0 * u[IDX(i,j,4)]
                              + 114.0 * u[IDX(i,j,5)]
                              -  56.0 * u[IDX(i,j,6)]
                              +  11.0 * u[IDX(i,j,7)]
                            ) * idz_sqrd_by_12;
        
        DzDzu[IDX(i,j,4)] = (    11.0 * u[IDX(i,j,3)]
                              -  20.0 * u[IDX(i,j,4)]
                              +   6.0 * u[IDX(i,j,5)]
                              +   4.0 * u[IDX(i,j,6)]
                              -   1.0 * u[IDX(i,j,7)]
                            ) * idz_sqrd_by_12;

        DzDzu[IDX(i,j,5)] = ( -   1.0 * u[IDX(i,j,3)]
                              +  16.0 * u[IDX(i,j,4)]
                              -  30.0 * u[IDX(i,j,5)]
                              +  16.0 * u[IDX(i,j,6)]
                              -   1.0 * u[IDX(i,j,7)]
                            ) * idz_sqrd_by_12;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {

        
        DzDzu[IDX(i,j,ke-3)] = ( -   1.0 * u[IDX(i,j,ke-5)]
                                 +  16.0 * u[IDX(i,j,ke-4)]
                                 -  30.0 * u[IDX(i,j,ke-3)]
                                 +  16.0 * u[IDX(i,j,ke-2)]
                                 -   1.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;


        DzDzu[IDX(i,j,ke-2)] = ( -   1.0 * u[IDX(i,j,ke-5)]
                                 +   4.0 * u[IDX(i,j,ke-4)]
                                 +   6.0 * u[IDX(i,j,ke-3)]
                                 -  20.0 * u[IDX(i,j,ke-2)]
                                 +  11.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;

        DzDzu[IDX(i,j,ke-1)] = ( +  11.0 * u[IDX(i,j,ke-5)]
                                 -  56.0 * u[IDX(i,j,ke-4)]
                                 + 114.0 * u[IDX(i,j,ke-3)]
                                 - 104.0 * u[IDX(i,j,ke-2)]
                                 +  35.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;


      }
    }
  }

  #ifdef DEBUG_DERIVS_COMP
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64adv_x(double * const  Dxu, const double * const  u,
                  const double dx, const unsigned int *sz,
                  const double * const betax, unsigned bflag)
{

    // hardcoded for not to use adv derivs, 
    deriv64_x(Dxu,u,dx,sz,bflag);
  

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64adv_y(double * const  Dyu, const double * const  u,
                  const double dy, const unsigned int *sz,
                  const double * const betay, unsigned bflag)
{

    // hardcoded for not to use adv derivs,
    deriv64_y(Dyu,u,dy,sz,bflag);

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv64adv_z(double * const  Dzu, const double * const  u,
                  const double dz, const unsigned int *sz,
                  const double * const betaz, unsigned bflag)
{

    // hardcoded for not to use adv derivs,
    deriv64_z(Dzu,u,dz,sz,bflag);
}
