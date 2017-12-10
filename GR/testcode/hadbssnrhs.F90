!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
      subroutine cal_had_bssn_rhs(
     &  dt_Alpha3D, dt_shiftx, dt_shifty, dt_shiftz, dt_chi3D,
     &  dt_trK3D, dt_gtxx, dt_gtxy, dt_gtxz, dt_gtyy, dt_gtyz,
     &  dt_gtzz, dt_Atxx, dt_Atxy, dt_Atxz, dt_Atyy, dt_Atyz,
     &  dt_Atzz, dt_Gamtx, dt_Gamty, dt_Gamtz, dt_gbx, dt_gby, dt_gbz,
     &  Alpha3D, shiftx, shifty, shiftz, chi3D, trK3D, 
     &  gtxx, gtxy, gtxz, gtyy, gtyz, gtzz, 
     &  Atxx, Atxy, Atxz, Atyy, Atyz,
     &  Atzz, Gamtx, Gamty, Gamtz, gbx, gby, gbz,
     &  dx_alpha,dy_alpha,dz_alpha,dx_shiftx,
     &  dy_shiftx,dz_shiftx,dx_shifty,dy_shifty,
     &  dz_shifty,dx_shiftz,dy_shiftz,dz_shiftz,
     &  dx_gbx,dy_gbx,dz_gbx,dx_gby,
     &  dy_gby,dz_gby,dx_gbz,dy_gbz,
     &  dz_gbz,dx_chi,dy_chi,dz_chi,
     &  dx_Gamtx,dy_Gamtx,dz_Gamtx,dx_Gamty,
     &  dy_Gamty,dz_Gamty,dx_Gamtz,dy_Gamtz,
     &  dz_Gamtz,dx_trK,dy_trK,dz_trK,
     &  dx_gtxx,dy_gtxx,dz_gtxx,dx_gtxy,
     &  dy_gtxy,dz_gtxy,dx_gtxz,dy_gtxz,
     &  dz_gtxz,dx_gtyy,dy_gtyy,dz_gtyy,
     &  dx_gtyz,dy_gtyz,dz_gtyz,dx_gtzz,
     &  dy_gtzz,dz_gtzz,dx_Atxx,dy_Atxx,
     &  dz_Atxx,dx_Atxy,dy_Atxy,dz_Atxy,
     &  dx_Atxz,dy_Atxz,dz_Atxz,dx_Atyy,
     &  dy_Atyy,dz_Atyy,dx_Atyz,dy_Atyz,
     &  dz_Atyz,dx_Atzz,dy_Atzz,dz_Atzz,
     &  dxx_gtxx,dxy_gtxx,dxz_gtxx,dyy_gtxx,
     &  dyz_gtxx,dzz_gtxx,dxx_gtxy,dxy_gtxy,
     &  dxz_gtxy,dyy_gtxy,dyz_gtxy,dzz_gtxy,
     &  dxx_gtxz,dxy_gtxz,dxz_gtxz,dyy_gtxz,
     &  dyz_gtxz,dzz_gtxz,dxx_gtyy,dxy_gtyy,
     &  dxz_gtyy,dyy_gtyy,dyz_gtyy,dzz_gtyy,
     &  dxx_gtyz,dxy_gtyz,dxz_gtyz,dyy_gtyz,
     &  dyz_gtyz,dzz_gtyz,dxx_gtzz,dxy_gtzz,
     &  dxz_gtzz,dyy_gtzz,dyz_gtzz,dzz_gtzz,
     &  dxx_chi,dxy_chi,dxz_chi,dyy_chi,
     &  dyz_chi,dzz_chi,dxx_alpha,dxy_alpha,
     &  dxz_alpha,dyy_alpha,dyz_alpha,dzz_alpha,
     &  dxx_shiftx,dxy_shiftx,dxz_shiftx,dyy_shiftx,
     &  dyz_shiftx,dzz_shiftx,dxx_shifty,dxy_shifty,
     &  dxz_shifty,dyy_shifty,dyz_shifty,dzz_shifty,
     &  dxx_shiftz,dxy_shiftz,dxz_shiftz,dyy_shiftz,
     &  dyz_shiftz,dzz_shiftz,
     &  lambda, lambda_f,
     &  eta, eta_damping, eta_damping_exp, R_0,
     &  trk0, chi_floor,
     &  x1d, y1d, z1d, nx, ny, nz)
        implicit none

        integer                                          :: nx, ny, nz
        real(kind=8), dimension(nx)                          :: x1d
        real(kind=8), dimension(ny)                          :: y1d
        real(kind=8), dimension(nz)                          :: z1d

      ! declare derivatives {{{
      real(kind=8), dimension(nx,ny,nz) :: dx_alpha
      real(kind=8), dimension(nx,ny,nz) :: dy_alpha
      real(kind=8), dimension(nx,ny,nz) :: dz_alpha
      real(kind=8), dimension(nx,ny,nz) :: dx_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dy_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dz_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dx_shifty
      real(kind=8), dimension(nx,ny,nz) :: dy_shifty
      real(kind=8), dimension(nx,ny,nz) :: dz_shifty
      real(kind=8), dimension(nx,ny,nz) :: dx_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dy_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dz_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dx_gbx
      real(kind=8), dimension(nx,ny,nz) :: dy_gbx
      real(kind=8), dimension(nx,ny,nz) :: dz_gbx
      real(kind=8), dimension(nx,ny,nz) :: dx_gby
      real(kind=8), dimension(nx,ny,nz) :: dy_gby
      real(kind=8), dimension(nx,ny,nz) :: dz_gby
      real(kind=8), dimension(nx,ny,nz) :: dx_gbz
      real(kind=8), dimension(nx,ny,nz) :: dy_gbz
      real(kind=8), dimension(nx,ny,nz) :: dz_gbz
      real(kind=8), dimension(nx,ny,nz) :: dx_chi
      real(kind=8), dimension(nx,ny,nz) :: dy_chi
      real(kind=8), dimension(nx,ny,nz) :: dz_chi
      real(kind=8), dimension(nx,ny,nz) :: dx_Gamtx
      real(kind=8), dimension(nx,ny,nz) :: dy_Gamtx
      real(kind=8), dimension(nx,ny,nz) :: dz_Gamtx
      real(kind=8), dimension(nx,ny,nz) :: dx_Gamty
      real(kind=8), dimension(nx,ny,nz) :: dy_Gamty
      real(kind=8), dimension(nx,ny,nz) :: dz_Gamty
      real(kind=8), dimension(nx,ny,nz) :: dx_Gamtz
      real(kind=8), dimension(nx,ny,nz) :: dy_Gamtz
      real(kind=8), dimension(nx,ny,nz) :: dz_Gamtz
      real(kind=8), dimension(nx,ny,nz) :: dx_trK
      real(kind=8), dimension(nx,ny,nz) :: dy_trK
      real(kind=8), dimension(nx,ny,nz) :: dz_trK
      real(kind=8), dimension(nx,ny,nz) :: dx_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dy_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dz_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dx_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dy_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dz_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dx_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dy_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dz_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dx_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dy_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dz_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dx_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dy_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dz_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dx_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dy_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dz_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dx_Atxx
      real(kind=8), dimension(nx,ny,nz) :: dy_Atxx
      real(kind=8), dimension(nx,ny,nz) :: dz_Atxx
      real(kind=8), dimension(nx,ny,nz) :: dx_Atxy
      real(kind=8), dimension(nx,ny,nz) :: dy_Atxy
      real(kind=8), dimension(nx,ny,nz) :: dz_Atxy
      real(kind=8), dimension(nx,ny,nz) :: dx_Atxz
      real(kind=8), dimension(nx,ny,nz) :: dy_Atxz
      real(kind=8), dimension(nx,ny,nz) :: dz_Atxz
      real(kind=8), dimension(nx,ny,nz) :: dx_Atyy
      real(kind=8), dimension(nx,ny,nz) :: dy_Atyy
      real(kind=8), dimension(nx,ny,nz) :: dz_Atyy
      real(kind=8), dimension(nx,ny,nz) :: dx_Atyz
      real(kind=8), dimension(nx,ny,nz) :: dy_Atyz
      real(kind=8), dimension(nx,ny,nz) :: dz_Atyz
      real(kind=8), dimension(nx,ny,nz) :: dx_Atzz
      real(kind=8), dimension(nx,ny,nz) :: dy_Atzz
      real(kind=8), dimension(nx,ny,nz) :: dz_Atzz
      real(kind=8), dimension(nx,ny,nz) :: dxx_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dxy_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dxz_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dyy_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dyz_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dzz_gtxx
      real(kind=8), dimension(nx,ny,nz) :: dxx_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dxy_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dxz_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dyy_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dyz_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dzz_gtxy
      real(kind=8), dimension(nx,ny,nz) :: dxx_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dxy_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dxz_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dyy_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dyz_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dzz_gtxz
      real(kind=8), dimension(nx,ny,nz) :: dxx_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dxy_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dxz_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dyy_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dyz_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dzz_gtyy
      real(kind=8), dimension(nx,ny,nz) :: dxx_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dxy_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dxz_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dyy_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dyz_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dzz_gtyz
      real(kind=8), dimension(nx,ny,nz) :: dxx_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dxy_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dxz_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dyy_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dyz_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dzz_gtzz
      real(kind=8), dimension(nx,ny,nz) :: dxx_chi
      real(kind=8), dimension(nx,ny,nz) :: dxy_chi
      real(kind=8), dimension(nx,ny,nz) :: dxz_chi
      real(kind=8), dimension(nx,ny,nz) :: dyy_chi
      real(kind=8), dimension(nx,ny,nz) :: dyz_chi
      real(kind=8), dimension(nx,ny,nz) :: dzz_chi
      real(kind=8), dimension(nx,ny,nz) :: dxx_alpha
      real(kind=8), dimension(nx,ny,nz) :: dxy_alpha
      real(kind=8), dimension(nx,ny,nz) :: dxz_alpha
      real(kind=8), dimension(nx,ny,nz) :: dyy_alpha
      real(kind=8), dimension(nx,ny,nz) :: dyz_alpha
      real(kind=8), dimension(nx,ny,nz) :: dzz_alpha
      real(kind=8), dimension(nx,ny,nz) :: dxx_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dxy_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dxz_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dyy_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dyz_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dzz_shiftx
      real(kind=8), dimension(nx,ny,nz) :: dxx_shifty
      real(kind=8), dimension(nx,ny,nz) :: dxy_shifty
      real(kind=8), dimension(nx,ny,nz) :: dxz_shifty
      real(kind=8), dimension(nx,ny,nz) :: dyy_shifty
      real(kind=8), dimension(nx,ny,nz) :: dyz_shifty
      real(kind=8), dimension(nx,ny,nz) :: dzz_shifty
      real(kind=8), dimension(nx,ny,nz) :: dxx_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dxy_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dxz_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dyy_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dyz_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dzz_shiftz
      real(kind=8), dimension(nx,ny,nz) :: dx_Ex
      real(kind=8), dimension(nx,ny,nz) :: dy_Ex
      real(kind=8), dimension(nx,ny,nz) :: dz_Ex
      real(kind=8), dimension(nx,ny,nz) :: dx_Ey
      real(kind=8), dimension(nx,ny,nz) :: dy_Ey
      real(kind=8), dimension(nx,ny,nz) :: dz_Ey
      real(kind=8), dimension(nx,ny,nz) :: dx_Ez
      real(kind=8), dimension(nx,ny,nz) :: dy_Ez
      real(kind=8), dimension(nx,ny,nz) :: dz_Ez
      real(kind=8), dimension(nx,ny,nz) :: dx_Bx
      real(kind=8), dimension(nx,ny,nz) :: dy_Bx
      real(kind=8), dimension(nx,ny,nz) :: dz_Bx
      real(kind=8), dimension(nx,ny,nz) :: dx_By
      real(kind=8), dimension(nx,ny,nz) :: dy_By
      real(kind=8), dimension(nx,ny,nz) :: dz_By
      real(kind=8), dimension(nx,ny,nz) :: dx_Bz
      real(kind=8), dimension(nx,ny,nz) :: dy_Bz
      real(kind=8), dimension(nx,ny,nz) :: dz_Bz
      real(kind=8), dimension(nx,ny,nz) :: dx_Phi_em
      real(kind=8), dimension(nx,ny,nz) :: dy_Phi_em
      real(kind=8), dimension(nx,ny,nz) :: dz_Phi_em
      real(kind=8), dimension(nx,ny,nz) :: dx_Psi_em
      real(kind=8), dimension(nx,ny,nz) :: dy_Psi_em
      real(kind=8), dimension(nx,ny,nz) :: dz_Psi_em
      ! }}}

        ! declare rhs and vars {{{
        real(kind=8), dimension(nx,ny,nz) :: Alpha3D
        real(kind=8), dimension(nx,ny,nz) :: shiftx
        real(kind=8), dimension(nx,ny,nz) :: shifty
        real(kind=8), dimension(nx,ny,nz) :: shiftz
        real(kind=8), dimension(nx,ny,nz) :: gtxx
        real(kind=8), dimension(nx,ny,nz) :: gtxy
        real(kind=8), dimension(nx,ny,nz) :: gtxz
        real(kind=8), dimension(nx,ny,nz) :: gtyy
        real(kind=8), dimension(nx,ny,nz) :: gtyz
        real(kind=8), dimension(nx,ny,nz) :: gtzz
        real(kind=8), dimension(nx,ny,nz) :: Atxx 
        real(kind=8), dimension(nx,ny,nz) :: Atxy 
        real(kind=8), dimension(nx,ny,nz) :: Atxz
        real(kind=8), dimension(nx,ny,nz) :: Atyy 
        real(kind=8), dimension(nx,ny,nz) :: Atyz 
        real(kind=8), dimension(nx,ny,nz) :: Atzz
        real(kind=8), dimension(nx,ny,nz) :: chi3D
        real(kind=8), dimension(nx,ny,nz) :: trK3D
        real(kind=8), dimension(nx,ny,nz) :: Gamtx
        real(kind=8), dimension(nx,ny,nz) :: Gamty
        real(kind=8), dimension(nx,ny,nz) :: Gamtz
        real(kind=8), dimension(nx,ny,nz) :: gxx 
        real(kind=8), dimension(nx,ny,nz) :: gxy 
        real(kind=8), dimension(nx,ny,nz) :: gxz
        real(kind=8), dimension(nx,ny,nz) :: gyy 
        real(kind=8), dimension(nx,ny,nz) :: gyz 
        real(kind=8), dimension(nx,ny,nz) :: gzz
        real(kind=8), dimension(nx,ny,nz) :: gbx 
        real(kind=8), dimension(nx,ny,nz) :: gby 
        real(kind=8), dimension(nx,ny,nz) :: gbz

        real(kind=8), dimension(nx,ny,nz) :: dt_Alpha3D
        real(kind=8), dimension(nx,ny,nz) :: dt_shiftx
        real(kind=8), dimension(nx,ny,nz) :: dt_shifty 
        real(kind=8), dimension(nx,ny,nz) :: dt_shiftz
        real(kind=8), dimension(nx,ny,nz) :: dt_gtxx
        real(kind=8), dimension(nx,ny,nz) :: dt_gtxy
        real(kind=8), dimension(nx,ny,nz) :: dt_gtxz
        real(kind=8), dimension(nx,ny,nz) :: dt_gtyy
        real(kind=8), dimension(nx,ny,nz) :: dt_gtyz
        real(kind=8), dimension(nx,ny,nz) :: dt_gtzz
        real(kind=8), dimension(nx,ny,nz) :: dt_Atxx
        real(kind=8), dimension(nx,ny,nz) :: dt_Atxy
        real(kind=8), dimension(nx,ny,nz) :: dt_Atxz
        real(kind=8), dimension(nx,ny,nz) :: dt_Atyy
        real(kind=8), dimension(nx,ny,nz) :: dt_Atyz
        real(kind=8), dimension(nx,ny,nz) :: dt_Atzz
        real(kind=8), dimension(nx,ny,nz) :: dt_chi3D
        real(kind=8), dimension(nx,ny,nz) :: dt_trK3D
        real(kind=8), dimension(nx,ny,nz) :: dt_Gamtx
        real(kind=8), dimension(nx,ny,nz) :: dt_Gamty
        real(kind=8), dimension(nx,ny,nz) :: dt_Gamtz
        real(kind=8), dimension(nx,ny,nz) :: dt_gbx
        real(kind=8), dimension(nx,ny,nz) :: dt_gby
        real(kind=8), dimension(nx,ny,nz) :: dt_gbz
        ! }}}

        real(kind=8)            :: chi_floor, trK0, R_0
        real(kind=8)            :: eta
        integer                 :: eta_damping, eta_damping_exp
        integer, dimension(4)   :: lambda
        integer, dimension(2)   :: lambda_f


        integer                       :: i, j, k, ii
        integer                       :: rc, err

        real(kind=8)                      :: dx, dy, dz
        real(kind=8)                      :: trK, Alpha, chi, G
        real(kind=8)                      :: rho_ADM, trPsi1, tr_pT
        real(kind=8)                      :: Alpha_rhs, trK_rhs, chi_rhs
        real(kind=8),  dimension(3)       :: Betau, Bu_rhs, Betau_rhs
        real(kind=8),  dimension(3)       :: Jtd_ADM, CalGamt, Gamt 
        real(kind=8),  dimension(3)       :: d_Alpha
        real(kind=8),  dimension(3)       :: d_chi, d_trK
        real(kind=8),  dimension(3)       :: Gamt_rhs
        real(kind=8),  dimension(3)       :: Bu
        real(kind=8),  dimension(3,3)     :: gd, gu, gtd, gtu
        real(kind=8),  dimension(3,3)     :: Atd, Atu, Atud, dd_chi
        real(kind=8),  dimension(3,3)     :: Rpd, Rtd, Rpd_1
        real(kind=8),  dimension(3,3)     :: pTtd_ADM
        real(kind=8),  dimension(3,3)     :: dd_Alpha, d_Gamt 
        real(kind=8),  dimension(3,3)     :: d_Betau, Psi1, Psi1TF
        real(kind=8),  dimension(3,3)     :: Atd_rhs, gtd_rhs
        real(kind=8),  dimension(3,3,3)   :: dd_Betau, d_gtd, d_Atd
        real(kind=8),  dimension(3,3,3,3) :: dd_gtd 
        real(kind=8),  dimension(0:3,0:3) :: gd4, gu4, Tu
        real(kind=8)                      :: detgd, idetgd
        real(kind=8)                      :: detgtd, idetgtd
        real(kind=8),  dimension(3,3,3)   :: Ctd, Ct, Chr
        real(kind=8)                      :: twelve_Pi_G
        real(kind=8)                      :: third_trPsi1
        real(kind=8)                      :: inv_chi, sqrt_chi 

        real(kind=8)                      :: trK2_rhs, trK3_rhs
        real(kind=8) :: t32

        real(kind=8)                      :: div_Beta  
        real(kind=8), dimension(3)        :: d_div_Beta
 
        real(kind=8), dimension(3)        :: xpt 
        real(kind=8)                      :: r_coord
 
        real(kind=8)                      :: det_gtd, idet_gtd
        real(kind=8)                      :: det_gtd_to_third
        real(kind=8)                      :: det_gtd_to_neg_third
        real(kind=8)                      :: trace_Atd
        real(kind=8)                      :: neg_one_third_trace_Atd

        integer                       :: lambda_1,  lambda_2
        integer                       :: lambda_3,  lambda_4
        integer                       :: lambda_f0, lambda_f1

        real(kind=8)                      :: feta
  
        ! "Fluff" variable
        !    ala Giddings
        !    or other weird stuff ala
        !    http://arxiv.org/pdf/1608.05974v1.pdf
        real(kind=8)                      :: tmp, radsq, radsq1, radsq2
        ! Center of the first black hole:
        real(kind=8)                      :: xcen1, ycen1, zcen1
        ! Center of the second black hole:
        real(kind=8)                      :: xcen2, ycen2, zcen2
        real(kind=8)                      :: time

        integer                       :: force_free 

        integer                       :: evolve_geometry

        real(kind=8), parameter    :: twelve = 12.0d0
        real(kind=8), parameter    :: six    = 6.0d0
        real(kind=8), parameter    :: four   = 4.0d0
        real(kind=8), parameter    :: three  = 3.0d0
        real(kind=8), parameter    :: two    = 2.0d0
        real(kind=8), parameter    :: one    = 1.0d0
        real(kind=8), parameter    :: threehalves  = 1.5d0
        real(kind=8), parameter    :: threefourths = 0.75d0
        real(kind=8), parameter    :: twothirds = 0.6666666666666667d0
        real(kind=8), parameter    :: half   = 0.5d0
        real(kind=8), parameter    :: third  = 0.3333333333333333d0
        real(kind=8), parameter    :: fourth = 0.25d0 
        real(kind=8), parameter    :: sixth  = 0.1666666666666667d0
        real(kind=8), parameter    :: four_pi_G  = 12.5663706143591729d0 
        real(kind=8), parameter    :: eight_pi_G = 25.1327412287183459d0 

        character(len=*), parameter    :: FMT3='(3E12.3)'
        character(len=*), parameter    :: FMT6='(6E12.3)'

        logical, parameter             :: ltrace      = .false.
        logical, parameter             :: ltrace2     = .false.

        ! include declarations for Maple temp variables
#include "bssn_emtest_temp_vars.h"

        ! G is Newton"s gravitation constant
        G = 1.0d0
        lambda_1  = lambda(1)
        lambda_2  = lambda(2)
        lambda_3  = lambda(3)
        lambda_4  = lambda(4)

        lambda_f0 = lambda_f(1)
        lambda_f1 = lambda_f(2)

        dx = x1d(2) - x1d(1)
        dy = y1d(2) - y1d(1)
        dz = z1d(2) - z1d(1)

        evolve_geometry = 1

        rho_ADM = 0.0d0
        tr_pT   = 0.0d0
        Jtd_ADM(1)   = 0.0d0
        Jtd_ADM(2)   = 0.0d0
        Jtd_ADM(3)   = 0.0d0
        pTtd_ADM(1,1)   = 0.0d0
        pTtd_ADM(2,1)   = 0.0d0
        pTtd_ADM(3,1)   = 0.0d0
        pTtd_ADM(1,2)   = 0.0d0
        pTtd_ADM(2,2)   = 0.0d0
        pTtd_ADM(3,2)   = 0.0d0
        pTtd_ADM(1,3)   = 0.0d0
        pTtd_ADM(2,3)   = 0.0d0
        pTtd_ADM(3,3)   = 0.0d0


        do k = 1, nz
        do j = 1, ny
        do i = 1, nx

          xpt(1) = x1d(i)
          xpt(2) = y1d(j)
          xpt(3) = z1d(k)

        ! Construct a radial coordinate from the Cartesian coordinates 
        ! in the usual way.  We need this if we damp the eta factor in 
        ! the gauge conditions.  (See arxiv:1003.0859) 
          r_coord = sqrt(    xpt(1)*xpt(1)
     &                     + xpt(2)*xpt(2)
     &                     + xpt(3)*xpt(3)  )

! define somehow the value of feta, the damped eta
	  feta = eta
	  if(r_coord.ge.R_0) feta=eta*(R_0/r_coord)**eta_damping_exp

          gd(1,1) = gxx(i,j,k)
          gd(1,2) = gxy(i,j,k)
          gd(1,3) = gxz(i,j,k)
          gd(2,1) = gxy(i,j,k)
          gd(2,2) = gyy(i,j,k)
          gd(2,3) = gyz(i,j,k)
          gd(3,1) = gxz(i,j,k)
          gd(3,2) = gyz(i,j,k)
          gd(3,3) = gzz(i,j,k)

          gtd(1,1) = gtxx(i,j,k)
          gtd(1,2) = gtxy(i,j,k)
          gtd(1,3) = gtxz(i,j,k)
          gtd(2,1) = gtxy(i,j,k)
          gtd(2,2) = gtyy(i,j,k)
          gtd(2,3) = gtyz(i,j,k)
          gtd(3,1) = gtxz(i,j,k)
          gtd(3,2) = gtyz(i,j,k)
          gtd(3,3) = gtzz(i,j,k)

          Atd(1,1) = Atxx(i,j,k)
          Atd(1,2) = Atxy(i,j,k)
          Atd(1,3) = Atxz(i,j,k)
          Atd(2,1) = Atxy(i,j,k)
          Atd(2,2) = Atyy(i,j,k)
          Atd(2,3) = Atyz(i,j,k)
          Atd(3,1) = Atxz(i,j,k)
          Atd(3,2) = Atyz(i,j,k)
          Atd(3,3) = Atzz(i,j,k)

          trK      = trK3D(i,j,k)
          chi      = chi3D(i,j,k)

! enforce a floor on chi 
          chi = max( chi3D(i,j,k) , chi_floor )

! and take the square root of chi (this is needed in the E&M equations)  
          if ( ltrace ) then
            write(0,*)' cal_bssn_rhs:  Checking that chi is nonnegative'
            if ( chi .lt. 0.d0 ) then
               write(0,*)' chi = ', chi
               write(0,*)' and the square root of chi will be undefined'
            endif
          endif

          sqrt_chi = sqrt( chi )

          Gamt(1) = Gamtx(i,j,k)
          Gamt(2) = Gamty(i,j,k)
          Gamt(3) = Gamtz(i,j,k)

          Alpha    = Alpha3D(i,j,k)
          Betau(1) = shiftx(i,j,k)
          Betau(2) = shifty(i,j,k)
          Betau(3) = shiftz(i,j,k)

          d_Alpha(1)    = dx_alpha(i,j,k)
          d_Alpha(2)    = dy_alpha(i,j,k)
          d_Alpha(3)    = dz_alpha(i,j,k)

          dd_Alpha(1,1) = dxx_alpha(i,j,k)
          dd_Alpha(1,2) = dxy_alpha(i,j,k)
          dd_Alpha(1,3) = dxz_alpha(i,j,k)
          dd_Alpha(2,2) = dyy_alpha(i,j,k)
          dd_Alpha(2,3) = dyz_alpha(i,j,k)
          dd_Alpha(3,3) = dzz_alpha(i,j,k)

          d_chi(1)      = dx_chi(i,j,k)
          d_chi(2)      = dy_chi(i,j,k)
          d_chi(3)      = dz_chi(i,j,k)

          dd_chi(1,1)   = dxx_chi(i,j,k)
          dd_chi(1,2)   = dxy_chi(i,j,k)
          dd_chi(1,3)   = dxz_chi(i,j,k)
          dd_chi(2,2)   = dyy_chi(i,j,k)
          dd_chi(2,3)   = dyz_chi(i,j,k)
          dd_chi(3,3)   = dzz_chi(i,j,k)

          d_trK(1)        = dx_trK(i,j,k)
          d_trK(2)        = dy_trK(i,j,k)
          d_trK(3)        = dz_trK(i,j,k)

          d_Gamt(1,1)   = dx_Gamtx(i,j,k)
          d_Gamt(2,1)   = dy_Gamtx(i,j,k)
          d_Gamt(3,1)   = dz_Gamtx(i,j,k)
          d_Gamt(1,2)   = dx_Gamty(i,j,k)
          d_Gamt(2,2)   = dy_Gamty(i,j,k)
          d_Gamt(3,2)   = dz_Gamty(i,j,k)
          d_Gamt(1,3)   = dx_Gamtz(i,j,k)
          d_Gamt(2,3)   = dy_Gamtz(i,j,k)
          d_Gamt(3,3)   = dz_Gamtz(i,j,k)

          d_Betau(1,1)   = dx_shiftx(i,j,k)
          d_Betau(2,1)   = dy_shiftx(i,j,k)
          d_Betau(3,1)   = dz_shiftx(i,j,k)
          d_Betau(1,2)   = dx_shifty(i,j,k)
          d_Betau(2,2)   = dy_shifty(i,j,k)
          d_Betau(3,2)   = dz_shifty(i,j,k)
          d_Betau(1,3)   = dx_shiftz(i,j,k)
          d_Betau(2,3)   = dy_shiftz(i,j,k)
          d_Betau(3,3)   = dz_shiftz(i,j,k)

          dd_Betau(1,1,1)   = dxx_shiftx(i,j,k)
          dd_Betau(1,2,1)   = dxy_shiftx(i,j,k)
          dd_Betau(1,3,1)   = dxz_shiftx(i,j,k)
          dd_Betau(2,1,1)   = dd_Betau(1,2,1)
          dd_Betau(2,2,1)   = dyy_shiftx(i,j,k)
          dd_Betau(2,3,1)   = dyz_shiftx(i,j,k)
          dd_Betau(3,1,1)   = dd_Betau(1,3,1)
          dd_Betau(3,2,1)   = dd_Betau(2,3,1)
          dd_Betau(3,3,1)   = dzz_shiftx(i,j,k)

          dd_Betau(1,1,2)   = dxx_shifty(i,j,k)
          dd_Betau(1,2,2)   = dxy_shifty(i,j,k)
          dd_Betau(1,3,2)   = dxz_shifty(i,j,k)
          dd_Betau(2,1,2)   = dd_Betau(1,2,2)
          dd_Betau(2,2,2)   = dyy_shifty(i,j,k)
          dd_Betau(2,3,2)   = dyz_shifty(i,j,k)
          dd_Betau(3,1,2)   = dd_Betau(1,3,2)
          dd_Betau(3,2,2)   = dd_Betau(2,3,2)
          dd_Betau(3,3,2)   = dzz_shifty(i,j,k)

          dd_Betau(1,1,3)   = dxx_shiftz(i,j,k)
          dd_Betau(1,2,3)   = dxy_shiftz(i,j,k)
          dd_Betau(1,3,3)   = dxz_shiftz(i,j,k)
          dd_Betau(2,1,3)   = dd_Betau(1,2,3)
          dd_Betau(2,2,3)   = dyy_shiftz(i,j,k)
          dd_Betau(2,3,3)   = dyz_shiftz(i,j,k)
          dd_Betau(3,1,3)   = dd_Betau(1,3,3)
          dd_Betau(3,2,3)   = dd_Betau(2,3,3)
          dd_Betau(3,3,3)   = dzz_shiftz(i,j,k)

          d_gtd(1,1,1)      = dx_gtxx(i,j,k)
          d_gtd(2,1,1)      = dy_gtxx(i,j,k)
          d_gtd(3,1,1)      = dz_gtxx(i,j,k)

          d_gtd(1,1,2)      = dx_gtxy(i,j,k)
          d_gtd(2,1,2)      = dy_gtxy(i,j,k)
          d_gtd(3,1,2)      = dz_gtxy(i,j,k)

          d_gtd(1,1,3)      = dx_gtxz(i,j,k)
          d_gtd(2,1,3)      = dy_gtxz(i,j,k)
          d_gtd(3,1,3)      = dz_gtxz(i,j,k)

          d_gtd(1,2,1)      = d_gtd(1,1,2)
          d_gtd(2,2,1)      = d_gtd(2,1,2)
          d_gtd(3,2,1)      = d_gtd(3,1,2)

          d_gtd(1,2,2)      = dx_gtyy(i,j,k)
          d_gtd(2,2,2)      = dy_gtyy(i,j,k)
          d_gtd(3,2,2)      = dz_gtyy(i,j,k)

          d_gtd(1,2,3)      = dx_gtyz(i,j,k)
          d_gtd(2,2,3)      = dy_gtyz(i,j,k)
          d_gtd(3,2,3)      = dz_gtyz(i,j,k)

          d_gtd(1,3,1)      = d_gtd(1,1,3)
          d_gtd(2,3,1)      = d_gtd(2,1,3)
          d_gtd(3,3,1)      = d_gtd(3,1,3)

          d_gtd(1,3,2)      = d_gtd(1,2,3)
          d_gtd(2,3,2)      = d_gtd(2,2,3)
          d_gtd(3,3,2)      = d_gtd(3,2,3)

          d_gtd(1,3,3)      = dx_gtzz(i,j,k)
          d_gtd(2,3,3)      = dy_gtzz(i,j,k)
          d_gtd(3,3,3)      = dz_gtzz(i,j,k)


          d_Atd(1,1,1)      = dx_Atxx(i,j,k)
          d_Atd(2,1,1)      = dy_Atxx(i,j,k)
          d_Atd(3,1,1)      = dz_Atxx(i,j,k)

          d_Atd(1,1,2)      = dx_Atxy(i,j,k)
          d_Atd(2,1,2)      = dy_Atxy(i,j,k)
          d_Atd(3,1,2)      = dz_Atxy(i,j,k)

          d_Atd(1,1,3)      = dx_Atxz(i,j,k)
          d_Atd(2,1,3)      = dy_Atxz(i,j,k)
          d_Atd(3,1,3)      = dz_Atxz(i,j,k)

          d_Atd(1,2,1)      = d_Atd(1,1,2)
          d_Atd(2,2,1)      = d_Atd(2,1,2)
          d_Atd(3,2,1)      = d_Atd(3,1,2)

          d_Atd(1,2,2)      = dx_Atyy(i,j,k)
          d_Atd(2,2,2)      = dy_Atyy(i,j,k)
          d_Atd(3,2,2)      = dz_Atyy(i,j,k)

          d_Atd(1,2,3)      = dx_Atyz(i,j,k)
          d_Atd(2,2,3)      = dy_Atyz(i,j,k)
          d_Atd(3,2,3)      = dz_Atyz(i,j,k)

          d_Atd(1,3,1)      = d_Atd(1,1,3)
          d_Atd(2,3,1)      = d_Atd(2,1,3)
          d_Atd(3,3,1)      = d_Atd(3,1,3)

          d_Atd(1,3,2)      = d_Atd(1,2,3)
          d_Atd(2,3,2)      = d_Atd(2,2,3)
          d_Atd(3,3,2)      = d_Atd(3,2,3)

          d_Atd(1,3,3)      = dx_Atzz(i,j,k)
          d_Atd(2,3,3)      = dy_Atzz(i,j,k)
          d_Atd(3,3,3)      = dz_Atzz(i,j,k)


          dd_gtd(1,1,1,1)      = dxx_gtxx(i,j,k)
          dd_gtd(1,2,1,1)      = dxy_gtxx(i,j,k)
          dd_gtd(1,3,1,1)      = dxz_gtxx(i,j,k)
          dd_gtd(2,2,1,1)      = dyy_gtxx(i,j,k)
          dd_gtd(2,3,1,1)      = dyz_gtxx(i,j,k)
          dd_gtd(3,3,1,1)      = dzz_gtxx(i,j,k)

          dd_gtd(1,1,1,2)      = dxx_gtxy(i,j,k)
          dd_gtd(1,2,1,2)      = dxy_gtxy(i,j,k)
          dd_gtd(1,3,1,2)      = dxz_gtxy(i,j,k)
          dd_gtd(2,2,1,2)      = dyy_gtxy(i,j,k)
          dd_gtd(2,3,1,2)      = dyz_gtxy(i,j,k)
          dd_gtd(3,3,1,2)      = dzz_gtxy(i,j,k)

          dd_gtd(1,1,1,3)      = dxx_gtxz(i,j,k)
          dd_gtd(1,2,1,3)      = dxy_gtxz(i,j,k)
          dd_gtd(1,3,1,3)      = dxz_gtxz(i,j,k)
          dd_gtd(2,2,1,3)      = dyy_gtxz(i,j,k)
          dd_gtd(2,3,1,3)      = dyz_gtxz(i,j,k)
          dd_gtd(3,3,1,3)      = dzz_gtxz(i,j,k)

          dd_gtd(1,1,2,2)      = dxx_gtyy(i,j,k)
          dd_gtd(1,2,2,2)      = dxy_gtyy(i,j,k)
          dd_gtd(1,3,2,2)      = dxz_gtyy(i,j,k)
          dd_gtd(2,2,2,2)      = dyy_gtyy(i,j,k)
          dd_gtd(2,3,2,2)      = dyz_gtyy(i,j,k)
          dd_gtd(3,3,2,2)      = dzz_gtyy(i,j,k)

          dd_gtd(1,1,2,3)      = dxx_gtyz(i,j,k)
          dd_gtd(1,2,2,3)      = dxy_gtyz(i,j,k)
          dd_gtd(1,3,2,3)      = dxz_gtyz(i,j,k)
          dd_gtd(2,2,2,3)      = dyy_gtyz(i,j,k)
          dd_gtd(2,3,2,3)      = dyz_gtyz(i,j,k)
          dd_gtd(3,3,2,3)      = dzz_gtyz(i,j,k)

          dd_gtd(1,1,3,3)      = dxx_gtzz(i,j,k)
          dd_gtd(1,2,3,3)      = dxy_gtzz(i,j,k)
          dd_gtd(1,3,3,3)      = dxz_gtzz(i,j,k)
          dd_gtd(2,2,3,3)      = dyy_gtzz(i,j,k)
          dd_gtd(2,3,3,3)      = dyz_gtzz(i,j,k)
          dd_gtd(3,3,3,3)      = dzz_gtzz(i,j,k)

          Bu(1)                = gbx(i,j,k)
          Bu(2)                = gby(i,j,k)
          Bu(3)                = gbz(i,j,k)
 
#include "bssn_emtest.h" 

            dt_Alpha3D(i,j,k)   = Alpha_rhs
            dt_shiftx(i,j,k)    = Betau_rhs(1)
            dt_shifty(i,j,k)    = Betau_rhs(2)
            dt_shiftz(i,j,k)    = Betau_rhs(3)

            dt_chi3D(i,j,k)     = chi_rhs
            dt_trK3D(i,j,k)     = trK_rhs

            dt_gtxx(i,j,k)      = gtd_rhs(1,1)
            dt_gtxy(i,j,k)      = gtd_rhs(1,2)
            dt_gtxz(i,j,k)      = gtd_rhs(1,3)
            dt_gtyy(i,j,k)      = gtd_rhs(2,2)
            dt_gtyz(i,j,k)      = gtd_rhs(2,3)
            dt_gtzz(i,j,k)      = gtd_rhs(3,3)

            dt_Atxx(i,j,k)      = Atd_rhs(1,1)
            dt_Atxy(i,j,k)      = Atd_rhs(1,2)
            dt_Atxz(i,j,k)      = Atd_rhs(1,3)
            dt_Atyy(i,j,k)      = Atd_rhs(2,2)
            dt_Atyz(i,j,k)      = Atd_rhs(2,3)
            dt_Atzz(i,j,k)      = Atd_rhs(3,3)

            dt_Gamtx(i,j,k)     = Gamt_rhs(1)
            dt_Gamty(i,j,k)     = Gamt_rhs(2)
            dt_Gamtz(i,j,k)     = Gamt_rhs(3)
  
            dt_gbx(i,j,k)       = Bu_rhs(1)
            dt_gby(i,j,k)       = Bu_rhs(2)
            dt_gbz(i,j,k)       = Bu_rhs(3)

        end do
        end do
        end do



        return
      end  subroutine

!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
      subroutine cal_had_bssn_rhs_lie(
     &  dt_Alpha3D, dt_shiftx, dt_shifty, dt_shiftz, dt_chi3D,
     &  dt_trK3D, dt_gtxx, dt_gtxy, dt_gtxz, dt_gtyy, dt_gtyz,
     &  dt_gtzz, dt_Atxx, dt_Atxy, dt_Atxz, dt_Atyy, dt_Atyz,
     &  dt_Atzz, dt_Gamtx, dt_Gamty, dt_Gamtz, dt_gbx, dt_gby, dt_gbz,
     &  shiftx, shifty, shiftz,
     &  adx_gtxx,ady_gtxx,adz_gtxx,adx_gtxy,
     &  ady_gtxy,adz_gtxy,adx_gtxz,ady_gtxz,
     &  adz_gtxz,adx_gtyy,ady_gtyy,adz_gtyy,
     &  adx_gtyz,ady_gtyz,adz_gtyz,adx_gtzz,
     &  ady_gtzz,adz_gtzz,adx_Atxx,ady_Atxx,
     &  adz_Atxx,adx_Atxy,ady_Atxy,adz_Atxy,
     &  adx_Atxz,ady_Atxz,adz_Atxz,adx_Atyy,
     &  ady_Atyy,adz_Atyy,adx_Atyz,ady_Atyz,
     &  adz_Atyz,adx_Atzz,ady_Atzz,adz_Atzz,
     &  adx_alpha,ady_alpha,adz_alpha,adx_shiftx,
     &  ady_shiftx,adz_shiftx,adx_shifty,ady_shifty,
     &  adz_shifty,adx_shiftz,ady_shiftz,adz_shiftz,
     &  adx_chi,ady_chi,adz_chi,adx_Gamtx,
     &  ady_Gamtx,adz_Gamtx,adx_Gamty,ady_Gamty,
     &  adz_Gamty,adx_Gamtz,ady_Gamtz,adz_Gamtz,
     &  adx_trK,ady_trK,adz_trK,adx_gbx,
     &  ady_gbx,adz_gbx,adx_gby,ady_gby,
     &  adz_gby,adx_gbz,ady_gbz,adz_gbz,
     &  lambda, lambda_f,
     &  x1d, y1d, z1d, nx, ny, nz)
        implicit none

        integer   :: nx, ny, nz
        real(kind=8), dimension(nx)                          :: x1d
        real(kind=8), dimension(ny)                          :: y1d
        real(kind=8), dimension(nz)                          :: z1d

      ! declare derivatives {{{
      real(kind=8), dimension(nx,ny,nz) :: adx_gtxx
      real(kind=8), dimension(nx,ny,nz) :: ady_gtxx
      real(kind=8), dimension(nx,ny,nz) :: adz_gtxx
      real(kind=8), dimension(nx,ny,nz) :: adx_gtxy
      real(kind=8), dimension(nx,ny,nz) :: ady_gtxy
      real(kind=8), dimension(nx,ny,nz) :: adz_gtxy
      real(kind=8), dimension(nx,ny,nz) :: adx_gtxz
      real(kind=8), dimension(nx,ny,nz) :: ady_gtxz
      real(kind=8), dimension(nx,ny,nz) :: adz_gtxz
      real(kind=8), dimension(nx,ny,nz) :: adx_gtyy
      real(kind=8), dimension(nx,ny,nz) :: ady_gtyy
      real(kind=8), dimension(nx,ny,nz) :: adz_gtyy
      real(kind=8), dimension(nx,ny,nz) :: adx_gtyz
      real(kind=8), dimension(nx,ny,nz) :: ady_gtyz
      real(kind=8), dimension(nx,ny,nz) :: adz_gtyz
      real(kind=8), dimension(nx,ny,nz) :: adx_gtzz
      real(kind=8), dimension(nx,ny,nz) :: ady_gtzz
      real(kind=8), dimension(nx,ny,nz) :: adz_gtzz
      real(kind=8), dimension(nx,ny,nz) :: adx_Atxx
      real(kind=8), dimension(nx,ny,nz) :: ady_Atxx
      real(kind=8), dimension(nx,ny,nz) :: adz_Atxx
      real(kind=8), dimension(nx,ny,nz) :: adx_Atxy
      real(kind=8), dimension(nx,ny,nz) :: ady_Atxy
      real(kind=8), dimension(nx,ny,nz) :: adz_Atxy
      real(kind=8), dimension(nx,ny,nz) :: adx_Atxz
      real(kind=8), dimension(nx,ny,nz) :: ady_Atxz
      real(kind=8), dimension(nx,ny,nz) :: adz_Atxz
      real(kind=8), dimension(nx,ny,nz) :: adx_Atyy
      real(kind=8), dimension(nx,ny,nz) :: ady_Atyy
      real(kind=8), dimension(nx,ny,nz) :: adz_Atyy
      real(kind=8), dimension(nx,ny,nz) :: adx_Atyz
      real(kind=8), dimension(nx,ny,nz) :: ady_Atyz
      real(kind=8), dimension(nx,ny,nz) :: adz_Atyz
      real(kind=8), dimension(nx,ny,nz) :: adx_Atzz
      real(kind=8), dimension(nx,ny,nz) :: ady_Atzz
      real(kind=8), dimension(nx,ny,nz) :: adz_Atzz
      real(kind=8), dimension(nx,ny,nz) :: adx_alpha
      real(kind=8), dimension(nx,ny,nz) :: ady_alpha
      real(kind=8), dimension(nx,ny,nz) :: adz_alpha
      real(kind=8), dimension(nx,ny,nz) :: adx_shiftx
      real(kind=8), dimension(nx,ny,nz) :: ady_shiftx
      real(kind=8), dimension(nx,ny,nz) :: adz_shiftx
      real(kind=8), dimension(nx,ny,nz) :: adx_shifty
      real(kind=8), dimension(nx,ny,nz) :: ady_shifty
      real(kind=8), dimension(nx,ny,nz) :: adz_shifty
      real(kind=8), dimension(nx,ny,nz) :: adx_shiftz
      real(kind=8), dimension(nx,ny,nz) :: ady_shiftz
      real(kind=8), dimension(nx,ny,nz) :: adz_shiftz
      real(kind=8), dimension(nx,ny,nz) :: adx_chi
      real(kind=8), dimension(nx,ny,nz) :: ady_chi
      real(kind=8), dimension(nx,ny,nz) :: adz_chi
      real(kind=8), dimension(nx,ny,nz) :: adx_Gamtx
      real(kind=8), dimension(nx,ny,nz) :: ady_Gamtx
      real(kind=8), dimension(nx,ny,nz) :: adz_Gamtx
      real(kind=8), dimension(nx,ny,nz) :: adx_Gamty
      real(kind=8), dimension(nx,ny,nz) :: ady_Gamty
      real(kind=8), dimension(nx,ny,nz) :: adz_Gamty
      real(kind=8), dimension(nx,ny,nz) :: adx_Gamtz
      real(kind=8), dimension(nx,ny,nz) :: ady_Gamtz
      real(kind=8), dimension(nx,ny,nz) :: adz_Gamtz
      real(kind=8), dimension(nx,ny,nz) :: adx_trK
      real(kind=8), dimension(nx,ny,nz) :: ady_trK
      real(kind=8), dimension(nx,ny,nz) :: adz_trK
      real(kind=8), dimension(nx,ny,nz) :: adx_gbx
      real(kind=8), dimension(nx,ny,nz) :: ady_gbx
      real(kind=8), dimension(nx,ny,nz) :: adz_gbx
      real(kind=8), dimension(nx,ny,nz) :: adx_gby
      real(kind=8), dimension(nx,ny,nz) :: ady_gby
      real(kind=8), dimension(nx,ny,nz) :: adz_gby
      real(kind=8), dimension(nx,ny,nz) :: adx_gbz
      real(kind=8), dimension(nx,ny,nz) :: ady_gbz
      real(kind=8), dimension(nx,ny,nz) :: adz_gbz
      real(kind=8), dimension(nx,ny,nz) :: adx_Ex
      real(kind=8), dimension(nx,ny,nz) :: ady_Ex
      real(kind=8), dimension(nx,ny,nz) :: adz_Ex
      real(kind=8), dimension(nx,ny,nz) :: adx_Ey
      real(kind=8), dimension(nx,ny,nz) :: ady_Ey
      real(kind=8), dimension(nx,ny,nz) :: adz_Ey
      real(kind=8), dimension(nx,ny,nz) :: adx_Ez
      real(kind=8), dimension(nx,ny,nz) :: ady_Ez
      real(kind=8), dimension(nx,ny,nz) :: adz_Ez
      real(kind=8), dimension(nx,ny,nz) :: adx_Bx
      real(kind=8), dimension(nx,ny,nz) :: ady_Bx
      real(kind=8), dimension(nx,ny,nz) :: adz_Bx
      real(kind=8), dimension(nx,ny,nz) :: adx_By
      real(kind=8), dimension(nx,ny,nz) :: ady_By
      real(kind=8), dimension(nx,ny,nz) :: adz_By
      real(kind=8), dimension(nx,ny,nz) :: adx_Bz
      real(kind=8), dimension(nx,ny,nz) :: ady_Bz
      real(kind=8), dimension(nx,ny,nz) :: adz_Bz
      real(kind=8), dimension(nx,ny,nz) :: adx_Phi_em
      real(kind=8), dimension(nx,ny,nz) :: ady_Phi_em
      real(kind=8), dimension(nx,ny,nz) :: adz_Phi_em
      real(kind=8), dimension(nx,ny,nz) :: adx_Psi_em
      real(kind=8), dimension(nx,ny,nz) :: ady_Psi_em
      real(kind=8), dimension(nx,ny,nz) :: adz_Psi_em
      ! }}}

        ! declare rhs functions {{{
        real(kind=8), dimension(nx,ny,nz) :: dt_Alpha3D
        real(kind=8), dimension(nx,ny,nz) :: dt_shiftx
        real(kind=8), dimension(nx,ny,nz) :: dt_shifty
        real(kind=8), dimension(nx,ny,nz) :: dt_shiftz
        real(kind=8), dimension(nx,ny,nz) :: dt_gtxx
        real(kind=8), dimension(nx,ny,nz) :: dt_gtxy
        real(kind=8), dimension(nx,ny,nz) :: dt_gtxz
        real(kind=8), dimension(nx,ny,nz) :: dt_gtyy
        real(kind=8), dimension(nx,ny,nz) :: dt_gtyz
        real(kind=8), dimension(nx,ny,nz) :: dt_gtzz
        real(kind=8), dimension(nx,ny,nz) :: dt_Atxx
        real(kind=8), dimension(nx,ny,nz) :: dt_Atxy
        real(kind=8), dimension(nx,ny,nz) :: dt_Atxz
        real(kind=8), dimension(nx,ny,nz) :: dt_Atyy
        real(kind=8), dimension(nx,ny,nz) :: dt_Atyz
        real(kind=8), dimension(nx,ny,nz) :: dt_Atzz
        real(kind=8), dimension(nx,ny,nz) :: dt_chi3D
        real(kind=8), dimension(nx,ny,nz) :: dt_trK3D
        real(kind=8), dimension(nx,ny,nz) :: dt_Gamtx
        real(kind=8), dimension(nx,ny,nz) :: dt_Gamty
        real(kind=8), dimension(nx,ny,nz) :: dt_Gamtz
        real(kind=8), dimension(nx,ny,nz) :: dt_gbx
        real(kind=8), dimension(nx,ny,nz) :: dt_gby
        real(kind=8), dimension(nx,ny,nz) :: dt_gbz

        real(kind=8), dimension(nx,ny,nz) :: shiftx
        real(kind=8), dimension(nx,ny,nz) :: shifty
        real(kind=8), dimension(nx,ny,nz) :: shiftz
        ! }}}

        integer, dimension(4)   :: lambda
        integer, dimension(2)   :: lambda_f
 
        real(kind=8),  dimension(3)       :: Betau
        real(kind=8),  dimension(3,3)     :: gtd_rhs, Atd_rhs 
        real(kind=8)                      :: chi_rhs, trK_rhs, Alpha_rhs
        real(kind=8),  dimension(3)       :: Bu_rhs, Betau_rhs, Gamt_rhs
        real(kind=8),  dimension(3)       :: adv_d_chi, adv_d_trK
        real(kind=8),  dimension(3)       :: adv_d_Alpha
        real(kind=8),  dimension(3,3)     :: adv_d_Betau, adv_d_Bu
        real(kind=8),  dimension(3,3)     :: adv_d_Gamt  
        real(kind=8),  dimension(3,3,3)   :: adv_d_gtd, adv_d_Atd

        integer                       :: lambda_1,  lambda_2
        integer                       :: lambda_3,  lambda_4
        integer                       :: evolve_geometry
        integer                       :: evolve_em_field

        real(kind=8)                      :: t1, t2, t4, t7  ! maple work vars
        integer                       ::  i, j, k, ii
        real(kind=8)                      ::  dx, dy, dz
        integer                       :: rc, err
        logical, parameter             :: error_check = .false.
        logical, parameter             :: ltrace      = .true.
        logical, parameter             :: ltrace2     = .true.

        if (ltrace2) then
          write(0,*) '### Begin cal_bssn_rhs_lie'
        end if

        dx = x1d(2) - x1d(1)
        dy = y1d(2) - y1d(1)
        dz = z1d(2) - z1d(1)

        lambda_1 = lambda(1)
        lambda_2 = lambda(2)
        lambda_3 = lambda(3)
        lambda_4 = lambda(4)

        evolve_geometry = 1

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          !print *,'i,j,k', i, j, k
          Betau(1) = shiftx(i,j,k)
          Betau(2) = shifty(i,j,k)
          Betau(3) = shiftz(i,j,k)

          !print *,'gtd: i,j,k', i, j, k
          gtd_rhs(1,1) = dt_gtxx(i,j,k)
          gtd_rhs(1,2) = dt_gtxy(i,j,k)
          gtd_rhs(1,3) = dt_gtxz(i,j,k)
          gtd_rhs(2,1) = dt_gtxy(i,j,k)
          gtd_rhs(2,2) = dt_gtyy(i,j,k)
          gtd_rhs(2,3) = dt_gtyz(i,j,k)
          gtd_rhs(3,1) = dt_gtxz(i,j,k)
          gtd_rhs(3,2) = dt_gtyz(i,j,k)
          gtd_rhs(3,3) = dt_gtzz(i,j,k)

          Atd_rhs(1,1) = dt_Atxx(i,j,k)
          Atd_rhs(1,2) = dt_Atxy(i,j,k)
          Atd_rhs(1,3) = dt_Atxz(i,j,k)
          Atd_rhs(2,1) = dt_Atxy(i,j,k)
          Atd_rhs(2,2) = dt_Atyy(i,j,k)
          Atd_rhs(2,3) = dt_Atyz(i,j,k)
          Atd_rhs(3,1) = dt_Atxz(i,j,k)
          Atd_rhs(3,2) = dt_Atyz(i,j,k)
          Atd_rhs(3,3) = dt_Atzz(i,j,k)

          trK_rhs      = dt_trK3D(i,j,k)
          chi_rhs      = dt_chi3D(i,j,k)

          Gamt_rhs(1) = dt_Gamtx(i,j,k)
          Gamt_rhs(2) = dt_Gamty(i,j,k)
          Gamt_rhs(3) = dt_Gamtz(i,j,k)

          Alpha_rhs    = dt_Alpha3D(i,j,k)
          Betau_rhs(1) = dt_shiftx(i,j,k)
          Betau_rhs(2) = dt_shifty(i,j,k)
          Betau_rhs(3) = dt_shiftz(i,j,k)

          Bu_rhs(1)    = dt_gbx(i,j,k)
          Bu_rhs(2)    = dt_gby(i,j,k)
          Bu_rhs(3)    = dt_gbz(i,j,k)

          !print *,'Gamt: i,j,k', i, j, k
          adv_d_Gamt(1,1) = adx_Gamtx(i,j,k)
          adv_d_Gamt(1,2) = adx_Gamty(i,j,k)
          adv_d_Gamt(1,3) = adx_Gamtz(i,j,k)
          adv_d_Gamt(2,1) = ady_Gamtx(i,j,k)
          adv_d_Gamt(2,2) = ady_Gamty(i,j,k)
          adv_d_Gamt(2,3) = ady_Gamtz(i,j,k)
          adv_d_Gamt(3,1) = adz_Gamtx(i,j,k)
          adv_d_Gamt(3,2) = adz_Gamty(i,j,k)
          adv_d_Gamt(3,3) = adz_Gamtz(i,j,k)

          adv_d_Atd(1,1,1) = adx_Atxx(i,j,k)
          adv_d_Atd(2,1,1) = ady_Atxx(i,j,k)
          adv_d_Atd(3,1,1) = adz_Atxx(i,j,k)
          adv_d_Atd(1,1,2) = adx_Atxy(i,j,k)
          adv_d_Atd(2,1,2) = ady_Atxy(i,j,k)
          adv_d_Atd(3,1,2) = adz_Atxy(i,j,k)
          adv_d_Atd(1,1,3) = adx_Atxz(i,j,k)
          adv_d_Atd(2,1,3) = ady_Atxz(i,j,k)
          adv_d_Atd(3,1,3) = adz_Atxz(i,j,k)
          adv_d_Atd(1,2,1) = adv_d_Atd(1,1,2)
          adv_d_Atd(2,2,1) = adv_d_Atd(2,1,2)
          adv_d_Atd(3,2,1) = adv_d_Atd(3,1,2)
          adv_d_Atd(1,2,2) = adx_Atyy(i,j,k)
          adv_d_Atd(2,2,2) = ady_Atyy(i,j,k)
          adv_d_Atd(3,2,2) = adz_Atyy(i,j,k)
          adv_d_Atd(1,2,3) = adx_Atyz(i,j,k)
          adv_d_Atd(2,2,3) = ady_Atyz(i,j,k)
          adv_d_Atd(3,2,3) = adz_Atyz(i,j,k)
          adv_d_Atd(1,3,1) = adv_d_Atd(1,1,3)
          adv_d_Atd(2,3,1) = adv_d_Atd(2,1,3)
          adv_d_Atd(3,3,1) = adv_d_Atd(3,1,3)
          adv_d_Atd(1,3,2) = adv_d_Atd(1,2,3)
          adv_d_Atd(2,3,2) = adv_d_Atd(2,2,3)
          adv_d_Atd(3,3,2) = adv_d_Atd(3,2,3)
          adv_d_Atd(1,3,3) = adx_Atzz(i,j,k)
          adv_d_Atd(2,3,3) = ady_Atzz(i,j,k)
          adv_d_Atd(3,3,3) = adz_Atzz(i,j,k)

          adv_d_gtd(1,1,1) = adx_gtxx(i,j,k)
          adv_d_gtd(2,1,1) = ady_gtxx(i,j,k)
          adv_d_gtd(3,1,1) = adz_gtxx(i,j,k)
          adv_d_gtd(1,1,2) = adx_gtxy(i,j,k)
          adv_d_gtd(2,1,2) = ady_gtxy(i,j,k)
          adv_d_gtd(3,1,2) = adz_gtxy(i,j,k)
          adv_d_gtd(1,1,3) = adx_gtxz(i,j,k)
          adv_d_gtd(2,1,3) = ady_gtxz(i,j,k)
          adv_d_gtd(3,1,3) = adz_gtxz(i,j,k)
          adv_d_gtd(1,2,1) = adv_d_gtd(1,1,2)
          adv_d_gtd(2,2,1) = adv_d_gtd(2,1,2)
          adv_d_gtd(3,2,1) = adv_d_gtd(3,1,2)
          adv_d_gtd(1,2,2) = adx_gtyy(i,j,k)
          adv_d_gtd(2,2,2) = ady_gtyy(i,j,k)
          adv_d_gtd(3,2,2) = adz_gtyy(i,j,k)
          adv_d_gtd(1,2,3) = adx_gtyz(i,j,k)
          adv_d_gtd(2,2,3) = ady_gtyz(i,j,k)
          adv_d_gtd(3,2,3) = adz_gtyz(i,j,k)
          adv_d_gtd(1,3,1) = adv_d_gtd(1,1,3)
          adv_d_gtd(2,3,1) = adv_d_gtd(2,1,3)
          adv_d_gtd(3,3,1) = adv_d_gtd(3,1,3)
          adv_d_gtd(1,3,2) = adv_d_gtd(1,2,3)
          adv_d_gtd(2,3,2) = adv_d_gtd(2,2,3)
          adv_d_gtd(3,3,2) = adv_d_gtd(3,2,3)
          adv_d_gtd(1,3,3) = adx_gtzz(i,j,k)
          adv_d_gtd(2,3,3) = ady_gtzz(i,j,k)
          adv_d_gtd(3,3,3) = adz_gtzz(i,j,k)

          adv_d_chi(1)     = adx_chi(i,j,k)
          adv_d_chi(2)     = ady_chi(i,j,k)
          adv_d_chi(3)     = adz_chi(i,j,k)

          adv_d_trK(1)     = adx_trK(i,j,k)
          adv_d_trK(2)     = ady_trK(i,j,k)
          adv_d_trK(3)     = adz_trK(i,j,k)

          adv_d_Alpha(1)   = adx_alpha(i,j,k)
          adv_d_Alpha(2)   = ady_alpha(i,j,k)
          adv_d_Alpha(3)   = adz_alpha(i,j,k)

          !print *,'Betau: i,j,k', i, j, k
          adv_d_Betau(1,1) = adx_shiftx(i,j,k)
          adv_d_Betau(2,1) = ady_shiftx(i,j,k)
          adv_d_Betau(3,1) = adz_shiftx(i,j,k)
          adv_d_Betau(1,2) = adx_shifty(i,j,k)
          adv_d_Betau(2,2) = ady_shifty(i,j,k)
          adv_d_Betau(3,2) = adz_shifty(i,j,k)
          adv_d_Betau(1,3) = adx_shiftz(i,j,k)
          adv_d_Betau(2,3) = ady_shiftz(i,j,k)
          adv_d_Betau(3,3) = adz_shiftz(i,j,k)

          !print *,'Bu: i,j,k', i, j, k
          adv_d_Bu(1,1)    = adx_gbx(i,j,k)
          adv_d_Bu(2,1)    = ady_gbx(i,j,k)
          adv_d_Bu(3,1)    = adz_gbx(i,j,k)
          adv_d_Bu(1,2)    = adx_gby(i,j,k)
          adv_d_Bu(2,2)    = ady_gby(i,j,k)
          adv_d_Bu(3,2)    = adz_gby(i,j,k)
          adv_d_Bu(1,3)    = adx_gbz(i,j,k)
          adv_d_Bu(2,3)    = ady_gbz(i,j,k)
          adv_d_Bu(3,3)    = adz_gbz(i,j,k)

#include "bssn_emtest_adv.h"

            dt_gtxx(i,j,k) = gtd_rhs(1,1)
            dt_gtxy(i,j,k) = gtd_rhs(1,2)
            dt_gtxz(i,j,k) = gtd_rhs(1,3)
            dt_gtyy(i,j,k) = gtd_rhs(2,2)
            dt_gtyz(i,j,k) = gtd_rhs(2,3)
            dt_gtzz(i,j,k) = gtd_rhs(3,3)
  
            dt_Atxx(i,j,k) = Atd_rhs(1,1)
            dt_Atxy(i,j,k) = Atd_rhs(1,2)
            dt_Atxz(i,j,k) = Atd_rhs(1,3)
            dt_Atyy(i,j,k) = Atd_rhs(2,2)
            dt_Atyz(i,j,k) = Atd_rhs(2,3)
            dt_Atzz(i,j,k) = Atd_rhs(3,3)
  
            dt_trK3D(i,j,k) = trK_rhs
            dt_chi3D(i,j,k) = chi_rhs
  
            dt_Gamtx(i,j,k) = Gamt_rhs(1)
            dt_Gamty(i,j,k) = Gamt_rhs(2)
            dt_Gamtz(i,j,k) = Gamt_rhs(3)
            
            dt_Alpha3D(i,j,k) = Alpha_rhs
            dt_shiftx(i,j,k)  = Betau_rhs(1)
            dt_shifty(i,j,k)  = Betau_rhs(2)
            dt_shiftz(i,j,k)  = Betau_rhs(3)
            
            !print *,'update rhs gb: i,j,k', i, j, k
            dt_gbx(i,j,k) = Bu_rhs(1)
            dt_gby(i,j,k) = Bu_rhs(2)
            dt_gbz(i,j,k) = Bu_rhs(3)

        end do
        end do
        end do


        return
      end  subroutine

!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
