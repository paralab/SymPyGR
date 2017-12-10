!  #include "cctk.h"
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
      subroutine cal_bssn_rhs(
     &  dt_Alpha3D, dt_shiftx, dt_shifty, dt_shiftz, dt_chi3D,
     &  dt_trK3D, dt_gtxx, dt_gtxy, dt_gtxz, dt_gtyy, dt_gtyz,
     &  dt_gtzz, dt_Atxx, dt_Atxy, dt_Atxz, dt_Atyy, dt_Atyz,
     &  dt_Atzz, dt_Gamtx, dt_Gamty, dt_Gamtz, dt_gbx, dt_gby, dt_gbz,
     &  Alpha3D, shiftx, shifty, shiftz, chi3D, trK3D, 
     &  gtxx, gtxy, gtxz, gtyy, gtyz, gtzz, 
     &  Atxx, Atxy, Atxz, Atyy, Atyz,
     &  Atzz, Gamtx, Gamty, Gamtz, gbx, gby, gbz,
     &  lambda, lambda_f,
     &  eta, eta_damping, eta_damping_exp, R_0,
     &  trk0, chi_floor,
     &  x1d, y1d, z1d, nx, ny, nz )

        implicit none

        integer, intent(in)                       :: nx, ny, nz
        real(kind=8), dimension(nx), intent(in)       :: x1d(nx)
        real(kind=8), dimension(ny), intent(in)       :: y1d(ny)
        real(kind=8), dimension(nz), intent(in)       :: z1d(nz)

        real(kind=8), dimension(nx,ny,nz), intent(in) :: Alpha3D
        real(kind=8), dimension(nx,ny,nz), intent(in) :: shiftx
        real(kind=8), dimension(nx,ny,nz), intent(in) :: shifty
        real(kind=8), dimension(nx,ny,nz), intent(in) :: shiftz
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gtxx
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gtxy
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gtxz
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gtyy
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gtyz
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gtzz
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Atxx 
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Atxy 
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Atxz
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Atyy 
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Atyz 
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Atzz
        real(kind=8), dimension(nx,ny,nz), intent(in) :: chi3D
        real(kind=8), dimension(nx,ny,nz), intent(in) :: trK3D
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Gamtx
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Gamty
        real(kind=8), dimension(nx,ny,nz), intent(in) :: Gamtz
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gbx
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gby
        real(kind=8), dimension(nx,ny,nz), intent(in) :: gbz

        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Alpha3D
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_shiftx
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_shifty 
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_shiftz
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gtxx
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gtxy
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gtxz
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gtyy
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gtyz
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gtzz
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Atxx
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Atxy
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Atxz
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Atyy
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Atyz
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Atzz
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_chi3D
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_trK3D
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Gamtx
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Gamty
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_Gamtz
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gbx
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gby
        real(kind=8), dimension(nx,ny,nz), intent(out) :: dt_gbz

        real(kind=8)            :: chi_floor, trK0 
        real(kind=8)            :: eta
        integer                 :: eta_damping, eta_damping_exp
        integer, dimension(4)   :: lambda
        integer, dimension(2)   :: lambda_f


        ! local variables
        integer                       :: i, j, k
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
!        real(kind=8),  dimension(3,3)     :: DelDel_Alpha, pTtd_ADM
        real(kind=8),  dimension(3,3)     :: pTtd_ADM
        real(kind=8),  dimension(3,3)     :: dd_Alpha, d_Gamt 
        real(kind=8),  dimension(3,3)     :: d_Betau, Psi1, Psi1TF, d_Bu
        real(kind=8),  dimension(3,3)     :: Atd_rhs, gtd_rhs
        real(kind=8),  dimension(3,3,3)   :: dd_Betau, d_gtd, d_Atd
        real(kind=8),  dimension(3,3,3,3) :: dd_gtd 
        real(kind=8),  dimension(0:3,0:3) :: gd4, gu4, Tu
        real(kind=8)                      :: detgd, idetgd
        real(kind=8)                      :: detgtd, idetgtd
        real(kind=8),  dimension(3,3,3)   :: Ctd, Ct, Chr
        real(kind=8)                      :: twelve_Pi_G
        real(kind=8)                      :: third_trPsi1
!        real(kind=8)                      :: expphi,  exp2phi,  exp4phi  ! set?
!        real(kind=8)                      :: expmphi, expm2phi, expm4phi ! set?
        real(kind=8)                      :: inv_chi, sqrt_chi 
!        real(kind=8)                      :: inv_sqrt_chi  

        real(kind=8)                      :: divE, divB, Bsq, inv_Bsq 
        real(kind=8)                      :: div_Beta  
        real(kind=8), dimension(3)        :: d_div_Beta
 
        real(kind=8)                      :: r_coord
 
        real(kind=8)                      :: det_gtd, idet_gtd
        real(kind=8)                      :: det_gtd_to_third
        real(kind=8)                      :: det_gtd_to_neg_third
        real(kind=8)                      :: trace_Atd
        real(kind=8)                      :: neg_one_third_trace_Atd

        integer                       :: lambda_1,  lambda_2
        integer                       :: lambda_3,  lambda_4
        integer                       :: lambda_f0, lambda_f1

        real(kind=8)                      :: R_0
        real(kind=8)                      :: feta
  
        real(kind=8)                      :: kappa_1, kappa_2 !! for the EM damping
        integer                       :: force_free 

        integer                       :: evolve_geometry
        real(kind=8)                      :: evolve_gr

        ! debugging
        real(kind=8),  dimension(3)       :: adv_d_chi, adv_d_trK
        real(kind=8),  dimension(3)       :: adv_d_Alpha
        real(kind=8),  dimension(3,3)     :: adv_d_Betau, adv_d_Bu
        real(kind=8),  dimension(3,3)     :: adv_d_Gamt
        real(kind=8),  dimension(3,3,3)   :: adv_d_gtd, adv_d_Atd
        real(kind=8)                      :: trK2_rhs, trK3_rhs
        real(kind=8) :: t32

!       real(kind=8), dimension(3)        :: adv_d_gtxx
!       real(kind=8), dimension(3)        :: adv_d_gtxy
!       real(kind=8), dimension(3)        :: adv_d_gtxz
!       real(kind=8), dimension(3)        :: adv_d_gtyy
!       real(kind=8), dimension(3)        :: adv_d_gtyz
!       real(kind=8), dimension(3)        :: adv_d_gtzz
!       real(kind=8), dimension(3)        :: adv_d_Atxx
!       real(kind=8), dimension(3)        :: adv_d_Atxy
!       real(kind=8), dimension(3)        :: adv_d_Atxz
!       real(kind=8), dimension(3)        :: adv_d_Atyy
!       real(kind=8), dimension(3)        :: adv_d_Atyz
!       real(kind=8), dimension(3)        :: adv_d_Atzz
!       real(kind=8), dimension(3)        :: adv_d_shiftx
!       real(kind=8), dimension(3)        :: adv_d_shifty
!       real(kind=8), dimension(3)        :: adv_d_shiftz
!       real(kind=8), dimension(3)        :: adv_d_Gamtz
!       real(kind=8), dimension(3)        :: adv_d_Gamtx
!       real(kind=8), dimension(3)        :: adv_d_Gamty
!       real(kind=8), dimension(3)        :: adv_d_gbx
!       real(kind=8), dimension(3)        :: adv_d_gby
!       real(kind=8), dimension(3)        :: adv_d_gbz
!       real(kind=8), dimension(3,3)      :: dd_shiftz, dd_shifty, dd_shiftx
!       real(kind=8), dimension(3,3)      :: dd_gtxx, dd_gtxy, dd_gtxz
!       real(kind=8), dimension(3,3)      :: dd_gtyy, dd_gtyz, dd_gtzz
!#include "include/new_derivs_defn.h"


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
        logical, parameter             :: error_check = .false.

        real(kind=8), parameter           :: eps = 1.0d-9
        logical, parameter             :: debug = .true.

!       This variable should be:  turn_off = .false.
        logical, parameter             :: turn_off    = .false.

        ! include declarations for Maple temp variables
#include "bssn_emtest_temp_vars.h"


!       real(kind=8)                      ::  ufuncs(7,7,7)
        integer                       ::  ii, jj, kk, mm
        real(kind=8)                      ::  idx, idxsq

        if (turn_off) then
          return
        end if
 
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

        idx = 1.0d0/dx
        idxsq = idx*idx

        evolve_geometry = 1

        if (evolve_geometry .eq. 1) then
          evolve_gr = 1.0d0
        else
          evolve_gr = 0.0d0
        end if

        if ( ltrace ) then
            write(0,*)' cal_bssn_rhs:  Checking parameters '
            write(0,*)'           These should be integers: '
            write(0,*)'                lambda_1  = ', lambda_1
            write(0,*)'                lambda_2  = ', lambda_2
            write(0,*)'                lambda_3  = ', lambda_3
            write(0,*)'                lambda_4  = ', lambda_4
            write(0,*)'                lambda_f0 = ', lambda_f0
            write(0,*)'                lambda_f1 = ', lambda_f1
            write(0,*)'           These should be reals: '
            write(0,*)'                eta       = ', eta
            write(0,*)'                trK0      = ', trK0
            write(0,*)'                kappa_1   = ', kappa_1
            write(0,*)'                kappa_2   = ', kappa_1
            write(0,*)'                chi_floor = ', chi_floor
        endif

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

!       Jd(1)   = 0.0d0
!       Jd(2)   = 0.0d0
!       Jd(3)   = 0.0d0
!       Sd(1,1)   = 0.0d0
!       Sd(2,1)   = 0.0d0
!       Sd(3,1)   = 0.0d0
!       Sd(1,2)   = 0.0d0
!       Sd(2,2)   = 0.0d0
!       Sd(3,2)   = 0.0d0
!       Sd(1,3)   = 0.0d0
!       Sd(2,3)   = 0.0d0
!       Sd(3,3)   = 0.0d0

!       x3d         => w(G_XPHYS)%d 
!       y3d         => w(G_YPHYS)%d 
!       z3d         => w(G_ZPHYS)%d 

!       Alpha3D     => u(G_ALPHA)%d
!       shiftx      => u(G_SHIFT1)%d
!       shifty      => u(G_SHIFT2)%d
!       shiftz      => u(G_SHIFT3)%d
!       chi3D       => u(G_CHI)%d
!       trK3D       => u(G_TRK)%d
!       gtxx        => u(G_GT11)%d
!       gtxy        => u(G_GT12)%d
!       gtxz        => u(G_GT13)%d
!       gtyy        => u(G_GT22)%d
!       gtyz        => u(G_GT23)%d
!       gtzz        => u(G_GT33)%d
!       Atxx        => u(G_A11)%d
!       Atxy        => u(G_A12)%d
!       Atxz        => u(G_A13)%d
!       Atyy        => u(G_A22)%d
!       Atyz        => u(G_A23)%d
!       Atzz        => u(G_A33)%d
!       Gamtx       => u(G_GAM1)%d
!       Gamty       => u(G_GAM2)%d
!       Gamtz       => u(G_GAM3)%d
!       gbx         => u(G_GB1)%d
!       gby         => u(G_GB2)%d
!       gbz         => u(G_GB3)%d

!       dt_Alpha3D     => dtu(G_ALPHA)%d
!       dt_shiftx      => dtu(G_SHIFT1)%d
!       dt_shifty      => dtu(G_SHIFT2)%d
!       dt_shiftz      => dtu(G_SHIFT3)%d
!       dt_chi3D       => dtu(G_CHI)%d
!       dt_trK3D       => dtu(G_TRK)%d
!       dt_gtxx        => dtu(G_GT11)%d
!       dt_gtxy        => dtu(G_GT12)%d
!       dt_gtxz        => dtu(G_GT13)%d
!       dt_gtyy        => dtu(G_GT22)%d
!       dt_gtyz        => dtu(G_GT23)%d
!       dt_gtzz        => dtu(G_GT33)%d
!       dt_Atxx        => dtu(G_A11)%d
!       dt_Atxy        => dtu(G_A12)%d
!       dt_Atxz        => dtu(G_A13)%d
!       dt_Atyy        => dtu(G_A22)%d
!       dt_Atyz        => dtu(G_A23)%d
!       dt_Atzz        => dtu(G_A33)%d
!       dt_Gamtx       => dtu(G_GAM1)%d
!       dt_Gamty       => dtu(G_GAM2)%d
!       dt_Gamtz       => dtu(G_GAM3)%d
!       dt_gbx         => dtu(G_GB1)%d
!       dt_gby         => dtu(G_GB2)%d
!       dt_gbz         => dtu(G_GB3)%d

!       gxx         => v(G_G11)%d
!       gxy         => v(G_G12)%d
!       gxz         => v(G_G13)%d
!       gyy         => v(G_G22)%d
!       gyz         => v(G_G23)%d
!       gzz         => v(G_G33)%d

        ! EM fields 
!       Ex          =>  u(G_EX)%d 
!       Ey          =>  u(G_EY)%d 
!       Ez          =>  u(G_EZ)%d 
!       Bx          =>  u(G_BX)%d 
!       By          =>  u(G_BY)%d 
!       Bz          =>  u(G_BZ)%d 
!       Psi_em      =>  u(G_PSI_EM)%d 
!       Phi_em      =>  u(G_PHI_EM)%d 

!       dt_Ex          =>  dtu(G_EX)%d 
!       dt_Ey          =>  dtu(G_EY)%d 
!       dt_Ez          =>  dtu(G_EZ)%d 
!       dt_Bx          =>  dtu(G_BX)%d 
!       dt_By          =>  dtu(G_BY)%d 
!       dt_Bz          =>  dtu(G_BZ)%d 
!       dt_Psi_em      =>  dtu(G_PSI_EM)%d 
!       dt_Phi_em      =>  dtu(G_PHI_EM)%d 
        

!----------------------------------------------------
! L O O P   B E G I N S   H E R E
!----------------------------------------------------

!!     dir$ simd
        do k = 4, nz-3
        do j = 4, ny-3
        do i = 4, nx-3

!         do mm = 1, NU_G
!         do kk = -3, 3
!         do jj = -3, 3
!         do ii = -3, 3
!           ufuncs(4+ii,4+jj,4+kk,mm) = 
!    &                         u(mm)%d(i+ii,j+jj,k+kk)
!         end do
!         end do
!         end do
!         end do

#if 0
          do mm = 1, NU_G

            if ((mm .ge. G_GT11 .and. mm .le. G_GT33) .or.
     &           mm .eq. G_ALPHA .or.
     &           mm .eq. G_CHI .or.
     &          (mm .ge. G_SHIFT1 .and. mm .le. G_SHIFT3)) then

!           load the x-y plane
            kk = 0
            do jj = -2, 2
            do ii = -2, 2
              ufuncs(4+ii,4+jj,4+kk,mm) =
     &                           u(mm)%d(i+ii,j+jj,k+kk)
            end do
            end do
            ufuncs(4,1,4,mm) = u(mm)%d(i,j-3,k)
            ufuncs(4,7,4,mm) = u(mm)%d(i,j+3,k)

!           load the y-z plane
            ii = 0
            do kk = -2, 2
            do jj = -2, 2
              ufuncs(4+ii,4+jj,4+kk,mm) =
     &                           u(mm)%d(i+ii,j+jj,k+kk)
            end do
            end do
            ufuncs(4,4,1,mm) = u(mm)%d(i,j,k-3)
            ufuncs(4,4,7,mm) = u(mm)%d(i,j,k+3)

!           load the x-z plane
            jj = 0
            do kk = -2, 2
            do ii = -2, 2
              ufuncs(4+ii,4+jj,4+kk,mm) =
     &                           u(mm)%d(i+ii,j+jj,k+kk)
            end do
            end do
            ufuncs(1,4,4,mm) = u(mm)%d(i-3,j,k)
            ufuncs(7,4,4,mm) = u(mm)%d(i+3,j,k)

          else

!           load the x-line
            do ii = -2, 2
              ufuncs(4+ii,4,4,mm) =
     &                           u(mm)%d(i+ii,j,k)
            end do
            ufuncs(1,4,4,mm) = u(mm)%d(i-3,j,k)
            ufuncs(7,4,4,mm) = u(mm)%d(i+3,j,k)

!           load the y-line
            do jj = -2, 2
              ufuncs(4,4+jj,4,mm) =
     &                           u(mm)%d(i,j+jj,k)
            end do
            ufuncs(4,1,4,mm) = u(mm)%d(i,j-3,k)
            ufuncs(4,7,4,mm) = u(mm)%d(i,j+3,k)

!           load the z-line
            do kk = -2, 2
              ufuncs(4,4,4+kk,mm) =
     &                           u(mm)%d(i,j,k+kk)
            end do
            ufuncs(4,4,1,mm) = u(mm)%d(i,j,k-3)
            ufuncs(4,4,7,mm) = u(mm)%d(i,j,k+3)
          end if
          end do
#endif

#include "all_derivs.h"

        ! Construct a radial coordinate from the Cartesian coordinates 
        ! in the usual way.  We need this if we damp the eta factor in 
        ! the gauge conditions.  (See arxiv:1003.0859) 
          r_coord = sqrt(x1d(i)*x1d(i) + y1d(j)*y1d(j) + z1d(k)*z1d(k))

          ! define somehow the value of feta, the damped eta
          feta = eta
          ! DWN - removing for vectorization
          if(r_coord.ge.R_0) feta=eta*(R_0/r_coord)**eta_damping_exp

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


#include "bssn_emtest.h" 

#include "bssn_emtest_adv.h"

          dt_Alpha3D(i,j,k)   = evolve_gr * Alpha_rhs
          dt_shiftx(i,j,k)    = evolve_gr * Betau_rhs(1)
          dt_shifty(i,j,k)    = evolve_gr * Betau_rhs(2)
          dt_shiftz(i,j,k)    = evolve_gr * Betau_rhs(3)

          dt_chi3D(i,j,k)     = evolve_gr * chi_rhs
          dt_trK3D(i,j,k)     = evolve_gr * trK_rhs

          dt_gtxx(i,j,k)      = evolve_gr * gtd_rhs(1,1)
          dt_gtxy(i,j,k)      = evolve_gr * gtd_rhs(1,2)
          dt_gtxz(i,j,k)      = evolve_gr * gtd_rhs(1,3)
          dt_gtyy(i,j,k)      = evolve_gr * gtd_rhs(2,2)
          dt_gtyz(i,j,k)      = evolve_gr * gtd_rhs(2,3)
          dt_gtzz(i,j,k)      = evolve_gr * gtd_rhs(3,3)

          dt_Atxx(i,j,k)      = evolve_gr * Atd_rhs(1,1)
          dt_Atxy(i,j,k)      = evolve_gr * Atd_rhs(1,2)
          dt_Atxz(i,j,k)      = evolve_gr * Atd_rhs(1,3)
          dt_Atyy(i,j,k)      = evolve_gr * Atd_rhs(2,2)
          dt_Atyz(i,j,k)      = evolve_gr * Atd_rhs(2,3)
          dt_Atzz(i,j,k)      = evolve_gr * Atd_rhs(3,3)

          dt_Gamtx(i,j,k)     = evolve_gr * Gamt_rhs(1)
          dt_Gamty(i,j,k)     = evolve_gr * Gamt_rhs(2)
          dt_Gamtz(i,j,k)     = evolve_gr * Gamt_rhs(3)
  
          dt_gbx(i,j,k)       = evolve_gr * Bu_rhs(1)
          dt_gby(i,j,k)       = evolve_gr * Bu_rhs(2)
          dt_gbz(i,j,k)       = evolve_gr * Bu_rhs(3)


          if (error_check) then
            err = 1
            call check_isfinite(dt_Alpha3D(i,j,k), 'dt_Alpha3d', rc)
            err = err*rc
            call check_isfinite(dt_shiftx(i,j,k), 'dt_shiftx', rc)
            err = err*rc
            call check_isfinite(dt_shifty(i,j,k), 'dt_shifty', rc)
            err = err*rc
            call check_isfinite(dt_shiftz(i,j,k), 'dt_shiftz', rc)
            err = err*rc
  
            call check_isfinite(dt_chi3D(i,j,k), 'dt_chi3D', rc)
            err = err*rc
            call check_isfinite(dt_trK3D(i,j,k), 'dt_trK3D', rc)
            err = err*rc
  
            call check_isfinite(dt_gtxx(i,j,k), 'dt_gtxx', rc)
            err = err*rc
            call check_isfinite(dt_gtxy(i,j,k), 'dt_gtxy', rc)
            err = err*rc
            call check_isfinite(dt_gtxz(i,j,k), 'dt_gtxz', rc)
            err = err*rc
            call check_isfinite(dt_gtyy(i,j,k), 'dt_gtyy', rc)
            err = err*rc
            call check_isfinite(dt_gtyz(i,j,k), 'dt_gtyz', rc)
            err = err*rc
            call check_isfinite(dt_gtzz(i,j,k), 'dt_gtzz', rc)
            err = err*rc
  
            call check_isfinite(dt_Atxx(i,j,k), 'dt_Atxx', rc)
            err = err*rc
            call check_isfinite(dt_Atxy(i,j,k), 'dt_Atxy', rc)
            err = err*rc
            call check_isfinite(dt_Atxz(i,j,k), 'dt_Atxz', rc)
            err = err*rc
            call check_isfinite(dt_Atyy(i,j,k), 'dt_Atyy', rc)
            err = err*rc
            call check_isfinite(dt_Atyz(i,j,k), 'dt_Atyz', rc)
            err = err*rc
            call check_isfinite(dt_Atzz(i,j,k), 'dt_Atzz', rc)
            err = err*rc
  
            call check_isfinite(dt_Gamtx(i,j,k), 'dt_Gamtx', rc)
            err = err*rc
            call check_isfinite(dt_Gamty(i,j,k), 'dt_Gamty', rc)
            err = err*rc
            call check_isfinite(dt_Gamtz(i,j,k), 'dt_Gamtz', rc)
            err = err*rc
  
            call check_isfinite(dt_gbx(i,j,k), 'dt_gbx', rc)
            err = err*rc
            call check_isfinite(dt_gby(i,j,k), 'dt_gby', rc)
            err = err*rc
            call check_isfinite(dt_gbz(i,j,k), 'dt_gbz', rc)
            err = err*rc
  
            if (err .eq. 0) then
              write(*,*)'### Problem with update'
              write(*,*)' location i,j,k: ',i,j,k
              write(*,*)' location nx,ny,nz: ',nx,ny,nz
  
              write(*,*)' gtd ', gtd(1,1), gtd(1,2), gtd(1,3)
              write(*,*)' gtd ', gtd(2,1), gtd(2,2), gtd(2,3)
              write(*,*)' gtd ', gtd(3,1), gtd(3,2), gtd(3,3)
  
              write(*,*)' Atd ', Atd(1,1), Atd(1,2), Atd(1,3)
              write(*,*)' Atd ', Atd(2,1), Atd(2,2), Atd(2,3)
              write(*,*)' Atd ', Atd(3,1), Atd(3,2), Atd(3,3)
  
              write(*,*)' trK,chi', trK, chi
  
              write(*,*)' Gamt ', Gamt(1), Gamt(2), Gamt(3)
  
              write(*,*)' Alpha ', Alpha
              write(*,*)' Betau ', Betau(1), Betau(2), Betau(3)
  
              write(*,*)' Tu ', Tu(0,0), Tu(0,1), Tu(0,2)
              write(*,*)' Tu ', Tu(0,3), Tu(1,0), Tu(1,1)
              write(*,*)' Tu ', Tu(1,2), Tu(1,3), Tu(2,0)
              write(*,*)' Tu ', Tu(2,1), Tu(2,2), Tu(2,3)
              write(*,*)' Tu ', Tu(3,0), Tu(3,1), Tu(3,2)
              write(*,*)' Tu ', Tu(3,3)
  
       
              write(*,*) ' d_Alpha ', d_Alpha(1), d_Alpha(2), d_Alpha(3)
  
              write(*,*) ' dd_Alpha ', dd_Alpha(1,1), dd_Alpha(1,2), 
     &                                 dd_Alpha(1,3)
              write(*,*)' dd_Alpha ', dd_Alpha(2,2), dd_Alpha(2,3),
     &                                dd_Alpha(3,3)
  
              write(*,*)' d_chi ', d_chi(1), d_chi(2), d_chi(3)
  
              write(*,*)' dd_chi ', dd_chi(1,1),dd_chi(1,2), dd_chi(1,3)
              write(*,*)' dd_chi ', dd_chi(2,2),dd_chi(2,3), dd_chi(3,3)

              write(*,*)' d_trK ', d_trK(1), d_trK(2), d_trK(3)

              write(*,*)' d_Gamt ', d_Gamt(1,1),d_Gamt(2,1), d_Gamt(3,1)
              write(*,*)' d_Gamt ', d_Gamt(1,2),d_Gamt(2,2), d_Gamt(3,2)
              write(*,*)' d_Gamt ', d_Gamt(1,3),d_Gamt(2,3), d_Gamt(3,3)

              write(*,*)' d_Betau 1'
              write(*,FMT3) d_Betau(1,1), d_Betau(2,1), d_Betau(3,1)
              write(*,*)' d_Betau 2'
              write(*,FMT3) d_Betau(1,2), d_Betau(2,2), d_Betau(3,2)
              write(*,*)' d_Betau 3'
              write(*,FMT3) d_Betau(1,3), d_Betau(2,3), d_Betau(3,3)

              write(*,*)' dd_Betau 1'
              write(*,FMT6) dd_Betau(1,1,1), dd_Betau(1,2,1),
     &                    dd_Betau(1,3,1), dd_Betau(2,1,1),
     &                    dd_Betau(2,2,1), dd_Betau(2,3,1)

              write(*,FMT3) dd_Betau(3,1,1), dd_Betau(3,2,1),
     &                    dd_Betau(3,3,1)

              write(*,*)' dd_Betau 2'
              write(*,FMT6) dd_Betau(1,1,2), dd_Betau(1,2,2),
     &                    dd_Betau(1,3,2), dd_Betau(2,1,2),
     &                    dd_Betau(2,2,2), dd_Betau(2,3,2)
              write(*,FMT3) dd_Betau(3,1,2), dd_Betau(3,2,2),
     &                    dd_Betau(3,3,2)

              write(*,*)' dd_Betau 2'
              write(*,FMT6) dd_Betau(1,1,3), dd_Betau(1,2,3),
     &                    dd_Betau(1,3,3), dd_Betau(2,1,3),
     &                    dd_Betau(2,2,3), dd_Betau(2,3,3)
              write(*,FMT3) dd_Betau(3,1,3), dd_Betau(3,2,3),
     &                    dd_Betau(3,3,3)

              write(*,*)' d_gtd 11'
              write(*,FMT3) d_gtd(1,1,1), d_gtd(2,1,1), d_gtd(3,1,1) 
              write(*,*)' d_gtd 12'
              write(*,FMT3) d_gtd(1,1,2), d_gtd(2,1,2), d_gtd(3,1,2)
              write(*,*)' d_gtd 13'
              write(*,FMT3) d_gtd(1,1,3), d_gtd(2,1,3), d_gtd(3,1,3)
              write(*,*)' d_gtd 21'
              write(*,FMT3) d_gtd(1,2,1), d_gtd(2,2,1), d_gtd(3,2,1)
              write(*,*)' d_gtd 22'
              write(*,FMT3) d_gtd(1,2,2), d_gtd(2,2,2), d_gtd(3,2,2)
              write(*,*)' d_gtd 23'
              write(*,FMT3) d_gtd(1,2,3), d_gtd(2,2,3), d_gtd(3,2,3)
              write(*,*)' d_gtd 31'
              write(*,FMT3) d_gtd(1,3,1), d_gtd(2,3,1), d_gtd(3,3,1)
              write(*,*)' d_gtd 32'
              write(*,FMT3) d_gtd(1,3,2), d_gtd(2,3,2), d_gtd(3,3,2)
              write(*,*)' d_gtd 33'
              write(*,FMT3) d_gtd(1,3,3), d_gtd(2,3,3), d_gtd(3,3,3)

              write(*,*)' d_Atd 11'
              write(*,FMT3) d_Atd(1,1,1), d_Atd(2,1,1), d_Atd(3,1,1) 
              write(*,*)' d_Atd 12'
              write(*,FMT3) d_Atd(1,1,2), d_Atd(2,1,2), d_Atd(3,1,2)
              write(*,*)' d_Atd 13'
              write(*,FMT3) d_Atd(1,1,3), d_Atd(2,1,3), d_Atd(3,1,3)
              write(*,*)' d_Atd 21'
              write(*,FMT3) d_Atd(1,2,1), d_Atd(2,2,1), d_Atd(3,2,1)
              write(*,*)' d_Atd 22'
              write(*,FMT3) d_Atd(1,2,2), d_Atd(2,2,2), d_Atd(3,2,2)
              write(*,*)' d_Atd 23'
              write(*,FMT3) d_Atd(1,2,3), d_Atd(2,2,3), d_Atd(3,2,3)
              write(*,*)' d_Atd 31'
              write(*,FMT3) d_Atd(1,3,1), d_Atd(2,3,1), d_Atd(3,3,1)
              write(*,*)' d_Atd 32'
              write(*,FMT3) d_Atd(1,3,2), d_Atd(2,3,2), d_Atd(3,3,2)
              write(*,*)' d_Atd 33'
              write(*,FMT3) d_Atd(1,3,3), d_Atd(2,3,3), d_Atd(3,3,3)

              write(*,*)' dd_gtd 11'
              write(*,FMT6) dd_gtd(1,1,1,1), dd_gtd(1,2,1,1),
     &                    dd_gtd(1,3,1,1), dd_gtd(2,2,1,1),
     &                    dd_gtd(2,3,1,1), dd_gtd(3,3,1,1)

              write(*,*)' dd_gtd 12'
              write(*,FMT6) dd_gtd(1,1,1,2), dd_gtd(1,2,1,2),
     &                    dd_gtd(1,3,1,2), dd_gtd(2,2,1,2),
     &                    dd_gtd(2,3,1,2), dd_gtd(3,3,1,2)

              write(*,*)' dd_gtd 13'
              write(*,FMT6) dd_gtd(1,1,1,3), dd_gtd(1,2,1,3),
     &                    dd_gtd(1,3,1,3), dd_gtd(2,2,1,3),
     &                    dd_gtd(2,3,1,3), dd_gtd(3,3,1,3)

              write(*,*)' dd_gtd 22'
              write(*,FMT6) dd_gtd(1,1,2,2), dd_gtd(1,2,2,2),
     &                    dd_gtd(1,3,2,2), dd_gtd(2,2,2,2),
     &                    dd_gtd(2,3,2,2), dd_gtd(3,3,2,2)

              write(*,*)' dd_gtd 23'
              write(*,FMT6) dd_gtd(1,1,2,3), dd_gtd(1,2,2,3),
     &                    dd_gtd(1,3,2,3), dd_gtd(2,2,2,3),
     &                    dd_gtd(2,3,2,3), dd_gtd(3,3,2,3)

              write(*,*)' dd_gtd 33'
              write(*,FMT6) dd_gtd(1,1,3,3), dd_gtd(1,2,3,3),
     &                    dd_gtd(1,3,3,3), dd_gtd(2,2,3,3),
     &                    dd_gtd(2,3,3,3), dd_gtd(3,3,3,3)

              write(*,*)' Bu'
              write(*,FMT3) Bu(1), Bu(2), Bu(3)

            end if
          end if

          if (ltrace2) then
          if (i .eq. 46 .and. j .eq. 10 .and. k .eq. 60) then
            print *,'gtu(11)  = ',gtu(1,1)
            print *,'gtu(22)  = ',gtu(2,2)
            print *,'d_Gamt(11)   = ',d_Gamt(1,1)
            print *,'d_Gamt(21)   = ',d_Gamt(2,1)
            print *,'d_Gamt(22)   = ',d_Gamt(2,2)
            print *,'d_Gamt(33)   = ',d_Gamt(3,3)
            print *,'dd_gtd(1111) = ',dd_gtd(1,1,1,1)
            print *,'dd_gtd(1122) = ',dd_gtd(1,1,2,2)
            print *,'dd_gtd(3333) = ',dd_gtd(3,3,3,3)
            print *,'Ct(111)      = ',0.5d0*Ct(1,1,1)
            print *,'Ct(112)      = ',0.5d0*Ct(1,1,2)
            print *,'Ct(211)      = ',0.5d0*Ct(2,1,1)
            print *,'Ct(222)      = ',0.5d0*Ct(2,2,2)
            print *,'Ctd(222)     = ',0.5d0*Ctd(2,2,2)
            print *,'rho_AMD      = ',rho_ADM
            print *,'Atu(11)      = ',Atu(1,1)
            print *,'Atu(22)      = ',Atu(2,2)
            print *,'Atu(33)      = ',Atu(3,3)
            print *,'Atd_rhs(11)  = ',Atd_rhs(1,1)
            print *,'Atd_rhs(12)  = ',Atd_rhs(1,2)
            print *,'Atd_rhs(22)  = ',Atd_rhs(2,2)
            print *,'Atd_rhs(33)  = ',Atd_rhs(3,3)
            print *,'Gamt_rhs(1)  = ',Gamt_rhs(1)
            print *,'Gamt_rhs(2)  = ',Gamt_rhs(2)
            print *,'Gamt_rhs(3)  = ',Gamt_rhs(3)
!           do ii = 1, NU_G
!             print *,'RHS:  ',ii , dtu(ii)%d(i,j,k)
!           end do
!           stop
          end if
          end if


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
