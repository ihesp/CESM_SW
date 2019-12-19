module pbl_utils
!-----------------------------------------------------------------------!
! Module to hold PBL-related subprograms that may be used with multiple !
! different vertical diffusion schemes.                                 !
!                                                                       !
! Public subroutines:                                                   !
!     
!     calc_obklen                                                       !
!                                                                       !
!------------------ History --------------------------------------------!
! Created: Apr. 2012, by S. Santos                                      !
!-----------------------------------------------------------------------!

use shr_kind_mod, only: r8 => shr_kind_r8

implicit none
private

! Public Procedures
!----------------------------------------------------------------------!
! Excepting the initialization procedure, these are elemental
! procedures, so they can accept scalars or any dimension of array as
! arguments, as long as all arguments have the same number of
! elements.
public pbl_utils_init
public calc_ustar
public calc_obklen
public virtem
public compute_radf

real(r8), parameter :: ustar_min = 0.01_r8

real(r8) :: g         ! acceleration of gravity
real(r8) :: vk        ! Von Karman's constant
real(r8) :: cpair     ! specific heat of dry air
real(r8) :: rair      ! gas constant for dry air
real(r8) :: zvir      ! rh2o/rair - 1

contains

subroutine pbl_utils_init(g_in,vk_in,cpair_in,rair_in,zvir_in)

  !-----------------------------------------------------------------------!
  ! Purpose: Set constants to be used in calls to later functions         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: g_in       ! acceleration of gravity
  real(r8), intent(in) :: vk_in      ! Von Karman's constant
  real(r8), intent(in) :: cpair_in   ! specific heat of dry air
  real(r8), intent(in) :: rair_in    ! gas constant for dry air
  real(r8), intent(in) :: zvir_in    ! rh2o/rair - 1

  g = g_in
  vk = vk_in
  cpair = cpair_in
  rair = rair_in
  zvir = zvir_in

end subroutine pbl_utils_init

elemental subroutine calc_ustar( t,    pmid, taux, tauy, &
                                 rrho, ustar)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate ustar and bottom level density (necessary for      !
  !  Obukhov length calculation).                                         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t         ! surface temperature
  real(r8), intent(in) :: pmid      ! midpoint pressure (bottom level)
  real(r8), intent(in) :: taux      ! surface u stress [N/m2]
  real(r8), intent(in) :: tauy      ! surface v stress [N/m2]

  real(r8), intent(out) :: rrho     ! 1./bottom level density
  real(r8), intent(out) :: ustar    ! surface friction velocity [m/s]

  rrho = rair * t / pmid
  ustar = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )
  
end subroutine calc_ustar

elemental subroutine calc_obklen( ths,  thvs, qflx, shflx, rrho, ustar, &
                                  khfs, kqfs, kbfs, obklen)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate Obukhov length and kinematic fluxes.               !
  !-----------------------------------------------------------------------!

  real(r8), intent(in)  :: ths           ! potential temperature at surface [K]
  real(r8), intent(in)  :: thvs          ! virtual potential temperature at surface
  real(r8), intent(in)  :: qflx          ! water vapor flux (kg/m2/s)
  real(r8), intent(in)  :: shflx         ! surface heat flux (W/m2)

  real(r8), intent(in)  :: rrho          ! 1./bottom level density [ m3/kg ]
  real(r8), intent(in)  :: ustar         ! Surface friction velocity [ m/s ]
  
  real(r8), intent(out) :: khfs          ! sfc kinematic heat flux [mK/s]
  real(r8), intent(out) :: kqfs          ! sfc kinematic water vapor flux [m/s]
  real(r8), intent(out) :: kbfs          ! sfc kinematic buoyancy flux [m^2/s^3]
  real(r8), intent(out) :: obklen        ! Obukhov length
  
  ! Need kinematic fluxes for Obukhov:
  khfs = shflx*rrho/cpair
  kqfs = qflx*rrho
  kbfs = khfs + zvir*ths*kqfs
  
  ! Compute Obukhov length:
  obklen = -thvs * ustar**3 / (g*vk*(kbfs + sign(1.e-10_r8,kbfs)))

end subroutine calc_obklen

elemental real(r8) function virtem(t,q)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate virtual temperature from temperature and specific  !
  !  humidity.                                                            !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t, q

  virtem = t * (1.0_r8 + zvir*q)

end function virtem

subroutine compute_radf( choice_radf, i, pcols, pver, ncvmax, ncvfin, ktop, qmin, &
                         ql, pi, qrlw, g, cldeff, zi, chs, lwp_CL, opt_depth_CL,  &
                         radinvfrac_CL, radf_CL )
  ! -------------------------------------------------------------------------- !
  ! Purpose:                                                                   !
  ! Calculate cloud-top radiative cooling contribution to buoyancy production. !
  ! Here,  'radf' [m2/s3] is additional buoyancy flux at the CL top interface  !
  ! associated with cloud-top LW cooling being mainly concentrated near the CL !
  ! top interface ( just below CL top interface ).  Contribution of SW heating !
  ! within the cloud is not included in this radiative buoyancy production     !
  ! since SW heating is more broadly distributed throughout the CL top layer.  !
  ! -------------------------------------------------------------------------- !
  
  !-----------------!
  ! Input variables !
  !-----------------!
  character(len=6), intent(in) :: choice_radf  ! Method for calculating radf
  integer,  intent(in)  :: i                   ! Index of current column
  integer,  intent(in)  :: pcols               ! Number of atmospheric columns
  integer,  intent(in)  :: pver                ! Number of atmospheric layers
  integer,  intent(in)  :: ncvmax              ! Max numbers of CLs (perhaps equal to pver)
  integer,  intent(in)  :: ncvfin(pcols)       ! Total number of CL in column
  integer,  intent(in)  :: ktop(pcols, ncvmax) ! ktop for current column
  real(r8), intent(in)  :: qmin                ! Minimum grid-mean LWC counted as clouds [kg/kg]
  real(r8), intent(in)  :: ql(pcols, pver)     ! Liquid water specific humidity [ kg/kg ]
  real(r8), intent(in)  :: pi(pcols, pver+1)   ! Interface pressures [ Pa ]
  real(r8), intent(in)  :: qrlw(pcols, pver)   ! Input grid-mean LW heating rate : [ K/s ] * cpair * dp = [ W/kg*Pa ]
  real(r8), intent(in)  :: g                   ! Gravitational acceleration
  real(r8), intent(in)  :: cldeff(pcols,pver)  ! Effective Cloud Fraction [fraction]
  real(r8), intent(in)  :: zi(pcols, pver+1)   ! Interface heights [ m ]
  real(r8), intent(in)  :: chs(pcols, pver+1)  ! Buoyancy coeffi. saturated sl (heat) coef. at all interfaces.

  !------------------!
  ! Output variables !
  !------------------!
  real(r8), intent(out) :: lwp_CL(ncvmax)         ! LWP in the CL top layer [ kg/m2 ]
  real(r8), intent(out) :: opt_depth_CL(ncvmax)   ! Optical depth of the CL top layer
  real(r8), intent(out) :: radinvfrac_CL(ncvmax)  ! Fraction of LW radiative cooling confined in the top portion of CL
  real(r8), intent(out) :: radf_CL(ncvmax)        ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]

  !-----------------!
  ! Local variables !
  !-----------------!
  integer :: kt, ncv
  real(r8) :: lwp, opt_depth, radinvfrac, radf


  !-----------------!
  ! Begin main code !
  !-----------------!
  lwp_CL        = 0._r8
  opt_depth_CL  = 0._r8
  radinvfrac_CL = 0._r8
  radf_CL       = 0._r8

  ! ---------------------------------------- !
  ! Perform do loop for individual CL regime !
  ! ---------------------------------------- !
  do ncv = 1, ncvfin(i)
    kt = ktop(i,ncv)
    !-----------------------------------------------------!
    ! Compute radf for each CL regime and for each column !
    !-----------------------------------------------------!
    if( choice_radf .eq. 'orig' ) then
      if( ql(i,kt) .gt. qmin .and. ql(i,kt-1) .lt. qmin ) then 
        lwp       = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
        opt_depth = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        ! Approximate LW cooling fraction concentrated at the inversion by using
        ! polynomial approx to exact formula 1-2/opt_depth+2/(exp(opt_depth)-1))

        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
        radf        = qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) ! Cp*radiative cooling = [ W/kg ] 
        radf        = max( radinvfrac * radf * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 ) * chs(i,kt)
        ! We can disable cloud LW cooling contribution to turbulence by uncommenting:
        ! radf = 0._r8
      end if

    elseif( choice_radf .eq. 'ramp' ) then

      lwp         = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
      opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
      radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
      radinvfrac  = max(cldeff(i,kt)-cldeff(i,kt-1),0._r8) * radinvfrac 
      radf        = qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) ! Cp*radiative cooling [W/kg] 
      radf        = max( radinvfrac * radf * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 ) * chs(i,kt)

    elseif( choice_radf .eq. 'maxi' ) then

      ! Radiative flux divergence both in 'kt' and 'kt-1' layers are included 
      ! 1. From 'kt' layer
        lwp         = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
        opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
        radf        = max( radinvfrac * qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 )
      ! 2. From 'kt-1' layer and add the contribution from 'kt' layer
        lwp         = ql(i,kt-1) * ( pi(i,kt) - pi(i,kt-1) ) / g
        opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth) + opt_depth**2 )
        radf        = radf + max( radinvfrac * qrlw(i,kt-1) / ( pi(i,kt-1) - pi(i,kt) ) * ( zi(i,kt-1) - zi(i,kt) ), 0._r8 )
        radf        = max( radf, 0._r8 ) * chs(i,kt) 

    endif

    lwp_CL(ncv)        = lwp
    opt_depth_CL(ncv)  = opt_depth
    radinvfrac_CL(ncv) = radinvfrac
    radf_CL(ncv)       = radf 
  end do ! ncv = 1, ncvfin(i)
end subroutine compute_radf

end module pbl_utils
