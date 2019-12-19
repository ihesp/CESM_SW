module wv_saturation_cpe_xk
#define PARASIZE 12
!--------------------------------------------------------------------!
! Module Overview:                                                   !
!                                                                    !
! This module provides an interface to wv_sat_methods, providing     !
! saturation vapor pressure and related calculations to CAM.         !
!                                                                    !
! The original wv_saturation codes were introduced by J. J. Hack,    !
! February 1990. The code has been extensively rewritten since then, !
! including a total refactoring in Summer 2012.                      !
!                                                                    !
!--------------------------------------------------------------------!
! Methods:                                                           !
!                                                                    !
! Pure water/ice saturation vapor pressures are calculated on the    !
! fly, with the specific method determined by a runtime option.      !
! Mixed phase SVP is interpolated from the internal table, estbl,    !
! which is created during initialization.                            !
!                                                                    !
! The default method for calculating SVP is determined by a namelist !
! option, and used whenever svp_water/ice or qsat are called.        !
!                                                                    !
!--------------------------------------------------------------------!

!use shr_kind_mod, only: r8 => shr_kind_r8
!use physconst,    only: epsilo, &
!                        latvap, &
!                        latice, &
!                        rh2o,   &
!                        cpair,  &
!                        tmelt,  &
!                        h2otrip
!
use wv_sat_methods_cpe_xk, only: &
     svp_to_qsat => wv_sat_svp_to_qsat

implicit none
private
save

! Public interfaces
! Namelist, initialization, finalization
!public wv_sat_readnl
!public wv_sat_init
!public wv_sat_final

! Saturation vapor pressure calculations
public svp_water
public svp_ice
  
! Mixed phase (water + ice) saturation vapor pressure table lookup
public estblf

public svp_to_qsat

! Subroutines that return both SVP and humidity
! Optional arguments do temperature derivatives
public qsat           ! Mixed phase
public qsat_water     ! SVP over water only
public qsat_ice       ! SVP over ice only

! Wet bulb temperature solver
public :: findsp_vc, findsp

! Data
!xukai++++++++++++++++++++++++++++++++
!real(8) ::              epsilo, &
!                        latvap, &
!                        latice, &
!                        rh2o,   &
!                        cpair,  &
!                        tmelt,  &
!                        h2otrip

! 1. epsilo
! 2. latvap
! 3. latice 
! 4. rh2o
! 5. cpair
! 6. tmelt
! 7. h2otrip
! 8. omeps ! 1.0_r8 - epsilo
! 9. c3  !  c3 = 287.04d0*(7.5d0*log(10.d0))/cpair
!xukai--------------------------------

! This value is slightly high, but it seems to be the value for the
! steam point of water originally (and most frequently) used in the
! Goff & Gratch scheme.

! Table of saturation vapor pressure values (estbl) from tmin to
! tmax+1 Kelvin, in one degree increments.  ttrice defines the
! transition region, estbl contains a combination of ice & water
! values.
! Make these public parameters in case another module wants to see the
! extent of the table.

  !real(8), parameter :: tboil = 373.16d0
  !real(8), public, parameter :: tmin = 127.16d0
  !real(8), public, parameter :: tmax = 375.16d0

  !real(8), parameter :: ttrice = 20.00d0  ! transition range from es over H2O to es over ice

  !integer :: plenest                             ! length of estbl
  !real(8), allocatable :: estbl(:)              ! table values of saturation vapor pressure

  !real(8) :: omeps      ! 1.0_r8 - epsilo

  !real(8) :: c3         ! parameter used by findsp

  ! Set coefficients for polynomial approximation of difference
  ! between saturation vapor press over water and saturation pressure
  ! over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial is
  ! valid in the range -40 < t < 0 (degrees C).
!  real(8) :: pcf(5) = (/ &
!       5.04469588506d-01, &
!       -5.47288442819d+00, &
!       -3.67471858735d-01, &
!       -8.95963532403d-03, &
!       -7.78053686625d-05 /)
!
!   --- Degree 6 approximation ---
!  real(8) :: pcf(6) = (/ &
!       7.63285250063e-02, &
!       5.86048427932e+00, &
!       4.38660831780e-01, &
!       1.37898276415e-02, &
!       2.14444472424e-04, &
!       1.36639103771e-06 /)

contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

!subroutine wv_sat_readnl(nlfile)
!  !------------------------------------------------------------------!
!  ! Purpose:                                                         !
!  !   Get runtime options for wv_saturation.                         !
!  !------------------------------------------------------------------!
!
!  use wv_sat_methods, only: wv_sat_get_scheme_idx, &
!                            wv_sat_valid_idx, &
!                            wv_sat_set_default
!
!  use spmd_utils,      only: masterproc
!  use namelist_utils,  only: find_group_name
!  use units,           only: getunit, freeunit
!  use mpishorthand
!  use cam_abortutils,  only: endrun
!
!  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
!   
!  ! Local variables
!  integer :: unitn, ierr
!
!  character(len=32) :: wv_sat_scheme = "GoffGratch"
!
!  character(len=*), parameter :: subname = 'wv_sat_readnl'
!
!  namelist /wv_sat_nl/ wv_sat_scheme
!  !-----------------------------------------------------------------------------
!
!  if (masterproc) then
!     unitn = getunit()
!     open( unitn, file=trim(nlfile), status='old' )
!     call find_group_name(unitn, 'wv_sat_nl', status=ierr)
!     if (ierr == 0) then
!        read(unitn, wv_sat_nl, iostat=ierr)
!        if (ierr /= 0) then
!           call endrun(subname // ':: ERROR reading namelist')
!           return
!        end if
!     end if
!     close(unitn)
!     call freeunit(unitn)
!
!  end if
!
!#ifdef SPMD
!  call mpibcast(wv_sat_scheme, len(wv_sat_scheme) , mpichar, 0, mpicom)
!#endif
!
!  if (.not. wv_sat_set_default(wv_sat_scheme)) then
!     call endrun('wv_sat_readnl :: Invalid wv_sat_scheme.')
!     return
!  end if
!
!end subroutine wv_sat_readnl
!
!subroutine wv_sat_init
!  !------------------------------------------------------------------!
!  ! Purpose:                                                         !
!  !   Initialize module (e.g. setting parameters, initializing the   !
!  !   SVP lookup table).                                             !
!  !------------------------------------------------------------------!
!
!  use wv_sat_methods, only: wv_sat_methods_init, &
!                            wv_sat_get_scheme_idx, &
!                            wv_sat_valid_idx
!  use spmd_utils,     only: masterproc
!  use cam_logfile,    only: iulog
!  use cam_abortutils, only: endrun
!  use shr_assert_mod, only: shr_assert_in_domain
!  use error_messages, only: handle_errmsg
!
!  integer :: status
!
!  ! For wv_sat_methods error reporting.
!  character(len=256) :: errstring
!
!  ! For generating internal SVP table.
!  real(8) :: t         ! Temperature
!  integer  :: i         ! Increment counter
!
!  ! Precalculated because so frequently used.
!  omeps  = 1.0d0 - epsilo
!
!  ! Transition range method is only valid for transition temperatures at:
!  ! -40 deg C < T < 0 deg C
!  call shr_assert_in_domain(ttrice, ge=0.d0, le=40.d0, varname="ttrice",&
!       msg="wv_sat_init: Invalid transition temperature range.")
!
!! This parameter uses a hardcoded 287.04_r8?
!  c3 = 287.04d0*(7.5d0*log(10.d0))/cpair
!
!! Init "methods" module containing actual SVP formulae.
!
!  call wv_sat_methods_init(r8, tmelt, h2otrip, tboil, ttrice, &
!       epsilo, errstring)
!
!  call handle_errmsg(errstring, subname="wv_sat_methods_init")
!
!  ! Add two to make the table slightly too big, just in case.
!  plenest = ceiling(tmax-tmin) + 2
!
!  ! Allocate SVP table.
!  allocate(estbl(plenest), stat=status)
!  if (status /= 0) then 
!     call endrun('wv_sat_init :: ERROR allocating saturation vapor pressure table')
!     return
!  end if
!
!  do i = 1, plenest
!    estbl(i) = svp_trans(tmin + real(i-1,r8))
!  end do
!
!  if (masterproc) then
!     write(iulog,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
!  end if
!
!end subroutine wv_sat_init
!
!subroutine wv_sat_final
!  !------------------------------------------------------------------!
!  ! Purpose:                                                         !
!  !   Deallocate global variables in module.                         !
!  !------------------------------------------------------------------!
!  use cam_abortutils, only: endrun
!
!  integer :: status
!
!  if (allocated(estbl)) then
!
!     deallocate(estbl, stat=status)
!
!     if (status /= 0) then
!        call endrun('wv_sat_final :: ERROR deallocating table')
!        return
!     end if
!
!  end if
!
!end subroutine wv_sat_final

!---------------------------------------------------------------------
! DEFAULT SVP FUNCTIONS
!---------------------------------------------------------------------

! Compute saturation vapor pressure over water
 function svp_water(t, wv_para) result(es)

  use wv_sat_methods_cpe_xk, only: &
       wv_sat_svp_water
   implicit none

  real(8), intent(in) :: t ! Temperature (K)
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: es            ! SVP (Pa)

  es = wv_sat_svp_water(T, wv_para)

end function svp_water

! Compute saturation vapor pressure over ice
 function svp_ice(t, wv_para) result(es)

  use wv_sat_methods_cpe_xk, only: &
       wv_sat_svp_ice

   implicit none
  real(8), intent(in) :: t ! Temperature (K)
  real(8) :: es            ! SVP (Pa)
  real(8), intent(in) :: wv_para(PARASIZE)

  es = wv_sat_svp_ice(T, wv_para)

end function svp_ice

! Compute saturation vapor pressure with an ice-water transition
 function svp_trans(t, wv_para) result(es)

  use wv_sat_methods_cpe_xk, only: &
       wv_sat_svp_trans

   implicit none
  real(8), intent(in) :: t ! Temperature (K)
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: es            ! SVP (Pa)

  es = wv_sat_svp_trans(T, wv_para)

end function svp_trans

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Does linear interpolation from nearest values found
! in the table (estbl).
 function estblf(t, estbl) result(es)

   implicit none
  real(8), intent(in) :: t ! Temperature 
  real(8) :: es            ! SVP (Pa)

  integer  :: i         ! Index for t in the table
  real(8) :: t_tmp     ! intermediate temperature for es look-up

  real(8) :: weight ! Weight for interpolation
  real(8),  intent(in) :: estbl(:)
  real(8),  parameter :: tmin = 127.16d0
  real(8),  parameter :: tmax = 375.16d0

  t_tmp = max(min(t,tmax)-tmin, 0.d0)   ! Number of table entries above tmin
  i = int(t_tmp) + 1                     ! Corresponding index.
  weight = t_tmp - aint(t_tmp, 8)       ! Fractional part of t_tmp (for interpolation).
  es = (1.d0 - weight)*estbl(i) + weight*estbl(i+1)

end function estblf

! Get enthalpy based only on temperature
! and specific humidity.
 function tq_enthalpy(t, q, hltalt, wv_para) result(enthalpy)

   implicit none
  real(8), intent(in) :: t      ! Temperature
  real(8), intent(in) :: q      ! Specific humidity
  real(8), intent(in) :: hltalt ! Modified hlat for T derivatives

  real(8) :: enthalpy
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: cpair 
  cpair = wv_para(5)

  enthalpy = cpair * t + hltalt * q
  
end function tq_enthalpy

!---------------------------------------------------------------------
! LATENT HEAT OF VAPORIZATION CORRECTIONS
!---------------------------------------------------------------------

 subroutine no_ip_hltalt(t, hltalt, wv_para)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate latent heat of vaporization of pure liquid water at  !
  !   a given temperature.                                           !
  !------------------------------------------------------------------!

   implicit none
  ! Inputs
  real(8), intent(in) :: t        ! Temperature
  ! Outputs
  real(8), intent(out) :: hltalt  ! Appropriately modified hlat
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: latvap
  real(8) :: tmelt
  latvap = wv_para(2)
  tmelt = wv_para(6)

  hltalt = latvap

  ! Account for change of latvap with t above freezing where
  ! constant slope is given by -2369 j/(kg c) = cpv - cw
  if (t >= tmelt) then
     hltalt = hltalt - 2369.0d0*(t-tmelt)
  end if

end subroutine no_ip_hltalt

 subroutine calc_hltalt(t, hltalt, wv_para, tterm)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate latent heat of vaporization of water at a given      !
  !   temperature, taking into account the ice phase if temperature  !
  !   is below freezing.                                             !
  !   Optional argument also calculates a term used to calculate     !
  !   d(es)/dT within the water-ice transition range.                !
  !------------------------------------------------------------------!

   implicit none
  ! Inputs
  real(8), intent(in) :: t        ! Temperature
  ! Outputs
  real(8), intent(out) :: hltalt  ! Appropriately modified hlat
  ! Term to account for d(es)/dT in transition region.
  real(8), intent(out), optional :: tterm

  ! Local variables
  real(8) :: tc      ! Temperature in degrees C
  real(8) :: weight  ! Weight for es transition from water to ice
  ! Loop iterator
  integer :: i
  real(8), parameter :: ttrice = 20.00d0  ! transition range from es over H2O to es over ice
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: tmelt
  real(8) :: latice
  real(8) :: pcf(5) 
  DATA pcf / &
       5.04469588506d-01, &
       -5.47288442819d+00, &
       -3.67471858735d-01, &
       -8.95963532403d-03, &
       -7.78053686625d-05 /

  tmelt = wv_para(6) 
  latice = wv_para(3)

  if (present(tterm)) tterm = 0.0d0

  call no_ip_hltalt(t,hltalt, wv_para)
  if (t < tmelt) then
     ! Weighting of hlat accounts for transition from water to ice.
     tc = t - tmelt

     if (tc >= -ttrice) then
        weight = -tc/ttrice

        ! polynomial expression approximates difference between es
        ! over water and es over ice from 0 to -ttrice (C) (max of
        ! ttrice is 40): required for accurate estimate of es
        ! derivative in transition range from ice to water
        if (present(tterm)) then
           do i = size(pcf), 1, -1
              tterm = pcf(i) + tc*tterm
           end do
           tterm = tterm/ttrice
        end if

     else
        weight = 1.0d0
     end if

     hltalt = hltalt + weight*latice

  end if

end subroutine calc_hltalt

!---------------------------------------------------------------------
! OPTIONAL OUTPUTS
!---------------------------------------------------------------------

! Temperature derivative outputs, for qsat_*
 subroutine deriv_outputs(t, p, es, qs, hltalt, tterm, wv_para,&
     gam, dqsdt)

   implicit none
  ! Inputs
  real(8), intent(in) :: t      ! Temperature
  real(8), intent(in) :: p      ! Pressure
  real(8), intent(in) :: es     ! Saturation vapor pressure
  real(8), intent(in) :: qs     ! Saturation specific humidity
  real(8), intent(in) :: hltalt ! Modified latent heat
  real(8), intent(in) :: tterm  ! Extra term for d(es)/dT in
                                 ! transition region.
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: omeps
  real(8) :: rh2o 
  real(8) :: cpair 

  ! Outputs
  real(8), intent(out), optional :: gam      ! (hltalt/cpair)*(d(qs)/dt)
  real(8), intent(out), optional :: dqsdt    ! (d(qs)/dt)

  ! Local variables
  real(8) :: desdt        ! d(es)/dt
  real(8) :: dqsdt_loc    ! local copy of dqsdt
  rh2o = wv_para(4)
  omeps = wv_para(8)
  cpair = wv_para(5)

  if (qs == 1.0d0) then
     dqsdt_loc = 0.d0
  else
     desdt = hltalt*es/(rh2o*t*t) + tterm
     dqsdt_loc = qs*p*desdt/(es*(p-omeps*es))
  end if

  if (present(dqsdt)) dqsdt = dqsdt_loc
  if (present(gam))   gam   = dqsdt_loc * (hltalt/cpair)

end subroutine deriv_outputs

!---------------------------------------------------------------------
! QSAT (SPECIFIC HUMIDITY) PROCEDURES
!---------------------------------------------------------------------

 subroutine qsat(t, p, es, qs, estbl, wv_para, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Look up and return saturation vapor pressure from precomputed  !
  !   table, then calculate and return saturation specific humidity. !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

   implicit none
  ! Inputs
  real(8), intent(in) :: t    ! Temperature
  real(8), intent(in) :: p    ! Pressure
  ! Outputs
  real(8), intent(out) :: es  ! Saturation vapor pressure
  real(8), intent(out) :: qs  ! Saturation specific humidity
  real(8), intent(in) :: estbl(:)
  real(8), intent(in) :: wv_para(PARASIZE)

  real(8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(8) :: hltalt       ! Modified latent heat for T derivatives
  real(8) :: tterm        ! Account for d(es)/dT in transition region

  es = estblf(t, estbl)

  qs = svp_to_qsat(es, p, wv_para)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

  ! Calculate optional arguments.
  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! "generalized" analytic expression for t derivative of es
     ! accurate to within 1 percent for 173.16 < t < 373.16
     call calc_hltalt(t, hltalt, wv_para, tterm)

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt, wv_para)

     call deriv_outputs(t, p, es, qs, hltalt, tterm, wv_para, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat

 subroutine qsat_water(t, p, es, qs, wv_para, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  use wv_sat_methods_cpe_xk, only: wv_sat_qsat_water

   implicit none
  ! Inputs
  real(8), intent(in) :: t    ! Temperature
  real(8), intent(in) :: p    ! Pressure
  ! Outputs
  real(8), intent(out) :: es  ! Saturation vapor pressure
  real(8), intent(out) :: qs  ! Saturation specific humidity
  real(8), intent(in) :: wv_para(PARASIZE)

  real(8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(8) :: hltalt       ! Modified latent heat for T derivatives

  !call put_tmp1(2, t)
  !call put_tmp2(2, p) !tmp2(tmpind(1)) = p

        !call lwpf_start_test_p6() 
  call wv_sat_qsat_water(t, p, es, qs, wv_para)
  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! "generalized" analytic expression for t derivative of es
     ! accurate to within 1 percent for 173.16 < t < 373.16
     call no_ip_hltalt(t, hltalt, wv_para)

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt, wv_para)

     ! For pure water/ice transition term is 0.
     call deriv_outputs(t, p, es, qs, hltalt, 0.d0, wv_para, &
          gam=gam, dqsdt=dqsdt)

  end if

        !call lwpf_stop_test_p6() 
end subroutine qsat_water

 subroutine qsat_ice(t, p, es, qs, wv_para, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  use wv_sat_methods_cpe_xk, only: wv_sat_qsat_ice
   implicit none

  ! Inputs
  real(8), intent(in) :: t    ! Temperature
  real(8), intent(in) :: p    ! Pressure
  ! Outputs
  real(8), intent(out) :: es  ! Saturation vapor pressure
  real(8), intent(out) :: qs  ! Saturation specific humidity
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: latvap, latice

  real(8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(8) :: hltalt       ! Modified latent heat for T derivatives
  latvap = wv_para(2) 
  latice = wv_para(3)

  call wv_sat_qsat_ice(t, p, es, qs, wv_para)

  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! For pure ice, just add latent heats.
     hltalt = latvap + latice

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt, wv_para)

     ! For pure water/ice transition term is 0.
     call deriv_outputs(t, p, es, qs, hltalt, 0.d0, wv_para, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat_ice

!---------------------------------------------------------------------
! FINDSP (WET BULB TEMPERATURE) PROCEDURES
!---------------------------------------------------------------------

subroutine findsp_vc(q, t, p, use_ice, tsp, qsp, estbl, wv_para)

  !use cam_logfile,      only: iulog
  !use cam_abortutils,   only: endrun

  ! Wrapper for findsp which is 1D and handles the output status.
  ! Changing findsp to elemental restricted debugging output.
  ! If that output is needed again, it's preferable *not* to copy findsp,
  ! but to change the existing version.

   implicit none
  ! input arguments
!  real(8), intent(in) :: q(:)        ! water vapor (kg/kg)
!  real(8), intent(in) :: t(:)        ! temperature (K)
!  real(8), intent(in) :: p(:)        ! pressure    (Pa)
!  logical,  intent(in) :: use_ice     ! flag to include ice phase in calculations
!
!  ! output arguments
!  real(8), intent(out) :: tsp(:)     ! saturation temp (K)
!  real(8), intent(out) :: qsp(:)     ! saturation mixing ratio (kg/kg)
  real(8), intent(in) :: q        ! water vapor (kg/kg)
  real(8), intent(in) :: t        ! temperature (K)
  real(8), intent(in) :: p        ! pressure    (Pa)
  logical,  intent(in) :: use_ice     ! flag to include ice phase in calculations

  ! output arguments
  real(8), intent(out) :: tsp     ! saturation temp (K)
  real(8), intent(out) :: qsp     ! saturation mixing ratio (kg/kg)
 
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8), intent(in) :: estbl(:)

  !integer :: status(size(q))   ! flag representing state of output
  integer :: status            ! flag representing state of output
                               ! 0 => Successful convergence
                               ! 1 => No calculation done: pressure or specific
                               !      humidity not within usable range
                               ! 2 => Run failed to converge
                               ! 4 => Temperature fell below minimum
                               ! 8 => Enthalpy not conserved

  integer :: n, i

  !n = size(q)

        !call lwpf_start_test_p6() 
  call findsp(q, t, p, use_ice, tsp, qsp, status, estbl, wv_para)
        !call lwpf_stop_test_p6() 
  ! Currently, only 2 and 8 seem to be treated as fatal errors.
  !do i = 1,n
  !   if (status(i) == 2) then
  !      write(iulog,*) ' findsp not converging at i = ', i
  !      write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
  !      write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
  !      call endrun ('wv_saturation::FINDSP -- not converging')
  !   else if (status(i) == 8) then
  !      write(iulog,*) ' the enthalpy is not conserved at i = ', i
  !      write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
  !      write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
  !      call endrun ('wv_saturation::FINDSP -- enthalpy is not conserved')
  !   endif
  !end do

end subroutine findsp_vc

 subroutine findsp (q, t, p, use_ice, tsp, qsp, status, estbl, wv_para)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!     find the wet bulb temperature for a given t and q
!     in a longitude height section
!     wet bulb temp is the temperature and spec humidity that is 
!     just saturated and has the same enthalpy
!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method: 
! a Newton method is used
! first guess uses an algorithm provided by John Petch from the UKMO
! we exclude points where the physical situation is unrealistic
! e.g. where the temperature is outside the range of validity for the
!      saturation vapor pressure, or where the water vapor pressure
!      exceeds the ambient pressure, or the saturation specific humidity is 
!      unrealistic
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
!
!     input arguments
!

   implicit none
  real(8), intent(in) :: q        ! water vapor (kg/kg)
  real(8), intent(in) :: t        ! temperature (K)
  real(8), intent(in) :: p        ! pressure    (Pa)
  logical,  intent(in) :: use_ice  ! flag to include ice phase in calculations
!
! output arguments
!
  real(8), intent(out) :: tsp      ! saturation temp (K)
  real(8), intent(out) :: qsp      ! saturation mixing ratio (kg/kg)
  integer,  intent(out) :: status   ! flag representing state of output
                                    ! 0 => Successful convergence
                                    ! 1 => No calculation done: pressure or specific
                                    !      humidity not within usable range
                                    ! 2 => Run failed to converge
                                    ! 4 => Temperature fell below minimum
                                    ! 8 => Enthalpy not conserved
    real(8), intent(in) :: estbl(:) 
    real(8), intent(in) :: wv_para(PARASIZE)
    real(8) :: cpair
    real(8) :: c3
!
! local variables
!
  integer, parameter :: iter = 8    ! max number of times to iterate the calculation
  integer :: l                      ! iterator

  real(8) es                   ! sat. vapor pressure
  real(8) gam                  ! change in sat spec. hum. wrt temperature (times hltalt/cpair)
  real(8) dgdt                 ! work variable
  real(8) g                    ! work variable
  real(8) hltalt               ! lat. heat. of vap.
  real(8) qs                   ! spec. hum. of water vapor

! work variables
  real(8) t1, q1, dt, dq
  real(8) qvd
  real(8) r1b, c1, c2

!  real(8), parameter :: dttol = 1.d-4 ! the relative temp error tolerance required to quit the iteration
!  real(8), parameter :: dqtol = 1.d-4 ! the relative moisture error tolerance required to quit the iteration
!  real(8) enin, enout
!  real(8),  parameter :: tmin = 127.16d0
!  real(8),  parameter :: tmax = 375.16d0

  real(8) :: dttol ! the relative temp error tolerance required to quit the iteration
  real(8) :: dqtol ! the relative moisture error tolerance required to quit the iteration
  real(8) enin, enout
  real(8) :: tmin
  real(8) :: tmax

  !      call lwpf_start_test_p6() 
  dttol = 1.d-4
  dqtol = 1.d-4
  tmin = 127.16d0
  tmax = 375.16d0
  cpair = wv_para(5)
  c3 = wv_para(9)
  ! Saturation specific humidity at this temperature
  if (use_ice) then
     call qsat(t, p, es, qs, estbl, wv_para)
  else
     call qsat_water(t, p, es, qs, wv_para)
  end if

  !      call lwpf_stop_test_p6() 
  ! make sure a meaningful calculation is possible
  if (p <= 5.d0*es .or. qs <= 0.d0 .or. qs >= 0.5d0 &
       .or. t < tmin .or. t > tmax) then
     status = 1
     ! Keep initial parameters when conditions aren't suitable
     tsp = t
     qsp = q
     enin = 1.d0
     enout = 1.d0

     return
  end if

  ! Prepare to iterate
  status = 2

  ! Get initial enthalpy
  if (use_ice) then
     call calc_hltalt(t,hltalt, wv_para)
  else
     call no_ip_hltalt(t,hltalt, wv_para)
  end if
  enin = tq_enthalpy(t, q, hltalt, wv_para)

  ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
  c1 = hltalt*c3
  c2 = (t + 36.d0)**2
  r1b = c2/(c2 + c1*qs)
  qvd = r1b * (q - qs)
  tsp = t + ((hltalt/cpair)*qvd)

  ! Generate qsp, gam, and enout from tsp.
  if (use_ice) then
     call qsat(tsp, p, es, qsp, estbl, wv_para, gam=gam, enthalpy=enout)
  else
     call qsat_water(tsp, p, es, qsp, wv_para, gam=gam, enthalpy=enout)
  end if

  ! iterate on first guess
  do l = 1, iter

     g = enin - enout
     dgdt = -cpair * (1 + gam)

     ! New tsp
     t1 = tsp - g/dgdt
     dt = abs(t1 - tsp)/t1
     tsp = t1

     ! bail out if past end of temperature range
     if ( tsp < tmin ) then
        tsp = tmin
        ! Get latent heat and set qsp to a value
        ! that preserves enthalpy.
        if (use_ice) then
           call calc_hltalt(tsp,hltalt, wv_para)
        else
           call no_ip_hltalt(tsp,hltalt, wv_para)
        end if
        qsp = (enin - cpair*tsp)/hltalt
        enout = tq_enthalpy(tsp, qsp, hltalt, wv_para)
        status = 4
        exit
     end if

     ! Re-generate qsp, gam, and enout from new tsp.
     if (use_ice) then
        call qsat(tsp, p, es, q1, estbl, wv_para, gam=gam, enthalpy=enout)
     else
        call qsat_water(tsp, p, es, q1, wv_para, gam=gam, enthalpy=enout)
     end if
     dq = abs(q1 - qsp)/max(q1,1.d-12)
     qsp = q1

     ! if converged at this point, exclude it from more iterations
     if (dt < dttol .and. dq < dqtol) then
        status = 0
        exit
     endif
  end do

  ! Test for enthalpy conservation
  if (abs((enin-enout)/(enin+enout)) > 1.d-4) status = 8

end subroutine findsp

end module wv_saturation_cpe_xk
