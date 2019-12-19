module wv_sat_methods_cpe_xk

#define PARASIZE 12
! This portable module contains all CAM methods for estimating
! the saturation vapor pressure of water.
!
! wv_saturation provides CAM-specific interfaces and utilities
! based on these formulae.
!
! Typical usage of this module:
!
! Init:
! call wv_sat_methods_init(r8, <constants>, errstring)
!
! Get scheme index from a name string:
! scheme_idx = wv_sat_get_scheme_idx(scheme_name)
! if (.not. wv_sat_valid_idx(scheme_idx)) <throw some error>
!
! Get pressures:
! es = wv_sat_svp_water(t, scheme_idx)
! es = wv_sat_svp_ice(t, scheme_idx)
!
! Use ice/water transition range:
! es = wv_sat_svp_trice(t, ttrice, scheme_idx)
!
! Note that elemental functions cannot be pointed to, nor passed
! as arguments. If you need to do either, it is recommended to
! wrap the function so that it can be given an explicit (non-
! elemental) interface.

implicit none
private
save
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
! 8. omeps
! 9. c3  !  c3 = 287.04d0*(7.5d0*log_c(10.d0))/cpair
!xukai--------------------------------


!integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real

!real(8) :: tmelt   ! Melting point of water at 1 atm (K)
!real(8) :: h2otrip ! Triple point temperature of water (K)
!real(8) :: tboil   ! Boiling point of water at 1 atm (K)
!
!real(8) :: ttrice  ! Ice-water transition range
!
!real(8) :: epsilo  ! Ice-water transition range
!real(8) :: omeps   ! 1._r8 - epsilo

! Indices representing individual schemes

!public wv_sat_methods_init
!public wv_sat_get_scheme_idx
!public wv_sat_valid_idx
!
!public wv_sat_set_default
!public wv_sat_reset_default

public wv_sat_svp_water
public wv_sat_svp_ice
public wv_sat_svp_trans

! pressure -> humidity conversion
public wv_sat_svp_to_qsat

! Combined qsat operations
public wv_sat_qsat_water
public wv_sat_qsat_ice
public wv_sat_qsat_trans
!public tmp1 
!public tmp2 
!public tmpind
!public wvflag  
!
!    real(8), pointer :: tmp1(:)
!    real(8), pointer :: tmp2(:)
!    integer, pointer :: tmpind(:)
!    integer :: wvflag 
 
interface 
    real(8) function log10_c(in1)
    real(8), intent(in) :: in1
    end function log10_c

    real(8) function log_c(in1)
    real(8), intent(in) :: in1
    end function log_c

    real(8) function pow_c(in1, in2)
    real(8), intent(in) :: in1, in2
    end function pow_c
   

end interface

!interface operator(**)
!              module procedure pow2
!end interface


contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

! Get physical constants
!subroutine wv_sat_methods_init(kind, tmelt_in, h2otrip_in, tboil_in, &
!     ttrice_in, epsilo_in, errstring)
!  integer, intent(in) :: kind
!  real(8), intent(in) :: tmelt_in
!  real(8), intent(in) :: h2otrip_in
!  real(8), intent(in) :: tboil_in
!  real(8), intent(in) :: ttrice_in
!  real(8), intent(in) :: epsilo_in
!  character(len=*), intent(out)  :: errstring
!
!  errstring = ' '
!
!  if (kind /= r8) then
!     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
!          kind,' was input kind but ',r8,' is internal kind.'
!     return
!  end if
!
!  if (ttrice_in < 0.d0) then
!     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
!          ttrice_in,' was input for ttrice, but negative range is invalid.'
!     return
!  end if
!
!  tmelt = tmelt_in
!  h2otrip = h2otrip_in
!  tboil = tboil_in
!  ttrice = ttrice_in
!  epsilo = epsilo_in
!
!  omeps = 1.d0 - epsilo
!
!end subroutine wv_sat_methods_init
!real(8) function pow2(in1, in2)
!real(8), intent(in) :: in1, in2 
!real(8) :: out1
!    call math_agent_pow_c(in1, in2, out1)
!    pow2 = out1
!end function
!
   

! Look up index by name.
!pure function wv_sat_get_scheme_idx(name) result(idx)
!  implicit none
!  character(len=*), intent(in) :: name
!  integer :: idx
!  
!  select case (name)
!  case("GoffGratch")
!     idx = GoffGratch_idx
!  case("MurphyKoop")
!     idx = MurphyKoop_idx
!  case("OldGoffGratch")
!     idx = OldGoffGratch_idx
!  case("Bolton")
!     idx = Bolton_idx
!  case default
!     idx = Invalid_idx
!  end select
!
!end function wv_sat_get_scheme_idx

! Check validity of an index from the above routine.
!pure function wv_sat_valid_idx(idx) result(status)
!  implicit none
!  integer, intent(in) :: idx
!  logical :: status
!
!  status = (idx /= Invalid_idx)
!
!end function wv_sat_valid_idx

! Set default scheme (otherwise, Goff & Gratch is default)
! Returns a logical representing success (.true.) or
! failure (.false.).
!function wv_sat_set_default(name) result(status)
!  implicit none
!  character(len=*), intent(in) :: name
!  logical :: status
!
!  ! Don't want to overwrite valid default with invalid,
!  ! so assign to temporary and check it first.
!  integer :: tmp_idx
!
!  tmp_idx = wv_sat_get_scheme_idx(name)
!
!  status = wv_sat_valid_idx(tmp_idx)
!
!  if (status) default_idx = tmp_idx
!
!end function wv_sat_set_default

! Reset default scheme to initial value.
! The same thing can be accomplished with wv_sat_set_default;
! the real reason to provide this routine is to reset the
! module for testing purposes.
!subroutine wv_sat_reset_default()
!  implicit none
!
!  default_idx = initial_default_idx
!
!end subroutine wv_sat_reset_default

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Get saturation specific humidity given pressure and SVP.
! Specific humidity is limited to range 0-1.
 function wv_sat_svp_to_qsat(es, p, wv_para) result(qs)

  implicit none
  real(8), intent(in) :: es  ! SVP
  real(8), intent(in) :: p   ! Current pressure.
  real(8), intent(in) :: wv_para(PARASIZE) 
  real(8) :: epsilo, omeps
  real(8) :: qs

   epsilo = wv_para(1)
   omeps = wv_para(8)
  ! If pressure is less than SVP, set qs to maximum of 1.
  if ( (p - es) <= 0.d0 ) then
     qs = 1.0d0
  else
     qs = epsilo*es / (p - omeps*es)
  end if

end function wv_sat_svp_to_qsat

 subroutine wv_sat_qsat_water(t, p, es, qs, wv_para, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!
  implicit none
  ! Inputs
  real(8), intent(in) :: t    ! Temperature
  real(8), intent(in) :: p    ! Pressure
  ! Outputs
  real(8), intent(out) :: es  ! Saturation vapor pressure
  real(8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
  real(8), intent(in) :: wv_para(PARASIZE)

  es = wv_sat_svp_water(t, wv_para, idx)
  !call put_tmp1(1, es)

  qs = wv_sat_svp_to_qsat(es, p, wv_para)
  !call put_tmp2(1, qs)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_water

 subroutine wv_sat_qsat_ice(t, p, es, qs, wv_para, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!
 implicit none
  ! Inputs
  real(8), intent(in) :: t    ! Temperature
  real(8), intent(in) :: p    ! Pressure
  ! Outputs
  real(8), intent(out) :: es  ! Saturation vapor pressure
  real(8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
  real(8), intent(in) :: wv_para(PARASIZE) 

  es = wv_sat_svp_ice(t, wv_para, idx)

  qs = wv_sat_svp_to_qsat(es, p, wv_para)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_ice

 subroutine wv_sat_qsat_trans(t, p, es, qs, wv_para, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  implicit none
  ! Inputs
  real(8), intent(in) :: t    ! Temperature
  real(8), intent(in) :: p    ! Pressure
  ! Outputs
  real(8), intent(out) :: es  ! Saturation vapor pressure
  real(8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
  real(8), intent(in) :: wv_para(PARASIZE) 

  es = wv_sat_svp_trans(t, wv_para, idx)

  qs = wv_sat_svp_to_qsat(es, p, wv_para)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_trans

!---------------------------------------------------------------------
! SVP INTERFACE FUNCTIONS
!---------------------------------------------------------------------

 function wv_sat_svp_water(t, wv_para, idx) result(es)
  implicit none
  real(8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(8) :: es
  real(8), intent(in) :: wv_para(PARASIZE)

  integer :: use_idx
integer :: Invalid_idx 
integer :: OldGoffGratch_idx 
integer :: GoffGratch_idx
integer :: MurphyKoop_idx 
integer :: Bolton_idx 

! Index representing the current default scheme.
!integer, parameter :: initial_default_idx = GoffGratch_idx
integer :: default_idx 
!call lwpf_start_test_p1()
    default_idx = 1 

   Invalid_idx = -1
   OldGoffGratch_idx = 0
   GoffGratch_idx = 1
   MurphyKoop_idx = 2
   Bolton_idx = 3



  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(1)
     es = GoffGratch_svp_water(t)
     !write(*,*) "GoffGratch_idx"
  case(2)
     es = MurphyKoop_svp_water(t)
     !write(*,*) "MurphyKoop_idx"
  case(0)
     es = OldGoffGratch_svp_water(t)
     !write(*,*) "OldGoffGratch_idx"
  case(3)
     ! write(*,*) "Bolton_idx" 
     es = Bolton_svp_water(t, wv_para)
  end select

!call lwpf_stop_test_p1()
end function wv_sat_svp_water

 function wv_sat_svp_ice(t, wv_para, idx) result(es)
  implicit none
  real(8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: es

  integer :: use_idx

integer :: Invalid_idx 
integer :: OldGoffGratch_idx 
integer :: GoffGratch_idx
integer :: MurphyKoop_idx 
integer :: Bolton_idx 

! Index representing the current default scheme.
!integer, parameter :: initial_default_idx = GoffGratch_idx
integer :: default_idx 
    default_idx = 1 

   Invalid_idx = -1
   OldGoffGratch_idx = 0
   GoffGratch_idx = 1
   MurphyKoop_idx = 2
   Bolton_idx = 3


  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(1)
     es = GoffGratch_svp_ice(t, wv_para)
  case(2)
     es = MurphyKoop_svp_ice(t)
  case(0)
     es = OldGoffGratch_svp_ice(t, wv_para)
  case(3)
     es = Bolton_svp_water(t, wv_para)
  end select

end function wv_sat_svp_ice

 function wv_sat_svp_trans(t, wv_para, idx) result (es)
 implicit none

  real(8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(8), intent(in) :: wv_para(PARASIZE)

  real(8) :: es

  real(8) :: esice      ! Saturation vapor pressure over ice
  real(8) :: weight     ! Intermediate scratch variable for es transition
  !real(8), parameter :: ttrice = 20.00d0  ! transition range from es over H2O to es over ice
  real(8) :: tmelt

  real(8) :: ttrice  ! transition range from es over H2O to es over ice

  ttrice = 20.00d0 
!
! Water
!
  tmelt = wv_para(6)
  if (t >= (tmelt - ttrice)) then
     es = wv_sat_svp_water(t, wv_para, idx)
  else
     es = 0.0d0
  end if

!
! Ice
!
  if (t < tmelt) then

     esice = wv_sat_svp_ice(t, wv_para, idx)

     if ( (tmelt - t) > ttrice ) then
        weight = 1.0d0
     else
        weight = (tmelt - t)/ttrice
     end if

     es = weight*esice + (1.0d0 - weight)*es
  end if

end function wv_sat_svp_trans

!---------------------------------------------------------------------
! SVP METHODS
!---------------------------------------------------------------------

! Goff & Gratch (1946)

 function GoffGratch_svp_water(t) result(es)
  implicit none
  real(8), intent(in) :: t  ! Temperature in Kelvin
  real(8) :: es             ! SVP in Pa
  !real(8), parameter :: tboil = 373.16d0
  real(8) :: tboil 
  real(8) :: pow_tmp
  real(8) :: base_tmp
  real(8) :: log_tmp

  tboil = 373.16d0

!pow_tmp = (-7.90298d0*(tboil/t-1.d0)+ &
!       5.02808d0*log10_c(tboil/t)- &
!       1.3816d-7*(pow_c(10.d0, (11.344d0*(1.d0-t/tboil)))-1.d0)+ &
!       8.1328d-3*(pow_c(10.d0, (-3.49149d0*(tboil/t-1.d0)))-1.d0)+ &
!       log10_c(1013.246d0))

  ! uncertain below -70 C
  !es = 10.d0
  base_tmp = 10.d0
  log_tmp = 1013.246d0
  es = pow_c( base_tmp, (-7.90298d0*(tboil/t-1.d0)+ &
       5.02808d0*log10_c(tboil/t)- &
       1.3816d-7*(pow_c(base_tmp, (11.344d0*(1.d0-t/tboil)))-1.d0)+ &
       8.1328d-3*(pow_c(base_tmp, (-3.49149d0*(tboil/t-1.d0)))-1.d0)+ &
       log10_c(log_tmp))) * 100.d0
!  es = 10._r8**(-7.90298_r8*(tboil/t-1._r8)+ &
!       5.02808_r8*log10(tboil/t)- &
!       1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/tboil))-1._r8)+ &
!       8.1328e-3_r8*(10._r8**(-3.49149_r8*(tboil/t-1._r8))-1._r8)+ &
!       log10(1013.246_r8))*100._r8



end function GoffGratch_svp_water

 function GoffGratch_svp_ice(t, wv_para) result(es)
  implicit none
  real(8), intent(in) :: t  ! Temperature in Kelvin
  real(8) :: es             ! SVP in Pa
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: h2otrip
  real(8) :: base_tmp
  real(8) :: log_tmp
  h2otrip = wv_para(7)


  base_tmp = 10.d0
  log_tmp = 6.1071d0

  ! good down to -100 C

  es = pow_c(base_tmp, (-9.09718d0*(h2otrip/t-1.d0)-3.56654d0* &
       log10_c(h2otrip/t)+0.876793d0*(1.d0-t/h2otrip)+ &
       log10_c(log_tmp)))*100.d0
  !es = pow_c(10.d0, (-9.09718d0*(h2otrip/t-1.d0)-3.56654d0* &
  !     log10_c(h2otrip/t)+0.876793d0*(1.d0-t/h2otrip)+ &
  !     log10_c(6.1071d0)))*100.d0

end function GoffGratch_svp_ice

! Murphy & Koop (2005)

 function MurphyKoop_svp_water(t) result(es)
  implicit none
  real(8), intent(in) :: t  ! Temperature in Kelvin
  real(8) :: es             ! SVP in Pa

  ! (good for 123 < T < 332 K)
  es = exp(54.842763d0 - (6763.22d0 / t) - (4.210d0 * log_c(t)) + &
       (0.000367d0 * t) + (tanh(0.0415d0 * (t - 218.8d0)) * &
       (53.878d0 - (1331.22d0 / t) - (9.44523d0 * log_c(t)) + &
       0.014025d0 * t)))

end function MurphyKoop_svp_water

 function MurphyKoop_svp_ice(t) result(es)
  implicit none
  real(8), intent(in) :: t  ! Temperature in Kelvin
  real(8) :: es             ! SVP in Pa

  ! (good down to 110 K)
  es = exp(9.550426d0 - (5723.265d0 / t) + (3.53068d0 * log_c(t)) &
       - (0.00728332d0 * t))

end function MurphyKoop_svp_ice

! Old CAM implementation, also labelled Goff & Gratch (1946)

! The water formula differs only due to compiler-dependent order of
! operations, so differences are roundoff level, usually 0.

! The ice formula gives fairly close answers to the current
! implementation, but has been rearranged, and uses the
! 1 atm melting point of water as the triple point.
! Differences are thus small but above roundoff.

! A curious fact: although using the melting point of water was
! probably a mistake, it mildly improves accuracy for ice svp,
! since it compensates for a systematic error in Goff & Gratch.

 function OldGoffGratch_svp_water(t) result(es)
  implicit none
  real(8), intent(in) :: t
  real(8) :: es
  real(8) :: ps, e1, e2, f1, f2, f3, f4, f5, f

  real(8), parameter :: tboil = 373.16d0

  ps = 1013.246d0
  e1 = 11.344d0*(1.0d0 - t/tboil)
  e2 = -3.49149d0*(tboil/t - 1.0d0)
  f1 = -7.90298d0*(tboil/t - 1.0d0)
  f2 = 5.02808d0*log10_c(tboil/t)
  f3 = -1.3816d0*(10.0d0**e1 - 1.0d0)/10000000.0d0
  f4 = 8.1328d0*(10.0d0**e2 - 1.0d0)/1000.0d0
  f5 = log10_c(ps)
  f  = f1 + f2 + f3 + f4 + f5

  es = (10.0d0**f)*100.0d0
  
end function OldGoffGratch_svp_water

 function OldGoffGratch_svp_ice(t, wv_para) result(es)
  implicit none
  real(8), intent(in) :: t
  real(8), intent(in) :: wv_para(PARASIZE)

  real(8) :: es
  real(8) :: term1, term2, term3
  real(8) :: tmelt 
  tmelt = wv_para(6)

  term1 = 2.01889049d0/(tmelt/t)
  term2 = 3.56654d0*log_c(tmelt/t)
  term3 = 20.947031d0*(tmelt/t)

  es = 575.185606d10*exp(-(term1 + term2 + term3))
  
end function OldGoffGratch_svp_ice

! Bolton (1980)
! zm_conv deep convection scheme contained this SVP calculation.
! It appears to be from D. Bolton, 1980, Monthly Weather Review.
! Unlike the other schemes, no distinct ice formula is associated
! with it. (However, a Bolton ice formula exists in CLUBB.)

! The original formula used degrees C, but this function
! takes Kelvin and internally converts.

 function Bolton_svp_water(t, wv_para) result(es)
  implicit none
  real(8),parameter :: c1 = 611.2d0
  real(8),parameter :: c2 = 17.67d0
  real(8),parameter :: c3 = 243.5d0

  real(8), intent(in) :: t  ! Temperature in Kelvin
  real(8), intent(in) :: wv_para(PARASIZE)
  real(8) :: es             ! SVP in Pa
  real(8) :: tmelt
  tmelt = wv_para(6)

  es = c1*exp( (c2*(t - tmelt))/((t - tmelt)+c3) )

end function Bolton_svp_water

end module wv_sat_methods_cpe_xk
