!/*************************************************************************
!	> File Name: wv_sat_methods_cpe.F90
!	> Author: Xu Kai
!	> Created Time: 2018年11月29日 星期四 09时35分30秒
! ************************************************************************/

        module wv_sat_methods_cpe
        
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
        
        !integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
        
        real(8) :: tmelt   ! Melting point of water at 1 atm (K)
        real(8) :: h2otrip ! Triple point temperature of water (K)
        real(8) :: tboil   ! Boiling point of water at 1 atm (K)
        
        real(8) :: ttrice  ! Ice-water transition range
        
        real(8) :: epsilo  ! Ice-water transition range
        real(8) :: omeps   ! 1._r8 - epsilo
        
        ! Indices representing individual schemes
        integer, parameter :: Invalid_idx = -1
        integer, parameter :: OldGoffGratch_idx = 0
        integer, parameter :: GoffGratch_idx = 1
        integer, parameter :: MurphyKoop_idx = 2
        integer, parameter :: Bolton_idx = 3
        
        ! Index representing the current default scheme.
        integer, parameter :: initial_default_idx = GoffGratch_idx
        integer :: default_idx = initial_default_idx
        
        public wv_sat_methods_init
        public wv_sat_get_scheme_idx
        public wv_sat_valid_idx
        
        public wv_sat_set_default
        public wv_sat_reset_default
        
        public wv_sat_svp_water
        public wv_sat_svp_ice
        public wv_sat_svp_trans
        
        ! pressure -> humidity conversion
        public wv_sat_svp_to_qsat
        
        ! Combined qsat operations
        public wv_sat_qsat_water
        public wv_sat_qsat_ice
        public wv_sat_qsat_trans

        !public shr_spfn_gamma_nonintrinsic_r8
        
        contains

        !---------------------------------------------------------------------
        ! ADMINISTRATIVE FUNCTIONS
        !---------------------------------------------------------------------
        
        ! Get physical constants
        subroutine wv_sat_methods_init(tmelt_in, h2otrip_in, tboil_in, &
             ttrice_in, epsilo_in)
          !integer, intent(in) :: kind
          real(8), intent(in) :: tmelt_in
          real(8), intent(in) :: h2otrip_in
          real(8), intent(in) :: tboil_in
          real(8), intent(in) :: ttrice_in
          real(8), intent(in) :: epsilo_in
          !character(len=*), intent(out)  :: errstring
        
          !errstring = ' '
        
          !if (kind /= 8) then
          !   !write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
          !   !     kind,' was input kind but ',r8,' is internal kind.'
          !   return
          !end if
        
          !if (ttrice_in < 0.d0) then
          !   !write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
          !   !     ttrice_in,' was input for ttrice, but negative range is invalid.'
          !   return
          !end if
        
          tmelt = tmelt_in
          h2otrip = h2otrip_in
          tboil = tboil_in
          ttrice = ttrice_in
          epsilo = epsilo_in
        
          omeps = 1.d0 - epsilo
        
        end subroutine wv_sat_methods_init
        
        ! Look up index by name.
        pure function wv_sat_get_scheme_idx(name) result(idx)
          character(len=*), intent(in) :: name
          integer :: idx
          
          select case (name)
          case("GoffGratch")
             idx = GoffGratch_idx
          case("MurphyKoop")
             idx = MurphyKoop_idx
          case("OldGoffGratch")
             idx = OldGoffGratch_idx
          case("Bolton")
             idx = Bolton_idx
          case default
             idx = Invalid_idx
          end select
        
        end function wv_sat_get_scheme_idx
        
        ! Check validity of an index from the above routine.
        pure function wv_sat_valid_idx(idx) result(status)
          integer, intent(in) :: idx
          logical :: status
        
          status = (idx /= Invalid_idx)
        
        end function wv_sat_valid_idx
        
        ! Set default scheme (otherwise, Goff & Gratch is default)
        ! Returns a logical representing success (.true.) or
        ! failure (.false.).
        function wv_sat_set_default(name) result(status)
          character(len=*), intent(in) :: name
          logical :: status
        
          ! Don't want to overwrite valid default with invalid,
          ! so assign to temporary and check it first.
          integer :: tmp_idx
        
          tmp_idx = wv_sat_get_scheme_idx(name)
        
          status = wv_sat_valid_idx(tmp_idx)
        
          if (status) default_idx = tmp_idx
        
        end function wv_sat_set_default
        
        ! Reset default scheme to initial value.
        ! The same thing can be accomplished with wv_sat_set_default;
        ! the real reason to provide this routine is to reset the
        ! module for testing purposes.
        subroutine wv_sat_reset_default()
        
          default_idx = initial_default_idx
        
        end subroutine wv_sat_reset_default
        
        !---------------------------------------------------------------------
        ! UTILITIES
        !---------------------------------------------------------------------
        
        ! Get saturation specific humidity given pressure and SVP.
        ! Specific humidity is limited to range 0-1.
        !elemental function wv_sat_svp_to_qsat(es, p) result(qs)
        function wv_sat_svp_to_qsat(es, p) result(qs)
        
          real(8), intent(in) :: es  ! SVP
          real(8), intent(in) :: p   ! Current pressure.
          real(8) :: qs
        
          ! If pressure is less than SVP, set qs to maximum of 1.
          if ( (p - es) <= 0.d0 ) then
             qs = 1.0d0
          else
             qs = epsilo*es / (p - omeps*es)
          end if
        
        end function wv_sat_svp_to_qsat
        
        !elemental subroutine wv_sat_qsat_water(t, p, es, qs, idx)
        subroutine wv_sat_qsat_water(t, p, es, qs, idx)
          !------------------------------------------------------------------!
          ! Purpose:                                                         !
          !   Calculate SVP over water at a given temperature, and then      !
          !   calculate and return saturation specific humidity.             !
          !------------------------------------------------------------------!
        
          ! Inputs
          real(8), intent(in) :: t    ! Temperature
          real(8), intent(in) :: p    ! Pressure
          ! Outputs
          real(8), intent(out) :: es  ! Saturation vapor pressure
          real(8), intent(out) :: qs  ! Saturation specific humidity
        
          integer,  intent(in), optional :: idx ! Scheme index
        
          es = wv_sat_svp_water(t, idx)
        
          qs = wv_sat_svp_to_qsat(es, p)
        
          ! Ensures returned es is consistent with limiters on qs.
          es = min(es, p)
        
        end subroutine wv_sat_qsat_water
        
        !elemental subroutine wv_sat_qsat_ice(t, p, es, qs, idx)
        subroutine wv_sat_qsat_ice(t, p, es, qs, idx)
          !------------------------------------------------------------------!
          ! Purpose:                                                         !
          !   Calculate SVP over ice at a given temperature, and then        !
          !   calculate and return saturation specific humidity.             !
          !------------------------------------------------------------------!
        
          ! Inputs
          real(8), intent(in) :: t    ! Temperature
          real(8), intent(in) :: p    ! Pressure
          ! Outputs
          real(8), intent(out) :: es  ! Saturation vapor pressure
          real(8), intent(out) :: qs  ! Saturation specific humidity
        
          integer,  intent(in), optional :: idx ! Scheme index
        
          es = wv_sat_svp_ice(t, idx)
        
          qs = wv_sat_svp_to_qsat(es, p)
        
          ! Ensures returned es is consistent with limiters on qs.
          es = min(es, p)
        
        end subroutine wv_sat_qsat_ice
        
        !elemental subroutine wv_sat_qsat_trans(t, p, es, qs, idx)
        subroutine wv_sat_qsat_trans(t, p, es, qs, idx)
          !------------------------------------------------------------------!
          ! Purpose:                                                         !
          !   Calculate SVP over ice at a given temperature, and then        !
          !   calculate and return saturation specific humidity.             !
          !------------------------------------------------------------------!
        
          ! Inputs
          real(8), intent(in) :: t    ! Temperature
          real(8), intent(in) :: p    ! Pressure
          ! Outputs
          real(8), intent(out) :: es  ! Saturation vapor pressure
          real(8), intent(out) :: qs  ! Saturation specific humidity
        
          integer,  intent(in), optional :: idx ! Scheme index
        
          es = wv_sat_svp_trans(t, idx)
        
          qs = wv_sat_svp_to_qsat(es, p)
        
          ! Ensures returned es is consistent with limiters on qs.
          es = min(es, p)
        
        end subroutine wv_sat_qsat_trans
        
        !---------------------------------------------------------------------
        ! SVP INTERFACE FUNCTIONS
        !---------------------------------------------------------------------
        
        !elemental function wv_sat_svp_water(t, idx) result(es)
        function wv_sat_svp_water(t, idx) result(es)
          real(8), intent(in) :: t
          integer,  intent(in), optional :: idx
          real(8) :: es
        
          integer :: use_idx
        
          if (present(idx)) then
             use_idx = idx
          else
             use_idx = default_idx
          end if
        
          select case (use_idx)
          case(GoffGratch_idx)
             es = GoffGratch_svp_water(t)
          case(MurphyKoop_idx)
             es = MurphyKoop_svp_water(t)
          case(OldGoffGratch_idx)
             es = OldGoffGratch_svp_water(t)
          case(Bolton_idx)
             es = Bolton_svp_water(t)
          end select
        
        end function wv_sat_svp_water
        
        !elemental function wv_sat_svp_ice(t, idx) result(es)
        function wv_sat_svp_ice(t, idx) result(es)
          real(8), intent(in) :: t
          integer,  intent(in), optional :: idx
          real(8) :: es
        
          integer :: use_idx
        
          if (present(idx)) then
             use_idx = idx
          else
             use_idx = default_idx
          end if
        
          select case (use_idx)
          case(GoffGratch_idx)
             es = GoffGratch_svp_ice(t)
          case(MurphyKoop_idx)
             es = MurphyKoop_svp_ice(t)
          case(OldGoffGratch_idx)
             es = OldGoffGratch_svp_ice(t)
          case(Bolton_idx)
             es = Bolton_svp_water(t)
          end select
        
        end function wv_sat_svp_ice
        
        !elemental function wv_sat_svp_trans(t, idx) result (es)
        function wv_sat_svp_trans(t, idx) result (es)
        
          real(8), intent(in) :: t
          integer,  intent(in), optional :: idx
        
          real(8) :: es
        
          real(8) :: esice      ! Saturation vapor pressure over ice
          real(8) :: weight     ! Intermediate scratch variable for es transition
        
        !
        ! Water
        !
          if (t >= (tmelt - ttrice)) then
             es = wv_sat_svp_water(t,idx)
          else
             es = 0.0d0
          end if
        
        !
        ! Ice
        !
          if (t < tmelt) then
        
             esice = wv_sat_svp_ice(t,idx)
        
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
        
        !elemental function GoffGratch_svp_water(t) result(es)
        function GoffGratch_svp_water(t) result(es)
          real(8), intent(in) :: t  ! Temperature in Kelvin
          real(8) :: es             ! SVP in Pa

          integer(8) :: lg10, power
          real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5
        
        
          call math_agent_log10_c(lg10)
          call math_agent_pow_c(power)

          tmp1 = 10.d0
          tmp2 = tboil/t

          call math_agent_1i1o(lg10, tmp2, tmp3)
          tmp4 = 11.344d0*(1.d0-t/tboil)
          call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          tmp4 = (-3.49149d0*(tboil/t-1.d0))
          call math_agent_2i1o(power, tmp1, tmp4, tmp5)
          tmp4 = 1013.246d0
          call math_agent_1i1o(lg10, tmp4, tmp1)

          tmp4 = -7.90298d0*(tboil/t-1.d0)+ &
               5.02808d0*tmp3- &
               1.3816d-7*(tmp2-1.d0)+ &
               8.1328d-3*(tmp5-1.d0)+ &
               tmp1

          tmp1 = 10.d0
          call math_agent_2i1o(power, tmp1, tmp4, tmp2)

          es = tmp2 * 100.d0

          !es = 10.d0

          ! uncertain below -70 C
          !es = 10.d0**(-7.90298d0*(tboil/t-1.d0)+ &
          !     5.02808d0*log10(tboil/t)- &
          !     1.3816d-7*(10.d0**(11.344d0*(1.d0-t/tboil))-1.d0)+ &
          !     8.1328d-3*(10.d0**(-3.49149d0*(tboil/t-1.d0))-1.d0)+ &
          !     log10(1013.246d0))*100.d0

        
        end function GoffGratch_svp_water
        
        !elemental function GoffGratch_svp_ice(t) result(es)
        function GoffGratch_svp_ice(t) result(es)
          real(8), intent(in) :: t  ! Temperature in Kelvin
          real(8) :: es             ! SVP in Pa

          integer(8) :: lg10, power
          real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5
        
          call math_agent_log10_c(lg10)
          call math_agent_pow_c(power)
        
          tmp1 = h2otrip/t
          call math_agent_1i1o(lg10, tmp1, tmp2)
          tmp1 = 6.1071d0
          call math_agent_1i1o(lg10, tmp1, tmp3)
          tmp1 = -9.09718d0*(h2otrip/t-1.d0)-3.56654d0* &
               tmp2+0.876793d0*(1.d0-t/h2otrip)+ &
               tmp3
          tmp2 = 10.d0
          call math_agent_2i1o(power, tmp2, tmp1, tmp3)
          es = tmp3*100.d0

          !! good down to -100 C
          !es = 10.d0**(-9.09718d0*(h2otrip/t-1.d0)-3.56654d0* &
          !     log10(h2otrip/t)+0.876793d0*(1.d0-t/h2otrip)+ &
          !     log10(6.1071d0))*100.d0
        
        end function GoffGratch_svp_ice
        
        ! Murphy & Koop (2005)
        
        !elemental function MurphyKoop_svp_water(t) result(es)
        function MurphyKoop_svp_water(t) result(es)
          real(8), intent(in) :: t  ! Temperature in Kelvin
          real(8) :: es             ! SVP in Pa
        
          ! (good for 123 < T < 332 K)
          es = exp(54.842763d0 - (6763.22d0 / t) - (4.210d0 * log(t)) + &
               (0.000367d0 * t) + (tanh(0.0415d0 * (t - 218.8d0)) * &
               (53.878d0 - (1331.22d0 / t) - (9.44523d0 * log(t)) + &
               0.014025d0 * t)))
        
        end function MurphyKoop_svp_water
        
        !elemental function MurphyKoop_svp_ice(t) result(es)
        function MurphyKoop_svp_ice(t) result(es)
          real(8), intent(in) :: t  ! Temperature in Kelvin
          real(8) :: es             ! SVP in Pa
        
          ! (good down to 110 K)
          es = exp(9.550426d0 - (5723.265d0 / t) + (3.53068d0 * log(t)) &
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
        
        !elemental function OldGoffGratch_svp_water(t) result(es)
        function OldGoffGratch_svp_water(t) result(es)
          real(8), intent(in) :: t
          real(8) :: es
          real(8) :: ps, e1, e2, f1, f2, f3, f4, f5, f
        
          ps = 1013.246d0
          e1 = 11.344d0*(1.0d0 - t/tboil)
          e2 = -3.49149d0*(tboil/t - 1.0d0)
          f1 = -7.90298d0*(tboil/t - 1.0d0)
          f2 = 5.02808d0*log10(tboil/t)
          f3 = -1.3816d0*(10.0d0**e1 - 1.0d0)/10000000.0d0
          f4 = 8.1328d0*(10.0d0**e2 - 1.0d0)/1000.0d0
          f5 = log10(ps)
          f  = f1 + f2 + f3 + f4 + f5
        
          es = (10.0d0**f)*100.0d0
          
        end function OldGoffGratch_svp_water
        
        !elemental function OldGoffGratch_svp_ice(t) result(es)
        function OldGoffGratch_svp_ice(t) result(es)
          real(8), intent(in) :: t
          real(8) :: es
          real(8) :: term1, term2, term3
        
          term1 = 2.01889049d0/(tmelt/t)
          term2 = 3.56654d0*log(tmelt/t)
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
        
        !elemental function Bolton_svp_water(t) result(es)
        function Bolton_svp_water(t) result(es)
          real(8),parameter :: c1 = 611.2d0
          real(8),parameter :: c2 = 17.67d0
          real(8),parameter :: c3 = 243.5d0
        
          real(8), intent(in) :: t  ! Temperature in Kelvin
          real(8) :: es             ! SVP in Pa
        
          es = c1*exp( (c2*(t - tmelt))/((t - tmelt)+c3) )
        
        end function Bolton_svp_water
        
        end module wv_sat_methods_cpe
