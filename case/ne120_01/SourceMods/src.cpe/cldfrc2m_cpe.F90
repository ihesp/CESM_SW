module cldfrc2m_cpe
#define PARASIZE 12
! cloud fraction calculations

!use shr_kind_mod,     only: r8=>shr_kind_r8
!use spmd_utils,       only: masterproc
!use ppgrid,           only: pcols
!use physconst,        only: rair
use wv_saturation_cpe_xk,    only: qsat_water, svp_water, svp_ice
!use cam_logfile,      only: iulog
!use cam_abortutils,   only: endrun

implicit none
private
save

public ::           &
   astG_PDF_single, &
   astG_PDF,        &
   astG_RHU_single, &
   astG_RHU,        &
   aist_single,     &
   aist_vector
interface 
    real(8) function exp_c(in1)
    real(8), intent(in) :: in1
    end function

    real(8) function acos_c(in1)
    real(8), intent(in) :: in1
    end function

    real(8) function cos_c(in1)
    real(8), intent(in) :: in1
    end function
    real(8) function pow_c(in1, in2)
    real(8), intent(in) :: in1, in2
    end function

    !real(8) function sqrt_c(in1)
    !real(8), intent(in) :: in1
    !end function


end interface




!   rhmini_const,    &
!   rhmaxi_const


!! Namelist variables
!real(8) :: cldfrc2m_rhmini            ! Minimum rh for ice cloud fraction > 0.
!real(8) :: cldfrc2m_rhmaxi
!
!! -------------------------- !
!! Parameters for Ice Stratus !
!! -------------------------- !
!real(8) :: rhmini_const                 ! Minimum rh for ice cloud fraction > 0.
!real(8) :: rhmaxi_const
!
!real(8),  parameter :: qist_min     = 1.d-7      ! Minimum in-stratus ice IWC constraint [ kg/kg ]
!real(8),  parameter :: qist_max     = 5.d-3      ! Maximum in-stratus ice IWC constraint [ kg/kg ]
!
!! ----------------------------- !
!! Parameters for Liquid Stratus !
!! ----------------------------- !
!
!logical,  parameter  :: CAMstfrac    = .false.    ! If .true. (.false.),
!                                                  ! use Slingo (triangular PDF-based) liquid stratus fraction
!logical,  parameter  :: freeze_dry   = .false.    ! If .true., use 'freeze dry' in liquid stratus fraction formula
!real(8)             :: rhminl_const              ! Critical RH for low-level  liquid stratus clouds
!real(8)             :: rhminl_adj_land_const     ! rhminl adjustment for snowfree land
!real(8)             :: rhminh_const              ! Critical RH for high-level liquid stratus clouds
!real(8)             :: premit                    ! Top    height for mid-level liquid stratus fraction
!real(8)             :: premib                    ! Bottom height for mid-level liquid stratus fraction
!integer              :: iceopt                    ! option for ice cloud closure 
!                                                  ! 1=wang & sassen 2=schiller (iciwc)  
!                                                  ! 3=wood & field, 4=Wilson (based on smith)
!                                                  ! 5=modified slingo (ssat & empyt cloud)        
!real(8)             :: icecrit                   ! Critical RH for ice clouds in Wilson & Ballard closure
                                                  ! ( smaller = more ice clouds )

!================================================================================================
contains

!================================================================================================


subroutine astG_PDF_single(U, p, qv, landfrac, snowh, a, Ga, &
                            premib, premit, & 
                            orhmin, &
                           rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! analytical formulation of triangular PDF.                 !
   ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', !
   ! so using constant 'dV' assume that width is proportional  !
   ! to the saturation specific humidity.                      !
   !    dV ~ 0.1.                                              !
   !    cldrh : RH of in-stratus( = 1 if no supersaturation)   !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.  In fact, it does not    !
   ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that   !
   ! they will produce the same results.                       !
   ! --------------------------------------------------------- !

   implicit none
   real(8), intent(in)  :: U                     ! Relative humidity
   real(8), intent(in)  :: p                     ! Pressure [Pa]
   real(8), intent(in)  :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(8), intent(in)  :: landfrac              ! Land fraction
   real(8), intent(in)  :: snowh                 ! Snow depth (liquid water equivalent)

   real(8), intent(out) :: a                     ! Stratus fraction
   real(8), intent(out) :: Ga                    ! dU/da
   real(8), optional, intent(out) :: orhmin      ! Critical RH

   real(8), optional, intent(in)  :: rhminl_in          ! Critical relative humidity for low-level  liquid stratus
   real(8), optional, intent(in)  :: rhminl_adj_land_in ! Adjustment drop of rhminl over the land
   real(8), optional, intent(in)  :: rhminh_in          ! Critical relative humidity for high-level liquid stratus

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(8) dV                                    ! Width of triangular PDF
   real(8) cldrh                                 ! RH of stratus cloud
   real(8) rhmin                                 ! Critical RH
   real(8) rhwght
                            
   real(8) :: rhminl
   real(8) :: rhminl_adj_land
   real(8) :: rhminh

   ! Statement functions
   logical land
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real(8), intent(in)  :: premib    
 real(8), intent(in)  :: premit 
!logical,  parameter  :: freeze_dry   = .false.    ! If .true., use 'freeze dry' in liquid stratus fraction formula
logical  :: freeze_dry    ! If .true., use 'freeze dry' in liquid stratus fraction formula
real(8) ::  tmp1
!, tmp2, tmp3
freeze_dry   = .false.
!xukai-------------------------------------------------------------------------------------
   land = nint(landfrac) == 1

   ! ---------- !
   ! Parameters !
   ! ---------- !

   cldrh  = 1.0d0
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !rhminl          = rhminl_const
   !rhminl_adj_land = rhminl_adj_land_const
   !rhminh          = rhminh_const
!xukai-------------------------------------------------------------------------------------
 
   !rhminl = rhminl_const
   if (present(rhminl_in)) rhminl = rhminl_in
   !rhminl_adj_land = rhminl_adj_land_const
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
   !rhminh = rhminh_const
   if (present(rhminh_in)) rhminh = rhminh_in

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   if( p .ge. premib ) then

       if( land .and. (snowh.le.0.000001d0) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif

       dV = cldrh - rhmin

       if( U .ge. 1.d0 ) then
           a  = 1.d0
           Ga = 1.d10
       elseif( U .gt. (cldrh-dV/6.d0) .and. U .lt. 1.d0 ) then
            tmp1 = 2.d0 / 3.d0
           a  = 1.d0 - pow_c((-3.d0/sqrt(2.d0)*(U-cldrh)/dV), tmp1)
           Ga = dV/sqrt(2.d0)*sqrt(1.d0-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6.d0) ) then
            tmp1 = 2.d0
           a  = 4.d0*pow_c((cos_c((1.d0/3.d0)*(acos_c((3.d0/2.d0/sqrt(2.d0))* & 
                      (1.d0+(U-cldrh)/dV))-2.d0*3.141592d0))), tmp1)
           Ga = dV/sqrt(2.d0)*(1.d0/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0.d0
           Ga = 1.d10
       endif

       if( freeze_dry ) then
           a  = a *max(0.15d0,min(1.0d0,qv/0.0030d0))
           Ga = Ga/max(0.15d0,min(1.0d0,qv/0.0030d0)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       dV    = cldrh - rhmin

       if( U .ge. 1.d0 ) then
           a  = 1.d0
           Ga = 1.d10
       elseif( U .gt. (cldrh-dV/6.d0) .and. U .lt. 1.d0 ) then
           tmp1 = 2.d0/3.d0 
           a  = 1.d0 - pow_c((-3.d0/sqrt(2.d0)*(U-cldrh)/dV), tmp1)
           Ga = dV/sqrt(2.d0)*sqrt(1.d0-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6.d0) ) then
           tmp1 = 2.d0
           a  = 4.d0*pow_c((cos_c((1.d0/3.d0)*(acos_c((3.d0/2.d0/sqrt(2.d0))* & 
                      (1.d0+(U-cldrh)/dV))-2.d0*3.141592d0))), tmp1)
           Ga = dV/sqrt(2.d0)*(1.d0/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0.d0
           Ga = 1.d10
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land .and. (snowh.le.0.000001d0) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0d0-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0d0-rhwght)
     ! endif

       dV    = cldrh - rhmin

       if( U .ge. 1.d0 ) then
           a  = 1.d0
           Ga = 1.d10
       elseif( U .gt. (cldrh-dV/6.d0) .and. U .lt. 1.d0 ) then
           tmp1 = 2.d0/3.d0
           a  = 1.d0 - pow_c((-3.d0/sqrt(2.d0)*(U-cldrh)/dV), tmp1)
           Ga = dV/sqrt(2.d0)*sqrt(1.d0-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6.d0) ) then
           tmp1 = 2.d0
           a  = 4.d0*pow_c((cos_c((1.d0/3.d0)*(acos_c((3.d0/2.d0/sqrt(2.d0))* & 
                         (1.d0+(U-cldrh)/dV))-2.d0*3.141592d0))), tmp1)
           Ga = dV/sqrt(2.d0)*(1.d0/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0.d0
           Ga = 1.d10
       endif

   endif

   if (present(orhmin)) orhmin = rhmin

end subroutine astG_PDF_single

!================================================================================================

subroutine astG_PDF(U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol, &
                    premib, premit, &
                    rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! analytical formulation of triangular PDF.                 !
   ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', !
   ! so using constant 'dV' assume that width is proportional  !
   ! to the saturation specific humidity.                      !
   !    dV ~ 0.1.                                              !
   !    cldrh : RH of in-stratus( = 1 if no supersaturation)   !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.  In fact, it does not    !
   ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that   !
   ! they will produce the same results.                       !
   ! --------------------------------------------------------- !

   implicit none
   real(8), intent(in)  :: U_in           ! Relative humidity
   real(8), intent(in)  :: p_in           ! Pressure [Pa]
   real(8), intent(in)  :: qv_in          ! Grid-mean water vapor specific humidity [kg/kg]
   real(8), intent(in)  :: landfrac_in    ! Land fraction
   real(8), intent(in)  :: snowh_in       ! Snow depth (liquid water equivalent)

   real(8), intent(out) :: a_out          ! Stratus fraction
   real(8), intent(out) :: Ga_out         ! dU/da
   integer,  intent(in)  :: ncol

   real(8), optional, intent(in)  :: rhminl_in                ! Critical relative humidity for low-level  liquid stratus
   real(8), optional, intent(in)  :: rhminl_adj_land_in       ! Adjustment drop of rhminl over the land
   real(8), optional, intent(in)  :: rhminh_in                ! Critical relative humidity for high-level liquid stratus

   real(8)              :: rhminl                ! Critical relative humidity for low-level  liquid stratus
   real(8)              :: rhminl_adj_land       ! Adjustment drop of rhminl over the land
   real(8)              :: rhminh                ! Critical relative humidity for high-level liquid stratus

   real(8)              :: U                     ! Relative humidity
   real(8)              :: p                     ! Pressure [Pa]
   real(8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(8)              :: landfrac              ! Land fraction
   real(8)              :: snowh                 ! Snow depth (liquid water equivalent)

   real(8)              :: a                     ! Stratus fraction
   real(8)              :: Ga                    ! dU/da

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(8) dV                                    ! Width of triangular PDF
   real(8) cldrh                                 ! RH of stratus cloud
   real(8) rhmin                                 ! Critical RH
   real(8) rhwght
                            
   ! Statement functions
   logical land
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real(8), intent(in)  :: premib    
 real(8), intent(in)  :: premit 
!logical,  parameter  :: freeze_dry   = .false.    ! If .true., use 'freeze dry' in liquid stratus fraction formula
logical  :: freeze_dry  ! If .true., use 'freeze dry' in liquid stratus fraction formula
real(8) :: tmp1
!, tmp2, tmp3
freeze_dry   = .false. 
!xukai-------------------------------------------------------------------------------------
 
   land = nint(landfrac_in) == 1

   ! ---------- !
   ! Parameters !
   ! ---------- !

   cldrh  = 1.0d0
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !rhminl          = rhminl_const
   !rhminl_adj_land = rhminl_adj_land_const
   !rhminh          = rhminh_const
!xukai-------------------------------------------------------------------------------------
   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   a_out  = 0.d0
   Ga_out = 0.d0

   !do i = 1, ncol

   U        = U_in      
   p        = p_in        
   qv       = qv_in       
   landfrac = landfrac_in 
   snowh    = snowh_in    

   if (present(rhminl_in))          rhminl          = rhminl_in      
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
   if (present(rhminh_in))          rhminh          = rhminh_in

   if( p .ge. premib ) then

       if( land .and. (snowh.le.0.000001d0) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif

       dV = cldrh - rhmin

       if( U .ge. 1.d0 ) then
           a  = 1.d0
           Ga = 1.d10
       elseif( U .gt. (cldrh-dV/6.d0) .and. U .lt. 1.d0 ) then
            tmp1 = 2.d0/3.d0
           a  = 1.d0 - pow_c((-3.d0/sqrt(2.d0)*(U-cldrh)/dV), tmp1)
           Ga = dV/sqrt(2.d0)*sqrt(1.d0-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6.d0) ) then
            tmp1 = 2.d0
           a  = 4.d0*pow_c((cos_c((1.d0/3.d0)*(acos_c((3.d0/2.d0/sqrt(2.d0))* & 
                      (1.d0+(U-cldrh)/dV))-2.d0*3.141592d0))), tmp1)
           Ga = dV/sqrt(2.d0)*(1.d0/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0.d0
           Ga = 1.d10
       endif

       if( freeze_dry ) then
           a  = a *max(0.15d0,min(1.0d0,qv/0.0030d0))
           Ga = Ga/max(0.15d0,min(1.0d0,qv/0.0030d0)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       dV    = cldrh - rhmin

       if( U .ge. 1.d0 ) then
           a  = 1.d0
           Ga = 1.d10
       elseif( U .gt. (cldrh-dV/6.d0) .and. U .lt. 1.d0 ) then
           tmp1 = 2.d0/3.d0
           a  = 1.d0 - pow_c((-3.d0/sqrt(2.d0)*(U-cldrh)/dV), tmp1)
           Ga = dV/sqrt(2.d0)*sqrt(1.d0-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6.d0) ) then
           tmp1 = 2.d0
           a  = 4.d0*pow_c((cos_c((1.d0/3.d0)*(acos_c((3.d0/2.d0/sqrt(2.d0))* & 
                      (1.d0+(U-cldrh)/dV))-2.d0*3.141592d0))), tmp1)
           Ga = dV/sqrt(2.d0)*(1.d0/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0.d0
           Ga = 1.d10
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land .and. (snowh.le.0.000001d0) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0d0-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0d0-rhwght)
     ! endif

       dV    = cldrh - rhmin

       if( U .ge. 1.d0 ) then
           a  = 1.d0
           Ga = 1.d10
       elseif( U .gt. (cldrh-dV/6.d0) .and. U .lt. 1.d0 ) then
            tmp1 = 2.d0/3.d0
           a  = 1.d0 - pow_c((-3.d0/sqrt(2.d0)*(U-cldrh)/dV), tmp1)
           Ga = dV/sqrt(2.d0)*sqrt(1.d0-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6.d0) ) then
            tmp1 = 2.d0
           a  = 4.d0*pow_c((cos_c((1.d0/3.d0)*(acos_c((3.d0/2.d0/sqrt(2.d0))* & 
                         (1.d0+(U-cldrh)/dV))-2.d0*3.141592d0))), tmp1)
           Ga = dV/sqrt(2.d0)*(1.d0/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0.d0
           Ga = 1.d10
       endif

   endif

   a_out  = a
   Ga_out = Ga 

   !enddo

end subroutine astG_PDF
!================================================================================================

subroutine astG_RHU_single(U, p, qv, landfrac, snowh, a, Ga, &
                            premib, premit, & 
                                orhmin, &
                              rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! CAM35 cloud fraction formula.                             !
   ! Below is valid only for CAMUW at 1.9x2.5 fv dynamics core !  
   ! For the other cases, I should re-define 'rhminl,rhminh' & !
   ! 'premib,premit'.                                          !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.                          !
   ! --------------------------------------------------------- !
   implicit none

   real(8), intent(in)  :: U               ! Relative humidity
   real(8), intent(in)  :: p               ! Pressure [Pa]
   real(8), intent(in)  :: qv              ! Grid-mean water vapor specific humidity [kg/kg]
   real(8), intent(in)  :: landfrac        ! Land fraction
   real(8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real(8), intent(out) :: a               ! Stratus fraction
   real(8), intent(out) :: Ga              ! dU/da
   real(8), optional, intent(out) :: orhmin ! Critical RH

   real(8), optional, intent(in)  :: rhminl_in          ! Critical relative humidity for low-level  liquid stratus
   real(8), optional, intent(in)  :: rhminl_adj_land_in ! Adjustment drop of rhminl over the land
   real(8), optional, intent(in)  :: rhminh_in          ! Critical relative humidity for high-level liquid stratus

   ! Local variables
   real(8) rhmin                                 ! Critical RH
   real(8) rhdif                                 ! Factor for stratus fraction
   real(8) rhwght

   real(8) :: rhminl
   real(8) :: rhminl_adj_land
   real(8) :: rhminh

   ! Statement functions
   logical land
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real(8), intent(in) :: premib 
    real(8), intent(in) :: premit 
!logical,  parameter  :: freeze_dry   = .false.    ! If .true., use 'freeze dry' in liquid stratus fraction formula
logical :: freeze_dry    ! If .true., use 'freeze dry' in liquid stratus fraction formula
real(8) :: tmp1
!, tmp2, tmp3
freeze_dry   = .false. 
!xukai-------------------------------------------------------------------------------------
 
   land = nint(landfrac) == 1

!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !rhminl          = rhminl_const
   !rhminl_adj_land = rhminl_adj_land_const
   !rhminh          = rhminh_const
!xukai-------------------------------------------------------------------------------------
 
   !rhminl = rhminl_const
   if (present(rhminl_in)) rhminl = rhminl_in
   !rhminl_adj_land = rhminl_adj_land_const
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
   !rhminh = rhminh_const
   if (present(rhminh_in)) rhminh = rhminh_in

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   if( p .ge. premib ) then

       if( land .and. (snowh.le.0.000001d0) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif
       rhdif = (U-rhmin)/(1.0d0-rhmin)
       !a  = min(1.d0,pow_c((max(rhdif,0.0d0)), 2.d0)) 
       tmp1 = 2.d0
       a  = min(1.d0,pow_c((max(rhdif,0.0d0)), tmp1)) 
       if( (U.ge.1.d0) .or. (U.le.rhmin) ) then
            Ga = 1.d20
       else          
            Ga = 0.5d0*(1.d0-rhmin)*((1.d0-rhmin)/(U-rhmin))
       endif
       if( freeze_dry ) then
           a  = a*max(0.15d0,min(1.0d0,qv/0.0030d0))
           Ga = Ga/max(0.15d0,min(1.0d0,qv/0.0030d0)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       rhdif = (U-rhmin)/(1.0d0-rhmin)
       !a  = min(1.d0,pow_c((max(rhdif,0.d0)), 2.d0))
       tmp1 = 2.d0
       a  = min(1.d0,pow_c((max(rhdif,0.d0)), tmp1))
       if( (U.ge.1.d0) .or. (U.le.rhmin) ) then
            Ga = 1.d20
       else          
            Ga = 0.5d0*(1.d0-rhmin)*((1.d0-rhmin)/(U-rhmin))
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land .and. (snowh.le.0.000001d0) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0d0-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0d0-rhwght)
     ! endif

       rhdif = (U-rhmin)/(1.0d0-rhmin)
       !a  = min(1.d0,pow_c((max(rhdif,0.d0)), 2.d0))
       tmp1 = 2.d0
       a  = min(1.d0,pow_c((max(rhdif,0.d0)), tmp1))
       if( (U.ge.1.d0) .or. (U.le.rhmin) ) then
            Ga = 1.d10
       else          
            Ga = 0.5d0*(1.d0-rhmin)*((1.d0-rhmin)/(U-rhmin))
       endif

   endif

   if (present(orhmin)) orhmin = rhmin

end subroutine astG_RHU_single

!================================================================================================

subroutine astG_RHU(U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol, &
                    premib, premit, &
                    rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! CAM35 cloud fraction formula.                             !
   ! Below is valid only for CAMUW at 1.9x2.5 fv dynamics core !  
   ! For the other cases, I should re-define 'rhminl,rhminh' & !
   ! 'premib,premit'.                                          !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.                          !
   ! --------------------------------------------------------- !

   implicit none
   real(8), intent(in)  :: U_in           ! Relative humidity
   real(8), intent(in)  :: p_in           ! Pressure [Pa]
   real(8), intent(in)  :: qv_in          ! Grid-mean water vapor specific humidity [kg/kg]
   real(8), intent(in)  :: landfrac_in    ! Land fraction
   real(8), intent(in)  :: snowh_in       ! Snow depth (liquid water equivalent)

   real(8), intent(out) :: a_out          ! Stratus fraction
   real(8), intent(out) :: Ga_out         ! dU/da
   integer,  intent(in)  :: ncol

   real(8), optional, intent(in)  :: rhminl_in          ! Critical relative humidity for low-level  liquid stratus
   real(8), optional, intent(in)  :: rhminl_adj_land_in ! Adjustment drop of rhminl over the land
   real(8), optional, intent(in)  :: rhminh_in          ! Critical relative humidity for high-level liquid stratus

   real(8)              :: U                     ! Relative humidity
   real(8)              :: p                     ! Pressure [Pa]
   real(8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(8)              :: landfrac              ! Land fraction
   real(8)              :: snowh                 ! Snow depth (liquid water equivalent)

   real(8)              :: rhminl                ! Critical relative humidity for low-level  liquid stratus
   real(8)              :: rhminl_adj_land       ! Adjustment drop of rhminl over the land
   real(8)              :: rhminh                ! Critical relative humidity for high-level liquid stratus

   real(8)              :: a                     ! Stratus fraction
   real(8)              :: Ga                    ! dU/da

   ! Local variables
   integer  i
   real(8) rhmin                                 ! Critical RH
   real(8) rhdif                                 ! Factor for stratus fraction
   real(8) rhwght

   ! Statement functions
   logical land
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real(8), intent(in) :: premib 
    real(8), intent(in) :: premit 
!logical,  parameter  :: freeze_dry   = .false.    ! If .true., use 'freeze dry' in liquid stratus fraction formula
logical :: freeze_dry 
    real(8) :: tmp1
    !, tmp2, tmp3
 freeze_dry = .false.    ! If .true., use 'freeze dry' in liquid stratus fraction formula
!xukai-------------------------------------------------------------------------------------
 
   land = nint(landfrac_in) == 1
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !rhminl          = rhminl_const
   !rhminl_adj_land = rhminl_adj_land_const
   !rhminh          = rhminh_const
!xukai-------------------------------------------------------------------------------------
 
   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   a_out = 0.d0
   Ga_out = 0.d0

   !do i = 1, ncol

   U        = U_in      
   p        = p_in        
   qv       = qv_in       
   landfrac = landfrac_in 
   snowh    = snowh_in    

   if (present(rhminl_in))          rhminl          = rhminl_in      
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
   if (present(rhminh_in))          rhminh          = rhminh_in

   if( p .ge. premib ) then

       if( land .and. (snowh.le.0.000001d0) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif
       rhdif = (U-rhmin)/(1.0d0-rhmin)
       !a  = min(1.d0,pow_c((max(rhdif,0.0d0)), 2.d0)) 
       tmp1 = 2.d0 
       a  = min(1.d0,pow_c((max(rhdif,0.0d0)), tmp1)) 
       if( (U.ge.1.d0) .or. (U.le.rhmin) ) then
            Ga = 1.d20
       else          
            Ga = 0.5d0*(1.d0-rhmin)*((1.d0-rhmin)/(U-rhmin))
       endif
       if( freeze_dry ) then
           a  = a*max(0.15d0,min(1.0d0,qv/0.0030d0))
           Ga = Ga/max(0.15d0,min(1.0d0,qv/0.0030d0)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       rhdif = (U-rhmin)/(1.0d0-rhmin)
       !a  = min(1.d0,pow_c((max(rhdif,0.d0)), 2.d0))
       tmp1 = 2.d0
       a  = min(1.d0,pow_c((max(rhdif,0.d0)), tmp1))
       if( (U.ge.1.d0) .or. (U.le.rhmin) ) then
            Ga = 1.d20
       else          
            Ga = 0.5d0*(1.d0-rhmin)*((1.d0-rhmin)/(U-rhmin))
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land .and. (snowh.le.0.000001d0) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0d0-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0d0-rhwght)
     ! endif

       rhdif = (U-rhmin)/(1.0d0-rhmin)
       !a  = min(1.d0,pow_c((max(rhdif,0.d0)), 2.d0))
       tmp1 = 2.d0
       a  = min(1.d0,pow_c((max(rhdif,0.d0)), tmp1))
       if( (U.ge.1.d0) .or. (U.le.rhmin) ) then
            Ga = 1.d10
       else          
            Ga = 0.5d0*(1.d0-rhmin)*((1.d0-rhmin)/(U-rhmin))
       endif

   endif

   a_out  = a
   Ga_out = Ga 

   !enddo

end subroutine astG_RHU

!================================================================================================

subroutine aist_single(qv, T, p, qi, landfrac, snowh, aist, &
                       wv_para, premib, premit, icecrit, iceopt, &
                       rhmaxi_in, rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in)

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   implicit none
   real(8), intent(in)  :: qv              ! Grid-mean water vapor[kg/kg]
   real(8), intent(in)  :: T               ! Temperature
   real(8), intent(in)  :: p               ! Pressure [Pa]
   real(8), intent(in)  :: qi              ! Grid-mean ice water content [kg/kg]
   real(8), intent(in)  :: landfrac        ! Land fraction
   real(8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real(8), intent(out) :: aist            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )
   real(8), intent(in)  :: wv_para(PARASIZE)

   real(8), optional, intent(in)  :: rhmaxi_in
   real(8), optional, intent(in)  :: rhmini_in          ! Critical relative humidity for               ice stratus
   real(8), optional, intent(in)  :: rhminl_in          ! Critical relative humidity for low-level  liquid stratus
   real(8), optional, intent(in)  :: rhminl_adj_land_in ! Adjustment drop of rhminl over the land
   real(8), optional, intent(in)  :: rhminh_in          ! Critical relative humidity for high-level liquid stratus

   ! Local variables
   real(8) rhmin                           ! Critical RH
   real(8) rhwght

   real(8) a,b,c,as,bs,cs                  ! Fit parameters
   real(8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(8) ttmp                            ! Limited temperature
   real(8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(8) rho                             ! Local air density
   real(8) esl                             ! Liq sat vapor pressure
   real(8) esi                             ! Ice sat vapor pressure
   real(8) ncf,phi                         ! Wilson and Ballard parameters
   real(8) es, qs

   real(8) rhi                             ! grid box averaged relative humidity over ice
   real(8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(8) mincld                          ! minimum ice cloud fraction threshold
   real(8) icimr                           ! in cloud ice mixing ratio
   real(8) rhdif                           ! working variable for slingo scheme

   real(8) :: rhmaxi
   real(8) :: rhmini
   real(8) :: rhminl
   real(8) :: rhminl_adj_land
   real(8) :: rhminh

   ! Statement functions
   logical land
   real(8) :: rair
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real(8), intent(in) :: premib 
    real(8), intent(in) :: premit 
    integer, intent(in) :: iceopt 
    real(8), intent(in) :: icecrit 
!real(8),  parameter :: qist_min     = 1.d-7      ! Minimum in-stratus ice IWC constraint [ kg/kg ]
!real(8),  parameter :: qist_max     = 5.d-3      ! Maximum in-stratus ice IWC constraint [ kg/kg ]
real(8) :: qist_min    ! Minimum in-stratus ice IWC constraint [ kg/kg ]
real(8) :: qist_max    ! Maximum in-stratus ice IWC constraint [ kg/kg ]
real(8) :: tmp1
!, tmp2,  tmp3
qist_min     = 1.d-7
qist_max     = 5.d-3 
!xukai-------------------------------------------------------------------------------------
   rair = wv_para(10)
   land = nint(landfrac) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87d0
     b = 0.569d0
     c = 0.002892d0
   ! Schiller parameters ( Option.2 )
     as = -68.4202d0
     bs = 0.983917d0
     cs = 2.81795d0
   ! Wood and Field parameters ( Option.3 )
     Kc = 75.d0
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.d-12
     mincld = 1.d-4
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !rhmaxi = rhmaxi_const
   !rhminl          = rhminl_const
   !rhminl_adj_land = rhminl_adj_land_const
   !rhminh          = rhminh_const
!xukai-------------------------------------------------------------------------------------
 
   !rhmaxi = rhmaxi_const
   if (present(rhmaxi_in)) rhmaxi = rhmaxi_in
   !rhmini = rhmini_const
   if (present(rhmini_in)) rhmini = rhmini_in
   !rhminl = rhminl_const
   if (present(rhminl_in)) rhminl = rhminl_in
   !rhminl_adj_land = rhminl_adj_land_const
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
   !rhminh = rhminh_const
   if (present(rhminh_in)) rhminh = rhminh_in

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     call qsat_water(T, p, es, qs, wv_para)
     esl = svp_water(T, wv_para)
     esi = svp_ice(T, wv_para)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195.d0,min(T,253.d0)) - 273.16d0
             !icicval = a + b * ttmp + c * pow_c(ttmp, 2.d0)
             tmp1 = 2.d0
             icicval = a + b * ttmp + c * pow_c(ttmp, tmp1)
             rho = p/(rair*T)
             icicval = icicval * 1.d-6 / rho 
         else
             ttmp = max(190.d0,min(T,273.16d0))
             !icicval = pow_c(10.d0 , (as * pow_c(bs, ttmp) + cs))
             tmp1 = 10.d0
             icicval = pow_c(tmp1 , (as * pow_c(bs, ttmp) + cs))
             icicval = icicval * 1.d-6 * 18.d0 / 28.97d0
         endif
         aist =  max(0.d0,min(qi/icicval,1.d0)) 
     elseif( iceopt.eq.3 ) then
         aist = 1.d0 - exp_c(-Kc*qi/(qs*(esi/esl)))
         aist = max(0.d0,min(aist,1.d0))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land .and. (snowh.le.0.000001d0) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land .and. (snowh.le.0.000001d0) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0d0-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0d0-rhwght)
           ! endif
         endif
         ncf = qi/((1.d0 - icecrit)*qs)
         if( ncf.le.0.d0 ) then 
             aist = 0.d0
         elseif( ncf.gt.0.d0 .and. ncf.le.1.d0/6.d0 ) then 
             !aist = 0.5d0*pow_c((6.d0 * ncf), (2.d0/3.d0))
             tmp1 = 2.d0/3.d0
             aist = 0.5d0*pow_c((6.d0 * ncf), tmp1)
         elseif( ncf.gt.1.d0/6.d0 .and. ncf.lt.1.d0 ) then
             !phi = (acos_c(pow_c(3.d0*(1.d0-ncf)/2.d0, (3.d0/2.d0)))+4.d0*3.1415927d0)/3.d0
             tmp1 = 3.d0/2.d0
             phi = (acos_c(pow_c(3.d0*(1.d0-ncf)/2.d0, tmp1))+4.d0*3.1415927d0)/3.d0
             aist = (1.d0 - 4.d0 * cos_c(phi) * cos_c(phi))
         else
             aist = 1.d0
         endif
             aist = max(0.d0,min(aist,1.d0))
     elseif (iceopt.eq.5) then 
        ! set rh ice cloud fraction
        rhi= (qv+qi)/qs * (esl/esi)
        rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
        !aist = min(1.0d0, pow_c(max(rhdif,0.d0), 2.d0))
        tmp1 = 2.d0
        aist = min(1.0d0, pow_c(max(rhdif,0.d0), tmp1))


        ! limiter to remove empty cloud and ice with no cloud
        ! and set icecld fraction to mincld if ice exists

        if (qi.lt.minice) then
           aist=0.d0
        else
           aist=max(mincld,aist)
        endif

        ! enforce limits on icimr
        if (qi.ge.minice) then
           icimr=qi/aist

           !minimum
           if (icimr.lt.qist_min) then
              aist = max(0.d0,min(1.d0,qi/qist_min))
           endif
           !maximum
           if (icimr.gt.qist_max) then
              aist = max(0.d0,min(1.d0,qi/qist_max))
           endif

        endif
     endif 

   ! 0.999d0 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0.d0,min(aist,0.999d0))

end subroutine aist_single

!================================================================================================

subroutine aist_vector(qv_in, T_in, p_in, qi_in, ni_in, landfrac_in, snowh_in, aist_out, ncol, &
                        wv_para, premib, premit, icecrit, iceopt,&
                       rhmaxi_in, rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   implicit none
   real(8), intent(in)  :: qv_in       ! Grid-mean water vapor[kg/kg]
   real(8), intent(in)  :: T_in        ! Temperature
   real(8), intent(in)  :: p_in        ! Pressure [Pa]
   real(8), intent(in)  :: qi_in       ! Grid-mean ice water content [kg/kg]
   real(8), intent(in)  :: ni_in       ! Grid-mean ice water number concentration [#/kg]
   real(8), intent(in)  :: landfrac_in ! Land fraction
   real(8), intent(in)  :: snowh_in    ! Snow depth (liquid water equivalent)

   real(8), intent(out) :: aist_out    ! Non-physical ice stratus fraction ( 0<= aist <= 1 )
   real(8), intent(in)  :: wv_para(PARASIZE)
   integer,  intent(in)  :: ncol 

   real(8), optional, intent(in)  :: rhmaxi_in
   real(8), optional, intent(in)  :: rhmini_in          ! Critical relative humidity for               ice stratus
   real(8), optional, intent(in)  :: rhminl_in          ! Critical relative humidity for low-level  liquid stratus
   real(8), optional, intent(in)  :: rhminl_adj_land_in ! Adjustment drop of rhminl over the land
   real(8), optional, intent(in)  :: rhminh_in          ! Critical relative humidity for high-level liquid stratus

   ! Local variables

   real(8) qv                              ! Grid-mean water vapor[kg/kg]
   real(8) T                               ! Temperature
   real(8) p                               ! Pressure [Pa]
   real(8) qi                              ! Grid-mean ice water content [kg/kg]
   real(8) ni
   real(8) landfrac                        ! Land fraction
   real(8) snowh                           ! Snow depth (liquid water equivalent)

   real(8) rhmaxi                          ! Critical relative humidity for               ice stratus
   real(8) rhmini                          ! Critical relative humidity for               ice stratus
   real(8) rhminl                          ! Critical relative humidity for low-level  liquid stratus
   real(8) rhminl_adj_land                 ! Adjustment drop of rhminl over the land
   real(8) rhminh                          ! Critical relative humidity for high-level liquid stratus

   real(8) aist                            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   real(8) rhmin                           ! Critical RH
   real(8) rhwght

   real(8) a,b,c,as,bs,cs,ah,bh,ch         ! Fit parameters
   real(8) nil
   real(8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(8) ttmp                            ! Limited temperature
   real(8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(8) rho                             ! Local air density
   real(8) esl                             ! Liq sat vapor pressure
   real(8) esi                             ! Ice sat vapor pressure
   real(8) ncf,phi                         ! Wilson and Ballard parameters
   real(8) qs
   real(8) esat_in
   real(8) qsat_in

   real(8) rhi                             ! grid box averaged relative humidity over ice
   real(8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(8) mincld                          ! minimum ice cloud fraction threshold
   real(8) icimr                           ! in cloud ice mixing ratio
   real(8) rhdif                           ! working variable for slingo scheme

   integer i


   ! Statement functions
   logical land
   real(8) :: rair
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real(8), intent(in) :: premib 
    real(8), intent(in) :: premit 
    integer, intent(in) :: iceopt 
    real(8), intent(in) :: icecrit 
!real(8),  parameter :: qist_min     = 1.d-7      ! Minimum in-stratus ice IWC constraint [ kg/kg ]
!real(8),  parameter :: qist_max     = 5.d-3      ! Maximum in-stratus ice IWC constraint [ kg/kg ]
real(8) :: qist_min       ! Minimum in-stratus ice IWC constraint [ kg/kg ]
real(8) :: qist_max       ! Maximum in-stratus ice IWC constraint [ kg/kg ]
real(8) :: tmp1, tmp2, tmp3
qist_min     = 1.d-7
qist_max     = 5.d-3
!xukai-------------------------------------------------------------------------------------
 
   rair = wv_para(10)
   land = nint(landfrac_in) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87d0
     b = 0.569d0
     c = 0.002892d0
   ! Schiller parameters ( Option.2 )
     as = -68.4202d0
     bs = 0.983917d0
     cs = 2.81795d0
   ! Wood and Field parameters ( Option.3 )
     Kc = 75.d0
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.d-12
     mincld = 1.d-4
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !rhmaxi = rhmaxi_const
   !rhminl          = rhminl_const
   !rhminl_adj_land = rhminl_adj_land_const
   !rhminh          = rhminh_const
!xukai-------------------------------------------------------------------------------------
 
     !rhmaxi          = rhmaxi_const
     if (present(rhmaxi_in))          rhmaxi          = rhmaxi_in

     !rhmini          = rhmini_const
     !rhminl          = rhminl_const
     !rhminl_adj_land = rhminl_adj_land_const
     !rhminh          = rhminh_const

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     aist_out = 0.d0
     esat_in  = 0.d0
     qsat_in  = 0.d0

     call qsat_water(T_in, p_in, &
          esat_in, qsat_in, wv_para)
     !call lwpf_stop_test_p6()
!    call set_wvflag(4)     
!    call put_tmp1(4, esat_in)
!    call put_tmp2(4, qsat_in)
!    call set_wvflag(0)     
     !do i = 1, ncol

     landfrac = landfrac_in     
     snowh = snowh_in   
     T = T_in
     qv = qv_in
     p = p_in
     qi = qi_in
     ni = ni_in
     qs = qsat_in
     esl = svp_water(T, wv_para)
     esi = svp_ice(T, wv_para)
     

 
     if (present(rhmini_in))          rhmini          = rhmini_in      
     if (present(rhminl_in))          rhminl          = rhminl_in      
     if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
     if (present(rhminh_in))          rhminh          = rhminh_in

!     call lwpf_start_test_p6()
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195.d0,min(T,253.d0)) - 273.16d0
             !icicval = a + b * ttmp + c * pow_c(ttmp, 2.d0)
             tmp1 = 2.d0
             icicval = a + b * ttmp + c * pow_c(ttmp, tmp1)
             rho = p/(rair*T)
             icicval = icicval * 1.d-6 / rho 
         else
             ttmp = max(190.d0,min(T,273.16d0))
             !icicval = pow_c(10.d0, (as * pow_c(bs, ttmp) + cs))
             tmp1 = 10.d0
             icicval = pow_c(tmp1, (as * pow_c(bs, ttmp) + cs))
             icicval = icicval * 1.d-6 * 18.d0 / 28.97d0
         endif
         aist =  max(0.d0,min(qi/icicval,1.d0)) 
     elseif( iceopt.eq.3 ) then
         aist = 1.d0 - exp_c(-Kc*qi/(qs*(esi/esl)))
         aist = max(0.d0,min(aist,1.d0))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land .and. (snowh.le.0.000001d0) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land .and. (snowh.le.0.000001d0) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0d0-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0d0-rhwght)
           ! endif
         endif
         ncf = qi/((1.d0 - icecrit)*qs)
         if( ncf.le.0.d0 ) then 
             aist = 0.d0
         elseif( ncf.gt.0.d0 .and. ncf.le.1.d0/6.d0 ) then 
             !aist = 0.5d0*pow_c((6.d0 * ncf), (2.d0/3.d0))
            tmp1 = 2.d0/3.d0
             aist = 0.5d0*pow_c((6.d0 * ncf), (2.d0/3.d0))
         elseif( ncf.gt.1.d0/6.d0 .and. ncf.lt.1.d0 ) then
             !phi = (acos_c(3.d0*(1.d0-ncf)/pow_c(2.d0, (3.d0/2.d0)))+4.d0*3.1415927d0)/3.d0
             tmp1 = 2.d0
             tmp2 = 3.d0/2.d0
             phi = (acos_c(3.d0*(1.d0-ncf)/pow_c(tmp1, tmp2))+4.d0*3.1415927d0)/3.d0
             aist = (1.d0 - 4.d0 * cos_c(phi) * cos_c(phi))
         else
             aist = 1.d0
         endif
             aist = max(0.d0,min(aist,1.d0))
     elseif (iceopt.eq.5) then 
        ! set rh ice cloud fraction
        rhi= (qv+qi)/qs * (esl/esi)
        rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
        !aist = min(1.0d0, pow_c(max(rhdif,0.d0), 2.d0))
        tmp1 = 2.d0
        aist = min(1.0d0, pow_c(max(rhdif,0.d0), tmp1))
        !tmp2 = max(rhdif,0.d0)
        !if( rhdif > 0.d0) then
        !    tmp2 = rhdif
        !else
        !    tmp2 = 0.d0
        !endif
     !call addr('tmp1', tmp1)
     !call addr('tmp2', tmp2)
      !call printf_double_f(tmp1)
      !call printf_double_f(tmp2)
     !call lwpf_start_test_p6()
     !   tmp3 = pow_c(tmp2, tmp1)
     !call lwpf_stop_test_p6()
      !call printf_double_f(tmp3)
      !  aist = min(1.0d0, tmp3)

     elseif (iceopt.eq.6) then
        !----- ICE CLOUD OPTION 6: fit based on T and Number (Gettelman: based on Heymsfield obs)
        ! Use observations from Heymsfield et al 2012 of IWC and Ni v. Temp
        ! Multivariate fit follows form of Boudala 2002: ICIWC = a * exp_c(b*T) * N^c
        ! a=6.73e-8, b=0.05, c=0.349
        ! N is #/L, so need to convert Ni_L=N*rhoa/1000.
        ah= 6.73834d-08
        bh= 0.0533110d0
        ch= 0.3493813d0
        rho=p/(rair*T)
        nil=ni*rho/1000.d0
        icicval = ah * exp_c(bh*T) * pow_c(nil, ch)
        !result is in g m-3, convert to kg H2O / kg air (icimr...)
        icicval = icicval / rho / 1000.d0
        aist =  max(0.d0,min(qi/icicval,1.d0))
        aist =  min(aist,1.d0)

     endif     

     if (iceopt.eq.5 .or. iceopt.eq.6) then


     !call lwpf_start_test_p7()
        ! limiter to remove empty cloud and ice with no cloud
        ! and set icecld fraction to mincld if ice exists

        if (qi.lt.minice) then
           aist=0.d0
        else
           aist=max(mincld,aist)
        endif

        ! enforce limits on icimr
        if (qi.ge.minice) then
           icimr=qi/aist

           !minimum
           if (icimr.lt.qist_min) then
              aist = max(0.d0,min(1.d0,qi/qist_min))
           endif
           !maximum
           if (icimr.gt.qist_max) then
              aist = max(0.d0,min(1.d0,qi/qist_max))
           endif

        endif

     !call lwpf_stop_test_p7()
     endif 

     !call lwpf_stop_test_p6()
   ! 0.999d0 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0.d0,min(aist,0.999d0))

     aist_out = aist

    !call set_wvflag(4)     
    !call put_tmp1(4, aist)
    !call put_tmp2(4, aist)
    !call set_wvflag(0)     
     !enddo

end subroutine aist_vector

!================================================================================================

end module cldfrc2m_cpe
