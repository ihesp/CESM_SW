!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module vmix_kpp

!BOP
! !MODULE: vmix_kpp
!
! !DESCRIPTION:
!  This module contains routines for initializing and computing
!  vertical mixing coefficients for the KPP parameterization
!  (see Large, McWilliams and Doney, Reviews of Geophysics, 32, 363 
!  November 1994).
!
! !REVISION HISTORY:
!  SVN:$Id: vmix_kpp.F90 65408 2014-11-13 20:08:35Z mlevy@ucar.edu $

! !USES:

   use kinds_mod
   use blocks
   use distribution
   use domain_size
   use domain
   use constants
   use grid
   use broadcast
   use io
   use state_mod
! state in kpp need state_mod be zengyh:
!  state, state_itype, state_type_mwjf, state_range_iopt, &
!  state_range_enforce, state_range_ignore, tmax, tmin, smax, smin, pressz

   use exit_mod
   use sw_absorption
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
   use io_types, only: stdout
   use communicate, only: my_task, master_task
   use niw_mixing
   use tidal_mixing, only: TIDAL_COEF, tidal_mix_max, ltidal_mixing
   use registry
   use prognostic
   use time_management
   use cvmix_kinds_and_types
   use cvmix_put_get
   use cvmix_kpp
   use cvmix_math

   implicit none
   private


   save


! !PUBLIC MEMBER FUNCTIONS:
   public :: init_vmix_kpp,   &
             vmix_coeffs_kpp, &
             add_kpp_sources, &
             smooth_hblt,     &
             linertial

! !PUBLIC DATA MEMBERS:

   real (r8), dimension(:,:,:), allocatable, public :: & 
      HMXL,               &! mixed layer depth
      KPP_HBLT,           &! boundary layer depth
      BOLUS_SP             ! scaled eddy-induced (bolus) speed used in inertial
                           !  mixing parameterization

   real (r8), public ::   &
      bckgrnd_vdc2         ! variation in diffusivity

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  mixing constants
!
!-----------------------------------------------------------------------

   real (r8) :: &
      rich_mix         ! coefficient for rich number term

   real (r8), dimension(:,:,:,:), allocatable :: &
      bckgrnd_vvc,    &! background value for viscosity
      bckgrnd_vdc      ! background value for diffusivity

   logical (log_kind) :: &
      lrich,             &! flag for computing Ri-dependent mixing
      ldbl_diff,         &! flag for computing double-diffusive mixing
      lshort_wave,       &! flag for computing short-wave forcing
      lcheckekmo,        &! check Ekman, Monin-Obhukov depth limit
      llangmuir,         &! flag for using Langmuir parameterization
      linertial,         &! flag for using inertial mixing parameterization
      lcvmix,            &! flag for using CVMix KPP routine instead of POP
      lccsm_control_compatible !flag for backwards compatibility with ccsm4 control

   integer (int_kind) :: & 
      num_v_smooth_Ri     ! num of times to vertically smooth Ri

   real (r8), parameter :: &
      epssfc = 0.1_r8       ! non-dimensional extent of sfc layer

   real (r8) ::           &
      Prandtl              ! Prandtl number

   real (r8), dimension(:,:,:), allocatable :: &
      FSTOKES        ! ratio of stokes velocity to ustar used in Langmuir 
                     !  parameterization

!-----------------------------------------------------------------------
!
!  non-local mixing (counter-gradient mixing), treated as source term
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:,:), allocatable :: & 
      KPP_SRC              ! non-local mixing (treated as source term)

!-----------------------------------------------------------------------
!
!  parameters for subroutine bldepth: computes bndy layer depth 
!
!  concv   = ratio of interior buoyancy frequency to 
!            buoyancy frequency at entrainment depth
!            parameter statement sets the minimum value.
!
!-----------------------------------------------------------------------

  real (r8), dimension(km) ::  &
      Ricr                   ! crit Rich no. for bndy layer depth
                             ! as a function of vertical resolution

   real (r8), parameter :: &
      cekman = 0.7_r8,      &! coefficient for Ekman depth
      cmonob = 1.0_r8,      &! coefficient for Monin-Obukhov depth
      concv  = 1.7_r8,      &! minimum allowed value
      hbf    = 1.0_r8        ! frac of layer affected by sol forcing

!-----------------------------------------------------------------------
!
!  parameters for subroutine ri_iwmix which computes
!  vertical mixing coefficients below boundary layer due to shear
!  instabilities, internal waves and convection
!
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      Riinfty = 0.8_r8,    &! Rich. no. limit for shear instab.
      BVSQcon = c0          ! Brunt-Vaisala square cutoff(s**-2)

!-----------------------------------------------------------------------
!
!  parameters for subroutine ddmix (double-diffusive mixing)
!
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      Rrho0  = 2.55_r8,     &! limit for double-diff density ratio
      dsfmax = 1.0_r8        ! max diffusivity for salt fingering

!-----------------------------------------------------------------------
!
!  parameters for subroutine blmix: mixing within boundary layer
!
!  cstar   = proportionality coefficient for nonlocal transport
!  cg      = non-dimensional coefficient for counter-gradient term
!
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      cstar = 10.0_r8       ! coeff for nonlocal transport

   real (r8) :: &
      cg,       &! coefficient for counter-gradient term
      Vtc        ! resolution and buoyancy independent part of the
                 ! turbulent velocity shear coefficient (for bulk Ri no)

!-----------------------------------------------------------------------
!
!  parameters for velocity scale function (from Large et al.)
!
!-----------------------------------------------------------------------

   real (r8), parameter ::   &
      zeta_m = -0.2_r8,      &
      zeta_s = -1.0_r8,      &
      c_m    =  8.38_r8,     &
      c_s    =  98.96_r8,    &
      a_m    =  1.26_r8,     &
      a_s    = -28.86_r8

!-----------------------------------------------------------------------
!
!  common vertical grid arrays used by KPP mixing
!
!-----------------------------------------------------------------------

   real (r8), dimension(:), allocatable :: & 
      zgrid,               &! depth at cell interfaces
      hwide                 ! layer thickness at interfaces

!-----------------------------------------------------------------------
!
!  ids for tavg diagnostics computed in this module
!
!-----------------------------------------------------------------------


   integer (int_kind) ::   &
      tavg_QSW_HBL,        &! tavg id for solar short-wave heat flux in bndry layer
      tavg_VDC_BCK,        &! tavg id for bckgrnd vertical tracer diffusivity
      tavg_VVC_BCK,        &! tavg id for bckgrnd vertical momentum viscosity
      tavg_KVMIX,          &! tavg id for tidal+bckgrnd vertical tracer diffusivity
      tavg_KVMIX_M,        &! tavg id for tidal+bckgrnd vertical momentum viscosity
      tavg_TPOWER

   integer (int_kind), dimension(nt) :: &
      tavg_KPP_SRC          ! tavg id for KPP_SRC for each tracer

   real (r8), dimension(:,:,:,:), allocatable :: &
      TIDAL_DIFF            ! diffusivity due to tidal mixing 

   integer (int_kind) ::   &
      tavg_KE_BL,    &! tavg id for boundary layer kinetic energy at mix time
      tavg_En         ! tavg id for boundary layer kinetic energy En



!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_vmix_kpp
! !INTERFACE:

 subroutine init_vmix_kpp(CVmix_vars, VDC, VVC)

! !DESCRIPTION:
!  Initializes constants and reads options for the KPP parameterization.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type(cvmix_data_type), dimension(:,:), intent(inout) :: CVmix_vars

   real (r8), dimension(:,:,:,:), intent(inout) :: &
      VVC        ! viscosity for momentum diffusion

   real (r8), dimension(:,:,0:,:,:),intent(inout) :: &
      VDC        ! diffusivity for tracer diffusion

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables and namelist inputs
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      n,                  &! local dummy index for tracer
      k,                  &! local dummy index for vertical lvl
      i, j, iblock,       &! local dummy indexes
      nml_error            ! namelist i/o error flag

   integer (int_kind) :: bid, jny, inx, ic, nlev

   real (r8), dimension(:), allocatable :: tmp_array

   real (r8) ::           &
      bckgrnd_vdc1,       &! background diffusivity (Ledwell)  
      bckgrnd_vdc_eq,     &! equatorial diffusivity (Gregg)
      bckgrnd_vdc_psim,   &! Max. PSI induced diffusivity (MacKinnon)
      bckgrnd_vdc_psin,   &! PSI diffusivity in northern hemisphere
      bckgrnd_vdc_psis,   &! PSI diffusivity in southern hemisphere
      bckgrnd_vdc_ban,    &! Banda Sea diffusivity (Gordon)
      bckgrnd_vdc_dpth,   &! depth at which diff equals vdc1
      bckgrnd_vdc_linv     ! inverse length for transition region

   logical (log_kind) ::  &
      lhoriz_varying_bckgrnd,  &
      larctic_bckgrnd_vdc  ! historically only used as a suboption of lniw_mixing

   namelist /vmix_kpp_nml/bckgrnd_vdc1, bckgrnd_vdc2,           &
                          bckgrnd_vdc_eq, bckgrnd_vdc_psim,     &
                          bckgrnd_vdc_ban,                      &
                          bckgrnd_vdc_dpth, bckgrnd_vdc_linv,   &
                          Prandtl, rich_mix,                    &
                          num_v_smooth_Ri, lrich, ldbl_diff,    &
                          lshort_wave, lcheckekmo,              &
                          larctic_bckgrnd_vdc,                  &
                          lhoriz_varying_bckgrnd, llangmuir,    &
                          linertial, lcvmix

   character (16), parameter :: &
      fmt_real = '(a30,2x,1pe12.5)'

   character (11), parameter :: &
      fmt_log  = '(a30,2x,l7)'

   character (11), parameter :: &
      fmt_int  = '(a30,2x,i5)'

   character (char_len) :: &
      string, string2

!-----------------------------------------------------------------------
!
!  set defaults for mixing coefficients, then read them from namelist
!
!-----------------------------------------------------------------------

   bckgrnd_vdc1           = 0.1_r8
   bckgrnd_vdc2           = c0
   bckgrnd_vdc_eq         = 0.01_r8
   bckgrnd_vdc_psim       = 0.13_r8
   bckgrnd_vdc_ban        = c1
   bckgrnd_vdc_dpth       = 2500.0e02_r8
   bckgrnd_vdc_linv       = 4.5e-05_r8
   Prandtl                = 10.0_r8
   rich_mix               = 50.0_r8
   lrich                  = .true.
   ldbl_diff              = .false.
   lshort_wave            = .false.
   lcheckekmo             = .false.
   larctic_bckgrnd_vdc    = .false.
   lhoriz_varying_bckgrnd = .false.
   llangmuir              = .false.
   linertial              = .false.
   lcvmix                 = .false.
   num_v_smooth_Ri        = 1

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=vmix_kpp_nml, iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR (init_vmix_kpp) reading vmix_kpp_nml')
   endif

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' vmix_kpp_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,vmix_kpp_nml)
      write(stdout,blank_fmt)

      write(stdout,fmt_real) '  bckgrnd_vdc1              =', bckgrnd_vdc1
      write(stdout,fmt_real) '  bckgrnd_vdc2              =', bckgrnd_vdc2
      write(stdout,fmt_real) '  bckgrnd_vdc_dpth          =', bckgrnd_vdc_dpth
      write(stdout,fmt_real) '  bckgrnd_vdc_linv          =', bckgrnd_vdc_linv
      write(stdout,fmt_real) '  bckgrnd_vdc_eq            =', bckgrnd_vdc_eq
      write(stdout,fmt_real) '  bckgrnd_vdc_psim          =', bckgrnd_vdc_psim
      write(stdout,fmt_real) '  bckgrnd_vdc_ban           =', bckgrnd_vdc_ban
      write(stdout,fmt_real) '  Prandtl                   =', Prandtl
      write(stdout,fmt_real) '  rich_mix                  =', rich_mix
      write(stdout,fmt_log ) '  Ri mixing                 =', lrich
      write(stdout,fmt_log ) '  double-diff               =', ldbl_diff
      write(stdout,fmt_log ) '  short_wave                =', lshort_wave
      write(stdout,fmt_log ) '  lcheckekmo                =', lcheckekmo
      write(stdout,fmt_int ) '  num_smooth_Ri             =', num_v_smooth_Ri
      write(stdout,fmt_log ) '  larctic_bckgrnd_vdc       =', larctic_bckgrnd_vdc
      write(stdout,fmt_log ) '  lhoriz_varying_bckgrnd    =', lhoriz_varying_bckgrnd
      write(stdout,fmt_log ) '  langmuir parameterization =', llangmuir
      write(stdout,fmt_log ) '  inertial mixing param.    =', linertial
      write(stdout,fmt_log ) '  use CVMix KPP routines    =', lcvmix
   endif

   call broadcast_scalar(bckgrnd_vdc1,          master_task)
   call broadcast_scalar(bckgrnd_vdc2,          master_task)
   call broadcast_scalar(bckgrnd_vdc_dpth,      master_task)
   call broadcast_scalar(bckgrnd_vdc_linv,      master_task)
   call broadcast_scalar(bckgrnd_vdc_eq  ,      master_task)
   call broadcast_scalar(bckgrnd_vdc_psim,      master_task)
   call broadcast_scalar(bckgrnd_vdc_ban ,      master_task)
   call broadcast_scalar(Prandtl,               master_task)
   call broadcast_scalar(rich_mix,              master_task)
   call broadcast_scalar(num_v_smooth_Ri,       master_task)
   call broadcast_scalar(lrich,                 master_task)
   call broadcast_scalar(ldbl_diff,             master_task)
   call broadcast_scalar(lshort_wave,           master_task)
   call broadcast_scalar(lcheckekmo,            master_task)
   call broadcast_scalar(larctic_bckgrnd_vdc,   master_task)
   call broadcast_scalar(lhoriz_varying_bckgrnd,master_task)
   call broadcast_scalar(llangmuir,             master_task)
   call broadcast_scalar(linertial,             master_task)
   call broadcast_scalar(lcvmix,                master_task)

!-----------------------------------------------------------------------
!
!  determine if this case must be backwards compatible with ccsm4 control
!
!-----------------------------------------------------------------------

   lccsm_control_compatible = registry_match('lccsm_control_compatible')

!-----------------------------------------------------------------------
!
!  define some non-dimensional constants 
!
!-----------------------------------------------------------------------

   Vtc     = sqrt(0.2_r8/c_s/epssfc)/vonkar**2
   cg      = cstar*vonkar*(c_s*vonkar*epssfc)**p33

!-----------------------------------------------------------------------
!
!  define vertical grid coordinates and cell widths
!  compute horizontally or vertically varying background 
!  (internal wave) diffusivity and viscosity
!
!  the vertical profile has the functional form
!
!  BCKGRND_VDC(z) = vdc1 + vdc2*atan((|z| - dpth)/L) or
!                 = vdc1 + vdc2*atan((|z| - dpth)*linv)
!
!    where
!
!  vdc1 = vertical diffusivity at |z|=D              (cm^2/s)
!  vdc2 = amplitude of variation                     (cm^2/s)
!  linv = inverse length scale for transition region (1/cm)
!  dpth = depth where diffusivity is vdc1            (cm)
!
!  the viscosity has the same form but multiplied by a constant
!  Prandtl number
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  if tidal mixing is on, use vertically constant (internal wave)
!  background diffusivity and viscosity (see consistency testing in POP_check)
!
!-----------------------------------------------------------------------

   if (bckgrnd_vdc2 /= c0 .and. my_task == master_task) then
      write (stdout,blank_fmt)
      write (stdout,'(a43)') &
            'Vertical Profile for background diffusivity'
   endif

!-----------------------------------------------------------------------
!
!  error checking
!
!-----------------------------------------------------------------------

   if (lhoriz_varying_bckgrnd .and. bckgrnd_vdc2 /= c0) then
      call exit_POP (sigAbort,  &
       'ERROR (init_vmix_kpp): lhoriz_varying_bckgrnd .and. bckgrnd_vdc2 /= c0')
   endif

!-----------------------------------------------------------------------
!
!  initialize grid info (only need one block since the vertical grid is
!  the same across blocks)
!
!-----------------------------------------------------------------------

   allocate  (zgrid(0:km+1), hwide(0:km+1))

   zgrid(0) = eps
   hwide(0) = eps
   do k=1,km
      zgrid(k) = -zt(k)
      hwide(k) =  dz(k)
   enddo
   zgrid(km+1) = -zw(km)
   hwide(km+1) = eps

   allocate  (bckgrnd_vvc(nx_block,ny_block,km,nblocks_clinic))
   allocate  (bckgrnd_vdc(nx_block,ny_block,km,nblocks_clinic))

   if (lhoriz_varying_bckgrnd) then
   
     k = 1
     do iblock=1,nblocks_clinic
     do j=1,ny_block
     do i=1,nx_block

      bckgrnd_vdc_psis= bckgrnd_vdc_psim*exp(-(0.4_r8*(TLATD(i,j,iblock)+28.9_r8))**c2)  
      bckgrnd_vdc_psin= bckgrnd_vdc_psim*exp(-(0.4_r8*(TLATD(i,j,iblock)-28.9_r8))**c2)  

      bckgrnd_vdc(i,j,k,iblock)=bckgrnd_vdc_eq+bckgrnd_vdc_psin+bckgrnd_vdc_psis

      if ( TLATD(i,j,iblock) .lt. -10.0_r8 ) then
        bckgrnd_vdc(i,j,k,iblock) = bckgrnd_vdc(i,j,k,iblock) + bckgrnd_vdc1
      elseif  ( TLATD(i,j,iblock) .le. 10.0_r8 ) then
        bckgrnd_vdc(i,j,k,iblock) = bckgrnd_vdc(i,j,k,iblock) +         &
        bckgrnd_vdc1 * (TLATD(i,j,iblock)/10.0_r8)**c2
      else
        bckgrnd_vdc(i,j,k,iblock)=bckgrnd_vdc(i,j,k,iblock) + bckgrnd_vdc1       
      endif
      
      !----------------
      ! North Banda Sea
      !----------------

      if ( (TLATD(i,j,iblock) .lt. -1.0_r8)  .and. (TLATD(i,j,iblock) .gt. -4.0_r8)  .and.  & 
           (TLOND(i,j,iblock) .gt. 103.0_r8) .and. (TLOND(i,j,iblock) .lt. 134.0_r8)) then
           bckgrnd_vdc(i,j,k,iblock) = bckgrnd_vdc_ban
      endif

      !-----------------
      ! Middle Banda Sea
      !-----------------

      if ( (TLATD(i,j,iblock) .le. -4.0_r8)  .and. (TLATD(i,j,iblock) .gt. -7.0_r8)  .and.  & 
           (TLOND(i,j,iblock) .gt. 106.0_r8) .and. (TLOND(i,j,iblock) .lt. 140.0_r8)) then
           bckgrnd_vdc(i,j,k,iblock) = bckgrnd_vdc_ban
      endif

      !----------------
      ! South Banda Sea
      !----------------

      if ( (TLATD(i,j,iblock) .le. -7.0_r8)  .and. (TLATD(i,j,iblock) .gt. -8.3_r8)  .and.  & 
           (TLOND(i,j,iblock) .gt. 111.0_r8) .and. (TLOND(i,j,iblock) .lt. 142.0_r8)) then
           bckgrnd_vdc(i,j,k,iblock) = bckgrnd_vdc_ban
      endif

      !----------------
      ! Arctic
      !----------------

      if (larctic_bckgrnd_vdc) then   ! historically only used as a suboption of lniw_mixing
      if (TLATD(i,j,iblock)  .ge. 70.0_r8) then
       bckgrnd_vdc(i,j,k,iblock) = bckgrnd_vdc_eq
      endif
      endif

      bckgrnd_vvc(i,j,k,iblock) = Prandtl*bckgrnd_vdc(i,j,k,iblock)

     end do ! i
     end do ! j
     end do ! iblock

     do k=2,km
      bckgrnd_vdc(:,:,k,:) = bckgrnd_vdc(:,:,1,:)
      bckgrnd_vvc(:,:,k,:) = bckgrnd_vvc(:,:,1,:)
     enddo

   else

!-----------------------------------------------------------------------
!
!  only need one block since the vertical grid is the same across blocks
!
!-----------------------------------------------------------------------
    do k=1,km
      bckgrnd_vdc(:,:,k,:) = bckgrnd_vdc1 + bckgrnd_vdc2* &
                       atan(bckgrnd_vdc_linv*       &
                            (zw(k)-bckgrnd_vdc_dpth))
      bckgrnd_vvc(:,:,k,:) = Prandtl*bckgrnd_vdc(:,:,k,:)

      if (bckgrnd_vdc2 /= c0 .and. my_task == master_task) then
        write (stdout,'(2x,e12.6)') bckgrnd_vdc(1,1,k,1)
      endif
    end do
 
   endif ! lhoriz_varying_bckgrnd

!-----------------------------------------------------------------------
!
!  compute crit Rich number for bndy layer depth as a function
!  of model vertical resolution. note that hwide is in cm.
!
!-----------------------------------------------------------------------

      Ricr(1:km) = 0.3_r8

!-----------------------------------------------------------------------
!
!  allocate and initialize KPP-specific arrays
!
!-----------------------------------------------------------------------

   allocate (HMXL     (nx_block,ny_block,nblocks_clinic), &
             KPP_HBLT (nx_block,ny_block,nblocks_clinic), &
             KPP_SRC  (nx_block,ny_block,km,nt,nblocks_clinic))

   HMXL     = c0
   KPP_HBLT = c0
   KPP_SRC  = c0
   VDC      = c0
   VVC      = c0

   if ( ltidal_mixing ) then
     allocate ( TIDAL_DIFF(nx_block,ny_block,km,nblocks_clinic) ) 
     TIDAL_DIFF = c0
   endif

!-----------------------------------------------------------------------
!
!  Set up CVMix variables
!
!-----------------------------------------------------------------------

   if ( lcvmix ) then
     call cvmix_init_kpp(lEkman=lcheckekmo,                                     &
                         lMonOb=lcheckekmo,                                     &
                         lnoDGat1=.false.,                                      &
                         surf_layer_ext = epssfc,                               &
                         interp_type2="POP",                                    &
                         MatchTechnique="MatchBoth")
     call cvmix_put_kpp("a_m", a_m)
     call cvmix_put_kpp("a_s", a_s)
     call cvmix_put_kpp("c_m", c_m)
     call cvmix_put_kpp("c_s", c_s)

     do bid=1,nblocks_clinic
       do jny=1,ny_block
         do inx=1,nx_block
           ic = (jny-1)*nx_block + inx
           nlev = KMT(inx, jny, bid)
           call cvmix_put(CVmix_vars(ic,bid), 'nlev', nlev)
           if (nlev.gt.0) then
             call cvmix_put(CVmix_vars(ic,bid), 'Mdiff', 0._r8)
             call cvmix_put(CVmix_vars(ic,bid), 'Tdiff', 0._r8)
             call cvmix_put(CVmix_vars(ic,bid), 'Sdiff', 0._r8)

             allocate(tmp_array(nlev+1))
             tmp_array = 0._r8
             call cvmix_put(CVmix_vars(ic,bid), 'Ri_bulk', tmp_array(1:nlev))
             call cvmix_put(CVmix_vars(ic,bid), 'zw', tmp_array)
             call cvmix_put(CVmix_vars(ic,bid), 'zt', -zt(1:nlev)*1e-2_r8)
             call cvmix_put(CVmix_vars(ic,bid), 'ocn_depth', zw(nlev)*1e-2_r8)
             CVmix_vars(ic,bid)%zw_iface(2:nlev+1) = -zw(1:nlev)*1e-2_r8
             deallocate(tmp_array)
           end if
         end do
       end do
     end do
   end if

!-----------------------------------------------------------------------
!
!  allocate and initialize the eddy-induced speed array
!
!-----------------------------------------------------------------------

   if ( linertial ) then

     allocate ( BOLUS_SP(nx_block,ny_block,nblocks_clinic) )

     BOLUS_SP = c0

   endif

!-----------------------------------------------------------------------
!
!  allocate and initialize ratio of stokes velocity to ustar 
!
!-----------------------------------------------------------------------

   allocate (FSTOKES(nx_block,ny_block,nblocks_clinic))

   do iblock=1,nblocks_clinic
     FSTOKES(:,:,iblock) = 11._r8 - MAX( c5*cos(c3*TLAT(:,:,iblock)) , c0 )
   enddo

!-----------------------------------------------------------------------
!
!  define tavg fields computed from vmix_kpp module routines
!
!-----------------------------------------------------------------------

   string = 'Solar Short-Wave Heat Flux in bndry layer'
   call define_tavg_field(tavg_QSW_HBL,'QSW_HBL',2,           &
                          long_name=trim(string),             &
                          units='watt/m^2', grid_loc='2110',  &
                          coordinates='TLONG TLAT time')

   if (ltidal_mixing) then
     string = 'Vertical diabatic diffusivity due to Tidal Mixing + background'
   else
     string = 'Vertical diabatic diffusivity due to background'
   endif
   call define_tavg_field(tavg_KVMIX,'KVMIX',3,               &
                          long_name=trim(string),             &
                          units='centimeter^2/s',             &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' ) 

   string = 'Vertical diabatic diffusivity due to background'
   call define_tavg_field(tavg_VDC_BCK,'VDC_BCK',3,               &
                          long_name=trim(string),             &
                          units='centimeter^2/s',             &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' ) 

   if (ltidal_mixing) then
     string = 'Vertical viscosity due to Tidal Mixing + background'
   else
     string = 'Vertical viscosity due to background'
   endif
   call define_tavg_field(tavg_KVMIX_M,'KVMIX_M',3,               &
                          long_name=trim(string),             &
                          units='centimeter^2/s',             &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' )

   string = 'Vertical viscosity due to background'
   call define_tavg_field(tavg_VVC_BCK,'VVC_BCK',3,               &
                          long_name=trim(string),             &
                          units='centimeter^2/s',             &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' )

   string = 'Energy Used by Vertical Mixing'
   call define_tavg_field(tavg_TPOWER,'TPOWER',3,             &
                          long_name=trim(string),             &
                          units='erg/s/cm^3',                 &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' ) 

   do n = 1,nt
     string  = 'KPP_SRC_' /&
               &/ trim(tracer_d(n)%short_name)
     string2 = trim(tracer_d(n)%short_name) /&
               &/ ' tendency from KPP non local mixing term'
     call define_tavg_field(tavg_KPP_SRC(n),trim(string),3,     &
                            long_name=trim(string2),            &
                            units=trim(tracer_d(n)%tend_units), &
                            scale_factor=tracer_d(n)%scale_factor,&
                            grid_loc='3111',                    &
                            coordinates  ='TLONG TLAT z_t time' ) 
   enddo

   if (lniw_mixing) then
     call define_tavg_field(tavg_KE_BL,'KE_BL',2,                 &
                            long_name='Boundary Layer KE',        &
                            units='ergs/centimeter**2',           &
                            grid_loc='2221')
     call define_tavg_field(tavg_En,'En',2,                       &
                            long_name='En for Boundary Layer KE ',&
                            units='Watts/meter^2',                &
                            grid_loc='2221')
   endif

!-----------------------------------------------------------------------
!EOC



!-----------defined by zengyh 2018-09-03----------
!   if(init_flag.eq.0) then
!      call athread_init()
!      init_flag=1
!   end if




 call flushm (stdout)

 end subroutine init_vmix_kpp








!***********************************************************************
!BOP
! !IROUTINE: vmix_coeffs_kpp
! !INTERFACE:

 subroutine vmix_coeffs_kpp(CVmix_vars, VDC, VVC, TRCR, UUU, VVV, UCUR, VCUR, &
                            RHOMIX, STF, SHF_QSW, this_block, convect_diff,   &
                            convect_visc, SMF, SMFT)

! !DESCRIPTION:
!  This is the main driver routine which calculates the vertical
!  mixing coefficients for the KPP mixing scheme as outlined in 
!  Large, McWilliams and Doney, Reviews of Geophysics, 32, 363 
!  (November 1994).  The non-local mixing is also computed here, but
!  is treated as a source term in baroclinic.
!
!----------------------------------------------------------------------------------
!  Updated late 2010/early 2011 to include near inertial wave parameterization.
!  Final code description is as follows. We use the diffusivity k to describe 
!  the order of computation.
!
!  k is diffusivity, k_w = background, k_n = near inertial wave, k_t = tidal, 
!  k_s = shear (Richardson), k_d = double diffusion, and k_c = convection. 
!  HBLT is boundary layer depth. The diffusivities in the description below
!  show total combination and limits after each routine is called. k_n_max is 
!  the maximum near inertial wave diffusivity allowed, and k_t_max is the 
!  corresponding maximum for the tidal diffusivity.
!
!    routine                                 description
!
!      buoydiff                (computes buoyancy difference with surface)
!      bldepth                    (computes boundary layer depth HBLT)
!      iw_reset                        initialize k:    k = k_w
!      niw_mix               (compute niwm diffusivity k_n below HBLT,
!                        and extend upward the top value into the BL for the
!                         matching slope condition on k required by blmix)
!                               k_w' = min(max(k_w,k_n),k_n_max)
!      ri_iwmix                 k = min((k_w' + k_t),k_t_max) + k_s
!      if( dbl_diff) ddmix      k = min((k_w' + k_t),k_t_max) + k_s + k_d
!      blmix            (computes k in boundary layer using present interior k)
!      .....                   (computes interior convective mixing)
!      .....    k = min((min(max(k_w,k_n),k_n_max) + k_t),k_t_max) + k_s + k_d + k_c
!  
!  Thus, in words, the background diffusivity k_w is initialized first. Then the
!  near inertial wave diffusivity is set, which is then combined with the background 
!  such that when k_n is greater than the background, it is chosen; otherwise, the
!  background is unchanged. Then, the tidal is evaluated and combined with the
!  modified background (including near inertial wave) in a similar manner. Finally,
!  the shear, double diffusion (if present) and the convection diffusivities are
!  added to the modified background/near inertial wave/tidal diffusivity.
!
!  Bruce P. Briegleb  April 2011
!----------------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TRCR                ! tracers at current time

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UUU,VVV,           &! velocities at mix time
      UCUR,VCUR           ! velocities at current time

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      RHOMIX              ! density at mix time

   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      STF                 ! surface forcing for all tracers

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SHF_QSW             ! short-wave forcing

   real (r8), dimension(nx_block,ny_block,2), intent(in), optional :: &
      SMF,               &! surface momentum forcing at U points
      SMFT                ! surface momentum forcing at T points
                         ! *** either one or the other (not
                         ! *** both) should be passed

   real (r8), intent(in) :: &
      convect_diff,      &! diffusivity to mimic convection
      convect_visc        ! viscosity   to mimic convection

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !INPUT/OUTPUT PARAMETERS:

   type(cvmix_data_type), dimension(:), intent(inout) :: CVmix_vars

   real (r8), dimension(nx_block,ny_block,km), intent(inout) ::      &
      VVC        ! viscosity for momentum diffusion

   real (r8), dimension(nx_block,ny_block,0:km+1,2),intent(inout) :: &
      VDC        ! diffusivity for tracer diffusion

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) ::  &
      error_string

   integer (int_kind) :: &
      k,                 &! vertical level index 
      i,j,               &! horizontal loop indices
      n,                 &! tracer index
      mt2,               &! index for separating temp from other trcrs
      bid                 ! local block address for this block

   integer (int_kind), dimension(nx_block,ny_block) :: &
      KBL                   ! index of first lvl below hbl

   real (r8), dimension(nx_block,ny_block) :: &
      USTAR,      &! surface friction velocity
      BFSFC,      &! surface buoyancy forcing
      WORK1,WORK2,&! temporary storage
      FCON,       &! convection temporary
      STABLE,     &! = 1 for stable forcing; = 0 for unstable forcing
      KE_mix,     &! kinetic energy at mix time
      KE_cur,     &! kinetic energy at cur time
      En           ! En for boundary layer kinetic energy
 
   real (r8), dimension(nx_block,ny_block,km) :: &
      DBLOC,      &! buoyancy difference between adjacent levels
      DBSFC,      &! buoyancy difference between level and surface
      GHAT         ! non-local mixing coefficient

   real (r8), dimension(nx_block,ny_block,km,2) :: &
      GHAT2         ! non-local mixing coefficient (CVMix)

   real (r8), dimension(nx_block,ny_block,0:km+1) :: &
      VISC        ! local temp for viscosity

   real (r8) ::  &
      factor

!-----------------------------------------------------------------------
!
!  initialize  and consistency checks
!
!-----------------------------------------------------------------------
! added by zengyh 2019-01-08:
   interface
     subroutine m_kpp_last_sw(KBL, KMT, KMU,            &
          dz, zt, zgrid,                                &
          AU0, AUN, AUE, AUNE, BFSFC,  HMXL,            &
          STF, GHAT, DBSFC, DBLOC, DZT,                 &
          VISC, VVC, VDC, KPP_SRC,                      &
          myproc, nx_block, ny_block, km, nt, bid,      &
          convect_visc, convect_diff, BVSQcon           &



          ) bind(C, name='m_kpp_last_sw')
         integer,value:: myproc, nx_block, ny_block, km, nt, bid
         integer:: KBL(*), KMT(*), KMU(*)
         real(8),value:: convect_visc, convect_diff, BVSQcon
         real(8):: dz(*), zt(*), zgrid(*),              &
          AU0(*), AUN(*), AUE(*), AUNE(*), BFSFC(*),  HMXL(*),   &
          STF(*), GHAT(*), DBSFC(*), DBLOC(*), DZT(*),  &
          VISC(*), VVC(*), VDC(*), KPP_SRC(*)



     end subroutine m_kpp_last_sw
   end interface

   bid = this_block%local_id

   if (.not. present(SMF) .and. .not. present(SMFT)) then
      error_string = 'ERROR KPP: must supply either SMF or SMFT'
      call document ('vmix_coeffs_kpp',  trim(error_string))
      call exit_POP(sigAbort, trim(error_string))
   endif

!-----------------------------------------------------------------------
!
!  compute buoyancy differences at each vertical level.
!
!-----------------------------------------------------------------------

   call buoydiff(DBLOC, DBSFC, TRCR, this_block)

   if (lniw_mixing) then
!-----------------------------------------------------------------------
!
!     when lniw_mixing, compute boundary layer depth now
!
!-----------------------------------------------------------------------

     if (present(SMFT)) then
        call bldepth (CVmix_vars, DBLOC, DBSFC, TRCR, UUU, VVV, UCUR, VCUR,   &
                      STF, SHF_QSW, KPP_HBLT(:,:,bid), USTAR, BFSFC, STABLE,  &
                      KBL, this_block, SMFT=SMFT)
     else
        call bldepth (CVmix_vars, DBLOC, DBSFC, TRCR, UUU, VVV, UCUR, VCUR,   &
                      STF, SHF_QSW, KPP_HBLT(:,:,bid), USTAR, BFSFC, STABLE,  &
                      KBL, this_block, SMF=SMF)
     endif


     call compute_niw_energy_flux(VISC,VDC,UUU,VVV,KE_mix,UCUR,VCUR,KE_cur,  &
                                  DBLOC, KPP_HBLT(:,:,bid),KBL,En,this_block)

   endif ! if (lniw_mixing) then

!-----------------------------------------------------------------------
!
!  compute mixing due to shear instability, internal waves and
!  convection
!
!-----------------------------------------------------------------------

   call ri_iwmix(DBLOC, VISC, VDC, UUU, VVV, RHOMIX,&
                 convect_diff, convect_visc, this_block)

!-----------------------------------------------------------------------
!
!  compute double diffusion if desired
!
!-----------------------------------------------------------------------

   if (ldbl_diff) call ddmix(VDC, TRCR, this_block)

!-----------------------------------------------------------------------
!
!     compute boundary layer depth (when no niw mixing)
!
!-----------------------------------------------------------------------

   if (.not. lniw_mixing) then
     if (present(SMFT)) then
#if defined(POPCHKJN)
         allocate(zyhKBL(nx_block,ny_block))
         allocate(tmpKBL(nx_block,ny_block))
         allocate(zyhHBLT(nx_block,ny_block))
         allocate(tmpHBLT(nx_block,ny_block))
         allocate(zyhrHBLT(nx_block,ny_block))
         allocate(zyhBFSFC(nx_block,ny_block))
         allocate(tmpBFSFC(nx_block,ny_block))
         allocate(zyhrBFSFC(nx_block,ny_block))
         allocate(zyhSTABLE(nx_block,ny_block))
         allocate(tmpSTABLE(nx_block,ny_block))
         allocate(zyhUSTAR(nx_block,ny_block))
         allocate(tmpUSTAR(nx_block,ny_block))

!         if (my_task == master_task) then
!            WRITE (*,*) 'step into call bldepth from vmix_kpp.F90:'
!         endif

         tmpKBL    = KBL
         tmpHBLT   = KPP_HBLT(:,:,bid)
         tmpBFSFC  = BFSFC
         tmpSTABLE = STABLE
         tmpUSTAR  = USTAR

       call getlooptime(start_kpp)

       call bldepth (CVmix_vars, DBLOC, DBSFC, TRCR, UUU, VVV, UCUR, VCUR,   &
                      STF, SHF_QSW, KPP_HBLT(:,:,bid), USTAR, BFSFC, STABLE,  &
                      KBL, this_block, SMFT=SMFT)

       call getlooptime(stop_kpp)

       if (my_task == master_task) then
            WRITE (*,*) 'initial kpp_bldepth costs zengyh:',stop_kpp-start_kpp,'(us);'
       endif

         zyhKBL    = KBL
         zyhHBLT   = KPP_HBLT(:,:,bid)
         zyhBFSFC  = BFSFC
         zyhSTABLE = STABLE
         zyhUSTAR  = USTAR

         KBL       = tmpKBL
         KPP_HBLT(:,:,bid) = tmpHBLT
         BFSFC     = tmpBFSFC
         STABLE    = tmpSTABLE
         USTAR     = tmpUSTAR

       call getlooptime(start_kpp)

#endif

 ! conditionally athreaded subroutine by zengyh 2018-11-01:
       if(.not. lcheckekmo .and. .not. linertial .and. .not. llangmuir &
         .and. .not. l1Ddyn .and. sw_absorption_type .eq. 'chlorophyll') then
         call bldepth_sw (CVmix_vars, DBLOC, DBSFC, TRCR, UUU, VVV, UCUR, VCUR,   &
                      STF, SHF_QSW, KPP_HBLT(:,:,bid), USTAR, BFSFC, STABLE,  &
                      KBL, this_block, SMFT=SMFT)
       else
         call bldepth (CVmix_vars, DBLOC, DBSFC, TRCR, UUU, VVV, UCUR, VCUR,   &
                      STF, SHF_QSW, KPP_HBLT(:,:,bid), USTAR, BFSFC, STABLE,  &
                      KBL, this_block, SMFT=SMFT)
       endif

#if defined(POPCHKJN)
       call getlooptime(stop_kpp)

       if (my_task == master_task) then
            WRITE (*,*) 'athreaded kpp_bldepth costs zengyh:',stop_kpp-start_kpp,'(us);'
       endif

         zyhKBL    = abs(zyhKBL    - KBL)
         zyhHBLT   = abs(zyhHBLT   - KPP_HBLT(:,:,bid))
         zyhrHBLT   = merge(zyhHBLT/KPP_HBLT(:,:,bid), c0, abs(KPP_HBLT(:,:,bid))>c0)
         zyhBFSFC  = abs(zyhBFSFC  - BFSFC)
         zyhrBFSFC  = merge(zyhBFSFC/BFSFC, c0, abs(BFSFC)>c0)
         zyhSTABLE = abs(zyhSTABLE - STABLE)
         zyhUSTAR  = abs(zyhUSTAR  - USTAR)

         zyhaveerrKBL = sum(zyhKBL)/(nx_block*ny_block)            
         zyhmaxerrKBL = maxval(zyhKBL)                
         zyhminerrKBL = minval(zyhKBL)                

         zyhaveerrHBLT = sum(zyhHBLT)/(nx_block*ny_block)            
         zyhavererrHBLT = sum(zyhrHBLT)/(nx_block*ny_block)            
         zyhmaxerrHBLT = maxval(zyhHBLT)                
         zyhmaxrerrHBLT = maxval(abs(zyhrHBLT))                
         zyhminerrHBLT = minval(zyhHBLT)                

         zyhaveerrBFSFC = sum(zyhBFSFC)/(nx_block*ny_block)            
         zyhavererrBFSFC = sum(zyhrBFSFC)/(nx_block*ny_block)            
         zyhmaxerrBFSFC = maxval(zyhBFSFC)                
         zyhmaxrerrBFSFC = maxval(abs(zyhrBFSFC))                
         zyhminerrBFSFC = minval(zyhBFSFC)                

         zyhaveerrSTABLE = sum(zyhSTABLE)/(nx_block*ny_block)            
         zyhmaxerrSTABLE = maxval(zyhSTABLE)                
         zyhminerrSTABLE = minval(zyhSTABLE)                

         zyhaveerrUSTAR = sum(zyhUSTAR)/(nx_block*ny_block)            
         zyhmaxerrUSTAR = maxval(zyhUSTAR)                
         zyhminerrUSTAR = minval(zyhUSTAR)                

         if ((my_task==0 .or. my_task==256 .or. my_task==512 .or. my_task==768) .and. (num_steps<=512 .or. mod(num_steps,64)==0)) then
            write(*,*) 'PROC: ', my_task,'; error of KBL: ', zyhmaxerrKBL,zyhaveerrKBL,zyhminerrKBL
            write(*,*) 'PROC: ', my_task,'; error of HBLT: ', zyhmaxrerrHBLT,zyhmaxerrHBLT,zyhavererrHBLT,zyhaveerrHBLT,zyhminerrHBLT
            write(*,*) 'PROC: ', my_task,'; error of BFSFC: ', zyhmaxrerrBFSFC,zyhmaxerrBFSFC,zyhavererrBFSFC,zyhaveerrBFSFC,zyhminerrBFSFC
!            write(*,*) 'error of STABLE: ', zyhmaxerrSTABLE,zyhaveerrSTABLE,zyhminerrSTABLE
!            write(*,*) 'error of USTAR: ', zyhmaxerrUSTAR,zyhaveerrUSTAR,zyhminerrUSTAR
         endif


         deallocate(zyhKBL)
         deallocate(zyhHBLT)
         deallocate(zyhrHBLT)
         deallocate(zyhBFSFC)
         deallocate(zyhrBFSFC)
         deallocate(zyhSTABLE)
         deallocate(zyhUSTAR)
#endif
     else
        call bldepth (CVmix_vars, DBLOC, DBSFC, TRCR, UUU, VVV, UCUR, VCUR,   &
                      STF, SHF_QSW, KPP_HBLT(:,:,bid), USTAR, BFSFC, STABLE,  &
                      KBL, this_block, SMF=SMF)
     endif
   endif ! .not. lniw_mixing

!-----------------------------------------------------------------------
!
!  compute boundary layer diffusivities and match to interior values
!
!-----------------------------------------------------------------------

   call blmix(CVmix_vars, VISC, VDC, KPP_HBLT(:,:,bid), USTAR, BFSFC, STABLE, &
              KBL, GHAT, GHAT2, this_block)





!-----------------------------------------------------------------------
!
!  consider interior convection:
!
!    compute function of Brunt-Vaisala squared for convection.
!
!  use either a smooth    
!
!    WORK1 = N**2,  FCON is function of N**2
!    FCON = 0 for N**2 > 0
!    FCON = [1-(1-WORK1/BVSQcon)**2]**3 for BVSQcon < N**2 < 0
!    FCON = 1 for N**2 < BVSQcon
!
!  or a step function. The smooth function has been used with
!  BVSQcon = -0.2e-4_dbl_kind.
!
!  after convection, average viscous coeffs to U-grid and reset sea 
!  floor values
!
!-----------------------------------------------------------------------

!----added by zengyh 2018-12-26:-------------
if(partial_bottom_cells .and. (.not. l1Ddyn) .and. (.not. lcvmix)) then

!initial kpp_last costs zengyh: 74046.-77272. (us);


   call m_kpp_last_sw(KBL, KMT(:,:,bid), KMU(:,:,bid),    &
          dz, zt, zgrid,                                &
          AU0(:,:,bid), AUN(:,:,bid), AUE(:,:,bid), AUNE(:,:,bid), &
          BFSFC,  HMXL(:,:,bid),                        &
          STF, GHAT, DBSFC, DBLOC, DZT(:,:,:,bid),      &
          VISC, VVC, VDC, KPP_SRC(:,:,:,:,bid),         &
          my_task, nx_block, ny_block, km, nt, bid,      &
          convect_visc, convect_diff, BVSQcon           &



          )
else


   do k=1,km-1           

      if (partial_bottom_cells) then
         WORK1 = DBLOC(:,:,k)/(p5*(DZT(:,:,k  ,bid) + &
                                   DZT(:,:,k+1,bid)))
      else
         WORK1 = DBLOC(:,:,k)/(zgrid(k) - zgrid(k+1))
      end if

      if (BVSQcon /= c0) then
         WORK2 = min(c1-(max(WORK1,BVSQcon))/BVSQcon, c1)
         FCON  = (c1 - WORK2*WORK2)**3
      else
         where (WORK1 > c0)
            FCON = c0
         elsewhere
            FCON = c1
         end where
      endif

      !*** add convection and reset sea floor values to zero
      

      do j=1,ny_block
      do i=1,nx_block
         if ( k >= KBL(i,j) ) then
            VISC(i,j,k)  = VISC(i,j,k)  + convect_visc * FCON(i,j)
            VDC(i,j,k,1) = VDC(i,j,k,1) + convect_diff * FCON(i,j)
            VDC(i,j,k,2) = VDC(i,j,k,2) + convect_diff * FCON(i,j)
         endif
         if (k >= KMT(i,j,bid)) then
         	  VISC(i,j,k  ) = c0
            VDC (i,j,k,1) = c0
            VDC (i,j,k,2) = c0
         endif
      end do
      end do

      !*** now average visc to U grid 
      if (l1Ddyn) then
        VVC(:,:,k) = merge(VISC(:,:,k), c0, (k < KMU(:,:,bid)))
      else
        call tgrid_to_ugrid(WORK2,VISC(:,:,k),bid)
        VVC(:,:,k) = merge(WORK2, c0, (k < KMU(:,:,bid)))
      end if

   enddo

   VDC(:,:,km,:) = c0
   VVC(:,:,km)   = c0

!-----------------------------------------------------------------------
!
!  add ghatp term from previous computation to right-hand-side 
!  source term on current row
!
!-----------------------------------------------------------------------

   do n=1,nt
      mt2=min(n,2)
      if (lcvmix) then
        KPP_SRC(:,:,1,n,bid) = STF(:,:,n)/dz(1)           &
                               *(-GHAT2(:,:,1,mt2))
        if (partial_bottom_cells) then
           do k=2,km
              KPP_SRC(:,:,k,n,bid) = STF(:,:,n)/DZT(:,:,k,bid)         &
                                   *( GHAT2(:,:,k-1,mt2) - GHAT2(:,:,k,mt2))
           enddo
        else
           do k=2,km
              KPP_SRC(:,:,k,n,bid) = STF(:,:,n)/dz(k)                  &
                                   *( GHAT2(:,:,k-1,mt2) - GHAT2(:,:,k,mt2))
           enddo
        endif
      else ! no CVMix
        KPP_SRC(:,:,1,n,bid) = STF(:,:,n)/dz(1)           &
                               *(-VDC(:,:,1,mt2)*GHAT(:,:,1))
        if (partial_bottom_cells) then
           do k=2,km
              KPP_SRC(:,:,k,n,bid) = STF(:,:,n)/DZT(:,:,k,bid)         &
                                   *( VDC(:,:,k-1,mt2)*GHAT(:,:,k-1)   &
                                     -VDC(:,:,k  ,mt2)*GHAT(:,:,k  ))
           enddo
        else
           do k=2,km
              KPP_SRC(:,:,k,n,bid) = STF(:,:,n)/dz(k)                  &
                                   *( VDC(:,:,k-1,mt2)*GHAT(:,:,k-1)   &
                                     -VDC(:,:,k  ,mt2)*GHAT(:,:,k  ))
           enddo
        endif
      end if
   enddo

!-----------------------------------------------------------------------
!
!  compute diagnostic mixed layer depth (cm) using a max buoyancy 
!  gradient criterion.  Use USTAR and BFSFC as temps.
!
!-----------------------------------------------------------------------

   USTAR = c0
   where (KMT(:,:,bid) == 1)
      HMXL(:,:,bid) = zt(1)
   elsewhere
      HMXL(:,:,bid) = c0
   endwhere

   if (partial_bottom_cells) then
      do k=2,km
         where (k <= KMT(:,:,bid))
            STABLE = zt(k-1) + p5*(DZT(:,:,k-1,bid) + DZT(:,:,k,bid))
            USTAR = max(DBSFC(:,:,k)/STABLE,USTAR)
            HMXL(:,:,bid) = STABLE
         endwhere
      enddo

      VISC(:,:,1) = c0
      do k=2,km
         where (USTAR > c0 )
            VISC(:,:,k) = (DBSFC(:,:,k)-DBSFC(:,:,k-1))/ &
                          (p5*(DZT(:,:,k,bid) + DZT(:,:,k-1,bid)))
         end where
         where ( VISC(:,:,k) >= USTAR .and.              &
                (VISC(:,:,k)-VISC(:,:,k-1)) /= c0 .and.  &
                 USTAR > c0 )   ! avoid divide by zero
            BFSFC = (VISC(:,:,k) - USTAR)/ &
                    (VISC(:,:,k)-VISC(:,:,k-1))
! tqian
!            HMXL(:,:,bid) =   (zt(k-1) + p5*DZT(:,:,k-1,bid))*(c1-BFSFC) &
!                            + (zt(k-1) - p5*DZT(:,:,k-1,bid))*BFSFC
             HMXL(:,:,bid) =   (zt(k-1) + p25*(DZT(:,:,k-1,bid)+DZT(:,:,k,bid)))*(c1-BFSFC) &
                             + (zt(k-1) - p25*(DZT(:,:,k-2,bid)+DZT(:,:,k-1,bid)))*BFSFC

            USTAR(:,:) = c0
         endwhere
      enddo
   else
      do k=2,km
         where (k <= KMT(:,:,bid))
            USTAR = max(DBSFC(:,:,k)/zt(k),USTAR)
            HMXL(:,:,bid) = zt(k)
         endwhere
      enddo

      VISC(:,:,1) = c0
      do k=2,km
         where (USTAR > c0 )
            VISC(:,:,k) = (DBSFC(:,:,k)-DBSFC(:,:,k-1))/ &
                          (zt(k) - zt(k-1))
         end where
         where ( VISC(:,:,k) >= USTAR .and.              &
                (VISC(:,:,k)-VISC(:,:,k-1)) /= c0 .and.  &
                 USTAR > c0 )   ! avoid divide by zero
            BFSFC = (VISC(:,:,k) - USTAR)/ &
                    (VISC(:,:,k)-VISC(:,:,k-1))
            HMXL(:,:,bid) = -p5*(zgrid(k  ) + zgrid(k-1))*(c1-BFSFC) &
                            -p5*(zgrid(k-1) + zgrid(k-2))*BFSFC
            USTAR(:,:) = c0
         endwhere
      enddo
   endif
endif
!-----------------------------------------------------------------------
!EOC

 end subroutine vmix_coeffs_kpp













!!!------- need to be athreaded, guoq.-------- 
!***********************************************************************
!BOP
! !IROUTINE: ri_iwmix
! !INTERFACE:


 subroutine ri_iwmix(DBLOC, VISC, VDC, UUU, VVV, RHOMIX, &
                     convect_diff, convect_visc, this_block)

! !DESCRIPTION:
!  Computes viscosity and diffusivity coefficients for the interior
!  ocean due to shear instability (richardson number dependent),
!  internal wave activity, and to static instability (Ri < 0).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UUU               ! U velocities at current time

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: & 
      VVV,             &! V velocities at current time
      DBLOC             ! buoyancy difference between adjacent levels

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      RHOMIX            ! density at mix time

   real (r8), intent(in) :: &
      convect_diff,         &! diffusivity to mimic convection
      convect_visc           ! viscosity   to mimic convection

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,0:km+1,2), intent(inout) :: & 
      VDC        ! diffusivity for tracer diffusion

   real (r8), dimension(nx_block,ny_block,0:km+1), intent(inout) :: & 
      VISC       ! viscosity

   logical (log_kind), parameter :: &
      prnt = .false. ! if true, diagnostic prints made of tidal, rich

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k,                 &! index for vertical levels
      i,j,               &! horizontal loop indices
      bid,               &! local block index
      n                   ! vertical smoothing index

   real (r8), dimension(nx_block,ny_block) :: &
      KVMIX, 		 &! vertical diffusivity
      KVMIX_M		  ! vertical viscosity

   real (r8), dimension(nx_block,ny_block,0:km+1) :: &
      WORK0               ! work array 

   real (r8), dimension(nx_block,ny_block) :: &
      VSHEAR,            &! (local velocity shear)^2
      RI_LOC,            &! local Richardson number 
      FRI,               &! function of Ri for shear
      FCON,              &
      WORK1,             &! local work array
      WORKN               ! local work array
  integer :: lniw,lr,par,Ddyn !guoq
!-----------------------------------------------------------------------
!
!  compute mixing at each level
!
!-----------------------------------------------------------------------
   bid = this_block%local_id

   KVMIX       = c0
   KVMIX_M     = c0
   WORK0(:,:,0)= c0

!-----------------------------------------------------------------------
!
!     compute velocity shear squared and average to T points:
!     VSHEAR = (UUU(k)-UUU(k+1))**2+(VVV(k)-VVV(k+1))**2
!     Use FRI here as a temporary.
!
!-----------------------------------------------------------------------
   if (partial_bottom_cells) then
       par=1
   else
       par=0
   endif
   if (l1Ddyn) then
       Ddyn=1
   else
       Ddyn=0
   endif
      call ri_loop01_fun(nx_block,ny_block,bid,km,par,Ddyn,p5,c0,eps,KMT,UUU,VVV,FRI,DZU,DZT,&
                      VSHEAR,DBLOC,zgrid,RI_LOC,WORK0,AT0,ATS,ATW,ATSW)
      call ri_loop01_sw()

#ifdef master 
    do k=1,km
      if (k < km) then

         FRI = (UUU(:,:,k)-UUU(:,:,k+1))**2 + &
               (VVV(:,:,k)-VVV(:,:,k+1))**2

         if (partial_bottom_cells) then
            FRI = FRI/(p5*(DZU(:,:,k  ,bid) + DZU(:,:,k+1,bid)))**2
         endif
       
         if (l1Ddyn) then   !l1Ddyn=F
           VSHEAR = FRI
         else
              
          call ugrid_to_tgrid(VSHEAR,FRI,bid)
         
      !  do j=2,ny_block
      !  do i=2,nx_block

      !      VSHEAR(i,j) =   AT0 (i,j,bid)*FRI(i  ,j  ) + &
      !                      ATS (i,j,bid)*FRI(i  ,j-1) + &
      !                      ATW (i,j,bid)*FRI(i-1,j  ) + &
      !                      ATSW(i,j,bid)*FRI(i-1,j-1)

     !   end do
     !   end do

     !      VSHEAR(:,1) = c0
     !      VSHEAR(1,:) = c0
       end if

    else

         VSHEAR = c0

    endif

!-----------------------------------------------------------------------
!
!     compute local richardson number
!
!-----------------------------------------------------------------------

      if (partial_bottom_cells) then
         if (k < km) then
            RI_LOC = DBLOC(:,:,k)/                              &
                     (VSHEAR + eps/(p5*(DZT(:,:,k  ,bid) +      &
                                        DZT(:,:,k+1,bid)))**2)/ &
                     (p5*(DZT(:,:,k,bid) + DZT(:,:,k+1,bid)))
         else
            RI_LOC = DBLOC(:,:,k)/(VSHEAR +              &
                     eps/(p5*DZT(:,:,k,bid))**2)/(p5*DZT(:,:,k,bid))
         end if
      else
         RI_LOC = DBLOC(:,:,k)*(zgrid(k)-zgrid(k+1))/(VSHEAR + eps)
      end if

      WORK0(:,:,k)   = merge(RI_LOC, WORK0(:,:,k-1), k <= KMT(:,:,bid))
 
   enddo
#endif
!-----------------------------------------------------------------------
!
!  vertically smooth Ri num_v_smooth_Ri times with 1-2-1 weighting
!  result again stored in WORK0 and use RI_LOC and FRI
!  as temps
!
!-----------------------------------------------------------------------
 
   do n = 1,num_v_smooth_Ri
      call ri_loop02_fun(nx_block,ny_block,bid,k,km,p5,p25,KMT,WORK0,FRI,RI_LOC)
      call ri_loop02_sw()

#ifdef master
 
      FRI            =  p25 * WORK0(:,:,1)
      WORK0(:,:,km+1) =       WORK0(:,:,km)
 
      do k=1,km
       !  !DIR$ NODEP
         do j=1,ny_block
       !  !DIR$ NODEP
         do i=1,nx_block
            RI_LOC(i,j) = WORK0(i,j,k)
            if (KMT(i,j,bid) >= 3) then
               WORK0(i,j,k) = FRI(i,j) + p5*RI_LOC(i,j) + p25*WORK0(i,j,k+1)
            endif
            FRI(i,j) = p25*RI_LOC(i,j)
         enddo
         enddo
      enddo
#endif
   enddo
!-----------------------------------------------------------------------
!
!  now that we have a smoothed Ri field, finish computing coeffs
!  at each level
!
!-----------------------------------------------------------------------

   if ( ltidal_mixing )  TIDAL_DIFF(:,:,:,bid) = c0   !ltidal_mixing=F   
!   do k = 1,km

!-----------------------------------------------------------------------
!
!     if Ri-number mixing requested,
!     evaluate function of Ri for shear instability:
!       for 0 < Ri < Riinfty, function = (1 - (Ri/Riinfty)**2)**3
!       for     Ri > Riinfty, function = 0
!       for     Ri < 0      , function = 1
!     compute contribution due to shear instability
!     WORK0 holds smoothed Ri at k
!
!     otherwise only use iw
!     convection is added later
!
!-----------------------------------------------------------------------
      if ( ltidal_mixing ) then   
        do k = 1,km 
!-----------------------------------------------------------------------
!
!  consider the internal wave mixing first. rich_mix is used as the
!  upper limit for internal wave mixing coefficient.
!
!  NOTE: no partial_bottom_cell implementation at this time 
!
!-----------------------------------------------------------------------

        WORK1 = DBLOC(:,:,k)/(zgrid(k) - zgrid(k+1))

        where (WORK1 > c0)
          TIDAL_DIFF(:,:,k,bid) = TIDAL_COEF(:,:,k,bid)/WORK1
        endwhere

        ! Notes:
        ! (1) this step breaks backwards compatibility
        ! (2) check for k>2 was added to if statement to avoid
        !     out of bounds access
        if ((.not. lccsm_control_compatible).and.(k.gt.2)) then
        where ( k == KMT(:,:,bid)-1  .or.  k == KMT(:,:,bid)-2 )
          TIDAL_DIFF(:,:,k,bid) = max( TIDAL_DIFF(:,:,k,  bid),  &
                                       TIDAL_DIFF(:,:,k-1,bid) )
        endwhere 
        endif

        if (lniw_mixing) then
         WORKN = VISC(:,:,k)
        else
         WORKN = bckgrnd_vvc(:,:,k,bid)
        endif

        WORK1 = Prandtl*min(WORKN(:,:)/Prandtl  &
                            + TIDAL_DIFF(:,:,k,bid), tidal_mix_max)
        if( prnt ) then             !prnt=F

! use DS product point to test, global (i,j) = (5,103)
           do j=1,ny_block
             if( this_block%j_glob(j) .eq. 103 ) then
             do i=1,nx_block
               if( this_block%i_glob(i) .eq. 5 ) then
       if( k < KMT(i,j,bid) ) then
 write(stdout,100) this_block%i_glob(i),this_block%j_glob(j),k, &
                         TIDAL_DIFF(i,j,k,bid)
 100 format(' tidal  i,j,k TIDAL_DIFF =',3(i3,1x),1pe11.4,1x)
       endif
               endif
             enddo ! i
             endif
           enddo ! j

         endif

        if ( k < km ) then
          KVMIX_M(:,:) = WORK1(:,:)
        endif

        if ( k < km ) then
          if (lniw_mixing) then
            VDC(:,:,k,2) = min(VDC(:,:,k,2) + TIDAL_DIFF(:,:,k,bid),  &
                               tidal_mix_max)
          else
            VDC(:,:,k,2) = min(bckgrnd_vdc(:,:,k,bid) + TIDAL_DIFF(:,:,k,bid),  &
                               tidal_mix_max)
          endif
          KVMIX(:,:) = VDC(:,:,k,2)
        endif
        if (lrich) then           !lrich=F
          FRI    = min((max(WORK0(:,:,k),c0))/Riinfty, c1)

          VISC(:,:,k) = WORK1 + rich_mix*(c1 - FRI*FRI)**3

          if( prnt ) then

! use DS product point to test, global (i,j) = (5,103)
            do j=1,ny_block
              if( this_block%j_glob(j) .eq. 103 ) then
              do i=1,nx_block
                if( this_block%i_glob(i) .eq. 5 ) then
        if( k < KMT(i,j,bid) ) then
  write(stdout,200) this_block%i_glob(i),this_block%j_glob(j),k, &
                          rich_mix*(c1 - FRI(i,j)*FRI(i,j))**3
  200 format(' Rich  i,j,k DIFF =',3(i3,1x),1pe11.4,1x)
        endif
                endif
              enddo ! i
              endif
            enddo ! j

          endif
          if ( k < km ) then
            VDC(:,:,k,2) = VDC(:,:,k,2) + rich_mix*(c1 - FRI*FRI)**3
            VDC(:,:,k,1) = VDC(:,:,k,2)
          endif
        else
          VISC(:,:,k) = WORK1 

          if ( k < km ) then
            VDC(:,:,k,1) = VDC(:,:,k,2)
          endif
        endif

      enddo  !guoq
   else ! .not. ltidal_mixin

   !    do k=1,km     
      if(lniw_mixing) then
        lniw=1
      else 
        lniw=0
      endif

      if(lrich) then
        lr=1
      else 
        lr=0
      endif
        
    call ri_loop03_fun(nx_block,ny_block,bid,km,lniw,lr,c0,c1,Riinfty,rich_mix,VDC, &
                          VISC,KVMIX,KVMIX_M,bckgrnd_vdc,bckgrnd_vvc,WORK0,FRI)
    call ri_loop03_sw()   
  
#ifdef master          
         do k=1,km
        if ( k < km ) then
          if (lniw_mixing) then
            KVMIX(:,:) = VDC(:,:,k,2)
            KVMIX_M(:,:) = VISC(:,:,k)
          else
            KVMIX(:,:) = bckgrnd_vdc(:,:,k,bid)
            KVMIX_M(:,:) = bckgrnd_vvc(:,:,k,bid)
          endif
        endif

        if (lniw_mixing) then
        if (lrich) then
            
           FRI    = min((max(WORK0(:,:,k),c0))/Riinfty, c1)

           VISC(:,:,k  ) = VISC(:,:,k) + &
                           rich_mix*(c1 - FRI*FRI)**3

           if ( k < km ) then
              VDC (:,:,k,2) = VDC(:,:,k,2) + &
                              rich_mix*(c1 - FRI*FRI)**3
              VDC(:,:,k,1) = VDC(:,:,k,2)
           endif
        endif
        else 
        if (lrich) then
           FRI    = min((max(WORK0(:,:,k),c0))/Riinfty, c1)

           VISC(:,:,k  ) = bckgrnd_vvc(:,:,k,bid) + &
                           rich_mix*(c1 - FRI*FRI)**3

           if ( k < km ) then
              VDC (:,:,k,2) = bckgrnd_vdc(:,:,k,bid) + &
                              rich_mix*(c1 - FRI*FRI)**3
              VDC(:,:,k,1) = VDC(:,:,k,2)
           endif
        else
           
           VISC(:,:,k  ) = bckgrnd_vvc(:,:,k,bid)

           if ( k < km ) then
              VDC (:,:,k,2) = bckgrnd_vdc(:,:,k,bid)
              VDC(:,:,k,1) = VDC(:,:,k,2)
           endif
        endif
        endif ! lniw_mixing

      enddo !guoq

#endif
    endif ! ltidal_mixing
!-----------------------------------------------------------------------
!
!     set seafloor values to zero
!
!-----------------------------------------------------------------------

     ! !DIR$ NODEP
     ! !DIR$ COLLAPSE
      do k=1,km
      do j=1,ny_block
     ! !DIR$ NODEP
      do i=1,nx_block
         if ( k >= KMT(i,j,bid) ) then
            VISC(i,j,k  ) = c0
            VDC (i,j,k,1) = c0
            VDC (i,j,k,2) = c0
         endif
      end do
      end do

      ! k index shifted because KVMIX and KVMIX_M are at cell bottom
      ! while output axis is at cell top
      call accumulate_tavg_field(KVMIX,tavg_KVMIX,bid,k)
      call accumulate_tavg_field(KVMIX_M,tavg_KVMIX_M,bid,k)
     
      if (lniw_mixing) then
      !*** accumulated in iw_reset
      else
      ! k index shifted because bckgrnd_vdc and bckgrnd_vvc are at cell bottom
      ! while output axis is at cell top
        call accumulate_tavg_field(bckgrnd_vdc(:,:,k,bid),tavg_VDC_BCK,bid,k)
        call accumulate_tavg_field(bckgrnd_vvc(:,:,k,bid),tavg_VVC_BCK,bid,k)
      endif
       
      if (accumulate_tavg_now(tavg_TPOWER)) then
         
         WORK1(:,:) = KVMIX(:,:)*RHOMIX(:,:,k)*DBLOC(:,:,k)/ &
            (zgrid(k) - zgrid(k+1))
          
         call accumulate_tavg_field(WORK1,tavg_TPOWER,bid,k)
      endif

!-----------------------------------------------------------------------
!
!     move to next level.
!
!-----------------------------------------------------------------------

   end do
!-----------------------------------------------------------------------
!
!  fill extra coefficients for blmix
!
!-----------------------------------------------------------------------

   if (lniw_mixing) then
   !*** see iw_reset
   else
     VISC(:,:,0  ) = c0
     VDC (:,:,0,:) = c0
     VISC(:,:,km+1  ) = c0
     VDC (:,:,km+1,:) = c0
   endif

!-----------------------------------------------------------------------
!EOC
 
 end subroutine ri_iwmix












!!!------- need to be athreaded, zengyh.-------- 

!***********************************************************************
!BOP
! !IROUTINE: bldepth
! !INTERFACE:

 subroutine bldepth (CVmix_vars, DBLOC, DBSFC, TRCR, UUU, VVV, UCUR, VCUR,     &
                     STF, SHF_QSW, HBLT, USTAR, BFSFC, STABLE, KBL, this_block,&
                     SMF, SMFT)

! !DESCRIPTION:
!  This routine computes the ocean boundary layer depth defined as
!  the shallowest depth where the bulk Richardson number is equal to
!  a critical value, Ricr.
!
!  NOTE: bulk richardson numbers are evaluated by computing 
!        differences between values at zgrid(kl) $< 0$ and surface
!        reference values. currently, the reference values are equal 
!        to the values in the surface layer.  when using higher 
!        vertical grid resolution, these reference values should be 
!        computed as the vertical averages from the surface down to 
!        epssfc*zgrid(kl).
!
!  This routine also computes where surface forcing is stable 
!  or unstable (STABLE)
!
! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:

   type(cvmix_data_type), dimension(:), intent(inout) :: CVmix_vars

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TRCR                ! tracers at current time

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UUU,VVV,       &! velocities at mix     time
      UCUR,VCUR,     &! velocities at current time    Markus (lniw_mixing)
      DBLOC,         &! buoyancy difference between adjacent levels
      DBSFC           ! buoyancy difference between level and surface


   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      STF                 ! surface forcing for all tracers

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SHF_QSW             ! short-wave forcing

   real (r8), dimension(nx_block,ny_block,2), intent(in), optional :: &
      SMF,               &! surface momentum forcing at U points
      SMFT                ! surface momentum forcing at T points
                         ! *** either one or the other (not
                         ! *** both) should be passed

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(nx_block,ny_block), intent(out) :: &
      KBL                    ! index of first lvl below hbl

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      HBLT,               &! boundary layer depth
      BFSFC,              &! Bo+radiation absorbed to d
      STABLE,             &! =1 stable forcing; =0 unstab
      USTAR                ! surface friction velocity

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) :: error_string

   integer (int_kind) :: &
      i,j,               &! loop indices
      bid,               &! local block index
      kupper, kup, kdn, ktmp, kl  ! vertical level indices

   real (r8), dimension(nx_block,ny_block) :: &
      VSHEAR,            &! (velocity shear re sfc)^2
      SIGMA,             &! d/hbl
      WM, WS,            &! turb vel scale functions
      BO,                &! surface buoyancy forcing
      BOSOL,             &! radiative buoyancy forcing
      TALPHA,            &! temperature expansion coeff
      SBETA,             &! salinity    expansion coeff
      RHO1,              &! density at the surface
      WORK,              &! temp array
      ZKL,               &! depth at current z level
      B_FRQNCY,          &! buoyancy frequency
      RSH_HBLT,          &! resolved shear contribution to HBLT (fraction)
      HLANGM,            &! Langmuir depth
      HEKMAN,            &! Eckman depth limit
      HLIMIT              ! limit to mixed-layer depth
                          ! (= min(HEKMAN,HMONOB))

   real (r8), dimension(nx_block,ny_block) :: &
      niuel,             &! unresolved NIW part
      nivel               ! unresolved NIW part

   real (r8), dimension(nx_block,ny_block,3) :: &
      RI_BULK,           &! Bulk Ri number at 3 lvls
      HMONOB              ! Monin-Obukhov depth limit

   real (r8) ::          &
      absorb_frac,       &! shortwave absorption frac
      perio,             &! inertial period
      sqrt_arg,          &! dummy sqrt argument
      z_upper, z_up       ! upper depths for RI_BULK interpolation

   real (r8) :: &
      a_co, b_co, c_co    ! coefficients of the quadratic equation
                          ! $(a_{co}z^2+b_{co}|z|+c_{co}=Ri_b) used to 
                          ! find the boundary layer depth. when
                          ! finding the roots, c_co = c_co - Ricr

   real (r8) :: &
      slope_up            ! slope of the above quadratic equation
                          ! at zup. this is used as a boundary
                          ! condition to determine the coefficients.

   real (r8) :: &
      factor              ! temporary scalar factor

   real (r8) :: &
      ni_obs_factor = 0.8_r8 ! scaling factor for obs vs model

!-----------------------------------------------------------------------
!
!  compute friction velocity USTAR.  compute on U-grid and average
!  to T-grid.
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   if (present(SMFT)) then
      USTAR = sqrt(sqrt(SMFT(:,:,1)**2 + SMFT(:,:,2)**2))
   else
      WORK = sqrt(sqrt(SMF(:,:,1)**2 + SMF(:,:,2)**2))
      call ugrid_to_tgrid(USTAR,WORK,bid)
   endif

!-----------------------------------------------------------------------
!
!  compute density and expansion coefficients at surface
!
!-----------------------------------------------------------------------

   WORK = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)

   call state(1,1,WORK,TRCR(:,:,1,2),this_block, &
                  RHOFULL=RHO1, DRHODT=TALPHA, DRHODS=SBETA)

!-----------------------------------------------------------------------
!
!  compute turbulent and radiative sfc buoyancy forcing
!
!-----------------------------------------------------------------------

   do j=1,ny_block
   do i=1,nx_block
      if (RHO1(i,j) /= c0) then
         BO   (i,j) = grav*(-TALPHA(i,j)*STF(i,j,1) - &
                             SBETA (i,j)*STF(i,j,2))/RHO1(i,j)

         BOSOL(i,j) = -grav*TALPHA(i,j)*SHF_QSW(i,j)/RHO1(i,j)
      else
         BO   (i,j) = c0
         BOSOL(i,j) = c0
      endif
   end do
   end do

!-----------------------------------------------------------------------
!
!  Find bulk Richardson number at every grid level until > Ricr
!  max values when Ricr never satisfied are KBL = KMT and
!  HBLT = -zgrid(KMT)
!
!  NOTE: the reference depth is -epssfc/2.*zgrid(i,k), but the 
!        reference u,v,t,s values are simply the surface layer 
!        values and not the averaged values from 0 to 2*ref.depth,
!        which is necessary for very fine grids(top layer < 2m 
!        thickness)
!
!
!  Initialize hbl and kbl to bottomed out values
!  Initialize HEKMAN and HLIMIT (= HMONOB until reset) to model bottom
!  Initialize Monin Obukhov depth to value at z_up
!  Set HMONOB=-zgrid(km) if unstable
!
!-----------------------------------------------------------------------

   kupper = 1
   kup = 2
   kdn = 3
   z_upper = c0
   z_up    = zgrid(1)

   RI_BULK(:,:,kupper) = c0
   RI_BULK(:,:,kup) = c0 
   KBL = merge(KMT(:,:,bid), 1, (KMT(:,:,bid) > 1))

   HLANGM = c0

   do kl=1,km
      if (partial_bottom_cells) then
      	 if (kl > 1) then
      	 	  ZKL = -zgrid(kl-1) + p5*(DZT(:,:,kl  ,bid) + &
                                     DZT(:,:,kl-1,bid))
         else
            ZKL = -zgrid(1)
         endif
      else
         ZKL = -zgrid(kl)
      endif

      !DIR$ COLLAPSE
      do j=1,ny_block
      do i=1,nx_block
         if (kl == KBL(i,j)) HBLT(i,j) = ZKL(i,j)
      end do
      end do
   enddo

   if ( lcheckekmo ) then

      HEKMAN = -zgrid(km) + eps
      HLIMIT = -zgrid(km) + eps

      if ( lshort_wave ) then
         select case (sw_absorption_type)

         case ('top-layer')

            BFSFC = BO + BOSOL
         
         case ('jerlov')

            call sw_absorb_frac(-z_up,absorb_frac)
            BFSFC = BO + BOSOL * (c1 - absorb_frac)

         case ('chlorophyll')

            call sw_trans_chl(1,this_block)
            BFSFC = BO + BOSOL*(c1-TRANS(:,:,bid))

         end select

      else
         BFSFC = BO
      endif

      STABLE = merge(c1, c0, BFSFC >= c0)

      BFSFC  = BFSFC + STABLE*eps

      WORK =   STABLE * cmonob*USTAR*USTAR*USTAR/vonkar/BFSFC &
            + (STABLE -c1)*zgrid(km)
      HMONOB(:,:,kup) = merge( -z_up+eps, WORK, WORK <= -z_up )
   endif

   RSH_HBLT = c0

!-----------------------------------------------------------------------
!
!  compute velocity shear squared on U-grid and use the maximum
!  of the four surrounding U-grid values for the T-grid.
!
!-----------------------------------------------------------------------

! account for the unresolved part of the NI velocity (equal to the resolved); Markus 9/27/11
! the 0.05 accounts like in En for the average over the abs(cycle), loosely based on
! 1/T * integral(cos), also tuned to get En = u'Tau' from model
! Model resolves half NI energy, so the unresolved part, nivel is added again 
! x4 because for 16 hours 1/f a time step of 1 hours gets only a quarter of max. increase

   if ( lniw_mixing .and. linertial ) then

     do j=1,ny_block
       do i=1,nx_block
         niuel(i,j)=c0
         nivel(i,j)=c0

         if( TLATD(i,j,bid) > 5.0_r8 ) then
 
          factor = ni_obs_factor*abs(12.0_r8*3600.0_r8/sin(TLAT(i,j,bid)))/pi2/dtt
          niuel(i,j) = - factor * (VCUR(i,j,1) - VVV(i,j,1))     
          nivel(i,j) =   factor * (UCUR(i,j,1) - UUU(i,j,1))     
    
           if( TLATD(i,j,bid) < 10.0_r8 ) then
            niuel(i,j) = niuel(i,j) * NIW_COS_FACTOR(i,j,bid)
            nivel(i,j) = nivel(i,j) * NIW_COS_FACTOR(i,j,bid)
           endif

         endif


         if( TLATD(i,j,bid) < -5.0_r8 ) then

          factor = ni_obs_factor*abs(12.0_r8*3600.0_r8/sin(TLAT(i,j,bid)))/pi2/dtt
          niuel(i,j) =   factor * (VCUR(i,j,1) - VVV(i,j,1))     
          nivel(i,j) = - factor * (UCUR(i,j,1) - UUU(i,j,1))     

           if( TLATD(i,j,bid) > -10.0_r8 ) then
            niuel(i,j) = niuel(i,j) * NIW_COS_FACTOR(i,j,bid)
            nivel(i,j) = nivel(i,j) * NIW_COS_FACTOR(i,j,bid)
           endif

         endif

       enddo
     enddo

   endif  ! if ( lniw_mixing .and. linertial ) then

   do kl = 2,km
 
      if ( lniw_mixing ) then
        WORK = (UUU(:,:,1) + niuel(:,:) - UUU(:,:,kl))**2 + &
               (VVV(:,:,1) + nivel(:,:) - VVV(:,:,kl))**2 
      else
        WORK = (UUU(:,:,1) - UUU(:,:,kl))**2 + &
               (VVV(:,:,1) - VVV(:,:,kl))**2
      endif
     
      if (partial_bottom_cells) then
         WORK = WORK/(-zgrid(kl-1) + & 
                      p5*(DZU(:,:,kl  ,bid) + &
                          DZU(:,:,kl-1,bid) - &
                          DZU(:,:,1   ,bid)))**2

      	 ZKL = -zgrid(kl-1) + p5*(DZT(:,:,kl  ,bid) + &
                                  DZT(:,:,kl-1,bid))
      else
         ZKL = -zgrid(kl)
      endif

      VSHEAR(:,1) = c0
      do j=2,ny_block
         VSHEAR(1,j) = c0
         do i=2,nx_block
            VSHEAR(i,j) = max(WORK(i,j  ), WORK(i-1,j  ),   &
                              WORK(i,j-1), WORK(i-1,j-1))
         enddo
      enddo

!-----------------------------------------------------------------------
!
!     compute bfsfc= Bo + radiative contribution down to hbf * hbl
!     add epsilon to BFSFC to ensure never = 0
!
!-----------------------------------------------------------------------

      if (lshort_wave) then

         select case (sw_absorption_type)

         case ('top-layer')

           BFSFC = BO + BOSOL

         case ('jerlov')

           do j=1,ny_block
           do i=1,nx_block
              call sw_absorb_frac(ZKL(i,j), absorb_frac)
              BFSFC(i,j) = BO(i,j) + BOSOL(i,j)*(c1 - absorb_frac) 
           enddo
           enddo

         case ('chlorophyll')

           call sw_trans_chl(2*kl-1,this_block)
           BFSFC = BO + BOSOL*(c1-TRANS(:,:,bid))

         end select


      else
         BFSFC = BO
      endif

      STABLE = MERGE(c1, c0, BFSFC >= c0)

      BFSFC  = BFSFC + STABLE*eps

!-----------------------------------------------------------------------
!
!     compute the Ekman and Monin Obukhov depths using above stability
!
!-----------------------------------------------------------------------

      if (lcheckekmo) then

         !DIR$ COLLAPSE
         do j=1,ny_block
         do i=1,nx_block
            if ( STABLE(i,j) > p5 .and. HEKMAN(i,j) >= -zgrid(km) ) then
               HEKMAN(i,j) = max(ZKL(i,j), &
                           cekman*USTAR(i,j)/(abs(FCORT(i,j,bid))+eps))
            endif
         end do
         end do

         HMONOB(:,:,kdn) = STABLE*cmonob*USTAR*USTAR*USTAR/vonkar/BFSFC + &
                          (STABLE-c1)*zgrid(km)

         !DIR$ COLLAPSE
         do j=1,ny_block
         do i=1,nx_block
            if (HMONOB(i,j,kdn) <= ZKL(i,j) .and. &
                HMONOB(i,j,kup) >  -z_up) then
               WORK(i,j) = (HMONOB(i,j,kdn) - HMONOB(i,j,kup))/ &
                           (z_up + ZKL(i,j))
               HLIMIT(i,j) = (HMONOB(i,j,kdn) - WORK(i,j)*ZKL(i,j))/ &
                             (c1 - WORK(i,j))

            endif
         end do
         end do
      endif

!-----------------------------------------------------------------------
!
!     compute velocity scales at sigma, for hbl = -zgrid(kl)
!
!-----------------------------------------------------------------------

      SIGMA = epssfc

      call wscale(SIGMA, ZKL, USTAR, BFSFC, 2, WM, WS)

!-----------------------------------------------------------------------
!
!     compute the turbulent shear contribution to RI_BULK and store
!     in WM.
!
!-----------------------------------------------------------------------

      if (partial_bottom_cells) then
         if (kl < km) then
            B_FRQNCY = sqrt( &
                       p5*(DBLOC(:,:,kl) + abs(DBLOC(:,:,kl)) + eps2)/  &
                      (p5*(DZT(:,:,kl,bid) + DZT(:,:,kl+1,bid))) )
         else
            B_FRQNCY = sqrt( &
                       p5*(DBLOC(:,:,kl) + abs(DBLOC(:,:,kl)) + eps2)/  &
                       DZT(:,:,kl,bid) )
         end if
      else
         B_FRQNCY = sqrt( &
                    p5*(DBLOC(:,:,kl) + abs(DBLOC(:,:,kl)) + eps2)/  &
                    (zgrid(kl)-zgrid(kl+1)) )
      endif

      WM = ZKL*WS*B_FRQNCY* &
          ( (Vtc/Ricr(kl))*max(2.1_r8 - 200.0_r8*B_FRQNCY,concv) )

!-----------------------------------------------------------------------
! 
!     compute bulk Richardson number at new level
!
!-----------------------------------------------------------------------

      if (partial_bottom_cells) then
         WORK = merge( DBSFC(:,:,kl)/(-zgrid(kl-1)+            &
                                      p5*(DZT(:,:,kl-1,bid) +  &
                                          DZT(:,:,kl  ,bid) -  &
                                          DZT(:,:,1   ,bid))), & 
                       c0, KMT(:,:,bid) >= kl)
         WM = WM/(-zgrid(kl-1) +          &
                  p5*(DZT(:,:,kl-1,bid) + &
                      DZT(:,:,kl  ,bid) - &
                      DZT(:,:,1   ,bid)))**2
         RI_BULK(:,:,kdn) = WORK/(VSHEAR+WM+eps/(-zgrid(kl-1)+   &
                                   p5*(DZU(:,:,kl,bid) +         &
                                       DZU(:,:,kl-1,bid) -       &
                                       DZU(:,:,1,bid)))**2)
      else
         WORK = MERGE( (zgrid(1)-zgrid(kl))*DBSFC(:,:,kl), &
                      c0, KMT(:,:,bid) >= kl)

         if (lniw_mixing) then
           RI_BULK(:,:,kdn) = WORK/(VSHEAR+WM+eps)
         else
           if ( linertial ) then
             RI_BULK(:,:,kdn) = WORK/(VSHEAR+WM+USTAR*BOLUS_SP(:,:,bid)+eps)
           else
             RI_BULK(:,:,kdn) = WORK/(VSHEAR+WM+eps)
           endif
         endif
      endif

!-----------------------------------------------------------------------
!
!       find hbl where Rib = Ricr. if possible, use a quadratic
!       interpolation. if not, linearly interpolate. the quadratic
!       equation coefficients are determined using the slope and
!       Ri_bulk at z_up and Ri_bulk at zgrid(kl). the slope at
!       z_up is computed linearly between z_upper and z_up.
!
!       compute Langmuir depth always 
!-----------------------------------------------------------------------

      do j=1,ny_block
      do i=1,nx_block
         if ( KBL(i,j) == KMT(i,j,bid) .and.  &
              RI_BULK(i,j,kdn) > Ricr(kl) ) then

            slope_up =  (RI_BULK(i,j,kupper) - RI_BULK(i,j,kup))/ &
                        (z_up - z_upper)
            a_co = (RI_BULK(i,j,kdn) - RI_BULK(i,j,kup) -         &
                    slope_up*(ZKL(i,j) + z_up) )/(z_up + ZKL(i,j))**2
            b_co = slope_up + c2 * a_co * z_up
            c_co = RI_BULK(i,j,kup) + &
                   z_up*(a_co*z_up + slope_up) - Ricr(kl)
            sqrt_arg = b_co**2 - c4*a_co*c_co

            if ( ( abs(b_co) > eps .and. abs(a_co)/abs(b_co) <= eps ) &
                 .or. sqrt_arg <= c0 ) then
                 	
               HBLT(i,j) = -z_up + (z_up + ZKL(i,j)) *               &
                           (Ricr(kl)         - RI_BULK(i,j,kup))/    &
                           (RI_BULK(i,j,kdn) - RI_BULK(i,j,kup))
            else
               HBLT(i,j) = (-b_co + sqrt(sqrt_arg)) / (c2*a_co)
            endif
            
            KBL(i,j) = kl
            RSH_HBLT(i,j) =  (VSHEAR(i,j)*Ricr(kl)/ &
                              (DBSFC(i,j,kl)+eps))/HBLT(i,j)

            HLANGM(i,j) = USTAR(i,j) * SQRT( FSTOKES(i,j,bid)*ZKL(i,j)/(DBSFC(i,j,kl)+eps) ) &
                          / 0.9_r8

         endif
      enddo
      enddo

!-----------------------------------------------------------------------
!
!     swap klevel indices and move to next level
!
!-----------------------------------------------------------------------

      ktmp   = kupper
      kupper = kup
      kup    = kdn
      kdn    = ktmp
      z_upper = z_up
      z_up    = zgrid(kl)

   end do

!-----------------------------------------------------------------------
!
!     apply Langmuir parameterization if requested 
!
!-----------------------------------------------------------------------

   if ( llangmuir ) then
     do kl = km,2,-1
        where ( HLANGM > HBLT          .and.   &
                HLANGM >  -zgrid(kl-1) .and.   &
                HLANGM <= ZKL                  )
           HBLT  = HLANGM
           KBL   = kl
        end where
     enddo
   endif

!-----------------------------------------------------------------------
!
!  first combine Ekman and Monin-Obukhov depth limits. then apply
!  these restrictions to HBLT. note that HLIMIT is set to -zgrid(km)
!  in unstable forcing.
!
!-----------------------------------------------------------------------

   if ( lcheckekmo ) then

      where ( HEKMAN < HLIMIT )  HLIMIT = HEKMAN

      do kl = 2,km
         where ( HLIMIT < HBLT          .and.   &
                 HLIMIT >  -zgrid(kl-1) .and.   &
                 HLIMIT <= ZKL                  )
            HBLT = HLIMIT
            KBL = kl
         end where
      enddo

   endif

!-----------------------------------------------------------------------
!
!  apply a Gaussian filter
!
!-----------------------------------------------------------------------

   if (.not.l1Ddyn) then
     call smooth_hblt (.true., .false., bid, HBLT=HBLT, KBL=KBL)
   end if

!-----------------------------------------------------------------------
!
!  correct stability and buoyancy forcing for SW up to boundary layer
!
!-----------------------------------------------------------------------

   if (lshort_wave) then
      select case (sw_absorption_type)

      case ('top-layer')

        BFSFC   = BO + BOSOL
!       QSW_HBL = SHF_QSW
        if (accumulate_tavg_now(tavg_QSW_HBL)) then
           WORK = SHF_QSW/hflux_factor
           call accumulate_tavg_field(WORK,tavg_QSW_HBL,bid,1)
        endif


      case ('jerlov')

         do j = 1,ny_block
         do i = 1,nx_block
            call sw_absorb_frac(HBLT(i,j),absorb_frac)
            BFSFC(i,j)  = BO(i,j) + BOSOL(i,j)*(c1 - absorb_frac) 
         enddo
         enddo

         if (accumulate_tavg_now(tavg_QSW_HBL)) then
           !QSW_HBL(i,j) = SHF_QSW(i,j)*(c1-absorb_frac)  ! boundary layer sw
            WORK = SHF_QSW*(c1-absorb_frac)/hflux_factor
            call accumulate_tavg_field(WORK,tavg_QSW_HBL,bid,1)
         endif

      case ('chlorophyll')

         ZTRANS(:,:,bid) = HBLT(:,:)
         call sw_trans_chl(0,this_block)
         BFSFC   = BO + BOSOL*(c1-TRANS(:,:,bid))

         if (accumulate_tavg_now(tavg_QSW_HBL)) then
           !QSW_HBL = SHF_QSW   *(c1-TRANS(:,:,bid)) ! boundary layer sw heating
            WORK = SHF_QSW*(c1-TRANS(:,:,bid))/hflux_factor
            call accumulate_tavg_field(WORK,tavg_QSW_HBL,bid,1)
         endif

      end select


   endif

   STABLE = MERGE(c1, c0, BFSFC >= c0)


   BFSFC  = BFSFC + STABLE * eps ! ensures bfsfc never=0

!-----------------------------------------------------------------------
!EOC

 end subroutine bldepth


!***********************************************************************
!BOP
! Athreaded subroutine added by zengyh 2018-11-01:
! !IROUTINE: bldepth_sw
! !INTERFACE:

 subroutine bldepth_sw (CVmix_vars, DBLOC, DBSFC, TRCR, UUU, VVV, UCUR, VCUR,     &
                     STF, SHF_QSW, HBLT, USTAR, BFSFC, STABLE, KBL, this_block,&
                     SMF, SMFT)

! !DESCRIPTION:
!  This routine computes the ocean boundary layer depth defined as
!  the shallowest depth where the bulk Richardson number is equal to
!  a critical value, Ricr.
!
!  NOTE: bulk richardson numbers are evaluated by computing 
!        differences between values at zgrid(kl) $< 0$ and surface
!        reference values. currently, the reference values are equal 
!        to the values in the surface layer.  when using higher 
!        vertical grid resolution, these reference values should be 
!        computed as the vertical averages from the surface down to 
!        epssfc*zgrid(kl).
!
!  This routine also computes where surface forcing is stable 
!  or unstable (STABLE)
!
! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:

   type(cvmix_data_type), dimension(:), intent(inout) :: CVmix_vars

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TRCR                ! tracers at current time

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UUU,VVV,       &! velocities at mix     time
      UCUR,VCUR,     &! velocities at current time    Markus (lniw_mixing)
      DBLOC,         &! buoyancy difference between adjacent levels
      DBSFC           ! buoyancy difference between level and surface


   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      STF                 ! surface forcing for all tracers

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SHF_QSW             ! short-wave forcing

   real (r8), dimension(nx_block,ny_block,2), intent(in), optional :: &
      SMF,               &! surface momentum forcing at U points
      SMFT                ! surface momentum forcing at T points
                         ! *** either one or the other (not
                         ! *** both) should be passed

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(nx_block,ny_block), intent(out) :: &
      KBL                    ! index of first lvl below hbl

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      HBLT,               &! boundary layer depth
      BFSFC,              &! Bo+radiation absorbed to d
      STABLE,             &! =1 stable forcing; =0 unstab
      USTAR                ! surface friction velocity

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) :: error_string

   integer (int_kind) :: &
      i,j,               &! loop indices
      bid,               &! local block index
      kupper, kup, kdn, ktmp, kl  ! vertical level indices

   real (r8), dimension(nx_block,ny_block) :: &
      VSHEAR,            &! (velocity shear re sfc)^2
      SIGMA,             &! d/hbl
      WM, WS,            &! turb vel scale functions
      BO,                &! surface buoyancy forcing
      BOSOL,             &! radiative buoyancy forcing
      TALPHA,            &! temperature expansion coeff
      SBETA,             &! salinity    expansion coeff
      RHO1,              &! density at the surface
      WORK,              &! temp array
      ZKL,               &! depth at current z level
      B_FRQNCY,          &! buoyancy frequency
      RSH_HBLT,          &! resolved shear contribution to HBLT (fraction)
      HLANGM,            &! Langmuir depth
      HEKMAN,            &! Eckman depth limit
! variables added by zengyh 2018-11-01:
      ZETA,              &! d/L or sigma*hbl/L(monin-obk)
      ZETAH,             &! sigma*hbl*vonkar*BFSFC or ZETA = ZETAH/USTAR**3
! variables added by zengyh 2018-11-01.
      HLIMIT              ! limit to mixed-layer depth
                          ! (= min(HEKMAN,HMONOB))

   real (r8), dimension(nx_block,ny_block) :: &
      niuel,             &! unresolved NIW part
      nivel               ! unresolved NIW part

   real (r8), dimension(nx_block,ny_block,3) :: &
      RI_BULK,           &! Bulk Ri number at 3 lvls
      HMONOB              ! Monin-Obukhov depth limit

   real (r8) ::          &
      absorb_frac,       &! shortwave absorption frac
      perio,             &! inertial period
      sqrt_arg,          &! dummy sqrt argument
      z_upper, z_up       ! upper depths for RI_BULK interpolation

   real (r8) :: &
      a_co, b_co, c_co    ! coefficients of the quadratic equation
                          ! $(a_{co}z^2+b_{co}|z|+c_{co}=Ri_b) used to 
                          ! find the boundary layer depth. when
                          ! finding the roots, c_co = c_co - Ricr

   real (r8) :: &
      slope_up            ! slope of the above quadratic equation
                          ! at zup. this is used as a boundary
                          ! condition to determine the coefficients.

   real (r8) :: &
      factor              ! temporary scalar factor

   real (r8) :: &
      ni_obs_factor = 0.8_r8 ! scaling factor for obs vs model


! variables added by zengyh 2018-11-3:
   integer (int_kind) :: &
      locbid, loc_short_wave, loc_partial_bottom_cells, loc_lniw_mixing, loc_linertial

#if defined(POPINSJN)
   real(8)::absVDCerr=0.0
   real(8)::tmpVDC
   real(8),allocatable, dimension(:,:,:) :: &
      MWORK21,SWORK21,             &
      MWORK11,SWORK11,             &
      MWORK1KM,SWORK1KM,MWORK2KM,SWORK2KM,         &
      MWORK1K,MWORK2K,SWORK1K,SWORK2K
   real(8)::avevdcerr=0.0
#endif




! added by zengyh 2018-11-03:
   interface
     subroutine m_kpp_bldepth_loop2(TRCR, DBLOC, DBSFC, &
          UUU, VVV, UCUR, VCUR,                         &
          HBLT, USTAR, BFSFC, STABLE,                   &
          KMT, KBL,                                     &
          DZT, DZU, niuel, nivel,                       &
          Tr, CHLINDX, BO, BOSOL, ZKL,                  &
          myproc, nx_block, ny_block, km, nt, bid,           &
          lshort_wave, partial_bottom_cells, lniw_mixing, linertial, &
          zgrid, Ricr,                                  & 
          Vtc, zeta_s, a_s, c_s, vonkar,                &
          eps, eps2, epssfc, concv                      &
#if defined(POPCHKJN) || defined(POPTIMEJN)
          ,steps                                        &
#endif
#if defined(POPINSJN)
          ,SWORK11,SWORK21,SWORK1K,SWORK2K,SWORK1KM,SWORK2KM &
#endif
          ) bind(C, name='m_kpp_bldepth_loop2') 
         real(8):: TRCR(*), DBLOC(*), DBSFC(*), UUU(*), VVV(*), UCUR(*), VCUR(*),  &
               HBLT(*), USTAR(*), BFSFC(*), STABLE(*)
         real(8):: DZT(*), DZU(*), niuel(*), nivel(*), Tr(*), &
               BO(*), BOSOL(*), ZKL(*)
         integer,value:: myproc, nx_block, ny_block, km, nt, bid 
         integer,value:: lshort_wave, partial_bottom_cells, lniw_mixing, linertial  
         integer:: KMT(*), KBL(*), CHLINDX(*)
         real(8):: zgrid(*), Ricr(*)
         real(8),value:: Vtc, zeta_s, a_s, c_s, vonkar, eps, eps2, epssfc, concv
#if defined(POPCHKJN) || defined(POPTIMEJN)
         integer,value:: steps
#endif
#if defined(POPINSJN)
         real(8),intent(out):: SWORK11(*), SWORK1K(*), SWORK1KM(*) 
         real(8),intent(out):: SWORK21(*), SWORK2K(*), SWORK2KM(*) 
#endif
     end subroutine m_kpp_bldepth_loop2
   end interface



  
!-----------------------------------------------------------------------
!
!  compute friction velocity USTAR.  compute on U-grid and average
!  to T-grid.
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   if (present(SMFT)) then
      USTAR = sqrt(sqrt(SMFT(:,:,1)**2 + SMFT(:,:,2)**2))
   else
      WORK = sqrt(sqrt(SMF(:,:,1)**2 + SMF(:,:,2)**2))
! kpp_call(63)      
      call ugrid_to_tgrid(USTAR,WORK,bid)
   endif

!-----------------------------------------------------------------------
!
!  compute density and expansion coefficients at surface
!
!-----------------------------------------------------------------------

   WORK = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)

! kpp_call(64): 636(us), called times:  1(s);  
   call state(1,1,WORK,TRCR(:,:,1,2),this_block, &
                  RHOFULL=RHO1, DRHODT=TALPHA, DRHODS=SBETA)

!-----------------------------------------------------------------------
!
!  compute turbulent and radiative sfc buoyancy forcing
!
!-----------------------------------------------------------------------

! kpp_loop(18):  204(us), called times:  1(s);
   do j=1,ny_block
   do i=1,nx_block
      if (RHO1(i,j) /= c0) then
         BO   (i,j) = grav*(-TALPHA(i,j)*STF(i,j,1) - &
                             SBETA (i,j)*STF(i,j,2))/RHO1(i,j)

         BOSOL(i,j) = -grav*TALPHA(i,j)*SHF_QSW(i,j)/RHO1(i,j)
      else
         BO   (i,j) = c0
         BOSOL(i,j) = c0
      endif
   end do
   end do

!-----------------------------------------------------------------------
!
!  Find bulk Richardson number at every grid level until > Ricr
!  max values when Ricr never satisfied are KBL = KMT and
!  HBLT = -zgrid(KMT)
!
!  NOTE: the reference depth is -epssfc/2.*zgrid(i,k), but the 
!        reference u,v,t,s values are simply the surface layer 
!        values and not the averaged values from 0 to 2*ref.depth,
!        which is necessary for very fine grids(top layer < 2m 
!        thickness)
!
!
!  Initialize hbl and kbl to bottomed out values
!  Initialize HEKMAN and HLIMIT (= HMONOB until reset) to model bottom
!  Initialize Monin Obukhov depth to value at z_up
!  Set HMONOB=-zgrid(km) if unstable
!
!-----------------------------------------------------------------------



 
! start to athread codes by zengyh 2018-10-31: 

if (lshort_wave) then
  loc_short_wave=1
else
  loc_short_wave=0
endif
if (partial_bottom_cells) then 
  loc_partial_bottom_cells=1
else
  loc_partial_bottom_cells=0
endif
if (lniw_mixing) then
  loc_lniw_mixing=1
else
  loc_lniw_mixing=0
endif
if (linertial) then
  loc_linertial=1
else
  loc_linertial=0
endif

locbid = bid-1
#if defined(POPTESTJN)
if((my_task==0 .or. my_task==256 .or. my_task==512 .or. my_task==768) .and. num_steps<=512) then
   write(*,*) 'PROC: ', my_task,'; bid: ', bid
   write(*,*) 'ny_block: ', ny_block, '; nx_block: ', nx_block, '; km: ', km, '; nt: ', nt
endif
#endif

#if defined(POPTESTJN)
if(my_task==0 .and. num_steps==0) then
   write(*,*) 'zgrid: '
   write(*,*) zgrid(0:km+1)
endif

if(my_task==0 .and. num_steps<=3) then
   write(*,*) 'UUU(:,8,1):'
   write(*,*) UUU(:,8,1) 
   write(*,*) 'UUU(:,9,1):'
   write(*,*) UUU(:,9,1)   
endif
#endif



call m_kpp_bldepth_loop2(TRCR, DBLOC, DBSFC,  &
      UUU, VVV, UCUR(:,:,1), VCUR(:,:,1),  &
      HBLT, USTAR, BFSFC, STABLE,&
      KMT(:,:,bid), KBL, &
      DZT(:,:,:,bid), DZU(:,:,:,bid), niuel, nivel, &
      Tr, CHLINDX(:,:,bid), BO, BOSOL, ZKL,&
!      KMT(:,:,1), KBL, &
!      DZT(:,:,:,1), DZU(:,:,:,1), niuel, nivel, &
!      Tr, CHLINDX(:,:,1), BO, BOSOL, ZKL,&
      my_task, nx_block, ny_block, km, nt, locbid, &
      loc_short_wave, loc_partial_bottom_cells, loc_lniw_mixing, loc_linertial, &
      zgrid, Ricr, & 
      Vtc, zeta_s, a_s, c_s, vonkar, &
      eps, eps2, epssfc, concv &
#if defined(POPCHKJN) || defined(POPTIMEJN)
      , num_steps &
#endif
#if defined(POPINSJN)
      , SWORK11,SWORK21,SWORK1K,SWORK2K,SWORK1KM,SWORK2KM &
#endif
      )

! end of athread codes by zengyh 2018-10-31.



#if defined(POPINSJN)
   if (my_task==0 .and. num_steps<=16) then
      write(*,*) 'lshort_wave: ', lshort_wave   
      write(*,*) 'partial_bottom_cells: ', partial_bottom_cells   
      write(*,*) 'lniw_mixing: ', lniw_mixing   
      write(*,*) 'linertial: ', linertial   
      write(*,*) 'lcheckekmo: ', lcheckekmo   
      write(*,*) 'llangmuir: ', llangmuir   
      write(*,*) 'l1Ddyn: ', l1Ddyn
   endif
#endif

!! llangmuir : False.
!! omit if ( llangmuir ) here.
!
!!-----------------------------------------------------------------------
!!
!!     apply Langmuir parameterization if requested 
!!
!!-----------------------------------------------------------------------



!! lcheckekmo : False.
!! omit if ( lcheckekmo ) here.

!!-----------------------------------------------------------------------
!!
!!  first combine Ekman and Monin-Obukhov depth limits. then apply
!!  these restrictions to HBLT. note that HLIMIT is set to -zgrid(km)
!!  in unstable forcing.
!!
!!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!
!  apply a Gaussian filter
!
!-----------------------------------------------------------------------

! l1Ddyn :  False.
   if (.not.l1Ddyn) then
   
! kpp_call(70):  6135(us), called times:  1(s);   
     call smooth_hblt (.true., .false., bid, HBLT=HBLT, KBL=KBL)
   end if

!-----------------------------------------------------------------------
!
!  correct stability and buoyancy forcing for SW up to boundary layer
!
!-----------------------------------------------------------------------

   if (lshort_wave) then
      select case (sw_absorption_type)
 
      case ('chlorophyll')

         ZTRANS(:,:,bid) = HBLT(:,:)
 
! kpp_call(74):  2359(us), called times:  1(s);		 
         call sw_trans_chl(0,this_block)
         BFSFC   = BO + BOSOL*(c1-TRANS(:,:,bid))

! noted out by zengyh 2018-11-01:
!         if (accumulate_tavg_now(tavg_QSW_HBL)) then
!           !QSW_HBL = SHF_QSW   *(c1-TRANS(:,:,bid)) ! boundary layer sw heating
!            WORK = SHF_QSW*(c1-TRANS(:,:,bid))/hflux_factor
!! kpp_call(75)			
!            call accumulate_tavg_field(WORK,tavg_QSW_HBL,bid,1)
!         endif

      end select

   endif

   STABLE = MERGE(c1, c0, BFSFC >= c0)

   BFSFC  = BFSFC + STABLE * eps ! ensures bfsfc never=0

   
   
! modified by zengyh 2018-11-01:
   if (lshort_wave) then
      select case (sw_absorption_type)
 
      case ('chlorophyll')
  
         if (accumulate_tavg_now(tavg_QSW_HBL)) then
           !QSW_HBL = SHF_QSW   *(c1-TRANS(:,:,bid)) ! boundary layer sw heating
            WORK = SHF_QSW*(c1-TRANS(:,:,bid))/hflux_factor
! kpp_call(75)			
            call accumulate_tavg_field(WORK,tavg_QSW_HBL,bid,1)
         endif

      end select

   endif   
  
#if defined(POPCHKJN) 
   num_steps = num_steps + 1
#endif
! modified by zengyh 2018-11-01.   

!-----------------------------------------------------------------------
!EOC

 end subroutine bldepth_sw










!***********************************************************************
!BOP
! !IROUTINE: blmix
! !INTERFACE:

 subroutine blmix(CVmix_vars, VISC, VDC, HBLT, USTAR, BFSFC, STABLE, &
                  KBL, GHAT, GHAT2, this_block) 

! !DESCRIPTION:
!  This routine computes mixing coefficients within boundary layer 
!  which depend on surface forcing and the magnitude and gradient 
!  of interior mixing below the boundary layer (matching).  These
!  quantities have been computed in other routines.
!
!  Caution: if mixing bottoms out at hbl = -zgrid(km) then
!  fictitious layer at km+1 is needed with small but finite width 
!  hwide(km+1).
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type(cvmix_data_type), dimension(:), intent(inout) :: CVmix_vars

   real (r8), dimension(nx_block,ny_block,0:km+1), intent(inout) :: & 
      VISC               ! interior mixing coeff on input
                         ! combined interior/bndy layer coeff output

   real (r8), dimension(nx_block,ny_block,0:km+1,2),intent(inout) :: &
      VDC        ! diffusivity for tracer diffusion

! !INPUT PARAMETERS:

   integer (int_kind), dimension(nx_block,ny_block), intent(in) ::  &
      KBL                    ! index of first lvl below hbl

   real (r8), dimension(nx_block,ny_block), intent(in) ::  &
      HBLT,                 & ! boundary layer depth
      BFSFC,                & ! surface buoyancy forcing
      STABLE,               & ! =1 stable forcing; =0 unstab
      USTAR                   ! surface friction velocity

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km),intent(out) :: &
      GHAT                ! non-local mixing coefficient
   real (r8), dimension(nx_block,ny_block,km,2),intent(inout) :: &
      GHAT2                ! non-local mixing coefficient (CVMix)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k,kp1,             &! dummy k level index
      i,j,               &! horizontal indices
      bid,               &! local block index
      ic, nlev

   integer (int_kind), dimension(nx_block,ny_block) :: &
      KN                  ! klvl closest to HBLT

   real (r8), dimension(nx_block,ny_block,km,3) :: &
      BLMC                ! bndy layer mixing coefs

   real (r8), dimension(nx_block,ny_block,3) :: &
      GAT1,              &! shape function at sigma=1
      DAT1,              &! derivative of shape function 
      DKM1                ! bndy layer difs at kbl-1 lvl

   real (r8), dimension(nx_block,ny_block) :: &
      WM,WS,             &! turbulent velocity scales
      CASEA,             &! =1 in case A, =0 in case B     
      SIGMA,             &! normalized depth (d/hbl)
      VISCH,             &! viscosity at hbl
      DIFTH,             &! temp diffusivity at hbl
      DIFSH,             &! tracer diffusivity at hbl
      DELHAT, R, DVDZUP, DVDZDN, &
      VISCP, DIFTP, DIFSP, F1, &
      WORK1,WORK2

!-----------------------------------------------------------------------
!
!  compute velocity scales at hbl
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   SIGMA = epssfc

   if (lcvmix) then
     GHAT = 0.0_r8
     GHAT2 = 0.0_r8
     do i=1,nx_block
       do j=1,ny_block
         k = KBL(i,j)
         CASEA(i,j)  = p5 + SIGN(p5, -zgrid(k)-p5*hwide(k)-HBLT(i,j))
         KN(i,j) = NINT(CASEA(i,j))*(k-1) + (1-NINT(CASEA(i,j)))*k
       end do
     end do

     do i=1,nx_block
       do j=1,ny_block
         ic = (j-1)*nx_block + i
         if (KMT(i,j,bid).gt.0) then
           nlev = CVmix_vars(ic)%nlev

           ! Pack into CVMix data type (convert from cgs to mks)
           CVmix_vars(ic)%Mdiff_iface(2:nlev+1) = VISC(i,j,1:nlev)*1e-4_r8
           CVmix_vars(ic)%Tdiff_iface(2:nlev+1) =  VDC(i,j,1:nlev,1)*1e-4_r8
           CVmix_vars(ic)%Sdiff_iface(2:nlev+1) =  VDC(i,j,1:nlev,2)*1e-4_r8
           CVmix_vars(ic)%BoundaryLayerDepth    = HBLT(i,j)*1e-2_r8
           CVmix_vars(ic)%kOBL_depth = real(KBL(i,j),r8)+p5*(p5-CASEA(i,j))
           CVmix_vars(ic)%SurfaceFriction = USTAR(i,j)*1e-2_r8
           CVmix_vars(ic)%SurfaceBuoyancyForcing = BFSFC(i,j)*1e-4_r8

           call cvmix_coeffs_kpp(CVmix_vars(ic))

           ! Unpack CVMix data type (convert from mks to cgs)
           VISC( i,j,1:nlev)   = CVmix_vars(ic)%Mdiff_iface(2:nlev+1)*1e4_r8
           VDC(  i,j,1:nlev,1) = CVmix_vars(ic)%Tdiff_iface(2:nlev+1)*1e4_r8
           VDC(  i,j,1:nlev,2) = CVmix_vars(ic)%Sdiff_iface(2:nlev+1)*1e4_r8
           GHAT2(i,j,1:nlev,1) = CVmix_vars(ic)%kpp_Tnonlocal_iface(2:nlev+1)
           GHAT2(i,j,1:nlev,2) = CVmix_vars(ic)%kpp_Snonlocal_iface(2:nlev+1)
         end if
       end do
     end do

   else ! No cvmix
#ifdef NOATH
     call wscale(SIGMA, HBLT, USTAR, BFSFC, 3, WM, WS)
#endif
     call loop1_blmix(HBLT,BFSFC,USTAR,WM,WS,nx_block,ny_block,my_task)

!-----------------------------------------------------------------------
!
!  determine caseA = 0 if closer to KBL than KBL-1
!  KN is then the closest klevel to HBLT
!
!-----------------------------------------------------------------------

     if (partial_bottom_cells) then
        !DIR$ COLLAPSE
        do j=1,ny_block
        do i=1,nx_block
           k = KBL(i,j)
           if (k == 1) then
              CASEA(i,j)  = p5 + SIGN(p5, -zgrid(0)-HBLT(i,j))
           else
              CASEA(i,j)  = p5 + SIGN(p5, &
                             -zgrid(k-1)+p5*DZT(i,j,k-1,bid)-HBLT(i,j))
           endif
        enddo
        enddo
     else
        !DIR$ COLLAPSE
        do j=1,ny_block
        do i=1,nx_block
           k = KBL(i,j)
           CASEA(i,j)  = p5 + SIGN(p5, -zgrid(k)-p5*hwide(k)-HBLT(i,j))
        enddo
        enddo
     endif

     KN = NINT(CASEA)*(KBL-1) + (1-NINT(CASEA))*KBL

!-----------------------------------------------------------------------
!
!  find the interior viscosities and derivatives at hbl by 
!  interpolating derivative values at vertical interfaces.  compute
!  matching conditions for shape function.  
!
!-----------------------------------------------------------------------

     F1 = STABLE*c5*BFSFC/(USTAR**4+eps)

     do k=1,km

        if (partial_bottom_cells) then

           if (k == 1) then
              WORK1 = c0
           else            
              WORK1 = DZT(:,:,k-1,bid)
           end if
           if (k == km) then
              WORK2 = eps
           else
              WORK2 = DZT(:,:,k+1,bid)
           end if

           !DIR$ COLLAPSE
           do j=1,ny_block
           do i=1,nx_block
              if (k == KN(i,j)) then

              DELHAT(i,j) = - zgrid(k-1) + DZT(i,j,k,bid) + &
                              p5*WORK1(i,j) - HBLT(i,j)  
              R     (i,j) = c1 - DELHAT(i,j) /DZT(i,j,k,bid)

              DVDZUP(i,j) = (VISC(i,j,k-1) - VISC(i,j,k  ))/DZT(i,j,k,bid)
              DVDZDN(i,j) = (VISC(i,j,k  ) - VISC(i,j,k+1))/WORK2(i,j)
              VISCP (i,j) = p5*( (c1-R(i,j))* &
                                 (DVDZUP(i,j) + abs(DVDZUP(i,j))) + &
                                     R(i,j) * &
                                 (DVDZDN(i,j) + abs(DVDZDN(i,j))) )

              DVDZUP(i,j) = (VDC(i,j,k-1,2) - VDC(i,j,k  ,2))/DZT(i,j,k,bid)
              DVDZDN(i,j) = (VDC(i,j,k  ,2) - VDC(i,j,k+1,2))/WORK2(i,j)
              DIFSP (i,j) = p5*( (c1-R(i,j))* &
                                 (DVDZUP(i,j) + abs(DVDZUP(i,j))) + &
                                     R(i,j) * &
                                 (DVDZDN(i,j) + abs(DVDZDN(i,j))) )

              DVDZUP(i,j) = (VDC(i,j,k-1,1) - VDC(i,j,k  ,1))/DZT(i,j,k,bid)
              DVDZDN(i,j) = (VDC(i,j,k  ,1) - VDC(i,j,k+1,1))/WORK2(i,j)
              DIFTP (i,j) = p5*( (c1-R(i,j))* &
                                 (DVDZUP(i,j) + abs(DVDZUP(i,j))) + &
                                     R(i,j) * &
                                 (DVDZDN(i,j) + abs(DVDZDN(i,j))) )

              VISCH(i,j) = VISC(i,j,k)  + VISCP(i,j)*DELHAT(i,j)
              DIFSH(i,j) = VDC(i,j,k,2) + DIFSP(i,j)*DELHAT(i,j)
              DIFTH(i,j) = VDC(i,j,k,1) + DIFTP(i,j)*DELHAT(i,j)

              GAT1(i,j,1) = VISCH(i,j) / HBLT(i,j) /(WM(i,j)+eps)
              DAT1(i,j,1) = -VISCP(i,j)/(WM(i,j)+eps) + F1(i,j)*VISCH(i,j)

              GAT1(i,j,2) = DIFSH(i,j) / HBLT(i,j) /(WS(i,j)+eps)
              DAT1(i,j,2) = -DIFSP(i,j)/(WS(i,j)+eps) + F1(i,j)*DIFSH(i,j)

              GAT1(i,j,3) = DIFTH(i,j) / HBLT(i,j) /(WS(i,j)+eps)
              DAT1(i,j,3) = -DIFTP(i,j)/(WS(i,j)+eps) + F1(i,j)*DIFTH(i,j)

              endif
           end do
           end do

        else

           !DIR$ COLLAPSE
           do j=1,ny_block
           do i=1,nx_block
              if (k == KN(i,j)) then

              DELHAT(i,j) = p5*hwide(k) - zgrid(k) - HBLT(i,j)        
              R     (i,j) = c1 - DELHAT(i,j) / hwide(k)

              DVDZUP(i,j) = (VISC(i,j,k-1) - VISC(i,j,k  ))/hwide(k)
              DVDZDN(i,j) = (VISC(i,j,k  ) - VISC(i,j,k+1))/hwide(k+1)
              VISCP (i,j) = p5*( (c1-R(i,j))* &
                                 (DVDZUP(i,j) + abs(DVDZUP(i,j))) + &
                                     R(i,j) * &
                                 (DVDZDN(i,j) + abs(DVDZDN(i,j))) )

              DVDZUP(i,j) = (VDC(i,j,k-1,2) - VDC(i,j,k  ,2))/hwide(k)
              DVDZDN(i,j) = (VDC(i,j,k  ,2) - VDC(i,j,k+1,2))/hwide(k+1)
              DIFSP (i,j) = p5*( (c1-R(i,j))* &
                                 (DVDZUP(i,j) + abs(DVDZUP(i,j))) + &
                                     R(i,j) * &
                                 (DVDZDN(i,j) + abs(DVDZDN(i,j))) )

              DVDZUP(i,j) = (VDC(i,j,k-1,1) - VDC(i,j,k  ,1))/hwide(k)
              DVDZDN(i,j) = (VDC(i,j,k  ,1) - VDC(i,j,k+1,1))/hwide(k+1)
              DIFTP (i,j) = p5*( (c1-R(i,j))* &
                                 (DVDZUP(i,j) + abs(DVDZUP(i,j))) + &
                                     R(i,j) * &
                                 (DVDZDN(i,j) + abs(DVDZDN(i,j))) )
   
              VISCH(i,j) = VISC(i,j,k)  + VISCP(i,j)*DELHAT(i,j)
              DIFSH(i,j) = VDC(i,j,k,2) + DIFSP(i,j)*DELHAT(i,j)
              DIFTH(i,j) = VDC(i,j,k,1) + DIFTP(i,j)*DELHAT(i,j)

              GAT1(i,j,1) = VISCH(i,j) / HBLT(i,j) /(WM(i,j)+eps)
              DAT1(i,j,1) = -VISCP(i,j)/(WM(i,j)+eps) + &
                             F1(i,j)*VISCH(i,j)

              GAT1(i,j,2) = DIFSH(i,j) / HBLT(i,j) /(WS(i,j)+eps)
              DAT1(i,j,2) = -DIFSP(i,j)/(WS(i,j)+eps) + &
                             F1(i,j)*DIFSH(i,j)

              GAT1(i,j,3) = DIFTH(i,j) / HBLT(i,j) /(WS(i,j)+eps)
              DAT1(i,j,3) = -DIFTP(i,j)/(WS(i,j)+eps) + &
                             F1(i,j)*DIFTH(i,j)

              endif
           end do
           end do

        endif                   ! pbc

     enddo

     DAT1 = min(DAT1,c0)

!-----------------------------------------------------------------------
!
!  compute the dimensionless shape functions and diffusivities
!  at the grid interfaces.  also compute function for non-local
!  transport term (GHAT).
!
!-----------------------------------------------------------------------
#ifdef NOATH
     do k = 1,km       

        if (partial_bottom_cells) then
           if (k > 1) then
              SIGMA = (-zgrid(k-1) + p5*DZT(:,:,k-1,bid) +  &
                       DZT(:,:,k,bid)) / HBLT 
           else
              SIGMA = (-zgrid(k) + p5*hwide(k)) / HBLT     
           end if
        else
           SIGMA = (-zgrid(k) + p5*hwide(k)) / HBLT     
        endif
        F1 = min(SIGMA,epssfc)

        call wscale(F1, HBLT, USTAR, BFSFC, 3, WM, WS)

        !DIR$ COLLAPSE
        do j=1,ny_block
        do i=1,nx_block
           BLMC(i,j,k,1) = HBLT(i,j)*WM(i,j)*SIGMA(i,j)*       &
                           (c1 + SIGMA(i,j)*((SIGMA(i,j)-c2) + &
                           (c3-c2*SIGMA(i,j))*GAT1(i,j,1) +    &
                           (SIGMA(i,j)-c1)*DAT1(i,j,1))) 
           BLMC(i,j,k,2) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*       &
                           (c1 + SIGMA(i,j)*((SIGMA(i,j)-c2) + &
                           (c3-c2*SIGMA(i,j))*GAT1(i,j,2) +    &
                           (SIGMA(i,j)-c1)*DAT1(i,j,2)))
           BLMC(i,j,k,3) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*       &
                           (c1 + SIGMA(i,j)*((SIGMA(i,j)-c2) + &    
                           (c3-c2*SIGMA(i,j))*GAT1(i,j,3) +    &
                           (SIGMA(i,j)-c1)*DAT1(i,j,3)))

           GHAT(i,j,k) = (c1-STABLE(i,j))* cg/(WS(i,j)*HBLT(i,j) +eps)
        end do
        end do

     end do
#endif
      call loop5_blmix(km,ny_block,nx_block,zgrid(1),DZT(1,1,0,bid),&
                       HBLT,hwide(1),BFSFC,USTAR,&
           BLMC(1,1,1,1),GAT1(1,1,1),DAT1(1,1,1),STABLE,GHAT,cg,my_task,p33)

!-----------------------------------------------------------------------
!
!  find diffusivities at kbl-1 grid level
!
!-----------------------------------------------------------------------

     !DIR$ COLLAPSE
     do j=1,ny_block
     do i=1,nx_block
        k = KBL(i,j) - 1
        SIGMA(i,j) = -zgrid(k)/HBLT(i,j)          
     enddo
     enddo

     F1 = min(SIGMA,epssfc)        
     call wscale(F1, HBLT, USTAR, BFSFC, 3, WM, WS)

     !DIR$ COLLAPSE
     do j=1,ny_block
     do i=1,nx_block
        DKM1(i,j,1) = HBLT(i,j)*WM(i,j)*SIGMA(i,j)*     &
                      (c1+SIGMA(i,j)*((SIGMA(i,j)-c2) + &
                      (c3-c2*SIGMA(i,j))*GAT1(i,j,1) +  &
                      (SIGMA(i,j)-c1)*DAT1(i,j,1)))
        DKM1(i,j,2) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*     &
                      (c1+SIGMA(i,j)*((SIGMA(i,j)-c2) + &
                      (c3-c2*SIGMA(i,j))*GAT1(i,j,2) +  &
                      (SIGMA(i,j)-c1)*DAT1(i,j,2)))
        DKM1(i,j,3) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*     &
                      (c1+SIGMA(i,j)*((SIGMA(i,j)-c2) + &       
                      (c3-c2*SIGMA(i,j))*GAT1(i,j,3) +  &
                      (SIGMA(i,j)-c1)*DAT1(i,j,3)))
     end do
     end do

!-----------------------------------------------------------------------
!
!  compute the enhanced mixing
!
!-----------------------------------------------------------------------

     !DIR$ NOVECTOR
     do k=1,km-1

        if (partial_bottom_cells) then
           if (k == 1) then 
              WORK1 = -p5*DZT(:,:,k,bid)
           else
              WORK1 = zgrid(k-1) - p5*(DZT(:,:,k-1,bid) + &
                                       DZT(:,:,k  ,bid))
           end if
           where (k == (KBL - 1)) &
              DELHAT = (HBLT + WORK1)/(p5*(DZT(:,:,k  ,bid) + &
                                           DZT(:,:,k+1,bid)))
        else
           where (k == (KBL - 1)) &
              DELHAT = (HBLT + zgrid(k))/(zgrid(k)-zgrid(k+1))
        endif

        !DIR$ COLLAPSE
        do j=1,ny_block
        do i=1,nx_block
           if (k == (KBL(i,j) - 1)) then

              BLMC(i,j,k,1) = (c1-DELHAT(i,j))*VISC(i,j,k) +             &
                          DELHAT(i,j) *(                                 &
                          (c1-DELHAT(i,j))**2 *DKM1(i,j,1) +             &
                              DELHAT(i,j)**2  *(CASEA(i,j)*VISC(i,j,k) + &
                          (c1-CASEA(i,j))*BLMC(i,j,k,1)))

              BLMC(i,j,k,2) = (c1-DELHAT(i,j))*VDC(i,j,k,2) +            &
                          DELHAT(i,j)*(                                  &
                          (c1-DELHAT(i,j))**2 *DKM1(i,j,2) +             &
                              DELHAT(i,j)**2 *(CASEA(i,j)*VDC(i,j,k,2) + &
                          (c1-CASEA(i,j))*BLMC(i,j,k,2)))

              BLMC(i,j,k,3) = (c1-DELHAT(i,j))*VDC(i,j,k,1) +            &
                          DELHAT(i,j) *(                                 &
                          (c1-DELHAT(i,j))**2 *DKM1(i,j,3) +             &
                              DELHAT(i,j)**2 *(CASEA(i,j)*VDC(i,j,k,1) + &
                          (c1-CASEA(i,j))*BLMC(i,j,k,3)))

              GHAT(i,j,k) = (c1-CASEA(i,j)) * GHAT(i,j,k)

           endif
        end do
        end do
     end do

!-----------------------------------------------------------------------
!
!  combine interior and boundary layer coefficients and nonlocal term
!
!-----------------------------------------------------------------------

     !DIR$ NOVECTOR
     do k=1,km
        !DIR$ NODEP
        do j=1,ny_block
        !DIR$ NODEP
        do i=1,nx_block
           if (k < KBL(i,j)) then 
              VISC(i,j,k)  = BLMC(i,j,k,1)
              VDC(i,j,k,2) = BLMC(i,j,k,2)
              VDC(i,j,k,1) = BLMC(i,j,k,3)
           else
              GHAT(i,j,k) = c0
           endif
        end do
        end do
     enddo
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine blmix









!***********************************************************************
!BOP
! !IROUTINE: wscale
! !INTERFACE:

 subroutine wscale(SIGMA, HBL, USTAR, BFSFC, m_or_s, WM, WS)

! !DESCRIPTION:
!  Computes turbulent velocity scales.
!
!  For $\zeta \geq 0, 
!    w_m = w_s = \kappa U^\star/(1+5\zeta)$
!
!  For $\zeta_m \leq \zeta < 0, 
!    w_m = \kappa U^\star (1-16\zeta)^{1\over 4}$
!
!  For $\zeta_s \leq \zeta < 0, 
!    w_s = \kappa U^\star (1-16\zeta)^{1\over 2}$
!
!  For $\zeta < \zeta_m, 
!    w_m = \kappa U^\star (a_m - c_m\zeta)^{1\over 3}$
!
!  For $\zeta < \zeta_s, 
!    w_s = \kappa U^\star (a_s - c_s\zeta)^{1\over 3}$
!
!  where $\kappa$ is the von Karman constant.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      m_or_s              ! flag =1 for wm only, 2 for ws, 3 for both

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SIGMA,             &! normalized depth (d/hbl)
      HBL,               &! boundary layer depth
      BFSFC,             &! surface buoyancy forcing
      USTAR               ! surface friction velocity

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      WM,                &! turb velocity scales: momentum
      WS                  ! turb velocity scales: tracer

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j  ! dummy loop indices

   real (r8), dimension(nx_block,ny_block) :: &
      ZETA,           &! d/L or sigma*hbl/L(monin-obk)
      ZETAH            ! sigma*hbl*vonkar*BFSFC or ZETA = ZETAH/USTAR**3

!-----------------------------------------------------------------------
!
!  compute zetah and zeta - surface layer is special case
!
!-----------------------------------------------------------------------

   ZETAH = SIGMA*HBL*vonkar*BFSFC
   ZETA  = ZETAH/(USTAR**3 + eps)

!-----------------------------------------------------------------------
!
!  compute velocity scales for momentum
!
!-----------------------------------------------------------------------

   if (m_or_s == 1 .or. m_or_s == 3) then
      do j=1,ny_block
      do i=1,nx_block
         if (ZETA(i,j) >= c0) then ! stable region
            WM(i,j) = vonkar*USTAR(i,j)/(c1 + c5*ZETA(i,j))
         else if (ZETA(i,j) >= zeta_m) then
            WM(i,j) = vonkar*USTAR(i,j)*(c1 - c16*ZETA(i,j))**p25
         else
            WM(i,j) = vonkar*(a_m*(USTAR(i,j)**3)-c_m*ZETAH(i,j))**p33
         endif
      end do
      end do
   endif

!-----------------------------------------------------------------------
!
!  compute velocity scales for tracers
!
!-----------------------------------------------------------------------

   if (m_or_s == 2 .or. m_or_s == 3) then
      do j=1,ny_block
      do i=1,nx_block
         if (ZETA(i,j) >= c0) then
            WS(i,j) = vonkar*USTAR(i,j)/(c1 + c5*ZETA(i,j))
         else if (ZETA(i,j) >= zeta_s) then
            WS(i,j) = vonkar*USTAR(i,j)*SQRT(c1 - c16*ZETA(i,j))
         else
            WS(i,j) = vonkar*(a_s*(USTAR(i,j)**3)-c_s*ZETAH(i,j))**p33
         endif
      end do
      end do
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine wscale











!***********************************************************************
!BOP
! !IROUTINE: ddmix
! !INTERFACE:

 subroutine ddmix(VDC, TRCR, this_block)

! !DESCRIPTION:
!  $R_\rho$ dependent interior flux parameterization.
!  Add double-diffusion diffusivities to Ri-mix values at blending
!  interface and below.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TRCR                ! tracers at current time

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,0:km+1,2),intent(inout) :: &
      VDC        ! diffusivity for tracer diffusion

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  k,kup,knxt

   real (r8), dimension(nx_block,ny_block) :: &
      ALPHADT,           &! alpha*DT  across interfaces
      BETADS,            &! beta *DS  across interfaces
      RRHO,              &! dd density ratio
      DIFFDD,            &! dd diffusivity scale
      PRANDTL             ! prandtl number

   real (r8), dimension(nx_block,ny_block,2) :: &
      TALPHA,            &! temperature expansion coeff
      SBETA               ! salinity    expansion coeff

!-----added for athread by zengyh 2018-10-23----
   integer::myknxt,mykup,mynx_block,myny_block,mykm  


!-----test for athread validation by zengyh 2018-10-23----
#ifdef POPCHKJN
   integer::i,j,n
   real(8),allocatable,dimension(:,:,:,:)::tmp_VDC
   real(8)::absVDCerr=0.0
   real(8)::tmpVDC
   real(8),allocatable, dimension(:)::aveVDC

   real(8),allocatable, dimension(:,:,:,:) :: &
      MWORK21,SWORK21
   real(8),allocatable, dimension(:,:,:) :: &
      MWORK11,SWORK11,             &
      MWORK1KM,SWORK1KM,MWORK2KM,SWORK2KM,         &
      MWORK1K,MWORK2K,SWORK1K,SWORK2K           
   real(8)::maxrho1err=1.0e-16,maxaerr=1.0e-16,maxberr=1.0e-16
   real(8)::maxdterr=1.0e-16,maxdserr=1.0e-16
   real(8)::maxvdcerr=1.0e-16,maxvdc2err=1.0e-16,tmperr,tmp2err
   real(8)::avevdcerr=0.0
   real(8)::posrho1me,posame,posbme,posdtme,posdsme,posvdcme,posvdc2me
#endif




   Interface
     subroutine m_kpp_ddmix_loop1(TRCR, VDC, myproc, nx_block, ny_block, km, nt,&
            knxt, kup, tmax, tmin, smax, smin, pressz, grav,&
            Rrho0, dsfmax &
#if defined(POPCHKJN) || defined(POPTIMEJN)
            ,steps&
#endif
#if defined(POPCHKJN)
            ,SWORK11,SWORK21,SWORK1K,SWORK2K,SWORK1KM,SWORK2KM&
#endif
            ) bind(C, name='m_kpp_ddmix_loop1') 
!         use iso_c_binding
         real(8):: TRCR(*)
         integer,value:: myproc, nx_block, ny_block, km, nt, knxt, kup
         real(8):: tmax(*), tmin(*), smax(*), smin(*), pressz(*)
         real(8),value:: grav, Rrho0, dsfmax
!         real(8),intent(inout):: VDC(*)
         real(8):: VDC(*)
#if defined(POPCHKJN) || defined(POPTIMEJN)
         integer,value:: steps
#endif
#if defined(POPCHKJN)
         real(8),intent(out):: SWORK11(*), SWORK1K(*), SWORK1KM(*) 
         real(8),intent(out):: SWORK21(*), SWORK2K(*), SWORK2KM(*) 
#endif
     end subroutine m_kpp_ddmix_loop1
   end interface


!-----------------------------------------------------------------------
!
!  compute alpha*DT and beta*DS at interfaces.  use RRHO and
!  PRANDTL for temporary storage for call to state
!
!-----------------------------------------------------------------------

#ifdef POPCHKJN
  allocate(aveVDC(0:km+1))
  allocate(tmp_VDC(nx_block,ny_block,0:km+1,2))

  allocate(MWORK11(nx_block,ny_block,km))
  allocate(MWORK1K(nx_block,ny_block,km))
  allocate(MWORK1KM(nx_block,ny_block,km))
  allocate(MWORK21(nx_block,ny_block,0:km+1,2))
  allocate(MWORK2K(nx_block,ny_block,km))
  allocate(MWORK2KM(nx_block,ny_block,km))

  allocate(SWORK11(km,nx_block,ny_block))
  allocate(SWORK1K(km,nx_block,ny_block))
  allocate(SWORK1KM(km,nx_block,ny_block))
  allocate(SWORK21(2,0:km+1,nx_block,ny_block))
  allocate(SWORK2K(km,nx_block,ny_block))
  allocate(SWORK2KM(km,nx_block,ny_block))
  
  tmp_VDC=VDC
#endif


   kup  = 1
   knxt = 2

#if defined(POPTIMEJN) || defined(POPCHKJN)    
!------- timing_kpp instrumented by zengyh 2018-10-23:----
  call getlooptime(start_kpp_state)

   PRANDTL = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)

   call state(1, 1, PRANDTL, TRCR(:,:,1,2), this_block, &
                    RHOFULL=RRHO, &
                    DRHODT=TALPHA(:,:,kup), DRHODS=SBETA(:,:,kup))

#if defined(POPCHKJN) 
   MWORK11(:,:,1)=RRHO
   MWORK1K(:,:,1)=TALPHA(:,:,kup)
   MWORK1KM(:,:,1)=SBETA(:,:,kup)
#endif

   do k=1,km

      if ( k < km ) then

         PRANDTL = merge(-c2,TRCR(:,:,k+1,1),TRCR(:,:,k+1,1) < -c2)

         call state(k+1, k+1, PRANDTL, TRCR(:,:,k+1,2),              &
                              this_block,                            &
                              RHOFULL=RRHO, DRHODT=TALPHA(:,:,knxt), &
                                            DRHODS= SBETA(:,:,knxt))

#if defined(POPCHKJN) 
!         MWORK11(:,:,k+1)=RRHO
         MWORK1K(:,:,k+1)=TALPHA(:,:,knxt)
         MWORK1KM(:,:,k+1)=SBETA(:,:,knxt)
#endif

         ALPHADT = -p5*(TALPHA(:,:,kup) + TALPHA(:,:,knxt)) &
                      *(TRCR(:,:,k,1) - TRCR(:,:,k+1,1))

         BETADS  = p5*( SBETA(:,:,kup) +  SBETA(:,:,knxt)) &
                     *(TRCR(:,:,k,2) - TRCR(:,:,k+1,2))

         kup  = knxt
         knxt = 3 - kup

      else

         ALPHADT = c0
         BETADS  = c0

      endif       

#if defined(POPCHKJN) 
!      MWORK2K(:,:,k)=ALPHADT
!      MWORK2KM(:,:,k)=BETADS
#endif

!-----------------------------------------------------------------------
!
!     salt fingering case
!
!-----------------------------------------------------------------------

#if defined(POPCHKJN) 
      DIFFDD=c0
#endif

      where ( ALPHADT > BETADS .and. BETADS > c0 )

         RRHO       = MIN(ALPHADT/BETADS, Rrho0)
         DIFFDD     = dsfmax*(c1-(RRHO-c1)/(Rrho0-c1))**3
         VDC(:,:,k,1) = VDC(:,:,k,1) + 0.7_r8*DIFFDD
         VDC(:,:,k,2) = VDC(:,:,k,2) + DIFFDD

      endwhere

#if defined(POPCHKJN) 
      MWORK11(:,:,k)=DIFFDD
      MWORK1K(:,:,k)=c0
      MWORK1KM(:,:,k)=c0
      MWORK2KM(:,:,k)=c0
#endif


!-----------------------------------------------------------------------
!
!     diffusive convection
!
!-----------------------------------------------------------------------

      where ( ALPHADT < c0 .and. BETADS < c0 .and. ALPHADT > BETADS )
         RRHO    = ALPHADT / BETADS
         DIFFDD  = 1.5e-2_r8*0.909_r8* &
                   exp(4.6_r8*exp(-0.54_r8*(c1/RRHO-c1)))
         PRANDTL = 0.15_r8*RRHO

#if defined(POPCHKJN) 
      MWORK1K(:,:,k)=exp(-0.54_r8*(c1/rrho-c1))
      MWORK1KM(:,:,k)=exp(4.6_r8*MWORK1K(:,:,k))
      MWORK2K(:,:,k)=DIFFDD
!      MWORK2KM(:,:,k)=PRANDTL
      MWORK2KM(:,:,k)=-0.54_r8*(c1/rrho-c1)
#endif

      elsewhere
         RRHO    = c0
         DIFFDD  = c0
         PRANDTL = c0
      endwhere

      where (RRHO > p5) PRANDTL = (1.85_r8 - 0.85_r8/RRHO)*RRHO

      VDC(:,:,k,1) = VDC(:,:,k,1) + DIFFDD
      VDC(:,:,k,2) = VDC(:,:,k,2) + PRANDTL*DIFFDD

   end do

!------- timing_kpp instrumented by zengyh 2018-10-23:----
  call getlooptime(stop_kpp_state)
  if (my_task == master_task) then
     WRITE (*,*) 'kpp_ddmix costs zengyh: ',stop_kpp_state-start_kpp_state,' (us);'
  endif
#endif  


#ifdef POPCHKJN    
!----------store in tmp_VDC for validation by zengyh 2018-10-23:-----!
   MWORK21=VDC
   VDC=tmp_VDC
   tmp_VDC=MWORK21

   num_steps=num_steps+1
#endif  


!--------Athreaded by zengyh 2018-10-23:-------

#ifdef POPTIMEJN    
!------- timing_kpp instrumented by zengyh 2018-09-03:----
  call getlooptime(start_kpp_state)
#endif  

   myknxt = knxt-1
   mykup = kup-1
   mynx_block = nx_block
   myny_block = ny_block
   mykm = km

!  if (my_task == master_task .and. num_steps == 0) then
!     WRITE (*,*) 'some loc paras: ',myknxt,mykup,mynx_block,myny_block,mykm,nt,';'
!     WRITE (*,*) 'size of k type:', sizeof(k),';'
!     WRITE (*,*) 'press:', pressz,';'
!  endif

!-----------------------------------------------------------------------
!
!  calculate VDC for all other levels
!
!-----------------------------------------------------------------------

if (state_itype==state_type_mwjf.and.state_range_iopt==state_range_enforce) then
! apply new athreaded codes:

! state called with present(RHOFULL).
  call m_kpp_ddmix_loop1(TRCR,VDC,my_task,&
          mynx_block,myny_block,mykm,2,myknxt,mykup, &
          tmax,tmin,smax,smin,pressz,grav, &
          Rrho0, dsfmax &
#if defined(POPCHKJN) || defined(POPTIMEJN)
          ,num_steps &
#endif
#if defined(POPCHKJN) 
          ,SWORK11,SWORK21,SWORK1K,SWORK2K,SWORK1KM,SWORK2KM &
#endif    
          ) 

else		  
! apply initial codes:	

   kup  = 1
   knxt = 2

   PRANDTL = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)

   call state(1, 1, PRANDTL, TRCR(:,:,1,2), this_block, &
                    RHOFULL=RRHO, &
                    DRHODT=TALPHA(:,:,kup), DRHODS=SBETA(:,:,kup))

   do k=1,km

      if ( k < km ) then

         PRANDTL = merge(-c2,TRCR(:,:,k+1,1),TRCR(:,:,k+1,1) < -c2)

         call state(k+1, k+1, PRANDTL, TRCR(:,:,k+1,2),              &
                              this_block,                            &
                              RHOFULL=RRHO, DRHODT=TALPHA(:,:,knxt), &
                                            DRHODS= SBETA(:,:,knxt))

         ALPHADT = -p5*(TALPHA(:,:,kup) + TALPHA(:,:,knxt)) &
                      *(TRCR(:,:,k,1) - TRCR(:,:,k+1,1))

         BETADS  = p5*( SBETA(:,:,kup) +  SBETA(:,:,knxt)) &
                     *(TRCR(:,:,k,2) - TRCR(:,:,k+1,2))

         kup  = knxt
         knxt = 3 - kup

      else

         ALPHADT = c0
         BETADS  = c0

      endif       

!-----------------------------------------------------------------------
!
!     salt fingering case
!
!-----------------------------------------------------------------------

      where ( ALPHADT > BETADS .and. BETADS > c0 )

         RRHO       = MIN(ALPHADT/BETADS, Rrho0)
         DIFFDD     = dsfmax*(c1-(RRHO-c1)/(Rrho0-c1))**3
         VDC(:,:,k,1) = VDC(:,:,k,1) + 0.7_r8*DIFFDD
         VDC(:,:,k,2) = VDC(:,:,k,2) + DIFFDD

      endwhere

!-----------------------------------------------------------------------
!
!     diffusive convection
!
!-----------------------------------------------------------------------

      where ( ALPHADT < c0 .and. BETADS < c0 .and. ALPHADT > BETADS )
         RRHO    = ALPHADT / BETADS
         DIFFDD  = 1.5e-2_r8*0.909_r8* &
                   exp(4.6_r8*exp(-0.54_r8*(c1/RRHO-c1)))
         PRANDTL = 0.15_r8*RRHO
      elsewhere
         RRHO    = c0
         DIFFDD  = c0
         PRANDTL = c0
      endwhere

      where (RRHO > p5) PRANDTL = (1.85_r8 - 0.85_r8/RRHO)*RRHO

      VDC(:,:,k,1) = VDC(:,:,k,1) + DIFFDD
      VDC(:,:,k,2) = VDC(:,:,k,2) + PRANDTL*DIFFDD

   end do
  
	
endif


#ifdef POPTIMEJN    
!------- timing_kpp instrumented by zengyh 2018-09-03:----
  call getlooptime(stop_kpp_state)

!  if (my_task == 0 .and. mod(num_steps,8) == 1) then
  if (my_task == 0) then
     WRITE (*,*) 'kpp_ddmix costs w athread zengyh: ',stop_kpp_state-start_kpp_state,' (us);'
  endif
#endif




#ifdef POPCHKJN    
 if (my_task == 0 .and. mod(num_steps,8)==1) then
    maxrho1err=1.0e-16
    maxaerr=1.0e-16
    maxberr=1.0e-16
    maxdterr=1.0e-16
    maxdserr=1.0e-16
    maxvdcerr=1.0e-16
    maxvdc2err=1.0e-16
    avevdcerr=0.0
    do k=1,km
       do j=1,ny_block
          do i=1,nx_block
             tmperr=abs(MWORK11(i,j,k)-SWORK11(k,i,j))
             if (tmperr>maxrho1err) then 
                maxrho1err=tmperr
                posrho1me=MWORK11(i,j,k) 
             endif
             tmperr=abs(MWORK1K(i,j,k)-SWORK1K(k,i,j))
             if (tmperr>maxaerr) then 
                maxaerr=tmperr
                posame=MWORK1K(i,j,k) 
             endif
             tmperr=abs(MWORK1KM(i,j,k)-SWORK1KM(k,i,j))
             if (tmperr>maxberr) then 
                maxberr=tmperr
                posbme=MWORK1KM(i,j,k) 
             endif
             tmperr=abs(MWORK2K(i,j,k)-SWORK2K(k,i,j))
             if (tmperr>maxdterr) then 
                maxdterr=tmperr
                posdtme=MWORK2K(i,j,k) 
             endif
             tmperr=abs(MWORK2KM(i,j,k)-SWORK2KM(k,i,j))
             if (tmperr>maxdserr) then 
                maxdserr=tmperr
                posdsme=MWORK2KM(i,j,k) 
             endif
          enddo
      enddo
    enddo
    do k=0,km+1
       do j=1,ny_block
          do i=1,nx_block
             tmperr =abs(MWORK21(i,j,k,1)-SWORK21(1,k,i,j))
             if (tmperr>maxvdcerr) then 
                maxvdcerr=tmperr
                posvdcme=MWORK21(i,j,k,1) 
             endif
             avevdcerr=avevdcerr+tmperr
             tmp2err=abs(MWORK21(i,j,k,2)-SWORK21(2,k,i,j))
             if (tmp2err>maxvdc2err) then 
                maxvdc2err=tmp2err
                posvdc2me=MWORK21(i,j,k,2) 
             endif
             avevdcerr=avevdcerr+tmp2err
          enddo
      enddo
    enddo
    avevdcerr=avevdcerr/(nx_block*ny_block*2*(km+2))
!    WRITE (*,*) 'RHO1 max absolute err: ',maxrho1err,'@',posrho1me, '(kpp_VDC);'
!    WRITE (*,*) 'ALPHA max absolute err: ',maxaerr,'@',posame, '(kpp_VDC);'
!    WRITE (*,*) 'BETA max absolute err: ',maxberr,'@',posbme, '(kpp_VDC);'
    WRITE (*,*) 'exp1 max absolute err: ',maxaerr,'@',posame, '(kpp_VDC);'
    WRITE (*,*) 'exp2 max absolute err: ',maxberr,'@',posbme, '(kpp_VDC);'
!    WRITE (*,*) 'ALPHADT max absolute err: ',maxdterr,'@',posdtme, '(kpp_VDC);'
!    WRITE (*,*) 'BETADS max absolute err: ',maxdserr,'@',posdsme, '(kpp_VDC);'
    WRITE (*,*) 'SALTDD max absolute err: ',maxrho1err,'@',posrho1me, '(kpp_VDC);'
    WRITE (*,*) 'DIFFDD max absolute err: ',maxdterr,'@',posdtme, '(kpp_VDC);'
!    WRITE (*,*) 'PRANDTL max absolute err: ',maxdserr,'@',posdsme, '(kpp_VDC);'
    WRITE (*,*) 'exp1_x max absolute err: ',maxdserr,'@',posdsme, '(kpp_VDC);'
    WRITE (*,*) 'VDC max absolute err: ',maxvdcerr,'@',posvdcme, '(kpp_VDC);'
    WRITE (*,*) 'VDC2 max absolute err: ',maxvdc2err,'@',posvdc2me, '(kpp_VDC);'
    WRITE (*,*) 'VDC ave absolute err: ',avevdcerr, '(kpp_VDC);'
 endif



!--------test absolute error of kpp_VDC for validation by zengyh 2018-10-26:----!
  absVDCerr=0.0
  do n=1,2
  do k=0,km+1
     aveVDC(k)=sum(abs(tmp_VDC(:,:,k,n)))/(nx_block*ny_block)
     do j=1,ny_block
        do i=1,nx_block
           tmpVDC=tmp_VDC(i,j,k,n)
           tmp_VDC(i,j,k,n)=abs(VDC(i,j,k,n)-tmp_VDC(i,j,k,n))
           if(tmpVDC/=c0) then
             tmp_VDC(i,j,k,n)=tmp_VDC(i,j,k,n)/abs(tmpVDC)
             absVDCerr=absVDCerr+tmp_VDC(i,j,k,n)
           endif
        end do
     end do
  end do
  end do

   if (my_task == master_task .and. mod(num_steps,8) == 1) then
      WRITE (*,*) 'kpp_state mean of VDC: ',sum(aveVDC)/((km+2)*2),'(VDC)'
      WRITE (*,*) 'kpp_state absolute sum relative error of athread w zengyh: ',absVDCerr,'(VDC)'
      WRITE (*,*) 'max abs relative error value:',maxval(tmp_VDC(:,:,:,:)),'(VDC)'
      WRITE (*,*) 'min abs relatibe error:',minval(tmp_VDC(:,:,:,:)),'(VDC)'
      WRITE (*,*) 'ave abs relative error:',absVDCerr/(nx_block*ny_block*(km+2)*2),'(VDC)'
   endif


 deallocate(aveVDC)
 deallocate(tmp_VDC)

 deallocate(MWORK11)
 deallocate(MWORK1K)
 deallocate(MWORK1KM)
 deallocate(MWORK21)
 deallocate(MWORK2K)
 deallocate(MWORK2KM)

 deallocate(SWORK11)
 deallocate(SWORK1K)
 deallocate(SWORK1KM)
 deallocate(SWORK21)
 deallocate(SWORK2K)
 deallocate(SWORK2KM)
#endif


!-----------------------------------------------------------------------
!EOC

 end subroutine ddmix












!***********************************************************************
!BOP
! !IROUTINE: buoydiff
! !INTERFACE:

 subroutine buoydiff(DBLOC, DBSFC, TRCR, this_block)

! !DESCRIPTION:
!  This routine calculates the buoyancy differences at model levels.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TRCR                ! tracers at current time

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km), intent(out) :: & 
      DBLOC,         &! buoyancy difference between adjacent levels
      DBSFC           ! buoyancy difference between level and surface

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k,                 &! vertical level index
      i,j,               &! horizontal indices
      kprev, klvl, ktmp, &! indices for 2-level TEMPK array
      bid                 ! local block index

   real (r8), dimension(nx_block,ny_block) :: &
      RHO1,              &! density of sfc t,s displaced to k
      RHOKM,             &! density of t(k-1),s(k-1) displaced to k
      RHOK,              &! density at level k
      TEMPSFC            ! adjusted temperature at surface
!      TALPHA,            &! temperature expansion coefficient
!      SBETA               ! salinity    expansion coefficient

   real (r8), dimension(nx_block,ny_block,2) :: &
      TEMPK               ! temp adjusted for freeze at levels k,k-1



	  

!-----added for athread by zengyh 2018-09-03----
   integer::myklvl,mykprev,mynx_block,myny_block,mykm  
   integer,allocatable,dimension(:,:)::LOCKMT    


!-----test for athread validation by zengyh 2018-09-03----
#ifdef POPCHKJN
   real(8),allocatable,dimension(:,:,:)::tmp_DBLOC, tmp_DBSFC
   real(8)::absLOCerr=0.0, absSFCerr=0.0
   real(8)::tmpLOC,tmpSFC
   real(8),allocatable, dimension(:)::aveLOC,aveSFC

   real(8),allocatable, dimension(:,:,:) :: &
      MWORK11,MWORK21,SWORK11,SWORK21,             &
      MWORK1KM,SWORK1KM,           &
      MWORK1K,MWORK2K,SWORK1K,SWORK2K           
   real(8)::maxrho1err=1.0e-16,maxrhokerr=1.0e-16,maxrhokmerr=1.0e-16
   real(8)::maxlocerr=1.0e-16,maxsfcerr=1.0e-16,tmperr
   real(8)::avelocerr=0.0,avesfcerr=0.0
   real(8)::posrho1me,posrhokme,posrhokmme,poslocme,possfcme
#endif




   Interface
     subroutine m_kpp_state_loop2(TRCR, KMT, myproc, nx_block, ny_block, km, nt,&
            klvl, kprev, tmax, tmin, smax, smin, pressz, grav, DBLOC, DBSFC&
#if defined(POPCHKJN) || defined(POPTIMEJN)
            ,steps&
#endif
#if defined(POPCHKJN)
            ,SWORK11,SWORK21,SWORK1K,SWORK2K,SWORK1KM&
#endif
            ) bind(C, name='m_kpp_state_loop2') 
!         use iso_c_binding
         real(8):: TRCR(*)
         integer,value:: myproc, nx_block, ny_block, km, nt, klvl, kprev
         integer:: KMT(*) 
         real(8):: tmax(*), tmin(*), smax(*), smin(*), pressz(*)
         real(8),value:: grav
         real(8),intent(out):: DBLOC(*), DBSFC(*)
#if defined(POPCHKJN) || defined(POPTIMEJN)
         integer,value:: steps
#endif
#if defined(POPCHKJN)
         real(8),intent(out):: SWORK11(*), SWORK1K(*), SWORK1KM(*) 
         real(8),intent(out):: SWORK21(*), SWORK2K(*) 
#endif
     end subroutine m_kpp_state_loop2
   end interface



  allocate(LOCKMT(nx_block,ny_block))    

#ifdef POPCHKJN
  allocate(aveLOC(km))
  allocate(aveSFC(km))
  allocate(tmp_DBLOC(nx_block,ny_block,km))
  allocate(tmp_DBSFC(nx_block,ny_block,km))

  allocate(MWORK11(nx_block,ny_block,km))
  allocate(MWORK1K(nx_block,ny_block,km))
  allocate(MWORK1KM(nx_block,ny_block,km))
  allocate(MWORK21(nx_block,ny_block,km))
  allocate(MWORK2K(nx_block,ny_block,km))
  allocate(SWORK11(km,nx_block,ny_block))
  allocate(SWORK1K(km,nx_block,ny_block))
  allocate(SWORK1KM(km,nx_block,ny_block))
  allocate(SWORK21(km,nx_block,ny_block))
  allocate(SWORK2K(km,nx_block,ny_block))
#endif




!-----------------------------------------------------------------------
!
!  calculate density and buoyancy differences at surface
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   klvl  = 2
   kprev = 1




#if defined(POPTIMEJN) || defined(POPCHKJN)    
!------- timing_kpp instrumented by zengyh 2018-09-03:----
  call getlooptime(start_kpp_state)

   TEMPSFC = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)
   TEMPK(:,:,kprev) = TEMPSFC
   DBSFC(:,:,1) = c0

   do k = 2,km

      TEMPK(:,:,klvl) = merge(-c2,TRCR(:,:,k,1),TRCR(:,:,k,1) < -c2)

      call state(k, k, TEMPSFC,          TRCR(:,:,1  ,2), &
        this_block, RHOFULL=RHO1)
      call state(k, k, TEMPK(:,:,kprev), TRCR(:,:,k-1,2), &
        this_block, RHOFULL=RHOKM)
      call state(k, k, TEMPK(:,:,klvl),  TRCR(:,:,k  ,2), &
        this_block, RHOFULL=RHOK)


      do j=1,ny_block
      do i=1,nx_block
         if (RHOK(i,j) /= c0) then
            DBSFC(i,j,k)   = grav*(c1 - RHO1 (i,j)/RHOK(i,j))
            DBLOC(i,j,k-1) = grav*(c1 - RHOKM(i,j)/RHOK(i,j))
         else
            DBSFC(i,j,k)   = c0
            DBLOC(i,j,k-1) = c0
         endif

         if (k-1 >= KMT(i,j,bid)) DBLOC(i,j,k-1) = c0
      end do
      end do
#endif  

#ifdef POPCHKJN    
      MWORK11(:,:,k)=RHO1(:,:)
      MWORK1K(:,:,k)=RHOK(:,:)
      MWORK1KM(:,:,k)=RHOKM(:,:)
#endif  

#if defined(POPTIMEJN) || defined(POPCHKJN)    
      ktmp  = klvl
      klvl  = kprev
      kprev = ktmp

   enddo

   DBLOC(:,:,km) = c0
 
!------- timing_kpp instrumented by zengyh 2018-09-03:----
  call getlooptime(stop_kpp_state)
  if (my_task == master_task) then
     WRITE (*,*) 'kpp_state costs zengyh: ',stop_kpp_state-start_kpp_state,' (us);'
  endif
#endif  


#ifdef POPCHKJN    
   MWORK21(:,:,:)=DBLOC(:,:,:)
   MWORK2K(:,:,:)=DBSFC(:,:,:)

!----------store in tmp_LOC & SFC for validation by zengyh 2018-09-03:-----!
   tmp_DBLOC(:,:,:)=DBLOC(:,:,:)
   tmp_DBSFC(:,:,:)=DBSFC(:,:,:)
#endif  





!--------Athreaded by zengyh 2018-09-03:-------

#ifdef POPTIMEJN    
!------- timing_kpp instrumented by zengyh 2018-09-03:----
  call getlooptime(start_kpp_state)
#endif  

   klvl  = 2
   kprev = 1


   LOCKMT = KMT(:,:,bid)-1
   myklvl = klvl-1
   mykprev = kprev-1
   mynx_block = nx_block
   myny_block = ny_block
   mykm = km


!  if (my_task == master_task .and. num_steps == 1) then
!     WRITE (*,*) 'some loc paras: ',myklvl,mykprev,mynx_block,myny_block,mykm,nt,';'
!     WRITE (*,*) 'size of k type:', sizeof(k),';'
!     WRITE (*,*) 'press:', pressz,';'
!  endif



!-----------------------------------------------------------------------
!
!  calculate DBLOC and DBSFC for all other levels
!
!-----------------------------------------------------------------------


if (state_itype==state_type_mwjf.and.state_range_iopt==state_range_enforce) then
! apply new athreaded codes:

! state called with present(RHOFULL).
  call m_kpp_state_loop2(TRCR,LOCKMT,my_task,&
          mynx_block,myny_block,mykm,2,myklvl,mykprev, &
          tmax,tmin,smax,smin,pressz,grav,DBLOC,DBSFC &
#if defined(POPCHKJN) || defined(POPTIMEJN)
          ,num_steps &
#endif
#if defined(POPCHKJN) 
          ,SWORK11,SWORK21,SWORK1K,SWORK2K,SWORK1KM &
#endif    
          ) 

else
! apply initial codes:

   TEMPSFC = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)
   TEMPK(:,:,kprev) = TEMPSFC
   DBSFC(:,:,1) = c0

   do k = 2,km
      TEMPK(:,:,klvl) = merge(-c2,TRCR(:,:,k,1),TRCR(:,:,k,1) < -c2)
      call state(k, k, TEMPSFC,          TRCR(:,:,1  ,2), &
                       this_block, RHOFULL=RHO1)
      call state(k, k, TEMPK(:,:,kprev), TRCR(:,:,k-1,2), &
                       this_block, RHOFULL=RHOKM)
      call state(k, k, TEMPK(:,:,klvl),  TRCR(:,:,k  ,2), &
                       this_block, RHOFULL=RHOK)
!      subroutine state(k, kk, TEMPK, SALTK, this_block, &
!                         RHOOUT, RHOFULL, DRHODT, DRHODS)

      do j=1,ny_block
      do i=1,nx_block
         if (RHOK(i,j) /= c0) then
            DBSFC(i,j,k)   = grav*(c1 - RHO1 (i,j)/RHOK(i,j))
            DBLOC(i,j,k-1) = grav*(c1 - RHOKM(i,j)/RHOK(i,j))
         else
            DBSFC(i,j,k)   = c0
            DBLOC(i,j,k-1) = c0
         endif
         if (k-1 >= KMT(i,j,bid)) DBLOC(i,j,k-1) = c0
      end do
      end do
      ktmp  = klvl
      klvl  = kprev
      kprev = ktmp
   enddo

   DBLOC(:,:,km) = c0

endif


#ifdef POPTIMEJN    
!------- timing_kpp instrumented by zengyh 2018-09-03:----
  call getlooptime(stop_kpp_state)

!  if (my_task == 0 .and. mod(num_steps,8) == 1) then
  if (my_task == 0) then
     WRITE (*,*) 'kpp_state costs w athread zengyh: ',stop_kpp_state-start_kpp_state,' (us);'
  endif
#endif




#ifdef POPCHKJN    
 if (my_task == -1 .and. num_steps<=2) then
    open(2,file="WORK_f90.txt",access='append')
    do k=2,km
         do j=1,ny_block
            do i=1,nx_block
!               write(2,"('c=',i1,' j=',i3,' i=',i2,' k=',i2,f20.16,f20.16,f20.16)") &  
! left alligned
               write(2,"(i4,i4,i4,i4,f21.16,f21.16,f21.16,f21.16,f21.16)") &
!               write(2,"(1x,'c=',i1,' j=',i3,' i=',i2,' k=',i2,f20.16,f20.16,f20.16,f20.16,f20.16,f20.16)") &   
! right alligned
                num_steps,j-1,i-1,k-1, &
                MWORK11(i,j,k)-SWORK11(k,i,j),MWORK21(i,j,k)-SWORK21(k,i,j),&
                MWORK1K(i,j,k)-SWORK1K(k,i,j),MWORK2K(i,j,k)-SWORK2K(k,i,j),&
                MWORK1KM(i,j,k)-SWORK1KM(k,i,j)
!               MWORK21(i,j,k),MWORK21(i,j,k)-SWORK21(k,i,j),&
!               MWORK2K(i,j,k),MWORK2K(i,j,k)-SWORK2K(k,i,j),&
!               MWORK1KM(i,j,k)-SWORK1KM(k,i,j)
            enddo
         enddo
    enddo
    close(2)
 endif


 if (my_task == 0 .and. mod(num_steps,8)==0) then
    maxrho1err=1.0e-16
    maxrhokerr=1.0e-16
    maxrhokmerr=1.0e-16
    maxlocerr=1.0e-16
    maxsfcerr=1.0e-16
    avelocerr=0.0
    avesfcerr=0.0
    do k=2,km
       do j=1,ny_block
          do i=1,nx_block
             tmperr=abs(MWORK11(i,j,k)-SWORK11(k,i,j))
             if (tmperr>maxrho1err) then 
                maxrho1err=tmperr
                posrho1me=MWORK11(i,j,k) 
             endif
             tmperr=abs(MWORK1K(i,j,k)-SWORK1K(k,i,j))
             if (tmperr>maxrhokerr) then
                maxrhokerr=tmperr 
                posrhokme=MWORK1K(i,j,k) 
             endif
             tmperr=abs(MWORK1KM(i,j,k)-SWORK1KM(k,i,j))
             if (tmperr>maxrhokmerr) then
                maxrhokmerr=tmperr 
                posrhokmme=MWORK1KM(i,j,k) 
             endif
             tmperr=abs(MWORK21(i,j,k)-SWORK21(k,i,j))
             avelocerr=avelocerr+tmperr
             if (tmperr>maxlocerr) then
                maxlocerr=tmperr 
                poslocme=MWORK21(i,j,k) 
             endif 
             tmperr=abs(MWORK2K(i,j,k)-SWORK2K(k,i,j))
             avesfcerr=avesfcerr+tmperr
             if (tmperr>maxsfcerr) then
                maxsfcerr=tmperr 
                possfcme=MWORK2k(i,j,k) 
             endif
          enddo
      enddo
    enddo
    avelocerr=avelocerr/(nx_block*ny_block*km)
    avesfcerr=avesfcerr/(nx_block*ny_block*km)
    WRITE (*,*) 'RHO1 max absolute err: ',maxrho1err,'@',posrho1me, '(kpp_RHO);'
    WRITE (*,*) 'RHOk max absolute err: ',maxrhokerr,'@',posrhokme, '(kpp_RHO);'
    WRITE (*,*) 'RHOkm max absolute err: ',maxrhokmerr,'@',posrhokmme, '(kpp_RHO);'
    WRITE (*,*) 'DBLOC max absolute err: ',maxlocerr,'@',poslocme, '(kpp_RHO);'
    WRITE (*,*) 'DBLOC ave absolute err: ',avelocerr, '(kpp_RHO);'
    WRITE (*,*) 'DBSFC max absolute err: ',maxsfcerr,'@',possfcme, '(kpp_RHO);'
    WRITE (*,*) 'DBSFC ave absolute err: ',avesfcerr, '(kpp_RHO);'
 endif



!--------test absolute error of kpp_LOC & SFC for validation by zengyh 2018-09-03:-------!
  absLOCerr=0.0
  absSFCerr=0.0
  do k=1,km
     aveLOC(k)=sum(abs(tmp_DBLOC(:,:,k)))/(nx_block*ny_block)
     aveSFC(k)=sum(abs(tmp_DBSFC(:,:,k)))/(nx_block*ny_block)
     do j=1,ny_block
        do i=1,nx_block
           tmpLOC=tmp_DBLOC(i,j,k)
           tmpSFC=tmp_DBSFC(i,j,k)
           tmp_DBLOC(i,j,k)=abs(DBLOC(i,j,k)-tmp_DBLOC(i,j,k))
           tmp_DBSFC(i,j,k)=abs(DBSFC(i,j,k)-tmp_DBSFC(i,j,k)) 
           if(tmpLOC/=c0) then
             tmp_DBLOC(i,j,k)=tmp_DBLOC(i,j,k)/abs(tmpLOC)
             absLOCerr=absLOCerr+tmp_DBLOC(i,j,k)
           endif
           if(tmpSFC/=c0) then
             tmp_DBSFC(i,j,k)=tmp_DBSFC(i,j,k)/abs(tmpSFC) 
             absSFCerr=absSFCerr+tmp_DBSFC(i,j,k)
           endif
        end do
     end do
  end do

   if (my_task == master_task .and. mod(num_steps,8) == 1) then
      WRITE (*,*) 'kpp_state mean of DBLOC: ',sum(aveLOC)/km,'(DBLOC)'
      WRITE (*,*) 'kpp_state absolute sum relative error of athread w zengyh: ',absLOCerr,'(DBLOC)'
      WRITE (*,*) 'max abs relative error value:',maxval(tmp_DBLOC(:,:,:)),'(DBLOC)'
      WRITE (*,*) 'min abs relatibe error:',minval(tmp_DBLOC(:,:,:)),'(DBLOC)'
      WRITE (*,*) 'ave abs relative error:',absLOCerr/(nx_block*ny_block*km),'(DBLOC)'
      WRITE (*,*) 'kpp_state mean of DBSFC: ',sum(aveSFC)/km,'(DBSFC)'
      WRITE (*,*) 'kpp_state absolute sum relative error of athread w zengyh: ',absSFCerr,'(DBSFC)'
      WRITE (*,*) 'max abs relative error value:',maxval(tmp_DBSFC(:,:,:)),'(DBSFC)'
      WRITE (*,*) 'min abs relative error:',minval(tmp_DBSFC(:,:,:)),'(DBSFC)'
      WRITE (*,*) 'ave abs relative error:',absSFCerr/(nx_block*ny_block*km),'(DBSFC)'
   endif


 deallocate(aveLOC)
 deallocate(aveSFC)

 deallocate(tmp_DBLOC)
 deallocate(tmp_DBSFC)

 deallocate(MWORK11)
 deallocate(MWORK1K)
 deallocate(MWORK1KM)
 deallocate(MWORK21)
 deallocate(MWORK2K)
 deallocate(SWORK11)
 deallocate(SWORK1K)
 deallocate(SWORK1KM)
 deallocate(SWORK21)
 deallocate(SWORK2K)
#endif

#if defined(POPCHKJN) || defined(POPTIMEJN)
 num_steps = num_steps + 1
#endif


 deallocate(LOCKMT)    


!-----------------------------------------------------------------------
!EOC

 end subroutine buoydiff











!***********************************************************************
!BOP
! !IROUTINE: add_kpp_sources
! !INTERFACE:

 subroutine add_kpp_sources(SRCARRAY, k, this_block)

! !DESCRIPTION:
!  This routine adds KPP non local mixing term to the tracer source 
!  tendencies.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                   ! vertical level index

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      SRCARRAY                ! array of tracer sources

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! local dummy index for tracer
      bid                 ! local block address

!-----------------------------------------------------------------------
!
!  if KPP not chosen, return
!
!-----------------------------------------------------------------------

   if (.not. allocated(KPP_SRC)) return

!-----------------------------------------------------------------------
!
!  determine location and add sources
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   SRCARRAY = SRCARRAY + KPP_SRC(:,:,k,:,bid)

   do n=1,nt
      call accumulate_tavg_field(KPP_SRC(:,:,k,n,bid),tavg_KPP_SRC(n),bid,k)
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine add_kpp_sources






!***********************************************************************
!BOP
! !IROUTINE: smooth_hblt
! !INTERFACE:

 subroutine smooth_hblt (overwrite_hblt, use_hmxl, &
                         bid, HBLT, KBL, SMOOTH_OUT)

! !DESCRIPTION:
!  This subroutine uses a 1-1-4-1-1 Laplacian filter one time
!  on HBLT or HMXL to reduce any horizontal two-grid-point noise.
!  If HBLT is overwritten, KBL is adjusted after smoothing.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in) :: &
      overwrite_hblt,   &    ! if .true.,  HBLT is overwritten
                             ! if .false., the result is returned in
                             !  a dummy array
      use_hmxl               ! if .true., smooth HMXL
                             ! if .false., smooth HBLT

   integer (int_kind), intent(in) :: &
      bid                    ! local block address

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), optional, intent(inout) :: &
      HBLT                   ! boundary layer depth

   integer (int_kind), dimension(nx_block,ny_block), optional, intent(inout) :: &
      KBL                    ! index of first lvl below hbl

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), optional, intent(out) ::  &
      SMOOTH_OUT              ! optional output array containing the
                              !  smoothened field if overwrite_hblt is false

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   character (char_len) ::  &
      message

   integer (int_kind) :: &
      i, j,              &  ! horizontal loop indices
      k                     ! vertical level index

   real (r8), dimension(nx_block,ny_block) ::  &
      WORK1, WORK2

   real (r8) ::  &
     cc, cw, ce, cn, cs, &  ! averaging weights
     ztmp                   ! temp for level depth

!-----------------------------------------------------------------------
!
!     consistency checks 
!
!-----------------------------------------------------------------------

   if ( overwrite_hblt  .and.  ( .not.present(KBL)  .or.        &
                                 .not.present(HBLT) ) ) then      
     message = 'incorrect subroutine arguments for smooth_hblt, error # 1'
     call exit_POP (sigAbort, trim(message))
   endif

   if ( .not.overwrite_hblt  .and.  .not.present(SMOOTH_OUT) ) then 
     message = 'incorrect subroutine arguments for smooth_hblt, error # 2'
     call exit_POP (sigAbort, trim(message))
   endif

   if ( use_hmxl .and. .not.present(SMOOTH_OUT) ) then          
     message = 'incorrect subroutine arguments for smooth_hblt, error # 3'
     call exit_POP (sigAbort, trim(message))
   endif

   if ( overwrite_hblt  .and.  use_hmxl ) then                  
     message = 'incorrect subroutine arguments for smooth_hblt, error # 4'
     call exit_POP (sigAbort, trim(message))
   endif

!-----------------------------------------------------------------------
!
!     perform one smoothing pass since we cannot do the necessary 
!     boundary updates for multiple passes.
!
!-----------------------------------------------------------------------

   if ( use_hmxl ) then
     WORK2 = HMXL(:,:,bid)
   else
     WORK2 = HBLT
   endif

   WORK1 = WORK2
   do j=2,ny_block-1
     do i=2,nx_block-1
       if ( KMT(i,j,bid) /= 0 ) then
         cw = p125
         ce = p125
         cn = p125
         cs = p125
         cc = p5
         if ( KMT(i-1,j,bid) == 0 ) then
           cc = cc + cw
           cw = c0
         endif
         if ( KMT(i+1,j,bid) == 0 ) then
           cc = cc + ce
           ce = c0
         endif
         if ( KMT(i,j-1,bid) == 0 ) then
           cc = cc + cs
           cs = c0
         endif
         if ( KMT(i,j+1,bid) == 0 ) then
           cc = cc + cn
           cn = c0
         endif
         WORK2(i,j) =  cw * WORK1(i-1,j)   &
                     + ce * WORK1(i+1,j)   &
                     + cs * WORK1(i,j-1)   &
                     + cn * WORK1(i,j+1)   &
                     + cc * WORK1(i,j)
       endif
     enddo
   enddo

   do k=1,km
     do j=2,ny_block-1
       do i=2,nx_block-1

         if (partial_bottom_cells) then
           ztmp = -zgrid(k-1) + p5*(DZT(i,j,k-1,bid) + &
                                    DZT(i,j,k  ,bid))
         else
           ztmp = -zgrid(k)
         endif

         if ( k == KMT(i,j,bid)  .and.  WORK2(i,j) > ztmp ) then
           WORK2(i,j) = ztmp
         endif

       enddo
     enddo
   enddo

   if ( overwrite_hblt  .and.  .not.use_hmxl ) then

     HBLT = WORK2

     do k=1,km
       do j=2,ny_block-1
         do i=2,nx_block-1

           if (partial_bottom_cells) then
             ztmp = -zgrid(k-1) + p5*(DZT(i,j,k-1,bid) + &
                                      DZT(i,j,k  ,bid))
           else
             ztmp = -zgrid(k)
           endif

           if ( KMT(i,j,bid) /= 0            .and.  &
                ( HBLT(i,j) >  -zgrid(k-1) ) .and.  &
                ( HBLT(i,j) <= ztmp        ) ) KBL(i,j) = k
     
         enddo
       enddo
     enddo

   else

     SMOOTH_OUT = WORK2

   endif

!-----------------------------------------------------------------------

 end subroutine smooth_hblt





!***********************************************************************
!BOP
! !IROUTINE: compute_niw_energy_flux
! !INTERFACE:

 subroutine compute_niw_energy_flux (VISC,VDC,UUU,VVV,KE_mix,UCUR,VCUR,KE_cur,  &
                                     DBLOC, KPP_HBLT,KBL,En,this_block)

! !DESCRIPTION:
!  Compute niw energy flux
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UUU,VVV,    &! velocities at mix time
      UCUR,VCUR    ! velocities at current time

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      KPP_HBLT   

   type (block), intent(in) :: &
      this_block   ! block information for current block

! !INPUT/OUTPUT PARAMETERS:
   integer (int_kind), dimension(nx_block,ny_block), intent(inout) :: &
      KBL          ! index of first lvl below hbl

   real (r8), dimension(nx_block,ny_block,0:km+1), intent(inout) :: &
      VISC         ! viscosity

    real (r8), dimension(nx_block,ny_block,0:km+1,2),intent(inout) :: &
      VDC          ! diffusivity for tracer diffusion

    real (r8), dimension(nx_block,ny_block), intent(inout) :: &
      KE_mix,     &! kinetic energy at mix time
      KE_cur,     &! kinetic energy at cur time
      En           ! En for boundary layer kinetic energy

   real (r8), dimension(nx_block,ny_block,km), intent(inout) :: &
      DBLOC        ! buoyancy difference between adjacent levels


! !OUTPUT PARAMETERS:


!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   character (char_len) ::  &
      message

   integer (int_kind) :: &
      i,j,               &! horizontal indices
      k,                 &! vertical level index
      bid                 ! block id

   real (r8) ::  &
      factor              ! temporary scalar factor

     bid = this_block%local_id

!-----------------------------------------------------------------------
!  compute niw energy flux internally using blke method; blke means
!  "boundary layer kinetic energy".
!-----------------------------------------------------------------------

     KE_mix = c0

     if( niw_energy_type .eq. 'blke' ) then

!-----------------------------------------------------------------------
!  compute boundary layer kinetic energy for mix (mixing) time
!-----------------------------------------------------------------------

       call blke(UUU,VVV,KBL,KE_mix)

!-----------------------------------------------------------------------
!  compute boundary layer kinetic energy for cur (current) time
!-----------------------------------------------------------------------

       call blke(UCUR,VCUR,KBL,KE_cur)

!-----------------------------------------------------------------------
!  compute En for boundary layer KE, enforcing positivity. En is the
!  energy flux for generating near-inertial waves. Note, the factor
!  multiplying the kinetic energy (KE) tendency is an empirical scaling
!  factor to map the total time-step change in boundary layer kinetic 
!  energy into that fraction which generates near-inertial waves.
!-----------------------------------------------------------------------

       En(:,:) = 0.05_r8 * (KE_cur(:,:)-KE_mix(:,:))/dtt
       where( En < c0 )
         En(:,:) = -c1 * En(:,:)
       endwhere


!-----------------------------------------------------------------------
!  limit En to > 10 and < -10 degrees latitude. Near-inertial waves do
!  not generate mixing close to the equator.
!-----------------------------------------------------------------------

       do j=1,ny_block
         do i=1,nx_block

            if( TLATD(i,j,bid) > -5.0_r8 .and. TLATD(i,j,bid) < 5.0_r8 ) then
             En(i,j) = c0
            endif

            if( TLATD(i,j,bid) > -10.0_r8 .and. TLATD(i,j,bid) < 10.0_r8 ) then
              En(i,j) = En(i,j) * NIW_COS_FACTOR(i,j,bid)

            endif

         enddo
       enddo

!-----------------------------------------------------------------------
!  write En for boundary layer KE to history file, converting 
!  erg/cm2/sec to Watts/m2
!-----------------------------------------------------------------------

       call accumulate_tavg_field(En/c1000,tavg_En,bid,1)

!-----------------------------------------------------------------------
!  include other factors for niw parameterization in final definition 
!  of En. These other factors account for various fractions of En that 
!  result in diffusivity forming mixing in the column.
!-----------------------------------------------------------------------

       En(:,:) = NIW_COEF(:,:,bid) * En(:,:)

   else ! niw_energy_type

!-----------------------------------------------------------------------
!  write En for boundary layer KE to history file,
!  converting erg/cm2/sec to Watts/m2
!-----------------------------------------------------------------------

       En(:,:) = NIW_ENERGY_FLUX(:,:,bid)
       call accumulate_tavg_field(En/c1000,tavg_En,bid,1)

!-----------------------------------------------------------------------
!  use external pre-defined niw energy source flux, including other
!  factors for niw parameterization in final definition of En
!-----------------------------------------------------------------------

       En(:,:) = NIW_COEF(:,:,bid) * NIW_ENERGY_FLUX(:,:,bid)

     endif ! niw_energy_type

!-----------------------------------------------------------------------
!  write mix time KE to history file
!-----------------------------------------------------------------------

     call accumulate_tavg_field(KE_mix,tavg_KE_BL,bid,1)


!-----------------------------------------------------------------------
!
!  set background diffusivities and viscosities
!
!-----------------------------------------------------------------------

     call iw_reset(VISC, VDC, this_block)

!-----------------------------------------------------------------------
!
!  compute mixing due to near inertial waves and combine with original
!  background values.
!
!-----------------------------------------------------------------------

     call niw_mix(DBLOC, KPP_HBLT, KBL, En, VISC, VDC, this_block)

!-----------------------------------------------------------------------

 end subroutine compute_niw_energy_flux





!***********************************************************************
!BOP
! !IROUTINE: blke
! !INTERFACE:

 subroutine blke(UUU, VVV, KBL, KE)

! !DESCRIPTION:
!  Computes boundary layer kinetic energy (per unit area)
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UUU               ! U velocities at current time

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      VVV               ! V velocities at current time

   integer (int_kind), dimension(nx_block,ny_block), intent(in) :: &
      KBL               ! index of first lvl below hbl

   real (r8), dimension(nx_block,ny_block) :: &
      KE                  ! kinetic energy


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k                   ! index for vertical levels


!-----------------------------------------------------------------------
!
!  compute total ke above mixed layer
!
!-----------------------------------------------------------------------

   KE = c0
   do k = 1,km
     where ( k <= KBL )
       KE(:,:) = KE(:,:) + &
         0.5_r8 * rho_sw * (UUU(:,:,k)**2 + VVV(:,:,k)**2) * dz(k)
     endwhere
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine blke





!***********************************************************************
!BOP
! !IROUTINE: iw_reset
! !INTERFACE:

 subroutine iw_reset(VISC, VDC, this_block)

! !DESCRIPTION:
!  Initialize viscosity and diffusivity coefficients to background values.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,0:km+1), intent(inout) :: & 
      VISC       ! viscosity

   real (r8), dimension(nx_block,ny_block,0:km+1,2), intent(inout) :: & 
      VDC        ! diffusivity for tracer diffusion

   type (block), intent(in) :: &
      this_block          ! block information for current block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k,                 &! index for vertical levels
      i,j,               &! horizontal loop indices
      bid                 ! local block index

!-----------------------------------------------------------------------
!
!  reset mixing to background values
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

    do k=1,km

      VISC(:,:,k  ) = bckgrnd_vvc(:,:,k,bid)
      VDC (:,:,k,1) = bckgrnd_vdc(:,:,k,bid)
      VDC (:,:,k,2) = bckgrnd_vdc(:,:,k,bid)

      ! k index shifted because bckgrnd_vdc and bckgrnd_vvc are at cell bottom
      ! while output axis is at cell top
      call accumulate_tavg_field(bckgrnd_vdc(:,:,k,bid),tavg_VDC_BCK,bid,k)
      call accumulate_tavg_field(bckgrnd_vvc(:,:,k,bid),tavg_VVC_BCK,bid,k)

    end do ! k

!-----------------------------------------------------------------------
!
!  fill extra coefficients for blmix
!
!-----------------------------------------------------------------------

   VISC(:,:,0  ) = c0
   VDC (:,:,0,:) = c0
   VISC(:,:,km+1  ) = c0 
   VDC (:,:,km+1,:) = c0 

!-----------------------------------------------------------------------
!EOC
 
 end subroutine iw_reset

!***********************************************************************

 end module vmix_kpp

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
