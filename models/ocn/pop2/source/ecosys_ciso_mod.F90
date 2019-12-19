!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ecosys_ciso_mod

!BOP
! !MODULE: ecosys_ciso_mod
!
! !DESCRIPTION:
!
!  Carbon 13 module and biotic 14C module
!  13C code is based on code form G. Xavier, ETH, 2010, which
!  was written for pop1 (CCSM3)
!  This code needs the ecosystem model to run, as it uses several
!  variables computed there. Data is shared using ecosys_share_mod.F90
!  This module adds 7 carbon pools for 13C and another 7 for 14C
!
!  Developer: Alexandra Jahn, NCAR, Started Nov 2012, last edited July 2014
!  A first version of the 13C code for POP1, which was used to develop 
!  the code for the current model version, was written by Xavier Giraud 
!
! !REVISION HISTORY:

!  SVN:$Id:  $

!-----------------------------------------------------------------------
!  variables/subroutines/function used from other modules
!  The following are used extensively in this module, so they are used at
!  the module level. The use statements for variables that are only needed
!  locally are located at the module subprogram level.
!-----------------------------------------------------------------------

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_GridHorzMod
   use POP_FieldMod
   use POP_HaloMod

   use kinds_mod
   use constants
   use communicate
   use broadcast
   use global_reductions
   use blocks
   use domain_size
   use domain
   use exit_mod
   use prognostic
   use grid
   use io
   use io_types
   use io_tools
   use tavg
   use timers
   use passive_tracer_tools
   use named_field_mod
   use forcing_tools
   use time_management
   use ecosys_parms
   use registry
   use named_field_mod
   use ecosys_share
#ifdef CCSMCOUPLED
   use POP_MCT_vars_mod
   use shr_strdata_mod
#endif

! !INPUT PARAMETERS:
!-----------------------------------------------------------------------
!  include ecosystem parameters
!  all variables from this modules have a parm_ prefix
!-----------------------------------------------------------------------

   implicit none
   save
   private


!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
      ecosys_ciso_tracer_cnt,            &
      ecosys_ciso_init,                  &
      ecosys_ciso_tracer_ref_val,        &
      ecosys_ciso_set_sflux,             &
      ecosys_ciso_tavg_forcing,          &
      ecosys_ciso_set_interior,          &
      ecosys_ciso_write_restart

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by forcing_passive_tracer
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      ecosys_ciso_tracer_cnt = 14

!-----------------------------------------------------------------------
!  flags controlling which portion of code are executed
!  usefull for debugging
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      ciso_lsource_sink

   logical (log_kind), dimension(:,:,:), allocatable :: &
      LAND_MASK

!-----------------------------------------------------------------------
!  non-autotroph relative tracer indices
!  autotroph relative tracer indices are in autotroph derived type and
!  are determined at run time
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      di13c_ind          =  1,  & ! dissolved inorganic carbon 13
      do13c_ind          =  2,  & ! dissolved organic carbon 13
      zoo13C_ind         =  3,  & ! zooplankton carbon 13
      di14c_ind          =  4,  & ! dissolved inorganic carbon 14
      do14c_ind          =  5,  & ! dissolved organic carbon 14
      zoo14C_ind         =  6     ! zooplankton carbon 14

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(ecosys_ciso_tracer_cnt) :: &
      ciso_ind_name_table

!-----------------------------------------------------------------------
!  tavg ids and buffer indices (into ECO_CISO_SFLUX_TAVG) for 2d fields
!  related to surface fluxes. Suplicates, which are used for placing fields
!  into multiple tavg streams, do not need separate buffer indices.
!  fields that are recoverable from the STF field do not need separate
!  buffer indices
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_CISO_DI13C_GAS_FLUX,       buf_ind_CISO_FG_13CO2,     & ! di13c flux
      tavg_CISO_DI13C_AS_GAS_FLUX,    buf_ind_CISO_FG_as_13CO2,  & ! air-sea di13c flux
      tavg_CISO_DI13C_SA_GAS_FLUX,    buf_ind_CISO_FG_sa_13CO2,  & ! sea-air di13c flux
      tavg_CISO_d13C_GAS_FLUX,        buf_ind_CISO_FG_d13C,      & ! surface ocean delta 13C
      tavg_CISO_R13C_DIC_surf,        buf_ind_CISO_R13C_DIC_surf,& ! 13C/12C ratio in total DIC
      tavg_CISO_R13C_atm,             buf_ind_CISO_R13C_atm,     & ! atmospheric ratio of 13C/12C
      tavg_CISO_D13C_atm,             buf_ind_CISO_D13C_atm,     & ! atmospheric delta13C in permil
      tavg_CISO_DI14C_GAS_FLUX,       buf_ind_CISO_FG_14CO2,     & ! di14c flux
      tavg_CISO_DI14C_AS_GAS_FLUX,    buf_ind_CISO_FG_as_14CO2,  & ! air-sea di14c flux
      tavg_CISO_DI14C_SA_GAS_FLUX,    buf_ind_CISO_FG_sa_14CO2,  & ! sea-air di14c flux
      tavg_CISO_d14C_GAS_FLUX,        buf_ind_CISO_FG_d14C,      & ! surface ocean delta 14C
      tavg_CISO_R14C_DIC_surf,        buf_ind_CISO_R14C_DIC_surf,& ! 14C/12C ratio in total DIC
      tavg_CISO_R14C_atm,             buf_ind_CISO_R14C_atm,     & ! atmospheric ratio of 14C/12C
      tavg_CISO_D14C_atm,             buf_ind_CISO_D14C_atm,     & ! atmospheric delta14C in permil
      tavg_CISO_DI13C_RIV_FLUX,       buf_ind_DI13C_RIV_FLUX,    & ! river input of DI13C
      tavg_CISO_DO13C_RIV_FLUX,       buf_ind_DO13C_RIV_FLUX,    & ! river input of DO13C
      tavg_CISO_DI14C_RIV_FLUX,       buf_ind_DI14C_RIV_FLUX,    & ! river input of DI14C
      tavg_CISO_DO14C_RIV_FLUX,       buf_ind_DO14C_RIV_FLUX,    & ! river input of DO14C
      tavg_CISO_eps_aq_g_surf,        buf_ind_CISO_eps_aq_g_surf,& ! tavg id for eps_aq_g_surf
      tavg_CISO_eps_dic_g_surf,       buf_ind_CISO_eps_dic_g_surf  ! tavg id for eps_dic_g_surf

! for debugging
   integer (int_kind) ::       &
      tavg_CISO_GLOBAL_D14C,   & ! tavg id for the global averaged atmos. D14C
      buf_ind_CISO_GLOBAL_D14C   ! buffer index for the global averaged atmos. D14C

!-----------------------------------------------------------------------
!  define tavg id for nonstandard 3d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_CISO_PO13C_FLUX_IN,    &! tavg id for po13c flux into cell
      tavg_CISO_PO13C_PROD,       &! tavg id for po13c production
      tavg_CISO_PO13C_REMIN,      &! tavg id for po13c remineralization
      tavg_CISO_Ca13CO3_FLUX_IN,  &! tavg id for ca13co3 flux into cell
      tavg_CISO_Ca13CO3_PROD,     &! tavg id for ca13co3 production
      tavg_CISO_Ca13CO3_REMIN,    &! tavg id for ca13co3 remineralization
      tavg_CISO_PO14C_FLUX_IN,    &! tavg id for po14c flux into cell
      tavg_CISO_PO14C_PROD,       &! tavg id for po14c production
      tavg_CISO_PO14C_REMIN,      &! tavg id for po14c remineralization
      tavg_CISO_Ca14CO3_FLUX_IN,  &! tavg id for ca14co3 flux into cell
      tavg_CISO_Ca14CO3_PROD,     &! tavg id for ca14co3 production
      tavg_CISO_Ca14CO3_REMIN,    &! tavg id for ca14co3 remineralization
      tavg_CISO_photo13C_TOT,     &! tavg id for total 13C fixation
      tavg_CISO_photo13C_TOT_zint,&! tavg id for total 13C fixation vertical integral
      tavg_CISO_photo14C_TOT,     &! tavg id for total 14C fixation
      tavg_CISO_photo14C_TOT_zint  ! tavg id for total 14C fixation vertical integral

!-----------------------------------------------------------------------
!  define tavg id for MORE nonstandard 3d fields
!-----------------------------------------------------------------------
   integer (int_kind) :: &
      tavg_CISO_eps_aq_g,        & ! tavg id for eps_aq_g
      tavg_CISO_eps_dic_g          ! tavg id for eps_dic_g

   integer (int_kind), dimension(autotroph_cnt) :: &
      tavg_CISO_eps_autotroph,   & ! tavg id for epsilon for each autotroph
      tavg_CISO_mui_to_co2star     ! tavg id for mui_to_co2star for each autotroph

!-----------------------------------------------------------------------
!  define tavg id for MORE nonstandard 3d fields
!-----------------------------------------------------------------------

   integer (int_kind), dimension(autotroph_cnt) :: &
      tavg_CISO_Ca13CO3_form,       &! tavg id for Ca13CO3 formation
      tavg_CISO_Ca14CO3_form,       &! tavg id for Ca14CO3 formation
      tavg_CISO_Ca13CO3_form_zint,  &! tavg id for Ca13CO3 formation vertical integral 0-100 m
      tavg_CISO_Ca14CO3_form_zint,  &! tavg id for Ca14CO3 formation vertical integral 0-100 m
      tavg_CISO_photo13C,           &! tavg id for 13C fixation
      tavg_CISO_photo13C_zint,      &! tavg id for 13C fixation vertical integral
      tavg_CISO_photo14C,           &! tavg id for 14C fixation
      tavg_CISO_photo14C_zint,      &! tavg id for 14C fixation vertical integral
      tavg_CISO_d13C,               &! tavg if for d13C of autotroph carbon
      tavg_CISO_d14C,               &! tavg if for d14C of autotroph carbon
      tavg_CISO_autotrophCaCO3_d14C,&! tavg if for d14C of autotrophCaCO3
      tavg_CISO_autotrophCaCO3_d13C  ! tavg if for d13C of autotrophCaCO3


   integer (int_kind) :: &
      tavg_CISO_DO13C_prod,         &! tavg id for do13c production
      tavg_CISO_DO13C_remin,        &! tavg id for do13c remineralization
      tavg_CISO_Jint_13Ctot,        &! tavg id for vertically integrated source sink term, 13Ctot
      tavg_CISO_Jint_100m_13Ctot,   &! tavg id for vertically integrated source sink term, 0-100m, 13Ctot
      tavg_CISO_DO14C_prod,         &! tavg id for do14c production
      tavg_CISO_DO14C_remin,        &! tavg id for do14c remineralization
      tavg_CISO_Jint_14Ctot,        &! tavg id for vertically integrated source sink term, 14Ctot
      tavg_CISO_Jint_100m_14Ctot,   &! tavg id for vertically integrated source sink term, 0-100m, 14Ctot
      tavg_CISO_zooC_d14C,          &! tavg if for d14C of zooC
      tavg_CISO_DOC_d14C,           &! tavg if for d14C of DOC
      tavg_CISO_DIC_d14C,           &! tavg if for d14C of DIC
      tavg_CISO_zooC_d13C,          &! tavg if for d13C of zooC
      tavg_CISO_DOC_d13C,           &! tavg if for d13C of DOC
      tavg_CISO_DIC_d13C             ! tavg if for d13C of DIC

   integer (int_kind) :: &
      tavg_calcToSed_13C,      &! tavg id for calcite flux sedimentary burial
      tavg_pocToSed_13C,       &! tavg id for poc burial flux to sediments
      tavg_calcToSed_14C,      &! tavg id for calcite flux sedimentary burial
      tavg_pocToSed_14C         ! tavg id for poc burial flux to sediments

!-----------------------------------------------------------------------
!  define array for holding flux-related quantities that need to be time-averaged
!  this is necessary since the forcing routines are called before tavg flags
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable :: &
      ECO_CISO_SFLUX_TAVG

!-----------------------------------------------------------------------
!  ciso_data_ind_d13c is the index for the D13C data for the
!  current timestep
!  Note that ciso_data_ind_d13c is always less than ciso_atm_d13c_data_nbval.
!  To enable OpenMP parallelism, duplicating data_ind for each block
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:), allocatable :: &
      ciso_data_ind_d13c, &      ! data index for D13C data
      ciso_data_ind_d14c         ! data index for D14C data

!-----------------------------------------------------------------------
!  average surface tracer value related variables
!  used as reference value for virtual flux computations
!-----------------------------------------------------------------------

   logical (log_kind), dimension(ecosys_ciso_tracer_cnt) :: &
      ciso_vflux_flag                ! which tracers get virtual fluxes applied

   integer (int_kind) :: &
      ciso_comp_surf_avg_flag        ! time flag id for computing average
                                     ! surface tracer values

   real (r8), dimension(ecosys_ciso_tracer_cnt) :: &
      ciso_surf_avg                  ! average surface tracer values

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ecosys_ciso_interior_timer,                   &
      ecosys_ciso_sflux_timer

!-----------------------------------------------------------------------
!  forcing related variables
!-----------------------------------------------------------------------

   integer (int_kind) ::          &
      ciso_atm_model_year,        &  ! arbitrary model year
      ciso_atm_data_year,         &  ! year in atmospheric ciso data that corresponds to ciso_atm_model_year
      ciso_atm_d13c_data_nbval,   &  ! number of values in ciso_atm_d13c_filename
      ciso_atm_d14c_data_nbval       ! number of values in ciso_atm_d14c_filename

   real (r8), dimension(:) , allocatable :: &
      ciso_atm_d13c_data,    &    !  atmospheric D13C values in datafile
      ciso_atm_d13c_data_yr       !  date of atmospheric D13C values in datafile

   real (r8), dimension(:,:) , allocatable :: &
      ciso_atm_d14c_data,    &    !  atmospheric D14C values in datafile (sh, eq, nh, in permil)
      ciso_atm_d14c_data_yr       !  date of atmospheric D14C values in datafile (sh, eq, nh)

   real (r8) :: &
      ciso_atm_d13c_const, &      !  atmospheric D13C constant [permil]
      ciso_atm_d14c_const         !  atmospheric D14C constant [permil]

   character(char_len) ::    &
      ciso_atm_d13c_opt,     &    ! option for CO2 and D13C varying or constant forcing
      ciso_atm_d13c_filename,&    ! filenames for varying atm D13C
      ciso_atm_d14c_opt           ! option for CO2 and D13C varying or constant forcing

   character (char_len), dimension(3) :: &
      ciso_atm_d14c_filename      ! filenames for varying atm D14C (one each for NH, SH, EQ)

!-----------------------------------------------------------------------
!  fractionation related variable
!-----------------------------------------------------------------------
   character(char_len) ::    &
      ciso_fract_factors          ! option for which biological fractionation calculation to use

!-----------------------------------------------------------------------
!  scalar constants for 14C decay calculation
!-----------------------------------------------------------------------

   real (r8), parameter :: c14_halflife_years = 5730.0_r8 !C14 half file
   real (r8) :: c14_lambda_inv_sec           ! Decay variable in seconds

!---------------------------------------------------------------------
!     Isotope standards
!---------------------------------------------------------------------

! Using scaled isotopic carbon pools, so Rstd =1
   real(r8), parameter :: &
      R13C_std = 1.0_r8,    & ! actual 13C/12C PDB standard ratio (Craig, 1957) = 1123.72e-5_r8
      R14C_std = 1.0_r8       ! actual 14C/12C NOSAMS standard ratio = 11.76e-13_r8

!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: ecosys_ciso_init
! !INTERFACE:

 subroutine ecosys_ciso_init(init_ts_file_fmt, read_restart_filename, &
                        tracer_d_module, TRACER_MODULE, lmarginal_seas, &
                        errorCode)

! !DESCRIPTION:
!  Initialize ecosys_ciso tracer module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

   logical (kind=log_kind), intent(in) :: &
      lmarginal_seas               

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(:), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(:,:,:,:,:,:), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'ciso_mod:ecosys_ciso_init'

   character(char_len) :: &
      ciso_init_ecosys_option,        & ! option for initialization of bgc
      ciso_init_ecosys_init_file,     & ! filename for option 'file'
      ciso_init_ecosys_init_file_fmt, & ! file format for option 'file'
      ciso_comp_surf_avg_freq_opt

   type(tracer_read), dimension(ecosys_ciso_tracer_cnt) :: &
      ciso_tracer_init_ext              ! namelist variable for initializing tracers


   logical (log_kind) :: &
      default                      ! arg to init_time_flag

   integer (int_kind) :: &
      non_autotroph_ecosys_ciso_tracer_cnt, & ! number of non-autotroph ecosystem tracers
      auto_ind,                        & ! autotroph functional group index
      n,                               & ! index for looping over tracers
      k,                               & ! index for looping over depth levels
      l,                               & ! index for looping over time levels
      ind,                             & ! tracer index for tracer name from namelist
      iblock,                          & ! index for looping over blocks
      nml_error                          ! namelist i/o error flag

   integer (int_kind) :: &
      freq_opt, freq,              & ! args for init_time_flag
      ciso_comp_surf_avg_freq_iopt,& ! choice for freq of comp_surf_avg
      ciso_comp_surf_avg_freq        ! choice for freq of comp_surf_avg

   logical (log_kind) :: &
      ciso_use_nml_surf_vals         ! do namelist surf values override values from restart file

   logical (log_kind) :: &
      ciso_lecovars_full_depth_tavg  ! should ecosystem vars be written full depth

!-----------------------------------------------------------------------
!  values to be used when comp_surf_avg_freq_opt==never
!-----------------------------------------------------------------------

   real (r8) :: &
      ciso_surf_avg_di13c_const, &
      ciso_surf_avg_di14c_const
!-----------------------------------------------------------------------
! ecosys_ciso_nml namelist
!-----------------------------------------------------------------------

   namelist /ecosys_ciso_nml/ &
      ciso_init_ecosys_option, ciso_init_ecosys_init_file, &
      ciso_init_ecosys_init_file_fmt, ciso_tracer_init_ext, &
      ciso_comp_surf_avg_freq_opt, ciso_comp_surf_avg_freq,  &
      ciso_use_nml_surf_vals, ciso_surf_avg_di13c_const, &
      ciso_surf_avg_di14c_const, &
      ciso_lsource_sink, &
      ciso_lecovars_full_depth_tavg, &
      ciso_atm_d13c_opt, ciso_atm_d13c_const, ciso_atm_d13c_filename, &
      ciso_atm_d14c_opt, ciso_atm_d14c_const, ciso_atm_d14c_filename, &
      ciso_fract_factors, ciso_atm_model_year, ciso_atm_data_year

   character (char_len) :: &
      ecosys_ciso_restart_filename  ! modified file name for restart file


!-----------------------------------------------------------------------
!  initialize name table
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!  initialize non-autotroph tracer_d values
!  accumulate non_autotroph_ecosys_ciso_tracer_cnt
!-----------------------------------------------------------------------
   non_autotroph_ecosys_ciso_tracer_cnt = 0

   tracer_d_module(di13c_ind)%short_name='DI13C'
   tracer_d_module(di13c_ind)%long_name='Dissolved Inorganic Carbon-13'
   non_autotroph_ecosys_ciso_tracer_cnt = non_autotroph_ecosys_ciso_tracer_cnt + 1

   tracer_d_module(do13c_ind)%short_name='DO13C'
   tracer_d_module(do13c_ind)%long_name='Dissolved Organic Carbon-13'
   non_autotroph_ecosys_ciso_tracer_cnt = non_autotroph_ecosys_ciso_tracer_cnt + 1

   tracer_d_module(zoo13C_ind)%short_name='zoo13C'
   tracer_d_module(zoo13C_ind)%long_name='Zooplankton Carbon-13'
   non_autotroph_ecosys_ciso_tracer_cnt = non_autotroph_ecosys_ciso_tracer_cnt + 1

   tracer_d_module(di14c_ind)%short_name='DI14C'
   tracer_d_module(di14c_ind)%long_name='Dissolved Inorganic Carbon-14'
   non_autotroph_ecosys_ciso_tracer_cnt = non_autotroph_ecosys_ciso_tracer_cnt + 1

   tracer_d_module(do14c_ind)%short_name='DO14C'
   tracer_d_module(do14c_ind)%long_name='Dissolved Organic Carbon-14'
   non_autotroph_ecosys_ciso_tracer_cnt = non_autotroph_ecosys_ciso_tracer_cnt + 1

   tracer_d_module(zoo14C_ind)%short_name='zoo14C'
   tracer_d_module(zoo14C_ind)%long_name='Zooplankton Carbon-14'
   non_autotroph_ecosys_ciso_tracer_cnt = non_autotroph_ecosys_ciso_tracer_cnt + 1

   do n = 1, non_autotroph_ecosys_ciso_tracer_cnt
      tracer_d_module(n)%units      = 'mmol/m^3'
      tracer_d_module(n)%tend_units = 'mmol/m^3/s'
      tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
   end do

!-----------------------------------------------------------------------
!  confirm that ecosys_ciso_tracer_cnt is consistent with autotroph declarations
!-----------------------------------------------------------------------

   n = non_autotroph_ecosys_ciso_tracer_cnt
   do auto_ind = 1, autotroph_cnt
      n = n + 2 ! C13, C14 tracers
      if (autotrophs(auto_ind)%imp_calcifier .or. &
          autotrophs(auto_ind)%exp_calcifier) n = n + 2 ! Ca13CO3 & Ca14CO3 tracers
   end do

   if (ecosys_ciso_tracer_cnt /= n) then
      call document(subname, 'actual ecosys_ciso_tracer_cnt', ecosys_ciso_tracer_cnt)
      call document(subname, 'computed ecosys_ciso_tracer_cnt', n)
      call exit_POP(sigAbort, 'inconsistency between actual ecosys_ciso_tracer_cnt and computed ecosys_ciso_tracer_cnt')
   endif

!-----------------------------------------------------------------------
!  initialize autotroph tracer_d values and tracer indices
!-----------------------------------------------------------------------

   n = non_autotroph_ecosys_ciso_tracer_cnt + 1

   do auto_ind = 1, autotroph_cnt
      tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // '13C'
      tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Carbon-13'
      tracer_d_module(n)%units      = 'mmol/m^3'
      tracer_d_module(n)%tend_units = 'mmol/m^3/s'
      tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
      autotrophs(auto_ind)%C13_ind = n
      n = n + 1


      tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // '14C'
      tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Carbon-14'
      tracer_d_module(n)%units      = 'mmol/m^3'
      tracer_d_module(n)%tend_units = 'mmol/m^3/s'
      tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
      autotrophs(auto_ind)%C14_ind = n
      n = n + 1

      if (autotrophs(auto_ind)%imp_calcifier .or. &
         autotrophs(auto_ind)%exp_calcifier) then
         tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Ca13CO3'
         tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Ca13CO3'
         tracer_d_module(n)%units      = 'mmol/m^3'
         tracer_d_module(n)%tend_units = 'mmol/m^3/s'
         tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
         autotrophs(auto_ind)%Ca13CO3_ind = n
         n = n + 1
         tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Ca14CO3'
         tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Ca14CO3'
         tracer_d_module(n)%units      = 'mmol/m^3'
         tracer_d_module(n)%tend_units = 'mmol/m^3/s'
         tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
         autotrophs(auto_ind)%Ca14CO3_ind = n
         n = n + 1
      else
         autotrophs(auto_ind)%Ca13CO3_ind = 0
         autotrophs(auto_ind)%Ca14CO3_ind = 0
      endif
   end do

   if (my_task == master_task) THEN
      write (stdout,*) '----- autotroph tracer indices -----'
      do auto_ind = 1, autotroph_cnt
         write (stdout,*) 'C13_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%C13_ind
         write (stdout,*) 'C14_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%C14_ind
         write (stdout,*) 'Ca13CO3_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Ca13CO3_ind
         write (stdout,*) 'Ca14CO3_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Ca14CO3_ind
         write (stdout,*) 'autotroph_cnt =',autotroph_cnt
      end do
      write (stdout,*) '------------------------------------'
   endif


!-----------------------------------------------------------------------
!  initialize ind_name_table
!-----------------------------------------------------------------------

   do n = 1, ecosys_ciso_tracer_cnt
      ciso_ind_name_table(n) = ind_name_pair(n, tracer_d_module(n)%short_name)
   end do

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   ciso_init_ecosys_option                 = 'unknown'
   ciso_init_ecosys_init_file              = 'unknown'
   ciso_init_ecosys_init_file_fmt          = 'bin'
   do n = 1,ecosys_ciso_tracer_cnt
      ciso_tracer_init_ext(n)%mod_varname  = 'unknown'
      ciso_tracer_init_ext(n)%filename     = 'unknown'
      ciso_tracer_init_ext(n)%file_varname = 'unknown'
      ciso_tracer_init_ext(n)%scale_factor = c1
      ciso_tracer_init_ext(n)%default_val  = c0
      ciso_tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   ciso_lsource_sink                       = .true.

   ciso_comp_surf_avg_freq_opt             = 'never'
   ciso_comp_surf_avg_freq                 = 1
   ciso_use_nml_surf_vals                  = .false.
   ciso_surf_avg_di13c_const               = 1944.0_r8
   ciso_surf_avg_di14c_const               = 1944.0_r8

   ciso_atm_d13c_opt                       = 'const'
   ciso_atm_d13c_const                     = -6.379_r8
   ciso_atm_d13c_filename                  = 'unknown'

   ciso_atm_d14c_opt                       = 'const'
   ciso_atm_d14c_const                     = 0.0_r8
   ciso_atm_d14c_filename(1)               = 'unknown'
   ciso_atm_d14c_filename(2)               = 'unknown'
   ciso_atm_d14c_filename(3)               = 'unknown'

   ciso_fract_factors                      = 'Rau'

   ciso_lecovars_full_depth_tavg           = .false.

   ciso_atm_model_year                     = 1
   ciso_atm_data_year                      = 1

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=ecosys_ciso_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'ecosys_ciso_nml not found')
      call exit_POP(sigAbort, 'ERROR : stopping in '/&
                           &/ subname)
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' ecosys_ciso:'
      write(stdout,blank_fmt)
      write(stdout,*) ' ecosys_ciso_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,ecosys_ciso_nml)
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(ciso_init_ecosys_option, master_task)
   call broadcast_scalar(ciso_init_ecosys_init_file, master_task)
   call broadcast_scalar(ciso_init_ecosys_init_file_fmt, master_task)

   call broadcast_scalar(ciso_atm_d13c_opt, master_task)
   call broadcast_scalar(ciso_atm_d13c_const, master_task)
   call broadcast_scalar(ciso_atm_d13c_filename, master_task)

   call broadcast_scalar(ciso_atm_d14c_opt, master_task)
   call broadcast_scalar(ciso_atm_d14c_const, master_task)
   call broadcast_scalar(ciso_atm_d14c_filename(1), master_task)
   call broadcast_scalar(ciso_atm_d14c_filename(2), master_task)
   call broadcast_scalar(ciso_atm_d14c_filename(3), master_task)

   do n = 1,ecosys_ciso_tracer_cnt
      call broadcast_scalar(ciso_tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(ciso_tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(ciso_tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(ciso_tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(ciso_tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(ciso_tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(ciso_comp_surf_avg_freq_opt, master_task)
   call broadcast_scalar(ciso_comp_surf_avg_freq, master_task)
   call broadcast_scalar(ciso_use_nml_surf_vals, master_task)
   call broadcast_scalar(ciso_surf_avg_di13c_const, master_task)
   call broadcast_scalar(ciso_surf_avg_di14c_const, master_task)


   call broadcast_scalar(ciso_lsource_sink, master_task)

   call broadcast_scalar(ciso_lecovars_full_depth_tavg, master_task)

   call broadcast_scalar(ciso_fract_factors, master_task)

   call broadcast_scalar(ciso_atm_model_year, master_task)
   call broadcast_scalar(ciso_atm_data_year, master_task)

!-----------------------------------------------------------------------
!  set variables immediately dependent on namelist variables
!-----------------------------------------------------------------------

   select case (ciso_comp_surf_avg_freq_opt)
   case ('never')
      ciso_comp_surf_avg_freq_iopt = freq_opt_never
   case ('nyear')
      ciso_comp_surf_avg_freq_iopt = freq_opt_nyear
   case ('nmonth')
      ciso_comp_surf_avg_freq_iopt = freq_opt_nmonth
   case default
      call document(subname, 'ciso_comp_surf_avg_freq_opt', ciso_comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'unknown ciso_comp_surf_avg_freq_opt')
   end select

   call init_time_flag('ciso_ecosys_comp_surf_avg', ciso_comp_surf_avg_flag, &
     default=.false., freq_opt=ciso_comp_surf_avg_freq_iopt,  &
     freq=ciso_comp_surf_avg_freq, owner='ciso_ecosys_init')


!-----------------------------------------------------------------------
!  namelist consistency checking
!-----------------------------------------------------------------------

   if (ciso_use_nml_surf_vals .and. ciso_comp_surf_avg_freq_iopt /= freq_opt_never) then
      call document(subname, 'ciso_use_nml_surf_vals', ciso_use_nml_surf_vals)
      call document(subname, 'ciso_comp_surf_avg_freq_opt', ciso_comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'ciso_use_nml_surf_vals can only be .true. if ' /&
                           &/ ' ciso_comp_surf_avg_freq_opt is never')
   endif

!-----------------------------------------------------------------------
!  initialize virtual flux flag array
!-----------------------------------------------------------------------

   ciso_vflux_flag = .false.
   ciso_vflux_flag(di13c_ind) = .true.
   ciso_vflux_flag(di14c_ind) = .true.

!-----------------------------------------------------------------------
!  allocate and initialize LAND_MASK
!-----------------------------------------------------------------------

   allocate( LAND_MASK(nx_block,ny_block,nblocks_clinic) )

   if (lmarginal_seas) then
      LAND_MASK = REGION_MASK /= c0
   else
      LAND_MASK = REGION_MASK > c0
   endif

!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   select case (ciso_init_ecosys_option)

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid')

      ecosys_ciso_restart_filename = char_blank

      if (ciso_init_ecosys_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read ciso vars from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         ecosys_ciso_restart_filename = read_restart_filename
         ciso_init_ecosys_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         ecosys_ciso_restart_filename = trim(ciso_init_ecosys_init_file)

      endif

      call rest_read_tracer_block(ciso_init_ecosys_init_file_fmt, &
                                  ecosys_ciso_restart_filename,   &
                                  tracer_d_module,           &
                                  TRACER_MODULE)

      if (ciso_use_nml_surf_vals) then

         ciso_surf_avg = c0
         ciso_surf_avg(di13c_ind) = ciso_surf_avg_di13c_const
         ciso_surf_avg(di14c_ind) = ciso_surf_avg_di14c_const
      else
         call extract_surf_avg(ciso_init_ecosys_init_file_fmt, &
                               ecosys_ciso_restart_filename,   &
                               ecosys_ciso_tracer_cnt,         &
                               ciso_vflux_flag,                &
                               ciso_ind_name_table,ciso_surf_avg)
      endif

      call eval_time_flag(ciso_comp_surf_avg_flag) ! evaluates time_flag(ciso_comp_surf_avg_flag)%value via time_to_do

      if (check_time_flag(ciso_comp_surf_avg_flag)) &
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                                 TRACER_MODULE(:,:,1,:,curtime,:), &
                                 ecosys_ciso_tracer_cnt,     &
                                 ciso_vflux_flag,ciso_surf_avg)

   case ('file', 'ccsm_startup')
      call document(subname, 'ciso vars being read from separate files')

      call file_read_tracer_block(ciso_init_ecosys_init_file_fmt, &
                                  ciso_init_ecosys_init_file,     &
                                  tracer_d_module,           &
                                  ciso_ind_name_table,            &
                                  ciso_tracer_init_ext,           &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, ecosys_ciso_tracer_cnt
            do k=1,km
               call fill_points(k,TRACER_MODULE(:,:,k,n,oldtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'ecosys_init: error in fill points for tracers(oldtime)')
                  return
               endif

               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'ecosys_init: error in fill points for tracers(newtime)')
                  return
               endif

            enddo
         enddo
      endif

      if (ciso_use_nml_surf_vals) then
         ciso_surf_avg = c0
         ciso_surf_avg(di13c_ind) = ciso_surf_avg_di13c_const
         ciso_surf_avg(di14c_ind) = ciso_surf_avg_di14c_const
      else
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                            TRACER_MODULE(:,:,1,:,curtime,:), &
                            ecosys_ciso_tracer_cnt,           &
                            ciso_vflux_flag,ciso_surf_avg)
      endif

   case default
      call document(subname, 'ciso_init_ecosys_option', ciso_init_ecosys_option)
      call exit_POP(sigAbort, 'unknown ciso_init_ecosys_option')

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,n,k)
   do iblock=1,nblocks_clinic
      do n = 1,ecosys_ciso_tracer_cnt
         do k = 1,km
            where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            endwhere
         end do
      end do
   enddo
   !$OMP END PARALLEL DO


!-----------------------------------------------------------------------
!  timer init
!-----------------------------------------------------------------------

   call get_timer(ecosys_ciso_interior_timer, 'ECOSYS_CISO_INTERIOR', &
                  nblocks_clinic, distrb_clinic%nprocs)
   call get_timer(ecosys_ciso_sflux_timer, 'ECOSYS_CISO_SFLUX',1, &
                  distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  Define decay variable for DI14C, using earlier defined half-life of 14C
!-----------------------------------------------------------------------

   c14_lambda_inv_sec = log(c2) / (c14_halflife_years * seconds_in_year)

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call ecosys_ciso_init_tavg
   call ecosys_ciso_init_sflux

!-----------------------------------------------------------------------
!  set lfull_depth_tavg flag for short-lived ecosystem tracers
!-----------------------------------------------------------------------

   tracer_d_module(zoo13C_ind   )%lfull_depth_tavg = ciso_lecovars_full_depth_tavg
   tracer_d_module(zoo14C_ind   )%lfull_depth_tavg = ciso_lecovars_full_depth_tavg

   do auto_ind = 1, autotroph_cnt
      n = autotrophs(auto_ind)%C13_ind
      tracer_d_module(n)%lfull_depth_tavg = ciso_lecovars_full_depth_tavg

      n = autotrophs(auto_ind)%C14_ind
      tracer_d_module(n)%lfull_depth_tavg = ciso_lecovars_full_depth_tavg

      n = autotrophs(auto_ind)%Ca13CO3_ind
      if (n > 0) then
         tracer_d_module(n)%lfull_depth_tavg = ciso_lecovars_full_depth_tavg
      endif
      n = autotrophs(auto_ind)%Ca14CO3_ind
      if (n > 0) then
         tracer_d_module(n)%lfull_depth_tavg = ciso_lecovars_full_depth_tavg
      endif
   end do


!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_ciso_init

!***********************************************************************
!BOP
! !IROUTINE: ecosys_ciso_init_tavg
! !INTERFACE:

 subroutine ecosys_ciso_init_tavg

! !DESCRIPTION:
!  call define_tavg_field for nonstandard tavg fields
!
! !REVISION HISTORY:
!  same as module
!

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'ecosys_ciso_mod:ecosys_ciso_init_tavg'

   integer (int_kind) :: &
      auto_ind,       & ! autotroph functional group index
      buf_len           ! how many surface flux fields are stored in ECO_CISO_SFLUX_TAVG

   character(char_len) :: &
      sname             ! short-name of tavg variable


!-----------------------------------------------------------------------
!  2D fields related to surface fluxes
!-----------------------------------------------------------------------

   buf_len = 0


!13C
   call define_tavg_field(tavg_CISO_DI13C_GAS_FLUX,'CISO_FG_13CO2',2,            &
                          long_name='DI13C Surface Gas Flux',            &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_FG_13CO2 = buf_len

   call define_tavg_field(tavg_CISO_DI13C_AS_GAS_FLUX,'CISO_FG_as_13CO2',2,            &
                          long_name='DI13C Surface Air-Sea Gas Flux',            &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_FG_as_13CO2 = buf_len

   call define_tavg_field(tavg_CISO_DI13C_SA_GAS_FLUX,'CISO_FG_sa_13CO2',2,            &
                          long_name='DI13C Surface Sea-Air Gas Flux',            &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_FG_sa_13CO2 = buf_len

   call define_tavg_field(tavg_CISO_d13C_GAS_FLUX,'CISO_FG_d13C',2,                 &
                          long_name='D13C Surface GAS FLUX',&
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_FG_d13C = buf_len

   call define_tavg_field(tavg_CISO_D13C_atm,'CISO_D13C_atm',2,                 &
                          long_name='Atmospheric Delta 13C in permil',&
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_D13C_atm = buf_len

   call define_tavg_field(tavg_CISO_R13C_DIC_surf,'CISO_R13C_DIC_surf',2,                 &
                          long_name='13C/12C ratio in total DIC',&
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_R13C_DIC_surf = buf_len

   call define_tavg_field(tavg_CISO_R13C_atm,'CISO_R13C_atm',2,                 &
                          long_name='13C/12C ratio in atmosphere',&
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_R13C_atm = buf_len

   call define_tavg_field(tavg_CISO_DI13C_RIV_FLUX,'CISO_DI13C_RIV_FLUX',2,          &
                          long_name='Flux of DI13C from rivers',         &
                          units='nmol/cm^2/s', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_DI13C_RIV_FLUX = buf_len

   call define_tavg_field(tavg_CISO_DO13C_RIV_FLUX,'CISO_DO13C_RIV_FLUX',2,          &
                          long_name='Flux of DO13C from rivers',         &
                          units='nmol/cm^2/s', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_DO13C_RIV_FLUX = buf_len

! Fractionation (for 14C and 13C)

   call define_tavg_field(tavg_CISO_eps_aq_g_surf,'CISO_eps_aq_g_surf',2,              &
                          long_name='Surface equilibrium fractionation (CO2_gaseous <-> CO2_aq)', &
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')

   buf_len = buf_len+1
   buf_ind_CISO_eps_aq_g_surf = buf_len

   call define_tavg_field(tavg_CISO_eps_dic_g_surf,'CISO_eps_dic_g_surf',2,              &
                          long_name='Surface equilibrium fractionation between total DIC and gaseous CO2', &
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')

   buf_len = buf_len+1
   buf_ind_CISO_eps_dic_g_surf = buf_len

!14C
   call define_tavg_field(tavg_CISO_DI14C_GAS_FLUX,'CISO_FG_14CO2',2,            &
                          long_name='DI14C Surface Gas Flux',            &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_FG_14CO2 = buf_len

   call define_tavg_field(tavg_CISO_DI14C_AS_GAS_FLUX,'CISO_FG_as_14CO2',2,            &
                          long_name='DI14C Surface Air-Sea Gas Flux',            &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_FG_as_14CO2 = buf_len

   call define_tavg_field(tavg_CISO_DI14C_SA_GAS_FLUX,'CISO_FG_sa_14CO2',2,            &
                          long_name='DI14C Surface Sea-Air Gas Flux',            &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_FG_sa_14CO2 = buf_len

   call define_tavg_field(tavg_CISO_d14C_GAS_FLUX,'CISO_FG_d14C',2,                 &
                          long_name='D14C Surface GAS FLUX',&
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_FG_d14C = buf_len

   call define_tavg_field(tavg_CISO_D14C_atm,'CISO_D14C_atm',2,                 &
                          long_name='Atmospheric Delta 14C in permil',&
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_D14C_atm = buf_len

   call define_tavg_field(tavg_CISO_R14C_DIC_surf,'CISO_R14C_DIC_surf',2,                 &
                          long_name='14C/12C ratio in total DIC',&
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_R14C_DIC_surf = buf_len

   call define_tavg_field(tavg_CISO_R14C_atm,'CISO_R14C_atm',2,                 &
                          long_name='14C/12C ratio in atmosphere',&
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_R14C_atm = buf_len

   call define_tavg_field(tavg_CISO_DI14C_RIV_FLUX,'CISO_DI14C_RIV_FLUX',2,          &
                          long_name='Flux of DI14C from rivers',         &
                          units='nmol/cm^2/s', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_DI14C_RIV_FLUX = buf_len

   call define_tavg_field(tavg_CISO_DO14C_RIV_FLUX,'CISO_DO14C_RIV_FLUX',2,          &
                          long_name='Flux of DO14C from rivers',         &
                          units='nmol/cm^2/s', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_DO14C_RIV_FLUX = buf_len


! for debugging

   call define_tavg_field(tavg_CISO_GLOBAL_D14C,'CISO_GLOBAL_D14C',2,            &
                          long_name='GLOBAL_D14C',            &
                          units='permil', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   buf_len = buf_len+1
   buf_ind_CISO_GLOBAL_D14C = buf_len

!-----------------------------------------------------------------------
!  allocate array for holding flux-related quantities that need to be time-averaged
!  this is necessary since the forcing routines are called before tavg flags
!-----------------------------------------------------------------------

   allocate(ECO_CISO_SFLUX_TAVG(nx_block,ny_block,buf_len,max_blocks_clinic))
   ECO_CISO_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
!  nonstandard 3D fields
!-----------------------------------------------------------------------
   call define_tavg_field(tavg_CISO_PO13C_FLUX_IN,'CISO_PO13C_FLUX_IN',3,&
                          long_name='PO13C Flux into Cell',              &
                          units='mmol/m^3 cm/s', grid_loc='3111',        &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_PO13C_PROD,'CISO_PO13C_PROD',3,      &
                          long_name='PO13C Production',                  &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_PO13C_REMIN,'CISO_PO13C_REMIN',3,    &
                          long_name='PO13C Remineralization',            &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_DO13C_prod,'CISO_DO13C_prod',3,      &
                          long_name='DO13C Production',                  &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_DO13C_remin,'CISO_DO13C_remin',3,    &
                          long_name='DO13C Remineralization',            &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_Ca13CO3_FLUX_IN,'CISO_Ca13CO3_FLUX_IN',3, &
                          long_name='Ca13CO3 flux into cell',                 &
                          units='mmol/m^3 cm/s', grid_loc='3111',             &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_Ca13CO3_PROD,'CISO_Ca13CO3_PROD',3,  &
                          long_name='Ca13CO3 Production',                &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_Ca13CO3_REMIN,'CISO_Ca13CO3_REMIN',3,&
                          long_name='Ca13CO3 Remineralization',          &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_photo13C_TOT,'CISO_photo13C_TOT',3, &
                          long_name='Total 13C Fixation',               &
                          units='mmol/m^3/s', grid_loc='3114',          &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_CISO_DIC_d13C,'CISO_DIC_d13C',3,   &
                          long_name='d13C of DIC',                &
                          units='permil', grid_loc='3111',        &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_DOC_d13C,'CISO_DOC_d13C',3,   &
                          long_name='d13C of DOC',                &
                          units='permil', grid_loc='3111',        &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_zooC_d13C,'CISO_zooC_d13C',3,  &
                          long_name='d13C of zooC',                &
                          units='permil', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

!14C

   call define_tavg_field(tavg_CISO_PO14C_FLUX_IN,'CISO_PO14C_FLUX_IN',3,&
                          long_name='PO14C Flux into Cell',              &
                          units='mmol/m^3 cm/s', grid_loc='3111',        &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_PO14C_PROD,'CISO_PO14C_PROD',3,      &
                          long_name='PO14C Production',                  &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_PO14C_REMIN,'CISO_PO14C_REMIN',3,    &
                          long_name='PO14C Remineralization',            &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_DO14C_prod,'CISO_DO14C_prod',3,      &
                          long_name='DO14C Production',                  &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_DO14C_remin,'CISO_DO14C_remin',3,    &
                          long_name='DO14C Remineralization',            &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_Ca14CO3_FLUX_IN,'CISO_Ca14CO3_FLUX_IN',3, &
                          long_name='Ca14CO3 flux into cell',                 &
                          units='mmol/m^3 cm/s', grid_loc='3111',             &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_Ca14CO3_PROD,'CISO_Ca14CO3_PROD',3,  &
                          long_name='Ca14CO3 Production',                &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_Ca14CO3_REMIN,'CISO_Ca14CO3_REMIN',3,&
                          long_name='Ca14CO3 Remineralization',          &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')


   call define_tavg_field(tavg_CISO_photo14C_TOT,'CISO_photo14C_TOT',3,  &
                          long_name='Total 14C Fixation',                &
                          units='mmol/m^3/s', grid_loc='3114',           &
                          coordinates='TLONG TLAT z_t_150m time')


   call define_tavg_field(tavg_CISO_DIC_d14C,'CISO_DIC_d14C',3,   &
                          long_name='d14C of DIC',                &
                          units='permil', grid_loc='3111',        &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_DOC_d14C,'CISO_DOC_d14C',3,   &
                          long_name='d14C of DOC',                &
                          units='permil', grid_loc='3111',        &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_zooC_d14C,'CISO_zooC_d14C',3, &
                          long_name='d14C of zooC',               &
                          units='permil', grid_loc='3111',        &
                          coordinates='TLONG TLAT z_t time')



!-----------------------------------------------------------------------
!  Nonstandard 2D fields
!-----------------------------------------------------------------------


   call define_tavg_field(tavg_CISO_photo13C_TOT_zint,'CISO_photo13C_TOT_zint',2,&
                          long_name='Total 13C Fixation Vertical Integral',      &
                          units='mmol/m^3 cm/s', grid_loc='2110',                &
                          coordinates='TLONG TLAT time')


   call define_tavg_field(tavg_CISO_Jint_13Ctot,'CISO_Jint_13Ctot',2,            &
                          long_name='13Ctot Source Sink Term Vertical Integral', &
                          units='mmol/m^3 cm/s', grid_loc='2110',                &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_CISO_Jint_100m_13Ctot,'CISO_Jint_100m_13Ctot',2,          &
                          long_name='13Ctot Source Sink Term Vertical Integral, 0-100m', &
                          units='mmol/m^3 cm/s', grid_loc='2110',                        &
                          coordinates='TLONG TLAT time')

! 14C

   call define_tavg_field(tavg_CISO_photo14C_TOT_zint,'CISO_photo14C_TOT_zint',2,&
                          long_name='Total 14C Fixation Vertical Integral',      &
                          units='mmol/m^3 cm/s', grid_loc='2110',                &
                          coordinates='TLONG TLAT time')


   call define_tavg_field(tavg_CISO_Jint_14Ctot,'CISO_Jint_14Ctot',2,            &
                          long_name='14Ctot Source Sink Term Vertical Integral', &
                          units='mmol/m^3 cm/s', grid_loc='2110',                &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_CISO_Jint_100m_14Ctot,'CISO_Jint_100m_14Ctot',2,          &
                          long_name='14Ctot Source Sink Term Vertical Integral, 0-100m', &
                          units='mmol/m^3 cm/s', grid_loc='2110',                        &
                          coordinates='TLONG TLAT time')


!-----------------------------------------------------------------------
!  Nonstandard autotroph 2D and 3D fields for each autotroph
!-----------------------------------------------------------------------


   do auto_ind = 1, autotroph_cnt


      if (autotrophs(auto_ind)%Ca13CO3_ind > 0) then
         sname = 'CISO_' // trim(autotrophs(auto_ind)%sname) // '_Ca13CO3_form'
         call define_tavg_field(tavg_CISO_Ca13CO3_form(auto_ind), sname, 3, &
                                long_name=trim(autotrophs(auto_ind)%lname) // ' Ca13CO3 Formation', &
                                units='mmol/m^3/s', grid_loc='3114', &
                                coordinates='TLONG TLAT z_t_150m time')

         sname = trim(sname) // '_zint'
         call define_tavg_field(tavg_CISO_Ca13CO3_form_zint(auto_ind), sname, 2, &
                                long_name=trim(autotrophs(auto_ind)%lname) // ' Ca13CO3 Formation Vertical Integral', &
                                units='mmol/m^3 cm/s', grid_loc='2110', &
                                coordinates='TLONG TLAT time')
      endif

      if (autotrophs(auto_ind)%Ca14CO3_ind > 0) then
         sname = 'CISO_' // trim(autotrophs(auto_ind)%sname) // '_Ca14CO3_form'
         call define_tavg_field(tavg_CISO_Ca14CO3_form(auto_ind), sname, 3, &
                                long_name=trim(autotrophs(auto_ind)%lname) // ' Ca14CO3 Formation', &
                                units='mmol/m^3/s', grid_loc='3114', &
                                coordinates='TLONG TLAT z_t_150m time')

         sname = trim(sname) // '_zint'
         call define_tavg_field(tavg_CISO_Ca14CO3_form_zint(auto_ind), sname, 2, &
                                long_name=trim(autotrophs(auto_ind)%lname) // ' Ca14CO3 Formation Vertical Integral', &
                                units='mmol/m^3 cm/s', grid_loc='2110', &
                                coordinates='TLONG TLAT time')
      endif


      call define_tavg_field(tavg_CISO_autotrophCaCO3_d13C(auto_ind), &
                            'CISO_autotrophCaCO3_d13C_' // trim(autotrophs(auto_ind)%sname), 3, &
                            long_name=trim(autotrophs(auto_ind)%lname) // ' d13C of CaCO3', &
                            units='mmol/m^3/s', grid_loc='3111', &
                            coordinates='TLONG TLAT z_t time')

      call define_tavg_field(tavg_CISO_autotrophCaCO3_d14C(auto_ind), &
                            'CISO_autotrophCaCO3_d14C_' // trim(autotrophs(auto_ind)%sname), 3, &
                            long_name=trim(autotrophs(auto_ind)%lname) // ' d14C of CaCO3', &
                            units='mmol/m^3/s', grid_loc='3111', &
                            coordinates='TLONG TLAT z_t time')

      call define_tavg_field(tavg_CISO_photo13C(auto_ind), &
                             'CISO_photo13C_' // trim(autotrophs(auto_ind)%sname), 3, &
                             long_name=trim(autotrophs(auto_ind)%lname) // ' 13C Fixation', &
                             units='mmol/m^3/s', grid_loc='3114', &
                             coordinates='TLONG TLAT z_t_150m time')

      call define_tavg_field(tavg_CISO_photo13C_zint(auto_ind), &
                             'CISO_photo13C_' // trim(autotrophs(auto_ind)%sname) // '_zint', 2, &
                             long_name=trim(autotrophs(auto_ind)%lname) // ' 13C Fixation Vertical Integral', &
                             units='mmol/m^3 cm/s', grid_loc='2110', &
                             coordinates='TLONG TLAT time')

      call define_tavg_field(tavg_CISO_eps_autotroph(auto_ind), &
                             'CISO_eps_autotroph_' // trim(autotrophs(auto_ind)%sname), 3, &
                             long_name=trim(autotrophs(auto_ind)%lname) // ' discrimination factor (eps)', &
                             units='permil', grid_loc='3111', &
                             coordinates='TLONG TLAT z_t time')

      call define_tavg_field(tavg_CISO_mui_to_co2star(auto_ind), &
                             'CISO_mui_to_co2star_' // trim(autotrophs(auto_ind)%sname), 3, &
                             long_name=trim(autotrophs(auto_ind)%lname) // ' instanteous growth rate over [CO2*]', &
                             units='m^3/mmol C/s', grid_loc='3111', &
                             coordinates='TLONG TLAT z_t time')


      call define_tavg_field(tavg_CISO_photo14C(auto_ind), &
                             'CISO_photo14C_' // trim(autotrophs(auto_ind)%sname), 3, &
                             long_name=trim(autotrophs(auto_ind)%lname) // ' 14C Fixation', &
                             units='mmol/m^3/s', grid_loc='3114', &
                             coordinates='TLONG TLAT z_t_150m time')

      call define_tavg_field(tavg_CISO_d14C(auto_ind), &
                             'CISO_d14C_' // trim(autotrophs(auto_ind)%sname), 3, &
                             long_name=trim(autotrophs(auto_ind)%lname) // ' d14C', &
                             units='permil', grid_loc='3111', &
                             coordinates='TLONG TLAT z_t time')

       call define_tavg_field(tavg_CISO_d13C(auto_ind), &
                             'CISO_d13C_' // trim(autotrophs(auto_ind)%sname), 3, &
                             long_name=trim(autotrophs(auto_ind)%lname) // ' d13C', &
                             units='permil', grid_loc='3111', &
                             coordinates='TLONG TLAT z_t time')

      call define_tavg_field(tavg_CISO_photo14C_zint(auto_ind), &
                             'CISO_photo14C_' // trim(autotrophs(auto_ind)%sname) // '_zint', 2, &
                             long_name=trim(autotrophs(auto_ind)%lname) // ' 14C Fixation Vertical Integral', &
                             units='mmol/m^3 cm/s', grid_loc='2110', &
                             coordinates='TLONG TLAT time')


   end do


!-----------------------------------------------------------------------
!  More nonstandard 3D fields
!-----------------------------------------------------------------------
   call define_tavg_field(tavg_CISO_eps_aq_g,'CISO_eps_aq_g',3,              &
                          long_name='Equilibrium fractionation (CO2_gaseous <-> CO2_aq)',  &
                          units='permil', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CISO_eps_dic_g,'CISO_eps_dic_g',3,              &
                          long_name='Equilibrium fractionation between total DIC and gaseous CO2',&
                          units='permil', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')



!-----------------------------------------------------------------------
!  Vars to sum up burial in sediments
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_calcToSed_13C,'calcToSed_13C',2,                &
                          long_name='Ca13CO3 Flux to Sediments',         &
                          units='nmolC/cm^2/s', grid_loc='2110',       &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_pocToSed_13C,'pocToSed_13C',2,                  &
                          long_name='PO13C Flux to Sediments',           &
                          units='nmolC/cm^2/s', grid_loc='2110',       &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_calcToSed_14C,'calcToSed_14C',2,                &
                          long_name='Ca14CO3 Flux to Sediments',         &
                          units='nmolC/cm^2/s', grid_loc='2110',       &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_pocToSed_14C,'pocToSed_14C',2,                  &
                          long_name='PO14C Flux to Sediments',           &
                          units='nmolC/cm^2/s', grid_loc='2110',       &
                          coordinates='TLONG TLAT time')


!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_ciso_init_tavg


!***********************************************************************
!BOP
! !IROUTINE: ecosys_ciso_init_sflux
! !INTERFACE:

 subroutine ecosys_ciso_init_sflux

! !USES:
   use named_field_mod, only: named_field_get_index
   use registry, only: registry_match


! !DESCRIPTION:
!  Initialize surface flux computations for the ecosys_ciso tracer module.
!  Includes reading CO2 and D13C data from file if option file is used
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-------------------------------------------------------------------------
!     Set D13C data source
!-------------------------------------------------------------------------

   select case (ciso_atm_d13c_opt)

   case ('const')
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*)'ciso: Using constant D13C values of ',ciso_atm_d13c_const
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
   endif

!-----------------------------------------------------------------------
!     READ in D13C data from file
!-----------------------------------------------------------------------

   case('file')
   call ciso_read_atm_D13C_data

   case default
   call exit_POP(sigAbort, 'unknown ciso_atm_d13c_opt in ecosys_ciso_init_sflux')

   end select
!-------------------------------------------------------------------------
!     Set D14C data source
!-------------------------------------------------------------------------

   select case (ciso_atm_d14c_opt)

   case ('const')
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*)'ciso: Using constant D14C values of ',ciso_atm_d14c_const
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
   endif

!-----------------------------------------------------------------------
!     READ in D14C data from files
!-----------------------------------------------------------------------

   case('file')
   call ciso_read_atm_D14C_data

   case default
   call exit_POP(sigAbort, 'unknown ciso_atm_d14c_opt in ecosys_ciso_init_sflux')

   end select


!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_ciso_init_sflux


!*****************************************************************************
!BOP
! !IROUTINE: ecosys_ciso_set_sflux
!
! !INTERFACE:

 subroutine ecosys_ciso_set_sflux(SST,SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)


! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      SST           ! sea surface temperature (C)


   real (r8), dimension(:,:,:,:), intent(in) :: &
           SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:,:), intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   character(*), parameter :: subname = 'ecosys_ciso_mod:ecosys_ciso_set_sflux'

   logical (log_kind), save :: &
      first = .true.  ! Logical for first iteration test

   type (block) :: &
      this_block      ! block info for the current block

   integer (int_kind) :: &
      i,j,iblock,n, & ! loop indices
      auto_ind,     & ! autotroph functional group index
      mcdate,sec,   & ! date vals for shr_strdata_advance
      errorCode       ! errorCode from HaloUpdate call

   real (r8), dimension(nx_block) :: &
      DI13C_ROW,    & ! row of DI13C values for solver
      DI14C_ROW,    & ! row of DI13C values for solver
      DIC_ROW,      & ! row of DIC values for solver
      DCO2STAR_ROW, & ! DCO2STAR from solver
      CO2STAR_ROW,  & ! CO2STAR from solver
      PV_ROW


   real (r8), dimension(nx_block,ny_block) :: &
      D13C,          &   ! atm 13co2 value
      R13C_DIC_surf, &   ! 13C/12C ratio in DIC
      R13C_atm,      &   ! 13C/12C ratio in atmospheric CO2
      FLUX13,        &   ! gas flux of 13CO2 (nmol/cm^2/s)
      FLUX13_as,     &   ! air-to-sea gas flux of 13CO2 (nmol/cm^2/s)
      FLUX13_sa,     &   ! sea-to-air gas flux of 13CO2 (nmol/cm^2/s)
      D14C,          &   ! atm 14co2 value
      R14C_DIC_surf, &   ! 14C/12C ratio in total DIC
      R14C_atm,      &   ! 14C/12C ratio in atmospheric CO2
      FLUX14,        &   ! gas flux of 14CO2 (nmol/cm^2/s)
      FLUX14_as,     &   ! air-to-sea gas flux of 14CO2 (nmol/cm^2/s)
      FLUX14_sa,     &   ! sea-to-air gas flux of 14CO2 (nmol/cm^2/s)
      FLUX,          &   ! gas flux of CO2 (nmol/cm^2/s)
      FLUX_as,       &   ! air-to-sea gas flux of CO2 (nmol/cm^2/s)
      FLUX_sa            ! sea-to-air gas flux of CO2 (nmol/cm^2/s)


   real (r8), dimension(nx_block,ny_block) :: &
      eps_aq_g_surf,       & ! equilibrium fractionation (CO2_gaseous <-> CO2_aq)
      alpha_aq_g_surf,     & ! alpha_xxx_g_surf => eps = ( alpa -1 ) * 1000
      eps_dic_g_surf,      & ! equilibrium fractionation between total DIC and gaseous CO2
      alpha_dic_g_surf,    & ! alpha_xxx_g_surf => eps = ( alpa -1 ) * 1000
      eps_hco3_g_surf,     & ! equilibrium fractionation between bicarbonate and gaseous CO2
      eps_co3_g_surf,      & ! equilibrium fractionation between carbonate and gaseous CO2
      frac_co3,            & ! carbonate fraction fCO3 = [CO3--]/DIC
      alpha_aq_g_surf_14c, & ! for 14C, with fractionation being twice as large for 14C than for 13C
      alpha_dic_g_surf_14c   ! for 14C, with fractionation being twice as large for 14C than for 13C


   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      di13c_riv_flux,   & ! River input of DI13C
      do13c_riv_flux,   & ! River input of DO13C
      di14c_riv_flux,   & ! River input of DI14C
      do14c_riv_flux      ! River input of DO14C

! For making global average of atmospheric D14C
   real (r8), dimension(nx_block,ny_block) :: &
      WORK1, &! local work space
      TFACT   ! factor for normalizing sums

   integer (int_kind) :: &
      ib,ie,jb,je

   real (r8), dimension(max_blocks_clinic) :: &
      D14C_local_sums, & ! array for holding block sums when calculating global D14C
      TAREA_local_sums   ! array for holding block sums of TAREA when calculating global D14C

   real (r8) :: &
      D14C_sum_tmp,  & ! temp for local sum of D14C
      TAREA_sum_tmp    ! temp for local sum of TAREA

   real (r8) :: &
      D14C_glo_avg  ! global average D14C over the ocean, computed from current D14C field

!-----------------------------------------------------------------------
!     local parameters for 13C
!     Zhang et al, 1995, Geochim. et Cosmochim. Acta, 59 (1), 107-114
!-----------------------------------------------------------------------
   real(r8) :: &
      alpha_k, &             ! eps = ( alpa -1 ) * 1000
      alpha_k_14c            ! for 14C, with fractionation being twice as large for 14C than for 13C


   real(r8), parameter :: &
      eps_k     = -0.81_r8  ! kinetic fraction during gas
                             ! transfert (per mil) (air-sea CO2
                             ! exchange) at 21C, Zhang et al 1995,
                             ! eps_k = -0.95 at 5C
 !-----------------------------------------------------------------------
 ! The following variables from ecosys_share are used directly, without
 ! defining pointers:
 !   CO2STAR_SURF_fields
 !   DCO2STAR_SURF_fields
 !   PV_SURF_fields
 !   DIC_SURF_fields
 !   CO3_SURF_fields
 !-----------------------------------------------------------------------

!-----------------------------------------------------------------------

   call timer_start(ecosys_ciso_sflux_timer)
!-----------------------------------------------------------------------

   if (check_time_flag(ciso_comp_surf_avg_flag))     &
      call comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR,&
                         ecosys_ciso_tracer_cnt,     &
                         ciso_vflux_flag,ciso_surf_avg)



   if (first) then
      allocate( ciso_data_ind_d13c(max_blocks_clinic) )
      allocate( ciso_data_ind_d14c(max_blocks_clinic) )
      ciso_data_ind_d13c = -1
      ciso_data_ind_d14c = -1
      first = .false.
   endif

!-----------------------------------------------------------------------
!  fluxes initially set to 0
!-----------------------------------------------------------------------
   WORK1            = c0
   D14C_local_sums  = c0
   TAREA_local_sums = c0

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
      STF_MODULE(:,:,:,iblock) = c0
   end do
  !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
   !$OMP PARALLEL DO PRIVATE(iblock, j, this_block, ib, ie, jb, je, D13C, R13C_atm,   &
   !$OMP                     DI13C_ROW, R13C_DIC_surf, FLUX13, FLUX13_as, &
   !$OMP                     FLUX13_sa, D14C, R14C_atm, DI14C_ROW, &
   !$OMP                     R14C_DIC_surf, FLUX14, FLUX14_as, FLUX14_sa, &
   !$OMP                     DIC_ROW, CO2STAR_ROW, DCO2STAR_ROW, FLUX, &
   !$OMP                     FLUX_as, FLUX_sa, PV_ROW, eps_aq_g_surf, &
   !$OMP                     eps_dic_g_surf, frac_co3, alpha_k, &
   !$OMP                     alpha_aq_g_surf, alpha_dic_g_surf, &
   !$OMP                     alpha_k_14c, alpha_aq_g_surf_14c, &
   !$OMP                     alpha_dic_g_surf_14c,TFACT,WORK1)
!-----------------------------------------------------------------------


   do iblock = 1, nblocks_clinic

!-----------------------------------------------------------------------
!  Set D13C and D14C (constant or from files read in _init) and put on global grid
!-----------------------------------------------------------------------
      select case (ciso_atm_d13c_opt)

      case ('const')
      D13C = ciso_atm_d13c_const

      case ('file')
      call ciso_comp_varying_D13C(iblock, ciso_data_ind_d13c(iblock),D13C)

      case default
      call exit_POP(sigAbort, 'unknown ciso_atm_d13c_opt in ecosys_ciso_set_sflux')

      end select

      select case (ciso_atm_d14c_opt)

      case ('const')
      D14C = ciso_atm_d14c_const

      case ('file')
      call ciso_comp_varying_D14C(iblock, ciso_data_ind_d14c(iblock),D14C)

      case default
      call exit_POP(sigAbort, 'unknown ciso_atm_d14c_opt in ecosys_ciso_set_sflux')

      end select
!-----------------------------------------------------------------------
! Save local D14C field for making global mean after end of iblock loop
!-----------------------------------------------------------------------

      this_block = get_block(blocks_clinic(iblock),iblock)
      ib = this_block%ib
      ie = this_block%ie
      jb = this_block%jb
      je = this_block%je

      where (LAND_MASK(:,:,iblock))
        TFACT = TAREA(:,:,iblock)
      elsewhere
        TFACT = 0.0_r8
      endwhere

      WORK1 = D14C * TFACT
      D14C_local_sums(iblock) = sum(WORK1(ib:ie,jb:je))
      TAREA_local_sums(iblock) = sum(TFACT(ib:ie,jb:je))


!-----------------------------------------------------------------------
!     initialize R13C_atm  and R14C_atm
!-----------------------------------------------------------------------

      R13C_atm = R13C_std * ( c1 + D13C / c1000 )

      R14C_atm = R14C_std * ( c1 + D14C / c1000 )

!-----------------------------------------------------------------------
!     compute 13C02 flux, based on CO2 flux calculated in ecosystem model
!     Zhang et al, 1995, Geochim. et Cosmochim. Acta, 59 (1), 107-114
!-----------------------------------------------------------------------

      do j = 1,ny_block
!-----------------------------------------------------------------------
!     compute R13C_DIC in surface ocean (assuming that DIC_ROW is 12C)
!-----------------------------------------------------------------------
         DI13C_ROW = p5*(SURF_VALS_OLD(:,j,di13c_ind,iblock) + &
                     SURF_VALS_CUR(:,j,di13c_ind,iblock))
         DI14C_ROW = p5*(SURF_VALS_OLD(:,j,di14c_ind,iblock) + &
                     SURF_VALS_CUR(:,j,di14c_ind,iblock))

         DIC_ROW   = DIC_SURF_fields(:,j,iblock)

         where ( DIC_ROW /= c0 )
            R13C_DIC_surf(:,j) = DI13C_ROW /DIC_ROW
         elsewhere
            R13C_DIC_surf(:,j) = c0
         endwhere

         where ( DIC_ROW /= c0 )
            R14C_DIC_surf(:,j) = DI14C_ROW / DIC_ROW
         elsewhere
            R14C_DIC_surf(:,j) = c0
         endwhere


!-----------------------------------------------------------------------
!     individal discrimination factor of each species with respect to
!     gaseous CO2, temperature dependent, based on Zhang et al. 95
!-----------------------------------------------------------------------
         eps_aq_g_surf(:,j)   = 0.0049_r8 * SST(:,j,iblock) - 1.31_r8
!!         eps_hco3_g_surf(:,j) = -0.141_r8  * SST(:,j,iblock) + 10.78_r8
!!         eps_co3_g_surf(:,j)  = -0.052_r8  * SST(:,j,iblock) + 7.22_r8


!-----------------------------------------------------------------------
!     compute the equilibrium discrimination factor between DIC and
!     gaseous CO2
!-----------------------------------------------------------------------
!     solution 1 : from individual species.
!     Not used as Zhang et al. 95
!     concluded that eps_dic_g_surf can not be calculated from the sum of
!     the three individual species
!----------------------------------------------------------------------
!         eps_dic_g_surf(:,j) = eps_aq_g_surf(:,j) + eps_hco3_g_surf(:,j) &
!                                + eps_co3_g_surf(:,j)
!-----------------------------------------------------------------------
!     solution 2: function of T and carbonate fraction (frac_co3)
!     Using this one, which is based on the empirical relationship from
!     the measured e_dic_g_surf of Zhang et al. 1995
!---------------------------------------------------------------------

         where (.not. LAND_MASK(:,j,iblock))
            frac_co3(:,j) = c0
         elsewhere
            frac_co3(:,j) = CO3_SURF_fields(:,j,iblock) / DIC_ROW
         end where

         eps_dic_g_surf(:,j)  = 0.014_r8 * SST(:,j,iblock) * frac_co3(:,j) - &
                                0.105_r8 * SST(:,j,iblock) + 10.53_r8

!-----------------------------------------------------------------------
!     compute alpha coefficients from eps :  eps = ( alpha -1 ) * 1000
!     => alpha = 1 + eps / 1000
!-----------------------------------------------------------------------

         alpha_k               = c1 + eps_k               / c1000
         alpha_aq_g_surf(:,j)  = c1 + eps_aq_g_surf(:,j)  / c1000
         alpha_dic_g_surf(:,j) = c1 + eps_dic_g_surf(:,j) / c1000

! Fractionation is twice as large for 14C than for 13C, so eps needs to be multiplied by 2 for 14C
         alpha_k_14c               = c1 + eps_k * 2.0_r8              / c1000
         alpha_aq_g_surf_14c(:,j)  = c1 + eps_aq_g_surf(:,j) *2.0_r8  / c1000
         alpha_dic_g_surf_14c(:,j) = c1 + eps_dic_g_surf(:,j) *2.0_r8 / c1000

!-----------------------------------------------------------------------
!     compute 13C flux and C flux
!-----------------------------------------------------------------------
         CO2STAR_ROW  = CO2STAR_SURF_fields(:,j,iblock)
         DCO2STAR_ROW = DCO2STAR_SURF_fields(:,j,iblock)
         PV_ROW       = PV_SURF_fields(:,j,iblock)

         FLUX13(:,j) = PV_ROW * alpha_k * alpha_aq_g_surf(:,j) * &
                      (( CO2STAR_ROW + DCO2STAR_ROW ) * R13C_atm(:,j) - &
                      CO2STAR_ROW * R13C_DIC_surf(:,j) / alpha_dic_g_surf(:,j) )

         FLUX14(:,j) = PV_ROW * alpha_k_14c * alpha_aq_g_surf_14c(:,j) * &
                      (( CO2STAR_ROW + DCO2STAR_ROW ) * R14C_atm(:,j) - &
                      CO2STAR_ROW * R14C_DIC_surf(:,j) / alpha_dic_g_surf_14C(:,j) )

         FLUX(:,j)   = PV_ROW * DCO2STAR_ROW



!-----------------------------------------------------------------------
!     compute fluxes in and out
!-----------------------------------------------------------------------

         FLUX_as(:,j)   = PV_ROW * ( DCO2STAR_ROW + CO2STAR_ROW )
         FLUX_sa(:,j)   = PV_ROW * CO2STAR_ROW

         FLUX13_as(:,j) = PV_ROW * alpha_k * alpha_aq_g_surf(:,j) * &
                          (( CO2STAR_ROW + DCO2STAR_ROW ) * R13C_atm(:,j))

         FLUX13_sa(:,j) = PV_ROW * alpha_k * alpha_aq_g_surf(:,j) * &
                          ( CO2STAR_ROW * R13C_DIC_surf(:,j) / alpha_dic_g_surf(:,j) )

         FLUX14_as(:,j) = PV_ROW * alpha_k_14c * alpha_aq_g_surf_14c(:,j) * &
                          (( CO2STAR_ROW + DCO2STAR_ROW ) * R14C_atm(:,j))

         FLUX14_sa(:,j) = PV_ROW * alpha_k_14c * alpha_aq_g_surf_14c(:,j) * &
                          ( CO2STAR_ROW * R14C_DIC_surf(:,j) / alpha_dic_g_surf_14c(:,j) )

!-----------------------------------------------------------------------
!     end of 13C computation for gass exchange
!-----------------------------------------------------------------------

      end do !j loop

!-----------------------------------------------------------------------
!     Adding 13C FLux to total DI13C
!-----------------------------------------------------------------------

      STF_MODULE(:,:,di13c_ind,iblock) = STF_MODULE(:,:,di13c_ind,iblock) + FLUX13
      STF_MODULE(:,:,di14c_ind,iblock) = STF_MODULE(:,:,di14c_ind,iblock) + FLUX14

!-----------------------------------------------------------------------
!    Tavg variables
!-----------------------------------------------------------------------

      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_13CO2,iblock)    = FLUX13
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_as_13CO2,iblock) = FLUX13_as
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_sa_13CO2,iblock) = FLUX13_sa

      where ( FLUX /= c0 )
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d13C,iblock) = ( FLUX13 / FLUX  &
                                              / R13C_std - c1 ) * c1000
      elsewhere
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d13C,iblock) = c0
      endwhere

      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R13C_DIC_surf,iblock)  = R13C_DIC_surf
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R13C_atm,iblock)  = R13C_atm
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_D13C_atm,iblock) = D13C
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_eps_aq_g_surf,iblock) = eps_aq_g_surf
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_eps_dic_g_surf,iblock) = eps_dic_g_surf
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_14CO2,iblock) = FLUX14
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_as_14CO2,iblock) = FLUX14_as
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_sa_14CO2,iblock) = FLUX14_sa
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R14C_DIC_surf,iblock) = R14C_DIC_surf
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R14C_atm,iblock) = R14C_atm
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_D14C_atm,iblock) = D14C

      where ( FLUX /= c0 )
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d14C,iblock) = ( FLUX14 / FLUX &
                                               / R14C_std - c1 ) * c1000
      elsewhere
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d14C,iblock) = c0
      endwhere

   enddo ! end of i-lopp
!-----------------------------------------------------------------------
!$OMP END PARALLEL DO
!-----------------------------------------------------------------------


!-------------------------------------------------------------------------
! River input of isotopic DIC and DOC.
! River input of BGC tracers in ecosys_mod is currently constant and from file
! So the isotopic carbon input is also done very simplified with one value
! globally, even though data shows it should vary from river to river.
!
! Using constant delta values of
! D13C=-10 permil for DIC (Mook 1986, Raymond et al 2004)
! D13C=-27.6 permil for DOC (Raymond et al 2004)
! D14C=-50 permil for DOC (Raymond et al 2004), Gruber et al
! D14C= atmos_D14C - 50 permil for DIC (based on very few data points and 
!       discussion with N. Gruber)
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!  Make average global DI14C
!-----------------------------------------------------------------------------

   D14C_sum_tmp  = sum(D14C_local_sums)
   TAREA_sum_tmp = sum(TAREA_local_sums)


   D14C_glo_avg  = global_sum(D14C_sum_tmp,distrb_clinic) / &
                   global_sum(TAREA_sum_tmp,distrb_clinic)


   di13c_riv_flux = dic_riv_flux_fields * (-10.0_r8/c1000 +c1) * R13C_std
   di14c_riv_flux = dic_riv_flux_fields * ((D14C_glo_avg - 50.0_r8)/c1000 +c1) * R14C_std

   do13c_riv_flux = doc_riv_flux_fields * (-27.6_r8/c1000 +c1) * R13C_std
   do14c_riv_flux = doc_riv_flux_fields * (-50.0_r8/c1000 +c1) * R14C_std

   STF_MODULE(:,:,di13c_ind,:) = STF_MODULE(:,:,di13c_ind,:) + di13c_riv_flux
   STF_MODULE(:,:,do13c_ind,:) = STF_MODULE(:,:,do13c_ind,:) + do13c_riv_flux

   STF_MODULE(:,:,di14c_ind,:) = STF_MODULE(:,:,di14c_ind,:) + di14c_riv_flux
   STF_MODULE(:,:,do14c_ind,:) = STF_MODULE(:,:,do14c_ind,:) + do14c_riv_flux

! write to tavg
   ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DI13C_RIV_FLUX,:) = di13c_riv_flux
   ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DO13C_RIV_FLUX,:) = do13c_riv_flux

   ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DI14C_RIV_FLUX,:) = di14c_riv_flux
   ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DO14C_RIV_FLUX,:) = do14c_riv_flux

   ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_GLOBAL_D14C,:) = D14C_glo_avg

 call timer_stop(ecosys_ciso_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_ciso_set_sflux

!***********************************************************************
!BOP
! !IROUTINE: ecosys_ciso_set_interior
! !INTERFACE:

 subroutine ecosys_ciso_set_interior(k, TEMP_OLD, TEMP_CUR, &
    TRACER_MODULE_OLD, TRACER_MODULE_CUR, DTRACER_MODULE, this_block)

! !DESCRIPTION:
!  Compute time derivatives for 13C and 14C state variables.
!
! !REVISION HISTORY:
!  13C code is based on code from X. Giraud, ETH Zürich, 2008, for pop1
!  Adapted to pop2 and new ecosystem model code and added biotic 14C
!  by A. Jahn, NCAR, 2012-2013
!
!


! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                   ! vertical level index

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      TEMP_OLD,          &! old potential temperature (C)
      TEMP_CUR            ! current potential temperature (C)

   real (r8), dimension(:,:,:,:), intent(in) :: &
      TRACER_MODULE_OLD, &! old tracer values
      TRACER_MODULE_CUR   ! current tracer values

   type (block), intent(in) :: &
      this_block          ! block info for the current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(out) :: &
      DTRACER_MODULE      ! computed source/sink terms

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'ecosys_ciso_mod:ecosys_ciso_set_interior'

   real (r8), dimension(nx_block,ny_block) :: &
      TEMP              ! local copy of model TEMP

   logical (log_kind), dimension(nx_block,ny_block) :: ZERO_MASK

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1,WORK2     ! temporaries

   real (r8) :: &
      ztop              ! depth of top of cell

   integer (int_kind) :: &
      bid,            & ! local_block id
      n,m,            & ! tracer index
      auto_ind,       & ! autotroph functional group index
      kk                ! index for looping over k levels

!-------------------------------------------------------------
! 13C Variables added
!-------------------------------------------------------------
   type(sinking_particle), save :: &
      PO13C,          & ! base units = nmol 13C
      P_Ca13CO3,      & ! base units = nmol CaCO3 13C
      PO14C,          & ! base units = nmol 14C
      P_Ca14CO3         ! base units = nmol CaCO3 14C
!     POC               ! base units = nmol C -> Defined in ecosys_share
!     P_CaCO3           ! base units = nmol CaCO3 -> Defined in ecosys_share

   real (r8), dimension(nx_block,ny_block) :: &
      DO13C_loc,         & ! local copy of model DO13C
      DI13C_loc,         & ! local copy of model DI13C
      zoo13C_loc,        & ! local copy of model zoo13C
      DO14C_loc,         & ! local copy of model DO14C
      DI14C_loc,         & ! local copy of model DI14C
      zoo14C_loc,        & ! local copy of model zoo14C
      R13C_CaCO3_PROD,   & ! 13C/12C in CaCO3 production of small phyto
      R13C_CO2STAR,      & ! 13C/12C in CO2* water
      R13C_DIC,          & ! 13C/12C in total DIC
      R13C_DOC,          & ! 13C/12C in total DOC
      R13C_zooC,         & ! 13C/12C in total zooplankton
      R14C_CaCO3_PROD,   & ! 14C/12C in CaCO3 production of small phyto
      R14C_CO2STAR,      & ! 14C/12C in CO2* water
      R14C_DIC,          & ! 14C/12C in total DIC
      R14C_DOC,          & ! 14C/12C in total DOC
      R14C_zooC            ! 14C/12C in total zooplankton


   real (r8), dimension(nx_block,ny_block,autotroph_cnt) :: &
      Ca13CO3_PROD,        & ! prod. of 13C CaCO3 by small phyto (mmol CaCO3/m^3/sec)
      Ca14CO3_PROD,        & ! prod. of 13C CaCO3 by small phyto (mmol CaCO3/m^3/sec)
      eps_autotroph,       & ! Permil fractionation (or discrimination factor) for Carbon autotroph types sp, diat, diaz
      mui_to_co2star,      & ! Carbon autotroph instanteous growth rate over [CO2*] (m^3 /mmol C /s)
      R13C_photoC,         & ! 13C/12C in Carbon autotroph C-fixation (mmol C/m^3/sec)
      R13C_autotroph,      & ! 13C/12C in total small phytoplankton
      photo13C,            & ! Carbon autotroph 13C-fixation (mmol C/m^3/sec)
      R14C_photoC,         & ! 14C/12C in Carbon autotroph C-fixation (mmol C/m^3/sec)
      R14C_autotroph,      & ! 14C/12C in total small phytoplankton
      photo14C,            & ! Carbon autotroph 14C-fixation (mmol C/m^3/sec)
      autotrophCaCO3_d13C, & ! d13C of autotrophCaCO3
      autotrophCaCO3_d14C, & ! d14C of autotrophCaCO3
      autotroph_d13C,      & ! d13C of autotroph C
      autotroph_d14C,      & ! d14C of autotroph C
      autotroph13C_loc,    &  ! local copy of model autotroph 13C
      autotroph14C_loc,    & ! local copy of model autotroph 14C
      autotrophCa13CO3_loc,& ! local copy of model autotrophCa13CO3
      autotrophCa14CO3_loc,& ! local copy of model autotrophCa14CO3
      R13C_autotrophCaCO3, & ! 13C/12C in total small phytoplankton carbonate
      R14C_autotrophCaCO3         ! 14C/12C in total small phytoplankton carbonate


   real (r8), dimension(nx_block,ny_block) :: &
      mui_to_co2star_loc,& ! Carbon autotroph instanteous growth rate over [CO2*] (m^3 /mol C /s)
      frac_co3,          & ! carbonate fraction fCO3 = [CO3--]/DIC
      CO2STAR_int,       & ! [CO2*] water (mmol/m^3) in interior domain (not only surface)
      DO13C_prod,        & ! production of 13C DOC (mmol C/m^3/sec)
      DO13C_remin,       & ! remineralization of 13C DOC (mmol C/m^3/sec)
      eps_aq_g,          & ! equilibrium fractionation (CO2_gaseous <-> CO2_aq)
      eps_dic_g,         & ! equilibrium fractionation between total DIC and gaseous CO2
      alpha_aq_g,        & ! eps = ( alpa -1 ) * 1000
      alpha_dic_g,       & ! eps = ( alpa -1 ) * 1000
      delta_C13_Corg,    & ! deltaC13 of Net Primary Production
      delta_C13_CO2STAR, & ! deltaC13 of CO2*
      DO14C_prod,        & ! production of 13C DOC (mmol C/m^3/sec)
      DO14C_remin,       & ! remineralization of 13C DOC (mmol C/m^3/sec)
      alpha_aq_g_14c,    & ! alpha for 14C, with fractionation twice as large as for 13C
      alpha_dic_g_14c,   & ! alpha for 14C, with fractionation twice as large as for 13C
      delta_C14_Corg,    & ! deltaC14 of Net Primary Production
      delta_C14_CO2STAR, & ! deltaC14 of CO2*
      DIC_d13C,          & ! d13C of DIC
      DOC_d13C,          & ! d13C of DOC
      zooC_d13C,         & ! d13C of zooC
      DIC_d14C,          & ! d14C of DIC
      DOC_d14C,          & ! d14C of DOC
      zooC_d14C            ! d14C of zooC


   real (r8), dimension(autotroph_cnt) :: &
      cell_active_C_uptake,  & ! ratio of active carbon uptake to carbon fixation
      cell_active_C,         & ! ratio of active carbon uptake to carbon fixation
      cell_surf,             & ! surface areas of cells ( m2 )
      cell_carb_cont,        & ! cell carbon content ( mol C cell-1 )
      cell_radius,           & ! cell radius ( um )
      cell_permea,           & ! cell wall permeability to CO2(aq) (m/s)
      cell_eps_fix             ! fractionation effect of carbon fixation



   real(r8), parameter :: &
      eps_carb = -2.0_r8          ! eps_carb = d13C(CaCO3) - d13C(DIC)  Ziveri et al., 2003

!---------------------------------------------------------------
! Define pointer variables, used to share values with other modules
! Target variables are defined in ecosys_share
! Below pointers are used to point to the right part of the
! global array in ecosys_share
!---------------------------------------------------------------

   real (r8), dimension(:,:), pointer :: DIC_loc      ! local copy of model DIC
   real (r8), dimension(:,:), pointer :: DOC_loc      ! local copy of model DOC
   real (r8), dimension(:,:), pointer :: O2_loc       ! local copy of model O2
   real (r8), dimension(:,:), pointer :: NO3_loc      ! local copy of model NO3
   real (r8), dimension(:,:), pointer :: CO3          ! carbonate ion
   real (r8), dimension(:,:), pointer :: HCO3         ! bicarbonate ion
   real (r8), dimension(:,:), pointer :: H2CO3        ! carbonic acid
   real (r8), dimension(:,:), pointer :: DOC_remin    ! remineralization of 13C DOC (mmol C/m^3/sec)

   real (r8), dimension(:,:,:), pointer :: zooC_loc     ! local copy of model zooC
   real (r8), dimension(:,:,:), pointer :: zoo_loss     ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: zoo_loss_poc ! zoo_loss routed to large detrital pool (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: zoo_loss_doc ! zoo_loss routed to doc (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: zoo_loss_dic ! zoo_loss routed to dic (mmol C/m^3/sec)

   real (r8), dimension(:,:,:), pointer :: QCaCO3               ! small phyto CaCO3/C ratio (mmol CaCO3/mmol C)
   real (r8), dimension(:,:,:), pointer :: autotrophCaCO3_loc  ! local copy of model autotroph CaCO3
   real (r8), dimension(:,:,:), pointer :: autotrophChl_loc    ! local copy of model autotroph Chl
   real (r8), dimension(:,:,:), pointer :: autotrophC_loc      ! local copy of model autotroph C
   real (r8), dimension(:,:,:), pointer :: autotrophFe_loc     ! local copy of model autotroph Fe
   real (r8), dimension(:,:,:), pointer :: autotrophSi_loc     ! local copy of model autotroph Si
   real (r8), dimension(:,:,:), pointer :: auto_graze          ! autotroph grazing rate (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_graze_zoo      ! auto_graze routed to zoo (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_graze_poc      ! auto_graze routed to poc (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_graze_doc      ! auto_graze routed to doc (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_graze_dic      ! auto_graze routed to dic (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_loss           ! autotroph non-grazing mort (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_loss_poc       ! auto_loss routed to poc (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_loss_doc       ! auto_loss routed to doc (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_loss_dic       ! auto_loss routed to dic (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: auto_agg            ! autotroph aggregation (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: photoC              ! C-fixation (mmol C/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: CaCO3_PROD          ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
   real (r8), dimension(:,:,:), pointer :: PCphoto             ! C-specific rate of photosynth. (1/sec)


!-------------------------------------------------------------

   bid = this_block%local_id

!-----------------------------------------------------------------------

   call timer_start(ecosys_ciso_interior_timer, block_id=bid)

   DTRACER_MODULE = c0

!-----------------------------------------------------------------------
!  exit immediately if computations are not to be performed
!-----------------------------------------------------------------------

   if (.not. ciso_lsource_sink) then
      call timer_stop(ecosys_ciso_interior_timer, block_id=bid)
      return
   endif

!-------------------------------------------------------------
! Assign locally used variables to pointer variables
! => use pointers to point to the right part of the global array
! in ecosys_share
!---------------------------------------------------------------

   DIC_loc => DIC_loc_fields(:,:,bid)
   DOC_loc => DOC_loc_fields(:,:,bid)
   O2_loc  => O2_loc_fields(:,:,bid)
   NO3_loc  => NO3_loc_fields(:,:,bid)
   CO3 => CO3_fields(:,:,bid)
   HCO3 => HCO3_fields(:,:,bid)
   H2CO3 => H2CO3_fields(:,:,bid)
   DOC_remin => DOC_remin_fields(:,:,bid)

   zooC_loc => zooC_loc_fields(:,:,:,bid)
   zoo_loss => zoo_loss_fields(:,:,:,bid)
   zoo_loss_poc => zoo_loss_poc_fields(:,:,:,bid)
   zoo_loss_doc => zoo_loss_doc_fields(:,:,:,bid)
   zoo_loss_dic => zoo_loss_dic_fields(:,:,:,bid)

   QCaCO3 => QCaCO3_fields(:,:,:,bid)
   autotrophCaCO3_loc => autotrophCaCO3_loc_fields(:,:,:,bid)
   autotrophChl_loc => autotrophChl_loc_fields(:,:,:,bid)
   autotrophC_loc => autotrophC_loc_fields(:,:,:,bid)
   autotrophFe_loc => autotrophFe_loc_fields(:,:,:,bid)
   autotrophSi_loc => autotrophSi_loc_fields(:,:,:,bid)
   auto_graze => auto_graze_fields(:,:,:,bid)
   auto_graze_zoo => auto_graze_zoo_fields(:,:,:,bid)
   auto_graze_poc => auto_graze_poc_fields(:,:,:,bid)
   auto_graze_doc => auto_graze_doc_fields(:,:,:,bid)
   auto_graze_dic => auto_graze_dic_fields(:,:,:,bid)
   auto_loss => auto_loss_fields(:,:,:,bid)
   auto_loss_poc => auto_loss_poc_fields(:,:,:,bid)
   auto_loss_doc => auto_loss_doc_fields(:,:,:,bid)
   auto_loss_dic => auto_loss_dic_fields(:,:,:,bid)
   auto_agg => auto_agg_fields(:,:,:,bid)
   photoC => photoC_fields(:,:,:,bid)
   CaCO3_PROD => CaCO3_PROD_fields(:,:,:,bid)
   PCphoto => PCphoto_fields(:,:,:,bid)

!----------------------------------------------------------------------------------------
! For Keller and Morel, set cell attributes based on autotroph type (from observations)
!----------------------------------------------------------------------------------------
   select case (ciso_fract_factors)
      case ('KellerMorel')
         do auto_ind = 1, autotroph_cnt
            if (autotrophs(auto_ind)%kSiO3 > c0) then
               !----------------------------------------------------------------------------------------
               ! Diatom based on P. tricornumtum ( Keller and morel, 1999; Popp et al., 1998 )
               !----------------------------------------------------------------------------------------
               cell_active_C_uptake(auto_ind) = 2.3_r8       ! ratio of active carbon uptake to carbon fixation
               cell_surf(auto_ind)            = 100.6e-12_r8 ! surface areas of cells ( m2 )
               cell_carb_cont(auto_ind)       = 63.3e-14_r8  ! cell carbon content ( mol C cell-1 )
               cell_radius(auto_ind)          = 14.2_r8      ! cell radius ( um )
               cell_permea(auto_ind)          = 3.3e-5_r8    ! cell wall permeability to CO2(aq) (m/s)
               cell_eps_fix(auto_ind)         = 26.6_r8      ! fractionation effect of carbon fixation

            else if (autotrophs(auto_ind)%Nfixer) then
               !----------------------------------------------------------------------------------------
               ! Diazotroph based on  Standard Phyto of Rau et al., (1996)
               !----------------------------------------------------------------------------------------
               !cell_active_C_uptake(auto_ind) = 0.0_r8        ! ratio of active carbon uptake to carbon fixation
               !cell_surf(auto_ind)            = -99.9_r8      ! surface areas of cells ( m2 ) - not used -
               !cell_carb_cont(auto_ind)       = -99.9_r8      ! cell carbon content ( mol C cell-1 ) - not used -
               !cell_radius(auto_ind)          = 10.0_r8       ! cell radius ( um )
               !cell_permea(auto_ind)          = 10.0e-5_r8    ! cell wall permeability to CO2(aq) (m/s)
               !cell_eps_fix(auto_ind)         = 25.0_r8       ! fractionation effect of carbon fixation
               !----------------------------------------------------------------------------------------
               ! Diazotroph based on Synechococcus sp. ( Keller and morel, 1999; Popp et al., 1998 )
               !----------------------------------------------------------------------------------------
               cell_active_C_uptake(auto_ind) = 7.5_r8        ! ratio of active carbon uptake to carbon fixation
               cell_surf(auto_ind)            = 5.8e-12_r8    ! surface areas of cells ( m2 )
               cell_carb_cont(auto_ind)       = 3e-14_r8      ! cell carbon content ( mol C cell-1 )
               cell_radius(auto_ind)          = 0.68_r8       ! cell radius ( um )
               cell_permea(auto_ind)          = 3.0e-8_r8     ! cell wall permeability to CO2(aq) (m/s)
               cell_eps_fix(auto_ind)         = 30.0_r8       ! fractionation effect of carbon fixation
            !else if (autotrophs(auto_ind)%exp_calcifier) then
            !Currently not set up to separate exp_calcifiers, needs cell_radius value from data
               !----------------------------------------------------------------------------------------
               ! Calcifier based on P. glacialis ( Keller and morel, 1999; Popp et al., 1998 )
               ! Popp et al express cell carbon content in ( pg C cell-1 ), here we convert in (mol C cell-1)
               ! convert pgC to molC : ! Mc = 12 g / mol ! Mc = 12 e12 pg / mol
               !----------------------------------------------------------------------------------------
               !   cell_active_C_uptake(auto_ind) = 0.0_r9       ! ratio of active carbon uptake to carbon fixation
               !   cell_surf(auto_ind)            = 3886.0_r8    ! surface areas of cells ( m2 )
               !   cell_carb_cont(auto_ind)       = 1.68e-10_r8  ! cell carbon content ( mol C cell-1 )
               !   cell_radius(auto_ind)          =        ! cell radius ( um )
               !   cell_permea(auto_ind)          = 1.1e-5_r8     ! cell wall permeability to CO2(aq) (m/s)
               !   cell_eps_fix(auto_ind)         = 23.0_r8       ! fractionation effect of carbon fixation

            else if (autotrophs(auto_ind)%Nfixer .and. autotrophs(auto_ind)%kSiO3 > c0) then
               call exit_POP(sigAbort, 'ciso: Currently Keller and Morel fractionation does not work for Diatoms-Diazotrophs')

            else
               !----------------------------------------------------------------------------------------
               ! Small phytoplankton based on E. huxleyi ( Keller and morel, 1999; Popp et al., 1998 )
               ! Popp et al express cell carbon content in ( pg C cell-1 ), here we convert in (mol C cell-1)
               ! convert pgC to molC : ! Mc = 12 g / mol ! Mc = 12 e12 pg / mol
               !----------------------------------------------------------------------------------------
               cell_active_C_uptake(auto_ind) = 2.2_r8      ! ratio of active carbon uptake to carbon fixation
               cell_surf(auto_ind)            = 87.6e-12_r8 ! surface areas of cells ( m2 )
               cell_carb_cont(auto_ind)       = 69.2e-14_r8 ! cell carbon content ( mol C cell-1 )
               cell_radius(auto_ind)          = 2.6_r8      ! cell radius ( um )
               cell_permea(auto_ind)          = 1.8e-5_r8   ! cell wall permeability to CO2(aq) (m/s)
               cell_eps_fix(auto_ind)         = 25.3_r8         ! fractionation effect of carbon fixation

         endif
      end do
   end select
!-----------------------------------------------------------------------
!  create local copies of model tracers
!  treat negative values as zero
!  apply mask to local copies
!-----------------------------------------------------------------------

   TEMP           = p5*(TEMP_OLD + TEMP_CUR)

   DI13C_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,di13c_ind) + &
                                TRACER_MODULE_CUR(:,:,k,di13c_ind)))

   DO13C_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,do13c_ind) + &
                              TRACER_MODULE_CUR(:,:,k,do13c_ind)))

   zoo13C_loc     = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,zoo13C_ind) + &
                              TRACER_MODULE_CUR(:,:,k,zoo13C_ind)))

   DI14C_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,di14c_ind) + &
                                TRACER_MODULE_CUR(:,:,k,di14c_ind)))

   DO14C_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,do14c_ind) + &
                              TRACER_MODULE_CUR(:,:,k,do14c_ind)))

   zoo14C_loc     = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,zoo14C_ind) + &
                              TRACER_MODULE_CUR(:,:,k,zoo14C_ind)))

   where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid))
      DIC_loc      = c0
      DI13C_loc    = c0
      DO13C_loc    = c0
      zoo13C_loc   = c0
      DI14C_loc    = c0
      DO14C_loc    = c0
      zoo14C_loc   = c0
   endwhere

   do auto_ind = 1, autotroph_cnt

      n = autotrophs(auto_ind)%C13_ind
      m = autotrophs(auto_ind)%C14_ind

      autotroph13C_loc(:,:,auto_ind) = max(c0, &
         p5*(TRACER_MODULE_OLD(:,:,k,n) + TRACER_MODULE_CUR(:,:,k,n)))

      autotroph14C_loc(:,:,auto_ind) = max(c0, &
         p5*(TRACER_MODULE_OLD(:,:,k,m) + TRACER_MODULE_CUR(:,:,k,m)))


      n = autotrophs(auto_ind)%Ca13CO3_ind
      m = autotrophs(auto_ind)%Ca14CO3_ind
      if (n > 0) then
         autotrophCa13CO3_loc(:,:,auto_ind) = max(c0, &
            p5*(TRACER_MODULE_OLD(:,:,k,n) + TRACER_MODULE_CUR(:,:,k,n)))
      else
         autotrophCa13CO3_loc(:,:,auto_ind) = c0
      endif

      if (m > 0) then
         autotrophCa14CO3_loc(:,:,auto_ind) = max(c0, &
            p5*(TRACER_MODULE_OLD(:,:,k,m) + TRACER_MODULE_CUR(:,:,k,m)))
      else
         autotrophCa14CO3_loc(:,:,auto_ind) = c0
      endif

      where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid))
         autotroph13C_loc(:,:,auto_ind)     = c0
         autotroph14C_loc(:,:,auto_ind)     = c0
         autotrophCa13CO3_loc(:,:,auto_ind) = c0
         autotrophCa14CO3_loc(:,:,auto_ind) = c0
      endwhere

   end do

!-----------------------------------------------------------------------
!  If any ecosys phyto box is zero, set others to zeros
!  (ZERO_MASK in ecosys_mod is equal to ZERO_MASK here)
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt

      ZERO_MASK = autotrophChl_loc(:,:,auto_ind) == c0 .or. &
                  autotrophC_loc(:,:,auto_ind)   == c0 .or. &
                  autotrophFe_loc(:,:,auto_ind)  == c0
      ! Zero_mask=true for any zero in Chl, C, Fw, it is false only if all are false
      ! Add si to zero mask... if it is present (ind > 0)
      if (autotrophs(auto_ind)%Si_ind > 0) &
         ZERO_MASK = ZERO_MASK .or. autotrophSi_loc(:,:,auto_ind) == c0

      where (ZERO_MASK)
         autotroph13C_loc(:,:,auto_ind) = c0
         autotroph14C_loc(:,:,auto_ind) = c0
      endwhere
      if (autotrophs(auto_ind)%Ca13CO3_ind > 0) &
         where (ZERO_MASK) autotrophCa13CO3_loc(:,:,auto_ind) = c0
      if (autotrophs(auto_ind)%Ca14CO3_ind > 0) &
         where (ZERO_MASK) autotrophCa14CO3_loc(:,:,auto_ind) = c0

   end do

!-----------------------------------------------------------------------
!  set local 13C/12C ratios, assuming ecosystem carries 12C (C=C12+C13+C14)
!  If any Carbon boxes are zero, set corresponding 13C to zeros.
!-----------------------------------------------------------------------

   where (DOC_loc > c0)
      R13C_DOC = DO13C_loc / DOC_loc
      R14C_DOC = DO14C_loc / DOC_loc
   elsewhere
      R13C_DOC = c0
      R14C_DOC = c0
   endwhere

   where (DIC_loc > c0)
      R13C_DIC = DI13C_loc / DIC_loc
      R14C_DIC = DI14C_loc / DIC_loc
   elsewhere
      R13C_DIC = c0
      R14C_DIC = c0
   endwhere

   WORK1 = sum(zooC_loc,dim=3)
   where (WORK1 > c0)
      R13C_zooC = zoo13C_loc / WORK1
      R14C_zooC = zoo14C_loc / WORK1
   elsewhere
      R13C_zooC = c0
      R14C_zooC = c0
   endwhere


   do auto_ind = 1, autotroph_cnt

      where (autotrophC_loc(:,:,auto_ind) > c0)
         R13C_autotroph(:,:,auto_ind)  = autotroph13C_loc(:,:,auto_ind) / &
                                         autotrophC_loc(:,:,auto_ind)
         R14C_autotroph(:,:,auto_ind)  = autotroph14C_loc(:,:,auto_ind) / &
                                         autotrophC_loc(:,:,auto_ind)
      elsewhere
         R13C_autotroph(:,:,auto_ind)  = c0
         R14C_autotroph(:,:,auto_ind)  = c0
      endwhere


      where (autotrophCaCO3_loc(:,:,auto_ind) > c0)
         R13C_autotrophCaCO3(:,:,auto_ind) = autotrophCa13CO3_loc(:,:,auto_ind) / &
                                             autotrophCaCO3_loc(:,:,auto_ind)
         R14C_autotrophCaCO3(:,:,auto_ind) = autotrophCa14CO3_loc(:,:,auto_ind) / &
                                             autotrophCaCO3_loc(:,:,auto_ind)
      elsewhere
         R13C_autotrophCaCO3(:,:,auto_ind) = c0
         R14C_autotrophCaCO3(:,:,auto_ind) = c0
      endwhere

   end do
!-----------------------------------------------------------------------
!  Initialize Particulate terms for k=1
!-----------------------------------------------------------------------

   if (k == 1) then
      call ciso_init_particulate_terms(PO13C, P_Ca13CO3, this_block)
      call ciso_init_particulate_terms(PO14C, P_Ca14CO3, this_block)
   endif
!-----------------------------------------------------------------------
! Calculate fraction of CO3
!-----------------------------------------------------------------------

   where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid))
      frac_co3 = c0
   elsewhere
      frac_co3 = CO3 / DIC_loc
   end where

!-----------------------------------------------------------------------
!   discrimination factors of carbone chemistry based on
!   Zhang et al, 1995, Geochim. et Cosmochim. Acta, 59 (1), 107-114
!
!   eps = permil fractionation and alpha is the fractionation factor
!   with eps =(alpha - 1) *1000
!
!   Fractionation is twice as large for 14C compared to 13C
!-----------------------------------------------------------------------

   eps_aq_g   = 0.0049_r8 * TEMP - 1.31_r8
   eps_dic_g  = 0.014_r8 * TEMP * frac_co3 - 0.105_r8 * TEMP + 10.53_r8

   alpha_aq_g  = c1 + eps_aq_g  / c1000
   alpha_dic_g = c1 + eps_dic_g / c1000

!fractionation is twice as large for 14C compared to 13C
   alpha_aq_g_14c  = c1 + eps_aq_g * 2.0_r8  / c1000
   alpha_dic_g_14c = c1 + eps_dic_g * 2.0_r8 / c1000

!-----------------------------------------------------------------------
!  13C/12C ratios of CO2* (CO2STAR)
!-----------------------------------------------------------------------

   R13C_CO2STAR = R13C_DIC * alpha_aq_g / alpha_dic_g

!-----------------------------------------------------------------------
!  delta_13C of CO2* (CO2STAR)
!-----------------------------------------------------------------------

   delta_C13_CO2STAR = ( R13C_CO2STAR / R13C_std - c1 ) * c1000

!-----------------------------------------------------------------------
!  14C/12C ratios of CO2* (CO2STAR)
!-----------------------------------------------------------------------

   R14C_CO2STAR = R14C_DIC * alpha_aq_g_14c / alpha_dic_g_14c

!-----------------------------------------------------------------------
!  delta_14C of CO2* (CO2STAR)
!-----------------------------------------------------------------------

   delta_C14_CO2STAR = ( R14C_CO2STAR / R14C_std - c1 ) * c1000

!-----------------------------------------------------------------------
!  [CO2STAR]  = [CO2*] = [CO2(aq)] + [H2CO3]
!  (this is eq 1.1.1 in Zeebe and Wolf-Gladrow, CO2 in seawater:
!  equilibrium, kinetics, isotopes, Elseview Oceanography Series 65)
!
!  DIC= [CO3] + [HCO3] + [CO2*] (eq 1.1.7)
!
!  => CO2STAR_int = DIC_loc - HCO3 - CO3 !
!-----------------------------------------------------------------------

   CO2STAR_int = DIC_loc - HCO3 - CO3

!------------------------------------------------------------------------
!  Loop over autotrophe types sp, diat, diaz and calculate fractionation
!  for each type
!------------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt

!------------------------------------------------------------------------
!   mu(i) / [CO2*]  of small phytoplankton ( m^3 / mmol C /s )
!-----------------------------------------------------------------------

      where ( CO2STAR_int /= c0 )
         mui_to_co2star(:,:,auto_ind) =  PCphoto(:,:,auto_ind) / CO2STAR_int
      elsewhere
         mui_to_co2star(:,:,auto_ind) = c0
      endwhere

!-----------------------------------------------------------------------
!  fractionation factors for 13C fixation  against CO2* in
!  authotrophe types (sp, diaz, diat)
!-----------------------------------------------------------------------
      select case (ciso_fract_factors)

!-----------------------------------------------------------------------
!   Rau et al., 1989 ( see Gruber et al., 1998 )
!   with restriction between -19 and -31 permil (see Rau et al., 1989)
!-----------------------------------------------------------------------
      case ('Rau')
         delta_C13_Corg = -0.8_r8 * CO2STAR_int - 12.6_r8

         delta_C13_Corg = min( delta_C13_Corg , -18.0_r8 )
         delta_C13_Corg = max( delta_C13_Corg , -32.0_r8 )


         eps_autotroph(:,:,auto_ind) = c1000 * (delta_C13_CO2STAR - delta_C13_Corg) &
                                      /(c1000 + delta_C13_Corg)


!-----------------------------------------------------------------------
!   Laws et al, 1995
!   with restriction between 10 and 26 for size effect (Tagliabue and Bopp, 2008)
!   convert mui_to_co2star from m3/mmol/s to kg/mumol/d
!-----------------------------------------------------------------------
      case ('Laws')
         eps_autotroph(:,:,auto_ind) = (( mui_to_co2star(:,:,auto_ind) * &
                                       seconds_in_day  ) - 0.371_r8 )    &
                                       / (-0.015_r8)

!--------------------------------------------------------------------------
! uncomment the following two lines to restrict eps_sp  between 10 and 26
!--------------------------------------------------------------------------

!         eps_autotroph(:,:,auto_ind) = min( eps_autotroph(:,:,auto_ind), 26.0_r8 )
!         eps_autotroph(:,:,auto_ind) = max( eps_autotroph(:,:,auto_ind), 10.0_r8 )


!-----------------------------------------------------------------------
!   Keller and morel, 1999
!-----------------------------------------------------------------------
      case ('KellerMorel')
!-----------------------------------------------------------------
! convert mui_to_co2start from m3/mmol/s to m3/mol/s
!-----------------------------------------------------------------

         mui_to_co2star_loc = mui_to_co2star(:,:,auto_ind) * c1000


         call fract_keller_morel( LAND_MASK(:,:,bid),        &
                 mui_to_co2star_loc, &
                 cell_active_C_uptake(auto_ind),   &
                 cell_surf(auto_ind),              &
                 cell_carb_cont(auto_ind),         &
                 cell_radius(auto_ind),            &
                 cell_permea(auto_ind),            &
                 cell_eps_fix(auto_ind),           &
                 eps_autotroph(:,:,auto_ind) )

!-----------------------------------------------------------------------
      end select
!-----------------------------------------------------------------------

      where (eps_autotroph(:,:,auto_ind) /= -c1000 )
         R13C_photoC(:,:,auto_ind) = R13C_CO2STAR *c1000 / &
                                    (eps_autotroph(:,:,auto_ind) + c1000)
         R14C_photoC(:,:,auto_ind) = R14C_CO2STAR *c1000 / &
                                    (2.0_r8* eps_autotroph(:,:,auto_ind) + c1000)
      elsewhere
         R13C_photoC(:,:,auto_ind) = c1
         R14C_photoC(:,:,auto_ind) = c1
      endwhere


!-----------------------------------------------------------------------
!     Use R13/14C_photoC to determine small phytoplankton, Diatom, and
!     Diaztroph 13C and 14C fixation
!-----------------------------------------------------------------------

      photo13C(:,:,auto_ind) = photoC(:,:,auto_ind) * R13C_photoC(:,:,auto_ind)


      photo14C(:,:,auto_ind) = photoC(:,:,auto_ind) * R14C_photoC(:,:,auto_ind)

!-----------------------------------------------------------------------
!     C13 & C14 CaCO3 production
!-----------------------------------------------------------------------
      if (autotrophs(auto_ind)%imp_calcifier) then

         R13C_CaCO3_PROD = R13C_DIC + R13C_std * eps_carb / c1000 

         R14C_CaCO3_PROD = R14C_DIC + R14C_std * eps_carb * 2.0_r8 / c1000

         Ca13CO3_PROD(:,:,auto_ind) = CaCO3_PROD(:,:,auto_ind) * R13C_CaCO3_PROD

         Ca14CO3_PROD(:,:,auto_ind) = CaCO3_PROD(:,:,auto_ind) * R14C_CaCO3_PROD

         call accumulate_tavg_field(Ca13CO3_PROD(:,:,auto_ind), tavg_CISO_Ca13CO3_form(auto_ind),bid,k)
         call accumulate_tavg_field(Ca14CO3_PROD(:,:,auto_ind), tavg_CISO_Ca14CO3_form(auto_ind),bid,k)

         if (accumulate_tavg_now(tavg_CISO_Ca13CO3_form_zint(auto_ind))) then
            if (partial_bottom_cells) then
               WORK1 = merge(DZT(:,:,k,bid) * Ca13CO3_PROD(:,:,auto_ind), c0,k<=KMT(:,:,bid))
            else
               WORK1 = merge(dz(k) * Ca13CO3_PROD(:,:,auto_ind), c0,k<=KMT(:,:,bid))
            endif
            call accumulate_tavg_field(WORK1, tavg_CISO_Ca13CO3_form_zint(auto_ind),bid,k)
         endif
         if (accumulate_tavg_now(tavg_CISO_Ca14CO3_form_zint(auto_ind))) then
            if (partial_bottom_cells) then
               WORK1 = merge(DZT(:,:,k,bid) * Ca14CO3_PROD(:,:,auto_ind), c0,k<=KMT(:,:,bid))
            else
               WORK1 = merge(dz(k) * Ca14CO3_PROD(:,:,auto_ind), c0,k<=KMT(:,:,bid))
            endif
            call accumulate_tavg_field(WORK1, tavg_CISO_Ca14CO3_form_zint(auto_ind),bid,k)
         endif
      endif

   end do ! end loop over autotroph types

!-----------------------------------------------------------------------
!  compute terms for DO13C and DO14C
!-----------------------------------------------------------------------

   DO13C_prod = sum(zoo_loss_doc, dim=3) *R13C_zooC +  &
                sum( (auto_loss_doc + auto_graze_doc) * R13C_autotroph, dim=3)

   DO13C_remin = DOC_remin * R13C_DOC


   DO14C_prod = sum(zoo_loss_doc,dim=3) *R14C_zooC +  &
                sum( (auto_loss_doc + auto_graze_doc) * R14C_autotroph, dim=3)

   DO14C_remin = DOC_remin * R14C_DOC



!-----------------------------------------------------------------------
!  large detritus 13C and 14C
!-----------------------------------------------------------------------

   PO13C%prod(:,:,bid) = sum( zoo_loss_poc, dim=3 ) * R13C_zooC + &
                        sum( ( (auto_graze_poc + auto_agg + auto_loss_poc) &
                        * R13C_autotroph), dim=3)


   PO14C%prod(:,:,bid) = sum( zoo_loss_poc, dim=3 ) * R14C_zooC + &
                         sum( ( (auto_graze_poc + auto_agg + auto_loss_poc) &
                         * R14C_autotroph), dim=3)


!-----------------------------------------------------------------------
!  large detrital Ca13CO3 and Ca14CCO3
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%CaCO3_ind > 0) then
         P_Ca13CO3%prod(:,:,bid) = P_CaCO3%prod(:,:,bid) * R13C_autotrophCaCO3(:,:,auto_ind)
         P_Ca14CO3%prod(:,:,bid) = P_CaCO3%prod(:,:,bid) * R14C_autotrophCaCO3(:,:,auto_ind)
      endif
   end do


!-----------------------------------------------------------------------
! Compute particulate terms and write to tavg
!-----------------------------------------------------------------------


   call ciso_compute_particulate_terms(k, POC, P_CaCO3, PO13C, P_Ca13CO3, &
                                  O2_loc, NO3_loc, this_block)


   call ciso_compute_particulate_terms(k, POC, P_CaCO3, PO14C, P_Ca14CO3, &
                                  O2_loc, NO3_loc, this_block)


   call ciso_tavg_particulate_terms(k, POC, P_CaCO3, PO13C, P_Ca13CO3, &
                                  PO14C, P_Ca14CO3, this_block)

!-----------------------------------------------------------------------
! Update DTRACER_MODULE for the 7 carbon pools for each Carbon isotope
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  autotroph Carbon (3 carbon pools)
!  autotroph Ca13CO3 and Ca14CO3
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      WORK1 = auto_graze(:,:,auto_ind) + auto_loss(:,:,auto_ind) + auto_agg(:,:,auto_ind)

      n = autotrophs(auto_ind)%C13_ind
      m = autotrophs(auto_ind)%C14_ind

      DTRACER_MODULE(:,:,n) = photo13C(:,:,auto_ind) - WORK1 * R13C_autotroph(:,:,auto_ind)
      DTRACER_MODULE(:,:,m) = photo14C(:,:,auto_ind) - WORK1 * R14C_autotroph(:,:,auto_ind) - &
                              c14_lambda_inv_sec * autotroph14C_loc(:,:,auto_ind)

      n = autotrophs(auto_ind)%Ca13CO3_ind
      if (n > 0) then
         DTRACER_MODULE(:,:,n) = Ca13CO3_PROD(:,:,auto_ind) - QCaCO3(:,:,auto_ind) &
                                 * WORK1 * R13C_autotrophCaCO3(:,:,auto_ind)
      endif
      n = autotrophs(auto_ind)%Ca14CO3_ind
      if (n > 0) then
         DTRACER_MODULE(:,:,n) = Ca14CO3_PROD(:,:,auto_ind) - QCaCO3(:,:,auto_ind) &
                                 * WORK1 * R14C_autotrophCaCO3(:,:,auto_ind)      &
                                 - c14_lambda_inv_sec * autotrophCa14CO3_loc(:,:,auto_ind)
      endif
   end do



!-----------------------------------------------------------------------
!  zoo 13 and 14 Carbon
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,zoo13C_ind) = sum(auto_graze_zoo * R13C_autotroph, dim=3) - &
                                    sum(zoo_loss,dim=3) *R13C_zooC

   DTRACER_MODULE(:,:,zoo14C_ind) = sum(auto_graze_zoo * R14C_autotroph, dim=3) - &
                                    sum(zoo_loss,dim=3) *R14C_zooC -                &
                                    c14_lambda_inv_sec * zoo14C_loc

!-----------------------------------------------------------------------
!  dissolved organic Matter 13C and 14C
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,do13c_ind) = DO13C_prod - DO13C_remin

   DTRACER_MODULE(:,:,do14c_ind) = DO14C_prod - DO14C_remin  -        &
                                   c14_lambda_inv_sec * DO14C_loc


!-----------------------------------------------------------------------
!   dissolved inorganic Carbon 13 and 14
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,di13c_ind) =        &
       sum( (auto_loss_dic+ auto_graze_dic) * R13C_autotroph, dim=3)   &
       - sum(photo13C, dim=3)                                          &
       + DO13C_remin + PO13C%remin(:,:,bid)                            &
       + sum(zoo_loss_dic,dim=3) * R13C_zooC                                      &
       + P_Ca13CO3%remin(:,:,bid)


   DTRACER_MODULE(:,:,di14c_ind) =                                     &
      sum( (auto_loss_dic+ auto_graze_dic) * R14C_autotroph, dim=3)    &
      - sum(photo14C, dim=3)                                           &
      + DO14C_remin + PO14C%remin(:,:,bid)                             &
      + sum(zoo_loss_dic,dim=3) * R14C_zooC                                       &
      + P_Ca14CO3%remin(:,:,bid)                                       &
      - c14_lambda_inv_sec * DI14C_loc

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%Ca14CO3_ind > 0) then
         DTRACER_MODULE(:,:,di14c_ind) = DTRACER_MODULE(:,:,di14c_ind)     &
                + f_graze_CaCO3_REMIN * auto_graze(:,:,auto_ind)           &
                * QCaCO3(:,:,auto_ind) * R14C_autotrophCaCO3(:,:,auto_ind) &
                - Ca14CO3_PROD(:,:,auto_ind)
     endif

      if (autotrophs(auto_ind)%Ca13CO3_ind > 0) then
         DTRACER_MODULE(:,:,di13c_ind) = DTRACER_MODULE(:,:,di13c_ind)     &
                + f_graze_CaCO3_REMIN * auto_graze(:,:,auto_ind)           &
                * QCaCO3(:,:,auto_ind) * R13C_autotrophCaCO3(:,:,auto_ind) &
                - Ca13CO3_PROD(:,:,auto_ind)
      endif
   end do


!-----------------------------------------------------------------------
!   Calculate oceanic D14C and D13C of carbon pools
!-----------------------------------------------------------------------

   DIC_d13C =  ( R13C_DIC / R13C_std - c1 ) * c1000
   DIC_d14C =  ( R14C_DIC / R14C_std - c1 ) * c1000

   DOC_d13C =  ( R13C_DOC / R13C_std - c1 ) * c1000
   DOC_d14C =  ( R14C_DOC / R14C_std - c1 ) * c1000

   zooC_d13C =  ( R13C_zooC / R13C_std - c1 ) * c1000
   zooC_d14C =  ( R14C_zooC / R14C_std - c1 ) * c1000

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%CaCO3_ind > 0) then
         autotrophCaCO3_d13C(:,:,auto_ind) =  ( R13C_autotrophCaCO3(:,:,auto_ind) / R13C_std - c1 ) * c1000
         autotrophCaCO3_d14C(:,:,auto_ind) =  ( R14C_autotrophCaCO3(:,:,auto_ind) / R14C_std - c1 ) * c1000
      else
         autotrophCaCO3_d13C(:,:,auto_ind) =  c0
         autotrophCaCO3_d14C(:,:,auto_ind) =  c0
      endif

      autotroph_d13C(:,:,auto_ind) =  ( R13C_autotroph(:,:,auto_ind) / R13C_std - c1 ) * c1000
      autotroph_d14C(:,:,auto_ind) =  ( R14C_autotroph(:,:,auto_ind) / R14C_std - c1 ) * c1000
   end do



!-----------------------------------------------------------------------
!  various tavg/history variables
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      call accumulate_tavg_field(autotroph_d13C(:,:,auto_ind), tavg_CISO_d13C(auto_ind),bid,k)
      call accumulate_tavg_field(autotroph_d14C(:,:,auto_ind), tavg_CISO_d14C(auto_ind),bid,k)
      call accumulate_tavg_field(autotrophCaCO3_d14C(:,:,auto_ind), tavg_CISO_autotrophCaCO3_d14C(auto_ind),bid,k)
      call accumulate_tavg_field(autotrophCaCO3_d13C(:,:,auto_ind), tavg_CISO_autotrophCaCO3_d13C(auto_ind),bid,k)
   end do

   call accumulate_tavg_field(DIC_d13C, tavg_CISO_DIC_d13C,bid,k)
   call accumulate_tavg_field(DOC_d13C, tavg_CISO_DOC_d13C,bid,k)
   call accumulate_tavg_field(zooC_d13C, tavg_CISO_zooC_d13C,bid,k)

   call accumulate_tavg_field(DIC_d14C, tavg_CISO_DIC_d14C,bid,k)
   call accumulate_tavg_field(DOC_d14C, tavg_CISO_DOC_d14C,bid,k)
   call accumulate_tavg_field(zooC_d14C, tavg_CISO_zooC_d14C,bid,k)


   if (accumulate_tavg_now(tavg_CISO_photo13C_TOT)) then
      WORK1 = sum(photo13C, dim=3)
      call accumulate_tavg_field(WORK1, tavg_CISO_photo13C_TOT,bid,k)
   endif

   if (accumulate_tavg_now(tavg_CISO_photo14C_TOT)) then
      WORK1 = sum(photo14C, dim=3)
      call accumulate_tavg_field(WORK1, tavg_CISO_photo14C_TOT,bid,k)
   endif

   if (accumulate_tavg_now(tavg_CISO_photo13C_TOT_zint)) then
      if (partial_bottom_cells) then
         WORK1 = merge(DZT(:,:,k,bid) * sum(photo13C, dim=3), c0,k<=KMT(:,:,bid))
      else
         WORK1 = merge(dz(k) * sum(photo13C, dim=3), c0,k<=KMT(:,:,bid))
      endif
      call accumulate_tavg_field(WORK1, tavg_CISO_photo13C_TOT_zint,bid,k)
   endif

    if (accumulate_tavg_now(tavg_CISO_photo14C_TOT_zint)) then
      if (partial_bottom_cells) then
         WORK1 = merge(DZT(:,:,k,bid) * sum(photo14C, dim=3), c0,k<=KMT(:,:,bid))
      else
         WORK1 = merge(dz(k) * sum(photo14C, dim=3), c0,k<=KMT(:,:,bid))
      endif
      call accumulate_tavg_field(WORK1, tavg_CISO_photo14C_TOT_zint,bid,k)
   endif

   do auto_ind = 1, autotroph_cnt
      call accumulate_tavg_field(photo13C(:,:,auto_ind), tavg_CISO_photo13C(auto_ind),bid,k)
      call accumulate_tavg_field(eps_autotroph(:,:,auto_ind), tavg_CISO_eps_autotroph(auto_ind),bid,k)
      call accumulate_tavg_field(mui_to_co2star(:,:,auto_ind), tavg_CISO_mui_to_co2star(auto_ind),bid,k)
      call accumulate_tavg_field(photo14C(:,:,auto_ind), tavg_CISO_photo14C(auto_ind),bid,k)
   end do


   do auto_ind = 1, autotroph_cnt
      if (accumulate_tavg_now(tavg_CISO_photo13C_zint(auto_ind))) then
         if (partial_bottom_cells) then
            WORK1 = merge(DZT(:,:,k,bid) * photo13C(:,:,auto_ind), c0,k<=KMT(:,:,bid))
         else
            WORK1 = merge(dz(k) * photo13C(:,:,auto_ind), c0,k<=KMT(:,:,bid))
         endif
         call accumulate_tavg_field(WORK1, tavg_CISO_photo13C_zint(auto_ind),bid,k)
      endif
   end do

   do auto_ind = 1, autotroph_cnt
      if (accumulate_tavg_now(tavg_CISO_photo14C_zint(auto_ind))) then
         if (partial_bottom_cells) then
            WORK1 = merge(DZT(:,:,k,bid) * photo14C(:,:,auto_ind), c0,k<=KMT(:,:,bid))
         else
            WORK1 = merge(dz(k) * photo14C(:,:,auto_ind), c0,k<=KMT(:,:,bid))
         endif
         call accumulate_tavg_field(WORK1, tavg_CISO_photo14C_zint(auto_ind),bid,k)
      endif
   end do


   call accumulate_tavg_field(DO13C_prod, tavg_CISO_DO13C_prod,bid,k)

   call accumulate_tavg_field(DO14C_prod, tavg_CISO_DO14C_prod,bid,k)

   call accumulate_tavg_field(DO13C_remin, tavg_CISO_DO13C_remin,bid,k)

   call accumulate_tavg_field(DO14C_remin, tavg_CISO_DO14C_remin,bid,k)

   call accumulate_tavg_field(eps_aq_g, tavg_CISO_eps_aq_g,bid,k)

   call accumulate_tavg_field(eps_dic_g, tavg_CISO_eps_dic_g,bid,k)

   ztop = c0
   if (k > 1) ztop = zw(k-1)

   if (accumulate_tavg_now(tavg_CISO_Jint_13Ctot) .or. &
       (accumulate_tavg_now(tavg_CISO_Jint_100m_13Ctot) .and. (ztop < 100.0e2_r8))) then
      WORK1 = DTRACER_MODULE(:,:,di13c_ind) + DTRACER_MODULE(:,:,do13c_ind) &
              + DTRACER_MODULE(:,:,zoo13C_ind)  &
              + sum(DTRACER_MODULE(:,:,autotrophs(:)%C13_ind), dim=3)
      do auto_ind = 1, autotroph_cnt
         n = autotrophs(auto_ind)%Ca13CO3_ind
         if (n > 0) then
            WORK1 = WORK1 + DTRACER_MODULE(:,:,n)
         endif
      end do
      if (accumulate_tavg_now(tavg_CISO_Jint_13Ctot)) then
         if (partial_bottom_cells) then
            WORK2 = merge(DZT(:,:,k,bid) * WORK1, c0, k<=KMT(:,:,bid))
         else
            WORK2 = merge(dz(k) * WORK1, c0, k<=KMT(:,:,bid))
         endif
         ! add back loss to sediments
         WORK2 = WORK2 + merge(PO13C%sed_loss(:,:,bid) + P_Ca13CO3%sed_loss(:,:,bid), &
                               c0, k<=KMT(:,:,bid))
         call accumulate_tavg_field(WORK2,tavg_CISO_Jint_13Ctot,bid,k)
      endif
      if (accumulate_tavg_now(tavg_CISO_Jint_100m_13Ctot) .and. (ztop < 100.0e2_r8)) then
         if (partial_bottom_cells) then
            WORK2 = merge(min(100.0e2_r8 - ztop, DZT(:,:,k,bid)) * WORK1, c0, k<=KMT(:,:,bid))
            ! add back loss to sediments
            WORK2 = WORK2 + merge(PO13C%sed_loss(:,:,bid) + P_Ca13CO3%sed_loss(:,:,bid), &
                                  c0, ztop + DZT(:,:,k,bid) <= 100.0e2_r8 .and. k<=KMT(:,:,bid))
         else
            WORK2 = merge(min(100.0e2_r8 - ztop, dz(k)) * WORK1, c0, k<=KMT(:,:,bid))
            ! add back loss to sediments
            WORK2 = WORK2 + merge(PO13C%sed_loss(:,:,bid) + P_Ca13CO3%sed_loss(:,:,bid), &
                                  c0, ztop + dz(k) <= 100.0e2_r8 .and. k<=KMT(:,:,bid))
         endif
         call accumulate_tavg_field(WORK2,tavg_CISO_Jint_100m_13Ctot,bid,k)
      endif
   endif

   if (accumulate_tavg_now(tavg_CISO_Jint_14Ctot) .or. &
       (accumulate_tavg_now(tavg_CISO_Jint_100m_14Ctot) .and. (ztop < 100.0e2_r8))) then
      WORK1 = DTRACER_MODULE(:,:,di14c_ind) + DTRACER_MODULE(:,:,do14c_ind) &
              + DTRACER_MODULE(:,:,zoo14C_ind)  &
              + sum(DTRACER_MODULE(:,:,autotrophs(:)%C14_ind), dim=3)
      do auto_ind = 1, autotroph_cnt
         n = autotrophs(auto_ind)%Ca14CO3_ind
         if (n > 0) then
            WORK1 = WORK1 + DTRACER_MODULE(:,:,n)
         endif
      end do
      if (accumulate_tavg_now(tavg_CISO_Jint_14Ctot)) then
         if (partial_bottom_cells) then
            WORK2 = merge(DZT(:,:,k,bid) * WORK1, c0, k<=KMT(:,:,bid))
         else
            WORK2 = merge(dz(k) * WORK1, c0, k<=KMT(:,:,bid))
         endif
         ! add back loss to sediments
         WORK2 = WORK2 + merge(PO14C%sed_loss(:,:,bid) + P_Ca14CO3%sed_loss(:,:,bid), &
                               c0, k<=KMT(:,:,bid))
         call accumulate_tavg_field(WORK2,tavg_CISO_Jint_14Ctot,bid,k)
      endif
      if (accumulate_tavg_now(tavg_CISO_Jint_100m_14Ctot) .and. (ztop < 100.0e2_r8)) then
         if (partial_bottom_cells) then
            WORK2 = merge(min(100.0e2_r8 - ztop, DZT(:,:,k,bid)) * WORK1, c0, k<=KMT(:,:,bid))
            ! add back loss to sediments
            WORK2 = WORK2 + merge(PO14C%sed_loss(:,:,bid) + P_Ca14CO3%sed_loss(:,:,bid), &
                                  c0, ztop + DZT(:,:,k,bid) <= 100.0e2_r8 .and. k<=KMT(:,:,bid))
         else
            WORK2 = merge(min(100.0e2_r8 - ztop, dz(k)) * WORK1, c0, k<=KMT(:,:,bid))
            ! add back loss to sediments
            WORK2 = WORK2 + merge(PO14C%sed_loss(:,:,bid) + P_Ca14CO3%sed_loss(:,:,bid), &
                                  c0, ztop + dz(k) <= 100.0e2_r8 .and. k<=KMT(:,:,bid))
         endif
         call accumulate_tavg_field(WORK2,tavg_CISO_Jint_100m_14Ctot,bid,k)
      endif
   endif

!-----------------------------------------------------------------------
 call timer_stop(ecosys_ciso_interior_timer, block_id=bid)


!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_ciso_set_interior
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fract_keller_morel
! !INTERFACE:


 subroutine fract_keller_morel( mask, mui_to_co2star, &
                  cell_active_C_uptake,               &
                  cell_surf,                          &
                  cell_carb_cont,                     &
                  cell_radius,                        &
                  cell_permea,                        &
                  cell_eps_fix,                       &
                  eps_p)

!---------------------------------------------------------------------------
! DESCRIPTION: Calculate the carbon isotopic fractionation of phytoplankton
!              photosynthesis : eps_p
!
! COMPUTATION : based on Keller and Morel, 1999
!
!   eps_p = d13C(co2(aq)) - d13C(phyto)
!
!   eps_p = eps_diff + (cell_active_C_uptake/(cell_active_C_uptake + 1/var)) * delta_d13C
!         + ( (1 + (cell_active_C_uptake-1)*var )/(1+cell_active_C_uptake*var) )
!         * ( cell_eps_fix - eps_diff )
!
!   delta_d13C = d13C(CO2) - d13C(source)
!              = 9 per mil
!
!   var = mui_to_co2star * cell_carb_cont / ( cell_permea * cell_surf )
!
!   mui_to_co2star = mu_i / [CO2*]
!
!  --------
!
!   Developed by X. Giraud, ETH Zürich, 21.07.2008
!   Converted to Fortran 90, adapted to changes since POP1, and included in
!   POP2 ciso code by A. Jahn, NCAR, in 10/2012
!---------------------------------------------------------------------------
! !USES:
!    use domain, only : nx_block, ny_block
!    use kinds_mod
!    use constants, only : c0, c1

! !INPUT PARAMETERS:

    logical (log_kind), dimension (nx_block,ny_block), intent(in) :: mask

    real (r8), dimension (nx_block,ny_block), intent(in) :: &
      mui_to_co2star    ! mui_to_co2star = mu_i / [CO2*] (m3 / mol C / s)

    real (r8), intent(in) :: &
      cell_active_C_uptake,  & ! ratio of active carbon uptake to carbon fixation
      cell_eps_fix,          & ! fractionation effect of carbon fixation
      cell_surf,             & ! surface areas of cells ( m2 )
      cell_carb_cont,        & ! cell carbon content ( mol C )
      cell_radius,           & ! cell radius ( um )
      cell_permea              ! cell wall permeability to CO2(aq) ( m /s )

! !OUTPUT PARAMETERS:

    real (r8), dimension (nx_block,ny_block), intent(out) :: &
      eps_p                 ! = d13C(co2(aq)) - d13C(phyto)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables and parameters
!-----------------------------------------------------------------------

    real (r8), dimension (nx_block,ny_block) :: &
      var, theta, eps_up

    real (r8) :: &
      pi, Vol, Qc, Surf, radius_m


   real (r8) :: &
      eps_diff   = 0.7_r8, &  ! fractionation by diffusion, O'Leary, 1984
      delta_d13C = -9.0_r8    ! = d13C(CO2) - d13C(source), difference between the
                              ! isotopic compositions of the external CO2 and the
                              ! organic matter pools (Goericke et al. 1994).
                              ! For active HCO3- uptake, the substrate for
                              ! the carbon uptake mechanism has an isotopic
                              ! composition which is around 9 permil higher
                              ! than the external CO2, so D13CO2 -D13C_source
                              ! is -9 permil ((Mook et al. 197)

!---------------------------------------------------------------------------
!   check for existence of ocean points
!---------------------------------------------------------------------------

    if (count(mask) == 0) then
       eps_p = c0
       return
    endif

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------

    pi  = 4.0_r8 * atan( 1.0_r8 )

!---------------------------------------------------------------------
!  cell surface in m^2
!---------------------------------------------------------------------

    radius_m = cell_radius * 1e-6_r8 ! convert radius from um to m

    if ( cell_surf > c0 ) then
      Surf = cell_surf
    else
      Surf = 4.0_r8 * pi * (radius_m ** 2)
    endif

!---------------------------------------------------------------------
!     cellular carbon content ( mol C )
!     volume in um^3
!---------------------------------------------------------------------

    if ( cell_carb_cont > c0 ) then
      Qc = cell_carb_cont
    else
      Vol = 4.0_r8 * pi * (cell_radius ** 3) / 3.0_r8
      Qc = 3.154e-14_r8 * (Vol ** (0.758_r8 ))
    endif

!---------------------------------------------------------------------
!     final expression of eps_p
!---------------------------------------------------------------------

    where ( mui_to_co2star /= c0 )

          var = mui_to_co2star * Qc / ( cell_permea * Surf )

          theta = c1 + ( cell_active_C_uptake - c1 ) * var
          theta = theta / ( c1 + cell_active_C_uptake * var )

          eps_up = eps_diff + ( cell_active_C_uptake /  &
                 ( cell_active_C_uptake + c1 / var ) ) * delta_d13C

          eps_p = eps_up + theta * ( cell_eps_fix - eps_diff )



    elsewhere

    eps_p = cell_eps_fix

    endwhere

    where (.not.mask)

    eps_p = c0

    endwhere
!-----------------------------------------------------------------------
!EOC
 end subroutine fract_keller_morel
!-----------------------------------------------------------------------


!***********************************************************************
!BOP
! !IROUTINE: ciso_init_particulate_terms
! !INTERFACE:

 subroutine ciso_init_particulate_terms(POC_ciso, P_CaCO3_ciso, this_block)

! !DESCRIPTION:
!  Set incoming fluxes (put into outgoing flux for first level usage).
!  Set dissolution length, production fraction and mass terms.
!
!  The first 2 arguments are intent(inout) in
!  order to preserve contents on other blocks.

! !INPUT/OUTPUT PARAMETERS:


   type(sinking_particle), intent(inout) :: &
      POC_ciso,        & ! base units = nmol C_ciso
      P_CaCO3_ciso       ! base units = nmol C_ciso


! !INPUT PARAMETERS:


   type (block), intent(in) :: &
      this_block      ! block info for the current block

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bid                 ! local_block id
!-----------------------------------------------------------------------

   POC_ciso%diss      = c0       ! not used
   POC_ciso%gamma     = c0       ! not used
   POC_ciso%mass      = c0       ! not used
   POC_ciso%rho       = c0       ! not used

   P_CaCO3_ciso%diss  = c0       ! not used
   P_CaCO3_ciso%gamma = c0       ! not used
   P_CaCO3_ciso%mass  = c0       ! not used
   P_CaCO3_ciso%rho   = c0       ! not used


!-----------------------------------------------------------------------
!  Set incoming fluxes
!-----------------------------------------------------------------------

   bid = this_block%local_id

   P_CaCO3_ciso%sflux_out(:,:,bid) = c0
   P_CaCO3_ciso%hflux_out(:,:,bid) = c0

!-----------------------------------------------------------------------
!  Hard POC is QA flux and soft POC is excess POC.
!-----------------------------------------------------------------------

   POC_ciso%sflux_out(:,:,bid) = c0
   POC_ciso%hflux_out(:,:,bid) = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine ciso_init_particulate_terms

!***********************************************************************
!BOP
! !IROUTINE: ciso_compute_particulate_terms
! !INTERFACE:
 subroutine ciso_compute_particulate_terms(k, POC, P_CaCO3, POC_ciso, P_CaCO3_ciso, &
              O2_loc, NO3_loc, this_block)

! !DESCRIPTION:
!  Compute outgoing fluxes and remineralization terms for Carbon isotopes.
!  Assumes that production terms have been set and that fluxes and remineralization
!  for Carbon 12 has already been computed.
!
!  Incoming fluxes are assumed to be the outgoing fluxes from the previous level.
!  For other comments, see compute_particulate_terms in ecosys_mod
!
!
!  Alex Jahn, Nov 2012
! !USES:

#ifdef CCSMCOUPLED
   use shr_sys_mod, only: shr_sys_abort
#endif

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k ! vertical model level

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      O2_loc,       & ! dissolved oxygen used to modify POC%diss, Sed fluxes
      NO3_loc         ! dissolved nitrate used to modify sed fluxes

   type (block), intent(in) :: &
      this_block      ! block info for the current block



   type(sinking_particle), intent(in) :: &
      POC,          & ! base units = nmol C
      P_CaCO3         ! base units = nmol CaCO3

! !INPUT/OUTPUT PARAMETERS:

   type(sinking_particle), intent(inout) :: &
      POC_ciso,       &  ! base units = nmol particulate organic Carbon isotope
      P_CaCO3_ciso       ! base units = nmol CaCO3 Carbon isotope




!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'ecosys_ciso_mod:ciso_compute_particulate_terms'

   real (r8) :: &
      dz_loc,              & ! dz at a particular i,j location
      dzr_loc                ! dzr at a particular i,j location

   integer (int_kind) ::   &
      i, j,                & ! loop indices
      bid                    ! local_block id

   logical (log_kind) ::   &
      poc_error              ! POC error flag

   real (r8) :: &
      flux, flux_alt,       & ! temp variables used to update sinking flux
      POC_ciso_PROD_avail,  & ! 13C POC production available for excess POC flux
      Rciso_POC_hflux_out, & ! ciso/12C of outgoing flux of hard POC
      Rciso_POC_in,        & ! ciso/12C of total POC ingoing component
      Rciso_CaCO3_in         ! ciso/12C of total CaCO3 ingoing component

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      SED_DENITRIF, & ! sedimentary denitrification (umolN/cm^2/s)
      OTHER_REMIN     ! sedimentary remin not due to oxic or denitrification

!-----------------------------------------------------------------------
!  Pointer variables (targets are defined in ecosys_share)
!-----------------------------------------------------------------------
   real (r8), dimension(:,:),pointer :: decay_CaCO3           ! scaling factor for dissolution of CaCO3
   real (r8), dimension(:,:),pointer :: DECAY_Hard            ! scaling factor for dissolution of Hard Ballast
   real (r8), dimension(:,:),pointer :: decay_POC_E           ! scaling factor for dissolution of excess POC
   real (r8), dimension(:,:),pointer :: POC_PROD_avail        ! scaling factor for dissolution of Hard Ballast
   real (r8), dimension(:,:),pointer :: poc_diss              ! diss. length used (cm)
   real (r8), dimension(:,:),pointer :: caco3_diss            ! CaCO3 diss. length used (cm)
   real (r8), dimension(:,:),pointer :: P_CaCO3_sflux_out
   real (r8), dimension(:,:),pointer :: P_CaCO3_hflux_out
   real (r8), dimension(:,:),pointer :: POC_sflux_out
   real (r8), dimension(:,:),pointer :: POC_hflux_out
   real (r8), dimension(:,:),pointer :: POC_remin
   real (r8), dimension(:,:),pointer :: P_CaCO3_remin

!-----------------------------------------------------------------------
!  this_block index as integer
!-----------------------------------------------------------------------

   bid = this_block%local_id

!---------------------------------------------------------------
! Need the following variables from ecosys_mod --> use pointers
! to point to the right part of the global array in ecosys_share
!---------------------------------------------------------------

   decay_CaCO3      => decay_CaCO3_fields(:,:,bid)
   DECAY_Hard       => DECAY_Hard_fields(:,:,bid)
   decay_POC_E      => decay_POC_E_fields(:,:,bid)
   POC_PROD_avail   => POC_PROD_avail_fields(:,:,bid)
   poc_diss         => poc_diss_fields(:,:,bid)
   caco3_diss       => caco3_diss_fields(:,:,bid)
   P_CaCO3_sflux_out=> P_CaCO3_sflux_out_fields(:,:,bid)
   P_CaCO3_hflux_out=> P_CaCO3_hflux_out_fields(:,:,bid)
   POC_sflux_out    => POC_sflux_out_fields(:,:,bid)
   POC_hflux_out    => POC_hflux_out_fields(:,:,bid)
   POC_remin        => POC_remin_fields(:,:,bid)
   P_CaCO3_remin    => P_CaCO3_remin_fields(:,:,bid)

!-----------------------------------------------------------------------
!  incoming fluxes are outgoing fluxes from previous level
!-----------------------------------------------------------------------

   POC_ciso%sflux_in(:,:,bid) = POC_ciso%sflux_out(:,:,bid)
   POC_ciso%hflux_in(:,:,bid) = POC_ciso%hflux_out(:,:,bid)

   P_CaCO3_ciso%sflux_in(:,:,bid) = P_CaCO3_ciso%sflux_out(:,:,bid)
   P_CaCO3_ciso%hflux_in(:,:,bid) = P_CaCO3_ciso%hflux_out(:,:,bid)

!-----------------------------------------------------------------------
!  initialize loss to sediments = 0 and local copy of percent sed
!-----------------------------------------------------------------------

   POC_ciso%sed_loss(:,:,bid)     = c0
   P_CaCO3_ciso%sed_loss(:,:,bid) = c0
   SED_DENITRIF(:,:,bid)          = c0
   OTHER_REMIN(:,:,bid)           = c0

!-----------------------------------------------------------------------
!     if any incoming Carbon flux is zero, set 13C flux to zero
!-----------------------------------------------------------------------

   where (POC%sflux_in(:,:,bid) == c0)
      POC_ciso%sflux_in(:,:,bid) = c0
   endwhere

   where (POC%hflux_in(:,:,bid) == c0)
      POC_ciso%hflux_in(:,:,bid) = c0
   endwhere


   where (P_CaCO3%sflux_in(:,:,bid) == c0)
      P_CaCO3_ciso%sflux_in(:,:,bid) = c0
   endwhere

   where (P_CaCO3%hflux_in(:,:,bid) == c0)
      P_CaCO3_ciso%hflux_in(:,:,bid) = c0
   endwhere


   dz_loc = dz(k)

   do j = 1,ny_block
      do i = 1,nx_block

         if (LAND_MASK(i,j,bid) .and. k <= KMT(i,j,bid)) then

            if (partial_bottom_cells) then
               dz_loc = DZT(i,j,k,bid)
            endif
            dzr_loc = c1 / dz_loc

!-----------------------------------------------------------------------
!  P_CaCO3_ciso sflux and hflux out
!-----------------------------------------------------------------------

            P_CaCO3_ciso%sflux_out(i,j,bid) = P_CaCO3_ciso%sflux_in(i,j,bid) * decay_CaCO3(i,j) + &
                P_CaCO3_ciso%prod(i,j,bid) * ((c1 - P_CaCO3%gamma) * (c1 - decay_CaCO3(i,j))  &
                * caco3_diss(i,j))

            P_CaCO3_ciso%hflux_out(i,j,bid) = P_CaCO3_ciso%hflux_in(i,j,bid) * DECAY_Hard(i,j) + &
                 P_CaCO3_ciso%prod(i,j,bid) * (P_CaCO3%gamma * dz_loc)

!-----------------------------------------------------------------
!   Compute how much 13C POC_PROD is available for deficit
!   reduction and excess POC flux
!-----------------------------------------------------------------

            if (POC%prod(i,j,bid) > c0 ) then
                POC_ciso_PROD_avail = POC_PROD_avail(i,j) * POC_ciso%prod(i,j,bid) / POC%prod(i,j,bid)
            else
                POC_ciso_PROD_avail = c0
            endif

!-----------------------------------------------------------------
!   Compute outgoing 13C POC fluxes of soft POC
!-----------------------------------------------------------------

            POC_ciso%sflux_out(i,j,bid) = POC_ciso%sflux_in(i,j,bid) *decay_POC_E(i,j) + &
                POC_ciso_PROD_avail *((c1 - decay_POC_E(i,j)) * (poc_diss(i,j)) )

!-----------------------------------------------------------------
!   Compute outgoing 13C POC fluxes of hard POC
!-----------------------------------------------------------------
            if (POC_ciso%hflux_in(i,j,bid) == c0 .and. POC_ciso%prod(i,j,bid) == c0) then
               POC_ciso%hflux_out(i,j,bid) = c0
            else

               Rciso_POC_hflux_out = POC%prod(i,j,bid) + &
                                     ( POC%sflux_in(i,j,bid) - POC_sflux_out(i,j) + &
                                     POC%hflux_in(i,j,bid) ) * dzr_loc

               if (Rciso_POC_hflux_out  /= c0) then
                  Rciso_POC_hflux_out = ( POC_ciso%prod(i,j,bid) + ( POC_ciso%sflux_in(i,j,bid) - &
                        POC_ciso%sflux_out(i,j,bid) + POC_ciso%hflux_in(i,j,bid) ) * dzr_loc )    &
                        / Rciso_POC_hflux_out
                else
                  Rciso_POC_hflux_out = c0
                endif

                POC_ciso%hflux_out(i,j,bid) = POC_hflux_out(i,j) * Rciso_POC_hflux_out
                POC_ciso%hflux_out(i,j,bid) = max(POC_ciso%hflux_out(i,j,bid), c0)

            endif

!-----------------------------------------------------------------------
!  Compute remineralization terms. It is assumed that there is no
!  sub-surface dust production.
!-----------------------------------------------------------------------

            P_CaCO3_ciso%remin(i,j,bid) = P_CaCO3_ciso%prod(i,j,bid) + &
               ((P_CaCO3_ciso%sflux_in(i,j,bid) - P_CaCO3_ciso%sflux_out(i,j,bid)) + &
               (P_CaCO3_ciso%hflux_in(i,j,bid) - P_CaCO3_ciso%hflux_out(i,j,bid))) *  dzr_loc

            POC_ciso%remin(i,j,bid) = POC_ciso%prod(i,j,bid)  &
                + ((POC_ciso%sflux_in(i,j,bid) - POC_ciso%sflux_out(i,j,bid)) &
                + (POC_ciso%hflux_in(i,j,bid) - POC_ciso%hflux_out(i,j,bid))) &
                * dzr_loc

!-----------------------------------------------------------------
!   Option to force the 13C/12C ratio of outgoing P_CaCO3 fluxes
!   to equal the rate of total incoming flux
!-----------------------------------------------------------------
!
!             Rciso_CaCO3_in = P_CaCO3%prod(i,j,bid) + &
!                 ( P_CaCO3%sflux_in(i,j,bid) + P_CaCO3%hflux_in(i,j,bid) ) *  dzr_loc
!
!             if ( Rciso_CaCO3_in > c0 ) then
!                 Rciso_CaCO3_in = ( P_CaCO3_ciso%prod(i,j,bid) + &
!                   ( P_CaCO3_ciso%sflux_in(i,j,bid) + P_CaCO3_ciso%hflux_in(i,j,bid) ) * dzr_loc) &
!                    / Rciso_CaCO3_in
!             else
!                     Rciso_CaCO3_in = c0
!             endif
!
!              P_CaCO3_ciso%sflux_out(i,j,bid) = P_CaCO3_sflux_out(i,j) * Rciso_CaCO3_in
!              P_CaCO3_ciso%hflux_out(i,j,bid) = P_CaCO3_hflux_out(i,j) * Rciso_CaCO3_in
!              P_CaCO3_ciso%remin(i,j,bid) = P_CaCO3_remin(i,j) * Rciso_CaCO3_in
!
!-----------------------------------------------------------------
!   Option to force the 13C/12C ratio of outgoing POC fluxes
!   to equal the rate of total incoming flux
!-----------------------------------------------------------------
!              Rciso_POC_in = POC%prod(i,j,bid) + ( POC%sflux_in(i,j,bid) + POC%hflux_in(i,j,bid) ) &
!                            * dzr_loc
!              if ( Rciso_POC_in > c0 ) then
!                 Rciso_POC_in = ( POC_ciso%prod(i,j,bid) + ( POC_ciso%sflux_in(i,j,bid) + &
!                    POC_ciso%hflux_in(i,j,bid) ) * dzr_loc ) / Rciso_POC_in
!              else
!                 Rciso_POC_in = c0
!              endif

!              POC_ciso%sflux_out(i,j,bid) = POC_sflux_out(i,j) *Rciso_POC_in
!              POC_ciso%hflux_out(i,j,bid) = POC_hflux_out(i,j) *Rciso_POC_in
!              POC_ciso%remin(i,j,bid) = POC_remin(i,j) * Rciso_POC_in
!-----------------------------------------------------------------

         else
            P_CaCO3_ciso%sflux_out(i,j,bid) = c0
            P_CaCO3_ciso%hflux_out(i,j,bid) = c0
            P_CaCO3_ciso%remin(i,j,bid)     = c0

            POC_ciso%sflux_out(i,j,bid) = c0
            POC_ciso%hflux_out(i,j,bid) = c0
            POC_ciso%remin(i,j,bid)     = c0
         endif


!-----------------------------------------------------------------------
!  Bottom Sediments Cell?
!  If so compute sedimentary burial and denitrification N losses.
!  Using empirical relations from Bohlen et al., 2012 (doi:10.1029/2011GB004198) for Sed Denitrification
!  OTHER_REMIN estimates organic matter remineralized in the sediments
!      by the processes other than oxic remin and denitrification (SO4 and CO2,
!      etc..)
!      based on Soetaert et al., 1996, varies between 10% and 50%
!      0.4_r8 is a coefficient with units mmolC/cm2/yr sinking flux,
!      OTHER_REMIN is 50% above this high flux value,
!      In special case where bottom O2 has been depleted to < 1.0 uM,
!               all sedimentary remin is due to DENITRIFICATION + OTHER_REMIN
!  POC burial from Dunne et al. 2007 (doi:10.1029/2006GB002907), maximum of 80% burial efficiency imposed
!  Bsi preservation in sediments = 0.3*sinkBsi - 0.06 mmol/m2/day
!     Ragueneau et al. 2000 (doi:10.1016/S0921-8181(00)00052-7)
!  Calcite is preserved in sediments above the lysocline, dissolves below.
!       Here a constant depth is used for lysocline.
!-----------------------------------------------------------------------

         if (LAND_MASK(i,j,bid) .and. (k == KMT(i,j,bid))) then

            flux = POC_ciso%sflux_out(i,j,bid)+POC_ciso%hflux_out(i,j,bid)

            if (flux > c0) then
               flux_alt = flux*mpercm*spd ! convert to mmol/m^2/day

               POC_ciso%sed_loss(i,j,bid) = flux * min(0.8_r8, parm_POMbury &
                  * (0.013_r8 + 0.53_r8 * flux_alt*flux_alt / (7.0_r8 + flux_alt)**2))

               SED_DENITRIF(i,j,bid) = dzr_loc * flux &
                  * (0.06_r8 + 0.19_r8 * 0.99_r8**(O2_loc(i,j)-NO3_loc(i,j)))

               flux_alt = flux*1.0e-6_r8*spd*365.0_r8 ! convert to mmol/cm^2/year
               OTHER_REMIN(i,j,bid) = dzr_loc &
                  * min(min(0.1_r8 + flux_alt,0.5_r8) * (flux - POC_ciso%sed_loss(i,j,bid)), &
                        (flux - POC_ciso%sed_loss(i,j,bid) - (SED_DENITRIF(i,j,bid)*dz_loc*denitrif_C_N)))

!----------------------------------------------------------------------------------
!              if bottom water O2 is depleted, assume all remin is denitrif + other
!----------------------------------------------------------------------------------

               if (O2_loc(i,j) < c1) then
                  OTHER_REMIN(i,j,bid) = dzr_loc * &
                                        (flux - POC_ciso%sed_loss(i,j,bid) - &
                                        (SED_DENITRIF(i,j,bid)*dz_loc*denitrif_C_N))
               endif

            endif

            if (zw(k) < 3300.0e2_r8) then
               flux = P_CaCO3_ciso%sflux_out(i,j,bid) + P_CaCO3_ciso%hflux_out(i,j,bid)
               P_CaCO3_ciso%sed_loss(i,j,bid) = flux
            endif

!----------------------------------------------------------------------------------
!  Update sinking fluxes and remin fluxes, accounting for sediments.
!  flux used to hold sinking fluxes before update.
!----------------------------------------------------------------------------------

            flux = P_CaCO3_ciso%sflux_out(i,j,bid) + P_CaCO3_ciso%hflux_out(i,j,bid)
            if (flux > c0) then
               P_CaCO3_ciso%remin(i,j,bid) = P_CaCO3_ciso%remin(i,j,bid) &
                     + ((flux - P_CaCO3_ciso%sed_loss(i,j,bid)) * dzr_loc)
            endif

            flux = POC_ciso%sflux_out(i,j,bid) + POC_ciso%hflux_out(i,j,bid)
            if (flux > c0) then
               POC_ciso%remin(i,j,bid) = POC_ciso%remin(i,j,bid) &
                     + ((flux - POC_ciso%sed_loss(i,j,bid)) * dzr_loc)
            endif

!-----------------------------------------------------------------------
!   Set all outgoing fluxes to 0.0
!-----------------------------------------------------------------------

            if (k == KMT(i,j,bid)) then
               P_CaCO3_ciso%sflux_out(i,j,bid) = c0
               P_CaCO3_ciso%hflux_out(i,j,bid) = c0

               POC_ciso%sflux_out(i,j,bid) = c0
               POC_ciso%hflux_out(i,j,bid) = c0

            endif

         endif
      end do
   end do





!-----------------------------------------------------------------------
!EOC

end subroutine ciso_compute_particulate_terms

!***********************************************************************
!BOP
! !IROUTINE: ciso_tavg_particulate_terms
! !INTERFACE:
 subroutine ciso_tavg_particulate_terms(k, POC, P_CaCO3, PO13C, P_Ca13CO3, &
               PO14C, P_Ca14CO3, this_block)

! !DESCRIPTION:
!  Writes tavg for particulate terms calculated in ciso_compute_particulate_terms
!
!
!  Alex Jahn, Nov 2012
! !USES:

#ifdef CCSMCOUPLED
   use shr_sys_mod, only: shr_sys_abort
#endif

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k ! vertical model level


   type (block), intent(in) :: &
      this_block      ! block info for the current block


   type(sinking_particle), intent(in) :: &
      POC,          &  ! base units = nmol C
      P_CaCO3,      &  ! base units = nmol CaCO3
      PO13C,        &  ! base units = nmol 13C
      P_Ca13CO3,    &  ! base units = nmol 13C CaCO3
      PO14C,        &  ! base units = nmol 14C
      P_Ca14CO3        ! base units = nmol 14C CaCO3





!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


   character(*), parameter :: &
      subname = 'ecosys_ciso_mod:ciso_tavg_particulate_terms'

   real (r8), dimension(nx_block,ny_block) :: &
      WORK                  ! temporary for summed quantities to be averaged

   integer (int_kind) :: &
      bid                   ! local_block id

!-----------------------------------------------------------------------
!  this_block index as integer
!-----------------------------------------------------------------------

   bid = this_block%local_id

!-----------------------------------------------------------------------
!  Set tavg variables.
!-----------------------------------------------------------------------

   if (accumulate_tavg_now(tavg_CISO_PO13C_FLUX_IN)) then
      WORK = PO13C%sflux_in(:,:,bid) + PO13C%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_CISO_PO13C_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(PO13C%prod(:,:,bid), tavg_CISO_PO13C_PROD,bid,k)

   call accumulate_tavg_field(PO13C%remin(:,:,bid), tavg_CISO_PO13C_REMIN,bid,k)

   if (accumulate_tavg_now(tavg_CISO_Ca13CO3_FLUX_IN)) then
      WORK = P_Ca13CO3%sflux_in(:,:,bid) + P_Ca13CO3%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_CISO_Ca13CO3_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(P_Ca13CO3%prod(:,:,bid), tavg_CISO_Ca13CO3_PROD,bid,k)

   call accumulate_tavg_field(P_Ca13CO3%remin(:,:,bid), tavg_CISO_Ca13CO3_REMIN,bid,k)


!14C

   if (accumulate_tavg_now(tavg_CISO_PO14C_FLUX_IN)) then
      WORK = PO14C%sflux_in(:,:,bid) + PO14C%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_CISO_PO14C_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(PO14C%prod(:,:,bid), tavg_CISO_PO14C_PROD,bid,k)

   call accumulate_tavg_field(PO14C%remin(:,:,bid), tavg_CISO_PO14C_REMIN,bid,k)

   if (accumulate_tavg_now(tavg_CISO_Ca14CO3_FLUX_IN)) then
      WORK = P_Ca14CO3%sflux_in(:,:,bid) + P_Ca14CO3%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_CISO_Ca14CO3_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(P_Ca14CO3%prod(:,:,bid), tavg_CISO_Ca14CO3_PROD,bid,k)

   call accumulate_tavg_field(P_Ca14CO3%remin(:,:,bid), tavg_CISO_Ca14CO3_REMIN,bid,k)

! ***********************************************************************
! - Accumulte losses of BGC tracers to sediments
! ***********************************************************************


   call accumulate_tavg_field(P_Ca13CO3%sed_loss(:,:,bid), tavg_calcToSed_13C,bid,k)

   call accumulate_tavg_field(PO13C%sed_loss(:,:,bid), tavg_pocToSed_13C,bid,k)


   call accumulate_tavg_field(P_Ca14CO3%sed_loss(:,:,bid), tavg_calcToSed_14C,bid,k)

   call accumulate_tavg_field(PO14C%sed_loss(:,:,bid), tavg_pocToSed_14C,bid,k)

!-----------------------------------------------------------------------
!EOC

 end subroutine ciso_tavg_particulate_terms

!***********************************************************************
!BOP
! !IROUTINE: ciso_read_atm_D13C_data
! !INTERFACE:

 subroutine ciso_read_atm_D13C_data

! !DESCRIPTION:
!  Read atmospheric D13C [permil] data from file
!
!  Have the master_task do the following :
!     1) get length of data
!     2) allocate memory for data
!     3) read in data, checking for consistent lengths
!  Then, outside master_task conditional
!     1) broadcast length of data
!     2) have non-mastertasks allocate memory for data
!     3) broadcast data
!
! !REVISION HISTORY:
!  same as module
!
! !USES:
!
!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   character(*), parameter :: sub_name = 'ecosys_ciso_mod:ciso_read_atm_D13C_data'

   integer (int_kind) ::    &
      stat,                 &  ! i/o status code
      irec,                 &  ! counter for looping
      skiplines,            &  ! number of comment lines at beginning of ascii file
      il                       ! looping index

   character (char_len) :: &
      sglchr                   ! variable to read characters from file into


!-----------------------------------------------------------------------
!     READ in D13C data from file
!-----------------------------------------------------------------------
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*)'ciso: Using varying D13C values from file ',trim(ciso_atm_d13c_filename)
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      open (nml_in, file=ciso_atm_d13c_filename, status='old',iostat=stat)
      if (stat /= 0) then
         write(stdout,fmt=*) 'open failed'
         go to 99
      endif
      read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d13c_data_nbval
      if (stat /= 0) then
         write(stdout,fmt=*) '1st line read failed'
         go to 99
      endif
      allocate(ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval))
      allocate(ciso_atm_d13c_data(ciso_atm_d13c_data_nbval))
      do irec=1,skiplines
         read(nml_in,FMT=*,iostat=stat) sglchr
         if (stat /= 0) then
            write(stdout,fmt=*) 'skipline read failed'
            go to 99
         endif
      enddo
      do irec=1,ciso_atm_d13c_data_nbval
         read(nml_in,FMT=*,iostat=stat) ciso_atm_d13c_data_yr(irec), ciso_atm_d13c_data(irec)
         if (stat /= 0) then
            write(stdout,fmt=*) 'data read failed'
            go to 99
         endif
      enddo
      close(nml_in)
   endif

99 call broadcast_scalar(stat, master_task)
   if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
                                                          &/ sub_name)

!---------------------------------------------------------------------
!     Need to allocate and broadcast the variables to other tasks beside master-task
!---------------------------------------------------------------------

   call broadcast_scalar(ciso_atm_d13c_data_nbval,master_task)

   if (my_task /= master_task) then
      allocate(ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval))
      allocate(ciso_atm_d13c_data(ciso_atm_d13c_data_nbval))
   endif


   call broadcast_array(ciso_atm_d13c_data, master_task)
   call broadcast_array(ciso_atm_d13c_data_yr, master_task)


!-----------------------------------------------------------------------
!EOC

 end subroutine ciso_read_atm_D13C_data

!***********************************************************************
!BOP
! !IROUTINE: ciso_read_atm_D14C_data
! !INTERFACE:

 subroutine ciso_read_atm_D14C_data

! !DESCRIPTION:
!  Read atmospheric D14C data from file
!
!  Have the master_task do the following :
!     1) get length of data
!     2) allocate memory for data
!     3) read in data, checking for consistent lengths
!  Then, outside master_task conditional
!     1) broadcast length of data
!     2) have non-mastertasks allocate memory for data
!     3) broadcast data
!
! !REVISION HISTORY:
!  same as module
!
! !USES:
!
!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   character(*), parameter :: sub_name = 'ciso_read_atm_D14C_data:ciso_read_atm_D14C_data'

   integer (int_kind) ::      &
      stat,                   &  ! i/o status code
      irec,                   &  ! counter for looping
      skiplines,              &  ! number of comment lines at beginning of ascii file
      il                         ! looping index

   character (char_len) ::  &
      sglchr                     ! variable to read characters from file into

   integer (int_kind) :: &
      ciso_atm_d14c_data_nbval_tmp

   logical (log_kind) :: &
      nbval_mismatch

!-----------------------------------------------------------------------
!     ensure that three datafiles have same number of entries
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,*)'ciso DIC14 calculation: Using varying C14 values from files'
      do il=1,3
         write(stdout,*) trim(ciso_atm_d14c_filename(il))
      enddo
      nbval_mismatch = .false.
      do il=1,3
         open (nml_in,file=ciso_atm_d14c_filename(il),status='old',iostat=stat)
         if (stat /= 0) then
            write(stdout,*) 'open failed for ', trim(ciso_atm_d14c_filename(il))
            go to 99
         endif
         read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d14c_data_nbval_tmp
         if (stat /= 0) then
            write(stdout,*) '1st line read failed for ', trim(ciso_atm_d14c_filename(il))
            go to 99
         endif
         close(nml_in)
         if (il == 1) then
            ciso_atm_d14c_data_nbval = ciso_atm_d14c_data_nbval_tmp
         else
            if (ciso_atm_d14c_data_nbval /= ciso_atm_d14c_data_nbval_tmp) nbval_mismatch = .true.
         endif
      enddo
   endif

   call broadcast_scalar(nbval_mismatch, master_task)
   if (nbval_mismatch) then
      call document(sub_name, 'D14C data files must all have the same number of values')
      call exit_POP(sigAbort, 'stopping in ' /&
                             &/ sub_name)
   endif

   call broadcast_scalar(ciso_atm_d14c_data_nbval, master_task)
   allocate(ciso_atm_d14c_data_yr(ciso_atm_d14c_data_nbval,3))
   allocate(ciso_atm_d14c_data(ciso_atm_d14c_data_nbval,3))

!-----------------------------------------------------------------------
!     READ in C14 data from files - three files, for SH, EQ, NH
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      do il=1,3
         open (nml_in,file=ciso_atm_d14c_filename(il),status='old',iostat=stat)
         if (stat /= 0) then
            write(stdout,*) 'open failed for ', trim(ciso_atm_d14c_filename(il))
            go to 99
         endif
         read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d14c_data_nbval_tmp
         if (stat /= 0) then
            write(stdout,*) '1st line read failed for ', trim(ciso_atm_d14c_filename(il))
            go to 99
         endif
         do irec=1,skiplines
            read(nml_in,FMT=*,iostat=stat) sglchr
            if (stat /= 0) then
               write(stdout,fmt=*) 'skipline read failed for ', trim(ciso_atm_d14c_filename(il))
               go to 99
            endif
         enddo
         do irec=1,ciso_atm_d14c_data_nbval
            read(nml_in,FMT=*,iostat=stat) ciso_atm_d14c_data_yr(irec,il), ciso_atm_d14c_data(irec,il)
            if (stat /= 0) then
               write(stdout,fmt=*) 'data read failed for ', trim(ciso_atm_d14c_filename(il))
               go to 99
            endif
         enddo
         close(nml_in)
      enddo
   endif

99 call broadcast_scalar(stat, master_task)
   if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
                                                          &/ sub_name)

!---------------------------------------------------------------------
! Broadcast the variables to other tasks beside master_task
!---------------------------------------------------------------------

   call broadcast_array(ciso_atm_d14c_data, master_task)
   call broadcast_array(ciso_atm_d14c_data_yr, master_task)

!-----------------------------------------------------------------------
!EOC

 end subroutine ciso_read_atm_D14C_data



!***********************************************************************
!BOP
! !IROUTINE: ciso_comp_varying_D13C
! !INTERFACE:

 subroutine ciso_comp_varying_D13C(iblock, ciso_data_ind_d13c, D13C)

! !DESCRIPTION:
!  Compute atmospheric mole fractions of d13c when temporarily
!  varying data is read from files
!  1. Linearly interpolate data values to current model time step
!  2. Spatial patern of D13Cis the same everywhere (90 S - 90 N)
!
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   integer (int_kind) :: &
      iblock           ! block index

   integer (int_kind) :: &
      ciso_data_ind_d13c  ! ciso_data_ind_d13cis the index for the data for current timestep,
                          ! note that ciso_data_ind_d13c is always strictly less than the length
                          ! of the data and is initialized to -1 before the first call

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      D13C             ! atmospheric D13C (permil)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j              ! loop indices

   real (r8) :: &
      model_date,     & ! date of current model timestep
      mapped_date,    & ! model_date mapped to data timeline
      weight            ! weighting for temporal interpolation


!-----------------------------------------------------------------------
!  Generate mapped_date and check to see if it is too large.
!-----------------------------------------------------------------------

   model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
   mapped_date = model_date - ciso_atm_model_year + ciso_atm_data_year

   if (mapped_date >= ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval)) then
      call exit_POP(sigAbort, 'ciso: Model date maps to date after end of D13C data in file.')
   endif

!--------------------------------------------------------------------------------------------------------------
!  Set atmospheric D13C to first value in record for years before record begins
!--------------------------------------------------------------------------------------------------------------

   if (mapped_date < ciso_atm_d13c_data_yr(1)) then
      D13C = ciso_atm_d13c_data(1)
      ciso_data_ind_d13c = 1
         if(my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,ndelim_fmt)
         write(stdout,blank_fmt)
         write(stdout,*)'ciso: Mapped date less than start of D13C data --> using first value in D13C data file'
         write(stdout,blank_fmt)
         write(stdout,ndelim_fmt)
         write(stdout,blank_fmt)
      endif
      return
   endif

!-----------------------------------------------------------------------
!  On first time step, perform linear search to find data_ind_d13c
!-----------------------------------------------------------------------

   if (ciso_data_ind_d13c == -1) then
      do ciso_data_ind_d13c = ciso_atm_d13c_data_nbval-1,1,-1
         if (mapped_date >= ciso_atm_d13c_data_yr(ciso_data_ind_d13c)) exit
      end do
   endif

!-----------------------------------------------------------------------
!  See if ciso_data_ind_d13c needs to be updated,
!  but do not set it to atm_d13c_data_nbval.
!-----------------------------------------------------------------------

   if (ciso_data_ind_d13c < ciso_atm_d13c_data_nbval-1) then
      if (mapped_date >= ciso_atm_d13c_data_yr(ciso_data_ind_d13c+1)) ciso_data_ind_d13c = ciso_data_ind_d13c + 1
   endif


!-----------------------------------------------------------------------
!  Generate hemisphere values for current time step.
!-----------------------------------------------------------------------

   weight = (mapped_date - ciso_atm_d13c_data_yr(ciso_data_ind_d13c)) &
            / (ciso_atm_d13c_data_yr(ciso_data_ind_d13c+1) - ciso_atm_d13c_data_yr(ciso_data_ind_d13c))

   D13C = weight * ciso_atm_d13c_data(ciso_data_ind_d13c+1) + (c1-weight) * ciso_atm_d13c_data(ciso_data_ind_d13c)

!-----------------------------------------------------------------------
!EOC

 end subroutine ciso_comp_varying_D13C

!***********************************************************************
!BOP
! !IROUTINE: ciso_comp_varying_D14C
! !INTERFACE:

 subroutine ciso_comp_varying_D14C(iblock, ciso_data_ind_d14c, D14C)

! !DESCRIPTION:
!  Compute atmospheric mole fractions of CO2 when temporarily
!  varying data is read from files
!  1. Linearly interpolate hemispheric values to current time step
!  2. Make global field of D14C, determined by:
!   -Northern Hemisphere value is used for 20N - 90 N
!   -Southern Hemisphere value is used for 20 S - 90 S
!   -Equator value is used for 20 S- 20 N
!
!
! !REVISION HISTORY:
!  same as module

! !USES:

!   use constants, only : radian
!   use grid, only : TLAT

! !INPUT PARAMETERS:

   integer (int_kind) :: &
      iblock          ! block index

   integer (int_kind) :: &
      ciso_data_ind_d14c   ! data_ind_d14c is the index into data for current timestep,
                      !  note that data_ind is always strictly less than the length of D14C data
                      !  and is initialized to -1 before the first call


! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      D14C            ! atmospheric delta C14 in permil on global grid

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, il        ! loop indices

   real (r8) :: &
      model_date,      & ! date of current model timestep
      mapped_date,     & ! model_date mapped to data timeline
      weight,          & ! weighting for temporal interpolation
      d14c_curr_sh,    & ! current atmospheric D14C value for SH (interpolated from data to model date)
      d14c_curr_nh,    & ! current atmospheric D14C value for NH (interpolated from data to model date)
      d14c_curr_eq       ! current atmospheric D14C value for EQ (interpolated from data to model date)

!-----------------------------------------------------------------------
!  Generate mapped_date and check to see if it is too large.
!-----------------------------------------------------------------------

   model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
   mapped_date = model_date - ciso_atm_model_year + ciso_atm_data_year
   do il=1,3
   if (mapped_date >= ciso_atm_d14c_data_yr(ciso_atm_d14c_data_nbval,il)) then
      call exit_POP(sigAbort, 'ciso: model date maps to date after end of D14C data in files.')
   endif
   enddo

!--------------------------------------------------------------------------------------------------------------
!  Set atmospheric D14C concentrations to zero before D14C record begins
!--------------------------------------------------------------------------------------------------------------

   if (mapped_date < ciso_atm_d14c_data_yr(1,1)) then
      D14C = c0
      ciso_data_ind_d14c = 1
      if(my_task == master_task) then
         write(stdout,*)'ciso: Model date less than start of D14C data --> D14C=0'
      endif
      return
   endif

!-----------------------------------------------------------------------
!  On first time step, perform linear search to find data_ind_d14c.
!-----------------------------------------------------------------------

   if (ciso_data_ind_d14c == -1) then
      do ciso_data_ind_d14c = ciso_atm_d14c_data_nbval-1,1,-1
         if (mapped_date >= ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1)) exit
      end do
   endif

!-----------------------------------------------------------------------
!  See if data_ind_d14c need to be updated,
!  but do not set it to atm_co2_data_nbval.
!-----------------------------------------------------------------------

   if (ciso_data_ind_d14c < ciso_atm_d14c_data_nbval-1) then
      if (mapped_date >= ciso_atm_d14c_data_yr(ciso_data_ind_d14c+1,1))  then
         ciso_data_ind_d14c = ciso_data_ind_d14c + 1
      endif
   endif
!
!-----------------------------------------------------------------------
!  Generate hemisphere values for current time step.
!-----------------------------------------------------------------------

   weight = (mapped_date - ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1)) &
            / (ciso_atm_d14c_data_yr(ciso_data_ind_d14c+1,1) - ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1))

   d14c_curr_sh = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,1) + &
                  (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,1)
   d14c_curr_eq = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,2) + &
                 (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,2)
   d14c_curr_nh = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,3) + &
                  (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,3)


!-----------------------------------------------------------------------
!  Merge hemisphere values for D14C
!      -Northern Hemisphere value is used for >20N - 90 N
!      -Southern Hemisphere value is used for >20 S - 90 S
!      -Equatorial value is used for 20 S to 20 N
!-----------------------------------------------------------------------

   do j = 1, ny_block
      do i = 1, nx_block
         if (TLATD(i,j,iblock) < -20.0_r8) then
            D14C(i,j) = d14c_curr_sh
         else if (TLATD(i,j,iblock) > 20.0_r8) then
            D14C(i,j) = d14c_curr_nh
         else
            D14C(i,j) = d14c_curr_eq
         endif
      end do
   end do
!-----------------------------------------------------------------------
!EOC

 end subroutine ciso_comp_varying_D14C


!*****************************************************************************
!BOP
! !IROUTINE: ecosys_ciso_tavg_forcing
! !INTERFACE:

 subroutine ecosys_ciso_tavg_forcing(STF_MODULE)

! !DESCRIPTION:
!  Accumulate non-standard forcing related tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real (r8), dimension(:,:,:,:), &
     intent(in) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

!-----------------------------------------------------------------------
!  Do components of flux calculations saved in hack array
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1,nblocks_clinic
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_13CO2,iblock), tavg_CISO_DI13C_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_as_13CO2,iblock), tavg_CISO_DI13C_AS_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_sa_13CO2,iblock), tavg_CISO_DI13C_SA_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d13C,iblock), tavg_CISO_d13C_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R13C_DIC_surf,iblock), tavg_CISO_R13C_DIC_surf,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R13C_atm,iblock), tavg_CISO_R13C_atm,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_D13C_atm,iblock), tavg_CISO_D13C_atm,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_eps_aq_g_surf,iblock), tavg_CISO_eps_aq_g_surf,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_eps_dic_g_surf,iblock), tavg_CISO_eps_dic_g_surf,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_14CO2,iblock), tavg_CISO_DI14C_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_as_14CO2,iblock), tavg_CISO_DI14C_AS_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_sa_14CO2,iblock), tavg_CISO_DI14C_SA_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R14C_DIC_surf,iblock), tavg_CISO_R14C_DIC_surf,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R14C_atm,iblock), tavg_CISO_R14C_atm,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_D14C_atm,iblock), tavg_CISO_D14C_atm,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d14C,iblock), tavg_CISO_d14C_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DI13C_RIV_FLUX,iblock),tavg_CISO_DI13C_RIV_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DO13C_RIV_FLUX,iblock),tavg_CISO_DO13C_RIV_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DI14C_RIV_FLUX,iblock),tavg_CISO_DI14C_RIV_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DO14C_RIV_FLUX,iblock),tavg_CISO_DO14C_RIV_FLUX,iblock,1)
 !debugging
      call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_GLOBAL_D14C,iblock), tavg_CISO_GLOBAL_D14C,iblock,1)
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_ciso_tavg_forcing

 !*****************************************************************************
!BOP
! !IROUTINE: ecosys_ciso_tracer_ref_val
! !INTERFACE:

 function ecosys_ciso_tracer_ref_val(ind)

! !DESCRIPTION:
!  return reference value for tracers using virtual fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: ind

! !OUTPUT PARAMETERS:

   real (r8) :: ecosys_ciso_tracer_ref_val

!EOP
!BOC
!-----------------------------------------------------------------------

   if (ciso_vflux_flag(ind)) then
      ecosys_ciso_tracer_ref_val = ciso_surf_avg(ind)
   else
      ecosys_ciso_tracer_ref_val = c0
   endif

!-----------------------------------------------------------------------
!EOC

 end function ecosys_ciso_tracer_ref_val


!*****************************************************************************
!BOP
! !IROUTINE: ecosys_ciso_write_restart
! !INTERFACE:

 subroutine ecosys_ciso_write_restart(restart_file, action)

! !DESCRIPTION:
!  write auxiliary fields & scalars to restart files
!
! !REVISION HISTORY:
!  same as module
!  use constants, only: char_blank, field_loc_center, field_type_scalar

! !INPUT PARAMETERS:

   character(*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character (char_len) :: &
      short_name   ! tracer name temporaries

   type (io_dim) :: &
      i_dim, j_dim ! dimension descriptors

   integer (int_kind) :: n

!-----------------------------------------------------------------------

   if (trim(action) == 'add_attrib_file') then
      short_name = char_blank
      do n=1,ecosys_ciso_tracer_cnt
         if (ciso_vflux_flag(n)) then
            short_name = 'surf_avg_' /&
                      &/ ciso_ind_name_table(n)%name
            call add_attrib_file(restart_file,trim(short_name),ciso_surf_avg(n))
         endif
      end do
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_ciso_write_restart
!***********************************************************************

 end module ecosys_ciso_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


