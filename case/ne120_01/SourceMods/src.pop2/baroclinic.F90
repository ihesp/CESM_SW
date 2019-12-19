!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!!! added_start by wangli
#ifndef NOATH
#define WLATH
#define ZYUATH
#endif
!!! added_end by wangli

 module baroclinic

!BOP
! !MODULE: baroclinic
!
! !DESCRIPTION:
!  Contains main driver routines and variables for computing the 
!  baroclinic velocities and tracer fields.
!
! !REVISION HISTORY:
!  SVN:$Id: baroclinic.F90 87926 2017-12-12 00:31:20Z nanr $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod 
   use POP_GridHorzMod
   use POP_HaloMod 

!   use perf_mod

   use kinds_mod, only: int_kind, r8, log_kind, r4, rtavg
   use blocks, only: nx_block, ny_block, block, get_block
!   use distribution, only: 
   use domain_size
   use domain, only: nblocks_clinic, blocks_clinic, POP_haloClinic
   use constants, only: delim_fmt, blank_fmt, p5, field_loc_center,p25,       &
       field_type_scalar, c0, c1, c2, grav, ndelim_fmt,                      &
       hflux_factor, salinity_factor, salt_to_ppt
!! modified by wangli to expand gradp and grad

   use prognostic, only: TRACER, UVEL, VVEL, max_blocks_clinic, km, mixtime, &
       RHO, newtime, oldtime, curtime, PSURF, nt
   use broadcast, only: broadcast_scalar
   use communicate, only: my_task, master_task
   use grid, only: FCOR, DZU, HUR, KMU, KMT, sfc_layer_type, l1Ddyn,         &
       sfc_layer_varthick, partial_bottom_cells, dz, DZT, CALCT, dzw, &
       dzr,dxur,dyur
!! modified by wangli /above to expand gradp and its inline sub grad

   use advection, only: advu, advt, advthdifft,comp_flux_vel_ghost
   use pressure_grad, only: lpressure_avg, gradp,bouss
   use horizontal_mix, only: hdiffu, advuhdiffu,hdifft, iso_impvmixt_tavg
   use vertical_mix, only: vmix_coeffs, implicit_vertical_mix, vdiffu,       &
       vdifft, impvmixt, impvmixu, impvmixt_correct, convad, impvmixt_tavg
   use vmix_kpp, only: add_kpp_sources
   use diagnostics, only: ldiag_cfl, cfl_check, ldiag_global,                &
       DIAG_KE_ADV_2D, DIAG_KE_PRESS_2D, DIAG_KE_HMIX_2D, DIAG_KE_VMIX_2D,   &
       DIAG_TRACER_HDIFF_2D, DIAG_PE_2D, DIAG_TRACER_ADV_2D,                 &
       DIAG_TRACER_SFC_FLX, DIAG_TRACER_VDIFF_2D, DIAG_TRACER_SOURCE_2D
   use movie, only: define_movie_field, movie_requested, update_movie_field
   use state_mod, only: state
   use state_mod 
   use ice, only: liceform, ice_formation, increment_tlast_ice
   use time_management, only: mix_pass, leapfrogts, impcor, c2dtu, beta,     &
       gamma, c2dtt
   use io_types, only: nml_in, nml_filename, stdout
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now, &
       tavg_method_max, tavg_method_min  &
       , TAVG_BUF_2D, TAVG_BUF_3D       & ! added by zengyh 2019-02-14.
       , dtavg, zt_150m_levs            &
       , avail_tavg_fields        ! added by zengyh 2019-02-14.	   
   use forcing_fields, only: STF, SMF, lsmft_avail, SMFT, TFW
   use forcing_shf, only: SHF_QSW
   use forcing_sfwf, only: lfw_as_salt_flx
   use sw_absorption, only:  add_sw_absorb
   use forcing_pt_interior, only: set_pt_interior
   use forcing_s_interior, only: set_s_interior
   use passive_tracers, only: set_interior_passive_tracers,  &
       reset_passive_tracers, tavg_passive_tracers, &
       tavg_passive_tracers_baroclinic_correct, &
       set_interior_passive_tracers_3D
   use exit_mod, only: sigAbort, exit_pop, flushm
   use overflows
   use overflow_type
!!   use operators, only: grad !!added by wangli


! added by zengyh 2018-11-28:
   use timers
   ! use step_mod
! added by zengyh 2019-2-15:   
   use io_tools, only: document
   use communicate, only: my_task, master_task   

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_baroclinic,          &
             baroclinic_driver,        &
             baroclinic_correct_adjust

! !PRIVATE DATA MEMBERS:

   logical (log_kind) :: &
      reset_to_freezing   ! flag to prevent very cold water

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  ids for tavg diagnostics computed from baroclinic
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_UDP,          &! tavg id for pressure grad work
      tavg_TEMP,         &! tavg id for temperature
      tavg_TEMP_2,       &! tavg id for temperature
      tavg_TEMP_MAX,     &! tavg id for maximum temperature
      tavg_TEMP_MIN,     &! tavg id for maximum temperature
      tavg_dTEMP_POS_3D, &! tavg id for positive temperature timestep difference
      tavg_dTEMP_POS_2D, &! tavg id for positive temperature timestep difference
      tavg_dTEMP_NEG_3D, &! tavg id for negative temperature timestep difference
      tavg_dTEMP_NEG_2D, &! tavg id for negative temperature timestep difference
      tavg_SST,          &! tavg id for surface temperature
      tavg_SST2,         &! tavg id for surface temperature squared
      tavg_SALT,         &! tavg id for salinity
      tavg_SALT_2,       &! tavg id for salinity	
      tavg_SALT_MAX,     &! tavg id for maximum salinity
      tavg_SALT_MIN,     &! tavg id for minimum salinity
      tavg_TEMP2,        &! tavg id for temperature squared
      tavg_SALT2,        &! tavg id for salinity    squared
      tavg_UVEL,         &! tavg id for U velocity
      tavg_UVEL_2,       &! tavg id for U velocity	
      tavg_UVEL2,        &! tavg id for U velocity squared
      tavg_VVEL,         &! tavg id for V velocity
      tavg_VVEL_2,       &! tavg id for V velocity
      tavg_VVEL2,        &! tavg id for V velocity squared
      tavg_KE,           &! tavg id for kinetic energy
      tavg_KE_2,           &! tavg id for kinetic energy
      tavg_ST,           &! tavg id for salt*temperature
      tavg_RHO,          &! tavg id for in-situ density
      tavg_RHO_VINT,     &! tavg id for vertical integral of in-situ density
      tavg_UV,           &! tavg id for u times v
      tavg_T1_8,         &! tavg id for temperature in top 8 lvls
      tavg_S1_8,         &! tavg id for salinity    in top 8 lvls
      tavg_U1_8,         &! tavg id for U           in top 8 lvls
      tavg_V1_8,         &! tavg id for V           in top 8 lvls
      tavg_U1_1,         &! tavg id for U           in top 1 lvl
      tavg_V1_1,         &! tavg id for V           in top 1 lvl
      tavg_RESID_T,      &! free-surface residual flux (T)
      tavg_RESID_S        ! free-surface residual flux (S)

!-----------------------------------------------------------------------
!
!  ids for movie diagnostics computed from baroclinic
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(km) :: &
      movie_TEMP,         &! movie id for temperature
      movie_SALT,         &! movie id for salinity
      movie_UVEL,         &! movie id for U velocity
      movie_VVEL,         &! movie id for V velocity
      movie_RHO            ! movie id for in-situ density


!!! added-start by wangli to expand gradp and its inline subroutine grad
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      SUMX, SUMY,        &! incremental k sum of Grad{x,y}(P(k))
      RHOKMX, RHOKMY      ! x,y gradient of density at k-1 level
!!! added-end by wangli


!!! defined by zengyh @ 2019-01-31: 
   character (128) :: exit_string
#if defined(POPCHKJN) || defined(POPTIMEJN)
   integer,save :: numsteps=0
#if defined(TAVG_R8)
   real(r8), parameter :: &
#else
   real(r4), parameter :: &
#endif
        f0=0.0, f1=1.0
!-----------timing_var_def_zengyh----------
   real(8)::start_accum,stop_accum
#endif
!!! end_defined by zengyh @ 2019-01-31.


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_baroclinic
! !INTERFACE:

 subroutine init_baroclinic

! !DESCRIPTION:
!  Initializes some baroclinic options.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error          ! namelist i/o error flag
   integer (int_kind) :: k

   namelist /baroclinic_nml/ reset_to_freezing

!-----------------------------------------------------------------------
!
!  read options from namelist and broadcast
!
!-----------------------------------------------------------------------
!   !call timer_start(clinic_init_sec1)
   reset_to_freezing = .true.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=baroclinic_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading baroclinic_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Baroclinic:'
      write(stdout,blank_fmt)
      write(stdout,*) ' baroclinic_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,baroclinic_nml)
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)

      write(stdout,'(a18)') 'Baroclinic options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      if (reset_to_freezing) then
         write(stdout,'(a40)') &
                          'Surface temperature reset to freezing on'
      else
         write(stdout,'(a41)') &
                          'Surface temperature reset to freezing off'
      endif
   endif

   call broadcast_scalar(reset_to_freezing, master_task)

!   !call timer_stop(clinic_init_sec1)
!-----------------------------------------------------------------------
!
!  define tavg fields computed from baroclinic driver routines
!
!-----------------------------------------------------------------------

!   !call timer_start(clinic_init_sec2)
   call define_tavg_field(tavg_UDP,'UDP',3,                            &
                          long_name='Pressure work',                   &
                          units='erg', grid_loc='3221')

   call define_tavg_field(tavg_U1_1,'U1_1',2,                          &
                          long_name='Zonal Velocity lvls 1-1',         &
                          units='centimeter/s', grid_loc='2221')

   call define_tavg_field(tavg_V1_1,'V1_1',2,                          &
                          long_name='Meridional Velocity lvls 1-1',    &
                          units='centimeter/s', grid_loc='2221')

   call define_tavg_field(tavg_U1_8,'U1_8',2,                          &
                          long_name='Zonal Velocity lvls 1-8',         &
                          units='centimeter/s', grid_loc='2221')

   call define_tavg_field(tavg_V1_8,'V1_8',2,                          &
                          long_name='Meridional Velocity lvls 1-8',    &
                          units='centimeter/s', grid_loc='2221')

   call define_tavg_field(tavg_T1_8,'T1_8',2,                          &
                          long_name='Potential Temperature lvls 1-8',  &
                          units='degC', grid_loc='2111')

   call define_tavg_field(tavg_S1_8,'S1_8',2,                          &
                          long_name='Salinity lvls 1-8',               &
                          scale_factor=1000.0_rtavg,                      &
                          units='gram/kilogram', grid_loc='2111')

   call define_tavg_field(tavg_UVEL,'UVEL',3,                          &
                          long_name='Velocity in grid-x direction',    &
                          units='centimeter/s', grid_loc='3221',       &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_UVEL_2,'UVEL_2',3,                          &
                          long_name='Velocity in grid-x direction',    &
                          units='centimeter/s', grid_loc='3221',       &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_UVEL2,'UVEL2',3,                        &
                          long_name='Velocity**2 in grid-x direction', &
                          units='centimeter^2/s^2', grid_loc='3221',   &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_VVEL,'VVEL',3,                          &
                          long_name='Velocity in grid-y direction',    &
                          units='centimeter/s', grid_loc='3221',       &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_VVEL_2,'VVEL_2',3,                          &
                          long_name='Velocity in grid-y direction',    &
                          units='centimeter/s', grid_loc='3221',       &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_VVEL2,'VVEL2',3,                        &
                          long_name='Velocity**2 in grid-y direction', &
                          units='centimeter^2/s^2', grid_loc='3221',   &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_KE,'KE',3,                              &
                          long_name='Horizontal Kinetic Energy',       &
                          units='centimeter^2/s^2', grid_loc='3221',   &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_KE_2,'KE_2',3,                              &
                          long_name='Horizontal Kinetic Energy',       &
                          units='centimeter^2/s^2', grid_loc='3221',   &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_TEMP,'TEMP',3,                          &
                          long_name='Potential Temperature',           &
                          units='degC', grid_loc='3111',               &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_TEMP_2,'TEMP_2',3,                          &
                          long_name='Potential Temperature',           &
                          units='degC', grid_loc='3111',               &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_TEMP_MAX,'TEMP_MAX',3,                  &
                          tavg_method=tavg_method_max,                 &
                          long_name='Maximum Potential Temperature',   &
                          units='degC', grid_loc='3111',               &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_TEMP_MIN,'TEMP_MIN',3,                  &
                          tavg_method=tavg_method_min,                 &
                          long_name='Minimum Potential Temperature',   &
                          units='degC', grid_loc='3111',               &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_dTEMP_POS_3D,'dTEMP_POS_3D',3,          &
                          tavg_method=tavg_method_max,                 &
                          long_name='max pos temperature timestep diff', &
                          units='degC', grid_loc='3111',               &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_dTEMP_POS_2D,'dTEMP_POS_2D',2,          &
                          tavg_method=tavg_method_max,                 &
                          long_name='max pos column temperature timestep diff', &
                          units='degC', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_dTEMP_NEG_3D,'dTEMP_NEG_3D',3,          &
                          tavg_method=tavg_method_min,                 &
                          long_name='min neg temperature timestep diff', &
                          units='degC', grid_loc='3111',               &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_dTEMP_NEG_2D,'dTEMP_NEG_2D',2,          &
                          tavg_method=tavg_method_min,                 &
                          long_name='min neg column temperature timestep diff', &
                          units='degC', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_SST,'SST',2,                            &
                          long_name='Surface Potential Temperature',   &
                          units='degC', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_SST2,'SST2',2,                          &
                          long_name='Surface Potential Temperature**2',&
                          units='degC', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_SALT,'SALT',3,                          &
                          long_name='Salinity',                        &
                          units='gram/kilogram', grid_loc='3111',      &
                          scale_factor=1000.0_rtavg,                      &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_SALT_2,'SALT_2',3,                          &
                          long_name='Salinity',                        &
                          units='gram/kilogram', grid_loc='3111',      &
                          scale_factor=1000.0_rtavg,                      &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_SALT_MAX,'SALT_MAX',3,                  &
                          tavg_method=tavg_method_max,                 &
                          long_name='Maximum Salinity',                &
                          units='gram/kilogram', grid_loc='3111',      &
                          scale_factor=1000.0_rtavg,                   &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_SALT_MIN,'SALT_MIN',3,                  &
                          tavg_method=tavg_method_min,                 &
                          long_name='Minimum Salinity',                &
                          units='gram/kilogram', grid_loc='3111',      &
                          scale_factor=1000.0_rtavg,                   &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_TEMP2,'TEMP2',3,                        &
                          long_name='Temperature**2',                  &
                          units='degC^2', grid_loc='3111')

   call define_tavg_field(tavg_SALT2,'SALT2',3,                        &
                          long_name='Salinity**2 ',                    &
                          units='(gram/gram)^2', grid_loc='3111',      &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_ST,'ST',3,                              &
                          long_name='Temperature*Salinity',            &
                          units='degC*gram/gram', grid_loc='3111',     &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_RHO,'RHO',3,                            &
                          long_name='In-Situ Density',                 &
                          units='gram/centimeter^3', grid_loc='3111',  &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_RHO_VINT,'RHO_VINT',2,                  &
                          long_name='Vertical Integral of In-Situ Density', &
                          units='gram/centimeter^2', grid_loc='2110',  &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_UV,'UV',3,                              &
                          long_name='UV velocity product',             &
                          units='centimeter^2/s^2', grid_loc='3221')

   call define_tavg_field(tavg_RESID_T,'RESID_T',2,                           &
                    long_name='Free-Surface Residual Flux (T)',               &
                          units='watt/m^2', grid_loc='2110',&
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_RESID_S,'RESID_S',2,                           &
                    long_name='Free-Surface Residual Flux (S)',               &
                          units='kg/m^2/s', grid_loc='2110',&
                          coordinates='TLONG TLAT time')

!   !call timer_stop(clinic_init_sec2)
!-----------------------------------------------------------------------
!
!  define movie fields computed from baroclinic driver routines
!
!-----------------------------------------------------------------------

!   !call timer_start(clinic_init_sec3)
   do k = 1, km

      call define_movie_field(movie_UVEL(k),'UVEL',k,                  &
                          long_name='Zonal Velocity',                  &
                          units='cm/s', grid_loc='3221')

      call define_movie_field(movie_VVEL(k),'VVEL',k,                  &
                          long_name='Meridional Velocity',             &
                          units='cm/s', grid_loc='3221')

      call define_movie_field(movie_TEMP(k),'TEMP',k,                  &
                          long_name='Potential Temperature',           &
                          units='degC', grid_loc='3111')

      call define_movie_field(movie_SALT(k),'SALT',k,                  &
                          long_name='Salinity',                        &
                          units='psu', grid_loc='3111')

      call define_movie_field(movie_RHO(k),'RHO',k,                    &
                          long_name='In-situ density',                 &
                          units='sigma units', grid_loc='3111')

   enddo

!   !call timer_stop(clinic_init_sec3)
!-----------------------------------------------------------------------
!EOC

 call flushm (stdout)

 end subroutine init_baroclinic

!***********************************************************************
!BOP
! !IROUTINE: baroclinic_driver
! !INTERFACE:

 subroutine baroclinic_driver(ZX,ZY,DH,DHU, errorCode)

! !DESCRIPTION:
!  This routine is the main driver for the explicit time integration of
!  baroclinic velocities $(u',v')$ and tracer fields $T$. 
!
!  Tracer equations:
!  \begin{equation}
!     (T^{n+1}-T^{n-1})/(2 \Delta t) = -L(T^n) + D_H(T^{n-1}) + 
!                                      D_V(T^{n-1}) + S
!  \end{equation}
!  where $S$ are source terms, $L$ is the advection operator and
!  $D$ is the diffusion operator in the horizontal and vertical
!  directions.
!
!  Momentum equations:
!  \begin{eqnarray}
!     (u'-u^{n-1})/(2 \Delta t) - f*\alpha*(v'-v^{n-1}) &=& F_x \\
!     (v'-v^{n-1})/(2 \Delta t) + f*\alpha*(u'-u^{n-1}) &=& F_y
!  \end{eqnarray}
!  \begin{eqnarray}
!     \tilde{u}' &=& u' - {1\over{H_U}}\sum_{k=1}^{km}dz_k u_k' \\
!     \tilde{v}' &=& v' - {1\over{H_U}}\sum_{k=1}^{km}dz_k v_k'
!  \end{eqnarray}
!
!  This routine calculates baroclinic velocities and tracers at
!  the new time level and stores them in arrays UVEL,VVEL,and TRACER
!  with time index newtime.  The arrays UVEL,VVEL, and TRACER
!  with time indices curtime and oldtime are not updated for the next
!  timestep until near the end of routine step.
!
!  The above equations are written for the case of (leapfrog)
!  implicit treatment of the coriolis terms.  if these terms are
!  treated explicitly then the coriolis terms appear only in the 
!  forcing terms $(F_x,F_y)$, which are calculated in clinic. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(in) :: &
      DH, DHU              ! change in surface height at T,U points

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(out) :: &
      ZX, ZY               ! vertical integrals of forcing

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      i,j,                &! dummy indices for horizontal directions
      n,k,                &! dummy indices for vertical level, tracer
      iblock,             &! counter for block loops
      kp1,km1,partial             ! level index for k+1, k-1 levels

   real (r8), dimension(nx_block,ny_block) :: & 
      FX,FY,              &! sum of r.h.s. forcing terms
      WORK1,WORK2,        &! local work space
      WUK,                &! vertical velocity at top of U box
      WTK                  ! vertical velocity at top of T box

   real (r8)    ::        &
      factor
   type (block) ::        &
      this_block           ! block information for current block
#ifndef NOWL_ATH
   integer(int_kind):: iparam(8)
   real(r8)::          rparam(8)
#endif


! added by zengyh 2019-02-15:
#if !defined(NOATH_TAVG_CLINIC_ZYH)

   integer (int_kind) ::  &
      itmp,               &
      bufloc,             &! location of field in tavg buffer
      ndims              ! rank of field (2=2d,3=3d)
   integer (int_kind), save ::  &
      part_bott_cells,    &! transform logical into int 
      numtavg            ! number of all tavg_ids
   integer (int_kind), dimension(48), save ::  &
      tavgids=0,              & ! tavg_id
      tavglocs=-1,            & ! bufloc of tavg_id
      tavgdims=0,             & ! ndims of tavg_id
      tavgmeths=-1               ! tavg method of tavg_id
   integer (int_kind), dimension(48*km) ::  &
      tavgflags               ! accum flag of tavg_id
   integer (int_kind), dimension(nx_block,ny_block) ::  &
      intCALCT                  ! transform logical CALCT(:,:) into int 
  
#if defined(POPCHKJN) 
#if defined(TAVG_R8)
   real (r8), dimension(nx_block,ny_block,km,48) :: &
      initavg3d=c0, tmptavg3d=c0, tmprtavg3d=c0
   real (r8), dimension(48) :: &
      errtavg3d, rerrtavg3d,             &
      maxerrtavg3d, maxrerrtavg3d,       &
      minerrtavg3d
#else
   real (r4), dimension(nx_block,ny_block,km,48) :: &
      initavg3d=f0, tmptavg3d=f0, tmprtavg3d=f0
   real (r4), dimension(48) :: &
      errtavg3d, rerrtavg3d,             &
      maxerrtavg3d, maxrerrtavg3d,       &
      minerrtavg3d
#endif  
   integer (int_kind) ::     &
      mykm = -1
   integer (int_kind), dimension(48,5) ::&
      cstat
#endif  


   interface
     subroutine m_clinic_tavg_sw(TAVG_BUF_2D,TAVG_BUF_3D, &
            UVEL,VVEL, &
            TRACER1cur,TRACER1old,TRACER2cur, &
            RHO,PSURF,DZT,DH, &
            dtavg,dzkm,grav,hflux_factor,salinity_factor, &
            KMT,intCALCT, &
            tavgids,tavgflags,tavglocs,tavgdims,tavgmeths, &
            myproc,km,nx_block,ny_block,iblock,niblock,&
            sfc_layer_type,sfc_layer_varthick,numtavg,part_bott_cells &
          ) bind(C, name='m_clinic_tavg_sw')
         integer,value:: myproc,km,nx_block,ny_block,iblock,niblock,sfc_layer_type,sfc_layer_varthick,numtavg,part_bott_cells
         integer:: KMT(*),intCALCT(*),tavgids(*),tavgflags(*),tavglocs(*),tavgdims(*),tavgmeths(*)
         real(8),value:: dtavg,grav,hflux_factor,salinity_factor
#ifdef TAVG_R8
         real(8):: TAVG_BUF_2D(*),TAVG_BUF_3D(*)
#else
         real(4):: TAVG_BUF_2D(*),TAVG_BUF_3D(*)
#endif
         real(8):: UVEL(*),VVEL(*),TRACER1cur(*),TRACER1old(*),TRACER2cur(*),RHO(*),PSURF(*),dzkm(*),DZT(*),DH(*)
     end subroutine m_clinic_tavg_sw
   end interface

#endif
! end_added by zengyh 2019-02-15.



!-----------------------------------------------------------------------
!
!  compute flux velocities in ghost cells
!
!-----------------------------------------------------------------------

!   !call timer_start(driver_velghost)
   errorCode = POP_Success

   if (.not.l1Ddyn) then

     call comp_flux_vel_ghost(DH, errorCode)

     if (errorCode /= POP_Success) then
        call POP_ErrorSet(errorCode, &
           'baroclinic_driver: error in comp_flux_vel_ghost')
        return
     endif

   endif

!   !call timer_stop(driver_velghost)
!-----------------------------------------------------------------------
!
!  prior to tracer update, compute 3D passive_tracer terms
!
!-----------------------------------------------------------------------

   !call timer_start(driver_tracer3D)
!   call set_interior_passive_tracers_3D(                            &
!           TRACER (:,:,:,:,oldtime,:), TRACER (:,:,:,:,curtime,:)  )
!!!! commented by wangli, null calls /cost over 110s
   !call timer_stop(driver_tracer3D)
!-----------------------------------------------------------------------
!
!  first block loop to update tracers
!
!-----------------------------------------------------------------------


   !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,kp1,km1,WTK,WORK1,factor)

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)  

! added by zengyh 2019-02-15:
!----------------------------
#if !defined(NOATH_TAVG_CLINIC_ZYH)
      do j=1,ny_block
         do i=1,nx_block
           if (CALCT(i,j,iblock)) then
              intCALCT(i,j) = 1
           else 
              intCALCT(i,j) = 0 
           endif
        enddo
      enddo
#endif
!----------------------------
! end_added by zengyh 2019-02-15.


      do k = 1,km 

   !call timer_start(driver_kpp)
         kp1 = k+1
         km1 = k-1
         if (k == 1) km1 = 1
         if (k == km) kp1 = km

!-----------------------------------------------------------------------
!
!        compute vertical viscosity and diffusion coeffs
!
!-----------------------------------------------------------------------


         if (lsmft_avail) then
            call vmix_coeffs(k,TRACER (:,:,:,:,mixtime,iblock), &
                               UVEL   (:,:,:  ,mixtime,iblock), &
                               VVEL   (:,:,:  ,mixtime,iblock), &
                               UVEL   (:,:,:  ,curtime,iblock), &
                               VVEL   (:,:,:  ,curtime,iblock), &
                               RHO    (:,:,:  ,mixtime,iblock), &
                               STF    (:,:,:          ,iblock), &
                               SHF_QSW(:,:            ,iblock), &
                               this_block, SMFT=SMFT(:,:,:,iblock))
         else
            call vmix_coeffs(k,TRACER (:,:,:,:,mixtime,iblock), &
                               UVEL   (:,:,:  ,mixtime,iblock), &
                               VVEL   (:,:,:  ,mixtime,iblock), &
                               UVEL   (:,:,:  ,curtime,iblock), &
                               VVEL   (:,:,:  ,curtime,iblock), &
                               RHO    (:,:,:  ,mixtime,iblock), &
                               STF    (:,:,:          ,iblock), &
                               SHF_QSW(:,:            ,iblock), &
                               this_block, SMF=SMF(:,:,:,iblock))
         endif

   !call timer_stop(driver_kpp)
!-----------------------------------------------------------------------
!
!        calculate level k tracers at new time
!
!-----------------------------------------------------------------------


   !call timer_start(driver_tracer_up)
         call tracer_update(k, WTK,                             &
                               TRACER (:,:,:,:,newtime,iblock), &
                               TRACER (:,:,:,:,oldtime,iblock), &
                               TRACER (:,:,:,:,mixtime,iblock), &
                               TRACER (:,:,:,:,curtime,iblock), &
                               UVEL   (:,:,:  ,curtime,iblock), &
                               VVEL   (:,:,:  ,curtime,iblock), &
                               UVEL   (:,:,:  ,mixtime,iblock), &
                               VVEL   (:,:,:  ,mixtime,iblock), &
                               RHO    (:,:,:  ,curtime,iblock), &
                               STF    (:,:,:          ,iblock), &
                               TFW    (:,:,:          ,iblock), &
                               SHF_QSW(:,:            ,iblock), &
                               DH     (:,:            ,iblock), &
                               PSURF  (:,:    ,oldtime,iblock), &
                               PSURF  (:,:    ,curtime,iblock), &
                               this_block)
         

   !call timer_stop(driver_tracer_up)

!-----------------------------------------------------------------------
!
!        accumulate some tavg diagnostics if requested; testing internal to
!          accumulate_tavg_field
!
!-----------------------------------------------------------------------

   !call timer_start(driver_accum_tavg)
         if (mix_pass /= 1) then


! added by zengyh 2019-02-15:
!----------------------------

#if defined(POPTIMEJN)
         call getlooptime(start_accum)
#endif

#if !defined(NOATH_TAVG_CLINIC_ZYH)

         tavgids(1) = tavg_UVEL
         tavgids(2) = tavg_UVEL_2
         tavgids(3) = tavg_UVEL2
         tavgids(4) = tavg_U1_1
         tavgids(5) = tavg_U1_8

         tavgids(6) = tavg_VVEL
         tavgids(7) = tavg_VVEL_2
         tavgids(8) = tavg_VVEL2
         tavgids(9) = tavg_V1_1
         tavgids(10) = tavg_V1_8

         tavgids(11) = tavg_KE
         tavgids(12) = tavg_KE_2
         tavgids(13) = tavg_UV

         tavgids(14) = tavg_TEMP
         tavgids(15) = tavg_TEMP_2
         tavgids(16) = tavg_TEMP_MAX
         tavgids(17) = tavg_TEMP_MIN

         tavgids(18) = tavg_dTEMP_POS_3D
         tavgids(19) = tavg_dTEMP_POS_2D
         tavgids(20) = tavg_dTEMP_NEG_3D
         tavgids(21) = tavg_dTEMP_NEG_2D

         tavgids(22) = tavg_T1_8
         tavgids(23) = tavg_TEMP2
         tavgids(24) = tavg_SST
         tavgids(25) = tavg_SST2

         tavgids(26) = tavg_SALT
         tavgids(27) = tavg_SALT_2
         tavgids(28) = tavg_SALT_MAX
         tavgids(29) = tavg_SALT_MIN

         tavgids(30) = tavg_S1_8
         tavgids(31) = tavg_SALT2

         tavgids(32) = tavg_ST

         tavgids(33) = tavg_RHO
         tavgids(34) = tavg_RHO_VINT

         tavgids(35) = tavg_RESID_T
         tavgids(36) = tavg_RESID_S

         numtavg = 36
 
         if(partial_bottom_cells) part_bott_cells = 1


! step into manycore tavg codes:
         do i=1,numtavg
            tavgflags((i-1)*km+k) = 1
         enddo


         if(k>1) then 
            tavgflags((4 -1)*km+k) = 0
            tavgflags((9 -1)*km+k) = 0
            tavgflags((24-1)*km+k) = 0
            tavgflags((25-1)*km+k) = 0
         endif

         if(k>8) then
            tavgflags((5 -1)*km+k) = 0
            tavgflags((10-1)*km+k) = 0
            tavgflags((22-1)*km+k) = 0
            tavgflags((30-1)*km+k) = 0
         endif



         do i=1,numtavg 
            if (.not. accumulate_tavg_now(tavgids(i))) then
              tavgflags((i-1)*km+k)=0    
            !--------------------------------------------------------
            ! get buffer location and field info from avail_tavg_field array
            !---------------------------------------------------------
            else
              bufloc = avail_tavg_fields(tavgids(i))%buf_loc
              if (bufloc <= 0) then
                  exit_string = 'FATAL ERROR: attempt to accumulate bad tavg field'
                  call document ('accumulate_tavg_field', exit_string)
                  call exit_POP (sigAbort,exit_string,out_unit=stdout)
              else
                  tavglocs(i) = bufloc
                  ndims = avail_tavg_fields(tavgids(i))%ndims
                  tavgdims(i)=ndims
                  tavgmeths(i)=avail_tavg_fields(tavgids(i))%method
              endif
              if(i /= 19 .and. i /= 21) then
                if ((ndims == 3) .and. (k > zt_150m_levs) .and. &
                    (avail_tavg_fields(tavgids(i))%grid_loc(4:4) == '4')) then
                  tavgflags((i-1)*km+k)=0
                endif
              else
                if ((ndims == 3) .and. (1 > zt_150m_levs) .and. &
                    (avail_tavg_fields(tavgids(i))%grid_loc(4:4) == '4')) then
                  tavgflags((i-1)*km+k)=0
                endif
              endif
            endif
         enddo ! end of loop i 
#endif



#if defined(POPTIMEJN)
         call getlooptime(stop_accum)
         if(my_task==0 .and. mod(numsteps,64)==0 .and. k<=1) then  
            WRITE (*,*) 'pre_athreaded accum_tavg of clinic_driver costs zengyh:',stop_accum-start_accum,'(us);'
         endif
#endif


#ifdef POPTESTJN
         if(my_task==0 .and. numsteps<1 .and. k<=2) then  
            do i=1,36 
               write(*,*) 'set no. tavgid, bufloc, ndims, method, flags: ', i, tavgids(i),tavglocs(i),tavgdims(i),tavgmeths(i),tavgflags((i-1)*km+k)
            enddo
         endif
#endif


#if defined(POPCHKJN)
         do i=1,36
           if(tavgdims(i)==3) then
             initavg3d(:,:,k,i) = TAVG_BUF_3D(:,:,k,iblock,tavglocs(i))
           elseif(tavgdims(i)==2 .and. k==1) then
             initavg3d(:,:,k,i) = TAVG_BUF_2D(:,:,iblock,tavglocs(i))
           endif
         enddo

#if defined(POPTIMEJN)
         call getlooptime(start_accum)
#endif
#endif

#if defined(POPCHKJN) || defined(NOATH_TAVG_CLINIC_ZYH)
! original tavg codes:

! end_added by zengyh 2019-02-15.
!----------------------------


         call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock),tavg_UVEL,iblock,k)

         call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock),tavg_UVEL_2,iblock,k)

         call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock)**2,tavg_UVEL2,iblock,k)

         if (k <= 1)  &
            call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock), &
                                       tavg_U1_1,iblock,k)
         if (k <= 8)  &
         call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock),tavg_U1_8,iblock,k)

         call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_VVEL,iblock,k)

         call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_VVEL_2,iblock,k)

         call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock)**2,tavg_VVEL2,iblock,k)

         if (k <= 1)  &
         call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_V1_1,iblock,k)

         if (k <= 8)  &
         call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_V1_8,iblock,k)

         call accumulate_tavg_field(p5*(UVEL(:,:,k,curtime,iblock)**2 + &
                                        VVEL(:,:,k,curtime,iblock)**2),tavg_KE,iblock,k)

         call accumulate_tavg_field(p5*(UVEL(:,:,k,curtime,iblock)**2 + &
                                        VVEL(:,:,k,curtime,iblock)**2),tavg_KE_2,iblock,k)

         call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock)* &
                                    VVEL(:,:,k,curtime,iblock), tavg_UV,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), &
                                    tavg_TEMP,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), &
                                    tavg_TEMP_2,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), &
                                    tavg_TEMP_MAX,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), &
                                    tavg_TEMP_MIN,iblock,k)

         call accumulate_tavg_field(max(TRACER(:,:,k,1,curtime,iblock) - &
                                        TRACER(:,:,k,1,oldtime,iblock), c0), &
                                    tavg_dTEMP_POS_3D,iblock,k)

         call accumulate_tavg_field(max(TRACER(:,:,k,1,curtime,iblock) - &
                                        TRACER(:,:,k,1,oldtime,iblock), c0), &
                                        tavg_dTEMP_POS_2D,iblock,1)

         call accumulate_tavg_field(min(TRACER(:,:,k,1,curtime,iblock) - &
                                        TRACER(:,:,k,1,oldtime,iblock), c0), &
                                    tavg_dTEMP_NEG_3D,iblock,k)

         call accumulate_tavg_field(min(TRACER(:,:,k,1,curtime,iblock) - &
                                        TRACER(:,:,k,1,oldtime,iblock), c0), &
                                        tavg_dTEMP_NEG_2D,iblock,1)

         if (k <= 8)  &
         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), &
                                    tavg_T1_8,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)**2, &
                                    tavg_TEMP2,iblock,k)

         if (k == 1)  then
            call accumulate_tavg_field(TRACER(:,:,1,1,curtime,iblock), &
                                       tavg_SST,iblock,1)

            call accumulate_tavg_field(TRACER(:,:,1,1,curtime,iblock)**2, &
                                       tavg_SST2,iblock,1)
         endif

         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), &
                                    tavg_SALT,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), &
                                    tavg_SALT_2,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), &
                                    tavg_SALT_MAX,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), &
                                    tavg_SALT_MIN,iblock,k)

         if (k <= 8)  &
         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), &
                                    tavg_S1_8,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock)**2, &
                                    tavg_SALT2,iblock,k)

         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)* &
                                    TRACER(:,:,k,2,curtime,iblock), &
                                    tavg_ST,iblock,k)

         call accumulate_tavg_field(RHO(:,:,k,curtime,iblock), &
                                    tavg_RHO,iblock,k)

!         if (partial_bottom_cells) then 
!            WORK1 = RHO(:,:,k,curtime,iblock) * DZT(:,:,k,iblock)
!         else
!            WORK1 = RHO(:,:,k,curtime,iblock) * dz(k)
!         endif
         if (sfc_layer_type == sfc_layer_varthick .and. &
                   k == 1) then
           where (k <= KMT(:,:,iblock))
             WORK1 = RHO(:,:,k,curtime,iblock) * ( dz(k) + PSURF(:,:,curtime,iblock)/grav )
           elsewhere
             WORK1 = c0
           endwhere
         else
           if (partial_bottom_cells) then
             where (k <= KMT(:,:,iblock))
               WORK1 = RHO(:,:,k,curtime,iblock) * DZT(:,:,k,iblock)
             elsewhere
               WORK1 = c0
             endwhere
           else
             where (k <= KMT(:,:,iblock))
               WORK1 = RHO(:,:,k,curtime,iblock) * dz(k)
             elsewhere
               WORK1 = c0
             endwhere
           endif
         endif 
         call accumulate_tavg_field(WORK1,tavg_RHO_VINT,iblock,k)

        if ( sfc_layer_type /= sfc_layer_varthick .and. k == 1) then
          if (accumulate_tavg_now(tavg_RESID_T)) then
              WORK1 = c0
              factor = c1/hflux_factor  ! converts to W/m^2
              where (CALCT(:,:,iblock))  &
                WORK1=DH(:,:,iblock)*TRACER(:,:,1,1,curtime,iblock)*factor
              call accumulate_tavg_field(WORK1,tavg_RESID_T,iblock,k)
          endif

          if (accumulate_tavg_now(tavg_RESID_S)) then
              WORK1 = c0
              factor = c1/salinity_factor  ! converts to kg(freshwater)/m^2/s
              where (CALCT(:,:,iblock)) &
                WORK1 = DH(:,:,iblock)*TRACER(:,:,k,2,curtime,iblock)*factor
              call accumulate_tavg_field(WORK1,tavg_RESID_S,iblock,k)
          endif
        endif  ! sfc_layer_type



! added by zengyh 2019-02-15:
!----------------------------
#endif


#if defined(POPCHKJN)
#if defined(POPTIMEJN)
         call getlooptime(stop_accum)
         if(my_task==0 .and. mod(numsteps,64)==0 .and. k<=1) then  
            WRITE (*,*) 'original accum_tavg of clinic_driver costs zengyh:',stop_accum-start_accum,'(us);'
         endif
#endif


         do i=1,36
           if(tavgdims(i)==3) then
             tmptavg3d(:,:,k,i) = TAVG_BUF_3D(:,:,k,iblock,tavglocs(i))
             TAVG_BUF_3D(:,:,k,iblock,tavglocs(i)) = initavg3d(:,:,k,i)
           elseif(tavgdims(i)==2 .and. k==1) then
             tmptavg3d(:,:,k,i) = TAVG_BUF_2D(:,:,iblock,tavglocs(i))
             TAVG_BUF_2D(:,:,iblock,tavglocs(i)) = initavg3d(:,:,k,i)
           endif
         enddo
#endif

! end_added by zengyh 2019-02-15.
!----------------------------



         if (nt > 2) call tavg_passive_tracers(iblock,k)

   !call timer_stop(driver_accum_tavg)
!-----------------------------------------------------------------------
!
!        update movie fields if requested
!
!-----------------------------------------------------------------------

   !call timer_start(driver_up_movie)
         if (movie_requested(movie_UVEL(k))) then
            call update_movie_field(UVEL(:,:,k,curtime,iblock), &
                                       movie_UVEL(k),iblock,k)
         endif

         if (movie_requested(movie_VVEL(k))) then
            call update_movie_field(VVEL(:,:,k,curtime,iblock), &
                                       movie_VVEL(k),iblock,k)
         endif

         if (movie_requested(movie_TEMP(k))) then
            call update_movie_field(TRACER(:,:,k,1,curtime,iblock), &
                                       movie_TEMP(k),iblock,k)
         endif

         if (movie_requested(movie_SALT(k))) then
            call update_movie_field(  &
                 TRACER(:,:,k,2,curtime,iblock)*salt_to_ppt, &
                                       movie_SALT(k),iblock,k)  !  convert to psu
         endif

         if (movie_requested(movie_RHO(k))) then
            call update_movie_field(  &
               salt_to_ppt*(RHO(:,:,k,curtime,iblock)-c1), movie_RHO(k),iblock,k)

         endif

         endif ! mix_pass

!-----------------------------------------------------------------------

   !call timer_stop(driver_up_movie)
      enddo  ! k loop




! added by zengyh 2019-02-15:
!----------------------------

#if defined(POPTIMEJN)
      numsteps = numsteps + 1
#endif


#if !defined(NOATH_TAVG_CLINIC_ZYH)

      if (mix_pass /= 1) then

#if defined(POPTIMEJN)
         call getlooptime(start_accum)
#endif

         call m_clinic_tavg_sw(TAVG_BUF_2D,TAVG_BUF_3D, &
            UVEL(:,:,1:km,curtime,iblock),VVEL(:,:,1:km,curtime,iblock), &
            TRACER(:,:,1:km,1,curtime,iblock),TRACER(:,:,1:km,1,oldtime,iblock),TRACER(:,:,1:km,2,curtime,iblock), &
            RHO(:,:,1:km,curtime,iblock),PSURF(:,:,curtime,iblock),DZT(:,:,1:km,iblock),DH(:,:,iblock), &
            dtavg,dz,grav,hflux_factor,salinity_factor, &
            KMT(:,:,iblock),intCALCT(:,:), &
            tavgids(:),tavgflags(:),tavglocs(:),tavgdims(:),tavgmeths(:), &
            my_task,km,nx_block,ny_block,iblock,nblocks_clinic,&
            sfc_layer_type,sfc_layer_varthick,numtavg,part_bott_cells)

#if defined(POPTIMEJN)
         call getlooptime(stop_accum)
         if(my_task==0 .and. mod(numsteps,64)==0) then  
            WRITE (*,*) 'athreaded accum_tavg of clinic_driver costs zengyh:',stop_accum-start_accum,'(us);'
         endif
#endif


#if defined(POPCHKJN)
         do i=1,36
          if(tavgdims(i)==3) then
             mykm = km
          elseif(tavgdims(i)==2) then
             mykm = 1
          endif
          do k=1,mykm
           if(tavgids(i)>0 .and. tavglocs(i)>0) then
!           if(tavgflags((i-1)*km+k)>0) then
! otherwise, tavgid(i) may be null.
             cstat(i,:) = 0
             if(tavgdims(i)==3) then
               initavg3d(:,:,k,i) = abs(tmptavg3d(:,:,k,i) - TAVG_BUF_3D(:,:,k,iblock,tavglocs(i)))
             elseif(tavgdims(i)==2) then
               initavg3d(:,:,k,i) = abs(tmptavg3d(:,:,k,i) - TAVG_BUF_2D(:,:,iblock,tavglocs(i)))
             endif
             tmprtavg3d(:,:,k,i)= merge(initavg3d(:,:,k,i)/abs(tmptavg3d(:,:,k,i)),initavg3d(:,:,k,i),abs(tmptavg3d(:,:,k,i))>f0)
             maxerrtavg3d(i) = f0
             maxerrtavg3d(i) = maxval(initavg3d(:,:,k,i))
             maxrerrtavg3d(i) = maxval(tmprtavg3d(:,:,k,i))
             minerrtavg3d(i) = f1
             minerrtavg3d(i) = minval(initavg3d(:,:,k,i))
             errtavg3d(i) = sum(initavg3d(:,:,k,i)/(nx_block*ny_block))
             rerrtavg3d(i) = sum(tmprtavg3d(:,:,k,i)/(nx_block*ny_block))
             if(my_task==0 .and. (rerrtavg3d(i)>=1.0e-6 .or. maxrerrtavg3d(i)>=1.0e-4)) then
               write(*,*) 'rerr-6-4_tavg: ',i, tavgids(i), ', k: ', k, ', flag: ', tavgflags(i), ', meth: ', tavgmeths(i), '; error: ', maxerrtavg3d(i),maxrerrtavg3d(i),errtavg3d(i),rerrtavg3d(i),minerrtavg3d(i)
             elseif(my_task==0 .and. (rerrtavg3d(i)>=1.0e-8 .or. maxrerrtavg3d(i)>=1.0e-6)) then
               write(*,*) 'rerr-8-6_tavg: ',i, tavgids(i), ', k: ', k, ', flag: ', tavgflags(i), ', meth: ', tavgmeths(i), '; error: ', maxerrtavg3d(i),maxrerrtavg3d(i),errtavg3d(i),rerrtavg3d(i),minerrtavg3d(i)
             elseif(my_task==0 .and. (rerrtavg3d(i)>=1.0e-10 .or. maxrerrtavg3d(i)>=1.0e-8)) then
               write(*,*) 'rerr-10-8_tavg: ',i, tavgids(i), ', k: ', k, ', flag: ', tavgflags(i), ', meth: ', tavgmeths(i), '; error: ', maxerrtavg3d(i),maxrerrtavg3d(i),errtavg3d(i),rerrtavg3d(i),minerrtavg3d(i)
             endif
#ifdef POPTESTJN 
             do j=1,ny_block
             do itmp=1,nx_block
               if(tmprtavg3d(itmp,j,k,i)>=1.0e-2) then
                  cstat(i,1) = cstat(i,1) + 1
               else if(tmprtavg3d(itmp,j,k,i)>=1.0e-4) then
                  cstat(i,2) = cstat(i,2) + 1
               else if(tmprtavg3d(itmp,j,k,i)>=1.0e-6) then
                  cstat(i,3) = cstat(i,3) + 1
               else if(tmprtavg3d(itmp,j,k,i)>=1.0e-8) then
                  cstat(i,4) = cstat(i,4) + 1
               else
                  cstat(i,5) = cstat(i,5) + 1
               endif
             enddo
             enddo
             if(my_task==0 .and. cstat(i,5)<(nx_block*ny_block)) then
               write(*,*) 'tavg_stat: ',i, tavgids(i), ', k: ', k, ', >=E-2: ', cstat(i,1), ', >=E-4: ', cstat(i,2), ', >=E-6: ', cstat(i,3), ', >=E-8: ', cstat(i,4), ', <E-8: ', cstat(i,5)
             endif

             if(my_task==0) then
               if(tavgdims(i)==3) then
                 initavg3d(:,:,k,i) = TAVG_BUF_3D(:,:,k,iblock,tavglocs(i))
               elseif(tavgdims(i)==2) then
                 initavg3d(:,:,k,i) = TAVG_BUF_2D(:,:,iblock,tavglocs(i))
               endif
               do j=1,ny_block
               do itmp=1,nx_block
                 if(tmprtavg3d(itmp,j,k,i)>1.0e-4) then
                   write(*,*) 'e-4(i,j) ',i, tavgids(i), ', k: ', k, ', j: ', j, ', i: ', itmp, '; m_value: ', tmptavg3d(itmp,j,k,i), ', s_value: ', initavg3d(itmp,j,k,i), ', rel_error: ', tmprtavg3d(itmp,j,k,i)
                 else if(tmprtavg3d(itmp,j,k,i)>1.0e-6) then
                   write(*,*) 'e-6(i,j) ',i, tavgids(i), ', k: ', k, ', j: ', j, ', i: ', itmp, '; m_value: ', tmptavg3d(itmp,j,k,i), ', s_value: ', initavg3d(itmp,j,k,i), ', rel_error: ',tmprtavg3d(itmp,j,k,i)
                 else if(tmprtavg3d(itmp,j,k,i)>1.0e-8) then
                   write(*,*) 'e-8(i,j) ',i, tavgids(i), ', k: ', k, ', j: ', j, ', i: ', itmp, '; m_value: ', tmptavg3d(itmp,j,k,i), ', s_value: ', initavg3d(itmp,j,k,i), ', rel_error: ',tmprtavg3d(itmp,j,k,i)
                 endif
               enddo
               enddo
             endif
#endif
           endif   ! end of tavgflags
         enddo  ! end of do k
         enddo  ! end of do i
#endif


      endif ! mix_pass
!----------------------------
#endif

! end_added by zengyh 2019-02-15.
!----------------------------




!-----------------------------------------------------------------------
!
!     if no implicit vertical mixing, we now have updated tracers
!     if using implicit vertical mixing and rigid lid or old free
!        surface form, update all the tracers now using impvmix
!     if using implicit vertical mixing and a variable thickness
!        surface layer, only update T,S as predicted values to
!        use for pressure averaging - the full update will
!        occur after the barotropic solver
!
!-----------------------------------------------------------------------

      if (implicit_vertical_mix) then
         if (sfc_layer_type /= sfc_layer_varthick) then
            call impvmixt(TRACER(:,:,:,:,newtime,iblock), &
                          TRACER(:,:,:,:,oldtime,iblock), &
                          PSURF (:,:,    curtime,iblock), &
                          1, nt, this_block)

         else if (lpressure_avg .and. leapfrogts) then

            !*** predictor update of T,S
            !*** with PSURF(curtime) on the LHS at k=1
 
            call impvmixt(TRACER(:,:,:,:,newtime,iblock), &
                          TRACER(:,:,:,:,oldtime,iblock), &
                          PSURF (:,:,curtime,iblock),     &
                          1, 2, this_block)

         endif
      endif

!-----------------------------------------------------------------------
!
!     end of first block loop
!
!-----------------------------------------------------------------------

   enddo ! first block loop

   !$OMP END PARALLEL DO


!-----------------------------------------------------------------------
!
!  update tracer ghost cells here outside the block loops (it
!  requires synchronization) - only update T,S ghost cells
!  for use in pressure averaging.
!
!-----------------------------------------------------------------------

   !call timer_start(timer_3dupdclinic)

   if (lpressure_avg .and. leapfrogts) then

      call POP_HaloUpdate(TRACER(:,:,:,1,newtime,:),          &
                              POP_haloClinic,                 &
                              POP_gridHorzLocCenter,          &
                              POP_fieldKindScalar, errorCode, &
                              fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'baroclinic_driver: error updating halo for PT')
         return
      endif

      call POP_HaloUpdate(TRACER(:,:,:,2,newtime,:),          &
                              POP_haloClinic,                 &
                              POP_gridHorzLocCenter,          &
                              POP_fieldKindScalar, errorCode, &
                              fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'baroclinic_driver: error updating halo for salinity')
         return
      endif

   endif

   !call timer_stop(timer_3dupdclinic)

!-----------------------------------------------------------------------
!
!  now loop over blocks to do momentum equations
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,km1,kp1,n, &
   !$OMP                     WUK,FX,FY,WORK1,WORK2)

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)  

!-----------------------------------------------------------------------
!
!     initialize arrays for vertical sums.
!
!-----------------------------------------------------------------------

      ZX(:,:,iblock) = c0
      ZY(:,:,iblock) = c0

      do k = 1,km

         kp1 = k+1
         km1 = k-1
         if (k == 1) km1 = 1
         if (k == km) kp1 = km

!-----------------------------------------------------------------------
!
!        if pressure averaging is on and it is a leapfrog time step
!        we need the updated density for the pressure averaging
!
!-----------------------------------------------------------------------

        !call timer_start(clinic_state1)
         if (lpressure_avg .and. leapfrogts) then

!!! added_start by wangli
#ifdef ZYUATH
        call loop2_state_baroclinic(nx_block,ny_block,km,TRACER(1,1,k,1,newtime,iblock), &
                                    TRACER(1,1,k,2,newtime,iblock),&
                       RHO(1,1,k,newtime,iblock),tmax(k),tmin(k),smax(k),smin(k),pressz(k))
#else
!!! added_end by wangli

            call state(k,k,TRACER(:,:,k,1,newtime,iblock), &
                           TRACER(:,:,k,2,newtime,iblock), &
                           this_block, RHOOUT=RHO(:,:,k,newtime,iblock))

!!! added_start by wangli
#endif
!!! added_end by wangli

         endif
        !call timer_stop(clinic_state1)

!-----------------------------------------------------------------------
!
!        calculate forcing terms (Fx,Fy) at level k.
!
!-----------------------------------------------------------------------

        !call timer_start(clinic_s)
         call clinic(k, FX, FY, WUK,                &
                        UVEL(:,:,:,curtime,iblock), &
                        VVEL(:,:,:,curtime,iblock), &
                        UVEL(:,:,:,oldtime,iblock), &
                        VVEL(:,:,:,oldtime,iblock), &
                        UVEL(:,:,k,mixtime,iblock), &
                        VVEL(:,:,k,mixtime,iblock), &
                        RHO (:,:,k,oldtime,iblock), &
                        RHO (:,:,k,curtime,iblock), &
                        RHO (:,:,k,newtime,iblock), &
                        SMF (:,:,:,iblock),         &
                        DHU (:,:,iblock),           &
                        this_block)

        !call timer_stop(clinic_s)
!-----------------------------------------------------------------------
!
!        store forces temporarily in UVEL(newtime),VVEL(newtime).
!
!-----------------------------------------------------------------------

        !call timer_start(clinic_storeforce)

#ifndef NOWL_ATH_ST
        iparam(1)=nx_block
        iparam(2)=ny_block
        iparam(3)=k
        iparam(4)=0     !impcor
        iparam(5)=0     !l1ddyn
        iparam(6)=0     !partial
        if(impcor) iparam(4)=1
        if(l1ddyn) iparam(5)=1
        if (partial_bottom_cells) iparam(6)=1
        rparam(1)=c2dtu
        rparam(2)=beta
        rparam(3)=c1
        rparam(4)=dz(k)
        call storeforce_c(iparam,rparam,fcor(:,:,iblock),fx,fy,&
                          uvel(:,:,:,newtime,iblock),vvel(:,:,:,newtime,iblock),&
                          zx(:,:,iblock),zy(:,:,iblock),dzu(:,:,:,iblock))

#else
         if (impcor) then   ! implicit treatment
            WORK1 = c2dtu*beta*FCOR(:,:,iblock)
            WORK2 = c2dtu/(c1 + WORK1**2)
            UVEL(:,:,k,newtime,iblock) = (FX + WORK1*FY)*WORK2 
            VVEL(:,:,k,newtime,iblock) = (FY - WORK1*FX)*WORK2 
         else               ! explicit treatment
            UVEL(:,:,k,newtime,iblock) = c2dtu*FX
            VVEL(:,:,k,newtime,iblock) = c2dtu*FY
         endif
!-----------------------------------------------------------------------
!
!        increment sum for vertically-averaged forcing ([Fx],[Fy]).
!
!-----------------------------------------------------------------------
         if (.not.l1Ddyn) then
           if (partial_bottom_cells) then
              ZX(:,:,iblock) = ZX(:,:,iblock) + FX*DZU(:,:,k,iblock)
              ZY(:,:,iblock) = ZY(:,:,iblock) + FY*DZU(:,:,k,iblock)
           else
              ZX(:,:,iblock) = ZX(:,:,iblock) + FX*dz(k)
              ZY(:,:,iblock) = ZY(:,:,iblock) + FY*dz(k)
           endif
        endif
#endif

        !call timer_stop(clinic_storeforce)
      enddo ! vertical (k) loop

!-----------------------------------------------------------------------
!
!     normalize sums for vertical averages ([Fx],[Fy]) by dividing
!     by depth at U points.
!
!-----------------------------------------------------------------------

      ZX(:,:,iblock) = ZX(:,:,iblock)*HUR(:,:,iblock)
      ZY(:,:,iblock) = ZY(:,:,iblock)*HUR(:,:,iblock)

!-----------------------------------------------------------------------
!
!     solve tridiagonal system with implicit treatment of vertical 
!     diffusion of velocity.
!
!-----------------------------------------------------------------------

      if (implicit_vertical_mix)                   &
         call impvmixu(UVEL(:,:,:,newtime,iblock), &
                       VVEL(:,:,:,newtime,iblock), & 
                       UVEL(:,:,:,oldtime,iblock), & 
                       VVEL(:,:,:,oldtime,iblock), & 
                       this_block)

!-----------------------------------------------------------------------
!
!     calculate unnormalized baroclinic velocities (Upp,Vpp)
!
!-----------------------------------------------------------------------

!   !call timer_start(timer_sec1clinic)
!      UVEL(:,:,:,newtime,iblock) = UVEL(:,:,:,oldtime,iblock) + &
!                                   UVEL(:,:,:,newtime,iblock)  ! holds c2dtu*Fx
!      VVEL(:,:,:,newtime,iblock) = VVEL(:,:,:,oldtime,iblock) + &
!                                   VVEL(:,:,:,newtime,iblock)  ! holds c2dtu*Fy
!   !call timer_stop(timer_sec1clinic)

      if ( overflows_on .and. overflows_interactive ) then
         call ovf_Utlda(iblock)
      endif

!-----------------------------------------------------------------------
!
!     find vertical averages ([Upp],[Vpp]).
!
!-----------------------------------------------------------------------
   !call timer_start(timer_sec2clinic)

!!! added-start by wangli
#ifdef WLATH
   if(partial_bottom_cells) then
     partial=1
   else
     partial=0
   endif
      call clinic_vertav_normalv_c(nx_block,ny_block,km,partial,c0,&
                                   uvel(:,:,:,newtime,iblock),&
                                   vvel(:,:,:,newtime,iblock),dzu(:,:,1,iblock),&
                                   hur(:,:,iblock),dz,kmu(:,:,iblock),l1Ddyn)
#else
!!! added-end by wangli

      if (.not.l1Ddyn) then

        WORK1 = c0  ! initialize sums
        WORK2 = c0

        if (partial_bottom_cells) then
           do k = 1,km
              WORK1 = WORK1 + UVEL(:,:,k,newtime,iblock)*DZU(:,:,k,iblock)
              WORK2 = WORK2 + VVEL(:,:,k,newtime,iblock)*DZU(:,:,k,iblock)
           enddo
        else
           do k = 1,km
              WORK1 = WORK1 + UVEL(:,:,k,newtime,iblock)*dz(k)
              WORK2 = WORK2 + VVEL(:,:,k,newtime,iblock)*dz(k)
           enddo
        endif

        WORK1 = WORK1*HUR(:,:,iblock)  ! normalize by dividing by depth
        WORK2 = WORK2*HUR(:,:,iblock)

!-----------------------------------------------------------------------
!
!     normalize baroclinic velocites by subtracting vertical mean:
!     (Up,Vp) = (Upp,Vpp) - ([Upp],[Vpp]), zero velocities at land pts.
!
!-----------------------------------------------------------------------

        do k = 1,km
           where (k <= KMU(:,:,iblock))
              UVEL(:,:,k,newtime,iblock) = &
              UVEL(:,:,k,newtime,iblock) - WORK1
              VVEL(:,:,k,newtime,iblock) = &
              VVEL(:,:,k,newtime,iblock) - WORK2
           elsewhere 
              UVEL(:,:,k,newtime,iblock) = c0
              VVEL(:,:,k,newtime,iblock) = c0
           endwhere
        enddo

      else

        do k = 1,km
          where (k > KMU(:,:,iblock))
            UVEL(:,:,k,newtime,iblock) = c0
            VVEL(:,:,k,newtime,iblock) = c0
          endwhere
        enddo

      endif

!!! added-start by wangli
#endif
!!! added-end by wangli

   !call timer_stop(timer_sec2clinic)
!-----------------------------------------------------------------------
!
!     note:  at this point UVEL(newtime) and VVEL(newtime) contain only 
!     the baroclinic velocities (Up,Vp) at the new time.  they are later
!     updated to the full velocities in step after barotropic is
!     is called, which calculates the barotropic velocites ([U],[V])
!     at the new time.  UVEL and VVEL at time levels oldtime ond curtime 
!     always contain the full velocites.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  end of second block loop over vertical levels.
!
!-----------------------------------------------------------------------

   enddo ! second block loop
   !$OMP END PARALLEL DO

#if drifter_particles
!-----------------------------------------------------------------------
!
!  advect drifters with old time level velocities.  then if its
!  time to sample ocean variables at drifter locations, do it.
!
!-----------------------------------------------------------------------

   if (mod(iday,10).eq.0 .and. newday) then
      if((ndrifters + array_size) .gt. ndrifters_total) then
         write(stdout,*)' '
         write(stdout,*) &
                 ' No more drifter arrays remaining for deployment'
         write(stdout,*)' '
      else
         arrays_deployed = arrays_deployed + 1
         ndrifters = ndrifters + array_size
         write(stdout,*)' '
         write(stdout,*)' Deploying drifter array # ',arrays_deployed
         write(stdout,*)' '
      endif
   endif

   call drifter_move
   if(mod(ihour,n_write_drifters).eq.0 .and. newhour) then
      call drifter_prop
   endif
#endif

!-----------------------------------------------------------------------
!
!  check cfl limits
!
!-----------------------------------------------------------------------

   if (ldiag_cfl) then
      call cfl_check
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine baroclinic_driver

!***********************************************************************
!BOP
! !IROUTINE: baroclinic_correct_adjust
! !INTERFACE:

 subroutine baroclinic_correct_adjust

! !DESCRIPTION:
!  This subroutine finishes updating tracers by performing
!  adjustment-like physics (convection and ice) and completes
!  the corrector step for tracers using the new surface pressure
!  to update tracers.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  & 
      k,                  &! vertical level index
      n,                  &! tracer index
      iblock               ! block index

   real (r8), dimension(nx_block,ny_block,nt) :: &
      RHS1                 ! r.h.s. for impvmix on corrector step

   type (block) ::       &
      this_block           ! block information for current block

!-----------------------------------------------------------------------
!
!  increment tlast_ice outside of threaded region
!
!-----------------------------------------------------------------------

   call increment_tlast_ice

!-----------------------------------------------------------------------
!
!  do everything for each sub block
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block,n,RHS1)

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)  

!-----------------------------------------------------------------------
!
!     do the corrector step for variable thickness surface layer
!
!-----------------------------------------------------------------------

      if (sfc_layer_type == sfc_layer_varthick) then

         if (implicit_vertical_mix) then

            !*** if implicit vertical mixing and pressure averaging:
            !*** correct new T and S and update remaining passive
            !*** tracers with tridiagonal solves

            if (lpressure_avg .and. leapfrogts) then
        
               do n = 1,2   ! corrector for T,S only
                  where (KMT(:,:,iblock) > 0)  ! corrector RHS at k=1
                     RHS1(:,:,n)=((c2*TRACER(:,:,1,n,curtime,iblock) - &
                                      TRACER(:,:,1,n,oldtime,iblock))  &
                                 *(PSURF(:,:,curtime,iblock) -         &
                                   PSURF(:,:,oldtime,iblock))          &
                                - TRACER(:,:,1,n,newtime,iblock)*      &
                                  (PSURF(:,:,newtime,iblock) -         &
                                   PSURF(:,:,curtime,iblock)))/        &
                                   (grav*dz(1)) 
                  elsewhere
                     RHS1(:,:,n) = c0
                  endwhere
               enddo

               !*** T,S update on corrector step
               call impvmixt_correct(TRACER(:,:,:,:,newtime,iblock), &
                                     PSURF (:,:,    newtime,iblock), &
                                     RHS1, 1, 2, this_block)

               do n = 3,nt
                  !*** surface RHS for passive tracers with press avg
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) =               &
                                   TRACER(:,:,1,n,newtime,iblock) - &
                                   TRACER(:,:,1,n,oldtime,iblock)   & 
                                   *(PSURF(:,:,newtime,iblock) -    &
                                     PSURF(:,:,oldtime,iblock))/    &
                                     (grav*dz(1))
                  endwhere
               enddo

               !*** standard update of all passive tracers
               !***  n=3,nt with PSURF(newtime) on LHS at k=1

               call impvmixt(TRACER(:,:,:,:,newtime,iblock), &
                             TRACER(:,:,:,:,oldtime,iblock), &
                             PSURF (:,:,    newtime,iblock), &
                             3, nt, this_block)

            !*** if implicit vertical mixing but no pressure averaging
            !*** update all tracers with tridiagonal solves

            else  ! no pressure averaging or not leapfrog

               do n = 1,nt
                  !***  surface RHS for tracers with pressure avg
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) = &
                     TRACER(:,:,1,n,newtime,iblock) - &
                     TRACER(:,:,1,n,oldtime,iblock)   & 
                     *(PSURF(:,:,newtime,iblock) -    &
                       PSURF(:,:,mixtime,iblock))/(grav*dz(1))
                  endwhere
               enddo

               !*** standard update all tracers:  n=1,nt
               !*** with PSURF(newtime) on the LHS at k=1
               call impvmixt(TRACER(:,:,:,:,newtime,iblock), &
                             TRACER(:,:,:,:,oldtime,iblock), &
                             PSURF (:,:,    newtime,iblock), &
                             1, nt, this_block)

            endif  ! pressure averaging and leapfrog

         else ! no implicit_vertical_mix

            !*** if explicit vertical mixing and pressure averaging:
            !*** correct new tracers level k=1

            if (lpressure_avg .and. leapfrogts) then

               do n = 1,2
                  !*** correct surface T and S with pressure avg
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) =           &
                       (TRACER(:,:,1,n,newtime,iblock)*         &
                       (dz(1) + PSURF(:,:,curtime,iblock)/grav) &
                       +(c2*TRACER(:,:,1,n,curtime,iblock) -    &
                            TRACER(:,:,1,n,oldtime,iblock))     &
                       *(PSURF(:,:,curtime,iblock) -            &
                         PSURF(:,:,oldtime,iblock))/grav)       &
                       /(dz(1) + PSURF(:,:,newtime,iblock)/grav) 
                  elsewhere
                     TRACER(:,:,1,n,newtime,iblock) = c0 ! zero on land pts
                  endwhere
               enddo

               do n = 3,nt
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) =             &
                        (TRACER(:,:,1,n,oldtime,iblock)*          &
                        (dz(1) + PSURF(:,:,oldtime,iblock)/grav)  &
                        + dz(1)*TRACER(:,:,1,n,newtime,iblock))   &
                        /(dz(1) + PSURF(:,:,newtime,iblock)/grav)
                  elsewhere
                     TRACER(:,:,1,n,newtime,iblock) = c0 ! zero on land pts
                  endwhere
               enddo

            else  ! no pressure avg or leapfrog

               do n = 1,nt

                  !*** exact update of all tracers in surface layer
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) =              &
                        (TRACER(:,:,1,n,oldtime,iblock)*           &
                        (dz(1) + PSURF(:,:,mixtime,iblock)/grav)   &
                        + dz(1)*TRACER(:,:,1,n,newtime,iblock))/   &
                        (dz(1) + PSURF(:,:,newtime,iblock)/grav)
                  elsewhere
                     TRACER(:,:,1,n,newtime,iblock) = c0 ! zero on land pts
                  endwhere

               enddo

            endif  ! pressure avg and leapfrog

         endif ! implicit_vertical_mix

      endif ! variable thickness surface layer
        !call timer_start(impt_tavg)
      if (mix_pass /= 1) then
         call impvmixt_tavg(TRACER(:,:,:,:,newtime,iblock), iblock)
         call iso_impvmixt_tavg(TRACER(:,:,:,:,newtime,iblock), iblock)
      endif

        !call timer_stop(impt_tavg)
!-----------------------------------------------------------------------
!
!     check for surface temperatures below freezing
!     do not reset if ice formation option is on
!
!-----------------------------------------------------------------------

      if (reset_to_freezing .and. .not. liceform) then
         TRACER(:,:,1,1,newtime,iblock) = &
            max(TRACER(:,:,1,1,newtime,iblock),-c2)
      endif

!-----------------------------------------------------------------------
!
!     convective adjustment of tracers - 
!     convad routine does nothing if convective adjustment not chosen
!     otherwise it performs convective adjustment and recomputes
!     density
!
!-----------------------------------------------------------------------

        !call timer_start(tm_convad)
      call convad(TRACER(:,:,:,:,newtime,iblock), &
                  RHO(:,:,:,newtime,iblock), this_block, iblock)
        !call timer_stop(tm_convad)

!-----------------------------------------------------------------------
!
!     compute ice formation and adjust temperature due to ice formation
!     if this option was requested
!
!-----------------------------------------------------------------------

        !call timer_start(posttracer)
      if (liceform .and. mix_pass /= 1) then
         call ice_formation(TRACER(:,:,:,:,newtime,iblock),          &
                            STF(:,:,1,iblock) + SHF_QSW(:,:,iblock), &
                            iblock,this_block,lfw_as_salt_flx)
      endif

      if (mix_pass /= 1) then
         if (nt > 2) call tavg_passive_tracers_baroclinic_correct(iblock)
      endif

!-----------------------------------------------------------------------
!
!     call passive tracer reset subroutines
!
!-----------------------------------------------------------------------

      if (nt > 2) call reset_passive_tracers(  &
         TRACER(:,:,:,:,newtime,iblock), iblock)


        !call timer_stop(posttracer)
!-----------------------------------------------------------------------
!
!     compute new density based on new tracers
!
!-----------------------------------------------------------------------

        !call timer_start(state2)

!!! added-start by wangli
#ifdef ZYUATH
     call loop_state_baroclinic(nx_block,ny_block,km,TRACER(:,:,:,:,newtime,iblock), &
                      RHO(:,:,:,newtime,iblock),tmax,tmin,smax,smin,pressz)
#else
!!! added-end by wangli

      do k = 1,km  ! recalculate new density

         call state(k,k,TRACER(:,:,k,1,newtime,iblock), &
                        TRACER(:,:,k,2,newtime,iblock), &
                        this_block, RHOOUT=RHO(:,:,k,newtime,iblock))

      enddo

!!! added-start by wangli
#endif
!!! added-end by wangli

        !call timer_stop(state2)
!-----------------------------------------------------------------------
!
!  end of block loop
!
!-----------------------------------------------------------------------

   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine baroclinic_correct_adjust

!***********************************************************************
!BOP
! !IROUTINE: clinic
! !INTERFACE:

 subroutine clinic(k, FX, FY, WUK, UCUR, VCUR, UOLD, VOLD,      &
                      UMIXK, VMIXK, RHOKOLD, RHOKCUR, RHOKNEW,  &
                      SMF_BLOCK, DHU_BLOCK, this_block)

! !DESCRIPTION:
!  Calculates forcing terms on r.h.s. of baroclinic momentum eqns.
!  \begin{eqnarray}
!     F_x &=& -L(u) + fv - \nabla p + D_H(u^{n-1}) + D_V(u^{n-1}) \\
!     F_y &=& -L(v) - fu - \nabla p + D_H(v^{n-1}) + D_V(v^{n-1})
!  \end{eqnarray}
!
!  The above equations are written for the case of explicit
!  treatment of the Coriolis terms.  If these terms are treated
!  implicitly, then the coriolis terms above should be replaced by:
!  \begin{eqnarray}
!       +fv &\rightarrow& +f(\gamma v + (1-\gamma)v^{n-1}) \\
!       -fu &\rightarrow& -f(\gamma u + (1-\gamma)u^{n-1})
!  \end{eqnarray}
!  on leapfrog timesteps and
!  \begin{eqnarray}
!       +fv &\rightarrow& +fv^{n-1} \\
!       -fu &\rightarrow& -fu^{n-1}
!  \end{eqnarray}
!  on Matsuno timesteps, where $\gamma$ is a parameter used to vary 
!  the time-centering of the Coriolis and pressure gradient terms on 
!  leapfrog steps.
!
!  The small metric terms for advection and diffusion of the
!  velocity field are calculated in the advection and horizontal 
!  diffusion routines.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k    ! depth level index

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UCUR, VCUR,           &! U,V for block at current time
      UOLD, VOLD             ! U,V for block at old     time

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      UMIXK, VMIXK,         &! U,V at level k and mix time level
      RHOKOLD,              &! density at level k and mix time level
      RHOKCUR,              &! density at level k and cur time level
      RHOKNEW,              &! density at level k and new time level
      DHU_BLOCK              ! change in surface height at U pts

   real (r8), dimension(nx_block,ny_block,2), intent(in) :: &
      SMF_BLOCK              ! surface momentum forcing for this block

   type (block), intent(in) :: &
      this_block             ! block info for the current block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(inout) :: &
      WUK             ! vertical velocity at top of U box

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      FX,       &! sum of terms contributing to Fx at level k 
      FY         ! sum of terms contributing to Fy at level k 

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bid                ! local block address

  real (r8), dimension(nx_block,ny_block) :: &
      WORKX,WORKY    &     ! local work space to hold forcing terms
      ,WORKX1,WORKY1 &     ! local work space to hold hdiiff.
      ,UUW, UUE, VUS, VUN, WUKB & ! local space for tmp varibles.
      ,CC,CN,CS,CE,CW & 
      ,D2UK, D2VK

#if defined(adhdu_chk)
   real (r8), dimension(nx_block,ny_block) :: &
      WUK_ini,WUK_s,     &   ! local work space to hold master forcing terms
      WORKX_s,WORKY_s    &   ! local work space to hold master forcing terms
      ,WORKX1_s,WORKY1_s     ! local work space to hold master hdiiff.
   real (r8) :: erreps         

   integer (int_kind) :: i,j
   integer (int_kind) :: advux_flg,advuy_flg,advuk_flg
   integer (int_kind) :: hdiffux_flg,hdiffuy_flg
   integer (int_kind), save :: numsteps=0
#endif
   real (r8) :: factor  ! temporary factor

!!! added-start by wangli, to expand gradp/grad
#ifndef WLATH
   real (r8), dimension(nx_block,ny_block) :: &
      RHOKX, RHOKY,      &! x,y gradients of level k density
      RHOAVG              ! avg density when pressure avg on
   integer (int_kind)::  i,j
#else
   integer (int_kind)::  iparam(16)
   real (r8) ::  rparam(12)
#endif
!!! added-end by wangli

!-----------------------------------------------------------------------
!
!  advection L(U),L(V)
!  set vertical velocity at surface
!
!-----------------------------------------------------------------------
        !call timer_start(call_advu)
   bid = this_block%local_id
 
   if (k == 1) WUK = DHU_BLOCK  ! free surface

   ! new combined advu and hdiffu:
!
   if (.not.l1Ddyn) then
#if defined(adhdu_chk)
     WUK_ini = WUK
#endif

#if !defined(NOATH_adhdu_GT)
     call advuhdiffu(k, WORKX, WORKY, WORKX1, WORKY1, WUK, UCUR, VCUR, UMIXK, VMIXK, this_block)
#endif


#if defined(adhdu_chk)
     WORKX_s  = WORKX
     WORKY_s  = WORKY
     WUK_s    = WUK
     WORKX1_s = WORKX1
     WORKY1_s = WORKY1
     WUK = WUK_ini
#endif




#if defined(adhdu_chk) || defined(NOATH_adhdu_GT) 
! original calls:
     call advu(k, WORKX, WORKY, WUK, UCUR, VCUR, this_block)
     call hdiffu(k, WORKX1, WORKY1, UMIXK, VMIXK, this_block)
#endif


#if defined(adhdu_chk)
    if(my_task==0 .and. numsteps<=768) then
      advux_flg = 1
      advuy_flg = 1
      advuk_flg = 1
      hdiffux_flg = 1
      hdiffuy_flg = 1
      erreps = 1.0E-10
      do j=1,ny_block
      do i=1,nx_block
        if((ABS(WORKX_s(i,j)-WORKX(i,j))/WORKX(i,j))>erreps) then
           print*, "advu_diff_LUK=", ABS(WORKX_s(i,j)-WORKX(i,j))/WORKX(i,j),WORKX_s(i,j),WORKX(i,j),i,j,k
           advux_flg = 0
        endif
        if((ABS(WORKY_s(i,j)-WORKY(i,j))/WORKY(i,j))>erreps) then
           print*, "advu_diff_LVK=", ABS(WORKY_s(i,j)-WORKY(i,j))/WORKY(i,j),WORKY_s(i,j),WORKY(i,j),i,j,k
           advuy_flg = 0
        endif

        if((ABS(WUK_s(i,j)-WUK(i,j))/WUK(i,j))>erreps) then
           print*, "advu_diff_WUK=", ABS(WUK_s(i,j)-WUK(i,j))/WUK(i,j),WUK_s(i,j),WUK(i,j),i,j,k
           advuk_flg = 0
        endif

        if((ABS(WORKX1_s(i,j)-WORKX1(i,j))/WORKX1(i,j))>erreps) then
           print*, "hdiffu_diff_HDUK=", ABS(WORKX1_s(i,j)-WORKX1(i,j))/WORKX1(i,j),WORKX1_s(i,j),WORKX1(i,j),i,j,k
           hdiffux_flg = 0
        endif
        if((ABS(WORKY1_s(i,j)-WORKY1(i,j))/WORKY1(i,j))>erreps) then
           print*, "hdiffu_diff_HDVK=", ABS(WORKY1_s(i,j)-WORKY1(i,j))/WORKY1(i,j),WORKY1_s(i,j),WORKY1(i,j),i,j,k
           hdiffuy_flg = 0
        endif
      enddo
      enddo
      print* ,"advu LUK with K= ", k, "flag= ", advux_flg
      print* ,"advu LVK with K= ", k, "flag= ", advuy_flg
      print* ,"advu WUK with K= ", k, "flag= ", advuk_flg
      print* ,"hdiffu HDUK with K= ", k, "flag= ", hdiffux_flg
      print* ,"hdiffu HDVK with K= ", k, "flag= ", hdiffuy_flg
      numsteps = numsteps + 1
    endif
#endif

     FX =  -WORKX   ! advu returns WORKX = +L(U) 
     FY =  -WORKY   ! advu returns WORKY = +L(V)
	 
	 
     FX = FX + WORKX1
     FY = FY + WORKY1

   else

     FX = c0
     FY = c0

     WORKX = c0
     WORKY = c0

   end if
        !call timer_stop(call_advu)



        !call timer_start(post_advu)

!!! added-start by wangli
#ifdef WLATH
   iparam(1)=nx_block
   iparam(2)=ny_block
   iparam(3)=k
   iparam(4)=0
   iparam(5)=0
   iparam(6)=0
   if(ldiag_global) iparam(4)=1
   if(partial_bottom_cells) iparam(5)=1
   if(l1ddyn) iparam(6)=1
   if(impcor .and. leapfrogts) then
        iparam(7)=1
   elseif(.not.impcor .and. leapfrogts) then
        iparam(7)=2
   else
        iparam(7)=3
   endif 
   if (lpressure_avg .and. leapfrogts) then
        iparam(8)=1
   else
        iparam(8)=0
   endif
   rparam(1)=c1
   rparam(2)=p5
   rparam(3)=p25
   rparam(4)=dz(k)
   rparam(5)=dzw(k)
   rparam(6)=bouss(k)
   rparam(7)=grav
   rparam(8)=gamma
   rparam(9)=c0
   rparam(10)= dzw(k-1)
 
   call post_advu_c(iparam,rparam,&
                     diag_ke_adv_2d(:,:,bid),dzu(:,:,1,bid),ucur,vcur,&
                     uold,vold,fcor(:,:,bid),fx,fy,&
                     RHOKNEW,RHOKCUR,RHOKOLD,SUMX(:,:,bid),SUMY(:,:,bid),&
                     RHOKMX(:,:,bid),RHOKMY(:,:,bid),dxur(:,:,bid),&
                     dyur(:,:,bid),workx,kmu(:,:,bid))
!!                     dyur(:,:,bid),diag_ke_press_2d(:,:,bid),kmu(:,:,bid))
#else
!!! added-end by wangli


   if (ldiag_global) then
      if (partial_bottom_cells) then
         DIAG_KE_ADV_2D(:,:,bid) = DIAG_KE_ADV_2D(:,:,bid) -           &
                                   DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                                   VCUR(:,:,k)*WORKY)
      else
         DIAG_KE_ADV_2D(:,:,bid) = DIAG_KE_ADV_2D(:,:,bid) -  &
                                   dz(k)*(UCUR(:,:,k)*WORKX + & 
                                          VCUR(:,:,k)*WORKY)
      endif
   endif

!-----------------------------------------------------------------------

   if (impcor .and. leapfrogts) then          ! implicit, leapfrog

      FX = FX + FCOR(:,:,bid)*(      gamma* VCUR(:,:,k) + &
                               (c1 - gamma)*VOLD(:,:,k))
      FY = FY - FCOR(:,:,bid)*(      gamma* UCUR(:,:,k) + & 
                               (c1 - gamma)*UOLD(:,:,k))

   elseif(.not.impcor .and. leapfrogts) then  ! explicit, leapfrog

      FX = FX + FCOR(:,:,bid)*VCUR(:,:,k)
      FY = FY - FCOR(:,:,bid)*UCUR(:,:,k)

   else                                  ! matsuno or foward euler

      FX = FX + FCOR(:,:,bid)*VOLD(:,:,k)
      FY = FY - FCOR(:,:,bid)*UOLD(:,:,k)

   endif

!-----------------------------------------------------------------------
!
!  hydrostatic pressure gradients
!
!-----------------------------------------------------------------------

   if (.not.l1Ddyn) then

!!     call gradp(k,WORKX, WORKY, RHOKOLD, RHOKCUR, RHOKNEW, this_block)



     if (lpressure_avg .and. leapfrogts) then
        RHOAVG = p25*(RHOKNEW + c2*RHOKCUR + RHOKOLD)*bouss(k)
     else
        RHOAVG = RHOKCUR*bouss(k)
     endif

!!     call grad(k,RHOKX,RHOKY,RHOAVG,this_block)
     rhokx = c0
     rhoky = c0
     
     do j=1,ny_block-1
     do i=1,nx_block-1
       if(k<=kmu(i,j,bid)) then
         rhokx(i,j) = dxur(i,j,bid)*p5*(rhoavg(i+1,j+1)-rhoavg(i,j)-&
                      rhoavg(i,j+1)+rhoavg(i+1,j))
         rhoky(i,j) = dyur(i,j,bid)*p5*(rhoavg(i+1,j+1)-rhoavg(i,j)+&
                      rhoavg(i,j+1)-rhoavg(i+1,j))
       endif
     enddo
     enddo

     if (k == 1) then
        RHOKMX(:,:,bid) = RHOKX
        RHOKMY(:,:,bid) = RHOKY
        SUMX  (:,:,bid) = c0
        SUMY  (:,:,bid) = c0
     endif

!-----------------------------------------------------------------------
!
!  obtain pressure gradient by incrementing sum of density gradients
!  from top level down to level k.
!
!-----------------------------------------------------------------------

     factor = dzw(k-1)*grav*p5
     SUMX(:,:,bid) = SUMX(:,:,bid) + factor*(RHOKX + RHOKMX(:,:,bid))
     SUMY(:,:,bid) = SUMY(:,:,bid) + factor*(RHOKY + RHOKMY(:,:,bid))

     WORKX  = SUMX(:,:,bid)
     WORKY  = SUMY(:,:,bid)

!-----------------------------------------------------------------------
!
!  overwrite level k-1 with level k density gradients for next pass
!
!-----------------------------------------------------------------------

     RHOKMX(:,:,bid) = RHOKX
     RHOKMY(:,:,bid) = RHOKY

!-----------------------------------------------------------------------
     FX = FX - WORKX   ! gradp returns WORKX as +Gradx(p)
     FY = FY - WORKY   ! gradp returns WORKY as +Grady(p)

   else

     WORKX = c0
     WORKY = c0

   end if ! not 1D POP

   if (partial_bottom_cells) then
      WORKX =  -DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                VCUR(:,:,k)*WORKY)
   else
      WORKX =  -dz(k)*(UCUR(:,:,k)*WORKX + & 
                       VCUR(:,:,k)*WORKY)
   endif

!!!added-start by wangli
#endif
!!!added-end by wangli

      !call timer_stop(post_advu)
   call accumulate_tavg_field(WORKX,tavg_UDP,bid,k)

   if (ldiag_global) then
      DIAG_KE_PRESS_2D(:,:,bid) = DIAG_KE_PRESS_2D(:,:,bid) + WORKX
   endif

!-----------------------------------------------------------------------
!
!  horizontal diffusion HDiff(Ub),HDiff(Vb)
!
!-----------------------------------------------------------------------
        !call timer_start(call_hdiffu)

!   if (.not.l1Ddyn) then

!     call hdiffu(k, WORKX, WORKY, UMIXK, VMIXK, this_block)

!     FX = FX + WORKX
!     FY = FY + WORKY

!   else

!     WORKX = c0
 !    WORKY = c0

 !  end if
        !call timer_stop(call_hdiffu)
        !call timer_start(post_hdiffu)

   if (ldiag_global) then
      if (partial_bottom_cells) then
         DIAG_KE_HMIX_2D(:,:,bid) = DIAG_KE_HMIX_2D(:,:,bid) +        &
                                  DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                                  VCUR(:,:,k)*WORKY)
      else
         DIAG_KE_HMIX_2D(:,:,bid) = DIAG_KE_HMIX_2D(:,:,bid) + &
                                    dz(k)*(UCUR(:,:,k)*WORKX + & 
                                           VCUR(:,:,k)*WORKY)
      endif
   endif

        !call timer_stop(post_hdiffu)
!-----------------------------------------------------------------------
!
!  vertical diffusion VDiff(Ub),VDiff(Vb)
!
!-----------------------------------------------------------------------

        !call timer_start(call_vdiffu)
   call vdiffu(k, WORKX, WORKY, UOLD, VOLD, SMF_BLOCK, this_block)

        !call timer_stop(call_vdiffu)
        !call timer_start(post_vdiffu)
   FX = FX + WORKX
   FY = FY + WORKY

   if (ldiag_global) then
      if (partial_bottom_cells) then
         DIAG_KE_VMIX_2D(:,:,bid) = DIAG_KE_VMIX_2D(:,:,bid) + &
                           DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                           VCUR(:,:,k)*WORKY)
      else
         DIAG_KE_VMIX_2D(:,:,bid) = DIAG_KE_VMIX_2D(:,:,bid) + &
                                    dz(k)*(UCUR(:,:,k)*WORKX + & 
                                           VCUR(:,:,k)*WORKY)
      endif
   endif

!-----------------------------------------------------------------------
!
!  zero forces (and hence velocities) at land points
!
!-----------------------------------------------------------------------

   where (k > KMU(:,:,bid))
      FX = c0
      FY = c0
   endwhere
        !call timer_stop(post_vdiffu)

!-----------------------------------------------------------------------
!EOC

 end subroutine clinic

!***********************************************************************
!BOP
! !IROUTINE: tracer_update
! !INTERFACE:

 subroutine tracer_update(k, WTK, TNEW, TOLD, TMIX, TCUR,            &
                             UCUR, VCUR, UMIX, VMIX, RHOCUR,         &
                             STF_IN, TFW_IN, QSW, DH_IN, POLD, PCUR, &
                             this_block)

! !DESCRIPTION:
!  Computes explicit forcing for tracer equations:
!  \begin{equation}
!     (T^{n+1}-T^{n-1})/(2 \Delta t) = -L(T) + D_H(T^{n-1}) + 
!                                              D_V(T^{n-1}) + S
!  \end{equation}
!  where $L$ is the advection operator, $D_{H,V}$ are the diffusion
!  operators in the horizontal and vertical and $S$ are source terms.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k     ! depth level index

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TCUR,                 &! tracers at current time level
      TOLD,                 &! tracers at old     time level
      TMIX                   ! tracers at mix     time level

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UCUR, VCUR,           &! U,V for block at current time
      UMIX, VMIX,           &! U,V at mix time level
      RHOCUR                 ! density at current time level

   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      STF_IN,               &! surface tracer fluxes
      TFW_IN                 ! tracer concentration in fresh water

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      DH_IN,                &! sfc height change at tracer points
      POLD,                 &! sfc pressure at old     time
      PCUR,                 &! sfc pressure at current time
      QSW                    ! short-wave heating

   type (block), intent(in) :: &
      this_block             ! block info for the current block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(inout) :: &
      WTK          ! on  input, vertical velocity at top    of T box
                   ! on output, vertical velocity at bottom of T box

   real (r8), dimension(nx_block,ny_block,km,nt), intent(inout) :: &
      TNEW                   ! tracers at new time level

!EOP
!BOC
!-----------------------------------------------------------------------
!         
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! dummy tracer index
      bid                 ! local_block id

   real (r8), dimension(nx_block,ny_block,nt) :: &
      FT,                &! sum of terms in dT/dt for the nth tracer
      WORKN,WORKN1               ! work array used for various dT/dt terms 

   real (r8), dimension(nx_block,ny_block) :: &
      WORKSW



#if defined(adhdt_chk)
   real (r8), dimension(nx_block,ny_block) :: &
      WTK_ini,WTK_s        ! local work space to hold master forcing terms
   real (r8), dimension(nx_block,ny_block,nt) :: &
      WORKN_s,           &   ! local work space to hold master forcing terms
      WORKN1_s               ! local work space to hold master hdiiff.
   real (r8) :: erreps, erreps0, erreps2         

   integer (int_kind) :: i,j
   integer (int_kind) :: advt_flg,advtk_flg
   integer (int_kind) :: hdifft_flg
   integer (int_kind), save :: numsteps=0
#endif
 



!-----------------------------------------------------------------------
!
!  initialize some arrays
!
!-----------------------------------------------------------------------


   FT    = c0

   bid = this_block%local_id


! especially, for advt:
   if (k == 1) WTK = DH_IN

   if (ldiag_global) then

      if (k == 1) then
         where (CALCT(:,:,bid)) &
            DIAG_PE_2D(:,:,bid) = dzw(0)*grav*RHOCUR(:,:,1)*WTK

         if (sfc_layer_type /= sfc_layer_varthick) then
            do n=1,nt
               where (CALCT(:,:,bid))
                  DIAG_TRACER_ADV_2D (:,:,n,bid) = WTK*TCUR(:,:,1,n)
                  DIAG_TRACER_SFC_FLX(:,:,n,bid) = WTK*TCUR(:,:,1,n)
               elsewhere
                  DIAG_TRACER_ADV_2D (:,:,n,bid) = c0
                  DIAG_TRACER_SFC_FLX(:,:,n,bid) = c0
               end where
            end do
         endif

      else

         !*** For energetic consistency, we use dzw even for 
         !*** partial bottom cell case

         where (k <= KMT(:,:,bid))
            DIAG_PE_2D(:,:,bid) = DIAG_PE_2D(:,:,bid) + &
                                  dzw(k-1)*WTK*grav*p5* &
                                  (RHOCUR(:,:,k-1) + RHOCUR(:,:,k))
         endwhere

      endif

   endif




! new combined advt and hdifft:
   if (.not.l1Ddyn) then
#if defined(adhdt_chk)
     WTK_ini = WTK
#endif

#if !defined(NOATH_adhdt_GT)
     call advthdifft(k,WORKN,WORKN1,WTK,TMIX,UMIX,VMIX,TCUR,UCUR,VCUR,this_block)
#endif




#if defined(adhdt_chk)
     WORKN_s  = WORKN
     WTK_s    = WTK
     WORKN1_s = WORKN1
     WTK = WTK_ini
#endif




#if defined(adhdt_chk) || defined(NOATH_adhdt_GT) 
! original calls:
     call hdifft(k, WORKN, TMIX, UMIX, VMIX, this_block)
     call advt(k,WORKN1,WTK,TMIX,TCUR,UCUR,VCUR,this_block)
#endif


#if defined(adhdt_chk)
    if(my_task==0 .and. numsteps<=768) then
      advt_flg = 1
      advtk_flg = 1
      hdifft_flg = 1
      erreps  = 1.0E-10
      erreps2 = 1.0E-20
      erreps0 = 1.0E-20
      do n=1,nt
      do j=1,ny_block
      do i=1,nx_block
!        if((ABS(WORKN_s(i,j,n)-WORKN(i,j,n))/WORKN(i,j,n))>erreps .and. ABS(WORKN(i,j,n))>erreps2) then
        if((ABS(WORKN_s(i,j,n)-WORKN(i,j,n))/WORKN(i,j,n)) * ABS(WORKN(i,j,n)) >erreps0) then
           print*, "hdifft_diff_HDTK=", ABS(WORKN_s(i,j,n)-WORKN(i,j,n))/WORKN(i,j,n),WORKN_s(i,j,n),WORKN(i,j,n),i,j,n,k
           hdifft_flg = 0
        endif

!        if((ABS(WORKN1_s(i,j,n)-WORKN1(i,j,n))/WORKN1(i,j,n))>erreps .and. ABS(WORKN1(i,j,n))>erreps2) then
        if((ABS(WORKN1_s(i,j,n)-WORKN1(i,j,n))/WORKN1(i,j,n)) * ABS(WORKN1(i,j,n)) >erreps0) then
           print*, "advt_diff_LTK=", ABS(WORKN1_s(i,j,n)-WORKN1(i,j,n))/WORKN1(i,j,n),WORKN1_s(i,j,n),WORKN1(i,j,n),i,j,n,k
           advt_flg = 0
        endif

        if(n==1) then
        if((ABS(WTK_s(i,j)-WTK(i,j))/WTK(i,j))>erreps) then
           print*, "advt_diff_WTK=", ABS(WTK_s(i,j)-WTK(i,j))/WTK(i,j),WTK_s(i,j),WTK(i,j),i,j,n,k
           advtk_flg = 0
        endif
        endif
      enddo
      enddo
      enddo
      print* ,"advt LTK with K= ", k, "flag= ", advt_flg
      print* ,"advt WTK with K= ", k, "flag= ", advtk_flg
      print* ,"hdifft HDTK with K= ", k, "flag= ", hdifft_flg
      numsteps = numsteps + 1
    endif
#endif


     FT = FT + WORKN
     FT = FT - WORKN1   ! advt returns WORKN = +L(T) 

   else

     WORKN = c0
     WORKN1 = c0

   end if






!-----------------------------------------------------------------------
!
!  horizontal diffusion HDiff(T)
!
!-----------------------------------------------------------------------
! original:
!   if (.not.l1Ddyn) then
!
!     call hdifft(k, WORKN, TMIX, UMIX, VMIX, this_block)
!     FT = FT + WORKN
!
!   else
!
!     WORKN = c0
!
!   end if
!

   if (ldiag_global) then

      if (partial_bottom_cells) then
         do n=1,nt
            where (k <= KMT(:,:,bid))            &
               DIAG_TRACER_HDIFF_2D(:,:,n,bid) = &
               DIAG_TRACER_HDIFF_2D(:,:,n,bid) + &
               WORKN(:,:,n)*DZT(:,:,k,bid)
         end do
      else
         do n=1,nt
            where (k <= KMT(:,:,bid))            &
               DIAG_TRACER_HDIFF_2D(:,:,n,bid) = & 
               DIAG_TRACER_HDIFF_2D(:,:,n,bid) + &
               WORKN(:,:,n)*dz(k)
         end do
      endif
   endif

!-----------------------------------------------------------------------
!
!  advection L(T)
!  set vertical velocity at surface
!  w = DH = dh/dt or dh/dt - Fw depending on surface type
!
!-----------------------------------------------------------------------

! original:
!   if (k == 1) WTK = DH_IN
!
!   if (ldiag_global) then
!
!      if (k == 1) then
!         where (CALCT(:,:,bid)) &
!            DIAG_PE_2D(:,:,bid) = dzw(0)*grav*RHOCUR(:,:,1)*WTK
!
!         if (sfc_layer_type /= sfc_layer_varthick) then
!            do n=1,nt
!               where (CALCT(:,:,bid))
!                  DIAG_TRACER_ADV_2D (:,:,n,bid) = WTK*TCUR(:,:,1,n)
!                  DIAG_TRACER_SFC_FLX(:,:,n,bid) = WTK*TCUR(:,:,1,n)
!               elsewhere
!                  DIAG_TRACER_ADV_2D (:,:,n,bid) = c0
!                  DIAG_TRACER_SFC_FLX(:,:,n,bid) = c0
!               end where
!            end do
!         endif
!
!      else
!
!         !*** For energetic consistency, we use dzw even for 
!         !*** partial bottom cell case
!
!         where (k <= KMT(:,:,bid))
!            DIAG_PE_2D(:,:,bid) = DIAG_PE_2D(:,:,bid) + &
!                                  dzw(k-1)*WTK*grav*p5* &
!                                  (RHOCUR(:,:,k-1) + RHOCUR(:,:,k))
!         endwhere
!
!      endif
!
!   endif
!
!   if (.not.l1Ddyn) then
!
!     call advt(k,WORKN1,WTK,TMIX,TCUR,UCUR,VCUR,this_block)
!     FT = FT - WORKN1   ! advt returns WORKN = +L(T) 
!
!   else
!
!     WORKN1 = c0
!
!   end if
!

   if (ldiag_global) then
     if (partial_bottom_cells) then
       do n=1,nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_ADV_2D(:,:,n,bid) = & 
                                   DIAG_TRACER_ADV_2D(:,:,n,bid) - &
                                   WORKN1(:,:,n)*DZT(:,:,k,bid)
       end do
     else
       do n=1,nt

         where (k <= KMT(:,:,bid)) DIAG_TRACER_ADV_2D(:,:,n,bid) = &
                                   DIAG_TRACER_ADV_2D(:,:,n,bid) - &
                                   WORKN1(:,:,n)*dz(k)
       end do
     endif 
   endif

!-----------------------------------------------------------------------
!
!  vertical diffusion VDiff(T)
!
!-----------------------------------------------------------------------

   call vdifft(k, WORKN, TOLD, STF_IN, this_block)

   FT = FT + WORKN

   if (ldiag_global) then
     if (partial_bottom_cells) then
       do n=1,nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_VDIFF_2D(:,:,n,bid) = &
                                   DIAG_TRACER_VDIFF_2D(:,:,n,bid) + &
                                   WORKN(:,:,n)*DZT(:,:,k,bid)
       end do
     else
       do n=1,nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_VDIFF_2D(:,:,n,bid) = &
                                   DIAG_TRACER_VDIFF_2D(:,:,n,bid) + &
                                   WORKN(:,:,n)*dz(k)
       end do
     endif
   endif

!-----------------------------------------------------------------------
!
!  add tracer change in surface layer due to freshwater flux
!  if using variable thickness surface layer
!
!-----------------------------------------------------------------------

   if (k == 1 .and. sfc_layer_type == sfc_layer_varthick) then
      do n = 1,nt
         FT(:,:,n) = FT(:,:,n) + dzr(1)*TFW_IN(:,:,n)
      enddo

      if (ldiag_global) then
         do n = 1,nt
            DIAG_TRACER_SFC_FLX(:,:,n,bid) = TFW_IN(:,:,n)
         enddo
      endif
   endif

!-----------------------------------------------------------------------
!
!  add source terms
!
!-----------------------------------------------------------------------

   WORKN = c0

   call set_pt_interior(k,this_block,WORKN(:,:,1))
   call set_s_interior (k,this_block,WORKN(:,:,2))


   if (nt > 2) call set_interior_passive_tracers(k, this_block, WORKN)
   

!-----------------------------------------------------------------------
!
!  add source terms from KPP and from shortwave solar absorption
!    if necessary.
!  NOTE:  this is here instead of in set_{pt,s}_interior in case
!    KPP and/or shortwave solar absorption are turned on but 
!    bulk restoring is not.
!
!-----------------------------------------------------------------------

   !*** does nothing if kpp not chosen - otherwise adds kpp sources
   call add_kpp_sources(WORKN, k, this_block)

   !*** if sw flux available, add penetrative shortwave
   call add_sw_absorb(WORKN, SHF_QSW(:,:,bid), k, this_block)


   FT = FT + WORKN


   if (ldiag_global) then
     if (partial_bottom_cells) then
       do n=1,nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_SOURCE_2D(:,:,n,bid) = &
                                   DIAG_TRACER_SOURCE_2D(:,:,n,bid) + &
                                   WORKN(:,:,n)*DZT(:,:,k,bid)
       end do
     else
       do n=1,nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_SOURCE_2D(:,:,n,bid) = &
                                   DIAG_TRACER_SOURCE_2D(:,:,n,bid) + &
                                   WORKN(:,:,n)*dz(k)
       end do
     endif
   endif


!-----------------------------------------------------------------------
!
!  save the explicit part of the RHS in TRACER(newtime)
!  if there is implicit vertical mixing
!
!  with pressure averaging and variable thickness surface layer, 
!  the RHS contains the surface height contribution for the 
!  predictor step (for T,S at k=1 only)
!
!-----------------------------------------------------------------------

   if (implicit_vertical_mix) then

      if (sfc_layer_type == sfc_layer_varthick .and. k == 1 .and. &
          lpressure_avg .and. leapfrogts) then

         do n = 1,2
            where (KMT(:,:,bid) > 0)  ! RHS for predictor
               TNEW(:,:,1,n) = c2dtt(1)*FT(:,:,n) - c2*TCUR(:,:,1,n)* &
                               (PCUR - POLD)/(grav*dz(1))
            endwhere
         enddo

         do n = 3,nt

            TNEW(:,:,k,n) = merge(c2dtt(k)*FT(:,:,n), &
                                  c0, k <= KMT(:,:,bid))
         enddo

      else

         do n = 1,nt
            TNEW(:,:,k,n) = merge(c2dtt(k)*FT(:,:,n), &
                                  c0, k <= KMT(:,:,bid))
         enddo

      endif

!-----------------------------------------------------------------------
!
!  for variable thickness surface layer, update all but surface
!    layers. 
!  at the surface:
!    if explicit vertical mixing and pressure averaging:
!      predict surface T,S and store RHS for all other tracers
!    otherwise
!      store RHS for all tracers for later update with new Psurf
!
!  if not a variable thickness surface layer, update tracers here
!     
!-----------------------------------------------------------------------

   else ! no implicit_vertical_mix

      if (sfc_layer_type == sfc_layer_varthick .and. k == 1) then
         if (lpressure_avg .and. leapfrogts) then

            !*** predict surface T and S with pressure avg

            do n = 1,2
               where (KMT(:,:,bid) > 0)
                  TNEW(:,:,1,n) = TOLD(:,:,1,n)                        &
                           + (c1/(c1 + PCUR/(grav*dz(1))))*            &
                             (c2dtt(1)*FT(:,:,n) - c2*TCUR(:,:,1,n)*   &
                             (PCUR-POLD)/(grav*dz(1)))
               elsewhere
                  TNEW(:,:,1,n) = c0  ! zero tracers on land pts
               endwhere
            end do

            !*** store RHS for other tracer surface layers

            do n = 3,nt
               TNEW(:,:,k,n) = merge(c2dtt(k)*FT(:,:,n), &
                                     c0, k <= KMT(:,:,bid))
            enddo

         else

            !*** store RHS for all tracer surface layers

            do n = 1,nt
               TNEW(:,:,k,n) = merge(c2dtt(k)*FT(:,:,n), &
                                     c0, k <= KMT(:,:,bid))
            enddo

         endif

      else   !*** update all tracers to new time

         do n = 1,nt
            where (k <= KMT(:,:,bid))
               TNEW(:,:,k,n) = TOLD(:,:,k,n) + c2dtt(k)*FT(:,:,n)
            elsewhere
               TNEW(:,:,k,n) = c0  ! zero tracers at land points
            endwhere
         enddo

      endif
   endif ! implicit_vertical_mix


!-----------------------------------------------------------------------
!EOC

 end subroutine tracer_update

!***********************************************************************

 end module baroclinic

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
