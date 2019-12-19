!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ecosys_driver

!BOP
! !MODULE: ecosys_driver

! !DESCRIPTION:
!  This module provides support for the ecosystem module and dependend tracer modules
!  The base model calls subroutines in passive_tracers, which then call
!  this module if ecosys_on is true. Ecosys_driver then calls subroutines
!  in individual ecosystem modules (so far ecosys_mod and ecosys_ciso_mod)
!
!  Written by: Alexandra Jahn, NCAR, Nov/Dec 2012


! !REVISION HISTORY:
!  SVN:$Id: passive_tracers.F90 28439 2011-05-18 21:40:58Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod

   use kinds_mod, only: r8, int_kind, log_kind, char_len
   use blocks, only: block, nx_block, ny_block
   use domain_size, only: max_blocks_clinic, km, nt
   use communicate, only: my_task, master_task
   use prognostic, only: TRACER, tracer_d, oldtime, curtime, newtime
   use forcing_shf, only: SHF_QSW_RAW, SHF_QSW
   use io_types, only: stdout, nml_in, nml_filename, io_field_desc, &
       datafile
   use exit_mod, only: sigAbort, exit_pop
   use io_tools, only: document
   use prognostic, only: tracer_field
   use constants, only: c0, c1, p5, delim_fmt, char_blank, ndelim_fmt
   use passive_tracer_tools, only: set_tracer_indices
   use ecosys_share, only: lmarginal_seas, ecosys_tadvect_ctype
   use broadcast, only: broadcast_scalar


   use ecosys_constants, only : ecosys_tracer_cnt

   use ecosys_mod, only:           &
       ecosys_init,                &
       ecosys_tracer_ref_val,      &
       ecosys_set_sflux,           &
       ecosys_tavg_forcing,        &
       ecosys_set_interior,        &
       ecosys_write_restart

   use ecosys_ciso_mod, only:       &
       ecosys_ciso_tracer_cnt,      &
       ecosys_ciso_init,            &
       ecosys_ciso_tracer_ref_val,  &
       ecosys_ciso_set_sflux,       &
       ecosys_ciso_tavg_forcing,    &
       ecosys_ciso_set_interior,    &
       ecosys_ciso_write_restart

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public ::                          &
      ecosys_driver_tracer_cnt_init,  &
      ecosys_driver_init,             &
      ecosys_driver_set_interior,     &
      ecosys_driver_set_sflux,        &
      ecosys_driver_tracer_cnt,       &
      ecosys_driver_tracer_ref_val,   &
      ecosys_driver_tavg_forcing,     &
      ecosys_driver_write_restart

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by forcing_passive_tracer
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ecosys_driver_tracer_cnt


!-----------------------------------------------------------------------
!     index bounds of passive tracer module variables in TRACER
!-----------------------------------------------------------------------

   integer (kind=int_kind) :: &
      ecosys_driver_ind_begin,  ecosys_driver_ind_end,  &
      ecosys_ind_begin,         ecosys_ind_end,         &
      ecosys_ciso_ind_begin,    ecosys_ciso_ind_end


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: ecosys_driver_tracer_cnt_init
! !INTERFACE:

 subroutine ecosys_driver_tracer_cnt_init(ciso_on)

! !DESCRIPTION:
!  Zero-level initialization of ecosys_driver,
!  which involves setting the ecosys_driver_tracer_cnt
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (kind=log_kind), intent(in)  ::  &
     ciso_on                 ! ecosys_ciso on

!EOP
!BOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine ecosys_driver_tracer_cnt, depending on whether only ecosys
!  or also other modules are on
!-----------------------------------------------------------------------

    ecosys_driver_tracer_cnt = ecosys_tracer_cnt

   if (ciso_on) then
       ecosys_driver_tracer_cnt = ecosys_driver_tracer_cnt + ecosys_ciso_tracer_cnt
   end if


!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_driver_tracer_cnt_init



!***********************************************************************
!BOP
! !IROUTINE: ecosys_driver_init
! !INTERFACE:

 subroutine ecosys_driver_init(ciso_on,init_ts_file_fmt,     &
                        read_restart_filename,                         &
                        tracer_d_module, TRACER_MODULE, tadvect_ctype, &
                        errorCode)


! !DESCRIPTION:
!  Initialize ecosys_driver passive tracers. This involves:
!  1) setting ecosys and ecosys_ciso module index bounds
!  2) calling ecosys and ecosys_ciso module init subroutine
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file


   logical (kind=log_kind), intent(in)  ::  &
     ciso_on                 ! ecosys_ciso on

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(:), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(:,:,:,:,:,:), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   character (char_len), dimension(:), intent(out) :: &
      tadvect_ctype     ! advection method for ecosys tracers

   integer (POP_i4), intent(out) :: &
      errorCode

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'ecosys_driver:ecosys_driver_init'

   integer (int_kind) :: cumulative_nt, n, &
      nml_error,        &! error flag for nml read
      iostat             ! io status flag

   character (char_len) :: sname, lname, units, coordinates
   character (4) :: grid_loc


!-----------------------------------------------------------------------
!  read in ecosys_driver namelist, to set namelist parameters that
!  should be the same for all ecosystem-related modules
!-----------------------------------------------------------------------

  namelist /ecosys_driver_nml/ &
      lmarginal_seas, ecosys_tadvect_ctype

   errorCode = POP_Success

   lmarginal_seas        = .true.
   ecosys_tadvect_ctype  = 'base_model'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
        nml_error = -1
      else
        nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
        read(nml_in, nml=ecosys_driver_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading ecosys_driver namelist')
   endif

    if (my_task == master_task) then
      write(stdout,*)
      write(stdout,ndelim_fmt)
      write(stdout,*)
      write(stdout,*) ' ecosys_driver:'
      write(stdout,*)
      write(stdout,*) ' ecosys_driver_nml namelist settings:'
      write(stdout,*)
      write(stdout,ecosys_driver_nml)
      write(stdout,*)
      write(stdout,delim_fmt)
   endif

   call broadcast_scalar(ecosys_tadvect_ctype, master_task)
   call broadcast_scalar(lmarginal_seas, master_task)

!-----------------------------------------------------------------------
!  set up indices for ecosys modules that are on
!  These indices are relative to ecosys_driver_ind_begin/end
!  This means they start at 1 and end at ecosys_driver_tracer_cnt
!  ecosys_driver then passes the ecosys module tracers back to passive
!  tracers within TRACER(:,:,:,ecosys_driver_ind_beg,ecosys_driver_ind_end)
!-----------------------------------------------------------------------

   cumulative_nt = 0

   call set_tracer_indices('ECOSYS', ecosys_tracer_cnt, cumulative_nt,  &
                                   ecosys_ind_begin, ecosys_ind_end)


   if (ciso_on) then
      call set_tracer_indices('CISO', ecosys_ciso_tracer_cnt, cumulative_nt,  &
                              ecosys_ciso_ind_begin, ecosys_ciso_ind_end)
   end if


   if (cumulative_nt /= ecosys_driver_tracer_cnt) then
      call document(subname, 'ecosys_driver_tracer_cnt', ecosys_driver_tracer_cnt)
      call document(subname, 'cumulative_nt', cumulative_nt)
      call exit_POP(sigAbort, &
         'ERROR in ecosys_driver_init: ecosys_driver_tracer_cnt does not match cumulative nt')
   end if


!-----------------------------------------------------------------------
!  ECOSYS block
!-----------------------------------------------------------------------

   tadvect_ctype(ecosys_ind_begin:ecosys_ind_end) = ecosys_tadvect_ctype
   call ecosys_init(init_ts_file_fmt, read_restart_filename,                  &
                    tracer_d_module(ecosys_ind_begin:ecosys_ind_end),         &
                    TRACER_MODULE(:,:,:,ecosys_ind_begin:ecosys_ind_end,:,:), &
                    lmarginal_seas,                                           &
                    errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_init')
      return
   endif


!-----------------------------------------------------------------------
!  ECOSYS CISO block
!-----------------------------------------------------------------------
   if (ciso_on) then
      tadvect_ctype(ecosys_ciso_ind_begin:ecosys_ciso_ind_end) = ecosys_tadvect_ctype
      call ecosys_ciso_init(init_ts_file_fmt, read_restart_filename,                       &
                       tracer_d_module(ecosys_ciso_ind_begin:ecosys_ciso_ind_end),         &
                       TRACER_MODULE(:,:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:,:), &
                       lmarginal_seas,                               &
                       errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_ecosys_driver: error in ecosys_ciso_init')
         return
      endif

   end if


!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_driver_init

!***********************************************************************
!BOP
! !IROUTINE: ecosys_driver_set_interior
! !INTERFACE:

 subroutine ecosys_driver_set_interior(k, ciso_on, TEMP_OLD, &
    TEMP_CUR, SALT_OLD, SALT_CUR, &
    TRACER_MODULE_OLD, TRACER_MODULE_CUR, DTRACER_MODULE, this_block)

! !DESCRIPTION:
!  call subroutines for each tracer module that compute source-sink terms
!  accumulate commnon tavg fields related to source-sink terms
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                   ! vertical level index

   logical (kind=log_kind), intent(in)  ::  &
     ciso_on                 ! ecosys_ciso on

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      TEMP_OLD,          &! old potential temperature (C)
      TEMP_CUR,          &! current potential temperature (C)
      SALT_OLD,          &! old salinity (msu)
      SALT_CUR            ! current salinity (msu)

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
!  ECOSYS block
!-----------------------------------------------------------------------

   call ecosys_set_interior(k,                                 &
      TEMP_OLD, TEMP_CUR,                                      &
      SALT_OLD, SALT_CUR,                                      &
      TRACER_MODULE_OLD(:,:,:,ecosys_ind_begin:ecosys_ind_end),&
      TRACER_MODULE_CUR(:,:,:,ecosys_ind_begin:ecosys_ind_end),&
      DTRACER_MODULE(:,:,ecosys_ind_begin:ecosys_ind_end),     &
      ciso_on,                                                 &
      this_block)

!-----------------------------------------------------------------------
!  ECOSYS_CISO block
!-----------------------------------------------------------------------

   if (ciso_on) then
      call ecosys_ciso_set_interior(k,                                      &
         TEMP_OLD, TEMP_CUR,                                                &
         TRACER_MODULE_OLD(:,:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end),&
         TRACER_MODULE_CUR(:,:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end),&
         DTRACER_MODULE(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end),     &
         this_block)
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_driver_set_interior

!***********************************************************************
!BOP
! !IROUTINE: ecosys_driver_set_sflux
! !INTERFACE:

 subroutine ecosys_driver_set_sflux(ciso_on,SHF_QSW_RAW, SHF_QSW, &
                             U10_SQR,IFRAC,PRESS,SST,SSS, &
                             SURFACE_VALS_OLD,SURFACE_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  call subroutines for each tracer module that compute surface fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  logical (kind=log_kind), intent(in)  ::  &
     ciso_on                 ! ecosys_ciso on

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      SHF_QSW_RAW,  &! penetrative solar heat flux, from coupler (degC*cm/s)
      SHF_QSW,      &! SHF_QSW used by physics, may have diurnal cylce imposed (degC*cm/s)
      U10_SQR,      &! 10m wind speed squared (cm/s)**2
      IFRAC,        &! sea ice fraction (non-dimensional)
      PRESS,        &! sea level atmospheric pressure (dyne/cm**2)
      SST,          &! sea surface temperature (C)
      SSS            ! sea surface salinity (psu)

   real (r8), dimension(:,:,:,:), &
      intent(in) :: SURFACE_VALS_OLD, SURFACE_VALS_CUR ! module tracers

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:,:), &
      intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  ECOSYS block
!-----------------------------------------------------------------------

   call ecosys_set_sflux(                                       &
      SHF_QSW_RAW, SHF_QSW,                                     &
      U10_SQR, IFRAC, PRESS,                                    &
      SST, SSS,                                                 &
      SURFACE_VALS_OLD(:,:,ecosys_ind_begin:ecosys_ind_end,:),  &
      SURFACE_VALS_CUR(:,:,ecosys_ind_begin:ecosys_ind_end,:),  &
      STF_MODULE(:,:,ecosys_ind_begin:ecosys_ind_end,:),        &
      ciso_on)

!-----------------------------------------------------------------------
!  ECOSYSC_CISO block
!-----------------------------------------------------------------------

   if (ciso_on) then
      call ecosys_ciso_set_sflux(                                         &
         SST, &
         SURFACE_VALS_OLD(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:), &
         SURFACE_VALS_CUR(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:), &
         STF_MODULE(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:))
   end if


!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_driver_set_sflux

!***********************************************************************
!BOP
! !IROUTINE: ecosys_driver_write_restart
! !INTERFACE:

 subroutine ecosys_driver_write_restart(ciso_on,restart_file, action)

! !DESCRIPTION:
!  call restart routines for each tracer module that
!  write fields besides the tracers themselves
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   logical (kind=log_kind), intent(in)  ::  &
     ciso_on                 ! ecosys_ciso on

   character(*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC

!-----------------------------------------------------------------------
!  ECOSYS block
!-----------------------------------------------------------------------

   call ecosys_write_restart(restart_file, action)

!-----------------------------------------------------------------------
!  ECOSYS_CISO block
!-----------------------------------------------------------------------
   if (ciso_on) then
      call ecosys_ciso_write_restart(restart_file, action)
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_driver_write_restart


!***********************************************************************
!BOP
! !IROUTINE: ecosys_driver_tavg_forcing
! !INTERFACE:

 subroutine ecosys_driver_tavg_forcing(ciso_on,STF_MODULE)

! !DESCRIPTION:
!  accumulate common tavg fields for tracer surface fluxes
!  call accumation subroutines for tracer modules that have additional
!     tavg fields related to surface fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
 logical (kind=log_kind), intent(in)  ::  &
     ciso_on                 ! ecosys_ciso on

  real (r8), dimension(:,:,:,:), &
     intent(in) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  call routines from modules that have additional sflux tavg fields
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ECOSYS block
!-----------------------------------------------------------------------

   call ecosys_tavg_forcing(STF_MODULE(:,:,ecosys_ind_begin:ecosys_ind_end,:))

!-----------------------------------------------------------------------
!  ECOSYS_CISO block
!-----------------------------------------------------------------------

   if (ciso_on) then
      call ecosys_ciso_tavg_forcing( &
         STF_MODULE(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:))
   end if


!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_driver_tavg_forcing



!***********************************************************************
!BOP
 !IROUTINE: ecosys_driver_tracer_ref_val
! !INTERFACE:
!
 function ecosys_driver_tracer_ref_val(ciso_on,ind)
!
! !DESCRIPTION:
!  return reference value for tracer with global tracer index ind
!  this is used in virtual flux computations
!
! !REVISION HISTORY:
!  same as module

 !INPUT PARAMETERS:
   logical (kind=log_kind), intent(in)  ::  &
     ciso_on                 ! ecosys_ciso on

   integer(int_kind), intent(in) :: ind

 !OUTPUT PARAMETERS:

   real(r8) :: ecosys_driver_tracer_ref_val

!EOP
!BOC

!-----------------------------------------------------------------------
!  default value for reference value is 0
!-----------------------------------------------------------------------

   ecosys_driver_tracer_ref_val = c0

!-----------------------------------------------------------------------
!  ECOSYS block
!-----------------------------------------------------------------------

   if (ind >= ecosys_ind_begin .and. ind <= ecosys_ind_end) then
      ecosys_driver_tracer_ref_val = ecosys_tracer_ref_val(ind-ecosys_ind_begin+1)
   endif

!-----------------------------------------------------------------------
!  ECOSYS_CISO block
!-----------------------------------------------------------------------

   if (ciso_on) then
      if (ind >= ecosys_ciso_ind_begin .and. ind <= ecosys_ciso_ind_end) then
         ecosys_driver_tracer_ref_val = ecosys_ciso_tracer_ref_val(ind-ecosys_ciso_ind_begin+1)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end function ecosys_driver_tracer_ref_val

!***********************************************************************

 end module ecosys_driver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
