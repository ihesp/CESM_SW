module ecosys_restore_mod
  !
  ! Module to generalize restoring any non-autotroph tracer
  !

  use kinds_mod, only : r8, log_kind, int_kind
  use domain_size, only : km
  use ecosys_constants, only : ecosys_tracer_cnt
  use passive_tracer_tools, only : tracer_read, tracer_read_init
  use ecosys_restore_timescale_file, only : ecosys_restore_timescale_file_type
  use ecosys_restore_timescale_interp, only : ecosys_restore_timescale_interp_type

  implicit none

  private

  type, private :: restore_info
     logical(log_kind) :: restore ! flag indicating if this tracer should be restored
     type(tracer_read) :: restore_file_info ! info about file containing restoring field
     real (r8), dimension(:, :, :, :), allocatable :: data ! restoring field
     integer(int_kind) :: tavg_restore_index ! index for tavg output
  end type restore_info

  type, public :: ecosys_restore_type
     logical(log_kind), public :: restore_any_tracer ! true if we are restoring any tracer
     type(restore_info), allocatable, private :: tracers(:)

     ! true if geographically varying nutrient restoring is read from
     ! a file (formally lnutr_variable_restore)
     logical(log_kind), private :: spatial_variability_from_file
     type(ecosys_restore_timescale_file_type), private :: timescale_file
     type(ecosys_restore_timescale_interp_type), private :: timescale_interp

   contains
     procedure, public :: init
     procedure, public :: read_restoring_fields
     procedure, public :: restore_variable
     procedure, public :: define_tavg_fields
     procedure, public :: accumulate_tavg
     procedure, public :: initialize_restoring_timescale
     procedure, private :: read_namelist
     procedure, private :: initialize_restore_read_vars

  end type ecosys_restore_type

contains

!*****************************************************************************

subroutine init(this, nml_filename, nml_in, ind_name_table)
  ! initialize ecosys_restore instance to default values, then read
  ! namelist and setup tracers that need to be restored

  use kinds_mod, only : char_len, int_kind, i4, log_kind
  use constants, only : c0, c2, c1000
  use passive_tracer_tools, only : tracer_read, tracer_read_init, &
       ind_name_pair, name_to_ind

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type) :: this
  character(len=*), intent(in) :: nml_filename
  integer(i4), intent(in) :: nml_in
  type(ind_name_pair), dimension(:), intent(in) :: ind_name_table

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer (int_kind) :: t
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_short_names
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_filenames
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_file_varnames

  !-----------------------------------------------------------------------

  ! initialize class data to default values
  this%restore_any_tracer = .false.
  allocate(this%tracers(ecosys_tracer_cnt))
  ! NOTE(bja, 2014-10) don't allocate tracers(t)%data here
  ! because we don't know if it is needed yet!
  do t = 1, ecosys_tracer_cnt
     call tracer_read_init(this%tracers(t)%restore_file_info)
     this%tracers(t)%restore = .false.
     this%tracers(t)%tavg_restore_index = -1
  end do

  call this%read_namelist(nml_filename, nml_in, &
       restore_short_names, restore_filenames, restore_file_varnames)

  call this%initialize_restore_read_vars(restore_short_names, restore_filenames, &
       restore_file_varnames, ind_name_table)

end subroutine Init


!*****************************************************************************
subroutine read_namelist(this, nml_filename, nml_in, &
     restore_short_names, restore_filenames, restore_file_varnames)

  ! Read the ecosys_restore namelist and broadcast to all
  ! processes. Store results in the ecosys_restore_vars

  use kinds_mod, only : char_len, char_len, int_kind, i4
  use communicate, only : master_task, my_task
  use broadcast, only : broadcast_scalar
  use io_types, only : stdout

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type) :: this
  character(len=*), intent(in) :: nml_filename
  integer(i4), intent(in) :: nml_in
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(out) :: restore_short_names
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(out) :: restore_filenames
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(out) :: restore_file_varnames

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: nml_error, t
  logical(log_kind) :: spatial_variability_from_file

  !-----------------------------------------------------------------------

  namelist /ecosys_restore_nml/ &
       restore_short_names, &
       restore_filenames, &
       restore_file_varnames, &
       spatial_variability_from_file

  ! initialize namelist variables to default values
  restore_short_names = ''
  restore_filenames = ''
  restore_file_varnames = ''
  spatial_variability_from_file = .false.

  if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     if (nml_error /= 0) then
        nml_error = -1
     else
        nml_error =  1
     endif
     do while (nml_error > 0)
        read(nml_in, nml=ecosys_restore_nml, iostat=nml_error)
     end do
     if (nml_error == 0) then
        close(nml_in)
     end if
     
     write(stdout, *) 'ecosys_restore_nml :'
     write(stdout, ecosys_restore_nml)
  endif

  if (my_task == master_task) then
     ! FIXME(bja, 2014-10) assert(len(restore_short_names) == len(restore_filenames))
     write(stdout, *) "Found restore variables : "
     do t = 1, size(restore_short_names)
        if (len(trim(restore_short_names(t))) > 0) then
           write(stdout, *) "my_task = ", my_task, "      ", trim(restore_short_names(t)), " --> ", &
                trim(restore_filenames(t)), " [ ", trim(restore_file_varnames(t)), " ]"
        end if
     end do
  end if

  ! NOTE(bja, 2014-10-15) broadcast_array_char_1d looks like it
  ! should work but corrupts the data...
  do t = 1, size(restore_short_names)
     call broadcast_scalar(restore_short_names(t), master_task)
     call broadcast_scalar(restore_filenames(t), master_task)
     call broadcast_scalar(restore_file_varnames(t), master_task)
  end do

  call broadcast_scalar(spatial_variability_from_file, master_task)

  ! assign namelist variables to corresponding instance variables
  this%spatial_variability_from_file = spatial_variability_from_file

end subroutine read_namelist

!*****************************************************************************

subroutine initialize_restore_read_vars(this, restore_short_names, restore_filenames, &
     restore_file_varnames, ind_name_table)
  !
  ! Read the ecosys_restore namelist and broadcast to all
  ! processes. Store results in the ecosys_restore_vars%tracers(i)%restore_file_info
  !
  ! FIXME(bja, 2014-10) this%restore_any_tracer is set as a
  ! side-effect of this function! This isn't clear and is a Bad Thing (tm)
  !
  ! NOTE(bja, 2014-10) assumes that restore file is ALWAYS netcdf!
  use kinds_mod, only : char_len, int_kind, log_kind
  use passive_tracer_tools, only : ind_name_pair, name_to_ind
  use io_types, only : stdout
  use communicate, only : master_task, my_task
  use exit_mod, only : exit_POP, sigAbort

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type) :: this
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(in) :: restore_short_names
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(in) :: restore_filenames
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(in) :: restore_file_varnames
  type(ind_name_pair), dimension(:), intent(in) :: ind_name_table

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: t, tracer_index
  logical(log_kind) :: unknown_user_request
  character(len=char_len) :: file_format, message

  !-----------------------------------------------------------------------

  file_format = 'nc'

  this%restore_any_tracer = .false.
  unknown_user_request = .false.
  ! reinitialize any restore vars requested by the user
  do t = 1, size(restore_short_names)
     if (len(trim(restore_short_names(t))) > 0) then
        ! attempt to map the user specified tracer name to tracer index
        tracer_index = name_to_ind(trim(restore_short_names(t)), ind_name_table)
        if (tracer_index > 0) then
           this%restore_any_tracer = .true.
           call tracer_read_init(this%tracers(tracer_index)%restore_file_info, &
                mod_varname=restore_short_names(t), filename=restore_filenames(t), &
                file_varname=restore_file_varnames(t), file_fmt=file_format)
           this%tracers(tracer_index)%restore = .true.
           if (my_task == master_task) then
              write(stdout, *) "Setting up restoring for '", trim(restore_short_names(t)), "'"
           end if
        else
           unknown_user_request = .true.
           if (my_task == master_task) then
              write(stdout, *) "ERROR: Could not find user requested restore variable '", &
                   trim(restore_short_names(t)), "'"
           end if
        end if
     end if
  end do

  if (unknown_user_request) then
     write(message, *) "ecosys_restore%initialize_restore_read_vars - ",      &
                       "encountered unknown user specified tracer names.",    &
                       "See ocn.log for details."
     call exit_POP(sigAbort, message)

  end if

end subroutine initialize_restore_read_vars

!*****************************************************************************

subroutine define_tavg_fields(this)
  !
  ! define the tavg fields for restoring data
  !
  use kinds_mod, only : char_len, int_kind
  use tavg, only : define_tavg_field

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type) :: this

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer (int_kind) :: n
  character (char_len) :: short_name
  character (char_len) :: long_name

  !-----------------------------------------------------------------------

  do n = 1, size(this%tracers)
     if (this%tracers(n)%restore) then
        short_name = trim(this%tracers(n)%restore_file_info%file_varname) // "_RESTORE"
        long_name = trim(this%tracers(n)%restore_file_info%file_varname) // " Restoring"

        call define_tavg_field(this%tracers(n)%tavg_restore_index, short_name, 3, &
             long_name=long_name, &
             units='mmol/m^3', &
             grid_loc='3111', &
             coordinates='TLONG TLAT z_t time')
     end if
  end do

end subroutine define_tavg_fields

!*****************************************************************************

subroutine initialize_restoring_timescale(this, nml_filename, nml_in, zt)
  !
  ! Initialize the spatially varying restoring timescale by 
  !
  use kinds_mod, only : i4
  use POP_KindsMod, only : POP_r8
  use io_types, only : stdout
  use communicate, only : master_task, my_task

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type) :: this
  character(len=*), intent(in) :: nml_filename
  integer(i4), intent(in) :: nml_in
  real (POP_r8), dimension(km), intent(in) :: zt

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------

  if (this%restore_any_tracer) then
     if (this%spatial_variability_from_file) then
        call this%timescale_file%init(nml_filename, nml_in)
     else
        call this%timescale_interp%init(nml_filename, nml_in, zt)
     end if
     if (my_task == master_task) then
        write(stdout, *) "Restoring timescale set."
     end if
  end if

end subroutine initialize_restoring_timescale

!*****************************************************************************

subroutine read_restoring_fields(this, LAND_MASK)
  !
  !  load restoring fields if required
  !
  use kinds_mod, only : int_kind
  use constants, only : c0
  use blocks, only : nx_block, ny_block
  use domain_size, only : km
  use domain, only : nblocks_clinic
  use prognostic, only : max_blocks_clinic
  use grid, only : KMT
  use passive_tracer_tools, only : read_field
  use io_types, only : stdout

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type) :: this
  logical (log_kind), dimension(:,:,:), intent(in) :: LAND_MASK

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer (int_kind) :: i,j    ! indoces for looping over horizontal grid
  integer (int_kind) :: k      ! index for looping over levels
  integer (int_kind) :: iblock ! index for looping over blocks
  integer (int_kind) :: tracer_index

  !-----------------------------------------------------------------------

  do tracer_index = 1, size(this%tracers)
     associate(tracer => this%tracers(tracer_index))

       if (tracer%restore) then
!!$          write(stdout, *) "Reading restore data for ", trim(tracer%restore_file_info%mod_varname)
          allocate(tracer%data(nx_block, ny_block, km, max_blocks_clinic))

          call read_field(tracer%restore_file_info%file_fmt, &
               tracer%restore_file_info%filename, &
               tracer%restore_file_info%file_varname, &
               tracer%data)

          do iblock = 1, nblocks_clinic
             do k = 1, km
                do j=1,ny_block
                   do i=1,nx_block
                      if (LAND_MASK(i,j,iblock) .and. (k.le.KMT(i,j,iblock))) then
                         tracer%data(i,j,k,iblock) = tracer%data(i,j,k,iblock) * tracer%restore_file_info%scale_factor
                      else
                         tracer%data(i,j,k,iblock) = c0
                      end if
                   end do
                end do
             end do
          end do
       endif
     end associate
  end do

end subroutine read_restoring_fields

!*****************************************************************************

subroutine restore_variable(this, tracer_index, vert_level, block_id, local_data, &
     restore_data)
  !
  !  restore a variable if required
  !
  use kinds_mod, only : r8, int_kind, log_kind
  use constants, only : c0
  use blocks, only : nx_block, ny_block
  use domain_size, only : km
  use passive_tracer_tools, only : tracer_read

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type) :: this
  integer(int_kind), intent(in) :: tracer_index
  integer(int_kind), intent(in) :: vert_level
  integer (int_kind), intent(in) :: block_id
  real(r8), dimension(nx_block, ny_block), intent(in) :: local_data
  real(r8), dimension(nx_block, ny_block), intent(out) :: restore_data
  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------

  associate( &
       tracer => this%tracers(tracer_index), &
       rtau => this%timescale_file%restore_rtau, &
       restore_max_level => this%timescale_file%restore_max_level, &
       inv_restoring_time_scale => this%timescale_interp%inv_restoring_time_scale &
       )

    if (tracer%restore) then
       if (this%spatial_variability_from_file) then
          restore_data = rtau(:, :, block_id) * &
               merge((tracer%data(:, :, vert_level, block_id) - local_data), &
               c0, vert_level <= restore_max_level(:, :, block_id))
       else
          restore_data = (tracer%data(:, :, vert_level, block_id) - local_data) * inv_restoring_time_scale(vert_level)
       endif
    else
       restore_data = c0
    endif

  end associate

end subroutine restore_variable

!*****************************************************************************

subroutine accumulate_tavg(this, tracer_index, vert_level, block_id, restore_local)
  !
  ! accumulate tavg field info for the specified tracer if approprate
  !
  use kinds_mod, only : r8, int_kind, log_kind
  use tavg, only : accumulate_tavg_field

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type) :: this
  integer(int_kind), intent(in) :: tracer_index
  integer(int_kind), intent(in) :: vert_level
  integer (int_kind), intent(in) :: block_id
  real (r8), intent(in) :: restore_local(:, :)

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------

  if (this%tracers(tracer_index)%tavg_restore_index > 0) then
     call accumulate_tavg_field(restore_local, &
          this%tracers(tracer_index)%tavg_restore_index, &
          block_id, vert_level)
  end if

end subroutine accumulate_tavg

!*****************************************************************************

end module ecosys_restore_mod
