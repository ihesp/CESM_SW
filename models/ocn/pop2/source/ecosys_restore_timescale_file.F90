module ecosys_restore_timescale_file
  !
  ! Control spatial variability of restoring timescale with an input file
  !
  use kinds_mod, only : r8, int_kind

  implicit none

  private

  type, public :: ecosys_restore_timescale_file_type
     ! inverse restoring timescale for variable interior restoring
     real (r8), dimension(:,:,:), allocatable, public :: restore_rtau
     ! maximum level for applying variable interior restoring
     integer (int_kind), dimension(:,:,:), allocatable, public :: restore_max_level

   contains
     procedure, public :: init
     procedure, private :: read_namelist
     procedure, private :: read_restoring_timescale_from_file

  end type ecosys_restore_timescale_file_type

contains

!*****************************************************************************

subroutine init(this, nml_filename, nml_in)

  use kinds_mod, only : char_len, i4

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_file_type) :: this
  character(len=*), intent(in) :: nml_filename
  integer(i4), intent(in) :: nml_in

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  character(char_len) :: file_format
  character(char_len) :: file_name

  !-----------------------------------------------------------------------
  call this%read_namelist(nml_filename, nml_in, file_name, file_format)
  call this%read_restoring_timescale_from_file(file_name, file_format)
  
end subroutine init

!*****************************************************************************

subroutine read_namelist(this, nml_filename, nml_in, file_name, file_format)

  use kinds_mod, only : int_kind, i4, char_len
  use io_types, only : stdout
  use communicate, only : master_task, my_task
  use broadcast, only : broadcast_scalar

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_file_type) :: this
  character(len=*), intent(in) :: nml_filename
  integer(i4), intent(in) :: nml_in
  character(char_len), intent(out) :: file_format
  character(char_len), intent(out) :: file_name

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: nml_error
  character(char_len) :: restore_timescale_file_format
  character(char_len) :: restore_timescale_file_name

  !-----------------------------------------------------------------------
  namelist /ecosys_restore_timescale_file_nml/ &
       restore_timescale_file_format, &
       restore_timescale_file_name

  restore_timescale_file_format = ''
  restore_timescale_file_name = ''

  if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     if (nml_error /= 0) then
        nml_error = -1
     else
        nml_error =  1
     endif
     do while (nml_error > 0)
        read(nml_in, nml=ecosys_restore_timescale_file_nml, iostat=nml_error)
     end do
     if (nml_error == 0) then
        close(nml_in)
     end if
     
     write(stdout, *) 'ecosys_restore_timescale_file_nml :'
     write(stdout, ecosys_restore_timescale_file_nml)
  endif

  call broadcast_scalar(restore_timescale_file_format, master_task)
  call broadcast_scalar(restore_timescale_file_name, master_task)

  file_format = restore_timescale_file_format
  file_name = restore_timescale_file_name

end subroutine read_namelist

!*****************************************************************************

subroutine read_restoring_timescale_from_file(this, file_format, file_name)
  !
  ! Initialize the spatially variable restoring timescale from the the
  ! user specified file
  !
  use kinds_mod, only : char_len
  use blocks, only : nx_block, ny_block
  use prognostic, only : max_blocks_clinic
  use time_management, only : seconds_in_day
  use passive_tracer_tools, only : read_field

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_file_type) :: this
  character(len=char_len), intent(in) :: file_format
  character(len=char_len), intent(in) :: file_name

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------

    allocate(this%restore_rtau(nx_block, ny_block, max_blocks_clinic))
    allocate(this%restore_max_level(nx_block, ny_block, max_blocks_clinic))

    call read_field(file_format, file_name, &
         'NUTR_RESTORE_MAX_LEVEL', this%restore_rtau)

    this%restore_max_level = nint(this%restore_rtau)

    call read_field(file_format, file_name, &
         'NUTR_RESTORE_RTAU', this%restore_rtau)

    this%restore_rtau = this%restore_rtau / seconds_in_day ! convert days to secs

  end subroutine read_restoring_timescale_from_file

!*****************************************************************************

end module
