module ecosys_restore_timescale_interp
  !
  ! Control spatial variability of restoring timescale with
  ! interpolation based on depth
  !

  use kinds_mod, only : r8, log_kind, int_kind
  use domain_size, only : km
  use constants, only : c0, c2, c1000

  implicit none

  private
  public :: interp_init
  private :: read_namelist
  private :: interpolate_restoring_timescale

  type, public :: ecosys_restore_timescale_interp_type
     real (r8), dimension(km) :: &
          inv_restoring_time_scale ! inverse restoring time scale for nutrients (1/secs)

     real (r8) :: rest_time_inv_surf ! inverse restoring timescale at surface
     real (r8) :: rest_time_inv_deep ! inverse restoring timescale at depth
     real (r8) :: rest_z0 ! shallow end of transition regime
     real (r8) :: rest_z1 ! deep end of transition regime

!   contains
!     procedure, public :: init
!     procedure, private :: read_namelist
!     procedure, private :: interpolate_restoring_timescale

  end type ecosys_restore_timescale_interp_type

  real (r8), parameter, private :: default_rest_time_inv_surf = c0
  real (r8), parameter, private :: default_rest_time_inv_deep = c0
  real (r8), parameter, private :: default_rest_z0 = c1000
  real (r8), parameter, private :: default_rest_z1 = c2 * c1000

contains


!*****************************************************************************

subroutine interp_init(this, nml_filename, nml_in, zt)

  use kinds_mod, only : i4
  use POP_KindsMod, only : POP_r8

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  type(ecosys_restore_timescale_interp_type) ,intent(inout):: this
  character(len=*), intent(in) :: nml_filename
  integer(i4), intent(in) :: nml_in
  real (POP_r8), dimension(km), intent(in) :: zt

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  call read_namelist(this,nml_filename, nml_in)
  call interpolate_restoring_timescale(this,zt)
  
end subroutine interp_init

!*****************************************************************************

subroutine read_namelist(this, nml_filename, nml_in)

  use kinds_mod, only : r8, i4
  use io_types, only : stdout
  use communicate, only : master_task, my_task
  use broadcast, only : broadcast_scalar
  use exit_mod, only : exit_POP, sigAbort
  use io_tools, only : document

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  type(ecosys_restore_timescale_interp_type) ,intent(inout):: this
  character(len=*), intent(in) :: nml_filename
  integer(i4), intent(in) :: nml_in

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: nml_error
  real(r8) :: rest_time_inv_surf, rest_time_inv_deep, rest_z0, rest_z1

  !-----------------------------------------------------------------------
  namelist /ecosys_restore_timescale_interp_nml/ &
       rest_time_inv_surf, &
       rest_time_inv_deep, &
       rest_z0, rest_z1

  rest_time_inv_surf = default_rest_time_inv_surf
  rest_time_inv_deep = default_rest_time_inv_deep
  rest_z0 = default_rest_z0
  rest_z1 = default_rest_z1

  if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     if (nml_error /= 0) then
        nml_error = -1
     else
        nml_error =  1
     endif
     do while (nml_error > 0)
        read(nml_in, nml=ecosys_restore_timescale_interp_nml, iostat=nml_error)
     end do
     if (nml_error == 0) then
        close(nml_in)
     end if
     
     write(stdout, *) 'ecosys_restore_timescale_interp_nml :'
     write(stdout, ecosys_restore_timescale_interp_nml)
  endif

  call broadcast_scalar(nml_error, master_task)
  if (nml_error /= 0) then
     call document("ecosys_restore_timescale_interp::read_namelist", &
          'ecosys_restore_timescale_interp_nml not found')
     call exit_POP(sigAbort, 'stopping in ' /&
          &/ "ecosys_restore_timescale_interp::read_namelist")
  endif
  call broadcast_scalar(rest_time_inv_surf, master_task)
  call broadcast_scalar(rest_time_inv_deep, master_task)
  call broadcast_scalar(rest_z0, master_task)
  call broadcast_scalar(rest_z1, master_task)

  this%rest_time_inv_surf = rest_time_inv_surf
  this%rest_time_inv_deep = rest_time_inv_deep
  this%rest_z0 = rest_z0
  this%rest_z1 = rest_z1

end subroutine read_namelist

!*****************************************************************************

subroutine interpolate_restoring_timescale(this, zt)
  !
  ! Initialize the spatial variability of the restoring time scale
  ! with an interpolation.
  !
  use kinds_mod, only : int_kind
  use POP_KindsMod, only : POP_r8
  use constants, only : p5
  use domain_size, only : km

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  type(ecosys_restore_timescale_interp_type) :: this
  real (POP_r8), dimension(km) :: zt

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: k

  do k = 1, km
     if (zt(k) < this%rest_z0) then
        this%inv_restoring_time_scale(k) = this%rest_time_inv_surf
     else if (zt(k) > this%rest_z1) then
        this%inv_restoring_time_scale(k) = this%rest_time_inv_deep
     else if (this%rest_z1 == this%rest_z0) then
        this%inv_restoring_time_scale(k) = this%rest_time_inv_surf + p5 * &
             (this%rest_time_inv_deep - this%rest_time_inv_surf)
     else
        this%inv_restoring_time_scale(k) = this%rest_time_inv_surf + &
             (zt(k) - this%rest_z0) / (this%rest_z1 - this%rest_z0) * &
             (this%rest_time_inv_deep - this%rest_time_inv_surf)
     endif
  end do
end subroutine interpolate_restoring_timescale

!*****************************************************************************

end module ecosys_restore_timescale_interp
