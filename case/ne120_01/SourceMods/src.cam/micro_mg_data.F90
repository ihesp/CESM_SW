module micro_mg_data

!
! Packing and time averaging for the MG interface.
!
! Use is as follows:
!
! 1) Figure out which columns will do averaging (mgncol) and the number of
!    levels where the microphysics will run (nlev).
!
! 2) Create an MGPacker object and assign it as follows:
!
!      packer = MGPacker(pcols, pver, mgcols, top_lev)
!
!    Where [pcols, pver] is the shape of the ultimate input/output arrays
!    that are defined at level midpoints.
!
! 3) Create a post-processing array of type MGPostProc:
!
!      post_proc = MGPostProc(packer)
!
! 4) Add pairs of pointers for packed and unpacked representations, already
!    associated with buffers of the correct dimensions:
!
!      call post_proc%add_field(unpacked_pointer, packed_pointer, &
!             fillvalue, accum_mean)
!
!    The third value is the default value used to "unpack" for points with
!    no "packed" part, and the fourth value is the method used to
!    accumulate values over time steps. These two arguments can be omitted,
!    in which case the default value will be 0 and the accumulation method
!    will take the mean.
!
! 5) Use the packed fields in MG, and for each MG iteration, do:
!
!      call post_proc%accumulate()
!
! 6) Perform final accumulation and scatter values into the unpacked arrays:
!
!      call post_proc%process_and_unpack()
!
! 7) Destroy the object when complete:
!
!      call post_proc%finalize()
!
! Caveat: MGFieldPostProc will hit a divide-by-zero error if you try to
!         take the mean over 0 steps.
!

! This include header defines CPP macros that only have an effect for debug
! builds.
#include "shr_assert.h"

use shr_kind_mod, only: r8 => shr_kind_r8
use shr_log_mod, only: &
     errMsg => shr_log_errMsg, &
     OOBMsg => shr_log_OOBMsg
use shr_sys_mod, only: shr_sys_abort

implicit none
private

public :: MGPacker,chen_new_MGPacker,chen_new_MGPostProc
public :: MGFieldPostProc
public :: accum_null
public :: accum_mean
public :: MGPostProc
public :: accumulate,process_and_unpack,finalize
public :: add_field

type :: MGPacker
   ! Unpacked array dimensions.
   integer :: pcols
   integer :: pver
   ! Calculated packed dimensions, stored for convenience.
   integer :: mgncol
   integer :: nlev
   ! Which columns are packed.
   integer, allocatable :: mgcols(:)
   ! Topmost level to copy into the packed array.
   integer :: top_lev
end type MGPacker
 
public:: xpack,pack_interface,xunpack

interface xpack
  module procedure pack_1D
  module procedure pack_2D
  module procedure pack_3D
end interface xpack

interface xunpack
  module procedure unpack_1D
  module procedure unpack_1D_array_fill
  module procedure unpack_2D
  module procedure unpack_2D_array_fill
  module procedure unpack_3D
  module procedure unpack_3D_array_fill
end interface xunpack

!interface finalize
!  module procedure MGPacker_finalize
!  module procedure MGFieldPostProc_finalize
!end interface finalize
  


!interface MGPacker
!   module procedure new_MGPacker
!end interface

! Enum for time accumulation/averaging methods.
integer, parameter :: accum_null = 0
integer, parameter :: accum_mean = 1

type :: MGFieldPostProc
   integer :: accum_method = -1
   integer :: rank = -1
   integer :: num_steps = 0
   real(r8) :: fillvalue = 0._r8
   real(r8), pointer :: unpacked_1D(:) => null()
   real(r8), pointer :: packed_1D(:) => null()
   real(r8), allocatable :: buffer_1D(:)
   real(r8), pointer :: unpacked_2D(:,:) => null()
   real(r8), pointer :: packed_2D(:,:) => null()
   real(r8), allocatable :: buffer_2D(:,:)
end type MGFieldPostProc
! contains
!   procedure :: accumulate => MGFieldPostProc_accumulate
!   procedure :: process_and_unpack => MGFieldPostProc_process_and_unpack
!   procedure :: unpack_only => MGFieldPostProc_unpack_only
!   procedure :: finalize => MGFieldPostProc_finalize
!end type MGFieldPostProc

interface accumulate
  module procedure MGFieldPostProc_accumulate
  module procedure MGPostProc_accumulate
end interface accumulate

interface process_and_unpack
  module procedure MGFieldPostProc_process_and_unpack
  module procedure MGPostProc_process_and_unpack
end interface process_and_unpack

interface  unpack_only
  module procedure  MGFieldPostProc_unpack_only
  module procedure MGPostProc_unpack_only
end interface   unpack_only

interface MGFieldPostProcs
   module procedure MGFieldPostProc_1D
   module procedure MGFieldPostProc_2D
end interface MGFieldPostProcs

#define VECTOR_NAME MGFieldPostProcVec
#define TYPE_NAME type(MGFieldPostProc)
#define THROW(string) call shr_sys_abort(string)

public :: VECTOR_NAME

#include "dynamic_vector_typedef.inc"

type MGPostProc
   type(MGPacker) :: packer
   type(MGFieldPostProcVec) :: field_procs
end type MGPostProc
 interface add_field
    module procedure add_field_1D
    module procedure add_field_2D
 end interface
    
!! contains
!!   procedure, private :: add_field_1D
!!   procedure, private :: add_field_2D
!!   generic :: add_field => add_field_1D, add_field_2D
!!   procedure :: accumulate => MGPostProc_accumulate
! interface  process_and_unpack
! end interface  process_and_unpack
!!   procedure :: process_and_unpack => MGPostProc_process_and_unpack
!!   procedure :: unpack_only => MGPostProc_unpack_only
 interface finalize
  module procedure MGPostProc_finalize
 end interface finalize
!!   procedure :: finalize => MGPostProc_finalize

!!   procedure, private :: MGPostProc_copy
!!   generic :: assignment(=) => MGPostProc_copy
!!end type MGPostProc


!!interface MGPostProc
!!   module procedure new_MGPostProc
!!end interface MGPostProc

interface assignment (=)
   module procedure MGPostProc_copy
end interface assignment (=)

interface assignment (=)
   module procedure MGFieldPostProc_copy
end interface assignment (=)

contains

subroutine chen_new_MGPacker(pcols, pver, mgcols, top_lev,new_MGPacker)
  integer, intent(in) :: pcols, pver
  integer, intent(in) :: mgcols(:)
  integer, intent(in) :: top_lev

  type(MGPacker),intent(inout) :: new_MGPacker

  new_MGPacker%pcols = pcols
  new_MGPacker%pver = pver
  new_MGPacker%mgncol = size(mgcols)
  new_MGPacker%nlev = pver - top_lev + 1

  allocate(new_MGPacker%mgcols(new_MGPacker%mgncol))
  new_MGPacker%mgcols = mgcols
  new_MGPacker%top_lev = top_lev

end subroutine chen_new_MGPacker

! Rely on the fact that intent(out) forces the compiler to deallocate all
! allocatable components and restart the type from scratch. Although
! compiler support for finalization varies, this seems to be one of the few
! cases where all major compilers are reliable, and humans are not.
subroutine MGPacker_finalize(self)
!  class(MGPacker), intent(out) :: self
  type(MGPacker), intent(out) :: self
end subroutine MGPacker_finalize

function pack_1D(self, unpacked) result(packed)
  type(MGPacker), intent(in) :: self
!  class(MGPacker), intent(in) :: self
  real(r8), intent(in) :: unpacked(:)

  real(r8) :: packed(self%mgncol)

  SHR_ASSERT(size(unpacked) == self%pcols, errMsg(__FILE__, __LINE__))
  !write(*,*) "printed by Asher in pack_1D",self%mgcols
  packed = unpacked(self%mgcols)
  !write(*,*) "printed by Asher at the end of pack_1D"

end function pack_1D

! Separation of pack and pack_interface is to workaround a PGI bug.
function pack_2D(self, unpacked) result(packed)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: unpacked(:,:)

  real(r8) :: packed(self%mgncol,self%nlev)

  SHR_ASSERT(size(unpacked, 1) == self%pcols, errMsg(__FILE__, __LINE__))

  packed = unpacked(self%mgcols,self%top_lev:)

end function pack_2D

function pack_interface(self, unpacked) result(packed)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: unpacked(:,:)

  real(r8) :: packed(self%mgncol,self%nlev+1)

  packed = unpacked(self%mgcols,self%top_lev:)

end function pack_interface

function pack_3D(self, unpacked) result(packed)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: unpacked(:,:,:)

  real(r8) :: packed(self%mgncol,self%nlev,size(unpacked, 3))

  SHR_ASSERT(size(unpacked,1) == self%pcols, errMsg(__FILE__, __LINE__))

  packed = unpacked(self%mgcols,self%top_lev:,:)

end function pack_3D

function unpack_1D(self, packed, fill) result(unpacked)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:)
  real(r8), intent(in) :: fill

  real(r8) :: unpacked(self%pcols)

  SHR_ASSERT(size(packed) == self%mgncol, errMsg(__FILE__, __LINE__))

  unpacked = fill
  unpacked(self%mgcols) = packed

end function unpack_1D

function unpack_1D_array_fill(self, packed, fill) result(unpacked)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:)
  real(r8), intent(in) :: fill(:)

  real(r8) :: unpacked(self%pcols)

  SHR_ASSERT(size(packed) == self%mgncol, errMsg(__FILE__, __LINE__))

  unpacked = fill
  unpacked(self%mgcols) = packed

end function unpack_1D_array_fill

function unpack_2D(self, packed, fill) result(unpacked)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:,:)
  real(r8), intent(in) :: fill

  real(r8) :: unpacked(self%pcols,self%pver+size(packed, 2)-self%nlev)

  SHR_ASSERT(size(packed, 1) == self%mgncol, errMsg(__FILE__, __LINE__))

  unpacked = fill
  unpacked(self%mgcols,self%top_lev:) = packed

end function unpack_2D

function unpack_2D_array_fill(self, packed, fill) result(unpacked)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:,:)
  real(r8), intent(in) :: fill(:,:)

  real(r8) :: unpacked(self%pcols,self%pver+size(packed, 2)-self%nlev)

  SHR_ASSERT(size(packed, 1) == self%mgncol, errMsg(__FILE__, __LINE__))

  unpacked = fill
  unpacked(self%mgcols,self%top_lev:) = packed

end function unpack_2D_array_fill

function unpack_3D(self, packed, fill) result(unpacked)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:,:,:)
  real(r8), intent(in) :: fill

  real(r8) :: unpacked(self%pcols,self%pver,size(packed, 3))

  SHR_ASSERT(size(packed, 1) == self%mgncol, errMsg(__FILE__, __LINE__))

  unpacked = fill
  unpacked(self%mgcols,self%top_lev:,:) = packed

end function unpack_3D

function unpack_3D_array_fill(self, packed, fill) result(unpacked)
!  class(MGPacker), intent(in) :: self
  type(MGPacker), intent(in) :: self
  real(r8), intent(in) :: packed(:,:,:)
  real(r8), intent(in) :: fill(:,:,:)

  real(r8) :: unpacked(self%pcols,self%pver,size(packed, 3))

  SHR_ASSERT(size(packed, 1) == self%mgncol, errMsg(__FILE__, __LINE__))

  unpacked = fill
  unpacked(self%mgcols,self%top_lev:,:) = packed

end function unpack_3D_array_fill

function MGFieldPostProc_1D(unpacked_ptr, packed_ptr, fillvalue, &
     accum_method) result(field_proc)
  real(r8), pointer, intent(in) :: unpacked_ptr(:)
  real(r8), pointer, intent(in) :: packed_ptr(:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  type(MGFieldPostProc) :: field_proc

  field_proc%rank = 1
  field_proc%unpacked_1D => unpacked_ptr
  field_proc%packed_1D => packed_ptr
  if (present(fillvalue)) then
     field_proc%fillvalue = fillvalue
  else
     field_proc%fillvalue = 0._r8
  end if
  if (present(accum_method)) then
     field_proc%accum_method = accum_method
  else
     field_proc%accum_method = accum_mean
  end if

end function MGFieldPostProc_1D

function MGFieldPostProc_2D(unpacked_ptr, packed_ptr, fillvalue, &
     accum_method) result(field_proc)
  real(r8), pointer, intent(in) :: unpacked_ptr(:,:)
  real(r8), pointer, intent(in) :: packed_ptr(:,:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  type(MGFieldPostProc) :: field_proc
  !print *, "xduan: postproc_2d", loc(unpacked_ptr), loc(packed_ptr), loc(field_proc)
  field_proc%rank = 2
  field_proc%unpacked_2D => unpacked_ptr
  field_proc%packed_2D => packed_ptr
  !print *, "xduan: postproc_2d, outer done"
  if (present(fillvalue)) then
     !print *, "xduan: postproc_2d fillval", fillvalue
     field_proc%fillvalue = fillvalue
  else
     field_proc%fillvalue = 0._r8
  end if
  if (present(accum_method)) then
     !print *, "xduan: postproc_2d accum", accum_method
     field_proc%accum_method = accum_method
  else
     field_proc%accum_method = accum_mean
  end if
  !print *, "xduan: postproc_2d, inner done"
end function MGFieldPostProc_2D


! Use the same intent(out) trick as for MGPacker, which is actually more
! useful here.
subroutine MGFieldPostProc_finalize(self)
!  class(MGFieldPostProc), intent(out) :: self
  type(MGFieldPostProc), intent(out) :: self
end subroutine MGFieldPostProc_finalize

subroutine MGFieldPostProc_accumulate(self)
!  class(MGFieldPostProc), intent(inout) :: self
  type(MGFieldPostProc), intent(inout) :: self

  select case (self%accum_method)
  case (accum_null)
     ! "Null" method does nothing.
  case (accum_mean)
     ! Allocation is done on the first accumulation step to allow the
     ! MGFieldPostProc to be copied after construction without copying the
     ! allocated array (until this function is first called).
     self%num_steps = self%num_steps + 1
     select case (self%rank)
     case (1)
        SHR_ASSERT(associated(self%packed_1D), errMsg(__FILE__, __LINE__))
        if (.not. allocated(self%buffer_1D)) then
           allocate(self%buffer_1D(size(self%packed_1D)))
           self%buffer_1D = 0._r8
        end if
        self%buffer_1D = self%buffer_1D + self%packed_1D
     case (2)
        SHR_ASSERT(associated(self%packed_2D), errMsg(__FILE__, __LINE__))
        if (.not. allocated(self%buffer_2D)) then
           ! Awkward; in F2008 can be replaced by source/mold.
           allocate(self%buffer_2D(&
                size(self%packed_2D, 1),size(self%packed_2D, 2)))
           self%buffer_2D = 0._r8
        end if
        self%buffer_2D = self%buffer_2D + self%packed_2D
     case default
        call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
             " Unsupported rank for MGFieldPostProc accumulation.")
     end select
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unrecognized MGFieldPostProc accumulation method.")
  end select

end subroutine MGFieldPostProc_accumulate

subroutine MGFieldPostProc_process_and_unpack(self, packer)
!  class(MGFieldPostProc), intent(inout) :: self
  type(MGFieldPostProc), intent(inout) :: self
!  class(MGPacker), intent(in) :: packer
  type(MGPacker), intent(in) :: packer

  select case (self%accum_method)
  case (accum_null)
     ! "Null" method just leaves the value as the last time step, so don't
     ! actually need to do anything.
  case (accum_mean)
     select case (self%rank)
     case (1)
        SHR_ASSERT(associated(self%packed_1D), errMsg(__FILE__, __LINE__))
        self%packed_1D = self%buffer_1D/self%num_steps
     case (2)
        SHR_ASSERT(associated(self%packed_2D), errMsg(__FILE__, __LINE__))
        self%packed_2D = self%buffer_2D/self%num_steps
     case default
        call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
             " Unsupported rank for MGFieldPostProc accumulation.")
     end select
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unrecognized MGFieldPostProc accumulation method.")
  end select

!  call self%unpack_only(packer)
  call unpack_only(self,packer)

end subroutine MGFieldPostProc_process_and_unpack

subroutine MGFieldPostProc_unpack_only(self, packer)
!  class(MGFieldPostProc), intent(inout) :: self
  type(MGFieldPostProc), intent(inout) :: self
!  class(MGPacker), intent(in) :: packer
  type(MGPacker), intent(in) :: packer

  select case (self%rank)
  case (1)
     SHR_ASSERT(associated(self%unpacked_1D), errMsg(__FILE__, __LINE__))
     self%unpacked_1D = xunpack(packer,self%packed_1D, self%fillvalue)
 !    self%unpacked_1D = packer%xunpack(self%packed_1D, self%fillvalue)
  case (2)
     SHR_ASSERT(associated(self%unpacked_2D), errMsg(__FILE__, __LINE__))
 !    self%unpacked_2D = packer%unpack(self%packed_2D, self%fillvalue)
     self%unpacked_2D = xunpack(packer,self%packed_2D, self%fillvalue)
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unsupported rank for MGFieldPostProc unpacking.")
  end select

end subroutine MGFieldPostProc_unpack_only

#include "dynamic_vector_procdef.inc"

subroutine chen_new_MGPostProc(packer,post_proc)
  type(MGPacker), intent(in) :: packer

  type(MGPostProc),intent(inout) :: post_proc

  post_proc%packer = packer
!  call post_proc%field_procs%clear()
  call clear(post_proc%field_procs)

end subroutine chen_new_MGPostProc

! Can't use the same intent(out) trick, because PGI doesn't get the
! recursive deallocation right.
subroutine MGPostProc_finalize(self)
!  class(MGPostProc), intent(inout) :: self
  type(MGPostProc), intent(inout) :: self

  integer :: i

!  call self%packer%finalize()
  !call finalize(self%packer)
  call MGPacker_finalize(self%packer)
!  do i = 1, self%field_procs%vsize()
  do i = 1, vsize(self%field_procs)
!     call self%field_procs%data(i)%finalize()
     !call finalize(self%field_procs%data(i))
     call MGFieldPostProc_finalize(self%field_procs%data(i))
  end do
!  call self%field_procs%clear()
  call clear(self%field_procs)
!  call self%field_procs%shrink_to_fit()
  call shrink_to_fit(self%field_procs)

end subroutine MGPostProc_finalize


subroutine add_field_1D(self, unpacked_ptr, packed_ptr, fillvalue, &
     accum_method)
!  class(MGPostProc), intent(inout) :: self
  type(MGPostProc), intent(inout) :: self
  real(r8), pointer, intent(in) :: unpacked_ptr(:)
  real(r8), pointer, intent(in) :: packed_ptr(:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  !print *, "xduan: add field 1d"
!  call self%field_procs%push_back(MGFieldPostProc(unpacked_ptr, &
  call push_back(self%field_procs,MGFieldPostProcs(unpacked_ptr, &
       packed_ptr, fillvalue, accum_method))

end subroutine add_field_1D

subroutine MGFieldPostProc_2D_routine(field_proc, unpacked_ptr, packed_ptr, &
     fillvalue, accum_method)
  real(r8), pointer, intent(in) :: unpacked_ptr(:,:)
  real(r8), pointer, intent(in) :: packed_ptr(:,:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  type(MGFieldPostProc), intent(inout) :: field_proc
  !print *, "xduan: postproc_2d", loc(unpacked_ptr), loc(packed_ptr), loc(field_proc)
  field_proc%rank = 2
  field_proc%unpacked_2D => unpacked_ptr
  field_proc%packed_2D => packed_ptr
  !print *, "xduan: postproc_2d, outer done"
  if (present(fillvalue)) then
     !print *, "xduan: postproc_2d fillval", fillvalue
     field_proc%fillvalue = fillvalue
  else
     field_proc%fillvalue = 0._r8
  end if
  if (present(accum_method)) then
     !print *, "xduan: postproc_2d accum", accum_method
     field_proc%accum_method = accum_method
  else
     field_proc%accum_method = accum_mean
  end if
  !print *, "xduan: postproc_2d, inner done"
end subroutine MGFieldPostProc_2D_routine

subroutine add_field_2D(self, unpacked_ptr, packed_ptr, fillvalue, &
     accum_method)
!  class(MGPostProc), intent(inout) :: self
  type(MGPostProc), intent(inout) :: self
  real(r8), pointer, intent(in) :: unpacked_ptr(:,:)
  real(r8), pointer, intent(in) :: packed_ptr(:,:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  type(MGFieldPostProc) :: field_proc
  !print *, "xduan: add field 2d"
!  call self%field_procs%push_back(MGFieldPostProc(unpacked_ptr, &
  !field_proc = MGFieldPostProcs(unpacked_ptr, &
  call MGFieldPostProc_2D_routine(field_proc, unpacked_ptr, &
       packed_ptr, fillvalue, accum_method)
  !print *, "field_proc_done"
  ! call resize(self%field_procs, self%field_procs%vec_size + 1)
  ! self%field_procs%data(self%field_procs%vec_size) = field_proc
  call push_back(self%field_procs, field_proc)
  !print *, "push back done"
end subroutine add_field_2D

subroutine MGPostProc_accumulate(self)
!  class(MGPostProc), intent(inout) :: self
  type(MGPostProc), intent(inout) :: self

  integer :: i

!  do i = 1, self%field_procs%vsize()
  do i = 1, vsize(self%field_procs)
!     call self%field_procs%data(i)%accumulate()
     call accumulate(self%field_procs%data(i))
  end do

end subroutine MGPostProc_accumulate

subroutine MGPostProc_process_and_unpack(self)
!  class(MGPostProc), intent(inout) :: self
  type(MGPostProc), intent(inout) :: self

  integer :: i

!  do i = 1, self%field_procs%vsize()
  do i = 1, vsize(self%field_procs)
!     call self%field_procs%data(i)%process_and_unpack(self%packer)
     call process_and_unpack(self%field_procs%data(i),self%packer)
  end do

end subroutine MGPostProc_process_and_unpack

subroutine MGPostProc_unpack_only(self)
!  class(MGPostProc), intent(inout) :: self
  type(MGPostProc), intent(inout) :: self

  integer :: i

!  do i = 1, self%field_procs%vsize()
  do i = 1, vsize(self%field_procs)
!     call self%field_procs%data(i)%unpack_only(self%packer)
     call unpack_only(self%field_procs%data(i),self%packer)
  end do

end subroutine MGPostProc_unpack_only

! This is necessary only to work around Intel/PGI bugs.
subroutine MGPostProc_copy(lhs, rhs)
!  class(MGPostProc), intent(out) :: lhs
  type(MGPostProc), intent(out) :: lhs
  type(MGPostProc), intent(in) :: rhs
  !print *, "MGPostProc copy"
  lhs%packer = rhs%packer
  lhs%field_procs = rhs%field_procs
  !print *, "MGPostProc copy exitting"
end subroutine MGPostProc_copy

subroutine MGFieldPostProc_copy(lhs, rhs)
  !  class(MGFieldPostProc), intent(out) :: lhs
  type(MGFieldPostProc), intent(out) :: lhs
  type(MGFieldPostProc), intent(in) :: rhs
  ! integer :: accum_method = -1
  ! integer :: rank = -1
  ! integer :: num_steps = 0
  ! real(r8) :: fillvalue = 0._r8
  ! real(r8), pointer :: unpacked_1D(:) => null()
  ! real(r8), pointer :: packed_1D(:) => null()
  ! real(r8), allocatable :: buffer_1D(:)
  ! real(r8), pointer :: unpacked_2D(:,:) => null()
  ! real(r8), pointer :: packed_2D(:,:) => null()
  ! real(r8), allocatable :: buffer_2D(:,:)

  !print *, "MGFieldPostProc copy"
  ! lhs%packer = rhs%packer
  ! lhs%field_procs = rhs%field_procs
  lhs%accum_method = rhs%accum_method
  lhs%rank         = rhs%rank
  lhs%num_steps    = rhs%num_steps
  lhs%fillvalue    = rhs%fillvalue
  lhs%unpacked_1D  => rhs%unpacked_1D
  lhs%packed_1D    => rhs%packed_1D
  !if (allocated(rhs%buffer_1D)) print *, "xduan: buffer_1D allocated"
  !lhs%buffer_1D    => rhs%buffer_1D

  lhs%unpacked_2D  => rhs%unpacked_2D
  lhs%packed_2D    => rhs%packed_2D
  !lhs%buffer_2D    => rhs%buffer_2D
  !if (allocated(rhs%buffer_2D)) print *, "xduan: buffer_2D allocated"
end subroutine MGFieldPostProc_copy

end module micro_mg_data
