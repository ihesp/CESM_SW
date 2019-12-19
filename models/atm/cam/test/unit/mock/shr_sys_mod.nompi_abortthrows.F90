module shr_sys_mod

! This is a mock version of shr_sys_mod.
! It contains only an abort method that throws a pFUnit exception instead
! of actually aborting.

use shr_kind_mod, only: &
     shr_kind_in

implicit none
private
save

! Fake abort
public :: shr_sys_abort

contains

subroutine shr_sys_abort(string, rc)
  use pfunit_mod, only: throw

  character(*), optional :: string
  integer(shr_kind_in), optional :: rc

  ! Prevent compiler spam about unused variables.
  if (.false.) rc = rc

  call throw("ABORTED: "//trim(string))

end subroutine shr_sys_abort

end module shr_sys_mod
