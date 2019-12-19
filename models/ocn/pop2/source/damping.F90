!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module damping

!BOP
! !MODULE: damping
!
! !DESCRIPTION:
!  This is the main driver module for UVEL and VVEL damping.
!
! !REVISION HISTORY:
!  SVN:$Id: damping.F90 65485 2014-11-16 22:26:36Z mlevy@ucar.edu $
!
! !USES:

   use kinds_mod,       only : r8, log_kind, int_kind
   use constants,       only : c1
   use time_management, only : steps_per_year
   use domain_size,     only : km
   use blocks,          only : nx_block, ny_block
   use communicate,     only : my_task, master_task
   use broadcast,       only : broadcast_scalar
   use exit_mod,        only : exit_POP, sigAbort
   use io_types,        only : nml_in, nml_filename
#ifdef CCSMCOUPLED
   use shr_sys_mod, only: shr_sys_abort
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_damping,                     &
             damping_uv

  logical(log_kind), public :: ldamp_uv        ! flag to damp UVEL and VVEL

!EOP
!***********************************************************************

 contains

!***********************************************************************

!BOP
! !IROUTINE: init_damping
! !INTERFACE:

  subroutine init_damping()

! !DESCRIPTION:
!  Initializes damping by reading in damping_nml namelist and setting
!  module-wide variables
!
! !REVISION HISTORY:
!  same as module
    integer(int_kind) :: nml_error

    namelist /damping_nml/ldamp_uv

    ldamp_uv      = .false.

    if (my_task.eq.master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
         read(nml_in, nml=damping_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
      print*, "Done setting default values!"
    end if

!!!! BROADCAST NAMELIST OPTIONS OUT
    call broadcast_scalar(nml_error, master_task)
    if (nml_error.ne.0) then
      call exit_POP(sigAbort,'ERROR: reading damping_nml')
    end if

    call broadcast_scalar(ldamp_uv, master_task)
!EOC

  end subroutine init_damping

!***********************************************************************

!BOP
! !IROUTINE: damping_uv
! !INTERFACE:

  subroutine damping_uv(UVEL, VVEL)

! !DESCRIPTION:
!  Apply damping to UVEL and VVEL fields
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:
    real(r8), dimension(nx_block, ny_block, km), intent(inout) :: UVEL, VVEL

! !LOCAL VARIABLES:
    real(r8), dimension(nx_block, ny_block, km) :: URESTORE, VRESTORE

!EOP

!BOC
    URESTORE = min(0.99_r8, abs(UVEL)/steps_per_year)
    VRESTORE = min(0.99_r8, abs(VVEL)/steps_per_year)

    UVEL = UVEL * (c1 - URESTORE)
    VVEL = VVEL * (c1 - VRESTORE)

!EOC

  end subroutine damping_uv

end module damping
