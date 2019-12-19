




      module mo_phtadj

      private
      public :: phtadj

      contains

      subroutine phtadj( p_rate, inv, m, ncol )

      use chem_mods, only : nfs, phtcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in) :: inv(:,:,:)
      real(r8), intent(in) :: m(:,:)
      real(r8), intent(inout) :: p_rate(:,:,:)

!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      integer :: k
      real(r8) :: im(ncol)

      do k = 1,pver
         im(:ncol) = 1._r8 / m(:ncol,k)
         p_rate(:,k, 97) = p_rate(:,k, 97) * inv(:,k, 2) * im(:)
         p_rate(:,k,101) = p_rate(:,k,101) * inv(:,k, 2) * im(:)
         p_rate(:,k,102) = p_rate(:,k,102) * inv(:,k, 2) * im(:)
         p_rate(:,k,104) = p_rate(:,k,104) * inv(:,k, 2) * im(:)
         p_rate(:,k,109) = p_rate(:,k,109) * inv(:,k, 2) * im(:)
         p_rate(:,k,113) = p_rate(:,k,113) * inv(:,k, 2) * im(:)
         p_rate(:,k,114) = p_rate(:,k,114) * inv(:,k, 2) * im(:)
         p_rate(:,k,116) = p_rate(:,k,116) * inv(:,k, 2) * im(:)
      end do

      end subroutine phtadj

      end module mo_phtadj
