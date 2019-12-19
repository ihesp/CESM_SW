




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
         p_rate(:,k, 94) = p_rate(:,k, 94) * inv(:,k, 2) * im(:)
         p_rate(:,k, 98) = p_rate(:,k, 98) * inv(:,k, 2) * im(:)
         p_rate(:,k, 99) = p_rate(:,k, 99) * inv(:,k, 2) * im(:)
         p_rate(:,k,101) = p_rate(:,k,101) * inv(:,k, 2) * im(:)
         p_rate(:,k,106) = p_rate(:,k,106) * inv(:,k, 2) * im(:)
         p_rate(:,k,110) = p_rate(:,k,110) * inv(:,k, 2) * im(:)
         p_rate(:,k,111) = p_rate(:,k,111) * inv(:,k, 2) * im(:)
         p_rate(:,k,113) = p_rate(:,k,113) * inv(:,k, 2) * im(:)
      end do

      end subroutine phtadj

      end module mo_phtadj
