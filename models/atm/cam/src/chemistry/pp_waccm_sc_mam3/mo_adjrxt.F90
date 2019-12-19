




      module mo_adjrxt

      private
      public :: adjrxt

      contains

      subroutine adjrxt( rate, inv, m, ncol )

      use ppgrid, only : pver
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : nfs, rxntot

      implicit none

!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in) :: inv(ncol,pver,nfs)
      real(r8), intent(in) :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      real(r8) :: im(ncol,pver)


         rate(:,:, 8) = rate(:,:, 8) * inv(:,:, 5)
         rate(:,:, 9) = rate(:,:, 9) * inv(:,:, 5)
         rate(:,:, 10) = rate(:,:, 10) * inv(:,:, 5)
         rate(:,:, 11) = rate(:,:, 11) * inv(:,:, 5)
         rate(:,:, 12) = rate(:,:, 12) * inv(:,:, 6)
         im(:,:) = 1._r8 / m(:,:)
         rate(:,:, 7) = rate(:,:, 7) * inv(:,:, 7) * inv(:,:, 7) * im(:,:)

      end subroutine adjrxt

      end module mo_adjrxt
