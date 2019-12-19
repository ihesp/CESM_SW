subroutine qneg3_f(idx, ncol, ncold, lver, lconst_beg, lconst_end, qmin, q, nvals, iw, kw, found, worst)
!-----------------------------------------------------------------------
!
! Purpose:
! Check moisture and tracers for minimum value, reset any below
! minimum value to minimum value and return information to allow
! warning message to be printed. The global average is NOT preserved.
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: J. Rosinski
!
!-----------------------------------------------------------------------
  implicit none
    integer, parameter :: R8 = selected_real_kind(12)

!------------------------------Arguments--------------------------------
!
! Input arguments
!
!   character*(*), intent(in) :: subnam ! name of calling routine

    integer, intent(in) :: idx          ! chunk/latitude index
    integer, intent(in) :: ncol         ! number of atmospheric columns
    integer, intent(in) :: ncold        ! declared number of atmospheric columns
    integer, intent(in) :: lver         ! number of vertical levels in column
    integer, intent(in) :: lconst_beg   ! beginning constituent
    integer, intent(in) :: lconst_end   ! ending    constituen
    integer, intent(inout) :: nvals(lconst_beg:lconst_end)      ! Global minimum constituent concentration
    integer, intent(inout) :: iw(lconst_beg:lconst_end)      ! Global minimum constituent concentration
    integer, intent(inout) :: kw(lconst_beg:lconst_end)      ! Global minimum constituent concentration
    integer, intent(inout) :: found(lconst_beg:lconst_end)      ! Global minimum constituent concentration
    real(r8), intent(inout) :: worst(lconst_beg:lconst_end) ! moisture/tracer field
    real(r8), intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

! Input/Output arguments
    real(r8), intent(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
!---------------------------Local workspace-----------------------------
    integer indx(ncol,lver)  ! array of indices of points < qmin
    integer nval(lver)       ! number of points < qmin for 1 level
    !integer nvals            ! number of values found < qmin
    integer nn
    integer iwtmp
    integer i,ii,k           ! longitude, level indices
    integer m                ! constituent index

!-----------------------------------------------------------------------
    !print *, iw(m), kw(m), found(m), worst(m)
    do m=lconst_beg,lconst_end
      nvals(m) = 0
      found(m) = -1
      worst(m) = 1.e35_r8
      iw(m) = -1

      do k=1,lver
        nval(k) = 0
        nn = 0
        do i=1,ncol
          if (q(i,k,m) < qmin(m)) then
             nn = nn + 1
             indx(nn,k) = i
          end if
        end do
        nval(k) = nn
      end do

      do k=1,lver
        if (nval(k) > 0) then
          found(m) = 1
          nvals(m) = nvals(m) + nval(k)
          iwtmp = -1
!cdir noe p,altcode=loopcnt
          do ii=1,nval(k)
             i = indx(ii,k)
             if (q(i,k,m) < worst(m)) then
                worst(m) = q(i,k,m)
                iwtmp = ii
             end if
          end do
          if (iwtmp /= -1 ) kw(m) = k
          if (iwtmp /= -1 ) iw(m) = indx(iwtmp,k)
!cdir noep,altcode=loopcnt
          do ii=1,nval(k)
           i = indx(ii,k)
           q(i,k,m) = qmin(m)
          end do
        end if
      end do
    end do
    return
end subroutine qneg3_f
