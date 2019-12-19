!--------------------------------------------------------------------------------
! Manages writing reaction rates to history
!--------------------------------------------------------------------------------
module rate_diags

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use shr_kind_mod,     only : CL => SHR_KIND_CL, CX => SHR_KIND_CX
  use cam_history,      only : fieldname_len
  use cam_history,      only : addfld,phys_decomp
  use cam_history,      only : outfld
  use chem_mods,        only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
  use ppgrid,           only : pver
  use spmd_utils,       only : masterproc
  use cam_abortutils,   only : endrun

  implicit none
  private 
  public :: rate_diags_init
  public :: rate_diags_calc
  public :: rate_diags_readnl

  character(len=fieldname_len) :: rate_names(rxt_tag_cnt)

  type rate_grp_t
    character(len=24) :: name
    integer :: nmembers = 0
    integer, allocatable :: map(:)
    real(r8), allocatable :: multipler(:)
  endtype rate_grp_t

  integer :: ngrps = 0
  type(rate_grp_t), allocatable :: grps(:)  

  integer, parameter :: maxsums = 100
  character(len=CX) :: rxn_rate_sums(maxsums) = ' '

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine rate_diags_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    ! args 
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    namelist /rxn_rate_diags_nl/ rxn_rate_sums

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'rxn_rate_diags_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, rxn_rate_diags_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('rate_diags_readnl:: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(rxn_rate_sums,len(rxn_rate_sums(1))*maxsums, mpichar, 0, mpicom)
#endif
  end subroutine rate_diags_readnl
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_init

    integer :: i, len, pos

    character(len=64) :: name

    do i = 1,rxt_tag_cnt
       pos = 0
       pos = index(rxt_tag_lst(i),'tag_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'usr_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'cph_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'ion_')
       if (pos>0) then
          name = 'r_'//trim(rxt_tag_lst(i)(5:))
       else
          name = 'r_'//trim(rxt_tag_lst(i)(1:))
       endif
       len = min(fieldname_len,len_trim(name))
       rate_names(i) = trim(name(1:len))
       call addfld(rate_names(i), 'molecules/cm3/sec', pver,'A','reaction rate', phys_decomp)
    enddo

    call parse_rate_sums()

    do i = 1, ngrps
       call addfld( grps(i)%name, 'molecules/cm3/sec', pver,'A','reaction rate group', phys_decomp)
    enddo

  end subroutine rate_diags_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_calc( rxt_rates, vmr, m, ncol, lchnk )

    use mo_rxt_rates_conv, only: set_rates

    real(r8), intent(inout) :: rxt_rates(:,:,:) ! 'molec/cm3/sec'
    real(r8), intent(in)    :: vmr(:,:,:)
    real(r8), intent(in)    :: m(:,:)           ! air density (molecules/cm3)
    integer,  intent(in)    :: ncol, lchnk

    integer :: i, j
    real(r8) :: group_rate(ncol,pver)

    call set_rates( rxt_rates, vmr, ncol )

    ! output individual tagged rates    
    do i = 1, rxt_tag_cnt
       ! convert from vmr/sec to molecules/cm3/sec
       rxt_rates(:ncol,:,rxt_tag_map(i)) = rxt_rates(:ncol,:,rxt_tag_map(i)) * m(:ncol,:)
       call outfld( rate_names(i), rxt_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    ! output rate groups ( or families )
    do i = 1, ngrps
       group_rate(:,:) = 0._r8
       do j = 1, grps(i)%nmembers
         group_rate(:ncol,:) = group_rate(:ncol,:) + grps(i)%multipler(j)*rxt_rates(:ncol,:,grps(i)%map(j))
       enddo 
       call outfld( grps(i)%name, group_rate(:ncol,:), ncol, lchnk )       
    end do

  end subroutine rate_diags_calc

!-------------------------------------------------------------------
! Private routines :
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  
  subroutine parse_rate_sums

    integer :: ndxs(512)
    integer :: nelem, spc_len, i,j,k, rxt_ndx
    character(len=CL) :: tmp_str, tmp_name

    character(len=8) :: xchr ! multipler
    real(r8) :: xdbl

    character(len=CX) :: sum_string

    ! a group is  a sum of reaction rates 

    ! count the numger of sums (or groups)
    sumcnt: do i = 1,maxsums
       spc_len=len_trim(rxn_rate_sums(i))
       if ( spc_len > 0 ) then
          ngrps = ngrps+1
       else
          exit sumcnt
       endif
    enddo sumcnt

    ! parse the individual sum strings...  and form the groupings
    has_grps: if (ngrps>0) then

       allocate( grps(ngrps) )

       ! from shr_megan_mod ... should be generalized and shared...
       grploop: do i = 1,ngrps

          ! parse out the rxn names and multipliers
          ! from first parsing out the terms in the summation equation ("+" separates the terms)

          sum_string = rxn_rate_sums(i)
          j = scan( sum_string, '=' )
          nelem = 1
          ndxs(nelem) = j ! ndxs stores the index of each term of the equation

          ! find indices of all the terms in the equation
          tmp_str = trim( sum_string(j+1:) )
          j = scan( tmp_str, '+' )
          do while(j>0)
             nelem = nelem+1
             ndxs(nelem) = ndxs(nelem-1) + j
             tmp_str = tmp_str(j+1:)
             j = scan( tmp_str, '+' )
          enddo
          ndxs(nelem+1) = len(sum_string)+1

          grps(i)%nmembers = nelem ! number of terms 
          grps(i)%name =  trim(adjustl( sum_string(:ndxs(1)-1))) ! thing to the left of the "=" is used as the name of the group

          ! now that we have the number of terms in the summation allocate memory for the map (reaction indices) and multipliers
          allocate(grps(i)%map(nelem)) 
          allocate(grps(i)%multipler(nelem))

          ! now parse out the  rxn names and multiplers from the terms 
          elmloop: do k = 1,nelem
             grps(i)%multipler(k) = 1._r8
             ! get the rxn name which follows the '*' operator if the is one
             tmp_name = adjustl(sum_string(ndxs(k)+1:ndxs(k+1)-1))
             j = scan( tmp_name, '*' )
             if (j>0) then
                xchr = tmp_name(1:j-1) ! get the multipler (left of the '*')
                read( xchr, * ) xdbl   ! convert the string to a real
                grps(i)%multipler(k) = xdbl ! store the multiplier
                tmp_name = adjustl(tmp_name(j+1:)) ! get the rxn name (right of the '*')
             endif
             ! look up the corresponding reaction index ...
             rxt_ndx = lookup_tag_ndx( tmp_name )
             if ( rxt_ndx > 0 ) then
                grps(i)%map(k) = rxt_ndx
             else
                call endrun('rate_diags::parse_rate_sums rate name not found : '//trim(tmp_name))
             endif
          enddo elmloop
       enddo grploop
    endif has_grps

  end subroutine parse_rate_sums

!-------------------------------------------------------------------
! finds the index corresponging to a given reacton name
!-------------------------------------------------------------------
  function lookup_tag_ndx( name ) result( ndx )
    character(len=*) :: name
    integer :: ndx

    integer :: i

    ndx = -1

    findloop: do i = 1,rxt_tag_cnt
       if (trim(name) .eq. trim(rate_names(i)(3:))) then
          ndx = i
          return
       endif
    end do findloop

  end function lookup_tag_ndx

end module rate_diags
