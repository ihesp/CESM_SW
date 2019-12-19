!---------------------------------------------------------------------------------
! Manages the CFC11* for radiation 
!  4 Dec 2009 -- Francis Vitt created
!  8 Mar 2013 -- expanded for waccm_tsmlt -- fvitt
!---------------------------------------------------------------------------------
module cfc11star

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  
  use physics_buffer, only : pbuf_add_field, dtype_r8
  use cam_abortutils, only : endrun
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use spmd_utils,   only : masterproc
  use constituents, only : cnst_get_ind

  implicit none
  save 

  private
  public :: register_cfc11star
  public :: update_cfc11star
  public :: init_cfc11star

  logical :: do_cfc11star
  character(len=16), parameter :: pbufname = 'CFC11STAR'
  integer :: pbf_idx = -1
  integer, parameter :: ncfcs = 13

  integer, target :: indices(ncfcs)
  
  real(r8) :: rel_rf(ncfcs)

contains

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine register_cfc11star

    implicit none

    integer :: m

    character(len=8), parameter :: species(ncfcs) = &
      (/ 'CFC11   ','CFC113  ','CFC114  ','CFC115  ','CCL4    ','CH3CCL3 ','CH3CL   ','HCFC22  ',&
         'HCFC141B','HCFC142B','CF2CLBR ','CF3BR   ','H2402   ' /)
    real(r8), parameter :: cfc_rf(ncfcs) = &
      (/  0.25_r8,   0.30_r8,   0.31_r8,   0.18_r8,   0.13_r8,   0.06_r8,   0.01_r8,   0.20_r8,  &
          0.14_r8,   0.20_r8,   0.30_r8,   0.32_r8,   0.33_r8 /) ! W/m2/ppb

    do m = 1, ncfcs 
       call cnst_get_ind(species(m), indices(m), abort=.false.)
    enddo

    do_cfc11star = any(indices(:)>0)
    if (.not.do_cfc11star) return

    call pbuf_add_field(pbufname,'global',dtype_r8,(/pcols,pver/),pbf_idx)

    rel_rf(:) = cfc_rf(:) / cfc_rf(1)

  endsubroutine register_cfc11star

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine init_cfc11star(pbuf2d)
    use cam_history,  only : addfld, phys_decomp
    use infnan,       only : nan, assignment(=)
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field

    implicit none

    real(r8) :: real_nan
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (.not.do_cfc11star) return

    real_nan = nan
    call pbuf_set_field(pbuf2d, pbf_idx, real_nan)

    call addfld(pbufname,'kg/kg',pver,'A','cfc11star for radiation', phys_decomp )
    
    if (masterproc) then
       write(iulog,*) 'init_cfc11star: CFC11STAR is added to pbuf2d for radiation'
    endif
  end subroutine init_cfc11star

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine update_cfc11star( pbuf2d, phys_state )
    use cam_history,  only : outfld
    use physics_types,only : physics_state
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    integer :: lchnk, ncol
    integer :: c, m
    real(r8), pointer :: cf11star(:,:)

    if (.not.do_cfc11star) return
    
    do c = begchunk,endchunk
       lchnk = phys_state(c)%lchnk
       ncol = phys_state(c)%ncol

       call pbuf_get_field(pbuf_get_chunk(pbuf2d, lchnk), pbf_idx, cf11star)

       cf11star(:ncol,:) = 0._r8
       do m = 1, ncfcs 
          if ( indices(m)>0 ) then
             cf11star(:ncol,:) = cf11star(:ncol,:) &
                               + phys_state(c)%q(:ncol,:,indices(m)) * rel_rf(m) 
          endif
       enddo

       call outfld( pbufname, cf11star(:ncol,:), ncol, lchnk) 

    enddo

  endsubroutine update_cfc11star

end module cfc11star
