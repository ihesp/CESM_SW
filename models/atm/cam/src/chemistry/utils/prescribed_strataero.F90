!-------------------------------------------------------------------
! manages reading and interpolation of prescribed stratospheric aerosols
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_strataero

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use tracer_data,      only : trfld, trfile
  use cam_logfile,      only : iulog

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_strataero_readnl
  public :: prescribed_strataero_register
  public :: prescribed_strataero_init
  public :: prescribed_strataero_adv
  public :: write_prescribed_strataero_restart
  public :: read_prescribed_strataero_restart
  public :: has_prescribed_strataero
  public :: init_prescribed_strataero_restart

  logical :: has_prescribed_strataero = .false.
  character(len=16), parameter :: mmr_name = 'VOLC_MMR'
  character(len=16), parameter :: rad_name = 'VOLC_RAD_GEOM'
  character(len=16), parameter :: sad_name = 'VOLC_SAD'
  character(len=16), parameter :: mass_name = 'VOLC_MASS'
  character(len=16), parameter :: mass_column_name = 'VOLC_MASS_C'
  character(len=16), parameter :: dens_name = 'VOLC_DENS'

  ! These variables are settable via the namelist (with longer names)
  character(len=32)  :: specifier(3) = (/'VOLC_MMR:H2SO4_mass             ', &
                                         'VOLC_RAD_GEOM:rmode             ', &
                                         'VOLC_SAD:sad                    ' /)
  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  integer            :: rad_ndx = -1
  integer            :: sad_ndx = -1
  integer            :: mmr_ndx = -1


contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_strataero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_strataero_readnl'

   character(len=32)  :: prescribed_strataero_specifier(3)
   character(len=256) :: prescribed_strataero_file
   character(len=256) :: prescribed_strataero_filelist
   character(len=256) :: prescribed_strataero_datapath
   character(len=32)  :: prescribed_strataero_type
   logical            :: prescribed_strataero_rmfile
   integer            :: prescribed_strataero_cycle_yr
   integer            :: prescribed_strataero_fixed_ymd
   integer            :: prescribed_strataero_fixed_tod

   namelist /prescribed_strataero_nl/ &
      prescribed_strataero_specifier, &
      prescribed_strataero_file,      &
      prescribed_strataero_filelist,  &
      prescribed_strataero_datapath,  &
      prescribed_strataero_type,      &
      prescribed_strataero_rmfile,    &
      prescribed_strataero_cycle_yr,  &
      prescribed_strataero_fixed_ymd, &
      prescribed_strataero_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_strataero_specifier= specifier
   prescribed_strataero_file     = filename
   prescribed_strataero_filelist = filelist
   prescribed_strataero_datapath = datapath
   prescribed_strataero_type     = data_type
   prescribed_strataero_rmfile   = rmv_file
   prescribed_strataero_cycle_yr = cycle_yr
   prescribed_strataero_fixed_ymd= fixed_ymd
   prescribed_strataero_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_strataero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_strataero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_strataero_specifier,len(prescribed_strataero_specifier)*2, mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_file,     len(prescribed_strataero_file),        mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_filelist, len(prescribed_strataero_filelist),    mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_datapath, len(prescribed_strataero_datapath),    mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_type,     len(prescribed_strataero_type),        mpichar, 0, mpicom)
   call mpibcast(prescribed_strataero_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_strataero_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_strataero_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_strataero_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier(:) = prescribed_strataero_specifier(:)
   filename   = prescribed_strataero_file
   filelist   = prescribed_strataero_filelist
   datapath   = prescribed_strataero_datapath
   data_type  = prescribed_strataero_type
   rmv_file   = prescribed_strataero_rmfile
   cycle_yr   = prescribed_strataero_cycle_yr
   fixed_ymd  = prescribed_strataero_fixed_ymd
   fixed_tod  = prescribed_strataero_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_strataero = .true.

end subroutine prescribed_strataero_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_strataero_register()
    use ppgrid,         only: pver,pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer :: idx

    if (has_prescribed_strataero) then
       call pbuf_add_field(mmr_name, 'physpkg', dtype_r8,(/pcols,pver/), mmr_ndx)
       call pbuf_add_field(rad_name, 'physpkg', dtype_r8,(/pcols,pver/), rad_ndx)
       call pbuf_add_field(sad_name, 'physpkg', dtype_r8,(/pcols,pver/), sad_ndx)
    endif

  endsubroutine prescribed_strataero_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_strataero_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    
    use physics_buffer, only: physics_buffer_desc, pbuf_get_index

    implicit none

    integer :: ndx, istat
    integer :: errcode
    
    if ( has_prescribed_strataero ) then
       if ( masterproc ) then
          write(iulog,*) 'stratospheric aerosol is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .true.
    file%geop_alt = .true.

    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

    call addfld(dens_name,'molecules/cm3', pver, 'I', 'prescribed volcanic aerosol number density', phys_decomp )
    call addfld(mmr_name,'kg/kg', pver, 'I', 'prescribed volcanic aerosol dry mass mixing ratio', phys_decomp )
    call addfld(rad_name,'m', pver, 'I', 'volcanic aerosol geometric-mean radius', phys_decomp )
    call addfld(mass_name,'kg/m^2', pver, 'I', 'volcanic aerosol vertical mass path in layer', phys_decomp )
    call addfld(mass_column_name,'kg/m^2', 1, 'I', 'volcanic aerosol column mass', phys_decomp )
    call addfld(sad_name,'cm2/cm3', pver, 'I', 'stratospheric aerosol density', phys_decomp )

  end subroutine prescribed_strataero_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_strataero_adv( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz, gravit        ! J/K/molecule
    use tropopause,   only : tropopause_find, TROP_ALG_TWMO, TROP_ALG_CLIMATE
    
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    
    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer :: c,ncol,i,k
    real(r8) :: to_mmr(pcols,pver)
    real(r8), parameter :: molmass = 4._r8/3._r8*98.0_r8 !convert dry mass to wet mass of h2so4 
    real(r8) :: ptrop
    real(r8) :: concvolc ! micrograms of wetted aerosol per cubic centimeter
    real(r8) :: volcmass(pcols,pver)
    real(r8) :: columnmass(pcols)
    real(r8) :: mmrvolc
    integer  :: tropLev(pcols)
    real(r8) :: area_fact

    real(r8) :: outdata(pcols,pver)
    real(r8), pointer :: mass(:,:)
    real(r8), pointer :: area(:,:)
    real(r8), pointer :: radius(:,:)

    !WACCM-derived relation between mass concentration and wet aerosol radius in meters
    real(r8),parameter :: radius_conversion = 1.9e-4_r8

    if( .not. has_prescribed_strataero ) return

    call advance_trcdata( fields, file, state, pbuf2d )

    ! copy prescribed tracer fields into state svariable with the correct units
    do c = begchunk,endchunk

       pbuf_chnk => pbuf_get_chunk(pbuf2d, c)

       ncol = state(c)%ncol

       select case ( to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))) )
       case ("molecules/cm3air", "molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
          to_mmr(:ncol,:) = (molmass*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
       case ('kg/kg','mmr','kg kg-1')
          to_mmr(:ncol,:) = 1._r8
       case ('mol/mol','mole/mole','vmr','fraction')
          to_mmr(:ncol,:) = molmass/mwdry
       case default
          write(iulog,*) 'prescribed_strataero_adv: mass units = ',trim(fields(1)%units) ,' are not recognized'
          call endrun('prescribed_strataero_adv: mass units are not recognized')
       end select

       call pbuf_get_field(pbuf_chnk, mmr_ndx, mass)

       call outfld( dens_name,         mass(:,:),     pcols, state(c)%lchnk)

       mass(:ncol,:) = to_mmr(:ncol,:) * mass(:ncol,:) ! mmr

       call pbuf_get_field(pbuf_chnk, rad_ndx, radius)
       call pbuf_get_field(pbuf_chnk, sad_ndx, area)

       select case ( to_lower(trim(fields(3)%units(:7))) )
       case ("um2/cm3")
          area_fact = 1.e-8_r8
       case default
          write(iulog,*) 'prescribed_strataero_adv: surface area density units = ',trim(fields(3)%units) ,' are not recognized'
          call endrun('prescribed_strataero_adv: surface area density are not recognized')
       end select
       area(:ncol,:) = area_fact*area(:ncol,:)

       ! this definition of tropopause is consistent with what is used in chemistry
       call tropopause_find(state(c), tropLev, primary=TROP_ALG_TWMO, backup=TROP_ALG_CLIMATE)

       do i = 1,ncol
          do k = 1,pver
             ! set to zero at and below tropopause
             if ( k >= tropLev(i) ) then
                mass(i,k) = 0._r8
                radius(i,k) = 0._r8
                area(i,k) = 0._r8
             endif
          enddo
       enddo

       volcmass(:ncol,:) = mass(:ncol,:)*state(c)%pdel(:ncol,:)/gravit
       columnmass(:ncol) = sum(volcmass(:ncol,:), 2)

       call outfld( mmr_name,         mass(:,:),     pcols, state(c)%lchnk)
       call outfld( mass_name,        volcmass(:,:), pcols, state(c)%lchnk)
       call outfld( mass_column_name, columnmass(:), pcols, state(c)%lchnk)
       call outfld( rad_name,         radius(:,:),   pcols, state(c)%lchnk)
       call outfld( sad_name,         area(:,:),     pcols, state(c)%lchnk)

    enddo

  end subroutine prescribed_strataero_adv

!-------------------------------------------------------------------
  subroutine init_prescribed_strataero_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_strataero', piofile, file )

  end subroutine init_prescribed_strataero_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_strataero_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_strataero_restart
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_strataero_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_strataero', piofile, file )

  end subroutine read_prescribed_strataero_restart

end module prescribed_strataero
