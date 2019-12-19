module mo_extfrc
  !---------------------------------------------------------------
  ! 	... insitu forcing module
  !---------------------------------------------------------------

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use ppgrid,        only : pcols, begchunk, endchunk, pver, pverp
  use chem_mods,     only : gas_pcnst, extcnt, extfrc_lst, frc_from_dataset
  use spmd_utils,    only : masterproc,iam
  use cam_abortutils,only : endrun
  use cam_history,   only : addfld, outfld, phys_decomp, add_default
  use cam_logfile,   only : iulog
  use tracer_data,   only : trfld,trfile

  implicit none

  type :: forcing
     integer           :: frc_ndx
     real(r8)          :: scalefactor
     character(len=265):: filename
     character(len=16) :: species
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type forcing

  private
  public  :: extfrc_inti
  public  :: extfrc_set
  public  :: extfrc_timestep_init

  save

  integer, parameter :: time_span = 1

  character(len=256) ::   filename

  type(forcing), allocatable  :: forcings(:)
  integer :: n_frc_files = 0

contains

  subroutine extfrc_inti( extfrc_specifier, extfrc_type_in, extfrc_cycle_yr, extfrc_fixed_ymd, extfrc_fixed_tod)

    !-----------------------------------------------------------------------
    ! 	... initialize the surface forcings
    !-----------------------------------------------------------------------
    use cam_pio_utils, only : cam_pio_openfile
    use pio,           only : pio_inq_dimid, pio_inquire, pio_inq_varndims, pio_closefile
    use pio,           only : pio_inq_varname, pio_nowrite, file_desc_t
    use pio,           only : pio_get_att, PIO_NOERR, PIO_GLOBAL
    use pio,           only : pio_seterrorhandling, PIO_BCAST_ERROR,PIO_INTERNAL_ERROR
    use mo_tracname,   only : solsym
    use mo_chem_utls,  only : get_extfrc_ndx, get_spc_ndx
    use chem_mods,     only : frc_from_dataset
    use tracer_data,   only : trcdata_init
    use phys_control,  only : phys_getopts
    use physics_buffer,only : physics_buffer_desc
    use string_utils,  only : GLC

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), dimension(:), intent(in) :: extfrc_specifier
    character(len=*), intent(in) :: extfrc_type_in
    integer  , intent(in)        :: extfrc_cycle_yr
    integer  , intent(in)        :: extfrc_fixed_ymd
    integer  , intent(in)        :: extfrc_fixed_tod

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: astat
    integer :: j, l, m, n, i,mm                          ! Indices
    character(len=16)  :: spc_name
    character(len=256) :: locfn
    character(len=256) :: frc_fnames(gas_pcnst)
    real(r8)           :: frc_scalefactor(gas_pcnst)
    character(len=16)  :: frc_species(gas_pcnst)
    integer            :: frc_indexes(gas_pcnst)

    integer ::  vid, ndims, nvars, isec, ierr
    type(file_desc_t) :: ncid
    character(len=32)  :: varname

    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    logical         , parameter :: rmv_file = .false.
    logical  :: history_aerosol      ! Output the MAM aerosol tendencies

    character(len=32) :: extfrc_type = ' '
    character(len=80) :: file_interp_type = ' '
    character(len=256) :: tmp_string = ' '
    character(len=32) :: xchr = ' '
    real(r8) :: xdbl

    !-----------------------------------------------------------------------
 
    call phys_getopts( history_aerosol_out = history_aerosol   )

    !-----------------------------------------------------------------------
    ! 	... species has insitu forcing ?
    !-----------------------------------------------------------------------

    !write(iulog,*) 'Species with insitu forcings'
    mm = 0

    count_emis: do n=1,gas_pcnst

       if ( len_trim(extfrc_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(extfrc_specifier(n),'->')
       spc_name = trim(adjustl(extfrc_specifier(n)(:i-1)))
       
       ! need to parse out scalefactor ...
       tmp_string = adjustl(extfrc_specifier(n)(i+2:))
       j = scan( tmp_string, '*' )
       if (j>0) then
          xchr = tmp_string(1:j-1) ! get the multipler (left of the '*')
          read( xchr, * ) xdbl   ! convert the string to a real
          tmp_string = adjustl(tmp_string(j+1:)) ! get the filepath name (right of the '*')
       else
          xdbl = 1._r8
       endif
       filename = trim(tmp_string)

       m = get_extfrc_ndx( spc_name )

       if ( m < 1 ) then
          call endrun('extfrc_inti: '//trim(spc_name)// ' does not have an external source')
       endif

       if ( .not. frc_from_dataset(m) ) then
          call endrun('extfrc_inti: '//trim(spc_name)//' cannot have external forcing from additional dataset')
       endif

       mm = mm+1
       frc_species(mm) = spc_name
       frc_fnames(mm) = filename
       frc_indexes(mm) = m
       frc_scalefactor(mm) = xdbl

    enddo count_emis

    n_frc_files = mm

    if( n_frc_files < 1 ) then
       if (masterproc) write(iulog,*) 'There are no species with insitu forcings'
       return
    end if

    if (masterproc) write(iulog,*) ' '

    !-----------------------------------------------------------------------
    ! 	... allocate forcings type array
    !-----------------------------------------------------------------------
    allocate( forcings(n_frc_files), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'extfrc_inti: failed to allocate forcings array; error = ',astat
       call endrun
    end if

    !-----------------------------------------------------------------------
    ! 	... setup the forcing type array
    !-----------------------------------------------------------------------
    do m=1,n_frc_files 
       forcings(m)%frc_ndx     = frc_indexes(m)
       forcings(m)%species     = frc_species(m)
       forcings(m)%filename    = frc_fnames(m)
       forcings(m)%scalefactor = frc_scalefactor(m)
    enddo
    
    do n= 1,extcnt 
       if (frc_from_dataset(n)) then
          spc_name = extfrc_lst(n)
          call addfld( trim(spc_name)//'_XFRC',  'molec/cm3/s', pver, 'A', &
               'external forcing for '//trim(spc_name),   phys_decomp )
          call addfld( trim(spc_name)//'_CLXF',  'molec/cm2/s', 1, 'A', &
               'vertically intergrated external forcing for '//trim(spc_name),   phys_decomp )
          if ( history_aerosol ) then 
             call add_default( trim(spc_name)//'_XFRC', 1, ' ' )
             call add_default( trim(spc_name)//'_CLXF', 1, ' ' )
          endif
       endif
    enddo

    if (masterproc) then
       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'extfrc_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'extfrc timing specs'
       write(iulog,*) 'type = ',extfrc_type
       if( extfrc_type == 'FIXED' ) then
          write(iulog,*) ' fixed date = ', extfrc_fixed_ymd
          write(iulog,*) ' fixed time = ', extfrc_fixed_tod
       else if( extfrc_type == 'CYCLICAL' ) then
          write(iulog,*) ' cycle year = ',extfrc_cycle_yr
       end if
       write(iulog,*) ' '
       write(iulog,*) 'there are ',n_frc_files,' species with external forcing files'
       do m = 1,n_frc_files
          write(iulog,*) ' '
          write(iulog,*) 'forcing type ',m
          write(iulog,*) 'species = ',trim(forcings(m)%species)
          write(iulog,*) 'frc ndx = ',forcings(m)%frc_ndx
          write(iulog,*) 'filename= ',trim(forcings(m)%filename)
       end do
       write(iulog,*) ' '
    endif

    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    frcing_loop: do m = 1, n_frc_files

       forcings(m)%nsectors = 0

       call cam_pio_openfile ( ncid, trim(forcings(m)%filename), PIO_NOWRITE)
       ierr = pio_inquire (ncid, nVariables=nvars)

       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, ndims)

          if( ndims < 4 ) then
             cycle
          elseif( ndims > 4 ) then
             ierr = pio_inq_varname (ncid, vid, varname)
             write(iulog,*) 'extfrc_inti: Skipping variable ', trim(varname),', ndims = ',ndims, &
                  ' , species=',trim(forcings(m)%species)
             cycle
          end if

          forcings(m)%nsectors = forcings(m)%nsectors+1

       enddo

       allocate( forcings(m)%sectors(forcings(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
         write(iulog,*) 'extfrc_inti: failed to allocate forcings(m)%sectors array; error = ',astat
         call endrun
       end if

       isec = 1
       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, ndims)
          if( ndims == 4 ) then
             ierr = pio_inq_varname(ncid, vid, forcings(m)%sectors(isec))
             isec = isec+1
          endif

       enddo

       ! Global attribute 'input_method' overrides the ext_frc_type namelist setting on
       ! a file-by-file basis.  If the ext_frc file does not contain the 'input_method' 
       ! attribute then the ext_frc_type namelist setting is used.
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ierr = pio_get_att(ncid, PIO_GLOBAL, 'input_method', file_interp_type)
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
       if ( ierr == PIO_NOERR) then
          l = GLC(file_interp_type)
          extfrc_type(1:l) = file_interp_type(1:l)
          extfrc_type(l+1:) = ' '
       else
          extfrc_type = trim(extfrc_type_in)
       endif

       call pio_closefile (ncid)

       allocate(forcings(m)%file%in_pbuf(size(forcings(m)%sectors)))
       forcings(m)%file%in_pbuf(:) = .false.
       call trcdata_init( forcings(m)%sectors, &
                          forcings(m)%filename, filelist, datapath, &
                          forcings(m)%fields,  &
                          forcings(m)%file, &
                          rmv_file, extfrc_cycle_yr, extfrc_fixed_ymd, extfrc_fixed_tod, trim(extfrc_type) )

    enddo frcing_loop


  end subroutine extfrc_inti

  subroutine extfrc_timestep_init( pbuf2d, state )
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use tracer_data,  only : advance_trcdata
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    do m = 1,n_frc_files
       call advance_trcdata( forcings(m)%fields, forcings(m)%file, state, pbuf2d  )
    end do

  end subroutine extfrc_timestep_init

  subroutine extfrc_set( lchnk, zint, frcing, ncol )

    !--------------------------------------------------------
    !	... form the external forcing
    !--------------------------------------------------------

    implicit none

    !--------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    :: ncol                  ! columns in chunk
    integer,  intent(in)    :: lchnk                 ! chunk index
    real(r8), intent(in)    :: zint(ncol, pverp)                  ! interface geopot above surface (km)
    real(r8), intent(inout) :: frcing(ncol,pver,extcnt)   ! insitu forcings (molec/cm^3/s)

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n
    character(len=16) :: xfcname
    real(r8) :: frcing_col(1:ncol)
    integer  :: k, isec
    real(r8),parameter :: km_to_cm = 1.e5_r8

    if( n_frc_files < 1 .or. extcnt < 1 ) then
       return
    end if

    frcing(:,:,:) = 0._r8

    !--------------------------------------------------------
    !	... set non-zero forcings
    !--------------------------------------------------------
    file_loop : do m = 1,n_frc_files

       n = forcings(m)%frc_ndx

       do isec = 1,forcings(m)%nsectors
          frcing(:ncol,:,n) = frcing(:ncol,:,n) + forcings(m)%fields(isec)%data(:ncol,:,lchnk)
       enddo

       frcing(:ncol,:,n) = forcings(m)%scalefactor*frcing(:ncol,:,n)

    enddo file_loop

    frc_loop : do n = 1,extcnt
       if (frc_from_dataset(n)) then

          xfcname = trim(extfrc_lst(n))//'_XFRC'
          call outfld( xfcname, frcing(:ncol,:,n), ncol, lchnk )

          frcing_col(:ncol) = 0._r8
          do k = 1,pver
             frcing_col(:ncol) = frcing_col(:ncol) + frcing(:ncol,k,n)*(zint(:ncol,k)-zint(:ncol,k+1))*km_to_cm
          enddo

          xfcname = trim(extfrc_lst(n))//'_CLXF'
          call outfld( xfcname, frcing_col(:ncol), ncol, lchnk )
       endif
    end do frc_loop

  end subroutine extfrc_set


end module mo_extfrc
