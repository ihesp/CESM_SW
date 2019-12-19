!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module passive_tracer_tools

!BOP
! !MODULE: passive_tracer_tools

! !DESCRIPTION:
!  This module provides support for passive tracers.
!  Subroutines in this module can be called by individual tracer modules.

! !REVISION HISTORY:
!  SVN:$Id:$

! !USES:

   use kinds_mod
   use constants, only: c0, c1, p5, char_blank, blank_fmt
   use domain_size, only: nx_global, ny_global, max_blocks_clinic,nt
   use domain, only: nblocks_clinic,distrb_clinic, blocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use constants, only: char_blank, field_loc_center, field_type_scalar,delim_fmt
   use prognostic, only: tracer_field
   use io_tools, only: document
   use io, only: data_set
   use io_types, only: datafile, io_dim, io_field_desc, rec_type_dbl, &
       construct_file, construct_io_dim, construct_io_field, &
       destroy_file, destroy_io_field, stdout,add_attrib_file,&
       extract_attrib_file
   use prognostic, only: curtime, oldtime
   use grid, only: TAREA, RCALCT, area_t
   use global_reductions, only: global_sum
   use blocks

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: ind_name_pair
      integer(int_kind) :: ind
      character(char_len) :: name
   end type ind_name_pair

!-----------------------------------------------------------------------
!  derived type for reading tracers from a file
!-----------------------------------------------------------------------

   type, public :: tracer_read
      character(char_len) :: mod_varname, filename, file_varname, file_fmt
      real(r8) :: scale_factor, default_val

   end type tracer_read

!-----------------------------------------------------------------------
!  monthly forcing variables
!-----------------------------------------------------------------------

   type, public :: forcing_monthly_every_ts
      type(tracer_read) :: input
      logical(log_kind) :: has_data
      real(r8), dimension(:,:,:,:,:), pointer :: DATA
      character(char_len) :: &
         interp_type,          & ! = 'linear'
         data_type,            & ! = 'monthly-calendar'
         interp_freq,          & ! = 'every-timestep'
         filename,             & ! = 'not-used-for-monthly'
         data_label              ! = 'not-used-for-monthly'
    real(r8), dimension(12) :: &
         data_time               ! times where DATA is given
    real(r8), dimension(20) :: &
         data_renorm             ! not used for monthly
    real(r8) :: &
         data_inc,             & ! not used for monthly data
         data_next,            & ! time that will be used for the next
                                 ! value of forcing data that is needed
         data_update,          & ! time when the a new forcing value
                                 ! needs to be added to interpolation set
         interp_inc,           & ! not used for 'every-timestep' interp
         interp_next,          & ! not used for 'every-timestep' interp
         interp_last             ! not used for 'every-timestep' interp
    integer(int_kind) :: &
         data_time_min_loc       ! index of the third dimension of data_time
                                 ! containing the minimum forcing time
  end type forcing_monthly_every_ts

! !PUBLIC MEMBER FUNCTIONS:

   public ::                                 &
      init_forcing_monthly_every_ts,         &
      field_exists_in_file,                  &
      file_read_tracer_block,                &
      rest_read_tracer_block,                &
      read_field,                            &
      name_to_ind,                           &
      extract_surf_avg,                      &
      comp_surf_avg,                         &
      set_tracer_indices,                    &
      tracer_read_init

!EOP
!BOC

   interface read_field
      module procedure read_field_2D,        &
                       read_field_3D
   end interface

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_forcing_monthly_every_ts(var)
! !INTERFACE:

 subroutine init_forcing_monthly_every_ts(var)

! !DESCRIPTION:
!  initialize fields of a variable of type forcing_monthly_every_ts
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   type(forcing_monthly_every_ts), intent(out) :: var

!EOP
!BOC
!-----------------------------------------------------------------------

   var%interp_type = 'linear'
   var%data_type   = 'monthly-calendar'
   var%interp_freq = 'every-timestep'
   var%filename    = 'not-used-for-monthly'
   var%data_label  = 'not-used-for-monthly'

!-----------------------------------------------------------------------
!EOC

 end subroutine init_forcing_monthly_every_ts

!***********************************************************************
!BOP
! !IROUTINE: name_to_ind
! !INTERFACE:

 function name_to_ind(name, ind_name_table)

! !DESCRIPTION:
!  translate string into a tracer index
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in) :: name

   type(ind_name_pair), dimension(:), intent(in) :: ind_name_table

! !OUTPUT PARAMETERS:

   integer(int_kind) :: name_to_ind

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer(int_kind) :: i

!-----------------------------------------------------------------------

   do i = 1,size(ind_name_table)
      if (trim(name) == trim(ind_name_table(i)%name)) then
         name_to_ind = ind_name_table(i)%ind
         return
      end if
   end do

   name_to_ind = 0

!-----------------------------------------------------------------------
!EOC

 end function name_to_ind

!***********************************************************************
!BOP
! !IROUTINE: field_exists_in_file
! !INTERFACE:

 function field_exists_in_file(fmt, filename, fieldname)

! !DESCRIPTION:
!  This routine determines if a field exists in a file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      fmt,                 & ! format (bin or nc)
      filename,            & ! file to read from
      fieldname              ! field to be read

! !OUTPUT PARAMETERS:

   logical (log_kind) :: field_exists_in_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'passive_tracer_tools:field_exists_in_file'

   type (datafile) :: &
      in_file         ! data file type for init ts file

!-----------------------------------------------------------------------

   in_file = construct_file(fmt, full_name=trim(filename), &
                            record_length=rec_type_dbl, &
                            recl_words=nx_global*ny_global)

   call data_set(in_file, 'open_read')

   call data_set (in_file, 'field_exists', &
                  fieldname=fieldname, field_exists=field_exists_in_file)

   call data_set (in_file, 'close')

   call destroy_file (in_file)

!-----------------------------------------------------------------------
!EOC

 end function field_exists_in_file

!***********************************************************************
!BOP
! !IROUTINE: rest_read_tracer_block
! !INTERFACE:

 subroutine rest_read_tracer_block(fmt, filename, tracer_d_module, &
    TRACER_MODULE)

! !DESCRIPTION:
!  read from a restart file all tracers for a tracer module
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      fmt,                 & ! format (bin or nc)
      filename               ! file name for restart file

   type (tracer_field), dimension(:), intent(in) :: &
      tracer_d_module   ! descriptors for each tracer

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(:,:,:,:,:,:), intent(inout) :: &
      TRACER_MODULE     ! tracers to be read in

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'passive_tracer_tools:rest_read_tracer_block'

   integer (int_kind) :: &
      n               ! tracer index

   character (char_len) ::  &
      fieldname       ! tracer name temporaries

!-----------------------------------------------------------------------

   call document(subname, 'reading tracer block from ' /&
                       &/ trim(filename))

   do n=1,size(tracer_d_module)

      fieldname = char_blank
      fieldname = trim(tracer_d_module(n)%short_name) /&
               &/ '_CUR'

      call read_field(fmt, filename, fieldname, TRACER_MODULE(:,:,:,n,curtime,:))

      fieldname = char_blank
      fieldname = trim(tracer_d_module(n)%short_name) /&
               &/ '_OLD'

      call read_field(fmt, filename, fieldname, TRACER_MODULE(:,:,:,n,oldtime,:))

   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine rest_read_tracer_block

!***********************************************************************
!BOP
! !IROUTINE: file_read_tracer_block
! !INTERFACE:

 subroutine file_read_tracer_block(init_filename_fmt, init_filename, &
   tracer_d_module, ind_name_table, tracer_init_ext, TRACER_MODULE)

! !DESCRIPTION:
!  read from a file all tracers for a tracer module
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_filename_fmt,   & ! format (bin or nc)
      init_filename          ! default filename

   type (tracer_field), dimension(:), intent(in) :: &
      tracer_d_module        ! descriptors for each tracer

   type(tracer_read), dimension(:), intent(in) :: &
      tracer_init_ext        ! namelist variable for initializing tracers

   type(ind_name_pair), dimension(:) :: &
      ind_name_table

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(:,:,:,:,:,:), intent(inout) :: &
      TRACER_MODULE          ! tracers to be read in

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'passive_tracer_tools:file_read_tracer_block'

   integer (int_kind) :: &
      tracer_cnt,   & ! number of tracers in module
      ind,          & ! tracer index for tracer name from namelist
      iblock,       & ! index for looping over blocks
      n               ! tracer index

   type(tracer_read), dimension(:), allocatable :: &
      tracer_init_int ! namelist variable for initializing tracers

   type (datafile) ::    &
      in_file         ! data file type for init ts file

!-----------------------------------------------------------------------

   tracer_cnt = size(tracer_d_module)

!-----------------------------------------------------------------------
!  initialize internal tracer_init array
!-----------------------------------------------------------------------

   allocate(tracer_init_int(tracer_cnt))
   do n = 1,tracer_cnt
      tracer_init_int(n)%mod_varname  = tracer_d_module(n)%short_name
      tracer_init_int(n)%filename     = init_filename
      tracer_init_int(n)%file_varname = tracer_d_module(n)%short_name
      tracer_init_int(n)%scale_factor = c1
      tracer_init_int(n)%default_val  = c0
      tracer_init_int(n)%file_fmt     = init_filename_fmt
   end do

!-----------------------------------------------------------------------
!  copy non-default values from external tracer_init array
!-----------------------------------------------------------------------

   do n = 1,tracer_cnt
      if (trim(tracer_init_ext(n)%mod_varname) /= 'unknown') then
         ind = name_to_ind(tracer_init_ext(n)%mod_varname, ind_name_table)
         if (ind == 0) then
            call document(subname, 'unknown external varname = ', &
                          trim(tracer_init_ext(n)%mod_varname))
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         end if

         if (trim(tracer_init_ext(n)%filename) /= 'unknown') &
              tracer_init_int(ind)%filename = &
              tracer_init_ext(n)%filename

         if (trim(tracer_init_ext(n)%file_varname) /= 'unknown') &
              tracer_init_int(ind)%file_varname = &
              tracer_init_ext(n)%file_varname

         if (tracer_init_ext(n)%scale_factor /= c1) &
              tracer_init_int(ind)%scale_factor = &
              tracer_init_ext(n)%scale_factor

         if (tracer_init_ext(n)%default_val /= c1) &
              tracer_init_int(ind)%default_val = &
              tracer_init_ext(n)%default_val
      end if
   end do

!-----------------------------------------------------------------------
!  process internal tracer_init array
!-----------------------------------------------------------------------

   do n = 1,tracer_cnt
      if (trim(tracer_init_int(n)%filename) == 'none' .or. &
          trim(tracer_init_int(n)%filename) == 'unknown') then
         call document(subname, 'initializing ' /&
                             &/ trim(tracer_init_int(n)%mod_varname) /&
                             &/ ' to default_val')
         do iblock = 1,nblocks_clinic
            TRACER_MODULE(:,:,:,n,curtime,iblock) =  &
               tracer_init_int(n)%default_val
         enddo
      else
         call document(subname, 'initializing ' /&
                             &/ trim(tracer_init_int(n)%mod_varname) /&
                             &/ ' with ' /&
                             &/ trim(tracer_init_int(n)%file_varname) /&
                             &/ ' from ' /&
                             &/ trim(tracer_init_int(n)%filename))

         call read_field(tracer_init_int(n)%file_fmt,        &
                         tracer_init_int(n)%filename,        &
                         tracer_init_int(n)%file_varname,    &
                         TRACER_MODULE(:,:,:,n,curtime,:))

         do iblock=1,nblocks_clinic
            TRACER_MODULE(:,:,:,n,curtime,iblock) = &
               TRACER_MODULE(:,:,:,n,curtime,iblock)*tracer_init_int(n)%scale_factor
            where (TRACER_MODULE(:,:,:,n,curtime,iblock) < c0) &
               TRACER_MODULE(:,:,:,n,curtime,iblock) = c0
         end do

         if (my_task == master_task) then
            write(stdout,blank_fmt)
            write(stdout,'(a12,a)') ' file read: ', &
               trim(tracer_init_int(n)%filename)
         endif

      end if
      do iblock=1,nblocks_clinic
         TRACER_MODULE(:,:,:,n,oldtime,iblock) = &
            TRACER_MODULE(:,:,:,n,curtime,iblock)
      enddo
   end do

   deallocate(tracer_init_int)

!-----------------------------------------------------------------------
!EOC

 end subroutine file_read_tracer_block

!***********************************************************************
!BOP
! !IROUTINE: read_field_2D
! !INTERFACE:

 subroutine read_field_2D(fmt, filename, fieldname, FIELD, record_length)

! !DESCRIPTION:
!  read 2D field from a file
!  Assumes the field is (nx_global,ny_global), cell centered, and scalar.
!  For binary files, the default external precision is double precision.
!  This can be overridden by passing the desired precision into record_length.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      fmt,                 & ! format (bin or nc)
      filename,            & ! file to read from
      fieldname              ! field to be read

   integer(int_kind), intent(in), optional :: &
      record_length          ! record length type for binary files

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(:,:,:), intent(inout), target :: &
      FIELD                  ! field to be read in

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'passive_tracer_tools:read_field_2D'

   integer(int_kind) :: &
      record_length_loc    ! record length type for binary files

   type (io_field_desc) :: &
      FIELD_DESC           ! IO field descriptors for FIELD

   type (datafile) :: &
      restart_file         ! io file descriptor

   type (io_dim) :: &
      i_dim, j_dim         ! dimension descriptors

!-----------------------------------------------------------------------

   call document(subname, 'reading ' /&
                       &/ trim(fieldname) /&
                       &/ ' from ' /&
                       &/ trim(filename))

   if (present(record_length)) then
      record_length_loc = record_length
   else
      record_length_loc = rec_type_dbl
   endif

   restart_file =                                     &
      construct_file(fmt,                             &
                     full_name=trim(filename),        &
                     record_length=record_length_loc, &
                     recl_words=nx_global*ny_global)

   call data_set(restart_file, 'open_read')

   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)

   FIELD_DESC =                                       &
      construct_io_field(trim(fieldname),             &
                         dim1=i_dim,                  &
                         dim2=j_dim,                  &
                         d2d_array = FIELD)

   call data_set (restart_file, 'define', FIELD_DESC)

   call data_set (restart_file, 'read', FIELD_DESC)

   call destroy_io_field (FIELD_DESC)

   call data_set (restart_file, 'close')

   call destroy_file (restart_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine read_field_2D

!***********************************************************************
!BOP
! !IROUTINE: read_field_3D
! !INTERFACE:

 subroutine read_field_3D(fmt, filename, fieldname, FIELD, record_length)

! !DESCRIPTION:
!  read 3D field from a file
!  Assumes the field is (nx_global,ny_global), cell centered, and scalar.
!  The length of the 3rd dimension is determined by the dimension of FIELD.
!  For binary files, the default external precision is double precision.
!  This can be overridden by passing the desired precision into record_length.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      fmt,                 & ! format (bin or nc)
      filename,            & ! file to read from
      fieldname              ! field to be read

   integer(int_kind), intent(in), optional :: &
      record_length          ! record length type for binary files

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(:,:,:,:), intent(inout), target :: &
      FIELD                  ! field to be read in

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'passive_tracer_tools:read_field_3D'

   integer(int_kind) :: &
      record_length_loc    ! record length type for binary files

   type (io_field_desc) :: &
      FIELD_DESC           ! IO field descriptors for FIELD

   type (datafile) :: &
      restart_file         ! io file descriptor

   type (io_dim) :: &
      i_dim, j_dim, k_dim  ! dimension descriptors

!-----------------------------------------------------------------------

   call document(subname, 'reading ' /&
                       &/ trim(fieldname) /&
                       &/ ' from ' /&
                       &/ trim(filename))

   if (present(record_length)) then
      record_length_loc = record_length
   else
      record_length_loc = rec_type_dbl
   endif

   restart_file =                                     &
      construct_file(fmt,                             &
                     full_name=trim(filename),        &
                     record_length=record_length_loc, &
                     recl_words=nx_global*ny_global)

   call data_set(restart_file, 'open_read')

   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)
   k_dim = construct_io_dim('k', size(FIELD,3))

   FIELD_DESC =                                       &
      construct_io_field(trim(fieldname),             &
                         dim1=i_dim,                  &
                         dim2=j_dim,                  &
                         dim3=k_dim,                  &
                         d3d_array = FIELD)

   call data_set (restart_file, 'define', FIELD_DESC)

   call data_set (restart_file, 'read', FIELD_DESC)

   call destroy_io_field (FIELD_DESC)

   call data_set (restart_file, 'close')

   call destroy_file (restart_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine read_field_3D

!***********************************************************************
!BOP
! !IROUTINE: extract_surf_avg
! !INTERFACE:

 subroutine extract_surf_avg(init_file_fmt,restart_filename,&
                             tracer_cnt,vflux_flag,ind_name_table,&
                             surf_avg)

! !DESCRIPTION:
!  Extract average surface values from restart file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_file_fmt, & ! file format (bin or nc)
      restart_filename ! file name for restart file

   integer (int_kind),intent(in) :: &
      tracer_cnt

   logical (log_kind), dimension(tracer_cnt),intent(in) :: &
      vflux_flag                ! which tracers get virtual fluxes applied

   type(ind_name_pair), dimension(tracer_cnt),intent(in) :: &
      ind_name_table

! !Output PARAMETERS:

   real (r8), dimension(tracer_cnt),intent(out) :: &
      surf_avg         ! average surface tracer values

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   type (datafile) ::&
      restart_file    ! io file descriptor

   integer (int_kind) :: &
      n               ! tracer index

   character (char_len) :: &
      short_name      ! tracer name temporaries


!-----------------------------------------------------------------------

   surf_avg = c0

   restart_file = construct_file(init_file_fmt, &
                                 full_name=trim(restart_filename), &
                                 record_length=rec_type_dbl, &
                                 recl_words=nx_global*ny_global)

   do n = 1, tracer_cnt
      if (vflux_flag(n)) then
         short_name = 'surf_avg_' /&
                   &/ ind_name_table(n)%name
         call add_attrib_file(restart_file, trim(short_name), surf_avg(n))
      endif
   end do

   call data_set(restart_file, 'open_read')

   do n = 1, tracer_cnt
      if (vflux_flag(n)) then
         short_name = 'surf_avg_' /&
                   &/ ind_name_table(n)%name
         call extract_attrib_file(restart_file, trim(short_name), surf_avg(n))
      endif
   end do

   call data_set (restart_file, 'close')

   call destroy_file (restart_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine extract_surf_avg
!*****************************************************************************
!BOP
! !IROUTINE: comp_surf_avg
! !INTERFACE:

 subroutine comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR,tracer_cnt,&
                         vflux_flag,surf_avg)

! !DESCRIPTION:
!  compute average surface tracer values
!
!  avg = sum(SURF_VAL*TAREA) / sum(TAREA)
!  with the sum taken over ocean points only
!
! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:
!real (r8), dimension(nx_block,ny_block,tracer_cnt,max_blocks_clinic), &
   real (r8), dimension(:,:,:,:), &
      intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR

   integer (int_kind),intent(in) :: &
      tracer_cnt                ! number of tracers

   logical (log_kind), dimension(tracer_cnt),intent(in) :: &
      vflux_flag                ! which tracers get virtual fluxes applied

! !Output PARAMETERS:

   real (r8), dimension(tracer_cnt),intent(out) :: &
      surf_avg         ! average surface tracer values


!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,       & ! tracer index
      iblock,  & ! block index
      dcount,  & ! diag counter
      ib,ie,jb,je

   real (r8), dimension(max_blocks_clinic,tracer_cnt) :: &
      local_sums ! array for holding block sums of each diagnostic

   real (r8) :: &
      sum_tmp    ! temp for local sum

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1, &   ! local work space
      TFACT      ! factor for normalizing sums

   type (block) :: &
      this_block ! block information for current block

!-----------------------------------------------------------------------

   local_sums = c0

!jw   !$OMP PARALLEL DO PRIVATE(iblock,this_block,ib,ie,jb,je,TFACT,n,WORK1)
   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      ib = this_block%ib
      ie = this_block%ie
      jb = this_block%jb
      je = this_block%je
      TFACT = TAREA(:,:,iblock)*RCALCT(:,:,iblock)

      do n = 1, tracer_cnt
         if (vflux_flag(n)) then
            WORK1 = p5*(SURF_VALS_OLD(:,:,n,iblock) + &
                        SURF_VALS_CUR(:,:,n,iblock))*TFACT
            local_sums(iblock,n) = sum(WORK1(ib:ie,jb:je))
         endif
      end do
   end do
!jw   !$OMP END PARALLEL DO

   do n = 1, tracer_cnt
      if (vflux_flag(n)) then
         sum_tmp = sum(local_sums(:,n))
         surf_avg(n) = global_sum(sum_tmp,distrb_clinic)/area_t
      endif
   end do

   if(my_task == master_task) then
      write(stdout,*)' Calculating surface tracer averages'
      do n = 1, tracer_cnt
         if (vflux_flag(n)) then
            write(stdout,*) n, surf_avg(n)
         endif
      end do
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_surf_avg
!*****************************************************************************
!BOP
! !IROUTINE: set_tracer_indices
! !INTERFACE:

 subroutine set_tracer_indices(module_string, module_nt,  &
        cumulative_nt, ind_begin, ind_end)

! !DESCRIPTION:
!  set the index bounds of a single passive tracer module
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      module_string

   integer (kind=int_kind), intent(in) ::  &
      module_nt

! !INPUT/OUTPUT PARAMETERS:

   integer (kind=int_kind), intent(inout) ::  &
      cumulative_nt

   integer (kind=int_kind), intent(out) ::  &
      ind_begin, &
      ind_end

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'passive_tracer_tools:set_tracer_indices'

   character (char_len) ::  &
      error_string

!-----------------------------------------------------------------------

   ind_begin = cumulative_nt + 1
   ind_end = ind_begin + module_nt - 1
   cumulative_nt = ind_end

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,*) module_string /&
         &/ ' ind_begin = ', ind_begin
      write(stdout,*) module_string /&
         &/ ' ind_end   = ', ind_end
      write(stdout,delim_fmt)
   end if

   if (cumulative_nt > nt) then
      call document(subname, 'nt', nt)
      call document(subname, 'cumulative_nt', cumulative_nt)
      error_string = 'nt too small for module ' /&
         &/ module_string
      call exit_POP(sigAbort, error_string)
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine set_tracer_indices
!***********************************************************************

!*****************************************************************************
!BOP
! !IROUTINE: tracer_read_init
! !INTERFACE:

 subroutine tracer_read_init(tracer, mod_varname, filename, file_varname, file_fmt, &
      scale_factor, default_val)

! !DESCRIPTION:
!  initialize a tracer_read type to common default values.
!

! !OUTPUT PARAMETERS:
   type(tracer_read), intent (out) :: tracer
! !INPUT PARAMETERS:
   character(char_len), optional, intent(in) :: mod_varname
   character(char_len), optional,  intent(in) :: filename
   character(char_len), optional,  intent(in) :: file_varname
   character(char_len), optional,  intent(in) :: file_fmt
   real(r8), optional,  intent(in) :: scale_factor
   real(r8), optional,  intent(in) :: default_val

!EOP
!BOC

   tracer%mod_varname = 'unknown'
   if (present(mod_varname)) then
      tracer%mod_varname = mod_varname
   end if

   tracer%filename = 'unknown'
   if (present(filename)) then
      tracer%filename = filename
   end if

   tracer%file_varname = 'unknown'
   if (present(file_varname)) then
      tracer%file_varname = file_varname
   end if

   tracer%file_fmt = 'bin'
   if (present(file_fmt)) then
      tracer%file_fmt = file_fmt
   end if

   tracer%scale_factor = c1
   if (present(scale_factor)) then
      tracer%scale_factor = scale_factor
   end if

   tracer%default_val = c0
   if (present(default_val)) then
      tracer%default_val = default_val
   end if

!-----------------------------------------------------------------------
!EOC
 end subroutine tracer_read_init
!***********************************************************************


 end module passive_tracer_tools

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
