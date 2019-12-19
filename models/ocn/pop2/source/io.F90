!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module io

!BOP
! !MODULE: io
!
! !DESCRIPTION:
!  This module provides a generic parallel input/output interface
!  for writing arrays.
!
! !REVISION HISTORY:
!  SVN:$Id: io.F90 63715 2014-09-22 22:47:53Z klindsay $

! !USES:

   use kinds_mod
   use blocks
   use communicate
   use broadcast
   use exit_mod
   use domain
   use constants
   use io_netcdf
   use io_binary
   use io_types

   implicit none
   public  ! to get io_types without having to explicitly use io_types
           ! module directly
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: data_set

!EOP
!BOC
!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: data_set
! !INTERFACE:

 subroutine data_set (data_file, operation, io_field, fieldname, field_exists)

! !DESCRIPTION:
!  This routine is the main interface for array and file io functions,
!  including read, write, open, close.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent (in)   :: operation
   character (*), intent (in), optional :: fieldname

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: data_file
   type (io_field_desc), intent (inout), optional :: io_field

! !OUTPUT PARAMETERS:

   logical (log_kind), intent (out), optional :: field_exists

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  select operation to perform
!
!-----------------------------------------------------------------------

   select case (trim(operation))

!-----------------------------------------------------------------------
!
!  open for reading
!
!-----------------------------------------------------------------------

   case ('open_read')

      if (data_file%data_format=='bin') then
         call open_read_binary(data_file)
      else if (data_file%data_format=='nc') then
         call open_read_netcdf(data_file)
      endif

!-----------------------------------------------------------------------
!
!  Open means open for write.  We also at this time write any global
!  attributes.
!
!-----------------------------------------------------------------------

   case ('open')

      if (data_file%data_format=='bin') then
         call open_binary(data_file)
      else if (data_file%data_format=='nc') then
         call open_netcdf(data_file)
      endif

!-----------------------------------------------------------------------
!
!  close a data file
!
!-----------------------------------------------------------------------

   case ('close')

      if (data_file%data_format=='bin') then
         call close_binary(data_file)
      else if (data_file%data_format=='nc') then
         call close_netcdf(data_file)
      endif

!-----------------------------------------------------------------------
!
!  flush data to a data file
!
!-----------------------------------------------------------------------

   case ('flush')

      if (data_file%data_format=='bin') then
        !*** tbd
      else if (data_file%data_format=='nc') then
         call sync_netcdf(data_file)
      endif

!-----------------------------------------------------------------------
!
!  determine if a field exists
!
!-----------------------------------------------------------------------

   case ('field_exists')

      if (.not.present(fieldname)) then
         call exit_POP(sigAbort,'data_file field_exists: missing fieldname arg')
      end if

      if (.not.present(field_exists)) then
         call exit_POP(sigAbort,'data_file field_exists: missing field_exists arg')
      end if

      if (data_file%data_format=='bin') then
         call field_exists_binary(data_file,fieldname,field_exists)
      else if (data_file%data_format=='nc') then
         call field_exists_netcdf(data_file,fieldname,field_exists)
      endif

!-----------------------------------------------------------------------
!
!  define an io field
!
!-----------------------------------------------------------------------

   case ('define')

      if (.not.present(io_field)) then
         call exit_POP(sigAbort, &
                       'data_file define: missing io_field arg')
      end if

      if (data_file%data_format=='bin') then
         call define_field_binary(data_file,io_field)
      else if (data_file%data_format=='nc') then
         call define_field_netcdf(data_file,io_field)
      endif

!-----------------------------------------------------------------------
!
!  write an io field
!
!-----------------------------------------------------------------------

   case ('write')

      if (.not.present(io_field)) then
         call exit_POP(sigAbort,'data_file write: missing io_field arg')
      end if

      if (data_file%data_format=='bin') then
         call write_field_binary(data_file,io_field)
      else if (data_file%data_format=='nc') then
         call write_field_netcdf(data_file,io_field)
      endif

!-----------------------------------------------------------------------
!
!  read an io field
!
!-----------------------------------------------------------------------

   case ('read')

      if (.not.present(io_field)) then
         call exit_POP(sigAbort,'data_file read: missing io_field arg')
      end if

      if (data_file%data_format=='bin') then
         call read_field_binary(data_file,io_field)
      else if (data_file%data_format=='nc') then
         call read_field_netcdf(data_file,io_field)
      endif

!-----------------------------------------------------------------------
!
!  unknown operation
!
!-----------------------------------------------------------------------

   case default

      if (my_task == master_task) &
         write(stdout,*) 'data_set operation: ',trim(operation)
      call exit_POP(sigAbort,'data_set: Unknown operation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine data_set

!***********************************************************************


 end module io

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
