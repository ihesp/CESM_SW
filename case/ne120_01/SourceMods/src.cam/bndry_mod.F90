#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module bndry_mod
  use parallel_mod, only : abortmp,iam, HME_BNDRY_P2P, HME_BNDRY_MASHM, HME_BNDRY_A2A, HME_BNDRY_A2AO, &
                           HME_BNDRY_GET1, HME_BNDRY_GET2, HME_BNDRY_PUT1, HME_BNDRY_PUT2
  use edgetype_mod, only : Ghostbuffer3D_t
  use thread_mod, only : omp_in_parallel, omp_get_thread_num
  use gbarrier_mod, only: gbarrier
  use perf_mod, only : t_startf, t_stopf ! EXTERNAL
  implicit none
  private
#ifdef _MPI
#include <mpif.h>
#endif

  public :: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation
  public :: ghost_exchangeV
  public :: bndry_exchangeS
  public :: bndry_exchangeS_start
  public :: bndry_exchangeS_finish
!  public :: ghost_exchangev3d
  public :: sort_neighbor_buffer_mapping



  interface bndry_exchangeV
     module procedure bndry_exchangeV_threaded
     module procedure bndry_exchangeV_nonthreaded
     module procedure long_bndry_exchangeV_nonth
  end interface

  interface bndry_exchangeS
     module procedure bndry_exchangeS_threaded 
     module procedure bndry_exchangeS_nonthreaded
  end interface

  interface bndry_exchangeS_finish
     module procedure bndry_exchangeS_threaded_finish
     module procedure bndry_exchangeS_nonthreaded_finish
  end interface

  interface bndry_exchangeS_start
     module procedure bndry_exchangeS_threaded_start
     module procedure bndry_exchangeS_nonthreaded_start
  end interface

contains 

  subroutine bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version,hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location
!
!    type (Schedule_t),pointer         :: pSchedule
!    type (Cycle_t),pointer            :: pCycle
!    integer                           :: icycle,ierr
!    integer                           :: length
!    integer                           :: iptr,source,nlyr
!    integer                           :: nSendCycles,nRecvCycles
    integer                           :: ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_a2a'
    character(len=80)                 :: locstring
!
!    logical(kind=log_kind), parameter :: Debug=.FALSE.
!    logical :: ompthreadMissmatch
!
!    integer                           :: i,j
!    integer :: lenMovePtr
    integer :: request
    integer :: lstatus(HME_status_size)

! Neighborhood collectives are only in MPI3 and up
#ifdef _MPI3
   if(ithr == 0) then 

      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsFull,buffer%sdisplsFull,MPIreal_t, &
                     buffer%receive,buffer%rcountsFull,buffer%rdisplsFull,MPIreal_t,par%commGraphFull,request,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,subname,': Error after call to MPI_neighbor_alltoallv: ',errorstring
      endif

      locstring = TRIM(subname) // ': ' // TRIM(location)
      ! location 1 for copyBuffer
      call copyBuffer(nthreads,ithr,buffer,locstring)

      call MPI_wait(request,lstatus,ierr)
   else

      locstring = TRIM(subname) // ': ' // TRIM(location)
      call copyBuffer(nthreads,ithr,buffer,locstring)

   endif
#else
    call abortmp('bndry_exchange_a2a requires MPI-3 feature support')
#endif

  end subroutine bndry_exchange_a2a
#define SW_SLAVE_COPYBUFFER
  subroutine copyBuffer(nthreads,ithr,buffer,location)
    use edgetype_mod, only : Edgebuffer_t
    use spmd_utils, only: iam
    integer :: nthreads
    integer :: ithr
    type (EdgeBuffer_t)          :: buffer
    character(len=80)            :: location
    logical ::  ompThreadMissmatch
    integer lenMovePtr, iptr,length,i,j
#ifdef SW_SLAVE_COPYBUFFER
	integer,parameter ::M1=3200,M2=256000
	type str_par
		integer(kind=8)	:: ptr_a,ptr_b
		integer					:: n1,n2
	end type str_par
	type(str_par) :: par 
	external		:: slave_copybuffer
	integer :: idx,nn,mm,jump
#endif
		call t_startf("copyBuffer")
		ompThreadMissmatch = .false.
    lenMovePtr = size(buffer%moveptr)
    if ( lenMOveptr .ne. nthreads) then
      ompthreadMissmatch = .true.
      write(*,30) TRIM(location), lenMoveptr, nthreads
    endif
#ifndef SW_SLAVE_COPYBUFFER
    if (.not. ompthreadMissmatch) then
      iptr   = buffer%moveptr(ithr+1)
      length = buffer%moveLength(ithr+1)
			!print*,length
			!print*,"id=   ",iam," length=    ",length
      if(length>0) then
        do i=0,length-1
           buffer%receive(iptr+i) = buffer%buf(iptr+i)
        enddo
      endif
    else if(ompthreadMissmatch .and. ithr == 0) then
       do j=1,lenMovePtr
          iptr   = buffer%moveptr(j)
          length = buffer%moveLength(j)
					!print*,length
          if(length>0) then
             do i=0,length-1
                buffer%receive(iptr+i) = buffer%buf(iptr+i)
             enddo
          endif
       enddo
    endif
#else
	  if (.not. ompthreadMissmatch) then
      iptr   = buffer%moveptr(ithr+1)
      length = buffer%moveLength(ithr+1)
      if(length>0) then
				!  do i=0,length-1
				!     buffer%receive(iptr+i) = buffer%buf(iptr+i)
				!  enddo
				if (length<=M1) then
					do i=0,length-1
						buffer%receive(iptr+i) = buffer%buf(iptr+i)
					enddo
				else if(length>M1 .and. length<=M2) then
					par%ptr_a=loc(buffer%buf(iptr))
					par%ptr_b=loc(buffer%receive(iptr))
					par%n1=length/64
					par%n2=mod(length,par%n1*64)
					call athread_spawn(slave_copybuffer,par)
					call athread_join()
				else
					nn=length/M2
					mm=mod(length,M2)
					if (mm<M1) then
						!并行部分
						jump=M2
						do idx=0,nn-1
							par%ptr_a=loc(buffer%buf(iptr+jump*idx))
							par%ptr_b=loc(buffer%receive(iptr+jump*idx))
							par%n1=4000
							par%n2=0
							call athread_spawn(slave_copybuffer,par)
							call athread_join()
						enddo
						!串行部分
						do i=iptr+nn*jump,mm+iptr+nn*jump
							buffer%receive(i) = buffer%buf(i)
						enddo
					else  
						!满载并行
						jump=M2
						do idx=0,nn-1
							par%ptr_a=loc(buffer%buf(iptr+jump*idx))
							par%ptr_b=loc(buffer%receive(iptr+jump*idx))
							par%n1=4000
							par%n2=0
							call athread_spawn(slave_copybuffer,par)
							call athread_join()
						enddo
						!非满载并行	
						par%ptr_a=loc(buffer%buf(iptr+jump*nn))
						par%ptr_b=loc(buffer%receive(iptr+jump*nn))
						par%n1=(length-jump*nn)/64
						par%n2=mod((length-jump*nn),par%n1*64)
						call athread_spawn(slave_copybuffer,par)
						call athread_join()
					endif 
				endif
			endif
    else if(ompthreadMissmatch .and. ithr == 0) then
       do j=1,lenMovePtr
          iptr   = buffer%moveptr(j)
          length = buffer%moveLength(j)
          if(length>0) then
            ! do i=0,length-1
            !    buffer%receive(iptr+i) = buffer%buf(iptr+i)
            ! enddo
						if (length<=M1) then
							do i=0,length-1
								buffer%receive(iptr+i) = buffer%buf(iptr+i)
							enddo
						else if(length>M1 .and. length<=M2) then
							par%ptr_a=loc(buffer%buf(iptr))
							par%ptr_b=loc(buffer%receive(iptr))
							par%n1=length/64
							par%n2=mod(length,par%n1*64)
							call athread_spawn(slave_copybuffer,par)
							call athread_join()
						else
							nn=length/M2
							mm=mod(length,M2)
							if (mm<M1) then
								!并行部分
								jump=M2
								do idx=0,nn-1
									par%ptr_a=loc(buffer%buf(iptr+jump*idx))
									par%ptr_b=loc(buffer%receive(iptr+jump*idx))
									par%n1=4000
									par%n2=0
									call athread_spawn(slave_copybuffer,par)
									call athread_join()
								enddo
								!串行部分
								do i=iptr+nn*jump,mm+iptr+nn*jump
									buffer%receive(i) = buffer%buf(i)
								enddo
							else  
								!满载并行
								jump=M2
								do idx=0,nn-1
									par%ptr_a=loc(buffer%buf(iptr+jump*idx))
									par%ptr_b=loc(buffer%receive(iptr+jump*idx))
									par%n1=4000
									par%n2=0
									call athread_spawn(slave_copybuffer,par)
									call athread_join()
								enddo
								!非满载并行	
								par%ptr_a=loc(buffer%buf(iptr+jump*nn))
								par%ptr_b=loc(buffer%receive(iptr+jump*nn))
								par%n1=(length-jump*nn)/64
								par%n2=mod((length-jump*nn),par%n1*64)
								call athread_spawn(slave_copybuffer,par)
								call athread_join()
							endif 
						endif
          endif
       enddo
    endif
#endif

30  format(a,'Potential performance issue: ',a,'LenMoveptr,nthreads: ',2(i3))
		call t_stopf("copyBuffer")
  end subroutine copyBuffer

  subroutine bndry_exchange_put1(par,nthreads,ithr,buffer,location)

    use kinds, only : log_kind,long_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version, mpi_lock_shared, hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: icycle,ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_put1'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.TRUE.

    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)
    integer :: iptr,tag,length,source,dest,nRecvCycles, nSendCycles
    integer(kind=long_kind) :: disp

   pSchedule => schedule(1)
   if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    call MPI_Win_fence(0,buffer%win,ierr)
    do icycle=1,nSendCycles
       pCycle         => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle)
       ! disp            = 0  !FIXME: This needs to point to the correct spot in the remove ranks memory
       disp            = INT(buffer%putDisplsFull(icycle),kind=long_kind)  !FIXME: This needs to point to the correct spot in the remove ranks memory
!       if(Debug .and. (source == 103 .or. iam == 104)) print *,'IAM: ',iam, subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'PTR: ',iptr, 'TAG: ',tag
!       call MPI_Win_lock(MPI_LOCK_SHARED,source,0,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_win_lock'
       call MPI_Put(buffer%buf(iptr+1),length,MPIreal_t,dest,disp,length,MPIreal_t,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_Get'
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Get: ',errorstring
       endif
    end do    ! icycle

    ! location 3 for copyBuffer
    call copyBuffer(nthreads,ithr,buffer,locstring)

    call MPI_Win_fence(0,buffer%win,ierr)

   else

    locstring = TRIM(subname) // ': ' // TRIM(location) 
    !Copy buffer for ithr!=0
    call copyBuffer(nthreads,ithr,buffer,locstring)

   endif

!   stop 'At end of bndry_exchange_put1'

  end subroutine bndry_exchange_put1

  subroutine bndry_exchange_put2(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind,long_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version, mpi_lock_shared, hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: icycle,ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_put2'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.TRUE.

    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)
    integer :: iptr,tag,length,source,dest, nRecvCycles, nSendCycles
    integer(kind=long_kind) :: disp

   pSchedule => schedule(1)
   if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    call MPI_Win_post(par%groupGraphFull,0,buffer%win,ierr)
    !==================================================
    ! Get the data from the remote process
    !==================================================
    call MPI_win_start(par%groupGraphFull,0,buffer%win,ierr)
    do icycle=1,nSendCycles
       pCycle         => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle)
       disp            = INT(buffer%putDisplsFull(icycle),kind=long_kind)
       call MPI_Put(buffer%buf(iptr+1),length,MPIreal_t,dest,disp,length,MPIreal_t,buffer%win,ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Get: ',errorstring
       endif
    end do    ! icycle

    ! location 3 for copyBuffer
    call copyBuffer(nthreads,ithr,buffer,locstring)

    call MPI_Win_complete(buffer%win,ierr)
    call MPI_Win_wait(buffer%win,ierr)

   else

    locstring = TRIM(subname) // ': ' // TRIM(location) 
    !Copy buffer for ithr!=0
    call copyBuffer(nthreads,ithr,buffer,locstring)

   endif

  end subroutine bndry_exchange_put2

  subroutine bndry_exchange_get1(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind,long_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version, mpi_lock_shared, hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: icycle,ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_get1'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.TRUE.

    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)
    integer :: iptr,tag,length,source,nRecvCycles, nSendCycles
    integer(kind=long_kind) :: disp

   pSchedule => schedule(1)
   if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    call MPI_Win_fence(0,buffer%win,ierr)
!    call MPI_Win_post(par%commGraphfull,0,buffer%win,ierr)
    !==================================================
    ! Get the data from the remote process
    !==================================================
!    call MPI_win_start(par%commGraphfull,0,buffer%win,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle)
       disp            = INT(buffer%getDisplsFull(icycle),kind=long_kind)
       call MPI_Get(buffer%receive(iptr+1),length,MPIreal_t,source,disp,length,MPIreal_t,buffer%win,ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Get: ',errorstring
       endif
    end do    ! icycle

    ! location 3 for copyBuffer
    call copyBuffer(nthreads,ithr,buffer,locstring)

    call MPI_Win_fence(0,buffer%win,ierr)

   else

    locstring = TRIM(subname) // ': ' // TRIM(location) 
    !Copy buffer for ithr!=0
    call copyBuffer(nthreads,ithr,buffer,locstring)

   endif

!   stop 'At end of bndry_exchange_get1'

  end subroutine bndry_exchange_get1

  subroutine bndry_exchange_get2(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind,long_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version, mpi_lock_shared, hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: icycle,ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_get2'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.TRUE.

    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)
    integer :: iptr,tag,length,source,nRecvCycles, nSendCycles
    integer(kind=long_kind) :: disp

   pSchedule => schedule(1)
   if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

!    call MPI_Win_fence(0,buffer%win,ierr)
    call MPI_Win_post(par%groupGraphFull,0,buffer%win,ierr)
    !==================================================
    ! Get the data from the remote process
    !==================================================
    call MPI_win_start(par%groupGraphFull,0,buffer%win,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle)
       ! disp            = 0  !FIXME: This needs to point to the correct spot in the remove ranks memory
       disp            = INT(buffer%getDisplsFull(icycle),kind=long_kind)  !FIXME: This needs to point to the correct spot in the remove ranks memory
!       if(Debug .and. (source == 103 .or. iam == 104)) print *,'IAM: ',iam, subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'PTR: ',iptr, 'TAG: ',tag
!       call MPI_Win_lock(MPI_LOCK_SHARED,source,0,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_win_lock'
       call MPI_Get(buffer%receive(iptr+1),length,MPIreal_t,source,disp,length,MPIreal_t,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_Get'
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Get: ',errorstring
       endif
!       call MPI_win_unlock(source,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_win_unlock'
    end do    ! icycle

!    if(Debug) print *,'IAM: ', iam, 'Before copyBuffer'
    ! location 3 for copyBuffer
    call copyBuffer(nthreads,ithr,buffer,locstring)

    call MPI_Win_complete(buffer%win,ierr)
    call MPI_Win_wait(buffer%win,ierr)

   else

    locstring = TRIM(subname) // ': ' // TRIM(location) 
    !Copy buffer for ithr!=0
    call copyBuffer(nthreads,ithr,buffer,locstring)

   endif

!   stop 'At end of bndry_exchange_get2'

  end subroutine bndry_exchange_get2

  subroutine bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version,hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

!    type (Schedule_t),pointer         :: pSchedule
!    type (Cycle_t),pointer            :: pCycle
!    integer                           :: icycle,ierr
!    integer                           :: length
!    integer                           :: iptr,source,nlyr
!    integer                           :: nSendCycles,nRecvCycles
    integer                           :: ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_a2ao'
    character(len=80)                 :: locstring
!
!    logical(kind=log_kind), parameter :: Debug=.FALSE.
!    logical :: ompthreadMissmatch
!
!    integer                           :: i,j
!    integer :: lenMovePtr
    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)

! Neighborhood collectives are only in MPI3 and up
#ifdef _MPI3

   if(ithr == 0) then 

      ! Start Inter-node communication
      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsInter,buffer%sdisplsInter,MPIreal_t, &
                     buffer%receive,buffer%rcountsInter,buffer%rdisplsInter,MPIreal_t,par%commGraphInter,requestInter,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,subname,': Error after call to MPI_neighbor_alltoallv: ',errorstring
      endif
      locstring = TRIM(subname) // ': ' // TRIM(location) 

      ! Start Intra-node communication
      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsIntra,buffer%sdisplsIntra,MPIreal_t, &
                     buffer%receive,buffer%rcountsIntra,buffer%rdisplsIntra,MPIreal_t,par%commGraphIntra,requestIntra,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,subname,': Error after call to MPI_neighbor_alltoallv: ',errorstring
      endif

      ! Finish the Intra-node communication 
      call MPI_wait(requestIntra,lstatus,ierr)

      ! location 3 for copyBuffer
      call copyBuffer(nthreads,ithr,buffer,locstring)

      ! Finish the Inter-node communication 
      call MPI_wait(requestInter,lstatus,ierr)

   else

      locstring = TRIM(subname) // ': ' // TRIM(location) 
      !Copy buffer for ithr!=0
      call copyBuffer(nthreads,ithr,buffer,locstring)

   endif
#else
    call abortmp('bndry_exchange_a2ao requires MPI-3 feature support')
#endif

  end subroutine bndry_exchange_a2ao

  subroutine bndry_exchangeV_p2p(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
!    integer                           :: imptr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchangeV_p2p'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.FALSE.
!    logical :: ompthreadMissmatch
!
!    integer                           :: i,j
!    integer :: lenMovePtr

    pSchedule => Schedule(1)
!    nlyr = buffer%nlyr
!    ompthreadMissmatch = .false. 

!    lenMovePtr = size(buffer%moveptr)

  if(ithr == 0) then 
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle)
       if(Debug) print *,'IAM: ',iam, subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'PTR: ',iptr, 'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr+1),length,MPIreal_t, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Fire off the sends
    !==================================================
    do icycle=1,nSendCycles
       pCycle         => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle)
       if(Debug) print *,'IAM: ',iam, subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'PTR: ',iptr, 'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr+1),length,MPIreal_t,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle
    
      !locstring = TRIM(subname) // ': ' // TRIM(location) 
      !add by chenyuhu
      locstring = TRIM(subname)
      call copyBuffer(nthreads,ithr,buffer,locstring)

      call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
      call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
  else  ! else for non-master threads

      !locstring = TRIM(subname) // ': ' // TRIM(location) 
      !add by chenyuhu
      locstring = TRIM(subname)
      !Copy buffer for ithr!=0
      call copyBuffer(nthreads,ithr,buffer,locstring)

  endif

  end subroutine bndry_exchangeV_p2p

  subroutine bndry_exchangeS_p2p(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif

    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*), parameter       :: subname = 'bndry_exchangeS_p2p'
    character(len=80)                 :: locstring
    logical(kind=log_kind), parameter :: Debug=.FALSE.

    integer                           :: i,j
    logical :: ompthreadMissmatch
    integer :: lenMovePtr

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr
    ompthreadMissmatch = .FALSE.

    lenMovePtr = size(buffer%moveptr)

  if(ithr == 0) then 
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles


    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle) + 1
       if(Debug) print *,subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr),length,MPIreal_t,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle) + 1
       if(Debug) print *,subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr),length,MPIreal_t, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle

    locstring = TRIM(subname) // ': ' // TRIM(location)
    call copyBuffer(nthreads,ithr,buffer,locstring)
    call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
  else
    locstring = TRIM(subname) // ': ' // TRIM(location)
    call copyBuffer(nthreads,ithr,buffer,locstring)
  endif

  end subroutine bndry_exchangeS_p2p

  subroutine bndry_exchangeS_p2p_start(par,nthreads,ithr,buffer,location)

    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character (len=*), optional :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchangeS_p2p_start'
    logical(kind=log_kind), parameter :: Debug=.FALSE.

    integer                           :: i,j, lenMovePtr
    logical :: ompthreadMissmatch

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr
    ompthreadMissmatch = .FALSE.

    lenMovePtr = size(buffer%moveptr)

  if(ithr == 0) then 
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle) + 1
       if(Debug) print *,subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr),length,MPIreal_t,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle) + 1
       if(Debug) print *,subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr),length,MPIreal_t, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle
  endif
    
  end subroutine bndry_exchangeS_p2p_start

  subroutine bndry_exchangeS_p2p_finish(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchangeS_p2p_finish'
    character(len=80)                 :: locstring
    logical(kind=log_kind), parameter :: Debug=.FALSE.

    integer                           :: i,j
    logical :: ompthreadMissmatch
    integer        :: lenMovePtr

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr

  locstring = TRIM(subname) // ': ' // TRIM(location)
  call copyBuffer(nthreads,ithr,buffer,locstring)

  if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
  endif

  end subroutine bndry_exchangeS_p2p_finish

  subroutine long_bndry_exchangeV_nonth(par,buffer)
    use kinds, only : log_kind
    use edgetype_mod, only : LongEdgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    type (LongEdgeBuffer_t)           :: buffer

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*),       parameter :: subname = 'long_bndry_exchangeV_nonth'
    logical(kind=log_kind), parameter :: Debug=.FALSE.

    integer                           :: i

#ifdef _MPI
    if(omp_in_parallel()) then
       print *,subname,': Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
       print *,'                Results are not predictable!!            '
    endif

    ! Setup the pointer to proper Schedule
    pSchedule => Schedule(1)
    nlyr = buffer%nlyr

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles


    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP
       !DBG if(Debug) print *,subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(1,iptr),length,MPIinteger_t,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP
       !DBG if(Debug) print *,subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(1,iptr),length,MPIinteger_t, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle


    !==================================================
    !  Wait for all the receives to complete
    !==================================================

    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       length             = pCycle%lengthP
       iptr            = pCycle%ptrP
       do i=0,length-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
    end do   ! icycle

#endif

  end subroutine long_bndry_exchangeV_nonth
  !********************************************************************************
  !
  !********************************************************************************
 subroutine bndry_exchangeV_threaded(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t, PrintHybrid
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (hybrid_t)               :: hybrid
    type (EdgeBuffer_t)           :: buffer
    character(len=*),   parameter :: subname = 'bndry_exchangeV_threaded'
    character(len=*), optional        :: location
    integer :: localsense

!    call t_adj_detailf(+2)
    call gbarrier(buffer%gbarrier, hybrid%ithr)
    if(buffer%bndry_type == HME_BNDRY_A2A) then 
       call bndry_exchange_a2a(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_A2AO) then 
       call bndry_exchange_a2ao(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_GET1) then 
       call bndry_exchange_get1(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_GET2) then 
       call bndry_exchange_get2(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_PUT1) then 
       call bndry_exchange_put1(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_PUT2) then 
       call bndry_exchange_put2(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else
       call bndry_exchangeV_p2p(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    endif
    call gbarrier(buffer%gbarrier, hybrid%ithr)
!    call t_adj_detailf(-2)

  end subroutine bndry_exchangeV_threaded

  subroutine bndry_exchangeV_nonthreaded(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    ! local
    character(len=*), parameter :: subname = 'bndry_exchangeV_nonthreaded'
    integer                     :: ithr
    integer                     :: nthreads
    integer                     :: localsense

!    call t_adj_detailf(+2)
    ithr=0
    nthreads = 1
!    localsense = 0
    if(buffer%bndry_type == HME_BNDRY_A2A) then 
       call bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then 
       call bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET1) then 
       call bndry_exchange_get1(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET2) then 
       call bndry_exchange_get2(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT1) then 
       call bndry_exchange_put1(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT2) then 
       call bndry_exchange_put2(par,nthreads,ithr,buffer,location)
    else
       call bndry_exchangeV_p2p(par,nthreads,ithr,buffer,location)
    endif
!    call t_adj_detailf(-2)

  end subroutine bndry_exchangeV_nonthreaded

 subroutine bndry_exchangeS_threaded(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchangeS_threaded'

!    call t_adj_detailf(+2)
    call gbarrier(buffer%gbarrier, hybrid%ithr)
    if(buffer%bndry_type == HME_BNDRY_A2A) then
       call bndry_exchange_a2a(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then 
       call bndry_exchange_a2ao(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET1) then 
       call bndry_exchange_get1(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET2) then 
       call bndry_exchange_get2(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT1) then 
       call bndry_exchange_put1(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT2) then 
       call bndry_exchange_put2(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else
       call bndry_exchangeS_p2p(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    endif
    call gbarrier(buffer%gbarrier, hybrid%ithr)
!    call t_adj_detailf(-2)

 end subroutine bndry_exchangeS_threaded

 subroutine bndry_exchangeS_threaded_start(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchangeS_threaded_start'

!    call t_adj_detailf(+2)
    call gbarrier(buffer%gbarrier, hybrid%ithr)
    call bndry_exchangeS_p2p_start(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
!    call t_adj_detailf(-2)

 end subroutine bndry_exchangeS_threaded_start

 subroutine bndry_exchangeS_threaded_finish(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchangeS_threaded_finish'

!    call t_adj_detailf(+2)
    call bndry_exchangeS_p2p_finish(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    call gbarrier(buffer%gbarrier, hybrid%ithr)
!    call t_adj_detailf(-2)

 end subroutine bndry_exchangeS_threaded_finish

 subroutine bndry_exchangeS_nonthreaded(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    integer                     :: ithr
    integer                     :: nthreads
!    integer                     :: localsense
    character(len=*), parameter :: subname = 'bndry_exchangeS_nonthreaded'

!    call t_adj_detailf(+2)
    ithr=0
    nthreads = 1
    if(buffer%bndry_type == HME_BNDRY_A2A) then 
       call bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then 
       call bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET1) then 
       call bndry_exchange_get1(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET2) then 
       call bndry_exchange_get2(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT1) then 
       call bndry_exchange_put1(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT2) then 
       call bndry_exchange_put2(par,nthreads,ithr,buffer,location)
    else
       call bndry_exchangeS_p2p(par,nthreads,ithr,buffer,location)
    endif
!    call t_adj_detailf(-2)

  end subroutine bndry_exchangeS_nonthreaded

 subroutine bndry_exchangeS_nonthreaded_start(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character (len=*), optional :: location

    integer                     :: ithr
    integer                     :: nthreads
    character(len=*), parameter :: subname = 'bndry_exchangeS_nonthreaded_start'

!    call t_adj_detailf(+2)
!    !$OMP BARRIER
    ithr=0
    nthreads=1
    call bndry_exchangeS_p2p_start(par,nthreads,ithr,buffer,location)
!    call t_adj_detailf(-2)

 end subroutine bndry_exchangeS_nonthreaded_start

 subroutine bndry_exchangeS_nonthreaded_finish(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)                 :: par
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character (len=*), optional :: location
    integer :: nthreads

    character(len=*), parameter :: subname = 'bndry_exchangeS_nonthreaded_finish'

!    call t_adj_detailf(+2)
    ithr=0
    nthreads=1
    call bndry_exchangeS_p2p_finish(par,nthreads,ithr,buffer,location)
!    call t_adj_detailf(-2)

  end subroutine bndry_exchangeS_nonthreaded_finish

  subroutine ghost_exchangeVfull(par,ithr,buffer)
!
!   MT 2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, parallel_t
#else
    use parallel_mod, only : abortmp, parallel_t
#endif
    implicit none
    type (parallel_t)                :: par
    integer                          :: ithr     ! hybrid%ithr 0 if called outside threaded region

!    type (hybrid_t)                 :: hybrid
    type (GhostBuffer3D_t)           :: buffer

    type (Schedule_t),pointer        :: pSchedule
    type (Cycle_t),pointer           :: pCycle
    integer                          :: dest,length,tag
    integer                          :: icycle,ierr
    integer                          :: iptr,source,nlyr
    integer                          :: nSendCycles,nRecvCycles
    integer                          :: errorcode,errorlen
    character(len=*), parameter      :: subname = 'ghost_exchangeVfull'
    character*(80) errorstring

    integer                          :: i,i1,i2
    logical(kind=log_kind),parameter :: Debug = .FALSE.

    !$OMP BARRIER
    if(ithr == 0) then 


#ifdef _MPI
       ! Setup the pointer to proper Schedule
       pSchedule => Schedule(1)
       nlyr = buffer%nlyr

       nSendCycles = pSchedule%nSendCycles
       nRecvCycles = pSchedule%nRecvCycles

       !==================================================
       !  Fire off the sends
       !==================================================
       do icycle=1,nSendCycles
          pCycle      => pSchedule%SendCycle(icycle)
          dest            = pCycle%dest - 1
          length      = nlyr * pCycle%lengthP_ghost * buffer%elem_size
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(buffer%buf(1,1,1,iptr),length,MPIreal_t,dest,tag,par%comm,Srequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,subname,': Error after call to MPI_Isend: ',errorstring
          endif
       end do    ! icycle

       !==================================================
       !  Post the Receives 
       !==================================================
       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          source          = pCycle%source - 1
          length      = nlyr * pCycle%lengthP_ghost * buffer%elem_size
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(buffer%receive(1,1,1,iptr),length,MPIreal_t, &
               source,tag,par%comm,Rrequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,subname,': Error after call to MPI_Irecv: ',errorstring
          endif
       end do    ! icycle


       !==================================================
       !  Wait for all the receives to complete
       !==================================================

       call MPI_Waitall(nSendCycles,Srequest,status,ierr)
       call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          length             = pCycle%lengthP_ghost
          iptr            = pCycle%ptrP_ghost
          do i=0,length-1
             buffer%buf(:,:,1:nlyr,iptr+i) = buffer%receive(:,:,1:nlyr,iptr+i)
          enddo
       end do   ! icycle


#endif
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER

  end subroutine ghost_exchangeVfull

  ! ===========================================
  !  GHOST_EXCHANGEV:
  !  Author: Christoph Erath
  !  derived from bndry_exchange, but copies an entire
  !             element of ghost cell information, including corner
  !             elements.  Requres cubed-sphere grid
  ! =========================================
 subroutine ghost_exchangeV(hybrid,buffer,nhc,npoints,ntrac)
!
!   2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use edgetype_mod, only : Ghostbuffertr_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : abortmp
#endif
    implicit none

    type (hybrid_t)                  :: hybrid
    type (GhostBuffertr_t)           :: buffer
    integer                          :: nhc,npoints,ntrac

    type (Schedule_t),pointer        :: pSchedule
    type (Cycle_t),pointer           :: pCycle
    integer                          :: dest,length,tag
    integer                          :: icycle,ierr
    integer                          :: iptr,source,nlyr
    integer                          :: nSendCycles,nRecvCycles
    integer                          :: errorcode,errorlen
    character(len=*), parameter      :: subname = 'ghost_exchangeV'
    character*(80) errorstring

    integer                          :: i,i1,i2
    logical(kind=log_kind),parameter :: Debug = .FALSE.


    !$OMP BARRIER
    if(hybrid%ithr == 0) then 

#ifdef _MPI
       ! Setup the pointer to proper Schedule
       pSchedule => Schedule(1)
       nlyr = buffer%nlyr
              
       nSendCycles = pSchedule%nSendCycles
       nRecvCycles = pSchedule%nRecvCycles

       !==================================================
       !  Fire off the sends
       !==================================================
       do icycle=1,nSendCycles
          pCycle      => pSchedule%SendCycle(icycle)
          dest            = pCycle%dest - 1
          length      = nlyr * ntrac * pCycle%lengthP_ghost*nhc*npoints
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(buffer%buf(1,1,1,1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,subname,': Error after call to MPI_Isend: ',errorstring
          endif
       end do    ! icycle

       !==================================================
       !  Post the Receives 
       !==================================================
       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          source          = pCycle%source - 1
          length      = nlyr * ntrac * pCycle%lengthP_ghost*nhc*npoints
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(buffer%receive(1,1,1,1,iptr),length,MPIreal_t, &
               source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,subname,': Error after call to MPI_Irecv: ',errorstring
          endif
       end do    ! icycle


       !==================================================
       !  Wait for all the receives to complete
       !==================================================

       call MPI_Waitall(nSendCycles,Srequest,status,ierr)
       call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          length             = pCycle%lengthP_ghost
          iptr            = pCycle%ptrP_ghost
          do i=0,length-1
             do i2=1,nhc
                do i1=1,npoints
                   buffer%buf(i1,i2,1:nlyr,1:ntrac,iptr+i) = buffer%receive(i1,i2,1:nlyr,1:ntrac,iptr+i)
                enddo
             enddo
          enddo
       end do   ! icycle


#endif
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER

  end subroutine ghost_exchangeV

  subroutine compute_ghost_corner_orientation(hybrid,elem,nets,nete)
!
!  this routine can NOT be called in a threaded region because then each thread
!  will have its on ghostbuffer.   initghostbufer3D() should detect this and abort.
!
  use kinds, only : real_kind
  use dimensions_mod, only: nelemd, np
  use parallel_mod, only : syncmp
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use edgetype_mod, only : ghostbuffer3D_t
  use edge_mod, only : ghostvpackfull, ghostvunpackfull, &
       initghostbuffer3D,freeghostbuffer3D
  use control_mod, only : north,south,east,west,neast, nwest, seast, swest

  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer :: nets,nete
  type (ghostBuffer3D_t)   :: ghostbuf_cv

  real (kind=real_kind) :: cin(2,2,1,nets:nete)  !CE: fvm tracer
  real (kind=real_kind) :: cout(-1:4,-1:4,1,nets:nete)  !CE: fvm tracer
  integer :: i,j,ie,kptr,np1,np2,nc,nc1,nc2,k,nlev
  logical :: fail,fail1,fail2
  real (kind=real_kind) :: tol=.1
  call syncmp(hybrid%par)
!   if (hybrid%par%masterproc) print *,'computing ghost cell corner orientations'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! first test on the Gauss Grid with same number of ghost cells:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nc=2   ! test using GLL interior points
  nc1=-1
  nc2=4

  nlev=1
  call initghostbuffer3D(ghostbuf_cv,nlev,nc)


  do ie=nets,nete
     cin(1,1,1,ie)=  elem(ie)%gdofp(1,1)
     cin(nc,nc,1,ie)=  elem(ie)%gdofp(np,np)
     cin(1,nc,1,ie)=   elem(ie)%gdofp(1,np)
     cin(nc,1,1,ie)=  elem(ie)%gdofp(np,1)
  enddo
  cout=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  run ghost exchange on c array to get corner orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     kptr=0
     call ghostVpackfull(ghostbuf_cv, cin(:,:,:,ie),1,nc,nc,nlev,kptr,elem(ie)%desc)
  end do
  call ghost_exchangeVfull(hybrid%par,hybrid%ithr,ghostbuf_cv)
  do ie=nets,nete
     kptr=0
     call ghostVunpackfull(ghostbuf_cv, cout(:,:,:,ie), nc1,nc2,nc,nlev, kptr, elem(ie)%desc)
  enddo

!       nc +--------+ 
!        ^ | nw  ne |    
!     j  | |        |
!        1 | sw  se |
!          +--------+
!           1 --> nc
!              i

! check SW corner
  do ie=nets,nete
     fail1=.false.
     fail2=.false.
     if ( elem(ie)%desc%putmapP_ghost(swest) /= -1) then
        if (abs(cout(nc1,1,1,ie)-cout(nc1,0,1,ie)) .gt. tol )  fail1=.true.
        if (abs(cout(1,nc1,1,ie)-cout(0,nc1,1,ie)).gt.tol) fail2=.true.
     endif
     if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange SW orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(swest)=.true.
        !print *,'reversion sw orientation ie',ie
        !print *,elem(ie)%desc%reverse(nwest),elem(ie)%desc%reverse(north),elem(ie)%desc%reverse(neast)
        !print *,elem(ie)%desc%reverse(west),' ',elem(ie)%desc%reverse(east)
        !print *,elem(ie)%desc%reverse(swest),elem(ie)%desc%reverse(south),elem(ie)%desc%reverse(seast)
     endif
  enddo
! check SE corner
  do ie=nets,nete
     fail1=.false.
     fail2=.false.
     if ( elem(ie)%desc%putmapP_ghost(seast) /= -1) then
        if (abs(cout(nc2,1,1,ie)-cout(nc2,0,1,ie)) .gt. tol )  fail1=.true.
        if (abs(cout(nc+1,nc1,1,ie)-cout(nc,nc1,1,ie)).gt.tol) fail2=.true.
     endif
     if (fail1 .neqv. fail2 ) call abortmp('ghost exchange SE orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(seast)=.true.
     endif
  enddo
! check NW corner
  do ie=nets,nete
     fail1=.false.
     fail2=.false.
     if ( elem(ie)%desc%putmapP_ghost(nwest) /= -1) then
        if (abs(cout(nc1,nc+1,1,ie)-cout(nc1,nc,1,ie)) .gt. tol )  fail1=.true.
        if (abs(cout(1,nc2,1,ie)-cout(0,nc2,1,ie)).gt.tol) fail2=.true.
     endif
     if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange NW orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(nwest)=.true.
     endif
  enddo
! check NE corner
  do ie=nets,nete
     fail1=.false.
     fail2=.false.
     if ( elem(ie)%desc%putmapP_ghost(neast) /= -1) then
        if (abs(cout(nc2,nc+1,1,ie)-cout(nc2,nc,1,ie)) .gt. tol )  fail1=.true.
        if (abs(cout(nc+1,nc2,1,ie)-cout(nc,nc2,1,ie)).gt.tol) fail2=.true.
     endif
     if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange NE orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(neast)=.true.
     endif
  enddo
  call freeghostbuffer3D(ghostbuf_cv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  end ghost exchange corner orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine


  subroutine sort_neighbor_buffer_mapping(par,elem,nets,nete)
!
!  gather global ID's of all neighbor elements.  Then create sorted (in global ID numbering)
!  mapping between edge buffer for each neighbor and a local map.  
!
!  this routine can NOT be called in a threaded region because then each thread
!  will have its on ghostbuffer.   initghostbufer3D() should detect this and abort.
!
!  also return num_neigh(ie) = number of neighbors (including onself) for element ie
!  
!
  use kinds, only : real_kind
  use dimensions_mod, only: nelemd, np, max_neigh_edges
  use parallel_mod, only : syncmp, parallel_t
  use element_mod, only : element_t
  use edgetype_mod, only : ghostbuffer3D_t
  use edge_mod, only : ghostvpack_unoriented, ghostvunpack_unoriented, &
       initghostbuffer3D,freeghostbuffer3D
  use control_mod, only : north,south,east,west,neast, nwest, seast, swest
  use coordinate_systems_mod, only: cartesian3D_t
  implicit none

  type (parallel_t)      , intent(in) :: par
  type (element_t)     , intent(inout), target :: elem(:)
  integer :: nets,nete
  type (ghostBuffer3D_t)   :: ghostbuf_cv

  real (kind=real_kind) :: cin(2,2,4,nets:nete)                    ! 1x1 element input data
  real (kind=real_kind) :: cout(2,2,4,max_neigh_edges+1,nets:nete)   ! 1x1 element output data
  real (kind=real_kind) :: u   (2,2,4)   
  integer :: i,j,ie,kptr,np1,np2,nc,k,nlev,patch_size,l,l2,sum1,sum2,m
  logical :: fail,fail1,fail2
  real (kind=real_kind) :: tol=.1


  if (par%masterproc) print *,'creating sorted ghost cell neigbor map...' 
  if (par%masterproc) print *,'checking ghost cell neighbor buffer sorting...' 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! first test on the Gauss Grid with same number of ghost cells:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nc=2
  nlev=4
  call initghostbuffer3D(ghostbuf_cv,nlev,nc)

  call syncmp(par)

  do ie=nets,nete
     cin(:,:,nlev,ie)=  elem(ie)%GlobalID
     k=0
     do i=1,nc
     do j=1,nc
        k=k+1
        cin(i,j,1,ie) = elem(ie)%corners3D(k)%x
        cin(i,j,2,ie) = elem(ie)%corners3D(k)%y
        cin(i,j,3,ie) = elem(ie)%corners3D(k)%z
     enddo
     enddo
  enddo
  cout=-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  run ghost exchange to get global ID of all neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     kptr=0
     call ghostVpack_unoriented(ghostbuf_cv, cin(:,:,:,ie),nc,nlev,kptr,elem(ie)%desc)
  end do

  ! check for array out of bouds overwriting 
  if (int(maxval(  cout(:,:,:,:,:))) /= -1 ) then
     call abortmp('ghost excchange unoriented failure ob1')
  endif
  call ghost_exchangeVfull(par,0,ghostbuf_cv)
  if (int(maxval(  cout(:,:,:,:,:))) /= -1 ) then
     call abortmp('ghost excchange unoriented failure ob2')
  endif

  do ie=nets,nete
     kptr=0

     call ghostVunpack_unoriented(ghostbuf_cv, cout(:,:,:,:,ie),nc,nlev, kptr, elem(ie)%desc, elem(ie)%GlobalId,cin(:,:,:,ie))

     ! check that we get the count of real neighbors correct
     patch_size=0

     do l=1,max_neigh_edges+1
        if (int(cout(1,1,nlev,l,ie)) /= -1 ) then
           patch_size = patch_size + 1
        endif
     enddo

     if (elem(ie)%desc%actual_neigh_edges+1 /= patch_size) then
        print *,'desc  actual_neigh_edges: ',elem(ie)%desc%actual_neigh_edges
        print *,'check patch_size: ',patch_size
        call abortmp( 'ghost exchange unoriented failure 1')
     endif

     ! check that all non-neighbors stayed -1
     do l=patch_size+1,max_neigh_edges+1
     if (int(cout(1,1,nlev,l,ie)) /= -1 ) then
        call abortmp( 'ghost exchange unoriented failure 2')
     endif
     enddo

     ! i am too lazy to check if all id's are identical since they are in
     ! different order.  check if there sum is identical
     sum1 = sum(int(cout(1,1,nlev,1:patch_size,ie)))
     sum2 = elem(ie)%globalID
     do l=1,max_neigh_edges
        if (elem(ie)%desc%globalID(l)>0) sum2 = sum2 + elem(ie)%desc%globalID(l)
     enddo
     if (sum1 /= sum2 ) then
        print *,int(cin(1,1,nlev,ie)),elem(ie)%desc%actual_neigh_edges,patch_size
        write(*,'(a,99i5)') 'ghost=',int(cout(1,1,nlev,1:patch_size,ie))
        write(*,'(a,99i5)') 'desc =',elem(ie)%desc%globalID(:)

        print *,'cout sum of all neighbor global ids:',sum1 
        print *,'desc sum of all neighbor global ids:',sum2  
        call abortmp( 'ghost exchange unoriented failure 3')        
     endif

     ALLOCATE(elem(ie)%desc%neigh_corners(4,patch_size))
     ! unpack corner data into array
     do l=1,patch_size
        k=0
        do i=1,nc
        do j=1,nc
           k=k+1
           elem(ie)%desc%neigh_corners(k,l)%x = cout(i,j,1,l,ie) 
           elem(ie)%desc%neigh_corners(k,l)%y = cout(i,j,2,l,ie) 
           elem(ie)%desc%neigh_corners(k,l)%z = cout(i,j,3,l,ie) 
        enddo
        enddo
     enddo
  enddo

  call freeghostbuffer3D(ghostbuf_cv)
  if (par%masterproc) print *,'passed.'
  end subroutine




end module bndry_mod
