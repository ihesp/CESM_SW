module static_part_mod
  use parallel_mod
  use gridgraph_mod
contains
  subroutine genstaticpart(GridEdge, GridVertex, par)
    implicit none
#include <mpif.h>
    type (parallel_t),   intent(in)    :: par
    type (GridVertex_t), intent(inout) :: GridVertex(:)
    type (GridEdge_t),   intent(inout) :: GridEdge(:)
    integer(kind=4), allocatable  :: procs(:)
    integer               :: nelem,nelem_edge,nelemd
    integer               :: head_part,tail_part
    integer               :: j,k,tmp1,id,s1,extra, ierr
    
    !print *, 'npart=', npart
    nelem      = SIZE(GridVertex(:))
    allocate(procs(nelem))
    if (par%masterproc) then
       open(1416, file="static-assignment.txt", form='formatted')
       read(1416, *), procs
       close(1416)
    endif
    call mpi_bcast(procs, nelem, MPI_INT, par%root, par%comm, ierr)
    do k=1,nelem
       id = GridVertex(k)%Number
       if (id <= 0 .or. id > nelem) print *, "Impossible id"
       !print *, "k=", k, "id=", id, "npart=", npart, "proc=", mod(id, npart) + 1
       GridVertex(k)%processor_number = procs(id) + 1
    end do
    deallocate(procs)
  end subroutine genstaticpart
end module static_part_mod
