#include <mpi.h>


/*
// mpi init function
int __wrap_MPI_Init(int *argc, char ***argv) {
  int Mpi_Init_Return = 0 ; 
  set_queue_malloc_mask_();
  Mpi_Init_Return = __real_MPI_Init(argc , argv) ;
  reset_queue_malloc_mask_();
  return Mpi_Init_Return ;
}
// blocking commucation function
int __wrap_MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  int Mpi_Send_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Send_Return = __real_MPI_Send(buf , count, datatype, dest, tag, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Send_Return ; 
}

int __wrap_MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm , MPI_Status *status) {
  int Mpi_Recv_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Recv_Return = __real_MPI_Recv(buf , count, datatype, source, tag, comm , status) ;
  reset_queue_malloc_mask_();
  return Mpi_Recv_Return ;
}
int __wrap_MPI_Sendrecv(void *sendbuf , int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  int Mpi_Sendrecv_Return  = 0 ;
  set_queue_malloc_mask_();
  Mpi_Sendrecv_Return = __real_MPI_Sendrecv(sendbuf , sendcount , sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status) ;
  reset_queue_malloc_mask_();
  return Mpi_Sendrecv_Return ;
}


int __wrap_MPI_Sendrecv_replace(void *buf , int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  int Mpi_Sendrecv_replace_Return  = 0 ;
  set_queue_malloc_mask_();
  Mpi_Sendrecv_replace_Return = __real_MPI_Sendrecv_replace(buf , count, datatype, dest, sendtag, source, recvtag, comm, status) ;
  reset_queue_malloc_mask_();
  return Mpi_Sendrecv_replace_Return ; 
}
// non_blocking communication function
int __wrap_MPI_Isend(void *buf , int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  int Mpi_Isend_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Isend_Return = __real_MPI_Isend(buf , count, datatype, dest, tag, comm, request) ;
  reset_queue_malloc_mask_();
  return Mpi_Isend_Return ; 
}

#ifdef NONBLOCAKING_BUFFERMODE
int __wrap_MPI_Ibsend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  int Mpi_Ibsend_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Ibsend_Return = __real_MPI_Ibsend(buf , count, datatype, dest, tag, comm, request) ;
  reset_queue_malloc_mask_();
  return Mpi_Isend_Return ;
}
#endif

#ifdef NONBLOCAKING_SYNCHRONOUSMODE
int __wrap_MPI_Issend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  int Mpi_Issend_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Issend_Return = __real_MPI_Issend(buf , count , datatype, dest, tag, comm, request) ;
  reset_queue_malloc_mask_();
  return Mpi_Issend_Return ;
}
#endif

#ifdef NONBLOCAKING_READYMODE 
int __wrap_MPI_Irsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  int Mpi_Irsend_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Irsend_Return = __real_MPI_Irsend(buf , count , datatype , dest , tag , comm , request) ;
  reset_queue_malloc_mask_();
  return Mpi_Irsend_Return ;
}
#endif

int __wrap_MPI_Irecv(void *buf , int count , MPI_Datatype datatype , int source , int tag , MPI_Comm comm , MPI_Request *request) {
  int Mpi_Irecv_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Irecv_Return = __real_MPI_Irecv(buf , count, datatype, source, tag, comm, request) ;
  reset_queue_malloc_mask_();
  return Mpi_Irecv_Return ;
}

int __wrap_MPI_Wait(MPI_Request *request , MPI_Status *status) {
  int Mpi_Wait_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Wait_Return = __real_MPI_Wait(request , status) ;  
  reset_queue_malloc_mask_();
  return Mpi_Wait_Return ;  
}

int __wrap_MPI_Waitany(int count , MPI_Request *Request_Array, int *index , MPI_Status *status) {
  int Mpi_Waitany_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Waitany_Return = __real_MPI_Waitany(count , Request_Array , index , status) ; 
  reset_queue_malloc_mask_();
  return Mpi_Waitany_Return ; 
}

int __wrap_MPI_Waitsome(int incount , MPI_Request *Request_Array , int *outcount , int *Indices_Array , MPI_Status *Status_Array) {
  int Mpi_Waitsome_Return  = 0 ;
  set_queue_malloc_mask_();
  Mpi_Waitsome_Return = __real_MPI_Waitsome(incount , Request_Array , outcount , Indices_Array , Status_Array) ;
  reset_queue_malloc_mask_();
  return Mpi_Waitsome_Return ;
}
int __wrap_MPI_Waitall(int count , MPI_Request *Request_Array , MPI_Status *Status_Array) {
  int Mpi_Waitall_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Waitall_Return = __real_MPI_Waitall(count , Request_Array , Status_Array) ;
  reset_queue_malloc_mask_();
  return Mpi_Waitall_Return ;
}

// collective communication function
int __wrap_MPI_Barrier(MPI_Comm comm) {
  int Mpi_Barrier_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Barrier_Return = __real_MPI_Barrier(comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Barrier_Return ;
}

int __wrap_MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm ){
  int Mpi_Bcast_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Bcast_Return = __real_MPI_Bcast(buffer , count , datatype , root , comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Bcast_Return ; 
}

int __wrap_MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,MPI_Comm comm) {
  int Mpi_Gather_Return  = 0 ;
  set_queue_malloc_mask_();
  Mpi_Gather_Return = __real_MPI_Gather(sendbuf , sendcount , sendtype , recvbuf , recvcount , recvtype , root , comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Gather_Return ;
}

int __wrap_MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  int Mpi_Gatherv_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Gatherv_Return = __real_MPI_Gatherv(sendbuf , sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Gatherv_Return ;
}

int __wrap_MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  int Mpi_Scatter_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Scatter_Return = __real_MPI_Scatter(sendbuf , sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Scatter_Return ;
}

int __wrap_MPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) { 
  int Mpi_Scatterv_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Scatterv_Return = __real_MPI_Scatterv(sendbuf , sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm) ; 
  reset_queue_malloc_mask_();
  return Mpi_Scatterv_Return ;
}

int __wrap_MPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  int Mpi_Allgather_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Allgather_Return = __real_MPI_Allgather(sendbuf , sendcount, sendtype, recvbuf, recvcount, recvtype, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Allgather_Return ; 
}

int __wrap_MPI_Allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm) {
  int Mpi_Allgatherv_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Allgatherv_Return = __real_MPI_Allgatherv(sendbuf , sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Allgatherv_Return ;
}

int __wrap_MPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  int Mpi_Alltoall_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Alltoall_Return = __real_MPI_Alltoall(sendbuf , sendcount, sendtype, recvbuf, recvcount, recvtype, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Alltoall_Return ;
}

int __wrap_MPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {
  int Mpi_Alltoallv_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Alltoallv_Return = __real_MPI_Alltoallv(sendbuf , sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Alltoallv_Return ;
}

int __wrap_MPI_Alltoallw(void *sendbuf, int sendcounts[], int sdispls[], MPI_Datatype sendtypes[], void *recvbuf, int recvcounts[], int rdispls[], MPI_Datatype recvtypes[], MPI_Comm comm) {
  int Mpi_Alltoallw_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Alltoallw_Return = __real_MPI_Alltoallw(sendbuf , sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Alltoallw_Return ;
}
int __wrap_MPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  int Mpi_Reduce_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Reduce_Return = __real_MPI_Reduce(sendbuf , recvbuf, count, datatype, op, root, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Reduce_Return ;
}

int __wrap_MPI_Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  int Mpi_Allreduce_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Allreduce_Return = __real_MPI_Allreduce(sendbuf , recvbuf, count, datatype, op, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Allreduce_Return ;
}

int __wrap_MPI_Reduce_scatter(void* sendbuf, void* recvbuf, int *recvcounts, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  int Mpi_Reduce_scatter_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Reduce_scatter_Return = __real_MPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Reduce_scatter_Return ;
}


int __wrap_MPI_Reduce_scatter_block(void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  int Mpi_Reduce_scatter_block_Return = 0 ;
  set_queue_malloc_mask_();
  Mpi_Reduce_scatter_block_Return = __real_MPI_Reduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op, comm) ;
  reset_queue_malloc_mask_();
  return Mpi_Reduce_scatter_block_Return ;
}
// mpi finish function
int __wrap_MPI_Finalize() {
  set_queue_malloc_mask_();
  int Mpi_Finalize_Return = 0 ;
  reset_queue_malloc_mask_();
} 
// */
