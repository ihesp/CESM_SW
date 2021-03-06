
module edgetype_mod 

  use kinds, only : int_kind, long_kind, log_kind, real_kind
  use coordinate_systems_mod, only : cartesian3D_t
  use gbarriertype_mod, only: gbarrier_t

  implicit none 
  private 
  save 

  integer, public :: initedgebuffer_callid = 0

  type, public :: rotation_t
     integer  :: nbr                                        ! nbr direction: north south east west
     integer  :: reverse                                    ! 0 = do not reverse order
     ! 1 = reverse order
     real (kind=real_kind), dimension(:,:,:), pointer :: R => null()  !  rotation matrix
  end type rotation_t

  type, public :: EdgeDescriptor_t
     integer(kind=int_kind)  :: use_rotation
     integer(kind=int_kind)  :: padding
     integer(kind=int_kind), pointer  :: putmapP(:) => null()
     integer(kind=int_kind), pointer  :: getmapP(:) => null()
     integer(kind=int_kind), pointer  :: putmapP_ghost(:) => null()
     integer(kind=int_kind), pointer  :: getmapP_ghost(:) => null()
     integer(kind=int_kind), pointer  :: putmapS(:) => null()
     integer(kind=int_kind), pointer  :: getmapS(:) => null()
     integer(kind=int_kind), pointer  :: globalID(:) => null()
     integer(kind=int_kind), pointer  :: loc2buf(:) => null()
     type (cartesian3D_t)  , pointer  :: neigh_corners(:,:) => null()
     integer                          :: actual_neigh_edges
     logical(kind=log_kind), pointer  :: reverse(:) => null()
     type (rotation_t), dimension(:), pointer :: rot => null() !  Identifies list of edges
     !  that must be rotated, and how
  end type EdgeDescriptor_t

  type, public :: EdgeBuffer_t
     real (kind=real_kind), dimension(:), allocatable :: buf
     real (kind=real_kind), dimension(:), allocatable :: receive
     integer(kind=int_kind), pointer :: putmap(:,:) => null()
     integer(kind=int_kind), pointer :: getmap(:,:) => null()
     logical(kind=log_kind), pointer :: reverse(:,:) => null()
     integer(kind=int_kind), pointer :: moveLength(:) => null()
     integer(kind=int_kind), pointer :: movePtr(:) => null()

     integer(kind=int_kind), pointer :: rcountsFull(:) => null()
     integer(kind=int_kind), pointer :: scountsFull(:) => null()
     integer(kind=int_kind), pointer :: sdisplsFull(:) => null()
     integer(kind=int_kind), pointer :: rdisplsFull(:) => null()

     integer(kind=int_kind), pointer :: rcountsInter(:) => null()
     integer(kind=int_kind), pointer :: scountsInter(:) => null()
     integer(kind=int_kind), pointer :: sdisplsInter(:) => null()
     integer(kind=int_kind), pointer :: rdisplsInter(:) => null()

     integer(kind=int_kind), pointer :: rcountsIntra(:) => null()
     integer(kind=int_kind), pointer :: scountsIntra(:) => null()
     integer(kind=int_kind), pointer :: sdisplsIntra(:) => null()
     integer(kind=int_kind), pointer :: rdisplsIntra(:) => null()

     integer(kind=int_kind), pointer :: getDisplsFull(:) => null()
     integer(kind=int_kind), pointer :: putDisplsFull(:) => null()

     integer(kind=int_kind), dimension(:), allocatable :: Rrequest,Srequest
     integer(kind=int_kind), dimension(:,:), allocatable :: status
     type (gbarrier_t) :: gbarrier
     integer :: nlyr ! Number of layers
     integer :: nbuf ! total size of message passing buffer, includes vertical levels
     integer :: nInter, nIntra
     integer :: id
     integer :: bndry_type
     integer :: tag
     integer :: win
     integer(kind=long_kind) :: winSize
  end type EdgeBuffer_t

  type, public :: LongEdgeBuffer_t
     integer :: nlyr
     integer :: nbuf
     integer (kind=int_kind), dimension(:,:), pointer :: buf => null()
     integer (kind=int_kind), dimension(:,:), pointer :: receive => null()
  end type LongEdgeBuffer_t


  type, public :: GhostBufferTR_t
     real (kind=real_kind), dimension(:,:,:,:,:), pointer :: buf => null()
     real (kind=real_kind), dimension(:,:,:,:,:), pointer :: receive => null()
     integer :: nlyr ! Number of layers
     integer :: nbuf ! size of the horizontal dimension of the buffers.
  end type GhostBufferTR_t

  type, public :: GhostBuffer3D_t
     real (kind=real_kind), dimension(:,:,:,:), pointer :: buf => null()
     real (kind=real_kind), dimension(:,:,:,:), pointer :: receive => null()
     integer :: nlyr ! Number of layers
     integer :: nhc  ! Number of layers of ghost cells
     integer :: np   ! Number of points in a cell
     integer :: nbuf ! size of the horizontal dimension of the buffers.
     integer :: elem_size ! size of 2D array (first two dimensions of buf())
  end type GhostBuffer3D_t

 
end module edgetype_mod
