#include <stdarg.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#include "edge_slave.h"

void edgevunpack3p1_parallel_(edge3p1_param_t *global_param) {
  edge3p1_param_t sl_param;
  int *sl_iwesn, *sl_getmap;
  double *sl_elem_v1, *sl_elem_v2, *sl_elem_v3; 
  double *sl_receive_east, *sl_receive_west, *sl_receive_south, *sl_receive_north;

  int i,j,k,ie;

  dma_init();
  if (_MYID == 0 ) {
    bcast_get(global_param, &sl_param, sizeof(edge3p1_param_t));
    dma_syn();
  }
  athread_syn(ARRAY_SCOPE, 0xFFFF);
  
  ldm_alloc_init();

  ldm_alloc_pointer(sl_getmap, 4*sl_param.max_corner_elem*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
  ldm_alloc_pointer(sl_iwesn, 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));

  ldm_alloc_pointer(sl_elem_v1, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));  
  ldm_alloc_pointer(sl_elem_v2, 2*sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
  if (sl_param.rsplit > 0){
    ldm_alloc_pointer(sl_elem_v3, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
  }

  if ( _MYID == 0 ) {
    bcast_get(sl_param.getmap_addr, sl_getmap, 4*sl_param.max_corner_elem*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
    bcast_get(sl_param.iwesn_addr, sl_iwesn, 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
    dma_syn();
  }
  athread_syn(ARRAY_SCOPE, 0xFFFF);

  int start_nlev, start_2nlev;
  int block_nlev = (sl_param.nlev+7) >> 3; //block_nlev=37/8=4
  int block_2nlev = (2*sl_param.nlev+7) >> 3; //block_2nlev=67/8=8
  int remainder_nlev = sl_param.nlev & 0x7; //remainder_nlev=6
  int remainder_2nlev = (2*sl_param.nlev) & 0x7; //remainder_2nlev=4
  if ( _COL < remainder_nlev ) {
    start_nlev = _COL * block_nlev;
  }else{
    start_nlev = remainder_nlev*block_nlev + (_COL-remainder_nlev)*(sl_param.nlev >> 3);
    block_nlev = sl_param.nlev >> 3;
  }
  if ( _COL < remainder_2nlev ) {
    start_2nlev = _COL * block_2nlev;
  }else{
    start_2nlev =  remainder_2nlev*block_2nlev + (_COL-remainder_2nlev)*((2*sl_param.nlev) >> 3);
    block_2nlev = (2*sl_param.nlev) >> 3;
  }

  ldm_alloc_pointer(sl_receive_west, 4*sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_receive_east, 4*sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_receive_south, 4*sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_receive_north, 4*sl_param.nlev*sl_param.np*sizeof(FTYPE));

  for ( ie=_ROW; ie<(sl_param.nete-sl_param.nets+1); ie+=8 ) {
    if ( _COL == 0 ) {
      bcast_set_mask(bcast_get, (1<<_ROW));

      bcast_get(sl_param.addrs[0]+ie*sl_param.offset, sl_elem_v1
          , sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
      bcast_get(sl_param.addrs[1]+ie*sl_param.offset, sl_elem_v2
          , 2*sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
	if (sl_param.rsplit > 0) {
        bcast_get(sl_param.addrs[2]+ie*sl_param.offset, sl_elem_v3
            , sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
	}
    
      bcast_get(sl_param.receive+sl_param.np+sl_iwesn[ie*4], sl_receive_west, 4*sl_param.nlev*sl_param.np*sizeof(FTYPE));
      bcast_get(sl_param.receive+sl_param.np+sl_iwesn[ie*4+1], sl_receive_east, 4*sl_param.nlev*sl_param.np*sizeof(FTYPE));
      bcast_get(sl_param.receive+sl_param.np+sl_iwesn[ie*4+2], sl_receive_south, 4*sl_param.nlev*sl_param.np*sizeof(FTYPE));
      bcast_get(sl_param.receive+sl_param.np+sl_iwesn[ie*4+3], sl_receive_north, 4*sl_param.nlev*sl_param.np*sizeof(FTYPE));

      dma_syn();
    }
    athread_syn(ROW_SCOPE, 0x00FF);

    for (k=0; k<block_2nlev; k++) {
      for ( i=0; i<sl_param.np; i++) {
	  //EAST
        sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1] =  sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1]+sl_receive_east[(sl_param.nlev+k+start_2nlev)*sl_param.np+i];
	  //SOUTH
        sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+i] = sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+i] + sl_receive_south[(sl_param.nlev+k+start_2nlev)*sl_param.np+i];
	  //NORTH
        sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i] =   sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i] + sl_receive_north[(sl_param.nlev+k+start_2nlev)*sl_param.np+i];
	  //WEST
	  sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+i*sl_param.np] = sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+i*sl_param.np] + sl_receive_west[(sl_param.nlev+k+start_2nlev)*sl_param.np+i];
	}
    }

    for (k=0; k<block_nlev; k++) {
      for (i=0; i<sl_param.np; i++) {
	  //EAST
	  sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1] = 
	      sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1]+ sl_receive_east[(k+start_nlev)*sl_param.np+i];
	  //SOUTH
	  sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i] = sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i]
            + sl_receive_south[(k+start_nlev)*sl_param.np+i];
	  //NORTH
	  sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i] = 
	      sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i] + sl_receive_north[(k+start_nlev)*sl_param.np+i];
	  //WEST
	  sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np] = sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np]
            + sl_receive_west[(k+start_nlev)*sl_param.np+i];
	}
    }
    if ( sl_param.rsplit > 0 ) {
      for (k=0; k<block_nlev; k++) {
        for (i=0; i<sl_param.np; i++) {
	  //EAST
	    sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1] = 
	        sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1]+ sl_receive_east[(3*sl_param.nlev+k+start_nlev)*sl_param.np+i];
	  //SOUTH
	    sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+i] = sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+i]
	        + sl_receive_south[(3*sl_param.nlev+k+start_nlev)*sl_param.np+i];
	  //NORTH
	    sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i] = 
	        sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i] + sl_receive_north[(3*sl_param.nlev+k+start_nlev)*sl_param.np+i];
	  //WEST
	    sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np] = sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np]
              + sl_receive_west[(3*sl_param.nlev+k+start_nlev)*sl_param.np+i];
	  }
      }
    }
    //SWEST
      for (i=0; i<sl_param.max_corner_elem; i++){
	  if (sl_getmap[ie*4*sl_param.max_corner_elem+i] != -1) {
	    pe_get(sl_param.receive+1+sl_getmap[ie*4*sl_param.max_corner_elem+i], sl_receive_west, 4*sl_param.nlev*sizeof(FTYPE));
	    dma_syn();
	    for (k =0; k<block_nlev; k++){
	      sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np] = sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np] + sl_receive_west[start_nlev+k];
	    }
	    for (k =0; k<block_2nlev; k++){
	      sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np] = sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np] + sl_receive_west[sl_param.nlev+k+start_2nlev];
	    }
	    if (sl_param.rsplit > 0) {
	      for (k =0; k<block_nlev; k++){
	        sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np] = sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np] + sl_receive_west[3*sl_param.nlev+start_nlev+k];
	      }
	    }
	  }
	}
    
    //SEAST
      for (i=sl_param.max_corner_elem; i<2*sl_param.max_corner_elem; i++) {
	  if ( sl_getmap[ie*4*sl_param.max_corner_elem+i] != -1) {
	    pe_get(sl_param.receive+1+sl_getmap[ie*4*sl_param.max_corner_elem+i], sl_receive_west, 4*sl_param.nlev*sizeof(FTYPE));
	    dma_syn();
	    for (k=0; k<block_nlev; k++) {
	      sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+sl_param.np-1] = sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+sl_param.np-1] + sl_receive_west[start_nlev+k];
	    }
	    for (k=0; k<block_2nlev; k++) {
	      sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+sl_param.np-1] = sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+sl_param.np-1] + sl_receive_west[sl_param.nlev+start_2nlev+k];
	    }
	    if ( sl_param.rsplit > 0) {
	      for (k=0; k<block_nlev; k++) {
	        sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+sl_param.np-1] = sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+sl_param.np-1] + sl_receive_west[3*sl_param.nlev+start_nlev+k];
	      }
	    }
	  }
	}
    //NEAST
    for (i=3*sl_param.max_corner_elem; i<4*sl_param.max_corner_elem; i++) {
	  if ( sl_getmap[ie*4*sl_param.max_corner_elem+i] != -1) {
	    pe_get(sl_param.receive+1+sl_getmap[ie*4*sl_param.max_corner_elem+i], sl_receive_west, 4*sl_param.nlev*sizeof(FTYPE));
	    dma_syn();
	    for (k=0; k<block_nlev; k++) {
	      sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1] = 
		    sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1] + sl_receive_west[k+start_nlev];
	    }
	    for (k=0; k<block_2nlev; k++) {
	      sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1] = 
		    sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1] + sl_receive_west[sl_param.nlev+k+start_2nlev];
	    }
	    if ( sl_param.rsplit > 0 ) {
	      for (k=0; k<block_nlev; k++) {
	        sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1] =  
		      sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1] + sl_receive_west[3*sl_param.nlev+start_nlev+k];
	      }
	    }
	  }
	}
    //NWEST
    for (i=2*sl_param.max_corner_elem; i<3*sl_param.max_corner_elem; i++) {
	  if ( sl_getmap[ie*4*sl_param.max_corner_elem+i] != -1) {
	    pe_get(sl_param.receive+1+sl_getmap[ie*4*sl_param.max_corner_elem+i], sl_receive_west, 4*sl_param.nlev*sizeof(FTYPE));
	    dma_syn();
	    for (k=0; k<block_nlev; k++) {
	      sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np] = 
		    sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np] + sl_receive_west[k+start_nlev];
	    }
	    for (k=0; k<block_2nlev; k++) {
	      sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np] = 
		    sl_elem_v2[(k+start_2nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np] + sl_receive_west[sl_param.nlev+k+start_2nlev];
	    }
	    if ( sl_param.rsplit > 0) {
	      for (k=0; k<block_nlev; k++) {
	        sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np] =
		      sl_elem_v3[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np] + sl_receive_west[3*sl_param.nlev+k+start_nlev];
	      }
	    }
	  }
      }  

    pe_put(sl_param.addrs[0]+ie*sl_param.offset+start_nlev*sl_param.np*sl_param.np
        , sl_elem_v1+start_nlev*sl_param.np*sl_param.np, block_nlev*sl_param.np*sl_param.np*sizeof(FTYPE) );
    pe_put(sl_param.addrs[1]+ie*sl_param.offset+start_2nlev*sl_param.np*sl_param.np
        , sl_elem_v2+start_2nlev*sl_param.np*sl_param.np, block_2nlev*sl_param.np*sl_param.np*sizeof(FTYPE) );
    if ( sl_param.rsplit > 0) {
      pe_put(sl_param.addrs[2]+ie*sl_param.offset+start_nlev*sl_param.np*sl_param.np
          , sl_elem_v3+start_nlev*sl_param.np*sl_param.np, block_nlev*sl_param.np*sl_param.np*sizeof(FTYPE) );
    }
    dma_syn();
    athread_syn(ROW_SCOPE, 0x00FF);
  }
  ldm_dealloc_after(sl_getmap);
}
