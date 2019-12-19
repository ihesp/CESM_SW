#include <stdio.h>
#include <stdarg.h>
#include "ldm_alloc.h"
#include "dma_macros.h"
#include "edge_slave.h"

void edgevunpackq_parallel_(edgeq_param_t *global_param) {
  edgeq_param_t sl_param;

  int *sl_iwesn, *sl_getmap;
  double *sl_elem_v; 
  double *sl_receive_east, *sl_receive_west, *sl_receive_south, *sl_receive_north;
  int kptr, iptr;
  int i,j,k,ie,q;

  dma_init();
  if (_MYID == 0 ) {
    bcast_get(global_param, &sl_param, sizeof(edgeq_param_t));
    dma_syn();
  }
  athread_syn(ARRAY_SCOPE, 0xFFFF);
  
  int start_qsize;
  int block_qsize = (sl_param.qsize+1+7) >> 3; //25+7/8=4
  int remainder_qsize = (sl_param.qsize+1) & 0x7; //remainder_qsize=25%8
  if ( _COL < remainder_qsize ) {
    start_qsize = _COL * block_qsize;
  }else{
    start_qsize = remainder_qsize*block_qsize + (_COL-remainder_qsize)*((sl_param.qsize+1) >> 3);
    block_qsize = (sl_param.qsize+1) >> 3;
  }

  ldm_alloc_init();

  ldm_alloc_pointer(sl_getmap, 4*sl_param.max_corner_elem*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
  ldm_alloc_pointer(sl_iwesn, 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));

  ldm_alloc_pointer(sl_elem_v, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));

  if ( _MYID == 0 ) {
    bcast_get(sl_param.getmap_addr, sl_getmap, 4*sl_param.max_corner_elem*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
    bcast_get(sl_param.iwesn_addr, sl_iwesn, 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
    dma_syn();
  }
  athread_syn(ARRAY_SCOPE, 0xFFFF);

  ldm_alloc_pointer(sl_receive_west, sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_receive_east, sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_receive_south, sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_receive_north, sl_param.nlev*sl_param.np*sizeof(FTYPE));

  for (ie=_ROW; ie<(sl_param.nete-sl_param.nets+1); ie+=8) {
    for (q=_COL; q<(sl_param.qsize); q+=8) {
      kptr = q*sl_param.nlev;
      iptr = kptr*sl_param.np;
      pe_get(sl_param.receive+iptr+sl_iwesn[ie*4], sl_receive_west, sl_param.nlev*sl_param.np*sizeof(FTYPE));
      pe_get(sl_param.receive+iptr+sl_iwesn[ie*4+1], sl_receive_east, sl_param.nlev*sl_param.np*sizeof(FTYPE));
      pe_get(sl_param.receive+iptr+sl_iwesn[ie*4+2], sl_receive_south, sl_param.nlev*sl_param.np*sizeof(FTYPE));
      pe_get(sl_param.receive+iptr+sl_iwesn[ie*4+3], sl_receive_north, sl_param.nlev*sl_param.np*sizeof(FTYPE));
	pe_get(sl_param.addrs[0]+ie*sl_param.qsize*sl_param.nlev*sl_param.np*sl_param.np+q*sl_param.nlev*sl_param.np*sl_param.np
	    , sl_elem_v, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
	dma_syn();

      for (k=0; k<sl_param.nlev; k++) {
        for ( i=0; i<sl_param.np; i++) {
	  //EAST
          sl_elem_v[k*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1] =  sl_elem_v[k*sl_param.np*sl_param.np
	        +i*sl_param.np+sl_param.np-1]+sl_receive_east[k*sl_param.np+i];
	  //SOUTH
          sl_elem_v[k*sl_param.np*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+i] + sl_receive_south[k*sl_param.np+i];
	  //NORTH
          sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i] = 
	        sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i] + sl_receive_north[k*sl_param.np+i];
	  //WEST
	    sl_elem_v[k*sl_param.np*sl_param.np+i*sl_param.np] = sl_elem_v[k*sl_param.np*sl_param.np+i*sl_param.np]
	        + sl_receive_west[k*sl_param.np+i];
	  }
      }

    //SWEST
      for (i=0; i<sl_param.max_corner_elem; i++){
	  if (sl_getmap[ie*4*sl_param.max_corner_elem+i] != -1) {
	    pe_get(sl_param.receive+kptr+sl_getmap[ie*4*sl_param.max_corner_elem+i], sl_receive_west, sl_param.nlev*sizeof(FTYPE));
	    dma_syn();
	    for (k =0; k<sl_param.nlev; k++){
	      sl_elem_v[k*sl_param.np*sl_param.np] = sl_elem_v[k*sl_param.np*sl_param.np] + sl_receive_west[k];
	    }
	  }
	}
    
    //SEAST
      for (i=sl_param.max_corner_elem; i<2*sl_param.max_corner_elem; i++) {
	  if ( sl_getmap[ie*4*sl_param.max_corner_elem+i] != -1) {
	    pe_get(sl_param.receive+kptr+sl_getmap[ie*4*sl_param.max_corner_elem+i], sl_receive_west, sl_param.nlev*sizeof(FTYPE));
	    dma_syn();
	    for (k=0; k<sl_param.nlev; k++) {
	      sl_elem_v[k*sl_param.np*sl_param.np+sl_param.np-1] = sl_elem_v[k*sl_param.np*sl_param.np+sl_param.np-1] + sl_receive_west[k];
	    }
	  }
	}
    //NEAST
      for (i=3*sl_param.max_corner_elem; i<4*sl_param.max_corner_elem; i++) {
	  if ( sl_getmap[ie*4*sl_param.max_corner_elem+i] != -1) {
	    pe_get(sl_param.receive+kptr+sl_getmap[ie*4*sl_param.max_corner_elem+i], sl_receive_west, sl_param.nlev*sizeof(FTYPE));
	    dma_syn();
	    for (k=0; k<sl_param.nlev; k++) {
	      sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1] = 
		    sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1] + sl_receive_west[k];
	    }
	  }
	}
    //NWEST
      for (i=2*sl_param.max_corner_elem; i<3*sl_param.max_corner_elem; i++) {
	  if ( sl_getmap[ie*4*sl_param.max_corner_elem+i] != -1) {
	    pe_get(sl_param.receive+kptr+sl_getmap[ie*4*sl_param.max_corner_elem+i], sl_receive_west, sl_param.nlev*sizeof(FTYPE));
	    dma_syn();
	    for (k=0; k<sl_param.nlev; k++) {
	      sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np] = 
		    sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np] + sl_receive_west[k];
	    }
	  }
      }  
   
    pe_put(sl_param.addrs[0]+ie*sl_param.qsize*sl_param.nlev*sl_param.np*sl_param.np+q*sl_param.nlev*sl_param.np*sl_param.np
        , sl_elem_v, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));

    dma_syn();
    }
//    athread_syn(ROW_SCOPE, 0x00FF);
  } 
  ldm_dealloc_after(sl_getmap);
}
