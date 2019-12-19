#include <stdio.h>
#include <stdarg.h>
#include "ldm_alloc.h"
#include "dma_macros.h"
#include "edge_slave.h"
void edgevpackadvp1_parallel_(edgeAdvp1_param_t *global_param) {
  edgeAdvp1_param_t sl_param;
  int *sl_reverse, *sl_putmap, *sl_iwesn;
  double *sl_elem_v, *sl_elem_v1;
  double *sl_elem_v_west, *sl_elem_v_east, *sl_elem_v_south, *sl_elem_v_north;
  
  double *sl_elem_swest, *sl_elem_seast, *sl_elem_nwest, *sl_elem_neast;
  
  int edgeptr, offset;
  int i, j, k, ie, q, kptr, iptr;

  dma_init();
  if ( _MYID == 0 ) {
    bcast_get(global_param, &sl_param, sizeof(edgeAdvp1_param_t));
    dma_syn();
  }
  athread_syn(ARRAY_SCOPE, 0xFFFF);
  
  ldm_alloc_init();
  ldm_alloc_pointer(sl_reverse, 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
  ldm_alloc_pointer(sl_iwesn, 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
  ldm_alloc_pointer(sl_putmap, 4*sl_param.max_corner_elem*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
  
  if ( _MYID == 0 ) {
    bcast_get(sl_param.reverse_addr, sl_reverse, 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
    bcast_get(sl_param.iwesn_addr, sl_iwesn, 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
    bcast_get(sl_param.putmap_addr, sl_putmap, 4*sl_param.max_corner_elem*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE));
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

  ldm_alloc_pointer(sl_elem_v, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_v1, sl_param.np*sl_param.np*sizeof(double));

  ldm_alloc_pointer(sl_elem_v_west, sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_v_east, sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_v_south, sl_param.nlev*sl_param.np*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_v_north, sl_param.nlev*sl_param.np*sizeof(FTYPE));
  
  ldm_alloc_pointer(sl_elem_swest, sl_param.nlev*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_seast, sl_param.nlev*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_nwest, sl_param.nlev*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_neast, sl_param.nlev*sizeof(FTYPE));
  

 
  for ( ie=_ROW; ie<(sl_param.nete-sl_param.nets+1); ie+=8 ) {

    for (q=_COL; q<(sl_param.qsize+1); q+=8) {
      if ( q != sl_param.qsize) {
        pe_get(sl_param.addrs[0]+ie*sl_param.state_size+q*sl_param.nlev*sl_param.np*sl_param.np, sl_elem_v, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
	  dma_syn();
	}else{
        pe_get(sl_param.addrs[1]+ie*sl_param.derived_size, sl_elem_v, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
        pe_get(sl_param.addrs[2]+ie*sl_param.elem_size, sl_elem_v1, sl_param.np*sl_param.np*sizeof(FTYPE));
        dma_syn();

	  for ( k=0; k<sl_param.nlev; k++) {
	    for ( j=0; j<4; j++) {
	      for (i=0; i<4; i++) {
	        sl_elem_v[k*sl_param.np*sl_param.np+j*sl_param.np+i] = 
	             sl_elem_v1[j*sl_param.np+i]*sl_elem_v[k*sl_param.np*sl_param.np+j*sl_param.np+i];
	      }
	    }
	  }

	pe_put(sl_param.addrs[1]+ie*sl_param.derived_size, sl_elem_v, sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE));
	}

      kptr = sl_param.nlev*q;
	iptr = sl_param.np*kptr;

	for (k=0; k<sl_param.nlev; k++) {
	  for (i=0; i<sl_param.np; i++) {
	    sl_elem_v_east[k*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1];
	    sl_elem_v_south[k*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+i];
	    sl_elem_v_north[k*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i];
	    sl_elem_v_west[k*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+i*sl_param.np];
	  }
	}
    //EAST
      if ( sl_reverse[ie*4+1] ) {
        for ( k=0; k < sl_param.nlev; k++ ) {
          for ( i=0; i< sl_param.np; i++ ) {
            sl_elem_v_east[k*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-i-1)*sl_param.np+sl_param.np-1];
          }
        }
	}
      pe_put(sl_param.buf+iptr+sl_iwesn[ie*4+1], sl_elem_v_east, sl_param.nlev*sl_param.np*sizeof(FTYPE));
  //SOUTH
      if ( sl_reverse[ie*4+2] ) {
        for ( k=0; k < sl_param.nlev; k++ ) {
          for ( i=0; i< sl_param.np; i++ ) {
            sl_elem_v_south[k*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+sl_param.np-i-1];
          }
        }
	}
      pe_put(sl_param.buf+iptr+sl_iwesn[ie*4+2], sl_elem_v_south, sl_param.nlev*sl_param.np*sizeof(FTYPE));
    //NORTH
	if (sl_reverse[ie*4+3]) {
        for ( k=0; k<sl_param.nlev; k++ ) {
          for ( i=0; i< sl_param.np; i++ ) {
            sl_elem_v_north[k*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1-i];
          }
        }
	}
      pe_put(sl_param.buf+iptr+sl_iwesn[ie*4+3], sl_elem_v_north, sl_param.nlev*sl_param.np*sizeof(FTYPE));
    //WEST
      if ( sl_reverse[ie*4] ) {
        for ( k=0; k < sl_param.nlev; k++ ) {
          for ( i=0; i< sl_param.np; i++ ) {
            sl_elem_v_west[k*sl_param.np+i] = sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-i-1)*sl_param.np];
          }
        }
      }
      pe_put(sl_param.buf+iptr+sl_iwesn[ie*4], sl_elem_v_west, sl_param.nlev*sl_param.np*sizeof(FTYPE));
	dma_syn();

    // SWEST
      for ( j=0; j<sl_param.max_corner_elem; j++ ) {
        if ( sl_putmap[ie*4*sl_param.max_corner_elem+j] != -1) {
          edgeptr = sl_putmap[ie*4*sl_param.max_corner_elem+j];
          for (k=0; k<sl_param.nlev; k++){
            sl_elem_swest[k] = sl_elem_v[k*sl_param.np*sl_param.np];
          }
          pe_put(sl_param.buf+kptr+edgeptr, sl_elem_swest, sl_param.nlev*sizeof(FTYPE));
          dma_syn();
        }
      }

    // SEAST
      for ( j=sl_param.max_corner_elem; j<2*sl_param.max_corner_elem; j++ ) {
        if ( sl_putmap[ie*4*sl_param.max_corner_elem+j] != -1) {
          edgeptr = sl_putmap[ie*4*sl_param.max_corner_elem+j];
          for ( k =0; k<sl_param.nlev; k++ ) {
            sl_elem_seast[k] = sl_elem_v[k*sl_param.np*sl_param.np+sl_param.np-1];
          }
          pe_put(sl_param.buf+kptr+edgeptr, sl_elem_seast, sl_param.nlev*sizeof(FTYPE));
          dma_syn();
        }
      }
    // NEAST
      for ( j =3*sl_param.max_corner_elem; j < 4*sl_param.max_corner_elem; j++) {
        if ( sl_putmap[ie*4*sl_param.max_corner_elem+j] != -1 ) {
          edgeptr = sl_putmap[ie*4*sl_param.max_corner_elem+j];
	    for ( k=0; k<sl_param.nlev; k++ ) {
            sl_elem_neast[k] = sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1];
	    }
          pe_put(sl_param.buf+kptr+edgeptr, sl_elem_neast, sl_param.nlev*sizeof(FTYPE));
          dma_syn();
        }
      }
    // NWEST
      for ( j=2*sl_param.max_corner_elem; j<3*sl_param.max_corner_elem; j++ ) {
        if(sl_putmap[ie*4*sl_param.max_corner_elem+j] != -1 ) {
          edgeptr = sl_putmap[ie*4*sl_param.max_corner_elem+j];
	    for ( k=0; k<sl_param.nlev; k++ ){
            sl_elem_nwest[k] = sl_elem_v[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np];
	    }
          pe_put(sl_param.buf+kptr+edgeptr, sl_elem_nwest, sl_param.nlev*sizeof(FTYPE));
          dma_syn();
        }
      }
    }
  }
  ldm_dealloc_after(sl_reverse);
}


