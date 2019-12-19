#include <stdio.h>
#include <stdarg.h>
#include "ldm_alloc.h"
#include "dma_macros.h"
#include "edge_slave.h"

void edgevpack3p1_parallel_(edge3p1_param_t *global_param) {
  edge3p1_param_t sl_param;
  int *sl_putmap, *sl_iwesn, *sl_reverse;
  double *sl_elem_v1, *sl_elem_v2, *sl_elem_v3;
  double *sl_elem_v1_west, *sl_elem_v1_east, *sl_elem_v1_south, *sl_elem_v1_north;
  double *sl_elem_v2_west, *sl_elem_v2_east, *sl_elem_v2_south, *sl_elem_v2_north;
  double *sl_elem_v3_west, *sl_elem_v3_east, *sl_elem_v3_south, *sl_elem_v3_north;
  double *sl_elem_swest, *sl_elem_seast, *sl_elem_nwest, *sl_elem_neast;
  int edgeptr, kptr, vptr, v1_ptr, v2_ptr, v3_ptr, ieptr;
  int blockbytes, block2bytes, netbytes, vbytes;
  int i,j,k,ie;


  dma_init();
  if (_MYID == 0){
    bcast_get(global_param, &sl_param, sizeof(edge3p1_param_t));
    dma_syn();
  }
  athread_syn(ARRAY_SCOPE, 0xFFFF);
  
  ldm_alloc_init();
  netbytes = 4*(sl_param.nete-sl_param.nets+1)*sizeof(ITYPE);
  ldm_alloc_pointer(sl_putmap, sl_param.max_corner_elem*netbytes);
  ldm_alloc_pointer(sl_iwesn, netbytes);
  ldm_alloc_pointer(sl_reverse, netbytes);
  if (_MYID == 0) {
    bcast_get(sl_param.putmap_addr, sl_putmap, sl_param.max_corner_elem*netbytes);
    bcast_get(sl_param.iwesn_addr, sl_iwesn, netbytes);
    bcast_get(sl_param.reverse_addr, sl_reverse, netbytes);
    dma_syn();
  }
  athread_syn(ARRAY_SCOPE, 0xFFFF);
  
  vbytes = sl_param.nlev*sl_param.np*sl_param.np*sizeof(FTYPE);
  ldm_alloc_pointer(sl_elem_v1, vbytes);
  ldm_alloc_pointer(sl_elem_v2, 2*vbytes);
  ldm_alloc_pointer(sl_elem_v3, vbytes);
  
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
  
  blockbytes = block_nlev*sl_param.np*sizeof(FTYPE);
  ldm_alloc_pointer(sl_elem_v1_west, blockbytes);
  ldm_alloc_pointer(sl_elem_v1_east, blockbytes);
  ldm_alloc_pointer(sl_elem_v1_south, blockbytes);
  ldm_alloc_pointer(sl_elem_v1_north, blockbytes);
  
  ldm_alloc_pointer(sl_elem_v3_west, blockbytes);
  ldm_alloc_pointer(sl_elem_v3_east, blockbytes);
  ldm_alloc_pointer(sl_elem_v3_south, blockbytes);
  ldm_alloc_pointer(sl_elem_v3_north, blockbytes);  
   
  block2bytes = block_2nlev*sl_param.np*sizeof(FTYPE);
  ldm_alloc_pointer(sl_elem_v2_west, block2bytes);
  ldm_alloc_pointer(sl_elem_v2_east, block2bytes);
  ldm_alloc_pointer(sl_elem_v2_south, block2bytes);
  ldm_alloc_pointer(sl_elem_v2_north, block2bytes);
  
  ldm_alloc_pointer(sl_elem_swest, 4*sl_param.nlev*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_seast, 4*sl_param.nlev*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_nwest, 4*sl_param.nlev*sizeof(FTYPE));
  ldm_alloc_pointer(sl_elem_neast, 4*sl_param.nlev*sizeof(FTYPE));
  
  v1_ptr = (1+start_nlev)*sl_param.np;
  v2_ptr = (start_2nlev+sl_param.nlev+1)*sl_param.np;
  v3_ptr = (start_nlev+3*sl_param.nlev+1)*sl_param.np;
  for ( ie=_ROW; ie<(sl_param.nete-sl_param.nets+1); ie+=8) {

    ieptr = ie*sl_param.offset;
    if ( _COL == 0) {
      bcast_set_mask(bcast_get, (1<<_ROW));
	bcast_get(sl_param.addrs[0]+ieptr, sl_elem_v1,  vbytes);
	bcast_get(sl_param.addrs[1]+ieptr, sl_elem_v2, 2*vbytes);
	if ( sl_param.rsplit > 0) {
        bcast_get(sl_param.addrs[2]+ieptr, sl_elem_v3, vbytes);
	}
	dma_syn();
	if (_MYID == 0){
	  
	}
    }
    athread_syn(ROW_SCOPE, 0x00FF);
    
    
    for ( k=0; k<block_2nlev; k++ ) {
      kptr = k*sl_param.np;
	vptr = (k+start_2nlev)*sl_param.np*sl_param.np;
      for ( i=0; i<sl_param.np; i++ ) {
        //WEST
        sl_elem_v2_west[kptr+i] = sl_elem_v2[vptr+i*sl_param.np];
	  //EAST
        sl_elem_v2_east[kptr+i] = sl_elem_v2[vptr+i*sl_param.np+sl_param.np-1];
	  //SOUTH
        sl_elem_v2_south[kptr+i] = sl_elem_v2[vptr+i];
	  //NORTH
        sl_elem_v2_north[kptr+i] = sl_elem_v2[vptr+(sl_param.np-1)*sl_param.np+i];	
      }
    }

  
    for ( k=0; k < block_nlev; k++ ) {
      for ( i=0; i< sl_param.np; i++ ) {
	  //WEST
        sl_elem_v1_west[k*sl_param.np+i] = sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np];
	  //EAST
        sl_elem_v1_east[k*sl_param.np+i] = sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i*sl_param.np+sl_param.np-1];
	  //SOUTH
        sl_elem_v1_south[k*sl_param.np+i] = sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+i];
	  //NORTH
        sl_elem_v1_north[k*sl_param.np+i] = sl_elem_v1[(k+start_nlev)*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+i];
      }
    }

    //WEST
    if ( sl_reverse[ie*4] ) {
      for ( k=0; k<block_2nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_2nlev)*sl_param.np*sl_param.np;
        for ( i=0; i<sl_param.np; i++) {
          sl_elem_v2_west[kptr+i] = sl_elem_v2[vptr+(sl_param.np-i-1)*sl_param.np];	
        }
      }
      for ( k=0; k < block_nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_nlev)*sl_param.np*sl_param.np;
        for ( i=0; i< sl_param.np; i++ ) {
          sl_elem_v1_west[kptr+i] = sl_elem_v1[vptr+(sl_param.np-i-1)*sl_param.np];
        }
      }
    }

    pe_put(sl_param.buf+(1+start_nlev)*sl_param.np+sl_iwesn[ie*4], sl_elem_v1_west, block_nlev*sl_param.np*sizeof(FTYPE));
    pe_put(sl_param.buf+v2_ptr+sl_iwesn[ie*4], sl_elem_v2_west, block2bytes);

    //EAST
    if ( sl_reverse[ie*4+1] ) {
      for ( k=0; k<block_2nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_2nlev)*sl_param.np*sl_param.np;
        for ( i=0; i< sl_param.np; i++ ) {
          sl_elem_v2_east[kptr+i] = sl_elem_v2[vptr+(sl_param.np-i-1)*sl_param.np+sl_param.np-1];
        }
      }
      for ( k=0; k < block_nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_nlev)*sl_param.np*sl_param.np;
        for ( i=0; i< sl_param.np; i++ ) {
          sl_elem_v1_east[kptr+i] = sl_elem_v1[vptr+(sl_param.np-i-1)*sl_param.np+sl_param.np-1];
        }
      }
    }

    pe_put(sl_param.buf+v1_ptr+sl_iwesn[ie*4+1], sl_elem_v1_east, blockbytes);
    pe_put(sl_param.buf+v2_ptr+sl_iwesn[ie*4+1], sl_elem_v2_east, block2bytes);
    
  //SOUTH
    if ( sl_reverse[ie*4+2] ) {
      for ( k=0; k<block_2nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_2nlev)*sl_param.np*sl_param.np;
        for (i=0; i<sl_param.np; i++){
          sl_elem_v2_south[kptr+i] = sl_elem_v2[vptr+sl_param.np-i-1];
        }
      }
      for ( k=0; k < block_nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_nlev)*sl_param.np*sl_param.np;
        for ( i=0; i< sl_param.np; i++ ) {
          sl_elem_v1_south[kptr+i] = sl_elem_v1[vptr+sl_param.np-i-1];
        }
      }
    }
    
    pe_put(sl_param.buf+v1_ptr+sl_iwesn[ie*4+2], sl_elem_v1_south, blockbytes);
    pe_put(sl_param.buf+v2_ptr+sl_iwesn[ie*4+2], sl_elem_v2_south, block2bytes);
    
    //NORTH
    if (sl_reverse[ie*4+3]) {
      for ( k=0; k<block_2nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_2nlev)*sl_param.np*sl_param.np;
        for ( i=0; i< sl_param.np; i++ ) {
          sl_elem_v2_north[kptr+i] = sl_elem_v2[vptr+(sl_param.np-1)*sl_param.np+sl_param.np-1-i];
        }
      }
      for ( k=0; k<block_nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_nlev)*sl_param.np*sl_param.np;
        for ( i=0; i< sl_param.np; i++ ) {
          sl_elem_v1_north[kptr+i] = sl_elem_v1[vptr+(sl_param.np-1)*sl_param.np+sl_param.np-1-i];
        }
      }
    }

    pe_put(sl_param.buf+v1_ptr+sl_iwesn[ie*4+3], sl_elem_v1_north, blockbytes);
    pe_put(sl_param.buf+v2_ptr+sl_iwesn[ie*4+3], sl_elem_v2_north, block2bytes);
    
    if ( sl_param.rsplit > 0 ) {
    
      for ( k=0; k < block_nlev; k++ ) {
        kptr = k*sl_param.np;
	  vptr = (k+start_nlev)*sl_param.np*sl_param.np;
        for ( i=0; i< sl_param.np; i++ ) {
	  //WEST
          sl_elem_v3_west[kptr+i] = sl_elem_v3[vptr+i*sl_param.np];
	  //EAST
          sl_elem_v3_east[kptr+i] = sl_elem_v3[vptr+i*sl_param.np+sl_param.np-1]; 	
	  //SOUTH
          sl_elem_v3_south[kptr+i] = sl_elem_v3[vptr+i]; 	
	  //NORTH
          sl_elem_v3_north[kptr+i] = sl_elem_v3[vptr+(sl_param.np-1)*sl_param.np+i]; 	
        }
      }

	
	//WEST
	if (sl_reverse[ie*4] ) {
        for ( k=0; k < block_nlev; k++ ) {
          kptr = k*sl_param.np;
	    vptr = (k+start_nlev)*sl_param.np*sl_param.np;
          for ( i=0; i< sl_param.np; i++ ) {
            sl_elem_v3_west[kptr+i] = sl_elem_v3[vptr+(sl_param.np-i-1)*sl_param.np]; 	
          }
        } 
	}
	
      pe_put(sl_param.buf+v3_ptr+sl_iwesn[ie*4], sl_elem_v3_west, blockbytes);
	
    //EAST
	if ( sl_reverse[ie*4+1] ) {
        for ( k=0; k < block_nlev; k++ ) {
          kptr = k*sl_param.np;
	    vptr = (k+start_nlev)*sl_param.np*sl_param.np;
          for ( i=0; i< sl_param.np; i++ ) {
            sl_elem_v3_east[kptr+i] = sl_elem_v3[vptr+(sl_param.np-i-1)*sl_param.np+sl_param.np-1];
          }
        }
      }
	
     pe_put(sl_param.buf+v3_ptr+sl_iwesn[ie*4+1], sl_elem_v3_east, blockbytes);
	
	//SOUTH
	if ( sl_reverse[ie*4+2] ) {
        for ( k=0; k < block_nlev; k++ ) {
          kptr = k*sl_param.np;
	    vptr = (k+start_nlev)*sl_param.np*sl_param.np;
          for ( i=0; i< sl_param.np; i++ ) {
            sl_elem_v1_south[kptr+i] = sl_elem_v1[vptr+sl_param.np-i-1];
            sl_elem_v3_south[kptr+i] = sl_elem_v3[vptr+sl_param.np-i-1];	
          }
        }
      }
	
      pe_put(sl_param.buf+v3_ptr+sl_iwesn[ie*4+2], sl_elem_v3_south, blockbytes);
	
	//NORTH
      if ( sl_reverse[ie*4+3] ) {
        for ( k=0; k<block_nlev; k++ ) {
          kptr = k*sl_param.np;
	    vptr = (k+start_nlev)*sl_param.np*sl_param.np;
          for ( i=0; i< sl_param.np; i++ ) {
            sl_elem_v3_north[kptr+i] = sl_elem_v3[vptr+(sl_param.np-1)*sl_param.np+sl_param.np-1-i]; 	
          }
        }
	}
	
      pe_put(sl_param.buf+v3_ptr+sl_iwesn[ie*4+3], sl_elem_v3_north, blockbytes);
    }
    dma_syn();

    // SWEST
    if ( _COL == 0) {
    for ( j=0; j<sl_param.max_corner_elem; j++ ) {
      if ( sl_putmap[ie*4*sl_param.max_corner_elem+j] != -1) {
        edgeptr = sl_putmap[ie*4*sl_param.max_corner_elem+j];
        for (k=0; k<sl_param.nlev; k++){
          sl_elem_swest[k] = sl_elem_v1[k*sl_param.np*sl_param.np];
          sl_elem_swest[k+3*sl_param.nlev] = sl_elem_v3[k*sl_param.np*sl_param.np];
        }
        for ( k=0; k<2*sl_param.nlev; k++) {
          sl_elem_swest[k+sl_param.nlev] = sl_elem_v2[k*sl_param.np*sl_param.np];
        }
	  if ( sl_param.rsplit <= 0) {
          pe_put(sl_param.buf+1+edgeptr, sl_elem_swest, 3*sl_param.nlev*sizeof(FTYPE));
	  }else{
          pe_put(sl_param.buf+1+edgeptr, sl_elem_swest, 4*sl_param.nlev*sizeof(FTYPE));
	  }
        dma_syn();
      }
    }
    }
 
    // SEAST
    if( _COL == 1) {
    for ( j=sl_param.max_corner_elem; j<2*sl_param.max_corner_elem; j++ ) {
      if ( sl_putmap[ie*4*sl_param.max_corner_elem+j] != -1) {
        edgeptr = sl_putmap[ie*4*sl_param.max_corner_elem+j];
        for ( k =0; k<sl_param.nlev; k++ ) {
          sl_elem_seast[k] = sl_elem_v1[k*sl_param.np*sl_param.np+sl_param.np-1];
          sl_elem_seast[k+3*sl_param.nlev] = sl_elem_v3[k*sl_param.np*sl_param.np+sl_param.np-1];
        }
        for ( k=0; k<2*sl_param.nlev; k++){
          sl_elem_seast[k+sl_param.nlev] = sl_elem_v2[k*sl_param.np*sl_param.np+sl_param.np-1];
        }
	  if ( sl_param.rsplit <= 0) {
          pe_put(sl_param.buf+1+edgeptr, sl_elem_seast, 3*sl_param.nlev*sizeof(FTYPE));
	  }else{
          pe_put(sl_param.buf+1+edgeptr, sl_elem_seast, 4*sl_param.nlev*sizeof(FTYPE));
	  }
        dma_syn();
      }
    }
    }
    
    // NWEST
   if ( _COL == 2 ) {
    for ( j=2*sl_param.max_corner_elem; j<3*sl_param.max_corner_elem; j++ ) {
      if(sl_putmap[ie*4*sl_param.max_corner_elem+j] != -1 ) {
        edgeptr = sl_putmap[ie*4*sl_param.max_corner_elem+j];
	  for ( k=0; k<sl_param.nlev; k++ ){
          sl_elem_nwest[k] = sl_elem_v1[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np];
          sl_elem_nwest[k+3*sl_param.nlev] = sl_elem_v3[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np];
	  }
        for ( k=0; k<2*sl_param.nlev; k++ ) {
          sl_elem_nwest[k+sl_param.nlev] = sl_elem_v2[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np];
	  }
	  if (sl_param.rsplit <=0 ) {
          pe_put(sl_param.buf+1+edgeptr, sl_elem_nwest, 3*sl_param.nlev*sizeof(double));
	  }else{
          pe_put(sl_param.buf+1+edgeptr, sl_elem_nwest, 4*sl_param.nlev*sizeof(double));
	  }
        dma_syn();
      }
    }
    }
    // NEAST
    if ( _COL == 3) {
    for ( j =3*sl_param.max_corner_elem; j < 4*sl_param.max_corner_elem; j++) {
      if ( sl_putmap[ie*4*sl_param.max_corner_elem+j] != -1 ) {
        edgeptr = sl_putmap[ie*4*sl_param.max_corner_elem+j];
	  for ( k=0; k<sl_param.nlev; k++ ) {
          sl_elem_neast[k] = sl_elem_v1[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1];
          sl_elem_neast[k+3*sl_param.nlev] = sl_elem_v3[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1];
	  }
	  for (k=0; k<2*sl_param.nlev; k++ ) {
          sl_elem_neast[k+sl_param.nlev] = sl_elem_v2[k*sl_param.np*sl_param.np+(sl_param.np-1)*sl_param.np+sl_param.np-1];
	  }
	  if ( sl_param.rsplit <= 0) {
          pe_put(sl_param.buf+1+edgeptr, sl_elem_neast, 3*sl_param.nlev*sizeof(FTYPE));
	  }else{
          pe_put(sl_param.buf+1+edgeptr, sl_elem_neast, 4*sl_param.nlev*sizeof(FTYPE));
	  }
        dma_syn();
      }
    }
    }

    athread_syn(ROW_SCOPE, 0x00FF);
  }
  ldm_dealloc_after(sl_putmap);
}

