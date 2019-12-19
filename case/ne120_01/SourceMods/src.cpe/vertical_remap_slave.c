//#define LDM_FAST
//#define DMA_FAST
#include<stdio.h>
#include<slave.h>
#include<simd.h>
#include"dma_macros.h"
#include"ldm_alloc.h"
#include"vertical_remap_macro.h"

static inline void setftz(){
  unsigned long fpcr;
  asm volatile("rfpcr %0": "=r"(fpcr));
  fpcr |= 1L << 60 | 1L << 2;
  fpcr &= ~(1L << 48);
  asm volatile("wfpcr %0": : "r"(fpcr));
}

static inline void unsetftz(){
  unsigned long fpcr;
  asm volatile("rfpcr %0": "=r"(fpcr));
  fpcr |= 1L << 2 ;
  fpcr &= ~(1L << 60);
  asm volatile("wfpcr %0": : "r"(fpcr));
}


void vertical_remap_parallel_(vertical_remap_param_t *vertical_remap_ptr) { 
  setftz();
  vertical_remap_param_t vertical_remap_param ;
  dma_init() ;
  ldm_alloc_init() ;
 
  // bcast vertical_remap_param
  if(_MYID == 0) {
    bcast_get(vertical_remap_ptr , &vertical_remap_param , sizeof(vertical_remap_param_t)) ;
    dma_syn();	
  }
  athread_syn(ARRAY_SCOPE , 0xFFFF) ;
  
  // vertical_remap_param assignment
  double *hyai_ptr , *hybi_ptr ;
  unsigned long ie_address_diff ;	
  double ps0 ;
  int net_len , vert_remap_q_alg ; 
 
  hyai_ptr = vertical_remap_param.hyai_ptr ;
  hybi_ptr = vertical_remap_param.hybi_ptr ;
  ie_address_diff = vertical_remap_param.ie_address_diff ;
  ps0 = vertical_remap_param.ps0 ;
  net_len = vertical_remap_param.net_len ;
  vert_remap_q_alg = vertical_remap_param.vert_remap_q_alg ;

 
  // loop_len
  int nlev = NLEV ;
  // ldm malloc size
  int array_sizev ;
  int array_size0 , array_size1 , array_size2 , array_size3 ;
  int array_size4 , array_size5 , array_size6 , array_size7 ;
  array_size0 = 256 ;
  array_size1 = 512 ;
  array_size2 = 576 ;
  array_size3 = 1024 ;
  array_size4 = 2048 ;
  array_size5 = 2304 ;
  array_size6 = 5120 ;
  array_size7 = 3072 ;
  
  // dma trans len
  int dma_trans_len_v , dma_trans_len0 , dma_trans_len1 , dma_trans_len2 , dma_trans_len3;
  dma_trans_len_v = 16 ;
  dma_trans_len0 = 256 ;
  dma_trans_len1 = (_ROW != 7) ? 512 : 256 ;
  dma_trans_len2 = (_ROW != 7) ? 1024 : 512 ; 
  
  // ldm malloc
  // input and output array
  doublev4 *hyai_local_v4_ptr , *hybi_local_v4_ptr ;
  doublev4 *hyai_diff_local_v4_ptr , *hybi_diff_local_v4_ptr ;
  doublev4 *dp3d_local_v4_ptr , *dp3d_trans_local_v4_ptr , *dpo_local_v4_ptr;
  doublev4 *array_local_v4_ptr , *array_trans_local_v4_ptr ;
  // tmp array  
  doublev4 *dpn_local_v4_ptr , *dpn_trans_local_v4_ptr ;
  doublev4 *pin_local_v4_ptr , *pio_local_v4_ptr;
  doublev4 *ppmdx_local_v4_ptr ;
  doublev4 *masso_local_v4_ptr , *dma_local_v4_ptr , *ai_local_v4_ptr ;
  doublev4 *coefs_local_v4_ptr ;
  double *z_local_ptr ;
  int *kid_local_ptr ;
  ldm_alloc_pointer(hyai_local_v4_ptr , array_size0) ;
  ldm_alloc_pointer(hybi_local_v4_ptr , array_size0) ;
  ldm_alloc_pointer(hyai_diff_local_v4_ptr , array_size0) ;
  ldm_alloc_pointer(hybi_diff_local_v4_ptr , array_size0) ;
  ldm_alloc_pointer(kid_local_ptr , array_size0) ;
  
  ldm_alloc_pointer(dpo_local_v4_ptr , array_size2) ;
  
  ldm_alloc_pointer(dp3d_local_v4_ptr , array_size1) ; \
  ldm_alloc_pointer(dp3d_trans_local_v4_ptr , array_size1) ;
  ldm_alloc_pointer(dpn_local_v4_ptr , array_size1) ;
  ldm_alloc_pointer(dpn_trans_local_v4_ptr , array_size1) ;
  ldm_alloc_pointer(pin_local_v4_ptr , array_size1) ;
  ldm_alloc_pointer(pio_local_v4_ptr , array_size1) ;
  ldm_alloc_pointer(z_local_ptr , array_size1) ;
  
  ldm_alloc_pointer(dma_local_v4_ptr , array_size3) ;
  ldm_alloc_pointer(ai_local_v4_ptr , array_size3) ; 
  ldm_alloc_pointer(masso_local_v4_ptr , array_size3) ;
 
  ldm_alloc_pointer(array_trans_local_v4_ptr , array_size4) ; 
  ldm_alloc_pointer(array_local_v4_ptr , array_size5) ;
  ldm_alloc_pointer(ppmdx_local_v4_ptr , array_size6) ;
  ldm_alloc_pointer(coefs_local_v4_ptr , array_size7) ;
  // bcast hyai and hybi's array
  if (_MYID == 0 ) {
    bcast_get(hyai_ptr , hyai_local_v4_ptr , dma_trans_len0) ;
    bcast_get(hybi_ptr , hybi_local_v4_ptr , dma_trans_len0) ;
    dma_syn() ;
  }
  athread_syn(ARRAY_SCOPE , 0xffff);
  // compute offset of send buffer and recv buffer 
  intv8 send_offset_v8_dp3d , recv_offset_v8_dp3d , send_id_v8 ;
  intv8 send_offset_v8_dpn , recv_offset_v8_dpn ;
  intv8 send_offset_v8_forward_qdp , recv_offset_v8_forward_qdp ;
  intv8 send_offset_v8_post_qdp , recv_offset_v8_post_qdp ;
  COMPUTE_SEND_RECV_OFFSET(send_id_v8 , \
		           send_offset_v8_dp3d , recv_offset_v8_dp3d , \
  	                   send_offset_v8_dpn , recv_offset_v8_dpn , \
                           send_offset_v8_forward_qdp , recv_offset_v8_forward_qdp , \
                           send_offset_v8_post_qdp , recv_offset_v8_post_qdp , \
                            _ROW) ;
 
  // compute hyai(1) * ps0
  double *hyai_local_ptr ;
  double hyai_var , hyai_ps0_var ;
  doublev4 ps0_var_v4 ;
  doublev4 hyai_ps0_var_v4 ;
  hyai_local_ptr = (double *)hyai_local_v4_ptr ;
  hyai_var = hyai_local_ptr[0] ;
  hyai_ps0_var = hyai_var*ps0 ;
  SET_ZERO_DOUBLEV4_HI(hyai_ps0_var_v4 , hyai_ps0_var , 0x0) ;
  VSHFF(ps0_var_v4 , ps0 , ps0 , 0x0) ;
       
 
  // compute difference between hyai(k+1) and hyai(k) 
  COMPUTE_HYAI_DIFFERENCE(hyai_diff_local_v4_ptr , hyai_local_v4_ptr , ps0_var_v4) ; 
  // compute difference between hybi(k+1) and hybi(k)
  COMPUTE_HYBI_DIFFERENCE(hybi_diff_local_v4_ptr , hybi_local_v4_ptr) ;
 
  int ie ;
  // when q != 0 , Qdp_q_ptr is the ptr of Qdp array 
  remap_q_ppm_param_t remap_q_ppm_param ;
  remap_q_ppm_param.hyai_ps0_var_v4 = hyai_ps0_var_v4 ;
  remap_q_ppm_param.send_id_v8 = send_id_v8 ;
  remap_q_ppm_param.send_offset_v8_dp3d = send_offset_v8_dp3d ;
  remap_q_ppm_param.recv_offset_v8_dp3d = recv_offset_v8_dp3d ;
  remap_q_ppm_param.send_offset_v8_dpn = send_offset_v8_dpn ;
  remap_q_ppm_param.recv_offset_v8_dpn = recv_offset_v8_dpn ;
  remap_q_ppm_param.send_offset_v8_forward_qdp = send_offset_v8_forward_qdp ;
  remap_q_ppm_param.recv_offset_v8_forward_qdp = recv_offset_v8_forward_qdp ;
  remap_q_ppm_param.send_offset_v8_post_qdp = send_offset_v8_post_qdp ;
  remap_q_ppm_param.recv_offset_v8_post_qdp = recv_offset_v8_post_qdp ;
  remap_q_ppm_param.hyai_diff_local_v4_ptr = hyai_diff_local_v4_ptr ; 
  remap_q_ppm_param.hybi_diff_local_v4_ptr = hybi_diff_local_v4_ptr ;
  remap_q_ppm_param.dp3d_local_v4_ptr = dp3d_local_v4_ptr ; 
  remap_q_ppm_param.dp3d_trans_local_v4_ptr = dp3d_trans_local_v4_ptr ;
  remap_q_ppm_param.dpo_local_v4_ptr = dpo_local_v4_ptr ;
  remap_q_ppm_param.dpn_local_v4_ptr = dpn_local_v4_ptr ;
  remap_q_ppm_param.dpn_trans_local_v4_ptr = dpn_trans_local_v4_ptr ; 
  remap_q_ppm_param.pio_local_v4_ptr = pio_local_v4_ptr ;
  remap_q_ppm_param.pin_local_v4_ptr = pin_local_v4_ptr ;
  remap_q_ppm_param.ppmdx_local_v4_ptr = ppmdx_local_v4_ptr ;
  remap_q_ppm_param.array_local_v4_ptr = array_local_v4_ptr ;
  remap_q_ppm_param.array_trans_local_v4_ptr = array_trans_local_v4_ptr ;
  remap_q_ppm_param.masso_local_v4_ptr = masso_local_v4_ptr ;
  remap_q_ppm_param.coefs_local_v4_ptr = coefs_local_v4_ptr ;
  remap_q_ppm_param.ai_local_v4_ptr = ai_local_v4_ptr ;
  remap_q_ppm_param.dma_local_v4_ptr = dma_local_v4_ptr ;  
  remap_q_ppm_param.z_local_ptr = z_local_ptr ;
  remap_q_ppm_param.kid_local_ptr = kid_local_ptr ; 
  for (ie = (_COL >> 2) ; ie < net_len ; ie = ie+2 ) {
        double *dp3d_ptr = (double *)((unsigned long)vertical_remap_param.dp3d_ptr + ie*ie_address_diff) ;
        double *ps_v_ptr = (double *)((unsigned long)vertical_remap_param.ps_v_ptr + ie*ie_address_diff) ;
        double *v_ptr = (double *)((unsigned long)vertical_remap_param.v_ptr + ie*ie_address_diff) ;
        double *t_ptr = (double *)((unsigned long)vertical_remap_param.t_ptr + ie*ie_address_diff) ;
        double *Qdp_ptr = (double *)((unsigned long)vertical_remap_param.Qdp_ptr + ie*ie_address_diff) ;
        double *Qdp_q_ptr = (double *)(Qdp_ptr + 480) ; 
  
 
        // dma trans dp3d 
	pe_get((dp3d_ptr + (_ROW << 6)) , dp3d_local_v4_ptr , dma_trans_len1) ;
	dma_syn() ;
        TRANSPOSE_DP3D(&remap_q_ppm_param) ;
        VERTICAL_REMAP_INIT(&remap_q_ppm_param , ie ) ;
        TRANSPOSE_DPN(&remap_q_ppm_param) ;
        pe_put(ps_v_ptr + (_ROW << 1) , &(remap_q_ppm_param.ps_v_var_v4) , dma_trans_len_v) ;
        dma_syn() ; 
        int q ;
	for (q = -4+((_COL - ((_COL >> 2) <<2))<<2) ; q < 24 ; q = q + 16) {
          double *array_ptr_0 , *array_ptr_1 , *array_ptr_2 , *array_ptr_3 ;
	  doublev4 *array_local_v4_ptr_0 , *array_local_v4_ptr_1 , *array_local_v4_ptr_2 , *array_local_v4_ptr_3 ;
	  if (q == -4) {
            array_local_v4_ptr_0 = array_local_v4_ptr ;
            array_local_v4_ptr_1 = &array_local_v4_ptr[8] ;
            array_local_v4_ptr_2 = &array_local_v4_ptr[32] ;
            array_local_v4_ptr_3 = &array_local_v4_ptr[48] ;
	    array_ptr_0 = v_ptr + (_ROW << 7) ;
	    array_ptr_2 = t_ptr + (_ROW << 6) ;
	    array_ptr_3 = Qdp_ptr + (_ROW << 6) ;

	    // v1 , v2 , t , Qdp
            pe_get(array_ptr_0 , array_local_v4_ptr_0 , dma_trans_len2) ;
            pe_get(array_ptr_2 , array_local_v4_ptr_2 , dma_trans_len1) ;
            pe_get(array_ptr_3 , array_local_v4_ptr_3 , dma_trans_len1) ;
            dma_syn() ;
          

            TRANSPOSE_ARRAY_FORWARD(&remap_q_ppm_param) ;
            REMAP_Q_PPM(&remap_q_ppm_param , ie , q) ;
            TRANSPOSE_ARRAY_POST(&remap_q_ppm_param) ;
            pe_put(array_ptr_0 , array_local_v4_ptr_0 , dma_trans_len2) ;
            pe_put(array_ptr_2 , array_local_v4_ptr_2 , dma_trans_len1) ;
            pe_put(array_ptr_3 , array_local_v4_ptr_3 , dma_trans_len1) ;
            dma_syn() ; 
          } else {
	    array_ptr_0 = Qdp_q_ptr + q*480 + (_ROW << 6) ;
	    array_ptr_1 = Qdp_q_ptr + (q+1)*480 + (_ROW << 6) ;
	    array_ptr_2 = Qdp_q_ptr + (q+2)*480 + (_ROW << 6) ;
	    array_ptr_3 = Qdp_q_ptr + (q+3)*480 + (_ROW << 6) ;
            array_local_v4_ptr_0 = array_local_v4_ptr ;
            array_local_v4_ptr_1 = &array_local_v4_ptr[16] ;
            array_local_v4_ptr_2 = &array_local_v4_ptr[32] ;
            array_local_v4_ptr_3 = &array_local_v4_ptr[48] ;
		// v1 , v2 , t , Qdp
            pe_get(array_ptr_0 , array_local_v4_ptr_0 , dma_trans_len1) ;
            pe_get(array_ptr_1 , array_local_v4_ptr_1 , dma_trans_len1) ;
            pe_get(array_ptr_2 , array_local_v4_ptr_2 , dma_trans_len1) ;
            pe_get(array_ptr_3 , array_local_v4_ptr_3 , dma_trans_len1) ;
            dma_syn() ;
            TRANSPOSE_ARRAY_FORWARD_QDP(&remap_q_ppm_param) ;
            REMAP_Q_PPM(&remap_q_ppm_param , ie , q) ;
            TRANSPOSE_ARRAY_POST_QDP(&remap_q_ppm_param) ;
            pe_put(array_ptr_0 , array_local_v4_ptr_0 , dma_trans_len1) ;
            pe_put(array_ptr_1 , array_local_v4_ptr_1 , dma_trans_len1) ;
            pe_put(array_ptr_2 , array_local_v4_ptr_2 , dma_trans_len1) ;
            pe_put(array_ptr_3 , array_local_v4_ptr_3 , dma_trans_len1) ;
            dma_syn() ;
 
	 }
       }   
    }
    unsetftz();
}  
