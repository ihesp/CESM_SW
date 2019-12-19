#include<math.h>
// 结构体定义;
typedef struct {
  double *hyai_ptr , *hybi_ptr ;
  double *dp3d_ptr , *ps_v_ptr ;
  double *v_ptr , *t_ptr , *Qdp_ptr ;
  unsigned long ie_address_diff ; 
  double ps0 ;
  int net_len , vert_remap_q_alg ;
} vertical_remap_param_t ;

typedef struct {
  doublev4 hyai_ps0_var_v4 ;
  doublev4 ps_v_var_v4 ;
  intv8 send_id_v8 ;
  intv8 send_offset_v8_dp3d ; 
  intv8 recv_offset_v8_dp3d ;
  intv8 send_offset_v8_dpn ;
  intv8 recv_offset_v8_dpn ;
  intv8 send_offset_v8_forward_qdp ;
  intv8 recv_offset_v8_forward_qdp ;
  intv8 send_offset_v8_post_qdp ;
  intv8 recv_offset_v8_post_qdp ;
  doublev4 *hyai_diff_local_v4_ptr ; 
  doublev4 *hybi_diff_local_v4_ptr ; 
  doublev4 *dp3d_local_v4_ptr ;
  doublev4 *dp3d_trans_local_v4_ptr ;
  doublev4 *dpo_local_v4_ptr ;
  doublev4 *dpn_local_v4_ptr ; 
  doublev4 *dpn_trans_local_v4_ptr ; 
  doublev4 *pio_local_v4_ptr ; 
  doublev4 *pin_local_v4_ptr ;
  doublev4 *ppmdx_local_v4_ptr ;
  doublev4 *array_local_v4_ptr ;
  doublev4 *array_trans_local_v4_ptr ;
  doublev4 *masso_local_v4_ptr ;
  doublev4 *coefs_local_v4_ptr ;
  doublev4 *ai_local_v4_ptr ; 
  doublev4 *dma_local_v4_ptr ;  
  double *z_local_ptr ;
  int *kid_local_ptr ;
} remap_q_ppm_param_t ;
// 宏定义文件
#define printf cpe_printf
#define NLEV_1 29
#define NLEV 30
#define NLEV_ALIGN 32
#define VSHFF(dst_v4 , src_hi_v4 , src_lo_v4 , mask) asm("vshff %1 , %2 , %3 , %0" : "=r"(dst_v4) : "r"(src_hi_v4) , "r"(src_lo_v4) , "r"(mask))
#define VSHFW(dst_v8 , src_hi_v8 , src_lo_v8 , mask) asm("vshfw %1 , %2 , %3 , %0" : "=r"(dst_v8) : "r"(src_hi_v8) , "r"(src_lo_v8) , "r"(mask))
#define SET_ZERO_DOUBLEV4_HI(dst_v4 , src_lo_v4 , mask) asm("vshff $31 , %1 , %2 , %0" : "=r"(dst_v4) : "r"(src_lo_v4) , "r"(mask))
#define SET_ZERO_DOUBLEV4_LO(dst_v4 , src_hi_v4 , mask) asm("vshff %1 , $31 , %2 , %0" : "=r"(dst_v4) : "r"(src_hi_v4) , "r"(mask))
// 寄存器通信_列发送
#define REGISTER_PUTC(send_v4 , send_id) asm("putc %0 , %1" : : "r"(send_v4) , "r"(send_id))

// 寄存器通信_列接受
#define REGISTER_GETC(recv_v4) asm("getc %0" : "=r"(recv_v4))

// 列同步
#define REGISTER_SYNC(mask) asm("sync %0" : : "r"(mask))
// abs
#define VABS(abs_var , var) asm volatile("vcpys $31 , %1 , %0" : "=r"(abs_var) : "r"(var))
// sign 
#define VSIGN(sign_var , sign , var) asm("vcpys %1 , %2 , %0" : "=r"(sign_var) : "r"(sign) , "r"(var))

#define VMIN(sellt_rslt , test_var , sellt_var1 , sellt_var2) asm("vsellt %1 , %2 , %3 , %0" : "=r"(sellt_rslt) : "r"(test_var) , "r"(sellt_var1) , "r"(sellt_var2)) 

#define VMINEZ(sellt_rslt , test_var , sellt_var) asm("vselle %1 , $31 , %2 , %0" : "=r"(sellt_rslt) : "r"(test_var) , "r"(sellt_var))

#define VMINE(sellt_rslt , test_var , sellt_var1 , sellt_var2) asm("vselle %1 , %2 , %3 , %0" : "=r"(sellt_rslt) : "r"(test_var) , "r"(sellt_var1) , "r"(sellt_var2))
// hyai[k+1] - hyai[k]
// hybi[k+1] - hybi[k]
#define COMPUTE_HYAI_DIFFERENCE(hyai_diff_local_v4_ptr , hyai_local_v4_ptr , ps0_var_v4) { \
    doublev4 hyai_var_v4_00 , hyai_var_v4_01 ;                        \
    doublev4 hyai_var_v4_02 , hyai_var_v4_03 ;                        \
    doublev4 hyai_var_v4_04 , hyai_var_v4_05 ;                        \
    doublev4 hyai_var_v4_06 , hyai_var_v4_07 ;                        \
                                                                        \
    hyai_var_v4_00 = hyai_local_v4_ptr[0] ;                           \
    hyai_var_v4_01 = hyai_local_v4_ptr[1] ;                           \
    hyai_var_v4_02 = hyai_local_v4_ptr[2] ;                           \
    hyai_var_v4_03 = hyai_local_v4_ptr[3] ;                           \
    hyai_var_v4_04 = hyai_local_v4_ptr[4] ;                           \
    hyai_var_v4_05 = hyai_local_v4_ptr[5] ;                           \
    hyai_var_v4_06 = hyai_local_v4_ptr[6] ;                           \
    hyai_var_v4_07 = hyai_local_v4_ptr[7] ;                           \
                                                                        \
    doublev4 tmpv4_0 , tmpv4_1 , tmpv4_2 , tmpv4_3 ;                    \
    doublev4 tmpv4_4 , tmpv4_5 , tmpv4_6 , tmpv4_7 ;                    \
                                                                        \
    tmpv4_0 = simd_vshff(hyai_var_v4_01 , hyai_var_v4_00 , 0x4e) ;    \
    tmpv4_1 = simd_vshff(hyai_var_v4_02 , hyai_var_v4_01 , 0x4e) ;    \
    tmpv4_2 = simd_vshff(hyai_var_v4_03 , hyai_var_v4_02 , 0x4e) ;    \
    tmpv4_3 = simd_vshff(hyai_var_v4_04 , hyai_var_v4_03 , 0x4e) ;    \
    tmpv4_4 = simd_vshff(hyai_var_v4_05 , hyai_var_v4_04 , 0x4e) ;    \
    tmpv4_5 = simd_vshff(hyai_var_v4_06 , hyai_var_v4_05 , 0x4e) ;    \
    tmpv4_6 = simd_vshff(hyai_var_v4_07 , hyai_var_v4_06 , 0x4e) ;    \
    tmpv4_7 = simd_vshff(hyai_var_v4_07 , hyai_var_v4_07 , 0x4e) ;    \
                                                                        \
    doublev4 hyai_var_v4_10 , hyai_var_v4_11 ;                        \
    doublev4 hyai_var_v4_12 , hyai_var_v4_13 ;                        \
    doublev4 hyai_var_v4_14 , hyai_var_v4_15 ;                        \
    doublev4 hyai_var_v4_16 , hyai_var_v4_17 ;                        \
                                                                        \
    hyai_var_v4_10 = simd_vshff(tmpv4_0 , hyai_var_v4_00 , 0x99) ;    \
    hyai_var_v4_11 = simd_vshff(tmpv4_1 , hyai_var_v4_01 , 0x99) ;    \
    hyai_var_v4_12 = simd_vshff(tmpv4_2 , hyai_var_v4_02 , 0x99) ;    \
    hyai_var_v4_13 = simd_vshff(tmpv4_3 , hyai_var_v4_03 , 0x99) ;    \
    hyai_var_v4_14 = simd_vshff(tmpv4_4 , hyai_var_v4_04 , 0x99) ;    \
    hyai_var_v4_15 = simd_vshff(tmpv4_5 , hyai_var_v4_05 , 0x99) ;    \
    hyai_var_v4_16 = simd_vshff(tmpv4_6 , hyai_var_v4_06 , 0x99) ;    \
    hyai_var_v4_17 = simd_vshff(tmpv4_7 , hyai_var_v4_07 , 0x99) ;    \
                                                                        \
    doublev4 hyai_diff_v4_0 , hyai_diff_v4_1 ;                        \
    doublev4 hyai_diff_v4_2 , hyai_diff_v4_3 ;                        \
    doublev4 hyai_diff_v4_4 , hyai_diff_v4_5 ;                        \
    doublev4 hyai_diff_v4_6 , hyai_diff_v4_7 ;                        \
                                                                        \
    hyai_diff_v4_0 = hyai_var_v4_10 - hyai_var_v4_00 ;               \
    hyai_diff_v4_1 = hyai_var_v4_11 - hyai_var_v4_01 ;               \
    hyai_diff_v4_2 = hyai_var_v4_12 - hyai_var_v4_02 ;               \
    hyai_diff_v4_3 = hyai_var_v4_13 - hyai_var_v4_03 ;               \
    hyai_diff_v4_4 = hyai_var_v4_14 - hyai_var_v4_04 ;               \
    hyai_diff_v4_5 = hyai_var_v4_15 - hyai_var_v4_05 ;               \
    hyai_diff_v4_6 = hyai_var_v4_16 - hyai_var_v4_06 ;               \
    hyai_diff_v4_7 = hyai_var_v4_17 - hyai_var_v4_07 ;               \
                                                                        \
    SET_ZERO_DOUBLEV4_HI(hyai_diff_v4_7 , hyai_diff_v4_7 , 0x44) ;    \
                                                                        \
    hyai_diff_v4_0 = hyai_diff_v4_0*ps0_var_v4 ;                      \
    hyai_diff_v4_1 = hyai_diff_v4_1*ps0_var_v4 ;                      \
    hyai_diff_v4_2 = hyai_diff_v4_2*ps0_var_v4 ;                      \
    hyai_diff_v4_3 = hyai_diff_v4_3*ps0_var_v4 ;                      \
    hyai_diff_v4_4 = hyai_diff_v4_4*ps0_var_v4 ;                      \
    hyai_diff_v4_5 = hyai_diff_v4_5*ps0_var_v4 ;                      \
    hyai_diff_v4_6 = hyai_diff_v4_6*ps0_var_v4 ;                      \
    hyai_diff_v4_7 = hyai_diff_v4_7*ps0_var_v4 ;                      \
    									\
    hyai_diff_local_v4_ptr[0] = hyai_diff_v4_0 ;                      \
    hyai_diff_local_v4_ptr[1] = hyai_diff_v4_1 ;                      \
    hyai_diff_local_v4_ptr[2] = hyai_diff_v4_2 ;                      \
    hyai_diff_local_v4_ptr[3] = hyai_diff_v4_3 ;                      \
    hyai_diff_local_v4_ptr[4] = hyai_diff_v4_4 ;                      \
    hyai_diff_local_v4_ptr[5] = hyai_diff_v4_5 ;                      \
    hyai_diff_local_v4_ptr[6] = hyai_diff_v4_6 ;                      \
    hyai_diff_local_v4_ptr[7] = hyai_diff_v4_7 ;                      \
  }

  #define COMPUTE_HYBI_DIFFERENCE(hybi_diff_local_v4_ptr , hybi_local_v4_ptr) { \
    doublev4 hybi_var_v4_00 , hybi_var_v4_01 ;                        \
    doublev4 hybi_var_v4_02 , hybi_var_v4_03 ;                        \
    doublev4 hybi_var_v4_04 , hybi_var_v4_05 ;                        \
    doublev4 hybi_var_v4_06 , hybi_var_v4_07 ;                        \
                                                                        \
    hybi_var_v4_00 = hybi_local_v4_ptr[0] ;                           \
    hybi_var_v4_01 = hybi_local_v4_ptr[1] ;                           \
    hybi_var_v4_02 = hybi_local_v4_ptr[2] ;                           \
    hybi_var_v4_03 = hybi_local_v4_ptr[3] ;                           \
    hybi_var_v4_04 = hybi_local_v4_ptr[4] ;                           \
    hybi_var_v4_05 = hybi_local_v4_ptr[5] ;                           \
    hybi_var_v4_06 = hybi_local_v4_ptr[6] ;                           \
    hybi_var_v4_07 = hybi_local_v4_ptr[7] ;                           \
                                                                        \
    doublev4 tmpv4_0 , tmpv4_1 , tmpv4_2 , tmpv4_3 ;                    \
    doublev4 tmpv4_4 , tmpv4_5 , tmpv4_6 , tmpv4_7 ;                    \
                                                                        \
    tmpv4_0 = simd_vshff(hybi_var_v4_01 , hybi_var_v4_00 , 0x4e) ;    \
    tmpv4_1 = simd_vshff(hybi_var_v4_02 , hybi_var_v4_01 , 0x4e) ;    \
    tmpv4_2 = simd_vshff(hybi_var_v4_03 , hybi_var_v4_02 , 0x4e) ;    \
    tmpv4_3 = simd_vshff(hybi_var_v4_04 , hybi_var_v4_03 , 0x4e) ;    \
    tmpv4_4 = simd_vshff(hybi_var_v4_05 , hybi_var_v4_04 , 0x4e) ;    \
    tmpv4_5 = simd_vshff(hybi_var_v4_06 , hybi_var_v4_05 , 0x4e) ;    \
    tmpv4_6 = simd_vshff(hybi_var_v4_07 , hybi_var_v4_06 , 0x4e) ;    \
    tmpv4_7 = simd_vshff(hybi_var_v4_07 , hybi_var_v4_07 , 0x4e) ;    \
                                                                        \
    doublev4 hybi_var_v4_10 , hybi_var_v4_11 ;                        \
    doublev4 hybi_var_v4_12 , hybi_var_v4_13 ;                        \
    doublev4 hybi_var_v4_14 , hybi_var_v4_15 ;                        \
    doublev4 hybi_var_v4_16 , hybi_var_v4_17 ;                        \
                                                                        \
    hybi_var_v4_10 = simd_vshff(tmpv4_0 , hybi_var_v4_00 , 0x99) ;    \
    hybi_var_v4_11 = simd_vshff(tmpv4_1 , hybi_var_v4_01 , 0x99) ;    \
    hybi_var_v4_12 = simd_vshff(tmpv4_2 , hybi_var_v4_02 , 0x99) ;    \
    hybi_var_v4_13 = simd_vshff(tmpv4_3 , hybi_var_v4_03 , 0x99) ;    \
    hybi_var_v4_14 = simd_vshff(tmpv4_4 , hybi_var_v4_04 , 0x99) ;    \
    hybi_var_v4_15 = simd_vshff(tmpv4_5 , hybi_var_v4_05 , 0x99) ;    \
    hybi_var_v4_16 = simd_vshff(tmpv4_6 , hybi_var_v4_06 , 0x99) ;    \
    hybi_var_v4_17 = simd_vshff(tmpv4_7 , hybi_var_v4_07 , 0x99) ;    \
                                                                        \
    doublev4 hybi_diff_v4_0 , hybi_diff_v4_1 ;                        \
    doublev4 hybi_diff_v4_2 , hybi_diff_v4_3 ;                        \
    doublev4 hybi_diff_v4_4 , hybi_diff_v4_5 ;                        \
    doublev4 hybi_diff_v4_6 , hybi_diff_v4_7 ;                        \
    							                \
    hybi_diff_v4_0 = hybi_var_v4_10 - hybi_var_v4_00 ;               \
    hybi_diff_v4_1 = hybi_var_v4_11 - hybi_var_v4_01 ;               \
    hybi_diff_v4_2 = hybi_var_v4_12 - hybi_var_v4_02 ;               \
    hybi_diff_v4_3 = hybi_var_v4_13 - hybi_var_v4_03 ;               \
    hybi_diff_v4_4 = hybi_var_v4_14 - hybi_var_v4_04 ;               \
    hybi_diff_v4_5 = hybi_var_v4_15 - hybi_var_v4_05 ;               \
    hybi_diff_v4_6 = hybi_var_v4_16 - hybi_var_v4_06 ;               \
    hybi_diff_v4_7 = hybi_var_v4_17 - hybi_var_v4_07 ;               \
                                                                        \
    SET_ZERO_DOUBLEV4_HI(hybi_diff_v4_7 , hybi_diff_v4_7 , 0x44) ;    \
                                                                        \
    hybi_diff_local_v4_ptr[0] = hybi_diff_v4_0 ;                      \
    hybi_diff_local_v4_ptr[1] = hybi_diff_v4_1 ;                      \
    hybi_diff_local_v4_ptr[2] = hybi_diff_v4_2 ;                      \
    hybi_diff_local_v4_ptr[3] = hybi_diff_v4_3 ;                      \
    hybi_diff_local_v4_ptr[4] = hybi_diff_v4_4 ;                      \
    hybi_diff_local_v4_ptr[5] = hybi_diff_v4_5 ;                      \
    hybi_diff_local_v4_ptr[6] = hybi_diff_v4_6 ;                      \
    hybi_diff_local_v4_ptr[7] = hybi_diff_v4_7 ;                      \
  } 
   

  #define TRANSPOSE_ARRAY_DP3D(dp3d_trans_local_v4_ptr , dp3d_local_v4_ptr) { \
  doublev4 dp3d_var_v4_0 , dp3d_var_v4_1 , dp3d_var_v4_2 , dp3d_var_v4_3 ;     \
  doublev4 dp3d_var_v4_4 , dp3d_var_v4_5 , dp3d_var_v4_6 , dp3d_var_v4_7 ;     \
  doublev4 dp3d_var_v4_8 , dp3d_var_v4_9 , dp3d_var_v4_10 , dp3d_var_v4_11 ;   \
  doublev4 dp3d_var_v4_12 , dp3d_var_v4_13 , dp3d_var_v4_14 , dp3d_var_v4_15 ; \
  \
  dp3d_var_v4_0 = dp3d_local_v4_ptr[0] ;  \
  dp3d_var_v4_1 = dp3d_local_v4_ptr[1] ;  \
  dp3d_var_v4_4 = dp3d_local_v4_ptr[4] ;  \
  dp3d_var_v4_5 = dp3d_local_v4_ptr[5] ;  \
  dp3d_var_v4_8 = dp3d_local_v4_ptr[8] ;  \
  dp3d_var_v4_9 = dp3d_local_v4_ptr[9] ;  \
  dp3d_var_v4_12 = dp3d_local_v4_ptr[12] ; \
  dp3d_var_v4_13 = dp3d_local_v4_ptr[13] ; \
  \
  doublev4 dtmpv4_0 , dtmpv4_1 , dtmpv4_2 , dtmpv4_3 ; \
  doublev4 dtmpv4_4 , dtmpv4_5 , dtmpv4_6 , dtmpv4_7 ; \
  \
  dtmpv4_0 = simd_vshff(dp3d_var_v4_4 , dp3d_var_v4_0 , 0x44) ; \
  dtmpv4_1 = simd_vshff(dp3d_var_v4_4 , dp3d_var_v4_0 , 0xee) ; \
  dtmpv4_2 = simd_vshff(dp3d_var_v4_5 , dp3d_var_v4_1 , 0x44) ; \
  dtmpv4_3 = simd_vshff(dp3d_var_v4_5 , dp3d_var_v4_1 , 0xee) ; \
  dtmpv4_4 = simd_vshff(dp3d_var_v4_12 , dp3d_var_v4_8 , 0x44) ; \
  dtmpv4_5 = simd_vshff(dp3d_var_v4_12 , dp3d_var_v4_8 , 0xee) ; \
  dtmpv4_6 = simd_vshff(dp3d_var_v4_13 , dp3d_var_v4_9 , 0x44) ; \
  dtmpv4_7 = simd_vshff(dp3d_var_v4_13 , dp3d_var_v4_9 , 0xee) ; \
  \
  doublev4 dp3d_trans_var_v4_0 , dp3d_trans_var_v4_1 , dp3d_trans_var_v4_2 , dp3d_trans_var_v4_3 ; \
  doublev4 dp3d_trans_var_v4_4 , dp3d_trans_var_v4_5 , dp3d_trans_var_v4_6 , dp3d_trans_var_v4_7 ; \
  \
  dp3d_trans_var_v4_0 = simd_vshff(dtmpv4_4 , dtmpv4_0 , 0x88) ; \
  dp3d_trans_var_v4_1 = simd_vshff(dtmpv4_4 , dtmpv4_0 , 0xdd) ; \
  dp3d_trans_var_v4_2 = simd_vshff(dtmpv4_5 , dtmpv4_1 , 0x88) ; \
  dp3d_trans_var_v4_3 = simd_vshff(dtmpv4_5 , dtmpv4_1 , 0xdd) ; \
  dp3d_trans_var_v4_4 = simd_vshff(dtmpv4_6 , dtmpv4_2 , 0x88) ; \
  dp3d_trans_var_v4_5 = simd_vshff(dtmpv4_6 , dtmpv4_2 , 0xdd) ; \
  dp3d_trans_var_v4_6 = simd_vshff(dtmpv4_7 , dtmpv4_3 , 0x88) ; \
  dp3d_trans_var_v4_7 = simd_vshff(dtmpv4_7 , dtmpv4_3 , 0xdd) ; \
  \
  dp3d_trans_local_v4_ptr[0] = dp3d_trans_var_v4_0 ; \
  dp3d_trans_local_v4_ptr[1] = dp3d_trans_var_v4_1 ; \
  dp3d_trans_local_v4_ptr[2] = dp3d_trans_var_v4_2 ; \
  dp3d_trans_local_v4_ptr[3] = dp3d_trans_var_v4_3 ; \
  dp3d_trans_local_v4_ptr[4] = dp3d_trans_var_v4_4 ; \
  dp3d_trans_local_v4_ptr[5] = dp3d_trans_var_v4_5 ; \
  dp3d_trans_local_v4_ptr[6] = dp3d_trans_var_v4_6 ; \
  dp3d_trans_local_v4_ptr[7] = dp3d_trans_var_v4_7 ; \
  \
  dp3d_var_v4_2 = dp3d_local_v4_ptr[2] ;  \
  dp3d_var_v4_3 = dp3d_local_v4_ptr[3] ;  \
  dp3d_var_v4_6 = dp3d_local_v4_ptr[6] ;  \
  dp3d_var_v4_7 = dp3d_local_v4_ptr[7] ;  \
  dp3d_var_v4_10 = dp3d_local_v4_ptr[10] ;  \
  dp3d_var_v4_11 = dp3d_local_v4_ptr[11] ;  \
  dp3d_var_v4_14 = dp3d_local_v4_ptr[14] ; \
  dp3d_var_v4_15 = dp3d_local_v4_ptr[15] ; \
  \
  dtmpv4_0 = simd_vshff(dp3d_var_v4_6 , dp3d_var_v4_2 , 0x44) ; \
  dtmpv4_1 = simd_vshff(dp3d_var_v4_6 , dp3d_var_v4_2 , 0xee) ; \
  dtmpv4_2 = simd_vshff(dp3d_var_v4_7 , dp3d_var_v4_3 , 0x44) ; \
  dtmpv4_3 = simd_vshff(dp3d_var_v4_7 , dp3d_var_v4_3 , 0xee) ; \
  dtmpv4_4 = simd_vshff(dp3d_var_v4_14 , dp3d_var_v4_10 , 0x44) ; \
  dtmpv4_5 = simd_vshff(dp3d_var_v4_14 , dp3d_var_v4_10 , 0xee) ; \
  dtmpv4_6 = simd_vshff(dp3d_var_v4_15 , dp3d_var_v4_11 , 0x44) ; \
  dtmpv4_7 = simd_vshff(dp3d_var_v4_15 , dp3d_var_v4_11 , 0xee) ; \
  \
  dp3d_trans_var_v4_0 = simd_vshff(dtmpv4_4 , dtmpv4_0 , 0x88) ; \
  dp3d_trans_var_v4_1 = simd_vshff(dtmpv4_4 , dtmpv4_0 , 0xdd) ; \
  dp3d_trans_var_v4_2 = simd_vshff(dtmpv4_5 , dtmpv4_1 , 0x88) ; \
  dp3d_trans_var_v4_3 = simd_vshff(dtmpv4_5 , dtmpv4_1 , 0xdd) ; \
  dp3d_trans_var_v4_4 = simd_vshff(dtmpv4_6 , dtmpv4_2 , 0x88) ; \
  dp3d_trans_var_v4_5 = simd_vshff(dtmpv4_6 , dtmpv4_2 , 0xdd) ; \
  dp3d_trans_var_v4_6 = simd_vshff(dtmpv4_7 , dtmpv4_3 , 0x88) ; \
  dp3d_trans_var_v4_7 = simd_vshff(dtmpv4_7 , dtmpv4_3 , 0xdd) ; \
  \
  dp3d_trans_local_v4_ptr[8] = dp3d_trans_var_v4_0 ; \
  dp3d_trans_local_v4_ptr[9] = dp3d_trans_var_v4_1 ; \
  dp3d_trans_local_v4_ptr[10] = dp3d_trans_var_v4_2 ; \
  dp3d_trans_local_v4_ptr[11] = dp3d_trans_var_v4_3 ; \
  dp3d_trans_local_v4_ptr[12] = dp3d_trans_var_v4_4 ; \
  dp3d_trans_local_v4_ptr[13] = dp3d_trans_var_v4_5 ; \
  dp3d_trans_local_v4_ptr[14] = dp3d_trans_var_v4_6 ; \
  dp3d_trans_local_v4_ptr[15] = dp3d_trans_var_v4_7 ; \
 }
// compute offset of send and recv buffer
#define COMPUTE_SEND_RECV_OFFSET(send_id_v8 , \
                                 send_offset_v8_dp3d , recv_offset_v8_dp3d , \
                                 send_offset_v8_dpn , recv_offset_v8_dpn , \
                                 send_offset_v8_forward_qdp , recv_offset_v8_forward_qdp , \
                                 send_offset_v8_post_qdp , recv_offset_v8_post_qdp , id) { \
  intv8 recv_index_v8 , send_index_v8 ; \
  intv8 recv_offset_v8 , send_offset_v8 ; \
  recv_index_v8 = simd_set_intv8(8 , 7 , 6 , 5 , 4 , 3 , 2 , 1) ; \
  send_index_v8 = simd_set_intv8(0 , 1 , 2 , 3 , 4 , 5 , 6 , 7) ; \
  VSHFW(send_id_v8 , id , id , 0x0) ; \
  send_offset_v8 = ((send_index_v8 + send_id_v8) & 0x7) ; \
  recv_offset_v8 = ((recv_index_v8 + send_id_v8) & 0x7) ; \
  send_id_v8 = send_offset_v8 ; \
  send_offset_v8_dp3d = send_offset_v8 << 1 ; \
  recv_offset_v8_dp3d = recv_offset_v8 ; \
  send_offset_v8_dpn = send_offset_v8 ; \
  recv_offset_v8_dpn =  recv_offset_v8 << 1 ; \
  send_offset_v8_forward_qdp = send_offset_v8 << 3 ; \
  recv_offset_v8_forward_qdp = recv_offset_v8 << 2 ; \
  send_offset_v8_post_qdp = send_offset_v8 << 2 ; \
  recv_offset_v8_post_qdp = recv_offset_v8 << 3 ; \
}
//zero commucation: copy module 
#define COPY_MEMERY_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset) { \
  doublev4 array_var_v4_0, array_var_v4_1 ;                            \
  \
  array_var_v4_0 = array_send_v4_ptr[send_offset] ; \
  array_var_v4_1 = array_send_v4_ptr[send_offset + 1] ; \
  \
  array_recv_v4_ptr[recv_offset + 1] = array_var_v4_0 ; \
  array_recv_v4_ptr[recv_offset + 10] = array_var_v4_1 ; \
}
// register ring commuication
#define REGISTER_RING_COMMUICATION_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) {\
  doublev4 array_var_v4_0, array_var_v4_1 ;                           \
  \
  array_var_v4_0 = array_send_v4_ptr[send_offset] ; \
  array_var_v4_1 = array_send_v4_ptr[send_offset + 1] ; \
  \
  REGISTER_PUTC(array_var_v4_0 , send_id) ; \
  REGISTER_PUTC(array_var_v4_1 , send_id) ; \
  \
  REGISTER_GETC(array_var_v4_0) ; \
  REGISTER_GETC(array_var_v4_1) ; \
  REGISTER_SYNC(0xff) ; \
  \
  array_recv_v4_ptr[recv_offset + 1] = array_var_v4_0 ; \
  array_recv_v4_ptr[recv_offset + 10] = array_var_v4_1 ; \
 \
}
// register full commuication
#define REGISTER_FULL_COMMUICATION_DP3D(array_recv_v4_ptr , array_send_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8) { \
  int send_id ; \
  int send_offset , recv_offset ; \
  \
  send_offset = simd_vextw0(send_offset_v8) ; \
  recv_offset = simd_vextw0(recv_offset_v8) ; \
  COPY_MEMERY_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset) ; \
  \
  send_id = simd_vextw1(send_id_v8) ; \
  send_offset = simd_vextw1(send_offset_v8) ; \
  recv_offset = simd_vextw1(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw2(send_id_v8) ; \
  send_offset = simd_vextw2(send_offset_v8) ; \
  recv_offset = simd_vextw2(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw3(send_id_v8) ; \
  send_offset = simd_vextw3(send_offset_v8) ; \
  recv_offset = simd_vextw3(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw4(send_id_v8) ; \
  send_offset = simd_vextw4(send_offset_v8) ; \
  recv_offset = simd_vextw4(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw5(send_id_v8) ; \
  send_offset = simd_vextw5(send_offset_v8) ; \
  recv_offset = simd_vextw5(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw6(send_id_v8) ; \
  send_offset = simd_vextw6(send_offset_v8) ; \
  recv_offset = simd_vextw6(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw7(send_id_v8) ; \
  send_offset = simd_vextw7(send_offset_v8) ; \
  recv_offset = simd_vextw7(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DP3D(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  doublev4 array_var_v4_8 , array_var_v4_17 ; \
  array_var_v4_8 = array_recv_v4_ptr[8] ; \
  array_var_v4_17 = array_recv_v4_ptr[17] ; \
  SET_ZERO_DOUBLEV4_HI(array_var_v4_8 , array_var_v4_8 , 0x44) ; \
  SET_ZERO_DOUBLEV4_HI(array_var_v4_17 , array_var_v4_17 , 0x44) ; \
  array_recv_v4_ptr[8] = array_var_v4_8 ; \
  array_recv_v4_ptr[17] = array_var_v4_17 ; \
}
// compute ps_v
#define COMPUTE_PS_V_SERIAL(ps_v_var_v4 , dpo_local_v4_ptr , hyai_ps0_var_v4) { \
  doublev4 dpo_var_v4_00 , dpo_var_v4_01 ; \
  doublev4 dpo_var_v4_10 , dpo_var_v4_11 ; \
  dpo_var_v4_00 = dpo_local_v4_ptr[1] ; \
  dpo_var_v4_01 = dpo_local_v4_ptr[2] ; \
  dpo_var_v4_10 = dpo_local_v4_ptr[10] ; \
  dpo_var_v4_11 = dpo_local_v4_ptr[11] ; \
  \
  doublev4 dpo_var_v4_0 , dpo_var_v4_1 , dpo_var_v4_2 , dpo_var_v4_3 ; \
  doublev4 dpo_var_v4_4 , dpo_var_v4_5 , dpo_var_v4_6 , dpo_var_v4_7 ; \
  VSHFF(dpo_var_v4_0 , dpo_var_v4_10 , dpo_var_v4_00 , 0x0) ; \
  VSHFF(dpo_var_v4_1 , dpo_var_v4_10 , dpo_var_v4_00 , 0x55) ; \
  VSHFF(dpo_var_v4_2 , dpo_var_v4_10 , dpo_var_v4_00 , 0xaa) ; \
  VSHFF(dpo_var_v4_3 , dpo_var_v4_10 , dpo_var_v4_00 , 0xff) ; \
  VSHFF(dpo_var_v4_4 , dpo_var_v4_11 , dpo_var_v4_01 , 0x0) ; \
  VSHFF(dpo_var_v4_5 , dpo_var_v4_11 , dpo_var_v4_01 , 0x55) ; \
  VSHFF(dpo_var_v4_6 , dpo_var_v4_11 , dpo_var_v4_01 , 0xaa) ; \
  VSHFF(dpo_var_v4_7 , dpo_var_v4_11 , dpo_var_v4_01 , 0xff) ; \
  \
  ps_v_var_v4 = dpo_var_v4_0 + dpo_var_v4_1 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_2 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_3 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_4 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_5 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_6 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_7 ; \
  \
  dpo_var_v4_00 = dpo_local_v4_ptr[3] ; \
  dpo_var_v4_01 = dpo_local_v4_ptr[4] ; \
  dpo_var_v4_10 = dpo_local_v4_ptr[12] ; \
  dpo_var_v4_11 = dpo_local_v4_ptr[13] ; \
  \
  VSHFF(dpo_var_v4_0 , dpo_var_v4_10 , dpo_var_v4_00 , 0x0) ; \
  VSHFF(dpo_var_v4_1 , dpo_var_v4_10 , dpo_var_v4_00 , 0x55) ; \
  VSHFF(dpo_var_v4_2 , dpo_var_v4_10 , dpo_var_v4_00 , 0xaa) ; \
  VSHFF(dpo_var_v4_3 , dpo_var_v4_10 , dpo_var_v4_00 , 0xff) ; \
  VSHFF(dpo_var_v4_4 , dpo_var_v4_11 , dpo_var_v4_01 , 0x0) ; \
  VSHFF(dpo_var_v4_5 , dpo_var_v4_11 , dpo_var_v4_01 , 0x55) ; \
  VSHFF(dpo_var_v4_6 , dpo_var_v4_11 , dpo_var_v4_01 , 0xaa) ; \
  VSHFF(dpo_var_v4_7 , dpo_var_v4_11 , dpo_var_v4_01 , 0xff) ; \
  \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_0 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_1 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_2 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_3 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_4 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_5 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_6 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_7 ; \
  \
  dpo_var_v4_00 = dpo_local_v4_ptr[5] ; \
  dpo_var_v4_01 = dpo_local_v4_ptr[6] ; \
  dpo_var_v4_10 = dpo_local_v4_ptr[14] ; \
  dpo_var_v4_11 = dpo_local_v4_ptr[15] ; \
  \
  VSHFF(dpo_var_v4_0 , dpo_var_v4_10 , dpo_var_v4_00 , 0x0) ; \
  VSHFF(dpo_var_v4_1 , dpo_var_v4_10 , dpo_var_v4_00 , 0x55) ; \
  VSHFF(dpo_var_v4_2 , dpo_var_v4_10 , dpo_var_v4_00 , 0xaa) ; \
  VSHFF(dpo_var_v4_3 , dpo_var_v4_10 , dpo_var_v4_00 , 0xff) ; \
  VSHFF(dpo_var_v4_4 , dpo_var_v4_11 , dpo_var_v4_01 , 0x0) ; \
  VSHFF(dpo_var_v4_5 , dpo_var_v4_11 , dpo_var_v4_01 , 0x55) ; \
  VSHFF(dpo_var_v4_6 , dpo_var_v4_11 , dpo_var_v4_01 , 0xaa) ; \
  VSHFF(dpo_var_v4_7 , dpo_var_v4_11 , dpo_var_v4_01 , 0xff) ; \
  \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_0 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_1 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_2 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_3 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_4 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_5 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_6 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_7 ; \
  \
  dpo_var_v4_00 = dpo_local_v4_ptr[7] ; \
  dpo_var_v4_01 = dpo_local_v4_ptr[8] ; \
  dpo_var_v4_10 = dpo_local_v4_ptr[16] ; \
  dpo_var_v4_11 = dpo_local_v4_ptr[17] ; \
  \
  VSHFF(dpo_var_v4_0 , dpo_var_v4_10 , dpo_var_v4_00 , 0x0) ; \
  VSHFF(dpo_var_v4_1 , dpo_var_v4_10 , dpo_var_v4_00 , 0x55) ; \
  VSHFF(dpo_var_v4_2 , dpo_var_v4_10 , dpo_var_v4_00 , 0xaa) ; \
  VSHFF(dpo_var_v4_3 , dpo_var_v4_10 , dpo_var_v4_00 , 0xff) ; \
  VSHFF(dpo_var_v4_4 , dpo_var_v4_11 , dpo_var_v4_01 , 0x0) ; \
  VSHFF(dpo_var_v4_5 , dpo_var_v4_11 , dpo_var_v4_01 , 0x55) ; \
  \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_0 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_1 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_2 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_3 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_4 ; \
  ps_v_var_v4 = ps_v_var_v4 + dpo_var_v4_5 ; \
  \
  VSHFF(ps_v_var_v4 , ps_v_var_v4 , ps_v_var_v4 , 0x88) ; \
  ps_v_var_v4 = ps_v_var_v4 + hyai_ps0_var_v4 ; \
}
// compute ps_v by vector
#define COMPUTE_PS_V(ps_v_var_v4 , dpo_local_v4_ptr , hyai_ps0_var_v4) { \
 doublev4 dpo_var_v4_0 , dpo_var_v4_1 , dpo_var_v4_2 , dpo_var_v4_3 ; \
 doublev4 dpo_var_v4_4 , dpo_var_v4_5 , dpo_var_v4_6 , dpo_var_v4_7 ; \
 doublev4 dpo_var_v4_8 , dpo_var_v4_9 , dpo_var_v4_10 , dpo_var_v4_11 ; \
 doublev4 dpo_var_v4_12 , dpo_var_v4_13 , dpo_var_v4_14 , dpo_var_v4_15 ; \
 \
 dpo_var_v4_0 = dpo_local_v4_ptr[1] ; \
 dpo_var_v4_1 = dpo_local_v4_ptr[2] ; \
 dpo_var_v4_2 = dpo_local_v4_ptr[3] ; \
 dpo_var_v4_3 = dpo_local_v4_ptr[4] ; \
 dpo_var_v4_4 = dpo_local_v4_ptr[5] ; \
 dpo_var_v4_5 = dpo_local_v4_ptr[6] ; \
 dpo_var_v4_6 = dpo_local_v4_ptr[7] ; \
 dpo_var_v4_7 = dpo_local_v4_ptr[8] ; \
 \
 dpo_var_v4_0 = dpo_var_v4_0 + dpo_var_v4_1 ; \
 dpo_var_v4_2 = dpo_var_v4_2 + dpo_var_v4_3 ; \
 dpo_var_v4_4 = dpo_var_v4_4 + dpo_var_v4_5 ; \
 dpo_var_v4_6 = dpo_var_v4_6 + dpo_var_v4_7 ; \
 \
 dpo_var_v4_8 = dpo_local_v4_ptr[10] ; \
 dpo_var_v4_9 = dpo_local_v4_ptr[11] ; \
 dpo_var_v4_10 = dpo_local_v4_ptr[12] ; \
 dpo_var_v4_11 = dpo_local_v4_ptr[13] ; \
 dpo_var_v4_12 = dpo_local_v4_ptr[14] ; \
 dpo_var_v4_13 = dpo_local_v4_ptr[15] ; \
 dpo_var_v4_14 = dpo_local_v4_ptr[16] ; \
 dpo_var_v4_15 = dpo_local_v4_ptr[17] ; \
 \
 dpo_var_v4_8 = dpo_var_v4_8 + dpo_var_v4_9 ; \
 dpo_var_v4_10 = dpo_var_v4_10 + dpo_var_v4_11 ; \
 dpo_var_v4_12 = dpo_var_v4_12 + dpo_var_v4_13 ; \
 dpo_var_v4_14 = dpo_var_v4_14 + dpo_var_v4_15 ; \
 \
 dpo_var_v4_0 = dpo_var_v4_0 + dpo_var_v4_2 ; \
 dpo_var_v4_4 = dpo_var_v4_4 + dpo_var_v4_6 ; \
 dpo_var_v4_8 = dpo_var_v4_8 + dpo_var_v4_10 ; \
 dpo_var_v4_12 = dpo_var_v4_12 + dpo_var_v4_14 ; \
 \
 dpo_var_v4_0 = dpo_var_v4_0 + dpo_var_v4_4 ; \
 dpo_var_v4_8 = dpo_var_v4_8 + dpo_var_v4_12 ; \
 \
 VSHFF(dpo_var_v4_1 , dpo_var_v4_8 , dpo_var_v4_0 , 0x44) ; \
 VSHFF(dpo_var_v4_11 , dpo_var_v4_8 , dpo_var_v4_0 , 0xee) ; \
 \
 dpo_var_v4_1 = dpo_var_v4_1 + dpo_var_v4_11 ; \
 \
 SET_ZERO_DOUBLEV4_HI(dpo_var_v4_2 , dpo_var_v4_1 , 0x88) ; \
 SET_ZERO_DOUBLEV4_HI(dpo_var_v4_12 , dpo_var_v4_1 , 0xdd) ; \
 \
 ps_v_var_v4 = dpo_var_v4_2 + dpo_var_v4_12 ; \
 ps_v_var_v4 = ps_v_var_v4 + hyai_ps0_var_v4 ; \
}

#define COMPUTE_DPN(dpn_local_v4_ptr , hyai_diff_local_v4_ptr , hybi_diff_local_v4_ptr , ps_v_var_v4) { \
  doublev4 hyai_diff_var_v4_0 , hybi_diff_var_v4_0 ; \
  doublev4 hyai_diff_var_v4_1 , hybi_diff_var_v4_1 ; \
  doublev4 hyai_diff_var_v4_2 , hybi_diff_var_v4_2 ; \
  doublev4 hyai_diff_var_v4_3 , hybi_diff_var_v4_3 ; \
  doublev4 hyai_diff_var_v4_4 , hybi_diff_var_v4_4 ; \
  doublev4 hyai_diff_var_v4_5 , hybi_diff_var_v4_5 ; \
  doublev4 hyai_diff_var_v4_6 , hybi_diff_var_v4_6 ; \
  doublev4 hyai_diff_var_v4_7 , hybi_diff_var_v4_7 ; \
  doublev4 ps_v_var_v4_0 , ps_v_var_v4_1 ; \
  \
  hyai_diff_var_v4_0 = hyai_diff_local_v4_ptr[0] ; \
  hyai_diff_var_v4_1 = hyai_diff_local_v4_ptr[1] ; \
  hyai_diff_var_v4_2 = hyai_diff_local_v4_ptr[2] ; \
  hyai_diff_var_v4_3 = hyai_diff_local_v4_ptr[3] ; \
  hyai_diff_var_v4_4 = hyai_diff_local_v4_ptr[4] ; \
  hyai_diff_var_v4_5 = hyai_diff_local_v4_ptr[5] ; \
  hyai_diff_var_v4_6 = hyai_diff_local_v4_ptr[6] ; \
  hyai_diff_var_v4_7 = hyai_diff_local_v4_ptr[7] ; \
  \
  VSHFF(ps_v_var_v4_0 , ps_v_var_v4 , ps_v_var_v4 , 0x0) ; \
  VSHFF(ps_v_var_v4_1 , ps_v_var_v4 , ps_v_var_v4 , 0x55) ; \
  \
  hybi_diff_var_v4_0 = hybi_diff_local_v4_ptr[0] ; \
  hybi_diff_var_v4_1 = hybi_diff_local_v4_ptr[1] ; \
  hybi_diff_var_v4_2 = hybi_diff_local_v4_ptr[2] ; \
  hybi_diff_var_v4_3 = hybi_diff_local_v4_ptr[3] ; \
  hybi_diff_var_v4_4 = hybi_diff_local_v4_ptr[4] ; \
  hybi_diff_var_v4_5 = hybi_diff_local_v4_ptr[5] ; \
  hybi_diff_var_v4_6 = hybi_diff_local_v4_ptr[6] ; \
  hybi_diff_var_v4_7 = hybi_diff_local_v4_ptr[7] ; \
  \
  doublev4 dpn_var_v4_0 , dpn_var_v4_1 , dpn_var_v4_2 , dpn_var_v4_3 ; \
  doublev4 dpn_var_v4_4 , dpn_var_v4_5 , dpn_var_v4_6 , dpn_var_v4_7 ; \
  \
  dpn_var_v4_0 = hyai_diff_var_v4_0 + hybi_diff_var_v4_0*ps_v_var_v4_0 ; \
  dpn_var_v4_1 = hyai_diff_var_v4_1 + hybi_diff_var_v4_1*ps_v_var_v4_0 ; \
  dpn_var_v4_2 = hyai_diff_var_v4_2 + hybi_diff_var_v4_2*ps_v_var_v4_0 ; \
  dpn_var_v4_3 = hyai_diff_var_v4_3 + hybi_diff_var_v4_3*ps_v_var_v4_0 ; \
  dpn_var_v4_4 = hyai_diff_var_v4_4 + hybi_diff_var_v4_4*ps_v_var_v4_0 ; \
  dpn_var_v4_5 = hyai_diff_var_v4_5 + hybi_diff_var_v4_5*ps_v_var_v4_0 ; \
  dpn_var_v4_6 = hyai_diff_var_v4_6 + hybi_diff_var_v4_6*ps_v_var_v4_0 ; \
  dpn_var_v4_7 = hyai_diff_var_v4_7 + hybi_diff_var_v4_7*ps_v_var_v4_0 ; \
  \
  dpn_local_v4_ptr[0] = dpn_var_v4_0 ; \
  dpn_local_v4_ptr[1] = dpn_var_v4_1 ; \
  dpn_local_v4_ptr[2] = dpn_var_v4_2 ; \
  dpn_local_v4_ptr[3] = dpn_var_v4_3 ; \
  dpn_local_v4_ptr[4] = dpn_var_v4_4 ; \
  dpn_local_v4_ptr[5] = dpn_var_v4_5 ; \
  dpn_local_v4_ptr[6] = dpn_var_v4_6 ; \
  dpn_local_v4_ptr[7] = dpn_var_v4_7 ; \
  \
  hyai_diff_var_v4_0 = hyai_diff_local_v4_ptr[0] ; \
  hyai_diff_var_v4_1 = hyai_diff_local_v4_ptr[1] ; \
  hyai_diff_var_v4_2 = hyai_diff_local_v4_ptr[2] ; \
  hyai_diff_var_v4_3 = hyai_diff_local_v4_ptr[3] ; \
  hyai_diff_var_v4_4 = hyai_diff_local_v4_ptr[4] ; \
  hyai_diff_var_v4_5 = hyai_diff_local_v4_ptr[5] ; \
  hyai_diff_var_v4_6 = hyai_diff_local_v4_ptr[6] ; \
  hyai_diff_var_v4_7 = hyai_diff_local_v4_ptr[7] ; \
  \
  hybi_diff_var_v4_0 = hybi_diff_local_v4_ptr[0] ; \
  hybi_diff_var_v4_1 = hybi_diff_local_v4_ptr[1] ; \
  hybi_diff_var_v4_2 = hybi_diff_local_v4_ptr[2] ; \
  hybi_diff_var_v4_3 = hybi_diff_local_v4_ptr[3] ; \
  hybi_diff_var_v4_4 = hybi_diff_local_v4_ptr[4] ; \
  hybi_diff_var_v4_5 = hybi_diff_local_v4_ptr[5] ; \
  hybi_diff_var_v4_6 = hybi_diff_local_v4_ptr[6] ; \
  hybi_diff_var_v4_7 = hybi_diff_local_v4_ptr[7] ; \
  \
  dpn_var_v4_0 = hyai_diff_var_v4_0 + hybi_diff_var_v4_0*ps_v_var_v4_1 ; \
  dpn_var_v4_1 = hyai_diff_var_v4_1 + hybi_diff_var_v4_1*ps_v_var_v4_1 ; \
  dpn_var_v4_2 = hyai_diff_var_v4_2 + hybi_diff_var_v4_2*ps_v_var_v4_1 ; \
  dpn_var_v4_3 = hyai_diff_var_v4_3 + hybi_diff_var_v4_3*ps_v_var_v4_1 ; \
  dpn_var_v4_4 = hyai_diff_var_v4_4 + hybi_diff_var_v4_4*ps_v_var_v4_1 ; \
  dpn_var_v4_5 = hyai_diff_var_v4_5 + hybi_diff_var_v4_5*ps_v_var_v4_1 ; \
  dpn_var_v4_6 = hyai_diff_var_v4_6 + hybi_diff_var_v4_6*ps_v_var_v4_1 ; \
  dpn_var_v4_7 = hyai_diff_var_v4_7 + hybi_diff_var_v4_7*ps_v_var_v4_1 ; \
  \
  dpn_local_v4_ptr[8] = dpn_var_v4_0 ; \
  dpn_local_v4_ptr[9] = dpn_var_v4_1 ; \
  dpn_local_v4_ptr[10] = dpn_var_v4_2 ; \
  dpn_local_v4_ptr[11] = dpn_var_v4_3 ; \
  dpn_local_v4_ptr[12] = dpn_var_v4_4 ; \
  dpn_local_v4_ptr[13] = dpn_var_v4_5 ; \
  dpn_local_v4_ptr[14] = dpn_var_v4_6 ; \
  dpn_local_v4_ptr[15] = dpn_var_v4_7 ; \
}
#define COMPUTE_PIO_PIN_2X2X2X4(pio_local_v4_ptr , pin_local_v4_ptr , dpo_local_v4_ptr , dpn_local_v4_ptr , pio_pin_var_v4_0 , ie) {\
  doublev4 dpo_var_v4_0 , dpn_var_v4_0 ; \
  doublev4 dpo_var_v4_1 , dpn_var_v4_1 ; \
  doublev4 dpo_var_v4_2 , dpn_var_v4_2 ; \
  doublev4 dpo_var_v4_3 , dpn_var_v4_3 ; \
  \
  dpo_var_v4_0 = dpo_local_v4_ptr[1] ; \
  dpo_var_v4_1 = dpo_local_v4_ptr[2] ; \
  dpo_var_v4_2 = dpo_local_v4_ptr[10] ; \
  dpo_var_v4_3 = dpo_local_v4_ptr[11] ; \
  dpn_var_v4_0 = dpn_local_v4_ptr[0] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[1] ; \
  dpn_var_v4_2 = dpn_local_v4_ptr[8] ; \
  dpn_var_v4_3 = dpn_local_v4_ptr[9] ; \
  \
  doublev4 dtmpv4_0 , dtmpv4_1 , dtmpv4_2 , dtmpv4_3 ; \
  doublev4 dtmpv4_4 , dtmpv4_5 , dtmpv4_6 , dtmpv4_7 ; \
  \
  dtmpv4_0 = simd_vshff(dpo_var_v4_2 , dpo_var_v4_0 , 0x44) ; \
  dtmpv4_1 = simd_vshff(dpo_var_v4_2 , dpo_var_v4_0 , 0xee) ; \
  dtmpv4_2 = simd_vshff(dpo_var_v4_3 , dpo_var_v4_1 , 0x44) ; \
  dtmpv4_3 = simd_vshff(dpo_var_v4_3 , dpo_var_v4_1 , 0xee) ; \
  dtmpv4_4 = simd_vshff(dpn_var_v4_2 , dpn_var_v4_0 , 0x44) ; \
  dtmpv4_5 = simd_vshff(dpn_var_v4_2 , dpn_var_v4_0 , 0xee) ; \
  dtmpv4_6 = simd_vshff(dpn_var_v4_3 , dpn_var_v4_1 , 0x44) ; \
  dtmpv4_7 = simd_vshff(dpn_var_v4_3 , dpn_var_v4_1 , 0xee) ; \
  \
  doublev4 pio_pin_var_v4_1 , pio_pin_var_v4_2 , pio_pin_var_v4_3 , pio_pin_var_v4_4 ; \
  doublev4 pio_pin_var_v4_5 , pio_pin_var_v4_6 , pio_pin_var_v4_7 , pio_pin_var_v4_8 ; \
  \
  pio_pin_var_v4_1 = simd_vshff(dtmpv4_4 , dtmpv4_0 , 0x88) ; \
  pio_pin_var_v4_2 = simd_vshff(dtmpv4_4 , dtmpv4_0 , 0xdd) ; \
  pio_pin_var_v4_3 = simd_vshff(dtmpv4_5 , dtmpv4_1 , 0x88) ; \
  pio_pin_var_v4_4 = simd_vshff(dtmpv4_5 , dtmpv4_1 , 0xdd) ; \
  pio_pin_var_v4_5 = simd_vshff(dtmpv4_6 , dtmpv4_2 , 0x88) ; \
  pio_pin_var_v4_6 = simd_vshff(dtmpv4_6 , dtmpv4_2 , 0xdd) ; \
  pio_pin_var_v4_7 = simd_vshff(dtmpv4_7 , dtmpv4_3 , 0x88) ; \
  pio_pin_var_v4_8 = simd_vshff(dtmpv4_7 , dtmpv4_3 , 0xdd) ; \
  \
  pio_pin_var_v4_1 = pio_pin_var_v4_0 + pio_pin_var_v4_1 ; \
  pio_pin_var_v4_2 = pio_pin_var_v4_1 + pio_pin_var_v4_2 ; \
  pio_pin_var_v4_3 = pio_pin_var_v4_2 + pio_pin_var_v4_3 ; \
  pio_pin_var_v4_4 = pio_pin_var_v4_3 + pio_pin_var_v4_4 ; \
  pio_pin_var_v4_5 = pio_pin_var_v4_4 + pio_pin_var_v4_5 ; \
  pio_pin_var_v4_6 = pio_pin_var_v4_5 + pio_pin_var_v4_6 ; \
  pio_pin_var_v4_7 = pio_pin_var_v4_6 + pio_pin_var_v4_7 ; \
  pio_pin_var_v4_8 = pio_pin_var_v4_7 + pio_pin_var_v4_8 ; \
  \
  dtmpv4_0 = simd_vshff(pio_pin_var_v4_1 , pio_pin_var_v4_0 , 0x44) ; \
  dtmpv4_1 = simd_vshff(pio_pin_var_v4_1 , pio_pin_var_v4_0 , 0xee) ; \
  dtmpv4_2 = simd_vshff(pio_pin_var_v4_3 , pio_pin_var_v4_2 , 0x44) ; \
  dtmpv4_3 = simd_vshff(pio_pin_var_v4_3 , pio_pin_var_v4_2 , 0xee) ; \
  dtmpv4_4 = simd_vshff(pio_pin_var_v4_5 , pio_pin_var_v4_4 , 0x44) ; \
  dtmpv4_5 = simd_vshff(pio_pin_var_v4_5 , pio_pin_var_v4_4 , 0xee) ; \
  dtmpv4_6 = simd_vshff(pio_pin_var_v4_7 , pio_pin_var_v4_6 , 0x44) ; \
  dtmpv4_7 = simd_vshff(pio_pin_var_v4_7 , pio_pin_var_v4_6 , 0xee) ; \
  \
  doublev4 pin_var_v4_0 , pin_var_v4_1 , pin_var_v4_2 , pin_var_v4_3 ; \
  doublev4 pio_var_v4_0 , pio_var_v4_1 , pio_var_v4_2 , pio_var_v4_3 ; \
  \
  pio_var_v4_0 = simd_vshff(dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  pio_var_v4_1 = simd_vshff(dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  pio_var_v4_2 = simd_vshff(dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  pio_var_v4_3 = simd_vshff(dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  pin_var_v4_0 = simd_vshff(dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  pin_var_v4_1 = simd_vshff(dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  pin_var_v4_2 = simd_vshff(dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  pin_var_v4_3 = simd_vshff(dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  pio_local_v4_ptr[0] = pio_var_v4_0 ; \
  pio_local_v4_ptr[1] = pio_var_v4_2 ; \
  pio_local_v4_ptr[8] = pio_var_v4_1 ; \
  pio_local_v4_ptr[9] = pio_var_v4_3 ; \
  pin_local_v4_ptr[0] = pin_var_v4_0 ; \
  pin_local_v4_ptr[1] = pin_var_v4_2 ; \
  pin_local_v4_ptr[8] = pin_var_v4_1 ; \
  pin_local_v4_ptr[9] = pin_var_v4_3 ; \
  \
  pio_pin_var_v4_0 = pio_pin_var_v4_8 ; \
}
#define COMPUTE_PIO_PIN(pio_local_v4_ptr , pin_local_v4_ptr , dpo_local_v4_ptr , dpn_local_v4_ptr , ie) { \
  doublev4 pio_pin_var_v4_0 = 0.0 ; \
  doublev4 *pio_local_v4_ptr_0 , *pio_local_v4_ptr_1 , *pio_local_v4_ptr_2 , *pio_local_v4_ptr_3 ; \
  doublev4 *pin_local_v4_ptr_0 , *pin_local_v4_ptr_1 , *pin_local_v4_ptr_2 , *pin_local_v4_ptr_3 ; \
  doublev4 *dpo_local_v4_ptr_0 , *dpo_local_v4_ptr_1 , *dpo_local_v4_ptr_2 , *dpo_local_v4_ptr_3 ; \
  doublev4 *dpn_local_v4_ptr_0 , *dpn_local_v4_ptr_1 , *dpn_local_v4_ptr_2 , *dpn_local_v4_ptr_3 ; \
  pio_local_v4_ptr_0 = &pio_local_v4_ptr[0] ; \
  pio_local_v4_ptr_1 = &pio_local_v4_ptr[2] ; \
  pio_local_v4_ptr_2 = &pio_local_v4_ptr[4] ; \
  pio_local_v4_ptr_3 = &pio_local_v4_ptr[6] ; \
  \
  pin_local_v4_ptr_0 = &pin_local_v4_ptr[0] ; \
  pin_local_v4_ptr_1 = &pin_local_v4_ptr[2] ; \
  pin_local_v4_ptr_2 = &pin_local_v4_ptr[4] ; \
  pin_local_v4_ptr_3 = &pin_local_v4_ptr[6] ; \
  \
  dpo_local_v4_ptr_0 = &dpo_local_v4_ptr[0] ; \
  dpo_local_v4_ptr_1 = &dpo_local_v4_ptr[2] ; \
  dpo_local_v4_ptr_2 = &dpo_local_v4_ptr[4] ; \
  dpo_local_v4_ptr_3 = &dpo_local_v4_ptr[6] ; \
  \
  dpn_local_v4_ptr_0 = &dpn_local_v4_ptr[0] ; \
  dpn_local_v4_ptr_1 = &dpn_local_v4_ptr[2] ; \
  dpn_local_v4_ptr_2 = &dpn_local_v4_ptr[4] ; \
  dpn_local_v4_ptr_3 = &dpn_local_v4_ptr[6] ; \
  \
  COMPUTE_PIO_PIN_2X2X2X4(pio_local_v4_ptr_0 , pin_local_v4_ptr_0 , dpo_local_v4_ptr_0 , dpn_local_v4_ptr_0 , pio_pin_var_v4_0 , ie) ; \
  COMPUTE_PIO_PIN_2X2X2X4(pio_local_v4_ptr_1 , pin_local_v4_ptr_1 , dpo_local_v4_ptr_1 , dpn_local_v4_ptr_1 , pio_pin_var_v4_0 , ie) ; \
  COMPUTE_PIO_PIN_2X2X2X4(pio_local_v4_ptr_2 , pin_local_v4_ptr_2 , dpo_local_v4_ptr_2 , dpn_local_v4_ptr_2 , pio_pin_var_v4_0 , ie) ; \
  COMPUTE_PIO_PIN_2X2X2X4(pio_local_v4_ptr_3 , pin_local_v4_ptr_3 , dpo_local_v4_ptr_3 , dpn_local_v4_ptr_3 , pio_pin_var_v4_0 , ie) ; \
  \
  double *pio_local_ptr , *pin_local_ptr ; \
  pio_local_ptr = (double *)pio_local_v4_ptr ; \
  pin_local_ptr = (double *)pin_local_v4_ptr ; \
  \
  pin_local_ptr[30] = pio_local_ptr[30] ; \
  pin_local_ptr[62] = pio_local_ptr[62] ; \
  pio_local_ptr[31] = pio_local_ptr[30] + 1.0 ; \
  pio_local_ptr[63] = pio_local_ptr[62] + 1.0 ; \
}
#define MEMCOPY_DP3D_2X36(dpo_local_v4_ptr) { \
  doublev4 dpo_var_v4_0 , dpo_var_v4_1 , dpo_var_v4_8 ; \
  doublev4 dpo_var_v4_9 , dpo_var_v4_10 , dpo_var_v4_17 ; \
  dpo_var_v4_1 = dpo_local_v4_ptr[1] ; \
  dpo_var_v4_8 = dpo_local_v4_ptr[8] ; \
  dpo_var_v4_10 = dpo_local_v4_ptr[10] ; \
  dpo_var_v4_17 = dpo_local_v4_ptr[17] ; \
  SET_ZERO_DOUBLEV4_LO(dpo_var_v4_0 , dpo_var_v4_1 , 0x10) ; \
  SET_ZERO_DOUBLEV4_LO(dpo_var_v4_9 , dpo_var_v4_10 , 0x10) ; \
  VSHFF(dpo_var_v4_8 , dpo_var_v4_8 , dpo_var_v4_8 , 0x14) ; \
  VSHFF(dpo_var_v4_17 , dpo_var_v4_17 , dpo_var_v4_17 , 0x14) ; \
  dpo_local_v4_ptr[0] = dpo_var_v4_0 ; \
  dpo_local_v4_ptr[8] = dpo_var_v4_8 ; \
  dpo_local_v4_ptr[9] = dpo_var_v4_9 ; \
  dpo_local_v4_ptr[17] = dpo_var_v4_17 ; \
}
#define COMPUTE_KID_Z(kid_local_ptr , z_local_ptr , pio_local_ptr , pin_local_ptr , dpo_local_ptr , ie ) { \
  \
  int k ; \
  for (k = 0 ; k < NLEV ; k = k + 2) { \
    double pin_var_0 , pin_var_1 , pin_var_2 , pin_var_3; \
    \
    pin_var_0 = pin_local_ptr[k+1] ; \
    pin_var_1 = pin_local_ptr[k+2] ; \
    pin_var_2 = pin_local_ptr[k+33] ; \
    pin_var_3 = pin_local_ptr[k+34] ; \
    \
    double pio_var_00 , pio_var_01 , pio_var_02 , pio_var_03 ; \
    int kk0 , kk1 , kk2 , kk3 ; \
    \
    for (kk0 = k ; kk0 < 32 ; kk0 ++) { \
      pio_var_00 = pio_local_ptr[kk0] ; \
      if (pio_var_00 > pin_var_0) { \
        kk0 = kk0 - 1 ; \
        break ; \
      } \
    }  \
    \
    for (kk1 = k+1 ; kk1 < 32 ; kk1 ++) { \
      pio_var_01 = pio_local_ptr[kk1] ; \
      if (pio_var_01 > pin_var_1) { \
        kk1 = kk1 - 1 ; \
        break ; \
      } \
    }  \
    \
    for (kk2 = k ; kk2 < 32 ; kk2 ++) { \
      pio_var_02 = pio_local_ptr[kk2 + 32] ; \
      if (pio_var_02 > pin_var_2) { \
        kk2 = kk2 - 1 ; \
        break ; \
      } \
    }  \
    \
    for (kk3 = k+1 ; kk3 < 32 ; kk3 ++) { \
      pio_var_03 = pio_local_ptr[kk3 + 32] ; \
      if (pio_var_03 > pin_var_3) { \
        kk3 = kk3 - 1 ; \
        break ; \
      } \
    }  \
    \
    kk0 = (kk0 == 30) ? 29 : kk0 ; \
    kk1 = (kk1 == 30) ? 29 : kk1 ; \
    kk2 = (kk2 == 30) ? 29 : kk2 ; \
    kk3 = (kk3 == 30) ? 29 : kk3 ; \
    \
    kid_local_ptr[k] = kk0 ; \
    kid_local_ptr[k + 1] = kk1 ; \
    kid_local_ptr[k + 32] = kk2 ; \
    kid_local_ptr[k + 33] = kk3 ; \
    \
    double pio_var_10 , pio_var_11 , pio_var_12 , pio_var_13 ; \
    \
    pio_var_00 = pio_local_ptr[kk0] ; \
    pio_var_10 = pio_local_ptr[kk0 + 1] ; \
    pio_var_01 = pio_local_ptr[kk1] ; \
    pio_var_11 = pio_local_ptr[kk1 + 1] ; \
    pio_var_02 = pio_local_ptr[kk2 + 32] ; \
    pio_var_12 = pio_local_ptr[kk2 + 33] ; \
    pio_var_03 = pio_local_ptr[kk3 + 32] ; \
    pio_var_13 = pio_local_ptr[kk3 + 33] ; \
    \
    double z_var_0 , z_var_1 , z_var_2 , z_var_3 ; \
    \
    z_var_0 = (pio_var_00 + pio_var_10)*0.5 ; \
    z_var_1 = (pio_var_01 + pio_var_11)*0.5 ; \
    z_var_2 = (pio_var_02 + pio_var_12)*0.5 ; \
    z_var_3 = (pio_var_03 + pio_var_13)*0.5 ; \
    \
    z_var_0 = pin_var_0 - z_var_0 ; \
    z_var_1 = pin_var_1 - z_var_1 ; \
    z_var_2 = pin_var_2 - z_var_2 ; \
    z_var_3 = pin_var_3 - z_var_3 ; \
    \
    double dpo_var_0 , dpo_var_1 , dpo_var_2 , dpo_var_3 ; \
    \
    dpo_var_0 = dpo_local_ptr[kk0 + 4] ; \
    dpo_var_1 = dpo_local_ptr[kk1 + 4] ; \
    dpo_var_2 = dpo_local_ptr[kk2 + 40] ; \
    dpo_var_3 = dpo_local_ptr[kk3 + 40] ; \
    \
    z_var_0 = z_var_0 / dpo_var_0 ; \
    z_var_1 = z_var_1 / dpo_var_1 ; \
    z_var_2 = z_var_2 / dpo_var_2 ; \
    z_var_3 = z_var_3 / dpo_var_3 ; \
    \
    z_local_ptr[k] = z_var_0 ; \
    z_local_ptr[k + 1] = z_var_1 ; \
    z_local_ptr[k + 32] = z_var_2 ; \
    z_local_ptr[k + 33] = z_var_3 ; \
    \
  } \
}
#define COMPUTE_PPM_GRIDS_2X2X4(ppmdx_local_v4_ptr , dpo_local_v4_ptr , dpo_var_v4_0_0 , condition ) { \
  doublev4 dpo_var_v4_1 , dpo_var_v4_2 , dpo_var_v4_3; \
  doublev4 dpo_var_v4_10 , dpo_var_v4_11 , dpo_var_v4_12 ; \
  dpo_var_v4_1 = dpo_local_v4_ptr[1] ; \
  dpo_var_v4_2 = dpo_local_v4_ptr[2] ; \
  dpo_var_v4_3 = dpo_local_v4_ptr[3] ; \
  dpo_var_v4_10 = dpo_local_v4_ptr[10] ; \
  dpo_var_v4_11 = dpo_local_v4_ptr[11] ; \
  dpo_var_v4_12 = dpo_local_v4_ptr[12] ; \
  \
  doublev4 dpo_var_v4_0_1 , dpo_var_v4_0_2 , dpo_var_v4_0_3 , dpo_var_v4_0_4 , dpo_var_v4_0_5; \
  VSHFF(dpo_var_v4_0_1 , dpo_var_v4_10 , dpo_var_v4_1 , 0x44) ; \
  VSHFF(dpo_var_v4_0_2 , dpo_var_v4_10 , dpo_var_v4_1 , 0xee) ; \
  VSHFF(dpo_var_v4_0_3 , dpo_var_v4_11 , dpo_var_v4_2 , 0x44) ; \
  VSHFF(dpo_var_v4_0_4 , dpo_var_v4_11 , dpo_var_v4_2 , 0xee) ; \
  VSHFF(dpo_var_v4_0_5 , dpo_var_v4_12 , dpo_var_v4_3 , 0x44) ; \
  \
  if(condition == 0 ) { \
    dpo_var_v4_0_5 = 0.0 ; \
  } \
  \
  doublev4 dpo_var_v4_1_0 , dpo_var_v4_1_1 , dpo_var_v4_1_2 , dpo_var_v4_1_3 , dpo_var_v4_1_4; \
  VSHFF(dpo_var_v4_1_0 , dpo_var_v4_0_1 , dpo_var_v4_0_0 , 0x8d) ; \
  VSHFF(dpo_var_v4_1_1 , dpo_var_v4_0_2 , dpo_var_v4_0_1 , 0x8d) ; \
  VSHFF(dpo_var_v4_1_2 , dpo_var_v4_0_3 , dpo_var_v4_0_2 , 0x8d) ; \
  VSHFF(dpo_var_v4_1_3 , dpo_var_v4_0_4 , dpo_var_v4_0_3 , 0x8d) ; \
  VSHFF(dpo_var_v4_1_4 , dpo_var_v4_0_5 , dpo_var_v4_0_4 , 0x8d) ; \
  \
  VSHFF(dpo_var_v4_1_0 , dpo_var_v4_1_0 , dpo_var_v4_1_0 , 0xd8) ; \
  VSHFF(dpo_var_v4_1_1 , dpo_var_v4_1_1 , dpo_var_v4_1_1 , 0xd8) ; \
  VSHFF(dpo_var_v4_1_2 , dpo_var_v4_1_2 , dpo_var_v4_1_2 , 0xd8) ; \
  VSHFF(dpo_var_v4_1_3 , dpo_var_v4_1_3 , dpo_var_v4_1_3 , 0xd8) ; \
  VSHFF(dpo_var_v4_1_4 , dpo_var_v4_1_4 , dpo_var_v4_1_4 , 0xd8) ; \
  \
  doublev4 dpo_sum_var_v4_12_0 , dpo_sum_var_v4_12_1 , dpo_sum_var_v4_12_2 , dpo_sum_var_v4_12_3 ; \
  dpo_sum_var_v4_12_0 = dpo_var_v4_1_0 + dpo_var_v4_0_1 ; \
  dpo_sum_var_v4_12_1 = dpo_var_v4_1_1 + dpo_var_v4_0_2 ; \
  dpo_sum_var_v4_12_2 = dpo_var_v4_1_2 + dpo_var_v4_0_3 ; \
  dpo_sum_var_v4_12_3 = dpo_var_v4_1_3 + dpo_var_v4_0_4 ; \
  \
  doublev4 dpo_sum_var_v4_01_0 , dpo_sum_var_v4_01_1 , dpo_sum_var_v4_01_2 , dpo_sum_var_v4_01_3 ; \
  dpo_sum_var_v4_01_0 = dpo_var_v4_0_0 + dpo_var_v4_1_0 ; \
  dpo_sum_var_v4_01_1 = dpo_var_v4_0_1 + dpo_var_v4_1_1 ; \
  dpo_sum_var_v4_01_2 = dpo_var_v4_0_2 + dpo_var_v4_1_2 ; \
  dpo_sum_var_v4_01_3 = dpo_var_v4_0_3 + dpo_var_v4_1_3 ; \
  \
  doublev4 ppmdx_var_v4_3_0 , ppmdx_var_v4_3_1 , ppmdx_var_v4_3_2 , ppmdx_var_v4_3_3 ; \
  ppmdx_var_v4_3_0 = dpo_var_v4_1_0 / dpo_sum_var_v4_12_0 ; \
  ppmdx_var_v4_3_1 = dpo_var_v4_1_1 / dpo_sum_var_v4_12_1 ; \
  ppmdx_var_v4_3_2 = dpo_var_v4_1_2 / dpo_sum_var_v4_12_2 ; \
  ppmdx_var_v4_3_3 = dpo_var_v4_1_3 / dpo_sum_var_v4_12_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_3_0 , ppmdx_var_v4_lo_3_1 ; \
  doublev4 ppmdx_var_v4_hi_3_0 , ppmdx_var_v4_hi_3_1 ; \
  VSHFF(ppmdx_var_v4_lo_3_0 , ppmdx_var_v4_3_1 , ppmdx_var_v4_3_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_3_0 , ppmdx_var_v4_3_1 , ppmdx_var_v4_3_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_3_1 , ppmdx_var_v4_3_3 , ppmdx_var_v4_3_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_3_1 , ppmdx_var_v4_3_3 , ppmdx_var_v4_3_2 , 0xee) ; \
  \
  doublev4 dpo_sum_var_v4_012_0 , dpo_sum_var_v4_012_1 , dpo_sum_var_v4_012_2 , dpo_sum_var_v4_012_3 ; \
  dpo_sum_var_v4_012_0 = dpo_sum_var_v4_01_0 + dpo_var_v4_0_1 ; \
  dpo_sum_var_v4_012_1 = dpo_sum_var_v4_01_1 + dpo_var_v4_0_2 ; \
  dpo_sum_var_v4_012_2 = dpo_sum_var_v4_01_2 + dpo_var_v4_0_3 ; \
  dpo_sum_var_v4_012_3 = dpo_sum_var_v4_01_3 + dpo_var_v4_0_4 ; \
  \
  ppmdx_local_v4_ptr[48] = ppmdx_var_v4_lo_3_0 ; \
  ppmdx_local_v4_ptr[49] = ppmdx_var_v4_lo_3_1 ; \
  ppmdx_local_v4_ptr[56] = ppmdx_var_v4_hi_3_0 ; \
  ppmdx_local_v4_ptr[57] = ppmdx_var_v4_hi_3_1 ; \
  \
  doublev4 ppmdx_var_v4_0_0 , ppmdx_var_v4_0_1 , ppmdx_var_v4_0_2 , ppmdx_var_v4_0_3 ; \
  ppmdx_var_v4_0_0 = dpo_var_v4_1_0 / dpo_sum_var_v4_012_0 ; \
  ppmdx_var_v4_0_1 = dpo_var_v4_1_1 / dpo_sum_var_v4_012_1 ; \
  ppmdx_var_v4_0_2 = dpo_var_v4_1_2 / dpo_sum_var_v4_012_2 ; \
  ppmdx_var_v4_0_3 = dpo_var_v4_1_3 / dpo_sum_var_v4_012_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_0_0 , ppmdx_var_v4_lo_0_1 ; \
  doublev4 ppmdx_var_v4_hi_0_0 , ppmdx_var_v4_hi_0_1 ; \
  VSHFF(ppmdx_var_v4_lo_0_0 , ppmdx_var_v4_0_1 , ppmdx_var_v4_0_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_0_0 , ppmdx_var_v4_0_1 , ppmdx_var_v4_0_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_0_1 , ppmdx_var_v4_0_3 , ppmdx_var_v4_0_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_0_1 , ppmdx_var_v4_0_3 , ppmdx_var_v4_0_2 , 0xee) ; \
  \
  doublev4 dpo_sum_var_v4_0123_0 , dpo_sum_var_v4_0123_1 , dpo_sum_var_v4_0123_2 , dpo_sum_var_v4_0123_3 ; \
  dpo_sum_var_v4_0123_0 = dpo_sum_var_v4_012_0 + dpo_var_v4_1_1 ; \
  dpo_sum_var_v4_0123_1 = dpo_sum_var_v4_012_1 + dpo_var_v4_1_2 ; \
  dpo_sum_var_v4_0123_2 = dpo_sum_var_v4_012_2 + dpo_var_v4_1_3 ; \
  dpo_sum_var_v4_0123_3 = dpo_sum_var_v4_012_3 + dpo_var_v4_1_4 ; \
  \
  ppmdx_local_v4_ptr[0] = ppmdx_var_v4_lo_0_0 ; \
  ppmdx_local_v4_ptr[1] = ppmdx_var_v4_lo_0_1 ; \
  ppmdx_local_v4_ptr[8] = ppmdx_var_v4_hi_0_0 ; \
  ppmdx_local_v4_ptr[9] = ppmdx_var_v4_hi_0_1 ; \
  \
  doublev4 ppmdx_var_v4_4_0 , ppmdx_var_v4_4_1 , ppmdx_var_v4_4_2 , ppmdx_var_v4_4_3 ; \
  ppmdx_var_v4_4_0 = 1.0 / dpo_sum_var_v4_0123_0 ; \
  ppmdx_var_v4_4_1 = 1.0 / dpo_sum_var_v4_0123_1 ; \
  ppmdx_var_v4_4_2 = 1.0 / dpo_sum_var_v4_0123_2 ; \
  ppmdx_var_v4_4_3 = 1.0 / dpo_sum_var_v4_0123_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_4_0 , ppmdx_var_v4_lo_4_1 ; \
  doublev4 ppmdx_var_v4_hi_4_0 , ppmdx_var_v4_hi_4_1 ; \
  VSHFF(ppmdx_var_v4_lo_4_0 , ppmdx_var_v4_4_1 , ppmdx_var_v4_4_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_4_0 , ppmdx_var_v4_4_1 , ppmdx_var_v4_4_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_4_1 , ppmdx_var_v4_4_3 , ppmdx_var_v4_4_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_4_1 , ppmdx_var_v4_4_3 , ppmdx_var_v4_4_2 , 0xee) ; \
  \
  doublev4 dpo_sum_var_v4_001_0 ,  dpo_sum_var_v4_001_1 ,  dpo_sum_var_v4_001_2 ,  dpo_sum_var_v4_001_3 ; \
  dpo_sum_var_v4_001_0 = dpo_sum_var_v4_01_0 + dpo_var_v4_0_0 ; \
  dpo_sum_var_v4_001_1 = dpo_sum_var_v4_01_1 + dpo_var_v4_0_1 ; \
  dpo_sum_var_v4_001_2 = dpo_sum_var_v4_01_2 + dpo_var_v4_0_2 ; \
  dpo_sum_var_v4_001_3 = dpo_sum_var_v4_01_3 + dpo_var_v4_0_3 ; \
  \
  ppmdx_local_v4_ptr[64] = ppmdx_var_v4_lo_4_0 ; \
  ppmdx_local_v4_ptr[65] = ppmdx_var_v4_lo_4_1 ; \
  ppmdx_local_v4_ptr[72] = ppmdx_var_v4_hi_4_0 ; \
  ppmdx_local_v4_ptr[73] = ppmdx_var_v4_hi_4_1 ; \
  \
  doublev4 ppmdx_var_v4_1_0 , ppmdx_var_v4_1_1 , ppmdx_var_v4_1_2 , ppmdx_var_v4_1_3 ; \
  ppmdx_var_v4_1_0 = dpo_sum_var_v4_001_0 / dpo_sum_var_v4_12_0 ; \
  ppmdx_var_v4_1_1 = dpo_sum_var_v4_001_1 / dpo_sum_var_v4_12_1 ; \
  ppmdx_var_v4_1_2 = dpo_sum_var_v4_001_2 / dpo_sum_var_v4_12_2 ; \
  ppmdx_var_v4_1_3 = dpo_sum_var_v4_001_3 / dpo_sum_var_v4_12_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_1_0 , ppmdx_var_v4_lo_1_1 ; \
  doublev4 ppmdx_var_v4_hi_1_0 , ppmdx_var_v4_hi_1_1 ; \
  VSHFF(ppmdx_var_v4_lo_1_0 , ppmdx_var_v4_1_1 , ppmdx_var_v4_1_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_1_0 , ppmdx_var_v4_1_1 , ppmdx_var_v4_1_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_1_1 , ppmdx_var_v4_1_3 , ppmdx_var_v4_1_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_1_1 , ppmdx_var_v4_1_3 , ppmdx_var_v4_1_2 , 0xee) ; \
  \
  doublev4 dpo_sum_var_v4_112_0 , dpo_sum_var_v4_112_1 , dpo_sum_var_v4_112_2 , dpo_sum_var_v4_112_3 ; \
  dpo_sum_var_v4_112_0 = dpo_sum_var_v4_12_0 + dpo_var_v4_1_0 ; \
  dpo_sum_var_v4_112_1 = dpo_sum_var_v4_12_1 + dpo_var_v4_1_1 ; \
  dpo_sum_var_v4_112_2 = dpo_sum_var_v4_12_2 + dpo_var_v4_1_2 ; \
  dpo_sum_var_v4_112_3 = dpo_sum_var_v4_12_3 + dpo_var_v4_1_3 ; \
  \
  ppmdx_local_v4_ptr[16] = ppmdx_var_v4_lo_1_0 ; \
  ppmdx_local_v4_ptr[17] = ppmdx_var_v4_lo_1_1 ; \
  ppmdx_local_v4_ptr[24] = ppmdx_var_v4_hi_1_0 ; \
  ppmdx_local_v4_ptr[25] = ppmdx_var_v4_hi_1_1 ; \
  \
  doublev4 ppmdx_var_v4_6_0 , ppmdx_var_v4_6_1 , ppmdx_var_v4_6_2 , ppmdx_var_v4_6_3 ; \
  ppmdx_var_v4_6_0 = dpo_sum_var_v4_01_0 / dpo_sum_var_v4_112_0 ; \
  ppmdx_var_v4_6_1 = dpo_sum_var_v4_01_1 / dpo_sum_var_v4_112_1 ; \
  ppmdx_var_v4_6_2 = dpo_sum_var_v4_01_2 / dpo_sum_var_v4_112_2 ; \
  ppmdx_var_v4_6_3 = dpo_sum_var_v4_01_3 / dpo_sum_var_v4_112_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_6_0 , ppmdx_var_v4_lo_6_1 ; \
  doublev4 ppmdx_var_v4_hi_6_0 , ppmdx_var_v4_hi_6_1 ; \
  VSHFF(ppmdx_var_v4_lo_6_0 , ppmdx_var_v4_6_1 , ppmdx_var_v4_6_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_6_0 , ppmdx_var_v4_6_1 , ppmdx_var_v4_6_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_6_1 , ppmdx_var_v4_6_3 , ppmdx_var_v4_6_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_6_1 , ppmdx_var_v4_6_3 , ppmdx_var_v4_6_2 , 0xee) ; \
  \
  doublev4 ppmdx_var_v4_8_0 , ppmdx_var_v4_8_1 , ppmdx_var_v4_8_2 , ppmdx_var_v4_8_3 ; \
  ppmdx_var_v4_8_0 = dpo_var_v4_1_0 * ppmdx_var_v4_6_0 ; \
  ppmdx_var_v4_8_1 = dpo_var_v4_1_1 * ppmdx_var_v4_6_1 ; \
  ppmdx_var_v4_8_2 = dpo_var_v4_1_2 * ppmdx_var_v4_6_2 ; \
  ppmdx_var_v4_8_3 = dpo_var_v4_1_3 * ppmdx_var_v4_6_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_8_0 , ppmdx_var_v4_lo_8_1 ; \
  doublev4 ppmdx_var_v4_hi_8_0 , ppmdx_var_v4_hi_8_1 ; \
  VSHFF(ppmdx_var_v4_lo_8_0 , ppmdx_var_v4_8_1 , ppmdx_var_v4_8_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_8_0 , ppmdx_var_v4_8_1 , ppmdx_var_v4_8_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_8_1 , ppmdx_var_v4_8_3 , ppmdx_var_v4_8_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_8_1 , ppmdx_var_v4_8_3 , ppmdx_var_v4_8_2 , 0xee) ; \
  \
  ppmdx_local_v4_ptr[96] = ppmdx_var_v4_lo_6_0 ; \
  ppmdx_local_v4_ptr[97] = ppmdx_var_v4_lo_6_1 ; \
  ppmdx_local_v4_ptr[104] = ppmdx_var_v4_hi_6_0 ; \
  ppmdx_local_v4_ptr[105] = ppmdx_var_v4_hi_6_1 ; \
  \
  doublev4 dpo_mul_var_v4_12_0 , dpo_mul_var_v4_12_1 , dpo_mul_var_v4_12_2 , dpo_mul_var_v4_12_3 ; \
  dpo_mul_var_v4_12_0 = dpo_var_v4_1_0 * dpo_var_v4_0_1 ; \
  dpo_mul_var_v4_12_1 = dpo_var_v4_1_1 * dpo_var_v4_0_2 ; \
  dpo_mul_var_v4_12_2 = dpo_var_v4_1_2 * dpo_var_v4_0_3 ; \
  dpo_mul_var_v4_12_3 = dpo_var_v4_1_3 * dpo_var_v4_0_4 ; \
  \
  ppmdx_local_v4_ptr[128] = ppmdx_var_v4_lo_8_0 ; \
  ppmdx_local_v4_ptr[129] = ppmdx_var_v4_lo_8_1 ; \
  ppmdx_local_v4_ptr[136] = ppmdx_var_v4_hi_8_0 ; \
  ppmdx_local_v4_ptr[137] = ppmdx_var_v4_hi_8_1 ; \
  \
  dpo_mul_var_v4_12_0 = dpo_mul_var_v4_12_0 * 2.0 ; \
  dpo_mul_var_v4_12_1 = dpo_mul_var_v4_12_1 * 2.0 ; \
  dpo_mul_var_v4_12_2 = dpo_mul_var_v4_12_2 * 2.0 ; \
  dpo_mul_var_v4_12_3 = dpo_mul_var_v4_12_3 * 2.0 ; \
  \
  doublev4 ppmdx_var_v4_5_0 , ppmdx_var_v4_5_1 , ppmdx_var_v4_5_2 , ppmdx_var_v4_5_3 ; \
  ppmdx_var_v4_5_0 = dpo_mul_var_v4_12_0 / dpo_sum_var_v4_12_0 ; \
  ppmdx_var_v4_5_1 = dpo_mul_var_v4_12_1 / dpo_sum_var_v4_12_1 ; \
  ppmdx_var_v4_5_2 = dpo_mul_var_v4_12_2 / dpo_sum_var_v4_12_2 ; \
  ppmdx_var_v4_5_3 = dpo_mul_var_v4_12_3 / dpo_sum_var_v4_12_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_5_0 , ppmdx_var_v4_lo_5_1 ; \
  doublev4 ppmdx_var_v4_hi_5_0 , ppmdx_var_v4_hi_5_1 ; \
  VSHFF(ppmdx_var_v4_lo_5_0 , ppmdx_var_v4_5_1 , ppmdx_var_v4_5_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_5_0 , ppmdx_var_v4_5_1 , ppmdx_var_v4_5_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_5_1 , ppmdx_var_v4_5_3 , ppmdx_var_v4_5_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_5_1 , ppmdx_var_v4_5_3 , ppmdx_var_v4_5_2 , 0xee) ; \
  \
  doublev4 dpo_sum_var_v4_122_0 , dpo_sum_var_v4_122_1 , dpo_sum_var_v4_122_2 , dpo_sum_var_v4_122_3 ; \
  dpo_sum_var_v4_122_0 = dpo_sum_var_v4_12_0 + dpo_var_v4_0_1 ; \
  dpo_sum_var_v4_122_1 = dpo_sum_var_v4_12_1 + dpo_var_v4_0_2 ; \
  dpo_sum_var_v4_122_2 = dpo_sum_var_v4_12_2 + dpo_var_v4_0_3 ; \
  dpo_sum_var_v4_122_3 = dpo_sum_var_v4_12_3 + dpo_var_v4_0_4 ; \
  \
  ppmdx_local_v4_ptr[80] = ppmdx_var_v4_lo_5_0 ; \
  ppmdx_local_v4_ptr[81] = ppmdx_var_v4_lo_5_1 ; \
  ppmdx_local_v4_ptr[88] = ppmdx_var_v4_hi_5_0 ; \
  ppmdx_local_v4_ptr[89] = ppmdx_var_v4_hi_5_1 ; \
  \
  doublev4 ppmdx_var_v4_2_0 , ppmdx_var_v4_2_1 , ppmdx_var_v4_2_2 , ppmdx_var_v4_2_3 ; \
  ppmdx_var_v4_2_0 = dpo_sum_var_v4_122_0 / dpo_sum_var_v4_01_0 ; \
  ppmdx_var_v4_2_1 = dpo_sum_var_v4_122_1 / dpo_sum_var_v4_01_1 ; \
  ppmdx_var_v4_2_2 = dpo_sum_var_v4_122_2 / dpo_sum_var_v4_01_2 ; \
  ppmdx_var_v4_2_3 = dpo_sum_var_v4_122_3 / dpo_sum_var_v4_01_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_2_0 , ppmdx_var_v4_lo_2_1 ; \
  doublev4 ppmdx_var_v4_hi_2_0 , ppmdx_var_v4_hi_2_1 ; \
  VSHFF(ppmdx_var_v4_lo_2_0 , ppmdx_var_v4_2_1 , ppmdx_var_v4_2_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_2_0 , ppmdx_var_v4_2_1 , ppmdx_var_v4_2_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_2_1 , ppmdx_var_v4_2_3 , ppmdx_var_v4_2_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_2_1 , ppmdx_var_v4_2_3 , ppmdx_var_v4_2_2 , 0xee) ; \
  \
  doublev4 dpo_sum_var_v4_23_0 , dpo_sum_var_v4_23_1 , dpo_sum_var_v4_23_2 , dpo_sum_var_v4_23_3; \
  dpo_sum_var_v4_23_0 = dpo_var_v4_0_1 + dpo_var_v4_1_1 ; \
  dpo_sum_var_v4_23_1 = dpo_var_v4_0_2 + dpo_var_v4_1_2 ; \
  dpo_sum_var_v4_23_2 = dpo_var_v4_0_3 + dpo_var_v4_1_3 ; \
  dpo_sum_var_v4_23_3 = dpo_var_v4_0_4 + dpo_var_v4_1_4 ; \
  \
  ppmdx_local_v4_ptr[32] = ppmdx_var_v4_lo_2_0 ; \
  ppmdx_local_v4_ptr[33] = ppmdx_var_v4_lo_2_1 ; \
  ppmdx_local_v4_ptr[40] = ppmdx_var_v4_hi_2_0 ; \
  ppmdx_local_v4_ptr[41] = ppmdx_var_v4_hi_2_1 ; \
  \
  doublev4 ppmdx_var_v4_7_0 , ppmdx_var_v4_7_1 , ppmdx_var_v4_7_2, ppmdx_var_v4_7_3 ; \
  ppmdx_var_v4_7_0 = dpo_sum_var_v4_23_0 / dpo_sum_var_v4_122_0 ; \
  ppmdx_var_v4_7_1 = dpo_sum_var_v4_23_1 / dpo_sum_var_v4_122_1 ; \
  ppmdx_var_v4_7_2 = dpo_sum_var_v4_23_2 / dpo_sum_var_v4_122_2 ; \
  ppmdx_var_v4_7_3 = dpo_sum_var_v4_23_3 / dpo_sum_var_v4_122_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_7_0 , ppmdx_var_v4_lo_7_1 ; \
  doublev4 ppmdx_var_v4_hi_7_0 , ppmdx_var_v4_hi_7_1 ; \
  VSHFF(ppmdx_var_v4_lo_7_0 , ppmdx_var_v4_7_1 , ppmdx_var_v4_7_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_7_0 , ppmdx_var_v4_7_1 , ppmdx_var_v4_7_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_7_1 , ppmdx_var_v4_7_3 , ppmdx_var_v4_7_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_7_1 , ppmdx_var_v4_7_3 , ppmdx_var_v4_7_2 , 0xee) ; \
  \
  ppmdx_local_v4_ptr[112] = ppmdx_var_v4_lo_7_0 ; \
  ppmdx_local_v4_ptr[113] = ppmdx_var_v4_lo_7_1 ; \
  ppmdx_local_v4_ptr[120] = ppmdx_var_v4_hi_7_0 ; \
  ppmdx_local_v4_ptr[121] = ppmdx_var_v4_hi_7_1 ; \
  \
  doublev4 ppmdx_var_v4_9_0 , ppmdx_var_v4_9_1 , ppmdx_var_v4_9_2 , ppmdx_var_v4_9_3 ; \
  ppmdx_var_v4_9_0 = dpo_var_v4_0_1 * ppmdx_var_v4_7_0 ; \
  ppmdx_var_v4_9_1 = dpo_var_v4_0_2 * ppmdx_var_v4_7_1 ; \
  ppmdx_var_v4_9_2 = dpo_var_v4_0_3 * ppmdx_var_v4_7_2 ; \
  ppmdx_var_v4_9_3 = dpo_var_v4_0_4 * ppmdx_var_v4_7_3 ; \
  \
  doublev4 ppmdx_var_v4_lo_9_0 , ppmdx_var_v4_lo_9_1 ; \
  doublev4 ppmdx_var_v4_hi_9_0 , ppmdx_var_v4_hi_9_1 ; \
  VSHFF(ppmdx_var_v4_lo_9_0 , ppmdx_var_v4_9_1 , ppmdx_var_v4_9_0 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_9_0 , ppmdx_var_v4_9_1 , ppmdx_var_v4_9_0 , 0xee) ; \
  VSHFF(ppmdx_var_v4_lo_9_1 , ppmdx_var_v4_9_3 , ppmdx_var_v4_9_2 , 0x44) ; \
  VSHFF(ppmdx_var_v4_hi_9_1 , ppmdx_var_v4_9_3 , ppmdx_var_v4_9_2 , 0xee) ; \
  \
  ppmdx_local_v4_ptr[144] = ppmdx_var_v4_lo_9_0 ; \
  ppmdx_local_v4_ptr[145] = ppmdx_var_v4_lo_9_1 ; \
  ppmdx_local_v4_ptr[152] = ppmdx_var_v4_hi_9_0 ; \
  ppmdx_local_v4_ptr[153] = ppmdx_var_v4_hi_9_1 ; \
  \
  dpo_var_v4_0_0 = dpo_var_v4_0_4 ; \
}  
#define COMPUTE_PPM_GRIDS(ppmdx_local_v4_ptr , dpo_local_v4_ptr) { \
  doublev4 dpo_var_v4_0 , dpo_var_v4_9 ; \
  dpo_var_v4_0 = dpo_local_v4_ptr[0] ; \
  dpo_var_v4_9 = dpo_local_v4_ptr[9] ; \
  \
  doublev4 dpo_var_v4_0_0 ; \
  VSHFF(dpo_var_v4_0_0 , dpo_var_v4_9 , dpo_var_v4_0 , 0xee) ; \
  \
  doublev4 *ppmdx_local_v4_ptr_0 , *ppmdx_local_v4_ptr_1 , *ppmdx_local_v4_ptr_2 , *ppmdx_local_v4_ptr_3 ; \
  doublev4 *dpo_local_v4_ptr_0 , *dpo_local_v4_ptr_1 , *dpo_local_v4_ptr_2 , *dpo_local_v4_ptr_3 ;\
  ppmdx_local_v4_ptr_0 = &ppmdx_local_v4_ptr[0] ; \
  ppmdx_local_v4_ptr_1 = &ppmdx_local_v4_ptr[2] ; \
  ppmdx_local_v4_ptr_2 = &ppmdx_local_v4_ptr[4] ; \
  ppmdx_local_v4_ptr_3 = &ppmdx_local_v4_ptr[6] ; \
  \
  dpo_local_v4_ptr_0 = &dpo_local_v4_ptr[0] ; \
  dpo_local_v4_ptr_1 = &dpo_local_v4_ptr[2] ; \
  dpo_local_v4_ptr_2 = &dpo_local_v4_ptr[4] ; \
  dpo_local_v4_ptr_3 = &dpo_local_v4_ptr[6] ; \
  \
  COMPUTE_PPM_GRIDS_2X2X4(ppmdx_local_v4_ptr_0 , dpo_local_v4_ptr_0 , dpo_var_v4_0_0 , 1) ; \
  COMPUTE_PPM_GRIDS_2X2X4(ppmdx_local_v4_ptr_1 , dpo_local_v4_ptr_1 , dpo_var_v4_0_0 , 1) ; \
  COMPUTE_PPM_GRIDS_2X2X4(ppmdx_local_v4_ptr_2 , dpo_local_v4_ptr_2 , dpo_var_v4_0_0 , 1) ; \
  COMPUTE_PPM_GRIDS_2X2X4(ppmdx_local_v4_ptr_3 , dpo_local_v4_ptr_3 , dpo_var_v4_0_0 , 0) ; \
}

#define COPY_MEMERY_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset) { \
  doublev4 array_var_v4_0, array_var_v4_1 ;                            \
  \
  array_var_v4_0 = array_send_v4_ptr[send_offset] ; \
  array_var_v4_1 = array_send_v4_ptr[send_offset + 8] ; \
  \
  array_recv_v4_ptr[recv_offset] = array_var_v4_0 ; \
  array_recv_v4_ptr[recv_offset + 1] = array_var_v4_1 ; \
  \
}


#define REGISTER_RING_COMMUICATION_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) {\
  doublev4 array_var_v4_0, array_var_v4_1 ;                           \
  \
  array_var_v4_0 = array_send_v4_ptr[send_offset] ; \
  array_var_v4_1 = array_send_v4_ptr[send_offset + 8] ; \
  \
  REGISTER_PUTC(array_var_v4_0 , send_id) ; \
  REGISTER_PUTC(array_var_v4_1 , send_id) ; \
  \
  REGISTER_GETC(array_var_v4_0) ; \
  REGISTER_GETC(array_var_v4_1) ; \
  REGISTER_SYNC(0xff) ; \
  \
  array_recv_v4_ptr[recv_offset] = array_var_v4_0 ; \
  array_recv_v4_ptr[recv_offset + 1] = array_var_v4_1 ; \
}

#define REGISTER_FULL_COMMUICATION_DPN(array_recv_v4_ptr , array_send_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8) { \
  int send_id ; \
  int send_offset , recv_offset ; \
  \
  send_offset = simd_vextw0(send_offset_v8) ; \
  recv_offset = simd_vextw0(recv_offset_v8) ; \
  COPY_MEMERY_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset) ; \
  \
  send_id = simd_vextw1(send_id_v8) ; \
  send_offset = simd_vextw1(send_offset_v8) ; \
  recv_offset = simd_vextw1(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw2(send_id_v8) ; \
  send_offset = simd_vextw2(send_offset_v8) ; \
  recv_offset = simd_vextw2(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw3(send_id_v8) ; \
  send_offset = simd_vextw3(send_offset_v8) ; \
  recv_offset = simd_vextw3(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw4(send_id_v8) ; \
  send_offset = simd_vextw4(send_offset_v8) ; \
  recv_offset = simd_vextw4(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw5(send_id_v8) ; \
  send_offset = simd_vextw5(send_offset_v8) ; \
  recv_offset = simd_vextw5(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw6(send_id_v8) ; \
  send_offset = simd_vextw6(send_offset_v8) ; \
  recv_offset = simd_vextw6(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw7(send_id_v8) ; \
  send_offset = simd_vextw7(send_offset_v8) ; \
  recv_offset = simd_vextw7(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_DPN(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
}

 #define TRANSPOSE_ARRAY_DPN(dpn_local_v4_ptr , dpn_trans_local_v4_ptr) { \
  doublev4 dpn_trans_var_v4_0 , dpn_trans_var_v4_1 , dpn_trans_var_v4_2 , dpn_trans_var_v4_3 ;     \
  doublev4 dpn_trans_var_v4_4 , dpn_trans_var_v4_5 , dpn_trans_var_v4_6 , dpn_trans_var_v4_7 ;     \
  \
  dpn_trans_var_v4_0 = dpn_trans_local_v4_ptr[0] ;  \
  dpn_trans_var_v4_1 = dpn_trans_local_v4_ptr[1] ;  \
  dpn_trans_var_v4_2 = dpn_trans_local_v4_ptr[2] ;  \
  dpn_trans_var_v4_3 = dpn_trans_local_v4_ptr[3] ;  \
  dpn_trans_var_v4_4 = dpn_trans_local_v4_ptr[4] ;  \
  dpn_trans_var_v4_5 = dpn_trans_local_v4_ptr[5] ;  \
  dpn_trans_var_v4_6 = dpn_trans_local_v4_ptr[6] ; \
  dpn_trans_var_v4_7 = dpn_trans_local_v4_ptr[7] ; \
  \
  doublev4 dtmpv4_0 , dtmpv4_1 , dtmpv4_2 , dtmpv4_3 ; \
  doublev4 dtmpv4_4 , dtmpv4_5 , dtmpv4_6 , dtmpv4_7 ; \
  \
  dtmpv4_0 = simd_vshff(dpn_trans_var_v4_1 , dpn_trans_var_v4_0 , 0x44) ; \
  dtmpv4_1 = simd_vshff(dpn_trans_var_v4_1 , dpn_trans_var_v4_0 , 0xee) ; \
  dtmpv4_2 = simd_vshff(dpn_trans_var_v4_3 , dpn_trans_var_v4_2 , 0x44) ; \
  dtmpv4_3 = simd_vshff(dpn_trans_var_v4_3 , dpn_trans_var_v4_2 , 0xee) ; \
  dtmpv4_4 = simd_vshff(dpn_trans_var_v4_5 , dpn_trans_var_v4_4 , 0x44) ; \
  dtmpv4_5 = simd_vshff(dpn_trans_var_v4_5 , dpn_trans_var_v4_4 , 0xee) ; \
  dtmpv4_6 = simd_vshff(dpn_trans_var_v4_7 , dpn_trans_var_v4_6 , 0x44) ; \
  dtmpv4_7 = simd_vshff(dpn_trans_var_v4_7 , dpn_trans_var_v4_6 , 0xee) ; \
  \
  doublev4 dpn_var_v4_0 , dpn_var_v4_1 , dpn_var_v4_2 , dpn_var_v4_3 ; \
  doublev4 dpn_var_v4_4 , dpn_var_v4_5 , dpn_var_v4_6 , dpn_var_v4_7 ; \
  \
  dpn_var_v4_0 = simd_vshff(dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  dpn_var_v4_1 = simd_vshff(dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  dpn_var_v4_2 = simd_vshff(dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  dpn_var_v4_3 = simd_vshff(dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  dpn_var_v4_4 = simd_vshff(dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  dpn_var_v4_5 = simd_vshff(dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  dpn_var_v4_6 = simd_vshff(dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  dpn_var_v4_7 = simd_vshff(dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_local_v4_ptr[0] = dpn_var_v4_0 ; \
  dpn_local_v4_ptr[1] = dpn_var_v4_4 ; \
  dpn_local_v4_ptr[4] = dpn_var_v4_1 ; \
  dpn_local_v4_ptr[5] = dpn_var_v4_5 ; \
  dpn_local_v4_ptr[8] = dpn_var_v4_2 ; \
  dpn_local_v4_ptr[9] = dpn_var_v4_6 ; \
  dpn_local_v4_ptr[12] = dpn_var_v4_3 ; \
  dpn_local_v4_ptr[13] = dpn_var_v4_7 ; \
  \
  dpn_trans_var_v4_0 = dpn_trans_local_v4_ptr[8] ;  \
  dpn_trans_var_v4_1 = dpn_trans_local_v4_ptr[9] ;  \
  dpn_trans_var_v4_2 = dpn_trans_local_v4_ptr[10] ;  \
  dpn_trans_var_v4_3 = dpn_trans_local_v4_ptr[11] ;  \
  dpn_trans_var_v4_4 = dpn_trans_local_v4_ptr[12] ;  \
  dpn_trans_var_v4_5 = dpn_trans_local_v4_ptr[13] ;  \
  dpn_trans_var_v4_6 = dpn_trans_local_v4_ptr[14] ; \
  dpn_trans_var_v4_7 = dpn_trans_local_v4_ptr[15] ; \
  \
  dtmpv4_0 = simd_vshff(dpn_trans_var_v4_1 , dpn_trans_var_v4_0 , 0x44) ; \
  dtmpv4_1 = simd_vshff(dpn_trans_var_v4_1 , dpn_trans_var_v4_0 , 0xee) ; \
  dtmpv4_2 = simd_vshff(dpn_trans_var_v4_3 , dpn_trans_var_v4_2 , 0x44) ; \
  dtmpv4_3 = simd_vshff(dpn_trans_var_v4_3 , dpn_trans_var_v4_2 , 0xee) ; \
  dtmpv4_4 = simd_vshff(dpn_trans_var_v4_5 , dpn_trans_var_v4_4 , 0x44) ; \
  dtmpv4_5 = simd_vshff(dpn_trans_var_v4_5 , dpn_trans_var_v4_4 , 0xee) ; \
  dtmpv4_6 = simd_vshff(dpn_trans_var_v4_7 , dpn_trans_var_v4_6 , 0x44) ; \
  dtmpv4_7 = simd_vshff(dpn_trans_var_v4_7 , dpn_trans_var_v4_6 , 0xee) ; \
  \
  dpn_var_v4_0 = simd_vshff(dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  dpn_var_v4_1 = simd_vshff(dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  dpn_var_v4_2 = simd_vshff(dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  dpn_var_v4_3 = simd_vshff(dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  dpn_var_v4_4 = simd_vshff(dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  dpn_var_v4_5 = simd_vshff(dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  dpn_var_v4_6 = simd_vshff(dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  dpn_var_v4_7 = simd_vshff(dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_local_v4_ptr[2] = dpn_var_v4_0 ; \
  dpn_local_v4_ptr[3] = dpn_var_v4_4 ; \
  dpn_local_v4_ptr[6] = dpn_var_v4_1 ; \
  dpn_local_v4_ptr[7] = dpn_var_v4_5 ; \
  dpn_local_v4_ptr[10] = dpn_var_v4_2 ; \
  dpn_local_v4_ptr[11] = dpn_var_v4_6 ; \
  dpn_local_v4_ptr[14] = dpn_var_v4_3 ; \
  dpn_local_v4_ptr[15] = dpn_var_v4_7 ; \
  \
}  
//vertical_remap fisrt part
void TRANSPOSE_DP3D(remap_q_ppm_param_t *remap_q_ppm_ptr) {
  intv8 send_id_v8 ;
  intv8 send_offset_v8 ;
  intv8 recv_offset_v8 ;
  doublev4 *dp3d_local_v4_ptr ;
  doublev4 *dp3d_trans_local_v4_ptr ;
  doublev4 *dpo_local_v4_ptr ;
  send_id_v8 = remap_q_ppm_ptr->send_id_v8 ;
  send_offset_v8 = remap_q_ppm_ptr->send_offset_v8_dp3d ;
  recv_offset_v8 = remap_q_ppm_ptr->recv_offset_v8_dp3d ;
  dp3d_local_v4_ptr = remap_q_ppm_ptr->dp3d_local_v4_ptr ;
  dp3d_trans_local_v4_ptr = remap_q_ppm_ptr->dp3d_trans_local_v4_ptr ;
  dpo_local_v4_ptr = remap_q_ppm_ptr->dpo_local_v4_ptr ;
  // transpose dp3d dimension(k,j,i) -> dimension(j,i,k)
  TRANSPOSE_ARRAY_DP3D(dp3d_trans_local_v4_ptr , dp3d_local_v4_ptr) ;
  // dp3d array's register commucation
  REGISTER_FULL_COMMUICATION_DP3D(dpo_local_v4_ptr , dp3d_trans_local_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8);
  
}
// compute dpo , z , kid ,array 
void VERTICAL_REMAP_INIT(remap_q_ppm_param_t *remap_q_ppm_ptr  ,int ie  ) {
  doublev4 hyai_ps0_var_v4 ; 
  doublev4 ps_v_var_v4 ;
  doublev4 *hyai_diff_local_v4_ptr ;
  doublev4 *hybi_diff_local_v4_ptr ;
  doublev4 *dpo_local_v4_ptr ;
  doublev4 *dpn_local_v4_ptr ;
  doublev4 *pio_local_v4_ptr ;
  doublev4 *pin_local_v4_ptr ;
  doublev4 *ppmdx_local_v4_ptr ;
  double *z_local_ptr ;
  int *kid_local_ptr ;
  hyai_ps0_var_v4 = remap_q_ppm_ptr->hyai_ps0_var_v4 ;
  hyai_diff_local_v4_ptr = remap_q_ppm_ptr->hyai_diff_local_v4_ptr ; 
  hybi_diff_local_v4_ptr = remap_q_ppm_ptr->hybi_diff_local_v4_ptr ;
  dpo_local_v4_ptr = remap_q_ppm_ptr->dpo_local_v4_ptr ;
  dpn_local_v4_ptr = remap_q_ppm_ptr->dpn_local_v4_ptr ;
  pio_local_v4_ptr = remap_q_ppm_ptr->pio_local_v4_ptr ; 
  pin_local_v4_ptr = remap_q_ppm_ptr->pin_local_v4_ptr ;
  ppmdx_local_v4_ptr =  remap_q_ppm_ptr->ppmdx_local_v4_ptr ;
  z_local_ptr =  remap_q_ppm_ptr->z_local_ptr ;
  kid_local_ptr = remap_q_ppm_ptr->kid_local_ptr ;
  // compute ps_v array
  COMPUTE_PS_V(ps_v_var_v4 , dpo_local_v4_ptr , hyai_ps0_var_v4) ;
  remap_q_ppm_ptr->ps_v_var_v4 = ps_v_var_v4 ;
  // compute dpn array
  COMPUTE_DPN(dpn_local_v4_ptr , hyai_diff_local_v4_ptr , hybi_diff_local_v4_ptr , ps_v_var_v4) ;
  // compute pio and pin array
  COMPUTE_PIO_PIN(pio_local_v4_ptr , pin_local_v4_ptr , dpo_local_v4_ptr , dpn_local_v4_ptr , ie) ;
  // dp3d special deal
  // Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence
  MEMCOPY_DP3D_2X36(dpo_local_v4_ptr) ;
  // compute kid array and z array
  double *pio_local_ptr , *pin_local_ptr , *dpo_local_ptr ;
  pio_local_ptr = (double *) pio_local_v4_ptr ;
  pin_local_ptr = (double *) pin_local_v4_ptr ;
  dpo_local_ptr = (double *) dpo_local_v4_ptr ;
  COMPUTE_KID_Z(kid_local_ptr , z_local_ptr , pio_local_ptr , pin_local_ptr , dpo_local_ptr , ie  ) ;
  // compute_ppm_grids
  COMPUTE_PPM_GRIDS(ppmdx_local_v4_ptr , dpo_local_v4_ptr) ;

}
void TRANSPOSE_DPN(remap_q_ppm_param_t *remap_q_ppm_ptr) {
  intv8 send_id_v8 ;
  intv8 send_offset_v8 ;
  intv8 recv_offset_v8 ;
  doublev4 *dpn_local_v4_ptr ;
  doublev4 *dpn_trans_local_v4_ptr ;
  send_id_v8 = remap_q_ppm_ptr->send_id_v8 ;
  send_offset_v8 = remap_q_ppm_ptr->send_offset_v8_dpn ;
  recv_offset_v8 = remap_q_ppm_ptr->recv_offset_v8_dpn ;
  dpn_local_v4_ptr = remap_q_ppm_ptr->dpn_local_v4_ptr ;
  dpn_trans_local_v4_ptr = remap_q_ppm_ptr->dpn_trans_local_v4_ptr ;
 
  REGISTER_FULL_COMMUICATION_DPN(dpn_trans_local_v4_ptr , dpn_local_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8);
  TRANSPOSE_ARRAY_DPN(dpn_local_v4_ptr , dpn_trans_local_v4_ptr) ;
}
#define TRANSPOSE_ARRAY_FORWARD_2X4X4(array_trans_local_v4_ptr , array_local_v4_ptr , dp3d_local_v4_ptr) { \
  doublev4 array_var_v4_0 , array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
  array_var_v4_0 = array_local_v4_ptr[0] ; \
  array_var_v4_1 = array_local_v4_ptr[8] ; \
  array_var_v4_2 = array_local_v4_ptr[4] ; \
  array_var_v4_3 = array_local_v4_ptr[12] ; \
  array_var_v4_4 = array_local_v4_ptr[32] ; \
  array_var_v4_5 = array_local_v4_ptr[36] ; \
  array_var_v4_6 = array_local_v4_ptr[48] ; \
  array_var_v4_7 = array_local_v4_ptr[52] ; \
  \
  doublev4 dp3d_var_v4_0 , dp3d_var_v4_1 ; \
  dp3d_var_v4_0 = dp3d_local_v4_ptr[0] ; \
  dp3d_var_v4_1 = dp3d_local_v4_ptr[4] ; \
  \
  array_var_v4_0 = array_var_v4_0 * dp3d_var_v4_0 ; \
  array_var_v4_1 = array_var_v4_1 * dp3d_var_v4_1 ; \
  array_var_v4_2 = array_var_v4_2 * dp3d_var_v4_0 ; \
  array_var_v4_3 = array_var_v4_3 * dp3d_var_v4_1 ; \
  array_var_v4_4 = array_var_v4_4 * dp3d_var_v4_0 ; \
  array_var_v4_5 = array_var_v4_5 * dp3d_var_v4_1 ; \
  \
  doublev4 dtmpv4_0 , dtmpv4_1 , dtmpv4_2 , dtmpv4_3 ; \
  doublev4 dtmpv4_4 , dtmpv4_5 , dtmpv4_6 , dtmpv4_7 ; \
  VSHFF(dtmpv4_0 , array_var_v4_2 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_2 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_1 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_1 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_6 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_6 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_5 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_5 , 0xee) ; \
  \
  doublev4 array_trans_var_v4_0 , array_trans_var_v4_1 , array_trans_var_v4_2 , array_trans_var_v4_3 ; \
  doublev4 array_trans_var_v4_4 , array_trans_var_v4_5 , array_trans_var_v4_6 , array_trans_var_v4_7 ; \
  VSHFF(array_trans_var_v4_0 , dtmpv4_4 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_4 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_5 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_5 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_2 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_2 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_3 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_3 , 0xdd) ; \
  \
  array_trans_local_v4_ptr[0] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[1] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[4] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[5] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[8] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[9] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[12] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[13] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[16] ; \
  array_var_v4_1 = array_local_v4_ptr[24] ; \
  array_var_v4_2 = array_local_v4_ptr[20] ; \
  array_var_v4_3 = array_local_v4_ptr[28] ; \
  array_var_v4_4 = array_local_v4_ptr[40] ; \
  array_var_v4_5 = array_local_v4_ptr[44] ; \
  array_var_v4_6 = array_local_v4_ptr[56] ; \
  array_var_v4_7 = array_local_v4_ptr[60] ; \
  \
  dp3d_var_v4_0 = dp3d_local_v4_ptr[8] ; \
  dp3d_var_v4_1 = dp3d_local_v4_ptr[12] ; \
  \
  array_var_v4_0 = array_var_v4_0 * dp3d_var_v4_0 ; \
  array_var_v4_1 = array_var_v4_1 * dp3d_var_v4_1 ; \
  array_var_v4_2 = array_var_v4_2 * dp3d_var_v4_0 ; \
  array_var_v4_3 = array_var_v4_3 * dp3d_var_v4_1 ; \
  array_var_v4_4 = array_var_v4_4 * dp3d_var_v4_0 ; \
  array_var_v4_5 = array_var_v4_5 * dp3d_var_v4_1 ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_2 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_2 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_1 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_1 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_6 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_6 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_5 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_5 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_4 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_4 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_5 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_5 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_2 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_2 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_3 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_3 , 0xdd) ; \
  \
  array_trans_local_v4_ptr[2] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[3] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[6] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[7] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[10] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[11] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[14] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[15] = array_trans_var_v4_7 ; \
}  

#define TRANSPOSE_ARRAY_FORWARD_QDP_2X4X4(array_trans_local_v4_ptr , array_local_v4_ptr) { \
  doublev4 array_var_v4_0 , array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
  array_var_v4_0 = array_local_v4_ptr[0] ; \
  array_var_v4_1 = array_local_v4_ptr[4] ; \
  array_var_v4_2 = array_local_v4_ptr[16] ; \
  array_var_v4_3 = array_local_v4_ptr[20] ; \
  array_var_v4_4 = array_local_v4_ptr[32] ; \
  array_var_v4_5 = array_local_v4_ptr[36] ; \
  array_var_v4_6 = array_local_v4_ptr[48] ; \
  array_var_v4_7 = array_local_v4_ptr[52] ; \
  \
  doublev4 dtmpv4_0 , dtmpv4_1 , dtmpv4_2 , dtmpv4_3 ; \
  doublev4 dtmpv4_4 , dtmpv4_5 , dtmpv4_6 , dtmpv4_7 ; \
  VSHFF(dtmpv4_0 , array_var_v4_2 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_2 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_1 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_1 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_6 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_6 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_5 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_5 , 0xee) ; \
  \
  doublev4 array_trans_var_v4_0 , array_trans_var_v4_1 , array_trans_var_v4_2 , array_trans_var_v4_3 ; \
  doublev4 array_trans_var_v4_4 , array_trans_var_v4_5 , array_trans_var_v4_6 , array_trans_var_v4_7 ; \
  VSHFF(array_trans_var_v4_0 , dtmpv4_4 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_4 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_5 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_5 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_2 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_2 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_3 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_3 , 0xdd) ; \
  \
  array_trans_local_v4_ptr[0] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[1] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[4] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[5] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[8] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[9] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[12] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[13] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[8] ; \
  array_var_v4_1 = array_local_v4_ptr[12] ; \
  array_var_v4_2 = array_local_v4_ptr[24] ; \
  array_var_v4_3 = array_local_v4_ptr[28] ; \
  array_var_v4_4 = array_local_v4_ptr[40] ; \
  array_var_v4_5 = array_local_v4_ptr[44] ; \
  array_var_v4_6 = array_local_v4_ptr[56] ; \
  array_var_v4_7 = array_local_v4_ptr[60] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_2 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_2 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_1 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_1 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_6 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_6 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_5 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_5 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_4 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_4 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_5 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_5 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_2 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_2 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_3 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_3 , 0xdd) ; \
   \
  array_trans_local_v4_ptr[2] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[3] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[6] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[7] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[10] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[11] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[14] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[15] = array_trans_var_v4_7 ; \
  \
}    
#define TRANSPOSE_ARRAY_FORWARD_4X4X4X4(array_trans_local_v4_ptr , array_local_v4_ptr , dp3d_local_v4_ptr) {\
  doublev4 *array_trans_local_v4_ptr_0 , *array_trans_local_v4_ptr_1 , *array_trans_local_v4_ptr_2 , *array_trans_local_v4_ptr_3 ; \
  doublev4 *array_local_v4_ptr_0 , *array_local_v4_ptr_1 , *array_local_v4_ptr_2 , *array_local_v4_ptr_3 ; \
  doublev4 *dp3d_local_v4_ptr_0 , *dp3d_local_v4_ptr_1 , *dp3d_local_v4_ptr_2 ,  *dp3d_local_v4_ptr_3 ; \
  array_trans_local_v4_ptr_0 = &array_trans_local_v4_ptr[0] ; \
  array_trans_local_v4_ptr_1 = &array_trans_local_v4_ptr[16] ; \
  array_trans_local_v4_ptr_2 = &array_trans_local_v4_ptr[32] ; \
  array_trans_local_v4_ptr_3 = &array_trans_local_v4_ptr[48] ; \
  \
  array_local_v4_ptr_0 = &array_local_v4_ptr[0] ; \
  array_local_v4_ptr_1 = &array_local_v4_ptr[1] ; \
  array_local_v4_ptr_2 = &array_local_v4_ptr[2] ; \
  array_local_v4_ptr_3 = &array_local_v4_ptr[3] ; \
  \
  dp3d_local_v4_ptr_0 = &dp3d_local_v4_ptr[0] ; \
  dp3d_local_v4_ptr_1 = &dp3d_local_v4_ptr[1] ; \
  dp3d_local_v4_ptr_2 = &dp3d_local_v4_ptr[2] ; \
  dp3d_local_v4_ptr_3 = &dp3d_local_v4_ptr[3] ; \
  \
  TRANSPOSE_ARRAY_FORWARD_2X4X4(array_trans_local_v4_ptr_0 , array_local_v4_ptr_0 , dp3d_local_v4_ptr_0) ; \
  TRANSPOSE_ARRAY_FORWARD_2X4X4(array_trans_local_v4_ptr_1 , array_local_v4_ptr_1 , dp3d_local_v4_ptr_1) ; \
  TRANSPOSE_ARRAY_FORWARD_2X4X4(array_trans_local_v4_ptr_2 , array_local_v4_ptr_2 , dp3d_local_v4_ptr_2) ; \
  TRANSPOSE_ARRAY_FORWARD_2X4X4(array_trans_local_v4_ptr_3 , array_local_v4_ptr_3 , dp3d_local_v4_ptr_3) ; \
  \
}
//
#define TRANSPOSE_ARRAY_FORWARD_QDP_4X4X4X4(array_trans_local_v4_ptr , array_local_v4_ptr) {\
  doublev4 *array_trans_local_v4_ptr_0 , *array_trans_local_v4_ptr_1 , *array_trans_local_v4_ptr_2 , *array_trans_local_v4_ptr_3 ; \
  doublev4 *array_local_v4_ptr_0 , *array_local_v4_ptr_1 , *array_local_v4_ptr_2 , *array_local_v4_ptr_3 ; \
  array_trans_local_v4_ptr_0 = &array_trans_local_v4_ptr[0] ; \
  array_trans_local_v4_ptr_1 = &array_trans_local_v4_ptr[16] ; \
  array_trans_local_v4_ptr_2 = &array_trans_local_v4_ptr[32] ; \
  array_trans_local_v4_ptr_3 = &array_trans_local_v4_ptr[48] ; \
  \
  array_local_v4_ptr_0 = &array_local_v4_ptr[0] ; \
  array_local_v4_ptr_1 = &array_local_v4_ptr[1] ; \
  array_local_v4_ptr_2 = &array_local_v4_ptr[2] ; \
  array_local_v4_ptr_3 = &array_local_v4_ptr[3] ; \
  \
  TRANSPOSE_ARRAY_FORWARD_QDP_2X4X4(array_trans_local_v4_ptr_0 , array_local_v4_ptr_0) ; \
  TRANSPOSE_ARRAY_FORWARD_QDP_2X4X4(array_trans_local_v4_ptr_1 , array_local_v4_ptr_1) ; \
  TRANSPOSE_ARRAY_FORWARD_QDP_2X4X4(array_trans_local_v4_ptr_2 , array_local_v4_ptr_2) ; \
  TRANSPOSE_ARRAY_FORWARD_QDP_2X4X4(array_trans_local_v4_ptr_3 , array_local_v4_ptr_3) ; \
}

#define COPY_MEMERY_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset) { \
  doublev4 array_var_v4_0, array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
  \
  array_var_v4_0 = array_send_v4_ptr[send_offset] ; \
  array_var_v4_1 = array_send_v4_ptr[send_offset + 1] ; \
  array_var_v4_2 = array_send_v4_ptr[send_offset + 2] ; \
  array_var_v4_3 = array_send_v4_ptr[send_offset + 3] ; \
  array_var_v4_4 = array_send_v4_ptr[send_offset + 4] ; \
  array_var_v4_5 = array_send_v4_ptr[send_offset + 5] ; \
  array_var_v4_6 = array_send_v4_ptr[send_offset + 6] ; \
  array_var_v4_7 = array_send_v4_ptr[send_offset + 7] ; \
  \
  array_recv_v4_ptr[recv_offset + 2] = array_var_v4_0 ; \
  array_recv_v4_ptr[recv_offset + 3] = array_var_v4_1 ; \
  array_recv_v4_ptr[recv_offset + 4] = array_var_v4_2 ; \
  array_recv_v4_ptr[recv_offset + 5] = array_var_v4_3 ; \
  array_recv_v4_ptr[recv_offset + 36] = array_var_v4_4 ; \
  array_recv_v4_ptr[recv_offset + 37] = array_var_v4_5 ; \
  array_recv_v4_ptr[recv_offset + 38] = array_var_v4_6 ; \
  array_recv_v4_ptr[recv_offset + 39] = array_var_v4_7 ; \
  \
}

#define REGISTER_RING_COMMUICATION_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) {\
  doublev4 array_var_v4_0, array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
  \
  array_var_v4_0 = array_send_v4_ptr[send_offset] ; \
  array_var_v4_1 = array_send_v4_ptr[send_offset + 1] ; \
  array_var_v4_2 = array_send_v4_ptr[send_offset + 2] ; \
  array_var_v4_3 = array_send_v4_ptr[send_offset + 3] ; \
  array_var_v4_4 = array_send_v4_ptr[send_offset + 4] ; \
  array_var_v4_5 = array_send_v4_ptr[send_offset + 5] ; \
  array_var_v4_6 = array_send_v4_ptr[send_offset + 6] ; \
  array_var_v4_7 = array_send_v4_ptr[send_offset + 7] ; \
  \
  REGISTER_PUTC(array_var_v4_0 , send_id) ; \
  REGISTER_PUTC(array_var_v4_1 , send_id) ; \
  REGISTER_PUTC(array_var_v4_2 , send_id) ; \
  REGISTER_PUTC(array_var_v4_3 , send_id) ; \
  REGISTER_PUTC(array_var_v4_4 , send_id) ; \
  REGISTER_PUTC(array_var_v4_5 , send_id) ; \
  REGISTER_PUTC(array_var_v4_6 , send_id) ; \
  REGISTER_PUTC(array_var_v4_7 , send_id) ; \
  \
  REGISTER_GETC(array_var_v4_0) ; \
  REGISTER_GETC(array_var_v4_1) ; \
  REGISTER_GETC(array_var_v4_2) ; \
  REGISTER_GETC(array_var_v4_3) ; \
  REGISTER_GETC(array_var_v4_4) ; \
  REGISTER_GETC(array_var_v4_5) ; \
  REGISTER_GETC(array_var_v4_6) ; \
  REGISTER_GETC(array_var_v4_7) ; \
  REGISTER_SYNC(0xff) ; \
  \
  array_recv_v4_ptr[recv_offset + 2] = array_var_v4_0 ; \
  array_recv_v4_ptr[recv_offset + 3] = array_var_v4_1 ; \
  array_recv_v4_ptr[recv_offset + 4] = array_var_v4_2 ; \
  array_recv_v4_ptr[recv_offset + 5] = array_var_v4_3 ; \
  array_recv_v4_ptr[recv_offset + 36] = array_var_v4_4 ; \
  array_recv_v4_ptr[recv_offset + 37] = array_var_v4_5 ; \
  array_recv_v4_ptr[recv_offset + 38] = array_var_v4_6 ; \
  array_recv_v4_ptr[recv_offset + 39] = array_var_v4_7 ; \
}
 
#define REGISTER_FULL_COMMUICATION_FORWARD(array_recv_v4_ptr , array_send_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8) { \
  int send_id ; \
  int send_offset , recv_offset ; \
  \
  send_offset = simd_vextw0(send_offset_v8) ; \
  recv_offset = simd_vextw0(recv_offset_v8) ; \
  COPY_MEMERY_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset) ; \
  \
  send_id = simd_vextw1(send_id_v8) ; \
  send_offset = simd_vextw1(send_offset_v8) ; \
  recv_offset = simd_vextw1(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw2(send_id_v8) ; \
  send_offset = simd_vextw2(send_offset_v8) ; \
  recv_offset = simd_vextw2(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw3(send_id_v8) ; \
  send_offset = simd_vextw3(send_offset_v8) ; \
  recv_offset = simd_vextw3(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw4(send_id_v8) ; \
  send_offset = simd_vextw4(send_offset_v8) ; \
  recv_offset = simd_vextw4(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw5(send_id_v8) ; \
  send_offset = simd_vextw5(send_offset_v8) ; \
  recv_offset = simd_vextw5(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw6(send_id_v8) ; \
  send_offset = simd_vextw6(send_offset_v8) ; \
  recv_offset = simd_vextw6(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw7(send_id_v8) ; \
  send_offset = simd_vextw7(send_offset_v8) ; \
  recv_offset = simd_vextw7(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_FORWARD(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
}

#define COMPUTE_AO_MASSO_2X4X4(masso_local_v4_ptr , array_local_v4_ptr , dpo_local_v4_ptr , masso_var_v4_0) { \
   doublev4 array_var_v4_0 , array_var_v4_1 , array_var_v4_2 ,array_var_v4_3 ; \
   doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
   array_var_v4_0 = array_local_v4_ptr[2] ; \
   array_var_v4_1 = array_local_v4_ptr[3] ; \
   array_var_v4_2 = array_local_v4_ptr[4] ; \
   array_var_v4_3 = array_local_v4_ptr[5] ; \
   array_var_v4_4 = array_local_v4_ptr[6] ; \
   array_var_v4_5 = array_local_v4_ptr[7] ; \
   array_var_v4_6 = array_local_v4_ptr[8] ; \
   array_var_v4_7 = array_local_v4_ptr[9] ; \
   \
   doublev4 masso_var_v4_1 , masso_var_v4_2 , masso_var_v4_3 , masso_var_v4_4 ; \
   doublev4 masso_var_v4_5 , masso_var_v4_6 , masso_var_v4_7, masso_var_v4_8 ; \
   \
   masso_var_v4_1 = masso_var_v4_0 + array_var_v4_0 ; \
   masso_var_v4_2 = masso_var_v4_1 + array_var_v4_1 ; \
   masso_var_v4_3 = masso_var_v4_2 + array_var_v4_2 ; \
   masso_var_v4_4 = masso_var_v4_3 + array_var_v4_3 ; \
   masso_var_v4_5 = masso_var_v4_4 + array_var_v4_4 ; \
   masso_var_v4_6 = masso_var_v4_5 + array_var_v4_5 ; \
   masso_var_v4_7 = masso_var_v4_6 + array_var_v4_6 ; \
   masso_var_v4_8 = masso_var_v4_7 + array_var_v4_7 ; \
   \
   doublev4 dpo_var_v4_0 , dpo_var_v4_1 ; \
   dpo_var_v4_0 = dpo_local_v4_ptr[0] ; \
   dpo_var_v4_1 = dpo_local_v4_ptr[1] ; \
   \
   doublev4 dtmpv4_0 , dtmpv4_1 , dtmpv4_2 , dtmpv4_3 ; \
   doublev4 dtmpv4_4 , dtmpv4_5 , dtmpv4_6 , dtmpv4_7 ; \
   VSHFF(dtmpv4_0 , dpo_var_v4_0 , dpo_var_v4_0 , 0x0) ; \
   VSHFF(dtmpv4_1 , dpo_var_v4_0 , dpo_var_v4_0 , 0x55) ; \
   VSHFF(dtmpv4_2 , dpo_var_v4_0 , dpo_var_v4_0 , 0xaa) ; \
   VSHFF(dtmpv4_3 , dpo_var_v4_0 , dpo_var_v4_0 , 0xff) ; \
   VSHFF(dtmpv4_4 , dpo_var_v4_1 , dpo_var_v4_1 , 0x0) ; \
   VSHFF(dtmpv4_5 , dpo_var_v4_1 , dpo_var_v4_1 , 0x55) ; \
   VSHFF(dtmpv4_6 , dpo_var_v4_1 , dpo_var_v4_1 , 0xaa) ; \
   VSHFF(dtmpv4_7 , dpo_var_v4_1 , dpo_var_v4_1 , 0xff) ; \
   \
   masso_local_v4_ptr[0] = masso_var_v4_0 ; \
   masso_local_v4_ptr[1] = masso_var_v4_1 ; \
   masso_local_v4_ptr[2] = masso_var_v4_2 ; \
   masso_local_v4_ptr[3] = masso_var_v4_3 ; \
   masso_local_v4_ptr[4] = masso_var_v4_4 ; \
   masso_local_v4_ptr[5] = masso_var_v4_5 ; \
   masso_local_v4_ptr[6] = masso_var_v4_6 ; \
   masso_local_v4_ptr[7] = masso_var_v4_7 ; \
   \
   array_var_v4_0 = array_var_v4_0 / dtmpv4_0 ; \
   array_var_v4_1 = array_var_v4_1 / dtmpv4_1 ; \
   array_var_v4_2 = array_var_v4_2 / dtmpv4_2 ; \
   array_var_v4_3 = array_var_v4_3 / dtmpv4_3 ; \
   array_var_v4_4 = array_var_v4_4 / dtmpv4_4 ; \
   array_var_v4_5 = array_var_v4_5 / dtmpv4_5 ; \
   array_var_v4_6 = array_var_v4_6 / dtmpv4_6 ; \
   array_var_v4_7 = array_var_v4_7 / dtmpv4_7 ; \
   \
   array_local_v4_ptr[2] = array_var_v4_0 ; \
   array_local_v4_ptr[3] = array_var_v4_1 ; \
   array_local_v4_ptr[4] = array_var_v4_2 ; \
   array_local_v4_ptr[5] = array_var_v4_3 ; \
   array_local_v4_ptr[6] = array_var_v4_4 ; \
   array_local_v4_ptr[7] = array_var_v4_5 ; \
   array_local_v4_ptr[8] = array_var_v4_6 ; \
   array_local_v4_ptr[9] = array_var_v4_7 ; \
   \
   masso_var_v4_0 = masso_var_v4_8 ; \
}

#define COMPUTE_AO_MASSO(masso_local_v4_ptr , array_local_v4_ptr , dpo_local_v4_ptr) { \
  doublev4 masso_var_v4_0 ; \
  masso_var_v4_0 = 0.0 ; \
  \
  doublev4 *masso_local_v4_ptr_0 , *masso_local_v4_ptr_1 , *masso_local_v4_ptr_2 , *masso_local_v4_ptr_3 ; \
  doublev4 *array_local_v4_ptr_0 , *array_local_v4_ptr_1 , *array_local_v4_ptr_2 , *array_local_v4_ptr_3 ; \
  doublev4 *dpo_local_v4_ptr_0 , *dpo_local_v4_ptr_1 , *dpo_local_v4_ptr_2 , *dpo_local_v4_ptr_3 ; \
  masso_local_v4_ptr_0 = &masso_local_v4_ptr[0] ; \
  masso_local_v4_ptr_1 = &masso_local_v4_ptr[8] ; \
  masso_local_v4_ptr_2 = &masso_local_v4_ptr[16] ; \
  masso_local_v4_ptr_3 = &masso_local_v4_ptr[24] ; \
  array_local_v4_ptr_0 = &array_local_v4_ptr[0] ; \
  array_local_v4_ptr_1 = &array_local_v4_ptr[8] ; \
  array_local_v4_ptr_2 = &array_local_v4_ptr[16] ; \
  array_local_v4_ptr_3 = &array_local_v4_ptr[24] ; \
  dpo_local_v4_ptr_0 = &dpo_local_v4_ptr[1] ; \
  dpo_local_v4_ptr_1 = &dpo_local_v4_ptr[3] ; \
  dpo_local_v4_ptr_2 = &dpo_local_v4_ptr[5] ; \
  dpo_local_v4_ptr_3 = &dpo_local_v4_ptr[7] ; \
  \
  COMPUTE_AO_MASSO_2X4X4(masso_local_v4_ptr_0 , array_local_v4_ptr_0 , dpo_local_v4_ptr_0 , masso_var_v4_0) ; \
  COMPUTE_AO_MASSO_2X4X4(masso_local_v4_ptr_1 , array_local_v4_ptr_1 , dpo_local_v4_ptr_1 , masso_var_v4_0) ; \
  COMPUTE_AO_MASSO_2X4X4(masso_local_v4_ptr_2 , array_local_v4_ptr_2 , dpo_local_v4_ptr_2 , masso_var_v4_0) ; \
  COMPUTE_AO_MASSO_2X4X4(masso_local_v4_ptr_3 , array_local_v4_ptr_3 , dpo_local_v4_ptr_3 , masso_var_v4_0) ; \
  \
}

#define MEMCOPY_2X2X4(array_local_v4_ptr) {\
  doublev4 array_var_v4_0 , array_var_v4_1 , array_var_v4_32 , array_var_v4_33 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[3] ; \
  array_var_v4_1 = array_local_v4_ptr[2] ; \
  array_var_v4_32 = array_local_v4_ptr[31] ; \
  array_var_v4_33 = array_local_v4_ptr[30] ; \
  \
  array_local_v4_ptr[0] = array_var_v4_0 ; \
  array_local_v4_ptr[1] = array_var_v4_1 ; \
  array_local_v4_ptr[32] = array_var_v4_32 ; \
  array_local_v4_ptr[33] = array_var_v4_33 ; \
}

#define COMPUTE_ARRAY_DIFFERNECE(array_diff_local_v4_ptr , array_local_v4_ptr) { \
  doublev4 array_var_v4_0 , array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 , array_var_v4_8; \
  array_var_v4_0 = array_local_v4_ptr[0] ; \
  array_var_v4_1 = array_local_v4_ptr[1] ; \
  array_var_v4_2 = array_local_v4_ptr[2] ; \
  array_var_v4_3 = array_local_v4_ptr[3] ; \
  array_var_v4_4 = array_local_v4_ptr[4] ; \
  array_var_v4_5 = array_local_v4_ptr[5] ; \
  array_var_v4_6 = array_local_v4_ptr[6] ; \
  array_var_v4_7 = array_local_v4_ptr[7] ; \
  array_var_v4_8 = array_local_v4_ptr[8] ; \
  \
  doublev4 array_diff_var_v4_0 , array_diff_var_v4_1 , array_diff_var_v4_2 , array_diff_var_v4_3 ; \
  doublev4 array_diff_var_v4_4 , array_diff_var_v4_5 , array_diff_var_v4_6 , array_diff_var_v4_7 ; \
  array_diff_var_v4_0 = array_var_v4_1 - array_var_v4_0 ; \
  array_diff_var_v4_1 = array_var_v4_2 - array_var_v4_1 ; \
  array_diff_var_v4_2 = array_var_v4_3 - array_var_v4_2 ; \
  array_diff_var_v4_3 = array_var_v4_4 - array_var_v4_3 ; \
  array_diff_var_v4_4 = array_var_v4_5 - array_var_v4_4 ; \
  array_diff_var_v4_5 = array_var_v4_6 - array_var_v4_5 ; \
  array_diff_var_v4_6 = array_var_v4_7 - array_var_v4_6 ; \
  array_diff_var_v4_7 = array_var_v4_8 - array_var_v4_7 ; \
  \
  array_diff_local_v4_ptr[0] = array_diff_var_v4_0 ; \
  array_diff_local_v4_ptr[1] = array_diff_var_v4_1 ; \
  array_diff_local_v4_ptr[2] = array_diff_var_v4_2 ; \
  array_diff_local_v4_ptr[3] = array_diff_var_v4_3 ; \
  array_diff_local_v4_ptr[4] = array_diff_var_v4_4 ; \
  array_diff_local_v4_ptr[5] = array_diff_var_v4_5 ; \
  array_diff_local_v4_ptr[6] = array_diff_var_v4_6 ; \
  array_diff_local_v4_ptr[7] = array_diff_var_v4_7 ; \
  \
  array_var_v4_0 = array_var_v4_8 ; \
  array_var_v4_1 = array_local_v4_ptr[9] ; \
  array_var_v4_2 = array_local_v4_ptr[10] ; \
  array_var_v4_3 = array_local_v4_ptr[11] ; \
  array_var_v4_4 = array_local_v4_ptr[12] ; \
  array_var_v4_5 = array_local_v4_ptr[13] ; \
  array_var_v4_6 = array_local_v4_ptr[14] ; \
  array_var_v4_7 = array_local_v4_ptr[15] ; \
  array_var_v4_8 = array_local_v4_ptr[16] ; \
  \
  array_diff_var_v4_0 = array_var_v4_1 - array_var_v4_0 ; \
  array_diff_var_v4_1 = array_var_v4_2 - array_var_v4_1 ; \
  array_diff_var_v4_2 = array_var_v4_3 - array_var_v4_2 ; \
  array_diff_var_v4_3 = array_var_v4_4 - array_var_v4_3 ; \
  array_diff_var_v4_4 = array_var_v4_5 - array_var_v4_4 ; \
  array_diff_var_v4_5 = array_var_v4_6 - array_var_v4_5 ; \
  array_diff_var_v4_6 = array_var_v4_7 - array_var_v4_6 ; \
  array_diff_var_v4_7 = array_var_v4_8 - array_var_v4_7 ; \
  \
  array_diff_local_v4_ptr[8] = array_diff_var_v4_0 ; \
  array_diff_local_v4_ptr[9] = array_diff_var_v4_1 ; \
  array_diff_local_v4_ptr[10] = array_diff_var_v4_2 ; \
  array_diff_local_v4_ptr[11] = array_diff_var_v4_3 ; \
  array_diff_local_v4_ptr[12] = array_diff_var_v4_4 ; \
  array_diff_local_v4_ptr[13] = array_diff_var_v4_5 ; \
  array_diff_local_v4_ptr[14] = array_diff_var_v4_6 ; \
  array_diff_local_v4_ptr[15] = array_diff_var_v4_7 ; \
  \
  array_var_v4_0 = array_var_v4_8 ; \
  array_var_v4_1 = array_local_v4_ptr[17] ; \
  array_var_v4_2 = array_local_v4_ptr[18] ; \
  array_var_v4_3 = array_local_v4_ptr[19] ; \
  array_var_v4_4 = array_local_v4_ptr[20] ; \
  array_var_v4_5 = array_local_v4_ptr[21] ; \
  array_var_v4_6 = array_local_v4_ptr[22] ; \
  array_var_v4_7 = array_local_v4_ptr[23] ; \
  array_var_v4_8 = array_local_v4_ptr[24] ; \
  \
  array_diff_var_v4_0 = array_var_v4_1 - array_var_v4_0 ; \
  array_diff_var_v4_1 = array_var_v4_2 - array_var_v4_1 ; \
  array_diff_var_v4_2 = array_var_v4_3 - array_var_v4_2 ; \
  array_diff_var_v4_3 = array_var_v4_4 - array_var_v4_3 ; \
  array_diff_var_v4_4 = array_var_v4_5 - array_var_v4_4 ; \
  array_diff_var_v4_5 = array_var_v4_6 - array_var_v4_5 ; \
  array_diff_var_v4_6 = array_var_v4_7 - array_var_v4_6 ; \
  array_diff_var_v4_7 = array_var_v4_8 - array_var_v4_7 ; \
  \
  array_diff_local_v4_ptr[16] = array_diff_var_v4_0 ; \
  array_diff_local_v4_ptr[17] = array_diff_var_v4_1 ; \
  array_diff_local_v4_ptr[18] = array_diff_var_v4_2 ; \
  array_diff_local_v4_ptr[19] = array_diff_var_v4_3 ; \
  array_diff_local_v4_ptr[20] = array_diff_var_v4_4 ; \
  array_diff_local_v4_ptr[21] = array_diff_var_v4_5 ; \
  array_diff_local_v4_ptr[22] = array_diff_var_v4_6 ; \
  array_diff_local_v4_ptr[23] = array_diff_var_v4_7 ; \
  \
  array_var_v4_0 = array_var_v4_8 ; \
  array_var_v4_1 = array_local_v4_ptr[25] ; \
  array_var_v4_2 = array_local_v4_ptr[26] ; \
  array_var_v4_3 = array_local_v4_ptr[27] ; \
  array_var_v4_4 = array_local_v4_ptr[28] ; \
  array_var_v4_5 = array_local_v4_ptr[29] ; \
  array_var_v4_6 = array_local_v4_ptr[30] ; \
  array_var_v4_7 = array_local_v4_ptr[31] ; \
  array_var_v4_8 = array_local_v4_ptr[32] ; \
  \
  array_diff_var_v4_0 = array_var_v4_1 - array_var_v4_0 ; \
  array_diff_var_v4_1 = array_var_v4_2 - array_var_v4_1 ; \
  array_diff_var_v4_2 = array_var_v4_3 - array_var_v4_2 ; \
  array_diff_var_v4_3 = array_var_v4_4 - array_var_v4_3 ; \
  array_diff_var_v4_4 = array_var_v4_5 - array_var_v4_4 ; \
  array_diff_var_v4_5 = array_var_v4_6 - array_var_v4_5 ; \
  array_diff_var_v4_6 = array_var_v4_7 - array_var_v4_6 ; \
  array_diff_var_v4_7 = array_var_v4_8 - array_var_v4_7 ; \
  \
  array_diff_local_v4_ptr[24] = array_diff_var_v4_0 ; \
  array_diff_local_v4_ptr[25] = array_diff_var_v4_1 ; \
  array_diff_local_v4_ptr[26] = array_diff_var_v4_2 ; \
  array_diff_local_v4_ptr[27] = array_diff_var_v4_3 ; \
  array_diff_local_v4_ptr[28] = array_diff_var_v4_4 ; \
  array_diff_local_v4_ptr[29] = array_diff_var_v4_5 ; \
  array_diff_local_v4_ptr[30] = array_diff_var_v4_6 ; \
  array_diff_local_v4_ptr[31] = array_diff_var_v4_7 ; \
  \
  array_var_v4_0 = array_var_v4_8 ; \
  array_var_v4_1 = array_local_v4_ptr[33] ; \
  \
  array_diff_var_v4_0 = array_var_v4_1 - array_var_v4_0 ; \
  array_diff_var_v4_1 = array_var_v4_2 - array_var_v4_1 ; \
  \
  array_diff_local_v4_ptr[32] = array_diff_var_v4_0 ; \
  array_diff_local_v4_ptr[33] = array_diff_var_v4_1 ; \
}

void COMPUTE_PPM(doublev4 *coefs_local_v4_ptr , doublev4 *ai_local_v4_ptr , doublev4 *array_local_v4_ptr , doublev4 *array_diff_local_v4_ptr , doublev4 *ppmdx_local_v4_ptr , doublev4 *dma_local_v4_ptr , int ie , int q) { 
  int i , k ; 
  COMPUTE_ARRAY_DIFFERNECE(array_diff_local_v4_ptr , array_local_v4_ptr) ; 
  for (k = 0 ; k < 32 ; k = k + 4) { 
    doublev4 ppmdx_var_v4_0 , ppmdx_var_v4_1 , ppmdx_var_v4_2 ; 
    ppmdx_var_v4_0 = ppmdx_local_v4_ptr[(k>>2)] ; 
    ppmdx_var_v4_1 = ppmdx_local_v4_ptr[(k>>2) + 16] ; 
    ppmdx_var_v4_2 = ppmdx_local_v4_ptr[(k>>2) + 32] ; 
   
 
    doublev4 dx_var_v4_0_0 , dx_var_v4_0_1 , dx_var_v4_0_2 , dx_var_v4_0_3 ; 
    doublev4 dx_var_v4_1_0 , dx_var_v4_1_1 , dx_var_v4_1_2 , dx_var_v4_1_3 ; 
    doublev4 dx_var_v4_2_0 , dx_var_v4_2_1 , dx_var_v4_2_2 , dx_var_v4_2_3 ; 
    VSHFF(dx_var_v4_2_0 , ppmdx_var_v4_2 , ppmdx_var_v4_2 , 0x0) ; 
    VSHFF(dx_var_v4_2_1 , ppmdx_var_v4_2 , ppmdx_var_v4_2 , 0x55) ; 
    VSHFF(dx_var_v4_2_2 , ppmdx_var_v4_2 , ppmdx_var_v4_2 , 0xaa) ; 
    VSHFF(dx_var_v4_2_3 , ppmdx_var_v4_2 , ppmdx_var_v4_2 , 0xff) ; 
    
    doublev4 ao_diff_var_v4_0 , ao_diff_var_v4_1 , ao_diff_var_v4_2 , ao_diff_var_v4_3 , ao_diff_var_v4_4 ; 
    ao_diff_var_v4_0 = array_diff_local_v4_ptr[k] ; 
    ao_diff_var_v4_1 = array_diff_local_v4_ptr[k + 1] ; 
    ao_diff_var_v4_2 = array_diff_local_v4_ptr[k + 2] ; 
    ao_diff_var_v4_3 = array_diff_local_v4_ptr[k + 3] ; 
    ao_diff_var_v4_4 = array_diff_local_v4_ptr[k + 4] ; 
    
    VSHFF(dx_var_v4_1_0 , ppmdx_var_v4_1 , ppmdx_var_v4_1 , 0x0) ; 
    VSHFF(dx_var_v4_1_1 , ppmdx_var_v4_1 , ppmdx_var_v4_1 , 0x55) ; 
    VSHFF(dx_var_v4_1_2 , ppmdx_var_v4_1 , ppmdx_var_v4_1 , 0xaa) ; 
    VSHFF(dx_var_v4_1_3 , ppmdx_var_v4_1 , ppmdx_var_v4_1 , 0xff) ; 
    
    doublev4 dx_ao_diff_var_v4_0 , dx_ao_diff_var_v4_1 , dx_ao_diff_var_v4_2 , dx_ao_diff_var_v4_3 ; 
    dx_ao_diff_var_v4_0 = dx_var_v4_2_0 * ao_diff_var_v4_0 ; 
    dx_ao_diff_var_v4_1 = dx_var_v4_2_1 * ao_diff_var_v4_1 ; 
    dx_ao_diff_var_v4_2 = dx_var_v4_2_2 * ao_diff_var_v4_2 ; 
    dx_ao_diff_var_v4_3 = dx_var_v4_2_3 * ao_diff_var_v4_3 ; 
   

 
    VSHFF(dx_var_v4_0_0 , ppmdx_var_v4_0 , ppmdx_var_v4_0 , 0x0) ; 
    VSHFF(dx_var_v4_0_1 , ppmdx_var_v4_0 , ppmdx_var_v4_0 , 0x55) ; 
    VSHFF(dx_var_v4_0_2 , ppmdx_var_v4_0 , ppmdx_var_v4_0 , 0xaa) ; 
    VSHFF(dx_var_v4_0_3 , ppmdx_var_v4_0 , ppmdx_var_v4_0 , 0xff) ; 
    
    dx_ao_diff_var_v4_0 = dx_var_v4_1_0 * ao_diff_var_v4_1 + dx_ao_diff_var_v4_0 ; 
    dx_ao_diff_var_v4_1 = dx_var_v4_1_1 * ao_diff_var_v4_2 + dx_ao_diff_var_v4_1 ; 
    dx_ao_diff_var_v4_2 = dx_var_v4_1_2 * ao_diff_var_v4_3 + dx_ao_diff_var_v4_2 ; 
    dx_ao_diff_var_v4_3 = dx_var_v4_1_3 * ao_diff_var_v4_4 + dx_ao_diff_var_v4_3 ; 
    
    dx_ao_diff_var_v4_0 = dx_var_v4_0_0 * dx_ao_diff_var_v4_0 ; 
    dx_ao_diff_var_v4_1 = dx_var_v4_0_1 * dx_ao_diff_var_v4_1 ; 
    dx_ao_diff_var_v4_2 = dx_var_v4_0_2 * dx_ao_diff_var_v4_2 ; 
    dx_ao_diff_var_v4_3 = dx_var_v4_0_3 * dx_ao_diff_var_v4_3 ; 
    
    doublev4 dx_ao_diff_abs_var_v4_0 , dx_ao_diff_abs_var_v4_1 , dx_ao_diff_abs_var_v4_2 , dx_ao_diff_abs_var_v4_3 ; 
    VABS(dx_ao_diff_abs_var_v4_0 , dx_ao_diff_var_v4_0) ; 
    VABS(dx_ao_diff_abs_var_v4_1 , dx_ao_diff_var_v4_1) ; 
    VABS(dx_ao_diff_abs_var_v4_2 , dx_ao_diff_var_v4_2) ; 
    VABS(dx_ao_diff_abs_var_v4_3 , dx_ao_diff_var_v4_3) ; 
    
    doublev4 ao_diff_abs_var_v4_0 , ao_diff_abs_var_v4_1 , ao_diff_abs_var_v4_2 , ao_diff_abs_var_v4_3 , ao_diff_abs_var_v4_4 ; 
    VABS(ao_diff_abs_var_v4_0 , ao_diff_var_v4_0) ; 
    VABS(ao_diff_abs_var_v4_1 , ao_diff_var_v4_1) ; 
    VABS(ao_diff_abs_var_v4_2 , ao_diff_var_v4_2) ; 
    VABS(ao_diff_abs_var_v4_3 , ao_diff_var_v4_3) ; 
    VABS(ao_diff_abs_var_v4_4 , ao_diff_var_v4_4) ; 
    
    doublev4 ao_diff_mul_var_v4_0 , ao_diff_mul_var_v4_1 , ao_diff_mul_var_v4_2 , ao_diff_mul_var_v4_3 ; 
    ao_diff_mul_var_v4_0 = ao_diff_var_v4_0 * ao_diff_var_v4_1 ; 
    ao_diff_mul_var_v4_1 = ao_diff_var_v4_1 * ao_diff_var_v4_2 ; 
    ao_diff_mul_var_v4_2 = ao_diff_var_v4_2 * ao_diff_var_v4_3 ; 
    ao_diff_mul_var_v4_3 = ao_diff_var_v4_3 * ao_diff_var_v4_4 ; 
    
    doublev4 ao_diff_min_var_v4_0 , ao_diff_min_var_v4_1 , ao_diff_min_var_v4_2 , ao_diff_min_var_v4_3 ; 
    ao_diff_min_var_v4_0 = ao_diff_abs_var_v4_1 - ao_diff_abs_var_v4_0 ; 
    ao_diff_min_var_v4_1 = ao_diff_abs_var_v4_2 - ao_diff_abs_var_v4_1 ; 
    ao_diff_min_var_v4_2 = ao_diff_abs_var_v4_3 - ao_diff_abs_var_v4_2 ; 
    ao_diff_min_var_v4_3 = ao_diff_abs_var_v4_4 - ao_diff_abs_var_v4_3 ; 
    
    VMIN(ao_diff_min_var_v4_0 , ao_diff_min_var_v4_0 , ao_diff_abs_var_v4_1 , ao_diff_abs_var_v4_0) ; 
    VMIN(ao_diff_min_var_v4_1 , ao_diff_min_var_v4_1 , ao_diff_abs_var_v4_2 , ao_diff_abs_var_v4_1) ; 
    VMIN(ao_diff_min_var_v4_2 , ao_diff_min_var_v4_2 , ao_diff_abs_var_v4_3 , ao_diff_abs_var_v4_2) ; 
    VMIN(ao_diff_min_var_v4_3 , ao_diff_min_var_v4_3 , ao_diff_abs_var_v4_4 , ao_diff_abs_var_v4_3) ; 
    
    ao_diff_min_var_v4_0 = 2.0 *ao_diff_min_var_v4_0 ; 
    ao_diff_min_var_v4_1 = 2.0 *ao_diff_min_var_v4_1 ; 
    ao_diff_min_var_v4_2 = 2.0 *ao_diff_min_var_v4_2 ; 
    ao_diff_min_var_v4_3 = 2.0 *ao_diff_min_var_v4_3 ; 
    
    doublev4 dma_var_v4_0 , dma_var_v4_1 , dma_var_v4_2 , dma_var_v4_3 ; 
    dma_var_v4_0 = dx_ao_diff_abs_var_v4_0 - ao_diff_min_var_v4_0 ; 
    dma_var_v4_1 = dx_ao_diff_abs_var_v4_1 - ao_diff_min_var_v4_1 ; 
    dma_var_v4_2 = dx_ao_diff_abs_var_v4_2 - ao_diff_min_var_v4_2 ; 
    dma_var_v4_3 = dx_ao_diff_abs_var_v4_3 - ao_diff_min_var_v4_3 ; 
    
    VSIGN(ao_diff_min_var_v4_0 , dx_ao_diff_var_v4_0 , ao_diff_min_var_v4_0) ; 
    VSIGN(ao_diff_min_var_v4_1 , dx_ao_diff_var_v4_1 , ao_diff_min_var_v4_1) ; 
    VSIGN(ao_diff_min_var_v4_2 , dx_ao_diff_var_v4_2 , ao_diff_min_var_v4_2) ; 
    VSIGN(ao_diff_min_var_v4_3 , dx_ao_diff_var_v4_3 , ao_diff_min_var_v4_3) ; 
    
    VMIN(dma_var_v4_0 , dma_var_v4_0 , dx_ao_diff_var_v4_0 , ao_diff_min_var_v4_0) ; 
    VMIN(dma_var_v4_1 , dma_var_v4_1 , dx_ao_diff_var_v4_1 , ao_diff_min_var_v4_1) ; 
    VMIN(dma_var_v4_2 , dma_var_v4_2 , dx_ao_diff_var_v4_2 , ao_diff_min_var_v4_2) ; 
    VMIN(dma_var_v4_3 , dma_var_v4_3 , dx_ao_diff_var_v4_3 , ao_diff_min_var_v4_3) ; 
   

 
    VMINEZ(dma_var_v4_0 , ao_diff_mul_var_v4_0 , dma_var_v4_0) ; 
    VMINEZ(dma_var_v4_1 , ao_diff_mul_var_v4_1 , dma_var_v4_1) ; 
    VMINEZ(dma_var_v4_2 , ao_diff_mul_var_v4_2 , dma_var_v4_2) ; 
    VMINEZ(dma_var_v4_3 , ao_diff_mul_var_v4_3 , dma_var_v4_3) ; 
   
 
    dma_local_v4_ptr[k] = dma_var_v4_0 ; 
    dma_local_v4_ptr[k + 1] = dma_var_v4_1 ; 
    dma_local_v4_ptr[k + 2] = dma_var_v4_2 ; 
    dma_local_v4_ptr[k + 3] = dma_var_v4_3 ; 
    
  } 
  
  for (k = 0 ; k < 32 ; k = k + 4) { 
    doublev4 ppmdx_var_v4_9 , ppmdx_var_v4_8 ; 
    ppmdx_var_v4_9 = ppmdx_local_v4_ptr[(k>>2) + 144] ; 
    ppmdx_var_v4_8 = ppmdx_local_v4_ptr[(k>>2) + 128] ; 
    

    doublev4 dx_var_v4_9_0 , dx_var_v4_9_1 , dx_var_v4_9_2 , dx_var_v4_9_3 ; 
    doublev4 dx_var_v4_8_0 , dx_var_v4_8_1 , dx_var_v4_8_2 , dx_var_v4_8_3 ; 
    VSHFF(dx_var_v4_9_0 , ppmdx_var_v4_9 , ppmdx_var_v4_9 , 0x0) ; 
    VSHFF(dx_var_v4_9_1 , ppmdx_var_v4_9 , ppmdx_var_v4_9 , 0x55) ; 
    VSHFF(dx_var_v4_9_2 , ppmdx_var_v4_9 , ppmdx_var_v4_9 , 0xaa) ; 
    VSHFF(dx_var_v4_9_3 , ppmdx_var_v4_9 , ppmdx_var_v4_9 , 0xff) ; 
    
    VSHFF(dx_var_v4_8_0 , ppmdx_var_v4_8 , ppmdx_var_v4_8 , 0x0) ; 
    VSHFF(dx_var_v4_8_1 , ppmdx_var_v4_8 , ppmdx_var_v4_8 , 0x55) ; 
    VSHFF(dx_var_v4_8_2 , ppmdx_var_v4_8 , ppmdx_var_v4_8 , 0xaa) ; 
    VSHFF(dx_var_v4_8_3 , ppmdx_var_v4_8 , ppmdx_var_v4_8 , 0xff) ; 
    
    doublev4 dma_var_v4_0 , dma_var_v4_1 ,dma_var_v4_2 , dma_var_v4_3 ,dma_var_v4_4 ; 
    dma_var_v4_0 = dma_local_v4_ptr[k] ; 
    dma_var_v4_1 = dma_local_v4_ptr[k + 1] ; 
    dma_var_v4_2 = dma_local_v4_ptr[k + 2] ; 
    dma_var_v4_3 = dma_local_v4_ptr[k + 3] ; 
    dma_var_v4_4 = dma_local_v4_ptr[k + 4] ; 
  
   
    doublev4 dx_dma_var_v4_0 , dx_dma_var_v4_1 , dx_dma_var_v4_2 , dx_dma_var_v4_3 ; 
    dx_dma_var_v4_0 = dx_var_v4_9_0 * dma_var_v4_0 ; 
    dx_dma_var_v4_1 = dx_var_v4_9_1 * dma_var_v4_1 ; 
    dx_dma_var_v4_2 = dx_var_v4_9_2 * dma_var_v4_2 ; 
    dx_dma_var_v4_3 = dx_var_v4_9_3 * dma_var_v4_3 ; 
    
    dx_dma_var_v4_0 = dx_dma_var_v4_0 - dx_var_v4_8_0 * dma_var_v4_1 ; 
    dx_dma_var_v4_1 = dx_dma_var_v4_1 - dx_var_v4_8_1 * dma_var_v4_2 ; 
    dx_dma_var_v4_2 = dx_dma_var_v4_2 - dx_var_v4_8_2 * dma_var_v4_3 ; 
    dx_dma_var_v4_3 = dx_dma_var_v4_3 - dx_var_v4_8_3 * dma_var_v4_4 ; 
    
    doublev4 ppmdx_var_v4_5 , ppmdx_var_v4_6 , ppmdx_var_v4_7 ; 
    ppmdx_var_v4_7 = ppmdx_local_v4_ptr[(k>>2) + 112] ; 
    ppmdx_var_v4_6 = ppmdx_local_v4_ptr[(k>>2) + 96] ; 
    ppmdx_var_v4_5 = ppmdx_local_v4_ptr[(k>>2) + 80] ; 
   

  
    doublev4 dx_var_v4_7_0 , dx_var_v4_7_1 , dx_var_v4_7_2 , dx_var_v4_7_3 ; 
    doublev4 dx_var_v4_6_0 , dx_var_v4_6_1 , dx_var_v4_6_2 , dx_var_v4_6_3 ; 
    doublev4 dx_var_v4_5_0 , dx_var_v4_5_1 , dx_var_v4_5_2 , dx_var_v4_5_3 ; 
    VSHFF(dx_var_v4_7_0 , ppmdx_var_v4_7 , ppmdx_var_v4_7 , 0x0) ; 
    VSHFF(dx_var_v4_7_1 , ppmdx_var_v4_7 , ppmdx_var_v4_7 , 0x55) ; 
    VSHFF(dx_var_v4_7_2 , ppmdx_var_v4_7 , ppmdx_var_v4_7 , 0xaa) ; 
    VSHFF(dx_var_v4_7_3 , ppmdx_var_v4_7 , ppmdx_var_v4_7 , 0xff) ; 
    
    VSHFF(dx_var_v4_6_0 , ppmdx_var_v4_6 , ppmdx_var_v4_6 , 0x0) ; 
    VSHFF(dx_var_v4_6_1 , ppmdx_var_v4_6 , ppmdx_var_v4_6 , 0x55) ; 
    VSHFF(dx_var_v4_6_2 , ppmdx_var_v4_6 , ppmdx_var_v4_6 , 0xaa) ; 
    VSHFF(dx_var_v4_6_3 , ppmdx_var_v4_6 , ppmdx_var_v4_6 , 0xff) ; 
    
    VSHFF(dx_var_v4_5_0 , ppmdx_var_v4_5 , ppmdx_var_v4_5 , 0x0) ; 
    VSHFF(dx_var_v4_5_1 , ppmdx_var_v4_5 , ppmdx_var_v4_5 , 0x55) ; 
    VSHFF(dx_var_v4_5_2 , ppmdx_var_v4_5 , ppmdx_var_v4_5 , 0xaa) ; 
    VSHFF(dx_var_v4_5_3 , ppmdx_var_v4_5 , ppmdx_var_v4_5 , 0xff) ; 
    
    doublev4 dx_diff_var_v4_0 , dx_diff_var_v4_1 , dx_diff_var_v4_2 , dx_diff_var_v4_3 ; 
    dx_diff_var_v4_0 = dx_var_v4_6_0 - dx_var_v4_7_0 ; 
    dx_diff_var_v4_1 = dx_var_v4_6_1 - dx_var_v4_7_1 ; 
    dx_diff_var_v4_2 = dx_var_v4_6_2 - dx_var_v4_7_2 ; 
    dx_diff_var_v4_3 = dx_var_v4_6_3 - dx_var_v4_7_3 ; 
    
    doublev4 dx_diff_mul_var_v4_0 , dx_diff_mul_var_v4_1 , dx_diff_mul_var_v4_2 , dx_diff_mul_var_v4_3 ; 
    dx_diff_mul_var_v4_0 = dx_var_v4_5_0 * dx_diff_var_v4_0 ; 
    dx_diff_mul_var_v4_1 = dx_var_v4_5_1 * dx_diff_var_v4_1 ; 
    dx_diff_mul_var_v4_2 = dx_var_v4_5_2 * dx_diff_var_v4_2 ; 
    dx_diff_mul_var_v4_3 = dx_var_v4_5_3 * dx_diff_var_v4_3 ; 
    
    doublev4 ao_diff_var_v4_1 , ao_diff_var_v4_2 , ao_diff_var_v4_3 , ao_diff_var_v4_4 ; 
    ao_diff_var_v4_1 = array_diff_local_v4_ptr[k + 1] ; 
    ao_diff_var_v4_2 = array_diff_local_v4_ptr[k + 2] ; 
    ao_diff_var_v4_3 = array_diff_local_v4_ptr[k + 3] ; 
    ao_diff_var_v4_4 = array_diff_local_v4_ptr[k + 4] ; 
    
    dx_dma_var_v4_0 = dx_dma_var_v4_0 + ao_diff_var_v4_1 * dx_diff_mul_var_v4_0 ; 
    dx_dma_var_v4_1 = dx_dma_var_v4_1 + ao_diff_var_v4_2 * dx_diff_mul_var_v4_1 ; 
    dx_dma_var_v4_2 = dx_dma_var_v4_2 + ao_diff_var_v4_3 * dx_diff_mul_var_v4_2 ; 
    dx_dma_var_v4_3 = dx_dma_var_v4_3 + ao_diff_var_v4_4 * dx_diff_mul_var_v4_3 ; 
    
    doublev4 ppmdx_var_v4_3 , ppmdx_var_v4_4 ; 
    ppmdx_var_v4_4 = ppmdx_local_v4_ptr[(k>>2) + 64] ; 
    ppmdx_var_v4_3 = ppmdx_local_v4_ptr[(k>>2) + 48] ; 
   
    
    doublev4 dx_var_v4_4_0 , dx_var_v4_4_1 , dx_var_v4_4_2 , dx_var_v4_4_3 ; 
    doublev4 dx_var_v4_3_0 , dx_var_v4_3_1 , dx_var_v4_3_2 , dx_var_v4_3_3 ; 
    VSHFF(dx_var_v4_4_0 , ppmdx_var_v4_4 , ppmdx_var_v4_4 , 0x0) ; 
    VSHFF(dx_var_v4_4_1 , ppmdx_var_v4_4 , ppmdx_var_v4_4 , 0x55) ; 
    VSHFF(dx_var_v4_4_2 , ppmdx_var_v4_4 , ppmdx_var_v4_4 , 0xaa) ; 
    VSHFF(dx_var_v4_4_3 , ppmdx_var_v4_4 , ppmdx_var_v4_4 , 0xff) ; 
    
    VSHFF(dx_var_v4_3_0 , ppmdx_var_v4_3 , ppmdx_var_v4_3 , 0x0) ; 
    VSHFF(dx_var_v4_3_1 , ppmdx_var_v4_3 , ppmdx_var_v4_3 , 0x55) ; 
    VSHFF(dx_var_v4_3_2 , ppmdx_var_v4_3 , ppmdx_var_v4_3 , 0xaa) ; 
    VSHFF(dx_var_v4_3_3 , ppmdx_var_v4_3 , ppmdx_var_v4_3 , 0xff) ; 
    
    doublev4 ao_var_v4_1 , ao_var_v4_2 , ao_var_v4_3 , ao_var_v4_4 ; 
    ao_var_v4_1 = array_local_v4_ptr[k + 1] ; 
    ao_var_v4_2 = array_local_v4_ptr[k + 2] ; 
    ao_var_v4_3 = array_local_v4_ptr[k + 3] ; 
    ao_var_v4_4 = array_local_v4_ptr[k + 4] ; 
   
    
    doublev4 ai_var_v4_0 , ai_var_v4_1 , ai_var_v4_2 , ai_var_v4_3 ; 
    ai_var_v4_0 = ao_var_v4_1 + dx_var_v4_3_0 * ao_diff_var_v4_1 ;
    ai_var_v4_1 = ao_var_v4_2 + dx_var_v4_3_1 * ao_diff_var_v4_2 ;
    ai_var_v4_2 = ao_var_v4_3 + dx_var_v4_3_2 * ao_diff_var_v4_3 ;
    ai_var_v4_3 = ao_var_v4_4 + dx_var_v4_3_3 * ao_diff_var_v4_4 ;

   
    ai_var_v4_0 = ai_var_v4_0 + dx_dma_var_v4_0 * dx_var_v4_4_0 ; 
    ai_var_v4_1 = ai_var_v4_1 + dx_dma_var_v4_1 * dx_var_v4_4_1 ; 
    ai_var_v4_2 = ai_var_v4_2 + dx_dma_var_v4_2 * dx_var_v4_4_2 ; 
    ai_var_v4_3 = ai_var_v4_3 + dx_dma_var_v4_3 * dx_var_v4_4_3 ; 
   
    
    ai_local_v4_ptr[k] = ai_var_v4_0 ; 
    ai_local_v4_ptr[k + 1] = ai_var_v4_1 ; 
    ai_local_v4_ptr[k + 2] = ai_var_v4_2 ; 
    ai_local_v4_ptr[k + 3] = ai_var_v4_3 ; 
    
  } 
  for (k = 1 ; k < 31 ; k++) { 
    doublev4 al_var_v4 , ar_var_v4 ; 
    al_var_v4 = ai_local_v4_ptr[k - 1] ; 
    ar_var_v4 = ai_local_v4_ptr[k] ; 
    

    doublev4 ao_var_v4 ; 
    ao_var_v4 = array_local_v4_ptr[k + 1] ; 
    
    doublev4 ar_ao_diff_var_v4 ; 
    doublev4 ao_al_diff_var_v4 ; 
    
    ar_ao_diff_var_v4 = ar_var_v4 - ao_var_v4 ; 
    ao_al_diff_var_v4 = ao_var_v4 - al_var_v4 ; 
    
    doublev4 ar_ao_al_mul_var ; 
    ar_ao_al_mul_var = ar_ao_diff_var_v4 * ao_al_diff_var_v4 ; 
    VMINE(al_var_v4 , ar_ao_al_mul_var , ao_var_v4 , al_var_v4) ; 
    VMINE(ar_var_v4 , ar_ao_al_mul_var , ao_var_v4 , ar_var_v4) ; 

   
 
    doublev4 ar_al_diff_var_v4 , ar_al_sum_var_v4 ; 
    ar_al_diff_var_v4 = ar_var_v4 - al_var_v4 ; 
    ar_al_sum_var_v4 = ar_var_v4 + al_var_v4 ; 
    
    doublev4 ao_al_ar_var_v4 ; 
    ao_al_ar_var_v4 = ao_var_v4 - 0.5*ar_al_sum_var_v4 ; 
    
    doublev4 ar_al_pow_var_v4 ; 
    ar_al_pow_var_v4 = ar_al_diff_var_v4 * ar_al_diff_var_v4 * 0.166666666666666666666666666666666666; 
    

    doublev4 judgement_var_v4_0 , judgement_var_v4_1 ; 
    judgement_var_v4_0 = ar_al_pow_var_v4 - ar_al_diff_var_v4 * ao_al_ar_var_v4 ; 
    judgement_var_v4_1 = ar_al_diff_var_v4 * ao_al_ar_var_v4 + ar_al_pow_var_v4 ; 
    
    doublev4 al_tmp_var_v4 , ar_tmp_var_v4 ; 
    al_tmp_var_v4 = (3.0*ao_var_v4 - 2.0*ar_var_v4) ; 
    ar_tmp_var_v4 = (3.0*ao_var_v4 - 2.0*al_var_v4) ; 
   
 
    VMINE(al_var_v4 , judgement_var_v4_0 , al_tmp_var_v4 , al_var_v4) ; 
    VMINE(ar_var_v4 , judgement_var_v4_1 , ar_tmp_var_v4 , ar_var_v4) ; 
    
 
    ar_al_sum_var_v4 = ar_var_v4 + al_var_v4 ; 
    ar_al_diff_var_v4 = ar_var_v4 - al_var_v4 ; 
    
    doublev4 coefs_var_v4_0 , coefs_var_v4_1 , coefs_var_v4_2 ; 
    coefs_var_v4_0 = 1.5*ao_var_v4 - 0.25*ar_al_sum_var_v4 ; 
    coefs_var_v4_1 = ar_al_diff_var_v4 ; 
    coefs_var_v4_2 = -6.0*ao_var_v4 + 3.0*ar_al_sum_var_v4 ; 
    
    coefs_local_v4_ptr[k] = coefs_var_v4_0 ; 
    coefs_local_v4_ptr[k + 32] = coefs_var_v4_1 ; 
    coefs_local_v4_ptr[k + 64] = coefs_var_v4_2 ; 
    
  } 
}
void COMPUTE_MASSN(doublev4 *array_local_v4_ptr , doublev4 *coefs_local_v4_ptr , doublev4 *masso_local_v4_ptr , double *z_local_ptr , int *kid_local_ptr , double *dpo_local_ptr , int ie ,int q) {
  doublev4 massn_var_v4_1 , massn_var_v4_2 ; 
  massn_var_v4_1 = 0.0 ; 
  int k ; 
  for(k = 0 ; k < 30 ; k ++) { 
    int kk ; 
    kk = kid_local_ptr[k] ; 
    
    double dpo_var ; 
    dpo_var = dpo_local_ptr[kk + 4] ; 
    
    double z_var ; 
    z_var = z_local_ptr[k] ; 
    doublev4 dpo_var_v4 , z_var_v4 ; 
    VSHFF(dpo_var_v4 , dpo_var , dpo_var , 0x0) ; 
    VSHFF(z_var_v4 , z_var , z_var , 0x0) ; 
    
    doublev4 coefs_var_v4_0 , coefs_var_v4_1 , coefs_var_v4_2 ; 
    coefs_var_v4_0 = coefs_local_v4_ptr[kk + 1] ; 
    coefs_var_v4_1 = coefs_local_v4_ptr[kk + 33] ; 
    coefs_var_v4_2 = coefs_local_v4_ptr[kk + 65] ; 
    
    doublev4 masso_var_v4 ; 
    masso_var_v4 = masso_local_v4_ptr[kk] ; 
    
    doublev4 z_square_var_v4 ; 
    z_square_var_v4 = z_var_v4 * z_var_v4 ;
    coefs_var_v4_1 = coefs_var_v4_1 * 0.5 ; 
    coefs_var_v4_2 = coefs_var_v4_2 * 0.3333333333333333333333 ; 
    
    doublev4 ztmpv4_0 , ztmpv4_1 , ztmpv4_2 ; 
    ztmpv4_0 = z_var_v4 + 0.5 ; 
    ztmpv4_1 = z_var_v4*z_var_v4  - 0.25 ; 
    ztmpv4_2 = z_var_v4*z_square_var_v4 + 0.125 ; 
    doublev4  mass_var_v4 ;
    mass_var_v4 = coefs_var_v4_0 * ztmpv4_0 ; 
    mass_var_v4 = mass_var_v4 + coefs_var_v4_1 * ztmpv4_1 ;
    mass_var_v4 = mass_var_v4 + coefs_var_v4_2 * ztmpv4_2 ;
    massn_var_v4_2 = masso_var_v4 + dpo_var_v4 * mass_var_v4 ; 
    array_local_v4_ptr[k + 2] = massn_var_v4_2 - massn_var_v4_1 ; 
    
    massn_var_v4_1 = massn_var_v4_2 ;
  
  } 
}
#define COPY_MEMERY_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset) { \
  doublev4 array_var_v4_0, array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
  \
  array_var_v4_0 = array_send_v4_ptr[send_offset + 2] ; \
  array_var_v4_1 = array_send_v4_ptr[send_offset + 3] ; \
  array_var_v4_2 = array_send_v4_ptr[send_offset + 4] ; \
  array_var_v4_3 = array_send_v4_ptr[send_offset + 5] ; \
  array_var_v4_4 = array_send_v4_ptr[send_offset + 36] ; \
  array_var_v4_5 = array_send_v4_ptr[send_offset + 37] ; \
  array_var_v4_6 = array_send_v4_ptr[send_offset + 38] ; \
  array_var_v4_7 = array_send_v4_ptr[send_offset + 39] ; \
  \
  array_recv_v4_ptr[recv_offset + 0] = array_var_v4_0 ; \
  array_recv_v4_ptr[recv_offset + 1] = array_var_v4_1 ; \
  array_recv_v4_ptr[recv_offset + 2] = array_var_v4_2 ; \
  array_recv_v4_ptr[recv_offset + 3] = array_var_v4_3 ; \
  array_recv_v4_ptr[recv_offset + 4] = array_var_v4_4 ; \
  array_recv_v4_ptr[recv_offset + 5] = array_var_v4_5 ; \
  array_recv_v4_ptr[recv_offset + 6] = array_var_v4_6 ; \
  array_recv_v4_ptr[recv_offset + 7] = array_var_v4_7 ; \
  \
}
#define REGISTER_RING_COMMUICATION_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) {\
  doublev4 array_var_v4_0, array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
  \
  array_var_v4_0 = array_send_v4_ptr[send_offset + 2] ; \
  array_var_v4_1 = array_send_v4_ptr[send_offset + 3] ; \
  array_var_v4_2 = array_send_v4_ptr[send_offset + 4] ; \
  array_var_v4_3 = array_send_v4_ptr[send_offset + 5] ; \
  array_var_v4_4 = array_send_v4_ptr[send_offset + 36] ; \
  array_var_v4_5 = array_send_v4_ptr[send_offset + 37] ; \
  array_var_v4_6 = array_send_v4_ptr[send_offset + 38] ; \
  array_var_v4_7 = array_send_v4_ptr[send_offset + 39] ; \
  \
  REGISTER_PUTC(array_var_v4_0 , send_id) ; \
  REGISTER_PUTC(array_var_v4_1 , send_id) ; \
  REGISTER_PUTC(array_var_v4_2 , send_id) ; \
  REGISTER_PUTC(array_var_v4_3 , send_id) ; \
  REGISTER_PUTC(array_var_v4_4 , send_id) ; \
  REGISTER_PUTC(array_var_v4_5 , send_id) ; \
  REGISTER_PUTC(array_var_v4_6 , send_id) ; \
  REGISTER_PUTC(array_var_v4_7 , send_id) ; \
  \
  REGISTER_GETC(array_var_v4_0) ; \
  REGISTER_GETC(array_var_v4_1) ; \
  REGISTER_GETC(array_var_v4_2) ; \
  REGISTER_GETC(array_var_v4_3) ; \
  REGISTER_GETC(array_var_v4_4) ; \
  REGISTER_GETC(array_var_v4_5) ; \
  REGISTER_GETC(array_var_v4_6) ; \
  REGISTER_GETC(array_var_v4_7) ; \
  REGISTER_SYNC(0xff) ; \
  \
  array_recv_v4_ptr[recv_offset] = array_var_v4_0 ; \
  array_recv_v4_ptr[recv_offset + 1] = array_var_v4_1 ; \
  array_recv_v4_ptr[recv_offset + 2] = array_var_v4_2 ; \
  array_recv_v4_ptr[recv_offset + 3] = array_var_v4_3 ; \
  array_recv_v4_ptr[recv_offset + 4] = array_var_v4_4 ; \
  array_recv_v4_ptr[recv_offset + 5] = array_var_v4_5 ; \
  array_recv_v4_ptr[recv_offset + 6] = array_var_v4_6 ; \
  array_recv_v4_ptr[recv_offset + 7] = array_var_v4_7 ; \
}
#define REGISTER_FULL_COMMUICATION_POST(array_recv_v4_ptr , array_send_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8) { \
  int send_id ; \
  int send_offset , recv_offset ; \
  \
  send_offset = simd_vextw0(send_offset_v8) ; \
  recv_offset = simd_vextw0(recv_offset_v8) ; \
  COPY_MEMERY_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset) ; \
  \
  send_id = simd_vextw1(send_id_v8) ; \
  send_offset = simd_vextw1(send_offset_v8) ; \
  recv_offset = simd_vextw1(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw2(send_id_v8) ; \
  send_offset = simd_vextw2(send_offset_v8) ; \
  recv_offset = simd_vextw2(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw3(send_id_v8) ; \
  send_offset = simd_vextw3(send_offset_v8) ; \
  recv_offset = simd_vextw3(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw4(send_id_v8) ; \
  send_offset = simd_vextw4(send_offset_v8) ; \
  recv_offset = simd_vextw4(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw5(send_id_v8) ; \
  send_offset = simd_vextw5(send_offset_v8) ; \
  recv_offset = simd_vextw5(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw6(send_id_v8) ; \
  send_offset = simd_vextw6(send_offset_v8) ; \
  recv_offset = simd_vextw6(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
  send_id = simd_vextw7(send_id_v8) ; \
  send_offset = simd_vextw7(send_offset_v8) ; \
  recv_offset = simd_vextw7(recv_offset_v8) ; \
  REGISTER_RING_COMMUICATION_POST(array_recv_v4_ptr , recv_offset , array_send_v4_ptr , send_offset , send_id) ; \
  \
}

#define TRANSPOSE_ARRAY_POST_4X4X4X4(array_trans_local_v4_ptr , array_local_v4_ptr , dpn_local_v4_ptr) { \
  doublev4 array_var_v4_0 , array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
  array_var_v4_0 = array_local_v4_ptr[0] ; \
  array_var_v4_1 = array_local_v4_ptr[4] ; \
  array_var_v4_2 = array_local_v4_ptr[8] ; \
  array_var_v4_3 = array_local_v4_ptr[12] ; \
  array_var_v4_4 = array_local_v4_ptr[16] ; \
  array_var_v4_5 = array_local_v4_ptr[20] ; \
  array_var_v4_6 = array_local_v4_ptr[24] ; \
  array_var_v4_7 = array_local_v4_ptr[28] ; \
  \
  doublev4 dtmpv4_0 , dtmpv4_1 , dtmpv4_2 , dtmpv4_3 ; \
  doublev4 dtmpv4_4 , dtmpv4_5 , dtmpv4_6 , dtmpv4_7 ; \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  doublev4 array_trans_var_v4_0 , array_trans_var_v4_1 , array_trans_var_v4_2 , array_trans_var_v4_3 ; \
  doublev4 array_trans_var_v4_4 , array_trans_var_v4_5 , array_trans_var_v4_6 , array_trans_var_v4_7 ; \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  doublev4 dpn_var_v4_0 , dpn_var_v4_1 ; \
  dpn_var_v4_0 = dpn_local_v4_ptr[0] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[1] ; \
  \
  array_trans_var_v4_0 = array_trans_var_v4_0 / dpn_var_v4_0 ; \
  array_trans_var_v4_1 = array_trans_var_v4_1 / dpn_var_v4_0 ; \
  array_trans_var_v4_2 = array_trans_var_v4_2 / dpn_var_v4_0 ; \
  array_trans_var_v4_4 = array_trans_var_v4_4 / dpn_var_v4_1 ; \
  array_trans_var_v4_5 = array_trans_var_v4_5 / dpn_var_v4_1 ; \
  array_trans_var_v4_6 = array_trans_var_v4_6 / dpn_var_v4_1 ; \
  \
  array_trans_local_v4_ptr[0] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[1] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[4] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[5] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[32] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[33] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[48] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[49] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[32] ; \
  array_var_v4_1 = array_local_v4_ptr[36] ; \
  array_var_v4_2 = array_local_v4_ptr[40] ; \
  array_var_v4_3 = array_local_v4_ptr[44] ; \
  array_var_v4_4 = array_local_v4_ptr[48] ; \
  array_var_v4_5 = array_local_v4_ptr[52] ; \
  array_var_v4_6 = array_local_v4_ptr[56] ; \
  array_var_v4_7 = array_local_v4_ptr[60] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_var_v4_0 = dpn_local_v4_ptr[2] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[3] ; \
  \
  array_trans_var_v4_0 = array_trans_var_v4_0 / dpn_var_v4_0 ; \
  array_trans_var_v4_1 = array_trans_var_v4_1 / dpn_var_v4_0 ; \
  array_trans_var_v4_2 = array_trans_var_v4_2 / dpn_var_v4_0 ; \
  array_trans_var_v4_4 = array_trans_var_v4_4 / dpn_var_v4_1 ; \
  array_trans_var_v4_5 = array_trans_var_v4_5 / dpn_var_v4_1 ; \
  array_trans_var_v4_6 = array_trans_var_v4_6 / dpn_var_v4_1 ; \
  \
  array_trans_local_v4_ptr[2] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[3] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[6] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[7] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[34] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[35] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[50] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[51] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[1] ; \
  array_var_v4_1 = array_local_v4_ptr[5] ; \
  array_var_v4_2 = array_local_v4_ptr[9] ; \
  array_var_v4_3 = array_local_v4_ptr[13] ; \
  array_var_v4_4 = array_local_v4_ptr[17] ; \
  array_var_v4_5 = array_local_v4_ptr[21] ; \
  array_var_v4_6 = array_local_v4_ptr[25] ; \
  array_var_v4_7 = array_local_v4_ptr[29] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_var_v4_0 = dpn_local_v4_ptr[4] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[5] ; \
  \
  array_trans_var_v4_0 = array_trans_var_v4_0 / dpn_var_v4_0 ; \
  array_trans_var_v4_1 = array_trans_var_v4_1 / dpn_var_v4_0 ; \
  array_trans_var_v4_2 = array_trans_var_v4_2 / dpn_var_v4_0 ; \
  array_trans_var_v4_4 = array_trans_var_v4_4 / dpn_var_v4_1 ; \
  array_trans_var_v4_5 = array_trans_var_v4_5 / dpn_var_v4_1 ; \
  array_trans_var_v4_6 = array_trans_var_v4_6 / dpn_var_v4_1 ; \
  \
  array_trans_local_v4_ptr[8] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[9] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[12] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[13] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[36] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[37] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[52] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[53] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[33] ; \
  array_var_v4_1 = array_local_v4_ptr[37] ; \
  array_var_v4_2 = array_local_v4_ptr[41] ; \
  array_var_v4_3 = array_local_v4_ptr[45] ; \
  array_var_v4_4 = array_local_v4_ptr[49] ; \
  array_var_v4_5 = array_local_v4_ptr[53] ; \
  array_var_v4_6 = array_local_v4_ptr[57] ; \
  array_var_v4_7 = array_local_v4_ptr[61] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_var_v4_0 = dpn_local_v4_ptr[6] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[7] ; \
  \
  array_trans_var_v4_0 = array_trans_var_v4_0 / dpn_var_v4_0 ; \
  array_trans_var_v4_1 = array_trans_var_v4_1 / dpn_var_v4_0 ; \
  array_trans_var_v4_2 = array_trans_var_v4_2 / dpn_var_v4_0 ; \
  array_trans_var_v4_4 = array_trans_var_v4_4 / dpn_var_v4_1 ; \
  array_trans_var_v4_5 = array_trans_var_v4_5 / dpn_var_v4_1 ; \
  array_trans_var_v4_6 = array_trans_var_v4_6 / dpn_var_v4_1 ; \
  \
  array_trans_local_v4_ptr[10] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[11] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[14] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[15] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[38] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[39] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[54] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[55] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[2] ; \
  array_var_v4_1 = array_local_v4_ptr[6] ; \
  array_var_v4_2 = array_local_v4_ptr[10] ; \
  array_var_v4_3 = array_local_v4_ptr[14] ; \
  array_var_v4_4 = array_local_v4_ptr[18] ; \
  array_var_v4_5 = array_local_v4_ptr[22] ; \
  array_var_v4_6 = array_local_v4_ptr[26] ; \
  array_var_v4_7 = array_local_v4_ptr[30] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_var_v4_0 = dpn_local_v4_ptr[8] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[9] ; \
  \
  array_trans_var_v4_0 = array_trans_var_v4_0 / dpn_var_v4_0 ; \
  array_trans_var_v4_1 = array_trans_var_v4_1 / dpn_var_v4_0 ; \
  array_trans_var_v4_2 = array_trans_var_v4_2 / dpn_var_v4_0 ; \
  array_trans_var_v4_4 = array_trans_var_v4_4 / dpn_var_v4_1 ; \
  array_trans_var_v4_5 = array_trans_var_v4_5 / dpn_var_v4_1 ; \
  array_trans_var_v4_6 = array_trans_var_v4_6 / dpn_var_v4_1 ; \
  \
  array_trans_local_v4_ptr[16] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[17] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[20] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[21] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[40] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[41] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[56] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[57] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[34] ; \
  array_var_v4_1 = array_local_v4_ptr[38] ; \
  array_var_v4_2 = array_local_v4_ptr[42] ; \
  array_var_v4_3 = array_local_v4_ptr[46] ; \
  array_var_v4_4 = array_local_v4_ptr[50] ; \
  array_var_v4_5 = array_local_v4_ptr[54] ; \
  array_var_v4_6 = array_local_v4_ptr[58] ; \
  array_var_v4_7 = array_local_v4_ptr[62] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_var_v4_0 = dpn_local_v4_ptr[10] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[11] ; \
  \
  array_trans_var_v4_0 = array_trans_var_v4_0 / dpn_var_v4_0 ; \
  array_trans_var_v4_1 = array_trans_var_v4_1 / dpn_var_v4_0 ; \
  array_trans_var_v4_2 = array_trans_var_v4_2 / dpn_var_v4_0 ; \
  array_trans_var_v4_4 = array_trans_var_v4_4 / dpn_var_v4_1 ; \
  array_trans_var_v4_5 = array_trans_var_v4_5 / dpn_var_v4_1 ; \
  array_trans_var_v4_6 = array_trans_var_v4_6 / dpn_var_v4_1 ; \
  \
  array_trans_local_v4_ptr[18] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[19] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[22] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[23] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[42] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[43] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[58] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[59] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[3] ; \
  array_var_v4_1 = array_local_v4_ptr[7] ; \
  array_var_v4_2 = array_local_v4_ptr[11] ; \
  array_var_v4_3 = array_local_v4_ptr[15] ; \
  array_var_v4_4 = array_local_v4_ptr[19] ; \
  array_var_v4_5 = array_local_v4_ptr[23] ; \
  array_var_v4_6 = array_local_v4_ptr[27] ; \
  array_var_v4_7 = array_local_v4_ptr[31] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_var_v4_0 = dpn_local_v4_ptr[12] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[13] ; \
  \
  array_trans_var_v4_0 = array_trans_var_v4_0 / dpn_var_v4_0 ; \
  array_trans_var_v4_1 = array_trans_var_v4_1 / dpn_var_v4_0 ; \
  array_trans_var_v4_2 = array_trans_var_v4_2 / dpn_var_v4_0 ; \
  array_trans_var_v4_4 = array_trans_var_v4_4 / dpn_var_v4_1 ; \
  array_trans_var_v4_5 = array_trans_var_v4_5 / dpn_var_v4_1 ; \
  array_trans_var_v4_6 = array_trans_var_v4_6 / dpn_var_v4_1 ; \
  \
  array_trans_local_v4_ptr[24] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[25] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[28] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[29] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[44] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[45] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[60] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[61] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[35] ; \
  array_var_v4_1 = array_local_v4_ptr[39] ; \
  array_var_v4_2 = array_local_v4_ptr[43] ; \
  array_var_v4_3 = array_local_v4_ptr[47] ; \
  array_var_v4_4 = array_local_v4_ptr[51] ; \
  array_var_v4_5 = array_local_v4_ptr[55] ; \
  array_var_v4_6 = array_local_v4_ptr[59] ; \
  array_var_v4_7 = array_local_v4_ptr[63] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  dpn_var_v4_0 = dpn_local_v4_ptr[14] ; \
  dpn_var_v4_1 = dpn_local_v4_ptr[15] ; \
  \
  array_trans_var_v4_0 = array_trans_var_v4_0 / dpn_var_v4_0 ; \
  array_trans_var_v4_1 = array_trans_var_v4_1 / dpn_var_v4_0 ; \
  array_trans_var_v4_2 = array_trans_var_v4_2 / dpn_var_v4_0 ; \
  array_trans_var_v4_4 = array_trans_var_v4_4 / dpn_var_v4_1 ; \
  array_trans_var_v4_5 = array_trans_var_v4_5 / dpn_var_v4_1 ; \
  array_trans_var_v4_6 = array_trans_var_v4_6 / dpn_var_v4_1 ; \
  \
  array_trans_local_v4_ptr[26] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[27] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[30] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[31] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[46] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[47] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[62] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[63] = array_trans_var_v4_7 ; \
}

#define TRANSPOSE_ARRAY_POST_QDP_2X4X4(array_trans_local_v4_ptr , array_local_v4_ptr) { \
  doublev4 array_var_v4_0 , array_var_v4_1 , array_var_v4_2 , array_var_v4_3 ; \
  doublev4 array_var_v4_4 , array_var_v4_5 , array_var_v4_6 , array_var_v4_7 ; \
  array_var_v4_0 = array_local_v4_ptr[0] ; \
  array_var_v4_1 = array_local_v4_ptr[4] ; \
  array_var_v4_2 = array_local_v4_ptr[8] ; \
  array_var_v4_3 = array_local_v4_ptr[12] ; \
  array_var_v4_4 = array_local_v4_ptr[16] ; \
  array_var_v4_5 = array_local_v4_ptr[20] ; \
  array_var_v4_6 = array_local_v4_ptr[24] ; \
  array_var_v4_7 = array_local_v4_ptr[28] ; \
  \
  doublev4 dtmpv4_0 , dtmpv4_1 , dtmpv4_2 , dtmpv4_3 ; \
  doublev4 dtmpv4_4 , dtmpv4_5 , dtmpv4_6 , dtmpv4_7 ; \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  doublev4 array_trans_var_v4_0 , array_trans_var_v4_1 , array_trans_var_v4_2 , array_trans_var_v4_3 ; \
  doublev4 array_trans_var_v4_4 , array_trans_var_v4_5 , array_trans_var_v4_6 , array_trans_var_v4_7 ; \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  array_trans_local_v4_ptr[0] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[1] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[16] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[17] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[32] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[33] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[48] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[49] = array_trans_var_v4_7 ; \
  \
  array_var_v4_0 = array_local_v4_ptr[32] ; \
  array_var_v4_1 = array_local_v4_ptr[36] ; \
  array_var_v4_2 = array_local_v4_ptr[40] ; \
  array_var_v4_3 = array_local_v4_ptr[44] ; \
  array_var_v4_4 = array_local_v4_ptr[48] ; \
  array_var_v4_5 = array_local_v4_ptr[52] ; \
  array_var_v4_6 = array_local_v4_ptr[56] ; \
  array_var_v4_7 = array_local_v4_ptr[60] ; \
  \
  VSHFF(dtmpv4_0 , array_var_v4_1 , array_var_v4_0 , 0x44) ; \
  VSHFF(dtmpv4_1 , array_var_v4_1 , array_var_v4_0 , 0xee) ; \
  VSHFF(dtmpv4_2 , array_var_v4_3 , array_var_v4_2 , 0x44) ; \
  VSHFF(dtmpv4_3 , array_var_v4_3 , array_var_v4_2 , 0xee) ; \
  VSHFF(dtmpv4_4 , array_var_v4_5 , array_var_v4_4 , 0x44) ; \
  VSHFF(dtmpv4_5 , array_var_v4_5 , array_var_v4_4 , 0xee) ; \
  VSHFF(dtmpv4_6 , array_var_v4_7 , array_var_v4_6 , 0x44) ; \
  VSHFF(dtmpv4_7 , array_var_v4_7 , array_var_v4_6 , 0xee) ; \
  \
  VSHFF(array_trans_var_v4_0 , dtmpv4_2 , dtmpv4_0 , 0x88) ; \
  VSHFF(array_trans_var_v4_1 , dtmpv4_2 , dtmpv4_0 , 0xdd) ; \
  VSHFF(array_trans_var_v4_2 , dtmpv4_3 , dtmpv4_1 , 0x88) ; \
  VSHFF(array_trans_var_v4_3 , dtmpv4_3 , dtmpv4_1 , 0xdd) ; \
  VSHFF(array_trans_var_v4_4 , dtmpv4_6 , dtmpv4_4 , 0x88) ; \
  VSHFF(array_trans_var_v4_5 , dtmpv4_6 , dtmpv4_4 , 0xdd) ; \
  VSHFF(array_trans_var_v4_6 , dtmpv4_7 , dtmpv4_5 , 0x88) ; \
  VSHFF(array_trans_var_v4_7 , dtmpv4_7 , dtmpv4_5 , 0xdd) ; \
  \
  array_trans_local_v4_ptr[2] = array_trans_var_v4_0 ; \
  array_trans_local_v4_ptr[3] = array_trans_var_v4_4 ; \
  array_trans_local_v4_ptr[18] = array_trans_var_v4_1 ; \
  array_trans_local_v4_ptr[19] = array_trans_var_v4_5 ; \
  array_trans_local_v4_ptr[34] = array_trans_var_v4_2 ; \
  array_trans_local_v4_ptr[35] = array_trans_var_v4_6 ; \
  array_trans_local_v4_ptr[50] = array_trans_var_v4_3 ; \
  array_trans_local_v4_ptr[51] = array_trans_var_v4_7 ; \
  \
}
#define TRANSPOSE_ARRAY_POST_QDP_4X4X4X4(array_trans_local_v4_ptr , array_local_v4_ptr) {\
  doublev4 *array_trans_local_v4_ptr_0 , *array_trans_local_v4_ptr_1 , *array_trans_local_v4_ptr_2 , *array_trans_local_v4_ptr_3 ; \
  doublev4 *array_local_v4_ptr_0 , *array_local_v4_ptr_1 , *array_local_v4_ptr_2 , *array_local_v4_ptr_3 ; \
  array_trans_local_v4_ptr_0 = &array_trans_local_v4_ptr[0] ; \
  array_trans_local_v4_ptr_1 = &array_trans_local_v4_ptr[4] ; \
  array_trans_local_v4_ptr_2 = &array_trans_local_v4_ptr[8] ; \
  array_trans_local_v4_ptr_3 = &array_trans_local_v4_ptr[12] ; \
  \
  array_local_v4_ptr_0 = &array_local_v4_ptr[0] ; \
  array_local_v4_ptr_1 = &array_local_v4_ptr[1] ; \
  array_local_v4_ptr_2 = &array_local_v4_ptr[2] ; \
  array_local_v4_ptr_3 = &array_local_v4_ptr[3] ; \
  \
  TRANSPOSE_ARRAY_POST_QDP_2X4X4(array_trans_local_v4_ptr_0 , array_local_v4_ptr_0) ; \
  TRANSPOSE_ARRAY_POST_QDP_2X4X4(array_trans_local_v4_ptr_1 , array_local_v4_ptr_1) ; \
  TRANSPOSE_ARRAY_POST_QDP_2X4X4(array_trans_local_v4_ptr_2 , array_local_v4_ptr_2) ; \
  TRANSPOSE_ARRAY_POST_QDP_2X4X4(array_trans_local_v4_ptr_3 , array_local_v4_ptr_3) ; \
}


void  TRANSPOSE_ARRAY_FORWARD(remap_q_ppm_param_t *remap_q_ppm_ptr) {
  doublev4 *array_local_v4_ptr ;
  doublev4 *array_trans_local_v4_ptr ;
  doublev4 *dp3d_local_v4_ptr ;
  intv8 send_id_v8 ;
  intv8 send_offset_v8 ;
  intv8 recv_offset_v8 ;
  send_id_v8 = remap_q_ppm_ptr->send_id_v8 ;
  send_offset_v8 = remap_q_ppm_ptr->send_offset_v8_forward_qdp ;
  recv_offset_v8 = remap_q_ppm_ptr->recv_offset_v8_forward_qdp ;
  array_local_v4_ptr = remap_q_ppm_ptr->array_local_v4_ptr ;
  array_trans_local_v4_ptr = remap_q_ppm_ptr->array_trans_local_v4_ptr ;
  dp3d_local_v4_ptr = remap_q_ppm_ptr->dp3d_local_v4_ptr ;
  TRANSPOSE_ARRAY_FORWARD_4X4X4X4(array_trans_local_v4_ptr , array_local_v4_ptr , dp3d_local_v4_ptr) ;
  REGISTER_FULL_COMMUICATION_FORWARD(array_local_v4_ptr , array_trans_local_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8) ;
}
void TRANSPOSE_ARRAY_FORWARD_QDP(remap_q_ppm_param_t *remap_q_ppm_ptr) {
  doublev4 *array_local_v4_ptr ;
  doublev4 *array_trans_local_v4_ptr ;
  intv8 send_id_v8 ;
  intv8 send_offset_v8 ;
  intv8 recv_offset_v8 ;
  send_id_v8 = remap_q_ppm_ptr->send_id_v8 ;
  send_offset_v8 = remap_q_ppm_ptr->send_offset_v8_forward_qdp ;
  recv_offset_v8 = remap_q_ppm_ptr->recv_offset_v8_forward_qdp ;
  array_local_v4_ptr = remap_q_ppm_ptr->array_local_v4_ptr ;
  array_trans_local_v4_ptr = remap_q_ppm_ptr->array_trans_local_v4_ptr ;
  TRANSPOSE_ARRAY_FORWARD_QDP_4X4X4X4(array_trans_local_v4_ptr , array_local_v4_ptr) ;
  REGISTER_FULL_COMMUICATION_FORWARD(array_local_v4_ptr , array_trans_local_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8) ;
}


void REMAP_Q_PPM(remap_q_ppm_param_t *remap_q_ppm_ptr , int ie , int q) {
  doublev4 *array_local_v4_ptr ;
  doublev4 *array_trans_local_v4_ptr ;
  doublev4 *masso_local_v4_ptr ;
  doublev4 *ppmdx_local_v4_ptr ; 
  doublev4 *coefs_local_v4_ptr ;
  doublev4 *ai_local_v4_ptr ;
  doublev4 *dma_local_v4_ptr ;
  doublev4 *dpo_local_v4_ptr ;
  double *dpo_local_ptr ;
  double *z_local_ptr ;
  int *kid_local_ptr ;
  array_local_v4_ptr = remap_q_ppm_ptr->array_local_v4_ptr ; 
  array_trans_local_v4_ptr = remap_q_ppm_ptr->array_trans_local_v4_ptr ;
  masso_local_v4_ptr = remap_q_ppm_ptr->masso_local_v4_ptr ;
  ppmdx_local_v4_ptr = remap_q_ppm_ptr->ppmdx_local_v4_ptr ;
  coefs_local_v4_ptr = remap_q_ppm_ptr->coefs_local_v4_ptr ;
  ai_local_v4_ptr = remap_q_ppm_ptr->ai_local_v4_ptr ;
  dma_local_v4_ptr = remap_q_ppm_ptr->dma_local_v4_ptr ;
  dpo_local_v4_ptr = remap_q_ppm_ptr->dpo_local_v4_ptr ; 
  dpo_local_ptr = (double *)remap_q_ppm_ptr->dpo_local_v4_ptr ;
  z_local_ptr = remap_q_ppm_ptr->z_local_ptr ; 
  kid_local_ptr =  remap_q_ppm_ptr->kid_local_ptr ;
  COMPUTE_AO_MASSO(masso_local_v4_ptr , array_local_v4_ptr , dpo_local_v4_ptr) ;
  MEMCOPY_2X2X4(array_local_v4_ptr) ;
  COMPUTE_PPM(coefs_local_v4_ptr , ai_local_v4_ptr , array_local_v4_ptr , array_trans_local_v4_ptr , ppmdx_local_v4_ptr , dma_local_v4_ptr , ie , q ) ;
  COMPUTE_MASSN(array_local_v4_ptr , coefs_local_v4_ptr , masso_local_v4_ptr , z_local_ptr , kid_local_ptr , dpo_local_ptr , ie , q) ;
  array_local_v4_ptr = &array_local_v4_ptr[34] ; 
  ppmdx_local_v4_ptr = &ppmdx_local_v4_ptr[8] ;
  dpo_local_v4_ptr = &dpo_local_v4_ptr[9] ; 
  dpo_local_ptr = (double *)dpo_local_v4_ptr ; 
  z_local_ptr = &z_local_ptr[32] ; 
  kid_local_ptr = &kid_local_ptr[32] ;
  COMPUTE_AO_MASSO(masso_local_v4_ptr , array_local_v4_ptr , dpo_local_v4_ptr) ;
  MEMCOPY_2X2X4(array_local_v4_ptr) ;
  COMPUTE_PPM(coefs_local_v4_ptr , ai_local_v4_ptr , array_local_v4_ptr , array_trans_local_v4_ptr , ppmdx_local_v4_ptr , dma_local_v4_ptr , ie , q) ;
  COMPUTE_MASSN(array_local_v4_ptr , coefs_local_v4_ptr , masso_local_v4_ptr , z_local_ptr , kid_local_ptr , dpo_local_ptr , ie ,q) ;

}

void TRANSPOSE_ARRAY_POST(remap_q_ppm_param_t *remap_q_ppm_ptr) {
  doublev4 *array_local_v4_ptr ;
  doublev4 *array_trans_local_v4_ptr ;
  doublev4 *dpn_local_v4_ptr ;
  intv8 send_id_v8 ;
  intv8 send_offset_v8 ; 
  intv8 recv_offset_v8 ; 
  array_local_v4_ptr = remap_q_ppm_ptr->array_local_v4_ptr ; 
  array_trans_local_v4_ptr = remap_q_ppm_ptr->array_trans_local_v4_ptr ;
  dpn_local_v4_ptr = remap_q_ppm_ptr->dpn_local_v4_ptr ;
  send_id_v8 = remap_q_ppm_ptr->send_id_v8 ; 
  send_offset_v8 = remap_q_ppm_ptr->send_offset_v8_post_qdp ; 
  recv_offset_v8 = remap_q_ppm_ptr->recv_offset_v8_post_qdp ;
  REGISTER_FULL_COMMUICATION_POST(array_trans_local_v4_ptr , array_local_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8) ;
  TRANSPOSE_ARRAY_POST_4X4X4X4(array_local_v4_ptr , array_trans_local_v4_ptr , dpn_local_v4_ptr) ;
 
}

void TRANSPOSE_ARRAY_POST_QDP(remap_q_ppm_param_t *remap_q_ppm_ptr) {
  doublev4 *array_local_v4_ptr ;
  doublev4 *array_trans_local_v4_ptr ;
  intv8 send_id_v8 ;
  intv8 send_offset_v8 ;
  intv8 recv_offset_v8 ;
  array_local_v4_ptr = remap_q_ppm_ptr->array_local_v4_ptr ;
  array_trans_local_v4_ptr = remap_q_ppm_ptr->array_trans_local_v4_ptr ;
  send_id_v8 = remap_q_ppm_ptr->send_id_v8 ;
  send_offset_v8 = remap_q_ppm_ptr->send_offset_v8_post_qdp ;
  recv_offset_v8 = remap_q_ppm_ptr->recv_offset_v8_post_qdp ;
  REGISTER_FULL_COMMUICATION_POST(array_trans_local_v4_ptr , array_local_v4_ptr , send_id_v8 , send_offset_v8 , recv_offset_v8) ;
  TRANSPOSE_ARRAY_POST_QDP_4X4X4X4(array_local_v4_ptr , array_trans_local_v4_ptr) ;
}
