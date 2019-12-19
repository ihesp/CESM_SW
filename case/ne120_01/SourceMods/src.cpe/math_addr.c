/*************************************************************************
	> File Name: SourceMods/src.cpe/math_addr.c
	> Author: Xu Kai
	> Created Time: 2019年01月16日 星期三 12时14分29秒
 ************************************************************************/
	
	
#include<stdio.h>
#include<stdlib.h>
#include "dma_macros.h"
#include "cpe_print.h"
#include "math_data.h"
void math_addr_(){
    if(_MYID == 0)
        cpe_printf("%p\n", log_data);

    double log_local_data[log_data_len] ;
    double log10_local_data[log10_data_len] ;
    double pow_local_data[pow_data_len] ;

    //log_data_local_ptr = log_local_data ;
    //log10_data_local_ptr = log10_local_data ;
    //pow_data_local_ptr = pow_local_data ;
    dma_init();
    if (_MYID == 0) {
        
        bcast_get(log_data , log_local_data , log_data_bytes) ;
        bcast_get(log10_data , log10_local_data , log10_data_bytes) ;
        bcast_get(pow_data , pow_local_data , pow_data_bytes) ;
        dma_syn() ;
    }
    athread_syn(ARRAY_SCOPE ,0xffff) ;
    //pe_get(log_data, log_local_data, log_data_bytes);
    //dma_syn() ;
    //memcpy(log_local_data, log_data, log_data_bytes);

}

