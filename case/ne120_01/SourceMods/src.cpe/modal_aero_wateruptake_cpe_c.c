/*************************************************************************
	> Created Time: 2019年01月08日 星期二 15时22分23秒
 ************************************************************************/
	
	
#include<stdio.h>
#include<stdlib.h>


#include"dma_macros.h"
#include<stdbool.h>
#include<stdarg.h>
#include "cpe_print.h"
#include"math_data.h"


#define REAL double
#define L1 16
#define L2 30

#define COL(x) (x & 0x07)
#define ROW(x) ((x & 0x38) >> 3)

typedef struct modal_aero_wateruptake_sub_para_c
{
    REAL pi                 ;
    REAL third              ;
    REAL pi43               ;
    int rank                ;
    int pcols               ; 
    int pver                ; 
    int top_lev             ;   
    int nmodes              ;   
    int ncol                ;    
    int modal_strat_sulfate ;
    REAL so4specdens        ;
    void * troplev          ;
    void * rhcrystal        ;
    void * rhdeliques       ;
    void * dryrad           ;        
    void * hygro            ;
    void * dryvol           ;
    void * so4dryvol        ;
    void * wtpct            ;
    void * sulden           ;
    void * wetrad           ;
    void * wetvol           ;
    void * wtrvol           ;
    void * rh               ; 
} modal_aero_wateruptake_sub_args_cc;

void pow_zzf_parallel_(void *ptr)
{
    if ( _MYID == 0 )
        pow_zzf_(ptr);
}
void get_double_(double *in, double *out){
	*out = in[1];
}

void modal_aero_wateruptake_sub_parallel_(void *modal_aero_wateruptake_sub_para)
{
    //set_func_ptr_();
    dma_init();
    modal_aero_wateruptake_sub_args_cc modal_aero_wateruptake_sub_spara;

    pe_get(modal_aero_wateruptake_sub_para, &modal_aero_wateruptake_sub_spara, sizeof(modal_aero_wateruptake_sub_args_cc));
    dma_syn();

    // in
    int  troplev   ; // pcols
    REAL rhcrystal [L1]; // nmodes
    REAL rhdeliques[L1]; // nmodes
    REAL dryrad    [L1][L2]; //nmodes, pver
    REAL hygro     [L1][L2]; 
    REAL dryvol    [L1][L2]; 
    REAL so4dryvol [L1][L2];
    REAL wtpct     [L1][L2]; 
    REAL sulden    [L1][L2]; 
    REAL rh        [L2];
    // out
    REAL wetrad    [L1][L2]; 
    REAL wetvol    [L1][L2]; 
    REAL wtrvol    [L1][L2]; 

    int rank                 = modal_aero_wateruptake_sub_spara.rank               ; 
    int pcols                = modal_aero_wateruptake_sub_spara.pcols              ; 
    int pver                 = modal_aero_wateruptake_sub_spara.pver               ; 
    int top_lev              = modal_aero_wateruptake_sub_spara.top_lev            ; 
    int nmodes               = modal_aero_wateruptake_sub_spara.nmodes             ; 
    int ncol                 = modal_aero_wateruptake_sub_spara.ncol               ; 
    int modal_strat_sulfate  = modal_aero_wateruptake_sub_spara.modal_strat_sulfate; 
    REAL so4specdens         = modal_aero_wateruptake_sub_spara.so4specdens;
    REAL pi                  = modal_aero_wateruptake_sub_spara.pi   ;
    REAL third               = modal_aero_wateruptake_sub_spara.third;
    REAL pi43                = modal_aero_wateruptake_sub_spara.pi43 ;

    //if ( rank == 0 && _MYID == 0 ) cpe_printf("CPE: %d %d %d %d %d %d %.20lf %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n", pcols, pver, top_lev, nmodes, ncol, modal_strat_sulfate, so4specdens,modal_aero_wateruptake_sub_spara.troplev,modal_aero_wateruptake_sub_spara.rhcrystal,modal_aero_wateruptake_sub_spara.rhdeliques,modal_aero_wateruptake_sub_spara.dryrad   ,modal_aero_wateruptake_sub_spara.hygro    ,modal_aero_wateruptake_sub_spara.dryvol   ,modal_aero_wateruptake_sub_spara.so4dryvol,modal_aero_wateruptake_sub_spara.wtpct    ,modal_aero_wateruptake_sub_spara.sulden   ,modal_aero_wateruptake_sub_spara.wetrad,modal_aero_wateruptake_sub_spara.wetvol,modal_aero_wateruptake_sub_spara.wtrvol,modal_aero_wateruptake_sub_spara.rh);

    int myid = COL(_MYID)*8 + ROW(_MYID);

    int m, k;
    troplev = ((int*)(modal_aero_wateruptake_sub_spara.troplev))[myid];
    //pe_get(modal_aero_wateruptake_sub_spara.troplev   , troplev   , pcols*sizeof(int));
    pe_get(modal_aero_wateruptake_sub_spara.rhcrystal , rhcrystal , nmodes*sizeof(REAL));
    pe_get(modal_aero_wateruptake_sub_spara.rhdeliques, rhdeliques, nmodes*sizeof(REAL));
    dma_syn();

    double sincos_local_data[sincos_data_len] ;
    double exp_local_data[exp_data_len] ;
    //double asin_local_data[asin_data_len] ;
    //double acos_local_data[acos_data_len] ;
    double log_local_data[log_data_len] ;
    double log10_local_data[log10_data_len] ;
    double pow_local_data[pow_data_len] ; 

    cos_data_local_ptr = sincos_local_data ;
    sin_data_local_ptr = sincos_local_data ; 
    exp_data_local_ptr = exp_local_data ;
    //asin_data_local_ptr = asin_local_data ;
    //acos_data_local_ptr = acos_local_data ;
    log_data_local_ptr = log_local_data ; 
    log10_data_local_ptr = log10_local_data ;  
    pow_data_local_ptr = pow_local_data ; 

    if (_MYID == 0) { 
        bcast_get(sincos_data , sincos_local_data , sincos_data_bytes) ; 
        bcast_get(exp_data , exp_local_data , exp_data_bytes) ;
        //bcast_get(asin_data , asin_local_data , asin_data_bytes) ;
        //bcast_get(acos_data , acos_local_data , acos_data_bytes) ; 
        bcast_get(log_data , log_local_data , log_data_bytes) ;
        bcast_get(log10_data , log10_local_data , log10_data_bytes) ;
        bcast_get(pow_data , pow_local_data , pow_data_bytes) ;
        dma_syn() ; 
    } 
    athread_syn(ARRAY_SCOPE ,0xffff) ;


#define OFFSET(dim1, dim2) (myid * (dim1) * (dim2) * sizeof(REAL))
    if ( myid < ncol )
    {
        pe_get(modal_aero_wateruptake_sub_spara.dryrad    + OFFSET(nmodes, pver), dryrad    , nmodes*pver*sizeof(REAL));
        pe_get(modal_aero_wateruptake_sub_spara.hygro     + OFFSET(nmodes, pver), hygro     , nmodes*pver*sizeof(REAL));
        pe_get(modal_aero_wateruptake_sub_spara.dryvol    + OFFSET(nmodes, pver), dryvol    , nmodes*pver*sizeof(REAL));
        pe_get(modal_aero_wateruptake_sub_spara.so4dryvol + OFFSET(nmodes, pver), so4dryvol , nmodes*pver*sizeof(REAL));
        pe_get(modal_aero_wateruptake_sub_spara.wtpct     + OFFSET(nmodes, pver), wtpct     , nmodes*pver*sizeof(REAL));
        pe_get(modal_aero_wateruptake_sub_spara.sulden    + OFFSET(nmodes, pver), sulden    , nmodes*pver*sizeof(REAL));
        pe_get(modal_aero_wateruptake_sub_spara.rh        + OFFSET(1, pver)     , rh        , pver*sizeof(REAL));
        //pe_get(modal_aero_wateruptake_sub_spara.wetrad    + OFFSET(nmodes, pver), wetrad    , nmodes*pver*sizeof(REAL));
        //pe_get(modal_aero_wateruptake_sub_spara.wetvol    + OFFSET(nmodes, pver), wetvol    , nmodes*pver*sizeof(REAL));
        //pe_get(modal_aero_wateruptake_sub_spara.wtrvol    + OFFSET(nmodes, pver), wtrvol    , nmodes*pver*sizeof(REAL));

        dma_syn();

        //for ( k = 0; k < pver; k ++ )
        //{
        //    for ( m = 0; m < nmodes; m ++ )
        //    {
        //        if ( abs(dryrad   [m][k] - ((REAL*)modal_aero_wateruptake_sub_spara.dryrad   )[myid*nmodes*pver + m*pver + k]) > 1e-15 ) cpe_printf("dryrad   wrong, %.20lf, %.20lf\n", dryrad   [m][k], ((REAL*)modal_aero_wateruptake_sub_spara.dryrad   )[myid*nmodes*pver + m*pver + k]);
        //        if ( abs(hygro    [m][k] - ((REAL*)modal_aero_wateruptake_sub_spara.hygro    )[myid*nmodes*pver + m*pver + k]) > 1e-15 ) cpe_printf("hygro    wrong, %.20lf, %.20lf\n", hygro    [m][k], ((REAL*)modal_aero_wateruptake_sub_spara.hygro    )[myid*nmodes*pver + m*pver + k]);
        //        if ( abs(dryvol   [m][k] - ((REAL*)modal_aero_wateruptake_sub_spara.dryvol   )[myid*nmodes*pver + m*pver + k]) > 1e-15 ) cpe_printf("dryvol   wrong, %.20lf, %.20lf\n", dryvol   [m][k], ((REAL*)modal_aero_wateruptake_sub_spara.dryvol   )[myid*nmodes*pver + m*pver + k]);
        //        if ( abs(so4dryvol[m][k] - ((REAL*)modal_aero_wateruptake_sub_spara.so4dryvol)[myid*nmodes*pver + m*pver + k]) > 1e-15 ) cpe_printf("so4dryvolwrong, %.20lf, %.20lf\n", so4dryvol[m][k], ((REAL*)modal_aero_wateruptake_sub_spara.so4dryvol)[myid*nmodes*pver + m*pver + k]);
        //        if ( abs(wtpct    [m][k] - ((REAL*)modal_aero_wateruptake_sub_spara.wtpct    )[myid*nmodes*pver + m*pver + k]) > 1e-15 ) cpe_printf("wtpct    wrong, %.20lf, %.20lf\n", wtpct    [m][k], ((REAL*)modal_aero_wateruptake_sub_spara.wtpct    )[myid*nmodes*pver + m*pver + k]);
        //        if ( abs(sulden   [m][k] - ((REAL*)modal_aero_wateruptake_sub_spara.sulden   )[myid*nmodes*pver + m*pver + k]) > 1e-15 ) cpe_printf("sulden   wrong, %.20lf, %.20lf\n", sulden   [m][k], ((REAL*)modal_aero_wateruptake_sub_spara.sulden   )[myid*nmodes*pver + m*pver + k]);
        //    }
        //    if ( abs(rh       [k] - ((REAL*)modal_aero_wateruptake_sub_spara.rh       )[myid*pver + k]) > 1e-15 ) cpe_printf("rh      wrong, %.20lf, %.20lf\n", rh[k], ((REAL*)modal_aero_wateruptake_sub_spara.rh       )[myid*pver + k]);
        //}

        modal_aero_wateruptake_sub_compute_(\
                &myid, &pcols, &pver, &top_lev, &modal_strat_sulfate, \
                &ncol, &nmodes, rhcrystal, rhdeliques, dryrad, \
                hygro, rh, dryvol, so4dryvol, &so4specdens, &troplev, \
                wetrad, wetvol, wtrvol, sulden, wtpct, &pi, &third, &pi43);

        pe_put(modal_aero_wateruptake_sub_spara.wetrad    + OFFSET(nmodes, pver), wetrad    , nmodes*pver*sizeof(REAL));
        pe_put(modal_aero_wateruptake_sub_spara.wetvol    + OFFSET(nmodes, pver), wetvol    , nmodes*pver*sizeof(REAL));
        pe_put(modal_aero_wateruptake_sub_spara.wtrvol    + OFFSET(nmodes, pver), wtrvol    , nmodes*pver*sizeof(REAL));
        dma_syn();
    }

    cos_data_local_ptr = sincos_data ;
    sin_data_local_ptr = sincos_data ; 
    exp_data_local_ptr = exp_data ;
    asin_data_local_ptr = asin_data ;
    acos_data_local_ptr = acos_data ;
    log_data_local_ptr = log_data ; 
    log10_data_local_ptr = log10_data ;  
    pow_data_local_ptr = pow_data ; 

}
