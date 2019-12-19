/*************************************************************************
	> File Name: micro_mg1_0_cpe_c.c
	> Author: Xu Kai
	> Created Time: 2018年11月28日 星期三 15时17分35秒
 ************************************************************************/
	
	
#define CPE
#define LWPF_KERNELS K(A) K(B) K(C) K(D)
#define LWPF_UNIT U(TEST)

#include<stdio.h>
#include<stdlib.h>
#include"dma_macros.h"
#include<stdbool.h>
#include<stdarg.h>
#include "cpe_print.h"
#include"/home/export/online1/cesm06/dxh/workspace/ldm_math/math_data.h"

#include "/home/export/online1/swmore/opensource/lwpf2/lwpf2.h"  
//#include "/home/export/online1/swmore/opensource/lwpf2/lwpf2_undef.h"  
//#include "/home/export/online1/swmore/opensource/lwpf2/pcrdef.h"  
//#include "/home/export/online1/swmore/opensource/lwpf2/sw5_pcr.h"  


#define REAL double
#define LAYER 15

#define COL(x) (x & 0x07)
#define ROW(x) ((x & 0x38) >> 3)

//extern "C" double __wv_sat_methods_cpe__shr_spfn_gamma_nonintrinsic_r8(double *x);
void  get_sp_(long long int *spa){
    long long int spt;
    asm volatile("mov $sp, %0 \n\t": "=r"(spt));
    *spa = spt; 
}

typedef struct micro_mg1_0_para_c
{
    REAL   svp_tmelt          ;           
    REAL   svp_h2otrip        ; 
    REAL   svp_tboil          ; 
    REAL   svp_ttrice         ; 
    REAL   svp_epsilo         ;  
    REAL   g                  ;  
    REAL   r                  ;
    REAL   rv                 ;
    REAL   cpp                ;
    REAL   rhow               ;
    REAL   tmelt              ;
    REAL   xxlv               ;
    REAL   xlf                ;
    REAL   xxls               ;
    REAL   rhosn              ;
    REAL   rhoi               ;
    REAL   ac                 ;
    REAL   bc                 ;
    REAL   as                 ;
    REAL   bs                 ;
    REAL   ai                 ;
    REAL   bi                 ;
    REAL   ar                 ;
    REAL   br                 ;
    REAL   ci                 ;
    REAL   di                 ;
    REAL   cs                 ;
    REAL   ds                 ;
    REAL   cr                 ;
    REAL   dr                 ;
    REAL   f1s                ;
    REAL   f2s                ;
    REAL   Eii                ;
    REAL   Ecr                ;
    REAL   f1r                ;
    REAL   f2r                ;
    REAL   DCS                ;
    REAL   qsmall             ;
    REAL   bimm               ;
    REAL   aimm               ;
    REAL   rhosu              ;
    REAL   mi0                ;
    REAL   rin                ;
    REAL   pi                 ;
    REAL   cons1              ;
    REAL   cons4              ;
    REAL   cons5              ;
    REAL   cons6              ;
    REAL   cons7              ;
    REAL   cons8              ;
    REAL   cons11             ;
    REAL   cons13             ;
    REAL   cons14             ;
    REAL   cons16             ;
    REAL   cons17             ;
    REAL   cons22             ;
    REAL   cons23             ;
    REAL   cons24             ;
    REAL   cons25             ;
    REAL   cons27             ;
    REAL   cons28             ;
    REAL   lammini            ;
    REAL   lammaxi            ;
    REAL   lamminr            ;
    REAL   lammaxr            ;
    REAL   lammins            ;
    REAL   lammaxs            ;
    REAL   tmax_fsnow         ;
    REAL   tmin_fsnow         ;
    REAL   tt0                ;
    REAL   csmin              ;
    REAL   csmax              ;
    REAL   minrefl            ;
    REAL   mindbz             ;
    REAL   rhmini             ;
    int    use_hetfrz_classnuc;

    int microp_uniform;
    int do_cldice;
    int pcols;
    int pver;
    int ncol;
    int top_lev;
    REAL deltatin;

    void *tn                        ;  
    void *qn                        ;
    void *relvar                    ;
    void *accre_enhan               ;
    void *p                         ;
    void *pdel                      ;
    void *cldn                      ;
    void *icecldf                   ;
    void *liqcldf                   ;
    void *naai                      ;
    void *npccnin                   ;
    void *rndst                     ;
    void *nacon                     ;

    void *tnd_qsnow                 ;
    void *tnd_nsnow                 ;
    void *re_ice                    ;
    void *frzimm                    ;
    void *frzcnt                    ;
    void *frzdep                    ;

    void *rate1ord_cw2pr_st         ;
    void *tlat                      ;
    void *qvlat                     ;
    void *qctend                    ;
    void *qitend                    ;
    void *nctend                    ;
    void *nitend                    ;
    void *effc                      ;
    void *effc_fn                   ;
    void *effi                      ;
    void *prect                     ;
    void *preci                     ;
    void *nevapr                    ;
    void *evapsnow                  ;
    void *am_evp_st                 ;
    void *prain                     ;
    void *prodsnow                  ;
    void *cmeout                    ;
    void *deffi                     ;
    void *pgamrad                   ;
    void *lamcrad                   ;
    void *qsout                     ;
    void *dsout                     ;
    void *rflx                      ;
    void *sflx                      ;
    void *qrout                     ;
    void *qcsevap                   ;
    void *qisevap                   ;
    void *qvres                     ;
    void *cmeiout                   ;
    void *vtrmc                     ;
    void *vtrmi                     ;
    void *qcsedten                  ;
    void *qisedten                  ;
    void *prao                      ;
    void *prco                      ;
    void *mnuccco                   ;
    void *mnuccto                   ;
    void *msacwio                   ;
    void *psacwso                   ;
    void *bergso                    ;
    void *bergo                     ;
    void *melto                     ;
    void *homoo                     ;
    void *qcreso                    ;
    void *prcio                     ;
    void *praio                     ;
    void *qireso                    ;
    void *mnuccro                   ;
    void *pracso                    ;
    void *meltsdt                   ;
    void *frzrdt                    ;
    void *mnuccdo                   ;
    void *nrout                     ;
    void *nsout                     ;
    void *refl                      ;
    void *arefl                     ;
    void *areflz                    ;
    void *frefl                     ;
    void *csrfl                     ;
    void *acsrfl                    ;
    void *fcsrfl                    ;
    void *rercld                    ;
    void *ncai                      ;
    void *ncal                      ;
    void *qrout2                    ;
    void *qsout2                    ;
    void *nrout2                    ;
    void *nsout2                    ;
    void *drout2                    ;
    void *dsout2                    ;
    void *freqs                     ;
    void *freqr                     ;
    void *nfice                     ;
    void *prer_evap                 ;
    void *qc                        ;
    void *qi                        ;
    void *nc                        ;
    void *ni                        ;
    void *reff_rain                 ;
    void *reff_snow                 ;
} micro_mg1_0_args_cc;

void log10_parallel_(REAL *log10_para)
{
    if ( _MYID == 0 ) 
    {
        long long int func_log10;
        REAL res;
        math_agent_log10_c_(&func_log10);
        math_agent_1i1o_(&func_log10, log10_para, &res);
        *log10_para = res;
    }
}
void exp_parallel_(REAL *exp_para)
{
    if ( _MYID == 0 ) 
    {
        long long int func_exp;
        REAL res;
        math_agent_exp_c_(&func_exp);
        math_agent_1i1o_(&func_exp, exp_para, &res);
        *exp_para = res;
    }
}
void pow_parallel_(REAL *pow_para)
{
    if ( _MYID == 0 ) 
    {
        long long int func_pow;
        REAL in1, in2, res;
        in1 = pow_para[0];
        in2 = pow_para[1];
        math_agent_pow_c_(&func_pow);
        math_agent_2i1o_(&func_pow, &in1, &in2, &res);
        pow_para[0] = res;
    }
}
void gamma_parallel_(REAL* x)
{
    if ( _MYID == 0 )
    {
        REAL in = *x;
        REAL res;
        //res = shr_spfn_gamma_nonintrinsic_r8_(&in);
        shr_spfn_gamma_nonintrinsic_r8_v2_(&in, &res);
        *x = res;
    }
}

#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_BCASTR(var) REG_PUTR(var, 8)
#define REG_BCASTC(var) REG_PUTC(var, 8)
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define ROWSYN  athread_syn(ROW_SCOPE,0xff)
#define COLSYN  athread_syn(COL_SCOPE,0xff)
#define ALLSYN  athread_syn(ARRAY_SCOPE,0xffff)
void ltrue_syn_(int *ltrue_, int *myid_, int *pcols_)
{
    int ltrue = *ltrue_;
    int myid = *myid_;
    int pcols = *pcols_;
    int get_ltrue;
    if ( myid >= 32 ) 
    {
        REG_PUTR(ltrue, (myid-32)/8);
        REG_GETR(get_ltrue);
    }
    else 
    {
        REG_PUTR(ltrue, (myid+32)/8);
        REG_GETR(get_ltrue);
    }
    if ( ltrue || get_ltrue ) *ltrue_ = 1;
}
void get_k_8_(REAL *array)
{
    int i = 0;
    doublev4 *a = (double*)array;
    for ( i = 0; i < 8; i ++ )
    {
        REG_GETR(a[i]);
    }
}
void get_k_3_(REAL *array)
{
    int i = 0;
    doublev4 *a = (double*)array;
    for ( i = 0; i < 3; i ++ )
    {
        REG_GETR(a[i]);
    }
}
void put_k_8_(REAL *array)
{
    int i = 0;
    doublev4 *a = (double*)array;
    int dst = COL(_MYID) + 4;
    for ( i = 0; i < 8; i ++ )
    {
        REG_PUTR(a[i], dst);
    }
}
void put_k_3_(REAL *array)
{
    int i = 0;
    doublev4 *a = (double*)array;
    int dst = COL(_MYID) + 4;
    for ( i = 0; i < 3; i ++ )
    {
        REG_PUTR(a[i], dst);
    }
}
void get_k_1_(REAL *array, int *n_, int *myid_)
{
    int i = 0;
    int n = *n_;
    int myid = *myid_;
    if ( myid >= 32 )
    {
        for ( i = 0; i < n; i ++ )
        {
            REG_GETR(array[i]);
        }
    }
}
void put_k_1_(REAL *array, int *n_, int *myid_)
{
    int i = 0;
    int n = *n_;
    int myid = *myid_;
    if ( myid < 32 )
    {
        for ( i = 0; i < n; i ++ )
        {
            REG_PUTR(array[i]  , (myid+32)/8);
        }
    }
}
void para_syn_dimk_(REAL *para, int *myid_)
{
    int i;
    int myid = *myid_;
    double tmp;
    if ( myid >= 32 ) REG_PUTR(*para, (myid-32)/8);
    else { REG_GETR(tmp);*para = tmp; }
}
void para_syn_max_dimk_(int *para, int *myid_)
{
    int i;
    int myid = *myid_;
    int tmp;
    if ( myid >= 32 ) 
    {
        REG_PUTR(*para, (myid-32)/8);
        REG_GETR(*para);
    }
    else 
    { 
        REG_GETR(tmp);
        *para = *para > tmp ? *para : tmp;
        REG_PUTR(*para, (myid+32)/8);
    }
}
void micro_mg1_0_parallel_(void *micro_mg1_0_para)
{
    //lwpf_enter(TEST);
//    if ( _MYID == 0 ) cpe_printf("myid:%d\n", _MYID);

    //set_func_ptr_();
    micro_mg1_0_args_cc micro_mg1_0_spara;
    
    dma_init();

    //micro_mg1_0_spara = *((micro_mg1_0_args_cc*)micro_mg1_0_para); 

    pe_get(micro_mg1_0_para, &micro_mg1_0_spara, sizeof(micro_mg1_0_args_cc));
    dma_syn();

    double sincos_local_data[sincos_data_len] ;
    double exp_local_data[exp_data_len] ;
    //double asin_local_data[asin_data_len] ;
    //double acos_local_data[acos_data_len] ;
    //double log_local_data[log_data_len] ;
    double log10_local_data[log10_data_len] ;
    double pow_local_data[pow_data_len] ; 

    cos_data_local_ptr = sincos_local_data ;
    sin_data_local_ptr = sincos_local_data ; 
    exp_data_local_ptr = exp_local_data ;
    //asin_data_local_ptr = asin_local_data ;
    //acos_data_local_ptr = acos_local_data ;
    //log_data_local_ptr = log_local_data ; 
    log10_data_local_ptr = log10_local_data ;  
    pow_data_local_ptr = pow_local_data ; 

    if (_MYID == 0) { 
        bcast_get(sincos_data , sincos_local_data , sincos_data_bytes) ; 
        bcast_get(exp_data , exp_local_data , exp_data_bytes) ;
        //bcast_get(asin_data , asin_local_data , asin_data_bytes) ;
        //bcast_get(acos_data , acos_local_data , acos_data_bytes) ; 
        //bcast_get(log_data , log_local_data , log_data_bytes) ;
        bcast_get(log10_data , log10_local_data , log10_data_bytes) ;
        bcast_get(pow_data , pow_local_data , pow_data_bytes) ;
        dma_syn() ; 
    } 
    athread_syn(ARRAY_SCOPE ,0xffff) ;
    // in
    REAL tn               [LAYER];
    REAL qn               [LAYER];
    REAL relvar           [LAYER];
    REAL accre_enhan      [LAYER];
    REAL p                [LAYER];
    REAL pdel             [LAYER];
    REAL cldn             [LAYER];
    REAL icecldf          [LAYER];
    REAL liqcldf          [LAYER];
    REAL naai             [LAYER];
    REAL npccnin          [LAYER];
    REAL rndst            [4][LAYER];
    REAL nacon            [4][LAYER];
    REAL tnd_qsnow        [LAYER];         
    REAL tnd_nsnow        [LAYER];
    REAL re_ice           [LAYER];
    REAL frzimm           [LAYER];
    REAL frzcnt           [LAYER];
    REAL frzdep           [LAYER];

    // out
    REAL rate1ord_cw2pr_st[LAYER];
    REAL tlat             [LAYER];
    REAL qvlat            [LAYER];
    REAL qctend           [LAYER];
    REAL qitend           [LAYER];
    REAL nctend           [LAYER];
    REAL nitend           [LAYER];
    REAL effc             [LAYER];
    REAL effc_fn          [LAYER];
    REAL effi             [LAYER];
    REAL prect            ;
    REAL preci            ;
    REAL nevapr           [LAYER];
    REAL evapsnow         [LAYER];
    REAL am_evp_st        [LAYER];
    REAL prain            [LAYER];
    REAL prodsnow         [LAYER];
    REAL cmeout           [LAYER];
    REAL deffi            [LAYER];
    REAL pgamrad          [LAYER];
    REAL lamcrad          [LAYER];
    REAL qsout            [LAYER];
    REAL dsout            [LAYER];
    REAL rflx             [LAYER+1];
    REAL sflx             [LAYER+1];
    REAL qrout            [LAYER];
    REAL qcsevap          [LAYER];
    REAL qisevap          [LAYER];
    REAL qvres            [LAYER];
    REAL cmeiout          [LAYER];
    REAL vtrmc            [LAYER];
    REAL vtrmi            [LAYER];
    REAL qcsedten         [LAYER];
    REAL qisedten         [LAYER];
    REAL prao             [LAYER];
    REAL prco             [LAYER];
    REAL mnuccco          [LAYER];
    REAL mnuccto          [LAYER];
    REAL msacwio          [LAYER];
    REAL psacwso          [LAYER];
    REAL bergso           [LAYER];
    REAL bergo            [LAYER];
    REAL melto            [LAYER];
    REAL homoo            [LAYER];
    REAL qcreso           [LAYER];
    REAL prcio            [LAYER];
    REAL praio            [LAYER];
    REAL qireso           [LAYER];
    REAL mnuccro          [LAYER];
    REAL pracso           [LAYER];
    REAL meltsdt          [LAYER];
    REAL frzrdt           [LAYER];
    REAL mnuccdo          [LAYER];
    REAL nrout            [LAYER];
    REAL nsout            [LAYER];
    REAL refl             [LAYER];
    REAL arefl            [LAYER];
    REAL areflz           [LAYER];
    REAL frefl            [LAYER];
    REAL csrfl            [LAYER];
    REAL acsrfl           [LAYER];
    REAL fcsrfl           [LAYER];
    REAL rercld           [LAYER];
    REAL ncai             [LAYER];
    REAL ncal             [LAYER];
    REAL qrout2           [LAYER];
    REAL qsout2           [LAYER];
    REAL nrout2           [LAYER];
    REAL nsout2           [LAYER];
    REAL drout2           [LAYER];
    REAL dsout2           [LAYER];
    REAL freqs            [LAYER];
    REAL freqr            [LAYER];
    REAL nfice            [LAYER];
    REAL prer_evap        [LAYER];

    // inout
    REAL qc               [LAYER];
    REAL qi               [LAYER];
    REAL nc               [LAYER];
    REAL ni               [LAYER];
    REAL reff_rain        [LAYER];
    REAL reff_snow        [LAYER];


    REAL   svp_tmelt           = micro_mg1_0_spara.svp_tmelt          ;     
    REAL   svp_h2otrip         = micro_mg1_0_spara.svp_h2otrip        ;   
    REAL   svp_tboil           = micro_mg1_0_spara.svp_tboil          ;   
    REAL   svp_ttrice          = micro_mg1_0_spara.svp_ttrice         ;  
    REAL   svp_epsilo          = micro_mg1_0_spara.svp_epsilo         ;  
    REAL   g                   = micro_mg1_0_spara.g                  ;     
    REAL   r                   = micro_mg1_0_spara.r                  ;     
    REAL   rv                  = micro_mg1_0_spara.rv                 ;     
    REAL   cpp                 = micro_mg1_0_spara.cpp                ;     
    REAL   rhow                = micro_mg1_0_spara.rhow               ;     
    REAL   tmelt               = micro_mg1_0_spara.tmelt              ;     
    REAL   xxlv                = micro_mg1_0_spara.xxlv               ;     
    REAL   xlf                 = micro_mg1_0_spara.xlf                ;     
    REAL   xxls                = micro_mg1_0_spara.xxls               ;     
    REAL   rhosn               = micro_mg1_0_spara.rhosn              ;     
    REAL   rhoi                = micro_mg1_0_spara.rhoi               ;     
    REAL   ac                  = micro_mg1_0_spara.ac                 ;     
    REAL   bc                  = micro_mg1_0_spara.bc                 ;     
    REAL   as                  = micro_mg1_0_spara.as                 ;     
    REAL   bs                  = micro_mg1_0_spara.bs                 ;     
    REAL   ai                  = micro_mg1_0_spara.ai                 ;     
    REAL   bi                  = micro_mg1_0_spara.bi                 ;     
    REAL   ar                  = micro_mg1_0_spara.ar                 ;     
    REAL   br                  = micro_mg1_0_spara.br                 ;     
    REAL   ci                  = micro_mg1_0_spara.ci                 ;     
    REAL   di                  = micro_mg1_0_spara.di                 ;     
    REAL   cs                  = micro_mg1_0_spara.cs                 ;     
    REAL   ds                  = micro_mg1_0_spara.ds                 ;     
    REAL   cr                  = micro_mg1_0_spara.cr                 ;     
    REAL   dr                  = micro_mg1_0_spara.dr                 ;     
    REAL   f1s                 = micro_mg1_0_spara.f1s                ;     
    REAL   f2s                 = micro_mg1_0_spara.f2s                ;     
    REAL   Eii                 = micro_mg1_0_spara.Eii                ;     
    REAL   Ecr                 = micro_mg1_0_spara.Ecr                ;     
    REAL   f1r                 = micro_mg1_0_spara.f1r                ;     
    REAL   f2r                 = micro_mg1_0_spara.f2r                ;     
    REAL   DCS                 = micro_mg1_0_spara.DCS                ;     
    REAL   qsmall              = micro_mg1_0_spara.qsmall             ;     
    REAL   bimm                = micro_mg1_0_spara.bimm               ;     
    REAL   aimm                = micro_mg1_0_spara.aimm               ;     
    REAL   rhosu               = micro_mg1_0_spara.rhosu              ;     
    REAL   mi0                 = micro_mg1_0_spara.mi0                ;     
    REAL   rin                 = micro_mg1_0_spara.rin                ;     
    REAL   pi                  = micro_mg1_0_spara.pi                 ;     
    REAL   cons1               = micro_mg1_0_spara.cons1              ;     
    REAL   cons4               = micro_mg1_0_spara.cons4              ;     
    REAL   cons5               = micro_mg1_0_spara.cons5              ;     
    REAL   cons6               = micro_mg1_0_spara.cons6              ;     
    REAL   cons7               = micro_mg1_0_spara.cons7              ;     
    REAL   cons8               = micro_mg1_0_spara.cons8              ;     
    REAL   cons11              = micro_mg1_0_spara.cons11             ;     
    REAL   cons13              = micro_mg1_0_spara.cons13             ;     
    REAL   cons14              = micro_mg1_0_spara.cons14             ;     
    REAL   cons16              = micro_mg1_0_spara.cons16             ;     
    REAL   cons17              = micro_mg1_0_spara.cons17             ;     
    REAL   cons22              = micro_mg1_0_spara.cons22             ;     
    REAL   cons23              = micro_mg1_0_spara.cons23             ;     
    REAL   cons24              = micro_mg1_0_spara.cons24             ;     
    REAL   cons25              = micro_mg1_0_spara.cons25             ;     
    REAL   cons27              = micro_mg1_0_spara.cons27             ;     
    REAL   cons28              = micro_mg1_0_spara.cons28             ;     
    REAL   lammini             = micro_mg1_0_spara.lammini            ;     
    REAL   lammaxi             = micro_mg1_0_spara.lammaxi            ;     
    REAL   lamminr             = micro_mg1_0_spara.lamminr            ;     
    REAL   lammaxr             = micro_mg1_0_spara.lammaxr            ;     
    REAL   lammins             = micro_mg1_0_spara.lammins            ;     
    REAL   lammaxs             = micro_mg1_0_spara.lammaxs            ;     
    REAL   tmax_fsnow          = micro_mg1_0_spara.tmax_fsnow         ;     
    REAL   tmin_fsnow          = micro_mg1_0_spara.tmin_fsnow         ;     
    REAL   tt0                 = micro_mg1_0_spara.tt0                ;     
    REAL   csmin               = micro_mg1_0_spara.csmin              ;     
    REAL   csmax               = micro_mg1_0_spara.csmax              ;     
    REAL   minrefl             = micro_mg1_0_spara.minrefl            ;     
    REAL   mindbz              = micro_mg1_0_spara.mindbz             ;     
    REAL   rhmini              = micro_mg1_0_spara.rhmini             ;     
    int   use_hetfrz_classnuc = micro_mg1_0_spara.use_hetfrz_classnuc;   

    int  microp_uniform        = micro_mg1_0_spara.microp_uniform;
    int  do_cldice             = micro_mg1_0_spara.do_cldice;
    int  pcols                 = micro_mg1_0_spara.pcols;
    int  pver                  = micro_mg1_0_spara.pver;
    int  ncol                  = micro_mg1_0_spara.ncol;
    int  top_lev               = micro_mg1_0_spara.top_lev;
    REAL deltatin              = micro_mg1_0_spara.deltatin;

    int myid = COL(_MYID)*8 + ROW(_MYID);

    int len = pver / 2;
	int i;
    int idx, idy;


    //prect            = 0.0;
    //preci            = 0.0;
	//rflx[0] = 0.0;
	//sflx[0] = 0.0;
	//for ( i = 0; i < LAYER; i ++ )
	//{
	//	rate1ord_cw2pr_st[i]       = 0.0;
	//	tlat             [i]       = 0.0;
	//	qvlat            [i]       = 0.0;
	//	qctend           [i]       = 0.0;
	//	qitend           [i]       = 0.0;
	//	nctend           [i]       = 0.0;
	//	nitend           [i]       = 0.0;
	//	effc             [i]       = 0.0;
	//	effc_fn          [i]       = 0.0;
	//	effi             [i]       = 0.0;
	//	nevapr           [i]       = 0.0;
	//	evapsnow         [i]       = 0.0;
	//	am_evp_st        [i]       = 0.0;
	//	prain            [i]       = 0.0;
	//	prodsnow         [i]       = 0.0;
	//	cmeout           [i]       = 0.0;
	//	deffi            [i]       = 0.0;
	//	pgamrad          [i]       = 0.0;
	//	lamcrad          [i]       = 0.0;
	//	qsout            [i]       = 0.0;
	//	dsout            [i]       = 0.0;
	//	rflx             [i+1]     = 0.0;
	//	sflx             [i+1]     = 0.0;
	//	qrout            [i]       = 0.0;
	//	qcsevap          [i]       = 0.0;
	//	qisevap          [i]       = 0.0;
	//	qvres            [i]       = 0.0;
	//	cmeiout          [i]       = 0.0;
	//	vtrmc            [i]       = 0.0;
	//	vtrmi            [i]       = 0.0;
	//	qcsedten         [i]       = 0.0;
	//	qisedten         [i]       = 0.0;
	//	prao             [i]       = 0.0;
	//	prco             [i]       = 0.0;
	//	mnuccco          [i]       = 0.0;
	//	mnuccto          [i]       = 0.0;
	//	msacwio          [i]       = 0.0;
	//	psacwso          [i]       = 0.0;
	//	bergso           [i]       = 0.0;
	//	bergo            [i]       = 0.0;
	//	melto            [i]       = 0.0;
	//	homoo            [i]       = 0.0;
	//	qcreso           [i]       = 0.0;
	//	prcio            [i]       = 0.0;
	//	praio            [i]       = 0.0;
	//	qireso           [i]       = 0.0;
	//	mnuccro          [i]       = 0.0;
	//	pracso           [i]       = 0.0;
	//	meltsdt          [i]       = 0.0;
	//	frzrdt           [i]       = 0.0;
	//	mnuccdo          [i]       = 0.0;
	//	nrout            [i]       = 0.0;
	//	nsout            [i]       = 0.0;
	//	refl             [i]       = 0.0;
	//	arefl            [i]       = 0.0;
	//	areflz           [i]       = 0.0;
	//	frefl            [i]       = 0.0;
	//	csrfl            [i]       = 0.0;
	//	acsrfl           [i]       = 0.0;
	//	fcsrfl           [i]       = 0.0;
	//	rercld           [i]       = 0.0;
	//	ncai             [i]       = 0.0;
	//	ncal             [i]       = 0.0;
	//	qrout2           [i]       = 0.0;
	//	qsout2           [i]       = 0.0;
	//	nrout2           [i]       = 0.0;
	//	nsout2           [i]       = 0.0;
	//	drout2           [i]       = 0.0;
	//	dsout2           [i]       = 0.0;
	//	freqs            [i]       = 0.0;
	//	freqr            [i]       = 0.0;
	//	nfice            [i]       = 0.0;
	//	prer_evap        [i]       = 0.0;
	//}

	//if ( _MYID == 0 ) cpe_printf("ni address:%lld\n", micro_mg1_0_spara.ni);

	//if ( _MYID == 8 )  cpe_printf("myid:%d pcols:%d pver:%d top_lev:%d\n", myid, pcols, pver, top_lev);

#define OFFSET(dim) (idy * (dim) * sizeof(REAL) + idx * (dim/2) * sizeof(REAL))
	//lwpf_start(A);
	if ( myid%32 < pcols )
	{
		idx = myid / 32;
		idy = myid % 32;
		//lwpf_start(B);
		pe_get(micro_mg1_0_spara.tn + OFFSET(pver), tn, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.qn + OFFSET(pver), qn, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.relvar + OFFSET(pver), relvar, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.accre_enhan + OFFSET(pver), accre_enhan, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.p + OFFSET(pver), p, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.pdel + OFFSET(pver), pdel, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.cldn + OFFSET(pver), cldn, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.icecldf + OFFSET(pver), icecldf, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.liqcldf + OFFSET(pver), liqcldf, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.naai + OFFSET(pver), naai, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.npccnin + OFFSET(pver), npccnin, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.tnd_qsnow + OFFSET(pver), tnd_qsnow, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.tnd_nsnow + OFFSET(pver), tnd_nsnow, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.re_ice    + OFFSET(pver), re_ice   , len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.frzimm    + OFFSET(pver), frzimm   , len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.frzcnt    + OFFSET(pver), frzcnt   , len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.frzdep    + OFFSET(pver), frzdep   , len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.rndst + OFFSET(pver) + 0*pcols*pver*sizeof(REAL), &rndst[0][0], len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.nacon + OFFSET(pver) + 0*pcols*pver*sizeof(REAL), &nacon[0][0], len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.rndst + OFFSET(pver) + 1*pcols*pver*sizeof(REAL), &rndst[1][0], len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.nacon + OFFSET(pver) + 1*pcols*pver*sizeof(REAL), &nacon[1][0], len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.rndst + OFFSET(pver) + 2*pcols*pver*sizeof(REAL), &rndst[2][0], len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.nacon + OFFSET(pver) + 2*pcols*pver*sizeof(REAL), &nacon[2][0], len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.rndst + OFFSET(pver) + 3*pcols*pver*sizeof(REAL), &rndst[3][0], len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.nacon + OFFSET(pver) + 3*pcols*pver*sizeof(REAL), &nacon[3][0], len*sizeof(REAL));

		pe_get(micro_mg1_0_spara.qc + OFFSET(pver), qc, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.qi + OFFSET(pver), qi, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.nc + OFFSET(pver), nc, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.ni + OFFSET(pver), ni, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.reff_rain + OFFSET(pver), reff_rain, len*sizeof(REAL));
		pe_get(micro_mg1_0_spara.reff_snow + OFFSET(pver), reff_snow, len*sizeof(REAL));

		dma_syn();
		//lwpf_stop(B);

		//for ( i = 0; i < len; i ++ )
		//{
		//	ni[i] = ((REAL*)micro_mg1_0_spara.ni)[myid*pver+i];
		//}

		//for ( i = 0; i < len; i ++ )
		//{
		//    if ( abs(tn         [i] - ((REAL*)micro_mg1_0_spara.tn         )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("tn         wrong, %.20lf \n", myid);
		//    if ( abs(qn         [i] - ((REAL*)micro_mg1_0_spara.qn         )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("qn         wrong, %.20lf \n", myid);
		//    if ( abs(relvar     [i] - ((REAL*)micro_mg1_0_spara.relvar     )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("relvar     wrong, %.20lf \n", myid);
		//    if ( abs(accre_enhan[i] - ((REAL*)micro_mg1_0_spara.accre_enhan)[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("accre_enhanwrong, %.20lf \n", myid);
		//    if ( abs(p          [i] - ((REAL*)micro_mg1_0_spara.p          )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("p          wrong, %.20lf \n", myid);
		//    if ( abs(pdel       [i] - ((REAL*)micro_mg1_0_spara.pdel       )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("pdel       wrong, %.20lf \n", myid);
		//    if ( abs(cldn       [i] - ((REAL*)micro_mg1_0_spara.cldn       )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("cldn       wrong, %.20lf \n", myid);
		//    if ( abs(icecldf    [i] - ((REAL*)micro_mg1_0_spara.icecldf    )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("icecldf    wrong, %.20lf \n", myid);
		//    if ( abs(liqcldf    [i] - ((REAL*)micro_mg1_0_spara.liqcldf    )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("liqcldf    wrong, %.20lf \n", myid);
		//    if ( abs(naai       [i] - ((REAL*)micro_mg1_0_spara.naai       )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("naai       wrong, %.20lf \n", myid);
		//    if ( abs(npccnin    [i] - ((REAL*)micro_mg1_0_spara.npccnin    )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("npccnin    wrong, %.20lf \n", myid);
		//    if ( abs(rndst   [0][i] - ((REAL*)micro_mg1_0_spara.rndst      )[idy*pver+i+idx*pver/2+0*pcols*pver]) > 1e-15 ) cpe_printf("rndst 0    wrong, %.20lf \n", myid);
		//    if ( abs(nacon   [0][i] - ((REAL*)micro_mg1_0_spara.nacon      )[idy*pver+i+idx*pver/2+0*pcols*pver]) > 1e-15 ) cpe_printf("nacon 0    wrong, %.20lf \n", myid);
		//    if ( abs(rndst   [1][i] - ((REAL*)micro_mg1_0_spara.rndst      )[idy*pver+i+idx*pver/2+1*pcols*pver]) > 1e-15 ) cpe_printf("rndst 1    wrong, %.20lf \n", myid);
		//    if ( abs(nacon   [1][i] - ((REAL*)micro_mg1_0_spara.nacon      )[idy*pver+i+idx*pver/2+1*pcols*pver]) > 1e-15 ) cpe_printf("nacon 1    wrong, %.20lf \n", myid);
		//    if ( abs(rndst   [2][i] - ((REAL*)micro_mg1_0_spara.rndst      )[idy*pver+i+idx*pver/2+2*pcols*pver]) > 1e-15 ) cpe_printf("rndst 2    wrong, %.20lf \n", myid);
		//    if ( abs(nacon   [2][i] - ((REAL*)micro_mg1_0_spara.nacon      )[idy*pver+i+idx*pver/2+2*pcols*pver]) > 1e-15 ) cpe_printf("nacon 2    wrong, %.20lf \n", myid);
		//    if ( abs(rndst   [3][i] - ((REAL*)micro_mg1_0_spara.rndst      )[idy*pver+i+idx*pver/2+3*pcols*pver]) > 1e-15 ) cpe_printf("rndst 3    wrong, %.20lf \n", myid);
		//    if ( abs(nacon   [3][i] - ((REAL*)micro_mg1_0_spara.nacon      )[idy*pver+i+idx*pver/2+3*pcols*pver]) > 1e-15 ) cpe_printf("nacon 3    wrong, %.20lf \n", myid);
		//    if ( abs(qc         [i] - ((REAL*)micro_mg1_0_spara.qc         )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("qc         wrong, %.20lf \n", myid);
		//    if ( abs(qi         [i] - ((REAL*)micro_mg1_0_spara.qi         )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("qi         wrong, %.20lf \n", myid);
		//    if ( abs(nc         [i] - ((REAL*)micro_mg1_0_spara.nc         )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("nc         wrong, %.20lf \n", myid);
		//    if ( abs(ni         [i] - ((REAL*)micro_mg1_0_spara.ni         )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("ni         wrong, %.20lf \n", myid);
		//    if ( abs(reff_rain  [i] - ((REAL*)micro_mg1_0_spara.reff_rain  )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("reff_rain  wrong, %.20lf \n", myid);
		//    if ( abs(reff_snow  [i] - ((REAL*)micro_mg1_0_spara.reff_snow  )[idy*pver+i+idx*pver/2])              > 1e-15 ) cpe_printf("reff_snow  wrong, %.20lf \n", myid);
		//}

		//for ( i = 0; i < len; i ++ )
		//{
		//	qrout2[i] = nacon[0][i];
		//	qsout2[i] = nacon[1][i];
		//	drout2[i] = nacon[2][i];
		//	dsout2[i] = nacon[3][i];
		//}
		//test_(nacon_t, qrout2, qsout2, drout2, dsout2, len)
		//lwpf_start(C);
		micro_mg1_0_compute_( \
				&g, &r, &rv, &cpp, &rhow, &tmelt, &xxlv, &xlf, \
				&xxls, &rhosn, &rhoi, &ac, &bc, &as, &bs, &ai, &bi, &ar, \
				&br, &ci, &di, &cs, &ds, &cr, &dr, &f1s, &f2s, &Eii, &Ecr, \
				&f1r, &f2r, &DCS, &qsmall, &bimm, &aimm, &rhosu, &mi0, \
				&rin, &pi, &cons1, &cons4, &cons5, &cons6, &cons7, \   
				&cons8, &cons11, &cons13, &cons14, &cons16, &cons17, \
				&cons22, &cons23, &cons24, &cons25, &cons27, &cons28, \
				&lammini, &lammaxi, &lamminr, &lammaxr, \
				&lammins, &lammaxs, &tmax_fsnow, &tmin_fsnow, &tt0, &csmin, \
				&csmax, &minrefl, &mindbz, &rhmini, &use_hetfrz_classnuc, \
				&microp_uniform, &pcols, &len, &ncol, &top_lev, &deltatin,\
				tn, qn, qc, qi, nc,                                        \
				ni, p, pdel, cldn, liqcldf,                                \
				relvar, accre_enhan,                                       \
				icecldf, rate1ord_cw2pr_st, naai, npccnin,                 \
				rndst, nacon, tlat, qvlat, qctend,                         \
				qitend, nctend, nitend, effc, effc_fn,                     \
				effi, &prect, &preci, nevapr, evapsnow, am_evp_st,         \
				prain, prodsnow, cmeout, deffi, pgamrad,                   \
				lamcrad, qsout, dsout, rflx, sflx,                         \
				qrout, reff_rain, reff_snow, qcsevap, qisevap,             \
				qvres, cmeiout, vtrmc, vtrmi, qcsedten,                    \
				qisedten, prao, prco, mnuccco, mnuccto,                    \
				msacwio, psacwso, bergso, bergo, melto,                    \
				homoo, qcreso, prcio, praio, qireso,                       \
				mnuccro, pracso, meltsdt, frzrdt, mnuccdo,                 \
				nrout, nsout, refl, arefl, areflz,                         \
				frefl, csrfl, acsrfl, fcsrfl, rercld,                      \
				ncai, ncal, qrout2, qsout2, nrout2,                        \
				nsout2, drout2, dsout2, freqs, freqr,                      \
				nfice, prer_evap, &do_cldice,                              \
				tnd_qsnow, tnd_nsnow, re_ice,                              \
				frzimm, frzcnt, frzdep,                                    \
				&svp_tmelt, &svp_h2otrip, &svp_tboil,                      \
				&svp_ttrice, &svp_epsilo, &myid                      );
		//lwpf_stop(C);

		//for ( i = 0; i < len; i ++ )
		//{
		//	qrout2[i] = nacon[0][i];
		//	qsout2[i] = nacon[1][i];
		//	drout2[i] = nacon[2][i];
		//	dsout2[i] = nacon[3][i];
		////    qrout2[i] = ((REAL*)(micro_mg1_0_spara.ni))[myid*pver+i];
		//}

		//lwpf_start(A);
		pe_put(micro_mg1_0_spara.rate1ord_cw2pr_st + OFFSET(pver), rate1ord_cw2pr_st, len*sizeof(REAL));

		pe_put(micro_mg1_0_spara.tlat        + OFFSET(pver), tlat       , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qvlat       + OFFSET(pver), qvlat      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qctend      + OFFSET(pver), qctend     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qitend      + OFFSET(pver), qitend     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nctend      + OFFSET(pver), nctend     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nitend      + OFFSET(pver), nitend     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.effc        + OFFSET(pver), effc       , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.effc_fn     + OFFSET(pver), effc_fn    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.effi        + OFFSET(pver), effi       , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nevapr      + OFFSET(pver), nevapr     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.evapsnow    + OFFSET(pver), evapsnow   , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.am_evp_st   + OFFSET(pver), am_evp_st  , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.prain       + OFFSET(pver), prain      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.prodsnow    + OFFSET(pver), prodsnow   , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.cmeout      + OFFSET(pver), cmeout     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.deffi       + OFFSET(pver), deffi      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.pgamrad     + OFFSET(pver), pgamrad    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.lamcrad     + OFFSET(pver), lamcrad    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qsout       + OFFSET(pver), qsout      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.dsout       + OFFSET(pver), dsout      , len*sizeof(REAL));
#define OFFSET1(dim) (idy * (dim) * sizeof(REAL) + idx * ((dim)/2+1) * sizeof(REAL))
		if ( myid < 32 )
		{
			pe_put(micro_mg1_0_spara.rflx        + OFFSET1(pver+1), rflx       , (1+len)*sizeof(REAL));
			pe_put(micro_mg1_0_spara.sflx        + OFFSET1(pver+1), sflx       , (1+len)*sizeof(REAL));
		}
		else
		{
			pe_put(micro_mg1_0_spara.rflx        + OFFSET1(pver+1), rflx+1       , (len)*sizeof(REAL));
			pe_put(micro_mg1_0_spara.sflx        + OFFSET1(pver+1), sflx+1       , (len)*sizeof(REAL));
			pe_put(micro_mg1_0_spara.prect       + idy*sizeof(REAL)   , &prect     , 1*sizeof(REAL));
			pe_put(micro_mg1_0_spara.preci       + idy*sizeof(REAL)   , &preci     , 1*sizeof(REAL));
		}
		pe_put(micro_mg1_0_spara.qrout       + OFFSET(pver), qrout      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qcsevap     + OFFSET(pver), qcsevap    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qisevap     + OFFSET(pver), qisevap    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qvres       + OFFSET(pver), qvres      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.cmeiout     + OFFSET(pver), cmeiout    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.vtrmc       + OFFSET(pver), vtrmc      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.vtrmi       + OFFSET(pver), vtrmi      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qcsedten    + OFFSET(pver), qcsedten   , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qisedten    + OFFSET(pver), qisedten   , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.prao        + OFFSET(pver), prao       , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.prco        + OFFSET(pver), prco       , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.mnuccco     + OFFSET(pver), mnuccco    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.mnuccto     + OFFSET(pver), mnuccto    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.msacwio     + OFFSET(pver), msacwio    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.psacwso     + OFFSET(pver), psacwso    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.bergso      + OFFSET(pver), bergso     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.bergo       + OFFSET(pver), bergo      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.melto       + OFFSET(pver), melto      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.homoo       + OFFSET(pver), homoo      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qcreso      + OFFSET(pver), qcreso     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.prcio       + OFFSET(pver), prcio      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.praio       + OFFSET(pver), praio      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qireso      + OFFSET(pver), qireso     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.mnuccro     + OFFSET(pver), mnuccro    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.pracso      + OFFSET(pver), pracso     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.meltsdt     + OFFSET(pver), meltsdt    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.frzrdt      + OFFSET(pver), frzrdt     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.mnuccdo     + OFFSET(pver), mnuccdo    , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nrout       + OFFSET(pver), nrout      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nsout       + OFFSET(pver), nsout      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.refl        + OFFSET(pver), refl       , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.arefl       + OFFSET(pver), arefl      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.areflz      + OFFSET(pver), areflz     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.frefl       + OFFSET(pver), frefl      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.csrfl       + OFFSET(pver), csrfl      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.acsrfl      + OFFSET(pver), acsrfl     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.fcsrfl      + OFFSET(pver), fcsrfl     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.rercld      + OFFSET(pver), rercld     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.ncai        + OFFSET(pver), ncai       , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.ncal        + OFFSET(pver), ncal       , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qrout2      + OFFSET(pver), qrout2     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qsout2      + OFFSET(pver), qsout2     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nrout2      + OFFSET(pver), nrout2     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nsout2      + OFFSET(pver), nsout2     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.drout2      + OFFSET(pver), drout2     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.dsout2      + OFFSET(pver), dsout2     , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.freqs       + OFFSET(pver), freqs      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.freqr       + OFFSET(pver), freqr      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nfice       + OFFSET(pver), nfice      , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.prer_evap   + OFFSET(pver), prer_evap  , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qc          + OFFSET(pver), qc         , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.qi          + OFFSET(pver), qi         , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.nc          + OFFSET(pver), nc         , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.ni          + OFFSET(pver), ni         , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.reff_rain   + OFFSET(pver), reff_rain  , len*sizeof(REAL));
		pe_put(micro_mg1_0_spara.reff_snow   + OFFSET(pver), reff_snow  , len*sizeof(REAL));
		dma_syn();
		//lwpf_stop(A);
	}                              
	cos_data_local_ptr = sincos_data ;
	sin_data_local_ptr = sincos_data ; 
	exp_data_local_ptr = exp_data ;
	log_data_local_ptr = log_data ; 
	log10_data_local_ptr = log10_data ;  
	pow_data_local_ptr = pow_data ; 

	//long long sp;
	//get_sp_(&sp);
	//if ( _MYID == 0 ) printf("sp:%lld\n", sp);

	//lwpf_stop(A);
	//lwpf_exit(TEST);
}                                  
