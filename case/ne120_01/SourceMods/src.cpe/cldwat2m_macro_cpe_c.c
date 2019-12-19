/*************************************************************************
  > File Name: src.cpe/cldwat2m_macro_cpe_c.c
  > Author: Xu Kai
  > Created Time: 2018年12月12日 星期三 13时58分32秒
 ************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include "dma_macros.h"
#include <stdarg.h>
#include "cpe_print.h"
//#include "lwpf_cpe_c.h"
#include"math_data.h"
//#include<regdef.h>
#define REAL double
typedef struct mmacro_pcond_args{

    void* T0;               
    void* qv0;              
    void* ql0;              
    void* qi0;              
    void* nl0;              
    void* ni0;              

    void* p;                
    void* dp;               

    void* A_T;              
    void* A_qv;             
    void* A_ql;             
    void* A_qi;             
    void* A_nl;             
    void* A_ni;             

    void* C_T;              
    void* C_qv;             
    void* C_ql;             
    void* C_qi;             
    void* C_nl;             
    void* C_ni;             
    void* C_qlst;           


    void* D_T;              
    void* D_qv;             
    void* D_ql;             
    void* D_qi;             
    void* D_nl;             
    void* D_ni;             

    void* a_cud;            
    void* a_cu0;            

    void* clrw_old;         
    void* clri_old;         

    void* tke;                     
    void* qtl_flx;                 
    void* qti_flx;                 
    void* cmfr_det;                
    void* qlr_det;                 
    void* qir_det;                 


    void* s_tendout;        
    void* qv_tendout;       
    void* ql_tendout;       
    void* qi_tendout;       
    void* nl_tendout;       
    void* ni_tendout;       

    void* qme;              
    void* qvadj;            
    void* qladj;            
    void* qiadj;            
    void* qllim;            
    void* qilim;            

    void* cld;              
    void* al_st_star;       
    void* ai_st_star;       
    void* ql_st_star;       
    void* qi_st_star;       

    void* d_rhmin_liq_PBL;
    void* d_rhmin_ice_PBL;
    void* d_rhmin_liq_det;
    void* d_rhmin_ice_det;

    void* landfrac;              
    void* snowh;                 
    void* wv_para;
    void* endrunflag;

    REAL dt;

    REAL rhmini_const;
    REAL rhmaxi_const;
    REAL rhminl_const;
    REAL rhminl_adj_land_const;
    REAL rhminh_const;
    REAL qmin1;
    REAL qmin2;
    REAL qmin3;
    REAL premib;
    REAL premit;
    REAL icecrit;

    //long long int spa1;
    //long long int spa2;

    int pcols;
    int pver;
    int pverp;
    int top_lev;
    int ncol;
    int do_cldice_int;
    int i_rhmini;
    int i_rhminl;
    int iceopt;


} mmacro_pcond_args_c;

#define LAYER 16
#define PARASIZE 12
extern void printf_double_c(double); 
static inline void setftz(){
    long fpcr;
    asm volatile("rfpcr %0": "=r"(fpcr));
    fpcr |= 1L << 60 | 1L << 2;
    asm volatile("wfpcr %0": : "r"(fpcr));
}

static inline void unsetftz(){
    long fpcr;
    asm volatile("rfpcr %0": "=r"(fpcr));
    fpcr |= 1L << 2;
    fpcr &= ~(1L << 60);
    asm volatile("wfpcr %0": : "r"(fpcr));
}
void mmacro_pcond_cpe_c_(void *args){

    //setftz();

    dma_init(); 
    double log_local_data[log_data_len] ;
    double log10_local_data[log10_data_len] ;
    double pow_local_data[pow_data_len] ; 
    double sincos_local_data[sincos_data_len] ;
    double acos_local_data[acos_data_len] ;

    log_data_local_ptr = log_local_data ; 
    log10_data_local_ptr = log10_local_data ;  
    pow_data_local_ptr = pow_local_data ; 
    cos_data_local_ptr = sincos_local_data ; 
    acos_data_local_ptr = acos_local_data ; 

    if (_MYID == 0) { 
        bcast_get(log_data , log_local_data , log_data_bytes) ;
        bcast_get(log10_data , log10_local_data , log10_data_bytes) ;
        bcast_get(pow_data , pow_local_data , pow_data_bytes) ;
        bcast_get(sincos_data , sincos_local_data , sincos_data_bytes) ;
        bcast_get(acos_data , acos_local_data , acos_data_bytes) ;
        dma_syn() ; 
    } 
    athread_syn(ARRAY_SCOPE ,0xffff) ;



    mmacro_pcond_args_c spara; 
    //spara.qmin1 = 0.0;

    pe_get(args, &spara, sizeof(mmacro_pcond_args_c)); 
    //memcpy(&spara, args, sizeof(mmacro_pcond_args_c)); 

    //volatile mmacro_pcond_args_c *sparaa = args; 

    dma_syn();
    //double qminslave = spara.qmin1; 

    //if(_MYID == 0){
    //    //printf_double_c(1.0);
    //    if(qminslave == 0)
    //        printf_double_c(qminslave);
    //}
#define GET_COL(col) asm volatile("rcsr    %0, 2" : "=r"(col))
#define GET_ROW(row) asm volatile("rcsr    %0, 1" : "=r"(row))
    int col_id;
    int row_id;
    GET_COL(col_id);
    GET_ROW(row_id);
    int core_id = col_id * 8 + row_id; 
    // inout
    REAL T0[LAYER];               
    REAL qv0[LAYER];              
    REAL ql0[LAYER];              
    REAL qi0[LAYER];              
    REAL nl0[LAYER];              
    REAL ni0[LAYER];              

    //in
    REAL p[LAYER];                
    REAL dp[LAYER];               

    REAL A_T[LAYER];              
    REAL A_qv[LAYER];             
    REAL A_ql[LAYER];             
    REAL A_qi[LAYER];             
    REAL A_nl[LAYER];             
    REAL A_ni[LAYER];             

    REAL C_T[LAYER];              
    REAL C_qv[LAYER];             
    REAL C_ql[LAYER];             
    REAL C_qi[LAYER];             
    REAL C_nl[LAYER];             
    REAL C_ni[LAYER];             
    REAL C_qlst[LAYER];           


    REAL D_T[LAYER];              
    REAL D_qv[LAYER];             
    REAL D_ql[LAYER];             
    REAL D_qi[LAYER];             
    REAL D_nl[LAYER];             
    REAL D_ni[LAYER];             

    REAL a_cud[LAYER];            
    REAL a_cu0[LAYER];            

    REAL clrw_old[LAYER];         
    REAL clri_old[LAYER];         

    REAL tke[LAYER];                     
    REAL qtl_flx[LAYER];                 
    REAL qti_flx[LAYER];                 
    REAL cmfr_det[LAYER];                
    REAL qlr_det[LAYER];                 
    REAL qir_det[LAYER];                 


    REAL s_tendout[LAYER];        
    REAL qv_tendout[LAYER];       
    REAL ql_tendout[LAYER];       
    REAL qi_tendout[LAYER];       
    REAL nl_tendout[LAYER];       
    REAL ni_tendout[LAYER];       

    REAL qme[LAYER];              
    REAL qvadj[LAYER];            
    REAL qladj[LAYER];            
    REAL qiadj[LAYER];            
    REAL qllim[LAYER];            
    REAL qilim[LAYER];            

    REAL cld[LAYER];              
    REAL al_st_star[LAYER];       
    REAL ai_st_star[LAYER];       
    REAL ql_st_star[LAYER];       
    REAL qi_st_star[LAYER];       

    REAL d_rhmin_liq_PBL[LAYER];
    REAL d_rhmin_ice_PBL[LAYER];
    REAL d_rhmin_liq_det[LAYER];
    REAL d_rhmin_ice_det[LAYER];

    REAL landfrac;
    REAL snowh;
    REAL wv_para[PARASIZE];
    //REAL aa[4];
    //REAL bb[4];

    REAL dt;
    REAL rhmini_const;
    REAL rhmaxi_const;
    REAL rhminl_const;
    REAL rhminl_adj_land_const;
    REAL rhminh_const;
    REAL qmin1;
    REAL qmin2;
    REAL qmin3;
    REAL premib;
    REAL premit;
    REAL icecrit;

    int pcols;
    int pver;
    int pverp;
    int top_lev;
    int ncol;
    int do_cldice_int;
    int i_rhmini;
    int i_rhminl;
    int iceopt;
    int endrunflag;
    long long int spa2;

    //int dim1 = 2;
    //int dim2 = 1;


    dt    = spara.dt;
    rhmini_const = spara.rhmini_const;
    rhmaxi_const = spara.rhmaxi_const;
    rhminl_const = spara.rhminl_const;
    rhminl_adj_land_const = spara.rhminl_adj_land_const;
    rhminh_const = spara.rhminh_const;
    qmin1       = spara.qmin1;
    qmin2       = spara.qmin2;
    qmin3       = spara.qmin3;
    premib = spara.premib;
    premit = spara.premit;
    icecrit = spara.icecrit;

    pcols = spara.pcols;
    pver  = spara.pver;
    pverp = spara.pverp;
    top_lev = spara.top_lev;
    ncol  = spara.ncol;
    do_cldice_int = spara.do_cldice_int;   
    i_rhmini = spara.i_rhmini; 
    i_rhminl = spara.i_rhminl;
    iceopt = spara.iceopt;
    int dimsize = pver;
    int dimnum = 0;
    int offsetnum = 0;
    int pver_org = pver;
    int pverp_org = pverp;


#define OFFSETARR(dim) (((core_id / 2) * dim + (core_id%2)*(dim / 2)) * sizeof(REAL)) 

#define OFFSET(dim) ((core_id/2) * (dim) * sizeof(REAL)) 

#define getmem(arr, num) {\
    offsetnum = OFFSETARR(num);\
    dimnum = (((core_id%2) == 0 )? (num / 2):(num - num/2)); \
    pe_get(spara.arr + offsetnum, arr,  dimnum*sizeof(REAL));  \
}

#define getmemp(arr, num) {\
    offsetnum = OFFSETARR(num);\
    dimnum = (((core_id%2) == 0) ? ((num / 2) + 1):(num - num/2));\
    pe_get(spara.arr + offsetnum, arr,  dimnum*sizeof(REAL));\
}

#define getpara(arr, num) \
        pe_get(spara.arr, arr, (num)*sizeof(REAL))

#define putpara(arr, num) \
        pe_put(spara.arr, arr, (num)*sizeof(REAL))

#define getreal(arr) \
        pe_get(spara.arr + OFFSET(1), &arr,  (1)*sizeof(REAL))

#define putmem(arr, num) {\
    offsetnum = OFFSETARR(num); \
    dimnum = (((core_id%2) == 0 )? (num / 2):(num - num/2)); \
    pe_put(spara.arr + offsetnum, arr,  dimnum*sizeof(REAL)); \
}

#define putint(arr) \
        pe_put(spara.arr + core_id * sizeof(int) , &arr,  1 * sizeof(int))

    if(core_id < ncol*2){
        // inout
        getmem(T0,  pver);               
        getmem(qv0, pver);              
        getmem(ql0, pver);              
        getmem(qi0, pver);              
        getmem(nl0, pver);              
        getmem(ni0, pver);              

        //in
        getmem(p,   pver);                
        getmem(dp,  pver);               

        getmem(A_T, pver);              
        getmem(A_qv,pver);             
        getmem(A_ql,pver);             
        getmem(A_qi,pver);             
        getmem(A_nl,pver);             
        getmem(A_ni,pver);             

        getmem(C_T,pver);              
        getmem(C_qv,pver);             
        getmem(C_ql,pver);             
        getmem(C_qi,pver);             
        getmem(C_nl,pver);             
        getmem(C_ni,pver);             
        getmem(C_qlst,pver);           


        getmem(D_T,pver);              
        getmem(D_qv,pver);             
        getmem(D_ql,pver);             
        getmem(D_qi,pver);             
        getmem(D_nl,pver);             
        getmem(D_ni,pver);             

        getmem(a_cud,pver);            
        getmem(a_cu0,pver);            

        getmem(clrw_old,pver);         
        getmem(clri_old,pver);         

        getmemp(tke,pverp);                     
        getmemp(qtl_flx,pverp);                 
        getmemp(qti_flx,pverp);                 

        getmem(cmfr_det,pver);                
        getmem(qlr_det,pver);                 
        getmem(qir_det,pver);                 

        getreal(landfrac);
        getreal(snowh);

        getpara(wv_para, PARASIZE);

        dma_syn();
        //spara.spa1 = spa1;

        //lwpf_start(kernel_b);

        top_lev = ((core_id % 2) == 0) ? top_lev : pver / 2 + 1; 
        pver = ((core_id % 2) == 0) ? pver / 2 : pver; 
        pverp = ((core_id % 2) == 0) ? pverp / 2 + 1: pverp; 

        //printf_int_(&top_lev);
        //printf_int_(&pver);
        //printf_int_(&pverp);

        wv_para[11] = 0;
        endrunflag = 0;

        mmacro_pcond_cpe_( &ncol, &dt, p, dp, 
                T0, qv0, ql0, qi0, nl0, ni0, 
                A_T, A_qv, A_ql, A_qi, A_nl, A_ni, 
                C_T, C_qv, C_ql, C_qi, C_nl, C_ni, C_qlst, 
                D_T, D_qv, D_ql, D_qi, D_nl, D_ni, 
                a_cud, a_cu0, clrw_old, clri_old, &landfrac, &snowh, 
                tke, qtl_flx, qti_flx, cmfr_det, qlr_det, qir_det, 
                s_tendout, qv_tendout, ql_tendout, qi_tendout, nl_tendout, ni_tendout, 
                qme, qvadj, qladj, qiadj, qllim, qilim, 
                cld, al_st_star, ai_st_star, ql_st_star, qi_st_star, 
                d_rhmin_liq_PBL, d_rhmin_ice_PBL, d_rhmin_liq_det, d_rhmin_ice_det,
                wv_para, &pver, &pverp, &top_lev,
                &i_rhmini, &i_rhminl, &rhmini_const, &rhmaxi_const, &rhminl_const, &rhminl_adj_land_const, &rhminh_const,
                &qmin1, &qmin2, &qmin3, &premib, &premit, 
                &icecrit, &iceopt,  &do_cldice_int );


        //lwpf_stop(kernel_b);

        pver = pver_org;
        pverp = pverp_org;

        //int  i;
        //for(i = 0; i < LAYER; i ++)
        //    qme[i] = spara.qmin1;

        ////------------------------------
        //
        putmem( T0,pver);              
        putmem(qv0,pver);            
        putmem(ql0,pver);            
        putmem(qi0,pver);            
        putmem(nl0,pver);            
        putmem(ni0,pver);            

        putmem( s_tendout,pver);              
        putmem(qv_tendout,pver);            
        putmem(ql_tendout,pver);            
        putmem(qi_tendout,pver);            
        putmem(nl_tendout,pver);            
        putmem(ni_tendout,pver);            

        putmem(qme,  pver);              
        putmem(qvadj,pver);            
        putmem(qladj,pver);            
        putmem(qiadj,pver);            
        putmem(qllim,pver);            
        putmem(qilim,pver);            

        putmem(cld,pver);              
        putmem(al_st_star,pver);       
        putmem(ai_st_star,pver);       
        putmem(ql_st_star,pver);       
        putmem(qi_st_star,pver);       

        putmem(d_rhmin_liq_PBL,pver);       
        putmem(d_rhmin_ice_PBL,pver);       
        putmem(d_rhmin_liq_det,pver);       
        putmem(d_rhmin_ice_det,pver);       

        if(wv_para[11] != 0){
            endrunflag = (int)wv_para[11];
            if(endrunflag == 0) 
                endrunflag = 5;
            if(endrunflag > 8) 
                endrunflag = 10;
        }
        putint(endrunflag);


        //spara.spa2 = spa2;
        //if(core_id == 0)
        //    pe_put(args, &spara, sizeof(mmacro_pcond_args_c)); 

        dma_syn();

    }

    log_data_local_ptr = log_data ; 
    log10_data_local_ptr = log10_data ;  
    pow_data_local_ptr = pow_data ; 
    cos_data_local_ptr = sincos_data ; 
    acos_data_local_ptr = acos_data ; 

    //unsetftz();

//lwpf_stop(kernel_a);
//lwpf_exit(test);

#ifdef XKOPT
#endif

}


