/*************************************************************************
  > File Name: src.cam/convect_deep_cpe_c.c
  > Author: Xu Kai
  > Created Time: 2018年11月07日 星期三 14时04分04秒
 ************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>


#include<math.h>
#include "dma_macros.h"
//#include "lwpf_cpe.h"
#include "cpe_print.h"

#define REAL double
#define LAYER 32
typedef struct zm_convr_args_c{
    void *t_trans; 
    void *qh_trans;
    void *pap_trans;     
    void *paph_trans;
    void *dpp_trans;
    void *zm_trans;
    void *zi_trans;

    void *qtnd_trans;
    void *heat_trans;
    void *mcon_trans;
    void *dlf_trans;
    void *pflx_trans;
    void *cme_trans;
    void *zdu_trans;
    void *rprd_trans;
    void *mu_trans;
    void *eu_trans;
    void *du_trans;
    void *md_trans;
    void *ed_trans;
    void *dp_trans;       
    //inout
    void *ql_trans;
    void *org_trans;
    void *orgt_trans;
    void *org2d_trans;
    //in 
    void *geos;   
    void *pblh;  
    void *tpert;
    void *landfrac;
    //out
    void *cape;
    void *dsubcld;
    void *jctop;
    void *jcbot;
    void *prec;
    void *rliq;

    //inout
    //integer
    void *ideep;
    void *maxg;
    void *jt;
    void *lengath;

    void *lcl;
    void *lel;
    void *lon;
    void *mx;

    REAL delt;
    REAL grav;
    REAL rgrav;
    REAL cpres;
    REAL rl;
    REAL rgas;
    REAL cpliq;
    REAL cpwv;
    REAL tfreez; 
    REAL eps1;
    REAL rh2o;
    REAL latice;
    REAL c0_ocn;
    REAL c0_lnd;

    
    //global
    int pcols;
    int pver;
    int pverp;
    //in
    int lchnk;
    int ncol;
    int zm_org;
    int limcnv;
    int no_deep_pbl;
} zm_convr_args_cc;

#include "math_data.h"
void zm_convr_parallel_(void *zm_convr_para)
{

    dma_init(); 
    double log_local_data[log_data_len] ;
    double log10_local_data[log10_data_len] ;
    double pow_local_data[pow_data_len] ;

    if (_MYID == 0) {
        
        bcast_get(log_data , log_local_data , log_data_bytes) ;
        bcast_get(log10_data , log10_local_data , log10_data_bytes) ;
        bcast_get(pow_data , pow_local_data , pow_data_bytes) ;
        dma_syn() ;
        //cpe_printf("%p\n", log_data); 
    }
    athread_syn(ARRAY_SCOPE ,0xffff) ;
    log_data_local_ptr = log_local_data ;
    log10_data_local_ptr = log10_local_data ;
    pow_data_local_ptr = pow_local_data ;



    ////lwpf_ZM_CONVR_enter();
    //lwpf_enter(ZM_CONVR);
    ////lwpf_A_start();
    //lwpf_start(A);
    zm_convr_args_cc zm_convr_spara; 


    pe_get(zm_convr_para, &zm_convr_spara, sizeof(zm_convr_args_cc)); 
    dma_syn();

    int core_id = athread_get_core(-1);
    core_id = (core_id % 8) * 8 + core_id /  8;

    REAL t[LAYER]; 
    REAL qh[LAYER];
    REAL pap[LAYER];     
    REAL paph[LAYER];
    REAL dpp[LAYER];
    REAL zm[LAYER];
    REAL zi[LAYER];

    REAL qtnd[LAYER];
    REAL heat[LAYER];
    REAL mcon[LAYER];
    REAL dlf[LAYER];
    REAL pflx[LAYER];
    REAL cme[LAYER];
    REAL zdu[LAYER];
    REAL rprd[LAYER];
    REAL mu[LAYER];
    REAL eu[LAYER];
    REAL du[LAYER];
    REAL md[LAYER];
    REAL ed[LAYER];
    REAL dp[LAYER];
    REAL ql[LAYER];
    REAL org[LAYER];
    REAL orgt[LAYER];
    REAL org2d[LAYER];

    //in 
    REAL geos;   
    REAL pblh;  
    REAL tpert;
    REAL landfrac;
    //out
    REAL cape;
    REAL dsubcld;
    REAL jctop;
    REAL jcbot;
    REAL prec;
    REAL rliq;

    int ideep;
    int maxg;
    int jt;
    int lengath;
    int lcl;
    int lel; 
    int lon;
    int mx;

    REAL delt;
    REAL grav;
    REAL rgrav;
    REAL cpres;
    REAL rl;
    REAL rgas;
    REAL cpliq;
    REAL cpwv;
    REAL tfreez;
    REAL eps1;
    REAL rh2o;
    REAL latice;
    REAL c0_ocn;
    REAL c0_lnd;
    

    //global
    int pcols;
    int pver;
    int pverp;

    //in
    int lchnk;
    int ncol;
    int zm_org;
    int limcnv;
    int no_deep_pbl;
   //pe_get(log_data, log_local_data, log_data_bytes);
    //dma_syn() ;
    //memcpy(log_local_data, log_data, log_data_bytes);

#define OFFSET(dim) (core_id * (dim) * sizeof(REAL)) 
#define OFFSETINT(dim) (core_id * (dim) * sizeof(int)) 

    pcols = zm_convr_spara.pcols;
    pver = zm_convr_spara.pver;
    pverp = zm_convr_spara.pverp;
    lchnk = zm_convr_spara.lchnk;
    ncol = zm_convr_spara.ncol;
    zm_org = zm_convr_spara.zm_org;
    limcnv = zm_convr_spara.limcnv;
    ideep = 0;
    int tttt[4];
    //if(pver != 30)
    //    printf("wrong pver != 30");

    ////lwpf_A_stop();
    //lwpf_stop(A);

    if(core_id < ncol){
        //lwpf_start(B);
        ////lwpf_B_start();
        pe_get(zm_convr_spara.t_trans + OFFSET(pver),  t, pver * sizeof(REAL)); 
        pe_get(zm_convr_spara.qh_trans + OFFSET(pver),  qh, pver * sizeof(REAL)); 
        pe_get(zm_convr_spara.pap_trans + OFFSET(pver),  pap, pver * sizeof(REAL)); 
        pe_get(zm_convr_spara.paph_trans + OFFSET(pver + 1),  paph, (pver + 1) * sizeof(REAL)); 
        pe_get(zm_convr_spara.dpp_trans + OFFSET(pver),  dpp, pver * sizeof(REAL)); 
        pe_get(zm_convr_spara.zm_trans + OFFSET(pver),  zm, pver * sizeof(REAL)); 
        pe_get(zm_convr_spara.zi_trans + OFFSET(pver + 1),  zi, (pver + 1) * sizeof(REAL)); 
        pe_get(zm_convr_spara.ql_trans + OFFSET(pver),  ql, pver * sizeof(REAL)); 

        //pe_get(zm_convr_spara.eu_trans + OFFSET(pver),  eu, pver * sizeof(REAL)); 
        //pe_get(zm_convr_spara.mu_trans + OFFSET(pver),  mu, pver * sizeof(REAL)); 
        //pe_get(zm_convr_spara.du_trans + OFFSET(pver),  du, pver * sizeof(REAL)); 
        //pe_get(zm_convr_spara.md_trans + OFFSET(pver),  md, pver * sizeof(REAL)); 
        //pe_get(zm_convr_spara.ed_trans + OFFSET(pver),  ed, pver * sizeof(REAL)); 
        dma_syn();
        if(zm_org == 1){
            pe_get(zm_convr_spara.org_trans + OFFSET(pver),  org, pver * sizeof(REAL)); 
            pe_get(zm_convr_spara.orgt_trans + OFFSET(pver),  orgt, pver * sizeof(REAL)); 
            pe_get(zm_convr_spara.org2d_trans + OFFSET(pver),  org2d, pver * sizeof(REAL)); 
        }

       pe_get(zm_convr_spara.geos + OFFSET(1),  &geos, sizeof(REAL)); 
       pe_get(zm_convr_spara.pblh + OFFSET(1),  &pblh, sizeof(REAL)); 
       pe_get(zm_convr_spara.tpert + OFFSET(1),  &tpert, sizeof(REAL)); 
       pe_get(zm_convr_spara.landfrac + OFFSET(1),  &landfrac, sizeof(REAL)); 

       //pe_get(zm_convr_spara.ideep + OFFSETINT(1),  &ideep, sizeof(int)); 
       pe_get(zm_convr_spara.maxg + OFFSETINT(1),  &maxg, sizeof(int)); 
       pe_get(zm_convr_spara.jt + OFFSETINT(1),  &jt, sizeof(int)); 
       dma_syn();

        delt = zm_convr_spara.delt;
        grav = zm_convr_spara.grav;
        rgrav = zm_convr_spara.rgrav;
        cpres = zm_convr_spara.cpres;
        rl  = zm_convr_spara.rl;
        rgas = zm_convr_spara.rgas;
        cpliq = zm_convr_spara.cpliq;
        cpwv = zm_convr_spara.cpwv;
        tfreez = zm_convr_spara.tfreez;
        eps1 = zm_convr_spara.eps1;
        rh2o = zm_convr_spara.rh2o;
        latice = zm_convr_spara.latice;
        c0_ocn = zm_convr_spara.c0_ocn;
        c0_lnd = zm_convr_spara.c0_lnd;

        //lengath = zm_convr_spara.lengath;
        no_deep_pbl = zm_convr_spara.no_deep_pbl;

        ////lwpf_B_stop();
        //lwpf_stop(B);
        ////lwpf_C_start();
        //lwpf_start(C);
        zm_convr_compute_(&lchnk   ,&ncol    , 
                t       ,qh      ,&prec    ,&jctop   ,&jcbot   , 
                &pblh    ,zm      ,&geos    ,zi      ,qtnd    , 
                heat    ,pap     ,paph    ,dpp     , 
                &delt    ,mcon    ,cme     ,&cape    , 
                &tpert   ,dlf     ,pflx    ,zdu     ,rprd    , 
                mu      ,md      ,du      ,eu      ,ed      , 
                dp      ,&dsubcld ,&jt      ,&maxg    ,&ideep   , 
                &lengath ,ql      ,&rliq    ,&landfrac,          
                &pver   ,&pcols, &pverp, &zm_org, &limcnv, org, 
                orgt, org2d, &grav, &rgrav, &cpres, &rl, &rgas,&cpliq, 
                &cpwv, &tfreez, &eps1, &rh2o, &latice, &c0_ocn, &c0_lnd, 
                &no_deep_pbl, tttt);


#define putmem(arr, num) \
    pe_put(zm_convr_spara.arr##_trans + OFFSET(num), arr, (num)*sizeof(REAL))

#define putmemreal(sig) \
    pe_put(zm_convr_spara.sig + OFFSET(1), &sig, (1)*sizeof(REAL))
#define putmemint(sig) \
    pe_put(zm_convr_spara.sig + OFFSETINT(1), &sig, (1)*sizeof(int))
        //lwpf_stop(C);
        ////lwpf_C_stop();
        //
      putmem(qtnd,  pver);
      putmem(heat,  pver);
      putmem(dlf,   pver);
      putmem(cme,   pver);

      //putmem(mcon,  pverp);
      //putmem(pflx,  pverp);
      //putmem(zdu,   pver);
      //putmem(rprd,  pver);
      //putmem(mu,    pver);
      //putmem(eu,    pver);
      //putmem(du,    pver);
      //putmem(md,    pver);
      //putmem(ed,    pver);
      //putmem(dp,    pver);
      //putmem(ql,    pver);

      //putmem(t,    pver);
      //putmem(qh,    pver);
      //putmem(pap,    pver);
      //putmem(dpp,    pver);
      //putmem(zm,    pver);

      dma_syn();

      //putmemreal(jctop);
      //putmemreal(jcbot);
      //putmemreal(prec);
      //putmemreal(rliq);
      if(zm_org == 1){
        putmem(org,   pver);
        putmem(orgt,  pver);
        putmem(org2d, pver);
      }
      //putmemint(ideep);
      lcl = tttt[0];
      lel = tttt[1];
      lon = tttt[2];
      mx = tttt[3];
      putmemint(lcl);
      putmemint(lel);
      putmemint(lon);
      putmemint(mx);

      ////putmemint(maxg);
      ////putmemint(jt);
      //lengath = ideep;
      //putmemint(lengath);
      dma_syn();
    }else{
#if 0
      int kcount_sig = 0;
      lengath = 0;
      dma_sum_(&kcount_sig, &lengath);

      int lel = pverp;  
      dma_min_(&lel); 
      lel = 1;  
      dma_max_(&lel); 

      int i = 0;
      int msg = limcnv - 1; 
      int kcount  = 0;
      kcount_sig = 0;

      for(i = pver; i >= msg + 2; i --){
          dma_sum_(&kcount_sig, &kcount);
           if(kcount >= lengath)
               break;
      }
#endif

    }
      //cpe_printf("begin\n");
      //dma_sum_(&kcount_sig, &lengath);
      //cpe_printf("end\n");
      
    log_data_local_ptr = log_data ;
    //log10_data_local_ptr = log10_data ;
    pow_data_local_ptr = pow_data ;



    ////lwpf_ZM_CONVR_exit();
    //lwpf_exit(ZM_CONVR);

#if 0 
#endif
}
#undef REAL
#undef LAYER
#undef OFFSET
#undef OFFSETINT
#undef putmem
#undef putmemreal
#undef putmemint
