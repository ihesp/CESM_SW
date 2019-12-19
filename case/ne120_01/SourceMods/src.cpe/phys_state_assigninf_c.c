/*************************************************************************
	> File Name: src.cpe/phys_state_assigninf_c.c
	> Author: Xu Kai
	> Created Time: 2018年12月20日 星期四 16时50分04秒
 ************************************************************************/
	
	
#include<stdio.h>
#include<stdlib.h>
#include"dma_macros.h"

typedef struct assigninf_para_c
{
    int psetcols;
    int pver;
    int pcnst;
    
    double infr8;
    void *lat       ; 
    void *lon       ; 
    void *ulat      ; 
    void *ulon      ; 
    void *ps        ; 
    void *psdry     ; 
    void *phis      ; 
    void *t         ; 
    void *u         ; 
    void *v         ; 
    void *s         ; 
    void *omega     ; 
    void *pmid      ; 
    void *pmiddry   ; 
    void *pdel      ; 
    void *pdeldry   ; 
    void *rpdel     ; 
    void *rpdeldry  ; 
    void *lnpmid    ; 
    void *lnpmiddry ; 
    void *exner     ; 
    void *zm        ; 
    void *q         ; 
    void *pint      ; 
    void *pintdry   ; 
    void *lnpint    ; 
    void *lnpintdry ; 
    void *zi        ; 
    void *te_ini    ; 
    void *te_cur    ; 
    void *tw_ini    ; 
    void *tw_cur    ; 
} assigninf_args_cc;

void assigninf_parallel_(void *assigninf_para)
{
    assigninf_args_cc assigninf_spara;

    dma_init();
    pe_get(assigninf_para, &assigninf_spara, sizeof(assigninf_args_cc));
    dma_syn();

    int psetcols = assigninf_spara.psetcols;
    int pver     = assigninf_spara.pver;
    int pcnst    = assigninf_spara.pcnst;
    int myid = _COL * 8 + _ROW;

    void *args[32];
    args[0] =  (void *)assigninf_spara.t;
    args[1] =  (void *)assigninf_spara.u;
    args[2] =  (void *)assigninf_spara.v;
    args[3] =  (void *)assigninf_spara.s;
    args[4] =  (void *)assigninf_spara.omega;
    args[5] =  (void *)assigninf_spara.pmid;
    args[6] =  (void *)assigninf_spara.pmiddry;
    args[7] =  (void *)assigninf_spara.pdel;
    args[8] =  (void *)assigninf_spara.pdeldry;
    args[9] =  (void *)assigninf_spara.rpdel;
    args[10] = (void *)assigninf_spara.rpdeldry;
    args[11] = (void *)assigninf_spara.lnpmid;
    args[12] = (void *)assigninf_spara.lnpmiddry;
    args[13] = (void *)assigninf_spara.exner;
    args[14] = (void *)assigninf_spara.zm;
    args[15] = (void *)assigninf_spara.pint;
    args[16] = (void *)assigninf_spara.pintdry;
    args[17] = (void *)assigninf_spara.lnpint;
    args[18] = (void *)assigninf_spara.lnpintdry;
    args[19] = (void *)assigninf_spara.zi;

    int i;
    if ( myid < 20 )
    {
        double ldm_array[1024];
        for ( i = 0; i < psetcols*pver; i ++ )
            ldm_array[i] = assigninf_spara.infr8;
        int psize;
        if ( myid < 15 ) psize = psetcols*pver*sizeof(double);
        else psize = psetcols*(pver+1)*sizeof(double);
        pe_put((void*)(args[myid])   , ldm_array, psize);
        dma_syn();
    }
    if ( myid >= 20 && myid < 20+pcnst )
    {
        double ldm_array[4096];
        for ( i = 0; i < psetcols*pver; i ++ )
            ldm_array[i] =  assigninf_spara.infr8;
        int psize = psetcols*pver*sizeof(double);
        int offset = (myid - 20) * psetcols*pver*sizeof(double);
        pe_put(assigninf_spara.q + offset   , ldm_array, psize);
        dma_syn();
    }
    if ( myid == 20+pcnst )
    {
        double ldm_array[1024];
        for ( i = 0; i < psetcols; i ++ )
            ldm_array[i] =  assigninf_spara.infr8;
        int psize = psetcols * sizeof(double);
        pe_put(assigninf_spara.lat   , ldm_array, psize);
        pe_put(assigninf_spara.lon   , ldm_array, psize);
        pe_put(assigninf_spara.ulat  , ldm_array, psize);
        pe_put(assigninf_spara.ulon  , ldm_array, psize);
        pe_put(assigninf_spara.ps    , ldm_array, psize);
        pe_put(assigninf_spara.psdry , ldm_array, psize);
        pe_put(assigninf_spara.phis  , ldm_array, psize);
        pe_put(assigninf_spara.te_ini, ldm_array, psize);
        pe_put(assigninf_spara.te_cur, ldm_array, psize);
        pe_put(assigninf_spara.tw_ini, ldm_array, psize);
        pe_put(assigninf_spara.tw_cur, ldm_array, psize);
        dma_syn();
    }
}
        



