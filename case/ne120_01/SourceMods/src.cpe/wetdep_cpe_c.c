/*************************************************************************
	> File Name: wetdep_cpe_c.c
	> Author: Xu Kai
	> Created Time: 2019年01月12日 星期六 15时52分43秒
 ************************************************************************/
	
	
#include<stdio.h>
#include<stdlib.h>
#include"dma_macros.h"
#include<stdbool.h>
#include<stdarg.h>
#include "cpe_print.h"
#include<math.h>


#define REAL double
#define L1 32
#define L2 30

#define COL(x) (x & 0x07)
#define ROW(x) ((x & 0x38) >> 3)

typedef struct wetdepa_v2_para_c
{
    int version            ;
    int ncol               ;
    int pcols              ;
    int pver               ;
    int is_strat_cloudborne;
    REAL deltat             ;
    REAL sol_fact           ;
    REAL sol_facti_in       ;
    REAL gravit             ;

    void * p            ;
    void * q            ;
    void * pdel         ;
    void * cldt         ;
    void * cldc         ;
    void * cmfdqr       ;
    void * evapc        ;
    void * conicw       ;
    void * cwat         ;
    void * precs        ;
    void * conds        ;
    void * evaps        ;
    void * cldvcu       ;
    void * cldvst       ;
    void * dlf          ;
    void * tracer       ;
    void * scavcoef     ;
    void * sol_factic_in;
    void * qqcw         ;
    void * f_act_conv   ;
    void * scavt        ;
    void * iscavt       ;
    void * fracis       ;
    void * icscavt      ;
    void * isscavt      ;
    void * bcscavt      ;
    void * bsscavt      ;
} wetdepa_v2_args_cc;

typedef struct transpose_para_c
{
    int version            ;
    int pcols              ;
    int pver               ;

    void * p            ;
    void * q            ;
    void * pdel         ;
    void * cldt         ;
    void * cldc         ;
    void * cmfdqr       ;
    void * evapc        ;
    void * conicw       ;
    void * cwat         ;
    void * precs        ;
    void * conds        ;
    void * evaps        ;
    void * cldvcu       ;
    void * cldvst       ;
    void * dlf          ;
    void * tracer       ;
    void * scavcoef     ;
    void * sol_factic_in;
    void * qqcw         ;
    void * f_act_conv   ;
    void * scavt        ;
    void * iscavt       ;
    void * fracis       ;
    void * icscavt      ;
    void * isscavt      ;
    void * bcscavt      ;
    void * bsscavt      ;

    void * p_permute            ;
    void * q_permute            ;
    void * pdel_permute         ;
    void * cldt_permute         ;
    void * cldc_permute         ;
    void * cmfdqr_permute       ;
    void * evapc_permute        ;
    void * conicw_permute       ;
    void * cwat_permute         ;
    void * precs_permute        ;
    void * conds_permute        ;
    void * evaps_permute        ;
    void * cldvcu_permute       ;
    void * cldvst_permute       ;
    void * dlf_permute          ;
    void * tracer_permute       ;
    void * scavcoef_permute     ;
    void * sol_factic_in_permute;
    void * qqcw_permute         ;
    void * f_act_conv_permute   ;
    void * scavt_permute        ;
    void * iscavt_permute       ;
    void * fracis_permute       ;
    void * icscavt_permute      ;
    void * isscavt_permute      ;
    void * bcscavt_permute      ;
    void * bsscavt_permute      ;
} transpose_args_cc;
void transpose_parallel_(void *transpose_para)
{
    transpose_args_cc transpose_spara;
    dma_init();

    pe_get(transpose_para, &transpose_spara, sizeof(transpose_args_cc));
    dma_syn();

    REAL *ptr[27];
    REAL *ptr_[27];

    ptr[ 0] = (REAL*)transpose_spara.p            ;
    ptr[ 1] = (REAL*)transpose_spara.q            ;
    ptr[ 2] = (REAL*)transpose_spara.pdel         ;
    ptr[ 3] = (REAL*)transpose_spara.cldt         ;
    ptr[ 4] = (REAL*)transpose_spara.cldc         ;
    ptr[ 5] = (REAL*)transpose_spara.cmfdqr       ;
    ptr[ 6] = (REAL*)transpose_spara.evapc        ;
    ptr[ 7] = (REAL*)transpose_spara.conicw       ;
    ptr[ 8] = (REAL*)transpose_spara.cwat         ;
    ptr[ 9] = (REAL*)transpose_spara.precs        ;
    ptr[10] = (REAL*)transpose_spara.conds        ;
    ptr[11] = (REAL*)transpose_spara.evaps        ;
    ptr[12] = (REAL*)transpose_spara.cldvcu       ;
    ptr[13] = (REAL*)transpose_spara.cldvst       ;
    ptr[14] = (REAL*)transpose_spara.dlf          ;
    ptr[15] = (REAL*)transpose_spara.tracer       ;
    ptr[16] = (REAL*)transpose_spara.scavcoef     ;
    ptr[17] = (REAL*)transpose_spara.sol_factic_in;
    ptr[18] = (REAL*)transpose_spara.qqcw         ;
    ptr[19] = (REAL*)transpose_spara.f_act_conv   ;
    ptr[20] = (REAL*)transpose_spara.scavt        ;
    ptr[21] = (REAL*)transpose_spara.iscavt       ;
    ptr[22] = (REAL*)transpose_spara.fracis       ;
    ptr[23] = (REAL*)transpose_spara.icscavt      ;
    ptr[24] = (REAL*)transpose_spara.isscavt      ;
    ptr[25] = (REAL*)transpose_spara.bcscavt      ;
    ptr[26] = (REAL*)transpose_spara.bsscavt      ;
    ptr_[ 0] = (REAL*)transpose_spara.p_permute            ;
    ptr_[ 1] = (REAL*)transpose_spara.q_permute            ;
    ptr_[ 2] = (REAL*)transpose_spara.pdel_permute         ;
    ptr_[ 3] = (REAL*)transpose_spara.cldt_permute         ;
    ptr_[ 4] = (REAL*)transpose_spara.cldc_permute         ;
    ptr_[ 5] = (REAL*)transpose_spara.cmfdqr_permute       ;
    ptr_[ 6] = (REAL*)transpose_spara.evapc_permute        ;
    ptr_[ 7] = (REAL*)transpose_spara.conicw_permute       ;
    ptr_[ 8] = (REAL*)transpose_spara.cwat_permute         ;
    ptr_[ 9] = (REAL*)transpose_spara.precs_permute        ;
    ptr_[10] = (REAL*)transpose_spara.conds_permute        ;
    ptr_[11] = (REAL*)transpose_spara.evaps_permute        ;
    ptr_[12] = (REAL*)transpose_spara.cldvcu_permute       ;
    ptr_[13] = (REAL*)transpose_spara.cldvst_permute       ;
    ptr_[14] = (REAL*)transpose_spara.dlf_permute          ;
    ptr_[15] = (REAL*)transpose_spara.tracer_permute       ;
    ptr_[16] = (REAL*)transpose_spara.scavcoef_permute     ;
    ptr_[17] = (REAL*)transpose_spara.sol_factic_in_permute;
    ptr_[18] = (REAL*)transpose_spara.qqcw_permute         ;
    ptr_[19] = (REAL*)transpose_spara.f_act_conv_permute   ;
    ptr_[20] = (REAL*)transpose_spara.scavt_permute        ;
    ptr_[21] = (REAL*)transpose_spara.iscavt_permute       ;
    ptr_[22] = (REAL*)transpose_spara.fracis_permute       ;
    ptr_[23] = (REAL*)transpose_spara.icscavt_permute      ;
    ptr_[24] = (REAL*)transpose_spara.isscavt_permute      ;
    ptr_[25] = (REAL*)transpose_spara.bcscavt_permute      ;
    ptr_[26] = (REAL*)transpose_spara.bsscavt_permute      ;
    int version = transpose_spara.version;
    int pcols   = transpose_spara.pcols;
    int pver    = transpose_spara.pver;

    int myid = COL(_MYID)*8 + ROW(_MYID);

    REAL array1[1024];
    REAL array2[1024];

    int i, k;
    int len;
    if ( version == 1 ) len = 19;
    else len = 17;
    if ( version != 3 )
    {
        if ( myid <= len )
        {
            pe_get(ptr[myid], array1, pcols*pver*sizeof(REAL));
            dma_syn();

            for ( i = 0; i < pcols; i ++ )
            {
                for ( k = 0; k < pver; k ++ )
                {
                    array2[i*pver+k] = array1[k*pcols+i];
                }
            }
            pe_put(ptr_[myid], array2, pcols*pver*sizeof(REAL));
            dma_syn();
        }
    }
    else
    {
        if ( myid >= 20 && myid <= 26 )
        {
            pe_get(ptr_[myid], array1, pcols*pver*sizeof(REAL));
            dma_syn();

            for ( i = 0; i < pcols; i ++ )
            {
                for ( k = 0; k < pver; k ++ )
                {
                    array2[k*pcols+i] = array1[i*pver+k];
                }
            }
            pe_put(ptr[myid], array2, pcols*pver*sizeof(REAL));
            dma_syn();
        }
    }
}

void wetdepa_v2_parallel_(void *wetdepa_v2_para)
{
    //set_func_ptr_();
    wetdepa_v2_args_cc wetdepa_v2_spara;
    dma_init();

    pe_get(wetdepa_v2_para, &wetdepa_v2_spara, sizeof(wetdepa_v2_args_cc));
    dma_syn();

    REAL p            [L2];
    REAL q            [L2];
    REAL pdel         [L2];
    REAL cldt         [L2];
    REAL cldc         [L2];
    REAL cmfdqr       [L2];
    REAL evapc        [L2];
    REAL conicw       [L2];
    REAL cwat         [L2];
    REAL precs        [L2];
    REAL conds        [L2];
    REAL evaps        [L2];
    REAL cldvcu       [L2];
    REAL cldvst       [L2];
    REAL dlf          [L2];
    REAL tracer       [L2];
    REAL scavcoef     [L2];
    REAL qqcw         [L2];
    REAL f_act_conv   [L2];
    REAL sol_factic_in[L2];
    REAL scavt        [L2];
    REAL iscavt       [L2];
    REAL fracis       [L2];
    REAL icscavt      [L2];
    REAL isscavt      [L2];
    REAL bcscavt      [L2];
    REAL bsscavt      [L2];

    int  version            = wetdepa_v2_spara.version            ; 
    int  ncol               = wetdepa_v2_spara.ncol               ; 
    int  pcols              = wetdepa_v2_spara.pcols              ; 
    int  pver               = wetdepa_v2_spara.pver               ; 
    int  is_strat_cloudborne= wetdepa_v2_spara.is_strat_cloudborne; 
    REAL deltat             = wetdepa_v2_spara.deltat             ; 
    REAL sol_fact           = wetdepa_v2_spara.sol_fact           ; 
    REAL sol_facti_in       = wetdepa_v2_spara.sol_facti_in       ; 
    REAL gravit             = wetdepa_v2_spara.gravit             ; 

    int myid = COL(_MYID)*8 + ROW(_MYID);

#define OFFSET(dim) (myid * (dim) * sizeof(REAL))
    if ( myid < ncol )
    {
        pe_get(wetdepa_v2_spara.p            + OFFSET(pver), p            , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.q            + OFFSET(pver), q            , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.pdel         + OFFSET(pver), pdel         , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.cldt         + OFFSET(pver), cldt         , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.cldc         + OFFSET(pver), cldc         , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.cmfdqr       + OFFSET(pver), cmfdqr       , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.evapc        + OFFSET(pver), evapc        , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.conicw       + OFFSET(pver), conicw       , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.cwat         + OFFSET(pver), cwat         , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.precs        + OFFSET(pver), precs        , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.conds        + OFFSET(pver), conds        , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.evaps        + OFFSET(pver), evaps        , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.cldvcu       + OFFSET(pver), cldvcu       , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.cldvst       + OFFSET(pver), cldvst       , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.dlf          + OFFSET(pver), dlf          , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.tracer       + OFFSET(pver), tracer       , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.scavcoef     + OFFSET(pver), scavcoef     , pver*sizeof(REAL));
        pe_get(wetdepa_v2_spara.sol_factic_in+ OFFSET(pver), sol_factic_in, pver*sizeof(REAL));
        if ( version == 1 )
        {
            pe_get(wetdepa_v2_spara.qqcw      + OFFSET(pver), qqcw      , pver*sizeof(REAL));
            pe_get(wetdepa_v2_spara.f_act_conv+ OFFSET(pver), f_act_conv, pver*sizeof(REAL));
        }
        dma_syn();

        //if ( version == 1 )
        //{
        //    wetdepa_v2_compute_(
        //            &pcols, &pver, &gravit,                                 \
        //            p, q, pdel, cldt, cldc,                                 \
        //            cmfdqr, evapc, conicw, precs, conds,                    \
        //            evaps, cwat, tracer, &deltat, scavt,                    \
        //            iscavt, cldvcu, cldvst, dlf, fracis,                    \
        //            &sol_fact, &ncol, scavcoef, &is_strat_cloudborne, qqcw, \
        //            f_act_conv, icscavt, isscavt, bcscavt, bsscavt,         \
        //            &sol_facti_in, sol_factic_in );
        //}
        //else
        //{
            wetdepa_v2_compute_(
                    &pcols, &pver, &gravit,                                 \
                    p, q, pdel, cldt, cldc,                                 \
                    cmfdqr, evapc, conicw, precs, conds,                    \
                    evaps, cwat, tracer, &deltat, scavt,                    \
                    iscavt, cldvcu, cldvst, dlf, fracis,                    \
                    &sol_fact, &ncol, scavcoef, &is_strat_cloudborne, qqcw, \
                    f_act_conv, icscavt, isscavt, bcscavt, bsscavt,         \
                    &sol_facti_in, sol_factic_in );
        //}

        pe_put(wetdepa_v2_spara.scavt  + OFFSET(pver), scavt  , pver*sizeof(REAL));
        pe_put(wetdepa_v2_spara.iscavt + OFFSET(pver), iscavt , pver*sizeof(REAL));
        pe_put(wetdepa_v2_spara.fracis + OFFSET(pver), fracis , pver*sizeof(REAL));
        pe_put(wetdepa_v2_spara.icscavt+ OFFSET(pver), icscavt, pver*sizeof(REAL));
        pe_put(wetdepa_v2_spara.isscavt+ OFFSET(pver), isscavt, pver*sizeof(REAL));
        pe_put(wetdepa_v2_spara.bcscavt+ OFFSET(pver), bcscavt, pver*sizeof(REAL));
        pe_put(wetdepa_v2_spara.bsscavt+ OFFSET(pver), bsscavt, pver*sizeof(REAL));
        dma_syn();


    }

}



