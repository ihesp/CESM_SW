/*************************************************************************
  > File Name: src.cpe/dma_comm.c
  > Author: Xu Kai
  > Created Time: 2018年11月29日 星期四 19时23分22秒
 ************************************************************************/


#include<stdio.h>
//#include<stdlib.h>
#include<slave.h>
//#include "cpe_print.h"
#define GET_COL(col) asm volatile("rcsr    %0, 2" : "=r"(col))
#define GET_ROW(row) asm volatile("rcsr    %0, 1" : "=r"(row))
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_BCASTR(var) REG_PUTR(var, 8)
#define REG_BCASTC(var) REG_PUTC(var, 8)
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define ROWSYN  athread_syn(ROW_SCOPE,0xff)
#define COLSYN  athread_syn(COL_SCOPE,0xff)
#define ALLSYN  athread_syn(ARRAY_SCOPE,0xffff)

void dma_min_(int *lel){
    int cid;// = (*core_id) % 8;
    int rid;// = (*core_id) / 8;
    GET_COL(cid);
    GET_ROW(rid);
    int kmin[8];
    int  i;
    kmin[0] = *lel;

    //printf("id %d\n", *core_id);
    if(rid == 0){
        
        int mink = kmin[0];

        int rnum = 8; 
        // get same column data
        //printf("get %d\n", *core_id);
        for(i = 1; i < rnum; i ++){
            REG_GETC(kmin[i]);
        }

        for(i = 1; i < 8; i ++) 
            if(kmin[0] > kmin[i])
                kmin[0] = kmin[i];

        //printf("syn %d\n", *core_id);
        ROWSYN;
        if(cid == 0){
            // get same row 0  data
            for(i = 1; i < rnum; i ++){
                REG_GETR(kmin[i]);
            }

            for(i = 1; i < 8; i ++) 
                if(kmin[0] > kmin[i])
                    kmin[0] = kmin[i];

            REG_BCASTR(kmin[0]);
        }else{
            REG_PUTR(kmin[0], 0);
            REG_GETR(kmin[0]);
        }
        ROWSYN;
        REG_BCASTC(kmin[0]);
    }else{
        //printf("put %d\n", *core_id);
        REG_PUTC(kmin[0], 0);
        REG_GETC(kmin[0]);
    }
    ALLSYN;
    *lel = kmin[0]; 

}

void dma_max_(int *lel){
    int cid;// = (*core_id) % 8;
    int rid;// = (*core_id) / 8;
    GET_COL(cid);
    GET_ROW(rid);
 
    int kmax[8];
    int  i;
    kmax[0] = *lel;

    if(rid == 0){
        
        int maxk = kmax[0];

        int rnum = 8; 
        // get same column data
        for(i = 1; i < rnum; i ++){
            REG_GETC(kmax[i]);
        }

        for(i = 1; i < 8; i ++) 
            if(kmax[0] < kmax[i])
                kmax[0] = kmax[i];

        if(cid == 0){
            // get same row 0  data
            for(i = 1; i < rnum; i ++){
                REG_GETR(kmax[i]);
            }

            for(i = 1; i < 8; i ++) 
                if(kmax[0] < kmax[i])
                    kmax[0] = kmax[i];

            REG_BCASTR(kmax[0]);
        }else{
            REG_PUTR(kmax[0], 0);
            REG_GETR(kmax[0]);
        }
        REG_BCASTC(kmax[0]);
    }else{
        REG_PUTC(kmax[0], 0);
        REG_GETC(kmax[0]);
    }
    ALLSYN;
    *lel = kmax[0]; 

}
void dma_sum_(int *in, int *out){
    int cid;// = (*core_id) % 8;
    int rid;// = (*core_id) / 8;
    GET_COL(cid);
    GET_ROW(rid);
 
    int sum = *in;
    int sig = 0;
    int  i;

    if(rid == 0){
        int rnum = 8; 
        // get same column data
        for(i = 1; i < rnum; i ++){
            REG_GETC(sig);
            sum += sig;
        }
        if(cid == 0){
            // get same row 0  data
            for(i = 1; i < rnum; i ++){
                REG_GETR(sig);
                sum += sig;
            }
            REG_BCASTR(sum);
        }else{
            REG_PUTR(sum, 0);
            REG_GETR(sum);
        }
        REG_BCASTC(sum);
    }else{
        REG_PUTC(sum, 0);
        REG_GETC(sum);
    }
    ALLSYN;
    *out = *out +  sum; 

}

