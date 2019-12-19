/*************************************************************************
  > File Name: src.cpe/reg_comm.c
  > Author: Xu Kai
  > Created Time: 2019年01月08日 星期二 14时55分54秒
 ************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<slave.h>
#include "dma_macros.h"
#include "cpe_print.h"
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
void get_core_id_(int *core_id){
    int col_id;
    int row_id;
    GET_COL(col_id);
    GET_ROW(row_id);
    *core_id = col_id * 8 + row_id; 

    *core_id = row_id % 2;
}
void get_core_oe_(int *core_ov){
    int row_id;
    GET_ROW(row_id);
    *core_ov = row_id % 2;
}
void regcomm_gotoflag_(int *f1, int *f2){
    int col_id;
    int row_id;
    GET_COL(col_id);
    GET_ROW(row_id);
    int core_id = col_id * 8 + row_id; 
    int t1 = *f1;
    int t2 = *f2;
    int dst; 
    if(core_id % 2 == 0){
        dst = row_id + 1;
        //cpe_printf("put");
        REG_PUTC(t1, dst);
        //cpe_printf("get");
        REG_GETC(t2);
        //cpe_printf("done %d", core_id);
    }else{
        dst = row_id - 1;
        REG_GETC(t2);
        REG_PUTC(t1, dst);
    }
    *f2 = t2;
}
void regcomm_d021_(double *d1){
    int col_id;
    int row_id;
    GET_COL(col_id);
    GET_ROW(row_id);
    int core_id = col_id * 8 + row_id; 
    double t1 = *d1;
    if(core_id % 2 == 0){
        int dst = row_id + 1; 
        REG_PUTC(t1, dst);
    }else{
        REG_GETC(t1);
        *d1 = t1;
    }

}
void regcomm_d120_(double *d1){
    int col_id;
    int row_id;
    GET_COL(col_id);
    GET_ROW(row_id);
    int core_id = col_id * 8 + row_id; 
    double t1 = *d1;

    if(core_id % 2 == 1){
        int dst = row_id - 1; 
        REG_PUTC(t1, dst );
    }else{
        REG_GETC(t1);
        *d1 = t1;
    }

}
void regcomm_sum_(double *sum){
    int col_id;
    int row_id;
    GET_COL(col_id);
    GET_ROW(row_id);
    int core_id = col_id * 8 + row_id; 

    double sum_tmp = *sum;
    double sum_t1 = *sum;
    int dst;
    if(core_id % 2 == 0){
        dst = row_id + 1;
        REG_PUTC(sum_t1, dst);
        REG_GETC(sum_tmp);
    }else{
        dst = row_id - 1;
        REG_GETC(sum_tmp);
        REG_PUTC(sum_t1, dst);
    }
    *sum = sum_t1 + sum_tmp;

}
