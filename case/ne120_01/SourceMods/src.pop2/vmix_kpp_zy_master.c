#include<stdio.h>
#include<athread.h>
#include "vmix_kpp_zy_struct.h"

extern SLAVE_FUN(fun1_blmix)(),SLAVE_FUN(fun5_blmix)(),SLAVE_FUN(fun8_blmix)();

void loop1_blmix_(double* HBLT ,double* BFSFC,
                  double* USTAR,double* WM   ,double* WS   ,int *nx_block,int* ny_block,int*my_task)
{   
    struct kpp_ddmix s2;
 
     s2.param[0]=*nx_block;
     s2.param[1]=*ny_block;
     s2.param[2]=athread_get_max_threads();
     
     s2.arr[0]=HBLT;
     s2.arr[1]=BFSFC;
     s2.arr[2]=USTAR;
     s2.arr[3]=WM;
     s2.arr[4]=WS; 

     athread_spawn(fun1_blmix,&s2);
     athread_join();     
        
}
void   loop5_blmix_(int *km,int *ny_block,int *nx_block,double *zgrid,double* DZT,double*HBLT,double *hwide,double *BFSFC,
                    double* USTAR,double*BLMC,double* GAT1,double* DAT1,double* STABLE,double* GHAT,double *cg,int* proc,double* p33)
{
       double *pkg;
       int k,km_size=*km;
        
       pkg=(double *)malloc(8*(2*(km_size +2)));

       for(k=0;k < km_size+1;k++)
       {
          pkg[k]             =zgrid[k];
          pkg[km_size+1+k]   =hwide[k];
       }

      struct kpp_ddmix s2;
     s2.param[0]=*nx_block;
     s2.param[1]=*ny_block;
     s2.param[2]=*km;
     s2.param[3]= athread_get_max_threads();
     s2.param[4]=*proc;

     s2.dparam[0]=*cg;
     s2.dparam[1]=*p33;

     s2.arr[0]=DZT;
     s2.arr[1]=HBLT;
     s2.arr[2]=BFSFC;
     s2.arr[3]=USTAR;
     s2.arr[4]=BLMC;
     s2.arr[5]=GAT1;
     s2.arr[6]=DAT1;
     s2.arr[7]=STABLE;
     s2.arr[8]=GHAT;
   
     s2.std[0]=pkg;


     athread_spawn(fun5_blmix,&s2);
     athread_join();

     free(pkg);
}


void  loop8_blmix(int* km,int* ny_block,int* nx_block,int* KBL,double *VISC,double *BLMC,double *GHAT,double *VDC)
{
/*
    struct kpp_ddmix s2; 
     s2.param[0]=*nx_block;
     s2.param[1]=*ny_block;
     s2.param[2]=*km;
     s2.arr[0]=VISC;
     s2.arr[1]=BLMC;
     s2.arr[2]=GHAT;
     s2.arr[3]=VDC;
     s2.arr_int[0]=KBL;
     athread_spawn(fun8_blmix,&s2);
     athread_join();     
*/     
}
