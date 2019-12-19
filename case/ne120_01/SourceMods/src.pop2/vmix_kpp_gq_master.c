#include<stdio.h>
#include<athread.h>
#include<math.h>
#include "vmix_kpp_gq_struct.h"

struct ri_Loop01_Para Loop01_Para;
struct ri_Loop02_Para Loop02_Para;
struct ri_Loop03_Para Loop03_Para;

extern SLAVE_FUN(ri_loop01_fun)();
extern SLAVE_FUN(ri_loop02_fun)();
extern SLAVE_FUN(ri_loop03_fun)();


void ri_loop01_fun_(int* nx_block, int* ny_block,int* bid,int*km,int*par,int*Ddyn,
                      double* p5,double* c0, double*eps,int*KMT,
                      double*UUU, double*VVV, double*FRI,double*DZU,double*DZT,double*VSHEAR,
                      double*DBLOC,double*zgrid,double*RI_LOC,double*WORK0,
                      double *AT0,double *ATS,double *ATW,double *ATSW)
                     // int *my_task)

{
    Loop01_Para.param[0]=*nx_block;
    Loop01_Para.param[1]=*ny_block;
    Loop01_Para.param[2]=*bid;
    Loop01_Para.param[3]=*km;
    Loop01_Para.param[4]=*par;
    Loop01_Para.param[5]=*Ddyn;
    Loop01_Para.param[6]=athread_get_max_threads();
    
   // Loop01_Para.param[7]=*my_task;
    

    Loop01_Para.dparam[0]=*p5;
    Loop01_Para.dparam[1]=*c0;
    Loop01_Para.dparam[2]=*eps;

    Loop01_Para.duse_arr[0]=KMT;
    
    Loop01_Para.use_arr[0]=UUU;
    Loop01_Para.use_arr[1]=VVV;
    Loop01_Para.use_arr[2]=FRI;
    Loop01_Para.use_arr[3]=DZU;
    Loop01_Para.use_arr[4]=DZT;
    Loop01_Para.use_arr[5]=VSHEAR;
    Loop01_Para.use_arr[6]=DBLOC;
    Loop01_Para.use_arr[7]=zgrid;
    Loop01_Para.use_arr[8]=RI_LOC;
    Loop01_Para.use_arr[9]=WORK0;

    Loop01_Para.use_arr[10]=AT0;
    Loop01_Para.use_arr[11]=ATS;
    Loop01_Para.use_arr[12]=ATW;
    Loop01_Para.use_arr[13]=ATSW;
}
   
 void ri_loop01_sw_()   
{
    athread_spawn(ri_loop01_fun,&Loop01_Para);
    athread_join();
}



void ri_loop02_fun_(int* nx_block, int* ny_block,int* bid,int*k,int*km,
                      double* p5,double* p25,int *KMT,
                      double* WORK0, double* FRI, double* RI_LOC)

{
    Loop02_Para.param[0]=*nx_block;
    Loop02_Para.param[1]=*ny_block;
    Loop02_Para.param[2]=*bid;
    Loop02_Para.param[3]=*k;
    Loop02_Para.param[4]=*km;
    Loop02_Para.param[5]=athread_get_max_threads();

    Loop02_Para.dparam[0]=*p5;
    Loop02_Para.dparam[1]=*p25;

    Loop02_Para.duse_arr[0]=KMT;
    
    Loop02_Para.use_arr[0]=WORK0;
    Loop02_Para.use_arr[1]=FRI;
    Loop02_Para.use_arr[2]=RI_LOC;
}
   
 void ri_loop02_sw_()   
{
    athread_spawn(ri_loop02_fun,&Loop02_Para);
    athread_join();
}


void ri_loop03_fun_(int* nx_block, int* ny_block,int* bid,int*km,int*lniw,int*lr,
                    double* c0,double* c1,double *Riinfty,double*rich_mix,
                    double* VDC, double*VISC, double*KVMIX,double *KVMIX_M,
                    double *bckgrnd_vdc,double*bckgrnd_vvc,double*WORK0,double*FRI)
{
  
    Loop03_Para.param[0]=*nx_block;
    Loop03_Para.param[1]=*ny_block;
    Loop03_Para.param[2]=*bid;
    Loop03_Para.param[3]=*km;
    Loop03_Para.param[4]=*lniw;
    Loop03_Para.param[5]=*lr;
    Loop03_Para.param[6]=athread_get_max_threads();

    Loop03_Para.dparam[0]=*c0;
    Loop03_Para.dparam[1]=*c1;
    Loop03_Para.dparam[2]=*Riinfty;
    Loop03_Para.dparam[3]=*rich_mix;   

    Loop03_Para.use_arr[0]=VDC;
    Loop03_Para.use_arr[1]=VISC;
    Loop03_Para.use_arr[2]=KVMIX;
    Loop03_Para.use_arr[3]=KVMIX_M;
    Loop03_Para.use_arr[4]=bckgrnd_vdc;
    Loop03_Para.use_arr[5]=bckgrnd_vvc;
    Loop03_Para.use_arr[6]=WORK0;
    Loop03_Para.use_arr[7]=FRI;
   
}
   
 void ri_loop03_sw_()   
{
    athread_spawn(ri_loop03_fun,&Loop03_Para);
    athread_join();
 }
                      
