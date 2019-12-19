#include<athread.h>
#include<stdio.h>
#include "step_mod_struct.h"

extern SLAVE_FUN(fun_state_step_mod)();
extern SLAVE_FUN(fun_state_baroclinic)();
extern SLAVE_FUN(fun2_state_baroclinic)();


void loop_step_mod_(int *nx_block, int *ny_block,int* km,int * oldtime,int* curtime,
                    double *TRACER,double *TRACER2,double *RHO,double* RHO2,double *tmax,double *tmin,double*smax,double*smin,double* pressz )
{
       double *pkg=NULL;
       int  k;
       struct kpp_ddmix s1;
       pkg=(double *)malloc(8*(5*(*km)));

       for(k=0;k < *km;k++)
       {
          pkg[k]       =pressz[k];
          pkg[*km+k]   =tmax[k];
          pkg[2* *km+k]=tmin[k];
          pkg[3* *km+k]=smax[k];
          pkg[4* *km+k]=smin[k];
       }
       
       s1.param[0]=*nx_block; 
       s1.param[1]=*ny_block;
       s1.param[2]=athread_get_max_threads();
       s1.param[3]=*km;
       s1.param[4]=*oldtime;
       s1.param[5]=*curtime;

       s1.arr[0]=TRACER;
       s1.arr[1]=RHO;
       s1.arr[2]=TRACER2;
       s1.arr[3]=RHO2;

       s1.std[0]=pkg;

       athread_spawn(fun_state_step_mod,&s1);
       athread_join();
       free(pkg);
}


void loop_state_baroclinic_(int *nx_block, int *ny_block,int* km,
                    double *TRACER,double *RHO,double *tmax,double *tmin,double*smax,double*smin,double* pressz )
{
    
       double *pkg=NULL;
       int  k;

       pkg=(double *)malloc(8*(5*(*km)));

       for(k=0;k < *km;k++)
       {
          pkg[k]       =pressz[k];
          pkg[*km+k]   =tmax[k];
          pkg[2* *km+k]=tmin[k];
          pkg[3* *km+k]=smax[k];
          pkg[4* *km+k]=smin[k];
       }

       struct kpp_ddmix s1;
       s1.param[0]=*nx_block; 
       s1.param[1]=*ny_block;
       s1.param[2]=athread_get_max_threads();
       s1.param[3]=*km;
       s1.arr[0]=TRACER;
       s1.arr[1]=RHO;
       s1.std[0]=pkg;

       athread_spawn(fun_state_baroclinic,&s1);
       athread_join();
       free(pkg);
}

void loop2_state_baroclinic_(int *nx_block, int *ny_block,int* km,
                    double *TRACER,double *TRACER2,double *RHO,double *tmax,double *tmin,double*smax,double*smin,double* pressz )
{
    

       struct kpp_ddmix s1;
       s1.param[0]=*nx_block; 
       s1.param[1]=*ny_block;
       s1.param[2]=athread_get_max_threads();
       s1.param[3]=*km;

       s1.dparam[0]=*pressz;
       s1.dparam[1]=*tmax;
       s1.dparam[2]=*tmin;
       s1.dparam[3]=*smax;
       s1.dparam[4]=*smin;


       s1.arr[0]=TRACER;
       s1.arr[1]=RHO;
       s1.arr[2]=TRACER2;

       athread_spawn(fun2_state_baroclinic,&s1);
       athread_join();
}

