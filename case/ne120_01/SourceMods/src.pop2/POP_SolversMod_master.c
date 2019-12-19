#include <athread.h>
#include<math.h>
#include "POP_SolversMod_struct.h"
#include<stdio.h> 
extern SLAVE_FUN(fun_pop_solvers_operator)(),SLAVE_FUN(fun_pop_solvers_preconditioner)(),SLAVE_FUN(fun_pop_solvers_step5)(),
       SLAVE_FUN(fun_pop_solvers_step5_2)();


void loop_pop_solvers_operator_(double *btropWgtCenter,double *btropWgtNorth,
		   double *btropWgtEast,double *btropWgtNE,double *AX,double *X,int *nx_block,int* ny_block)
  {

   struct param_zy s1;
   s1.param[0]=*nx_block;
   s1.param[1]=*ny_block;

   s1.use_arr[0]=btropWgtCenter;
   s1.use_arr[1]=btropWgtNorth;
   s1.use_arr[2]=btropWgtEast;
   s1.use_arr[3]=btropWgtNE;
   s1.use_arr[4]=AX;
   s1.use_arr[5]=X;

   athread_spawn(fun_pop_solvers_operator,&s1);
   athread_join();
  }

void loop_pop_solvers_step5_(int* nx,int *ny,double *csomga,double *csy,double *Q,double *R,double *X,double *S,double *B, 
                         double *btropWgtCenter,double *btropWgtNorth,double *btropWgtEast,double *btropWgtNE)
{
   struct param_zy s1;
   s1.param[0]=*nx;
   s1.param[1]=*ny;

   s1.dparam[0]=*csomga;
   s1.dparam[1]=*csy;
   s1.use_arr[0] = Q;
   s1.use_arr[1] = R;
   s1.use_arr[2] = X;
   s1.use_arr[3] = S; 
   s1.use_arr[4]=btropWgtCenter;
   s1.use_arr[5]=btropWgtNorth;
   s1.use_arr[6]=btropWgtEast;
   s1.use_arr[7]=btropWgtNE;
   s1.use_arr[8] = B;

   athread_spawn(fun_pop_solvers_step5,&s1);
   athread_join();
  }

void loop_pop_solvers_step5_2_(int* nx,int *ny,double *csomga,double *csy,double *Q,double *R,double *X,double *S,double *B, 
                         double *btropWgtCenter,double *btropWgtNorth,double *btropWgtEast,double *btropWgtNE,double* work0)
{
   struct param_zy s1;
   s1.param[0]=*nx;
   s1.param[1]=*ny;

   s1.dparam[0]=*csomga;
   s1.dparam[1]=*csy;
   s1.use_arr[0] = Q;
   s1.use_arr[1] = R;
   s1.use_arr[2] = X;
   s1.use_arr[3] = S; 
   s1.use_arr[4]=btropWgtCenter;
   s1.use_arr[5]=btropWgtNorth;
   s1.use_arr[6]=btropWgtEast;
   s1.use_arr[7]=btropWgtNE;
   s1.use_arr[8] = B;
   s1.use_arr[9] = work0;

   athread_spawn(fun_pop_solvers_step5_2,&s1);
   athread_join();
  }


    
void loop_preconditioner_(int* nx,int* ny,int* EvpYnb,int* EvpXnb,int* EvpYbidx,int* EvpXbidx,int* landIndx,double* R,double* InvEvpCenterWgt,double* EvpCenterWgt,double* EvpNeWgt,double* EvpRinv,double* InvEvpNeWgt,int* m)  
{
    struct param_zy s1;

    s1.param[0]=*nx;
    s1.param[1]=*ny;
    s1.param[2]=*EvpYnb;
    s1.param[3]=*EvpXnb;
    s1.param[4]=*m;
    

    s1.use_arr_int[0]=EvpYbidx;
    s1.use_arr_int[1]=EvpXbidx;
    s1.use_arr_int[2]=landIndx;
    s1.use_arr[0]=R;
    s1.use_arr[1]=InvEvpCenterWgt;
    s1.use_arr[2]=EvpCenterWgt;
    s1.use_arr[3]=EvpNeWgt;
    s1.use_arr[4]=EvpRinv;
    s1.use_arr[5]=InvEvpNeWgt;
    
    athread_spawn(fun_pop_solvers_preconditioner,&s1);
    athread_join();
}

