#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include "advt_tavg_accum.h"

#define locmin(a,b)  ((a)<(b)?(a):(b))

extern SLAVE_FUN(s_advt_tavg)();


void m_advt_tavg_sw(
#ifdef TAVG_R8
    double *TAVG_BUF_2D,double *TAVG_BUF_3D, 
#else
    float *TAVG_BUF_2D,float *TAVG_BUF_3D,
#endif
    double *VTN,double *UTE,double *WTK,
    double *TR_E,double *TR_N,double *TRCR,
    double *DZT,double *TAREA_R, 
    double dtavg,
    double *tmax,double *tmin,double *smax,double *smin,double *pressz, 
    int *tavgids,int *tavgflags,int *tavglocs,int *tavgdims,int *tavgmeths, 
    int* tadvtype,int *tavgintpara
     )

  {
  struct param_tavg_accum st2;

  int i, j, nt, numtavg, threads_num=64;

  threads_num = athread_get_max_threads();

  st2.ipts[0] = tavgintpara[0];   //nx_block
  st2.ipts[1] = tavgintpara[1];   //ny_block
  st2.ipts[2] = tavgintpara[15];  //my_task
  st2.ipts[3] = tavgintpara[2];   //km

  st2.ipts[6] = 100*tavgintpara[5]+tavgintpara[4];  //100*bid+nblocks

  numtavg = tavgintpara[6];   //numtavg
  st2.ipts[7] = numtavg;
  nt = tavgintpara[8];   //nt
  st2.ipts[4] = nt;

  st2.ipts[8] = tavgintpara[7];   //part_bott_cells
  st2.ipts[9] = tavgintpara[9];   //tadvect_centered
  for(i=0;i<nt;i++)
     st2.ipts[10+i] = tadvtype[i];  
  st2.ipts[15] = tavgintpara[3];  //k


  
  st2.dpts[0] = dtavg;  



#ifdef TAVG_R8
  st2.darr[0] = TAVG_BUF_2D;
  st2.darr[1] = TAVG_BUF_3D;
#else
  st2.farr[0] = TAVG_BUF_2D;
  st2.farr[1] = TAVG_BUF_3D;
#endif
  st2.darr[2] = VTN;  
  st2.darr[3] = UTE;  
  st2.darr[4] = WTK;  
  st2.darr[5] = TR_E;  
  st2.darr[6] = TR_N; 
  st2.darr[7] = TRCR; 
  st2.darr[9] = DZT;   
  st2.darr[10] = TAREA_R;  
  
  st2.darr[11] = tmax;   
  st2.darr[12] = tmin;   
  st2.darr[13] = smax;   
  st2.darr[14] = smin;   
  st2.darr[15] = pressz;   


  for(i=0;i<numtavg;i++)
     {
     st2.tavgindx[i*12+0] = tavgids[i]-1;
     st2.tavgindx[i*12+1] = tavgflags[i];      
     st2.tavgindx[i*12+2] = tavglocs[i];
     st2.tavgindx[i*12+3] = tavgdims[i];	  
     st2.tavgindx[i*12+4] = tavgmeths[i];	  
     }
 
  athread_spawn(s_advt_tavg, &st2);
  athread_join();

  return;
  }
