#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include "clinic_tavg_accum.h"

#define locmin(a,b)  ((a)<(b)?(a):(b))

extern SLAVE_FUN(s_clinic_tavg)();


void m_clinic_tavg_sw(
#ifdef TAVG_R8
    double *TAVG_BUF_2D,double *TAVG_BUF_3D, 
#else
    float *TAVG_BUF_2D,float *TAVG_BUF_3D,
#endif
    double *UVEL,double *VVEL, 
    double *TRACER1cur,double *TRACER1old,double *TRACER2cur,
    double *RHO,double *PSURF,double *DZT,double *DH, 
    double dtavg,double *dzkm,double grav,double hflux_factor,double salinity_factor, 
    int *KMT,int *intCALCT, 
    int *tavgids,int *tavgflags,int *tavglocs,int *tavgdims,int *tavgmeths, 
    int myproc,int km,int nx,int ny,int iblock,int niblock,
    int sfc_layer_type,int sfc_layer_varthick,int numtavg,int part_bott_cells
     )

  {
  struct param_tavg_accum st2;

  int i, j, k,ktmp,kdiv,kmod, threads_num=64;

  threads_num = athread_get_max_threads();

  st2.ipts[0] = nx;
  st2.ipts[1] = ny;
  st2.ipts[2] = sfc_layer_type;
  st2.ipts[3] = sfc_layer_varthick;
  st2.ipts[4] = myproc;
  st2.ipts[6] = km;
  st2.ipts[7] = iblock;
  st2.ipts[8] = niblock; 
  st2.ipts[9] = part_bott_cells;  
  st2.ipts[10] = numtavg;  
  st2.ipts[11] = threads_num;  
  
  st2.iarr[0] = KMT;
  st2.iarr[1] = intCALCT;  

  st2.dpts[0] = grav;
  st2.dpts[1] = hflux_factor;
  st2.dpts[2] = salinity_factor;  
  st2.dpts[3] = dtavg;  
  for(k=0;k<km;k++)
     st2.dpts[4+k] = dzkm[k];

#ifdef TAVG_R8
  st2.darr[0] = TAVG_BUF_2D;
  st2.darr[1] = TAVG_BUF_3D;
#else
  st2.farr[0] = TAVG_BUF_2D;
  st2.farr[1] = TAVG_BUF_3D;
#endif
  st2.darr[2] = UVEL;  
  st2.darr[3] = VVEL;  
  st2.darr[4] = TRACER1cur;  
  st2.darr[5] = TRACER1old;  
  st2.darr[6] = TRACER2cur; 
  st2.darr[7] = RHO; 
  st2.darr[8] = PSURF; 
  st2.darr[9] = DZT;   
  st2.darr[10] = DH;  
  

  for(i=0;i<numtavg;i++)
     {
     st2.tavgindx[i*12+0] = tavgids[i]-1;
     st2.tavgindx[i*12+2] = tavglocs[i];
     st2.tavgindx[i*12+3] = tavgdims[i];	  
     st2.tavgindx[i*12+4] = tavgmeths[i];	  

     st2.tavgindx[i*12+8] = 0;	  
     st2.tavgindx[i*12+9] = 0;	  
     st2.tavgindx[i*12+10] = 0;	  
     st2.tavgindx[i*12+11] = 0;	  
     for(k=0;k<km;k++)
        {
        kdiv = k/16;
        kmod = k%16;
        ktmp = i*12+(8+kdiv);
        st2.tavgindx[ktmp] = st2.tavgindx[ktmp] + (int)(tavgflags[i*km+k]<<kmod);      
        }
     }
 
  athread_spawn(s_clinic_tavg, &st2);
  athread_join();

  return;
  }


