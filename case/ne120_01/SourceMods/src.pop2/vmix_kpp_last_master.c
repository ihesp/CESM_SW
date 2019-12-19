#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include "vmix_kpp_last.h"

extern SLAVE_FUN(s_kpp_last1)();
extern SLAVE_FUN(s_kpp_last2)();


void m_kpp_last_sw(int *KBL, int *KMT, int *KMU,  
     double *dz, double *zt, double *zgrid, 
     double *AU0, double *AUN, double *AUE, double *AUNE, 
     double *BFSFC, double *HMXL,
     double *STF, double *GHAT, double *DBSFC, double *DBLOC, double *DZT,
     double *VISC, double *VVC, double *VDC, double *KPP_SRC, 
     int myproc, int nx_block, int ny_block, int km, int nt, int bid,
     double convect_visc, double convect_diff, double BVSQcon
#if defined(POPCHKJN) 
     , int steps
#endif
     )

 {
  double locdbl[10*80*60];
  
  struct param_kpp_last_s1 st2;

  int locint[3*80*60],*m2di;  
  double *m1dd, *m2dd;

  int nm,nxy,pij;
  int i,j,k,n,threads_num=64;

#ifdef POPINSJN 
  if(myproc==0 && steps<=2)
     printf("step into kpp_bldepth_master.\n"); 
#endif


  threads_num=athread_get_max_threads();

  st2.ipts[0] = nx_block;
  st2.ipts[1] = ny_block;
  st2.ipts[2] = km;
  st2.ipts[3] = nt;
  
  st2.ipts[4] = myproc;
  st2.ipts[5] = threads_num;


  st2.dpts[0] = convect_visc;
  st2.dpts[1] = convect_diff;
  st2.dpts[2] = BVSQcon;


//--------------------------------  
  nxy = nx_block*ny_block;

  m2di = NULL;
//  m2di = (int *)malloc(SIZEINT*3*nxy);  
  m2di = locint;  
  for(j=0;j<ny_block;j++)
  for(i=0;i<nx_block;i++)
     {
     pij = j*nx_block + i;
     m2di[pij      ] = KBL[pij];
     m2di[pij+  nxy] = KMT[pij];
     m2di[pij+2*nxy] = KMU[pij];
     } 
  st2.iarr[0]=m2di;
  

  m1dd = NULL;
//  m1dd = (double *)malloc(SIZEDBL*(3*km+2));  
  m1dd = locdbl;  
  for(i=0;i<km;i++)
     {
     m1dd[i     ] = dz[i];
     m1dd[i+  km] = zt[i];
     m1dd[i+2*km] = zgrid[i];
     } 
  m1dd[3*km] = zgrid[km];
  m1dd[3*km+1] = zgrid[km+1];
  st2.darr[0]=m1dd;
 
#if defined(POPTESTJN) 
  if(myproc==0 && steps<=2)
    {
    printf("master- nx: %d, km: %d, nt:%d. \n",nx_block,km,nt);
    printf("master- dz[0]: %f, dz[15]: %f. \n",m1dd[0],m1dd[15]);
    printf("master- zt[0]: %f, zt[15]: %f. \n",m1dd[km],m1dd[km+15]);
    printf("master- zgrid[0]: %f, zgrid[km+1]: %f. \n",m1dd[2*km],m1dd[3*km+1]);
    }
#endif
 

  m2dd = NULL;
//  m2dd = (double *)malloc(SIZEDBL*(6+nt)*nxy);  
  m2dd = locdbl + 4*km;  
  for(j=0;j<ny_block;j++)
  for(i=0;i<nx_block;i++)
     {
     pij = j*nx_block + i;
     m2dd[pij      ] = AU0[pij];
     m2dd[pij+  nxy] = AUN[pij];
     m2dd[pij+2*nxy] = AUE[pij];
     m2dd[pij+3*nxy] = AUNE[pij];

     m2dd[pij+4*nxy] = BFSFC[pij];
     m2dd[pij+5*nxy] = HMXL[pij];

     for(k=0;k<nt;k++)
        m2dd[pij+(6+k)*nxy] = STF[pij+k*nxy];
     } 
  st2.darr[1]=m2dd;


  st2.darr[2] = GHAT;
  st2.darr[3] = DBSFC;
  st2.darr[4] = DBLOC;
  st2.darr[5] = DZT;
  st2.darr[6] = VISC;
  st2.darr[7] = VVC;
  st2.darr[8] = VDC;
  st2.darr[9] = KPP_SRC;

  
// transfer the addr of struct st2!!!

  athread_spawn(s_kpp_last1, &st2);
  athread_join();

//  c_kpp_last1(&st2);



  athread_spawn(s_kpp_last2, &st2);
  athread_join();


//  c_kpp_last2(&st2);


  for(j=0;j<ny_block;j++)
  for(i=0;i<nx_block;i++)
     {
     pij = j*nx_block + i;
     HMXL[pij]  = m2dd[pij+5*nxy];
     }


//  free(m2di);
//  free(m1dd);
//  free(m2dd);

  return;
}

