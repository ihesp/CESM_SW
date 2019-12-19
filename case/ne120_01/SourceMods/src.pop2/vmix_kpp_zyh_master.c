#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include <math.h>
#include <simd.h>
#include "vmix_kpp_zyh_struct.h"


extern SLAVE_FUN(kpp_ddmix_loop1)();

extern SLAVE_FUN(kpp_state_loop2)();

extern SLAVE_FUN(s_kpp_bldepth_loop2)();

void m_kpp_ddmix_loop1(double *TRCR, double *VDC, int myproc, int nx_block, int ny_block, int km, int nt, 
     int knxt, int kup, double *tmax, double *tmin, double *smax, double *smin, double *pressz, 
     double grav, double Rrho0, double dsfmax
#if defined(POPCHKJN) || defined(POPTESTJN) 
     , int steps
#endif
#if defined(POPCHKJN)
     , double *swork11, double *swork21, double *swork1k
     , double *swork2k, double *swork1km, double *swork2km
#endif
     )

 {
 
  int k=0;
  int i,j,cid,threads_num;

#ifdef POPTESTJN 
  FILE *svar;
#endif

  double *pts;

  struct param_kpp_state_s2 st2;

#ifdef POPTESTJN 
  double start_tm,stop_tm;
#endif



  threads_num=athread_get_max_threads();

  pts=(double *)malloc(sizeof(double)*(5*km+3));

       
  st2.nblock[0]=nx_block;
  st2.nblock[1]=ny_block;
  st2.nblock[2]=km;
  st2.nblock[3]=nt;
  st2.nblock[4]=knxt;
  st2.nblock[5]=kup;
  
  st2.core_info[0]=myproc;
  st2.core_info[1]=threads_num;


  st2.local_arr[0]=TRCR;
  st2.local_arr[1]=VDC;


  for(k=0;k<km;k++)
    {
     pts[k]=pressz[k];
     pts[km+k]=tmax[k];
     pts[2*km+k]=tmin[k];
     pts[3*km+k]=smax[k];
     pts[4*km+k]=smin[k];
    }
  pts[5*km]=grav;
  pts[5*km+1]=Rrho0;
  pts[5*km+2]=dsfmax;
  st2.local_pts[0]=pts;

#ifdef POPCHKJN 
  st2.sdiag[0]=swork11;
  st2.sdiag[1]=swork1k;
  st2.sdiag[2]=swork1km;
  st2.sdiag[3]=swork21;
  st2.sdiag[4]=swork2k;
  st2.sdiag[5]=swork2km;
#endif

#ifdef POPTESTJN 
  if(myproc==0 && steps<=2)
    {
    printf("nx_block_master: %d \n",nx_block);
    printf("grav_master: %lf; \n",grav);
    }
#endif


#ifdef POPTESTJN 
  getlooptime(&start_tm);
#endif

  athread_spawn(kpp_ddmix_loop1,&st2);
  athread_join();

#ifdef POPTESTJN 
  getlooptime(&stop_tm);
  if(myproc==0 && steps<=2)
    {
    printf("Computing VDC costs %lf (us); \n",stop_tm-start_tm);
    }
#endif

  free(pts);

  return;
}




void m_kpp_state_loop2(double *TRCR, int *KMT, int myproc, int nx_block, int ny_block, int km, int nt,
     int klvl, int kprev, double *tmax, double *tmin, double *smax, double *smin, double *pressz, 
     double grav, double *DBLOC, double *DBSFC
#if defined(POPCHKJN) || defined(POPTIMEJN) 
     , int steps
#endif
#if defined(POPCHKJN)
     , double *swork11, double *swork21, double *swork1k
     , double *swork2k, double *swork1km
#endif
     )

 {
 
  int k=0;
  int i,j,cid,lnx,lny,lkm,lnt,threads_num;
  int idxn,idxo;

#ifdef POPTESTJN 
  FILE *svar;
#endif

  double *pts;

  struct param_kpp_state_s2 st2;

#ifdef POPTIMEJN 
  double start_tm,stop_tm;
#endif



  threads_num=athread_get_max_threads();

  lnx=nx_block;
  lny=ny_block;
  lkm=km;
  lnt=nt;

  pts=(double *)malloc(sizeof(double)*(5*lkm+1));

       
  st2.nblock[0]=nx_block;
  st2.nblock[1]=ny_block;
  st2.nblock[2]=km;
  st2.nblock[3]=nt;
  st2.nblock[4]=klvl;
  st2.nblock[5]=kprev;
  
  st2.core_info[0]=myproc;
  st2.core_info[1]=threads_num;

  st2.local_iarr[0]=KMT;

  st2.local_arr[0]=TRCR;
  st2.local_arr[1]=DBLOC;
  st2.local_arr[2]=DBSFC;
  


  for(k=0;k<lkm;k++)
    {
     pts[k]=pressz[k];
     pts[lkm+k]=tmax[k];
     pts[2*lkm+k]=tmin[k];
     pts[3*lkm+k]=smax[k];
     pts[4*lkm+k]=smin[k];
    }
  pts[5*lkm]=grav;
  st2.local_pts[0]=pts;

#ifdef POPCHKJN 
  st2.sdiag[0]=swork11;
  st2.sdiag[1]=swork1k;
  st2.sdiag[2]=swork1km;
  st2.sdiag[3]=swork21;
  st2.sdiag[4]=swork2k;
#endif


#ifdef POPTIMEJN 
  getlooptime(&start_tm);
#endif

  athread_spawn(kpp_state_loop2,&st2);
  athread_join();

#ifdef POPTIMEJN 
  getlooptime(&stop_tm);
  if(myproc==0 && steps<=2)
    {
    printf("Computing LOC and SFC costs %lf (us); \n",stop_tm-start_tm);
    }
#endif


/*
 if(myproc==0 && steps<=2)
  {
  svar=fopen("WORK_slave.txt","a");
  for(k=1;k<lkm;k++)
     for(j=0;j<lny;j++)
        for(i=0;i<lnx;i++)
           {
//           fprintf(svar,"c=%1d j=%-3d i=%-2d k=%-2d %20.16lf %20.16lf %20.16lf \n", steps,j,i,k,swork11[j*lnx*lkm+i*lkm+k], swork1k[j*lnx*lkm+i*lkm+k], swork1km[j*lnx*lkm+i*lkm+k]);    //left alligned
           idxn=j*lnx*lkm+i*lkm+k;
           fprintf(svar,"c=%1d j=%-3d i=%-2d k=%-2d %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf \n", steps,j,i,k,swork11[idxn], swork1k[idxn], swork1km[idxn], swork21[idxn], swork2k[idxn]);
//           fprintf(svar,"c=%1d j=%3d i=%2d k=%2d %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf \n", steps,j,i,k,swork11[idxn], swork1k[idxn], swork1km[idxn], swork21[idxn], swork2k[idxn]);   //right alligned
           }
  fclose(svar);
  }
*/

  free(pts);

  return;
}




void m_kpp_bldepth_loop2(double *TRCR, double *DBLOC, double *DBSFC, 
     double *UUU, double *VVV, double *UCUR, double *VCUR, 
     double *HBLT, double *USTAR, double *BFSFC, double *STABLE,
     int *KMT, int *KBL, 
     double *DZT, double *DZU, double *niuel, double *nivel, 
     double *Tr, int *CHLINDX, double *BO, double *BOSOL, double *ZKL,
     int myproc, int nx_block, int ny_block, int km, int nt, int bid,
     int loc_short_wave, int loc_partial_bottom_cells, int loc_lniw_mixing, int loc_linertial,
     double *zgrid, double *Ricr,
     double Vtc, double zeta_s, double a_s, double c_s, double vonkar, 
     double eps, double eps2, double epssfc, double concv
#if defined(POPCHKJN) || defined(POPTIMEJN) 
     , int steps
#endif
#if defined(POPINSJN)
     , double *swork11, double *swork21, double *swork1k
     , double *swork2k, double *swork1km, double *swork2km
#endif
     )

 {
 
  int k=0,ksol;
  int i,j,cid,threads_num=64;
  int idxn,idxo;

//  double mTrans[SZswtr], m1din[SZ1ddbl];  
//  double m2din[SZ2ddbl];  
  double *mTrans, *m1din;  
  double *m1dout, *m2din, *m2dout;
  
  struct param_kpp_bldepth_s1 st2;

#ifdef POPINSJN 
  FILE *svar;
#endif

#ifdef POPTIMEJN 
  double start_tm,stop_tm;
#endif


#ifdef POPINSJN 
  if(myproc==0 && steps<=2)
     printf("step into kpp_bldepth_master.\n"); 
#endif


  threads_num=athread_get_max_threads();

  st2.nblock[0]=nx_block;
  st2.nblock[1]=ny_block;
  st2.nblock[2]=km;
  st2.nblock[3]=nt;
  
  st2.core_info[0]=myproc;
  st2.core_info[1]=threads_num;

  st2.local_iarr[0] = KMT;
//  st2.local_iarr[0] = KMT + bid*ny_block*nx_block;
  st2.local_iarr[1] = KBL;
#ifdef POPINSJN 
  if(myproc==0 && steps<=3)
     printf("KMT in kpp_bldepth_master: %ld\n", KMT); 
#endif


//--------------------------------
  st2.global_logic[0]=loc_short_wave;
  st2.global_logic[1]=loc_partial_bottom_cells;
  st2.global_logic[2]=loc_lniw_mixing;
  st2.global_logic[3]=loc_linertial;




//--------------------------------
  st2.local_arr[0] = TRCR;
  st2.local_arr[1] = DBLOC;
  st2.local_arr[2] = DBSFC;
  st2.local_arr[3] = UUU;
  st2.local_arr[4] = VVV;
#ifdef POPINSJN 
  if(myproc==0 && steps<=3)
     printf("UUU in kpp_bldepth_master: %ld\n", UUU); 
#endif


////  STF;
////  SHF_QSW;
////  STABLE;
////  SMFT;

////  BOLUS_SP;

  
  
//--------------------------------  
  st2.global_arr[0]=DZT;
  st2.global_arr[1]=DZU;
//  st2.global_arr[0]=DZT + bid*(km+2)*ny_block*nx_block;
//  st2.global_arr[1]=DZU + bid*(km+2)*ny_block*nx_block;
#ifdef POPINSJN 
  if(myproc==0 && steps<=3)
     printf("DZT in kpp_bldepth_master: %ld\n", DZT); 
#endif

  

#ifdef POPINSJN 
  if(myproc==0 && steps<=2)
     printf("malloc_1din in kpp_bldepth_master.\n"); 
#endif




//--------------------------------  

  m1din = NULL;
  m1din = (double *)malloc(8*(2*km+2));  
#ifdef POPINSJN 
  if(m1din==NULL)
    {
    if(myproc==0 && steps<=3)
       printf("malloc m1din of kpp_bldep_master failed! \n");
    return;
    }
#endif

  for(i=0;i<km;i++)
     {
     m1din[i] = Ricr[i];
     m1din[i+km] = zgrid[i];
     } 
  m1din[2*km] = zgrid[km];
  m1din[2*km+1] = zgrid[km+1];
  st2.global_arr[2]=m1din;
#ifdef POPINSJN 
  if(myproc==0 && steps<=3)
     printf("m1din in kpp_bldepth_master: %ld\n", m1din); 
#endif

  
  
#ifdef POPINSJN 
  if(myproc==0 && steps<=2)
     printf("malloc_mTrans in kpp_bldepth_master:\n"); 
#endif

//--------------------------------  
//  st2.global_arr[6]=Tr;

// computing mTrans with Tr(0:2*km,0:400) with ksol=2*km, nsub=400,
// and CHLINDX(nx_block,ny_block,max_bclinic):  

  mTrans = NULL;
  mTrans = (double *)malloc(8*km*nx_block*ny_block);
#ifdef POPINSJN 
  if(mTrans==NULL)
    {
    if(myproc==0 && steps<=3)
       printf("malloc mTrans of kpp_bldep_master failed! \n");
    return;
    }
#endif

  ksol = 2*km;
  for(j=0;j<ny_block;j++)
     for(i=0;i<nx_block;i++)
        {
        idxn = CHLINDX[j*nx_block+i];
//        idxn = CHLINDX[j*nx_block+i + bid*ny_block*nx_block];
        mTrans[(0*ny_block+j)*nx_block+i] = Tr[idxn*(ksol+1)+1];
        for(k=1;k<km;k++)
           mTrans[(k*ny_block+j)*nx_block+i] = Tr[idxn*(ksol+1)+(2*(k+1)-1)];  
        }  
  st2.global_arr[3]=mTrans;

#ifdef POPINSJN 
  if(myproc==0 && steps<=2)
     printf("malloc_mTrans in kpp_bldepth_master.\n"); 
#endif




  
#ifdef POPINSJN 
  if(myproc==0 && steps<=2)
     printf("malloc_2din in kpp_bldepth_master:\n"); 
#endif

  m2din = NULL;
  m2din = (double *)malloc(8*7*nx_block*ny_block);  
#ifdef POPINSJN 
  if(m2din==NULL)
    {
    if(myproc==0 && steps<=3)
       printf("malloc m2din of kpp_bldep_master failed! \n");
    return;
    }
#endif


  for(j=0;j<ny_block;j++)
  for(i=0;i<nx_block;i++)
     {
     m2din[j*nx_block+i] = USTAR[j*nx_block+i];
     m2din[(1*ny_block+j)*nx_block+i] = BO[j*nx_block+i];
     m2din[(2*ny_block+j)*nx_block+i] = BOSOL[j*nx_block+i];
     m2din[(3*ny_block+j)*nx_block+i] = UCUR[j*nx_block+i];	 
     m2din[(4*ny_block+j)*nx_block+i] = VCUR[j*nx_block+i];	 	 
/*
     niuel[j*nx_block+i];
     nivel[j*nx_block+i];
     BOLUS_SP[j*nx_block+i]; // reserved!
*/
     } 
  st2.local_arr[5]=m2din;
#ifdef POPINSJN 
  if(myproc==0 && steps<=3)
     printf("m2din in kpp_bldepth_master: %ld\n", m2din); 
#endif


  
#ifdef POPINSJN 
  if(myproc==0 && steps<=2)
     printf("malloc_2din in kpp_bldepth_master.\n"); 
#endif



  
#ifdef POPINSJN 
  if(myproc==0 && steps<=4)
     printf("pack others into struct in kpp_bldepth_master.\n"); 
#endif

//--------------------------------  
  st2.local_pts[0]=Vtc;
  st2.local_pts[1]=a_s;
  st2.local_pts[2]=c_s;
  st2.local_pts[3]=vonkar; 
  st2.local_pts[4]=eps;
  st2.local_pts[5]=eps2;
  st2.local_pts[6]=epssfc;
  st2.local_pts[7]=concv;
  st2.local_pts[8]=zeta_s;

#ifdef POPINSJN 
  st2.sdiag[0]=swork11;
  st2.sdiag[1]=swork1k;
  st2.sdiag[2]=swork1km;
  st2.sdiag[3]=swork21;
  st2.sdiag[4]=swork2k;
  st2.sdiag[5]=swork2km;
#endif


#ifdef POPINSJN 
  if(myproc==0 && steps<=4)
     printf("call s_kpp_bldepth from kpp_master:\n"); 
#endif


#ifdef POPINSJN 
  getlooptime(&start_tm);
#endif

// transfer the addr of struct st2!!!
//  s_kpp_bldepth_loop2(&st2);

  athread_spawn(s_kpp_bldepth_loop2, &st2);
  athread_join();

#ifdef POPINSJN 
  getlooptime(&stop_tm);
  if(myproc==0 && steps<=16)
    {
    printf("Computing bldepth slave costs %lf (us); \n",stop_tm-start_tm);
    }
#endif

  for(j=0;j<ny_block;j++)
  for(i=0;i<nx_block;i++)
     {
     HBLT[j*nx_block+i]  = m2din[(5*ny_block+j)*nx_block+i];
     BFSFC[j*nx_block+i] = m2din[(6*ny_block+j)*nx_block+i];
     }


#ifdef POPINSJN 
  if(myproc==0 && steps<=2)
     printf("call s_kpp_bldepth back to kpp_master.\n"); 
#endif


  free(mTrans);

  free(m1din);
  free(m2din);

  return;
}

