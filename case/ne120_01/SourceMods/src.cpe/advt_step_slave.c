#include<slave.h>
#include<stdio.h>
#include<math.h>
#include"advt_tavg_accum.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))


void advt_state_s(struct param_tavg_accum *s1)
{

      double  mwjfnums0t0,  mwjfnums0t1,  mwjfnums0t2, mwjfnums0t3,
              mwjfnums1t0,  mwjfnums1t1,  mwjfnums2t0,
              mwjfdens0t0,  mwjfdens0t1,  mwjfdens0t2,  mwjfdens0t3, mwjfdens0t4,
              mwjfdens1t0,  mwjfdens1t1,  mwjfdens1t3,
              mwjfdensqt0,  mwjfdensqt2;
      const  double
      mwjfdp0s0t0 =   1.0e+0,
      mwjfdp0s0t1 =   7.28606739e-3,
      mwjfdp0s0t2 =  -4.60835542e-5,
      mwjfdp0s0t3 =   3.68390573e-7,
      mwjfdp0s0t4 =   1.80809186e-10,
      mwjfdp0s1t0 =   2.14691708e-3,
      mwjfdp0s1t1 =  -9.27062484e-6,
      mwjfdp0s1t3 =  -1.78343643e-10,
      mwjfdp0sqt0 =   4.76534122e-6,
      mwjfdp0sqt2 =   1.63410736e-9,
      mwjfdp1s0t0 =   5.30848875e-6,
      mwjfdp2s0t3 =  -3.03175128e-16,
      mwjfdp3s0t1 =  -1.27934137e-17;
      const  double
      mwjfnp0s0t0 =   9.99843699e-1,
      mwjfnp0s0t1 =   7.35212840e-3,
      mwjfnp0s0t2 =  -5.45928211e-5,
      mwjfnp0s0t3 =   3.98476704e-7,
      mwjfnp0s1t0 =   2.96938239e-3,
      mwjfnp0s1t1 =  -7.23268813e-6,
      mwjfnp0s2t0 =   2.12382341e-6,
      mwjfnp1s0t0 =   1.04004591e-5,
      mwjfnp1s0t2 =   1.03970529e-10,
      mwjfnp1s1t0 =   5.18761880e-9,
      mwjfnp2s0t0 =  -3.24041825e-11,
      mwjfnp2s0t2 =  -1.23869360e-14;
      const double
      c1=1.0,
      c10=10.0,
      c1000=1000.0;

     int myid,nx,ny,km,tn,k,quotient,i,quotient1,task_size1,koffset;
     int usize,usize1,usize3,usize_tracer,max_usize,j_size,i_size,task_offset,task_num,task_num1,task,task_size,tail_size;
     int tile,task_tracer,rho_koffset,rho_toffset,task_rho;
     double tmax,tmin,smax,smin,pressz;
     double *p_TRACER=NULL,*p_RHO=NULL;
     double p;
     double TEMPK,SALTK,TQ,SQ,SQR,ZWORK1,ZWORK2,DENOMK; 


//     nx       = s1->ipts[0];
//     ny       = s1->ipts[1];
//     tn       = s1->ipts[14];
//     km       = s1->ipts[3];
     task_size1 = s1->ipts[12];
    
     p_TRACER   = s1->darr[6];
     p_RHO      = s1->darr[8];

     pressz = s1->darr[15][0];
     tmax   = s1->darr[11][0];
     tmin   = s1->darr[12][0];
     smax   = s1->darr[13][0];
     smin   = s1->darr[14][0]; 

    
     mwjfnums0t1 = mwjfnp0s0t1;
     mwjfnums0t3 = mwjfnp0s0t3;
     mwjfnums1t1 = mwjfnp0s1t1;
     mwjfnums2t0 = mwjfnp0s2t0;

     mwjfdens0t2 = mwjfdp0s0t2;
     mwjfdens0t4 = mwjfdp0s0t4;
     mwjfdens1t0 = mwjfdp0s1t0;
     mwjfdens1t1 = mwjfdp0s1t1;
     mwjfdens1t3 = mwjfdp0s1t3;
     mwjfdensqt0 = mwjfdp0sqt0;
     mwjfdensqt2 = mwjfdp0sqt2;
        
     p=c10*pressz;

     mwjfnums0t0 = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0);
     mwjfnums0t2 = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2);
     mwjfnums1t0 = mwjfnp0s1t0 + p*mwjfnp1s1t0;

     mwjfdens0t0 = mwjfdp0s0t0 + p*mwjfdp1s0t0;
     mwjfdens0t1 = mwjfdp0s0t1 + p*p*p* mwjfdp3s0t1;
     mwjfdens0t3 = mwjfdp0s0t3 + p*p*   mwjfdp2s0t3;

  
     for(i=0;i<task_size1;i++)
       {
       koffset = task_size1+i;
                
       TEMPK = p_TRACER[i];
       if(TEMPK>tmax) TQ=tmax;
       else if(TEMPK<tmin) TQ=tmin;
       else TQ=TEMPK;
        	   
       SALTK = p_TRACER[koffset];         // s with j=1
       if(SALTK>smax) SQ=smax;
       else if(SALTK<smin) SQ=smin;
       else SQ=SALTK;
 
       SQ  = c1000*SQ;
       SQR =(double)sqrt((double)SQ);

       ZWORK1 = mwjfnums0t0   + TQ * (mwjfnums0t1      + TQ * (mwjfnums0t2 +
                mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0 +
                mwjfnums1t1 * TQ   +       mwjfnums2t0      * SQ);
                   
       ZWORK2 = mwjfdens0t0 + TQ * (mwjfdens0t1 + TQ * (mwjfdens0t2 +
                TQ * (mwjfdens0t3 +   mwjfdens0t4 * TQ))) +
                SQ * (mwjfdens1t0  + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
                SQR* (mwjfdensqt0      + TQ*TQ*mwjfdensqt2));

       DENOMK = c1/ZWORK2;
                  
       p_RHO[i] = ZWORK1*DENOMK;
       }   

return;
}




