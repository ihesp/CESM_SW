#include <stdio.h>
#include <math.h>
#include <slave.h>
#include "vmix_kpp_zyh_struct.h"

#ifdef POPTESTJN
#include "cpe_print.h"
#endif

#define SZldmint 128
#define SZldmdbl 5600
//#define SZldmdbl 6000
//#define SZldmdbl 2600
#define locmax(a,b)  ((a)>(b)?(a):(b))


#define get_myrid(row)  \
asm volatile("rcsr   %0, 1" : "=r"(row))
#define get_mycid(col)  \
asm volatile("rcsr   %0, 2" : "=r"(col))
#define sync_mem   \
asm volatile("memb")
//asm volatile("memb\n\t":::"memory")


double *ldm_malloc(); 

#include "math_data.h"
//------------------------------------------------------------------
//void slave_kpp_ddmix_loop1(struct param_kpp_state_s2 *st2_slave)
void kpp_ddmix_loop1(struct param_kpp_state_s2 *st2_slave)
{
	
// BOP:
 volatile int get_reply, get_reply2, put_reply;
#ifdef POPCHKJN
 volatile int num_calls=0;
#endif

// Common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,task_pos,task_size,task_tile,task_all,last_size;
 int nx_block,ny_block,km,nt,lx_block;

 int offset,offset1,offset2,offset3;
 int offsetg,offsetg1,offsetg2,offsetg3;
 int offsetp,offsetp1,offsetp2,offsetp3;
 int max_ldm,mx_usize,usize,usize1,usize2,usize3;
 int i,j,k,n,m,l,ni,nj,nk,nm,nitr;

 double *a_buf,*b_buf,*c_buf,*e_buf,*f_buf,*g_buf;
 double *a0_buf,*b0_buf,*c0_buf,*e0_buf,*f0_buf,*g0_buf;

 int *d_buf,*d0_buf;

 int itmp,i1tmp,jtmp,j1tmp,ktmp,k1tmp;
 double etmp,stmp,tmpd,tmpres;


// Input variables: 
 double *TRCR;
 int knxt0,kup0;


// Input & Output variables:
 double *VDC;


// Temp variables:
 int knxt,kup;
 double 
      TQ, SQ,             // adjusted T,S
      TEMPK, SALTK, SQR,     // temporary variables
      p,                  // temporary pressure scalar
      DENOMK,WORK1,WORK2,WORK3,WORK4,
#ifdef POPCHKJN
      *mwork11,*mwork1k,*mwork1km,
      *mwork21,*mwork2k,*mwork2km,
      *swork11,*swork1k,*swork1km,
      *swork21,*swork2k,*swork2km,
#endif
      tmpw1,tmpw2,tmprho,
      *TEMPSFC,*TEMPKLP,
      *ALPHADT,*BETADS,PRANDTL,DIFFDD,
      *RRHO,*TALPHA,*SBETA,
      *RHO1,*RHOK,*RHOKM;


// Constants && parameters:
 double Rrho0,dsfmax;
//      Rrho0  = 2.55,     // limit for double-diff density ratio
//      dsfmax = 1.0;      // max diffusivity for salt fingering

 double 
      c0=0.0,
      c1=1.0,   
      c1p85=1.85,   
      c2=2.0,
      c3=3.0,
      c4=4.0,
      c4p6=4.6,
      c10=10.0,   
      c1000=1000.0,   
      c1p5=1.5,
      p001=0.001, 
      p015=0.015, 
      p15=0.15,
      p25=0.25,
      p5=0.5,
      p54=0.54,
      p85=0.85,
      p909=0.909,
      leps=1.0e-16;

//    these constants will be used to construct the numerator
//    factor unit change (kg/m^3 -> g/cm^3) into numerator terms

 double 
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

//    these constants will be used to construct the denomin
 double      
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

 double      
      *mwjfnums0t0,  mwjfnums0t1, *mwjfnums0t2, mwjfnums0t3, 
      *mwjfnums1t0,  mwjfnums1t1,  mwjfnums2t0,      
      *mwjfdens0t0, *mwjfdens0t1,  mwjfdens0t2, *mwjfdens0t3, mwjfdens0t4, 
       mwjfdens1t0,  mwjfdens1t1,  mwjfdens1t3,       
       mwjfdensqt0,  mwjfdensqt2;

 double *pts,*loc_pts,*tmax,*tmin,*smax,*smin,*pressz,grav;
// EOP.

 double exp_local_data[exp_data_len];
 athread_syn(ARRAY_SCOPE, 0xffff);
 if (_MYID == 0){
   get_reply = 0;
   athread_get(BCAST_MODE, exp_data, exp_local_data, exp_data_bytes, &get_reply, 0xff, 0, 0);
   while (get_reply != 1);
   sync_mem;
 }
 athread_syn(ARRAY_SCOPE, 0xffff);
 exp_data_local_ptr = exp_local_data;




#ifdef POPCHKJN
 num_calls = num_calls+1;      
#endif
    
 myid = athread_get_id(-1);
 myproc = st2_slave->core_info[0];
 threads_num = st2_slave->core_info[1];

 nx_block = st2_slave->nblock[0];
 ny_block = st2_slave->nblock[1];
 km = st2_slave->nblock[2];
 nt = st2_slave->nblock[3];
 knxt0 = st2_slave->nblock[4];
 kup0 = st2_slave->nblock[5];


 TRCR = st2_slave->local_arr[0];
 VDC = st2_slave->local_arr[1];


 pts = st2_slave->local_pts[0];
 loc_pts=(double *)ldm_malloc(8*(5*km+3));

 get_reply=0;
 athread_get(PE_MODE,pts,loc_pts,8*(5*km+3),(void*)&get_reply,0,0,0);

#ifdef POPCHKJN
 mwork11 = st2_slave->sdiag[0];  
 mwork1k = st2_slave->sdiag[1];  
 mwork1km = st2_slave->sdiag[2];  
 mwork21 = st2_slave->sdiag[3];  
 mwork2k = st2_slave->sdiag[4];  
 mwork2km = st2_slave->sdiag[5];  
#endif


 mwjfnums0t0=(double*)ldm_malloc(8*km);
 mwjfnums0t2=(double*)ldm_malloc(8*km);
 mwjfnums1t0=(double*)ldm_malloc(8*km);
 mwjfdens0t0=(double*)ldm_malloc(8*km);
 mwjfdens0t1=(double*)ldm_malloc(8*km);
 mwjfdens0t3=(double*)ldm_malloc(8*km);

 //max_ldm=get_allocatable_size();
 
 mx_usize=MAXUSZ1;
 if(km*nx_block<=mx_usize)
   {
    ni=1;
    nj=1;
    usize=2*nx_block*km;
    usize1=nx_block;
    usize2=2*nx_block*(km+2);
    usize3=nx_block*km;
   }
 else if(km<=mx_usize)  // general case!
   {
    ni=1;
    nj=nx_block;
    usize=2*km;
    usize1=1;
    usize2=2*(km+2);
    usize3=km;
   }
 else                   // extreme case
   {
    // MUST to modify! 
    // when km is larger than MAXUSZ.
   }
                


 task_all=ny_block*nj*ni;
 task_num0=task_all/threads_num;
 if(myid<task_all%threads_num)
   {
    task_num=task_num0+1;
    task_pos=task_num*myid;
   }
 else
   {
    task_num=task_num0;
    task_pos=(task_all%threads_num)+task_num0*myid;
   }


 while(get_reply!=1);
 pressz=loc_pts;
 tmax=loc_pts+km;
 tmin=loc_pts+2*km;
 smax=loc_pts+3*km;
 smin=loc_pts+4*km;
 grav=loc_pts[5*km];
 Rrho0=loc_pts[5*km+1];
 dsfmax=loc_pts[5*km+2];


#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0)
             {
             cpe_printf("grav_slave:%.18lf \n",grav);
             }
#endif



if(task_num>0)
 {

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


 for(k=0;k<km;k++) 
    {
    p=c10*pressz[k];

    mwjfnums0t0[k] = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0);
    mwjfnums0t2[k] = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2);
    mwjfnums1t0[k] = mwjfnp0s1t0 + p*mwjfnp1s1t0;

    mwjfdens0t0[k] = mwjfdp0s0t0 + p*mwjfdp1s0t0;
    mwjfdens0t1[k] = mwjfdp0s0t1 + p*p * mwjfdp3s0t1*p;
    mwjfdens0t3[k] = mwjfdp0s0t3 + p * mwjfdp2s0t3*p;
    }


 task_tile=task_num/(mx_usize/usize3);       // use usize3 as beseline
 // per tile has (task_num/task_tile) tasks with each task ~(8*usize3*4).
 if(task_tile>0)
   {
   task_size = mx_usize/usize3;
   last_size = task_num%task_size;
   }
 else
   {
   task_size = task_num;
   last_size = task_size;
   }
 
#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0)
             {
             itmp = get_allocatable_size();
             cpe_printf("max0 ldm_size: %d \n",itmp);
             }
#endif

 
 a_buf=(double*)ldm_malloc(8*usize*task_size);      //TRCR
 RRHO=(double*)ldm_malloc(8*usize3*task_size);
 ALPHADT=(double*)ldm_malloc(8*usize3*task_size);
 BETADS=(double*)ldm_malloc(8*usize3*task_size);

 TALPHA=(double*)ldm_malloc(8*usize1*2);
 SBETA=(double*)ldm_malloc(8*usize1*2);

 e_buf=(double*)ldm_malloc(8*usize2*task_size);     //VDC
 
#ifdef POPCHKJN
 swork11=(double*)ldm_malloc(8*usize3*task_size);
 swork1k=(double*)ldm_malloc(8*usize3*task_size);
 swork1km=(double*)ldm_malloc(8*usize3*task_size);
 swork21=(double*)ldm_malloc(8*usize2*task_size);
 swork2k=(double*)ldm_malloc(8*usize3*task_size);
 swork2km=(double*)ldm_malloc(8*usize3*task_size);
 for(k=0;k<(usize3*task_size);k++)
    {
    swork11[k]=0.0;
    swork1k[k]=0.0;
    swork1km[k]=0.0;
    swork2k[k]=0.0;
    swork2km[k]=0.0;
    }
 for(k=0;k<(usize2*task_size);k++)
    {
    swork21[k]=0.0;
    }
#endif



#ifdef POPCHKJN
 put_reply=7;
#else
 put_reply=1;
#endif

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0)
             {
             itmp = get_allocatable_size();
             cpe_printf("max1 ldm_size: %d \n",itmp);
             }
#endif


 for(l=0;l<task_tile;l++)
    {
    offset=task_pos*usize;
    offset2=task_pos*usize2;
    offset3=task_pos*usize3;


    get_reply=0;
    offsetg=task_pos;
    athread_get(PE_MODE,TRCR+offsetg,a_buf,8*usize*task_size,(void*)&get_reply,0,8*(ny_block*nx_block-task_size*usize1),8*task_size*usize1);
//  a_buf saved as (nt,km,task_size,usize1) with nt=2, indexed by ((j*km+k)*task_size+i)*usize1+nm,
//     while used as (task_size,usize1,km,nt), indexed by ((i*usize1+nm)*km+k)*nt+j.

    get_reply2=0;
    offsetg=task_pos;
//  e_buf saved as (nt,km,task_size,usize1) with nt=2, indexed by ((j*km+k)*task_size+i)*usize1+nm,
//     while used as (task_size,usize1,km,nt), indexed by ((i*usize1+nm)*km+k)*nt+j.
    athread_get(PE_MODE,VDC+offsetg,e_buf,8*usize2*task_size,(void*)&get_reply2,0,8*(ny_block*nx_block-task_size*usize1),8*task_size*usize1);

    while(get_reply!=1);

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0 && l==1)
             cpe_printf("tile1_get_TRCR: ok! \n");
#endif

    for(i=0;i<task_size;i++)
       {

       for(nm=0;nm<usize1;nm++)     
        {
        knxt  = knxt0;
        kup = kup0;


/*
!  compute alpha*DT and beta*DS at interfaces.  use RRHO and
!  PRANDTL for temporary storage for call to state
*/

        nk=nm*km;
  
        //  PRANDTL = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)
/*
        if(a_buf[usize*i+2*nm*km]>-c2)
         PRANDTL[i*usize1+nm]=a_buf[usize*i+2*nm*km];
*/
        // transpose to shape (km,task_size,usize1)
        jtmp=i*usize1+nm;              // t with j=0, k=0
        if(a_buf[jtmp]>-c2)
         PRANDTL=a_buf[jtmp];
        else
         PRANDTL=-c2;

         // call state(1, 1, PRANDTL,  TRCR(:,:,1  ,2),  this_block, RHOFULL=RRHO,DRHODT=TALPHA(:,:,kup), DRHODS=SBETA(:,:,kup))
	 TEMPK = PRANDTL;    
         if(TEMPK>tmax[0]) TQ=tmax[0];
         else if(TEMPK<tmin[0]) TQ=tmin[0];
         else TQ=TEMPK;      

//	 SALTK = a_buf[1+usize*i+2*nm*km];
         SALTK = a_buf[(km*task_size+i)*usize1+nm];   // s with j=1, k=0
	 if(SALTK>smax[0]) SQ=smax[0];
         else if(SALTK<smin[0]) SQ=smin[0];
         else SQ=SALTK;      

         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[0] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[0] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[0] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[0] + TQ * (mwjfdens0t1[0] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[0] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));


         DENOMK = c1/WORK2;
//         RHO1[nk+usize3*i] = WORK1*DENOMK;
         RRHO[jtmp] = WORK1*DENOMK;

#ifdef POPCHKJN
//         swork11[nk+usize3*i]=RRHO[jtmp];
#endif


         WORK3 =           // dP_1/dT
                 mwjfnums0t1 + TQ * (c2*mwjfnums0t2[0] +  
                 c3*mwjfnums0t3 * TQ) + mwjfnums1t1 * SQ;

         WORK4 =           // dP_2/dT
                 mwjfdens0t1[0] + SQ * mwjfdens1t1 +   
                 TQ * (c2*(mwjfdens0t2 + SQ*SQR*mwjfdensqt2) + 
                 TQ * (c3*(mwjfdens0t3[0] + SQ * mwjfdens1t3) +
                 TQ *  c4*mwjfdens0t4));

         TALPHA[2*nm+kup] = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK;
#ifdef POPCHKJN
         swork1k[nk+usize3*i]=TALPHA[2*nm+kup];
#endif


         WORK3 =           // dP_1/dS
                 mwjfnums1t0[0] + mwjfnums1t1 * TQ + c2*mwjfnums2t0 * SQ;

         WORK4 = mwjfdens1t0 +   // dP_2/dS
                 TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3) +  
                 c1p5*SQR*(mwjfdensqt0 + TQ*TQ*mwjfdensqt2);

         SBETA[2*nm+kup] = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK * c1000;
#ifdef POPCHKJN
         swork1km[nk+usize3*i]=SBETA[2*nm+kup];
#endif

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0)
             cpe_printf("computing state(1,1,P,T) @ tile %d ok! \n",l);
#endif



        for(k=0;k<km-1;k++)       //do k=1,km in fortran
         {
         nk=(k+1)+nm*km;    //!!! 2018-10-26
         
      
         // PRANDTL = merge(-c2,TRCR(:,:,k+1,1),TRCR(:,:,k+1,1) < -c2)
/*
         if(a_buf[usize*i+2*nk]>-c2) 
            TEMPKLP[klvl+2*(i*usize1+nm)]=a_buf[usize*i+2*nk];
*/
         // transpose to shape (km,task_size,usize1)
         jtmp=(k*task_size+i)*usize1+nm;              // t with j=0
         j1tmp=((k+1)*task_size+i)*usize1+nm;              // t with j=0
         if(a_buf[j1tmp]>-c2) 
            PRANDTL=a_buf[j1tmp];
         else
            PRANDTL=-c2; 


         // call state(k+1, k+1, PRANDTL,  TRCR(:,:,k+1,2),  this_block, RHOFULL=RRHO,DRHODT=TALPHA(:,:,knxt), DRHODS=SBETA(:,:,knxt)))
	 TEMPK = PRANDTL;       
         if(TEMPK>tmax[k+1]) TQ=tmax[k+1];
         else if(TEMPK<tmin[k+1]) TQ=tmin[k+1];
         else TQ=TEMPK;      

//	 SALTK = a_buf[1+usize*i+2*(nm*km+k)];
	 ktmp=((km+k)*task_size+i)*usize1+nm;
	 k1tmp=((km+k+1)*task_size+i)*usize1+nm;
         SALTK = a_buf[k1tmp];         // s with j=1
	 if(SALTK>smax[k+1]) SQ=smax[k+1];
         else if(SALTK<smin[k+1]) SQ=smin[k+1];
         else SQ=SALTK;      

         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[k+1] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[k+1] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k+1] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[k+1] + TQ * (mwjfdens0t1[k+1] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[k+1] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));


         DENOMK = c1/WORK2;
//         RHO1[nk+usize3*i] = WORK1*DENOMK;
         RRHO[j1tmp] = WORK1*DENOMK;

#ifdef POPCHKJN
//         swork11[nk+usize3*i]=RRHO[j1tmp];
#endif


         WORK3 =           // dP_1/dT
                 mwjfnums0t1 + TQ * (c2*mwjfnums0t2[k+1] + 
                 c3*mwjfnums0t3 * TQ) + mwjfnums1t1 * SQ;

         WORK4 =           // dP_2/dT
                 mwjfdens0t1[k+1] + SQ * mwjfdens1t1 +   
                 TQ * (c2*(mwjfdens0t2 + SQ*SQR*mwjfdensqt2) + 
                 TQ * (c3*(mwjfdens0t3[k+1] + SQ * mwjfdens1t3) +  
                 TQ *  c4*mwjfdens0t4));

         TALPHA[2*nm+knxt] = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK;
#ifdef POPCHKJN
         swork1k[nk+usize3*i]=TALPHA[2*nm+knxt];
#endif


         WORK3 =           // dP_1/dS
                 mwjfnums1t0[k+1] + mwjfnums1t1 * TQ + c2*mwjfnums2t0 * SQ;

         WORK4 = mwjfdens1t0 +   // dP_2/dS
                 TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3) +  
                 c1p5*SQR*(mwjfdensqt0 + TQ*TQ*mwjfdensqt2);

         SBETA[2*nm+knxt] = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK * c1000;
#ifdef POPCHKJN
         swork1km[nk+usize3*i]=SBETA[2*nm+knxt];
#endif


/*
         ALPHADT[nk+usize3*i]=-p5*(TALPHA[kup]+TALPHA[knxt])
                             *(a_buf[jtmp]-a_buf[j1tmp]);
         BETADS [nk+usize3*i]=-p5*(SBETA[kup]+SBETA[knxt])
                             *(a_buf[ktmp]-a_buf[k1tmp]);
*/
//       transpose to shape (km,task_size,usize1)
//       itmp=(k*task_size+i)*usize1+nm;
         ALPHADT[jtmp]=-p5*(TALPHA[2*nm+kup]+TALPHA[2*nm+knxt])
                             *(a_buf[jtmp]-a_buf[j1tmp]);
         BETADS [jtmp]=p5*(SBETA[2*nm+kup]+SBETA[2*nm+knxt])
                             *(a_buf[ktmp]-a_buf[k1tmp]);
 

         kup = knxt;
         knxt  = 1-kup;

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0)
             cpe_printf("computing state @ tile %d ok! \n",l);
#endif

        }        // end of for(k)

/*
        ALPHADT[nm*km+km-1+usize3*i]=c0;
        BETADS [nm*km+km-1+usize3*i]=c0;  // with k=0,..,km-1
*/
        ALPHADT[((km-1)*task_size+i)*usize1+nm]=c0;
        BETADS [((km-1)*task_size+i)*usize1+nm]=c0;

       }       // end of for(nm)


#if defined(POPTESTJN)
         if(num_calls%8==1 && j==0 && i==0 && k==1 && myproc==0 && myid==0)
             {
             cpe_printf("9.99843699e-1: %.16e \n",mwjfnp0s0t0);
             cpe_printf("7.35212840e-3: %.16e \n",mwjfnp0s0t1);
             cpe_printf("-5.45928211e-5: %.16e \n",mwjfnp0s0t2);
             cpe_printf("3.98476704e-7: %.16e \n",mwjfnp0s0t3);
             cpe_printf("2.96938239e-3: %.16e \n",mwjfnp0s1t0);
             cpe_printf("-7.23268813e-6: %.16e \n",mwjfnp0s1t1);
             cpe_printf("2.12382341e-6: %.16e \n",mwjfnp0s2t0);
             cpe_printf("1.04004591e-5: %.16e \n",mwjfnp1s0t0);
             cpe_printf("1.03970529e-10: %.16e \n",mwjfnp1s0t2);
             cpe_printf("5.18761880e-9: %.16e \n",mwjfnp1s1t0);
             cpe_printf("-3.24041825e-11: %.16e \n",mwjfnp2s0t0);
             cpe_printf("-1.23869360e-14: %.16e \n",mwjfnp2s0t2);

             cpe_printf("1.0e+0: %.16e \n",mwjfdp0s0t0);
             cpe_printf("7.28606739e-3: %.16e \n",mwjfdp0s0t1);
             cpe_printf("-4.60835542e-5: %.16e \n",mwjfdp0s0t2);
             cpe_printf("3.68390573e-7: %.16e \n",mwjfdp0s0t3);
             cpe_printf("1.80809186e-10: %.16e \n",mwjfdp0s0t4);
             cpe_printf("2.14691708e-3: %.16e \n",mwjfdp0s1t0);
             cpe_printf("-9.27062484e-6: %.16e \n",mwjfdp0s1t1);
             cpe_printf("-1.78343643e-10: %.16e \n",mwjfdp0s1t3);
             cpe_printf("4.76534122e-6: %.16e \n",mwjfdp0sqt0);
             cpe_printf("1.63410736e-9: %.16e \n",mwjfdp0sqt2);
             cpe_printf("5.30848875e-6: %.16e \n",mwjfdp1s0t0);
             cpe_printf("-3.03175128e-16: %.16e \n",mwjfdp2s0t3);
             cpe_printf("-1.27934137e-17: %.16e \n",mwjfdp3s0t1);
             }
#endif

    }    // end of for(i)


    while(get_reply2!=1);

#if defined(POPTESTJN)
    if(num_calls<=2 && myproc==0 && myid==0)
             cpe_printf("tile_get_VDC %d:  ok! \n",l);
#endif



    // put with stride at slave: 
#ifdef POPCHKJN
    while(put_reply!=7);
#else
    while(put_reply!=1);
#endif


    for(k=0;k<km;k++)
       for(i=0;i<task_size;i++)
          for(nm=0;nm<usize1;nm++)
             {
             itmp =(  k      *task_size+i)*usize1+nm; 
             ktmp =((km+2+k)   *task_size+i)*usize1+nm;
             i1tmp =((k+1)   *task_size+i)*usize1+nm; 
             k1tmp =((km+2+k+1)*task_size+i)*usize1+nm;
#ifdef POPCHKJN
             DIFFDD=c0;
             tmpw1=c0;
             tmpw2=c0;
#endif

//!     salt fingering case

             if(ALPHADT[itmp]>BETADS[itmp] && BETADS[itmp]>c0)
               {
               tmprho=ALPHADT[itmp]/BETADS[itmp];
               if(tmprho>Rrho0) tmprho=Rrho0;
               etmp=c1-(tmprho-c1)/(Rrho0-c1);
               DIFFDD=dsfmax*etmp*etmp*etmp;
               e_buf[i1tmp]=e_buf[i1tmp]+0.7*DIFFDD;
               e_buf[k1tmp]=e_buf[k1tmp]+DIFFDD;
               }
#ifdef POPCHKJN
             swork11[(i*usize1+nm)*km+k]=DIFFDD;
             swork2km[(i*usize1+nm)*km+k]=c0;
#endif
 
//!     diffusive convection

//             if(ALPHADT[itmp]>BETADS[itmp] && BETADS[itmp]<c0 && ALPHADT[itmp]<c0)
             if(ALPHADT[itmp]>BETADS[itmp] && ALPHADT[itmp]<c0)
               {
               tmprho=ALPHADT[itmp]/BETADS[itmp];
               DIFFDD=p015*p909*exp(c4p6*exp(-p54*(c1/tmprho-c1)));

#ifdef POPCHKJN
        swork2km[(i*usize1+nm)*km+k]=-p54*(c1/tmprho-c1);
               tmpw1=exp(-p54*(c1/tmprho-c1));
               tmpw2=exp(c4p6*tmpw1);
               DIFFDD=p015*p909*tmpw2;
#endif
               PRANDTL=0.15*tmprho;            
               }
             else
               {
               tmprho=c0;
               DIFFDD=c0;
               PRANDTL=c0;
               }

//             if(tmprho>p5) PRANDTL=(1.85-0.85/tmprho)*tmprho;
             if(tmprho>p5) PRANDTL=1.85*tmprho-0.85;

#ifdef POPCHKJN
        swork1k[(i*usize1+nm)*km+k]=tmpw1;
        swork1km[(i*usize1+nm)*km+k]=tmpw2;
        swork2k[(i*usize1+nm)*km+k]=DIFFDD;
//        swork2km[(i*usize1+nm)*km+k]=PRANDTL;
#endif
             e_buf[i1tmp]=e_buf[i1tmp]+DIFFDD;
             e_buf[k1tmp]=e_buf[k1tmp]+PRANDTL*DIFFDD;
             }

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0)
             cpe_printf("tile_computing_VDC @ tile %d ok! \n",l);
#endif


#ifdef POPCHKJN
    for(i=0;i<task_size;i++)
    for(nm=0;nm<usize1;nm++)
      {
      for(k=0;k<=km+1;k++)
        {
        swork21[(i*usize1+nm)*2*(km+2)+2*k]=e_buf[(k*task_size+i)*usize1+nm];
        swork21[(i*usize1+nm)*2*(km+2)+2*k+1]=e_buf[((k+km+2)*task_size+i)*usize1+nm];
        }
/*
      for(k=0;k<km;k++)
        {
        swork2k[(i*usize1+nm)*km+k]=ALPHADT[(k*task_size+i)*usize1+nm];
        swork2km[(i*usize1+nm)*km+k]=BETADS[(k*task_size+i)*usize1+nm];
        }
*/
      }
#endif


#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0)
             cpe_printf("tile_computing_VDC&&swork @ tile %d ok! \n",l);
#endif


/*     
     put_reply=0;
     athread_put(PE_MODE,e_buf,VDC+offset3,8*2*usize2*task_size,(void*)&put_reply,0,0);
     while(put_reply!=1);
*/

     // put with stride at slave: 
     put_reply=0;
     offsetp3=task_pos;
     athread_put(PE_MODE,e_buf,VDC+offsetp3,8*usize2*task_size,(void*)&put_reply,8*(nx_block*ny_block-task_size*usize1),8*task_size*usize1);

#if defined(POPTESTJN)
     while(put_reply!=1);
         if(num_calls<=2 && myproc==0 && myid==0)
             cpe_printf("tile_put_VDC ok! \n");
#endif



#ifdef POPCHKJN
     athread_put(PE_MODE,swork11, mwork11+offset3, 8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork1k, mwork1k+offset3, 8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork1km,mwork1km+offset3,8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork21, mwork21+offset2, 8*usize2*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork2k, mwork2k+offset3, 8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork2km, mwork2km+offset3, 8*usize3*task_size,(void*)&put_reply,0,0);
#endif


#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0 && l==1)
             cpe_printf("tile1_put_VDC&&swork: ok! \n");
#endif

     
     task_pos = task_pos + task_size;
   }
   // end of for(l)



#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==0)
             cpe_printf("step_into_last_part: \n");
#endif

 // last_size
 if(last_size>0)
   {
   offset=task_pos*usize;
   offset2=task_pos*usize2;
   offset3=task_pos*usize3;

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("step_into_last_part! \n");
#endif


    get_reply=0;
    offsetg=task_pos;
    athread_get(PE_MODE,TRCR+offsetg,a_buf,8*usize*last_size,(void*)&get_reply,0,8*(ny_block*nx_block-last_size*usize1),8*last_size*usize1);
//  a_buf saved as (nt,km,task_size,usize1) with nt=2, indexed by ((j*km+k)*task_size+i)*usize1+nm,
//  while used as (task_size,usize1,km,nt), indexed by ((i*usize1+nm)*km+k)*nt+j.

    get_reply2=0;
    offsetg=task_pos;
    athread_get(PE_MODE,VDC+offsetg,e_buf,8*usize2*last_size,(void*)&get_reply2,0,8*(ny_block*nx_block-last_size*usize1),8*last_size*usize1);

    while(get_reply!=1);
#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("last_part: get_TRCR&&VDC ok! \n");
#endif



    for(i=0;i<last_size;i++)
       {

       for(nm=0;nm<usize1;nm++)     
        {
        knxt  = knxt0;
        kup = kup0;
        

/*
!  compute alpha*DT and beta*DS at interfaces.  use RRHO and
!  PRANDTL for temporary storage for call to state
*/

        nk=nm*km;
        
        //  PRANDTL = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)
/*
        if(a_buf[usize*i+2*nm*km]>-c2)
         PRANDTL[i*usize1+nm]=a_buf[usize*i+2*nm*km];
*/
        jtmp=i*usize1+nm;              // t with j=0, k=0
        if(a_buf[jtmp]>-c2)
         PRANDTL=a_buf[jtmp];
        else
         PRANDTL=-c2;

         // call state(1, 1, PRANDTL,  TRCR(:,:,1  ,2),  this_block, RHOFULL=RRHO,DRHODT=TALPHA(:,:,kup), DRHODS=SBETA(:,:,kup))
	 TEMPK = PRANDTL;    
         if(TEMPK>tmax[0]) TQ=tmax[0];
         else if(TEMPK<tmin[0]) TQ=tmin[0];
         else TQ=TEMPK;      

//	 SALTK = a_buf[1+usize*i+2*nm*km];
         SALTK = a_buf[(km*last_size+i)*usize1+nm];   // s with j=1, k=0
	 if(SALTK>smax[0]) SQ=smax[0];
         else if(SALTK<smin[0]) SQ=smin[0];
         else SQ=SALTK;      

         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[0] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[0] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[0] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[0] + TQ * (mwjfdens0t1[0] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[0] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));


         DENOMK = c1/WORK2;
//         RHO1[nk+usize3*i] = WORK1*DENOMK;
         RRHO[jtmp] = WORK1*DENOMK;

#ifdef POPCHKJN
//         swork11[nk+usize3*i]=RRHO[jtmp];
#endif


         WORK3 =           // dP_1/dT
                 mwjfnums0t1 + TQ * (c2*mwjfnums0t2[0] +  
                 c3*mwjfnums0t3 * TQ) + mwjfnums1t1 * SQ;

         WORK4 =           // dP_2/dT
                 mwjfdens0t1[0] + SQ * mwjfdens1t1 +   
                 TQ * (c2*(mwjfdens0t2 + SQ*SQR*mwjfdensqt2) + 
                 TQ * (c3*(mwjfdens0t3[0] + SQ * mwjfdens1t3) +
                 TQ *  c4*mwjfdens0t4));

         TALPHA[2*nm+kup] = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK;
#ifdef POPCHKJN
         swork1k[nk+usize3*i]=TALPHA[2*nm+kup];
#endif


         WORK3 =           // dP_1/dS
                 mwjfnums1t0[0] + mwjfnums1t1 * TQ + c2*mwjfnums2t0 * SQ;

         WORK4 = mwjfdens1t0 +   // dP_2/dS
                 TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3) +  
                 c1p5*SQR*(mwjfdensqt0 + TQ*TQ*mwjfdensqt2);

         SBETA[2*nm+kup] = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK * c1000;
#ifdef POPCHKJN
         swork1km[nk+usize3*i]=SBETA[2*nm+kup];
#endif

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("computing state(1,1,P,T) @ last_part task ok! \n");
#endif



        for(k=0;k<km-1;k++)       //do k=1,km in fortran
         {
         nk=(k+1)+nm*km;    //!!! 2018-10-26
      
         // PRANDTL = merge(-c2,TRCR(:,:,k+1,1),TRCR(:,:,k+1,1) < -c2)
/*
         if(a_buf[usize*i+2*nk]>-c2) 
            TEMPKLP[klvl+2*(i*usize1+nm)]=a_buf[usize*i+2*nk];
*/
         jtmp=(k*last_size+i)*usize1+nm;              // t with j=0
         j1tmp=((k+1)*last_size+i)*usize1+nm;              // t with j=0
         if(a_buf[j1tmp]>-c2) 
            PRANDTL=a_buf[j1tmp];
         else
            PRANDTL=-c2; 


         // call state(k+1, k+1, PRANDTL,  TRCR(:,:,k+1,2),  this_block, RHOFULL=RRHO,DRHODT=TALPHA(:,:,knxt), DRHODS=SBETA(:,:,knxt)))
	 TEMPK = PRANDTL;       
         if(TEMPK>tmax[k+1]) TQ=tmax[k+1];
         else if(TEMPK<tmin[k+1]) TQ=tmin[k+1];
         else TQ=TEMPK;      

//	 SALTK = a_buf[1+usize*i+2*(nm*km+k)];
	 ktmp=((km+k)*last_size+i)*usize1+nm;
	 k1tmp=((km+k+1)*last_size+i)*usize1+nm;
         SALTK = a_buf[k1tmp];         // s with j=1
	 if(SALTK>smax[k+1]) SQ=smax[k+1];
         else if(SALTK<smin[k+1]) SQ=smin[k+1];
         else SQ=SALTK;      

         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[k+1] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[k+1] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k+1] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[k+1] + TQ * (mwjfdens0t1[k+1] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[k+1] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));


         DENOMK = c1/WORK2;
//         RHO1[nk+usize3*i] = WORK1*DENOMK;
         RRHO[j1tmp] = WORK1*DENOMK;

#ifdef POPCHKJN
//         swork11[nk+usize3*i]=RRHO[j1tmp];
#endif


         WORK3 =           // dP_1/dT
                 mwjfnums0t1 + TQ * (c2*mwjfnums0t2[k+1] + 
                 c3*mwjfnums0t3 * TQ) + mwjfnums1t1 * SQ;

         WORK4 =           // dP_2/dT
                 mwjfdens0t1[k+1] + SQ * mwjfdens1t1 +   
                 TQ * (c2*(mwjfdens0t2 + SQ*SQR*mwjfdensqt2) + 
                 TQ * (c3*(mwjfdens0t3[k+1] + SQ * mwjfdens1t3) +  
                 TQ *  c4*mwjfdens0t4));

         TALPHA[2*nm+knxt] = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK;
#ifdef POPCHKJN
         swork1k[nk+usize3*i]=TALPHA[2*nm+knxt];
#endif


         WORK3 =           // dP_1/dS
                 mwjfnums1t0[k+1] + mwjfnums1t1 * TQ + c2*mwjfnums2t0 * SQ;

         WORK4 = mwjfdens1t0 +   // dP_2/dS
                 TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3) +  
                 c1p5*SQR*(mwjfdensqt0 + TQ*TQ*mwjfdensqt2);

         SBETA[2*nm+knxt] = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK * c1000;
#ifdef POPCHKJN
         swork1km[nk+usize3*i]=SBETA[2*nm+knxt];
#endif


/*
         ALPHADT[nk+usize3*i]=-p5*(TALPHA[kup]+TALPHA[knxt])
                             *(a_buf[jtmp]-a_buf[j1tmp]);
         BETADS [nk+usize3*i]=p5*(SBETA[kup]+SBETA[knxt])
                             *(a_buf[ktmp]-a_buf[k1tmp]);
*/
//       transpose to shape (km,task_size,usize1)
//       itmp=(k*task_size+i)*usize1+nm;
         ALPHADT[jtmp]=-p5*(TALPHA[2*nm+kup]+TALPHA[2*nm+knxt])
                             *(a_buf[jtmp]-a_buf[j1tmp]);
         BETADS [jtmp]=p5*(SBETA[2*nm+kup]+SBETA[2*nm+knxt])
                             *(a_buf[ktmp]-a_buf[k1tmp]);
 

         kup = knxt;
         knxt  = 1-kup;

        }        // end of for(k)

/*
        ALPHADT[nm*km+km-1+usize3*i]=c0;
        BETADS [nm*km+km-1+usize3*i]=c0;
*/
        ALPHADT[((km-1)*last_size+i)*usize1+nm]=c0;
        BETADS [((km-1)*last_size+i)*usize1+nm]=c0;

       }       // end of for(nm)

    }    // end of for(i)

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("computing_state@last_part ok! \n");
#endif



    while(get_reply2!=1);


    // put with stride at slave: 
#ifdef POPCHKJN
    while(put_reply!=7);
#else
    while(put_reply!=1);
#endif



    for(k=0;k<km;k++)
       for(i=0;i<last_size;i++)
          for(nm=0;nm<usize1;nm++)
             {
             itmp =(  k      *last_size+i)*usize1+nm;
             i1tmp =((k+1)   *last_size+i)*usize1+nm;
             ktmp =((km+2+k)   *last_size+i)*usize1+nm;
             k1tmp =((km+2+k+1)*last_size+i)*usize1+nm;
#ifdef POPCHKJN
             DIFFDD=c0;
             tmpw1=c0;
             tmpw2=c0;
#endif


//!     salt fingering case
             if(ALPHADT[itmp]>BETADS[itmp] && BETADS[itmp]>c0)
               {
               tmprho=ALPHADT[itmp]/BETADS[itmp];
               if(tmprho>Rrho0) tmprho=Rrho0;
               etmp=c1-(tmprho-c1)/(Rrho0-c1);
               DIFFDD=dsfmax*etmp*etmp*etmp;
               e_buf[i1tmp]=e_buf[i1tmp]+0.7*DIFFDD;
               e_buf[k1tmp]=e_buf[k1tmp]+DIFFDD;
               }

#ifdef POPCHKJN
        swork11[(i*usize1+nm)*km+k]=DIFFDD;
        swork2km[(i*usize1+nm)*km+k]=c0;
#endif

 //!     diffusive convection
             if(ALPHADT[itmp]>BETADS[itmp] && BETADS[itmp]<c0 && ALPHADT[itmp]<c0)
               {
               tmprho=ALPHADT[itmp]/BETADS[itmp];
//               DIFFDD=1.5e-2*0.909*exp(4.6*exp(-0.54*(c1/tmprho-c1)));
               DIFFDD=p015*p909*exp(c4p6*exp(-p54*(c1/tmprho-c1)));
#ifdef POPCHKJN
               swork2km[(i*usize1+nm)*km+k]=-p54*(c1/tmprho-c1);
               tmpw1=exp(-p54*(c1/tmprho-c1));
               tmpw2=exp(c4p6*tmpw1);
               DIFFDD=p015*p909*tmpw2;
#endif
               PRANDTL=0.15*tmprho;
               }
             else
               {
               tmprho=c0;
               DIFFDD=c0;
               PRANDTL=c0;
               }

             if(tmprho>p5) PRANDTL=(1.85-0.85/tmprho)*tmprho;

#ifdef POPCHKJN
        swork1k[(i*usize1+nm)*km+k]=tmpw1;
        swork1km[(i*usize1+nm)*km+k]=tmpw2;
        swork2k[(i*usize1+nm)*km+k]=DIFFDD;
//        swork2km[(i*usize1+nm)*km+k]=PRANDTL;
#endif
             e_buf[i1tmp]=e_buf[i1tmp]+DIFFDD;
             e_buf[k1tmp]=e_buf[k1tmp]+PRANDTL*DIFFDD;

             }

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("computing_VDC@last_part task ok! \n");
#endif


#ifdef POPCHKJN
    for(i=0;i<last_size;i++)
    for(nm=0;nm<usize1;nm++)
      {
      for(k=0;k<=km+1;k++)
        {
        swork21[(i*usize1+nm)*2*(km+2)+2*k]=e_buf[(k*last_size+i)*usize1+nm];
        swork21[(i*usize1+nm)*2*(km+2)+2*k+1]=e_buf[((k+km+2)*last_size+i)*usize1+nm];
        }
/*
      for(k=0;k<km;k++)
        {
        swork2k[(i*usize1+nm)*km+k]=ALPHADT[(k*last_size+i)*usize1+nm];
        swork2km[(i*usize1+nm)*km+k]=BETADS[(k*last_size+i)*usize1+nm];
        }
*/
      }
#endif

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("computing_VDC&swork @last_part task ok! \n");
#endif


     // put with stride at slave: 
     put_reply=0;
     offsetp3=task_pos;
     athread_put(PE_MODE,e_buf,VDC+offsetp3,8*usize2*last_size,(void*)&put_reply,8*(nx_block*ny_block-last_size*usize1),8*last_size*usize1);

#if defined(POPTESTJN)
     while(put_reply!=1);
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("put_VDC@last_part ok! \n");
#endif



#ifdef POPCHKJN
     athread_put(PE_MODE,swork11, mwork11+offset3, 8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork1k, mwork1k+offset3, 8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork1km,mwork1km+offset3,8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork21, mwork21+offset2, 8*usize2*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork2k, mwork2k+offset3, 8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork2km, mwork2km+offset3, 8*usize3*last_size,(void*)&put_reply,0,0);
#endif

     
#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("put_VDC&swork @last_part ok! \n");
#endif


   }
   // end of if(last_size>0)
    

 ldm_free(a_buf,8*usize*task_size);

 ldm_free(TALPHA,8*usize1*2);
 ldm_free(SBETA,8*usize1*2);

 ldm_free(RRHO,8*usize3*task_size);
 ldm_free(ALPHADT,8*usize3*task_size);
 ldm_free(BETADS,8*usize3*task_size);


 // put with stride at slave: 
#ifdef POPCHKJN
 while(put_reply!=7);
#else
 while(put_reply!=1);
#endif
 ldm_free(e_buf,8*usize2*task_size);


#ifdef POPCHKJN
 ldm_free(swork11,8*usize3*task_size);
 ldm_free(swork1k,8*usize3*task_size);
 ldm_free(swork1km,8*usize3*task_size);
 ldm_free(swork21,8*usize2*task_size);
 ldm_free(swork2k,8*usize3*task_size);
 ldm_free(swork2km,8*usize3*task_size);
#endif

}
// end of if(task_num>0)


 ldm_free(loc_pts,8*(5*km+3));

 ldm_free(mwjfnums0t0,8*km);
 ldm_free(mwjfnums0t2,8*km);
 ldm_free(mwjfnums1t0,8*km);
 ldm_free(mwjfdens0t0,8*km);
 ldm_free(mwjfdens0t1,8*km);
 ldm_free(mwjfdens0t3,8*km);

#if defined(POPTESTJN)
         if(num_calls<=2 && myproc==0 && myid==50)
             cpe_printf("ddmix_slave finished! \n");
#endif
  exp_data_local_ptr = exp_data;

}
//------------------------------------------------------------






//------------------------------------------------------------------
void slave_kpp_state_loop2(struct param_kpp_state_s2 *st2_slave)
{
	
 volatile int get_reply, put_reply;
#ifdef POPCHKJN
 volatile long num_calls=0;
#endif

 int myproc,myid,threads_num;
 int task_num0,task_num,task_pos,task_size,task_tile,task_all,last_size;
 int nx_block,ny_block,km,nt,klvl,kprev,klvl0,kprev0,lx_block;
 int *KMT;

 double *TRCR,*TEMPSFC,*TEMPKLP,*DBLOC,*DBSFC;

 int offset,offset1,offset2,offset3;
 long offsetg,offsetg1,offsetg2,offsetg3;
 long offsetp,offsetp1,offsetp2,offsetp3;
 int max_ldm,mx_usize,usize,usize1,usize2,usize3;
 int i,j,k,n,m,l,ni,nj,nk,nm,nitr,ktmp;

 double *a_buf,*b_buf,*c_buf,*e_buf,*f_buf,*g_buf;
 double *a0_buf,*b0_buf,*c0_buf,*e0_buf,*f0_buf,*g0_buf;

 int *d_buf,*d0_buf;

 int itmp,i1tmp,jtmp,j1tmp;
 double etmp;


 double 
      TQ, SQ,             // adjusted T,S
      TEMPK, SALTK, SQR,     // temporary variables
      p,                  // temporary pressure scalar
      DENOMK,WORK1,WORK2,
#ifdef POPCHKJN
      *mwork11,*mwork1k,*mwork1km,
      *mwork21,*mwork2k,
      *swork11,*swork1k,*swork1km,
      *swork21,*swork2k,
#endif
      tmpw1,tmpw2,
      *RHO1,*RHOK,*RHOKM;

 double 
      c0=0.0,
      c1=1.0,   
      c2=2.0,
      nc1=-1.0,
      nc2=-2.0,
      c10=10.0,   
      c1000=1000.0,   
      p001=0.001, 
      p25=0.25,
      p5=0.5,
      leps=1.0e-16;

//    these constants will be used to construct the numerator
//    factor unit change (kg/m^3 -> g/cm^3) into numerator terms

 double 
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

//    these constants will be used to construct the denomin
 double      
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

 double      
      *mwjfnums0t0,  mwjfnums0t1, *mwjfnums0t2, mwjfnums0t3, 
      *mwjfnums1t0,  mwjfnums1t1,  mwjfnums2t0,      
      *mwjfdens0t0, *mwjfdens0t1,  mwjfdens0t2, *mwjfdens0t3, mwjfdens0t4, 
       mwjfdens1t0,  mwjfdens1t1,  mwjfdens1t3,       
       mwjfdensqt0,  mwjfdensqt2;

 double *pts,*loc_pts,*tmax,*tmin,*smax,*smin,*pressz,grav;
	
	
#ifdef POPCHKJN
 num_calls = num_calls+1;      
#endif
    
 myid = athread_get_id(-1);
 myproc = st2_slave->core_info[0];
 threads_num = st2_slave->core_info[1];

 nx_block = st2_slave->nblock[0];
 ny_block = st2_slave->nblock[1];
 km = st2_slave->nblock[2];
 nt = st2_slave->nblock[3];
 klvl0 = st2_slave->nblock[4];
 kprev0 = st2_slave->nblock[5];

 KMT = st2_slave->local_iarr[0];

 TRCR = st2_slave->local_arr[0];
 DBLOC = st2_slave->local_arr[1];
 DBSFC = st2_slave->local_arr[2];


 pts = st2_slave->local_pts[0];
 loc_pts=(double *)ldm_malloc(8*(5*km+1));

 get_reply=0;
 athread_get(PE_MODE,pts,loc_pts,8*(5*km+1),(void*)&get_reply,0,0,0);

#ifdef POPCHKJN
 mwork11 = st2_slave->sdiag[0];  
 mwork1k = st2_slave->sdiag[1];  
 mwork1km = st2_slave->sdiag[2];  
 mwork21 = st2_slave->sdiag[3];  
 mwork2k = st2_slave->sdiag[4];  
#endif


 mwjfnums0t0=(double*)ldm_malloc(8*km);
 mwjfnums0t2=(double*)ldm_malloc(8*km);
 mwjfnums1t0=(double*)ldm_malloc(8*km);
 mwjfdens0t0=(double*)ldm_malloc(8*km);
 mwjfdens0t1=(double*)ldm_malloc(8*km);
 mwjfdens0t3=(double*)ldm_malloc(8*km);

 //max_ldm=get_allocatable_size();
 
 mx_usize=MAXUSZ;
 if(km*nx_block<=mx_usize)
   {
    ni=1;
    nj=1;
    usize=2*nx_block*km;
    usize1=nx_block;
    usize3=nx_block*km;
   }
 else if(km<=mx_usize)  // general case!
   {
    ni=1;
    nj=nx_block;
    usize=2*km;
    usize1=1;
    usize3=km;
   }
 else                   // extreme case
   {
    // MUST to modify! 
    // when km is larger than MAXUSZ.
   }
                


 task_all=ny_block*nj*ni;
 task_num0=task_all/threads_num;
 if(myid<task_all%threads_num)
   {
    task_num=task_num0+1;
    task_pos=task_num*myid;
   }
 else
   {
    task_num=task_num0;
    task_pos=(task_all%threads_num)+task_num0*myid;
   }


 while(get_reply!=1);
 pressz=loc_pts;
 tmax=loc_pts+km;
 tmin=loc_pts+2*km;
 smax=loc_pts+3*km;
 smin=loc_pts+4*km;
 grav=loc_pts[5*km];




if(task_num>0)
 {

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


 for(k=0;k<km;k++) 
    {
    p=c10*pressz[k];

    mwjfnums0t0[k] = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0);
    mwjfnums0t2[k] = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2);
    mwjfnums1t0[k] = mwjfnp0s1t0 + p*mwjfnp1s1t0;

    mwjfdens0t0[k] = mwjfdp0s0t0 + p*mwjfdp1s0t0;
    mwjfdens0t1[k] = mwjfdp0s0t1 + p*p * mwjfdp3s0t1*p;
    mwjfdens0t3[k] = mwjfdp0s0t3 + p * mwjfdp2s0t3*p;
    }


 task_tile=task_num/(mx_usize/usize3);       // use usize3 as beseline
 // per tile has (task_num/task_tile) tasks with each task ~(8*usize3*4).
 if(task_tile>0)
   {
   task_size = mx_usize/usize3;
   last_size = task_num%task_size;
   }
 else
   {
   task_size = task_num;
   last_size = task_size;
   }
 
 
 
 a_buf=(double*)ldm_malloc(8*usize*task_size);      //TRCR
 d_buf=(int*)ldm_malloc(SIZEINT*usize1*task_size);     //KMT

 TEMPSFC=(double*)ldm_malloc(8*usize1*task_size);
 TEMPKLP=(double*)ldm_malloc(8*usize1*2*task_size);

 RHO1=(double*)ldm_malloc(8*usize3*task_size);
 RHOK=(double*)ldm_malloc(8*usize3*task_size);
 RHOKM=(double*)ldm_malloc(8*usize3*task_size);

 e_buf=(double*)ldm_malloc(8*usize3*task_size);     //DBLOC
 f_buf=(double*)ldm_malloc(8*usize3*task_size);     //DBSFC
 
#ifdef POPCHKJN
 swork11=(double*)ldm_malloc(8*usize3*task_size);
 swork1k=(double*)ldm_malloc(8*usize3*task_size);
 swork1km=(double*)ldm_malloc(8*usize3*task_size);
 swork21=(double*)ldm_malloc(8*usize3*task_size);
 swork2k=(double*)ldm_malloc(8*usize3*task_size);
 for(k=0;k<(usize3*task_size);k++)
    {
    swork11[k]=0.0;
    swork1k[k]=0.0;
    swork1km[k]=0.0;
    swork21[k]=0.0;
    swork2k[k]=0.0;
    }
#endif



#ifdef POPCHKJN
 put_reply=7;
#else
 put_reply=2;
#endif

 for(l=0;l<task_tile;l++)
    {
    offset=task_pos*usize;
    offset1=task_pos*usize1;
    offset3=task_pos*usize3;


    get_reply=0;
    offsetg=task_pos;
    athread_get(PE_MODE,TRCR+offsetg,a_buf,8*usize*task_size,(void*)&get_reply,0,8*(ny_block*nx_block-task_size*usize1),8*task_size*usize1);
//  a_buf saved as (nt,km,task_size,usize1) with nt=2, indexed by ((j*km+k)*task_size+i)*usize1+nm,
//     while used as (task_size,usize1,km,nt), indexed by ((i*usize1+nm)*km+k)*nt+j.
    athread_get(PE_MODE,KMT+offset1,d_buf,SIZEINT*usize1*task_size,(void*)&get_reply,0,0,0);
    while(get_reply!=2);


    for(i=0;i<task_size;i++)
       {

       for(nm=0;nm<usize1;nm++)     
        {
        klvl  = klvl0;
        kprev = kprev0;
        
        //  TEMPSFC = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)
/*
        if(a_buf[usize*i+2*nm*km]>-c2)
         TEMPSFC[i*usize1+nm]=a_buf[usize*i+2*nm*km];
*/
        jtmp=i*usize1+nm;              // t with j=0, k=0
        if(a_buf[jtmp]>-c2)
         TEMPSFC[i*usize1+nm]=a_buf[jtmp];
        else
         TEMPSFC[i*usize1+nm]=-c2;

        TEMPKLP[kprev+2*(i*usize1+nm)]=TEMPSFC[i*usize1+nm];  //TEMPK(:,:,kprev) = TEMPSFC

        for(k=1;k<km;k++)       //do k=2,km in fortran
         {
         nk=k+nm*km;
         itmp=(k*task_size+i)*usize1+nm;   // transpose to shape (km,task_size,usize1)
         
#if defined(POPTESTJN)
         if(num_calls==1 && j==0 && i==0 && myproc==0 && myid==0)
             cpe_printf("pressz @ %d: %.18lf \n",k,pressz[k]);
#endif

      
         // TEMPK(:,:,klvl) = merge(-c2,TRCR(:,:,k,1),TRCR(:,:,k,1) < -c2)
/*
         if(a_buf[usize*i+2*nk]>-c2) 
            TEMPKLP[klvl+2*(i*usize1+nm)]=a_buf[usize*i+2*nk];
*/
         jtmp=(k*task_size+i)*usize1+nm;              // t with j=0
         if(a_buf[jtmp]>-c2) 
            TEMPKLP[klvl+2*(i*usize1+nm)]=a_buf[jtmp];
         else
            TEMPKLP[klvl+2*(i*usize1+nm)]=-c2; 


         // call state(k, k, TEMPSFC,  TRCR(:,:,1  ,2),  this_block, RHOFULL=RHO1)
	 TEMPK = TEMPSFC[i*usize1+nm];       
         if(TEMPK>tmax[k]) TQ=tmax[k];
         else if(TEMPK<tmin[k]) TQ=tmin[k];
         else TQ=TEMPK;      

//	 SALTK = a_buf[1+usize*i+2*nm*km];
         SALTK = a_buf[(km*task_size+i)*usize1+nm];   // s with j=1, k=0
	 if(SALTK>smax[k]) SQ=smax[k];
         else if(SALTK<smin[k]) SQ=smin[k];
         else SQ=SALTK;      

         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[k] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[k] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[k] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));


         DENOMK = c1/WORK2;
//         RHO1[nk+usize3*i] = WORK1*DENOMK;
         RHO1[itmp] = WORK1*DENOMK;

#ifdef POPCHKJN
         swork11[nk+usize3*i]=RHO1[itmp];
#endif

         // call state(k, k, TEMPK(:,:,kprev), TRCR(:,:,k-1,2), this_block, RHOFULL=RHOKM)
	 TEMPK = TEMPKLP[kprev+2*(i*usize1+nm)];       
         if(TEMPK>tmax[k]) TQ=tmax[k];
         else if(TEMPK<tmin[k]) TQ=tmin[k];
         else TQ=TEMPK;      

//         SALTK = a_buf[1+2*(nk-1)+usize*i];
         SALTK = a_buf[((km+k-1)*task_size+i)*usize1+nm];    // s with j=1
	 if(SALTK>smax[k]) SQ=smax[k];
         else if(SALTK<smin[k]) SQ=smin[k];
         else SQ=SALTK;      

         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[k] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[k] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[k] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));

         DENOMK = c1/WORK2;
         RHOKM[itmp] = WORK1*DENOMK;

#ifdef POPCHKJN
         swork1km[nk+usize3*i]=RHOKM[itmp];
#endif




         // call state(k, k, TEMPK(:,:,klvl),  TRCR(:,:,k  ,2), this_block, RHOFULL=RHOK)
	 TEMPK = TEMPKLP[klvl+2*(i*usize1+nm)];       
         if(TEMPK>tmax[k]) TQ=tmax[k];
         else if(TEMPK<tmin[k]) TQ=tmin[k];
         else TQ=TEMPK;      

//         SALTK = a_buf[1+2*nk+usize*i];
         SALTK = a_buf[((km+k)*task_size+i)*usize1+nm];         // s with j=1
         if(SALTK>smax[k]) SQ=smax[k];
         else if(SALTK<smin[k]) SQ=smin[k];
         else SQ=SALTK;      


         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[k] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[k] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[k] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));

         DENOMK = c1/WORK2;
         RHOK[itmp] = WORK1*DENOMK;

#ifdef POPCHKJN
         swork1k[nk+usize3*i]=RHOK[itmp];
#endif

         ktmp  = klvl;
         klvl  = kprev;
         kprev = ktmp;

        }        // end of for(k)

       }       // end of for(nm)

    }    // end of for(i)



    // put with stride at slave: 
#ifdef POPCHKJN
    while(put_reply!=7);
#else
    while(put_reply!=2);
#endif



    for(i=0;i<task_size;i++)
       for(nm=0;nm<usize1;nm++)
          f_buf[i*usize1+nm]=c0;      //DBSFC(:,:,1) = c0 
    for(k=1;k<km;k++)
       for(i=0;i<task_size;i++)
          for(nm=0;nm<usize1;nm++)
             {
             itmp =(k    *task_size+i)*usize1+nm; 
             i1tmp=((k-1)*task_size+i)*usize1+nm; 
             if(RHOK[itmp]!=c0)
               {
               f_buf[itmp] =grav*(c1-RHO1[itmp] /RHOK[itmp]);
               e_buf[i1tmp]=grav*(c1-RHOKM[itmp]/RHOK[itmp]);
               }
             else
               {
               f_buf[itmp] =c0;
               e_buf[i1tmp]=c0;
               }
             if((k-1)>=d_buf[i*usize1+nm]) e_buf[i1tmp]=c0;
             }
    for(i=0;i<task_size;i++)
       for(nm=0;nm<usize1;nm++)
          e_buf[((km-1)*task_size+i)*usize1+nm]=c0;          //DBLOC(:,:,km) = c0


#ifdef POPCHKJN
    for(i=0;i<task_size;i++)
    for(nm=0;nm<usize1;nm++)
    for(k=0;k<km;k++)
        {
        swork21[(i*usize1+nm)*km+k]=e_buf[(k*task_size+i)*usize1+nm];
        swork2k[(i*usize1+nm)*km+k]=f_buf[(k*task_size+i)*usize1+nm];
        }
#endif


/*     
     put_reply=0;
     athread_put(PE_MODE,e_buf,DBLOC+offset3,8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,f_buf,DBSFC+offset3,8*usize3*task_size,(void*)&put_reply,0,0);
     while(put_reply!=2);
*/

     // put with stride at slave: 
     put_reply=0;
     offsetp3=task_pos;
     athread_put(PE_MODE,e_buf,DBLOC+offsetp3,8*usize3*task_size,(void*)&put_reply,8*(nx_block*ny_block-task_size*usize1),8*task_size*usize1);
     athread_put(PE_MODE,f_buf,DBSFC+offsetp3,8*usize3*task_size,(void*)&put_reply,8*(nx_block*ny_block-task_size*usize1),8*task_size*usize1);



#ifdef POPCHKJN
     athread_put(PE_MODE,swork11, mwork11+offset3, 8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork1k, mwork1k+offset3, 8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork1km,mwork1km+offset3,8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork21, mwork21+offset3, 8*usize3*task_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork2k, mwork2k+offset3, 8*usize3*task_size,(void*)&put_reply,0,0);
#endif
     
     task_pos = task_pos + task_size;
   }
   // end of for(l)


 // last_size
 if(last_size>0)
   {
   offset=task_pos*usize;
   offset1=task_pos*usize1;
   offset3=task_pos*usize3;


    get_reply=0;
    offsetg=task_pos;
    athread_get(PE_MODE,TRCR+offsetg,a_buf,8*usize*last_size,(void*)&get_reply,0,8*(ny_block*nx_block-last_size*usize1),8*last_size*usize1);
//  a_buf saved as (nt,km,last_size,usize1) with nt=2, indexed by ((j*km+k)*last_size+i)*usize1+nm,
//     while used as (last_size,usize1,km,nt), indexed by ((i*usize1+nm)*km+k)*nt+j.
    athread_get(PE_MODE,KMT+offset1,d_buf,SIZEINT*usize1*last_size,(void*)&get_reply,0,0,0);
    while(get_reply!=2);


    for(i=0;i<last_size;i++)
       {

       for(nm=0;nm<usize1;nm++)     
        {
        klvl  = klvl0;
        kprev = kprev0;
        
        //  TEMPSFC = merge(-c2,TRCR(:,:,1,1),TRCR(:,:,1,1) < -c2)
/*
        if(a_buf[usize*i+2*nm*km]>-c2)
         TEMPSFC[i*usize1+nm]=a_buf[usize*i+2*nm*km];
*/
        jtmp=i*usize1+nm;              // t with j=0, k=0
        if(a_buf[jtmp]>-c2)
         TEMPSFC[i*usize1+nm]=a_buf[jtmp];
        else
         TEMPSFC[i*usize1+nm]=-c2;

        TEMPKLP[kprev+2*(i*usize1+nm)]=TEMPSFC[i*usize1+nm];  //TEMPK(:,:,kprev) = TEMPSFC

        for(k=1;k<km;k++)       //do k=2,km in fortran
         {
         nk=k+nm*km;
         itmp=(k*last_size+i)*usize1+nm;   // transpose to shape (km,last_size,usize1)
 
         // TEMPK(:,:,klvl) = merge(-c2,TRCR(:,:,k,1),TRCR(:,:,k,1) < -c2)
/*
         if(a_buf[usize*i+2*nk]>-c2) 
            TEMPKLP[klvl+2*(i*usize1+nm)]=a_buf[usize*i+2*nk];
*/
         jtmp=(k*last_size+i)*usize1+nm;              // t with j=0
         if(a_buf[jtmp]>-c2) 
            TEMPKLP[klvl+2*(i*usize1+nm)]=a_buf[jtmp];
         else
            TEMPKLP[klvl+2*(i*usize1+nm)]=-c2; 


         // call state(k, k, TEMPSFC,  TRCR(:,:,1  ,2),  this_block, RHOFULL=RHO1)
	 TEMPK = TEMPSFC[i*usize1+nm];       
         if(TEMPK>tmax[k]) TQ=tmax[k];
         else if(TEMPK<tmin[k]) TQ=tmin[k];
         else TQ=TEMPK;      

//	 SALTK = a_buf[1+usize*i+2*nm*km];
         SALTK = a_buf[(km*last_size+i)*usize1+nm];   // s with j=1, k=0
	 if(SALTK>smax[k]) SQ=smax[k];
         else if(SALTK<smin[k]) SQ=smin[k];
         else SQ=SALTK;      

         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[k] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[k] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[k] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));


         DENOMK = c1/WORK2;
         RHO1[itmp] = WORK1*DENOMK;

#ifdef POPCHKJN
         swork11[nk+usize3*i]=RHO1[itmp];
#endif

         // call state(k, k, TEMPK(:,:,kprev), TRCR(:,:,k-1,2), this_block, RHOFULL=RHOKM)
	 TEMPK = TEMPKLP[kprev+2*(i*usize1+nm)];       
         if(TEMPK>tmax[k]) TQ=tmax[k];
         else if(TEMPK<tmin[k]) TQ=tmin[k];
         else TQ=TEMPK;      

//         SALTK = a_buf[1+2*(nk-1)+usize*i];
         SALTK = a_buf[((km+k-1)*last_size+i)*usize1+nm];    // s with j=1
	 if(SALTK>smax[k]) SQ=smax[k];
         else if(SALTK<smin[k]) SQ=smin[k];
         else SQ=SALTK;      

         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[k] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[k] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[k] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));

         DENOMK = c1/WORK2;
         RHOKM[itmp] = WORK1*DENOMK;

#ifdef POPCHKJN
         swork1km[nk+usize3*i]=RHOKM[itmp];
#endif




         // call state(k, k, TEMPK(:,:,klvl),  TRCR(:,:,k  ,2), this_block, RHOFULL=RHOK)
	 TEMPK = TEMPKLP[klvl+2*(i*usize1+nm)];       
         if(TEMPK>tmax[k]) TQ=tmax[k];
         else if(TEMPK<tmin[k]) TQ=tmin[k];
         else TQ=TEMPK;      

//         SALTK = a_buf[1+2*nk+usize*i];
         SALTK = a_buf[((km+k)*last_size+i)*usize1+nm];         // s with j=1
         if(SALTK>smax[k]) SQ=smax[k];
         else if(SALTK<smin[k]) SQ=smin[k];
         else SQ=SALTK;      


         SQ  = c1000*SQ;
         SQR = (double)sqrt((double)SQ);

         //      !*** calculate numerator of MWJF density [P_1(S,T,p)]
         WORK1 = mwjfnums0t0[k] + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2[k] +
                 mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                 mwjfnums1t1 * TQ + mwjfnums2t0 * SQ);

         //      !*** calculate denominator of MWJF density [P_2(S,T,p)]
         WORK2 = mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
              TQ * (mwjfdens0t3[k] + mwjfdens0t4 * TQ))) +
              SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
              SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2));

         DENOMK = c1/WORK2;
         RHOK[itmp] = WORK1*DENOMK;

#ifdef POPCHKJN
         swork1k[nk+usize3*i]=RHOK[itmp];
#endif

         ktmp  = klvl;
         klvl  = kprev;
         kprev = ktmp;

        }        // end of for(k)

       }       // end of for(nm)

    }    // end of for(i)




    // put with stride at slave: 
#ifdef POPCHKJN
    while(put_reply!=7);
#else
    while(put_reply!=2);
#endif



    for(i=0;i<last_size;i++)
       for(nm=0;nm<usize1;nm++)
          f_buf[i*usize1+nm]=c0;      //DBSFC(:,:,1) = c0 
    for(k=1;k<km;k++)
       for(i=0;i<last_size;i++)
          for(nm=0;nm<usize1;nm++)
             {
             itmp =(k    *last_size+i)*usize1+nm; 
             i1tmp=((k-1)*last_size+i)*usize1+nm; 
             if(RHOK[itmp]!=c0)
               {
               f_buf[itmp] =grav*(c1-RHO1[itmp] /RHOK[itmp]);
               e_buf[i1tmp]=grav*(c1-RHOKM[itmp]/RHOK[itmp]);
               }
             else
               {
               f_buf[itmp] =c0;
               e_buf[i1tmp]=c0;
               }
             if((k-1)>=d_buf[i*usize1+nm]) e_buf[i1tmp]=c0;
             }
    for(i=0;i<last_size;i++)
       for(nm=0;nm<usize1;nm++)
          e_buf[((km-1)*last_size+i)*usize1+nm]=c0;          //DBLOC(:,:,km) = c0


#ifdef POPCHKJN
    for(i=0;i<last_size;i++)
    for(nm=0;nm<usize1;nm++)
    for(k=0;k<km;k++)
        {
        swork21[(i*usize1+nm)*km+k]=e_buf[(k*last_size+i)*usize1+nm];
        swork2k[(i*usize1+nm)*km+k]=f_buf[(k*last_size+i)*usize1+nm];
        }
#endif


/*     
     put_reply=0;
     athread_put(PE_MODE,e_buf,DBLOC+offset3,8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,f_buf,DBSFC+offset3,8*usize3*last_size,(void*)&put_reply,0,0);
     while(put_reply!=2);
*/

     // put with stride at slave: 
     put_reply=0;
     offsetp3=task_pos;
     athread_put(PE_MODE,e_buf,DBLOC+offsetp3,8*usize3*last_size,(void*)&put_reply,8*(nx_block*ny_block-last_size*usize1),8*last_size*usize1);
     athread_put(PE_MODE,f_buf,DBSFC+offsetp3,8*usize3*last_size,(void*)&put_reply,8*(nx_block*ny_block-last_size*usize1),8*last_size*usize1);



#ifdef POPCHKJN
     athread_put(PE_MODE,swork11, mwork11+offset3, 8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork1k, mwork1k+offset3, 8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork1km,mwork1km+offset3,8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork21, mwork21+offset3, 8*usize3*last_size,(void*)&put_reply,0,0);
     athread_put(PE_MODE,swork2k, mwork2k+offset3, 8*usize3*last_size,(void*)&put_reply,0,0);
#endif
     

   }
   // end of if(last_size>0)
    

 ldm_free(a_buf,8*usize*task_size);
 ldm_free(d_buf,SIZEINT*usize1*task_size);
 ldm_free(TEMPSFC,8*usize1*task_size);
 ldm_free(TEMPKLP,8*usize1*2*task_size);

 ldm_free(RHO1,8*usize3*task_size);
 ldm_free(RHOK,8*usize3*task_size);
 ldm_free(RHOKM,8*usize3*task_size);

 // put with stride at slave: 
#ifdef POPCHKJN
 while(put_reply!=7);
#else
 while(put_reply!=2);
#endif
 ldm_free(e_buf,8*usize3*task_size);
 ldm_free(f_buf,8*usize3*task_size);


#ifdef POPCHKJN
 ldm_free(swork11,8*usize3*task_size);
 ldm_free(swork1k,8*usize3*task_size);
 ldm_free(swork1km,8*usize3*task_size);
 ldm_free(swork21,8*usize3*task_size);
 ldm_free(swork2k,8*usize3*task_size);
#endif

}
// end of if(task_num>0)


 ldm_free(loc_pts,8*(5*km+1));

 ldm_free(mwjfnums0t0,8*km);
 ldm_free(mwjfnums0t2,8*km);
 ldm_free(mwjfnums1t0,8*km);
 ldm_free(mwjfdens0t0,8*km);
 ldm_free(mwjfdens0t1,8*km);
 ldm_free(mwjfdens0t3,8*km);

}
//------------------------------------------------------------







//#include "math_data.h"
//------------------------------------------------------------------
//void slave_s_kpp_bldepth_loop2(struct param_kpp_bldepth_s1 *st2_master)
void s_kpp_bldepth_loop2(struct param_kpp_bldepth_s1 *st2_master)
{
	
// BOP:

 int locint[SZldmint], *s1di,*s2di;
 double locdbl[SZldmdbl], *s1dd, *s2dd, *s3dd;

// double crt();
 volatile int get_reply, get_reply2, put_reply, put_reply2, retnp,retng;
#ifdef POPCHKJN
 volatile int num_calls=0;
#endif

 double pow_local_data[pow_data_len];
 athread_syn(ARRAY_SCOPE, 0xffff);
 if (_MYID == 0){
   get_reply = 0;
   athread_get(BCAST_MODE, pow_data, pow_local_data, pow_data_bytes, &get_reply, 0xff, 0, 0);
   while (get_reply != 1);
   sync_mem; 
 }
 athread_syn(ARRAY_SCOPE, 0xffff);
 pow_data_local_ptr = pow_local_data;

// Common variables:
 int myproc,myid,myrid,mycid,threads_num;
 int task_num0,task_num,task_numj,task_numj0,task_pos,task_posj,task_post, 
     task_size,task_tile,task_all,last_size,task_tile1,task_size1;
 int nx_block,ny_block,km,nt,lx_block;

 int offset,offset1;
 int max_ldm,mx_usize,usize0,usize,usize1,usize2,usize3;
 int i,j,k,n,m,l,ni,nj,nj1,nk,nm,kl,pij;

 int itmp,jtmp,ktmp,kltmp;


// Input variables: 
 int *KMT; 
 int lshort_wave, partial_bottom_cells, lniw_mixing, linertial;
 double *TRCR, *DBLOC, *DBSFC, *UUU, *VVV, *UCUR, *VCUR, 
     *DZT, *DZU, 
//     *niuel, *nivel,  // used only if lniw_mixing .eq. TRUE!
     *mTrans, *BO, *BOSOL, 
     *m1din, *m2din;
 double Vtc, zeta_s, a_s, c_s, vonkar, eps, eps2, epssfc, concv;


// Input & Output variables:
 int *KBL; 
 double *HBLT, *USTAR, *BFSFC, *STABLE, *ZKL;


// Slave & Temp variables:
 int    kup, kdn, kupper;
 double z_up, z_upper, SIGMA, 
        slope_up, a_co, b_co, c_co, sqrt_arg;
 double *RI_BULK, *WORK, *VSHEAR,  
       *WM, *WS, *B_FRQNCY, *RSH_HBLT,  
#ifdef POPINSJN
      *mwork11,*mwork1k,*mwork1km,
      *mwork21,*mwork2k,*mwork2km,
      *swork11,*swork1k,*swork1km,
      *swork21,*swork2k,*swork2km,
#endif
      tmpw1,tmpw2,ZETA,ZETAH,
//      raup, radn, ratmp, valtmp, inval, retval,    // for pow(tmp,p33) 
      *U1, *V1, *DZU1, *DZT1;
 struct param_kpp_bldepth_s1 st2_slave;

 int *s_KMT, *s_KBL; 
 double *s_DBLOC, *s_DBSFC, *s_UUU, *s_VVV, *s_UCUR, *s_VCUR,  
     *s_DZT, *s_DZU, 
//     *s_niuel, *s_nivel,  
     *s_mTrans, *s_BO, *s_BOSOL,  
     *s_zgrid, *s_Ricr;
 double *s_HBLT, *s_USTAR, *s_BFSFC, *s_ZKL;


// Constants && parameters:
 double 
      c0=0.0,
      c1=1.0,   
      c16=16.0,   
      c1p85=1.85,   
      c2=2.0,
      c200=200.0,
      c2p1=2.1,
      c3=3.0,
      c4=4.0,
      c4p6=4.6,
      c5=5.0,   
      c8=8.0,   
      c10=10.0,   
      c1000=1000.0,   
      c1p5=1.5,
      p001=0.001, 
      p015=0.015, 
      p125=0.125,
      p15=0.15,
      p25=0.25,
      p33=1.0/3.0,
      p5=0.5,
      p54=0.54,
      p85=0.85,
      p9=0.9,
      p909=0.909,
      loceps=1.0e-14;

// EOP.


 get_reply=0;
 athread_get(PE_MODE,st2_master,&st2_slave,sizeof(st2_slave),(void*)&get_reply,0,0,0);

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;

 get_myrid(myrid);
 get_mycid(mycid);

 if(myrid%2>0)		// convient for register communication!
   myid = myrid*8 + (7-mycid);	 
 else
   myid = myrid*8 + mycid; 	 


#ifdef POPINSJN
 num_calls = num_calls+1;      
 if(myid==51 && num_calls<=4)
    cpe_printf("step into s_kpp_bldep_slave: \n"); 
#endif

// sync_mem;
 while(get_reply!=1);
 sync_mem;


 
 myproc      = st2_slave.core_info[0];

 threads_num = st2_slave.core_info[1];


 nx_block = st2_slave.nblock[0];
 ny_block = st2_slave.nblock[1];
 km       = st2_slave.nblock[2];
 nt       = st2_slave.nblock[3];

 
#ifdef POPINSJN
if(myproc==0 && myid==51 && num_calls<=3)
 {
 cpe_printf("get parameter from struct st2_slave: myproc= %d, myid= %d(%d,%d)\n", myproc,myid,myrid,mycid);
 cpe_printf("nx_block= %d, ny_block= %d, km= %d, threads_num=%d.\n", nx_block, ny_block, km, threads_num);
 }
#endif

 

 //max_ldm=get_allocatable_size();

 task_numj0 = ny_block/threads_num;
 if(myid < ny_block%threads_num)
    {
    task_numj=task_numj0+1;
    task_posj=task_numj*myid*nx_block;
    }
 else
    {
    task_numj=task_numj0;
    task_posj=(ny_block%threads_num + task_numj0*myid) *nx_block;
    }
 usize0 = task_numj*nx_block;    // total size for 2d-array.


#ifdef POPINSJN
if(myproc==0 && myid==51 && num_calls<=3)
 {
 cpe_printf("nx_block= %d, ny_block= %d, km= %d, task_numj=%d.\n", nx_block, ny_block, km, task_numj);
 cpe_printf("usize0= %d.\n", usize0);
 }
#endif


 
if(task_numj>0)
{	

 m1din = st2_slave.global_arr[2];
 s1dd = locdbl;			//2*km+2

 get_reply = 0;
 retng = 0;
 retng = athread_get(PE_MODE,m1din,s1dd,8*(2*km+2),(void*)&get_reply,0,0,0); 


 KMT = st2_slave.local_iarr[0];		// int(nx,ny), in
 KBL = st2_slave.local_iarr[1];		// int(nx,ny), out
  
//-------------------------------- 
 lshort_wave          = st2_slave.global_logic[0];
 partial_bottom_cells = st2_slave.global_logic[1];
 lniw_mixing = st2_slave.global_logic[2];
 linertial   = st2_slave.global_logic[3];
 
//--------------------------------
// TRCR    // reserved, double(nx,ny,1,nt), in
 DBLOC   = st2_slave.local_arr[1];      // double(nx,ny,km), in 
 DBSFC   = st2_slave.local_arr[2];      // double(nx,ny,km), in 
 UUU     = st2_slave.local_arr[3];      // double(nx,ny,km), in 
 VVV     = st2_slave.local_arr[4];      // double(nx,ny,km), in 
 
//// STF        //reserved for future athread work. 
//// SHF_QSW 

//--------------------------------
// UCUR    // double(nx,ny), in 
// VCUR    // double(nx,ny), in 
// HBLT    // double(nx,ny), out 
// USTAR   // double(nx,ny), in
// BFSFC   // double(nx,ny), inout

//// STABLE    //reserved for future athread work.
//// SMFT   
 
//// niuel   // double(nx,ny), in, reserved
//// nivel   // double(nx,ny), in
// BO      // double(nx,ny), in
// BOSOL   // double(nx,ny), in
// HLANGM  // double(nx,ny), out for other cases
// ZKL     // double(nx,ny), out for other cases
 
//// BOLUS_SP     // reserved! double(nx,ny), in
// FSTOKES // double(nx,ny), in 

 m2din = st2_slave.local_arr[5];
 
//-------------------------------- 
 DZT   = st2_slave.global_arr[0];     // double(nx,ny,0:km+1), in
 DZU   = st2_slave.global_arr[1];     // double(nx,ny,0:km+1), in
  
//-------------------------------- 
// zgrid   // double(0:km+1), in
// Ricr    // double(km), in

// sync_mem;
 while(get_reply!=1);
 sync_mem;

#ifdef POPINSJN
 if(myproc==0 && myid==51 && num_calls<=3)
    {
    if(retng != 8*(2*km+2))
       cpe_printf("get m1din from master of s_kpp_bldep failed! \n"); 
    cpe_printf("err_s1dd[0]= %lf, err_s1dd[km]= %lf, err_s1dd[2*km+1]= %lf. \n", m1din[0]-s1dd[0],m1din[km]-s1dd[km],m1din[2*km+1]-s1dd[2*km+1]);
    }
#endif

 s_Ricr = s1dd;
 s_zgrid= s1dd+km;     // zgrid(0:km+1)



//--------------------------------  
 mTrans  = st2_slave.global_arr[3];     // double(nx,ny,km), in
 

 Vtc    = st2_slave.local_pts[0];
 a_s    = st2_slave.local_pts[1];
 c_s    = st2_slave.local_pts[2];
 vonkar = st2_slave.local_pts[3];
 eps    = st2_slave.local_pts[4];
 eps2   = st2_slave.local_pts[5];
 epssfc = st2_slave.local_pts[6];
 concv  = st2_slave.local_pts[7];
 zeta_s = st2_slave.local_pts[8];


#ifdef POPINSJN
 mwork11 = st2_slave.sdiag[0];  
 mwork1k = st2_slave.sdiag[1];  
 mwork1km = st2_slave.sdiag[2];  
 mwork21 = st2_slave.sdiag[3];  
 mwork2k = st2_slave.sdiag[4];  
 mwork2km = st2_slave.sdiag[5];  
#endif

 
 mx_usize = MAXUSZBLD;    // MAXUSZBLD SHOULD be larger than usize0! 
 if(km*usize0<=mx_usize)
   {
    ni = 1;
    nj = 1;
    nj1= 1;
    usize  = km*usize0;
    usize1 = (km+2)*usize0;	
   }
 else if(usize0<=mx_usize)  // general case!
   {
    ni = 1;
    nj = km; 
    nj1= km+2;
    usize  = usize0;	
    usize1 = usize0;
   }
 else                   // extreme case
   {
    // MUST to modify! 
    // when km is larger than MAXUSZ.
   }

 
 
 s2dd = locdbl + 4*km; 		//28*usize0

 get_reply = 0;
 retng = 0;
 retng = athread_get(PE_MODE, m2din+task_posj, s2dd, 8*5*usize0, (void*)&get_reply, 0, 8*(ny_block-task_numj)*nx_block, 8*usize0); 

 
 task_num = nj*ni;
// task_num0 = nj1*ni;

 task_tile=task_num/(mx_usize/usize);       // use usize as beseline
 // per tile has (task_num/task_tile) tasks with each task ~(8*usize).
 if(task_tile>0)
   { 
   task_size = mx_usize/usize;
   last_size = task_num%task_size;
   }
 else
   {
   task_size = mx_usize/usize;
   last_size = task_num;
   }  



 if(last_size>0)
   task_tile1=task_tile+1;
 else 
   task_tile1=task_tile;
    
   


 while(get_reply!=1);
 sync_mem;

#ifdef POPINSJN
 if(myproc==0 && myid==51 && num_calls<=3)
    {
    if(retng != 8*5*usize0)
       cpe_printf("get m2din from master of s_kpp_bldep failed! \n"); 
    }
#endif

#ifdef POPINSJN
 if(myproc==0 && myid==51 && num_calls<=3)
    {
    cpe_printf("err_s2dd[0]= %lf, err_s2dd[usize0]= %lf, err_s2dd[2*usize0]= %lf, err_s2dd[3*usize0]= %lf, err_s2dd[4*usize0]= %lf, err_s2dd[5*usize0-1]= %lf. \n", m2din[task_posj+0]-s2dd[0],m2din[task_posj+ny_block*nx_block]-s2dd[usize0],m2din[task_posj+2*ny_block*nx_block]-s2dd[2*usize0],m2din[task_posj+3*ny_block*nx_block]-s2dd[3*usize0],m2din[task_posj+4*ny_block*nx_block]-s2dd[4*usize0],m2din[task_posj+5*ny_block*nx_block-1]-s2dd[5*usize0-1]);
    }
#endif

// get in: 
 s_USTAR = s2dd;
 s_BO    = s2dd + 1*usize0;
 s_BOSOL = s2dd + 2*usize0;
 s_UCUR  = s2dd + 3*usize0;
 s_VCUR  = s2dd + 4*usize0; 
// s_niuel   //reserved! 
// s_nivel   //reserved! 
// s_BOLUS_SP   //reserved! 



 s2di  = locint;			//2*usize0
 s_KMT = s2di;
 s_KBL = s2di + usize0; 

 get_reply = 0;
 retng = 0;
 retng = athread_get(PE_MODE,KMT+task_posj,s_KMT,4*usize0,(void*)&get_reply,0,0,0);




// put out:
 s_HBLT  = s2dd + 5*usize0;
 s_BFSFC = s2dd + 6*usize0;
 s_ZKL   = s2dd + 7*usize0;

// tmp vars:
 RI_BULK = s2dd + 8*usize0;     // double(nx,ny,3)
 WORK    = s2dd + 11*usize0;     // double(nx,ny)
 VSHEAR  = s2dd + 12*usize0;     // double(nx,ny)
 WM      = s2dd + 13*usize0;     // double(nx,ny)
 WS      = s2dd + 14*usize0;     // double(nx,ny)
 B_FRQNCY= s2dd + 15*usize0;     // double(nx,ny)
 RSH_HBLT= s2dd + 16*usize0;     // double(nx,ny)
 U1      = s2dd + 17*usize0;     // double(nx,ny)
 V1      = s2dd + 18*usize0;     // double(nx,ny) 
 DZU1    = s2dd + 19*usize0;     // double(nx,ny) 
 DZT1    = s2dd + 20*usize0;     // double(nx,ny) 


 
#ifdef POPINSJN
 s3dd = locdbl + 4*km + 21*usize0; //15*task_size*usize1
 offset = SZldmdbl;
 swork11 = s3dd + offset - 1*usize1*task_size;
 swork1k = s3dd + offset - 2*usize1*task_size;
 swork1km= s3dd + offset - 3*usize1*task_size;
 swork21 = s3dd + offset - 4*usize1*task_size;
 swork2k = s3dd + offset - 5*usize1*task_size;
 swork2km= s3dd + offset - 6*usize1*task_size;

#else
 s3dd = locdbl + 4*km + 21*usize0; //9*task_size*usize1

#endif


 s_DBLOC = s3dd;
 s_DBSFC = s3dd +   usize1*task_size;
 s_UUU   = s3dd + 2*usize1*task_size;
 s_VVV   = s3dd + 3*usize1*task_size;
 s_mTrans= s3dd + 4*usize1*task_size;
 s_DZT   = s3dd + 5*usize1*task_size;
 s_DZU   = s3dd + 6*usize1*task_size + 2*usize1;




// sync_mem;
 while(get_reply!=1);
 sync_mem;

#ifdef POPINSJN
 if(myproc==0 && myid==51 && num_calls<=3)
    {
    if(retng != 4*usize0)
       cpe_printf("get KMT from master of s_kpp_bldep failed! \n"); 
    cpe_printf("err_sKMT[0]= %d, err_sKMT[usize0-1]= %d. \n", KMT[task_posj+0]-s_KMT[0],KMT[task_posj+usize0-1]-s_KMT[usize0-1]);
    }
#endif





#ifdef POPINSJN
if(myproc==0 && myid==51 && num_calls<=3)
 cpe_printf("init of s_kpp_bldep ok. \n");
#endif


//-------------------------	  
      
   kupper = 0;  // kupper = 1
   kup = 1;     // kup = 2
   kdn = 2;     // kdn = 3

   z_upper = c0;
   z_up    = s_zgrid[1];    // z_up = zgrid(1) with zgrid(0:km+1)
   
  
    
   for(j=0;j<task_numj;j++)
     for(i=0;i<nx_block;i++)
       {
       pij = j*nx_block+i;
       RI_BULK[(kupper*task_numj+j)*nx_block+i] = c0;
       RI_BULK[(kup   *task_numj+j)*nx_block+i] = c0; 
       RSH_HBLT[pij] = c0;
//     KBL = merge(KMT(:,:,bid), 1, (KMT(:,:,bid) > 1))
       if(s_KMT[pij]>1)
          s_KBL[pij] = s_KMT[pij];
       else
          s_KBL[pij] = 1;
       }   


#ifdef POPINSJN
//#ifdef POPTESTJN
if(myproc==0 && myid==1 && num_calls<=3)
 cpe_printf("set s_KBL of s_kpp_bldep ok: KBL(16)=%d; KBL(48)=%d. \n",s_KBL[16],s_KBL[48]);
#endif




// partial_bottom_cells:  True.


// kpp_loop(19):  3545(us), called times:  1(s);   
   if (partial_bottom_cells>0)
     {
     task_pos = task_posj;
//     kltmp = -1;
     for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
  	task_size1 = task_size;
      else  
        task_size1 = last_size;

      get_reply = 0;
      retng = 0;
      retng = athread_get(PE_MODE,DZT+task_pos,s_DZT,8*usize0*(task_size1+1),(void*)&get_reply,0,8*(ny_block-task_numj)*nx_block,8*usize0);
      while(get_reply!=1);
      sync_mem;

#ifdef POPINSJN
      if(myproc==0 && myid==51 && num_calls<=3 && l%8==0)
        {
        if(retng!=8*usize*(task_size1+1))
           cpe_printf("get DZT in loop19 from master of s_kpp_bldep failed: l=%d \n",l);
        cpe_printf("err_DZT[0]= %lf, err_DZT[usize0-1]= %lf. \n", DZT[task_pos+0]-s_DZT[0],DZT[task_pos+usize0-1]-s_DZT[usize0-1]);
        }
#endif
   		 
      for(kl=0;kl<task_size1;kl++)
       {  
//       kltmp = kltmp + 1;
       kltmp = l*task_size + kl;
       for(j=0;j<task_numj;j++)
         for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
            if(kltmp > 0)
// DZT(:,:,0:km+1,bid)
               s_ZKL[pij] = -s_zgrid[(kltmp+1)-1] + p5*(s_DZT[((kl+1)*task_numj+j)*nx_block+i] + s_DZT[(((kl+1)-1)*task_numj+j)*nx_block+i]);
            else
               s_ZKL[pij] = -s_zgrid[1];    //zgrid(0:km+1)
            }
       for(j=0;j<task_numj;j++)
         for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
            if ((kltmp+1) == s_KBL[pij])
               s_HBLT[pij] = s_ZKL[pij];
            }   //end of for(i,j)
       }  //end of for(kl)

      task_pos = task_pos + task_size1*ny_block*nx_block;	
      }   
      //end of for(l)
     }    
     //end of if(partial_bottom_cells>0)


   else    //else of if(partial_bottom_cells>0)
     {
     for(kl=0;kl<km;kl++)
       for(j=0;j<task_numj;j++)
         for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
            s_ZKL[pij] = -s_zgrid[kl+1];
            if((kl+1) == s_KBL[pij])
               s_HBLT[pij] = s_ZKL[pij];
            }
     }    //end of else of if(partial_bottom_cells>0)

#ifdef POPINSJN
if(myproc==0 && myid==1 && num_calls<=3)
 {
 cpe_printf("loop19 of s_kpp_bldep: HBLT(16)=%lf; HBLT(48)=%lf. \n", s_HBLT[15], s_HBLT[47]);
// cpe_printf("task_numj=%d; task_posj=%d; \n", task_numj, task_posj);
// cpe_printf("task_size=%d, last_size=%d, task_tile=%d. \n", task_size, last_size, task_tile);
 }
#endif



//-----------------------------------------------------------------------
//!  compute velocity shear squared on U-grid and use the maximum
//!  of the four surrounding U-grid values for the T-grid.
//!-----------------------------------------------------------------------

//! account for the unresolved part of the NI velocity (equal to the resolved); Markus 9/27/11
//! the 0.05 accounts like in En for the average over the abs(cycle), loosely based on
//! 1/T * integral(cos), also tuned to get En = u'Tau' from model
//! Model resolves half NI energy, so the unresolved part, nivel is added again 
//! x4 because for 16 hours 1/f a time step of 1 hours gets only a quarter of max. increase






#ifdef POPINSJN
if(myproc==0 && myid==51 && num_calls<=3)
 cpe_printf("step into loop21 of s_kpp_bldep: \n");
#endif


//! kpp_loop(21):  114787(us), called times:  1(s);  



//   kltmp = 0;
   task_pos = task_posj;	


   for(l=0;l<task_tile1;l++)
    {
    if(l<task_tile)
       task_size1 = task_size;
    else   
       task_size1 = last_size;
	 
 
    get_reply = 0;
    retng = 0;
    retng = athread_get(PE_MODE,UUU+task_pos,s_UUU,8*usize0*task_size1,(void*)&get_reply,0,8*(ny_block-task_numj)*nx_block,8*usize0);

#ifdef POPINSJN
    if(myproc==0 && myid==51 && num_calls<=3 && l%8==0)
      {
      if(retng != 8*usize*task_size1)
         cpe_printf("get UUU from master of s_kpp_bldep failed! \n");
//       return;
      cpe_printf("err_UUU[0]= %lf, err_UUU[usize0-1]= %lf. \n", UUU[task_pos+0]-s_UUU[0],UUU[task_pos+usize0-1]-s_UUU[usize0-1]);
        }
#endif

    retng = 0;
    retng = athread_get(PE_MODE,VVV+task_pos,s_VVV,8*usize0*task_size1,(void*)&get_reply,0,8*(ny_block-task_numj)*nx_block,8*usize0);

#ifdef POPINSJN
    if(myproc==0 && myid==51 && num_calls<=3 && l%8==0)
        {
        if(retng != 8*usize*task_size1)
          cpe_printf("get VVV from master of s_kpp_bldep failed! \n");
        cpe_printf("err_VVV[0]= %lf, err_VVV[usize0-1]= %lf. \n", VVV[task_pos+0]-s_VVV[0],VVV[task_pos+usize0-1]-s_VVV[usize0-1]);
        }
#endif

    retng = 0;
    retng = athread_get(PE_MODE,DZU+task_pos,s_DZU,8*usize0*(1+task_size1+1),(void*)&get_reply,0,8*(ny_block-task_numj)*nx_block,8*usize0);

#ifdef POPINSJN
    if(myproc==0 && myid==51 && num_calls<=3 && l%8==0)
        {
        if(retng != 8*usize*(task_size1+1))
          cpe_printf("get DZU from master of s_kpp_bldep failed! \n");
        cpe_printf("err_DZU[0]= %lf, err_DZU[usize0-1]= %lf. \n", DZU[task_pos+0]-s_DZU[0],DZU[task_pos+usize0-1]-s_DZU[usize0-1]);
        }
#endif

    retng = 0;
    retng = athread_get(PE_MODE,DZT+task_pos,s_DZT,8*usize0*(1+task_size1+1),(void*)&get_reply,0,8*(ny_block-task_numj)*nx_block,8*usize0);	

#ifdef POPINSJN
    if(myproc==0 && myid==51 && num_calls<=3 && l%8==0)
        {
        if(retng != 8*usize*(task_size1+1))
          cpe_printf("get DZT from master of s_kpp_bldep failed! \n");
        cpe_printf("err_DZT[0]= %lf, err_DZT[usize0-1]= %lf. \n", DZT[task_pos+0]-s_DZT[0],DZT[task_pos+usize0-1]-s_DZT[usize0-1]);
        }
#endif


    get_reply2 = 0;	
    retng = 0;
    retng = athread_get(PE_MODE,DBSFC+task_pos,s_DBSFC,8*usize0*task_size1,(void*)&get_reply2,0,8*(ny_block-task_numj)*nx_block,8*usize0);

#ifdef POPINSJN
    if(myproc==0 && myid==51 && num_calls<=3 && l%8==0)
        {
        if(retng != 8*usize*task_size1)
          cpe_printf("get DBSFC from master of s_kpp_bldep failed! \n");
        cpe_printf("err_DBSFC[0]= %lf, err_DBSFC[usize0-1]= %lf. \n", DBSFC[task_pos+0]-s_DBSFC[0],DBSFC[task_pos+usize0-1]-s_DBSFC[usize0-1]);
        }
#endif

    retng = 0;
    retng = athread_get(PE_MODE,DBLOC+task_pos,s_DBLOC,8*usize0*task_size1,(void*)&get_reply2,0,8*(ny_block-task_numj)*nx_block,8*usize0);

#ifdef POPINSJN
    if(myproc==0 && myid==51 && num_calls<=3 && l%8==0)
        {
        if(retng != 8*usize*task_size1)
          cpe_printf("get DBLOC from master of s_kpp_bldep failed! \n");
        cpe_printf("err_DBLOC[0]= %lf, err_DBLOC[usize0-1]= %lf. \n", DBLOC[task_pos+0]-s_DBLOC[0],DBLOC[task_pos+usize0-1]-s_DBLOC[usize0-1]);
        }
#endif

    retng = 0;
    retng = athread_get(PE_MODE,mTrans+task_pos,s_mTrans,8*usize0*task_size1,(void*)&get_reply2,0,8*(ny_block-task_numj)*nx_block,8*usize0);	

#ifdef POPINSJN
    if(myproc==0 && myid==51 && num_calls<=3 && l%8==0)
        {
        if(retng != 8*usize*task_size1)
          cpe_printf("get mTrans from master of s_kpp_bldep failed! \n");
        cpe_printf("err_mTrans[0]= %lf, err_ mTrans[usize0-1]= %lf. \n", mTrans[task_pos+0]-s_mTrans[0],mTrans[task_pos+usize0-1]-s_mTrans[usize0-1]);
        }
#endif



//    sync_mem;
    while(get_reply!=4);
    sync_mem;


//    sync_mem;
    while(get_reply2!=3);
    sync_mem;



	
    for(kl=0;kl<task_size1;kl++)
     {  
//     kltmp = kltmp + 1;
     kltmp = l*task_size + kl;
	
     if(kltmp>0) 
     //  do kl = 2,km
      { 
      
#ifdef POPINSJN
if(myproc==0 && myid%8==0 && num_calls<=3 && kltmp%8==0)
  cpe_printf("computing WORK of loop21 of s_kpp_bldep: l=%d, kl=%d \n",l,kl);
#endif

      if (lniw_mixing>0)
        {
/*
        for(j=0;j<task_numj;j++)
          for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
            tmpw1 = U1[pij] + s_niuel[pij] - s_UUU[(kl*task_numj+j)*nx_block+i];
            tmpw2 = V1[pij] + s_nivel[pij] - s_VVV[(kl*task_numj+j)*nx_block+i];
            WORK[pij] = tmpw1*tmpw1 + tmpw2*tmpw2; 
            }
*/

        }

      else
        {
        for(j=0;j<task_numj;j++)
          for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
            tmpw1 = U1[pij] - s_UUU[(kl*task_numj+j)*nx_block+i];
            tmpw2 = V1[pij] - s_VVV[(kl*task_numj+j)*nx_block+i];
            WORK[pij] = tmpw1*tmpw1 + tmpw2*tmpw2; 
            }
        }
     
#ifdef POPINSJN
if(myproc==0 && myid==8 && num_calls<=3 && kltmp%16==0)
 {
 cpe_printf("kl=%d: U1(16) = %lf, UUU(16) =  %lf; U1(36) = %lf, UUU(36) = %lf \n", kltmp, U1[15], s_UUU[15], U1[35], s_UUU[35]);
 cpe_printf("kl=%d: V1(16) = %lf, VVV(16) =  %lf; V1(36) = %lf, VVV(36) = %lf \n", kltmp, V1[15], s_VVV[15], V1[35], s_VVV[35]);
 cpe_printf("kl=%d: WORK(16) =  %lf; WORK(36) = %lf \n", kltmp, WORK[15], WORK[35]);
 }
#endif


      if(partial_bottom_cells>0)
        {
        for(j=0;j<task_numj;j++)
          for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
// DZU(:,:,0:km+1,bid)
            tmpw1 = -s_zgrid[(kltmp+1)-1] + p5*(s_DZU[((kl+1)*task_numj+j)*nx_block+i] + s_DZU[(((kl+1)-1)*task_numj+j)*nx_block+i] - DZU1[pij]);
            WORK[pij] = WORK[pij]/(tmpw1*tmpw1); 
// DZT(:,:,0:km+1,bid)
            s_ZKL[pij] = -s_zgrid[(kltmp+1)-1] + p5*(s_DZT[((kl+1)*task_numj+j)*nx_block+i] + s_DZT[(((kl+1)-1)*task_numj+j)*nx_block+i]);
            }
        }
      else
        {
        for(j=0;j<task_numj;j++)
          for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
            s_ZKL[pij] = -s_zgrid[kltmp+1];
            }
        }

#ifdef POPINSJN
if(myproc==0 && myid==8 && num_calls<=3 && kltmp%16==0)
 {
 cpe_printf("kl=%d: ZKL(16) =  %lf; ZKL(36) = %lf \n", kltmp, s_ZKL[15], s_ZKL[35]);
 cpe_printf("kl=%d: WORK(16) =  %lf; WORK(36) = %lf \n", kltmp, WORK[15], WORK[35]);
 }
#endif


//      if(task_posj==0)           // the global 1st row of vshear.
      if(myid==0)
	{  
        for(i=0;i<nx_block;i++)
          VSHEAR[0*nx_block+i] = c0;
	}

      for(j=0;j<task_numj;j++)
        VSHEAR[j*nx_block+0] = c0;

      m4pmax_vshear_s(VSHEAR,WORK,WM,WS,nx_block,ny_block,task_numj,threads_num,myrid,mycid,myid);	 

	  
//!-----------------------------------------------------------------------
//!     compute bfsfc= Bo + radiative contribution down to hbf * hbl
//!     add epsilon to BFSFC to ensure never = 0
//!-----------------------------------------------------------------------

#ifdef POPINSJN
if(myproc==0 && myid==8 && num_calls<=3 && kltmp%16==0)
 {
 cpe_printf("kl=%d: VSHEAR(0,16)=%lf; VSHEAR(0,24)=%lf; VSHEAR(0,36)=%lf; VSHEAR(0,40)=%lf; VSHEAR(0,48)=%lf; VSHEAR(0,76)=%lf. \n", kltmp, VSHEAR[15], VSHEAR[23], VSHEAR[35], VSHEAR[39], VSHEAR[47], VSHEAR[75]);
 }
#endif



      if(lshort_wave>0)
         {
//         select case (sw_absorption_type)
//
//         case ('chlorophyll')   

//!           call sw_trans_chl(2*kl-1,this_block)	
//! start of sw_trans_chl codes expanded by zengyh 2018-11-01: 
         for(j=0;j<task_numj;j++)
            for(i=0;i<nx_block;i++)
               {
               pij = j*nx_block+i;
//              TRANS[pij] = Tr[(CHLINDX[pij])*(ksol+1)+(2*(kl+1)-1)];
               tmpw1 = s_mTrans[(kl*task_numj+j)*nx_block+i];
//! end of sw_trans_chl codes expanded by zengyh 2018-11-01.
               s_BFSFC[pij] = s_BO[pij] + s_BOSOL[pij]*(c1-tmpw1); 
               }
         }

      else
         {
         for(j=0;j<task_numj;j++)
            for(i=0;i<nx_block;i++)
               {
               pij = j*nx_block+i;
               s_BFSFC[pij] = s_BO[pij];
               }
         }

      for(j=0;j<task_numj;j++)
         for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
            if(s_BFSFC[pij]>=c0)
               {
               s_BFSFC[pij] = s_BFSFC[pij] + eps;
               }
            }


#ifdef POPINSJN
if(myproc==0 && myid==1 && num_calls<=3 && kltmp%16==0)
 {
 cpe_printf("kl=%d: s_BFSFC(0,16)=%lf; BO(0,16)=%lf, BOSOL(0,16)=%lf. \n", kltmp, s_BFSFC[15], s_BO[15], s_BOSOL[15]);
 cpe_printf("kl=%d: s_BFSFC(0,36)=%lf; BO(0,36)=%lf, BOSOL(0,36)=%lf. \n", kltmp, s_BFSFC[35], s_BO[35], s_BOSOL[35]);
 }
#endif


//!-----------------------------------------------------------------------
//!     compute velocity scales at sigma, for hbl = -zgrid(kl)
//!-----------------------------------------------------------------------

      SIGMA = epssfc;

#ifdef POPINSJN
if(myproc==0 && myid==8 && num_calls<=3 && kl<=2)
 cpe_printf("step into wscale of loop21 of s_kpp_bldep: \n");
#endif
 
//!      call wscale(SIGMA, ZKL, USTAR, BFSFC, 2, WM, WS)
 
//! start of wscale codes expanded by zengyh 2018-11-01: 
//! subroutine wscale(SIGMA, HBL, USTAR, BFSFC, m_or_s, WM, WS)
      for(j=0;j<task_numj;j++)
         for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
//!----------------------------------------------------------------------
//!  compute zetah and zeta - surface layer is special case
//!-----------------------------------------------------------------------
            ZETAH = SIGMA*s_ZKL[pij]*vonkar*s_BFSFC[pij];
            tmpw1 = s_USTAR[pij];
            ZETA  = ZETAH/(tmpw1*tmpw1*tmpw1 + eps);

//!-----------------------------------------------------------------------
//!  compute velocity scales for tracers
//!-----------------------------------------------------------------------
            if(ZETA >= c0)
               WS[pij] = vonkar*tmpw1/(c1 + c5*ZETA);
            else if(ZETA >= zeta_s)
               WS[pij] = vonkar*tmpw1*sqrt(c1 - c16*ZETA);
            else
               {
               tmpw2 = a_s*(tmpw1*tmpw1*tmpw1)-c_s*ZETAH;
               WS[pij] = vonkar*pow(tmpw2,p33);
//               WS[pij] = vonkar*crt(tmpw2);
// MUST use ldm_math with pow! 2018-11-29
               }
            }
//!  end of wscale codes expanded by zengyh 2018-11-01. 
  
  
//!-----------------------------------------------------------------------
//!     compute the turbulent shear contribution to RI_BULK and store in WM.
//!-----------------------------------------------------------------------

#ifdef POPINSJN
if(myproc==0 && myid==1 && num_calls<=3 && kltmp%16==0)
 {
 cpe_printf("kl=%d: WS(0,16)=%lf; WS(0,36)=%lf. \n", kltmp, WS[15], WS[35]);
 }
#endif



      if(partial_bottom_cells>0)
         {
         if(kltmp<km-1)
            nk = kl+1;
         else
            nk = kl;

         for(j=0;j<task_numj;j++)
            for(i=0;i<nx_block;i++)
               {
               pij = j*nx_block+i;
               tmpw1 = s_DBLOC[(kl*task_numj+j)*nx_block+i];
// apply macro to replace abs;
// else, use 'fabs';
// otherwise, float exception!
               tmpw1 = tmpw1+fabs(tmpw1)+eps2;

// DZT(:,:,0:km+1,bid)
               tmpw2 = s_DZT[((kl+1)*task_numj+j)*nx_block+i]+s_DZT[((nk+1)*task_numj+j)*nx_block+i];

               B_FRQNCY[pij] = sqrt(tmpw1/tmpw2);
               }
         }

      else
            {
            for(j=0;j<task_numj;j++)
               for(i=0;i<nx_block;i++)
                  {
                  pij = j*nx_block+i;
                  tmpw1 = s_DBLOC[(kl*task_numj+j)*nx_block+i]+fabs(s_DBLOC[(kl*task_numj+j)*nx_block+i])+eps2;
                  tmpw2 = s_zgrid[kltmp+1]-s_zgrid[(kltmp+1)+1];
                  B_FRQNCY[pij] = sqrt(p5*tmpw1/tmpw2);
                  }
            }

#ifdef POPINSJN
if(myproc==1 && myid==1 && num_calls<=3 && kltmp%16==0)
   {
   cpe_printf("kl=%d: B_FRQNCY(0,16)=%lf; B_FRQNCY(0,36)=%lf. \n", kltmp, B_FRQNCY[15], B_FRQNCY[35]);
   }
fflush(0);
#endif



      for(j=0;j<task_numj;j++)
         for(i=0;i<nx_block;i++)
            {
            pij = j*nx_block+i;
            tmpw1 = c2p1 - c200*B_FRQNCY[pij];
            tmpw2 = (Vtc/s_Ricr[kltmp])*locmax(tmpw1,concv);
            WM[pij] = s_ZKL[pij]*WS[pij]*B_FRQNCY[pij]* tmpw2;
            }

//!-----------------------------------------------------------------------
//!     compute bulk Richardson number at new level
//!-----------------------------------------------------------------------

#ifdef POPINSJN
if(myproc==0 && myid==1 && num_calls<=3 && kltmp%16==0)
   {
//   cpe_printf("step into Rich_number of loop21 of s_kpp_bldep: %lf, %lf \n", tmpw1, tmpw2);
   cpe_printf("kl=%d: WM(0,16)=%lf; WM(0,36)=%lf. \n", kltmp, WM[15], WM[35]);
   }
//fflush(0);
#endif


      if (partial_bottom_cells>0) 
         {
         for(j=0;j<task_numj;j++)
            for(i=0;i<nx_block;i++)
               {
               pij = j*nx_block+i;
// DZT(:,:,0:km+1,bid)
               tmpw1 = -s_zgrid[(kltmp+1)-1] + p5*(s_DZT[(((kl+1)-1)*task_numj+j)*nx_block+i]+s_DZT[((kl+1)*task_numj+j)*nx_block+i]-DZT1[pij]);

               if(s_KMT[pij]>=kltmp+1)
                  WORK[pij] = s_DBSFC[(kl*task_numj+j)*nx_block+i]/tmpw1;
               else
                  WORK[pij] = c0;

               WM[pij] = WM[pij]/(tmpw1*tmpw1);

               tmpw2 = -s_zgrid[(kltmp+1)-1] + p5*(s_DZU[((kl+1)*task_numj+j)*nx_block+i]+s_DZU[(((kl+1)-1)*task_numj+j)*nx_block+i]-DZU1[pij]);

//               RI_BULK[(kdn*ny_block+j)*nx_block+i] = WORK[pij]/(VSHEAR[pij]+WM[pij]+eps/(tmpw2*tmpw2));
               tmpw1 = eps/(tmpw2*tmpw2);
               tmpw1 = VSHEAR[pij]+WM[pij]+tmpw1;
               RI_BULK[(kdn*task_numj+j)*nx_block+i] = WORK[pij]/tmpw1;
               }
         }

      else
         {
         for(j=0;j<task_numj;j++)
            for(i=0;i<nx_block;i++)
               {
               pij = j*nx_block+i;
               if(s_KMT[pij]>=kltmp+1)
                  WORK[pij] = (s_zgrid[1]-s_zgrid[kltmp+1])*s_DBSFC[(kl*task_numj+j)*nx_block+i];
               else
                  WORK[pij] = c0;
               }

         if (lniw_mixing>0)
           {
           for(j=0;j<task_numj;j++)
             for(i=0;i<nx_block;i++)
               {
               pij = j*nx_block+i;
               RI_BULK[(kdn*task_numj+j)*nx_block+i] = WORK[pij]/(VSHEAR[pij]+WM[pij]+eps);
               }
           }
         else
           {
           if (linertial>0)
             {
// linertial = .False.
//reserved!
//               for(j=0;j<task_numj;j++)
//               for(i=0;i<nx_block;i++)
//                 RI_BULK[(kdn*task_numj+j)*nx_block+i] = WORK[pij]/(VSHEAR[pij]+WM[pij]+s_USTAR[pij]*s_BOLUS_SP[pij]+eps);

             }
           else
             {
             for(j=0;j<task_numj;j++)
               for(i=0;i<nx_block;i++)
                 {
                 pij = j*nx_block+i;
                 RI_BULK[(kdn*task_numj+j)*nx_block+i] = WORK[pij]/(VSHEAR[pij]+WM[pij]+eps);
                 }
             }
           }

         }

//!-----------------------------------------------------------------------
//!       find hbl where Rib = Ricr. if possible, use a quadratic
//!       interpolation. if not, linearly interpolate. the quadratic
//!       equation coefficients are determined using the slope and
//!       Ri_bulk at z_up and Ri_bulk at zgrid(kl). the slope at
//!       z_up is computed linearly between z_upper and z_up.
//!
//!       compute Langmuir depth always 
//!-----------------------------------------------------------------------

#ifdef POPINSJN
if(myproc==0 && myid==1 && num_calls<=3 && kltmp%16==0)
 {
 cpe_printf("kl=%d: RI_BULK(kdn,0,16)=%lf; RI_BULK(kdn,0,36)=%lf. \n", kltmp, RI_BULK[(kdn*task_numj+0)*nx_block+15], RI_BULK[(kdn*task_numj+0)*nx_block+35]);
 }
#endif



#ifdef POPINSJN
if(myproc==0 && myid%8==0 && num_calls<=3 && kl%5==0)
 cpe_printf("step into Langmuir_depth of loop21 of s_kpp_bldep: \n");
#endif

      for(j=0;j<task_numj;j++)
      for(i=0;i<nx_block;i++)
         { 
         pij = j*nx_block+i;
         if(s_KBL[pij] == s_KMT[pij] && RI_BULK[(kdn*task_numj+j)*nx_block+i] > s_Ricr[kltmp])
            {
            slope_up =  (RI_BULK[(kupper*task_numj+j)*nx_block+i] - RI_BULK[(kup*task_numj+j)*nx_block+i])/(z_up - z_upper);
            a_co = (RI_BULK[(kdn*task_numj+j)*nx_block+i] - RI_BULK[(kup*task_numj+j)*nx_block+i] -
                    slope_up*(s_ZKL[pij] + z_up))/((z_up + s_ZKL[pij])*(z_up + s_ZKL[pij]));
            b_co = slope_up + c2 * a_co * z_up;
            c_co = RI_BULK[(kup*task_numj+j)*nx_block+i] + z_up*(a_co*z_up + slope_up) - s_Ricr[kltmp];
            sqrt_arg = b_co*b_co - c4*a_co*c_co;

            if ((fabs(b_co)>eps && fabs(a_co/b_co)<=eps) || sqrt_arg<=c0)
               s_HBLT[pij] = -z_up + (z_up + s_ZKL[pij]) *
                   (s_Ricr[kltmp]         - RI_BULK[(kup*task_numj+j)*nx_block+i])/
                   (RI_BULK[(kdn*task_numj+j)*nx_block+i] - RI_BULK[(kup*task_numj+j)*nx_block+i]);
            else
               s_HBLT[pij] = (-b_co + sqrt(sqrt_arg))/(c2*a_co);
            
            s_KBL[pij] = kltmp+1;

            RSH_HBLT[pij] = (VSHEAR[pij]*s_Ricr[kltmp]/(s_DBSFC[(kl*task_numj+j)*nx_block+i]+eps))/s_HBLT[pij];
            }
         }

#ifdef POPINSJN
if(myproc==0 && myid==8 && num_calls<=3 && kltmp%16==0)
 {
 cpe_printf("kl=%d: KBL(0,16)=%d; KBL(0,36)=%d. \n", kltmp, s_KBL[15], s_KBL[35]);
 cpe_printf("kl=%d: HBLT(0,16)=%lf; HBLT(0,36)=%lf. \n", kltmp, s_HBLT[15], s_HBLT[35]);
 }
#endif

      //!! swap klevel indices and move to next level
      ktmp   = kupper;
      kupper = kup;
      kup    = kdn;
      kdn    = ktmp;
      z_upper = z_up;
      z_up    = s_zgrid[kltmp+1];
      }

     else       
     // else of if(kltmp>0), i.e., kltmp=0.
      {
      for(j=0;j<task_numj;j++)
       for(i=0;i<nx_block;i++)
	   {
            pij = j*nx_block+i;
            U1[pij] = s_UUU[pij];
            V1[pij] = s_VVV[pij];
            DZU1[pij] = s_DZU[pij+usize0];      //usize0: task_numj*nx_block
            DZT1[pij] = s_DZT[pij+usize0];
	   }
      }

     } 
     // end of for(kl)

    task_pos = task_pos + task_size*ny_block*nx_block;	
   }
    // end of for(l)

//-------------------------
//! end of do loop(21).  



 put_reply = 0;
 retnp = 0;
 retnp = athread_put(PE_MODE, s_KBL, KBL+task_posj, 4*usize0, (void*)&put_reply, 0, 0); 

#ifdef POPINSJN
 if(myproc==0 && myid==51 && num_calls<=3)
    {
    if(retnp != 4*usize0)
       cpe_printf("put KBL to master of s_kpp_bldep failed! \n");
    }
#endif
#ifdef POPTESTJN
 if(myproc==0 && myid==51 && num_calls<=9)
    {
    cpe_printf("err_KBL[0]= %d, err_KBL[usize0-1]= %d. \n", KBL[task_posj+0]-s_KBL[0],KBL[task_posj+usize0-1]-s_KBL[usize0-1]);
    }
#endif



#ifdef POPINSJN
if(myproc==0 && myid==51 && num_calls<=3)
 {
 cpe_printf("HBLT(0,16)=%lf; HBLT(0,36)=%lf. \n", s_HBLT[15], s_HBLT[35]);
 cpe_printf("BFSFC(0,16)=%lf; BFSFC(0,36)=%lf. \n", s_BFSFC[15], s_BFSFC[35]);
 }
#endif



 put_reply2 = 0;
 retnp = 0;
 HBLT = m2din + 5*ny_block*nx_block;
 retnp = athread_put(PE_MODE, s_HBLT, HBLT+task_posj, 8*2*usize0, (void*)&put_reply2, 8*(ny_block-task_numj)*nx_block, 8*usize0); 

#ifdef POPINSJN
 if(myproc==0 && myid==51 && num_calls<=3)
    {
    if(retnp != 8*2*usize0)
       cpe_printf("put HBLT+BFSFC to master of s_kpp_bldep failed! \n");
    }
#endif

#ifdef POPTESTJN
 if(myproc==0 && myid==51 && num_calls<=9)
    {
    cpe_printf("err_HBLT[0]= %lf, err_HBLT[usize0-1]= %lf. \n", HBLT[task_posj+0]-s_HBLT[0],HBLT[task_posj+usize0-1]-s_HBLT[usize0-1]);
    cpe_printf("err_BFSFC[0]= %lf, err_BFSFC[usize0-1]= %lf. \n", BFSFC[task_posj+0]-s_BFSFC[0],BFSFC[task_posj+usize0-1]-s_BFSFC[usize0-1]);
    }
#endif



// sync_mem;
 while(put_reply!=1);
 sync_mem;


// sync_mem;
 while(put_reply2!=1);
 sync_mem;


#ifdef POPINSJN
if(myproc==0 && myid==51 && num_calls<=3)
 cpe_printf("step out of s_kpp_bldep. \n");
#endif

}
// end of if(task_numj>0)
  pow_data_local_ptr = pow_data;
return;
}

