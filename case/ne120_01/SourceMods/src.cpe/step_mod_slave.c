#include<slave.h>
#include"step_mod_struct.h"
#include<stdio.h>
#include"cpe_print.h"
#include<math.h>

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

double *ldm_malloc();

void fun_state_step_mod(struct kpp_ddmix *ss1)
{
     //total 52~53k ldm 
      double  *mwjfnums0t0,  mwjfnums0t1, *mwjfnums0t2, mwjfnums0t3,
              *mwjfnums1t0,  mwjfnums1t1,  mwjfnums2t0,
              *mwjfdens0t0, *mwjfdens0t1,  mwjfdens0t2, *mwjfdens0t3, mwjfdens0t4,
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

     int myid,mpe_id,nx,ny,km,tn,k,quotient,iblock,i,oldtime,curtime;
     int usize,usize1,usize3,usize_tracer,max_usize,j_size,i_size,k_size,task_offset,task_num,task_num1,task,task_size,tail_size;
     int tile,task_tracer,koffset,toffset,rho_koffset,rho_toffset,task_rho;
     volatile unsigned long get_reply, put_reply,get_reply2;
     struct kpp_ddmix local_s1;
     double *RHO=NULL,*TRACER=NULL,*pkg=NULL,*pkg_buf=NULL,*tmax=NULL,*tmin=NULL,*smax=NULL,*smin=NULL,*pressz=NULL,*TRACER2,*RHO2;
     double *p_TRACER=NULL,*p_RHO=NULL,*p_RHO2=NULL,*p_TRACER2;
     double p;
     double TEMPK,SALTK,TQ,SQ,SQR,ZWORK1,ZWORK2,DENOMK; 

     
     get_reply=0;
     athread_get(0,ss1,&local_s1,sizeof(local_s1),(void *)&get_reply,0,0,0);//kpp_ddmix ->local_s1
     myid = athread_get_id(-1);   
     while(get_reply!=1);
     asm volatile("memb");
     nx       = local_s1.param[0];
     ny       = local_s1.param[1];
     tn       = local_s1.param[2];
     km       = local_s1.param[3];
     oldtime  = local_s1.param[4];
     curtime  = local_s1.param[5];
     TRACER   = local_s1.arr[0];
     RHO      = local_s1.arr[1];
     TRACER2  = local_s1.arr[2];
     RHO2     = local_s1.arr[3]; 
     pkg      = local_s1.std[0];
  
     pkg_buf=(double *)ldm_malloc(8*km*5);
     get_reply=0;
     athread_get(0,pkg,pkg_buf,8*km*5,(void *)&get_reply,0,0,0);

     mwjfnums0t0=(double*)ldm_malloc(8*km);
     mwjfnums0t2=(double*)ldm_malloc(8*km);
     mwjfnums1t0=(double*)ldm_malloc(8*km);
     mwjfdens0t0=(double*)ldm_malloc(8*km);
     mwjfdens0t1=(double*)ldm_malloc(8*km);
     mwjfdens0t3=(double*)ldm_malloc(8*km);

     max_usize=1000;    
     
       j_size =1;
       i_size =nx; 
       usize  =2*km;  //t nt km combine   TRACER  //3*2*km  //3*3*km
       usize3 =km;     //t    km combine   RHO  //3*km
       usize1 =1;
                      // km nt t task_size uszie1 
    

    task= ny * i_size * j_size;    
    task_num1  =  task/tn;
    if(myid<task%tn)  
    {
       task_num=task_num1+1; 
       task_offset=task_num*myid;
    }
    else
    {
       task_num = task_num1;
       task_offset= ( task%tn ) + task_num1*myid;
    }  
 
    while(get_reply!=1);
    pressz=pkg_buf;
    tmax  =pkg_buf+km;
    tmin  =pkg_buf+2*km;
    smax  =pkg_buf+3*km;
    smin  =pkg_buf+4*km;

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
            mwjfdens0t1[k] = mwjfdp0s0t1 + p*p*p* mwjfdp3s0t1;
            mwjfdens0t3[k] = mwjfdp0s0t3 + p*p*   mwjfdp2s0t3;
        }        

        quotient =task_num/(max_usize/usize3);      
        if( quotient >0) 
        {
            task_size = max_usize/usize3;
            tail_size = task_num%task_size;
        }
        else  
        {
            task_size = task_num;
            tail_size = task_size;     
        }
     

       p_TRACER= ldm_malloc(8* task_size * usize);
       p_TRACER2= ldm_malloc(8* task_size * usize);
       p_RHO   = ldm_malloc(8* task_size * usize3);
       p_RHO2  = ldm_malloc(8* task_size * usize3);


       for(tile =0;tile <quotient;tile++)  //jmix
       {
            task_tracer = task_offset;
            get_reply=0;
            get_reply2=0;
            athread_get(0,TRACER2 +task_tracer  ,p_TRACER2  ,8*usize  *task_size, (void *)&get_reply2,
                         0,8*(ny*nx-task_size*usize1),8*task_size*usize1);
            athread_get(0,TRACER +task_tracer  ,p_TRACER    ,8*usize  *task_size, (void *)&get_reply,
                         0,8*(ny*nx-task_size*usize1),8*task_size*usize1);
           // while(get_reply!=2);
            for(i=0;i<task_size;i++)
                   for(k=0;k<km;k++)      
                   {

                   koffset    =((k)   *task_size+i)*usize1;      //nt 1 old         // t with j=0
                   toffset    =((km+k)*task_size+i)*usize1; //nt 2 old 
                   rho_koffset=(k     *task_size+i)*usize1;
   
 
                   while(get_reply!=1);
                   TEMPK = p_TRACER[koffset];
                   if(TEMPK>tmax[k]) TQ=tmax[k];
                   else if(TEMPK<tmin[k]) TQ=tmin[k];
                   else TQ=TEMPK;
        	   
                   SALTK = p_TRACER[toffset];         // s with j=1
        	   if(SALTK>smax[k]) SQ=smax[k];
        	   else if(SALTK<smin[k]) SQ=smin[k];
       		   else SQ=SALTK;


 
                   SQ  = c1000*SQ;
                   SQR =(double)sqrt((double)SQ);

                   ZWORK1 = mwjfnums0t0[k]   + TQ * (mwjfnums0t1      + TQ * (mwjfnums0t2[k] +
                            mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                            mwjfnums1t1 * TQ   +       mwjfnums2t0      * SQ);
                   
                   ZWORK2 =      mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
                           TQ * (mwjfdens0t3[k] +       mwjfdens0t4 * TQ))) +
                           SQ * (mwjfdens1t0      + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
                           SQR* (mwjfdensqt0      + TQ*TQ*mwjfdensqt2));
                   DENOMK = c1/ZWORK2;
                  
                   p_RHO[rho_koffset] = ZWORK1*DENOMK;
                
                   while(get_reply2!=1);
                   TEMPK = p_TRACER2[koffset];
                   if(TEMPK>tmax[k]) TQ=tmax[k];
                   else if(TEMPK<tmin[k]) TQ=tmin[k];
                   else TQ=TEMPK;
        	   
                   SALTK = p_TRACER2[toffset];         // s with j=1
        	   if(SALTK>smax[k]) SQ=smax[k];
        	   else if(SALTK<smin[k]) SQ=smin[k];
       		   else SQ=SALTK;

 
                   SQ  = c1000*SQ;
                   SQR =(double)sqrt((double)SQ);

                   ZWORK1 = mwjfnums0t0[k]   + TQ * (mwjfnums0t1      + TQ * (mwjfnums0t2[k] +
                            mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                            mwjfnums1t1 * TQ   +       mwjfnums2t0      * SQ);
                   
                   ZWORK2 =      mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
                           TQ * (mwjfdens0t3[k] +       mwjfdens0t4 * TQ))) +
                           SQ * (mwjfdens1t0      + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
                           SQR* (mwjfdensqt0      + TQ*TQ*mwjfdensqt2));
                   DENOMK = c1/ZWORK2;
                   p_RHO2[rho_koffset] = ZWORK1*DENOMK;
                 }   
            

           put_reply=0;
           
           task_rho = task_offset;
           
           athread_put(0,p_RHO   ,RHO     +task_rho   ,8*usize3*task_size,(void *)&put_reply,8*(nx*ny-task_size*usize1),8*task_size*usize1);
           athread_put(0,p_RHO2  ,RHO2    +task_rho   ,8*usize3*task_size,(void *)&put_reply,8*(nx*ny-task_size*usize1),8*task_size*usize1); 
           task_offset = task_offset + task_size;       
           while(put_reply!=2);
        
        }
        
      if(tail_size >0)
       {
            task_tracer = task_offset;
            get_reply=0;
            get_reply2=0;
            athread_get(0,TRACER2 +task_tracer  ,p_TRACER2  ,8*usize  *tail_size, (void *)&get_reply2,
                        0,8*(ny*nx-tail_size*usize1),8*tail_size*usize1);
            athread_get(0,TRACER +task_tracer  ,p_TRACER  ,8*usize  *tail_size, (void *)&get_reply,
                        0,8*(ny*nx-tail_size*usize1),8*tail_size*usize1);
           // while(get_reply!=2);
            for(i=0;i<tail_size;i++)
                   for(k=0;k<km;k++)      
                   {
                   koffset    =((k)   *tail_size+i)*usize1;      //nt 1 old         // t with j=0
                   toffset    =((km+k)*tail_size+i)*usize1; //nt 2 old
                   rho_koffset=((k)   *tail_size+i)*usize1;
                   while(get_reply!=1);
                   TEMPK = p_TRACER[koffset];
                   if(TEMPK>tmax[k]) TQ=tmax[k];
                   else if(TEMPK<tmin[k]) TQ=tmin[k];
                   else TQ=TEMPK;
        	   
                   SALTK = p_TRACER[toffset];         // s with j=1
        	   if(SALTK>smax[k]) SQ=smax[k];
        	   else if(SALTK<smin[k]) SQ=smin[k];
       		   else SQ=SALTK;

        
 
                   SQ  = c1000*SQ;
                   SQR =(double)sqrt((double)SQ);

                   ZWORK1 = mwjfnums0t0[k]   + TQ * (mwjfnums0t1      + TQ * (mwjfnums0t2[k] +
                            mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                            mwjfnums1t1 * TQ   +       mwjfnums2t0      * SQ);
                   
                   ZWORK2 =      mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
                           TQ * (mwjfdens0t3[k] +       mwjfdens0t4 * TQ))) +
                           SQ * (mwjfdens1t0      + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
                           SQR* (mwjfdensqt0      + TQ*TQ*mwjfdensqt2));
                   DENOMK = c1/ZWORK2;
                   p_RHO[rho_koffset] = ZWORK1*DENOMK;
                

                   while(get_reply2!=1);
                   TEMPK = p_TRACER2[koffset];
                   if(TEMPK>tmax[k]) TQ=tmax[k];
                   else if(TEMPK<tmin[k]) TQ=tmin[k];
                   else TQ=TEMPK;
        	   
                   SALTK = p_TRACER2[toffset];         // s with j=1
        	   if(SALTK>smax[k]) SQ=smax[k];
        	   else if(SALTK<smin[k]) SQ=smin[k];
       		   else SQ=SALTK;

        

 
                   SQ  = c1000*SQ;
                   SQR =(double)sqrt((double)SQ);

                   ZWORK1 = mwjfnums0t0[k]   + TQ * (mwjfnums0t1      + TQ * (mwjfnums0t2[k] +
                            mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                            mwjfnums1t1 * TQ   +       mwjfnums2t0      * SQ);
                   
                   ZWORK2 =      mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
                           TQ * (mwjfdens0t3[k] +       mwjfdens0t4 * TQ))) +
                           SQ * (mwjfdens1t0      + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
                           SQR* (mwjfdensqt0      + TQ*TQ*mwjfdensqt2));
                   DENOMK = c1/ZWORK2;
                   p_RHO2[rho_koffset] = ZWORK1*DENOMK;
                 }   //k

           put_reply=0;
          
           task_rho = task_offset;
          
           athread_put(0,p_RHO,RHO       +task_rho   ,8*usize3*tail_size,(void *)&put_reply,8*(nx*ny-tail_size*usize1),8*tail_size*usize1);
           athread_put(0,p_RHO2,RHO2       +task_rho   ,8*usize3*tail_size,(void *)&put_reply,8*(nx*ny-tail_size*usize1),8*tail_size*usize1); 
         
            while(put_reply!=2); 
        }//if
     ldm_free(p_RHO2       ,8*usize3*task_size);
     ldm_free(p_RHO        ,8*usize3*task_size);
     ldm_free(p_TRACER     ,8*usize*task_size);
     ldm_free(p_TRACER2     ,8*usize*task_size);
     
   }//if task_num>0

 ldm_free(pkg_buf,8*(5*km));
 ldm_free(mwjfnums0t0,8*km);
 ldm_free(mwjfnums0t2,8*km);
 ldm_free(mwjfnums1t0,8*km);
 ldm_free(mwjfdens0t0,8*km);
 ldm_free(mwjfdens0t1,8*km);
 ldm_free(mwjfdens0t3,8*km);
}

void fun_state_baroclinic(struct kpp_ddmix *s1)
{

      double  *mwjfnums0t0,  mwjfnums0t1, *mwjfnums0t2, mwjfnums0t3,
              *mwjfnums1t0,  mwjfnums1t1,  mwjfnums2t0,
              *mwjfdens0t0, *mwjfdens0t1,  mwjfdens0t2, *mwjfdens0t3, mwjfdens0t4,
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

     int myid,nx,ny,km,tn,k,quotient,iblock,i;
     int usize,usize1,usize3,usize_tracer,max_usize,j_size,i_size,k_size,task_offset,task_num,task_num1,task,task_size,tail_size;
     int tile,task_tracer,koffset,toffset,rho_koffset,rho_toffset,task_rho;
     volatile unsigned long get_reply, put_reply;
     struct kpp_ddmix local_s1;
     double *RHO=NULL,*TRACER=NULL,*pkg=NULL,*pkg_buf=NULL,*tmax=NULL,*tmin=NULL,*smax=NULL,*smin=NULL,*pressz=NULL;
     double *p_TRACER=NULL,*p_RHO=NULL;
     double p;
     double TEMPK,SALTK,TQ,SQ,SQR,ZWORK1,ZWORK2,DENOMK; 

     
     get_reply=0;
     athread_get(0,s1,&local_s1,sizeof(local_s1),(void *)&get_reply,0,0,0);//kpp_ddmix ->local_s1
     myid = athread_get_id(-1);   
     while(get_reply!=1);
     asm volatile("memb");
     nx       = local_s1.param[0];
     ny       = local_s1.param[1];
     tn       = local_s1.param[2];
     km       = local_s1.param[3];
     TRACER   = local_s1.arr[0];
     RHO      = local_s1.arr[1];
     pkg      = local_s1.std[0];
  
     pkg_buf=(double *)ldm_malloc(8*km*5);
     get_reply=0;
     athread_get(0,pkg,pkg_buf,8*km*5,(void *)&get_reply,0,0,0);

     mwjfnums0t0=(double*)ldm_malloc(8*km);
     mwjfnums0t2=(double*)ldm_malloc(8*km);
     mwjfnums1t0=(double*)ldm_malloc(8*km);
     mwjfdens0t0=(double*)ldm_malloc(8*km);
     mwjfdens0t1=(double*)ldm_malloc(8*km);
     mwjfdens0t3=(double*)ldm_malloc(8*km);

     max_usize=2000;    
     
       j_size =1;
       i_size =nx; 
       usize  =2*km;  //t nt km combine   TRACER  //3*2*km  //3*3*km
       usize3 =km;     //t    km combine   RHO  //3*km
       usize1 =1;
    

    task= ny * i_size * j_size;    
    task_num1  =  task/tn;
    if(myid<task%tn)  
    {
       task_num=task_num1+1; 
       task_offset=task_num*myid;
    }
    else
    {
       task_num = task_num1;
       task_offset= ( task%tn ) + task_num1*myid;
    }  
 
    while(get_reply!=1);
    pressz=pkg_buf;
    tmax  =pkg_buf+km;
    tmin  =pkg_buf+2*km;
    smax  =pkg_buf+3*km;
    smin  =pkg_buf+4*km;

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
            mwjfdens0t1[k] = mwjfdp0s0t1 + p*p*p* mwjfdp3s0t1;
            mwjfdens0t3[k] = mwjfdp0s0t3 + p*p*   mwjfdp2s0t3;
        }        

        quotient =task_num/(max_usize/usize3);      
        if( quotient >0) 
        {
            task_size = max_usize/usize3;
            tail_size = task_num%task_size;
        }
        else  
        {
            task_size = task_num;
            tail_size = task_size;     
        }
     

       p_TRACER= ldm_malloc(8* task_size * usize);
       p_RHO   = ldm_malloc(8* task_size * usize3);


       for(tile =0;tile <quotient;tile++)  //jmix
       {
            task_tracer = task_offset;
            get_reply=0;
            athread_get(0,TRACER +task_tracer  ,p_TRACER    ,8*usize  *task_size, (void *)&get_reply,
                         0,8*(ny*nx-task_size*usize1),8*task_size*usize1);
         //  while(get_reply!=1);
            for(i=0;i<task_size;i++)
                   for(k=0;k<km;k++)      
                   {

                   koffset    =(k   *task_size+i)*usize1;      //nt 1 old         // t with j=0
                   toffset    =((km+k)*task_size+i)*usize1; //nt 2 old 
                   rho_koffset=(k*task_size+i)*usize1;
 
                   while(get_reply!=1);
                   TEMPK = p_TRACER[koffset];
                   if(TEMPK>tmax[k]) TQ=tmax[k];
                   else if(TEMPK<tmin[k]) TQ=tmin[k];
                   else TQ=TEMPK;
        	   
                   SALTK = p_TRACER[toffset];         // s with j=1
        	   if(SALTK>smax[k]) SQ=smax[k];
        	   else if(SALTK<smin[k]) SQ=smin[k];
       		   else SQ=SALTK;

 
                   SQ  = c1000*SQ;
                   SQR =(double)sqrt((double)SQ);

                   ZWORK1 = mwjfnums0t0[k]   + TQ * (mwjfnums0t1      + TQ * (mwjfnums0t2[k] +
                            mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                            mwjfnums1t1 * TQ   +       mwjfnums2t0      * SQ);
                   
                   ZWORK2 =      mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
                           TQ * (mwjfdens0t3[k] +       mwjfdens0t4 * TQ))) +
                           SQ * (mwjfdens1t0      + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
                           SQR* (mwjfdensqt0      + TQ*TQ*mwjfdensqt2));
                   DENOMK = c1/ZWORK2;
                  
                   p_RHO[rho_koffset] = ZWORK1*DENOMK;
                 }   

           put_reply=0;
           task_rho = task_offset;
           athread_put(0,p_RHO   ,RHO     +task_rho   ,8*usize3*task_size,(void *)&put_reply,8*(nx*ny-task_size*usize1),8*task_size*usize1);        
           task_offset = task_offset + task_size;
           while(put_reply!=1);
        }
        
      if(tail_size >0)
       {
            task_tracer = task_offset;
            get_reply=0;
            athread_get(0,TRACER +task_tracer  ,p_TRACER  ,8*usize  *tail_size, (void *)&get_reply,
                        0,8*(ny*nx-tail_size*usize1),8*tail_size*usize1);
           // while(get_reply!=1);
            for(i=0;i<tail_size;i++)
                   for(k=0;k<km;k++)      
                   {
                   koffset    =(k   *tail_size+i)*usize1;      //nt 1 old         // t with j=0
                   toffset    =((km+k)*tail_size+i)*usize1; //nt 2 old
                   rho_koffset=(k   *tail_size+i)*usize1;
                   
                   while(get_reply!=1);
                   TEMPK = p_TRACER[koffset];
                   if(TEMPK>tmax[k]) TQ=tmax[k];
                   else if(TEMPK<tmin[k]) TQ=tmin[k];
                   else TQ=TEMPK;
        	   
                   SALTK = p_TRACER[toffset];         // s with j=1
        	   if(SALTK>smax[k]) SQ=smax[k];
        	   else if(SALTK<smin[k]) SQ=smin[k];
       		   else SQ=SALTK;

                   SQ  = c1000*SQ;
                   SQR =(double)sqrt((double)SQ);

                   ZWORK1 = mwjfnums0t0[k]   + TQ * (mwjfnums0t1      + TQ * (mwjfnums0t2[k] +
                            mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0[k] +
                            mwjfnums1t1 * TQ   +       mwjfnums2t0      * SQ);
                   
                   ZWORK2 =      mwjfdens0t0[k] + TQ * (mwjfdens0t1[k] + TQ * (mwjfdens0t2 +
                           TQ * (mwjfdens0t3[k] +       mwjfdens0t4 * TQ))) +
                           SQ * (mwjfdens1t0      + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
                           SQR* (mwjfdensqt0      + TQ*TQ*mwjfdensqt2));
                   DENOMK = c1/ZWORK2;
                   p_RHO[rho_koffset] = ZWORK1*DENOMK;
                
                }   //k

           put_reply=0;
           task_rho = task_offset;
           athread_put(0,p_RHO,RHO       +task_rho   ,8*usize3*tail_size,(void *)&put_reply,8*(nx*ny-tail_size*usize1),8*tail_size*usize1);
         
            while(put_reply!=1); 
        }//if
     ldm_free(p_RHO        ,8*usize3*task_size);
     ldm_free(p_TRACER     ,8*usize*task_size);     
   }//if task_num>0

 ldm_free(pkg_buf,8*(5*km));
 ldm_free(mwjfnums0t0,8*km);
 ldm_free(mwjfnums0t2,8*km);
 ldm_free(mwjfnums1t0,8*km);
 ldm_free(mwjfdens0t0,8*km);
 ldm_free(mwjfdens0t1,8*km);
 ldm_free(mwjfdens0t3,8*km);
}


void fun2_state_baroclinic(struct kpp_ddmix *s1)
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

     int myid,nx,ny,km,tn,k,quotient,i,quotient1,task_size1;
     int usize,usize1,usize3,usize_tracer,max_usize,j_size,i_size,task_offset,task_num,task_num1,task,task_size,tail_size;
     int tile,task_tracer,rho_koffset,rho_toffset,task_rho;
     volatile int  get_reply, put_reply;
     struct kpp_ddmix local_s1;
     double *RHO=NULL,*TRACER=NULL,*TRACER2=NULL,tmax,tmin,smax,smin,pressz;
     double *p_TRACER=NULL,*p_RHO=NULL,*p_TRACER2=NULL;
     double p;
     double TEMPK,SALTK,TQ,SQ,SQR,ZWORK1,ZWORK2,DENOMK; 

     double ldm_arr[5000];
  
     get_reply=0;
     athread_get(0,s1,&local_s1,sizeof(local_s1),(void *)&get_reply,0,0,0);//kpp_ddmix ->local_s1
     myid = athread_get_id(-1);   
     while(get_reply!=1);
     asm volatile("memb");

     nx       = local_s1.param[0];
     ny       = local_s1.param[1];
     tn       = local_s1.param[2];
     km       = local_s1.param[3];
    
     TRACER   = local_s1.arr[0];
     RHO      = local_s1.arr[1];
     TRACER2  = local_s1.arr[2]; 

     pressz = local_s1.dparam[0];
     tmax   = local_s1.dparam[1];
     tmin   = local_s1.dparam[2];
     smax   = local_s1.dparam[3];
     smin   = local_s1.dparam[4]; 

     max_usize=1000;    
     
       j_size =1;
       i_size =nx; 
       usize  =1;  //t nt km combine   TRACER  //3*2*km  //3*3*km
       usize3 =1;     //t    km combine   RHO  //3*km
       usize1 =1;
    

    task= ny * i_size * j_size;    
    task_num1  =  task/tn;
    if(myid<task%tn)  
    {
       task_num=task_num1+1; 
       task_offset=task_num*myid;
    }
    else
    {
       task_num = task_num1;
       task_offset= ( task%tn ) + task_num1*myid;
    }  
 

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
        
        p=c10*pressz;

        mwjfnums0t0 = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0);
        mwjfnums0t2 = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2);
        mwjfnums1t0 = mwjfnp0s1t0 + p*mwjfnp1s1t0;

        mwjfdens0t0 = mwjfdp0s0t0 + p*mwjfdp1s0t0;
        mwjfdens0t1 = mwjfdp0s0t1 + p*p*p* mwjfdp3s0t1;
        mwjfdens0t3 = mwjfdp0s0t3 + p*p*   mwjfdp2s0t3;

        quotient =task_num/(max_usize/usize3);      
        if( quotient >0) 
        {
            task_size = max_usize/usize3;
            tail_size = task_num%task_size;
        }
        else  
        {
            task_size = task_num;
            tail_size = task_size;     
        }
     
       p_TRACER2= ldm_arr;
       p_TRACER = p_TRACER2 + task_size*usize;
       p_RHO    = p_TRACER  + task_size*usize;

       if(tail_size>0)
            quotient1 =quotient+1;
       else
            quotient1 =quotient;


       for(tile =0;tile <quotient1;tile++)  //jmix
       {
            if(tile<quotient)
                   task_size1=task_size;
            else
                   task_size1=tail_size;


            task_tracer = task_offset;
            get_reply=0;
            athread_get(0,TRACER  +task_tracer  ,p_TRACER    ,8*usize   *task_size1, (void *)&get_reply,0,0,0);
            athread_get(0,TRACER2 +task_tracer  ,p_TRACER2    ,8*usize  *task_size1, (void *)&get_reply,0,0,0);
  
            while(get_reply!=2);
            asm volatile("memb");
            for(i=0;i<task_size1;i++)
            {
                   TEMPK = p_TRACER[i];
                   if(TEMPK>tmax) TQ=tmax;
                   else if(TEMPK<tmin) TQ=tmin;
                   else TQ=TEMPK;
        	   
                   SALTK = p_TRACER2[i];         // s with j=1
        	   if(SALTK>smax) SQ=smax;
        	   else if(SALTK<smin) SQ=smin;
       		   else SQ=SALTK;

 
                   SQ  = c1000*SQ;
                   SQR =(double)sqrt((double)SQ);

                   ZWORK1 = mwjfnums0t0   + TQ * (mwjfnums0t1      + TQ * (mwjfnums0t2 +
                            mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0 +
                            mwjfnums1t1 * TQ   +       mwjfnums2t0      * SQ);
                   
                   ZWORK2 =      mwjfdens0t0 + TQ * (mwjfdens0t1 + TQ * (mwjfdens0t2 +
                           TQ * (mwjfdens0t3 +       mwjfdens0t4 * TQ))) +
                           SQ * (mwjfdens1t0      + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+
                           SQR* (mwjfdensqt0      + TQ*TQ*mwjfdensqt2));
                   DENOMK = c1/ZWORK2;
                  
                   p_RHO[i] = ZWORK1*DENOMK;
                 }   

               put_reply=0;
               task_rho = task_offset;
               athread_put(0,p_RHO   ,RHO     +task_rho   ,8*usize3*task_size1,(void *)&put_reply,0,0);
               while(put_reply!=1);
               asm volatile("memb");
               task_offset = task_offset + task_size1;        
  
            }//tile
        
   }//if task_num>0
}


