#include <stdio.h>
#include "slave.h"
#include <math.h>
#include "vmix_kpp_gq_struct.h"
#include "cpe_print.h"

#define s_double sizeof(double)
#define s_int sizeof(int)

#define max(x,y) ( x>y?x:y )
#define min(x,y) ( x<y?x:y )

#define get_myrid(row)  \
asm volatile("rcsr   %0, 1" : "=r"(row))
#define get_mycid(col)  \
asm volatile("rcsr   %0, 2" : "=r"(col))
__thread_local int i,j,k,km,nx_block,ny_block,bid,lniw,lr,myid,thread_num,size_i,size_j,iblock;
__thread_local int Ddyn,par,mycid,myrid;
__thread_local double c0,c1,Riinfty,rich_mix,p5,p25;
__thread_local int *KMT,*KMT_slave,*KMT_S;
__thread_local double*AT0,*ATS,*ATW,*ATSW,eps;
__thread_local double*AT0_S,*ATS_S,*ATW_S,*ATSW_S,*FRI1_S;
__thread_local double*WORK0,*FRI,*RI_LOC,*UUU,*VVV,*DZU,*DZT,*VSHEAR,*DBLOC,*zgrid;
__thread_local double*WORK0_slave,*FRI_slave,*RI_LOC_slave,*WORK0_S,*FRI_S,*RI_LOC_S;
__thread_local double *WORK0_1_slave,*WORK0_km_slave,*WORK0_km1_slave,*WORK0_right_slave,
                      *UUU_right_S,*VVV_right_S,*DZT_right_S,*DZU_right_S,*WORK0_left_S;
__thread_local double *VDC,*VISC,*KVMIX,*KVMIX_M,*bckgrnd_vdc,*bckgrnd_vvc,*WORK0,*FRI;
__thread_local double *VDC1_S,*VDC2_S,*VISC_S,*KVMIX_S,*KVMIX_M_S,*bckgrnd_vdc_S,
      *bckgrnd_vvc_S,*WORK0_S,*FRI_S,*UUU_S,*VVV_S,*DZU_S,*DZT_S,*VSHEAR_S,*DBLOC_S,*zgrid_S;
__thread_local int volatile get_reply,put_reply;


double *ldm_malloc();

void ri_loop01_fun(struct ri_Loop01_Para *Loop01_Para)
{
/* int mycid,myrid,i,j,k,nx_block,ny_block,thread_num,size_i,size_j,bid,km,par,Ddyn,
      my_task,myid;
 double p5,c0,eps;
 int *KMT,*KMT_S;
 int volatile get_reply,put_reply;
 double *UUU,*VVV,*FRI,*DZU,*DZT,*VSHEAR,*DBLOC,*zgrid,*WORK0,*RI_LOC,*AT0,*ATS,*ATW,*ATSW;
 double *UUU_S,*UUU_right_S,*VVV_S,*VVV_right_S,*FRI_S,*DZU_S,*DZU_right_S,*DZT_S,
        *DZT_right_S,*VSHEAR_S,*WORK0_left_S,*DBLOC_S,*zgrid_S,*RI_LOC_S,*WORK0_S,
        *AT0_S,*ATS_S,*ATW_S,*ATSW_S;
*/
     get_myrid(myrid);
     get_mycid(mycid);
     if(myrid%2>0)
      myid=myrid*8+(7-mycid);
     else
      myid=myrid*8+mycid;

     nx_block=Loop01_Para->param[0];
     ny_block=Loop01_Para->param[1];
     bid=Loop01_Para->param[2];
     km=Loop01_Para->param[3];
     par=Loop01_Para->param[4];
     Ddyn=Loop01_Para->param[5];
     thread_num=Loop01_Para->param[6];
   //  my_task=Loop01_Para->param[7];

     p5=Loop01_Para->dparam[0];
     c0=Loop01_Para->dparam[1];
     eps=Loop01_Para->dparam[2];

     KMT=Loop01_Para->duse_arr[0];
        
     UUU=Loop01_Para->use_arr[0];
     VVV=Loop01_Para->use_arr[1];
     FRI=Loop01_Para->use_arr[2];
     DZU=Loop01_Para->use_arr[3];
     DZT=Loop01_Para->use_arr[4];
     VSHEAR=Loop01_Para->use_arr[5];
     DBLOC=Loop01_Para->use_arr[6];
     zgrid=Loop01_Para->use_arr[7];
     RI_LOC=Loop01_Para->use_arr[8];
     WORK0=Loop01_Para->use_arr[9];
     AT0=Loop01_Para->use_arr[10];
     ATS=Loop01_Para->use_arr[11];
     ATW=Loop01_Para->use_arr[12];
     ATSW=Loop01_Para->use_arr[13];
       
         
         size_i=nx_block-1+1;
       KMT_S=(int*)ldm_malloc(s_int*size_i);
       UUU_S=ldm_malloc(s_double*size_i);
       UUU_right_S=ldm_malloc(s_double*size_i);
       VVV_S=ldm_malloc(s_double*size_i);
       VVV_right_S=ldm_malloc(s_double*size_i);
       FRI_S=ldm_malloc(s_double*size_i);
       DZU_S=ldm_malloc(s_double*size_i);
       DZU_right_S=ldm_malloc(s_double*size_i);
       DZT_S=ldm_malloc(s_double*size_i);
       DZT_right_S=ldm_malloc(s_double*size_i);
       VSHEAR_S=ldm_malloc(s_double*size_i);
       DBLOC_S=ldm_malloc(s_double*size_i);
       zgrid_S=ldm_malloc(s_double*(km+2));
       RI_LOC_S=ldm_malloc(s_double*size_i);
       WORK0_S=ldm_malloc(s_double*size_i);
       WORK0_left_S=ldm_malloc(s_double*size_i);
       AT0_S=ldm_malloc(s_double*size_i);
       ATS_S=ldm_malloc(s_double*size_i);
       ATW_S=ldm_malloc(s_double*size_i);
       ATSW_S=ldm_malloc(s_double*size_i);
   
       for(i=0;i<nx_block;i++)
	   WORK0_left_S[i]=c0;
   
 for(k=1;k<=km;k++)
{     
 for(j=0;j<ny_block;j++)
{  
  if(j%thread_num == myid) 
 
  {
      get_reply=0;
      athread_get(0,KMT+j*nx_block+(bid-1)*nx_block*ny_block,KMT_S,
                                     s_int*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,zgrid,zgrid_S,s_double*(km+2),(void*)&get_reply,0,0,0);
      athread_get(0,DBLOC+j*nx_block+(k-1)*nx_block*ny_block,
                              DBLOC_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,UUU+j*nx_block+(k-1)*nx_block*ny_block,
                                   UUU_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,UUU+j*nx_block+k*nx_block*ny_block,
                                   UUU_right_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,VVV+j*nx_block+(k-1)*nx_block*ny_block,
                                   VVV_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,VVV+j*nx_block+k*nx_block*ny_block,
                                   VVV_right_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,DZT+j*nx_block+k*nx_block*ny_block+(bid-1)*nx_block*ny_block*(km+2),
                                   DZT_S,s_double*size_i,(void*)&get_reply,0,0,0); 
      athread_get(0,DZT+j*nx_block+(k+1)*nx_block*ny_block+(bid-1)*nx_block*ny_block*(km+2),
                                   DZT_right_S,s_double*size_i,(void*)&get_reply,0,0,0); 
      athread_get(0,DZU+j*nx_block+k*nx_block*ny_block+(bid-1)*nx_block*ny_block*(km+2),
                                   DZU_S,s_double*size_i,(void*)&get_reply,0,0,0); 
      athread_get(0,DZU+j*nx_block+(k+1)*nx_block*ny_block+(bid-1)*nx_block*ny_block*(km+2),
                                   DZU_right_S,s_double*size_i,(void*)&get_reply,0,0,0);
      //athread_get(0,WORK0+j*nx_block+(k-1)*nx_block*ny_block,
       //                            WORK0_left_S,s_double*size_i,(void*)&get_reply,0,0,0);
      
      athread_get(0,AT0+j*nx_block+(bid-1)*nx_block*ny_block,
                                   AT0_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,ATS+j*nx_block+(bid-1)*nx_block*ny_block,
                                   ATS_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,ATW+j*nx_block+(bid-1)*nx_block*ny_block,
                                   ATW_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,ATSW+j*nx_block+(bid-1)*nx_block*ny_block,
                                   ATSW_S,s_double*size_i,(void*)&get_reply,0,0,0);
     
     
      
      while(get_reply!=15); 
       asm volatile("memb");
      for(i=0;i<size_i;i++)
    {
      if (k < km)
       {
         FRI_S[i] = (UUU_S[i]-UUU_right_S[i])*(UUU_S[i]-UUU_right_S[i]) + 
                    (VVV_S[i]-VVV_right_S[i])*(VVV_S[i]-VVV_right_S[i]);

         if (par)
      FRI_S[i] = FRI_S[i]/((p5*(DZU_S[i] + DZU_right_S[i]))*(p5*(DZU_S[i] + DZU_right_S[i])));
       
          if (Ddyn)    //l1Ddyn=F
            VSHEAR_S[i] = FRI_S[i];
             
          else
          { 
        
           if(i==0)
               VSHEAR_S[0]=c0;    
           else
               
      ugrid_to_tgrid_slave(&VSHEAR_S[i],FRI_S[i],FRI_S[i-1],AT0_S[i],ATS_S[i],ATW_S[i],
                         ATSW_S[i],ny_block);
#ifdef men             
                 if(myid==1&&my_task==0&&i==1)
               { 
                 cpe_printf("k=%d\n",k);
               }
#endif           
          }
   
            
        }        
     else

       VSHEAR_S[i] = c0;
       if (par)
         if(k<km)
       RI_LOC_S[i] = DBLOC_S[i]/
     (VSHEAR_S[i] + eps/((p5*(DZT_S[i]+DZT_right_S[i]))*(p5*(DZT_S[i]+DZT_right_S[i]))))/ 
         (p5*(DZT_S[i] + DZT_right_S[i]));
         else
            RI_LOC_S[i] = DBLOC_S[i]/(VSHEAR_S[i] +
                     eps/((p5*DZT_S[i])*(p5*DZT_S[i])))/(p5*DZT_right_S[i]);
      else
         RI_LOC_S[i] = DBLOC_S[i]*(zgrid_S[k]-zgrid_S[k+1])/(VSHEAR_S[i]+ eps);
         
            
     // WORK0(:,:,k)   = merge(RI_LOC, WORK0(:,:,k-1), k <= KMT(:,:,bid))
       if(k<=KMT_S[i])
         WORK0_S[i]=RI_LOC_S[i];
       else
         WORK0_S[i]=WORK0_left_S[i];
	      
	 WORK0_left_S[i]=WORK0_S[i];
         
  }
    
           
       put_reply=0;
    athread_put(0,RI_LOC_S,RI_LOC+j*nx_block,s_double*size_i,(void*)&put_reply,0,0);
   // athread_put(0,VSHEAR_S,VSHEAR+j*nx_block,s_double*size_i,(void*)&put_reply,0,0);
  //  athread_put(0,FRI_S,FRI+j*nx_block,s_double*size_i,(void*)&put_reply,0,0);
    
    athread_put(0,WORK0_S,WORK0+j*nx_block+k*nx_block*ny_block,s_double*size_i,(void*)&put_reply,0,0);
     while(put_reply!=2);
     asm volatile("memb");
  }
 }    
}

       ldm_free(KMT_S,s_int*size_i);
       ldm_free(UUU_S,s_double*size_i);
       ldm_free(UUU_right_S,s_double*size_i);
       ldm_free(VVV_S,s_double*size_i);
       ldm_free(VVV_right_S,s_double*size_i);
       ldm_free(FRI_S,s_double*size_i);
       ldm_free(DZU_S,s_double*size_i);
       ldm_free(DZU_right_S,s_double*size_i);
       ldm_free(DZT_S,s_double*size_i);
       ldm_free(DZT_right_S,s_double*size_i);
       ldm_free(VSHEAR_S,s_double*size_i);
       ldm_free(DBLOC_S,s_double*size_i);
       ldm_free(zgrid_S,s_double*(km+2));
       ldm_free(RI_LOC_S,s_double*size_i);
       ldm_free(WORK0_S,s_double*size_i);
       ldm_free(WORK0_left_S,s_double*size_i);
       ldm_free(AT0_S,s_double*size_i);
       ldm_free(ATS_S,s_double*size_i);
       ldm_free(ATW_S,s_double*size_i);
       ldm_free(ATSW_S,s_double*size_i);
      
}

void ri_loop02_fun(struct ri_Loop02_Para *Loop02_Para)
{

/* int i,j,k,km,nx_block,ny_block,bid,myid,thread_num,size_i,size_j;
 double p5,p25;
 int *KMT,*KMT_slave;
 double *RI_LOC,*FRI,*WORK0;
 double *WORK0_slave,*FRI_slave,*RI_LOC_slave;
 double *WORK0_1_slave,*WORK0_km_slave,*WORK0_km1_slave,*WORK0_right_slave;
 int volatile get_reply,put_reply;
*/

	myid = athread_get_id(-1);

	nx_block=Loop02_Para->param[0];
	ny_block=Loop02_Para->param[1];
	bid=Loop02_Para->param[2];
	k=Loop02_Para->param[3];
	km=Loop02_Para->param[4];
	thread_num=Loop02_Para->param[5];

        p5=Loop02_Para->dparam[0];
        p25=Loop02_Para->dparam[1];
        
        KMT=Loop02_Para->duse_arr[0];

        WORK0=Loop02_Para->use_arr[0];
        FRI=Loop02_Para->use_arr[1];
        RI_LOC=Loop02_Para->use_arr[2];
        size_j=ny_block-1+1;
        size_i=nx_block-1+1;
       KMT_slave=(int*)ldm_malloc(s_int*size_i);   
       WORK0_slave=ldm_malloc(s_double*size_i);
       WORK0_1_slave=ldm_malloc(s_double*size_i);
       WORK0_km_slave=ldm_malloc(s_double*size_i);
       WORK0_km1_slave=ldm_malloc(s_double*size_i);
       WORK0_right_slave=ldm_malloc(s_double*size_i);
       FRI_slave=ldm_malloc(s_double*size_i);
       RI_LOC_slave=ldm_malloc(s_double*size_i);
      
      for(j=0;j<ny_block;j++)
 {  
    if(j%thread_num == myid)
   {   
    
      get_reply=0;
      athread_get(0,WORK0+j*nx_block+1*nx_block*ny_block,WORK0_1_slave,s_double*size_i,
                                               (void*)&get_reply,0,0,0);
      athread_get(0,WORK0+j*nx_block+km*nx_block*ny_block,
                 WORK0_km_slave,s_double*size_i,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      asm volatile("memb"); 
       for(i=0;i<size_i;i++)
      {
      FRI_slave[i]  =  p25 * WORK0_1_slave[i];
      WORK0_km1_slave[i] =WORK0_km_slave[i];
        } 
       put_reply=0;
       athread_put(0,WORK0_km1_slave,WORK0+j*nx_block+(km+1)*nx_block*ny_block,s_double*size_i,                                                                                                                    
                                                                  (void*)&put_reply,0,0);   
      while(put_reply!=1); 
      asm volatile("memb");   
     }
   
 }

 for(k=1;k<=km;k++)
{     
 for(j=0;j<ny_block;j++)
{  
  if(j%thread_num == myid) 
  { 
      get_reply=0;
      athread_get(0,KMT+j*nx_block+(bid-1)*nx_block*ny_block,
                 KMT_slave,s_int*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,WORK0+j*nx_block+k*nx_block*ny_block,
                 WORK0_slave,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,WORK0+j*nx_block+(k+1)*nx_block*ny_block,
                 WORK0_right_slave,s_double*size_i,(void*)&get_reply,0,0,0);        
      while(get_reply!=3); 
      asm volatile("memb");
      for(i=0;i<size_i;i++)
      {
          
          RI_LOC_slave[i] = WORK0_slave[i];
          if (KMT_slave[i]>=3) 
           {
              WORK0_slave[i] = FRI_slave[i] + p5*RI_LOC_slave[i] + p25*WORK0_right_slave[i];
           }
          FRI_slave[i] = p25*RI_LOC_slave[i];
          
        }
    put_reply=0;
   athread_put(0,WORK0_slave,WORK0+j*nx_block+k*nx_block*ny_block,s_double*size_i,
                                                                    (void*)&put_reply,0,0);
   athread_put(0,FRI_slave,FRI+j*nx_block,s_double*size_i,(void*)&put_reply,0,0);
   while(put_reply!=2);
    asm volatile("memb");
   
    }
  }
}      
       ldm_free(KMT_slave,s_int*size_i);
       ldm_free(WORK0_slave,s_double*size_i);
       ldm_free(WORK0_1_slave,s_double*size_i);
       ldm_free(WORK0_km_slave,s_double*size_i);
       ldm_free(WORK0_km1_slave,s_double*size_i);
       ldm_free(WORK0_right_slave,s_double*size_i);
       ldm_free(FRI_slave,s_double*size_i);
       ldm_free(RI_LOC_slave,s_double*size_i);
  
}


void ri_loop03_fun(struct ri_Loop03_Para *Loop03_Para)
{

 /*
 int i,j,k,km,nx_block,ny_block,bid,lniw,lr,myid,thread_num,size_i,size_j;
 double c0,c1,Riinfty,rich_mix;
 double *VDC,*VISC,*KVMIX,*KVMIX_M,*bckgrnd_vdc,*bckgrnd_vvc,*WORK0,*FRI;
 double *VDC1_S,*VDC2_S,*VISC_S,*KVMIX_S,*KVMIX_M_S,*bckgrnd_vdc_S,*bckgrnd_vvc_S,
                      *WORK0_S,*FRI_S;
 int volatile get_reply,put_reply;
 */

	myid = athread_get_id(-1);

	nx_block=Loop03_Para->param[0];
	ny_block=Loop03_Para->param[1];
	bid=Loop03_Para->param[2];
	km=Loop03_Para->param[3];
        lniw=Loop03_Para->param[4];
        lr=Loop03_Para->param[5];
	thread_num=Loop03_Para->param[6];

        c0=Loop03_Para->dparam[0];
        c1=Loop03_Para->dparam[1];
        Riinfty=Loop03_Para->dparam[2];
        rich_mix=Loop03_Para->dparam[3];
        
        VDC=Loop03_Para->use_arr[0];
        VISC=Loop03_Para->use_arr[1];
        KVMIX=Loop03_Para->use_arr[2];
        KVMIX_M=Loop03_Para->use_arr[3];
        bckgrnd_vdc=Loop03_Para->use_arr[4];
        bckgrnd_vvc=Loop03_Para->use_arr[5];
        WORK0=Loop03_Para->use_arr[6];
        FRI=Loop03_Para->use_arr[7];
         
         size_i=nx_block-1+1;
       VDC1_S=ldm_malloc(s_double*size_i);
       VDC2_S=ldm_malloc(s_double*size_i);
       VISC_S=ldm_malloc(s_double*size_i);
       KVMIX_S=ldm_malloc(s_double*size_i);
       KVMIX_M_S=ldm_malloc(s_double*size_i);
       bckgrnd_vdc_S=ldm_malloc(s_double*size_i);
       bckgrnd_vvc_S=ldm_malloc(s_double*size_i);
       WORK0_S=ldm_malloc(s_double*size_i);
       FRI_S=ldm_malloc(s_double*size_i);

 for(k=1;k<=km;k++)
{     
 for(j=0;j<ny_block;j++)
{  
  if(j%thread_num == myid) 
  {
      get_reply=0;
      athread_get(0,VDC+j*nx_block+k*nx_block*ny_block+(2-1)*nx_block*ny_block*(km+2),
                              VDC2_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,VISC+j*nx_block+k*nx_block*ny_block,VISC_S,s_double*size_i,
                              (void*)&get_reply,0,0,0);
      athread_get(0,bckgrnd_vvc+j*nx_block+(k-1)*nx_block*ny_block+(bid-1)*nx_block*ny_block*km,
                              bckgrnd_vvc_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,bckgrnd_vdc+j*nx_block+(k-1)*nx_block*ny_block+(bid-1)*nx_block*ny_block*km,
                              bckgrnd_vdc_S,s_double*size_i,(void*)&get_reply,0,0,0);
      athread_get(0,WORK0+j*nx_block+k*nx_block*ny_block,WORK0_S,s_double*size_i,(void*)&get_reply,0,0,0);
      while(get_reply!=5);
      asm volatile("memb"); 
      for(i=0;i<size_i;i++)
  {
       if(k<km) 
       {  
       
          if(lniw)                        //lniw=F
           {
              KVMIX_S[i]   = VDC2_S[i];
              KVMIX_M_S[i] = VISC_S[i];
           }
          else 
           {
              KVMIX_S[i]   = bckgrnd_vdc_S[i];
              KVMIX_M_S[i] = bckgrnd_vvc_S[i];
           }
       }

      if (lniw)                          //lniw=F
      {
        if (lr)
        { 
            
          FRI_S[i]      = min((max(WORK0_S[i],c0))/Riinfty, c1);

          VISC_S[i]     = VISC_S[i] +
                  rich_mix*(c1 - FRI_S[i]*FRI_S[i])*(c1 - FRI_S[i]*FRI_S[i])*(c1 - FRI_S[i]*FRI_S[i]);
         }
        if ( k < km ) 
         { 
          VDC2_S[i]     = VDC2_S[i] + 
                  rich_mix*(c1 - FRI_S[i]*FRI_S[i])*(c1 - FRI_S[i]*FRI_S[i])*(c1 - FRI_S[i]*FRI_S[i]);
       
          VDC1_S[i]     = VDC2_S[i];
         }
      }


      else
     {
        if (lr)                             //lrich=T
       {
          {  
            FRI_S[i]   = min((max(WORK0_S[i],c0))/Riinfty, c1);

            VISC_S[i] =  bckgrnd_vvc_S[i] + 
               rich_mix*(c1 - FRI_S[i]*FRI_S[i])*(c1 - FRI_S[i]*FRI_S[i])*(c1 - FRI_S[i]*FRI_S[i]);
          }
           if ( k < km ) 
          { 
              VDC2_S[i] = bckgrnd_vdc_S[i] + 
                  rich_mix*(c1 - FRI_S[i]*FRI_S[i])*(c1 - FRI_S[i]*FRI_S[i])*(c1 - FRI_S[i]*FRI_S[i]);
              VDC1_S[i] = VDC2_S[i];
          } 
       } 
       
       else
        {
           VISC_S[i] = bckgrnd_vvc_S[i];
            
           if ( k < km ) 
            {
              VDC2_S[i] = bckgrnd_vdc_S[i];
              VDC1_S[i] = VDC2_S[i];
           }
       }
     }
}
    
           
    put_reply=0;
    
   // athread_put(0,FRI_S,FRI+j*nx_block,s_double*size_i,(void*)&put_reply,0,0);
    athread_put(0,KVMIX_S,KVMIX+j*nx_block,s_double*size_i,(void*)&put_reply,0,0);
    athread_put(0,KVMIX_M_S,KVMIX_M+j*nx_block,s_double*size_i,(void*)&put_reply,0,0);
    athread_put(0,VISC_S,VISC+j*nx_block+k*nx_block*ny_block,s_double*size_i,(void*)&put_reply,0,0);
    athread_put(0,VDC1_S,VDC+j*nx_block+k*nx_block*ny_block+(1-1)*nx_block*ny_block*(km+2),
                                     s_double*size_i,(void*)&put_reply,0,0);
    athread_put(0,VDC2_S,VDC+j*nx_block+k*nx_block*ny_block+(2-1)*nx_block*ny_block*(km+2),
                                   s_double*size_i,(void*)&put_reply,0,0);
     while(put_reply!=5);
     asm volatile("memb");
   
    }
  }
}

       ldm_free(VDC1_S,s_double*size_i);
       ldm_free(VDC2_S,s_double*size_i);
       ldm_free(VISC_S,s_double*size_i);
       ldm_free(KVMIX_S,s_double*size_i);
       ldm_free(KVMIX_M_S,s_double*size_i);
       ldm_free(bckgrnd_vdc_S,s_double*size_i);
       ldm_free(bckgrnd_vvc_S,s_double*size_i);
       ldm_free(WORK0_S,s_double*size_i);
       ldm_free(FRI_S,s_double*size_i);
       
  
}
