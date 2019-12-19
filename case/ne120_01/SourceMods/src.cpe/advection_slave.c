
#include <stdio.h>
#include "slave.h"
#include "advection_struct.h"

#define TYI 400 
#define TYR 5600

#define s_double sizeof(double)
#define s_int sizeof(int)   

#define locmin(x,y) ((x)<(y)?(x) : (y)) 
#define locmax(x,y) ((x)>(y)?(x) : (y)) 



void slave_adhdt_loop11_fun();
void loop11_advt();
void loop16_hdifft();
void slave_adhdt_loop13_fun();
void loop13_advt();
void loop17_hdifft();



   void slave_adhdt_loop11_fun(struct adhdt_Loop11_Para* Loop11_Para)   
{
    int phyid;
	
    phyid = athread_get_id(-1);

    if(phyid<24)
	loop11_advt(Loop11_Para);

    else if(phyid<48)
	loop16_hdifft(Loop11_Para);
	
    return;
}
  
  
void loop11_advt(struct adhdt_Loop11_Para *Loop11_Para)
{
   int volatile get_reply,put_reply;
   int phyid,stid,myid,jbeg,jend,ibeg,iend,bid,nx_block,ny_block,k,thread_num,km,my_task;
   int advtk[TYI], *KMT;   
   double advt_11[TYR];
   double p5=0.5, c0=0.0;
   double *UUU,*VVV,*DZU,*DXU,*DYU,*UTE,*UTW,*VTN,*VTS,*WTK,*WTKB,*TAREA_R;
   double *UUU_slave,*VVV_slave,*DZU_slave,*DXU_slave,*DYU_slave,*UTE_slave,
          *UTW_slave,*VTN_slave,*VTS_slave,*WTK_slave,*WTKB_slave,
          *FC_slave,*TAREA_R_slave;
   int i,j,nxy,*KMT_slave;
   struct adhdt_Loop11_Para Loop11_Para_S;
   
   get_reply=0;
   athread_get(PE_MODE,Loop11_Para,&Loop11_Para_S,sizeof(Loop11_Para_S),
                                              (void*)&get_reply,0,0,0);
   asm volatile("memb"); 
   while(get_reply!=1);
 
    phyid = athread_get_id(-1);
    stid  = Loop11_Para_S.param[18];
    myid = phyid - stid;

    thread_num=Loop11_Para_S.param[17];	
    if(myid>=thread_num) return;
	
    jbeg = Loop11_Para_S.param[0];
    jend = Loop11_Para_S.param[1];
    ibeg = Loop11_Para_S.param[2];
    iend = Loop11_Para_S.param[3];
    nx_block = Loop11_Para_S.param[4];
    ny_block = Loop11_Para_S.param[5];
    bid  = Loop11_Para_S.param[7];
    k    = Loop11_Para_S.param[8];
    km   = Loop11_Para_S.param[9];

    my_task=Loop11_Para_S.param[10];
	
    nxy = nx_block*ny_block;

    KMT = Loop11_Para_S.duse_arr[0];
 
   
    UUU = Loop11_Para_S.use_arr[0];
    VVV = Loop11_Para_S.use_arr[1];
    DYU = Loop11_Para_S.use_arr[2];
    DXU = Loop11_Para_S.use_arr[3];
    DZU = Loop11_Para_S.use_arr[4];

    VTN = Loop11_Para_S.use_arr[5];       //UVNSEW
    VTS = Loop11_Para_S.use_arr[5]+  nxy;
    UTE = Loop11_Para_S.use_arr[5]+2*nxy;
    UTW = Loop11_Para_S.use_arr[5]+3*nxy;
	
    WTK = Loop11_Para_S.use_arr[12];
    WTKB    = Loop11_Para_S.use_arr[13];
    TAREA_R = Loop11_Para_S.use_arr[14];
	
	

    KMT_slave = advtk;
	
    UUU_slave= advt_11;
    VVV_slave= advt_11+2*nx_block;
    DYU_slave= advt_11+4*nx_block;
    DXU_slave= advt_11+6*nx_block;
    DZU_slave= advt_11+8*nx_block;

    VTN_slave= advt_11+10*nx_block;
    VTS_slave= advt_11+11*nx_block;  
    UTE_slave= advt_11+12*nx_block;
    UTW_slave= advt_11+13*nx_block;

    WTK_slave= advt_11+15*nx_block;
    WTKB_slave= advt_11+16*nx_block;  
    FC_slave = advt_11+17*nx_block;
    TAREA_R_slave = advt_11+18*nx_block;

 
 for(j=myid+1;j<=ny_block;j+=thread_num)
    {
    if(j>=jbeg && j<=jend)
      {   
      get_reply=0;
      
      athread_get(0,DYU+(j-2)*nx_block,DYU_slave,s_double*2*nx_block,(void*)&get_reply,0,0,0);
      athread_get(0,DXU+(j-2)*nx_block,DXU_slave,s_double*2*nx_block,(void*)&get_reply,0,0,0);
      athread_get(0,DZU+(j-2)*nx_block,DZU_slave,s_double*2*nx_block,(void*)&get_reply,0,0,0);  
      athread_get(0,UUU+(j-2)*nx_block,UUU_slave,s_double*2*nx_block,(void*)&get_reply,0,0,0);
      athread_get(0,VVV+(j-2)*nx_block,VVV_slave,s_double*2*nx_block,(void*)&get_reply,0,0,0);   
      while(get_reply!=5); 
      asm volatile("memb");
        
      for(i=0;i<ibeg-1;i++)
        {
        UTE_slave[i]=c0;
        UTW_slave[i]=c0;
        VTN_slave[i]=c0;
        VTS_slave[i]=c0;
        }		   

      for(i=ibeg-1;i<iend;i++)
        {
        UTE_slave[i] = p5*(UUU_slave[nx_block+i]  * DYU_slave[nx_block+i]*  
                                                    DZU_slave[nx_block+i] + 
                           UUU_slave[i]*            DYU_slave[i]*  
                                                    DZU_slave[i]);
        UTW_slave[i] = p5*(UUU_slave[nx_block+i-1]* DYU_slave[nx_block+i-1]*  
                                                    DZU_slave[nx_block+i-1] + 
                           UUU_slave[i-1]*          DYU_slave[i-1]*  
                                                    DZU_slave[i-1]);
        VTN_slave[i] = p5*(VVV_slave[nx_block+i]  * DXU_slave[nx_block+i]*  
                                                    DZU_slave[nx_block+i] + 
                           VVV_slave[nx_block+i-1]* DXU_slave[nx_block+i-1]*  
                                                    DZU_slave[nx_block+i-1]);
        VTS_slave[i] = p5*(VVV_slave[i]           * DXU_slave[i]*  
                                                    DZU_slave[i] + 
                           VVV_slave[i-1]*          DXU_slave[i-1]*  
                                                    DZU_slave[i-1]);
        }
		   
      for(i=iend;i<nx_block;i++)
        {		   
           UTE_slave[i]=c0;
           UTW_slave[i]=c0;
           VTN_slave[i]=c0;
           VTS_slave[i]=c0;   
        }            
  
     }  // if(j&)
 
   else
     {
       for(i=0;i<nx_block;i++)
         { 
         UTE_slave[i]=c0;
         UTW_slave[i]=c0;
         VTN_slave[i]=c0;
         VTS_slave[i]=c0;             
         }
     }

	
     if(k<km)
       {		
       get_reply=0;
       athread_get(0,KMT    +(j-1)*nx_block, KMT_slave,    s_int*nx_block,   (void*)&get_reply,0,0,0);
       athread_get(0,TAREA_R+(j-1)*nx_block, TAREA_R_slave,s_double*nx_block,(void*)&get_reply,0,0,0);
       athread_get(0,WTK    +(j-1)*nx_block, WTK_slave,    s_double*nx_block,(void*)&get_reply,0,0,0);
       while(get_reply!=3);
       asm volatile("memb");

       for(i=0;i<nx_block;i++)
         {
         FC_slave[i] = (VTN_slave[i] - VTS_slave[i] + UTE_slave[i] - UTW_slave[i])*TAREA_R_slave[i]; 

     // WTKB = merge(WTK+FC, c0, k < KMT(:,:,bid))
         if(k<KMT_slave[i])
            WTKB_slave[i]=WTK_slave[i]+FC_slave[i];
         else
            WTKB_slave[i]=c0;
         }
       } // if(k)

     else
       {
       for(i=0;i<nx_block;i++)  
          WTKB_slave[i] = c0;
       }   
   

     put_reply=0;
     athread_put(0,WTKB_slave, WTKB+(j-1)*nx_block,s_double*nx_block,  (void*)&put_reply,0,0);
     athread_put(0,VTN_slave,  VTN +(j-1)*nx_block,s_double*4*nx_block,(void*)&put_reply,s_double*(nxy-nx_block),s_double*nx_block);
     while(put_reply!=2); 
     asm volatile("memb");

     } // end of for(j)
  
     return;
}




void loop16_hdifft(struct adhdt_Loop11_Para *Loop16_Para)
{
   int phyid,stid,myid,jbegin,jend,ibegin,iend,bid,k,i,j;
   int array_x,array_y,array_z,thread_num,size_j,size_i,nxy;
   int volatile get_reply,put_reply; 
   double *DTN,*DTS,*DTE,*DTW,*DZT,*CN,*CS,*CE,*CW,*CC;
   int  *KMT, hdifk[TYI];
   double hdift[TYR];
   double *DTN_slave,*DTS_slave,*DTE_slave,*DTW_slave,*DZT_slave,
          *CN_slave,*CS_slave,*CE_slave,*CW_slave,*CC_slave;
   int *KMT_slave;
   double c0=0.0, c1=1.0; 

   struct adhdt_Loop11_Para  Loop16_Para_S;   

   get_reply=0; 
   athread_get(PE_MODE,Loop16_Para,&Loop16_Para_S,sizeof(Loop16_Para_S), (void*)&get_reply,0,0,0);
   asm volatile("memb");   
   while(get_reply!=1);  

    phyid = athread_get_id(-1); 
    stid  = Loop16_Para_S.param[20];	
    myid = phyid - stid;

    thread_num=Loop16_Para_S.param[19];
    if(myid>=thread_num)  return;

    jbegin=Loop16_Para_S.param[0];       
    jend  =Loop16_Para_S.param[1];
    ibegin=Loop16_Para_S.param[2];
    iend  =Loop16_Para_S.param[3];

    array_x=Loop16_Para_S.param[4];
    array_y=Loop16_Para_S.param[5];
    array_z=Loop16_Para_S.param[6];   
	
    nxy = array_x*array_y;

    bid=Loop16_Para_S.param[7];
    k=Loop16_Para_S.param[8];   

    KMT=Loop16_Para_S.duse_arr[0];      
	
    DTN=Loop16_Para_S.use_arr[6];
    DTS=Loop16_Para_S.use_arr[7];
    DTE=Loop16_Para_S.use_arr[8];
    DTW=Loop16_Para_S.use_arr[9];

    DZT=Loop16_Para_S.use_arr[10];

    CN=Loop16_Para_S.use_arr[11];
    CS=Loop16_Para_S.use_arr[11]+  nxy;
    CE=Loop16_Para_S.use_arr[11]+2*nxy;
    CW=Loop16_Para_S.use_arr[11]+3*nxy;
    CC=Loop16_Para_S.use_arr[11]+4*nxy;
          
 
    CN_slave= hdift;   // ldm_malloc(s_double*5*size_i);
    CS_slave = CN_slave +   array_x;
    CE_slave = CN_slave + 2*array_x; 
    CW_slave = CN_slave + 3*array_x;
    CC_slave = CN_slave + 4*array_x;

    DTN_slave= hdift + 6*array_x;   //ldm_malloc(s_double*size_i);
    DTS_slave= hdift + 7*array_x;   //ldm_malloc(s_double*size_i);
    DTE_slave= hdift + 8*array_x;   //ldm_malloc(s_double*size_i);
    DTW_slave= hdift + 9*array_x;   //ldm_malloc(s_double*size_i);
    DZT_slave= hdift +10*array_x;   //ldm_malloc(s_double*(2*array_x+size_i+2));
  
    KMT_slave= hdifk;    //(int*)ldm_malloc(s_int*(2*array_x+size_i+2));        

    for(j=myid+1;j<=array_y;j+=thread_num)
      {
      if(j>=jbegin && j<=jend)
        {
        get_reply = 0;         
        athread_get(0,DZT+(j-1-1)*array_x,DZT_slave,s_double*3*array_x,(void*)&get_reply,0,0,0);  
       
        athread_get(0,DTN+(j-1)*array_x,DTN_slave,s_double*array_x,(void*)&get_reply,0,0,0);
        athread_get(0,DTS+(j-1)*array_x,DTS_slave,s_double*array_x,(void*)&get_reply,0,0,0);
        athread_get(0,DTE+(j-1)*array_x,DTE_slave,s_double*array_x,(void*)&get_reply,0,0,0);
        athread_get(0,DTW+(j-1)*array_x,DTW_slave,s_double*array_x,(void*)&get_reply,0,0,0);

        athread_get(0,KMT+(j-1-1)*array_x,KMT_slave,s_int*3*array_x,(void*)&get_reply,0,0,0);
        asm volatile("memb");
        while(get_reply!=6);

        for(i=0;i<ibegin-1;i++)
          {
          CN_slave[i] = c0;
          CS_slave[i] = c0;
          CE_slave[i] = c0;
          CW_slave[i] = c0;
          CC_slave[i] = c0;
          }

        for(i=ibegin-1;i<iend;i++)
          {
          CN_slave[i]=DTN_slave[i] * locmin( DZT_slave[array_x+i], DZT_slave[2*array_x+i]) * (c1/ DZT_slave[array_x+i]);
          CS_slave[i]=DTS_slave[i] * locmin( DZT_slave[array_x+i], DZT_slave[i]          ) * (c1/ DZT_slave[array_x+i]);
          CE_slave[i]=DTE_slave[i] * locmin( DZT_slave[array_x+i], DZT_slave[array_x+1+i]) * (c1/ DZT_slave[array_x+i]);
          CW_slave[i]=DTW_slave[i] * locmin( DZT_slave[array_x+i], DZT_slave[array_x-1+i]) * (c1/ DZT_slave[array_x+i]);

          if ( k >  KMT_slave[array_x+i] )  //KMT_slave[nx_block+i+1]
               {
               CN_slave[i] = c0;
               CS_slave[i] = c0;
               CE_slave[i] = c0;
               CW_slave[i] = c0;  
               }
          else
               {
               if ( k > KMT_slave[2*array_x+i])
                   { // KMT_slave[2*nx_block+i+1]
                   CN_slave[i] = c0;
                   }
               if ( k > KMT_slave[i])
                   {     // KMT_slave[i+1]
                   CS_slave[i] = c0;
                   }
               if ( k > KMT_slave[array_x+i+1])
                   {   //KMT_slave[nx_block+i+1+1]
                   CE_slave[i] = c0;
                   } 
               if ( k > KMT_slave[array_x+i-1])
                   {   // KMT_slave[nx_block+i]
                   CW_slave[i] = c0;
                   }
               }  // if(k) 
			  
          CC_slave[i] = -(CN_slave[i]+ CS_slave[i]+ CE_slave[i]+ CW_slave[i]);     
          } // for(i)  
		
        for(i=iend;i<array_x;i++)
          {
          CN_slave[i] = c0;
          CS_slave[i] = c0;
          CE_slave[i] = c0;
          CW_slave[i] = c0;
          CC_slave[i] = c0;
          }

        }  // if(j&)

  
     else
        {
        for(i=0;i<array_x;i++)
          {
          CN_slave[i] = c0;
          CS_slave[i] = c0;
          CE_slave[i] = c0;
          CW_slave[i] = c0;
          CC_slave[i] = c0;
          }
        }
      

     put_reply=0;
     athread_put(0,CN_slave,CN+(j-1)*array_x,s_double*5*array_x,(void*)&put_reply,s_double*(nxy-array_x),s_double*array_x);
     asm volatile("memb");
     while(put_reply!=1);
   } // end of for(j)

 
   return;
}







   void slave_adhdt_loop13_fun(struct adhdt_Loop13_Para* Loop13_Para)   
{
    int phyid;
	
    phyid = athread_get_id(-1);

    if(phyid<24)
		loop13_advt(Loop13_Para);

    else if(phyid<48)
		loop17_hdifft(Loop13_Para);
	
    return;
}
  



void loop13_advt(struct adhdt_Loop13_Para *Loop13_Para)
{   
    int volatile get_reply,put_reply,get_reply1;
    int phyid,stid,myid,jb,je,ib,ie,nx_block,ny_block,k,km,bid,nt,ns,thread_num,size_i,size_j,i,j,n;
    int trmask[6],nxy,sfc_layer_type,sfc_layer_varthick=1;
    double *VTN,*UTE,*TRCR,*TAREA_R,*DZT,*LTK,*WTK,*WTKB;
    double *VTN_slave,*UTE_slave,*TRCR_slave,*TAREA_R_slave,*DZT_slave,*LTK_slave,*WTK_slave,*WTKB_slave;
    double advt_13[TYR],dzrk,p5=0.5,c1=1.0,c0=0.0;
    struct adhdt_Loop13_Para Loop13_Para_S;
	
    get_reply=0;
    athread_get(PE_MODE,Loop13_Para,&Loop13_Para_S,sizeof(Loop13_Para_S),
                                              (void*)&get_reply,0,0,0);
    asm volatile("memb");
    while(get_reply!=1);
    
    phyid = athread_get_id(-1);
    stid  = Loop13_Para_S.param[18];
    myid = phyid - stid;

    thread_num=Loop13_Para_S.param[17];	
    if(myid>=thread_num) return;
    
    jb = Loop13_Para_S.param[13];
    je = Loop13_Para_S.param[14];
    ib = Loop13_Para_S.param[15];
    ie = Loop13_Para_S.param[16];
    nx_block = Loop13_Para_S.param[4];
    ny_block = Loop13_Para_S.param[5];
    bid= Loop13_Para_S.param[7];
    k  = Loop13_Para_S.param[8];
    km = Loop13_Para_S.param[9];

    nt = Loop13_Para_S.param[11];
    ns = Loop13_Para_S.param[12]-1;
	
//    tr_mask = Loop13_Para_S.param[12];
//    if(tr_mask>0) return;
    for(n=ns;n<nt;n++)
       trmask[n]=Loop13_Para_S.param[25+n];

    sfc_layer_type = Loop13_Para_S.param[21];


    dzrk = Loop13_Para_S.dparam[1];

	
    nxy = nx_block*ny_block;  

    VTN = Loop13_Para_S.use_arr[0];  //UVNSEW
    UTE = Loop13_Para_S.use_arr[0]+2*nxy;
    TRCR= Loop13_Para_S.use_arr[1];
    TAREA_R=Loop13_Para_S.use_arr[2];
    DZT = Loop13_Para_S.use_arr[3];
    LTK = Loop13_Para_S.use_arr[4];
    WTK = Loop13_Para_S.use_arr[10];
    WTKB= Loop13_Para_S.use_arr[11];	

    VTN_slave = advt_13;
    UTE_slave = advt_13+2*nx_block;

    TAREA_R_slave =advt_13+6*nx_block;
    DZT_slave = advt_13+7*nx_block;
    WTK_slave = advt_13+9*nx_block; 
    WTKB_slave= advt_13+10*nx_block; 	
	
    TRCR_slave= advt_13+12*nx_block;	
    LTK_slave = advt_13+24*nx_block; 
  
 
// for(j=1;j<=ny_block;j++)
//  {   
//  if(j%thread_num == myid) 
 for(j=myid+1;j<=ny_block;j+=thread_num)
   {

   if(j>=jb && j<=je)
      {  
      get_reply=0;
      athread_get(0,VTN+(j-2)*nx_block,VTN_slave,s_double*2*nx_block,(void *)&get_reply,0,0,0);
      athread_get(0,UTE+(j-1)*nx_block,UTE_slave,s_double*nx_block,(void *)&get_reply,0,0,0);
      athread_get(0,TRCR+(j-2)*nx_block+(k-1)*nxy+ns*km*nxy, TRCR_slave,s_double*3*(nt-ns)*nx_block,(void *)&get_reply,0,s_double*(nxy*km-3*nx_block), s_double*3*nx_block);
      athread_get(0,TAREA_R+(j-1)*nx_block,TAREA_R_slave,s_double*nx_block,(void *)&get_reply,0,0,0);
      athread_get(0,DZT+(j-1)*nx_block+k*nxy+(bid-1)*(km+2)*nxy, DZT_slave,s_double*nx_block,(void *)&get_reply,0,0,0);
      asm volatile("memb");
      while(get_reply!=5); 
	  
  
     for(n=0;n<(nt-ns);n++)
      {	 
      for(i=0;i<(ib-1);i++)
     	 LTK_slave[i+n*nx_block] = c0;	  
		 
      for(i=ib-1;i<ie;i++)
         {
         LTK_slave[i+n*nx_block] = p5*((VTN_slave[i+nx_block]-VTN_slave[i]
                            +UTE_slave[i]-UTE_slave[i-1])  
                                               *TRCR_slave[i  +(3*n+1)*nx_block] +           
                          VTN_slave[i+nx_block]*TRCR_slave[i  +(3*n+2)*nx_block] -           
                          VTN_slave[i]*         TRCR_slave[i  +(3*n  )*nx_block] +           
                          UTE_slave[i]*         TRCR_slave[i+1+(3*n+1)*nx_block] -           
                          UTE_slave[i-1]*       TRCR_slave[i-1+(3*n+1)*nx_block])*           
                          TAREA_R_slave[i] * (c1/ DZT_slave[i]);
         }
		
      for(i=ie;i<nx_block;i++)
     	 LTK_slave[i+n*nx_block] = c0;	
      } // for(n)
	 
      }  // if(j&)

   else  // else (j&)
      {
      for(n=0;n<(nt-ns);n++)
        {		  
        for(i=0;i<nx_block;i++)
          LTK_slave[i+n*nx_block] = c0;
	}
      } 



    get_reply1=0;
    athread_get(0,WTK +(j-1)*nx_block,WTK_slave, s_double*nx_block,(void *)&get_reply1,0,0,0);
    athread_get(0,WTKB+(j-1)*nx_block,WTKB_slave,s_double*nx_block,(void *)&get_reply1,0,0,0);
    asm volatile("memb");
    while(get_reply1!=2);	


    for(n=ns;n<nt;n++)
     {
     get_reply1=0;		 
     if(k>1 && k<km)
       {
       athread_get(0,TRCR+(j-1)*nx_block+(k-2)*nxy+n*km*nxy,TRCR_slave,s_double*3*nx_block,(void *)&get_reply1,0,s_double*(ny_block-1)*nx_block,s_double*nx_block);
       asm volatile("memb");
       while(get_reply1!=1);

       for(i=0;i<nx_block;i++)
          LTK_slave[i+(n-ns)*nx_block] = LTK_slave[i+(n-ns)*nx_block]
            + p5*(c1/DZT_slave[i])*WTK_slave[i]*(TRCR_slave[i] + TRCR_slave[nx_block+i])
            - p5*(c1/DZT_slave[i])*WTKB_slave[i]*(TRCR_slave[nx_block+i] + TRCR_slave[2*nx_block+i]);
       }

     else if(k==km)
       {
       athread_get(0,TRCR+(j-1)*nx_block+(k-2)*nxy+n*km*nxy,TRCR_slave,s_double*2*nx_block,(void *)&get_reply1,0,s_double*(ny_block-1)*nx_block,s_double*nx_block);
       asm volatile("memb");
       while(get_reply1!=1);

       for(i=0;i<nx_block;i++)
           LTK_slave[i+(n-ns)*nx_block] = LTK_slave[i+(n-ns)*nx_block]
         	  + p5*(c1/DZT_slave[i])*WTK_slave[i]*(TRCR_slave[i] + TRCR_slave[nx_block+i]);
       }

     else if(k==1)
       {
       athread_get(0,TRCR+(j-1)*nx_block+(k-1)*nxy+n*km*nxy,TRCR_slave,s_double*2*nx_block,(void *)&get_reply1,0,s_double*(ny_block-1)*nx_block,s_double*nx_block);
       asm volatile("memb");
       while(get_reply1!=1);

       if(sfc_layer_type == sfc_layer_varthick)
          {
          for(i=0;i<nx_block;i++)
             LTK_slave[i+(n-ns)*nx_block] = LTK_slave[i+(n-ns)*nx_block] - p5*(c1/DZT_slave[i])*WTKB_slave[i]*(TRCR_slave[i] + TRCR_slave[nx_block+i]);
          }
       else
          {
          for(i=0;i<nx_block;i++)
             LTK_slave[i+(n-ns)*nx_block] = LTK_slave[i+(n-ns)*nx_block] + dzrk*WTK_slave[i]*TRCR_slave[i]
                - p5*(c1/DZT_slave[i])*WTKB_slave[i]*(TRCR_slave[i] + TRCR_slave[nx_block+i]);
          }
       }

    if(trmask[n]>0) 
      { 	 
      put_reply=0;
      athread_put(0,LTK_slave+(n-ns)*nx_block,LTK+(j-1)*nx_block+n*nxy, s_double*nx_block,(void *)&put_reply,0,0);
      asm volatile("memb");
      while(put_reply!=1);
      }

    } // for(n)

//    } // if(j%)	 

  } // end of for(j)

  return;
}





void loop17_hdifft(struct adhdt_Loop13_Para *Loop17_Para)
{ 
   int phyid,stid,myid;
   int jbegin,jend,ibegin,iend,ib,ie,jb,je,nt,ns;
   int array_x,array_y,array_z,nxy,i,j,n,k,bid,thread_num;
   int volatile get_reply,put_reply; 
   double hdift[TYR], ah, c0=0.0;
   double *AHF,*D2TK,*HDTK,*TMIX,*CC,*CN,*CS,*CE,*CW;
   double *AHF_slave,*HDTK_slave,*D2TK_slave,*TMIX_slave,
          *CC_slave, *CN_slave,*CS_slave,*CE_slave,*CW_slave;

   struct adhdt_Loop13_Para  Loop17a_Para_S;


    get_reply=0;
    athread_get(PE_MODE,Loop17_Para,&Loop17a_Para_S,sizeof(Loop17a_Para_S), (void*)&get_reply,0,0,0);
    asm volatile("memb");
    while(get_reply!=1);
 
    phyid = athread_get_id(-1);
    stid   = Loop17a_Para_S.param[20];
    myid = phyid - stid;
	
    thread_num=Loop17a_Para_S.param[19];
    if(myid>=thread_num)  return;


    jbegin=Loop17a_Para_S.param[0];  
    jend=Loop17a_Para_S.param[1];   
    ibegin=Loop17a_Para_S.param[2];   
    iend=Loop17a_Para_S.param[3];

    jb=Loop17a_Para_S.param[13];  
    je=Loop17a_Para_S.param[14];   
    ib=Loop17a_Para_S.param[15];   
    ie=Loop17a_Para_S.param[16];


    array_x=Loop17a_Para_S.param[4];
    array_y=Loop17a_Para_S.param[5];
    array_z=Loop17a_Para_S.param[6];
	
    nxy = array_x*array_y; 
	
    k=Loop17a_Para_S.param[8];
    bid=Loop17a_Para_S.param[7]; 

//    n=Loop17a_Para_S.param[11];
    nt  = Loop17a_Para_S.param[11];
    ns  = Loop17a_Para_S.param[12]-1;
	

    ah=Loop17a_Para_S.dparam[0]; 
	
    AHF=Loop17a_Para_S.use_arr[5];       
    HDTK=Loop17a_Para_S.use_arr[6];
    D2TK=Loop17a_Para_S.use_arr[9];
   
    TMIX=Loop17a_Para_S.use_arr[7];
           
    CN=Loop17a_Para_S.use_arr[8];
    CS=Loop17a_Para_S.use_arr[8]+  nxy;
    CE=Loop17a_Para_S.use_arr[8]+2*nxy;
    CW=Loop17a_Para_S.use_arr[8]+3*nxy;
    CC=Loop17a_Para_S.use_arr[8]+4*nxy;


    CN_slave= hdift;    
    CS_slave= CN_slave + 3*array_x;
    CE_slave= CN_slave + 6*array_x;
    CW_slave= CN_slave + 9*array_x;
    CC_slave= CN_slave + 12*array_x;

    AHF_slave = hdift + 16*array_x;    
    TMIX_slave= hdift + 19*array_x;

    D2TK_slave= TMIX_slave + 5*(nt-ns)*array_x;
    HDTK_slave= D2TK_slave + 3*(nt-ns)*array_x;  

    for(j=1;j<=array_y;j++)
       {
       if(j%thread_num == myid)     
         {

         if(j>=jb && j<=je)
            {
            get_reply = 0; 
            athread_get(0,AHF+(j-2)*array_x+(bid-1)*nxy,AHF_slave,s_double*3*array_x,(void*)&get_reply,0,0,0);          
            athread_get(0,TMIX+(j-3)*array_x+(k-1)*nxy+ns*nxy*array_z,TMIX_slave,s_double*5*(nt-ns)*array_x,(void*)&get_reply,0,s_double*(nxy*array_z-5*array_x),s_double*5*array_x);
            athread_get(0,CN+(j-2)*array_x, CN_slave,s_double*15*array_x,(void*)&get_reply,0,s_double*(nxy-3*array_x),s_double*3*array_x); 
            asm volatile("memb");
            while(get_reply!=3);


          for(n=0;n<(nt-ns);n++)
            {
//            for(i=0;i<array_x;i++) 
//              {
//              if(i>=(ibegin-1) && i<=(iend-1))
            for(i=0;i<(ibegin-1);i++)
                {
                D2TK_slave[i+          3*n*array_x] = c0;
                D2TK_slave[i+  array_x+3*n*array_x] = c0;
                D2TK_slave[i+2*array_x+3*n*array_x] = c0;
                }
            for(i=ibegin-1;i<iend;i++)
                {
                D2TK_slave[i          +3*n*array_x] = AHF_slave[i] *
                               (CC_slave[i] * TMIX_slave[i  +  array_x+5*n*array_x] +
                                CN_slave[i] * TMIX_slave[i  +2*array_x+5*n*array_x] + 
                                CS_slave[i] * TMIX_slave[i            +5*n*array_x] +
                                CE_slave[i] * TMIX_slave[i+1+  array_x+5*n*array_x] + 
                                CW_slave[i] * TMIX_slave[i-1+  array_x+5*n*array_x]); 

                D2TK_slave[i+  array_x+3*n*array_x] = AHF_slave[i+array_x] *
                               (CC_slave[i+array_x] * TMIX_slave[i  +2*array_x+5*n*array_x] +
                                CN_slave[i+array_x] * TMIX_slave[i  +3*array_x+5*n*array_x] + 
                                CS_slave[i+array_x] * TMIX_slave[i  +  array_x+5*n*array_x] +
                                CE_slave[i+array_x] * TMIX_slave[i+1+2*array_x+5*n*array_x] + 
                                CW_slave[i+array_x] * TMIX_slave[i-1+2*array_x+5*n*array_x]); 

                D2TK_slave[i+2*array_x+3*n*array_x] = AHF_slave[i+2*array_x] *
                               (CC_slave[i+2*array_x] * TMIX_slave[i  +3*array_x+5*n*array_x] +
                                CN_slave[i+2*array_x] * TMIX_slave[i  +4*array_x+5*n*array_x] + 
                                CS_slave[i+2*array_x] * TMIX_slave[i  +2*array_x+5*n*array_x] +
                                CE_slave[i+2*array_x] * TMIX_slave[i+1+3*array_x+5*n*array_x] + 
                                CW_slave[i+2*array_x] * TMIX_slave[i-1+3*array_x+5*n*array_x]); 
                }
            for(i=iend;i<array_x;i++)
                {
                D2TK_slave[i          +3*n*array_x] = c0;
                D2TK_slave[i+  array_x+3*n*array_x] = c0;
                D2TK_slave[i+2*array_x+3*n*array_x] = c0;
                }
//              } // end of for(i_D2TK)


//            for(i=0;i<array_x;i++) 
//              {
//              if(i>=(ib-1) && i<=(ie-1))
            for(i=0;i<(ib-1);i++)
                HDTK_slave[i+n*array_x] = c0;
            for(i=ib-1;i<ie;i++)
                HDTK_slave[i+n*array_x] = ah*(CC_slave[i+array_x] * D2TK_slave[i  +  array_x+3*n*array_x] +
                                    CN_slave[i+array_x] * D2TK_slave[i  +2*array_x+3*n*array_x] +
                                    CS_slave[i+array_x] * D2TK_slave[i            +3*n*array_x] +
                                    CE_slave[i+array_x] * D2TK_slave[i+1+  array_x+3*n*array_x] +
                                    CW_slave[i+array_x] * D2TK_slave[i-1+  array_x+3*n*array_x]);   
            for(i=ie;i<array_x;i++)
                HDTK_slave[i+n*array_x] = c0;
//              } // end of for(i_HDTK)
            } // for(n)	
				
           } // if(j&)


         else  //else if(j&)
           {
           for(n=0;n<(nt-ns);n++) 
            {				   
            for(i=0;i<array_x;i++) 
              HDTK_slave[i+n*array_x] = c0;
	    }
           }


         put_reply=0;
         athread_put(0,HDTK_slave,        HDTK+(j-1)*array_x+ns*nxy, s_double*(nt-ns)*array_x,(void*)&put_reply,s_double*(nxy-array_x),s_double*array_x);   
         athread_put(0,D2TK_slave+array_x,D2TK+(j-1)*array_x+ns*nxy, s_double*(nt-ns)*array_x,(void*)&put_reply,s_double*(nxy-array_x),s_double*array_x);   		 
         asm volatile("memb");
         while(put_reply!=2);
                
         } // end of if(j%)
     } // end of for(j)

   return;
}
