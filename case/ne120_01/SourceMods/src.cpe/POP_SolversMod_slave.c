#include<stdio.h>
#include<slave.h>
#include"POP_SolversMod_struct.h"
#include"cpe_print.h"
#include <unistd.h>


#define get_myrid(row)  asm volatile("rcsr   %0, 1" : "=r"(row))
#define get_mycid(col)  asm volatile("rcsr   %0, 2" : "=r"(col))

#define   CPE_N 64

//POP_Solversmod##############################################################

void fun_pop_solvers_operator(struct param_zy *s1) 
{ 
      int myid;
      myid           = athread_get_id(-1); 
        
      if(myid>=24){return ;}

      double *btropWgtCenter,*btropWgtNorth,*btropWgtEast,*btropWgtNE,*AX,*X,*p_btropWgtCenter,*p_btropWgtNorth,*p_btropWgtEast,
               *p_btropWgtNE,*p_AX,*p_X;

      volatile int  get_reply,put_reply;
      int offset1,offset2,i,j,k,istart,jstart,taskall,tasknum,tmp_rest,array_index,tindex;
      double s_zero=0.0;
      int nx,ny;

      double ldm_arr[4000];

      struct param_zy sldm;

      get_reply=0;
      athread_get(0,s1,&sldm,sizeof(sldm),(void*)&get_reply,0,0,0);
      while(get_reply!=1);asm volatile("memb");
      
      nx             = sldm.param[0];
      ny             = sldm.param[1];
      btropWgtCenter = sldm.use_arr[0];
      btropWgtNorth  = sldm.use_arr[1];
      btropWgtEast   = sldm.use_arr[2];
      btropWgtNE     = sldm.use_arr[3];
      AX             = sldm.use_arr[4];
      X              = sldm.use_arr[5]; 
 
        p_btropWgtCenter = ldm_arr ;                // ldm_malloc(8*size_i);
        p_X              = p_btropWgtCenter+ 256;// ldm_malloc(8*(2*nx+size_i+2));
        p_btropWgtNorth  = p_X             + 256 ; //ldm_malloc(8*(nx+size_i));
        p_btropWgtEast   = p_btropWgtNorth + 256;		//ldm_malloc(8*(size_i+1));
        p_btropWgtNE     = p_btropWgtEast  + 256;		//ldm_malloc(8*(nx+size_i+1));
        p_AX             = p_btropWgtNE    + 256;		//ldm_malloc(8*nx);

        istart=2;
        jstart=2;

        taskall=ny-2;
        tasknum= taskall/24;
        tmp_rest=taskall-tasknum*24;

     
        if(myid<tmp_rest)
        {
                tasknum ++;
                array_index=jstart+myid*tasknum;
        }
        else
        {
                array_index=jstart+tmp_rest*(tasknum+1)+(myid-tmp_rest)*tasknum;
        }

          tindex=0;
          for(;tindex<tasknum;tindex++)
          {
          get_reply=0;
          offset1= (tindex+array_index-1)*nx;    
          offset2= (tindex+array_index-1-1)*nx;;//(i,j-1)
          athread_get(0,X             +offset2     ,p_X               ,8*(3*nx),(void*)&get_reply,0,0,0); 
          
          offset1= (tindex+array_index-1)*nx;
          athread_get(0,btropWgtCenter+offset1 +1  ,p_btropWgtCenter  ,8*(nx-2)         ,(void*)&get_reply,0,0,0);
          athread_get(0,btropWgtNorth +offset2 +1  ,p_btropWgtNorth   ,8*(2*nx-2)    ,(void*)&get_reply,0,0,0); //j-1,j 
          athread_get(0,btropWgtEast  +offset1     ,p_btropWgtEast    ,8*(nx-1)     ,(void*)&get_reply,0,0,0);
          athread_get(0,btropWgtNE    +offset2     ,p_btropWgtNE      ,8*(2*nx-1)  ,(void*)&get_reply,0,0,0);
         
          put_reply=0;
          p_AX[0]   =s_zero;
          p_AX[nx-1]=s_zero;
          while(get_reply!=5);        
          asm volatile("memb"); 
          for(i=0;i<nx-2;i++)
          {
          p_AX[i+1]  = p_btropWgtCenter[i]     *    p_X[nx+1 +i]     +
                       p_btropWgtNorth[nx +i]  *    p_X[2*nx+1 +i]   +
                       p_btropWgtNorth[i]      *    p_X[1 +i]        +
                       p_btropWgtEast[1 +i]    *    p_X[nx+2 +i]     +
                       p_btropWgtEast[i]       *    p_X[nx +i]       +
                       p_btropWgtNE[nx+1 +i]   *    p_X[2*nx+2 +i]   +
                       p_btropWgtNE[1 +i]      *    p_X[2 +i]        +
                       p_btropWgtNE[nx +i]     *    p_X[2*nx +i]     +
                       p_btropWgtNE[i]         *    p_X[i] ;     
           
          }    
          athread_put(0,p_AX   ,AX+offset1    ,8*nx,(void *)&put_reply,0,0);            
          while(put_reply!=1);
          asm volatile("memb");
         }
}
void fun_pop_solvers_step5(struct param_zy *s1)
{
 
     int myid,my_rid,my_cid,last_slave,pre,pre2; 
           
     volatile int  get_reply,put_reply,get_reply1,get_reply2;
     int offset1,offset2,my_task,i,j,nx,ny,k,istart,jstart,taskall,tasknum,tmp_rest,array_index,tindex;
     double c1=1.0,c0=0.0;
     struct param_zy sldm;
     double ldm_arr[4096],csomga,csy,*Q,*R,*X,*p_Q,*p_R,*p_X,*S,*p_S,*B,*p_B;
     double *btropWgtCenter,*btropWgtNorth,*btropWgtEast,*btropWgtNE,*p_btropWgtCenter,*p_btropWgtNorth,*p_btropWgtEast,
               *p_btropWgtNE,*p_last,*left_x,*right_x,*right_nor,*right_ne;

     get_reply=0;
     athread_get(0,s1 ,&sldm,sizeof(param_zy),(void*)&get_reply,0,0,0);

     get_myrid(my_rid);
     get_mycid(my_cid);

     if(my_rid%2)
         myid = my_rid*8 + (7-my_cid);
     else
         myid = my_rid*8 + my_cid;

     while(get_reply!=1);asm volatile("memb");

    csomga      = sldm.dparam[0];
    csy         = sldm.dparam[1];
    nx          = sldm.param[0];
    ny          = sldm.param[1];


    if(myid >= ny){}
     else{
     
    Q           = sldm.use_arr[0];
    R           = sldm.use_arr[1];
    X           = sldm.use_arr[2];
    S           = sldm.use_arr[3];
    btropWgtCenter = sldm.use_arr[4];
    btropWgtNorth  = sldm.use_arr[5];
    btropWgtEast   = sldm.use_arr[6];
    btropWgtNE     = sldm.use_arr[7];
    B              = sldm.use_arr[8];

    p_R         = ldm_arr;
    p_Q         = p_R+128;
    p_B         = p_Q+128;

    p_btropWgtCenter = p_B+128;
    p_X              = p_btropWgtCenter+ 128;
    p_btropWgtNorth  = p_X             + 128;
    p_btropWgtEast   = p_btropWgtNorth + 256;
    p_btropWgtNE     = p_btropWgtEast  + 128;
    p_S              = p_btropWgtNE    + 256;
    left_x           = p_S             + 128;

    right_x=left_x+128;
    right_nor=right_x+128;
    right_ne =right_nor+128;

    taskall=ny;
    tasknum= taskall/CPE_N;
    tmp_rest=taskall-tasknum*CPE_N;
    if(myid<tmp_rest)
    {
         tasknum ++;//this thread need +1 task
         array_index=1+myid*tasknum ;
    }
    else    
         array_index=1+tmp_rest*(tasknum+1) + (myid-tmp_rest)*tasknum ;
    
    tindex=0;
    for(;tindex<tasknum;tindex++)
    {
                offset1=(tindex+array_index-1)*nx;                
                get_reply=0;
                athread_get(0,R             +offset1          ,p_R        ,8*nx    ,(void *)&get_reply,0,0,0);
		athread_get(0,Q             +offset1          ,p_Q        ,8*nx    ,(void *)&get_reply,0,0,0);
                athread_get(0,X             +offset1          ,p_X        ,8*nx    ,(void *)&get_reply,0,0,0);
                get_reply1=0;	
                offset2=(tindex+array_index-1-1)*nx;
                
                athread_get(0,btropWgtCenter+offset1   ,p_btropWgtCenter  ,8*nx    ,(void *)&get_reply1,0,0,0);
         	athread_get(0,btropWgtNorth +offset2   ,p_btropWgtNorth   ,8*nx*2  ,(void *)&get_reply1,0,0,0);
         	athread_get(0,btropWgtEast  +offset1   ,p_btropWgtEast    ,8*nx    ,(void *)&get_reply1,0,0,0);
         	athread_get(0,btropWgtNE    +offset2   ,p_btropWgtNE      ,8*nx*2  ,(void *)&get_reply1,0,0,0);
            
                put_reply =0;
                get_reply2=0;
                athread_get(0,B             +offset1          ,p_B        ,8*nx    ,(void *)&get_reply2,0,0,0);

                while(get_reply!=3);
//                asm volatile("memb");

                for(i=0;i<nx;i++)
                {
                    p_Q[i]=csomga * p_R[i] + (csy * csomga - c1) * p_Q[i];
                    p_X[i]=p_X[i] + p_Q[i];
                    p_S[i]=c0;   
                }
          
                athread_put(0,p_Q       ,Q  +offset1    ,8*nx ,(void *)&put_reply,0,0);
                athread_put(0,p_X       ,X  +offset1    ,8*nx ,(void *)&put_reply,0,0);
                while(get_reply1!=4); 
               
                //Reg_comm_zy(nx,ny,p_X,p_S,p_btropWgtCenter,p_btropWgtNorth,p_btropWgtEast,p_btropWgtNE,
                //         my_rid,my_cid,myid,left_x,right_x,right_nor,right_ne); 
 
                while(get_reply2!=1);
          
                for(i=0;i<nx;i++)
                {
                    p_R[i]=p_B[i] - p_S[i];
                }
                athread_put(0,p_R      ,R  +offset1    ,8*nx ,(void *)&put_reply,0,0);       
                while(put_reply!=3);
                asm volatile("memb");
    }   
  }//if >=ny
}

void fun_pop_solvers_step5_2(struct param_zy *s1)     //just of  m =60  or m % 10= 0
{
 
     int myid,my_rid,my_cid,last_slave,pre,pre2; 
           
     volatile int  get_reply,put_reply,get_reply1,get_reply2;
     int offset1,offset2,my_task,i,j,nx,ny,k,istart,jstart,taskall,tasknum,tmp_rest,array_index,tindex;
     double c1=1.0,c0=0.0;
     struct param_zy sldm;
     double ldm_arr[4096],csomga,csy,*Q,*R,*X,*p_Q,*p_R,*p_X,*S,*p_S,*B,*p_B,*work0,*p_work0;
     double *btropWgtCenter,*btropWgtNorth,*btropWgtEast,*btropWgtNE,*p_btropWgtCenter,*p_btropWgtNorth,*p_btropWgtEast,
               *p_btropWgtNE,*p_last,*left_x,*right_x,*right_nor,*right_ne;

     get_reply=0;
     athread_get(0,s1 ,&sldm,sizeof(param_zy),(void*)&get_reply,0,0,0);

     get_myrid(my_rid);
     get_mycid(my_cid);

     if(my_rid%2)
         myid = my_rid*8 + (7-my_cid);
     else
         myid = my_rid*8 + my_cid;

     while(get_reply!=1);asm volatile("memb");

    csomga      = sldm.dparam[0];
    csy         = sldm.dparam[1];
    nx          = sldm.param[0];
    ny          = sldm.param[1];


    if(myid >= ny){}
     else{
     
    Q              = sldm.use_arr[0];
    R              = sldm.use_arr[1];
    X              = sldm.use_arr[2];
    S              = sldm.use_arr[3];
    btropWgtCenter = sldm.use_arr[4];
    btropWgtNorth  = sldm.use_arr[5];
    btropWgtEast   = sldm.use_arr[6];
    btropWgtNE     = sldm.use_arr[7];
    B              = sldm.use_arr[8];
    work0          = sldm.use_arr[9];


    p_R         = ldm_arr;
    p_Q         = p_R+128;
    p_B         = p_Q+128;

    p_btropWgtCenter = p_B+128;
    p_X              = p_btropWgtCenter+ 128;
    p_btropWgtNorth  = p_X             + 128;
    p_btropWgtEast   = p_btropWgtNorth + 256;
    p_btropWgtNE     = p_btropWgtEast  + 128;
    p_S              = p_btropWgtNE    + 256;
    left_x           = p_S             + 128;
    p_work0          = left_x          + 128;
    right_x          = p_work0           +128;
    right_nor        = right_x+128;
    right_ne         = right_nor+128;

    taskall=ny;
    tasknum= taskall/CPE_N;
    tmp_rest=taskall-tasknum*CPE_N;
    if(myid<tmp_rest)
    {
         tasknum ++;//this thread need +1 task
         array_index=1+myid*tasknum ;
    }
    else    
         array_index=1+tmp_rest*(tasknum+1) + (myid-tmp_rest)*tasknum ;
    
    tindex=0;
    for(;tindex<tasknum;tindex++)
    {
                offset1=(tindex+array_index-1)*nx;                
                get_reply=0;
                athread_get(0,R             +offset1          ,p_R        ,8*nx    ,(void *)&get_reply,0,0,0);
		athread_get(0,Q             +offset1          ,p_Q        ,8*nx    ,(void *)&get_reply,0,0,0);
                athread_get(0,X             +offset1          ,p_X        ,8*nx    ,(void *)&get_reply,0,0,0);
                get_reply1=0;	
                offset2=(tindex+array_index-1-1)*nx;
                
                athread_get(0,btropWgtCenter+offset1   ,p_btropWgtCenter  ,8*nx    ,(void *)&get_reply1,0,0,0);
         	athread_get(0,btropWgtNorth +offset2   ,p_btropWgtNorth   ,8*nx*2  ,(void *)&get_reply1,0,0,0);
         	athread_get(0,btropWgtEast  +offset1   ,p_btropWgtEast    ,8*nx    ,(void *)&get_reply1,0,0,0);
         	athread_get(0,btropWgtNE    +offset2   ,p_btropWgtNE      ,8*nx*2  ,(void *)&get_reply1,0,0,0);
            
                put_reply =0;
                get_reply2=0;
                athread_get(0,B             +offset1          ,p_B        ,8*nx    ,(void *)&get_reply2,0,0,0);

                while(get_reply!=3);
//                asm volatile("memb");

                for(i=0;i<nx;i++)
                {
                    p_Q[i]=csomga * p_R[i] + (csy * csomga - c1) * p_Q[i];
                    p_X[i]=p_X[i] + p_Q[i];
                    p_S[i]=c0;   
                }
          
                athread_put(0,p_Q       ,Q  +offset1    ,8*nx ,(void *)&put_reply,0,0);
                athread_put(0,p_X       ,X  +offset1    ,8*nx ,(void *)&put_reply,0,0);
                while(get_reply1!=4); 
               
                //Reg_comm_zy(nx,ny,p_X,p_S,p_btropWgtCenter,p_btropWgtNorth,p_btropWgtEast,p_btropWgtNE,
                //         my_rid,my_cid,myid,left_x,right_x,right_nor,right_ne); 
 
                while(get_reply2!=1);
          
                for(i=0;i<nx;i++)
                {
                    p_R[i]=p_B[i] - p_S[i];
                }
           
                athread_put(0,p_R      ,R  +offset1    ,8*nx ,(void *)&put_reply,0,0);       
           
                for(i=0;i<nx;i++)
                {
                    p_work0[i]=p_R[i] * p_R[i];
                }  

               athread_put(0,p_work0   ,work0  +offset1    ,8*nx ,(void *)&put_reply,0,0);
               while(put_reply!=4);
                asm volatile("memb");
    }   
  }//if >=ny
}


void fun_pop_solvers_preconditioner(struct param_zy *s1)
{
     int myid,nx,ny,m,iblock,EvpYnb,EvpXnb,offset,offsetib,offsetrinv;
     volatile int get_reply=0,put_reply=0,get_reply1=0,get_reply2=0,get_reply3=0,get_replyf=0; 
     int* EvpYbidx,*EvpXbidx,*landIndx;
     int volatile  *p_landIndx,*p_EvpYbidx,*p_EvpXbidx;
     double*R,*f  ,*work1,*InvEvpCenterWgt,*EvpCenterWgt,*EvpNeWgt,*EvpRinv,
           *p_work1,*p_R,*p_f,*InvEvpNeWgt,*p_r,*p_y ;
     double volatile *p_InvEvpCenterWgt,*p_EvpCenterWgt,*p_EvpNeWgt,*p_EvpRinv,*p_InvEvpNeWgt;
     double  *p_InvEvpCenterWgt1,*p_EvpCenterWgt1,*p_EvpNeWgt1,*p_EvpRinv1,*p_InvEvpNeWgt1;
     struct param_zy sldm;
     int   ldm_arr[256];
     double ldm_arr2[3200];
     
     int i,j,je,js,lm,ln,l,ib,is,ie,Indxj,Indxi,ybs=8,xbs=8,kk,bsize=0,stripe=0;
     double s_zero=0.0;

     myid = athread_get_id(-1);
   
     get_reply=0;
     athread_get(0,s1,&sldm,sizeof(sldm),(void*)&get_reply,0,0,0);
     while(get_reply!=1);asm volatile("memb");

      nx        = sldm.param[0];
      ny        = sldm.param[1];
      EvpYnb    = sldm.param[2];
      EvpXnb    = sldm.param[3];
      m         = sldm.param[4];

      if(myid >= EvpXnb *EvpYnb){return ;}
     
      EvpYbidx  = sldm.use_arr_int[0];
      EvpXbidx  = sldm.use_arr_int[1];
      landIndx  = sldm.use_arr_int[2];

      R               = sldm.use_arr[0];
      InvEvpCenterWgt = sldm.use_arr[1];
      EvpCenterWgt    = sldm.use_arr[2];
      EvpNeWgt        = sldm.use_arr[3];
      EvpRinv         = sldm.use_arr[4];
      InvEvpNeWgt     = sldm.use_arr[5];
      
       p_EvpYbidx        =   ldm_arr;
       p_EvpXbidx        =   p_EvpYbidx        +64;
       p_landIndx        =   p_EvpXbidx        +64;

       p_R               =  ldm_arr2            ;
       p_f               =   p_R               +128;
       p_work1           =   p_f               +128;
       p_y               =   p_work1           +128;
       p_r               =   p_y               +128;    
 
       p_EvpNeWgt        =    ldm_arr2         +6*128;
       p_InvEvpCenterWgt =   p_EvpNeWgt        +128;
       p_InvEvpNeWgt     =   p_InvEvpCenterWgt +128;
       p_EvpCenterWgt    =   p_InvEvpNeWgt     +128;
       p_EvpRinv         =   p_EvpCenterWgt    +128; 
     
     if(m==1)
     {

       get_reply=0;
       athread_get(0,EvpYbidx         ,p_EvpYbidx          ,4*(EvpYnb+1), (void *)&get_reply,0,0,0);
       athread_get(0,EvpXbidx         ,p_EvpXbidx          ,4*(EvpXnb+1), (void *)&get_reply,0,0,0);
    
       j  = myid / EvpXnb ;
       i  = myid % EvpXnb ;
       ib = j    * EvpXnb+i;
       while(get_reply!=2);

       
       js = p_EvpYbidx[j];
       je = p_EvpYbidx[j+1]+1;
       lm = (je-js) +1;
       is = p_EvpXbidx[i];
       ie = p_EvpXbidx[i+1] +1;
       ln = (ie-is)+1;
       l  = ln + lm -5;
 
       offset    = js*nx + is;
       offsetib  = ib*(ybs+2)*(xbs+2);
       offsetrinv= ib*(ybs+xbs-1)*(ybs+xbs-1);
       
       stripe = 8*(xbs+2-ln);
       bsize  = 8*ln;

       get_reply1=0;
       athread_get(0,R               +offset    ,p_R                 ,8*(ln-2)*(lm-2) , (void*)&get_reply1,0,8*(nx-(ln-2)),8*(ln-2));
       athread_get(0,landIndx                   ,p_landIndx          ,4*EvpXnb*EvpYnb , (void*)&get_reply1,0,0,0);

       get_reply2=0; 
       athread_get(0,InvEvpCenterWgt +offsetib  ,p_InvEvpCenterWgt   ,8*ln*lm         , (void*)&get_reply2,0,stripe,bsize);     

       get_reply3=0;
       athread_get(0,EvpNeWgt        +offsetib  ,p_EvpNeWgt          ,8*ln*lm         , (void*)&get_reply3,0,stripe,bsize);
       athread_get(0,InvEvpNeWgt     +offsetib  ,p_InvEvpNeWgt       ,8*ln*lm         , (void*)&get_reply3,0,stripe,bsize); 
       athread_get(0,EvpCenterWgt    +offsetib  ,p_EvpCenterWgt      ,8*ln*lm         , (void*)&get_reply3,0,stripe,bsize);
       athread_get(0,EvpRinv         +offsetrinv,p_EvpRinv           ,8*l*l           , (void*)&get_reply3,0,8*(ybs+xbs-1-l),8*l);

       put_reply=0; 
       while(get_reply1!=2);     
    //   asm volatile("memb");      
    
        if(p_landIndx[ib]==1)
        {
                    while(get_reply2!=1);
                    for(Indxj=0;Indxj<lm-2;Indxj++)
                    for(Indxi=0;Indxi<ln-2;Indxi++)
                        p_work1[Indxj*(ln-2) + Indxi] =  p_R[Indxj*(ln-2)+Indxi]  * p_InvEvpCenterWgt[(Indxj+1)*ln+Indxi+1 ];


        }       
        else if(p_landIndx[ib]==0)
        {

                    for(Indxj=0;Indxj<lm;Indxj++)
                    for(Indxi=0;Indxi<ln;Indxi++)
                    {
                              p_f[Indxj*ln+Indxi]=s_zero;
                              p_y[Indxj*ln+Indxi]=s_zero;
                    }
                    for(Indxj=0;Indxj<lm-2;Indxj++)
                    for(Indxi=0;Indxi<ln-2;Indxi++)
                              p_f[(Indxj+1)*ln+Indxi+1] =  p_R[Indxj*(ln-2)+Indxi];
                    
                    while(get_reply2!=1);

                    while(get_reply3!=4);

                    for(Indxj=1;Indxj<lm-1;Indxj++)
                    for(Indxi=1;Indxi<ln-1;Indxi++)
                    {
                        p_y[(Indxj+1)*ln+1+Indxi] =
                                  (p_f       [ Indxj   *ln+Indxi  ] 
                             - p_EvpCenterWgt[ Indxj   *ln+Indxi  ] * p_y[ Indxj    *ln+Indxi  ]
                                  -p_EvpNeWgt[(Indxj-1)*ln+Indxi  ] * p_y[(Indxj-1) *ln+Indxi+1]
                                  -p_EvpNeWgt[ Indxj   *ln+Indxi-1] * p_y[(Indxj+1) *ln+Indxi-1]
                                  -p_EvpNeWgt[(Indxj-1)*ln+Indxi-1] * p_y[(Indxj-1) *ln+Indxi-1])*p_InvEvpNeWgt[Indxj*ln+Indxi];

                    }
                    for(Indxi=0;Indxi<ln-2;Indxi++)
                         p_r[Indxi]= p_y[ln*(lm-1)+ Indxi+2];
      

                    for(Indxi=ln-2;Indxi<l;Indxi++)
                          p_r[Indxi] = p_y[ln-1  + (lm -1 + (ln - 3 - Indxi) )*ln];  

                    for(Indxj=0;Indxj<lm-2;Indxj++)
                    for(kk=0;kk<l;kk++)
                          p_y[1+(lm-Indxj-2)*ln]  = p_y[1+(lm-Indxj-2)*ln] + p_EvpRinv[Indxj*l+kk] * p_r[kk]; 

                    for( Indxi = 0;Indxi<ln-3;Indxi++)
                    for( kk = 0;kk<l;kk++)
                          p_y[Indxi+2 + ln*1]     = p_y[Indxi+2+ ln*1]     + p_EvpRinv[kk+(lm-2+Indxi)*l] * p_r[kk];  

                    for(Indxj=1;Indxj<lm-2;Indxj++)
                    for(Indxi=1;Indxi<ln-2;Indxi++)
                      {   p_y[Indxi+1+ln*(Indxj+1)]  = (p_f[Indxi+ln*Indxj]- p_EvpCenterWgt[Indxi+ln*Indxj] * p_y[Indxi+ln  *Indxj    ]
                                - p_EvpNeWgt [Indxi+ln  *(Indxj-1)]   * p_y[Indxi+1+ln*(Indxj-1)]
                                - p_EvpNeWgt [Indxi-1+ln* Indxj   ]   * p_y[Indxi-1+ln*(Indxj+1)]
                                - p_EvpNeWgt [Indxi-1+ln*(Indxj-1)]   * p_y[Indxi-1+ln*(Indxj-1)] ) *p_InvEvpNeWgt[Indxi+ln*Indxj];
                      }
                    for(Indxj=0;Indxj<lm-2;Indxj++)
                    for(Indxi=0;Indxi<ln-2;Indxi++)
                         p_work1[Indxi+(ln-2)*Indxj] = p_y[(Indxj+1)*ln+1+Indxi ];
        }//ib
      }//m
     else if(m !=1)
      {
     //  asm volatile("memb");
       j  = myid / EvpXnb ;
       i  = myid % EvpXnb ;
       ib = j    * EvpXnb+i;
       
       js = p_EvpYbidx[j];
       je = p_EvpYbidx[j+1]+1;
       lm = (je-js) +1;
       is = p_EvpXbidx[i];
       ie = p_EvpXbidx[i+1] +1;
       ln = (ie-is)+1;
       l  = ln + lm -5;
 
       offset    = js*nx + is;

       get_reply1=0;
       athread_get(0,R               +offset    ,p_R                 ,8*(ln-2)*(lm-2) , (void*)&get_reply1,0,8*(nx-(ln-2)),8*(ln-2));
       while(get_reply1!=1);     
    //   asm volatile("memb"); 
        if(p_landIndx[ib]==1)
        {
                    for(Indxj=0;Indxj<lm-2;Indxj++)
                    for(Indxi=0;Indxi<ln-2;Indxi++)
                        p_work1[Indxj*(ln-2) + Indxi] =  p_R[Indxj*(ln-2)+Indxi]  * p_InvEvpCenterWgt[(Indxj+1)*ln+Indxi+1 ];
        }       
        else if(p_landIndx[ib]==0)
        {
                    for(Indxj=0;Indxj<lm;Indxj++)
                    for(Indxi=0;Indxi<ln;Indxi++)
                    {
                              p_f[Indxj*ln+Indxi]=s_zero;
                              p_y[Indxj*ln+Indxi]=s_zero;
                    }
                    for(Indxj=0;Indxj<lm-2;Indxj++)
                    for(Indxi=0;Indxi<ln-2;Indxi++)
                              p_f[(Indxj+1)*ln+Indxi+1] =  p_R[Indxj*(ln-2)+Indxi];
                    

                    for(Indxj=1;Indxj<lm-1;Indxj++)
                    for(Indxi=1;Indxi<ln-1;Indxi++)
                    {
                        p_y[(Indxj+1)*ln+1+Indxi] =
                                  (p_f       [ Indxj   *ln+Indxi  ] 
                             - p_EvpCenterWgt[ Indxj   *ln+Indxi  ] * p_y[ Indxj    *ln+Indxi  ]
                                  -p_EvpNeWgt[(Indxj-1)*ln+Indxi  ] * p_y[(Indxj-1) *ln+Indxi+1]
                                  -p_EvpNeWgt[ Indxj   *ln+Indxi-1] * p_y[(Indxj+1) *ln+Indxi-1]
                                  -p_EvpNeWgt[(Indxj-1)*ln+Indxi-1] * p_y[(Indxj-1) *ln+Indxi-1])*p_InvEvpNeWgt[Indxj*ln+Indxi];

                    }
                    for(Indxi=0;Indxi<ln-2;Indxi++)
                         p_r[Indxi]= p_y[ln*(lm-1)+ Indxi+2];
      

                    for(Indxi=ln-2;Indxi<l;Indxi++)
                          p_r[Indxi] = p_y[ln-1  + (lm -1 + (ln - 3 - Indxi) )*ln];  

                    for(Indxj=0;Indxj<lm-2;Indxj++)
                    for(kk=0;kk<l;kk++)
                          p_y[1+(lm-Indxj-2)*ln]  = p_y[1+(lm-Indxj-2)*ln] + p_EvpRinv[Indxj*l+kk] * p_r[kk]; 

                    for( Indxi = 0;Indxi<ln-3;Indxi++)
                    for( kk = 0;kk<l;kk++)
                          p_y[Indxi+2 + ln*1]     = p_y[Indxi+2+ ln*1]     + p_EvpRinv[kk+(lm-2+Indxi)*l] * p_r[kk];  


                    for(Indxj=1;Indxj<lm-2;Indxj++)
                    for(Indxi=1;Indxi<ln-2;Indxi++)
                      {   p_y[Indxi+1+ln*(Indxj+1)]  = (p_f[Indxi+ln*Indxj]- p_EvpCenterWgt[Indxi+ln*Indxj] * p_y[Indxi+ln  *Indxj    ]
                                - p_EvpNeWgt [Indxi+ln  *(Indxj-1)]   * p_y[Indxi+1+ln*(Indxj-1)]
                                - p_EvpNeWgt [Indxi-1+ln* Indxj   ]   * p_y[Indxi-1+ln*(Indxj+1)]
                                - p_EvpNeWgt [Indxi-1+ln*(Indxj-1)]   * p_y[Indxi-1+ln*(Indxj-1)] ) *p_InvEvpNeWgt[Indxi+ln*Indxj];
                      }
                    for(Indxj=0;Indxj<lm-2;Indxj++)
                    for(Indxi=0;Indxi<ln-2;Indxi++)
                         p_work1[Indxi+(ln-2)*Indxj] = p_y[(Indxj+1)*ln+1+Indxi ];
        }//bi
      }//else

        athread_put(0,p_work1       ,R  +offset      ,8*(ln-2)*(lm-2), (void *)&put_reply,8*(nx-(ln-2)),8*(ln-2));

        while(put_reply!=1);
        asm volatile("memb"); 
}
   
                  

