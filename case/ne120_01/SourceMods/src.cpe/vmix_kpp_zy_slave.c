#include<slave.h>
#include"vmix_kpp_zy_struct.h"
#include<math.h>
#include"cpe_print.h"
#include<stdio.h>
#include"math_data.h"
//#define   CPE_N   32
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))


#define ALIGNED(addr) ((((unsigned long)(addr)>>5)+1)<<5)


void fun1_blmix(struct kpp_ddmix *s2)
{

     int myid;   
     double       zeta_m = -0.2,
                  zeta_s = -1.0,
                  c_m    =  8.38,
                  c_s    =  98.96,
                  a_m    =  1.26,
                  a_s    = -28.86,
                  vonkar =  0.4,
                  eps    =  1.0e-10,
                  c0     =  0.0,                 
                  c1     =  1.0,
                  c3     =  3.0,  
                  c5     =  5.0,
                  c16    =  16.0,
                  p33    =  1.0/3.0,
                  p5     =  0.500,
                  p25    =  0.250,
                  p125   =  0.125,
                  epssfc =  0.1;    


     double *SIGMA=NULL  ,*HBLT=NULL  ,*BFSFC=NULL  , *ZETAH=NULL ,*USTAR=NULL  ,*ZETA=NULL  ,*WM=NULL  ,*WS=NULL;
     double *p_SIGMA=NULL,*p_HBLT=NULL,*p_BFSFC=NULL,*p_ZETAH=NULL,*p_USTAR,*p_ZETA=NULL,*p_WM=NULL,*p_WS=NULL; 

     int size,nx,ny,km,task,task_offset,offset,task_num,i,j,tindex,threads_num,task_pos,task_size;
     volatile int   get_reply,put_reply;
     struct kpp_ddmix local_s2;
     double tmp3;

     myid = athread_get_id(-1);

   //  if(myid >=CPE_N){ return;}
     

    double pow_local_data[pow_data_len];
     athread_syn(ARRAY_SCOPE, 0xffff);
     if (myid == 0){
      get_reply = 0;
      athread_get(BCAST_MODE, pow_data, pow_local_data, pow_data_bytes, (void *)&get_reply, 0xff, 0, 0);
     while (get_reply != 1);
     asm volatile("memb");
     }
    athread_syn(ARRAY_SCOPE, 0xffff);
    pow_data_local_ptr = pow_local_data;
     
     get_reply=0;
     athread_get(0,s2,&local_s2,sizeof(local_s2),(void *)&get_reply,0,0,0);//kpp_ddmix ->local_s1
     while(get_reply!=1);
     asm volatile("memb");
     nx       = local_s2.param[0];
     ny       = local_s2.param[1];
     threads_num   = local_s2.param[2];
    // threads_num = CPE_N; 

     HBLT     = local_s2.arr[0];
     BFSFC    = local_s2.arr[1];
     USTAR    = local_s2.arr[2];
     WM       = local_s2.arr[3]; 
     WS       = local_s2.arr[4];

     double  MEM[700];
     task= nx*ny;
     task_num=task/threads_num;
     if(myid<task%threads_num)
     {
        task_num=task_num+1;
        task_pos=task_num*myid;
     }
     else
     {
        task_pos=(task%threads_num)+task_num*myid;
     }
    
      
     if(task_num>0)
     {

        task_size= task_num;

        p_HBLT =MEM,
        p_BFSFC=MEM+task_size,
        p_USTAR=MEM+2*task_size,
        p_ZETAH=MEM+3*task_size,
        p_ZETA =MEM+4*task_size,
        p_WM   =MEM+5*task_size,
        p_WS   =MEM+6*task_size;
      
        get_reply=0;
        athread_get(0,HBLT  +task_pos   ,p_HBLT ,8*task_num,(void *)&get_reply,0,0,0);
        athread_get(0,BFSFC +task_pos   ,p_BFSFC,8*task_num,(void *)&get_reply,0,0,0);
        athread_get(0,USTAR +task_pos   ,p_USTAR,8*task_num,(void *)&get_reply,0,0,0);
        while(get_reply!=3);
        asm volatile("memb"); 
       

     for(i=0;i<task_size;i++)
     {
        
  
           tmp3  = *(p_USTAR+i) *   *(p_USTAR+i) *   *(p_USTAR+i);
           p_ZETAH[i] = epssfc  * p_HBLT[i] * vonkar * p_BFSFC[i];
           p_ZETA[i]  = p_ZETAH[i]  / (tmp3 + eps);
              
           if (p_ZETA[i] >= c0) 
                p_WM[i] = vonkar*p_USTAR[i]/(c1 + c5*p_ZETA[i]);
              else if (p_ZETA[i] >= zeta_m) 
                p_WM[i] = vonkar*p_USTAR[i]*pow((c1 - c16*p_ZETA[i]),p25);
           else
                p_WM[i] = vonkar*pow((a_m*tmp3-c_m*p_ZETAH[i]),p33) ;
       
           if (p_ZETA[i] >= c0) 
               p_WS[i] = vonkar*p_USTAR[i]/(c1 + c5*p_ZETA[i]);
           else if (p_ZETA[i] >= zeta_s) 
               p_WS[i] = vonkar*p_USTAR[i]*sqrt(c1 - c16*p_ZETA[i]);
           else 
               p_WS[i] = vonkar*pow((a_s*tmp3-c_s*p_ZETAH[i]),p33);    
   
      }
        put_reply=0;
        athread_put(0,p_WM   ,WM    +task_pos  ,8*task_size, (void *)&put_reply,0,0 );
        athread_put(0,p_WS   ,WS    +task_pos  ,8*task_size, (void *)&put_reply,0,0 ) ;
        while(put_reply!=2);
         asm volatile("memb");
    
}//if
    pow_data_local_ptr = pow_data;
}
   
void fun2_blmix(struct kpp_ddmix *s2)
{
/*
     int myid;
     double *VISC,*BLMC,*GHAT,*VDC;
     int *KBL, *p_KBL;
     double *p_VISC,*p_BLMC,*p_GHAT,*p_VDC; 
     int size,nx,ny_block,km,task,task_offset,offset,task_num,i,j,tindex,mpe_id,task_rest;
     volatile int   get_reply,put_reply;
     struct kpp_ddmix local_s2;

     myid = athread_get_id(-1);

     if(myid >=CPE_N){ return;}

     get_reply=0;
     athread_get(0,s2,&local_s2,sizeof(local_s2),(void *)&get_reply,0,0,0);//kpp_ddmix ->local_s1
     nx       = local_s2.param[0];
     ny_block = local_s2.param[1];
     km   = local_s2.param[2];

     VISC    = local_s2.arr[0];
     BLMC    = local_s2.arr[1];
     GHAT    = local_s2.arr[2];
     VDC     = local_s2.arr[3];
     KBL     = local_s2.arr_int[0];

     taskall=ny_block*nj*ni;
     task_num=task_all/threads_num;
     if(myid<task_all%threads_num)
     {
        task_num=task_num+1;
        task_pos=task_num*myid;
     }
     else
     {
        task_pos=(task_all%threads_num)+task_num0*myid;
     }

if(task_num>0)
{

     for(k=0;k<km;k++)
     {

     }

     task_tile=task_num/(mx_usize/usize3);
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

      for(l=0;l<task_tile;l++)
      {

       if(l<task_tile)
        task_size1 = task_size;
       else
        task_size1 = last_size;
     
         for(i=0;i<task_size;l+)
         { 
          
*/


/*
     do k=1,km
       do j=1,ny_block
        do i=1,nx_block
           if (k < KBL(i,j)) then
              VISC(i,j,k)  = BLMC(i,j,k,1)
              VDC(i,j,k,2) = BLMC(i,j,k,2)
              VDC(i,j,k,1) = BLMC(i,j,k,3)
           else
              GHAT(i,j,k) = c0
           endif
        end do
        end do
     enddo
*/

}

void fun5_blmix(struct kpp_ddmix *s2)
{

     double       zeta_m = -0.2,
                  zeta_s = -1.0,
                  c_m    =  8.38,
                  c_s    =  98.96,
                  a_m    =  1.26,
                  a_s    = -28.86,
                  vonkar =  0.4,
                  eps    =  1.0e-10,
                  c0     =  0.0,                 
                  c1     =  1.0,
                  c2     =  2.0,
                  c3     =  3.0,  
                  c5     =  5.0,
                  c16    =  16.0,
                  p33    =  1.0/3.0,
                  p5     =  0.500,
                  p25    =  0.250,
                  p125   =  0.125,
                  epssfc =  0.1;    


     int i,k,offset,koffset,k1offset,joffset,offset1,task_size1,task_tile,tile_num,tile_num1,usize_blmc,k2offset;
     int myid,ny,nx,bid,km,task_all,task_num,usize,usize1,tile,ni,nj,task_pos,task_num0,task_size,tail_size,l,usize2,joffset1,usize3;
     volatile int   get_reply,put_reply;
     struct kpp_ddmix local_s2;
     double *pkg=NULL,*p_DZT=NULL,*p_HBLT=NULL,*p_USTAR=NULL,*p_BFSFC=NULL,*p_GAT1=NULL,*p_DAT1=NULL,*p_GHAT=NULL,*p_STABLE=NULL;
     double *p_BLMC=NULL; 
     int threads_num,max_usize,proc,toffset1,toffset2,joffset2,joffset3; 
     double *BLMC,*DZT,*HBLT,*BFSFC,*USTAR,*BLMS,*GAT1,*DAT1,*STABLE,*GHAT;
     double *zgrid,*hwide;
     double SIGMA,F1,tmp3,ZETAH,ZETA,WS,WM,cg;
     double *p_test=NULL;
     myid = athread_get_id(-1);
     int ntimes=0;
 //    double ldm_arr[6000];
    double pow_local_data[pow_data_len];
//     pow_data_local_ptr = pow_local_data; 
     athread_syn(ARRAY_SCOPE, 0xffff);
     if (myid == 0){
      get_reply = 0;
      athread_get(BCAST_MODE, pow_data, pow_local_data, pow_data_bytes, (void *)&get_reply, 0xff, 0, 0);
     while (get_reply != 1);
     asm volatile("memb");
     }
    athread_syn(ARRAY_SCOPE, 0xffff);
    pow_data_local_ptr = pow_local_data;


     get_reply=0;
     athread_get(0,s2,&local_s2,sizeof(local_s2),(void *)&get_reply,0,0,0);//kpp_ddmix ->local_s1

     double ldm_arr[5900]; //32kb
     
     while(get_reply!=1);
     asm volatile("memb");
     
     nx       = local_s2.param[0];
     ny       = local_s2.param[1];
     km       = local_s2.param[2];
     threads_num = local_s2.param[3];
     proc     = local_s2.param[4];

     cg       = local_s2.dparam[0]; 
     p33       = local_s2.dparam[1];

    
     DZT   = local_s2.arr[0];
     HBLT  = local_s2.arr[1]; 
     
     BFSFC = local_s2.arr[2];
     USTAR = local_s2.arr[3];
     BLMC  = local_s2.arr[4];
     GAT1  = local_s2.arr[5];
     DAT1  = local_s2.arr[6];
     STABLE= local_s2.arr[7];
     GHAT  = local_s2.arr[8];

     pkg     = local_s2.std[0];


//   double ldm_arr[6100]; //48.8k
     
     get_reply=0;
     athread_get(0,pkg,ldm_arr,8*(km+1)*2,(void *)&get_reply,0,0,0);
//      while(get_reply!=1);
//     asm volatile("memb");
     max_usize=1000;//44-45k
     usize  = km+2;   //km+2
     usize2 = km;
     usize3 = 3;
     usize1 = 1;
     usize_blmc= 3*km;

     nj=ny;
     ni=nx;

     task_all= nj*ni;
     task_num=task_all/threads_num;
     if(myid<task_all%threads_num)
     {
        task_num=task_num+1;
        task_pos=task_num*myid;
     }
     else
     {
        task_pos=(task_all%threads_num)+task_num*myid;
     }
 

     while(get_reply!=1);
     asm volatile("memb");

     zgrid = ldm_arr;
     hwide = ldm_arr +km+1;
/*
     if(ntimes==0&&myid==0&&proc==0)
     {       ntimes++;
            for(i=0;i<km+1;i++)
            {
               cpe_printf("   zgird[%d]= %lf  hwide[%d]=%lf  ",i,*(zgrid+i),i,*(hwide+i)); 
            }  
     }
*/
     if(task_num>0)
     { 
     	 tile_num  =task_num/(max_usize/usize2);  //tile ~20   1.5k
        
         if( tile_num >0)
         {
            task_size = max_usize / usize2;
            tail_size = task_num  % task_size;
         }
         else
         {
            task_size = task_num;
            tail_size = task_size;
         }

         // 322* task_size+130     19 task_size =50k
         p_DZT   = ldm_arr + 2*(km+1) ;//usize
         p_HBLT  = p_DZT   + usize*(task_size);//usize1
      
         p_BFSFC = p_HBLT  + usize1*task_size;//usize1

         p_USTAR = p_BFSFC + usize1*task_size;//usize1
         p_DAT1  = p_USTAR + usize1*task_size; //usize3
         p_GAT1  = p_DAT1  + usize3*task_size; //usize3
         p_STABLE= p_GAT1  + usize3*task_size; //usize1
         p_GHAT  = p_STABLE+ usize1*task_size;  //usize2
         p_BLMC  = p_GHAT  + usize2*task_size; //usize_blmc
         p_test  = p_BLMC  + usize_blmc*task_size; 
         if(tail_size>0)
             tile_num1=tile_num+1;
          else
             tile_num1=tile_num;

 
         for(l=0;l<tile_num1;l++)
         {   
                if(l<tile_num)
                    task_size1=task_size;
                else
                    task_size1=tail_size;                 
   

        	offset = task_pos;
        	get_reply=0;

                athread_get(0,  DAT1+offset,  p_DAT1  ,8*usize3*task_size1,(void *)&get_reply,//0, 0, 0);
                            0,  8*(ny*nx-task_size1*usize1),8*task_size1*usize1);

                athread_get(0,  GAT1+offset,  p_GAT1  ,8*usize3*task_size1,(void *)&get_reply,//0, 0, 0);
                            0,  8*(ny*nx-task_size1*usize1),8*task_size1*usize1);
             
               athread_get(0,  DZT +offset,   p_DZT   ,8*usize *task_size1 ,(void *)&get_reply,
                            0,  8*(ny*nx-task_size1*usize1),8*task_size1*usize1); 

               athread_get(0,  HBLT  +offset,  p_HBLT   ,8*usize1*task_size1,(void *)&get_reply,0, 0, 0);
               athread_get(0,  BFSFC +offset,  p_BFSFC  ,8*usize1*task_size1,(void *)&get_reply,0, 0, 0);
               athread_get(0,  USTAR +offset,  p_USTAR  ,8*usize1*task_size1,(void *)&get_reply,0, 0, 0);
               athread_get(0,  STABLE+offset,  p_STABLE ,8*usize1*task_size1,(void *)&get_reply,0, 0, 0);

        	while(get_reply!=7);
        	asm volatile("memb");

        	for(i=0;i<task_size1;i++)
        	{

                        joffset  =  i*usize1;  //i,j
                        joffset1 = (  task_size1 +i)*usize1;  //i,j,2
                        joffset2 = (2*task_size1 +i)*usize1;  //i,j,3

			for(k=0;k<km;k++)
                 	{
				koffset   = ( k        * task_size1+i)*usize1;//i,j,k
                                k1offset  = ((k+1)     * task_size1+i)*usize1;//i,j,k-1 
                                toffset1  = ((k +km)    * task_size1+i)*usize1;
                                toffset2  = ((k +2*km)  * task_size1+i)*usize1;
                                
                                 if(k>0) 
                                     SIGMA  = (-zgrid[k-1] + p5 * p_DZT[koffset] + p_DZT[k1offset]) / p_HBLT[joffset];
                                 else
                                     SIGMA  = (-zgrid[k]   + p5 * hwide[k]) / p_HBLT[joffset] ;


                                 F1     = MIN(SIGMA,epssfc);
                  

                                 ZETAH  = F1 * p_HBLT[joffset] * vonkar * p_BFSFC[joffset];
                                 tmp3   = p_USTAR[joffset]  *  p_USTAR[joffset] * p_USTAR[joffset]; 
                                 ZETA   = ZETAH   / ( tmp3 + eps);
	

                                 if (ZETA >= c0) 
                                 	WM = vonkar * p_USTAR[joffset]/(c1 + c5*ZETA);
                                 else if (ZETA >= zeta_m) 
                                 	WM = vonkar * p_USTAR[joffset]*pow((c1 - c16*ZETA),p25);
                                 else
                			WM = vonkar * pow((a_m * tmp3-c_m*ZETAH),p33) ;
       
           			 if (ZETA >= c0) 
               			 	WS = vonkar * p_USTAR[joffset]/(c1 + c5*ZETA);
           			 else if (ZETA >= zeta_s) 
               			 	WS = vonkar * p_USTAR[joffset]*sqrt(c1 - c16*ZETA);
           			 else 
               			 	WS = vonkar * pow((a_s *tmp3 -c_s*ZETAH),p33);    

                                 p_BLMC[koffset] = p_HBLT[joffset]  *  WM  *  SIGMA    *       
                         		(c1  +  SIGMA  *  ((SIGMA-c2)        + 
                       	      	 	(c3  -  c2*SIGMA)  *  p_GAT1[joffset]    +    
                           		(SIGMA-c1)  *  p_DAT1[joffset]));
           			 p_BLMC[toffset1] = p_HBLT[joffset]  *  WS  *  SIGMA  *       
                           		(c1  +  SIGMA  *  ((SIGMA  -  c2)  +  
                           		(c3  -  c2  *  SIGMA)  *  p_GAT1[joffset1]  +    
                           		(SIGMA  -  c1)  *  p_DAT1[joffset1]));
               			 p_BLMC[toffset2] = p_HBLT[joffset]  *  WS  *  SIGMA  *       
                           		(c1  +  SIGMA  *  ((SIGMA  -  c2)   + 
                           		(c3  -  c2  *  SIGMA)  *  p_GAT1[joffset2] +    
                           		(SIGMA  -  c1)  *  p_DAT1[joffset2]));
                                 p_GHAT[koffset] = (c1 - p_STABLE[joffset])  *  cg  /  (WS  *  p_HBLT[joffset]  +  eps);
                        }//k

        	}//i
              
                  put_reply=0;
                  offset = task_pos;
                  athread_put(0,p_GHAT   ,GHAT    +offset  ,8*task_size1*usize2     ,(void*)&put_reply,8*(nx*ny-task_size1*usize1),8*task_size1*usize1);
                  athread_put(0,p_BLMC   ,BLMC    +offset  ,8*task_size1*usize_blmc ,(void*)&put_reply,8*(nx*ny-task_size1*usize1),8*task_size1*usize1) ;
                  while(put_reply!=2);
      
                  asm volatile("memb");
                  task_pos = task_pos + task_size1; 
    	}//l
   }//if

     pow_data_local_ptr = pow_data;
}//fun

/*
do k = 1,km      
        do j=1,ny_block
        do i=1,nx_block
           if (k > 1) then
              SIGMA(i,j) = (-zgrid(k-1) + p5*DZT(i,j,k-1,bid) +  &
                       DZT(i,j,k,bid)) / HBLT(i,j)
           else
              SIGMA(i,j) = (-zgrid(k) + p5*hwide(k)) / HBLT(i,j)
           end if
           F1(i,j) = min(SIGMA(i,j),epssfc)

        ZETAH(i,j) = F1(i,j)*HBLT(i,j)*vonkar*BFSFC(i,j)
        ZETA(i,j)  = ZETAH(i,j)/(USTAR(i,j)**3 + eps)


           if (ZETA(i,j) >= c0) then ! stable region
               WM(i,j) = vonkar*USTAR(i,j)/(c1 + c5*ZETA(i,j))
           else if (ZETA(i,j) >= zeta_m) then
               WM(i,j) = vonkar*USTAR(i,j)*(c1 - c16*ZETA(i,j))**p25
           else
               WM(i,j) = vonkar*(a_m*(USTAR(i,j)**3)-c_m*ZETAH(i,j))**p33
           endif
          
           if (ZETA(i,j) >= c0) then
            WS(i,j) = vonkar*USTAR(i,j)/(c1 + c5*ZETA(i,j))
          else if (ZETA(i,j) >= zeta_s) then
            WS(i,j) = vonkar*USTAR(i,j)*SQRT(c1 - c16*ZETA(i,j))
          else
            WS(i,j) = vonkar*(a_s*(USTAR(i,j)**3)-c_s*ZETAH(i,j))**p33
          endif
           
           BLMC(i,j,k,1) = HBLT(i,j)*WM(i,j)*SIGMA(i,j)*       &
                           (c1 + SIGMA(i,j)*((SIGMA(i,j)-c2) + &
                           (c3-c2*SIGMA(i,j))*GAT1(i,j,1) +    &
                           (SIGMA(i,j)-c1)*DAT1(i,j,1))) 
           BLMC(i,j,k,2) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*       &
                           (c1 + SIGMA(i,j)*((SIGMA(i,j)-c2) + &
                           (c3-c2*SIGMA(i,j))*GAT1(i,j,2) +    &
                           (SIGMA(i,j)-c1)*DAT1(i,j,2)))
           BLMC(i,j,k,3) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*       &
                           (c1 + SIGMA(i,j)*((SIGMA(i,j)-c2) + &    
                           (c3-c2*SIGMA(i,j))*GAT1(i,j,3) +    &
                           (SIGMA(i,j)-c1)*DAT1(i,j,3)))

           GHAT(i,j,k) = (c1-STABLE(i,j))* cg/(WS(i,j)*HBLT(i,j) +eps)
        end do
        end do

     end do
*/


void fun8_blmix(struct kpp_ddmix *s2)
{
/*
     int myid;
     double *VISC,*BLMC,*GHAT,*VDC;
     int *KBL, *p_KBL;
     double *p_VISC,*p_BLMC,*p_GHAT,*p_VDC; 
     int size,nx,ny_block,km,task,task_offset,offset,task_num,i,j,tindex,mpe_id,task_rest;
     volatile int   get_reply,put_reply;
     struct kpp_ddmix local_s2;

     myid = athread_get_id(-1);

     if(myid >=CPE_N){ return;}

     get_reply=0;
     athread_get(0,s2,&local_s2,sizeof(local_s2),(void *)&get_reply,0,0,0);//kpp_ddmix ->local_s1
     nx       = local_s2.param[0];
     ny_block = local_s2.param[1];
     km   = local_s2.param[2];

     VISC    = local_s2.arr[0];
     BLMC    = local_s2.arr[1];
     GHAT    = local_s2.arr[2];
     VDC     = local_s2.arr[3];
     KBL     = local_s2.arr_int[0];

     taskall=ny_block*nj*ni;
     task_num=task_all/threads_num;
     if(myid<task_all%threads_num)
     {
        task_num=task_num+1;
        task_pos=task_num*myid;
     }
     else
     {
        task_pos=(task_all%threads_num)+task_num0*myid;
     }

if(task_num>0)
{

     for(k=0;k<km;k++)
     {

     }

     task_tile=task_num/(mx_usize/usize3);
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

      for(l=0;l<task_tile;l++)
      {

       if(l<task_tile)
        task_size1 = task_size;
       else
        task_size1 = last_size;
     
         for(i=0;i<task_size;l+)
         { 
          
*/


/*
     do k=1,km
       do j=1,ny_block
        do i=1,nx_block
           if (k < KBL(i,j)) then
              VISC(i,j,k)  = BLMC(i,j,k,1)
              VDC(i,j,k,2) = BLMC(i,j,k,2)
              VDC(i,j,k,1) = BLMC(i,j,k,3)
           else
              GHAT(i,j,k) = c0
           endif
        end do
        end do
     enddo
*/

}

