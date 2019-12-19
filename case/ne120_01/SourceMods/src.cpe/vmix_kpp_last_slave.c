#include <stdio.h>
#include <math.h>
#include <slave.h>
#include <unistd.h>
#include "vmix_kpp_last.h"

#ifdef POPCHKJN
#include "cpe_print.h"
#endif

#define SZldmint 384
#define SZldmdbl 5600
#define locmax(a,b)  ((a)>(b)?(a):(b))
#define locmin(a,b)  ((a)<(b)?(a):(b))


#define get_myrid(row)  \
asm volatile("rcsr   %0, 1" : "=r"(row))
#define get_mycid(col)  \
asm volatile("rcsr   %0, 2" : "=r"(col))
#define sync_mem   \
asm volatile("memb")
//asm volatile("memb\n\t":::"memory")



//------------------------------------------------------------------
void s_kpp_last1(struct param_kpp_last_s1 *st2_master)
{
	
// BOP:

 double locdbl[SZldmdbl];

 int locint[SZldmint];

 double *s1dd, *s2dd, *s3dd;
 int *s2di;

 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(POPCHKJN)
 volatile int num_calls=0;
#endif

// Common variables:
 int myproc,myid,myrid,mycid,threads_num;
 int task_num0,task_num,task_numj,task_numj0,
     task_pos,task_posi,task_posj,task_post, 
     task_size,task_tile,task_all,last_size,task_tile1,task_size1;
 int nx,ny,km,nt;

 int offset,offset1,retng,retnp;
 int mx_usize,usize0,usize,usize1,usize2,usize3;
 int i,j,k,n,mt2,l,ni,nj,nj1,nk,nm,kl,pij,nxy;

 int itmp,jtmp,ktmp,kltmp;


// Input variables: 
 int *KBL, *KMT, *KMU;			// nx*ny
 double conv_visc, conv_diff, bvsqcon;

// double *dz, *zt;	// km
// double *zgrid;		// 0:km+1

// double *AU0,*AUN,*AUE,*AUNE;   // nx*ny
// double *STF;		// nx*ny*nt

 double *GHAT, *DBSFC, *DBLOC;   	// nx*ny*km
 double *DZT;		// nx*ny*(0:km+1)


// Output variables:
 double *HMXL;		// nx*ny
 double *VVC;		// nx*ny*km


// Input & Output variables:
 double *VISC;		// nx*ny*(0:km+1)
 double *VDC;		// nx*ny*(0:km+1)*2
 double *KPP_SRC;	// nx*ny*km*nt
 

// Slave & Temp variables:
 struct param_kpp_last_s1 st2_slave;
 double WORK1, *WORK2, FCON, USTAR, STABLE, BFSFC;  //nx*ny

 int *m2di;
 double *m1dd, *m2dd;
// unsigned long m1dd, m2dd;

 double xtmp,ytmp,ztmp;
 int *s_KMT, *s_KMU, *s_KBL; 
 double *s_dz, *s_zt, *s_zgrid;
 double *s_AU0, *s_AUN, *s_AUE, *s_AUNE, *s_STF;
 double *s_GHAT, *s_DBLOC, *s_DBSFC, *s_DZT;
 double *s_HMXL, *s_VVC, *s_VDC, *s_VISC, *s_KPP_SRC;
 double *u_buf, *d_buf;


// Constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;

// EOP.


#if defined(POPCHKJN)
 num_calls = num_calls + 1;
#endif

 get_reply=0;
 athread_get(PE_MODE,st2_master,&st2_slave,sizeof(st2_slave),(void*)&get_reply,0,0,0);
 while(get_reply!=1);
 sync_mem;


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


 nx = st2_slave.ipts[0];
 ny = st2_slave.ipts[1];
 km = st2_slave.ipts[2];
 nt = st2_slave.ipts[3];
 myproc      = st2_slave.ipts[4];
 threads_num = st2_slave.ipts[5];


 conv_visc = st2_slave.dpts[0];
 conv_diff = st2_slave.dpts[1];
 bvsqcon   = st2_slave.dpts[2];



#if defined(POPINSJN)
 if(myproc==0 && num_calls<8 && myid==0)
    {
    cpe_printf("nx: %d, km: %d, bvsqcon:%f. \n",nx,km,bvsqcon);
    cpe_printf("st2.darr[0]_s1dd: %ld, st2.darr[1]_s2dd: %ld. \n",st2_slave.darr[0], st2_slave.darr[1]);
    cpe_printf("st2.darr[5]_DZT: %ld, st2.darr[9]_KPP_SRC: %ld. \n",st2_slave.darr[5], st2_slave.darr[9]);
    }
#endif



 task_numj0 = ny/threads_num;
 if(myid < ny%threads_num)
    {
    task_numj=task_numj0+1;
    task_posi=task_numj*myid;
    task_posj=task_numj*myid*nx;
    }
 else
    {
    task_numj=task_numj0;
    task_posi=ny%threads_num + task_numj0*myid;
    task_posj=(ny%threads_num + task_numj0*myid) *nx;
    }
 usize0 = task_numj*nx;    // total size for 2d-array.






if(task_numj>0)
{	
 
 mx_usize = MAXUSZ;    // MAXUSZ SHOULD be larger than usize0! 
 if((km-1)*usize0<=mx_usize)
   {
    ni = 1;
    nj = 1;
    nj1= 1;
    usize  = (km-1)*usize0;
    usize1 = (km+2-1)*usize0;	
   }
 else if(usize0<=mx_usize)  // general case!
   {
    ni = 1;
    nj = km-1; 
    nj1= km+2-1;
    usize  = usize0;	
    usize1 = usize0;
   }
 else                   // extreme case
   {
    // MUST to modify! 
    // when km is larger than MAXUSZ.
   }

 
 
 task_num = nj*ni;

 task_tile=task_num/(mx_usize/usize);       // use usize as beseline
 // per tile has (task_num/task_tile) tasks with each task ~(8*usize).
 if(task_tile>0)
   { 
   task_size = mx_usize/usize;
   last_size = task_num%task_size;
   }
 else
   {
   task_size = task_num;
   last_size = task_num;
   }  


 if(last_size>0)
   task_tile1=task_tile+1;
 else 
   task_tile1=task_tile;
    



 m1dd = st2_slave.darr[0];
// m1dd = (unsigned long)st2_slave.darr[0];
// dz    =  m1dd;
// zt    =  m1dd +   km;
// zgrid =  m1dd + 2*km;

 m2di = st2_slave.iarr[0];
// KBL =  m2di;
// KMT =  m2di +   nxy;
// KMU =  m2di + 2*nxy;

 m2dd = st2_slave.darr[1];
// AU0  = m2dd;
// AUN  = m2dd +   nxy;
// AUE  = m2dd + 2*nxy;
// AUNE = m2dd + 3*nxy;
// // BFSFC= m2dd + 4*nxy;

 GHAT  = st2_slave.darr[2];
 DBSFC = st2_slave.darr[3];
 DBLOC = st2_slave.darr[4];
 DZT   = st2_slave.darr[5];
 VISC  = st2_slave.darr[6];
 VVC   = st2_slave.darr[7];
 VDC   = st2_slave.darr[8];



 s1dd  = locdbl; 

#if defined(POPINSJN)
 if(myproc==0 && num_calls<8 && myid==0)
    {
    cpe_printf("nx: %d, km: %d, size:%d. \n",nx,km,SIZEDBL*(3*km+2));
    cpe_printf("m1dd: %ld, s1dd: %ld. \n",m1dd, s1dd);
    }
#endif

 get_reply = 0;
 athread_get(PE_MODE,m1dd,s1dd,SIZEDBL*(3*km+2),(void*)&get_reply,0,0,0);
 sync_mem;
 while(get_reply!=1);

#if defined(POPINSJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("km: %d, retng:%d. \n",km,retng);
#endif


 s_dz  = s1dd; 
 s_zt  = s1dd + km; 
 s_zgrid  = s1dd + 2*km; 


 nxy = nx*ny;
 


 s2di = locint;

#if defined(POPINSJN)
 if(myproc==0 && num_calls<8 && myid==0)
    {
    cpe_printf("task_posj: %d, size:%d. \n",task_posj,SIZEINT*(3*usize0));
    cpe_printf("m2di: %ld, s2di: %ld. \n",m2di, s2di);
    }
#endif

 get_reply = 0;
 athread_get(PE_MODE,m2di+task_posj,s2di,SIZEINT*3*usize0,(void*)&get_reply,0,SIZEINT*(nxy-usize0),SIZEINT*usize0);
 sync_mem;
 while(get_reply!=1);

 s_KBL = s2di;
 s_KMT = s2di +   usize0;
 s_KMU = s2di + 2*usize0;


 

#if defined(POPINSJN)
 if(myproc==0 && num_calls<8 && myid==0)
    {
    cpe_printf("nx: %d, km: %d, nt:%d. \n",nx,km,nt);
    cpe_printf("dz[0]: %f, dz[15]: %f. \n",s_dz[0],s_dz[15]);
    cpe_printf("zt[0]: %f, zt[15]: %f. \n",s_zt[0],s_zt[15]);
    cpe_printf("zgrid[0]: %f, zgrid[15]: %f. \n",s_zgrid[0],s_zgrid[15]);
    }
#endif


#if defined(POPINSJN)
 if(myproc==0 && num_calls<8 && myid==0)
    {
    cpe_printf("nx: %d, km: %d, size:%d. \n",nx,km,SIZEDBL*(6+nt)*usize0);
    cpe_printf("m2dd: %ld, s2dd: %ld. \n",m2dd, s2dd);
    }
#endif


 s2dd = locdbl + 4*km;

 get_reply = 0;
 athread_get(PE_MODE,m2dd+task_posj,s2dd,SIZEDBL*4*usize0,(void*)&get_reply,0,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);
 while(get_reply!=1);
 sync_mem;

 s_AU0  = s2dd;
 s_AUN  = s2dd +   usize0;
 s_AUE  = s2dd + 2*usize0;
 s_AUNE = s2dd + 3*usize0;

 WORK2  = s2dd + 4*usize0;
 u_buf  = s2dd + 5*usize0;
 d_buf  = s2dd + 6*usize0;


 
 s3dd = locdbl + 4*km + 8*usize0; //14*task_size*usize1

 s_GHAT  = s3dd;
 s_DBSFC = s3dd +   usize1*task_size;
 s_DBLOC = s3dd + 2*usize1*task_size;
 s_VVC   = s3dd + 3*usize1*task_size;
 s_DZT   = s3dd + 4*usize1*task_size;
 s_VISC  = s3dd + 6*usize1*task_size;
 s_VDC   = s3dd + 8*usize1*task_size;
 




#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last1_ task_tile1: %d, task_size: %d, last_size:%d. \n",task_tile1,task_size,last_size);
#endif


   task_pos = task_posj;
//  for(k=0;k<km-1;k++)
   for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
  	task_size1 = task_size;
      else  
        task_size1 = last_size;


      // get: DBLOC, DZT, VISC, VDC. 
      get_reply = 0;
      athread_get(PE_MODE,DBLOC+task_pos,s_DBLOC,SIZEDBL*usize0*task_size1,(void*)&get_reply,0,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);

      athread_get(PE_MODE,DZT+task_pos,s_DZT,SIZEDBL*usize0*(task_size1+2),(void*)&get_reply,0,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);

      athread_get(PE_MODE,VISC+task_pos,s_VISC,SIZEDBL*usize0*(task_size1+2),(void*)&get_reply,0,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);

      athread_get(PE_MODE,VDC+task_pos,s_VDC,SIZEDBL*usize0*(task_size1+2),(void*)&get_reply,0,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);
      athread_get(PE_MODE,VDC+(km+2)*nxy+task_pos,s_VDC+usize0*(task_size1+2),SIZEDBL*usize0*(task_size1+2),(void*)&get_reply,0,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);
      while(get_reply!=5);
      sync_mem;

#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
   cpe_printf("last1 get ok with l: %d. \n", l);
#endif


		 
      for(kl=0;kl<task_size1;kl++)
        {  
        kltmp = l*task_size + kl;
/*
        for(j=0;j<task_numj;j++)
          for(i=0;i<nx;i++)
            {
            pij = j*nx+i;
*/
        for(pij=0;pij<usize0;pij++)
            {
            WORK1 = s_DBLOC[kl*usize0+pij]/(p5*(s_DZT[(kl+1)*usize0+pij]+s_DZT[(kl+1+1)*usize0+pij]));

            if(bvsqcon != c0)
              {
              xtmp = locmax(WORK1,bvsqcon); 
              xtmp = c1 - xtmp/bvsqcon;
              WORK1 = locmin(xtmp, c1);
              xtmp = c1 - WORK1*WORK1;
              FCON = xtmp*xtmp*xtmp;
              }
            else
              {
              if(WORK1>c0)
                FCON = c0;
              else
                FCON = c1;
              }

            if((kltmp+1) >= s_KBL[pij])
              {
              s_VISC[(kl+1)*usize0+pij] = s_VISC[(kl+1)*usize0+pij] + conv_visc*FCON;
              s_VDC[(kl+1)*usize0+pij]  = s_VDC[(kl+1)*usize0+pij]  + conv_diff*FCON;
              s_VDC[((task_size1+2)+(kl+1))*usize0+pij] = s_VDC[((task_size1+2)+(kl+1))*usize0+pij] + conv_diff*FCON;
              }

            if((kltmp+1) >= s_KMT[pij])
              {
              s_VISC[(kl+1)*usize0+pij] = c0;
              s_VDC[(kl+1)*usize0+pij]  = c0;
              s_VDC[((task_size1+2)+(kl+1))*usize0+pij] = c0;
              }
            }


        // call tgrid_to_ugrid(WORK2,VISC(:,:,k),bid)

        m4pwgt_work2_s(WORK2,s_VISC+(kl+1)*usize0,s_AU0,s_AUN,s_AUE,s_AUNE,u_buf,d_buf,nx,ny,task_numj,threads_num,myrid,mycid,myid);
        if((task_posi+task_numj)==ny)
          {
          for(i=0;i<nx;i++)
            WORK2[(task_numj-1)*nx+i] = c0;
          }
        for(j=0;j<task_numj;j++)
          WORK2[j*nx+nx-1] = c0;
     
        //!--expand tgrid_to_ugrid:
        /*
        for(j=0;j<ny-1;j++)
         for(i=0;i<nx-1;i++)
            {
            pij = j*nx+i;
            WORK2[pij] = AU0[pij]*VISC[(k+1)*nxy+pij]+ 
                     AUN[pij]*VISC[(k+1)*nxy+pij+nx] + 
                     AUE[pij]*VISC[(k+1)*nxy+pij+1]  + 
                     AUNE[pij]*VISC[(k+1)*nxy+pij+nx+1];
            }

        for(i=0;i<nx;i++)
          WORK2[(ny-1)*nx+i] = c0;
        for(j=0;j<ny;j++)
          WORK2[j*nx+nx-1] = c0;
        */
        //--expand tgrid_to_ugrid finished.


        for(pij=0;pij<usize0;pij++)
          {
          if((kltmp+1) < s_KMU[pij])
            s_VVC[kl*usize0+pij] = WORK2[pij];
          else
            s_VVC[kl*usize0+pij] = c0;
          }

       }
       //  end of for kl.

#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
   cpe_printf("last1 comp ok with l: %d. \n", l);
#endif



      // put: VVC, VDC.
      put_reply = 0;
      athread_put(PE_MODE,s_VVC,VVC+task_pos,SIZEDBL*usize0*task_size1,(void*)&put_reply,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);

      athread_put(PE_MODE,s_VDC+usize0,VDC+nxy+task_pos,SIZEDBL*usize0*task_size1,(void*)&put_reply,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);
      athread_put(PE_MODE,s_VDC+usize0+usize0*(task_size1+2),VDC+nxy+(km+2)*nxy+task_pos,SIZEDBL*usize0*task_size1,(void*)&put_reply,SIZEDBL*(nxy-usize0),SIZEDBL*usize0);
      while(put_reply!=3);
      sync_mem;



#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last1 put ok with l: %d. \n", l);
#endif


      task_pos = task_pos + task_size1*nxy;	
 
      }
      // end of for l.



   for(pij=0;pij<usize0;pij++)
        {
        s_VVC[pij] = c0;
        s_VDC[pij] = c0;
        s_VDC[usize0+pij] = c0;
        }
   // k=km.

   // put: VVC, VDC.
   task_pos = task_posj + (km-1)*nxy;

   put_reply = 0;
   athread_put(PE_MODE,s_VVC,VVC+task_pos,SIZEDBL*usize0,(void*)&put_reply,0,0);

   athread_put(PE_MODE,s_VDC,VDC+task_pos+nxy,SIZEDBL*usize0,(void*)&put_reply,0,0);
   athread_put(PE_MODE,s_VDC+usize0,VDC+task_pos+(km+2)*nxy+nxy,SIZEDBL*usize0,(void*)&put_reply,0,0);
   while(put_reply!=3);
   sync_mem;


#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last1 finished. \n");
#endif



}
// end of if(task_numj>0)



return;
}








//---------------:---------------------------------------------------
void s_kpp_last2(struct param_kpp_last_s1 *st2_master)
{
	
// BOP:

 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s1dd, *s2dd, *s3dd;

 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(POPTESTJN)
 volatile int num_calls=0;
#endif

// Common variables:
 int myproc,myid,myrid,mycid,threads_num;
 int task_num0,task_num,task_numj,task_numj0,
     task_pos,task_pos1,task_posi,task_posj,task_post, 
     task_size,task_tile,task_all,last_size,task_tile1,task_size1;
 int nx,ny,km,nt;

 int offset,offset1;
 int mx_usize,usize0,usize,usize1,usize2,usize3;
 int i,j,k,n,mt2,l,ni,nj,nj1,nk,nm,kl,pij,nxy;

 int itmp,jtmp,ktmp,kltmp;


// Input variables: 
 int *KBL, *KMT, *KMU;			// nx*ny
 double conv_visc, conv_diff, bvsqcon;

 double *dz, *zt;	// km
 double *zgrid;		// 0:km+1

 double *STF;		// nx*ny*nt

 double *GHAT, *DBSFC, *DBLOC;   	// nx*ny*km
 double *DZT;		// nx*ny*(0:km+1)


// Output variables:
 double *HMXL;		// nx*ny
 double *VVC;		// nx*ny*km


// Input & Output variables:
 double *VISC;		// nx*ny*(0:km+1)
 double *VDC;		// nx*ny*(0:km+1)*2
 double *KPP_SRC;	// nx*ny*km*nt
 

// Slave & Temp variables:
 struct param_kpp_last_s1 st2_slave;
 double WORK1, *WORK2, *FCON, USTAR, STABLE, BFSFC;  //nx*ny

 int *m2di;
 double *m1dd, *m2dd;

 double xtmp,ytmp,ztmp;
 int *s_KMT, *s_KMU, *s_KBL; 
 double *s_dz, *s_zt, *s_zgrid;
 double *s_GHAT, *s_DBLOC, *s_DBSFC, *s_DZT;
 double *s_HMXL, *s_VVC, *s_VDC, *s_VISC, *s_KPP_SRC;
 double *s_STF;


// Constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;

// EOP.


#if defined(POPTESTJN)
 num_calls = num_calls + 1;
#endif

 get_reply=0;
 athread_get(PE_MODE,st2_master,&st2_slave,sizeof(st2_slave),(void*)&get_reply,0,0,0);
 while(get_reply!=1);
 sync_mem;


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


 nx = st2_slave.ipts[0];
 ny = st2_slave.ipts[1];
 km = st2_slave.ipts[2];
 nt = st2_slave.ipts[3];
 myproc      = st2_slave.ipts[4];
 threads_num = st2_slave.ipts[5];


 conv_visc = st2_slave.dpts[0];
 conv_diff = st2_slave.dpts[1];
 bvsqcon   = st2_slave.dpts[2];




 task_numj0 = (nx*ny)/threads_num;
 if(myid < (nx*ny)%threads_num)
    {
    task_numj=task_numj0+1;
    task_posi=task_numj*myid;                    // for 2-D arr
    task_pos =task_posi;                         // for 3-D arr w. km
    }
 else
    {
    task_numj=task_numj0;
    task_posi=(nx*ny)%threads_num + task_numj0*myid;    // for 2-D arr
    task_pos =task_posi;                                // for 3-D arr w. km
    }
 usize0 = km; 



if(task_numj>0)
{	
 
 mx_usize = MAXUSZ;    // MAXUSZ SHOULD be larger than usize0! 
 if(task_numj*usize0<=mx_usize)
   {
    ni = 1;
    nj = 1;
    usize  = task_numj*usize0;
    usize1 = task_numj*(usize0+2);	
   }
 else if(usize0<=mx_usize)  // general case!
   {
    ni = 1;
    nj = task_numj; 
    usize  = usize0;	
    usize1 = usize0+2;
   }
 else                   // extreme case
   {
    // MUST to modify! 
    // when km is larger than MAXUSZ.
   }

 

 task_tile = nj/(mx_usize/usize);       // use usize as beseline
 // per tile has (task_num/task_tile) tasks with each task ~(8*usize).
 if(task_tile>0)
   { 
   task_size = mx_usize/usize;
   last_size = nj%task_size;
   }
 else
   {
   task_size = nj;
   last_size = nj;
   }  


 if(last_size>0)
   task_tile1=task_tile+1;
 else 
   task_tile1=task_tile;
    


 nxy = nx*ny;
 
 m2di = st2_slave.iarr[0];
// KBL =  m2di;
// KMT =  m2di +   nxy;
// KMU =  m2di + 2*nxy;



 s2di = locint;

 get_reply = 0;
 athread_get(PE_MODE,m2di+task_posi,s2di,SIZEINT*3*task_numj,(void*)&get_reply,0,SIZEINT*(nxy-task_numj),SIZEINT*task_numj);
 while(get_reply!=1);
 sync_mem;

 s_KBL = s2di;
 s_KMT = s2di +   task_numj;
 s_KMU = s2di + 2*task_numj;



 

 m1dd = st2_slave.darr[0];
// dz    =  m1dd;
// zt    =  m1dd +   km;
// zgrid =  m1dd + 2*km;


 s1dd  = locdbl; 

 get_reply = 0;
 athread_get(PE_MODE,m1dd,s1dd,SIZEDBL*(3*km+2),(void*)&get_reply,0,0,0);
 while(get_reply!=1);
 sync_mem;
 s_dz  = s1dd; 
 s_zt  = s1dd + km; 
 s_zgrid  = s1dd + 2*km; 




 m2dd = st2_slave.darr[1];
// AU0  = m2dd;
// AUN  = m2dd +   nxy;
// AUE  = m2dd + 2*nxy;
// AUNE = m2dd + 3*nxy;
// // BFSFC= m2dd + 4*nxy;
 HMXL = m2dd + 5*nxy;
// STF  = m2dd + 6*nxy;

 s2dd = locdbl + 4*km;

 get_reply = 0;
 athread_get(PE_MODE,m2dd+task_posi,s2dd,SIZEDBL*(6+nt)*task_numj,(void*)&get_reply,0,SIZEDBL*(nxy-task_numj),SIZEDBL*task_numj);
 while(get_reply!=1);
 sync_mem;

// s_AU0;
// s_AUN;
// s_AUE;
// s_AUNE;
 s_HMXL = s2dd + 5*task_numj;
 s_STF  = s2dd + 6*task_numj;





 GHAT  = st2_slave.darr[2];
 DBSFC = st2_slave.darr[3];
 DBLOC = st2_slave.darr[4];
 DZT   = st2_slave.darr[5];
 VISC  = st2_slave.darr[6];
 VVC   = st2_slave.darr[7];
 VDC   = st2_slave.darr[8];
 KPP_SRC = st2_slave.darr[9];

 
 s3dd = locdbl + 4*km + (9+nt)*task_numj; //12*task_size*usize1

 s_GHAT  = s3dd;
 s_DBSFC = s3dd +   task_size*usize;
 s_VVC   = s3dd + 2*task_size*usize;
 s_KPP_SRC = s3dd + 3*task_size*usize;
 s_DZT   = s3dd + (3+nt)*task_size*usize;
 s_VISC  = s3dd + (5+nt)*task_size*usize;
 s_VDC   = s3dd + (7+nt)*task_size*usize;
 




#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last2_ task_tile1: %d, task_size: %d, last_size:%d. \n",task_tile1,task_size,last_size);
#endif



 //!  add ghatp term from previous computation to right-hand-side 
 //!  source term on current row

 task_pos = task_posi;
 for(l=0;l<task_tile1;l++)
   {	 
   if(l<task_tile)
      task_size1 = task_size;
   else  
      task_size1 = last_size;


   // get: GHAT, DZT, VDC. 
   get_reply = 0;
   athread_get(PE_MODE,GHAT+task_pos,s_GHAT,SIZEDBL*usize*task_size1,(void*)&get_reply,0,SIZEDBL*(nxy-task_size1),SIZEDBL*task_size1);

   athread_get(PE_MODE,DZT+task_pos,s_DZT,SIZEDBL*usize1*task_size1,(void*)&get_reply,0,SIZEDBL*(nxy-task_size1),SIZEDBL*task_size1);

   athread_get(PE_MODE,VDC+task_pos,s_VDC,SIZEDBL*usize1*task_size1,(void*)&get_reply,0,SIZEDBL*(nxy-task_size1),SIZEDBL*task_size1);
   athread_get(PE_MODE,VDC+(km+2)*nxy+task_pos,s_VDC+usize1*task_size1,SIZEDBL*usize1*task_size1,(void*)&get_reply,0,SIZEDBL*(nxy-task_size1),SIZEDBL*task_size1);
   while(get_reply!=4);
   sync_mem;

#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
   cpe_printf("last2_KPP_SRC get ok with l: %d. \n", l);
#endif

   for(n=0;n<nt;n++)
      {
      mt2 = locmin(n+1,2)-1;

      for(j=0;j<task_size1;j++)       // for 3-D arr:
        {
        pij = l*task_size + j;        // for 2-D arr:
        s_KPP_SRC[n*usize*task_size1+j] = s_STF[n*task_numj+pij]/s_dz[0]
              *(-s_VDC[(mt2*usize1+1)*task_size1+j]*s_GHAT[j]);
        for(k=1;k<km;k++)
          s_KPP_SRC[(n*usize+k)*task_size1+j] = s_STF[n*task_numj+pij]/s_DZT[(k+1)*task_size1+j] 
              *(s_VDC[(mt2*usize1+(k+1)-1)*task_size1+j]*s_GHAT[(k-1)*task_size1+j]
               -s_VDC[(mt2*usize1+(k+1))*task_size1+j]*s_GHAT[k*task_size1+j]);
        }
        // end of for(j)

#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last2_KPP_SRC comp with l: %d. \n", l);
#endif
 
      }
      // end of for(n)

   // put: KPP_SRC.
   put_reply = 0;
   athread_put(PE_MODE,s_KPP_SRC,KPP_SRC+task_pos,SIZEDBL*nt*usize*task_size1,(void*)&put_reply,SIZEDBL*(nxy-task_size1),SIZEDBL*task_size1);
   while(put_reply!=1);
   sync_mem;

#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last2_KPP_SRC put with l: %d. \n", l);
#endif
         
   task_pos = task_pos + task_size1;
   }
   // end of for(l)
    



#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("step into last2_HMXL.\n");
#endif
 


//!  compute diagnostic mixed layer depth (cm) using a max buoyancy 
//!  gradient criterion.  Use USTAR and BFSFC as temps.

   task_pos = task_posi;
   for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;


      // get: DBSFC, DZT. 
      get_reply = 0;
      athread_get(PE_MODE,DBSFC+task_pos,s_DBSFC,SIZEDBL*usize*task_size1,(void*)&get_reply,0,SIZEDBL*(nxy-task_size1),SIZEDBL*task_size1);

      athread_get(PE_MODE,DZT+task_pos,s_DZT,SIZEDBL*usize1*task_size1,(void*)&get_reply,0,SIZEDBL*(nxy-task_size1),SIZEDBL*task_size1);
      while(get_reply!=2);
      sync_mem;

#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last2_HMXL get ok with l: %d. \n", l);
#endif
 


      for(j=0;j<task_size1;j++)
        {
        pij = l*task_size + j;
        USTAR = c0;
        s_VISC[task_size1+j] = c0;
        if(s_KMT[pij] == 1)
          s_HMXL[pij] = s_zt[0];
        else
          s_HMXL[pij] = c0;
      
        nk = locmin(km-1,s_KMT[pij]-1);
        for(k=1;k<=nk;k++)
          {
          STABLE = s_zt[k-1] + p5*(s_DZT[(k+1-1)*task_size1+j] + s_DZT[(k+1)*task_size1+j]); 
          USTAR = locmax(s_DBSFC[k*task_size1+j]/STABLE,USTAR);
          s_HMXL[pij] = STABLE;
          }
          //end of for(k) 
 
        for(k=1;k<km;k++)
          {
          if(USTAR > c0)
            s_VISC[(k+1)*task_size1+j] = (s_DBSFC[k*task_size1+j]-s_DBSFC[(k-1)*task_size1+j])/ 
                 (p5*(s_DZT[(k+1)*task_size1+j] + s_DZT[(k+1-1)*task_size1+j]));
          ztmp = s_VISC[(k+1)*task_size1+j]-s_VISC[(k+1-1)*task_size1+j];
          if(s_VISC[(k+1)*task_size1+j] >= USTAR && USTAR > c0 &&
              fabs(ztmp) > c0)  // avoid divide by zero
            {
            BFSFC = (s_VISC[(k+1)*task_size1+j] - USTAR)/ ztmp;
            s_HMXL[pij] = (s_zt[k-1] + p25*(s_DZT[(k+1-1)*task_size1+j] + s_DZT[(k+1)*task_size1+j]))*(c1-BFSFC) 
                 + (s_zt[k-1] - p25*(s_DZT[(k+1-2)*task_size1+j] + s_DZT[(k+1-1)*task_size1+j]))*BFSFC;
            USTAR = c0;
            }
          }
          //end of for(k) 
         
        }
        //end of for(j) 

       

#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last2_HMXL comp with l: %d. \n", l);
#endif
         
      task_pos = task_pos + task_size1;
      }
      // end of for(l)


   //put: HMXL.
   put_reply = 0;
   athread_put(PE_MODE,s_HMXL,HMXL+task_posi,SIZEDBL*task_numj,(void*)&put_reply,0,0);
   while(put_reply!=1);
   sync_mem;


#if defined(POPTESTJN)
 if(myproc==0 && num_calls<8 && myid==0)
    cpe_printf("last2 finished. \n");
#endif

}
// end of if(task_numj>0)


return;
}

