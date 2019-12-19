#include <slave.h>
#include "baroclinic_struct.h"
//#include "simd.h"
//#include "dma.h"
#include <math.h>
#define   s_int sizeof(int)
#define   s_double sizeof(double)
#include "cpe_print.h"
//
#define sync_mem   \
asm volatile("memb")
#define get_myrid(row)  \
asm volatile("rcsr   %0, 1" : "=r"(row))
#define get_mycid(col)  \
asm volatile("rcsr   %0, 2" : "=r"(col))

#define REG_PUTR(var,dest) \
asm volatile ("putr %0,%1\n"::"r"(var),"r"(dest):"memory")
#define REG_GETR(var) \
asm volatile ("getr %0\n":"=r"(var)::"memory")
#define REG_PUTC(var,dest) \
asm volatile ("putc %0,%1\n"::"r"(var),"r"(dest):"memory")
#define REG_GETC(var)  \
asm volatile ("getc %0\n":"=r"(var)::"memory")


#define max(x,y) ( x>y?x:y )
#define min(x,y) ( x<y?x:y )

//################  baroclinic-vert-average-nomal-vec  #########################
void  s_clsec1(struct param_clinic *clsec1)
{
   double *ldm_malloc();
   int volatile get_replya[2],get_reply,put_reply;
   int nx,ny,km,partial,i,k,l1ddyn;
   double c0;
   int myid,thds,current,next,lenmax;
   int nxy,tail,len,ijst,lenk,tijst,tlen,nblocks,block;
   int *kmu;
   double *unew,*vnew,*dzu,*dz,*hur;
   int *s_kmu;
   double *s_unew,*s_vnew,*s_dzu,*s_hur,*s_dz,*s_work1,*s_work2;
   int rank,dsize;

   param_clinic s_sec1;

    get_reply=0;
    athread_get(0,clsec1,&s_sec1,sizeof(s_sec1),(void *)&get_reply,0,0,0);
    while(get_reply!=1);

    nx           = s_sec1.param[0]; 
    ny           = s_sec1.param[1]; 
    km           = s_sec1.param[2];
    partial      = s_sec1.param[3];
    l1ddyn       = s_sec1.param[4];
    thds         = s_sec1.param[5];
    rank         = s_sec1.param[6];

    c0           = s_sec1.dparam[0];

    myid           = athread_get_id(-1);
     
    unew  = s_sec1.addr[0];
    vnew  = s_sec1.addr[1];
    dzu   = s_sec1.addr[2];
    hur   = s_sec1.addr[3];
    dz    = s_sec1.addr[4];

    kmu   = s_sec1.addr_int[0];

       nxy  = nx*ny;
       tail = nxy%thds;
       tlen  = nxy/thds+(tail-myid+thds-1)/thds;
       tijst = nxy/thds*myid+min(myid,tail);

       nblocks=1;
       lenmax = (tlen+nblocks-1)/nblocks;  //for malloc
       lenk=(lenmax*km+3)/4*4;
       lenmax=(lenmax+3)/4*4;
       dsize=lenk*2+lenmax*5+km;
       while(dsize > 5000)
       {
         nblocks++;
         lenmax = (tlen+nblocks-1)/nblocks;  //for malloc
         lenk=(lenmax*km+3)/4*4;
         lenmax=(lenmax+3)/4*4;
         dsize=lenk*2+lenmax*5+km;
       }

     s_unew  = ldm_malloc(s_double*lenk);
     s_vnew  = ldm_malloc(s_double*lenk);
     s_dzu   = ldm_malloc(s_double*lenmax*2);
     s_hur   = ldm_malloc(s_double*lenmax);
     s_work1 = ldm_malloc(s_double*lenmax);
     s_work2 = ldm_malloc(s_double*lenmax);
     s_dz    = ldm_malloc(s_double*km);
     s_kmu  = (int *)ldm_malloc(s_int*lenmax);
     

 for(block=0;block<nblocks;block++)
  { //block <
         tail=tlen%nblocks;
         len=tlen/nblocks+(tail-block+nblocks-1)/nblocks;
         ijst=tlen/nblocks*block+min(block,tail)+tijst;

        get_reply=0;
        athread_get(0,kmu+ijst, s_kmu, len*s_int,   (void *)&get_reply,0,0,0);
    if(l1ddyn !=1)
      {
        get_replya[0]=0;
        athread_get(0,unew+ijst,s_unew,len*s_double,(void *)&get_replya[0],0,0,0);
        athread_get(0,vnew+ijst,s_vnew,len*s_double,(void *)&get_replya[0],0,0,0);
        athread_get(0,hur+ijst, s_hur, len*s_double,(void *)&get_reply,0,0,0);
        if(partial==1) {athread_get(0,dzu+ijst,s_dzu,len*s_double,(void *)&get_replya[0],0,0,0);}
             else      {athread_get(0,dz, s_dz, km*s_double, (void *)&get_reply,0,0,0);}
        for(i=0;i<len;i++)
            {
               s_work1[i]=0;
               s_work2[i]=0;
            }
        if(partial==1)
          {        
            while(get_reply!=2);
            for(k=0;k<km;k++)
            {
              current=k%2;
              next   =(k+1)%2;
              if(k<(km-1)) 
                  {
                    get_replya[next]=0;
                    athread_get(0,unew+(k+1)*nxy+ijst,s_unew+k*len+len,len*s_double,(void *)&get_replya[next],0,0,0);
                    athread_get(0,vnew+k*nxy+nxy+ijst,s_vnew+k*len+len,len*s_double,(void *)&get_replya[next],0,0,0);
                    athread_get(0,dzu+k*nxy+nxy+ijst, s_dzu+next*len,len*s_double,(void *)&get_replya[next],0,0,0);
                  }
              while(get_replya[current]!=3);
              for(i=0;i<len;i++)
                 {
                   s_work1[i]+=s_unew[k*len+i]*s_dzu[current*len+i];
                   s_work2[i]+=s_vnew[k*len+i]*s_dzu[current*len+i];
                 }

            } // kloop right
          } //partial yes right
        else
          {

            while(get_reply!=3);
            for(k=0;k<km;k++)
            {
              current=k%2;
              next   =(k+1)%2;
              if(k<(km-1)) 
                  {
                    get_replya[next]=0;
                    athread_get(0,unew+k*nxy+nxy+ijst,s_unew+k*len+len,len*s_double,(void *)&get_replya[next],0,0,0);
                    athread_get(0,vnew+k*nxy+nxy+ijst,s_vnew+k*len+len,len*s_double,(void *)&get_replya[next],0,0,0);
                  }
              while(get_replya[current]!=2);
              for(i=0;i<len;i++)
                 {
                   s_work1[i]+=s_unew[k*len+i]*s_dz[k];
                   s_work2[i]+=s_vnew[k*len+i]*s_dz[k];
                 }
            } // kloop right
          } // partial no right

          for(i=0;i<len;i++)
             {
               s_work1[i]=s_work1[i]*s_hur[i];
               s_work2[i]=s_work2[i]*s_hur[i];
             }
         
          for(i=0;i<len;i++)
           {
             for(k=0;k<s_kmu[i];k++)
              {
               s_unew[k*len+i]=s_unew[k*len+i]-s_work1[i];
               s_vnew[k*len+i]=s_vnew[k*len+i]-s_work2[i];
              }
             for(k=s_kmu[i];k<km;k++)
              {
               s_unew[k*len+i]=c0;
               s_vnew[k*len+i]=c0;
              }
           }

      } // l1ddyn yes right
    else
      {
        while(get_reply!=1);
        for(i=0;i<len;i++)
           {
             for(k=s_kmu[i];k<km;k++)
              {
               s_unew[k*len+i]=c0;
               s_vnew[k*len+i]=c0;
              }
           }

      } // l1ddyn no right

    //   put
    put_reply=0;
    athread_put(0,s_unew,unew+ijst,len*s_double*km,(void *)&put_reply,(nxy-len)*8,len*8);
    athread_put(0,s_vnew,vnew+ijst,len*s_double*km,(void *)&put_reply,(nxy-len)*8,len*8);
    while(put_reply!=2);

  }  //nblocks right
      ldm_free(s_unew,s_double*lenk);
      ldm_free(s_vnew,s_double*lenk);
      ldm_free(s_dzu,s_double*lenmax*2);
      ldm_free(s_hur,s_double*lenmax);
      ldm_free(s_work1,s_double*lenmax);
      ldm_free(s_work2,s_double*lenmax);
      ldm_free(s_dz,s_double*km);
      ldm_free(s_kmu,s_int*lenmax);
}
/*
      if (.not.l1Ddyn) then
        WORK1 = c0  ! initialize sums
        if (partial_bottom_cells) then
           do k = 1,km
              WORK1 = WORK1 + UVEL(:,:,k,newtime,iblock)*DZU(:,:,k,iblock)
           enddo
        else
           do k = 1,km
              WORK1 = WORK1 + UVEL(:,:,k,newtime,iblock)*dz(k)
           enddo
        endif
        WORK1 = WORK1*HUR(:,:,iblock)  ! normalize by dividing by depth
        do k = 1,km
           where (k <= KMU(:,:,iblock))
              UVEL(:,:,k,newtime,iblock) = &
              UVEL(:,:,k,newtime,iblock) - WORK1
              VVEL(:,:,k,newtime,iblock) = &
              VVEL(:,:,k,newtime,iblock) - WORK2
           elsewhere 
              UVEL(:,:,k,newtime,iblock) = c0
              VVEL(:,:,k,newtime,iblock) = c0
           endwhere
        enddo

      else

        do k = 1,km
          where (k > KMU(:,:,iblock))
            UVEL(:,:,k,newtime,iblock) = c0
            VVEL(:,:,k,newtime,iblock) = c0
          endwhere
        enddo

      endif

   call timer_stop(timer_sec2clinic)
*/


void  s_postadvu(struct param_clinic *postadvu)
{
   double *ldm_malloc();
   int volatile get_reply1,get_reply2,get_reply3,get_reply,put_reply;
   int nx,ny,k,partial,l1ddyn,nxy,k1;
   int imp_leap,lpres_leap,ldiag;
   double c0,c1,p5,p25,dzk,dzwk,dzwkm1,boussk,grav,gamma,c2;
   int  threads,i,j,mycid,myrid,remap_id,pre_cid,next_cid,rank;

   double *diag_ke_adv_2d,*dzu,*ucur,*vcur,*fx,*fy,*fcor,*uold,*vold;
   double *s_diag_ke_adv_2d,*s_dzu,*s_ucur,*s_vcur,*s_fx,*s_fy,*s_fcor;
   double *s_workx,*s_worky,*s_uold,*s_vold,*s_right; //,*s_dkpress2d;
   double *rhoknew,*rhokcur,*rhokold;
   double *s_rhoknew,*s_rhokcur,*s_rhokold;
   double *sumx,*sumy,*rhokmx,*rhokmy;
   double *s_sumx,*s_sumy,*s_rhokmx,*s_rhokmy;
   double *dxur,*dyur,*s_dxur,*s_dyur;
   double *s_rhokx,*s_rhoky,*s_rhoavg;
   double factor,*workx,*s_xtmp,*s_ytmp; //,*dkpress2d;
//   doublev4 s_a[20],s_b[20];

//   static int times=0;
   int *kmu,*s_kmu,myid,ii;
   

   struct param_clinic s_postadu;
//    times++;
//       if(rank==0 && remap_id==0) cpe_printf("astart_%d!\n",times);
//    times--;

/*   athread_syn(ARRAY_SCOPE, 0xffff);
   if (_MYID == 0){
   get_reply = 0;
   athread_get(BCAST_MODE, postadvu, &s_postadu, sizeof(s_postadu), (void *)&get_reply, 0xff, 0, 0);
   while (get_reply != 1);
   sync_mem;
                  }
   athread_syn(ARRAY_SCOPE, 0xffff);
*/
    get_reply=0;
    athread_get(0,postadvu,&s_postadu,sizeof(s_postadu),(void *)&get_reply,0,0,0);
    while(get_reply!=1);
//    sync_mem;


   nx          = s_postadu.param[0];
   ny          = s_postadu.param[1];
   k           = s_postadu.param[2];
   ldiag       = s_postadu.param[3];
   partial     = s_postadu.param[4];
   l1ddyn      = s_postadu.param[5];
   imp_leap    = s_postadu.param[6];
   lpres_leap  = s_postadu.param[7];
   threads     = s_postadu.param[8];   
   rank        = s_postadu.param[9];   

   k1=k-1;
   nxy=nx*ny;
 
   c1          = s_postadu.dparam[0];
   p5          = s_postadu.dparam[1];
   p25         = s_postadu.dparam[2];
   dzk         = s_postadu.dparam[3];
   dzwk        = s_postadu.dparam[4];
   boussk      = s_postadu.dparam[5];
   grav        = s_postadu.dparam[6];
   gamma       = s_postadu.dparam[7];
   c0          = s_postadu.dparam[8];
   dzwkm1      = s_postadu.dparam[9];
   c2          = 2.0;

   get_myrid(myrid);
   get_mycid(mycid);
//   myid=myrid*8+mycid;

/*    times++;
      if(times<300){
      if(rank==0 && myid==0) cpe_printf("nx,k,imp,times %d %d %d %d\n",nx,k,imp_leap,times);
      if(rank==0 && myid==0) cpe_printf("dzwkm1=%lf, dzwk=%lf, grav=%lf, gamma=%lf, dzk=%lf, p25=%lf \n",dzwkm1,dzwk,grav,gamma,dzk,p25);
                   }
*/

   diag_ke_adv_2d  = s_postadu.addr[0];
   dzu             = s_postadu.addr[1];
   ucur            = s_postadu.addr[2];
   vcur            = s_postadu.addr[3];
   uold            = s_postadu.addr[4];
   vold            = s_postadu.addr[5];
   fcor            = s_postadu.addr[6];
   fx              = s_postadu.addr[7];
   fy              = s_postadu.addr[8];
   rhoknew         = s_postadu.addr[9];
   rhokcur         = s_postadu.addr[10];
   rhokold         = s_postadu.addr[11];
   sumx            = s_postadu.addr[12];
   sumy            = s_postadu.addr[13];
   rhokmx          = s_postadu.addr[14];
   rhokmy          = s_postadu.addr[15];
   dxur            = s_postadu.addr[16];
   dyur            = s_postadu.addr[17];
//   dkpress2d       = s_postadu.addr[18];
   workx           = s_postadu.addr[18];

   kmu             = s_postadu.addr_int[0];

   s_diag_ke_adv_2d = ldm_malloc(s_double*nx);
   s_dzu            = ldm_malloc(s_double*nx);
   s_ucur           = ldm_malloc(s_double*nx);
   s_vcur           = ldm_malloc(s_double*nx);
   s_uold           = ldm_malloc(s_double*nx);
   s_vold           = ldm_malloc(s_double*nx);
   s_fcor           = ldm_malloc(s_double*nx);
   s_fx             = ldm_malloc(s_double*nx);
   s_fy             = ldm_malloc(s_double*nx);
   s_rhoknew        = ldm_malloc(s_double*nx*2);
   s_rhokcur        = ldm_malloc(s_double*nx*2);
   s_rhokold        = ldm_malloc(s_double*nx*2);
   s_sumx           = ldm_malloc(s_double*nx);
   s_sumy           = ldm_malloc(s_double*nx);
   s_rhokmx         = ldm_malloc(s_double*nx);
   s_rhokmy         = ldm_malloc(s_double*nx);
   s_dxur           = ldm_malloc(s_double*nx);
   s_dyur           = ldm_malloc(s_double*nx);
   s_workx          = ldm_malloc(s_double*nx);
   s_worky          = ldm_malloc(s_double*nx);
   s_rhokx          = ldm_malloc(s_double*nx);
   s_rhoky          = ldm_malloc(s_double*nx);
   s_rhoavg         = ldm_malloc(s_double*nx*2);
//   s_right          = ldm_malloc(s_double*nx); //get rhoavg(:,j+1)
///   s_dkpress2d      = ldm_malloc(s_double*nx);
   s_xtmp           = ldm_malloc(s_double*nx); 
   s_ytmp           = ldm_malloc(s_double*nx); 
   

   s_kmu            = (int *) ldm_malloc(s_int*nx);

   if(myrid&1)          
     remap_id = myrid*8 + (7-mycid);
   else
     remap_id = myrid*8 + mycid;

//   if(myrid&1) next_cid = mycid-1;         
//   else        next_cid = mycid+1;
   if(myrid&1) pre_cid = mycid+1;         
   else        pre_cid = mycid-1;

   for(j=0;j<ny;j++)
    {
     if(j%threads == remap_id) 
      {

//        times++;
//       if(rank==0 && remap_id==0) cpe_printf("cstart_%d!\n",times);
//        times--;
        get_reply=0; get_reply1=0;
        get_reply2=0;get_reply3=0;
        athread_get(0,diag_ke_adv_2d+j*nx,s_diag_ke_adv_2d,s_double*nx,(void *)&get_reply,0,0,0);
        athread_get(0,dzu +k1*nxy+j*nx,   s_dzu,           s_double*nx,(void *)&get_reply,0,0,0);
        athread_get(0,ucur+k1*nxy+j*nx,   s_ucur,          s_double*nx,(void *)&get_reply,0,0,0);
        athread_get(0,vcur+k1*nxy+j*nx,   s_vcur,          s_double*nx,(void *)&get_reply,0,0,0);
        athread_get(0,fx+         j*nx,   s_fx,            s_double*nx,(void *)&get_reply,0,0,0);
        athread_get(0,fy+         j*nx,   s_fy,            s_double*nx,(void *)&get_reply,0,0,0);
        athread_get(0,fcor+       j*nx,   s_fcor,          s_double*nx,(void *)&get_reply1,0,0,0);
        athread_get(0,uold+k1*nxy+j*nx,   s_uold,          s_double*nx,(void *)&get_reply1,0,0,0);
        athread_get(0,vold+k1*nxy+j*nx,   s_vold,          s_double*nx,(void *)&get_reply1,0,0,0);
       if(l1ddyn == 0 && j<(ny-1))
        {
        athread_get(0,rhoknew+     j*nx,  s_rhoknew,     2*s_double*nx,(void *)&get_reply2,0,0,0);
        athread_get(0,rhokcur+     j*nx,  s_rhokcur,     2*s_double*nx,(void *)&get_reply2,0,0,0);
        athread_get(0,rhokold+     j*nx,  s_rhokold,     2*s_double*nx,(void *)&get_reply2,0,0,0);
        athread_get(0,dxur+        j*nx,  s_dxur,          s_double*nx,(void *)&get_reply2,0,0,0);
        athread_get(0,dyur+        j*nx,  s_dyur,          s_double*nx,(void *)&get_reply2,0,0,0);
        athread_get(0,kmu+         j*nx,  s_kmu,           s_int*nx,   (void *)&get_reply2,0,0,0);
        }
       if(l1ddyn ==0 && k>1)
        {
        athread_get(0,sumx+        j*nx,  s_sumx,          s_double*nx,(void *)&get_reply3,0,0,0);
        athread_get(0,sumy+        j*nx,  s_sumy,          s_double*nx,(void *)&get_reply3,0,0,0);
        athread_get(0,rhokmx+      j*nx,  s_rhokmx,        s_double*nx,(void *)&get_reply3,0,0,0);
        athread_get(0,rhokmy+      j*nx,  s_rhokmy,        s_double*nx,(void *)&get_reply3,0,0,0);
        }
//        athread_get(0,dkpress2d+   j*nx,  s_dkpress2d,     s_double*nx,(void *)&get_reply3,0,0,0);

        while(get_reply!=6);
        if(ldiag == 1) 
          {
            if(partial==1) 
              {
                for(i=0;i<nx;i++)
                  {
		    s_diag_ke_adv_2d[i]+=s_dzu[i]*(s_ucur[i]*s_fx[i]+s_vcur[i]*s_fy[i]);
                  }
              }
            else
              {
                for(i=0;i<nx;i++)
                  {
		    s_diag_ke_adv_2d[i]+=dzk*(s_ucur[i]*s_fx[i]+s_vcur[i]*s_fy[i]);
                  }
              }
          }  //if-ldiag

       put_reply=0;
       athread_put(0,s_diag_ke_adv_2d,diag_ke_adv_2d+j*nx,s_double*nx,(void *)&put_reply,0,0);

       while(get_reply1!=3);       

       if(imp_leap == 1)
         {
          for(i=0;i<nx;i++)
           {
             s_fx[i]+=s_fcor[i]*(gamma*s_vcur[i]+(c1-gamma)*s_vold[i]); 
             s_fy[i]-=s_fcor[i]*(gamma*s_ucur[i]+(c1-gamma)*s_uold[i]); 
           }
         }
       if(imp_leap == 2)
         {
          for(i=0;i<nx;i++)
           {
             s_fx[i]+=s_fcor[i]*s_vcur[i]; 
             s_fy[i]-=s_fcor[i]*s_ucur[i]; 
           }
         }
       if(imp_leap == 3)
         {
          for(i=0;i<nx;i++)
           {
             s_fx[i]+=s_fcor[i]*s_vold[i];
             s_fy[i]-=s_fcor[i]*s_uold[i];
           }
         }

       if(l1ddyn == 0) 
         {
            for(i=0;i<nx;i++)
             {
               s_rhokx[i]=c0;
               s_rhoky[i]=c0;
             }

           if(j<ny-1)
           {
             while(get_reply2!=6);
             if(lpres_leap == 1)
               {
                 for(i=0;i<(nx+nx);i++) s_rhoavg[i]=p25*(s_rhoknew[i]+c2*s_rhokcur[i]+s_rhokold[i])*boussk;
               }
             else
               {
                 for(i=0;i<(nx+nx);i++) s_rhoavg[i]=s_rhokcur[i]*boussk;
               }

             for(i=0;i<nx;i++)
              {
                s_xtmp[i]=s_rhoavg[nx+i]+s_rhoavg[i];
                s_ytmp[i]=s_rhoavg[nx+i]-s_rhoavg[i];
              }

             for(i=0;i<nx-1;i++)
               {
                 if(k<=s_kmu[i])
                   {
                     s_rhokx[i]=s_dxur[i]*p5*(s_xtmp[i+1]-s_xtmp[i]);
                     s_rhoky[i]=s_dyur[i]*p5*(s_ytmp[i+1]+s_ytmp[i]);
                   }
               } 
           } //j<ny-1 right

             if(k==1)
              {
                for(i=0;i<nx;i++)
                 {
                   s_rhokmx[i]=s_rhokx[i];
                   s_rhokmy[i]=s_rhoky[i];
                   s_sumx[i]=c0;
                   s_sumy[i]=c0;
                 }
              }

          factor = dzwkm1*grav*p5;
        if(k>1) while(get_reply3!=4);
          for(i=0;i<nx;i++)
           {
             s_sumx[i]+=factor*(s_rhokx[i]+s_rhokmx[i]);
             s_sumy[i]+=factor*(s_rhoky[i]+s_rhokmy[i]);
             s_workx[i]=s_sumx[i];
             s_worky[i]=s_sumy[i];
           }

          for(i=0;i<nx;i++)
           {
             s_rhokmx[i]=s_rhokx[i];
             s_rhokmy[i]=s_rhoky[i];
             s_fx[i]-=s_workx[i];
             s_fy[i]-=s_worky[i];
           
           }
           athread_put(0,s_sumx,sumx+j*nx,s_double*nx,(void*)&put_reply,0,0);
           athread_put(0,s_sumy,sumy+j*nx,s_double*nx,(void*)&put_reply,0,0);
           athread_put(0,s_rhokmx,rhokmx+j*nx,s_double*nx,(void*)&put_reply,0,0);
           athread_put(0,s_rhokmy,rhokmy+j*nx,s_double*nx,(void*)&put_reply,0,0);

         } //l1ddyn if
       else
         {
           for(i=0;i<nx;i++)
            {
              s_workx[i]=c0;
              s_worky[i]=c0;
            }


         }

       athread_put(0,s_fx,fx+j*nx,s_double*nx,(void*)&put_reply,0,0);
       athread_put(0,s_fy,fy+j*nx,s_double*nx,(void*)&put_reply,0,0);
       if(partial==1)
        {

            s_workx[i]=-s_dzu[i]*(s_ucur[i]*s_workx[i]+s_vcur[i]*s_worky[i]);
            
        }
       else
        {
          for(i=0;i<nx;i++)
            {
            s_workx[i]=-dzk*(s_ucur[i]*s_workx[i]+s_vcur[i]*s_worky[i]);
            }
        }

        athread_put(0,s_workx,workx+j*nx,s_double*nx,(void*)&put_reply,0,0);

//	if(rank==0 && remap_id==0) cpe_printf("putreply_bf %d\n",put_reply);
        if(l1ddyn==0) while(put_reply!=8);
        else  while(put_reply!=4);
//	if(rank==0 && remap_id==0) cpe_printf("putreply_af %d\n",put_reply);
          
/*   if (partial_bottom_cells) then
      WORKX =  -DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                VCUR(:,:,k)*WORKY)
   else
      WORKX =  -dz(k)*(UCUR(:,:,k)*WORKX + &
                       VCUR(:,:,k)*WORKY)
   endif

   call accumulate_tavg_field(WORKX,tavg_UDP,bid,k)

   if (ldiag_global) then
      DIAG_KE_PRESS_2D(:,:,bid) = DIAG_KE_PRESS_2D(:,:,bid) + WORKX
   endif
*/
      



      } //if mod-j
    } //j
 
   ldm_free(s_diag_ke_adv_2d,s_double*nx);
   ldm_free(s_dzu           ,s_double*nx);
   ldm_free(s_ucur          ,s_double*nx);
   ldm_free(s_vcur          ,s_double*nx);
   ldm_free(s_fcor          ,s_double*nx);
   ldm_free(s_fx            ,s_double*nx);
   ldm_free(s_fy            ,s_double*nx);
   ldm_free(s_workx         ,s_double*nx);
   ldm_free(s_worky         ,s_double*nx);
   ldm_free(s_uold          ,s_double*nx);
   ldm_free(s_vold          ,s_double*nx);
   ldm_free(s_rhoknew       ,s_double*nx*2);
   ldm_free(s_rhokcur       ,s_double*nx*2);
   ldm_free(s_rhokold       ,s_double*nx*2);
   ldm_free(s_rhoavg        ,s_double*nx*2);
   ldm_free(s_sumx          ,s_double*nx);
   ldm_free(s_sumy          ,s_double*nx);
   ldm_free(s_rhokmx        ,s_double*nx);
   ldm_free(s_rhokmy        ,s_double*nx);
   ldm_free(s_dxur          ,s_double*nx);
   ldm_free(s_dyur          ,s_double*nx);
//   ldm_free(s_right         ,s_double*nx);
   ldm_free(s_rhokx         ,s_double*nx);
   ldm_free(s_rhoky         ,s_double*nx);
   ldm_free(s_xtmp          ,s_double*nx);
   ldm_free(s_ytmp          ,s_double*nx);
   ldm_free(s_kmu           ,s_int*nx   );

//  if(rank==0 && remap_id==0)
//           { cpe_printf("all_finish!\n");}


}


void  s_grad(struct param_clinic *grad)
{
   double *ldm_malloc();
   int volatile get_reply1,get_reply2,get_reply3,get_reply,put_reply;
   int i,j,myrid,mycid,remap_id,pre_cid;

   double c0,p5;
   int    nx,ny,k,threads;
   double *rhokx,*rhoky,*rhoavg,*dxur,*dyur;
   double *s_rhokx,*s_rhoky,*s_rhoavg,*s_dxur,*s_dyur,*s_right;
   int    *kmu,*s_kmu;

   struct param_clinic s_grad;

   get_reply=0;
   athread_get(0,grad,&s_grad,sizeof(s_grad),(void *)&get_reply,0,0,0);
   while(get_reply!=1);

   nx=s_grad.param[0];
   ny=s_grad.param[1];
   k =s_grad.param[2];
   threads=s_grad.param[3];
   c0=s_grad.dparam[0];
   p5=s_grad.dparam[1];
   rhokx=s_grad.addr[0];
   rhoky=s_grad.addr[1];
   rhoavg=s_grad.addr[2];
   dxur=s_grad.addr[3];
   dyur=s_grad.addr[4];
   kmu=s_grad.addr_int[0];

   s_rhokx  = ldm_malloc(nx*s_double);
   s_rhoky  = ldm_malloc(nx*s_double);
   s_rhoavg = ldm_malloc(nx*s_double);
   s_dxur   = ldm_malloc(nx*s_double);
   s_dyur   = ldm_malloc(nx*s_double);
   s_right  = ldm_malloc(nx*s_double);
   s_kmu    = (int *) ldm_malloc(nx*s_int);

   get_myrid(myrid);
   get_mycid(mycid);
   if(myrid&1)
     remap_id = myrid*8 + (7-mycid);
   else
     remap_id = myrid*8 + mycid;
   if(myrid&1) pre_cid = mycid+1;
   else        pre_cid = mycid-1;

   for(j=0;j<ny-1;j++)
   {
     if(j%threads == remap_id)
     {
       get_reply=0;
       athread_get(0,rhoavg+j*nx,s_rhoavg,nx*s_double,(void *)&get_reply,0,0,0);
       athread_get(0,dxur+j*nx,  s_dxur,  nx*s_double,(void *)&get_reply,0,0,0);
       athread_get(0,dyur+j*nx,  s_dyur,  nx*s_double,(void *)&get_reply,0,0,0);
       athread_get(0,kmu+j*nx,   s_kmu,   nx*s_int,   (void *)&get_reply,0,0,0);

       while(get_reply!=4);
       sync_mem;

             if(remap_id%8 >0 )
              {
               for(i=0;i<nx;i++)
                {
                  REG_PUTR(s_rhoavg[i],pre_cid);
                  sync_mem;
                }
              }
             else
              {
               if(remap_id>0)
                {
                  for(i=0;i<nx;i++)
                   {
                    REG_PUTC(s_rhoavg[i],myrid-1);
                    sync_mem;
                   }
                }
              }

             if(remap_id%8 <7 && remap_id <(ny-1))
               {

                 REG_GETR(s_right[0]);
                 for(i=0;i<nx-1;i++)
                   {
                     REG_GETR(s_right[i+1]);
                     sync_mem;
                     if(k<=s_kmu[i])
                     {
                       s_rhokx[i]=s_dxur[i]*p5*(s_right[i+1]-s_rhoavg[i]-s_right[i]+s_rhoavg[i+1]);
                       s_rhoky[i]=s_dyur[i]*p5*(s_right[i+1]-s_rhoavg[i]+s_right[i]-s_rhoavg[i+1]);
                     }
                   }
               }
             else
               {
                 if(remap_id<(ny-1))
                   {
                     REG_GETC(s_right[0]);
                     for(i=0;i<nx-1;i++)
                      {
                        REG_GETC(s_right[i+1]);
                        sync_mem;
                        if(k<=s_kmu[i])
                        {
                          s_rhokx[i]=s_dxur[i]*p5*(s_right[i+1]-s_rhoavg[i]-s_right[i]+s_rhoavg[i+1]);
                          s_rhoky[i]=s_dyur[i]*p5*(s_right[i+1]-s_rhoavg[i]+s_right[i]-s_rhoavg[i+1]);
                        }

                      }
                   }
               }
       
        put_reply=0;
        athread_put(0,s_rhokx,rhokx+j*ny,nx*s_double,(void *)&put_reply,0,0);
        athread_put(0,s_rhoky,rhoky+j*ny,nx*s_double,(void *)&put_reply,0,0);

        while(put_reply!=2);


     }
   }


   ldm_free(s_rhokx,nx*s_double);
   ldm_free(s_rhoky,nx*s_double);
   ldm_free(s_rhoavg,nx*s_double);
   ldm_free(s_dxur,nx*s_double);
   ldm_free(s_dyur,nx*s_double);
   ldm_free(s_right,nx*s_double);
   ldm_free(s_kmu,nx*s_int);


}


void  s_storeforce(struct little_clinic *stforce)
{
   double *ldm_malloc();
   int volatile get_reply,put_reply,get_reply1;
   int nx,ny,k,i,thds,myid; 
   int impcor,l1ddyn,partial;
   int tail,len,ijst,nxy;
   double c2dtu,beta,c1,dzk;

   double *fcor,*fx,*fy,*zx,*zy,*uvel,*vvel,*dzu;
   double *s_fcor,*s_fx,*s_fy,*s_zx,*s_zy,*s_uvel,*s_vvel,*s_dzu;
   double *s_work1,*s_work2;
   
   struct little_clinic s_stforce;
    get_reply=0;
    athread_get(0,stforce,&s_stforce,sizeof(s_stforce),(void *)&get_reply,0,0,0);
    while(get_reply!=1);
    
    nx     =s_stforce.param[0];
    ny     =s_stforce.param[1];
    k      =s_stforce.param[2];
    impcor =s_stforce.param[3];
    l1ddyn =s_stforce.param[4];
    partial=s_stforce.param[5];
    thds   =s_stforce.param[6];
    myid   =athread_get_id(-1);
   
    c2dtu  =s_stforce.dparam[0];
    beta   =s_stforce.dparam[1];
    c1     =s_stforce.dparam[2];
    dzk    =s_stforce.dparam[3];

    fcor    =s_stforce.addr[0];
    fx      =s_stforce.addr[1];
    fy      =s_stforce.addr[2];
    uvel    =s_stforce.addr[3];
    vvel    =s_stforce.addr[4];
    zx      =s_stforce.addr[5];
    zy      =s_stforce.addr[6];
    dzu     =s_stforce.addr[7];

       nxy  = nx*ny;
       tail = nxy%thds;
       len  = nxy/thds+(tail-myid+thds-1)/thds;
       ijst = nxy/thds*myid+min(myid,tail);

    s_fcor = ldm_malloc(len*s_double);
    s_fx   = ldm_malloc(len*s_double);
    s_fy   = ldm_malloc(len*s_double);
    s_zx   = ldm_malloc(len*s_double);
    s_zy   = ldm_malloc(len*s_double);
    s_uvel = ldm_malloc(len*s_double);
    s_vvel = ldm_malloc(len*s_double);
    s_dzu  = ldm_malloc(len*s_double);
    s_work1= ldm_malloc(len*s_double);
    s_work2= ldm_malloc(len*s_double);

    get_reply=0;get_reply1=0;
    if(impcor==1) athread_get(0,fcor+ijst,s_fcor,len*s_double,(void *)&get_reply,0,0,0);
    athread_get(0,fx  +ijst,s_fx  ,len*s_double,(void *)&get_reply,0,0,0);
    athread_get(0,fy  +ijst,s_fy  ,len*s_double,(void *)&get_reply,0,0,0);
    if(l1ddyn==0) {
    if(partial==1) athread_get(0,dzu+k*nxy+ijst,s_dzu,len*s_double,(void *)&get_reply1,0,0,0);
    athread_get(0,zx  +ijst,     s_zx ,len*s_double,(void *)&get_reply1,0,0,0);
    athread_get(0,zy  +ijst,     s_zy ,len*s_double,(void *)&get_reply1,0,0,0);
                  } 

    if(impcor==1)
      {
        while(get_reply!=3);
        for(i=0;i<len;i++)
          {
            s_work1[i]=c2dtu*beta*s_fcor[i];
            s_work2[i]=c2dtu/(c1+s_work1[i]*s_work1[i]);
            s_uvel[i]=(s_fx[i]+s_work1[i]*s_fy[i])*s_work2[i];
            s_vvel[i]=(s_fy[i]-s_work1[i]*s_fx[i])*s_work2[i];
          }
      }
    else
      {
        while(get_reply!=2);
        for(i=0;i<len;i++)
          {
            s_uvel[i]=c2dtu*s_fx[i];
            s_vvel[i]=c2dtu*s_fy[i];
          }
      }  //impcor else >
     put_reply=0;
     athread_put(0,s_uvel,uvel+(k-1)*nxy+ijst,len*s_double,(void *)&put_reply,0,0);
     athread_put(0,s_vvel,vvel+(k-1)*nxy+ijst,len*s_double,(void *)&put_reply,0,0);

     if(l1ddyn == 0)
       {
         if(partial==1) 
           {
             while(get_reply1!=3);
             for(i=0;i<len;i++)
               {              
                 s_zx[i]+=s_fx[i]*s_dzu[i];
                 s_zy[i]+=s_fy[i]*s_dzu[i];
               }
           } //partial >
         else
           {
             for(i=0;i<len;i++)
               {
                  while(get_reply1!=2);
                  s_zx[i]+=s_fx[i]*dzk;
                  s_zy[i]+=s_fy[i]*dzk;
               }

           }//partial else>
         athread_put(0,s_zx,zx+ijst,len*s_double,(void *)&put_reply,0,0);
         athread_put(0,s_zy,zy+ijst,len*s_double,(void *)&put_reply,0,0);
         while(put_reply!=4);
       } //l1dynn >
     else
       {
          while(put_reply!=2);
       }

    ldm_free(s_fcor,len*s_double);
    ldm_free(s_fx,len*s_double);
    ldm_free(s_fy,len*s_double);
    ldm_free(s_zx,len*s_double);
    ldm_free(s_zy,len*s_double);
    ldm_free(s_uvel,len*s_double);
    ldm_free(s_vvel,len*s_double);
    ldm_free(s_dzu,len*s_double);
    ldm_free(s_work1,len*s_double);
    ldm_free(s_work2,len*s_double);
}
