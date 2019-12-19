
#include <stdio.h>
#include <slave.h>
//#include "simd.h"
//#include "dma.h"
#include <math.h>
#define   s_int sizeof(int)
#define   s_double sizeof(double)
#include "vertical_mix_struct.h"

#define max(x,y) ( x>y?x:y )
#define min(x,y) ( x<y?x:y )

//################  vdiffu  #########################
void  s_vmix(struct param_vmix *v1)
{
    double s_dw[2048];
    int    s_iw[128];

    int volatile get_reply,get_reply0,get_reply1,put_reply;

    int nx,ny,k,km,ist,ied,jst,jed,partial,myid; //input
    int thds_p,remap_id;
    double p5,drag,dzwrk,dzrk,c0; //input
    double *uold,*vold,*vvc,*dzu,*smf,*vuf,*vvf,*vduk,*vdvk; //input
    int *kmu; //input

    int nxy,tail,len,le4,le24,ijst,i,io,jo; //index vars
    double *s_uold,*s_vold,*s_vvc,*s_dzu,*s_smf,*s_vuf,*s_vvf,*s_vduk,*s_vdvk; //slave
    double *s_vufb,*s_vvfb; //slave
    int *s_kmu; //slave

    struct param_vmix s_v2;
    double vmag,dtmp; //local
    double one=1.0,zero=0.0;


    get_reply=0;
    athread_get(0,v1,&s_v2,sizeof(s_v2),(void *)&get_reply,0,0,0);
    while(get_reply!=1);


    nx           = s_v2.param[0]; 
    ny           = s_v2.param[1]; 
    k            = s_v2.param[2];
    km           = s_v2.param[3];
    ist          = s_v2.param[4]; 
    ied          = s_v2.param[5]; 
    jst          = s_v2.param[6];
    jed          = s_v2.param[7];
    partial      = s_v2.param[8];
    thds_p       = s_v2.param[9]/2;   //0-31u, 32-63v

    p5           = s_v2.dparam[0];
    drag         = s_v2.dparam[1];
    dzwrk        = s_v2.dparam[2];
    dzrk         = s_v2.dparam[3];
    c0           = s_v2.dparam[4];

    myid           = athread_get_id(-1);
    remap_id       =myid%thds_p;
     
    uold = s_v2.addr[0];
    vold = s_v2.addr[1];
    vvc  = s_v2.addr[2];
    dzu  = s_v2.addr[3];
    smf  = s_v2.addr[4];
    vuf  = s_v2.addr[5];
    vvf  = s_v2.addr[6];
    vduk = s_v2.addr[7];
    vdvk = s_v2.addr[8];
 
    kmu  = s_v2.addr_int[0];

       nxy  = nx*ny;
//       tail = nxy%thds;
       tail = nxy%thds_p;
//       len  = nxy/thds+(tail-myid+thds-1)/thds;
       len  = nxy/thds_p+(tail-remap_id+thds_p-1)/thds_p;
//       ijst = nxy/thds*myid+min(myid,tail);
       ijst = nxy/thds_p*remap_id+min(remap_id,tail);
       le4  = (len+3)/4*4;
       le24 = (len+len+3)/4*4;

     s_uold = s_dw;
     s_vold = s_dw+le24;
     s_dzu  = s_vold+le24; 
     s_smf  = s_dzu +le24; // 2*len dzu
//     s_vvc  = s_smf +le24;
     s_vvc  = s_smf +le4; //len for u, for v
     s_vufb = s_vvc +le4;
     s_vvfb = s_vvc +le4;
//     s_vvfb = s_vufb+le4;
     s_vduk = s_vvfb+le4;
     s_vdvk = s_vvfb+le4;
//     s_vdvk = s_vduk+le4;
     s_vuf  = s_vdvk+le4;
     s_vvf  = s_vdvk+le4;
//     s_vvf =  s_vuf +le4;

//     s_kmu  = (int*)ldm_malloc(s_int*len);
     s_kmu  = s_iw;

       get_reply =0;
       get_reply0=0;
       get_reply1=0;
       put_reply=0;

   if(myid < thds_p)
     {
       athread_get(0,uold+ijst,s_uold,s_double*len*2,(void *)&get_reply,0,(nxy-len)*8,len*8);
       athread_get(0,vold+ijst,s_vold,s_double*len*2,(void *)&get_reply,0,(nxy-len)*8,len*8);
       athread_get(0, dzu+ijst,s_dzu, s_double*len*2,(void *)&get_reply,0,(nxy-len)*8,len*8);
       athread_get(0, vvc+ijst,s_vvc, s_double*len,  (void *)&get_reply,0,0,0);
       if(k==1)
         {
//         athread_get(0,smf+ijst,s_smf,s_double*len*2,(void *)&get_reply1,0,(nxy-len)*8,len*8);
           athread_get(0,smf+ijst,s_smf,s_double*len,(void *)&get_reply1,0,0,0);
//         only the first half
         }
       else
         {
           athread_get(0,vuf+ijst,s_vuf,s_double*len,(void *)&get_reply1,0,0,0);
//         athread_get(0,vvf+ijst,s_vvf,s_double*len,(void *)&get_reply1,0,0,0);
         }
       athread_get(0,kmu+ijst,s_kmu,s_int*len,(void *)&get_reply0,0,0,0);
       while(get_reply!=4);

// ifblock2
       if(k<km) 
         { //if k-block2
           if(partial==1) 
             { //if partial-block
               for (i=0;i<len;i++)
                 { // i loop left
//                 *(s_work+i)=p5*( *(s_dzu+i) + *(s_dzu+i+len));
                   dtmp=one/(p5*( *(s_dzu+i) + *(s_dzu+i+len)));
                   //s_vufb[i]=s_vvc[i]*(s_uold[i]-s_uold[i+len])/s_work[i];
                   *(s_vufb+i)= *(s_vvc+i) * ( *(s_uold+i) - *(s_uold+i+len))*dtmp;
                   //s_vvfb[i]=s_vvc[i]*(s_vold[i]-s_vold[i+len])/s_work[i];
//                 *(s_vvfb+i)= *(s_vvc+i)*( *(s_vold+i)- *(s_vold+i+len))*dtmp;
                 } //i loop right
             } //if partial right
          else
             { //else pa left
               for (i=0;i<len;i++)
                 { //i left
                    *(s_vufb+i)= *(s_vvc+i)* ( *(s_uold+i)- *(s_uold+i+len))*dzwrk;
 //                 *(s_vvfb+i)= *(s_vvc+i)* ( *(s_vold+i)- *(s_vold+i+len))*dzwrk;
                 } //i right
             } // else par right 
         } //if k block2 right                        
       else 
         { // else k left
            for (i=0;i<len;i++)
              {
              *(s_vufb+i)=zero;
//            *(s_vvfb+i)=0.0;
              }
         }  //if k-block2 else right

//   do j=this_block%jb,this_block%je
//   do i=this_block%ib,this_block%ie
//      if (k == KMU(i,j,bid)) then
//         vmag = bottom_drag*sqrt(UOLD(i,j,k)**2 + VOLD(i,j,k)**2)
//         VUFB(i,j) = vmag*UOLD(i,j,k)
//         VVFB(i,j) = vmag*VOLD(i,j,k)
//      endif
//   end do
//   end do

    
    while(get_reply0!=1);
    for(i=0;i<len;i++)
      {
        io=(i+ijst)%nx;
        jo=(i+ijst)/nx;
        if(k==s_kmu[i] && io>=(ist-1) && io<ied && jo>=(jst-1) && jo<jed) 
          {
            vmag=drag*sqrt( *(s_uold+i) * *(s_uold+i)+ *(s_vold+i)* *(s_vold+i));
            *(s_vufb+i) = vmag * *(s_uold+i); 
//          *(s_vvfb+i) = vmag * *(s_vold+i); 
          }
      }//i> 


//    if (partial_bottom_cells) then
//      VDUK = merge((VUF(:,:,bid) - VUFB)/DZU(:,:,k,bid), &
//                   c0, k <= KMU(:,:,bid))
//      VDVK = merge((VVF(:,:,bid) - VVFB)/DZU(:,:,k,bid), &
//                   c0, k <= KMU(:,:,bid))
//    else
//      VDUK = merge((VUF(:,:,bid) - VUFB)*dzr(k), &
//                        c0, k <= KMU(:,:,bid))
//      VDVK = merge((VVF(:,:,bid) - VVFB)*dzr(k), &
//                   c0, k <= KMU(:,:,bid))
//    endif
//
//
//   if (k == 1) then
//         VUF(:,:,bid) = merge(SMF(:,:,1), c0, KMU(:,:,bid) >= 1)
//         VVF(:,:,bid) = merge(SMF(:,:,2), c0, KMU(:,:,bid) >= 1)
//   endif
//
     if(k==1)
       {
          while(get_reply1!=1);
          for(i=0;i<len;i++)
            {
              if(s_kmu[i] >=1)
                {
                  s_vuf[i]=s_smf[i];
//                s_vvf[i]=s_smf[i+len];
                }
            }
       } //k==1>
     else
       {
          while(get_reply1!=1);
       } 

     if(partial==1) 
       {
          for(i=0;i<len;i++)
            {
               if(k<=s_kmu[i])
                 {
                   *(s_vduk+i)=(*(s_vuf+i)- *(s_vufb+i))/ *(s_dzu+i);
                   // *(s_vdvk+i)=(*(s_vvf+i)- *(s_vvfb+i))/ *(s_dzu+i);
                 }
               else
                 {
                   *(s_vduk+i)=c0;//*(s_vdvk+i)=c0;
                 }
            }  // i >
//        for(i=0;i<len;i++){
//           if(k>s_kmu[i]) {*(s_vduk+i)=c0;//*(s_vdvk+i)=c0;
//                          }
//                          }
       } // partial >
     else
       {
           for(i=0;i<len;i++)
             {
               if(k<=s_kmu[i])
                 {
                   *(s_vduk+i)=(*(s_vuf+i)- *(s_vufb+i))*dzrk;
                   // *(s_vdvk+i)=(*(s_vvf+i)- *(s_vvfb+i))*dzrk;
                 }
               else
                 {
                   *(s_vduk+i)=c0;//*(s_vdvk+i)=c0;
                 }
             } //i >
//                      for(i=0;i<len;i++){
//                         if(k>s_kmu[i]) {*(s_vduk+i)=c0;//*(s_vdvk+i)=c0;
//                                        }
//                                        }
       } //partial no >

      athread_put(0,s_vduk,vduk+ijst,s_double*len,(void *)&put_reply,0,0);
//    athread_put(0,s_vdvk,vdvk+ijst,s_double*len,(void *)&put_reply,0,0);
      athread_put(0,s_vufb, vuf+ijst,s_double*len,(void *)&put_reply,0,0);
//    athread_put(0,s_vvfb, vvf+ijst,s_double*len,(void *)&put_reply,0,0);
//    while(put_reply!=4); 
      while(put_reply!=2); 
    } //first half >
  else
    { // second half <

       athread_get(0,uold+ijst,s_uold,s_double*len*2,(void *)&get_reply,0,(nxy-len)*8,len*8);
       athread_get(0,vold+ijst,s_vold,s_double*len*2,(void *)&get_reply,0,(nxy-len)*8,len*8);
       athread_get(0, dzu+ijst,s_dzu, s_double*len*2,(void *)&get_reply,0,(nxy-len)*8,len*8);
       athread_get(0, vvc+ijst,s_vvc, s_double*len,  (void *)&get_reply,0,0,0);
       if(k==1)
         {
           athread_get(0,smf+ijst+nxy,s_smf,s_double*len,(void *)&get_reply1,0,0,0);
         }
       else
         {
//           athread_get(0,vuf+ijst,s_vuf,s_double*len,(void *)&get_reply1,0,0,0);
           athread_get(0,vvf+ijst,s_vvf,s_double*len,(void *)&get_reply1,0,0,0);
         }
       athread_get(0,kmu+ijst,s_kmu,s_int*len,(void *)&get_reply0,0,0,0);
       while(get_reply!=4);

// ifblock2
       if(k<km) 
         { //if k-block2
           if(partial==1) 
             { //if partial-block
               for (i=0;i<len;i++)
                 { // i loop left
//                 *(s_work+i)=p5*( *(s_dzu+i) + *(s_dzu+i+len));
                   dtmp=one/(p5*( *(s_dzu+i) + *(s_dzu+i+len)));
                   //s_vufb[i]=s_vvc[i]*(s_uold[i]-s_uold[i+len])/s_work[i];
//                 *(s_vufb+i)= *(s_vvc+i) * ( *(s_uold+i) - *(s_uold+i+len))*dtmp;
                   //s_vvfb[i]=s_vvc[i]*(s_vold[i]-s_vold[i+len])/s_work[i];
                   *(s_vvfb+i)= *(s_vvc+i)*( *(s_vold+i)- *(s_vold+i+len))*dtmp;
                  } //i loop right
              } //if partial right
          else   
              { //else pa left
                for (i=0;i<len;i++)
                  { //i left
 //                 *(s_vufb+i)= *(s_vvc+i)* ( *(s_uold+i)- *(s_uold+i+len))*dzwrk;
                    *(s_vvfb+i)= *(s_vvc+i)* ( *(s_vold+i)- *(s_vold+i+len))*dzwrk;
                  } //i right
              } // else par right 
         } //if k block2 right                        
       else 
         { // else k left
            for (i=0;i<len;i++)
              {
//              *(s_vufb+i)=0.0;
                *(s_vvfb+i)=zero;
              }
         }  //if k-block2 else right

//   do j=this_block%jb,this_block%je
//   do i=this_block%ib,this_block%ie
//      if (k == KMU(i,j,bid)) then
//         vmag = bottom_drag*sqrt(UOLD(i,j,k)**2 + VOLD(i,j,k)**2)
//         VUFB(i,j) = vmag*UOLD(i,j,k)
//         VVFB(i,j) = vmag*VOLD(i,j,k)
//      endif
//   end do
//   end do

    
    while(get_reply0!=1);
    for(i=0;i<len;i++)
      {
        io=(i+ijst)%nx;
        jo=(i+ijst)/nx;
        if(k==s_kmu[i] && io>=(ist-1) && io<ied && jo>=(jst-1) && jo<jed) 
          {
             vmag=drag*sqrt( *(s_uold+i) * *(s_uold+i)+ *(s_vold+i)* *(s_vold+i));
 //         *(s_vufb+i) = vmag * *(s_uold+i); 
            *(s_vvfb+i) = vmag * *(s_vold+i); 
          } //if >
      }  //i >


//    if (partial_bottom_cells) then
//      VDUK = merge((VUF(:,:,bid) - VUFB)/DZU(:,:,k,bid), &
//                   c0, k <= KMU(:,:,bid))
//      VDVK = merge((VVF(:,:,bid) - VVFB)/DZU(:,:,k,bid), &
//                   c0, k <= KMU(:,:,bid))
//    else
//      VDUK = merge((VUF(:,:,bid) - VUFB)*dzr(k), &
//                        c0, k <= KMU(:,:,bid))
//      VDVK = merge((VVF(:,:,bid) - VVFB)*dzr(k), &
//                   c0, k <= KMU(:,:,bid))
//    endif
//
//
//   if (k == 1) then
//         VUF(:,:,bid) = merge(SMF(:,:,1), c0, KMU(:,:,bid) >= 1)
//         VVF(:,:,bid) = merge(SMF(:,:,2), c0, KMU(:,:,bid) >= 1)
//   endif
//
     if(k==1)
       {
          while(get_reply1!=1);
          for(i=0;i<len;i++)
            {
              if(s_kmu[i] >=1)
                {
//                s_vuf[i]=s_smf[i];
//                s_vvf[i]=s_smf[i+len];
                  s_vvf[i]=s_smf[i];
                }
            }
       } //k>
     else
       {
          while(get_reply1!=1);
       } 

     if(partial==1) 
       {
         for(i=0;i<len;i++)
           {
             if(k<=s_kmu[i])
               {
//               *(s_vduk+i)=(*(s_vuf+i)- *(s_vufb+i))/ *(s_dzu+i);
                 *(s_vdvk+i)=(*(s_vvf+i)- *(s_vvfb+i))/ *(s_dzu+i);
               }           
             else
               {
                 *(s_vdvk+i)=c0;//*(s_vduk+i)=c0;
               }
           } //i >
//             for(i=0;i<len;i++){
//                if(k>s_kmu[i]) {*(s_vdvk+i)=c0;//*(s_vduk+i)=c0;
//                               }
//                               }
       } //par >
     else
       {
         for(i=0;i<len;i++)
            {
               if(k<=s_kmu[i])
                 {
//                 *(s_vduk+i)=(*(s_vuf+i)- *(s_vufb+i))*dzrk;
                   *(s_vdvk+i)=(*(s_vvf+i)- *(s_vvfb+i))*dzrk;
                 }
               else
                 {
                   *(s_vdvk+i)=c0;//*(s_vduk+i)=c0;
                 }
            } // i >
//              for(i=0;i<len;i++){
//                 if(k>s_kmu[i]) {*(s_vdvk+i)=c0;//*(s_vduk+i)=c0;
//                                }
//                                }
        } //par else >>

//    athread_put(0,s_vduk,vduk+ijst,s_double*len,(void *)&put_reply,0,0);
      athread_put(0,s_vdvk,vdvk+ijst,s_double*len,(void *)&put_reply,0,0);
//    athread_put(0,s_vufb, vuf+ijst,s_double*len,(void *)&put_reply,0,0);
      athread_put(0,s_vvfb, vvf+ijst,s_double*len,(void *)&put_reply,0,0);
      while(put_reply!=2); 
    } // second half >

} 
// s_impt
//

void  s_impt(struct param_vmix *vit)
{
    double *ldm_malloc();
    int volatile get_reply,get_reply0,get_reply1,put_reply;
   
    int nfirst,nlast,nx,ny,nt,km,mt2_r,sfc_ju,partial,thds,ist,ied,jst,jed,myid;
    double grav,p5,c0,aidif;

    double *dz,*c2dtt,*afact,*psfc,*vdc,*dzt,*tnew,*told;
    int *kmt;

    double *p_dz,*p_c2dtt,*p_afact,*p_psfc,*p_vdc,*p_dzt,*p_tnew,*p_told;
    int *p_kmt,*p_kif;
 
    double *p_h1,*p_a,*p_b,*p_e,*p_f;

    int nxy,tail,len,ijst,tlen,tijst,nblocks,block,n,i,k,lenmax,km4,len2,lenk;
    int i1,i2,j1,j2,mt2,off0,off1,len4,dsize,isize;
    double ctmp,dtmp,one,hfact1,hfactk;
    struct param_vmix s_vit;
    one=1.0;


    get_reply=0;
    athread_get(0,vit,&s_vit,sizeof(s_vit),(void *)&get_reply,0,0,0);
    while(get_reply!=1);

    nfirst = s_vit.param[0];
    nlast  = s_vit.param[1];
    nx     = s_vit.param[2];
    ny     = s_vit.param[3];
    nt     = s_vit.param[4];
    km     = s_vit.param[5];
    mt2_r  = s_vit.param[6];
    sfc_ju = s_vit.param[7];
    partial= s_vit.param[8];
    thds   = s_vit.param[9];
    ist    = s_vit.param[10];
    ied    = s_vit.param[11];
    jst    = s_vit.param[12];
    jed    = s_vit.param[13];
    
    grav   = s_vit.dparam[0];
    p5     = s_vit.dparam[1];
    c0     = s_vit.dparam[2];
    aidif  = s_vit.dparam[3];

    dz     = s_vit.addr[0];
    c2dtt  = s_vit.addr[1];
    afact  = s_vit.addr[2];
    psfc   = s_vit.addr[3];
    vdc    = s_vit.addr[4];
    dzt    = s_vit.addr[5];
    tnew   = s_vit.addr[6];
    told   = s_vit.addr[7];
//    ff     = vit->addr[8];

    kmt    = s_vit.addr_int[0];

    myid           = athread_get_id(-1);


       nxy  = nx*ny;
       tail = nxy%thds;
       tlen  = nxy/thds+(tail-myid+thds-1)/thds;  //total-len for myid
       tijst = nxy/thds*myid+min(myid,tail);     //start /myid

     nblocks=(tlen*8*6+tlen*62*4*11)/(60*1024-62*8*4)+1;
     nblocks=1;
#ifdef SMALL
     nblocks=4;
#endif
#ifdef LARGE
     nblocks=1;
#endif

     lenmax = (tlen+nblocks-1)/nblocks;  //for malloc
     km4=km;
     while(km4%4 !=0) km4++;
     len4=(lenmax+3)/4*4;
     len2=(lenmax+lenmax+3)/4*4;
     lenk=(lenmax*km+3)/4*4;
     dsize=km4*3+len4*5+lenk*5;
     isize=len4+lenk;
     while(dsize > 5200|| isize>1024)
     {
      nblocks++;
     lenmax = (tlen+nblocks-1)/nblocks;  //for malloc
     km4=km;
     while(km4%4 !=0) km4++;
     len4=(lenmax+3)/4*4;
     len2=(lenmax+lenmax+3)/4*4;
     lenk=(lenmax*km+3)/4*4;
     dsize=km4*3+len4*5+lenk*5;
     isize=len4+lenk;
     }

     p_dz       = (double*)ldm_malloc(s_double*km4);   // dimension(km)
     p_c2dtt    = (double*)ldm_malloc(s_double*km4);   // dimension(km)
     p_afact    = (double*)ldm_malloc(s_double*km4);
     p_psfc     = (double*)ldm_malloc(s_double*len4);  // (nx_block,ny_block)
     p_vdc      = (double*)ldm_malloc(s_double*len4);  // (nx,ny,0:km+1,2,bid)
     p_dzt      = (double*)ldm_malloc(s_double*lenk);  // (nx,ny,0:km+1,bid)   k,k+1
     p_tnew     = (double*)ldm_malloc(s_double*lenk);  // (nx,ny,km,nt)
     p_told     = (double*)ldm_malloc(s_double*lenk);  // (nx,ny,km,nt)

     p_h1       = (double*)ldm_malloc(s_double*len4);
     p_a        = (double*)ldm_malloc(s_double*len4);
     p_b        = (double*)ldm_malloc(s_double*len4);
//     p_c        = (double*)ldm_malloc(s_double*lenmax);
//     p_d        = (double*)ldm_malloc(s_double*lenmax);
     p_e        = (double*)ldm_malloc(s_double*lenk); //i,j,k               ~32KB
     p_f        = (double*)ldm_malloc(s_double*lenk); //i,j,k  76x52/64x8*62~32KB
     p_kif      = (int *)ldm_malloc(s_int*lenk);
     p_kmt      = (int *)ldm_malloc(s_int*len4);     // (nx,ny,bid)

       get_reply=0;
       athread_get(0,dz,p_dz,s_double*km,(void *)&get_reply,0,0,0);
       athread_get(0,c2dtt,p_c2dtt,s_double*km,(void *)&get_reply,0,0,0);
       athread_get(0,afact,p_afact,s_double*km,(void *)&get_reply,0,0,0);
       while (get_reply!=3);

 for(block=0;block<nblocks;block++)
  { //block <
         tail=tlen%nblocks;
         len=tlen/nblocks+(tail-block+nblocks-1)/nblocks;
         ijst=tlen/nblocks*block+min(block,tail)+tijst;
                   
       get_reply=0;
       get_reply1=0;
       athread_get(0,psfc+ijst,p_psfc,s_double*len,(void *)&get_reply,0,0,0);
       athread_get(0,kmt+ijst,p_kmt,s_int*len,(void *)&get_reply1,0,0,0);

       hfact1 = *(p_dz)/ *(p_c2dtt);

//       if(l_type==l_vart) {
       if(sfc_ju==1) {
                      while (get_reply!=1);
                      for(i=0;i<len;i++){
                       *(p_h1+i) = hfact1 + *(p_psfc+i)/(grav* *(p_c2dtt));
                                        }
                     }
       else          {
                      for(i=0;i<len;i++){
                       *(p_h1+i) = hfact1;
                                        }
                     } 

       while (get_reply1!=1);
           for(i=0;i<len;i++)
            {
              for (k=0;k< (p_kmt[i]-1);k++) { *(p_kif+i+k*len) = 1;}
              for (k=(p_kmt[i]-1);k<km;k++) { *(p_kif+i+k*len) = 0;}
            }

       for(n=(nfirst-1);n<nlast;n++)
       { //n loop left
         mt2 = min (n+1, mt2_r)-1;
         off0=ijst+nxy+nxy*(km+2)*mt2;   //:,:,0,1,bid -> :,:,1,mt2,bid(fortran index)
         off1=ijst+nxy*km*n;             //:,:,1,bid   -> :,:,1,n_c,bid
//         off2=ijst+nxy*km*n;
         get_reply=0; 
         get_reply1=0; 
         get_reply0=0; 
         athread_get(0,vdc+off0,p_vdc,s_double*len,(void *)&get_reply,0,0,0);
         athread_get(0,tnew+off1,p_tnew,s_double*len,(void *)&get_reply,0,0,0);
         if(n==(nfirst-1)) athread_get(0,dzt+ijst,p_dzt,s_double*len,(void *)&get_reply0,0,0,0);
         athread_get(0,told+off1,p_told,s_double*len*km,(void *)&get_reply1,0,(nxy-len)*8,len*8);
//         athread_get(0,dzt+ijst,p_dzt,s_double*len*km,&get_reply1,0,(nxy-len)*8,len*8); 
                            //get dzt(ijst:ijed,2:km+1,bid)

          while(get_reply!=2);
//         if(n==(nfirst-1)) { while(get_reply!=3);}
//                        else { while(get_reply!=2);}

         for(i=0;i<len;i++)
         { //i <
           *(p_a+i) = *(p_afact) * *(p_vdc+i);           
           dtmp = one/(*(p_h1+i) + *(p_a+i));
           *(p_e+i) = *(p_a+i) *dtmp;
           *(p_b+i) = *(p_h1+i) * *(p_e+i);
           *(p_f+i) = hfact1 * *(p_tnew+i)*dtmp; 
         } //i >
        
         if(n==(nfirst-1)) { while(get_reply0!=1);}
       if(partial==1)
         { //if par <
         for(k=1;k<km;k++)
           {   //k <
            get_reply=0;
            athread_get(0,vdc+off0+k*nxy,p_vdc,s_double*len,(void *)&get_reply,0,0,0);
            athread_get(0,tnew+off1+k*nxy,p_tnew+k*len,s_double*len,(void *)&get_reply,0,0,0);
            if(n==(nfirst-1)) athread_get(0,dzt+ijst+k*nxy,p_dzt+k*len,s_double*len,(void *)&get_reply,0,0,0); 

         if(n==(nfirst-1)) { while(get_reply!=3);}
                        else { while(get_reply!=2);}

             for(i=0;i<len;i++)
              {
                ctmp = *(p_a+i);
                *(p_a+i) = aidif* *(p_vdc+i)/(p5 * (*(p_dzt+i+(k-1)*len)+ *(p_dzt+i+k*len)));
                hfactk   = *(p_dzt+i+(k-1)*len)/ *(p_c2dtt+k);        
                dtmp   = one/(hfactk + *(p_kif+i+len*k)* *(p_a+i) + *(p_b+i));
                *(p_e+k*len+i)= *(p_a+i) *dtmp;
                *(p_b+i)   = (hfactk + *(p_b+i))* *(p_e+k*len+i);
                *(p_f+k*len+i)= (hfactk* *(p_tnew+i+k*len)+ ctmp* *(p_f+(k-1)*len+i))*dtmp;
              }
           } // k >
             for(i=0;i<len;i++)
              {
                for(k=p_kmt[i];k<km;k++) { *(p_f+k*len+i)=c0;}
              }
         } //if par >
       else
         { //if par else <
         for(k=1;k<km;k++)
           {   //k <
            get_reply=0;
//            get_reply1=0;
            athread_get(0,vdc+off0+k*nxy,p_vdc,s_double*len,(void *)&get_reply,0,0,0);
            athread_get(0,tnew+off1+k*nxy,p_tnew+k*len,s_double*len,(void *)&get_reply,0,0,0);
//            if(n==(nfirst-1)) athread_get(0,dzt+ijst+k*nxy,p_dzt+k*len,s_double*len,(void *)&get_reply,0,0,0); 

//         if (n==(nfirst-1)) {while(get_reply!=3);}
//                else { while(get_reply!=2);}       
             while(get_reply!=2);       

             hfactk = *(p_dz+k)/ *(p_c2dtt+k);
             for(i=0;i<len;i++)
               {
                ctmp = *(p_a+i);
                *(p_a+i) = *(p_afact+k)* *(p_vdc+i);
                dtmp = one/(hfactk + *(p_kif+i+k*len) * *(p_a+i)+ *(p_b+i));
                *(p_e+k*len+i) = *(p_a+i)*dtmp;
                *(p_b+i) = (hfactk + *(p_b+i))* *(p_e+k*len+i);
                *(p_f+k*len+i) = (hfactk * *(p_tnew+k*len+i) +ctmp* *(p_f+(k-1)*len+i))*dtmp;
               }
           }   //k >
             for(i=0;i<len;i++)
               {
                for(k=p_kmt[i];k<km;k++) { *(p_f+k*len+i)=c0;}
               }
         } //if par else >

      
         while(get_reply1!=1);

         put_reply=0;
//          athread_put(0,p_f+km*len-len,ff+km*nxy-nxy+n*nxy*km+ijst,s_double*len,&put_reply,0,0);

//           for(i=0;i<len;i++)
//            {
//                for(k=(p_kmt[i]-2);k>-1;k--)
//                {
//                 *(p_f+k*len+i)+=  *(p_e+k*len+i)* *(p_f+k*len+len+i);
//                }
//            }
            
           for(k=km-2;k>-1;k--)
             {
                for(i=0;i<len;i++)
                 {
                   *(p_f+k*len+i)+=  *(p_e+k*len+i)* *(p_f+k*len+len+i)*p_kif[k*len+i];
                 }
             } 

          for(k=0;k<km;k++)
          {       // only adapt to len<=nx
              j1=ijst/nx;
              j2=(ijst+len-1)/nx;
              if(j1==j2)
                {
                   if( j1 >= jst-1 && j1<jed)
                     {
                        i1=j1*nx+ist-1-ijst;i2=j1*nx+ied-ijst;
                        for(i=max(i1,0);i<min(i2,len);i++)
                        {
                          p_tnew[k*len+i]=p_told[k*len+i]+p_f[k*len+i]; 
                        } 

                     }
                }
              else   // j2=j1+1
                {
                   if( j1 >= jst-1  )
                     {
                       if(j2 <= jed)
                         { 
                           i1=j1*nx+ist-1-ijst;i2=j1*nx+ied-ijst;
                           for(i=max(i1,0);i<min(i2,len);i++) p_tnew[k*len+i]=p_told[k*len+i]+p_f[k*len+i]; 
                         } 
                       if(j2<jed)
                         {
                           i1=j2*nx+ist-1-ijst;i2=j2*nx+ied-ijst;
                           for(i=max(i1,0);i<min(i2,len);i++) p_tnew[k*len+i]=p_told[k*len+i]+p_f[k*len+i]; 
                         }
		     } 
                   else
                     {
                       if( j1=jst-2)
                         {
                           i1=j2*nx+ist-1-ijst;i2=j2*nx+ied-ijst;
                           for(i=max(i1,0);i<min(i2,len);i++) p_tnew[k*len+i]=p_told[k*len+i]+p_f[k*len+i]; 
                         }
                     }

                }

//           j1=ijst/nx;j2=(ijst+len-1)/nx;
//           for(j=max(j1,jst-1);j<min(j2+1,jed);j++)
//             {
//                i1=j*nx+ist-1-ijst;i2=j*nx+ied-ijst;
//                  for(i=max(i1,0);i<min(i2,len);i++) p_tnew[k*len+i]=p_told[k*len+i]+p_f[k*len+i]; 
//             }

//            for(i=0;i<len;i++)
//            { 
//              io=(ijst+i)%nx; jo=(ijst+i)/nx;
//              if(io>=(ist-1) && io<ied && jo>=(jst-1) && jo<jed) p_tnew[k*len+i]=p_told[k*len+i]+p_f[k*len+i];
//            }
            athread_put(0,p_tnew+k*len,tnew+k*nxy+n*nxy*km+ijst,s_double*len,(void *)&put_reply,0,0);
          }  //k>
         
//       for(k=0;k<km;k++)
//         {
//          for(i=0;i<len;i++)
//           {
//             p_tnew[i]=p_told[i]+ *(p_f+k*len+i);
//           }
//         }
//         athread_put(0,p_f,ff+off1,s_double*len*km,&put_reply,(nxy-len)*8,len*8);
//         while (put_reply !=1);

//        for(k=0;k<km;k++)
//         {
//          athread_put(0,p_f+k*len,ff+k*nxy+n*nxy*km+ijst,s_double*len,&put_reply,0,0);
//         }
        while(put_reply !=km);

       } //n loop right

  } //blcok >

     ldm_free(p_dz,s_double*km4);   // dimension(km)
     ldm_free(p_c2dtt,s_double*km4);   // dimension(km)
     ldm_free(p_psfc,s_double*len4);   // dimension(km)
     ldm_free(p_vdc,s_double*len4);   
     ldm_free(p_dzt,s_double*lenk);   
     ldm_free(p_tnew,s_double*lenk);   
     ldm_free(p_told,s_double*lenk);   
     ldm_free(p_kmt,s_int*len4);   
     ldm_free(p_afact,s_double*km4);   

// temp local
     ldm_free(p_h1,s_double*len4);   
     ldm_free(p_a,s_double*len4);   
     ldm_free(p_b,s_double*len4);   
//     ldm_free(p_c,s_double*lenmax);   
//     ldm_free(p_d,s_double*lenmax);   
     ldm_free(p_e,s_double*lenk);   
     ldm_free(p_f,s_double*lenk);   

     ldm_free(p_kif,s_int*lenk);   

}




//################  vdifft  #########################
void  s_vdift(struct param_vmix *vdift)
{
    double *ldm_malloc();
    int volatile get_reply,get_reply0,get_reply1,put_reply;

    int nx,ny,k,km,kvdc,mt2_r,partial,lbtmhf,nt,thds,myid;
    double c0,p5,ztk,dzwrk,zwk,dzrk,bm_h_flx,bm_h_flxdp;

    double *stf,*vdc,*dzt,*told,*vtf,*vdtk;
    int *kmt;

    double *s_stf,*s_vdc,*s_dzt,*s_told,*s_vtf,*s_vdtk;
    int *s_kmt;

    double *s_vtfb,*s_work;
    int nxy,tail,len,ijst,le4;
    double zero;
    int rank,mt2,n,i,off0,off1,off2;
    struct param_vmix s_vdift;
    zero=0.0;

    get_reply=0;
    athread_get(0,vdift,&s_vdift,sizeof(s_vdift),(void *)&get_reply,0,0,0);

    while(get_reply!=1);

    nx           = s_vdift.param[0]; 
    ny           = s_vdift.param[1]; 
    k            = s_vdift.param[2]; 
    km           = s_vdift.param[3]; 
    kvdc         = s_vdift.param[4]; 
    mt2_r        = s_vdift.param[5]; 
    partial      = s_vdift.param[6]; 
    lbtmhf       = s_vdift.param[7]; 
    nt           = s_vdift.param[8]; 
    thds         = s_vdift.param[9];
    rank         = s_vdift.param[10];

    myid         = athread_get_id(-1);

    c0           = s_vdift.dparam[0];
    p5           = s_vdift.dparam[1];
    ztk          = s_vdift.dparam[2];
    dzwrk        = s_vdift.dparam[3];
    zwk          = s_vdift.dparam[4];
    dzrk         = s_vdift.dparam[5];
    bm_h_flx     = s_vdift.dparam[6];
    bm_h_flxdp   = s_vdift.dparam[7];

    stf          = s_vdift.addr[0];
    vdc          = s_vdift.addr[1];
    dzt          = s_vdift.addr[2];
    told         = s_vdift.addr[3];
    vtf          = s_vdift.addr[4];
    vdtk         = s_vdift.addr[5];

    kmt          = s_vdift.addr_int[0];

    nxy=nx*ny;
    tail=nxy%thds;
    len=nxy/thds+(tail-myid+thds-1)/thds;
    ijst=nxy/thds*myid+min(myid,tail);
    le4=((len+3)>>2)<<2;

    s_stf  =  (double*)ldm_malloc(s_double*le4);
    s_vtf  =  (double*)ldm_malloc(s_double*le4);
    s_vdc  =  (double*)ldm_malloc(s_double*le4);
    s_dzt  =  (double*)ldm_malloc(s_double*le4*2);
    s_told =  (double*)ldm_malloc(s_double*le4*2);
    s_work =  (double*)ldm_malloc(s_double*le4);
    s_vdtk =  (double*)ldm_malloc(s_double*le4);
    s_vtfb =  (double*)ldm_malloc(s_double*le4);  //for vtfb
    s_kmt  =  (int *)ldm_malloc(s_int*le4);

      get_reply=0;
      get_reply0=0;
      off2=k*nxy+ijst;
      athread_get(0,kmt+ijst,s_kmt,len*s_int,(void *)&get_reply,0,0,0);
//      athread_get(0,dzt+off2,s_dzt,len*s_double,&get_reply1,0,0,0);
//      if(k<km)   athread_get(0,dzt+off2+nxy,s_dzt+len,len*s_double,&get_reply1,0,0,0);
      if(k<km) 
        {
           athread_get(0,dzt+off2,s_dzt,len*s_double*2,(void *)&get_reply0,0,(nxy-len)*s_double,len*s_double);
        }
      else
        {
           athread_get(0,dzt+off2,s_dzt,len*s_double,(void *)&get_reply0,0,0,0);
        }
          while(get_reply!=1);  // for   kmt

     for(n=0;n<nt;n++)
       {
//      io=n+2;  // for while kmt(1) + 1*(n-0+1) stf or vtf
      get_reply=0;   //for dzt indepof n vars
      get_reply1=0;
      mt2 = min(n+1,mt2_r)-1;
//!      if (k == 1) then
//!         VTF(:,:,n,bid) = merge(STF(:,:,n), c0, KMT(:,:,bid) >= 1)
//!      endif
//!  need STF(:,:,n),KMT(:,:,bid),c0, malloc vtf,kmt
//
      if (k == 1) 
       { // k<
    
          athread_get(0,stf+n*nxy+ijst,s_stf,len*s_double,(void *)&get_reply,0,0,0);
          while(get_reply!=1);
          for(i=0;i<len;i++) 
           {
             if(s_kmt[i]>=1) {s_vtf[i] = s_stf[i];}
             else { s_vtf[i]=c0;}
           }
//           {
//             s_vtf[i] = s_stf[i];
//           }
//          for(i=0;i<len;i++) 
//           {
//             if(s_kmt[i] < 1) s_vtf[i]=c0;
//           }
        } // k=1 >
       else  // get vtf for k>1
        {
          athread_get(0,vtf+n*nxy+ijst,s_vtf,len*s_double,(void *)&get_reply,0,0,0);
        }       


//     if(myid==0&& rank==0) cpe_printf("a2 get_reply= k=io n %d %d  %d %d\n",get_reply,k,io,n);
     
//!-----------------------------------------------------------------------
//!
//!     vertical tracer flux vdc*Dz(T) at bottom  of T box
//!     calculate Dz(vdc*Dz(T))
//!
//!-----------------------------------------------------------------------
//!!     need vdc(:,:,kvdc,mt2,bid) told(,,k,n kp1,n), p5,
//!      dzt(k,bid)(kp1,bid),zt(k)--consants,
//!      bottom_heat_flx,bottom_heat_flx_depth
//        if(k==km)
//         {
//           block=1;  // use for reply 
//           for(i=0;i<len;i++) s_a[i]=zero;
//         }
        if(k<km) 
         {
//           block=3*(1+n)+2;   //use for the following reply
//         off0=mt2*(km+2)*nxy+kvdc*nxy+ijst;
           off0=mt2*(km+2)*nxy+kvdc*nxy+ijst; 
           off1=n*km*nxy+(k-1)*nxy+ijst;
           athread_get(0,vdc+off0,s_vdc,len*s_double,(void *)&get_reply1,0,0,0);
           athread_get(0,told+off1,s_told,len*s_double*2,(void *)&get_reply1,0,(nxy-len)*8,len*8);
//           athread_get(0,told+off1,s_told,len*s_double,&get_reply1,0,0,0);
//           athread_get(0,told+off1+nxy,s_told+len,len*s_double,&get_reply1,0,0,0);
         }

//         if(myid==0 && rank==0) cpe_printf("a3 get_reply1= k=  block %d %d %d\n",get_reply1,block,k,n);

    if (partial==1) 
      {
//      if (partial_bottom_cells) then
//!CDIR COLLAPSE
//!         VTFB = merge(VDC(:,:,kvdc,mt2,bid)*                    &
//!                      (TOLD(:,:,k  ,n) - TOLD(:,:,kp1,n))/      &
//!                      (p5*(DZT(:,:,k,bid) + DZT(:,:,kp1,bid)))  &
//!                      ,c0, KMT(:,:,bid) > k)
//         if(myid==0 && rank==0) cpe_printf("a3 get_reply1= k=  block %d %d %d\n",get_reply1,k,block);
        if (n==0) while(get_reply0 !=1);
        if (k<km) while(get_reply1 !=2);
//         if(myid==0 && rank==0) {
//          cpe_printf("b3 get_reply1= k=  block %d %d %d\n",get_reply1,k,block);
//                               }

        if(k==km)
         {
           for(i=0;i<len;i++) s_vtfb[i]=zero;  //for vtfb
         }
        else
         {
           for(i=0;i<len;i++)
            {
              if(s_kmt[i]>k)
                {
                  s_vtfb[i]=s_vdc[i]*(s_told[i]-s_told[i+len])/(p5*(s_dzt[i]+s_dzt[i+len])); 
                }
               else{s_vtfb[i]=c0;}
            } 
         }
//            {
//              s_a[i]=s_vdc[i]*(s_told[i]-s_told[i+len])/(p5*(s_dzt[i]+s_dzt[i+len])); 
//            } 
//         } 
//         for(i=0;i<len;i++)
//           {
//             if(s_kmt[i] <=k) s_a[i]=c0;
//           }

//       if(k==km) then
//         vtfb=0.d0
//       else 
//         VTFB= VDC(:,:,kvdc,mt2,bid)*                    &
//               (TOLD(:,:,k  ,n) - TOLD(:,:,kp1,n))/      &
//               (p5*(DZT(:,:,k,bid) + DZT(:,:,kp1,bid)))  
//       endif
//         do i=1,nx_block*ny_block
//            if(kmt(i,1,bid) <=k) then
//              VTFB(i,1)=c0
//            endif
//         enddo
//
         if(lbtmhf==1 && n==0) 
          {
            for(i=0;i<len;i++)
             {
               s_work[i]=ztk+p5*s_dzt[i];
             }
            for(i=0;i<len;i++)
             {
               if(k==s_kmt[i] && s_work[i]>=bm_h_flxdp) s_vtfb[i]=-bm_h_flx;
             }
          }
//         if (lbottom_heat_flx .and. n == 1) then
//!CDIR COLLAPSE
//!            VTFB = merge( -bottom_heat_flx, VTFB,       &
//!                         k == KMT(:,:,bid) .and.        &
//!                         (zt(k) + p5*DZT(:,:,k,bid)) >= &
//!                         bottom_heat_flx_depth)
//
//            work=zt(k) + p5*DZT(:,:,k,bid)
//            do i=1,nx_block*ny_block
//              if(k == kmt(i,1,bid) .and. work(i,1)>=bottom_heat_flx_depth) then
//                vtfb(i,1)=-bottom_heat_flx
//              endif
//            enddo
//         endif
//!    need vtf,
    if(k>1) while(get_reply!=1);
       for(i=0;i<len;i++)
        {
          if(k<=s_kmt[i]) 
            {
              s_vdtk[i]=(s_vtf[i]-s_vtfb[i])/s_dzt[i]; 
            }
          else
            {
              s_vdtk[i]=c0;
            }
        }
//       for(i=0;i<len;i++)
//        {
//          s_vdtk[i]=(s_vtf[i]-s_a[i])/s_dzt[i]; 
//        }
//       for(i=0;i<len;i++)
//        {
//         if(k>s_kmt[i]) s_vdtk[i]=c0;
//        }
//!CDIR COLLAPSE
//!         VDTK(:,:,n) = merge((VTF(:,:,n,bid) - VTFB)/DZT(:,:,k,bid), &
//!                             c0, k <= KMT(:,:,bid))
//         VDTK(:,:,n) = (VTF(:,:,n,bid) - VTFB)/DZT(:,:,k,bid)
//         do i=1,nx_block*ny_block
//           if(k>kmt(i,1,bid)) vdtk(i,1,n)=c0
//         enddo
      }  // if partial >
//      else
    else
      {   // < partial else 
        if (n==0) while(get_reply0 !=1);
        if (k<km) while(get_reply1 !=2);
//!CDIR COLLAPSE
//!need dzwrk,zwk
//!         VTFB = merge(VDC(:,:,kvdc,mt2,bid)*                      &
//!                      (TOLD(:,:,k  ,n) - TOLD(:,:,kp1,n))*dzwr(k) &
//!                      ,c0, KMT(:,:,bid) > k)
        if(k==km)
         {
          for(i=0;i<len;i++) s_vtfb[i]=zero;
         }
        else
         {
          for(i=0;i<len;i++) s_vtfb[i]=s_vdc[i]*(s_told[i]-s_told[i+len])*dzwrk;
         }
         for(i=0;i<len;i++)
          {
            if(s_kmt[i]<=k) s_vtfb[i]=c0;
          }
//
//         if(k==km) then
//           vtfb=0.d0
//         else
//           VTFB = VDC(:,:,kvdc,mt2,bid)*                      &
//                  (TOLD(:,:,k  ,n) - TOLD(:,:,kp1,n))*dzwr(k)
//         endif
//         do i=1,nx_block*ny_block
//           if(kmt(i,1,bid) <=k) vtfb(i,1)=c0
//         enddo
         if(lbtmhf==1 && n==0)
           {
            for(i=0;i<len;i++)
             {
               if(k==s_kmt[i] && zwk>=bm_h_flxdp) s_vtfb[i]=-bm_h_flx;
             }
           }
          
//         if (lbottom_heat_flx .and. n == 1) then
//!CDIR COLLAPSE
//!            VTFB = merge( -bottom_heat_flx, VTFB,      &
//!                         k == KMT(:,:,bid) .and.       &
//!                         zw(k) >= bottom_heat_flx_depth)
//           do i=1,nx_block*ny_block
//             if(k==kmt(i,1,bid) .and. zw(k)>=bottom_heat_flx_depth) then
//               VTFB(i,1)=-bottom_heat_flx
//             endif
//           enddo
//         endif

//! put out
//!CDIR COLLAPSE
//!         VDTK(:,:,n) = merge((VTF(:,:,n,bid) - VTFB)*dzr(k), &
//!                             c0, k <= KMT(:,:,bid))
          if(k>1) while(get_reply!=1);
          for(i=0;i<len;i++) s_vdtk[i]=(s_vtf[i]-s_vtfb[i])*dzrk;
          for(i=0;i<len;i++)
           {
            if(k>s_kmt[i]) s_vdtk[i]=c0;
           } 


//         VDTK(:,:,n) = (VTF(:,:,n,bid) - VTFB)*dzr(k)
//         do i=1,nx_block*ny_block
//           if(k>kmt(i,1,bid)) vdtk(i,1,n)=c0
//         enddo
//      endif
       } // partial end >

    put_reply=0;
    athread_put(0,s_vdtk,vdtk+n*nxy+ijst,len*s_double,(void *)&put_reply,0,0);
    athread_put(0,s_vtfb,vtf+n*nxy+ijst,len*s_double,(void *)&put_reply,0,0);
    while(put_reply!=2);

//!-----------------------------------------------------------------------
//!
//!     set top value of VTF to bottom value for next pass at level k+1
//!
//!-----------------------------------------------------------------------
//!!   put vtfb to vtf
//!CDIR COLLAPSE
//      VTF(:,:,n,bid) = VTFB  
//
//   enddo

   } // n>

//      if(myid==0 && rank==0) cpe_printf("end loopn,nx,ny,nt,km % %d %d %d\n",nx,ny,nt,km);

    ldm_free(s_stf,s_double*le4);
    ldm_free(s_kmt,s_int*le4);
    ldm_free(s_vtf,s_double*le4);
    ldm_free(s_vdc,s_double*le4);
    ldm_free(s_told,s_double*le4*2);
    ldm_free(s_dzt,s_double*le4*2);
    ldm_free(s_work,s_double*le4);
    ldm_free(s_vdtk,s_double*le4);
    ldm_free(s_vtfb,s_double*le4);

}


void  s_impu(struct param_vmix *vdiu)
{
   double *ldm_malloc();
   int volatile get_reply,get_reply0,get_reply1,put_reply;

   int nx,ny,km,partial,ist,ied,jst,jed,thds;
   double c0,p5,aidif,c2dtu;
   double *dz,*vvc,*unew,*vnew,*dzu,*afacu,*uold,*vold;
   int *kmu;

   double *s_dz,*s_vvc,*s_unew,*s_vnew,*s_dzu,*s_afacu;
   int *s_kmu,*s_kif;

   int nxy,len,tail,ijst,tlen,tijst,lenmax,nblocks,block,k,i;

   double *s_a,*s_b,*s_e,*s_f1,*s_f2,hfacu1,hfacuk;
   int j,i1,j1,i2,j2,rank,myid,current,next,km4,dsize;
   double dtmp,one,ctmp;
   struct param_vmix s_vdiu;
   one=1.0;


    get_reply=0;
    athread_get(0,vdiu,&s_vdiu,sizeof(s_vdiu),(void *)&get_reply,0,0,0);
    while(get_reply!=1);

    nx           = s_vdiu.param[0]; 
    ny           = s_vdiu.param[1]; 
    km           = s_vdiu.param[2]; 
    partial      = s_vdiu.param[3]; 
    ist          = s_vdiu.param[4]; 
    ied          = s_vdiu.param[5]; 
    jst          = s_vdiu.param[6]; 
    jed          = s_vdiu.param[7]; 
    thds         = s_vdiu.param[8];
    rank         = s_vdiu.param[9];
  
    c0           = s_vdiu.dparam[0]; 
    p5           = s_vdiu.dparam[1]; 
    aidif        = s_vdiu.dparam[2]; 
    c2dtu        = s_vdiu.dparam[3]; 

    dz           = s_vdiu.addr[0];
    vvc          = s_vdiu.addr[1];
    unew         = s_vdiu.addr[2];
    vnew         = s_vdiu.addr[3];
    uold         = s_vdiu.addr[4];
    vold         = s_vdiu.addr[5];
    dzu          = s_vdiu.addr[6];
    afacu        = s_vdiu.addr[7];   //use as afac_u

    kmu          = s_vdiu.addr_int[0];
    
    myid           = athread_get_id(-1);

#ifdef POPchechJP    
    if(myid==0&& rank==0) 
     {
       cpe_printf("nx,ny,km,partial %d %d  %d %d\n",nx,ny,km,partial);
     }
#endif
    nxy    = nx*ny;
    tail   = nxy%thds;
    tlen   = nxy/thds+(tail-myid+thds-1)/thds;  //total-len for myid
    tijst  = nxy/thds*myid+min(myid,tail);     //start /myid
    nblocks = (nxy/thds*8*km*11/2+nxy/thds*8*10)/(60*1024-62*8*7/2)+1;
    nblocks=1;
#ifdef SMALL
     nblocks=4;
#endif
#ifdef LARGE
     nblocks=1;
#endif

    lenmax = (tlen+nblocks-1)/nblocks;
    while(lenmax%4 != 0) lenmax++;
    km4=km;
    while(km4%4 !=0) km4++;
    dsize=km4+km4+lenmax*3+lenmax*km*5;

    while (dsize > 5500) 
    {
      nblocks++;
      lenmax = (tlen+nblocks-1)/nblocks;
      while(lenmax%4 != 0) lenmax++;
      km4=km;
      while(km4%4 !=0) km4++;
      dsize=km4+km4+lenmax*5+lenmax*km*5;
    }

     s_vvc     = ldm_malloc(s_double*lenmax);     // (nx,ny,k,bid) (:,:,km,max)
     s_unew    = ldm_malloc(s_double*lenmax*km);     // (nx,ny,k,bid)
     s_vnew    = ldm_malloc(s_double*lenmax*km);     // (nx,ny,k,bid)
     s_dzu     = ldm_malloc(s_double*lenmax*2);     // (nx,ny,k:k+1,bid) (:,:,0:km+1,max_c)
     s_a       = ldm_malloc(s_double*lenmax);     // (nx,ny)
     s_b       = ldm_malloc(s_double*lenmax);     // (nx,ny)
//     s_c       = ldm_malloc(s_double*lenmax);     // (nx,ny)
//     s_d       = ldm_malloc(s_double*lenmax);     // (nx,ny)
     s_e       = ldm_malloc(s_double*lenmax*km);     // (nx,ny,km)
     s_f1      = ldm_malloc(s_double*lenmax*km);     // (nx,ny,km)
     s_f2      = ldm_malloc(s_double*lenmax*km);     // (nx,ny,km)
     s_dz      = ldm_malloc(s_double*km4);
     s_afacu   = ldm_malloc(s_double*km4);     //afacu (km)
//     s_hfacu   = ldm_malloc(s_double*km);     //hfac_u (km)
     s_kif     = (int *)ldm_malloc(s_int*lenmax*km);     
     s_kmu     = (int *)ldm_malloc(s_int*lenmax);     

       get_reply=0;
       get_reply0=0;
       athread_get(0,dz,s_dz,s_double*km,(void *)&get_reply,0,0,0);
       athread_get(0,afacu,s_afacu,s_double*km,(void *)&get_reply0,0,0,0); // afacu(km)
       while (get_reply!=1);
//       for(k=0;k<km;k++) {s_hfacu[k]=s_dz[k]/c2dtu;}
         hfacu1=s_dz[0]/c2dtu;
//   do k=1,km
//      hfac_u(k) = dz(k)/c2dtu
//   end do


 for(block=0;block<nblocks;block++)
  { //block <
         tail=tlen%nblocks;
         len=tlen/nblocks+(tail-block+nblocks-1)/nblocks;
         ijst=tlen/nblocks*block+min(block,tail)+tijst;
          
         get_reply=0;
         athread_get(0,vvc+ijst,s_vvc,s_double*len,(void *)&get_reply,0,0,0);
         athread_get(0,unew+ijst,s_unew,s_double*len,(void *)&get_reply,0,0,0);
         athread_get(0,vnew+ijst,s_vnew,s_double*len,(void *)&get_reply,0,0,0);
         get_reply1=0;
//         athread_get(0,dzu+ijst+2*nxy,s_dzu,s_double*len*km,&get_reply,0,(nxy-len)*8,len*8);
         athread_get(0,dzu+ijst+2*nxy,s_dzu,s_double*len,(void *)&get_reply1,0,0,0);
         athread_get(0,kmu+ijst,s_kmu,s_int*len,(void *)&get_reply1,0,0,0);
         if(block==0) while (get_reply0!=1);
         while (get_reply!=3);
         for(i=0;i<len;i++)
           {
              s_a[i]=s_afacu[0]*s_vvc[i];
//              dtmp=one/(s_hfacu[0]+s_a[i]);
              dtmp=one/(hfacu1+s_a[i]);
              s_e[i]=s_a[i]*dtmp;
//              s_b[i]=s_hfacu[0]*s_e[i];
              s_b[i]=hfacu1*s_e[i];
//              s_f1[i]=s_hfacu[0]*s_unew[i]*dtmp;
              s_f1[i]=hfacu1*s_unew[i]*dtmp;
//              s_f2[i]=s_hfacu[0]*s_vnew[i]*dtmp;
              s_f2[i]=hfacu1*s_vnew[i]*dtmp;
           }

         while (get_reply1!=2);
           for(i=0;i<len;i++)
            {
              for (k=1;k< (s_kmu[i]-1);k++) { *(s_kif+i+k*len) = 1;}
              for (k=(s_kmu[i]-1);k<km;k++) { *(s_kif+i+k*len) = 0;}
            }

         if(partial==1) 
           {
              for(k=1;k<km;k++)
                {
                  current=(k-1)%2; next=k%2;get_reply=0;
                  athread_get(0,dzu+ijst+(k+2)*nxy,s_dzu+len*next,s_double*len,(void *)&get_reply,0,0,0);
                  athread_get(0,vvc+ijst+k*nxy,s_vvc,s_double*len,(void *)&get_reply,0,0,0);
                  athread_get(0,unew+ijst+k*nxy,s_unew+k*len,s_double*len,(void *)&get_reply,0,0,0);
                  athread_get(0,vnew+ijst+k*nxy,s_vnew+k*len,s_double*len,(void *)&get_reply,0,0,0);
                  while (get_reply!=4);
                  for(i=0;i<len;i++)
                    {
                      ctmp=s_a[i];
//		      s_hfacu[k]=s_dzu[len*current+i]/c2dtu;
		      hfacuk=s_dzu[len*current+i]/c2dtu;
                      s_a[i]=aidif*s_vvc[i]/(p5*(s_dzu[len*current+i]+s_dzu[len*next+i]));
//                      dtmp=one/(s_hfacu[k]+s_kif[k*len+i]*s_a[i]+s_b[i]);
                      dtmp=one/(hfacuk+s_kif[k*len+i]*s_a[i]+s_b[i]);
                      s_e[k*len+i]=s_a[i]*dtmp;
//                      s_b[i]=(s_hfacu[k]+s_b[i])*s_e[k*len+i];
                      s_b[i]=(hfacuk+s_b[i])*s_e[k*len+i];
//                      s_f1[k*len+i]=(s_hfacu[k]*s_unew[i+k*len]+ctmp*s_f1[len*k-len+i])*dtmp;
                      s_f1[k*len+i]=(hfacuk*s_unew[i+k*len]+ctmp*s_f1[len*k-len+i])*dtmp;
//                      s_f2[k*len+i]=(s_hfacu[k]*s_vnew[i+k*len]+ctmp*s_f2[len*k-len+i])*dtmp;
                      s_f2[k*len+i]=(hfacuk*s_vnew[i+k*len]+ctmp*s_f2[len*k-len+i])*dtmp;
                    }
                }  //k >
              for(i=0;i<len;i++)
                {
                  for(k=s_kmu[i];k<km;k++)
                    {s_f1[k*len+i]=c0; s_f2[k*len+i]=c0;
                    }
                }
           } // if partial>



         else //if no partial
           {
              for(k=1;k<km;k++)
                {
                  hfacuk=s_dz[k]/c2dtu;
                  current=(k-1)%2; next=k%2;get_reply=0;
                  athread_get(0,dzu+ijst+k*nxy,s_dzu+len*next,s_double*len,(void *)&get_reply,0,0,0);
                  athread_get(0,vvc+ijst+k*nxy,s_vvc,s_double*len,(void *)&get_reply,0,0,0);
                  athread_get(0,unew+ijst+k*nxy,s_unew+k*len,s_double*len,(void *)&get_reply,0,0,0);
                  athread_get(0,vnew+ijst+k*nxy,s_vnew+k*len,s_double*len,(void *)&get_reply,0,0,0);
                  while (get_reply!=4);
                  for(i=0;i<len;i++)
                    {
                      ctmp=s_a[i];
                      s_a[i]=s_afacu[k]*s_vvc[i];
//                      dtmp=one/(s_hfacu[k]+s_a[i]*s_kif[k*len+i]+s_b[i]);
                      dtmp=one/(hfacuk+s_a[i]*s_kif[k*len+i]+s_b[i]);
                      s_e[i]=s_a[i]*dtmp;
                      s_b[i]=(hfacuk+s_b[i])*s_e[len*k+i];
                      s_f1[i]=(hfacuk*s_unew[i+k*len]+ctmp*s_f1[k*len-len+i])*dtmp;
                      s_f2[i]=(hfacuk*s_vnew[i+k*len]+ctmp*s_f2[k*len-len+i])*dtmp;
                    }
                }  // k >
        
              for(i=0;i<len;i++)
                {
                  for(k=s_kmu[i];k<km;k++)
                    {
                    s_f1[k*len+i]=c0; s_f2[k*len+i]=c0;
                    }
                }
           } //  > else




       for(i=0;i<len;i++)
         {
           for(k=s_kmu[i]-2;k>-1;k--)
            {
              s_f1[k*len+i]+= s_e[i+k*len]*s_f1[k*len+len+i];
              s_f2[k*len+i]+= s_e[i+k*len]*s_f2[k*len+len+i];
            }
         }

       
       j1=ijst/nx;j2=(ijst+len-1)/nx;
       put_reply=0;
       for(k=0;k<km;k++)
         {
//    for uvew=new+old
           get_reply=0;
           athread_get(0,uold+ijst+k*nxy,s_a,s_double*len,(void *)&get_reply,0,0,0);
           athread_get(0,vold+ijst+k*nxy,s_b,s_double*len,(void *)&get_reply,0,0,0);
//    old+new >-

           for(j=max(j1,jst-1);j<min(j2+1,jed);j++)
             {
              i1=j*nx+ist-1-ijst;i2=j*nx+ied-ijst;
              for(i=max(i1,0);i<min(i2,len);i++) 
               {
                   s_unew[k*len+i]=s_f1[k*len+i];
                   s_vnew[k*len+i]=s_f2[k*len+i];
               }
             }
//   for new+old
           while(get_reply != 2);
           for(i=0;i<len;i++)
             {
               s_unew[k*len+i]+=s_a[i];
               s_vnew[k*len+i]+=s_b[i];
             }
//   new+old >-
         athread_put(0,s_unew+k*len,unew+ijst+k*nxy,len*s_double,(void *)&put_reply,0,0);
         athread_put(0,s_vnew+k*len,vnew+ijst+k*nxy,len*s_double,(void *)&put_reply,0,0);
         }

      while(put_reply!=km*2);

  }  //block >
     ldm_free(s_dz,s_double*km4);
     ldm_free(s_vvc,s_double*lenmax);     // (nx,ny,k,bid) (:,:,km,max)
     ldm_free(s_unew,s_double*lenmax*km);     // (nx,ny,k,bid)
     ldm_free(s_vnew,s_double*lenmax*km);     // (nx,ny,k,bid)
     ldm_free(s_dzu,s_double*lenmax*2);     // (nx,ny,k:k+1,bid) (:,:,0:km+1,max_c)
     ldm_free(s_a,s_double*lenmax);     // (nx,ny)
     ldm_free(s_b,s_double*lenmax);     // (nx,ny)
//     ldm_free(s_c,s_double*lenmax);     // (nx,ny)
//     ldm_free(s_c,s_double*lenmax);     // (nx,ny)
//     ldm_free(s_d,s_double*lenmax);     // (nx,ny)
     ldm_free(s_e,s_double*lenmax*km);     // (nx,ny,km)
     ldm_free(s_f1,s_double*lenmax*km);     // (nx,ny,km)
     ldm_free(s_f2,s_double*lenmax*km);     // (nx,ny,km)
//     ldm_free(s_hfacu,s_double*km);     //hfac_u (km)
     ldm_free(s_afacu,s_double*km4);     //afacu (km)
     ldm_free(s_kif,s_int*lenmax*km);
     ldm_free(s_kmu,s_int*lenmax);

}

void  s_impt_cr(struct param_vmix *vitc)
{
   double *ldm_malloc();
   int volatile get_reply,get_reply0,get_reply1,put_reply;
   int dsize;

   int nt,nfirst,nlast,nx,ny,km,ist,ied,jst,jed,mt2_r,sfc_ju,partial,thds;
   double c0,p5,aidif,grav;
   double *dz,*c2dtt,*afact,*psfc,*vdc,*dzt,*rhs,*tnew;
   int *kmt;

   double *p_dz,*p_c2dtt,*p_afact,*p_psfc,*p_vdc,*p_dzt,*p_rhs,*p_tnew;
   int *p_kmt,*p_kif;

   double *p_a,*p_b,*p_e,*p_f,*p_h1;
   int nxy,len,tail,ijst,tlen,tijst,lenmax,nblocks,block,i,k;
   int j,i1,i2,j1,j2,rank,myid,n,mt2,off0,off1,len2,lenk,km4;
   double dtmp,one,ctmp,hfact1,hfactk;
   struct param_vmix s_vitc;
   one=1.0;

    get_reply=0;
    athread_get(0,vitc,&s_vitc,sizeof(s_vitc),(void *)&get_reply,0,0,0);
    while(get_reply!=1);

    nt           = s_vitc.param[0]; 
    nfirst       = s_vitc.param[1]; 
    nlast        = s_vitc.param[2]; 
    nx           = s_vitc.param[3]; 
    ny           = s_vitc.param[4]; 
    km           = s_vitc.param[5]; 
    ist          = s_vitc.param[6]; 
    ied          = s_vitc.param[7]; 
    jst          = s_vitc.param[8];
    jed          = s_vitc.param[9];
    mt2_r        = s_vitc.param[10];
    sfc_ju       = s_vitc.param[11];
    partial      = s_vitc.param[12];
    rank         = s_vitc.param[13];
    thds         = s_vitc.param[14];
  
    c0           = s_vitc.dparam[0]; 
    p5           = s_vitc.dparam[1]; 
    aidif        = s_vitc.dparam[2]; 
    grav         = s_vitc.dparam[3]; 

    dz           = s_vitc.addr[0];
    c2dtt        = s_vitc.addr[1];
    afact        = s_vitc.addr[2];
    psfc         = s_vitc.addr[3];
    vdc          = s_vitc.addr[4];
    dzt          = s_vitc.addr[5];
    rhs          = s_vitc.addr[6];
    tnew         = s_vitc.addr[7];

    kmt          = s_vitc.addr_int[0];
    myid         = athread_get_id(-1);

#ifdef POPcheckJP
    if(rank==0 && myid==0)

    {
      cpe_printf("nx,ny,km,nt,nfirst,nlast,thds %d %d %d %d %d %d %d\n",nx,ny,km,nt,nfirst,nlast, thds);
      cpe_printf("ist ied jst jed partial sfc_ju mt2_r%d %d %d %d %d %d %d\n",ist,ied,jst,jed,partial,sfc_ju,mt2_r);

    }
#endif

       nxy  = nx*ny;
       tail = nxy%thds;
       tlen  = nxy/thds+(tail-myid+thds-1)/thds;  //total-len for myid
       tijst = nxy/thds*myid+min(myid,tail);     //start /myid

     nblocks=(tlen*8*7+tlen*km*4*9)/(60*1024-km*8*4)+1; // assume available ldm=60k
     nblocks=1;
#ifdef SMALL
     nblocks=3;
#endif
#ifdef MED
          nblocks=2;
#endif
#ifdef LARGE
        nblocks=1;
#endif
     lenmax = (tlen+nblocks-1)/nblocks;  //for malloc
     km4=(km+3)/4*4;
     len2=(lenmax+lenmax+3)/4*4;
     lenk=(lenmax*km+3)/4*4;
     km4=(km+3)/4*4;
     lenmax=(lenmax+3)/4*4;
     dsize=km4*3+lenmax*6+lenk*4;

     while(dsize>5800)
     {
     nblocks++;
     lenmax = (tlen+nblocks-1)/nblocks;  //for malloc
     km4=(km+3)/4*4;
     len2=(lenmax+lenmax+3)/4*4;
     lenk=(lenmax*km+3)/4*4;
     km4=(km+3)/4*4;
     lenmax=(lenmax+3)/4*4;
     dsize=km4*3+lenmax*6+lenk*4;
     }

     p_dz       = (double*)ldm_malloc(s_double*km4);   // dimension(km)
     p_c2dtt    = (double*)ldm_malloc(s_double*km4);   // dimension(km)
     p_afact    = (double*)ldm_malloc(s_double*km4);
     p_psfc     = (double*)ldm_malloc(s_double*lenmax);  // (nx_block,ny_block)
     p_vdc      = (double*)ldm_malloc(s_double*lenmax);  // (nx,ny,0:km+1,2,bid)
     p_dzt      = (double*)ldm_malloc(s_double*lenk);  // (nx,ny,0:km+1,bid)   k,k+1
     p_tnew     = (double*)ldm_malloc(s_double*lenk);  // (nx,ny,km,nt)
     p_rhs      = (double*)ldm_malloc(s_double*lenmax);  // (nx,ny,nt)

     p_h1       = (double*)ldm_malloc(s_double*lenmax);
     p_a        = (double*)ldm_malloc(s_double*lenmax);
     p_b        = (double*)ldm_malloc(s_double*lenmax);
//     p_c        = (double*)ldm_malloc(s_double*lenmax);
//     p_d        = (double*)ldm_malloc(s_double*lenmax);
     p_e        = (double*)ldm_malloc(s_double*lenk); //i,j,k               ~32KB
     p_f        = (double*)ldm_malloc(s_double*lenk); //i,j,k  76x52/64x8*62~32KB
     p_kif      = (int *)ldm_malloc(s_int*lenk);
     p_kmt      = (int *)ldm_malloc(s_int*lenmax);     // (nx,ny,bid)

       get_reply =0;
       athread_get(0,dz,p_dz,s_double*km,(void *)&get_reply,0,0,0);
       athread_get(0,c2dtt,p_c2dtt,s_double*km,(void *)&get_reply,0,0,0);
       athread_get(0,afact,p_afact,s_double*km,(void *)&get_reply,0,0,0);
       while (get_reply!=3);

 for(block=0;block<nblocks;block++)
  { //block <
         tail=tlen%nblocks;
         len=tlen/nblocks+(tail-block+nblocks-1)/nblocks;
         ijst=tlen/nblocks*block+min(block,tail)+tijst;
                   
       get_reply  =0;
       get_reply0 =0;
       athread_get(0,psfc+ijst,p_psfc,s_double*len,(void *)&get_reply,0,0,0);
       athread_get(0,kmt+ijst,p_kmt,s_int*len,(void *)&get_reply0,0,0,0);

       hfact1 = *(p_dz)/ *(p_c2dtt);
       while (get_reply!=1);

//       if(l_type==l_vart) {
       if(sfc_ju==1) {
                      for(i=0;i<len;i++){
                       *(p_h1+i) = hfact1 + *(p_psfc+i)/(grav* *(p_c2dtt));
                                        }
                     }
       else          {
                      for(i=0;i<len;i++){
                       *(p_h1+i) = hfact1;
                                        }
                     } 

       while (get_reply0!=1);
           for(i=0;i<len;i++)
            {
              for (k=0;k< (p_kmt[i]-1);k++) { *(p_kif+i+k*len) = 1;}
              for (k=(p_kmt[i]-1);k<km;k++) { *(p_kif+i+k*len) = 0;}
            }

       for(n=(nfirst-1);n<nlast;n++)
       { //n loop left
         mt2 = min (n+1, mt2_r)-1;
         off0=ijst+nxy+nxy*(km+2)*mt2;   //:,:,0,1,bid -> :,:,1,mt2,bid(fortran index)
         off1=ijst+nxy*km*n;             //:,:,1,bid   -> :,:,1,n_c,bid
         get_reply=0; 
         get_reply1=0; 
         get_reply0=0; 
         athread_get(0,vdc+off0,p_vdc,s_double*len,(void *)&get_reply,0,0,0);
         athread_get(0,tnew+off1,p_tnew,s_double*len,(void *)&get_reply1,0,0,0);
         if(n==(nfirst-1)) athread_get(0,dzt+2*nxy+ijst,p_dzt,s_double*len,(void *)&get_reply0,0,0,0);
         athread_get(0,rhs+n*nxy+ijst,p_rhs,s_double*len,(void *)&get_reply,0,0,0);
//         athread_get(0,dzt+ijst,p_dzt,s_double*len*km,&get_reply1,0,(nxy-len)*8,len*8); 
                            //get dzt(ijst:ijed,2:km+1,bid)

          while(get_reply!=2);
//         if(n==(nfirst-1)) { while(get_reply!=3);}
//                        else { while(get_reply!=2);}

         for(i=0;i<len;i++)
         { //i <
           *(p_a+i) = *(p_afact) * *(p_vdc+i);           
           dtmp = one/(*(p_h1+i) + *(p_a+i));
           *(p_e+i) = *(p_a+i) *dtmp;
           *(p_b+i) = *(p_h1+i) * *(p_e+i);
           *(p_f+i) = hfact1 * *(p_rhs+i)*dtmp; 
         } //i >
        
        if(n==(nfirst-1)) { while(get_reply0!=1);}

       if(partial==1)
         { //if par <
         for(k=1;k<km;k++)
           {   //k <
            get_reply=0;
            athread_get(0,vdc+off0+k*nxy,p_vdc,s_double*len,(void *)&get_reply,0,0,0);
            athread_get(0,tnew+off1+k*nxy,p_tnew+k*len,s_double*len,(void *)&get_reply1,0,0,0);
            if(n==(nfirst-1)) athread_get(0,dzt+ijst+(k+2)*nxy,p_dzt+k*len,s_double*len,(void *)&get_reply,0,0,0); 

         if(n==(nfirst-1)) { while(get_reply!=2);}
                        else { while(get_reply!=1);}

             for(i=0;i<len;i++)
              {
                ctmp = *(p_a+i);
                *(p_a+i) = aidif* *(p_vdc+i)/(p5 * (*(p_dzt+i+(k-1)*len)+ *(p_dzt+i+k*len)));
                hfactk   = *(p_dzt+i+(k-1)*len)/ *(p_c2dtt+k);        
                dtmp   = one/(hfactk + *(p_kif+i+len*k)* *(p_a+i) + *(p_b+i));
                *(p_e+k*len+i)= *(p_a+i) *dtmp;
                *(p_b+i)   = (hfactk + *(p_b+i))* *(p_e+k*len+i);
                *(p_f+k*len+i)= ctmp* *(p_f+(k-1)*len+i)*dtmp;
              }
           } // k >
             for(i=0;i<len;i++)
              {
                for(k=p_kmt[i];k<km;k++) { *(p_f+k*len+i)=c0;}
              }
         } //if par >
       else
         { //if par else <
         for(k=1;k<km;k++)
           {   //k <
            get_reply=0;
//            get_reply1=0;
            athread_get(0,vdc+off0+k*nxy,p_vdc,s_double*len,(void *)&get_reply,0,0,0);
            athread_get(0,tnew+off1+k*nxy,p_tnew+k*len,s_double*len,(void *)&get_reply1,0,0,0);
//            if(n==(nfirst-1)) athread_get(0,dzt+ijst+(k+2)*nxy,p_dzt+k*len,s_double*len,(void *)&get_reply,0,0,0); 

//         if (n==(nfirst-1)) {while(get_reply!=2);}
//                else { while(get_reply!=1);}       
             while(get_reply!=1);       

             hfactk = *(p_dz+k)/ *(p_c2dtt+k);
             for(i=0;i<len;i++)
               {
                ctmp = *(p_a+i);
                *(p_a+i) = *(p_afact+k)* *(p_vdc+i);
                dtmp =one/( hfactk + *(p_kif+i+k*len) * *(p_a+i)+ *(p_b+i));
                *(p_e+k*len+i) = *(p_a+i)*dtmp;
                *(p_b+i) = (hfactk + *(p_b+i))* *(p_e+k*len+i);
                *(p_f+k*len+i) = ctmp* *(p_f+(k-1)*len+i)*dtmp;
               }
           }   //k >
             for(i=0;i<len;i++)
               {
                for(k=p_kmt[i];k<km;k++) { *(p_f+k*len+i)=c0;}
               }
         } //if par else >

      
//         while(get_reply1!=1);

         put_reply=0;
//          athread_put(0,p_f+km*len-len,ff+km*nxy-nxy+n*nxy*km+ijst,s_double*len,&put_reply,0,0);

//           for(i=0;i<len;i++)
//            {
//                for(k=(p_kmt[i]-2);k>-1;k--)
//                {
//                 *(p_f+k*len+i)+=  *(p_e+k*len+i)* *(p_f+k*len+len+i);
//                }
//            }
            
           for(k=km-2;k>-1;k--)
             {
                for(i=0;i<len;i++)
                 {
                   *(p_f+k*len+i)+=  *(p_e+k*len+i)* *(p_f+k*len+len+i)*p_kif[k*len+i];
                 }
             } 

         while(get_reply1!=km);
          for(k=0;k<km;k++)
          {       // only adapt to len<=nx
              j1=ijst/nx;
              j2=(ijst+len-1)/nx;
              if(j1==j2)
                {
                   if( j1 >= jst-1 && j1<jed)
                     {
                        i1=j1*nx+ist-1-ijst;i2=j1*nx+ied-ijst;
                        for(i=max(i1,0);i<min(i2,len);i++)
                        {
                          p_tnew[k*len+i]+=p_f[k*len+i]; 
                        } 

                     }
                }
              else   // j2=j1+1
                {
                   if( j1 >= jst-1  )
                     {
                       if(j2 <= jed)
                         { 
                           i1=j1*nx+ist-1-ijst;i2=j1*nx+ied-ijst;
                           for(i=max(i1,0);i<min(i2,len);i++) p_tnew[k*len+i]+=p_f[k*len+i]; 
                         } 
                       if(j2<jed)
                         {
                           i1=j2*nx+ist-1-ijst;i2=j2*nx+ied-ijst;
                           for(i=max(i1,0);i<min(i2,len);i++) p_tnew[k*len+i]+=p_f[k*len+i]; 
                         }
		     } 
                   else
                     {
                       if( j1=jst-2)
                         {
                           i1=j2*nx+ist-1-ijst;i2=j2*nx+ied-ijst;
                           for(i=max(i1,0);i<min(i2,len);i++) p_tnew[k*len+i]+=p_f[k*len+i]; 
                         }
                     }

                }

//           j1=ijst/nx;j2=(ijst+len-1)/nx;
//           for(j=max(j1,jst-1);j<min(j2+1,jed);j++)
//             {
//                i1=j*nx+ist-1-ijst;i2=j*nx+ied-ijst;
//                  for(i=max(i1,0);i<min(i2,len);i++) p_tnew[k*len+i]=p_told[k*len+i]+p_f[k*len+i]; 
//             }

//            for(i=0;i<len;i++)
//            { 
//              io=(ijst+i)%nx; jo=(ijst+i)/nx;
//              if(io>=(ist-1) && io<ied && jo>=(jst-1) && jo<jed) p_tnew[k*len+i]=p_told[k*len+i]+p_f[k*len+i];
//            }
            athread_put(0,p_tnew+k*len,tnew+k*nxy+n*nxy*km+ijst,s_double*len,(void *)&put_reply,0,0);
          }  //k>
         
//       for(k=0;k<km;k++)
//         {
//          for(i=0;i<len;i++)
//           {
//             p_tnew[i]=p_told[i]+ *(p_f+k*len+i);
//           }
//         }
//         athread_put(0,p_f,ff+off1,s_double*len*km,&put_reply,(nxy-len)*8,len*8);
//         while (put_reply !=1);

//        for(k=0;k<km;k++)
//         {
//          athread_put(0,p_f+k*len,ff+k*nxy+n*nxy*km+ijst,s_double*len,&put_reply,0,0);
//         }
        while(put_reply !=km);

       } //n loop right

  } //blcok >

     ldm_free(p_dz,s_double*km4);   // dimension(km)
     ldm_free(p_c2dtt,s_double*km4);   // dimension(km)
     ldm_free(p_psfc,s_double*lenmax);   // dimension(km)
     ldm_free(p_vdc,s_double*lenmax);   
     ldm_free(p_dzt,s_double*lenk);   
     ldm_free(p_tnew,s_double*lenk);   
     ldm_free(p_rhs,s_double*lenmax);   
     ldm_free(p_kmt,s_int*lenmax);   
     ldm_free(p_afact,s_double*km4);   

// temp local
     ldm_free(p_h1,s_double*lenmax);   
     ldm_free(p_a,s_double*lenmax);   
     ldm_free(p_b,s_double*lenmax);   
//     ldm_free(p_c,s_double*lenmax);   
//     ldm_free(p_d,s_double*lenmax);   
     ldm_free(p_e,s_double*lenk);   
     ldm_free(p_f,s_double*lenk);   
     ldm_free(p_kif,s_int*lenk);   

}





