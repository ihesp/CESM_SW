#include <stdio.h>
#include "slave.h"
#include "hmix_del4_struct.h"
#include "cpe_print.h"
#define s_double sizeof(double)
#define s_int sizeof(int)

#define max(x,y) ( x>y?x:y )
#define min(x,y) ( x<y?x:y )


void s_advuhdiffu(struct param_hdiffu *hdiffu1)
{    
    int volatile get_reply,put_reply,get_reply1,get_reply2,get_reply3;
    int myid,nx,ny,k,thds,partial,lhmixu,ldiag;
    double hdiffu_d[5680];
    double c0,c16,am;
    int jlen,jst,jed,i,j,tail,jlen1,jlen2,j1,j2,jt1,jt2;
    int ib,ie,jb,je;

    double *s_duc,*s_dum,*s_dzu;
    double *duc,*dum,*dzu;
 
    double *dun,*dus,*due,*duw,*dmc,*dmn,*dms,*dme,*dmw; 
    double *s_dun,*s_dus,*s_due,*s_duw,*s_dmc,*s_dmn,*s_dms,*s_dme,*s_dmw; 
    double *amf, *s_amf,*umixk,*vmixk,*s_umixk,*s_vmixk;

    double *hduk,*hdvk,*s_hduk,*s_hdvk;

    double *d2uk,*d2vk,*s_d2uk,*s_d2vk;
    double *s_cc,*s_cn,*s_cs,*s_ce,*s_cw;
    double one,sdzu,ddzu;

    int    *kmu;
    int    s_kmu[240];

    struct param_hdiffu s_hdiffu1;
/*zyu  added   */

    int  phid,km,tavgyes,rank,jbeg,ibeg,jend,iend;
    double p25,p5,p125,p5c2dzk,dzrk,dz2rk;

    double *uuu,*vvv,*dyu,*dxu,*uvewns,*kxyu,*uue,*vun  ;

    double *s_uuu,*s_vvv,*s_dxu,*s_dyu,*s_uuw,*s_uue,*s_vus,*s_vun,*s_utmp,*s_vtmp;
    double *wuk,*wukb,*uarear,*s_wukb,*s_wuk,*s_uarear;
    double *luk,*lvk,*s_luk,*s_lvk,*s_uuub,*s_vvvb;
    double tmp1,*s_uuua,*s_vvva;

    double *kxu,*kyu,*s_kxu,*s_kyu;
    
 
//zyu
//   struct param_hdiffu  s_advu1;
   
/*   end      */


    one=1.;
    get_reply=0;
    athread_get(PE_MODE,hdiffu1,&s_hdiffu1,sizeof(s_hdiffu1),
                                             (void*)&get_reply,0,0,0);
    asm volatile("memb");
    while(get_reply!=1);

    phid = athread_get_id(-1);
       
    nx     =s_hdiffu1.param[0];
    ny     =s_hdiffu1.param[1];
    k      =s_hdiffu1.param[2];
    km     =s_hdiffu1.param[3];
    ib     =s_hdiffu1.param[4];
    ie     =s_hdiffu1.param[5];
    jb     =s_hdiffu1.param[6];
    je     =s_hdiffu1.param[7];
    partial=s_hdiffu1.param[8];
    tavgyes=s_hdiffu1.param[9];
    lhmixu =s_hdiffu1.param[10];
    rank   =s_hdiffu1.param[11];
    thds   =s_hdiffu1.param[12]/2;

    ibeg =ib;
    iend =ie;
    jbeg =jb;
    jend =je;


if(phid < thds)
{
    myid=phid;
    c0   =s_hdiffu1.dparam[0];
    c16  =s_hdiffu1.dparam[1];
    am   =s_hdiffu1.dparam[2];

    dzu  =s_hdiffu1.addr[4];
    duc  =s_hdiffu1.addr[12];
    dum  =s_hdiffu1.addr[13];
    dun  =s_hdiffu1.addr[14];
    dus  =s_hdiffu1.addr[15];
    due  =s_hdiffu1.addr[16];
    duw  =s_hdiffu1.addr[17];
    dmc  =s_hdiffu1.addr[18];
    dmn  =s_hdiffu1.addr[19];
    dms  =s_hdiffu1.addr[20];
    dme  =s_hdiffu1.addr[21];
    dmw  =s_hdiffu1.addr[22];
    umixk=s_hdiffu1.addr[23];
    vmixk=s_hdiffu1.addr[24];
//    d2uk =s_hdiffu1.addr[14];
//    d2vk =s_hdiffu1.addr[15];
    hduk =s_hdiffu1.addr[25];
    hdvk =s_hdiffu1.addr[26];
    amf  =s_hdiffu1.addr[27];
    

    kmu  =s_hdiffu1.addr_int[0];

        tail=ny%thds;
        jlen=ny/thds+(tail-myid+thds-1)/thds;
        jst =ny/thds*myid+min(myid,tail);
        jed =jst+jlen-1;

if(jlen <=0) return;

    s_cc    = hdiffu_d;
    s_duc   = hdiffu_d+  nx*jlen + 2*nx;
    s_dum   = hdiffu_d+2*nx*jlen + 4*nx;
    s_dun   = hdiffu_d+3*nx*jlen + 6*nx;
    s_dus   = hdiffu_d+4*nx*jlen + 8*nx;
    s_due   = hdiffu_d+5*nx*jlen +10*nx;
    s_duw   = hdiffu_d+6*nx*jlen +12*nx;
    s_dzu   = hdiffu_d+7*nx*jlen +14*nx;
    s_umixk = hdiffu_d+8*nx*jlen +18*nx;
    s_vmixk = hdiffu_d+9*nx*jlen +22*nx;
    s_d2uk  = hdiffu_d+10*nx*jlen +26*nx;
    s_d2vk  = hdiffu_d+11*nx*jlen +28*nx;
    s_cn    = hdiffu_d+12*nx*jlen +30*nx;
    s_cs    = hdiffu_d+13*nx*jlen +32*nx;
    s_ce    = hdiffu_d+14*nx*jlen +34*nx;
    s_cw    = hdiffu_d+15*nx*jlen +36*nx;
    s_dmn   = hdiffu_d+16*nx*jlen +38*nx;
    s_dms   = hdiffu_d+17*nx*jlen +40*nx;
    s_dme   = hdiffu_d+18*nx*jlen +42*nx;
    s_dmw   = hdiffu_d+19*nx*jlen +44*nx;
    s_dmc   = hdiffu_d+20*nx*jlen +46*nx;
    s_hduk  = hdiffu_d+21*nx*jlen +48*nx;
    s_hdvk  = hdiffu_d+22*nx*jlen +48*nx;
    s_amf   = hdiffu_d+23*nx*jlen +48*nx;
    

//   [jst,min[jed,jb-3]   c0      
//   [max(jst,jb-2),min(jed,je)]  for d loop 
//   [max(jst,je+1),jed] c0
//
//   [jst,min[jed,jb-2]   c0      
//   [max(jst,jb-1),min(jed,je-1)]  for H 
//   [max(jst,je),jed] c0
    put_reply=0;
    if(jst >= je || jed < (jb-1))
      {
        for(i=0;i<jlen*nx;i++)
          {
            s_hduk[i]=c0;
            s_hdvk[i]=c0;
          }
      } // j field1 >

    else if(jst >=(jb-1) && jed < je )
      { // j field 2 <
        get_reply=0;
        get_reply1=0;
        jt1=(jst-1)*nx;
        athread_get(0,duc+jt1,s_duc,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
        athread_get(0,dum+jt1,s_dum,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
        athread_get(0,dzu+jt1-nx,  s_dzu,(jlen+4)*nx*s_double,  (void *)&get_reply1,0,0,0);
        athread_get(0,umixk+jt1-nx,s_umixk,(jlen+4)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,vmixk+jt1-nx,s_vmixk,(jlen+4)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,dun+jt1,s_dun,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,dus+jt1,s_dus,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,due+jt1,s_due,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,duw+jt1,s_duw,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,amf+jt1,s_amf,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,kmu+jt1,s_kmu,(jlen+2)*nx*s_int,   (void *)&get_reply1,0,0,0);
        athread_get(0,dmc+jt1,s_dmc,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,dmn+jt1,s_dmn,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,dms+jt1,s_dms,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,dme+jt1,s_dme,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,dmw+jt1,s_dmw,(jlen+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        while(get_reply!=2);
        for(i=0;i<(jlen+2)*nx;i++)
          {
            s_cc[i]=s_duc[i]+s_dum[i];
          }

        while(get_reply1!=14);
        if(partial==1)
          { //partial <
            for(j=0;j<(jlen+2);j++)  //[jst-1,jed+1]
              { //j <
                for(i=ib-2;i<ie+1;i++)
                  {
                    jt2=j*nx+i;
                    if(k>s_kmu[jt2])
                      {s_d2uk[jt2]=c0;
                       s_d2vk[jt2]=c0;
                      }
                    else
                      {
                        sdzu=s_dzu[jt2+nx];
                        ddzu=one/sdzu;
                        s_cn[jt2]=s_dun[jt2]*min(s_dzu[jt2+nx*2],sdzu)*ddzu;
                        s_cs[jt2]=s_dus[jt2]*min(s_dzu[jt2     ],sdzu)*ddzu;
                        s_ce[jt2]=s_due[jt2]*min(s_dzu[jt2+nx+1],sdzu)*ddzu;
                        s_cw[jt2]=s_duw[jt2]*min(s_dzu[jt2+nx-1],sdzu)*ddzu;
                        s_d2uk[jt2]=(s_cc[jt2]*s_umixk[jt2+nx]
                                   +s_cn[jt2]*s_umixk[jt2+nx+nx]
                                   +s_cs[jt2]*s_umixk[jt2]
                                   +s_ce[jt2]*s_umixk[jt2+nx+1]
                                   +s_cw[jt2]*s_umixk[jt2+nx-1])
                                  +(s_dmc[jt2]*s_vmixk[jt2+nx]
                                   +s_dmn[jt2]*s_vmixk[jt2+nx+nx]
                                   +s_dms[jt2]*s_vmixk[jt2]
                                   +s_dme[jt2]*s_vmixk[jt2+nx+1]
                                   +s_dmw[jt2]*s_vmixk[jt2+nx-1]);
                        s_d2vk[jt2]=(s_cc[jt2]*s_vmixk[jt2+nx]
                                   +s_cn[jt2]*s_vmixk[jt2+nx+nx]
                                   +s_cs[jt2]*s_vmixk[jt2]
                                   +s_ce[jt2]*s_vmixk[jt2+nx+1]
                                   +s_cw[jt2]*s_vmixk[jt2+nx-1])
                                   +(s_dmc[jt2]*s_umixk[jt2+nx]
                                   +s_dmn[jt2]*s_umixk[jt2+nx+nx]
                                   +s_dms[jt2]*s_umixk[jt2]
                                   +s_dme[jt2]*s_umixk[jt2+nx+1]
                                   +s_dmw[jt2]*s_umixk[jt2+nx-1]);
                       if(lhmixu==1)
                         {
                           s_d2uk[jt2]*=s_amf[jt2];
                           s_d2vk[jt2]*=s_amf[jt2];
                         }
                     } //k else >
                  } //i >
                for(i=0;i<ib-2;i++)
                  {
                     s_d2uk[j*nx+i]=c0;
                     s_d2vk[j*nx+i]=c0;
                  } //i>
                for(i=ie+1;i<nx;i++)
                  {
                     s_d2uk[j*nx+i]=c0;
                     s_d2vk[j*nx+i]=c0;
                  } //i>
           } //j loop


       for(j=0;j<jlen;j++)  //[jst,jed]
         {
                for(i=ib-1;i<ie;i++)
                  {
                    jt2=j*nx+i;
                    if(k<=s_kmu[jt2+nx])
                      {
                         s_hduk[jt2]=am*((s_cc[jt2+nx]*s_d2uk[jt2+nx]
                                       +s_cn[jt2+nx]*s_d2uk[jt2+nx+nx]
                                       +s_cs[jt2+nx]*s_d2uk[jt2]
                                       +s_ce[jt2+nx]*s_d2uk[jt2+nx+1]
                                       +s_cw[jt2+nx]*s_d2uk[jt2+nx-1])
                                      +(s_dmc[jt2+nx]*s_d2vk[jt2+nx]
                                       +s_dmn[jt2+nx]*s_d2vk[jt2+nx+nx]
                                       +s_dms[jt2+nx]*s_d2vk[jt2]
                                       +s_dme[jt2+nx]*s_d2vk[jt2+nx+1]
                                       +s_dmw[jt2+nx]*s_d2vk[jt2+nx-1]));
                         s_hdvk[jt2]=am*((s_cc[jt2+nx]*s_d2vk[jt2+nx]
                                       +s_cn[jt2+nx]*s_d2vk[jt2+nx+nx]
                                       +s_cs[jt2+nx]*s_d2vk[jt2]
                                       +s_ce[jt2+nx]*s_d2vk[jt2+nx+1]
                                       +s_cw[jt2+nx]*s_d2vk[jt2+nx-1])
                                      +(s_dmc[jt2+nx]*s_d2uk[jt2+nx]
                                       +s_dmn[jt2+nx]*s_d2uk[jt2+nx+nx]
                                       +s_dms[jt2+nx]*s_d2uk[jt2]
                                       +s_dme[jt2+nx]*s_d2uk[jt2+nx+1]
                                       +s_dmw[jt2+nx]*s_d2uk[jt2+nx-1]));
                      }
                    else
                      {
                        s_hduk[jt2]=c0; 
                        s_hdvk[jt2]=c0;
                      }
                  } //i>
                for(i=0;i<ib-1;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
                for(i=ie;i<nx;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
         } // j >
      } //partial if>
    else 
      { //partial else <
         for(j=0;j<(jlen+2);j++)  //[jst-1,jed+1]
           {
              for(i=ib-2;i<ie+1;i++)
                {
                  if(k>s_kmu[j*nx+i]) 
                    {
                      s_d2uk[j*nx+i]=c0;
                      s_d2vk[j*nx+i]=c0;
                    }
                  else
                    {
                  s_d2uk[j*nx+i]=(s_cc[j*nx+i]*s_umixk[j*nx+nx+i]
                                 +s_dun[j*nx+i]*s_umixk[(j+2)*nx+i]
                                 +s_dus[j*nx+i]*s_umixk[j*nx+i]
                                 +s_due[j*nx+i]*s_umixk[j*nx+nx+i+1]
                                 +s_duw[j*nx+i]*s_umixk[j*nx+nx+i-1])
                                 +(s_dmc[j*nx+i]*s_vmixk[j*nx+nx+i]
                                 +s_dmn[j*nx+i]*s_vmixk[(j+2)*nx+i]
                                 +s_dms[j*nx+i]*s_vmixk[j*nx+i]
                                 +s_dme[j*nx+i]*s_vmixk[j*nx+nx+i+1]
                                 +s_dmw[j*nx+i]*s_vmixk[j*nx+nx+i-1]);
                  s_d2vk[j*nx+i]=(s_cc[j*nx+i]*s_vmixk[j*nx+nx+i]
                                 +s_dun[j*nx+i]*s_vmixk[(j+2)*nx+i]
                                 +s_dus[j*nx+i]*s_vmixk[j*nx+i]
                                 +s_due[j*nx+i]*s_vmixk[j*nx+nx+i+1]
                                 +s_duw[j*nx+i]*s_vmixk[j*nx+nx+i-1])
                                 +(s_dmc[j*nx+i]*s_umixk[j*nx+nx+i]
                                 +s_dmn[j*nx+i]*s_umixk[(j+2)*nx+i]
                                 +s_dms[j*nx+i]*s_umixk[j*nx+i]
                                 +s_dme[j*nx+i]*s_umixk[j*nx+nx+i+1]
                                 +s_dmw[j*nx+i]*s_umixk[j*nx+nx+i-1]);
                       if(lhmixu==1)
                         {
                           s_d2uk[j*nx+i]=s_d2uk[j*nx+i]*s_amf[j*nx+i];
                           s_d2vk[j*nx+i]=s_d2vk[j*nx+i]*s_amf[j*nx+i];

                         }
                     }

              } //i >
            for(i=0;i<ib-2;i++)
              {
                 s_d2uk[j*nx+i]=c0;
                 s_d2vk[j*nx+i]=c0;
              } //i>
            for(i=ie+1;i<nx;i++)
              {
                 s_d2uk[j*nx+i]=c0;
                 s_d2vk[j*nx+i]=c0;
              } //i>
           } //j loop


       for(j=0;j<jlen;j++)  //[jst,jed]
         { // j <
                for(i=ib-1;i<ie;i++)
                  {
                    if(k<=s_kmu[i+j*nx+nx])
                      {
                         s_hduk[i+j*nx]=am*(( s_cc[i+j*nx+nx]*s_d2uk[j*nx+nx+i]
                                       +s_dun[i+j*nx+nx]*s_d2uk[(j+2)*nx+i]
                                       +s_dus[i+j*nx+nx]*s_d2uk[j*nx+i]
                                       +s_due[i+j*nx+nx]*s_d2uk[j*nx+nx+i+1]
                                       +s_duw[i+j*nx+nx]*s_d2uk[j*nx+nx+i-1])
                                      +(s_dmc[i+j*nx+nx]*s_d2vk[j*nx+nx+i]
                                       +s_dmn[i+j*nx+nx]*s_d2vk[(j+2)*nx+i]
                                       +s_dms[i+j*nx+nx]*s_d2vk[j*nx+i]
                                       +s_dme[i+j*nx+nx]*s_d2vk[j*nx+nx+i+1]
                                       +s_dmw[i+j*nx+nx]*s_d2vk[j*nx+nx+i-1]));
                        s_hdvk[i+j*nx]=am*(( s_cc[i+j*nx+nx]*s_d2vk[j*nx+nx+i]
                                       +s_dun[i+j*nx+nx]*s_d2vk[(j+2)*nx+i]
                                       +s_dus[i+j*nx+nx]*s_d2vk[j*nx+i]
                                       +s_due[i+j*nx+nx]*s_d2vk[j*nx+nx+i+1]
                                       +s_duw[i+j*nx+nx]*s_d2vk[j*nx+nx+i-1])
                                      +(s_dmc[i+j*nx+nx]*s_d2uk[j*nx+nx+i]
                                       +s_dmn[i+j*nx+nx]*s_d2uk[(j+2)*nx+i]
                                       +s_dms[i+j*nx+nx]*s_d2uk[j*nx+i]
                                       +s_dme[i+j*nx+nx]*s_d2uk[j*nx+nx+i+1]
                                       +s_dmw[i+j*nx+nx]*s_d2uk[j*nx+nx+i-1]));
                     }
                   else
                     {
                       s_hduk[i+j*nx]=c0; s_hdvk[i+j*nx]=c0;
                     }
                  }//i>
                for(i=0;i<ib-1;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
                for(i=ie;i<nx;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
          } //j >
        } //partial if else>
      } //j field2
//
    else if(jed >=(jb-1) && jst < (jb-1))
      { // j field3
        // [jst, jb-1)     [jb-1,jed]
        // the first part [jst,jb-1) h=c0;
         jlen1=jb-1-jst;
         for(j=0;j<jb-1-jst;j++)
           {
             for(i=0;i<nx;i++)
              {
                 hduk[j*nx+i]=c0;
                 hdvk[j*nx+i]=c0;
              }
           } //j >
         // the second part [jb-1,jed] for h, and [jb-2,jed+1 for d2]
         jlen2=jlen-jlen1;
         get_reply=0;
         get_reply1=0;
         athread_get(0,duc+(jst+jlen1-1)*nx,s_duc,(jlen2+2)*nx*s_double,(void *)&get_reply,0,0,0);
         athread_get(0,dum+(jst+jlen1-1)*nx,s_dum,(jlen2+2)*nx*s_double,(void *)&get_reply,0,0,0);
         athread_get(0,dzu+(jst+jlen1-2)*nx,s_dzu,(jlen2+4)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,umixk+(jst+jlen1-2)*nx,s_umixk,(jlen2+4)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,vmixk+(jst+jlen1-2)*nx,s_vmixk,(jlen2+4)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dun+(jst+jlen1-1)*nx,s_dun,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dus+(jst+jlen1-1)*nx,s_dus,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,due+(jst+jlen1-1)*nx,s_due,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,duw+(jst+jlen1-1)*nx,s_duw,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,amf+(jst+jlen1-1)*nx,s_amf,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,kmu+(jst+jlen1-1)*nx,s_kmu,(jlen2+2)*nx*s_int,   (void *)&get_reply1,0,0,0);
         athread_get(0,dmc+(jst+jlen1-1)*nx,s_dmc,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dmn+(jst+jlen1-1)*nx,s_dmn,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dms+(jst+jlen1-1)*nx,s_dms,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dme+(jst+jlen1-1)*nx,s_dme,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dmw+(jst+jlen1-1)*nx,s_dmw,(jlen2+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         while(get_reply!=2);
         for(i=0;i<(jlen2+2)*nx;i++)
           {
             s_cc[i+(jlen1-1)*nx]=s_duc[i+(jlen1-1)*nx]+s_dum[i+(jlen1-1)*nx];
           }
         while(get_reply1!=14);
 
        if(partial==1)
          { //partial <
            for(j=0;j<(jlen2+2);j++)  //[jb-1,jed+1]
              { //j <
                for(i=ib-2;i<ie+1;i++)
                  { //i <
                    if(k>s_kmu[(j+jlen1-1)*nx+i])
                      {
                        s_d2uk[(j+jlen1-1)*nx+i]=c0;
                        s_d2vk[(j+jlen1-1)*nx+i]=c0;
                      }
                    else
                      { //k else <
                    sdzu=s_dzu[(jlen1+j)*nx+i];
                    ddzu=one/sdzu;
                    s_cn[(j+jlen1-1)*nx+i]=s_dun[(j+jlen1-1)*nx+i]*min(s_dzu[(j+jlen1+1)*nx+i],sdzu)*ddzu;
                    s_cs[(j+jlen1-1)*nx+i]=s_dus[(j+jlen1-1)*nx+i]*min(s_dzu[(j+jlen1-1)*nx+i],sdzu)*ddzu;
                    s_ce[(j+jlen1-1)*nx+i]=s_due[(j+jlen1-1)*nx+i]*min(s_dzu[(j+jlen1  )*nx+i+1],sdzu)*ddzu;
                    s_cw[(j+jlen1-1)*nx+i]=s_duw[(j+jlen1-1)*nx+i]*min(s_dzu[(j+jlen1  )*nx+i-1],sdzu)*ddzu;
                    s_d2uk[(j+jlen1-1)*nx+i]=(s_cc[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1)*nx+i]
                                   +s_cn[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1+1)*nx+i]
                                   +s_cs[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1-1)*nx+i]
                                   +s_ce[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1  )*nx+i+1]
                                   +s_cw[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1  )*nx+i-1])
                                  +(s_dmc[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1)*nx+i]
                                   +s_dmn[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1+1)*nx+i]
                                   +s_dms[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1-1)*nx+i]
                                   +s_dme[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1  )*nx+i+1]
                                   +s_dmw[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1  )*nx+i-1]);
                    s_d2vk[j*nx+i]=(s_cc[j*nx+i]*s_vmixk[j*nx+nx+i]
                                   +s_cn[j*nx+i]*s_vmixk[(j+2)*nx+i]
                                   +s_cs[j*nx+i]*s_vmixk[j*nx+i]
                                   +s_ce[j*nx+i]*s_vmixk[j*nx+nx+i+1]
                                   +s_cw[j*nx+i]*s_vmixk[j*nx+nx+i-1])
                                   +(s_dmc[j*nx+i]*s_umixk[j*nx+nx+i]
                                   +s_dmn[j*nx+i]*s_umixk[(j+2)*nx+i]
                                   +s_dms[j*nx+i]*s_umixk[j*nx+i]
                                   +s_dme[j*nx+i]*s_umixk[j*nx+nx+i+1]
                                   +s_dmw[j*nx+i]*s_umixk[j*nx+nx+i-1]);
                         if(lhmixu==1)
                           {
                              s_d2uk[(j+jlen1-1)*nx+i]*=s_amf[(j+jlen1-1)*nx+i];
                              s_d2vk[(j+jlen1-1)*nx+i]*=s_amf[(j+jlen1-1)*nx+i];
                           }
                       } // k else>
                    } //i >
              } //j <

       for(j=0;j<jlen2;j++)  //[jb-1,jed]
         {
                for(i=ib-1;i<ie;i++)
                  {
                    if(k<=s_kmu[i+j*nx+nx])
                      {
                         s_hduk[i+(j+jlen1)*nx]=am*((s_cc[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i]
                                       +s_cn[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1+1)*nx+i]
                                       +s_cs[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1-1)*nx+i]
                                       +s_ce[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i+1]
                                       +s_cw[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i-1])
                                      +(s_dmc[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i]
                                       +s_dmn[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1+1)*nx+i]
                                       +s_dms[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1-1)*nx+i]
                                       +s_dme[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i+1]
                                       +s_dmw[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i-1]));
                         s_hdvk[i+j*nx]=am*((s_cc[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i]
                                       +s_cn[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1+1)*nx+i]
                                       +s_cs[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1-1)*nx+i]
                                       +s_ce[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i+1]
                                       +s_cw[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i-1])
                                      +(s_dmc[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i]
                                       +s_dmn[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1+1)*nx+i]
                                       +s_dms[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1-1)*nx+i]
                                       +s_dme[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i+1]
                                       +s_dmw[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i-1]));
                      }
                    else
                      {
                        s_hduk[i+j*nx]=c0; 
                        s_hdvk[i+j*nx]=c0;
                      }
                  } //i>
                for(i=0;i<ib-1;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
                for(i=ie;i<nx;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
         } // j >


          } //partial >
          { //partial else<
            for(j=0;j<(jlen2+2);j++)  //[jb-1,jed+1]
              { //j <
                for(i=ib-2;i<ie+1;i++)
                  { //i <
                    if(k>s_kmu[(j+jlen1-1)*nx+i])
                      {
                        s_d2uk[(j+jlen1-1)*nx+i]=c0;
                        s_d2vk[(j+jlen1-1)*nx+i]=c0;
                      }
                    else
                      { //k else <
                    s_d2uk[(j+jlen1-1)*nx+i]=(s_cc[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1)*nx+i]
                                   +s_dun[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1+1)*nx+i]
                                   +s_dus[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1-1)*nx+i]
                                   +s_due[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1  )*nx+i+1]
                                   +s_duw[(j+jlen1-1)*nx+i]*s_umixk[(j+jlen1  )*nx+i-1])
                                  +(s_dmc[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1)*nx+i]
                                   +s_dmn[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1+1)*nx+i]
                                   +s_dms[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1-1)*nx+i]
                                   +s_dme[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1  )*nx+i+1]
                                   +s_dmw[(j+jlen1-1)*nx+i]*s_vmixk[(j+jlen1  )*nx+i-1]);
                    s_d2vk[j*nx+i]=(s_cc[j*nx+i]*s_vmixk[j*nx+nx+i]
                                   +s_dun[j*nx+i]*s_vmixk[(j+2)*nx+i]
                                   +s_dus[j*nx+i]*s_vmixk[j*nx+i]
                                   +s_due[j*nx+i]*s_vmixk[j*nx+nx+i+1]
                                   +s_duw[j*nx+i]*s_vmixk[j*nx+nx+i-1])
                                   +(s_dmc[j*nx+i]*s_umixk[j*nx+nx+i]
                                   +s_dmn[j*nx+i]*s_umixk[(j+2)*nx+i]
                                   +s_dms[j*nx+i]*s_umixk[j*nx+i]
                                   +s_dme[j*nx+i]*s_umixk[j*nx+nx+i+1]
                                   +s_dmw[j*nx+i]*s_umixk[j*nx+nx+i-1]);
                         if(lhmixu==1)
                           {
                              s_d2uk[(j+jlen1-1)*nx+i]*=s_amf[(j+jlen1-1)*nx+i];
                              s_d2vk[(j+jlen1-1)*nx+i]*=s_amf[(j+jlen1-1)*nx+i];
                           }
                       } // k else>
                    } //i >
              } //j <

       for(j=0;j<jlen2;j++)  //[jb-1,jed]
         {
                for(i=ib-1;i<ie;i++)
                  {
                    if(k<=s_kmu[i+j*nx+nx])
                      {
                         s_hduk[i+(j+jlen1)*nx]=am*((s_cc[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i]
                                       +s_dun[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1+1)*nx+i]
                                       +s_dus[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1-1)*nx+i]
                                       +s_due[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i+1]
                                       +s_duw[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i-1])
                                      +(s_dmc[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i]
                                       +s_dmn[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1+1)*nx+i]
                                       +s_dms[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1-1)*nx+i]
                                       +s_dme[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i+1]
                                       +s_dmw[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i-1]));
                         s_hdvk[i+j*nx]=am*((s_cc[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i]
                                       +s_dun[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1+1)*nx+i]
                                       +s_dus[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1-1)*nx+i]
                                       +s_due[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i+1]
                                       +s_duw[i+(j+jlen1)*nx]*s_d2vk[(j+jlen1)*nx+i-1])
                                      +(s_dmc[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i]
                                       +s_dmn[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1+1)*nx+i]
                                       +s_dms[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1-1)*nx+i]
                                       +s_dme[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i+1]
                                       +s_dmw[i+(j+jlen1)*nx]*s_d2uk[(j+jlen1)*nx+i-1]));
                      }
                    else
                      {
                        s_hduk[i+j*nx]=c0; 
                        s_hdvk[i+j*nx]=c0;
                      }
                  } //i>
                for(i=0;i<ib-1;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
                for(i=ie;i<nx;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
         } // j >


          } //partial else>
      }//j field3
//
    else if(jst <=(je-1) && jed > (je-1) )
      { //j field4  [jst,je-1] [je,jed]
         jlen1=je-jst;
         get_reply=0;
         get_reply1=0;
         athread_get(0,duc+(jst-1)*nx,s_duc,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
         athread_get(0,dum+(jst-1)*nx,s_dum,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
         athread_get(0,dzu+(jst-2)*nx,s_dzu,(jlen1+4)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,umixk+(jst-2)*nx,s_umixk,(jlen1+4)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,vmixk+(jst-2)*nx,s_vmixk,(jlen1+4)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dun+(jst-1)*nx,s_dun,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dus+(jst-1)*nx,s_dus,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,due+(jst-1)*nx,s_due,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,duw+(jst-1)*nx,s_duw,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,amf+(jst-1)*nx,s_amf,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,kmu+(jst-1)*nx,s_kmu,(jlen1+2)*nx*s_int,   (void *)&get_reply1,0,0,0);
         athread_get(0,dmc+(jst-1)*nx,s_dmc,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dmn+(jst-1)*nx,s_dmn,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dms+(jst-1)*nx,s_dms,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dme+(jst-1)*nx,s_dme,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
         athread_get(0,dmw+(jst-1)*nx,s_dmw,(jlen1+2)*nx*s_double,(void *)&get_reply1,0,0,0);
        while(get_reply!=2);
        for(i=0;i<(jlen1+2)*nx;i++)
          {
            s_cc[i]=s_duc[i]+s_dum[i];
          }
//    first part [jst,je-1] for h, [jst-1,je] for d2
        while(get_reply1!=14);
        if(partial==1)
          { //partial <
            for(j=0;j<(jlen1+2);j++)  //[jst-1,je]
              { //j <
                for(i=ib-2;i<ie+1;i++)
                  {
                    if(k>s_kmu[i+j*nx])
                      {
                         s_d2uk[j*nx+i]=c0;
                         s_d2vk[j*nx+i]=c0;
                      }
                    else
                      {
                    sdzu=s_dzu[j*nx+nx+i];
                    ddzu=one/sdzu;
                    s_cn[j*nx+i]=s_dun[j*nx+i]*min(s_dzu[(j+2)*nx+i],sdzu)*ddzu;
                    s_cs[j*nx+i]=s_dus[j*nx+i]*min(s_dzu[j*nx+i],sdzu)*ddzu;
                    s_ce[j*nx+i]=s_due[j*nx+i]*min(s_dzu[(j+1)*nx+i+1],sdzu)*ddzu;
                    s_cw[j*nx+i]=s_duw[j*nx+i]*min(s_dzu[(j+1)*nx+i-1],sdzu)*ddzu;
                    s_d2uk[j*nx+i]=(s_cc[j*nx+i]*s_umixk[j*nx+nx+i]
                                   +s_cn[j*nx+i]*s_umixk[(j+2)*nx+i]
                                   +s_cs[j*nx+i]*s_umixk[j*nx+i]
                                   +s_ce[j*nx+i]*s_umixk[j*nx+nx+i+1]
                                   +s_cw[j*nx+i]*s_umixk[j*nx+nx+i-1])
                                  +(s_dmc[j*nx+i]*s_vmixk[j*nx+nx+i]
                                   +s_dmn[j*nx+i]*s_vmixk[(j+2)*nx+i]
                                   +s_dms[j*nx+i]*s_vmixk[j*nx+i]
                                   +s_dme[j*nx+i]*s_vmixk[j*nx+nx+i+1]
                                   +s_dmw[j*nx+i]*s_vmixk[j*nx+nx+i-1]);
                    s_d2vk[j*nx+i]=(s_cc[j*nx+i]*s_vmixk[j*nx+nx+i]
                                   +s_cn[j*nx+i]*s_vmixk[(j+2)*nx+i]
                                   +s_cs[j*nx+i]*s_vmixk[j*nx+i]
                                   +s_ce[j*nx+i]*s_vmixk[j*nx+nx+i+1]
                                   +s_cw[j*nx+i]*s_vmixk[j*nx+nx+i-1])
                                   +(s_dmc[j*nx+i]*s_umixk[j*nx+nx+i]
                                   +s_dmn[j*nx+i]*s_umixk[(j+2)*nx+i]
                                   +s_dms[j*nx+i]*s_umixk[j*nx+i]
                                   +s_dme[j*nx+i]*s_umixk[j*nx+nx+i+1]
                                   +s_dmw[j*nx+i]*s_umixk[j*nx+nx+i-1]);
                         if(lhmixu==1)
                           {
                             s_d2uk[j*nx+i]*=s_amf[j*nx+i];
                             s_d2vk[j*nx+i]*=s_amf[j*nx+i];
                           }
                       }
                  } //i >
                for(i=0;i<ib-2;i++)
                  {
                     s_d2uk[j*nx+i]=c0;
                     s_d2vk[j*nx+i]=c0;
                  } //i>
                for(i=ie+1;i<nx;i++)
                  {
                     s_d2uk[j*nx+i]=c0;
                     s_d2vk[j*nx+i]=c0;
                  } //i>
           } //j loop

       jlen1=je-jst;
       for(j=0;j<jlen1;j++)  //[jst,je)
         {
                for(i=ib-1;i<ie;i++)
                  {
                    if(k<=s_kmu[i+j*nx+nx])
                      {
                         s_hduk[i+j*nx]=am*((s_cc[i+j*nx+nx]*s_d2uk[j*nx+nx+i]
                                       +s_cn[i+j*nx+nx]*s_d2uk[(j+2)*nx+i]
                                       +s_cs[i+j*nx+nx]*s_d2uk[j*nx+i]
                                       +s_ce[i+j*nx+nx]*s_d2uk[j*nx+nx+i+1]
                                       +s_cw[i+j*nx+nx]*s_d2uk[j*nx+nx+i-1])
                                      +(s_dmc[i+j*nx+nx]*s_d2vk[j*nx+nx+i]
                                       +s_dmn[i+j*nx+nx]*s_d2vk[(j+2)*nx+i]
                                       +s_dms[i+j*nx+nx]*s_d2vk[j*nx+i]
                                       +s_dme[i+j*nx+nx]*s_d2vk[j*nx+nx+i+1]
                                       +s_dmw[i+j*nx+nx]*s_d2vk[j*nx+nx+i-1]));
                         s_hdvk[i+j*nx]=am*((s_cc[i+j*nx+nx]*s_d2vk[j*nx+nx+i]
                                       +s_cn[i+j*nx+nx]*s_d2vk[(j+2)*nx+i]
                                       +s_cs[i+j*nx+nx]*s_d2vk[j*nx+i]
                                       +s_ce[i+j*nx+nx]*s_d2vk[j*nx+nx+i+1]
                                       +s_cw[i+j*nx+nx]*s_d2vk[j*nx+nx+i-1])
                                      +(s_dmc[i+j*nx+nx]*s_d2uk[j*nx+nx+i]
                                       +s_dmn[i+j*nx+nx]*s_d2uk[(j+2)*nx+i]
                                       +s_dms[i+j*nx+nx]*s_d2uk[j*nx+i]
                                       +s_dme[i+j*nx+nx]*s_d2uk[j*nx+nx+i+1]
                                       +s_dmw[i+j*nx+nx]*s_d2uk[j*nx+nx+i-1]));
                      }
                    else
                      {
                        s_hduk[i+j*nx]=c0; 
                        s_hdvk[i+j*nx]=c0;
                      }
                  } //i>
                for(i=0;i<ib-1;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
                for(i=ie;i<nx;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
         } // j >

        } //partial if>
    else 
        { // partial else<
         for(j=0;j<(jlen1+2);j++)  //[jst-1,jed+1]
           {
              for(i=ib-2;i<ie+1;i++)
                {
                  if(k>s_kmu[j*nx+i])
                    {
                      s_d2uk[j*nx+i]=c0;
                      s_d2vk[j*nx+i]=c0;
                    }
                  else
                    {
                  s_d2uk[j*nx+i]=(s_cc[j*nx+i]*s_umixk[j*nx+nx+i]
                                 +s_dun[j*nx+i]*s_umixk[(j+2)*nx+i]
                                 +s_dus[j*nx+i]*s_umixk[j*nx+i]
                                 +s_due[j*nx+i]*s_umixk[j*nx+nx+i+1]
                                 +s_duw[j*nx+i]*s_umixk[j*nx+nx+i-1])
                                 +(s_dmc[j*nx+i]*s_vmixk[j*nx+nx+i]
                                 +s_dmn[j*nx+i]*s_vmixk[(j+2)*nx+i]
                                 +s_dms[j*nx+i]*s_vmixk[j*nx+i]
                                 +s_dme[j*nx+i]*s_vmixk[j*nx+nx+i+1]
                                 +s_dmw[j*nx+i]*s_vmixk[j*nx+nx+i-1]);
                  s_d2vk[j*nx+i]=(s_cc[j*nx+i]*s_vmixk[j*nx+nx+i]
                                 +s_dun[j*nx+i]*s_vmixk[(j+2)*nx+i]
                                 +s_dus[j*nx+i]*s_vmixk[j*nx+i]
                                 +s_due[j*nx+i]*s_vmixk[j*nx+nx+i+1]
                                 +s_duw[j*nx+i]*s_vmixk[j*nx+nx+i-1])
                                 +(s_dmc[j*nx+i]*s_umixk[j*nx+nx+i]
                                 +s_dmn[j*nx+i]*s_umixk[(j+2)*nx+i]
                                 +s_dms[j*nx+i]*s_umixk[j*nx+i]
                                 +s_dme[j*nx+i]*s_umixk[j*nx+nx+i+1]
                                 +s_dmw[j*nx+i]*s_umixk[j*nx+nx+i-1]);
                        if(lhmixu==1)
                          {
                            s_d2uk[j*nx+i]*=s_amf[j*nx+i];
                            s_d2vk[j*nx+i]*=s_amf[j*nx+i];
                          }
                     }
              } //i >
            for(i=0;i<ib-2;i++)
              {
                 s_d2uk[j*nx+i]=c0;
                 s_d2vk[j*nx+i]=c0;
              } //i>
            for(i=ie+1;i<nx;i++)
              {
                 s_d2uk[j*nx+i]=c0;
                 s_d2vk[j*nx+i]=c0;
              } //i>
            } //j loop

       jlen1=je-jst;
       for(j=0;j<jlen1;j++)  //[jst,je)
         {
                for(i=ib-1;i<ie;i++)
                  {
                    if(k<=s_kmu[i+j*nx+nx])
                      {
                         s_hduk[i+j*nx]=am*((s_cc[i+j*nx+nx]*s_d2uk[j*nx+nx+i]
                                       +s_dun[i+j*nx+nx]*s_d2uk[(j+2)*nx+i]
                                       +s_dus[i+j*nx+nx]*s_d2uk[j*nx+i]
                                       +s_due[i+j*nx+nx]*s_d2uk[j*nx+nx+i+1]
                                       +s_duw[i+j*nx+nx]*s_d2uk[j*nx+nx+i-1])
                                      +(s_dmc[i+j*nx+nx]*s_d2vk[j*nx+nx+i]
                                       +s_dmn[i+j*nx+nx]*s_d2vk[(j+2)*nx+i]
                                       +s_dms[i+j*nx+nx]*s_d2vk[j*nx+i]
                                       +s_dme[i+j*nx+nx]*s_d2vk[j*nx+nx+i+1]
                                       +s_dmw[i+j*nx+nx]*s_d2vk[j*nx+nx+i-1]));
                         s_hdvk[i+j*nx]=am*((s_cc[i+j*nx+nx]*s_d2vk[j*nx+nx+i]
                                       +s_dun[i+j*nx+nx]*s_d2vk[(j+2)*nx+i]
                                       +s_dus[i+j*nx+nx]*s_d2vk[j*nx+i]
                                       +s_due[i+j*nx+nx]*s_d2vk[j*nx+nx+i+1]
                                       +s_duw[i+j*nx+nx]*s_d2vk[j*nx+nx+i-1])
                                      +(s_dmc[i+j*nx+nx]*s_d2uk[j*nx+nx+i]
                                       +s_dmn[i+j*nx+nx]*s_d2uk[(j+2)*nx+i]
                                       +s_dms[i+j*nx+nx]*s_d2uk[j*nx+i]
                                       +s_dme[i+j*nx+nx]*s_d2uk[j*nx+nx+i+1]
                                       +s_dmw[i+j*nx+nx]*s_d2uk[j*nx+nx+i-1]));
                      }
                    else
                      {
                        s_hduk[i+j*nx]=c0; 
                        s_hdvk[i+j*nx]=c0;
                      }
                  } //i>
                for(i=0;i<ib-1;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
                for(i=ie;i<nx;i++)
                  {
                    s_hduk[i+j*nx]=c0;
                    s_hdvk[i+j*nx]=c0;
                  }
         } // j >
         } // partial else>
//  part2 of field4
     for(j=je-jst;j<jlen;j++)
       {
          s_hduk[j*nx+i]=c0;
          s_hdvk[j*nx+i]=c0;

       }
     } //j field4  [jst,je-1] [je,jed]

   put_reply=0;
   athread_put(0,s_hduk,hduk+jst*nx,nx*jlen*s_double,(void *)&put_reply,0,0);
   athread_put(0,s_hdvk,hdvk+jst*nx,nx*jlen*s_double,(void *)&put_reply,0,0);
   while(put_reply!=2);
} // if h>
else
{
   myid=phid-thds;

/*    int volatile get_reply,put_reply,get_reply1,get_reply2,get_reply3;
    int myid,jbeg,jend,ibeg,iend,nx,ny,k,km,thds,my_task,partial,tavgyes;
    double advu_d[3320];
    double c0,p25,p125,p5,p5c2dzk,dzrk,dz2rk;
    int jlen,jst,jed,i,j,tail,jlen1,j1,j2;
    double *uuu,*vvv,*dyu,*dxu,*dzu,*uuw,*uue,*vus,*vun;
    double *s_uuu,*s_vvv,*s_dxu,*s_dyu,*s_dzu,*s_uuw,*s_uue,*s_vus,*s_vun,*s_utmp,*s_vtmp;
    double *wuk,*wukb,*uarear,*s_wukb,*s_wuk,*s_uarear;
    double *luk,*lvk,*s_luk,*s_lvk,*s_cc,*s_uuub,*s_vvvb;
    double tmp1,*s_uuua,*s_vvva;

    double *kxu,*kyu,*s_kxu,*s_kyu;
    int    *kmu;
    int    s_kmu[160];

    struct param_advu s_advu1;
    get_reply=0;
    athread_get(PE_MODE,advu1,&s_advu1,sizeof(s_advu1),
                                             (void*)&get_reply,0,0,0);
    asm volatile("memb");
    while(get_reply!=1);

    myid = athread_get_id(-1);
*/       

    c0     =s_hdiffu1.dparam[0];
    p25    =s_hdiffu1.dparam[4];
    p125   =s_hdiffu1.dparam[5];
    p5     =s_hdiffu1.dparam[3];
    p5c2dzk=p5*s_hdiffu1.dparam[6];
    dzrk   =s_hdiffu1.dparam[7];
    dz2rk  =s_hdiffu1.dparam[8];

    uuu    =s_hdiffu1.addr[0];
    vvv    =s_hdiffu1.addr[1];
    dxu    =s_hdiffu1.addr[2];
    dyu    =s_hdiffu1.addr[3];
    dzu    =s_hdiffu1.addr[4];
    uvewns =s_hdiffu1.addr[5];
    wukb   =s_hdiffu1.addr[6];
    wuk    =s_hdiffu1.addr[7];
    uarear =s_hdiffu1.addr[8];
    luk    =s_hdiffu1.addr[9];
    lvk    =s_hdiffu1.addr[10];
//zyu
    kxyu   =s_hdiffu1.addr[11];

    kmu    =s_hdiffu1.addr_int[0];

        tail=ny%thds;
        jlen=ny/thds+(tail-myid+thds-1)/thds;
        jst =ny/thds*myid+min(myid,tail);
        jed =jst+jlen-1;

if(jlen <=0) return;

    s_uuw = hdiffu_d;
    s_uue = hdiffu_d+  nx*jlen;
    s_vus = hdiffu_d+2*nx*jlen;
    s_vun = hdiffu_d+3*nx*jlen+nx;
    s_uuu = hdiffu_d+4*nx*jlen+nx;
    s_vvv = hdiffu_d+5*nx*jlen+ 3*nx;
    s_dxu = hdiffu_d+6*nx*jlen+ 5*nx;
    s_dyu = hdiffu_d+7*nx*jlen+ 7*nx;
    s_dzu = hdiffu_d+8*nx*jlen+ 9*nx;
    s_utmp= hdiffu_d+9*nx*jlen+11*nx;
    s_vtmp= hdiffu_d+10*nx*jlen+13*nx;
    s_wukb= hdiffu_d+11*nx*jlen+15*nx;
    s_wuk = hdiffu_d+12*nx*jlen+15*nx;
  s_uarear= hdiffu_d+13*nx*jlen+15*nx;
    s_luk = hdiffu_d+14*nx*jlen+15*nx;
    s_lvk = hdiffu_d+15*nx*jlen+15*nx;
    s_cc  = hdiffu_d+16*nx*jlen+15*nx;
    s_uuub= hdiffu_d+17*nx*jlen+15*nx;
    s_vvvb= hdiffu_d+18*nx*jlen+15*nx;
    s_uuua= hdiffu_d+19*nx*jlen+15*nx;
    s_vvva= hdiffu_d+20*nx*jlen+15*nx;
    s_kxu = hdiffu_d+21*nx*jlen+15*nx;
    s_kyu = hdiffu_d+22*nx*jlen+15*nx;
   

//   [max(0,jst),min[jed,jbeg-3]       c0
//   [max(jst,jbeg-2),min(jed,jend)]   -->utmp vtmp
//   [max(jst,jend+1),min(jed,ny-1)]   c0
    get_reply=0;
    get_reply1=0;
    get_reply2=0;
    get_reply3=0;
    put_reply=0;
    athread_get(0,uuu+(jst-1)*nx,s_uuu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
    athread_get(0,vvv+(jst-1)*nx,s_vvv,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
    athread_get(0,dxu+(jst-1)*nx,s_dxu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
    athread_get(0,dyu+(jst-1)*nx,s_dyu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);

    athread_get(0,wuk+jst*nx,s_wuk,jlen*nx*s_double,(void *)&get_reply1,0,0,0);
    athread_get(0,uarear+jst*nx,s_uarear,jlen*nx*s_double,(void *)&get_reply1,0,0,0);
    if(k>1)
      {
         athread_get(0,uuu+jst*nx-nx*ny,s_uuub,jlen*nx*s_double,(void *)&get_reply2,0,0,0);
         athread_get(0,vvv+jst*nx-nx*ny,s_vvvb,jlen*nx*s_double,(void *)&get_reply2,0,0,0);
      }
    if(k<km)
      {
         athread_get(0,uuu+jst*nx+nx*ny,s_uuua,jlen*nx*s_double,(void *)&get_reply3,0,0,0);
         athread_get(0,vvv+jst*nx+nx*ny,s_vvva,jlen*nx*s_double,(void *)&get_reply3,0,0,0);
      }

  if(partial==1)
   {
      athread_get(0,dzu+(jst-1)*nx,s_dzu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
    if(jed < jbeg-2 || jst > jend) // do j=jbeg-1,jend+1 in F90, [jbeg-2,jend] in c   
      {
        for(i=0;i<jlen*nx;i++) 
          {
            s_uuw[i]=c0;
            s_uue[i]=c0;
            s_vus[i]=c0;
            s_vun[i]=c0;
          }
        while(get_reply!=5);
      } // j field 1>

    else if(jst < jbeg-2 && jed >=jbeg-2)
      {
        jlen1=jbeg-2-jst; // [jst,jbeg-3] 
        for(i=0;i<jlen1*nx;i++)
          {
            s_uuw[i]=c0;
            s_uue[i]=c0;
            s_vus[i]=c0;
            s_vun[i]=c0;
          }
//        athread_get(0,uuu+(jbeg-3)*nx,s_uuu,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,vvv+(jbeg-3)*nx,s_vvv,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dxu+(jbeg-3)*nx,s_dxu,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dyu+(jbeg-3)*nx,s_dyu,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dzu+(jbeg-3)*nx,s_dzu,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
        while(get_reply!=5);
        for(i=0;i<(jlen-jlen1+2)*nx;i++)
          {
            s_utmp[i]=s_uuu[i+jlen1*nx]*s_dyu[i+jlen1*nx]*s_dzu[i+jlen1*nx];
            s_vtmp[i]=s_vvv[i+jlen1*nx]*s_dxu[i+jlen1*nx]*s_dzu[i+jlen1*nx];
          }
        for(j=0;j<jlen-jlen1;j++)
          {
            for(i=ibeg-2;i<iend+1;i++)
              {
                s_uuw[i+jlen1*nx+j*nx]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i-1]) //i,j  i-1,j
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i-1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i-1]);
                               //i,j-1           i-1,j-1        i,j+1             i-1,j+1            
                s_uue[i+jlen1*nx+j*nx]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i+1])  //i+1,j i,j
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i+1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i+1]);
                    //            i,j-1                  i+1,j-1        i,j+1              i+1,j+1
                s_vus[i+jlen1*nx+j*nx]=p25*(s_vtmp[j*nx+nx+i]+s_vtmp[j*nx+i]) //i,j  i,j-1
                    +p125*(s_vtmp[j*nx+i-1]+s_vtmp[(j+1)*nx+i-1]+s_vtmp[j*nx+i+1]+s_vtmp[(j+1)*nx+i+1]);
                    //            i-1,j-1             i-1,j          i+1,j-1             i+1,j
                s_vun[i+jlen1*nx+j*nx]=p25*(s_vtmp[(j+1)*nx+i]+s_vtmp[(j+2)*nx+i]) //i,j   i,j+1
                    +p125*(s_vtmp[(j+1)*nx+i-1]+s_vtmp[(j+2)*nx+i-1]+s_vtmp[(j+1)*nx+i+1]+s_vtmp[(j+2)*nx+i+1]); 
                    //         i-1,j               i-1,j+1                     i+1,j             i+1,j+1    
              }
            for(i=0;i<ibeg-2;i++) 
              {
                s_uuw[i+jlen1*nx+j*nx]=c0;
                s_uue[i+jlen1*nx+j*nx]=c0;
                s_vus[i+jlen1*nx+j*nx]=c0;
                s_vun[i+jlen1*nx+j*nx]=c0;
              }
            for(i=iend+1;i<nx;i++)
              {
                s_uuw[i+jlen1*nx+j*nx]=c0;
                s_uue[i+jlen1*nx+j*nx]=c0;
                s_vus[i+jlen1*nx+j*nx]=c0;
                s_vun[i+jlen1*nx+j*nx]=c0;
              }
          } // j loop>
          for(i=ibeg-2;i<iend+1;i++)
            {
              s_vus[jlen*nx+i]=p25*(s_vtmp[jlen*nx+nx+i]+s_vtmp[jlen*nx+i])
                              +p125*(s_vtmp[jlen*nx+i-1]+s_vtmp[jlen*nx+nx+i-1]+
                                     s_vtmp[jlen*nx+i+1]+s_vtmp[jlen*nx+nx+i+1]);
            }
          for(i=0;i<ibeg-2;i++)
            {
              s_vus[jlen*nx+i]=c0;
            }
          for(i=iend+1;i<nx;i++)
            {
              s_vus[jlen*nx+i]=c0;
            }
      }  //  j field 2 >
    else if(jst >= jbeg-2 && jed <=jend)
      {
//        athread_get(0,uuu+(jst-1)*nx,s_uuu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,vvv+(jst-1)*nx,s_vvv,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dxu+(jst-1)*nx,s_dxu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dyu+(jst-1)*nx,s_dyu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dzu+(jst-1)*nx,s_dzu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
        while(get_reply!=5);
        for(i=0;i<(jlen+2)*nx;i++)
          {
            s_utmp[i]=s_uuu[i]*s_dyu[i]*s_dzu[i];
            s_vtmp[i]=s_vvv[i]*s_dxu[i]*s_dzu[i];
          }

        for(j=0;j<jlen;j++)
          {
            for(i=ibeg-2;i<iend+1;i++)
              {
                s_uuw[j*nx+i]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i-1]) //i,j  i-1,j
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i-1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i-1]);              
                     //               i,j-1             i-1,j-1          i,j+1                i-1,j+1
                s_uue[j*nx+i]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i+1])
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i+1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i+1]);
                s_vus[j*nx+i]=p25*(s_vtmp[j*nx+nx+i]+s_vtmp[j*nx+i])
                    +p125*(s_vtmp[j*nx+i-1]+s_vtmp[(j+1)*nx+i-1]+s_vtmp[j*nx+i+1]+s_vtmp[(j+1)*nx+i+1]);
                s_vun[j*nx+i]=p25*(s_vtmp[(j+1)*nx+i]+s_vtmp[(j+2)*nx+i])
                    +p125*(s_vtmp[(j+1)*nx+i-1]+s_vtmp[(j+2)*nx+i-1]+s_vtmp[(j+1)*nx+i+1]+s_vtmp[(j+2)*nx+i+1]); 
              }
            for(i=0;i<ibeg-2;i++) 
              {
                s_uuw[j*nx+i]=c0;
                s_uue[j*nx+i]=c0;
                s_vus[j*nx+i]=c0;
                s_vun[j*nx+i]=c0;
              }
            for(i=iend+1;i<nx;i++)
              {
                s_uuw[j*nx+i]=c0;
                s_uue[j*nx+i]=c0;
                s_vus[j*nx+i]=c0;
                s_vun[j*nx+i]=c0;
              }
          } //j
          for(i=ibeg-2;i<iend+1;i++) //vus(j+1) by d2dk=f(  vus(j+1))
            {
              s_vus[jlen*nx+i]=p25*(s_vtmp[jlen*nx+nx+i]+s_vtmp[jlen*nx+i])
                              +p125*(s_vtmp[jlen*nx+i-1]+s_vtmp[jlen*nx+nx+i-1]+
                                     s_vtmp[jlen*nx+i+1]+s_vtmp[jlen*nx+nx+i+1]);
            }
          for(i=0;i<ibeg-2;i++) s_vus[jlen*nx+i]=c0;
          for(i=iend+1;i<nx;i++) s_vus[jlen*nx+i]=c0;

      } // j field 3>
    else if(jst <= jend && jed > jend)
      {
        jlen1=jend-jst+1;
//        athread_get(0,uuu+(jst-1)*nx,s_uuu,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,vvv+(jst-1)*nx,s_vvv,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dxu+(jst-1)*nx,s_dxu,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dyu+(jst-1)*nx,s_dyu,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dzu+(jst-1)*nx,s_dzu,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
        while(get_reply!=5);
        for(i=0;i<(jlen1+2)*nx;i++)
          {
            s_utmp[i]=s_uuu[i]*s_dyu[i]*s_dzu[i];
            s_vtmp[i]=s_vvv[i]*s_dxu[i]*s_dzu[i];
          }
        for(j=0;j<jlen1;j++)
          {
            for(i=ibeg-2;i<iend+1;i++)
              {
                s_uuw[i+j*nx]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i-1])
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i-1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i-1]);
                s_uue[i+j*nx]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i+1])
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i+1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i+1]);
                s_vus[i+j*nx]=p25*(s_vtmp[j*nx+nx+i]+s_vtmp[j*nx+i])
                    +p125*(s_vtmp[j*nx+i-1]+s_vtmp[(j+1)*nx+i-1]+s_vtmp[j*nx+i+1]+s_vtmp[(j+1)*nx+i+1]);
                s_vun[i+j*nx]=p25*(s_vtmp[(j+1)*nx+i]+s_vtmp[(j+2)*nx+i])
                    +p125*(s_vtmp[(j+1)*nx+i-1]+s_vtmp[(j+2)*nx+i-1]+s_vtmp[(j+1)*nx+i+1]+s_vtmp[(j+2)*nx+i+1]);
              }
            for(i=0;i<ibeg-2;i++)
              {
                s_uuw[i+j*nx]=c0;
                s_uue[i+j*nx]=c0;
                s_vus[i+j*nx]=c0;
                s_vun[i+j*nx]=c0;
              }
            for(i=iend+1;i<nx;i++)
              {
                s_uuw[i+j*nx]=c0;
                s_uue[i+j*nx]=c0;
                s_vus[i+j*nx]=c0;
                s_vun[i+j*nx]=c0;
              }
          }
        for(i=0;i<(jlen-jlen1)*nx;i++)
          {
                s_uuw[i+jlen1*nx]=c0;
                s_uue[i+jlen1*nx]=c0;
                s_vus[i+jlen1*nx]=c0;
                s_vun[i+jlen1*nx]=c0;
          }
      } // j field4 >
      if(tavgyes==1)
        {
        athread_put(0,s_uue,uue+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
        athread_put(0,s_vun,vun+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
        }


        while(get_reply1!=2);
        for(i=0;i<jlen*nx;i++)
         {
           s_wukb[i]=s_wuk[i]+(s_vun[i]-s_vus[i]+s_uue[i]-s_uuw[i])*s_uarear[i];
         }


        get_reply1=0;
        athread_get(0,kxyu+jst*nx,      s_kxu,jlen*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,kxyu+jst*nx+nx*ny,s_kyu,jlen*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,kmu+jst*nx,s_kmu,jlen*nx*s_int,   (void *)&get_reply1,0,0,0);



     if(tavgyes==1)
      {
        athread_put(0,s_wukb,wukb+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
      }
     else
      {
        athread_put(0,s_wukb,wuk+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
      }

//   luk lvk
    if(jed < jbeg-1 || jst > jend-1) // do j=jbeg,jend in F90, [jbeg-1,jend-1] in c   
      {
        for(i=0;i<jlen*nx;i++)
          {
            s_luk[i]=c0; s_lvk[i]=c0;
          }
      }//j field1>
    else if(jst <jbeg-1 && jed >=jbeg-1)
      {
        jlen1=jbeg-2-jst+1; //[jst,jbeg-2]
        for(i=0;i<jlen1*nx;i++)
          {
            s_luk[i]=c0; s_lvk[i]=c0;
          }
        for(j=0;j<jlen-jlen1;j++)
         {
           for(i=ibeg-1;i<iend;i++)
             {
              s_cc[i+jlen1*nx+j*nx]=s_vus[i+jlen1*nx+j*nx+nx]-s_vus[i+jlen1*nx+j*nx]
                                    +s_uuw[i+jlen1*nx+j*nx+1] -s_uuw[i+jlen1*nx+j*nx];
              s_luk[i+jlen1*nx+j*nx]=p5*(s_cc[i+jlen1*nx+j*nx]*s_uuu[i+jlen1*nx+j*nx+nx]+
                                         s_vus[i+jlen1*nx+(j+1)*nx]*s_uuu[i+jlen1*nx+(j+2)*nx]
                                        -s_vus[i+jlen1*nx+j*nx]*s_uuu[i+jlen1*nx+j*nx]
                                        +s_uuw[i+1+jlen1*nx+j*nx]*s_uuu[i+1+jlen1*nx+j*nx+nx]
                                        -s_uuw[i+jlen1*nx+j*nx]*s_uuu[i-1+jlen1*nx+j*nx+nx])
                                        *s_uarear[i+jlen1*nx+j*nx]/s_dzu[i+jlen1*nx+j*nx+nx];
              s_lvk[i+jlen1*nx+j*nx]=p5*(s_cc[i+jlen1*nx+j*nx]*s_vvv[i+jlen1*nx+j*nx+nx]+
                                        s_vus[i+jlen1*nx+j*nx+nx]*s_vvv[i+jlen1*nx+(j+2)*nx]
                                       -s_vus[i+jlen1*nx+j*nx]*s_vvv[i+jlen1*nx+j*nx]
                                       +s_uuw[i+1+jlen1*nx+j*nx]*s_vvv[i+1+jlen1*nx+j*nx+nx]
                                       -s_uuw[i+jlen1*nx+j*nx]*s_vvv[i-1+jlen1*nx+j*nx+nx])
                                       *s_uarear[i+jlen1*nx+j*nx]/s_dzu[i+jlen1*nx+j*nx+nx];


             }
           for(i=0;i<ibeg-1;i++)
             {
               s_luk[i]=c0; s_lvk[i]=0;
             }
           for(i=iend;i<nx;i++)
             {
               s_luk[i]=c0; s_lvk[i]=0;
             }

         }// j> 

      } //j field2 >
    else if(jst >=jbeg-1 && jed < jend) //[jbeg-1,jend-1]
      {
        for(j=0;j<jlen;j++)
          {
            for (i=ibeg-1;i<iend;i++)
              {
                s_cc[j*nx+i]=s_vus[j*nx+nx+i]-s_vus[j*nx+i]+s_uuw[i+1+j*nx]-s_uuw[j*nx+i];
                s_luk[j*nx+i]=p5*(s_cc[i+j*nx]*s_uuu[i+j*nx+nx]
                                 +s_vus[j*nx+nx+i]*s_uuu[i+(j+2)*nx]
                                 -s_vus[j*nx+i]   *s_uuu[i+j*nx]
                                 +s_uuw[i+1+j*nx]*s_uuu[i+1+j*nx+nx]
                                 -s_uuw[i+j*nx]*s_uuu[i-1+j*nx+nx])
                                 *s_uarear[j*nx+i]/s_dzu[i+j*nx+nx];
                s_lvk[j*nx+i]=p5*(s_cc[j*nx+i]*s_vvv[i+j*nx+nx]
                                 +s_vus[i+j*nx+nx]*s_vvv[i+(j+2)*nx]
                                 -s_vus[i+j*nx]*s_vvv[i+j*nx]
                                 +s_uuw[i+1+j*nx]*s_vvv[i+1+j*nx+nx]
                                 -s_uuw[i+j*nx]*s_vvv[i-1+j*nx+nx])
                                 *s_uarear[j*nx+i]/s_dzu[i+j*nx+nx]; 
              } // i>
            for(i=0;i<ibeg-1;i++)
              {
                s_luk[j*nx+i]=c0; s_lvk[j*nx+i]=c0;
              }
            for(i=iend;i<nx;i++)
              {
                s_luk[j*nx+i]=c0; s_lvk[j*nx+i]=c0;
              }
          }// j>
      } //j field 3>
    else if(jst <= jend-1 && jed > jend-1)
      {
         jlen1=jend-jst;
         for(j=0;j<jlen1;j++)          
          {
            for (i=ibeg-1;i<iend;i++)
              {
                s_cc[j*nx+i]=s_vus[j*nx+nx+i]-s_vus[j*nx+i]+s_uuw[i+1+j*nx]-s_uuw[j*nx+i];
                s_luk[j*nx+i]=p5*(s_cc[i+j*nx]*s_uuu[i+j*nx+nx]
                                 +s_vus[j*nx+nx+i]*s_uuu[i+(j+2)*nx]
                                 -s_vus[j*nx+i]   *s_uuu[i+j*nx]
                                 +s_uuw[i+1+j*nx]*s_uuu[i+1+j*nx+nx]
                                 -s_uuw[i+j*nx]*s_uuu[i-1+j*nx+nx])
                                 *s_uarear[j*nx+i]/s_dzu[i+j*nx+nx];
                s_lvk[j*nx+i]=p5*(s_cc[j*nx+i]*s_vvv[i+j*nx+nx]
                                 +s_vus[i+j*nx+nx]*s_vvv[i+(j+2)*nx]
                                 -s_vus[i+j*nx]*s_vvv[i+j*nx]
                                 +s_uuw[i+1+j*nx]*s_vvv[i+1+j*nx+nx]
                                 -s_uuw[i+j*nx]*s_vvv[i-1+j*nx+nx])
                                 *s_uarear[j*nx+i]/s_dzu[i+j*nx+nx];
              } // i>
            for(i=0;i<ibeg-1;i++)
              {
                s_luk[j*nx+i]=c0; s_lvk[j*nx+i]=c0;
              }
            for(i=iend;i<nx;i++)
              {
                s_luk[j*nx+i]=c0; s_lvk[j*nx+i]=c0;
              }
          }// j>
        for(i=0;j<(jlen-jlen1)*nx;i++)
          {
            s_luk[jlen1*nx+i]=c0;
            s_lvk[jlen1*nx+i]=c0;

          }
      } //j field4 >

     if(k==1) 
       {
         for(i=0;i<jlen*nx;i++)
           {
             s_luk[i]+=dzrk*s_wuk[i]*s_uuu[i+nx];
             s_lvk[i]+=dzrk*s_wuk[i]*s_vvv[i+nx];
           }
       }
     else
       {
         while(get_reply2!=2);
         for(i=0;i<jlen*nx;i++)
           {
             tmp1=p5/s_dzu[i+nx]*s_wuk[i];
             s_luk[i]+=tmp1*(s_uuub[i]+s_uuu[i+nx]);
             s_lvk[i]+=tmp1*(s_vvvb[i]+s_vvv[i+nx]);
           }
       }

    if(k<km)
      {
        while(get_reply3!=2);
        for(i=0;i<jlen*nx;i++)
          {
            tmp1=p5/s_dzu[i+nx]*s_wukb[i];
            s_luk[i]-=tmp1*(s_uuu[i+nx]+s_uuua[i]); 
            s_lvk[i]-=tmp1*(s_vvv[i+nx]+s_vvva[i]); 
          }
      }
   

   } // partial>
 else
   { //partial else <
    if(jed < jbeg-2 || jst > jend) // do j=jbeg-1,jend+1 in F90, [jbeg-2,jend] in c   
      {
        for(i=0;i<jlen*nx;i++) 
          {
            s_uuw[i]=c0;
            s_uue[i]=c0;
            s_vus[i]=c0;
            s_vun[i]=c0;
          }
        while(get_reply!=4);
      } // j field 1>

    else if(jst < jbeg-2 && jed >=jbeg-2)
      {
        jlen1=jbeg-2-jst; // [jst,jbeg-3] 
        for(i=0;i<jlen1*nx;i++)
          {
            s_uuw[i]=c0;
            s_uue[i]=c0;
            s_vus[i]=c0;
            s_vun[i]=c0;
          }
//        athread_get(0,uuu+(jbeg-3)*nx,s_uuu,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,vvv+(jbeg-3)*nx,s_vvv,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dxu+(jbeg-3)*nx,s_dxu,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dyu+(jbeg-3)*nx,s_dyu,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dzu+(jbeg-3)*nx,s_dzu,(jlen-jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
        while(get_reply!=4);
        for(i=0;i<(jlen-jlen1+2)*nx;i++)
          {
            s_utmp[i]=s_uuu[i+jlen1*nx]*s_dyu[i+jlen1*nx];
            s_vtmp[i]=s_vvv[i+jlen1*nx]*s_dxu[i+jlen1*nx];
          }
        for(j=0;j<jlen-jlen1;j++)
          {
            for(i=ibeg-2;i<iend+1;i++)
              {
                s_uuw[i+jlen1*nx+j*nx]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i-1]) //i,j  i-1,j
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i-1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i-1]);
                               //i,j-1           i-1,j-1        i,j+1             i-1,j+1            
                s_uue[i+jlen1*nx+j*nx]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i+1])  //i+1,j i,j
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i+1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i+1]);
                    //            i,j-1                  i+1,j-1        i,j+1              i+1,j+1
                s_vus[i+jlen1*nx+j*nx]=p25*(s_vtmp[j*nx+nx+i]+s_vtmp[j*nx+i]) //i,j  i,j-1
                    +p125*(s_vtmp[j*nx+i-1]+s_vtmp[(j+1)*nx+i-1]+s_vtmp[j*nx+i+1]+s_vtmp[(j+1)*nx+i+1]);
                    //            i-1,j-1             i-1,j          i+1,j-1             i+1,j
                s_vun[i+jlen1*nx+j*nx]=p25*(s_vtmp[(j+1)*nx+i]+s_vtmp[(j+2)*nx+i]) //i,j   i,j+1
                    +p125*(s_vtmp[(j+1)*nx+i-1]+s_vtmp[(j+2)*nx+i-1]+s_vtmp[(j+1)*nx+i+1]+s_vtmp[(j+2)*nx+i+1]); 
                    //         i-1,j               i-1,j+1                     i+1,j             i+1,j+1    
              }
            for(i=0;i<ibeg-2;i++) 
              {
                s_uuw[i+jlen1*nx+j*nx]=c0;
                s_uue[i+jlen1*nx+j*nx]=c0;
                s_vus[i+jlen1*nx+j*nx]=c0;
                s_vun[i+jlen1*nx+j*nx]=c0;
              }
            for(i=iend+1;i<nx;i++)
              {
                s_uuw[i+jlen1*nx+j*nx]=c0;
                s_uue[i+jlen1*nx+j*nx]=c0;
                s_vus[i+jlen1*nx+j*nx]=c0;
                s_vun[i+jlen1*nx+j*nx]=c0;
              }
          } //j loop>
          for(i=ibeg-2;i<iend+1;i++)
            {
              s_vus[jlen*nx+i]=p25*(s_vtmp[jlen*nx+nx+i]+s_vtmp[jlen*nx+i])
                              +p125*(s_vtmp[jlen*nx+i-1]+s_vtmp[jlen*nx+nx+i-1]+
                                     s_vtmp[jlen*nx+i+1]+s_vtmp[jlen*nx+nx+i+1]);
            }
      }  //  j field 2 >
    else if(jst >= jbeg-2 && jed <=jend)
      {
//        athread_get(0,uuu+(jst-1)*nx,s_uuu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,vvv+(jst-1)*nx,s_vvv,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dxu+(jst-1)*nx,s_dxu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dyu+(jst-1)*nx,s_dyu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dzu+(jst-1)*nx,s_dzu,(jlen+2)*nx*s_double,(void *)&get_reply,0,0,0);
        while(get_reply!=4);
        for(i=0;i<(jlen+2)*nx;i++)
          {
            s_utmp[i]=s_uuu[i]*s_dyu[i];
            s_vtmp[i]=s_vvv[i]*s_dxu[i];
          }

        for(j=0;j<jlen;j++)
          {
            for(i=ibeg-2;i<iend+1;i++)
              {
                s_uuw[j*nx+i]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i-1]) //i,j  i-1,j
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i-1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i-1]);              
                     //               i,j-1             i-1,j-1          i,j+1                i-1,j+1
                s_uue[j*nx+i]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i+1])
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i+1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i+1]);
                s_vus[j*nx+i]=p25*(s_vtmp[j*nx+nx+i]+s_vtmp[j*nx+i])
                    +p125*(s_vtmp[j*nx+i-1]+s_vtmp[(j+1)*nx+i-1]+s_vtmp[j*nx+i+1]+s_vtmp[(j+1)*nx+i+1]);
                s_vun[j*nx+i]=p25*(s_vtmp[(j+1)*nx+i]+s_vtmp[(j+2)*nx+i])
                    +p125*(s_vtmp[(j+1)*nx+i-1]+s_vtmp[(j+2)*nx+i-1]+s_vtmp[(j+1)*nx+i+1]+s_vtmp[(j+2)*nx+i+1]); 
              }
            for(i=0;i<ibeg-2;i++) 
              {
                s_uuw[j*nx+i]=c0;
                s_uue[j*nx+i]=c0;
                s_vus[j*nx+i]=c0;
                s_vun[j*nx+i]=c0;
              }
            for(i=iend+1;i<nx;i++)
              {
                s_uuw[j*nx+i]=c0;
                s_uue[j*nx+i]=c0;
                s_vus[j*nx+i]=c0;
                s_vun[j*nx+i]=c0;
              }
          } //j
          for(i=ibeg-2;i<iend+1;i++)
            {
              s_vus[jlen*nx+i]=p25*(s_vtmp[jlen*nx+nx+i]+s_vtmp[jlen*nx+i])
                              +p125*(s_vtmp[jlen*nx+i-1]+s_vtmp[jlen*nx+nx+i-1]+
                                     s_vtmp[jlen*nx+i+1]+s_vtmp[jlen*nx+nx+i+1]);
            }
   
      } // j field 3>
    else if(jst <= jend && jed > jend)
      {
        jlen1=jend-jst+1;
//        athread_get(0,uuu+(jst-1)*nx,s_uuu,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,vvv+(jst-1)*nx,s_vvv,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dxu+(jst-1)*nx,s_dxu,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dyu+(jst-1)*nx,s_dyu,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
//        athread_get(0,dzu+(jst-1)*nx,s_dzu,(jlen1+2)*nx*s_double,(void *)&get_reply,0,0,0);
        while(get_reply!=4);
        for(i=0;i<(jlen1+2)*nx;i++)
          {
            s_utmp[i]=s_uuu[i]*s_dyu[i];
            s_vtmp[i]=s_vvv[i]*s_dxu[i];
          }
        for(j=0;j<jlen1;j++)
          {
            for(i=ibeg-2;i<iend+1;i++)
              {
                s_uuw[i+j*nx]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i-1])
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i-1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i-1]);
                s_uue[i+j*nx]=p25*(s_utmp[j*nx+nx+i]+s_utmp[j*nx+nx+i+1])
                    +p125*(s_utmp[j*nx+i]+s_utmp[j*nx+i+1]+s_utmp[(j+2)*nx+i]+s_utmp[(j+2)*nx+i+1]);
                s_vus[i+j*nx]=p25*(s_vtmp[j*nx+nx+i]+s_vtmp[j*nx+i])
                    +p125*(s_vtmp[j*nx+i-1]+s_vtmp[(j+1)*nx+i-1]+s_vtmp[j*nx+i+1]+s_vtmp[(j+1)*nx+i+1]);
                s_vun[i+j*nx]=p25*(s_vtmp[(j+1)*nx+i]+s_vtmp[(j+2)*nx+i])
                    +p125*(s_vtmp[(j+1)*nx+i-1]+s_vtmp[(j+2)*nx+i-1]+s_vtmp[(j+1)*nx+i+1]+s_vtmp[(j+2)*nx+i+1]);
              }
            for(i=0;i<ibeg-2;i++)
              {
                s_uuw[i+j*nx]=c0;
                s_uue[i+j*nx]=c0;
                s_vus[i+j*nx]=c0;
                s_vun[i+j*nx]=c0;
              }
            for(i=iend+1;i<nx;i++)
              {
                s_uuw[i+j*nx]=c0;
                s_uue[i+j*nx]=c0;
                s_vus[i+j*nx]=c0;
                s_vun[i+j*nx]=c0;
              }
          }
        for(i=0;i<(jlen-jlen1)*nx;i++)
          {
                s_uuw[i+jlen1*nx]=c0;
                s_uue[i+jlen1*nx]=c0;
                s_vus[i+jlen1*nx]=c0;
                s_vun[i+jlen1*nx]=c0;
          }
      } // j field4 >
      if(tavgyes==1)
        {
//        athread_put(0,s_uue,uue+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
        athread_put(0,s_uue,uvewns+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
        athread_put(0,s_vun,uvewns+jst*nx+2*nx*ny,s_double*nx*jlen,(void *)&put_reply,0,0);
        }

        while(get_reply1!=2);
        for(i=0;i<jlen*nx;j++)
         {
           s_wukb[i]=s_wuk[i]+p5c2dzk*(s_vun[i]-s_vus[i]+s_uue[i]-s_uuw[i])*s_uarear[i];
         }

        get_reply1=0;
        athread_get(0,kxu+jst*nx,s_kxu,jlen*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,kyu+jst*nx,s_kxu,jlen*nx*s_double,(void *)&get_reply1,0,0,0);
        athread_get(0,kmu+jst*nx,s_kxu,jlen*nx*s_int,   (void *)&get_reply1,0,0,0);

     if(tavgyes==1)
      {
        athread_put(0,s_wukb,wukb+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
      }
     else
      {
        athread_put(0,s_wukb,wuk+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
      }

//   luk lvk
    if(jed < jbeg-1 || jst > jend-1) // do j=jbeg,jend in F90, [jbeg-1,jend-1] in c   
      {
        for(i=0;i<jlen*nx;i++)
          {
            s_luk[i]=c0; s_lvk[i]=c0;
          }
      }//j field1>
    else if(jst <jbeg-1 && jed >=jbeg-1)
      {
        jlen1=jbeg-2-jst+1; //[jst,jbeg-2]
        for(i=0;i<jlen1*nx;i++)
          {
            s_luk[i]=c0; s_lvk[i]=c0;
          }
        for(j=0;j<jlen-jlen1;j++)
         {
           for(i=ibeg-1;i<iend;i++)
             {
               s_cc[i+jlen1*nx+j*nx]=s_vus[i+jlen1*nx+j*nx+nx]-s_vus[i+jlen1*nx+j*nx]
                                    +s_uuw[i+jlen1*nx+j*nx+1] -s_uuw[i+jlen1*nx+j*nx];
              s_luk[i+jlen1*nx+j*nx]=p5*(s_cc[i+jlen1*nx+j*nx]*s_uuu[i+jlen1*nx+j*nx+nx]+
                                         s_vus[i+jlen1*nx+(j+1)*nx]*s_uuu[i+jlen1*nx+(j+2)*nx]
                                        -s_vus[i+jlen1*nx+j*nx]*s_uuu[i+jlen1*nx+j*nx]
                                        +s_uuw[i+1+jlen1*nx+j*nx]*s_uuu[i+1+jlen1*nx+j*nx+nx]
                                        -s_uuw[i+jlen1*nx+j*nx]*s_uuu[i-1+jlen1*nx+j*nx+nx])
                                        *s_uarear[i+jlen1*nx+j*nx];
             s_lvk[i+jlen1*nx+j*nx]=p5*(s_cc[i+jlen1*nx+j*nx]*s_vvv[i+jlen1*nx+j*nx+nx]+
                                        s_vus[i+jlen1*nx+j*nx+nx]*s_vvv[i+jlen1*nx+(j+2)*nx]
                                       -s_vus[i+jlen1*nx+j*nx]*s_vvv[i+jlen1*nx+j*nx]
                                       +s_uuw[i+1+jlen1*nx+j*nx]*s_vvv[i+1+jlen1*nx+j*nx+nx]
                                       -s_uuw[i+jlen1*nx+j*nx]*s_vvv[i-1+jlen1*nx+j*nx+nx])
                                       *s_uarear[i+jlen1*nx+j*nx];


             }
           for(i=0;i<ibeg-1;i++)
             {
               s_luk[i]=c0; s_lvk[i]=0;
             }
           for(i=iend;i<nx;i++)
             {
               s_luk[i]=c0; s_lvk[i]=0;
             }

         }// j> 

      } //j field2 >
    else if(jst >=jbeg-1 && jed < jend) //[jbeg-1,jend-1]
      {
        for(j=0;j<jlen;j++)
          {
            for (i=ibeg-1;i<iend;i++)
              {
                s_cc[j*nx+i]=s_vus[j*nx+nx+i]-s_vus[j*nx+i]+s_uuw[i+1+j*nx]-s_uuw[j*nx+i];
                s_luk[j*nx+i]=p5*(s_cc[i+j*nx]*s_uuu[i+j*nx+nx]
                                 +s_vus[j*nx+nx+i]*s_uuu[i+(j+2)*nx]
                                 -s_vus[j*nx+i]   *s_uuu[i+j*nx]
                                 +s_uuw[i+1+j*nx]*s_uuu[i+1+j*nx+nx]
                                 -s_uuw[i+j*nx]*s_uuu[i-1+j*nx+nx])*s_uarear[j*nx+i];
                s_lvk[j*nx+i]=p5*(s_cc[j*nx+i]*s_vvv[i+j*nx+nx]
                                 +s_vus[i+j*nx+nx]*s_vvv[i+(j+2)*nx]
                                 -s_vus[i+j*nx]*s_vvv[i+j*nx]
                                 +s_uuw[i+1+j*nx]*s_vvv[i+1+j*nx+nx]
                                 -s_uuw[i+j*nx]*s_vvv[i-1+j*nx+nx])*s_uarear[j*nx+i]; 
              } // i>
            for(i=0;i<ibeg-1;i++)
              {
                s_luk[j*nx+i]=c0; s_lvk[j*nx+i]=c0;
              }
            for(i=iend;i<nx;i++)
              {
                s_luk[j*nx+i]=c0; s_lvk[j*nx+i]=c0;
              }
          }// j>
      } //j field 3>
    else if(jst <= jend-1 && jed > jend-1)
      {
         jlen1=jend-jst;
         for(j=0;j<jlen1;j++)          
          {
            for (i=ibeg-1;i<iend;i++)
              {
                s_cc[j*nx+i]=s_vus[j*nx+nx+i]-s_vus[j*nx+i]+s_uuw[i+1+j*nx]-s_uuw[j*nx+i];
                s_luk[j*nx+i]=p5*(s_cc[i+j*nx]*s_uuu[i+j*nx+nx]
                                 +s_vus[j*nx+nx+i]*s_uuu[i+(j+2)*nx]
                                 -s_vus[j*nx+i]   *s_uuu[i+j*nx]
                                 +s_uuw[i+1+j*nx]*s_uuu[i+1+j*nx+nx]
                                 -s_uuw[i+j*nx]*s_uuu[i-1+j*nx+nx])*s_uarear[j*nx+i];
                s_lvk[j*nx+i]=p5*(s_cc[j*nx+i]*s_vvv[i+j*nx+nx]
                                 +s_vus[i+j*nx+nx]*s_vvv[i+(j+2)*nx]
                                 -s_vus[i+j*nx]*s_vvv[i+j*nx]
                                 +s_uuw[i+1+j*nx]*s_vvv[i+1+j*nx+nx]
                                 -s_uuw[i+j*nx]*s_vvv[i-1+j*nx+nx])*s_uarear[j*nx+i];
              } // i>
            for(i=0;i<ibeg-1;i++)
              {
                s_luk[j*nx+i]=c0; s_lvk[j*nx+i]=c0;
              }
            for(i=iend;i<nx;i++)
              {
                s_luk[j*nx+i]=c0; s_lvk[j*nx+i]=c0;
              }
          }// j>
        for(i=0;j<(jlen-jlen1)*nx;i++)
          {
            s_luk[jlen1*nx+i]=c0;
            s_lvk[jlen1*nx+i]=c0;

          }

      } //j field4 >


    if(k==1)
      {
        for(i=0;i<jlen*nx;i++)
          {
            s_luk[i]+=dzrk*s_wuk[i]*s_uuu[i+nx];
            s_lvk[i]+=dzrk*s_wuk[i]*s_vvv[i+nx];
          }
      }
    else
      {
        while(get_reply2!=2);
        for(i=0;i<jlen*nx;i++)
          {
            s_luk[i]+=dz2rk*s_wuk[i]*(s_uuub[i]+s_uuu[i+nx]);
            s_lvk[i]+=dz2rk*s_wuk[i]*(s_vvvb[i]+s_vvv[i+nx]);

          }
      }

    if(k<km)
      {
        while(get_reply3!=2);
        for(i=0;i<jlen*nx;i++)
          {
            s_luk[i]-=dz2rk*s_wukb[i]*(s_uuu[i+nx]+s_uuua[i]);
            s_lvk[i]-=dz2rk*s_wukb[i]*(s_vvv[i+nx]+s_vvva[i]);

          }


      }


   } // partial else>

    while(get_reply1!=3);
    j1=max(0,ibeg-1-jst);
    j2=min(jlen,jend-jst); 
    for(j=j1;j<j2;j++)
      {
        for(i=ibeg-1;i<iend;i++)
          {
            if(k<=s_kmu[j*nx+i])
              {
                s_luk[j*nx+i]+=s_uuu[j*nx+nx+i]*s_vvv[j*nx+nx+i]*s_kyu[j*nx+i]
                              -s_vvv[j*nx+nx+i]*s_vvv[j*nx+nx+i]*s_kxu[j*nx+i]; 
                s_lvk[j*nx+i]+=s_uuu[j*nx+nx+i]*s_vvv[j*nx+nx+i]*s_kxu[j*nx+i]
                              -s_uuu[j*nx+nx+i]*s_uuu[j*nx+nx+i]*s_kyu[j*nx+i]; 
              }
            else
              {
                s_luk[j*nx+i]=c0;
                s_lvk[j*nx+i]=c0;
              }
          }

      }// j>


//     athread_put(0,s_uuw,uuw+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
//     athread_put(0,s_uue,uue+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
//     athread_put(0,s_vus,vus+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
//     athread_put(0,s_vun,vun+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
     athread_put(0,s_luk,luk+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
     athread_put(0,s_lvk,lvk+jst*nx,s_double*nx*jlen,(void *)&put_reply,0,0);
     if(tavgyes==0) while(put_reply!=3);
     else while(put_reply!=5);
 }




} // a >


//} 
