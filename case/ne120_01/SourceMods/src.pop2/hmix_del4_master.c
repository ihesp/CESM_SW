#include <stdio.h>
#include <athread.h>
#include <math.h>
#include "hmix_del4_struct.h"

extern SLAVE_FUN(s_advuhdiffu)();


void advuhdiffu_c_(int *iparam, double *rparam, 
                 double *uuu, double *vvv, double *dxu, double *dyu, double *dzu,
                 double *uvewns, double *wukb, double *wuk, double *uarea_r,
                 double *luk, double *lvk, double *kxyu,
                 double *duc,double *dum,
                 double *dun, double *dus, double *due, double *duw, 
                 double *dmc, double *dmn, double *dms, double *dme, double *dmw,
                 double *umixk, double *vmixk, 
                 double *hduk, double *hdvk, double *amf, int *kmu)


{
     int i;
     struct param_hdiffu hdiffu;

     for(i=0;i<12;i++) hdiffu.param[i]=*(iparam+i);
     hdiffu.param[12]=athread_get_max_threads();     
     for(i=0;i<9;i++) hdiffu.dparam[i]=*(rparam+i); 

     hdiffu.addr[0]=uuu;
     hdiffu.addr[1]=vvv;
     hdiffu.addr[2]=dxu;
     hdiffu.addr[3]=dyu;
     hdiffu.addr[4]=dzu;
     hdiffu.addr[5]=uvewns;
     hdiffu.addr[6]=wukb;
     hdiffu.addr[7]=wuk;
     hdiffu.addr[8]=uarea_r;
     hdiffu.addr[9]=luk;
     hdiffu.addr[10]=lvk;
     hdiffu.addr[11]=kxyu;

     hdiffu.addr[12]=duc;
     hdiffu.addr[13]=dum;
     hdiffu.addr[14]=dun;
     hdiffu.addr[15]=dus;
     hdiffu.addr[16]=due;
     hdiffu.addr[17]=duw;
     hdiffu.addr[18]=dmc;
     hdiffu.addr[19]=dmn;
     hdiffu.addr[20]=dms;
     hdiffu.addr[21]=dme;
     hdiffu.addr[22]=dmw;
     hdiffu.addr[23]=umixk;
     hdiffu.addr[24]=vmixk;
     hdiffu.addr[25]=hduk;
     hdiffu.addr[26]=hdvk;
     hdiffu.addr[27]=amf;

     hdiffu.addr_int[0]=kmu;
     
     athread_spawn(s_advuhdiffu,&hdiffu);         
     athread_join();      
}
