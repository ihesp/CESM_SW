#include<athread.h>
#include<math.h>
//#include<stdio.h>
#include "baroclinic_struct.h"
#include<mpi.h>
 
/************ master  *************/
extern SLAVE_FUN(s_clsec1)();
extern SLAVE_FUN(s_postadvu)();
extern SLAVE_FUN(s_storeforce)();
//extern SLAVE_FUN(s_grad)();



//############ clinic_sec2_  ###################### 
void clinic_vertav_normalv_c_(int *nx, int *ny, int *km,int *partial, double *c0, double *unew,
                              double *vnew, double *dzu, double *hur, double *dz, int *kmu,int *l1ddyn)
 {
   struct param_clinic clsec1;
   int rank;
   clsec1.param[0]=*nx;
   clsec1.param[1]=*ny;
   clsec1.param[2]=*km;
   clsec1.param[3]=*partial;
   clsec1.param[4]=*l1ddyn;
   clsec1.param[5]=athread_get_max_threads();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   clsec1.param[6]=rank;

   clsec1.dparam[0]=*c0; 

   clsec1.addr[0]=unew;
   clsec1.addr[1]=vnew;
   clsec1.addr[2]=dzu;
   clsec1.addr[3]=hur;
   clsec1.addr[4]=dz;

   clsec1.addr_int[0]=kmu;

   athread_spawn(s_clsec1,&clsec1);
   athread_join();
 }
//############ post_advu ####################### 
void post_advu_c_(int *iparam, double *rparam, 
                  double *diag_ke_adv_2d,double *dzu,double *ucur,double *vcur,
                  double *uold, double *vold, double *fcor,double *fx,double *fy,
                  double *rhoknew,double *rhokcur,double *rhokold,
                  double *sumx, double *sumy,double *rhokmx, double *rhokmy,
                  double *dxur, double *dyur,double *workx,int *kmu)
//                  double *dxur, double *dyur,double *dkpress2d,int *kmu)
 {
   struct param_clinic postadvu;
   int rank,i;
   for(i=0;i<8;i++) postadvu.param[i]=*(iparam+i);
   for(i=0;i<10;i++) postadvu.dparam[i]=*(rparam+i); 
   postadvu.param[8] = athread_get_max_threads();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   postadvu.param[9] = rank;
   postadvu.addr[0]=diag_ke_adv_2d;
   postadvu.addr[1]=dzu;
   postadvu.addr[2]=ucur;
   postadvu.addr[3]=vcur;
   postadvu.addr[4]=uold;
   postadvu.addr[5]=vold;
   postadvu.addr[6]=fcor;
   postadvu.addr[7]=fx;
   postadvu.addr[8]=fy;
   postadvu.addr[9]=rhoknew;
   postadvu.addr[10]=rhokcur;
   postadvu.addr[11]=rhokold;
   postadvu.addr[12]=sumx;
   postadvu.addr[13]=sumy;
   postadvu.addr[14]=rhokmx;
   postadvu.addr[15]=rhokmy;
   postadvu.addr[16]=dxur;
   postadvu.addr[17]=dyur;
   postadvu.addr[18]=workx;
//   postadvu.addr[18]=dkpress2d;
   postadvu.addr_int[0]=kmu;
   athread_spawn(s_postadvu,&postadvu);
   athread_join();
 }

    void storeforce_c_(int *iparam, double *rparam, double *fcor, 
                       double *fx, double *fy, double *uvel, double *vvel,
                       double *zx, double *zy,
                       double *dzu)
{
       int i;
       struct little_clinic stforce;
       for(i=0;i<6;i++) stforce.param[i]=*(iparam+i);
       for(i=0;i<4;i++) stforce.dparam[i]=*(rparam+i);
       stforce.param[6]= athread_get_max_threads();
       stforce.addr[0]=fcor;
       stforce.addr[1]=fx;
       stforce.addr[2]=fy;
       stforce.addr[3]=uvel;
       stforce.addr[4]=vvel;
       stforce.addr[5]=zx;
       stforce.addr[6]=zy;
       stforce.addr[7]=dzu;
       athread_spawn(s_storeforce,&stforce);
       athread_join();

}

/*
     void grad_c_(int *nx,int *ny,int *k,double *c0,double *p5,double *rhokx,double *rhoky,
                 double *rhoavg,double *dxur,double *dyur,int *kmu)
{
     struct param_clinic grad;
     grad.param[0]=*nx;
     grad.param[1]=*ny;
     grad.param[2]=*k;
     grad.param[3]=athread_get_max_threads();
     grad.dparam[0]=*c0;
     grad.dparam[1]=*p5;
     grad.addr[0]=rhokx;
     grad.addr[1]=rhoky;
     grad.addr[2]=rhoavg;
     grad.addr[3]=dxur;
     grad.addr[3]=dyur;
     grad.addr_int[0]=kmu;
     athread_spawn(s_grad,&grad);
     athread_join();

}     
*/
