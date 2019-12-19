#include<athread.h>
#include<math.h>
//#include<stdio.h>
#include "vertical_mix_struct.h"
#include <mpi.h>
 
/************ master  *************/
extern SLAVE_FUN(s_vmix)();
extern SLAVE_FUN(s_impt)();
extern SLAVE_FUN(s_vdift)();
extern SLAVE_FUN(s_impu)();
extern SLAVE_FUN(s_impt_cr)();

int rank;
struct param_vmix v1;     //vdiffu 
struct param_vmix vit;    //impvmixt
struct param_vmix vdift;  //vdifft
struct param_vmix vdiu;   //impvmixu
struct param_vmix vitc;   //impvmix_correct


//############ vertical_mix ####################### 
void vdiffu_c_(int *nx, int *ny, int *k, int *km, int *ist, int *ied, int *jst, int *jed, double *p5, 
              double *drag, double *dzwrk, double *dzrk, double *uold, double *vold, double *vvc, 
              double *dzu, double *smf, double *vuf,double *vvf,double *vduk, double *vdvk, 
              int *kmu, double *c0, int *partial)
 {
   v1.param[0]=*nx;
   v1.param[1]=*ny;
   v1.param[2]=*k;
   v1.param[3]=*km;
   v1.param[4]=*ist;
   v1.param[5]=*ied;
   v1.param[6]=*jst;
   v1.param[7]=*jed;
   v1.param[8]=*partial;
   v1.param[9]=athread_get_max_threads();

   v1.dparam[0]=*p5; 
   v1.dparam[1]=*drag; 
   v1.dparam[2]=*dzwrk;
   v1.dparam[3]=*dzrk; 
   v1.dparam[4]=*c0;

   v1.addr[0]=uold;
   v1.addr[1]=vold;
   v1.addr[2]=vvc;
   v1.addr[3]=dzu;
   v1.addr[4]=smf;
   v1.addr[5]=vuf;
   v1.addr[6]=vvf;
   v1.addr[7]=vduk;
   v1.addr[8]=vdvk;

   
   v1.addr_int[0]=kmu;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);


//   if (*k<2&& rank==0) {
//   printf("nx,ny,k,km %d %d %d\n",*nx,*ny,*k,*km);
//   printf("k, uold vold vvc dzu vufb vvfb dzu%d %d %d %d %d %d %d\n",*k,uold, 
//        vold,vvc,dzu,vufb,vvfb); }

   athread_spawn(s_vmix,&v1);
   athread_join();
  }


/*void c_addr2_(int *i, int *j,int *k, double *z, double *t)
 {  
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if(rank ==0){
     printf("k,mt2,bid,vdc::kbid,vdc %d %d %d %d\n", *i,*j,*k,z,t);
     }
 }
*/


void impvmixt_c_(int *nfirst,int *nlast,int *nx,int *ny,int *nt,int *km,int *mt2_r, int *sfc_ju,
                 double *grav,double *p5, double *c0, double *aidif,
                 double *dz, double *c2dtt, double *afact,  
                 double *psfc, double *vdc, double *dzt, double *tnew, double *told,int *kmt, 
                 int *partial,int *ib, int *ie, int *jb, int *je)
{
      vit.param[0]  =  *nfirst;
      vit.param[1]  =  *nlast;
      vit.param[2]  =  *nx;
      vit.param[3]  =  *ny;
      vit.param[4]  =  *nt;
      vit.param[5]  =  *km;
      vit.param[6]  =  *mt2_r;
      vit.param[7]  =  *sfc_ju;
      vit.param[8]  =  *partial;
      vit.param[9]  =  athread_get_max_threads();
      vit.param[10]  = *ib;
      vit.param[11]  = *ie;
      vit.param[12]  = *jb;
      vit.param[13]  = *je;

      vit.dparam[0] =  *grav; 
      vit.dparam[1] =  *p5; 
      vit.dparam[2] =  *c0; 
      vit.dparam[3] =  *aidif;

      vit.addr[0]   =  dz;
      vit.addr[1]   =  c2dtt;
      vit.addr[2]   =  afact;
      vit.addr[3]   =  psfc;
      vit.addr[4]   =  vdc;
      vit.addr[5]   =  dzt;
      vit.addr[6]   =  tnew;
      vit.addr[7]   =  told;
//      vit.addr[8]   =  ff;

      vit.addr_int[0]  = kmt;

//              MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      
      
//         if ( rank==0) {
//         printf("nx,ny,k,km,nfirst,nlast %d %d %d %d %d %d\n",*nfirst, *nlast, *nx,*ny,*nt,*km);
//         printf("cgrav,p5,c0,aidif  %f %f %f %f\n",*grav, *p5,*c0,*aidif);
//         printf("cdz%f %f %f %f %f\n",*dz, *(dz+1),*(dz+2),*(dz+3),*(dz+4));
//         printf("cc2dtt%f %f %f %f %f\n",*c2dtt, *(c2dtt+1),*(c2dtt+2),*(c2dtt+3),*(c2dtt+4));
//         printf("cafact%f %f %f %f %f\n",*afact, *(afact+1),*(afact+2),*(afact+3),*(afact+4));
//         printf("told,tnew,ff,afact %d %d %d %d \n",told,tnew,ff,afact);
//         printf("k, uold vold vvc dzu vufb vvfb dzu %d %d\n",*sfc_judge,*partial);
//               }
      

      athread_spawn(s_impt,&vit);
      athread_join();
     
}

void vdifft_c_(int *iparam, double *rparam, double *stf, double *vdc, double *dzt, double *told, double *vtf, double *vdtk, int *kmt)
 {
   int i;
   for (i=0;i<9;i++){vdift.param[i]=*(iparam+i);}
   for (i=0;i<8;i++){vdift.dparam[i]=*(rparam+i);}

   vdift.param[9]=athread_get_max_threads();
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   vdift.param[10]=rank;

   vdift.addr[0]=stf;
   vdift.addr[1]=vdc;
   vdift.addr[2]=dzt;
   vdift.addr[3]=told;
   vdift.addr[4]=vtf;
   vdift.addr[5]=vdtk;

   vdift.addr_int[0]=kmt;
#ifdef POPcheckJP
   if(rank==0) {
printf("xxxxxx%d %d %d %d %d %d %d %d\n",*iparam,*(iparam+1),*(iparam+2),*(iparam+3),*(iparam+4),*(iparam+5),*(iparam+6),*(iparam+7));
printf("xxxxxx%f %f %f %f %f %f %f %f\n",*rparam,*(rparam+1),*(rparam+2),*(rparam+3),*(rparam+4),*(rparam+5),*(rparam+6),*(rparam+7));
              }
#endif 

   athread_spawn(s_vdift,&vdift);
   athread_join();
 }



void impvu_c_(int *iparam,double *rparam, double *dz,double *vvc,
              double *unew,double *vnew,double *uold, double *vold, double *dzu,
              double *afacu, int *kmu)
{
  int i,rank;
  for (i=0;i<8;i++){vdiu.param[i]=*(iparam+i);}
  for (i=0;i<4;i++){vdiu.dparam[i]=*(rparam+i);}
  vdiu.param[8]=athread_get_max_threads();
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  vdiu.param[9]=rank;
  vdiu.addr[0]=dz;
  vdiu.addr[1]=vvc;
  vdiu.addr[2]=unew;
  vdiu.addr[3]=vnew;
  vdiu.addr[4]=uold;
  vdiu.addr[5]=vold;
  vdiu.addr[6]=dzu;
  vdiu.addr[7]=afacu;
  
  vdiu.addr_int[0]=kmu;

  athread_spawn(s_impu,&vdiu);
  athread_join();
}

void impt_cr_c_(int *iparam, double *rparam,double *dz, double *c2dtt,
                double *afact,double *psfc, double *vdc,double *dzt,
                double *rhs, double *tnew, int *kmt)
{
  int i;
  for(i=0;i<14;i++) {vitc.param[i]=*(iparam+i);}
  for(i=0;i<4;i++)  {vitc.dparam[i]=*(rparam+i);}
    vitc.param[14]=athread_get_max_threads();
  
    vitc.addr[0]=dz;
    vitc.addr[1]=c2dtt;
    vitc.addr[2]=afact;
    vitc.addr[3]=psfc;
    vitc.addr[4]=vdc;
    vitc.addr[5]=dzt;
    vitc.addr[6]=rhs;
    vitc.addr[7]=tnew;

    vitc.addr_int[0]=kmt;
    
    athread_spawn(s_impt_cr,&vitc);
    athread_join();
}
