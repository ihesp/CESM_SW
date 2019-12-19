#include <stdio.h>
#include <athread.h>
#include <math.h>
#include "advection_struct.h"

#define locmin(x,y) ((x)<(y)?(x) : (y))     

extern SLAVE_FUN(adhdt_loop11_fun)();
extern SLAVE_FUN(adhdt_loop13_fun)();


void adhdt_loop11_fun_(int* intpara, 
      double* UUU, double* VVV, double* DYU, double* DXU,double* DZU,
      double* UVNSEW, double* DTN, double* DTS, double* DTE, double* DTW, double* DZT,
      double* CNSEW, int* KMT,double *WTK,double *WTKB,double *TAREA_R)

{
    int nsc, nscad, j;
    int nxy;
    
    struct adhdt_Loop11_Para Loop11_Para;

    Loop11_Para.param[0]= intpara[0];  //*jb-1;
    Loop11_Para.param[1]= intpara[1];  //*je+1;
    Loop11_Para.param[2]= intpara[2];  //*ib-1;
    Loop11_Para.param[3]= intpara[3];  //*ie+1;
    Loop11_Para.param[4]= intpara[4];  //*nx_block;
    Loop11_Para.param[5]= intpara[5];  //*ny_block;
    Loop11_Para.param[6]= intpara[6];  //*array_z;	
    Loop11_Para.param[7]= intpara[7];  //*bid;
    Loop11_Para.param[8]= intpara[8];  //*k;
    Loop11_Para.param[9]= intpara[9];  //*km;
    Loop11_Para.param[10]= intpara[10]; //*my_task;
    Loop11_Para.param[13]= intpara[14];  //*jb;
    Loop11_Para.param[14]= intpara[15];  //*je;
    Loop11_Para.param[15]= intpara[16];  //*ib;
    Loop11_Para.param[16]= intpara[17];  //*ie;	
    nxy = intpara[4]*intpara[5];  	

//advt:	
    nsc = locmin(24,intpara[5]);
    Loop11_Para.param[17]=nsc;
    Loop11_Para.param[18] = 0;	     // ser. No. of pthread		
   
    Loop11_Para.use_arr[0]= UUU;
    Loop11_Para.use_arr[1]= VVV;
    Loop11_Para.use_arr[2]= DYU;
    Loop11_Para.use_arr[3]= DXU;
    Loop11_Para.use_arr[4]= DZU;
    Loop11_Para.use_arr[5]= UVNSEW;
    Loop11_Para.use_arr[6]= DTN;
    Loop11_Para.use_arr[7]= DTS;
    Loop11_Para.use_arr[8]= DTE;
    Loop11_Para.use_arr[9]= DTW;
    Loop11_Para.use_arr[10]= DZT;
    Loop11_Para.use_arr[11]= CNSEW;
    Loop11_Para.use_arr[12]= WTK;
    Loop11_Para.use_arr[13]= WTKB;
    Loop11_Para.use_arr[14]= TAREA_R;

	
    Loop11_Para.duse_arr[0]= KMT;	

//hdifft:
    nsc = locmin(24,intpara[5]);
    Loop11_Para.param[19]=nsc;
    Loop11_Para.param[20] = 24;	     // ser. No. of pthread	

    athread_spawn(adhdt_loop11_fun,&Loop11_Para);
    athread_join();
 }





   
void adhdt_loop13_fun_(int* intpara, 
   double* UVNSEW, double* TRCR, double* TAREA_R, double* DZT, double* LTK,double* WTK,double* WTKB,
   double* ah, double* dzrk, double* AHF, double* D2TK, 
   double* HDTK, double* TMIX, double* CNSEW)

{
    int nsc, n;
   	
    struct adhdt_Loop13_Para Loop13_Para;
  
    Loop13_Para.param[0]= intpara[0];  //*jb-1;
    Loop13_Para.param[1]= intpara[1];  //*je+1;
    Loop13_Para.param[2]= intpara[2];  //*ib-1;
    Loop13_Para.param[3]= intpara[3];  //*ie+1;
    Loop13_Para.param[4]= intpara[4];  //*nx_block;
    Loop13_Para.param[5]= intpara[5];  //*ny_block;
    Loop13_Para.param[6]= intpara[6];  //*array_z;	
    Loop13_Para.param[7]= intpara[7];  //*bid;
    Loop13_Para.param[8]= intpara[8];  //*k;
    Loop13_Para.param[9]= intpara[9];  //*km;
    Loop13_Para.param[10]= intpara[10]; //*my_task;
    Loop13_Para.param[11]= intpara[11]; //*nt;
    Loop13_Para.param[12]= intpara[12]; //*ns;
    Loop13_Para.param[13]= intpara[14];  //*jb;
    Loop13_Para.param[14]= intpara[15];  //*je;
    Loop13_Para.param[15]= intpara[16];  //*ib;
    Loop13_Para.param[16]= intpara[17];  //*ie;	
    Loop13_Para.param[21]= intpara[19];  //*sfc_layer_type;		

//advt:	
    nsc = locmin(24,intpara[5]);
    Loop13_Para.param[17]=nsc;
    Loop13_Para.param[18] = 0;	     // ser. No. of pthread		
   
//    Loop13_Para.param[12]= intpara[12]; //*tr_mask(n);	
    for(n=intpara[12]-1;n<intpara[11];n++)
       Loop13_Para.param[25+n]= intpara[20+n]; //*tr_mask(n);
    
    Loop13_Para.use_arr[0]= UVNSEW;
    Loop13_Para.use_arr[1]= TRCR;
    Loop13_Para.use_arr[2]= TAREA_R;
    Loop13_Para.use_arr[3]= DZT;
    Loop13_Para.use_arr[4]= LTK;
    Loop13_Para.use_arr[5]= AHF;
    Loop13_Para.use_arr[6]= HDTK;
    Loop13_Para.use_arr[7]= TMIX;
    Loop13_Para.use_arr[8]= CNSEW;
    Loop13_Para.use_arr[9]= D2TK;
    Loop13_Para.use_arr[10]= WTK;
    Loop13_Para.use_arr[11]= WTKB;	

//hdifft:
    nsc = locmin(24,intpara[5]);
    Loop13_Para.param[19] = nsc;
    Loop13_Para.param[20] = 24;	     // ser. No. of pthread	

    Loop13_Para.dparam[0] = *ah;
    Loop13_Para.dparam[1] = *dzrk;	

 
    athread_spawn(adhdt_loop13_fun,&Loop13_Para);
    athread_join();
 }


