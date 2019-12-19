#include <stdio.h>
#include <math.h>
#include <slave.h>
#include <unistd.h>
#include "clinic_tavg_accum.h"

#ifdef POPTESTJN
#include "cpe_print.h"
#endif

#define SZldmint 2200
#define SZldmbuf 1100
#define SZldmdbl 4400
#define locmax(a,b)  ((a)>(b)?(a):(b))
#define locmin(a,b)  ((a)<(b)?(a):(b))


#define sync_mem   \
asm volatile("memb")
//asm volatile("memb\n\t":::"memory")



void s_clinic_tavg_0();
//accumulate_tavg_field(UVEL(:,:,k,curtime,iblock),tavg_UVEL,iblock,k)

void s_clinic_tavg_2();
//accumulate_tavg_field(UVEL(:,:,k,curtime,iblock)**2,tavg_UVEL2,iblock,k)

void s_clinic_tavg_5();
//ccumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_VVEL,iblock,k)

void s_clinic_tavg_7();
//accumulate_tavg_field(VVEL(:,:,k,curtime,iblock)**2,tavg_VVEL2,iblock,k)

void s_clinic_tavg_10();
//accumulate_tavg_field(p5*(UVEL(:,:,k,curtime,iblock)**2 + &
//                          VVEL(:,:,k,curtime,iblock)**2),tavg_KE,iblock,k)

void s_clinic_tavg_12();
//accumulate_tavg_field(UVEL(:,:,k,curtime,iblock)* &
//                      VVEL(:,:,k,curtime,iblock), tavg_UV,iblock,k)

void s_clinic_tavg_13();
//accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), tavg_TEMP,iblock,k)

void s_clinic_tavg_17();
//accumulate_tavg_field(max(TRACER(:,:,k,1,curtime,iblock) - &
//            TRACER(:,:,k,1,oldtime,iblock), c0), tavg_dTEMP_POS_3D,iblock,k)

void s_clinic_tavg_18();
//accumulate_tavg_field(max(TRACER(:,:,k,1,curtime,iblock) - &
//            TRACER(:,:,k,1,oldtime,iblock), c0), tavg_dTEMP_POS_2D,iblock,1)

void s_clinic_tavg_19();
//accumulate_tavg_field(min(TRACER(:,:,k,1,curtime,iblock) - &
//            TRACER(:,:,k,1,oldtime,iblock), c0), tavg_dTEMP_NEG_3D,iblock,k)

void s_clinic_tavg_20();
//accumulate_tavg_field(min(TRACER(:,:,k,1,curtime,iblock) - &
//            TRACER(:,:,k,1,oldtime,iblock), c0), tavg_dTEMP_NEG_2D,iblock,1)

void s_clinic_tavg_22();
//accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)**2,tavg_TEMP2,iblock,k)

void s_clinic_tavg_25();
//accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock),tavg_SALT,iblock,k)

void s_clinic_tavg_30();
//accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock)**2,tavg_SALT2,iblock,k)

void s_clinic_tavg_31();
//accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)* &
//            TRACER(:,:,k,2,curtime,iblock), tavg_ST,iblock,k)

void s_clinic_tavg_32();
//accumulate_tavg_field(RHO(:,:,k,curtime,iblock), tavg_RHO,iblock,k)

void s_clinic_tavg_33();
/*
if (sfc_layer_type == sfc_layer_varthick .and. k == 1) then
    where (k <= KMT(:,:,iblock))
       WORK1 = RHO(:,:,k,curtime,iblock) * ( dz(k) + PSURF(:,:,curtime,iblock)/grav )
    elsewhere
       WORK1 = c0
    endwhere
else
  if (partial_bottom_cells) then
    where (k <= KMT(:,:,iblock))
       WORK1 = RHO(:,:,k,curtime,iblock) * DZT(:,:,k,iblock)
    elsewhere
       WORK1 = c0
    endwhere
  else
    where (k <= KMT(:,:,iblock))
       WORK1 = RHO(:,:,k,curtime,iblock) * dz(k)
    elsewhere
       WORK1 = c0
    endwhere
  endif
endif
accumulate_tavg_field(WORK1,tavg_RHO_VINT,iblock,k)
*/

void s_clinic_tavg_34();
/*
if ( sfc_layer_type /= sfc_layer_varthick .and. k == 1) then
  if (accumulate_tavg_now(tavg_RESID_T)) then
     WORK1 = c0
     factor = c1/hflux_factor  ! converts to W/m^2
     where (CALCT(:,:,iblock))  &
        WORK1=DH(:,:,iblock)*TRACER(:,:,1,1,curtime,iblock)*factor
        accumulate_tavg_field(WORK1,tavg_RESID_T,iblock,k)
  endif
endif
*/

void s_clinic_tavg_35();
/*
if ( sfc_layer_type /= sfc_layer_varthick .and. k == 1) then
  if (accumulate_tavg_now(tavg_RESID_S)) then
     WORK1 = c0
     factor = c1/salinity_factor  ! converts to kg(freshwater)/m^2/s
     where (CALCT(:,:,iblock)) &
        WORK1 = DH(:,:,iblock)*TRACER(:,:,k,2,curtime,iblock)*factor
        accumulate_tavg_field(WORK1,tavg_RESID_S,iblock,k)
     endif
endif
*/




//---------------:---------------------------------------------------
void s_clinic_tavg(struct param_tavg_accum *st2_master)
 {
	
// BOP:
 volatile int get_reply, put_reply;

// Common variables:
 int myproc,myid,threads_num;

// Input variables: 
 int numtavg;

// Output variables:

// Slave & Temp variables:
 struct param_tavg_accum st2_slave, *st_slave;

// Constants && parameters:

// EOP.


 get_reply=0;
 athread_get(PE_MODE,st2_master,&st2_slave,sizeof(st2_slave),(void*)&get_reply,0,0,0);

 myid = athread_get_id(-1);

 st_slave = &st2_slave;

 while(get_reply!=1);
 sync_mem;


  numtavg= st2_slave.ipts[10];
  threads_num= st2_slave.ipts[11];

  switch(myid) 
     {
     case 0:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 8;
        st2_slave.tavgindx[6] = 0;
        s_clinic_tavg_0(st_slave);
//call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock),tavg_UVEL,iblock,k)
        break;

     case 1:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 8;
        st2_slave.tavgindx[6] = 100;
        s_clinic_tavg_0(st_slave);
        break;
     
     case 2:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 8;
        st2_slave.tavgindx[6] = 200;
        s_clinic_tavg_0(st_slave);
        break;
     
     case 3:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 8;
        st2_slave.tavgindx[6] = 300;
        s_clinic_tavg_0(st_slave);
        break;
     
     case 4:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 8;
        st2_slave.tavgindx[6] = 400;
        s_clinic_tavg_0(st_slave);
        break;
     
     case 5:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 8;
        st2_slave.tavgindx[6] = 500;
        s_clinic_tavg_0(st_slave);
        break;
     
     case 6:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 8;
        st2_slave.tavgindx[6] = 600;
        s_clinic_tavg_0(st_slave);
        break;
     
     case 7:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 8;
        st2_slave.tavgindx[6] = 700;
        s_clinic_tavg_0(st_slave);
        break;
     


     case 8:
        st2_slave.ipts[5] = myid*100 + 1;
        st2_slave.tavgindx[12+5] = 1;
        st2_slave.tavgindx[12+6] = 1;
        s_clinic_tavg_0(st_slave);
//call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock),tavg_UVEL_2,iblock,k)

        st2_slave.ipts[5] = myid*100 + 2;
        st2_slave.tavgindx[2*12+5] = 1;
        st2_slave.tavgindx[2*12+6] = 2;
        s_clinic_tavg_2(st_slave);
//call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock)**2,tavg_UVEL2,iblock,k)
        break;

     case 9:
        st2_slave.ipts[5] = myid*100 + 3;
        st2_slave.tavgindx[3*12+5] = 1;
        st2_slave.tavgindx[3*12+6] = 3;
        s_clinic_tavg_0(st_slave);
//if (k <= 1)  &
//     call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock), &
//              tavg_U1_1,iblock,k)

        st2_slave.ipts[5] = myid*100 + 4;
        st2_slave.tavgindx[4*12+5] = 1;
        st2_slave.tavgindx[4*12+6] = 4;
        s_clinic_tavg_0(st_slave);
//if (k <= 8)  &
//     call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock),tavg_U1_8,iblock,k)
        break;



     case 10:
        st2_slave.ipts[5] = myid*100 + 5;
        st2_slave.tavgindx[5*12+5] = 8;
        st2_slave.tavgindx[5*12+6] = 5;
        s_clinic_tavg_5(st_slave);
//call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_VVEL,iblock,k)
        break;

     case 11:
        st2_slave.ipts[5] = myid*100 + 5;
        st2_slave.tavgindx[5*12+5] = 8;
        st2_slave.tavgindx[5*12+6] = 105;
        s_clinic_tavg_5(st_slave);
        break;

     case 12:
        st2_slave.ipts[5] = myid*100 + 5;
        st2_slave.tavgindx[5*12+5] = 8;
        st2_slave.tavgindx[5*12+6] = 205;
        s_clinic_tavg_5(st_slave);
        break;

     case 13:
        st2_slave.ipts[5] = myid*100 + 5;
        st2_slave.tavgindx[5*12+5] = 8;
        st2_slave.tavgindx[5*12+6] = 305;
        s_clinic_tavg_5(st_slave);
        break;

     case 14:
        st2_slave.ipts[5] = myid*100 + 5;
        st2_slave.tavgindx[5*12+5] = 8;
        st2_slave.tavgindx[5*12+6] = 405;
        s_clinic_tavg_5(st_slave);
        break;

     case 15:
        st2_slave.ipts[5] = myid*100 + 5;
        st2_slave.tavgindx[5*12+5] = 8;
        st2_slave.tavgindx[5*12+6] = 505;
        s_clinic_tavg_5(st_slave);
        break;

     case 16:
        st2_slave.ipts[5] = myid*100 + 5;
        st2_slave.tavgindx[5*12+5] = 8;
        st2_slave.tavgindx[5*12+6] = 605;
        s_clinic_tavg_5(st_slave);
        break;

     case 17:
        st2_slave.ipts[5] = myid*100 + 5;
        st2_slave.tavgindx[5*12+5] = 8;
        st2_slave.tavgindx[5*12+6] = 705;
        s_clinic_tavg_5(st_slave);
        break;

     

     case 18:
        st2_slave.ipts[5] = myid*100 + 6;
        st2_slave.tavgindx[6*12+5] = 1;
        st2_slave.tavgindx[6*12+6] = 6;
        s_clinic_tavg_5(st_slave);
//call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_VVEL_2,iblock,k)

        st2_slave.ipts[5] = myid*100 + 7;
        st2_slave.tavgindx[7*12+5] = 1;
        st2_slave.tavgindx[7*12+6] = 7;
        s_clinic_tavg_7(st_slave);
//call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock)**2,tavg_VVEL2,iblock,k)
        break;

     case 19:
        st2_slave.ipts[5] = myid*100 + 8;
        st2_slave.tavgindx[8*12+5] = 1;
        st2_slave.tavgindx[8*12+6] = 8;
        s_clinic_tavg_5(st_slave);
//if (k <= 1)  &
//     call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_V1_1,iblock,k)
      
        st2_slave.ipts[5] = myid*100 + 9;
        st2_slave.tavgindx[9*12+5] = 1;
        st2_slave.tavgindx[9*12+6] = 9;
        s_clinic_tavg_5(st_slave);
//if (k <= 8)  &
//     call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock),tavg_V1_8,iblock,k)
        break;




     case 20:
        st2_slave.ipts[5] = myid*100 + 10;
        st2_slave.tavgindx[10*12+5] = 8;
        st2_slave.tavgindx[10*12+6] = 10;
        s_clinic_tavg_10(st_slave);
//         call accumulate_tavg_field(p5*(uvel(:,:,k,curtime,iblock)**2 + &
//                                        vvel(:,:,k,curtime,iblock)**2),tavg_KE,iblock,k)
        break;

     case 21:
        st2_slave.ipts[5] = myid*100 + 10;
        st2_slave.tavgindx[10*12+5] = 8;
        st2_slave.tavgindx[10*12+6] = 110;
        s_clinic_tavg_10(st_slave);
        break;

     case 22:
        st2_slave.ipts[5] = myid*100 + 10;
        st2_slave.tavgindx[10*12+5] = 8;
        st2_slave.tavgindx[10*12+6] = 210;
        s_clinic_tavg_10(st_slave);
        break;

     case 23:
        st2_slave.ipts[5] = myid*100 + 10;
        st2_slave.tavgindx[10*12+5] = 8;
        st2_slave.tavgindx[10*12+6] = 310;
        s_clinic_tavg_10(st_slave);
        break;

     case 24:
        st2_slave.ipts[5] = myid*100 + 10;
        st2_slave.tavgindx[10*12+5] = 8;
        st2_slave.tavgindx[10*12+6] = 410;
        s_clinic_tavg_10(st_slave);
        break;

     case 25:
        st2_slave.ipts[5] = myid*100 + 10;
        st2_slave.tavgindx[10*12+5] = 8;
        st2_slave.tavgindx[10*12+6] = 510;
        s_clinic_tavg_10(st_slave);
        break;

     case 26:
        st2_slave.ipts[5] = myid*100 + 10;
        st2_slave.tavgindx[10*12+5] = 8;
        st2_slave.tavgindx[10*12+6] = 610;
        s_clinic_tavg_10(st_slave);
        break;

     case 27:
        st2_slave.ipts[5] = myid*100 + 10;
        st2_slave.tavgindx[10*12+5] = 8;
        st2_slave.tavgindx[10*12+6] = 710;
        s_clinic_tavg_10(st_slave);
        break;




     case 28:
        st2_slave.ipts[5] = myid*100 + 11;
        st2_slave.tavgindx[11*12+5] = 1;
        st2_slave.tavgindx[11*12+6] = 11;
        s_clinic_tavg_10(st_slave);
//         call accumulate_tavg_field(p5*(uvel(:,:,k,curtime,iblock)**2 + &
//                                        vvel(:,:,k,curtime,iblock)**2),tavg_KE_2,iblock,k)

        st2_slave.ipts[5] = myid*100 + 12;
        st2_slave.tavgindx[12*12+5] = 1;
        st2_slave.tavgindx[12*12+6] = 12;
        s_clinic_tavg_12(st_slave);
//         call accumulate_tavg_field(uvel(:,:,k,curtime,iblock)* &
//                                    vvel(:,:,k,curtime,iblock), tavg_UV,iblock,k)
        break;





     case 29:
        st2_slave.ipts[5] = myid*100 + 13;
        st2_slave.tavgindx[13*12+5] = 8;
        st2_slave.tavgindx[13*12+6] = 13;
        s_clinic_tavg_13(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), tavg_TEMP,iblock,k)
        break;

     case 30:
        st2_slave.ipts[5] = myid*100 + 13;
        st2_slave.tavgindx[13*12+5] = 8;
        st2_slave.tavgindx[13*12+6] = 113;
        s_clinic_tavg_13(st_slave);
        break;

     case 31:
        st2_slave.ipts[5] = myid*100 + 13;
        st2_slave.tavgindx[13*12+5] = 8;
        st2_slave.tavgindx[13*12+6] = 213;
        s_clinic_tavg_13(st_slave);
        break;

     case 32:
        st2_slave.ipts[5] = myid*100 + 13;
        st2_slave.tavgindx[13*12+5] = 8;
        st2_slave.tavgindx[13*12+6] = 313;
        s_clinic_tavg_13(st_slave);
        break;

     case 33:
        st2_slave.ipts[5] = myid*100 + 13;
        st2_slave.tavgindx[13*12+5] = 8;
        st2_slave.tavgindx[13*12+6] = 413;
        s_clinic_tavg_13(st_slave);
        break;

     case 34:
        st2_slave.ipts[5] = myid*100 + 13;
        st2_slave.tavgindx[13*12+5] = 8;
        st2_slave.tavgindx[13*12+6] = 513;
        s_clinic_tavg_13(st_slave);
        break;

     case 35:
        st2_slave.ipts[5] = myid*100 + 13;
        st2_slave.tavgindx[13*12+5] = 8;
        st2_slave.tavgindx[13*12+6] = 613;
        s_clinic_tavg_13(st_slave);
        break;

     case 36:
        st2_slave.ipts[5] = myid*100 + 13;
        st2_slave.tavgindx[13*12+5] = 8;
        st2_slave.tavgindx[13*12+6] = 713;
        s_clinic_tavg_13(st_slave);
        break;




     case 37:
        st2_slave.ipts[5] = myid*100 + 14;
        st2_slave.tavgindx[14*12+5] = 1;
        st2_slave.tavgindx[14*12+6] = 14;
        s_clinic_tavg_13(st_slave);
//        call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), tavg_TEMP_2,iblock,k)

        st2_slave.ipts[5] = myid*100 + 15;
        st2_slave.tavgindx[15*12+5] = 1;
        st2_slave.tavgindx[15*12+6] = 15;
        s_clinic_tavg_13(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), tavg_TEMP_MAX,iblock,k)
        break;


     case 38:
        st2_slave.ipts[5] = myid*100 + 16;
        st2_slave.tavgindx[16*12+5] = 1;
        st2_slave.tavgindx[16*12+6] = 16;
        s_clinic_tavg_13(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), tavg_TEMP_MIN,iblock,k)

        st2_slave.ipts[5] = myid*100 + 17;
        st2_slave.tavgindx[17*12+5] = 1;
        st2_slave.tavgindx[17*12+6] = 17;
        s_clinic_tavg_17(st_slave);
//         call accumulate_tavg_field(max(TRACER(:,:,k,1,curtime,iblock) - TRACER(:,:,k,1,oldtime,iblock), c0), tavg_dTEMP_POS_3D,iblock,k)
        break;



     case 39:
        st2_slave.ipts[5] = myid*100 + 18;
        st2_slave.tavgindx[18*12+5] = 1;
        st2_slave.tavgindx[18*12+6] = 18;
        s_clinic_tavg_18(st_slave);
//         call accumulate_tavg_field(max(TRACER(:,:,k,1,curtime,iblock) - TRACER(:,:,k,1,oldtime,iblock), c0), tavg_dTEMP_POS_2D,iblock,1)

        st2_slave.ipts[5] = myid*100 + 19;
        st2_slave.tavgindx[19*12+5] = 1;
        st2_slave.tavgindx[19*12+6] = 19;
        s_clinic_tavg_19(st_slave);
//         call accumulate_tavg_field(min(TRACER(:,:,k,1,curtime,iblock) - TRACER(:,:,k,1,oldtime,iblock), c0), tavg_dTEMP_NEG_3D,iblock,k)
        break;

     case 40:
        st2_slave.ipts[5] = myid*100 + 20;
        st2_slave.tavgindx[20*12+5] = 1;
        st2_slave.tavgindx[20*12+6] = 20;
        s_clinic_tavg_20(st_slave);
//         call accumulate_tavg_field(min(TRACER(:,:,k,1,curtime,iblock) - TRACER(:,:,k,1,oldtime,iblock), c0), tavg_dTEMP_NEG_2D,iblock,1)

        st2_slave.ipts[5] = myid*100 + 21;
        st2_slave.tavgindx[21*12+5] = 1;
        st2_slave.tavgindx[21*12+6] = 21;
        s_clinic_tavg_13(st_slave);
//         if (k <= 8) call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), tavg_T1_8,iblock,k)
        break;



     case 41:
        st2_slave.ipts[5] = myid*100 + 22;
        st2_slave.tavgindx[22*12+5] = 1;
        st2_slave.tavgindx[22*12+6] = 22;
        s_clinic_tavg_22(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)**2, tavg_TEMP2,iblock,k)
        break;



  
     case 42:
        st2_slave.ipts[5] = myid*100 + 23;
        st2_slave.tavgindx[23*12+5] = 8;
        st2_slave.tavgindx[23*12+6] = 23;
        s_clinic_tavg_13(st_slave);
// if(k==1)
//   accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock),tavg_SST,iblock,k)
        break;

     case 43:
        st2_slave.ipts[5] = myid*100 + 23;
        st2_slave.tavgindx[23*12+5] = 8;
        st2_slave.tavgindx[23*12+6] = 123;
        s_clinic_tavg_13(st_slave);
        break;

     case 44:
        st2_slave.ipts[5] = myid*100 + 23;
        st2_slave.tavgindx[23*12+5] = 8;
        st2_slave.tavgindx[23*12+6] = 223;
        s_clinic_tavg_13(st_slave);
        break;

     case 45:
        st2_slave.ipts[5] = myid*100 + 23;
        st2_slave.tavgindx[23*12+5] = 8;
        st2_slave.tavgindx[23*12+6] = 323;
        s_clinic_tavg_13(st_slave);
        break;

     case 46:
        st2_slave.ipts[5] = myid*100 + 23;
        st2_slave.tavgindx[23*12+5] = 8;
        st2_slave.tavgindx[23*12+6] = 423;
        s_clinic_tavg_13(st_slave);
        break;

     case 47:
        st2_slave.ipts[5] = myid*100 + 23;
        st2_slave.tavgindx[23*12+5] = 8;
        st2_slave.tavgindx[23*12+6] = 523;
        s_clinic_tavg_13(st_slave);
        break;

     case 48:
        st2_slave.ipts[5] = myid*100 + 23;
        st2_slave.tavgindx[23*12+5] = 8;
        st2_slave.tavgindx[23*12+6] = 623;
        s_clinic_tavg_13(st_slave);
        break;

     case 49:
        st2_slave.ipts[5] = myid*100 + 23;
        st2_slave.tavgindx[23*12+5] = 8;
        st2_slave.tavgindx[23*12+6] = 723;
        s_clinic_tavg_13(st_slave);
        break;





     case 50:
        st2_slave.ipts[5] = myid*100 + 24;
        st2_slave.tavgindx[24*12+5] = 1;
        st2_slave.tavgindx[24*12+6] = 24;
        s_clinic_tavg_22(st_slave);
// if(k==1)
//   accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)**2,tavg_SST2,iblock,k)
        break;



     case 51:
        st2_slave.ipts[5] = myid*100 + 25;
        st2_slave.tavgindx[25*12+5] = 8;
        st2_slave.tavgindx[25*12+6] = 25;
        s_clinic_tavg_25(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), tavg_SALT,iblock,k)
        break;

     case 52:
        st2_slave.ipts[5] = myid*100 + 25;
        st2_slave.tavgindx[25*12+5] = 8;
        st2_slave.tavgindx[25*12+6] = 125;
        s_clinic_tavg_25(st_slave);
        break;

     case 53:
        st2_slave.ipts[5] = myid*100 + 25;
        st2_slave.tavgindx[25*12+5] = 8;
        st2_slave.tavgindx[25*12+6] = 225;
        s_clinic_tavg_25(st_slave);
        break;

     case 54:
        st2_slave.ipts[5] = myid*100 + 25;
        st2_slave.tavgindx[25*12+5] = 8;
        st2_slave.tavgindx[25*12+6] = 325;
        s_clinic_tavg_25(st_slave);
        break;

     case 55:
        st2_slave.ipts[5] = myid*100 + 25;
        st2_slave.tavgindx[25*12+5] = 8;
        st2_slave.tavgindx[25*12+6] = 425;
        s_clinic_tavg_25(st_slave);
        break;

     case 56:
        st2_slave.ipts[5] = myid*100 + 25;
        st2_slave.tavgindx[25*12+5] = 8;
        st2_slave.tavgindx[25*12+6] = 525;
        s_clinic_tavg_25(st_slave);
        break;

     case 57:
        st2_slave.ipts[5] = myid*100 + 25;
        st2_slave.tavgindx[25*12+5] = 8;
        st2_slave.tavgindx[25*12+6] = 625;
        s_clinic_tavg_25(st_slave);
        break;

     case 58:
        st2_slave.ipts[5] = myid*100 + 25;
        st2_slave.tavgindx[25*12+5] = 8;
        st2_slave.tavgindx[25*12+6] = 725;
        s_clinic_tavg_25(st_slave);
        break;





     case 59:
        st2_slave.ipts[5] = myid*100 + 26;
        st2_slave.tavgindx[26*12+5] = 1;
        st2_slave.tavgindx[26*12+6] = 26;
        s_clinic_tavg_25(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), tavg_SALT_2,iblock,k)

        st2_slave.ipts[5] = myid*100 + 27;
        st2_slave.tavgindx[27*12+5] = 1;
        st2_slave.tavgindx[27*12+6] = 27;
        s_clinic_tavg_25(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), tavg_SALT_MAX,iblock,k)
        break;


     case 60:
        st2_slave.ipts[5] = myid*100 + 28;
        st2_slave.tavgindx[28*12+5] = 1;
        st2_slave.tavgindx[28*12+6] = 28;
        s_clinic_tavg_25(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), tavg_SALT_MIN,iblock,k)

        st2_slave.ipts[5] = myid*100 + 29;
        st2_slave.tavgindx[29*12+5] = 1;
        st2_slave.tavgindx[29*12+6] = 29;
        s_clinic_tavg_25(st_slave);
//         if (k <= 8) call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), tavg_S1_8,iblock,k)
        break;


     case 61:
        st2_slave.ipts[5] = myid*100+ 30;
        st2_slave.tavgindx[30*12+5] = 1;
        st2_slave.tavgindx[30*12+6] = 30;
        s_clinic_tavg_30(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock)**2, tavg_SALT2,iblock,k)

        st2_slave.ipts[5] = myid*100 + 31;
        st2_slave.tavgindx[31*12+5] = 1;
        st2_slave.tavgindx[31*12+6] = 31;
        s_clinic_tavg_31(st_slave);
//         call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)* TRACER(:,:,k,2,curtime,iblock), tavg_ST,iblock,k)
        break;


     case 62:
        st2_slave.ipts[5] = myid*100 + 32;
        st2_slave.tavgindx[32*12+5] = 1;
        st2_slave.tavgindx[32*12+6] = 32;
        s_clinic_tavg_32(st_slave);
//         call accumulate_tavg_field(RHO(:,:,k,curtime,iblock), tavg_RHO,iblock,k)

        st2_slave.ipts[5] = myid*100 + 33;
        st2_slave.tavgindx[33*12+5] = 1;
        st2_slave.tavgindx[33*12+6] = 33;
        s_clinic_tavg_33(st_slave);
        break;


     case 63:
        st2_slave.ipts[5] = myid*100 + 34;
        st2_slave.tavgindx[34*12+5] = 1;
        st2_slave.tavgindx[34*12+6] = 34;
        s_clinic_tavg_34(st_slave);

        st2_slave.ipts[5] = myid*100 + 35;
        st2_slave.tavgindx[35*12+5] = 1;
        st2_slave.tavgindx[35*12+6] = 35;
        s_clinic_tavg_35(st_slave);
        break;



     default:
        break;
      
     }

 return;
 }
 // end of s_clinic_tavg.









//---------------:---------------------------------------------------
void s_clinic_tavg_0(struct param_tavg_accum *st2_slave)
 {
	
// BOP:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef TAVG_R8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif



 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(POPTESTJN)
 volatile int num_calls=0;
#endif

// Common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// Input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, nsc, part_bott_cells;
 int *KMT, *CALCT;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *uvel, *vvel;
 double *TRACER1cur, *TRACER1old, *TRACER2cur;
 double *RHO, *PSURF, *DZT, *DH; 

// Output variables:

// Input & Output variables:
#ifdef TAVG_R8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// Slave & Temp variables:
 int *s_KMT, *s_CALCT, *s_flag; 
 int i,j,k,l,ni,nj,nj1,pij;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_uvel, *s_vvel, *s_dzkm;
#ifdef TAVG_R8
 double *m_bufloc, *s_bufloc;
#else
 float *m_bufloc, *s_bufloc;
#endif


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

/*
 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/

 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];

 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;



 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 nsc      = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef TAVG_R8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 uvel = st2_slave->darr[2];
 vvel = st2_slave->darr[3];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];


//for(k=0;k<km;k++)
//  s_dzkm[k] = st2_slave->dpts[4+k];



for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/nsc;
  if(myid < nxy%nsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%nsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
      }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;

    s_uvel = locdbl;
    s_vvel = locdbl + usize;
    s_bufloc = locbuf;
//  s_dzkm = locdbl + 2*usize;
    

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

#if defined(POPINSJN)
 if(myproc==0 && myid==1 && l==1 && num_calls<8)
    {
    cpe_printf("k: %d, j: 40, i: 37, slave_master_value:%f. \n",myk,s_bufloc[0],m_bufloc[0]);
//    while(1);
    }
#endif

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,uvel+k*nxy+task_pos1,s_uvel,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef TAVG_R8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif

      athread_get(PE_MODE,m_bufloc + task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      sync_mem;
      while(get_reply!=2);
 

      accum_tavg_field_s(s_uvel,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc + task_pos1,pgsize,(void*)&put_reply,0,0);
      sync_mem;
      while(put_reply!=1);
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)

 return;
 }
 // end of s_clinic_tavg_0.








//---------------:---------------------------------------------------
void s_clinic_tavg_2(struct param_tavg_accum *st2_slave)
 {
	
// BOP:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef TAVG_R8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(POPTESTJN)
 volatile int num_calls=0;
#endif

// Common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// Input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *KMT, *CALCT;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *uvel, *vvel;
 double *TRACER1cur, *TRACER1old, *TRACER2cur;
 double *RHO, *PSURF, *DZT, *DH; 

// Output variables:

// Input & Output variables:
#ifdef TAVG_R8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// Slave & Temp variables:
 int *s_KMT, *s_CALCT, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_uvel, *s_vvel, *s_uvel2, *s_dzkm;
#ifdef TAVG_R8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


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


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef TAVG_R8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 uvel = st2_slave->darr[2];
 vvel = st2_slave->darr[3];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];




for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
      }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_uvel = locdbl;
    s_vvel = locdbl + usize;
    s_uvel2= locdbl + 2*usize;
    s_bufloc = locbuf;
//  s_dzkm = locdbl + 2*usize;

    if(tavgdim == 2)
        {
//      call accumulate_tavg_field(uvel(:,:,k,curtime,iblock)**2,tavg_uvel2,iblock,k)
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,uvel+k*nxy+task_pos1,s_uvel,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef TAVG_R8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_uvel2[i] = s_uvel[i]*s_uvel[i];


      accum_tavg_field_s(s_uvel2,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 
 return;
 }
 // end of s_clinic_tavg_2.







//---------------:---------------------------------------------------
void s_clinic_tavg_5(struct param_tavg_accum *st2_slave)
 {
	
// BOP:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef TAVG_R8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(POPTESTJN)
 volatile int num_calls=0;
#endif

// Common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// Input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *KMT, *CALCT;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *uvel, *vvel;
 double *TRACER1cur, *TRACER1old, *TRACER2cur;
 double *RHO, *PSURF, *DZT, *DH; 

// Output variables:

// Input & Output variables:
#ifdef TAVG_R8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// Slave & Temp variables:
 int *s_KMT, *s_CALCT, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_uvel, *s_vvel, *s_dzkm;
#ifdef TAVG_R8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


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


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef TAVG_R8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 uvel = st2_slave->darr[2];
 vvel = st2_slave->darr[3];
    
 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];


for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
      }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
 
    s_uvel = locdbl;
    s_vvel = locdbl + usize;
    s_bufloc = locbuf;
    
    if(tavgdim == 2)
        {
//!  call accumulate_tavg_field(vvel(:,:,curtime,iblock),tavg_vvel,iblock)
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
// ref line 1105 of tavg.F90:
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,vvel+k*nxy+task_pos1,s_vvel,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef TAVG_R8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      sync_mem;
 

      accum_tavg_field_s(s_vvel,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 
 return;
 }
 // end of s_clinic_tavg_5.








//---------------:---------------------------------------------------
void s_clinic_tavg_7(struct param_tavg_accum *st2_slave)
 {
	
// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *uvel, *vvel;
 double *tracer1cur, *tracer1old, *tracer2cur;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_uvel, *s_vvel, *s_vvel2, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 uvel = st2_slave->darr[2];
 vvel = st2_slave->darr[3];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_uvel = locdbl;
    s_vvel = locdbl + usize;
    s_vvel2= locdbl + 3*usize;
    s_bufloc = locbuf;
    // s_dzkm = locdbl + 2*usize;

    if(tavgdim == 2)
        {
//      call accumulate_tavg_field(vvel(:,:,k,curtime,iblock)**2,tavg_vvel2,iblock,k)
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,vvel+k*nxy+task_pos1,s_vvel,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_vvel2[i] = s_vvel[i]*s_vvel[i];


      accum_tavg_field_s(s_vvel2,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 
 return;
 }
 // end of s_clinic_tavg_7.





//---------------:---------------------------------------------------
void s_clinic_tavg_10(struct param_tavg_accum *st2_slave)
 {
	
// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *uvel, *vvel;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_uvel, *s_vvel, *s_uv2, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;

 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 uvel = st2_slave->darr[2];
 vvel = st2_slave->darr[3];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_uvel = locdbl;
    s_vvel = locdbl + usize;
    s_uv2  = locdbl + 2*usize;
    s_bufloc = locbuf;
    // s_dzkm = locdbl + 2*usize;

    if(tavgdim == 2)
        {
//      call accumulate_tavg_field(uvel(:,:,k,curtime,iblock)**2+vvel(:,:,k,curtime,iblock)**2,tavg_ke,iblock,k)
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,uvel+k*nxy+task_pos1,s_uvel,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,vvel+k*nxy+task_pos1,s_vvel,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=3);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_uv2[i] = p5*(s_uvel[i]*s_uvel[i] + s_vvel[i]*s_vvel[i]);


      accum_tavg_field_s(s_uv2,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 
 return;
 }
 // end of s_clinic_tavg_10.




//---------------:---------------------------------------------------
void s_clinic_tavg_12(struct param_tavg_accum *st2_slave)
 {

// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *uvel, *vvel;
 double *tracer1cur, *tracer1old, *tracer2cur;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_uvel, *s_vvel, *s_uv, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 uvel = st2_slave->darr[2];
 vvel = st2_slave->darr[3];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];



for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
      s_uvel = locdbl;
      s_vvel = locdbl + usize;
      s_uv   = locdbl + 2*usize;
      s_bufloc = locbuf;
      // s_dzkm = locdbl + 2*usize;

    if(tavgdim == 2)
        {
//      call accumulate_tavg_field(vvel(:,:,k,curtime,iblock)**2,tavg_vvel2,iblock,k)
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,uvel+k*nxy+task_pos1,s_uvel,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,vvel+k*nxy+task_pos1,s_vvel,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=3);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_uv[i] = s_uvel[i]*s_vvel[i];


      accum_tavg_field_s(s_uv,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 
 }
 // end of s_clinic_tavg_12.




//---------------:---------------------------------------------------
void s_clinic_tavg_13(struct param_tavg_accum *st2_slave)
 {
// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr1cur;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr1cur, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];



for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur = locdbl;
    s_bufloc = locbuf;
    // s_dzkm = locdbl + 2*usize;

    if(tavgdim == 2)
        {
//      call accumulate_tavg_field(vvel(:,:,k,curtime,iblock)**2,tavg_vvel2,iblock,k)
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      sync_mem;
 

      accum_tavg_field_s(s_tr1cur,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 

 }
 // end of s_clinic_tavg_13.




//---------------:---------------------------------------------------
void s_clinic_tavg_17(struct param_tavg_accum *st2_slave)
 {
// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr1cur, *tr1old, *tr1;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr1cur, *s_tr1old, *s_tr1, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];
 tr1old = st2_slave->darr[5];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur = locdbl;
    s_tr1old = locdbl + usize;
    s_tr1    = locdbl + 2*usize;
    s_bufloc = locbuf;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,tr1old+k*nxy+task_pos1,s_tr1old,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=3);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_tr1[i] = locmax(s_tr1cur[i]-s_tr1old[i],c0);



      accum_tavg_field_s(s_tr1,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 

 }
 // end of s_clinic_tavg_17.




//---------------:---------------------------------------------------
void s_clinic_tavg_18(struct param_tavg_accum *st2_slave)
 {
// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr1cur, *tr1old, *tr1;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr1cur, *s_tr1old, *s_tr1, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];
 tr1old = st2_slave->darr[5];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = 1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur = locdbl;
    s_tr1old = locdbl + usize;
    s_tr1    = locdbl + 2*usize;
    s_bufloc = locbuf;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1old,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=3);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_tr1[i] = locmax(s_tr1cur[i]-s_tr1old[i],c0);


      accum_tavg_field_s(s_tr1,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 

 }
 // end of s_clinic_tavg_18.




//---------------:---------------------------------------------------
void s_clinic_tavg_19(struct param_tavg_accum *st2_slave)
 {
// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr1cur, *tr1old, *tr1;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr1cur, *s_tr1old, *s_tr1, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];
 tr1old = st2_slave->darr[5];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur = locdbl;
    s_tr1old = locdbl + usize;
    s_tr1    = locdbl + 2*usize;
    s_bufloc = locbuf;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,tr1old+k*nxy+task_pos1,s_tr1old,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=3);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_tr1[i] = locmin(s_tr1cur[i]-s_tr1old[i],c0);


      accum_tavg_field_s(s_tr1,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 

 }
 // end of s_clinic_tavg_19.




//---------------:---------------------------------------------------
void s_clinic_tavg_20(struct param_tavg_accum *st2_slave)
 {
// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr1cur, *tr1old, *tr1;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr1cur, *s_tr1old, *s_tr1, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];
 tr1old = st2_slave->darr[5];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];



for(k=0;k<km;k++)
 {
 myk = 1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur = locdbl;
    s_tr1old = locdbl + usize;
    s_tr1    = locdbl + 2*usize;
    s_bufloc = locbuf;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,tr1old+k*nxy+task_pos1,s_tr1old,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=3);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_tr1[i] = locmin(s_tr1cur[i]-s_tr1old[i],c0);


      accum_tavg_field_s(s_tr1,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 

 }
 // end of s_clinic_tavg_20.





//---------------:---------------------------------------------------
void s_clinic_tavg_22(struct param_tavg_accum *st2_slave)
 {
// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr1cur, *tr1cur2;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr1cur, *s_tr1cur2, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];



for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur  = locdbl;
    s_tr1cur2 = locdbl + usize;
    s_bufloc = locbuf;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_tr1cur2[i] = s_tr1cur[i]*s_tr1cur[i];


      accum_tavg_field_s(s_tr1cur2,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 


 }
 // end of s_clinic_tavg_22.








//---------------:---------------------------------------------------
void s_clinic_tavg_25(struct param_tavg_accum *st2_slave)
 {

// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr2cur, *tr2cur2;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr2cur, *s_tr2cur2, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr2cur = st2_slave->darr[6];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr2cur  = locdbl;
    s_tr2cur2 = locdbl + usize;
    s_bufloc = locbuf;
    // s_dzkm = locdbl + 2*usize;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr2cur+k*nxy+task_pos1,s_tr2cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      sync_mem;
 


      accum_tavg_field_s(s_tr2cur,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 

 }
 // end of s_clinic_tavg_25.




//---------------:---------------------------------------------------
void s_clinic_tavg_30(struct param_tavg_accum *st2_slave)
 {

// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr2cur, *tr2cur2;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr2cur, *s_tr2cur2, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr2cur = st2_slave->darr[6];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr2cur  = locdbl;
    s_tr2cur2 = locdbl + usize;
    s_bufloc = locbuf;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr2cur+k*nxy+task_pos1,s_tr2cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_tr2cur2[i] = s_tr2cur[i]*s_tr2cur[i];


      accum_tavg_field_s(s_tr2cur2,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 

 }
 // end of s_clinic_tavg_30.




//---------------:---------------------------------------------------
void s_clinic_tavg_31(struct param_tavg_accum *st2_slave)
 {

// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *tr1cur, *tr2cur;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_tr1cur, *s_tr2cur, *s_tr12, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];
 tr2cur = st2_slave->darr[6];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur  = locdbl;
    s_tr2cur  = locdbl + usize;
    s_tr12   = locdbl + 2*usize;
    s_bufloc = locbuf;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,tr2cur+k*nxy+task_pos1,s_tr2cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=3);
      sync_mem;
 

      for(i=0;i<task_size1;i++)
         s_tr12[i] = s_tr1cur[i]*s_tr2cur[i];


      accum_tavg_field_s(s_tr12,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)
 

 }
 // end of s_clinic_tavg_31.




//---------------:---------------------------------------------------
void s_clinic_tavg_32(struct param_tavg_accum *st2_slave)
 {

// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply2, put_reply, put_reply2;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_rho, *s_dzkm;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 rho = st2_slave->darr[7];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];





for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_rho  = locdbl;
    s_bufloc = locbuf;

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      get_reply = 0;
      athread_get(PE_MODE,rho+k*nxy+task_pos1,s_rho,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=2);
      sync_mem;
 


      accum_tavg_field_s(s_rho,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)

 }
 // end of s_clinic_tavg_32.







//---------------:---------------------------------------------------
void s_clinic_tavg_33(struct param_tavg_accum *st2_slave)
 {

// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply1, put_reply, put_reply1;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg,dzk, grav, hflux_factor, salinity_factor;
 double *rho, *psurf, *dzt, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_rho, *s_psurf, *s_dzt, *s_dzkm, *s_work;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 rho = st2_slave->darr[7];
 psurf = st2_slave->darr[8];
 dzt = st2_slave->darr[9];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];


 kmt = st2_slave->iarr[0];
 calct = st2_slave->iarr[1];


 s_dzkm = locdbl;
 for(k=0;k<km;k++)
   {
   s_dzkm[k] = st2_slave->dpts[4+k];
   }


for(k=0;k<km;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_rho  = locdbl + km;
    s_psurf  = locdbl + km + usize;
    s_dzt   = locdbl + km + 2*usize;
    s_work   = locdbl + km + 3*usize;
    s_bufloc = locbuf;

    s_kmt = locint + 8;
    s_calct = locint + 8 +usize;
   

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      
      get_reply1 = 0;
      athread_get(PE_MODE,kmt+k*nxy+task_pos1,s_kmt,SIZEINT*task_size1,(void*)&get_reply1,0,0,0);
      while(get_reply1!=1);
      sync_mem;
      
      
      get_reply = 0;
      athread_get(PE_MODE,rho+k*nxy+task_pos1,s_rho,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,psurf+k*nxy+task_pos1,s_psurf,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,dzt+k*nxy+task_pos1,s_dzt,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=4);
      sync_mem;
 

      if(sfc_layer_type == sfc_layer_varthick && myk==1)
        {
        for(i=0;i<task_size1;i++)
          {
          if(myk <= s_kmt[i])
            s_work[i] = s_rho[i] * (s_dzkm[k] + s_psurf[i]/grav);
          else
            s_work[i] = c0;
          }
        }
      else
        {
        if(part_bott_cells > 0)
          {
          for(i=0;i<task_size1;i++)
            {
            if(myk <= s_kmt[i])
              s_work[i] = s_rho[i] * s_dzt[i];
            else
              s_work[i] = c0;
            }
          }
        else
          {
          for(i=0;i<task_size1;i++)
            {
            if(myk <= s_kmt[i])
              s_work[i] = s_rho[i] * s_dzkm[k];
            else
              s_work[i] = c0;
            }
          }
        }

      accum_tavg_field_s(s_work,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)

 }
 // end of s_clinic_tavg_33.






//---------------:---------------------------------------------------
void s_clinic_tavg_34(struct param_tavg_accum *st2_slave)
 {

// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply1, put_reply, put_reply1;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg, grav, hflux_factor, salinity_factor;
 double *tr1cur, *tr2cur, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,factor;
 double *s_tr1cur, *s_tr2cur, *s_dh, *s_work;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*

 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];
 tr2cur = st2_slave->darr[6];
 dh = st2_slave->darr[10];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];


 calct = st2_slave->iarr[1];



// for(k=0;k<km;k++)
for(k=0;k<1;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur  = locdbl;
    s_tr2cur  = locdbl + usize;
    s_dh  = locdbl + 2*usize;
    s_work   = locdbl + 3*usize;
    s_bufloc = locbuf;

    s_calct = locint + 8 + usize;
   

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      
      if(k<1)
        {
        get_reply1 = 0;
        athread_get(PE_MODE,calct+k*nxy+task_pos1,s_calct,SIZEINT*task_size1,(void*)&get_reply1,0,0,0);
        while(get_reply1!=1);
        sync_mem;
        }
      
      
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,tr2cur+k*nxy+task_pos1,s_tr2cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,dh+k*nxy+task_pos1,s_dh,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=4);
      sync_mem;
 

      if(sfc_layer_type != sfc_layer_varthick && myk==1)
        {
        factor = c1/hflux_factor;
        for(i=0;i<task_size1;i++)
          {
          s_work[i] = c0;
          if(s_calct[i] > 0)
            s_work[i] = s_dh[i] * s_tr1cur[k] * factor;
          }
        }

      accum_tavg_field_s(s_work,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)

 }
 // end of s_clinic_tavg_34.






//---------------:---------------------------------------------------
void s_clinic_tavg_35(struct param_tavg_accum *st2_slave)
 {

// bop:
 int locint[SZldmint], *s2di;
 double locdbl[SZldmdbl], *s2dd;
#ifdef tavg_r8
 double locbuf[SZldmbuf];
#else
 float locbuf[SZldmbuf];
#endif


 volatile int get_reply, get_reply1, put_reply, put_reply1;
#if defined(poptestjn)
 volatile int num_calls=0;
#endif

// common variables:
 int myproc,myid,threads_num;
 int task_num0,task_num,
     task_pos,task_pos1, 
     task_size,task_tile,last_size,task_tile1,task_size1;
 int nx,ny,myk,nsc,nxy,km;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 

// input variables: 
 int iblock, niblock, sfc_layer_type, sfc_layer_varthick;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, lnsc, part_bott_cells;
 int *kmt, *calct;			// nx*ny
 double dtavg, grav, hflux_factor, salinity_factor;
 double *tr1cur, *tr2cur, *dh; 

// output variables:

// input & output variables:
#ifdef tavg_r8
 double *tavg_buf2d;		// nx*ny
 double *tavg_buf3d; 
#else
 float *tavg_buf2d;
 float *tavg_buf3d;
#endif


// slave & temp variables:
 int *s_kmt, *s_calct, *s_flag; 
 int i,j,k,l;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,factor;
 double *s_tr1cur, *s_tr2cur, *s_dh, *s_work;
#ifdef tavg_r8
 double  *m_bufloc, *s_bufloc;
#else
 float  *m_bufloc, *s_bufloc;
#endif


// constants && parameters:
 double 
      p25= 0.25,
      p5 = 0.5,
      c0 = 0.0,
      c1 = 1.0;


// eop.


#if defined(poptestjn)
 num_calls = num_calls + 1;
#endif


/*
 for(i=0;i<SZldmint;i++)
    locint[i] = 0;
 for(i=0;i<SZldmbuf;i++)
    locbuf[i] = 0.0f;
 for(i=0;i<SZldmdbl;i++)
    locdbl[i] = c0;
*/



 nx = st2_slave->ipts[0];
 ny = st2_slave->ipts[1];
 sfc_layer_type = st2_slave->ipts[2]; 
 sfc_layer_varthick = st2_slave->ipts[3];
 myproc = st2_slave->ipts[4];
 km     = st2_slave->ipts[6];
 iblock = st2_slave->ipts[7];
 niblock= st2_slave->ipts[8]; 
 part_bott_cells= st2_slave->ipts[9];

 grav = st2_slave->dpts[0];
 hflux_factor    = st2_slave->dpts[1];
 salinity_factor = st2_slave->dpts[2];
 dtavg  = st2_slave->dpts[3];


 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;


 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 lnsc     = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 


#ifdef tavg_r8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 tr1cur = st2_slave->darr[4];
 tr2cur = st2_slave->darr[6];
 dh = st2_slave->darr[10];

 s_flag = locint;
 s_flag[0]=st2_slave->tavgindx[jtmp*12+8];
 s_flag[1]=st2_slave->tavgindx[jtmp*12+9];
 s_flag[2]=st2_slave->tavgindx[jtmp*12+10];
 s_flag[3]=st2_slave->tavgindx[jtmp*12+11];


 calct = st2_slave->iarr[1];



// for(k=0;k<km;k++)
for(k=0;k<1;k++)
 {
 myk = k+1;
 tavgflag = ((s_flag[k/16])>>(k%16))%2;

 if(tavgflag > 0)
  {	
  nxy = nx*ny;
  task_num0 = nxy/lnsc;
  if(myid < nxy%lnsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid;  
    }
  else
    {
    task_num = task_num0;
    task_pos = nxy%lnsc + task_num0*myid; 
    }


  if(task_num > 0)
    {	
    mx_usize = MAXUSZ;  
    if(task_num <= mx_usize)
      {	
      usize  = task_num;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = mx_usize;	
      task_tile = task_num/usize;
      task_size = usize;
      last_size = task_num%usize;
       }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;
    
    s_tr1cur  = locdbl;
    s_tr2cur  = locdbl + usize;
    s_dh  = locdbl + 2*usize;
    s_work   = locdbl + 3*usize;
    s_bufloc = locbuf;

    s_calct = locint + 8 + usize;
   

    if(tavgdim == 2)
        {
        offset2d = ((tavgloc-1)*niblock+(iblock-1))*nxy;
        m_bufloc = tavg_buf2d + offset2d;
        }
    else
        {
        offset3d = (((tavgloc-1)*niblock+(iblock-1))*km+(myk-1))*nxy;
        m_bufloc = tavg_buf3d + offset3d;
        }

    task_pos1 = task_pos;
    for(l=0;l<task_tile1;l++)
      {	 
      if(l<task_tile)
         task_size1 = task_size;
      else  
         task_size1 = last_size;

      // get:  
      
      if(k<1)
        {
        get_reply1 = 0;
        athread_get(PE_MODE,calct+k*nxy+task_pos1,s_calct,SIZEINT*task_size1,(void*)&get_reply1,0,0,0);
        while(get_reply1!=1);
        sync_mem;
        }
      
      
      get_reply = 0;
      athread_get(PE_MODE,tr1cur+k*nxy+task_pos1,s_tr1cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,tr2cur+k*nxy+task_pos1,s_tr2cur,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,dh+k*nxy+task_pos1,s_dh,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef tavg_r8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc+task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      while(get_reply!=4);
      sync_mem;
 

      if(sfc_layer_type != sfc_layer_varthick && myk==1)
        {
        factor = c1/salinity_factor;
        for(i=0;i<task_size1;i++)
          {
          s_work[i] = c0;
          if(s_calct[i] > 0)
            s_work[i] = s_dh[i] * s_tr2cur[k] * factor;
          }
        }

      accum_tavg_field_s(s_work,s_bufloc,task_size1,tavgmeth,dtavg);


      // put:  
      put_reply = 0;
      athread_put(PE_MODE,s_bufloc,m_bufloc+task_pos1,pgsize,(void*)&put_reply,0,0);
      while(put_reply!=1);
      sync_mem;
     
      task_pos1 = task_pos1 + task_size1;
      }
      // end of for(l)

    }
    // end of if(task_num>0)

  }
  // end of if(tavgflag > 0)

}
// end of for(k)

 }
 // end of s_clinic_tavg_35.



/*
         if (sfc_layer_type == sfc_layer_varthick .and. k == 1) then
           where (k <= KMT(:,:,iblock))
             WORK1 = RHO(:,:,k,curtime,iblock) * ( dz(k) + PSURF(:,:,curtime,iblock)/grav )
           elsewhere
             WORK1 = c0
           endwhere
         else
           if (partial_bottom_cells) then
             where (k <= KMT(:,:,iblock))
               WORK1 = RHO(:,:,k,curtime,iblock) * DZT(:,:,k,iblock)
             elsewhere
               WORK1 = c0
             endwhere
           else
             where (k <= KMT(:,:,iblock))
               WORK1 = RHO(:,:,k,curtime,iblock) * dz(k)
             elsewhere
               WORK1 = c0
             endwhere
           endif
         endif 
         call accumulate_tavg_field(WORK1,tavg_RHO_VINT,iblock,k)

         if ( sfc_layer_type /= sfc_layer_varthick .and. k == 1) then
           if (accumulate_tavg_now(tavg_RESID_T)) then
              WORK1 = c0
              factor = c1/hflux_factor  ! converts to W/m^2
              where (CALCT(:,:,iblock))  &
                WORK1=DH(:,:,iblock)*TRACER(:,:,1,1,curtime,iblock)*factor
              call accumulate_tavg_field(WORK1,tavg_RESID_T,iblock,k)
           endif

           if (accumulate_tavg_now(tavg_RESID_S)) then
              WORK1 = c0
              factor = c1/salinity_factor  ! converts to kg(freshwater)/m^2/s
              where (CALCT(:,:,iblock)) &
                WORK1 = DH(:,:,iblock)*TRACER(:,:,k,2,curtime,iblock)*factor
              call accumulate_tavg_field(WORK1,tavg_RESID_S,iblock,k)
           endif
         endif  ! sfc_layer_type

*/


