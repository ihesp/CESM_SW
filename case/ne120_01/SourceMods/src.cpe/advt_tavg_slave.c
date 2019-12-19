#include <stdio.h>
#include <math.h>
#include <slave.h>
#include <unistd.h>
#include "advt_tavg_accum.h"

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



void s_advt_tavg_0();
// call accumulate_tavg_field(WTK,tavg_WVEL,bid,k)

void s_advt_tavg_3();
// call state(k,1,TRCR(:,:,k,1), TRCR(:,:,k,2), this_block, RHOFULL=RHOK1)
// call accumulate_tavg_field(RHOK1,tavg_PD,bid,k)

void s_advt_tavg_17();
// call accumulate_tavg_field(WORK,tavg_UE_TRACER(1),bid,k)

void s_advt_tavg_18();
// call accumulate_tavg_field(WORK,tavg_VN_TRACER(1),bid,k)




//---------------:---------------------------------------------------
void s_advt_tavg(struct param_tavg_accum *st2_master)
 {
	
// BOP:
 volatile int get_reply, put_reply;

// Common variables:
 int myproc,myid,nt,threads_num;

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


 nt = st2_slave.ipts[4];
//  numtavg= st2_slave.ipts[7];

// For current case, only 6 tavg variables are activated:
// 1   tavgid: 129,  tavgname: tavg_WVEL 
// 4   tavgid: 152,  tavgname: tavg_PD
// 18  tavgid: 138,  tavgname: tavg_UE_TRACER(1)
// 19  tavgid: 139,  tavgname: tavg_VN_TRACER(1)
// 22  tavgid: 141,  tavgname: tavg_UE_TRACER(2)
// 23  tavgid: 142,  tavgname: tavg_VN_TRACER(2)

  switch(myid) 
     {
     case 0:
        st2_slave.ipts[5] = myid*100;
        st2_slave.tavgindx[5] = 4;
        st2_slave.tavgindx[6] = 0;
        s_advt_tavg_0(st_slave);
        // call accumulate_tavg_field(WTK,tavg_WVEL,bid,k)
        break;

     case 1:
        st2_slave.ipts[5] = myid*100+0;
        st2_slave.tavgindx[5] = 4;
        st2_slave.tavgindx[6] = 100;
        s_advt_tavg_0(st_slave);
        break;
     
     case 2:
        st2_slave.ipts[5] = myid*100+0;
        st2_slave.tavgindx[5] = 4;
        st2_slave.tavgindx[6] = 200;
        s_advt_tavg_0(st_slave);
        break;
     
     case 3:
        st2_slave.ipts[5] = myid*100+0;
        st2_slave.tavgindx[5] = 4;
        st2_slave.tavgindx[6] = 300;
        s_advt_tavg_0(st_slave);
        break;
//------------------------    

	
     case 4:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 3;
        s_advt_tavg_3(st_slave);
        // call accumulate_tavg_field(RHOK1,tavg_PD,bid,k)		
        break;		

     case 5:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 103;
        s_advt_tavg_3(st_slave);
        break;

     case 6:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 203;
        s_advt_tavg_3(st_slave);
        break;

     case 7:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 303;
        s_advt_tavg_3(st_slave);
        break;

     case 8:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 403;
        s_advt_tavg_3(st_slave);
        break;

     case 9:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 503;
        s_advt_tavg_3(st_slave);
        break;

     case 10:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 603;
        s_advt_tavg_3(st_slave);
        break;

     case 11:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 703;
        s_advt_tavg_3(st_slave);
        break;

     case 12:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 803;
        s_advt_tavg_3(st_slave);
        break;

     case 13:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 903;
        s_advt_tavg_3(st_slave);
        break;

     case 14:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1003;
        s_advt_tavg_3(st_slave);
        break;

     case 15:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1103;
        s_advt_tavg_3(st_slave);
        break;

     case 16:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1203;
        s_advt_tavg_3(st_slave);
        break;

     case 17:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1303;
        s_advt_tavg_3(st_slave);
        break;

     case 18:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1403;
        s_advt_tavg_3(st_slave);
        break;
     case 19:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1503;
        s_advt_tavg_3(st_slave);
        break;

     case 20:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1603;
        s_advt_tavg_3(st_slave);
        break;

     case 21:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1703;
        s_advt_tavg_3(st_slave);
        break;
     case 22:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1803;
        s_advt_tavg_3(st_slave);
        break;

     case 23:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 1903;
        s_advt_tavg_3(st_slave);
        break;

     case 24:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 2003;
        s_advt_tavg_3(st_slave);
        break;
		
     case 25:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 2103;
        s_advt_tavg_3(st_slave);
        break;

     case 26:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 2203;
        s_advt_tavg_3(st_slave);
        break;

     case 27:
        st2_slave.ipts[5] = myid*100+3;
        st2_slave.tavgindx[3*12+5] = 24;
        st2_slave.tavgindx[3*12+6] = 2303;
        s_advt_tavg_3(st_slave);
        break;		
//------------------------    

     case 28:
        st2_slave.ipts[5] = myid*100+17;
        st2_slave.tavgindx[17*12+5] = 8;
        st2_slave.tavgindx[17*12+6] = 17;
        s_advt_tavg_17(st_slave);
// call accumulate_tavg_field(WORK,tavg_UE_TRACER(1),bid,k)
        break;
		
     case 29:
        st2_slave.ipts[5] = myid*100+17;
        st2_slave.tavgindx[17*12+5] = 8;
        st2_slave.tavgindx[17*12+6] = 117;
        s_advt_tavg_17(st_slave);
        break;
		
     case 30:
        st2_slave.ipts[5] = myid*100+17;
        st2_slave.tavgindx[17*12+5] = 8;
        st2_slave.tavgindx[17*12+6] = 217;
        s_advt_tavg_17(st_slave);
        break;

     case 31:
        st2_slave.ipts[5] = myid*100+17;
        st2_slave.tavgindx[17*12+5] = 8;
        st2_slave.tavgindx[17*12+6] = 317;
        s_advt_tavg_17(st_slave);
        break;

     case 32:
        st2_slave.ipts[5] = myid*100+17;
        st2_slave.tavgindx[17*12+5] = 8;
        st2_slave.tavgindx[17*12+6] = 417;
        s_advt_tavg_17(st_slave);
        break;

     case 33:
        st2_slave.ipts[5] = myid*100+17;
        st2_slave.tavgindx[17*12+5] = 8;
        st2_slave.tavgindx[17*12+6] = 517;
        s_advt_tavg_17(st_slave);
        break;

     case 34:
        st2_slave.ipts[5] = myid*100+17;
        st2_slave.tavgindx[17*12+5] = 8;
        st2_slave.tavgindx[17*12+6] = 617;
        s_advt_tavg_17(st_slave);
        break;

     case 35:
        st2_slave.ipts[5] = myid*100+17;
        st2_slave.tavgindx[17*12+5] = 8;
        st2_slave.tavgindx[17*12+6] = 717;
        s_advt_tavg_17(st_slave);
        break;
//------------------------    



     case 36:
        st2_slave.ipts[5] = myid*100+18;
        st2_slave.tavgindx[18*12+5] = 8;
        st2_slave.tavgindx[18*12+6] = 18;
        s_advt_tavg_18(st_slave);
// call accumulate_tavg_field(WORK,tavg_VN_TRACER(1),bid,k)
        break;

     case 37:
        st2_slave.ipts[5] = myid*100+18;
        st2_slave.tavgindx[18*12+5] = 8;
        st2_slave.tavgindx[18*12+6] = 118;
        s_advt_tavg_18(st_slave);
        break;

     case 38:
        st2_slave.ipts[5] = myid*100+18;
        st2_slave.tavgindx[18*12+5] = 8;
        st2_slave.tavgindx[18*12+6] = 218;
        s_advt_tavg_18(st_slave);
        break;

     case 39:
        st2_slave.ipts[5] = myid*100+18;
        st2_slave.tavgindx[18*12+5] = 8;
        st2_slave.tavgindx[18*12+6] = 318;
        s_advt_tavg_18(st_slave);
        break;		

     case 40:
        st2_slave.ipts[5] = myid*100+18;
        st2_slave.tavgindx[18*12+5] = 8;
        st2_slave.tavgindx[18*12+6] = 418;
        s_advt_tavg_18(st_slave);
        break;

     case 41:
        st2_slave.ipts[5] = myid*100+18;
        st2_slave.tavgindx[18*12+5] = 8;
        st2_slave.tavgindx[18*12+6] = 518;
        s_advt_tavg_18(st_slave);
        break;		

     case 42:
        st2_slave.ipts[5] = myid*100+18;
        st2_slave.tavgindx[18*12+5] = 8;
        st2_slave.tavgindx[18*12+6] = 618;
        s_advt_tavg_18(st_slave);
        break;

     case 43:
        st2_slave.ipts[5] = myid*100+18;
        st2_slave.tavgindx[18*12+5] = 8;
        st2_slave.tavgindx[18*12+6] = 718;
        s_advt_tavg_18(st_slave);
        break;				
//------------------------    



     case 44:
        st2_slave.ipts[4] = 100+nt;	   // n=1 
        st2_slave.ipts[5] = myid*100+21;
        st2_slave.tavgindx[21*12+5] = 8;
        st2_slave.tavgindx[21*12+6] = 21;
        s_advt_tavg_17(st_slave);
// call accumulate_tavg_field(WORK,tavg_UE_TRACER(2),bid,k)
        break;

     case 45:
        st2_slave.ipts[4] = 100+nt;	
        st2_slave.ipts[5] = myid*100+21;
        st2_slave.tavgindx[21*12+5] = 8;
        st2_slave.tavgindx[21*12+6] = 121;
        s_advt_tavg_17(st_slave);
        break;				

     case 46:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+21;
        st2_slave.tavgindx[21*12+5] = 8;
        st2_slave.tavgindx[21*12+6] = 221;
        s_advt_tavg_17(st_slave);
        break;

     case 47:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+21;
        st2_slave.tavgindx[21*12+5] = 8;
        st2_slave.tavgindx[21*12+6] = 321;
        s_advt_tavg_17(st_slave);
        break;				

     case 48:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+21;
        st2_slave.tavgindx[21*12+5] = 8;
        st2_slave.tavgindx[21*12+6] = 421;
        s_advt_tavg_17(st_slave);
        break;

     case 49:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+21;
        st2_slave.tavgindx[21*12+5] = 8;
        st2_slave.tavgindx[21*12+6] = 521;
        s_advt_tavg_17(st_slave);
        break;				

     case 50:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+21;
        st2_slave.tavgindx[21*12+5] = 8;
        st2_slave.tavgindx[21*12+6] = 621;
        s_advt_tavg_17(st_slave);
        break;

     case 51:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+21;
        st2_slave.tavgindx[21*12+5] = 8;
        st2_slave.tavgindx[21*12+6] = 721;
        s_advt_tavg_17(st_slave);
        break;				
//------------------------    



     case 52:
        st2_slave.ipts[4] = 100+nt;	   // n=1 	 
        st2_slave.ipts[5] = myid*100+22;
        st2_slave.tavgindx[22*12+5] = 8;
        st2_slave.tavgindx[22*12+6] = 22;
        s_advt_tavg_18(st_slave);
// call accumulate_tavg_field(WORK,tavg_VN_TRACER(2),bid,k)
        break;

     case 53:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+22;
        st2_slave.tavgindx[22*12+5] = 8;
        st2_slave.tavgindx[22*12+6] = 122;
        s_advt_tavg_18(st_slave);
        break;				

     case 54:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+22;
        st2_slave.tavgindx[22*12+5] = 8;
        st2_slave.tavgindx[22*12+6] = 222;
        s_advt_tavg_18(st_slave);
        break;

     case 55:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+22;
        st2_slave.tavgindx[22*12+5] = 8;
        st2_slave.tavgindx[22*12+6] = 322;
        s_advt_tavg_18(st_slave);
        break;	
		
     case 56:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+22;
        st2_slave.tavgindx[22*12+5] = 8;
        st2_slave.tavgindx[22*12+6] = 422;
        s_advt_tavg_18(st_slave);
        break;

     case 57:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+22;
        st2_slave.tavgindx[22*12+5] = 8;
        st2_slave.tavgindx[22*12+6] = 522;
        s_advt_tavg_18(st_slave);
        break;				

     case 58:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+22;
        st2_slave.tavgindx[22*12+5] = 8;
        st2_slave.tavgindx[22*12+6] = 622;
        s_advt_tavg_18(st_slave);
        break;

     case 59:
        st2_slave.ipts[4] = 100+nt;		 
        st2_slave.ipts[5] = myid*100+22;
        st2_slave.tavgindx[22*12+5] = 8;
        st2_slave.tavgindx[22*12+6] = 722;
        s_advt_tavg_18(st_slave);
        break;				


     default:
        break;
      
     }

 return;
 }
 // end of s_advt_tavg.









//---------------:---------------------------------------------------
void s_advt_tavg_0(struct param_tavg_accum *st2_slave)
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
 int nx,ny,myk,nxy,km,nt;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 
// call accumulate_tavg_field(WTK,tavg_WVEL,bid,k)
// Input variables: 
 int iblock, niblock;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, nsc, part_bott_cells;
 double dtavg;
 double *wtk3d;

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
 int i,j,k,l,ni,nj,nj1,pij;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_wtk;
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
 myproc = st2_slave->ipts[2];
 km     = st2_slave->ipts[3];
 nt     = st2_slave->ipts[4];
 iblock = (st2_slave->ipts[6])/100; 
 niblock = (st2_slave->ipts[6])%100;  
 myk     = st2_slave->ipts[15];
 

 dtavg  = st2_slave->dpts[0];

 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;

 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgflag = st2_slave->tavgindx[jtmp*12+1];
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

 wtk3d = st2_slave->darr[4];
 

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

    s_wtk = locdbl;
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
      athread_get(PE_MODE,wtk3d+task_pos1,s_wtk,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
#ifdef TAVG_R8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif

      athread_get(PE_MODE,m_bufloc + task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);
      sync_mem;
      while(get_reply!=2);
 

      accum_tavg_field_s(s_wtk,s_bufloc,task_size1,tavgmeth,dtavg);


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


 return;
 }
 // end of s_advt_tavg_0.






//---------------:---------------------------------------------------
void s_advt_tavg_3(struct param_tavg_accum *st2_slave)
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
 int nx,ny,myk,nxy,km,nt;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 
// call accumulate_tavg_field(RHOK1,tavg_PD,bid,k)
// Input variables: 
 int iblock, niblock;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, nsc, part_bott_cells;
 double dtavg, *trcr;

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
 int i,j,k,l,ni,nj,nj1,pij;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,ztmp;
 double *s_rhok1,*s_trcr;
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
 myproc = st2_slave->ipts[2];
 km     = st2_slave->ipts[3];
 nt     = st2_slave->ipts[4];
 iblock = (st2_slave->ipts[6])/100; 
 niblock = (st2_slave->ipts[6])%100;  
 myk     = st2_slave->ipts[15];
 
 dtavg  = st2_slave->dpts[0];

 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;

 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgflag = st2_slave->tavgindx[jtmp*12+1];
 tavgloc  = st2_slave->tavgindx[jtmp*12+2];
 tavgdim  = st2_slave->tavgindx[jtmp*12+3];
 tavgmeth = st2_slave->tavgindx[jtmp*12+4];
 nsc      = st2_slave->tavgindx[jtmp*12+5]; 
 myid     = (st2_slave->tavgindx[jtmp*12+6])/100; 

 st2_slave->ipts[14] = nsc;

#ifdef TAVG_R8
 tavg_buf2d = st2_slave->darr[0];
 tavg_buf3d = st2_slave->darr[1];
#else
 tavg_buf2d = st2_slave->farr[0];
 tavg_buf3d = st2_slave->farr[1];
#endif

 trcr  = st2_slave->darr[7];



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

    s_trcr = locdbl;
    s_rhok1 = locdbl + 3*usize;

    s_bufloc = locbuf;

    st2_slave->darr[6] = s_trcr;
    st2_slave->darr[8] = s_rhok1;
    

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
      get_reply2= 0;

      athread_get(PE_MODE,trcr+(myk-1)*nxy+task_pos1,s_trcr,2*SIZEDBL*task_size1,(void*)&get_reply2,0,SIZEDBL*(nxy*km-task_size1),SIZEDBL*task_size1);

#ifdef TAVG_R8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif
      athread_get(PE_MODE,m_bufloc + task_pos1,s_bufloc,pgsize,(void*)&get_reply,0,0,0);

      sync_mem; 
      while(get_reply2!=1);

      st2_slave->ipts[12] = task_size1;
      advt_state_s(st2_slave);


      sync_mem;
      while(get_reply!=1);
 

      accum_tavg_field_s(s_rhok1,s_bufloc,task_size1,tavgmeth,dtavg);


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


 return;
 }
 // end of s_advt_tavg_3.





//---------------:---------------------------------------------------
void s_advt_tavg_17(struct param_tavg_accum *st2_slave)
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
 int nx,ny,myk,nxy,km,nt;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 
// Input variables: 
 int iblock, niblock;
 int tadvt_cent, tadvtype;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, nsc, part_bott_cells;
 double dtavg;
 double *tarea_r;
 double *ute3d, *trcr, *trcr_e, *dzt; 

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
 int i,j,k,l,ni,n;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,*sztmp;
 double *s_tarea_r, *s_ute, *s_trcr, *s_trcr_e, *s_dzt;
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
      c1 = 1.0,
      c2 = 2.0;


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
 myproc = st2_slave->ipts[2];
 km     = st2_slave->ipts[3];
 nt = (st2_slave->ipts[4])%100;
 n  = (st2_slave->ipts[4])/100;
 iblock = (st2_slave->ipts[6])/100; 
 niblock = (st2_slave->ipts[6])%100;  
 myk     = st2_slave->ipts[15];

 part_bott_cells= st2_slave->ipts[8];

 tadvt_cent = st2_slave->ipts[9];
 tadvtype   = st2_slave->ipts[10+n];

 dtavg  = st2_slave->dpts[0];

 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;

 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgflag = st2_slave->tavgindx[jtmp*12+1];
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

 ute3d = st2_slave->darr[3];
 trcr  = st2_slave->darr[7];
 trcr_e= st2_slave->darr[5];
 dzt   = st2_slave->darr[9];
 tarea_r = st2_slave->darr[10];


 if(tavgflag > 0)
  {	
  task_num0 = ny/nsc;
  if(myid < ny%nsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid*nx;  
    }
  else
    {
    task_num = task_num0;
    task_pos = (ny%nsc + task_num0*myid)*nx; 
    }


  if(task_num > 0)
    {	
    nxy = nx*ny;
    mx_usize = MAXUSZ;  
    if((task_num*nx) <= mx_usize)
      {	
      usize  = task_num*nx;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = (mx_usize/nx)*nx;	
      task_tile = task_num/(mx_usize/nx);
      task_size = usize;
      last_size = (task_num%(mx_usize/nx))*nx;
      }
 
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;

    s_ute  = locdbl;
    s_trcr = locdbl + usize;
    s_dzt  = locdbl + 3*usize;	
    s_tarea_r = locdbl + 4*usize;		
    sztmp  = locdbl + 5*usize;
    if(tadvtype != tadvt_cent)
      s_trcr_e  = locdbl + 2*usize;	
	
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
      get_reply2= 0;	  
	  
      athread_get(PE_MODE,ute3d+task_pos1,s_ute,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,dzt+task_pos1,s_dzt,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
	  
      athread_get(PE_MODE,tarea_r+task_pos1,s_tarea_r,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);	  
#ifdef TAVG_R8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif

      athread_get(PE_MODE,m_bufloc + task_pos1,s_bufloc,pgsize,(void*)&get_reply2,0,0,0);
	  
      while(get_reply!=3);	  



      if(tadvtype == tadvt_cent)
        {
        athread_get(PE_MODE,trcr+(n*km+myk-1)*nxy+task_pos1,s_trcr,SIZEDBL*task_size1,(void*)&get_reply2,0,0,0);
        sync_mem; 
        while(get_reply2!=2);
        }
      else
        {
        athread_get(PE_MODE,trcr+(n*km+myk-1)*nxy+task_pos1,s_trcr,SIZEDBL*task_size1,(void*)&get_reply2,0,0,0);
        athread_get(PE_MODE,trcr_e+n*nxy+task_pos1,s_trcr_e,SIZEDBL*task_size1,(void*)&get_reply2,0,0,0);
        sync_mem; 
        while(get_reply2!=3);
        }
 
 
/*
      if (partial_bottom_cells) then
         FUE =  p5*UTE*TAREA_R(:,:,bid)/DZT(:,:,k,bid)
      else
         FUE =  p5*UTE*TAREA_R(:,:,bid)
      endif

      do n=1,nt
         if (tadvect_itype(n) == tadvect_centered) then
            WORK = FUE*(TRCR(:,:,k,n) + eoshift(TRCR(:,:,k,n),dim=1,shift=1))  // by row.
         else
            WORK = c2*FUE*TRACER_E(:,:,n)
         endif
         call accumulate_tavg_field(WORK,tavg_UE_TRACER(n),bid,k)
      enddo

where
     TRCR(nx,ny,km,nt)
     TRACER_E(nx,ny,nt) 
*/
 
      if(part_bott_cells>0)
        { 
        for(i=0;i<task_size1;i++)
          sztmp[i] = p5*s_ute[i]*s_tarea_r[i]/s_dzt[i];
	}
      else
	{ 
        for(i=0;i<task_size1;i++)
	  sztmp[i] = p5*s_ute[i]*s_tarea_r[i];
        }
		 
      if(tadvtype == tadvt_cent)
        { 
        for(i=0;i<task_size1;i++)
          {
          jtmp = i%nx;
          if(jtmp<(nx-1))
            sztmp[i] = sztmp[i]*(s_trcr[i]+s_trcr[i+1]);
          else
            sztmp[i] = sztmp[i]*s_trcr[i];     // the last row is set to 0.
          }
        }
      else
	{ 
        for(i=0;i<task_size1;i++)
	  sztmp[i] = c2*sztmp[i]*s_trcr_e[i];
        }
 

      accum_tavg_field_s(sztmp,s_bufloc,task_size1,tavgmeth,dtavg);


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


 return;
 }
 // end of s_advt_tavg_17.









//---------------:---------------------------------------------------
void s_advt_tavg_18(struct param_tavg_accum *st2_slave)
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
 int nx,ny,myk,nxy,km,nt;

 int offset,offset1,offset2d,offset3d;
 int mx_usize,usize0,usize,usize1,pgsize;
 
// Input variables: 
 int iblock, niblock;
 int tadvt_cent, tadvtype;
 int tavgid, tavgflag, tavgloc, tavgdim, tavgmeth, nsc, part_bott_cells;
 double dtavg;
 double *tarea_r;
 double *vtn3d, *trcr, *trcr_n, *dzt; 

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
 int *s_flag; 
 int i,j,k,l,ni,n;
 int itmp,jtmp,ktmp;
 double xtmp,ytmp,*sztmp;
 double *s_tarea_r, *s_vtn, *s_trcr, *s_trcr_n, *s_dzt;
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
      c1 = 1.0,
      c2 = 2.0;


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
 myproc = st2_slave->ipts[2];
 km     = st2_slave->ipts[3];
 nt = (st2_slave->ipts[4])%100;
 n  = (st2_slave->ipts[4])/100;
 iblock = (st2_slave->ipts[6])/100; 
 niblock = (st2_slave->ipts[6])%100;  
 myk     = st2_slave->ipts[15];

 part_bott_cells= st2_slave->ipts[8];

 tadvt_cent = st2_slave->ipts[9];
 tadvtype   = st2_slave->ipts[10+n];

 dtavg  = st2_slave->dpts[0];

 itmp = st2_slave->ipts[5];
 jtmp = itmp%100;

 tavgid   = st2_slave->tavgindx[jtmp*12+0];
 tavgflag = st2_slave->tavgindx[jtmp*12+1];
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

 vtn3d = st2_slave->darr[2];
 trcr  = st2_slave->darr[7];
 trcr_n= st2_slave->darr[6];
 dzt   = st2_slave->darr[9];
 tarea_r = st2_slave->darr[10];


 if(tavgflag > 0)
  {	
  task_num0 = ny/nsc;
  if(myid < ny%nsc)
    {
    task_num = task_num0 + 1;
    task_pos = task_num*myid*nx;  
    }
  else
    {
    task_num = task_num0;
    task_pos = (ny%nsc + task_num0*myid)*nx; 
    }


  if(task_num > 0)
    {	
    nxy = nx*ny;
    mx_usize = MAXUSZ;  
    if((task_num*nx) <= mx_usize)
      {	
      usize  = task_num*nx;
      task_tile = 0;
      last_size = usize;
      }
    else                 // general case!
      {
      usize  = (mx_usize/nx)*nx;	
      task_tile = task_num/(mx_usize/nx);
      task_size = usize;
      last_size = (task_num%(mx_usize/nx))*nx;
      }
  
    if(last_size>0)
      task_tile1 = task_tile + 1;
    else 
      task_tile1 = task_tile;

    s_vtn  = locdbl;
    s_trcr = locdbl + usize;
    s_dzt  = locdbl + 3*usize;	
    s_tarea_r = locdbl + 4*usize;		
    sztmp  = locdbl + 5*usize;
    if(tadvtype != tadvt_cent)
      s_trcr_n  = locdbl + 2*usize;	
	
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
      get_reply2= 0;	  
	  
      athread_get(PE_MODE,vtn3d+task_pos1,s_vtn,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
      athread_get(PE_MODE,dzt+task_pos1,s_dzt,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);
	  
      athread_get(PE_MODE,tarea_r+task_pos1,s_tarea_r,SIZEDBL*task_size1,(void*)&get_reply,0,0,0);	  
	  
	  
#ifdef TAVG_R8 
      pgsize = SIZEDBL*task_size1;
#else
      pgsize = SIZEFLT*task_size1;
#endif

      athread_get(PE_MODE,m_bufloc + task_pos1,s_bufloc,pgsize,(void*)&get_reply2,0,0,0);
	  
      while(get_reply!=3);	  


      if(tadvtype == tadvt_cent)
        {
        if((task_pos1+task_size1)<nxy)
          athread_get(PE_MODE,trcr+(n*km+myk-1)*nxy+task_pos1,s_trcr,SIZEDBL*(task_size1+nx),(void*)&get_reply2,0,0,0);
        else
          {
          athread_get(PE_MODE,trcr+(n*km+myk-1)*nxy+task_pos1,s_trcr,SIZEDBL*task_size1,(void*)&get_reply2,0,0,0);
          for(i=0;i<nx;i++)
            s_trcr[task_size1+i] = c0;   // the last column is set to 0.
          }
        sync_mem; 
        while(get_reply2!=2);
        }
      else
        {
        athread_get(PE_MODE,trcr+(n*km+myk-1)*nxy+task_pos1,s_trcr,SIZEDBL*task_size1,(void*)&get_reply2,0,0,0);
        athread_get(PE_MODE,trcr_n+n*nxy+task_pos1,s_trcr_n,SIZEDBL*task_size1,(void*)&get_reply2,0,0,0);
        while(get_reply2!=3);
        sync_mem;
        }

 

/*
      if (partial_bottom_cells) then
         FVN =  p5*VTN*TAREA_R(:,:,bid)/DZT(:,:,k,bid)
      else
         FVN =  p5*VTN*TAREA_R(:,:,bid)
      endif

      do n=1,nt
         if (tadvect_itype(n) == tadvect_centered) then
            WORK = FVN*(TRCR(:,:,k,n) +  eoshift(TRCR(:,:,k,n),dim=2,shift=1))  // by column.
         else
            WORK = c2*FVN*TRACER_N(:,:,n)
         endif
         call accumulate_tavg_field(WORK,tavg_VN_TRACER(n),bid,k)
      enddo

where
     TRCR(nx,ny,km,nt)
     TRACER_N(nx,ny,nt) 
*/
 
      if(part_bott_cells>0)
         { 
         for(i=0;i<task_size1;i++)
           sztmp[i] = p5*s_vtn[i]*s_tarea_r[i]/s_dzt[i];
	 }
      else
	 { 
         for(i=0;i<task_size1;i++)
	   sztmp[i] = p5*s_vtn[i]*s_tarea_r[i];
	 }
		 
      if(tadvtype == tadvt_cent)
         { 
         for(i=0;i<task_size1;i++)
           sztmp[i] = sztmp[i]*(s_trcr[i]+s_trcr[i+nx]);   // the right column.
	 }
      else
	 { 
         for(i=0;i<task_size1;i++)
	   sztmp[i] = c2*sztmp[i]*s_trcr_n[i];
	 }
 

      accum_tavg_field_s(sztmp,s_bufloc,task_size1,tavgmeth,dtavg);


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


 return;
 }
 // end of s_advt_tavg_18.




