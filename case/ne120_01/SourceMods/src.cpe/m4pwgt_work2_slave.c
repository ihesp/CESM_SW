#include <unistd.h>     // MUST include this head file!

#if defined(POPCHKJN) || defined(POPTESTJN)
//#ifdef POPINSJN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cpe_print.h"
#endif


#define REG_PUTR(var,dest) \
asm volatile ("putr %0,%1\n"::"r"(var),"r"(dest):"memory")
#define REG_GETR(var) \
asm volatile ("getr %0\n":"=r"(var)::"memory")
#define REG_PUTC(var,dest) \
asm volatile ("putc %0,%1\n"::"r"(var),"r"(dest):"memory")
#define REG_GETC(var)  \
asm volatile ("getc %0\n":"=r"(var)::"memory")

#define sync_mem   \
asm volatile("memb")


void m4pwgt_work2_s(double *swork, double *svisc, double *au0, double *aun, double *aue, double *aune, double *u_buf, double *d_buf, int nx, int ny, int numj, int nthrds, int lrid, int lcid, int lid)
{
// BOP:

  int myid,myrid,mycid;
  int i,j,lny, pij;
  int nthreads, task_inum, task_jnum; 
 
  double tmpd, tmpw1, tmpw2;

  task_inum = nx;
  task_jnum = numj;
  lny = ny;
  nthreads = nthrds;

  myrid = lrid;
  mycid = lcid;  
  myid  = lid;

#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("step into swork2_s %d:\n", myid);
#endif
  

/*
 for(j=0;j<ny-1;j++)
    for(i=0;i<nx-1;i++)
       {
         WORK2(i,j) = AU0 (i,j)*VISC(i  ,j  ) + &
                      AUN (i,j)*VISC(i  ,j+1) + &
                      AUE (i,j)*VISC(i+1,j  ) + &
                      AUNE(i,j)*VISC(i+1,j+1)
       }
*/


// put first line to upside:
  if(myid > 0)
     {	   
     if(myid%8==0)
        {
        for(i=0;i<task_inum;i++)
           {
           u_buf[i] = svisc[(task_jnum-1)*task_inum+i];  
           REG_PUTC(u_buf[i],myrid-1);
           sync_mem;
           }
        }
     else 
        {
        if(mycid>0 && myrid%2==0 )
           { 
           for(i=0;i<task_inum;i++)
              {
              u_buf[i] = svisc[(task_jnum-1)*task_inum+i];  
	      REG_PUTR(u_buf[i],mycid-1);
              sync_mem;
	      }      
           }
        else if(mycid<7 && myrid%2==1)
           { 
           for(i=0;i<task_inum;i++)
              {
              u_buf[i] = svisc[(task_jnum-1)*task_inum+i];  
	      REG_PUTR(u_buf[i],mycid+1);
              sync_mem;
	      }      
           }
        }
     }    
// end of if(myid > 0) 
	  
#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("put visc @%d to up ok!\n", myid);
#endif


  
// get down line from downside:
  if(myid < (lny-1))
     {

     if(myid%8==7)
        {
        REG_GETC(d_buf[0]);
        sync_mem;
        for(i=0;i<(task_inum-1);i++)	
           {   
    	   REG_GETC(d_buf[i+1]);   
           sync_mem;

/*
 for(j=0;j<ny-1;j++)
    for(i=0;i<nx-1;i++)
       {
         WORK2(i,j) = AU0 (i,j)*VISC(i  ,j  ) + &
                      AUN (i,j)*VISC(i  ,j+1) + &
                      AUE (i,j)*VISC(i+1,j  ) + &
                      AUNE(i,j)*VISC(i+1,j+1)
       }
*/

// computing the down line:	  
           pij = (task_jnum-1)*task_inum+i; 
           tmpw1 = aun[pij] * d_buf[i];
           tmpw2 = aune[pij]* d_buf[i+1];
           tmpd  = tmpw1 + tmpw2;
           tmpw1 = au0[pij] * svisc[pij];
           tmpw2 = aue[pij] * svisc[pij+1];
           swork[pij] = tmpw1 + tmpw2 + tmpd;
           }
        }

     else
        {
        if((mycid<7 && myrid%2==0) || (mycid>0 && myrid%2==1))
           {
           REG_GETR(d_buf[0]);
           sync_mem;
           for(i=0;i<(task_inum-1);i++)	
              {   
    	      REG_GETR(d_buf[i+1]);   
              sync_mem;

#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("get visc @%d from down ok!\n", myid);
#endif

// computing the down line:	  
              pij = (task_jnum-1)*task_inum+i; 
              tmpw1 = aun[pij] * d_buf[i];
              tmpw2 = aune[pij]* d_buf[i+1];
              tmpd  = tmpw1 + tmpw2;
              tmpw1 = au0[pij] * svisc[pij];
              tmpw2 = aue[pij] * svisc[pij+1];
              swork[pij] = tmpw1 + tmpw2 + tmpd;
	      }
           }		
        }		

#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("computing down border block @%d ok!\n", myid);
#endif
     }



// computing inner lines:
  for(j=0;j<(task_jnum-1);j++)
     for(i=0;i<(task_inum-1);i++)	
        {
        pij = j*task_inum+i; 
        tmpw1 = au0[pij]*svisc[pij]          + aue[pij] *svisc[pij+1]; 
        tmpw2 = aun[pij]*svisc[pij+task_inum]+ aune[pij]*svisc[pij+task_inum+1];
        swork[pij] = tmpw1 + tmpw2;
        }   



#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("step out of swork2_s %d:\n", myid);
#endif
  
// EOP 
}
