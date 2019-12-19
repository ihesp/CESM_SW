#include <unistd.h>     // MUST include this head file!

//#if defined(POPCHKJN) || defined(POPTESTJN)
#ifdef POPINSJN
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

#define locdmax(a,b) ( (a)>(b)?(a):(b) )


void m4pmax_vshear_s(double *sshear, double *swork, double *u_buf, double *d_buf, int nx, int ny, int numj, int nthrds, int lrid, int lcid, int lid)
{
// BOP:

  int myid,myrid,mycid;
  int i,j,lny;
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
    cpe_printf("step into vshear_s %d:\n", myid);
#endif
  
  


// put last line to downside:
  if(myid < (lny-1))
     {	   
     if(myid%8==7 && myrid<7)
        {
        for(i=0;i<task_inum;i++)
           {
           d_buf[i] = swork[(task_jnum-1)*task_inum+i];  
           REG_PUTC(d_buf[i],myrid+1);
           sync_mem;
           }
        }
     else 
        {
        if(mycid<7 && myrid%2==0 )
           { 
           for(i=0;i<task_inum;i++)
              {
              d_buf[i] = swork[(task_jnum-1)*task_inum+i];  
	      REG_PUTR(d_buf[i],mycid+1);
              sync_mem;
	      }      
           }
        else if(mycid>0 && myrid%2==1)
           { 
           for(i=0;i<task_inum;i++)
              {
              d_buf[i] = swork[(task_jnum-1)*task_inum+i];  
	      REG_PUTR(d_buf[i],mycid-1);
              sync_mem;
	      }      
           }
        }
     }    
// end of if(myid<ny) 
	  
#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("put work @%d to down ok!\n", myid);
#endif


  
// get up line from upside:
  if(myid>0)
     {

     if(myid%8==0)
        {
        REG_GETC(u_buf[0]);
        sync_mem;
        for(i=1;i<task_inum;i++)	
           {   
    	   REG_GETC(u_buf[i]);   
           sync_mem;

// computing the upper line:	  
           tmpw1 = locdmax(swork[0*task_inum+i], swork[0*task_inum+(i-1)]); 
           tmpw2 = locdmax(u_buf[i]            , u_buf[i-1]);
           sshear[0*task_inum+i] = locdmax(tmpw1, tmpw2);
           }
        }

     else
        {
        if((mycid>0 && myrid%2==0) || (mycid<7 && myrid%2==1))
           {
           REG_GETR(u_buf[0]);
           sync_mem;
           for(i=1;i<task_inum;i++)	
              {   
    	      REG_GETR(u_buf[i]);   
              sync_mem;

#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("get work @%d from upper ok!\n", myid);
#endif

// computing the upper line:	  
              tmpw1 = locdmax(swork[0*task_inum+i], swork[0*task_inum+(i-1)]); 
              tmpw2 = locdmax(u_buf[i]            , u_buf[i-1]);
              sshear[0*task_inum+i] = locdmax(tmpw1, tmpw2);
	      }
           }		
        }		

#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("computing upper border block @%d ok!\n", myid);
#endif
     }



// computing inner lines:
  for(j=1;j<task_jnum;j++)
     for(i=1;i<task_inum;i++)	
        {
        tmpw1 = locdmax(swork[j*task_inum+i]    ,swork[j*task_inum+(i-1)]); 
        tmpw2 = locdmax(swork[(j-1)*task_inum+i],swork[(j-1)*task_inum+(i-1)]);
        sshear[j*task_inum+i] = locdmax(tmpw1, tmpw2);
        }   



#ifdef POPINSJN
 if(myid>10 && myid<15)
    cpe_printf("step out of vshear_s %d:\n", myid);
#endif
  
// EOP 
}
