#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>     // MUST include this head file!
#include "slave.h"
//#include "cpe_print.h"
#include "vmix_kpp_gq_struct.h"
#define s_double sizeof(double)


#define REG_PUTR(var,dest) \
asm volatile ("putr %0,%1\n"::"r"(var),"r"(dest):"memory")
#define REG_GETR(var) \
asm volatile ("getr %0\n":"=r"(var)::"memory")
#define REG_PUTC(var,dest) \
asm volatile ("putc %0,%1\n"::"r"(var),"r"(dest):"memory")
#define REG_GETC(var)  \
asm volatile ("getc %0\n":"=r"(var)::"memory")

#define get_myrid(row)  \
asm volatile("rcsr   %0, 1" : "=r"(row))
#define get_mycid(col)  \
asm volatile("rcsr   %0, 2" : "=r"(col))



double *ldm_malloc(); 


void ugrid_to_tgrid_slave(double *TGRID,double UGRID,double UGRID1,double AT0_U,double ATS_U,
                  double ATW_U,double ATSW_U,int ny)
{
	

 int i,j,k,myid,myrid,mycid;

 double u_buf,d_buf;
 double c0=0.0;

 get_myrid(myrid);
 get_mycid(mycid);
   
 
   if(myrid%2>0)
      myid=myrid*8+(7-mycid);
   else
     myid=myrid*8+mycid;  

     
    if(myid<(ny-1))   
       {
		   
	   if(myrid%2==0&&mycid<7)
	    {					  
	       REG_PUTR(UGRID,mycid+1);
	       REG_PUTR(UGRID1,mycid+1);
               
           }
#ifdef mem    
		 if(mycid==1&&my_task==0)
	         cpe_printf("mycid=\n",mycid);
#endif
           if(myrid%2==1&&mycid>0) 
	     {		
	            REG_PUTR(UGRID,mycid-1);
	            REG_PUTR(UGRID1,mycid-1);
                    
	      }
#ifdef mem
              if(mycid==1&&my_task==0)
               cpe_printf("REG_PUTR\n");
#endif
           if(myid%8==7)
			   
	      {
	            REG_PUTC(UGRID,myrid+1);
	            REG_PUTC(UGRID1,myrid+1);
	          
      	    }
	   }	
           		
#ifdef mem    
            if(my_task==0&&mycid==1)
              cpe_printf("REG_PUTC\n");
#endif
                     
                                 
       if(myid>0&&myid<ny)
       {
        if(myid%8==0)
        {			
          
             
	     REG_GETC(u_buf);
             REG_GETC(d_buf);		  		   

             *TGRID= AT0_U*UGRID + 
                      ATS_U*u_buf + 
                      ATW_U*UGRID1 + 
                      ATSW_U*d_buf;
          
             
    	   
         }
        else   		 
        {   
        				   
	     REG_GETR(u_buf);
             REG_GETR(d_buf);		  		   
	     		   
            
       
          *TGRID = AT0_U*UGRID + 
                      ATS_U*u_buf + 
                      ATW_U*UGRID1 + 
                      ATSW_U*d_buf;
                 
    	    
          }
#ifdef mem    
           if(my_task==0&&mycid==1)
           cpe_printf("REG_GET\n");
#endif       
       }
      
    
      	  
      else
     	  	   
	 
          *TGRID = c0; 
         
	  
    
#ifdef mem    
      if(my_task==0&&mycid==1)
          cpe_printf("ldm_free\n");
#endif
 
}  

