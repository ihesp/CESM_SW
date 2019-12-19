#include<stdio.h>
#include<slave.h>
#include"cpe_print.h"
#include <unistd.h>
#include <string.h>


#define REG_PUTR(var,dest) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dest):"memory")
#define REG_GETR(var)      asm volatile ("getr %0\n":"=r"(var)::"memory")
#define REG_PUTC(var,dest) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dest):"memory")
#define REG_GETC(var)      asm volatile ("getc %0\n":"=r"(var)::"memory")


void Reg_comm_zy(int nx,int ny,double *p_X,double* p_S,double* p_btropWgtCenter,double* p_btropWgtNorth,double* p_btropWgtEast,double*p_btropWgtNE,
             int my_rid,int my_cid,int myid,double *left_x,double *right_x,double *left_nor,double *left_ne)
{
           int i,j,last_slave,pre,pre2;
           

           if(myid >= ny ){}
           else {         
           if(my_rid & 1) {pre =my_cid-1;pre2=my_cid+1; } //1,3
           else           {pre =my_cid+1;pre2=my_cid-1; } //0,2  
 
           last_slave = ny;
           for(i=0;i<nx;i++)  //px -1
           {
               if(myid < (last_slave-1))
               {
                  if((myid % 8 == 7) && (my_rid < 7))                  
                  {

                      REG_PUTC(p_X[i],my_rid+1 );
                      asm volatile("memb");

                  }                           
                  else        
                  {

                      REG_PUTR( p_X[i],pre );                          
                      asm volatile("memb");

                  }                                         
               }      
              if(myid>0)
              {
                if(myid %8==0)               
                {

                     REG_GETC(left_x[i]); 
                     asm volatile("memb");

                }    
                else
                {

                     REG_GETR(left_x[i]);
                     asm volatile("memb");

                }                     
             }  
         }

      for(i=0;i<nx;i++)
       {
             if(myid >0)  //px+1
             {
               if(myid % 8 == 0 )
               {

                   REG_PUTC( p_X[i],my_rid-1 );
                   asm volatile("memb");
               }
               else
               {

                   REG_PUTR(p_X[i],pre2 );
                   asm volatile("memb");
               }
            }
            if(myid<last_slave-1)
            {
                 if(myid %8==7)
                 {

                     REG_GETC(right_x[i]);
                     asm volatile("memb");
                 }
                 else
                 {

                      REG_GETR(right_x[i]);
                      asm volatile("memb");
                 }
            }
       }
      if(myid >0   && myid < last_slave-1)   //!!!!!j +-1
      {
         for(i=1 ;i<nx-1; i++)
         {       
                     p_S[i] =   p_btropWgtCenter[i]     *    p_X[i]       +
                                p_btropWgtNorth[i+nx]   *    right_x[i]   +
                                p_btropWgtNorth[i]      *    left_x[i]    +
                                p_btropWgtEast[i]       *    p_X[i+1]     +
                                p_btropWgtEast[i-1]     *    p_X[i-1]     +
                                p_btropWgtNE[i+nx]      *    right_x[i+1] +
                                p_btropWgtNE[i]         *    left_x[i+1]  +
                                p_btropWgtNE[i-1+nx]    *    right_x[i-1] +
                                p_btropWgtNE[i-1]       *    left_x[i-1] ;
     }
    }          
   } 
}

