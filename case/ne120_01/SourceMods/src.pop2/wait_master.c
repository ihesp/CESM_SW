#include<athread.h>

void wait_fun_(int *pst, int *ped)
        
{ 
   
   int i;
   for(i=*pst-1;i<*ped;i++)
      athread_wait(i);
	
 }

// may directly in .F90:
// call athread_wait()


