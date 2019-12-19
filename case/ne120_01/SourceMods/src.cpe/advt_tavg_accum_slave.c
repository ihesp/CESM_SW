#include <unistd.h>  

//#if defined(POPCHKJN) || defined(POPTESTJN)
#ifdef POPINSJN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cpe_print.h"
#endif



 void accum_tavg_field_s(double *ARRAY, 
#ifdef TAVG_R8
     double *TAVGBUF
#else
     float  *TAVGBUF
#endif
     , int task_size, int tavgmeth, double dtavg)
/*
! !DESCRIPTION:
!  This routine updates a tavg field.  If the time average of the
!  field is requested, it accumulates a time sum of a field by 
!  multiplying by the time step and accumulating the sum into the 
!  tavg buffer array.  If the min or max of a field is requested, it
!  checks the current value and replaces the min, max if the current
!  value is less than or greater than the stored value.

   integer (int_kind), intent(in) :: &
      task_size       & ! size of double array ARRAY
      tavgmeth          ! method of tavg accumulation

   real (r8), dimension(task_size) :: &
      ARRAY           & ! array of data to add to accumulated sum in tavg buffer
      TAVGBUF           ! tavg buffer

   real (r8), intent(in) ::  &
      dtavg             ! time interval used in avg method
*/

 {
   int i,j;

//!  update the field into the tavg buffer

   switch(tavgmeth)
      {
      case 1:   // tavg_method_avg: accumulate running time sum for time avg
         for(i=0;i<task_size;i++)
           {
#ifdef TAVG_R8
           TAVGBUF[i] = TAVGBUF[i] + dtavg*ARRAY[i];
#else
           TAVGBUF[i] = TAVGBUF[i] + dtavg*ARRAY[i];          // ok!
//           TAVGBUF[i] = (float)((double)TAVGBUF[i] + (dtavg*ARRAY[i]));  // ok!
#endif
           }
         break;
         
      case 2:   // tavg_method_min: replace with current minimum value
         for(i=0;i<task_size;i++)
           {
#ifdef TAVG_R8
           if(ARRAY[i] < TAVGBUF[i]) 
              TAVGBUF[i] = ARRAY[i];
#else
           if((float)ARRAY[i] < TAVGBUF[i]) 
//              TAVGBUF[i] = (float)ARRAY[i];
              TAVGBUF[i] = ARRAY[i];
#endif
           }
         break;

      case 3:   // tavg_method_max: replace with current minimum value
         for(i=0;i<task_size;i++)
           {
#ifdef TAVG_R8
           if(ARRAY[i] > TAVGBUF[i])
              TAVGBUF[i] = ARRAY[i];
#else
           if((float)ARRAY[i] > TAVGBUF[i])
//              TAVGBUF[i] = (float)ARRAY[i];
              TAVGBUF[i] = ARRAY[i];
#endif
           }
         break;

      case 5:   // tavg_method_constant: overwrite with current value; intended for time-invariant fields
         for(i=0;i<task_size;i++)
           {
#ifdef TAVG_R8
           TAVGBUF[i] = ARRAY[i];
#else
//           TAVGBUF[i] = (float)ARRAY[i];
           TAVGBUF[i] = ARRAY[i];
#endif
           }
         break;

//      case 4: // tavg_method_qflux:  
//         TAVG_BUF_2D(:,:,block,bufloc) =  &
//         TAVG_BUF_2D(:,:,block,bufloc) + const*max (c0,ARRAY);
//         break;

      default:
         break;
      }

 return;
 }
 // end of accum_tavg_field_s.

