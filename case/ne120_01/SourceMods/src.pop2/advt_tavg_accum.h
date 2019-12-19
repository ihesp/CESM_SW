#include <stdlib.h>

#define SIZEINT sizeof(int)
#define SIZEFLT sizeof(float)
#define SIZEDBL sizeof(double)

#if defined(POPCHKJN) 
#define MAXUSZ 512      //4KB
#else
#define MAXUSZ 768      //6KB
#endif


typedef struct param_tavg_accum{
   int tavgindx[32*12];
   int ipts[16];
   int *iarr[8];
   double dpts[8];
   float *farr[4];
   double *darr[16];
}param_tavg_accum;


