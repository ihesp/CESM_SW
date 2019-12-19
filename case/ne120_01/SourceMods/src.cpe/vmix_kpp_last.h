#include <stdlib.h>

#define SIZEINT sizeof(int)
#define SIZEDBL sizeof(double)

#if defined(POPCHKJN) 
#define MAXUSZ 256      //2KB
#else
#define MAXUSZ 400      //4KB
#endif


typedef struct param_kpp_last_s1{
   int logic[8];
   int ipts[8];
   int *iarr[8];
   double dpts[8];
   double *darr[16];
/*
#if defined(POPTESTJN)||defined(POPCHKJN) 
   int *liin[64];
   double *ldin[64];
   double *lval[64];
   double *sdiag[8];
#endif
*/
}param_kpp_last_s1;


