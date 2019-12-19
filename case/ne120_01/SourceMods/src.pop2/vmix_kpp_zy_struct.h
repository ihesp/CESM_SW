#include <stdlib.h>

typedef struct kpp_ddmix{
   int     param[16];
   double  dparam[16];
   double  *arr[16];
   int     *arr_int[16];
   double  *std[16];
}kpp_ddmix;
