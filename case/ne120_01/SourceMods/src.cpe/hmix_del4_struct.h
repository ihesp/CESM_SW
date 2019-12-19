#ifndef HDIFFU_STRUCT_H
#define HDIFFU_STRUCT_H
#include <stdlib.h>

typedef struct param_hdiffu
{ 
   int param[16];
   double dparam[12];
   double *addr[32];
   int *addr_int[2];
}param_hdiffu;

#endif
