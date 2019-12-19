#ifndef VMIX_KPP_STRUCT_H 
#define VMIX_KPP_STRUCT_H
#include <stdlib.h>


typedef struct ri_Loop01_Para
{ 
   int param[10];
   double dparam[4];
   int* duse_arr[2];
   double* use_arr[16];

}ri_Loop01_Para;

typedef struct ri_Loop02_Para
{ 
   int param[8];
   double dparam[2];
   int* duse_arr[2];
   double* use_arr[5];

}ri_Loop02_Para;

typedef struct ri_Loop03_Para
{ 
   int param[10];
   double dparam[10];
   double* use_arr[10];

}ri_Loop03_Para;
#endif
