#ifndef C_STRUCT_H
#define C_STRUCT_H
#include <stdlib.h>

typedef struct param_vmix
{
        int    param[16];
        double dparam[8];
        double *addr[12];
        int    *addr_int[4];
}param_vmix;
#endif
