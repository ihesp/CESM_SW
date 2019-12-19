#ifndef C_STRUCT_H
#define C_STRUCT_H
#include <stdlib.h>

typedef struct param_clinic
{
        int    param[20];
        double dparam[12];
        double *addr[20];
        int    *addr_int[4];
}param_clinic;
typedef struct little_clinic
{
        int    param[8];
        double dparam[6];
        double *addr[8];
        int    *addr_int[4];
}little_clinic;
#endif
