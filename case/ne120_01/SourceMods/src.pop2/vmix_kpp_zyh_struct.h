#include <stdlib.h>

#define SIZEINT 4

#if defined(POPCHKJN) 
#define MAXUSZ 256      //2KB
#define MAXUSZ1 128      //1KB
//#define MAXUSZ 64        //0.5KB
#else
#define MAXUSZ 512      //6KB
#define MAXUSZ1 256      //6KB
//#define MAXUSZ1 256      //2KB
#endif


typedef struct param_kpp_state_s2{
   int nblock[8];
   int core_info[3];
   int *local_iarr[5];
   double *local_arr[8];
   double *local_pts[8];
   int *liin[64];
   double *ldin[64];
   double *lval[64];
   double *sdiag[8];
}param_kpp_state_s2;

//#include <stdlib.h>

//#define SIZEINT 4

//#if defined(POPCHKJN) 
#if defined(POPTESTJN) 
#define MAXUSZ2 256      //2KB
#define MAXUSZBLD 260      //2.0KB
//#define MAXUSZBLD 160      //1.6KB
//#define MAXUSZBLD 80      //0.6KB
//#define MAXUSZ1 128      //1KB
//#define MAXUSZ 64        //0.5KB
#else
#define MAXUSZ2 760      //6KB
//#define MAXUSZBLD 80      //0.6KB
#define MAXUSZBLD 480      //3.8KB
//#define MAXUSZ1 760      //6KB
//#define MAXUSZ1 256      //2KB
#endif


typedef struct param_kpp_bldepth_s1{
   int nblock[8];
   int core_info[4];
   int *local_iarr[2];
   int global_logic[4];
   double *local_arr[6];
   double *global_arr[4];
   double local_pts[10];
/*
#if defined(POPTESTJN)||defined(POPCHKJN) 
   int *liin[64];
   double *ldin[64];
   double *lval[64];
   double *sdiag[8];
#endif
*/
}param_kpp_bldepth_s1;


