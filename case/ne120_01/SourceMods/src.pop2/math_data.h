#include <slave.h>
//数据表大小
#define sincos_data_len 570
#define acos_data_len 554
#define asin_data_len 554
#define exp_data_len 528
#define log_data_len 537
#define log10_data_len 540
#define pow_data_len 1054

#define sincos_data_bytes 4560
#define acos_data_bytes 4432
#define asin_data_bytes 4432
#define exp_data_bytes 4224
#define log_data_bytes 4320
#define log10_data_bytes 4320
#define pow_data_bytes 8432 
typedef union
{
        unsigned long LTYPE;
        double  DTYPE;
} LDTYPE;
extern LDTYPE sincos_data[sincos_data_len] ;
extern LDTYPE exp_data[exp_data_len] ;
extern LDTYPE acos_data[acos_data_len] ;
extern LDTYPE asin_data[asin_data_len] ;
extern LDTYPE log_data[log_data_len] ;
extern LDTYPE log10_data[log10_data_len] ;
extern LDTYPE pow_data[pow_data_len] ;

__thread_local extern double *sin_data_local_ptr ;
__thread_local extern double *cos_data_local_ptr ;
__thread_local extern double *acos_data_local_ptr ;
__thread_local extern double *asin_data_local_ptr ;
__thread_local extern double *exp_data_local_ptr ;
__thread_local extern double *log_data_local_ptr ;
__thread_local extern double *log10_data_local_ptr ;
__thread_local extern double *pow_data_local_ptr ;


