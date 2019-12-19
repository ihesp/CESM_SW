#include <math.h>
#include <slave.h>
void  get_sp_(long long int *spa){
    long long int spt;
    asm volatile("mov $sp, %0 \n\t"
            : "=r"(spt)
            );
    *spa = spt; 
}

__thread_local long long int *log10_addr;
__thread_local long long int *log_addr;
__thread_local long long int *pow_addr;

__thread_local long long int *exp_addr;
__thread_local long long int *acos_addr;
__thread_local long long int *cos_addr;
//__thread_local long long int *sqrt_addr;

//__thread_local double *tmp1_addr;
//__thread_local double *tmp2_addr;

//__thread_local pver_c;
//__thread_local k_c;
//__thread_local wvflag_c;

//void init_tmp_addr_( int *pver, double *in1, double *in2){
//    pver_c = *pver;
//    tmp1_addr = in1;
//    tmp2_addr = in2;
//}
//void set_wvflag0_(){
//    wvflag_c = 0;
//}
//
//void set_wvflag1_(){
//    wvflag_c = 1;
//}
//
//void set_wvflag2_(){
//    wvflag_c = 2;
//}
//void set_wvflag_(int *flag){
//    wvflag_c = *flag;
//}
//void put_k_(int *k){
//    k_c = *k;
//}
//
//void put_tmp1_(int *flag, double *in1){
//    if(wvflag_c == *flag)
//        tmp1_addr[k_c - 1] = *in1; 
//}
//
//void put_tmp2_(int *flag, double *in1){
//    if(wvflag_c == *flag)
//        tmp2_addr[k_c - 1] = *in1; 
//}
double log10_c_(double *in1){
    double (*func_addr)(double);
    func_addr = log10_addr;
    //= (void*)*log10_addr;
    double out = func_addr(*in1);
    return out;
}

double log_c_(double *in1){
    double (*func_addr)(double);
    func_addr = log_addr;
    double out = func_addr(*in1);
    return out;

}
double pow_c_(double *in1, double *in2){
    double (*func_addr)(double, double);
    func_addr = pow_addr;
    double out = func_addr(*in1, *in2);
    return out;
}
double exp_c_(double *in1){
    double (*func_addr)(double);
    func_addr = exp_addr;
    double out = func_addr(*in1);
    return out;

}
double acos_c_(double *in1){
    double (*func_addr)(double);
    func_addr = acos_addr;
    double out = func_addr(*in1);
    return out;
}
double cos_c_(double *in1){
    double (*func_addr)(double);
    func_addr = cos_addr;
    double out = func_addr(*in1);
    return out;
}
//double sqrt_c_(double *in1){
//    double (*func_addr)(double);
//    func_addr = sqrt_addr;
//    double out = func_addr(*in1);
//    return out;
//}
void  math_agent_log10_c_(long long int *func_addr_addr){
    *func_addr_addr = log10;
    log10_addr = &log10;
}
void  math_agent_abs_c_(long long int *func_addr_addr){
  *func_addr_addr = fabs;
}
void  math_agent_log_c_(long long int *func_addr_addr){
    *func_addr_addr = log;
    log_addr = &log;
}

void  math_agent_pow_c_(long long int *func_addr_addr){
    *func_addr_addr = pow;
    pow_addr = &pow;
}

void  math_agent_exp_c_(long long int *func_addr_addr){
    *func_addr_addr = exp;
    exp_addr = &exp;
}
void  math_agent_sin_c_(long long int *func_addr_addr){
  *func_addr_addr = sin;
}

void  math_agent_sqrt_c_(long long int *func_addr_addr){
  *func_addr_addr = sqrt;
}

void  math_agent_acos_c_(){
    acos_addr = &acos;
}

void  math_agent_cos_c_(){
    cos_addr = &cos;
}


void  math_agent_1i1o_(long long int *func_addr_addr, double *in, double *out){
  double (*func_addr)(double) = (void*)*func_addr_addr;
  *out = func_addr(*in);
}

void  math_agent_2i1o_(long long int *func_addr_addr, double *in1, double *in2, double *out){
  double (*func_addr)(double, double) = (void*)*func_addr_addr;
  *out = func_addr(*in1,  *in2);
}
