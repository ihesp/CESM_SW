#include <slave.h>
#include <stdio.h>
#include <stdarg.h>

#include "cpe_print.h"

void test_(){
  if (_MYID == 0)
    cpe_printf("dxh cpe test!!!\n");
}
void printf_double_c(double in) {
    cpe_printf("cpe: %lf\n", in);
}
void printf_out_(int *in1, int *in2, int *in3){
    if(_MYID == 0)
        cpe_printf("out : %d %d %d\n", *in1, *in2, *in3);

}
void printf_c(int in1, int in2, int in3){
    if(_MYID == 0)
        cpe_printf("cout : %d %d %d\n", in1, in2, in3);

}
void printf_double_f_(double in1){
    if(_MYID == 0)
        cpe_printf("cout : %lf \n", in1);

}

void printf_double2_c(double in1, double in2){
    if(_MYID == 0)
        cpe_printf("cout : %lf %lf \n", in1, in2);

}
void printf_int_(int *in1){
    if(_MYID < 2)
        cpe_printf("cout : %d \n", *in1);

}
