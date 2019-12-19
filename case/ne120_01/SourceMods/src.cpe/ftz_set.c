/*************************************************************************
  > File Name: src.cam/ftz_set.c
  > Author: Xu Kai
  > Created Time: 2019年01月17日 星期四 13时42分53秒
 ************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<slave.h>
void setftz_f_(){
    long fpcr;
    asm volatile("rfpcr %0": "=r"(fpcr));
    fpcr |= 1L << 60 | 1L << 2;
    asm volatile("wfpcr %0": : "r"(fpcr));
}

void unsetftz_f_(){
    long fpcr;
    asm volatile("rfpcr %0": "=r"(fpcr));
    fpcr |= 1L << 2;
    fpcr &= ~(1L << 60);
    asm volatile("wfpcr %0": : "r"(fpcr));
}
