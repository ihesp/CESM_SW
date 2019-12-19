#ifndef SW5_PCR_H_
#define SW5_PCR_H_
#include <stdlib.h>
#include <string.h>

#define PCR_DEF_BEGIN(name) enum name ## _OPTS {
#define PCR_OPT_DEF(opt, value) opt = value,
#define PCR_DEF_END(name) name ## _MAX};
#include "pcrdef.h"
#undef PCR_DEF_BEGIN
#undef PCR_OPT_DEF
#undef PCR_DEF_END

#define PCR_DEF_BEGIN(name) char *name ##_NAMES[] = {
#define PCR_OPT_DEF(opt, value) #opt,
#define PCR_DEF_END(name) #name "_END"};
#include "pcrdef.h"
#undef PCR_DEF_BEGIN
#undef PCR_OPT_DEF
#undef PCR_DEF_END

#define PCRC_ALL_FLOP          (0x3f << 16)
#define PCRC_ALL_PC            (0x1f)
#define PCRC_ALL               (PCRC_ALL_FLOP | PCRC_ALL_PC)
#define PCRC_VDIV_VSQRT        (1 << 21) //浮点向量除法平方根指令计数使能。本位为1，从核性能计数器1的计数使能为1，且计数事件为5’h15（浮点除法、平方根类指令等效操作计数）时，浮点向量除法平方指令参与计数（每次+4）。
#define PCRC_FDIV_FSQRT        (1 << 20) //浮点标量除法平方根指令计数使能。本位为1，从核性能计数器1的计数使能为1，且计数事件为5’h15（浮点除法、平方根类指令等效操作计数）时，浮点标量除法平方指令参与计数（每次+1）。
#define PCRC_VMA               (1 << 19) //浮点向量乘加类指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点向量乘加类指令参与计数（每次+8）。
#define PCRC_VADD_VSUB_VMUL    (1 << 18) //浮点向量加减乘指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点向量加减乘指令参与计数（每次+4）。
#define PCRC_FMA               (1 << 17) //浮点标量乘加类指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点标量乘加类指令参与计数（每次+2）。
#define PCRC_FADD_FSUB_FMUL    (1 << 16) //浮点标量加减乘指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点标量加减乘指令参与计数（每次+1）。
#define PCRC_PC2_OVERFLOW      (1 <<  5) //TA性能计数器的溢出中断使能。为1时使能。   
#define PCRC_PC2               (1 <<  4) //TA性能计数器的计数使能。为1时使能。        
#define PCRC_PC1_OVERFLOW      (1 <<  3) //从核性能计数器1的溢出中断使能。为1时使能。
#define PCRC_PC1               (1 <<  2) //从核性能计数器1的计数使能。为1时使能。    
#define PCRC_PC0_OVERFLOW      (1 <<  1) //从核性能计数器0的溢出中断使能。为1时使能。
#define PCRC_PC0               (1 <<  0) //从核性能计数器0的计数使能。为1时使能。

#define HAVE_ALL(super, sub) (((super) & (sub)) == (sub))
static inline void decode_pcrc(char *ret, long pcrc_val){
  strcpy(ret, "(");
  if (HAVE_ALL(pcrc_val, PCRC_ALL)) {
    strcat(ret, "ALL, ");
  } else {
    if (HAVE_ALL(pcrc_val, PCRC_ALL_FLOP)) {
      strcat(ret, "ALL_FLOP, ");
    } else {
      if (HAVE_ALL(pcrc_val, PCRC_VDIV_VSQRT)) strcat(ret, "VDIV_VSQRT, ");
      if (HAVE_ALL(pcrc_val, PCRC_FDIV_FSQRT)) strcat(ret, "FDIV_FSQRT, ");
      if (HAVE_ALL(pcrc_val, PCRC_VMA)) strcat(ret, "VMA, ");
      if (HAVE_ALL(pcrc_val, PCRC_VADD_VSUB_VMUL)) strcat(ret, "VADD_VSUB_VMUL, ");
      if (HAVE_ALL(pcrc_val, PCRC_FMA)) strcat(ret, "FMA, ");
      if (HAVE_ALL(pcrc_val, PCRC_FADD_FSUB_FMUL)) strcat(ret, "FADD_FSUB_FMUL, ");
    }
    if (HAVE_ALL(pcrc_val, PCRC_ALL_PC)) {
      strcat(ret, "ALL_PC, ");
    } else {
      if (HAVE_ALL(pcrc_val, PCRC_PC0)) strcat(ret, "PC0, ");
      if (HAVE_ALL(pcrc_val, PCRC_PC1)) strcat(ret, "PC1, ");
      if (HAVE_ALL(pcrc_val, PCRC_PC2)) strcat(ret, "PC2, ");
    }
  }
  int retlen = strlen(ret);
  if (ret[retlen - 1] == ' ' && ret[retlen - 2] == ',')
    ret[retlen - 2] = 0;
  strcat(ret, ")");
}

#endif
