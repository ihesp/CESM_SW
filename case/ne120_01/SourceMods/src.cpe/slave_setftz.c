#include <slave.h>

static inline void setftz(){
  long fpcr;
  asm volatile("rfpcr %0": "=r"(fpcr));
  fpcr |= 1L << 60 | 1L << 2;
  asm volatile("wfpcr %0": : "r"(fpcr));
}

void slave_setftz_() {
  setftz();
}
