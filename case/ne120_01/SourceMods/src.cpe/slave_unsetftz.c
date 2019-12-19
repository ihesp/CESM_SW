#include <slave.h>

static inline void unsetftz(){
  long fpcr;
  asm volatile("rfpcr %0": "=r"(fpcr));
  fpcr |= 1L << 2;
  fpcr &= ~(1L << 60);
  asm volatile("wfpcr %0": : "r"(fpcr));
}

void slave_unsetftz_() {
  unsetftz();
}
