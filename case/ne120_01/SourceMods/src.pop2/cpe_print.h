#include <stdarg.h>
static void cpe_printf(const char *fmt ,...){
  volatile long vprintf_addr = (long)vprintf;

  int (*vprintf_ptr)(const char *, va_list) = (void*)vprintf_addr;
  va_list vlist;
  va_start(vlist, fmt);
  vprintf_ptr(fmt, vlist);
  va_end(vlist);
}

