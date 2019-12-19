#include <string.h>

#define deftransfer(ftype, ctype) \
void transfer_ ## ftype ##_c_ (void *src, void *dst) \
{ \
  *(ctype*) dst = *(ctype*) src; \
}

deftransfer(r4, float)
deftransfer(r8, double)
