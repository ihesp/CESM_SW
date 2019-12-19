#include <stdio.h>
void print_binary_(double *a, double *b) {
  printf("%016lx, %016lx \n", *((long *) a), *((long *) b));
}
