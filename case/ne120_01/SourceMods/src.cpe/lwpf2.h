#ifndef __LWPF__
#define __LWPF__

#ifndef __LWPF_SHARED__
#define __LWPF_SHARED__
typedef struct perf_config_t{
  long pcr0, pcr1, pcr2, pcrc;
} perf_config_t;
#define LWPF_WEAK_DATA __attribute__((weak)) __attribute__((section(".data")))
#define LWPF_WEAK __attribute__((weak))

#ifdef CPE
#include <simd.h>
#include <slave.h>

extern int256 lwpf_local_counter[];
static inline void set_perf_mode(long pcr0, long pcr1, long pcr2, long pcrc){
  pcr0 <<= 59;
  pcr1 <<= 59;
  pcr2 <<= 59;

  asm volatile("wcsr %0, 5\n"
               "wcsr %1, 6\n"
               "wcsr %2, 7\n"
               "wcsr %3, 8\n": :
               "r"(pcr0), "r"(pcr1), "r"(pcr2), "r"(pcrc));
}

static void lwpf_sync_counters_m2c(void *lwpf_cpe_counter, int kernel_count) {
  volatile int reply = 0;
  athread_get(PE_MODE, lwpf_cpe_counter, lwpf_local_counter, sizeof(int256) * kernel_count, (void*)&reply, 0, 0, 0);
  while (reply != 1);
}

static void lwpf_sync_counters_c2m(void *lwpf_cpe_counter, int kernel_count) {
  volatile int reply = 0;
  athread_put(PE_MODE, lwpf_local_counter, lwpf_cpe_counter, sizeof(int256) * kernel_count, (void*)&reply, 0, 0);
  while (reply != 1);
}

#define lwpf_enter(x) lwpf_sync_counters_m2c(lwpf_global_counter_ ## x[_MYID], lwpf_kernel_count_ ## x);

#define lwpf_exit(x) lwpf_sync_counters_c2m(lwpf_global_counter_ ## x[_MYID], lwpf_kernel_count_ ## x);

#define lwpf_start(kernel)						\
  {									\
    int256 cntrs;                                                       \
    asm volatile("rcsr %0, 7\n\t"					\
		 "vshff %0, %0, 0x24, %0\n\t"				\
		 "rcsr %0, 6\n\t"					\
		 "vshff %0, %0, 0xc4, %0\n\t"				\
		 "rcsr %0, 5\n\t"					\
		 "vshff %0, %0, 0xe0, %0\n\t"				\
		 "rcsr %0, 4\n\t"                                       \
                 "vsubl %2, %0, %1\n\t"                                 \
                 : "=&r"(cntrs), "=r"(lwpf_local_counter[kernel])       \
                 : "r"(lwpf_local_counter[kernel]));                    \
									\
  }									\

#define lwpf_stop(kernel)						\
  {									\
    int256 cntrs;                                                       \
    asm volatile("rcsr %0, 7\n\t"					\
		 "vshff %0, %0, 0x24, %0\n\t"				\
		 "rcsr %0, 6\n\t"					\
		 "vshff %0, %0, 0xc4, %0\n\t"				\
		 "rcsr %0, 5\n\t"					\
		 "vshff %0, %0, 0xe0, %0\n\t"				\
		 "rcsr %0, 4\n\t"                                       \
                 "vaddl %2, %0, %1\n\t"                                 \
                 : "=&r"(cntrs), "=r"(lwpf_local_counter[kernel])       \
                 : "r"(lwpf_local_counter[kernel]));                    \
									\
  }									\

#endif
#endif

#ifdef CPE

#define K(x) x,
typedef enum {
  LWPF_KERNELS
#define U(x) LWPF_KERNELS_END_ ## x
  LWPF_UNIT
#undef U
#define U(x) lwpf_kernel_ ## x
} LWPF_UNIT;
#undef U
#undef K

#define K(x) #x,
#define U(x) lwpf_kernel_names_ ## x
LWPF_WEAK_DATA char *LWPF_UNIT[] = {LWPF_KERNELS "LWPF_KERNELS_END"};
#undef K
#undef U

#define U(x) lwpf_kernel_count_ ## x
LWPF_WEAK const long LWPF_UNIT = 
#undef U
#define U(x) LWPF_KERNELS_END_ ## x
  LWPF_UNIT;
#undef U

#define U(x) lwpf_global_counter_ ## x
long LWPF_UNIT[64][
#undef U
#define U(x) LWPF_KERNELS_END_ ## x
                   LWPF_UNIT
                   ][4];


__thread_local int256 lwpf_local_counter[LWPF_UNIT] LWPF_WEAK;
#undef U

#define U(x) lwpf_init_ ## x
LWPF_WEAK void LWPF_UNIT(perf_config_t *conf){
#undef U
#define U(x) LWPF_KERNELS_END_ ## x
  int256 lwpf_local_counter[LWPF_UNIT];
  set_perf_mode(conf->pcr0, conf->pcr1, conf->pcr2, conf->pcrc);
  int256 v0 = 0;
  int i;
  for (i = 0; i < LWPF_UNIT; i ++){
    lwpf_local_counter[i] = v0;
  }
#undef U
#define U(x) lwpf_exit(x);
  LWPF_UNIT
#undef U
}
#endif

#ifdef MPE
#include <athread.h>
#include <stdio.h>
#include "sw5_pcr.h"

#define U(x) extern char *lwpf_kernel_names_ ## x[];
LWPF_UNITS
#undef U

#define U(x) extern long lwpf_global_counter_ ## x[];
LWPF_UNITS
#undef U
#define U(x) extern long lwpf_kernel_count_ ## x;
LWPF_UNITS
#undef U

#define U(x) extern void slave_lwpf_init_ ## x(perf_config_t *);
LWPF_UNITS
#undef U

static inline void lwpf_init(perf_config_t *conf){
#define U(x) athread_spawn(lwpf_init_ ## x, conf); athread_join();
  LWPF_UNITS
#undef U
}
//GMK1
#ifndef LWPF_NOCOLOR
#define LWPF_KG   "\x1B[31mG\x1B[0m"
#define LWPF_KM   "\x1B[33mM\x1B[0m"
#define LWPF_KK   "\x1B[32mK\x1B[0m"
#define LWPF_K1   "\x1B[30m."
#define LWPF_KNRM "\x1B[0m"
#else
#define LWPF_KG   "G"
#define LWPF_KM   "M"
#define LWPF_KK   "K"
#define LWPF_K1   "."
#define LWPF_KNRM ""
#endif
#define LWPF_DG   1000000000
#define LWPF_OG    500000000
#define LWPF_DM   1000000
#define LWPF_OM    500000
#define LWPF_DK   1000
#define LWPF_OK    500
#define LWPF_D1   1
#define LWPF_O1   0
static inline void print_unit_num_w5(FILE *output, long value){
  char *unit = LWPF_K1;
  int div = 1, off = LWPF_O1;
  if (value > LWPF_DK * 10L) {div = LWPF_DK; unit = LWPF_KK; off = LWPF_OK;}
  if (value > LWPF_DM * 10L) {div = LWPF_DM; unit = LWPF_KM; off = LWPF_OM;}
  if (value > LWPF_DG * 10L) {div = LWPF_DG; unit = LWPF_KG; off = LWPF_OG;}
  fprintf(output, "%4d%s%s", (value + off) / div, unit, LWPF_KNRM);
}

static inline void print_unit_num_w7(FILE *output, long value){
  char *unit = LWPF_K1;
  int div = 1, off = LWPF_O1;
  if (value > LWPF_DK * 1000L) {div = LWPF_DK; unit = LWPF_KK; off = LWPF_OK;}
  if (value > LWPF_DM * 1000L) {div = LWPF_DM; unit = LWPF_KM; off = LWPF_OM;}
  if (value > LWPF_DG * 1000L) {div = LWPF_DG; unit = LWPF_KG; off = LWPF_OG;}
  fprintf(output, "%6d%s%s", (value + off) / div, unit, LWPF_KNRM);
}

static inline void print_unit_num_w16(FILE *output, long value){
  char *unit = LWPF_K1;
  long div = 1, off = LWPF_O1;
  if (value > LWPF_DK * 100000000000L) {div = LWPF_DK; unit = LWPF_KK; off = LWPF_OK;}
  if (value > LWPF_DM * 100000000000L) {div = LWPF_DM; unit = LWPF_KM; off = LWPF_OM;}
  fprintf(output, "%15d%s%s", (value + off) / div, unit, LWPF_KNRM);
}

void center_print(FILE *output, const char *s, int width) {
  int length = strlen(s);
  int i;
  for (i=0; i< (width-length) / 2; i++) {
    fputs(" ", output);
  }
  fputs(s, output);
  i += length;
  for (; i < width; i++) {
    fputs(" ", output);
  }
}

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

static inline void lwpf_report_summary_one(FILE *output, char *uname, long *lwpf_global_counter, char **lwpf_kernel_names, long lwpf_kernel_count, perf_config_t *conf) {
  fprintf(output, "LWPF kernel summary for unit %s:\n", uname);
  //fprintf(output, "%11s|%17s|%17s|%17s|%17s\n", "KERNEL", "TC", PCR0_NAMES[conf->pcr0], PCR1_NAMES[conf->pcr1], PCR2_NAMES[conf->pcr2]);

  fprintf(output, "|");
  center_print(output, "KERNEL", 12);
  fprintf(output, "|");
  center_print(output, "TC", 17);
  fprintf(output, "|");
  center_print(output, PCR0_NAMES[conf->pcr0] + 4, 17);
  fprintf(output, "|");
  center_print(output, PCR1_NAMES[conf->pcr1] + 4, 17);
  fprintf(output, "|");
  center_print(output, PCR2_NAMES[conf->pcr2] + 4, 17);
  fprintf(output, "|\n");


  char *avgminmax = " avg | min | max ";
  fprintf(output, "|%12s|%17s|%17s|%17s|%17s|\n", "", avgminmax, avgminmax, avgminmax, avgminmax);
  //"(avg/min/max)", "(avg/min/max)", "(avg/min/max)", "(avg/min/max)");
  char *dashes12 = "------------";
  char *plus5 = "-----+-----+-----";
  fprintf(output, "|%s+%s+%s+%s+%s|\n", dashes12, plus5, plus5, plus5, plus5);
  int i, j;
  for (i = 0; i < lwpf_kernel_count; i ++){
    //fprintf(output, "%20s|", lwpf_kernel_names[i]);
    fprintf(output, "|");
    center_print(output, lwpf_kernel_names[i], 12);
    fprintf(output, "|");
    long rpcc_sum = 0, rpcc_max = 0, rpcc_min = 1L << 59;
    long pcr0_sum = 0, pcr0_max = 0, pcr0_min = 1L << 59;
    long pcr1_sum = 0, pcr1_max = 0, pcr1_min = 1L << 59;
    long pcr2_sum = 0, pcr2_max = 0, pcr2_min = 1L << 59;
    for (j = 0; j < 64; j ++){
      rpcc_sum += lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 0];
      pcr0_sum += lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 1];
      pcr1_sum += lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 2];
      pcr2_sum += lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 3];
      rpcc_max = max(rpcc_max, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 0]);
      pcr0_max = max(pcr0_max, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 1]);
      pcr1_max = max(pcr1_max, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 2]);
      pcr2_max = max(pcr2_max, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 3]);
      rpcc_min = min(rpcc_min, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 0]);
      pcr0_min = min(pcr0_min, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 1]);
      pcr1_min = min(pcr1_min, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 2]);
      pcr2_min = min(pcr2_min, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 3]);
    }

    long rpcc_avg = rpcc_sum / 64;
    long pcr0_avg = pcr0_sum / 64;
    long pcr1_avg = pcr1_sum / 64;
    long pcr2_avg = pcr2_sum / 64;
    print_unit_num_w5(output, rpcc_avg);
    fprintf(output, "|");
    print_unit_num_w5(output, rpcc_min);
    fprintf(output, "|");
    print_unit_num_w5(output, rpcc_max);
    fprintf(output, "|");
    
    print_unit_num_w5(output, pcr0_avg);
    fprintf(output, "|");
    print_unit_num_w5(output, pcr0_min);
    fprintf(output, "|");
    print_unit_num_w5(output, pcr0_max);
    fprintf(output, "|");

    print_unit_num_w5(output, pcr1_avg);
    fprintf(output, "|");
    print_unit_num_w5(output, pcr1_min);
    fprintf(output, "|");
    print_unit_num_w5(output, pcr1_max);
    fprintf(output, "|");

    print_unit_num_w5(output, pcr2_avg);
    fprintf(output, "|");
    print_unit_num_w5(output, pcr2_min);
    fprintf(output, "|");
    print_unit_num_w5(output, pcr2_max);
    fprintf(output, "|\n");
    //fprintf(output, "%6d/%6d/%6d|", 
  }
}

static inline void lwpf_report_summary_one_wide(FILE *output, char *uname, long *lwpf_global_counter, char **lwpf_kernel_names, long lwpf_kernel_count, perf_config_t *conf) {
  fprintf(output, "LWPF kernel summary for unit %s:\n", uname);

  fprintf(output, "|");
  center_print(output, "KERNEL", 12);
  fprintf(output, "|");
  center_print(output, "TC", 23);
  fprintf(output, "|");
  center_print(output, PCR0_NAMES[conf->pcr0] + 4, 23);
  fprintf(output, "|");
  center_print(output, PCR1_NAMES[conf->pcr1] + 4, 23);
  fprintf(output, "|");
  center_print(output, PCR2_NAMES[conf->pcr2] + 4, 23);
  fprintf(output, "|\n");


  char *avgminmax = "  avg  |  min  |  max  ";
  fprintf(output, "|%12s|%23s|%23s|%23s|%23s|\n", "", avgminmax, avgminmax, avgminmax, avgminmax);
  //"(avg/min/max)", "(avg/min/max)", "(avg/min/max)", "(avg/min/max)");
  char *dashes12 = "------------";
  char *plus5 = "-------+-------+-------";
  fprintf(output, "|%s+%s+%s+%s+%s|\n", dashes12, plus5, plus5, plus5, plus5);
  int i, j;
  for (i = 0; i < lwpf_kernel_count; i ++){
    //fprintf(output, "%20s|", lwpf_kernel_names[i]);
    fprintf(output, "|");
    center_print(output, lwpf_kernel_names[i], 12);
    fprintf(output, "|");
    long rpcc_sum = 0, rpcc_max = 0, rpcc_min = 1L << 59;
    long pcr0_sum = 0, pcr0_max = 0, pcr0_min = 1L << 59;
    long pcr1_sum = 0, pcr1_max = 0, pcr1_min = 1L << 59;
    long pcr2_sum = 0, pcr2_max = 0, pcr2_min = 1L << 59;
    for (j = 0; j < 64; j ++){
      rpcc_sum += lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 0];
      pcr0_sum += lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 1];
      pcr1_sum += lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 2];
      pcr2_sum += lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 3];
      rpcc_max = max(rpcc_max, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 0]);
      pcr0_max = max(pcr0_max, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 1]);
      pcr1_max = max(pcr1_max, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 2]);
      pcr2_max = max(pcr2_max, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 3]);
      rpcc_min = min(rpcc_min, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 0]);
      pcr0_min = min(pcr0_min, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 1]);
      pcr1_min = min(pcr1_min, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 2]);
      pcr2_min = min(pcr2_min, lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 3]);
    }

    long rpcc_avg = rpcc_sum / 64;
    long pcr0_avg = pcr0_sum / 64;
    long pcr1_avg = pcr1_sum / 64;
    long pcr2_avg = pcr2_sum / 64;
    print_unit_num_w7(output, rpcc_avg);
    fprintf(output, "|");
    print_unit_num_w7(output, rpcc_min);
    fprintf(output, "|");
    print_unit_num_w7(output, rpcc_max);
    fprintf(output, "|");
    
    print_unit_num_w7(output, pcr0_avg);
    fprintf(output, "|");
    print_unit_num_w7(output, pcr0_min);
    fprintf(output, "|");
    print_unit_num_w7(output, pcr0_max);
    fprintf(output, "|");

    print_unit_num_w7(output, pcr1_avg);
    fprintf(output, "|");
    print_unit_num_w7(output, pcr1_min);
    fprintf(output, "|");
    print_unit_num_w7(output, pcr1_max);
    fprintf(output, "|");

    print_unit_num_w7(output, pcr2_avg);
    fprintf(output, "|");
    print_unit_num_w7(output, pcr2_min);
    fprintf(output, "|");
    print_unit_num_w7(output, pcr2_max);
    fprintf(output, "|\n");
    //fprintf(output, "%6d/%6d/%6d|", 
  }
}

static inline void lwpf_report_detail_one(FILE *output, char *uname, long *lwpf_global_counter, char **lwpf_kernel_names, long lwpf_kernel_count, perf_config_t *conf) {
  fprintf(output, "LWPF kernel details for unit %s:\n", uname);

  fprintf(output, "|");
  center_print(output, "KERNEL", 12);
  fprintf(output, "|");
  center_print(output, "PE", 4);
  fprintf(output, "|");
  center_print(output, "TC", 16);
  fprintf(output, "|");
  center_print(output, PCR0_NAMES[conf->pcr0] + 4, 16);
  fprintf(output, "|");
  center_print(output, PCR1_NAMES[conf->pcr1] + 4, 16);
  fprintf(output, "|");
  center_print(output, PCR2_NAMES[conf->pcr2] + 4, 16);
  fprintf(output, "|\n");

  fprintf(output, "|------------+----+----------------+----------------+----------------+----------------|\n");
  int i, j;
  for (i = 0; i < lwpf_kernel_count; i ++){
    for (j = 0; j < 64; j ++){
      fprintf(output, "|");
      if (j == 0)
        center_print(output, lwpf_kernel_names[i], 12);
      else
        center_print(output, "", 12);
      fprintf(output, "|");
      fprintf(output, "%4d|", j);
      long rpcc = lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 0];
      long pcr0 = lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 1];
      long pcr1 = lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 2];
      long pcr2 = lwpf_global_counter[(j * lwpf_kernel_count + i) * 4 + 3];
      print_unit_num_w16(output, rpcc);
      fprintf(output, "|");
      print_unit_num_w16(output, pcr0);
      fprintf(output, "|");
      print_unit_num_w16(output, pcr1);
      fprintf(output, "|");
      print_unit_num_w16(output, pcr2);
      fprintf(output, "|");
      fprintf(output, "\n");
    }
    if (i < lwpf_kernel_count - 1)
      fprintf(output, "|------------+----+----------------+----------------+----------------+----------------|\n");
  }
}
#undef max
#undef min

static inline void lwpf_report_summary(FILE *output, perf_config_t *conf){
#define U(x) lwpf_report_summary_one(output, #x, lwpf_global_counter_ ## x, lwpf_kernel_names_ ## x, lwpf_kernel_count_ ## x, conf);
  LWPF_UNITS
#undef U
}

static inline void lwpf_report_summary_wide(FILE *output, perf_config_t *conf){
#define U(x) lwpf_report_summary_one_wide(output, #x, lwpf_global_counter_ ## x, lwpf_kernel_names_ ## x, lwpf_kernel_count_ ## x, conf);
  LWPF_UNITS
#undef U
}

static inline void lwpf_report_detail(FILE *output, perf_config_t *conf){
#define U(x) lwpf_report_detail_one(output, #x, lwpf_global_counter_ ## x, lwpf_kernel_names_ ## x, lwpf_kernel_count_ ## x, conf);
  LWPF_UNITS
#undef U
}

#endif
#endif
