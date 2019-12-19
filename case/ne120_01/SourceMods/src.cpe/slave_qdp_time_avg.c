#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define KBLK 12 // KBLB is unit block size which along k axis per cpe
#define NC 8
#define NR 8
#define qstep_Qdp (NLEV*NP*NP) // stripe in q axis of Qtens_biharmonic array
#define block  (KBLK*NP*NP)


#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0":"=r"(var))
#define MAX(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define MIN(val1, val2) ((val1) < (val2) ? (val1) : (val2))


typedef struct {
  double *gl_qdp;
  int rkstage, n0_qdp, np1_qdp, nets, nete, qsize        \
      , qsize_d, step_elem;
} param_t;

void slave_qdp_time_avg_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
  get_col_id(cid);
  get_row_id(rid);
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();
  double *gl_qdp = param_d.gl_qdp;
  int rkstage = param_d.rkstage;
  int n0_qdp = param_d.n0_qdp;
  int np1_qdp = param_d.np1_qdp;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;
  int qsize_d = param_d.qsize_d;
  int step_elem = param_d.step_elem;
  int sum_k = qsize*NLEV;

  double Qdp_n0[block];
  double Qdp_np1[block];

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block + NLEV*NP*NP*3 + block + NLEV*2       \
        + 2*UC*NLEV);
    printf("slave_euler_step size_tol:%dk\n", size_tol/1024);
  }
#endif

  double *gl_qdp_n0 = gl_qdp + (n0_qdp - 1)*qsize_d*NLEV*NP*NP;
  double *gl_qdp_np1 = gl_qdp + (np1_qdp - 1)*qsize_d*NLEV*NP*NP;
  double *src_qdp_n0, *src_qdp_np1;

  int i, j, k, ie, k_beg, k_end, k_n;

  for (ie = 0; ie < nete - nets + 1; ie++) {
    k_beg = id*KBLK;
    k_end = (id + 1)*KBLK;
    k_end = k_end < sum_k ? k_end : sum_k;
    k_n = k_end - k_beg;

    //if (id == 63) printf("k_beg:%d, k_n:%d\n", k_beg, k_n);
    if (k_n > 0) {
      src_qdp_n0 = gl_qdp_n0 + ie*step_elem + k_beg*NP*NP;
      src_qdp_np1 = gl_qdp_np1 + ie*step_elem + k_beg*NP*NP;
      pe_get(src_qdp_n0, Qdp_n0, k_n*NP*NP*sizeof(double));
      pe_get(src_qdp_np1, Qdp_np1, k_n*NP*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos = k*NP*NP + j*NP + i;
            Qdp_np1[pos] = (Qdp_n0[pos] + (rkstage - 1)*Qdp_np1[pos])/rkstage;
          }
        }
      }
      pe_put(src_qdp_np1, Qdp_np1, k_n*NP*NP*sizeof(double));
      dma_syn();
    }
  }
}
