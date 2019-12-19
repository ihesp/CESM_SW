#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define KBLK  12


#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0":"=r"(var))
#define max(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define min(val1, val2) ((val1) < (val2) ? (val1) : (val2))


typedef struct {
  double *Qtens_biharmonic, *dpdiss_ave, *dp0;
  int nets, nete, qsize, step_elem;
} param_t;

void slave_rhs_multiplier_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();
  double *gl_Qtens_biharmonic = param_d.Qtens_biharmonic;
  double *gl_dpdiss_ave = param_d.dpdiss_ave;
  double *gl_dp0 = param_d.dp0;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;
  int step_elem = param_d.step_elem;

  double Qtens_biharmonic[KBLK*NP*NP];
  double dpdiss_ave[KBLK*NP*NP];
  double dp0[NLEV];

  int sum_k = qsize*NLEV;

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block + NLEV*NP*NP*3 + block + NLEV*2       \
        + 2*UC*NLEV);
    printf("slave_euler_step size_tol:%dk\n", size_tol/1024);
  }
#endif

  double *src_Qtens_biharmonic, *src_dpdiss_ave, *src_dp0;

  int istep_Qtens = qsize*NLEV*NP*NP;
  int c, r, i, j, k, q, ie, k_beg, k_end, q_beg, q_end, k_n, q_n, pos_dp
      , pos_qdp, pos_Qtens_bi, pos_qmax;

  pe_get(gl_dp0, dp0, (NLEV*sizeof(double)));
  dma_syn();

  for (ie = 0; ie < nete - nets + 1; ie++) {
    // roll the q and k into the id of cpe,
    // other words is that (q, k) has relationship with id
    // be careful that the edge of possbilly exception
    q_beg = KBLK*id/NLEV;
    q_end = (KBLK*(id + 1))/NLEV;
    q_end = q_end < qsize ? q_end : qsize;
    k_beg = KBLK*id - q_beg*NLEV;
    k_end = KBLK*(id + 1) - q_beg*NLEV;
    k_end = k_end < NLEV ? k_end : NLEV;
    k_n = k_end - k_beg;
    q_n = (k_n == 6 && q_beg != (qsize - 1)) ? 1 : 0;
    q_n = (id*KBLK > sum_k) ? -1 : q_n;
    for (q = 0; q <= q_n; q++) {
      src_dpdiss_ave = gl_dpdiss_ave + ie*step_elem                            \
          + (k_beg + q*(q*NLEV - k_beg) - q*NLEV)*NP*NP;
      src_Qtens_biharmonic = gl_Qtens_biharmonic + ie*istep_Qtens              \
          + q_beg*NLEV*NP*NP + (k_beg + q*(q*NLEV - k_beg))*NP*NP;
      pe_get(src_dpdiss_ave, dpdiss_ave, (k_n*NP*NP*sizeof(double)));
      pe_get(src_Qtens_biharmonic, Qtens_biharmonic, (k_n*NP*NP*sizeof(double)));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos = k*NP*NP + j*NP + i;
            int pos_dpdiss = k*NP*NP + j*NP + i;
            pos_dp = (((k_beg + q*(q*NLEV - k_beg)) + k) >= 30)             \
                ? k : (k_beg + q*(q*NLEV - k_beg)) + k;
            Qtens_biharmonic[pos] = Qtens_biharmonic[pos]*dpdiss_ave[pos_dpdiss]/dp0[pos_dp];
          }
        }
      }
      pe_put(src_Qtens_biharmonic, Qtens_biharmonic, (k_n*NP*NP*sizeof(double)));
      dma_syn();
    }
  }
}
