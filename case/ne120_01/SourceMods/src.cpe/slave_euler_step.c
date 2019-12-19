#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define KBLK 12
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
  double *qdp_s, *qdp_leap, *dp, *divdp_proj, *Qtens_biharmonic, *qmax, *qmin;
  double dt;
  int  nets, nete, np1_qdp, n0_qdp, DSSopt, rhs_multiplier, qsize, qsize_d;
} param_t;

void slave_euler_step_(param_t *param_s) {
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
  double *gl_qdp_s = param_d.qdp_s;
  double *gl_qdp_leap = param_d.qdp_leap;
  double *gl_dp = param_d.dp;
  double *gl_divdp_proj = param_d.divdp_proj;
  double *gl_Qtens_bi = param_d.Qtens_biharmonic;
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
  double dt = param_d.dt;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int np1_qdp = param_d.np1_qdp;
  int n0_qdp = param_d.n0_qdp;
  int DSSopt = param_d.DSSopt;
  int rhs_multiplier = param_d.rhs_multiplier;
  int qsize = param_d.qsize;
  int qsize_d = param_d.qsize_d;

  int istep_Qten = qsize*NLEV*NP*NP; // stripe in ie axis of Qtens_biharmonic array
  int istep_qmax = qsize*NLEV;
  int sum_k = qsize*NLEV;

  double Qdp[block];
  double dp[block];
  double divdp_proj[block];
  double dp_l[block];
  double Qtens_biharmonic[block];
  double qmax_val[KBLK];
  double qmin_val[KBLK];
  double qmax[KBLK];
  double qmin[KBLK];

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block + NLEV*NP*NP*3 + block + NLEV*2       \
        + 2*UC*NLEV);
    printf("slave_euler_step size_tol:%dk\n", size_tol/1024);
  }
#endif

  int slice_qdp = (int)(gl_qdp_leap - gl_qdp_s);  //stripe in ie axis of elem.state.Qdp
  double *src_np1_qdp = (double *)(gl_qdp_s) + (np1_qdp - 1)*qsize_d*NLEV*NP*NP;
  double *src_n0_qdp = (double *)(gl_qdp_s) + (n0_qdp - 1)*qsize_d*NLEV*NP*NP;

  double *src_qdp_ptr, *src_dp_ptr,  *src_divdp_proj_ptr, *src_Qtens_bi_ptr    \
      , *src_qmax_ptr, *src_qmin_ptr;

  int i, j, k, q, ie, k_beg, k_end, k_n, q_beg, q_end, q_n, pos;

  for (ie = 0; ie < nete - nets + 1; ie++) {
    q_beg = KBLK*id/NLEV;
    q_end = (KBLK*(id + 1))/NLEV;
    q_end = q_end < qsize ? q_end : qsize;
    k_beg = KBLK*id - q_beg*NLEV;
    k_end = KBLK*(id + 1) - q_beg*NLEV;
    k_end = k_end < NLEV ? k_end : NLEV;
    k_n = k_end - k_beg;
    q_n = (k_n == 6 && q_beg != (qsize - 1)) ? 1 : 0;
    q_n = (id*KBLK > sum_k) ? -1 : q_n;
    //if (id == 62) printf("q_beg:%d, q_end:%d, k_beg:%d, k_end:%d, k_n:%d\n", q_beg, q_end, k_beg, k_end, k_n);
    for (q = 0; q <= q_n; q++) {
      src_dp_ptr = gl_dp + ie*slice_qdp + (k_beg + q*(q*NLEV - k_beg) - q*NLEV)*NP*NP; // if k over NLEV, k should be reset pointer to start
      src_divdp_proj_ptr = gl_divdp_proj + ie*slice_qdp                        \
          + (k_beg + q*(q*NLEV - k_beg) - q*NLEV)*NP*NP; // if k over NLEV, k should be reset pointer to start
      src_qdp_ptr = src_n0_qdp + ie*slice_qdp + q_beg*NLEV*NP*NP               \
          + (k_beg + q*(q*NLEV - k_beg))*NP*NP;
      src_Qtens_bi_ptr = gl_Qtens_bi + ie*istep_Qten + q_beg*NLEV*NP*NP        \
          + (k_beg + q*(q*NLEV - k_beg))*NP*NP;
      src_qmax_ptr = gl_qmax + ie*istep_qmax + q_beg*NLEV                      \
          + (k_beg + q*(q*NLEV - k_beg));
      src_qmin_ptr = gl_qmin + ie*istep_qmax + q_beg*NLEV                      \
          + (k_beg + q*(q*NLEV - k_beg));
      //if (id == 59) printf("%d, %d, %d, %d, %d\n", q_beg, k_beg, (k_beg + q*(q*NLEV - k_beg)), q_n, k_n);
      pe_get(src_dp_ptr, dp, k_n*NP*NP*sizeof(double));
      pe_get(src_divdp_proj_ptr, divdp_proj, k_n*NP*NP*sizeof(double));
      pe_get(src_qdp_ptr, Qdp, k_n*NP*NP*sizeof(double));
      pe_get(src_Qtens_bi_ptr, Qtens_biharmonic, k_n*NP*NP*sizeof(double));
      pe_get(src_qmax_ptr, qmax, k_n*sizeof(double));
      pe_get(src_qmin_ptr, qmin, k_n*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            pos = k*NP*NP + j*NP + i;
            dp_l[pos] = dp[pos] - rhs_multiplier*dt*divdp_proj[pos];
          }
        }
      }
      for (k = 0; k < k_n; k++) {
        qmin_val[k] = +1.0e+24;
        qmax_val[k] = -1.0e+24;
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            pos = k*NP*NP + j*NP + i;
            Qtens_biharmonic[pos] = Qdp[pos]/dp_l[pos];
            qmin_val[k] = MIN(qmin_val[k], Qtens_biharmonic[pos]);
            qmax_val[k] = MAX(qmax_val[k], Qtens_biharmonic[pos]);
          }
        }
        if (rhs_multiplier == 1) {
          qmin[k] = MIN(qmin[k], qmin_val[k]);
          qmin[k] = MAX(qmin[k], 0.0);
          qmax[k] = MAX(qmax[k], qmax_val[k]);
        } else {
          qmin[k] = MAX(qmin_val[k], 0.0);
          qmax[k] = qmax_val[k];
        }
      }
      pe_put(src_Qtens_bi_ptr, Qtens_biharmonic, k_n*NP*NP*sizeof(double));
      pe_put(src_qmax_ptr, qmax, k_n*sizeof(double));
      pe_put(src_qmin_ptr, qmin, k_n*sizeof(double));
      dma_syn();
    }
  }
}
