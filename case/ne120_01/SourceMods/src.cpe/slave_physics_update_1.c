#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"

#define KBLK 12

#define MAX(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define MIN(val1, val2) ((val1) < (val2) ? (val1) : (val2))

typedef struct {
  double *state_q, *ptend_q, *qmin, *worst;
  int *nvals, *iw, *kw, *found, *ptend_lq;
  double dt;
  int pcnst, top_level, bot_level, lchnk, psetcols, pver, ncol, ixnumice   \
      , ixnumliq, ixnumsnow, ixnumrain;
} param_t;

void slave_physics_update_1_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();
  double *gl_state_q = param_d.state_q;
  double *gl_ptend_q = param_d.ptend_q;
  double *qmin = param_d.qmin;
  double *worst = param_d.worst;
  int *nvals = param_d.nvals;
  int *iw = param_d.iw;
  int *kw = param_d.kw;
  int *found = param_d.found;
  int *ptend_lq = param_d.ptend_lq;
  double dt = param_d.dt;
  int pcnst = param_d.pcnst;
  int top_level = param_d.top_level;
  int bot_level = param_d.bot_level;
  int lchnk = param_d.lchnk;
  int psetcols = param_d.psetcols;
  int pver = param_d.pver;
  int ncol = param_d.ncol;
  int ixnumice = param_d.ixnumice;
  int ixnumliq = param_d.ixnumliq;
  int ixnumsnow = param_d.ixnumsnow;
  int ixnumrain = param_d.ixnumrain;

  double state_q[pver*psetcols];
  double ptend_q[pver*psetcols];

  int i, k;


#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block + NLEV*NP*NP*3 + block + NLEV*2       \
        + 2*UC*NLEV);
    printf("slave_euler_step size_tol:%dk\n", size_tol/1024);
  }
#endif

  if (id == -1) {
    printf("bot_level:%d, ncol:%d, dt:%lf, ixnumrain:%d\n",bot_level, ncol, dt, ixnumrain);
  }

  if (id < pcnst) {
    int m = id;
    double *src_state_q = gl_state_q + m*pver*psetcols;
    if (ptend_lq[m] == 1) {
      pe_get(gl_state_q, state_q, pver*psetcols*sizeof(double));
      pe_get(gl_ptend_q, ptend_q, pver*psetcols*sizeof(double));
      dma_syn();
      for (k = top_level - 1; k < bot_level; k++) {
        for (i = 0; i < ncol; i++) {
          int pos = k*psetcols + i;
          state_q[pos] = state_q[pos] + ptend_q[pos]*dt;
        }
      }
      if (m != ixnumice && m != ixnumliq && m != ixnumrain && m != ixnumsnow) {
        if (id == 0)
        qneg3_f_(&lchnk, &ncol, &psetcols, &pver, &m, &m, &qmin[m], state_q \
            , &nvals[m], &iw[m], &kw[m], &found[m], &worst[m]);
      } else {
        for (k = top_level -1; k < bot_level; k++) {
          for (i = 0; i < ncol; i++) {
            int pos = k*psetcols + i;
            state_q[pos] = MAX(1.0e-12, state_q[pos]);
            state_q[pos] = MIN(1.0e10, state_q[pos]);
          }
        }
      }
      pe_put(gl_state_q, state_q, pver*psetcols*sizeof(double));
      pe_put(gl_ptend_q, ptend_q, pver*psetcols*sizeof(double));
      dma_syn();
    }
  }
}
