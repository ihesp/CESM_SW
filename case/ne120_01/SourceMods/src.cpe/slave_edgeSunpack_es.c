#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define KBLK  12
#define west  1
#define east  2
#define south 3
#define north 4
#define swest 5
#define seast 6
#define nwest 7
#define neast 8
#define max_neigh_edges   8
#define max_corner_elem   1

#define MAX(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define MIN(val1, val2) ((val1) < (val2) ? (val1) : (val2))

#define MAX5(a, a1, a2, a3, a4) { \
  double max = a; \
  max = MAX(max, a); \
  max = MAX(max, a1); \
  max = MAX(max, a2); \
  max = MAX(max, a3); \
  max = MAX(max, a4); \
  a = max; \
}

#define MIN5(a, a1, a2, a3, a4) { \
  double min = a; \
  min = MIN(min, a); \
  min = MIN(min, a1); \
  min = MIN(min, a2); \
  min = MIN(min, a3); \
  min = MIN(min, a4); \
  a = min; \
}

typedef struct {
  double *qmax, *qmin, *receive;
  int *getmap;
  int nets, nete, qsize;
} param_t;

void slave_edgesunpack_es_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
  double *gl_receive = param_d.receive;
  int *gl_getmap = param_d.getmap;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;

  double qmax[KBLK];
  double qmin[KBLK];
  double receive[KBLK];
  double receive_is[KBLK];
  double receive_ie[KBLK];
  double receive_in[KBLK];
  double receive_iw[KBLK];
  int getmap[max_neigh_edges];

  int sum_k = qsize*NLEV;
  int istep_qmax = qsize*NLEV;
  
  int k, n, ie, q, k_beg, k_end, k_n, q_beg, q_end, q_n, pos, kptr, iptr, ll    \
      , is, iee, in, iw, ir;

  double *src_receive, *src_qmax, *src_qmin;
  int *src_getmap;

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
      src_getmap = gl_getmap + ie*max_neigh_edges;
      src_qmax = gl_qmax + ie*istep_qmax + q_beg*NLEV                          \
          + (k_beg + q*(q*NLEV - k_beg)); // (k_beg + q*(q*NLEV - k_beg)) could be 0 12 24 0 6 18 0
      src_qmin = gl_qmin + ie*istep_qmax + q_beg*NLEV                          \
          + (k_beg + q*(q*NLEV - k_beg));
      pe_get(src_getmap, getmap, max_neigh_edges*sizeof(int));
      pe_get(src_qmax, qmax, k_n*sizeof(double));
      pe_get(src_qmin, qmin, k_n*sizeof(double));
      dma_syn();

      is  = getmap[south - 1];
      iee = getmap[east - 1];
      in  = getmap[north - 1];
      iw  = getmap[west - 1];

      src_receive = gl_receive + iee + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
      pe_get(src_receive, receive_ie, k_n*sizeof(double));
      dma_syn();
      src_receive = gl_receive + is + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
      pe_get(src_receive, receive_is, k_n*sizeof(double));
      dma_syn();
      src_receive = gl_receive + iw + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
      pe_get(src_receive, receive_iw, k_n*sizeof(double));
      dma_syn();
      src_receive = gl_receive + in + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
      pe_get(src_receive, receive_in, k_n*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        MIN5(qmin[k], receive_ie[k], receive_is[k], receive_iw[k], receive_in[k]);
      }
      src_receive = gl_receive + iee + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
      pe_get(src_receive, receive_ie, k_n*sizeof(double));
      dma_syn();
      src_receive = gl_receive + is + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
      pe_get(src_receive, receive_is, k_n*sizeof(double));
      dma_syn();
      src_receive = gl_receive + iw + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
      pe_get(src_receive, receive_iw, k_n*sizeof(double));
      dma_syn();
      src_receive = gl_receive + in + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
      pe_get(src_receive, receive_in, k_n*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        MAX5(qmax[k], receive_ie[k], receive_is[k], receive_iw[k], receive_in[k]);
      }
      // ---------------------------- SWEST ----------------------------------//
      for (ll = swest - 1; ll < (swest + max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV                   \
              + (k_beg + q*(q*NLEV - k_beg));
          pe_get(src_receive, receive, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            qmin[k] = MIN(qmin[k], receive[k]);
          }
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV                   \
              + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
          pe_get(src_receive, receive, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            qmax[k] = MAX(qmax[k], receive[k]);
          }
        }
      }  // end loop ll
      // ---------------------------- SEAST ----------------------------------//
      for (ll = swest + max_corner_elem - 1; ll < (swest + 2*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV                   \
              + (k_beg + q*(q*NLEV - k_beg));
          pe_get(src_receive, receive, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            qmin[k] = MIN(qmin[k], receive[k]);
          }
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV                   \
              + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
          pe_get(src_receive, receive, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            qmax[k] = MAX(qmax[k], receive[k]);
          }
        }
      } // end loop ll
      // ----------------------- NEAST -------------------------------- //
      for (ll = swest + 3*max_corner_elem - 1; ll < (swest + 4*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV                   \
              + (k_beg + q*(q*NLEV - k_beg));
          pe_get(src_receive, receive, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            qmin[k] = MIN(qmin[k], receive[k]);
          }
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV                   \
              + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
          pe_get(src_receive, receive, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            qmax[k] = MAX(qmax[k], receive[k]);
          }
        }
      } // end loop ll
      // ----------------------- NWEST -------------------------------- //
      for (ll = swest + 2*max_corner_elem - 1; ll < (swest + 3*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV                   \
              + (k_beg + q*(q*NLEV - k_beg));
          pe_get(src_receive, receive, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            qmin[k] = MIN(qmin[k], receive[k]);
          }
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV                   \
              + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
          pe_get(src_receive, receive, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            qmax[k] = MAX(qmax[k], receive[k]);
          }
        }
      } // end loop ll
      for (k = 0; k < k_n; k++) {
        qmin[k] = MAX(qmin[k], 0.0);
      }
      pe_put(src_qmin, qmin, k_n*sizeof(double));
      pe_put(src_qmax, qmax, k_n*sizeof(double));
      dma_syn();
    }
  }
}
