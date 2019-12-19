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

typedef struct {
  double *qmax, *qmin, *buf;
  int *putmap;
  int nets, nete, qsize;
} param_t;

void slave_edgespack_es_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
  double *gl_buf = param_d.buf;
  int *gl_putmap = param_d.putmap;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;

  double qmax[KBLK];
  double qmin[KBLK];
  int putmap[max_neigh_edges];

  int i, j, k, q, n, ie, k_beg, k_end, k_n, q_beg, q_end, q_n, pos, id_map     \
      , kptr, iptr, ll, is, iee, in, iw, ir;
  int sum_k = qsize*NLEV;
  int istep_qmax = qsize*NLEV;

  double *src_buf, *src_qmax, *src_qmin;
  int *src_putmap;

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
      src_putmap = gl_putmap + ie*max_neigh_edges;
      src_qmax = gl_qmax + ie*istep_qmax + q_beg*NLEV                          \
          + (k_beg + q*(q*NLEV - k_beg)); // (k_beg + q*(q*NLEV - k_beg)) could be 0 12 24 0 6 18 0
      src_qmin = gl_qmin + ie*istep_qmax + q_beg*NLEV                          \
          + (k_beg + q*(q*NLEV - k_beg));
      pe_get(src_putmap, putmap, max_neigh_edges*sizeof(int));
      pe_get(src_qmax, qmax, k_n*sizeof(double));
      pe_get(src_qmin, qmin, k_n*sizeof(double));
      dma_syn();
      is = putmap[south - 1];
      iee = putmap[east - 1];
      in = putmap[north - 1];
      iw = putmap[west - 1];
      src_buf = gl_buf + iee + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
      pe_put(src_buf, qmin, k_n*sizeof(double)); // East for qmin
      dma_syn();
      src_buf = gl_buf + iee + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
      pe_put(src_buf, qmax, k_n*sizeof(double)); // East for qmax
      dma_syn();
      src_buf = gl_buf + is + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
      pe_put(src_buf, qmin, k_n*sizeof(double)); // South for qmin
      dma_syn();
      src_buf = gl_buf + is + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
      pe_put(src_buf, qmax, k_n*sizeof(double)); // South for qmax
      dma_syn();
      src_buf = gl_buf + in + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
      pe_put(src_buf, qmin, k_n*sizeof(double)); // North for qmin
      dma_syn();
      src_buf = gl_buf + in + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
      pe_put(src_buf, qmax, k_n*sizeof(double)); // North for qmax
      dma_syn();
      src_buf = gl_buf + iw + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
      pe_put(src_buf, qmin, k_n*sizeof(double)); // West for qmin
      dma_syn();
      src_buf = gl_buf + iw + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
      pe_put(src_buf, qmax, k_n*sizeof(double)); // West for qmax
      dma_syn();
      // ---------------------------- SWEST ----------------------------------//
      for (ll = swest - 1; ll < (swest + max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
          pe_put(src_buf, qmin, k_n*sizeof(double));
          dma_syn();
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
          pe_put(src_buf, qmax, k_n*sizeof(double));
          dma_syn();
        }
      }  // end loop ll
      // ---------------------------- SEAST ----------------------------------//
      for (ll = swest + max_corner_elem - 1; ll < (swest + 2*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
          pe_put(src_buf, qmin, k_n*sizeof(double));
          dma_syn();
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
          pe_put(src_buf, qmax, k_n*sizeof(double));
          dma_syn();
        }
      }
      // ----------------------------- NEAST -------------------------------- //
      for (ll = swest + 3*max_corner_elem - 1; ll < (swest + 4*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
          pe_put(src_buf, qmin, k_n*sizeof(double));
          dma_syn();
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
          pe_put(src_buf, qmax, k_n*sizeof(double));
          dma_syn();
        }
      } // end loop ll
      // ----------------------- NWEST -------------------------------- //
      for (ll = swest + 2*max_corner_elem - 1; ll < (swest + 3*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
          pe_put(src_buf, qmin, k_n*sizeof(double));
          dma_syn();
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg)) + qsize*NLEV;
          pe_put(src_buf, qmax, k_n*sizeof(double));
          dma_syn();
        }
      } // end loop ll
    } // end loop q
  } // end loop ie
}
