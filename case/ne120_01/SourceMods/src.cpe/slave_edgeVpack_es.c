#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define KBLK  15
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

#define abs(value, ret) asm volatile ("fcpys $31, %1, %0" : "=r"(ret) : "r"(value))

typedef struct {
  double *v1, *v2, *spheremp, *buf;
  int *putmap, *reverse;
  int nets, nete, qsize, step_elem;
} param_t;

void slave_edgevpack_es_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();

  double *gl_v1 = param_d.v1;
  double *gl_v2 = param_d.v2;
  double *gl_spheremp = param_d.spheremp;
  double *gl_buf = param_d.buf;
  int    *gl_putmap = param_d.putmap;
  int    *gl_reverse = param_d.reverse;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;
  int step_elem = param_d.step_elem;

  double v[KBLK*NP*NP];
  //double v1[KBLK*NP*NP];
  //double v2[KBLK*NP*NP];
  double spheremp[NP*NP];
  int putmap[max_neigh_edges];
  int reverse[4];
  double buf[KBLK*NP];
  double buf_1[KBLK];

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block_Dinv + block_tensor + block_qtens     \
        + 11*NP*NP + max_neigh_edges*2 + NP*NLEV + NLEV);
    printf("slave_euler_step size_tol:%dk\n", size_tol/1024);
  }
#endif

  // local variables
  // local variables of edgeVpack
  int ll, kptr, iptr, edgeptr, is, iee, in, iw, ir;
  // local variables of laplace_sphere_wk
  double *src_v, *src_v1, *src_v2, *src_spheremp, *src_buf;
  int *src_putmap, *src_reverse;
  int sum_k = (qsize + 1)*NLEV;
  int sum_sub_k = qsize*NLEV;
  int istep_qmax = qsize*NLEV;
  int istep_qtens = qsize*NLEV*NP*NP;
  int c, r, i, j, d, l, m, n, k, q, ie, k_beg, k_end, q_beg, q_end, k_n, q_n   \
      , pos_dp, pos_qdp, pos_Qtens_bi, pos_qmax;

  for (ie = 0; ie < nete - nets + 1; ie++) {
    // roll the q and k into the id of cpe,
    // other words is that (q, k) has relationship with id
    // be careful that the edge of possbilly exception
    q_beg = KBLK*id/NLEV;
    q_end = (KBLK*(id + 1))/NLEV;
    q_end = q_end < qsize + 1 ? q_end : qsize + 1;
    k_beg = KBLK*id - q_beg*NLEV;
    k_end = KBLK*(id + 1) - q_beg*NLEV;
    k_end = k_end < NLEV ? k_end : NLEV;
    k_n = k_end - k_beg;
    q_n = (id*KBLK >= sum_k) ? -1 : 0;

    src_spheremp = gl_spheremp + ie*step_elem;
    src_putmap = gl_putmap + ie*max_neigh_edges;
    src_reverse = gl_reverse + ie*4;
    pe_get(src_spheremp, spheremp, NP*NP*sizeof(double));
    pe_get(src_putmap, putmap, max_neigh_edges*sizeof(int));
    pe_get(src_reverse, reverse, 4*sizeof(int));
    dma_syn();
    // is iee in iw is offset of edgebuf which has relationship with four direcion
    // is iee in iw only has ie dimension
    is = putmap[south - 1];
    iee = putmap[east - 1];
    in = putmap[north - 1];
    iw = putmap[west - 1];

    for (q = 0; q <= q_n; q++) {
      src_v = gl_v1 + ie*step_elem + q_beg*NLEV*NP*NP + k_beg*NP*NP;
      if (KBLK*id >= sum_sub_k) src_v = gl_v2 + ie*step_elem + k_beg*NP*NP;
      pe_get(src_v, v, k_n*NP*NP*sizeof(double));
      dma_syn();
      if (KBLK*id >= sum_sub_k) {
        for (k = 0; k < k_n; k++)
          for (j = 0; j < NP; j++)
            for (i = 0; i < NP; i++) {
              int pos_v = k*NP*NP + j*NP + i;
              int pos_spheremp = j*NP + i;
              v[pos_v] = spheremp[pos_spheremp]*v[pos_v];
            }
        pe_put(src_v, v, k_n*NP*NP*sizeof(double));
        dma_syn();
      }
      // -------------------------- start edgeVpack --------------------------//
      // ----------------------------- South ---------------------------------//
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_v = k*NP*NP + i;
          int pos_buf = k*NP + i;
          buf[pos_buf] = v[pos_v];
        }
      }
      // reverse the south
      if(reverse[south - 1]) {
        for (k = 0; k < k_n; k++) {
          for (i = 0; i < NP; i++) {
            ir = NP - 1 - i;
            int pos_buf = k*NP + ir;
            int pos_v = k*NP*NP + i;
            buf[pos_buf] = v[pos_v];
          }
        }
      }
      src_buf = gl_buf + is + q_beg*NLEV*NP + k_beg*NP;
      pe_put(src_buf, buf, k_n*NP*sizeof(double));
      dma_syn();
      // ---------------------------- East -----------------------------------//
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++)  {
          int pos_v = k*NP*NP + i*NP + NP - 1;
          int pos_buf = k*NP + i;
          buf[pos_buf] = v[pos_v];
        }
      }
      // reverse the east
      if(reverse[east - 1]) {
        for (k = 0; k < k_n; k++) {
          for (i = 0; i < NP; i++) {
            ir = NP - 1 - i;
            int pos_buf = k*NP + ir;
            int pos_v = k*NP*NP + i*NP + NP - 1;
            buf[pos_buf] = v[pos_v];
          }
        }
      }
      src_buf = gl_buf + iee + q_beg*NLEV*NP + k_beg*NP;
      pe_put(src_buf, buf, k_n*NP*sizeof(double));
      dma_syn();
      // ------------------------------ North --------------------------------//
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++)  {
          int pos_v = k*NP*NP + (NP - 1)*NP + i;
          int pos_buf = k*NP + i;
          buf[pos_buf] = v[pos_v];
        }
      }
      // reverse the north
      if(reverse[north - 1]) {
        for (k = 0; k < k_n; k++) {
          for (i = 0; i < NP; i++) {
            ir = NP - 1 - i;
            int pos_buf = k*NP + ir;
            int pos_v = k*NP*NP + (NP - 1)*NP + i;
            buf[pos_buf] = v[pos_v];
          }
        }
      }
      src_buf = gl_buf + in + q_beg*NLEV*NP + k_beg*NP;
      pe_put(src_buf, buf, k_n*NP*sizeof(double));
      dma_syn();
      // ------------------------------ West ---------------------------------//
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++)  {
          int pos_v = k*NP*NP + i*NP;
          int pos_buf = k*NP + i;
          buf[pos_buf] = v[pos_v];
        }
      }
      // reverse the west
      if(reverse[west - 1]) {
        for (k = 0; k < k_n; k++) {
          for (i = 0; i < NP; i++) {
            ir = NP - 1 - i;
            int pos_buf = k*NP + ir;
            int pos_v = k*NP*NP + i*NP;
            buf[pos_buf] = v[pos_v];
          }
        }
      }
      src_buf = gl_buf + iw + q_beg*NLEV*NP + k_beg*NP;
      pe_put(src_buf, buf, k_n*NP*sizeof(double));
      dma_syn();
      // ---------------------------- SWEST -------------------------------- //
      for (ll = swest - 1; ll < (swest + max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          for (k = 0; k < k_n; k++) {
            int pos_v = k*NP*NP;
            buf_1[k] = v[pos_v];
          }
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV+ k_beg;
          pe_put(src_buf, buf_1, k_n*sizeof(double));
          dma_syn();
        }
      }  // end loop ll
      // ----------------------------- SEAST -------------------------------- //
      for (ll = swest + max_corner_elem - 1; ll < (swest + 2*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          for (k = 0; k < k_n; k++) {
            int pos_v = k*NP*NP + (NP - 1);
            buf_1[k] = v[pos_v];
          }
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV+ k_beg;
          pe_put(src_buf, buf_1, k_n*sizeof(double));
          dma_syn();
        }
      }  // end loop ll
      // ----------------------------- NEAST -------------------------------- //
      for (ll = swest + 3*max_corner_elem - 1; ll < (swest + 4*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          for (k = 0; k < k_n; k++) {
            int pos_v = k*NP*NP + (NP - 1)*NP + (NP - 1);
            buf_1[k] = v[pos_v];
          }
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV+ k_beg;
          pe_put(src_buf, buf_1, k_n*sizeof(double));
          dma_syn();
        }
      }  // end loop ll
      // ----------------------- NWEST -------------------------------- //
      for (ll = swest + 2*max_corner_elem - 1; ll < (swest + 3*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          for (k = 0; k < k_n; k++) {
            int pos_v = k*NP*NP + (NP - 1)*NP;
            buf_1[k] = v[pos_v];
          }
          src_buf = gl_buf + putmap[ll] + q_beg*NLEV+ k_beg;
          pe_put(src_buf, buf_1, k_n*sizeof(double));
          dma_syn();
        }
      }  // end loop ll
      // -------------------- end of edgeVpack -------------------------//
    } // end loop q
  } // end loop ie

#if 0
  if (id == 0)
    printf("qsize:%d, var_coef:%d\n", qsize, var_coef);
#endif
}
