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
  double *v1, *v2, *rspheremp, *receive;
  int *getmap;
  int nets, nete, qsize, step_elem;
} param_t;

void slave_edgevunpack_es_(param_t *param_s) {
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
  double *gl_rspheremp = param_d.rspheremp;
  double *gl_receive = param_d.receive;
  int    *gl_getmap = param_d.getmap;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;
  int step_elem = param_d.step_elem;

  double v[KBLK*NP*NP];
  double rspheremp[NP*NP];
  int getmap[max_neigh_edges];
  double receive[KBLK*NP];
  double receive_1[KBLK];

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
  double *src_v, *src_v1, *src_v2, *src_rspheremp, *src_receive;
  int *src_getmap;
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

    src_rspheremp = gl_rspheremp + ie*step_elem;
    src_getmap = gl_getmap + ie*max_neigh_edges;
    pe_get(src_rspheremp, rspheremp, NP*NP*sizeof(double));
    pe_get(src_getmap, getmap, max_neigh_edges*sizeof(int));
    dma_syn();
    // is iee in iw is offset of edgebuf which has relationship with four direcion
    // is iee in iw only has ie dimension
    is = getmap[south - 1];
    iee = getmap[east - 1];
    in = getmap[north - 1];
    iw = getmap[west - 1];

    for (q = 0; q <= q_n; q++) {
      src_v = gl_v1 + ie*step_elem + q_beg*NLEV*NP*NP + k_beg*NP*NP;
      if (KBLK*id >= sum_sub_k) src_v = gl_v2 + ie*step_elem + k_beg*NP*NP;
      pe_get(src_v, v, k_n*NP*NP*sizeof(double));
      dma_syn();
      // ------------------------- start of edgeVunpack ----------------------//
      // ------------------------------- East --------------------------------//
      src_receive = gl_receive + iee + q_beg*NLEV*NP + k_beg*NP;
      pe_get(src_receive, receive, k_n*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_v = k*NP*NP + i*NP + (NP - 1);
          int pos_receive = k*NP + i;
          v[pos_v] = v[pos_v] + receive[pos_receive];
        }
      }
      // ----------------------------- South ---------------------------------//
      src_receive = gl_receive + is + q_beg*NLEV*NP + k_beg*NP;
      pe_get(src_receive, receive, k_n*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_v = k*NP*NP + i;
          int pos_receive = k*NP + i;
          v[pos_v] = v[pos_v] + receive[pos_receive];
        }
      }
      // ----------------------------- North ---------------------------------//
      src_receive = gl_receive + in + q_beg*NLEV*NP + k_beg*NP;
      pe_get(src_receive, receive, k_n*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_v = k*NP*NP + (NP - 1)*NP + i;
          int pos_receive = k*NP + i;
          v[pos_v] = v[pos_v] + receive[pos_receive];
        }
      }
      // ----------------------------- West ----------------------------------//
      src_receive = gl_receive + iw + q_beg*NLEV*NP + k_beg*NP;
      pe_get(src_receive, receive, k_n*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_v = k*NP*NP + i*NP;
          int pos_receive = k*NP + i;
          v[pos_v] = v[pos_v] + receive[pos_receive];
        }
      }
      // ----------------------- SWEST -------------------------------------- //
      for (ll = swest - 1; ll < (swest + max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV + k_beg;
          pe_get(src_receive, receive_1, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            int pos_v = k*NP*NP;
            v[pos_v] = v[pos_v] + receive_1[k];
          }
        }
      }  // end loop ll
      // ----------------------- SEAST -------------------------------------- //
      for (ll = swest + max_corner_elem - 1; ll < (swest + 2*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV + k_beg;
          pe_get(src_receive, receive_1, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            int pos_v = k*NP*NP + NP - 1;
            v[pos_v] = v[pos_v] + receive_1[k];
          }
        }
      }  // end loop ll
      // ----------------------- NEAST -------------------------------------- //
      for (ll = swest + 3*max_corner_elem - 1; ll < (swest + 4*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV + k_beg;
          pe_get(src_receive, receive_1, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            int pos_v = k*NP*NP + NP*(NP - 1) + NP - 1;
            v[pos_v] = v[pos_v] + receive_1[k];
          }
        }
      }  // end loop ll
      // ----------------------------- NWEST -------------------------------- //
      for (ll = swest + 2*max_corner_elem - 1; ll < (swest + 3*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV + k_beg;
          pe_get(src_receive, receive_1, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            int pos_v = k*NP*NP + NP*(NP - 1);
            v[pos_v] = v[pos_v] + receive_1[k];
          }
        }
      }  // end loop ll ----- end of edgeVunpack ----
      for (k = 0; k < k_n; k++)
        for (j = 0; j < NP; j++)
          for (i = 0; i < NP; i++) {
            int pos_v = k*NP*NP + j*NP + i;
            int pos_spheremp = j*NP + i;
            v[pos_v] = rspheremp[pos_spheremp]*v[pos_v];
          }
      pe_put(src_v, v, k_n*NP*NP*sizeof(double));
      dma_syn();
    } // end loop q
  } // end loop ie
}
