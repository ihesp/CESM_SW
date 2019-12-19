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

#define abs(value, ret) asm volatile ("fcpys $31, %1, %0" : "=r"(ret) : "r"(value))

typedef struct {
  double *Dinv, *rspheremp, *spheremp, *variable_hyper, *tensorVisc, *qtens, *Dvv, *receive;
  int *getmap;
  double rrearth, hypervis_scaling, hypervis_power;
  int nets, nete, qsize, step_elem;
} param_t;

void slave_biharmonic_wk_scalar_2_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();

  double *gl_Dinv = param_d.Dinv;
  double *gl_rspheremp = param_d.rspheremp;
  double *gl_spheremp = param_d.spheremp;
  double *gl_variable_hyper = param_d.variable_hyper;
  double *gl_tensorVisc = param_d.tensorVisc;
  double *gl_qtens = param_d.qtens;
  double *gl_Dvv = param_d.Dvv;
  double *gl_receive = param_d.receive;
  int    *gl_getmap = param_d.getmap;
  double rrearth = param_d.rrearth;
  double hypervis_scaling = param_d.hypervis_scaling;
  double hypervis_power = param_d.hypervis_power;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;
  int step_elem = param_d.step_elem;

  double Dinv[2*2*NP*NP];
  double Dvv[NP*NP];
  double variable_hyper[NP*NP];
  double tensorVisc[2*2*NP*NP];
  double qtens[KBLK*NP*NP];
  double rspheremp[NP*NP];
  double spheremp[NP*NP];
  double lap_p[NP*NP];
  double grads[2*NP*NP];
  double v1[NP*NP];
  double v2[NP*NP];
  double vtemp[NP*NP*2];
  int getmap[max_neigh_edges];
  double receive[KBLK*NP];
  double receive_1[KBLK];

  double dsdx00, dsdy00; int var_coef = 1;
  double oldgrads_1, oldgrads_2, tensor_1, tensor_2;
  if (hypervis_scaling > 0) var_coef = 0;

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block_Dinv + block_tensor + block_qtens     \
        + 11*NP*NP);
    printf("slave_euler_step size_tol:%dk\n", size_tol/1024);
  }
#endif

  pe_get(gl_Dvv, Dvv, NP*NP*sizeof(double));
  dma_syn();
  // local variables
  // local variables of edgeVunpack
  int ll, kptr, iptr, edgeptr, is, iee, in, iw, ir;
  // local variables of laplace_sphere_wk
  double *src_Dinv, *src_rspheremp, *src_spheremp, *src_variable_hyper,        \
      *src_tensorVisc, *src_qtens, *src_receive;
  int *src_getmap;

  int sum_k = qsize*NLEV;
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
    q_end = q_end < qsize ? q_end : qsize;
    k_beg = KBLK*id - q_beg*NLEV;
    k_end = KBLK*(id + 1) - q_beg*NLEV;
    k_end = k_end < NLEV ? k_end : NLEV;
    k_n = k_end - k_beg;
    q_n = (k_n == 6 && q_beg != (qsize - 1)) ? 1 : 0;
    q_n = (id*KBLK > sum_k) ? -1 : q_n;

    src_Dinv = gl_Dinv + ie*step_elem;
    src_rspheremp = gl_rspheremp + ie*step_elem;
    src_spheremp = gl_spheremp + ie*step_elem;
    src_variable_hyper = gl_variable_hyper + ie*step_elem;
    src_tensorVisc = gl_tensorVisc + ie*step_elem;
    src_getmap = gl_getmap + ie*max_neigh_edges;
    pe_get(src_Dinv, Dinv, 2*2*NP*NP*sizeof(double));
    pe_get(src_rspheremp, rspheremp, NP*NP*sizeof(double));
    pe_get(src_spheremp, spheremp, NP*NP*sizeof(double));
    pe_get(src_tensorVisc, tensorVisc, 2*2*NP*NP*sizeof(double));
    pe_get(src_variable_hyper, variable_hyper, NP*NP*sizeof(double));
    pe_get(src_getmap, getmap, max_neigh_edges*sizeof(int));
    dma_syn();
    // is iee in iw is offset of edgebuf which has relationship with four direcion
    // is iee in iw only has ie dimension
    is = getmap[south - 1];
    iee = getmap[east - 1];
    in = getmap[north - 1];
    iw = getmap[west - 1];

    for (q = 0; q <= q_n; q++) {
      src_qtens = gl_qtens + ie*istep_qtens + q_beg*NLEV*NP*NP                 \
          + (k_beg + q*(q*NLEV - k_beg))*NP*NP;
      pe_get(src_qtens, qtens, k_n*NP*NP*sizeof(double));
      dma_syn();
      // ------------------------- start of edgeVunpack ----------------------//
      // ------------------------------- East --------------------------------//
      src_receive = gl_receive + iee + q_beg*NLEV*NP + (k_beg + q*(q*NLEV - k_beg))*NP;
      pe_get(src_receive, receive, k_n*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_qtens = k*NP*NP + i*NP + (NP - 1);
          int pos_receive = k*NP + i;
          qtens[pos_qtens] = qtens[pos_qtens] + receive[pos_receive];
        }
      }
      // ----------------------------- South ---------------------------------//
      src_receive = gl_receive + is + q_beg*NLEV*NP + (k_beg + q*(q*NLEV - k_beg))*NP;
      pe_get(src_receive, receive, k_n*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_qtens = k*NP*NP + i;
          int pos_receive = k*NP + i;
          qtens[pos_qtens] = qtens[pos_qtens] + receive[pos_receive];
        }
      }
      // ----------------------------- North ---------------------------------//
      src_receive = gl_receive + in + q_beg*NLEV*NP + (k_beg + q*(q*NLEV - k_beg))*NP;
      pe_get(src_receive, receive, k_n*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_qtens = k*NP*NP + (NP - 1)*NP + i;
          int pos_receive = k*NP + i;
          qtens[pos_qtens] = qtens[pos_qtens] + receive[pos_receive];
        }
      }
      // ----------------------------- West ----------------------------------//
      src_receive = gl_receive + iw + q_beg*NLEV*NP + (k_beg + q*(q*NLEV - k_beg))*NP;
      pe_get(src_receive, receive, k_n*NP*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (i = 0; i < NP; i++) {
          int pos_qtens = k*NP*NP + i*NP;
          int pos_receive = k*NP + i;
          qtens[pos_qtens] = qtens[pos_qtens] + receive[pos_receive];
        }
      }
      // ----------------------- SWEST -------------------------------------- //
      for (ll = swest - 1; ll < (swest + max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
          pe_get(src_receive, receive_1, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            int pos_qtens = k*NP*NP;
            qtens[pos_qtens] = qtens[pos_qtens] + receive_1[k];
          }
        }
      }  // end loop ll
      // ----------------------- SEAST -------------------------------------- //
      for (ll = swest + max_corner_elem - 1; ll < (swest + 2*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
          pe_get(src_receive, receive_1, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            int pos_qtens = k*NP*NP + NP - 1;
            qtens[pos_qtens] = qtens[pos_qtens] + receive_1[k];
          }
        }
      }  // end loop ll
      // ----------------------- NEAST -------------------------------------- //
      for (ll = swest + 3*max_corner_elem - 1; ll < (swest + 4*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
          pe_get(src_receive, receive_1, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            int pos_qtens = k*NP*NP + NP*(NP - 1) + NP - 1;
            qtens[pos_qtens] = qtens[pos_qtens] + receive_1[k];
          }
        }
      }  // end loop ll
      // ----------------------------- NWEST -------------------------------- //
      for (ll = swest + 2*max_corner_elem - 1; ll < (swest + 3*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q_beg*NLEV + (k_beg + q*(q*NLEV - k_beg));
          pe_get(src_receive, receive_1, k_n*sizeof(double));
          dma_syn();
          for (k = 0; k < k_n; k++) {
            int pos_qtens = k*NP*NP + NP*(NP - 1);
            qtens[pos_qtens] = qtens[pos_qtens] + receive_1[k];
          }
        }
      }  // end loop ll ----- end of edgeVunpack ----
      // --------------- start gradient_sphere divergence_sphere -------------//
      for (k = 0; k < k_n; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos = j*NP + i;
            int pos_qtens = k*NP*NP + j*NP + i;
            lap_p[pos] = rspheremp[pos]*qtens[pos_qtens];
          }
        }
        for (j = 0; j < NP; j++) {  // start gradient_sphere
          for (l = 0; l < NP; l++) {
            dsdx00 = 0.0;
            dsdy00 = 0.0;
            for (i = 0; i < NP; i++) {
              int pos_Dvv = l*NP + i;
              int pos_lap_1 = j*NP + i;
              int pos_lap_2 = i*NP + j;
              dsdx00 = dsdx00 + Dvv[pos_Dvv]*lap_p[pos_lap_1];
              dsdy00 = dsdy00 + Dvv[pos_Dvv]*lap_p[pos_lap_2];
            }
            int pos_v1 = j*NP + l;
            int pos_v2 = l*NP + j;
            v1[pos_v1] = dsdx00*rrearth;
            v2[pos_v2] = dsdy00*rrearth;
          }
        }
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos_grads_1 = j*NP + i;
            int pos_grads_2 = NP*NP + j*NP + i;
            int pos_v = j*NP + i;
            int pos_Dinv_1 = j*NP + i;
            int pos_Dinv_2 = NP*NP + j*NP + i;
            int pos_Dinv_3 = 2*NP*NP + j*NP + i;
            int pos_Dinv_4 = 3*NP*NP + j*NP + i;
            grads[pos_grads_1] = Dinv[pos_Dinv_1]*v1[pos_v]                    \
                + Dinv[pos_Dinv_2]*v2[pos_v];
            grads[pos_grads_2] = Dinv[pos_Dinv_3]*v1[pos_v]                    \
                + Dinv[pos_Dinv_4]*v2[pos_v];
          }
        }
        if (var_coef) {
          //if (id == 0)
          //  printf("var_coef_in:%d\n", var_coef);
          if (hypervis_power != 0) {
            for (j = 0; j < NP; j++) {
              for (i = 0; i < NP; i++) {
                int pos_grads_1 = j*NP + i;
                int pos_grads_2 = NP*NP + j*NP + i;
                int pos = j*NP + i;
                grads[pos_grads_1] = grads[pos_grads_1]*variable_hyper[pos];
                grads[pos_grads_2] = grads[pos_grads_2]*variable_hyper[pos];
              }
            }
          } else if (hypervis_scaling != 0) {
            for (j = 0; j < NP; j++) {
              for (i = 0; i < NP; i++) {
                int pos_grads_1 = j*NP + i;
                int pos_grads_2 = NP*NP + j*NP + i;
                int pos_tensor_1 = j*NP + i;
                int pos_tensor_2 = NP*NP + j*NP + i;
                int pos_tensor_3 = 2*NP*NP + j*NP + i;
                int pos_tensor_4 = 3*NP*NP + j*NP + i;
                oldgrads_1 = grads[pos_grads_1];
                oldgrads_2 = grads[pos_grads_2];
                grads[pos_grads_1] = oldgrads_1*tensorVisc[pos_tensor_1]       \
                    + oldgrads_2*tensorVisc[pos_tensor_3];
                grads[pos_grads_2] = oldgrads_1*tensorVisc[pos_tensor_2]       \
                    + oldgrads_2*tensorVisc[pos_tensor_4];
              }
            }
          }
        }  // end of gradient_sphere
        // start divergence_sphere_wk
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos_1 = j*NP + i;
            int pos_2 = NP*NP + j*NP + i;
            int pos_Dinv_1 = j*NP + i;
            int pos_Dinv_2 = NP*NP + j*NP + i;
            int pos_Dinv_3 = 2*NP*NP + j*NP + i;
            int pos_Dinv_4 = 3*NP*NP + j*NP + i;
            vtemp[pos_1] = Dinv[pos_Dinv_1]*grads[pos_1]                       \
                + Dinv[pos_Dinv_3]*grads[pos_2];
            vtemp[pos_2] = Dinv[pos_Dinv_2]*grads[pos_1]                       \
                + Dinv[pos_Dinv_4]*grads[pos_2];
          }
        }
        for (n = 0; n < NP; n++) {
          for (m = 0; m < NP; m++) {
            int pos_qtens = k*NP*NP + n*NP + m;
            qtens[pos_qtens] = 0;
            for (j = 0; j < NP; j++) {
              int pos_spheremp_1 = n*NP + j;
              int pos_spheremp_2 = j*NP + m;
              int pos_vtemp_1 = n*NP + j;
              int pos_vtemp_2 = NP*NP + j*NP + m;
              int pos_Dvv_1 = j*NP + m;
              int pos_Dvv_2 = j*NP + n;
              qtens[pos_qtens] = qtens[pos_qtens]
                  - (spheremp[pos_spheremp_1]*vtemp[pos_vtemp_1]*Dvv[pos_Dvv_1]  \
                  + spheremp[pos_spheremp_2]*vtemp[pos_vtemp_2]*Dvv[pos_Dvv_2])  \
                  *rrearth;
            }
          }
        }  // end divergence_sphere_wk
      }  // end loop k
      pe_put(src_qtens, qtens, k_n*NP*NP*sizeof(double));
      dma_syn();
    }  // end loop q
  }  // end loop ie
#if 0
  if (id == 0)
    printf("qsize:%d, var_coef:%d\n", qsize, var_coef);
#endif
}
