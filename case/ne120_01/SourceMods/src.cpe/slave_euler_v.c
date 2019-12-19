
#include <slave.h>
#include "dma_macros.h"
//#include "ldm_alloc.h"

#define NP    4
#define NLEV  30
#define KBLK  12
#define maxiter           (NP*NP - 1)
#define stripe_qdp        (NLEV*NP*NP)
#define istep_dp          (NLEV*NP*NP)
#define qstep_Qten        (NLEV*NP*NP)       // stripe in q axis of Qtens_biharmonic array
#define qstep_dp_star     (NLEV*NP*NP)
#define qstep_qmax        NLEV
#define qstep_Qtens       (NLEV*NP*NP)
#define block             (KBLK*NP*NP)
#define block_dp          (KBLK*NP*NP)
#define block_vn0         (KBLK*2*NP*NP)
#define block_gradQ       (2*NP*NP)
#define block_dp_star     (KBLK*NP*NP)
#define block_Qtens       (KBLK*NP*NP)
#define block_Dinv        (2*2*NP*NP)
#define block_det         (NP*NP)

#define abs(value, ret) asm volatile ("fcpys $31, %1, %0" : "=r"(ret) : "r"(value))

#define sum_array(sum, array, len) {   \
  int i;   \
  sum = 0.0;   \
  for (i = 0; i < len; i++)     \
    sum = sum + array[i];  \
}

#define sum_array_multiply_1(sum, a, b, len, type) { \
  int i;            \
  sum = 0.0;         \
  type _array[len];    \
  for (i = 0; i < len; i++) { \
    _array[i] = a[i] * b[i];   \
    sum = sum + _array[i]; \
  }  \
}

#define sum_array_multiply(sum, a, b) { \
  int i;  \
  double _array[16];   \
  sum = 0;     \
  for (i = 0; i < 16; i++) {  \
    _array[i] = a[i] * b[i]; \
    sum = sum + _array[i];  \
  } \
}

typedef struct {
  double *gl_qdp, *gl_qdp_leap, *divdp_proj, *dp, *vn0, *Dvv, *Dinv, *metdet,  \
      *rmetdet, *Qtens_biharmonic, *divdp, *dpdiss_biharmonic,                 \
      *spheremp, *qmax, *qmin;
  double dt, rrearth, nu_p, nu_q;
  int nets, nete, rhs_multiplier, qsize, qsize_d, n0_qdp, np1_qdp, limiter_option       \
      , rhs_viss;
} param_t;

void slave_euler_v_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
	dma_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();

  double *gl_qdp = param_d.gl_qdp;
  double *gl_qdp_leap = param_d.gl_qdp_leap;
  double *gl_divdp_proj = param_d.divdp_proj;
  double *gl_dp = param_d.dp;
  double *gl_vn0 = param_d.vn0;
  double *gl_Dvv = param_d.Dvv;
  double *gl_Dinv = param_d.Dinv;
  double *gl_metdet = param_d.metdet;
  double *gl_rmetdet = param_d.rmetdet;
  double *gl_Qtens_biharmonic = param_d.Qtens_biharmonic;
  double *gl_divdp = param_d.divdp;
  double *gl_dpdiss_biharmonic = param_d.dpdiss_biharmonic;
  double *gl_spheremp = param_d.spheremp;
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
  double dt = param_d.dt;
  double rrearth = param_d.rrearth;
  double nu_p = param_d.nu_p;
  double nu_q = param_d.nu_q;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int rhs_multiplier = param_d.rhs_multiplier;
  int qsize = param_d.qsize;
  int qsize_d = param_d.qsize_d;
  int n0_qdp = param_d.n0_qdp;
  int np1_qdp = param_d.np1_qdp;
  int limiter_option = param_d.limiter_option;
  int rhs_viss = param_d.rhs_viss;

  int slice_qdp = (int)(gl_qdp_leap - gl_qdp);
  //int istep_dp_star = qsize*NLEV*NP*NP;
  int istep_Qtens = qsize*NLEV*NP*NP;
  int istep_qmax = qsize*NLEV;
  int i, j, k, q, ie, l, k_beg, k_end, k_n, q_beg, q_end, q_n, pos;

  double Qdp[KBLK*NP*NP];
  double Qdp_np1[KBLK*NP*NP];
  double dp[KBLK*NP*NP];
  double dp_temp[KBLK*NP*NP];
  double dp_star[KBLK*NP*NP];
  double vn0[KBLK*2*NP*NP];
  double Vstar[KBLK*2*NP*NP];                 // same size as vn0
  double divdp_proj[KBLK*NP*NP];             // same size as dp_tmp
  double gradQ[2*NP*NP];
  double Dvv[NP*NP];
  double Dinv[2*2*NP*NP];
  double metdet[NP*NP];
  double rmetdet[NP*NP];
  double Qtens_temp[KBLK*NP*NP];
  double Qtens_biharmonic[KBLK*NP*NP];
  double dpdiss_biharmonic[KBLK*NP*NP];     // same size as dp
  double divdp[KBLK*NP*NP];                 // same size as dp
  double dpdiss[NP*NP];
  double spheremp[NP*NP];
  double qmax[KBLK];
  double qmin[KBLK];

  int sum_k = qsize*NLEV;

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block + block + block_dp + block_dp_star    \
        + block_vn0 + block_vn0 + block_dp + block_gradQ + NP*NP + 4*NP*NP     \
        + NP*NP + NP*NP + block + block + block_dp + block_dp + NP*NP          \
        + NP*NP + NP*NP + NP*NP + UC*NLEV + UC*NLEV + 6*NP*NP);
    printf("size_tol:%dk\n", size_tol/1024);
  }
#endif

  pe_get(gl_Dvv, Dvv, NP*NP*sizeof(double));
  dma_syn();

  /* local Variables of deivergence_sphere */
  double dudx00, dvdy00;
  double gv[2*NP*NP], vvtemp[NP*NP];
  /* local Variables of deivergence_sphere */
  double cc[NP*NP], xx[NP*NP];
  double addmass, weightssum;
  double mass, sumc;
  double tol_limiter = 5e-14;

  double *src_qdp, *src_qdp_np1, *src_dp, *src_vn0, *src_divdp_proj, *src_Dinv,\
     *src_metdet, *src_rmetdet, *src_Qtens_biharmonic, *src_divdp,             \
     *src_dpdiss_biharmonic, *src_spheremp, *src_qmax, *src_qmin;              \
  double *gl_n0_qdp = gl_qdp + (n0_qdp - 1)*qsize_d*stripe_qdp;
  double *gl_np1_qdp = gl_qdp + (np1_qdp - 1)*qsize_d*stripe_qdp;

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
    //if (id == 62) printf("q_beg:%d, q_end:%d, k_beg:%d, k_end:%d, k_n:%d\n", q_beg, q_end, k_beg, k_end, k_n);
    for (q = 0; q <= q_n; q++) {
      src_dp = gl_dp + ie*slice_qdp + (k_beg + q*(q*NLEV - k_beg) - q*NLEV)*NP*NP;
      src_vn0 = gl_vn0 + ie*slice_qdp + (k_beg + q*(q*NLEV - k_beg) - q*NLEV)*2*NP*NP;
      src_divdp_proj = gl_divdp_proj + ie*slice_qdp                            \
          + (k_beg + q*(q*NLEV - k_beg) - q*NLEV)*NP*NP;
      src_divdp = gl_divdp + ie*slice_qdp + (k_beg + q*(q*NLEV - k_beg) - q*NLEV)*NP*NP;
      src_dpdiss_biharmonic = gl_dpdiss_biharmonic + ie*slice_qdp              \
          + (k_beg + q*(q*NLEV - k_beg) - q*NLEV)*NP*NP;
      src_Dinv = gl_Dinv + ie*slice_qdp;
      src_metdet = gl_metdet + ie*slice_qdp;
      src_rmetdet = gl_rmetdet + ie*slice_qdp;
      src_spheremp = gl_spheremp + ie*slice_qdp;
      src_qdp = gl_n0_qdp + ie*slice_qdp + q_beg*NLEV*NP*NP                    \
          + (k_beg + q*(q*NLEV - k_beg))*NP*NP;
      src_qdp_np1 = gl_np1_qdp + ie*slice_qdp + q_beg*NLEV*NP*NP               \
          + (k_beg + q*(q*NLEV - k_beg))*NP*NP;
      src_Qtens_biharmonic = gl_Qtens_biharmonic + ie*istep_Qtens + q_beg*NLEV*NP*NP \
          + (k_beg + q*(q*NLEV - k_beg))*NP*NP;
      src_qmax = gl_qmax + ie*istep_qmax + q_beg*NLEV                          \
          + (k_beg + q*(q*NLEV - k_beg));
      src_qmin = gl_qmin + ie*istep_qmax + q_beg*NLEV                          \
          + (k_beg + q*(q*NLEV - k_beg));
      pe_get(src_dp, dp, k_n*NP*NP*sizeof(double));
      pe_get(src_vn0, vn0, k_n*2*NP*NP*sizeof(double));
      pe_get(src_divdp_proj, divdp_proj, k_n*NP*NP*sizeof(double));
      pe_get(src_divdp, divdp, k_n*NP*NP*sizeof(double));
      pe_get(src_dpdiss_biharmonic, dpdiss_biharmonic, k_n*NP*NP*sizeof(double));
      pe_get(src_Dinv, Dinv, 2*2*NP*NP*sizeof(double));
      pe_get(src_metdet, metdet, NP*NP*sizeof(double));
      pe_get(src_rmetdet, rmetdet, NP*NP*sizeof(double));
      pe_get(src_spheremp, spheremp, NP*NP*sizeof(double));
      pe_get(src_qdp, Qdp, k_n*NP*NP*sizeof(double));
      pe_get(src_Qtens_biharmonic, Qtens_biharmonic, k_n*NP*NP*sizeof(double));
      pe_get(src_qmax, qmax, k_n*sizeof(double));
      pe_get(src_qmin, qmin, k_n*sizeof(double));
      dma_syn();
      for (k = 0; k < k_n; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos_dp = k*NP*NP + j*NP + i;
            int pos_vn0_1 = k*2*NP*NP + j*NP + i;
            int pos_vn0_2 = k*2*NP*NP + NP*NP + j*NP + i;
            dp_temp[pos_dp] = dp[pos_dp] - rhs_multiplier*dt*divdp_proj[pos_dp];
            Vstar[pos_vn0_1] = vn0[pos_vn0_1]/dp_temp[pos_dp];
            Vstar[pos_vn0_2] = vn0[pos_vn0_2]/dp_temp[pos_dp];
          }
        }
      }
      for (k = 0; k < k_n; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos_qdp = k*NP*NP + j*NP + i;
            int pos_vn0_1 = k*2*NP*NP + j*NP + i;
            int pos_vn0_2 = k*2*NP*NP + NP*NP + j*NP + i;
            int pos_gradQ_1 = j*NP + i;
            int pos_gradQ_2 = NP*NP + j*NP + i;
            int pos_Dinv_1 = j*NP + i;
            int pos_Dinv_2 = NP*NP + j*NP + i;
            int pos_Dinv_3 = 2*NP*NP + j*NP + i;
            int pos_Dinv_4 = 3*NP*NP + j*NP + i;
            int pos_gv_1 = j*NP + i;
            int pos_gv_2 = NP*NP + j*NP + i;
            int pos_det = j*NP + i;
            gradQ[pos_gradQ_1] = Vstar[pos_vn0_1]*Qdp[pos_qdp];
            gradQ[pos_gradQ_2] = Vstar[pos_vn0_2]*Qdp[pos_qdp];
            gv[pos_gv_1] = metdet[pos_det]*(Dinv[pos_Dinv_1]*gradQ[pos_gradQ_1]\
                + Dinv[pos_Dinv_3]*gradQ[pos_gradQ_2]);
            gv[pos_gv_2] = metdet[pos_det]*(Dinv[pos_Dinv_2]*gradQ[pos_gradQ_1]\
                + Dinv[pos_Dinv_4]*gradQ[pos_gradQ_2]);
          }
        }
        for (j = 0; j < NP; j++) {
          for (l = 0; l < NP; l++) {
            dudx00 = 0.0;
            dvdy00 = 0.0;
            for (i = 0; i < NP; i++) {
              int pos_Dvv = l*NP + i;
              int pos_gv_1 = j*NP + i;
              int pos_gv_2 = NP*NP + i*NP + j;
              dudx00 = dudx00 + Dvv[pos_Dvv]*gv[pos_gv_1];
              dvdy00 = dvdy00 + Dvv[pos_Dvv]*gv[pos_gv_2];
            }
            int pos_dp_star = k*NP*NP + j*NP + l;
            int pos_vvtemp = l*NP + j;
            dp_star[pos_dp_star] = dudx00;
            vvtemp[pos_vvtemp] = dvdy00;
          }
        }
        for ( j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos_dp_star = k*NP*NP + j*NP + i;
            int pos_vvtemp = j*NP + i;
            int pos_det = j*NP + i;
            dp_star[pos_dp_star] = (dp_star[pos_dp_star] + vvtemp[pos_vvtemp]) \
                *(rmetdet[pos_det]*rrearth);
          }
        } // end of divergence_sphere
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos = k*NP*NP + j*NP + i;
            int pos_dp_star = k*NP*NP + j*NP + i;
            Qtens_temp[pos] = Qdp[pos] - dt*dp_star[pos_dp_star];
          }
        }
        if (rhs_viss != 0) {
          for (j = 0; j < NP; j++) {
            for (i = 0; i < NP; i++) {
              int pos = k*NP*NP + j*NP + i;
              Qtens_temp[pos] = Qtens_temp[pos] + Qtens_biharmonic[pos];
            }
          }
        }
      } // end loop k
      if (limiter_option == 8) {
        for (k = 0; k < k_n; k++) {
          for (j = 0; j < NP; j++) {
            for (i = 0; i < NP; i++) {
              int pos = k*NP*NP + j*NP + i;
              dp_star[pos] = dp_temp[pos] - dt*divdp[pos];
            }
          }
          if (nu_p > 0 && rhs_viss != 0) {
            for (j = 0; j < NP; j++) {
              for (i = 0; i < NP; i++) {
                int pos_dpdiss = j*NP + i;
                int pos_spheremp = j*NP + i;
                int pos_dpdiss_bi = k*NP*NP + j*NP + i;
                int pos_dp_star = k*NP*NP + j*NP + i;
                dpdiss[pos_dpdiss] = dpdiss_biharmonic[pos_dpdiss_bi];
                dp_star[pos_dp_star] = dp_star[pos_dp_star]                    \
                    - rhs_viss*dt*nu_q*dpdiss[pos_dpdiss]/spheremp[pos_spheremp];
              }
            }
          }
        }
        limiter_optim_iter_full_f_(Qtens_temp, spheremp, qmin, qmax, dp_star, &k_n);
      } // end if limiter_option = 8
      for ( k = 0; k < k_n; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos = k*NP*NP + j*NP + i;
            int pos_spheremp = j*NP + i;
            Qdp_np1[pos] = spheremp[pos_spheremp]*Qtens_temp[pos];
          }
        }
      } // end loop k
      pe_put(src_qdp_np1, Qdp_np1, k_n*NP*NP*sizeof(double));
      pe_put(src_qmax, qmax, k_n*sizeof(double));
      pe_put(src_qmin, qmin, k_n*sizeof(double));
      dma_syn();
    } // end loop q
  } // end loop ie
#if 0
  if (id == 0)
    printf("value:%lf\n", (tol_limiter*5E-3*5.82989072416270755E-290));
#endif
}
