#include <slave.h>
#include "dma_macros.h"
//#include "ldm_alloc.h"

#define NP    4
#define NLEV  30
#define UC    1         // a unit that divides nete by column direction
#define UR    2         // a unit that divides qsize by row direcion
#define NC    8
#define NR    8
#define maxiter           (NP*NP - 1)
#define stripe_qdp        (NLEV*NP*NP)
#define istep_dp          (NLEV*NP*NP)
#define qstep_Qten        (NLEV*NP*NP)       // stripe in q axis of Qtens_biharmonic array
#define qstep_dp_star     (NLEV*NP*NP)
#define qstep_qmax        NLEV
#define qstep_Qtens       (NLEV*NP*NP)
#define block             (UC*NLEV*NP*NP)
#define block_dp          (NLEV*NP*NP)
#define block_vn0         (NLEV*2*NP*NP)
#define block_gradQ       (2*NP*NP)
#define block_dp_star     (NLEV*NP*NP)
#define block_Qtens       (NLEV*NP*NP)
#define block_Dinv        (2*2*NP*NP)
#define block_det         (NP*NP)

#define abs(value, ret) asm volatile ("fcpys $31, %1, %0" : "=r"(ret) : "r"(value))

#define sum_array(sum, array) {   \
  int i;   \
  sum = 0.0;   \
  for (i = 0; i < 16; i++)     \
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
  sum = 0.0;     \
  for (i = 0; i < 16; i++) {  \
    _array[i] = a[i] * b[i]; \
    sum = sum + _array[i];  \
  } \
}

typedef struct {
  double *spheremp, *qmax, *qmin, *Qtens_temp, *Qtens_temp_test, *dp_star_temp;
  int nets, nete, qsize, limiter_option, step_elem;
} param_t;

extern void limiter_optim_iter_full_f_(double *Qtens_temp, double *spheremp    \
    , double *qmin, double *qmax, double *dp_star);

void slave_euler_limit_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
	dma_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();

  double *gl_spheremp = param_d.spheremp;
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
  double *gl_Qtens_temp = param_d.Qtens_temp;
  double *gl_Qtens_temp_test = param_d.Qtens_temp_test;
  double *gl_dp_star_temp = param_d.dp_star_temp;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;
  int limiter_option = param_d.limiter_option;
  int step_elem = param_d.step_elem;

  rid = id / NC;
  cid = id % NC;
  int loop_r = ((nete - nets + 1) + UR*NR - 1)/(UR*NR);
  int loop_c = (qsize + UC*NC - 1)/(UC*NC);
  int istep_dp_star = qsize*NLEV*NP*NP;
  int istep_Qtens = qsize*NLEV*NP*NP;
  int istep_qmax = qsize*NLEV;
  int c, r, i, j, l, k, k1, iter, q, ie, cbeg, cend, rbeg, rend, cn, rn;

  double spheremp[NP*NP];
  double qmax[UC*NLEV];
  double qmin[UC*NLEV];
  double dp_star[UC*NLEV*NP*NP];
  double Qtens_temp[UC*NLEV*NP*NP];

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block + block + block_dp + block_dp_star    \
        + block_vn0 + block_vn0 + block_dp + block_gradQ + NP*NP + 4*NP*NP     \
        + NP*NP + NP*NP + block + block + block_dp + block_dp + NP*NP          \
        + NP*NP + NP*NP + NP*NP + UC*NLEV + UC*NLEV + 6*NP*NP);
    printf("size_tol:%dk\n", size_tol/1024);
  }
#endif

  /* local Variables of deivergence_sphere */
  double dudx00, dvdy00;
  double gv[2*NP*NP], vvtemp[NP*NP];
  /* local Variables of deivergence_sphere */
  double cc[NP*NP], xx[NP*NP];
  double addmass, weightssum;
  double mass, sumc;
  double tol_limiter = 5e-14;

  double *src_spheremp, *src_qmax, *src_qmin;
  // for Debug
  double *src_Qtens_temp, *src_dp_star_temp, *src_Qtens_temp_test;
  // Divide ie-axis data on the row cpe with loop_r
  // the value of loop_r rely on NR, UR; NP is the number of cloumn cpe,
  // UR is the unit that divides ie size by cloumn direcion,
  // Divide q-axis data on the cloumn cpe with loop_c
  // same as loop_r, the loop_c rely on NC, UC
  // UC is the unit that divides q size by row direction
  for (r = 0; r < loop_r; r++) {
    rbeg = r*NR*UR + rid*UR;
    rend = r*NR*UR + (rid + 1)*UR;
    rend = rend < (nete - nets + 1) ? rend : (nete - nets + 1);
    rn = rend - rbeg;
    rn = rn < 0 ? 0 : rn;   // handling boundary issues, removing the case where rn < 0
    for (ie = 0; ie < rn; ie++) {
      src_spheremp = gl_spheremp + (rbeg + ie)*step_elem;
      pe_get(src_spheremp, spheremp, NP*NP*sizeof(double));
      dma_syn();
      for (c = 0; c < loop_c; c++) {
        cbeg = c*NC*UC + cid*UC;
        cend = c*NC*UC + (cid + 1)*UC;
        cend = cend < qsize ? cend : qsize;
        cn = cend - cbeg;
        if (cn > 0) {   // if cn < 0, the dma will get exceptional contribution
          src_qmax = gl_qmax + (rbeg + ie)*istep_qmax + cbeg*qstep_qmax;
          src_qmin = gl_qmin + (rbeg + ie)*istep_qmax + cbeg*qstep_qmax;
          src_Qtens_temp = gl_Qtens_temp + (rbeg + ie)*qsize*NLEV*NP*NP + cbeg*NLEV*NP*NP;
          src_Qtens_temp_test = gl_Qtens_temp_test + (rbeg + ie)*qsize*NLEV*NP*NP + cbeg*NLEV*NP*NP;
          src_dp_star_temp = gl_dp_star_temp + (rbeg + ie)*qsize*NLEV*NP*NP + cbeg*NLEV*NP*NP;
          pe_get(src_qmax, qmax, UC*NLEV*sizeof(double));
          pe_get(src_qmin, qmin, UC*NLEV*sizeof(double));
          pe_get(src_Qtens_temp, Qtens_temp, cn*NLEV*NP*NP*sizeof(double));
          pe_get(src_dp_star_temp, dp_star, cn*NLEV*NP*NP*sizeof(double));
          dma_syn();
          for (q = 0; q < cn; q++) {
            if (limiter_option == 8) {
              limiter_optim_iter_full_f_(Qtens_temp, spheremp    \
                  , &qmin[q*NLEV], &qmax[q*NLEV], dp_star);
            }   // end if limiter_option = 8
            //for (k = 0; k < NLEV; k++) {
            //  for (j = 0; j < NP; j++) {
            //    for (i = 0; i < NP; i++) {
            //      int pos_Qtens = k*NP*NP + j*NP + i;
            //      int pos = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
            //      int pos_spheremp = j*NP + i;
            //      Qdp_np1[pos] = spheremp[pos_spheremp]*Qtens_temp[pos_Qtens];
            //    }
            //  }
            //}
          }     // end loop q
          //pe_put(src_qmax, qmax, cn*NLEV*sizeof(double));
          //pe_put(src_qmin, qmin, cn*NLEV*sizeof(double));
          //pe_put(src_Qtens_temp, Qtens_temp, cn*NLEV*NP*NP*sizeof(double));
          //pe_put(src_dp_star_temp, dp_star, cn*NLEV*NP*NP*sizeof(double));
          pe_put(src_Qtens_temp_test, Qtens_temp, cn*NLEV*NP*NP*sizeof(double));
          dma_syn();
        }
      }
    }
  }
}
