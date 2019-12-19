#include <athread.h>

#define NLAY 31
#define NGPTSW 224
#define MXMOL 38
#define NBNDSW 14

struct param_t{
  int lchnk, iplon, nlayers, istart, iend; 
  double *palbd, *palbp;
  double *pcldfmc, *ptaucmc, *pasycmc, *pomgcmc;
  double *ptaua, *pasya, *pomga;
  double prmu0;
  double *adjflux;

  double *pbbfd, *pbbfu, *pbbcd, *pbbcu;
  double *puvfd, *puvcd, *pnifd, *pnicd;
  double *pnifu, *pnicu;
  double *pbbfddir, *pbbcddir, *puvfddir, *puvcddir;
  double *pnifddir, *pnicddir;
  double *pbbfsu, *pbbfsd;
  double *zsflxzen, *ztaug, *ztaur;
  int *ngc;
  double od_lo, tblint, bpade, *exp_tbl;
};
//extern SLAVE_FUN(spcvmc_sw_lyx_trns_slv)();
//extern SLAVE_FUN(spcvmc_sw_lyx_trns_slv)(struct param_t *gl_pm);
//extern void slave_spcvmc_sw_lyx_trns_slv(struct param_t *gl_pm);
void spcvmc_sw_parallel_(int *lchnk, int *iplon, int *nlayers, int *istart, int *iend, int *icpr,
                         int *idelm, int *iout,
                         double *pavel, double *tavel, double *pz, double *tz,
                         double *tbound, double *palbd, double *palbp,
                         double (*pcldfmc)[NLAY], double (*ptaucmc)[NLAY], double (*pasycmc)[NLAY],
                         double (*pomgcmc)[NLAY], double (*ptaormc)[NLAY],
                         double (*ptaua)[NLAY], double (*pasya)[NLAY], double (*pomga)[NLAY],
                         double *prmu0,
                         double *coldry, double (*wkl)[MXMOL], double *adjflux,
                         int *laytrop, int *layswtch, int *laylow,
                         int *jp, int *jt, int *jt1,
                         double *co2mult, double *colch4, double *colco2, double *colh2o,
                         double *colmol, double *coln2o, double *colo2, double *colo3,
                         double *fac00, double *fac01, double *fac10, double *fac11,
                         double *selffac, double *selffrac, int *indself,
                         double *forfac, double *forfrac, int *indfor,
                         double *pbbfd, double *pbbfu, double *pbbcd, double *pbbcu,
                         double *puvfd, double *puvcd, double *pnifd, double *pnicd,
                         double *pnifu, double *pnicu,
                         double *pbbfddir, double *pbbcddir, double *puvfddir, double *puvcddir,
                         double *pnifddir, double *pnicddir,
                         double (*pbbfsu)[NBNDSW], double (*pbbfsd)[NBNDSW],
                         double *zsflxzen, double (*ztaug)[NLAY], double (*ztaur)[NLAY],
                         int *ngc, int *ngs,
                         double *od_lo, double *tblint, double *bpade, double *exp_tbl){

  /*
    output :
    pbbcd, pbbcu, pbbfd, pbbfu, pbbfddir, pbbcddir, puvcd, puvfd, puvcddir, puvfddir,
    pnicd, pnifd, pnicddir, pnifddir, pnicu, pnifu, pbbfsu, pbbfsd
  */
  struct param_t pm;
  extern void slave_spcvmc_sw_lyx_trns_slv(struct param_t *);

  pm.lchnk = *lchnk, pm.iplon = *iplon, pm.nlayers = *nlayers, pm.istart = *istart, pm.iend = *iend; 
  pm.pcldfmc = (double*)pcldfmc, pm.ptaucmc = (double*)ptaucmc, pm.pasycmc = (double*)pasycmc, pm.pomgcmc = (double*)pomgcmc;
  pm.palbd = palbd, pm.palbp = palbp;
  pm.ptaua = (double*)ptaua, pm.pasya = (double*)pasya, pm.pomga = (double*)pomga;
  pm.prmu0 = *prmu0;
  pm.adjflux = adjflux;

  pm.pbbfd = pbbfd, pm.pbbfu = pbbfu, pm.pbbcd = pbbcd, pm.pbbcu = pbbcu;
  pm.puvfd = puvfd, pm.puvcd = puvcd, pm.pnifd = pnifd, pm.pnicd = pnicd;
  pm.pnifu = pnifu, pm.pnicu = pnicu;
  pm.pbbfddir = pbbfddir, pm.pbbcddir = pbbcddir, pm.puvfddir = puvfddir, pm.puvcddir = puvcddir;
  pm.pnifddir = pnifddir, pm.pnicddir = pnicddir;
  pm.pbbfsu = (double*)pbbfsu, pm.pbbfsd = (double*)pbbfsd;
  pm.zsflxzen = zsflxzen, pm.ztaug = (double*)ztaug, pm.ztaur = (double*)ztaur;
  pm.ngc = ngc;
  pm.od_lo = *od_lo, pm.tblint = *tblint, pm.bpade = *bpade, pm.exp_tbl = exp_tbl;
  
  //extern void slave_spcvmc_sw_lyx_trns_slv();
  athread_spawn(spcvmc_sw_lyx_trns_slv, &pm);
  athread_join();

}
