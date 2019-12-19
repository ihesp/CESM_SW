#include <slave.h>
#include <dma.h>
#include <math.h>

#define get_row_id(row) asm volatile("rcsr %0, 1" : "=r"(row));
#define get_col_id(col) asm volatile("rcsr %0, 2" : "=r"(col));
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define memb() asm volatile("memb")
#define my_dma_wait(reply, count) while ((*(reply)) != count) memb()

#define dma_set_size_wb(desc, size) {			\
    dma_set_size(desc, size);				\
    dma_set_bsize(desc, size);				\
  }

#define dma_init_get(desc) {				\
    memset(desc, 0, sizeof(dma_desc));			\
    dma_set_stepsize(desc, 0);				\
    dma_set_op(desc, DMA_GET);				\
    dma_set_mode(desc, PE_MODE);			\
  }

#define dma_init_put(desc) {				\
    memset(desc, 0, sizeof(dma_desc));			\
    dma_set_stepsize(desc, 0);				\
    dma_set_op(desc, DMA_PUT);				\
    dma_set_mode(desc, PE_MODE);			\
  }

#define dma_supp(x, y, z) {			\
    asm volatile("memb");			\
    dma(x, (long long)(y), (long long)(z));	\
    asm volatile("memb");			\
  }

#define min(x, y) ((x) < (y) ? (x) : (y))
#define dma_get(m_ptr, s_ptr) {                    \
  get_reply = 0;                                   \
  dma_supp(get_desc, m_ptr, s_ptr);                \
  my_dma_wait(&get_reply, 1);                      \
}
#define dma_put(m_ptr, s_ptr) {                    \
  put_reply = 0;                                   \
  dma_supp(put_desc, m_ptr, s_ptr);                \
  my_dma_wait(&put_reply, 1);                      \
}

void reftra_sw_dxh_trns_(dma_desc *desc, volatile int *reply, int *nlay, int *lrtchk, 
                         double *pgg, double *prmuz, double *ptau, double *pw,
                         double *pref, double *prefd, double *ptra, double *ptrad,
                         double *od_lo, double *tblint, double *bpade, double *exp_tbl){
  double zsr3 = 1.7320508075688772; //sqrt(3)
  double zwcrit = 0.9999995;
  double eps = 1e-8;
  int jk;
  int klev = *nlay;

  //printf("%8x %8x\n", lrtchk[0], lrtchk[1]);
  for (jk = 0; jk < klev; jk ++){
    int ikl = klev - jk - 1;
    if (!lrtchk[ikl]){
      pref[jk] = 0.0;
      ptra[jk] = 1.0;
      prefd[jk] = 0.0;
      ptrad[jk] = 1.0;
    } else {
      // use diff reverse from pref
      double zto1 = ptau[ikl];
      double zw = pw[ikl];
      double zg = pgg[ikl];
      double zg3 = 3.0 * zg;

      double zgamma1, zgamma2, zgamma3, zgamma4;
      zgamma1 = (8.0 - zw * (5.0 + zg3)) * 0.25;
      zgamma2 = 3.0 * (zw * (1.0 - zg)) * 0.25;
      zgamma3 = (2.0 - zg3 * (*prmuz)) * 0.25;
      zgamma4 = 1.0 - zgamma3;
      double zwo = zw / (1.0 - (1.0 - zw) * (zg / (1.0 - zg)) * (zg / (1.0 - zg)));
      //double za, za1, zgt, ze1, ze2;
      //if (isnan(zwo))
      //printf("%f %f\n", zw, zg);
      if (zwo >= zwcrit){
        double za = zgamma1 * (*prmuz);
        double za1 = za - zgamma3;
        double zgt = zgamma1 * zto1;
        double ze1 = min(zto1 / (*prmuz), 500.0);
        double ze2;
        if (ze1 <= (*od_lo))
          ze2 = 1.0 - ze1 + 0.5 * ze1 * ze1;
        else{
          double tblind = ze1 / ((*bpade) + ze1);
          int itind = (*tblint) * tblind + 0.5;
          //double exp_tbl_itind; 
          //dma_set_size_wb(desc, sizeof(double));
          //*reply = 0;
          //dma_supp(*desc, exp_tbl + itind, &exp_tbl_itind);
          //my_dma_wait(reply, 1);
          //ze2 = exp_tbl_itind;
          ze2 = exp_tbl[itind];
        }
        //ze2 = exp(-ze1);
        pref[jk] = (zgt - za1 * (1.0 - ze2)) / (1.0 + zgt);
        ptra[jk] = 1.0 - pref[jk];

        prefd[jk] = zgt / (1.0 + zgt);
        ptrad[jk] = 1.0 - prefd[jk];

        if (ze2 == 1.0){
          pref[jk] = 0.0;
          ptra[jk] = 1.0;
          prefd[jk] = 0.0;
          ptrad[jk] = 1.0;
        }
      } else {
        double za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3;
        double za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4;
        double zrk = sqrt(zgamma1 * zgamma1 - zgamma2 * zgamma2);
        double zrp = zrk * (*prmuz);
        double zrp1 = 1.0 + zrp;
        double zrm1 = 1.0 - zrp;
        double zrk2 = 2.0 * zrk;
        double zrpp = 1.0 - zrp * zrp;
        double zrkg = zrk + zgamma1;
        double zr1 = zrm1 * (za2 + zrk * zgamma3);
        double zr2 = zrp1 * (za2 - zrk * zgamma3);
        double zr3 = zrk2 * (zgamma3 - za2 * (*prmuz));
        double zr4 = zrpp * zrkg;
        double zr5 = zrpp * (zrk - zgamma1);
        double zt1 = zrp1 * (za1 + zrk * zgamma4);
        double zt2 = zrm1 * (za1 - zrk * zgamma4);
        double zt3 = zrk2 * (zgamma4 + za1 * (*prmuz));
        double zt4 = zr4;
        double zt5 = zr5;
        double zbeta = (zgamma1 - zrk) / zrkg;

        double zem1, zem2, zep1, zep2;
        double ze1 = min(zrk * zto1, 500.0);
        double ze2 = min(zto1 / (*prmuz), 500.0);
        if (ze1 <= (*od_lo))
          zem1 = 1.0 - ze1 + 0.5 * ze1 * ze1;
        else{
          double tblind = ze1 / ((*bpade) + ze1);
          int itind = (*tblint) * tblind + 0.5;
          //double exp_tbl_itind; 
          //dma_set_size_wb(desc, sizeof(double));
          //*reply = 0;
          //dma_supp(*desc, exp_tbl + itind, &exp_tbl_itind);
          //my_dma_wait(reply, 1);
          //zem1 = exp_tbl_itind;
          zem1 = exp_tbl[itind];
        }
        //zem1 = exp(-ze1);
        zep1 = 1.0 / zem1;
        if (ze2 <= (*od_lo))
          zem2 = 1.0 - ze2 + 0.5 * ze2 * ze2;
        else {
          double tblind = ze2 / ((*bpade) + ze2);
          int itind = (*tblint) * tblind + 0.5;
          //double exp_tbl_itind; 
          //dma_set_size_wb(desc, sizeof(double));
          //*reply = 0;
          //dma_supp(*desc, exp_tbl + itind, &exp_tbl_itind);
          //my_dma_wait(reply, 1);
          //zem2 = exp_tbl_itind;
          zem2 = exp_tbl[itind];
        }
        //zem2 = exp(-ze2);
        zep2 = 1.0 / zem2;

        double zdenr = zr4 * zep1 + zr5 * zem1;
        double zdent = zt4 * zep1 + zt5 * zem1;

        if (zdenr >= -eps && zdenr <= eps){
          pref[jk] = eps;
          ptra[jk] = zem2;
        } else {
          pref[jk] = zw * (zr1 * zep1 - zr2 * zem1 - zr3 * zem2) / zdenr;
          ptra[jk] = zem2 - zem2 * zw * (zt1 * zep1 - zt2 * zem1 - zt3 * zep2) / zdent;
        }
        double zemm = zem1 * zem1;
        double zdend = 1.0 / ((1.0 - zbeta * zemm) * zrkg);
        prefd[jk] = zgamma2 * (1.0 - zemm) * zdend;
        ptrad[jk] = zrk2 * zem1 * zdend;
      }
    }
  }
}

void vrtqdr_sw_c_(int *klev_p, double *pref, double *prefd, double *ptra, double *ptrad, 
               double *pdbt, double *prdnd, double *prup, double *prupd, double *ptdbt, 
               double *pfd, double *pfu){
  //integer, parameter :: ngptsw = 112
  //klev is a variable, ne30 = 31

  int klev = *klev_p;
  int jk;
  int klevm1 = klev - 1;
  double ztdn[32];

  double zreflect = 1. / (1. - prefd[klev] * prefd[klevm1]);
  prup[1] = pref[klevm1] + (ptrad[klevm1] * ((ptra[klevm1] - pdbt[1]) * prefd[klev  ] + pdbt[1] * pref[klev  ])) * zreflect;
  prupd[1] = prefd[klevm1] + ptrad[klevm1] * ptrad[klevm1] * prefd[klev  ] * zreflect;

  // can't vectorization
  for(jk = 0; jk < klev - 1; ++ jk){
    // jk : 0 ~ 29 ikp : 30 ~ 1 ikx : 29 ~ 0
    int ikp = klev - 1 - jk;
    int ikx = ikp - 1;
    int jk1 = jk + 1, jk2 = jk + 2;
    zreflect = 1. / (1. - prupd[jk1] * prefd[ikx]);
    prup[jk2] = pref[ikx] + (ptrad[ikx] * ((ptra[ikx] - pdbt[klev - ikx]) * prupd[jk1] + pdbt[jk2] * prup[jk1])) * zreflect;
    prupd[jk2] = prefd[ikx] + ptrad[ikx] * ptrad[ikx] * prupd[jk1] * zreflect;
  }

  ztdn[31] = 1.;
  prdnd[31] = 0.;
  ztdn[30] = ptra[0];
  prdnd[30] = prefd[0];

  // can't vectorization
  for(jk = 1; jk < klev; ++ jk){
    // jk : 1 ~ 30 ikp : 2 ~ 31
    int ikl1 = klev - jk - 1;
    int ikl = klev - jk;
    zreflect = 1. / (1. - prefd[jk] * prdnd[ikl]);
    ztdn[ikl1] = ptdbt[ikl] * ptra[jk] + (ptrad[jk] * ((ztdn[ikl] - ptdbt[ikl]) + ptdbt[ikl] * pref[jk] * prdnd[ikl])) * zreflect;
    prdnd[ikl1] = prefd[jk] + ptrad[jk] * ptrad[jk] * prdnd[ikl] * zreflect;
  }

  // vectorization
  for(jk = 0; jk < 32; ++ jk){
    // jk : 0 ~ 31
    zreflect = 1. / (1. - prdnd[jk] * prupd[jk]);
    pfu[jk] = (ptdbt[jk] * prup[jk] + (ztdn[jk] - ptdbt[jk]) * prupd[jk]) * zreflect;
    pfd[jk] = ptdbt[jk] + (ztdn[jk] - ptdbt[jk] + ptdbt[jk] * prup[jk] * prdnd[jk]) * zreflect;
  }
}

#define NLAY 31
#define NGPTSW 224
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

void slave_spcvmc_sw_lyx_trns_slv(struct param_t *gl_pm){

  /*
    output :
    pbbcd, pbbcu, pbbfd, pbbfu, pbbfddir, pbbcddir, puvcd, puvfd, puvcddir, puvfddir,
    pnicd, pnifd, pnicddir, pnifddir, pnicu, pnifu, pbbfsu, pbbfsd
  */
  int i, j;

  double zclear, zcloud;
  double zdbt[NLAY + 1], zdbtc[NLAY + 1];
  double zgc[NLAY], zgcc[NLAY], zgco[NLAY];
  double zomc[NLAY], zomcc[NLAY], zomco[NLAY];
  double zrdnd[NLAY + 1], zrdndc[NLAY + 1];
  double zref[NLAY + 1], zrefc[NLAY + 1], zrefo[NLAY + 1];
  double zrefd[NLAY + 1], zrefdc[NLAY + 1], zrefdo[NLAY + 1];
  double zrup[NLAY + 1], zrupd[NLAY + 1];
  double zrupc[NLAY + 1], zrupdc[NLAY + 1];
  double zs1[NLAY + 1];
  double ztauc[NLAY], ztauo[NLAY];
  double ztdn[NLAY + 1], ztdnd[NLAY + 1];
  double ztoc[NLAY], ztor[NLAY];
  double ztra[NLAY + 1], ztrac[NLAY + 1], ztrao[NLAY + 1];
  double ztrad[NLAY + 1], ztradc[NLAY + 1], ztrado[NLAY + 1];
  doublev4 ztdbt[8], ztdbtc[8];
  doublev4 zcd[8], zcu[8];
  doublev4 zfd[8], zfu[8];
  int lrtchkclr[NLAY], lrtchkcld[NLAY];

  int iw2ibm[4];   // NGPTSW / 64 : num of iw for every slave
  int iw2jb[4];

  // ============================ param ================================
  // === in ===
  // every iw
  double pcldfmc[32], ptaucmc[32], pasycmc[32], pomgcmc[32];
  double ztaug[32], ztaur[32];
  // every ibm
  double ptaua[32], pasya[32], pomga[32];
  // once
  double adjflux[29];       // jbband == 29
  double zsflxzen[NGPTSW];  // TODO NGPTSW to 4
  double palbd[NBNDSW], palbp[NBNDSW];
  int ngc[NBNDSW];
  // exp_tbl host

  // === out ===
  doublev4 pbbfd[8], pbbfu[8], pbbcd[8], pbbcu[8];
  doublev4 puvfd[8], puvcd[8], pnifd[8], pnicd[8];
  doublev4 pnifu[8], pnicu[8];
  doublev4 pbbfddir[8], pbbcddir[8], puvfddir[8], puvcddir[8];
  doublev4 pnifddir[8], pnicddir[8];
  doublev4 pbbfsu_T[NBNDSW][8], pbbfsd_T[NBNDSW][8];
  doublev4 pbbfsu[8 * NBNDSW], pbbfsd[8 * NBNDSW]; // vectorization for 2d

  dma_desc get_desc, put_desc;
  volatile int get_reply, put_reply;

  struct param_t pm;

  int rid, cid ,sid;
  get_row_id(rid);
  get_col_id(cid);
  sid = (rid << 3) + cid;

  dma_init_get(&get_desc);
  dma_init_put(&put_desc);

  dma_set_reply(&get_desc, &get_reply);
  dma_set_reply(&put_desc, &put_reply);

  dma_set_size_wb(&get_desc, sizeof(struct param_t));
  dma_get(gl_pm, &pm);

  dma_set_size_wb(&get_desc, sizeof(double) * 29);
  dma_get(pm.adjflux, adjflux);
  dma_set_size_wb(&get_desc, sizeof(double) * NGPTSW);
  dma_get(pm.zsflxzen, zsflxzen);
  dma_set_size_wb(&get_desc, sizeof(double) * NBNDSW);
  dma_get(pm.palbd, palbd);
  dma_get(pm.palbp, palbp);
  dma_set_size_wb(&get_desc, sizeof(int) * NBNDSW);
  dma_get(pm.ngc, ngc);

  zclear = zcloud = 0;

  int ib1 = pm.istart;
  int ib2 = pm.iend;
  int klev = pm.nlayers;  // klev == 31

  double repclc = 1e-12;
  int jb, jg, jk;
  doublev4 v4_0 = 0.0;

  for (i = 0; i < 8; ++ i){
    pbbcd   [i] = v4_0;
    pbbcu   [i] = v4_0;
    pbbfd   [i] = v4_0;
    pbbfu   [i] = v4_0;
    pbbcddir[i] = v4_0;
    pbbfddir[i] = v4_0;
    puvcd   [i] = v4_0;
    puvfd   [i] = v4_0;
    puvcddir[i] = v4_0;
    puvfddir[i] = v4_0;
    pnicd   [i] = v4_0;
    pnifd   [i] = v4_0;
    pnicddir[i] = v4_0;
    pnifddir[i] = v4_0;
    pnicu   [i] = v4_0;
    pnifu   [i] = v4_0;
  }
  for (i = 0; i < NBNDSW * 8; ++ i){
    pbbfsu[i] = v4_0;
    pbbfsd[i] = v4_0;
  }

  int iw, iwn;
  int iwb = 4 * ((sid + 1) >> 1) + 3 * (sid >> 1);
  int iwl = 3 + (~sid & 1);
  int iwe = iwb + iwl;
  int ibm_last = -1;

  for (iw = 0, jb = ib1 - 1; jb < ib2; ++ jb){
    int ibm = jb - 15;
    int igt = ngc[ibm];
    // delete code for iout > 0
    for(jg = 0; jg < igt; ++ jg, ++ iw){
      if(iwb <= iw && iw < iwe){
        iw2ibm[iw - iwb] = ibm;
        iw2jb [iw - iwb] = jb;
        if(ibm_last != ibm){
          ibm_last = ibm;
          for (i = 0; i < 8; ++ i){
            pbbfsu_T[ibm][i] = v4_0;
            pbbfsd_T[ibm][i] = v4_0;
          }
        }
      }
    }
  }
  iwn = iw;               //  iwn == 224
  iwe = (iwe < iwn? iwe : iwn);


  ibm_last = -1;
  for(iw = iwb; iw < iwe; ++ iw){
    int jb = iw2jb[iw - iwb];
    int ibm = iw2ibm[iw - iwb];
    double zincflx = adjflux[jb] * zsflxzen[iw] * pm.prmu0;
    
    dma_set_size_wb(&get_desc, sizeof(double) * NLAY);
    dma_get(pm.pcldfmc + iw * NLAY, pcldfmc); 
    dma_get(pm.ptaucmc + iw * NLAY, ptaucmc); 
    dma_get(pm.pasycmc + iw * NLAY, pasycmc); 
    dma_get(pm.pomgcmc + iw * NLAY, pomgcmc); 
    dma_get(pm.ztaug + iw * NLAY, ztaug); 
    dma_get(pm.ztaur + iw * NLAY, ztaur); 

    if(ibm != ibm_last){
      ibm_last = ibm;
      dma_get(pm.ptaua + ibm * NLAY, ptaua); 
      dma_get(pm.pasya + ibm * NLAY, pasya); 
      dma_get(pm.pomga + ibm * NLAY, pomga); 
    }
  // extract and no check
  {
    ((double*)ztdbtc)[klev] = 1.0;
    zdbtc[0]  = 0.0;
    ztrac[klev]  = 0.0;
    ztradc[klev] = 0.0;
  
    ztrao[klev]  = 0.0;
    ztrado[klev] = 0.0;

    ((double*)ztdbt)[klev] = 1.0;
    zdbt[0]  = 0.0;
    ztra[klev]  = 0.0;
    ztrad[klev] = 0.0;
  }
    
    zrefc[klev]  = palbp[ibm];
    zrefdc[klev] = palbd[ibm];
    zrupc[0]  = palbp[ibm];
    zrupdc[0] = palbd[ibm];

    zrefo[klev]  = palbp[ibm];
    zrefdo[klev] = palbd[ibm];
    
    zref[klev]  = palbp[ibm];
    zrefd[klev] = palbd[ibm];
    zrup[0]  = palbp[ibm];
    zrupd[0] = palbd[ibm];

    int jk;
    // vectorization
    for (jk = 0; jk < klev; jk ++){
      // jk : 0 ~ 30 
      lrtchkclr[jk] = 1;
      lrtchkcld[jk] = (pcldfmc[jk] > repclc);
      
      ztauc[jk] = ztaur[jk] + ztaug[jk] + ptaua[jk];
      zomcc[jk] = ztaur[jk] + ptaua[jk] * pomga[jk];   // * 1.
      zgcc[jk]  = pasya[jk] * pomga[jk] * ptaua[jk] / zomcc[jk];
      zomcc[jk] = zomcc[jk] / ztauc[jk];

      double zf = zgcc[jk] * zgcc[jk];
      double zwf = zomcc[jk] * zf;
      ztauc[jk] = (1.0 - zwf) * ztauc[jk];
      zomcc[jk] = (zomcc[jk] - zwf) / (1.0 - zwf);
      zgcc[jk] = (zgcc[jk] - zf) / (1.0 - zf);

      ztauo[jk] = ztauc[jk] + ptaucmc[jk];
      zomco[jk] = ztauc[jk] * zomcc[jk] + ptaucmc[jk] * pomgcmc[jk];
      zgco[jk] = (ptaucmc[jk] * pomgcmc[jk] * pasycmc[jk] +
                  ztauc[jk] * zomcc[jk] * zgcc[jk]) / zomco[jk];
      zomco[jk] = zomco[jk] / ztauo[jk];
    }

    reftra_sw_dxh_trns_(&get_desc, &get_reply, &klev, lrtchkclr, 
                        zgcc, &pm.prmu0, ztauc, zomcc, 
                        zrefc, zrefdc, ztrac, ztradc,
                        &pm.od_lo, &pm.tblint, &pm.bpade, pm.exp_tbl);
    reftra_sw_dxh_trns_(&get_desc, &get_reply, &klev, lrtchkcld, 
                        zgco, &pm.prmu0, ztauo, zomco,
                        zrefo, zrefdo, ztrao, ztrado,
                        &pm.od_lo, &pm.tblint, &pm.bpade, pm.exp_tbl);

    // vectorization
    for (jk = 0; jk < klev; jk ++){
      // jk : 0 ~ 30 ikl : 30 ~ 0
      int ikl = klev - jk - 1;
      double zclear = 1.0 - pcldfmc[ikl];
      double zcloud = pcldfmc[ikl];

      zref [jk] = zclear * zrefc [jk] + zcloud *  zrefo[jk];
      zrefd[jk] = zclear * zrefdc[jk] + zcloud * zrefdo[jk];
      ztra [jk] = zclear * ztrac [jk] + zcloud *  ztrao[jk];
      ztrad[jk] = zclear * ztradc[jk] + zcloud * ztrado[jk];

      double ze1 = ztauc[ikl] / pm.prmu0;
      double zdbtmc;
      if (ze1 <= pm.od_lo)
        zdbtmc = 1.0 - ze1 + 0.5 * ze1 * ze1;
      else {
        double tblind = ze1 / (pm.bpade + ze1);
        int itind = pm.tblint * tblind + 0.5;
        //double exp_tbl_itind; 
        //dma_set_size_wb(&get_desc, sizeof(double));
        //dma_get(pm.exp_tbl + itind, &exp_tbl_itind);
        //zdbtmc = exp_tbl_itind;
        zdbtmc = pm.exp_tbl[itind];
      }

      zdbtc[ikl + 1] = zdbtmc;
      ((double*)ztdbtc)[ikl] = zdbtc[ikl + 1] * ((double*)ztdbtc)[ikl + 1];

      ze1 = ztauo[ikl] / pm.prmu0;
      double zdbtmo;
      if (ze1 <= pm.od_lo)
        zdbtmo = 1.0 - ze1 + 0.5 * ze1 * ze1;
      else{
        // can't vectorization
        double tblind = ze1 / (pm.bpade + ze1);
        int itind = pm.tblint * tblind + 0.5;
        //double exp_tbl_itind; 
        //dma_set_size_wb(&get_desc, sizeof(double));
        //dma_get(pm.exp_tbl + itind, &exp_tbl_itind);
        //zdbtmo = exp_tbl_itind;
        zdbtmo = pm.exp_tbl[itind];
      }
      zdbt[ikl + 1] = zclear * zdbtmc + zcloud * zdbtmo;
      ((double*)ztdbt)[ikl] = zdbt[ikl + 1] * ((double*)ztdbt)[ikl + 1];
    }

    vrtqdr_sw_c_(&klev, zrefc, zrefdc, ztrac, ztradc, zdbtc, zrdndc, zrupc, zrupdc, ((double*)ztdbtc), (double*)zcd, (double*)zcu);
    vrtqdr_sw_c_(&klev, zref, zrefd, ztra, ztrad, zdbt, zrdnd, zrup, zrupd, ((double*)ztdbt), (double*)zfd, (double*)zfu);

    // vectorization : shuffle zfu, zfd etc.
    doublev4 zincflx_v4 = zincflx;
    doublev4 zincflxd5_v4 = zincflx_v4 * 0.5;
    for (jk = 0; jk < 8; jk ++){
      // jk : 0 ~ 31
      pbbfsu_T[ibm][jk] += zincflx_v4 * zfu[jk];
      pbbfsd_T[ibm][jk] += zincflx_v4 * zfd[jk];
      
      pbbfu[jk] += zincflx_v4 * zfu[jk];
      pbbfd[jk] += zincflx_v4 * zfd[jk];
      pbbcu[jk] += zincflx_v4 * zcu[jk];
      pbbcd[jk] += zincflx_v4 * zcd[jk];

      pbbfddir[jk] += zincflx_v4 * ztdbt[jk];
      pbbcddir[jk] += zincflx_v4 * ztdbtc[jk];

      if (ibm >= 9 && ibm <= 12) {
        puvcd[jk] += zincflx_v4 * zcd[jk];
        puvfd[jk] += zincflx_v4 * zfd[jk];

        puvfddir[jk] += zincflx_v4 * ztdbt[jk];
        puvcddir[jk] += zincflx_v4 * ztdbtc[jk];
      } else if (ibm == 8){
        puvcd[jk] += zincflxd5_v4 * zcd[jk];
        puvfd[jk] += zincflxd5_v4 * zfd[jk];
        pnicd[jk] += zincflxd5_v4 * zcd[jk];
        pnifd[jk] += zincflxd5_v4 * zfd[jk];
        
        puvfddir[jk] += zincflxd5_v4 * ztdbt[jk];
        puvcddir[jk] += zincflxd5_v4 * ztdbtc[jk];
        pnifddir[jk] += zincflxd5_v4 * ztdbt[jk];
        pnicddir[jk] += zincflxd5_v4 * ztdbtc[jk];

        pnicu[jk] += zincflxd5_v4 * zcu[jk];
        pnifu[jk] += zincflxd5_v4 * zfu[jk];
      } else if (ibm == 13 || ibm <= 7){
        pnicd[jk] += zincflx_v4 * zcd[jk];
        pnifd[jk] += zincflx_v4 * zfd[jk];

        pnifddir[jk] += zincflx_v4 * ztdbt[jk];
        pnicddir[jk] += zincflx_v4 * ztdbtc[jk];

        pnicu[jk] += zincflx_v4 * zcu[jk];
        pnifu[jk] += zincflx_v4 * zfu[jk];
      }
    }
  }
  
  // shuffle for specfic ibm
  ibm_last = -1;
  for(iw = iwb; iw < iwe; ++ iw){
    int ibm = iw2ibm[iw - iwb];
    if(ibm == ibm_last) continue;
    ibm_last = ibm;
    for(i = 0; i < 32; ++ i){
      ((double*)pbbfsu)[i * NBNDSW + ibm] = ((double*)pbbfsu_T)[ibm * 32 + i];
      ((double*)pbbfsd)[i * NBNDSW + ibm] = ((double*)pbbfsd_T)[ibm * 32 + i];
    }
  }
  //return;
  athread_syn(ARRAY_SCOPE, 0xffff);

  doublev4 get_buf[16];
  // register conmunication
  volatile int stplen;
  for(stplen = 1; stplen < 8; stplen <<= 1){
    volatile int lm = stplen - 1, nm = stplen;
    athread_syn(COL_SCOPE, 0xff);
    if(rid & lm) continue;
    if(rid & nm){
      for(i = 0; i < 8; ++ i){
        REG_PUTC(pbbcd[i]   , rid ^ nm); 
        REG_PUTC(pbbcu[i]   , rid ^ nm); 
        REG_PUTC(pbbfd[i]   , rid ^ nm); 
        REG_PUTC(pbbfu[i]   , rid ^ nm); 
        REG_PUTC(pbbcddir[i], rid ^ nm); 
        REG_PUTC(pbbfddir[i], rid ^ nm); 
        REG_PUTC(puvcd[i]   , rid ^ nm); 
        REG_PUTC(puvfd[i]   , rid ^ nm); 
        REG_PUTC(puvcddir[i], rid ^ nm); 
        REG_PUTC(puvfddir[i], rid ^ nm); 
        REG_PUTC(pnicd[i]   , rid ^ nm); 
        REG_PUTC(pnifd[i]   , rid ^ nm); 
        REG_PUTC(pnicddir[i], rid ^ nm); 
        REG_PUTC(pnifddir[i], rid ^ nm); 
        REG_PUTC(pnicu[i]   , rid ^ nm); 
        REG_PUTC(pnifu[i]   , rid ^ nm); 
      }
      for(i = 0; i < 8 * NBNDSW ; ++ i){
        REG_PUTC(pbbfsu[i], rid ^ nm); 
        REG_PUTC(pbbfsd[i], rid ^ nm); 
      }
    }
    else {
      for(i = 0; i < 8; ++ i){
        REG_GETC(get_buf[0]); pbbcd[i] += get_buf[0];
        REG_GETC(get_buf[1]); pbbcu[i] += get_buf[1];
        REG_GETC(get_buf[2]); pbbfd[i] += get_buf[2];
        REG_GETC(get_buf[3]); pbbfu[i] += get_buf[3];
        REG_GETC(get_buf[4]); pbbcddir[i] += get_buf[4];
        REG_GETC(get_buf[5]); pbbfddir[i] += get_buf[5];
        REG_GETC(get_buf[6]); puvcd[i] += get_buf[6];
        REG_GETC(get_buf[7]); puvfd[i] += get_buf[7];
        REG_GETC(get_buf[8]); puvcddir[i] += get_buf[8];
        REG_GETC(get_buf[9]); puvfddir[i] += get_buf[9];
        REG_GETC(get_buf[10]); pnicd[i] += get_buf[10];
        REG_GETC(get_buf[11]); pnifd[i] += get_buf[11];
        REG_GETC(get_buf[12]); pnicddir[i] += get_buf[12];
        REG_GETC(get_buf[13]); pnifddir[i] += get_buf[13];
        REG_GETC(get_buf[14]); pnicu[i] += get_buf[14];
        REG_GETC(get_buf[15]); pnifu[i] += get_buf[15];
      }
      for(i = 0; i < 8 * NBNDSW ; ++ i){
        REG_GETC(get_buf[0]); pbbfsu[i] += get_buf[0];
        REG_GETC(get_buf[1]); pbbfsd[i] += get_buf[1];
      }
    }
  }
  if(rid) return;

  for(stplen = 1; stplen < 8; stplen <<= 1){
    volatile int lm = stplen - 1, nm = stplen;  
    athread_syn(ROW_SCOPE, 0xff);
    if(cid & lm) continue;
    if(cid & nm){
      for(i = 0; i < 8; ++ i){
        REG_PUTR(pbbcd[i]   , cid ^ nm); 
        REG_PUTR(pbbcu[i]   , cid ^ nm); 
        REG_PUTR(pbbfd[i]   , cid ^ nm); 
        REG_PUTR(pbbfu[i]   , cid ^ nm); 
        REG_PUTR(pbbcddir[i], cid ^ nm); 
        REG_PUTR(pbbfddir[i], cid ^ nm); 
        REG_PUTR(puvcd[i]   , cid ^ nm); 
        REG_PUTR(puvfd[i]   , cid ^ nm); 
        REG_PUTR(puvcddir[i], cid ^ nm); 
        REG_PUTR(puvfddir[i], cid ^ nm); 
        REG_PUTR(pnicd[i]   , cid ^ nm); 
        REG_PUTR(pnifd[i]   , cid ^ nm); 
        REG_PUTR(pnicddir[i], cid ^ nm); 
        REG_PUTR(pnifddir[i], cid ^ nm); 
        REG_PUTR(pnicu[i]   , cid ^ nm); 
        REG_PUTR(pnifu[i]   , cid ^ nm); 
      }
      for(i = 0; i < 8 * NBNDSW ; ++ i){
        REG_PUTR(pbbfsu[i], cid ^ nm); 
        REG_PUTR(pbbfsd[i], cid ^ nm); 
      }
    }
    else {
      for(i = 0; i < 8; ++ i){
        REG_GETR(get_buf[0]); pbbcd[i] += get_buf[0];
        REG_GETR(get_buf[1]); pbbcu[i] += get_buf[1];
        REG_GETR(get_buf[2]); pbbfd[i] += get_buf[2];
        REG_GETR(get_buf[3]); pbbfu[i] += get_buf[3];
        REG_GETR(get_buf[4]); pbbcddir[i] += get_buf[4];
        REG_GETR(get_buf[5]); pbbfddir[i] += get_buf[5];
        REG_GETR(get_buf[6]); puvcd[i] += get_buf[6];
        REG_GETR(get_buf[7]); puvfd[i] += get_buf[7];
        REG_GETR(get_buf[8]); puvcddir[i] += get_buf[8];
        REG_GETR(get_buf[9]); puvfddir[i] += get_buf[9];
        REG_GETR(get_buf[10]); pnicd[i] += get_buf[10];
        REG_GETR(get_buf[11]); pnifd[i] += get_buf[11];
        REG_GETR(get_buf[12]); pnicddir[i] += get_buf[12];
        REG_GETR(get_buf[13]); pnifddir[i] += get_buf[13];
        REG_GETR(get_buf[14]); pnicu[i] += get_buf[14];
        REG_GETR(get_buf[15]); pnifu[i] += get_buf[15];
      }
      for(i = 0; i < 8 * NBNDSW ; ++ i){
        REG_GETR(get_buf[0]); pbbfsu[i] += get_buf[0];
        REG_GETR(get_buf[1]); pbbfsd[i] += get_buf[1];
      }
    }
  }
  if(cid) return;

  dma_set_size_wb(&put_desc, sizeof(double) * 32);
  dma_put(pm.pbbcd, pbbcd);
  dma_put(pm.pbbcu, pbbcu);
  dma_put(pm.pbbfd, pbbfd);
  dma_put(pm.pbbfu, pbbfu);
  dma_put(pm.pbbcddir, pbbcddir);
  dma_put(pm.pbbfddir, pbbfddir);
  dma_put(pm.puvcd, puvcd);
  dma_put(pm.puvfd, puvfd);
  dma_put(pm.puvcddir, puvcddir);
  dma_put(pm.puvfddir, puvfddir);
  dma_put(pm.pnicd, pnicd);
  dma_put(pm.pnifd, pnifd);
  dma_put(pm.pnicddir, pnicddir);
  dma_put(pm.pnifddir, pnifddir);
  dma_put(pm.pnicu, pnicu);
  dma_put(pm.pnifu, pnifu);
  dma_set_size_wb(&put_desc, sizeof(double) * 32 * NBNDSW);
  dma_put(pm.pbbfsu, pbbfsu);
  dma_put(pm.pbbfsd, pbbfsd);
}
