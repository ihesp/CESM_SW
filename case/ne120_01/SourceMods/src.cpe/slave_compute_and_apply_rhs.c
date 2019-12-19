#include "subFunForSlave.h"
/*kernel-CAAR-1*/
void slave_compute_and_apply_rhs_kernelcaar1a_(str_caar1* dat){
	/*local index variables*/
	volatile unsigned int id=athread_get_id(-1);
	volatile unsigned int cid,rid;GET_COL(cid);GET_ROW(rid);
	volatile int elemId=rid/2;
	volatile int nlevId=(id-elemId*16)*N;
	volatile unsigned long get_reply,put_reply;	
	int i,j,k,l,m,n,ie;

	/*get parameter:pointers & parameter*/
	str_caar1 dat0;
	dma_init();
	if(id==0){
		bcast_get(dat,&dat0,sizeof(str_caar1));
		dma_syn();
	}
	athread_syn(ARRAY_SCOPE,0xFFFF);

	double *ptr_e_v,*ptr_e_ps_v,*ptr_e_vn0,*ptr_e_Dinv,*ptr_e_metdet, 
				 *ptr_e_rmetdet,*ptr_e_D,*ptr_e_dp3d,*ptr_e_T,*ptr_e_Qdp, 
				 *ptr_e_phis,*ptr_e_phi,*ptr_e_eta_dot_dpdn,*ptr_e_omega_p, 
				 *ptr_e_spheremp,
				 *ptr_e_fcor,*ptr_e_pecnd,
				 *ptr_e_v_np1,*ptr_e_T_np1,*ptr_e_dp3d_np1,*ptr_e_ps_v_np1, 
				 *ptr_e_v_nm1,*ptr_e_T_nm1,*ptr_e_dp3d_nm1,*ptr_e_ps_v_nm1, 
				 *ptr_d_Dvv,*ptr_hvcoord;
	ptr_e_v            = dat0.ptr_e_v           ;
	ptr_e_ps_v         = dat0.ptr_e_ps_v        ;
	ptr_e_vn0          = dat0.ptr_e_vn0         ;
	ptr_e_Dinv         = dat0.ptr_e_Dinv        ;
	ptr_e_metdet       = dat0.ptr_e_metdet      ;
	ptr_e_rmetdet      = dat0.ptr_e_rmetdet     ;
	ptr_e_D            = dat0.ptr_e_D           ;
	ptr_e_dp3d         = dat0.ptr_e_dp3d        ;
	ptr_e_T            = dat0.ptr_e_T           ;
	ptr_e_Qdp          = dat0.ptr_e_Qdp         ;
	ptr_e_phi          = dat0.ptr_e_phi         ;
	ptr_e_phis         = dat0.ptr_e_phis        ;
	ptr_e_eta_dot_dpdn = dat0.ptr_e_eta_dot_dpdn;
	ptr_e_omega_p      = dat0.ptr_e_omega_p     ;
	ptr_e_spheremp     = dat0.ptr_e_spheremp    ;
	ptr_d_Dvv          = dat0.ptr_d_Dvv         ;
	ptr_hvcoord        = dat0.ptr_hvcoord       ;
	ptr_e_fcor         = dat0.ptr_e_fcor        ;
	ptr_e_pecnd        = dat0.ptr_e_pecnd       ;

	ptr_e_v_np1        = dat0.ptr_e_v_np1       ;
	ptr_e_T_np1        = dat0.ptr_e_T_np1       ;
	ptr_e_dp3d_np1     = dat0.ptr_e_dp3d_np1    ;
	ptr_e_ps_v_np1     = dat0.ptr_e_ps_v_np1    ;

	ptr_e_v_nm1        = dat0.ptr_e_v_nm1       ;
	ptr_e_T_nm1        = dat0.ptr_e_T_nm1       ;
	ptr_e_dp3d_nm1     = dat0.ptr_e_dp3d_nm1    ;
	ptr_e_ps_v_nm1     = dat0.ptr_e_ps_v_nm1    ;

	/*local index variables*/
	int len_elem;      
	len_elem           = dat0.len_elem          ;      
	int nets,nete,qn0,
			use_cpstar,qsplit,rsplit,n0,np1,nm1     ;
	nets               = dat0.nets              ;
	nete               = dat0.nete              ;
	qn0                = dat0.qn0               ;
	use_cpstar         = dat0.use_cpstar        ;
	qsplit             = dat0.qsplit            ;
	rsplit             = dat0.rsplit            ;
	n0                 = dat0.n0                ;
	np1                = dat0.np1               ; 
	nm1                = dat0.nm1               ;
	double eta_ave_w,rrearth,dt2                ;
	eta_ave_w          = dat0.eta_ave_w         ;
	rrearth            = dat0.rrearth           ;
	dt2                = dat0.dt2               ;
	double Cp,Rgas,Cpwater_vapor,Rwater_vapor,kappa;
	Cp                  = dat0.Cp               ;
	Rgas                = dat0.Rgas             ;
	Cpwater_vapor       = dat0.Cpwater_vapor    ;
	Rwater_vapor        = dat0.Rwater_vapor     ;
	kappa               = dat0.kappa            ;


	/*local variable */
	/*element variabels*/
	double e_dp3d[NP*NP*N],e_v[NP*NP*2*N],e_vn0[NP*NP*2*N];
	double e_ps_v[NP*NP],e_Dinv[NP*NP*2*2];
	double e_metdet[NP*NP],e_rmetdet[NP*NP];
	double e_D[NP*NP*2*2],e_T[NP*NP*N],e_Qdp[NP*NP*N];
	double e_omega_p[NP*NP*N];
	double e_spheremp[NP*NP];
	double e_pecnd[NP*NP*N],e_fcor[NP*NP];
	double e_phi[NP*NP*N];
	double e_phis[NP*NP];
	double e_ps_v_nm1[NP*NP],e_ps_v_np1[NP*NP];
	double e_T_np1[NP*NP*N],e_T_nm1[NP*NP*N];
	double e_v_np1[NP*NP*2*N],e_v_nm1[NP*NP*2*N];
	double e_dp3d_np1[NP*NP*N],e_dp3d_nm1[NP*NP*N]; 
	/*deriv variabels*/
	double d_Dvv[NP*NP];
	/*hvcoord variabels*/
	double hvcVar[2];
	/*
	 *     ps0           hvcVar[0  ]   1
	 *     hyai(plevp)   hvcVar[1  ]   31
	 *     hyam(plev)    hvcVar[32 ]   30
	 *     hybi(plevp)   hvcVar[62 ]   31
	 *     hybm(plev)    hvcVar[93 ]   30
	 *     hybd(plev)    hvcVar[123]   30
	 *     prsfac        hvcVar[153]   1
	 *     etam(plev)    hvcVar[154]   30
	 *     etai(plevp)   hvcVar[184]   31
	 */
	/*local tmp variabels*/
	double omega_p[NP*NP*N],vgrad_T[NP*NP],Ephi[NP*NP],grad_ps[NP*NP*2];
	double vtens1[NP*NP*N],vtens2[NP*NP*N],ttens[NP*NP*N];
	double E,glnps1,glnps2,gpterm;
	double grad_p[NP*NP*N*N],vgrad_p[NP*NP*N],vtemp[NP*NP*2],T_v[NP*NP*N];
	double shared_dp[NP*NP*L*N],shared_p[NP*NP*L*N],dp[NP*NP*N],p[NP*NP*N];
	double p1[NP*NP],dp1[NP*NP],rdp[NP*NP*N],vort[NP*NP*N],kappa_star[NP*NP*N];
	double v1,v2,Qt;
	double shared_phi[NP*NP*L*2],shared_T_v[NP*NP*L*2],divdp[NP*NP*N],shared_divdp[NP*NP*NLEV];

	double phii[NP*NP*L*2],phiiExchage[NP*NP];
	double hkk,hk1;
	double sumX[NP*NP*NLEV];
	double suml[NP*NP];
	/*for register communication*/
	doublev4  v4grad_ps;
	doublev4 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
	doublev4 v4C1,v4C2,v4C3,v4C4;
	v4C1=simd_set_doublev4(0.5,0.5,0.5,0.5);
	v4C2=simd_set_doublev4(2,2,2,2);
	v4C3=simd_set_doublev4(2*0.5,2*0.5,2*0.5,2*0.5);
	E=Rwater_vapor/Rgas-1.0;
	v4C4=simd_set_doublev4(E,E,E,E);

	/*BroadCast Global Variabels for each cpes*/ 
	if(id==0){
		get_reply=0;	
		bcast_get(ptr_hvcoord,hvcVar,sizeof(double)*2);
		bcast_get(ptr_d_Dvv,d_Dvv,sizeof(double)*NP*NP);
		dma_syn();
	}
	athread_syn(ARRAY_SCOPE,0xFFFF);
	int NUMELEM=nete-nets+1;
	int numEleLoop=(NUMELEM+3)/4;
	int mask[4]={0x03,0x0C,0x30,0xC0};
	elemId-=4;
	/*element loop*/
	for(ie=0;ie<numEleLoop;ie++){
		elemId+=4;
		/*Variabels     contains Nlev*/
		if(elemId<NUMELEM&&nlevId<29){
			pe_get(ptr_e_dp3d+elemId*len_elem+nlevId*NP*NP,e_dp3d,sizeof(double)*NP*NP*2  );
			pe_get(ptr_e_v   +elemId*len_elem+nlevId*NP*NP*2,e_v,   sizeof(double)*NP*NP*2*N);
			pe_get(ptr_e_vn0   +elemId*len_elem+nlevId*NP*NP*2,e_vn0,   sizeof(double)*NP*NP*2*N);
			pe_get(ptr_e_T   +elemId*len_elem+nlevId*NP*NP,e_T,   sizeof(double)*NP*NP*N);
			pe_get(ptr_e_Qdp    +elemId*len_elem+nlevId*NP*NP,e_Qdp,   sizeof(double)*NP*NP*N);
			pe_get(ptr_e_omega_p+elemId*len_elem+nlevId*NP*NP,e_omega_p,sizeof(double)*NP*NP*N);
			pe_get(ptr_e_pecnd+elemId*len_elem+nlevId*NP*NP,e_pecnd,sizeof(double)*NP*NP*N);
			pe_get(ptr_e_phi+elemId*len_elem+nlevId*NP*NP,e_phi,sizeof(double)*NP*NP*N);
			dma_syn();
		}

		/*Variabels Not contains Nlev*/
		if(nlevId==0&&elemId<NUMELEM){
			get_reply=0;
			athread_get(BCAST_MODE,ptr_e_Dinv+elemId*len_elem,e_Dinv,sizeof(double)*NP*NP*2*2,&get_reply,mask[elemId%4],0,0);
			athread_get(BCAST_MODE,ptr_e_metdet+elemId*len_elem,e_metdet,sizeof(double)*NP*NP,&get_reply,mask[elemId%4],0,0);
			athread_get(BCAST_MODE,ptr_e_rmetdet+elemId*len_elem,e_rmetdet,sizeof(double)*NP*NP,&get_reply,mask[elemId%4],0,0);
			athread_get(BCAST_MODE,ptr_e_fcor+elemId*len_elem,e_fcor,sizeof(double)*NP*NP,&get_reply,mask[elemId%4],0,0);
			athread_get(BCAST_MODE,ptr_e_D+elemId*len_elem,e_D,sizeof(double)*NP*NP*2*2,&get_reply,mask[elemId%4],0,0);
			athread_get(BCAST_MODE,ptr_e_spheremp+elemId*len_elem,e_spheremp,sizeof(double)*NP*NP,&get_reply,mask[elemId%4],0,0);
			athread_get(BCAST_MODE,ptr_e_phis+elemId*len_elem,e_phis,sizeof(double)*NP*NP,&get_reply,mask[elemId%4],0,0);
			while(get_reply!=7);
		}
		athread_syn(ARRAY_SCOPE,0xFFFF);
		/*Function Begin*/
		/*start of Gather--dp*/
		for(k=0;k<N;k++){
			for(i=0;i<NP;i++){
				simd_load(tmp1,e_dp3d+k*NP*NP+i*NP);
				simd_store(tmp1,dp+k*NP*NP+i*NP);
			}
		}

		if(cid==0){
			for(k=0;k<N;k++){
				for(i=0;i<NP;i++){
					simd_load(tmp1,        dp+k*NP*NP+i*NP);
					simd_store(tmp1,shared_dp+k*NP*NP+i*NP);
				}
			}
		}

		for(i=1;i<8;i++){
			if(cid==i){
				for(k=0;k<N*NP;k++){
					simd_load(v4grad_ps,dp+k*NP);
					REG_PUTR(v4grad_ps,0);
				}
			}
			if(cid==0){
				for(k=0;k<N*NP;k++){
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,shared_dp+i*NP*NP*N+k*NP);
				}
			}
			//athread_syn(ROW_SCOPE,0xFF);
			athread_syn(ARRAY_SCOPE,0xFFFF);
		}
		/*end of Gather--dp*/
		/*uppper part of each element*/
		if(nlevId==0){
			for(j=0;j<NP;j++){
				simd_load(tmp1,dp+j*NP);
				tmp2=hvcVar[1]*hvcVar[0]+tmp1/2.0;
				simd_store(tmp2,shared_p+j*NP);
			}

			for(k=1;k<2*L;k++){
				for(j=0;j<NP;j++){
					simd_load(tmp1,shared_p +(k-1)*NP*NP+j*NP);
					simd_load(tmp2,shared_dp+(k-1)*NP*NP+j*NP);
					tmp3=tmp1+tmp2/2.0;
					simd_load(tmp1,shared_dp+k*NP*NP+j*NP);
					tmp3+=tmp1/2.0;
					simd_store(tmp3,shared_p+k*NP*NP+j*NP);
				}
			}
		}
		/* register communication*/
		if(nlevId==0){
			for(j=0;j<NP;j++){
				simd_load(v4grad_ps,shared_p+(2*L-1)*NP*NP+j*NP);
				REG_PUTC(v4grad_ps,rid+1);
				simd_load(v4grad_ps,shared_dp+(2*L-1)*NP*NP+j*NP);
				REG_PUTC(v4grad_ps,rid+1);
			}
		}
		if(nlevId==16){
			for(j=0;j<NP;j++){
				REG_GETC(v4grad_ps);
				simd_store(v4grad_ps,p1+j*NP);
				REG_GETC(v4grad_ps);
				simd_store(v4grad_ps,dp1+j*NP);
			}
		}
		//athread_syn(ARRAY_SCOPE,0xFFFF);

		if(nlevId==16){
			//k=0;
			for(j=0;j<NP;j++){
				simd_load(tmp1,p1 +j*NP);
				simd_load(tmp2,dp1+j*NP);
				tmp3=tmp1+tmp2/2;
				simd_load(tmp1,shared_dp+j*NP);
				tmp3+=tmp1/2;
				simd_store(tmp3,shared_p+j*NP);
			}

			for(k=1;k<2*L;k++){
				for(j=0;j<NP;j++){
					simd_load(tmp1,shared_p+(k-1)*NP*NP+j*NP);
					simd_load(tmp2,shared_dp+(k-1)*NP*NP+j*NP);
					tmp3=tmp1+tmp2/2;
					simd_load(tmp1,shared_dp+k*NP*NP+j*NP);
					tmp3+=tmp1/2;
					simd_store(tmp3,shared_p+k*NP*NP+j*NP);
				}
			}
		}
		/*start of Scatter--p*/
		if(cid==0){
			for(k=0;k<N;k++){
				for(i=0;i<NP;i++){
					simd_load(tmp1,shared_p+k*NP*NP+i*NP);
					simd_store(tmp1,p+k*NP*NP+i*NP);
				}
			}
		}
		for(i=1;i<8;i++){
			if(cid==0){
				for(k=0;k<N*NP;k++){
					simd_load(v4grad_ps,shared_p+i*NP*NP*N+k*NP);
					REG_PUTR(v4grad_ps,i);
				}
			}
			if(cid==i){
				for(k=0;k<N*NP;k++){
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,p+k*NP);
				}
			}
			//athread_syn(ROW_SCOPE,0xFF);
			athread_syn(ARRAY_SCOPE,0xFFFF);
		}
		/*end of Scatter --p*/

		for(k=0;k<N;k++){
			gradient_sphere_fortran_(&p[k*NP*NP],&d_Dvv[0],&e_Dinv[0],&grad_p[k*NP*NP*2],&rrearth);
		}
		for(k=0;k<N;k++){
			for(i=0;i<NP;i++){
				simd_load(tmp1,dp+k*NP*NP+i*NP);
				tmp2=1.0/tmp1;
				simd_store(tmp1,rdp+k*NP*NP+i*NP);
			}
		}
		/* ============================
			 compute vgrad_lnps
			 ============================*/
		for(k=0;k<N;k++){
			for(i=0;i<NP;i++){
				simd_load(tmp1,e_v+k*NP*NP*2       +i*NP);
				simd_load(tmp2,e_v+k*NP*NP*2+NP*NP +i*NP);

				simd_load(tmp3,grad_p+k*NP*NP*2      +i*NP);
				tmp4=tmp1*tmp3;	
				simd_load(tmp3,grad_p+k*NP*NP*2+NP*NP+i*NP);
				tmp4+=(tmp2*tmp3);
				simd_store(tmp4,vgrad_p+k*NP*NP+i*NP);

				simd_load(tmp3,dp+k*NP*NP+i*NP);

				tmp4=tmp1*tmp3;
				simd_store(tmp4,vtemp+i*NP);
				simd_load(tmp5,e_vn0+k*NP*NP*2+i*NP);
				tmp5+=eta_ave_w*tmp4;
				simd_store(tmp5,e_vn0+k*NP*NP*2+i*NP);

				tmp4=tmp2*tmp3;
				simd_store(tmp4,vtemp+NP*NP+i*NP);
				simd_load(tmp5,e_vn0+k*NP*NP*2+NP*NP+i*NP);
				tmp5+=eta_ave_w*tmp4;
				simd_store(tmp5,e_vn0+k*NP*NP*2+NP*NP+i*NP);
			}
			/* =========================================
				 Compute relative vorticity and divergence
				 =========================================*/
			divergence_sphere_fortran_(vtemp,divdp+k*NP*NP,d_Dvv,e_metdet,e_Dinv,e_rmetdet,&rrearth);
			vorticity_sphere_fortran_(e_v+k*NP*NP*2,e_D,d_Dvv,e_rmetdet,vort+k*NP*NP,&rrearth);
		}

		for(k=0;k<N;k++){
			for(j=0;j<NP;j++){
				simd_load(tmp1,e_Qdp+k*NP*NP+j*NP);
				simd_load(tmp2,   dp+k*NP*NP+j*NP);
				simd_load(tmp3,  e_T+k*NP*NP+j*NP);

				tmp5=tmp3*(1.0+v4C4*(tmp1/tmp2));
				simd_store(tmp5,T_v+k*NP*NP+j*NP);

				tmp1=simd_set_doublev4(kappa,kappa,kappa,kappa);
				simd_store(tmp1,kappa_star+k*NP*NP+j*NP);
			}
		}
		/***************************/
		/***   preq_hydrostatic  ***/
		/*Broadcast grad_ps located on thread 0,16,32,4L to L,24,40,56 */
		/*start of Gather*/
		if(cid==0){
			for(k=0;k<N;k++){
				for(i=0;i<NP;i++){
					simd_load(tmp1,e_phi+k*NP*NP+i*NP);
					simd_load(tmp2,T_v  +k*NP*NP+i*NP);
					simd_load(tmp3,p    +k*NP*NP+i*NP);
					simd_load(tmp4,dp   +k*NP*NP+i*NP);

					simd_store(tmp1,shared_phi  +k*NP*NP+i*NP);
					simd_store(tmp2,shared_T_v  +k*NP*NP+i*NP);
					simd_store(tmp3,shared_p    +k*NP*NP+i*NP);
					simd_store(tmp4,shared_dp   +k*NP*NP+i*NP);
				}
			}
		}
		for(i=1;i<8;i++){
			if(cid==i){
				for(k=0;k<N*NP;k++){
					simd_load(v4grad_ps,e_phi+k*NP);
					REG_PUTR(v4grad_ps,0);
					simd_load(v4grad_ps,T_v+k*NP);
					REG_PUTR(v4grad_ps,0);
					simd_load(v4grad_ps,p+k*NP);
					REG_PUTR(v4grad_ps,0);
					simd_load(v4grad_ps,dp+k*NP);
					REG_PUTR(v4grad_ps,0);

				}
			}
			if(cid==0){
				for(k=0;k<N*NP;k++){
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,shared_phi+i*NP*NP*N+k*NP);
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,shared_T_v+i*NP*NP*N+k*NP);
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,shared_p+i*NP*NP*N+k*NP);
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,shared_dp+i*NP*NP*N+k*NP);

				}
			}
			athread_syn(ARRAY_SCOPE,0xFFFF);
		}
		/*end of Gather*/
		/*low part of each element */
		if(nlevId==16){
			for(j=0;j<NP;j++){
				//K=30
				simd_load(tmp1,shared_dp+(2*(L-2)+1)*NP*NP+j*NP);
				simd_load(tmp2,shared_p +(2*(L-2)+1)*NP*NP+j*NP);
				tmp3=tmp1*v4C1/tmp2;

				simd_load(tmp1,shared_T_v+(2*(L-2)+1)*NP*NP+j*NP);
				tmp2=Rgas*tmp1*v4C2*tmp3;
				simd_store(tmp2,phii+(2*(L-2)+1)*NP*NP      +j*NP);

				simd_load(tmp2,e_phis+j*NP);
				tmp3=tmp2+Rgas*tmp1*tmp3;
				simd_store(tmp3,shared_phi+(2*(L-2)+1)*NP*NP+j*NP);
			}
			for(j=0;j<NP;j++){
				//K=29->16
				for(k=12;k>=0;k--){
					simd_load(tmp1,shared_dp+k*NP*NP+j*NP);
					simd_load(tmp2,shared_p +k*NP*NP+j*NP);
					tmp3=tmp1*v4C1/tmp2;

					simd_load(tmp1,phii  +(k+1)*NP*NP+j*NP);
					simd_load(tmp2,shared_T_v+k*NP*NP+j*NP);
					tmp5=tmp1+Rgas*tmp2*v4C2*tmp3;
					simd_store(tmp5,phii+k*NP*NP+j*NP);

					simd_load(tmp4,e_phis+j*NP);
					tmp5=tmp4+tmp1+Rgas*tmp2*tmp3;
					simd_store(tmp5,shared_phi+k*NP*NP+j*NP);
				}
			}
		}
		for(i=0;i<NP;i++){
			//put phii[16] to neighbour along column
			if(nlevId==16){
				simd_load(v4grad_ps,phii+i*NP);
				REG_PUTC(v4grad_ps,rid-1);
			}
			/*uppper part of each element*/
			if(nlevId==0){
				/*receiver phii form ID=8 24 40 56 */
				REG_GETC(v4grad_ps);
				simd_store(v4grad_ps,phiiExchage+i*NP);
			}
			athread_syn(ARRAY_SCOPE,0xFFFF);
		}

		if(nlevId==0){
			//K=15
			for(j=0;j<NP;j++){
				simd_load(tmp1,shared_dp+(2*(L-1)+1)*NP*NP+j*NP);
				simd_load(tmp2,shared_p +(2*(L-1)+1)*NP*NP+j*NP);
				tmp3=tmp1*v4C1/tmp2;

				simd_load(tmp1,phiiExchage+j*NP);
				simd_load(tmp2,shared_T_v+(2*(L-1)+1)*NP*NP+j*NP);
				tmp5=tmp1+Rgas*tmp2*v4C2*tmp3;
				simd_store(tmp5,phii+(2*(L-1)+1)*NP*NP+j*NP);	

				simd_load(tmp4,e_phis+j*NP);
				tmp5=tmp4+tmp1+Rgas*tmp2*tmp3;
				simd_store(tmp5,shared_phi+(2*(L-1)+1)*NP*NP+j*NP);

			}
			//k=14->1
			for(j=0;j<NP;j++){
				for(k=14;k>0;k--){
					simd_load(tmp1,shared_dp+k*NP*NP+j*NP);
					simd_load(tmp2,shared_p +k*NP*NP+j*NP);
					tmp3=tmp1*v4C1/tmp2;

					simd_load(tmp1,phii+(k+1)*NP*NP+j*NP);
					simd_load(tmp2,shared_T_v+k*NP*NP+j*NP);
					tmp5=tmp1+Rgas*tmp2*v4C2*tmp3;
					simd_store(tmp5,phii+k*NP*NP+j*NP);

					simd_load(tmp4,e_phis+j*NP);
					tmp5=tmp4+tmp1+Rgas*tmp2*tmp3;
					simd_store(tmp5,shared_phi+k*NP*NP+j*NP);
				}
			}
			//K=0
			for(j=0;j<NP;j++){
				simd_load(tmp1,shared_dp+j*NP);
				simd_load(tmp2,shared_p +j*NP);
				tmp3=tmp1*v4C1/tmp2;

				simd_load(tmp1,phii+1*NP*NP+j*NP);
				simd_load(tmp2,shared_T_v+j*NP);
				simd_load(tmp4,e_phis+j*NP);

				tmp5=tmp4+tmp1+Rgas*tmp2*tmp3;
				simd_store(tmp5,shared_phi+j*NP);
			}
		}
		/*start of Scatter*/
		if(cid==0){
			for(k=0;k<N;k++){
				for(i=0;i<NP;i++){
					simd_load(tmp1,shared_phi+k*NP*NP+i*NP);
					simd_load(tmp2,shared_T_v+k*NP*NP+i*NP);
					simd_load(tmp3,shared_p+k*NP*NP+i*NP);
					simd_load(tmp4,shared_dp+k*NP*NP+i*NP);

					simd_store(tmp1,e_phi+k*NP*NP+i*NP);
					simd_store(tmp2,T_v+k*NP*NP+i*NP);
					simd_store(tmp3,p+k*NP*NP+i*NP);
					simd_store(tmp4,dp+k*NP*NP+i*NP);

				}
			}
		}
		for(k=0;k<N*NP;k++){
			for(i=1;i<8;i++){
				if(cid==0){
					simd_load(v4grad_ps,shared_phi+i*NP*NP*N+k*NP);
					REG_PUTR(v4grad_ps,i);
				}
				if(cid==i){
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,e_phi+k*NP);
				}
			}
			athread_syn(ARRAY_SCOPE,0xFFFF);
			//athread_syn(ROW_SCOPE,0xFF);
		}
		/*end of Scatter*/
		/**** preq_hydrostatic  ****/
		/***************************/

		/***************************/
		/******  preq_omega_ps  ****/ 
		/*start of Gather--devdp*/
		if(cid==0){
			for(k=0;k<N;k++){
				for(i=0;i<NP;i++){
					simd_load(tmp1,divdp+k*NP*NP+i*NP);
					simd_store(tmp1,shared_divdp+k*NP*NP+i*NP);
				}
			}
		}
		for(i=1;i<8;i++){
			if(cid==i){
				for(k=0;k<N*NP;k++){
					simd_load(v4grad_ps,divdp+k*NP);
					REG_PUTR(v4grad_ps,0);
				}
			}
			if(cid==0){
				for(k=0;k<N*NP;k++){
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,shared_divdp+i*NP*NP*N+k*NP);
				}
			}
			athread_syn(ARRAY_SCOPE,0xFFFF);
		}

		/*rid=1,3,5,7->rid=0,2,4,6*/
		if(nlevId==16){
			for(k=0;k<NP*L*N;k++){
				simd_load(v4grad_ps,shared_divdp+k*NP);
				REG_PUTC(v4grad_ps,rid-1);
			}
		}
		if(nlevId==0){
			for(k=0;k<NP*L*N;k++){
				REG_GETC(v4grad_ps);
				simd_store(v4grad_ps,shared_divdp+NP*NP*N*L+k*NP);
			}
		}
		//athread_syn(ARRAY_SCOPE,0xFFFF);

		/*end of Gather--divdp*/
		if(nlevId==0){
			for(j=0;j<NP;j++){
				/*  k=1  */
				simd_load(tmp1,shared_divdp+j*NP);
				simd_store(tmp1,suml       +j*NP);
				simd_store(tmp1,sumX       +j*NP);
				/* k=2->29 */
				for(k=1;k<NLEV-1;k++){
					simd_load(tmp1,suml+j*NP);
					simd_store(tmp1,sumX+k*NP*NP+j*NP);

					simd_load(tmp2,shared_divdp+k*NP*NP+j*NP);
					tmp1+=tmp2;
					simd_store(tmp1,suml+j*NP);
				}
				/*  k=30 */
				simd_load(tmp1,suml+j*NP);
				simd_store(tmp1,sumX+(NLEV-1)*NP*NP+j*NP);
			}
		}
		/* start of Scatter--sumX*/
		if(nlevId==0){
			for(k=0;k<NP*L*N;k++){
				simd_load(v4grad_ps,sumX+NP*NP*N*L+k*NP);
				REG_PUTC(v4grad_ps,rid+1);
			}
		}
		if(nlevId==16){
			for(k=0;k<NP*L*N;k++){
				REG_GETC(v4grad_ps);
				simd_store(v4grad_ps,sumX+k*NP);
			}
		}
		//athread_syn(ARRAY_SCOPE,0xFFFF);
		for(i=1;i<8;i++){
			if(cid==0){
				for(k=0;k<N*NP;k++){
					simd_load(v4grad_ps,sumX+i*NP*NP*N+k*NP);
					REG_PUTR(v4grad_ps,i);
				}
			}
			if(cid==i){
				for(k=0;k<N*NP;k++){
					REG_GETR(v4grad_ps);
					simd_store(v4grad_ps,sumX+k*NP);
				}
			}
			athread_syn(ARRAY_SCOPE,0xFFFF);
		}
		/* end of Scatter--sumX*/

		if(nlevId==0){
			/*  k=0  */
			for(j=0;j<NP;j++){
				simd_load(tmp1,vgrad_p+j*NP);
				simd_load(tmp2,p      +j*NP);
				simd_load(tmp3,divdp  +j*NP);
				tmp4=tmp1/tmp2-v4C1/tmp2*tmp3;
				simd_store(tmp4,omega_p+j*NP);	

				simd_load(tmp2,e_omega_p+j*NP);
				tmp2+=(eta_ave_w*tmp4);
				simd_store(tmp2,e_omega_p+j*NP);
			}
			/* k=1  */
			for(j=0;j<NP;j++){
				simd_load(tmp1,vgrad_p+NP*NP+j*NP);
				simd_load(tmp2,p      +NP*NP+j*NP);
				simd_load(tmp3,sumX   +NP*NP+j*NP);
				simd_load(tmp4,divdp  +NP*NP+j*NP);
				tmp5=tmp1/tmp2-tmp3*(v4C3/tmp2)-v4C1/tmp2*tmp4;
				simd_store(tmp5,omega_p+NP*NP+j*NP);
				simd_load(tmp2,e_omega_p+NP*NP+j*NP);
				tmp2+=(eta_ave_w*tmp5);
				simd_store(tmp2,e_omega_p+NP*NP+j*NP);
			}
		}
		else{
			for(k=0;k<N;k++){
				for(j=0;j<NP;j++){
					simd_load(tmp1,vgrad_p+k*NP*NP+j*NP);
					simd_load(tmp2,p      +k*NP*NP+j*NP);
					simd_load(tmp3,sumX   +k*NP*NP+j*NP);
					simd_load(tmp4,divdp  +k*NP*NP+j*NP);
					tmp5=tmp1/tmp2-tmp3*(v4C3/tmp2)-v4C1/tmp2*tmp4;
					simd_store(tmp5,omega_p+k*NP*NP+j*NP);
					simd_load(tmp2,e_omega_p+k*NP*NP+j*NP);
					tmp2+=(eta_ave_w*tmp5);
					simd_store(tmp2,e_omega_p+k*NP*NP+j*NP);

				}
			}
		}
		/******  preq_omega_ps  ****/ 
		/***************************/

		/*==============================================
			! Compute phi + kinetic energy term: 10*nv*nv Flops
			! ==============================================*/
		for(k=0;k<N;k++){
			for(j=0;j<NP;j++){
				simd_load(tmp1,e_v    +k*NP*NP*2+      j*NP);
				simd_load(tmp2,e_v    +k*NP*NP*2+NP*NP+j*NP);
				simd_load(tmp3,e_phi  +k*NP*NP+        j*NP);
				simd_load(tmp4,e_pecnd+k*NP*NP+        j*NP);
				tmp5=v4C1*(tmp1*tmp1+tmp2*tmp2)+tmp3+tmp4;
				simd_store(tmp5,Ephi+j*NP);
			}

			gradient_sphere_fortran_(e_T+k*NP*NP,d_Dvv,e_Dinv,vtemp,&rrearth);
			//gradient_sphere(e_T+k*NP*NP,d_Dvv,e_Dinv,vtemp,rrearth);


			for(j=0;j<NP;j++){
				simd_load(tmp1,e_v+k*NP*NP*2+      j*NP);
				simd_load(tmp2,e_v+k*NP*NP*2+NP*NP+j*NP);
				simd_load(tmp3,vtemp+j*NP);
				simd_load(tmp4,vtemp+NP*NP+j*NP);
				tmp5=tmp1*tmp3+tmp2*tmp4;
				simd_store(tmp5,vgrad_T+j*NP);
			}

			gradient_sphere_fortran_(Ephi,d_Dvv,e_Dinv,vtemp,&rrearth);
			//gradient_sphere(Ephi,d_Dvv,e_Dinv,vtemp,rrearth);
			for(j=0;j<NP;j++){
				simd_load(tmp1,T_v+k*NP*NP+j*NP);
				simd_load(tmp2,p  +k*NP*NP+j*NP);
				tmp1=tmp1/tmp2;//gpterm
				simd_load(tmp2,grad_p+k*NP*NP*2+j*NP);
				tmp2=Rgas*tmp1*tmp2;//glnps1
				//vtens1
				simd_load(tmp3,e_fcor+j*NP);
				simd_load(tmp4,vort+k*NP*NP+j*NP);
				simd_load(tmp5,e_v+k*NP*NP*2+NP*NP+j*NP);//v2
				tmp5=tmp5*(tmp3+tmp4);	
				//simd_load(tmp7,v_vadv+k*NP*NP*2+j*NP);
				tmp7=simd_set_doublev4(0,0,0,0);
				simd_load(tmp8,vtemp+j*NP);
				tmp5=-tmp7+tmp5-tmp8-tmp2;
				simd_store(tmp5,vtens1+k*NP*NP+j*NP);
				//vtens2
				simd_load(tmp2,grad_p+k*NP*NP*2+NP*NP+j*NP);
				tmp2=Rgas*tmp1*tmp2;//glnps2
				simd_load(tmp5,e_v+k*NP*NP*2+j*NP);//v1

				tmp5=tmp5*(tmp3+tmp4);	
				//simd_load(tmp3,v_vadv+k*NP*NP*2+NP*NP+j*NP);
				tmp3=simd_set_doublev4(0,0,0,0);
				simd_load(tmp4,vtemp+NP*NP+j*NP);
				tmp5=-tmp3-tmp5-tmp4-tmp2;
				simd_store(tmp5,vtens2+k*NP*NP+j*NP);
				//ttens
				//simd_load(tmp1,T_vadv +k*NP*NP+j*NP);
				tmp1=simd_set_doublev4(0,0,0,0);
				simd_load(tmp2,vgrad_T+        j*NP);
				simd_load(tmp3,kappa_star+k*NP*NP+j*NP);
				simd_load(tmp4,T_v+k*NP*NP+j*NP);
				simd_load(tmp5,omega_p+k*NP*NP+j*NP);
				tmp6=-tmp1-tmp2+tmp3*tmp4*tmp5;
				simd_store(tmp6,ttens+k*NP*NP+j*NP);
			}
		}
		/*update host data*/
		if(elemId<NUMELEM&&nlevId<29){
			put_reply=0;
			athread_put(PE_MODE,e_vn0,ptr_e_vn0   +elemId*len_elem+nlevId*NP*NP*2,
					sizeof(double)*NP*NP*2*N,&put_reply,0,0);	
			athread_put(PE_MODE,e_phi,ptr_e_phi+elemId*len_elem+nlevId*NP*NP,
					sizeof(double)*NP*NP*N,&put_reply,0,0);
			athread_put(PE_MODE,e_omega_p,ptr_e_omega_p+elemId*len_elem+nlevId*NP*NP,
					sizeof(double)*NP*NP*N,&put_reply,0,0);
			while(put_reply!=3);
		}

		/* =========================================================
			 local element timestep, store in np1.
			 note that we allow np1=n0 or nm1
			 apply mass matrix
			 =========================================================*/
		if(elemId<NUMELEM&&nlevId<29){
			if(nm1==n0){
				for(k=0;k<N;k++){
					for(j=0;j<NP;j++){
						simd_load(tmp1,e_v     +k*NP*NP*2      +j*NP);
						simd_store(tmp1,e_v_nm1+k*NP*NP*2      +j*NP);
						simd_load(tmp1,e_v     +k*NP*NP*2+NP*NP+j*NP);
						simd_store(tmp1,e_v_nm1+k*NP*NP*2+NP*NP+j*NP);

						simd_load(tmp1,e_T     +k*NP*NP+j*NP);
						simd_store(tmp1,e_T_nm1+k*NP*NP+j*NP);

						simd_load(tmp1,e_dp3d     +k*NP*NP+j*NP);
						simd_store(tmp1,e_dp3d_nm1+k*NP*NP+j*NP);
					}
				}
			}
			else{
				/*--1---update e_v_np1*/
				pe_get(ptr_e_v_nm1   +elemId*len_elem+nlevId*NP*NP*2,e_v_nm1,sizeof(double)*NP*NP*2*N);
				/*--2---update e_T_np1*/
				pe_get(ptr_e_T_nm1   +elemId*len_elem+nlevId*NP*NP,e_T_nm1,sizeof(double)*NP*NP*N);
				/*--3---update e_dp3d_np1*/
				pe_get(ptr_e_dp3d_nm1+elemId*len_elem+nlevId*NP*NP,e_dp3d_nm1,sizeof(double)*NP*NP*N);
				dma_syn();	
			}
			for(k=0;k<N;k++){
				for(j=0;j<NP;j++){
					/*--1---update e_v_np1*/
					simd_load(tmp1,e_spheremp+j*NP);
					simd_load(tmp2,e_v_nm1+k*NP*NP*2+j*NP);
					simd_load(tmp3,vtens1 +k*NP*NP  +j*NP);
					tmp4=tmp1*(tmp2+dt2*tmp3);
					simd_store(tmp4,e_v_np1+k*NP*NP*2+j*NP);

					simd_load(tmp2,e_v_nm1+k*NP*NP*2+NP*NP+j*NP);
					simd_load(tmp3,vtens2 +k*NP*NP  +j*NP);
					tmp4=tmp1*(tmp2+dt2*tmp3);
					simd_store(tmp4,e_v_np1+k*NP*NP*2+NP*NP+j*NP);
					/*--2---update e_T_np1*/
					simd_load(tmp2,e_T_nm1+k*NP*NP+j*NP);
					simd_load(tmp3,ttens +k*NP*NP  +j*NP);
					tmp4=tmp1*(tmp2+dt2*tmp3);
					simd_store(tmp4,e_T_np1+k*NP*NP+j*NP);

					/*--3---update e_dp3d_np1*/
					simd_load(tmp2,e_dp3d_nm1+k*NP*NP+j*NP);
					simd_load(tmp3,divdp +k*NP*NP  +j*NP);
					tmp4=tmp1*(tmp2-dt2*tmp3);
					simd_store(tmp4,e_dp3d_np1+k*NP*NP+j*NP);
				}
			}
			/*update e_v_np1,e_T_np1,e_dp3d_np1 to host*/
			put_reply=0;
			athread_put(PE_MODE,e_v_np1,ptr_e_v_np1+elemId*len_elem+nlevId*NP*NP*2,
					sizeof(double)*NP*NP*2*N,&put_reply,0,0);
			athread_put(PE_MODE,e_T_np1,  ptr_e_T_np1   +elemId*len_elem+nlevId*NP*NP,
					sizeof(double)*NP*NP*N,&put_reply,0,0);
			athread_put(PE_MODE,e_dp3d_np1,ptr_e_dp3d_np1+elemId*len_elem+nlevId*NP*NP,
					sizeof(double)*NP*NP*2  ,&put_reply,0,0);
			while(put_reply!=3);
		}

		if(elemId<NUMELEM&&nlevId==0){
			get_reply=0;
			athread_get(PE_MODE,ptr_e_ps_v_nm1+elemId*len_elem,e_ps_v_nm1,sizeof(double)*NP*NP,&get_reply,0,0,0);
			while(get_reply!=1);
			for(j=0;j<NP;j++){
				simd_load(tmp1,e_spheremp+j*NP);
				simd_load(tmp2,e_ps_v_nm1+j*NP);
				tmp3=tmp1*tmp2;
				simd_store(tmp3,e_ps_v_np1+j*NP);
			}
			put_reply=0;
			athread_put(PE_MODE,e_ps_v_np1,ptr_e_ps_v_np1+elemId*len_elem,
					sizeof(double)*NP*NP,&put_reply,0,0);
			while(put_reply!=1);
		}
		/*Function End*/
	}
}

