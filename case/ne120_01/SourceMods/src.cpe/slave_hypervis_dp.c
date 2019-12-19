#include<slave.h> 
#include"register.h"
#include<simd.h>
#define NP 4
#define N 8

#define ADDRESS2(p,a,len1,b,len2) p+a*len1+b*len2 //caculate memory address 

/* Fortran function for kernel a&d and g&h */
extern void vlaplace_sphere_wk_ad_fortran_(double* v,double* laplace,
		double* d_Dvv,double* e_metdet,double* e_Dinv,
		double* e_rmetdet,double* rrearth,double* e_D,
		double* e_mp,double* e_metinv,double* e_spheremp);
extern void  vlaplace_sphere_wk_gh_fortran_(double* v,double* laplace,
		double* e_vec_sphere2cart,double* hypervis_power,
		double* hypervis_scaling,int* var_coef,double* nu_ratio,double* rrearth,
		/*part-lap*/
		double* e_variable_hyperviscosity,double* e_tensorVisc,
		double* d_Dvv,double* e_Dinv,double* e_spheremp,
		/*part-Vlap*/
		double* e_metdet,double* e_rmetdet,
		double* e_D,double* e_mp,double* e_metinv);

/*Kernel-A*/
typedef struct{
	/*pointer-0*/
	double*v,*T,*e_spheremp,*lap_t,*lap_v;
	double nu_s,nu,dt;
	/*pointer-laplace_sphere_wk*/
	double* d_Dvv,*e_Dinv,rrearth;
	/*pointer-Vlaplace_sphere_wk*/
	double* e_metdet,*e_rmetdet,*e_D,*e_mp,*e_metinv;
	/*parameter-0*/
	int nets,nete,nlev,np,len_elem
}str_a;

/* Laplace-Function */
void gradient_sphere(int np ,double s[np][np],double d_Dvv[np][np],
		double Dinv[2][2][np][np],double ds[2][np][np],double rrearth){
	/*local variable*/
	int i,j,l;
	double dsdx00,dsdy00,v1[np][np],v2[np][np];
	for(j=0;j<np;j++){
		for(l=0;l<np;l++){
			dsdx00=0.0;
			dsdy00=0.0;
			for(i=0;i<np;i++){
				dsdx00 = dsdx00 +d_Dvv[l][i]*s[j][i];
				dsdy00 = dsdy00 +d_Dvv[l][i]*s[i][j];
			}
			v1[j][l]= dsdx00*rrearth;
			v2[l][j]= dsdy00*rrearth;
		}
	}
	for(j=0;j<np;j++){
		for(i=0;i<np;i++){
			ds[0][j][i]=Dinv[0][0][j][i]*v1[j][i]+Dinv[0][1][j][i]*v2[j][i];
			ds[1][j][i]=Dinv[1][0][j][i]*v1[j][i]+Dinv[1][1][j][i]*v2[j][i];
		}
	}
}
void divergence_sphere_wk(int np,double v[2][np][np],double d_Dvv[np][np],
		double e_Dinv[2][2][np][np],double e_rspheremp[np][np],double div[np][np],double rrearth){
	/*local variable*/
	double vtemp[2][np][np];
	int i,j,m,n;
	for(j=0;j<np;j++){
		for(i=0;i<np;i++)
		{
			vtemp[0][j][i]=(e_Dinv[0][0][j][i]*v[0][j][i] + e_Dinv[1][0][j][i]*v[1][j][i]);
			vtemp[1][j][i]=(e_Dinv[0][1][j][i]*v[0][j][i] + e_Dinv[1][1][j][i]*v[1][j][i]);
		}
	}
	for(n=0;n<np;n++)
		for(m=0;m<np;m++){
			div[n][m]=0;
			for(j=0;j<np;j++){
				div[n][m]=div[n][m]-(e_rspheremp[n][j]*vtemp[0][n][j]*d_Dvv[j][m]+
						e_rspheremp[j][m]*vtemp[1][j][m]*d_Dvv[j][n])*rrearth;
			}
		}
}
void laplace_sphere_wk(int np,double s[np][np],double d_Dvv[np][np],double e_Dinv[2][2][np][np],
		double e_rspheremp[np][np],double laplace[np][np],double rrearth){
	/*local variable*/
	double grads[2][np][np];

	gradient_sphere(np,s,d_Dvv,e_Dinv,grads,rrearth);
	divergence_sphere_wk(np,grads,d_Dvv,e_Dinv,e_rspheremp,laplace,rrearth);
}

void slave_advance_hypervis_dp_kernela_(str_a* dat){	
	/*receive  variable and pointer from Fortran*/
	double *ptr_v,*ptr_T,*ptr_e_spheremp,*ptr_lap_t,*ptr_lap_v;
	ptr_v=dat->v;ptr_T=dat->T;
	ptr_e_spheremp=dat->e_spheremp;
	ptr_lap_t=dat->lap_t;
	ptr_lap_v=dat->lap_v;
	double nu_s,nu,dt;
	nu_s=dat->nu_s;nu=dat->nu;dt=dat->dt;
	/*laplace_sphere_wk*/
	double* ptr_d_Dvv,*ptr_e_Dinv,rrearth;
	ptr_d_Dvv=dat->d_Dvv; ptr_e_Dinv=dat->e_Dinv;
	rrearth=dat->rrearth;
	/*Vlaplace_sphere_wk*/
	double* ptr_e_metdet,*ptr_e_rmetdet,*ptr_e_D,*ptr_e_mp,*ptr_e_metinv;
	ptr_e_metdet=dat->e_metdet;
	ptr_e_rmetdet=dat->e_rmetdet;
	ptr_e_D=dat->e_D;
	ptr_e_mp=dat->e_mp;
	ptr_e_metinv=dat->e_metinv;
	int nets,nete,nlev,np,len_elem;
	nets=dat->nets;nete=dat->nete;
	nlev=dat->nlev;np=dat->np;
	len_elem=dat->len_elem;
	/*get threads id*/
	volatile unsigned int id=athread_get_id(-1);
	/*local parameter */
	volatile unsigned long get_reply,put_reply;
	int GridX=64;
	int Nsub=nlev;
	int LoopNumX=((nete-nets+1)+GridX-1)/GridX;
	int colid=id%GridX;
	int globalx,globaly,i,j,n,iloopx;
	/*main LDM*/
	double v[Nsub][2][np][np],T[Nsub][np][np],e_spheremp[np][np],lap_v[2][np][np],lap_t[np][np];//main
	double d_Dvv[np][np],e_Dinv[2][2][np][np];//laplace_sphere_wk
	double e_metdet[np][np],e_rmetdet[np][np],//Vlaplace_sphere_wk
				 e_D[2][2][np][np],e_mp[np][np],
				 e_metinv[2][2][np][np];
	/*block loop*/
	for(iloopx=0;iloopx<LoopNumX;iloopx++){
		globalx=colid+iloopx*GridX;
		globaly=0;
		if(globalx<(nete-nets+1)){
			/*main*/
			get_reply=0;
			athread_get(PE_MODE,ADDRESS2(ptr_v,globalx,len_elem,globaly,np*np*2),
					v,sizeof(double)*np*np*2*Nsub,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_T,globalx,len_elem,globaly,np*np),T
					,sizeof(double)*np*np*Nsub,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_spheremp,globalx,len_elem,globaly
						,0),e_spheremp,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ptr_lap_t,lap_t,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ptr_lap_v,lap_v,sizeof(double)*np*np*2,&get_reply,0,0,0);
			/*lap*/
			athread_get(PE_MODE,ADDRESS2(ptr_e_Dinv,globalx,len_elem,0,0),e_Dinv,
					sizeof(double)*np*np*2*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ptr_d_Dvv,d_Dvv,sizeof(double)*np*np,&get_reply,0,0,0);
			/*Vlap*/
			athread_get(PE_MODE,ADDRESS2(ptr_e_metdet,globalx,len_elem,0
						,0),e_metdet,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_rmetdet,globalx,len_elem,0
						,0),e_rmetdet,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_D,globalx,len_elem,0
						,0),e_D,sizeof(double)*np*np*2*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_mp,globalx,len_elem,0
						,0),e_mp,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_metinv,globalx,len_elem,0
						,0),e_metinv,sizeof(double)*np*np*2*2,&get_reply,0,0,0);
			while(get_reply!=12);
			/*caculate*/
			for(n=0;n<Nsub;n++){
				laplace_sphere_wk(np,T[n],d_Dvv,e_Dinv,e_spheremp,lap_t,rrearth);
				vlaplace_sphere_wk_ad_fortran_(&v[n][0][0][0],&lap_v[0][0][0],&d_Dvv[0][0],&e_metdet[0][0],
						&e_Dinv[0][0][0][0],&e_rmetdet[0][0],&rrearth,&e_D,
						&e_mp[0][0],&e_metinv[0][0][0][0],&e_spheremp[0][0]);
				for(i=0;i<np;i++){
					for(j=0;j<np;j++){
						T[n][i][j]=e_spheremp[i][j]*T[n][i][j]+dt*nu_s*lap_t[i][j];
						v[n][0][i][j]=e_spheremp[i][j]*v[n][0][i][j]+dt*nu*lap_v[0][i][j];
						v[n][1][i][j]=e_spheremp[i][j]*v[n][1][i][j]+dt*nu*lap_v[1][i][j];
					}
				}
			}
			/*save data to host */
			put_reply=0;	
			athread_put(PE_MODE,v,ADDRESS2(ptr_v,globalx,len_elem,globaly,np*np*2),
					sizeof(double)*np*np*2*Nsub,&put_reply,0,0);
			athread_put(PE_MODE,T,ADDRESS2(ptr_T,globalx,len_elem,globaly,np*np),
					sizeof(double)*np*np*Nsub,&put_reply,0,0);
			while(put_reply!=2);
		}
	}
}

/*Kernel-B*/
typedef struct{
	double*v,*T,*rspheremp;
	int nets,nete,nlev,np,len_elem;
}str_b;

void slave_advance_hypervis_dp_kernelb_(str_b* dat){	
	/*receive  variable and pointer from Fortran*/
	double *ptr_v,*ptr_T,*ptr_rspheremp;
	ptr_v=dat->v;ptr_T=dat->T;
	ptr_rspheremp=dat->rspheremp;
	int nets,nete,nlev,np,len_elem;
	nets=dat->nets;nete=dat->nete;
	nlev=dat->nlev;np=dat->np;
	len_elem=dat->len_elem;
	/*get threads id*/
	volatile unsigned int id=athread_get_id(-1);
	/*local parameter */
	volatile unsigned long get_reply,put_reply;
	int GridX=64;
	int Nsub=nlev;
	int LoopNumX=((nete-nets+1)+GridX-1)/GridX;
	int colid=id%GridX;
	int globalx,globaly,i,j,n,iloopx;
	/*main LDM*/
	double v[Nsub][2][np][np],T[Nsub][np][np],rspheremp[np][np];
	/*block loop */
	for(iloopx=0;iloopx<LoopNumX;iloopx++){
		globalx=colid+iloopx*GridX;
		globaly=0;
		if(globalx<(nete-nets+1)){
			get_reply=0;
			athread_get(PE_MODE,ADDRESS2(ptr_v,globalx,len_elem,globaly,np*np*2),
					v,Nsub*sizeof(double)*np*np*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_T,globalx,len_elem,globaly,np*np),T
					,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_rspheremp,globalx,len_elem,globaly
						,0),rspheremp,sizeof(double)*np*np,&get_reply,0,0,0);
			while(get_reply!=3);
			/*caculate*/
			for(n=0;n<Nsub;n++){
				for(i=0;i<np;i++){
					for(j=0;j<np;j++){
						T[n][i][j]	 =rspheremp[i][j]*T[n][i][j];
						v[n][0][i][j]=rspheremp[i][j]*v[n][0][i][j];
						v[n][1][i][j]=rspheremp[i][j]*v[n][1][i][j];
					}
				}
			}
			/*save data to host */
			put_reply=0;	
			athread_put(PE_MODE,v,ADDRESS2(ptr_v,globalx,len_elem,globaly,np*np*2),
					sizeof(double)*np*np*2*Nsub,&put_reply,0,0);
			athread_put(PE_MODE,T,ADDRESS2(ptr_T,globalx,len_elem,globaly,np*np),
					sizeof(double)*np*np*Nsub,&put_reply,0,0);
			while(put_reply!=2);
		}
	}
}

/*Kernel-C*/
typedef struct{
	double *dp3d,*dpdiss_ave,*dpdiss_biharmonic,*dptens;
	double eta_ave_w;
	int hypervis_subcycle,nets,nete,nlev,np,len_elem,len_dptens;
}str_c;

void slave_advance_hypervis_dp_kernelc_(str_c* dat){	
	/*receive  variable and pointer from Fortran*/
	double*ptr_dp3d,*ptr_dpdiss_ave,*ptr_dpdiss_biharmonic,*ptr_dptens;
	ptr_dp3d=dat->dp3d;
	ptr_dpdiss_ave=dat->dpdiss_ave;
	ptr_dpdiss_biharmonic=dat->dpdiss_biharmonic;
	ptr_dptens=dat->dptens;
	double eta_ave_w;
	eta_ave_w=dat->eta_ave_w;
	int hypervis_subcycle,nets,nete,nlev,np,len_elem,len_dptens;
	hypervis_subcycle=dat->hypervis_subcycle;
	nets=dat->nets;nete=dat->nete;
	nlev=dat->nlev;np=dat->np;
	len_elem=dat->len_elem;
	len_dptens=dat->len_dptens;
	/*get threads id*/
	volatile unsigned int id=athread_get_id(-1);
	/*local parameter*/
	volatile unsigned long get_reply,put_reply;
	int GridX=64;
	int Nsub=nlev;
	int LoopNumX=((nete-nets+1)+GridX-1)/GridX;
	int colid=id%GridX;
	int globalx,globaly,i,j,n,iloopx;
	/*main LDM*/
	double dpdiss_ave[Nsub][np][np],
				 dp3d[Nsub][np][np],dpdiss_biharmonic[Nsub][np][np],
				 dptens[Nsub][np][np];
	/*block loop*/
	for(iloopx=0;iloopx<LoopNumX;iloopx++){
		globalx=colid+iloopx*GridX;
		globaly=0;
		if(globalx<(nete-nets+1)){
			get_reply=0;
			athread_get(PE_MODE,ADDRESS2(ptr_dpdiss_ave,globalx,len_elem,globaly,np*np),
					dpdiss_ave,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_dpdiss_biharmonic,globalx,len_elem,globaly,np*np),
					dpdiss_biharmonic,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_dp3d,globalx,len_elem,globaly,np*np),
					dp3d,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_dptens,globalx,len_dptens,globaly,np*np),
					dptens,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			while(get_reply!=4);
			/*caculate*/
			for(n=0;n<Nsub;n++){
				for(i=0;i<np;i++){
					for(j=0;j<np;j++){
						dpdiss_ave[n][i][j]=dpdiss_ave[n][i][j]+eta_ave_w*dp3d[n][i][j]/hypervis_subcycle;
						dpdiss_biharmonic[n][i][j]=dpdiss_biharmonic[n][i][j]+eta_ave_w*dptens[n][i][j]/hypervis_subcycle;
					}
				}
			}
			/*save data to host */
			put_reply=0;	
			athread_put(PE_MODE,dpdiss_ave,ADDRESS2(ptr_dpdiss_ave,globalx,len_elem,globaly,np*np),
					sizeof(double)*np*np*Nsub,&put_reply,0,0);
			athread_put(PE_MODE,dpdiss_biharmonic,ADDRESS2(ptr_dpdiss_biharmonic,globalx,len_elem,globaly,np*np),
					sizeof(double)*np*np*Nsub,&put_reply,0,0);
			while(put_reply!=2);
		}
	}
}

/*Kernel-D*/
typedef struct{
	/*pointer-0*/
	double*v,*T,*e_dp3d,*e_spheremp,*lap_t,*lap_v,*lap_dp;
	double nu_s,nu,dt;
	/*pointer-laplace_sphere_wk*/
	double* d_Dvv,*e_Dinv,rrearth,nu_p,nu_scale_top,nu_top;
	/*pointer-Vlaplace_sphere_wk*/
	double* e_metdet,*e_rmetdet,*e_D,*e_mp,*e_metinv;
	double *ttens,*dptens,*vtens;
	/*parameter-0*/
	int nets,nete,nlev,np,len_elem;
	/*parameter-Vlaplace_sphere_wk*/
	int len_ttens,len_dptens,len_vtens;
}str_d;

void slave_advance_hypervis_dp_kerneld_(str_d* dat)
{	
	/*receive  variable and pointer from Fortran*/
	double *ptr_v,*ptr_T,*ptr_e_dp3d,*ptr_e_spheremp,
				 *ptr_lap_t,*ptr_lap_v,*ptr_lap_dp;
	ptr_v=dat->v;ptr_T=dat->T;
	ptr_e_dp3d=dat->e_dp3d;
	ptr_e_spheremp=dat->e_spheremp;
	ptr_lap_t=dat->lap_t;
	ptr_lap_v=dat->lap_v;
	ptr_lap_dp=dat->lap_dp;
	double nu_s,nu,dt;
	nu_s=dat->nu_s;nu=dat->nu;dt=dat->dt;
	/*laplace_sphere_wk*/
	double* ptr_d_Dvv,*ptr_e_Dinv,rrearth;
	ptr_d_Dvv=dat->d_Dvv; ptr_e_Dinv=dat->e_Dinv;
	rrearth=dat->rrearth;
	/*Vlaplace_sphere_wk*/
	double* ptr_e_metdet,*ptr_e_rmetdet,*ptr_e_D,*ptr_e_mp,*ptr_e_metinv;
	ptr_e_metdet=dat->e_metdet;
	ptr_e_rmetdet=dat->e_rmetdet;
	ptr_e_D=dat->e_D;
	ptr_e_mp=dat->e_mp;
	ptr_e_metinv=dat->e_metinv;
	double* ptr_ttens,*ptr_dptens,*ptr_vtens;
	ptr_ttens=dat->ttens;
	ptr_dptens=dat->dptens;
	ptr_vtens=dat->vtens;
	double nu_p,nu_scale_top,nu_top;
	nu_p=dat->nu_p;
	nu_scale_top=dat->nu_scale_top;
	nu_top=dat->nu_top;
	int nets,nete,nlev,np;
	nets=dat->nets;
	nete=dat->nete;
	nlev=dat->nlev;
	np=dat->np;
	int len_elem,len_ttens,len_dptens,len_vtens;
	len_elem=dat->len_elem;
	len_ttens=dat->len_ttens;
	len_dptens=dat->len_dptens;
	len_vtens=dat->len_vtens;
	/*get threads id*/
	volatile unsigned int id=athread_get_id(-1);
	/*local parameter */
	volatile unsigned long get_reply,put_reply;
	int Nsub=nlev;
	int GridX=64;
	int LoopNumX=((nete-nets+1)+GridX-1)/GridX;
	int colid=id%GridX;
	int globalx,globaly,i,j,n,iloopx;
	/*main LDM*/
	double e_v[Nsub][2][np][np],e_T[Nsub][np][np],e_dp3d[Nsub][np][np],
				 e_spheremp[np][np],lap_v[2][np][np],lap_t[np][np],lap_dp[np][np];
	double d_Dvv[np][np],e_Dinv[2][2][np][np],e_metdet[np][np],e_rmetdet[np][np],
				 e_D[2][2][np][np],e_mp[np][np],e_metinv[2][2][np][np];
	double ttens[Nsub][np][np],dptens[Nsub][np][np],
				 vtens[Nsub][2][np][np];
	/*block loop*/
	for(iloopx=0;iloopx<LoopNumX;iloopx++){
		globalx=colid+iloopx*GridX;
		globaly=0;
		if(globalx<(nete-nets+1)){
			get_reply=0;
			athread_get(PE_MODE,ADDRESS2(ptr_v,globalx,len_elem,globaly,np*np*2),
					e_v,sizeof(double)*np*np*2*Nsub,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_T,globalx,len_elem,globaly,np*np),
					e_T,sizeof(double)*np*np*Nsub,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_dp3d,globalx,len_elem,globaly,np*np),
					e_dp3d,sizeof(double)*np*np*Nsub,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_spheremp,globalx,len_elem,globaly,0),
					e_spheremp,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ptr_lap_v,lap_v,sizeof(double)*np*np*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ptr_lap_t,lap_t,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ptr_lap_dp,lap_dp,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ptr_d_Dvv,d_Dvv,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_Dinv,globalx,len_elem,0,0),
					e_Dinv,sizeof(double)*np*np*2*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_metdet,globalx,len_elem,0,0),
					e_metdet,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_rmetdet,globalx,len_elem,0,0),
					e_rmetdet,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_D,globalx,len_elem,0,0),
					e_D,sizeof(double)*np*np*2*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_mp,globalx,len_elem,0,0),
					e_mp,sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_e_metinv,globalx,len_elem,0,0),
					e_metinv,sizeof(double)*np*np*2*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_ttens,globalx,len_ttens,globaly,np*np),
					ttens,sizeof(double)*np*np*Nsub,&get_reply,0,0,0); 
			athread_get(PE_MODE,ADDRESS2(ptr_dptens,globalx,len_dptens,globaly,np*np),
					dptens,sizeof(double)*np*np*Nsub,&get_reply,0,0,0); 
			athread_get(PE_MODE,ADDRESS2(ptr_vtens,globalx,len_vtens,globaly,np*np*2),
					vtens,sizeof(double)*np*np*2*Nsub,&get_reply,0,0,0); 
			while(get_reply!=17);

			for(n=0;n<Nsub;n++){
				nu_scale_top=1;
				if((globaly+n)==0)	nu_scale_top=4;
				if((globaly+n)==1)	nu_scale_top=2;
				if(nu_top>0&&(globaly+n)<3)
				{
					laplace_sphere_wk(np,e_T[n],d_Dvv,e_Dinv,e_spheremp,lap_t,rrearth);
					laplace_sphere_wk(np,e_dp3d[n],d_Dvv,e_Dinv,e_spheremp,lap_dp,rrearth);
					vlaplace_sphere_wk_ad_fortran_(&e_v[n][0][0][0],&lap_v[0][0][0],&d_Dvv[0][0],&e_metdet[0][0],
							&e_Dinv[0][0][0][0],&e_rmetdet[0][0],&rrearth,&e_D,&e_mp[0][0],&e_metinv[0][0][0][0],&e_spheremp[0][0]);
				}
				if(nu_top>0&&(globaly+n)<3)
				{
					for(j=0;j<np;j++)
						for(i=0;i<np;i++){
							ttens[n][j][i]=(-nu_s*ttens[n][j][i]+nu_scale_top*nu_top*lap_t[j][i]);
							dptens[n][j][i]=(-nu_p*dptens[n][j][i]+nu_scale_top*nu_top*lap_dp[j][i]);
							vtens[n][0][j][i]=(-nu*vtens[n][0][j][i]+nu_scale_top*nu_top*lap_v[0][j][i]);
							vtens[n][1][j][i]=(-nu*vtens[n][1][j][i]+nu_scale_top*nu_top*lap_v[1][j][i]);
						}
				}
				else{
					for(j=0;j<np;j++)
						for(i=0;i<np;i++){
							ttens[n][j][i]=-nu_s*ttens[n][j][i];
							dptens[n][j][i]=-nu_p*dptens[n][j][i];
							vtens[n][0][j][i]=-nu*vtens[n][0][j][i];
							vtens[n][1][j][i]=-nu*vtens[n][1][j][i];
						}
				}
			}
			/*save data to host*/
			put_reply=0;
			athread_put(PE_MODE,ttens,ADDRESS2(ptr_ttens,globalx,len_ttens,globaly,np*np),
					sizeof(double)*np*np*Nsub,&put_reply,0,0); 
			athread_put(PE_MODE,dptens,ADDRESS2(ptr_dptens,globalx,len_dptens,globaly,np*np),
					sizeof(double)*np*np*Nsub,&put_reply,0,0); 
			athread_put(PE_MODE,vtens,ADDRESS2(ptr_vtens,globalx,len_vtens,globaly,np*np*2),
					sizeof(double)*np*np*2*Nsub,&put_reply,0,0); 
			while(put_reply!=3);
		}
	}
}


/*Kernel-EF*/
typedef struct{
	double* vtens,*ttens,*dptens,*v,*T,*dp3d,*rspheremp;
	double cp,dt;
	int nets,nete,nlev,np,
			len_vtens,len_ttens,len_dptens,
			len_elem;
}str_ef;

void slave_advance_hypervis_dp_kernelef_(str_ef* dat){	
	/*receive  variable and pointer from Fortran*/
	double* ptr_vtens,*ptr_ttens,*ptr_dptens,*ptr_v,*ptr_T,*ptr_dp3d,*ptr_rspheremp;
	ptr_vtens=dat->vtens;ptr_ttens=dat->ttens;ptr_dptens=dat->dptens;
	ptr_rspheremp=dat->rspheremp;
	ptr_v=dat->v;ptr_T=dat->T;ptr_dp3d=dat->dp3d;
	double cp=dat->cp,dt=dat->dt;
	int nets,nete,nlev,np,
			len_vtens,len_ttens,len_dptens,
			len_elem;
	nets=dat->nets;
	nete=dat->nete;
	nlev=dat->nlev;
	np=dat->np;
	len_vtens=dat->len_vtens;
	len_ttens=dat->len_ttens;
	len_dptens=dat->len_dptens;
	len_elem=dat->len_elem;
	/*get threads id*/
	volatile unsigned int id=athread_get_id(-1);
	/*local parameter*/
	volatile unsigned long get_reply,put_reply;
	int GridX=64;
	int Nsub=nlev;
	int LoopNumX=((nete-nets+1)+GridX-1)/GridX;
	int colid=id%GridX;
	int globalx,globaly,i,j,n,iloopx;
	double heating;
	/*main LDM*/
	double vtens[Nsub][2][np][np],ttens[Nsub][np][np],dptens[Nsub][np][np],
				 rspheremp[np][np],v[Nsub][2][np][np],T[Nsub][np][np],
				 dp3d[Nsub][np][np];
	/*block loop*/
	for(iloopx=0;iloopx<LoopNumX;iloopx++){
		globalx=colid+iloopx*GridX;
		globaly=0;
		if(globalx<(nete-nets+1)){
			get_reply=0;
			athread_get(PE_MODE,ADDRESS2(ptr_vtens,globalx,len_vtens,globaly,np*np*2),
					vtens,Nsub*sizeof(double)*np*np*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_ttens,globalx,len_ttens,globaly,np*np),
					ttens,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_dptens,globalx,len_dptens,globaly,np*np),
					dptens,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_v,globalx,len_elem,globaly,np*np*2),
					v,Nsub*sizeof(double)*np*np*2,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_dp3d,globalx,len_elem,globaly,np*np),
					dp3d,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_T,globalx,len_elem,globaly,np*np),T
					,Nsub*sizeof(double)*np*np,&get_reply,0,0,0);
			athread_get(PE_MODE,ADDRESS2(ptr_rspheremp,globalx,len_elem,globaly,0),rspheremp,sizeof(double)*np*np,&get_reply,0,0,0);
			while(get_reply!=7);
			/*caculate*/
			for(n=0;n<Nsub;n++){
				for(i=0;i<np;i++){
					for(j=0;j<np;j++){
						vtens[n][0][i][j]=dt*vtens[n][0][i][j]*rspheremp[i][j];
						vtens[n][1][i][j]=dt*vtens[n][1][i][j]*rspheremp[i][j];
						ttens[n][i][j]=dt*ttens[n][i][j]*rspheremp[i][j];
						dptens[n][i][j]=dt*dptens[n][i][j]*rspheremp[i][j];
					}
				}
			}
			for(n=0;n<Nsub;n++){
				for(i=0;i<np;i++){
					for(j=0;j<np;j++){
						v[n][0][i][j]=v[n][0][i][j]+vtens[n][0][i][j];
						v[n][1][i][j]=v[n][1][i][j]+vtens[n][1][i][j];
						heating=vtens[n][0][i][j]*v[n][0][i][j]+vtens[n][1][i][j]*v[n][1][i][j];
						T[n][i][j]=T[n][i][j]+ttens[n][i][j]-heating/cp;
						dp3d[n][i][j]=dp3d[n][i][j]+dptens[n][i][j];
					}
				}
			}
			/*save data to host*/
			put_reply=0;	
			athread_put(PE_MODE,v,ADDRESS2(ptr_v,globalx,len_elem,globaly,np*np*2),
					sizeof(double)*np*np*2*Nsub,&put_reply,0,0);
			athread_put(PE_MODE,T,ADDRESS2(ptr_T,globalx,len_elem,globaly,np*np),
					sizeof(double)*np*np*Nsub,&put_reply,0,0);
			athread_put(PE_MODE,dp3d,ADDRESS2(ptr_dp3d,globalx,len_elem,globaly,np*np),
					sizeof(double)*np*np*Nsub,&put_reply,0,0);
			while(put_reply!=3);
		}
	}
}

/*Kernel-G-H*/
typedef struct {
	/*lap*/
	double*e_T,*e_dp3d,*ptens,*dptens,*e_variable_hyperviscosity,*e_tensorVisc,*e_Dinv,*d_Dvv,*e_spheremp;
	/*vlap*/
	double*e_v,*vtens,*e_vec_sphere2cart,*e_metdet,*e_rmetdet,*e_Dvv,*e_D,*e_mp,*e_metinv;
	double hypervis_power,hypervis_scaling,rrearth;
	/*lap*/
	double nu_ratio;
	int var_coef;
	int nets,nete,nlev,np,len_elem,len_ptens,len_dptens,len_vtens;
} str_g ;
/* Laplace */
void gradient_sphere_gh(int np ,double s[np][np],double d_Dvv[np][np],
		double Dinv[2][2][np][np],double ds[2][np][np],double rrearth){
	/*local variable*/	
	int i,j,l;
	double dsdx00,dsdy00,v1[np][np],v2[np][np];
	for(j=0;j<np;j++){
		for(l=0;l<np;l++){
			dsdx00=0.0;
			dsdy00=0.0;
			for(i=0;i<np;i++){
				dsdx00 = dsdx00 +d_Dvv[l][i]*s[j][i];
				dsdy00 = dsdy00 +d_Dvv[l][i]*s[i][j];
			}
			v1[j][l]= dsdx00*rrearth;
			v2[l][j]= dsdy00*rrearth;
		}
	}
	for(j=0;j<np;j++){
		for(i=0;i<np;i++){
			ds[0][j][i]=Dinv[0][0][j][i]*v1[j][i]+Dinv[0][1][j][i]*v2[j][i];
			ds[1][j][i]=Dinv[1][0][j][i]*v1[j][i]+Dinv[1][1][j][i]*v2[j][i];
		}
	}
}
void divergence_sphere_wk_gh(int np,double v[2][np][np],double d_Dvv[np][np],
		double e_Dinv[2][2][np][np],double e_rspheremp[np][np],double div[np][np],double rrearth){
	/*local variable*/
	double vtemp[2][np][np];
	int i,j,m,n;
	for(j=0;j<np;j++){
		for(i=0;i<np;i++)
		{
			vtemp[0][j][i]=(e_Dinv[0][0][j][i]*v[0][j][i] + e_Dinv[1][0][j][i]*v[1][j][i]);
			vtemp[1][j][i]=(e_Dinv[0][1][j][i]*v[0][j][i] + e_Dinv[1][1][j][i]*v[1][j][i]);
		}
	}

	for(n=0;n<np;n++)
		for(m=0;m<np;m++){
			div[n][m]=0;
			for(j=0;j<np;j++){
				div[n][m]=div[n][m]-(e_rspheremp[n][j]*vtemp[0][n][j]*d_Dvv[j][m]+
						e_rspheremp[j][m]*vtemp[1][j][m]*d_Dvv[j][n])*rrearth;
			}
		}
}
void laplace_sphere_wk_gh(int np,double s[np][np],double d_Dvv[np][np],double e_Dinv[2][2][np][np],double e_spheremp[np][np],double laplace[np][np],
		double e_variable_hyperviscosity[np][np],double e_tensorVisc[2][2][np][np],double rrearth,double hypervis_power,double hypervis_scaling,int var_coef){
	/*local variable*/
	double grads[2][np][np],oldgrads[2][np][np];
	int i,j;
	gradient_sphere_gh(np,s,d_Dvv,e_Dinv,grads,rrearth);
	if(var_coef){
		if(hypervis_power!=0){
			/*scalar viscosity with variable coefficient*/
			for(i=0;i<np;i++)
				for(j=0;j<np;j++){
					grads[0][i][j]=grads[0][i][j]*e_variable_hyperviscosity[i][j];
					grads[1][i][j]=grads[1][i][j]*e_variable_hyperviscosity[i][j];
				}
		}
		else if(hypervis_scaling!=0){
			for(i=0;i<np;i++)
				for(j=0;j<np;j++){
					oldgrads[0][i][j]=grads[0][i][j];
					oldgrads[1][i][j]=grads[1][i][j];
				}
			for(i=0;i<np;i++)
				for(j=0;j<np;j++){
					grads[0][i][j]=oldgrads[0][i][j]*e_tensorVisc[0][0][i][j]+
												 oldgrads[1][i][j]*e_tensorVisc[1][0][i][j];
					grads[1][i][j]=oldgrads[0][i][j]*e_tensorVisc[0][1][i][j]+
												 oldgrads[1][i][j]*e_tensorVisc[1][1][i][j];
				}

		}
	}
	divergence_sphere_wk_gh(np,grads,d_Dvv,e_Dinv,e_spheremp,laplace,rrearth);
}

typedef struct {
	double*e_rspheremp,*ptens,*dptens,*e_variable_hyperviscosity,*e_tensorVisc,*e_Dinv,*d_Dvv,*e_spheremp;
	double* vtens,*e_vec_sphere2cart,*e_metdet,	//vlap
				*e_rmetdet,*e_Dvv,*e_D,*e_mp,*e_metinv;
	double hypervis_power,hypervis_scaling,rrearth,nu_ratio2;
	int var_coef;
	int nets,nete,nlev,np,len_elem,len_ptens,len_dptens,len_vtens;
} str_h ;

/*kernel-G*/
void slave_advance_hypervis_dp_kernelg_(str_g* dat)
{
	/*get threads id*/
	volatile unsigned int id=athread_get_id(-1);
	volatile unsigned int colid,rowid;
	GET_ROW(rowid);GET_COL(colid);
	volatile unsigned long get_reply,put_reply;
	
	str_g dat0;
	if(id==0){
		get_reply=0;	
		athread_get(BCAST_MODE,dat,&dat0,sizeof(str_g),&get_reply,0xFF,0,0);
		while(get_reply!=1);
	}
	athread_syn(ARRAY_SCOPE,0xFFFF);
	
	//lap
	double* ptr_e_T,*ptr_e_dp3d,*ptr_ptens,*ptr_dptens,
		*ptr_e_variable_hyperviscosity,*ptr_e_tensorVisc,
		*ptr_e_Dinv,*ptr_d_Dvv,*ptr_e_spheremp;
	ptr_e_T=dat0.e_T;
	ptr_e_dp3d=dat0.e_dp3d;
	ptr_ptens=dat0.ptens;
	ptr_dptens=dat0.dptens;
	ptr_e_variable_hyperviscosity=dat0.e_variable_hyperviscosity;
	ptr_e_tensorVisc=dat0.e_tensorVisc;
	ptr_e_Dinv=dat0.e_Dinv;
	ptr_d_Dvv=dat0.d_Dvv;
	ptr_e_spheremp=dat0.e_spheremp;
	//vlap

	double*ptr_e_v,*ptr_vtens,*ptr_e_vec_sphere2cart,*ptr_e_metdet,
		*ptr_e_rmetdet,*ptr_e_Dvv,*ptr_e_D,*ptr_e_mp,*ptr_e_metinv;
	ptr_e_v       =dat0.e_v;
	ptr_vtens     =dat0.vtens ;
	ptr_e_vec_sphere2cart=dat0.e_vec_sphere2cart;
	ptr_e_metdet  =dat0.e_metdet ;
	ptr_e_rmetdet =dat0.e_rmetdet ;
	ptr_e_Dvv     =dat0.e_Dvv ;
	ptr_e_D       =dat0.e_D ;
	ptr_e_mp      =dat0.e_mp ;
	ptr_e_metinv =dat0.e_metinv;
	double hypervis_power,hypervis_scaling,rrearth;
	hypervis_power=dat0.hypervis_power;
	hypervis_scaling=dat0.hypervis_scaling;
	rrearth=dat0.rrearth;
	double nu_ratio;
	nu_ratio=dat0.nu_ratio;
	int var_coef;
	var_coef=dat0.var_coef;	
	int nets,nete,nlev,np,len_elem,len_ptens,len_dptens,len_vtens;
	nets=dat0.nets;
	nete=dat0.nete;
	nlev=dat0.nlev;
	np=dat0.np;
	len_elem=dat0.len_elem;
	len_ptens=dat0.len_ptens;
	len_dptens=dat0.len_dptens;
	len_vtens=dat0.len_vtens;

/*local parameter */
	int i,j,k,l,n,m;
	//NLEV
	double e_T[N*NP*NP],e_dp3d[N*NP*NP],
				 ptens[N*NP*NP],dptens[N*NP*NP];
	double e_v[N*2*NP*NP],vtens[N*2*NP*NP];
	//NO NLEV
	double e_variable_hyperviscosity[NP*NP],
				 e_tensorVisc[2*2*NP*NP],
				 e_Dinv[2*2*NP*NP],
				 e_spheremp[NP*NP];
	double e_vec_sphere2cart[2*3*NP*NP],
				 e_metdet[NP*NP],e_rmetdet[NP*NP],
				 e_D[2*2*NP*NP],e_mp[NP*NP],e_metinv[2*2*NP*NP];
	double	 d_Dvv[NP*NP];
	double tmpVal[23*NP*NP];
	doublev4 v4_e_T,v4_e_dp3d,v4_ptens,v4_dptens,v4_e_v,v4_vtens;
	/*BroadCast Global Variabels for each cpes*/
	if(id==0){
		get_reply=0;	
		athread_get(BCAST_MODE,ptr_d_Dvv,d_Dvv,sizeof(double)*np*np,&get_reply,0xFF,0,0);
		while(get_reply!=1);
	}
	for(l=0;l<(nete-nets+1+1)/2;l++){
		/*colid=0 Load variables that contains nlev */
		if(colid==0){
			if(rowid<4){//upper part cpes (0->30)
				//Laplace variables
				get_reply=0; 
				athread_get(PE_MODE,ptr_e_T+len_elem*2*l+rowid*N*np*np,
						e_T,sizeof(double)*np*np*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_e_dp3d+len_elem*2*l+rowid*N*np*np,
						e_dp3d,sizeof(double)*np*np*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_ptens+len_ptens*2*l+rowid*N*np*np,
						ptens,sizeof(double)*np*np*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_dptens+len_dptens*2*l+rowid*N*np*np,
						dptens,sizeof(double)*np*np*N,&get_reply,0,0,0);
				//Laplace variabels
				athread_get(PE_MODE,ptr_e_v+len_elem*2*l+rowid*N*np*np*2,
						e_v,sizeof(double)*np*np*2*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_vtens+len_vtens*2*l+rowid*N*np*np*2,
						vtens,sizeof(double)*np*np*2*N,&get_reply,0,0,0); 
				while(get_reply!=6);
			}
			else if((2*l+1)<(nete-nets+1)){ //lower part cpes (32->62)
				//Laplace variables
				get_reply=0; 
				athread_get(PE_MODE,ptr_e_T+len_elem*(2*l+1)+(rowid-4)*N*np*np,
						e_T,sizeof(double)*np*np*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_e_dp3d+len_elem*(2*l+1)+(rowid-4)*N*np*np,
						e_dp3d,sizeof(double)*np*np*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_ptens+len_ptens*(2*l+1)+(rowid-4)*N*np*np,
						ptens,sizeof(double)*np*np*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_dptens+len_dptens*(2*l+1)+(rowid-4)*N*np*np,
						dptens,sizeof(double)*np*np*N,&get_reply,0,0,0);
				//Laplace variabels
				athread_get(PE_MODE,ptr_e_v+len_elem*(2*l+1)+(rowid-4)*N*np*np*2,
						e_v,sizeof(double)*np*np*2*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_vtens+len_vtens*(2*l+1)+(rowid-4)*N*np*np*2,
						vtens,sizeof(double)*np*np*2*N,&get_reply,0,0,0); 
				while(get_reply!=6);
			}
		}
		/*Scatter data according to rowid by register communication*/
		for(n=1;n<8;n++){
			if(colid==0){
				for(m=0;m<np;m++){
					simd_load(v4_e_T,e_T+n*np*np+m*np);
					REG_PUTR(v4_e_T,n);	
					simd_load(v4_e_dp3d,e_dp3d+n*np*np+m*np);
					REG_PUTR(v4_e_dp3d,n);
					simd_load(v4_ptens,ptens+n*np*np+m*np);
					REG_PUTR(v4_ptens,n);
					simd_load(v4_dptens,dptens+n*np*np+m*np);
					REG_PUTR(v4_dptens,n);
					simd_load(v4_e_v,e_v+n*(np*np*2)+2*m*np);
					REG_PUTR(v4_e_v,n);
					simd_load(v4_e_v,e_v+n*(np*np*2)+(2*m+1)*np);
					REG_PUTR(v4_e_v,n);
					simd_load(v4_vtens,vtens+n*(np*np*2)+2*m*np);
					REG_PUTR(v4_vtens,n);
					simd_load(v4_vtens,vtens+n*(np*np*2)+(2*m+1)*np);
					REG_PUTR(v4_vtens,n);
				}
			}
			if(colid==n){
				for(m=0;m<np;m++){
					REG_GETR(v4_e_T);
					simd_store(v4_e_T,e_T+m*np);
					REG_GETR(v4_e_dp3d);
					simd_store(v4_e_dp3d,e_dp3d+m*np);
					REG_GETR(v4_ptens);
					simd_store(v4_ptens,ptens+m*np);
					REG_GETR(v4_dptens);
					simd_store(v4_dptens,dptens+m*np);
					REG_GETR(v4_e_v);
					simd_store(v4_e_v,e_v+2*m*np);
					REG_GETR(v4_e_v);
					simd_store(v4_e_v,e_v+(2*m+1)*np);
					REG_GETR(v4_vtens);
					simd_store(v4_vtens,vtens+2*m*np);
					REG_GETR(v4_vtens);
					simd_store(v4_vtens,vtens+(2*m+1)*np);
				}
			}
		//athread_syn(ROW_SCOPE,0xFF);
		}
		/*id=0 && id=32(if necessary) Broadcast variabels that don't contains nlev  */
		if(id==0){
			get_reply=0;
			athread_get(BCAST_MODE,ptr_e_metinv+(2*l)*len_elem,tmpVal,sizeof(double)*np*np*23,&get_reply,0x0F,0,0);
			athread_get(BCAST_MODE,ptr_e_variable_hyperviscosity+2*l*len_elem,e_variable_hyperviscosity,sizeof(double)*np*np,&get_reply,0x0F,0,0);
			athread_get(BCAST_MODE,ptr_e_tensorVisc+2*l*len_elem,e_tensorVisc,sizeof(double)*np*np*2*2,&get_reply,0x0F,0,0);

			//athread_get(BCAST_MODE,ptr_e_Dinv+(2*l)*len_elem,e_Dinv,sizeof(double)*np*np*2*2,&get_reply,0x0F,0,0);
			//athread_get(BCAST_MODE,ptr_e_spheremp+(2*l)*len_elem,e_spheremp,sizeof(double)*np*np,&get_reply,0x0F,0,0);
			///*vlaplace*/
			//athread_get(BCAST_MODE,ptr_e_vec_sphere2cart+(2*l)*len_elem,e_vec_sphere2cart,sizeof(double)*np*np*2*3,&get_reply,0x0F,0,0);
			//athread_get(BCAST_MODE,ptr_e_metdet+(2*l)*len_elem,e_metdet,sizeof(double)*np*np,&get_reply,0x0F,0,0);
			//athread_get(BCAST_MODE,ptr_e_rmetdet+(2*l)*len_elem,e_rmetdet,sizeof(double)*np*np,&get_reply,0x0F,0,0);
			//athread_get(BCAST_MODE,ptr_e_D+(2*l)*len_elem,e_D,sizeof(double)*np*np*2*2,&get_reply,0x0F,0,0);
			//athread_get(BCAST_MODE,ptr_e_mp+(2*l)*len_elem,e_mp,sizeof(double)*np*np,&get_reply,0x0F,0,0);
			//athread_get(BCAST_MODE,ptr_e_metinv+(2*l)*len_elem,e_metinv,sizeof(double)*np*np*2*2,&get_reply,0x0F,0,0);

			while(get_reply!=3);	
		}
		if(id==32&&(2*l+1)<(nete-nets+1)){
			get_reply=0;
			athread_get(BCAST_MODE,ptr_e_metinv+(2*l+1)*len_elem,tmpVal,sizeof(double)*np*np*23,&get_reply,0xF0,0,0);
			athread_get(BCAST_MODE,ptr_e_variable_hyperviscosity+(2*l+1)*len_elem,e_variable_hyperviscosity,sizeof(double)*np*np,&get_reply,0xF0,0,0);
			athread_get(BCAST_MODE,ptr_e_tensorVisc+(2*l+1)*len_elem,e_tensorVisc,sizeof(double)*np*np*2*2,&get_reply,0xF0,0,0);

			//athread_get(BCAST_MODE,ptr_e_Dinv+(2*l+1)*len_elem,e_Dinv,sizeof(double)*np*np*2*2,&get_reply,0xF0,0,0);
			//athread_get(BCAST_MODE,ptr_e_spheremp+(2*l+1)*len_elem,e_spheremp,sizeof(double)*np*np,&get_reply,0xF0,0,0);
			///*vlaplace*/
			//athread_get(BCAST_MODE,ptr_e_vec_sphere2cart+(2*l+1)*len_elem,e_vec_sphere2cart,sizeof(double)*np*np*2*3,&get_reply,0xF0,0,0);
			//athread_get(BCAST_MODE,ptr_e_metdet+(2*l+1)*len_elem,e_metdet,sizeof(double)*np*np,&get_reply,0xF0,0,0);
			//athread_get(BCAST_MODE,ptr_e_rmetdet+(2*l+1)*len_elem,e_rmetdet,sizeof(double)*np*np,&get_reply,0xF0,0,0);
			//athread_get(BCAST_MODE,ptr_e_D+(2*l+1)*len_elem,e_D,sizeof(double)*np*np*2*2,&get_reply,0xF0,0,0);
			//athread_get(BCAST_MODE,ptr_e_mp+(2*l+1)*len_elem,e_mp,sizeof(double)*np*np,&get_reply,0xF0,0,0);
			//athread_get(BCAST_MODE,ptr_e_metinv+(2*l+1)*len_elem,e_metinv,sizeof(double)*np*np*2*2,&get_reply,0xF0,0,0);
	
			while(get_reply!=3);	
		}
		athread_syn(ARRAY_SCOPE,0xFFFF);
	//#define NEW
#ifdef  NEW
		laplace_sphere_wk_gh(np,e_T,d_Dvv,e_Dinv,e_spheremp,ptens,e_variable_hyperviscosity,e_tensorVisc,rrearth,hypervis_power,hypervis_scaling,var_coef);
		laplace_sphere_wk_gh(np,e_dp3d,d_Dvv,e_Dinv,e_spheremp,dptens,e_variable_hyperviscosity	,e_tensorVisc,rrearth,hypervis_power,hypervis_scaling,var_coef);
		vlaplace_sphere_wk_gh_fortran_(&e_v[0],&vtens[0],
				&e_vec_sphere2cart[0],&hypervis_power,
				&hypervis_scaling,&var_coef,&nu_ratio,&rrearth,
				/*part-lap*/
				&e_variable_hyperviscosity[0],&e_tensorVisc[0],
				&d_Dvv[0],&e_Dinv[0],&e_spheremp[0],
				/*part-Vlap*/
				&e_metdet[0],&e_rmetdet[0],
				&e_D[0],&e_mp[0],&e_metinv[0]);
#else
		laplace_sphere_wk_gh(np,e_T,d_Dvv,tmpVal+10*np*np,tmpVal+22*np*np,ptens,e_variable_hyperviscosity,e_tensorVisc,rrearth,hypervis_power,hypervis_scaling,var_coef);
		laplace_sphere_wk_gh(np,e_dp3d,d_Dvv,tmpVal+10*np*np,tmpVal+22*np*np,dptens,e_variable_hyperviscosity	,e_tensorVisc,rrearth,hypervis_power,hypervis_scaling,var_coef);
		vlaplace_sphere_wk_gh_fortran_(&e_v[0],&vtens[0],tmpVal+14*np*np,&hypervis_power,&hypervis_scaling,&var_coef,&nu_ratio,&rrearth,
				//part-lap
				&e_variable_hyperviscosity[0],&e_tensorVisc[0],&d_Dvv[0],tmpVal+10*np*np,tmpVal+22*np*np,
				//part-Vlap
				tmpVal+4*np*np,tmpVal+5*np*np,tmpVal+6*np*np,tmpVal+20*np*np,tmpVal);
#endif
		if(rowid<4){
			if(id<30){
				put_reply=0;
				athread_put(PE_MODE,ptens,ptr_ptens+len_ptens*2*l+id*np*np,
						sizeof(double)*np*np,&put_reply,0,0); 
				athread_put(PE_MODE,dptens,ptr_dptens+len_dptens*2*l+id*np*np,
						sizeof(double)*np*np,&put_reply,0,0); 
				athread_put(PE_MODE,vtens,ptr_vtens+len_vtens*2*l+id*np*np*2,
						sizeof(double)*np*np*2,&put_reply,0,0); 
				while(put_reply!=3);
			}
		}
		else if((2*l+1)<(nete-nets+1)){
			if(id<62){
				put_reply=0;
				athread_put(PE_MODE,ptens,ptr_ptens+len_ptens*(2*l+1)+(id-32)*np*np,
						sizeof(double)*np*np,&put_reply,0,0); 
				athread_put(PE_MODE,dptens,ptr_dptens+len_dptens*(2*l+1)+(id-32)*np*np,
						sizeof(double)*np*np,&put_reply,0,0); 
				athread_put(PE_MODE,vtens,ptr_vtens+len_vtens*(2*l+1)+(id-32)*np*np*2,
						sizeof(double)*np*np*2,&put_reply,0,0); 
				while(put_reply!=3);
			}
		}
	}

}

/*kernel-H*/
void slave_advance_hypervis_dp_kernelh_(str_h* dat){

	/*get threads id*/
	volatile unsigned int id=athread_get_id(-1);
	volatile unsigned long get_reply,put_reply;
	str_h dat0;	
	if(id==0){
		get_reply=0;	
		athread_get(BCAST_MODE,dat,&dat0,sizeof(str_h),&get_reply,0xFF,0,0);
		while(get_reply!=1);
	}
	athread_syn(ARRAY_SCOPE,0xFFFF);
	double* ptr_e_rspheremp,*ptr_ptens,*ptr_dptens,
		*ptr_e_variable_hyperviscosity,*ptr_e_tensorVisc,
		*ptr_e_Dinv,*ptr_d_Dvv,*ptr_e_spheremp;
	ptr_e_rspheremp=dat0.e_rspheremp;
	ptr_ptens=dat0.ptens;
	ptr_dptens=dat0.dptens;
	ptr_e_variable_hyperviscosity=dat0.e_variable_hyperviscosity;
	ptr_e_tensorVisc=dat0.e_tensorVisc;
	ptr_e_Dinv=dat0.e_Dinv;
	ptr_d_Dvv=dat0.d_Dvv;
	ptr_e_spheremp=dat0.e_spheremp;

	double* ptr_vtens,*ptr_e_vec_sphere2cart,*ptr_e_metdet,
		*ptr_e_rmetdet,*ptr_e_Dvv,*ptr_e_D,*ptr_e_mp,*ptr_e_metinv;
	ptr_vtens     =dat0.vtens ;
	ptr_e_vec_sphere2cart=dat0.e_vec_sphere2cart;
	ptr_e_metdet  =dat0.e_metdet ;
	ptr_e_rmetdet =dat0.e_rmetdet ;
	ptr_e_Dvv     =dat0.e_Dvv ;
	ptr_e_D       =dat0.e_D ;
	ptr_e_mp      =dat0.e_mp ;
	ptr_e_metinv =dat0.e_metinv;

	double hypervis_power,hypervis_scaling,rrearth;
	hypervis_power=dat0.hypervis_power;
	hypervis_scaling=dat0.hypervis_scaling;
	rrearth=dat0.rrearth;
	double nu_ratio2=dat0.nu_ratio2;
	int var_coef;
	var_coef=dat0.var_coef;	
	int nets,nete,nlev,np,len_elem,len_ptens,len_dptens,len_vtens;
	nets=dat0.nets;
	nete=dat0.nete;
	nlev=dat0.nlev;
	np=dat0.np;
	len_elem=dat0.len_elem;
	len_ptens=dat0.len_ptens;
	len_dptens=dat0.len_dptens;
	len_vtens=dat0.len_vtens;

/*local parameter */
	int i,j,k,l,m,n;
	volatile unsigned int rowid;GET_ROW(rowid);
	volatile unsigned int colid;GET_COL(colid);
	//Laplace
	/*The variables contains Nlev -- Using REGISTER COMMUNICATION to share data */ 
	double ptens[N*NP*NP],dptens[N*NP*NP];
	/*Global variables -- Using BROADCAST to share data*/
	double d_Dvv[NP*NP];
	//Vlaplace
	/*The variables contains Nlev -- Using REGISTER COMMUNICATION to share data */ 
	double vtens[N*2*NP*NP];
	/*The variables don't contains Nlev -- Using ROW BROADCAST to share data*/
	double tmpVal[24*NP*NP];
	double e_variable_hyperviscosity[NP*NP],e_tensorVisc[2*2*NP*NP];
	/*local variables*/
	double tmp[NP*NP];
	double v[2*NP*NP];	
	doublev4 v4_ptens,v4_dptens,v4_vtens;
	/*BroadCast Global Variabels for each cpes*/
	if(id==0){
		get_reply=0;	
		athread_get(BCAST_MODE,ptr_d_Dvv,d_Dvv,sizeof(double)*np*np,&get_reply,0xFF,0,0);
		while(get_reply!=1);
	}
	for(l=0;l<(nete-nets+1+1)/2;l++){
		/*colid=0 Load variables that contains nlev */
		if(colid==0){
			if(rowid<4){//upper part cpes (0->30)
				//Laplace variables
				get_reply=0; 
				athread_get(PE_MODE,ptr_ptens+len_ptens*2*l+rowid*N*np*np,
						ptens,sizeof(double)*np*np*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_dptens+len_dptens*2*l+rowid*N*np*np,
						dptens,sizeof(double)*np*np*N,&get_reply,0,0,0);
				//Laplace variables
				athread_get(PE_MODE,ptr_vtens+len_vtens*2*l+rowid*N*np*np*2,
						vtens,sizeof(double)*np*np*2*N,&get_reply,0,0,0); 
				while(get_reply!=3);
			}
			else if((2*l+1)<(nete-nets+1)){ //lower part cpes (32->62)
				//Laplace variables
				get_reply=0; 
				athread_get(PE_MODE,ptr_ptens+len_ptens*(2*l+1)+(rowid-4)*N*np*np,
						ptens,sizeof(double)*np*np*N,&get_reply,0,0,0);
				athread_get(PE_MODE,ptr_dptens+len_dptens*(2*l+1)+(rowid-4)*N*np*np,
						dptens,sizeof(double)*np*np*N,&get_reply,0,0,0);
				//VLaplace variables
				athread_get(PE_MODE,ptr_vtens+len_vtens*(2*l+1)+(rowid-4)*N*np*np*2,
						vtens,sizeof(double)*np*np*2*N,&get_reply,0,0,0); 
				while(get_reply!=3);
			}
		}
		/*Scatter data according to rowid by register communication*/
		for(n=1;n<8;n++){
			if(colid==0){
				for(m=0;m<np;m++){
					simd_load(v4_ptens,ptens+n*np*np+m*np);
					REG_PUTR(v4_ptens,n);

					simd_load(v4_dptens,dptens+n*np*np+m*np);
					REG_PUTR(v4_dptens,n);

					simd_load(v4_vtens,vtens+n*(np*np*2)+2*m*np);
					REG_PUTR(v4_vtens,n);
					simd_load(v4_vtens,vtens+n*(np*np*2)+(2*m+1)*np);
					REG_PUTR(v4_vtens,n);
				}
			}
			if(colid==n){
				for(m=0;m<np;m++){
					REG_GETR(v4_ptens);
					simd_store(v4_ptens,ptens+m*np);

					REG_GETR(v4_dptens);
					simd_store(v4_dptens,dptens+m*np);

					REG_GETR(v4_vtens);
					simd_store(v4_vtens,vtens+2*m*np);
					REG_GETR(v4_vtens);
					simd_store(v4_vtens,vtens+(2*m+1)*np);
				}
			}
			//athread_syn(ROW_SCOPE,0xFF);
		}
		/*id=0 && id=32(if necessary) Broadcast variabels that don't contains nlev  */
		if(id==0){
			get_reply=0;
			athread_get(BCAST_MODE,ptr_e_metinv                 +(2*l)*len_elem,tmpVal,                   sizeof(double)*np*np*24,&get_reply,0x0F,0,0);
			athread_get(BCAST_MODE,ptr_e_variable_hyperviscosity+(2*l)*len_elem,e_variable_hyperviscosity,sizeof(double)*np*np,&get_reply,0x0F,0,0);
			athread_get(BCAST_MODE,ptr_e_tensorVisc             +(2*l)*len_elem,e_tensorVisc,             sizeof(double)*np*np*2*2,&get_reply,0x0F,0,0);
			while(get_reply!=3);	
		}
		if(id==32&&(2*l+1)<(nete-nets+1)){
			get_reply=0;
			athread_get(BCAST_MODE,ptr_e_metinv                 +(2*l+1)*len_elem,tmpVal,                   sizeof(double)*np*np*24,&get_reply,0xF0,0,0);
			athread_get(BCAST_MODE,ptr_e_variable_hyperviscosity+(2*l+1)*len_elem,e_variable_hyperviscosity,sizeof(double)*np*np,&get_reply,0xF0,0,0);
			athread_get(BCAST_MODE,ptr_e_tensorVisc             +(2*l+1)*len_elem,e_tensorVisc,             sizeof(double)*np*np*2*2,&get_reply,0xF0,0,0);
			while(get_reply!=3);	
		}	
		athread_syn(ARRAY_SCOPE,0xFFFF);
		/*Calculate ptens,dptens and vtens , and write them to host*/
		doublev4 tmp1,tmp2,tmp3;
		for(i=0;i<np;i++){
			simd_load(tmp1,tmpVal+23*np*np+i*np);
			simd_load(tmp2,ptens+i*np);
			tmp3=tmp1*tmp2;
			simd_store(tmp3,tmp+i*np);
		//	for(j=0;j<np;j++){
		//		tmp[i*np+j]=tmpVal[23*np*np+i*np+j]*ptens[i*np+j];
		//	}
		}
		laplace_sphere_wk_gh(np,tmp,d_Dvv,tmpVal+10*np*np,tmpVal+22*np*np,ptens,e_variable_hyperviscosity,
				e_tensorVisc,rrearth,hypervis_power,hypervis_scaling,var_coef);
		for(i=0;i<np;i++){
			simd_load(tmp1,tmpVal+23*np*np+i*np);
			simd_load(tmp2,dptens+i*np);
			tmp3=tmp1*tmp2;
			simd_store(tmp3,tmp+i*np);
	
		//	for(j=0;j<np;j++){
		//		tmp[i*np+j]=tmpVal[23*np*np+i*np+j]*dptens[i*np+j];
		//	}
		}
		laplace_sphere_wk_gh(np,tmp,d_Dvv,tmpVal+10*np*np,tmpVal+22*np*np,dptens,e_variable_hyperviscosity,
				e_tensorVisc,rrearth,hypervis_power,hypervis_scaling,var_coef);
		for(i=0;i<np;i++){
			simd_load(tmp1,tmpVal+23*np*np+i*np);
			simd_load(tmp2,vtens+i*np);
			tmp3=tmp1*tmp2;
			simd_store(tmp3,v+i*np);
			simd_load(tmp2,vtens+1*np*np+i*np);
			tmp3=tmp1*tmp2;
			simd_store(tmp3,v+1*np*np+i*np);
		//	for(j=0;j<np;j++){
		//		v[0*np*np+i*np+j]=tmpVal[23*np*np+i*np+j]*vtens[0*np*np+i*np+j];
		//		v[1*np*np+i*np+j]=tmpVal[23*np*np+i*np+j]*vtens[1*np*np+i*np+j];
		//	}
		}
		vlaplace_sphere_wk_gh_fortran_(v,vtens,tmpVal+14*np*np,&hypervis_power,
				&hypervis_scaling,&var_coef,&nu_ratio2,&rrearth,
				//part-lap
				&e_variable_hyperviscosity[0],&e_tensorVisc[0],d_Dvv,tmpVal+10*np*np,tmpVal+22*np*np,
				//part-Vlap
				tmpVal+4*np*np,tmpVal+5*np*np,tmpVal+6*np*np,tmpVal+20*np*np,tmpVal);

		if(rowid<4){
			if(id<30){
				put_reply=0;
				athread_put(PE_MODE,ptens,ptr_ptens+len_ptens*2*l+id*np*np,
						sizeof(double)*np*np,&put_reply,0,0); 
				athread_put(PE_MODE,dptens,ptr_dptens+len_dptens*2*l+id*np*np,
						sizeof(double)*np*np,&put_reply,0,0); 
				athread_put(PE_MODE,vtens,ptr_vtens+len_vtens*2*l+id*np*np*2,
						sizeof(double)*np*np*2,&put_reply,0,0); 
				while(put_reply!=3);
			}
		}
		if(rowid>3&&(2*l+1)<(nete-nets+1)){
			if(id<62){
				put_reply=0;
				athread_put(PE_MODE,ptens,ptr_ptens+len_ptens*(2*l+1)+(id-32)*np*np,
						sizeof(double)*np*np,&put_reply,0,0); 
				athread_put(PE_MODE,dptens,ptr_dptens+len_dptens*(2*l+1)+(id-32)*np*np,
						sizeof(double)*np*np,&put_reply,0,0); 
				athread_put(PE_MODE,vtens,ptr_vtens+len_vtens*(2*l+1)+(id-32)*np*np*2,
						sizeof(double)*np*np*2,&put_reply,0,0); 
				while(put_reply!=3);
			}
		}
	}
}

