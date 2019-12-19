/*************************************************************************
    > File Name: subFunForSlave.h
    > Author: wusihai
    > Mail: wusihai18@gmail.com 
    > Created Time: Fri Dec 21 16:38:59 2018
 ************************************************************************/

/*constants parameter*/
#define NLEV      30
#define NP        4
#define N         2
#define L         8

#include <slave.h>
#include <simd.h>
#include "register.h"
#include "dma_macros.h"

/*Struct of slave */
typedef struct {
	double *ptr_e_v,*ptr_e_ps_v,*ptr_e_vn0,*ptr_e_Dinv,*ptr_e_metdet, 
				 *ptr_e_rmetdet,*ptr_e_D,*ptr_e_dp3d,*ptr_e_T,*ptr_e_Qdp, 
				 *ptr_e_phis,*ptr_e_phi,*ptr_e_eta_dot_dpdn,*ptr_e_omega_p, 
				 *ptr_e_spheremp,
				 *ptr_e_fcor,*ptr_e_pecnd,
				 *ptr_e_v_np1,*ptr_e_T_np1,*ptr_e_dp3d_np1,*ptr_e_ps_v_np1, 
				 *ptr_e_v_nm1,*ptr_e_T_nm1,*ptr_e_dp3d_nm1,*ptr_e_ps_v_nm1, 
				 *ptr_d_Dvv,*ptr_hvcoord;
	double eta_ave_w,rrearth,dt2;
	double Cp,Rgas,Cpwater_vapor,Rwater_vapor,kappa;
	int len_elem;
	int nets,nete,qn0,use_cpstar,qsplit,rsplit,n0,np1,nm1;
} str_caar1;


/*Declaration of external fortran subroutine*/
extern void gradient_sphere_fortran_(double* s,double*Dvv,double*e_Dinv,
		double*ds,double* rrearth);
extern void divergence_sphere_fortran_(double* v,double* div,double* d_Dvv,
		double* e_metdet,double* e_Dinv,double* e_rmetdet,double* rrearth);
extern void vorticity_sphere_fortran_(double* v,double* e_D,double* d_Dvv,
		double* e_rmetdet,double* vort,double* rrearth);
