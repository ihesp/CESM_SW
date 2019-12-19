subroutine vlaplace_sphere_wk_ad_fortran(v,laplace, &
d_Dvv,e_metdet,e_Dinv, &
e_rmetdet,rrearth,e_D, &
e_mp,e_metinv,e_spheremp) 
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   logic:
!      tensorHV:     requires cartesian
!      nu_div/=nu:   requires contra formulatino
!
!   One combination NOT supported:  tensorHV and nu_div/=nu then abort
!
implicit none
integer,parameter  :: np=4,real_kind=8
real(kind=8),dimension(np,np,2)	:: v,laplace
real(kind=8),dimension(np,np)	::   d_Dvv,e_metdet
real(kind=8),dimension(np,np,2,2)	:: e_Dinv
real(kind=8),dimension(np,np)	::   e_rmetdet
real(kind=8)	::   rrearth
real(kind=8),dimension(np,np,2,2)	::   e_D
real(kind=8),dimension(np,np)	::   e_mp
real(kind=8),dimension(np,np,2,2)	::  e_metinv
real(kind=8),dimension(np,np)	::   e_spheremp

!-------------var_coef=.false.-----------------!
!    if (hypervis_scaling/=0 .and. var_coef) then
!       ! tensorHV is turned on - requires cartesian formulation
!       if (present(nu_ratio)) then
!          if (nu_ratio /= 1) then
!             call abortmp('ERROR: tensorHV can not be used with nu_div/=nu')
!          endif
!       endif
!       call vlaplace_sphere_wk_cartesian(v,deriv,elem,laplace,var_coef)
!    else  
!       ! all other cases, use contra formulation:
!       call vlaplace_sphere_wk_contra(v,deriv,elem,laplace,var_coef,nu_ratio)
!    endif
!----------------------------------------------!

!       ! all other cases, use contra formulation:
call vlaplace_sphere_wk_contra1(v,laplace,d_Dvv,e_metdet,e_Dinv,e_rmetdet,rrearth,e_D,e_mp,e_metinv,e_spheremp)
end subroutine vlaplace_sphere_wk_ad_fortran

subroutine vlaplace_sphere_wk_contra1(v,laplace, &
d_Dvv,e_metdet,e_Dinv, &
e_rmetdet,rrearth,e_D, &
e_mp,e_metinv,e_spheremp) 
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   du/dt = laplace(u) = grad(div) - curl(vor)
!   < PHI du/dt > = < PHI laplace(u) >        PHI = covariant, u = contravariant
!                 = < PHI grad(div) >  - < PHI curl(vor) >
!                 = grad_wk(div) - curl_wk(vor)               
!
integer,parameter  :: np=4,real_kind=8
real(kind=8),dimension(np,np,2)	:: v,laplace
real(kind=8),dimension(np,np)	::   d_Dvv,e_metdet
real(kind=8),dimension(np,np,2,2)	:: e_Dinv
real(kind=8),dimension(np,np)	::   e_rmetdet
real(kind=8)	::   rrearth
real(kind=8),dimension(np,np,2,2)	::   e_D
real(kind=8),dimension(np,np)	::   e_mp
real(kind=8),dimension(np,np,2,2)	::  e_metinv
real(kind=8),dimension(np,np)	::   e_spheremp

!local 
real(kind=real_kind) :: lap_tmp(np,np,2)
real(kind=real_kind) :: lap_tmp2(np,np,2)
integer i,j,l,m,n
real(kind=real_kind) :: vor(np,np),div(np,np)
real(kind=real_kind) :: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y

call divergence_sphere1(v,d_Dvv,e_metdet,e_Dinv,e_rmetdet,rrearth,div)
call vorticity_sphere1(v,e_D,d_Dvv,e_rmetdet,vor,rrearth)

!--------var_coef=.false. && nu_ratio is not present-----!
!if (var_coef .and. hypervis_power/=0 ) then
!	! scalar viscosity with variable coefficient
!	div = div*elem%variable_hyperviscosity(:,:)
!	vor = vor*elem%variable_hyperviscosity(:,:)
!endif
!
!if (present(nu_ratio)) div = nu_ratio*div
!---------------------------------------------------------!

call gradient_sphere_wk_testcov1(div,lap_tmp,e_mp,e_metinv,d_Dvv,e_metdet,e_D,rrearth)
call curl_sphere_wk_testcov1(vor,lap_tmp2,d_Dvv,e_D,e_mp,rrearth)
laplace = lap_tmp - lap_tmp2

do n=1,np
	do m=1,np
		! add in correction so we dont damp rigid rotation
#define UNDAMPRR
#ifdef UNDAMPRR
		laplace(m,n,1)=laplace(m,n,1) + 2*e_spheremp(m,n)*v(m,n,1)*(rrearth**2)
		laplace(m,n,2)=laplace(m,n,2) + 2*e_spheremp(m,n)*v(m,n,2)*(rrearth**2)
#endif
	enddo
enddo
end subroutine 


subroutine divergence_sphere1(v,d_Dvv,e_metdet,e_Dinv,e_rmetdet,rrearth,div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
!
integer,parameter  :: np=4,real_kind=8
real(kind=8), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
real(kind=8),dimension(np,np)	::   d_Dvv,e_metdet
real(kind=8),dimension(np,np,2,2)	:: e_Dinv
real(kind=8),dimension(np,np)	::   e_rmetdet
real(kind=8)	::   rrearth
real(kind=real_kind) :: div(np,np)

integer i
integer j
integer l

real(kind=real_kind) ::  dudx00
real(kind=real_kind) ::  dvdy00
real(kind=real_kind) ::  gv(np,np,2),vvtemp(np,np)

! convert to contra variant form and multiply by g
!OMP_COLLAPSE_SIMD
!DIR_VECTOR_ALIGNED
do j=1,np
	do i=1,np
		gv(i,j,1)=e_metdet(i,j)*(e_Dinv(i,j,1,1)*v(i,j,1) + e_Dinv(i,j,1,2)*v(i,j,2))
		gv(i,j,2)=e_metdet(i,j)*(e_Dinv(i,j,2,1)*v(i,j,1) + e_Dinv(i,j,2,2)*v(i,j,2))
	enddo
enddo

! compute d/dx and d/dy         
do j=1,np
	do l=1,np
		dudx00=0.0d0
		dvdy00=0.0d0
		!DIR$ UNROLL(NP)
		do i=1,np
			dudx00 = dudx00 + d_Dvv(i,l)*gv(i,j ,1)
			dvdy00 = dvdy00 + d_Dvv(i,l)*gv(j,i,2)
		end do
		div(l,j) = dudx00
		vvtemp(j,l) = dvdy00
	end do
end do

!OMP_COLLAPSE_SIMD
!DIR_VECTOR_ALIGNED
do j=1,np
	do i=1,np
		div(i,j)=(div(i,j)+vvtemp(i,j))*(e_rmetdet(i,j)*rrearth)
	end do
end do

end subroutine 


subroutine vorticity_sphere1(v,e_D,d_Dvv,e_rmetdet,vort,rrearth)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  spherical vorticity of v
!
integer,parameter  :: np=4,real_kind=8
real(kind=real_kind), intent(in) :: v(np,np,2)
real(kind=8),dimension(np,np,2,2)	::   e_D
real(kind=8),dimension(np,np)	::   d_Dvv
real(kind=8),dimension(np,np)	::   e_rmetdet
real(kind=real_kind) :: vort(np,np)
real(kind=8)	::   rrearth
!local 
integer i
integer j
integer l

real(kind=real_kind) ::  dvdx00
real(kind=real_kind) ::  dudy00
real(kind=real_kind) ::  vco(np,np,2)
real(kind=real_kind) ::  vtemp(np,np)

! convert to covariant form
do j=1,np
	do i=1,np
		vco(i,j,1)=(e_D(i,j,1,1)*v(i,j,1) + e_D(i,j,2,1)*v(i,j,2))
		vco(i,j,2)=(e_D(i,j,1,2)*v(i,j,1) + e_D(i,j,2,2)*v(i,j,2))
	enddo
enddo

do j=1,np
	do l=1,np

		dudy00=0.0d0
		dvdx00=0.0d0
		!DIR$ UNROLL(NP)
		do i=1,np
			dvdx00 = dvdx00 + d_Dvv(i,l  )*vco(i,j  ,2)
			dudy00 = dudy00 + d_Dvv(i,l  )*vco(j  ,i,1)
		enddo

		vort(l  ,j  ) = dvdx00
		vtemp(j  ,l  ) = dudy00
	enddo
enddo

do j=1,np
	do i=1,np
		vort(i,j)=(vort(i,j)-vtemp(i,j))*(e_rmetdet(i,j)*rrearth)
	end do
end do

end subroutine 

subroutine gradient_sphere_wk_testcov1(s,ds,e_mp,e_metinv,d_Dvv,e_metdet,e_D,rrearth)
!
!   integrated-by-parts gradient, w.r.t. COVARIANT test functions
!   input s:  scalar
!   output  ds: weak gradient, lat/lon coordinates
!   ds = - integral[ div(PHIcov) s ]
!
!     PHIcov1 = (PHI^mn,0)   
!     PHIcov2 = (0,PHI^mn)
!   div() acts on contra components, so convert test function to contra: 
!     PHIcontra1 =  metinv PHIcov1  = (a^mn,b^mn)*PHI^mn   
!                                     a = metinv(1,1)  b=metinv(2,1)
!
!   ds1 = sum wij g  s_ij 1/g ( g a PHI^mn)_x  + ( g b PHI^mn)_y ) 
!       = sum  wij s_ij  ag(m,n)  d/dx( PHI^mn ) + bg(m,n) d/dy( PHI^mn)
!          i,j 
! for d/dx component, only sum over j=n
!       = sum  w_in s_in  ag(m,n)  d( PHI^m)(i)
!          i
! for d/dy component, only sum over i=m
!       = sum  w_mj s_mj  bg(m,n)  d( PHI^n)(j)
!          j
!  
!
! This formula is identical to gradient_sphere_wk_testcontra, except that
!    g(m,n) is replaced by a(m,n)*g(m,n)   
!  and we have two terms for each componet of ds 
!
!
 
integer,parameter  :: np=4,real_kind=8
real(kind=real_kind), intent(in) :: s(np,np)
real(kind=real_kind) :: ds(np,np,2)
real(kind=8),dimension(np,np)	::   e_mp
real(kind=8),dimension(np,np,2,2)	::  e_metinv
real(kind=8),dimension(np,np)	::   d_Dvv,e_metdet
real(kind=8),dimension(np,np,2,2)	::   e_D
real(kind=8)	::   rrearth
!local 
integer i,j,l,m,n
real(kind=real_kind) ::  dscontra(np,np,2)

dscontra=0
do n=1,np
	do m=1,np
		!DIR$ UNROLL(NP)
		do j=1,np
			dscontra(m,n,1)=dscontra(m,n,1)-(&
			(e_mp(j,n)*e_metinv(m,n,1,1)*e_metdet(m,n)*s(j,n)*d_Dvv(m,j) ) +&
			(e_mp(m,j)*e_metinv(m,n,2,1)*e_metdet(m,n)*s(m,j)*d_Dvv(n,j) ) &
			) *rrearth

			dscontra(m,n,2)=dscontra(m,n,2)-(&
			(e_mp(j,n)*e_metinv(m,n,1,2)*e_metdet(m,n)*s(j,n)*d_Dvv(m,j) ) +&
			(e_mp(m,j)*e_metinv(m,n,2,2)*e_metdet(m,n)*s(m,j)*d_Dvv(n,j) ) &
			) *rrearth
		enddo
	enddo
enddo
!--------------------commnets--------------------!
!#if 0
!    ! slow form, for debugging
!    do m=1,np
!       do n=1,np
!          vcontra=0
!          vcontra(m,n,1)=1
!
!          ! contra->latlon:
!          v(:,:,1)=(elem%D(1,1,:,:)*vcontra(:,:,1) + elem%D(1,2,:,:)*vcontra(:,:,2))
!          v(:,:,2)=(elem%D(2,1,:,:)*vcontra(:,:,1) + elem%D(2,2,:,:)*vcontra(:,:,2))
!
!
!          ! compute div(metdet phivec) * s
!          div = divergence_sphere(v,deriv,elem)
!          ! compute integral[ div(phi) * s ]
!          ds(m,n,1)=0
!          do i=1,np
!             do j=1,np
!                ds(m,n,1)=ds(m,n,1) + div(i,j)*s(i,j)*elem%spheremp(i,j)
!             enddo
!          enddo
!
!          vcontra=0
!          vcontra(m,n,2)=1
!
!          ! contra->latlon:
!          v(:,:,1)=(elem%D(1,1,:,:)*vcontra(:,:,1) + elem%D(1,2,:,:)*vcontra(:,:,2))
!          v(:,:,2)=(elem%D(2,1,:,:)*vcontra(:,:,1) + elem%D(2,2,:,:)*vcontra(:,:,2))
!
!          ! compute div(metdet phivec) * s
!          div = divergence_sphere(v,deriv,elem)
!          ! compute integral[ div(phi) * s ]
!          ds(m,n,2)=0
!          do i=1,np
!             do j=1,np
!                ds(m,n,2)=ds(m,n,2) + div(i,j)*s(i,j)*elem%spheremp(i,j)
!             enddo
!          enddo
!       enddo
!    enddo
!    ! change sign 
!    ds=-ds
!    print *,'ds,dscov:1 ',ds(1,1,1),dscov(1,1,1),ds(1,1,1)/dscov(1,1,1)
!    print *,'ds,dscov:2 ',ds(1,1,2),dscov(1,1,2),ds(1,1,2)/dscov(1,1,2)
!
!    dscov=ds
!#endif
!------------------------------------------------!

! convert contra -> latlon 
do j=1,np
	do i=1,np
		ds(i,j,1)=(e_D(i,j,1,1)*dscontra(i,j,1) + e_D(i,j,1,2)*dscontra(i,j,2))
		ds(i,j,2)=(e_D(i,j,2,1)*dscontra(i,j,1) + e_D(i,j,2,2)*dscontra(i,j,2))
	enddo
enddo

end subroutine 


subroutine curl_sphere_wk_testcov1(s,ds,d_Dvv,e_D,e_mp,rrearth)
!
!   integrated-by-parts gradient, w.r.t. COVARIANT test functions
!   input s:  scalar  (assumed to be s*khat)
!   output  ds: weak curl, lat/lon coordinates
!   
! starting with: 
!   PHIcov1 = (PHI,0)  covariant vector 
!   PHIcov2 = (0,PHI)  covariant vector 
!
!   ds1 = integral[ PHIcov1 dot curl(s*khat) ] 
!   ds2 = integral[ PHIcov2 dot curl(s*khat) ] 
! integrate by parts: 
!   ds1 = integral[ vor(PHIcov1) * s ]       
!   ds2 = integral[ vor(PHIcov1) * s ]
!
!     PHIcov1 = (PHI^mn,0)   
!     PHIcov2 = (0,PHI^mn)
!  vorticity() acts on covariant vectors:
!   ds1 = sum wij g  s_ij 1/g (  (PHIcov1_2)_x  - (PHIcov1_1)_y ) 
!       = -sum wij s_ij  d/dy (PHI^mn )
! for d/dy component, only sum over i=m
!       = -sum  w_mj s_mj   d( PHI^n)(j)
!           j
!
!   ds2 = sum wij g  s_ij 1/g (  (PHIcov2_2)_x  - (PHIcov2_1)_y ) 
!       = +sum wij s_ij  d/dx (PHI^mn )
! for d/dx component, only sum over j=n
!       = +sum  w_in s_in  d( PHI^m)(i)
!           i
!
integer,parameter  :: np=4,real_kind=8
real(kind=real_kind), intent(in) :: s(np,np)
real(kind=real_kind) :: ds(np,np,2)
real(kind=8),dimension(np,np)	::   d_Dvv
real(kind=8),dimension(np,np,2,2)	::   e_D
real(kind=8),dimension(np,np)	::   e_mp
real(kind=8)	::   rrearth

!local
integer i,j,l,m,n
real(kind=real_kind) ::  dscontra(np,np,2)

dscontra=0
do n=1,np
	do m=1,np
		!DIR$ UNROLL(NP)
		do j=1,np
			! phi(n)_y  sum over second index, 1st index fixed at m
			dscontra(m,n,1)=dscontra(m,n,1)-(e_mp(m,j)*s(m,j)*d_Dvv(n,j) )*rrearth
			! phi(m)_x  sum over first index, second index fixed at n
			dscontra(m,n,2)=dscontra(m,n,2)+(e_mp(j,n)*s(j,n)*d_Dvv(m,j) )*rrearth
		enddo
	enddo
enddo

! convert contra -> latlon 
do j=1,np
	do i=1,np
		ds(i,j,1)=(e_D(i,j,1,1)*dscontra(i,j,1) + e_D(i,j,1,2)*dscontra(i,j,2))
		ds(i,j,2)=(e_D(i,j,2,1)*dscontra(i,j,1) + e_D(i,j,2,2)*dscontra(i,j,2))
	enddo
enddo
end subroutine 

subroutine vlaplace_sphere_wk_gh_fortran(v,laplace,e_vec_sphere2cart,hypervis_power,&
hypervis_scaling,var_coef_int,nu_ratio,rrearth,e_variable_hyperviscosity,e_tensorVisc,&
d_Dvv,e_Dinv,e_spheremp,e_metdet,e_rmetdet,e_D,e_mp,e_metinv)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   logic:
!      tensorHV:     requires cartesian
!      nu_div/=nu:   requires contra formulatino
!
!   One combination NOT supported:  tensorHV and nu_div/=nu then abort
!
implicit none
integer,parameter :: np=4,real_kind=8
real(kind=real_kind),dimension(np,np,2)	:: v
real(kind=real_kind),dimension(np,np,2)	:: laplace
real(kind=real_kind),dimension(np,np,3,2)	:: e_vec_sphere2cart
real(kind=real_kind)		:: hypervis_power,hypervis_scaling
integer		:: var_coef_int
real(kind=real_kind)		:: nu_ratio,rrearth
real(kind=real_kind),dimension(np,np)	:: e_variable_hyperviscosity
real(kind=real_kind),dimension(np,np,2,2)	:: e_tensorVisc
real(kind=real_kind),dimension(np,np)	:: d_Dvv
real(kind=real_kind),dimension(np,np,2,2)	:: e_Dinv
real(kind=real_kind),dimension(np,np)	:: e_spheremp
real(kind=real_kind),dimension(np,np)	:: e_metdet
real(kind=real_kind),dimension(np,np)	:: e_rmetdet
real(kind=real_kind),dimension(np,np,2,2)	:: e_D
real(kind=real_kind),dimension(np,np)	:: e_mp
real(kind=real_kind),dimension(np,np,2,2)	:: e_metinv
!local 
logical :: var_coef
if(var_coef_int==0) then
	var_coef=.false.
else
	var_coef=.true.
endif

if (hypervis_scaling/=0 .and. var_coef) then
	! abortmp function is not used!
	! tensorHV is turned on - requires cartesian formulation
	!	if (present(nu_ratio)) then
	!		if (nu_ratio /= 1) then
	!			call abortmp('ERROR: tensorHV can not be used with nu_div/=nu')
	!		endif
	!	endif
	call vlaplace_sphere_wk_cartesian2(v,laplace,&
	e_vec_sphere2cart,rrearth,hypervis_power,&
	hypervis_scaling,var_coef,e_variable_hyperviscosity,&
	e_tensorVisc,d_Dvv,e_Dinv,e_spheremp)
else  
	! all other cases, use contra formulation:
	call  vlaplace_sphere_wk_contra2(v,laplace,d_Dvv,e_metdet,e_Dinv,&
	e_rmetdet, rrearth,e_D,e_mp,e_metinv,e_spheremp,&
	e_variable_hyperviscosity, hypervis_power, nu_ratio,var_coef)
endif
end subroutine 

subroutine  vlaplace_sphere_wk_cartesian2(v,laplace,e_vec_sphere2cart,rrearth,hypervis_power,&
hypervis_scaling,var_coef,e_variable_hyperviscosity,e_tensorVisc,d_Dvv,e_Dinv,e_spheremp)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
integer,parameter :: np=4,real_kind=8
real(kind=real_kind),dimension(np,np,2)	:: v
real(kind=real_kind),dimension(np,np,2)	:: laplace
real(kind=real_kind),dimension(np,np,3,2)	:: e_vec_sphere2cart
real(kind=real_kind)		:: rrearth
real(kind=real_kind)		:: hypervis_power,hypervis_scaling
logical		:: var_coef
real(kind=real_kind),dimension(np,np)	:: e_variable_hyperviscosity
real(kind=real_kind),dimension(np,np,2,2)	:: e_tensorVisc
real(kind=real_kind),dimension(np,np)	:: d_Dvv
real(kind=real_kind),dimension(np,np,2,2)	:: e_Dinv
real(kind=real_kind),dimension(np,np)	:: e_spheremp

! Local
integer component
real(kind=real_kind) :: dum_cart(np,np,3)
real(kind=real_kind) :: dum_tmp(np,np)

! latlon -> cartesian
do component=1,3
	dum_cart(:,:,component)=sum( e_vec_sphere2cart(:,:,component,:)*v(:,:,:) ,3)
end do

! Do laplace on cartesian comps
do component=1,3
	call laplace_sphere_wk2(dum_cart(:,:,component),d_Dvv,e_Dinv, e_spheremp,dum_tmp, &
	e_variable_hyperviscosity, e_tensorVisc, rrearth, hypervis_power, hypervis_scaling,var_coef)

	dum_cart(:,:,component) = dum_tmp
enddo

! cartesian -> latlon
do component=1,2
	! vec_sphere2cart is its own pseudoinverse.
	laplace(:,:,component)=sum( dum_cart(:,:,:)*e_vec_sphere2cart(:,:,:,component) ,3)
end do 

!	
end subroutine 

subroutine  laplace_sphere_wk2(s,d_Dvv,e_Dinv, e_spheremp, laplace, &
e_variable_hyperviscosity, e_tensorVisc, rrearth, hypervis_power, hypervis_scaling,var_coef)
!
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operator, grad(s) does not need to be made C0
!
implicit none
integer,parameter :: np=4,real_kind=8
real(kind=real_kind), intent(in) :: s(np,np) 
real(kind=real_kind),dimension(np,np)	:: d_Dvv
real(kind=real_kind),dimension(np,np,2,2)	:: e_Dinv
real(kind=real_kind),dimension(np,np)	:: e_spheremp
real(kind=real_kind)  :: laplace(np,np)
real(kind=real_kind),dimension(np,np)	:: e_variable_hyperviscosity
real(kind=real_kind),dimension(np,np,2,2)	:: e_tensorVisc
real(kind=real_kind)		:: rrearth,hypervis_power,hypervis_scaling
logical, intent(in) :: var_coef

! Local
real(kind=real_kind) :: grads(np,np,2), oldgrads(np,np,2)
integer i,j
call gradient_sphere2(e_Dinv,s,grads,d_Dvv,rrearth)
if (var_coef) then
	if (hypervis_power/=0 ) then
		! scalar viscosity with variable coefficient
		grads(:,:,1) = grads(:,:,1)*e_variable_hyperviscosity(:,:)
		grads(:,:,2) = grads(:,:,2)*e_variable_hyperviscosity(:,:)
	else if (hypervis_scaling /=0 ) then
		! tensor hv, (3)
		oldgrads=grads
		do j=1,np
			do i=1,np
				grads(i,j,1) = sum(oldgrads(i,j,:)*e_tensorVisc(i,j,1,:))
				grads(i,j,2) = sum(oldgrads(i,j,:)*e_tensorVisc(i,j,2,:))
			end do
		end do
	else
		! do nothing: constant coefficient viscsoity
	endif
endif

! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  
call  divergence_sphere_wk2(grads,laplace,d_Dvv,e_Dinv,e_spheremp,rrearth)
end subroutine 

subroutine gradient_sphere2(Dinv,s,ds,d_Dvv,rrearth)
!
!   input s:  scalar
!   output  ds: spherical gradient of s, lat-lon coordinates
!
implicit none
integer,parameter :: np=4,real_kind=8
real(kind=real_kind), intent(in), dimension(np,np,2,2) :: Dinv
real(kind=real_kind), intent(in) :: s(np,np)
real(kind=real_kind) :: ds(np,np,2)
real(kind=real_kind), intent(in), dimension(np,np) :: d_DVV
real(kind=real_kind) :: rrearth
!local 
integer i
integer j
integer l
real(kind=real_kind) ::  dsdx00
real(kind=real_kind) ::  dsdy00
real(kind=real_kind) ::  v1(np,np),v2(np,np)

do j=1,np
	do l=1,np
		dsdx00=0.0d0
		dsdy00=0.0d0
		!DIR$ UNROLL(NP)
		do i=1,np
			dsdx00 = dsdx00 + d_Dvv(i,l  )*s(i,j  )
			dsdy00 = dsdy00 + d_Dvv(i,l  )*s(j  ,i)
		end do
		v1(l  ,j  ) = dsdx00*rrearth
		v2(j  ,l  ) = dsdy00*rrearth
	end do
end do
! convert covarient to latlon
!OMP_COLLAPSE_SIMD
!DIR_VECTOR_ALIGNED
do j=1,np
	do i=1,np
		ds(i,j,1)=Dinv(i,j,1,1)*v1(i,j) + Dinv(i,j,2,1)*v2(i,j)
		ds(i,j,2)=Dinv(i,j,1,2)*v1(i,j) + Dinv(i,j,2,2)*v2(i,j)
	enddo
enddo

end subroutine 

subroutine  divergence_sphere_wk2(v,div,d_Dvv,e_Dinv,e_spheremp,rrearth)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!
!   Computes  -< grad(psi) dot v > 
!   (the integrated by parts version of < psi div(v) > )
!
!   note: after DSS, divergence_sphere() and divergence_sphere_wk() 
!   are identical to roundoff, as theory predicts.
!
implicit none
integer,parameter :: np=4,real_kind=8
real(kind=real_kind), intent(in),dimension(np,np,2) :: v
real(kind=real_kind),dimension(np,np) :: div 
real(kind=real_kind),dimension(np,np)	:: d_Dvv
real(kind=real_kind), intent(in) ,dimension(np,np,2,2)	:: e_Dinv
real(kind=real_kind), intent(in) ,dimension(np,np)	:: e_spheremp
real(kind=real_kind)	:: rrearth

! Local
integer i,j,m,n
real(kind=real_kind) ::  vtemp(np,np,2)
real(kind=real_kind) ::  ggtemp(np,np,2)
real(kind=real_kind) ::  gtemp(np,np,2)
real(kind=real_kind) ::  psi(np,np)
real(kind=real_kind) :: xtmp

! latlon- > contra
!OMP_COLLAPSE_SIMD
!DIR_VECTOR_ALIGNED
do j=1,np
	do i=1,np
		vtemp(i,j,1)=(e_Dinv(i,j,1,1)*v(i,j,1) + e_Dinv(i,j,1,2)*v(i,j,2))
		vtemp(i,j,2)=(e_Dinv(i,j,2,1)*v(i,j,1) + e_Dinv(i,j,2,2)*v(i,j,2))
	enddo
enddo

do n=1,np
	do m=1,np

		div(m,n)=0
		!DIR$ UNROLL(NP)
		do j=1,np
			div(m,n)=div(m,n)-(e_spheremp(j,n)*vtemp(j,n,1)*d_Dvv(m,j) &
			+e_spheremp(m,j)*vtemp(m,j,2)*d_Dvv(n,j)) &
			* rrearth
		enddo
	end do
end do
end subroutine  

subroutine vlaplace_sphere_wk_contra2(v,laplace,d_Dvv,e_metdet,e_Dinv,&
e_rmetdet, rrearth,e_D,e_mp,e_metinv,e_spheremp,&
e_variable_hyperviscosity, hypervis_power, nu_ratio,var_coef)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   du/dt = laplace(u) = grad(div) - curl(vor)
!   < PHI du/dt > = < PHI laplace(u) >        PHI = covariant, u = contravariant
!                 = < PHI grad(div) >  - < PHI curl(vor) >
!                 = grad_wk(div) - curl_wk(vor)               
!
implicit none
integer,parameter :: np=4,real_kind=8
real(kind=real_kind), intent(in) :: v(np,np,2) 
real(kind=real_kind) :: laplace(np,np,2)
real(kind=real_kind),dimension(np,np)	:: d_Dvv
real(kind=real_kind),dimension(np,np)	:: e_metdet
real(kind=real_kind),dimension(np,np,2,2)	:: e_Dinv
real(kind=real_kind),dimension(np,np)	:: e_rmetdet
real(kind=real_kind) :: rrearth 
real(kind=real_kind),dimension(np,np,2,2)	:: e_D
real(kind=real_kind),dimension(np,np)	:: e_mp
real(kind=real_kind),dimension(np,np,2,2)	:: e_metinv 
real(kind=real_kind),dimension(np,np)	:: e_spheremp
real(kind=real_kind),dimension(np,np)	:: e_variable_hyperviscosity
real(kind=real_kind) :: hypervis_power,nu_ratio
logical, intent(in) :: var_coef
! Local
integer i,j,l,m,n
real(kind=real_kind) :: lap_tmp(np,np,2)
real(kind=real_kind) :: lap_tmp2(np,np,2)
real(kind=real_kind) :: vor(np,np),div(np,np)
real(kind=real_kind) :: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y

call  divergence_sphere2(v,div,d_Dvv,e_metdet,e_Dinv,e_rmetdet,rrearth)
call vorticity_sphere2(v,e_D,d_Dvv,e_rmetdet,vor,rrearth)
if (var_coef .and. hypervis_power/=0 ) then
	! scalar viscosity with variable coefficient
	div = div*e_variable_hyperviscosity(:,:)
	vor = vor*e_variable_hyperviscosity(:,:)
endif

!if (present(nu_ratio))
div = nu_ratio*div

call  gradient_sphere_wk_testcov2(div,lap_tmp,e_mp,e_metdet,e_metinv,d_Dvv,e_D,rrearth)
call  curl_sphere_wk_testcov2(vor,lap_tmp2,d_Dvv,e_mp,e_D,rrearth)
laplace = lap_tmp - lap_tmp2
do n=1,np
	do m=1,np
		! add in correction so we dont damp rigid rotation
#define UNDAMPRR
#ifdef UNDAMPRR
		laplace(m,n,1)=laplace(m,n,1) + 2*e_spheremp(m,n)*v(m,n,1)*(rrearth**2)
		laplace(m,n,2)=laplace(m,n,2) + 2*e_spheremp(m,n)*v(m,n,2)*(rrearth**2)
#endif
	enddo
enddo

end subroutine

subroutine divergence_sphere2(v,div,d_Dvv,e_metdet,e_Dinv,e_rmetdet,rrearth)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
!
implicit none
integer,parameter		::np=4,real_kind=8
real(kind=real_kind), intent(in),dimension(np,np,2) :: v  ! in lat-lon coordinates
real(kind=real_kind) :: div(np,np)
real(kind=real_kind), intent(in),dimension(np,np) :: d_Dvv  ! in lat-lon coordinates
real(kind=real_kind), intent(in),dimension(np,np) :: e_metdet  ! in lat-lon coordinates
real(kind=real_kind), intent(in),dimension(np,np,2,2) :: e_Dinv  ! in lat-lon coordinates
real(kind=real_kind), intent(in),dimension(np,np) :: e_rmetdet  ! in lat-lon coordinates
real(kind=real_kind) :: rrearth
! Local
integer i
integer j
integer l
real(kind=real_kind) ::  dudx00
real(kind=real_kind) ::  dvdy00
real(kind=real_kind) ::  gv(np,np,2),vvtemp(np,np)
! convert to contra variant form and multiply by g
do j=1,np
	do i=1,np
		gv(i,j,1)=e_metdet(i,j)*(e_Dinv(i,j,1,1)*v(i,j,1) + e_Dinv(i,j,1,2)*v(i,j,2))
		gv(i,j,2)=e_metdet(i,j)*(e_Dinv(i,j,2,1)*v(i,j,1) + e_Dinv(i,j,2,2)*v(i,j,2))
	enddo
enddo

! compute d/dx and d/dy         
do j=1,np
	do l=1,np
		dudx00=0.0d0
		dvdy00=0.0d0
		do i=1,np
			dudx00 = dudx00 + d_Dvv(i,l  )*gv(i,j  ,1)
			dvdy00 = dvdy00 + d_Dvv(i,l  )*gv(j  ,i,2)
		end do
		div(l  ,j  ) = dudx00
		vvtemp(j  ,l  ) = dvdy00
	end do
end do

do j=1,np
	do i=1,np
		div(i,j)=(div(i,j)+vvtemp(i,j))*(e_rmetdet(i,j)*rrearth)
	end do
end do
end subroutine 

subroutine vorticity_sphere2(v,e_D,d_Dvv,e_rmetdet,vort,rrearth)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  spherical vorticity of v
!
implicit none
integer,parameter		::np=4,real_kind=8
real(kind=real_kind), intent(in) :: v(np,np,2)
real(kind=real_kind) :: vort(np,np)
real(kind=real_kind), intent(in),dimension(np,np,2,2) :: e_D
real(kind=real_kind), intent(in),dimension(np,np) :: d_Dvv
real(kind=real_kind), intent(in),dimension(np,np) :: e_rmetdet
real(kind=real_kind), intent(in) :: rrearth
!local	
integer i
integer j
integer l

real(kind=real_kind) ::  dvdx00
real(kind=real_kind) ::  dudy00
real(kind=real_kind) ::  vco(np,np,2)
real(kind=real_kind) ::  vtemp(np,np)

! convert to covariant form
do j=1,np
	do i=1,np
		vco(i,j,1)=(e_D(i,j,1,1)*v(i,j,1) + e_D(i,j,2,1)*v(i,j,2))
		vco(i,j,2)=(e_D(i,j,1,2)*v(i,j,1) + e_D(i,j,2,2)*v(i,j,2))
	enddo
enddo

do j=1,np
	do l=1,np

		dudy00=0.0d0
		dvdx00=0.0d0
		!DIR$ UNROLL(NP)
		do i=1,np
			dvdx00 = dvdx00 + d_Dvv(i,l  )*vco(i,j  ,2)
			dudy00 = dudy00 + d_Dvv(i,l  )*vco(j  ,i,1)
		enddo

		vort(l  ,j  ) = dvdx00
		vtemp(j  ,l  ) = dudy00
	enddo
enddo

do j=1,np
	do i=1,np
		vort(i,j)=(vort(i,j)-vtemp(i,j))*(e_rmetdet(i,j)*rrearth)
	end do
end do
end subroutine 

subroutine gradient_sphere_wk_testcov2(s,ds,e_mp,e_metdet,e_metinv,d_Dvv,e_D,rrearth)
!
!   integrated-by-parts gradient, w.r.t. COVARIANT test functions
!   input s:  scalar
!   output  ds: weak gradient, lat/lon coordinates
!   ds = - integral[ div(PHIcov) s ]
!
!     PHIcov1 = (PHI^mn,0)   
!     PHIcov2 = (0,PHI^mn)
!   div() acts on contra components, so convert test function to contra: 
!     PHIcontra1 =  metinv PHIcov1  = (a^mn,b^mn)*PHI^mn   
!                                     a = metinv(1,1)  b=metinv(2,1)
!
!   ds1 = sum wij g  s_ij 1/g ( g a PHI^mn)_x  + ( g b PHI^mn)_y ) 
!       = sum  wij s_ij  ag(m,n)  d/dx( PHI^mn ) + bg(m,n) d/dy( PHI^mn)
!          i,j 
! for d/dx component, only sum over j=n
!       = sum  w_in s_in  ag(m,n)  d( PHI^m)(i)
!          i
! for d/dy component, only sum over i=m
!       = sum  w_mj s_mj  bg(m,n)  d( PHI^n)(j)
!          j
!  
!
! This formula is identical to gradient_sphere_wk_testcontra, except that
!    g(m,n) is replaced by a(m,n)*g(m,n)   
!  and we have two terms for each componet of ds 
!
!
implicit none
integer,parameter :: np=4,real_kind=8 
real(kind=real_kind), intent(in) :: s(np,np)
real(kind=real_kind) :: ds(np,np,2)
real(kind=real_kind),dimension(np,np) :: e_mp
real(kind=real_kind),dimension(np,np) :: e_metdet 
real(kind=real_kind),dimension(np,np,2,2) :: e_metinv
real(kind=real_kind),dimension(np,np) :: d_Dvv
real(kind=real_kind),dimension(np,np,2,2) :: e_D
real(kind=real_kind)	:: rrearth
!local 	
integer i,j,l,m,n
real(kind=real_kind) ::  dscontra(np,np,2)
dscontra=0
do n=1,np
	do m=1,np
		!DIR$ UNROLL(NP)
		do j=1,np
			dscontra(m,n,1)=dscontra(m,n,1)-(&
			(e_mp(j,n)*e_metinv(m,n,1,1)*e_metdet(m,n)*s(j,n)*d_Dvv(m,j) ) +&
			(e_mp(m,j)*e_metinv(m,n,2,1)*e_metdet(m,n)*s(m,j)*d_Dvv(n,j) ) &
			) *rrearth

			dscontra(m,n,2)=dscontra(m,n,2)-(&
			(e_mp(j,n)*e_metinv(m,n,1,2)*e_metdet(m,n)*s(j,n)*d_Dvv(m,j) ) +&
			(e_mp(m,j)*e_metinv(m,n,2,2)*e_metdet(m,n)*s(m,j)*d_Dvv(n,j) ) &
			) *rrearth
		enddo
	enddo
enddo
! convert contra -> latlon 
do j=1,np
	do i=1,np
		ds(i,j,1)=(e_D(i,j,1,1)*dscontra(i,j,1) + e_D(i,j,1,2)*dscontra(i,j,2))
		ds(i,j,2)=(e_D(i,j,2,1)*dscontra(i,j,1) + e_D(i,j,2,2)*dscontra(i,j,2))
	enddo
enddo

end subroutine 

subroutine curl_sphere_wk_testcov2(s,ds,d_Dvv,e_mp,e_D,rrearth)
!
!   integrated-by-parts gradient, w.r.t. COVARIANT test functions
!   input s:  scalar  (assumed to be s*khat)
!   output  ds: weak curl, lat/lon coordinates
!   
! starting with: 
!   PHIcov1 = (PHI,0)  covariant vector 
!   PHIcov2 = (0,PHI)  covariant vector 
!
!   ds1 = integral[ PHIcov1 dot curl(s*khat) ] 
!   ds2 = integral[ PHIcov2 dot curl(s*khat) ] 
! integrate by parts: 
!   ds1 = integral[ vor(PHIcov1) * s ]       
!   ds2 = integral[ vor(PHIcov1) * s ]
!
!     PHIcov1 = (PHI^mn,0)   
!     PHIcov2 = (0,PHI^mn)
!  vorticity() acts on covariant vectors:
!   ds1 = sum wij g  s_ij 1/g (  (PHIcov1_2)_x  - (PHIcov1_1)_y ) 
!       = -sum wij s_ij  d/dy (PHI^mn )
! for d/dy component, only sum over i=m
!       = -sum  w_mj s_mj   d( PHI^n)(j)
!           j
!
!   ds2 = sum wij g  s_ij 1/g (  (PHIcov2_2)_x  - (PHIcov2_1)_y ) 
!       = +sum wij s_ij  d/dx (PHI^mn )
! for d/dx component, only sum over j=n
!       = +sum  w_in s_in  d( PHI^m)(i)
!           i
!
implicit none
integer,parameter :: np=4,real_kind=8 
real(kind=real_kind), intent(in) :: s(np,np)
real(kind=real_kind) :: ds(np,np,2)
real(kind=real_kind),dimension(np,np) :: d_Dvv
real(kind=real_kind),dimension(np,np) :: e_mp
real(kind=real_kind),dimension(np,np,2,2) :: e_D
real(kind=real_kind) :: rrearth
!local 
integer i,j,l,m,n
real(kind=real_kind) ::  dscontra(np,np,2)

dscontra=0
do n=1,np
	do m=1,np
		!DIR$ UNROLL(NP)
		do j=1,np
			! phi(n)_y  sum over second index, 1st index fixed at m
			dscontra(m,n,1)=dscontra(m,n,1)-(e_mp(m,j)*s(m,j)*d_Dvv(n,j) )*rrearth
			! phi(m)_x  sum over first index, second index fixed at n
			dscontra(m,n,2)=dscontra(m,n,2)+(e_mp(j,n)*s(j,n)*d_Dvv(m,j) )*rrearth
		enddo
	enddo
enddo

! convert contra -> latlon 
do j=1,np
	do i=1,np
		ds(i,j,1)=(e_D(i,j,1,1)*dscontra(i,j,1) + e_D(i,j,1,2)*dscontra(i,j,2))
		ds(i,j,2)=(e_D(i,j,2,1)*dscontra(i,j,1) + e_D(i,j,2,2)*dscontra(i,j,2))
	enddo
enddo
end subroutine 


