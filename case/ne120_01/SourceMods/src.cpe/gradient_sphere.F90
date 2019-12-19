subroutine gradient_sphere_fortran(s,d_Dvv,Dinv,ds,rrearth)
		implicit none
		integer,parameter :: np=4,real_kind=8
		real(kind=real_kind) :: ds(np,np,2)
		real(kind=real_kind), intent(in), dimension(np,np,2,2) :: Dinv
		real(kind=real_kind), intent(in) :: s(np,np)
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

end subroutine gradient_sphere_fortran

subroutine divergence_sphere_fortran(v,div,d_Dvv,e_metdet,e_Dinv,e_rmetdet,rrearth)
				
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
end subroutine divergence_sphere_fortran

subroutine vorticity_sphere_fortran(v,e_D,d_Dvv,e_rmetdet,vort,rrearth)

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
end subroutine vorticity_sphere_fortran


