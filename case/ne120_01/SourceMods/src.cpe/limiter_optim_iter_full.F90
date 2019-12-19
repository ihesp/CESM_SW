subroutine limiter_optim_iter_full_f(ptens,sphweights,minp,maxp,dpmass,k_n)
implicit none
  integer, parameter :: np=4, nlev=30, real_kind = 8
  integer, intent(in) :: k_n
  real (kind=real_kind), dimension(k_n), intent(inout)   :: minp, maxp
  real (kind=real_kind), dimension(np*np,k_n), intent(inout)   :: ptens
  real (kind=real_kind), dimension(np*np,k_n), intent(in), optional  :: dpmass
  real (kind=real_kind), dimension(np*np), intent(in)   :: sphweights

  real (kind=real_kind), dimension(np,np) :: ptens_mass
  integer  k1, k, i, j, iter, weightsnum
  real (kind=real_kind) :: addmass, weightssum, mass, sumc
  real (kind=real_kind) :: x(np*np),c(np*np)
  integer :: maxiter = np*np-1
  real (kind=real_kind) :: tol_limiter = 5D-14

  do k = 1, k_n
    do k1=1,np*np
      c(k1)=sphweights(k1)*dpmass(k1,k)
      x(k1)=ptens(k1,k)/dpmass(k1,k)
    enddo

    sumc=sum(c)
    if (sumc <= 0 ) CYCLE   ! this should never happen, but if it does, dont limit
    mass=sum(c*x)

    ! relax constraints to ensure limiter has a solution:
    ! This is only needed if runnign with the SSP CFL>1 or
    ! due to roundoff errors
    if( mass < minp(k)*sumc ) then
      minp(k) = mass / sumc
    endif
    if( mass > maxp(k)*sumc ) then
      maxp(k) = mass / sumc
    endif

    do iter=1,maxiter

      addmass=0.0d0

      do k1=1,np*np
        if((x(k1)>maxp(k))) then
          addmass=addmass+(x(k1)-maxp(k))*c(k1)
          x(k1)=maxp(k)
        endif
        if((x(k1)<minp(k))) then
          addmass=addmass-(minp(k)-x(k1))*c(k1)
          x(k1)=minp(k)
        endif
      enddo !k1

      if(abs(addmass)<=tol_limiter*abs(mass)) exit

      weightssum=0.0d0
!        weightsnum=0
      if(addmass>0)then
        do k1=1,np*np
          if(x(k1)<maxp(k))then
            weightssum=weightssum+c(k1)
!              weightsnum=weightsnum+1
          endif
        enddo !k1
        do k1=1,np*np
          if(x(k1)<maxp(k))then
              x(k1)=x(k1)+addmass/weightssum
!                x(k1)=x(k1)+addmass/(c(k1)*weightsnum)
          endif
        enddo
      else
        do k1=1,np*np
          if(x(k1)>minp(k))then
            weightssum=weightssum+c(k1)
!              weightsnum=weightsnum+1
          endif
        enddo
        do k1=1,np*np
          if(x(k1)>minp(k))then
            x(k1)=x(k1)+addmass/weightssum
!             x(k1)=x(k1)+addmass/(c(k1)*weightsnum)
          endif
        enddo
      endif
    enddo!end of iteration

    do k1=1,np*np
       ptens(k1,k)=x(k1)
    enddo
  enddo

  do k = 1, k_n
    do k1=1,np*np
      ptens(k1,k)=ptens(k1,k)*dpmass(k1,k)
    enddo
  enddo
#if 0
  print *, tol_limiter*5D-3*5.82989072416270755D-293
#endif
end subroutine limiter_optim_iter_full_f
