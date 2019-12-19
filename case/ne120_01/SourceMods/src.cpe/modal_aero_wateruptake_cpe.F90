!/*************************************************************************
!	> File Name: modal_aero_wateruptake_cpe.F90
!	> Author: Xu Kai
!	> Created Time: 2019年01月08日 星期二 21时26分14秒
! ************************************************************************/

!      subroutine pow_zzf(n1)
!
!      complex*16, intent(inout)  :: n1(2)
!
!      n1(2) = n1(1) ** (1.d0/3.d0)
!      n1(1) = n1(2)
!      end subroutine
!-----------------------------------------------------------------------
      subroutine makoh_cubic( cx, p2, p1, p0, im )
!
!     solves  x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2 are real
!

      integer, parameter :: imx=5
      integer :: im
      real(8) :: p0(imx), p1(imx), p2(imx)
      complex*16 :: cx(3,imx)

      integer :: i
      real(8) :: eps, q(imx), r(imx), sqrt3, third
      complex*16 :: ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

      !save eps
      !data eps/1.e-20_r8/

      eps = 1.d-20
      third=1.d0/3.d0
      ci=dcmplx(0.d0,1.d0)
      !ci=cmplx(0.d0,1.d0,kind=16)
      sqrt3=sqrt(3.d0)
      cw=0.5d0*(-1+ci*sqrt3)
      cwsq=0.5d0*(-1-ci*sqrt3)

      do i=1,im
      if(p1(i).eq.0.d0)then
!        completely insoluble particle
         cx(1,i)=(-p0(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
      else
         q(i)=p1(i)/3.d0
         r(i)=p0(i)/2.d0
         crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
         crad(i)=sqrt(crad(i))

         cy(i)=r(i)-crad(i)
         if (abs(cy(i)).gt.eps) cy(i)=cy(i)**third
         cq=q(i)
         cz(i)=-cq/cy(i)

         cx(1,i)=-cy(i)-cz(i)
         cx(2,i)=-cw*cy(i)-cwsq*cz(i)
         cx(3,i)=-cwsq*cy(i)-cw*cz(i)
      endif
      enddo
      return
      end subroutine makoh_cubic


!-----------------------------------------------------------------------
      subroutine makoh_quartic( cx, p3, p2, p1, p0, im )

!     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2, p3 are real
!
      integer, parameter :: imx=5
      integer :: im
      real(8) :: p0(imx), p1(imx), p2(imx), p3(imx)
      complex*16 :: cx(4,imx)

      integer :: i
      real(8) :: third, q(imx), r(imx)
      complex*16 :: cb(imx), cb0(imx), cb1(imx),   &
                     crad(imx), cy(imx), czero
	  
	  real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
	  real(8) :: pi

	  !complex*16 :: haha 

	  pi = 3.14159265358979323846264338327950288419716939937510d0


      !czero=cmplx(0.0d0,0.0d0,kind=16)
      czero=dcmplx(0.0d0,0.0d0)
      third=1.d0/3.d0

	  !haha = dcmplx(third, 0.d0)

      do 10 i=1,im

      q(i)=-p2(i)*p2(i)/36.d0+(p3(i)*p1(i)-4*p0(i))/12.d0
      !r(i)=-(p2(i)/6)**3+p2(i)*(p3(i)*p1(i)-4*p0(i))/48._r8   &
      ! +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16
      r(i)=-(p2(i)/6)*(p2(i)/6)*(p2(i)/6)+p2(i)*(p3(i)*p1(i)-4*p0(i))/48.d0   &
       +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16

      crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
      crad(i)=sqrt(crad(i))

      cb(i)=r(i)-crad(i)
      if(cb(i).eq.czero)then
!        insoluble particle
         cx(1,i)=(-p1(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
         cx(4,i)=cx(1,i)
      else

         !cx(1,i)=cb(i) 
         cb(i)=cb(i)**third
		 !call get_double(third, tmp1)
		 !cx(2,i) = tmp1
		 
         !cb(i)=cb(i)**haha

		 !tmp1 = dreal(cb(i))
		 !tmp2 = dimag(cb(i))
		 !tmp3 = sqrt(tmp1*tmp1 + tmp2*tmp2)
		 !tmp4 = asin(tmp2/tmp3)
		 !tmp3 = tmp3 ** third


		 !tmp5 = cos(third*tmp4)
		 !tmp6 = sin(third*tmp4)
		 !!tmp6 = -1.d0 * sqrt(1 - tmp5*tmp5)

		 !tmp7 = tmp3*tmp5
		 !tmp8 = tmp3*tmp6

         !cx(1,i)=cb(i) 
		 !!cb(i) = dcmplx(tmp7, tmp8) 

		 !tmp1 = dreal(cb(i))
		 !tmp2 = dimag(cb(i))
		 !tmp3 = sqrt(tmp1*tmp1 + tmp2*tmp2)
		 !tmp4 = atan2(tmp2, tmp1) !logi
		 !!tmp4 = asin(tmp2/tmp3)

		 !tmp5 = log(tmp3) ! logr
		 !tmp6 = exp(tmp5*third - tmp4*tmp2)
		 !tmp7 = tmp5*tmp2 + tmp4*third
		 !tmp8 = tmp6 * cos(tmp7)
		 !tmp9 = tmp6 * sin(tmp7)

		 !cb(i) = dcmplx(tmp8, tmp9) 

         !cx(2,i)=cb(i) 
         !cx(3,i)=cb(i) 
         !cx(4,i)=cb(i) 


         cy(i)=-cb(i)+q(i)/cb(i)+p2(i)/6.d0

         cb0(i)=sqrt(cy(i)*cy(i)-p0(i))
         cb1(i)=(p3(i)*cy(i)-p1(i))/(2*cb0(i))

         cb(i)=p3(i)/2.d0+cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)+cb0(i))
         crad(i)=sqrt(crad(i))
         cx(1,i)=(-cb(i)+crad(i))/2.d0
         cx(2,i)=(-cb(i)-crad(i))/2.d0

         cb(i)=p3(i)/2-cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)-cb0(i))
         crad(i)=sqrt(crad(i))
         cx(3,i)=(-cb(i)+crad(i))/2.d0
         cx(4,i)=(-cb(i)-crad(i))/2.d0
      endif
   10 continue

         !cx(1,1)=p3(1) 
         !cx(2,1)=p2(1) 
         !cx(3,1)=p1(1) 
         !cx(4,1)=p0(1) 
      return
      end subroutine makoh_quartic
!-----------------------------------------------------------------------
      subroutine modal_aero_kohler(   &
          rdry_in, hygro, s, rwet_out, im)

! calculates equlibrium radius r of haze droplets as function of
! dry particle mass and relative humidity s using kohler solution
! given in pruppacher and klett (eqn 6-35)

! for multiple aerosol types, assumes an internal mixture of aerosols

      implicit none
      interface 
          subroutine makoh_cubic( cx, p2, p1, p0, im )
          complex*16 :: cx(3,5)
          real(8) :: p2(5), p1(5), p0(5)
          integer :: im
          end subroutine 
          subroutine makoh_quartic( cx, p3, p2, p1, p0, im )
          complex*16 :: cx(4,5)
          real(8) :: p3(5), p2(5), p1(5), p0(5)
          integer :: im
          end subroutine
      end interface
! arguments
      integer :: im         ! number of grid points to be processed
      real(8) :: rdry_in(:)    ! aerosol dry radius (m)
      real(8) :: hygro(:)      ! aerosol volume-mean hygroscopicity (--)
      real(8) :: s(:)          ! relative humidity (1 = saturated)
      real(8) :: rwet_out(:)   ! aerosol wet radius (m)

! local variables
      integer, parameter :: imax=5
      integer :: i, n, nsol

      real(8) :: a, b
      real(8) :: p40(imax),p41(imax),p42(imax),p43(imax) ! coefficients of polynomial
      real(8) :: p30(imax),p31(imax),p32(imax) ! coefficients of polynomial
      real(8) :: p
      real(8) :: r3, r4
      real(8) :: r(im)         ! wet radius (microns)
      real(8) :: rdry(imax)    ! radius of dry particle (microns)
      real(8) :: ss            ! relative humidity (1 = saturated)
      real(8) :: slog(imax)    ! log relative humidity
      real(8) :: vol(imax)     ! total volume of particle (microns**3)
      real(8) :: xi, xr

      complex*16 :: cx4(4,imax),cx3(3,imax)

      real(8) :: eps
      real(8) :: mw
      real(8) :: pi
      real(8) :: rhow
      real(8) :: surften
      real(8) :: tair
      real(8) :: third
      real(8) :: ugascon

      !integer(8) :: log_agent
      eps = 1.d-4
      mw = 18.d0
      pi = 3.14159d0
      rhow = 1.d0
      surften = 76.d0
      tair = 273.d0
      third = 1.d0/3.d0
      ugascon = 8.3d7   
      !call math_agent_log_c(log_agent)

!     effect of organics on surface tension is neglected
      a=2.d4*mw*surften/(ugascon*tair*rhow)

      do i=1,im
           rdry(i) = rdry_in(i)*1.0d6   ! convert (m) to (microns)
           vol(i) = rdry(i)*rdry(i)*rdry(i)          ! vol is r**3, not volume
           b = vol(i)*hygro(i)

!          quartic
           ss=min(s(i),1.d0-eps)
           ss=max(ss,1.d-10)
           slog(i)=log(ss)
           !call math_agent_1i1o(log_agent, ss, slog(i))
           p43(i)=-a/slog(i)
           p42(i)=0.d0
           p41(i)=b/slog(i)-vol(i)
           p40(i)=a*vol(i)/slog(i)
!          cubic for rh=1
           p32(i)=0.d0
           p31(i)=-b/a
           p30(i)=-vol(i)

      end do

	  r(1) = 0.d0
       do 100 i=1,im

!       if(vol(i).le.1.e-20)then
        if(vol(i).le.1.d-12)then
           r(i)=rdry(i)
           go to 100
        endif

        p=abs(p31(i))/(rdry(i)*rdry(i))
        if(p.lt.eps)then
!          approximate solution for small particles
           r(i)=rdry(i)*(1.d0+p*third/(1.d0-slog(i)*rdry(i)/a))
        else
           call makoh_quartic(cx4,p43,p42,p41,p40,1)
           !call makoh_quartic(cx4(1,i),p43(i),p42(i),p41(i),p40(i),1)
           !cx4(1,i) = p43(i)
           !cx4(2,i) = p42(i)
           !cx4(3,i) = p41(i)
           !cx4(4,i) = p40(i)
!          find smallest real(r8) solution
           r(i)=1000.d0*rdry(i)
           nsol=0
           do n=1,4
              xr=real(cx4(n,i))
              xi=aimag(cx4(n,i))
              if(abs(xi).gt.abs(xr)*eps) cycle  
              if(xr.gt.r(i)) cycle  
              if(xr.lt.rdry(i)*(1.d0-eps)) cycle  
              if(xr.ne.xr) cycle  
              r(i)=xr
              nsol=n
           end do  
           if(nsol.eq.0)then
              !write(iulog,*)   &
              ! 'ccm kohlerc - no real(r8) solution found (quartic)'
              !write(iulog,*)'roots =', (cx4(n,i),n=1,4)
              !write(iulog,*)'p0-p3 =', p40(i), p41(i), p42(i), p43(i)
              !write(iulog,*)'rh=',s(i)
              !write(iulog,*)'setting radius to dry radius=',rdry(i)
              r(i)=rdry(i)
!             stop
           endif
		   !r(1) = dreal(cx4(1,1))
		   !r(2) = dimag(cx4(1,1))
		   !r(3) = dreal(cx4(2,1))
		   !r(4) = dimag(cx4(2,1))
        endif

        if(s(i).gt.1.d0-eps)then
!          save quartic solution at s=1-eps
           r4=r(i)
!          cubic for rh=1
           p=abs(p31(i))/(rdry(i)*rdry(i))
           if(p.lt.eps)then
              r(i)=rdry(i)*(1.d0+p*third)
           else
              call makoh_cubic(cx3,p32,p31,p30,im)
              !call makoh_cubic(cx3(1,i),p32(i),p31(i),p30(i),im)
!             find smallest real(r8) solution
              r(i)=1000.d0*rdry(i)
              nsol=0
              do n=1,3
                 xr=real(cx3(n,i))
                 xi=aimag(cx3(n,i))
                 if(abs(xi).gt.abs(xr)*eps) cycle  
                 if(xr.gt.r(i)) cycle  
                 if(xr.lt.rdry(i)*(1.d0-eps)) cycle  
                 if(xr.ne.xr) cycle  
                 r(i)=xr
                 nsol=n
              end do  
              if(nsol.eq.0)then
                 !write(iulog,*)   &
                 ! 'ccm kohlerc - no real(r8) solution found (cubic)'
                 !write(iulog,*)'roots =', (cx3(n,i),n=1,3)
                 !write(iulog,*)'p0-p2 =', p30(i), p31(i), p32(i)
                 !write(iulog,*)'rh=',s(i)
                 !write(iulog,*)'setting radius to dry radius=',rdry(i)
                 r(i)=rdry(i)
!                stop
              endif
           endif
           r3=r(i)
!          now interpolate between quartic, cubic solutions
           r(i)=(r4*(1.d0-s(i))+r3*(s(i)-1.d0+eps))/eps
        endif

  100 continue

! bound and convert from microns to m
      do i=1,im
         r(i) = min(r(i),30.d0) ! upper bound based on 1 day lifetime
         rwet_out(i) = r(i)*1.d-6
      end do

!rwet_out(1) = r(1)
!rwet_out2(1) = r(3)
!rwet_out3(1) = r(4)


      return
      end subroutine modal_aero_kohler


subroutine modal_aero_wateruptake_sub_compute( &
    myid, pcols, pver, top_lev, modal_strat_sulfate, &
    ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
    hygro, rh, dryvol, so4dryvol, so4specdens, troplev, &
    wetrad, wetvol, wtrvol, sulden, wtpct, pi, third, pi43)

    implicit none

    interface 
        subroutine modal_aero_kohler(rdry_in, hygro, s, rwet_out, im)
        integer :: im
        real(8) :: rdry_in(:), hygro(:), s(:), rwet_out(:)
        end subroutine
    end interface

    ! Arguments
    integer, intent(in)  :: myid
    integer, intent(in)  :: ncol                    ! number of columns
    integer, intent(in)  :: nmodes
    integer, intent(in)  :: pcols
    integer, intent(in)  :: pver
    integer, intent(in)  :: top_lev
    logical, intent(in)  :: modal_strat_sulfate
    integer, intent(in)  :: troplev

    real(8), intent(in) :: pi
    real(8), intent(in) :: third
    real(8), intent(in) :: pi43 
    real(8), intent(in) :: rhcrystal(16)
    real(8), intent(in) :: rhdeliques(16)
    real(8), intent(in) :: dryrad(30,16)         ! dry volume mean radius of aerosol (m)
    real(8), intent(in) :: hygro(30,16)          ! volume-weighted mean hygroscopicity (--)
    real(8), intent(in) :: rh(30)               ! relative humidity (0-1)
    real(8), intent(in) :: dryvol(30,16)         ! dry volume of single aerosol (m3)
    real(8), intent(in) :: so4dryvol(30,16)      ! dry volume of sulfate in single aerosol (m3)
    real(8), intent(in) :: so4specdens           ! mass density sulfate in single aerosol (kg/m3)
    real(8), intent(in) :: wtpct(30,16)          ! sulfate aerosol composition, weight % H2SO4
    real(8), intent(in) :: sulden(30,16)         ! sulfate aerosol mass density (g/cm3)

    real(8), intent(inout) :: wetrad(30,16)        ! wet radius of aerosol (m)
    real(8), intent(inout) :: wetvol(30,16)        ! single-particle-mean wet volume (m3)
    real(8), intent(inout) :: wtrvol(30,16)        ! single-particle-mean water volume in wet aerosol (m3)

    ! local variables

    integer :: i, k, m

    real(8) :: hystfac                ! working variable for hysteresis

    integer(8) :: pow_agent
    real(8) :: tmp1, tmp2
    !-----------------------------------------------------------------------

    !call math_agent_pow_c(pow_agent)

    ! loop over all aerosol modes
    do m = 1, nmodes

       hystfac = 1.0d0 / max(1.0d-5, (rhdeliques(m) - rhcrystal(m)))

       do k = top_lev, pver

             if ( modal_strat_sulfate .and. (k<troplev)) then
                wetvol(k,m) = dryvol(k,m)-so4dryvol(k,m)
                wetvol(k,m) = wetvol(k,m)+so4dryvol(k,m)*so4specdens/sulden(k,m)/wtpct(k,m)/10.d0
                wetvol(k,m) = max(wetvol(k,m), dryvol(k,m))
                !tmp1 = (wetvol(k,m)/pi43)
                !call math_agent_2i1o(pow_agent, tmp1, third, wetrad(k,m))
                wetrad(k,m) = (wetvol(k,m)/pi43)**third
                wetrad(k,m) = max(wetrad(k,m), dryrad(k,m))
                wtrvol(k,m) = wetvol(k,m) - dryvol(k,m)
                wtrvol(k,m) = max(wtrvol(k,m), 0.0d0)
             else
               ! compute wet radius for each mode
               !call modal_aero_kohler(dryrad(k:k,m), hygro(k:k,m), rh(k:k), wetrad(k:k,m), 1)
               call modal_aero_kohler(dryrad(k:k,m), hygro(k:k,m), rh(k:k), wetrad(k:k,m), 1)

               !wetrad(k,m) = dryrad(k,m) + hygro(k,m) + rh(k)

               wetrad(k,m) = max(wetrad(k,m), dryrad(k,m))
               wetvol(k,m) = pi43*wetrad(k,m)*wetrad(k,m)*wetrad(k,m)
               wetvol(k,m) = max(wetvol(k,m), dryvol(k,m))
               wtrvol(k,m) = wetvol(k,m) - dryvol(k,m)
               wtrvol(k,m) = max(wtrvol(k,m), 0.0d0)

               ! apply simple treatment of deliquesence/crystallization hysteresis
               ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
               ! the "upper curve" value, and the fraction is a linear function of rh
               if (rh(k) < rhcrystal(m)) then
                  wetrad(k,m) = dryrad(k,m)
                  wetvol(k,m) = dryvol(k,m)
                  wtrvol(k,m) = 0.0d0
               else if (rh(k) < rhdeliques(m)) then
                  wtrvol(k,m) = wtrvol(k,m)*hystfac*(rh(k) - rhcrystal(m))
                  wtrvol(k,m) = max(wtrvol(k,m), 0.0d0)
                  wetvol(k,m) = dryvol(k,m) + wtrvol(k,m)
                  !tmp1 = (wetvol(k,m)/pi43)
                  !call math_agent_2i1o(pow_agent, tmp1, third, wetrad(k,m))
                  wetrad(k,m) = (wetvol(k,m)/pi43)**third
               end if
             end if

       end do     ! levels

    end do ! modes

    end subroutine
	
