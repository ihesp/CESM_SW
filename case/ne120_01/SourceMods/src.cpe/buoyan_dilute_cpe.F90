!*************************************************************************
! File Name: src.cpe/buoyan_dilute_cpe.f90
! Author: Xu Kai
! Created Time: 2018年11月12日 星期一 16时33分22秒
!************************************************************************/
        module buoyan_dilute_cpe
        public buoyan_dilute_compute
        public cldprp 
        public closure 
        public q1q2_pjr
        contains

        subroutine cldprp(lchnk   , &
                  q       ,t       ,u       ,v       ,p       , &
                  z       ,s       ,mu      ,eu      ,du      , &
                  md      ,ed      ,sd      ,qd      ,mc      , &
                  qu      ,su      ,zf      ,qst     ,hmn     , &
                  hsat    ,shat    ,ql      , &
                  cmeg    ,jb      ,lel     ,jt      ,jlcl    , &
                  mx      ,j0      ,jd      ,rl      ,il2g    , &
                  rd      ,grav    ,cp      ,msg     , &
                  pflx    ,evp     ,cu      ,rprd    ,limcnv  ,landfrac,&
                  pver, pverp, tfreez, eps1, rh2o, c0_ocn, c0_lnd, &
                  log_addr,log10_addr, pow_addr&
                  )
!               qtnd, dlf, heat, mcon, pflxo, cme, &
!               zdu,  rprdo)
!
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.
! 
! Author: See above, modified by P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
!-----------------------------------------------------------------------

        implicit none

!------------------------------------------------------------------------------
!
! Input arguments
!
       integer, intent(in) :: pver 
       integer, intent(in) :: pverp
       integer, intent(in) :: lchnk                  ! chunk identifier
    
       real(8), intent(in) :: q(pver)         ! spec. humidity of env
       real(8), intent(in) :: t(pver)         ! temp of env
       real(8), intent(in) :: p(pver)         ! pressure of env
       real(8), intent(in) :: z(pver)         ! height of env
       real(8), intent(in) :: s(pver)         ! normalized dry static energy of env
       real(8), intent(in) :: zf(pverp)       ! height of interfaces
       real(8), intent(in) :: u(pver)         ! zonal velocity of env
       real(8), intent(in) :: v(pver)         ! merid. velocity of env
    
       real(8), intent(in) :: landfrac ! RBN Landfrac
    
       integer, intent(in) :: jb              ! updraft base level
       integer, intent(in) :: lel             ! updraft launch level
       integer, intent(out) :: jt              ! updraft plume top
       integer, intent(out) :: jlcl            ! updraft lifting cond level
       integer, intent(in) :: mx              ! updraft base level (same is jb)
       integer, intent(out) :: j0              ! level where updraft begins detraining
       integer, intent(out) :: jd              ! level of downdraft
       integer, intent(in) :: limcnv                 ! convection limiting level
       integer, intent(in) :: il2g                   !CORE GROUP REMOVE
       integer, intent(in) :: msg                    ! missing moisture vals (always 0)
       real(8), intent(in) :: rl                    ! latent heat of vap
       real(8), intent(in) :: shat(pver)      ! interface values of dry stat energy
!
! output
!
       real(8), intent(out) :: rprd(pver)     ! rate of production of precip at that layer
       real(8), intent(out) :: du(pver)       ! detrainement rate of updraft
       real(8), intent(out) :: ed(pver)       ! entrainment rate of downdraft
       real(8), intent(out) :: eu(pver)       ! entrainment rate of updraft
       real(8), intent(out) :: hmn(pver)      ! moist stat energy of env
       real(8), intent(out) :: hsat(pver)     ! sat moist stat energy of env
       real(8), intent(out) :: mc(pver)       ! net mass flux
       real(8), intent(out) :: md(pver)       ! downdraft mass flux
       real(8), intent(out) :: mu(pver)       ! updraft mass flux
       real(8), intent(out) :: pflx(pverp)    ! precipitation flux thru layer
       real(8), intent(out) :: qd(pver)       ! spec humidity of downdraft
       real(8), intent(out) :: ql(pver)       ! liq water of updraft
       real(8), intent(out) :: qst(pver)      ! saturation mixing ratio of env.
       real(8), intent(out) :: qu(pver)       ! spec hum of updraft
       real(8), intent(out) :: sd(pver)       ! normalized dry stat energy of downdraft
       real(8), intent(out) :: su(pver)       ! normalized dry stat energy of updraft
    
    
       real(8) rd                   ! gas constant for dry air
       real(8) grav                 ! gravity
       real(8) cp                   ! heat capacity of dry air

!
! Local workspace
!
       real(8) gamma(pver)
       real(8) dz(pver)
       real(8) iprm(pver)
       real(8) hu(pver)
       real(8) hd(pver)
       real(8) eps(pver)
       real(8) f(pver)
       real(8) k1(pver)
       real(8) i2(pver)
       real(8) ihat(pver)
       real(8) i3(pver)
       real(8) idag(pver)
       real(8) i4(pver)
       real(8) qsthat(pver)
       real(8) hsthat(pver)
       real(8) gamhat(pver)
       real(8) cu(pver)
       real(8) evp(pver)
       real(8) cmeg(pver)
       real(8) qds(pver)
! RBN For c0mask
       real(8) c0mask
    
       real(8) hmin
       real(8) expdif
       real(8) expnum
       real(8) ftemp
       real(8) eps0
       real(8) rmue
       real(8) zuef
       real(8) zdef
       real(8) epsm
       real(8) ratmjb
       real(8) est
       real(8) totpcp
       real(8) totevp
       real(8) alfa
       real(8) ql1
       real(8) tu
       real(8) estu
       real(8) qstu
    
       real(8) small
       real(8) mdt
    
       integer khighest
       integer klowest
       integer kount
       integer i,k
    
       logical doit
       logical done
       real(8), intent(in) :: tfreez, eps1, rh2o

       real(8), intent(in) :: c0_ocn
       real(8), intent(in) :: c0_lnd
       real(8),parameter ::  tiedke_add = 5.0d-1 
       integer(8), intent(in) :: log_addr,log10_addr, pow_addr
       real(8) :: login, logout
       real(8) :: powbase
       real(8) :: powexp
       real(8) :: powout1, powout2, powout3, powout4, powout5,&
                    powout6, powout7

!      real(kind=8), intent(out) :: qtnd(pver) 
!      real(kind=8), intent(out) :: heat(pver)           
!      real(kind=8), intent(out) :: mcon(pverp)
!      real(kind=8), intent(out) :: dlf(pver)   
!      real(kind=8), intent(out) :: pflxo(pverp)
!      real(kind=8), intent(out) :: cme(pver)
!      real(kind=8), intent(out) :: zdu(pver)
!      real(kind=8), intent(out) :: rprdo(pver) 
      integer :: kcount_sig
      real(8) :: tmp1, tmp2

!

!------------------------------------------------------------------------------
!
!   do i = 1,il2g
      ftemp = 0.0d0
      expnum = 0.0d0
      expdif = 0.0d0
      c0mask  = c0_ocn * (1.0d0-landfrac) +   c0_lnd * landfrac 
!   end do
!
!jr Change from msg+1 to 1 to prevent blowup
!
        do k = 1,pver
      !do i = 1,il2g
         dz(k) = zf(k) - zf(k+1)
      !end do
        end do

!
! initialize many output and work variables to zero
!
       pflx(1) = 0

        do k = 1,pver
!      do i = 1,il2g
         k1(k) = 0.0d0
         i2(k) = 0.0d0
         i3(k) = 0.0d0
         i4(k) = 0.0d0
         mu(k) = 0.0d0
         f(k) = 0.0d0
         eps(k) = 0.0d0
         eu(k) = 0.0d0
         du(k) = 0.0d0
         ql(k) = 0.0d0
         cu(k) = 0.0d0
         evp(k) = 0.0d0
         cmeg(k) = 0.0d0
         qds(k) = q(k)
         md(k) = 0.0d0
         ed(k) = 0.0d0
         sd(k) = s(k)
         qd(k) = q(k)
         mc(k) = 0.0d0
         qu(k) = q(k)
         su(k) = s(k)
         call qsat_hPa(t(k), p(k), est, qst(k), tfreez, eps1, rl, rh2o, cp, log_addr,log10_addr, pow_addr)
!++bee
         if ( p(k)-est <= 0.0d0 ) then
            qst(k) = 1.0d0
         end if
!--bee

         gamma(k) = qst(k)*(1.0d0 + qst(k)/eps1)*eps1*rl/(rd*(t(k)*t(k)))*rl/cp
         hmn(k) = cp*t(k) + grav*z(k) + rl*q(k)
         hsat(k) = cp*t(k) + grav*z(k) + rl*qst(k)
         hu(k) = hmn(k)
         hd(k) = hmn(k)
         rprd(k) = 0.0d0
!      end do
       end do

!
!jr Set to zero things which make this routine blow up

       do k=1,msg
!      do i=1,il2g
         rprd(k) = 0.0d0
!      end do
       end do
!
! interpolate the layer values of qst, hsat and gamma to
! layer interfaces
!
   do k = 1, msg+1
!      do i = 1,il2g
         hsthat(k) = hsat(k)
         qsthat(k) = qst(k)
         gamhat(k) = gamma(k)
!      end do
   end do
!   do i = 1,il2g
      totpcp = 0.0d0
      totevp = 0.0d0
!   end do
   do k = msg + 2,pver
!      do i = 1,il2g
         if (abs(qst(k-1)-qst(k)) > 1.d-6) then
             login = qst(k-1)/qst(k)
             call math_agent_1i1o(log_addr, login, logout)
            qsthat(k) = logout*qst(k-1)*qst(k)/ (qst(k-1)-qst(k))
         else
            qsthat(k) = qst(k)
         end if
         hsthat(k) = cp*shat(k) + rl*qsthat(k)
         if (abs(gamma(k-1)-gamma(k)) > 1.d-6) then
             login = gamma(k-1)/gamma(k)
             call math_agent_1i1o(log_addr, login, logout)
            gamhat(k) = logout*gamma(k-1)*gamma(k)/ &
                                (gamma(k-1)-gamma(k))
         else
            gamhat(k) = gamma(k)
         end if


!      end do
   end do


!
! initialize cloud top to highest plume top.
!jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
!
   jt = pver
!   do i = 1,il2g
      jt = max(lel,limcnv+1)
      jt = min(jt,pver)
      jd = pver
      jlcl = lel
      hmin = 1.0d6
!   end do
!
! find the level of minimum hsat, where detrainment starts
!

   do k = msg + 1,pver
!      do i = 1,il2g
         if (hsat(k) <= hmin .and. k >= jt .and. k <= jb) then
            hmin = hsat(k)
            j0 = k
         end if
!      end do
   end do
!   do i = 1,il2g
      j0 = min(j0,jb-2)
      j0 = max(j0,jt+2)
!
! Fix from Guang Zhang to address out of bounds array reference
!
      j0 = min(j0,pver)
!   end do
!
! Initialize certain arrays inside cloud
!
   do k = msg + 1,pver
!      do i = 1,il2g
         if (k >= jt .and. k <= jb) then
            hu(k) = hmn(mx) + cp*tiedke_add
            su(k) = s(mx) + tiedke_add
         end if
!      end do
   end do

!
! *********************************************************
! compute taylor series for approximate eps(z) below
! *********************************************************
!
   do k = pver - 1,msg + 1,-1
!      do i = 1,il2g
         if (k < jb .and. k >= jt) then
            k1(k) = k1(k+1) + (hmn(mx)-hmn(k))*dz(k)
            ihat(k) = 0.5d0* (k1(k+1)+k1(k))
            i2(k) = i2(k+1) + ihat(k)*dz(k)
            idag(k) = 0.5d0* (i2(k+1)+i2(k))
            i3(k) = i3(k+1) + idag(k)*dz(k)
            iprm(k) = 0.5d0* (i3(k+1)+i3(k))
            i4(k) = i4(k+1) + iprm(k)*dz(k)
         end if
!      end do
   end do

!
! re-initialize hmin array for ensuing calculation.
!
!   do i = 1,il2g
      hmin = 1.d6
!   end do
   do k = msg + 1,pver
!      do i = 1,il2g
         if (k >= j0 .and. k <= jb .and. hmn(k) <= hmin) then
            hmin = hmn(k)
            expdif = hmn(mx) - hmin
         end if
!      end do
   end do
!
! *********************************************************
! compute approximate eps(z) using above taylor series
! *********************************************************
!
   do k = msg + 2,pver
!      do i = 1,il2g
         expnum = 0.d0
         ftemp = 0.d0
         if (k < jt .or. k >= jb) then
            k1(k) = 0.d0
            expnum = 0.d0
         else
            expnum = hmn(mx) - (hsat(k-1)*(zf(k)-z(k)) + &
                        hsat(k)* (z(k-1)-zf(k)))/(z(k-1)-z(k))
         end if
         if ((expdif > 100.d0 .and. expnum > 0.d0) .and. &
            k1(k) > expnum*dz(k)) then
            ftemp = expnum/k1(k)
            powbase = ftemp
            powexp = 2
            call math_agent_2i1o(pow_addr, powbase, powexp,  powout1) 
            powexp = 3
            call math_agent_2i1o(pow_addr, powbase, powexp,  powout2) 
            powexp = 4
            call math_agent_2i1o(pow_addr, powbase, powexp,  powout3) 
            powbase = k1(k)
            powexp = 2
            call math_agent_2i1o(pow_addr, powbase, powexp,  powout4) 
            powexp = 3
            call math_agent_2i1o(pow_addr, powbase, powexp,  powout5) 
 
            powbase =i2(k)
            powexp = 2
            call math_agent_2i1o(pow_addr, powbase, powexp,  powout6) 
            powexp = 3
            call math_agent_2i1o(pow_addr, powbase, powexp,  powout7) 

         f(k) = ftemp + i2(k)/k1(k)*powout1 + &
                     (2.d0*powout6-k1(k)*i3(k))/powout4* &
                     powout2 + (-5.d0*k1(k)*i2(k)*i3(k)+ &
                     5.d0*powout7+ powout4*i4(k))/ &
                     powout5*powout3
!   
            !f(k) = ftemp + i2(k)/k1(k)*(ftemp*ftemp) + &
            !         (2.d0*(i2(k)*i2(k))-k1(k)*i3(k))/(k1(k)*k1(k))* &
            !         (ftemp*ftemp*ftemp) + (-5.d0*k1(k)*i2(k)*i3(k)+ &
            !         5.d0*(i2(k)*i2(k)*i2(k))+(k1(k)*k1(k))*i4(k))/ &
            !         (k1(k)*k1(k)*k1(k))*(ftemp*ftemp*ftemp*ftemp)

!            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 + &
!                     (2._r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2* &
!                     ftemp(i)**3 + (-5._r8*k1(i,k)*i2(i,k)*i3(i,k)+ &
!                     5._r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/ &
!                     k1(i,k)**3*ftemp(i)**4
 
            f(k) = max(f(k),0.d0)
            f(k) = min(f(k),2.0d-4)
         end if


!      end do
   end do

!   do i = 1,il2g
      if (j0 < jb) then
         if (f(j0) < 1.d-6 .and. f(j0+1) > f(j0)) j0 = j0 + 1
      end if
!   end do
   do k = msg + 2,pver
!      do i = 1,il2g
         if (k >= jt .and. k <= j0) then
            f(k) = max(f(k),f(k-1))
         end if
!      end do
   end do
!   do i = 1,il2g
      eps0 = f(j0)
      eps(jb) = eps0
!   end do
!
! This is set to match the Rasch and Kristjansson paper
!
   do k = pver,msg + 1,-1
!      do i = 1,il2g
         if (k >= j0 .and. k <= jb) then
            eps(k) = f(j0)
         end if
!      end do
   end do
   do k = pver,msg + 1,-1
!      do i = 1,il2g
         if (k < j0 .and. k >= jt) eps(k) = f(k)
!      end do
   end do
!
! specify the updraft mass flux mu, entrainment eu, detrainment du
! and moist static energy hu.
! here and below mu, eu,du, md and ed are all normalized by mb
!
!   do i = 1,il2g
      if (eps0 > 0.d0) then
         mu(jb) = 1.d0
         eu(jb) = mu(jb)/dz(jb)
      end if
!   end do
   do k = pver,msg + 1,-1
!      do i = 1,il2g
         if (eps0 > 0.d0 .and. (k >= jt .and. k < jb)) then
            zuef = zf(k) - zf(jb)
            rmue = (1.d0/eps0)* (exp(eps(k+1)*zuef)-1.d0)/zuef
            mu(k) = (1.d0/eps0)* (exp(eps(k  )*zuef)-1.d0)/zuef
            eu(k) = (rmue-mu(k+1))/dz(k)
            du(k) = (rmue-mu(k))/dz(k)
         end if
!      end do
   end do

!
   khighest = pverp
   klowest = 1

      khighest = min(khighest,lel)
      call dma_min(khighest) 

      klowest = max(klowest,jb)
      call dma_max(klowest) 

       !qtnd(1)  = khighest 
        !qtnd(2)  = klowest 
        !return 
!   end do
   do k = klowest-1,khighest,-1
!      do i = 1,il2g
         if (k <= jb-1 .and. k >= lel .and. eps0 > 0.d0) then
            if (mu(k) < 2.0d-2) then
               hu(k) = hmn(k)
               mu(k) = 0.d0
               eu(k) = 0.d0
               du(k) = mu(k+1)/dz(k)
            else
               hu(k) = mu(k+1)/mu(k)*hu(k+1) + &
                         dz(k)/mu(k)* (eu(k)*hmn(k)- du(k)*hsat(k))
            end if
         end if
!      end do
   end do


!
! reset cloud top index beginning from two layers above the
! cloud base (i.e. if cloud is only one layer thick, top is not reset
!
!   do i=1,il2g
      doit = .true.
!   end do
   do k=klowest-2,khighest-1,-1
!      do i=1,il2g
         if (doit .and. k <= jb-2 .and. k >= lel-1) then
            if (hu(k) <= hsthat(k) .and. hu(k+1) > hsthat(k+1) &
               .and. mu(k) >= 2.d-2) then
               if (hu(k)-hsthat(k) < -2000.d0) then
                  jt = k + 1
                  doit = .false.
               else
                  jt = k
                  doit = .false.
               end if
            else if (hu(k) > hu(jb) .or. mu(k) < 2.d-2) then
               jt = k + 1
               doit = .false.
            end if
         end if
!      end do
   end do


   do k = pver,msg + 1,-1
!      do i = 1,il2g
         if (k >= lel .and. k <= jt .and. eps0 > 0.d0) then
            mu(k) = 0.d0
            eu(k) = 0.d0
            du(k) = 0.d0
            hu(k) = hmn(k)
         end if
         if (k == jt .and. eps0 > 0.d0) then
            du(k) = mu(k+1)/dz(k)
            eu(k) = 0.d0
            mu(k) = 0.d0
         end if
!      end do
   end do

  
!
! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.
!
!   do i = 1,il2g
!
! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
!
      alfa = 0.1d0
      jt = min(jt,jb-1)
      jd = max(j0,jt+1)
      jd = min(jd,jb)
      hd(jd) = hmn(jd-1)
      if (jd < jb .and. eps0 > 0.d0) then
         epsm = eps0
         md(jd) = -alfa*epsm/eps0
      end if
!   end do
   do k = msg + 1,pver
!      do i = 1,il2g
         if ((k > jd .and. k <= jb) .and. eps0 > 0.d0) then
            zdef = zf(jd) - zf(k)
            md(k) = -alfa/ (2.d0*eps0)*(exp(2.d0*epsm*zdef)-1.d0)/zdef
         end if
!      end do
   end do
   do k = msg + 1,pver
!      do i = 1,il2g
         if ((k >= jt .and. k <= jb) .and. eps0 > 0.d0 .and. jd < jb) then
            ratmjb = min(abs(mu(jb)/md(jb)),1.d0)
            md(k) = md(k)*ratmjb
         end if
!      end do
   end do

   small = 1.d-20
   do k = msg + 1,pver
!      do i = 1,il2g
         if ((k >= jt .and. k <= pver) .and. eps0 > 0.d0) then
            ed(k-1) = (md(k-1)-md(k))/dz(k-1)
            mdt = min(md(k),-small)
            hd(k) = (md(k-1)*hd(k-1) - dz(k-1)*ed(k-1)*hmn(k-1))/mdt
         end if
!      end do
   end do
!
! calculate updraft and downdraft properties.
!
   do k = msg + 2,pver
!      do i = 1,il2g
         if ((k >= jd .and. k <= jb) .and. eps0 > 0.d0 .and. jd < jb) then
            qds(k) = qsthat(k) + gamhat(k)*(hd(k)-hsthat(k))/ &
               (rl*(1.d0 + gamhat(k)))
         end if
!      end do
   end do

!
!   do i = 1,il2g
      done = .false.
!   end do
   kount = 0

   do k = pver,msg + 2,-1
!      do i = 1,il2g
         if (k == jb .and. eps0 > 0.d0) then
            qu(k) = q(mx)
            su(k) = (hu(k)-rl*qu(k))/cp
         end if
         if (( .not. done .and. k > jt .and. k < jb) .and. eps0 > 0.d0) then
            su(k) = mu(k+1)/mu(k)*su(k+1) + &
                      dz(k)/mu(k)* (eu(k)-du(k))*s(k)
            qu(k) = mu(k+1)/mu(k)*qu(k+1) + dz(k)/mu(k)* (eu(k)*q(k)- &
                            du(k)*qst(k))
            tu = su(k) - grav/cp*zf(k)
            call qsat_hPa(tu, (p(k)+p(k-1))/2.d0, estu, qstu, tfreez, eps1, rl, rh2o, cp, log_addr,log10_addr, pow_addr)
            if (qu(k) >= qstu) then
               jlcl = k
               !kount = kount + 1
               done = .true.
                kcount_sig = 1
                call dma_sum(kcount_sig, kount)
             else
                kcount_sig = 0
                call dma_sum(kcount_sig, kount)
            end if
         else
             kcount_sig = 0
             call dma_sum(kcount_sig, kount)
         end if
!      end do
      if (kount >= il2g) goto 690
   end do
690 continue
         
   do k = msg + 2,pver
!      do i = 1,il2g
         if ((k > jt .and. k <= jlcl) .and. eps0 > 0.d0) then
            su(k) = shat(k) + (hu(k)-hsthat(k))/(cp* (1.d0+gamhat(k)))
            qu(k) = qsthat(k) + gamhat(k)*(hu(k)-hsthat(k))/ &
                     (rl* (1.d0+gamhat(k)))
         end if
!      end do
   end do

 do k = pver,msg + 2,-1
!      do i = 1,il2g
         if (k >= jt .and. k < jb .and. eps0 > 0.d0) then
            cu(k) = ((mu(k)*su(k)-mu(k+1)*su(k+1))/ &
                      dz(k)- (eu(k)-du(k))*s(k))/(rl/cp)
            if (k == jt) cu(k) = 0.d0
            cu(k) = max(0.d0,cu(k))
        end if
!      end do
   end do
   
                 return 


! compute condensed liquid, rain production rate
! accumulate total precipitation (condensation - detrainment of liquid)
! Note ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
! The differencing is somewhat strange (e.g. du(i,k)*ql(i,k+1)) but is
! consistently applied.
!    mu, ql are interface quantities
!    cu, du, eu, rprd are midpoint quantites

   do k = pver,msg + 2,-1
!      do i = 1,il2g
         rprd(k) = 0.d0
         if (k >= jt .and. k < jb .and. eps0 > 0.d0 .and. mu(k) >= 0.0d0) then
            if (mu(k) > 0.d0) then
               ql1 = 1.d0/mu(k)* (mu(k+1)*ql(k+1)- &
                     dz(k)*du(k)*ql(k+1)+dz(k)*cu(k))
               ql(k) = ql1/ (1.d0+dz(k)*c0mask)
            else
               ql(k) = 0.d0
            end if
            totpcp = totpcp + dz(k)*(cu(k)-du(k)*ql(k+1))
            rprd(k) = c0mask*mu(k)*ql(k)
         end if
!      end do
   end do


!
!   do i = 1,il2g
      qd(jd) = qds(jd)
      sd(jd) = (hd(jd) - rl*qd(jd))/cp
!   end do
!
   do k = msg + 2,pver
!      do i = 1,il2g
         if (k >= jd .and. k < jb .and. eps0 > 0.d0) then
            qd(k+1) = qds(k+1)
            evp(k) = -ed(k)*q(k) + (md(k)*qd(k)-md(k+1)*qd(k+1))/dz(k)
            evp(k) = max(evp(k),0.d0)
            mdt = min(md(k+1),-small)
            sd(k+1) = ((rl/cp*evp(k)-ed(k)*s(k))*dz(k) + md(k)*sd(k))/mdt
            totevp = totevp - dz(k)*ed(k)*q(k)
         end if
!      end do
   end do
!   do i = 1,il2g
!*guang         totevp = totevp + md(i,jd)*q(i,jd-1) -
      totevp = totevp + md(jd)*qd(jd) - md(jb)*qd(jb)
!   end do
!!$   if (.true.) then
   if (.false.) then
!      do i = 1,il2g
         k = jb
         if (eps0 > 0.d0) then
            evp(k) = -ed(k)*q(k) + (md(k)*qd(k))/dz(k)
            evp(k) = max(evp(k),0.d0)
            totevp = totevp - dz(k)*ed(k)*q(k)
         end if
!      end do
   endif

!   do i = 1,il2g
      totpcp = max(totpcp,0.d0)
      totevp = max(totevp,0.d0)
!   end do
!
   do k = msg + 2,pver
!      do i = 1,il2g
         if (totevp > 0.d0 .and. totpcp > 0.d0) then
            md(k)  = md (k)*min(1.d0, totpcp/(totevp+totpcp))
            ed(k)  = ed (k)*min(1.d0, totpcp/(totevp+totpcp))
            evp(k) = evp(k)*min(1.d0, totpcp/(totevp+totpcp))
         else
            md(k) = 0.d0
            ed(k) = 0.d0
            evp(k) = 0.d0
         end if
! cmeg is the cloud water condensed - rain water evaporated
! rprd is the cloud water converted to rain - (rain evaporated)
         cmeg(k) = cu(k) - evp(k)
         rprd(k) = rprd(k)-evp(k)
!      end do
   end do

! compute the net precipitation flux across interfaces
   pflx(1) = 0.d0
   do k = 2,pverp
!      do i = 1,il2g
         pflx(k) = pflx(k-1) + rprd(k-1)*dz(k-1)
!      end do
   end do
!
   do k = msg + 1,pver
!      do i = 1,il2g
         mc(k) = mu(k) + md(k)
!      end do
   end do
!
   return
end subroutine cldprp
         subroutine closure(lchnk   , &
                   q       ,t       ,p       ,z       ,s       , &
                   tp      ,qs      ,qu      ,su      ,mc      , &
                   du      ,mu      ,md      ,qd      ,sd      , &
                   qhat    ,shat    ,dp      ,qstp    ,zf      , &
                   ql      ,dsubcld ,mb      ,cape    ,tl      , &
                   lcl     ,lel     ,jt      ,mx      ,il1g    , &
                   il2g    ,rd      ,grav    ,cp      ,rl      , &
                   msg     ,pver    ,eps1, &
                   log_addr, log10_addr, pow_addr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: G. Zhang and collaborators. CCM contact:P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
! 
!-----------------------------------------------------------------------
!   use dycore,    only: dycore_is, get_resolution

   implicit none

!
!-----------------------------Arguments---------------------------------
!
       integer, intent(in) :: pver
       integer, intent(in) :: lchnk                 ! chunk identifier
    
       real(8), intent(inout) :: q(pver)        ! spec humidity
       real(8), intent(inout) :: t(pver)        ! temperature
       real(8), intent(inout) :: p(pver)        ! pressure (mb)
       real(8), intent(inout) :: mb            ! cloud base mass flux
       real(8), intent(in) :: z(pver)        ! height (m)
       real(8), intent(in) :: s(pver)        ! normalized dry static energy
       real(8), intent(in) :: tp(pver)       ! parcel temp
       real(8), intent(in) :: qs(pver)       ! sat spec humidity
       real(8), intent(in) :: qu(pver)       ! updraft spec. humidity
       real(8), intent(in) :: su(pver)       ! normalized dry stat energy of updraft
       real(8), intent(in) :: mc(pver)       ! net convective mass flux
       real(8), intent(in) :: du(pver)       ! detrainment from updraft
       real(8), intent(in) :: mu(pver)       ! mass flux of updraft
       real(8), intent(in) :: md(pver)       ! mass flux of downdraft
       real(8), intent(in) :: qd(pver)       ! spec. humidity of downdraft
       real(8), intent(in) :: sd(pver)       ! dry static energy of downdraft
       real(8), intent(in) :: qhat(pver)     ! environment spec humidity at interfaces
       real(8), intent(in) :: shat(pver)     ! env. normalized dry static energy at intrfcs
       real(8), intent(in) :: dp(pver)       ! pressure thickness of layers
       real(8), intent(in) :: qstp(pver)     ! spec humidity of parcel
       real(8), intent(in) :: zf(pver+1)     ! height of interface levels
       real(8), intent(in) :: ql(pver)       ! liquid water mixing ratio
    
       real(8), intent(in) :: cape          ! available pot. energy of column
       real(8), intent(in) :: tl
       real(8), intent(in) :: dsubcld       ! thickness of subcloud layer
    
       integer, intent(in) :: lcl        ! index of lcl
       integer, intent(in) :: lel        ! index of launch leve
       integer, intent(in) :: jt         ! top of updraft
       integer, intent(in) :: mx         ! base of updraft
!
!--------------------------Local variables------------------------------
!
       real(8) dtpdt(pver)
       real(8) dqsdtp(pver)
       real(8) dtmdt(pver)
       real(8) dqmdt(pver)
       real(8) dboydt(pver)
       real(8) thetavp(pver)
       real(8) thetavm(pver)
    
       real(8) dtbdt,dqbdt,dtldt
       real(8) beta
       !real(8) capelmt
       real(8) cp
       real(8) dadt
       real(8) debdt
       real(8) dltaa
       real(8) eb
       real(8) grav
    
       integer i
       integer il1g
       integer il2g
       integer k, kmin, kmax
       integer msg
    
       real(8) rd
       real(8) rl
       real(8), intent(in) :: eps1 
       real(8), parameter :: tau = 3600.d0
       real(8), parameter :: capelmt = 70.d0
       real(8) :: tmpres
       integer(8), intent(in) :: log_addr,log10_addr, pow_addr
       real(8) :: login, login1, login2, logout, logout1, logout2
       real(8) :: powbase1, powexp1, powout1, powbase2, powexp2, powout2
   ! threshold value for cape for deep convection.
! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.
!
!   do i = il1g,il2g
      mb = 0.d0
      eb = p(mx)*q(mx)/ (eps1+q(mx))
      dtbdt = (1.d0/dsubcld)* (mu(mx)*(shat(mx)-su(mx))+ &
                  md(mx)* (shat(mx)-sd(mx)))
      dqbdt = (1.d0/dsubcld)* (mu(mx)*(qhat(mx)-qu(mx))+ &
                 md(mx)* (qhat(mx)-qd(mx)))
      debdt = eps1*p(mx)/ ((eps1+q(mx))*(eps1+q(mx)))*dqbdt
      login1 = t(mx)
      call math_agent_1i1o(log_addr, login1, logout1)
      login2 = eb 
      call math_agent_1i1o(log_addr, login2, logout2)
      tmpres = 3.5d0*logout1-logout2-4.805d0
      dtldt = -2840.d0* (3.5d0/t(mx)*dtbdt-debdt/eb)/ &
                 (tmpres*tmpres)
!   end do
!
!   dtmdt and dqmdt are cumulus heating and drying.
!
   do k = msg + 1,pver
!      do i = il1g,il2g
         dtmdt(k) = 0.d0
         dqmdt(k) = 0.d0
!      end do
   end do
!
   do k = msg + 1,pver - 1
!      do i = il1g,il2g
         if (k == jt) then
            dtmdt(k) = (1.d0/dp(k))*(mu(k+1)* (su(k+1)-shat(k+1)- &
                          rl/cp*ql(k+1))+md(k+1)* (sd(k+1)-shat(k+1)))
            dqmdt(k) = (1.d0/dp(k))*(mu(k+1)* (qu(k+1)- &
                         qhat(k+1)+ql(k+1))+md(k+1)*(qd(k+1)-qhat(k+1)))
         end if
!      end do
   end do
!
   beta = 0.d0
   do k = msg + 1,pver - 1
!      do i = il1g,il2g
         if (k > jt .and. k < mx) then
            dtmdt(k) = (mc(k)* (shat(k)-s(k))+mc(k+1)* (s(k)-shat(k+1)))/ &
                         dp(k) - rl/cp*du(k)*(beta*ql(k)+ (1-beta)*ql(k+1))
!          dqmdt(k)=(mc(k)*(qhat(k)-q(k))
!     1                +mc(k+1)*(q(k)-qhat(k+1)))/dp(k)
!     2                +du(k)*(qs(k)-q(k))
!     3                +du(k)*(beta*ql(k)+(1-beta)*ql(k+1))

            dqmdt(k) = (mu(k+1)* (qu(k+1)-qhat(k+1)+cp/rl* (su(k+1)-s(k)))- &
                          mu(k)* (qu(k)-qhat(k)+cp/rl*(su(k)-s(k)))+md(k+1)* &
                         (qd(k+1)-qhat(k+1)+cp/rl*(sd(k+1)-s(k)))-md(k)* &
                         (qd(k)-qhat(k)+cp/rl*(sd(k)-s(k))))/dp(k) + &
                          du(k)* (beta*ql(k)+(1-beta)*ql(k+1))
         end if
!      end do
   end do
!
            powbase1 =  (1000.d0/p(k))
            powexp1 =   (rd/cp)
            call math_agent_2i1o(pow_addr, powbase1, powexp1, powout1)


   do k = msg + 1,pver
!      do i = il1g,il2g
         if (k >= lel .and. k <= lcl) then
            thetavp(k) = tp(k)* powout1 *(1.d0+1.608d0*qstp(k)-q(mx))

            !thetavp(k) = tp(k)* (1000.d0/p(k))** (rd/cp)*(1.d0+1.608d0*qstp(k)-q(mx))

            thetavm(k) = t(k)* powout1 *(1.d0+0.608d0*q(k))
            !thetavm(k) = t(k)* (1000.d0/p(k))** (rd/cp)*(1.d0+0.608d0*q(k))
            dqsdtp(k) = qstp(k)* (1.d0+qstp(k)/eps1)*eps1*rl/(rd*(tp(k)*tp(k)))
!
! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.
!
            dtpdt(k) = tp(k)/ (1.d0+rl/cp* (dqsdtp(k)-qstp(k)/tp(k)))* &
                        (dtbdt/t(mx)+rl/cp* (dqbdt/tl-q(mx)/ &
                         (tl*tl)*dtldt))
!
! dboydt is the integrand of cape change.
!
            dboydt(k) = ((dtpdt(k)/tp(k)+1.d0/(1.d0+1.608d0*qstp(k)-q(mx))* &
                          (1.608d0 * dqsdtp(k) * dtpdt(k) -dqbdt)) - (dtmdt(k)/t(k)+0.608d0/ &
                          (1.d0+0.608d0*q(k))*dqmdt(k)))*grav*thetavp(k)/thetavm(k)
         end if
!      end do
   end do
!
   do k = msg + 1,pver
!      do i = il1g,il2g
         if (k > lcl .and. k < mx) then
            thetavp(k) = tp(k)* powout1*(1.d0+0.608d0*q(mx))
            !thetavp(k) = tp(k)* (1000.d0/p(k))** (rd/cp)*(1.d0+0.608d0*q(mx))
            thetavm(k) = t(k)* powout1*(1.d0+0.608d0*q(k))
            !thetavm(k) = t(k)* (1000.d0/p(k))** (rd/cp)*(1.d0+0.608d0*q(k))
!
! dboydt is the integrand of cape change.
!
            dboydt(k) = (dtbdt/t(mx)+0.608d0/ (1.d0+0.608d0*q(mx))*dqbdt- &
                          dtmdt(k)/t(k)-0.608d0/ (1.d0+0.608d0*q(k))*dqmdt(k))* &
                          grav*thetavp(k)/thetavm(k)
         end if
!      end do
   end do

!
! buoyant energy change is set to 2/3*excess cape per 3 hours
!
!!!!!xukai need to fix 
   dadt  = 0.d0
   !kmin = minval(lel)
   !kmax = maxval(mx) - 1
   kmin = lel
   kmax = mx - 1
   do k = kmin, kmax
!      do i = il1g,il2g
         if ( k >= lel .and. k <= mx - 1) then
            dadt = dadt + dboydt(k)* (zf(k)-zf(k+1))
         endif
!      end do
   end do
!   do i = il1g,il2g
      dltaa = -1.d0* (cape-capelmt)
      if (dadt /= 0.d0) mb = max(dltaa/tau/dadt,0.d0)
!   end do
!
   return
end subroutine closure

subroutine q1q2_pjr(lchnk   , &
                    dqdt    ,dsdt    ,q       ,qs      ,qu      , &
                    su      ,du      ,qhat    ,shat    ,dp      , &
                    mu      ,md      ,sd      ,qd      ,ql      , &
                    dsubcld ,jt      ,mx      ,il1g    ,il2g    , &
                    cp      ,rl      ,msg     ,          &
                    dl      ,evp     ,cu      ,pver, &
                    log_addr, log10_addr, pow_addr)


   implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: phil rasch dec 19 1995
! 
!-----------------------------------------------------------------------


   integer, intent(in)  :: pver
   real(8), intent(in) :: cp

   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: il1g
   integer, intent(in) :: il2g
   integer, intent(in) :: msg

   real(8), intent(in) :: q(pver)
   real(8), intent(in) :: qs(pver)
   real(8), intent(in) :: qu(pver)
   real(8), intent(in) :: su(pver)
   real(8), intent(in) :: du(pver)
   real(8), intent(in) :: qhat(pver)
   real(8), intent(in) :: shat(pver)
   real(8), intent(in) :: dp(pver)
   real(8), intent(in) :: mu(pver)
   real(8), intent(in) :: md(pver)
   real(8), intent(in) :: sd(pver)
   real(8), intent(in) :: qd(pver)
   real(8), intent(in) :: ql(pver)
   real(8), intent(in) :: evp(pver)
   real(8), intent(in) :: cu(pver)
   real(8), intent(in) :: dsubcld

   real(8),intent(out) :: dqdt(pver),dsdt(pver)
   real(8),intent(out) :: dl(pver)
   integer kbm
   integer ktm
   integer jt
   integer mx
!
! work fields:
!
   integer i
   integer k

   real(8) emc
   real(8) rl
   integer(8), intent(in) :: log_addr, log10_addr, pow_addr
  !-------------------------------------------------------------------
   do k = msg + 1,pver
!      do i = il1g,il2g
         dsdt(k) = 0.d0
         dqdt(k) = 0.d0
         dl(k) = 0.d0
!      end do
   end do
!
! find the highest level top and bottom levels of convection
!
   ktm = pver
   kbm = pver
!   do i = il1g, il2g
      ktm = min(ktm,jt)
      kbm = min(kbm,mx)
!   end do

   do k = ktm,pver-1
!      do i = il1g,il2g
         emc = -cu (k)               &         ! condensation in updraft
               +evp(k)                         ! evaporating rain in downdraft

         dsdt(k) = -rl/cp*emc &
                     + (+mu(k+1)* (su(k+1)-shat(k+1)) &
                        -mu(k)*   (su(k)-shat(k)) &
                        +md(k+1)* (sd(k+1)-shat(k+1)) &
                        -md(k)*   (sd(k)-shat(k)) &
                       )/dp(k)

         dqdt(k) = emc + &
                    (+mu(k+1)* (qu(k+1)-qhat(k+1)) &
                     -mu(k)*   (qu(k)-qhat(k)) &
                     +md(k+1)* (qd(k+1)-qhat(k+1)) &
                     -md(k)*   (qd(k)-qhat(k)) &
                    )/dp(k)

         dl(k) = du(k)*ql(k+1)

!      end do
   end do

!
!DIR$ NOINTERCHANGE!
   do k = kbm,pver
!      do i = il1g,il2g
         if (k == mx) then
            dsdt(k) = (1.d0/dsubcld)* &
                        (-mu(k)* (su(k)-shat(k)) &
                         -md(k)* (sd(k)-shat(k)) &
                        )
            dqdt(k) = (1.d0/dsubcld)* &
                        (-mu(k)*(qu(k)-qhat(k)) &
                         -md(k)*(qd(k)-qhat(k)) &
                        )
         else if (k > mx) then
            dsdt(k) = dsdt(k-1)
            dqdt(k) = dqdt(k-1)
         end if
!      end do
   end do
!
   return
end subroutine q1q2_pjr

!!!#define pver 30
!        subroutine buoyan_dilute_compute(lchnk   ,ncol  ,   &
!                  q       ,t       ,p       ,z       ,pf    )
        subroutine buoyan_dilute_compute(lchnk   ,ncol         , &
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,rl      ,cape    , &
                  pblt    ,lcl     ,lel     ,lon     ,mx      , &
                  rd      ,grav    ,cp      ,msg     , &
                  tpert   , org    , landfrac, pver, zm_org, cpliq, &
                  cpwv, tfreez, eps1, rh2o, latice, log_addr, log10_addr, &
                  pow_addr &
                 )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
! 
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for 
!            testing (dmpdp). 
! 08/04/05 - Swap to convert dmpdz to dmpdp  
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled 
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
! 
! References:
! Raymond and Blythe (1992) JAS 
! 
! Author:
! Richard Neale - September 2004
! 
!-----------------------------------------------------------------------
        implicit none
!-----------------------------------------------------------------------
!
! input arguments
        integer, intent(in) :: pver
!      integer, parameter :: pver = 30        
!
       integer, intent(in) :: lchnk                 ! chunk identifier
       integer, intent(in) :: ncol                  ! number of atmospheric columns
    
       real(8), intent(in) :: q(pver)        ! spec. humidity
       real(8), intent(in) :: t(pver)        ! temperature
       real(8), intent(in) :: p(pver)        ! pressure
       real(8), intent(in) :: z(pver)        ! height
       real(8), intent(in) :: pf(pver+1)     ! pressure at interfaces
       real(8), intent(in) :: pblt          ! index of pbl depth
       real(8), intent(in) :: tpert         ! perturbation temperature by pbl processes

!
! outpu!t arguments
!
       real(8), intent(out) :: tp(pver)       ! parcel temperature
       real(8), intent(out) :: qstp(pver)     ! saturation mixing ratio of parcel (only above lcl, just q below).
       real(8), intent(out) :: tl            ! parcel temperature at lcl
       real(8), intent(out) :: cape          ! convective aval. pot. energy.
       integer lcl !
       integer lel !
       integer lon ! level of onset of deep convection
       integer mx         ! level of max moist static energy
    
       real(8), intent(inout) :: org(pver)      ! organization parameter
       real(8), intent(in) :: landfrac
       logical, intent(in) :: zm_org

!--------------------------Local Variables------------------------------
!
       real(8) capeten(5)     ! provisional value of cape
       real(8) tv(pver)       !
       real(8) tpv(pver)      !
       real(8) buoy(pver)
    
       real(8) a1
       real(8) a2
       real(8) estp
       real(8) pl
       real(8) plexp
       real(8) hmax
       real(8) hmn
       real(8) y
    
       logical plge600
       integer knt
       integer lelten(5)
    
       real(8), intent(in) :: cp
       real(8) e
       real(8) grav
    
       integer i
       integer k
       integer msg
       integer n
    
       real(8), intent(in)::  rd
       real(8), intent(in) :: rl
       real(8), intent(in) :: cpliq
       real(8), intent(in) :: cpwv 
       real(8), intent(in) :: tfreez 
       real(8), intent(in) :: eps1 
       real(8), intent(in) :: rh2o 
       real(8), intent(in) :: latice 
       real(8),parameter ::  tiedke_add = 0.5d0  
       !real(8) ::  tiedke_add  
       integer(8), intent(in) :: log_addr, log10_addr, pow_addr
       real(8) :: logout, login
!#ifdef PERGRO
!        real(8) rhd
!#endif
!
!-----------------------------------------------------------------------
!
      !tiedke_add = 0.5d0 
       do n = 1,5
             lelten(n) = pver
             capeten(n) = 0.d0
       end do

          lon = pver
          knt = 0
          lel = pver
          mx = lon
          cape = 0.d0
          hmax = 0.d0
   
       tp(:) = t(:)
       qstp(:) = q(:)

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
       tv(:) = t(:) *(1.0d0+1.608d0*q(:))/ (1.0d0+q(:))
       tpv(:) = tv(:)
       buoy(:) = 0.0d0


! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
!
!#ifdef PERGRO
!       do k = pver,msg + 1,-1
!         hmn = cp*t(k) + grav*z(k) + rl*q(k)
!!
!! Reset max moist static energy level when relative difference exceeds 1.e-4
!!
!         rhd = (hmn - hmax)/(hmn + hmax)
!         if (k >= nint(pblt) .and. k <= lon .and. rhd > -1.0d-4) then
!            hmax = hmn
!            mx = k
!         end if
!      end do
!#else
       do k = pver,msg + 1,-1
         hmn = cp*t(k) + grav*z(k) + rl*q(k)
         if (k >= nint(pblt) .and. k <= lon .and. hmn > hmax) then
            hmax = hmn
            mx = k
         end if
       end do
!#endif

! LCL dilute calculation - initialize to mx(i)
! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
! Original code actually sets LCL as level above wher condensate forms.
! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.

        ! Initialise LCL variables.
          lcl = mx
          tl = t(mx)
          pl = p(mx)
        !do k = 1, pver
        !    qtnd(k) = t(k)
        !enddo
!
!
! main buoyancy calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
!!   RBN 9/9/04   !!!
      !mcon(1) = lcl


       call parcel_dilute(lchnk, ncol, msg, mx, p, t, q, &
       tpert, tp, tpv, qstp, pl, tl, lcl, &
       org, landfrac, zm_org, pver, grav, rd, &
       rl, cpliq, cpwv, tfreez, eps1, rh2o, cp, latice, &
       log_addr, log10_addr,  pow_addr &
       )
      
     ! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
        plge600 = pl.ge.600.0d0 ! Just change to always allow buoy calculation.
!
!
! Main buoyancy calculation.
!
       do k = pver,msg + 1,-1
      !do i=1,ncol
         if (k <= mx .and. plge600) then   ! Define buoy from launch level to cloud top.
            tv(k) = t(k)* (1.0d0+1.608d0*q(k))/ (1.0d0+q(k))
            buoy(k) = tpv(k) - tv(k) + tiedke_add  ! +0.5K or not?

            !dlf(k) = buoy(k) 
         else
            qstp(k) = q(k)
            tp(k)   = t(k)            
            tpv(k)  = tv(k)

            !dlf(k) = qstp(k) 
            !heat(k) = tp(k)
    

         endif
        !qtnd(k) = tv(k)
      !end do
       end do
!

  

!-------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
       do k = msg + 2,pver
!      do i = 1,ncol
         if (k < lcl .and. plge600) then
            if (buoy(k+1) > 0.0d0 .and. buoy(k) <= 0.0d0) then
               knt = min(5,knt + 1)
               lelten(knt) = k
            end if
         end if
!      end do
       end do


      !tnd(1) = lcl + 0.d0
!
! calculate convective available potential energy (cape).
!
       !call math_agent_log_c(log_addr)
      !do k =  1,pver
      !      qtnd(k) = pf(k) 
      ! enddo
       do n = 1,5
       !qtnd(n) = lelten(n) + 3.d0
      do k = msg + 1,pver
         !do i = 1,ncol
            if (plge600 .and. k <= mx .and. k > lelten(n)) then
                login = pf(k+1)/pf(k)
                call math_agent_1i1o(log_addr, login, logout)
               capeten(n) = capeten(n) + rd*buoy(k)*logout
               !capeten(n) = capeten(n) + rd*buoy(k)*log(pf(k+1)/pf(k))
            !else
            !   qtnd(k) = 2.0d0
            end if
!        qtnd(k) = pf(k+1)/pf(k) 
         !end do
      end do
     end do
     !do i = 1, ncol

! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
       do n = 1,5
      !do i = 1,ncol
         if (capeten(n) > cape) then
            cape = capeten(n)
            lel = lelten(n)
         end if
        !qtnd(n) = capeten(n) 
      !end do
     end do
!
! put lower bound on cape for diagnostic purposes.
!
       !do i = 1,ncol
          cape = max(cape, 0.0d0)
       !end do
!
!        qtnd(1) = cape 
        !do k = 1, pver
            !qtnd(k) = buoy(k) 
        !enddo


       return
      end subroutine buoyan_dilute_compute

        subroutine parcel_dilute (lchnk, ncol, msg, klaunch, p, t, q, &
        tpert, tp, tpv, qstp, pl, tl, lcl, &
        org, landfrac, zm_org, pver, grav, rgas, rl, cpliq, cpwv, tfreez, &
        eps1, rh2o, cp, latice, log_addr, log10_addr, pow_addr &
        )

! Routine  to determine 
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

!--------------------
        implicit none
!--------------------
        integer, intent(in) :: pver 
        integer, intent(in) :: lchnk
        integer, intent(in) :: ncol
        integer, intent(in) :: msg
        
        integer, intent(in) :: klaunch
        
        real(8), intent(in), dimension(pver) :: p
        real(8), intent(in), dimension(pver) :: t
        real(8), intent(in), dimension(pver) :: q
        real(8), intent(in) :: tpert ! PBL temperature perturbation.
        
        real(8), intent(inout), dimension(pver) :: tp    ! Parcel temp.
        real(8), intent(inout), dimension(pver) :: qstp  ! Parcel water vapour (sat value above lcl).
        real(8), intent(inout) :: tl         ! Actual temp of LCL.
        real(8), intent(inout) :: pl          ! Actual pressure of LCL. 
        
        integer, intent(inout) :: lcl ! Lifting condesation level (first model level with saturation).
        
        real(8), intent(out), dimension(pver) :: tpv   ! Define tpv within this routine.
        
        real(8), intent(inout), dimension(pver) :: org
        real(8), intent(in) :: landfrac
!--------------------

! Have to be careful as s is also dry static energy.


! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.


        real(8) tmix(pver)        ! Tempertaure of the entraining parcel.
        real(8) qtmix(pver)       ! Total water of the entraining parcel.
        real(8) qsmix(pver)       ! Saturated mixing ratio at the tmix.
        real(8) smix(pver)        ! Entropy of the entraining parcel.
        real(8) xsh2o(pver)       ! Precipitate lost from parcel.
        real(8) ds_xsh2o(pver)    ! Entropy change due to loss of condensate.
        real(8) ds_freeze(pver)   ! Entropy change sue to freezing of precip.
        real(8) dmpdz2d(pver)     ! variable detrainment rate
        
        real(8) mp    ! Parcel mass flux.
        real(8) qtp   ! Parcel total water.
        real(8) sp    ! Parcel entropy.
        
        real(8) sp0    ! Parcel launch entropy.
        real(8) qtp0   ! Parcel launch total water.
        real(8) mp0    ! Parcel launch relative mass flux.
        
        real(8) lwmax      ! Maximum condesate that can be held in cloud before rainout.
        real(8) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
!real(r8) dmpdpc     ! In cloud parcel mass entrainment rate (/mb).
        real(8) dmpdz      ! Parcel fractional mass entrainment rate (/m)
        real(8) dpdz,dzdp  ! Hydrstatic relation and inverse of.
        real(8) senv       ! Environmental entropy at each grid point.
        real(8) qtenv      ! Environmental total water "   "   ".
        real(8) penv       ! Environmental total pressure "   "   ".
        real(8) tenv       ! Environmental total temperature "   "   ".
        real(8) new_s      ! Hold value for entropy after condensation/freezing adjustments.
        real(8) new_q      ! Hold value for total water after condensation/freezing adjustments.
        real(8) dp         ! Layer thickness (center to center)
        real(8) tfguess    ! First guess for entropy inversion - crucial for efficiency!
        real(8) tscool     ! Super cooled temperature offset (in degC) (eg -35).

        real(8) qxsk, qxskp1        ! LCL excess water (k, k+1)
        real(8) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
        real(8) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.
        real(8) org2rkm, org2Tpert
        real(8) dmpdz_lnd, dmpdz_mask
        
        integer rcall       ! Number of ientropy call for errors recording
        integer nit_lheat     ! Number of iterations for condensation/freezing loop.

        integer i,k,ii   ! Loop counters.
!+++++++++++++++
        real(8), intent(in) :: grav
        logical, intent(in) :: zm_org
        real(8), intent(in) :: rgas
        real(8), intent(in) :: rl 
        real(8), intent(in) :: cpliq 
        real(8), intent(in) :: cpwv
        real(8), intent(in) :: tfreez
        real(8), intent(in) :: eps1 
        real(8), intent(in) :: rh2o
        real(8), intent(in) :: cp 
        real(8), intent(in) :: latice 
       integer(8), intent(in) :: log_addr, log10_addr, pow_addr
       real(8) :: logout, login
!-----------------

!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
!
!======================================================================
!
! Set some values that may be changed frequently.
!

        if (zm_org) then
           org2rkm = 10.0d0
           org2Tpert = 0.0d0
        endif

        nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.
        dmpdz=-1.0d-3       ! Entrainment rate. (-ve for /m)
        dmpdz_lnd=-1.0d-3
        !dmpdpc = 3.e-2_r8   ! In cloud entrainment rate (/mb).
        lwmax = 1.0d-3    ! Need to put formula in for this.
        tscool = 0.0d0   ! Temp at which water loading freezes in the cloud.
        
        qtmix=0.0d0
        smix=0.0d0
        
        qtenv = 0.0d0
        senv = 0.0d0
        tenv = 0.0d0
        penv = 0.0d0
        
        qtp0 = 0.0d0
        sp0  = 0.0d0
        mp0 = 0.0d0
        
        qtp = 0.0d0
        sp = 0.0d0
        mp = 0.0d0
        
        new_q = 0.0d0
        new_s = 0.0d0

! **** Begin loops ****

        do k = pver, msg+1, -1
        !do i=1,ncol 

! Initialize parcel values at launch level.

      if (k == klaunch) then 
         qtp0 = q(k)   ! Parcel launch total water (assuming subsaturated) - OK????.
         sp0  = entropy(t(k),p(k),qtp0,rl,cpliq, cpwv, tfreez,&
                        eps1, rh2o, cp, rgas, log_addr, log10_addr, pow_addr)  ! Parcel launch entropy.
         mp0  = 1.0d0       

        ! Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute). 
         smix(k)  = sp0
         qtmix(k) = qtp0
         tfguess = t(k)
         rcall = 1
         ! i need to determine
         call ientropy (rcall,i,lchnk,smix(k),p(k),qtmix(k),tmix(k),qsmix(k),tfguess,&
                rl,cpliq, cpwv, tfreez, eps1, rh2o, cp, rgas, log_addr, log10_addr, pow_addr)

      end if

! Entraining levels
      
      if (k < klaunch) then 

! Set environmental values for this level.                 
         
         dp = (p(k)-p(k+1)) ! In -ve mb as p decreasing with height - difference between center of layers.
         qtenv = 0.5d0*(q(k)+q(k+1))         ! Total water of environment.
         tenv  = 0.5d0*(t(k)+t(k+1)) 
         penv  = 0.5d0*(p(k)+p(k+1))

         senv  = entropy(tenv,penv,qtenv,rl, cpliq, cpwv, tfreez,&
                            eps1, rh2o, cp, rgas, log_addr, log10_addr, pow_addr)  ! Entropy of environment.   


! Determine fractional entrainment rate /pa given value /m.

         dpdz = -(penv*grav)/(rgas*tenv) ! in mb/m since  p in mb.
         dzdp = 1.0d0/dpdz                  ! in m/mb
         if (zm_org) then
            dmpdz_mask = landfrac * dmpdz_lnd + (1.0d0 - landfrac) * dmpdz
            dmpdp = (dmpdz_mask/(1.0d0+org(k)*org2rkm))*dzdp              ! /mb Fractional entrainment
         else
            dmpdp = dmpdz*dzdp
         endif

! Sum entrainment to current level
! entrains q,s out of intervening dp layers, in which linear variation is assumed
! so really it entrains the mean of the 2 stored values.

         sp  = sp  - dmpdp*dp*senv 
         qtp = qtp - dmpdp*dp*qtenv 
         mp  = mp - dmpdp*dp
            
! Entrain s and qt to next level.

         smix(k)  = (sp0  +  sp) / (mp0 + mp)

        !kai
        !qtnd(k) = smix(k) 
        !qtnd(k) = dmpdp 


         qtmix(k) = (qtp0 + qtp) / (mp0 + mp)

! Invert entropy from s and q to determine T and saturation-capped q of mixture.
! t(i,k) used as a first guess so that it converges faster.

         tfguess = tmix(k+1)
         rcall = 2
         ! i need to determine 
         call ientropy(rcall,i,lchnk,smix(k),p(k),qtmix(k),tmix(k),qsmix(k),tfguess,&
                rl,cpliq, cpwv, tfreez, eps1, rh2o, cp, rgas, log_addr, log10_addr, pow_addr)   
!
! Determine if this is lcl of this column if qsmix <= qtmix.
! FIRST LEVEL where this happens on ascending.

         if (qsmix(k) <= qtmix(k) .and. qsmix(k+1) > qtmix(k+1)) then
            lcl = k
            qxsk   = qtmix(k) - qsmix(k)
            qxskp1 = qtmix(k+1) - qsmix(k+1)
            dqxsdp = (qxsk - qxskp1)/dp
            pl  = p(k+1) - qxskp1/dqxsdp    ! pressure level of actual lcl.
            dsdp   = (smix(k)  - smix(k+1))/dp
            dqtdp  = (qtmix(k) - qtmix(k+1))/dp
            slcl   = smix(k+1)  +  dsdp* (pl-p(k+1))  
            qtlcl  = qtmix(k+1) +  dqtdp*(pl-p(k+1))

            tfguess = tmix(k)
            rcall = 3
            ! i need to determine, need to fix
            call ientropy (rcall,i,lchnk,slcl,pl,qtlcl,tl,qslcl,tfguess,&
                rl,cpliq, cpwv, tfreez, eps1, rh2o, cp, rgas, log_addr, log10_addr, pow_addr)

!            write(iulog,*)' '
!            write(iulog,*)' p',p(i,k+1),pl(i),p(i,lcl(i))
!            write(iulog,*)' t',tmix(i,k+1),tl(i),tmix(i,lcl(i))
!            write(iulog,*)' s',smix(i,k+1),slcl,smix(i,lcl(i))
!            write(iulog,*)'qt',qtmix(i,k+1),qtlcl,qtmix(i,lcl(i))
!            write(iulog,*)'qs',qsmix(i,k+1),qslcl,qsmix(i,lcl(i))

         endif
!         
      end if !  k < klaunch

 
!   end do ! Levels loop
        end do ! Columns loop


        !do k = 1, klaunch
        !do k = 1, pver
            !qtnd(k) = klaunch !tmix(k)
            !qtnd(k) = smix(k)
        !enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the 
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or 
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
!! provides latent heating to the mixed parcel and so this has to be added back 
!! to it. But does this also increase qsmix as well? Also freezing processes
 

        xsh2o = 0.0d0
        ds_xsh2o = 0.0d0
        ds_freeze = 0.0d0
        
!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
! Iterate solution twice for accuracy

       !call math_agent_log_c(log_addr)


        do k = pver, msg+1, -1
!        do i=1,ncol    
      
! Initialize variables at k=klaunch
      
      if (k == klaunch) then

! Set parcel values at launch level assume no liquid water.            

         tp(k)    = tmix(k)
         qstp(k)  = q(k) 
         if (zm_org) then
            tpv(k)   =  (tp(k) + (org2Tpert*org(k)+tpert)) * (1.0d0+1.608d0*qstp(k)) / (1.0d0+qstp(k))
         else
            tpv(k)   =  (tp(k) + tpert) * (1.d0+1.608d0*qstp(k)) / (1.d0+qstp(k))
         endif
         
      end if

      if (k < klaunch) then
            
! Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.

! Iterate nit_lheat times for s,qt changes.

         do ii=0,nit_lheat-1            

! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).

            xsh2o(k) = max (0.0d0, qtmix(k) - qsmix(k) - lwmax)

! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)                     
             login =  tmix(k)/tfreez         
            call math_agent_1i1o(log_addr, login, logout)
            ds_xsh2o(k) = ds_xsh2o(k+1) - cpliq * logout * max(0.0d0,(xsh2o(k)-xsh2o(k+1)))
            !ds_xsh2o(k) = ds_xsh2o(k+1) - cpliq * log (tmix(k)/tfreez) * max(0.0d0,(xsh2o(k)-xsh2o(k+1)))
!
! Entropy of freezing: latice times amount of water involved divided by T.
!
 
            if (tmix(k) <= tfreez+tscool .and. ds_freeze(k+1) == 0.0d0) then ! One off freezing of condensate. 
               ds_freeze(k) = (latice/tmix(k)) * max(0.0d0,qtmix(k)-qsmix(k)-xsh2o(k)) ! Gain of LH
            end if
            
            if (tmix(k) <= tfreez+tscool .and. ds_freeze(k+1) /= 0.d0) then ! Continual freezing of additional condensate.
               ds_freeze(k) = ds_freeze(k+1)+(latice/tmix(k)) * max(0.0d0,(qsmix(k+1)-qsmix(k)))
            end if
            
! Adjust entropy and accordingly to sum of ds (be careful of signs).

            new_s = smix(k) + ds_xsh2o(k) + ds_freeze(k) 

! Adjust liquid water and accordingly to xsh2o.

            new_q = qtmix(k) - xsh2o(k)

! Invert entropy to get updated Tmix and qsmix of parcel.

            tfguess = tmix(k)
            rcall =4
            call ientropy (rcall,i,lchnk,new_s, p(k), new_q, tmix(k), qsmix(k), tfguess, &
                rl,cpliq, cpwv, tfreez, eps1, rh2o, cp, rgas, log_addr, log10_addr, pow_addr)
            
         end do  ! Iteration loop for freezing processes.

! tp  - Parcel temp is temp of mixture.
! tpv - Parcel v. temp should be density temp with new_q total water. 

         tp(k)    = tmix(k)

! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)

         if (new_q > qsmix(k)) then  ! Super-saturated so condensate present - reduces buoyancy.
            qstp(k) = qsmix(k)
         else                          ! Just saturated/sub-saturated - no condensate virtual effects.
            qstp(k) = new_q
         end if

         if (zm_org) then
            tpv(k) = (tp(k)+(org2Tpert*org(k)+tpert))* (1.0d0+1.608d0*qstp(k)) / (1.0d0+ new_q) 
         else
            tpv(k) = (tp(k)+tpert)* (1.0d0+1.608d0*qstp(k)) / (1.0d0+ new_q) 
         endif

      end if ! k < klaunch
      
!end do ! Loop for columns
   
        end do  ! Loop for vertical levels.

return
        end subroutine parcel_dilute
!-----------------------------------------------------------------------------------------
        real(8) function entropy(TK,p,qtot, rl, cpliq, cpwv, tfreez, &
        eps1, rh2o, cpres, rgas, log_addr, log10_addr, pow_addr)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
!
         implicit none
     real(8), intent(in) :: p,qtot,TK, rl, cpliq, &
     cpwv, tfreez, eps1,  rh2o, cpres, rgas
     real(8) :: qv,qst,e,est,L
     real(8), parameter :: pref = 1000.0d0
       integer(8), intent(in) :: log_addr,log10_addr, pow_addr
       real(8) :: logout1,logout2, logout3, login

        L = rl - (cpliq - cpwv)*(TK-tfreez)         ! T IN CENTIGRADE

        call qsat_hPa(TK, p, est, qst, tfreez, eps1, rl, rh2o, cpres, log_addr,log10_addr, pow_addr)


        qv = min(qtot,qst)                         ! Partition qtot into vapor part only.
        e = qv*p / (eps1 +qv)

       !call math_agent_log_c(log_addr)

             login =  TK/tfreez
            call math_agent_1i1o(log_addr, login, logout1)
             login = (p-e)/pref 
            call math_agent_1i1o(log_addr, login, logout2)
             login = qv/qst
            call math_agent_1i1o(log_addr, login, logout3)
        !entropy = (cpres + qtot*cpliq)*( TK/tfreez) - rgas*( (p-e)/pref ) + &
        !L*qv/TK - qv*rh2o*(qv/qst)
        entropy = (cpres + qtot*cpliq)*logout1 - rgas*logout2 + &
        L*qv/TK - qv*rh2o*logout3
        !entropy = (cpres + qtot*cpliq)*log( TK/tfreez) - rgas*log( (p-e)/pref ) + &
        !L*qv/TK - qv*rh2o*log(qv/qst)

        end FUNCTION entropy

!
!-----------------------------------------------------------------------------------------
        SUBROUTINE ientropy (rcall,icol,lchnk,s,p,qt,T,qst,Tfg, rl, cpliq, cpwv, tfreez, &
        eps1, rh2o, cpres, rgas, log_addr,log10_addr, pow_addr)
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
! 

         implicit none
      !use phys_grid, only: get_rlon_p, get_rlat_p
    
      integer, intent(in) :: icol, lchnk, rcall
      real(8), intent(in)  :: s, p, Tfg, qt
      real(8), intent(in) :: rl, cpliq, cpwv, tfreez, &
      eps1, rh2o, cpres, rgas
      real(8), intent(out) :: qst, T
      real(8) :: est, this_lat,this_lon
      real(8) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm,tol
      integer :: i
    
      logical :: converged

       integer(8), intent(in) :: log_addr,log10_addr, pow_addr
! Max number of iteration loops.
        integer, parameter :: LOOPMAX = 100
        real(8), parameter :: EPS = 3.0d-8

        converged = .false.

! Invert the entropy equation -- use Brent's method
! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.

      T = Tfg                  ! Better first guess based on Tprofile from conv.
    
      a = Tfg-10    !low bracket
      b = Tfg+10    !high bracket
    
      fa = entropy(a, p, qt, rl, cpliq, cpwv, &
                     tfreez, eps1, rh2o, cpres, rgas, log_addr,log10_addr, pow_addr) - s
      fb = entropy(b, p, qt, rl, cpliq, cpwv, &
                     tfreez, eps1, rh2o, cpres, rgas, log_addr,log10_addr, pow_addr) - s
    
      c=b
      fc=fb
      tol=0.001d0
    
      converge: do i=0, LOOPMAX
         if ((fb > 0.0d0 .and. fc > 0.0d0) .or. &
              (fb < 0.0d0 .and. fc < 0.0d0)) then
            c=a
            fc=fa
            d=b-a
            ebr=d
         end if
         if (abs(fc) < abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         end if

         tol1=2.0d0*EPS*abs(b)+0.5d0*tol
         xm=0.5d0*(c-b)
         converged = (abs(xm) <= tol1 .or. fb == 0.0d0)
         if (converged) exit converge

         if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
            sbr=fb/fa
            if (a == c) then
               pbr=2.0d0*xm*sbr
               qbr=1.0d0-sbr
            else
               qbr=fa/fc
               rbr=fb/fc
               pbr=sbr*(2.0d0*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0d0))
               qbr=(qbr-1.0d0)*(rbr-1.0d0)*(sbr-1.0d0)
            end if
            if (pbr > 0.0d0) qbr=-qbr
            pbr=abs(pbr)
            if (2.0d0*pbr  <  min(3.0d0*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
               ebr=d
               d=pbr/qbr
            else
               d=xm
               ebr=d
            end if
         else
            d=xm
            ebr=d
         end if
         a=b
         fa=fb
         b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )

         fb = entropy(b, p, qt, rl, cpliq, cpwv, &
                     tfreez, eps1, rh2o, cpres, rgas, log_addr,log10_addr, pow_addr) - s

        end do converge

      T = b
      call qsat_hPa(T, p, est, qst, tfreez, eps1, rl, rh2o, cpres, log_addr,log10_addr, pow_addr)
    
!      if (.not. converged) then
!     this_lat = get_rlat_p(lchnk, icol)*57.296_r8
!     this_lon = get_rlon_p(lchnk, icol)*57.296_r8
!     write(iulog,*) '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
!     write(iulog,100) 'ZM_CONV: IENTROPY. Details: call#,lchnk,icol= ',rcall,lchnk,icol, &
!          ' lat: ',this_lat,' lon: ',this_lon, &
!          ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000._r8*qt, &
!          ' qst(g/kg) = ', 1000._r8*qst,', s(J/kg) = ',s
!     call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
!      end if
!
!    100 format (A,I1,I4,I4,7(A,F6.2))

        End SUBROUTINE ientropy

! Wrapper for qsat_water that does translation between Pa and hPa
! qsat_water uses Pa internally, so get it right, need to pass in Pa.
! Afterward, set es back to hPa.
         subroutine qsat_hPa(t, p, es, qm, tfreez, epsilo, latvap, rh2o, cpair, log_addr,log10_addr, pow_addr)
        !use wv_saturation, only: qsat_water

         implicit none
! Inputs
      real(8), intent(in) :: t    ! Temperature (K)
      real(8), intent(in) :: p    ! Pressure (hPa)
      ! Outputs
      real(8), intent(out) :: es  ! Saturation vapor pressure (hPa)
      real(8), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass, kg/kg)

        real(8), intent(in) :: tfreez
        real(8), intent(in) :: epsilo 
        real(8), intent(in) :: latvap
        real(8), intent(in) :: rh2o 
        real(8), intent(in) :: cpair 
        real(8) :: p_t
        integer(8), intent(in) :: log_addr,log10_addr, pow_addr
        !call qsat_water(t, p*100.0d0, es, qm, tfreez, epsilo, latvap, rh2o, cpair)
        p_t = p * 100.0d0
        call wv_sat_qsat_water(t, p_t, es, qm, tfreez, epsilo, log_addr,log10_addr, pow_addr)
        es = es*0.01d0

        End subroutine qsat_hPa
!         subroutine qsat_water(t, p, es, qs, tfreez, epsilo,latvap, rh2o, cpair, gam, dqsdt, enthalpy)
!!------------------------------------------------------------------!
!! Purpose:                                                         !
!!   Calculate SVP over water at a given temperature, and then      !
!!   calculate and return saturation specific humidity.             !
!!   Optionally return various temperature derivatives or enthalpy  !
!!   at saturation.                                                 !
!!------------------------------------------------------------------!
!
!!  use wv_sat_methods, only: wv_sat_qsat_water
!
!         implicit none
!! Inputs
!        real(8), intent(in) :: t    ! Temperature
!        real(8), intent(in) :: p    ! Pressure
!! Outputs
!        real(8), intent(out) :: es  ! Saturation vapor pressure
!        real(8), intent(out) :: qs  ! Saturation specific humidity
!
!        real(8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
!        real(8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
!        real(8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q
!
!! Local variables
!        real(8) :: hltalt       ! Modified latent heat for T derivatives
!        real(8), intent(in) :: tfreez
!        real(8), intent(in) :: epsilo 
!        real(8), intent(in) :: latvap 
!        real(8), intent(in) :: rh2o 
!        real(8), intent(in) :: cpair 
!
!        call wv_sat_qsat_water(t, p, es, qs, tfreez, epsilo)
!
!       ! if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then
!
!       !  ! "generalized" analytic expression for t derivative of es
!       !   ! accurate to within 1 percent for 173.16 < t < 373.16
!       !  call no_ip_hltalt(t, hltalt, tfreez, latvap)
!
!       ! !if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt,cpair)
!       ! if (present(enthalpy))  enthalpy = cpair * t + hltalt * qs
!
!       ! ! For pure water/ice transition term is 0.
!       ! call deriv_outputs(t, p, es, qs, hltalt, 0.0d0, epsilo, rh2o, cpair, &
!       !   gam=gam, dqsdt=dqsdt)
!
!       ! end if
!
!        end subroutine qsat_water
         subroutine wv_sat_qsat_water(t, p, es, qs, tmelt, epsilo, log_addr,log10_addr, pow_addr, idx)
          !------------------------------------------------------------------!
          ! Purpose:                                                         !
          !   Calculate SVP over water at a given temperature, and then      !
          !   calculate and return saturation specific humidity.             !
          !------------------------------------------------------------------!
        
         implicit none
          ! Inputs
          real(8), intent(in) :: t    ! Temperature
          real(8), intent(in) :: p    ! Pressure
          ! Outputs
          real(8), intent(out) :: es  ! Saturation vapor pressure
          real(8), intent(out) :: qs  ! Saturation specific humidity
        
          integer,  intent(in), optional :: idx ! Scheme index
          real(8), intent(in) :: tmelt   
          real(8), intent(in) :: epsilo   
        
            !real(8), parameter :: tboil = 373.16d0
            real(8) :: tboil
          !  integer(kind=8) :: log10_addr, pow_addr
          real(8) :: logout1, logout2, login, &
                    powbase, powbase2, powbase3, &
                    powexp1, powexp2, powexp3, &
                    powout1, powout2, powout3

        integer(8), intent(in) :: log_addr,log10_addr, pow_addr
        !integer(8) :: abs_addr
                  tboil = 373.16d0

          !es = wv_sat_svp_water(t, tmelt, idx)
          !es = GoffGratch_svp_water(t)
         !!es = 10.0d0**(-7.90298d0*(tboil/t-1.0d0)+ &
         !      5.02808d0*(tboil/t)- &
         !      1.3816d-7*(10.0d0**(11.344d0*(1.0d0-t/tboil))-1.0d0)+ &
         !      8.1328d-3*(10.0d0**(-3.49149d0*(tboil/t-1.0d0))-1.0d0)+ &
         !      (1013.246d0))*100.0d0

          !call  math_agent_log10_c(log10_addr)

          login  = tboil/t
          !call math_agent_abs_c(abs_addr)
          !!if(login < 0.d0) login = 0.d0  - login
          call math_agent_1i1o(log10_addr, login, logout1)
          login  = 1013.246d0 
          call math_agent_1i1o(log10_addr, login, logout2)

          !logout1 = 1.0d0 
          !logout2  = logout2 + logout1
          !call math_agent_pow_c(pow_addr)
          powbase = 10.0d0
          powexp1 = (11.344d0*(1.d0-t/tboil))
          call math_agent_2i1o(pow_addr, powbase, powexp1, powout1)
          powexp2 = (-3.49149d0*(tboil/t - 1.d0))
          call math_agent_2i1o(pow_addr, powbase, powexp2, powout2)
           powexp3 = (-7.90298d0*(tboil/t-1.d0)+ &
            5.02808d0*logout1 - &
            1.3816d-7*(powout1-1.d0)+ &
            8.1328d-3*(powout2-1.d0)+ &
            logout2)

          call math_agent_2i1o(pow_addr, powbase, powexp3, powout3)
          es = powout3 * 100.d0 
            !if(es > 0.0d0 ) es = 1.0d-16
            !if(es < 0.0d0 ) es = 1.0d-15
            !es = 1.0d-16


!          es = 10.0d0**(-7.90298d0*(tboil/t-1.0d0)+ &
!               5.02808d0*logout1- &
!               1.3816d-7*(10.0d0**(11.344d0*(1.0d0-t/tboil))-1.0d0)+ &
!               8.1328d-3*(10.0d0**(-3.49149d0*(tboil/t-1.0d0))-1.0d0)+ &
!               logout2)*100.0d0
!

          qs = wv_sat_svp_to_qsat(es, p, epsilo)
        
          ! Ensures returned es is consistent with limiters on qs.
          es = min(es, p)
        
        end subroutine wv_sat_qsat_water


         !function wv_sat_svp_water(t, tmelt, idx) result(es)
         !implicit none
         ! real(8), intent(in) :: t
         ! real(8), intent(in) :: tmelt
         ! integer,  intent(in), optional :: idx
         ! real(8) :: es
        
         ! integer :: use_idx
        
         ! integer, parameter :: OldGoffGratch_idx = 0
         ! integer, parameter :: GoffGratch_idx = 1
         ! integer, parameter :: MurphyKoop_idx = 2
         ! integer, parameter :: Bolton_idx = 3

         ! if (present(idx)) then
         !    use_idx = idx
         ! else
         !    use_idx = 1 
         ! end if
        
         ! es = GoffGratch_svp_water(t)
         ! !select case (use_idx)
         ! !case(GoffGratch_idx)
         ! !   es = GoffGratch_svp_water(t)
         ! !case(MurphyKoop_idx)
         ! !   es = MurphyKoop_svp_water(t)
         ! !case(OldGoffGratch_idx)
         ! !   es = OldGoffGratch_svp_water(t)
         ! !case(Bolton_idx)
         ! !   es = Bolton_svp_water(t, tmelt)
         ! !end select
        
         !end function wv_sat_svp_water
         function wv_sat_svp_to_qsat(es, p, epsilo) result(qs)
        
         implicit none
          real(8), intent(in) :: es  ! SVP
          real(8), intent(in) :: p   ! Current pressure.
          real(8) :: qs
        
          real(8), intent(in) :: epsilo  
          real(8) :: omeps
          omeps = 1.0d0 - epsilo

          ! If pressure is less than SVP, set qs to maximum of 1.
          if ( (p - es) <= 0.0d0 ) then
             qs = 1.0d0
          else
             qs = epsilo*es / (p - omeps*es)
          end if
        
        end function wv_sat_svp_to_qsat


         !subroutine no_ip_hltalt(t, hltalt, tmelt, latvap)
         ! !------------------------------------------------------------------!
         ! ! Purpose:                                                         !
         ! !   Calculate latent heat of vaporization of pure liquid water at  !
         ! !   a given temperature.                                           !
         ! !------------------------------------------------------------------!
         ! implicit none       
         ! ! Inputs
         ! real(8), intent(in) :: t        ! Temperature
         ! ! Outputs
         ! real(8), intent(out) :: hltalt  ! Appropriately modified hlat
         ! real(8), intent(in) :: tmelt
         ! real(8), intent(in) :: latvap 
        
         ! hltalt = latvap
        
         ! ! Account for change of latvap with t above freezing where
         ! ! constant slope is given by -2369 j/(kg c) = cpv - cw
         ! if (t >= tmelt) then
         !    hltalt = hltalt - 2369.0d0*(t-tmelt)
         ! end if

        !end subroutine no_ip_hltalt



        !! Temperature derivative outputs, for qsat_*
         !subroutine deriv_outputs(t, p, es, qs, hltalt, tterm, epsilo, rh2o, cpair,&
         !    gam, dqsdt)
        
         ! implicit none
         ! ! Inputs
         ! real(8), intent(in) :: t      ! Temperature
         ! real(8), intent(in) :: p      ! Pressure
         ! real(8), intent(in) :: es     ! Saturation vapor pressure
         ! real(8), intent(in) :: qs     ! Saturation specific humidity
         ! real(8), intent(in) :: hltalt ! Modified latent heat
         ! real(8), intent(in) :: tterm  ! Extra term for d(es)/dT in
         !                                ! transition region.
        
         ! ! Outputs
         ! real(8), intent(out), optional :: gam      ! (hltalt/cpair)*(d(qs)/dt)
         ! real(8), intent(out), optional :: dqsdt    ! (d(qs)/dt)
        
         ! ! Local variables
         ! real(8) :: desdt        ! d(es)/dt
         ! real(8) :: dqsdt_loc    ! local copy of dqsdt
        
         ! real(8), intent(in) :: epsilo   
         ! real(8), intent(in) :: rh2o   
         ! real(8), intent(in) :: cpair   
         ! real(8) :: omeps
         ! omeps = 1.0d0 - epsilo
         ! if (qs == 1.0d0) then
         !    dqsdt_loc = 0.0d0
         ! else
         !    desdt = hltalt*es/(rh2o*t*t) + tterm
         !    dqsdt_loc = qs*p*desdt/(es*(p-omeps*es))
         ! end if
        
         ! if (present(dqsdt)) dqsdt = dqsdt_loc
         ! if (present(gam))   gam   = dqsdt_loc * (hltalt/cpair)
        
         !end subroutine deriv_outputs
         !function tq_enthalpy(t, q, hltalt, cpair) result(enthalpy)
         ! implicit none
        
         ! real(8), intent(in) :: t      ! Temperature
         ! real(8), intent(in) :: q      ! Specific humidity
         ! real(8), intent(in) :: hltalt ! Modified hlat for T derivatives
         ! real(8), intent(in) :: cpair
        
         ! real(8) :: enthalpy
        
         ! enthalpy = cpair * t + hltalt * q
         ! 
         !end function tq_enthalpy

         !elemental function GoffGratch_svp_water(t) result(es)
         ! implicit none
         ! real(8), intent(in) :: t  ! Temperature in Kelvin
         ! real(8) :: es             ! SVP in Pa
        
         !   real(8), parameter :: tboil = 373.16d0
         !   integer(kind=8) :: log10_addr
         ! real(8) :: logout1, logout2, logout3, login
         ! ! uncertain below -70 C
         !es = 10.0d0**(-7.90298d0*(tboil/t-1.0d0)+ &
         !      5.02808d0*(tboil/t)- &
         !      1.3816d-7*(10.0d0**(11.344d0*(1.0d0-t/tboil))-1.0d0)+ &
         !      8.1328d-3*(10.0d0**(-3.49149d0*(tboil/t-1.0d0))-1.0d0)+ &
         !      (1013.246d0))*100.0d0

         ! !call  math_agent_log10_c(log10_addr)

         ! !login  = tboil/t
         ! !call math_agent_1i1o(log10_addr, login, logout1)

         ! !login  = 1013.246d0 
         ! !call math_agent_1i1o(log10_addr, login, logout2)
  
         ! !es = 10.0d0**(-7.90298d0*(tboil/t-1.0d0)+ &
         ! !     5.02808d0*logout1- &
         ! !     1.3816d-7*(10.0d0**(11.344d0*(1.0d0-t/tboil))-1.0d0)+ &
         ! !     8.1328d-3*(10.0d0**(-3.49149d0*(tboil/t-1.0d0))-1.0d0)+ &
         ! !     logout2)*100.0d0
!        !   es = 10.0d0**(-7.90298d0*(tboil/t-1.0d0)+ &
!        !       5.02808d0*log10(tboil/t)- &
!        !       1.3816d-7*(10.0d0**(11.344d0*(1.0d0-t/tboil))-1.0d0)+ &
!        !       8.1328d-3*(10.0d0**(-3.49149d0*(tboil/t-1.0d0))-1.0d0)+ &
!        !       log10(1013.246d0))*100.0d0
        
         !end function GoffGratch_svp_water

! Murphy & Koop (2005)

!         function MurphyKoop_svp_water(t) result(es)
!          implicit none
!          real(8), intent(in) :: t  ! Temperature in Kelvin
!          real(8) :: es             ! SVP in Pa
!        
!          ! (good for 123 < T < 332 K)
!          es = exp(54.842763d0 - (6763.22d0 / t) - (4.210d0 * log(t)) + &
!               (0.000367d0 * t) + (tanh(0.0415d0 * (t - 218.8d0)) * &
!               (53.878d0 - (1331.22d0 / t) - (9.44523d0 * log(t)) + &
!               0.014025d0 * t)))
!        
!        end function MurphyKoop_svp_water
!
!
!         function OldGoffGratch_svp_water(t) result(es)
!          implicit none
!          real(8), intent(in) :: t
!          real(8) :: es
!          real(8) :: ps, e1, e2, f1, f2, f3, f4, f5, f
!        
!        real(8), parameter :: tboil = 373.16d0
!          ps = 1013.246d0
!          e1 = 11.344d0*(1.0d0 - t/tboil)
!          e2 = -3.49149d0*(tboil/t - 1.0d0)
!          f1 = -7.90298d0*(tboil/t - 1.0d0)
!          f2 = 5.02808d0*log10(tboil/t)
!          f3 = -1.3816d0*(10.0d0**e1 - 1.0d0)/10000000.0d0
!          f4 = 8.1328d0*(10.0d0**e2 - 1.0d0)/1000.0d0
!          f5 = log10(ps)
!          f  = f1 + f2 + f3 + f4 + f5
!        
!          es = (10.0d0**f)*100.0d0
!          
!        end function OldGoffGratch_svp_water
!
!         function Bolton_svp_water(t, tmelt) result(es)
!          implicit none
!          real(8),parameter :: c1 = 611.2d0
!          real(8),parameter :: c2 = 17.67d0
!          real(8),parameter :: c3 = 243.5d0
!        
!          real(8), intent(in) :: t  ! Temperature in Kelvin
!          real(8) :: es             ! SVP in Pa
!        
!          real(8), intent(in) :: tmelt
!          es = c1*exp( (c2*(t - tmelt))/((t - tmelt)+c3) )
!        
!        end function Bolton_svp_water



        end module 
