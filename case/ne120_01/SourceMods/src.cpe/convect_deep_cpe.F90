!*************************************************************************
!File Name: src.cpe/convect_deep_cpe.f90
!Author: Xu Kai
!Created Time: 2018年11月07日 星期三 15时36分12秒
! ************************************************************************/

      subroutine zm_convr_compute(lchnk   ,ncol    , &
      t       ,qh      ,prec    ,jctop   ,jcbot   , &
      pblh    ,zm      ,geos    ,zi      ,qtnd    , &
      heat    ,pap     ,paph    ,dpp     , &
      delt    ,mcon    ,cme     ,cape    , &
      tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
      mu      ,md      ,du      ,eu      ,ed      , &
      dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
      lengath ,ql      ,rliq    ,landfrac,          &
      pver    ,pcols, pverp, zm_org_int, limcnv, org, orgt,&
      org2d   ,grav    ,rgrav   ,cpres, rl, rgas, cpliq, &
      cpwv, tfreez, eps1, rh2o, latice, c0_ocn, c0_lnd, no_deep_pbl_in,&
      tttt)
      use buoyan_dilute_cpe, only : buoyan_dilute_compute, cldprp, &
      closure, q1q2_pjr
      !      use buoyan_dilute_cpe, only : buoyan_dilute_compute

      implicit none
      integer, intent(in) :: lchnk    
      integer, intent(in) :: ncol        

      integer, intent(in) :: pver        
      !integer, parameter :: pver = 30        
      integer, intent(in) :: pverp        
      integer, intent(in) :: pcols        
      integer, intent(in) :: zm_org_int        
      integer, intent(in) :: limcnv       

      real(8), intent(inout) :: t(pver)  
      real(8), intent(inout) :: qh(pver) 
      real(8), intent(inout) :: pap(pver)     
      real(8), intent(inout) :: paph(pver+1)
      real(8), intent(inout) :: dpp(pver)  
      real(8), intent(inout) :: zm(pver)
      real(8), intent(in) :: geos
      real(8), intent(inout) :: zi(pver+1)
      real(8), intent(in) :: pblh
      real(8), intent(in) :: tpert
      real(8), intent(in) :: landfrac 

      real(8), intent(out) :: qtnd(pver) 
      real(8), intent(out) :: heat(pver)           
      real(8), intent(out) :: mcon(pverp)
      real(8), intent(out) :: dlf(pver)   
      real(8), intent(out) :: pflx(pverp)
      real(8), intent(out) :: cme(pver)
      real(8), intent(out) :: zdu(pver)
      real(8), intent(out) :: rprd(pver) 
      real(8), intent(out) :: mu(pver)
      real(8), intent(out) :: eu(pver)
      real(8), intent(out) :: du(pver)
      real(8), intent(out) :: md(pver)
      real(8), intent(out) :: ed(pver)
      real(8), intent(out) :: dp(pver)       
      real(8), intent(out) :: cape      
      real(8), intent(out) :: dsubcld
      real(8), intent(out) :: jctop  
      real(8), intent(out) :: jcbot  
      real(8), intent(out) :: prec
      real(8), intent(out) :: rliq 
      integer, intent(inout) :: maxg     
      real(8), intent(inout) ::  ql(pver)
      real(8), intent(inout) :: org(pver)
      real(8), intent(inout) :: orgt(pver)
      real(8), intent(inout) :: org2d(pver)


      real(8) zs
      real(8) dlg(pver)  
      real(8) pflxg(pverp) 
      real(8) cug(pver)   
      real(8) evpg(pver) 
      real(8) orgavg
      real(8) dptot
      real(8) mumax
      integer jt        
      integer ideep   
      integer lengath
      real(8) pblt  


      real(8) q(pver)
      real(8) p(pver)
      real(8) z(pver)
      real(8) s(pver)
      real(8) tp(pver)
      real(8) zf(pver+1)
      real(8) pf(pver+1)
      real(8) qstp(pver)

      real(8) tl      
      integer lcl     
      integer lel    
      integer lon 
      integer maxi
      integer index

      real(8) precip
      !real(8) qg(pver)  
      !real(8) tg(pver) 
      !real(8) pg(pver)
      !real(8) zg(pver)
      !real(8) sg(pver)
      !real(8) tpg(pver)
      !real(8) zfg(pver+1) 
      !real(8) qstpg(pver)
      !real(8) ug(pver)  
      !real(8) vg(pver) 
      !real(8) cmeg(pver)

      real(8) rprdg(pver) 
      real(8) capeg     
      real(8) tlg      
      real(8) landfracg 

      integer lclg      
      integer lelg
      !real(8) dqdt(pver)    
      !real(8) dsdt(pver)    
      !real(8) sd(pver)      
      !real(8) qd(pver)      
      !real(8) mc(pver)      
      !real(8) qhat(pver)    
      !real(8) qu(pver)      
      !real(8) su(pver)      
      !real(8) qs(pver)      
      !real(8) shat(pver)    
      !real(8) hmn(pver)     
      !real(8) hsat(pver)    
      !real(8) qlg(pver)
      !real(8) dudt(pver)    
      !real(8) dvdt(pver)    

      real(8) mb           

      integer jlcl
      integer j0           
      integer jd           


      real(8), intent(in) :: delt               
      real(8), intent(in) :: grav              
      real(8), intent(in) :: rgrav               
      real(8), intent(in) :: cpres               
      real(8), intent(in) :: rl               
      real(8), intent(in) :: rgas               

      integer i
      integer ii
      integer k
      integer msg                
      real(8) qdifr
      real(8) sdifr
      logical :: zm_org 
      real(8), intent(in) :: cpliq
      real(8), intent(in) :: cpwv 
      real(8), intent(in) :: tfreez 
      real(8), intent(in) :: eps1
      real(8), intent(in) :: rh2o
      real(8), intent(in) :: latice 
      real(8), intent(in) :: c0_ocn 
      real(8), intent(in) :: c0_lnd
!      integer, intent(in) :: no_deep_pbl
      integer, intent(in) :: no_deep_pbl_in

      logical :: no_deep_pbl
      real(8) :: abst
      integer(8) :: log_addr, log10_addr, pow_addr
      real(8), parameter :: capelmt = 70.0d0
      integer :: lelk
      integer :: kcount_sig 
      integer :: kcount
      integer, intent(out) :: tttt(4)

      call math_agent_log_c(log_addr)
      call math_agent_log10_c(log10_addr)
      call math_agent_pow_c(pow_addr)
      !integer, external :: lwpf_D_start
      if(no_deep_pbl_in .eq. 1) then
          no_deep_pbl = .true.
      else
          no_deep_pbl = .false.
      endif
!--
      if(zm_org_int .eq. 1) then
          zm_org = .true.
      else
          zm_org = .false.
      endif
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
      msg = limcnv - 1
!
! initialize necessary arrays.
! zero out variables not used in cam


     if (zm_org) then
      orgt(:) = 0d0
     end if

     qtnd(:) = 0.0d0
     heat(:) = 0.0d0
     mcon(:) = 0.0d0
     rliq   = 0.0d0
!
! initialize convective tendencies
!
      prec = 0d0
      do k = 1,pver
      !dqdt(k)  = 0.0d0
      !dsdt(k)  = 0.0d0
      !dudt(k)  = 0.0d0
      !dvdt(k)  = 0.0d0
      pflx(k)  = 0.0d0
      pflxg(k) = 0.0d0
      cme(k)   = 0.0d0
      rprd(k)  = 0.0d0
      zdu(k)   = 0.0d0
      ql(k)    = 0.0d0
      !qlg(k)   = 0.0d0
      dlf(k)   = 0.0d0
      dlg(k)   = 0.0d0
      end do
      pflx(pverp) = 0.0d0
      pflxg(pverp) = 0.0d0
      pblt = pver
      dsubcld = 0.0d0

      jctop = pver
      jcbot = 1


!       compute vertical average here
      if (zm_org) then
          orgavg = 0.0d0
          dptot = 0.0d0

          do k = 1, pver
          if (org(k) .gt. 0) then
              orgavg = orgavg+dpp(k)*org(k)
              dptot = dptot+dpp(k)
          endif
          enddo  

     if (dptot .gt. 0) then
          orgavg = orgavg/dptot
     endif

          do k = 1, pver
          org2d(k) = orgavg
          enddo

     endif

!
! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
!
          zs = geos*rgrav
          pf(pver+1) = paph(pver+1)*0.01d0
          zf(pver+1) = zi(pver+1) + zs
          do k = 1,pver
          p(k) = pap(k)*0.01d0
          pf(k) = paph(k)*0.01d0
          z(k) = zm(k) + zs
          zf(k) = zi(k) + zs
          end do
          !print *, loc(zf(:))

          do k = pver - 1, msg + 1,-1
          if (abs(z(k)-zs-pblh) < (zf(k)-zf(k+1))*0.5d0) pblt = k
          end do

          !call lwpf_d_start()
          !call lwpf_d_stop()

          do k = 1,pver
          q(k) = qh(k)
          s(k) = t(k) + (grav/cpres)*z(k)
          tp(k)=0.0d0
          !shat(k) = s(k)
          !qhat(k) = q(k)
          end do

          capeg = 0.0d0
          lclg = 1
          lelg = pver
          maxg = 1
          tlg = 400.0d0
          dsubcld = 0.0d0

!call lwpf_e_start()
          call buoyan_dilute_compute(lchnk   ,ncol    , &
          q       ,t       ,p       ,z       ,pf       , &
          tp      ,qstp    ,tl      ,rl      ,cape     , &
          pblt    ,lcl     ,lel     ,lon     ,maxi     , &
          rgas    ,grav    ,cpres   ,msg     , &
          tpert   , org2d  , landfrac, pver, zm_org, cpliq, &
          cpwv, tfreez, eps1, rh2o, latice, log_addr,   &
          log10_addr, pow_addr &
          )

          tttt(1) = lcl
          tttt(2) = lel
          tttt(3) = lon
          tttt(4) = maxi 

          do k = 1, pver
          qtnd(k) = tp(k)
          dlf(k) = qstp(k) 
          enddo
          heat(1) = tl
          cme(1) = cape
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          return 
!------------------------------------------------------------------
#ifdef CPEDEBUG
!call lwpf_f_start()

          lengath = 0
          if (cape > capelmt) then
              lengath = lengath + 1
              index = i
          end if

          kcount = 0
          call dma_sum(lengath, kcount)  

           lengath = kcount

if (lengath.eq.1) ideep =  1 

if (lengath.eq.0) then
               lelk = pverp 
               call dma_min(lelk)
               lelk = 1
               call dma_max(lelk)
               
               lengath = kcount
               kcount = 0
               kcount_sig = 0
               do k = pver, msg + 2, -1
                    call dma_sum(kcount_sig, kcount)
                    if(kcount >= lengath) then
                        goto 960
                    endif
               enddo
             960 continue 
else
     ideep = 1
     lengath = kcount
          
       !end do
       do k = 1,pver
!       do i = 1,lengath
         dp(k) = 0.01d0*dpp(k)
         qg(k) = q(k)
         tg(k) = t(k)
         pg(k) = p(k)
         zg(k) = z(k)
         sg(k) = s(k)
         tpg(k) = tp(k)
         zfg(k) = zf(k)
         qstpg(k) = qstp(k)
         ug(k) = 0.0d0
         vg(k) = 0.0d0
      end do
!        end do
!
!   do i = 1,lengath
      zfg(pver+1) = zf(pver+1)
!   end do
!    do i = 1,lengath
      capeg = cape
      lclg = lcl
      lelg = lel
      maxg = maxi
      tlg = tl
      landfracg = landfrac
!   end do
!
!
       do k = msg + 1,pver
!      do i = 1,lengath
         if (k >= maxg) then
            dsubcld = dsubcld + dp(k)
         end if
!      end do
       end do
!
       do k = msg + 2,pver
!      do i = 1,lengath
         sdifr = 0.0d0
         qdifr = 0.0d0
         if (sg(k) > 0.0d0 .or. sg(k-1) > 0.0d0) &
            sdifr = abs((sg(k)-sg(k-1))/max(sg(k-1),sg(k)))
         if (qg(k) > 0.0d0 .or. qg(k-1) > 0.0d0) &
            qdifr = abs((qg(k)-qg(k-1))/max(qg(k-1),qg(k)))
         if (sdifr > 1.0d-6) then
            shat(k) = log(sg(k-1)/sg(k))*sg(k-1)*sg(k)/(sg(k-1)-sg(k))
         else
            shat(k) = 0.5d0* (sg(k)+sg(k-1))
         end if
         if (qdifr > 1.0d-6) then
            qhat(k) = log(qg(k-1)/qg(k))*qg(k-1)*qg(k)/(qg(k-1)-qg(k))
         else
            qhat(k) = 0.5d0* (qg(k)+qg(k-1))
         end if
!      end do
       end do

      call cldprp(lchnk   , &
               qg      ,tg      ,ug      ,vg      ,pg      , &
               zg      ,sg      ,mu      ,eu      ,du      , &
               md      ,ed      ,sd      ,qd      ,mc      , &
               qu      ,su      ,zfg     ,qs      ,hmn     , &
               hsat    ,shat    ,qlg     , &
               cmeg    ,maxg    ,lelg    ,jt      ,jlcl    , &
               maxg    ,j0      ,jd      ,rl      ,lengath , &
               rgas    ,grav    ,cpres   ,msg     , &
               pflxg   ,evpg    ,cug     ,rprdg   ,limcnv  ,landfracg,&
               pver, pverp, tfreez, eps1, rh2o, c0_ocn, c0_lnd, &
               log_addr, log10_addr, pow_addr &
               )

               !qtnd, dlf, heat, mcon, pflx, cme, &
               !zdu, rprd)
!   real(r8), intent(out) :: rprd(pcols,pver)     ! rate of production of precip at that layer
!   real(r8), intent(out) :: pflx(pcols,pverp)    ! precipitation flux thru layer
!
!   real(r8), intent(out) :: du(pcols,pver)       ! detrainement rate of updraft
!   real(r8), intent(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
!   real(r8), intent(out) :: eu(pcols,pver)       ! entrainment rate of updraft
!   real(r8), intent(out) :: md(pcols,pver)       ! downdraft mass flux
!   real(r8), intent(out) :: ql(pcols,pver)       ! liq water of updraft
!   real(r8), intent(out) :: mu(pcols,pver)       ! updraft mass flux
       
!   real(r8), intent(out) :: hmn(pcols,pver)      ! moist stat energy of env
!   real(r8), intent(out) :: hsat(pcols,pver)     ! sat moist stat energy of env
!   real(r8), intent(out) :: mc(pcols,pver)       ! net mass flux
!   real(r8), intent(out) :: qd(pcols,pver)       ! spec humidity of downdraft
!   real(r8), intent(out) :: qst(pcols,pver)      ! saturation mixing ratio of env.
!   real(r8), intent(out) :: qu(pcols,pver)       ! spec hum of updraft
!   real(r8), intent(out) :: sd(pcols,pver)       ! normalized dry stat energy of downdraft
!   real(r8), intent(out) :: su(pcols,pver)       ! normalized dry stat energy of updraft
                  do k = 1, pver

                    mcon(k) = hmn(k)
                    zdu(k) = hsat(k)
                    dp(k) = mc(k)
                    t(k) = qd(k)
                    qh(k) = qs(k)
                    pap(k) = qu(k)
                    dpp(k) = sd(k)
                    zm(k) = su(k)

                    !rprd(k) = sd(k)
                  enddo
       !goto 886

       do k = msg + 1,pver
!      do i = 1,lengath
         du   (k) = du   (k)* (zfg(k)-zfg(k+1))/dp(k)
         eu   (k) = eu   (k)* (zfg(k)-zfg(k+1))/dp(k)
         ed   (k) = ed   (k)* (zfg(k)-zfg(k+1))/dp(k)
         cug  (k) = cug  (k)* (zfg(k)-zfg(k+1))/dp(k)
         cmeg (k) = cmeg (k)* (zfg(k)-zfg(k+1))/dp(k)
         rprdg(k) = rprdg(k)* (zfg(k)-zfg(k+1))/dp(k)
         evpg (k) = evpg (k)* (zfg(k)-zfg(k+1))/dp(k)
!      end do
       end do

       call closure(lchnk   , &
                qg      ,tg      ,pg      ,zg      ,sg      , &
                tpg     ,qs      ,qu      ,su      ,mc      , &
                du      ,mu      ,md      ,qd      ,sd      , &
                qhat    ,shat    ,dp      ,qstpg   ,zfg     , &
                qlg     ,dsubcld ,mb      ,capeg   ,tlg     , &
                lclg    ,lelg    ,jt      ,maxg    ,1       , &
                lengath ,rgas    ,grav    ,cpres   ,rl      , &
                msg     ,pver    ,eps1  ,&
                log_addr, log10_addr, pow_addr)

                  !do k = 1, pver

                  !  mcon(k) = hmn(k)
                  !  zdu(k) = hsat(k)
                  !  dp(k) = mc(k)
                  !  t(k) = qd(k)
                  !  qh(k) = qs(k)
                  !  pap(k) = qu(k)
                  !  dpp(k) = sd(k)
                  !  zm(k) = su(k)

                  !  !rprd(k) = sd(k)
                  !enddo


!   do i=1,lengath
      mumax = 0
!   end do
     do k=msg + 2,pver
!      do i=1,lengath
        mumax = max(mumax, mu(k)/dp(k))
!      end do
     end do

!   do i=1,lengath
      if (mumax > 0.d0) then
         mb = min(mb,0.5d0/(delt*mumax))
      else
         mb = 0.d0
      endif
!   end do
! If no_deep_pbl = .true., don't allow convection entirely 
! within PBL (suggestion of Bjorn Stevens, 8-2000)

     if (no_deep_pbl) then
!      do i=1,lengath
         if (zm(jt) < pblh) mb = 0
!      end do
     end if

     do k=msg+1,pver
!      do i=1,lengath
         mu   (k)  = mu   (k)*mb
         md   (k)  = md   (k)*mb
         mc   (k)  = mc   (k)*mb
         du   (k)  = du   (k)*mb
         eu   (k)  = eu   (k)*mb
         ed   (k)  = ed   (k)*mb
         cmeg (k)  = cmeg (k)*mb
         rprdg(k)  = rprdg(k)*mb
         cug  (k)  = cug  (k)*mb
         evpg (k)  = evpg (k)*mb
         pflxg(k+1)= pflxg(k+1)*mb*100.d0/grav
!      end do
     end do
!
! compute temperature and moisture changes due to convection.
!
        call q1q2_pjr(lchnk   , &
                 dqdt    ,dsdt    ,qg      ,qs      ,qu      , &
                 su      ,du      ,qhat    ,shat    ,dp      , &
                 mu      ,md      ,sd      ,qd      ,qlg     , &
                 dsubcld ,jt      ,maxg    ,1       ,lengath , &
                 cpres   ,rl      ,msg     ,          &
                 dlg     ,evpg    ,cug     ,pver,  &
                 log_addr, log10_addr, &
                 pow_addr)

!
!
!
!!!!!xukai need to fix
       do k = msg + 1,pver
!      do i = 1,lengath
            q(k) = qh(k) + 2.d0*delt*dqdt(k)
         qtnd(k) = dqdt (k)
         cme (k) = cmeg (k)
         rprd(k) = rprdg(k)
         zdu (k) = du   (k)
         mcon(k) = mc   (k)
         heat(k) = dsdt (k)*cpres
         dlf (k) = dlg  (k)
         pflx(k) = pflxg(k)
         ql  (k) = qlg  (k)

!      end do
     end do
!!
!!!!!!DIR$ CONCURRENT
!   do i = 1,lengath
      jctop = jt
!++bee
      jcbot = maxg
!--bee
      pflx(pverp) = pflxg(pverp)
!   end do


end if
  886 continue
      return
      if(lengath > 0) then

!!!!!!! loop change
     do k = pver,msg + 1,-1
!      do i = 1,ncol
         prec = prec - dpp(k)* (q(k)-qh(k)) - dpp(k)*dlf(k)*2*delt
!      end do
     end do

!   do i = 1,ncol
      prec = rgrav*max(prec,0.d0)/ (2.d0*delt)/1000.d0
!   end do

     do k = 1, pver
!      do i = 1, ncol
         rliq = rliq + dlf(k)*dpp(k)/grav
         !it
!      end do
     end do
     rliq = rliq /1000.d0

    endif
#endif

!call lwpf_f_stop()


          end subroutine

