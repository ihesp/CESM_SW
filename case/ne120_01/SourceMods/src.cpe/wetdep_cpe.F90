!/*************************************************************************
!	> File Name: wetdep_cpe.F90
!	> Author: Xu Kai
!	> Created Time: 2019年01月14日 星期一 13时38分30秒
! ************************************************************************/
	
subroutine wetdepa_v2_compute(                          &
   pcols, pver, gravit,                                 &
   p, q, pdel, cldt, cldc,                              &
   cmfdqr, evapc, conicw, precs, conds,                 &
   evaps, cwat, tracer, deltat, scavt,                  &
   iscavt, cldvcu, cldvst, dlf, fracis,                 &
   sol_fact, ncol, scavcoef, is_strat_cloudborne, qqcw, &
   f_act_conv, icscavt, isscavt, bcscavt, bsscavt,      &
   sol_facti_in, sol_factic_in )

   implicit none
   !----------------------------------------------------------------------- 
   !
   ! scavenging code for very soluble aerosols
   ! 
   !-----------------------------------------------------------------------

   integer, intent(in) :: pcols, pver
   real(8), intent(in) :: gravit
   real(8), intent(in) ::&
      p(pver),        &! pressure
      q(pver),        &! moisture
      pdel(pver),     &! pressure thikness
      cldt(pver),     &! total cloud fraction
      cldc(pver),     &! convective cloud fraction
      cmfdqr(pver),   &! rate of production of convective precip
      evapc(pver),    &! Evaporation rate of convective precipitation
      conicw(pver),   &! convective cloud water
      cwat(pver),     &! cloud water amount 
      precs(pver),    &! rate of production of stratiform precip
      conds(pver),    &! rate of production of condensate
      evaps(pver),    &! rate of evaporation of precip
      cldvcu(pver),   &! Convective precipitation area at the top interface of each layer
      cldvst(pver),   &! Stratiform precipitation area at the top interface of each layer
      dlf(pver),      &! Detrainment of convective condensate [kg/kg/s]
      deltat,               &! time step
      tracer(pver)     ! trace species

   ! If subroutine is called with just sol_fact:
   !    sol_fact is used for both in- and below-cloud scavenging
   ! If subroutine is called with optional argument sol_facti_in:
   !    sol_fact  is used for below cloud scavenging
   !    sol_facti is used for in cloud scavenging

   real(8), intent(in)  :: sol_fact
   integer,  intent(in)  :: ncol
   real(8), intent(in)  :: scavcoef(pver) ! Dana and Hales coefficient (/mm) (0.1 if not MODAL_AERO)
   real(8), intent(out) ::&
      scavt(pver),   &! scavenging tend 
      iscavt(pver),  &! incloud scavenging tends
      fracis(pver)    ! fraction of species not scavenged

   ! Setting is_strat_cloudborne=.true. indicates that tracer is stratiform-cloudborne aerosol.
   !   This is only used by MAM code.  The optional args qqcw and f_act_conv are not referenced
   !   in this case.
   ! Setting is_strat_cloudborne=.false. is being used to indicate that the tracers are the
   !   interstitial modal aerosols.  In this case the optional qqcw (the cloud borne mixing ratio
   !   corresponding to the interstitial aerosol) must be provided, as well as the optional f_act_conv.
   logical,  intent(in), optional :: is_strat_cloudborne   
   real(8), intent(in), optional :: qqcw(pver)
   real(8), intent(in), optional :: f_act_conv(pver)

   real(8), intent(in), optional :: sol_facti_in   ! solubility factor (frac of aerosol scavenged in cloud)
   real(8), intent(in), optional :: sol_factic_in(pver)  ! sol_facti_in for convective clouds
         

   real(8), intent(out), optional :: icscavt(pver)     ! incloud, convective
   real(8), intent(out), optional :: isscavt(pver)     ! incloud, stratiform
   real(8), intent(out), optional :: bcscavt(pver)     ! below cloud, convective
   real(8), intent(out), optional :: bsscavt(pver)     ! below cloud, stratiform

   ! local variables

   integer :: i, k

   real(8) :: omsm                 ! 1 - (a small number)
   real(8) :: clds          ! stratiform cloud fraction
   real(8) :: fracev        ! fraction of precip from above that is evaporating
   real(8) :: fracev_cu     ! Fraction of convective precip from above that is evaporating
   real(8) :: fracp         ! fraction of cloud water converted to precip
   real(8) :: pdog          ! work variable (pdel/gravit)
   real(8) :: rpdog         ! work variable (gravit/pdel)
   real(8) :: precabc       ! conv precip from above (work array)
   real(8) :: precabs       ! strat precip from above (work array)
   real(8) :: rat           ! ratio of amount available to amount removed
   real(8) :: scavab        ! scavenged tracer flux from above (work array)
   real(8) :: scavabc       ! scavenged tracer flux from above (work array)
   real(8) :: srcc          ! tend for convective rain
   real(8) :: srcs          ! tend for stratiform rain
   real(8) :: srct          ! work variable

   real(8) :: fins          ! fraction of rem. rate by strat rain
   real(8) :: finc          ! fraction of rem. rate by conv. rain
   real(8) :: conv_scav_ic  ! convective scavenging incloud
   real(8) :: conv_scav_bc  ! convective scavenging below cloud
   real(8) :: st_scav_ic    ! stratiform scavenging incloud
   real(8) :: st_scav_bc    ! stratiform scavenging below cloud

   real(8) :: odds          ! limit on removal rate (proportional to prec)
   real(8) :: dblchek
   logical :: found

   real(8) :: trac_qqcw
   real(8) :: tracer_incu
   real(8) :: tracer_mean

   ! For stratiform cloud, cloudborne aerosol is treated explicitly,
   !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
   ! For convective cloud, cloudborne aerosol is not treated explicitly,
   !    and sol_factic is 1.0 for both cloudborne and interstitial.

   real(8) :: sol_facti              ! in cloud fraction of aerosol scavenged
   real(8) :: sol_factb              ! below cloud fraction of aerosol scavenged
   real(8) :: sol_factic(pver) ! in cloud fraction of aerosol scavenged for convective clouds

   real(8) :: rdeltat
   ! ------------------------------------------------------------------------

   omsm = 1.d0-2*epsilon(1.d0) ! used to prevent roundoff errors below zero

   ! default (if other sol_facts aren't in call, set all to required sol_fact)
   sol_facti = sol_fact
   sol_factb = sol_fact

   if ( present(sol_facti_in) )  sol_facti = sol_facti_in

   sol_factic  = sol_facti
   if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in

   ! this section of code is for highly soluble aerosols,
   ! the assumption is that within the cloud that
   ! all the tracer is in the cloud water
   !
   ! for both convective and stratiform clouds, 
   ! the fraction of cloud water converted to precip defines
   ! the amount of tracer which is pulled out.

   precabs = 0.0d0
   precabc = 0.0d0
   scavab  = 0.0d0
   scavabc = 0.0d0

   do k = 1, pver

         clds  = cldt(k) - cldc(k)
         pdog  = pdel(k)/gravit
         rpdog = gravit/pdel(k)
         rdeltat  = 1.0d0/deltat

         ! ****************** Evaporation **************************
         ! calculate the fraction of strat precip from above 
         !                 which evaporates within this layer
         fracev = evaps(k)*pdog &
                     /max(1.d-12,precabs)

         ! trap to ensure reasonable ratio bounds
         fracev = max(0.d0,min(1.d0,fracev))

         ! Same as above but convective precipitation part
         fracev_cu = evapc(k)*pdog/max(1.d-12,precabc)
         fracev_cu = max(0.d0,min(1.d0,fracev_cu))

         ! ****************** Convection ***************************
         !
         ! set odds proportional to fraction of the grid box that is swept by the 
         ! precipitation =precabc/rhoh20*(area of sphere projected on plane
         !                                /volume of sphere)*deltat
         ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
         ! unless the fraction of the area that is cloud is less than odds, in which
         ! case use the cloud fraction (assumes precabs is in kg/m2/s)
         ! is really: precabs*3/4/1000./1e-3*deltat
         ! here I use .1 from Balkanski
         !
         ! use a local rate of convective rain production for incloud scav
         !
         ! Fraction of convective cloud water converted to rain.  This version is used
         ! in 2 of the 3 branches below before fracp is reused in the stratiform calc.
         ! NB: In below formula for fracp conicw is a LWC/IWC that has already
         !     precipitated out, i.e., conicw does not contain precipitation

         fracp = cmfdqr(k)*deltat / &
                    max( 1.d-12, cldc(k)*conicw(k) + (cmfdqr(k)+dlf(k))*deltat )
         fracp = max( min( 1.d0, fracp), 0.d0 )

         if ( present(is_strat_cloudborne) ) then

            if ( is_strat_cloudborne ) then

               ! convective scavenging

               conv_scav_ic = 0.d0

               conv_scav_bc = 0.d0

               ! stratiform scavenging

               fracp = precs(k)*deltat / &
                          max( 1.d-12, cwat(k) + precs(k)*deltat )
               fracp = max( 0.d0, min(1.d0, fracp) )
               st_scav_ic = sol_facti *fracp*tracer(k)*rdeltat

               st_scav_bc = 0.d0

            else

               ! convective scavenging

               trac_qqcw = min(qqcw(k), &
                                  tracer(k)*( clds/max( 0.01d0, 1.d0-clds ) ) )

               tracer_incu = f_act_conv(k)*(tracer(k) + trac_qqcw)

               conv_scav_ic = sol_factic(k)*cldc(k)*fracp*tracer_incu*rdeltat

               tracer_mean = tracer(k)*(1.d0 - cldc(k)*f_act_conv(k)) - &
                                cldc(k)*f_act_conv(k)*trac_qqcw
               tracer_mean = max(0.d0,tracer_mean)

               odds = precabc/max(cldvcu(k),1.d-5)*scavcoef(k)*deltat
               odds = max(min(1.d0,odds),0.d0)
               conv_scav_bc = sol_factb *cldvcu(k)*odds*tracer_mean*rdeltat


               ! stratiform scavenging

               st_scav_ic = 0.d0

               odds = precabs/max(cldvst(k),1.d-5)*scavcoef(k)*deltat
               odds = max(min(1.d0,odds),0.d0)
               st_scav_bc = sol_factb *cldvst(k)*odds*tracer_mean*rdeltat

            end if

         else

            ! convective scavenging

            conv_scav_ic = sol_factic(k)*cldc(k)*fracp*tracer(k)*rdeltat

            odds = precabc/max(cldvcu(k), 1.d-5)*scavcoef(k)*deltat
            odds = max( min(1.d0, odds), 0.d0)
            conv_scav_bc = sol_factb*cldvcu(k)*odds*tracer(k)*rdeltat

            ! stratiform scavenging

            ! fracp is the fraction of cloud water converted to precip
            ! NB: In below formula for fracp cwat is a LWC/IWC that has already
            !     precipitated out, i.e., cwat does not contain precipitation
            fracp = precs(k)*deltat / &
                       max( 1.d-12, cwat(k) + precs(k)*deltat )
            fracp = max( 0.d0, min( 1.d0, fracp ) )
            
            ! assume the corresponding amnt of tracer is removed
            st_scav_ic = sol_facti*clds*fracp*tracer(k)*rdeltat

            odds = precabs/max(cldvst(k),1.d-5)*scavcoef(k)*deltat
            odds = max(min(1.d0,odds),0.d0)
            st_scav_bc =sol_factb*(cldvst(k)*odds) *tracer(k)*rdeltat

         end if

         ! total convective scavenging
         srcc = conv_scav_ic + conv_scav_bc
         finc = conv_scav_ic/(srcc + 1.d-36)

         ! total stratiform scavenging
         srcs = st_scav_ic + st_scav_bc
         fins = st_scav_ic/(srcs + 1.d-36)

         ! make sure we dont take out more than is there
         ! ratio of amount available to amount removed
         rat = tracer(k)/max(deltat*(srcc+srcs),1.d-36)
         if (rat.lt.1.d0) then
            srcs = srcs*rat
            srcc = srcc*rat
         endif
         srct = (srcc+srcs)*omsm

            
         ! fraction that is not removed within the cloud
         ! (assumed to be interstitial, and subject to convective transport)
         fracp = deltat*srct/max(cldvst(k)*tracer(k),1.d-36)  ! amount removed
         fracp = max(0.d0,min(1.d0,fracp))
         fracis(k) = 1.d0 - fracp

         ! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
         ! Sungsu added cumulus contribution in the below 3 blocks
         scavt(k) = -srct + (fracev*scavab+fracev_cu*scavabc)*rpdog
         iscavt(k) = -(srcc*finc + srcs*fins)*omsm

         if ( present(icscavt) ) icscavt(k) = -(srcc*finc) * omsm
         if ( present(isscavt) ) isscavt(k) = -(srcs*fins) * omsm
         if ( present(bcscavt) ) bcscavt(k) = -(srcc * (1-finc)) * omsm +  &
            fracev_cu*scavabc*rpdog

         if ( present(bsscavt) ) bsscavt(k) = -(srcs * (1-fins)) * omsm +  &
            fracev*scavab*rpdog

         !icscavt(k) = -(srcc*finc) * omsm
         !isscavt(k) = -(srcs*fins) * omsm
         !bcscavt(k) = -(srcc * (1-finc)) * omsm +  &
         !   fracev_cu*scavabc*rpdog

         !bsscavt(k) = -(srcs * (1-fins)) * omsm +  &
         !   fracev*scavab*rpdog
         dblchek = tracer(k) + deltat*scavt(k)

         ! now keep track of scavenged mass and precip
         scavab = scavab*(1-fracev) + srcs*pdog
         precabs = precabs + (precs(k) - evaps(k))*pdog
         scavabc = scavabc*(1-fracev_cu) + srcc*pdog
         precabc = precabc + (cmfdqr(k) - evapc(k))*pdog

      !end do ! End of i = 1, ncol

      !found = .false.
      !do i = 1,ncol
      !   if ( dblchek < 0.d0 ) then
      !      found = .true.
      !      exit
      !   end if
      !end do

      !if ( found ) then
      !   do i = 1,ncol
      !      if (dblchek .lt. 0._r8) then
      !         write(iulog,*) ' wetdapa: negative value ', i, k, tracer(k), &
      !            dblchek, scavt(k), srct, rat, fracev
      !      endif
      !   end do
      !endif

   end do ! End of k = 1, pver

end subroutine wetdepa_v2_compute
	
