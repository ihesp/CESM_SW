
  module cldwat2m_macro_cpe
#define PARASIZE 12
  !--------------------------------------------------- !
  ! Purpose     : CAM Interface for Cloud Macrophysics !
  ! Author      : Sungsu Park                          !
  ! Description : Park et al. 2010.                    !
  ! For questions, contact Sungsu Park                 !
  !                        e-mail : sungsup@ucar.edu   !
  !                        phone  : 303-497-1375       !
  !--------------------------------------------------- !

   !use shr_kind_mod,     only: r8=>shr_kind_r8
   !use spmd_utils,       only: masterproc
   !use ppgrid,           only: pcols, pver, pverp
   !use cam_abortutils,   only: endrun
   !use physconst,        only: cpair, latvap, latice, rh2o, gravit, rair
   use wv_saturation_cpe_xk,    only: qsat_water, svp_water, svp_ice, qsat_ice
   !use cam_history,      only: addfld, phys_decomp, outfld, hist_fld_active
   !use cam_logfile,      only: iulog
   !use ref_pres,         only: top_lev=>trop_cloud_top_lev
   use cldfrc2m_cpe,         only: astG_PDF_single, astG_PDF, astG_RHU_single, &
                               astG_RHU, aist_single, aist_vector

    !use wv_sat_methods_cpe, only: tmp1, tmp2, tmpind, wvflag
   implicit none
   private
   save

   public ::           &
        rhcrit_calc,   &
        instratus_condensate, & 
        instratus_core,&
        funcd_instratus, &
        gridmean_RH, &
        positive_moisture, & 
        gaussj
        interface 
            real(8) function pow_c(in1, in2)
            real(8), intent(in) :: in1
            real(8), intent(in) :: in2
            end function
        end interface


  ! -------------- !
   ! Set Parameters !
   ! -------------- !

!   ! ------------------------------------------------------------------------------- !
!   ! Parameter used for selecting generalized critical RH for liquid and ice stratus !
!   ! ------------------------------------------------------------------------------- !
!
!   integer :: i_rhminl ! This is for liquid stratus fraction.
!                       ! If 0 : Original fixed critical RH from the namelist.
!                       ! If 1 : Add convective detrainment effect on the above '0' option. 
!                       !        In this case, 'tau_detw' [s] should be specified below.
!                       ! If 2 : Use fully scale-adaptive method.
!                       !        In this case, 'tau_detw' [s] and 'c_aniso' [no unit] should
!                       !        be specified below. 
!
!   integer :: i_rhmini ! This is for ice stratus fraction.
!                       ! If 0 : Original fixed critical RH from the namelist.
!                       ! If 1 : Add convective detrainment effect on the above '0' option. 
!                       !        In this case, 'tau_deti' [s] should be specified below.
!                       ! If 2 : Use fully scale-adaptive method.
!                       !        In this case, 'tau_deti' [s] and 'c_aniso' [no unit] should
!                       !        be specified below. 
!                       ! Note that 'micro_mg_cam' is using below 'rhmini_const', regardless
!                       ! of 'i_rhmini'.  This connection should be built in future.
!
!   real(8), parameter :: tau_detw =100.d0   ! Dissipation time scale of convective liquid condensate detrained
!                                              !  into the clear portion. [hr]. 0.5-3 hr is possible.
!   real(8), parameter :: tau_deti =  1.d0   ! Dissipation time scale of convective ice    condensate detrained
!                                              !  into the clear portion. [hr]. 0.5-3 hr is possible.
!   real(8), parameter :: c_aniso  =  1.d0   ! Inverse of anisotropic factor of PBL turbulence
!
!   ! ----------------------------- !
!   ! Parameters for Liquid Stratus !
!   ! ----------------------------- !
!
!   logical,  parameter  :: CAMstfrac    = .false.    ! If .true. (.false.),
!                                                     ! use Slingo (triangular PDF-based) liquid stratus fraction
!   real(8), parameter  :: qlst_min     = 2.d-5   ! Minimum in-stratus LWC constraint [ kg/kg ]
!   real(8), parameter  :: qlst_max     = 3.d-3   ! Maximum in-stratus LWC constraint [ kg/kg ]
!   real(8), parameter  :: cc           = 0.1d0     ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )
!   integer,  parameter  :: niter        = 2          ! For iterative computation of QQ with 'ramda' below.
!   real(8), parameter  :: ramda        = 0.5d0     ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
!   real(8), private    :: rhminl_const              ! Critical RH for low-level  liquid stratus clouds
!   real(8), private    :: rhminl_adj_land_const     ! rhminl adjustment for snowfree land
!   real(8), private    :: rhminh_const              ! Critical RH for high-level liquid stratus clouds
!   real(8), private    :: premit                    ! Top    height for mid-level liquid stratus fraction
!   real(8), private    :: premib                    ! Bottom height for mid-level liquid stratus fraction
!
!   real(8), parameter :: qsmall = 1.d-18         ! Smallest mixing ratio considered in the macrophysics

   contains

   ! -------------- !
   ! Initialization !
   ! -------------- !


!=======================================================================================================
subroutine endrun()
    integer :: endcesm
    endcesm = 1
end subroutine

subroutine rhcrit_calc( &
   ncol, dp, T0, p, &
   clrw_old, clri_old, tke, qtl_flx, &
   qti_flx, cmfr_det, qlr_det, qir_det, &
   rhmini_arr, rhminl_arr, rhminl_adj_land_arr, rhminh_arr, &
   d_rhmin_liq_PBL, d_rhmin_ice_PBL, d_rhmin_liq_det, d_rhmin_ice_det,&
   wv_para, pver, pverp, top_lev, &
   i_rhmini, i_rhminl, rhmini_const, rhmaxi_const, rhminl_const, rhminl_adj_land_const, rhminh_const)
   ! ------------------------------------------------- !
   ! Compute a drop of critical RH for stratus by      !
   ! (1) PBL turbulence, and                           !
   ! (2) convective detrainment.                       !
   ! Note that all of 'd_rhmin...' terms are positive. !
   ! ------------------------------------------------- !
   implicit none

   integer, intent(in)    :: pver
   integer, intent(in)    :: pverp
   integer, intent(in)    :: top_lev 

   integer,  intent(in) :: ncol                         ! Number of active columns
   real(8), intent(in) :: dp(top_lev:pver)               ! Pressure thickness [Pa] > 0
   real(8), intent(in) :: T0(top_lev:pver)               ! Temperature [K]
   real(8), intent(in) :: p(top_lev:pver)                ! Pressure at the layer mid-point [Pa]
   real(8), intent(in) :: clrw_old(top_lev:pver)         ! Clear sky fraction at the previous time step for liquid stratus process
   real(8), intent(in) :: clri_old(top_lev:pver)         ! Clear sky fraction at the previous time step for    ice stratus process
   real(8), intent(in) :: tke(top_lev:pverp)                     ! (top_lev:pverp) TKE from the PBL scheme
   real(8), intent(in) :: qtl_flx(top_lev:pverp)                 ! (top_lev:pverp) overbar(w'qtl') from PBL scheme where qtl = qv + ql
   real(8), intent(in) :: qti_flx(top_lev:pverp)                 ! (top_lev:pverp) overbar(w'qti') from PBL scheme where qti = qv + qi
   real(8), intent(in) :: cmfr_det(top_lev:pver)                ! (top_lev:pver)  Detrained mass flux from the convection scheme
   real(8), intent(in) :: qlr_det(top_lev:pver)                 ! (top_lev:pver)  Detrained        ql from the convection scheme
   real(8), intent(in) :: qir_det(top_lev:pver)                 ! (top_lev:pver)  Detrained        qi from the convection scheme

   real(8), intent(out) :: rhmini_arr(top_lev:pver)
   real(8), intent(out) :: rhminl_arr(top_lev:pver)
   real(8), intent(out) :: rhminl_adj_land_arr(top_lev:pver)
   real(8), intent(out) :: rhminh_arr(top_lev:pver) 
   real(8), intent(out) :: d_rhmin_liq_PBL(top_lev:pver)
   real(8), intent(out) :: d_rhmin_ice_PBL(top_lev:pver)
   real(8), intent(out) :: d_rhmin_liq_det(top_lev:pver)
   real(8), intent(out) :: d_rhmin_ice_det(top_lev:pver)


   ! local variables

   integer ::  k

   real(8) :: esat_tmp          ! Dummy for saturation vapor pressure calc.
   real(8) :: qsat_tmp          ! Saturation water vapor specific humidity [kg/kg]
   real(8) :: sig_tmp

   real(8) :: gravit
!xukai++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, intent(in) :: i_rhmini 
   integer, intent(in) :: i_rhminl

   real(8), intent(in) :: wv_para(PARASIZE) 

   real(8) :: tau_deti !!=  1.d0   ! Dissipation time scale of convective ice    condensate detrained
   real(8) :: qsmall != 1.d-18         ! Smallest mixing ratio considered in the macrophysics
   real(8) :: c_aniso  !=  1.d0   ! Inverse of anisotropic factor of PBL turbulence

   real(8), intent(in)    :: rhmini_const     
   real(8), intent(in)    :: rhmaxi_const     
   real(8), intent(in)    :: rhminl_const     
   real(8), intent(in)    :: rhminl_adj_land_const 
   real(8), intent(in)    :: rhminh_const         

   real(8) :: tau_detw !=100.d0   ! Dissipation time scale of convective liquid condensate detrained
   real(8)  :: qlst_min !     = 2.d-5   ! Minimum in-stratus LWC constraint [ kg/kg ]
   real(8)  :: qlst_max !    = 3.d-3   ! Maximum in-stratus LWC constraint [ kg/kg ]
   !integer(8), intent(out) :: spa2

    tau_deti =  1.d0   ! Dissipation time scale of convective ice    condensate detrained
    qsmall = 1.d-18         ! Smallest mixing ratio considered in the macrophysics
    c_aniso  =  1.d0   ! Inverse of anisotropic factor of PBL turbulence
    tau_detw =100.d0   ! Dissipation time scale of convective liquid condensate detrained
    qlst_min     = 2.d-5   ! Minimum in-stratus LWC constraint [ kg/kg ]
    qlst_max     = 3.d-3   ! Maximum in-stratus LWC constraint [ kg/kg ]
 

!xukai--------------------------------------------------------------------------------------
   !---------------------------------------------------------------------------------------------------
   gravit = wv_para(11)



   ! ---------------------------------- !
   ! Calc critical RH for ice stratus   !
   ! ---------------------------------- !

   rhmini_arr(top_lev:pver) = rhmini_const

   if (i_rhmini > 0) then

      ! Compute the drop of critical RH by convective detrainment of cloud condensate

      do k = top_lev, pver
         !do i = 1, ncol
            d_rhmin_ice_det(k) = tau_deti*(gravit/dp(k))*cmfr_det(k)*clri_old(k)*qir_det(k)*3.6d6 
            d_rhmin_ice_det(k) = max(0.d0,min(0.5d0,d_rhmin_ice_det(k)))
         !end do
      end do

      if (i_rhmini == 1) then
         rhmini_arr(top_lev:pver) = rhmini_const - d_rhmin_ice_det(top_lev:pver)
      end if

   end if

   if (i_rhmini == 2) then

      ! Compute the drop of critical RH by the variability induced by PBL turbulence

      do k = top_lev, pver
         call qsat_ice(T0(k), p(k), esat_tmp, qsat_tmp, wv_para)

         !do i = 1, ncol
            sig_tmp = 0.5d0 * ( qti_flx(k)   / sqrt(max(qsmall,tke(k))) + & 
                                 qti_flx(k+1) / sqrt(max(qsmall,tke(k+1))) )
            d_rhmin_ice_PBL(k) = c_aniso*sig_tmp/max(qsmall,qsat_tmp) 
            d_rhmin_ice_PBL(k) = max(0.d0,min(0.5d0,d_rhmin_ice_PBL(k)))

            rhmini_arr(k) = 1.d0 - d_rhmin_ice_PBL(k) - d_rhmin_ice_det(k)
         !end do
      end do
   end if

   if (i_rhmini > 0) then
      do k = top_lev, pver
         !do i = 1, ncol
            rhmini_arr(k) = max(0.d0,min(rhmaxi_const,rhmini_arr(k))) 
         !end do
      end do
   end if

   ! ------------------------------------- !
   ! Choose critical RH for liquid stratus !
   ! ------------------------------------- !

   rhminl_arr(top_lev:pver)          = rhminl_const
   rhminl_adj_land_arr(top_lev:pver) = rhminl_adj_land_const
   rhminh_arr(top_lev:pver)          = rhminh_const

   if (i_rhminl > 0) then

      ! Compute the drop of critical RH by convective detrainment of cloud condensate

      do k = top_lev, pver
         !do i = 1, ncol
            d_rhmin_liq_det(k) = tau_detw*(gravit/dp(k))*cmfr_det(k)*clrw_old(k)*qlr_det(k)*3.6d6 
            d_rhmin_liq_det(k) = max(0.d0,min(0.5d0,d_rhmin_liq_det(k)))
         !end do
      end do

      if (i_rhminl == 1) then
         rhminl_arr(top_lev:pver) = rhminl_const - d_rhmin_liq_det(top_lev:pver)
         rhminh_arr(top_lev:pver) = rhminh_const - d_rhmin_liq_det(top_lev:pver)
      end if

   end if

   if (i_rhminl == 2) then

      ! Compute the drop of critical RH by the variability induced by PBL turbulence

      do k = top_lev, pver
         call qsat_water(T0(k), p(k), esat_tmp, qsat_tmp, wv_para)

         !do i = 1, ncol
            sig_tmp = 0.5d0 * ( qtl_flx(k)   / sqrt(max(qsmall,tke(k))) + & 
                                 qtl_flx(k+1) / sqrt(max(qsmall,tke(k+1))) )
            d_rhmin_liq_PBL(k) = c_aniso*sig_tmp/max(qsmall,qsat_tmp) 
            d_rhmin_liq_PBL(k) = max(0.d0,min(0.5d0,d_rhmin_liq_PBL(k)))

            rhminl_arr(k) = 1.d0 - d_rhmin_liq_PBL(k) - d_rhmin_liq_det(k)
            rhminl_adj_land_arr(k) = 0.d0
            rhminh_arr(k) = rhminl_arr(k)
         !end do
      end do
   end if

   if (i_rhminl > 0) then
      do k = top_lev, pver
         !do i = 1, ncol
            rhminl_arr(k) = max(rhminl_adj_land_arr(k),min(1.d0,rhminl_arr(k))) 
            rhminh_arr(k) = max(0.d0,min(1.d0,rhminh_arr(k))) 
         !end do
      end do
   end if
   !call get_sp(spa2)

end subroutine rhcrit_calc

!=======================================================================================================

   subroutine instratus_condensate( ncol, k,                      &  
                                    p_in, T0_in, qv0_in, ql0_in, qi0_in, & 
                                    ni0_in,                              &
                                    a_dc_in, ql_dc_in, qi_dc_in,         &
                                    a_sc_in, ql_sc_in, qi_sc_in,         & 
                                    landfrac, snowh,                     &
                                    rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in, &
                                    T_out, qv_out, ql_out, qi_out,       &
                                    al_st_out, ai_st_out, ql_st_out, qi_st_out,&
                                    wv_para, rhmaxi_const, premib, premit, &
                                    icecrit, iceopt)

   ! ------------------------------------------------------- !
   ! Diagnostically force in-stratus condensate to be        ! 
   ! in the range of 'qlst_min < qc_st < qlst_max'           !
   ! whenever stratus exists in the equilibrium state        !
   ! ------------------------------------------------------- !

   implicit none
   !integer,  intent(in)  :: lchnk                ! Chunk identifier
   integer,  intent(in)  :: ncol                 ! Number of atmospheric columns
   integer,  intent(in)  :: k                    ! Layer index

   real(8), intent(in)  :: p_in          ! Pressure [Pa]
   real(8), intent(in)  :: T0_in         ! Temperature [K]
   real(8), intent(in)  :: qv0_in        ! Grid-mean water vapor [kg/kg]
   real(8), intent(in)  :: ql0_in        ! Grid-mean LWC [kg/kg]
   real(8), intent(in)  :: qi0_in        ! Grid-mean IWC [kg/kg]
   real(8), intent(in)  :: ni0_in

   real(8), intent(in)  :: a_dc_in       ! Deep cumulus cloud fraction
   real(8), intent(in)  :: ql_dc_in      ! In-deep cumulus LWC [kg/kg]
   real(8), intent(in)  :: qi_dc_in      ! In-deep cumulus IWC [kg/kg]
   real(8), intent(in)  :: a_sc_in       ! Shallow cumulus cloud fraction
   real(8), intent(in)  :: ql_sc_in      ! In-shallow cumulus LWC [kg/kg]
   real(8), intent(in)  :: qi_sc_in      ! In-shallow cumulus IWC [kg/kg]

   real(8), intent(in)  :: landfrac      ! Land fraction
   real(8), intent(in)  :: snowh         ! Snow depth (liquid water equivalent)

   real(8), intent(in)  :: rhmini_in 
   real(8), intent(in)  :: rhminl_in
   real(8), intent(in)  :: rhminl_adj_land_in
   real(8), intent(in)  :: rhminh_in     

   real(8), intent(out) :: T_out         ! Temperature [K]
   real(8), intent(out) :: qv_out        ! Grid-mean water vapor [kg/kg]
   real(8), intent(out) :: ql_out        ! Grid-mean LWC [kg/kg]
   real(8), intent(out) :: qi_out        ! Grid-mean IWC [kg/kg]

   real(8), intent(out) :: al_st_out     ! Liquid stratus fraction
   real(8), intent(out) :: ai_st_out     ! Ice stratus fraction
   real(8), intent(out) :: ql_st_out     ! In-stratus LWC [kg/kg]
   real(8), intent(out) :: qi_st_out     ! In-stratus IWC [kg/kg]
   real(8), intent(inout) ::  wv_para(PARASIZE) 

   ! Local variables

   integer i                                     ! Column    index

   real(8) p    
   real(8) T0   
   real(8) qv0    
   real(8) ql0    
   real(8) qi0    
   real(8) a_dc   
   real(8) ql_dc  
   real(8) qi_dc  
   real(8) a_sc   
   real(8) ql_sc  
   real(8) qi_sc  
   real(8) esat0  
   real(8) qsat0  
   real(8) U0     
   real(8) U0_nc  
   real(8) G0_nc
   real(8) al0_st_nc            
   real(8) al0_st
   real(8) ai0_st_nc            
   real(8) ai0_st               
   real(8) a0_st               
   real(8) ql0_nc
   real(8) qi0_nc
   real(8) qc0_nc
   real(8) ql0_st
   real(8) qi0_st
   real(8) qc0_st
   real(8) T   
   real(8) qv    
   real(8) ql    
   real(8) qi
   real(8) ql_st
   real(8) qi_st
   real(8) es  
   real(8) qs  
   real(8) esat_in  
   real(8) qsat_in  
   real(8) U0_in  
   real(8) al0_st_nc_in
   real(8) ai0_st_nc_in
   real(8) G0_nc_in
   integer  idxmod 
   real(8) U
   real(8) U_nc
   real(8) al_st_nc
   real(8) ai_st_nc
   real(8) G_nc
   real(8) a_st
   real(8) al_st
   real(8) ai_st
   real(8) Tmin0
   real(8) Tmax0
   real(8) Tmin
   real(8) Tmax
   integer caseid

   real(8) rhmini
   real(8) rhminl
   real(8) rhminl_adj_land
   real(8) rhminh
   real(8) latvap
   real(8) cpair
   real(8) latice 
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real(8), intent(in) :: rhmaxi_const
    real(8), intent(in) :: premib 
    real(8), intent(in) :: premit 
    real(8), intent(in) :: icecrit 
    integer, intent(in) :: iceopt 
   !logical,  parameter  :: CAMstfrac    = .false.    ! If .true. (.false.),
   !real(8), parameter  :: qlst_min     = 2.d-5   ! Minimum in-stratus LWC constraint [ kg/kg ]
   !real(8), parameter  :: qlst_max     = 3.d-3   ! Maximum in-stratus LWC constraint [ kg/kg ]

   logical CAMstfrac
   real(8) qlst_min
   real(8) qlst_max

   !call lwpf_start_test_p4() 
   CAMstfrac  = .false.    ! If .true. (.false.),
   qlst_min     = 2.d-5   ! Minimum in-stratus LWC constraint [ kg/kg ]
   qlst_max     = 3.d-3   ! Maximum in-stratus LWC constraint [ kg/kg ]
!xukai-------------------------------------------------------------------------------------
   
   latvap = wv_para(2)
   cpair  = wv_para(5)
   latice = wv_para(3)
   

   ! ---------------- !
   ! Main Computation ! 
   ! ---------------- !
    !tmpind(1) = k
    
    !call put_k(k)
   call qsat_water(T0_in, p_in, &
        esat_in, qsat_in, wv_para)


    !call put_tmp1(7, 0.d0)
   U0_in = qv0_in/qsat_in
   if( CAMstfrac ) then
       call astG_RHU(U0_in,p_in,qv0_in,landfrac,snowh,al0_st_nc_in,G0_nc_in,ncol,&
                     premib, premit, &
                     rhminl_in, rhminl_adj_land_in, rhminh_in)
   else
       call astG_PDF(U0_in,p_in,qv0_in,landfrac,snowh,al0_st_nc_in,G0_nc_in,ncol,&
                     premib, premit, &
                     rhminl_in, rhminl_adj_land_in, rhminh_in)
   endif
    !call lwpf_start_test_p6() 
   call aist_vector(qv0_in,T0_in,p_in,qi0_in,ni0_in,landfrac,snowh,ai0_st_nc_in,ncol,&
                    wv_para, premib, premit, icecrit, iceopt, &
                    rhmaxi_const, rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in)
!xukai delete
!   T_out = al0_st_nc_in 
!   qv_out = G0_nc_in
!   ql_out = ai0_st_nc_in
!   qi_out = 0.d0


 
    !call lwpf_stop_test_p6() 
    !call lwpf_stop_test_p6() 
   !do i = 1, ncol

      ! ---------------------- !
      ! Define local variables !
      ! ---------------------- !

      p   = p_in

      T0  = T0_in
      qv0 = qv0_in
      ql0 = ql0_in
      qi0 = qi0_in

      a_dc  = a_dc_in
      ql_dc = ql_dc_in
      qi_dc = qi_dc_in

      a_sc  = a_sc_in
      ql_sc = ql_sc_in
      qi_sc = qi_sc_in

      ql_dc = 0.d0
      qi_dc = 0.d0
      ql_sc = 0.d0
      qi_sc = 0.d0

      es  = esat_in 
      qs  = qsat_in 
 
      rhmini = rhmini_in     
      rhminl = rhminl_in     
      rhminl_adj_land = rhminl_adj_land_in     
      rhminh = rhminh_in     

      idxmod = 0
      caseid = -1

      ! ------------------------------------------------------------ !
      ! Force the grid-mean RH to be smaller than 1 if oversaturated !
      ! In order to be compatible with reduced 3x3 QQ, condensation  !
      ! should occur only into the liquid in gridmean_RH.            !
      ! ------------------------------------------------------------ !

      if( qv0 .gt. qs ) then
        !call put_tmp1(7, 1.d0)
          call gridmean_RH( i, k, p, T0, qv0, ql0, qi0,      &
                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, &
                            landfrac, snowh , &
                            wv_para)
          call qsat_water(T0, p, esat0, qsat0, wv_para)
          U0      = (qv0/qsat0)
          U0_nc   =  U0 
          if( CAMstfrac ) then
              call astG_RHU_single(U0_nc, p, qv0, landfrac, snowh, al0_st_nc, G0_nc, &
                premib, premit, & 
                 rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
          else
              call astG_PDF_single(U0_nc, p, qv0, landfrac, snowh, al0_st_nc, G0_nc, &
                premib, premit, & 
                 rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
          endif
          call aist_single(qv0,T0,p,qi0,landfrac,snowh,ai0_st_nc,&
                            wv_para, premib, premit, icecrit, iceopt, &
                           rhmaxi_const, rhmini, rhminl, rhminl_adj_land, rhminh)
          ai0_st  = (1.d0-a_dc-a_sc)*ai0_st_nc
          al0_st  = (1.d0-a_dc-a_sc)*al0_st_nc
          a0_st   = max(ai0_st,al0_st)         
          idxmod  = 1 


      else
          ai0_st  = (1.d0-a_dc-a_sc)*ai0_st_nc_in
          al0_st  = (1.d0-a_dc-a_sc)*al0_st_nc_in
            !call put_tmp1(7, 2.d0)
      endif    
      a0_st   = max(ai0_st,al0_st)         



  
      ! ----------------------- ! 
      ! Handling of input state !
      ! ----------------------- !

      ql0_nc  = max(0.d0,ql0-a_dc*ql_dc-a_sc*ql_sc)
      qi0_nc  = max(0.d0,qi0-a_dc*qi_dc-a_sc*qi_sc)
      qc0_nc  = ql0_nc + qi0_nc 

      Tmin0 = T0 - (latvap/cpair)*ql0
      Tmax0 = T0 + ((latvap+latice)/cpair)*qv0

      ! ------------------------------------------------------------- !
      ! Do nothing and just exit if generalized in-stratus condensate !
      ! condition is satisfied. This includes the case I.             !
      ! For 4x4 liquid stratus, a0_st --> al0_st.                     ! 
      ! ------------------------------------------------------------- !
      if( ( ql0_nc .ge. qlst_min*al0_st ) .and. ( ql0_nc .le. qlst_max*al0_st ) ) then

            !call put_tmp1(7, 3.d0)
          ! ------------------ !
          ! This is the case I !
          ! ------------------ ! 
             T = T0
             qv = qv0
             ql = ql0
             qi = qi0
             caseid = 0

           !call put_tmp2(7, ql)
             goto 10

      else
         ! ----------------------------- !
         ! This is case II : Dense Cloud !
         ! ----------------------------- !   
         if( al0_st .eq. 0.d0 .and. ql0_nc .gt. 0.d0 ) then
             ! ------------------------------------- !
             ! Compute hypothetical full evaporation !
             ! ------------------------------------- !
             T  = Tmin0
             qv = qv0 + ql0 
             call qsat_water(T, p, es, qs, wv_para)
             U  = qv/qs
             U_nc = U  
             if( CAMstfrac ) then
                 call astG_RHU_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
                premib, premit, & 
                    rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
             else
                 call astG_PDF_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
                premib, premit, & 
                    rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
             endif
             al_st = (1.d0-a_dc-a_sc)*al_st_nc  
             caseid = 0



             if( al_st .eq. 0.d0 ) then
                !call put_tmp1(7, 5.d0)
                 ql = 0.d0
                 qi = qi0
                 idxmod = 1
                 caseid = 1
            !call put_tmp2(7, ql)
                 goto 10
             else
                !call put_tmp1(7, 6.d0)
                 ! ------------------------------------------- !
                 ! Evaporate until qc_st decreases to qlst_max !
                 ! ------------------------------------------- !
                 Tmin = Tmin0
                 Tmax = T0 
                 call instratus_core( i, k, p,                              &
                                      T0, qv0, ql0, 0.d0,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_max, Tmin, Tmax, landfrac, snowh, &
                                      rhminl, rhminl_adj_land, rhminh,             &
                                      T, qv, ql, qi ,&
                                      wv_para, premib, premit)   
                 idxmod = 1
                 caseid = 2
            !call put_tmp2(7, ql)
                 goto 10
             endif
         ! ------------------------------ !
         ! This is case III : Empty Cloud !
         ! ------------------------------ !  
         elseif( al0_st .gt. 0.d0 .and. ql0_nc .eq. 0.d0 ) then
                !call put_tmp1(7, 7.d0)
              ! ------------------------------------------ ! 
              ! Condense until qc_st increases to qlst_min !
              ! ------------------------------------------ !
              Tmin = Tmin0
              Tmax = Tmax0  
              call instratus_core(  i, k, p,                              &
                                   T0, qv0, ql0, 0.d0,                         &
                                   a_dc, ql_dc, qi_dc,                          &
                                   a_sc, ql_sc, qi_sc, ai0_st,                  &
                                   qlst_min, Tmin, Tmax, landfrac, snowh, &
                                   rhminl, rhminl_adj_land, rhminh,             &
                                   T, qv, ql, qi , &
                                   wv_para, premib, premit)   
              idxmod = 1 
              caseid = 3
            !call put_tmp2(7, ql)
              goto 10
         ! --------------- !
         ! This is case IV !
         ! --------------- !   
         elseif( al0_st .gt. 0.d0 .and. ql0_nc .gt. 0.d0 ) then

             if( ql0_nc .gt. qlst_max*al0_st ) then
                !call put_tmp1(7, 8.d0)
                 ! --------------------------------------- !
                 ! Evaporate until qc_st drops to qlst_max !
                 ! --------------------------------------- !
                 Tmin = Tmin0
                 Tmax = Tmax0
                 call instratus_core( i, k, p,                              &
                                      T0, qv0, ql0, 0.d0,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_max, Tmin, Tmax, landfrac, snowh, &
                                      rhminl, rhminl_adj_land, rhminh,             &
                                      T, qv, ql, qi,&
                                      wv_para, premib, premit)   
                 idxmod = 1
                 caseid = 4
            !call put_tmp2(7, ql)
                 goto 10
             elseif( ql0_nc .lt. qlst_min*al0_st ) then
                !call put_tmp1(7, 9.d0)
                 ! -------------------------------------------- !
                 ! Condensate until qc_st increases to qlst_min !
                 ! -------------------------------------------- !
                 Tmin = Tmin0
                 Tmax = Tmax0 
                 call instratus_core( i, k, p,                              &
                                      T0, qv0, ql0, 0.d0,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_min, Tmin, Tmax, landfrac, snowh, & 
                                      rhminl, rhminl_adj_land, rhminh,             &
                                      T, qv, ql, qi ,&
                                      wv_para, premib, premit)   
                 idxmod = 1
                 caseid = 5
            !call put_tmp2(7, ql)
                 goto 10
             else
                !call put_tmp1(7, 10.d0)
            !call put_tmp2(7, ql)
                 ! ------------------------------------------------ !
                 ! This case should not happen. Issue error message !
                 ! ------------------------------------------------ !
                 !write(iulog,*) 'Impossible case1 in instratus_condensate' 
                 wv_para(PARASIZE) = 1.0d0
                 !call endrun()
             endif
         ! ------------------------------------------------ !                   
         ! This case should not happen. Issue error message !
         ! ------------------------------------------------ !    
         else
                !call put_tmp1(7, 11.d0)
            !call put_tmp2(7, i, k, ql)
             !write(iulog,*) 'Impossible case2 in instratus_condensate' 
             !write(iulog,*)  al0_st, a_sc, a_dc
             !write(iulog,*)  1000*ql0_nc, 1000*(ql0+qi0)
              wv_para(PARASIZE) = 2.0d0
             !call endrun()
         endif
      endif

10 continue   
!call set_wvflag(4)
!call put_tmp1(4,ql) 
!call put_tmp2(4,qi) 
!call set_wvflag(0)


   ! -------------------------------------------------- !
   ! Force final energy-moisture conserving consistency !
   ! -------------------------------------------------- !

     qi = qi0

     if( idxmod .eq. 1 ) then
         call aist_single(qv,T,p,qi,landfrac,snowh,ai_st_nc,&
                           wv_para, premib, premit, icecrit, iceopt, &
                          rhmaxi_const, rhmini, rhminl, rhminl_adj_land, rhminh)
         ai_st = (1.d0-a_dc-a_sc)*ai_st_nc
         call qsat_water(T, p, es, qs, wv_para)
         U     = (qv/qs)
         U_nc  =  U
         if( CAMstfrac ) then
             call astG_RHU_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
                premib, premit, & 
                rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
         else
             call astG_PDF_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
                premib, premit, & 
                rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
         endif
         al_st = (1.d0-a_dc-a_sc)*al_st_nc
           !call put_tmp1(7, 12.d0)
           !call put_tmp2(7, al_st)
 

     else
         ai_st  = (1.d0-a_dc-a_sc)*ai0_st_nc_in
         al_st  = (1.d0-a_dc-a_sc)*al0_st_nc_in
           !call put_tmp1(7, 13.d0)
           !call put_tmp2(7,  al_st)
     endif

     a_st  = max(ai_st,al_st)

     if( al_st .eq. 0.d0 ) then
         ql_st = 0.d0
           !call put_tmp1(7, 14.d0)
           !call put_tmp2(7,  ql_st)
     else
         ql_st = ql/al_st
         ql_st = min(qlst_max,max(qlst_min,ql_st)) ! PJR
           !call put_tmp1(7, 15.d0)
           !call put_tmp2(7,  ql_st)

     endif
     if( ai_st .eq. 0.d0 ) then
         qi_st = 0.d0
     else
         qi_st = qi/ai_st
     endif

     qi    = ai_st*qi_st
     ql    = al_st*ql_st

           !call put_tmp1(7, 16.d0)
           !call put_tmp2(7,  ql)
     T     = T0 - (latvap/cpair)*(ql0-ql) - ((latvap+latice)/cpair)*(qi0-qi)
     qv    = qv0 + ql0 - ql + qi0 - qi

   ! -------------- !
   ! Send to output !
   ! -------------- !


   T_out  = T
   qv_out = qv
   ql_out = ql
   qi_out = qi

   al_st_out = al_st
   ai_st_out = ai_st
   ql_st_out = ql_st
   qi_st_out = qi_st

!call set_wvflag(4)
!call put_tmp1(4,ql) 
!call put_tmp2(4,qi) 
!call set_wvflag(0)

   !enddo 

   return
   end subroutine instratus_condensate

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine instratus_core( icol, k, p,                      &
                              T0, qv0, ql0, qi0,                      &
                              a_dc, ql_dc, qi_dc,                     & 
                              a_sc, ql_sc, qi_sc, ai_st,              &
                              qcst_crit, Tmin, Tmax, landfrac, snowh, &
                              rhminl, rhminl_adj_land, rhminh,        &
                              T, qv, ql, qi ,&
                              wv_para, premib, premit)

   ! ------------------------------------------------------ !
   ! Subroutine to find saturation equilibrium state using  ! 
   ! a Newton iteration method, so that 'qc_st = qcst_crit' !
   ! is satisfied.                                          !
   ! ------------------------------------------------------ !

   implicit none
   !integer,  intent(in)  :: lchnk      ! Chunk identifier
   integer,  intent(in)  :: icol       ! Number of atmospheric columns
   integer,  intent(in)  :: k          ! Layer index

   real(8), intent(in)  :: p          ! Pressure [Pa]
   real(8), intent(in)  :: T0         ! Temperature [K]
   real(8), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
   real(8), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
   real(8), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]

   real(8), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
   real(8), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(8), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(8), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
   real(8), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(8), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(8), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)

   real(8), intent(in)  :: Tmin       ! Minimum temperature system can have [K]
   real(8), intent(in)  :: Tmax       ! Maximum temperature system can have [K]
   real(8), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
   real(8), intent(in)  :: landfrac   ! Land fraction
   real(8), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)

   real(8), intent(in)  :: rhminl
   real(8), intent(in)  :: rhminl_adj_land
   real(8), intent(in)  :: rhminh

   real(8), intent(out) :: T          ! Temperature [K]
   real(8), intent(out) :: qv         ! Grid-mean water vapor [kg/kg]
   real(8), intent(out) :: ql         ! Grid-mean LWC [kg/kg]
   real(8), intent(out) :: qi         ! Grid-mean IWC [kg/kg]
   real(8), intent(in) ::  wv_para(PARASIZE)

   ! Local variables

   integer i                           ! Iteration index

   real(8) muQ0, muQ
   real(8) ql_nc0, qi_nc0, qc_nc0, qc_nc    
   real(8) fice0, fice    
   real(8) ficeg0, ficeg   
   real(8) esat0
   real(8) qsat0
   real(8) dqcncdt, dastdt, dUdt
   real(8) alpha, beta
   real(8) U, U_nc
   real(8) al_st_nc, G_nc
   real(8) al_st

   ! Variables for root-finding algorithm

   integer j                          
   real(8)  x1, x2
   real(8)  rtsafe
   real(8)  df, dx, dxold, f, fh, fl, temp, xh, xl
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !real(8), parameter :: xacc = 1.d-3
    real(8), intent(in) :: premib
    real(8), intent(in) :: premit
 
   !real(8), parameter  :: qlst_min     = 2.d-5   ! Minimum in-stratus LWC constraint [ kg/kg ]
   !real(8), parameter  :: qlst_max     = 3.d-3   ! Maximum in-stratus LWC constraint [ kg/kg ]
   real(8) xacc 
   real(8) qlst_min
   real(8) qlst_max

   xacc = 1.d-3
   qlst_min     = 2.d-5   ! Minimum in-stratus LWC constraint [ kg/kg ]
   qlst_max     = 3.d-3   ! Maximum in-stratus LWC constraint [ kg/kg ]

!xukai-------------------------------------------------------------------------------------
   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   ql_nc0 = max(0.d0,ql0-a_dc*ql_dc-a_sc*ql_sc)
   qi_nc0 = max(0.d0,qi0-a_dc*qi_dc-a_sc*qi_sc)
   qc_nc0 = max(0.d0,ql0+qi0-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
   fice0  = 0.d0
   ficeg0 = 0.d0
   muQ0   = 1.d0

   ! ------------ !
   ! Root finding !
   ! ------------ !

   x1 = Tmin
   x2 = Tmax
   call funcd_instratus( x1, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                         qcst_crit, landfrac, snowh,                    &
                         rhminl, rhminl_adj_land, rhminh,               &
                         fl, df, qc_nc, fice, al_st ,&
                         wv_para, premib, premit)
   call funcd_instratus( x2, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                         qcst_crit, landfrac, snowh,                    &
                         rhminl, rhminl_adj_land, rhminh,               &
                         fh, df, qc_nc, fice, al_st ,&
                         wv_para, premib, premit)
   if((fl > 0.d0 .and. fh > 0.d0) .or. (fl < 0.d0 .and. fh < 0.d0)) then
       call funcd_instratus( T0, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                             a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                             qcst_crit, landfrac, snowh,                    &
                             rhminl, rhminl_adj_land, rhminh,               &
                             fl, df, qc_nc, fice, al_st , &
                             wv_para, premib, premit)
       rtsafe = T0 
       goto 10       
   endif
   if( fl == 0.d0) then
           rtsafe = x1
           goto 10
   elseif( fh == 0.d0) then
           rtsafe = x2
           goto 10
   elseif( fl < 0.d0) then
           xl = x1
           xh = x2
   else
           xh = x1
           xl = x2
   end if
   rtsafe = 0.5d0*(x1+x2)
   dxold = abs(x2-x1)
   dx = dxold
   call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
                         qcst_crit, landfrac, snowh,                        &
                         rhminl, rhminl_adj_land, rhminh,                   &
                         f, df, qc_nc, fice, al_st ,&
                         wv_para, premib, premit)
   do j = 1, 20
      if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0.d0 .or. abs(2.0d0*f) > abs(dxold*df) ) then
           dxold = dx
           dx = 0.5d0*(xh-xl)
           rtsafe = xl + dx
           if(xl == rtsafe) goto 10
      else
           dxold = dx
           dx = f/df
           temp = rtsafe
           rtsafe = rtsafe - dx
           if (temp == rtsafe) goto 10
      end if
    ! if(abs(dx) < xacc) goto 10
      call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
                            qcst_crit, landfrac, snowh,                        &
                            rhminl, rhminl_adj_land, rhminh,                   &
                            f, df, qc_nc, fice, al_st ,&
                            wv_para, premib, premit)
    ! Sep.21.2010. Sungsu modified to enhance convergence and guarantee 'qlst_min <  qlst < qlst_max'.
      if( qcst_crit < 0.5d0 * ( qlst_min + qlst_max ) ) then
          if( ( qc_nc*(1.d0-fice) .gt.          qlst_min*al_st .and. &
                qc_nc*(1.d0-fice) .lt. 1.1d0 * qlst_min*al_st ) ) goto 10
      else
          if( ( qc_nc*(1.d0-fice) .gt. 0.9d0 * qlst_max*al_st .and. &
                qc_nc*(1.d0-fice) .lt.          qlst_max*al_st ) ) goto 10
      endif
      if(f < 0.d0) then
          xl = rtsafe
      else
          xh = rtsafe
      endif

   enddo

10 continue

   ! ------------------------------------------- !
   ! Final safety check before sending to output !
   ! ------------------------------------------- !

   qc_nc = max(0.d0,qc_nc)

   T  = rtsafe
   ql = qc_nc*(1.d0-fice) + a_dc*ql_dc + a_sc*ql_sc
   qi = qc_nc*fice + a_dc*qi_dc + a_sc*qi_sc
   qv = qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc))
   qv = max(qv,1.d-12) 

   return
   end subroutine instratus_core

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine funcd_instratus( T, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0,   &
                               a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,  &
                               qcst_crit, landfrac, snowh,                     &
                               rhminl, rhminl_adj_land, rhminh,                &
                               f, fg, qc_nc, fice, al_st ,&
                               wv_para, premib, premit) 

   ! --------------------------------------------------- !
   ! Subroutine to find function value and gradient at T !
   ! --------------------------------------------------- !

   implicit none

   real(8), intent(in)  :: T          ! Iteration temperature [K]

   real(8), intent(in)  :: p          ! Pressure [Pa]
   real(8), intent(in)  :: T0         ! Initial temperature [K]
   real(8), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
   real(8), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
   real(8), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]
   real(8), intent(in)  :: fice0      ! 
   real(8), intent(in)  :: muQ0       ! 
   real(8), intent(in)  :: qc_nc0     ! 

   real(8), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
   real(8), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(8), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(8), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
   real(8), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(8), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(8), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)

   real(8), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
   real(8), intent(in)  :: landfrac   ! Land fraction
   real(8), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)

   real(8), intent(in)  :: rhminl
   real(8), intent(in)  :: rhminl_adj_land
   real(8), intent(in)  :: rhminh

   real(8), intent(out) :: f          ! Value of minimization function at T
   real(8), intent(out) :: fg         ! Gradient of minimization function 
   real(8), intent(out) :: qc_nc      !
   real(8), intent(out) :: al_st      !
   real(8), intent(out) :: fice       !
   real(8), intent(in)  :: wv_para(PARASIZE)

   ! Local variables

   real(8) es
   real(8) qs
   real(8) dqsdT
   real(8) dqcncdt
   real(8) alpha
   real(8) beta
   real(8) U
   real(8) U_nc
   real(8) al_st_nc
   real(8) G_nc
   real(8) dUdt
   real(8) dalstdt
   real(8) qv
   real(8) cpair
   real(8) latvap 
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real(8), intent(in) :: premib
    real(8), intent(in) :: premit
   !logical,  parameter  :: CAMstfrac    = .false.    ! If .true. (.false.),
   logical  :: CAMstfrac  
   real(8) :: pow_base_tmp
   CAMstfrac = .false.    ! If .true. (.false.),

!xukai-------------------------------------------------------------------------------------
   cpair = wv_para(5)
   latvap = wv_para(2)
   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   call qsat_water(T, p, es, qs, wv_para, dqsdt=dqsdt)

   fice    = fice0 
   qc_nc   = (cpair/latvap)*(T-T0)+muQ0*qc_nc0       
   dqcncdt = (cpair/latvap) 
   qv      = (qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc)))
   alpha   = (1.d0/qs)
   pow_base_tmp = 2.d0
   beta    = (qv/pow_c(qs, pow_base_tmp))*dqsdT 

   U      =  (qv/qs)
   U_nc   =   U
   if( CAMstfrac ) then
       call astG_RHU_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
                premib, premit, & 
          rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
   else
       call astG_PDF_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
                premib, premit, & 
          rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
   endif
   al_st   =  (1.d0-a_dc-a_sc)*al_st_nc 
   dUdt    = -(alpha*dqcncdt+beta)
   dalstdt =  (1.d0/G_nc)*dUdt
   if( U_nc .eq. 1.d0 ) dalstdt = 0.d0

   f  = qc_nc   - qcst_crit*al_st
   fg = dqcncdt - qcst_crit*dalstdt

   return
   end subroutine funcd_instratus

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine gridmean_RH( icol, k, p, T, qv, ql, qi,       &
                           a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, &
                           landfrac, snowh ,&
                           wv_para)

   ! ------------------------------------------------------------- !
   ! Subroutine to force grid-mean RH = 1 when RH > 1              !
   ! This is condensation process similar to instratus_condensate. !
   ! During condensation, we assume 'fice' is maintained in this   !
   ! verison for MG not for RK.                                    !
   ! ------------------------------------------------------------- !

   implicit none
   !integer,  intent(in)    :: lchnk      ! Chunk identifier
   integer,  intent(in)    :: icol       ! Number of atmospheric columns
   integer,  intent(in)    :: k          ! Layer index

   real(8), intent(in)    :: p          ! Pressure [Pa]
   real(8), intent(inout) :: T          ! Temperature [K]
   real(8), intent(inout) :: qv         ! Grid-mean water vapor [kg/kg]
   real(8), intent(inout) :: ql         ! Grid-mean LWC [kg/kg]
   real(8), intent(inout) :: qi         ! Grid-mean IWC [kg/kg]

   real(8), intent(in)    :: a_dc       ! Deep cumulus cloud fraction
   real(8), intent(in)    :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(8), intent(in)    :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(8), intent(in)    :: a_sc       ! Shallow cumulus cloud fraction
   real(8), intent(in)    :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(8), intent(in)    :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(8), intent(in)    :: landfrac   ! Land fraction
   real(8), intent(in)    :: snowh      ! Snow depth (liquid water equivalent)
   real(8), intent(in)    :: wv_para(PARASIZE)

   !integer, intent(in)    :: pver
   !integer, intent(in)    :: pverp

   ! Local variables

   integer m                             ! Iteration index

   real(8)  ql_nc0, qi_nc0, qc_nc0
   real(8)  Tscale
   real(8)  Tc, qt, qc, dqcdt, qc_nc    
   real(8)  es, qs, dqsdT
   real(8)  al_st_nc, G_nc
   real(8)  f, fg
   !real(8), parameter :: xacc = 1.d-3
   real(8) :: xacc
   real(8) cpair
   real(8) latvap 
   cpair = wv_para(5)
   latvap = wv_para(2)
   xacc = 1.d-3


   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   ql_nc0 = max(0.d0,ql-a_dc*ql_dc-a_sc*ql_sc)
   qi_nc0 = max(0.d0,qi-a_dc*qi_dc-a_sc*qi_sc)
   qc_nc0 = max(0.d0,ql+qi-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
   Tc    = T - (latvap/cpair)*ql
   qt    = qv + ql

   do m = 1, 20
      call qsat_water(T, p, es, qs, wv_para, dqsdt=dqsdT)
      Tscale = latvap/cpair
      qc     = (T-Tc)/Tscale
      dqcdt  = 1.d0/Tscale
      f      = qs + qc - qt 
      fg     = dqsdT + dqcdt
      fg     = sign(1.d0,fg)*max(1.d-10,abs(fg))
    ! Sungsu modified convergence criteria to speed up convergence and guarantee RH <= 1.
      if( qc .ge. 0.d0 .and. ( qt - qc ) .ge. 0.999d0*qs .and. ( qt - qc ) .le. 1.d0*qs ) then
          goto 10
      endif
      T = T - f/fg
   enddo
 ! write(iulog,*) 'Convergence in gridmean_RH is not reached. RH = ', ( qt - qc ) / qs
10 continue

   call qsat_water(T, p, es, qs, wv_para)
 ! Sungsu modified 'qv = qs' in consistent with the modified convergence criteria above.
   qv = min(qt,qs) ! Modified
   ql = qt - qv
   T  = Tc + (latvap/cpair)*ql

   return
   end subroutine gridmean_RH

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine positive_moisture( ncol, dt, qvmin, qlmin, qimin, dp, &
                                 qv,   ql, qi,    t,     qvten, &
                                 qlten,    qiten, tten,  do_cldice, &
                                 wv_para, pver, pverp, top_lev)

   ! ------------------------------------------------------------------------------- !
   ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
   ! force them to be larger than minimum value by (1) condensating water vapor      !
   ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
   ! layer. '2.d0' is multiplied to the minimum values for safety.                  !
   ! Update final state variables and tendencies associated with this correction.    !
   ! If any condensation happens, update (s,t) too.                                  !
   ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
   ! input tendencies.                                                               !
   ! Be careful the order of k : '1': top layer, 'pver' : near-surface layer         ! 
   ! ------------------------------------------------------------------------------- !

   implicit none
   integer,  intent(in)     :: pver 
   integer,  intent(in)     :: pverp 
   integer, intent(in)    :: top_lev
   integer,  intent(in)     :: ncol
   real(8), intent(in)     :: dt
   real(8), intent(in)     :: dp(top_lev:pver), qvmin(top_lev:pver), qlmin(top_lev:pver), qimin(top_lev:pver)
   real(8), intent(inout)  :: qv(top_lev:pver), ql(top_lev:pver), qi(top_lev:pver), t(top_lev:pver)
   real(8), intent(out)    :: qvten(top_lev:pver), qlten(top_lev:pver), qiten(top_lev:pver), tten(top_lev:pver)
   logical, intent(in)      :: do_cldice
   integer   i, k
   real(8)  dql, dqi, dqv, sum, aa, dum 
   real(8), intent(inout)    :: wv_para(PARASIZE)
   real(8) cpair
   real(8) latvap 
   real(8) latice
   integer :: gotoflag
   integer :: gotoflag_other
   integer :: core_oe
   real(8) :: dp_tmp 

   
   cpair = wv_para(5)
   latvap = wv_para(2)
   latice = wv_para(3)

   gotoflag = 0

!   tten(:pver)  = 0.d0
!   qvten(:pver) = 0.d0
!   qlten(:pver) = 0.d0
!   qiten(:pver) = 0.d0

   tten(top_lev:pver)  = 0.d0
   qvten(top_lev:pver) = 0.d0
   qlten(top_lev:pver) = 0.d0
   qiten(top_lev:pver) = 0.d0

    call get_core_oe(core_oe)
    !dp_tmp = dp(pver)
!
!need DMA commmunicate
   !do i = 1, ncol
      do k = top_lev, pver
         if( qv(k) .lt. qvmin(k) .or. ql(k) .lt. qlmin(k) .or. qi(k) .lt. qimin(k) ) then
            gotoflag = 10
!             goto 10
         endif
      enddo
      call regcomm_gotoflag(gotoflag, gotoflag_other)
      if(gotoflag == 10 .or. gotoflag_other == 10) then
          goto 10
      endif

      goto 11

   10 continue
    if(core_oe == 1) then
             call regcomm_d021(dp_tmp)
             call regcomm_d021(dqv)
             qv(top_lev)    = qv(top_lev)    - dqv*dp_tmp/dp(top_lev)
             qvten(top_lev) = qvten(top_lev) - dqv*dp_tmp/dp(top_lev)/dt
    endif

      do k = top_lev, pver    ! From the top to the 1st (lowest) layer from the surface
        !call put_k(k)
         dql = max(0.d0,1.d0*qlmin(k)-ql(k))

         if (do_cldice) then
         dqi = max(0.d0,1.d0*qimin(k)-qi(k))
            !call put_tmp1(6, 77.d0)
         else
           dqi = 0.d0
            !call put_tmp1(6, 66.d0)
         end if

         qlten(k) = qlten(k) +  dql/dt
         qiten(k) = qiten(k) +  dqi/dt
         qvten(k) = qvten(k) - (dql+dqi)/dt
         tten(k)  = tten(k)  + (latvap/cpair)*(dql/dt) + ((latvap+latice)/cpair)*(dqi/dt)
         ql(k)    = ql(k) + dql
         qi(k)    = qi(k) + dqi
         qv(k)    = qv(k) - dql - dqi
         t(k)     = t(k)  + (latvap * dql + (latvap+latice) * dqi)/cpair
         dqv        = max(0.d0,1.d0*qvmin(k)-qv(k))
         qvten(k) = qvten(k) + dqv/dt
         qv(k)    = qv(k)    + dqv
         if( k .ne. pver ) then 
             qv(k+1)    = qv(k+1)    - dqv*dp(k)/dp(k+1)
             qvten(k+1) = qvten(k+1) - dqv*dp(k)/dp(k+1)/dt
        else
            if(core_oe == 0) then
                  dp_tmp = dp(k)
                 call regcomm_d021(dp_tmp)
                 call regcomm_d021(dqv)
            endif
         endif

         qv(k) = max(qv(k),qvmin(k))
         ql(k) = max(ql(k),qlmin(k))
         qi(k) = max(qi(k),qimin(k))
      end do
      ! Extra moisture used to satisfy 'qv(pver)=qvmin' is proportionally 
      ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
      ! preserves column moisture. 

      call regcomm_d120(dqv)
      if( dqv .gt. 1.d-20 ) then
          sum = 0.d0
          do k = top_lev, pver
             if( qv(k) .gt. 2.d0*qvmin(k) ) sum = sum + qv(k)*dp(k)
          enddo
          call regcomm_sum(sum)
          aa = dqv*dp(pver)/max(1.d-20,sum)
          call regcomm_d120(aa)
          if( aa .lt. 0.5d0 ) then
              do k = top_lev, pver
                 if( qv(k) .gt. 2.d0*qvmin(k) ) then
                     dum        = aa*qv(k)
                     qv(k)    = qv(k) - dum
                     qvten(k) = qvten(k) - dum/dt
                 endif
              enddo 
          else 
			  wv_para(PARASIZE) = 10.d0
              
          endif
      endif 
11 continue
   !enddo
   return

   end subroutine positive_moisture

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !
      SUBROUTINE gaussj(a,n,np,b,m,mp, wv_para)
      implicit none
      INTEGER, intent(in) :: m,mp,n,np
      !,NMAX
      real(8), intent(inout) :: a(np,np),b(np,mp)
      real(8) aa(np,np),bb(np,mp)
      !PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,ii,jj,indxc(8),indxr(8),ipiv(8)
      real(8) big,dum,pivinv
      real(8), intent(inout) :: wv_para(PARASIZE)

      aa(:,:) = a(:,:)
      bb(:,:) = b(:,:)

      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                !write(iulog,*) 'singular matrix in gaussj 1'
                !do ii = 1, np
                !do jj = 1, np
                !   write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
                !end do
                !end do   
                !call endrun
                wv_para(PARASIZE) = 3.d0
                goto 25
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) then
            !write(iulog,*) 'singular matrix in gaussj 2'
            !do ii = 1, np
            !do jj = 1, np
            !   write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
            !end do
            !end do   
            !call endrun

            wv_para(PARASIZE) = 4.d0
            goto 25
        endif 
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue


25 continue
      return
      end subroutine gaussj

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !



      SUBROUTINE gaussjGG(a,n,np,b,m,mp, wv_para)
      INTEGER, intent(in) ::  m,mp,n,np
      !,NMAX
      real(8), intent(inout) ::  a(np,np),b(np,mp)
      real(8), intent(inout) ::  wv_para(PARASIZE) 
      real(8) aa(np,np),bb(np,mp)
      !PARAMETER (NMAX=8)
      INTEGER i,icol,irow,j,k,l,ll,ii,jj,indxc(8),indxr(8),ipiv(8)
      real(8) big,dum,pivinv

      aa(:,:) = a(:,:)
      bb(:,:) = b(:,:)

      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                !write(iulog,*) 'singular matrix in gaussj 1'
                !do ii = 1, np
                !do jj = 1, np
                !   write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
                !end do
                !end do   
                !call endrun
                wv_para(PARASIZE) = 3.d0
                goto 25
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) then
            !write(iulog,*) 'singular matrix in gaussj 2'
            !do ii = 1, np
            !do jj = 1, np
            !   write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
            !end do
            !end do   
            !call endrun
            wv_para(PARASIZE) = 4.d0
            goto 25
        endif 
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue

25 continue

      return
      end subroutine gaussjGG


  ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

end module cldwat2m_macro_cpe
