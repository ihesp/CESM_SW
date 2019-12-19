!*************************************************************************
! File Name: src.cpe/cldwat2m_macro_cpe.F90
! Author: Xu Kai
! Created Time: 2018年12月13日 星期四 10时30分38秒
! ************************************************************************
#define PARASIZE 12

       subroutine mmacro_pcond_cpe(  ncol, dt, p, dp, &
            T0, qv0, ql0, qi0, nl0, ni0, &
            A_T, A_qv, A_ql, A_qi, A_nl, A_ni, &
            C_T, C_qv, C_ql, C_qi, C_nl, C_ni, C_qlst, &
            D_T, D_qv, D_ql, D_qi, D_nl, D_ni, &
            a_cud, a_cu0, clrw_old, clri_old, landfrac, snowh, & 
            tke, qtl_flx, qti_flx, cmfr_det, qlr_det, qir_det, &
            s_tendout, qv_tendout, ql_tendout, qi_tendout, nl_tendout, ni_tendout, &
            qme, qvadj, qladj, qiadj, qllim, qilim, &
            cld, al_st_star, ai_st_star, ql_st_star, qi_st_star, &
            d_rhmin_liq_PBL, d_rhmin_ice_PBL, d_rhmin_liq_det, d_rhmin_ice_det,&
            wv_para, pver, pverp, top_lev, &
            i_rhmini, i_rhminl, rhmini_const, rhmaxi_const, rhminl_const, rhminl_adj_land_const, rhminh_const, &
            qmins1, qmins2, qmins3, premib, premit,  &
            icecrit, iceopt, do_cldice_int)
    
       !use constituents, only : qmin, cnst_get_ind
       use wv_saturation_cpe_xk,    only : findsp_vc, qsat_water
       use  cldwat2m_macro_cpe,  only : rhcrit_calc, positive_moisture, &
                                        instratus_condensate, gaussj
       use cldfrc2m_cpe,         only: astG_PDF_single, astG_PDF, astG_RHU_single, &
                               astG_RHU, aist_single, aist_vector


        !use wv_sat_methods_cpe, only: tmp1, tmp2, tmpind, wvflag
        implicit none
    
       integer   icol
       !integer,  intent(in)    :: lchnk                        ! Chunk number
       integer,  intent(in)    :: pver                         
       integer,  intent(in)    :: pverp                         
       integer,  intent(in)    :: ncol                         
        integer, intent(in) :: top_lev 
       !integer,  intent(in)    :: dim1                         
       !integer,  intent(in)    :: dim2                         
    
       ! Input-Output variables
    
       real(8), intent(inout) :: T0(top_lev:pver)               
       real(8), intent(inout) :: qv0(top_lev:pver)              
       real(8), intent(inout) :: ql0(top_lev:pver)              
       real(8), intent(inout) :: qi0(top_lev:pver)              
       real(8), intent(inout) :: nl0(top_lev:pver)              
       real(8), intent(inout) :: ni0(top_lev:pver)              
    
       ! Input variables
    
       real(8), intent(in)    :: dt                           
       real(8), intent(in)    :: p(top_lev:pver)                
       real(8), intent(in)    :: dp(top_lev:pver)               
    
       real(8), intent(in)    :: A_T(top_lev:pver)              
       real(8), intent(in)    :: A_qv(top_lev:pver)             
       real(8), intent(in)    :: A_ql(top_lev:pver)             
       real(8), intent(in)    :: A_qi(top_lev:pver)             
       real(8), intent(in)    :: A_nl(top_lev:pver)             
       real(8), intent(in)    :: A_ni(top_lev:pver)             
    
       real(8), intent(in)    :: C_T(top_lev:pver)              
       real(8), intent(in)    :: C_qv(top_lev:pver)             
       real(8), intent(in)    :: C_ql(top_lev:pver)             
       real(8), intent(in)    :: C_qi(top_lev:pver)             
       real(8), intent(in)    :: C_nl(top_lev:pver)             
       real(8), intent(in)    :: C_ni(top_lev:pver)             
       real(8), intent(in)    :: C_qlst(top_lev:pver)           
                                                               
    
       real(8), intent(in)    :: D_T(top_lev:pver)              
       real(8), intent(in)    :: D_qv(top_lev:pver)             
       real(8), intent(in)    :: D_ql(top_lev:pver)             
       real(8), intent(in)    :: D_qi(top_lev:pver)             
       real(8), intent(in)    :: D_nl(top_lev:pver)             
       real(8), intent(in)    :: D_ni(top_lev:pver)             
    
       real(8), intent(in)    :: a_cud(top_lev:pver)            
       real(8), intent(in)    :: a_cu0(top_lev:pver)            
    
       real(8), intent(in)    :: clrw_old(top_lev:pver)         
       real(8), intent(in)    :: clri_old(top_lev:pver)         
       real(8), intent(in)    :: tke(top_lev:pverp)                     
       real(8), intent(in)    :: qtl_flx(top_lev:pverp)                 
       real(8), intent(in)    :: qti_flx(top_lev:pverp)                 
       real(8), intent(in)    :: cmfr_det(top_lev:pver)                
       real(8), intent(in)    :: qlr_det(top_lev:pver)                 
       real(8), intent(in)    :: qir_det(top_lev:pver)                 
    
       real(8), intent(in)    :: landfrac              
       real(8), intent(in)    :: snowh                 
       real(8), intent(inout)    :: wv_para(PARASIZE)                 

       integer,  intent(in)    :: do_cldice_int                    
       logical                 :: do_cldice                    

    
       
    
       real(8), target, intent(out)   :: s_tendout(top_lev:pver)        
       real(8), target, intent(out)   :: qv_tendout(top_lev:pver)       
       real(8), intent(out)   :: ql_tendout(top_lev:pver)       
       real(8), intent(out)   :: qi_tendout(top_lev:pver)       
       real(8), intent(out)   :: nl_tendout(top_lev:pver)       
       real(8), intent(out)   :: ni_tendout(top_lev:pver)       
    
       real(8), intent(out)   :: qme  (top_lev:pver)            
       real(8), intent(out)   :: qvadj(top_lev:pver)            
       real(8), intent(out)   :: qladj(top_lev:pver)            
       real(8), intent(out)   :: qiadj(top_lev:pver)            
       real(8), intent(out)   :: qllim(top_lev:pver)            
       real(8), intent(out)   :: qilim(top_lev:pver)            
    
       real(8), intent(out)   :: cld(top_lev:pver)              
       real(8), intent(out)   :: al_st_star(top_lev:pver)       
       real(8), intent(out)   :: ai_st_star(top_lev:pver)       
       real(8), intent(out)   :: ql_st_star(top_lev:pver)       
       real(8), intent(out)   :: qi_st_star(top_lev:pver)       
       !integer(8), intent(out) :: spa2
       integer(8) :: spa2
        
       !real(8), intent(inout) :: aa(dim1,dim1)
       !real(8), intent(inout) :: bb(dim1,dim2)
    
 
       ! --------------- !
       ! Local variables !
       ! --------------- !
       integer :: ixcldliq, ixcldice
     
       integer :: i, j, k, iter, ii, jj                        
    
       
    
       real(8) T(top_lev:pver)                                  
                                                               
       real(8) T1(top_lev:pver)                                 
       real(8) T_0(top_lev:pver)                                
       real(8) T_05(top_lev:pver)                               
       real(8) T_prime0(top_lev:pver)                           
       real(8) T_dprime(top_lev:pver)                           
       real(8) T_star(top_lev:pver)                             
    
       real(8) qv(top_lev:pver)                                 
                                                               
       real(8) qv1(top_lev:pver)                                
       real(8) qv_0(top_lev:pver)                               
       real(8) qv_05(top_lev:pver)                              
       real(8) qv_prime0(top_lev:pver)                          
       real(8) qv_dprime(top_lev:pver)                          
       real(8) qv_star(top_lev:pver)                            
    
       real(8) ql(top_lev:pver)                                 
                                                               
       real(8) ql1(top_lev:pver)                                
       real(8) ql_0(top_lev:pver)                               
       real(8) ql_05(top_lev:pver)                              
       real(8) ql_prime0(top_lev:pver)                          
       real(8) ql_dprime(top_lev:pver)                          
       real(8) ql_star(top_lev:pver)                            
    
       real(8) qi(top_lev:pver)                                 
                                                               
       real(8) qi1(top_lev:pver)                                
       real(8) qi_0(top_lev:pver)                               
       real(8) qi_05(top_lev:pver)                              
       real(8) qi_prime0(top_lev:pver)                          
       real(8) qi_dprime(top_lev:pver)                          
       real(8) qi_star(top_lev:pver)                            
    
       real(8) nl(top_lev:pver)                                 
                                                               
       real(8) nl1(top_lev:pver)                                
       real(8) nl_0(top_lev:pver)                               
       real(8) nl_05(top_lev:pver)                              
       real(8) nl_prime0(top_lev:pver)                          
       real(8) nl_dprime(top_lev:pver)                          
       real(8) nl_star(top_lev:pver)                            
    
       real(8) ni(top_lev:pver)                                 
                                                               
       real(8) ni1(top_lev:pver)                                
       real(8) ni_0(top_lev:pver)                               
       real(8) ni_05(top_lev:pver)                              
       real(8) ni_prime0(top_lev:pver)                          
       real(8) ni_dprime(top_lev:pver)                          
       real(8) ni_star(top_lev:pver)                            
    
       real(8) a_st(top_lev:pver)                               
       real(8) a_st_0(top_lev:pver)                             
       real(8) a_st_star(top_lev:pver)                          
    
       real(8) al_st(top_lev:pver)                              
       real(8) al_st_0(top_lev:pver)                            
       real(8) al_st_nc(top_lev:pver)                           
    
       real(8) ai_st(top_lev:pver)                              
       real(8) ai_st_0(top_lev:pver)                            
       real(8) ai_st_nc(top_lev:pver)                           
    
       real(8) ql_st(top_lev:pver)                              
       real(8) ql_st_0(top_lev:pver)                            
    
       real(8) qi_st(top_lev:pver)                              
       real(8) qi_st_0(top_lev:pver)                            
    
     ! Cumulus properties 
    
       real(8) dacudt(top_lev:pver)
       real(8) a_cu(top_lev:pver)
    
     ! Adjustment tendency in association with 'positive_moisture'
    
       real(8) Tten_pwi1(top_lev:pver)                          
       real(8) qvten_pwi1(top_lev:pver)                         
       real(8) qlten_pwi1(top_lev:pver)                         
       real(8) qiten_pwi1(top_lev:pver)                         
       real(8) nlten_pwi1(top_lev:pver)                         
       real(8) niten_pwi1(top_lev:pver)                         
    
       real(8) Tten_pwi2(top_lev:pver)                          
       real(8) qvten_pwi2(top_lev:pver)                         
       real(8) qlten_pwi2(top_lev:pver)                         
       real(8) qiten_pwi2(top_lev:pver)                         
       real(8) nlten_pwi2(top_lev:pver)                         
       real(8) niten_pwi2(top_lev:pver)                         
    
       real(8) A_T_adj(top_lev:pver)                            
       real(8) A_qv_adj(top_lev:pver)                           
       real(8) A_ql_adj(top_lev:pver)                           
       real(8) A_qi_adj(top_lev:pver)                           
       real(8) A_nl_adj(top_lev:pver)                           
       real(8) A_ni_adj(top_lev:pver)                           
    
     
    
       real(8) QQw1(top_lev:pver)           
       real(8) QQi1(top_lev:pver)           
       real(8) QQw2(top_lev:pver)           
       real(8) QQi2(top_lev:pver)           
    
       real(8) QQnl1(top_lev:pver)          
       real(8) QQni1(top_lev:pver)          
       real(8) QQnl2(top_lev:pver)          
       real(8) QQni2(top_lev:pver)          
    
     
    
       real(8) QQ(top_lev:pver)             
       real(8) QQw(top_lev:pver)            
       real(8) QQi(top_lev:pver)            
       real(8) QQnl(top_lev:pver)           
       real(8) QQni(top_lev:pver)           
       real(8) ACnl(top_lev:pver)           
       real(8) ACni(top_lev:pver)           
    
       real(8) QQw_prev(top_lev:pver)   
       real(8) QQi_prev(top_lev:pver)   
       real(8) QQnl_prev(top_lev:pver)  
       real(8) QQni_prev(top_lev:pver)  
    
       real(8) QQw_prog(top_lev:pver)   
       real(8) QQi_prog(top_lev:pver)   
       real(8) QQnl_prog(top_lev:pver)  
       real(8) QQni_prog(top_lev:pver)  
    
       real(8) QQ_final(top_lev:pver)                           
       real(8) QQw_final(top_lev:pver)                           
       real(8) QQi_final(top_lev:pver)                           
       real(8) QQn_final(top_lev:pver)                           
       real(8) QQnl_final(top_lev:pver)                          
       real(8) QQni_final(top_lev:pver)                          
    
       real(8) QQ_all(top_lev:pver)         
       real(8) QQw_all(top_lev:pver)        
       real(8) QQi_all(top_lev:pver)        
       real(8) QQn_all(top_lev:pver)        
       real(8) QQnl_all(top_lev:pver)       
       real(8) QQni_all(top_lev:pver)       
    
     
    
       real(8) U(top_lev:pver)                                  
       real(8) U_nc(top_lev:pver)                               
       real(8) G_nc(top_lev:pver)                               
       real(8) F_nc(top_lev:pver)                               
       real(8) alpha                                          
       real(8) beta                                           
       real(8) betast                                         
       real(8) gammal                                         
       real(8) gammai                                         
       real(8) gammaQ                                         
       real(8) deltal                                         
       real(8) deltai                                         
       real(8) A_Tc                                           
       real(8) A_qt                                           
       real(8) C_Tc                                           
       real(8) C_qt                                           
       real(8) dTcdt                                          
       real(8) dqtdt                                          
       real(8) dqtstldt                                       
       real(8) dqidt                                          
    
       real(8) dqlstdt                                        
       real(8) dalstdt                                        
       real(8) dastdt                                         
    
       real(8) anic                                           
       real(8) GG                                             
    
       real(8) aa(2,2)
       real(8) bb(2,1)
    
       real(8) zeros(top_lev:pver)
    
       real(8) qmin1(top_lev:pver)
       real(8) qmin2(top_lev:pver)
       real(8) qmin3(top_lev:pver)
    
       real(8) esat_a                             
       real(8) qsat_a(top_lev:pver)                        
       real(8) Twb_aw                             
       real(8) qvwb_aw(top_lev:pver)                       
    
       real(8) esat_b                                 
       real(8) qsat_b                                 
       real(8) dqsdT_b                                 
    
       logical  land
       real(8) tmp
    
       real(8), intent(out) :: d_rhmin_liq_PBL(top_lev:pver)
       real(8), intent(out) :: d_rhmin_ice_PBL(top_lev:pver)
       real(8), intent(out) :: d_rhmin_liq_det(top_lev:pver)
       real(8), intent(out) :: d_rhmin_ice_det(top_lev:pver)
       real(8) rhmini_arr(top_lev:pver)
       real(8) rhminl_arr(top_lev:pver)
       real(8) rhminl_adj_land_arr(top_lev:pver)
       real(8) rhminh_arr(top_lev:pver) 
       real(8) rhmin_liq_diag(top_lev:pver)
       real(8) rhmin_ice_diag(top_lev:pver)
    
       real(8) QQmax,QQmin,QQwmin,QQimin                      
       real(8) cone                                           
!xukai+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!i_rhmini, rhminl_const, rhminl_adj_land_const, rhminh_const, i_rhminl
        integer, intent(in) :: i_rhmini
        integer, intent(in) :: i_rhminl

        real(8), intent(in) :: rhmini_const
        real(8), intent(in) :: rhmaxi_const
        real(8), intent(in) :: rhminl_const
        real(8), intent(in) :: rhminl_adj_land_const     ! rhminl adjustment for snowfree land
        real(8), intent(in) :: rhminh_const              ! Critical RH for high-level liquid stratus clouds

        real(8), intent(in) :: qmins1
        real(8), intent(in) :: qmins2
        real(8), intent(in) :: qmins3

        real(8), intent(in) :: premib 
        real(8), intent(in) :: premit 
        real(8), intent(in) :: icecrit
        integer, intent(in) :: iceopt
        
        real(8) ::  qsmall 
        integer(8) :: pow_addr, exp_addr, log_addr, log10_addr
        real(8) :: estbl(8)
        integer :: niter
        integer :: core_oe
        logical CAMstfrac 
        real(8) :: qmin
        real(8) :: latice
        real(8) :: latvap
        real(8) :: cpair 
        !real(8), parameter  :: cc           = 0.1d0     ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )
        !real(8), parameter  :: ramda        = 0.5d0    ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
        real(8)  :: cc       ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )
        real(8)  :: ramda    ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
 
        !integer, target :: tmppind(10)
        integer :: adim, bdim
        call  math_agent_pow_c(pow_addr)
        call  math_agent_exp_c(exp_addr)
        call  math_agent_log_c(log_addr)
        call  math_agent_log10_c(log10_addr)
        call  math_agent_cos_c()
        call  math_agent_acos_c()
        !call  math_agent_sqrt_c()
        adim = 2
        bdim = 1
        qsmall = 1.d-18 
        niter = 2
        CAMstfrac = .false.
        qmin = qmins1
        latvap = wv_para(2) 
        latice = wv_para(3) 
        cpair = wv_para(5)
        !tmp1 => s_tendout(:)
        !tmp2 => qv_tendout(:)
        !tmpind => tmppind(:)
        !call init_tmp_addr(pver, nl_tendout(1), ni_tendout(1)) 
        cc           = 0.1d0  
        ramda        = 0.5d0 
!xukai--------------------------------------------------------------------------------------
       if(do_cldice_int == 1) then
           do_cldice = .true.
        else
           do_cldice = .false.
        endif
        !call lwpf_start_test_p1() 

       cone            = 0.999d0
       zeros(top_lev:pver)  = 0.d0
    
       ! ------------------------------------ !
       ! Global initialization of main output !
       ! ------------------------------------ !
    
         s_tendout(top_lev:pver)     = 0.d0
         qv_tendout(top_lev:pver)    = 0.d0
         ql_tendout(top_lev:pver)    = 0.d0
         qi_tendout(top_lev:pver)    = 0.d0
         nl_tendout(top_lev:pver)    = 0.d0
         ni_tendout(top_lev:pver)    = 0.d0
    
         qme(top_lev:pver)           = 0.d0
    
         cld(top_lev:pver)           = 0.d0
         al_st_star(top_lev:pver)    = 0.d0
         ai_st_star(top_lev:pver)    = 0.d0
         ql_st_star(top_lev:pver)    = 0.d0
         qi_st_star(top_lev:pver)    = 0.d0
    
       ! --------------------------------------- !
       ! Initialization of internal 2D variables !
       ! --------------------------------------- !
    
         T(top_lev:pver)             = 0.d0
         T1(top_lev:pver)            = 0.d0
         T_0(top_lev:pver)           = 0.d0
         T_05(top_lev:pver)          = 0.d0
         T_prime0(top_lev:pver)      = 0.d0
         T_dprime(top_lev:pver)      = 0.d0
         T_star(top_lev:pver)        = 0.d0
    
         qv(top_lev:pver)            = 0.d0
         qv1(top_lev:pver)           = 0.d0
         qv_0(top_lev:pver)          = 0.d0
         qv_05(top_lev:pver)         = 0.d0
         qv_prime0(top_lev:pver)     = 0.d0
         qv_dprime(top_lev:pver)     = 0.d0
         qv_star(top_lev:pver)       = 0.d0
    
         ql(top_lev:pver)            = 0.d0
         ql1(top_lev:pver)           = 0.d0
         ql_0(top_lev:pver)          = 0.d0
         ql_05(top_lev:pver)         = 0.d0
         ql_prime0(top_lev:pver)     = 0.d0
         ql_dprime(top_lev:pver)     = 0.d0
         ql_star(top_lev:pver)       = 0.d0
    
         qi(top_lev:pver)            = 0.d0
         qi1(top_lev:pver)           = 0.d0
         qi_0(top_lev:pver)          = 0.d0
         qi_05(top_lev:pver)         = 0.d0
         qi_prime0(top_lev:pver)     = 0.d0
         qi_dprime(top_lev:pver)     = 0.d0
         qi_star(top_lev:pver)       = 0.d0
    
         nl(top_lev:pver)            = 0.d0
         nl1(top_lev:pver)           = 0.d0
         nl_0(top_lev:pver)          = 0.d0
         nl_05(top_lev:pver)         = 0.d0
         nl_prime0(top_lev:pver)     = 0.d0
         nl_dprime(top_lev:pver)     = 0.d0
         nl_star(top_lev:pver)       = 0.d0
    
         ni(top_lev:pver)            = 0.d0
         ni1(top_lev:pver)           = 0.d0
         ni_0(top_lev:pver)          = 0.d0
         ni_05(top_lev:pver)         = 0.d0
         ni_prime0(top_lev:pver)     = 0.d0
         ni_dprime(top_lev:pver)     = 0.d0
         ni_star(top_lev:pver)       = 0.d0
    
         a_st(top_lev:pver)          = 0.d0
         a_st_0(top_lev:pver)        = 0.d0
         a_st_star(top_lev:pver)     = 0.d0
    
         al_st(top_lev:pver)         = 0.d0
         al_st_0(top_lev:pver)       = 0.d0
         al_st_nc(top_lev:pver)      = 0.d0
    
         ai_st(top_lev:pver)         = 0.d0
         ai_st_0(top_lev:pver)       = 0.d0
         ai_st_nc(top_lev:pver)      = 0.d0
    
         ql_st(top_lev:pver)         = 0.d0
         ql_st_0(top_lev:pver)       = 0.d0
    
         qi_st(top_lev:pver)         = 0.d0
         qi_st_0(top_lev:pver)       = 0.d0
    
     ! Cumulus properties 
    
         dacudt(top_lev:pver)        = 0.d0
         a_cu(top_lev:pver)          = 0.d0
    
     ! Adjustment tendency in association with 'positive_moisture'
    
         Tten_pwi1(top_lev:pver)     = 0.d0
         qvten_pwi1(top_lev:pver)    = 0.d0
         qlten_pwi1(top_lev:pver)    = 0.d0
         qiten_pwi1(top_lev:pver)    = 0.d0
         nlten_pwi1(top_lev:pver)    = 0.d0
         niten_pwi1(top_lev:pver)    = 0.d0
    
         Tten_pwi2(top_lev:pver)     = 0.d0
         qvten_pwi2(top_lev:pver)    = 0.d0
         qlten_pwi2(top_lev:pver)    = 0.d0
         qiten_pwi2(top_lev:pver)    = 0.d0
         nlten_pwi2(top_lev:pver)    = 0.d0
         niten_pwi2(top_lev:pver)    = 0.d0
    
         A_T_adj(top_lev:pver)       = 0.d0
         A_qv_adj(top_lev:pver)      = 0.d0
         A_ql_adj(top_lev:pver)      = 0.d0
         A_qi_adj(top_lev:pver)      = 0.d0
         A_nl_adj(top_lev:pver)      = 0.d0
         A_ni_adj(top_lev:pver)      = 0.d0
    
         qvadj   (top_lev:pver)      = 0.d0
         qladj   (top_lev:pver)      = 0.d0
         qiadj   (top_lev:pver)      = 0.d0
    
     ! Adjustment tendency in association with 'instratus_condensate'
    
         QQw1(top_lev:pver)          = 0.d0
         QQi1(top_lev:pver)          = 0.d0
         QQw2(top_lev:pver)          = 0.d0
         QQi2(top_lev:pver)          = 0.d0
    
         QQnl1(top_lev:pver)         = 0.d0
         QQni1(top_lev:pver)         = 0.d0
         QQnl2(top_lev:pver)         = 0.d0
         QQni2(top_lev:pver)         = 0.d0
    
         QQnl(top_lev:pver)          = 0.d0
         QQni(top_lev:pver)          = 0.d0
    
     ! Macrophysical process tendency variables
    
         QQ(top_lev:pver)            = 0.d0
         QQw(top_lev:pver)           = 0.d0
         QQi(top_lev:pver)           = 0.d0
         QQnl(top_lev:pver)          = 0.d0
         QQni(top_lev:pver)          = 0.d0
         ACnl(top_lev:pver)          = 0.d0
         ACni(top_lev:pver)          = 0.d0
    
         QQw_prev(top_lev:pver)      = 0.d0
         QQi_prev(top_lev:pver)      = 0.d0
         QQnl_prev(top_lev:pver)     = 0.d0
         QQni_prev(top_lev:pver)     = 0.d0
    
         QQw_prog(top_lev:pver)      = 0.d0
         QQi_prog(top_lev:pver)      = 0.d0
         QQnl_prog(top_lev:pver)     = 0.d0
         QQni_prog(top_lev:pver)     = 0.d0
    
         QQ_final(top_lev:pver)      = 0.d0                        
         QQw_final(top_lev:pver)     = 0.d0                  
         QQi_final(top_lev:pver)     = 0.d0           
         QQn_final(top_lev:pver)     = 0.d0    
         QQnl_final(top_lev:pver)    = 0.d0
         QQni_final(top_lev:pver)    = 0.d0
    
         QQ_all(top_lev:pver)        = 0.d0
         QQw_all(top_lev:pver)       = 0.d0
         QQi_all(top_lev:pver)       = 0.d0
         QQn_all(top_lev:pver)       = 0.d0
         QQnl_all(top_lev:pver)      = 0.d0
         QQni_all(top_lev:pver)      = 0.d0
    
     ! Coefficient for computing QQ and related processes
    
         U(top_lev:pver)             = 0.d0
         U_nc(top_lev:pver)          = 0.d0
         G_nc(top_lev:pver)          = 0.d0
         F_nc(top_lev:pver)          = 0.d0
    
     ! Other
    
         qmin1(top_lev:pver)         = 0.d0
         qmin2(top_lev:pver)         = 0.d0
         qmin3(top_lev:pver)         = 0.d0

        !!call swpf_start_test_p1() 
        !do i = 1, 100
        !    do k = top_lev, pver
        !        ql_tendout(k) = a_ql(k)    
        !    enddo
        !enddo
        !!call swpf_stop_test_p1() 
   
  
       s_tendout(top_lev:pver)  = a_t(top_lev:pver)              
       qv_tendout(top_lev:pver) = a_qv(top_lev:pver)           
       ql_tendout(top_lev:pver) = a_ql(top_lev:pver)    
       qi_tendout(top_lev:pver) = a_qi(top_lev:pver)           
       nl_tendout(top_lev:pver) = a_nl(top_lev:pver)      
       ni_tendout(top_lev:pver) = a_ni(top_lev:pver)      


        !call swpf_start_test_p2() 
      call rhcrit_calc( &
          ncol, dp, T0, p, &
          clrw_old, clri_old, tke, qtl_flx, &
          qti_flx, cmfr_det, qlr_det, qir_det, &
          rhmini_arr, rhminl_arr, rhminl_adj_land_arr, rhminh_arr, &
          d_rhmin_liq_PBL, d_rhmin_ice_PBL, d_rhmin_liq_det, d_rhmin_ice_det, &
          wv_para, pver, pverp, top_lev, &
          i_rhmini, i_rhminl, rhmini_const, rhmaxi_const, rhminl_const, rhminl_adj_land_const, rhminh_const) 

!#define DIFFRES
#ifdef DIFFRES
       s_tendout(top_lev:pver)  = rhmini_arr(top_lev:pver)              
       qv_tendout(top_lev:pver) = rhminl_arr(top_lev:pver)           
       ql_tendout(top_lev:pver) = rhminl_adj_land_arr(top_lev:pver)    
       qi_tendout(top_lev:pver) = rhminh_arr(top_lev:pver)           

#endif
!#undef DIFFRES

       !nl_tendout(top_lev:pver) = d_rhmin_liq_PBL(top_lev:pver)      
       !ni_tendout(top_lev:pver) = d_rhmin_ice_PBL(top_lev:pver)      
       !qme  (top_lev:pver)  =               d_rhmin_liq_det(top_lev:pver)
       !qvadj(top_lev:pver)  =               d_rhmin_ice_det(top_lev:pver)

        !call printf_int(i_rhmini)
        !call printf_int(i_rhminl)
        !call swpf_stop_test_p2() 
   ! ---------------------------------- !
   ! Compute cumulus-related properties ! 
   ! ---------------------------------- !

   dacudt(top_lev:pver) = &
        (a_cu0(top_lev:pver) - a_cud(top_lev:pver))/dt

   ! ---------------------------------------------------------------------- !
   ! set to zero for levels above
   ! ---------------------------------------------------------------------- !
   ql0(:top_lev-1) = 0.d0
   qi0(:top_lev-1) = 0.d0
   nl0(:top_lev-1) = 0.d0
   ni0(:top_lev-1) = 0.d0
   
   ! ---------------------------------------------------------------------- !
   ! Check if input non-cumulus pixels satisfie a non-negative constraint.  !
   ! If not, force all water vapor substances to be positive in all layers. !
   ! We should use 'old' cumulus properties for this routine.               !                
   ! ---------------------------------------------------------------------- !

   T1(:)    =  T0(:) 
   qv1(:)   = qv0(:) 
   ql1(:)   = ql0(:) 
   qi1(:)   = qi0(:) 
   nl1(:)   = nl0(:) 
   ni1(:)   = ni0(:) 

 
!do k = top_lev, pver
   qmin1(:) = qmins1
   qmin2(:) = qmins2
   qmin3(:) = qmins3
!enddo

!#define DIFFRES
#ifdef DIFFRES
        do k = top_lev, pver
        qme(k)      = qmins1  
        enddo
        !qme(:)      = qmin1(:)   
        qvadj(:)    = qmin2(:)
        qladj(:)    = qmin3(:)
 
#endif
!#undef  DIFFRES


        !call swpf_start_test_p3() 
   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp, & 
                           qv1, ql1, qi1, T1, qvten_pwi1, qlten_pwi1, &
                           qiten_pwi1, Tten_pwi1, do_cldice, &
                           wv_para, pver, pverp, top_lev)

        !call swpf_stop_test_p3() 

!#define DIFFRES
#ifdef DIFFRES
       s_tendout(:)  = qv1(:)              
       qv_tendout(:) = ql1(:)           
       ql_tendout(:) = qi1(:)    
       qi_tendout(:) = t1(:)           
       nl_tendout(:) = qvten_pwi1(:)      
       ni_tendout(:) = qlten_pwi1(:)      
       qme  (:)  =      qiten_pwi1(:)
       qvadj(:)  =      tten_pwi1(:)
#endif
!#undef DIFFRES

   do k = top_lev, pver
   !do i = 1, ncol
      if( ql1(k) .lt. qsmall ) then
          nlten_pwi1(k) = -nl1(k)/dt
          nl1(k)        = 0.d0
      endif 
      if( qi1(k) .lt. qsmall ) then
          niten_pwi1(k) = -ni1(k)/dt
          ni1(k)        = 0.d0
      endif 
   !enddo
   enddo
    !call set_wvflag(0)
        !call swpf_start_test_p4() 
   do k = top_lev, pver
      call instratus_condensate(ncol, k,                                   &
                                 p(k), T1(k), qv1(k), ql1(k), qi1(k),    &
                                 ni1(k),                                         &
                                 a_cud(k), zeros(k), zeros(k),               &
                                 zeros(k), zeros(k), zeros(k),               &
                                 landfrac, snowh,                                  &
                                 rhmini_arr(k), rhminl_arr(k), rhminl_adj_land_arr(k), rhminh_arr(k), &
                                 T_0(k), qv_0(k), ql_0(k), qi_0(k),        & 
                                 al_st_0(k), ai_st_0(k), ql_st_0(k), qi_st_0(k) ,&
                                 wv_para, rhmaxi_const, premib, premit, &
                                 icecrit, iceopt)
      a_st_0(k) = max(al_st_0(k),ai_st_0(k))
      QQw1(k)   = (ql_0(k) - ql1(k))/dt
      QQi1(k)   = (qi_0(k) - qi1(k))/dt
      ! -------------------------------------------------- !
      ! Reduce droplet concentration if evaporation occurs !
      ! Set a limit such that negative state not happens.  ! 
      ! -------------------------------------------------- !
      !do i = 1, ncol
         if( QQw1(k) .le. 0.d0 ) then
             if( ql1(k) .gt. qsmall ) then
                 QQnl1(k) = QQw1(k)*nl1(k)/ql1(k)
                 QQnl1(k) = min(0.d0,cone*max(QQnl1(k),-nl1(k)/dt))
             else
                 QQnl1(k) = 0.d0
             endif  
         endif 
         if( QQi1(k) .le. 0.d0 ) then
             if( qi1(k) .gt. qsmall ) then
                 QQni1(k) = QQi1(k)*ni1(k)/qi1(k)
                 QQni1(k) = min(0.d0,cone*max(QQni1(k),-ni1(k)/dt))
             else
                 QQni1(k) = 0.d0
             endif  
         endif 
      !enddo
   enddo
   nl_0(top_lev:) = max(0.d0,nl1(top_lev:)+QQnl1(top_lev:)*dt) 
   ni_0(top_lev:) = max(0.d0,ni1(top_lev:)+QQni1(top_lev:)*dt)

        !call swpf_stop_test_p4() 

    !call set_wvflag(0)
!#define DIFFRES
#ifdef DIFFRES
       s_tendout(:)  = t_0(:)              
       qv_tendout(:) = qv_0(:)           
       ql_tendout(:) = ql_0(:)    
       qi_tendout(:) = qi_0(:)           
       nl_tendout(:) = al_st_0(:)      
       ni_tendout(:) = ai_st_0(:)      
       qme  (:)  =      ql_st_0(:)
       qvadj(:)  =      qi_st_0(:)
#endif
!#undef DIFFRES


!+xukai - precision different much between MPE and CPE
   T_05(top_lev:)  =  T_0(top_lev:) + (  A_T(top_lev:) +  C_T(top_lev:) ) * dt
   qv_05(top_lev:) = qv_0(top_lev:) + ( A_qv(top_lev:) + C_qv(top_lev:) ) * dt
   ql_05(top_lev:) = ql_0(top_lev:) + ( A_ql(top_lev:) + C_ql(top_lev:) ) * dt
   qi_05(top_lev:) = qi_0(top_lev:) + ( A_qi(top_lev:) + C_qi(top_lev:) ) * dt 
   nl_05(top_lev:) = max(0.d0, nl_0(top_lev:) + ( A_nl(top_lev:) + C_nl(top_lev:) ) * dt )
   ni_05(top_lev:) = max(0.d0, ni_0(top_lev:) + ( A_ni(top_lev:) + C_ni(top_lev:) ) * dt )

!-xukai - precision different much between MPE and CPE
!#define VERRESULT 
#ifdef  VERRESULT 
       s_tendout(:)  =  t_05(:)              
       qv_tendout(:) =  qv_05(:)         
       ql_tendout(:) =  ql_05(:)  
       qi_tendout(:) =  qi_05(:)        
       nl_tendout(:) = nl_05(:)      
       ni_tendout(:) = ni_05(:)      

#endif
!#undef VERRESULT


        !call lwpf_start_test_p5() 
!   do k = top_lev, pver
!   ql_05(k) = ql_0(k) + ( A_ql(k) + C_ql(k) ) * dt
!   qi_05(k) = qi_0(k) + ( A_qi(k) + C_qi(k) ) * dt 
!   enddo
   !ql_05(top_lev:) =  ( A_ql(top_lev:) + C_ql(top_lev:) ) * dt
   !qi_05(top_lev:) =  ( A_qi(top_lev:) + C_qi(top_lev:) ) * dt 
!ql_05(top_lev:) =  ql_05(top_lev:)  + ql_0(top_lev:)  
!qi_05(top_lev:) =  qi_05(top_lev:)  + qi_0(top_lev:)  
   !ql_05(top_lev:) =   C_ql(top_lev:)
   !qi_05(top_lev:) =   C_qi(top_lev:) 

   !ql_05(top_lev:) =   A_ql(top_lev:) + C_ql(top_lev:) 
   !qi_05(top_lev:) =   A_qi(top_lev:) + C_qi(top_lev:) 

   !ql_05(top_lev:) =  ( A_ql(top_lev:) + C_ql(top_lev:) ) * dt
   !qi_05(top_lev:) =  ( A_qi(top_lev:) + C_qi(top_lev:) ) * dt 

!       s_tendout(:)  = ql_0(:)              
!       qv_tendout(:) = ( A_ql(top_lev:) + C_ql(top_lev:) ) * dt
!       ql_tendout(:) = qi_0(:)    
!       qi_tendout(:) = ( A_qi(top_lev:) + C_qi(top_lev:) ) * dt           
!
!
!       nl_tendout(:) = qv_05(:)      
!
!       ni_tendout(:) = ql_05(:)      
!       qme  (:)  =      qi_05(:)
!
!       qvadj(:)  =      t_05(:)
   !call set_wvflag(0)

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp, & 
                           qv_05, ql_05, qi_05, T_05, A_qv_adj, &
                           A_ql_adj, A_qi_adj, A_T_adj, do_cldice, &
                           wv_para, pver, pverp, top_lev)
                           
   !call set_wvflag(0)

!#define DIFFRES
#ifdef DIFFRES
        qme(:)      = qv_05(:)   
        qvadj(:)    = ql_05(:)
        qladj(:)    = qi_05(:)
        qiadj(:)    = t_05(:)
        qllim(:)    =  a_qv_adj(:)
        qilim(:)    = a_ql_adj(:)
#endif
!#undef  DIFFRES

!#define DIFFRES
#ifdef DIFFRES
        qllim(:)    =  a_qi_adj(:)
        qilim(:)    = a_t_adj(:)

#endif
!#undef  DIFFRES

        !call lwpf_stop_test_p5() 

   T(top_lev:)     = T_0(top_lev:)
   qv(top_lev:)    = qv_0(top_lev:)
   ql(top_lev:)    = ql_0(top_lev:)
   qi(top_lev:)    = qi_0(top_lev:)
   al_st(top_lev:) = al_st_0(top_lev:)
   ai_st(top_lev:) = ai_st_0(top_lev:)
   a_st(top_lev:)  = a_st_0(top_lev:)
   ql_st(top_lev:) = ql_st_0(top_lev:)
   qi_st(top_lev:) = qi_st_0(top_lev:)
   nl(top_lev:)    = nl_0(top_lev:)
   ni(top_lev:)    = ni_0(top_lev:)

   !     !call lwpf_start_test_p6() 
  do k = top_lev, pver
      call findsp_vc(qv_05(k), T_05(k), p(k), .false., &
                     Twb_aw, qvwb_aw(k), estbl, wv_para)

      call qsat_water(T_05(k), p(k), &
                      esat_a, qsat_a(k), wv_para)
   enddo
!#define DIFFRES
#ifdef DIFFRES
        qvadj(:)    = qvwb_aw(:)
        qiadj(:)    = qsat_a(:)
#endif
!#undef  DIFFRES
   !     !call lwpf_stop_test_p6() 
 
   


        !call lwpf_start_test_p7() 

   do iter = 1, niter

      ! ------------------------------------------ !
      ! Initialize array within the iteration loop !
      ! ------------------------------------------ !

      QQ(:)         = 0.d0
      QQw(:)        = 0.d0
      QQi(:)        = 0.d0
      QQnl(:)       = 0.d0
      QQni(:)       = 0.d0 
      QQw2(:)       = 0.d0
      QQi2(:)       = 0.d0
      QQnl2(:)      = 0.d0
      QQni2(:)      = 0.d0
      nlten_pwi2(:) = 0.d0
      niten_pwi2(:) = 0.d0
      ACnl(:)       = 0.d0
      ACni(:)       = 0.d0 
      aa(:,:)         = 0.d0
      bb(:,:)         = 0.d0

      do k = top_lev, pver

      call qsat_water(T(k), p(k), &
                      esat_b, qsat_b, wv_para, dqsdt=dqsdT_b)

      if( iter .eq. 1 ) then
          a_cu(k) = a_cud(k)
      else
          a_cu(k) = a_cu0(k)
      endif
      !do i = 1, ncol
         U(k)    =  qv(k)/qsat_b
         U_nc(k) =  U(k)
      !enddo
      if( CAMstfrac ) then
          call astG_RHU(U_nc(k),p(k),qv(k),landfrac,snowh,al_st_nc(k),G_nc(k),ncol,&
                        premib, premit, &
                        rhminl_arr(k), rhminl_adj_land_arr(k), rhminh_arr(k))                          
      else
          call astG_PDF(U_nc(k),p(k),qv(k),landfrac,snowh,al_st_nc(k),G_nc(k),ncol,&
                        premib, premit, &
                        rhminl_arr(k), rhminl_adj_land_arr(k), rhminh_arr(k))
      endif
      call aist_vector(qv(k),T(k),p(k),qi(k),ni(k),landfrac,snowh,ai_st_nc(k),ncol,&
                       wv_para, premib, premit, icecrit, iceopt, & 
                       rhmaxi_const, rhmini_arr(k), rhminl_arr(k), rhminl_adj_land_arr(k), rhminh_arr(k))

      ai_st(k)  =  (1.d0 -a_cu(k))*ai_st_nc(k)
      al_st(k)  =  (1.d0 -a_cu(k))*al_st_nc(k)
      a_st(k)   =  max(al_st(k),ai_st(k))  

      !do i = 1, ncol

         ! -------------------------------------------------------- !
         ! Compute basic thermodynamic coefficients for computing Q !
         ! -------------------------------------------------------- !

         alpha  =  1.d0/qsat_b
         beta   =  dqsdT_b*(qv(k)/(qsat_b*qsat_b))
         betast =  alpha*dqsdT_b 
         gammal =  alpha + (latvap/cpair)*beta
         gammai =  alpha + ((latvap+latice)/cpair)*beta
         gammaQ =  alpha + (latvap/cpair)*beta
         deltal =  1.d0 + a_st(k)*(latvap/cpair)*(betast/alpha)
         deltai =  1.d0 + a_st(k)*((latvap+latice)/cpair)*(betast/alpha)
         A_Tc   =  A_T(k)+A_T_adj(k)-(latvap/cpair)*(A_ql(k)+A_ql_adj(k))-((latvap+latice)/cpair)*(A_qi(k)+A_qi_adj(k))
         A_qt   =  A_qv(k) + A_qv_adj(k) + A_ql(k) + A_ql_adj(k) + A_qi(k) + A_qi_adj(k)
         C_Tc   =  C_T(k) - (latvap/cpair)*C_ql(k) - ((latvap+latice)/cpair)*C_qi(k)
         C_qt   =  C_qv(k) + C_ql(k) + C_qi(k)
         dTcdt  =  A_Tc + C_Tc
         dqtdt  =  A_qt + C_qt
       ! dqtstldt = A_qt + C_ql(k)/max(1.e-2_r8,al_st(k))                             ! Original  
       ! dqtstldt = A_qt - A_qi(k) - A_qi_adj(k) + C_ql(k)/max(1.e-2_r8,al_st(k)) ! New 1 on Dec.30.2009.
         dqtstldt = A_qt - A_qi(k) - A_qi_adj(k) + C_qlst(k)                        ! New 2 on Dec.30.2009.
       ! dqtstldt = A_qt + C_qt                                                           ! Original Conservative treatment
       ! dqtstldt = A_qt - A_qi(k) - A_qi_adj(k) + C_qt - C_qi(k)            ! New Conservative treatment on Dec.30.2009
         dqidt = A_qi(k) + A_qi_adj(k) + C_qi(k) 

         anic    = max(1.d-8,(1.d0-a_cu(k)))
         GG      = G_nc(k)/anic
         aa(1,1) = gammal*al_st(k)
         aa(1,2) = GG + gammal*cc*ql_st(k)          
         aa(2,1) = alpha + (latvap/cpair)*betast*al_st(k)
         aa(2,2) = (latvap/cpair)*betast*cc*ql_st(k) 
         bb(1,1) = alpha*dqtdt - beta*dTcdt - gammai*dqidt - GG*al_st_nc(k)*dacudt(k) + F_nc(k) 
         bb(2,1) = alpha*dqtstldt - betast*(dTcdt + ((latvap+latice)/cpair)*dqidt) 

       !s_tendout(k)  = aa(1,1)              
       !qv_tendout(k) = aa(1,2)           
       !ql_tendout(k) = aa(2,1)    
       !qi_tendout(k) = aa(2,2)           

       !call gaussj(aa(1:2,1:2),2,2,bb(1:2,1),1,1, wv_para)
       call gaussj(aa(1:2,1:2),adim,adim,bb(1:2,1),bdim,bdim, wv_para)

       !qladj(k) = aa(1,1)
       !qiadj(k) = aa(1,2)
       !qllim(k) = aa(2,1)
       !qilim(k) = aa(2,2)
       !nl_tendout(k) = bb(1,1)      
       !ni_tendout(k) = bb(2,1)      


         dqlstdt = bb(1,1)
         dalstdt = bb(2,1)
         QQ(k) = al_st(k)*dqlstdt + cc*ql_st(k)*dalstdt - ( A_ql(k) + A_ql_adj(k) + C_ql(k) )

       ! ------------------------------------------------------------ !
       ! Limiter for QQ                                               !
       ! Here, 'fice' should be from the reference equilibrium state  !
       ! since QQ itself is computed from the reference state.        !
       ! From the assumption used for derivation of QQ(i), it must be !
       ! that QQw(i) = QQ(i)*(1._r8-fice(i)), QQi(i) = QQ(i)*fice(i)  !  
       ! ------------------------------------------------------------ !

         if( QQ(k) .ge. 0.d0 ) then
             QQmax    = (qv_05(k) - qmins1)/dt ! For ghost cumulus & semi-ghost ice stratus
             QQmax    = max(0.d0,QQmax) 
             QQ(k)  = min(QQ(k),QQmax)
             QQw(k) = QQ(k)
             QQi(k) = 0.d0 
         else
             QQmin  = 0.d0
             if( qv_05(k) .lt. qsat_a(k) ) QQmin = min(0.d0,cone*(qv_05(k)-qvwb_aw(k))/dt)
             QQ(k)  = max(QQ(k),QQmin)
             QQw(k) = QQ(k)
             QQi(k) = 0.d0
             QQwmin   = min(0.d0,-cone*ql_05(k)/dt)
             QQimin   = min(0.d0,-cone*qi_05(k)/dt)
             QQw(k) = min(0.d0,max(QQw(k),QQwmin))
             QQi(k) = min(0.d0,max(QQi(k),QQimin))
         endif

       ! -------------------------------------------------- !
       ! Reduce droplet concentration if evaporation occurs !
       ! Note 'QQnl1,QQni1' are computed from the reference !
       ! equilibrium state but limiter is from 'nl_05'.     !
       ! -------------------------------------------------- !

         if( QQw(k) .lt. 0.d0 ) then
             if( ql_05(k) .gt. qsmall ) then
                 QQnl(k) = QQw(k)*nl_05(k)/ql_05(k)
                 QQnl(k) = min(0.d0,cone*max(QQnl(k),-nl_05(k)/dt))
             else
                 QQnl(k) = 0.d0
             endif  
         endif 

         if( QQi(k) .lt. 0.d0 ) then
             if( qi_05(k) .gt. qsmall ) then
                 QQni(k) = QQi(k)*ni_05(k)/qi_05(k)
                 QQni(k) = min(0.d0,cone*max(QQni(k),-ni_05(k)/dt))
             else
                 QQni(k) = 0.d0
             endif  
         endif 

      !enddo
      enddo
   !if(iter == 1) then
   !    qladj(top_lev:)      = ai_st(top_lev:)      
   !    qiadj(top_lev:)      = al_st(top_lev:)
   !    qllim(top_lev:)      = a_st(top_lev:)
   !    qilim(top_lev:)      = qq(top_lev:)
   !    qme(top_lev:) = ql_05(top_lev:)       
   !    qvadj(top_lev:) = qqw(top_lev:)       
   !endif


    ! -------------------------------------------------------------------- !
    ! Until now, we have finished computing all necessary tendencies       ! 
    ! from the equilibrium input state (T_0).                              !
    ! If ramda = 0 : fully explicit scheme                                 !
    !    ramda = 1 : fully implicit scheme                                 !
    ! Note that 'ramda = 0.5 with niter = 2' can mimic                     !
    ! -------------------------------------------------------------------- !

      if( iter .eq. 1 ) then
          QQw_prev(top_lev:)  = QQw(top_lev:)       
          QQi_prev(top_lev:)  = QQi(top_lev:)   
          QQnl_prev(top_lev:) = QQnl(top_lev:)       
          QQni_prev(top_lev:) = QQni(top_lev:)   
      endif

      QQw_prog(top_lev:)   = ramda*QQw(top_lev:)   + (1.d0-ramda)*QQw_prev(top_lev:)
      QQi_prog(top_lev:)   = ramda*QQi(top_lev:)   + (1.d0-ramda)*QQi_prev(top_lev:)
      QQnl_prog(top_lev:)  = ramda*QQnl(top_lev:)  + (1.d0-ramda)*QQnl_prev(top_lev:)
      QQni_prog(top_lev:)  = ramda*QQni(top_lev:)  + (1.d0-ramda)*QQni_prev(top_lev:)

      QQw_prev(top_lev:)   = QQw_prog(top_lev:)
      QQi_prev(top_lev:)   = QQi_prog(top_lev:)
      QQnl_prev(top_lev:)  = QQnl_prog(top_lev:)
      QQni_prev(top_lev:)  = QQni_prog(top_lev:)

!#define DIFFRES
#ifdef DIFFRES 
   if(iter == 1) then
       qladj(top_lev:)      = T_star(top_lev:)      
       qiadj(top_lev:)      = qv_star(top_lev:)
       qllim(top_lev:)      = ql_star(top_lev:)
       qilim(top_lev:)      = qi_star(top_lev:)
       qme(top_lev:) = al_st_star(top_lev:)       
       qvadj(top_lev:) = ai_st_star(top_lev:)       


       !qiadj(top_lev:)      = T_dprime(top_lev:)
       !qllim(top_lev:)      =  qv_dprime(top_lev:)
       !qilim(top_lev:)      = ql_dprime(top_lev:)
       !qvadj(top_lev:)      = qi_dprime(top_lev:)       
       !qme(top_lev:)        =  ni_dprime(top_lev:)       

       !qiadj(top_lev:)      = qqw_prev(top_lev:)
       !qllim(top_lev:)      =  qqi_prev(top_lev:)
       !qilim(top_lev:)      =  qqnl_prev(top_lev:)
       !qvadj(top_lev:)      = QQni_prev(top_lev:)       
 
       !qiadj(top_lev:)      = A_nl(top_lev:)
       !qllim(top_lev:)      = C_nl(top_lev:)
       !qilim(top_lev:)      = A_ql_adj(top_lev:)
       !qvadj(top_lev:)      = A_qi_adj(top_lev:)       
       !qme(top_lev:)        = C_qv(top_lev:)       
 
   endif
#endif
!#undef DIFFRES


    ! -------------------------------------------------------- !
    ! Compute final prognostic state on which final diagnostic !
    ! in-stratus condensate adjustment is applied in the below.!
    ! Important : I must check whether there are any external  !  
    !             advective forcings of 'A_nl(i,k),A_ni(i,k)'. !
    !             Even they are (i.e., advection of aerosol),  !
    !             actual droplet activation will be performd   !
    !             in microphysics, so it will be completely    !
    !             reasonable to 'A_nl(i,k)=A_ni(i,k)=0'.       !
    ! -------------------------------------------------------- !
!+xukai - precision different much between MPE and CPE
    do k = top_lev, pver
    !do i = 1, ncol
       T_prime0(k)  = T_0(k)  + dt*( A_T(k)  +  A_T_adj(k) +  C_T(k) + &
            (latvap*QQw_prog(k)+(latvap+latice)*QQi_prog(k))/cpair )
       qv_prime0(k) = qv_0(k) + dt*( A_qv(k) + A_qv_adj(k) + C_qv(k) - QQw_prog(k) - QQi_prog(k) )
       ql_prime0(k) = ql_0(k) + dt*( A_ql(k) + A_ql_adj(k) + C_ql(k) + QQw_prog(k) )
       qi_prime0(k) = qi_0(k) + dt*( A_qi(k) + A_qi_adj(k) + C_qi(k) + QQi_prog(k) )
       nl_prime0(k) = max(0.d0,nl_0(k) + dt*( A_nl(k) + C_nl(k) + QQnl_prog(k) ))
       ni_prime0(k) = max(0.d0,ni_0(k) + dt*( A_ni(k) + C_ni(k) + QQni_prog(k) ))
       if( ql_prime0(k) .lt. qsmall ) nl_prime0(k) = 0.d0
       if( qi_prime0(k) .lt. qsmall ) ni_prime0(k) = 0.d0
    !enddo
    enddo
!-xukai - precision different much between MPE and CPE

   ! -------------------------------------------------- !
   ! Perform diagnostic 'positive_moisture' constraint. !
   ! -------------------------------------------------- !

   T_dprime(top_lev:)  =  T_prime0(top_lev:) 
   qv_dprime(top_lev:) = qv_prime0(top_lev:) 
   ql_dprime(top_lev:) = ql_prime0(top_lev:) 
   qi_dprime(top_lev:) = qi_prime0(top_lev:) 
   nl_dprime(top_lev:) = nl_prime0(top_lev:) 
   ni_dprime(top_lev:) = ni_prime0(top_lev:) 


   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp,          & 
                           qv_dprime, ql_dprime, qi_dprime, T_dprime,  &
                           qvten_pwi2, qlten_pwi2, qiten_pwi2, Tten_pwi2, do_cldice, &
                           wv_para, pver, pverp, top_lev)

   do k = top_lev, pver
   !do i = 1, ncol
      if( ql_dprime(k) .lt. qsmall ) then
          nlten_pwi2(k) = -nl_dprime(k)/dt
          nl_dprime(k)   = 0.d0
      endif 
      if( qi_dprime(k) .lt. qsmall ) then
          niten_pwi2(k) = -ni_dprime(k)/dt
          ni_dprime(k)   = 0.d0
      endif 
   !enddo
   enddo


   ! -------------------------------------------------------------- !
   ! Add tendency associated with detrainment of cumulus condensate !
   ! This tendency is not used in computing Q                       !
   ! Since D_ql,D_qi,D_nl,D_ni > 0, don't need to worry about       !
   ! negative scalar.                                               !
   ! This tendency is not reflected into Fzs2, which is OK.         !
   ! -------------------------------------------------------------- !

   T_dprime(top_lev:)   =  T_dprime(top_lev:)  + D_T(top_lev:) * dt 
   qv_dprime(top_lev:)  = qv_dprime(top_lev:) + D_qv(top_lev:) * dt 
   ql_dprime(top_lev:)  = ql_dprime(top_lev:) + D_ql(top_lev:) * dt
   qi_dprime(top_lev:)  = qi_dprime(top_lev:) + D_qi(top_lev:) * dt
   nl_dprime(top_lev:)  = nl_dprime(top_lev:) + D_nl(top_lev:) * dt 
   ni_dprime(top_lev:)  = ni_dprime(top_lev:) + D_ni(top_lev:) * dt

    !---------------------------------------------------------- !
    !Impose diagnostic upper and lower limits on the in-stratus !
    !condensate amount. This produces a final equilibrium state !
    !at the end of each iterative process.                      !
    !---------------------------------------------------------- !

   do k = top_lev, pver
      call instratus_condensate( ncol           , k              , p(k)        , &
                                 T_dprime(k)  , qv_dprime(k) , ql_dprime(k) , qi_dprime(k), &
                                 ni_dprime(k) ,                                                   &
                                 a_cu0(k)     , zeros(k)     , zeros(k)     ,                 & 
                                 zeros(k)     , zeros(k)     , zeros(k)     ,                 &
                                 landfrac       , snowh          ,                                  &
                                 rhmini_arr(k), rhminl_arr(k), rhminl_adj_land_arr(k), rhminh_arr(k), &
                                 T_star(k)    , qv_star(k)   , ql_star(k)   , qi_star(k)  , & 
                                 al_st_star(k), ai_st_star(k), ql_st_star(k), qi_st_star(k) , &
                                 wv_para, rhmaxi_const, premib, premit, &
                                 icecrit, iceopt)
      a_st_star(k)  = max(al_st_star(k),ai_st_star(k))
      QQw2(k) = (ql_star(k) - ql_dprime(k))/dt
      QQi2(k) = (qi_star(k) - qi_dprime(k))/dt
      ! -------------------------------------------------- !
      ! Reduce droplet concentration if evaporation occurs !
      ! -------------------------------------------------- !
      !do i = 1, ncol
         if( QQw2(k) .le. 0.d0 ) then
             if( ql_dprime(k) .ge. qsmall ) then
                 QQnl2(k) = QQw2(k)*nl_dprime(k)/ql_dprime(k)
                 QQnl2(k) = min(0.d0,cone*max(QQnl2(k),-nl_dprime(k)/dt))
             else
                 QQnl2(k) = 0.d0
             endif  
         endif 
         if( QQi2(k) .le. 0.d0 ) then
             if( qi_dprime(k) .gt. qsmall ) then
                 QQni2(k) = QQi2(k)*ni_dprime(k)/qi_dprime(k)
                 QQni2(k) = min(0.d0,cone*max(QQni2(k),-ni_dprime(k)/dt))
             else
                 QQni2(k) = 0.d0
             endif  
         endif 
      !enddo
   enddo
   nl_star(top_lev:) = max(0.d0,nl_dprime(top_lev:)+QQnl2(top_lev:)*dt) 
   ni_star(top_lev:) = max(0.d0,ni_dprime(top_lev:)+QQni2(top_lev:)*dt)

   !if(iter == 3) then
   !    qladj(top_lev:)      = T_star(top_lev:)      
   !    qiadj(top_lev:)      = qv_star(top_lev:)
   !    qllim(top_lev:)      = ql_star(top_lev:)
   !    qilim(top_lev:)      = qi_star(top_lev:)
   !    qme(top_lev:) = al_st_star(top_lev:)       
   !    qvadj(top_lev:) = ai_st_star(top_lev:)       


   !    !qiadj(top_lev:)      = a_st_star(top_lev:)
   !    !qllim(top_lev:)      = ql_st_star(top_lev:)
   !    !qilim(top_lev:)      = qi_st_star(top_lev:)
   !    !qvadj(top_lev:)      = nl_star(top_lev:)       
   !    !qme(top_lev:)        = ni_star(top_lev:)       
   !endif



   ! ------------------------------------------ !
   ! Final adjustment of droplet concentration. !
   ! Set # to zero if there is no cloud.        !
   ! ------------------------------------------ !

   do k = top_lev, pver
   !do i = 1, ncol 
      if( ql_star(k) .lt. qsmall ) then
          ACnl(k) = - nl_star(k)/dt
          nl_star(k) = 0.d0
      endif
      if( qi_star(k) .lt. qsmall ) then
          ACni(k) = - ni_star(k)/dt
          ni_star(k) = 0.d0
      endif
   !enddo
   enddo

   ! ----------------------------------------------------- !
   ! Define equilibrium reference state for next iteration !
   ! ----------------------------------------------------- !    


   T(top_lev:)     = T_star(top_lev:)
   qv(top_lev:)    = qv_star(top_lev:)
   ql(top_lev:)    = ql_star(top_lev:)
   qi(top_lev:)    = qi_star(top_lev:)
   al_st(top_lev:) = al_st_star(top_lev:)
   ai_st(top_lev:) = ai_st_star(top_lev:)
   a_st(top_lev:)  = a_st_star(top_lev:)
   ql_st(top_lev:) = ql_st_star(top_lev:)
   qi_st(top_lev:) = qi_st_star(top_lev:)
   nl(top_lev:)    = nl_star(top_lev:)
   ni(top_lev:)    = ni_star(top_lev:)
!#define DIFFRES
#ifdef DIFFRES 
   if(iter == 2) then
       qladj(top_lev:)      = T_star(top_lev:)      
       qiadj(top_lev:)      = qv_star(top_lev:)
       qllim(top_lev:)      = ql_star(top_lev:)
       qilim(top_lev:)      = qi_star(top_lev:)
       qme(top_lev:) = al_st_star(top_lev:)       
       qvadj(top_lev:) = ai_st_star(top_lev:)       


       !qiadj(top_lev:)      = a_st_star(top_lev:)
       !qllim(top_lev:)      = ql_st_star(top_lev:)
       !qilim(top_lev:)      = qi_st_star(top_lev:)
       !qvadj(top_lev:)      = nl_star(top_lev:)       
       !qme(top_lev:)        = ni_star(top_lev:)       
   endif
#endif
!#undef DIFFRES


   enddo ! End of 'iter' prognostic iterative computation

        !call lwpf_stop_test_p7() 

   ! ------------------------------------------------------------------------ !
   ! Compute final tendencies of main output variables and diagnostic outputs !
   ! Note that the very input state [T0,qv0,ql0,qi0] are                      !
   ! marched to [T_star,qv_star,ql_star,qi_star] with equilibrium             !
   ! stratus informations of [a_st_star,ql_st_star,qi_st_star] by             !
   ! below final tendencies and [A_T,A_qv,A_ql,A_qi]                          !
   ! ------------------------------------------------------------------------ !

   !t0(:) = nl_star(:) 
   !t0(top_lev:) = nl0(top_lev:) 
   !t0(top_lev:) = ( nl_star(top_lev:) - nl0(top_lev:) )/dt 
   !!qv0(:) = nl0(:)  
   !qv0(top_lev:) = (A_nl(top_lev:)+C_nl(top_lev:))


   !ql0(top_lev:) = ( nl_star(top_lev:) - nl0(top_lev:) )/dt - &
   !     (A_nl(top_lev:)+C_nl(top_lev:))

   !ql0(:) = a_nl(:)
   !qi0(:) = c_nl(:)

   ! ------------------ !
   ! Process tendencies !
   ! ------------------ !
    !call lwpf_start_test_p8() 

   QQw_final(top_lev:)  = QQw_prog(top_lev:)
   QQi_final(top_lev:)  = QQi_prog(top_lev:)
   QQ_final(top_lev:)   = QQw_final(top_lev:) + QQi_final(top_lev:)
   QQw_all(top_lev:)    = QQw_prog(top_lev:)  + QQw1(top_lev:) + QQw2(top_lev:) + &
        qlten_pwi1(top_lev:) + qlten_pwi2(top_lev:) + A_ql_adj(top_lev:)
   QQi_all(top_lev:)    = QQi_prog(top_lev:)  + QQi1(top_lev:) + QQi2(top_lev:) + &
        qiten_pwi1(top_lev:) + qiten_pwi2(top_lev:) + A_qi_adj(top_lev:)
   QQ_all(top_lev:)     = QQw_all(top_lev:)   + QQi_all(top_lev:)
   QQnl_final(top_lev:) = QQnl_prog(top_lev:)
   QQni_final(top_lev:) = QQni_prog(top_lev:)
   QQn_final(top_lev:)  = QQnl_final(top_lev:) + QQni_final(top_lev:)
   QQnl_all(top_lev:)   = QQnl_prog(top_lev:)  + QQnl1(top_lev:) + QQnl2(top_lev:) + &
        nlten_pwi1(top_lev:) + nlten_pwi2(top_lev:) + ACnl(top_lev:) + A_nl_adj(top_lev:)
   QQni_all(top_lev:)   = QQni_prog(top_lev:)  + QQni1(top_lev:) + QQni2(top_lev:) + &
        niten_pwi1(top_lev:) + niten_pwi2(top_lev:) + ACni(top_lev:) + A_ni_adj(top_lev:)
   QQn_all(top_lev:)    = QQnl_all(top_lev:)   + QQni_all(top_lev:)

   qme(top_lev:)        = QQ_final(top_lev:)   
   qvadj(top_lev:)      = qvten_pwi1(top_lev:) + qvten_pwi2(top_lev:) + A_qv_adj(top_lev:)
   qladj(top_lev:)      = qlten_pwi1(top_lev:) + qlten_pwi2(top_lev:) + A_ql_adj(top_lev:)
   qiadj(top_lev:)      = qiten_pwi1(top_lev:) + qiten_pwi2(top_lev:) + A_qi_adj(top_lev:)
   qllim(top_lev:)      = QQw1      (top_lev:) + QQw2      (top_lev:)
   qilim(top_lev:)      = QQi1      (top_lev:) + QQi2      (top_lev:)

   ! ----------------- !
   ! Output tendencies !
   ! ----------------- !

   s_tendout(top_lev:)  = cpair*( T_star(top_lev:)  -  T0(top_lev:) )/dt - &
        cpair*(A_T(top_lev:)+C_T(top_lev:))
   qv_tendout(top_lev:) =    ( qv_star(top_lev:) - qv0(top_lev:) )/dt - &
        (A_qv(top_lev:)+C_qv(top_lev:))
   ql_tendout(top_lev:) =    ( ql_star(top_lev:) - ql0(top_lev:) )/dt - &
        (A_ql(top_lev:)+C_ql(top_lev:))
   qi_tendout(top_lev:) =    ( qi_star(top_lev:) - qi0(top_lev:) )/dt - &
        (A_qi(top_lev:)+C_qi(top_lev:))
   nl_tendout(top_lev:) =    ( nl_star(top_lev:) - nl0(top_lev:) )/dt - &
        (A_nl(top_lev:)+C_nl(top_lev:))
   ni_tendout(top_lev:) =    ( ni_star(top_lev:) - ni0(top_lev:) )/dt - &
        (A_ni(top_lev:)+C_ni(top_lev:))

   if (.not. do_cldice) then
      do k = top_lev, pver
         do i = 1, ncol

            ! Don't want either qi or ni tendencies, but the code above is somewhat convoluted and
            ! is trying to adjust both (small numbers). Just force it to zero here.
            qi_tendout(k) = 0.d0
            ni_tendout(k) = 0.d0
          end do
      end do
   end if

   ! ------------------ !
   ! Net cloud fraction !
   ! ------------------ !

   cld(top_lev:) = a_st_star(top_lev:) + a_cu0(top_lev:)

   ! --------------------------------- !
   ! Updated grid-mean state variables !
   ! --------------------------------- !

   T0(top_lev:)  = T_star(top_lev:)
   qv0(top_lev:) = qv_star(top_lev:)
   ql0(top_lev:) = ql_star(top_lev:)
   qi0(top_lev:) = qi_star(top_lev:)
   nl0(top_lev:) = nl_star(top_lev:)
   ni0(top_lev:) = ni_star(top_lev:)

    !call lwpf_stop_test_p8() 
#ifdef CPERES
#endif

#ifdef XKOPT
#endif
        end subroutine
