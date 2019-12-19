module ecosys_share

! !MODULE: ecosys_share

!-----------------------------------------------------------------------------
!   This module contains definitions of variables, derived types, and
!   functions/subroutines that are used in ecocys_mod.F90 as well as by
!   other modules that make use of the ecosys_mod.
!
!   The variables are shared using threading with pointers, and need to be
!   pointed to in the code.
!   Note: So far the values of all of these fields are set in ecosys_mod
!   and are NOT modified in the other modules
!   A. Jahn, NCAR
!-----------------------------------------------------------------------------
  USE kinds_mod
  USE blocks, ONLY: nx_block, ny_block
  USE domain_size, ONLY: max_blocks_clinic
  USE constants, only: c0

  implicit none

  public
  save

!-----------------------------------------------------------------------------

  logical (log_kind) :: &
     lmarginal_seas               ! Is ecosystem active in marginal seas ?

  character(char_len) :: &
     ecosys_tadvect_ctype         ! advection method for ecosys tracers

!-----------------------------------------------------------------------------
! number of ecosystem constituents and grazing interactions
!-----------------------------------------------------------------------------
  INTEGER (KIND=int_kind), PARAMETER :: &
       zooplankton_cnt = ZOOPLANKTON_CNT, &
       autotroph_cnt   = AUTOTROPH_CNT,   &
       grazer_prey_cnt = GRAZER_PREY_CNT

!-----------------------------------------------------------------------------
!   derived type for grazers
!-----------------------------------------------------------------------------

 TYPE, PUBLIC :: zooplankton_type
     CHARACTER(char_len) :: sname, lname
     INTEGER (KIND=int_kind) :: &
          C_ind                    ! tracer indices for zooplankton carbon
     REAL(KIND=r8) :: &
          z_mort_0,         & ! zoo linear mort rate (1/sec)
          z_mort2_0,        & ! zoo quad mort rate (1/sec/((mmol C/m3))
          loss_thres      !zoo conc. where losses go to zero
  END TYPE

  TYPE(zooplankton_type), DIMENSION(zooplankton_cnt) :: zooplankton

!-----------------------------------------------------------------------------
!   derived type for functional group
!-----------------------------------------------------------------------------

  TYPE, PUBLIC :: autotroph_type
     CHARACTER(char_len) :: sname, lname
     LOGICAL(KIND=log_kind) :: &
        Nfixer,                             & ! flag set to true if this autotroph fixes N2
        imp_calcifier,                      & ! flag set to true if this autotroph implicitly handles calcification
        exp_calcifier                         ! flag set to true if this autotroph explicitly handles calcification
     INTEGER (KIND=int_kind) :: &
        Chl_ind, C_ind, Fe_ind,             & ! tracer indices for Chl, C, Fe content
        Si_ind, CaCO3_ind,                  & ! tracer indices for Si, CaCO3 content
        C13_ind, C14_ind,                   & ! tracer indices for 13C, 14C
        Ca13CO3_ind, Ca14CO3_ind              ! tracer indices for 13CaCO3, 14CaCO3
     REAL(KIND=r8) :: &
        kFe, kPO4, kDOP, kNO3, kNH4, kSiO3, & ! nutrient uptake half-sat constants
        Qp,                                 & ! P/C ratio
        gQfe_0, gQfe_min,                   & ! initial and minimum fe/C ratio
        alphaPI,                            & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
        PCref,                              & ! max C-spec. grth rate at tref (1/sec)
        thetaN_max,                         & ! max thetaN (Chl/N) (mg Chl/mmol N)
        loss_thres, loss_thres2,            & ! conc. where losses go to zero
        temp_thres,                         & ! Temp. where concentration threshold and photosynth. rate drops
        mort, mort2,                        & ! linear and quadratic mortality rates (1/sec), (1/sec/((mmol C/m3))
        agg_rate_max, agg_rate_min,         & ! max and min agg. rate (1/d)
        loss_poc                              ! routing of loss term

  END TYPE

  INTEGER (KIND=int_kind), PARAMETER :: &
     sp_ind          = 1, &  ! small phytoplankton
     diat_ind        = 2, &  ! diatoms
     diaz_ind        = 3     ! diazotrophs

  TYPE(autotroph_type), DIMENSION(autotroph_cnt) :: autotrophs


!-----------------------------------------------------------------------------
!   derived type for grazing
!-----------------------------------------------------------------------------

 TYPE, PUBLIC :: grazing_type
    CHARACTER(char_len) :: sname, lname
     INTEGER (KIND=int_kind) :: &
          grazing_function,           & ! functional form of grazing parameterization
          auto_ind_cnt,               & ! number of autotrophs in prey-clase auto_ind
          zoo_ind_cnt                   ! number of zooplankton in prey-clase zoo_ind
     REAL(KIND=r8) :: &
          z_umax_0,                           & ! max zoo growth rate at tref (1/sec)
          z_grz,                              & ! grazing coef. (mmol C/m^3)^2
          graze_zoo, graze_poc, graze_doc,    & ! routing of grazed term, remainder goes to dic
          f_zoo_detr                            ! fraction of zoo losses to detrital
     INTEGER (KIND=int_kind), DIMENSION(autotroph_cnt) :: &
          auto_ind
     INTEGER (KIND=int_kind), DIMENSION(zooplankton_cnt) :: &
          zoo_ind
  END TYPE

  TYPE(grazing_type), DIMENSION(grazer_prey_cnt,zooplankton_cnt) :: grazing

!-----------------------------------------------------------------------------
!  derived type for implicit handling of sinking particulate matter
!-----------------------------------------------------------------------------

  TYPE, PUBLIC :: sinking_particle
      REAL(KIND=r8)  :: &
         diss,        & ! dissolution length for soft subclass
         gamma,       & ! fraction of production -> hard subclass
         mass,        & ! mass of 1e9 base units in g
         rho            ! QA mass ratio of POC to this particle class

      REAL(KIND=r8) , DIMENSION(nx_block,ny_block,max_blocks_clinic) :: &
         sflux_in,    & ! incoming flux of soft subclass (base units/cm^2/sec)
         hflux_in,    & ! incoming flux of hard subclass (base units/cm^2/sec)
         prod,        & ! production term (base units/cm^3/sec)
         sflux_out,   & ! outgoing flux of soft subclass (base units/cm^2/sec)
         hflux_out,   & ! outgoing flux of hard subclass (base units/cm^2/sec)
         sed_loss,    & ! loss to sediments (base units/cm^s/sec)
         remin          ! remineralization term (base units/cm^3/sec)
    END TYPE

!*****************************************************************************
! Variable definitions
!*****************************************************************************


   type(sinking_particle) :: &
      POC,            & ! base units = nmol C
      P_CaCO3           ! base units = nmol CaCO3


  real (r8), dimension(nx_block,ny_block,max_blocks_clinic), target, public :: &
      DIC_SURF_fields,         & ! surface values of DIC for solver
      CO2STAR_SURF_fields,     & ! CO2STAR from solver
      DCO2STAR_SURF_fields,    & ! DCO2STAR from solver
      PV_SURF_fields,          & ! piston velocity (cm/s)
      CO3_fields,              & ! carbonate ion
      CO3_SURF_fields,         & ! Surface carbonate ion
      HCO3_fields,             & ! bicarbonate ion
      H2CO3_fields,            & ! carbonic acid
      DIC_loc_fields,          & ! local copy of model DIC
      DOC_loc_fields,          & ! local copy of model DOC
      O2_loc_fields,           & ! local copy of model O2
      NO3_loc_fields,          & ! local copy of model NO3
      decay_CaCO3_fields,      & ! scaling factor for dissolution of CaCO3
      DECAY_Hard_fields,       & ! scaling factor for dissolution of Hard Ballast
      decay_POC_E_fields,      & ! scaling factor for dissolution of excess POC
      poc_diss_fields,         & ! diss. length used (cm)
      caco3_diss_fields,       & ! caco3 diss. length used (cm)
      POC_PROD_avail_fields,   & ! POC production available for excess POC flux
      DOC_remin_fields,        & ! remineralization of 13C DOC (mmol C/m^3/sec)
      P_CaCO3_sflux_out_fields,& ! P_CaCO3 sflux_out from ecosys before getting set to zero for k=KMT
      P_CaCO3_hflux_out_fields,& ! P_CaCO3_hflux_out from ecosys before getting set to zero for k=KMT
      POC_sflux_out_fields,    & ! POC_sflux_out from ecosys before getting set to zero for k=KMT
      POC_hflux_out_fields,    & ! POC_hflux_out from ecosys before getting set to zero for k=KMT
      P_CaCO3_remin_fields,    & ! P_CaCO3 remin from ecosys before it gets modified for k=KMT
      POC_remin_fields,        & ! POC remin from ecosys before it gets modified for k=KMT
      dic_riv_flux_fields,     & ! River input of DIC in ecosystem (from file)
      doc_riv_flux_fields        ! River input of DOC in ecosystem (from file)

    real (r8), dimension(nx_block,ny_block,zooplankton_cnt,max_blocks_clinic), target, public :: &
      zooC_loc_fields,         & ! local copy of model zooC
      zoo_loss_fields,         & ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
      zoo_loss_poc_fields,     & ! zoo_loss routed to large detrital (mmol C/m^3/sec)
      zoo_loss_doc_fields,     & ! zoo_loss routed to doc (mmol C/m^3/sec)
      zoo_loss_dic_fields        ! zoo_loss routed to dic (mmol C/m^3/sec)

    real (r8), dimension(nx_block,ny_block,autotroph_cnt,max_blocks_clinic), target, public :: &
      CaCO3_PROD_fields,        & ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
      QCaCO3_fields,            & ! small phyto CaCO3/C ratio (mmol CaCO3/mmol C)
      autotrophCaCO3_loc_fields,& ! local copy of model autotroph CaCO3
      autotrophChl_loc_fields,  & ! local copy of model autotroph Chl
      autotrophC_loc_fields,    & ! local copy of model autotroph C
      autotrophFe_loc_fields,   & ! local copy of model autotroph Fe
      autotrophSi_loc_fields,   & ! local copy of model autotroph Si
      auto_graze_fields,        & ! autotroph grazing rate (mmol C/m^3/sec)
      auto_graze_zoo_fields,    & ! auto_graze routed to zoo (mmol C/m^3/sec)
      auto_graze_poc_fields,    & ! auto_graze routed to poc (mmol C/m^3/sec)
      auto_graze_doc_fields,    & ! auto_graze routed to doc (mmol C/m^3/sec)
      auto_graze_dic_fields,    & ! auto_graze routed to dic (mmol C/m^3/sec)
      auto_loss_fields,         & ! autotroph non-grazing mort (mmol C/m^3/sec)
      auto_loss_poc_fields,     & ! auto_loss routed to poc (mmol C/m^3/sec)
      auto_loss_doc_fields,     & ! auto_loss routed to doc (mmol C/m^3/sec)
      auto_loss_dic_fields,     & ! auto_loss routed to dic (mmol C/m^3/sec)
      auto_agg_fields,          & ! autotroph aggregation (mmol C/m^3/sec)
      photoC_fields,            & ! C-fixation (mmol C/m^3/sec)
      PCphoto_fields              ! C-specific rate of photosynth. (1/sec)

!***********************************************************************

 contains

!*****************************************************************************
! Functions and subroutines used by more than one ecosystem-related module
!*****************************************************************************

!*****************************************************************************
!BOP
! !IROUTINE: SCHMIDT_CO2
! !INTERFACE:

 function SCHMIDT_CO2(SST, LAND_MASK)

! !DESCRIPTION:
!  Compute Schmidt number of CO2 in seawater as function of SST
!  where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!  ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
!  pp. 7373-7382, May 15, 1992
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: SST

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: SCHMIDT_CO2

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a = 2073.1_r8, &
      b = 125.62_r8, &
      c = 3.6276_r8, &
      d = 0.043219_r8

   where (LAND_MASK)
      SCHMIDT_CO2 = a + SST * (-b + SST * (c + SST * (-d)))
   elsewhere
      SCHMIDT_CO2 = c0
   end where

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_CO2

!*****************************************************************************


 end module ecosys_share
