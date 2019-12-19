
module ionosphere

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to compute the relevant physics needed to improve and add a realistic
! simulation of the ionosphere.  Main routines are initialization, ion/electron
! temperature calculation, ambipolar diffusion
!
! Authors: Joe McInerney/Hanli Liu/Art Richmond
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,   only : r8 => shr_kind_r8            ! Real kind to declare variables
  use ppgrid,         only : pcols, pver, pverp           ! Dimensions and chunk bounds
  use cam_history,    only : outfld                       ! Routine to output fields to history files
  use physics_types,  only : physics_state, &             ! Structures containing physics state variables
                             physics_ptend, &             ! Structures containing physics tendency variables
                             physics_ptend_init           ! Routine to initialize physics tendency variables
  use physics_buffer, only : pbuf_add_field, &            ! 
                             pbuf_get_index,dtype_r8, &   !
			     physics_buffer_desc, &       !
			     pbuf_get_field, &            ! Needed to access physics buffer
			     pbuf_set_field, &            !
			     pbuf_get_chunk               !
  use cam_logfile,    only : iulog                        ! Output unit for run.out file
  use mo_jeuv,        only : nIonRates                    ! Number of ionization rates in mo_photo
  use shr_const_mod,  only : kboltz => shr_const_boltz, &
                             pi => shr_const_pi           ! Boltzmann constant and pi
  use time_manager,   only : is_first_step,get_nstep
  use chem_mods,      only : adv_mass                     ! Array holding mass values for short lived species
  use cam_abortutils, only : endrun
  use mo_chem_utls,   only : get_spc_ndx                  ! Routine to get index of adv_mass array for short lived species

  implicit none

  save
  
  private   ! Make default type private to the module

  !------------------------
  ! PUBLIC: interfaces 
  !------------------------
  public ionos_init     ! Initialization
  public ionos_register ! Registration of ionosphere variables in pbuf physics buffer
  public ionos_tend     ! Calculate tendencies for extended model ionosphere
  public tridag         ! Generic tridiagonal solver routine

  !------------------------------------------------------------------------
  ! PRIVATE: Rest of the data and interfaces are private to this module
  !------------------------------------------------------------------------   
  real(r8), parameter               :: kboltz_ev = 8.617E-5_r8 ! Boltzmann constant (eV/K)
  real(r8), parameter               :: temax = 7.0E3_r8        ! maximum electron temperature (K)
  real(r8), parameter               :: dayOPFlux = 2.0E8_r8    ! Daytime O+ flux at upper boundary (
  real(r8), parameter               :: nightOPFlux = -2.0E8_r8 ! Nighttime O+ flux at upper boundary (

  !------------------------------------------------------------------------------ 
  ! magnetic grid dimensions (assumed resolution of 2deg)
  !------------------------------------------------------------------------------ 
  integer, parameter                :: nmlat = 90	       ! mlat
  integer, parameter                :: nmlon = 180             ! mlon 
                
  integer                           :: teTiBot                  ! bottom of ionosphere calculations
  integer                           :: teTiBotP                 ! bottom of ionosphere calculations

  integer                           :: vAmbBot                   ! bottom of O+ transport calculations
  integer                           :: vAmbBotP                  ! bottom of O+ transport calculations

  real(r8)                          :: calDay                  ! current calendar day
  real(r8)                          :: rads2Degs               ! radians to degrees
  real(r8)                          :: degs2Rads               ! degrees to radians

  real(r8)                          :: f107                    ! 10.7 cm solar flux

  type ionos_state

    real(r8), dimension(pcols)      :: cosZenAngR              ! cosine of zenith angle (radians)
    real(r8), dimension(pcols)      :: zenAngD                 ! zenith angle (degrees)

    real(r8), dimension(pcols,pver) :: bNorth3d  ! northward component of magnetic field units?
    real(r8), dimension(pcols,pver) :: bEast3d   ! eastward component of magnetic field
    real(r8), dimension(pcols,pver) :: bDown3d   ! downward component of magnetic field

    real(r8), dimension(pcols,pver,nIonRates) :: ionPRates    ! ionization rates temporary array (s-1 cm-3)
    real(r8), dimension(pcols,pver)           :: sumIonPRates ! Sum of ionization rates for O+,O2+,N+,N2+,NO+ (s-2 cm-3)

    real(r8), dimension(pcols,pver)  :: ue2d  ! horizontal x(eastward) component of ExB drift with added vertical dimension
    real(r8), dimension(pcols,pver)  :: ve2d  ! horizontal y(northward) component of ExB drift with added vertical dimension
    real(r8), dimension(pcols,pver)  :: we2d  ! vertical z component of ExB drift with added vertical dimension - midpoints

    real(r8), dimension(pcols,pverp) :: wei2d    ! vertical z component of ExB drift with added vertical dimension - interfaces

    real(r8), dimension(pcols,pverp) :: omegai   ! vertical velocity on interface levels (Pa/s)
 
    real(r8), dimension(pcols,pver)  :: dipMag   ! dip angle for each column (radians)
    real(r8), dimension(pcols,pver)  :: dipMagD  ! dip angle for each column (degrees)

    real(r8), dimension(pcols,pverp) :: tNInt    ! Interface Temperature (K)

    real(r8), dimension(pcols,pver)  :: ndensN2  ! N2 number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO2  ! O2 number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO1  ! O number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensNO  ! NO number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensN1  ! N number density  (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensE   ! E electron number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensOp  ! O plus number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO2p ! O2 plus ion number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensNOp ! NO plus ion number density  (cm-3)

    real(r8), dimension(pcols,pver)  :: sourceg4 ! g4 source term for electron/ion temperature update

    real(r8), dimension(pcols,pverp) :: rairvi   ! Constituent dependent gas constant on interface levels

  end type ionos_state

contains

!==============================================================================

  subroutine ionos_init()
  
!-----------------------------------------------------------------------
! Time independent initialization for ionosphere simulation.
!-----------------------------------------------------------------------

    use cam_history,      only : addfld, add_default, phys_decomp ! Routines and variables for adding fields to history output
    use phys_control,     only : phys_getopts !Method used to get flag for waccmx ionosphere output variables     

    logical :: history_waccmx

    call phys_getopts(history_waccmx_out=history_waccmx)

    !-------------------------------------------------------------------------------
    !  Add history variables for ionosphere 
    !-------------------------------------------------------------------------------
    call addfld ('TE'           , 'K',pver, 'I','Electron Temperature', phys_decomp)
    call addfld ('TI'           , 'K',pver, 'I','Ion Temperature', phys_decomp)
    call addfld ('QIN'          , 'J/kg/s',pver, 'I','JOULE+IN Heating', phys_decomp)
    call addfld ('LOSS_g3'     , ' ',pver         , 'I','Loss Term g3',  phys_decomp)
    call addfld ('LOSS_EI'      , ' ',pver         , 'I','Loss Term EI',  phys_decomp)
    call addfld ('LOSS_IN'      , ' ',pver         , 'I','Loss Term IN',  phys_decomp)
    call addfld ('IONI_RATE'    , ' ',pver         , 'I','Total ionization',  phys_decomp)
    call addfld ('IONI_Eff'     , ' ',pver         , 'I','ionization efficiency',  phys_decomp)
    call addfld ('SOURCE_g4'    , ' ',pver         , 'I','SOURCE g4',  phys_decomp)
    call addfld ('AUR_IRATESUM' , ' ',pver         , 'I','Auroral ionization',  phys_decomp)

    call addfld ('OpI' , ' ',pver         , 'I','O+ Ionosphere',  phys_decomp)
    call addfld ('eI' , ' ',pver         , 'I','e Ionosphere',  phys_decomp)

    !-------------------------------------------------------------------------------
    !  Set default values for ionosphere history variables
    !-------------------------------------------------------------------------------
    if (history_waccmx) then
       call add_default ('TE'	 , 1, ' ')
       call add_default ('TI'	 , 1, ' ')
       call add_default ('QIN'   , 1, ' ')
       call add_default ('LOSS_g3'	, 1, ' ')
       call add_default ('LOSS_EI'	 , 1, ' ')
       call add_default ('LOSS_IN'	 , 1, ' ')
       call add_default ('IONI_RATE'	 , 1, ' ')
       call add_default ('IONI_Eff'	 , 1, ' ')
       call add_default ('SOURCE_g4'	 , 1, ' ')
       call add_default ('AUR_IRATESUM'  , 1, ' ')

       call add_default ('OpI'  , 1, ' ')
       call add_default ('eI'  , 1, ' ')    
    end if

  end subroutine ionos_init

!==============================================================================     

  subroutine ionos_register

    !-----------------------------------------------------------------------
    ! Register ionosphere variables with physics buffer:
    !
    ! Ion production rates pcols,pver,nIonRates,
    !   so firstdim = 1 middledim = pver lastdim = nIonRates.
    ! 
    ! pcols dimension and lchnk assumed here
    !
    !-----------------------------------------------------------------------

    integer :: idx
  
    !------------------------------------------------------------------------------
    ! Electron temperature in physics buffer (global so can write to history files) 
    !------------------------------------------------------------------------------
    call pbuf_add_field('ElecTemp','global',dtype_r8,(/pcols,pver/),idx)
    
    !--------------------------------------------------------------------------
    ! Ion temperature in physics buffer (global so can write to history files)
    !--------------------------------------------------------------------------
    call pbuf_add_field('IonTemp', 'global',dtype_r8,(/pcols,pver/),idx)

  end subroutine ionos_register

!==============================================================================

  subroutine ionos_tend(state, ptend, pbuf, ztodt)

    !-------------------------------------------------------------------------------------
    ! Calculate dry static energy and O+ tendency for extended ionosphere simulation 
    !-------------------------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    use time_manager,        only : get_step_size
    use short_lived_species, only : slvd_index,slvd_pbf_ndx => pbf_idx ! Routines to access short lived species in storage and pbuf
    use constituents,        only : pcnst, cnst_get_ind                ! Number of and routine to get indices for constituents
    use physics_types,       only : physics_ptend_sum
    
    type(physics_state), intent(in)    :: state               ! physics state structure
    type(physics_ptend), intent(inout) :: ptend               ! parameterization tendency structure
    type(physics_buffer_desc),pointer  :: pbuf(:)             ! physics buffer

    real(r8),            intent(in)    :: ztodt               ! Physics time step

!---------------------------Local storage-------------------------------

    type(physics_ptend)                :: ptend_loc           ! Local parameterization tendencies
    type(ionos_state)                  :: istate              ! ionosphere state structure

    integer :: lchnk      ! Chunk number 
    integer :: ncol       ! Number of columns in chunk 
                
    integer :: teTiBot			          ! bottom of ionosphere calculations
    integer :: vAmbBot			          ! bottom of O+ transport calculations

    integer :: indxTe                                   ! pbuf index for electron temperature
    integer :: indxTi                                   ! pbuf index for ion temperature

    integer :: sIndxOp                   ! pbuf index for O+ mixing ratio
    integer :: indxOp                    ! state%q index for O+ mixing ratio

    real(r8), dimension(:,:), pointer   :: tE           ! Pointer to electron temperature in pbuf (K) 
    real(r8), dimension(:,:), pointer   :: ti           ! Pointer to ion temperature in pbuf (K) 

    real(r8), dimension(:,:), pointer   :: mmrPOp       ! Pointer to access O+ in pbuf

    real(r8), dimension(pcols,pver)     :: opTend       ! O+ tendency
 
    logical :: ls         
    logical :: lq(pcnst)         

    !----------------------------------------------------------------
    !  Get number of this chunk
    !----------------------------------------------------------------
     lchnk = state%lchnk
     ncol = state%ncol
    
    !------------------------------------------------------------
    !  Initialize data needed in the ionosphere calculations
    !------------------------------------------------------------
    call ionos_timestep_init(state, pbuf, istate, teTiBot, vAmbBot)

    !-------------------------------------------------------------------------------------------------------------------
    !  Get electron temperature from physics buffer. 
    !-------------------------------------------------------------------------------------------------------------------
    indxTe = pbuf_get_index( 'ElecTemp' )
    call pbuf_get_field(pbuf, indxTe, tE)

    indxTi = pbuf_get_index( 'IonTemp' )
    call pbuf_get_field(pbuf, indxTi, tI)

   !------------------------------------------------------------------------------------------
   !  Attempt to get indices to determine if O+ is advected in state and short-lived in pbuf  
   !------------------------------------------------------------------------------------------
    call cnst_get_ind( 'Op',  indxOp, abort=.false. )
    sIndxOp = slvd_index( 'Op' )

    !--------------------------------------------------------------------
    ! Initialize local physics_ptend object or get Op from physics buffer
    !--------------------------------------------------------------------
    if (indxOp > 0) then
      ls = .TRUE.   
      lq(:) = .FALSE.
      lq(indxOp) = .TRUE. 
      call physics_ptend_init(ptend_loc, state%psetcols, 'ionosphere', ls=ls, lq=lq)
    else if (sIndxOp > 0) then 
      ls = .TRUE.   
      call physics_ptend_init(ptend_loc, state%psetcols, 'ionosphere', ls=ls)
      call pbuf_get_field(pbuf, slvd_pbf_ndx, mmrPOp, start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
    else 
      call endrun('ionosphere: Cannot find state or pbuf index for Op in ionos_tend')             
    endif    

    !-----------------------------------------------------------------
    !  Get electron temperature and update dry static energy tendency
    !-----------------------------------------------------------------
    call ionos_teti(state, ptend%s, ptend_loc%s, ztodt, istate, tE, tI, teTiBot)
       
    !------------------------------------
    !  Get vertical O+ transport tendency
    !------------------------------------
    call ionos_vamb(state, opTend, ztodt, pbuf, istate, tE, tI, mmrPOp, vAmbBot)

    !----------------------------------------------------------------------------
    !  Add O+ tendency to the output tendency
    ! ---------------------------------------------------------------------------- 
    if (indxOp > 0) then
      ptend_loc%q(:,:,indxOp) = opTend(:,:)
    endif

    call physics_ptend_sum(ptend_loc, ptend, ncol)

    !--------------------------------------------------------------
    !  Make Te and Ti fields available for output to history files
    !--------------------------------------------------------------
    call outfld ('TE'            , tE           , pcols, lchnk)
    call outfld ('TI'            , ti           , pcols, lchnk)
       
    return

  end subroutine ionos_tend

!===============================================================================

  subroutine ionos_timestep_init(state, pbuf, istate, teTiBot, vAmbBot)
  
    !---------------------------------------------------------------------------------------
    ! Time independent initialization for extended ionosphere simulation called in phys_init
    ! of physpkg module which is called in cam_comp module
    !---------------------------------------------------------------------------------------
    use dycore,           only : get_resolution
    use interpolate_data, only : lininterp
    use constituents,     only : cnst_get_ind, cnst_mw            ! Routines to get molecular weights for constituents
    use mo_apex,          only : bnorth, beast, bdown             ! Magnetic field components
    use time_manager,     only : get_curr_calday                  ! Routine to get current calendar day
    use mo_solar_parms,   only : solar_parms_get                  ! Routine to get solar parameters, i.e. f107
    use physconst,        only : rairv, mbarv                     ! Constituent dependent rair and mbar
    use ref_pres,         only : press_lim_idx                    
    
    use short_lived_species, only : slvd_index,slvd_pbf_ndx => pbf_idx ! Routines to access short lived species     

    type(physics_buffer_desc), pointer  :: pbuf(:)             ! physics buffer
    type(physics_state), intent(in),    target :: state        ! physics state structure
    type(ionos_state),   intent(inout), target :: istate       ! ionosphere state structure

    integer, intent(out) :: teTiBot  ! bottom of ionosphere calculations 
    integer, intent(out) :: vAmbBot   ! bottom of ionosphere calculations 

!---------------------------Local storage-------------------------------    
    integer,parameter :: nCnst = 9     ! Number of species needed from state%q or pbuf

    integer :: lchnk      ! Chunk number 
    integer :: ncol       ! Number of columns in current chunk 

    integer :: indxUe     ! pbuf index for eastward ExB drift
    integer :: indxVe     ! pbuf index for northward ExB drift
    integer :: indxWe     ! pbuf index for upward ExB drift

    integer :: indxIR     ! pbuf index for ionization rates
    integer :: indxAIPRS  ! pbuf index for aurora ion production rate sum

    integer :: indxCnst   ! Constituent index used in cslculating densities

    integer :: indxSLvd   ! index of pbuf to access short lived species
    integer :: sIndx      ! index of adv_mass for any short lived species to access constituent mass

    integer :: iVer       ! Counter for vertical loops
    integer :: iCol       ! Counter for column loops
    integer :: iIonR      ! Counter for ionization rates loops
    integer :: iCnst      ! Counter for constituent loop
   
    integer :: indxSP     ! pbuf index for Pedersen Conductivity
    integer :: indxSH     ! pbuf index for Hall Conductivity

    real(r8), parameter :: teTiBotPres   = 50._r8   ! Pressure above which electron/ion temperature are calculated in WACCM-X. (Pa)
    real(r8), parameter :: vAmbBotPres   = 0.001_r8 ! Pressure above which O+ is calculated in WACCM-X. (Pa)
    
    real(r8), dimension(:), pointer   :: ue  ! pointer to eastward ExB drift in pbuf (from module iondrag)
    real(r8), dimension(:), pointer   :: ve  ! pointer to northward ExB drift in pbuf
    real(r8), dimension(:), pointer   :: we  ! pointer to upward ExB drift in pbuf

    character(len = 3), dimension(nCnst) :: cCnst

    real(r8), dimension(:,:), pointer :: sigma_ped    ! Pointer to Pedersen Conductivity in pbuf (siemens/m) from module iondrag
    real(r8), dimension(:,:), pointer :: sigma_hall   ! Pointer to Hall Conductivity in pbuf (siemens/m)

    real(r8), dimension(:,:), pointer :: mmrP         ! Pointer to access short lived species in pbuf

    real(r8), dimension(:),pointer    :: geoLatR  ! Latitude (radians)  Make ncol because zenith aurora are ncol
    real(r8), dimension(:),pointer    :: geoLonR  ! Longitude (radians)

    real(r8), dimension(:,:),pointer  :: pMid     ! Midpoint pressure (Pa)
    real(r8), dimension(:,:),pointer  :: tN       ! Neutral temperature (K)

    real(r8), dimension(:,:),pointer  :: tNInt   ! Interface Temperture (K)

    real(r8), dimension(:),pointer    :: cosZenAngR            ! cosine of zenith angle (radians)
    real(r8), dimension(:),pointer    :: zenAngD               ! zenith angle (degrees)

    real(r8), dimension(:,:),pointer  :: bNorth3d  ! northward component of magnetic field units?
    real(r8), dimension(:,:),pointer  :: bEast3d   ! eastward component of magnetic field
    real(r8), dimension(:,:),pointer  :: bDown3d          ! downward component of magnetic field

    real(r8), dimension(:,:),pointer  :: ue2d  ! horizontal x(eastward) component of ExB drift with added vertical dimension
    real(r8), dimension(:,:),pointer  :: ve2d  ! horizontal y(northward) component of ExB drift with added vertical dimension
    real(r8), dimension(:,:),pointer  :: we2d  ! vertical z component of ExB drift with added vertical dimension - midpoints

    real(r8), dimension(pcols,pver)   :: omega         ! vertical velocity at midpoint levels (Pa/s)
    real(r8), dimension(pcols,pver)   :: sourceR       ! R term of source g4 calculation
    real(r8), dimension(pcols,pver)   :: sourceEff     ! Efficiency term of source g4 calculation

    real(r8), dimension(:,:),pointer  :: wei2d   ! vertical z component of ExB drift with added vertical dimension - interfaces

    real(r8), dimension(:,:),pointer  :: omegai  ! vertical velocity on interface levels (Pa/s)
    
    real(r8), dimension(:,:),pointer  :: rairvi        ! Constituent dependent gas constant
 
    real(r8), dimension(:,:),pointer  :: dipMag  ! dip angle for each column (radians)
    real(r8), dimension(:,:),pointer  :: dipMagD ! dip angle for each column (degrees)

    real(r8) :: rMassN2    ! N2 molecular weight kg/kmol

    real(r8) :: rMass      ! Constituent molecular weight kg/kmol

    real(r8), dimension(pcols,pver)   :: mmrN2    ! N2 mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver)   :: mmrO2    ! O2 mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver)   :: mmrO1    ! O mass mixing ratio kg/kg

    real(r8), dimension(pcols,pver)   :: mmr      ! Constituent mass mixing ratio kg/kg

    real(r8), dimension(:,:),pointer  :: ndensN2  ! N2 number density (cm-3)
    real(r8), dimension(pcols,pver)   :: ndensO2  ! O2 number density (cm-3)
    real(r8), dimension(pcols,pver)   :: ndensO1  ! O number density (cm-3)
    real(r8), dimension(pcols,pver)   :: ndensN1  ! N number density  (cm-3)
    real(r8), dimension(pcols,pver)   :: ndensE   ! E electron number density (cm-3)

    real(r8), dimension(:,:)  ,pointer :: ndens    ! Constituent number density  (cm-3)

    real(r8), dimension(:,:,:),pointer :: ionRates     ! Pointer to ionization rates for O+,O2+,N+,N2+,NO+ in pbuf (s-1) 
                                                       !                               (from modules mo_jeuv and mo_jshort)

    real(r8), dimension(:,:,:),pointer :: ionPRates    ! ionization rates temporary array (s-1 cm-3)
    real(r8), dimension(:,:)  ,pointer :: sumIonPRates ! Sum of ionization rates for O+,O2+,N+,N2+,NO+ (s-1 cm-3)

    real(r8), dimension(:,:)  ,pointer :: aurIPRateSum ! Auroral ion production sum for O2+,O+,N2+ (s-1 cm-3 from module mo_aurora)

    real(r8), dimension(:,:)  ,pointer :: sourceg4     ! g4 source term for electron/ion temperature update

    logical, dimension(pcols)          :: do_aurora    ! Logical to test if aurora ion production rates calculated in mo_aurora 
                                                       ! (set to zero)

!--------------------------------------------------------------------------------
     
    omega             = 0._r8
    sourceR           = 0._r8  
    sourceEff         = 0._r8
   
    mmrN2             = 0._r8
    mmrO2             = 0._r8
    mmrO1             = 0._r8
    mmr               = 0._r8

    ndensO2(:,:)      = 0._r8
    ndensO1(:,:)      = 0._r8
    ndensN1(:,:)      = 0._r8
    ndensE(:,:)       = 0._r8

    sourceR(:,:)      = 0._r8
    sourceEff(:,:)    = 0._r8
    
    !--------------------------------------------------------------------------------------
    !  Get lchnk from state 
    !--------------------------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

    !------------------------------------------------------------------------------------------------------
    !  Set the bottom of the ionosphere calculations at around 50 Pascals or 0.5 hectopascals(millibars).  
    !  teTiBotPres is in Pascals.
    !------------------------------------------------------------------------------------------------------
    teTiBot  = press_lim_idx(teTiBotPres, top=.false.)       

    !------------------------------------------------------------------------------------------------------
    !  Set the bottom of the O+ transport calculations at around .001 Pascals or 
    !  .00001 hectopascals(millibars) or ~140km.  
    !------------------------------------------------------------------------------------------------------
    vAmbBot  = press_lim_idx(vAmbBotPres, top=.false.)
    
    !--------------------------------------------------------------
    !  Set radians to degrees variable and get zenith angle
    !--------------------------------------------------------------
    rads2Degs   = 180._r8/pi
    degs2Rads   = pi/180._r8

    !----------------------------------------------------------------
    !  Get latitude and longitude of each column in this chunk
    !----------------------------------------------------------------    
    geoLatR => state%lat(1:ncol)
    geoLonR => state%lon(1:ncol)

    !-------------------------------------------------------------------------------------------------------
    !  Need to get midpoint and interface pressure and neutral temperature from state structure (pcols,pver)
    !-------------------------------------------------------------------------------------------------------
    pMid => state%pmid(1:ncol,1:pver)
    tN   => state%t(1:ncol,1:pver)

    tNInt => istate%tNInt(1:ncol,1:pverp)
    cosZenAngR    => istate%cosZenAngR(1:ncol)
    zenAngD       => istate%zenAngD(1:ncol)
    
    bNorth3d => istate%bNorth3d(1:ncol,1:pver)
    bEast3d  => istate%bEast3d(1:ncol,1:pver)
    bDown3d  => istate%bDown3d(1:ncol,1:pver)
 
    ue2d => istate%ue2d(1:ncol,1:pver)
    ve2d => istate%ve2d(1:ncol,1:pver)
    we2d => istate%we2d(1:ncol,1:pver)
    
    wei2d => istate%wei2d(1:ncol,1:pverp)

    omegai => istate%omegai(1:ncol,1:pverp)

    rairvi => istate%rairvi(1:ncol,1:pverp)
 
    dipMag  => istate%dipMag(1:ncol,1:pver)
    dipMagD => istate%dipMagD(1:ncol,1:pver)

    ndensN2 => istate%ndensN2(1:ncol,1:pver)

    ionPRates    => istate%ionPRates(1:ncol,1:pver,1:nIonRates)
    sumIonPRates => istate%sumIonPRates(1:ncol,1:pver)
         
    sourceg4 => istate%sourceg4(1:ncol,1:pver)

    !-------------------------------------------------------------------------------------
    !  Calculate neutral temperature on interface levels.  tN vertical dimension is pver
    !-------------------------------------------------------------------------------------   
    do iVer = 2, pver
 
      do iCol = 1, ncol

        tNInt(iCol,iVer) = 0.5_r8 * tN(iCol,iVer) + 0.5_r8 * tN(iCol,iVer-1)

      enddo
    enddo

    do iCol = 1, ncol
        tNInt(iCol,1) = 1.5_r8 * tNInt(iCol,2) - 0.5_r8 * tNInt(iCol,3) 
    enddo
    do iCol = 1, ncol
        tNInt(iCol,pverp) = 1.5_r8 * tNInt(iCol,pver) - 0.5_r8 * tNInt(iCol,pver-1) 
    enddo

    !--------------------------------------------------------------
    !  Get zenith angle
    !-------------------------------------------------------------- 
    calDay = get_curr_calday()    
    call zenith(calDay,geoLatR(1:ncol),geoLonR(1:ncol),cosZenAngR(1:ncol),ncol)

    do iCol = 1, ncol

      zenAngD(iCol) = ACOS(cosZenAngR(iCol)) * rads2Degs
    
    enddo

    !--------------------------------------------------------------
    !  Get F10.7 solar flux
    !--------------------------------------------------------------
    call solar_parms_get( f107_s = f107 )

    !---------------------------------------------------------------------------------------
    !  Expand magnetic field components in vertical to make 3D, pcols,pver,begchunk:endchunk
    !  These are used in calculation of magnetic dip angle and magnetic declination angle so 
    !  store in local ionosphere module structure.
    !---------------------------------------------------------------------------------------
    do iVer = 1, pver
    
      do iCol = 1, ncol

        bNorth3d(iCol,iVer) = bnorth(iCol,lchnk)
        bEast3d(iCol,iVer) = beast(iCol,lchnk)
        bDown3d(iCol,iVer) = bdown(iCol,lchnk)

      enddo
    
    enddo

    !---------------------------------------------------------------------------------
    ! Get ion ExB drift from physics buffer (they were defined by exbdrift module)
    ! Ion drifts are 2d output arrays, i.e., no vertical dimension but 1d here (pcols)
    !---------------------------------------------------------------------------------
    indxUe = pbuf_get_index( 'UE' )
    indxVe = pbuf_get_index( 'VE' )
    indxWe = pbuf_get_index( 'WE' )
    call pbuf_get_field(pbuf, indxUe, ue)
    call pbuf_get_field(pbuf, indxVe, ve)
    call pbuf_get_field(pbuf, indxWe, we)

    !-------------------------------------------------------------------------------
    !  Form 2D drifts needed for this module.  Vertical level goes to teTiBotP since
    !  needed below to get vertical drift on interface levels
    !-------------------------------------------------------------------------------
    do iVer = 1, pver
    
      do iCol = 1, ncol

        ue2d(iCol,iVer) = ue(iCol)
        ve2d(iCol,iVer) = ve(iCol)
        we2d(iCol,iVer) = we(iCol)

      enddo
    
    enddo
     
    !-------------------------------------------------------------------------------
    !  Need vertical ExB drift on interface levels
    !-------------------------------------------------------------------------------
    do iVer = 2, pver
      do iCol = 1, ncol
        wei2d(iCol,iVer) = 0.5_r8 * (we2d(iCol,iVer-1) + we2d(iCol,iVer))
      enddo    
    enddo
    do iCol = 1, ncol
        wei2d(iCol,1) = 1.5_r8 * wei2d(iCol,2) - 0.5_r8 * wei2d(iCol,3)
    enddo
    do iCol = 1, ncol
        wei2d(iCol,pverp) = 1.5_r8 * wei2d(iCol,pver) - 0.5_r8 * wei2d(iCol,pver-1)
    enddo

    !--------------------------------------------------------------------------------------
    !  Need to get vertical velocity on interface levels handling top level specifically
    !--------------------------------------------------------------------------------------
    omega(1:ncol,1:pver) = state%omega(1:ncol,1:pver)
    
    do iVer = 2, pver
      do iCol = 1, ncol
        omegai(iCol,iVer) = 0.5_r8 * (omega(iCol,iVer-1) + omega(iCol,iVer))
      enddo
    enddo
    do iCol = 1, ncol
        omegai(iCol,1) = 1.5_r8 * omegai(iCol,2) - 0.5_r8 * omegai(iCol,3)
    enddo
    do iCol = 1, ncol
        omegai(iCol,pverp) = 1.5_r8 * omegai(iCol,pver) - 0.5_r8 * omegai(iCol,pver-1)
    enddo

    !------------------------------------------------------------------------
    !  Get constituent dependent gas constant and derive on interface levels
    !------------------------------------------------------------------------        
    do iVer = 2, pver
      do iCol = 1, ncol
        rairvi(iCol,iVer) = 0.5_r8 * rairv(iCol,iVer-1,lchnk) + 0.5_r8 * rairv(iCol,iVer,lchnk)
      enddo
    enddo

    do iCol = 1, ncol
       rairvi(iCol,1) = 1.5_r8 * rairvi(iCol,2) - 0.5_r8 * rairvi(iCol,3)
    enddo
    do iCol = 1, ncol
       rairvi(iCol,pverp) = 1.5_r8 * rairvi(iCol,pver) - 0.5_r8 * rairvi(iCol,pver-1)
    enddo

    !-------------------------------------------------------------------------------
    !  Need to get dip angle from magnetic field components
    !-------------------------------------------------------------------------------     
    do iVer = 1, pver
      do iCol = 1, ncol
        dipMag(iCol,iVer) = ATAN(bDown3d(iCol,iVer) / SQRT(bNorth3d(iCol,iVer)**2 + bEast3d(iCol,iVer)**2))
        if (dipMag(iCol,iVer) < 0.17_r8 .and. dipMag(iCol,iVer) > 0._r8 ) dipMag(iCol,iVer) = 0.17_r8
        if (dipMag(iCol,iVer) > -0.17_r8 .and. dipMag(iCol,iVer) < 0._r8 ) dipMag(iCol,iVer) = 0.17_r8
        dipMagD(iCol,iVer) = dipMag(iCol,iVer) * rads2Degs
      enddo
    enddo

    !-------------------------------------------------------------------------------------------
    !  Set up constituents to be accessed here from pbuf or state%q. 
    !-------------------------------------------------------------------------------------------
    cCnst = (/'O  ','O2 ','NO ','H  ','N  ','e  ','Op ','O2p','NOp'/)

    do iCnst = 1, nCnst

      !--------------------------------------
      !  Assign density to istate array
      !-------------------------------------- 
      if (cCnst(iCnst) == 'O  ') ndens => istate%ndensO1(1:ncol,1:pver)
      if (cCnst(iCnst) == 'O2 ') ndens => istate%ndensO2(1:ncol,1:pver) 
      if (cCnst(iCnst) == 'NO ') ndens => istate%ndensNO(1:ncol,1:pver) 
      if (cCnst(iCnst) == 'N  ') ndens => istate%ndensN1(1:ncol,1:pver) 
      if (cCnst(iCnst) == 'e  ') ndens => istate%ndensE(1:ncol,1:pver)
      if (cCnst(iCnst) == 'Op ') ndens => istate%ndensOp(1:ncol,1:pver)
      if (cCnst(iCnst) == 'O2p') ndens => istate%ndensO2p(1:ncol,1:pver)
      if (cCnst(iCnst) == 'NOp') ndens => istate%ndensNOp(1:ncol,1:pver)

      !-------------------------------------------------------------------------------------------
      !  Set flag and get field mmr whether each constituent is short-lived(pbuf) or not(state%q). 
      !-------------------------------------------------------------------------------------------
      call cnst_get_ind( TRIM(cCnst(iCnst)), indxCnst, abort=.false. )
      if (indxCnst < 0) then
         indxSlvd = slvd_index( TRIM(cCnst(iCnst)) )
  	 if (indxSLvd > 0) then
  	    call pbuf_get_field(pbuf, slvd_pbf_ndx, mmrP, start=(/1,1,indxSLvd/), kount=(/pcols,pver,1/) )
  	    mmr(1:ncol,1:pver) = mmrP(1:ncol,1:pver)
            sIndx = get_spc_ndx( TRIM(cCnst(iCnst)) )
  	    rMass  = adv_mass(sIndx)
  	 endif
      else
  	 mmr(1:ncol,1:pver) = state%q(1:ncol,1:pver,indxCnst)
  	 rMass  = cnst_mw(indxCnst)
      endif

      !--------------------------------------------------------------------------------------------------------------
      !  Need to get number density (cgs units) from mass mixing ratio.  mbarv is kg/mole, same as rMass units
      !  kg/kg * (kg/mole)/(kg/mole) * (Pa or N/m*m)/((Joules/K or N*m/K) * (K)) = m-3 * 1E-06 = cm-3
      !--------------------------------------------------------------------------------------------------------------- 
      ndens(1:ncol,1:pver)  = mmr(1:ncol,1:pver) * mbarv(1:ncol,1:pver,lchnk) / rMass * &
                                   pMid(1:ncol,1:pver) / (kboltz * tN(1:ncol,1:pver)) * 1.E-06_r8
 
      if (cCnst(iCnst) == 'O  ') then
        mmrO1(1:ncol,1:pver) = mmr(1:ncol,1:pver)
	ndensO1(1:ncol,1:pver) = ndens(1:ncol,1:pver)
      endif
      if (cCnst(iCnst) == 'O2 ') then
        mmrO2(1:ncol,1:pver) = mmr(1:ncol,1:pver)
	ndensO2(1:ncol,1:pver) = ndens(1:ncol,1:pver)
      endif
      if (cCnst(iCnst) == 'N  ') ndensN1(1:ncol,1:pver) = ndens(1:ncol,1:pver)
      if (cCnst(iCnst) == 'e  ') ndensE(1:ncol,1:pver) = ndens(1:ncol,1:pver)

      !----------------------------------------------------------------------------
      !  Calculate N2 density from O2 and O and assign to istate array
      !----------------------------------------------------------------------------
      if (iCnst == nCnst) then

        rMassN2 = 28._r8     
        mmrN2(1:ncol,1:pver) = 1._r8 - (mmrO2(1:ncol,1:pver) + mmrO1(1:ncol,1:pver)) 
        mmrN2(1:ncol,1:pver) = MAX(1.e-20_r8,mmrN2(1:ncol,1:pver))
        ndensN2(1:ncol,1:pver) = mmrN2(1:ncol,1:pver) * mbarv(1:ncol,1:pver,lchnk) / rMassN2 * &
				            pMid(1:ncol,1:pver) / (kboltz * tN(1:ncol,1:pver)) * 1.E-06_r8		
	 
      endif
 
    enddo ! nCnst

    !------------------------------------------------------------------------------------
    ! Get ionization rates from physics buffer which were calculated in mo_jeuv and 
    ! mo_jshort modules.  Rates array dimensions are pcols, pver, nIonRates.  Units s-1
    !------------------------------------------------------------------------------------
    indxIR = pbuf_get_index( 'IonRates' )
    call pbuf_get_field(pbuf, indxIR, ionRates)

    !----------------------------------------------------------------------------------------------
    !  Need to convert these ionization rates to ion production rates by multiplying number density 
    !  of neutral species appropriate from reactions in mo_jeuv(jeuv) and mo_jshort(jshort)(for NO)  
    !----------------------------------------------------------------------------------------------         
    do iVer = 1, pver
      do iCol = 1, ncol
    
        do iIonR = 1, nIonRates
          IF (iIonR <= 3) ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensO1(iCol,iVer)
          IF (iIonR == 4) ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensN1(iCol,iVer)
          IF ((iIonR == 5) .OR. (iIonR >= 7 .AND. iIonR <= 9)) &
                                    ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensO2(iCol,iVer)
          IF (iIonR == 6 .OR. iIonR == 10 .OR. iIonR == 11) &
                                    ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensN2(iCol,iVer)
        enddo
                                    
        !----------------------------------------------
        !  Sum ion production rates all reactions
        !----------------------------------------------   
        sumIonPRates(iCol,iVer) = SUM(ionPRates(iCol,iVer,1:11))

       enddo
    enddo
        
    !-------------------------------------------------------------------------------------------
    ! Get aurora ion production rate sum from physics buffer which were calculated in mo_aurora 
    ! module.  Rate array dimensions are pcols, pver.  Units s-1 cm-3
    !-------------------------------------------------------------------------------------------
    indxAIPRS = pbuf_get_index( 'AurIPRateSum' )
    call pbuf_get_field(pbuf, indxAIPRS, aurIPRateSum)

    !-------------------------------------------------------------------------------------
    ! Check latitudes, and set aurora ion production rates for all columns and levels 
    ! to zero if all equatorward of 32.5 deg since not calculated in mo_aurora for those 
    ! latitudes (same as criteria in mo_aurora)
    !-------------------------------------------------------------------------------------
    do_aurora(1:ncol) = abs( geoLatR(1:ncol) ) > pi/6._r8
    if( all( .not. do_aurora(1:ncol) ) ) then
       aurIPRateSum(:,:) = 0._r8
    end if
 
    !-------------------------------------------------------------------------------------------------
    !  Calculate electron heating rate which is a source in electron/ion temperature derivation
    !-------------------------------------------------------------------------------------------------
    do iVer = 1, teTiBot
      do iCol = 1, ncol
        sourceR(iCol,iVer) = LOG( ndensE(iCol,iVer) / (ndensO2(iCol,iVer) + ndensN2(iCol,iVer) + &
                                                                                  0.1_r8 * ndensO1(iCol,iVer)) )
        sourceEff(iCol,iVer) = EXP( -(12.75_r8 + 6.941_r8 * sourceR(iCol,iVer) + 1.166_r8 * sourceR(iCol,iVer)**2 + &
                                        0.08043_r8 * sourceR(iCol,iVer)**3 + 0.001996_r8 * sourceR(iCol,iVer)**4) )

        !-------------------------------------------------------------------------------
        !  Calculate g4 source term for electron temperature update
        !-------------------------------------------------------------------------------         
        sourceg4(iCol,iVer) = (sumIonPRates(iCol,iVer) + aurIPRateSum(iCol,iVer)) * sourceEff(iCol,iVer)

      enddo

    enddo

    !----------------------------------------------------------------------------------------------
    ! Get Pedersen and Hall Conductivities from physics buffer which were calculated in iondrag 
    ! module.  Conductivity array dimensions are pcols, pver
    !-------------------------------------------------------------------------------
    indxSP = pbuf_get_index( 'PedConduct' )
    indxSH = pbuf_get_index( 'HallConduct' )
    call pbuf_get_field(pbuf, indxSP, sigma_ped)
    call pbuf_get_field(pbuf, indxSH, sigma_hall)

    return

  end subroutine ionos_timestep_init
!
!===============================================================================

  subroutine ionos_teti(state, dSETendIn, dSETendOut, ztodt, istate, tE, tI, teTiBot)

  !-----------------------------------------------------------------------
  ! Routine to compute the electron and ion temperature needed for
  ! ambipolar diffusion calculation.
  !-----------------------------------------------------------------------

    use phys_grid,       only : get_lat_p, get_lon_p, get_rlat_p,get_rlon_p
    use physconst,       only : gravit ! Gravity (m/s2)
    use physconst,       only : rairv, mbarv                     ! Constituent dependent rair and mbar

!------------------------------Arguments--------------------------------

    type(physics_buffer_desc), pointer           :: pbuf(:)  ! physics buffer
    type(physics_state),   intent(in), target    :: state    ! physics state structure
    type(ionos_state),     intent(in), target    :: istate   ! ionosphere state structure

    real(r8), dimension(pcols,pver),   intent(in)      :: dSETendIn    ! dry static energy tendency
    real(r8), dimension(pcols,pver),   intent(out)     :: dSETendOut   ! dry static energy tendency

    real(r8), intent(in)                               :: ztodt     ! physics time step

    real(r8), dimension(:,:), pointer, intent(inout)   :: tE        ! Pointer to electron temperature in pbuf (K) 
    real(r8), dimension(:,:), pointer, intent(inout)   :: ti        ! Pointer to ion temperature in pbuf (K) 

    integer, intent(in) :: teTiBot  ! bottom of ionosphere calculations 

!---------------------------Local storage-------------------------------
    integer :: lchnk                                    ! Chunk number 
    integer :: ncol                                     ! Number of atmospheric columns 
    integer :: teTiBotP                                 ! bottom of ionosphere calculations plus one more level 

    integer :: i,k                                      ! loop indexes
    integer :: iVer                                     ! Counter for vertical loops
    integer :: iCol                                     ! Counter for column loops
    integer :: iter                                     ! Counter for iteration loop

    integer, parameter  :: maxIter  = 6                 ! maximum number of iterations to solve for electron/ion temperature

    real(r8), parameter :: Kec1   = 7.5E5_r8            ! c1 constant for calculation of electron conductivity(Ke)
    real(r8), parameter :: Kec2   = 3.22E4_r8           ! c2 constant for calculation of electron conductivity(Ke)
    real(r8), parameter :: stepweight  = 1.0_r8         ! weight of previous and current times step for diagonals
    real(r8), parameter :: sToQConv  = 6.24E15_r8       ! Conversion from J/kg/s to ev/g/s

    real(r8), parameter :: lossc5  = 1.21E-4_r8         ! c5 constant needed for loss term g3 for electron temperature update
    real(r8), parameter :: lossc7  = 3.6E-2_r8          ! c7 constant needed for loss term g3 for electron temperature update
    real(r8), parameter :: lossc9  = 5.7E-4_r8          ! c9 constant needed for loss term g3 for electron temperature update
    real(r8), parameter :: lossc13 = 7.E-5_r8           ! c13 constant needed for loss term g3 for electron temperature update

    real(r8), parameter :: lossc4pCoef  = 1.77E-19_r8 
    real(r8), parameter :: lossc6pCoef  = 1.21E-18_r8 
    real(r8), parameter :: lossc8pCoef  = 7.9E-19_r8
    real(r8), parameter :: lossc10pCoef = 1.3E-4_r8 
    real(r8), parameter :: lossc11pCoef = 3.125E-21_r8
    real(r8), parameter :: lossc12pCoef = 3.4E-12_r8
    real(r8), parameter :: lossc14pCoef = 1.57E-12_r8 
    real(r8), parameter :: lossc15pCoef = 2.9E-14_r8
    real(r8), parameter :: lossc16pCoef = 6.9E-14_r8
    real(r8), parameter :: lossc3pC1 = 3.2E-8_r8
    real(r8), parameter :: lossc3pC2 = 15._r8
    real(r8), parameter :: lossc3pC3 = 0.53_r8

    real(r8), parameter :: losscinCoef1 = 6.6e-14_r8
    real(r8), parameter :: losscinCoef2 = 5.8e-14_r8
    real(r8), parameter :: losscinCoef3 = 0.21e-14_r8
    real(r8), parameter :: losscinCoef4 = 5.9e-14_r8
    real(r8), parameter :: losscinCoef5 = 5.45e-14_r8
    real(r8), parameter :: losscinCoef6 = 4.5e-14_r8
    real(r8), parameter :: losscinCoef7 = 5.8e-14_r8
    real(r8), parameter :: losscinCoef8 = 0.14e-14_r8
    real(r8), parameter :: losscinCoef9 = 4.4e-14_r8
    
    real(r8), parameter :: FeDCoef1 = -5.0E+7_r8
    real(r8), parameter :: FeDCoef2 = 4.0E+7_r8

    real(r8), parameter :: losscACoef1 = 5.71E-8_r8
    real(r8), parameter :: losscACoef2 = -3352.6_r8
    real(r8), parameter :: losscACoef3 = 2.0E-7_r8
    real(r8), parameter :: losscACoef4 = -4605.2_r8
    real(r8), parameter :: losscACoef5 = 2.53E-6_r8 
    real(r8), parameter :: losscACoef6 = -17620._r8

    real(r8), parameter :: loss10pCoef = 3200._r8
    real(r8), parameter :: lossc12pC1  = 0.4_r8
    real(r8), parameter :: lossc12pC2  = 150._r8

    real(r8), parameter :: losscf2dC1 = 2.4E+4_r8
    real(r8), parameter :: losscf2dC2 = 0.3_r8
    real(r8), parameter :: losscf2dC3 = 1500._r8
    real(r8), parameter :: losscf2dC4 = 1.947E-5_r8
    real(r8), parameter :: losscf2dC5 = 4000._r8

    real(r8), parameter :: losscf2C1 = 3000._r8

    real(r8), parameter :: losscf3c1 = -22713._r8

    real(r8), parameter :: f1Ted1C1 = 2.82E-17_r8
    real(r8), parameter :: f1Ted1C2 = 3.41E-21_r8

    real(r8), parameter :: f1Ted2C1 = 2.2E-16_r8
    real(r8), parameter :: f1Ted2C2 = 7.92E-18_r8

    real(r8), parameter :: f1Ted3C1 = 1.1E-16_r8
    real(r8), parameter :: f1Ted3C2 = 5.7E-4_r8

    real(r8) :: wrk1                                    ! 2/3/kboltz_ev
    real(r8) :: FeDB                                    ! B term of electron heat flux of UB
    real(r8) :: FeD                                     ! Day time flux
    real(r8) :: FeN                                     ! Night time flux
    real(r8) :: f1Ted1                                  ! d1 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Ted2                                  ! d2 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Ted3                                  ! d3 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Te

    real(r8), dimension(:,:), pointer  	:: pMid 	! Midpoint pressure (Pa)
    real(r8), dimension(:,:), pointer 	:: tN		! Neutral temperature (K)
    real(r8), dimension(pcols,pver)  	:: tE0  	! Electron temperature from last time step (K)
    real(r8), dimension(pcols,pver)  	:: tEPrevI	! Electron temperature from previous iteration (K)

    real(r8), dimension(:,:), pointer 	:: pInt 	! Interface pressure (Pa)
    real(r8), dimension(:,:), pointer 	:: tNInt	! Interface Temperture (K)
    real(r8), dimension(:,:), pointer	:: rairvi	! Constituent dependent gas constant on interface levels

    real(r8), dimension(:,:), pointer	 :: ndensN2	! N2 number density (cm-3)
    real(r8), dimension(:,:), pointer	 :: ndensO2	! O2 number density (cm-3)
    real(r8), dimension(:,:), pointer	 :: ndensO1	! O number density (cm-3)
    real(r8), dimension(:,:), pointer	 :: ndensNO	! NO number density (cm-3)
    real(r8), dimension(:,:), pointer	 :: ndensN1	! N number density  (cm-3)
    real(r8), dimension(:,:), pointer	 :: ndensE	! E electron number density (cm-3)
    real(r8), dimension(:,:), pointer	 :: ndensOp	! O plus number density (cm-3)
    real(r8), dimension(:,:), pointer	 :: ndensO2p	! O2 plus ion number density (cm-3)
    real(r8), dimension(:,:), pointer	 :: ndensNOp	! NO plus ion number density  (cm-3)

    real(r8), dimension(:,:), pointer  	:: sourceg4	! g4 source term for electron/ion temperature update
 
    real(r8), dimension(:,:), pointer	:: dipMag       ! dip angle for each column (radians)
    real(r8), dimension(:,:), pointer	:: dipMagD	! dip angle for each column (degrees)

    real(r8), dimension(:),   pointer   :: zenAngD	! zenith angle (degrees)

    real(r8), dimension(pcols)       	:: FeUB 	! electron heat flux at upper boundary
 
    real(r8), dimension(pver)   	:: sqrtTE	! Square root of electron temperature
 
    real(r8), dimension(pver)        	:: Ke		! electron conductivity

    real(r8), dimension(pverp)       	:: Kei  	! electron conductivity interface levels

    real(r8), dimension(pcols,pver)  	:: lossc4p	! c4 prime of Lc(eN2) component of loss term
    real(r8), dimension(pcols,pver)  	:: lossceN2	! Lc(eN2) component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc6p	! c6 prime of Lc(eO2) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO2	! Lc(eO2) component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc8p	! c8 prime of Lc(eO) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO1	! Lc(eO) component of loss term equation

    real(r8), dimension(pcols,pver)  	:: lossc10p	! c10 prime of Lc(eN2) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscA	! A of Lc(eN2)v component of loss term equation
    real(r8), dimension(pcols,pver)  	:: tENDiff	! Difference between electron and neutral temperatures
    real(r8), dimension(pcols,pver)  	:: lossceN2v	! Lc(eN2)v component of loss term equation

    real(r8), dimension(pcols,pver)  	:: lossc11p	! c11 prime of Lc(eO2)v component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO2v	! Lc(eO2)v component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc12p	! c12 prime of Lc(eO)f component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceOf	! Lc(eO)f component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc14p	! c14 prime of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscf2d	! d of f2 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscf2	! f2 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscf3	! f3 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO1D	! Lc(eO)1D component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc15p	! c15 prime of Lc(eN2)Rot component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceN2Rot  ! Lc(eN2)Rot component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc16p	! c16 prime of Lc(eO2)Rot component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO2Rot  ! Lc(eO2)Rot component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc3p	! c3 prime of Lc(ei) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscei	! Lc(ei) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscin	! ion-neutral heating coeff.

    real(r8), dimension(pcols,pver)  	:: lossg3	! g3 loss term for Te tendency

    real(r8), dimension(pcols,pver)  	:: delTEN	! Difference between electron and neutral temperatures from production/loss

    real(r8), dimension(pcols,pverp) 	:: delZi	! Delta z: interfaces
    real(r8), dimension(pcols,pver)  	:: delZ 	! Delta z: midpoints
 
    real(r8), dimension(pcols,pver)  	:: qjoule	! joule heating
    real(r8), dimension(pcols,pver)  	:: qen  	! electron-neutral heating
    real(r8), dimension(pcols,pver)  	:: qei  	! electron-ion Coulomb heating
    real(r8), dimension(pcols,pver)  	:: rho  	! mass density

    real(r8), dimension(pcols,pver)  	:: wrk2

    real(r8), dimension(teTiBot)      	:: subdiag	! subdiagonal values for Te tendency solving
    real(r8), dimension(teTiBot)      	:: superdiag	! superdiagonal values for Te tendency solving
    real(r8), dimension(teTiBot)      	:: diag 	! diagonal values for Te tendency solving
    real(r8), dimension(teTiBot)      	:: rHSInit	! initial RHS of electron temperature update
    real(r8), dimension(teTiBot)      	:: rHSH 	! h for RHS of electron temperature update
    real(r8), dimension(teTiBot)      	:: rHS  	! RHS of electron temperature update
    real(r8), dimension(teTiBot)      	:: tETemp	! temporary electron temperature array for input to tridag

    logical, dimension(pcols)        	:: colConv	! flag for column converging = 1 if converged otherwise = 0

    logical                          	:: converged	! Flag for convergence in electron temperature calculation iteration loop
        
!-----------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------- 
    !  Initialize arrays to zero and column convergence logical to .false.
    !---------------------------------------------------------------------------------------------------------
    sqrtTE(:)           = 0._r8
    Ke(:)               = 0._r8
    Kei(:)              = 0._r8
    lossc4p(:,:)        = 0._r8
    lossceN2(:,:)       = 0._r8
    lossc6p(:,:)        = 0._r8
    lossc6p(:,:)        = 0._r8
    lossceO2(:,:)       = 0._r8
    lossc8p(:,:)        = 0._r8
    lossceO1(:,:)       = 0._r8
    lossc10p(:,:)       = 0._r8
    losscA(:,:)         = 0._r8
    tENDiff(:,:)        = 0._r8
    lossceN2v(:,:)      = 0._r8
    lossc11p(:,:)       = 0._r8
    lossceO2v(:,:)      = 0._r8
    lossc12p(:,:)       = 0._r8
    lossceOf(:,:)       = 0._r8
    lossc14p(:,:)       = 0._r8
    losscf2d(:,:)       = 0._r8
    losscf2(:,:)        = 0._r8
    losscf3(:,:)        = 0._r8
    lossceO1D(:,:)      = 0._r8
    lossc15p(:,:)       = 0._r8
    lossceN2Rot(:,:)    = 0._r8
    lossc16p(:,:)       = 0._r8
    lossceO2Rot(:,:)    = 0._r8
    lossc3p(:,:)        = 0._r8
    losscei(:,:)        = 0._r8
    losscin(:,:)        = 0._r8
    lossg3(:,:)         = 0._r8
    delTEN(:,:)         = 0._r8 
    delZi(:,:)          = 0._r8
    delZ(:,:)           = 0._r8 
    subDiag(:)          = 0._r8         
    superDiag(:)        = 0._r8        
    diag(:)             = 0._r8 
    rHSInit(:)          = 0._r8 
    rHSH(:)             = 0._r8
    rHS(:)              = 0._r8 
    teTemp(:)           = 0._r8 
    qjoule(:,:)         = 0._r8
    qei(:,:)            = 0._r8
    qen(:,:)            = 0._r8
    rho(:,:)            = 0._r8
    dSETendOut          = 0._r8
    colConv(:)          = .false.
    
    !--------------------------------------------------------------------------------------
    !  Get lchnk and ncol from state
    !--------------------------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol
    
    !-------------------------------------------
    !  Calculate some commonly used variables
    !-------------------------------------------
    wrk1 = 2._r8 / 3._r8/ kboltz_ev
    teTiBotP = teTiBot + 1

    !-------------------------------------------------------------------------------------------------------
    !  Need to get midpoint and interface pressure and neutral temperature from state structure (ncol,teTiBot)
    !-------------------------------------------------------------------------------------------------------
    pMid  => state%pmid(1:ncol,1:pver)
    tN    => state%t(1:ncol,1:pver)
    rho(1:ncol,1:pver) = pMid(1:ncol,1:pver)/rairv(1:ncol,1:pver,lchnk)/tN(1:ncol,1:pver) * 1.E-3_r8     ! convert to g/cm3

    qjoule(1:ncol,1:teTiBot) = dSETendIn(1:ncol,1:teTiBot) * sToQConv     ! convert from J/kg/s to ev/g/s

    pInt    => state%pint(1:ncol,1:pverp)
    tNInt   => istate%tNInt(1:ncol,1:pverp)        
    rairvi  => istate%rairvi(1:ncol,1:pverp)

    !----------------------------------------------------------------
    !  Get variables needed from the ionosphere state structure
    !----------------------------------------------------------------
    ndensO2  => istate%ndensO2(1:ncol,1:pver) 
    ndensO1  => istate%ndensO1(1:ncol,1:pver)
    ndensNO  => istate%ndensNO(1:ncol,1:pver) 
    ndensN1  => istate%ndensN1(1:ncol,1:pver) 
    ndensE   => istate%ndensE(1:ncol,1:pver)  
    ndensOp  => istate%ndensOp(1:ncol,1:pver) 
    ndensO2p => istate%ndensO2p(1:ncol,1:pver)
    ndensNOp => istate%ndensNOp(1:ncol,1:pver)
    ndensN2  => istate%ndensN2(1:ncol,1:pver) 

    sourceg4 => istate%sourceg4(1:ncol,1:pver)

    dipMag   => istate%dipMag(1:ncol,1:pver)
    dipMagD  => istate%dipMagD(1:ncol,1:pver)

    zenAngD  => istate%zenAngD(1:ncol) 
      
    !-------------------------------------------------------------------------------------------------------------------
    !  Get electron temperature from physics buffer and if this is the first time calculated then initialize to neutral
    !  temperature.  is_first_step means this is the first time step of an initial run. 
    !-------------------------------------------------------------------------------------------------------------------
    if( is_first_step() ) then   

      tI(1:ncol,1:pver) = tN(1:ncol,1:pver)
      tE(1:ncol,1:pver) = tN(1:ncol,1:pver)
    
    else

      tE(1:ncol,1:pver) = MAX(tN(1:ncol,1:pver),tE(1:ncol,1:pver))
      tE(1:ncol,1:pver) = MIN(temax,tE(1:ncol,1:pver))

      tI(1:ncol,1:pver) = MAX(tN(1:ncol,1:pver),ti(1:ncol,1:pver))
      tI(1:ncol,1:pver) = MIN(ti(1:ncol,1:pver),tE(1:ncol,1:pver))

      tE(1:ncol,teTiBotP:pver) = tN(1:ncol,teTiBotP:pver)
      tI(1:ncol,teTiBotP:pver) = tN(1:ncol,teTiBotP:pver)

    endif

    tE0(1:ncol,1:pver) = tE(1:ncol,1:pver)

    wrk2(1:ncol,1:teTiBot) =  ndensE(1:ncol,1:teTiBot)/wrk1/(SIN(dipMag(1:ncol,1:teTiBot)))**2._r8
    
    !-----------------------------------------------------------------------------
    !  Get terms needed for loss term g3 for electron temperature update which do 
    !  not need to be updated in iteration loop.  
    !-----------------------------------------------------------------------------
    do iCol = 1, ncol

      if (.not. colConv(iCol)) then

        do iVer = 1, teTiBot

	  lossc4p(iCol,iVer)  = lossc4pCoef * ndensN2(iCol,iVer) * ndensE(iCol,iVer)
	  lossc6p(iCol,iVer)  = lossc6pCoef * ndensO2(iCol,iVer) * ndensE(iCol,iVer)
	  lossc8p(iCol,iVer)  = lossc8pCoef * ndensO1(iCol,iVer) * ndensE(iCol,iVer)
	  lossc10p(iCol,iVer) = lossc10pCoef * ndensN2(iCol,iVer) * ndensE(iCol,iVer)
	  lossc11p(iCol,iVer) = lossc11pCoef * ndensO2(iCol,iVer) * ndensE(iCol,iVer)
	  lossc12p(iCol,iVer) = lossc12pCoef * ndensO1(iCol,iVer) * ndensE(iCol,iVer)
	  lossc14p(iCol,iVer) = lossc14pCoef * ndensO1(iCol,iVer) * ndensE(iCol,iVer)
	  lossc15p(iCol,iVer) = lossc15pCoef * ndensN2(iCol,iVer) * ndensE(iCol,iVer)
	  lossc16p(iCol,iVer) = lossc16pCoef * ndensO2(iCol,iVer) * ndensE(iCol,iVer)
	  lossc3p(iCol,iVer)  = lossc3pC1 * lossc3pC2 * (ndensOP(iCol,iVer) + &
			      0.5_r8 * ndensO2P(iCol,iVer) + lossc3pC3 * ndensNOP(iCol,iVer)) * ndensE(iCol,iVer)

	  losscin(iCol,iVer) = (losscinCoef1*ndensN2(iCol,iVer) + losscinCoef2*ndensO2(iCol,iVer)		     &
			     + losscinCoef3*ndensO1(iCol,iVer)*SQRT(2._r8)*SQRT(tN(iCol,iVer)))*ndensOP(iCol,iVer)    &
			     +(losscinCoef4*ndensN2(iCol,iVer) + losscinCoef5*ndensO2(iCol,iVer)		   &
			     + losscinCoef6*ndensO1(iCol,iVer))*ndensOP(iCol,iVer)			       &
			     +(losscinCoef7*ndensN2(iCol,iVer) + losscinCoef8*ndensO2(iCol,iVer)*SQRT(tN(iCol,iVer)) &
			     + losscinCoef9*ndensO1(iCol,iVer)) * ndensO2P(iCol,iVer)


        enddo !iVer loop

        !----------------------------------------------------------------------------------
        !  Calculate upper boundary heat flux 
        !----------------------------------------------------------------------------------
        if (ABS(dipMagD(iCol,1)) < 40.0_r8) FeDB = 0.5_r8 * &
                                        (1._r8 + SIN(pi * (ABS(dipMagD(iCol,1)) - 20.0_r8) /40.0_r8))

        if (ABS(dipMagD(iCol,1)) >= 40.0_r8) FeDB = 1._r8 

        FeD = FeDCoef1 * f107 * FeDB - FeDCoef2 * f107
        FeN = .5_r8 * FeD

        !---------------------------------------------------
        !  Set upper boundary condition for right hand side
        !---------------------------------------------------
        if (zenAngD(iCol) <= 80.0_r8) FeUB(iCol) = FeD
        if (zenAngD(iCol) > 80.0_r8 .AND. zenAngD(iCol) < 100.0_r8) FeUB(iCol) = 0.5_r8 * (FeD + FeN) &
                                                                         + 0.5_r8 * (FeD - FeN) * &
                                                                 COS(pi * ((zenAngD(iCol) - 80.0_r8) / 20.0_r8))
        if (zenAngD(iCol) >= 100.0_r8) FeUB(iCol) = FeN

        !------------------------------------------------------------------------------------------
        !  Calculate thickness terms for vertical derivative
        !------------------------------------------------------------------------------------------
        do iVer = 1, teTiBot

          delZ(iCol,iVer) = (pInt(iCol,iVer+1) - pInt(iCol,iVer)) * rairv(iCol,iVer,lchnk) * &
	                                                               tN(iCol,iVer) / pMid(iCol,iVer) / gravit

        enddo

        do iVer = 2, teTiBotP		! Assuming teTiBotP < pverp
          delZi(iCol,iVer) = (pMid(iCol,iVer) - pMid(iCol,iVer-1)) * rairvi(iCol,iVer) * &
        							  tNInt(iCol,iVer) / pInt(iCol,iVer) / gravit
        enddo
        delZi(iCol,1) = 1.5_r8*delZi(iCol,2) - .5_r8*delZi(iCol,3)

        !----------------------------------------------------------
        !  Convert delZ variables from meters to centimeters
        !----------------------------------------------------------
        delZi(iCol,1:teTiBotP) = delZi(iCol,1:teTiBotP)*100._r8
        delZ(iCol,1:teTiBot) = delZ(iCol,1:teTiBot)*100._r8
  
      endif ! Column not converged

    enddo !iCol loop

    !-------------------------------------------------------------------------------------------------------
    !  Iterate to calculate new electron temperature. 
    !  Time splitting is used: first solve the heating/cooling equation, then solve the diffusion equations.
    !  Also, set convergence flag to false and iterate until true or 6 iterations, whichever comes first 
    !-------------------------------------------------------------------------------------------------------
    converged = .false. 
    iter = 0
    do while (.not. converged .and. iter < maxIter)
    
      !--------------------------------------------------------------------------------------------------------
      !  Increment iteration loop counter and save electron temperature from previous iteration for convergence 
      !  test at end of interation loop.  Also, take square root of electron temperature to be used later
      !--------------------------------------------------------------------------------------------------------        
      iter = iter + 1
      
      tEPrevI(1:ncol,1:teTiBot) = tE(1:ncol,1:teTiBot)

      !--------------------------------------------------------------------------------------------------------
      !  Loop over columns then vertical levels and call tridiagonal solver for each column to get electron 
      !  temperature
      !--------------------------------------------------------------------------------------------------------
      do iCol = 1, ncol

        if (.not. colConv(iCol)) then

          sqrtTE(1:teTiBot) = SQRT(tE(iCol,1:teTiBot))

          do iVer = 1, teTiBot
 
            !-----------------------------------------------------------------------------
            !  Get loss term g3 for electron temperature update.  Need to calculate 
            !  constituent dependent loss terms which make up g3
            !-----------------------------------------------------------------------------
            lossceN2(iCol,iVer) = lossc4p(iCol,iVer) * (1._r8 - lossc5 * tE(iCol,iVer)) * tE(iCol,iVer)
            lossceO2(iCol,iVer) = lossc6p(iCol,iVer) * (1._r8 + lossc7 * sqrtTE(iVer)) * sqrtTE(iVer)
            lossceO1(iCol,iVer) = lossc8p(iCol,iVer) * (1._r8 + lossc9 * tE(iCol,iVer)) * sqrtTE(iVer)

            if (tE(iCol,iVer) < 1000.0_r8) losscA(iCol,iVer) = losscACoef1 * EXP(losscACoef2 / tE(iCol,iVer)) 
            if (tE(iCol,iVer) >= 1000.0_r8 .AND. tE(iCol,iVer) <= 2000.0_r8) &
                                          losscA(iCol,iVer) = losscACoef3 * EXP(losscACoef4 / tE(iCol,iVer)) 
            if (tE(iCol,iVer) > 2000.0_r8) losscA(iCol,iVer) = losscACoef5 * sqrtTE(iVer) * &
                                                                          EXP(losscACoef6 / tE(iCol,iVer))

            tENDiff(iCol,iVer) = tE(iCol,iVer) -  tN(iCol,iVer)
            if (ABS(tENDiff(iCol,iVer)) < 0.1_r8) tENDiff(iCol,iVer) = 0.1_r8

            lossceN2v(iCol,iVer) = lossc10p(iCol,iVer) * losscA(iCol,iVer) * &
               (1._r8 - EXP(loss10pCoef * (1._r8 / tE(iCol,iVer) - 1._r8 / tN(iCol,iVer)))) / tENDiff(iCol,iVer)
            lossceO2v(iCol,iVer) = lossc11p(iCol,iVer) * tE(iCol,iVer) * tE(iCol,iVer)
            lossceOf(iCol,iVer) = lossc12p(iCol,iVer) * (1._r8 - lossc13 * tE(iCol,iVer)) * &
                                               (lossc12pC1 + lossc12pC2 / tE(iCol,iVer)) / tN(iCol,iVer)
            losscf2d(iCol,iVer) = losscf2dC1 + losscf2dC2 * (tE(iCol,iVer) - losscf2dC3) + &
                          losscf2dC4 * (tE(iCol,iVer) - losscf2dC3) * (tE(iCol,iVer) - losscf2dC5)
            losscf2(iCol,iVer) = losscf2d(iCol,iVer) * (1._r8 / losscf2C1 - 1._r8 / tE(iCol,iVer))
            losscf3(iCol,iVer) = losscf3c1 * (1._r8 / tN(iCol,iVer) - 1._r8 / tE(iCol,iVer))
            lossceO1D(iCol,iVer) = lossc14p(iCol,iVer) * EXP(losscf2(iCol,iVer)) * & 
                                                   (1._r8 - EXP(losscf3(iCol,iVer))) / tENDiff(iCol,iVer)
            lossceN2Rot(iCol,iVer) = lossc15p(iCol,iVer) / sqrtTE(iVer)
            lossceO2Rot(iCol,iVer) = lossc16p(iCol,iVer) / sqrtTE(iVer)
            losscei(iCol,iVer) = lossc3p(iCol,iVer) / tE(iCol,iVer)**1.5_r8

            !------------------------------------------------
            ! Loss term: lossg3*tE/sin(I)^2
            !------------------------------------------------
            lossg3(iCol,iVer) =  lossceN2(iCol,iVer) + lossceO2(iCol,iVer) + lossceO1(iCol,iVer) + lossceN2v(iCol,iVer)   &
                             + lossceO2v(iCol,iVer) + lossceOf(iCol,iVer) + lossceO1D(iCol,iVer)                        &
                             + lossceN2Rot(iCol,iVer) + lossceO2Rot(iCol,iVer)

          enddo !iVer loop

	endif ! Column not converged

      enddo ! End of column loop

      !-----------------------------------------------------
      !  Calculate thermal conductivity of electron gas   
      !-----------------------------------------------------
      do iCol = 1, ncol

        if (.not. colConv(iCol)) then

          sqrtTE(1:teTiBot) = SQRT(tE(iCol,1:teTiBot))
      
          do iVer = 1, teTiBot

            f1Ted1 = f1Ted1C1 * sqrtTE(iVer) - f1Ted1C2 * tE(iCol,iVer)**1.5_r8
            f1Ted2 = f1Ted2C1 + f1Ted2C2 * sqrtTE(iVer)
            f1Ted3 = f1Ted3C1 * (1._r8 + f1Ted3C2 * tE(iCol,iVer))

            f1Te = ndensN2(iCol,iVer) / ndensE(iCol,iVer) * f1Ted1 + ndensO2(iCol,iVer) / &
                        ndensE(iCol,iVer) * f1Ted2 + ndensO1(iCol,iVer) / ndensE(iCol,iVer) * f1Ted3

            !-----------------------------------------------------------------------------
            !  Calculate electron conductivity using parameters set in module and f1(Te)
            !-----------------------------------------------------------------------------
            Ke(iVer) = Kec1 * tE(iCol,iVer)**2.5_r8 / (1._r8 + Kec2 * tE(iCol,iVer)**2._r8 * f1Te)

          enddo !iVer loop

          !----------------------------------------------------------------------
          !  Get electron conductivity at interface levels to be used later
          !----------------------------------------------------------------------
          do iVer = 2,teTiBot
            Kei(iVer) = SQRT(Ke(iVer-1)*Ke(iVer))
          enddo
          Kei(1) = 1.5_r8*Kei(2)-.5_r8*Kei(3)
          Kei(teTiBotP) = 1.5_r8*Kei(teTiBot)-.5_r8*Kei(teTiBot-1)

          !------------------------------------------------------------------------------------------------------
          !  Derive subdiagonal, superdiagonal, and diagonal as input to solver for electron temperature tendency
          !------------------------------------------------------------------------------------------------------
          do iVer = 2, teTiBot-1
              subDiag(iVer) = -Kei(iVer) / delZi(iCol,iVer) / delZ(iCol,iVer)
              superDiag(iVer) = -Kei(iVer+1) / delZi(iCol,iVer+1) / delZ(iCol,iVer)
              diag(iVer) = wrk2(iCol,iVer)/ztodt + (lossg3(iCol,iVer)+losscei(iCol,iVer))/SIN(dipMag(iCol,iVer))**2._r8  &
                           -subDiag(iVer)-superDiag(iVer)
              rHS(iVer) = tE(iCol,iVer) * wrk2(iCol,iVer)/ztodt + sourceg4(iCol,iVer)/SIN(dipMag(iCol,iVer))**2._r8 &
                          +(lossg3(iCol,iVer)*tN(iCol,iVer)+losscei(iCol,iVer)*ti(iCol,iVer))/SIN(dipMag(iCol,iVer))**2._r8
          enddo !iVer loop

          !-------------------------------------------------------------------------------------
          !  Calculate diagonal, superdiagonal, and right hand side upper boundary values
          !-------------------------------------------------------------------------------------
          superDiag(1)  = -Kei(2) / delZi(iCol,2) / delZ(iCol,1)
          diag(1) = wrk2(iCol,1)/ztodt - superDiag(1)
          rHS(1) = tE(iCol,1) * wrk2(iCol,1) / ztodt - FeUB(iCol) / delZ(iCol,1)

          !---------------------------------------------------------------------------------------------
          !  Calculate subdiagonal, diagonal, superdiagonal, and right hand side lower boundary values
          !---------------------------------------------------------------------------------------------
          subDiag(teTiBot) = -Kei(teTiBot) / delZi(iCol,teTiBot) / delZ(iCol,teTiBot)
          superDiag(teTiBot) = -Kei(teTiBotP) / delZi(iCol,teTiBotP) / delZ(iCol,teTiBot)
          diag(teTiBot) = wrk2(iCol,teTiBot)/ztodt + (lossg3(iCol,teTiBot)+losscei(iCol,teTiBot))/SIN(dipMag(iCol,teTiBot))**2._r8 &
                           -subDiag(teTiBot)-superDiag(teTiBot)
          rHS(teTiBot) = tE(iCol,teTiBot) * wrk2(iCol,teTiBot)/ztodt + sourceg4(iCol,teTiBot)/SIN(dipMag(iCol,teTiBot))**2._r8  &
                         +(lossg3(iCol,teTiBot)*tN(iCol,teTiBot)+losscei(iCol,teTiBot)*ti(iCol,teTiBot))/ &
			 SIN(dipMag(iCol,teTiBot))**2._r8 - superDiag(teTiBot) * tN(iCol,teTiBotP)

          !-------------------------------------------------
          ! Call solver to get electron temperature update
          !-------------------------------------------------
          call tridag(subDiag,diag,superDiag,rHS,tETemp,teTiBot)

          tE(iCol,1:teTiBot) = tETemp(1:teTiBot)
          do iVer = 1,teTiBot
             tE(iCol,iVer) = min(temax,tE(iCol,iVer))
             tE(iCol,iVer) = max(tN(iCol,iVer),tE(iCol,iVer))
          enddo
	  
          !---------------------------------------------------------------------------------------------------------
          !  Calculate ion temperature from electron temperature, ion-neutral and electron-ion loss terms, neutral 
          !  temperature, mass density and joule heating.  Set minimum value to neutral temperature and maximum 
          !  value to electron temperature for each column and vertical level
          !---------------------------------------------------------------------------------------------------------
          do iVer = 1,teTiBot
              ti(iCol,iVer) = (losscei(iCol,iVer) * tE(iCol,iVer) + losscin(iCol,iVer) * tN(iCol,iVer) +   &
                               rho(iCol,iVer) * qjoule(iCol,iVer))/(losscei(iCol,iVer) + losscin(iCol,iVer))
              ti(iCol,iVer) = max(tN(iCol,iVer),ti(iCol,iVer))
              ti(iCol,iVer) = min(tE(iCol,iVer),ti(iCol,iVer))
          enddo
      
          !--------------------------------------------------------------------------------------------------------
          ! Check for convergence which is a change of electron temperature ratio to previous loop for all levels
          ! and columns of less than 0.05K.  Had to modify this to do convergence check on each column since 
          ! checking all columns in a chunk gives different answers depending on number of tasks and tasks per node.
          !--------------------------------------------------------------------------------------------------------
          if (ALL(ABS(tE(iCol,1:teTiBot) / tEPrevI(iCol,1:teTiBot) - 1._r8) < 0.05_r8)) then 
	
	    colConv(iCol) = .true.

          endif

        endif ! Column not converged

      enddo ! iCol loop
      !--------------------------------------------------------------
      !  Check to see if all columns have converged and set flag
      !--------------------------------------------------------------
      if (ALL(colConv(1:ncol))) converged = .true.

    enddo ! End of iteration loop
      
!    write(iulog,*)'lchnk, Number of tempei iterations, converged flag ', lchnk, iter, converged

    !--------------------------------------------------------------------------------------------------------
    ! Calculate electron-neutral heating and electron-ion Coulomb heating.  Then update dry static energy.
    !--------------------------------------------------------------------------------------------------------
    do iVer = 1, teTiBot
       do iCol = 1, ncol
          sqrtTE(iVer) = SQRT(tE(iCol,iVer))
          lossceN2(iCol,iVer) = lossc4p(iCol,iVer) * (1._r8 - lossc5 * tE(iCol,iVer)) * tE(iCol,iVer)
          lossceO2(iCol,iVer) = lossc6p(iCol,iVer) * (1._r8 + lossc7 * sqrtTE(iVer)) * sqrtTE(iVer)
          lossceO1(iCol,iVer) = lossc8p(iCol,iVer) * (1._r8 + lossc9 * tE(iCol,iVer)) * sqrtTE(iVer)
       enddo
    enddo
    qen(1:ncol,1:teTiBot) = (lossceN2(1:ncol,1:teTiBot)+lossceO2(1:ncol,1:teTiBot)+lossceO1(1:ncol,1:teTiBot)) *  &
                           (tE(1:ncol,1:teTiBot)-tN(1:ncol,1:teTiBot)) / rho(1:ncol,1:teTiBot)
    qei(1:ncol,1:teTiBot) = losscei(1:ncol,1:teTiBot) * (tE(1:ncol,1:teTiBot)-ti(1:ncol,1:teTiBot)) / rho(1:ncol,1:teTiBot)

    dSETendOut(1:ncol,1:teTiBot) = (qei(1:ncol,1:teTiBot)+qen(1:ncol,1:teTiBot)) / sToQConv    

    return

  end subroutine ionos_teti

!===============================================================================
 
  subroutine ionos_vamb(state, vAmbOpTend, ztodt, pbuf, istate, tE, tI, mmrPOp, vAmbBot)

!------------------------------------------------------------------------------------------
! Routine to compute O+ ion transport in the vertical dimension through ambipolar diffusion 
! Goal is to get diagonals and right hand side for tridiagonal O+ solver.
!------------------------------------------------------------------------------------------

    use cam_history,         only : addfld, add_default, phys_decomp
    use phys_grid,           only : get_lat_p, get_lon_p, get_rlat_p,get_rlon_p
    use physconst,           only : gravit, avogad                   ! Gravity (m/s2) and Avogadro's number (molecules/kmole)
    use time_manager,        only : get_step_size                    ! Routine to get current time step and time step size
    use physconst,           only : rairv, mbarv                     ! Constituent dependent rair and mbar
    use efield,              only : ed1, ed2                         ! x and y components of the electric field
    use mo_apex,             only : bnorth, beast, bmag, alatm       ! Magnetic field components and magnitude (nT) and 
                                                                     !                       geomagnetic latitude (radians)
    use exbdrift,            only : map_mag2geo, cal_mlt             ! Methods to convert e field from mag to geo coords
    use constituents,        only : cnst_get_ind, cnst_mw            ! Routines to get molecular weights for constituents
    use short_lived_species, only : slvd_index                       ! Routine to access short lived species in pbuf
    use charge_neutrality,   only : charge_fix
    use constituents,        only : pcnst                            ! Number of state%q constituents

    implicit none

!------------------------------Arguments--------------------------------

    type(physics_buffer_desc), pointer           :: pbuf(:)             ! physics buffer
    type(physics_state),   intent(in), target    :: state               ! physics state structure
    type(ionos_state),     intent(in), target    :: istate              ! ionosphere state structure

    real(r8), intent(in) :: ztodt   ! Physics time step

    real(r8), dimension(pcols,pver),   intent(out)    :: vAmbOpTend   ! O+ tendency

    real(r8), dimension(:,:), pointer, intent(in)     :: tE           ! Pointer to electron temperature in pbuf (K) 
    real(r8), dimension(:,:), pointer, intent(in)     :: tI           ! Pointer to ion temperature in pbuf (K) 

    real(r8), dimension(:,:), pointer, intent(inout)  :: mmrPOp       ! O plus ion mass mixing ratio kg/kg

    integer, intent(in)  :: vAmbBot   ! bottom of O+ transport calculations 

!---------------------------Local storage-------------------------------
    integer :: lchnk                                    ! Chunk number 
    integer :: ncol                                     ! Number of atmospheric columns
    integer :: vAmbBotP                                 ! bottom of O+ transport calculations plus one more level 

    integer :: iVer                                     ! Counter for vertical loops
    integer :: iCol                                     ! Counter for column loops

    integer :: indxOp,sIndxOp                           ! state%q or pbuf index for O+ mixing ratio

    real(r8) :: cBF                                     ! Burnside factor
    real(r8) :: rMassOp                                 ! Molecular mass of O+ (kg/kmole)
    real(r8) :: tR                                      ! Reduced temperature (K)
    real(r8) :: dAPrime                                 ! Factor needed for dA (cm*cm/s/K
    real(r8) :: delTpDelZ                               ! Partial derivative of Tp with z
    real(r8) :: delTpDelZ1                              ! Top level of partial derivative of Tp with z
    real(r8) :: pScaleHeight                            ! Pressure scale height
    real(r8) :: bMagT                                   ! Magnetic field magnitude in Teslas (tesla)

    real(r8) :: c31,c41,c43r,c5,c6                      ! Coefficients for upper boundary O+ calculations
    real(r8) :: aOpFluxCoeff                            ! Coefficient for O+ flux calculation
    real(r8) :: oPFlux                                  ! O+ flux (?)
    real(r8) :: ndensOpTop                              ! Top (zeroth) level O+ flux (?)
         
    real(r8), dimension(pcols)          :: ed1_geo      ! zonal electric field on geographic grid [V/m]
    real(r8), dimension(pcols)		:: ed2_geo      ! meridional electric field on geographic grid [V/m]
    real(r8), dimension(pcols)		:: epot_geo     ! electric potential on geographic grid
    real(r8), dimension(pcols)          :: mlt  	! mag.local time of WACCM geographic grid points
 
    real(r8), dimension(pcols)          :: decMag	! magnetic declination angle (radians)
    real(r8), dimension(pcols)          :: sinDec	! sine of the magnetic declination angle
    real(r8), dimension(pcols)          :: cosDec	! cosine of the magnetic declination angle

    real(r8), dimension(pcols)          :: geoMagLat    ! Geomagnetic latitude (radians)

    real(r8), dimension(pcols)          :: eHPB         ! Horizontal electric field perpendicular to magnetic field [V/m]
    
    real(r8), dimension(pcols,pver)  	:: tP		! Plasma temperature (K)

    real(r8), dimension(:,:), pointer	:: pMid	        ! Midpoint pressure (Pa)
    real(r8), dimension(:,:), pointer	:: tN  	        ! Neutral temperature (K)
    real(r8), dimension(:,:), pointer	:: uN     	! Neutral zonal wind velocity (m/s)
    real(r8), dimension(:,:), pointer	:: vN  	        ! Neutral meridional wind velocity (m/s)

    real(r8), dimension(:,:), pointer	:: pInt 	! Interface pressure (Pa)
    real(r8), dimension(:,:), pointer	:: tNInt	! Interface Neutral Temperture (K)
    real(r8), dimension(pcols,pverp) 	:: tIInt	! Interface Ion Temperture (K)
    real(r8), dimension(pcols,pverp) 	:: tEInt	! Interface Electron Temperture (K)
    real(r8), dimension(:,:), pointer 	:: rairvi	! Constituent dependent gas constant on interface levels

    real(r8), dimension(:,:), pointer  	:: omega 	! Vertical velocity (Pa/s?)
    real(r8), dimension(pcols,pver)  	:: wN		! Neutral vertical wind velocity (m/s)

    real(r8), dimension(:,:), pointer	:: ndensN2      ! N2 number density (cm-3)
    real(r8), dimension(:,:), pointer   :: ndensO2	! O2 number density (cm-3)
    real(r8), dimension(:,:), pointer	:: ndensO1	! O number density (cm-3)
    real(r8), dimension(:,:), pointer	:: ndensOp	! O+ number density (cm-3)

    real(r8), dimension(:,:), pointer  	:: dipMag	! dip angle for each column (radians)
    real(r8), dimension(pcols,pver)  	:: dipMagD	! dip angle for each column (degrees)

    real(r8), dimension(:), pointer     :: zenAngD	! zenith angle (degrees)

    real(r8), dimension(pcols)       	:: FeUB 	! electron heat flux at upper boundary

    real(r8), dimension(pcols,pver)  	:: delZ 	! Delta z: midpoints (cm)
    real(r8), dimension(pcols,pverp) 	:: delZi	! Delta z: interfaces (cm)

    real(r8), dimension(pcols,pver)     :: dA           ! Ambipolar diffusion coefficient for c1 and c2 coefficients (cm*cm/s)
    real(r8), dimension(pcols,pver)  	:: c1    	! c1 coefficient for matrix diagonals: midpoints
    real(r8), dimension(pcols,pver)  	:: c2    	! c1 coefficient for matrix diagonals: midpoints
    real(r8), dimension(pcols,pverp) 	:: c1i   	! c2 coefficient for matrix diagonals: interfaces
    real(r8), dimension(pcols,pverp) 	:: c2i   	! c2 coefficient for matrix diagonals: interfaces
 
    real(r8), dimension(vAmbBot)      	:: subdiag	! subdiagonal values for Te tendency solving
    real(r8), dimension(vAmbBot)      	:: superdiag	! superdiagonal values for Te tendency solving
    real(r8), dimension(vAmbBot)        :: diag 	! diagonal values for Te tendency solving
    real(r8), dimension(vAmbBot)        :: rHS  	! RHS of electron temperature update
    real(r8), dimension(vAmbBot)        :: newDensOp	! Op result from tridag solver
    real(r8), dimension(vAmbBot)        :: newMMROp	! O+ mass mixing ratio converted from O+ number density from solver
    real(r8), dimension(pcols,pver)     :: mmrSOp       ! O+ mass mixing ratio converted for all columns
    real(r8), dimension(pcols,pver)     :: opmmr_tmp
    real(r8), dimension(pcols,pver)     :: emmr_tmp

    logical :: lq(pcnst)         

!-----------------------------------------------------------------------------------------------------------------
      
!    write(iulog,*)'lchnk, Inside vamb '

    !--------------------------------------------------------------------------------------------------------- 
    !  Initialize arrays to zero
    !---------------------------------------------------------------------------------------------------------
    ed1_geo(:)          = 0._r8
    ed2_geo(:)          = 0._r8
    epot_geo(:)         = 0._r8
    mlt(:)              = 0._r8
    decMag(:)           = 0._r8
    sinDec(:)           = 0._r8
    cosDec(:)           = 0._r8
    geoMagLat(:)        = 0._r8
    eHPB(:)             = 0._r8
    tP(:,:)             = 0._r8
    tIInt(:,:)          = 0._r8
    tEInt(:,:)          = 0._r8
    wN(:,:)             = 0._r8
    dipMagD(:,:)        = 0._r8
    FeUB(:)             = 0._r8 
    delZ(:,:)           = 0._r8 
    delZi(:,:)          = 0._r8
    dA(:,:)             = 0._r8 
    c1(:,:)             = 0._r8 
    c2(:,:)             = 0._r8 
    c1i(:,:)            = 0._r8 
    c2i(:,:)            = 0._r8 
    subDiag(:)          = 0._r8         
    superDiag(:)        = 0._r8        
    diag(:)             = 0._r8 
    rHS(:)              = 0._r8 
    newDensOp(:)        = 0._r8 
    newMMROp(:)         = 0._r8 
    mmrSOp(:,:)         = 0._r8 
    vAmbOpTend(:,:)     = 0._r8

    !--------------------------------------------------------------------------------------
    !  Get lchnk and ncol from 
    !--------------------------------------------------------------------------------------
    lchnk = state%lchnk
    ncol = state%ncol

    !-------------------------------------------------------------------------------
    !  Set Burnside Factor and bottom level plus one to calculate O+ transport
    !-------------------------------------------------------------------------------
!    cBF = 1._r8
    cBF = 0.8_r8
    vAmbBotP = vAmbBot + 1
    
    !----------------------------------------------------------------
    !  Get variables needed from the ionosphere state structure
    !----------------------------------------------------------------
    dipMag   => istate%dipMag(1:ncol,1:pver)
     
    ndensO1  => istate%ndensO1(1:ncol,1:pver)
    ndensN2  => istate%ndensN2(1:ncol,1:pver)
    ndensO2  => istate%ndensO2(1:ncol,1:pver)
    ndensOp  => istate%ndensOp(1:ncol,1:pver)
    
    tNInt    => istate%tNInt(1:ncol,1:pverp)
    rairvi   => istate%rairvi(1:ncol,1:pverp)

    zenAngD  => istate%zenAngD(1:ncol) 

    !----------------------------------------------------------------
    !  Get variables needed from the physics state structure
    !----------------------------------------------------------------
    pInt  => state%pint(1:ncol,1:pverp)
    pMid  => state%pmid(1:ncol,1:pver)
    omega => state%omega(1:ncol,1:pver)
    tN    => state%t(1:ncol,1:pver)
    uN    => state%u(1:ncol,1:pver)
    vN    => state%v(1:ncol,1:pver)

    !----------------------------------------------------------------------------------------------------------------------
    !  If this is the first time calculated then initialize to neutral temperature. is_first_step means this is the first
    !  time step of an initial run. 
    !----------------------------------------------------------------------------------------------------------------------
    if( is_first_step() ) then   
      
      ti(1:ncol,1:pver) = tN(1:ncol,1:pver)
      tE(1:ncol,1:pver) = tN(1:ncol,1:pver)

    endif

    call cnst_get_ind( 'Op',  indxOp, abort=.false. )
    if (indxOp > 0) then
      rMassOp = cnst_mw(indxOP)
    else
      sIndxOp = get_spc_ndx( 'Op' )
      if (sIndxOp > 0) then
	rMassOp  = adv_mass(sIndxOp)
      else
        call endrun('ionosphere: Cannot find short-lived index for Op in ionos_vamb')         
      endif
    endif

    !-------------------------------------------------
    !  Get geomagnetic latitudes for this chunk
    !-------------------------------------------------
    geoMagLat(1:ncol) = alatm(1:ncol,lchnk)

    !--------------------------------------------------
    !  Calculate cosine/sine of declination angle
    !--------------------------------------------------
    decMag(1:ncol) = -atan2( beast(1:ncol,lchnk),bnorth(1:ncol,lchnk) )
    cosDec(1:ncol) = cos( decMag(1:ncol) )
    sinDec(1:ncol) = sin( decMag(1:ncol) )

    !------------------------------------------------------------------------------------------------------------
    !  Calculate the magnetic local time of WACCM geographic grid points needed for mag to geo conversion
    !------------------------------------------------------------------------------------------------------------
    call cal_mlt( mlt, lchnk, ncol )
    
    !------------------------------------------------------------------------------------------------------------------
    !  Convert e_field module zonal and meridional electric field components to geographic coords using exbdrift module
    !  map_mag2geo method
    !------------------------------------------------------------------------------------------------------------------ 
    call map_mag2geo( mlt, lchnk, ncol, ed1_geo, ed2_geo, epot_geo )       
    
    !-------------------------------------------------------------------------------------
    !  Calculate ion temperature on interface levels
    !-------------------------------------------------------------------------------------   
    do iVer = 2, vAmbBotP
 
      do iCol = 1, ncol

        tIInt(iCol,iVer) = 0.5_r8 * tI(iCol,iVer) + 0.5_r8 * tI(iCol,iVer-1)

      enddo
      
    enddo
    
    !--------------------------
    !  Get top interface level
    !--------------------------
    do iCol = 1, ncol
    
      tIInt(iCol,1) = 1.5_r8 * tIInt(iCol,2) - 0.5_r8 * tIInt(iCol,3)
       
    enddo
    
    !--------------------------
    !  Get bottom interface level
    !--------------------------
    do iCol = 1, ncol
    
      tIInt(iCol,pverp) = 1.5_r8 * tIInt(iCol,pver) - 0.5_r8 * tIInt(iCol,pver-1) 
      
    enddo
    
    !-------------------------------------------------------------------------------------
    !  Calculate electron temperature on interface levels
    !-------------------------------------------------------------------------------------   
    do iVer = 2, vAmbBotP
 
      do iCol = 1, ncol

        tEInt(iCol,iVer) = 0.5_r8 * tE(iCol,iVer) + 0.5_r8 * tE(iCol,iVer-1)

      enddo
      
    enddo

    !--------------------------
    !  Get top interface level
    !--------------------------
    do iCol = 1, ncol
    
      tEInt(iCol,1) = 1.5_r8 * tEInt(iCol,2) - 0.5_r8 * tEInt(iCol,3) 
      
    enddo
    
    !--------------------------------------------------------------------------------------------------------
    !  Get plasma temperature needed for delTp/delz
    !--------------------------------------------------------------------------------------------------------
    do iVer = 1, vAmbBot

      do iCol = 1, ncol

	tP(iCol,iVer) = ( tI(iCol,iVer) + tE(iCol,iVer) ) / 2._r8

      enddo
      
    enddo
    
    !--------------------------
    !  Get bottom interface level
    !--------------------------
    do iCol = 1, ncol
    
      tP(iCol,vAmbBotP) = 1.5_r8 * tP(iCol,vAmbBot) - 0.5_r8 * tP(iCol,vAmbBot-1)
       
    enddo

    !------------------------------------------------------------------------------------------
    !  Calculate thickness terms for vertical derivative
    !------------------------------------------------------------------------------------------
    do iVer = 1, vAmbBot

      do iCol = 1, ncol

        delZ(iCol,iVer) = (pInt(iCol,iVer+1) - pInt(iCol,iVer)) * rairv(iCol,iVer,lchnk) * tN(iCol,iVer) / pMid(iCol,iVer) / gravit

      enddo

    enddo

    do iVer = 2, vAmbBotP	    ! Assuming vAmbBotP < pverp

      do iCol = 1, ncol

        delZi(iCol,iVer) = (pMid(iCol,iVer) - pMid(iCol,iVer-1)) * rairvi(iCol,iVer) * &
    							      tNInt(iCol,iVer) / pInt(iCol,iVer) / gravit
      enddo

    enddo

    do iCol = 1, ncol

      delZi(iCol,1) = 1.5_r8 * delZi(iCol,2) - .5_r8 * delZi(iCol,3)

      !----------------------------------------------------------
      !  Convert delZ variables from meters to centimeters
      !----------------------------------------------------------
      delZi(iCol,1:vAmbBotP) = delZi(iCol,1:vAmbBotP) * 100._r8
      delZ(iCol,1:vAmbBot) = delZ(iCol,1:vAmbBot) * 100._r8

    enddo
    
    !-----------------------------------------------------------------------------
    !  Compute the horizontal electric field perpendicular to the magnetic field
    !-----------------------------------------------------------------------------
    eHPB(1:ncol) = ed1_geo(1:ncol) * cosDec(1:ncol) - ed2_geo(1:ncol) * sinDec(1:ncol) 
    
    !------------------------------------------------------------------
    ! Calculate pressure scale height and neutral vertical velocity
    !------------------------------------------------------------------
    do iVer = 2,vAmbBot
    
      do iCol = 1,ncol
      
        pScaleHeight = .5_r8*(rairv(iCol,iVer,lchnk)*tN(iCol,iVer)+rairv(iCol,iVer-1,lchnk)*tN(iCol,iVer-1))/gravit
      
        wN(iCol,iVer) = -omega(iCol,iVer) / pMid(iCol,iVer) * pScaleHeight    

     enddo
       
    enddo

    !--------------------------
    !  Get top midpoint level
    !--------------------------
    wN(1:ncol,1) = 1.5_r8 * wN(1:ncol,2) - 0.5_r8 * wN(1:ncol,3) 
	    
    !--------------------------------------------------------------------------------------------------------
    !  Loop over vertical levels then columns and call tridiagonal solver for each column to get O+
    !--------------------------------------------------------------------------------------------------------
    do iVer = 1, vAmbBot

      do iCol = 1, ncol

        !-------------------------------------------------------------------------------
        !  Calculate Da ambipolar diffusion coefficient for c1 and c2 coefficients (cgs)
        !-------------------------------------------------------------------------------
        tR = ( tI(iCol,iVer) + tN(iCol,iVer) ) *.5_r8
	
        dAPrime = 1.42E+17_r8 / (SQRT(tR) * (1._r8 - 0.064_r8 * LOG10(tR))**2 * ndensO1(iCol,iVer) * cBF + 18.6_r8 * &
								    ndensN2(iCol,iVer) + 18.1_r8 * ndensO2(iCol,iVer))

	dA(iCol,iVer) = MIN(2.0e11_r8, dAPrime * 2._r8 * tP(iCol,iVer))
	
        !-------------------------------------------------------------------------------
        !  Calculate derivative for c1 and c2 coefficients (cgs)
        !-------------------------------------------------------------------------------
	delTpDelZ = ( tP(iCol,iVer) - tP(iCol,iVer+1) ) / delZ(iCol,iVer)
	
	!-------------------------------------------------------------
	!  Need to save top value for upper boundary calculations
	!-------------------------------------------------------------
	if (iVer == 1) delTpDelZ1 = delTpDelZ	 

        !-----------------------------------------------------------------------------
        !  Convert magnetic field magnitude units from nanoTesla to Tesla
        !-----------------------------------------------------------------------------
        bMagT = bmag(iCol,lchnk)*1.e-9_r8
       
        !---------------------------------------------------------------------------------
        !  Compute C1 coefficient for matrix formation and RHS
	!  eHPB(V/m) bMagT(nT=V*s/m2) so eHPB/bMagT=(V/m/V/s*m2=m/s=.01*cm/s)
	!  Factors .01 and 100. is to convert from SI to cgs units (meters to centimeters)
        !---------------------------------------------------------------------------------
        c1(iCol,iVer) = istate%we2d(iCol,iVer) * 100._r8 + &
	                ( ( vN(iCol,iVer) * cos(decMag(iCol)) + uN(iCol,iVer) * sin(decMag(iCol)) ) * &
			cos(dipMag(iCol,iVer)) - wN(iCol,iVer) * sin(dipMag(iCol,iVer)) ) * &
			sin(dipMag(iCol,iVer)) * 100._r8 + &
			dA(iCol,iVer) * ( delTpDelZ  / tP(iCol,iVer) + &
			0.5_r8 * rMassOp / avogad * gravit / kboltz / tP(iCol,iVer) * .01_r8 ) * &
			sin(dipMag(iCol,iVer))**2
	          
        !--------------------------------------------------------------
        !  Compute C2 coefficient for matrix formation 
	!  dA(cm2/s)
        !--------------------------------------------------------------
        c2(iCol,iVer) = dA(iCol,iVer) * sin(dipMag(iCol,iVer))**2

      enddo
    
    enddo    

    !----------------------------------------------------
    !  Get c1 and c2 coefficients on interface levels
    !----------------------------------------------------
    do iVer = 2, vAmbBot	    ! Assuming vAmbBot < pverp

      do iCol = 1, ncol

        c1i(iCol,iVer) = 0.5_r8 * c1(iCol,iVer) + 0.5_r8 * c1(iCol,iVer-1)
        c2i(iCol,iVer) = 0.5_r8 * c2(iCol,iVer) + 0.5_r8 * c2(iCol,iVer-1)

      enddo

    enddo

    !-------------------------------------------------------
    !  Get top interface level for c1i and c2i coefficients
    !-------------------------------------------------------
    do iCol = 1, ncol
  
      c1i(iCol,1) = 1.5_r8 * c1i(iCol,2) - 0.5_r8 * c1i(iCol,3) 
      c2i(iCol,1) = 1.5_r8 * c2i(iCol,2) - 0.5_r8 * c2i(iCol,3) 
    
    enddo

    !---------------------------------------------------------
    !  Get bottom interface level for c1i and c2i coefficients
    !---------------------------------------------------------
    do iCol = 1, ncol
    
      c1i(iCol,vAmbBotP) = 1.5_r8 * c1i(iCol,vAmbBot) - 0.5_r8 * c1i(iCol,vAmbBot-1) 
      c2i(iCol,vAmbBotP) = 1.5_r8 * c2i(iCol,vAmbBot) - 0.5_r8 * c2i(iCol,vAmbBot-1) 
    
    enddo
   
    !--------------------------------------------------------------------------------------------------------
    !  Loop over columns then vertical levels and call tridiagonal solver for each column to get O+
    !--------------------------------------------------------------------------------------------------------	 
    do iCol = 1, ncol
      
      do iVer = 2, vAmbBot-1
 
        !-----------------------------------------------------------------------------------------------------
        !  Calculate sub-diagonal, diagonal, superdiagonal, and right hand side
        !-----------------------------------------------------------------------------------------------------
        subDiag(iVer) = -( 0.25_r8 * c1i(iCol,iVer) / delZ(iCol,iVer) + 0.5_r8 *c2i(iCol,iVer) / &
	                                                                  delZi(iCol,iVer) / delZ(iCol,iVer) )
		     
        diag(iVer) =  1._r8 / ztodt - 0.25_r8 * ( c1i(iCol,iVer) - c1i(iCol,iVer+1) ) / delZ(iCol,iVer) + 0.5_r8 * &
	             ( c2i(iCol,iVer) / delZi(iCol,iVer) /delZ(iCol,iVer) + c2i(iCol,iVer+1) / delZ(iCol,iVer) / &
		     delZi(iCol,iVer+1) ) 

        superDiag(iVer) =  0.25_r8 * c1i(iCol,iVer+1) / delZ(iCol,iVer) - 0.5_r8 * c2i(iCol,iVer+1) / &
	                                                                delZ(iCol,iVer) / delZi(iCol,iVer+1)
	
        rHS(iVer) = -subDiag(iVer) * ndensOp(iCol,iVer-1) + &
	            ( 1._r8 / ztodt + 0.25_r8 * ( c1i(iCol,iVer) - c1i(iCol,iVer+1) ) / delZ(iCol,iVer) - &
		    0.5_r8 * ( c2i(iCol,iVer) / delZi(iCol,iVer) /delZ(iCol,iVer) + &
		    c2i(iCol,iVer+1) / delZ(iCol,iVer) / delZi(iCol,iVer+1) ) ) * ndensOp(iCol,iVer) -  &
		    superDiag(iVer) * ndensOp(iCol,iVer+1)
				  
      enddo !iVer loop
     
      !-----------------------------------------------------------------
      !  Get upper boundary conditions for diagonal and RHS in cgs units
      !-----------------------------------------------------------------
      
      ! Can we make the flux a constant array, to save some computing each time we go through the iCol loop?  ?-HL

      !-----------------------------------------------------------------
      !  Calculate O+ flux 
      !-----------------------------------------------------------------      
      if (ABS(geoMagLat(iCol)) > pi/12._r8) aOPFluxCoeff = 1._r8
      if (ABS(geoMagLat(iCol)) <= pi/12._r8) aOPFluxCoeff = 0.5_r8 * ( 1._r8 + sin((ABS(geoMagLat(iCol)) - pi / 24._r8) * 12._r8) )
  
      if (zenAngD(iCol) <= 80._r8) oPFlux = dayOPFlux * aOPFluxCoeff
      if (zenAngD(iCol) > 80._r8 .and. zenAngD(iCol) < 100._r8) &
             oPFlux = ( 0.5_r8 * (dayOPFlux + nightOPFlux) + 0.5_r8 *(dayOPFlux - nightOPFlux) * & 
	     cos(pi * (zenAngD(iCol) - 80._r8) / 20._r8) ) * aOPFluxCoeff
      if (zenAngD(iCol) >= 100._r8) oPFlux = nightOPFlux * aOPFluxCoeff
     
      !-----------------------------------------------------------------
      !  Calculate coefficients needed for diagonal and RHS top values
      !-----------------------------------------------------------------      
      c31 = -dA(iCol,1) * sin(dipMag(iCol,1))**2
      c43r = ( 2._r8 * delTpDelZ1 + rMassOp / avogad * gravit / kboltz * .01_r8 )/2._r8/tP(iCol,1)
      c41 = c31 * c43r
      
      if (dipMag(iCol,1) == 0._r8) c5 = 0._r8
      if (dipMag(iCol,1) /= 0._r8) c5 = 2._r8 * delZ(iCol,1) / c31 * oPFlux
      
      if (dipMag(iCol,1) == 0._r8) c6 = 0._r8
      if (dipMag(iCol,1) /= 0._r8) c6 = - 2._r8 * delZ(iCol,1) * c43r
      
      !-----------------------------------------------------------------
      !  Calculate diagonals and RHS top level values
      !-----------------------------------------------------------------
      subDiag(1) = -( 0.25_r8 * c1i(iCol,1) / delZ(iCol,1) + 0.5_r8 *c2i(iCol,1) / &
	                                                                  delZi(iCol,1) / delZ(iCol,1) )
		     
      diag(1) =  1._r8 / ztodt - 0.25_r8 * ( c1i(iCol,1) - c1i(iCol,2) ) / delZ(iCol,1) + 0.5_r8 * &
	             ( c2i(iCol,1) / delZi(iCol,1) /delZ(iCol,1) + c2i(iCol,2) / delZ(iCol,1) / &
		     delZi(iCol,2) )

      superDiag(1) =  0.25_r8 * c1i(iCol,2) / delZ(iCol,1) - 0.5_r8 * c2i(iCol,2) / &
	                                                                delZ(iCol,1) / delZi(iCol,2)
      ndensOpTop = c5 + c6 * ndensOP(iCol,1) + ndensOP(iCol,2)

      rHS(1) = - subDiag(1) * ndensOpTop + &
	            ( 1._r8 / ztodt + 0.25_r8 * ( c1i(iCol,1) - c1i(iCol,2) ) / delZ(iCol,1) - &
		    0.5_r8 * ( c2i(iCol,1) / delZi(iCol,1) /delZ(iCol,1) + &
		    c2i(iCol,2) / delZ(iCol,1) / delZi(iCol,2) ) ) * ndensOp(iCol,1) - &
                    superDiag(1) * ndensOp(iCol,2)

		    
      !--------------------------------------------------------------
      !  Set upper boundary conditions for diagonals and RHS
      !--------------------------------------------------------------      
      diag(1) = subDiag(1) * c6 + diag(1)
      
      superDiag(1) = superDiag(1) + subDiag(1)
      
      rHS(1) = rHS(1) - subDiag(1) * c5

      !----------------------------------------------
      !  Calculate diagonal and RHS bottom values
      !----------------------------------------------
      subDiag(vAmbBot) = -( 0.25_r8 * c1i(iCol,vAmbBot) / delZ(iCol,vAmbBot) + 0.5_r8 *c2i(iCol,vAmbBot) / &
	                                                                  delZi(iCol,vAmbBot) / delZ(iCol,vAmbBot) )
		     
      diag(vAmbBot) =  1._r8 / ztodt - 0.25_r8 * ( c1i(iCol,vAmbBot) - c1i(iCol,vAmbBotP) ) / delZ(iCol,vAmbBot) + 0.5_r8 * &
	             ( c2i(iCol,vAmbBot) / delZi(iCol,vAmbBot) /delZ(iCol,vAmbBot) + c2i(iCol,vAmbBotP) / delZ(iCol,vAmbBot) / &
		     delZi(iCol,vAmbBotP) ) 

      superDiag(vAmbBot) = 0.25_r8 * c1i(iCol,vAmbBotP) / delZ(iCol,vAmbBot) - 0.5_r8 * c2i(iCol,vAmbBotP) / &
	                                                                delZ(iCol,vAmbBot) / delZi(iCol,vAmbBotP)

      rHS(vAmbBot) = -subDiag(vAmbBot) * ndensOp(iCol,vAmbBot-1) + & 
	            ( 1._r8 / ztodt + 0.25_r8 * ( c1i(iCol,vAmbBot) - c1i(iCol,vAmbBotP) ) / delZ(iCol,vAmbBot) - &
		    0.5_r8 * ( c2i(iCol,vAmbBot) / delZi(iCol,vAmbBot) /delZ(iCol,vAmbBot) + &
		    c2i(iCol,vAmbBotP) / delZ(iCol,vAmbBot) / delZi(iCol,vAmbBotP) ) ) * ndensOp(iCol,vAmbBot) - &
                    superDiag(vAmbBot) * ndensOp(iCol,vAmbBotP)

      !--------------------------------------------------------------
      !  Set lower boundary condition for RHS
      !--------------------------------------------------------------      

      rHS(vAmbBot) = rHS(vAmbBot) - superDiag(vAmbBot) * ndensOp(iCol,vAmbBotP)
         
      !-------------------------------------------------
      ! Call solver to get O+ number density update
      !-------------------------------------------------
      call tridag(subDiag,diag,superDiag,rHS,newDensOp,vAmbBot)

      !-----------------------------------------------------------------------------
      !!! Fix to avoid very small Op number density around the F peak
      !-----------------------------------------------------------------------------
      do iVer = 1,vAmbBot
         newDensOp(iVer) = MAX(100._r8,newDensOp(iVer))
      enddo

     !------------------------------------------------------------------------------------------------------------------------------
     !  Convert new O+ number density to mixing ratio and output to pbuf or vambtend%q depending on whether it is short lived or not  
     !------------------------------------------------------------------------------------------------------------------------------      
      newMMROp(1:vAmbBot) = newDensOp(1:vAmbBot) / mbarv(iCol,1:vAmbBot,lchnk) * rMassOp / pMid(iCol,1:vAmbBot) * &
                                                                                (kboltz * tN(iCol,1:vAmbBot)) / 1.E-06_r8

      if (indxOp > 0) then
          mmrSOp(iCol,1:vAmbBot) = newMMROp(1:vAmbBot)
      else
          mmrPOp(iCol,1:vAmbBot) = newMMROp(1:vAmbBot)
      endif      
 
    enddo ! iCol loop   

    !-----------------------------
    !  Calculate O+ tendency
    !-----------------------------    
    if (indxOp > 0) vAmbOpTend(1:ncol,1:vAmbBot) = ( ( mmrSOp(1:ncol,1:vAmbBot) - state%q(1:ncol,1:vAmbBot,indxOp) ) / ztodt )

    return

  end subroutine ionos_vamb

!===============================================================================

!-----------------------------------------------------------------------
! Simple tridiagonal solver routine
!-----------------------------------------------------------------------

  SUBROUTINE tridag(a,b,c,r,u,n)

    INTEGER,INTENT(IN)      :: n
    REAL(r8),INTENT(IN)     :: a(n),b(n),c(n),r(n)
    REAL(r8),INTENT(INOUT)  :: u(n)
    !------------------------------
    !  Local variables
    !------------------------------
    INTEGER j
    REAL(r8) :: bet,gam(n)

    if(b(1).eq.0._r8) call endrun('ionosphere: bt(1)=0 in tridag')
    bet=b(1)
    u(1)=r(1)/bet
    do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if(bet.eq.0._r8) call endrun('ionosphere: bet=0 in tridag')
      u(j)=(r(j)-a(j)*u(j-1))/bet
    end do

    do j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do

    return

  END SUBROUTINE tridag

end module ionosphere
