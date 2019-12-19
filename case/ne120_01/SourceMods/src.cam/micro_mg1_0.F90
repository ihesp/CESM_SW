module micro_mg1_0

!---------------------------------------------------------------------------------
! Purpose:
!   MG microphysics
!
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!
! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
! microphysics in cooperation with the MG liquid microphysics. This is
! controlled by the do_cldice variable.
!
! NOTE: If do_cldice is false, then MG microphysics should not update CLDICE
! or NUMICE; however, it is assumed that the other microphysics scheme will have
! updated CLDICE and NUMICE. The other microphysics should handle the following
! processes that would have been done by MG:
!   - Detrainment (liquid and ice)
!   - Homogeneous ice nucleation
!   - Heterogeneous ice nucleation
!   - Bergeron process
!   - Melting of ice
!   - Freezing of cloud drops
!   - Autoconversion (ice -> snow)
!   - Growth/Sublimation of ice
!   - Sedimentation of ice
!---------------------------------------------------------------------------------
! modification for sub-columns, HM, (orig 8/11/10)
! This is done using the logical 'microp_uniform' set to .true. = uniform for subcolumns
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure to specific humidity formula
! 3) svp over water
! 4) svp over ice

#ifndef HAVE_GAMMA_INTRINSICS
use shr_spfn_mod, only: gamma => shr_spfn_gamma
#endif

  use wv_sat_methods, only: &
       svp_water => wv_sat_svp_water, &
       svp_ice => wv_sat_svp_ice, &
       svp_to_qsat => wv_sat_svp_to_qsat, &
       default_idx

  use phys_control, only: phys_getopts
  use physconst, only: svp_epsilo  => epsilo, & 
                       svp_latvap  => latvap, & 
                       svp_latice  => latice, & 
                       svp_rh2o    => rh2o,   & 
                       svp_cpair   => cpair,  & 
                       svp_tmelt   => tmelt,  & 
                       svp_h2otrip => h2otrip
  use perf_mod

implicit none
private
save

! Note: The liu_in option has been removed, as there was a serious bug with this
! option being set to false. The code now behaves as if the default liu_in=.true.
! is always on. Addition/reinstatement of ice nucleation options will likely be
! done outside of this module.

public :: &
     micro_mg_init, &
     micro_mg_get_cols, &
     micro_mg_tend

integer, parameter :: r8 = selected_real_kind(12)      ! 8 byte real

real(r8) :: g              !gravity
real(r8) :: r              !Dry air Gas constant
real(r8) :: rv             !water vapor gas contstant
real(r8) :: cpp            !specific heat of dry air
real(r8) :: rhow           !density of liquid water
real(r8) :: tmelt          ! Freezing point of water (K)
real(r8) :: xxlv           ! latent heat of vaporization
real(r8) :: xlf            !latent heat of freezing
real(r8) :: xxls           !latent heat of sublimation

real(r8) :: rhosn  ! bulk density snow
real(r8) :: rhoi   ! bulk density ice

real(r8) :: ac,bc,as,bs,ai,bi,ar,br  !fall speed parameters 
real(r8) :: ci,di    !ice mass-diameter relation parameters
real(r8) :: cs,ds    !snow mass-diameter relation parameters
real(r8) :: cr,dr    !drop mass-diameter relation parameters
real(r8) :: f1s,f2s  !ventilation param for snow
real(r8) :: Eii      !collection efficiency aggregation of ice
real(r8) :: Ecr      !collection efficiency cloud droplets/rain
real(r8) :: f1r,f2r  !ventilation param for rain
real(r8) :: DCS      !autoconversion size threshold
real(r8) :: qsmall   !min mixing ratio 
real(r8) :: bimm,aimm !immersion freezing
real(r8) :: rhosu     !typical 850mn air density
real(r8) :: mi0       ! new crystal mass
real(r8) :: rin       ! radius of contact nuclei
real(r8) :: pi       ! pi

! Additional constants to help speed up code

real(r8) :: cons1
real(r8) :: cons4
real(r8) :: cons5
real(r8) :: cons6
real(r8) :: cons7
real(r8) :: cons8
real(r8) :: cons11
real(r8) :: cons13
real(r8) :: cons14
real(r8) :: cons16
real(r8) :: cons17
real(r8) :: cons22
real(r8) :: cons23
real(r8) :: cons24
real(r8) :: cons25
real(r8) :: cons27
real(r8) :: cons28

real(r8) :: lammini
real(r8) :: lammaxi
real(r8) :: lamminr
real(r8) :: lammaxr
real(r8) :: lammins
real(r8) :: lammaxs

! parameters for snow/rain fraction for convective clouds
real(r8) :: tmax_fsnow ! max temperature for transition to convective snow
real(r8) :: tmin_fsnow ! min temperature for transition to convective snow

!needed for findsp
real(r8) :: tt0       ! Freezing temperature

real(r8) :: csmin,csmax,minrefl,mindbz

real(r8) :: rhmini     ! Minimum rh for ice cloud fraction > 0.

logical :: use_hetfrz_classnuc ! option to use heterogeneous freezing

character(len=16)  :: micro_mg_precip_frac_method  ! type of precipitation fraction method
real(r8)           :: micro_mg_berg_eff_factor     ! berg efficiency factor


!===============================================================================
contains
!===============================================================================

subroutine micro_mg_init( &
     kind, gravit, rair, rh2o, cpair,  &
     rhoh2o, tmelt_in, latvap, latice, &
     rhmini_in, micro_mg_dcs, use_hetfrz_classnuc_in, &
     micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, errstring)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! initialize constants for the morrison microphysics
! 
! Author: Andrew Gettelman Dec 2005
! 
!-----------------------------------------------------------------------

integer,          intent(in)  :: kind            ! Kind used for reals
real(r8),         intent(in)  :: gravit
real(r8),         intent(in)  :: rair
real(r8),         intent(in)  :: rh2o
real(r8),         intent(in)  :: cpair
real(r8),         intent(in)  :: rhoh2o
real(r8),         intent(in)  :: tmelt_in        ! Freezing point of water (K)
real(r8),         intent(in)  :: latvap
real(r8),         intent(in)  :: latice
real(r8),         intent(in)  :: rhmini_in       ! Minimum rh for ice cloud fraction > 0.
real(r8),         intent(in)  :: micro_mg_dcs
logical,          intent(in)  :: use_hetfrz_classnuc_in
character(len=16),intent(in)  :: micro_mg_precip_frac_method_in  ! type of precipitation fraction method
real(r8),         intent(in)  :: micro_mg_berg_eff_factor_in     ! berg efficiency factor

character(128),   intent(out) :: errstring       ! Output status (non-blank for error return)

integer k

integer l,m, iaer
real(r8) surften       ! surface tension of water w/respect to air (N/m)
real(r8) arg
!-----------------------------------------------------------------------

!print*, "default_idx = ", default_idx
errstring = ' '

if( kind .ne. r8 ) then
   errstring = 'micro_mg_init: KIND of reals does not match'
   return
end if

!declarations for morrison codes (transforms variable names)

g= gravit                  !gravity
r= rair                    !Dry air Gas constant: note units(phys_constants are in J/K/kmol)
rv= rh2o                   !water vapor gas contstant
cpp = cpair                !specific heat of dry air
rhow = rhoh2o              !density of liquid water
tmelt = tmelt_in
rhmini = rhmini_in
micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
micro_mg_berg_eff_factor    = micro_mg_berg_eff_factor_in

! latent heats

xxlv = latvap         ! latent heat vaporization
xlf = latice          ! latent heat freezing
xxls = xxlv + xlf     ! latent heat of sublimation

! flags
use_hetfrz_classnuc = use_hetfrz_classnuc_in

! parameters for snow/rain fraction for convective clouds

tmax_fsnow = tmelt
tmin_fsnow = tmelt-5._r8

! parameters below from Reisner et al. (1998)
! density parameters (kg/m3)

rhosn = 250._r8    ! bulk density snow  (++ ceh)
rhoi = 500._r8     ! bulk density ice
rhow = 1000._r8    ! bulk density liquid


! fall speed parameters, V = aD^b
! V is in m/s

! droplets
ac = 3.e7_r8
bc = 2._r8

! snow
as = 11.72_r8
bs = 0.41_r8

! cloud ice
ai = 700._r8
bi = 1._r8

! rain
ar = 841.99667_r8
br = 0.8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d

pi= 3.1415927_r8

! cloud ice mass-diameter relationship

ci = rhoi*pi/6._r8
di = 3._r8

! snow mass-diameter relationship

cs = rhosn*pi/6._r8
ds = 3._r8

! drop mass-diameter relationship

cr = rhow*pi/6._r8
dr = 3._r8

! ventilation parameters for snow
! hall and prupacher

f1s = 0.86_r8
f2s = 0.28_r8

! collection efficiency, aggregation of cloud ice and snow

Eii = 0.1_r8

! collection efficiency, accretion of cloud water by rain

Ecr = 1.0_r8

! ventilation constants for rain

f1r = 0.78_r8
f2r = 0.32_r8

! autoconversion size threshold for cloud ice to snow (m)

Dcs = micro_mg_dcs

! smallest mixing ratio considered in microphysics

qsmall = 1.e-18_r8  

! immersion freezing parameters, bigg 1953

bimm = 100._r8
aimm = 0.66_r8

! typical air density at 850 mb

rhosu = 85000._r8/(rair * tmelt)

! mass of new crystal due to aerosol freezing and growth (kg)

mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)

! radius of contact nuclei aerosol (m)

rin = 0.1e-6_r8

! freezing temperature
tt0=273.15_r8

pi=4._r8*atan(1.0_r8)

!Range of cloudsat reflectivities (dBz) for analytic simulator
csmin= -30._r8
csmax= 26._r8
mindbz = -99._r8
!      minrefl = 10._r8**(mindbz/10._r8)
minrefl = 1.26e-10_r8

! Define constants to help speed up code (limit calls to gamma function)

cons1=gamma(1._r8+di)
cons4=gamma(1._r8+br)
cons5=gamma(4._r8+br)
cons6=gamma(1._r8+ds)
cons7=gamma(1._r8+bs)     
cons8=gamma(4._r8+bs)     
cons11=gamma(3._r8+bs)
cons13=gamma(5._r8/2._r8+br/2._r8)
cons14=gamma(5._r8/2._r8+bs/2._r8)
cons16=gamma(1._r8+bi)
cons17=gamma(4._r8+bi)
cons22=(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
cons23=dcs**3
cons24=dcs**2
cons25=dcs**bs
cons27=xxlv**2
cons28=xxls**2

lammaxi = 1._r8/10.e-6_r8
lammini = 1._r8/(2._r8*dcs)
lammaxr = 1._r8/20.e-6_r8
lamminr = 1._r8/500.e-6_r8
lammaxs = 1._r8/10.e-6_r8
lammins = 1._r8/2000.e-6_r8

end subroutine micro_mg_init

#define optimize

#ifdef optimize
!===============================================================================
!microphysics routine for each timestep goes here...

subroutine micro_mg_tend ( &
     microp_uniform, pcols, pver, ncol, top_lev, deltatin,&
     tn, qn, qc, qi, nc,                              &
     ni, p, pdel, cldn, liqcldf,                      &
     relvar, accre_enhan,                             &
     icecldf, rate1ord_cw2pr_st, naai, npccnin,       &
     rndst, nacon, tlat, qvlat, qctend,               &
     qitend, nctend, nitend, effc, effc_fn,           &
     effi, prect, preci, nevapr, evapsnow, am_evp_st, &
     prain, prodsnow, cmeout, deffi, pgamrad,         &
     lamcrad, qsout, dsout, rflx, sflx,               &
     qrout, reff_rain, reff_snow, qcsevap, qisevap,   &
     qvres, cmeiout, vtrmc, vtrmi, qcsedten,          &
     qisedten, prao, prco, mnuccco, mnuccto,          &
     msacwio, psacwso, bergso, bergo, melto,          &
     homoo, qcreso, prcio, praio, qireso,             &
     mnuccro, pracso, meltsdt, frzrdt, mnuccdo,       &
     nrout, nsout, refl, arefl, areflz,               &
     frefl, csrfl, acsrfl, fcsrfl, rercld,            &
     ncai, ncal, qrout2, qsout2, nrout2,              &
     nsout2, drout2, dsout2, freqs, freqr,            &
     nfice, prer_evap, do_cldice, errstring,          &
     tnd_qsnow, tnd_nsnow, re_ice,                    &
     frzimm, frzcnt, frzdep)

     include 'mpif.h'
! input arguments
logical,  intent(in) :: microp_uniform  ! True = configure uniform for sub-columns  False = use w/o sub-columns (standard)
integer,  intent(in) :: pcols                ! size of column (first) index
integer,  intent(in) :: pver                 ! number of layers in columns
integer,  intent(in) :: ncol                 ! number of columns
integer,  intent(in) :: top_lev              ! top level microphys is applied
real(r8), intent(in) :: deltatin             ! time step (s)
real(r8), intent(in) :: tn(pcols,pver)       ! input temperature (K)
real(r8), intent(in) :: qn(pcols,pver)       ! input h20 vapor mixing ratio (kg/kg)
real(r8), intent(in) :: relvar(pcols,pver)   ! relative variance of cloud water (-)
real(r8), intent(in) :: accre_enhan(pcols,pver) ! optional accretion enhancement factor (-)

! note: all input cloud variables are grid-averaged
real(r8), intent(inout) :: qc(pcols,pver)    ! cloud water mixing ratio (kg/kg)
real(r8), intent(inout) :: qi(pcols,pver)    ! cloud ice mixing ratio (kg/kg)
real(r8), intent(inout) :: nc(pcols,pver)    ! cloud water number conc (1/kg)
real(r8), intent(inout) :: ni(pcols,pver)    ! cloud ice number conc (1/kg)
real(r8), intent(in) :: p(pcols,pver)        ! air pressure (pa)
real(r8), intent(in) :: pdel(pcols,pver)     ! pressure difference across level (pa)
real(r8), intent(in) :: cldn(pcols,pver)     ! cloud fraction
real(r8), intent(in) :: icecldf(pcols,pver)  ! ice cloud fraction   
real(r8), intent(in) :: liqcldf(pcols,pver)  ! liquid cloud fraction
          
real(r8), intent(out) :: rate1ord_cw2pr_st(pcols,pver) ! 1st order rate for direct cw to precip conversion 
! used for scavenging
! Inputs for aerosol activation
real(r8), intent(in) :: naai(pcols,pver)      ! ice nulceation number (from microp_aero_ts) 
real(r8), intent(in) :: npccnin(pcols,pver)   ! ccn activated number tendency (from microp_aero_ts)
real(r8), intent(in) :: rndst(pcols,pver,4)   ! radius of 4 dust bins for contact freezing (from microp_aero_ts)
real(r8), intent(in) :: nacon(pcols,pver,4)   ! number in 4 dust bins for contact freezing  (from microp_aero_ts)

! Used with CARMA cirrus microphysics
! (or similar external microphysics model)
logical,  intent(in) :: do_cldice             ! Prognosing cldice

! output arguments

real(r8), intent(out) :: tlat(pcols,pver)    ! latent heating rate       (W/kg)
real(r8), intent(out) :: qvlat(pcols,pver)   ! microphysical tendency qv (1/s)
real(r8), intent(out) :: qctend(pcols,pver)  ! microphysical tendency qc (1/s) 
real(r8), intent(out) :: qitend(pcols,pver)  ! microphysical tendency qi (1/s)
real(r8), intent(out) :: nctend(pcols,pver)  ! microphysical tendency nc (1/(kg*s))
real(r8), intent(out) :: nitend(pcols,pver)  ! microphysical tendency ni (1/(kg*s))
real(r8), intent(out) :: effc(pcols,pver)    ! droplet effective radius (micron)
real(r8), intent(out) :: effc_fn(pcols,pver) ! droplet effective radius, assuming nc = 1.e8 kg-1
real(r8), intent(out) :: effi(pcols,pver)    ! cloud ice effective radius (micron)
real(r8), intent(out) :: prect(pcols)        ! surface precip rate (m/s)
real(r8), intent(out) :: preci(pcols)        ! cloud ice/snow precip rate (m/s)
real(r8), intent(out) :: nevapr(pcols,pver)  ! evaporation rate of rain + snow
real(r8), intent(out) :: evapsnow(pcols,pver)! sublimation rate of snow
real(r8), intent(out) :: am_evp_st(pcols,pver)! stratiform evaporation area
real(r8), intent(out) :: prain(pcols,pver)   ! production of rain + snow
real(r8), intent(out) :: prodsnow(pcols,pver)! production of snow
real(r8), intent(out) :: cmeout(pcols,pver)  ! evap/sub of cloud
real(r8), intent(out) :: deffi(pcols,pver)   ! ice effective diameter for optics (radiation)
real(r8), intent(out) :: pgamrad(pcols,pver) ! ice gamma parameter for optics (radiation)
real(r8), intent(out) :: lamcrad(pcols,pver) ! slope of droplet distribution for optics (radiation)
real(r8), intent(out) :: qsout(pcols,pver)   ! snow mixing ratio (kg/kg)
real(r8), intent(out) :: dsout(pcols,pver)   ! snow diameter (m)
real(r8), intent(out) :: rflx(pcols,pver+1)  ! grid-box average rain flux (kg m^-2 s^-1)
real(r8), intent(out) :: sflx(pcols,pver+1)  ! grid-box average snow flux (kg m^-2 s^-1)
real(r8), intent(out) :: qrout(pcols,pver)     ! grid-box average rain mixing ratio (kg/kg)
real(r8), intent(inout) :: reff_rain(pcols,pver) ! rain effective radius (micron)
real(r8), intent(inout) :: reff_snow(pcols,pver) ! snow effective radius (micron)
real(r8), intent(out) :: qcsevap(pcols,pver) ! cloud water evaporation due to sedimentation
real(r8), intent(out) :: qisevap(pcols,pver) ! cloud ice sublimation due to sublimation
real(r8), intent(out) :: qvres(pcols,pver) ! residual condensation term to ensure RH < 100%
real(r8), intent(out) :: cmeiout(pcols,pver) ! grid-mean cloud ice sub/dep
real(r8), intent(out) :: vtrmc(pcols,pver) ! mass-weighted cloud water fallspeed
real(r8), intent(out) :: vtrmi(pcols,pver) ! mass-weighted cloud ice fallspeed
real(r8), intent(out) :: qcsedten(pcols,pver) ! qc sedimentation tendency
real(r8), intent(out) :: qisedten(pcols,pver) ! qi sedimentation tendency
! microphysical process rates for output (mixing ratio tendencies)
real(r8), intent(out) :: prao(pcols,pver) ! accretion of cloud by rain 
real(r8), intent(out) :: prco(pcols,pver) ! autoconversion of cloud to rain
real(r8), intent(out) :: mnuccco(pcols,pver) ! mixing rat tend due to immersion freezing
real(r8), intent(out) :: mnuccto(pcols,pver) ! mixing ratio tend due to contact freezing
real(r8), intent(out) :: msacwio(pcols,pver) ! mixing ratio tend due to H-M splintering
real(r8), intent(out) :: psacwso(pcols,pver) ! collection of cloud water by snow
real(r8), intent(out) :: bergso(pcols,pver) ! bergeron process on snow
real(r8), intent(out) :: bergo(pcols,pver) ! bergeron process on cloud ice
real(r8), intent(out) :: melto(pcols,pver) ! melting of cloud ice
real(r8), intent(out) :: homoo(pcols,pver) ! homogeneos freezign cloud water
real(r8), intent(out) :: qcreso(pcols,pver) ! residual cloud condensation due to removal of excess supersat
real(r8), intent(out) :: prcio(pcols,pver) ! autoconversion of cloud ice to snow
real(r8), intent(out) :: praio(pcols,pver) ! accretion of cloud ice by snow
real(r8), intent(out) :: qireso(pcols,pver) ! residual ice deposition due to removal of excess supersat
real(r8), intent(out) :: mnuccro(pcols,pver) ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
real(r8), intent(out) :: pracso (pcols,pver) ! mixing ratio tendency due to accretion of rain by snow (1/s)
real(r8), intent(out) :: meltsdt(pcols,pver) ! latent heating rate due to melting of snow  (W/kg)
real(r8), intent(out) :: frzrdt (pcols,pver) ! latent heating rate due to homogeneous freezing of rain (W/kg)
real(r8), intent(out) :: mnuccdo(pcols,pver) ! mass tendency from ice nucleation
real(r8), intent(out) :: nrout(pcols,pver) ! rain number concentration (1/m3)
real(r8), intent(out) :: nsout(pcols,pver) ! snow number concentration (1/m3)
real(r8), intent(out) :: refl(pcols,pver)    ! analytic radar reflectivity        
real(r8), intent(out) :: arefl(pcols,pver)  !average reflectivity will zero points outside valid range
real(r8), intent(out) :: areflz(pcols,pver)  !average reflectivity in z.
real(r8), intent(out) :: frefl(pcols,pver)
real(r8), intent(out) :: csrfl(pcols,pver)   !cloudsat reflectivity 
real(r8), intent(out) :: acsrfl(pcols,pver)  !cloudsat average
real(r8), intent(out) :: fcsrfl(pcols,pver)
real(r8), intent(out) :: rercld(pcols,pver) ! effective radius calculation for rain + cloud
real(r8), intent(out) :: ncai(pcols,pver) ! output number conc of ice nuclei available (1/m3)
real(r8), intent(out) :: ncal(pcols,pver) ! output number conc of CCN (1/m3)
real(r8), intent(out) :: qrout2(pcols,pver)
real(r8), intent(out) :: qsout2(pcols,pver)
real(r8), intent(out) :: nrout2(pcols,pver)
real(r8), intent(out) :: nsout2(pcols,pver)
real(r8), intent(out) :: drout2(pcols,pver) ! mean rain particle diameter (m)
real(r8), intent(out) :: dsout2(pcols,pver) ! mean snow particle diameter (m)
real(r8), intent(out) :: freqs(pcols,pver)
real(r8), intent(out) :: freqr(pcols,pver)
real(r8), intent(out) :: nfice(pcols,pver)
real(r8), intent(out) :: prer_evap(pcols,pver)

real(r8) :: nevapr2(pcols,pver)

character(128),   intent(out) :: errstring       ! Output status (non-blank for error return)

! Tendencies calculated by external schemes that can replace MG's native
! process tendencies.

! Used with CARMA cirrus microphysics
! (or similar external microphysics model)
real(r8), intent(in), pointer :: tnd_qsnow(:,:) ! snow mass tendency (kg/kg/s)
real(r8), intent(in), pointer :: tnd_nsnow(:,:) ! snow number tendency (#/kg/s)
real(r8), intent(in), pointer :: re_ice(:,:)    ! ice effective radius (m)

! From external ice nucleation.
real(r8), intent(in), pointer :: frzimm(:,:) ! Number tendency due to immersion freezing (1/cm3)
real(r8), intent(in), pointer :: frzcnt(:,:) ! Number tendency due to contact freezing (1/cm3)
real(r8), intent(in), pointer :: frzdep(:,:) ! Number tendency due to deposition nucleation (1/cm3)

! local workspace
! all units mks unless otherwise stated

! Additional constants to help speed up code
real(r8) :: cons2
real(r8) :: cons3
real(r8) :: cons9
real(r8) :: cons10
real(r8) :: cons12
real(r8) :: cons15
real(r8) :: cons18
real(r8) :: cons19
real(r8) :: cons20

! temporary variables for sub-stepping 
real(r8) :: t1(pcols,pver)
real(r8) :: q1(pcols,pver)
real(r8) :: qc1(pcols,pver)
real(r8) :: qi1(pcols,pver)
real(r8) :: nc1(pcols,pver)
real(r8) :: ni1(pcols,pver)
real(r8) :: tlat1(pcols,pver)
real(r8) :: qvlat1(pcols,pver)
real(r8) :: qctend1(pcols,pver)
real(r8) :: qitend1(pcols,pver)
real(r8) :: nctend1(pcols,pver)
real(r8) :: nitend1(pcols,pver)
real(r8) :: prect1(pcols)
real(r8) :: preci1(pcols)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(r8) :: deltat        ! sub-time step (s)
real(r8) :: omsm    ! number near unity for round-off issues
real(r8) :: dto2    ! dt/2 (s)
real(r8) :: mincld  ! minimum allowed cloud fraction
real(r8) :: q(pcols,pver) ! water vapor mixing ratio (kg/kg)
real(r8) :: t(pcols,pver) ! temperature (K)
real(r8) :: rho(pcols,pver) ! air density (kg m-3)
real(r8) :: dv(pcols,pver)  ! diffusivity of water vapor in air
real(r8) :: mu(pcols,pver)  ! viscocity of air
real(r8) :: sc(pcols,pver)  ! schmidt number
real(r8) :: kap(pcols,pver) ! thermal conductivity of air
real(r8) :: rhof(pcols,pver) ! air density correction factor for fallspeed
real(r8) :: cldmax(pcols,pver) ! precip fraction assuming maximum overlap
real(r8) :: cldm(pcols,pver)   ! cloud fraction
real(r8) :: icldm(pcols,pver)   ! ice cloud fraction
real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction
real(r8) :: icwc(pcols)    ! in cloud water content (liquid+ice)
real(r8) :: calpha(pcols)  ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cbeta(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cbetah(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cgamma(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cgamah(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: rcgama(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec1(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec2(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec3(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec4(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: qtmp ! dummy qv 
real(r8) :: dum  ! temporary dummy variable

real(r8) :: cme(pcols,pver)  ! total (liquid+ice) cond/evap rate of cloud

real(r8) :: cmei(pcols,pver) ! dep/sublimation rate of cloud ice
real(r8) :: cwml(pcols,pver) ! cloud water mixing ratio
real(r8) :: cwmi(pcols,pver) ! cloud ice mixing ratio
real(r8) :: nnuccd(pver)   ! ice nucleation rate from deposition/cond.-freezing
real(r8) :: mnuccd(pver)   ! mass tendency from ice nucleation
real(r8) :: qcld              ! total cloud water
real(r8) :: lcldn(pcols,pver) ! fractional coverage of new liquid cloud
real(r8) :: lcldo(pcols,pver) ! fractional coverage of old liquid cloud
real(r8) :: nctend_mixnuc(pcols,pver)
real(r8) :: arg ! argument of erfc

! for calculation of rate1ord_cw2pr_st
real(r8) :: qcsinksum_rate1ord(pver)   ! sum over iterations of cw to precip sink
real(r8) :: qcsum_rate1ord(pver)    ! sum over iterations of cloud water       

real(r8) :: alpha

real(r8) :: dum1,dum2   !general dummy variables

real(r8) :: npccn(pver)     ! droplet activation rate
real(r8) :: qcic(pcols,pver) ! in-cloud cloud liquid mixing ratio
real(r8) :: qiic(pcols,pver) ! in-cloud cloud ice mixing ratio
real(r8) :: qniic(pcols,pver) ! in-precip snow mixing ratio
real(r8) :: qric(pcols,pver) ! in-precip rain mixing ratio
real(r8) :: ncic(pcols,pver) ! in-cloud droplet number conc
real(r8) :: niic(pcols,pver) ! in-cloud cloud ice number conc
real(r8) :: nsic(pcols,pver) ! in-precip snow number conc
real(r8) :: nric(pcols,pver) ! in-precip rain number conc
real(r8) :: lami(pver) ! slope of cloud ice size distr
real(r8) :: n0i(pver) ! intercept of cloud ice size distr
real(r8) :: lamc(pver) ! slope of cloud liquid size distr
real(r8) :: n0c(pver) ! intercept of cloud liquid size distr
real(r8) :: lams(pver) ! slope of snow size distr
real(r8) :: n0s(pver) ! intercept of snow size distr
real(r8) :: lamr(pver) ! slope of rain size distr
real(r8) :: n0r(pver) ! intercept of rain size distr
real(r8) :: cdist1(pver) ! size distr parameter to calculate droplet freezing
! combined size of precip & cloud drops
real(r8) :: arcld(pcols,pver) ! averaging control flag
real(r8) :: Actmp  !area cross section of drops
real(r8) :: Artmp  !area cross section of rain

real(r8) :: pgam(pver) ! spectral width parameter of droplet size distr
real(r8) :: lammax  ! maximum allowed slope of size distr
real(r8) :: lammin  ! minimum allowed slope of size distr
real(r8) :: nacnt   ! number conc of contact ice nuclei
real(r8) :: mnuccc(pver) ! mixing ratio tendency due to freezing of cloud water
real(r8) :: nnuccc(pver) ! number conc tendency due to freezing of cloud water

real(r8) :: mnucct(pver) ! mixing ratio tendency due to contact freezing of cloud water
real(r8) :: nnucct(pver) ! number conc tendency due to contact freezing of cloud water
real(r8) :: msacwi(pver) ! mixing ratio tendency due to HM ice multiplication
real(r8) :: nsacwi(pver) ! number conc tendency due to HM ice multiplication

real(r8) :: prc(pver) ! qc tendency due to autoconversion of cloud droplets
real(r8) :: nprc(pver) ! number conc tendency due to autoconversion of cloud droplets
real(r8) :: nprc1(pver) ! qr tendency due to autoconversion of cloud droplets
real(r8) :: nsagg(pver) ! ns tendency due to self-aggregation of snow
real(r8) :: dc0  ! mean size droplet size distr
real(r8) :: ds0  ! mean size snow size distr (area weighted)
real(r8) :: eci  ! collection efficiency for riming of snow by droplets
real(r8) :: psacws(pver) ! mixing rat tendency due to collection of droplets by snow
real(r8) :: npsacws(pver) ! number conc tendency due to collection of droplets by snow
real(r8) :: uni ! number-weighted cloud ice fallspeed
real(r8) :: umi ! mass-weighted cloud ice fallspeed
real(r8) :: uns(pver) ! number-weighted snow fallspeed
real(r8) :: ums(pver) ! mass-weighted snow fallspeed
real(r8) :: unr(pver) ! number-weighted rain fallspeed
real(r8) :: umr(pver) ! mass-weighted rain fallspeed
real(r8) :: unc ! number-weighted cloud droplet fallspeed
real(r8) :: umc ! mass-weighted cloud droplet fallspeed
real(r8) :: pracs(pver) ! mixing rat tendency due to collection of rain by snow
real(r8) :: npracs(pver) ! number conc tendency due to collection of rain by snow
real(r8) :: mnuccr(pver) ! mixing rat tendency due to freezing of rain
real(r8) :: nnuccr(pver) ! number conc tendency due to freezing of rain
real(r8) :: pra(pver) ! mixing rat tendnency due to accretion of droplets by rain
real(r8) :: npra(pver) ! nc tendnency due to accretion of droplets by rain
real(r8) :: nragg(pver) ! nr tendency due to self-collection of rain
real(r8) :: prci(pver) ! mixing rat tendency due to autoconversion of cloud ice to snow
real(r8) :: nprci(pver) ! number conc tendency due to autoconversion of cloud ice to snow
real(r8) :: prai(pver) ! mixing rat tendency due to accretion of cloud ice by snow
real(r8) :: nprai(pver) ! number conc tendency due to accretion of cloud ice by snow
real(r8) :: qvs ! liquid saturation vapor mixing ratio
real(r8) :: qvi ! ice saturation vapor mixing ratio
real(r8) :: dqsdt ! change of sat vapor mixing ratio with temperature
real(r8) :: dqsidt ! change of ice sat vapor mixing ratio with temperature
real(r8) :: ab ! correction factor for rain evap to account for latent heat
real(r8) :: qclr ! water vapor mixing ratio in clear air
real(r8) :: abi ! correction factor for snow sublimation to account for latent heat
real(r8) :: epss ! 1/ sat relaxation timescale for snow
real(r8) :: epsr ! 1/ sat relaxation timescale for rain
real(r8) :: pre(pver) ! rain mixing rat tendency due to evaporation
real(r8) :: prds(pver) ! snow mixing rat tendency due to sublimation
real(r8) :: qce ! dummy qc for conservation check
real(r8) :: qie ! dummy qi for conservation check
real(r8) :: nce ! dummy nc for conservation check
real(r8) :: nie ! dummy ni for conservation check
real(r8) :: ratio ! parameter for conservation check
real(r8) :: dumc(pcols,pver) ! dummy in-cloud qc
real(r8) :: dumnc(pcols,pver) ! dummy in-cloud nc
real(r8) :: dumi(pcols,pver) ! dummy in-cloud qi
real(r8) :: dumni(pcols,pver) ! dummy in-cloud ni
real(r8) :: dums(pcols,pver) ! dummy in-cloud snow mixing rat
real(r8) :: dumns(pcols,pver) ! dummy in-cloud snow number conc
real(r8) :: dumr(pcols,pver) ! dummy in-cloud rain mixing rat
real(r8) :: dumnr(pcols,pver) ! dummy in-cloud rain number conc
! below are parameters for cloud water and cloud ice sedimentation calculations
real(r8) :: fr(pver)
real(r8) :: fnr(pver)
real(r8) :: fc(pver)
real(r8) :: fnc(pver)
real(r8) :: fi(pver)
real(r8) :: fni(pver)
real(r8) :: fs(pver)
real(r8) :: fns(pver)
real(r8) :: faloutr(pver)
real(r8) :: faloutnr(pver)
real(r8) :: faloutc(pver)
real(r8) :: faloutnc(pver)
real(r8) :: falouti(pver)
real(r8) :: faloutni(pver)
real(r8) :: falouts(pver)
real(r8) :: faloutns(pver)
real(r8) :: faltndr
real(r8) :: faltndnr
real(r8) :: faltndc
real(r8) :: faltndnc
real(r8) :: faltndi
real(r8) :: faltndni
real(r8) :: faltnds
real(r8) :: faltndns
real(r8) :: faltndqie
real(r8) :: faltndqce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(r8) :: relhum(pcols,pver) ! relative humidity
real(r8) :: csigma(pcols) ! parameter for cond/evap of cloud water/ice
real(r8) :: rgvm ! max fallspeed for all species
real(r8) :: arn(pcols,pver) ! air density corrected rain fallspeed parameter
real(r8) :: asn(pcols,pver) ! air density corrected snow fallspeed parameter
real(r8) :: acn(pcols,pver) ! air density corrected cloud droplet fallspeed parameter
real(r8) :: ain(pcols,pver) ! air density corrected cloud ice fallspeed parameter
real(r8) :: nsubi(pver) ! evaporation of cloud ice number
real(r8) :: nsubc(pver) ! evaporation of droplet number
real(r8) :: nsubs(pver) ! evaporation of snow number
real(r8) :: nsubr(pver) ! evaporation of rain number
real(r8) :: mtime ! factor to account for droplet activation timescale
real(r8) :: dz(pcols,pver) ! height difference across model vertical level


!! add precip flux variables for sub-stepping
real(r8) :: rflx1(pcols,pver+1)
real(r8) :: sflx1(pcols,pver+1)

! returns from function/subroutine calls
real(r8) :: tsp(pcols,pver)      ! saturation temp (K)
real(r8) :: qsp(pcols,pver)      ! saturation mixing ratio (kg/kg)
real(r8) :: qsphy(pcols,pver)      ! saturation mixing ratio (kg/kg): hybrid rh
real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
real(r8) :: esl(pcols,pver)      ! liquid sat vapor pressure (pa)
real(r8) :: esi(pcols,pver)      ! ice sat vapor pressure (pa)

! sum of source/sink terms for diagnostic precip

real(r8) :: qnitend(pcols,pver) ! snow mixing ratio source/sink term
real(r8) :: nstend(pcols,pver)  ! snow number concentration source/sink term
real(r8) :: qrtend(pcols,pver) ! rain mixing ratio source/sink term
real(r8) :: nrtend(pcols,pver)  ! rain number concentration source/sink term
real(r8) :: qrtot ! vertically-integrated rain mixing rat source/sink term
real(r8) :: nrtot ! vertically-integrated rain number conc source/sink term
real(r8) :: qstot ! vertically-integrated snow mixing rat source/sink term
real(r8) :: nstot ! vertically-integrated snow number conc source/sink term

! new terms for Bergeron process

real(r8) :: dumnnuc ! provisional ice nucleation rate (for calculating bergeron)
real(r8) :: ninew  ! provisional cloud ice number conc (for calculating bergeron)
real(r8) :: qinew ! provisional cloud ice mixing ratio (for calculating bergeron)
real(r8) :: qvl  ! liquid sat mixing ratio   
real(r8) :: epsi ! 1/ sat relaxation timecale for cloud ice
real(r8) :: prd ! provisional deposition rate of cloud ice at water sat 
real(r8) :: berg(pcols,pver) ! mixing rat tendency due to bergeron process for cloud ice
real(r8) :: bergs(pver) ! mixing rat tendency due to bergeron process for snow

!bergeron terms
real(r8) :: bergtsf   !bergeron timescale to remove all liquid
real(r8) :: rhin      !modified RH for vapor deposition

! diagnostic rain/snow for output to history
! values are in-precip (local) !!!!

real(r8) :: drout(pcols,pver)     ! rain diameter (m)

!averageed rain/snow for history
real(r8) :: dumfice

!ice nucleation, droplet activation
real(r8) :: dum2i(pcols,pver) ! number conc of ice nuclei available (1/kg)
real(r8) :: dum2l(pcols,pver) ! number conc of CCN (1/kg)
real(r8) :: ncmax
real(r8) :: nimax

real(r8) :: qcvar     ! 1/relative variance of sub-grid qc

! loop array variables
integer i,k,nstep,n, l
integer ii,kk, m

! loop variables for sub-step solution
integer iter,it,ltrue(pcols)

! used in contact freezing via dust particles
real(r8)  tcnt, viscosity, mfp
real(r8)  slip1, slip2, slip3, slip4
!        real(r8)  dfaer1, dfaer2, dfaer3, dfaer4
!        real(r8)  nacon1,nacon2,nacon3,nacon4
real(r8)  ndfaer1, ndfaer2, ndfaer3, ndfaer4
real(r8)  nslip1, nslip2, nslip3, nslip4

! used in ice effective radius
real(r8)  bbi, cci, ak, iciwc, rvi

! used in Bergeron processe and water vapor deposition
real(r8)  Tk, deles, Aprpr, Bprpr, Cice, qi0, Crate, qidep

! mean cloud fraction over the time step
real(r8)  cldmw(pcols,pver)

! used in secondary ice production
real(r8) ni_secp

! variabels to check for RH after rain evap

real(r8) :: esn
real(r8) :: qsn
real(r8) :: ttmp



real(r8) :: rainrt(pcols,pver)  ! rain rate for reflectivity calculation
real(r8) :: rainrt1(pcols,pver)
real(r8) :: tmp

real(r8) dmc,ssmc,dstrn  ! variables for modal scheme.

real(r8), parameter :: cdnl    = 0.e6_r8    ! cloud droplet number limiter

! heterogeneous freezing
real(r8) :: mnudep(pver) ! mixing ratio tendency due to deposition of water vapor
real(r8) :: nnudep(pver) ! number conc tendency due to deposition of water vapor
real(r8) :: con1 ! work cnstant
real(r8) :: r3lx ! Mean volume radius (m)
real(r8) :: mi0l
real(r8) :: frztmp

logical  :: do_clubb_sgs

! Temporal variables for transposition -- Ping Xu
real(r8) :: tn_permute(pver, pcols)
real(r8) :: qn_permute(pver, pcols)
real(r8) :: relvar_permute(pver, pcols)
real(r8) :: accre_enhan_permute(pver, pcols)

real(r8) :: qc_permute(pver, pcols)
real(r8) :: qi_permute(pver, pcols)
real(r8) :: nc_permute(pver, pcols)
real(r8) :: ni_permute(pver, pcols)
real(r8) :: p_permute(pver, pcols)
real(r8) :: pdel_permute(pver, pcols)
real(r8) :: cldn_permute(pver, pcols)
real(r8) :: icecldf_permute(pver, pcols)
real(r8) :: liqcldf_permute(pver, pcols)
real(r8) :: rate1ord_cw2pr_st_permute(pver, pcols)

real(r8) :: naai_permute(pver, pcols)
real(r8) :: npccnin_permute(pver, pcols)
real(r8) :: rndst_permute(pver, pcols, 4)
real(r8) :: nacon_permute(pver, pcols, 4)

real(r8) :: tlat_permute(pver, pcols)
real(r8) :: qvlat_permute(pver, pcols)
real(r8) :: qctend_permute(pver, pcols)
real(r8) :: qitend_permute(pver, pcols)
real(r8) :: nctend_permute(pver, pcols)
real(r8) :: nitend_permute(pver, pcols)

real(r8) :: effc_permute(pver, pcols)
real(r8) :: effc_fn_permute(pver, pcols)
real(r8) :: effi_permute(pver, pcols)
!real(r8) :: prect_permute(pver, pcols)
!real(r8) :: preci_permute(pver, pcols)
real(r8) :: nevapr_permute(pver, pcols)
real(r8) :: evapsnow_permute(pver, pcols)
real(r8) :: am_evp_st_permute(pver, pcols)
real(r8) :: prain_permute(pver, pcols)
real(r8) :: prodsnow_permute(pver, pcols)
real(r8) :: cmeout_permute(pver, pcols)
real(r8) :: deffi_permute(pver, pcols)
real(r8) :: pgamrad_permute(pver, pcols)
real(r8) :: lamcrad_permute(pver, pcols)
real(r8) :: qsout_permute(pver, pcols)
real(r8) :: dsout_permute(pver, pcols)

real(r8) :: rflx_permute(pver+1, pcols)
real(r8) :: sflx_permute(pver+1, pcols)
real(r8) :: qrout_permute(pver, pcols)
real(r8) :: reff_rain_permute(pver, pcols)
real(r8) :: reff_snow_permute(pver, pcols)
real(r8) :: qcsevap_permute(pver, pcols)
real(r8) :: qisevap_permute(pver, pcols)
real(r8) :: qvres_permute(pver, pcols)
real(r8) :: cmeiout_permute(pver, pcols)
real(r8) :: vtrmc_permute(pver, pcols)
real(r8) :: vtrmi_permute(pver, pcols)
real(r8) :: qcsedten_permute(pver, pcols)
real(r8) :: qisedten_permute(pver, pcols)
real(r8) :: prao_permute(pver, pcols)
real(r8) :: prco_permute(pver, pcols)
real(r8) :: mnuccco_permute(pver, pcols)
real(r8) :: mnuccto_permute(pver, pcols)
real(r8) :: msacwio_permute(pver, pcols)
real(r8) :: psacwso_permute(pver, pcols)
real(r8) :: bergso_permute(pver, pcols)
real(r8) :: bergo_permute(pver, pcols)
real(r8) :: melto_permute(pver, pcols)
real(r8) :: homoo_permute(pver, pcols)
real(r8) :: qcreso_permute(pver, pcols)
real(r8) :: prcio_permute(pver, pcols)
real(r8) :: praio_permute(pver, pcols)
real(r8) :: qireso_permute(pver, pcols)
real(r8) :: mnuccro_permute(pver, pcols)
real(r8) :: pracso_permute(pver, pcols)
real(r8) :: meltsdt_permute(pver, pcols)
real(r8) :: frzrdt_permute(pver, pcols)
real(r8) :: mnuccdo_permute(pver, pcols)
real(r8) :: nrout_permute(pver, pcols)
real(r8) :: nsout_permute(pver, pcols)
real(r8) :: refl_permute(pver, pcols)
real(r8) :: arefl_permute(pver, pcols)
real(r8) :: areflz_permute(pver, pcols)
real(r8) :: frefl_permute(pver, pcols)
real(r8) :: csrfl_permute(pver, pcols)
real(r8) :: acsrfl_permute(pver, pcols)
real(r8) :: fcsrfl_permute(pver, pcols)
real(r8) :: rercld_permute(pver, pcols)
real(r8) :: ncai_permute(pver, pcols)
real(r8) :: ncal_permute(pver, pcols)
real(r8) :: qrout2_permute(pver, pcols)
real(r8) :: qsout2_permute(pver, pcols)
real(r8) :: nrout2_permute(pver, pcols)
real(r8) :: nsout2_permute(pver, pcols)
real(r8) :: drout2_permute(pver, pcols)
real(r8) :: dsout2_permute(pver, pcols)
real(r8) :: freqs_permute(pver, pcols)
real(r8) :: freqr_permute(pver, pcols)
real(r8) :: nfice_permute(pver, pcols)
real(r8) :: prer_evap_permute(pver, pcols)
real(r8) :: nevapr2_permute(pver, pcols)

real(r8) :: tnd_qsnow_permute(pver, pcols)
real(r8) :: tnd_nsnow_permute(pver, pcols)
real(r8) :: re_ice_permute(pver, pcols)
real(r8) :: frzimm_permute(pver, pcols)
real(r8) :: frzcnt_permute(pver, pcols)
real(r8) :: frzdep_permute(pver, pcols)

real(r8) :: t1_permute(pver, pcols)
real(r8) :: q1_permute(pver, pcols)
real(r8) :: qc1_permute(pver, pcols)
real(r8) :: qi1_permute(pver, pcols)
real(r8) :: nc1_permute(pver, pcols)
real(r8) :: ni1_permute(pver, pcols)
real(r8) :: tlat1_permute(pver, pcols)
real(r8) :: qvlat1_permute(pver, pcols)
real(r8) :: qctend1_permute(pver, pcols)
real(r8) :: qitend1_permute(pver, pcols)
real(r8) :: nctend1_permute(pver, pcols)
real(r8) :: nitend1_permute(pver, pcols)
real(r8) :: q_permute(pver, pcols)
real(r8) :: t_permute(pver, pcols)
real(r8) :: rho_permute(pver, pcols)
real(r8) :: dv_permute(pver, pcols)
real(r8) :: mu_permute(pver, pcols)
real(r8) :: sc_permute(pver, pcols)
real(r8) :: kap_permute(pver, pcols)
real(r8) :: rhof_permute(pver, pcols)
real(r8) :: cldmax_permute(pver, pcols)
real(r8) :: cldm_permute(pver, pcols)
real(r8) :: icldm_permute(pver, pcols)
real(r8) :: lcldm_permute(pver, pcols)
real(r8) :: cme_permute(pver, pcols)
real(r8) :: cmei_permute(pver, pcols)
real(r8) :: cwml_permute(pver, pcols)
real(r8) :: cwmi_permute(pver, pcols)
real(r8) :: lcldn_permute(pver, pcols)
real(r8) :: lcldo_permute(pver, pcols)
real(r8) :: nctend_mixnuc_permute(pver, pcols)

real(r8) :: qcic_permute(pver, pcols)
real(r8) :: qiic_permute(pver, pcols)
real(r8) :: qniic_permute(pver, pcols)
real(r8) :: qric_permute(pver, pcols)
real(r8) :: ncic_permute(pver, pcols)
real(r8) :: niic_permute(pver, pcols)
real(r8) :: nsic_permute(pver, pcols)
real(r8) :: nric_permute(pver, pcols)
real(r8) :: arcld_permute(pver, pcols)
real(r8) :: dumc_permute(pver, pcols)
real(r8) :: dumnc_permute(pver, pcols)
real(r8) :: dumi_permute(pver, pcols)
real(r8) :: dumni_permute(pver, pcols)
real(r8) :: dums_permute(pver, pcols)
real(r8) :: dumns_permute(pver, pcols)
real(r8) :: dumr_permute(pver, pcols)
real(r8) :: dumnr_permute(pver, pcols)
real(r8) :: relhum_permute(pver, pcols)
real(r8) :: arn_permute(pver, pcols)
real(r8) :: asn_permute(pver, pcols)
real(r8) :: acn_permute(pver, pcols)
real(r8) :: ain_permute(pver, pcols)
real(r8) :: dz_permute(pver, pcols)

real(r8) :: rflx1_permute(pver+1, pcols)
real(r8) :: sflx1_permute(pver+1, pcols)
real(r8) :: tsp_permute(pver, pcols)
real(r8) :: qsp_permute(pver, pcols)
real(r8) :: qsphy_permute(pver, pcols)
real(r8) :: esl_permute(pver, pcols)
real(r8) :: esi_permute(pver, pcols)
real(r8) :: qnitend_permute(pver, pcols)
real(r8) :: nstend_permute(pver, pcols)
real(r8) :: qrtend_permute(pver, pcols)
real(r8) :: nrtend_permute(pver, pcols)
real(r8) :: berg_permute(pver, pcols)
real(r8) :: drout_permute(pver, pcols)
real(r8) :: dum2i_permute(pver, pcols)
real(r8) :: dum2l_permute(pver, pcols)
real(r8) :: cldmw_permute(pver, pcols)
real(r8) :: rainrt_permute(pver, pcols)
real(r8) :: rainrt1_permute(pver, pcols)
real(r8) :: prect_permute(pcols)        ! surface precip rate (m/s)
real(r8) :: preci_permute(pcols)        ! cloud ice/snow precip rate (m/s)

real(r8), parameter :: svp_tboil = 373.16_r8
real(r8), parameter :: svp_ttrice = 20.00_r8  ! transition range from es over H2O to es over ice

type micro_mg1_0_args
    real(r8) :: svp_tmelt
    real(r8) :: svp_h2otrip
    real(r8) :: svp_tboil
    real(r8) :: svp_ttrice
    real(r8) :: svp_epsilo
    ! init
    real(r8) :: g              !gravity
    real(r8) :: r              !Dry air Gas constant
    real(r8) :: rv             !water vapor gas contstant
    real(r8) :: cpp            !specific heat of dry air
    real(r8) :: rhow           !density of liquid water
    real(r8) :: tmelt          ! Freezing point of water (K)
    real(r8) :: xxlv           ! latent heat of vaporization
    real(r8) :: xlf            !latent heat of freezing
    real(r8) :: xxls           !latent heat of sublimation
    real(r8) :: rhosn  ! bulk density snow
    real(r8) :: rhoi   ! bulk density ice
    real(r8) :: ac
    real(r8) :: bc
    real(r8) :: as
    real(r8) :: bs
    real(r8) :: ai
    real(r8) :: bi
    real(r8) :: ar
    real(r8) :: br  !fall speed parameters 
    real(r8) :: ci
    real(r8) :: di    !ice mass-diameter relation parameters
    real(r8) :: cs
    real(r8) :: ds    !snow mass-diameter relation parameters
    real(r8) :: cr
    real(r8) :: dr    !drop mass-diameter relation parameters
    real(r8) :: f1s
    real(r8) :: f2s  !ventilation param for snow
    real(r8) :: Eii      !collection efficiency aggregation of ice
    real(r8) :: Ecr      !collection efficiency cloud droplets/rain
    real(r8) :: f1r
    real(r8) :: f2r  !ventilation param for rain
    real(r8) :: DCS      !autoconversion size threshold
    real(r8) :: qsmall   !min mixing ratio 
    real(r8) :: bimm
    real(r8) :: aimm !immersion freezing
    real(r8) :: rhosu     !typical 850mn air density
    real(r8) :: mi0       ! new crystal mass
    real(r8) :: rin       ! radius of contact nuclei
    real(r8) :: pi       ! pi
    real(r8) :: cons1
    real(r8) :: cons4
    real(r8) :: cons5
    real(r8) :: cons6
    real(r8) :: cons7
    real(r8) :: cons8
    real(r8) :: cons11
    real(r8) :: cons13
    real(r8) :: cons14
    real(r8) :: cons16
    real(r8) :: cons17
    real(r8) :: cons22
    real(r8) :: cons23
    real(r8) :: cons24
    real(r8) :: cons25
    real(r8) :: cons27
    real(r8) :: cons28
    real(r8) :: lammini
    real(r8) :: lammaxi
    real(r8) :: lamminr
    real(r8) :: lammaxr
    real(r8) :: lammins
    real(r8) :: lammaxs
    real(r8) :: tmax_fsnow ! max temperature for transition to convective snow
    real(r8) :: tmin_fsnow ! min temperature for transition to convective snow
    real(r8) :: tt0       ! Freezing temperature
    real(r8) :: csmin
    real(r8) :: csmax
    real(r8) :: minrefl
    real(r8) :: mindbz
    real(r8) :: rhmini     ! Minimum rh for ice cloud fraction > 0.
    logical :: use_hetfrz_classnuc ! option to use heterogeneous freezing
    ! others
    logical :: microp_uniform
    logical :: do_cldice
    integer :: pcols
    integer :: pver
    integer :: ncol
    integer :: top_lev
    real(r8) :: deltatin
    ! in
    integer(8) :: tn
    integer(8) :: qn
    integer(8) :: relvar
    integer(8) :: accre_enhan
    integer(8) :: p
    integer(8) :: pdel
    integer(8) :: cldn
    integer(8) :: icecldf
    integer(8) :: liqcldf
    integer(8) :: naai
    integer(8) :: npccnin
    integer(8) :: rndst
    integer(8) :: nacon
    integer(8) :: tnd_qsnow
    integer(8) :: tnd_nsnow
    integer(8) :: re_ice
    integer(8) :: frzimm
    integer(8) :: frzcnt
    integer(8) :: frzdep   
    ! out
    integer(8) :: rate1ord_cw2pr_st
    integer(8) :: tlat
    integer(8) :: qvlat
    integer(8) :: qctend
    integer(8) :: qitend
    integer(8) :: nctend
    integer(8) :: nitend
    integer(8) :: effc
    integer(8) :: effc_fn
    integer(8) :: effi
    integer(8) :: prect
    integer(8) :: preci
    integer(8) :: nevapr
    integer(8) :: evapsnow
    integer(8) :: am_evp_st
    integer(8) :: prain
    integer(8) :: prodsnow
    integer(8) :: cmeout
    integer(8) :: deffi
    integer(8) :: pgamrad
    integer(8) :: lamcrad
    integer(8) :: qsout
    integer(8) :: dsout
    integer(8) :: rflx
    integer(8) :: sflx
    integer(8) :: qrout
    integer(8) :: qcsevap
    integer(8) :: qisevap
    integer(8) :: qvres
    integer(8) :: cmeiout
    integer(8) :: vtrmc
    integer(8) :: vtrmi
    integer(8) :: qcsedten
    integer(8) :: qisedten
    integer(8) :: prao
    integer(8) :: prco
    integer(8) :: mnuccco
    integer(8) :: mnuccto
    integer(8) :: msacwio
    integer(8) :: psacwso
    integer(8) :: bergso
    integer(8) :: bergo
    integer(8) :: melto
    integer(8) :: homoo
    integer(8) :: qcreso
    integer(8) :: prcio
    integer(8) :: praio
    integer(8) :: qireso
    integer(8) :: mnuccro
    integer(8) :: pracso 
    integer(8) :: meltsdt
    integer(8) :: frzrdt 
    integer(8) :: mnuccdo
    integer(8) :: nrout
    integer(8) :: nsout
    integer(8) :: refl
    integer(8) :: arefl
    integer(8) :: areflz
    integer(8) :: frefl
    integer(8) :: csrfl
    integer(8) :: acsrfl
    integer(8) :: fcsrfl
    integer(8) :: rercld
    integer(8) :: ncai
    integer(8) :: ncal
    integer(8) :: qrout2
    integer(8) :: qsout2
    integer(8) :: nrout2
    integer(8) :: nsout2
    integer(8) :: drout2
    integer(8) :: dsout2
    integer(8) :: freqs
    integer(8) :: freqr
    integer(8) :: nfice
    integer(8) :: prer_evap
    ! inout
    integer(8) :: qc
    integer(8) :: qi
    integer(8) :: nc
    integer(8) :: ni
    integer(8) :: reff_rain
    integer(8) :: reff_snow                
end type

type(micro_mg1_0_args) :: micro_mg1_0_para
integer, external :: slave_micro_mg1_0_parallel
integer, external :: slave_log10_parallel
integer, external :: slave_gamma_parallel
integer, external :: slave_pow_parallel
integer, external :: slave_exp_parallel
real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5
real(8) :: tmp_array(2)

integer :: ierr, rank

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!call set_consist_math_on()
call mpi_comm_rank(mpi_comm_world, rank, ierr)
! Return error message
errstring = ' '

if (.not. (do_cldice .or. &
     (associated(tnd_qsnow) .and. associated(tnd_nsnow) .and. associated(re_ice)))) then
   errstring = "MG's native cloud ice processes are disabled, but &
        &no replacement values were passed in."
end if

if (use_hetfrz_classnuc .and. (.not. &
     (associated(frzimm) .and. associated(frzcnt) .and. associated(frzdep)))) then
   errstring = "Hoose heterogeneous freezing is enabled, but the &
        &required tendencies were not all passed in."
end if

call phys_getopts(do_clubb_sgs_out = do_clubb_sgs)

!qrout2(:,:)=0._r8
!qsout2(:,:)=0._r8
!nrout2(:,:)=0._r8
!nsout2(:,:)=0._r8
!drout2(:,:)=0._r8
!dsout2(:,:)=0._r8
!pgam(:) = 0._r8
!pra(:) = 0.d0
!pre(:) = 0.d0
!lammin=0.d0
!lammax=0.d0
!

micro_mg1_0_para%svp_tmelt           = svp_tmelt
micro_mg1_0_para%svp_h2otrip         = svp_h2otrip
micro_mg1_0_para%svp_tboil           = svp_tboil
micro_mg1_0_para%svp_ttrice          = svp_ttrice
micro_mg1_0_para%svp_epsilo          = svp_epsilo
micro_mg1_0_para%g                   = g
micro_mg1_0_para%r                   = r                       
micro_mg1_0_para%rv                  = rv                      
micro_mg1_0_para%cpp                 = cpp                     
micro_mg1_0_para%rhow                = rhow                    
micro_mg1_0_para%tmelt               = tmelt                   
micro_mg1_0_para%xxlv                = xxlv                    
micro_mg1_0_para%xlf                 = xlf                     
micro_mg1_0_para%xxls                = xxls                    
micro_mg1_0_para%rhosn               = rhosn                   
micro_mg1_0_para%rhoi                = rhoi                    
micro_mg1_0_para%ac                  = ac                   
micro_mg1_0_para%bc                  = bc                   
micro_mg1_0_para%as                  = as                   
micro_mg1_0_para%bs                  = bs                   
micro_mg1_0_para%ai                  = ai                   
micro_mg1_0_para%bi                  = bi                   
micro_mg1_0_para%ar                  = ar                   
micro_mg1_0_para%br                  = br                   
micro_mg1_0_para%ci                  = ci                   
micro_mg1_0_para%di                  = di                   
micro_mg1_0_para%cs                  = cs                   
micro_mg1_0_para%ds                  = ds                   
micro_mg1_0_para%cr                  = cr                   
micro_mg1_0_para%dr                  = dr                   
micro_mg1_0_para%f1s                 = f1s                   
micro_mg1_0_para%f2s                 = f2s                   
micro_mg1_0_para%Eii                 = Eii                   
micro_mg1_0_para%Ecr                 = Ecr                   
micro_mg1_0_para%f1r                 = f1r                   
micro_mg1_0_para%f2r                 = f2r                   
micro_mg1_0_para%DCS                 = DCS                   
micro_mg1_0_para%qsmall              = qsmall                   
micro_mg1_0_para%bimm                = bimm                   
micro_mg1_0_para%aimm                = aimm                   
micro_mg1_0_para%rhosu               = rhosu                   
micro_mg1_0_para%mi0                 = mi0                          
micro_mg1_0_para%rin                 = rin                          
micro_mg1_0_para%pi                  = pi                          
micro_mg1_0_para%cons1               = cons1                        
micro_mg1_0_para%cons4               = cons4                        
micro_mg1_0_para%cons5               = cons5                        
micro_mg1_0_para%cons6               = cons6                        
micro_mg1_0_para%cons7               = cons7                        
micro_mg1_0_para%cons8               = cons8                        
micro_mg1_0_para%cons11              = cons11                       
micro_mg1_0_para%cons13              = cons13                       
micro_mg1_0_para%cons14              = cons14                       
micro_mg1_0_para%cons16              = cons16                       
micro_mg1_0_para%cons17              = cons17                       
micro_mg1_0_para%cons22              = cons22                       
micro_mg1_0_para%cons23              = cons23                       
micro_mg1_0_para%cons24              = cons24                       
micro_mg1_0_para%cons25              = cons25                       
micro_mg1_0_para%cons27              = cons27                       
micro_mg1_0_para%cons28              = cons28                       
micro_mg1_0_para%lammini             = lammini                      
micro_mg1_0_para%lammaxi             = lammaxi                      
micro_mg1_0_para%lamminr             = lamminr                      
micro_mg1_0_para%lammaxr             = lammaxr                      
micro_mg1_0_para%lammins             = lammins                      
micro_mg1_0_para%lammaxs             = lammaxs                      
micro_mg1_0_para%tmax_fsnow          = tmax_fsnow                   
micro_mg1_0_para%tmin_fsnow          = tmin_fsnow                   
micro_mg1_0_para%tt0                 = tt0                          
micro_mg1_0_para%csmin               = csmin                   
micro_mg1_0_para%csmax               = csmax                   
micro_mg1_0_para%minrefl             = minrefl                   
micro_mg1_0_para%mindbz              = mindbz                   
micro_mg1_0_para%rhmini              = rhmini                   
micro_mg1_0_para%use_hetfrz_classnuc = use_hetfrz_classnuc                   

micro_mg1_0_para%microp_uniform      = microp_uniform
micro_mg1_0_para%do_cldice           = do_cldice
micro_mg1_0_para%pcols               = pcols
micro_mg1_0_para%pver                = pver
micro_mg1_0_para%ncol                = ncol
micro_mg1_0_para%top_lev             = top_lev
micro_mg1_0_para%deltatin             = deltatin

micro_mg1_0_para%tn                  = loc(tn_permute                (:,:))
micro_mg1_0_para%qn                  = loc(qn_permute                (:,:))
micro_mg1_0_para%relvar              = loc(relvar_permute            (:,:))
micro_mg1_0_para%accre_enhan         = loc(accre_enhan_permute       (:,:))
micro_mg1_0_para%p                   = loc(p_permute                 (:,:))
micro_mg1_0_para%pdel                = loc(pdel_permute              (:,:))
micro_mg1_0_para%cldn                = loc(cldn_permute              (:,:))
micro_mg1_0_para%icecldf             = loc(icecldf_permute           (:,:))
micro_mg1_0_para%liqcldf             = loc(liqcldf_permute           (:,:))
micro_mg1_0_para%naai                = loc(naai_permute              (:,:))
micro_mg1_0_para%npccnin             = loc(npccnin_permute           (:,:))
micro_mg1_0_para%rndst               = loc(rndst_permute             (:,:,:))
micro_mg1_0_para%nacon               = loc(nacon_permute             (:,:,:))
micro_mg1_0_para%tnd_qsnow           = loc(tnd_qsnow_permute         (:,:))
micro_mg1_0_para%tnd_nsnow           = loc(tnd_nsnow_permute         (:,:))
micro_mg1_0_para%re_ice              = loc(re_ice_permute            (:,:))
micro_mg1_0_para%frzimm              = loc(frzimm_permute            (:,:))
micro_mg1_0_para%frzcnt              = loc(frzcnt_permute            (:,:))
micro_mg1_0_para%frzdep              = loc(frzdep_permute            (:,:))
                                                             
micro_mg1_0_para%rate1ord_cw2pr_st   = loc(rate1ord_cw2pr_st_permute (:,:))
micro_mg1_0_para%tlat                = loc(tlat_permute              (:,:))
micro_mg1_0_para%qvlat               = loc(qvlat_permute             (:,:))
micro_mg1_0_para%qctend              = loc(qctend_permute            (:,:))
micro_mg1_0_para%qitend              = loc(qitend_permute            (:,:))
micro_mg1_0_para%nctend              = loc(nctend_permute            (:,:))
micro_mg1_0_para%nitend              = loc(nitend_permute            (:,:))
micro_mg1_0_para%effc                = loc(effc_permute              (:,:))
micro_mg1_0_para%effc_fn             = loc(effc_fn_permute           (:,:))
micro_mg1_0_para%effi                = loc(effi_permute              (:,:))
micro_mg1_0_para%prect               = loc(prect_permute             (:))
micro_mg1_0_para%preci               = loc(preci_permute             (:))
micro_mg1_0_para%nevapr              = loc(nevapr_permute            (:,:))
micro_mg1_0_para%evapsnow            = loc(evapsnow_permute          (:,:))
micro_mg1_0_para%am_evp_st           = loc(am_evp_st_permute         (:,:))
micro_mg1_0_para%prain               = loc(prain_permute             (:,:))
micro_mg1_0_para%prodsnow            = loc(prodsnow_permute          (:,:))
micro_mg1_0_para%cmeout              = loc(cmeout_permute            (:,:))
micro_mg1_0_para%deffi               = loc(deffi_permute             (:,:))
micro_mg1_0_para%pgamrad             = loc(pgamrad_permute           (:,:))
micro_mg1_0_para%lamcrad             = loc(lamcrad_permute           (:,:))
micro_mg1_0_para%qsout               = loc(qsout_permute             (:,:))
micro_mg1_0_para%dsout               = loc(dsout_permute             (:,:))
micro_mg1_0_para%rflx                = loc(rflx_permute              (:,:))
micro_mg1_0_para%sflx                = loc(sflx_permute              (:,:))
micro_mg1_0_para%qrout               = loc(qrout_permute             (:,:))
micro_mg1_0_para%qcsevap             = loc(qcsevap_permute           (:,:))
micro_mg1_0_para%qisevap             = loc(qisevap_permute           (:,:))
micro_mg1_0_para%qvres               = loc(qvres_permute             (:,:))
micro_mg1_0_para%cmeiout             = loc(cmeiout_permute           (:,:))
micro_mg1_0_para%vtrmc               = loc(vtrmc_permute             (:,:))
micro_mg1_0_para%vtrmi               = loc(vtrmi_permute             (:,:))
micro_mg1_0_para%qcsedten            = loc(qcsedten_permute          (:,:))
micro_mg1_0_para%qisedten            = loc(qisedten_permute          (:,:))
micro_mg1_0_para%prao                = loc(prao_permute              (:,:))
micro_mg1_0_para%prco                = loc(prco_permute              (:,:))
micro_mg1_0_para%mnuccco             = loc(mnuccco_permute           (:,:))
micro_mg1_0_para%mnuccto             = loc(mnuccto_permute           (:,:))
micro_mg1_0_para%msacwio             = loc(msacwio_permute           (:,:))
micro_mg1_0_para%psacwso             = loc(psacwso_permute           (:,:))
micro_mg1_0_para%bergso              = loc(bergso_permute            (:,:))
micro_mg1_0_para%bergo               = loc(bergo_permute             (:,:))
micro_mg1_0_para%melto               = loc(melto_permute             (:,:))
micro_mg1_0_para%homoo               = loc(homoo_permute             (:,:))
micro_mg1_0_para%qcreso              = loc(qcreso_permute            (:,:))
micro_mg1_0_para%prcio               = loc(prcio_permute             (:,:))
micro_mg1_0_para%praio               = loc(praio_permute             (:,:))
micro_mg1_0_para%qireso              = loc(qireso_permute            (:,:))
micro_mg1_0_para%mnuccro             = loc(mnuccro_permute           (:,:))
micro_mg1_0_para%pracso              = loc(pracso_permute            (:,:))
micro_mg1_0_para%meltsdt             = loc(meltsdt_permute           (:,:))
micro_mg1_0_para%frzrdt              = loc(frzrdt_permute            (:,:))
micro_mg1_0_para%mnuccdo             = loc(mnuccdo_permute           (:,:))
micro_mg1_0_para%nrout               = loc(nrout_permute             (:,:))
micro_mg1_0_para%nsout               = loc(nsout_permute             (:,:))
micro_mg1_0_para%refl                = loc(refl_permute              (:,:))
micro_mg1_0_para%arefl               = loc(arefl_permute             (:,:))
micro_mg1_0_para%areflz              = loc(areflz_permute            (:,:))
micro_mg1_0_para%frefl               = loc(frefl_permute             (:,:))
micro_mg1_0_para%csrfl               = loc(csrfl_permute             (:,:))
micro_mg1_0_para%acsrfl              = loc(acsrfl_permute            (:,:))
micro_mg1_0_para%fcsrfl              = loc(fcsrfl_permute            (:,:))
micro_mg1_0_para%rercld              = loc(rercld_permute            (:,:))
micro_mg1_0_para%ncai                = loc(ncai_permute              (:,:))
micro_mg1_0_para%ncal                = loc(ncal_permute              (:,:))
micro_mg1_0_para%qrout2              = loc(qrout2_permute            (:,:))
micro_mg1_0_para%qsout2              = loc(qsout2_permute            (:,:))
micro_mg1_0_para%nrout2              = loc(nrout2_permute            (:,:))
micro_mg1_0_para%nsout2              = loc(nsout2_permute            (:,:))
micro_mg1_0_para%drout2              = loc(drout2_permute            (:,:))
micro_mg1_0_para%dsout2              = loc(dsout2_permute            (:,:))
micro_mg1_0_para%freqs               = loc(freqs_permute             (:,:))
micro_mg1_0_para%freqr               = loc(freqr_permute             (:,:))
micro_mg1_0_para%nfice               = loc(nfice_permute             (:,:))
micro_mg1_0_para%prer_evap           = loc(prer_evap_permute         (:,:))
                                                             
micro_mg1_0_para%qc                  = loc(qc_permute                (:,:))
micro_mg1_0_para%qi                  = loc(qi_permute                (:,:))
micro_mg1_0_para%nc                  = loc(nc_permute                (:,:))
micro_mg1_0_para%ni                  = loc(ni_permute                (:,:))
micro_mg1_0_para%reff_rain           = loc(reff_rain_permute         (:,:))
micro_mg1_0_para%reff_snow           = loc(reff_snow_permute         (:,:))

call t_startf("Transpose1")
call acc_transpose_c(tn         , tn_permute         , pcols, pver) 
call acc_transpose_c(qn         , qn_permute         , pcols, pver) 
call acc_transpose_c(relvar     , relvar_permute     , pcols, pver) 
call acc_transpose_c(accre_enhan, accre_enhan_permute, pcols, pver) 
call acc_transpose_c(qc         , qc_permute         , pcols, pver) 
call acc_transpose_c(qi         , qi_permute         , pcols, pver) 
call acc_transpose_c(nc         , nc_permute         , pcols, pver) 
call acc_transpose_c(ni         , ni_permute         , pcols, pver) 
call acc_transpose_c(p          , p_permute          , pcols, pver) 
call acc_transpose_c(pdel       , pdel_permute       , pcols, pver) 
call acc_transpose_c(cldn       , cldn_permute       , pcols, pver) 
call acc_transpose_c(icecldf    , icecldf_permute    , pcols, pver) 
call acc_transpose_c(liqcldf    , liqcldf_permute    , pcols, pver) 
call acc_transpose_c(naai       , naai_permute       , pcols, pver) 
call acc_transpose_c(npccnin    , npccnin_permute    , pcols, pver) 
call acc_transpose_c(reff_rain  , reff_rain_permute  , pcols, pver) 
call acc_transpose_c(reff_snow  , reff_snow_permute  , pcols, pver) 
if (associated(tnd_qsnow)) then 
    call acc_transpose_c(tnd_qsnow  , tnd_qsnow_permute, pcols, pver) 
endif
if (associated(tnd_nsnow)) then 
    call acc_transpose_c(tnd_nsnow  , tnd_nsnow_permute, pcols, pver) 
endif
if (associated(re_ice)) then 
    call acc_transpose_c(re_ice     , re_ice_permute, pcols, pver) 
endif
if (associated(frzimm)) then 
    call acc_transpose_c(frzimm     , frzimm_permute, pcols, pver) 
endif
if (associated(frzcnt)) then 
    call acc_transpose_c(frzcnt     , frzcnt_permute, pcols, pver) 
endif
if (associated(frzdep)) then 
    call acc_transpose_c(frzdep     , frzdep_permute, pcols, pver) 
endif
call acc_transpose_c_v2(rndst  , rndst_permute  , pcols, pver, 4) 
call acc_transpose_c_v2(nacon  , nacon_permute  , pcols, pver, 4) 
!do i=1,pcols
!    do k=1,pver
!        rndst_permute(k, i, :)           =    rndst(i, k, :)
!        nacon_permute(k, i, :)           =    nacon(i, k, :)
!    end do
!end do

!do i=1,pcols
!    !prect_permute(i)               =    prect(i) 
!    !preci_permute(i)               =    preci(i) 
!    do k=1,pver
!
!        tn_permute(k,i)                  =    tn(i,k)                 
!        qn_permute(k,i)                  =    qn(i,k) 
!        relvar_permute(k,i)              =    relvar(i,k) 
!        accre_enhan_permute(k,i)         =    accre_enhan(i,k) 
!                                             
!        qc_permute(k,i)                  =    qc(i,k) 
!        qi_permute(k,i)                  =    qi(i,k) 
!        nc_permute(k,i)                  =    nc(i,k) 
!        ni_permute(k,i)                  =    ni(i,k) 
!        p_permute(k,i)                   =    p(i,k) 
!        pdel_permute(k,i)                =    pdel(i,k) 
!        cldn_permute(k,i)                =    cldn(i,k) 
!        icecldf_permute(k,i)             =    icecldf(i,k) 
!        liqcldf_permute(k,i)             =    liqcldf(i,k) 
!                                             
!        naai_permute(k,i)                =    naai(i,k) 
!        npccnin_permute(k,i)             =    npccnin(i,k) 
!        reff_rain_permute(k,i)           =    reff_rain(i,k) 
!        reff_snow_permute(k,i)           =    reff_snow(i,k) 
!                                             
!        rndst_permute(k, i, :)           =    rndst(i, k, :)
!        nacon_permute(k, i, :)           =    nacon(i, k, :)
!
!        !rate1ord_cw2pr_st_permute(k,i)   =    rate1ord_cw2pr_st(i,k) 
!        !tlat_permute(k,i)                =    tlat(i,k) 
!        !qvlat_permute(k,i)               =    qvlat(i,k) 
!        !qctend_permute(k,i)              =    qctend(i,k) 
!        !qitend_permute(k,i)              =    qitend(i,k) 
!        !nctend_permute(k,i)              =    nctend(i,k) 
!        !nitend_permute(k,i)              =    nitend(i,k) 
!        !                                     
!        !effc_permute(k,i)                =    effc(i,k) 
!        !effc_fn_permute(k,i)             =    effc_fn(i,k) 
!        !effi_permute(k,i)                =    effi(i,k) 
!        !nevapr_permute(k,i)              =    nevapr(i,k) 
!        !evapsnow_permute(k,i)            =    evapsnow(i,k) 
!        !am_evp_st_permute(k,i)           =    am_evp_st(i,k) 
!        !prain_permute(k,i)               =    prain(i,k) 
!        !prodsnow_permute(k,i)            =    prodsnow(i,k) 
!        !cmeout_permute(k,i)              =    cmeout(i,k) 
!        !deffi_permute(k,i)               =    deffi(i,k) 
!        !pgamrad_permute(k,i)             =    pgamrad(i,k) 
!        !lamcrad_permute(k,i)             =    lamcrad(i,k) 
!        !qsout_permute(k,i)               =    qsout(i,k) 
!        !dsout_permute(k,i)               =    dsout(i,k) 
!        !                                     
!        !qrout_permute(k,i)               =    qrout(i,k) 
!        !qcsevap_permute(k,i)             =    qcsevap(i,k) 
!        !qisevap_permute(k,i)             =    qisevap(i,k) 
!        !qvres_permute(k,i)               =    qvres(i,k) 
!        !cmeiout_permute(k,i)             =    cmeiout(i,k) 
!        !vtrmc_permute(k,i)               =    vtrmc(i,k) 
!        !vtrmi_permute(k,i)               =    vtrmi(i,k) 
!        !qcsedten_permute(k,i)            =    qcsedten(i,k) 
!        !qisedten_permute(k,i)            =    qisedten(i,k) 
!        !prao_permute(k,i)                =    prao(i,k) 
!        !prco_permute(k,i)                =    prco(i,k) 
!        !mnuccco_permute(k,i)             =    mnuccco(i,k) 
!        !mnuccto_permute(k,i)             =    mnuccto(i,k) 
!        !msacwio_permute(k,i)             =    msacwio(i,k) 
!        !psacwso_permute(k,i)             =    psacwso(i,k) 
!        !bergso_permute(k,i)              =    bergso(i,k) 
!        !bergo_permute(k,i)               =    bergo(i,k) 
!        !melto_permute(k,i)               =    melto(i,k) 
!        !homoo_permute(k,i)               =    homoo(i,k) 
!        !qcreso_permute(k,i)              =    qcreso(i,k) 
!        !prcio_permute(k,i)               =    prcio(i,k) 
!        !praio_permute(k,i)               =    praio(i,k) 
!        !qireso_permute(k,i)              =    qireso(i,k) 
!        !mnuccro_permute(k,i)             =    mnuccro(i,k) 
!        !pracso_permute(k,i)              =    pracso(i,k) 
!        !meltsdt_permute(k,i)             =    meltsdt(i,k) 
!        !frzrdt_permute(k,i)              =    frzrdt(i,k) 
!        !mnuccdo_permute(k,i)             =    mnuccdo(i,k) 
!        !nrout_permute(k,i)               =    nrout(i,k) 
!        !nsout_permute(k,i)               =    nsout(i,k) 
!        !refl_permute(k,i)                =    refl(i,k) 
!        !arefl_permute(k,i)               =    arefl(i,k) 
!        !areflz_permute(k,i)              =    areflz(i,k) 
!        !frefl_permute(k,i)               =    frefl(i,k) 
!        !csrfl_permute(k,i)               =    csrfl(i,k) 
!        !acsrfl_permute(k,i)              =    acsrfl(i,k) 
!        !fcsrfl_permute(k,i)              =    fcsrfl(i,k) 
!        !rercld_permute(k,i)              =    rercld(i,k) 
!        !ncai_permute(k,i)                =    ncai(i,k) 
!        !ncal_permute(k,i)                =    ncal(i,k) 
!        !qrout2_permute(k,i)              =    qrout2(i,k) 
!        !qsout2_permute(k,i)              =    qsout2(i,k) 
!        !nrout2_permute(k,i)              =    nrout2(i,k) 
!        !nsout2_permute(k,i)              =    nsout2(i,k) 
!        !drout2_permute(k,i)              =    drout2(i,k) 
!        !dsout2_permute(k,i)              =    dsout2(i,k) 
!        !freqs_permute(k,i)               =    freqs(i,k) 
!        !freqr_permute(k,i)               =    freqr(i,k) 
!        !nfice_permute(k,i)               =    nfice(i,k) 
!        !prer_evap_permute(k,i)           =    prer_evap(i,k) 
!        !nevapr2_permute(k,i)             =    nevapr2(i,k) 
!        !                                    
!        !t1_permute(k,i)                  =    t1(i,k) 
!        !q1_permute(k,i)                  =    q1(i,k) 
!        !qc1_permute(k,i)                 =    qc1(i,k) 
!        !qi1_permute(k,i)                 =    qi1(i,k) 
!        !nc1_permute(k,i)                 =    nc1(i,k) 
!        !ni1_permute(k,i)                 =    ni1(i,k) 
!        !tlat1_permute(k,i)               =    tlat1(i,k) 
!        !qvlat1_permute(k,i)              =    qvlat1(i,k) 
!        !qctend1_permute(k,i)             =    qctend1(i,k) 
!        !qitend1_permute(k,i)             =    qitend1(i,k) 
!        !nctend1_permute(k,i)             =    nctend1(i,k) 
!        !nitend1_permute(k,i)             =    nitend1(i,k) 
!        !q_permute(k,i)                   =    q(i,k) 
!        !t_permute(k,i)                   =    t(i,k) 
!        !rho_permute(k,i)                 =    rho(i,k) 
!        !dv_permute(k,i)                  =    dv(i,k) 
!        !mu_permute(k,i)                  =    mu(i,k) 
!        !sc_permute(k,i)                  =    sc(i,k) 
!        !kap_permute(k,i)                 =    kap(i,k) 
!        !rhof_permute(k,i)                =    rhof(i,k) 
!        !cldmax_permute(k,i)              =    cldmax(i,k) 
!        !cldm_permute(k,i)                =    cldm(i,k) 
!        !icldm_permute(k,i)               =    icldm(i,k) 
!        !lcldm_permute(k,i)               =    lcldm(i,k) 
!        !cme_permute(k,i)                 =    cme(i,k) 
!        !cmei_permute(k,i)                =    cmei(i,k) 
!        !cwml_permute(k,i)                =    cwml(i,k) 
!        !cwmi_permute(k,i)                =    cwmi(i,k) 
!        !lcldn_permute(k,i)               =    lcldn(i,k) 
!        !lcldo_permute(k,i)               =    lcldo(i,k) 
!        !nctend_mixnuc_permute(k,i)       =    nctend_mixnuc(i,k) 
!        !                                    
!        !qcic_permute(k,i)                =    qcic(i,k) 
!        !qiic_permute(k,i)                =    qiic(i,k) 
!        !qniic_permute(k,i)               =    qniic(i,k) 
!        !qric_permute(k,i)                =    qric(i,k) 
!        !ncic_permute(k,i)                =    ncic(i,k) 
!        !niic_permute(k,i)                =    niic(i,k) 
!        !nsic_permute(k,i)                =    nsic(i,k) 
!        !nric_permute(k,i)                =    nric(i,k) 
!        !arcld_permute(k,i)               =    arcld(i,k) 
!        !dumc_permute(k,i)                =    dumc(i,k) 
!        !dumnc_permute(k,i)               =    dumnc(i,k) 
!        !dumi_permute(k,i)                =    dumi(i,k) 
!        !dumni_permute(k,i)               =    dumni(i,k) 
!        !dums_permute(k,i)                =    dums(i,k) 
!        !dumns_permute(k,i)               =    dumns(i,k) 
!        !dumr_permute(k,i)                =    dumr(i,k) 
!        !dumnr_permute(k,i)               =    dumnr(i,k) 
!        !relhum_permute(k,i)              =    relhum(i,k) 
!        !arn_permute(k,i)                 =    arn(i,k) 
!        !asn_permute(k,i)                 =    asn(i,k) 
!        !acn_permute(k,i)                 =    acn(i,k) 
!        !ain_permute(k,i)                 =    ain(i,k) 
!        !dz_permute(k,i)                  =    dz(i,k) 
!        !                                    
!        !tsp_permute(k,i)                 =    tsp(i,k) 
!        !qsp_permute(k,i)                 =    qsp(i,k) 
!        !qsphy_permute(k,i)               =    qsphy(i,k) 
!        !esl_permute(k,i)                 =    esl(i,k) 
!        !esi_permute(k,i)                 =    esi(i,k) 
!        !qnitend_permute(k,i)             =    qnitend(i,k) 
!        !nstend_permute(k,i)              =    nstend(i,k) 
!        !qrtend_permute(k,i)              =    qrtend(i,k) 
!        !nrtend_permute(k,i)              =    nrtend(i,k) 
!        !berg_permute(k,i)                =    berg(i,k) 
!        !drout_permute(k,i)               =    drout(i,k) 
!        !dum2i_permute(k,i)               =    dum2i(i,k) 
!        !dum2l_permute(k,i)               =    dum2l(i,k) 
!        !cldmw_permute(k,i)               =    cldmw(i,k) 
!        !rainrt_permute(k,i)              =    rainrt(i,k) 
!        !rainrt1_permute(k,i)             =    rainrt1(i,k)                
!
!    enddo
!enddo
!!do i=1,pcols
!!    do k=1,pver+1
!!        rflx1_permute(k,i) = rflx1(i,k)
!!        sflx1_permute(k,i) = sflx1(i,k)
!!        rflx_permute(k,i) = rflx(i,k)
!!        sflx_permute(k,i) = sflx(i,k)
!!    enddo
!!enddo
!
!if (associated(tnd_qsnow)) then 
!    do i = 1, pcols
!        do k = 1, pver
!            tnd_qsnow_permute(k,i)           =    tnd_qsnow(i,k) 
!        enddo
!    enddo 
!endif
!if (associated(tnd_nsnow)) then 
!    do i = 1, pcols
!        do k = 1, pver
!            tnd_nsnow_permute(k,i)           =    tnd_nsnow(i,k) 
!        enddo
!    enddo 
!endif
!if (associated(re_ice)) then 
!    do i = 1, pcols
!        do k = 1, pver
!            re_ice_permute(k,i)              =    re_ice(i,k) 
!        enddo
!    enddo 
!endif
!if (associated(frzimm)) then 
!    do i = 1, pcols
!        do k = 1, pver
!            frzimm_permute(k,i)              =    frzimm(i,k) 
!        enddo
!    enddo 
!endif
!if (associated(frzcnt)) then 
!    do i = 1, pcols
!        do k = 1, pver
!            frzcnt_permute(k,i)              =    frzcnt(i,k) 
!        enddo
!    enddo 
!endif
!if (associated(frzdep)) then 
!    do i = 1, pcols
!        do k = 1, pver
!            frzdep_permute(k,i)              =    frzdep(i,k) 
!        enddo
!    enddo 
!endif
call t_stopf("Transpose1")

call t_startf("mg_tend_computing")
!call lwpf_init_fortran()
call athread_spawn(slave_micro_mg1_0_parallel, micro_mg1_0_para)
call athread_join()
!call lwpf_report_fortran()
call t_stopf("mg_tend_computing")

call t_startf("Transpose2")
call acc_transpose_c(qc_permute               , qc               , pver, pcols) 
call acc_transpose_c(qi_permute               , qi               , pver, pcols) 
call acc_transpose_c(nc_permute               , nc               , pver, pcols) 
call acc_transpose_c(ni_permute               , ni               , pver, pcols) 
call acc_transpose_c(reff_rain_permute        , reff_rain        , pver, pcols) 
call acc_transpose_c(reff_snow_permute        , reff_snow        , pver, pcols) 
call acc_transpose_c(rate1ord_cw2pr_st_permute, rate1ord_cw2pr_st, pver, pcols) 
call acc_transpose_c(tlat_permute             , tlat             , pver, pcols) 
call acc_transpose_c(qvlat_permute            , qvlat            , pver, pcols) 
call acc_transpose_c(qctend_permute           , qctend           , pver, pcols) 
call acc_transpose_c(qitend_permute           , qitend           , pver, pcols) 
call acc_transpose_c(nctend_permute           , nctend           , pver, pcols) 
call acc_transpose_c(nitend_permute           , nitend           , pver, pcols) 
call acc_transpose_c(effc_permute             , effc             , pver, pcols) 
call acc_transpose_c(effc_fn_permute          , effc_fn          , pver, pcols) 
call acc_transpose_c(effi_permute             , effi             , pver, pcols) 
call acc_transpose_c(nevapr_permute           , nevapr           , pver, pcols) 
call acc_transpose_c(evapsnow_permute         , evapsnow         , pver, pcols) 
call acc_transpose_c(am_evp_st_permute        , am_evp_st        , pver, pcols) 
call acc_transpose_c(prain_permute            , prain            , pver, pcols) 
call acc_transpose_c(prodsnow_permute         , prodsnow         , pver, pcols) 
call acc_transpose_c(cmeout_permute           , cmeout           , pver, pcols) 
call acc_transpose_c(deffi_permute            , deffi            , pver, pcols) 
call acc_transpose_c(pgamrad_permute          , pgamrad          , pver, pcols) 
call acc_transpose_c(lamcrad_permute          , lamcrad          , pver, pcols) 
call acc_transpose_c(qsout_permute            , qsout            , pver, pcols) 
call acc_transpose_c(dsout_permute            , dsout            , pver, pcols) 
call acc_transpose_c(qrout_permute            , qrout            , pver, pcols) 
call acc_transpose_c(qcsevap_permute          , qcsevap          , pver, pcols) 
call acc_transpose_c(qisevap_permute          , qisevap          , pver, pcols) 
call acc_transpose_c(qvres_permute            , qvres            , pver, pcols) 
call acc_transpose_c(cmeiout_permute          , cmeiout          , pver, pcols) 
call acc_transpose_c(vtrmc_permute            , vtrmc            , pver, pcols) 
call acc_transpose_c(vtrmi_permute            , vtrmi            , pver, pcols) 
call acc_transpose_c(qcsedten_permute         , qcsedten         , pver, pcols) 
call acc_transpose_c(qisedten_permute         , qisedten         , pver, pcols) 
call acc_transpose_c(prao_permute             , prao             , pver, pcols) 
call acc_transpose_c(prco_permute             , prco             , pver, pcols) 
call acc_transpose_c(mnuccco_permute          , mnuccco          , pver, pcols) 
call acc_transpose_c(mnuccto_permute          , mnuccto          , pver, pcols) 
call acc_transpose_c(msacwio_permute          , msacwio          , pver, pcols) 
call acc_transpose_c(psacwso_permute          , psacwso          , pver, pcols) 
call acc_transpose_c(bergso_permute           , bergso           , pver, pcols) 
call acc_transpose_c(bergo_permute            , bergo            , pver, pcols) 
call acc_transpose_c(melto_permute            , melto            , pver, pcols) 
call acc_transpose_c(homoo_permute            , homoo            , pver, pcols) 
call acc_transpose_c(qcreso_permute           , qcreso           , pver, pcols) 
call acc_transpose_c(prcio_permute            , prcio            , pver, pcols) 
call acc_transpose_c(praio_permute            , praio            , pver, pcols) 
call acc_transpose_c(qireso_permute           , qireso           , pver, pcols) 
call acc_transpose_c(mnuccro_permute          , mnuccro          , pver, pcols) 
call acc_transpose_c(pracso_permute           , pracso           , pver, pcols) 
call acc_transpose_c(meltsdt_permute          , meltsdt          , pver, pcols) 
call acc_transpose_c(frzrdt_permute           , frzrdt           , pver, pcols) 
call acc_transpose_c(mnuccdo_permute          , mnuccdo          , pver, pcols) 
call acc_transpose_c(nrout_permute            , nrout            , pver, pcols) 
call acc_transpose_c(nsout_permute            , nsout            , pver, pcols) 
call acc_transpose_c(refl_permute             , refl             , pver, pcols) 
call acc_transpose_c(arefl_permute            , arefl            , pver, pcols) 
call acc_transpose_c(areflz_permute           , areflz           , pver, pcols) 
call acc_transpose_c(frefl_permute            , frefl            , pver, pcols) 
call acc_transpose_c(csrfl_permute            , csrfl            , pver, pcols) 
call acc_transpose_c(acsrfl_permute           , acsrfl           , pver, pcols) 
call acc_transpose_c(fcsrfl_permute           , fcsrfl           , pver, pcols) 
call acc_transpose_c(rercld_permute           , rercld           , pver, pcols) 
call acc_transpose_c(ncai_permute             , ncai             , pver, pcols) 
call acc_transpose_c(ncal_permute             , ncal             , pver, pcols) 
call acc_transpose_c(qrout2_permute           , qrout2           , pver, pcols) 
call acc_transpose_c(qsout2_permute           , qsout2           , pver, pcols) 
call acc_transpose_c(nrout2_permute           , nrout2           , pver, pcols) 
call acc_transpose_c(nsout2_permute           , nsout2           , pver, pcols) 
call acc_transpose_c(drout2_permute           , drout2           , pver, pcols) 
call acc_transpose_c(dsout2_permute           , dsout2           , pver, pcols) 
call acc_transpose_c(freqs_permute            , freqs            , pver, pcols) 
call acc_transpose_c(freqr_permute            , freqr            , pver, pcols) 
call acc_transpose_c(nfice_permute            , nfice            , pver, pcols) 
call acc_transpose_c(prer_evap_permute        , prer_evap        , pver, pcols) 
call acc_transpose_c(rflx_permute             , rflx             , pver+1, pcols) 
call acc_transpose_c(sflx_permute             , sflx             , pver+1, pcols) 
do i=1,pcols
    prect(i)               =    prect_permute(i) 
    preci(i)               =    preci_permute(i) 
enddo
!do i=1,pcols
!    prect(i)               =    prect_permute(i) 
!    preci(i)               =    preci_permute(i) 
!    do k=1,pver
!
!         qc(i,k)                      =    qc_permute(k,i)                
!         qi(i,k)                      =    qi_permute(k,i)                
!         nc(i,k)                      =    nc_permute(k,i)                
!         ni(i,k)                      =    ni_permute(k,i)                
!         rate1ord_cw2pr_st(i,k)       =    rate1ord_cw2pr_st_permute(k,i) 
!         tlat(i,k)                    =    tlat_permute(k,i)              
!         qvlat(i,k)                   =    qvlat_permute(k,i)             
!         qctend(i,k)                  =    qctend_permute(k,i)            
!         qitend(i,k)                  =    qitend_permute(k,i)            
!         nctend(i,k)                  =    nctend_permute(k,i)            
!         nitend(i,k)                  =    nitend_permute(k,i)            
!                                                                          
!         effc(i,k)                    =    effc_permute(k,i)              
!         effc_fn(i,k)                 =    effc_fn_permute(k,i)           
!         effi(i,k)                    =    effi_permute(k,i)              
!         nevapr(i,k)                  =    nevapr_permute(k,i)            
!         evapsnow(i,k)                =    evapsnow_permute(k,i)          
!         am_evp_st(i,k)               =    am_evp_st_permute(k,i)         
!         prain(i,k)                   =    prain_permute(k,i)             
!         prodsnow(i,k)                =    prodsnow_permute(k,i)          
!         cmeout(i,k)                  =    cmeout_permute(k,i)            
!         deffi(i,k)                   =    deffi_permute(k,i)             
!         pgamrad(i,k)                 =    pgamrad_permute(k,i)           
!         lamcrad(i,k)                 =    lamcrad_permute(k,i)           
!         qsout(i,k)                   =    qsout_permute(k,i)             
!         dsout(i,k)                   =    dsout_permute(k,i)             
!                                                                          
!         qrout(i,k)                   =    qrout_permute(k,i)             
!         reff_rain(i,k)               =    reff_rain_permute(k,i)         
!         reff_snow(i,k)               =    reff_snow_permute(k,i)         
!         qcsevap(i,k)                 =    qcsevap_permute(k,i)           
!         qisevap(i,k)                 =    qisevap_permute(k,i)           
!         qvres(i,k)                   =    qvres_permute(k,i)             
!         cmeiout(i,k)                 =    cmeiout_permute(k,i)           
!         vtrmc(i,k)                   =    vtrmc_permute(k,i)             
!         vtrmi(i,k)                   =    vtrmi_permute(k,i)             
!         qcsedten(i,k)                =    qcsedten_permute(k,i)          
!         qisedten(i,k)                =    qisedten_permute(k,i)          
!         prao(i,k)                    =    prao_permute(k,i)              
!         prco(i,k)                    =    prco_permute(k,i)              
!         mnuccco(i,k)                 =    mnuccco_permute(k,i)           
!         mnuccto(i,k)                 =    mnuccto_permute(k,i)           
!         msacwio(i,k)                 =    msacwio_permute(k,i)           
!         psacwso(i,k)                 =    psacwso_permute(k,i)           
!         bergso(i,k)                  =    bergso_permute(k,i)            
!         bergo(i,k)                   =    bergo_permute(k,i)             
!         melto(i,k)                   =    melto_permute(k,i)             
!         homoo(i,k)                   =    homoo_permute(k,i)             
!         qcreso(i,k)                  =    qcreso_permute(k,i)            
!         prcio(i,k)                   =    prcio_permute(k,i)             
!         praio(i,k)                   =    praio_permute(k,i)             
!         qireso(i,k)                  =    qireso_permute(k,i)            
!         mnuccro(i,k)                 =    mnuccro_permute(k,i)           
!         pracso(i,k)                  =    pracso_permute(k,i)            
!         meltsdt(i,k)                 =    meltsdt_permute(k,i)           
!         frzrdt(i,k)                  =    frzrdt_permute(k,i)            
!         mnuccdo(i,k)                 =    mnuccdo_permute(k,i)           
!         nrout(i,k)                   =    nrout_permute(k,i)             
!         nsout(i,k)                   =    nsout_permute(k,i)             
!         refl(i,k)                    =    refl_permute(k,i)              
!         arefl(i,k)                   =    arefl_permute(k,i)             
!         areflz(i,k)                  =    areflz_permute(k,i)            
!         frefl(i,k)                   =    frefl_permute(k,i)             
!         csrfl(i,k)                   =    csrfl_permute(k,i)             
!         acsrfl(i,k)                  =    acsrfl_permute(k,i)            
!         fcsrfl(i,k)                  =    fcsrfl_permute(k,i)            
!         rercld(i,k)                  =    rercld_permute(k,i)            
!         ncai(i,k)                    =    ncai_permute(k,i)              
!         ncal(i,k)                    =    ncal_permute(k,i)              
!         qrout2(i,k)                  =    qrout2_permute(k,i)            
!         qsout2(i,k)                  =    qsout2_permute(k,i)            
!         nrout2(i,k)                  =    nrout2_permute(k,i)            
!         nsout2(i,k)                  =    nsout2_permute(k,i)            
!         drout2(i,k)                  =    drout2_permute(k,i)            
!         dsout2(i,k)                  =    dsout2_permute(k,i)            
!         freqs(i,k)                   =    freqs_permute(k,i)             
!         freqr(i,k)                   =    freqr_permute(k,i)             
!         nfice(i,k)                   =    nfice_permute(k,i)             
!         prer_evap(i,k)               =    prer_evap_permute(k,i)         
!    enddo
!enddo
!do i=1,pcols
!    do k=1,pver+1
!        rflx(i,k)   =  rflx_permute(k,i) 
!        sflx(i,k)   =  sflx_permute(k,i) 
!    enddo
!enddo
call t_stopf("Transpose2")

!do i=1,pcols
!    prect(i)               =    0.d0
!    preci(i)               =    0.d0
!    do k=1,pver
!         !qc(i,k)                      =    0.d0
!         !qi(i,k)                      =    0.d0
!         !nc(i,k)                      =    0.d0
!         !ni(i,k)                      =    0.d0
!         rate1ord_cw2pr_st(i,k)       =    0.d0 
!         tlat(i,k)                    =    0.d0 
!         qvlat(i,k)                   =    0.d0 
!         qctend(i,k)                  =    0.d0 
!         qitend(i,k)                  =    0.d0 
!         nctend(i,k)                  =    0.d0 
!         nitend(i,k)                  =    0.d0 
!                                           
!         effc(i,k)                    =    0.d0 
!         effc_fn(i,k)                 =    0.d0 
!         effi(i,k)                    =    0.d0 
!         nevapr(i,k)                  =    0.d0 
!         evapsnow(i,k)                =    0.d0 
!         am_evp_st(i,k)               =    0.d0 
!         prain(i,k)                   =    0.d0 
!         prodsnow(i,k)                =    0.d0 
!         cmeout(i,k)                  =    0.d0 
!         deffi(i,k)                   =    0.d0 
!         pgamrad(i,k)                 =    0.d0 
!         lamcrad(i,k)                 =    0.d0 
!         qsout(i,k)                   =    0.d0 
!         dsout(i,k)                   =    0.d0 
!                                           
!         qrout(i,k)                   =    0.d0 
!         !reff_rain(i,k)               =    reff_rain_permute(k,i)
!         !reff_snow(i,k)               =    reff_rain_permute(k,i)
!         qcsevap(i,k)                 =    0.d0 
!         qisevap(i,k)                 =    0.d0 
!         qvres(i,k)                   =    0.d0 
!         cmeiout(i,k)                 =    0.d0 
!         vtrmc(i,k)                   =    0.d0 
!         vtrmi(i,k)                   =    0.d0 
!         qcsedten(i,k)                =    0.d0 
!         qisedten(i,k)                =    0.d0 
!         prao(i,k)                    =    0.d0 
!         prco(i,k)                    =    0.d0 
!         mnuccco(i,k)                 =    0.d0 
!         mnuccto(i,k)                 =    0.d0 
!         msacwio(i,k)                 =    0.d0 
!         psacwso(i,k)                 =    0.d0 
!         bergso(i,k)                  =    0.d0 
!         bergo(i,k)                   =    0.d0 
!         melto(i,k)                   =    0.d0 
!         homoo(i,k)                   =    0.d0 
!         qcreso(i,k)                  =    0.d0 
!         prcio(i,k)                   =    0.d0 
!         praio(i,k)                   =    0.d0 
!         qireso(i,k)                  =    0.d0 
!         mnuccro(i,k)                 =    0.d0 
!         pracso(i,k)                  =    0.d0 
!         meltsdt(i,k)                 =    0.d0 
!         frzrdt(i,k)                  =    0.d0 
!         mnuccdo(i,k)                 =    0.d0 
!         nrout(i,k)                   =    0.d0 
!         nsout(i,k)                   =    0.d0 
!         refl(i,k)                    =    0.d0 
!         arefl(i,k)                   =    0.d0 
!         areflz(i,k)                  =    0.d0 
!         frefl(i,k)                   =    0.d0 
!         csrfl(i,k)                   =    0.d0 
!         acsrfl(i,k)                  =    0.d0 
!         fcsrfl(i,k)                  =    0.d0 
!         rercld(i,k)                  =    0.d0 
!         ncai(i,k)                    =    0.d0 
!         ncal(i,k)                    =    0.d0 
!         qrout2(i,k)                  =    0.d0 
!         qsout2(i,k)                  =    0.d0 
!         nrout2(i,k)                  =    0.d0 
!         nsout2(i,k)                  =    0.d0 
!         drout2(i,k)                  =    0.d0 
!         dsout2(i,k)                  =    0.d0 
!         freqs(i,k)                   =    0.d0 
!         freqr(i,k)                   =    0.d0 
!         nfice(i,k)                   =    0.d0 
!         prer_evap(i,k)               =    0.d0 
!    end do
!end do
!call t_startf("Mg_computing_mpe")
!! initialize  output fields for number conc qand ice nucleation
!ncai(1:ncol,1:pver)=0._r8 
!ncal(1:ncol,1:pver)=0._r8  
!
!!Initialize rain size
!rercld(1:ncol,1:pver)=0._r8
!arcld(1:ncol,1:pver)=0._r8
!
!!initialize radiation output variables
!pgamrad(1:ncol,1:pver)=0._r8 ! liquid gamma parameter for optics (radiation)
!lamcrad(1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics (radiation)
!deffi  (1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics (radiation)
!!initialize radiation output variables
!!initialize water vapor tendency term output
!qcsevap(1:ncol,1:pver)=0._r8 
!qisevap(1:ncol,1:pver)=0._r8 
!qvres  (1:ncol,1:pver)=0._r8 
!cmeiout (1:ncol,1:pver)=0._r8
!vtrmc (1:ncol,1:pver)=0._r8
!vtrmi (1:ncol,1:pver)=0._r8
!qcsedten (1:ncol,1:pver)=0._r8
!qisedten (1:ncol,1:pver)=0._r8    
!
!prao(1:ncol,1:pver)=0._r8 
!prco(1:ncol,1:pver)=0._r8 
!mnuccco(1:ncol,1:pver)=0._r8 
!mnuccto(1:ncol,1:pver)=0._r8 
!msacwio(1:ncol,1:pver)=0._r8 
!psacwso(1:ncol,1:pver)=0._r8 
!bergso(1:ncol,1:pver)=0._r8 
!bergo(1:ncol,1:pver)=0._r8 
!melto(1:ncol,1:pver)=0._r8 
!homoo(1:ncol,1:pver)=0._r8 
!qcreso(1:ncol,1:pver)=0._r8 
!prcio(1:ncol,1:pver)=0._r8 
!praio(1:ncol,1:pver)=0._r8 
!qireso(1:ncol,1:pver)=0._r8 
!mnuccro(1:ncol,1:pver)=0._r8 
!pracso (1:ncol,1:pver)=0._r8 
!meltsdt(1:ncol,1:pver)=0._r8
!frzrdt (1:ncol,1:pver)=0._r8
!mnuccdo(1:ncol,1:pver)=0._r8
!
!rflx(:,:)=0._r8
!sflx(:,:)=0._r8
!effc(:,:)=0._r8
!effc_fn(:,:)=0._r8
!effi(:,:)=0._r8
!
!! assign variable deltat for sub-stepping...
!deltat=deltatin
!
!! parameters for scheme
!
!omsm=0.99999_r8
!dto2=0.5_r8*deltat
!mincld=0.0001_r8
!
!! initialize multi-level fields
!q(1:ncol,1:pver)=qn(1:ncol,1:pver)
!t(1:ncol,1:pver)=tn(1:ncol,1:pver)
!
!! initialize time-varying parameters
!
!do k=1,pver
!   do i=1,ncol
!      rho(i,k)=p(i,k)/(r*t(i,k))
!      dv(i,k) = 8.794E-5_r8*t(i,k)**1.81_r8/p(i,k)
!      mu(i,k) = 1.496E-6_r8*t(i,k)**1.5_r8/(t(i,k)+120._r8) 
!      sc(i,k) = mu(i,k)/(rho(i,k)*dv(i,k))
!      kap(i,k) = 1.414e3_r8*1.496e-6_r8*t(i,k)**1.5_r8/(t(i,k)+120._r8) 
!
!      ! air density adjustment for fallspeed parameters
!      ! includes air density correction factor to the
!      ! power of 0.54 following Heymsfield and Bansemer 2007
!
!      rhof(i,k)=(rhosu/rho(i,k))**0.54_r8
!
!      arn(i,k)=ar*rhof(i,k)
!      asn(i,k)=as*rhof(i,k)
!      acn(i,k)=ac*rhof(i,k)
!      ain(i,k)=ai*rhof(i,k)
!
!      ! get dz from dp and hydrostatic approx
!      ! keep dz positive (define as layer k-1 - layer k)
!
!      dz(i,k)= pdel(i,k)/(rho(i,k)*g)
!
!   end do
!end do
!
!! initialization
!qc(1:ncol,1:top_lev-1) = 0._r8
!qi(1:ncol,1:top_lev-1) = 0._r8
!nc(1:ncol,1:top_lev-1) = 0._r8
!ni(1:ncol,1:top_lev-1) = 0._r8
!t1(1:ncol,1:pver) = t(1:ncol,1:pver)
!q1(1:ncol,1:pver) = q(1:ncol,1:pver)
!qc1(1:ncol,1:pver) = qc(1:ncol,1:pver)
!qi1(1:ncol,1:pver) = qi(1:ncol,1:pver)
!nc1(1:ncol,1:pver) = nc(1:ncol,1:pver)
!ni1(1:ncol,1:pver) = ni(1:ncol,1:pver)
!
!! initialize tendencies to zero
!tlat1(1:ncol,1:pver)=0._r8
!qvlat1(1:ncol,1:pver)=0._r8
!qctend1(1:ncol,1:pver)=0._r8
!qitend1(1:ncol,1:pver)=0._r8
!nctend1(1:ncol,1:pver)=0._r8
!nitend1(1:ncol,1:pver)=0._r8
!
!! initialize precip output
!qrout(1:ncol,1:pver)=0._r8
!qsout(1:ncol,1:pver)=0._r8
!nrout(1:ncol,1:pver)=0._r8
!nsout(1:ncol,1:pver)=0._r8
!dsout(1:ncol,1:pver)=0._r8
!
!drout(1:ncol,1:pver)=0._r8
!
!reff_rain(1:ncol,1:pver)=0._r8
!reff_snow(1:ncol,1:pver)=0._r8
!
!! initialize variables for trop_mozart
!nevapr(1:ncol,1:pver) = 0._r8
!nevapr2(1:ncol,1:pver) = 0._r8
!evapsnow(1:ncol,1:pver) = 0._r8
!prain(1:ncol,1:pver) = 0._r8
!prodsnow(1:ncol,1:pver) = 0._r8
!cmeout(1:ncol,1:pver) = 0._r8
!
!am_evp_st(1:ncol,1:pver) = 0._r8
!
!! for refl calc
!rainrt1(1:ncol,1:pver) = 0._r8
!
!! initialize precip fraction and output tendencies
!cldmax(1:ncol,1:pver)=mincld
!
!!initialize aerosol number
!!        naer2(1:ncol,1:pver,:)=0._r8
!dum2l(1:ncol,1:pver)=0._r8
!dum2i(1:ncol,1:pver)=0._r8
!
!! initialize avg precip rate
!prect1(1:ncol)=0._r8
!preci1(1:ncol)=0._r8
!
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!Get humidity and saturation vapor pressures
!
!do k=top_lev,pver
!
!   do i=1,ncol
!
!      ! find wet bulk temperature and saturation value for provisional t and q without
!      ! condensation
!      
!      es(i) = svp_water(t(i,k))
!      qs(i) = svp_to_qsat(es(i), p(i,k))
!
!      ! Prevents negative values.
!      if (qs(i) < 0.0_r8) then
!         qs(i) = 1.0_r8
!         es(i) = p(i,k)
!      end if
!
!      esl(i,k)=svp_water(t(i,k))
!      esi(i,k)=svp_ice(t(i,k))
!
!      ! hm fix, make sure when above freezing that esi=esl, not active yet
!      if (t(i,k).gt.tmelt)esi(i,k)=esl(i,k)
!
!      relhum(i,k)=q(i,k)/qs(i)
!
!      ! get cloud fraction, check for minimum
!
!      cldm(i,k)=max(cldn(i,k),mincld)
!      cldmw(i,k)=max(cldn(i,k),mincld)
!
!      icldm(i,k)=max(icecldf(i,k),mincld)
!      lcldm(i,k)=max(liqcldf(i,k),mincld)
!
!      ! subcolumns, set cloud fraction variables to one
!      ! if cloud water or ice is present, if not present
!      ! set to mincld (mincld used instead of zero, to prevent
!      ! possible division by zero errors
!
!      if (microp_uniform) then
!
!         cldm(i,k)=mincld
!         cldmw(i,k)=mincld
!         icldm(i,k)=mincld
!         lcldm(i,k)=mincld
!
!         if (qc(i,k).ge.qsmall) then
!            lcldm(i,k)=1._r8           
!            cldm(i,k)=1._r8
!            cldmw(i,k)=1._r8
!         end if
!
!         if (qi(i,k).ge.qsmall) then             
!            cldm(i,k)=1._r8
!            icldm(i,k)=1._r8
!         end if
!
!      end if               ! sub-columns
!
!      ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
!
!      nfice(i,k)=0._r8
!      dumfice=qc(i,k)+qi(i,k)
!      if (dumfice.gt.qsmall .and. qi(i,k).gt.qsmall) then
!         nfice(i,k)=qi(i,k)/dumfice
!      endif
!
!      if (do_cldice .and. (t(i,k).lt.tmelt - 5._r8)) then
!
!         ! if aerosols interact with ice set number of activated ice nuclei
!         dum2=naai(i,k)
!
!         dumnnuc=(dum2-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
!         dumnnuc=max(dumnnuc,0._r8)
!         ! get provisional ni and qi after nucleation in order to calculate
!         ! Bergeron process below
!         ninew=ni(i,k)+dumnnuc*deltat
!         qinew=qi(i,k)+dumnnuc*deltat*mi0
!
!         !T>268
!      else
!         ninew=ni(i,k)
!         qinew=qi(i,k)
!      end if
!
!      ! Initialize CME components
!
!      cme(i,k) = 0._r8
!      cmei(i,k)=0._r8
!
!
!      !-------------------------------------------------------------------
!      !Bergeron process
!
!      ! make sure to initialize bergeron process to zero
!      berg(i,k)=0._r8
!      prd = 0._r8
!
!      !condensation loop.
!
!      ! get in-cloud qi and ni after nucleation
!      if (icldm(i,k) .gt. 0._r8) then 
!         qiic(i,k)=qinew/icldm(i,k)
!         niic(i,k)=ninew/icldm(i,k)
!      else
!         qiic(i,k)=0._r8
!         niic(i,k)=0._r8
!      endif
!
!      !if T < 0 C then bergeron.
!      if (do_cldice .and. (t(i,k).lt.273.15_r8)) then
!
!         !if ice exists
!         if (qi(i,k).gt.qsmall) then
!
!            bergtsf = 0._r8 ! bergeron time scale (fraction of timestep)
!
!            qvi = svp_to_qsat(esi(i,k), p(i,k))
!            qvl = svp_to_qsat(esl(i,k), p(i,k))
!
!            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
!            abi = 1._r8+dqsidt*xxls/cpp
!
!            ! get ice size distribution parameters
!
!            if (qiic(i,k).ge.qsmall) then
!               lami(k) = (cons1*ci* &
!                    niic(i,k)/qiic(i,k))**(1._r8/di)
!               n0i(k) = niic(i,k)*lami(k)
!
!               ! check for slope
!               ! adjust vars
!               if (lami(k).lt.lammini) then
!
!                  lami(k) = lammini
!                  n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!               else if (lami(k).gt.lammaxi) then
!                  lami(k) = lammaxi
!                  n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!               end if
!
!               epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))
!
!               !if liquid exists  
!               if (qc(i,k).gt. qsmall) then 
!
!                  !begin bergeron process
!                  !     do bergeron (vapor deposition with RHw=1)
!                  !     code to find berg (a rate) goes here
!
!                  ! calculate Bergeron process
!
!                  prd = epsi*(qvl-qvi)/abi
!
!               else
!                  prd = 0._r8
!               end if
!
!               ! multiply by cloud fraction
!
!               prd = prd*min(icldm(i,k),lcldm(i,k))
!
!               !     transfer of existing cloud liquid to ice
!
!               berg(i,k)=max(0._r8,prd)
!
!            end if  !end liquid exists bergeron
!
!            if (berg(i,k).gt.0._r8) then
!               bergtsf=max(0._r8,(qc(i,k)/berg(i,k))/deltat) 
!
!               if(bergtsf.lt.1._r8) berg(i,k) = max(0._r8,qc(i,k)/deltat)
!
!            endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            if (bergtsf.lt.1._r8.or.icldm(i,k).gt.lcldm(i,k)) then
!
!               if (qiic(i,k).ge.qsmall) then
!
!                  ! first case is for case when liquid water is present, but is completely depleted 
!                  ! in time step, i.e., bergrsf > 0 but < 1
!
!                  if (qc(i,k).ge.qsmall) then
!                     rhin  = (1.0_r8 + relhum(i,k)) / 2._r8
!                     if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then
!                        prd = epsi*(rhin*qvl-qvi)/abi
!
!                        ! multiply by cloud fraction assuming liquid/ice maximum overlap
!                        prd = prd*min(icldm(i,k),lcldm(i,k))
!
!                        ! add to cmei
!                        cmei(i,k) = cmei(i,k) + (prd * (1._r8- bergtsf))
!
!                     end if ! rhin 
!                  end if ! qc > qsmall
!
!                  ! second case is for pure ice cloud, either no liquid, or icldm > lcldm
!
!                  if (qc(i,k).lt.qsmall.or.icldm(i,k).gt.lcldm(i,k)) then
!
!                     ! note: for case of no liquid, need to set liquid cloud fraction to zero
!                     ! store liquid cloud fraction in 'dum'
!
!                     if (qc(i,k).lt.qsmall) then 
!                        dum=0._r8 
!                     else
!                        dum=lcldm(i,k)
!                     end if
!
!                     ! set RH to grid-mean value for pure ice cloud
!                     rhin = relhum(i,k)
!
!                     if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then
!
!                        prd = epsi*(rhin*qvl-qvi)/abi
!
!                        ! multiply by relevant cloud fraction for pure ice cloud
!                        ! assuming maximum overlap of liquid/ice
!                        prd = prd*max((icldm(i,k)-dum),0._r8)
!                        cmei(i,k) = cmei(i,k) + prd
!
!                     end if ! rhin
!                  end if ! qc or icldm > lcldm
!               end if ! qiic
!            end if ! bergtsf or icldm > lcldm
!
!            !     if deposition, it should not reduce grid mean rhi below 1.0
!            if(cmei(i,k) > 0.0_r8 .and. (relhum(i,k)*esl(i,k)/esi(i,k)) > 1._r8 ) &
!                 cmei(i,k)=min(cmei(i,k),(q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat)
!
!         end if            !end ice exists loop
!         !this ends temperature < 0. loop
!
!         !-------------------------------------------------------------------
!      end if  ! 
!      !..............................................................
!
!      ! evaporation should not exceed available water
!
!      if ((-berg(i,k)).lt.-qc(i,k)/deltat) berg(i,k) = max(qc(i,k)/deltat,0._r8)
!
!      !sublimation process...
!      if (do_cldice .and. ((relhum(i,k)*esl(i,k)/esi(i,k)).lt.1._r8 .and. qiic(i,k).ge.qsmall )) then
!
!         qvi = svp_to_qsat(esi(i,k), p(i,k))
!         qvl = svp_to_qsat(esl(i,k), p(i,k))
!         dqsidt =  xxls*qvi/(rv*t(i,k)**2)
!         abi = 1._r8+dqsidt*xxls/cpp
!
!         ! get ice size distribution parameters
!
!         lami(k) = (cons1*ci* &
!              niic(i,k)/qiic(i,k))**(1._r8/di)
!         n0i(k) = niic(i,k)*lami(k)
!
!         ! check for slope
!         ! adjust vars
!         if (lami(k).lt.lammini) then
!
!            lami(k) = lammini
!            n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!         else if (lami(k).gt.lammaxi) then
!            lami(k) = lammaxi
!            n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!         end if
!
!         epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))
!
!         ! modify for ice fraction below
!         prd = epsi*(relhum(i,k)*qvl-qvi)/abi * icldm(i,k)
!         cmei(i,k)=min(prd,0._r8)
!
!      endif
!
!      ! sublimation should not exceed available ice
!      if (cmei(i,k).lt.-qi(i,k)/deltat) cmei(i,k)=-qi(i,k)/deltat
!
!      ! sublimation should not increase grid mean rhi above 1.0 
!      if(cmei(i,k) < 0.0_r8 .and. (relhum(i,k)*esl(i,k)/esi(i,k)) < 1._r8 ) &
!           cmei(i,k)=min(0._r8,max(cmei(i,k),(q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat))
!
!      ! limit cmei due for roundoff error
!
!      cmei(i,k)=cmei(i,k)*omsm
!
!      ! conditional for ice nucleation 
!      if (do_cldice .and. (t(i,k).lt.(tmelt - 5._r8))) then 
!
!         ! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
!         ! ice nucleation rate (dum2) has already been calculated and read in (naai)
!
!         dum2i(i,k)=naai(i,k)
!      else
!         dum2i(i,k)=0._r8
!      end if
!
!   end do ! i loop
!end do ! k loop
!
!
!!! initialize sub-step precip flux variables
!do i=1,ncol
!   !! flux is zero at top interface, so these should stay as 0.
!   rflx1(i,1)=0._r8
!   sflx1(i,1)=0._r8
!   do k=top_lev,pver
!
!      ! initialize normal and sub-step precip flux variables
!      rflx1(i,k+1)=0._r8
!      sflx1(i,k+1)=0._r8
!   end do ! i loop
!end do ! k loop
!!! initialize final precip flux variables.
!do i=1,ncol
!   !! flux is zero at top interface, so these should stay as 0.
!   rflx(i,1)=0._r8
!   sflx(i,1)=0._r8
!   do k=top_lev,pver
!      ! initialize normal and sub-step precip flux variables
!      rflx(i,k+1)=0._r8
!      sflx(i,k+1)=0._r8
!   end do ! i loop
!end do ! k loop
!
!do i=1,ncol
!   ltrue(i)=0
!   do k=top_lev,pver
!      ! skip microphysical calculations if no cloud water
!
!      if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall.or.cmei(i,k).ge.qsmall) ltrue(i)=1
!   end do
!end do
!
!! assign number of sub-steps to iter
!! use 2 sub-steps, following tests described in MG2008
!iter = 2
!
!! get sub-step time step
!deltat=deltat/real(iter)
!
!! since activation/nucleation processes are fast, need to take into account
!! factor mtime = mixing timescale in cloud / model time step
!! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
!! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk
!
!!        note: mtime for bulk aerosols was set to: mtime=deltat/1200._r8
!
!mtime=1._r8
!rate1ord_cw2pr_st(:,:)=0._r8 ! rce 2010/05/01
!
!!!!! skip calculations if no cloud water
!do i=1,ncol
!   if (ltrue(i).eq.0) then
!      tlat(i,1:pver)=0._r8
!      qvlat(i,1:pver)=0._r8
!      qctend(i,1:pver)=0._r8
!      qitend(i,1:pver)=0._r8
!      qnitend(i,1:pver)=0._r8
!      qrtend(i,1:pver)=0._r8
!      nctend(i,1:pver)=0._r8
!      nitend(i,1:pver)=0._r8
!      nrtend(i,1:pver)=0._r8
!      nstend(i,1:pver)=0._r8
!      prect(i)=0._r8
!      preci(i)=0._r8
!      qniic(i,1:pver)=0._r8
!      qric(i,1:pver)=0._r8
!      nsic(i,1:pver)=0._r8
!      nric(i,1:pver)=0._r8
!      rainrt(i,1:pver)=0._r8
!      goto 300
!   end if
!
!   qcsinksum_rate1ord(1:pver)=0._r8 
!   qcsum_rate1ord(1:pver)=0._r8 
!
!
!!!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !.....................................................................................................
!   do it=1,iter
!
!      ! initialize sub-step microphysical tendencies
!
!      tlat(i,1:pver)=0._r8
!      qvlat(i,1:pver)=0._r8
!      qctend(i,1:pver)=0._r8
!      qitend(i,1:pver)=0._r8
!      qnitend(i,1:pver)=0._r8
!      qrtend(i,1:pver)=0._r8
!      nctend(i,1:pver)=0._r8
!      nitend(i,1:pver)=0._r8
!      nrtend(i,1:pver)=0._r8
!      nstend(i,1:pver)=0._r8
!
!      ! initialize diagnostic precipitation to zero
!
!      qniic(i,1:pver)=0._r8
!      qric(i,1:pver)=0._r8
!      nsic(i,1:pver)=0._r8
!      nric(i,1:pver)=0._r8
!
!      rainrt(i,1:pver)=0._r8
!
!
!      ! begin new i,k loop, calculate new cldmax after adjustment to cldm above
!
!      ! initialize vertically-integrated rain and snow tendencies
!
!      qrtot = 0._r8
!      nrtot = 0._r8
!      qstot = 0._r8
!      nstot = 0._r8
!
!      ! initialize precip at surface
!
!      prect(i)=0._r8
!      preci(i)=0._r8
!
!      do k=top_lev,pver
!      
!         qcvar=relvar(i,k)
!         cons2=gamma(qcvar+2.47_r8)
!         cons3=gamma(qcvar)
!         cons9=gamma(qcvar+2._r8)
!         cons10=gamma(qcvar+1._r8)
!         cons12=gamma(qcvar+1.15_r8) 
!         cons15=gamma(qcvar+bc/3._r8)
!         cons18=qcvar**2.47_r8
!         cons19=qcvar**2
!         cons20=qcvar**1.15_r8
!
!         ! set cwml and cwmi to current qc and qi
!
!         cwml(i,k)=qc(i,k)
!         cwmi(i,k)=qi(i,k)
!
!         ! initialize precip fallspeeds to zero
!
!         ums(k)=0._r8 
!         uns(k)=0._r8 
!         umr(k)=0._r8 
!         unr(k)=0._r8
!
!         ! calculate precip fraction based on maximum overlap assumption
!
!         ! for sub-columns cldm has already been set to 1 if cloud
!         ! water or ice is present, so cldmax will be correctly set below
!         ! and nothing extra needs to be done here
!
!         if (k.eq.top_lev) then
!            cldmax(i,k)=cldm(i,k)
!         else
!            ! if rain or snow mix ratio is smaller than
!            ! threshold, then set cldmax to cloud fraction at current level
!
!            if (do_clubb_sgs) then
!               if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall) then
!                  cldmax(i,k)=cldm(i,k)
!               else
!                  cldmax(i,k)=cldmax(i,k-1)
!               end if
!            else
!
!               if (qric(i,k-1).ge.qsmall.or.qniic(i,k-1).ge.qsmall) then
!                  cldmax(i,k)=max(cldmax(i,k-1),cldm(i,k))
!               else
!                  cldmax(i,k)=cldm(i,k)
!               end if
!            endif
!         end if
!
!         ! decrease in number concentration due to sublimation/evap
!         ! divide by cloud fraction to get in-cloud decrease
!         ! don't reduce Nc due to bergeron process
!
!         if (cmei(i,k) < 0._r8 .and. qi(i,k) > qsmall .and. cldm(i,k) > mincld) then
!            nsubi(k)=cmei(i,k)/qi(i,k)*ni(i,k)/cldm(i,k)
!         else
!            nsubi(k)=0._r8
!         end if
!         nsubc(k)=0._r8
!
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
!
!         if (do_cldice .and. dum2i(i,k).gt.0._r8.and.t(i,k).lt.(tmelt - 5._r8).and. &
!              relhum(i,k)*esl(i,k)/esi(i,k).gt. rhmini+0.05_r8) then
!
!            !if NCAI > 0. then set numice = ncai (as before)
!            !note: this is gridbox averaged
!
!            nnuccd(k)=(dum2i(i,k)-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
!            nnuccd(k)=max(nnuccd(k),0._r8)
!            nimax = dum2i(i,k)*icldm(i,k)
!
!            !Calc mass of new particles using new crystal mass...
!            !also this will be multiplied by mtime as nnuccd is...
!
!            mnuccd(k) = nnuccd(k) * mi0
!
!            !  add mnuccd to cmei....
!            cmei(i,k)= cmei(i,k) + mnuccd(k) * mtime
!
!            !  limit cmei
!
!            qvi = svp_to_qsat(esi(i,k), p(i,k))
!            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
!            abi = 1._r8+dqsidt*xxls/cpp
!            cmei(i,k)=min(cmei(i,k),(q(i,k)-qvi)/abi/deltat)
!
!            ! limit for roundoff error
!            cmei(i,k)=cmei(i,k)*omsm
!
!         else
!            nnuccd(k)=0._r8
!            nimax = 0._r8
!            mnuccd(k) = 0._r8
!         end if
!
!         !c............................................................................
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
!         ! for microphysical process calculations
!         ! units are kg/kg for mixing ratio, 1/kg for number conc
!
!         ! limit in-cloud values to 0.005 kg/kg
!
!         qcic(i,k)=min(cwml(i,k)/lcldm(i,k),5.e-3_r8)
!         qiic(i,k)=min(cwmi(i,k)/icldm(i,k),5.e-3_r8)
!         ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)
!         niic(i,k)=max(ni(i,k)/icldm(i,k),0._r8)
!
!         if (qc(i,k) - berg(i,k)*deltat.lt.qsmall) then
!            qcic(i,k)=0._r8
!            ncic(i,k)=0._r8
!            if (qc(i,k)-berg(i,k)*deltat.lt.0._r8) then
!               berg(i,k)=qc(i,k)/deltat*omsm
!            end if
!         end if
!
!         if (do_cldice .and. qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.qsmall) then
!            qiic(i,k)=0._r8
!            niic(i,k)=0._r8
!            if (qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.0._r8) then
!               cmei(i,k)=(-qi(i,k)/deltat-berg(i,k))*omsm
!            end if
!         end if
!
!         ! add to cme output
!
!         cmeout(i,k) = cmeout(i,k)+cmei(i,k)
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! droplet activation
!         ! calculate potential for droplet activation if cloud water is present
!         ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
!         ! number tendency (npccnin) is read in from companion routine
!
!         ! assume aerosols already activated are equal to number of existing droplets for simplicity
!         ! multiply by cloud fraction to obtain grid-average tendency
!
!         if (qcic(i,k).ge.qsmall) then   
!            npccn(k) = max(0._r8,npccnin(i,k))  
!            dum2l(i,k)=(nc(i,k)+npccn(k)*deltat)/lcldm(i,k)
!            dum2l(i,k)=max(dum2l(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3  
!            ncmax = dum2l(i,k)*lcldm(i,k)
!         else
!            npccn(k)=0._r8
!            dum2l(i,k)=0._r8
!            ncmax = 0._r8
!         end if
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! get size distribution parameters based on in-cloud cloud water/ice 
!         ! these calculations also ensure consistency between number and mixing ratio
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         !......................................................................
!         ! cloud ice
!
!         if (qiic(i,k).ge.qsmall) then
!
!            ! add upper limit to in-cloud number concentration to prevent numerical error
!            niic(i,k)=min(niic(i,k),qiic(i,k)*1.e20_r8)
!
!            lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**(1._r8/di)
!            n0i(k) = niic(i,k)*lami(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lami(k).lt.lammini) then
!
!               lami(k) = lammini
!               n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!               niic(i,k) = n0i(k)/lami(k)
!            else if (lami(k).gt.lammaxi) then
!               lami(k) = lammaxi
!               n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!               niic(i,k) = n0i(k)/lami(k)
!            end if
!
!         else
!            lami(k) = 0._r8
!            n0i(k) = 0._r8
!         end if
!
!         if (qcic(i,k).ge.qsmall) then
!
!            ! add upper limit to in-cloud number concentration to prevent numerical error
!            ncic(i,k)=min(ncic(i,k),qcic(i,k)*1.e20_r8)
!
!            ncic(i,k)=max(ncic(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm  
!
!            ! get pgam from fit to observations of martin et al. 1994
!
!            pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
!            pgam(k)=1._r8/(pgam(k)**2)-1._r8
!            pgam(k)=max(pgam(k),2._r8)
!            pgam(k)=min(pgam(k),15._r8)
!
!            ! calculate lamc
!
!            lamc(k) = (pi/6._r8*rhow*ncic(i,k)*gamma(pgam(k)+4._r8)/ &
!                 (qcic(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!
!            ! lammin, 50 micron diameter max mean size
!
!            lammin = (pgam(k)+1._r8)/50.e-6_r8
!            lammax = (pgam(k)+1._r8)/2.e-6_r8
!
!            if (lamc(k).lt.lammin) then
!               lamc(k) = lammin
!               ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
!                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
!            else if (lamc(k).gt.lammax) then
!               lamc(k) = lammax
!               ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
!                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
!            end if
!
!            ! parameter to calculate droplet freezing
!
!            cdist1(k) = ncic(i,k)/gamma(pgam(k)+1._r8) 
!
!         else
!            lamc(k) = 0._r8
!            cdist1(k) = 0._r8
!         end if
!
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! begin micropysical process calculations 
!         !.................................................................
!         ! autoconversion of cloud liquid water to rain
!         ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
!         ! minimum qc of 1 x 10^-8 prevents floating point error
!
!         if (qcic(i,k).ge.1.e-8_r8) then
!
!            ! nprc is increase in rain number conc due to autoconversion
!            ! nprc1 is decrease in cloud droplet conc due to autoconversion
!
!            ! assume exponential sub-grid distribution of qc, resulting in additional
!            ! factor related to qcvar below
!
!            ! hm switch for sub-columns, don't include sub-grid qc
!            if (microp_uniform) then
!
!               prc(k) = 1350._r8*qcic(i,k)**2.47_r8* &
!                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
!               nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
!               nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
!
!            else
!
!               prc(k) = cons2/(cons3*cons18)*1350._r8*qcic(i,k)**2.47_r8* &
!                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
!               nprc(k) = prc(k)/cons22
!               nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
!
!            end if               ! sub-column switch
!
!         else
!            prc(k)=0._r8
!            nprc(k)=0._r8
!            nprc1(k)=0._r8
!         end if
!
!         ! add autoconversion to precip from above to get provisional rain mixing ratio
!         ! and number concentration (qric and nric)
!
!         ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)
!
!         dum=0.45_r8
!         dum1=0.45_r8
!
!         if (k.eq.top_lev) then
!            qric(i,k)=prc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
!            nric(i,k)=nprc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
!         else
!            if (qric(i,k-1).ge.qsmall) then
!               dum=umr(k-1)
!               dum1=unr(k-1)
!            end if
!
!            ! no autoconversion of rain number if rain/snow falling from above
!            ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
!            ! by the existing rain/snow particles from above
!
!            if (qric(i,k-1).ge.1.e-9_r8.or.qniic(i,k-1).ge.1.e-9_r8) then
!               nprc(k)=0._r8
!            end if
!
!            qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
!                 (rho(i,k)*dz(i,k)*((pra(k-1)+prc(k))*lcldm(i,k)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax(i,k))))&
!                 /(dum*rho(i,k)*cldmax(i,k))
!            nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
!                 (rho(i,k)*dz(i,k)*(nprc(k)*lcldm(i,k)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax(i,k))))&
!                 /(dum1*rho(i,k)*cldmax(i,k))
!
!         end if
!
!         !.......................................................................
!         ! Autoconversion of cloud ice to snow
!         ! similar to Ferrier (1994)
!
!         if (do_cldice) then
!            if (t(i,k).le.273.15_r8.and.qiic(i,k).ge.qsmall) then
!
!               ! note: assumes autoconversion timescale of 180 sec
!               
!               nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)
!
!               prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
!                    (cons23/lami(k)+3._r8*cons24/lami(k)**2+ &
!                    6._r8*dcs/lami(k)**3+6._r8/lami(k)**4)*exp(-lami(k)*dcs)
!            else
!               prci(k)=0._r8
!               nprci(k)=0._r8
!            end if
!         else
!            ! Add in the particles that we have already converted to snow, and
!            ! don't do any further autoconversion of ice.
!            prci(k)  = tnd_qsnow(i, k) / cldm(i,k)
!            nprci(k) = tnd_nsnow(i, k) / cldm(i,k)
!         end if
!
!         ! add autoconversion to flux from level above to get provisional snow mixing ratio
!         ! and number concentration (qniic and nsic)
!
!         dum=(asn(i,k)*cons25)
!         dum1=(asn(i,k)*cons25)
!
!         if (k.eq.top_lev) then
!            qniic(i,k)=prci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
!            nsic(i,k)=nprci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
!         else
!            if (qniic(i,k-1).ge.qsmall) then
!               dum=ums(k-1)
!               dum1=uns(k-1)
!            end if
!
!            qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
!                 (rho(i,k)*dz(i,k)*((prci(k)+prai(k-1)+psacws(k-1)+bergs(k-1))*icldm(i,k)+(prds(k-1)+ &
!                 pracs(k-1)+mnuccr(k-1))*cldmax(i,k))))&
!                 /(dum*rho(i,k)*cldmax(i,k))
!
!            nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
!                 (rho(i,k)*dz(i,k)*(nprci(k)*icldm(i,k)+(nsubs(k-1)+nsagg(k-1)+nnuccr(k-1))*cldmax(i,k))))&
!                 /(dum1*rho(i,k)*cldmax(i,k))
!
!         end if
!
!         ! if precip mix ratio is zero so should number concentration
!
!         if (qniic(i,k).lt.qsmall) then
!            qniic(i,k)=0._r8
!            nsic(i,k)=0._r8
!         end if
!
!         if (qric(i,k).lt.qsmall) then
!            qric(i,k)=0._r8
!            nric(i,k)=0._r8
!         end if
!
!         ! make sure number concentration is a positive number to avoid 
!         ! taking root of negative later
!
!         nric(i,k)=max(nric(i,k),0._r8)
!         nsic(i,k)=max(nsic(i,k),0._r8)
!
!         !.......................................................................
!         ! get size distribution parameters for precip
!         !......................................................................
!         ! rain
!
!         if (qric(i,k).ge.qsmall) then
!            lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
!            n0r(k) = nric(i,k)*lamr(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lamr(k).lt.lamminr) then
!
!               lamr(k) = lamminr
!
!               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
!               nric(i,k) = n0r(k)/lamr(k)
!            else if (lamr(k).gt.lammaxr) then
!               lamr(k) = lammaxr
!               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
!               nric(i,k) = n0r(k)/lamr(k)
!            end if
!
!            ! provisional rain number and mass weighted mean fallspeed (m/s)
!
!            unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
!            umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(i,k))
!
!         else
!            lamr(k) = 0._r8
!            n0r(k) = 0._r8
!            umr(k) = 0._r8
!            unr(k) = 0._r8
!         end if
!
!         !......................................................................
!         ! snow
!
!         if (qniic(i,k).ge.qsmall) then
!            lams(k) = (cons6*cs*nsic(i,k)/qniic(i,k))**(1._r8/ds)
!            n0s(k) = nsic(i,k)*lams(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lams(k).lt.lammins) then
!               lams(k) = lammins
!               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
!               nsic(i,k) = n0s(k)/lams(k)
!
!            else if (lams(k).gt.lammaxs) then
!               lams(k) = lammaxs
!               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
!               nsic(i,k) = n0s(k)/lams(k)
!            end if
!
!            ! provisional snow number and mass weighted mean fallspeed (m/s)
!
!            ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(i,k))
!            uns(k) = min(asn(i,k)*cons7/lams(k)**bs,1.2_r8*rhof(i,k))
!
!         else
!            lams(k) = 0._r8
!            n0s(k) = 0._r8
!            ums(k) = 0._r8
!            uns(k) = 0._r8
!         end if
!
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         ! heterogeneous freezing of cloud water
!
!         if (.not. use_hetfrz_classnuc) then
!
!            if (do_cldice .and. qcic(i,k).ge.qsmall .and. t(i,k).lt.269.15_r8) then
!
!               ! immersion freezing (Bigg, 1953)
!
!
!               ! subcolumns
!
!               if (microp_uniform) then
!
!                  mnuccc(k) = &
!                     pi*pi/36._r8*rhow* &
!                     cdist1(k)*gamma(7._r8+pgam(k))* &
!                     bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
!                     lamc(k)**3/lamc(k)**3
!
!                  nnuccc(k) = &
!                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
!                     *bimm* &
!                     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3
!
!               else
!
!                  mnuccc(k) = cons9/(cons3*cons19)* &
!                     pi*pi/36._r8*rhow* &
!                     cdist1(k)*gamma(7._r8+pgam(k))* &
!                     bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
!                     lamc(k)**3/lamc(k)**3
!
!                  nnuccc(k) = cons10/(cons3*qcvar)* &
!                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
!                     *bimm* &
!                     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3
!               end if           ! sub-columns
!
!
!               ! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
!               ! dust size and number in 4 bins are read in from companion routine
!
!               tcnt=(270.16_r8-t(i,k))**1.3_r8
!               viscosity=1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
!               mfp=2.0_r8*viscosity/(p(i,k)  &                   ! Mean free path (m)
!                  *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k))))           
!
!               nslip1=1.0_r8+(mfp/rndst(i,k,1))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,1)/mfp))))! Slip correction factor
!               nslip2=1.0_r8+(mfp/rndst(i,k,2))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,2)/mfp))))
!               nslip3=1.0_r8+(mfp/rndst(i,k,3))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,3)/mfp))))
!               nslip4=1.0_r8+(mfp/rndst(i,k,4))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,4)/mfp))))
!
!               ndfaer1=1.381e-23_r8*t(i,k)*nslip1/(6._r8*pi*viscosity*rndst(i,k,1))  ! aerosol diffusivity (m2/s)
!               ndfaer2=1.381e-23_r8*t(i,k)*nslip2/(6._r8*pi*viscosity*rndst(i,k,2))
!               ndfaer3=1.381e-23_r8*t(i,k)*nslip3/(6._r8*pi*viscosity*rndst(i,k,3))
!               ndfaer4=1.381e-23_r8*t(i,k)*nslip4/(6._r8*pi*viscosity*rndst(i,k,4))
!
!
!               if (microp_uniform) then
!
!                  mnucct(k) = &
!                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
!                     cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
!
!                  nnucct(k) = (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
!                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
!
!               else
!
!                  mnucct(k) = gamma(qcvar+4._r8/3._r8)/(cons3*qcvar**(4._r8/3._r8))*  &
!                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
!                     cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
!
!                  nnucct(k) =  gamma(qcvar+1._r8/3._r8)/(cons3*qcvar**(1._r8/3._r8))*  &
!                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
!                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
!
!               end if      ! sub-column switch
!
!               ! make sure number of droplets frozen does not exceed available ice nuclei concentration
!               ! this prevents 'runaway' droplet freezing
!
!               if (nnuccc(k)*lcldm(i,k).gt.nnuccd(k)) then
!                  dum=(nnuccd(k)/(nnuccc(k)*lcldm(i,k)))
!                  ! scale mixing ratio of droplet freezing with limit
!                  mnuccc(k)=mnuccc(k)*dum
!                  nnuccc(k)=nnuccd(k)/lcldm(i,k)
!               end if
!
!            else
!               mnuccc(k)=0._r8
!               nnuccc(k)=0._r8
!               mnucct(k)=0._r8
!               nnucct(k)=0._r8
!            end if
!
!         else
!            if (do_cldice .and. qcic(i,k) >= qsmall) then
!               con1 = 1._r8/(1.333_r8*pi)**0.333_r8
!               r3lx = con1*(rho(i,k)*qcic(i,k)/(rhow*max(ncic(i,k)*rho(i,k), 1.0e6_r8)))**0.333_r8 ! in m
!               r3lx = max(4.e-6_r8, r3lx)
!               mi0l = 4._r8/3._r8*pi*rhow*r3lx**3_r8
!                
!               nnuccc(k) = frzimm(i,k)*1.0e6_r8/rho(i,k)
!               mnuccc(k) = nnuccc(k)*mi0l 
!
!               nnucct(k) = frzcnt(i,k)*1.0e6_r8/rho(i,k)
!               mnucct(k) = nnucct(k)*mi0l 
!
!               nnudep(k) = frzdep(i,k)*1.0e6_r8/rho(i,k)
!               mnudep(k) = nnudep(k)*mi0
!            else
!               nnuccc(k) = 0._r8
!               mnuccc(k) = 0._r8
!
!               nnucct(k) = 0._r8
!               mnucct(k) = 0._r8
!
!               nnudep(k) = 0._r8
!               mnudep(k) = 0._r8
!            end if
!         endif
!
!
!         !.......................................................................
!         ! snow self-aggregation from passarelli, 1978, used by reisner, 1998
!         ! this is hard-wired for bs = 0.4 for now
!         ! ignore self-collection of cloud ice
!
!         if (qniic(i,k).ge.qsmall .and. t(i,k).le.273.15_r8) then
!            nsagg(k) = -1108._r8*asn(i,k)*Eii* &
!                 pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)*rho(i,k)** &
!                 ((2._r8+bs)/3._r8)*qniic(i,k)**((2._r8+bs)/3._r8)* &
!                 (nsic(i,k)*rho(i,k))**((4._r8-bs)/3._r8)/ &
!                 (4._r8*720._r8*rho(i,k))
!         else
!            nsagg(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! accretion of cloud droplets onto snow/graupel
!         ! here use continuous collection equation with
!         ! simple gravitational collection kernel
!         ! ignore collisions between droplets/cloud ice
!         ! since minimum size ice particle for accretion is 50 - 150 micron
!
!         ! ignore collision of snow with droplets above freezing
!
!         if (qniic(i,k).ge.qsmall .and. t(i,k).le.tmelt .and. &
!              qcic(i,k).ge.qsmall) then
!
!            ! put in size dependent collection efficiency
!            ! mean diameter of snow is area-weighted, since
!            ! accretion is function of crystal geometric area
!            ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)
!
!            dc0 = (pgam(k)+1._r8)/lamc(k)
!            ds0 = 1._r8/lams(k)
!            dum = dc0*dc0*uns(k)*rhow/(9._r8*mu(i,k)*ds0)
!            eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))
!
!            eci = max(eci,0._r8)
!            eci = min(eci,1._r8)
!
!
!            ! no impact of sub-grid distribution of qc since psacws
!            ! is linear in qc
!
!            psacws(k) = pi/4._r8*asn(i,k)*qcic(i,k)*rho(i,k)* &
!                 n0s(k)*Eci*cons11/ &
!                 lams(k)**(bs+3._r8)
!            npsacws(k) = pi/4._r8*asn(i,k)*ncic(i,k)*rho(i,k)* &
!                 n0s(k)*Eci*cons11/ &
!                 lams(k)**(bs+3._r8)
!         else
!            psacws(k)=0._r8
!            npsacws(k)=0._r8
!         end if
!
!         ! add secondary ice production due to accretion of droplets by snow 
!         ! (Hallet-Mossop process) (from Cotton et al., 1986)
!
!         if (.not. do_cldice) then
!            ni_secp   = 0.0_r8
!            nsacwi(k) = 0.0_r8
!            msacwi(k) = 0.0_r8
!         else if((t(i,k).lt.270.16_r8) .and. (t(i,k).ge.268.16_r8)) then
!            ni_secp   = 3.5e8_r8*(270.16_r8-t(i,k))/2.0_r8*psacws(k)
!            nsacwi(k) = ni_secp
!            msacwi(k) = min(ni_secp*mi0,psacws(k))
!         else if((t(i,k).lt.268.16_r8) .and. (t(i,k).ge.265.16_r8)) then
!            ni_secp   = 3.5e8_r8*(t(i,k)-265.16_r8)/3.0_r8*psacws(k)
!            nsacwi(k) = ni_secp
!            msacwi(k) = min(ni_secp*mi0,psacws(k))
!         else
!            ni_secp   = 0.0_r8
!            nsacwi(k) = 0.0_r8
!            msacwi(k) = 0.0_r8
!         endif
!         psacws(k) = max(0.0_r8,psacws(k)-ni_secp*mi0)
!
!         !.......................................................................
!         ! accretion of rain water by snow
!         ! formula from ikawa and saito, 1991, used by reisner et al., 1998
!
!         if (qric(i,k).ge.1.e-8_r8 .and. qniic(i,k).ge.1.e-8_r8 .and. & 
!              t(i,k).le.273.15_r8) then
!
!            pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))**2+ &
!                 0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(i,k)* &
!                 n0r(k)*n0s(k)* &
!                 (5._r8/(lamr(k)**6*lams(k))+ &
!                 2._r8/(lamr(k)**5*lams(k)**2)+ &
!                 0.5_r8/(lamr(k)**4*lams(k)**3)))
!
!            npracs(k) = pi/2._r8*rho(i,k)*ecr*(1.7_r8*(unr(k)-uns(k))**2+ &
!                 0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
!                 (1._r8/(lamr(k)**3*lams(k))+ &
!                 1._r8/(lamr(k)**2*lams(k)**2)+ &
!                 1._r8/(lamr(k)*lams(k)**3))
!
!         else
!            pracs(k)=0._r8
!            npracs(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! heterogeneous freezing of rain drops
!         ! follows from Bigg (1953)
!
!         if (t(i,k).lt.269.15_r8 .and. qric(i,k).ge.qsmall) then
!
!            mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
!                 (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamr(k)**3 &
!                 /lamr(k)**3
!
!            nnuccr(k) = pi*nric(i,k)*bimm* &
!                 (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamr(k)**3
!         else
!            mnuccr(k)=0._r8
!            nnuccr(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! accretion of cloud liquid water by rain
!         ! formula from Khrouditnov and Kogan (2000)
!         ! gravitational collection kernel, droplet fall speed neglected
!
!         if (qric(i,k).ge.qsmall .and. qcic(i,k).ge.qsmall) then
!
!            ! include sub-grid distribution of cloud water
!
!            ! add sub-column switch
!
!            if (microp_uniform) then
!
!               pra(k) = 67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
!               npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
!
!            else
!
!               pra(k) = accre_enhan(i,k)*(cons12/(cons3*cons20)*67._r8*(qcic(i,k)*qric(i,k))**1.15_r8)
!               npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
!
!            end if               ! sub-column switch
!
!         else
!            pra(k)=0._r8
!            npra(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! Self-collection of rain drops
!         ! from Beheng(1994)
!
!         if (qric(i,k).ge.qsmall) then
!            nragg(k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
!         else
!            nragg(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! Accretion of cloud ice by snow
!         ! For this calculation, it is assumed that the Vs >> Vi
!         ! and Ds >> Di for continuous collection
!
!         if (do_cldice .and. qniic(i,k).ge.qsmall.and.qiic(i,k).ge.qsmall &
!              .and.t(i,k).le.273.15_r8) then
!
!            prai(k) = pi/4._r8*asn(i,k)*qiic(i,k)*rho(i,k)* &
!                 n0s(k)*Eii*cons11/ &
!                 lams(k)**(bs+3._r8)
!            nprai(k) = pi/4._r8*asn(i,k)*niic(i,k)* &
!                 rho(i,k)*n0s(k)*Eii*cons11/ &
!                 lams(k)**(bs+3._r8)
!         else
!            prai(k)=0._r8
!            nprai(k)=0._r8
!         end if
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! calculate evaporation/sublimation of rain and snow
!         ! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
!         ! in-cloud condensation/deposition of rain and snow is neglected
!         ! except for transfer of cloud water to snow through bergeron process
!
!         ! initialize evap/sub tendncies
!         pre(k)=0._r8
!         prds(k)=0._r8
!
!         ! evaporation of rain
!         ! only calculate if there is some precip fraction > cloud fraction
!
!         if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8.or.cldmax(i,k).gt.lcldm(i,k)) then
!
!            ! set temporary cloud fraction to zero if cloud water + ice is very small
!            ! this will ensure that evaporation/sublimation of precip occurs over
!            ! entire grid cell, since min cloud fraction is specified otherwise
!            if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8) then
!               dum=0._r8
!            else
!               dum=lcldm(i,k)
!            end if
!
!            ! saturation vapor pressure
!            esn=svp_water(t(i,k))
!            qsn=svp_to_qsat(esn, p(i,k))
!
!            ! recalculate saturation vapor pressure for liquid and ice
!            esl(i,k)=esn
!            esi(i,k)=svp_ice(t(i,k))
!            ! hm fix, make sure when above freezing that esi=esl, not active yet
!            if (t(i,k).gt.tmelt)esi(i,k)=esl(i,k)
!
!            ! calculate q for out-of-cloud region
!            qclr=(q(i,k)-dum*qsn)/(1._r8-dum)
!
!            if (qric(i,k).ge.qsmall) then
!
!               qvs=svp_to_qsat(esl(i,k), p(i,k))
!               dqsdt = xxlv*qvs/(rv*t(i,k)**2)
!               ab = 1._r8+dqsdt*xxlv/cpp
!               epsr = 2._r8*pi*n0r(k)*rho(i,k)*Dv(i,k)* &
!                    (f1r/(lamr(k)*lamr(k))+ &
!                    f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
!                    sc(i,k)**(1._r8/3._r8)*cons13/ &
!                    (lamr(k)**(5._r8/2._r8+br/2._r8)))
!
!               pre(k) = epsr*(qclr-qvs)/ab
!
!               ! only evaporate in out-of-cloud region
!               ! and distribute across cldmax
!               pre(k)=min(pre(k)*(cldmax(i,k)-dum),0._r8)
!               pre(k)=pre(k)/cldmax(i,k)
!               am_evp_st(i,k) = max(cldmax(i,k)-dum, 0._r8)
!            end if
!
!            ! sublimation of snow
!            if (qniic(i,k).ge.qsmall) then
!               qvi=svp_to_qsat(esi(i,k), p(i,k))
!               dqsidt =  xxls*qvi/(rv*t(i,k)**2)
!               abi = 1._r8+dqsidt*xxls/cpp
!               epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
!                    (f1s/(lams(k)*lams(k))+ &
!                    f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
!                    sc(i,k)**(1._r8/3._r8)*cons14/ &
!                    (lams(k)**(5._r8/2._r8+bs/2._r8)))
!               prds(k) = epss*(qclr-qvi)/abi
!
!               ! only sublimate in out-of-cloud region and distribute over cldmax
!               prds(k)=min(prds(k)*(cldmax(i,k)-dum),0._r8)
!               prds(k)=prds(k)/cldmax(i,k)
!               am_evp_st(i,k) = max(cldmax(i,k)-dum, 0._r8)
!            end if
!
!            ! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
!            ! get updated RH at end of time step based on cloud water/ice condensation/evap
!
!            qtmp=q(i,k)-(cmei(i,k)+(pre(k)+prds(k))*cldmax(i,k))*deltat
!            ttmp=t(i,k)+((pre(k)*cldmax(i,k))*xxlv+ &
!                 (cmei(i,k)+prds(k)*cldmax(i,k))*xxls)*deltat/cpp
!
!            !limit range of temperatures!
!            ttmp=max(180._r8,min(ttmp,323._r8))
!
!            esn=svp_water(ttmp)  ! use rhw to allow ice supersaturation
!            qsn=svp_to_qsat(esn, p(i,k))
!
!            ! modify precip evaporation rate if q > qsat
!            if (qtmp.gt.qsn) then
!               if (pre(k)+prds(k).lt.-1.e-20_r8) then
!                  dum1=pre(k)/(pre(k)+prds(k))
!                  ! recalculate q and t after cloud water cond but without precip evap
!                  qtmp=q(i,k)-(cmei(i,k))*deltat
!                  ttmp=t(i,k)+(cmei(i,k)*xxls)*deltat/cpp
!                  esn=svp_water(ttmp) ! use rhw to allow ice supersaturation
!                  qsn=svp_to_qsat(esn, p(i,k))
!                  dum=(qtmp-qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))
!                  dum=min(dum,0._r8)
!
!                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
!                  pre(k)=dum*dum1/deltat/cldmax(i,k)
!
!                  ! do separately using RHI for prds....
!                  esn=svp_ice(ttmp) ! use rhi to allow ice supersaturation
!                  qsn=svp_to_qsat(esn, p(i,k))
!                  dum=(qtmp-qsn)/(1._r8 + cons28*qsn/(cpp*rv*ttmp**2))
!                  dum=min(dum,0._r8)
!
!                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
!                  prds(k)=dum*(1._r8-dum1)/deltat/cldmax(i,k)
!               end if
!            end if
!         end if
!
!         ! bergeron process - evaporation of droplets and deposition onto snow
!
!         if (qniic(i,k).ge.qsmall.and.qcic(i,k).ge.qsmall.and.t(i,k).lt.tmelt) then
!            qvi=svp_to_qsat(esi(i,k), p(i,k))
!            qvs=svp_to_qsat(esl(i,k), p(i,k))
!            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
!            abi = 1._r8+dqsidt*xxls/cpp
!            epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
!                 (f1s/(lams(k)*lams(k))+ &
!                 f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
!                 sc(i,k)**(1._r8/3._r8)*cons14/ &
!                 (lams(k)**(5._r8/2._r8+bs/2._r8)))
!            bergs(k)=epss*(qvs-qvi)/abi
!         else
!            bergs(k)=0._r8
!         end if
!
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! conservation to ensure no negative values of cloud water/precipitation
!         ! in case microphysical process rates are large
!
!         ! make sure and use end-of-time step values for cloud water, ice, due
!         ! condensation/deposition
!
!         ! note: for check on conservation, processes are multiplied by omsm
!         ! to prevent problems due to round off error
!
!         ! include mixing timescale  (mtime)
!
!         qce=(qc(i,k) - berg(i,k)*deltat)
!         nce=(nc(i,k)+npccn(k)*deltat*mtime)
!         qie=(qi(i,k)+(cmei(i,k)+berg(i,k))*deltat)
!         nie=(ni(i,k)+nnuccd(k)*deltat*mtime)
!
!         ! conservation of qc
!
!         dum = (prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+ &
!              psacws(k)+bergs(k))*lcldm(i,k)*deltat
!
!         if (dum.gt.qce) then
!            ratio = qce/deltat/lcldm(i,k)/(prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+psacws(k)+bergs(k))*omsm 
!
!            prc(k) = prc(k)*ratio
!            pra(k) = pra(k)*ratio
!            mnuccc(k) = mnuccc(k)*ratio
!            mnucct(k) = mnucct(k)*ratio  
!            msacwi(k) = msacwi(k)*ratio  
!            psacws(k) = psacws(k)*ratio
!            bergs(k) = bergs(k)*ratio
!         end if
!
!         ! conservation of nc
!
!         dum = (nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
!              npsacws(k)-nsubc(k))*lcldm(i,k)*deltat
!
!         if (dum.gt.nce) then
!            ratio = nce/deltat/((nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+&
!                 npsacws(k)-nsubc(k))*lcldm(i,k))*omsm
!
!            nprc1(k) = nprc1(k)*ratio
!            npra(k) = npra(k)*ratio
!            nnuccc(k) = nnuccc(k)*ratio
!            nnucct(k) = nnucct(k)*ratio  
!            npsacws(k) = npsacws(k)*ratio
!            nsubc(k)=nsubc(k)*ratio
!         end if
!
!         ! conservation of qi
!
!         if (do_cldice) then
!
!            frztmp = -mnuccc(k) - mnucct(k) - msacwi(k)
!            if (use_hetfrz_classnuc) frztmp = -mnuccc(k)-mnucct(k)-mnudep(k)-msacwi(k)
!            dum = ( frztmp*lcldm(i,k) + (prci(k)+prai(k))*icldm(i,k) )*deltat
!
!            if (dum.gt.qie) then
!
!               frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
!               if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
!               ratio = (qie/deltat + frztmp*lcldm(i,k))/((prci(k)+prai(k))*icldm(i,k))*omsm 
!               prci(k) = prci(k)*ratio
!               prai(k) = prai(k)*ratio
!            end if
!
!            ! conservation of ni
!            frztmp = -nnucct(k) - nsacwi(k)
!            if (use_hetfrz_classnuc) frztmp = -nnucct(k) - nnuccc(k) - nnudep(k) - nsacwi(k)
!            dum = ( frztmp*lcldm(i,k) + (nprci(k)+nprai(k)-nsubi(k))*icldm(i,k) )*deltat
!
!            if (dum.gt.nie) then
!
!               frztmp = nnucct(k) + nsacwi(k)
!               if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
!               ratio = (nie/deltat + frztmp*lcldm(i,k))/ &  
!                     ((nprci(k)+nprai(k)-nsubi(k))*icldm(i,k))*omsm
!               nprci(k) = nprci(k)*ratio
!               nprai(k) = nprai(k)*ratio
!               nsubi(k) = nsubi(k)*ratio
!            end if
!         end if
!
!         ! for precipitation conservation, use logic that vertical integral 
!         ! of tendency from current level to top of model (i.e., qrtot) cannot be negative
!
!         ! conservation of rain mixing rat
!
!         if (((prc(k)+pra(k))*lcldm(i,k)+(-mnuccr(k)+pre(k)-pracs(k))*&
!              cldmax(i,k))*dz(i,k)*rho(i,k)+qrtot.lt.0._r8) then
!
!            if (-pre(k)+pracs(k)+mnuccr(k).ge.qsmall) then
!
!               ratio = (qrtot/(dz(i,k)*rho(i,k))+(prc(k)+pra(k))*lcldm(i,k))/&
!                    ((-pre(k)+pracs(k)+mnuccr(k))*cldmax(i,k))*omsm 
!
!               pre(k) = pre(k)*ratio
!               pracs(k) = pracs(k)*ratio
!               mnuccr(k) = mnuccr(k)*ratio
!            end if
!         end if
!
!         ! conservation of nr
!         ! for now neglect evaporation of nr
!         nsubr(k)=0._r8
!
!         if ((nprc(k)*lcldm(i,k)+(-nnuccr(k)+nsubr(k)-npracs(k)&
!              +nragg(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+nrtot.lt.0._r8) then
!
!            if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).ge.qsmall) then
!
!               ratio = (nrtot/(dz(i,k)*rho(i,k))+nprc(k)*lcldm(i,k))/&
!                    ((-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k))*cldmax(i,k))*omsm
!
!               nsubr(k) = nsubr(k)*ratio
!               npracs(k) = npracs(k)*ratio
!               nnuccr(k) = nnuccr(k)*ratio
!               nragg(k) = nragg(k)*ratio
!            end if
!         end if
!
!         ! conservation of snow mix ratio
!
!         if (((bergs(k)+psacws(k))*lcldm(i,k)+(prai(k)+prci(k))*icldm(i,k)+(pracs(k)+&
!              mnuccr(k)+prds(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+qstot.lt.0._r8) then
!
!            if (-prds(k).ge.qsmall) then
!
!               ratio = (qstot/(dz(i,k)*rho(i,k))+(bergs(k)+psacws(k))*lcldm(i,k)+(prai(k)+prci(k))*icldm(i,k)+&
!                    (pracs(k)+mnuccr(k))*cldmax(i,k))/(-prds(k)*cldmax(i,k))*omsm
!
!               prds(k) = prds(k)*ratio
!            end if
!         end if
!
!         ! conservation of ns
!
!         ! calculate loss of number due to sublimation
!         ! for now neglect sublimation of ns
!         nsubs(k)=0._r8
!
!         if ((nprci(k)*icldm(i,k)+(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax(i,k))*&
!              dz(i,k)*rho(i,k)+nstot.lt.0._r8) then
!
!            if (-nsubs(k)-nsagg(k).ge.qsmall) then
!
!               ratio = (nstot/(dz(i,k)*rho(i,k))+nprci(k)*icldm(i,k)+&
!                    nnuccr(k)*cldmax(i,k))/((-nsubs(k)-nsagg(k))*cldmax(i,k))*omsm
!
!               nsubs(k) = nsubs(k)*ratio
!               nsagg(k) = nsagg(k)*ratio
!            end if
!         end if
!
!         ! get tendencies due to microphysical conversion processes
!         ! note: tendencies are multiplied by appropaiate cloud/precip 
!         ! fraction to get grid-scale values
!         ! note: cmei is already grid-average values
!
!         qvlat(i,k) = qvlat(i,k)-(pre(k)+prds(k))*cldmax(i,k)-cmei(i,k) 
!
!         tlat(i,k) = tlat(i,k)+((pre(k)*cldmax(i,k)) &
!              *xxlv+(prds(k)*cldmax(i,k)+cmei(i,k))*xxls+ &
!              ((bergs(k)+psacws(k)+mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(i,k)+(mnuccr(k)+ &
!              pracs(k))*cldmax(i,k)+berg(i,k))*xlf)
!
!         qctend(i,k) = qctend(i,k)+ &
!              (-pra(k)-prc(k)-mnuccc(k)-mnucct(k)-msacwi(k)- & 
!              psacws(k)-bergs(k))*lcldm(i,k)-berg(i,k)
!
!         if (do_cldice) then
!
!            frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
!            if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
!            qitend(i,k) = qitend(i,k) + frztmp*lcldm(i,k) + &
!               (-prci(k)-prai(k))*icldm(i,k) + cmei(i,k) + berg(i,k)
!
!         end if
!
!         qrtend(i,k) = qrtend(i,k)+ &
!              (pra(k)+prc(k))*lcldm(i,k)+(pre(k)-pracs(k)- &
!              mnuccr(k))*cldmax(i,k)
!
!         qnitend(i,k) = qnitend(i,k)+ &
!              (prai(k)+prci(k))*icldm(i,k)+(psacws(k)+bergs(k))*lcldm(i,k)+(prds(k)+ &
!              pracs(k)+mnuccr(k))*cldmax(i,k)
!
!         ! add output for cmei (accumulate)
!         cmeiout(i,k) = cmeiout(i,k) + cmei(i,k)
!
!         ! assign variables for trop_mozart, these are grid-average
!         ! evaporation/sublimation is stored here as positive term
!
!         evapsnow(i,k) = evapsnow(i,k)-prds(k)*cldmax(i,k)
!         nevapr(i,k) = nevapr(i,k)-pre(k)*cldmax(i,k)
!         nevapr2(i,k) = nevapr2(i,k)-pre(k)*cldmax(i,k)
!
!         ! change to make sure prain is positive: do not remove snow from
!         ! prain used for wet deposition
!         prain(i,k) = prain(i,k)+(pra(k)+prc(k))*lcldm(i,k)+(-pracs(k)- &
!              mnuccr(k))*cldmax(i,k)
!         prodsnow(i,k) = prodsnow(i,k)+(prai(k)+prci(k))*icldm(i,k)+(psacws(k)+bergs(k))*lcldm(i,k)+(&
!              pracs(k)+mnuccr(k))*cldmax(i,k)
!
!         ! following are used to calculate 1st order conversion rate of cloud water
!         !    to rain and snow (1/s), for later use in aerosol wet removal routine
!         ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
!         !    used to calculate pra, prc, ... in this routine
!         ! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
!         !                      (no cloud ice or bergeron terms)
!         ! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }
!
!         qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k) + (pra(k)+prc(k)+psacws(k))*lcldm(i,k) 
!         qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc(i,k) 
!
!         ! microphysics output, note this is grid-averaged
!         prao(i,k)=prao(i,k)+pra(k)*lcldm(i,k)
!         prco(i,k)=prco(i,k)+prc(k)*lcldm(i,k)
!         mnuccco(i,k)=mnuccco(i,k)+mnuccc(k)*lcldm(i,k)
!         mnuccto(i,k)=mnuccto(i,k)+mnucct(k)*lcldm(i,k)
!         mnuccdo(i,k)=mnuccdo(i,k)+mnuccd(k)*lcldm(i,k)
!         msacwio(i,k)=msacwio(i,k)+msacwi(k)*lcldm(i,k)
!         psacwso(i,k)=psacwso(i,k)+psacws(k)*lcldm(i,k)
!         bergso(i,k)=bergso(i,k)+bergs(k)*lcldm(i,k)
!         bergo(i,k)=bergo(i,k)+berg(i,k)
!         prcio(i,k)=prcio(i,k)+prci(k)*icldm(i,k)
!         praio(i,k)=praio(i,k)+prai(k)*icldm(i,k)
!         mnuccro(i,k)=mnuccro(i,k)+mnuccr(k)*cldmax(i,k)
!         pracso (i,k)=pracso (i,k)+pracs (k)*cldmax(i,k)
!
!         ! multiply activation/nucleation by mtime to account for fast timescale
!
!         nctend(i,k) = nctend(i,k)+ npccn(k)*mtime+&
!              (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k) & 
!              -npra(k)-nprc1(k))*lcldm(i,k)      
!
!         if (do_cldice) then
!
!            frztmp = nnucct(k) + nsacwi(k)
!            if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
!            nitend(i,k) = nitend(i,k) + nnuccd(k)*mtime + & 
!                  frztmp*lcldm(i,k) + (nsubi(k)-nprci(k)-nprai(k))*icldm(i,k)
!
!         end if
!
!         nstend(i,k) = nstend(i,k)+(nsubs(k)+ &
!              nsagg(k)+nnuccr(k))*cldmax(i,k)+nprci(k)*icldm(i,k)
!
!         nrtend(i,k) = nrtend(i,k)+ &
!              nprc(k)*lcldm(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
!              +nragg(k))*cldmax(i,k)
!
!         ! make sure that nc and ni at advanced time step do not exceed
!         ! maximum (existing N + source terms*dt), which is possible due to
!         ! fast nucleation timescale
!
!         if (nctend(i,k).gt.0._r8.and.nc(i,k)+nctend(i,k)*deltat.gt.ncmax) then
!            nctend(i,k)=max(0._r8,(ncmax-nc(i,k))/deltat)
!         end if
!
!         if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax) then
!            nitend(i,k)=max(0._r8,(nimax-ni(i,k))/deltat)
!         end if
!
!         ! get final values for precipitation q and N, based on
!         ! flux of precip from above, source/sink term, and terminal fallspeed
!         ! see eq. 15-16 in MG2008
!
!         ! rain
!
!         if (qric(i,k).ge.qsmall) then
!            if (k.eq.top_lev) then
!               qric(i,k)=qrtend(i,k)*dz(i,k)/cldmax(i,k)/umr(k)
!               nric(i,k)=nrtend(i,k)*dz(i,k)/cldmax(i,k)/unr(k)
!            else
!               qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
!                    (rho(i,k)*dz(i,k)*qrtend(i,k)))/(umr(k)*rho(i,k)*cldmax(i,k))
!               nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
!                    (rho(i,k)*dz(i,k)*nrtend(i,k)))/(unr(k)*rho(i,k)*cldmax(i,k))
!
!            end if
!         else
!            qric(i,k)=0._r8
!            nric(i,k)=0._r8
!         end if
!
!         ! snow
!
!         if (qniic(i,k).ge.qsmall) then
!            if (k.eq.top_lev) then
!               qniic(i,k)=qnitend(i,k)*dz(i,k)/cldmax(i,k)/ums(k)
!               nsic(i,k)=nstend(i,k)*dz(i,k)/cldmax(i,k)/uns(k)
!            else
!               qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
!                    (rho(i,k)*dz(i,k)*qnitend(i,k)))/(ums(k)*rho(i,k)*cldmax(i,k))
!               nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
!                    (rho(i,k)*dz(i,k)*nstend(i,k)))/(uns(k)*rho(i,k)*cldmax(i,k))
!            end if
!         else
!            qniic(i,k)=0._r8
!            nsic(i,k)=0._r8
!         end if
!
!         ! calculate precipitation flux at surface
!         ! divide by density of water to get units of m/s
!
!         prect(i) = prect(i)+(qrtend(i,k)*dz(i,k)*rho(i,k)+&
!              qnitend(i,k)*dz(i,k)*rho(i,k))/rhow
!         preci(i) = preci(i)+qnitend(i,k)*dz(i,k)*rho(i,k)/rhow
!
!         ! convert rain rate from m/s to mm/hr
!
!         rainrt(i,k)=qric(i,k)*rho(i,k)*umr(k)/rhow*3600._r8*1000._r8
!
!         ! vertically-integrated precip source/sink terms (note: grid-averaged)
!
!         qrtot = max(qrtot+qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
!         qstot = max(qstot+qnitend(i,k)*dz(i,k)*rho(i,k),0._r8)
!         nrtot = max(nrtot+nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
!         nstot = max(nstot+nstend(i,k)*dz(i,k)*rho(i,k),0._r8)
!
!         ! calculate melting and freezing of precip
!
!         ! melt snow at +2 C
!
!         if (t(i,k)+tlat(i,k)/cpp*deltat > 275.15_r8) then
!            if (qstot > 0._r8) then
!
!               ! make sure melting snow doesn't reduce temperature below threshold
!               dum = -xlf/cpp*qstot/(dz(i,k)*rho(i,k))
!               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.275.15_r8) then
!                  dum = (t(i,k)+tlat(i,k)/cpp*deltat-275.15_r8)*cpp/xlf
!                  dum = dum/(xlf/cpp*qstot/(dz(i,k)*rho(i,k)))
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qric(i,k)=qric(i,k)+dum*qniic(i,k)
!               nric(i,k)=nric(i,k)+dum*nsic(i,k)
!               qniic(i,k)=(1._r8-dum)*qniic(i,k)
!               nsic(i,k)=(1._r8-dum)*nsic(i,k)
!               ! heating tendency 
!               tmp=-xlf*dum*qstot/(dz(i,k)*rho(i,k))
!               meltsdt(i,k)=meltsdt(i,k) + tmp
!
!               tlat(i,k)=tlat(i,k)+tmp
!               qrtot=qrtot+dum*qstot
!               nrtot=nrtot+dum*nstot
!               qstot=(1._r8-dum)*qstot
!               nstot=(1._r8-dum)*nstot
!               preci(i)=(1._r8-dum)*preci(i)
!            end if
!         end if
!
!         ! freeze all rain at -5C for Arctic
!
!         if (t(i,k)+tlat(i,k)/cpp*deltat < (tmelt - 5._r8)) then
!
!            if (qrtot > 0._r8) then
!
!               ! make sure freezing rain doesn't increase temperature above threshold
!               dum = xlf/cpp*qrtot/(dz(i,k)*rho(i,k))
!               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.(tmelt - 5._r8)) then
!                  dum = -(t(i,k)+tlat(i,k)/cpp*deltat-(tmelt-5._r8))*cpp/xlf
!                  dum = dum/(xlf/cpp*qrtot/(dz(i,k)*rho(i,k)))
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qniic(i,k)=qniic(i,k)+dum*qric(i,k)
!               nsic(i,k)=nsic(i,k)+dum*nric(i,k)
!               qric(i,k)=(1._r8-dum)*qric(i,k)
!               nric(i,k)=(1._r8-dum)*nric(i,k)
!               ! heating tendency 
!               tmp = xlf*dum*qrtot/(dz(i,k)*rho(i,k))
!               frzrdt(i,k)=frzrdt(i,k) + tmp
!
!               tlat(i,k)=tlat(i,k)+tmp
!               qstot=qstot+dum*qrtot
!               qrtot=(1._r8-dum)*qrtot
!               nstot=nstot+dum*nrtot
!               nrtot=(1._r8-dum)*nrtot
!               preci(i)=preci(i)+dum*(prect(i)-preci(i))
!            end if
!         end if
!
!         ! if rain/snow mix ratio is zero so should number concentration
!
!         if (qniic(i,k).lt.qsmall) then
!            qniic(i,k)=0._r8
!            nsic(i,k)=0._r8
!         end if
!
!         if (qric(i,k).lt.qsmall) then
!            qric(i,k)=0._r8
!            nric(i,k)=0._r8
!         end if
!
!         ! make sure number concentration is a positive number to avoid 
!         ! taking root of negative
!
!         nric(i,k)=max(nric(i,k),0._r8)
!         nsic(i,k)=max(nsic(i,k),0._r8)
!
!         !.......................................................................
!         ! get size distribution parameters for fallspeed calculations
!         !......................................................................
!         ! rain
!
!         if (qric(i,k).ge.qsmall) then
!            lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
!            n0r(k) = nric(i,k)*lamr(k)
!
!            ! check for slope
!            ! change lammax and lammin for rain and snow
!            ! adjust vars
!
!            if (lamr(k).lt.lamminr) then
!
!               lamr(k) = lamminr
!
!               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
!               nric(i,k) = n0r(k)/lamr(k)
!            else if (lamr(k).gt.lammaxr) then
!               lamr(k) = lammaxr
!               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
!               nric(i,k) = n0r(k)/lamr(k)
!            end if
!
!
!            ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
!
!            unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
!            umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(i,k))
!
!         else
!            lamr(k) = 0._r8
!            n0r(k) = 0._r8
!            umr(k)=0._r8
!            unr(k)=0._r8
!         end if
!
!         !calculate mean size of combined rain and snow
!
!         if (lamr(k).gt.0._r8) then
!            Artmp = n0r(k) * pi / (2._r8 * lamr(k)**3._r8)
!         else 
!            Artmp = 0._r8
!         endif
!
!         if (lamc(k).gt.0._r8) then
!            Actmp = cdist1(k) * pi * gamma(pgam(k)+3._r8)/(4._r8 * lamc(k)**2._r8)
!         else 
!            Actmp = 0._r8
!         endif
!
!         if (Actmp.gt.0_r8.or.Artmp.gt.0) then
!            rercld(i,k)=rercld(i,k) + 3._r8 *(qric(i,k) + qcic(i,k)) / (4._r8 * rhow * (Actmp + Artmp))
!            arcld(i,k)=arcld(i,k)+1._r8
!         endif
!
!         !......................................................................
!         ! snow
!
!         if (qniic(i,k).ge.qsmall) then
!            lams(k) = (cons6*cs*nsic(i,k)/ &
!                 qniic(i,k))**(1._r8/ds)
!            n0s(k) = nsic(i,k)*lams(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lams(k).lt.lammins) then
!               lams(k) = lammins
!               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
!               nsic(i,k) = n0s(k)/lams(k)
!
!            else if (lams(k).gt.lammaxs) then
!               lams(k) = lammaxs
!               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
!               nsic(i,k) = n0s(k)/lams(k)
!            end if
!
!            ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
!
!            ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(i,k))
!            uns(k) = min(asn(i,k)*cons7/lams(k)**bs,1.2_r8*rhof(i,k))
!
!         else
!            lams(k) = 0._r8
!            n0s(k) = 0._r8
!            ums(k) = 0._r8
!            uns(k) = 0._r8
!         end if
!
!         !c........................................................................
!         ! sum over sub-step for average process rates
!
!         ! convert rain/snow q and N for output to history, note, 
!         ! output is for gridbox average
!
!         qrout(i,k)=qrout(i,k)+qric(i,k)*cldmax(i,k)
!         qsout(i,k)=qsout(i,k)+qniic(i,k)*cldmax(i,k)
!         nrout(i,k)=nrout(i,k)+nric(i,k)*rho(i,k)*cldmax(i,k)
!         nsout(i,k)=nsout(i,k)+nsic(i,k)*rho(i,k)*cldmax(i,k)
!
!         tlat1(i,k)=tlat1(i,k)+tlat(i,k)
!         qvlat1(i,k)=qvlat1(i,k)+qvlat(i,k)
!         qctend1(i,k)=qctend1(i,k)+qctend(i,k)
!         qitend1(i,k)=qitend1(i,k)+qitend(i,k)
!         nctend1(i,k)=nctend1(i,k)+nctend(i,k)
!         nitend1(i,k)=nitend1(i,k)+nitend(i,k)
!
!         t(i,k)=t(i,k)+tlat(i,k)*deltat/cpp
!         q(i,k)=q(i,k)+qvlat(i,k)*deltat
!         qc(i,k)=qc(i,k)+qctend(i,k)*deltat
!         qi(i,k)=qi(i,k)+qitend(i,k)*deltat
!         nc(i,k)=nc(i,k)+nctend(i,k)*deltat
!         ni(i,k)=ni(i,k)+nitend(i,k)*deltat
!
!         rainrt1(i,k)=rainrt1(i,k)+rainrt(i,k)
!
!         !divide rain radius over substeps for average
!         if (arcld(i,k) .gt. 0._r8) then
!            rercld(i,k)=rercld(i,k)/arcld(i,k)
!         end if
!
!         !calculate precip fluxes and adding them to summing sub-stepping variables
!         !! flux is zero at top interface
!         rflx(i,1)=0.0_r8
!         sflx(i,1)=0.0_r8
!
!         !! calculating the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
!         rflx(i,k+1)=qrout(i,k)*rho(i,k)*umr(k)
!         sflx(i,k+1)=qsout(i,k)*rho(i,k)*ums(k)
!
!         !! add to summing sub-stepping variable
!         rflx1(i,k+1)=rflx1(i,k+1)+rflx(i,k+1)
!         sflx1(i,k+1)=sflx1(i,k+1)+sflx(i,k+1)
!
!         !c........................................................................
!
!      end do ! k loop
!
!      prect1(i)=prect1(i)+prect(i)
!      preci1(i)=preci1(i)+preci(i)
!
!   end do ! it loop, sub-step
!
!   do k = top_lev, pver
!      rate1ord_cw2pr_st(i,k) = qcsinksum_rate1ord(k)/max(qcsum_rate1ord(k),1.0e-30_r8) 
!   end do
!
!300 continue  ! continue if no cloud water
!end do ! i loop
!
!! convert dt from sub-step back to full time step
!deltat=deltat*real(iter)
!
!!c.............................................................................
!
!do i=1,ncol
!
!   ! skip all calculations if no cloud water
!   if (ltrue(i).eq.0) then
!
!      do k=1,top_lev-1
!         ! assign zero values for effective radius above 1 mbar
!         effc(i,k)=0._r8
!         effi(i,k)=0._r8
!         effc_fn(i,k)=0._r8
!         lamcrad(i,k)=0._r8
!         pgamrad(i,k)=0._r8
!         deffi(i,k)=0._r8
!      end do
!
!      do k=top_lev,pver
!         ! assign default values for effective radius
!         effc(i,k)=10._r8
!         effi(i,k)=25._r8
!         effc_fn(i,k)=10._r8
!         lamcrad(i,k)=0._r8
!         pgamrad(i,k)=0._r8
!         deffi(i,k)=0._r8
!      end do
!      goto 500
!   end if
!
!   ! initialize nstep for sedimentation sub-steps
!   nstep = 1
!
!   ! divide precip rate by number of sub-steps to get average over time step
!
!   prect(i)=prect1(i)/real(iter)
!   preci(i)=preci1(i)/real(iter)
!
!   do k=top_lev,pver
!
!      ! assign variables back to start-of-timestep values before updating after sub-steps 
!
!      t(i,k)=t1(i,k)
!      q(i,k)=q1(i,k)
!      qc(i,k)=qc1(i,k)
!      qi(i,k)=qi1(i,k)
!      nc(i,k)=nc1(i,k)
!      ni(i,k)=ni1(i,k)
!
!      ! divide microphysical tendencies by number of sub-steps to get average over time step
!
!      tlat(i,k)=tlat1(i,k)/real(iter)
!      qvlat(i,k)=qvlat1(i,k)/real(iter)
!      qctend(i,k)=qctend1(i,k)/real(iter)
!      qitend(i,k)=qitend1(i,k)/real(iter)
!      nctend(i,k)=nctend1(i,k)/real(iter)
!      nitend(i,k)=nitend1(i,k)/real(iter)
!
!      rainrt(i,k)=rainrt1(i,k)/real(iter)
!
!      ! divide by number of sub-steps to find final values
!      rflx(i,k+1)=rflx1(i,k+1)/real(iter)
!      sflx(i,k+1)=sflx1(i,k+1)/real(iter)
!
!      ! divide output precip q and N by number of sub-steps to get average over time step
!
!      qrout(i,k)=qrout(i,k)/real(iter)
!      qsout(i,k)=qsout(i,k)/real(iter)
!      nrout(i,k)=nrout(i,k)/real(iter)
!      nsout(i,k)=nsout(i,k)/real(iter)
!
!      ! divide trop_mozart variables by number of sub-steps to get average over time step 
!
!      nevapr(i,k) = nevapr(i,k)/real(iter)
!      nevapr2(i,k) = nevapr2(i,k)/real(iter)
!      evapsnow(i,k) = evapsnow(i,k)/real(iter)
!      prain(i,k) = prain(i,k)/real(iter)
!      prodsnow(i,k) = prodsnow(i,k)/real(iter)
!      cmeout(i,k) = cmeout(i,k)/real(iter)
!
!      cmeiout(i,k) = cmeiout(i,k)/real(iter)
!      meltsdt(i,k) = meltsdt(i,k)/real(iter)
!      frzrdt (i,k) = frzrdt (i,k)/real(iter)
!
!
!      ! microphysics output
!      prao(i,k)=prao(i,k)/real(iter)
!      prco(i,k)=prco(i,k)/real(iter)
!      mnuccco(i,k)=mnuccco(i,k)/real(iter)
!      mnuccto(i,k)=mnuccto(i,k)/real(iter)
!      msacwio(i,k)=msacwio(i,k)/real(iter)
!      psacwso(i,k)=psacwso(i,k)/real(iter)
!      bergso(i,k)=bergso(i,k)/real(iter)
!      bergo(i,k)=bergo(i,k)/real(iter)
!      prcio(i,k)=prcio(i,k)/real(iter)
!      praio(i,k)=praio(i,k)/real(iter)
!
!      mnuccro(i,k)=mnuccro(i,k)/real(iter)
!      pracso (i,k)=pracso (i,k)/real(iter)
!
!      mnuccdo(i,k)=mnuccdo(i,k)/real(iter)
!
!      ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
!      nevapr(i,k) = nevapr(i,k) + evapsnow(i,k)
!      prer_evap(i,k) = nevapr2(i,k)
!      prain(i,k) = prain(i,k) + prodsnow(i,k)
!
!      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      ! calculate sedimentation for cloud water and ice
!      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      ! update in-cloud cloud mixing ratio and number concentration 
!      ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
!      ! note: these are in-cloud values***, hence we divide by cloud fraction
!
!      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)
!      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)/icldm(i,k)
!      dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k),0._r8)
!      dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)/icldm(i,k),0._r8)
!
!      ! obtain new slope parameter to avoid possible singularity
!
!      if (dumi(i,k).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)
!
!         lami(k) = (cons1*ci* &
!              dumni(i,k)/dumi(i,k))**(1._r8/di)
!         lami(k)=max(lami(k),lammini)
!         lami(k)=min(lami(k),lammaxi)
!      else
!         lami(k)=0._r8
!      end if
!
!      if (dumc(i,k).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)
!         ! add lower limit to in-cloud number concentration
!         dumnc(i,k)=max(dumnc(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3 
!         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)**2)-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
!              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         lammax = (pgam(k)+1._r8)/2.e-6_r8
!         lamc(k)=max(lamc(k),lammin)
!         lamc(k)=min(lamc(k),lammax)
!      else
!         lamc(k)=0._r8
!      end if
!
!      ! calculate number and mass weighted fall velocity for droplets
!      ! include effects of sub-grid distribution of cloud water
!
!
!      if (dumc(i,k).ge.qsmall) then
!         unc = acn(i,k)*gamma(1._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+1._r8))
!         umc = acn(i,k)*gamma(4._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+4._r8))
!         ! fallspeed for output
!         vtrmc(i,k)=umc
!      else
!         umc = 0._r8
!         unc = 0._r8
!      end if
!
!      ! calculate number and mass weighted fall velocity for cloud ice
!
!      if (dumi(i,k).ge.qsmall) then
!         uni =  ain(i,k)*cons16/lami(k)**bi
!         umi = ain(i,k)*cons17/(6._r8*lami(k)**bi)
!         uni=min(uni,1.2_r8*rhof(i,k))
!         umi=min(umi,1.2_r8*rhof(i,k))
!
!         ! fallspeed
!         vtrmi(i,k)=umi
!      else
!         umi = 0._r8
!         uni = 0._r8
!      end if
!
!      fi(k) = g*rho(i,k)*umi
!      fni(k) = g*rho(i,k)*uni
!      fc(k) = g*rho(i,k)*umc
!      fnc(k) = g*rho(i,k)*unc
!
!      ! calculate number of split time steps to ensure courant stability criteria
!      ! for sedimentation calculations
!
!      rgvm = max(fi(k),fc(k),fni(k),fnc(k))
!      nstep = max(int(rgvm*deltat/pdel(i,k)+1._r8),nstep)
!
!      ! redefine dummy variables - sedimentation is calculated over grid-scale
!      ! quantities to ensure conservation
!
!      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
!      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
!      dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
!      dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)
!
!      if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
!      if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
!
!   end do       !!! vertical loop
!   do n = 1,nstep  !! loop over sub-time step to ensure stability
!
!      do k = top_lev,pver
!         if (do_cldice) then
!            falouti(k) = fi(k)*dumi(i,k)
!            faloutni(k) = fni(k)*dumni(i,k)
!         else
!            falouti(k)  = 0._r8
!            faloutni(k) = 0._r8
!         end if
!
!         faloutc(k) = fc(k)*dumc(i,k)
!         faloutnc(k) = fnc(k)*dumnc(i,k)
!      end do
!
!      ! top of model
!
!      k = top_lev
!      faltndi = falouti(k)/pdel(i,k)
!      faltndni = faloutni(k)/pdel(i,k)
!      faltndc = faloutc(k)/pdel(i,k)
!      faltndnc = faloutnc(k)/pdel(i,k)
!
!      ! add fallout terms to microphysical tendencies
!
!      qitend(i,k) = qitend(i,k)-faltndi/nstep
!      nitend(i,k) = nitend(i,k)-faltndni/nstep
!      qctend(i,k) = qctend(i,k)-faltndc/nstep
!      nctend(i,k) = nctend(i,k)-faltndnc/nstep
!
!      ! sedimentation tendencies for output
!      qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
!      qisedten(i,k)=qisedten(i,k)-faltndi/nstep
!
!      dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
!      dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
!      dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
!      dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep
!
!      do k = top_lev+1,pver
!
!         ! for cloud liquid and ice, if cloud fraction increases with height
!         ! then add flux from above to both vapor and cloud water of current level
!         ! this means that flux entering clear portion of cell from above evaporates
!         ! instantly
!
!         dum=lcldm(i,k)/lcldm(i,k-1)
!         dum=min(dum,1._r8)
!         dum1=icldm(i,k)/icldm(i,k-1)
!         dum1=min(dum1,1._r8)
!
!         faltndqie=(falouti(k)-falouti(k-1))/pdel(i,k)
!         faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(i,k)
!         faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(i,k)
!         faltndqce=(faloutc(k)-faloutc(k-1))/pdel(i,k)
!         faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(i,k)
!         faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(i,k)
!
!         ! add fallout terms to eulerian tendencies
!
!         qitend(i,k) = qitend(i,k)-faltndi/nstep
!         nitend(i,k) = nitend(i,k)-faltndni/nstep
!         qctend(i,k) = qctend(i,k)-faltndc/nstep
!         nctend(i,k) = nctend(i,k)-faltndnc/nstep
!
!         ! sedimentation tendencies for output
!         qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
!         qisedten(i,k)=qisedten(i,k)-faltndi/nstep
!
!         ! add terms to to evap/sub of cloud water
!
!         qvlat(i,k)=qvlat(i,k)-(faltndqie-faltndi)/nstep
!         ! for output
!         qisevap(i,k)=qisevap(i,k)-(faltndqie-faltndi)/nstep
!         qvlat(i,k)=qvlat(i,k)-(faltndqce-faltndc)/nstep
!         ! for output
!         qcsevap(i,k)=qcsevap(i,k)-(faltndqce-faltndc)/nstep
!
!         tlat(i,k)=tlat(i,k)+(faltndqie-faltndi)*xxls/nstep
!         tlat(i,k)=tlat(i,k)+(faltndqce-faltndc)*xxlv/nstep
!
!         dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
!         dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
!         dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
!         dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep
!
!         Fni(K)=MAX(Fni(K)/pdel(i,K),Fni(K-1)/pdel(i,K-1))*pdel(i,K)
!         FI(K)=MAX(FI(K)/pdel(i,K),FI(K-1)/pdel(i,K-1))*pdel(i,K)
!         fnc(k)=max(fnc(k)/pdel(i,k),fnc(k-1)/pdel(i,k-1))*pdel(i,k)
!         Fc(K)=MAX(Fc(K)/pdel(i,K),Fc(K-1)/pdel(i,K-1))*pdel(i,K)
!
!      end do   !! k loop
!
!      ! units below are m/s
!      ! cloud water/ice sedimentation flux at surface 
!      ! is added to precip flux at surface to get total precip (cloud + precip water)
!      ! rate
!
!      prect(i) = prect(i)+(faloutc(pver)+falouti(pver))/g/nstep/1000._r8  
!      preci(i) = preci(i)+(falouti(pver))/g/nstep/1000._r8
!
!   end do   !! nstep loop
!
!   ! end sedimentation
!   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   ! get new update for variables that includes sedimentation tendency
!   ! note : here dum variables are grid-average, NOT in-cloud
!
!   do k=top_lev,pver
!
!      dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)
!      dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)
!      dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)
!      dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)
!
!      if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
!      if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
!
!      ! calculate instantaneous processes (melting, homogeneous freezing)
!      if (do_cldice) then
!
!         if (t(i,k)+tlat(i,k)/cpp*deltat > tmelt) then
!            if (dumi(i,k) > 0._r8) then
!
!               ! limit so that melting does not push temperature below freezing
!               dum = -dumi(i,k)*xlf/cpp
!               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.tmelt) then
!                  dum = (t(i,k)+tlat(i,k)/cpp*deltat-tmelt)*cpp/xlf
!                  dum = dum/dumi(i,k)*xlf/cpp 
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qctend(i,k)=qctend(i,k)+dum*dumi(i,k)/deltat
!
!               ! for output
!               melto(i,k)=dum*dumi(i,k)/deltat
!
!               ! assume melting ice produces droplet
!               ! mean volume radius of 8 micron
!
!               nctend(i,k)=nctend(i,k)+3._r8*dum*dumi(i,k)/deltat/ &
!                    (4._r8*pi*5.12e-16_r8*rhow)
!
!               qitend(i,k)=((1._r8-dum)*dumi(i,k)-qi(i,k))/deltat
!               nitend(i,k)=((1._r8-dum)*dumni(i,k)-ni(i,k))/deltat
!               tlat(i,k)=tlat(i,k)-xlf*dum*dumi(i,k)/deltat
!            end if
!         end if
!
!         ! homogeneously freeze droplets at -40 C
!
!         if (t(i,k)+tlat(i,k)/cpp*deltat < 233.15_r8) then
!            if (dumc(i,k) > 0._r8) then
!
!               ! limit so that freezing does not push temperature above threshold
!               dum = dumc(i,k)*xlf/cpp
!               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.233.15_r8) then
!                  dum = -(t(i,k)+tlat(i,k)/cpp*deltat-233.15_r8)*cpp/xlf
!                  dum = dum/dumc(i,k)*xlf/cpp
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qitend(i,k)=qitend(i,k)+dum*dumc(i,k)/deltat
!               ! for output
!               homoo(i,k)=dum*dumc(i,k)/deltat
!
!               ! assume 25 micron mean volume radius of homogeneously frozen droplets
!               ! consistent with size of detrained ice in stratiform.F90
!               nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*1.563e-14_r8* &
!                    500._r8)/deltat
!               qctend(i,k)=((1._r8-dum)*dumc(i,k)-qc(i,k))/deltat
!               nctend(i,k)=((1._r8-dum)*dumnc(i,k)-nc(i,k))/deltat
!               tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)/deltat
!            end if
!         end if
!
!         ! remove any excess over-saturation, which is possible due to non-linearity when adding 
!         ! together all microphysical processes
!         ! follow code similar to old CAM scheme
!
!         qtmp=q(i,k)+qvlat(i,k)*deltat
!         ttmp=t(i,k)+tlat(i,k)/cpp*deltat
!
!         esn = svp_water(ttmp)  ! use rhw to allow ice supersaturation
!         qsn = svp_to_qsat(esn, p(i,k))
!
!         if (qtmp > qsn .and. qsn > 0) then
!            ! expression below is approximate since there may be ice deposition
!            dum = (qtmp-qsn)/(1._r8+cons27*qsn/(cpp*rv*ttmp**2))/deltat
!            ! add to output cme
!            cmeout(i,k) = cmeout(i,k)+dum
!            ! now add to tendencies, partition between liquid and ice based on temperature
!            if (ttmp > 268.15_r8) then
!               dum1=0.0_r8
!               ! now add to tendencies, partition between liquid and ice based on te
!            else if (ttmp < 238.15_r8) then
!               dum1=1.0_r8
!            else
!               dum1=(268.15_r8-ttmp)/30._r8
!            end if
!
!            dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
!                 *qsn/(cpp*rv*ttmp**2))/deltat
!            qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
!            ! for output
!            qcreso(i,k)=dum*(1._r8-dum1)
!            qitend(i,k)=qitend(i,k)+dum*dum1
!            qireso(i,k)=dum*dum1
!            qvlat(i,k)=qvlat(i,k)-dum
!            ! for output
!            qvres(i,k)=-dum
!            tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
!         end if
!      end if
!
!      !...............................................................................
!      ! calculate effective radius for pass to radiation code
!      ! if no cloud water, default value is 10 micron for droplets,
!      ! 25 micron for cloud ice
!
!      ! update cloud variables after instantaneous processes to get effective radius
!      ! variables are in-cloud to calculate size dist parameters
!
!      dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/lcldm(i,k)
!      dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/icldm(i,k)
!      dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/lcldm(i,k)
!      dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/icldm(i,k)
!
!      ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
!
!      dumc(i,k)=min(dumc(i,k),5.e-3_r8)
!      dumi(i,k)=min(dumi(i,k),5.e-3_r8)
!
!      !...................
!      ! cloud ice effective radius
!
!      if (dumi(i,k).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)
!         lami(k) = (cons1*ci*dumni(i,k)/dumi(i,k))**(1._r8/di)
!
!         if (lami(k).lt.lammini) then
!            lami(k) = lammini
!            n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*cons1)
!            niic(i,k) = n0i(k)/lami(k)
!            ! adjust number conc if needed to keep mean size in reasonable range
!            if (do_cldice) nitend(i,k)=(niic(i,k)*icldm(i,k)-ni(i,k))/deltat
!
!         else if (lami(k).gt.lammaxi) then
!            lami(k) = lammaxi
!            n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*cons1)
!            niic(i,k) = n0i(k)/lami(k)
!            ! adjust number conc if needed to keep mean size in reasonable range
!            if (do_cldice) nitend(i,k)=(niic(i,k)*icldm(i,k)-ni(i,k))/deltat
!         end if
!         effi(i,k) = 1.5_r8/lami(k)*1.e6_r8
!
!      else
!         effi(i,k) = 25._r8
!      end if
!
!      ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
!      ! radius has already been determined from the size distribution.
!      if (.not. do_cldice) then
!         effi(i,k) = re_ice(i,k) * 1e6_r8      ! m -> um
!      end if
!
!      !...................
!      ! cloud droplet effective radius
!
!      if (dumc(i,k).ge.qsmall) then
!
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! set tendency to ensure minimum droplet concentration
!         ! after update by microphysics, except when lambda exceeds bounds on mean drop
!         ! size or if there is no cloud water
!         if (dumnc(i,k).lt.cdnl/rho(i,k)) then   
!            nctend(i,k)=(cdnl/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat   
!         end if
!         dumnc(i,k)=max(dumnc(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)**2)-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
!              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         ! Multiply by omsm to fit within RRTMG's table.
!         lammax = (pgam(k)+1._r8)*omsm/2.e-6_r8
!         if (lamc(k).lt.lammin) then
!            lamc(k) = lammin
!            ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
!                 gamma(pgam(k)+1._r8)/ &
!                 (pi*rhow*gamma(pgam(k)+4._r8))
!            ! adjust number conc if needed to keep mean size in reasonable range
!            nctend(i,k)=(ncic(i,k)*lcldm(i,k)-nc(i,k))/deltat
!
!         else if (lamc(k).gt.lammax) then
!            lamc(k) = lammax
!            ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
!                 gamma(pgam(k)+1._r8)/ &
!                 (pi*rhow*gamma(pgam(k)+4._r8))
!            ! adjust number conc if needed to keep mean size in reasonable range
!            nctend(i,k)=(ncic(i,k)*lcldm(i,k)-nc(i,k))/deltat
!         end if
!
!         effc(i,k) = &
!              gamma(pgam(k)+4._r8)/ &
!              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
!         !assign output fields for shape here
!         lamcrad(i,k)=lamc(k)
!         pgamrad(i,k)=pgam(k)
!
!      else
!         effc(i,k) = 10._r8
!         lamcrad(i,k)=0._r8
!         pgamrad(i,k)=0._r8
!      end if
!
!      ! ice effective diameter for david mitchell's optics
!      if (do_cldice) then
!         deffi(i,k)=effi(i,k)*rhoi/917._r8*2._r8
!      else
!         deffi(i,k)=effi(i,k) * 2._r8
!      end if
!
!
!!!! recalculate effective radius for constant number, in order to separate
!      ! first and second indirect effects
!      ! assume constant number of 10^8 kg-1
!
!      dumnc(i,k)=1.e8_r8
!
!      if (dumc(i,k).ge.qsmall) then
!         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)**2)-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
!              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         lammax = (pgam(k)+1._r8)/2.e-6_r8
!         if (lamc(k).lt.lammin) then
!            lamc(k) = lammin
!         else if (lamc(k).gt.lammax) then
!            lamc(k) = lammax
!         end if
!         effc_fn(i,k) = &
!              gamma(pgam(k)+4._r8)/ &
!              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
!
!      else
!         effc_fn(i,k) = 10._r8
!      end if
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!
!
!   end do ! vertical k loop
!
!500 continue
!
!   do k=top_lev,pver
!      ! if updated q (after microphysics) is zero, then ensure updated n is also zero
!
!      if (qc(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-nc(i,k)/deltat
!      if (do_cldice .and. qi(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-ni(i,k)/deltat
!   end do
!
!end do ! i loop
!
!! add snow ouptut
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
!         dsout(i,k)=3._r8*rhosn/917._r8*(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
!      endif
!   end do
!end do
!
!!calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual
!do i = 1,ncol
!   do k=top_lev,pver
!      !! RAIN
!      if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
!         reff_rain(i,k)=1.5_r8*(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)*1.e6_r8
!      endif
!      !! SNOW
!      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
!         reff_snow(i,k)=1.5_r8*(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)*1.e6_r8
!      end if
!   end do
!end do
!
!! analytic radar reflectivity
!! formulas from Matthew Shupe, NOAA/CERES
!! *****note: radar reflectivity is local (in-precip average)
!! units of mm^6/m^3
!
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qc(i,k)+qctend(i,k)*deltat.ge.qsmall) then
!         dum=((qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)*1000._r8)**2 &
!              /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)/1.e6_r8)*lcldm(i,k)/cldmax(i,k)
!      else
!         dum=0._r8
!      end if
!      if (qi(i,k)+qitend(i,k)*deltat.ge.qsmall) then
!         dum1=((qi(i,k)+qitend(i,k)*deltat)*rho(i,k)/icldm(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(i,k)/cldmax(i,k)
!      else 
!         dum1=0._r8
!      end if
!
!      if (qsout(i,k).ge.qsmall) then
!         dum1=dum1+(qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
!      end if
!
!      refl(i,k)=dum+dum1
!
!      ! add rain rate, but for 37 GHz formulation instead of 94 GHz
!      ! formula approximated from data of Matrasov (2007)
!      ! rainrt is the rain rate in mm/hr
!      ! reflectivity (dum) is in DBz
!
!      if (rainrt(i,k).ge.0.001_r8) then
!         dum=log10(rainrt(i,k)**6._r8)+16._r8
!
!         ! convert from DBz to mm^6/m^3
!
!         dum = 10._r8**(dum/10._r8)
!      else
!         ! don't include rain rate in R calculation for values less than 0.001 mm/hr
!         dum=0._r8
!      end if
!
!      ! add to refl
!
!      refl(i,k)=refl(i,k)+dum
!
!      !output reflectivity in Z.
!      areflz(i,k)=refl(i,k)
!
!      ! convert back to DBz 
!
!      if (refl(i,k).gt.minrefl) then 
!         refl(i,k)=10._r8*log10(refl(i,k))
!      else
!         refl(i,k)=-9999._r8
!      end if
!
!      !set averaging flag
!      if (refl(i,k).gt.mindbz) then 
!         arefl(i,k)=refl(i,k)
!         frefl(i,k)=1.0_r8  
!      else
!         arefl(i,k)=0._r8
!         areflz(i,k)=0._r8
!         frefl(i,k)=0._r8
!      end if
!
!      ! bound cloudsat reflectivity
!
!      csrfl(i,k)=min(csmax,refl(i,k))
!
!      !set averaging flag
!      if (csrfl(i,k).gt.csmin) then 
!         acsrfl(i,k)=refl(i,k)
!         fcsrfl(i,k)=1.0_r8  
!      else
!         acsrfl(i,k)=0._r8
!         fcsrfl(i,k)=0._r8
!      end if
!
!   end do
!end do
!
!
!! averaging for snow and rain number and diameter
!
!qrout2(:,:)=0._r8
!qsout2(:,:)=0._r8
!nrout2(:,:)=0._r8
!nsout2(:,:)=0._r8
!drout2(:,:)=0._r8
!dsout2(:,:)=0._r8
!freqs(:,:)=0._r8
!freqr(:,:)=0._r8
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
!         qrout2(i,k)=qrout(i,k)
!         nrout2(i,k)=nrout(i,k)
!         drout2(i,k)=(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)
!         freqr(i,k)=1._r8
!      endif
!      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
!         qsout2(i,k)=qsout(i,k)
!         nsout2(i,k)=nsout(i,k)
!         dsout2(i,k)=(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
!         freqs(i,k)=1._r8
!      endif
!   end do
!end do
!
!! output activated liquid and ice (convert from #/kg -> #/m3)
!do i = 1,ncol
!   do k=top_lev,pver
!      ncai(i,k)=dum2i(i,k)*rho(i,k)
!      ncal(i,k)=dum2l(i,k)*rho(i,k)
!   end do
!end do
!
!
!!redefine fice here....
!nfice(:,:)=0._r8
!do k=top_lev,pver
!   do i=1,ncol
!      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
!      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
!      dumfice=qsout(i,k) + qrout(i,k) + dumc(i,k) + dumi(i,k)  
!
!      if (dumfice.gt.qsmall.and.(qsout(i,k)+dumi(i,k).gt.qsmall)) then
!         nfice(i,k)=(qsout(i,k) + dumi(i,k))/dumfice
!      endif
!
!      if (nfice(i,k).gt.1._r8) then
!         nfice(i,k)=1._r8
!      endif
!
!   enddo
!enddo
!call t_stopf("Mg_computing_mpe")

!! initialize  output fields for number conc qand ice nucleation
!ncai(1:ncol,1:pver)=0._r8 
!ncal(1:ncol,1:pver)=0._r8  
!
!!Initialize rain size
!rercld(1:ncol,1:pver)=0._r8
!arcld(1:ncol,1:pver)=0._r8
!
!!initialize radiation output variables
!pgamrad(1:ncol,1:pver)=0._r8 ! liquid gamma parameter for optics (radiation)
!lamcrad(1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics (radiation)
!deffi  (1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics (radiation)
!!initialize radiation output variables
!!initialize water vapor tendency term output
!qcsevap(1:ncol,1:pver)=0._r8 
!qisevap(1:ncol,1:pver)=0._r8 
!qvres  (1:ncol,1:pver)=0._r8 
!cmeiout (1:ncol,1:pver)=0._r8
!vtrmc (1:ncol,1:pver)=0._r8
!vtrmi (1:ncol,1:pver)=0._r8
!qcsedten (1:ncol,1:pver)=0._r8
!qisedten (1:ncol,1:pver)=0._r8    
!
!prao(1:ncol,1:pver)=0._r8 
!prco(1:ncol,1:pver)=0._r8 
!mnuccco(1:ncol,1:pver)=0._r8 
!mnuccto(1:ncol,1:pver)=0._r8 
!msacwio(1:ncol,1:pver)=0._r8 
!psacwso(1:ncol,1:pver)=0._r8 
!bergso(1:ncol,1:pver)=0._r8 
!bergo(1:ncol,1:pver)=0._r8 
!melto(1:ncol,1:pver)=0._r8 
!homoo(1:ncol,1:pver)=0._r8 
!qcreso(1:ncol,1:pver)=0._r8 
!prcio(1:ncol,1:pver)=0._r8 
!praio(1:ncol,1:pver)=0._r8 
!qireso(1:ncol,1:pver)=0._r8 
!mnuccro(1:ncol,1:pver)=0._r8 
!pracso (1:ncol,1:pver)=0._r8 
!meltsdt(1:ncol,1:pver)=0._r8
!frzrdt (1:ncol,1:pver)=0._r8
!mnuccdo(1:ncol,1:pver)=0._r8
!
!rflx(:,:)=0._r8
!sflx(:,:)=0._r8
!effc(:,:)=0._r8
!effc_fn(:,:)=0._r8
!effi(:,:)=0._r8
!
!! assign variable deltat for sub-stepping...
!deltat=deltatin
!
!! parameters for scheme
!
!omsm=0.99999_r8
!dto2=0.5_r8*deltat
!mincld=0.0001_r8
!
!! initialize multi-level fields
!q(1:ncol,1:pver)=qn(1:ncol,1:pver)
!t(1:ncol,1:pver)=tn(1:ncol,1:pver)
!
!! initialize time-varying parameters
!
!do k=1,pver
!   do i=1,ncol
!      rho(i,k)=p(i,k)/(r*t(i,k))
!
!      tmp_array(1) = t(i,k)
!      tmp_array(2) = 1.81_r8
!      call athread_spawn(slave_pow_parallel, tmp_array)
!      call athread_join()
!      tmp1 = tmp_array(1)
!
!      tmp_array(1) = t(i,k)
!      tmp_array(2) = 1.5_r8
!      call athread_spawn(slave_pow_parallel, tmp_array)
!      call athread_join()
!      tmp2 = tmp_array(1)
!
!      dv(i,k) = 8.794E-5_r8*tmp1/p(i,k)
!      mu(i,k) = 1.496E-6_r8*tmp2/(t(i,k)+120._r8) 
!      sc(i,k) = mu(i,k)/(rho(i,k)*dv(i,k))
!      kap(i,k) = 1.414e3_r8*1.496e-6_r8*tmp2/(t(i,k)+120._r8) 
!
!      !dv(i,k) = 8.794E-5_r8*t(i,k)**1.81_r8/p(i,k)
!      !mu(i,k) = 1.496E-6_r8*t(i,k)**1.5_r8/(t(i,k)+120._r8) 
!      !sc(i,k) = mu(i,k)/(rho(i,k)*dv(i,k))
!      !kap(i,k) = 1.414e3_r8*1.496e-6_r8*t(i,k)**1.5_r8/(t(i,k)+120._r8) 
!
!      ! air density adjustment for fallspeed parameters
!      ! includes air density correction factor to the
!      ! power of 0.54 following Heymsfield and Bansemer 2007
!
!      tmp_array(1) = rhosu/rho(i,k)
!      tmp_array(2) = 0.54_r8
!      call athread_spawn(slave_pow_parallel, tmp_array)
!      call athread_join()
!      tmp1 = tmp_array(1)
!      rhof(i,k) = tmp1
!
!      !rhof(i,k)=(rhosu/rho(i,k))**0.54_r8
!
!      arn(i,k)=ar*rhof(i,k)
!      asn(i,k)=as*rhof(i,k)
!      acn(i,k)=ac*rhof(i,k)
!      ain(i,k)=ai*rhof(i,k)
!
!      ! get dz from dp and hydrostatic approx
!      ! keep dz positive (define as layer k-1 - layer k)
!
!      dz(i,k)= pdel(i,k)/(rho(i,k)*g)
!
!      !! Check
!      !qrout2(i,k) = rho(i,k)
!      !qsout2(i,k) = dv(i,k)
!      !nrout2(i,k) = mu(i,k)
!      !nsout2(i,k) = sc(i,k)
!      !drout2(i,k) = kap(i,k)
!      !dsout2(i,k) = dz(i,k)
!
!   end do
!end do
!
!! initialization
!qc(1:ncol,1:top_lev-1) = 0._r8
!qi(1:ncol,1:top_lev-1) = 0._r8
!nc(1:ncol,1:top_lev-1) = 0._r8
!ni(1:ncol,1:top_lev-1) = 0._r8
!t1(1:ncol,1:pver) = t(1:ncol,1:pver)
!q1(1:ncol,1:pver) = q(1:ncol,1:pver)
!qc1(1:ncol,1:pver) = qc(1:ncol,1:pver)
!qi1(1:ncol,1:pver) = qi(1:ncol,1:pver)
!nc1(1:ncol,1:pver) = nc(1:ncol,1:pver)
!ni1(1:ncol,1:pver) = ni(1:ncol,1:pver)
!
!! initialize tendencies to zero
!tlat1(1:ncol,1:pver)=0._r8
!qvlat1(1:ncol,1:pver)=0._r8
!qctend1(1:ncol,1:pver)=0._r8
!qitend1(1:ncol,1:pver)=0._r8
!nctend1(1:ncol,1:pver)=0._r8
!nitend1(1:ncol,1:pver)=0._r8
!
!! initialize precip output
!qrout(1:ncol,1:pver)=0._r8
!qsout(1:ncol,1:pver)=0._r8
!nrout(1:ncol,1:pver)=0._r8
!nsout(1:ncol,1:pver)=0._r8
!dsout(1:ncol,1:pver)=0._r8
!
!drout(1:ncol,1:pver)=0._r8
!
!reff_rain(1:ncol,1:pver)=0._r8
!reff_snow(1:ncol,1:pver)=0._r8
!
!! initialize variables for trop_mozart
!nevapr(1:ncol,1:pver) = 0._r8
!nevapr2(1:ncol,1:pver) = 0._r8
!evapsnow(1:ncol,1:pver) = 0._r8
!prain(1:ncol,1:pver) = 0._r8
!prodsnow(1:ncol,1:pver) = 0._r8
!cmeout(1:ncol,1:pver) = 0._r8
!
!am_evp_st(1:ncol,1:pver) = 0._r8
!
!! for refl calc
!rainrt1(1:ncol,1:pver) = 0._r8
!
!! initialize precip fraction and output tendencies
!cldmax(1:ncol,1:pver)=mincld
!
!!initialize aerosol number
!!        naer2(1:ncol,1:pver,:)=0._r8
!dum2l(1:ncol,1:pver)=0._r8
!dum2i(1:ncol,1:pver)=0._r8
!
!! initialize avg precip rate
!prect1(1:ncol)=0._r8
!preci1(1:ncol)=0._r8
!
!lami(:) = 0._r8
!n0i(:) = 0._r8
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!Get humidity and saturation vapor pressures
!
!do k=top_lev,pver
!
!   do i=1,ncol
!
!      ! find wet bulk temperature and saturation value for provisional t and q without
!      ! condensation
!      
!      tmp1 = svp_tboil/t(i,k)
!      call athread_spawn(slave_log10_parallel, tmp1)
!      call athread_join()
!
!      es(i) = 10._r8**(-7.90298_r8*(svp_tboil/t(i,k)-1._r8)+ &
!           5.02808_r8*tmp1- &
!           1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t(i,k)/svp_tboil))-1._r8)+ &
!           8.1328e-3_r8*(10._r8**(-3.49149_r8*(svp_tboil/t(i,k)-1._r8))-1._r8)+ &
!           log10(1013.246_r8))*100._r8
!
!      !es(i) = svp_water(t(i,k))
!      qs(i) = svp_to_qsat(es(i), p(i,k))
!
!      ! Prevents negative values.
!      if (qs(i) < 0.0_r8) then
!         qs(i) = 1.0_r8
!         es(i) = p(i,k)
!      end if
!
!      !esl(i,k)=svp_water(t(i,k))
!      !esi(i,k)=svp_ice(t(i,k))
!
!
!      esl(i,k) = 10._r8**(-7.90298_r8*(svp_tboil/t(i,k)-1._r8)+ &
!           5.02808_r8*tmp1- &
!           1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t(i,k)/svp_tboil))-1._r8)+ &
!           8.1328e-3_r8*(10._r8**(-3.49149_r8*(svp_tboil/t(i,k)-1._r8))-1._r8)+ &
!           log10(1013.246_r8))*100._r8
!
!      tmp2 = svp_h2otrip/t(i,k)
!      call athread_spawn(slave_log10_parallel, tmp2)
!      call athread_join()
!
!      esi(i,k) = 10.d0**(-9.09718d0*(svp_h2otrip/t(i,k)-1.d0)-3.56654d0* &
!           tmp2+0.876793d0*(1.d0-t(i,k)/svp_h2otrip)+ &
!           log10(6.1071d0))*100.d0
!      
!      !qrout2(i,k) = esl(i,k)
!
!
!
!            !qrout2(i,k) = log10(t(i,k))
!            !qsout2(i,k) = 10.d0**(11.344d0*(1.d0-t(i,k)/svp_tboil))
!            !nrout2(i,k) = 10.d0**(-3.49149d0*(svp_tboil/t(i,k)-1.d0))
!            !nsout2(i,k) = log10(1013.246d0)
!            !drout2(i,k) = (-7.90298d0*(svp_tboil/t(i,k)-1.d0)+ &
!            !   5.02808d0*log10(svp_tboil/t(i,k))- &
!            !   1.3816d-7*(10.d0**(11.344d0*(1.d0-t(i,k)/svp_tboil))-1.d0)+ &
!            !   8.1328d-3*(10.d0**(-3.49149d0*(svp_tboil/t(i,k)-1.d0))-1.d0)+ &
!            !   log10(1013.246d0))
!            !!dsout2(i,k) = log10(svp_tboil)
!            !dsout2(i,k) = t(i,k)
!
!
!      ! hm fix, make sure when above freezing that esi=esl, not active yet
!      if (t(i,k).gt.tmelt)esi(i,k)=esl(i,k)
!
!      relhum(i,k)=q(i,k)/qs(i)
!
!      ! get cloud fraction, check for minimum
!
!      cldm(i,k)=max(cldn(i,k),mincld)
!      cldmw(i,k)=max(cldn(i,k),mincld)
!
!      icldm(i,k)=max(icecldf(i,k),mincld)
!      lcldm(i,k)=max(liqcldf(i,k),mincld)
!
!      ! subcolumns, set cloud fraction variables to one
!      ! if cloud water or ice is present, if not present
!      ! set to mincld (mincld used instead of zero, to prevent
!      ! possible division by zero errors
!
!      if (microp_uniform) then
!
!         cldm(i,k)=mincld
!         cldmw(i,k)=mincld
!         icldm(i,k)=mincld
!         lcldm(i,k)=mincld
!
!         if (qc(i,k).ge.qsmall) then
!            lcldm(i,k)=1._r8           
!            cldm(i,k)=1._r8
!            cldmw(i,k)=1._r8
!         end if
!
!         if (qi(i,k).ge.qsmall) then             
!            cldm(i,k)=1._r8
!            icldm(i,k)=1._r8
!         end if
!
!      end if               ! sub-columns
!
!      ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
!
!      nfice(i,k)=0._r8
!      dumfice=qc(i,k)+qi(i,k)
!      if (dumfice.gt.qsmall .and. qi(i,k).gt.qsmall) then
!         nfice(i,k)=qi(i,k)/dumfice
!      endif
!
!      if (do_cldice .and. (t(i,k).lt.tmelt - 5._r8)) then
!
!         ! if aerosols interact with ice set number of activated ice nuclei
!         dum2=naai(i,k)
!
!         dumnnuc=(dum2-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
!         dumnnuc=max(dumnnuc,0._r8)
!         ! get provisional ni and qi after nucleation in order to calculate
!         ! Bergeron process below
!         ninew=ni(i,k)+dumnnuc*deltat
!         qinew=qi(i,k)+dumnnuc*deltat*mi0
!
!         !T>268
!      else
!         ninew=ni(i,k)
!         qinew=qi(i,k)
!      end if
!
!      ! Initialize CME components
!
!      cme(i,k) = 0._r8
!      cmei(i,k)=0._r8
!
!
!      !-------------------------------------------------------------------
!      !Bergeron process
!
!      ! make sure to initialize bergeron process to zero
!      berg(i,k)=0._r8
!      prd = 0._r8
!
!      !condensation loop.
!
!      ! get in-cloud qi and ni after nucleation
!      if (icldm(i,k) .gt. 0._r8) then 
!         qiic(i,k)=qinew/icldm(i,k)
!         niic(i,k)=ninew/icldm(i,k)
!      else
!         qiic(i,k)=0._r8
!         niic(i,k)=0._r8
!      endif
!
!      !! Check
!      !qrout2(i,k) = relhum(i,k)
!      !qsout2(i,k) = qiic(i,k)
!      !nrout2(i,k) = niic(i,k)
!      !nsout2(i,k) = nfice(i,k)
!      !drout2(i,k) = qinew
!      !dsout2(i,k) = ninew
!
!
!      !if T < 0 C then bergeron.
!      if (do_cldice .and. (t(i,k).lt.273.15_r8)) then
!
!         !if ice exists
!         if (qi(i,k).gt.qsmall) then
!
!            bergtsf = 0._r8 ! bergeron time scale (fraction of timestep)
!
!            qvi = svp_to_qsat(esi(i,k), p(i,k))
!            qvl = svp_to_qsat(esl(i,k), p(i,k))
!
!            dqsidt =  xxls*qvi/(rv*t(i,k)*t(i,k))
!            abi = 1._r8+dqsidt*xxls/cpp
!
!            ! get ice size distribution parameters
!
!            if (qiic(i,k).ge.qsmall) then
!                tmp_array(1) = (cons1*ci* &
!                    niic(i,k)/qiic(i,k))
!                tmp_array(2) = 1.d0/di
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!               lami(k) = tmp1
!
!
!               !lami(k) = (cons1*ci* &
!               !     niic(i,k)/qiic(i,k))**(1._r8/di)
!               n0i(k) = niic(i,k)*lami(k)
!
!               ! check for slope
!               ! adjust vars
!               if (lami(k).lt.lammini) then
!
!                  lami(k) = lammini
!
!                tmp_array(1) = lami(k)
!                tmp_array(2) = di+1.d0
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!                  n0i(k) = tmp1*qiic(i,k)/(ci*cons1)
!
!                  !n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!               else if (lami(k).gt.lammaxi) then
!                  lami(k) = lammaxi
!
!                tmp_array(1) = lami(k)
!                tmp_array(2) = di+1.d0
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!                  n0i(k) = tmp1*qiic(i,k)/(ci*cons1)
!
!                  !n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!               end if
!
!               epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))
!
!               !if liquid exists  
!               if (qc(i,k).gt. qsmall) then 
!
!                  !begin bergeron process
!                  !     do bergeron (vapor deposition with RHw=1)
!                  !     code to find berg (a rate) goes here
!
!                  ! calculate Bergeron process
!
!                  prd = epsi*(qvl-qvi)/abi
!
!               else
!                  prd = 0._r8
!               end if
!
!               ! multiply by cloud fraction
!
!               prd = prd*min(icldm(i,k),lcldm(i,k))
!
!               !     transfer of existing cloud liquid to ice
!
!               berg(i,k)=max(0._r8,prd)
!
!      !qrout2(i,k) = qrout2(i,k)+epsi
!      !qsout2(i,k) = qsout2(i,k)+qvl-qvi
!      !nrout2(i,k) = nrout2(i,k)+abi
!      !nsout2(i,k) = nsout2(i,k)+prd
!      !drout2(i,k) = drout2(i,k)+berg(i,k)
!      !dsout2(i,k) = dsout2(i,k)+lami(k)
!
!
!            end if  !end liquid exists bergeron
!
!            if (berg(i,k).gt.0._r8) then
!               bergtsf=max(0._r8,(qc(i,k)/berg(i,k))/deltat) 
!
!               if(bergtsf.lt.1._r8) berg(i,k) = max(0._r8,qc(i,k)/deltat)
!
!            endif
!
!      !qrout2(i,k) = qrout2(i,k)+bergtsf
!      !qsout2(i,k) = qsout2(i,k)+qiic(i,k)
!      !nrout2(i,k) = nrout2(i,k)+qc(i,k)
!      !nsout2(i,k) = nsout2(i,k)+relhum(i,k)
!      !drout2(i,k) = drout2(i,k)+icldm(i,k)
!      !dsout2(i,k) = dsout2(i,k)+lcldm(i,k)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            if (bergtsf.lt.1._r8.or.icldm(i,k).gt.lcldm(i,k)) then
!
!               if (qiic(i,k).ge.qsmall) then
!
!                  ! first case is for case when liquid water is present, but is completely depleted 
!                  ! in time step, i.e., bergrsf > 0 but < 1
!
!                  if (qc(i,k).ge.qsmall) then
!                     rhin  = (1.0_r8 + relhum(i,k)) / 2._r8
!                     if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then
!                        prd = epsi*(rhin*qvl-qvi)/abi
!
!                        ! multiply by cloud fraction assuming liquid/ice maximum overlap
!                        prd = prd*min(icldm(i,k),lcldm(i,k))
!
!                        ! add to cmei
!                        cmei(i,k) = cmei(i,k) + (prd * (1._r8- bergtsf))
!
!
!                     end if ! rhin 
!                  end if ! qc > qsmall
!
!                  ! second case is for pure ice cloud, either no liquid, or icldm > lcldm
!
!                  if (qc(i,k).lt.qsmall.or.icldm(i,k).gt.lcldm(i,k)) then
!
!                     ! note: for case of no liquid, need to set liquid cloud fraction to zero
!                     ! store liquid cloud fraction in 'dum'
!
!                     if (qc(i,k).lt.qsmall) then 
!                        dum=0._r8 
!                     else
!                        dum=lcldm(i,k)
!                     end if
!
!                     ! set RH to grid-mean value for pure ice cloud
!                     rhin = relhum(i,k)
!
!                     if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then
!
!                        prd = epsi*(rhin*qvl-qvi)/abi
!
!                        ! multiply by relevant cloud fraction for pure ice cloud
!                        ! assuming maximum overlap of liquid/ice
!                        prd = prd*max((icldm(i,k)-dum),0._r8)
!                        cmei(i,k) = cmei(i,k) + prd
!
!                     end if ! rhin
!                  end if ! qc or icldm > lcldm
!               end if ! qiic
!            end if ! bergtsf or icldm > lcldm
!
!            !     if deposition, it should not reduce grid mean rhi below 1.0
!            if(cmei(i,k) > 0.0_r8 .and. (relhum(i,k)*esl(i,k)/esi(i,k)) > 1._r8 ) &
!                 cmei(i,k)=min(cmei(i,k),(q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat)
!
!         end if            !end ice exists loop
!         !this ends temperature < 0. loop
!
!         !-------------------------------------------------------------------
!      end if  ! 
!      !..............................................................
!
!      ! evaporation should not exceed available water
!
!      if ((-berg(i,k)).lt.-qc(i,k)/deltat) berg(i,k) = max(qc(i,k)/deltat,0._r8)
!
!      !sublimation process...
!      if (do_cldice .and. ((relhum(i,k)*esl(i,k)/esi(i,k)).lt.1._r8 .and. qiic(i,k).ge.qsmall )) then
!
!         qvi = svp_to_qsat(esi(i,k), p(i,k))
!         qvl = svp_to_qsat(esl(i,k), p(i,k))
!         dqsidt =  xxls*qvi/(rv*t(i,k)*t(i,k))
!         abi = 1._r8+dqsidt*xxls/cpp
!
!         ! get ice size distribution parameters
!
!                tmp_array(1) = (cons1*ci* &
!              niic(i,k)/qiic(i,k))
!                tmp_array(2) = 1.d0/di
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!         lami(k) = tmp1
!
!         !lami(k) = (cons1*ci* &
!         !     niic(i,k)/qiic(i,k))**(1._r8/di)
!         n0i(k) = niic(i,k)*lami(k)
!
!         ! check for slope
!         ! adjust vars
!         if (lami(k).lt.lammini) then
!
!            lami(k) = lammini
!
!                tmp_array(1) = lami(k)
!                tmp_array(2) = di+1.d0
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!                  n0i(k) = tmp1*qiic(i,k)/(ci*cons1)
!
!            !n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!         else if (lami(k).gt.lammaxi) then
!            lami(k) = lammaxi
!
!                tmp_array(1) = lami(k)
!                tmp_array(2) = di+1.d0
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!                  n0i(k) = tmp1*qiic(i,k)/(ci*cons1)
!
!            !n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!         end if
!
!         epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))
!
!         ! modify for ice fraction below
!         prd = epsi*(relhum(i,k)*qvl-qvi)/abi * icldm(i,k)
!         cmei(i,k)=min(prd,0._r8)
!
!      endif
!
!      ! sublimation should not exceed available ice
!      if (cmei(i,k).lt.-qi(i,k)/deltat) cmei(i,k)=-qi(i,k)/deltat
!
!      ! sublimation should not increase grid mean rhi above 1.0 
!      if(cmei(i,k) < 0.0_r8 .and. (relhum(i,k)*esl(i,k)/esi(i,k)) < 1._r8 ) &
!           cmei(i,k)=min(0._r8,max(cmei(i,k),(q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat))
!
!       !limit cmei due for roundoff error
!
!      cmei(i,k)=cmei(i,k)*omsm
!
!      ! conditional for ice nucleation 
!      if (do_cldice .and. (t(i,k).lt.(tmelt - 5._r8))) then 
!
!         ! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
!         ! ice nucleation rate (dum2) has already been calculated and read in (naai)
!
!         dum2i(i,k)=naai(i,k)
!      else
!         dum2i(i,k)=0._r8
!      end if
!
!   end do ! i loop
!end do ! k loop
!
!!do i=1,ncol
!!do k=top_lev,pver
!!      qrout2(i,k) = qiic(i,k)
!!      qsout2(i,k) = niic(i,k)
!!      nrout2(i,k) = nfice(i,k)
!!      nsout2(i,k) = berg(i,k)
!!      drout2(i,k) = cmei(i,k)
!!      dsout2(i,k) = dum2i(i,k)
!!end do
!!end do
!
!
!!! initialize sub-step precip flux variables
!do i=1,ncol
!   !! flux is zero at top interface, so these should stay as 0.
!   rflx1(i,1)=0._r8
!   sflx1(i,1)=0._r8
!   do k=top_lev,pver
!
!      ! initialize normal and sub-step precip flux variables
!      rflx1(i,k+1)=0._r8
!      sflx1(i,k+1)=0._r8
!   end do ! i loop
!end do ! k loop
!!! initialize final precip flux variables.
!do i=1,ncol
!   !! flux is zero at top interface, so these should stay as 0.
!   rflx(i,1)=0._r8
!   sflx(i,1)=0._r8
!   do k=top_lev,pver
!      ! initialize normal and sub-step precip flux variables
!      rflx(i,k+1)=0._r8
!      sflx(i,k+1)=0._r8
!   end do ! i loop
!end do ! k loop
!
!do i=1,ncol
!   ltrue(i)=0
!   do k=top_lev,pver
!      ! skip microphysical calculations if no cloud water
!
!      if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall.or.cmei(i,k).ge.qsmall) ltrue(i)=1
!   end do
!end do
!
!!do i=1,pcols
!!do k=1,pver
!!      qrout2(i,k) = qc(i,k)
!!      qsout2(i,k) = qi(i,k)
!!      nrout2(i,k) = cmei(i,k)
!!      nsout2(i,k) = qsmall
!!      drout2(i,k) = ltrue(i)
!!      dsout2(i,k) = deltat
!!end do
!!end do
!
!! assign number of sub-steps to iter
!! use 2 sub-steps, following tests described in MG2008
!iter = 2
!
!! get sub-step time step
!deltat=deltat/real(iter)
!
!
!! since activation/nucleation processes are fast, need to take into account
!! factor mtime = mixing timescale in cloud / model time step
!! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
!! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk
!
!!        note: mtime for bulk aerosols was set to: mtime=deltat/1200._r8
!
!mtime=1._r8
!rate1ord_cw2pr_st(:,:)=0._r8 ! rce 2010/05/01
!
!!!!! skip calculations if no cloud water
!do i=1,ncol
!
!!pra(:) = 0._r8
!!pre(:) = 0._r8
!!pracs(:) = 0._r8
!!mnuccr(:) = 0._r8
!
!   if (ltrue(i).eq.0) then
!      tlat(i,1:pver)=0._r8
!      qvlat(i,1:pver)=0._r8
!      qctend(i,1:pver)=0._r8
!      qitend(i,1:pver)=0._r8
!      qnitend(i,1:pver)=0._r8
!      qrtend(i,1:pver)=0._r8
!      nctend(i,1:pver)=0._r8
!      nitend(i,1:pver)=0._r8
!      nrtend(i,1:pver)=0._r8
!      nstend(i,1:pver)=0._r8
!      prect(i)=0._r8
!      preci(i)=0._r8
!      qniic(i,1:pver)=0._r8
!      qric(i,1:pver)=0._r8
!      nsic(i,1:pver)=0._r8
!      nric(i,1:pver)=0._r8
!      rainrt(i,1:pver)=0._r8
!      goto 300
!   end if
!
!   qcsinksum_rate1ord(1:pver)=0._r8 
!   qcsum_rate1ord(1:pver)=0._r8 
!
!            !do it=1,iter
!            !do k=1,pver
!            !      qrout2(i,k) = qrout2(i,k) + rflx1(i,k)
!            !      qsout2(i,k) = qsout2(i,k) + sflx1(i,k)
!            !      nrout2(i,k) = nrout2(i,k) + rflx(i,k)
!            !      nsout2(i,k) = nsout2(i,k) + sflx(i,k)
!            !      drout2(i,k) = drout2(i,k) + qc(i,k)
!            !      dsout2(i,k) = dsout2(i,k) + qi(i,k)
!            !end do
!            !end do
!
!!!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !.....................................................................................................
!   do it=1,iter
!
!      ! initialize sub-step microphysical tendencies
!
!      tlat(i,1:pver)=0._r8
!      qvlat(i,1:pver)=0._r8
!      qctend(i,1:pver)=0._r8
!      qitend(i,1:pver)=0._r8
!      qnitend(i,1:pver)=0._r8
!      qrtend(i,1:pver)=0._r8
!      nctend(i,1:pver)=0._r8
!      nitend(i,1:pver)=0._r8
!      nrtend(i,1:pver)=0._r8
!      nstend(i,1:pver)=0._r8
!
!      ! initialize diagnostic precipitation to zero
!
!      qniic(i,1:pver)=0._r8
!      qric(i,1:pver)=0._r8
!      nsic(i,1:pver)=0._r8
!      nric(i,1:pver)=0._r8
!
!      rainrt(i,1:pver)=0._r8
!
!
!      ! begin new i,k loop, calculate new cldmax after adjustment to cldm above
!
!      ! initialize vertically-integrated rain and snow tendencies
!
!      qrtot = 0._r8
!      nrtot = 0._r8
!      qstot = 0._r8
!      nstot = 0._r8
!
!      ! initialize precip at surface
!
!      prect(i)=0._r8
!      preci(i)=0._r8
!
!      do k=top_lev,pver
!      
!         qcvar=relvar(i,k)
!         cons2=gamma(qcvar+2.47_r8)
!         cons3=gamma(qcvar)
!         cons9=gamma(qcvar+2._r8)
!         cons10=gamma(qcvar+1._r8)
!         cons12=gamma(qcvar+1.15_r8) 
!         cons15=gamma(qcvar+bc/3._r8)
!
!         tmp_array(1) = qcvar
!         tmp_array(2) = 2.47d0
!         call athread_spawn(slave_pow_parallel, tmp_array)
!         call athread_join()
!         cons18 = tmp_array(1)
!
!         cons19=qcvar*qcvar
!
!         tmp_array(1) = qcvar
!         tmp_array(2) = 1.15d0
!         call athread_spawn(slave_pow_parallel, tmp_array)
!         call athread_join()
!         cons20 = tmp_array(1)
!
!
!         !cons18=qcvar**2.47_r8
!         !cons19=qcvar**2
!         !cons20=qcvar**1.15_r8
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + qcvar
!         !    qsout2(i,k) = qsout2(i,k) + cons2
!         !    nrout2(i,k) = nrout2(i,k) + cons15
!         !    nsout2(i,k) = nsout2(i,k) + cons18
!         !    drout2(i,k) = drout2(i,k) + qc(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + qi(i,k)
!         !end if
!
!         ! set cwml and cwmi to current qc and qi
!
!         cwml(i,k)=qc(i,k)
!         cwmi(i,k)=qi(i,k)
!
!         ! initialize precip fallspeeds to zero
!
!         ums(k)=0._r8 
!         uns(k)=0._r8 
!         umr(k)=0._r8 
!         unr(k)=0._r8
!
!         ! calculate precip fraction based on maximum overlap assumption
!
!         ! for sub-columns cldm has already been set to 1 if cloud
!         ! water or ice is present, so cldmax will be correctly set below
!         ! and nothing extra needs to be done here
!
!         if (k.eq.top_lev) then
!            cldmax(i,k)=cldm(i,k)
!         else
!            ! if rain or snow mix ratio is smaller than
!            ! threshold, then set cldmax to cloud fraction at current level
!
!            if (do_clubb_sgs) then
!               if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall) then
!                  cldmax(i,k)=cldm(i,k)
!               else
!                  cldmax(i,k)=cldmax(i,k-1)
!               end if
!            else
!
!               if (qric(i,k-1).ge.qsmall.or.qniic(i,k-1).ge.qsmall) then
!                  cldmax(i,k)=max(cldmax(i,k-1),cldm(i,k))
!               else
!                  cldmax(i,k)=cldm(i,k)
!               end if
!            endif
!         end if
!
!         !if ( it .eq. 1 .and. k > 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + cldmax(i,k-1) 
!         !    qsout2(i,k) = qsout2(i,k) + cldmax(i,k)
!         !    nrout2(i,k) = nrout2(i,k) + qniic(i,k-1)
!         !    nsout2(i,k) = nsout2(i,k) + qi(i,k-1)
!         !    drout2(i,k) = drout2(i,k) + cldm(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + qric(i,k-1)
!         !endif
!
!         ! decrease in number concentration due to sublimation/evap
!         ! divide by cloud fraction to get in-cloud decrease
!         ! don't reduce Nc due to bergeron process
!
!         if (cmei(i,k) < 0._r8 .and. qi(i,k) > qsmall .and. cldm(i,k) > mincld) then
!            nsubi(k)=cmei(i,k)/qi(i,k)*ni(i,k)/cldm(i,k)
!         else
!            nsubi(k)=0._r8
!         end if
!         nsubc(k)=0._r8
!
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
!
!         if (do_cldice .and. dum2i(i,k).gt.0._r8.and.t(i,k).lt.(tmelt - 5._r8).and. &
!              relhum(i,k)*esl(i,k)/esi(i,k).gt. rhmini+0.05_r8) then
!
!            !if NCAI > 0. then set numice = ncai (as before)
!            !note: this is gridbox averaged
!
!            nnuccd(k)=(dum2i(i,k)-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
!            nnuccd(k)=max(nnuccd(k),0._r8)
!            nimax = dum2i(i,k)*icldm(i,k)
!
!            !Calc mass of new particles using new crystal mass...
!            !also this will be multiplied by mtime as nnuccd is...
!
!            mnuccd(k) = nnuccd(k) * mi0
!
!            !  add mnuccd to cmei....
!            cmei(i,k)= cmei(i,k) + mnuccd(k) * mtime
!
!            !  limit cmei
!
!            qvi = svp_to_qsat(esi(i,k), p(i,k))
!            !dqsidt =  xxls*qvi/(rv*t(i,k)**2)
!            dqsidt =  xxls*qvi/(rv*t(i,k)*t(i,k))
!            abi = 1._r8+dqsidt*xxls/cpp
!            cmei(i,k)=min(cmei(i,k),(q(i,k)-qvi)/abi/deltat)
!
!            ! limit for roundoff error
!            cmei(i,k)=cmei(i,k)*omsm
!
!         else
!            nnuccd(k)=0._r8
!            nimax = 0._r8
!            mnuccd(k) = 0._r8
!         end if
!
!         !c............................................................................
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
!         ! for microphysical process calculations
!         ! units are kg/kg for mixing ratio, 1/kg for number conc
!
!         ! limit in-cloud values to 0.005 kg/kg
!
!         qcic(i,k)=min(cwml(i,k)/lcldm(i,k),5.e-3_r8)
!         qiic(i,k)=min(cwmi(i,k)/icldm(i,k),5.e-3_r8)
!         ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)
!         niic(i,k)=max(ni(i,k)/icldm(i,k),0._r8)
!
!         if (qc(i,k) - berg(i,k)*deltat.lt.qsmall) then
!            qcic(i,k)=0._r8
!            ncic(i,k)=0._r8
!            if (qc(i,k)-berg(i,k)*deltat.lt.0._r8) then
!               berg(i,k)=qc(i,k)/deltat*omsm
!            end if
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+qcic(i,k)
!         !    qsout2(i,k) = qsout2(i,k)+qiic(i,k)
!         !    nrout2(i,k) = nrout2(i,k)+ncic(i,k)
!         !    nsout2(i,k) = nsout2(i,k)+niic(i,k)
!         !    drout2(i,k) = drout2(i,k)+qc(i,k)
!         !    dsout2(i,k) = dsout2(i,k)+berg(i,k)
!         !endif
!
!         if (do_cldice .and. qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.qsmall) then
!            qiic(i,k)=0._r8
!            niic(i,k)=0._r8
!            if (qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.0._r8) then
!               cmei(i,k)=(-qi(i,k)/deltat-berg(i,k))*omsm
!            end if
!         end if
!
!         ! add to cme output
!
!         cmeout(i,k) = cmeout(i,k)+cmei(i,k)
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + qcic(i,k)
!         !    qsout2(i,k) = qsout2(i,k) + ncic(i,k)
!         !    nrout2(i,k) = nrout2(i,k) + berg(i,k)
!         !    nsout2(i,k) = nsout2(i,k) + niic(i,k)
!         !    drout2(i,k) = drout2(i,k) + qiic(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + cmeout(i,k)
!         !end if
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! droplet activation
!         ! calculate potential for droplet activation if cloud water is present
!         ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
!         ! number tendency (npccnin) is read in from companion routine
!
!         ! assume aerosols already activated are equal to number of existing droplets for simplicity
!         ! multiply by cloud fraction to obtain grid-average tendency
!
!         if (qcic(i,k).ge.qsmall) then   
!            npccn(k) = max(0._r8,npccnin(i,k))  
!            dum2l(i,k)=(nc(i,k)+npccn(k)*deltat)/lcldm(i,k)
!            dum2l(i,k)=max(dum2l(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3  
!            ncmax = dum2l(i,k)*lcldm(i,k)
!         else
!            npccn(k)=0._r8
!            dum2l(i,k)=0._r8
!            ncmax = 0._r8
!         end if
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! get size distribution parameters based on in-cloud cloud water/ice 
!         ! these calculations also ensure consistency between number and mixing ratio
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         !......................................................................
!         ! cloud ice
!
!         if (qiic(i,k).ge.qsmall) then
!
!            ! add upper limit to in-cloud number concentration to prevent numerical error
!            niic(i,k)=min(niic(i,k),qiic(i,k)*1.e20_r8)
!
!            tmp_array(1) = (cons1*ci*niic(i,k)/qiic(i,k))
!            tmp_array(2) = 1._r8/di
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp1 = tmp_array(1)
!            lami(k) = tmp1
!
!            !lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**(1._r8/di)
!            n0i(k) = niic(i,k)*lami(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lami(k).lt.lammini) then
!
!               lami(k) = lammini
!
!               tmp_array(1) = lami(k)
!               tmp_array(2) = di+1._r8
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               n0i(k) = tmp1*qiic(i,k)/(ci*cons1)
!
!               !n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!               niic(i,k) = n0i(k)/lami(k)
!            else if (lami(k).gt.lammaxi) then
!               lami(k) = lammaxi
!
!               tmp_array(1) = lami(k)
!               tmp_array(2) = di+1._r8
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               n0i(k) = tmp1*qiic(i,k)/(ci*cons1)
!
!               !n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
!               niic(i,k) = n0i(k)/lami(k)
!            end if
!
!         else
!            lami(k) = 0._r8
!            n0i(k) = 0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+npccn(k)
!         !    qsout2(i,k) = qsout2(i,k)+ncmax
!         !    nrout2(i,k) = nrout2(i,k)+dum2l(i,k)
!         !    nsout2(i,k) = nsout2(i,k)+lami(k)
!         !    drout2(i,k) = drout2(i,k)+qcic(i,k)
!         !    dsout2(i,k) = dsout2(i,k)+qsmall
!         !endif
!         pgam(k) = 0.d0
!         lammin = 0.d0
!         lammax = 0.d0
!         if (qcic(i,k).ge.qsmall) then
!
!            ! add upper limit to in-cloud number concentration to prevent numerical error
!            ncic(i,k)=min(ncic(i,k),qcic(i,k)*1.e20_r8)
!
!            ncic(i,k)=max(ncic(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm  
!
!            ! get pgam from fit to observations of martin et al. 1994
!
!            pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
!            pgam(k)=1._r8/(pgam(k)*pgam(k))-1._r8
!            !pgam(k)=1._r8/(pgam(k)**2)-1._r8
!            pgam(k)=max(pgam(k),2._r8)
!            pgam(k)=min(pgam(k),15._r8)
!
!            ! calculate lamc
!
!            !lamc(k) = (pi/6._r8*rhow*ncic(i,k)*gamma(pgam(k)+4._r8)/ &
!            !     (qcic(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!
!
!            ! lammin, 50 micron diameter max mean size
!
!            lammin = (pgam(k)+1._r8)/50.e-6_r8
!            lammax = (pgam(k)+1._r8)/2.e-6_r8
!            
!            tmp1 = pgam(k) + 4._r8
!            call athread_spawn(slave_gamma_parallel, tmp1)
!            call athread_join()
!            tmp2 = pgam(k) + 1._r8
!            call athread_spawn(slave_gamma_parallel, tmp2)
!            call athread_join()
!            tmp_array(1) = (pi/6._r8*rhow*ncic(i,k)*tmp1/ &
!                 (qcic(i,k)*tmp2))
!            tmp_array(2) = 1.d0/3.d0
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!
!            lamc(k) = tmp_array(1)
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + pgam(k)
!         !    qsout2(i,k) = qsout2(i,k) + lammin
!         !    nrout2(i,k) = nrout2(i,k) + lammax
!         !    !nsout2(i,k) = nsout2(i,k) + qcic(i,k)
!         !    !drout2(i,k) = drout2(i,k) + tmp1
!         !    !dsout2(i,k) = dsout2(i,k) + tmp2
!         !endif
!
!            if (lamc(k).lt.lammin) then
!               lamc(k) = lammin
!               !ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
!               !     gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
!               ncic(i,k) = 6._r8*lamc(k)*lamc(k)*lamc(k)*qcic(i,k)* &
!                    tmp2/(pi*rhow*tmp1)
!            else if (lamc(k).gt.lammax) then
!               lamc(k) = lammax
!               ncic(i,k) = 6._r8*lamc(k)*lamc(k)*lamc(k)*qcic(i,k)* &
!                    tmp2/(pi*rhow*tmp1)
!               !ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
!               !     gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
!            end if
!
!            ! parameter to calculate droplet freezing
!
!            !cdist1(k) = ncic(i,k)/gamma(pgam(k)+1._r8) 
!            cdist1(k) = ncic(i,k)/tmp2
!
!         else
!            lamc(k) = 0._r8
!            cdist1(k) = 0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + ncic(i,k) 
!         !    qsout2(i,k) = qsout2(i,k) + pgam(k)
!         !    nrout2(i,k) = nrout2(i,k) + lamc(k)
!         !    nsout2(i,k) = nsout2(i,k) + pgam(k)
!         !    drout2(i,k) = drout2(i,k) + lammin
!         !    dsout2(i,k) = dsout2(i,k) + lammax
!         !endif
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! begin micropysical process calculations 
!         !.................................................................
!         ! autoconversion of cloud liquid water to rain
!         ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
!         ! minimum qc of 1 x 10^-8 prevents floating point error
!
!         if (qcic(i,k).ge.1.e-8_r8) then
!
!            ! nprc is increase in rain number conc due to autoconversion
!            ! nprc1 is decrease in cloud droplet conc due to autoconversion
!
!            ! assume exponential sub-grid distribution of qc, resulting in additional
!            ! factor related to qcvar below
!
!            ! hm switch for sub-columns, don't include sub-grid qc
!            if (microp_uniform) then
!
!               tmp_array(1) = qcic(i,k)
!               tmp_array(2) = 2.47d0
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!
!               tmp_array(1) = (ncic(i,k)/1.e6_r8*rho(i,k))
!               tmp_array(2) = -1.79d0
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp2 = tmp_array(1)
!               prc(k) = 1350._r8*tmp1* &
!                   tmp2 
!
!
!               !prc(k) = 1350._r8*qcic(i,k)**2.47_r8* &
!               !     (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
!               nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)*(25.e-6_r8)*(25.e-6_r8))
!               !nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
!               nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
!
!            else
!
!               tmp_array(1) = qcic(i,k)
!               tmp_array(2) = 2.47d0
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!
!               tmp_array(1) = (ncic(i,k)/1.e6_r8*rho(i,k))
!               tmp_array(2) = -1.79d0
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp2 = tmp_array(1)
!               prc(k) = cons2/(cons3*cons18)*1350._r8*tmp1* &
!                   tmp2 
!
!               !prc(k) = cons2/(cons3*cons18)*1350._r8*qcic(i,k)**2.47_r8* &
!               !     (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
!               nprc(k) = prc(k)/cons22
!               nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
!
!            end if               ! sub-column switch
!
!         else
!            prc(k)=0._r8
!            nprc(k)=0._r8
!            nprc1(k)=0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + qcic(i,k) 
!         !    qsout2(i,k) = qsout2(i,k) + ncic(i,k)
!         !    nrout2(i,k) = nrout2(i,k) + rho(i,k)
!         !    nsout2(i,k) = nsout2(i,k) + prc(k)
!         !    drout2(i,k) = drout2(i,k) + qric(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + nric(i,k)
!         !endif
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + qcic(i,k) 
!         !    qsout2(i,k) = qsout2(i,k) + nprc1(k)
!         !    nrout2(i,k) = nrout2(i,k) + prc(k)
!         !    nsout2(i,k) = nsout2(i,k) + rho(i,k)
!         !    drout2(i,k) = drout2(i,k) + qric(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + nric(i,k)
!         !endif
!         ! add autoconversion to precip from above to get provisional rain mixing ratio
!         ! and number concentration (qric and nric)
!
!         ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)
!
!         dum=0.45_r8
!         dum1=0.45_r8
!
!         !if ( it .eq. 1 .and. k.eq.pver/2) then
!         !    qrout2(i,1) = qrout2(i,1) + umr(pver/2) 
!         !    qsout2(i,1) = qsout2(i,1) + unr(pver/2)
!         !    nrout2(i,1) = nrout2(i,1) + qniic(i,pver/2)
!         !    nsout2(i,1) = nsout2(i,1) + rho(i,pver/2)
!         !    drout2(i,1) = drout2(i,1) + nric(i,pver/2)
!         !    dsout2(i,1) = dsout2(i,1) + qric(i,pver/2)
!         !    qrout2(i,16) = qrout2(i,16) + umr(pver/2) 
!         !    qsout2(i,16) = qsout2(i,16) + unr(pver/2)
!         !    nrout2(i,16) = nrout2(i,16) + qniic(i,pver/2)
!         !    nsout2(i,16) = nsout2(i,16) + rho(i,pver/2)
!         !    drout2(i,16) = drout2(i,16) + nric(i,pver/2)
!         !    dsout2(i,16) = dsout2(i,16) + qric(i,pver/2)
!         !endif
!
!         !if ( it .eq. 1 .and. k <=8 .and. k > 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + unr(k-1) 
!         !    qsout2(i,k) = qsout2(i,k) + nric(k-1)
!         !    nrout2(i,k) = nrout2(i,k) + nsubr(k-1)
!         !    nsout2(i,k) = nsout2(i,k) + npracs(k-1)
!         !    drout2(i,k) = drout2(i,k) + nnuccr(k-1)
!         !    dsout2(i,k) = dsout2(i,k) + nragg(k-1)
!         !endif
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + umr(k) 
!         !    qsout2(i,k) = qsout2(i,k) + unr(k)
!         !    nrout2(i,k) = nrout2(i,k) + qniic(i,k)
!         !    nsout2(i,k) = nsout2(i,k) + pra(k)
!         !    drout2(i,k) = drout2(i,k) + pre(k)
!         !    dsout2(i,k) = dsout2(i,k) + qric(i,k)
!         !endif
!         !if ( it .eq. 1 .and. k > 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + umr(k-1) 
!         !    qsout2(i,k) = qsout2(i,k) + qric(i,k-1)
!         !    nrout2(i,k) = nrout2(i,k) + pra(k-1)
!         !    nsout2(i,k) = nsout2(i,k) + pre(k-1)
!         !    drout2(i,k) = drout2(i,k) + pracs(k-1)
!         !    dsout2(i,k) = dsout2(i,k) + mnuccr(k-1)
!         !endif
!         if (k.eq.top_lev) then
!             qric(i,k)=prc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
!             nric(i,k)=nprc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
!         else 
!             if (qric(i,k-1).ge.qsmall) then
!                 dum=umr(k-1)
!                 dum1=unr(k-1)
!             end if
!
!             ! no autoconversion of rain number if rain/snow falling from above
!             ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
!             ! by the existing rain/snow particles from above
!
!             if (qric(i,k-1).ge.1.e-9_r8.or.qniic(i,k-1).ge.1.e-9_r8) then
!                 nprc(k)=0._r8
!             end if
!
!             qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
!             (rho(i,k)*dz(i,k)*((pra(k-1)+prc(k))*lcldm(i,k)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax(i,k))))&
!             /(dum*rho(i,k)*cldmax(i,k))
!             nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
!             (rho(i,k)*dz(i,k)*(nprc(k)*lcldm(i,k)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax(i,k))))&
!             /(dum1*rho(i,k)*cldmax(i,k))
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    !qrout2(i,k) = qrout2(i,k) + umr(k-1)
!         !    !qsout2(i,k) = qsout2(i,k) + rho(i,k-1)
!         !    !nrout2(i,k) = nrout2(i,k) + pra(k-1)
!         !    !nsout2(i,k) = nsout2(i,k) + pre(k-1)
!         !    !drout2(i,k) = drout2(i,k) + pracs(k-1)
!         !    !dsout2(i,k) = dsout2(i,k) + mnuccr(k-1)
!
!         !    qrout2(i,k) = qrout2(i,k) + cldmax(i,k-1)
!         !    qsout2(i,k) = qsout2(i,k) + rho(i,k)
!         !    nrout2(i,k) = nrout2(i,k) + dz(i,k)
!         !    nsout2(i,k) = nsout2(i,k) + prc(k)
!         !    drout2(i,k) = drout2(i,k) + lcldm(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + cldmax(i,k)
!         !endif
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + qric(i,k)
!         !endif
!         !if ( it .eq. 1 .and. k > 1) then
!         !    qrout2(i,k) = qrout2(i,k) + pra(k-1)
!         !    qsout2(i,k) = qsout2(i,k) + pre(k-1)
!         !    nrout2(i,k) = nrout2(i,k) + pracs(k-1)
!         !    nsout2(i,k) = nsout2(i,k) + mnuccr(k-1)
!         !    drout2(i,k) = drout2(i,k) + qric(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + nric(i,k)
!         !endif
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + niic(i,k)
!         !    qsout2(i,k) = qsout2(i,k) + berg(i,k)
!         !    nrout2(i,k) = nrout2(i,k) + qcic(i,k)
!         !    nsout2(i,k) = nsout2(i,k) + ncic(i,k)
!         !    drout2(i,k) = drout2(i,k) + qric(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + nric(i,k)
!         !endif
!         !.......................................................................
!         ! Autoconversion of cloud ice to snow
!         ! similar to Ferrier (1994)
!
!         if (do_cldice) then
!            if (t(i,k).le.273.15_r8.and.qiic(i,k).ge.qsmall) then
!
!               ! note: assumes autoconversion timescale of 180 sec
!               
!               tmp1 = -lami(k)*dcs
!               call athread_spawn(slave_exp_parallel, tmp1)
!               call athread_join()
!
!               nprci(k) = n0i(k)/(lami(k)*180._r8)*tmp1
!
!               prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
!                    (cons23/lami(k)+3._r8*cons24/(lami(k)*lami(k))+ &
!                    6._r8*dcs/(lami(k)*lami(k)*lami(k))+6._r8/(lami(k)*lami(k)*lami(k)*lami(k)))*tmp1
!
!               !nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)
!
!               !!prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
!               !!     (cons23/lami(k)+3._r8*cons24/lami(k)**2+ &
!               !!     6._r8*dcs/lami(k)**3+6._r8/lami(k)**4)*exp(-lami(k)*dcs)
!               !prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
!               !     (cons23/lami(k)+3._r8*cons24/(lami(k)*lami(k))+ &
!               !     6._r8*dcs/(lami(k)*lami(k)*lami(k))+6._r8/(lami(k)*lami(k)*lami(k)*lami(k)))*exp(-lami(k)*dcs)
!            else
!               prci(k)=0._r8
!               nprci(k)=0._r8
!            end if
!         else
!            ! Add in the particles that we have already converted to snow, and
!            ! don't do any further autoconversion of ice.
!            prci(k)  = tnd_qsnow(i, k) / cldm(i,k)
!            nprci(k) = tnd_nsnow(i, k) / cldm(i,k)
!         end if
!
!         ! add autoconversion to flux from level above to get provisional snow mixing ratio
!         ! and number concentration (qniic and nsic)
!
!         dum=(asn(i,k)*cons25)
!         dum1=(asn(i,k)*cons25)
!
!         !if ( it .eq. 1  .and. k > 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + uns(k-1)
!         !    qsout2(i,k) = qsout2(i,k) + nprci(k)
!         !    nrout2(i,k) = nrout2(i,k) + nsubs(k-1)
!         !    nsout2(i,k) = nsout2(i,k) + nsagg(k-1)
!         !    drout2(i,k) = drout2(i,k) + nnuccr(k-1)
!         !    dsout2(i,k) = dsout2(i,k) + rho(i,k-1)
!         !endif
!
!         if (k.eq.top_lev) then
!            qniic(i,k)=prci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
!            nsic(i,k)=nprci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
!         else
!            if (qniic(i,k-1).ge.qsmall) then
!               dum=ums(k-1)
!               dum1=uns(k-1)
!            end if
!
!            qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
!                 (rho(i,k)*dz(i,k)*((prci(k)+prai(k-1)+psacws(k-1)+bergs(k-1))*icldm(i,k)+(prds(k-1)+ &
!                 pracs(k-1)+mnuccr(k-1))*cldmax(i,k))))&
!                 /(dum*rho(i,k)*cldmax(i,k))
!
!            nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
!                 (rho(i,k)*dz(i,k)*(nprci(k)*icldm(i,k)+(nsubs(k-1)+nsagg(k-1)+nnuccr(k-1))*cldmax(i,k))))&
!                 /(dum1*rho(i,k)*cldmax(i,k))
!
!         end if
!
!         !!if ( it .eq. 1 .and. k .eq. 1 ) then
!         !if ( it .eq. 1  ) then
!         !    qrout2(i,k) = qrout2(i,k) + dum
!         !    qsout2(i,k) = qsout2(i,k) + dum1
!         !    nrout2(i,k) = nrout2(i,k) + prci(k)
!         !    nsout2(i,k) = nsout2(i,k) + nprci(k)
!         !    drout2(i,k) = drout2(i,k) + qniic(i,k)
!         !    dsout2(i,k) = dsout2(i,k) + nsic(i,k)
!         !endif
!
!         ! if precip mix ratio is zero so should number concentration
!
!         if (qniic(i,k).lt.qsmall) then
!            qniic(i,k)=0._r8
!            nsic(i,k)=0._r8
!         end if
!
!         if (qric(i,k).lt.qsmall) then
!            qric(i,k)=0._r8
!            nric(i,k)=0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qsout2(i,k) = qsout2(i,k) + qric(i,k)
!         !endif
!
!         ! make sure number concentration is a positive number to avoid 
!         ! taking root of negative later
!
!         nric(i,k)=max(nric(i,k),0._r8)
!         nsic(i,k)=max(nsic(i,k),0._r8)
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+prci(k) 
!         !    qsout2(i,k) = qsout2(i,k)+nprci(k)
!         !    nrout2(i,k) = nrout2(i,k)+qniic(i,k)
!         !    nsout2(i,k) = nsout2(i,k)+qric(i,k) 
!         !    drout2(i,k) = drout2(i,k)+nric(i,k) 
!         !    dsout2(i,k) = dsout2(i,k)+nsic(i,k) 
!         !endif
!         !.......................................................................
!         ! get size distribution parameters for precip
!         !......................................................................
!         ! rain
!
!         if (qric(i,k).ge.qsmall) then
!             tmp_array(1) = (pi*rhow*nric(i,k)/qric(i,k))
!             tmp_array(2) = 1._r8/3._r8
!             call athread_spawn(slave_pow_parallel, tmp_array)
!             call athread_join()
!             tmp1 = tmp_array(1)
!             lamr(k) = tmp1
!
!            !lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
!            n0r(k) = nric(i,k)*lamr(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lamr(k).lt.lamminr) then
!
!               lamr(k) = lamminr
!
!               !n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
!               n0r(k) = lamr(k)*lamr(k)*lamr(k)*lamr(k)*qric(i,k)/(pi*rhow)
!               nric(i,k) = n0r(k)/lamr(k)
!            else if (lamr(k).gt.lammaxr) then
!               lamr(k) = lammaxr
!               !n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
!               n0r(k) = lamr(k)*lamr(k)*lamr(k)*lamr(k)*qric(i,k)/(pi*rhow)
!               nric(i,k) = n0r(k)/lamr(k)
!            end if
!
!            ! provisional rain number and mass weighted mean fallspeed (m/s)
!
!             tmp_array(1) = lamr(k)
!             tmp_array(2) = br
!             call athread_spawn(slave_pow_parallel, tmp_array)
!             call athread_join()
!             tmp1 = tmp_array(1)
!            unr(k) = min(arn(i,k)*cons4/tmp1,9.1_r8*rhof(i,k))
!            umr(k) = min(arn(i,k)*cons5/(6._r8*tmp1),9.1_r8*rhof(i,k))
!
!            !unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
!            !umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(i,k))
!
!         else
!            lamr(k) = 0._r8
!            n0r(k) = 0._r8
!            umr(k) = 0._r8
!            unr(k) = 0._r8
!         end if
!
!         !......................................................................
!         ! snow
!
!         if (qniic(i,k).ge.qsmall) then
!             tmp_array(1) =(cons6*cs*nsic(i,k)/qniic(i,k))
!             tmp_array(2) = 1._r8/ds
!             call athread_spawn(slave_pow_parallel, tmp_array)
!             call athread_join()
!             tmp1 = tmp_array(1)
!             lams(k) = tmp1
!
!            !lams(k) = (cons6*cs*nsic(i,k)/qniic(i,k))**(1._r8/ds)
!            n0s(k) = nsic(i,k)*lams(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lams(k).lt.lammins) then
!               lams(k) = lammins
!
!               tmp_array(1) = lams(k)
!               tmp_array(2) = ds+1._r8
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               n0s(k) = tmp1*qniic(i,k)/(cs*cons6)
!
!               !n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
!               nsic(i,k) = n0s(k)/lams(k)
!
!            else if (lams(k).gt.lammaxs) then
!               lams(k) = lammaxs
!
!               tmp_array(1) = lams(k)
!               tmp_array(2) = ds+1._r8
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               n0s(k) = tmp1*qniic(i,k)/(cs*cons6)
!
!               !n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
!               nsic(i,k) = n0s(k)/lams(k)
!            end if
!
!            ! provisional snow number and mass weighted mean fallspeed (m/s)
!
!            tmp_array(1) = lams(k)
!            tmp_array(2) = bs
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp1 = tmp_array(1)
!
!            ums(k) = min(asn(i,k)*cons8/(6._r8*tmp1),1.2_r8*rhof(i,k))
!            uns(k) = min(asn(i,k)*cons7/tmp1,1.2_r8*rhof(i,k))
!
!            !ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(i,k))
!            !uns(k) = min(asn(i,k)*cons7/lams(k)**bs,1.2_r8*rhof(i,k))
!
!         else
!            lams(k) = 0._r8
!            n0s(k) = 0._r8
!            ums(k) = 0._r8
!            uns(k) = 0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+umr(k)
!         !    qsout2(i,k) = qsout2(i,k)+unr(k)
!         !    nrout2(i,k) = nrout2(i,k)+lams(k)
!         !    nsout2(i,k) = nsout2(i,k)+n0s(k)
!         !    drout2(i,k) = drout2(i,k)+ums(k)
!         !    dsout2(i,k) = dsout2(i,k)+uns(k)
!         !endif
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+cdist1(k)
!         !    qsout2(i,k) = qsout2(i,k)+pgam(k)
!         !    nrout2(i,k) = nrout2(i,k)+t(i,k)
!         !    nsout2(i,k) = nsout2(i,k)+lamc(k)
!         !    !drout2(i,k) = drout2(i,k)+use_hetfrz_classnuc
!         !    !dsout2(i,k) = dsout2(i,k)+do_cldice
!         !endif
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         ! heterogeneous freezing of cloud water
!         nnuccc(k) = 0._r8
!         mnuccc(k) = 0._r8
!
!         nnucct(k) = 0._r8
!         mnucct(k) = 0._r8
!
!         nnudep(k) = 0._r8
!         mnudep(k) = 0._r8
!         if (.not. use_hetfrz_classnuc) then
!
!            if (do_cldice .and. qcic(i,k).ge.qsmall .and. t(i,k).lt.269.15_r8) then
!
!               ! immersion freezing (Bigg, 1953)
!
!
!               ! subcolumns
!
!               if (microp_uniform) then
!
!                  !mnuccc(k) = &
!                  !   pi*pi/36._r8*rhow* &
!                  !   cdist1(k)*gamma(7._r8+pgam(k))* &
!                  !   bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
!                  !   lamc(k)**3/lamc(k)**3
!
!                  mnuccc(k) = &
!                     pi*pi/36._r8*rhow* &
!                     cdist1(k)*gamma(7._r8+pgam(k))* &
!                     bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
!                     (lamc(k)*lamc(k)*lamc(k))/(lamc(k)*lamc(k)*lamc(k))
!
!                  !nnuccc(k) = &
!                  !   pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
!                  !   *bimm* &
!                  !   (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3
!
!                  nnuccc(k) = &
!                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
!                     *bimm* &
!                     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/(lamc(k)*lamc(k)*lamc(k))
!
!               else
!
!                  tmp1 = 7.d0+pgam(k)
!                  call athread_spawn(slave_gamma_parallel, tmp1)
!                  call athread_join()
!                  tmp2 = aimm*(273.15d0-t(i,k))
!                  call athread_spawn(slave_exp_parallel, tmp2)
!                  call athread_join()
!                  tmp3 = pgam(k)+4.d0
!                  call athread_spawn(slave_gamma_parallel, tmp3)
!                  call athread_join()
!
!                  mnuccc(k) = cons9/(cons3*cons19)* &
!                     pi*pi/36._r8*rhow* &
!                     cdist1(k)*tmp1* &
!                     bimm*(tmp2-1._r8)/ &
!                     (lamc(k)*lamc(k)*lamc(k))/(lamc(k)*lamc(k)*lamc(k))
!
!                  nnuccc(k) = cons10/(cons3*qcvar)* &
!                     pi/6._r8*cdist1(k)*tmp3 &
!                     *bimm* &
!                     (tmp2-1._r8)/(lamc(k)*lamc(k)*lamc(k))
!
!                  !mnuccc(k) = 1.d0
!                  !nnuccc(k) = 1.d0
!
!                  !mnuccc(k) = cons9/(cons3*cons19)* &
!                  !   pi*pi/36._r8*rhow* &
!                  !   cdist1(k)*gamma(7._r8+pgam(k))* &
!                  !   bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
!                  !   (lamc(k)*lamc(k)*lamc(k))/(lamc(k)*lamc(k)*lamc(k))
!
!                  !nnuccc(k) = cons10/(cons3*qcvar)* &
!                  !   pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
!                  !   *bimm* &
!                  !   (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/(lamc(k)*lamc(k)*lamc(k))
!
!                  !mnuccc(k) = cons9/(cons3*cons19)* &
!                  !   pi*pi/36._r8*rhow* &
!                  !   cdist1(k)*gamma(7._r8+pgam(k))* &
!                  !   bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
!                  !   lamc(k)**3/lamc(k)**3
!
!                  !nnuccc(k) = cons10/(cons3*qcvar)* &
!                  !   pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
!                  !   *bimm* &
!                  !   (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3
!               end if           ! sub-columns
!
!         !if ( it .eq. 1 ) then
!         !    !qrout2(k) = qrout2(k)+nnuccc(k)
!         !    !qsout2(k) = qsout2(k)+mnuccc(k)
!         !    !nrout2(k) = nrout2(k)+nnucct(k)
!         !    !nsout2(k) = nsout2(k)+mnucct(k)
!         !    drout2(i,k) = drout2(i,k)+nnuccc(k)
!         !    dsout2(i,k) = dsout2(i,k)+mnuccc(k)
!         !endif
!
!               ! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
!               ! dust size and number in 4 bins are read in from companion routine
!               tmp_array(1) = 270.16d0-t(i,k)
!               tmp_array(2) = 1.3d0
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tcnt = tmp_array(1)
!
!               !tcnt=(270.16_r8-t(i,k))**1.3_r8
!
!               tmp_array(1) =(t(i,k)/298.0_r8)
!               tmp_array(2) = 0.85_r8
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               viscosity=1.8e-5_r8*tmp1    ! Viscosity (kg/m/s)
!
!               !viscosity=1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
!               mfp=2.0_r8*viscosity/(p(i,k)  &                   ! Mean free path (m)
!                  *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k))))           
!
!               nslip1=1.0_r8+(mfp/rndst(i,k,1))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,1)/mfp))))! Slip correction factor
!               nslip2=1.0_r8+(mfp/rndst(i,k,2))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,2)/mfp))))
!               nslip3=1.0_r8+(mfp/rndst(i,k,3))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,3)/mfp))))
!               nslip4=1.0_r8+(mfp/rndst(i,k,4))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,4)/mfp))))
!
!               ndfaer1=1.381e-23_r8*t(i,k)*nslip1/(6._r8*pi*viscosity*rndst(i,k,1))  ! aerosol diffusivity (m2/s)
!               ndfaer2=1.381e-23_r8*t(i,k)*nslip2/(6._r8*pi*viscosity*rndst(i,k,2))
!               ndfaer3=1.381e-23_r8*t(i,k)*nslip3/(6._r8*pi*viscosity*rndst(i,k,3))
!               ndfaer4=1.381e-23_r8*t(i,k)*nslip4/(6._r8*pi*viscosity*rndst(i,k,4))
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+nslip1
!         !    qsout2(i,k) = qsout2(i,k)+nslip2
!         !    nrout2(i,k) = nrout2(i,k)+nslip3
!         !    nsout2(i,k) = nsout2(i,k)+ndfaer1
!         !    drout2(i,k) = drout2(i,k)+ndfaer2
!         !    dsout2(i,k) = dsout2(i,k)+ndfaer3
!         !endif
!
!
!               if (microp_uniform) then
!
!                  !mnucct(k) = &
!                  !   (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                  !   ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
!                  !   cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
!
!                  mnucct(k) = &
!                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
!                     cdist1(k)*gamma(pgam(k)+5._r8)/(lamc(k)*lamc(k)*lamc(k)*lamc(k))
!
!                  nnucct(k) = (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
!                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
!
!               else
!
!                  !mnucct(k) = gamma(qcvar+4._r8/3._r8)/(cons3*qcvar**(4._r8/3._r8))*  &
!                  !   (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                  !   ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
!                  !   cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
!
!                  tmp1 = qcvar+4.d0/3.d0
!                  call athread_spawn(slave_gamma_parallel, tmp1)
!                  call athread_join()
!                  tmp_array(1) = qcvar
!                  tmp_array(2) = 4.d0/3.d0
!                  call athread_spawn(slave_pow_parallel, tmp_array)
!                  call athread_join()
!                  tmp2 = pgam(k)+5.d0
!                  call athread_spawn(slave_gamma_parallel, tmp2)
!                  call athread_join()
!                  mnucct(k) = tmp1/(cons3*tmp_array(1))*  &
!                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
!                     cdist1(k)*tmp2/(lamc(k)*lamc(k)*lamc(k)*lamc(k))
!
!
!                  tmp1 = qcvar+1.d0/3.d0
!                  call athread_spawn(slave_gamma_parallel, tmp1)
!                  call athread_join()
!                  tmp_array(1) = qcvar
!                  tmp_array(2) = 1.d0/3.d0
!                  call athread_spawn(slave_pow_parallel, tmp_array)
!                  call athread_join()
!                  tmp2 = pgam(k)+2.d0
!                  call athread_spawn(slave_gamma_parallel, tmp2)
!                  call athread_join()
!                  nnucct(k) =  tmp1/(cons3*tmp_array(1))*  &
!                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
!                     cdist1(k)*tmp2/lamc(k)
!
!                  !mnucct(k) = gamma(qcvar+4._r8/3._r8)/(cons3*qcvar**(4._r8/3._r8))*  &
!                  !   (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                  !   ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
!                  !   cdist1(k)*gamma(pgam(k)+5._r8)/(lamc(k)*lamc(k)*lamc(k)*lamc(k))
!
!                  !nnucct(k) =  gamma(qcvar+1._r8/3._r8)/(cons3*qcvar**(1._r8/3._r8))*  &
!                  !   (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
!                  !   ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
!                  !   cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+nacon(i,k,1)
!         !    qsout2(i,k) = qsout2(i,k)+nacon(i,k,2)
!         !    nrout2(i,k) = nrout2(i,k)+mnucct(k) 
!         !    nsout2(i,k) = nsout2(i,k)+nnucct(k) 
!         !    drout2(i,k) = drout2(i,k)+nacon(i,k,3)
!         !    dsout2(i,k) = dsout2(i,k)+nacon(i,k,4)
!         !endif
!
!               end if      ! sub-column switch
!
!               ! make sure number of droplets frozen does not exceed available ice nuclei concentration
!               ! this prevents 'runaway' droplet freezing
!
!               if (nnuccc(k)*lcldm(i,k).gt.nnuccd(k)) then
!                  dum=(nnuccd(k)/(nnuccc(k)*lcldm(i,k)))
!                  ! scale mixing ratio of droplet freezing with limit
!                  mnuccc(k)=mnuccc(k)*dum
!                  nnuccc(k)=nnuccd(k)/lcldm(i,k)
!               end if
!
!            else
!               mnuccc(k)=0._r8
!               nnuccc(k)=0._r8
!               mnucct(k)=0._r8
!               nnucct(k)=0._r8
!            end if
!
!         else
!            if (do_cldice .and. qcic(i,k) >= qsmall) then
!                tmp_array(1) = (1.333_r8*pi)
!                tmp_array(2) = 0.333_r8
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!               con1 = 1._r8/tmp1
!
!                tmp_array(1) =(rho(i,k)*qcic(i,k)/(rhow*max(ncic(i,k)*rho(i,k), 1.0e6_r8)))
!                tmp_array(2) = 0.333_r8
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!
!               r3lx = con1*tmp1 ! in m
!
!               !con1 = 1._r8/(1.333_r8*pi)**0.333_r8
!               !r3lx = con1*(rho(i,k)*qcic(i,k)/(rhow*max(ncic(i,k)*rho(i,k), 1.0e6_r8)))**0.333_r8 ! in m
!               r3lx = max(4.e-6_r8, r3lx)
!               mi0l = 4._r8/3._r8*pi*rhow*r3lx*r3lx*r3lx
!               !mi0l = 4._r8/3._r8*pi*rhow*r3lx**3_r8
!                
!               nnuccc(k) = frzimm(i,k)*1.0e6_r8/rho(i,k)
!               mnuccc(k) = nnuccc(k)*mi0l 
!
!               nnucct(k) = frzcnt(i,k)*1.0e6_r8/rho(i,k)
!               mnucct(k) = nnucct(k)*mi0l 
!
!               nnudep(k) = frzdep(i,k)*1.0e6_r8/rho(i,k)
!               mnudep(k) = nnudep(k)*mi0
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+mnuccc(k)
!         !    qsout2(i,k) = qsout2(i,k)+nnuccc(k) 
!         !    nrout2(i,k) = nrout2(i,k)+mnucct(k) 
!         !    nsout2(i,k) = nsout2(i,k)+nnucct(k) 
!         !    drout2(i,k) = drout2(i,k)+mnudep(k) 
!         !    dsout2(i,k) = dsout2(i,k)+nnudep(k) 
!         !endif
!
!            else
!               nnuccc(k) = 0._r8
!               mnuccc(k) = 0._r8
!
!               nnucct(k) = 0._r8
!               mnucct(k) = 0._r8
!
!               nnudep(k) = 0._r8
!               mnudep(k) = 0._r8
!            end if
!         endif
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+nnuccc(k)
!         !    qsout2(i,k) = qsout2(i,k)+mnuccc(k)
!         !    nrout2(i,k) = nrout2(i,k)+nnucct(k)
!         !    nsout2(i,k) = nsout2(i,k)+mnucct(k)
!         !    drout2(i,k) = drout2(i,k)+nnudep(k)
!         !    dsout2(i,k) = dsout2(i,k)+mnudep(k)
!         !endif
!         !.......................................................................
!         ! snow self-aggregation from passarelli, 1978, used by reisner, 1998
!         ! this is hard-wired for bs = 0.4 for now
!         ! ignore self-collection of cloud ice
!
!         if (qniic(i,k).ge.qsmall .and. t(i,k).le.273.15_r8) then
!            tmp_array(1) = pi
!            tmp_array(2) = ((1._r8-bs)/3._r8)
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp1 = tmp_array(1)
!
!            tmp_array(1) = rhosn
!            tmp_array(2) = ((-2._r8-bs)/3._r8)
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp2 = tmp_array(1)
!
!            tmp_array(1) = rho(i,k)
!            tmp_array(2) = ((2._r8+bs)/3._r8)
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp3 = tmp_array(1)
!
!            tmp_array(1) = qniic(i,k)
!            tmp_array(2) = ((2._r8+bs)/3._r8)
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp4 = tmp_array(1)
!
!            tmp_array(1) = (nsic(i,k)*rho(i,k))
!            tmp_array(2) = ((4._r8-bs)/3._r8)
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp5 = tmp_array(1)
!
!            nsagg(k) = -1108._r8*asn(i,k)*Eii* &
!                 tmp1*tmp2*tmp3*tmp4* &
!                 tmp5/ &
!                 (4._r8*720._r8*rho(i,k))
!            !nsagg(k) = -1108._r8*asn(i,k)*Eii* &
!            !     pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)*rho(i,k)** &
!            !     ((2._r8+bs)/3._r8)*qniic(i,k)**((2._r8+bs)/3._r8)* &
!            !     (nsic(i,k)*rho(i,k))**((4._r8-bs)/3._r8)/ &
!            !     (4._r8*720._r8*rho(i,k))
!         else
!            nsagg(k)=0._r8
!         end if
!
!
!         !.......................................................................
!         ! accretion of cloud droplets onto snow/graupel
!         ! here use continuous collection equation with
!         ! simple gravitational collection kernel
!         ! ignore collisions between droplets/cloud ice
!         ! since minimum size ice particle for accretion is 50 - 150 micron
!
!         ! ignore collision of snow with droplets above freezing
!
!         if (qniic(i,k).ge.qsmall .and. t(i,k).le.tmelt .and. &
!              qcic(i,k).ge.qsmall) then
!
!            ! put in size dependent collection efficiency
!            ! mean diameter of snow is area-weighted, since
!            ! accretion is function of crystal geometric area
!            ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)
!
!            dc0 = (pgam(k)+1._r8)/lamc(k)
!            ds0 = 1._r8/lams(k)
!            dum = dc0*dc0*uns(k)*rhow/(9._r8*mu(i,k)*ds0)
!            eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))
!
!            eci = max(eci,0._r8)
!            eci = min(eci,1._r8)
!
!
!            ! no impact of sub-grid distribution of qc since psacws
!            ! is linear in qc
!            tmp_array(1) = lams(k)
!            tmp_array(2) = bs+3._r8
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp1 = tmp_array(1)
!            psacws(k) = pi/4._r8*asn(i,k)*qcic(i,k)*rho(i,k)* &
!                 n0s(k)*Eci*cons11/ &
!                 tmp1
!            npsacws(k) = pi/4._r8*asn(i,k)*ncic(i,k)*rho(i,k)* &
!                 n0s(k)*Eci*cons11/ &
!                 tmp1
!
!            !psacws(k) = pi/4._r8*asn(i,k)*qcic(i,k)*rho(i,k)* &
!            !     n0s(k)*Eci*cons11/ &
!            !     lams(k)**(bs+3._r8)
!            !npsacws(k) = pi/4._r8*asn(i,k)*ncic(i,k)*rho(i,k)* &
!            !     n0s(k)*Eci*cons11/ &
!            !     lams(k)**(bs+3._r8)
!         else
!            psacws(k)=0._r8
!            npsacws(k)=0._r8
!         end if
!
!         ! add secondary ice production due to accretion of droplets by snow 
!         ! (Hallet-Mossop process) (from Cotton et al., 1986)
!
!         if (.not. do_cldice) then
!            ni_secp   = 0.0_r8
!            nsacwi(k) = 0.0_r8
!            msacwi(k) = 0.0_r8
!         else if((t(i,k).lt.270.16_r8) .and. (t(i,k).ge.268.16_r8)) then
!            ni_secp   = 3.5e8_r8*(270.16_r8-t(i,k))/2.0_r8*psacws(k)
!            nsacwi(k) = ni_secp
!            msacwi(k) = min(ni_secp*mi0,psacws(k))
!         else if((t(i,k).lt.268.16_r8) .and. (t(i,k).ge.265.16_r8)) then
!            ni_secp   = 3.5e8_r8*(t(i,k)-265.16_r8)/3.0_r8*psacws(k)
!            nsacwi(k) = ni_secp
!            msacwi(k) = min(ni_secp*mi0,psacws(k))
!         else
!            ni_secp   = 0.0_r8
!            nsacwi(k) = 0.0_r8
!            msacwi(k) = 0.0_r8
!         endif
!         psacws(k) = max(0.0_r8,psacws(k)-ni_secp*mi0)
!
!         !.......................................................................
!         ! accretion of rain water by snow
!         ! formula from ikawa and saito, 1991, used by reisner et al., 1998
!
!         if (qric(i,k).ge.1.e-8_r8 .and. qniic(i,k).ge.1.e-8_r8 .and. & 
!              t(i,k).le.273.15_r8) then
!
!            !pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))**2+ &
!            !     0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(i,k)* &
!            !     n0r(k)*n0s(k)* &
!            !     (5._r8/(lamr(k)**6*lams(k))+ &
!            !     2._r8/(lamr(k)**5*lams(k)**2)+ &
!            !     0.5_r8/(lamr(k)**4*lams(k)**3)))
!
!            !npracs(k) = pi/2._r8*rho(i,k)*ecr*(1.7_r8*(unr(k)-uns(k))**2+ &
!            !     0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
!            !     (1._r8/(lamr(k)**3*lams(k))+ &
!            !     1._r8/(lamr(k)**2*lams(k)**2)+ &
!            !     1._r8/(lamr(k)*lams(k)**3))
!
!            tmp_array(1) = ((1.2_r8*umr(k)-0.95_r8*ums(k))*(1.2_r8*umr(k)-0.95_r8*ums(k))+ &
!                 0.08_r8*ums(k)*umr(k))
!            tmp_array(2) = 0.5d0
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp1 = tmp_array(1)
!            pracs(k) = pi*pi*ecr*(tmp1*rhow*rho(i,k)* &
!                 n0r(k)*n0s(k)* &
!                 (5._r8/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k))+ &
!                 2._r8/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k)*lams(k))+ &
!                 0.5_r8/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k)*lams(k)*lams(k))))
!            tmp_array(1) = (1.7_r8*(unr(k)-uns(k))*(unr(k)-uns(k))+ &
!                 0.3_r8*unr(k)*uns(k))
!            tmp_array(2) = 0.5d0
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp1 = tmp_array(1)
!            npracs(k) = pi/2._r8*rho(i,k)*ecr*tmp1*n0r(k)*n0s(k)* &
!                 (1._r8/(lamr(k)*lamr(k)*lamr(k)*lams(k))+ &
!                 1._r8/(lamr(k)*lamr(k)*lams(k)*lams(k))+ &
!                 1._r8/(lamr(k)*lams(k)*lams(k)*lams(k)))
!
!
!            !pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))*(1.2_r8*umr(k)-0.95_r8*ums(k))+ &
!            !     0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(i,k)* &
!            !     n0r(k)*n0s(k)* &
!            !     (5._r8/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k))+ &
!            !     2._r8/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k)*lams(k))+ &
!            !     0.5_r8/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k)*lams(k)*lams(k))))
!
!            !npracs(k) = pi/2._r8*rho(i,k)*ecr*(1.7_r8*(unr(k)-uns(k))*(unr(k)-uns(k))+ &
!            !     0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
!            !     (1._r8/(lamr(k)*lamr(k)*lamr(k)*lams(k))+ &
!            !     1._r8/(lamr(k)*lamr(k)*lams(k)*lams(k))+ &
!            !     1._r8/(lamr(k)*lams(k)*lams(k)*lams(k)))
!
!         else
!            pracs(k)=0._r8
!            npracs(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! heterogeneous freezing of rain drops
!         ! follows from Bigg (1953)
!
!         if (t(i,k).lt.269.15_r8 .and. qric(i,k).ge.qsmall) then
!
!            !!mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
!            !!     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamr(k)**3 &
!            !!     /lamr(k)**3
!
!            !mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
!            !     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/(lamr(k)*lamr(k)*lamr(k)) &
!            !     /(lamr(k)*lamr(k)*lamr(k))
!
!            !!nnuccr(k) = pi*nric(i,k)*bimm* &
!            !!     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamr(k)**3
!
!            !nnuccr(k) = pi*nric(i,k)*bimm* &
!            !     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/(lamr(k)*lamr(k)*lamr(k))
!
!
!            tmp1 = aimm*(273.15d0-t(i,k))
!            call athread_spawn(slave_exp_parallel, tmp1)
!            call athread_join()
!            mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
!                 (tmp1-1._r8)/(lamr(k)*lamr(k)*lamr(k)) &
!                 /(lamr(k)*lamr(k)*lamr(k))
!
!            nnuccr(k) = pi*nric(i,k)*bimm* &
!                 (tmp1-1._r8)/(lamr(k)*lamr(k)*lamr(k))
!         else
!            mnuccr(k)=0._r8
!            nnuccr(k)=0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+nsagg(k)
!         !    !qsout2(i,k) = qsout2(i,k)+npsacws(k)
!         !    !nrout2(i,k) = nrout2(i,k)+pracs(k)
!         !    !nsout2(i,k) = nsout2(i,k)+npracs(k)
!         !    !drout2(i,k) = drout2(i,k)+nnuccr(k)
!         !    dsout2(i,k) = dsout2(i,k)+mnuccr(k)
!         !endif
!         !.......................................................................
!         ! accretion of cloud liquid water by rain
!         ! formula from Khrouditnov and Kogan (2000)
!         ! gravitational collection kernel, droplet fall speed neglected
!
!         if (qric(i,k).ge.qsmall .and. qcic(i,k).ge.qsmall) then
!
!            ! include sub-grid distribution of cloud water
!
!            ! add sub-column switch
!
!            if (microp_uniform) then
!
!                tmp_array(1) = (qcic(i,k)*qric(i,k))
!                tmp_array(2) = 1.15_r8
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!               pra(k) = 67._r8*tmp1
!
!               !pra(k) = 67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
!               npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
!
!            else
!
!                tmp_array(1) = (qcic(i,k)*qric(i,k))
!                tmp_array(2) = 1.15_r8
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!               pra(k) = accre_enhan(i,k)*(cons12/(cons3*cons20)*67._r8*tmp1)
!
!               !pra(k) = accre_enhan(i,k)*(cons12/(cons3*cons20)*67._r8*(qcic(i,k)*qric(i,k))**1.15_r8)
!               npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
!
!            end if               ! sub-column switch
!
!         else
!            pra(k)=0._r8
!            npra(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! Self-collection of rain drops
!         ! from Beheng(1994)
!
!         if (qric(i,k).ge.qsmall) then
!            nragg(k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
!         else
!            nragg(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! Accretion of cloud ice by snow
!         ! For this calculation, it is assumed that the Vs >> Vi
!         ! and Ds >> Di for continuous collection
!
!         if (do_cldice .and. qniic(i,k).ge.qsmall.and.qiic(i,k).ge.qsmall &
!              .and.t(i,k).le.273.15_r8) then
!
!                tmp_array(1) = lams(k)
!                tmp_array(2) = bs+3._r8
!                call athread_spawn(slave_pow_parallel, tmp_array)
!                call athread_join()
!                tmp1 = tmp_array(1)
!            prai(k) = pi/4._r8*asn(i,k)*qiic(i,k)*rho(i,k)* &
!                 n0s(k)*Eii*cons11/ &
!                 tmp1
!            nprai(k) = pi/4._r8*asn(i,k)*niic(i,k)* &
!                 rho(i,k)*n0s(k)*Eii*cons11/ &
!                 tmp1
!
!            !prai(k) = pi/4._r8*asn(i,k)*qiic(i,k)*rho(i,k)* &
!            !     n0s(k)*Eii*cons11/ &
!            !     lams(k)**(bs+3._r8)
!            !nprai(k) = pi/4._r8*asn(i,k)*niic(i,k)* &
!            !     rho(i,k)*n0s(k)*Eii*cons11/ &
!            !     lams(k)**(bs+3._r8)
!         else
!            prai(k)=0._r8
!            nprai(k)=0._r8
!         end if
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! calculate evaporation/sublimation of rain and snow
!         ! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
!         ! in-cloud condensation/deposition of rain and snow is neglected
!         ! except for transfer of cloud water to snow through bergeron process
!
!         ! initialize evap/sub tendncies
!         pre(k)=0._r8
!         prds(k)=0._r8
!
!         ! evaporation of rain
!         ! only calculate if there is some precip fraction > cloud fraction
!
!         if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8.or.cldmax(i,k).gt.lcldm(i,k)) then
!
!            ! set temporary cloud fraction to zero if cloud water + ice is very small
!            ! this will ensure that evaporation/sublimation of precip occurs over
!            ! entire grid cell, since min cloud fraction is specified otherwise
!            if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8) then
!               dum=0._r8
!            else
!               dum=lcldm(i,k)
!            end if
!
!            ! saturation vapor pressure
!      tmp1 = svp_tboil/t(i,k)
!      call athread_spawn(slave_log10_parallel, tmp1)
!      call athread_join()
!
!      esn = 10._r8**(-7.90298_r8*(svp_tboil/t(i,k)-1._r8)+ &
!           5.02808_r8*tmp1- &
!           1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t(i,k)/svp_tboil))-1._r8)+ &
!           8.1328e-3_r8*(10._r8**(-3.49149_r8*(svp_tboil/t(i,k)-1._r8))-1._r8)+ &
!           log10(1013.246_r8))*100._r8
!
!            !esn=svp_water(t(i,k))
!            qsn=svp_to_qsat(esn, p(i,k))
!
!            ! recalculate saturation vapor pressure for liquid and ice
!            esl(i,k)=esn
!
!      tmp2 = svp_h2otrip/t(i,k)
!      call athread_spawn(slave_log10_parallel, tmp2)
!      call athread_join()
!
!      esi(i,k) = 10.d0**(-9.09718d0*(svp_h2otrip/t(i,k)-1.d0)-3.56654d0* &
!           tmp2+0.876793d0*(1.d0-t(i,k)/svp_h2otrip)+ &
!           log10(6.1071d0))*100.d0
!
!            !esi(i,k)=svp_ice(t(i,k))
!            ! hm fix, make sure when above freezing that esi=esl, not active yet
!            if (t(i,k).gt.tmelt)esi(i,k)=esl(i,k)
!
!            ! calculate q for out-of-cloud region
!            qclr=(q(i,k)-dum*qsn)/(1._r8-dum)
!
!            if (qric(i,k).ge.qsmall) then
!
!               qvs=svp_to_qsat(esl(i,k), p(i,k))
!               dqsdt = xxlv*qvs/(rv*t(i,k)*t(i,k))
!               ab = 1._r8+dqsdt*xxlv/cpp
!
!               tmp_array(1) = (arn(i,k)*rho(i,k)/mu(i,k))
!               tmp_array(2) = 0.5d0
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               tmp_array(1) = sc(i,k)
!               tmp_array(2) =(1._r8/3._r8)
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp2 = tmp_array(1)
!               tmp_array(1) = lamr(k)
!               tmp_array(2) = (5._r8/2._r8+br/2._r8)
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp3 = tmp_array(1)
!               epsr = 2._r8*pi*n0r(k)*rho(i,k)*Dv(i,k)* &
!                    (f1r/(lamr(k)*lamr(k))+ &
!                    f2r*tmp1* &
!                    tmp2*cons13/ &
!                    (tmp3))
!
!               !epsr = 2._r8*pi*n0r(k)*rho(i,k)*Dv(i,k)* &
!               !     (f1r/(lamr(k)*lamr(k))+ &
!               !     f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
!               !     sc(i,k)**(1._r8/3._r8)*cons13/ &
!               !     (lamr(k)**(5._r8/2._r8+br/2._r8)))
!
!               pre(k) = epsr*(qclr-qvs)/ab
!
!               ! only evaporate in out-of-cloud region
!               ! and distribute across cldmax
!               pre(k)=min(pre(k)*(cldmax(i,k)-dum),0._r8)
!               pre(k)=pre(k)/cldmax(i,k)
!               am_evp_st(i,k) = max(cldmax(i,k)-dum, 0._r8)
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+epsr
!         !    qsout2(i,k) = qsout2(i,k)+qclr
!         !    nrout2(i,k) = nrout2(i,k)+qvs
!         !    nsout2(i,k) = nsout2(i,k)+qclr-qvs
!         !    drout2(i,k) = drout2(i,k)+dum
!         !    dsout2(i,k) = dsout2(i,k)+pre(k)
!         !endif
!            end if
!
!            ! sublimation of snow
!            if (qniic(i,k).ge.qsmall) then
!               qvi=svp_to_qsat(esi(i,k), p(i,k))
!               dqsidt =  xxls*qvi/(rv*t(i,k)*t(i,k))
!               abi = 1._r8+dqsidt*xxls/cpp
!
!               tmp_array(1) = (asn(i,k)*rho(i,k)/mu(i,k))
!               tmp_array(2) = 0.5d0
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               tmp_array(1) = sc(i,k)
!               tmp_array(2) =(1._r8/3._r8)
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp2 = tmp_array(1)
!               tmp_array(1) = lams(k)
!               tmp_array(2) = (5._r8/2._r8+bs/2._r8)
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp3 = tmp_array(1)
!               epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
!                    (f1s/(lams(k)*lams(k))+ &
!                    f2s*tmp1* &
!                    tmp2*cons14/ &
!                    (tmp3))
!
!               !epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
!               !     (f1s/(lams(k)*lams(k))+ &
!               !     f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
!               !     sc(i,k)**(1._r8/3._r8)*cons14/ &
!               !     (lams(k)**(5._r8/2._r8+bs/2._r8)))
!               prds(k) = epss*(qclr-qvi)/abi
!
!               ! only sublimate in out-of-cloud region and distribute over cldmax
!               prds(k)=min(prds(k)*(cldmax(i,k)-dum),0._r8)
!               prds(k)=prds(k)/cldmax(i,k)
!               am_evp_st(i,k) = max(cldmax(i,k)-dum, 0._r8)
!            end if
!
!            ! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
!            ! get updated RH at end of time step based on cloud water/ice condensation/evap
!
!            qtmp=q(i,k)-(cmei(i,k)+(pre(k)+prds(k))*cldmax(i,k))*deltat
!            ttmp=t(i,k)+((pre(k)*cldmax(i,k))*xxlv+ &
!                 (cmei(i,k)+prds(k)*cldmax(i,k))*xxls)*deltat/cpp
!
!            !limit range of temperatures!
!            ttmp=max(180._r8,min(ttmp,323._r8))
!
!      tmp1 = svp_tboil/ttmp
!      call athread_spawn(slave_log10_parallel, tmp1)
!      call athread_join()
!
!      esn = 10._r8**(-7.90298_r8*(svp_tboil/ttmp-1._r8)+ &
!           5.02808_r8*tmp1- &
!           1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-ttmp/svp_tboil))-1._r8)+ &
!           8.1328e-3_r8*(10._r8**(-3.49149_r8*(svp_tboil/ttmp-1._r8))-1._r8)+ &
!           log10(1013.246_r8))*100._r8
!
!            !esn=svp_water(ttmp)  ! use rhw to allow ice supersaturation
!            qsn=svp_to_qsat(esn, p(i,k))
!
!            ! modify precip evaporation rate if q > qsat
!            if (qtmp.gt.qsn) then
!               if (pre(k)+prds(k).lt.-1.e-20_r8) then
!                  dum1=pre(k)/(pre(k)+prds(k))
!                  ! recalculate q and t after cloud water cond but without precip evap
!                  qtmp=q(i,k)-(cmei(i,k))*deltat
!                  ttmp=t(i,k)+(cmei(i,k)*xxls)*deltat/cpp
!
!      tmp1 = svp_tboil/ttmp
!      call athread_spawn(slave_log10_parallel, tmp1)
!      call athread_join()
!
!      esn = 10._r8**(-7.90298_r8*(svp_tboil/ttmp-1._r8)+ &
!           5.02808_r8*tmp1- &
!           1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-ttmp/svp_tboil))-1._r8)+ &
!           8.1328e-3_r8*(10._r8**(-3.49149_r8*(svp_tboil/ttmp-1._r8))-1._r8)+ &
!           log10(1013.246_r8))*100._r8
!
!                  !esn=svp_water(ttmp) ! use rhw to allow ice supersaturation
!                  qsn=svp_to_qsat(esn, p(i,k))
!                  dum=(qtmp-qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp*ttmp))
!                  dum=min(dum,0._r8)
!
!                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
!                  pre(k)=dum*dum1/deltat/cldmax(i,k)
!
!                  ! do separately using RHI for prds....
!
!      tmp2 = svp_h2otrip/ttmp
!      call athread_spawn(slave_log10_parallel, tmp2)
!      call athread_join()
!
!      esn = 10.d0**(-9.09718d0*(svp_h2otrip/ttmp-1.d0)-3.56654d0* &
!           tmp2+0.876793d0*(1.d0-ttmp/svp_h2otrip)+ &
!           log10(6.1071d0))*100.d0
!
!                  !esn=svp_ice(ttmp) ! use rhi to allow ice supersaturation
!                  qsn=svp_to_qsat(esn, p(i,k))
!                  dum=(qtmp-qsn)/(1._r8 + cons28*qsn/(cpp*rv*ttmp*ttmp))
!                  dum=min(dum,0._r8)
!
!                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
!                  prds(k)=dum*(1._r8-dum1)/deltat/cldmax(i,k)
!               end if
!            end if
!         end if
!
!         ! bergeron process - evaporation of droplets and deposition onto snow
!
!         if (qniic(i,k).ge.qsmall.and.qcic(i,k).ge.qsmall.and.t(i,k).lt.tmelt) then
!            qvi=svp_to_qsat(esi(i,k), p(i,k))
!            qvs=svp_to_qsat(esl(i,k), p(i,k))
!            dqsidt =  xxls*qvi/(rv*t(i,k)*t(i,k))
!            abi = 1._r8+dqsidt*xxls/cpp
!
!            tmp_array(1) = (asn(i,k)*rho(i,k)/mu(i,k))
!            tmp_array(2) = 0.5_r8
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp1 = tmp_array(1)
!            tmp_array(1) = sc(i,k)
!            tmp_array(2) = (1._r8/3._r8)
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp2 = tmp_array(1)
!            tmp_array(1) = lams(k)
!            tmp_array(2) = (5._r8/2._r8+bs/2._r8)
!            call athread_spawn(slave_pow_parallel, tmp_array)
!            call athread_join()
!            tmp3 = tmp_array(1)
!            epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
!                 (f1s/(lams(k)*lams(k))+ &
!                 f2s*tmp1* &
!                 tmp2*cons14/ &
!                 (tmp3))
!
!            !epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
!            !     (f1s/(lams(k)*lams(k))+ &
!            !     f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
!            !     sc(i,k)**(1._r8/3._r8)*cons14/ &
!            !     (lams(k)**(5._r8/2._r8+bs/2._r8)))
!            bergs(k)=epss*(qvs-qvi)/abi
!         else
!            bergs(k)=0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+npra(k)
!         !    qsout2(i,k) = qsout2(i,k)+nprai(k)
!         !    nrout2(i,k) = nrout2(i,k)+prai(k)
!         !    nsout2(i,k) = nsout2(i,k)+pre(k)
!         !    drout2(i,k) = drout2(i,k)+prds(k)
!         !    dsout2(i,k) = dsout2(i,k)+bergs(k)
!         !endif
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! conservation to ensure no negative values of cloud water/precipitation
!         ! in case microphysical process rates are large
!
!         ! make sure and use end-of-time step values for cloud water, ice, due
!         ! condensation/deposition
!
!         ! note: for check on conservation, processes are multiplied by omsm
!         ! to prevent problems due to round off error
!
!         ! include mixing timescale  (mtime)
!
!         qce=(qc(i,k) - berg(i,k)*deltat)
!         nce=(nc(i,k)+npccn(k)*deltat*mtime)
!         qie=(qi(i,k)+(cmei(i,k)+berg(i,k))*deltat)
!         nie=(ni(i,k)+nnuccd(k)*deltat*mtime)
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+qce
!         !    qsout2(i,k) = qsout2(i,k)+nce
!         !    nrout2(i,k) = nrout2(i,k)+qie
!         !    nsout2(i,k) = nsout2(i,k)+nie
!         !    !drout2(i,k) = drout2(i,k)+prds(k)
!         !    !dsout2(i,k) = dsout2(i,k)+bergs(k)
!         !endif
!         ! conservation of qc
!
!         dum = (prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+ &
!              psacws(k)+bergs(k))*lcldm(i,k)*deltat
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+prc(k)
!         !    qsout2(i,k) = qsout2(i,k)+pra(k)
!         !    nrout2(i,k) = nrout2(i,k)+mnuccc(k)
!         !    nsout2(i,k) = nsout2(i,k)+mnucct(k)
!         !    drout2(i,k) = drout2(i,k)+msacwi(k)
!         !    dsout2(i,k) = dsout2(i,k)+psacws(k)
!         !endif
!
!         if (dum.gt.qce) then
!            ratio = qce/deltat/lcldm(i,k)/(prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+psacws(k)+bergs(k))*omsm 
!
!            prc(k) = prc(k)*ratio
!            pra(k) = pra(k)*ratio
!            mnuccc(k) = mnuccc(k)*ratio
!            mnucct(k) = mnucct(k)*ratio  
!            msacwi(k) = msacwi(k)*ratio  
!            psacws(k) = psacws(k)*ratio
!            bergs(k) = bergs(k)*ratio
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+prc(k)
!         !    qsout2(i,k) = qsout2(i,k)+pra(k)
!         !    nrout2(i,k) = nrout2(i,k)+mnuccc(k)
!         !    nsout2(i,k) = nsout2(i,k)+mnucct(k)
!         !    drout2(i,k) = drout2(i,k)+msacwi(k)
!         !    dsout2(i,k) = dsout2(i,k)+psacws(k)
!         !endif
!
!         ! conservation of nc
!
!         dum = (nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
!              npsacws(k)-nsubc(k))*lcldm(i,k)*deltat
!
!         if (dum.gt.nce) then
!            ratio = nce/deltat/((nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+&
!                 npsacws(k)-nsubc(k))*lcldm(i,k))*omsm
!
!            nprc1(k) = nprc1(k)*ratio
!            npra(k) = npra(k)*ratio
!            nnuccc(k) = nnuccc(k)*ratio
!            nnucct(k) = nnucct(k)*ratio  
!            npsacws(k) = npsacws(k)*ratio
!            nsubc(k)=nsubc(k)*ratio
!         end if
!
!         ! conservation of qi
!
!         if (do_cldice) then
!
!            frztmp = -mnuccc(k) - mnucct(k) - msacwi(k)
!            if (use_hetfrz_classnuc) frztmp = -mnuccc(k)-mnucct(k)-mnudep(k)-msacwi(k)
!            dum = ( frztmp*lcldm(i,k) + (prci(k)+prai(k))*icldm(i,k) )*deltat
!
!            if (dum.gt.qie) then
!
!               frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
!               if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
!               ratio = (qie/deltat + frztmp*lcldm(i,k))/((prci(k)+prai(k))*icldm(i,k))*omsm 
!               prci(k) = prci(k)*ratio
!               prai(k) = prai(k)*ratio
!            end if
!
!            ! conservation of ni
!            frztmp = -nnucct(k) - nsacwi(k)
!            if (use_hetfrz_classnuc) frztmp = -nnucct(k) - nnuccc(k) - nnudep(k) - nsacwi(k)
!            dum = ( frztmp*lcldm(i,k) + (nprci(k)+nprai(k)-nsubi(k))*icldm(i,k) )*deltat
!
!            if (dum.gt.nie) then
!
!               frztmp = nnucct(k) + nsacwi(k)
!               if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
!               ratio = (nie/deltat + frztmp*lcldm(i,k))/ &  
!                     ((nprci(k)+nprai(k)-nsubi(k))*icldm(i,k))*omsm
!               nprci(k) = nprci(k)*ratio
!               nprai(k) = nprai(k)*ratio
!               nsubi(k) = nsubi(k)*ratio
!            end if
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+ratio
!         !    qsout2(i,k) = qsout2(i,k)+qnitend(i,k)
!         !    nrout2(i,k) = nrout2(i,k)+qrtot
!         !    nsout2(i,k) = nsout2(i,k)+nrtot
!         !    drout2(i,k) = drout2(i,k)+qstot
!         !    dsout2(i,k) = dsout2(i,k)+nstot
!         !endif
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+frztmp
!         !    qsout2(i,k) = qsout2(i,k)+mnucct(k)
!         !    nrout2(i,k) = nrout2(i,k)+nnucct(k)
!         !    nsout2(i,k) = nsout2(i,k)+prci(k)
!         !    drout2(i,k) = drout2(i,k)+prai(k)
!         !    dsout2(i,k) = dsout2(i,k)+nsubi(k)
!         !endif
!         ! for precipitation conservation, use logic that vertical integral 
!         ! of tendency from current level to top of model (i.e., qrtot) cannot be negative
!
!         ! conservation of rain mixing rat
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+pre(k)
!         !    qsout2(i,k) = qsout2(i,k)+pracs(k)
!         !    nrout2(i,k) = nrout2(i,k)+mnuccr(k)
!         !endif
!         if (((prc(k)+pra(k))*lcldm(i,k)+(-mnuccr(k)+pre(k)-pracs(k))*&
!              cldmax(i,k))*dz(i,k)*rho(i,k)+qrtot.lt.0._r8) then
!
!            if (-pre(k)+pracs(k)+mnuccr(k).ge.qsmall) then
!
!               ratio = (qrtot/(dz(i,k)*rho(i,k))+(prc(k)+pra(k))*lcldm(i,k))/&
!                    ((-pre(k)+pracs(k)+mnuccr(k))*cldmax(i,k))*omsm 
!
!               pre(k) = pre(k)*ratio
!               pracs(k) = pracs(k)*ratio
!               mnuccr(k) = mnuccr(k)*ratio
!            end if
!         end if
!         !if ( it .eq. 1 ) then
!         !    nsout2(i,k) = nsout2(i,k)+ qrtot
!         !    drout2(i,k) = drout2(i,k)+ pracs(k)
!         !    dsout2(i,k) = dsout2(i,k)+ mnuccr(k)
!         !endif
!
!         ! conservation of nr
!         ! for now neglect evaporation of nr
!         nsubr(k)=0._r8
!
!         if ((nprc(k)*lcldm(i,k)+(-nnuccr(k)+nsubr(k)-npracs(k)&
!              +nragg(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+nrtot.lt.0._r8) then
!
!            if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).ge.qsmall) then
!
!               ratio = (nrtot/(dz(i,k)*rho(i,k))+nprc(k)*lcldm(i,k))/&
!                    ((-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k))*cldmax(i,k))*omsm
!
!               nsubr(k) = nsubr(k)*ratio
!               npracs(k) = npracs(k)*ratio
!               nnuccr(k) = nnuccr(k)*ratio
!               nragg(k) = nragg(k)*ratio
!            end if
!         end if
!
!         ! conservation of snow mix ratio
!
!         if (((bergs(k)+psacws(k))*lcldm(i,k)+(prai(k)+prci(k))*icldm(i,k)+(pracs(k)+&
!              mnuccr(k)+prds(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+qstot.lt.0._r8) then
!
!            if (-prds(k).ge.qsmall) then
!
!               ratio = (qstot/(dz(i,k)*rho(i,k))+(bergs(k)+psacws(k))*lcldm(i,k)+(prai(k)+prci(k))*icldm(i,k)+&
!                    (pracs(k)+mnuccr(k))*cldmax(i,k))/(-prds(k)*cldmax(i,k))*omsm
!
!               prds(k) = prds(k)*ratio
!            end if
!         end if
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+ratio
!         !    qsout2(i,k) = qsout2(i,k)+prds(k)
!         !    nrout2(i,k) = nrout2(i,k)+qstot
!         !    nsout2(i,k) = nsout2(i,k)+psacws(k)
!         !    drout2(i,k) = drout2(i,k)+bergs(k)
!         !    dsout2(i,k) = dsout2(i,k)+dz(i,k)
!         !endif
!
!         ! conservation of ns
!
!         ! calculate loss of number due to sublimation
!         ! for now neglect sublimation of ns
!         nsubs(k)=0._r8
!
!         if ((nprci(k)*icldm(i,k)+(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax(i,k))*&
!              dz(i,k)*rho(i,k)+nstot.lt.0._r8) then
!
!            if (-nsubs(k)-nsagg(k).ge.qsmall) then
!
!               ratio = (nstot/(dz(i,k)*rho(i,k))+nprci(k)*icldm(i,k)+&
!                    nnuccr(k)*cldmax(i,k))/((-nsubs(k)-nsagg(k))*cldmax(i,k))*omsm
!
!               nsubs(k) = nsubs(k)*ratio
!               nsagg(k) = nsagg(k)*ratio
!            end if
!         end if
!
!         ! get tendencies due to microphysical conversion processes
!         ! note: tendencies are multiplied by appropaiate cloud/precip 
!         ! fraction to get grid-scale values
!         ! note: cmei is already grid-average values
!
!         qvlat(i,k) = qvlat(i,k)-(pre(k)+prds(k))*cldmax(i,k)-cmei(i,k) 
!
!         tlat(i,k) = tlat(i,k)+((pre(k)*cldmax(i,k)) &
!              *xxlv+(prds(k)*cldmax(i,k)+cmei(i,k))*xxls+ &
!              ((bergs(k)+psacws(k)+mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(i,k)+(mnuccr(k)+ &
!              pracs(k))*cldmax(i,k)+berg(i,k))*xlf)
!
!         qctend(i,k) = qctend(i,k)+ &
!              (-pra(k)-prc(k)-mnuccc(k)-mnucct(k)-msacwi(k)- & 
!              psacws(k)-bergs(k))*lcldm(i,k)-berg(i,k)
!
!         if (do_cldice) then
!
!            frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
!            if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
!            qitend(i,k) = qitend(i,k) + frztmp*lcldm(i,k) + &
!               (-prci(k)-prai(k))*icldm(i,k) + cmei(i,k) + berg(i,k)
!
!         end if
!
!         qrtend(i,k) = qrtend(i,k)+ &
!              (pra(k)+prc(k))*lcldm(i,k)+(pre(k)-pracs(k)- &
!              mnuccr(k))*cldmax(i,k)
!
!         qnitend(i,k) = qnitend(i,k)+ &
!              (prai(k)+prci(k))*icldm(i,k)+(psacws(k)+bergs(k))*lcldm(i,k)+(prds(k)+ &
!              pracs(k)+mnuccr(k))*cldmax(i,k)
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+qnitend(i,k)
!         !    qsout2(i,k) = qsout2(i,k)+icldm(i,k)
!         !    nrout2(i,k) = nrout2(i,k)+lcldm(i,k)
!         !    nsout2(i,k) = nsout2(i,k)+mnuccr(k)
!         !    drout2(i,k) = drout2(i,k)+psacws(k)
!         !    dsout2(i,k) = dsout2(i,k)+cldmax(i,k)
!         !    !qsout2(i,k) = qsout2(i,k)+prai(k)
!         !    !nrout2(i,k) = nrout2(i,k)+prci(k)
!         !    !nsout2(i,k) = nsout2(i,k)+bergs(k)
!         !    !drout2(i,k) = drout2(i,k)+prds(k)
!         !    !dsout2(i,k) = dsout2(i,k)+pracs(k)
!         !endif
!
!         ! add output for cmei (accumulate)
!         cmeiout(i,k) = cmeiout(i,k) + cmei(i,k)
!
!         ! assign variables for trop_mozart, these are grid-average
!         ! evaporation/sublimation is stored here as positive term
!
!         evapsnow(i,k) = evapsnow(i,k)-prds(k)*cldmax(i,k)
!         nevapr(i,k) = nevapr(i,k)-pre(k)*cldmax(i,k)
!         nevapr2(i,k) = nevapr2(i,k)-pre(k)*cldmax(i,k)
!
!         ! change to make sure prain is positive: do not remove snow from
!         ! prain used for wet deposition
!         prain(i,k) = prain(i,k)+(pra(k)+prc(k))*lcldm(i,k)+(-pracs(k)- &
!              mnuccr(k))*cldmax(i,k)
!         prodsnow(i,k) = prodsnow(i,k)+(prai(k)+prci(k))*icldm(i,k)+(psacws(k)+bergs(k))*lcldm(i,k)+(&
!              pracs(k)+mnuccr(k))*cldmax(i,k)
!
!         ! following are used to calculate 1st order conversion rate of cloud water
!         !    to rain and snow (1/s), for later use in aerosol wet removal routine
!         ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
!         !    used to calculate pra, prc, ... in this routine
!         ! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
!         !                      (no cloud ice or bergeron terms)
!         ! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }
!
!         qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k) + (pra(k)+prc(k)+psacws(k))*lcldm(i,k) 
!         qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc(i,k) 
!
!         ! microphysics output, note this is grid-averaged
!         prao(i,k)=prao(i,k)+pra(k)*lcldm(i,k)
!         prco(i,k)=prco(i,k)+prc(k)*lcldm(i,k)
!         mnuccco(i,k)=mnuccco(i,k)+mnuccc(k)*lcldm(i,k)
!         mnuccto(i,k)=mnuccto(i,k)+mnucct(k)*lcldm(i,k)
!         mnuccdo(i,k)=mnuccdo(i,k)+mnuccd(k)*lcldm(i,k)
!         msacwio(i,k)=msacwio(i,k)+msacwi(k)*lcldm(i,k)
!         psacwso(i,k)=psacwso(i,k)+psacws(k)*lcldm(i,k)
!         bergso(i,k)=bergso(i,k)+bergs(k)*lcldm(i,k)
!         bergo(i,k)=bergo(i,k)+berg(i,k)
!         prcio(i,k)=prcio(i,k)+prci(k)*icldm(i,k)
!         praio(i,k)=praio(i,k)+prai(k)*icldm(i,k)
!         mnuccro(i,k)=mnuccro(i,k)+mnuccr(k)*cldmax(i,k)
!         pracso (i,k)=pracso (i,k)+pracs (k)*cldmax(i,k)
!
!         ! multiply activation/nucleation by mtime to account for fast timescale
!
!         nctend(i,k) = nctend(i,k)+ npccn(k)*mtime+&
!              (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k) & 
!              -npra(k)-nprc1(k))*lcldm(i,k)      
!
!         if (do_cldice) then
!
!            frztmp = nnucct(k) + nsacwi(k)
!            if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
!            nitend(i,k) = nitend(i,k) + nnuccd(k)*mtime + & 
!                  frztmp*lcldm(i,k) + (nsubi(k)-nprci(k)-nprai(k))*icldm(i,k)
!
!         end if
!
!         nstend(i,k) = nstend(i,k)+(nsubs(k)+ &
!              nsagg(k)+nnuccr(k))*cldmax(i,k)+nprci(k)*icldm(i,k)
!
!
!         nrtend(i,k) = nrtend(i,k)+ &
!              nprc(k)*lcldm(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
!              +nragg(k))*cldmax(i,k)
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+nprc(k)*lcldm(i,k)
!         !    qsout2(i,k) = qsout2(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
!         !     +nragg(k))*cldmax(i,k)
!         !    nrout2(i,k) = nrout2(i,k)+nprc(k)*lcldm(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
!         !     +nragg(k))*cldmax(i,k)
!         !    nsout2(i,k) = nsout2(i,k)+nrtend(i,k)
!         !    !drout2(i,k) = drout2(i,k)+nsubr(k)-npracs(k)-nnuccr(k)
!         !    !dsout2(i,k) = dsout2(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
!         !    ! +nragg(k))*cldmax(i,k)
!         !endif
!
!         ! make sure that nc and ni at advanced time step do not exceed
!         ! maximum (existing N + source terms*dt), which is possible due to
!         ! fast nucleation timescale
!
!         if (nctend(i,k).gt.0._r8.and.nc(i,k)+nctend(i,k)*deltat.gt.ncmax) then
!            nctend(i,k)=max(0._r8,(ncmax-nc(i,k))/deltat)
!         end if
!
!         if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax) then
!            nitend(i,k)=max(0._r8,(nimax-ni(i,k))/deltat)
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k)+qrtend(i,k)
!         !    qsout2(i,k) = qsout2(i,k)+qnitend(i,k)
!         !    nrout2(i,k) = nrout2(i,k)+nstend(i,k)
!         !    nsout2(i,k) = nsout2(i,k)+nrtend(i,k)
!         !    drout2(i,k) = drout2(i,k)+nctend(i,k)
!         !    dsout2(i,k) = dsout2(i,k)+nitend(i,k)
!         !endif
!         ! get final values for precipitation q and N, based on
!         ! flux of precip from above, source/sink term, and terminal fallspeed
!         ! see eq. 15-16 in MG2008
!
!         ! rain
!
!         if (qric(i,k).ge.qsmall) then
!            if (k.eq.top_lev) then
!               qric(i,k)=qrtend(i,k)*dz(i,k)/cldmax(i,k)/umr(k)
!               nric(i,k)=nrtend(i,k)*dz(i,k)/cldmax(i,k)/unr(k)
!            else
!               qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
!                    (rho(i,k)*dz(i,k)*qrtend(i,k)))/(umr(k)*rho(i,k)*cldmax(i,k))
!               nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
!                    (rho(i,k)*dz(i,k)*nrtend(i,k)))/(unr(k)*rho(i,k)*cldmax(i,k))
!
!            end if
!         else
!            qric(i,k)=0._r8
!            nric(i,k)=0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    nrout2(i,k) = nrout2(i,k) + qric(i,k)
!         !endif
!
!         ! snow
!
!         if (qniic(i,k).ge.qsmall) then
!            if (k.eq.top_lev) then
!               qniic(i,k)=qnitend(i,k)*dz(i,k)/cldmax(i,k)/ums(k)
!               nsic(i,k)=nstend(i,k)*dz(i,k)/cldmax(i,k)/uns(k)
!            else
!               qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
!                    (rho(i,k)*dz(i,k)*qnitend(i,k)))/(ums(k)*rho(i,k)*cldmax(i,k))
!               nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
!                    (rho(i,k)*dz(i,k)*nstend(i,k)))/(uns(k)*rho(i,k)*cldmax(i,k))
!            end if
!         else
!            qniic(i,k)=0._r8
!            nsic(i,k)=0._r8
!         end if
!        !if ( it .eq. 1  ) then
!        !    qrout2(i,k) = qrout2(i,k)+qniic(i,k)
!        !    qsout2(i,k) = qsout2(i,k)+nsic(i,k)
!        !    !nrout2(i,k) = nrout2(i,k)+qc(i,k)
!        !    !nsout2(i,k) = nsout2(i,k)+qi(i,k)
!        !    !drout2(i,k) = drout2(i,k)+nc(i,k)
!        !    !dsout2(i,k) = dsout2(i,k)+ni(i,k)
!        !endif
!
!         ! calculate precipitation flux at surface
!         ! divide by density of water to get units of m/s
!
!                !qrout2(i,k) = qrtend(i,k)
!                !qsout2(i,k) = qnitend(i,k)
!                !nrout2(i,k) = prect(i)
!                !nsout2(i,k) = (qrtend(i,k)*dz(i,k)*rho(i,k)+qnitend(i,k)*dz(i,k)*rho(i,k))/rhow
!                !drout2(i,k) = qrtend(i,k)*dz(i,k)*rho(i,k)
!                !dsout2(i,k) = qnitend(i,k)*dz(i,k)*rho(i,k)
!
!         prect(i) = prect(i)+(qrtend(i,k)*dz(i,k)*rho(i,k)+&
!              qnitend(i,k)*dz(i,k)*rho(i,k))/rhow
!         preci(i) = preci(i)+qnitend(i,k)*dz(i,k)*rho(i,k)/rhow
!
!         ! convert rain rate from m/s to mm/hr
!
!         rainrt(i,k)=qric(i,k)*rho(i,k)*umr(k)/rhow*3600._r8*1000._r8
!
!         ! vertically-integrated precip source/sink terms (note: grid-averaged)
!
!         qrtot = max(qrtot+qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
!         qstot = max(qstot+qnitend(i,k)*dz(i,k)*rho(i,k),0._r8)
!         nrtot = max(nrtot+nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
!         nstot = max(nstot+nstend(i,k)*dz(i,k)*rho(i,k),0._r8)
!
!         ! calculate melting and freezing of precip
!
!         ! melt snow at +2 C
!
!         if (t(i,k)+tlat(i,k)/cpp*deltat > 275.15_r8) then
!            if (qstot > 0._r8) then
!
!               ! make sure melting snow doesn't reduce temperature below threshold
!               dum = -xlf/cpp*qstot/(dz(i,k)*rho(i,k))
!               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.275.15_r8) then
!                  dum = (t(i,k)+tlat(i,k)/cpp*deltat-275.15_r8)*cpp/xlf
!                  dum = dum/(xlf/cpp*qstot/(dz(i,k)*rho(i,k)))
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qric(i,k)=qric(i,k)+dum*qniic(i,k)
!               nric(i,k)=nric(i,k)+dum*nsic(i,k)
!               qniic(i,k)=(1._r8-dum)*qniic(i,k)
!               nsic(i,k)=(1._r8-dum)*nsic(i,k)
!               ! heating tendency 
!               tmp=-xlf*dum*qstot/(dz(i,k)*rho(i,k))
!               meltsdt(i,k)=meltsdt(i,k) + tmp
!
!               tlat(i,k)=tlat(i,k)+tmp
!               qrtot=qrtot+dum*qstot
!               nrtot=nrtot+dum*nstot
!               qstot=(1._r8-dum)*qstot
!               nstot=(1._r8-dum)*nstot
!               preci(i)=(1._r8-dum)*preci(i)
!            end if
!         end if
!         !if ( it .eq. 1 ) then
!         !    nsout2(i,k) = nsout2(i,k) + qric(i,k)
!         !endif
!
!         ! freeze all rain at -5C for Arctic
!
!         if (t(i,k)+tlat(i,k)/cpp*deltat < (tmelt - 5._r8)) then
!
!            if (qrtot > 0._r8) then
!
!               ! make sure freezing rain doesn't increase temperature above threshold
!               dum = xlf/cpp*qrtot/(dz(i,k)*rho(i,k))
!               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.(tmelt - 5._r8)) then
!                  dum = -(t(i,k)+tlat(i,k)/cpp*deltat-(tmelt-5._r8))*cpp/xlf
!                  dum = dum/(xlf/cpp*qrtot/(dz(i,k)*rho(i,k)))
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qniic(i,k)=qniic(i,k)+dum*qric(i,k)
!               nsic(i,k)=nsic(i,k)+dum*nric(i,k)
!               qric(i,k)=(1._r8-dum)*qric(i,k)
!               nric(i,k)=(1._r8-dum)*nric(i,k)
!               ! heating tendency 
!               tmp = xlf*dum*qrtot/(dz(i,k)*rho(i,k))
!               frzrdt(i,k)=frzrdt(i,k) + tmp
!
!               tlat(i,k)=tlat(i,k)+tmp
!               qstot=qstot+dum*qrtot
!               qrtot=(1._r8-dum)*qrtot
!               nstot=nstot+dum*nrtot
!               nrtot=(1._r8-dum)*nrtot
!               preci(i)=preci(i)+dum*(prect(i)-preci(i))
!            end if
!         end if
!
!         ! if rain/snow mix ratio is zero so should number concentration
!
!         if (qniic(i,k).lt.qsmall) then
!            qniic(i,k)=0._r8
!            nsic(i,k)=0._r8
!         end if
!
!         if (qric(i,k).lt.qsmall) then
!            qric(i,k)=0._r8
!            nric(i,k)=0._r8
!         end if
!
!         !if ( it .eq. 1 ) then
!         !    drout2(i,k) = drout2(i,k) + qric(i,k)
!         !endif
!         ! make sure number concentration is a positive number to avoid 
!         ! taking root of negative
!
!         nric(i,k)=max(nric(i,k),0._r8)
!         nsic(i,k)=max(nsic(i,k),0._r8)
!
!         !.......................................................................
!         ! get size distribution parameters for fallspeed calculations
!         !......................................................................
!         ! rain
!
!         if (qric(i,k).ge.qsmall) then
!             tmp_array(1) = (pi*rhow*nric(i,k)/qric(i,k))
!             tmp_array(2) = 1._r8/3._r8
!             call athread_spawn(slave_pow_parallel, tmp_array)
!             call athread_join()
!             tmp1 = tmp_array(1)
!            lamr(k) = tmp1
!
!            !lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
!            n0r(k) = nric(i,k)*lamr(k)
!
!            ! check for slope
!            ! change lammax and lammin for rain and snow
!            ! adjust vars
!
!            if (lamr(k).lt.lamminr) then
!
!               lamr(k) = lamminr
!
!               n0r(k) = lamr(k)*lamr(k)*lamr(k)*lamr(k)*qric(i,k)/(pi*rhow)
!               !n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
!               nric(i,k) = n0r(k)/lamr(k)
!            else if (lamr(k).gt.lammaxr) then
!               lamr(k) = lammaxr
!               n0r(k) = lamr(k)*lamr(k)*lamr(k)*lamr(k)*qric(i,k)/(pi*rhow)
!               !n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
!               nric(i,k) = n0r(k)/lamr(k)
!            end if
!
!
!            ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
!
!             tmp_array(1) = lamr(k)
!             tmp_array(2) = br
!             call athread_spawn(slave_pow_parallel, tmp_array)
!             call athread_join()
!             tmp1 = tmp_array(1)
!            unr(k) = min(arn(i,k)*cons4/tmp1,9.1_r8*rhof(i,k))
!            umr(k) = min(arn(i,k)*cons5/(6._r8*tmp1),9.1_r8*rhof(i,k))
!
!            !unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
!            !umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(i,k))
!
!         else
!            lamr(k) = 0._r8
!            n0r(k) = 0._r8
!            umr(k)=0._r8
!            unr(k)=0._r8
!         end if
!
!         !calculate mean size of combined rain and snow
!
!         if (lamr(k).gt.0._r8) then
!            Artmp = n0r(k) * pi / (2._r8 * lamr(k)*lamr(k)*lamr(k))
!            !Artmp = n0r(k) * pi / (2._r8 * lamr(k)**3._r8)
!         else 
!            Artmp = 0._r8
!         endif
!
!         if (lamc(k).gt.0._r8) then
!            Actmp = cdist1(k) * pi * gamma(pgam(k)+3._r8)/(4._r8 * lamc(k)*lamc(k))
!            !Actmp = cdist1(k) * pi * gamma(pgam(k)+3._r8)/(4._r8 * lamc(k)**2._r8)
!         else 
!            Actmp = 0._r8
!         endif
!
!         if (Actmp.gt.0_r8.or.Artmp.gt.0) then
!            rercld(i,k)=rercld(i,k) + 3._r8 *(qric(i,k) + qcic(i,k)) / (4._r8 * rhow * (Actmp + Artmp))
!            arcld(i,k)=arcld(i,k)+1._r8
!         endif
!
!         !......................................................................
!         ! snow
!
!         if (qniic(i,k).ge.qsmall) then
!             tmp_array(1) =(cons6*cs*nsic(i,k)/qniic(i,k))
!             tmp_array(2) = 1._r8/ds
!             call athread_spawn(slave_pow_parallel, tmp_array)
!             call athread_join()
!             tmp1 = tmp_array(1)
!             lams(k) = tmp1
!
!            !lams(k) = (cons6*cs*nsic(i,k)/ &
!            !     qniic(i,k))**(1._r8/ds)
!            n0s(k) = nsic(i,k)*lams(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lams(k).lt.lammins) then
!               lams(k) = lammins
!
!               tmp_array(1) = lams(k)
!               tmp_array(2) = ds+1._r8
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               n0s(k) = tmp1*qniic(i,k)/(cs*cons6)
!
!               !n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
!               nsic(i,k) = n0s(k)/lams(k)
!
!            else if (lams(k).gt.lammaxs) then
!               lams(k) = lammaxs
!
!               tmp_array(1) = lams(k)
!               tmp_array(2) = ds+1._r8
!               call athread_spawn(slave_pow_parallel, tmp_array)
!               call athread_join()
!               tmp1 = tmp_array(1)
!               n0s(k) = tmp1*qniic(i,k)/(cs*cons6)
!
!               !n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
!               nsic(i,k) = n0s(k)/lams(k)
!            end if
!
!            ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
!
!             tmp_array(1) = lams(k)
!             tmp_array(2) = bs
!             call athread_spawn(slave_pow_parallel, tmp_array)
!             call athread_join()
!             tmp1 = tmp_array(1)
!
!            ums(k) = min(asn(i,k)*cons8/(6._r8*tmp1),1.2_r8*rhof(i,k))
!            uns(k) = min(asn(i,k)*cons7/tmp1,1.2_r8*rhof(i,k))
!
!         else
!            lams(k) = 0._r8
!            n0s(k) = 0._r8
!            ums(k) = 0._r8
!            uns(k) = 0._r8
!         end if
!
!         !c........................................................................
!         ! sum over sub-step for average process rates
!
!         ! convert rain/snow q and N for output to history, note, 
!         ! output is for gridbox average
!
!         qrout(i,k)=qrout(i,k)+qric(i,k)*cldmax(i,k)
!         qsout(i,k)=qsout(i,k)+qniic(i,k)*cldmax(i,k)
!         nrout(i,k)=nrout(i,k)+nric(i,k)*rho(i,k)*cldmax(i,k)
!         nsout(i,k)=nsout(i,k)+nsic(i,k)*rho(i,k)*cldmax(i,k)
!
!         tlat1(i,k)=tlat1(i,k)+tlat(i,k)
!         qvlat1(i,k)=qvlat1(i,k)+qvlat(i,k)
!         qctend1(i,k)=qctend1(i,k)+qctend(i,k)
!         qitend1(i,k)=qitend1(i,k)+qitend(i,k)
!         nctend1(i,k)=nctend1(i,k)+nctend(i,k)
!         nitend1(i,k)=nitend1(i,k)+nitend(i,k)
!
!         t(i,k)=t(i,k)+tlat(i,k)*deltat/cpp
!         q(i,k)=q(i,k)+qvlat(i,k)*deltat
!         qc(i,k)=qc(i,k)+qctend(i,k)*deltat
!         qi(i,k)=qi(i,k)+qitend(i,k)*deltat
!         nc(i,k)=nc(i,k)+nctend(i,k)*deltat
!         ni(i,k)=ni(i,k)+nitend(i,k)*deltat
!
!         rainrt1(i,k)=rainrt1(i,k)+rainrt(i,k)
!
!         !divide rain radius over substeps for average
!         if (arcld(i,k) .gt. 0._r8) then
!            rercld(i,k)=rercld(i,k)/arcld(i,k)
!         end if
!
!        !if ( it .eq. 1  ) then
!        !    qrout2(i,k) = qrout2(i,k)+t(i,k)
!        !    qsout2(i,k) = qsout2(i,k)+q(i,k)
!        !    nrout2(i,k) = nrout2(i,k)+qc(i,k)
!        !    nsout2(i,k) = nsout2(i,k)+qi(i,k)
!        !    drout2(i,k) = drout2(i,k)+nc(i,k)
!        !    dsout2(i,k) = dsout2(i,k)+ni(i,k)
!        !endif
!
!
!         !calculate precip fluxes and adding them to summing sub-stepping variables
!         !! flux is zero at top interface
!         rflx(i,1)=0.0_r8
!         sflx(i,1)=0.0_r8
!
!         !! calculating the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
!         rflx(i,k+1)=qrout(i,k)*rho(i,k)*umr(k)
!         sflx(i,k+1)=qsout(i,k)*rho(i,k)*ums(k)
!
!         !! add to summing sub-stepping variable
!         rflx1(i,k+1)=rflx1(i,k+1)+rflx(i,k+1)
!         sflx1(i,k+1)=sflx1(i,k+1)+sflx(i,k+1)
!
!         !c........................................................................
!
!      end do ! k loop
!
!      prect1(i)=prect1(i)+prect(i)
!      preci1(i)=preci1(i)+preci(i)
!
!   end do ! it loop, sub-step
!
!   do k = top_lev, pver
!      rate1ord_cw2pr_st(i,k) = qcsinksum_rate1ord(k)/max(qcsum_rate1ord(k),1.0e-30_r8) 
!   end do
!
!300 continue  ! continue if no cloud water
!end do ! i loop
!
!! convert dt from sub-step back to full time step
!deltat=deltat*real(iter)
!
!
!!do i=1,ncol
!!do k=1,pver
!!!qrout2(i,k) = tlat1(i,k)
!!!qsout2(i,k) = qvlat1(i,k)
!!!nrout2(i,k) = qctend1(i,k)
!!!nsout2(i,k) = qitend1(i,k)
!!!drout2(i,k) = nctend1(i,k)
!!!dsout2(i,k) = nitend1(i,k)
!!qrout2(i,k) = prect1(i)
!!qsout2(i,k) = preci1(i)
!!nrout2(i,k) = sflx1(i,k+1)
!!nsout2(i,k) = rflx1(i,k+1)
!!drout2(i,k) = rainrt1(i,k)
!!dsout2(i,k) = nsout(i,k)
!!end do
!!end do
!
!!c.............................................................................
!
!do i=1,ncol
!
!   ! skip all calculations if no cloud water
!   if (ltrue(i).eq.0) then
!
!      do k=1,top_lev-1
!         ! assign zero values for effective radius above 1 mbar
!         effc(i,k)=0._r8
!         effi(i,k)=0._r8
!         effc_fn(i,k)=0._r8
!         lamcrad(i,k)=0._r8
!         pgamrad(i,k)=0._r8
!         deffi(i,k)=0._r8
!      end do
!
!      do k=top_lev,pver
!         ! assign default values for effective radius
!         effc(i,k)=10._r8
!         effi(i,k)=25._r8
!         effc_fn(i,k)=10._r8
!         lamcrad(i,k)=0._r8
!         pgamrad(i,k)=0._r8
!         deffi(i,k)=0._r8
!      end do
!      goto 500
!   end if
!
!   ! initialize nstep for sedimentation sub-steps
!   nstep = 1
!
!   ! divide precip rate by number of sub-steps to get average over time step
!
!   prect(i)=prect1(i)/real(iter)
!   preci(i)=preci1(i)/real(iter)
!
!   do k=top_lev,pver
!
!      ! assign variables back to start-of-timestep values before updating after sub-steps 
!
!      t(i,k)=t1(i,k)
!      q(i,k)=q1(i,k)
!      qc(i,k)=qc1(i,k)
!      qi(i,k)=qi1(i,k)
!      nc(i,k)=nc1(i,k)
!      ni(i,k)=ni1(i,k)
!
!      ! divide microphysical tendencies by number of sub-steps to get average over time step
!
!      tlat(i,k)=tlat1(i,k)/real(iter)
!      qvlat(i,k)=qvlat1(i,k)/real(iter)
!      qctend(i,k)=qctend1(i,k)/real(iter)
!      qitend(i,k)=qitend1(i,k)/real(iter)
!      nctend(i,k)=nctend1(i,k)/real(iter)
!      nitend(i,k)=nitend1(i,k)/real(iter)
!
!      rainrt(i,k)=rainrt1(i,k)/real(iter)
!
!      ! divide by number of sub-steps to find final values
!      rflx(i,k+1)=rflx1(i,k+1)/real(iter)
!      sflx(i,k+1)=sflx1(i,k+1)/real(iter)
!
!      ! divide output precip q and N by number of sub-steps to get average over time step
!
!      qrout(i,k)=qrout(i,k)/real(iter)
!      qsout(i,k)=qsout(i,k)/real(iter)
!      nrout(i,k)=nrout(i,k)/real(iter)
!      nsout(i,k)=nsout(i,k)/real(iter)
!
!      ! divide trop_mozart variables by number of sub-steps to get average over time step 
!
!      nevapr(i,k) = nevapr(i,k)/real(iter)
!      nevapr2(i,k) = nevapr2(i,k)/real(iter)
!      evapsnow(i,k) = evapsnow(i,k)/real(iter)
!      prain(i,k) = prain(i,k)/real(iter)
!      prodsnow(i,k) = prodsnow(i,k)/real(iter)
!      cmeout(i,k) = cmeout(i,k)/real(iter)
!
!      cmeiout(i,k) = cmeiout(i,k)/real(iter)
!      meltsdt(i,k) = meltsdt(i,k)/real(iter)
!      frzrdt (i,k) = frzrdt (i,k)/real(iter)
!
!
!      ! microphysics output
!      prao(i,k)=prao(i,k)/real(iter)
!      prco(i,k)=prco(i,k)/real(iter)
!      mnuccco(i,k)=mnuccco(i,k)/real(iter)
!      mnuccto(i,k)=mnuccto(i,k)/real(iter)
!      msacwio(i,k)=msacwio(i,k)/real(iter)
!      psacwso(i,k)=psacwso(i,k)/real(iter)
!      bergso(i,k)=bergso(i,k)/real(iter)
!      bergo(i,k)=bergo(i,k)/real(iter)
!      prcio(i,k)=prcio(i,k)/real(iter)
!      praio(i,k)=praio(i,k)/real(iter)
!
!      mnuccro(i,k)=mnuccro(i,k)/real(iter)
!      pracso (i,k)=pracso (i,k)/real(iter)
!
!      mnuccdo(i,k)=mnuccdo(i,k)/real(iter)
!
!      ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
!      nevapr(i,k) = nevapr(i,k) + evapsnow(i,k)
!      prer_evap(i,k) = nevapr2(i,k)
!      prain(i,k) = prain(i,k) + prodsnow(i,k)
!
!      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      ! calculate sedimentation for cloud water and ice
!      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      ! update in-cloud cloud mixing ratio and number concentration 
!      ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
!      ! note: these are in-cloud values***, hence we divide by cloud fraction
!
!      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)
!      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)/icldm(i,k)
!      dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k),0._r8)
!      dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)/icldm(i,k),0._r8)
!
!      ! obtain new slope parameter to avoid possible singularity
!
!      if (dumi(i,k).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)
!
!         lami(k) = (cons1*ci* &
!              dumni(i,k)/dumi(i,k))**(1._r8/di)
!         lami(k)=max(lami(k),lammini)
!         lami(k)=min(lami(k),lammaxi)
!      else
!         lami(k)=0._r8
!      end if
!
!      if (dumc(i,k).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)
!         ! add lower limit to in-cloud number concentration
!         dumnc(i,k)=max(dumnc(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3 
!         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)*pgam(k))-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
!              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         lammax = (pgam(k)+1._r8)/2.e-6_r8
!         lamc(k)=max(lamc(k),lammin)
!         lamc(k)=min(lamc(k),lammax)
!      else
!         lamc(k)=0._r8
!      end if
!
!      ! calculate number and mass weighted fall velocity for droplets
!      ! include effects of sub-grid distribution of cloud water
!
!
!      if (dumc(i,k).ge.qsmall) then
!         unc = acn(i,k)*gamma(1._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+1._r8))
!         umc = acn(i,k)*gamma(4._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+4._r8))
!         ! fallspeed for output
!         vtrmc(i,k)=umc
!      else
!         umc = 0._r8
!         unc = 0._r8
!      end if
!
!      ! calculate number and mass weighted fall velocity for cloud ice
!
!      if (dumi(i,k).ge.qsmall) then
!         uni =  ain(i,k)*cons16/lami(k)**bi
!         umi = ain(i,k)*cons17/(6._r8*lami(k)**bi)
!         uni=min(uni,1.2_r8*rhof(i,k))
!         umi=min(umi,1.2_r8*rhof(i,k))
!
!         ! fallspeed
!         vtrmi(i,k)=umi
!      else
!         umi = 0._r8
!         uni = 0._r8
!      end if
!
!      fi(k) = g*rho(i,k)*umi
!      fni(k) = g*rho(i,k)*uni
!      fc(k) = g*rho(i,k)*umc
!      fnc(k) = g*rho(i,k)*unc
!
!      ! calculate number of split time steps to ensure courant stability criteria
!      ! for sedimentation calculations
!
!      rgvm = max(fi(k),fc(k),fni(k),fnc(k))
!      nstep = max(int(rgvm*deltat/pdel(i,k)+1._r8),nstep)
!
!      ! redefine dummy variables - sedimentation is calculated over grid-scale
!      ! quantities to ensure conservation
!
!      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
!      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
!      dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
!      dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)
!
!      if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
!      if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
!
!   end do       !!! vertical loop
!
!   !do k=1,pver
!   !     qrout2(i,k) = nstep 
!   !     qsout2(i,k) = nstep 
!   !     nrout2(i,k) = nstep 
!   !     nsout2(i,k) = nstep 
!   !     drout2(i,k) = nstep 
!   !     dsout2(i,k) = nstep 
!   !end do
!
!   do n = 1,nstep  !! loop over sub-time step to ensure stability
!
!      do k = top_lev,pver
!         if (do_cldice) then
!            falouti(k) = fi(k)*dumi(i,k)
!            faloutni(k) = fni(k)*dumni(i,k)
!         else
!            falouti(k)  = 0._r8
!            faloutni(k) = 0._r8
!         end if
!
!         faloutc(k) = fc(k)*dumc(i,k)
!         faloutnc(k) = fnc(k)*dumnc(i,k)
!      end do
!
!      ! top of model
!
!      k = top_lev
!      faltndi = falouti(k)/pdel(i,k)
!      faltndni = faloutni(k)/pdel(i,k)
!      faltndc = faloutc(k)/pdel(i,k)
!      faltndnc = faloutnc(k)/pdel(i,k)
!
!      ! add fallout terms to microphysical tendencies
!
!      qitend(i,k) = qitend(i,k)-faltndi/nstep
!      nitend(i,k) = nitend(i,k)-faltndni/nstep
!      qctend(i,k) = qctend(i,k)-faltndc/nstep
!      nctend(i,k) = nctend(i,k)-faltndnc/nstep
!
!      ! sedimentation tendencies for output
!      qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
!      qisedten(i,k)=qisedten(i,k)-faltndi/nstep
!
!      dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
!      dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
!      dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
!      dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep
!
!      do k = top_lev+1,pver
!
!         ! for cloud liquid and ice, if cloud fraction increases with height
!         ! then add flux from above to both vapor and cloud water of current level
!         ! this means that flux entering clear portion of cell from above evaporates
!         ! instantly
!
!         dum=lcldm(i,k)/lcldm(i,k-1)
!         dum=min(dum,1._r8)
!         dum1=icldm(i,k)/icldm(i,k-1)
!         dum1=min(dum1,1._r8)
!
!         faltndqie=(falouti(k)-falouti(k-1))/pdel(i,k)
!         faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(i,k)
!         faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(i,k)
!         faltndqce=(faloutc(k)-faloutc(k-1))/pdel(i,k)
!         faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(i,k)
!         faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(i,k)
!
!         ! add fallout terms to eulerian tendencies
!
!         qitend(i,k) = qitend(i,k)-faltndi/nstep
!         nitend(i,k) = nitend(i,k)-faltndni/nstep
!         qctend(i,k) = qctend(i,k)-faltndc/nstep
!         nctend(i,k) = nctend(i,k)-faltndnc/nstep
!
!         !if ( nstep .eq. 1 ) then
!         !    qrout2(i,k) = qrout2(i,k) + qcsedten(i,k)
!         !    qsout2(i,k) = qsout2(i,k) + faltndc/nstep
!         !    nrout2(i,k) = nrout2(i,k) + faltndc
!         !    nsout2(i,k) = nsout2(i,k) + nstep
!         !endif
!         ! sedimentation tendencies for output
!         qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
!         qisedten(i,k)=qisedten(i,k)-faltndi/nstep
!
!         !if ( nstep .eq. 1 ) then
!         !    drout2(i,k) = drout2(i,k) + qcsedten(i,k)          
!         !    dsout2(i,k) = dsout2(i,k) + qisedten(i,k)
!         !endif
!
!         ! add terms to to evap/sub of cloud water
!
!         qvlat(i,k)=qvlat(i,k)-(faltndqie-faltndi)/nstep
!         ! for output
!         qisevap(i,k)=qisevap(i,k)-(faltndqie-faltndi)/nstep
!         qvlat(i,k)=qvlat(i,k)-(faltndqce-faltndc)/nstep
!
!         !if ( n .eq. 2 ) then
!         !    qrout2(i,k) = qrout2(i,k) + qcsevap(i,k) 
!         !    qsout2(i,k) = qsout2(i,k) + faltndqce
!         !    nrout2(i,k) = nrout2(i,k) + faltndc
!         !    nsout2(i,k) = nsout2(i,k) + nstep
!         !    drout2(i,k) = drout2(i,k) + (faltndqce-faltndc)/nstep          
!         !    dsout2(i,k) = dsout2(i,k) + qcsevap(i,k)-(faltndqce-faltndc)/nstep
!         !endif
!
!         ! for output
!         qcsevap(i,k)=qcsevap(i,k)-(faltndqce-faltndc)/nstep
!
!         tlat(i,k)=tlat(i,k)+(faltndqie-faltndi)*xxls/nstep
!         tlat(i,k)=tlat(i,k)+(faltndqce-faltndc)*xxlv/nstep
!
!         dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
!         dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
!         dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
!         dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep
!
!         Fni(K)=MAX(Fni(K)/pdel(i,K),Fni(K-1)/pdel(i,K-1))*pdel(i,K)
!         FI(K)=MAX(FI(K)/pdel(i,K),FI(K-1)/pdel(i,K-1))*pdel(i,K)
!         fnc(k)=max(fnc(k)/pdel(i,k),fnc(k-1)/pdel(i,k-1))*pdel(i,k)
!         Fc(K)=MAX(Fc(K)/pdel(i,K),Fc(K-1)/pdel(i,K-1))*pdel(i,K)
!
!      end do   !! k loop
!
!      ! units below are m/s
!      ! cloud water/ice sedimentation flux at surface 
!      ! is added to precip flux at surface to get total precip (cloud + precip water)
!      ! rate
!
!      prect(i) = prect(i)+(faloutc(pver)+falouti(pver))/g/nstep/1000._r8  
!      preci(i) = preci(i)+(falouti(pver))/g/nstep/1000._r8
!
!   end do   !! nstep loop
!
!           !do k=1,pver
!           !     qrout2(i,k) = qrout2(i,k) + qitend(i,k) 
!           !     qsout2(i,k) = qsout2(i,k) + nitend(i,k) 
!           !     nrout2(i,k) = nrout2(i,k) + qctend(i,k)
!           !     nsout2(i,k) = nsout2(i,k) + nctend(i,k) 
!           !     drout2(i,k) = drout2(i,k) + qcsedten(i,k)          
!           !     dsout2(i,k) = dsout2(i,k) + qisedten(i,k)
!           !end do        
!   ! end sedimentation
!   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   ! get new update for variables that includes sedimentation tendency
!   ! note : here dum variables are grid-average, NOT in-cloud
!
!   do k=top_lev,pver
!
!      dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)
!      dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)
!      dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)
!      dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)
!
!      if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
!      if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
!
!      ! calculate instantaneous processes (melting, homogeneous freezing)
!      if (do_cldice) then
!
!         if (t(i,k)+tlat(i,k)/cpp*deltat > tmelt) then
!            if (dumi(i,k) > 0._r8) then
!
!               ! limit so that melting does not push temperature below freezing
!               dum = -dumi(i,k)*xlf/cpp
!               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.tmelt) then
!                  dum = (t(i,k)+tlat(i,k)/cpp*deltat-tmelt)*cpp/xlf
!                  dum = dum/dumi(i,k)*xlf/cpp 
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qctend(i,k)=qctend(i,k)+dum*dumi(i,k)/deltat
!
!               ! for output
!               melto(i,k)=dum*dumi(i,k)/deltat
!
!               ! assume melting ice produces droplet
!               ! mean volume radius of 8 micron
!
!               nctend(i,k)=nctend(i,k)+3._r8*dum*dumi(i,k)/deltat/ &
!                    (4._r8*pi*5.12e-16_r8*rhow)
!
!               qitend(i,k)=((1._r8-dum)*dumi(i,k)-qi(i,k))/deltat
!               nitend(i,k)=((1._r8-dum)*dumni(i,k)-ni(i,k))/deltat
!               tlat(i,k)=tlat(i,k)-xlf*dum*dumi(i,k)/deltat
!            end if
!         end if
!
!         ! homogeneously freeze droplets at -40 C
!
!         if (t(i,k)+tlat(i,k)/cpp*deltat < 233.15_r8) then
!            if (dumc(i,k) > 0._r8) then
!
!               ! limit so that freezing does not push temperature above threshold
!               dum = dumc(i,k)*xlf/cpp
!               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.233.15_r8) then
!                  dum = -(t(i,k)+tlat(i,k)/cpp*deltat-233.15_r8)*cpp/xlf
!                  dum = dum/dumc(i,k)*xlf/cpp
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qitend(i,k)=qitend(i,k)+dum*dumc(i,k)/deltat
!               ! for output
!               homoo(i,k)=dum*dumc(i,k)/deltat
!
!               ! assume 25 micron mean volume radius of homogeneously frozen droplets
!               ! consistent with size of detrained ice in stratiform.F90
!               nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*1.563e-14_r8* &
!                    500._r8)/deltat
!               qctend(i,k)=((1._r8-dum)*dumc(i,k)-qc(i,k))/deltat
!               nctend(i,k)=((1._r8-dum)*dumnc(i,k)-nc(i,k))/deltat
!               tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)/deltat
!            end if
!         end if
!
!         ! remove any excess over-saturation, which is possible due to non-linearity when adding 
!         ! together all microphysical processes
!         ! follow code similar to old CAM scheme
!
!         qtmp=q(i,k)+qvlat(i,k)*deltat
!         ttmp=t(i,k)+tlat(i,k)/cpp*deltat
!
!         !esn = svp_water(ttmp)  ! use rhw to allow ice supersaturation
!
!         tmp1 = svp_tboil/ttmp
!         call athread_spawn(slave_log10_parallel, tmp1)
!         call athread_join()
!         esn = 10.d0**(-7.90298d0*(svp_tboil/ttmp-1.d0)+ &
!            5.02808d0*tmp1- &
!            1.3816d-7*(10.d0**(11.344d0*(1.d0-ttmp/svp_tboil))-1.d0)+ &
!            8.1328d-3*(10.d0**(-3.49149d0*(svp_tboil/ttmp-1.d0))-1.d0)+ &
!            log10(1013.246d0))*100.d0
!
!         qsn = svp_to_qsat(esn, p(i,k))
!
!         if (qtmp > qsn .and. qsn > 0) then
!            ! expression below is approximate since there may be ice deposition
!            dum = (qtmp-qsn)/(1._r8+cons27*qsn/(cpp*rv*ttmp*ttmp))/deltat
!            ! add to output cme
!            cmeout(i,k) = cmeout(i,k)+dum
!            ! now add to tendencies, partition between liquid and ice based on temperature
!            if (ttmp > 268.15_r8) then
!               dum1=0.0_r8
!               ! now add to tendencies, partition between liquid and ice based on te
!            else if (ttmp < 238.15_r8) then
!               dum1=1.0_r8
!            else
!               dum1=(268.15_r8-ttmp)/30._r8
!            end if
!
!            dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))*(xxls*dum1+xxlv*(1._r8-dum1)) &
!                 *qsn/(cpp*rv*ttmp*ttmp))/deltat
!            qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
!            ! for output
!            qcreso(i,k)=dum*(1._r8-dum1)
!            qitend(i,k)=qitend(i,k)+dum*dum1
!            qireso(i,k)=dum*dum1
!            qvlat(i,k)=qvlat(i,k)-dum
!            ! for output
!            qvres(i,k)=-dum
!            tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
!         end if
!      end if
!
!      !...............................................................................
!      ! calculate effective radius for pass to radiation code
!      ! if no cloud water, default value is 10 micron for droplets,
!      ! 25 micron for cloud ice
!
!      ! update cloud variables after instantaneous processes to get effective radius
!      ! variables are in-cloud to calculate size dist parameters
!
!      dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/lcldm(i,k)
!      dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/icldm(i,k)
!      dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/lcldm(i,k)
!      dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/icldm(i,k)
!
!      ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
!
!      dumc(i,k)=min(dumc(i,k),5.e-3_r8)
!      dumi(i,k)=min(dumi(i,k),5.e-3_r8)
!
!      !...................
!      ! cloud ice effective radius
!
!      if (dumi(i,k).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)
!         lami(k) = (cons1*ci*dumni(i,k)/dumi(i,k))**(1._r8/di)
!
!         if (lami(k).lt.lammini) then
!            lami(k) = lammini
!            n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*cons1)
!            niic(i,k) = n0i(k)/lami(k)
!            ! adjust number conc if needed to keep mean size in reasonable range
!            if (do_cldice) nitend(i,k)=(niic(i,k)*icldm(i,k)-ni(i,k))/deltat
!
!         else if (lami(k).gt.lammaxi) then
!            lami(k) = lammaxi
!            n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*cons1)
!            niic(i,k) = n0i(k)/lami(k)
!            ! adjust number conc if needed to keep mean size in reasonable range
!            if (do_cldice) nitend(i,k)=(niic(i,k)*icldm(i,k)-ni(i,k))/deltat
!         end if
!         effi(i,k) = 1.5_r8/lami(k)*1.e6_r8
!
!      else
!         effi(i,k) = 25._r8
!      end if
!
!      ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
!      ! radius has already been determined from the size distribution.
!      if (.not. do_cldice) then
!         effi(i,k) = re_ice(i,k) * 1e6_r8      ! m -> um
!      end if
!
!      !...................
!      ! cloud droplet effective radius
!
!      if (dumc(i,k).ge.qsmall) then
!
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! set tendency to ensure minimum droplet concentration
!         ! after update by microphysics, except when lambda exceeds bounds on mean drop
!         ! size or if there is no cloud water
!         if (dumnc(i,k).lt.cdnl/rho(i,k)) then   
!            nctend(i,k)=(cdnl/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat   
!         end if
!         dumnc(i,k)=max(dumnc(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)*pgam(k))-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
!              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         ! Multiply by omsm to fit within RRTMG's table.
!         lammax = (pgam(k)+1._r8)*omsm/2.e-6_r8
!         if (lamc(k).lt.lammin) then
!            lamc(k) = lammin
!            ncic(i,k) = 6._r8*lamc(k)*lamc(k)*lamc(k)*dumc(i,k)* &
!                 gamma(pgam(k)+1._r8)/ &
!                 (pi*rhow*gamma(pgam(k)+4._r8))
!            ! adjust number conc if needed to keep mean size in reasonable range
!            nctend(i,k)=(ncic(i,k)*lcldm(i,k)-nc(i,k))/deltat
!
!         else if (lamc(k).gt.lammax) then
!            lamc(k) = lammax
!            ncic(i,k) = 6._r8*lamc(k)*lamc(k)*lamc(k)*dumc(i,k)* &
!                 gamma(pgam(k)+1._r8)/ &
!                 (pi*rhow*gamma(pgam(k)+4._r8))
!            ! adjust number conc if needed to keep mean size in reasonable range
!            nctend(i,k)=(ncic(i,k)*lcldm(i,k)-nc(i,k))/deltat
!         end if
!
!         effc(i,k) = &
!              gamma(pgam(k)+4._r8)/ &
!              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
!         !assign output fields for shape here
!         lamcrad(i,k)=lamc(k)
!         pgamrad(i,k)=pgam(k)
!
!      else
!         effc(i,k) = 10._r8
!         lamcrad(i,k)=0._r8
!         pgamrad(i,k)=0._r8
!      end if
!
!      ! ice effective diameter for david mitchell's optics
!      if (do_cldice) then
!         deffi(i,k)=effi(i,k)*rhoi/917._r8*2._r8
!      else
!         deffi(i,k)=effi(i,k) * 2._r8
!      end if
!
!
!!!! recalculate effective radius for constant number, in order to separate
!      ! first and second indirect effects
!      ! assume constant number of 10^8 kg-1
!
!      dumnc(i,k)=1.e8_r8
!
!      if (dumc(i,k).ge.qsmall) then
!         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)*pgam(k))-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
!              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         lammax = (pgam(k)+1._r8)/2.e-6_r8
!         if (lamc(k).lt.lammin) then
!            lamc(k) = lammin
!         else if (lamc(k).gt.lammax) then
!            lamc(k) = lammax
!         end if
!         effc_fn(i,k) = &
!              gamma(pgam(k)+4._r8)/ &
!              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
!
!      else
!         effc_fn(i,k) = 10._r8
!      end if
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!
!
!   end do ! vertical k loop
!
!           !do k=1,pver
!           !     qrout2(i,k) = qcreso(i,k) 
!           !     qsout2(i,k) = qitend(i,k) 
!           !     nrout2(i,k) = qireso(i,k)
!           !     nsout2(i,k) = qvlat(i,k) 
!           !     drout2(i,k) = deffi(i,k)          
!           !     dsout2(i,k) = effc_fn(i,k)
!           !end do        
!
!500 continue
!
!   do k=top_lev,pver
!      ! if updated q (after microphysics) is zero, then ensure updated n is also zero
!
!      if (qc(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-nc(i,k)/deltat
!      if (do_cldice .and. qi(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-ni(i,k)/deltat
!   end do
!
!end do ! i loop
!
!! add snow ouptut
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
!         dsout(i,k)=3._r8*rhosn/917._r8*(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
!      endif
!   end do
!end do
!
!!calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual
!do i = 1,ncol
!   do k=top_lev,pver
!      !! RAIN
!      if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
!         reff_rain(i,k)=1.5_r8*(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)*1.e6_r8
!      endif
!      !! SNOW
!      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
!         reff_snow(i,k)=1.5_r8*(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)*1.e6_r8
!      end if
!   end do
!end do
!
!! analytic radar reflectivity
!! formulas from Matthew Shupe, NOAA/CERES
!! *****note: radar reflectivity is local (in-precip average)
!! units of mm^6/m^3
!
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qc(i,k)+qctend(i,k)*deltat.ge.qsmall) then
!         dum=((qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)*1000._r8)*((qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)*1000._r8) &
!              /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)/1.e6_r8)*lcldm(i,k)/cldmax(i,k)
!      else
!         dum=0._r8
!      end if
!      if (qi(i,k)+qitend(i,k)*deltat.ge.qsmall) then
!         dum1=((qi(i,k)+qitend(i,k)*deltat)*rho(i,k)/icldm(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(i,k)/cldmax(i,k)
!      else 
!         dum1=0._r8
!      end if
!
!      if (qsout(i,k).ge.qsmall) then
!         dum1=dum1+(qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
!      end if
!
!      refl(i,k)=dum+dum1
!
!      ! add rain rate, but for 37 GHz formulation instead of 94 GHz
!      ! formula approximated from data of Matrasov (2007)
!      ! rainrt is the rain rate in mm/hr
!      ! reflectivity (dum) is in DBz
!
!      if (rainrt(i,k).ge.0.001_r8) then
!         dum=log10(rainrt(i,k)**6._r8)+16._r8
!
!         ! convert from DBz to mm^6/m^3
!
!         dum = 10._r8**(dum/10._r8)
!      else
!         ! don't include rain rate in R calculation for values less than 0.001 mm/hr
!         dum=0._r8
!      end if
!
!      ! add to refl
!
!      refl(i,k)=refl(i,k)+dum
!
!      !output reflectivity in Z.
!      areflz(i,k)=refl(i,k)
!
!      ! convert back to DBz 
!
!      if (refl(i,k).gt.minrefl) then 
!         refl(i,k)=10._r8*log10(refl(i,k))
!      else
!         refl(i,k)=-9999._r8
!      end if
!
!      !set averaging flag
!      if (refl(i,k).gt.mindbz) then 
!         arefl(i,k)=refl(i,k)
!         frefl(i,k)=1.0_r8  
!      else
!         arefl(i,k)=0._r8
!         areflz(i,k)=0._r8
!         frefl(i,k)=0._r8
!      end if
!
!      ! bound cloudsat reflectivity
!
!      csrfl(i,k)=min(csmax,refl(i,k))
!
!      !set averaging flag
!      if (csrfl(i,k).gt.csmin) then 
!         acsrfl(i,k)=refl(i,k)
!         fcsrfl(i,k)=1.0_r8  
!      else
!         acsrfl(i,k)=0._r8
!         fcsrfl(i,k)=0._r8
!      end if
!
!   end do
!end do
!
!
!! averaging for snow and rain number and diameter
!
!qrout2(:,:)=0._r8
!qsout2(:,:)=0._r8
!nrout2(:,:)=0._r8
!nsout2(:,:)=0._r8
!drout2(:,:)=0._r8
!dsout2(:,:)=0._r8
!freqs(:,:)=0._r8
!freqr(:,:)=0._r8
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
!         qrout2(i,k)=qrout(i,k)
!         nrout2(i,k)=nrout(i,k)
!         drout2(i,k)=(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)
!         freqr(i,k)=1._r8
!      endif
!      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
!         qsout2(i,k)=qsout(i,k)
!         nsout2(i,k)=nsout(i,k)
!         dsout2(i,k)=(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
!         freqs(i,k)=1._r8
!      endif
!   end do
!end do
!
!! output activated liquid and ice (convert from #/kg -> #/m3)
!do i = 1,ncol
!   do k=top_lev,pver
!      ncai(i,k)=dum2i(i,k)*rho(i,k)
!      ncal(i,k)=dum2l(i,k)*rho(i,k)
!   end do
!end do
!
!
!!redefine fice here....
!nfice(:,:)=0._r8
!do k=top_lev,pver
!   do i=1,ncol
!      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
!      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
!      dumfice=qsout(i,k) + qrout(i,k) + dumc(i,k) + dumi(i,k)  
!
!      if (dumfice.gt.qsmall.and.(qsout(i,k)+dumi(i,k).gt.qsmall)) then
!         nfice(i,k)=(qsout(i,k) + dumi(i,k))/dumfice
!      endif
!
!      if (nfice(i,k).gt.1._r8) then
!         nfice(i,k)=1._r8
!      endif
!
!   enddo
!enddo


!if ( rank .eq. 0 ) then
!do i=1,pcols
!        !if ( abs(prect    (i) - prect_permute(i)) / max(abs(prect    (i)), abs(prect_permute(i))) > 1e-10_r8 ) then
!        !    print*, "err:", abs(prect    (i) - prect_permute(i)) / max(abs(prect    (i)), abs(prect_permute(i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", prect    (i), " prect_permute: ", prect_permute(i)
!        !endif
!        !if ( abs(preci    (i) - preci_permute(i)) / max(abs(preci    (i)), abs(preci_permute(i))) > 1e-10_r8 ) then
!        !    print*, "err:", abs(preci    (i) - preci_permute(i)) / max(abs(preci    (i)), abs(preci_permute(i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", preci    (i), " preci_permute: ", preci_permute(i)
!        !endif
!    do k=1,pver
!
!        !if ( abs(qrout2(i,k) - qrout2_permute(k,i)) / abs(qrout2(i,k)) > 1e-19_r8 ) then
!        !    print*, "err:", abs(qrout2(i,k) - qrout2_permute(k,i)) / abs(qrout2(i,k)), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "ncol: ", ncol, "pcols: ", pcols, "qrout2: ", qrout2(i,k), " qrout2_permute: ", qrout2_permute(k,i)
!        !endif
!        !if ( abs(qsout2(i,k) - qsout2_permute(k,i)) / abs(qsout2(i,k)) > 1e-19_r8 ) then
!        !    print*, "err:", abs(qsout2(i,k) - qsout2_permute(k,i)) / abs(qsout2(i,k)), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "ncol: ", ncol, "pcols: ", pcols, "qsout2: ", qsout2(i,k), " qsout2_permute: ", qsout2_permute(k,i)
!        !endif
!        !if ( abs(nrout2(i,k) - nrout2_permute(k,i)) / abs(nrout2(i,k)) > 1e-19_r8 ) then
!        !    print*, "err:", abs(nrout2(i,k) - nrout2_permute(k,i)) / abs(nrout2(i,k)), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "ncol: ", ncol, "pcols: ", pcols, "nrout2: ", nrout2(i,k), " nrout2_permute: ", nrout2_permute(k,i)
!        !endif
!        !if ( abs(nsout2(i,k) - nsout2_permute(k,i)) / abs(nsout2(i,k)) > 1e-19_r8 ) then
!        !    print*, "err:", abs(nsout2(i,k) - nsout2_permute(k,i)) / abs(nsout2(i,k)), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "ncol: ", ncol, "pcols: ", pcols, "nsout2: ", nsout2(i,k), " nsout2_permute: ", nsout2_permute(k,i)
!        !endif
!        !if ( abs(drout2(i,k) - drout2_permute(k,i)) / abs(drout2(i,k)) > 1e-19_r8 ) then
!        !    print*, "err:", abs(drout2(i,k) - drout2_permute(k,i)) / abs(drout2(i,k)), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "ncol: ", ncol, "pcols: ", pcols, "drout2: ", drout2(i,k), " drout2_permute: ", drout2_permute(k,i)
!        !endif
!        !if ( abs(dsout2(i,k) - dsout2_permute(k,i)) / abs(dsout2(i,k)) > 1e-19_r8 ) then
!        !    print*, "err:", abs(dsout2(i,k) - dsout2_permute(k,i)) / abs(dsout2(i,k)), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "dsout2: ", dsout2(i,k), " dsout2_permute: ", dsout2_permute(k,i)
!        !endif
!
!        if ( abs(qc(i,k) - qc_permute(k,i)) / max(abs(qc(i,k)), abs(qc_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qc(i,k) - qc_permute(k,i)) / max(abs(qc(i,k)), abs(qc_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "qc: ", qc(i,k), " qc_permute: ", qc_permute(k,i)
!        endif
!        if ( abs(qi(i,k) - qi_permute(k,i)) / max(abs(qi(i,k)), abs(qi_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qi(i,k) - qi_permute(k,i)) / max(abs(qi(i,k)), abs(qi_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "qi: ", qi(i,k), " qi_permute: ", qi_permute(k,i)
!        endif
!        if ( abs(nc(i,k) - nc_permute(k,i)) / max(abs(nc(i,k)), abs(nc_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nc(i,k) - nc_permute(k,i)) / max(abs(nc(i,k)), abs(nc_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "nc: ", nc(i,k), " nc_permute: ", nc_permute(k,i)
!        endif
!        if ( abs(ni(i,k) - ni_permute(k,i)) / max(abs(ni(i,k)), abs(ni_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(ni(i,k) - ni_permute(k,i)) / max(abs(ni(i,k)), abs(ni_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "ni: ", ni(i,k), " ni_permute: ", ni_permute(k,i)
!        endif
!
!        if ( abs(rate1ord_cw2pr_st(i,k) - rate1ord_cw2pr_st_permute(k,i)) / max(abs(rate1ord_cw2pr_st(i,k)), abs(rate1ord_cw2pr_st_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(rate1ord_cw2pr_st(i,k) - rate1ord_cw2pr_st_permute(k,i)) / max(abs(rate1ord_cw2pr_st(i,k)), abs(rate1ord_cw2pr_st_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, "rate1ord_cw2pr_st: ", rate1ord_cw2pr_st(i,k), " rate1ord_cw2pr_st_permute: ", rate1ord_cw2pr_st_permute(k,i)
!        endif
!
!        if ( abs(tlat     (i,k) - tlat_permute(k,i)) / max(abs(tlat     (i,k)), abs(tlat_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(tlat     (i,k) - tlat_permute(k,i)) / max(abs(tlat     (i,k)), abs(tlat_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", tlat     (i,k), " tlat_permute: ", tlat_permute(k,i)
!        endif
!        if ( abs(qvlat    (i,k) - qvlat_permute(k,i)) / max(abs(qvlat    (i,k)), abs(qvlat_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qvlat    (i,k) - qvlat_permute(k,i)) / max(abs(qvlat    (i,k)), abs(qvlat_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qvlat    (i,k), " qvlat_permute: ", qvlat_permute(k,i)
!        endif
!        if ( abs(qctend   (i,k) - qctend_permute(k,i)) / max(abs(qctend   (i,k)), abs(qctend_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qctend   (i,k) - qctend_permute(k,i)) / max(abs(qctend   (i,k)), abs(qctend_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qctend   (i,k), " qctend_permute: ", qctend_permute(k,i)
!        endif
!        if ( abs(qitend   (i,k) - qitend_permute(k,i)) / max(abs(qitend   (i,k)), abs(qitend_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qitend   (i,k) - qitend_permute(k,i)) / max(abs(qitend   (i,k)), abs(qitend_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qitend   (i,k), " qitend_permute: ", qitend_permute(k,i)
!        endif
!        if ( abs(nctend   (i,k) - nctend_permute(k,i)) / max(abs(nctend   (i,k)), abs(nctend_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nctend   (i,k) - nctend_permute(k,i)) / max(abs(nctend   (i,k)), abs(nctend_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", nctend   (i,k), " nctend_permute: ", nctend_permute(k,i)
!        endif
!        if ( abs(nitend   (i,k) - nitend_permute(k,i)) / max(abs(nitend   (i,k)), abs(nitend_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nitend   (i,k) - nitend_permute(k,i)) / max(abs(nitend   (i,k)), abs(nitend_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", nitend   (i,k), " nitend_permute: ", nitend_permute(k,i)
!        endif
!        if ( abs(effc     (i,k) - effc_permute(k,i)) / max(abs(effc     (i,k)), abs(effc_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(effc     (i,k) - effc_permute(k,i)) / max(abs(effc     (i,k)), abs(effc_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", effc     (i,k), " effc_permute: ", effc_permute(k,i)
!        endif
!        if ( abs(effc_fn  (i,k) - effc_fn_permute(k,i)) / max(abs(effc_fn  (i,k)), abs(effc_fn_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(effc_fn  (i,k) - effc_fn_permute(k,i)) / max(abs(effc_fn  (i,k)), abs(effc_fn_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", effc_fn  (i,k), " effc_fn_permute: ", effc_fn_permute(k,i)
!        endif
!        if ( abs(effi     (i,k) - effi_permute(k,i)) / max(abs(effi     (i,k)), abs(effi_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(effi     (i,k) - effi_permute(k,i)) / max(abs(effi     (i,k)), abs(effi_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", effi     (i,k), " effi_permute: ", effi_permute(k,i)
!        endif
!        if ( abs(nevapr   (i,k) - nevapr_permute(k,i)) / max(abs(nevapr   (i,k)), abs(nevapr_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nevapr   (i,k) - nevapr_permute(k,i)) / max(abs(nevapr   (i,k)), abs(nevapr_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", nevapr   (i,k), " nevapr_permute: ", nevapr_permute(k,i)
!        endif
!        if ( abs(evapsnow (i,k) - evapsnow_permute(k,i)) / max(abs(evapsnow (i,k)), abs(evapsnow_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(evapsnow (i,k) - evapsnow_permute(k,i)) / max(abs(evapsnow (i,k)), abs(evapsnow_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", evapsnow (i,k), " evapsnow_permute: ", evapsnow_permute(k,i)
!        endif
!        if ( abs(am_evp_st(i,k) - am_evp_st_permute(k,i)) / max(abs(am_evp_st(i,k)), abs(am_evp_st_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(am_evp_st(i,k) - am_evp_st_permute(k,i)) / max(abs(am_evp_st(i,k)), abs(am_evp_st_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", am_evp_st(i,k), " am_evp_st_permute: ", am_evp_st_permute(k,i)
!        endif
!        if ( abs(prain    (i,k) - prain_permute(k,i)) / max(abs(prain    (i,k)), abs(prain_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(prain    (i,k) - prain_permute(k,i)) / max(abs(prain    (i,k)), abs(prain_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", prain    (i,k), " prain_permute: ", prain_permute(k,i)
!        endif
!        if ( abs(prodsnow (i,k) - prodsnow_permute(k,i)) / max(abs(prodsnow (i,k)), abs(prodsnow_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(prodsnow (i,k) - prodsnow_permute(k,i)) / max(abs(prodsnow (i,k)), abs(prodsnow_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", prodsnow (i,k), " prodsnow_permute: ", prodsnow_permute(k,i)
!        endif
!        if ( abs(cmeout   (i,k) - cmeout_permute(k,i)) / max(abs(cmeout   (i,k)), abs(cmeout_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(cmeout   (i,k) - cmeout_permute(k,i)) / max(abs(cmeout   (i,k)), abs(cmeout_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", cmeout   (i,k), " cmeout_permute: ", cmeout_permute(k,i)
!        endif
!        if ( abs(deffi    (i,k) - deffi_permute(k,i)) / max(abs(deffi    (i,k)), abs(deffi_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(deffi    (i,k) - deffi_permute(k,i)) / max(abs(deffi    (i,k)), abs(deffi_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", deffi    (i,k), " deffi_permute: ", deffi_permute(k,i)
!        endif
!        if ( abs(pgamrad  (i,k) - pgamrad_permute(k,i)) / max(abs(pgamrad  (i,k)), abs(pgamrad_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(pgamrad  (i,k) - pgamrad_permute(k,i)) / max(abs(pgamrad  (i,k)), abs(pgamrad_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", pgamrad  (i,k), " pgamrad_permute: ", pgamrad_permute(k,i)
!        endif
!        if ( abs(lamcrad  (i,k) - lamcrad_permute(k,i)) / max(abs(lamcrad  (i,k)), abs(lamcrad_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(lamcrad  (i,k) - lamcrad_permute(k,i)) / max(abs(lamcrad  (i,k)), abs(lamcrad_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", lamcrad  (i,k), " lamcrad_permute: ", lamcrad_permute(k,i)
!        endif
!        if ( abs(qsout    (i,k) - qsout_permute(k,i)) / max(abs(qsout    (i,k)), abs(qsout_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qsout    (i,k) - qsout_permute(k,i)) / max(abs(qsout    (i,k)), abs(qsout_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qsout    (i,k), " qsout_permute: ", qsout_permute(k,i)
!        endif
!        if ( abs(dsout    (i,k) - dsout_permute(k,i)) / max(abs(dsout    (i,k)), abs(dsout_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(dsout    (i,k) - dsout_permute(k,i)) / max(abs(dsout    (i,k)), abs(dsout_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", dsout    (i,k), " dsout_permute: ", dsout_permute(k,i)
!        endif
!        if ( abs(rflx     (i,k) - rflx_permute(k,i)) / max(abs(rflx     (i,k)), abs(rflx_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(rflx     (i,k) - rflx_permute(k,i)) / max(abs(rflx     (i,k)), abs(rflx_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", rflx     (i,k), " rflx_permute: ", rflx_permute(k,i)
!        endif
!        if ( abs(sflx     (i,k) - sflx_permute(k,i)) / max(abs(sflx     (i,k)), abs(sflx_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(sflx     (i,k) - sflx_permute(k,i)) / max(abs(sflx     (i,k)), abs(sflx_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", sflx     (i,k), " sflx_permute: ", sflx_permute(k,i)
!        endif
!        if ( abs(qrout    (i,k) - qrout_permute(k,i)) / max(abs(qrout    (i,k)), abs(qrout_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qrout    (i,k) - qrout_permute(k,i)) / max(abs(qrout    (i,k)), abs(qrout_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qrout    (i,k), " qrout_permute: ", qrout_permute(k,i)
!        endif
!        if ( abs(qcsevap  (i,k) - qcsevap_permute(k,i)) / max(abs(qcsevap  (i,k)), abs(qcsevap_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qcsevap  (i,k) - qcsevap_permute(k,i)) / max(abs(qcsevap  (i,k)), abs(qcsevap_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qcsevap  (i,k), " qcsevap_permute: ", qcsevap_permute(k,i)
!        endif
!        if ( abs(qisevap  (i,k) - qisevap_permute(k,i)) / max(abs(qisevap  (i,k)), abs(qisevap_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qisevap  (i,k) - qisevap_permute(k,i)) / max(abs(qisevap  (i,k)), abs(qisevap_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qisevap  (i,k), " qisevap_permute: ", qisevap_permute(k,i)
!        endif
!        if ( abs(qvres    (i,k) - qvres_permute(k,i)) / max(abs(qvres    (i,k)), abs(qvres_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qvres    (i,k) - qvres_permute(k,i)) / max(abs(qvres    (i,k)), abs(qvres_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qvres    (i,k), " qvres_permute: ", qvres_permute(k,i)
!        endif
!        if ( abs(cmeiout  (i,k) - cmeiout_permute(k,i)) / max(abs(cmeiout  (i,k)), abs(cmeiout_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(cmeiout  (i,k) - cmeiout_permute(k,i)) / max(abs(cmeiout  (i,k)), abs(cmeiout_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", cmeiout  (i,k), " cmeiout_permute: ", cmeiout_permute(k,i)
!        endif
!        if ( abs(vtrmc    (i,k) - vtrmc_permute(k,i)) / max(abs(vtrmc    (i,k)), abs(vtrmc_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(vtrmc    (i,k) - vtrmc_permute(k,i)) / max(abs(vtrmc    (i,k)), abs(vtrmc_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", vtrmc    (i,k), " vtrmc_permute: ", vtrmc_permute(k,i)
!        endif
!        if ( abs(vtrmi    (i,k) - vtrmi_permute(k,i)) / max(abs(vtrmi    (i,k)), abs(vtrmi_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(vtrmi    (i,k) - vtrmi_permute(k,i)) / max(abs(vtrmi    (i,k)), abs(vtrmi_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", vtrmi    (i,k), " vtrmi_permute: ", vtrmi_permute(k,i)
!        endif
!        if ( abs(qcsedten (i,k) - qcsedten_permute(k,i)) / max(abs(qcsedten (i,k)), abs(qcsedten_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qcsedten (i,k) - qcsedten_permute(k,i)) / max(abs(qcsedten (i,k)), abs(qcsedten_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qcsedten (i,k), " qcsedten_permute: ", qcsedten_permute(k,i)
!        endif
!        if ( abs(qisedten (i,k) - qisedten_permute(k,i)) / max(abs(qisedten (i,k)), abs(qisedten_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qisedten (i,k) - qisedten_permute(k,i)) / max(abs(qisedten (i,k)), abs(qisedten_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qisedten (i,k), " qisedten_permute: ", qisedten_permute(k,i)
!        endif
!        if ( abs(prao     (i,k) - prao_permute(k,i)) / max(abs(prao     (i,k)), abs(prao_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(prao     (i,k) - prao_permute(k,i)) / max(abs(prao     (i,k)), abs(prao_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", prao     (i,k), " prao_permute: ", prao_permute(k,i)
!        endif
!        if ( abs(prco     (i,k) - prco_permute(k,i)) / max(abs(prco     (i,k)), abs(prco_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(prco     (i,k) - prco_permute(k,i)) / max(abs(prco     (i,k)), abs(prco_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", prco     (i,k), " prco_permute: ", prco_permute(k,i)
!        endif
!        if ( abs(mnuccco  (i,k) - mnuccco_permute(k,i)) / max(abs(mnuccco  (i,k)), abs(mnuccco_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(mnuccco  (i,k) - mnuccco_permute(k,i)) / max(abs(mnuccco  (i,k)), abs(mnuccco_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", mnuccco  (i,k), " mnuccco_permute: ", mnuccco_permute(k,i)
!        endif
!        if ( abs(mnuccto  (i,k) - mnuccto_permute(k,i)) / max(abs(mnuccto  (i,k)), abs(mnuccto_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(mnuccto  (i,k) - mnuccto_permute(k,i)) / max(abs(mnuccto  (i,k)), abs(mnuccto_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", mnuccto  (i,k), " mnuccto_permute: ", mnuccto_permute(k,i)
!        endif
!        if ( abs(msacwio  (i,k) - msacwio_permute(k,i)) / max(abs(msacwio  (i,k)), abs(msacwio_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(msacwio  (i,k) - msacwio_permute(k,i)) / max(abs(msacwio  (i,k)), abs(msacwio_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", msacwio  (i,k), " msacwio_permute: ", msacwio_permute(k,i)
!        endif
!        if ( abs(psacwso  (i,k) - psacwso_permute(k,i)) / max(abs(psacwso  (i,k)), abs(psacwso_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(psacwso  (i,k) - psacwso_permute(k,i)) / max(abs(psacwso  (i,k)), abs(psacwso_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", psacwso  (i,k), " psacwso_permute: ", psacwso_permute(k,i)
!        endif
!        if ( abs(bergso   (i,k) - bergso_permute(k,i)) / max(abs(bergso   (i,k)), abs(bergso_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(bergso   (i,k) - bergso_permute(k,i)) / max(abs(bergso   (i,k)), abs(bergso_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", bergso   (i,k), " bergso_permute: ", bergso_permute(k,i)
!        endif
!        if ( abs(bergo    (i,k) - bergo_permute(k,i)) / max(abs(bergo    (i,k)), abs(bergo_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(bergo    (i,k) - bergo_permute(k,i)) / max(abs(bergo    (i,k)), abs(bergo_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", bergo    (i,k), " bergo_permute: ", bergo_permute(k,i)
!        endif
!        if ( abs(melto    (i,k) - melto_permute(k,i)) / max(abs(melto    (i,k)), abs(melto_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(melto    (i,k) - melto_permute(k,i)) / max(abs(melto    (i,k)), abs(melto_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", melto    (i,k), " melto_permute: ", melto_permute(k,i)
!        endif
!        if ( abs(homoo    (i,k) - homoo_permute(k,i)) / max(abs(homoo    (i,k)), abs(homoo_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(homoo    (i,k) - homoo_permute(k,i)) / max(abs(homoo    (i,k)), abs(homoo_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", homoo    (i,k), " homoo_permute: ", homoo_permute(k,i)
!        endif
!        if ( abs(qcreso   (i,k) - qcreso_permute(k,i)) / max(abs(qcreso   (i,k)), abs(qcreso_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qcreso   (i,k) - qcreso_permute(k,i)) / max(abs(qcreso   (i,k)), abs(qcreso_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qcreso   (i,k), " qcreso_permute: ", qcreso_permute(k,i)
!        endif
!        if ( abs(prcio    (i,k) - prcio_permute(k,i)) / max(abs(prcio    (i,k)), abs(prcio_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(prcio    (i,k) - prcio_permute(k,i)) / max(abs(prcio    (i,k)), abs(prcio_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", prcio    (i,k), " prcio_permute: ", prcio_permute(k,i)
!        endif
!        if ( abs(praio    (i,k) - praio_permute(k,i)) / max(abs(praio    (i,k)), abs(praio_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(praio    (i,k) - praio_permute(k,i)) / max(abs(praio    (i,k)), abs(praio_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", praio    (i,k), " praio_permute: ", praio_permute(k,i)
!        endif
!        if ( abs(qireso   (i,k) - qireso_permute(k,i)) / max(abs(qireso   (i,k)), abs(qireso_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qireso   (i,k) - qireso_permute(k,i)) / max(abs(qireso   (i,k)), abs(qireso_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qireso   (i,k), " qireso_permute: ", qireso_permute(k,i)
!        endif
!        if ( abs(mnuccro  (i,k) - mnuccro_permute(k,i)) / max(abs(mnuccro  (i,k)), abs(mnuccro_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(mnuccro  (i,k) - mnuccro_permute(k,i)) / max(abs(mnuccro  (i,k)), abs(mnuccro_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", mnuccro  (i,k), " mnuccro_permute: ", mnuccro_permute(k,i)
!        endif
!        if ( abs(pracso   (i,k) - pracso_permute(k,i)) / max(abs(pracso   (i,k)), abs(pracso_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(pracso   (i,k) - pracso_permute(k,i)) / max(abs(pracso   (i,k)), abs(pracso_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", pracso   (i,k), " pracso_permute: ", pracso_permute(k,i)
!        endif
!        if ( abs(meltsdt  (i,k) - meltsdt_permute(k,i)) / max(abs(meltsdt  (i,k)), abs(meltsdt_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(meltsdt  (i,k) - meltsdt_permute(k,i)) / max(abs(meltsdt  (i,k)), abs(meltsdt_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", meltsdt  (i,k), " meltsdt_permute: ", meltsdt_permute(k,i)
!        endif
!        if ( abs(frzrdt   (i,k) - frzrdt_permute(k,i)) / max(abs(frzrdt   (i,k)), abs(frzrdt_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(frzrdt   (i,k) - frzrdt_permute(k,i)) / max(abs(frzrdt   (i,k)), abs(frzrdt_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", frzrdt   (i,k), " frzrdt_permute: ", frzrdt_permute(k,i)
!        endif
!        if ( abs(mnuccdo  (i,k) - mnuccdo_permute(k,i)) / max(abs(mnuccdo  (i,k)), abs(mnuccdo_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(mnuccdo  (i,k) - mnuccdo_permute(k,i)) / max(abs(mnuccdo  (i,k)), abs(mnuccdo_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", mnuccdo  (i,k), " mnuccdo_permute: ", mnuccdo_permute(k,i)
!        endif
!        if ( abs(nrout    (i,k) - nrout_permute(k,i)) / max(abs(nrout    (i,k)), abs(nrout_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nrout    (i,k) - nrout_permute(k,i)) / max(abs(nrout    (i,k)), abs(nrout_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", nrout    (i,k), " nrout_permute: ", nrout_permute(k,i)
!        endif
!        if ( abs(nsout    (i,k) - nsout_permute(k,i)) / max(abs(nsout    (i,k)), abs(nsout_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nsout    (i,k) - nsout_permute(k,i)) / max(abs(nsout    (i,k)), abs(nsout_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", nsout    (i,k), " nsout_permute: ", nsout_permute(k,i)
!        endif
!        if ( abs(refl     (i,k) - refl_permute(k,i)) / max(abs(refl     (i,k)), abs(refl_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(refl     (i,k) - refl_permute(k,i)) / max(abs(refl     (i,k)), abs(refl_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", refl     (i,k), " refl_permute: ", refl_permute(k,i)
!        endif
!        if ( abs(arefl    (i,k) - arefl_permute(k,i)) / max(abs(arefl    (i,k)), abs(arefl_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(arefl    (i,k) - arefl_permute(k,i)) / max(abs(arefl    (i,k)), abs(arefl_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", arefl    (i,k), " arefl_permute: ", arefl_permute(k,i)
!        endif
!        if ( abs(areflz   (i,k) - areflz_permute(k,i)) / max(abs(areflz   (i,k)), abs(areflz_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(areflz   (i,k) - areflz_permute(k,i)) / max(abs(areflz   (i,k)), abs(areflz_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", areflz   (i,k), " areflz_permute: ", areflz_permute(k,i)
!        endif
!        if ( abs(frefl    (i,k) - frefl_permute(k,i)) / max(abs(frefl    (i,k)), abs(frefl_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(frefl    (i,k) - frefl_permute(k,i)) / max(abs(frefl    (i,k)), abs(frefl_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", frefl    (i,k), " frefl_permute: ", frefl_permute(k,i)
!        endif
!        if ( abs(csrfl    (i,k) - csrfl_permute(k,i)) / max(abs(csrfl    (i,k)), abs(csrfl_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(csrfl    (i,k) - csrfl_permute(k,i)) / max(abs(csrfl    (i,k)), abs(csrfl_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", csrfl    (i,k), " csrfl_permute: ", csrfl_permute(k,i)
!        endif
!        if ( abs(acsrfl   (i,k) - acsrfl_permute(k,i)) / max(abs(acsrfl   (i,k)), abs(acsrfl_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(acsrfl   (i,k) - acsrfl_permute(k,i)) / max(abs(acsrfl   (i,k)), abs(acsrfl_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", acsrfl   (i,k), " acsrfl_permute: ", acsrfl_permute(k,i)
!        endif
!        if ( abs(fcsrfl   (i,k) - fcsrfl_permute(k,i)) / max(abs(fcsrfl   (i,k)), abs(fcsrfl_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(fcsrfl   (i,k) - fcsrfl_permute(k,i)) / max(abs(fcsrfl   (i,k)), abs(fcsrfl_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", fcsrfl   (i,k), " fcsrfl_permute: ", fcsrfl_permute(k,i)
!        endif
!        if ( abs(rercld   (i,k) - rercld_permute(k,i)) / max(abs(rercld   (i,k)), abs(rercld_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(rercld   (i,k) - rercld_permute(k,i)) / max(abs(rercld   (i,k)), abs(rercld_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", rercld   (i,k), " rercld_permute: ", rercld_permute(k,i)
!        endif
!        if ( abs(ncai     (i,k) - ncai_permute(k,i)) / max(abs(ncai     (i,k)), abs(ncai_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(ncai     (i,k) - ncai_permute(k,i)) / max(abs(ncai     (i,k)), abs(ncai_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", ncai     (i,k), " ncai_permute: ", ncai_permute(k,i)
!        endif
!        if ( abs(ncal     (i,k) - ncal_permute(k,i)) / max(abs(ncal     (i,k)), abs(ncal_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(ncal     (i,k) - ncal_permute(k,i)) / max(abs(ncal     (i,k)), abs(ncal_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", ncal     (i,k), " ncal_permute: ", ncal_permute(k,i)
!        endif
!        if ( abs(qrout2   (i,k) - qrout2_permute(k,i)) / max(abs(qrout2   (i,k)), abs(qrout2_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qrout2   (i,k) - qrout2_permute(k,i)) / max(abs(qrout2   (i,k)), abs(qrout2_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qrout2   (i,k), " qrout2_permute: ", qrout2_permute(k,i)
!        endif
!        if ( abs(qsout2   (i,k) - qsout2_permute(k,i)) / max(abs(qsout2   (i,k)), abs(qsout2_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(qsout2   (i,k) - qsout2_permute(k,i)) / max(abs(qsout2   (i,k)), abs(qsout2_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", qsout2   (i,k), " qsout2_permute: ", qsout2_permute(k,i)
!        endif
!        if ( abs(nrout2   (i,k) - nrout2_permute(k,i)) / max(abs(nrout2   (i,k)), abs(nrout2_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nrout2   (i,k) - nrout2_permute(k,i)) / max(abs(nrout2   (i,k)), abs(nrout2_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", nrout2   (i,k), " nrout2_permute: ", nrout2_permute(k,i)
!        endif
!        if ( abs(nsout2   (i,k) - nsout2_permute(k,i)) / max(abs(nsout2   (i,k)), abs(nsout2_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nsout2   (i,k) - nsout2_permute(k,i)) / max(abs(nsout2   (i,k)), abs(nsout2_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", nsout2   (i,k), " nsout2_permute: ", nsout2_permute(k,i)
!        endif
!        if ( abs(drout2   (i,k) - drout2_permute(k,i)) / max(abs(drout2   (i,k)), abs(drout2_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(drout2   (i,k) - drout2_permute(k,i)) / max(abs(drout2   (i,k)), abs(drout2_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", drout2   (i,k), " drout2_permute: ", drout2_permute(k,i)
!        endif
!        if ( abs(dsout2   (i,k) - dsout2_permute(k,i)) / max(abs(dsout2   (i,k)), abs(dsout2_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(dsout2   (i,k) - dsout2_permute(k,i)) / max(abs(dsout2   (i,k)), abs(dsout2_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", dsout2   (i,k), " dsout2_permute: ", dsout2_permute(k,i)
!        endif
!        if ( abs(freqs    (i,k) - freqs_permute(k,i)) / max(abs(freqs    (i,k)), abs(freqs_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(freqs    (i,k) - freqs_permute(k,i)) / max(abs(freqs    (i,k)), abs(freqs_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", freqs    (i,k), " freqs_permute: ", freqs_permute(k,i)
!        endif
!        if ( abs(freqr    (i,k) - freqr_permute(k,i)) / max(abs(freqr    (i,k)), abs(freqr_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(freqr    (i,k) - freqr_permute(k,i)) / max(abs(freqr    (i,k)), abs(freqr_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", freqr    (i,k), " freqr_permute: ", freqr_permute(k,i)
!        endif
!        if ( abs(nfice    (i,k) - nfice_permute(k,i)) / max(abs(nfice    (i,k)), abs(nfice_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(nfice    (i,k) - nfice_permute(k,i)) / max(abs(nfice    (i,k)), abs(nfice_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", nfice    (i,k), " nfice_permute: ", nfice_permute(k,i)
!        endif
!        if ( abs(prer_evap(i,k) - prer_evap_permute(k,i)) / max(abs(prer_evap(i,k)), abs(prer_evap_permute(k,i))) > 1e-10_r8 ) then
!            print*, "err:", abs(prer_evap(i,k) - prer_evap_permute(k,i)) / max(abs(prer_evap(i,k)), abs(prer_evap_permute(k,i))), "i: ", i, "k: ", k, "top_lev: ", top_lev, "pver: ", pver, ": ", prer_evap(i,k), " prer_evap_permute: ", prer_evap_permute(k,i)
!        endif
!        !qrout2(i,k) = qrout2_permute(k,i)
!        !qsout2(i,k) = qsout2_permute(k,i)
!
!    end do
!enddo
!endif
!call set_consist_math_off()


!! initialize  output fields for number conc qand ice nucleation
!ncai_permute(1:pver,1:ncol)=0._r8 
!ncal_permute(1:pver,1:ncol)=0._r8  
!
!!Initialize rain size
!rercld_permute(1:pver,1:ncol)=0._r8
!arcld_permute(1:pver,1:ncol)=0._r8
!
!!initialize radiation output variables
!pgamrad_permute(1:pver,1:ncol)=0._r8 ! liquid gamma parameter for optics (radiation)
!lamcrad_permute(1:pver,1:ncol)=0._r8 ! slope of droplet distribution for optics (radiation)
!deffi_permute  (1:pver,1:ncol)=0._r8 ! slope of droplet distribution for optics (radiation)
!!initialize radiation output variables
!!initialize water vapor tendency term output
!qcsevap_permute(1:pver,1:ncol)=0._r8 
!qisevap_permute(1:pver,1:ncol)=0._r8 
!qvres_permute  (1:pver,1:ncol)=0._r8 
!cmeiout_permute (1:pver,1:ncol)=0._r8
!vtrmc_permute (1:pver,1:ncol)=0._r8
!vtrmi_permute (1:pver,1:ncol)=0._r8
!qcsedten_permute (1:pver,1:ncol)=0._r8
!qisedten_permute (1:pver,1:ncol)=0._r8    
!
!prao_permute(1:pver,1:ncol)=0._r8 
!prco_permute(1:pver,1:ncol)=0._r8 
!mnuccco_permute(1:pver,1:ncol)=0._r8 
!mnuccto_permute(1:pver,1:ncol)=0._r8 
!msacwio_permute(1:pver,1:ncol)=0._r8 
!psacwso_permute(1:pver,1:ncol)=0._r8 
!bergso_permute(1:pver,1:ncol)=0._r8 
!bergo_permute(1:pver,1:ncol)=0._r8 
!melto_permute(1:pver,1:ncol)=0._r8 
!homoo_permute(1:pver,1:ncol)=0._r8 
!qcreso_permute(1:pver,1:ncol)=0._r8 
!prcio_permute(1:pver,1:ncol)=0._r8 
!praio_permute(1:pver,1:ncol)=0._r8 
!qireso_permute(1:pver,1:ncol)=0._r8 
!mnuccro_permute(1:pver,1:ncol)=0._r8 
!pracso_permute (1:pver,1:ncol)=0._r8 
!meltsdt_permute(1:pver,1:ncol)=0._r8
!frzrdt_permute (1:pver,1:ncol)=0._r8
!mnuccdo_permute(1:pver,1:ncol)=0._r8
!
!rflx_permute(:,:)=0._r8
!sflx_permute(:,:)=0._r8
!effc_permute(:,:)=0._r8
!effc_fn_permute(:,:)=0._r8
!effi_permute(:,:)=0._r8
!
!! assign variable deltat for sub-stepping...
!deltat=deltatin
!
!! parameters for scheme
!
!omsm=0.99999_r8
!dto2=0.5_r8*deltat
!mincld=0.0001_r8
!
!! initialize multi-level fields
!q_permute(1:pver,1:ncol)=qn_permute(1:pver,1:ncol)
!t_permute(1:pver,1:ncol)=tn_permute(1:pver,1:ncol)
!
!! initialize time-varying parameters
!
!do k=1,pver
!   do i=1,ncol
!      rho_permute(k,i)=p_permute(k,i)/(r*t_permute(k,i))
!      dv_permute(k,i) = 8.794E-5_r8*t_permute(k,i)**1.81_r8/p_permute(k,i)
!      mu_permute(k,i) = 1.496E-6_r8*t_permute(k,i)**1.5_r8/(t_permute(k,i)+120._r8) 
!      sc_permute(k,i) = mu_permute(k,i)/(rho_permute(k,i)*dv_permute(k,i))
!      kap_permute(k,i) = 1.414e3_r8*1.496e-6_r8*t_permute(k,i)**1.5_r8/(t_permute(k,i)+120._r8) 
!
!      ! air density adjustment for fallspeed parameters
!      ! includes air density correction factor to the
!      ! power of 0.54 following Heymsfield and Bansemer 2007
!
!      rhof_permute(k,i)=(rhosu/rho_permute(k,i))**0.54_r8
!
!      arn_permute(k,i)=ar*rhof_permute(k,i)
!      asn_permute(k,i)=as*rhof_permute(k,i)
!      acn_permute(k,i)=ac*rhof_permute(k,i)
!      ain_permute(k,i)=ai*rhof_permute(k,i)
!
!      ! get dz from dp and hydrostatic approx
!      ! keep dz positive (define as layer k-1 - layer k)
!
!      dz_permute(k,i)= pdel_permute(k,i)/(rho_permute(k,i)*g)
!
!   end do
!end do
!
!! initialization
!qc_permute(1:top_lev-1,1:ncol) = 0._r8
!qi_permute(1:top_lev-1,1:ncol) = 0._r8
!nc_permute(1:top_lev-1,1:ncol) = 0._r8
!ni_permute(1:top_lev-1,1:ncol) = 0._r8
!t1_permute(1:pver,1:ncol) = t_permute(1:pver,1:ncol)
!q1_permute(1:pver,1:ncol) = q_permute(1:pver,1:ncol)
!qc1_permute(1:pver,1:ncol) = qc_permute(1:pver,1:ncol)
!qi1_permute(1:pver,1:ncol) = qi_permute(1:pver,1:ncol)
!nc1_permute(1:pver,1:ncol) = nc_permute(1:pver,1:ncol)
!ni1_permute(1:pver,1:ncol) = ni_permute(1:pver,1:ncol)
!
!! initialize tendencies to zero
!tlat1_permute(1:pver,1:ncol)=0._r8
!qvlat1_permute(1:pver,1:ncol)=0._r8
!qctend1_permute(1:pver,1:ncol)=0._r8
!qitend1_permute(1:pver,1:ncol)=0._r8
!nctend1_permute(1:pver,1:ncol)=0._r8
!nitend1_permute(1:pver,1:ncol)=0._r8
!
!! initialize precip output
!qrout_permute(1:pver,1:ncol)=0._r8
!qsout_permute(1:pver,1:ncol)=0._r8
!nrout_permute(1:pver,1:ncol)=0._r8
!nsout_permute(1:pver,1:ncol)=0._r8
!dsout_permute(1:pver,1:ncol)=0._r8
!
!drout_permute(1:pver,1:ncol)=0._r8
!
!reff_rain_permute(1:pver,1:ncol)=0._r8
!reff_snow_permute(1:pver,1:ncol)=0._r8
!
!! initialize variables for trop_mozart
!nevapr_permute(1:pver,1:ncol) = 0._r8
!nevapr2_permute(1:pver,1:ncol) = 0._r8
!evapsnow_permute(1:pver,1:ncol) = 0._r8
!prain_permute(1:pver,1:ncol) = 0._r8
!prodsnow_permute(1:pver,1:ncol) = 0._r8
!cmeout_permute(1:pver,1:ncol) = 0._r8
!
!am_evp_st_permute(1:pver,1:ncol) = 0._r8
!
!! for refl calc
!rainrt1_permute(1:pver,1:ncol) = 0._r8
!
!! initialize precip fraction and output tendencies
!cldmax_permute(1:pver,1:ncol)=mincld
!
!!initialize aerosol number
!!        naer2(1:ncol,1:pver,:)=0._r8
!dum2l_permute(1:pver,1:ncol)=0._r8
!dum2i_permute(1:pver,1:ncol)=0._r8
!
!! initialize avg precip rate
!prect1(1:ncol)=0._r8
!preci1(1:ncol)=0._r8
!
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!Get humidity and saturation vapor pressures
!
!do k=top_lev,pver
!
!   do i=1,ncol
!
!      ! find wet bulk temperature and saturation value for provisional t and q without
!      ! condensation
!      
!      es(i) = svp_water(t_permute(k,i))
!      qs(i) = svp_to_qsat(es(i), p_permute(k,i))
!
!      ! Prevents negative values.
!      if (qs(i) < 0.0_r8) then
!         qs(i) = 1.0_r8
!         es(i) = p_permute(k,i)
!      end if
!
!      esl_permute(k,i)=svp_water(t_permute(k,i))
!      esi_permute(k,i)=svp_ice(t_permute(k,i))
!
!      ! hm fix, make sure when above freezing that esi=esl, not active yet
!      if (t_permute(k,i).gt.tmelt)esi_permute(k,i)=esl_permute(k,i)
!
!      relhum_permute(k,i)=q_permute(k,i)/qs(i)
!
!      ! get cloud fraction, check for minimum
!
!      cldm_permute(k,i)=max(cldn_permute(k,i),mincld)
!      cldmw_permute(k,i)=max(cldn_permute(k,i),mincld)
!
!      icldm_permute(k,i)=max(icecldf_permute(k,i),mincld)
!      lcldm_permute(k,i)=max(liqcldf(i,k),mincld)
!
!      ! subcolumns, set cloud fraction variables to one
!      ! if cloud water or ice is present, if not present
!      ! set to mincld (mincld used instead of zero, to prevent
!      ! possible division by zero errors
!
!      if (microp_uniform) then
!
!         cldm_permute(k,i)=mincld
!         cldmw_permute(k,i)=mincld
!         icldm_permute(k,i)=mincld
!         lcldm_permute(k,i)=mincld
!
!         if (qc_permute(k,i).ge.qsmall) then
!            lcldm_permute(k,i)=1._r8           
!            cldm_permute(k,i)=1._r8
!            cldmw_permute(k,i)=1._r8
!         end if
!
!         if (qi_permute(k,i).ge.qsmall) then             
!            cldm_permute(k,i)=1._r8
!            icldm_permute(k,i)=1._r8
!         end if
!
!      end if               ! sub-columns
!
!      ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
!
!      nfice_permute(k,i)=0._r8
!      dumfice=qc_permute(k,i)+qi_permute(k,i)
!      if (dumfice.gt.qsmall .and. qi_permute(k,i).gt.qsmall) then
!         nfice_permute(k,i)=qi_permute(k,i)/dumfice
!      endif
!
!      if (do_cldice .and. (t_permute(k,i).lt.tmelt - 5._r8)) then
!
!         ! if aerosols interact with ice set number of activated ice nuclei
!         dum2=naai_permute(k,i)
!
!         dumnnuc=(dum2 - ni_permute(k,i)/icldm_permute(k,i))/deltat*icldm_permute(k,i)
!         dumnnuc=max(dumnnuc,0._r8)
!         ! get provisional ni and qi after nucleation in order to calculate
!         ! Bergeron process below
!         ninew=ni_permute(k,i)+dumnnuc*deltat
!         qinew=qi_permute(k,i)+dumnnuc*deltat*mi0
!
!         !T>268
!      else
!         ninew=ni_permute(k,i)
!         qinew=qi_permute(k,i)
!      end if
!
!      ! Initialize CME components
!
!      cme_permute(k,i) = 0._r8
!      cmei_permute(k,i)=0._r8
!
!
!      !-------------------------------------------------------------------
!      !Bergeron process
!
!      ! make sure to initialize bergeron process to zero
!      berg_permute(k,i)=0._r8
!      prd = 0._r8
!
!      !condensation loop.
!
!      ! get in-cloud qi and ni after nucleation
!      if (icldm_permute(k,i) .gt. 0._r8) then 
!         qiic_permute(k,i)=qinew/icldm_permute(k,i)
!         niic_permute(k,i)=ninew/icldm_permute(k,i)
!      else
!         qiic_permute(k,i)=0._r8
!         niic_permute(k,i)=0._r8
!      endif
!
!      !if T < 0 C then bergeron.
!      if (do_cldice .and. (t_permute(k,i).lt.273.15_r8)) then
!
!         !if ice exists
!         if (qi_permute(k,i).gt.qsmall) then
!
!            bergtsf = 0._r8 ! bergeron time scale (fraction of timestep)
!
!            qvi = svp_to_qsat(esi_permute(k,i), p_permute(k,i))
!            qvl = svp_to_qsat(esl_permute(k,i), p_permute(k,i))
!
!            dqsidt =  xxls*qvi/(rv*t_permute(k,i)**2)
!            abi = 1._r8+dqsidt*xxls/cpp
!
!            ! get ice size distribution parameters
!
!            if (qiic_permute(k,i).ge.qsmall) then
!               lami(k) = (cons1*ci* &
!                    niic_permute(k,i)/qiic_permute(k,i))**(1._r8/di)
!               n0i(k) = niic_permute(k,i)*lami(k)
!
!               ! check for slope
!               ! adjust vars
!               if (lami(k).lt.lammini) then
!
!                  lami(k) = lammini
!                  n0i(k) = lami(k)**(di+1._r8)*qiic_permute(k,i)/(ci*cons1)
!               else if (lami(k).gt.lammaxi) then
!                  lami(k) = lammaxi
!                  n0i(k) = lami(k)**(di+1._r8)*qiic_permute(k,i)/(ci*cons1)
!               end if
!
!               epsi = 2._r8*pi*n0i(k)*rho_permute(k,i)*dv_permute(k,i)/(lami(k)*lami(k))
!
!               !if liquid exists  
!               if (qc_permute(k,i).gt. qsmall) then 
!
!                  !begin bergeron process
!                  !     do bergeron (vapor deposition with RHw=1)
!                  !     code to find berg (a rate) goes here
!
!                  ! calculate Bergeron process
!
!                  prd = epsi*(qvl-qvi)/abi
!
!               else
!                  prd = 0._r8
!               end if
!
!               ! multiply by cloud fraction
!
!               prd = prd*min(icldm_permute(k,i),lcldm_permute(k,i))
!
!               !     transfer of existing cloud liquid to ice
!
!               berg_permute(k,i)=max(0._r8,prd)
!
!            end if  !end liquid exists bergeron
!
!            if (berg_permute(k,i).gt.0._r8) then
!               bergtsf=max(0._r8,(qc_permute(k,i)/berg_permute(k,i))/deltat) 
!
!               if(bergtsf.lt.1._r8) berg_permute(k,i) = max(0._r8,qc_permute(k,i)/deltat)
!
!            endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            if (bergtsf.lt.1._r8.or.icldm_permute(k,i).gt.lcldm_permute(k,i)) then
!
!               if (qiic_permute(k,i).ge.qsmall) then
!
!                  ! first case is for case when liquid water is present, but is completely depleted 
!                  ! in time step, i.e., bergrsf > 0 but < 1
!
!                  if (qc_permute(k,i).ge.qsmall) then
!                     rhin  = (1.0_r8 + relhum_permute(k,i)) / 2._r8
!                     if ((rhin*esl_permute(k,i)/esi_permute(k,i)) > 1._r8) then
!                        prd = epsi*(rhin*qvl-qvi)/abi
!
!                        ! multiply by cloud fraction assuming liquid/ice maximum overlap
!                        prd = prd*min(icldm_permute(k,i),lcldm_permute(k,i))
!
!                        ! add to cmei
!                        cmei_permute(k,i) = cmei_permute(k,i) + (prd * (1._r8- bergtsf))
!
!                     end if ! rhin 
!                  end if ! qc > qsmall
!
!                  ! second case is for pure ice cloud, either no liquid, or icldm > lcldm
!
!                  if (qc_permute(k,i).lt.qsmall.or.icldm_permute(k,i).gt.lcldm_permute(k,i)) then
!
!                     ! note: for case of no liquid, need to set liquid cloud fraction to zero
!                     ! store liquid cloud fraction in 'dum'
!
!                     if (qc_permute(k,i).lt.qsmall) then 
!                        dum=0._r8 
!                     else
!                        dum=lcldm_permute(k,i)
!                     end if
!
!                     ! set RH to grid-mean value for pure ice cloud
!                     rhin = relhum_permute(k,i)
!
!                     if ((rhin*esl_permute(k,i)/esi_permute(k,i)) > 1._r8) then
!
!                        prd = epsi*(rhin*qvl-qvi)/abi
!
!                        ! multiply by relevant cloud fraction for pure ice cloud
!                        ! assuming maximum overlap of liquid/ice
!                        prd = prd*max((icldm_permute(k,i)-dum),0._r8)
!                        cmei_permute(k,i) = cmei_permute(k,i) + prd
!
!                     end if ! rhin
!                  end if ! qc or icldm > lcldm
!               end if ! qiic
!            end if ! bergtsf or icldm > lcldm
!
!            !     if deposition, it should not reduce grid mean rhi below 1.0
!            if(cmei_permute(k,i) > 0.0_r8 .and. (relhum_permute(k,i)*esl_permute(k,i)/esi_permute(k,i)) > 1._r8 ) &
!                 cmei_permute(k,i)=min(cmei_permute(k,i),(q_permute(k,i)-qs(i)*esi_permute(k,i)/esl_permute(k,i))/abi/deltat)
!
!         end if            !end ice exists loop
!         !this ends temperature < 0. loop
!
!         !-------------------------------------------------------------------
!      end if  ! 
!      !..............................................................
!
!      ! evaporation should not exceed available water
!
!      if ((-berg_permute(k,i)).lt.-qc_permute(k,i)/deltat) berg_permute(k,i) = max(qc_permute(k,i)/deltat,0._r8)
!
!      !sublimation process...
!      if (do_cldice .and. ((relhum_permute(k,i)*esl_permute(k,i)/esi_permute(k,i)).lt.1._r8 .and. qiic_permute(k,i).ge.qsmall )) then
!
!         qvi = svp_to_qsat(esi_permute(k,i), p_permute(k,i))
!         qvl = svp_to_qsat(esl_permute(k,i), p_permute(k,i))
!         dqsidt =  xxls*qvi/(rv*t_permute(k,i)**2)
!         abi = 1._r8+dqsidt*xxls/cpp
!
!         ! get ice size distribution parameters
!
!         lami(k) = (cons1*ci* &
!              niic_permute(k,i)/qiic_permute(k,i))**(1._r8/di)
!         n0i(k) = niic_permute(k,i)*lami(k)
!
!         ! check for slope
!         ! adjust vars
!         if (lami(k).lt.lammini) then
!
!            lami(k) = lammini
!            n0i(k) = lami(k)**(di+1._r8)*qiic_permute(k,i)/(ci*cons1)
!         else if (lami(k).gt.lammaxi) then
!            lami(k) = lammaxi
!            n0i(k) = lami(k)**(di+1._r8)*qiic_permute(k,i)/(ci*cons1)
!         end if
!
!         epsi = 2._r8*pi*n0i(k)*rho_permute(k,i)*dv_permute(k,i)/(lami(k)*lami(k))
!
!         ! modify for ice fraction below
!         prd = epsi*(relhum_permute(k,i)*qvl-qvi)/abi * icldm_permute(k,i)
!         cmei_permute(k,i)=min(prd,0._r8)
!
!      endif
!
!      ! sublimation should not exceed available ice
!      if (cmei_permute(k,i).lt.-qi_permute(k,i)/deltat) cmei_permute(k,i)=-qi_permute(k,i)/deltat
!
!      ! sublimation should not increase grid mean rhi above 1.0 
!      if(cmei_permute(k,i) < 0.0_r8 .and. (relhum_permute(k,i)*esl_permute(k,i)/esi_permute(k,i)) < 1._r8 ) &
!           cmei_permute(k,i)=min(0._r8,max(cmei_permute(k,i),(q_permute(k,i)-qs(i)*esi_permute(k,i)/esl_permute(k,i))/abi/deltat))
!
!      ! limit cmei due for roundoff error
!
!      cmei_permute(k,i)=cmei_permute(k,i)*omsm
!
!      ! conditional for ice nucleation 
!      if (do_cldice .and. (t_permute(k,i).lt.(tmelt - 5._r8))) then 
!
!         ! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
!         ! ice nucleation rate (dum2) has already been calculated and read in (naai)
!
!         dum2i_permute(k,i)=naai_permute(k,i)
!      else
!         dum2i_permute(k,i)=0._r8
!      end if
!
!   end do ! i loop
!end do ! k loop
!
!
!!! initialize sub-step precip flux variables
!do i=1,ncol
!   !! flux is zero at top interface, so these should stay as 0.
!   rflx1_permute(1,i)=0._r8
!   sflx1_permute(1,i)=0._r8
!   do k=top_lev,pver
!
!      ! initialize normal and sub-step precip flux variables
!      rflx1_permute(k+1,i)=0._r8
!      sflx1_permute(k+1,i)=0._r8
!   end do ! i loop
!end do ! k loop
!!! initialize final precip flux variables.
!do i=1,ncol
!   !! flux is zero at top interface, so these should stay as 0.
!   rflx_permute(1,i)=0._r8
!   sflx_permute(1,i)=0._r8
!   do k=top_lev,pver
!      ! initialize normal and sub-step precip flux variables
!      rflx_permute(k+1,i)=0._r8
!      sflx_permute(k+1,i)=0._r8
!   end do ! i loop
!end do ! k loop
!
!do i=1,ncol
!   ltrue(i)=0
!   do k=top_lev,pver
!      ! skip microphysical calculations if no cloud water
!
!      if (qc_permute(k,i).ge.qsmall.or.qi_permute(k,i).ge.qsmall.or.cmei_permute(k,i).ge.qsmall) ltrue(i)=1
!   end do
!end do
!
!! assign number of sub-steps to iter
!! use 2 sub-steps, following tests described in MG2008
!iter = 2
!
!! get sub-step time step
!deltat=deltat/real(iter)
!
!! since activation/nucleation processes are fast, need to take into account
!! factor mtime = mixing timescale in cloud / model time step
!! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
!! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk
!
!!        note: mtime for bulk aerosols was set to: mtime=deltat/1200._r8
!
!mtime=1._r8
!rate1ord_cw2pr_st_permute(:,:)=0._r8 ! rce 2010/05/01
!
!!!!! skip calculations if no cloud water
!do i=1,ncol
!   if (ltrue(i).eq.0) then
!      tlat_permute(1:pver,i)=0._r8
!      qvlat_permute(1:pver,i)=0._r8
!      qctend_permute(1:pver,i)=0._r8
!      qitend_permute(1:pver,i)=0._r8
!      qnitend_permute(1:pver,i)=0._r8
!      qrtend_permute(1:pver,i)=0._r8
!      nctend_permute(1:pver,i)=0._r8
!      nitend_permute(1:pver,i)=0._r8
!      nrtend_permute(1:pver,i)=0._r8
!      nstend_permute(1:pver,i)=0._r8
!      prect(i)=0._r8
!      preci(i)=0._r8
!      qniic_permute(1:pver,i)=0._r8
!      qric_permute(1:pver,i)=0._r8
!      nsic_permute(1:pver,i)=0._r8
!      nric_permute(1:pver,i)=0._r8
!      rainrt_permute(1:pver,i)=0._r8
!      goto 300
!   end if
!
!   qcsinksum_rate1ord(1:pver)=0._r8 
!   qcsum_rate1ord(1:pver)=0._r8 
!
!
!!!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !.....................................................................................................
!   do it=1,iter
!
!      ! initialize sub-step microphysical tendencies
!
!      tlat_permute(1:pver,i)=0._r8
!      qvlat_permute(1:pver,i)=0._r8
!      qctend_permute(1:pver,i)=0._r8
!      qitend_permute(1:pver,i)=0._r8
!      qnitend_permute(1:pver,i)=0._r8
!      qrtend_permute(1:pver,i)=0._r8
!      nctend_permute(1:pver,i)=0._r8
!      nitend_permute(1:pver,i)=0._r8
!      nrtend_permute(1:pver,i)=0._r8
!      nstend_permute(1:pver,i)=0._r8
!
!      ! initialize diagnostic precipitation to zero
!
!      qniic_permute(1:pver,i)=0._r8
!      qric_permute(1:pver,i)=0._r8
!      nsic_permute(1:pver,i)=0._r8
!      nric_permute(1:pver,i)=0._r8
!
!      rainrt_permute(1:pver,i)=0._r8
!
!
!      ! begin new i,k loop, calculate new cldmax after adjustment to cldm above
!
!      ! initialize vertically-integrated rain and snow tendencies
!
!      qrtot = 0._r8
!      nrtot = 0._r8
!      qstot = 0._r8
!      nstot = 0._r8
!
!      ! initialize precip at surface
!
!      prect(i)=0._r8
!      preci(i)=0._r8
!
!      do k=top_lev,pver
!      
!         qcvar=relvar_permute(k,i)
!         cons2=gamma(qcvar+2.47_r8)
!         cons3=gamma(qcvar)
!         cons9=gamma(qcvar+2._r8)
!         cons10=gamma(qcvar+1._r8)
!         cons12=gamma(qcvar+1.15_r8) 
!         cons15=gamma(qcvar+bc/3._r8)
!         cons18=qcvar**2.47_r8
!         cons19=qcvar**2
!         cons20=qcvar**1.15_r8
!
!         ! set cwml and cwmi to current qc and qi
!
!         cwml_permute(k,i)=qc_permute(k,i)
!         cwmi_permute(k,i)=qi_permute(k,i)
!
!         ! initialize precip fallspeeds to zero
!
!         ums(k)=0._r8 
!         uns(k)=0._r8 
!         umr(k)=0._r8 
!         unr(k)=0._r8
!
!         ! calculate precip fraction based on maximum overlap assumption
!
!         ! for sub-columns cldm has already been set to 1 if cloud
!         ! water or ice is present, so cldmax will be correctly set below
!         ! and nothing extra needs to be done here
!
!         if (k.eq.top_lev) then
!            cldmax_permute(k,i)=cldm_permute(k,i)
!         else
!            ! if rain or snow mix ratio is smaller than
!            ! threshold, then set cldmax to cloud fraction at current level
!
!            if (do_clubb_sgs) then
!               if (qc_permute(k,i).ge.qsmall.or.qi_permute(k,i).ge.qsmall) then
!                  cldmax_permute(k,i)=cldm_permute(k,i)
!               else
!                  cldmax_permute(k,i)=cldmax_permute(k-1,i)
!               end if
!            else
!
!               if (qric_permute(k-1,i).ge.qsmall.or.qniic_permute(k-1,i).ge.qsmall) then
!                  cldmax_permute(k,i)=max(cldmax_permute(k-1,i),cldm_permute(k,i))
!               else
!                  cldmax_permute(k,i)=cldm_permute(k,i)
!               end if
!            endif
!         end if
!
!         ! decrease in number concentration due to sublimation/evap
!         ! divide by cloud fraction to get in-cloud decrease
!         ! don't reduce Nc due to bergeron process
!
!         if (cmei_permute(k,i) < 0._r8 .and. qi_permute(k,i) > qsmall .and. cldm_permute(k,i) > mincld) then
!            nsubi(k)=cmei_permute(k,i)/qi_permute(k,i)*ni_permute(k,i)/cldm_permute(k,i)
!         else
!            nsubi(k)=0._r8
!         end if
!         nsubc(k)=0._r8
!
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
!
!         if (do_cldice .and. dum2i_permute(k,i).gt.0._r8.and.t_permute(k,i).lt.(tmelt - 5._r8).and. &
!              relhum_permute(k,i)*esl_permute(k,i)/esi_permute(k,i).gt. rhmini+0.05_r8) then
!
!            !if NCAI > 0. then set numice = ncai (as before)
!            !note: this is gridbox averaged
!
!            nnuccd(k)=(dum2i_permute(k,i)-ni_permute(k,i)/icldm_permute(k,i))/deltat*icldm_permute(k,i)
!            nnuccd(k)=max(nnuccd(k),0._r8)
!            nimax = dum2i_permute(k,i)*icldm_permute(k,i)
!
!            !Calc mass of new particles using new crystal mass...
!            !also this will be multiplied by mtime as nnuccd is...
!
!            mnuccd(k) = nnuccd(k) * mi0
!
!            !  add mnuccd to cmei....
!            cmei_permute(k,i)= cmei_permute(k,i) + mnuccd(k) * mtime
!
!            !  limit cmei
!
!            qvi = svp_to_qsat(esi_permute(k,i), p_permute(k,i))
!            dqsidt =  xxls*qvi/(rv*t_permute(k,i)**2)
!            abi = 1._r8+dqsidt*xxls/cpp
!            cmei_permute(k,i)=min(cmei_permute(k,i),(q_permute(k,i)-qvi)/abi/deltat)
!
!            ! limit for roundoff error
!            cmei_permute(k,i)=cmei_permute(k,i)*omsm
!
!         else
!            nnuccd(k)=0._r8
!            nimax = 0._r8
!            mnuccd(k) = 0._r8
!         end if
!
!         !c............................................................................
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
!         ! for microphysical process calculations
!         ! units are kg/kg for mixing ratio, 1/kg for number conc
!
!         ! limit in-cloud values to 0.005 kg/kg
!
!         qcic_permute(k,i)=min(cwml_permute(k,i)/lcldm_permute(k,i),5.e-3_r8)
!         qiic_permute(k,i)=min(cwmi_permute(k,i)/icldm_permute(k,i),5.e-3_r8)
!         ncic_permute(k,i)=max(nc_permute(k,i)/lcldm_permute(k,i),0._r8)
!         niic_permute(k,i)=max(ni_permute(k,i)/icldm_permute(k,i),0._r8)
!
!         if (qc_permute(k,i) - berg_permute(k,i)*deltat.lt.qsmall) then
!            qcic_permute(k,i)=0._r8
!            ncic_permute(k,i)=0._r8
!            if (qc_permute(k,i)-berg_permute(k,i)*deltat.lt.0._r8) then
!               berg_permute(k,i)=qc_permute(k,i)/deltat*omsm
!            end if
!         end if
!
!         if (do_cldice .and. qi_permute(k,i)+(cmei_permute(k,i)+berg_permute(k,i))*deltat.lt.qsmall) then
!            qiic_permute(k,i)=0._r8
!            niic_permute(k,i)=0._r8
!            if (qi_permute(k,i)+(cmei_permute(k,i)+berg_permute(k,i))*deltat.lt.0._r8) then
!               cmei_permute(k,i)=(-qi_permute(k,i)/deltat-berg_permute(k,i))*omsm
!            end if
!         end if
!
!         ! add to cme output
!
!         cmeout_permute(k,i) = cmeout_permute(k,i)+cmei_permute(k,i)
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! droplet activation
!         ! calculate potential for droplet activation if cloud water is present
!         ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
!         ! number tendency (npccnin) is read in from companion routine
!
!         ! assume aerosols already activated are equal to number of existing droplets for simplicity
!         ! multiply by cloud fraction to obtain grid-average tendency
!
!         if (qcic_permute(k,i).ge.qsmall) then   
!            npccn(k) = max(0._r8,npccnin_permute(k,i))  
!            dum2l_permute(k,i)=(nc_permute(k,i)+npccn(k)*deltat)/lcldm_permute(k,i)
!            dum2l_permute(k,i)=max(dum2l_permute(k,i),cdnl/rho_permute(k,i)) ! sghan minimum in #/cm3  
!            ncmax = dum2l_permute(k,i)*lcldm_permute(k,i)
!         else
!            npccn(k)=0._r8
!            dum2l_permute(k,i)=0._r8
!            ncmax = 0._r8
!         end if
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! get size distribution parameters based on in-cloud cloud water/ice 
!         ! these calculations also ensure consistency between number and mixing ratio
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         !......................................................................
!         ! cloud ice
!
!         if (qiic_permute(k,i).ge.qsmall) then
!
!            ! add upper limit to in-cloud number concentration to prevent numerical error
!            niic_permute(k,i)=min(niic_permute(k,i),qiic_permute(k,i)*1.e20_r8)
!
!            lami(k) = (cons1*ci*niic_permute(k,i)/qiic_permute(k,i))**(1._r8/di)
!            n0i(k) = niic_permute(k,i)*lami(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lami(k).lt.lammini) then
!
!               lami(k) = lammini
!               n0i(k) = lami(k)**(di+1._r8)*qiic_permute(k,i)/(ci*cons1)
!               niic_permute(k,i) = n0i(k)/lami(k)
!            else if (lami(k).gt.lammaxi) then
!               lami(k) = lammaxi
!               n0i(k) = lami(k)**(di+1._r8)*qiic_permute(k,i)/(ci*cons1)
!               niic_permute(k,i) = n0i(k)/lami(k)
!            end if
!
!         else
!            lami(k) = 0._r8
!            n0i(k) = 0._r8
!         end if
!
!         if (qcic_permute(k,i).ge.qsmall) then
!
!            ! add upper limit to in-cloud number concentration to prevent numerical error
!            ncic_permute(k,i)=min(ncic_permute(k,i),qcic_permute(k,i)*1.e20_r8)
!
!            ncic_permute(k,i)=max(ncic_permute(k,i),cdnl/rho_permute(k,i)) ! sghan minimum in #/cm  
!
!            ! get pgam from fit to observations of martin et al. 1994
!
!            pgam(k)=0.0005714_r8*(ncic_permute(k,i)/1.e6_r8*rho_permute(k,i))+0.2714_r8
!            pgam(k)=1._r8/(pgam(k)**2)-1._r8
!            pgam(k)=max(pgam(k),2._r8)
!            pgam(k)=min(pgam(k),15._r8)
!
!            ! calculate lamc
!
!            lamc(k) = (pi/6._r8*rhow*ncic_permute(k,i)*gamma(pgam(k)+4._r8)/ &
!                 (qcic_permute(k,i)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!
!            ! lammin, 50 micron diameter max mean size
!
!            lammin = (pgam(k)+1._r8)/50.e-6_r8
!            lammax = (pgam(k)+1._r8)/2.e-6_r8
!
!            if (lamc(k).lt.lammin) then
!               lamc(k) = lammin
!               ncic_permute(k,i) = 6._r8*lamc(k)**3*qcic_permute(k,i)* &
!                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
!            else if (lamc(k).gt.lammax) then
!               lamc(k) = lammax
!               ncic_permute(k,i) = 6._r8*lamc(k)**3*qcic_permute(k,i)* &
!                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
!            end if
!
!            ! parameter to calculate droplet freezing
!
!            cdist1(k) = ncic_permute(k,i)/gamma(pgam(k)+1._r8) 
!
!         else
!            lamc(k) = 0._r8
!            cdist1(k) = 0._r8
!         end if
!
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! begin micropysical process calculations 
!         !.................................................................
!         ! autoconversion of cloud liquid water to rain
!         ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
!         ! minimum qc of 1 x 10^-8 prevents floating point error
!
!         if (qcic_permute(k,i).ge.1.e-8_r8) then
!
!            ! nprc is increase in rain number conc due to autoconversion
!            ! nprc1 is decrease in cloud droplet conc due to autoconversion
!
!            ! assume exponential sub-grid distribution of qc, resulting in additional
!            ! factor related to qcvar below
!
!            ! hm switch for sub-columns, don't include sub-grid qc
!            if (microp_uniform) then
!
!               prc(k) = 1350._r8*qcic_permute(k,i)**2.47_r8* &
!                    (ncic_permute(k,i)/1.e6_r8*rho_permute(k,i))**(-1.79_r8)
!               nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
!               nprc1(k) = prc(k)/(qcic_permute(k,i)/ncic_permute(k,i))
!
!            else
!
!               prc(k) = cons2/(cons3*cons18)*1350._r8*qcic_permute(k,i)**2.47_r8* &
!                    (ncic_permute(k,i)/1.e6_r8*rho_permute(k,i))**(-1.79_r8)
!               nprc(k) = prc(k)/cons22
!               nprc1(k) = prc(k)/(qcic_permute(k,i)/ncic_permute(k,i))
!
!            end if               ! sub-column switch
!
!         else
!            prc(k)=0._r8
!            nprc(k)=0._r8
!            nprc1(k)=0._r8
!         end if
!
!         ! add autoconversion to precip from above to get provisional rain mixing ratio
!         ! and number concentration (qric and nric)
!
!         ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)
!
!         dum=0.45_r8
!         dum1=0.45_r8
!
!         if (k.eq.top_lev) then
!            qric_permute(k,i)=prc(k)*lcldm_permute(k,i)*dz_permute(k,i)/cldmax_permute(k,i)/dum
!            nric_permute(k,i)=nprc(k)*lcldm_permute(k,i)*dz_permute(k,i)/cldmax_permute(k,i)/dum
!         else
!            if (qric_permute(k-1,i).ge.qsmall) then
!               dum=umr(k-1)
!               dum1=unr(k-1)
!            end if
!
!            ! no autoconversion of rain number if rain/snow falling from above
!            ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
!            ! by the existing rain/snow particles from above
!
!            if (qric_permute(k-1,i).ge.1.e-9_r8.or.qniic_permute(k-1,i).ge.1.e-9_r8) then
!               nprc(k)=0._r8
!            end if
!
!            qric_permute(k,i) = (rho_permute(k-1,i)*umr(k-1)*qric_permute(k-1,i)*cldmax_permute(k-1,i)+ &
!                 (rho_permute(k,i)*dz_permute(k,i)*((pra(k-1)+prc(k))*lcldm_permute(k,i)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax_permute(k,i))))&
!                 /(dum*rho_permute(k,i)*cldmax_permute(k,i))
!            nric_permute(k,i) = (rho_permute(k-1,i)*unr(k-1)*nric_permute(k-1,i)*cldmax_permute(k-1,i)+ &
!                 (rho_permute(k,i)*dz_permute(k,i)*(nprc(k)*lcldm_permute(k,i)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax_permute(k,i))))&
!                 /(dum1*rho_permute(k,i)*cldmax_permute(k,i))
!
!         end if
!
!         !.......................................................................
!         ! Autoconversion of cloud ice to snow
!         ! similar to Ferrier (1994)
!
!         if (do_cldice) then
!            if (t_permute(k,i).le.273.15_r8.and.qiic_permute(k,i).ge.qsmall) then
!
!               ! note: assumes autoconversion timescale of 180 sec
!               
!               nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)
!
!               prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
!                    (cons23/lami(k)+3._r8*cons24/lami(k)**2+ &
!                    6._r8*dcs/lami(k)**3+6._r8/lami(k)**4)*exp(-lami(k)*dcs)
!            else
!               prci(k)=0._r8
!               nprci(k)=0._r8
!            end if
!         else
!            ! Add in the particles that we have already converted to snow, and
!            ! don't do any further autoconversion of ice.
!            prci(k)  = tnd_qsnow_permute(k,i) / cldm_permute(k,i)
!            nprci(k) = tnd_nsnow_permute(k,i) / cldm_permute(k,i)
!         end if
!
!         ! add autoconversion to flux from level above to get provisional snow mixing ratio
!         ! and number concentration (qniic and nsic)
!
!         dum=(asn_permute(k,i)*cons25)
!         dum1=(asn_permute(k,i)*cons25)
!
!         if (k.eq.top_lev) then
!            qniic_permute(k,i)=prci(k)*icldm_permute(k,i)*dz_permute(k,i)/cldmax_permute(k,i)/dum
!            nsic_permute(k,i)=nprci(k)*icldm_permute(k,i)*dz_permute(k,i)/cldmax_permute(k,i)/dum
!         else
!            if (qniic_permute(k-1,i).ge.qsmall) then
!               dum=ums(k-1)
!               dum1=uns(k-1)
!            end if
!
!            qniic_permute(k,i) = (rho_permute(k-1,i)*ums(k-1)*qniic_permute(k-1,i)*cldmax_permute(k-1,i)+ &
!                 (rho_permute(k,i)*dz_permute(k,i)*((prci(k)+prai(k-1)+psacws(k-1)+bergs(k-1))*icldm_permute(k,i)+(prds(k-1)+ &
!                 pracs(k-1)+mnuccr(k-1))*cldmax_permute(k,i))))&
!                 /(dum*rho_permute(k,i)*cldmax_permute(k,i))
!
!            nsic_permute(k,i) = (rho_permute(k-1,i)*uns(k-1)*nsic_permute(k-1,i)*cldmax_permute(k-1,i)+ &
!                 (rho_permute(k,i)*dz_permute(k,i)*(nprci(k)*icldm_permute(k,i)+(nsubs(k-1)+nsagg(k-1)+nnuccr(k-1))*cldmax_permute(k,i))))&
!                 /(dum1*rho_permute(k,i)*cldmax_permute(k,i))
!
!         end if
!
!         ! if precip mix ratio is zero so should number concentration
!
!         if (qniic_permute(k,i).lt.qsmall) then
!            qniic_permute(k,i)=0._r8
!            nsic_permute(k,i)=0._r8
!         end if
!
!         if (qric_permute(k,i).lt.qsmall) then
!            qric_permute(k,i)=0._r8
!            nric_permute(k,i)=0._r8
!         end if
!
!         ! make sure number concentration is a positive number to avoid 
!         ! taking root of negative later
!
!         nric_permute(k,i)=max(nric_permute(k,i),0._r8)
!         nsic_permute(k,i)=max(nsic_permute(k,i),0._r8)
!
!         !.......................................................................
!         ! get size distribution parameters for precip
!         !......................................................................
!         ! rain
!
!         if (qric_permute(k,i).ge.qsmall) then
!            lamr(k) = (pi*rhow*nric_permute(k,i)/qric_permute(k,i))**(1._r8/3._r8)
!            n0r(k) = nric_permute(k,i)*lamr(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lamr(k).lt.lamminr) then
!
!               lamr(k) = lamminr
!
!               n0r(k) = lamr(k)**4*qric_permute(k,i)/(pi*rhow)
!               nric_permute(k,i) = n0r(k)/lamr(k)
!            else if (lamr(k).gt.lammaxr) then
!               lamr(k) = lammaxr
!               n0r(k) = lamr(k)**4*qric_permute(k,i)/(pi*rhow)
!               nric_permute(k,i) = n0r(k)/lamr(k)
!            end if
!
!            ! provisional rain number and mass weighted mean fallspeed (m/s)
!
!            unr(k) = min(arn_permute(k,i)*cons4/lamr(k)**br,9.1_r8*rhof_permute(k,i))
!            umr(k) = min(arn_permute(k,i)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof_permute(k,i))
!
!         else
!            lamr(k) = 0._r8
!            n0r(k) = 0._r8
!            umr(k) = 0._r8
!            unr(k) = 0._r8
!         end if
!
!         !......................................................................
!         ! snow
!
!         if (qniic_permute(k,i).ge.qsmall) then
!            lams(k) = (cons6*cs*nsic_permute(k,i)/qniic_permute(k,i))**(1._r8/ds)
!            n0s(k) = nsic_permute(k,i)*lams(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lams(k).lt.lammins) then
!               lams(k) = lammins
!               n0s(k) = lams(k)**(ds+1._r8)*qniic_permute(k,i)/(cs*cons6)
!               nsic_permute(k,i) = n0s(k)/lams(k)
!
!            else if (lams(k).gt.lammaxs) then
!               lams(k) = lammaxs
!               n0s(k) = lams(k)**(ds+1._r8)*qniic_permute(k,i)/(cs*cons6)
!               nsic_permute(k,i) = n0s(k)/lams(k)
!            end if
!
!            ! provisional snow number and mass weighted mean fallspeed (m/s)
!
!            ums(k) = min(asn_permute(k,i)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof_permute(k,i))
!            uns(k) = min(asn_permute(k,i)*cons7/lams(k)**bs,1.2_r8*rhof_permute(k,i))
!
!         else
!            lams(k) = 0._r8
!            n0s(k) = 0._r8
!            ums(k) = 0._r8
!            uns(k) = 0._r8
!         end if
!
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!         ! heterogeneous freezing of cloud water
!
!         if (.not. use_hetfrz_classnuc) then
!
!            if (do_cldice .and. qcic_permute(k,i).ge.qsmall .and. t_permute(k,i).lt.269.15_r8) then
!
!               ! immersion freezing (Bigg, 1953)
!
!
!               ! subcolumns
!
!               if (microp_uniform) then
!
!                  mnuccc(k) = &
!                     pi*pi/36._r8*rhow* &
!                     cdist1(k)*gamma(7._r8+pgam(k))* &
!                     bimm*(exp(aimm*(273.15_r8-t_permute(k,i)))-1._r8)/ &
!                     lamc(k)**3/lamc(k)**3
!
!                  nnuccc(k) = &
!                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
!                     *bimm* &
!                     (exp(aimm*(273.15_r8-t_permute(k,i)))-1._r8)/lamc(k)**3
!
!               else
!
!                  mnuccc(k) = cons9/(cons3*cons19)* &
!                     pi*pi/36._r8*rhow* &
!                     cdist1(k)*gamma(7._r8+pgam(k))* &
!                     bimm*(exp(aimm*(273.15_r8-t_permute(k,i)))-1._r8)/ &
!                     lamc(k)**3/lamc(k)**3
!
!                  nnuccc(k) = cons10/(cons3*qcvar)* &
!                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
!                     *bimm* &
!                     (exp(aimm*(273.15_r8-t_permute(k,i)))-1._r8)/lamc(k)**3
!               end if           ! sub-columns
!
!
!               ! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
!               ! dust size and number in 4 bins are read in from companion routine
!
!               tcnt=(270.16_r8-t_permute(k,i))**1.3_r8
!               viscosity=1.8e-5_r8*(t_permute(k,i)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
!               mfp=2.0_r8*viscosity/(p_permute(k,i)  &                   ! Mean free path (m)
!                  *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t_permute(k,i))))           
!
!               nslip1=1.0_r8+(mfp/rndst_permute(k,i,1))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst_permute(k,i,1)/mfp))))! Slip correction factor
!               nslip2=1.0_r8+(mfp/rndst_permute(k,i,2))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst_permute(k,i,2)/mfp))))
!               nslip3=1.0_r8+(mfp/rndst_permute(k,i,3))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst_permute(k,i,3)/mfp))))
!               nslip4=1.0_r8+(mfp/rndst_permute(k,i,4))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst_permute(k,i,4)/mfp))))
!
!               ndfaer1=1.381e-23_r8*t_permute(k,i)*nslip1/(6._r8*pi*viscosity*rndst_permute(k,i,1))  ! aerosol diffusivity (m2/s)
!               ndfaer2=1.381e-23_r8*t_permute(k,i)*nslip2/(6._r8*pi*viscosity*rndst_permute(k,i,2))
!               ndfaer3=1.381e-23_r8*t_permute(k,i)*nslip3/(6._r8*pi*viscosity*rndst_permute(k,i,3))
!               ndfaer4=1.381e-23_r8*t_permute(k,i)*nslip4/(6._r8*pi*viscosity*rndst_permute(k,i,4))
!
!
!               if (microp_uniform) then
!
!                  mnucct(k) = &
!                     (ndfaer1*(nacon_permute(k,i,1)*tcnt)+ndfaer2*(nacon_permute(k,i,2)*tcnt)+ &
!                     ndfaer3*(nacon_permute(k,i,3)*tcnt)+ndfaer4*(nacon_permute(k,i,4)*tcnt))*pi*pi/3._r8*rhow* &
!                     cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
!
!                  nnucct(k) = (ndfaer1*(nacon_permute(k,i,1)*tcnt)+ndfaer2*(nacon_permute(k,i,2)*tcnt)+ &
!                     ndfaer3*(nacon_permute(k,i,3)*tcnt)+ndfaer4*(nacon_permute(k,i,4)*tcnt))*2._r8*pi*  &
!                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
!
!               else
!
!                  mnucct(k) = gamma(qcvar+4._r8/3._r8)/(cons3*qcvar**(4._r8/3._r8))*  &
!                     (ndfaer1*(nacon_permute(k,i,1)*tcnt)+ndfaer2*(nacon_permute(k,i,2)*tcnt)+ &
!                     ndfaer3*(nacon_permute(k,i,3)*tcnt)+ndfaer4*(nacon_permute(k,i,4)*tcnt))*pi*pi/3._r8*rhow* &
!                     cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
!
!                  nnucct(k) =  gamma(qcvar+1._r8/3._r8)/(cons3*qcvar**(1._r8/3._r8))*  &
!                     (ndfaer1*(nacon_permute(k,i,1)*tcnt)+ndfaer2*(nacon_permute(k,i,2)*tcnt)+ &
!                     ndfaer3*(nacon_permute(k,i,3)*tcnt)+ndfaer4*(nacon_permute(k,i,4)*tcnt))*2._r8*pi*  &
!                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
!
!               end if      ! sub-column switch
!
!               ! make sure number of droplets frozen does not exceed available ice nuclei concentration
!               ! this prevents 'runaway' droplet freezing
!
!               if (nnuccc(k)*lcldm_permute(k,i).gt.nnuccd(k)) then
!                  dum=(nnuccd(k)/(nnuccc(k)*lcldm_permute(k,i)))
!                  ! scale mixing ratio of droplet freezing with limit
!                  mnuccc(k)=mnuccc(k)*dum
!                  nnuccc(k)=nnuccd(k)/lcldm_permute(k,i)
!               end if
!
!            else
!               mnuccc(k)=0._r8
!               nnuccc(k)=0._r8
!               mnucct(k)=0._r8
!               nnucct(k)=0._r8
!            end if
!
!         else
!            if (do_cldice .and. qcic_permute(k,i) >= qsmall) then
!               con1 = 1._r8/(1.333_r8*pi)**0.333_r8
!               r3lx = con1*(rho_permute(k,i)*qcic_permute(k,i)/(rhow*max(ncic_permute(k,i)*rho_permute(k,i), 1.0e6_r8)))**0.333_r8 ! in m
!               r3lx = max(4.e-6_r8, r3lx)
!               mi0l = 4._r8/3._r8*pi*rhow*r3lx**3_r8
!                
!               nnuccc(k) = frzimm_permute(k,i)*1.0e6_r8/rho_permute(k,i)
!               mnuccc(k) = nnuccc(k)*mi0l 
!
!               nnucct(k) = frzcnt_permute(k,i)*1.0e6_r8/rho_permute(k,i)
!               mnucct(k) = nnucct(k)*mi0l 
!
!               nnudep(k) = frzdep_permute(k,i)*1.0e6_r8/rho_permute(k,i)
!               mnudep(k) = nnudep(k)*mi0
!            else
!               nnuccc(k) = 0._r8
!               mnuccc(k) = 0._r8
!
!               nnucct(k) = 0._r8
!               mnucct(k) = 0._r8
!
!               nnudep(k) = 0._r8
!               mnudep(k) = 0._r8
!            end if
!         endif
!
!
!         !.......................................................................
!         ! snow self-aggregation from passarelli, 1978, used by reisner, 1998
!         ! this is hard-wired for bs = 0.4 for now
!         ! ignore self-collection of cloud ice
!
!         if (qniic_permute(k,i).ge.qsmall .and. t_permute(k,i).le.273.15_r8) then
!            nsagg(k) = -1108._r8*asn_permute(k,i)*Eii* &
!                 pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)*rho_permute(k,i)** &
!                 ((2._r8+bs)/3._r8)*qniic_permute(k,i)**((2._r8+bs)/3._r8)* &
!                 (nsic_permute(k,i)*rho_permute(k,i))**((4._r8-bs)/3._r8)/ &
!                 (4._r8*720._r8*rho_permute(k,i))
!         else
!            nsagg(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! accretion of cloud droplets onto snow/graupel
!         ! here use continuous collection equation with
!         ! simple gravitational collection kernel
!         ! ignore collisions between droplets/cloud ice
!         ! since minimum size ice particle for accretion is 50 - 150 micron
!
!         ! ignore collision of snow with droplets above freezing
!
!         if (qniic_permute(k,i).ge.qsmall .and. t_permute(k,i).le.tmelt .and. &
!              qcic_permute(k,i).ge.qsmall) then
!
!            ! put in size dependent collection efficiency
!            ! mean diameter of snow is area-weighted, since
!            ! accretion is function of crystal geometric area
!            ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)
!
!            dc0 = (pgam(k)+1._r8)/lamc(k)
!            ds0 = 1._r8/lams(k)
!            dum = dc0*dc0*uns(k)*rhow/(9._r8*mu_permute(k,i)*ds0)
!            eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))
!
!            eci = max(eci,0._r8)
!            eci = min(eci,1._r8)
!
!
!            ! no impact of sub-grid distribution of qc since psacws
!            ! is linear in qc
!
!            psacws(k) = pi/4._r8*asn_permute(k,i)*qcic_permute(k,i)*rho_permute(k,i)* &
!                 n0s(k)*Eci*cons11/ &
!                 lams(k)**(bs+3._r8)
!            npsacws(k) = pi/4._r8*asn_permute(k,i)*ncic_permute(k,i)*rho_permute(k,i)* &
!                 n0s(k)*Eci*cons11/ &
!                 lams(k)**(bs+3._r8)
!         else
!            psacws(k)=0._r8
!            npsacws(k)=0._r8
!         end if
!
!         ! add secondary ice production due to accretion of droplets by snow 
!         ! (Hallet-Mossop process) (from Cotton et al., 1986)
!
!         if (.not. do_cldice) then
!            ni_secp   = 0.0_r8
!            nsacwi(k) = 0.0_r8
!            msacwi(k) = 0.0_r8
!         else if((t_permute(k,i).lt.270.16_r8) .and. (t_permute(k,i).ge.268.16_r8)) then
!            ni_secp   = 3.5e8_r8*(270.16_r8-t_permute(k,i))/2.0_r8*psacws(k)
!            nsacwi(k) = ni_secp
!            msacwi(k) = min(ni_secp*mi0,psacws(k))
!         else if((t_permute(k,i).lt.268.16_r8) .and. (t_permute(k,i).ge.265.16_r8)) then
!            ni_secp   = 3.5e8_r8*(t_permute(k,i)-265.16_r8)/3.0_r8*psacws(k)
!            nsacwi(k) = ni_secp
!            msacwi(k) = min(ni_secp*mi0,psacws(k))
!         else
!            ni_secp   = 0.0_r8
!            nsacwi(k) = 0.0_r8
!            msacwi(k) = 0.0_r8
!         endif
!         psacws(k) = max(0.0_r8,psacws(k)-ni_secp*mi0)
!
!         !.......................................................................
!         ! accretion of rain water by snow
!         ! formula from ikawa and saito, 1991, used by reisner et al., 1998
!
!         if (qric_permute(k,i).ge.1.e-8_r8 .and. qniic_permute(k,i).ge.1.e-8_r8 .and. & 
!              t_permute(k,i).le.273.15_r8) then
!
!            pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))**2+ &
!                 0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho_permute(k,i)* &
!                 n0r(k)*n0s(k)* &
!                 (5._r8/(lamr(k)**6*lams(k))+ &
!                 2._r8/(lamr(k)**5*lams(k)**2)+ &
!                 0.5_r8/(lamr(k)**4*lams(k)**3)))
!
!            npracs(k) = pi/2._r8*rho_permute(k,i)*ecr*(1.7_r8*(unr(k)-uns(k))**2+ &
!                 0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
!                 (1._r8/(lamr(k)**3*lams(k))+ &
!                 1._r8/(lamr(k)**2*lams(k)**2)+ &
!                 1._r8/(lamr(k)*lams(k)**3))
!
!         else
!            pracs(k)=0._r8
!            npracs(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! heterogeneous freezing of rain drops
!         ! follows from Bigg (1953)
!
!         if (t_permute(k,i).lt.269.15_r8 .and. qric_permute(k,i).ge.qsmall) then
!
!            mnuccr(k) = 20._r8*pi*pi*rhow*nric_permute(k,i)*bimm* &
!                 (exp(aimm*(273.15_r8-t_permute(k,i)))-1._r8)/lamr(k)**3 &
!                 /lamr(k)**3
!
!            nnuccr(k) = pi*nric_permute(k,i)*bimm* &
!                 (exp(aimm*(273.15_r8-t_permute(k,i)))-1._r8)/lamr(k)**3
!         else
!            mnuccr(k)=0._r8
!            nnuccr(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! accretion of cloud liquid water by rain
!         ! formula from Khrouditnov and Kogan (2000)
!         ! gravitational collection kernel, droplet fall speed neglected
!
!         if (qric_permute(k,i).ge.qsmall .and. qcic_permute(k,i).ge.qsmall) then
!
!            ! include sub-grid distribution of cloud water
!
!            ! add sub-column switch
!
!            if (microp_uniform) then
!
!               pra(k) = 67._r8*(qcic_permute(k,i)*qric_permute(k,i))**1.15_r8
!               npra(k) = pra(k)/(qcic_permute(k,i)/ncic_permute(k,i))
!
!            else
!
!               pra(k) = accre_enhan_permute(k,i)*(cons12/(cons3*cons20)*67._r8*(qcic_permute(k,i)*qric_permute(k,i))**1.15_r8)
!               npra(k) = pra(k)/(qcic_permute(k,i)/ncic_permute(k,i))
!
!            end if               ! sub-column switch
!
!         else
!            pra(k)=0._r8
!            npra(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! Self-collection of rain drops
!         ! from Beheng(1994)
!
!         if (qric_permute(k,i).ge.qsmall) then
!            nragg(k) = -8._r8*nric_permute(k,i)*qric_permute(k,i)*rho_permute(k,i)
!         else
!            nragg(k)=0._r8
!         end if
!
!         !.......................................................................
!         ! Accretion of cloud ice by snow
!         ! For this calculation, it is assumed that the Vs >> Vi
!         ! and Ds >> Di for continuous collection
!
!         if (do_cldice .and. qniic_permute(k,i).ge.qsmall.and.qiic_permute(k,i).ge.qsmall &
!              .and.t_permute(k,i).le.273.15_r8) then
!
!            prai(k) = pi/4._r8*asn_permute(k,i)*qiic_permute(k,i)*rho_permute(k,i)* &
!                 n0s(k)*Eii*cons11/ &
!                 lams(k)**(bs+3._r8)
!            nprai(k) = pi/4._r8*asn_permute(k,i)*niic_permute(k,i)* &
!                 rho_permute(k,i)*n0s(k)*Eii*cons11/ &
!                 lams(k)**(bs+3._r8)
!         else
!            prai(k)=0._r8
!            nprai(k)=0._r8
!         end if
!
!         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! calculate evaporation/sublimation of rain and snow
!         ! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
!         ! in-cloud condensation/deposition of rain and snow is neglected
!         ! except for transfer of cloud water to snow through bergeron process
!
!         ! initialize evap/sub tendncies
!         pre(k)=0._r8
!         prds(k)=0._r8
!
!         ! evaporation of rain
!         ! only calculate if there is some precip fraction > cloud fraction
!
!         if (qcic_permute(k,i)+qiic_permute(k,i).lt.1.e-6_r8.or.cldmax_permute(k,i).gt.lcldm_permute(k,i)) then
!
!            ! set temporary cloud fraction to zero if cloud water + ice is very small
!            ! this will ensure that evaporation/sublimation of precip occurs over
!            ! entire grid cell, since min cloud fraction is specified otherwise
!            if (qcic_permute(k,i)+qiic_permute(k,i).lt.1.e-6_r8) then
!               dum=0._r8
!            else
!               dum=lcldm_permute(k,i)
!            end if
!
!            ! saturation vapor pressure
!            esn=svp_water(t_permute(k,i))
!            qsn=svp_to_qsat(esn, p_permute(k,i))
!
!            ! recalculate saturation vapor pressure for liquid and ice
!            esl_permute(k,i)=esn
!            esi_permute(k,i)=svp_ice(t_permute(k,i))
!            ! hm fix, make sure when above freezing that esi=esl, not active yet
!            if (t_permute(k,i).gt.tmelt)esi_permute(k,i)=esl_permute(k,i)
!
!            ! calculate q for out-of-cloud region
!            qclr=(q_permute(k,i)-dum*qsn)/(1._r8-dum)
!
!            if (qric_permute(k,i).ge.qsmall) then
!
!               qvs=svp_to_qsat(esl_permute(k,i), p_permute(k,i))
!               dqsdt = xxlv*qvs/(rv*t_permute(k,i)**2)
!               ab = 1._r8+dqsdt*xxlv/cpp
!               epsr = 2._r8*pi*n0r(k)*rho_permute(k,i)*dv_permute(k,i)* &
!                    (f1r/(lamr(k)*lamr(k))+ &
!                    f2r*(arn_permute(k,i)*rho_permute(k,i)/mu_permute(k,i))**0.5_r8* &
!                    sc_permute(k,i)**(1._r8/3._r8)*cons13/ &
!                    (lamr(k)**(5._r8/2._r8+br/2._r8)))
!
!               pre(k) = epsr*(qclr-qvs)/ab
!
!               ! only evaporate in out-of-cloud region
!               ! and distribute across cldmax
!               pre(k)=min(pre(k)*(cldmax_permute(k,i)-dum),0._r8)
!               pre(k)=pre(k)/cldmax_permute(k,i)
!               am_evp_st_permute(k,i) = max(cldmax_permute(k,i)-dum, 0._r8)
!            end if
!
!            ! sublimation of snow
!            if (qniic_permute(k,i).ge.qsmall) then
!               qvi=svp_to_qsat(esi_permute(k,i), p_permute(k,i))
!               dqsidt =  xxls*qvi/(rv*t_permute(k,i)**2)
!               abi = 1._r8+dqsidt*xxls/cpp
!               epss = 2._r8*pi*n0s(k)*rho_permute(k,i)*dv_permute(k,i)* &
!                    (f1s/(lams(k)*lams(k))+ &
!                    f2s*(asn_permute(k,i)*rho_permute(k,i)/mu_permute(k,i))**0.5_r8* &
!                    sc_permute(k,i)**(1._r8/3._r8)*cons14/ &
!                    (lams(k)**(5._r8/2._r8+bs/2._r8)))
!               prds(k) = epss*(qclr-qvi)/abi
!
!               ! only sublimate in out-of-cloud region and distribute over cldmax
!               prds(k)=min(prds(k)*(cldmax_permute(k,i)-dum),0._r8)
!               prds(k)=prds(k)/cldmax_permute(k,i)
!               am_evp_st_permute(k,i) = max(cldmax_permute(k,i)-dum, 0._r8)
!            end if
!
!            ! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
!            ! get updated RH at end of time step based on cloud water/ice condensation/evap
!
!            qtmp=q_permute(k,i)-(cmei_permute(k,i)+(pre(k)+prds(k))*cldmax_permute(k,i))*deltat
!            ttmp=t_permute(k,i)+((pre(k)*cldmax_permute(k,i))*xxlv+ &
!                 (cmei_permute(k,i)+prds(k)*cldmax_permute(k,i))*xxls)*deltat/cpp
!
!            !limit range of temperatures!
!            ttmp=max(180._r8,min(ttmp,323._r8))
!
!            esn=svp_water(ttmp)  ! use rhw to allow ice supersaturation
!            qsn=svp_to_qsat(esn, p_permute(k,i))
!
!            ! modify precip evaporation rate if q > qsat
!            if (qtmp.gt.qsn) then
!               if (pre(k)+prds(k).lt.-1.e-20_r8) then
!                  dum1=pre(k)/(pre(k)+prds(k))
!                  ! recalculate q and t after cloud water cond but without precip evap
!                  qtmp=q_permute(k,i)-(cmei_permute(k,i))*deltat
!                  ttmp=t_permute(k,i)+(cmei_permute(k,i)*xxls)*deltat/cpp
!                  esn=svp_water(ttmp) ! use rhw to allow ice supersaturation
!                  qsn=svp_to_qsat(esn, p_permute(k,i))
!                  dum=(qtmp-qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))
!                  dum=min(dum,0._r8)
!
!                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
!                  pre(k)=dum*dum1/deltat/cldmax_permute(k,i)
!
!                  ! do separately using RHI for prds....
!                  esn=svp_ice(ttmp) ! use rhi to allow ice supersaturation
!                  qsn=svp_to_qsat(esn, p_permute(k,i))
!                  dum=(qtmp-qsn)/(1._r8 + cons28*qsn/(cpp*rv*ttmp**2))
!                  dum=min(dum,0._r8)
!
!                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
!                  prds(k)=dum*(1._r8-dum1)/deltat/cldmax_permute(k,i)
!               end if
!            end if
!         end if
!
!         ! bergeron process - evaporation of droplets and deposition onto snow
!
!         if (qniic_permute(k,i).ge.qsmall.and.qcic_permute(k,i).ge.qsmall.and.t_permute(k,i).lt.tmelt) then
!            qvi=svp_to_qsat(esi_permute(k,i), p_permute(k,i))
!            qvs=svp_to_qsat(esl_permute(k,i), p_permute(k,i))
!            dqsidt =  xxls*qvi/(rv*t_permute(k,i)**2)
!            abi = 1._r8+dqsidt*xxls/cpp
!            epss = 2._r8*pi*n0s(k)*rho_permute(k,i)*dv_permute(k,i)* &
!                 (f1s/(lams(k)*lams(k))+ &
!                 f2s*(asn_permute(k,i)*rho_permute(k,i)/mu_permute(k,i))**0.5_r8* &
!                 sc_permute(k,i)**(1._r8/3._r8)*cons14/ &
!                 (lams(k)**(5._r8/2._r8+bs/2._r8)))
!            bergs(k)=epss*(qvs-qvi)/abi
!         else
!            bergs(k)=0._r8
!         end if
!
!         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ! conservation to ensure no negative values of cloud water/precipitation
!         ! in case microphysical process rates are large
!
!         ! make sure and use end-of-time step values for cloud water, ice, due
!         ! condensation/deposition
!
!         ! note: for check on conservation, processes are multiplied by omsm
!         ! to prevent problems due to round off error
!
!         ! include mixing timescale  (mtime)
!
!         qce=(qc_permute(k,i) - berg_permute(k,i)*deltat)
!         nce=(nc_permute(k,i)+npccn(k)*deltat*mtime)
!         qie=(qi_permute(k,i)+(cmei_permute(k,i)+berg_permute(k,i))*deltat)
!         nie=(ni_permute(k,i)+nnuccd(k)*deltat*mtime)
!
!         ! conservation of qc
!
!         dum = (prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+ &
!              psacws(k)+bergs(k))*lcldm_permute(k,i)*deltat
!
!         if (dum.gt.qce) then
!            ratio = qce/deltat/lcldm_permute(k,i)/(prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+psacws(k)+bergs(k))*omsm 
!
!            prc(k) = prc(k)*ratio
!            pra(k) = pra(k)*ratio
!            mnuccc(k) = mnuccc(k)*ratio
!            mnucct(k) = mnucct(k)*ratio  
!            msacwi(k) = msacwi(k)*ratio  
!            psacws(k) = psacws(k)*ratio
!            bergs(k) = bergs(k)*ratio
!         end if
!
!         ! conservation of nc
!
!         dum = (nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
!              npsacws(k)-nsubc(k))*lcldm_permute(k,i)*deltat
!
!         if (dum.gt.nce) then
!            ratio = nce/deltat/((nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+&
!                 npsacws(k)-nsubc(k))*lcldm_permute(k,i))*omsm
!
!            nprc1(k) = nprc1(k)*ratio
!            npra(k) = npra(k)*ratio
!            nnuccc(k) = nnuccc(k)*ratio
!            nnucct(k) = nnucct(k)*ratio  
!            npsacws(k) = npsacws(k)*ratio
!            nsubc(k)=nsubc(k)*ratio
!         end if
!
!         ! conservation of qi
!
!         if (do_cldice) then
!
!            frztmp = -mnuccc(k) - mnucct(k) - msacwi(k)
!            if (use_hetfrz_classnuc) frztmp = -mnuccc(k)-mnucct(k)-mnudep(k)-msacwi(k)
!            dum = ( frztmp*lcldm_permute(k,i) + (prci(k)+prai(k))*icldm_permute(k,i) )*deltat
!
!            if (dum.gt.qie) then
!
!               frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
!               if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
!               ratio = (qie/deltat + frztmp*lcldm_permute(k,i))/((prci(k)+prai(k))*icldm_permute(k,i))*omsm 
!               prci(k) = prci(k)*ratio
!               prai(k) = prai(k)*ratio
!            end if
!
!            ! conservation of ni
!            frztmp = -nnucct(k) - nsacwi(k)
!            if (use_hetfrz_classnuc) frztmp = -nnucct(k) - nnuccc(k) - nnudep(k) - nsacwi(k)
!            dum = ( frztmp*lcldm_permute(k,i) + (nprci(k)+nprai(k)-nsubi(k))*icldm_permute(k,i) )*deltat
!
!            if (dum.gt.nie) then
!
!               frztmp = nnucct(k) + nsacwi(k)
!               if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
!               ratio = (nie/deltat + frztmp*lcldm_permute(k,i))/ &  
!                     ((nprci(k)+nprai(k)-nsubi(k))*icldm_permute(k,i))*omsm
!               nprci(k) = nprci(k)*ratio
!               nprai(k) = nprai(k)*ratio
!               nsubi(k) = nsubi(k)*ratio
!            end if
!         end if
!
!         ! for precipitation conservation, use logic that vertical integral 
!         ! of tendency from current level to top of model (i.e., qrtot) cannot be negative
!
!         ! conservation of rain mixing rat
!
!         if (((prc(k)+pra(k))*lcldm_permute(k,i)+(-mnuccr(k)+pre(k)-pracs(k))*&
!              cldmax_permute(k,i))*dz_permute(k,i)*rho_permute(k,i)+qrtot.lt.0._r8) then
!
!            if (-pre(k)+pracs(k)+mnuccr(k).ge.qsmall) then
!
!               ratio = (qrtot/(dz_permute(k,i)*rho_permute(k,i))+(prc(k)+pra(k))*lcldm_permute(k,i))/&
!                    ((-pre(k)+pracs(k)+mnuccr(k))*cldmax_permute(k,i))*omsm 
!
!               pre(k) = pre(k)*ratio
!               pracs(k) = pracs(k)*ratio
!               mnuccr(k) = mnuccr(k)*ratio
!            end if
!         end if
!
!         ! conservation of nr
!         ! for now neglect evaporation of nr
!         nsubr(k)=0._r8
!
!         if ((nprc(k)*lcldm_permute(k,i)+(-nnuccr(k)+nsubr(k)-npracs(k)&
!              +nragg(k))*cldmax_permute(k,i))*dz_permute(k,i)*rho_permute(k,i)+nrtot.lt.0._r8) then
!
!            if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).ge.qsmall) then
!
!               ratio = (nrtot/(dz_permute(k,i)*rho_permute(k,i))+nprc(k)*lcldm_permute(k,i))/&
!                    ((-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k))*cldmax_permute(k,i))*omsm
!
!               nsubr(k) = nsubr(k)*ratio
!               npracs(k) = npracs(k)*ratio
!               nnuccr(k) = nnuccr(k)*ratio
!               nragg(k) = nragg(k)*ratio
!            end if
!         end if
!
!         ! conservation of snow mix ratio
!
!         if (((bergs(k)+psacws(k))*lcldm_permute(k,i)+(prai(k)+prci(k))*icldm_permute(k,i)+(pracs(k)+&
!              mnuccr(k)+prds(k))*cldmax_permute(k,i))*dz_permute(k,i)*rho_permute(k,i)+qstot.lt.0._r8) then
!
!            if (-prds(k).ge.qsmall) then
!
!               ratio = (qstot/(dz_permute(k,i)*rho_permute(k,i))+(bergs(k)+psacws(k))*lcldm_permute(k,i)+(prai(k)+prci(k))*icldm_permute(k,i)+&
!                    (pracs(k)+mnuccr(k))*cldmax_permute(k,i))/(-prds(k)*cldmax_permute(k,i))*omsm
!
!               prds(k) = prds(k)*ratio
!            end if
!         end if
!
!         ! conservation of ns
!
!         ! calculate loss of number due to sublimation
!         ! for now neglect sublimation of ns
!         nsubs(k)=0._r8
!
!         if ((nprci(k)*icldm_permute(k,i)+(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax_permute(k,i))*&
!              dz_permute(k,i)*rho_permute(k,i)+nstot.lt.0._r8) then
!
!            if (-nsubs(k)-nsagg(k).ge.qsmall) then
!
!               ratio = (nstot/(dz_permute(k,i)*rho_permute(k,i))+nprci(k)*icldm_permute(k,i)+&
!                    nnuccr(k)*cldmax_permute(k,i))/((-nsubs(k)-nsagg(k))*cldmax_permute(k,i))*omsm
!
!               nsubs(k) = nsubs(k)*ratio
!               nsagg(k) = nsagg(k)*ratio
!            end if
!         end if
!
!         ! get tendencies due to microphysical conversion processes
!         ! note: tendencies are multiplied by appropaiate cloud/precip 
!         ! fraction to get grid-scale values
!         ! note: cmei is already grid-average values
!
!         qvlat_permute(k,i) = qvlat_permute(k,i)-(pre(k)+prds(k))*cldmax_permute(k,i)-cmei_permute(k,i) 
!
!         tlat_permute(k,i) = tlat_permute(k,i)+((pre(k)*cldmax_permute(k,i)) &
!              *xxlv+(prds(k)*cldmax_permute(k,i)+cmei_permute(k,i))*xxls+ &
!              ((bergs(k)+psacws(k)+mnuccc(k)+mnucct(k)+msacwi(k))*lcldm_permute(k,i)+(mnuccr(k)+ &
!              pracs(k))*cldmax_permute(k,i)+berg_permute(k,i))*xlf)
!
!         qctend_permute(k,i) = qctend_permute(k,i)+ &
!              (-pra(k)-prc(k)-mnuccc(k)-mnucct(k)-msacwi(k)- & 
!              psacws(k)-bergs(k))*lcldm_permute(k,i)-berg_permute(k,i)
!
!         if (do_cldice) then
!
!            frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
!            if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
!            qitend_permute(k,i) = qitend_permute(k,i) + frztmp*lcldm_permute(k,i) + &
!               (-prci(k)-prai(k))*icldm_permute(k,i) + cmei_permute(k,i) + berg_permute(k,i)
!
!         end if
!
!         qrtend_permute(k,i) = qrtend_permute(k,i)+ &
!              (pra(k)+prc(k))*lcldm_permute(k,i)+(pre(k)-pracs(k)- &
!              mnuccr(k))*cldmax_permute(k,i)
!
!         qnitend_permute(k,i) = qnitend_permute(k,i)+ &
!              (prai(k)+prci(k))*icldm_permute(k,i)+(psacws(k)+bergs(k))*lcldm_permute(k,i)+(prds(k)+ &
!              pracs(k)+mnuccr(k))*cldmax_permute(k,i)
!
!         ! add output for cmei (accumulate)
!         cmeiout_permute(k,i) = cmeiout_permute(k,i) + cmei_permute(k,i)
!
!         ! assign variables for trop_mozart, these are grid-average
!         ! evaporation/sublimation is stored here as positive term
!
!         evapsnow_permute(k,i) = evapsnow_permute(k,i)-prds(k)*cldmax_permute(k,i)
!         nevapr_permute(k,i) = nevapr_permute(k,i)-pre(k)*cldmax_permute(k,i)
!         nevapr2_permute(k,i) = nevapr2_permute(k,i)-pre(k)*cldmax_permute(k,i)
!
!         ! change to make sure prain is positive: do not remove snow from
!         ! prain used for wet deposition
!         prain_permute(k,i) = prain_permute(k,i)+(pra(k)+prc(k))*lcldm_permute(k,i)+(-pracs(k)- &
!              mnuccr(k))*cldmax_permute(k,i)
!         prodsnow_permute(k,i) = prodsnow_permute(k,i)+(prai(k)+prci(k))*icldm_permute(k,i)+(psacws(k)+bergs(k))*lcldm_permute(k,i)+(&
!              pracs(k)+mnuccr(k))*cldmax_permute(k,i)
!
!         ! following are used to calculate 1st order conversion rate of cloud water
!         !    to rain and snow (1/s), for later use in aerosol wet removal routine
!         ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
!         !    used to calculate pra, prc, ... in this routine
!         ! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
!         !                      (no cloud ice or bergeron terms)
!         ! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }
!
!         qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k) + (pra(k)+prc(k)+psacws(k))*lcldm_permute(k,i) 
!         qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc_permute(k,i) 
!
!         ! microphysics output, note this is grid-averaged
!         prao_permute(k,i)=prao_permute(k,i)+pra(k)*lcldm_permute(k,i)
!         prco_permute(k,i)=prco_permute(k,i)+prc(k)*lcldm_permute(k,i)
!         mnuccco_permute(k,i)=mnuccco_permute(k,i)+mnuccc(k)*lcldm_permute(k,i)
!         mnuccto_permute(k,i)=mnuccto_permute(k,i)+mnucct(k)*lcldm_permute(k,i)
!         mnuccdo_permute(k,i)=mnuccdo_permute(k,i)+mnuccd(k)*lcldm_permute(k,i)
!         msacwio_permute(k,i)=msacwio_permute(k,i)+msacwi(k)*lcldm_permute(k,i)
!         psacwso_permute(k,i)=psacwso_permute(k,i)+psacws(k)*lcldm_permute(k,i)
!         bergso_permute(k,i)=bergso_permute(k,i)+bergs(k)*lcldm_permute(k,i)
!         bergo_permute(k,i)=bergo_permute(k,i)+berg_permute(k,i)
!         prcio_permute(k,i)=prcio_permute(k,i)+prci(k)*icldm_permute(k,i)
!         praio_permute(k,i)=praio_permute(k,i)+prai(k)*icldm_permute(k,i)
!         mnuccro_permute(k,i)=mnuccro_permute(k,i)+mnuccr(k)*cldmax_permute(k,i)
!         pracso_permute(k,i)=pracso_permute(k,i)+pracs (k)*cldmax_permute(k,i)
!
!         ! multiply activation/nucleation by mtime to account for fast timescale
!
!         nctend_permute(k,i) = nctend_permute(k,i)+ npccn(k)*mtime+&
!              (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k) & 
!              -npra(k)-nprc1(k))*lcldm_permute(k,i)      
!
!         if (do_cldice) then
!
!            frztmp = nnucct(k) + nsacwi(k)
!            if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
!            nitend_permute(k,i) = nitend_permute(k,i) + nnuccd(k)*mtime + & 
!                  frztmp*lcldm_permute(k,i) + (nsubi(k)-nprci(k)-nprai(k))*icldm_permute(k,i)
!
!         end if
!
!         nstend_permute(k,i) = nstend_permute(k,i)+(nsubs(k)+ &
!              nsagg(k)+nnuccr(k))*cldmax_permute(k,i)+nprci(k)*icldm_permute(k,i)
!
!         nrtend_permute(k,i) = nrtend_permute(k,i)+ &
!              nprc(k)*lcldm_permute(k,i)+(nsubr(k)-npracs(k)-nnuccr(k) &
!              +nragg(k))*cldmax_permute(k,i)
!
!         ! make sure that nc and ni at advanced time step do not exceed
!         ! maximum (existing N + source terms*dt), which is possible due to
!         ! fast nucleation timescale
!
!         if (nctend_permute(k,i).gt.0._r8.and.nc_permute(k,i)+nctend_permute(k,i)*deltat.gt.ncmax) then
!            nctend_permute(k,i)=max(0._r8,(ncmax-nc_permute(k,i))/deltat)
!         end if
!
!         if (do_cldice .and. nitend_permute(k,i).gt.0._r8.and.ni_permute(k,i)+nitend_permute(k,i)*deltat.gt.nimax) then
!            nitend_permute(k,i)=max(0._r8,(nimax-ni_permute(k,i))/deltat)
!         end if
!
!         ! get final values for precipitation q and N, based on
!         ! flux of precip from above, source/sink term, and terminal fallspeed
!         ! see eq. 15-16 in MG2008
!
!         ! rain
!
!         if (qric_permute(k,i).ge.qsmall) then
!            if (k.eq.top_lev) then
!               qric_permute(k,i)=qrtend_permute(k,i)*dz_permute(k,i)/cldmax_permute(k,i)/umr(k)
!               nric_permute(k,i)=nrtend_permute(k,i)*dz_permute(k,i)/cldmax_permute(k,i)/unr(k)
!            else
!               qric_permute(k,i) = (rho_permute(k-1,i)*umr(k-1)*qric_permute(k-1,i)*cldmax_permute(k-1,i)+ &
!                    (rho_permute(k,i)*dz_permute(k,i)*qrtend_permute(k,i)))/(umr(k)*rho_permute(k,i)*cldmax_permute(k,i))
!               nric_permute(k,i) = (rho_permute(k-1,i)*unr(k-1)*nric_permute(k-1,i)*cldmax_permute(k-1,i)+ &
!                    (rho_permute(k,i)*dz_permute(k,i)*nrtend_permute(k,i)))/(unr(k)*rho_permute(k,i)*cldmax_permute(k,i))
!
!            end if
!         else
!            qric_permute(k,i)=0._r8
!            nric_permute(k,i)=0._r8
!         end if
!
!         ! snow
!
!         if (qniic_permute(k,i).ge.qsmall) then
!            if (k.eq.top_lev) then
!               qniic_permute(k,i)=qnitend_permute(k,i)*dz_permute(k,i)/cldmax_permute(k,i)/ums(k)
!               nsic_permute(k,i)=nstend_permute(k,i)*dz_permute(k,i)/cldmax_permute(k,i)/uns(k)
!            else
!               qniic_permute(k,i) = (rho_permute(k-1,i)*ums(k-1)*qniic_permute(k-1,i)*cldmax_permute(k-1,i)+ &
!                    (rho_permute(k,i)*dz_permute(k,i)*qnitend_permute(k,i)))/(ums(k)*rho_permute(k,i)*cldmax_permute(k,i))
!               nsic_permute(k,i) = (rho_permute(k-1,i)*uns(k-1)*nsic_permute(k-1,i)*cldmax_permute(k-1,i)+ &
!                    (rho_permute(k,i)*dz_permute(k,i)*nstend_permute(k,i)))/(uns(k)*rho_permute(k,i)*cldmax_permute(k,i))
!            end if
!         else
!            qniic_permute(k,i)=0._r8
!            nsic_permute(k,i)=0._r8
!         end if
!
!         ! calculate precipitation flux at surface
!         ! divide by density of water to get units of m/s
!
!         prect(i) = prect(i)+(qrtend_permute(k,i)*dz_permute(k,i)*rho_permute(k,i)+&
!              qnitend_permute(k,i)*dz_permute(k,i)*rho_permute(k,i))/rhow
!         preci(i) = preci(i)+qnitend_permute(k,i)*dz_permute(k,i)*rho_permute(k,i)/rhow
!
!         ! convert rain rate from m/s to mm/hr
!
!         rainrt_permute(k,i)=qric_permute(k,i)*rho_permute(k,i)*umr(k)/rhow*3600._r8*1000._r8
!
!         ! vertically-integrated precip source/sink terms (note: grid-averaged)
!
!         qrtot = max(qrtot+qrtend_permute(k,i)*dz_permute(k,i)*rho_permute(k,i),0._r8)
!         qstot = max(qstot+qnitend_permute(k,i)*dz_permute(k,i)*rho_permute(k,i),0._r8)
!         nrtot = max(nrtot+nrtend_permute(k,i)*dz_permute(k,i)*rho_permute(k,i),0._r8)
!         nstot = max(nstot+nstend_permute(k,i)*dz_permute(k,i)*rho_permute(k,i),0._r8)
!
!         ! calculate melting and freezing of precip
!
!         ! melt snow at +2 C
!
!         if (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat > 275.15_r8) then
!            if (qstot > 0._r8) then
!
!               ! make sure melting snow doesn't reduce temperature below threshold
!               dum = -xlf/cpp*qstot/(dz_permute(k,i)*rho_permute(k,i))
!               if (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat+dum.lt.275.15_r8) then
!                  dum = (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat-275.15_r8)*cpp/xlf
!                  dum = dum/(xlf/cpp*qstot/(dz_permute(k,i)*rho_permute(k,i)))
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qric_permute(k,i)=qric_permute(k,i)+dum*qniic_permute(k,i)
!               nric_permute(k,i)=nric_permute(k,i)+dum*nsic_permute(k,i)
!               qniic_permute(k,i)=(1._r8-dum)*qniic_permute(k,i)
!               nsic_permute(k,i)=(1._r8-dum)*nsic_permute(k,i)
!               ! heating tendency 
!               tmp=-xlf*dum*qstot/(dz_permute(k,i)*rho_permute(k,i))
!               meltsdt_permute(k,i)=meltsdt_permute(k,i) + tmp
!
!               tlat_permute(k,i)=tlat_permute(k,i)+tmp
!               qrtot=qrtot+dum*qstot
!               nrtot=nrtot+dum*nstot
!               qstot=(1._r8-dum)*qstot
!               nstot=(1._r8-dum)*nstot
!               preci(i)=(1._r8-dum)*preci(i)
!            end if
!         end if
!
!         ! freeze all rain at -5C for Arctic
!
!         if (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat < (tmelt - 5._r8)) then
!
!            if (qrtot > 0._r8) then
!
!               ! make sure freezing rain doesn't increase temperature above threshold
!               dum = xlf/cpp*qrtot/(dz_permute(k,i)*rho_permute(k,i))
!               if (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat+dum.gt.(tmelt - 5._r8)) then
!                  dum = -(t_permute(k,i)+tlat_permute(k,i)/cpp*deltat-(tmelt-5._r8))*cpp/xlf
!                  dum = dum/(xlf/cpp*qrtot/(dz_permute(k,i)*rho_permute(k,i)))
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qniic_permute(k,i)=qniic_permute(k,i)+dum*qric_permute(k,i)
!               nsic_permute(k,i)=nsic_permute(k,i)+dum*nric_permute(k,i)
!               qric_permute(k,i)=(1._r8-dum)*qric_permute(k,i)
!               nric_permute(k,i)=(1._r8-dum)*nric_permute(k,i)
!               ! heating tendency 
!               tmp = xlf*dum*qrtot/(dz_permute(k,i)*rho_permute(k,i))
!               frzrdt_permute(k,i)=frzrdt_permute(k,i) + tmp
!
!               tlat_permute(k,i)=tlat_permute(k,i)+tmp
!               qstot=qstot+dum*qrtot
!               qrtot=(1._r8-dum)*qrtot
!               nstot=nstot+dum*nrtot
!               nrtot=(1._r8-dum)*nrtot
!               preci(i)=preci(i)+dum*(prect(i)-preci(i))
!            end if
!         end if
!
!         ! if rain/snow mix ratio is zero so should number concentration
!
!         if (qniic_permute(k,i).lt.qsmall) then
!            qniic_permute(k,i)=0._r8
!            nsic_permute(k,i)=0._r8
!         end if
!
!         if (qric_permute(k,i).lt.qsmall) then
!            qric_permute(k,i)=0._r8
!            nric_permute(k,i)=0._r8
!         end if
!
!         ! make sure number concentration is a positive number to avoid 
!         ! taking root of negative
!
!         nric_permute(k,i)=max(nric_permute(k,i),0._r8)
!         nsic_permute(k,i)=max(nsic_permute(k,i),0._r8)
!
!         !.......................................................................
!         ! get size distribution parameters for fallspeed calculations
!         !......................................................................
!         ! rain
!
!         if (qric_permute(k,i).ge.qsmall) then
!            lamr(k) = (pi*rhow*nric_permute(k,i)/qric_permute(k,i))**(1._r8/3._r8)
!            n0r(k) = nric_permute(k,i)*lamr(k)
!
!            ! check for slope
!            ! change lammax and lammin for rain and snow
!            ! adjust vars
!
!            if (lamr(k).lt.lamminr) then
!
!               lamr(k) = lamminr
!
!               n0r(k) = lamr(k)**4*qric_permute(k,i)/(pi*rhow)
!               nric_permute(k,i) = n0r(k)/lamr(k)
!            else if (lamr(k).gt.lammaxr) then
!               lamr(k) = lammaxr
!               n0r(k) = lamr(k)**4*qric_permute(k,i)/(pi*rhow)
!               nric_permute(k,i) = n0r(k)/lamr(k)
!            end if
!
!
!            ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
!
!            unr(k) = min(arn_permute(k,i)*cons4/lamr(k)**br,9.1_r8*rhof_permute(k,i))
!            umr(k) = min(arn_permute(k,i)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof_permute(k,i))
!
!         else
!            lamr(k) = 0._r8
!            n0r(k) = 0._r8
!            umr(k)=0._r8
!            unr(k)=0._r8
!         end if
!
!         !calculate mean size of combined rain and snow
!
!         if (lamr(k).gt.0._r8) then
!            Artmp = n0r(k) * pi / (2._r8 * lamr(k)**3._r8)
!         else 
!            Artmp = 0._r8
!         endif
!
!         if (lamc(k).gt.0._r8) then
!            Actmp = cdist1(k) * pi * gamma(pgam(k)+3._r8)/(4._r8 * lamc(k)**2._r8)
!         else 
!            Actmp = 0._r8
!         endif
!
!         if (Actmp.gt.0_r8.or.Artmp.gt.0) then
!            rercld_permute(k,i)=rercld_permute(k,i) + 3._r8 *(qric_permute(k,i) + qcic_permute(k,i)) / (4._r8 * rhow * (Actmp + Artmp))
!            arcld_permute(k,i)=arcld_permute(k,i)+1._r8
!         endif
!
!         !......................................................................
!         ! snow
!
!         if (qniic_permute(k,i).ge.qsmall) then
!            lams(k) = (cons6*cs*nsic_permute(k,i)/ &
!                 qniic_permute(k,i))**(1._r8/ds)
!            n0s(k) = nsic_permute(k,i)*lams(k)
!
!            ! check for slope
!            ! adjust vars
!
!            if (lams(k).lt.lammins) then
!               lams(k) = lammins
!               n0s(k) = lams(k)**(ds+1._r8)*qniic_permute(k,i)/(cs*cons6)
!               nsic_permute(k,i) = n0s(k)/lams(k)
!
!            else if (lams(k).gt.lammaxs) then
!               lams(k) = lammaxs
!               n0s(k) = lams(k)**(ds+1._r8)*qniic_permute(k,i)/(cs*cons6)
!               nsic_permute(k,i) = n0s(k)/lams(k)
!            end if
!
!            ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
!
!            ums(k) = min(asn_permute(k,i)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof_permute(k,i))
!            uns(k) = min(asn_permute(k,i)*cons7/lams(k)**bs,1.2_r8*rhof_permute(k,i))
!
!         else
!            lams(k) = 0._r8
!            n0s(k) = 0._r8
!            ums(k) = 0._r8
!            uns(k) = 0._r8
!         end if
!
!         !c........................................................................
!         ! sum over sub-step for average process rates
!
!         ! convert rain/snow q and N for output to history, note, 
!         ! output is for gridbox average
!
!         qrout_permute(k,i)=qrout_permute(k,i)+qric_permute(k,i)*cldmax_permute(k,i)
!         qsout_permute(k,i)=qsout_permute(k,i)+qniic_permute(k,i)*cldmax_permute(k,i)
!         nrout_permute(k,i)=nrout_permute(k,i)+nric_permute(k,i)*rho_permute(k,i)*cldmax_permute(k,i)
!         nsout_permute(k,i)=nsout_permute(k,i)+nsic_permute(k,i)*rho_permute(k,i)*cldmax_permute(k,i)
!
!         tlat1_permute(k,i)=tlat1_permute(k,i)+tlat_permute(k,i)
!         qvlat1_permute(k,i)=qvlat1_permute(k,i)+qvlat_permute(k,i)
!         qctend1_permute(k,i)=qctend1_permute(k,i)+qctend_permute(k,i)
!         qitend1_permute(k,i)=qitend1_permute(k,i)+qitend_permute(k,i)
!         nctend1_permute(k,i)=nctend1_permute(k,i)+nctend_permute(k,i)
!         nitend1_permute(k,i)=nitend1_permute(k,i)+nitend_permute(k,i)
!
!         t_permute(k,i)=t_permute(k,i)+tlat_permute(k,i)*deltat/cpp
!         q_permute(k,i)=q_permute(k,i)+qvlat_permute(k,i)*deltat
!         qc_permute(k,i)=qc_permute(k,i)+qctend_permute(k,i)*deltat
!         qi_permute(k,i)=qi_permute(k,i)+qitend_permute(k,i)*deltat
!         nc_permute(k,i)=nc_permute(k,i)+nctend_permute(k,i)*deltat
!         ni_permute(k,i)=ni_permute(k,i)+nitend_permute(k,i)*deltat
!
!         rainrt1_permute(k,i)=rainrt1_permute(k,i)+rainrt_permute(k,i)
!
!         !divide rain radius over substeps for average
!         if (arcld_permute(k,i) .gt. 0._r8) then
!            rercld_permute(k,i)=rercld_permute(k,i)/arcld_permute(k,i)
!         end if
!
!         !calculate precip fluxes and adding them to summing sub-stepping variables
!         !! flux is zero at top interface
!         rflx_permute(1,i)=0.0_r8
!         sflx_permute(1,i)=0.0_r8
!
!         !! calculating the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
!         rflx_permute(k+1,i)=qrout_permute(k,i)*rho_permute(k,i)*umr(k)
!         sflx_permute(k+1,i)=qsout_permute(k,i)*rho_permute(k,i)*ums(k)
!
!         !! add to summing sub-stepping variable
!         rflx1_permute(k+1,i)=rflx1_permute(k+1,i)+rflx_permute(k+1,i)
!         sflx1_permute(k+1,i)=sflx1_permute(k+1,i)+sflx_permute(k+1,i)
!
!         !c........................................................................
!
!      end do ! k loop
!
!      prect1(i)=prect1(i)+prect(i)
!      preci1(i)=preci1(i)+preci(i)
!
!   end do ! it loop, sub-step
!
!   do k = top_lev, pver
!      rate1ord_cw2pr_st_permute(k,i) = qcsinksum_rate1ord(k)/max(qcsum_rate1ord(k),1.0e-30_r8) 
!   end do
!
!300 continue  ! continue if no cloud water
!end do ! i loop
!
!! convert dt from sub-step back to full time step
!deltat=deltat*real(iter)
!
!!c.............................................................................
!
!do i=1,ncol
!
!   ! skip all calculations if no cloud water
!   if (ltrue(i).eq.0) then
!
!      do k=1,top_lev-1
!         ! assign zero values for effective radius above 1 mbar
!         effc_permute(k,i)=0._r8
!         effi_permute(k,i)=0._r8
!         effc_fn_permute(k,i)=0._r8
!         lamcrad_permute(k,i)=0._r8
!         pgamrad_permute(k,i)=0._r8
!         deffi_permute(k,i)=0._r8
!      end do
!
!      do k=top_lev,pver
!         ! assign default values for effective radius
!         effc_permute(k,i)=10._r8
!         effi_permute(k,i)=25._r8
!         effc_fn_permute(k,i)=10._r8
!         lamcrad_permute(k,i)=0._r8
!         pgamrad_permute(k,i)=0._r8
!         deffi_permute(k,i)=0._r8
!      end do
!      goto 500
!   end if
!
!   ! initialize nstep for sedimentation sub-steps
!   nstep = 1
!
!   ! divide precip rate by number of sub-steps to get average over time step
!
!   prect(i)=prect1(i)/real(iter)
!   preci(i)=preci1(i)/real(iter)
!
!   do k=top_lev,pver
!
!      ! assign variables back to start-of-timestep values before updating after sub-steps 
!
!      t_permute(k,i)=t1_permute(k,i)
!      q_permute(k,i)=q1_permute(k,i)
!      qc_permute(k,i)=qc1_permute(k,i)
!      qi_permute(k,i)=qi1_permute(k,i)
!      nc_permute(k,i)=nc1_permute(k,i)
!      ni_permute(k,i)=ni1_permute(k,i)
!
!      ! divide microphysical tendencies by number of sub-steps to get average over time step
!
!      tlat_permute(k,i)=tlat1_permute(k,i)/real(iter)
!      qvlat_permute(k,i)=qvlat1_permute(k,i)/real(iter)
!      qctend_permute(k,i)=qctend1_permute(k,i)/real(iter)
!      qitend_permute(k,i)=qitend1_permute(k,i)/real(iter)
!      nctend_permute(k,i)=nctend1_permute(k,i)/real(iter)
!      nitend_permute(k,i)=nitend1_permute(k,i)/real(iter)
!
!      rainrt_permute(k,i)=rainrt1_permute(k,i)/real(iter)
!
!      ! divide by number of sub-steps to find final values
!      rflx_permute(k+1,i)=rflx1_permute(k+1,i)/real(iter)
!      sflx_permute(k+1,i)=sflx1_permute(k+1,i)/real(iter)
!
!      ! divide output precip q and N by number of sub-steps to get average over time step
!
!      qrout_permute(k,i)=qrout_permute(k,i)/real(iter)
!      qsout_permute(k,i)=qsout_permute(k,i)/real(iter)
!      nrout_permute(k,i)=nrout_permute(k,i)/real(iter)
!      nsout_permute(k,i)=nsout_permute(k,i)/real(iter)
!
!      ! divide trop_mozart variables by number of sub-steps to get average over time step 
!
!      nevapr_permute(k,i) = nevapr_permute(k,i)/real(iter)
!      nevapr2_permute(k,i) = nevapr2_permute(k,i)/real(iter)
!      evapsnow_permute(k,i) = evapsnow_permute(k,i)/real(iter)
!      prain_permute(k,i) = prain_permute(k,i)/real(iter)
!      prodsnow_permute(k,i) = prodsnow_permute(k,i)/real(iter)
!      cmeout_permute(k,i) = cmeout_permute(k,i)/real(iter)
!
!      cmeiout_permute(k,i) = cmeiout_permute(k,i)/real(iter)
!      meltsdt_permute(k,i) = meltsdt_permute(k,i)/real(iter)
!      frzrdt_permute(k,i) = frzrdt_permute(k,i)/real(iter)
!
!
!      ! microphysics output
!      prao_permute(k,i)=prao_permute(k,i)/real(iter)
!      prco_permute(k,i)=prco_permute(k,i)/real(iter)
!      mnuccco_permute(k,i)=mnuccco_permute(k,i)/real(iter)
!      mnuccto_permute(k,i)=mnuccto_permute(k,i)/real(iter)
!      msacwio_permute(k,i)=msacwio_permute(k,i)/real(iter)
!      psacwso_permute(k,i)=psacwso_permute(k,i)/real(iter)
!      bergso_permute(k,i)=bergso_permute(k,i)/real(iter)
!      bergo_permute(k,i)=bergo_permute(k,i)/real(iter)
!      prcio_permute(k,i)=prcio_permute(k,i)/real(iter)
!      praio_permute(k,i)=praio_permute(k,i)/real(iter)
!
!      mnuccro_permute(k,i)=mnuccro_permute(k,i)/real(iter)
!      pracso_permute(k,i)=pracso_permute(k,i)/real(iter)
!
!      mnuccdo_permute(k,i)=mnuccdo_permute(k,i)/real(iter)
!
!      ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
!      nevapr_permute(k,i) = nevapr_permute(k,i) + evapsnow_permute(k,i)
!      prer_evap_permute(k,i) = nevapr2_permute(k,i)
!      prain_permute(k,i) = prain_permute(k,i) + prodsnow_permute(k,i)
!
!      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      ! calculate sedimentation for cloud water and ice
!      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      ! update in-cloud cloud mixing ratio and number concentration 
!      ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
!      ! note: these are in-cloud values***, hence we divide by cloud fraction
!
!      dumc_permute(k,i) = (qc_permute(k,i)+qctend_permute(k,i)*deltat)/lcldm_permute(k,i)
!      dumi_permute(k,i) = (qi_permute(k,i)+qitend_permute(k,i)*deltat)/icldm_permute(k,i)
!      dumnc_permute(k,i) = max((nc_permute(k,i)+nctend_permute(k,i)*deltat)/lcldm_permute(k,i),0._r8)
!      dumni_permute(k,i) = max((ni_permute(k,i)+nitend_permute(k,i)*deltat)/icldm_permute(k,i),0._r8)
!
!      ! obtain new slope parameter to avoid possible singularity
!
!      if (dumi_permute(k,i).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumni_permute(k,i)=min(dumni_permute(k,i),dumi_permute(k,i)*1.e20_r8)
!
!         lami(k) = (cons1*ci* &
!              dumni_permute(k,i)/dumi_permute(k,i))**(1._r8/di)
!         lami(k)=max(lami(k),lammini)
!         lami(k)=min(lami(k),lammaxi)
!      else
!         lami(k)=0._r8
!      end if
!
!      if (dumc_permute(k,i).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumnc_permute(k,i)=min(dumnc_permute(k,i),dumc_permute(k,i)*1.e20_r8)
!         ! add lower limit to in-cloud number concentration
!         dumnc_permute(k,i)=max(dumnc_permute(k,i),cdnl/rho_permute(k,i)) ! sghan minimum in #/cm3 
!         pgam(k)=0.0005714_r8*(ncic_permute(k,i)/1.e6_r8*rho_permute(k,i))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)**2)-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc_permute(k,i)*gamma(pgam(k)+4._r8)/ &
!              (dumc_permute(k,i)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         lammax = (pgam(k)+1._r8)/2.e-6_r8
!         lamc(k)=max(lamc(k),lammin)
!         lamc(k)=min(lamc(k),lammax)
!      else
!         lamc(k)=0._r8
!      end if
!
!      ! calculate number and mass weighted fall velocity for droplets
!      ! include effects of sub-grid distribution of cloud water
!
!
!      if (dumc_permute(k,i).ge.qsmall) then
!         unc = acn_permute(k,i)*gamma(1._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+1._r8))
!         umc = acn_permute(k,i)*gamma(4._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+4._r8))
!         ! fallspeed for output
!         vtrmc_permute(k,i)=umc
!      else
!         umc = 0._r8
!         unc = 0._r8
!      end if
!
!      ! calculate number and mass weighted fall velocity for cloud ice
!
!      if (dumi_permute(k,i).ge.qsmall) then
!         uni =  ain_permute(k,i)*cons16/lami(k)**bi
!         umi = ain_permute(k,i)*cons17/(6._r8*lami(k)**bi)
!         uni=min(uni,1.2_r8*rhof_permute(k,i))
!         umi=min(umi,1.2_r8*rhof_permute(k,i))
!
!         ! fallspeed
!         vtrmi_permute(k,i)=umi
!      else
!         umi = 0._r8
!         uni = 0._r8
!      end if
!
!      fi(k) = g*rho_permute(k,i)*umi
!      fni(k) = g*rho_permute(k,i)*uni
!      fc(k) = g*rho_permute(k,i)*umc
!      fnc(k) = g*rho_permute(k,i)*unc
!
!      ! calculate number of split time steps to ensure courant stability criteria
!      ! for sedimentation calculations
!
!      rgvm = max(fi(k),fc(k),fni(k),fnc(k))
!      nstep = max(int(rgvm*deltat/pdel_permute(k,i)+1._r8),nstep)
!
!      ! redefine dummy variables - sedimentation is calculated over grid-scale
!      ! quantities to ensure conservation
!
!      dumc_permute(k,i) = (qc_permute(k,i)+qctend_permute(k,i)*deltat)
!      dumi_permute(k,i) = (qi_permute(k,i)+qitend_permute(k,i)*deltat)
!      dumnc_permute(k,i) = max((nc_permute(k,i)+nctend_permute(k,i)*deltat),0._r8)
!      dumni_permute(k,i) = max((ni_permute(k,i)+nitend_permute(k,i)*deltat),0._r8)
!
!      if (dumc_permute(k,i).lt.qsmall) dumnc_permute(k,i)=0._r8
!      if (dumi_permute(k,i).lt.qsmall) dumni_permute(k,i)=0._r8
!
!   end do       !!! vertical loop
!   do n = 1,nstep  !! loop over sub-time step to ensure stability
!
!      do k = top_lev,pver
!         if (do_cldice) then
!            falouti(k) = fi(k)*dumi_permute(k,i)
!            faloutni(k) = fni(k)*dumni_permute(k,i)
!         else
!            falouti(k)  = 0._r8
!            faloutni(k) = 0._r8
!         end if
!
!         faloutc(k) = fc(k)*dumc_permute(k,i)
!         faloutnc(k) = fnc(k)*dumnc_permute(k,i)
!      end do
!
!      ! top of model
!
!      k = top_lev
!      faltndi = falouti(k)/pdel_permute(k,i)
!      faltndni = faloutni(k)/pdel_permute(k,i)
!      faltndc = faloutc(k)/pdel_permute(k,i)
!      faltndnc = faloutnc(k)/pdel_permute(k,i)
!
!      ! add fallout terms to microphysical tendencies
!
!      qitend_permute(k,i) = qitend_permute(k,i)-faltndi/nstep
!      nitend_permute(k,i) = nitend_permute(k,i)-faltndni/nstep
!      qctend_permute(k,i) = qctend_permute(k,i)-faltndc/nstep
!      nctend_permute(k,i) = nctend_permute(k,i)-faltndnc/nstep
!
!      ! sedimentation tendencies for output
!      qcsedten_permute(k,i)=qcsedten_permute(k,i)-faltndc/nstep
!      qisedten_permute(k,i)=qisedten_permute(k,i)-faltndi/nstep
!
!      dumi_permute(k,i) = dumi_permute(k,i)-faltndi*deltat/nstep
!      dumni_permute(k,i) = dumni_permute(k,i)-faltndni*deltat/nstep
!      dumc_permute(k,i) = dumc_permute(k,i)-faltndc*deltat/nstep
!      dumnc_permute(k,i) = dumnc_permute(k,i)-faltndnc*deltat/nstep
!
!      do k = top_lev+1,pver
!
!         ! for cloud liquid and ice, if cloud fraction increases with height
!         ! then add flux from above to both vapor and cloud water of current level
!         ! this means that flux entering clear portion of cell from above evaporates
!         ! instantly
!
!         dum=lcldm_permute(k,i)/lcldm_permute(k-1,i)
!         dum=min(dum,1._r8)
!         dum1=icldm_permute(k,i)/icldm_permute(k-1,i)
!         dum1=min(dum1,1._r8)
!
!         faltndqie=(falouti(k)-falouti(k-1))/pdel_permute(k,i)
!         faltndi=(falouti(k)-dum1*falouti(k-1))/pdel_permute(k,i)
!         faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel_permute(k,i)
!         faltndqce=(faloutc(k)-faloutc(k-1))/pdel_permute(k,i)
!         faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel_permute(k,i)
!         faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel_permute(k,i)
!
!         ! add fallout terms to eulerian tendencies
!
!         qitend_permute(k,i) = qitend_permute(k,i)-faltndi/nstep
!         nitend_permute(k,i) = nitend_permute(k,i)-faltndni/nstep
!         qctend_permute(k,i) = qctend_permute(k,i)-faltndc/nstep
!         nctend_permute(k,i) = nctend_permute(k,i)-faltndnc/nstep
!
!         ! sedimentation tendencies for output
!         qcsedten_permute(k,i)=qcsedten_permute(k,i)-faltndc/nstep
!         qisedten_permute(k,i)=qisedten_permute(k,i)-faltndi/nstep
!
!         ! add terms to to evap/sub of cloud water
!
!         qvlat_permute(k,i)=qvlat_permute(k,i)-(faltndqie-faltndi)/nstep
!         ! for output
!         qisevap_permute(k,i)=qisevap_permute(k,i)-(faltndqie-faltndi)/nstep
!         qvlat_permute(k,i)=qvlat_permute(k,i)-(faltndqce-faltndc)/nstep
!         ! for output
!         qcsevap_permute(k,i)=qcsevap_permute(k,i)-(faltndqce-faltndc)/nstep
!
!         tlat_permute(k,i)=tlat_permute(k,i)+(faltndqie-faltndi)*xxls/nstep
!         tlat_permute(k,i)=tlat_permute(k,i)+(faltndqce-faltndc)*xxlv/nstep
!
!         dumi_permute(k,i) = dumi_permute(k,i)-faltndi*deltat/nstep
!         dumni_permute(k,i) = dumni_permute(k,i)-faltndni*deltat/nstep
!         dumc_permute(k,i) = dumc_permute(k,i)-faltndc*deltat/nstep
!         dumnc_permute(k,i) = dumnc_permute(k,i)-faltndnc*deltat/nstep
!
!         Fni(K)=MAX(Fni(K)/pdel_permute(K,i),Fni(K-1)/pdel_permute(K-1,i))*pdel_permute(K,i)
!         FI(K)=MAX(FI(K)/pdel_permute(K,i),FI(K-1)/pdel_permute(K-1,i))*pdel_permute(K,i)
!         fnc(k)=max(fnc(k)/pdel_permute(k,i),fnc(k-1)/pdel_permute(k-1,i))*pdel_permute(k,i)
!         Fc(K)=MAX(Fc(K)/pdel_permute(K,i),Fc(K-1)/pdel_permute(K-1,i))*pdel_permute(K,i)
!
!      end do   !! k loop
!
!      ! units below are m/s
!      ! cloud water/ice sedimentation flux at surface 
!      ! is added to precip flux at surface to get total precip (cloud + precip water)
!      ! rate
!
!      prect(i) = prect(i)+(faloutc(pver)+falouti(pver))/g/nstep/1000._r8  
!      preci(i) = preci(i)+(falouti(pver))/g/nstep/1000._r8
!
!   end do   !! nstep loop
!
!   ! end sedimentation
!   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   ! get new update for variables that includes sedimentation tendency
!   ! note : here dum variables are grid-average, NOT in-cloud
!
!   do k=top_lev,pver
!
!      dumc_permute(k,i) = max(qc_permute(k,i)+qctend_permute(k,i)*deltat,0._r8)
!      dumi_permute(k,i) = max(qi_permute(k,i)+qitend_permute(k,i)*deltat,0._r8)
!      dumnc_permute(k,i) = max(nc_permute(k,i)+nctend_permute(k,i)*deltat,0._r8)
!      dumni_permute(k,i) = max(ni_permute(k,i)+nitend_permute(k,i)*deltat,0._r8)
!
!      if (dumc_permute(k,i).lt.qsmall) dumnc_permute(k,i)=0._r8
!      if (dumi_permute(k,i).lt.qsmall) dumni_permute(k,i)=0._r8
!
!      ! calculate instantaneous processes (melting, homogeneous freezing)
!      if (do_cldice) then
!
!         if (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat > tmelt) then
!            if (dumi_permute(k,i) > 0._r8) then
!
!               ! limit so that melting does not push temperature below freezing
!               dum = -dumi_permute(k,i)*xlf/cpp
!               if (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat+dum.lt.tmelt) then
!                  dum = (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat-tmelt)*cpp/xlf
!                  dum = dum/dumi_permute(k,i)*xlf/cpp 
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qctend_permute(k,i)=qctend_permute(k,i)+dum*dumi_permute(k,i)/deltat
!
!               ! for output
!               melto_permute(k,i)=dum*dumi_permute(k,i)/deltat
!
!               ! assume melting ice produces droplet
!               ! mean volume radius of 8 micron
!
!               nctend_permute(k,i)=nctend_permute(k,i)+3._r8*dum*dumi_permute(k,i)/deltat/ &
!                    (4._r8*pi*5.12e-16_r8*rhow)
!
!               qitend_permute(k,i)=((1._r8-dum)*dumi_permute(k,i)-qi_permute(k,i))/deltat
!               nitend_permute(k,i)=((1._r8-dum)*dumni_permute(k,i)-ni_permute(k,i))/deltat
!               tlat_permute(k,i)=tlat_permute(k,i)-xlf*dum*dumi_permute(k,i)/deltat
!            end if
!         end if
!
!         ! homogeneously freeze droplets at -40 C
!
!         if (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat < 233.15_r8) then
!            if (dumc_permute(k,i) > 0._r8) then
!
!               ! limit so that freezing does not push temperature above threshold
!               dum = dumc_permute(k,i)*xlf/cpp
!               if (t_permute(k,i)+tlat_permute(k,i)/cpp*deltat+dum.gt.233.15_r8) then
!                  dum = -(t_permute(k,i)+tlat_permute(k,i)/cpp*deltat-233.15_r8)*cpp/xlf
!                  dum = dum/dumc_permute(k,i)*xlf/cpp
!                  dum = max(0._r8,dum)
!                  dum = min(1._r8,dum)
!               else
!                  dum = 1._r8
!               end if
!
!               qitend_permute(k,i)=qitend_permute(k,i)+dum*dumc_permute(k,i)/deltat
!               ! for output
!               homoo_permute(k,i)=dum*dumc_permute(k,i)/deltat
!
!               ! assume 25 micron mean volume radius of homogeneously frozen droplets
!               ! consistent with size of detrained ice in stratiform.F90
!               nitend_permute(k,i)=nitend_permute(k,i)+dum*3._r8*dumc_permute(k,i)/(4._r8*3.14_r8*1.563e-14_r8* &
!                    500._r8)/deltat
!               qctend_permute(k,i)=((1._r8-dum)*dumc_permute(k,i)-qc_permute(k,i))/deltat
!               nctend_permute(k,i)=((1._r8-dum)*dumnc_permute(k,i)-nc_permute(k,i))/deltat
!               tlat_permute(k,i)=tlat_permute(k,i)+xlf*dum*dumc_permute(k,i)/deltat
!            end if
!         end if
!
!         ! remove any excess over-saturation, which is possible due to non-linearity when adding 
!         ! together all microphysical processes
!         ! follow code similar to old CAM scheme
!
!         qtmp=q_permute(k,i)+qvlat_permute(k,i)*deltat
!         ttmp=t_permute(k,i)+tlat_permute(k,i)/cpp*deltat
!
!         esn = svp_water(ttmp)  ! use rhw to allow ice supersaturation
!         qsn = svp_to_qsat(esn, p_permute(k,i))
!
!         if (qtmp > qsn .and. qsn > 0) then
!            ! expression below is approximate since there may be ice deposition
!            dum = (qtmp-qsn)/(1._r8+cons27*qsn/(cpp*rv*ttmp**2))/deltat
!            ! add to output cme
!            cmeout_permute(k,i) = cmeout_permute(k,i)+dum
!            ! now add to tendencies, partition between liquid and ice based on temperature
!            if (ttmp > 268.15_r8) then
!               dum1=0.0_r8
!               ! now add to tendencies, partition between liquid and ice based on te
!            else if (ttmp < 238.15_r8) then
!               dum1=1.0_r8
!            else
!               dum1=(268.15_r8-ttmp)/30._r8
!            end if
!
!            dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
!                 *qsn/(cpp*rv*ttmp**2))/deltat
!            qctend_permute(k,i)=qctend_permute(k,i)+dum*(1._r8-dum1)
!            ! for output
!            qcreso_permute(k,i)=dum*(1._r8-dum1)
!            qitend_permute(k,i)=qitend_permute(k,i)+dum*dum1
!            qireso_permute(k,i)=dum*dum1
!            qvlat_permute(k,i)=qvlat_permute(k,i)-dum
!            ! for output
!            qvres_permute(k,i)=-dum
!            tlat_permute(k,i)=tlat_permute(k,i)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
!         end if
!      end if
!
!      !...............................................................................
!      ! calculate effective radius for pass to radiation code
!      ! if no cloud water, default value is 10 micron for droplets,
!      ! 25 micron for cloud ice
!
!      ! update cloud variables after instantaneous processes to get effective radius
!      ! variables are in-cloud to calculate size dist parameters
!
!      dumc_permute(k,i) = max(qc_permute(k,i)+qctend_permute(k,i)*deltat,0._r8)/lcldm_permute(k,i)
!      dumi_permute(k,i) = max(qi_permute(k,i)+qitend_permute(k,i)*deltat,0._r8)/icldm_permute(k,i)
!      dumnc_permute(k,i) = max(nc_permute(k,i)+nctend_permute(k,i)*deltat,0._r8)/lcldm_permute(k,i)
!      dumni_permute(k,i) = max(ni_permute(k,i)+nitend_permute(k,i)*deltat,0._r8)/icldm_permute(k,i)
!
!      ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
!
!      dumc_permute(k,i)=min(dumc_permute(k,i),5.e-3_r8)
!      dumi_permute(k,i)=min(dumi_permute(k,i),5.e-3_r8)
!
!      !...................
!      ! cloud ice effective radius
!
!      if (dumi_permute(k,i).ge.qsmall) then
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumni_permute(k,i)=min(dumni_permute(k,i),dumi_permute(k,i)*1.e20_r8)
!         lami(k) = (cons1*ci*dumni_permute(k,i)/dumi_permute(k,i))**(1._r8/di)
!
!         if (lami(k).lt.lammini) then
!            lami(k) = lammini
!            n0i(k) = lami(k)**(di+1._r8)*dumi_permute(k,i)/(ci*cons1)
!            niic_permute(k,i) = n0i(k)/lami(k)
!            ! adjust number conc if needed to keep mean size in reasonable range
!            if (do_cldice) nitend_permute(k,i)=(niic_permute(k,i)*icldm_permute(k,i)-ni_permute(k,i))/deltat
!
!         else if (lami(k).gt.lammaxi) then
!            lami(k) = lammaxi
!            n0i(k) = lami(k)**(di+1._r8)*dumi_permute(k,i)/(ci*cons1)
!            niic_permute(k,i) = n0i(k)/lami(k)
!            ! adjust number conc if needed to keep mean size in reasonable range
!            if (do_cldice) nitend_permute(k,i)=(niic_permute(k,i)*icldm_permute(k,i)-ni_permute(k,i))/deltat
!         end if
!         effi_permute(k,i) = 1.5_r8/lami(k)*1.e6_r8
!
!      else
!         effi_permute(k,i) = 25._r8
!      end if
!
!      ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
!      ! radius has already been determined from the size distribution.
!      if (.not. do_cldice) then
!         effi_permute(k,i) = re_ice_permute(k,i) * 1e6_r8      ! m -> um
!      end if
!
!      !...................
!      ! cloud droplet effective radius
!
!      if (dumc_permute(k,i).ge.qsmall) then
!
!         ! add upper limit to in-cloud number concentration to prevent numerical error
!         dumnc_permute(k,i)=min(dumnc_permute(k,i),dumc_permute(k,i)*1.e20_r8)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! set tendency to ensure minimum droplet concentration
!         ! after update by microphysics, except when lambda exceeds bounds on mean drop
!         ! size or if there is no cloud water
!         if (dumnc_permute(k,i).lt.cdnl/rho_permute(k,i)) then   
!            nctend_permute(k,i)=(cdnl/rho_permute(k,i)*lcldm_permute(k,i)-nc_permute(k,i))/deltat   
!         end if
!         dumnc_permute(k,i)=max(dumnc_permute(k,i),cdnl/rho_permute(k,i)) ! sghan minimum in #/cm3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         pgam(k)=0.0005714_r8*(ncic_permute(k,i)/1.e6_r8*rho_permute(k,i))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)**2)-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc_permute(k,i)*gamma(pgam(k)+4._r8)/ &
!              (dumc_permute(k,i)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         ! Multiply by omsm to fit within RRTMG's table.
!         lammax = (pgam(k)+1._r8)*omsm/2.e-6_r8
!         if (lamc(k).lt.lammin) then
!            lamc(k) = lammin
!            ncic_permute(k,i) = 6._r8*lamc(k)**3*dumc_permute(k,i)* &
!                 gamma(pgam(k)+1._r8)/ &
!                 (pi*rhow*gamma(pgam(k)+4._r8))
!            ! adjust number conc if needed to keep mean size in reasonable range
!            nctend_permute(k,i)=(ncic_permute(k,i)*lcldm_permute(k,i)-nc_permute(k,i))/deltat
!
!         else if (lamc(k).gt.lammax) then
!            lamc(k) = lammax
!            ncic_permute(k,i) = 6._r8*lamc(k)**3*dumc_permute(k,i)* &
!                 gamma(pgam(k)+1._r8)/ &
!                 (pi*rhow*gamma(pgam(k)+4._r8))
!            ! adjust number conc if needed to keep mean size in reasonable range
!            nctend_permute(k,i)=(ncic_permute(k,i)*lcldm_permute(k,i)-nc_permute(k,i))/deltat
!         end if
!
!         effc_permute(k,i) = &
!              gamma(pgam(k)+4._r8)/ &
!              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
!         !assign output fields for shape here
!         lamcrad_permute(k,i)=lamc(k)
!         pgamrad_permute(k,i)=pgam(k)
!
!      else
!         effc_permute(k,i) = 10._r8
!         lamcrad_permute(k,i)=0._r8
!         pgamrad_permute(k,i)=0._r8
!      end if
!
!      ! ice effective diameter for david mitchell's optics
!      if (do_cldice) then
!         deffi_permute(k,i)=effi_permute(k,i)*rhoi/917._r8*2._r8
!      else
!         deffi_permute(k,i)=effi_permute(k,i) * 2._r8
!      end if
!
!
!!!! recalculate effective radius for constant number, in order to separate
!      ! first and second indirect effects
!      ! assume constant number of 10^8 kg-1
!
!      dumnc_permute(k,i)=1.e8_r8
!
!      if (dumc_permute(k,i).ge.qsmall) then
!         pgam(k)=0.0005714_r8*(ncic_permute(k,i)/1.e6_r8*rho_permute(k,i))+0.2714_r8
!         pgam(k)=1._r8/(pgam(k)**2)-1._r8
!         pgam(k)=max(pgam(k),2._r8)
!         pgam(k)=min(pgam(k),15._r8)
!
!         lamc(k) = (pi/6._r8*rhow*dumnc_permute(k,i)*gamma(pgam(k)+4._r8)/ &
!              (dumc_permute(k,i)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam(k)+1._r8)/50.e-6_r8
!         lammax = (pgam(k)+1._r8)/2.e-6_r8
!         if (lamc(k).lt.lammin) then
!            lamc(k) = lammin
!         else if (lamc(k).gt.lammax) then
!            lamc(k) = lammax
!         end if
!         effc_fn_permute(k,i) = &
!              gamma(pgam(k)+4._r8)/ &
!              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
!
!      else
!         effc_fn_permute(k,i) = 10._r8
!      end if
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!
!
!   end do ! vertical k loop
!
!500 continue
!
!   do k=top_lev,pver
!      ! if updated q (after microphysics) is zero, then ensure updated n is also zero
!
!      if (qc_permute(k,i)+qctend_permute(k,i)*deltat.lt.qsmall) nctend_permute(k,i)=-nc_permute(k,i)/deltat
!      if (do_cldice .and. qi_permute(k,i)+qitend_permute(k,i)*deltat.lt.qsmall) nitend_permute(k,i)=-ni_permute(k,i)/deltat
!   end do
!
!end do ! i loop
!
!! add snow ouptut
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qsout_permute(k,i).gt.1.e-7_r8.and.nsout_permute(k,i).gt.0._r8) then
!         dsout_permute(k,i)=3._r8*rhosn/917._r8*(pi * rhosn * nsout_permute(k,i)/qsout_permute(k,i))**(-1._r8/3._r8)
!      endif
!   end do
!end do
!
!!calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual
!do i = 1,ncol
!   do k=top_lev,pver
!      !! RAIN
!      if (qrout_permute(k,i).gt.1.e-7_r8.and.nrout_permute(k,i).gt.0._r8) then
!         reff_rain_permute(k,i)=1.5_r8*(pi * rhow * nrout_permute(k,i)/qrout_permute(k,i))**(-1._r8/3._r8)*1.e6_r8
!      endif
!      !! SNOW
!      if (qsout_permute(k,i).gt.1.e-7_r8.and.nsout_permute(k,i).gt.0._r8) then
!         reff_snow_permute(k,i)=1.5_r8*(pi * rhosn * nsout_permute(k,i)/qsout_permute(k,i))**(-1._r8/3._r8)*1.e6_r8
!      end if
!   end do
!end do
!
!! analytic radar reflectivity
!! formulas from Matthew Shupe, NOAA/CERES
!! *****note: radar reflectivity is local (in-precip average)
!! units of mm^6/m^3
!
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qc_permute(k,i)+qctend_permute(k,i)*deltat.ge.qsmall) then
!         dum=((qc_permute(k,i)+qctend_permute(k,i)*deltat)/lcldm_permute(k,i)*rho_permute(k,i)*1000._r8)**2 &
!              /(0.109_r8*(nc_permute(k,i)+nctend_permute(k,i)*deltat)/lcldm_permute(k,i)*rho_permute(k,i)/1.e6_r8)*lcldm_permute(k,i)/cldmax_permute(k,i)
!      else
!         dum=0._r8
!      end if
!      if (qi_permute(k,i)+qitend_permute(k,i)*deltat.ge.qsmall) then
!         dum1=((qi_permute(k,i)+qitend_permute(k,i)*deltat)*rho_permute(k,i)/icldm_permute(k,i)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm_permute(k,i)/cldmax_permute(k,i)
!      else 
!         dum1=0._r8
!      end if
!
!      if (qsout_permute(k,i).ge.qsmall) then
!         dum1=dum1+(qsout_permute(k,i)*rho_permute(k,i)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
!      end if
!
!      refl_permute(k,i)=dum+dum1
!
!      ! add rain rate, but for 37 GHz formulation instead of 94 GHz
!      ! formula approximated from data of Matrasov (2007)
!      ! rainrt is the rain rate in mm/hr
!      ! reflectivity (dum) is in DBz
!
!      if (rainrt_permute(k,i).ge.0.001_r8) then
!         dum=log10(rainrt_permute(k,i)**6._r8)+16._r8
!
!         ! convert from DBz to mm^6/m^3
!
!         dum = 10._r8**(dum/10._r8)
!      else
!         ! don't include rain rate in R calculation for values less than 0.001 mm/hr
!         dum=0._r8
!      end if
!
!      ! add to refl
!
!      refl_permute(k,i)=refl_permute(k,i)+dum
!
!      !output reflectivity in Z.
!      areflz_permute(k,i)=refl_permute(k,i)
!
!      ! convert back to DBz 
!
!      if (refl_permute(k,i).gt.minrefl) then 
!         refl_permute(k,i)=10._r8*log10(refl_permute(k,i))
!      else
!         refl_permute(k,i)=-9999._r8
!      end if
!
!      !set averaging flag
!      if (refl_permute(k,i).gt.mindbz) then 
!         arefl_permute(k,i)=refl_permute(k,i)
!         frefl_permute(k,i)=1.0_r8  
!      else
!         arefl_permute(k,i)=0._r8
!         areflz_permute(k,i)=0._r8
!         frefl_permute(k,i)=0._r8
!      end if
!
!      ! bound cloudsat reflectivity
!
!      csrfl_permute(k,i)=min(csmax,refl_permute(k,i))
!
!      !set averaging flag
!      if (csrfl_permute(k,i).gt.csmin) then 
!         acsrfl_permute(k,i)=refl_permute(k,i)
!         fcsrfl_permute(k,i)=1.0_r8  
!      else
!         acsrfl_permute(k,i)=0._r8
!         fcsrfl_permute(k,i)=0._r8
!      end if
!
!   end do
!end do
!
!
!! averaging for snow and rain number and diameter
!
!qrout2_permute(:,:)=0._r8
!qsout2_permute(:,:)=0._r8
!nrout2_permute(:,:)=0._r8
!nsout2_permute(:,:)=0._r8
!drout2_permute(:,:)=0._r8
!dsout2_permute(:,:)=0._r8
!freqs_permute(:,:)=0._r8
!freqr_permute(:,:)=0._r8
!do i = 1,ncol
!   do k=top_lev,pver
!      if (qrout_permute(k,i).gt.1.e-7_r8.and.nrout_permute(k,i).gt.0._r8) then
!         qrout2_permute(k,i)=qrout_permute(k,i)
!         nrout2_permute(k,i)=nrout_permute(k,i)
!         drout2_permute(k,i)=(pi * rhow * nrout_permute(k,i)/qrout_permute(k,i))**(-1._r8/3._r8)
!         freqr_permute(k,i)=1._r8
!      endif
!      if (qsout_permute(k,i).gt.1.e-7_r8.and.nsout_permute(k,i).gt.0._r8) then
!         qsout2_permute(k,i)=qsout_permute(k,i)
!         nsout2_permute(k,i)=nsout_permute(k,i)
!         dsout2_permute(k,i)=(pi * rhosn * nsout_permute(k,i)/qsout_permute(k,i))**(-1._r8/3._r8)
!         freqs_permute(k,i)=1._r8
!      endif
!   end do
!end do
!
!! output activated liquid and ice (convert from #/kg -> #/m3)
!do i = 1,ncol
!   do k=top_lev,pver
!      ncai_permute(k,i)=dum2i_permute(k,i)*rho_permute(k,i)
!      ncal_permute(k,i)=dum2l_permute(k,i)*rho_permute(k,i)
!   end do
!end do
!
!
!!redefine fice here....
!nfice_permute(:,:)=0._r8
!do k=top_lev,pver
!   do i=1,ncol
!      dumc_permute(k,i) = (qc_permute(k,i)+qctend_permute(k,i)*deltat)
!      dumi_permute(k,i) = (qi_permute(k,i)+qitend_permute(k,i)*deltat)
!      dumfice=qsout_permute(k,i) + qrout_permute(k,i) + dumc_permute(k,i) + dumi_permute(k,i)  
!
!      if (dumfice.gt.qsmall.and.(qsout_permute(k,i)+dumi_permute(k,i).gt.qsmall)) then
!         nfice_permute(k,i)=(qsout_permute(k,i) + dumi_permute(k,i))/dumfice
!      endif
!
!      if (nfice_permute(k,i).gt.1._r8) then
!         nfice_permute(k,i)=1._r8
!      endif
!
!   enddo
!enddo


!do i=1,pcols
!    do k=1,pver
!
!         !relvar(i,k)                  =    relvar_permute(k,i)            
!         !accre_enhan(i,k)             =    accre_enhan_permute(k,i)       
!         !                                                                 
!         qc(i,k)                      =    qc_permute(k,i)                
!         qi(i,k)                      =    qi_permute(k,i)                
!         nc(i,k)                      =    nc_permute(k,i)                
!         ni(i,k)                      =    ni_permute(k,i)                
!         !p(i,k)                       =    p_permute(k,i)                 
!         !pdel(i,k)                    =    pdel_permute(k,i)              
!         !cldn(i,k)                    =    cldn_permute(k,i)              
!         !icecldf(i,k)                 =    icecldf_permute(k,i)           
!         !liqcldf(i,k)                 =    liqcldf_permute(k,i)           
!         rate1ord_cw2pr_st(i,k)       =    rate1ord_cw2pr_st_permute(k,i) 
!         !                                                                 
!         !naai(i,k)                    =    naai_permute(k,i)              
!         !npccnin(i,k)                 =    npccnin_permute(k,i)           
!         !                                                                 
!         tlat(i,k)                    =    tlat_permute(k,i)              
!         qvlat(i,k)                   =    qvlat_permute(k,i)             
!         qctend(i,k)                  =    qctend_permute(k,i)            
!         qitend(i,k)                  =    qitend_permute(k,i)            
!         nctend(i,k)                  =    nctend_permute(k,i)            
!         nitend(i,k)                  =    nitend_permute(k,i)            
!                                                                          
!         effc(i,k)                    =    effc_permute(k,i)              
!         effc_fn(i,k)                 =    effc_fn_permute(k,i)           
!         effi(i,k)                    =    effi_permute(k,i)              
!         !prect(i,k)                   =   !prect_permute(k,i)            
!         !preci(i,k)                   =   !preci_permute(k,i)            
!         nevapr(i,k)                  =    nevapr_permute(k,i)            
!         evapsnow(i,k)                =    evapsnow_permute(k,i)          
!         am_evp_st(i,k)               =    am_evp_st_permute(k,i)         
!         prain(i,k)                   =    prain_permute(k,i)             
!         prodsnow(i,k)                =    prodsnow_permute(k,i)          
!         cmeout(i,k)                  =    cmeout_permute(k,i)            
!         deffi(i,k)                   =    deffi_permute(k,i)             
!         pgamrad(i,k)                 =    pgamrad_permute(k,i)           
!         lamcrad(i,k)                 =    lamcrad_permute(k,i)           
!         qsout(i,k)                   =    qsout_permute(k,i)             
!         dsout(i,k)                   =    dsout_permute(k,i)             
!                                                                          
!         qrout(i,k)                   =    qrout_permute(k,i)             
!         reff_rain(i,k)               =    reff_rain_permute(k,i)         
!         reff_snow(i,k)               =    reff_snow_permute(k,i)         
!         qcsevap(i,k)                 =    qcsevap_permute(k,i)           
!         qisevap(i,k)                 =    qisevap_permute(k,i)           
!         qvres(i,k)                   =    qvres_permute(k,i)             
!         cmeiout(i,k)                 =    cmeiout_permute(k,i)           
!         vtrmc(i,k)                   =    vtrmc_permute(k,i)             
!         vtrmi(i,k)                   =    vtrmi_permute(k,i)             
!         qcsedten(i,k)                =    qcsedten_permute(k,i)          
!         qisedten(i,k)                =    qisedten_permute(k,i)          
!         prao(i,k)                    =    prao_permute(k,i)              
!         prco(i,k)                    =    prco_permute(k,i)              
!         mnuccco(i,k)                 =    mnuccco_permute(k,i)           
!         mnuccto(i,k)                 =    mnuccto_permute(k,i)           
!         msacwio(i,k)                 =    msacwio_permute(k,i)           
!         psacwso(i,k)                 =    psacwso_permute(k,i)           
!         bergso(i,k)                  =    bergso_permute(k,i)            
!         bergo(i,k)                   =    bergo_permute(k,i)             
!         melto(i,k)                   =    melto_permute(k,i)             
!         homoo(i,k)                   =    homoo_permute(k,i)             
!         qcreso(i,k)                  =    qcreso_permute(k,i)            
!         prcio(i,k)                   =    prcio_permute(k,i)             
!         praio(i,k)                   =    praio_permute(k,i)             
!         qireso(i,k)                  =    qireso_permute(k,i)            
!         mnuccro(i,k)                 =    mnuccro_permute(k,i)           
!         pracso(i,k)                  =    pracso_permute(k,i)            
!         meltsdt(i,k)                 =    meltsdt_permute(k,i)           
!         frzrdt(i,k)                  =    frzrdt_permute(k,i)            
!         mnuccdo(i,k)                 =    mnuccdo_permute(k,i)           
!         nrout(i,k)                   =    nrout_permute(k,i)             
!         nsout(i,k)                   =    nsout_permute(k,i)             
!         refl(i,k)                    =    refl_permute(k,i)              
!         arefl(i,k)                   =    arefl_permute(k,i)             
!         areflz(i,k)                  =    areflz_permute(k,i)            
!         frefl(i,k)                   =    frefl_permute(k,i)             
!         csrfl(i,k)                   =    csrfl_permute(k,i)             
!         acsrfl(i,k)                  =    acsrfl_permute(k,i)            
!         fcsrfl(i,k)                  =    fcsrfl_permute(k,i)            
!         rercld(i,k)                  =    rercld_permute(k,i)            
!         ncai(i,k)                    =    ncai_permute(k,i)              
!         ncal(i,k)                    =    ncal_permute(k,i)              
!         qrout2(i,k)                  =    qrout2_permute(k,i)            
!         qsout2(i,k)                  =    qsout2_permute(k,i)            
!         nrout2(i,k)                  =    nrout2_permute(k,i)            
!         nsout2(i,k)                  =    nsout2_permute(k,i)            
!         drout2(i,k)                  =    drout2_permute(k,i)            
!         dsout2(i,k)                  =    dsout2_permute(k,i)            
!         freqs(i,k)                   =    freqs_permute(k,i)             
!         freqr(i,k)                   =    freqr_permute(k,i)             
!         nfice(i,k)                   =    nfice_permute(k,i)             
!         prer_evap(i,k)               =    prer_evap_permute(k,i)         
!         !nevapr2(i,k)                 =    nevapr2_permute(k,i)           
!         !                                                                 
!         !tnd_qsnow(i,k)               =    tnd_qsnow_permute(k,i)         
!         !tnd_nsnow(i,k)               =    tnd_nsnow_permute(k,i)         
!         !re_ice(i,k)                  =    re_ice_permute(k,i)            
!         !frzimm(i,k)                  =    frzimm_permute(k,i)            
!         !frzcnt(i,k)                  =    frzcnt_permute(k,i)            
!         !frzdep(i,k)                  =    frzdep_permute(k,i)            
!         !                                                                 
!         !t1(i,k)                      =    t1_permute(k,i)                
!         !q1(i,k)                      =    q1_permute(k,i)                
!         !qc1(i,k)                     =    qc1_permute(k,i)               
!         !qi1(i,k)                     =    qi1_permute(k,i)               
!         !nc1(i,k)                     =    nc1_permute(k,i)               
!         !ni1(i,k)                     =    ni1_permute(k,i)               
!         !tlat1(i,k)                   =    tlat1_permute(k,i)             
!         !qvlat1(i,k)                  =    qvlat1_permute(k,i)            
!         !qctend1(i,k)                 =    qctend1_permute(k,i)           
!         !qitend1(i,k)                 =    qitend1_permute(k,i)           
!         !nctend1(i,k)                 =    nctend1_permute(k,i)           
!         !nitend1(i,k)                 =    nitend1_permute(k,i)           
!         !q(i,k)                       =    q_permute(k,i)                 
!         !t(i,k)                       =    t_permute(k,i)                 
!         !rho(i,k)                     =    rho_permute(k,i)               
!         !dv(i,k)                      =    dv_permute(k,i)                
!         !mu(i,k)                      =    mu_permute(k,i)                
!         !sc(i,k)                      =    sc_permute(k,i)                
!         !kap(i,k)                     =    kap_permute(k,i)               
!         !rhof(i,k)                    =    rhof_permute(k,i)              
!         !cldmax(i,k)                  =    cldmax_permute(k,i)            
!         !cldm(i,k)                    =    cldm_permute(k,i)              
!         !icldm(i,k)                   =    icldm_permute(k,i)             
!         !lcldm(i,k)                   =    lcldm_permute(k,i)             
!         !cme(i,k)                     =    cme_permute(k,i)               
!         !cmei(i,k)                    =    cmei_permute(k,i)              
!         !cwml(i,k)                    =    cwml_permute(k,i)              
!         !cwmi(i,k)                    =    cwmi_permute(k,i)              
!         !lcldn(i,k)                   =    lcldn_permute(k,i)             
!         !lcldo(i,k)                   =    lcldo_permute(k,i)             
!         !nctend_mixnuc(i,k)           =    nctend_mixnuc_permute(k,i)     
!         !                                                                 
!         !qcic(i,k)                    =    qcic_permute(k,i)              
!         !qiic(i,k)                    =    qiic_permute(k,i)              
!         !qniic(i,k)                   =    qniic_permute(k,i)             
!         !qric(i,k)                    =    qric_permute(k,i)              
!         !ncic(i,k)                    =    ncic_permute(k,i)              
!         !niic(i,k)                    =    niic_permute(k,i)              
!         !nsic(i,k)                    =    nsic_permute(k,i)              
!         !nric(i,k)                    =    nric_permute(k,i)              
!         !arcld(i,k)                   =    arcld_permute(k,i)             
!         !dumc(i,k)                    =    dumc_permute(k,i)              
!         !dumnc(i,k)                   =    dumnc_permute(k,i)             
!         !dumi(i,k)                    =    dumi_permute(k,i)              
!         !dumni(i,k)                   =    dumni_permute(k,i)             
!         !dums(i,k)                    =    dums_permute(k,i)              
!         !dumns(i,k)                   =    dumns_permute(k,i)             
!         !dumr(i,k)                    =    dumr_permute(k,i)              
!         !dumnr(i,k)                   =    dumnr_permute(k,i)             
!         !relhum(i,k)                  =    relhum_permute(k,i)            
!         !arn(i,k)                     =    arn_permute(k,i)               
!         !asn(i,k)                     =    asn_permute(k,i)               
!         !acn(i,k)                     =    acn_permute(k,i)               
!         !ain(i,k)                     =    ain_permute(k,i)               
!         !dz(i,k)                      =    dz_permute(k,i)                
!         !                                                                 
!         !tsp(i,k)                     =    tsp_permute(k,i)               
!         !qsp(i,k)                     =    qsp_permute(k,i)               
!         !qsphy(i,k)                   =    qsphy_permute(k,i)             
!         !esl(i,k)                     =    esl_permute(k,i)               
!         !esi(i,k)                     =    esi_permute(k,i)               
!         !qnitend(i,k)                 =    qnitend_permute(k,i)           
!         !nstend(i,k)                  =    nstend_permute(k,i)            
!         !qrtend(i,k)                  =    qrtend_permute(k,i)            
!         !nrtend(i,k)                  =    nrtend_permute(k,i)            
!         !berg(i,k)                    =    berg_permute(k,i)              
!         !drout(i,k)                   =    drout_permute(k,i)             
!         !dum2i(i,k)                   =    dum2i_permute(k,i)             
!         !dum2l(i,k)                   =    dum2l_permute(k,i)             
!         !cldmw(i,k)                   =    cldmw_permute(k,i)             
!         !rainrt(i,k)                  =    rainrt_permute(k,i)            
!         !rainrt1(i,k)                 =    rainrt1_permute(k,i)           
!         !                                                                 
!         !rndst(i, k, :)               =    rndst_permute(k, i, :)         
!         !nacon(i, k, :)               =    nacon_permute(k, i, :)           
!
!    enddo
!enddo
!
!do i=1,pcols
!    do k=1,pver+1
!        
!        !rflx1(i,k)  =  rflx1_permute(k,i)
!        !sflx1(i,k)  =  sflx1_permute(k,i)
!        rflx(i,k)   =  rflx_permute(k,i) 
!        sflx(i,k)   =  sflx_permute(k,i) 
!    enddo
!enddo

end subroutine micro_mg_tend

!========================================================================
!UTILITIES
!========================================================================

#else
!===============================================================================
!microphysics routine for each timestep goes here...

subroutine micro_mg_tend ( &
     microp_uniform, pcols, pver, ncol, top_lev, deltatin,&
     tn, qn, qc, qi, nc,                              &
     ni, p, pdel, cldn, liqcldf,                      &
     relvar, accre_enhan,                             &
     icecldf, rate1ord_cw2pr_st, naai, npccnin,       &
     rndst, nacon, tlat, qvlat, qctend,               &
     qitend, nctend, nitend, effc, effc_fn,           &
     effi, prect, preci, nevapr, evapsnow, am_evp_st, &
     prain, prodsnow, cmeout, deffi, pgamrad,         &
     lamcrad, qsout, dsout, rflx, sflx,               &
     qrout, reff_rain, reff_snow, qcsevap, qisevap,   &
     qvres, cmeiout, vtrmc, vtrmi, qcsedten,          &
     qisedten, prao, prco, mnuccco, mnuccto,          &
     msacwio, psacwso, bergso, bergo, melto,          &
     homoo, qcreso, prcio, praio, qireso,             &
     mnuccro, pracso, meltsdt, frzrdt, mnuccdo,       &
     nrout, nsout, refl, arefl, areflz,               &
     frefl, csrfl, acsrfl, fcsrfl, rercld,            &
     ncai, ncal, qrout2, qsout2, nrout2,              &
     nsout2, drout2, dsout2, freqs, freqr,            &
     nfice, prer_evap, do_cldice, errstring,          &
     tnd_qsnow, tnd_nsnow, re_ice,                    &
     frzimm, frzcnt, frzdep)

! input arguments
logical,  intent(in) :: microp_uniform  ! True = configure uniform for sub-columns  False = use w/o sub-columns (standard)
integer,  intent(in) :: pcols                ! size of column (first) index
integer,  intent(in) :: pver                 ! number of layers in columns
integer,  intent(in) :: ncol                 ! number of columns
integer,  intent(in) :: top_lev              ! top level microphys is applied
real(r8), intent(in) :: deltatin             ! time step (s)
real(r8), intent(in) :: tn(pcols,pver)       ! input temperature (K)
real(r8), intent(in) :: qn(pcols,pver)       ! input h20 vapor mixing ratio (kg/kg)
real(r8), intent(in) :: relvar(pcols,pver)   ! relative variance of cloud water (-)
real(r8), intent(in) :: accre_enhan(pcols,pver) ! optional accretion enhancement factor (-)

! note: all input cloud variables are grid-averaged
real(r8), intent(inout) :: qc(pcols,pver)    ! cloud water mixing ratio (kg/kg)
real(r8), intent(inout) :: qi(pcols,pver)    ! cloud ice mixing ratio (kg/kg)
real(r8), intent(inout) :: nc(pcols,pver)    ! cloud water number conc (1/kg)
real(r8), intent(inout) :: ni(pcols,pver)    ! cloud ice number conc (1/kg)
real(r8), intent(in) :: p(pcols,pver)        ! air pressure (pa)
real(r8), intent(in) :: pdel(pcols,pver)     ! pressure difference across level (pa)
real(r8), intent(in) :: cldn(pcols,pver)     ! cloud fraction
real(r8), intent(in) :: icecldf(pcols,pver)  ! ice cloud fraction   
real(r8), intent(in) :: liqcldf(pcols,pver)  ! liquid cloud fraction
          
real(r8), intent(out) :: rate1ord_cw2pr_st(pcols,pver) ! 1st order rate for direct cw to precip conversion 
! used for scavenging
! Inputs for aerosol activation
real(r8), intent(in) :: naai(pcols,pver)      ! ice nulceation number (from microp_aero_ts) 
real(r8), intent(in) :: npccnin(pcols,pver)   ! ccn activated number tendency (from microp_aero_ts)
real(r8), intent(in) :: rndst(pcols,pver,4)   ! radius of 4 dust bins for contact freezing (from microp_aero_ts)
real(r8), intent(in) :: nacon(pcols,pver,4)   ! number in 4 dust bins for contact freezing  (from microp_aero_ts)

! Used with CARMA cirrus microphysics
! (or similar external microphysics model)
logical,  intent(in) :: do_cldice             ! Prognosing cldice

! output arguments

real(r8), intent(out) :: tlat(pcols,pver)    ! latent heating rate       (W/kg)
real(r8), intent(out) :: qvlat(pcols,pver)   ! microphysical tendency qv (1/s)
real(r8), intent(out) :: qctend(pcols,pver)  ! microphysical tendency qc (1/s) 
real(r8), intent(out) :: qitend(pcols,pver)  ! microphysical tendency qi (1/s)
real(r8), intent(out) :: nctend(pcols,pver)  ! microphysical tendency nc (1/(kg*s))
real(r8), intent(out) :: nitend(pcols,pver)  ! microphysical tendency ni (1/(kg*s))
real(r8), intent(out) :: effc(pcols,pver)    ! droplet effective radius (micron)
real(r8), intent(out) :: effc_fn(pcols,pver) ! droplet effective radius, assuming nc = 1.e8 kg-1
real(r8), intent(out) :: effi(pcols,pver)    ! cloud ice effective radius (micron)
real(r8), intent(out) :: prect(pcols)        ! surface precip rate (m/s)
real(r8), intent(out) :: preci(pcols)        ! cloud ice/snow precip rate (m/s)
real(r8), intent(out) :: nevapr(pcols,pver)  ! evaporation rate of rain + snow
real(r8), intent(out) :: evapsnow(pcols,pver)! sublimation rate of snow
real(r8), intent(out) :: am_evp_st(pcols,pver)! stratiform evaporation area
real(r8), intent(out) :: prain(pcols,pver)   ! production of rain + snow
real(r8), intent(out) :: prodsnow(pcols,pver)! production of snow
real(r8), intent(out) :: cmeout(pcols,pver)  ! evap/sub of cloud
real(r8), intent(out) :: deffi(pcols,pver)   ! ice effective diameter for optics (radiation)
real(r8), intent(out) :: pgamrad(pcols,pver) ! ice gamma parameter for optics (radiation)
real(r8), intent(out) :: lamcrad(pcols,pver) ! slope of droplet distribution for optics (radiation)
real(r8), intent(out) :: qsout(pcols,pver)   ! snow mixing ratio (kg/kg)
real(r8), intent(out) :: dsout(pcols,pver)   ! snow diameter (m)
real(r8), intent(out) :: rflx(pcols,pver+1)  ! grid-box average rain flux (kg m^-2 s^-1)
real(r8), intent(out) :: sflx(pcols,pver+1)  ! grid-box average snow flux (kg m^-2 s^-1)
real(r8), intent(out) :: qrout(pcols,pver)     ! grid-box average rain mixing ratio (kg/kg)
real(r8), intent(inout) :: reff_rain(pcols,pver) ! rain effective radius (micron)
real(r8), intent(inout) :: reff_snow(pcols,pver) ! snow effective radius (micron)
real(r8), intent(out) :: qcsevap(pcols,pver) ! cloud water evaporation due to sedimentation
real(r8), intent(out) :: qisevap(pcols,pver) ! cloud ice sublimation due to sublimation
real(r8), intent(out) :: qvres(pcols,pver) ! residual condensation term to ensure RH < 100%
real(r8), intent(out) :: cmeiout(pcols,pver) ! grid-mean cloud ice sub/dep
real(r8), intent(out) :: vtrmc(pcols,pver) ! mass-weighted cloud water fallspeed
real(r8), intent(out) :: vtrmi(pcols,pver) ! mass-weighted cloud ice fallspeed
real(r8), intent(out) :: qcsedten(pcols,pver) ! qc sedimentation tendency
real(r8), intent(out) :: qisedten(pcols,pver) ! qi sedimentation tendency
! microphysical process rates for output (mixing ratio tendencies)
real(r8), intent(out) :: prao(pcols,pver) ! accretion of cloud by rain 
real(r8), intent(out) :: prco(pcols,pver) ! autoconversion of cloud to rain
real(r8), intent(out) :: mnuccco(pcols,pver) ! mixing rat tend due to immersion freezing
real(r8), intent(out) :: mnuccto(pcols,pver) ! mixing ratio tend due to contact freezing
real(r8), intent(out) :: msacwio(pcols,pver) ! mixing ratio tend due to H-M splintering
real(r8), intent(out) :: psacwso(pcols,pver) ! collection of cloud water by snow
real(r8), intent(out) :: bergso(pcols,pver) ! bergeron process on snow
real(r8), intent(out) :: bergo(pcols,pver) ! bergeron process on cloud ice
real(r8), intent(out) :: melto(pcols,pver) ! melting of cloud ice
real(r8), intent(out) :: homoo(pcols,pver) ! homogeneos freezign cloud water
real(r8), intent(out) :: qcreso(pcols,pver) ! residual cloud condensation due to removal of excess supersat
real(r8), intent(out) :: prcio(pcols,pver) ! autoconversion of cloud ice to snow
real(r8), intent(out) :: praio(pcols,pver) ! accretion of cloud ice by snow
real(r8), intent(out) :: qireso(pcols,pver) ! residual ice deposition due to removal of excess supersat
real(r8), intent(out) :: mnuccro(pcols,pver) ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
real(r8), intent(out) :: pracso (pcols,pver) ! mixing ratio tendency due to accretion of rain by snow (1/s)
real(r8), intent(out) :: meltsdt(pcols,pver) ! latent heating rate due to melting of snow  (W/kg)
real(r8), intent(out) :: frzrdt (pcols,pver) ! latent heating rate due to homogeneous freezing of rain (W/kg)
real(r8), intent(out) :: mnuccdo(pcols,pver) ! mass tendency from ice nucleation
real(r8), intent(out) :: nrout(pcols,pver) ! rain number concentration (1/m3)
real(r8), intent(out) :: nsout(pcols,pver) ! snow number concentration (1/m3)
real(r8), intent(out) :: refl(pcols,pver)    ! analytic radar reflectivity        
real(r8), intent(out) :: arefl(pcols,pver)  !average reflectivity will zero points outside valid range
real(r8), intent(out) :: areflz(pcols,pver)  !average reflectivity in z.
real(r8), intent(out) :: frefl(pcols,pver)
real(r8), intent(out) :: csrfl(pcols,pver)   !cloudsat reflectivity 
real(r8), intent(out) :: acsrfl(pcols,pver)  !cloudsat average
real(r8), intent(out) :: fcsrfl(pcols,pver)
real(r8), intent(out) :: rercld(pcols,pver) ! effective radius calculation for rain + cloud
real(r8), intent(out) :: ncai(pcols,pver) ! output number conc of ice nuclei available (1/m3)
real(r8), intent(out) :: ncal(pcols,pver) ! output number conc of CCN (1/m3)
real(r8), intent(out) :: qrout2(pcols,pver)
real(r8), intent(out) :: qsout2(pcols,pver)
real(r8), intent(out) :: nrout2(pcols,pver)
real(r8), intent(out) :: nsout2(pcols,pver)
real(r8), intent(out) :: drout2(pcols,pver) ! mean rain particle diameter (m)
real(r8), intent(out) :: dsout2(pcols,pver) ! mean snow particle diameter (m)
real(r8), intent(out) :: freqs(pcols,pver)
real(r8), intent(out) :: freqr(pcols,pver)
real(r8), intent(out) :: nfice(pcols,pver)
real(r8), intent(out) :: prer_evap(pcols,pver)

real(r8) :: nevapr2(pcols,pver)

character(128),   intent(out) :: errstring       ! Output status (non-blank for error return)

! Tendencies calculated by external schemes that can replace MG's native
! process tendencies.

! Used with CARMA cirrus microphysics
! (or similar external microphysics model)
real(r8), intent(in), pointer :: tnd_qsnow(:,:) ! snow mass tendency (kg/kg/s)
real(r8), intent(in), pointer :: tnd_nsnow(:,:) ! snow number tendency (#/kg/s)
real(r8), intent(in), pointer :: re_ice(:,:)    ! ice effective radius (m)

! From external ice nucleation.
real(r8), intent(in), pointer :: frzimm(:,:) ! Number tendency due to immersion freezing (1/cm3)
real(r8), intent(in), pointer :: frzcnt(:,:) ! Number tendency due to contact freezing (1/cm3)
real(r8), intent(in), pointer :: frzdep(:,:) ! Number tendency due to deposition nucleation (1/cm3)

! local workspace
! all units mks unless otherwise stated

! Additional constants to help speed up code
real(r8) :: cons2
real(r8) :: cons3
real(r8) :: cons9
real(r8) :: cons10
real(r8) :: cons12
real(r8) :: cons15
real(r8) :: cons18
real(r8) :: cons19
real(r8) :: cons20

! temporary variables for sub-stepping 
real(r8) :: t1(pcols,pver)
real(r8) :: q1(pcols,pver)
real(r8) :: qc1(pcols,pver)
real(r8) :: qi1(pcols,pver)
real(r8) :: nc1(pcols,pver)
real(r8) :: ni1(pcols,pver)
real(r8) :: tlat1(pcols,pver)
real(r8) :: qvlat1(pcols,pver)
real(r8) :: qctend1(pcols,pver)
real(r8) :: qitend1(pcols,pver)
real(r8) :: nctend1(pcols,pver)
real(r8) :: nitend1(pcols,pver)
real(r8) :: prect1(pcols)
real(r8) :: preci1(pcols)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(r8) :: deltat        ! sub-time step (s)
real(r8) :: omsm    ! number near unity for round-off issues
real(r8) :: dto2    ! dt/2 (s)
real(r8) :: mincld  ! minimum allowed cloud fraction
real(r8) :: q(pcols,pver) ! water vapor mixing ratio (kg/kg)
real(r8) :: t(pcols,pver) ! temperature (K)
real(r8) :: rho(pcols,pver) ! air density (kg m-3)
real(r8) :: dv(pcols,pver)  ! diffusivity of water vapor in air
real(r8) :: mu(pcols,pver)  ! viscocity of air
real(r8) :: sc(pcols,pver)  ! schmidt number
real(r8) :: kap(pcols,pver) ! thermal conductivity of air
real(r8) :: rhof(pcols,pver) ! air density correction factor for fallspeed
real(r8) :: cldmax(pcols,pver) ! precip fraction assuming maximum overlap
real(r8) :: cldm(pcols,pver)   ! cloud fraction
real(r8) :: icldm(pcols,pver)   ! ice cloud fraction
real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction
real(r8) :: icwc(pcols)    ! in cloud water content (liquid+ice)
real(r8) :: calpha(pcols)  ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cbeta(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cbetah(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cgamma(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cgamah(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: rcgama(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec1(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec2(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec3(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec4(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: qtmp ! dummy qv 
real(r8) :: dum  ! temporary dummy variable

real(r8) :: cme(pcols,pver)  ! total (liquid+ice) cond/evap rate of cloud

real(r8) :: cmei(pcols,pver) ! dep/sublimation rate of cloud ice
real(r8) :: cwml(pcols,pver) ! cloud water mixing ratio
real(r8) :: cwmi(pcols,pver) ! cloud ice mixing ratio
real(r8) :: nnuccd(pver)   ! ice nucleation rate from deposition/cond.-freezing
real(r8) :: mnuccd(pver)   ! mass tendency from ice nucleation
real(r8) :: qcld              ! total cloud water
real(r8) :: lcldn(pcols,pver) ! fractional coverage of new liquid cloud
real(r8) :: lcldo(pcols,pver) ! fractional coverage of old liquid cloud
real(r8) :: nctend_mixnuc(pcols,pver)
real(r8) :: arg ! argument of erfc

! for calculation of rate1ord_cw2pr_st
real(r8) :: qcsinksum_rate1ord(pver)   ! sum over iterations of cw to precip sink
real(r8) :: qcsum_rate1ord(pver)    ! sum over iterations of cloud water       

real(r8) :: alpha

real(r8) :: dum1,dum2   !general dummy variables

real(r8) :: npccn(pver)     ! droplet activation rate
real(r8) :: qcic(pcols,pver) ! in-cloud cloud liquid mixing ratio
real(r8) :: qiic(pcols,pver) ! in-cloud cloud ice mixing ratio
real(r8) :: qniic(pcols,pver) ! in-precip snow mixing ratio
real(r8) :: qric(pcols,pver) ! in-precip rain mixing ratio
real(r8) :: ncic(pcols,pver) ! in-cloud droplet number conc
real(r8) :: niic(pcols,pver) ! in-cloud cloud ice number conc
real(r8) :: nsic(pcols,pver) ! in-precip snow number conc
real(r8) :: nric(pcols,pver) ! in-precip rain number conc
real(r8) :: lami(pver) ! slope of cloud ice size distr
real(r8) :: n0i(pver) ! intercept of cloud ice size distr
real(r8) :: lamc(pver) ! slope of cloud liquid size distr
real(r8) :: n0c(pver) ! intercept of cloud liquid size distr
real(r8) :: lams(pver) ! slope of snow size distr
real(r8) :: n0s(pver) ! intercept of snow size distr
real(r8) :: lamr(pver) ! slope of rain size distr
real(r8) :: n0r(pver) ! intercept of rain size distr
real(r8) :: cdist1(pver) ! size distr parameter to calculate droplet freezing
! combined size of precip & cloud drops
real(r8) :: arcld(pcols,pver) ! averaging control flag
real(r8) :: Actmp  !area cross section of drops
real(r8) :: Artmp  !area cross section of rain

real(r8) :: pgam(pver) ! spectral width parameter of droplet size distr
real(r8) :: lammax  ! maximum allowed slope of size distr
real(r8) :: lammin  ! minimum allowed slope of size distr
real(r8) :: nacnt   ! number conc of contact ice nuclei
real(r8) :: mnuccc(pver) ! mixing ratio tendency due to freezing of cloud water
real(r8) :: nnuccc(pver) ! number conc tendency due to freezing of cloud water

real(r8) :: mnucct(pver) ! mixing ratio tendency due to contact freezing of cloud water
real(r8) :: nnucct(pver) ! number conc tendency due to contact freezing of cloud water
real(r8) :: msacwi(pver) ! mixing ratio tendency due to HM ice multiplication
real(r8) :: nsacwi(pver) ! number conc tendency due to HM ice multiplication

real(r8) :: prc(pver) ! qc tendency due to autoconversion of cloud droplets
real(r8) :: nprc(pver) ! number conc tendency due to autoconversion of cloud droplets
real(r8) :: nprc1(pver) ! qr tendency due to autoconversion of cloud droplets
real(r8) :: nsagg(pver) ! ns tendency due to self-aggregation of snow
real(r8) :: dc0  ! mean size droplet size distr
real(r8) :: ds0  ! mean size snow size distr (area weighted)
real(r8) :: eci  ! collection efficiency for riming of snow by droplets
real(r8) :: psacws(pver) ! mixing rat tendency due to collection of droplets by snow
real(r8) :: npsacws(pver) ! number conc tendency due to collection of droplets by snow
real(r8) :: uni ! number-weighted cloud ice fallspeed
real(r8) :: umi ! mass-weighted cloud ice fallspeed
real(r8) :: uns(pver) ! number-weighted snow fallspeed
real(r8) :: ums(pver) ! mass-weighted snow fallspeed
real(r8) :: unr(pver) ! number-weighted rain fallspeed
real(r8) :: umr(pver) ! mass-weighted rain fallspeed
real(r8) :: unc ! number-weighted cloud droplet fallspeed
real(r8) :: umc ! mass-weighted cloud droplet fallspeed
real(r8) :: pracs(pver) ! mixing rat tendency due to collection of rain by snow
real(r8) :: npracs(pver) ! number conc tendency due to collection of rain by snow
real(r8) :: mnuccr(pver) ! mixing rat tendency due to freezing of rain
real(r8) :: nnuccr(pver) ! number conc tendency due to freezing of rain
real(r8) :: pra(pver) ! mixing rat tendnency due to accretion of droplets by rain
real(r8) :: npra(pver) ! nc tendnency due to accretion of droplets by rain
real(r8) :: nragg(pver) ! nr tendency due to self-collection of rain
real(r8) :: prci(pver) ! mixing rat tendency due to autoconversion of cloud ice to snow
real(r8) :: nprci(pver) ! number conc tendency due to autoconversion of cloud ice to snow
real(r8) :: prai(pver) ! mixing rat tendency due to accretion of cloud ice by snow
real(r8) :: nprai(pver) ! number conc tendency due to accretion of cloud ice by snow
real(r8) :: qvs ! liquid saturation vapor mixing ratio
real(r8) :: qvi ! ice saturation vapor mixing ratio
real(r8) :: dqsdt ! change of sat vapor mixing ratio with temperature
real(r8) :: dqsidt ! change of ice sat vapor mixing ratio with temperature
real(r8) :: ab ! correction factor for rain evap to account for latent heat
real(r8) :: qclr ! water vapor mixing ratio in clear air
real(r8) :: abi ! correction factor for snow sublimation to account for latent heat
real(r8) :: epss ! 1/ sat relaxation timescale for snow
real(r8) :: epsr ! 1/ sat relaxation timescale for rain
real(r8) :: pre(pver) ! rain mixing rat tendency due to evaporation
real(r8) :: prds(pver) ! snow mixing rat tendency due to sublimation
real(r8) :: qce ! dummy qc for conservation check
real(r8) :: qie ! dummy qi for conservation check
real(r8) :: nce ! dummy nc for conservation check
real(r8) :: nie ! dummy ni for conservation check
real(r8) :: ratio ! parameter for conservation check
real(r8) :: dumc(pcols,pver) ! dummy in-cloud qc
real(r8) :: dumnc(pcols,pver) ! dummy in-cloud nc
real(r8) :: dumi(pcols,pver) ! dummy in-cloud qi
real(r8) :: dumni(pcols,pver) ! dummy in-cloud ni
real(r8) :: dums(pcols,pver) ! dummy in-cloud snow mixing rat
real(r8) :: dumns(pcols,pver) ! dummy in-cloud snow number conc
real(r8) :: dumr(pcols,pver) ! dummy in-cloud rain mixing rat
real(r8) :: dumnr(pcols,pver) ! dummy in-cloud rain number conc
! below are parameters for cloud water and cloud ice sedimentation calculations
real(r8) :: fr(pver)
real(r8) :: fnr(pver)
real(r8) :: fc(pver)
real(r8) :: fnc(pver)
real(r8) :: fi(pver)
real(r8) :: fni(pver)
real(r8) :: fs(pver)
real(r8) :: fns(pver)
real(r8) :: faloutr(pver)
real(r8) :: faloutnr(pver)
real(r8) :: faloutc(pver)
real(r8) :: faloutnc(pver)
real(r8) :: falouti(pver)
real(r8) :: faloutni(pver)
real(r8) :: falouts(pver)
real(r8) :: faloutns(pver)
real(r8) :: faltndr
real(r8) :: faltndnr
real(r8) :: faltndc
real(r8) :: faltndnc
real(r8) :: faltndi
real(r8) :: faltndni
real(r8) :: faltnds
real(r8) :: faltndns
real(r8) :: faltndqie
real(r8) :: faltndqce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(r8) :: relhum(pcols,pver) ! relative humidity
real(r8) :: csigma(pcols) ! parameter for cond/evap of cloud water/ice
real(r8) :: rgvm ! max fallspeed for all species
real(r8) :: arn(pcols,pver) ! air density corrected rain fallspeed parameter
real(r8) :: asn(pcols,pver) ! air density corrected snow fallspeed parameter
real(r8) :: acn(pcols,pver) ! air density corrected cloud droplet fallspeed parameter
real(r8) :: ain(pcols,pver) ! air density corrected cloud ice fallspeed parameter
real(r8) :: nsubi(pver) ! evaporation of cloud ice number
real(r8) :: nsubc(pver) ! evaporation of droplet number
real(r8) :: nsubs(pver) ! evaporation of snow number
real(r8) :: nsubr(pver) ! evaporation of rain number
real(r8) :: mtime ! factor to account for droplet activation timescale
real(r8) :: dz(pcols,pver) ! height difference across model vertical level


!! add precip flux variables for sub-stepping
real(r8) :: rflx1(pcols,pver+1)
real(r8) :: sflx1(pcols,pver+1)

! returns from function/subroutine calls
real(r8) :: tsp(pcols,pver)      ! saturation temp (K)
real(r8) :: qsp(pcols,pver)      ! saturation mixing ratio (kg/kg)
real(r8) :: qsphy(pcols,pver)      ! saturation mixing ratio (kg/kg): hybrid rh
real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
real(r8) :: esl(pcols,pver)      ! liquid sat vapor pressure (pa)
real(r8) :: esi(pcols,pver)      ! ice sat vapor pressure (pa)

! sum of source/sink terms for diagnostic precip

real(r8) :: qnitend(pcols,pver) ! snow mixing ratio source/sink term
real(r8) :: nstend(pcols,pver)  ! snow number concentration source/sink term
real(r8) :: qrtend(pcols,pver) ! rain mixing ratio source/sink term
real(r8) :: nrtend(pcols,pver)  ! rain number concentration source/sink term
real(r8) :: qrtot ! vertically-integrated rain mixing rat source/sink term
real(r8) :: nrtot ! vertically-integrated rain number conc source/sink term
real(r8) :: qstot ! vertically-integrated snow mixing rat source/sink term
real(r8) :: nstot ! vertically-integrated snow number conc source/sink term

! new terms for Bergeron process

real(r8) :: dumnnuc ! provisional ice nucleation rate (for calculating bergeron)
real(r8) :: ninew  ! provisional cloud ice number conc (for calculating bergeron)
real(r8) :: qinew ! provisional cloud ice mixing ratio (for calculating bergeron)
real(r8) :: qvl  ! liquid sat mixing ratio   
real(r8) :: epsi ! 1/ sat relaxation timecale for cloud ice
real(r8) :: prd ! provisional deposition rate of cloud ice at water sat 
real(r8) :: berg(pcols,pver) ! mixing rat tendency due to bergeron process for cloud ice
real(r8) :: bergs(pver) ! mixing rat tendency due to bergeron process for snow

!bergeron terms
real(r8) :: bergtsf   !bergeron timescale to remove all liquid
real(r8) :: rhin      !modified RH for vapor deposition

! diagnostic rain/snow for output to history
! values are in-precip (local) !!!!

real(r8) :: drout(pcols,pver)     ! rain diameter (m)

!averageed rain/snow for history
real(r8) :: dumfice

!ice nucleation, droplet activation
real(r8) :: dum2i(pcols,pver) ! number conc of ice nuclei available (1/kg)
real(r8) :: dum2l(pcols,pver) ! number conc of CCN (1/kg)
real(r8) :: ncmax
real(r8) :: nimax

real(r8) :: qcvar     ! 1/relative variance of sub-grid qc

! loop array variables
integer i,k,nstep,n, l
integer ii,kk, m

! loop variables for sub-step solution
integer iter,it,ltrue(pcols)

! used in contact freezing via dust particles
real(r8)  tcnt, viscosity, mfp
real(r8)  slip1, slip2, slip3, slip4
!        real(r8)  dfaer1, dfaer2, dfaer3, dfaer4
!        real(r8)  nacon1,nacon2,nacon3,nacon4
real(r8)  ndfaer1, ndfaer2, ndfaer3, ndfaer4
real(r8)  nslip1, nslip2, nslip3, nslip4

! used in ice effective radius
real(r8)  bbi, cci, ak, iciwc, rvi

! used in Bergeron processe and water vapor deposition
real(r8)  Tk, deles, Aprpr, Bprpr, Cice, qi0, Crate, qidep

! mean cloud fraction over the time step
real(r8)  cldmw(pcols,pver)

! used in secondary ice production
real(r8) ni_secp

! variabels to check for RH after rain evap

real(r8) :: esn
real(r8) :: qsn
real(r8) :: ttmp



real(r8) :: rainrt(pcols,pver)  ! rain rate for reflectivity calculation
real(r8) :: rainrt1(pcols,pver)
real(r8) :: tmp

real(r8) dmc,ssmc,dstrn  ! variables for modal scheme.

real(r8), parameter :: cdnl    = 0.e6_r8    ! cloud droplet number limiter

! heterogeneous freezing
real(r8) :: mnudep(pver) ! mixing ratio tendency due to deposition of water vapor
real(r8) :: nnudep(pver) ! number conc tendency due to deposition of water vapor
real(r8) :: con1 ! work cnstant
real(r8) :: r3lx ! Mean volume radius (m)
real(r8) :: mi0l
real(r8) :: frztmp

logical  :: do_clubb_sgs

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Return error message
errstring = ' '

if (.not. (do_cldice .or. &
     (associated(tnd_qsnow) .and. associated(tnd_nsnow) .and. associated(re_ice)))) then
   errstring = "MG's native cloud ice processes are disabled, but &
        &no replacement values were passed in."
end if

if (use_hetfrz_classnuc .and. (.not. &
     (associated(frzimm) .and. associated(frzcnt) .and. associated(frzdep)))) then
   errstring = "Hoose heterogeneous freezing is enabled, but the &
        &required tendencies were not all passed in."
end if

call phys_getopts(do_clubb_sgs_out = do_clubb_sgs)

! initialize  output fields for number conc qand ice nucleation
ncai(1:ncol,1:pver)=0._r8 
ncal(1:ncol,1:pver)=0._r8  

!Initialize rain size
rercld(1:ncol,1:pver)=0._r8
arcld(1:ncol,1:pver)=0._r8

!initialize radiation output variables
pgamrad(1:ncol,1:pver)=0._r8 ! liquid gamma parameter for optics (radiation)
lamcrad(1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics (radiation)
deffi  (1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics (radiation)
!initialize radiation output variables
!initialize water vapor tendency term output
qcsevap(1:ncol,1:pver)=0._r8 
qisevap(1:ncol,1:pver)=0._r8 
qvres  (1:ncol,1:pver)=0._r8 
cmeiout (1:ncol,1:pver)=0._r8
vtrmc (1:ncol,1:pver)=0._r8
vtrmi (1:ncol,1:pver)=0._r8
qcsedten (1:ncol,1:pver)=0._r8
qisedten (1:ncol,1:pver)=0._r8    

prao(1:ncol,1:pver)=0._r8 
prco(1:ncol,1:pver)=0._r8 
mnuccco(1:ncol,1:pver)=0._r8 
mnuccto(1:ncol,1:pver)=0._r8 
msacwio(1:ncol,1:pver)=0._r8 
psacwso(1:ncol,1:pver)=0._r8 
bergso(1:ncol,1:pver)=0._r8 
bergo(1:ncol,1:pver)=0._r8 
melto(1:ncol,1:pver)=0._r8 
homoo(1:ncol,1:pver)=0._r8 
qcreso(1:ncol,1:pver)=0._r8 
prcio(1:ncol,1:pver)=0._r8 
praio(1:ncol,1:pver)=0._r8 
qireso(1:ncol,1:pver)=0._r8 
mnuccro(1:ncol,1:pver)=0._r8 
pracso (1:ncol,1:pver)=0._r8 
meltsdt(1:ncol,1:pver)=0._r8
frzrdt (1:ncol,1:pver)=0._r8
mnuccdo(1:ncol,1:pver)=0._r8

rflx(:,:)=0._r8
sflx(:,:)=0._r8
effc(:,:)=0._r8
effc_fn(:,:)=0._r8
effi(:,:)=0._r8

! assign variable deltat for sub-stepping...
deltat=deltatin

! parameters for scheme

omsm=0.99999_r8
dto2=0.5_r8*deltat
mincld=0.0001_r8

! initialize multi-level fields
q(1:ncol,1:pver)=qn(1:ncol,1:pver)
t(1:ncol,1:pver)=tn(1:ncol,1:pver)

! initialize time-varying parameters

do k=1,pver
   do i=1,ncol
      rho(i,k)=p(i,k)/(r*t(i,k))
      dv(i,k) = 8.794E-5_r8*t(i,k)**1.81_r8/p(i,k)
      mu(i,k) = 1.496E-6_r8*t(i,k)**1.5_r8/(t(i,k)+120._r8) 
      sc(i,k) = mu(i,k)/(rho(i,k)*dv(i,k))
      kap(i,k) = 1.414e3_r8*1.496e-6_r8*t(i,k)**1.5_r8/(t(i,k)+120._r8) 

      ! air density adjustment for fallspeed parameters
      ! includes air density correction factor to the
      ! power of 0.54 following Heymsfield and Bansemer 2007

      rhof(i,k)=(rhosu/rho(i,k))**0.54_r8

      arn(i,k)=ar*rhof(i,k)
      asn(i,k)=as*rhof(i,k)
      acn(i,k)=ac*rhof(i,k)
      ain(i,k)=ai*rhof(i,k)

      ! get dz from dp and hydrostatic approx
      ! keep dz positive (define as layer k-1 - layer k)

      dz(i,k)= pdel(i,k)/(rho(i,k)*g)

   end do
end do

! initialization
qc(1:ncol,1:top_lev-1) = 0._r8
qi(1:ncol,1:top_lev-1) = 0._r8
nc(1:ncol,1:top_lev-1) = 0._r8
ni(1:ncol,1:top_lev-1) = 0._r8
t1(1:ncol,1:pver) = t(1:ncol,1:pver)
q1(1:ncol,1:pver) = q(1:ncol,1:pver)
qc1(1:ncol,1:pver) = qc(1:ncol,1:pver)
qi1(1:ncol,1:pver) = qi(1:ncol,1:pver)
nc1(1:ncol,1:pver) = nc(1:ncol,1:pver)
ni1(1:ncol,1:pver) = ni(1:ncol,1:pver)

! initialize tendencies to zero
tlat1(1:ncol,1:pver)=0._r8
qvlat1(1:ncol,1:pver)=0._r8
qctend1(1:ncol,1:pver)=0._r8
qitend1(1:ncol,1:pver)=0._r8
nctend1(1:ncol,1:pver)=0._r8
nitend1(1:ncol,1:pver)=0._r8

! initialize precip output
qrout(1:ncol,1:pver)=0._r8
qsout(1:ncol,1:pver)=0._r8
nrout(1:ncol,1:pver)=0._r8
nsout(1:ncol,1:pver)=0._r8
dsout(1:ncol,1:pver)=0._r8

drout(1:ncol,1:pver)=0._r8

reff_rain(1:ncol,1:pver)=0._r8
reff_snow(1:ncol,1:pver)=0._r8

! initialize variables for trop_mozart
nevapr(1:ncol,1:pver) = 0._r8
nevapr2(1:ncol,1:pver) = 0._r8
evapsnow(1:ncol,1:pver) = 0._r8
prain(1:ncol,1:pver) = 0._r8
prodsnow(1:ncol,1:pver) = 0._r8
cmeout(1:ncol,1:pver) = 0._r8

am_evp_st(1:ncol,1:pver) = 0._r8

! for refl calc
rainrt1(1:ncol,1:pver) = 0._r8

! initialize precip fraction and output tendencies
cldmax(1:ncol,1:pver)=mincld

!initialize aerosol number
!        naer2(1:ncol,1:pver,:)=0._r8
dum2l(1:ncol,1:pver)=0._r8
dum2i(1:ncol,1:pver)=0._r8

! initialize avg precip rate
prect1(1:ncol)=0._r8
preci1(1:ncol)=0._r8

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Get humidity and saturation vapor pressures

do k=top_lev,pver

   do i=1,ncol

      ! find wet bulk temperature and saturation value for provisional t and q without
      ! condensation
      
      es(i) = svp_water(t(i,k))
      qs(i) = svp_to_qsat(es(i), p(i,k))

      ! Prevents negative values.
      if (qs(i) < 0.0_r8) then
         qs(i) = 1.0_r8
         es(i) = p(i,k)
      end if

      esl(i,k)=svp_water(t(i,k))
      esi(i,k)=svp_ice(t(i,k))

      ! hm fix, make sure when above freezing that esi=esl, not active yet
      if (t(i,k).gt.tmelt)esi(i,k)=esl(i,k)

      relhum(i,k)=q(i,k)/qs(i)

      ! get cloud fraction, check for minimum

      cldm(i,k)=max(cldn(i,k),mincld)
      cldmw(i,k)=max(cldn(i,k),mincld)

      icldm(i,k)=max(icecldf(i,k),mincld)
      lcldm(i,k)=max(liqcldf(i,k),mincld)

      ! subcolumns, set cloud fraction variables to one
      ! if cloud water or ice is present, if not present
      ! set to mincld (mincld used instead of zero, to prevent
      ! possible division by zero errors

      if (microp_uniform) then

         cldm(i,k)=mincld
         cldmw(i,k)=mincld
         icldm(i,k)=mincld
         lcldm(i,k)=mincld

         if (qc(i,k).ge.qsmall) then
            lcldm(i,k)=1._r8           
            cldm(i,k)=1._r8
            cldmw(i,k)=1._r8
         end if

         if (qi(i,k).ge.qsmall) then             
            cldm(i,k)=1._r8
            icldm(i,k)=1._r8
         end if

      end if               ! sub-columns

      ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)

      nfice(i,k)=0._r8
      dumfice=qc(i,k)+qi(i,k)
      if (dumfice.gt.qsmall .and. qi(i,k).gt.qsmall) then
         nfice(i,k)=qi(i,k)/dumfice
      endif

      if (do_cldice .and. (t(i,k).lt.tmelt - 5._r8)) then

         ! if aerosols interact with ice set number of activated ice nuclei
         dum2=naai(i,k)

         dumnnuc=(dum2-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
         dumnnuc=max(dumnnuc,0._r8)
         ! get provisional ni and qi after nucleation in order to calculate
         ! Bergeron process below
         ninew=ni(i,k)+dumnnuc*deltat
         qinew=qi(i,k)+dumnnuc*deltat*mi0

         !T>268
      else
         ninew=ni(i,k)
         qinew=qi(i,k)
      end if

      ! Initialize CME components

      cme(i,k) = 0._r8
      cmei(i,k)=0._r8


      !-------------------------------------------------------------------
      !Bergeron process

      ! make sure to initialize bergeron process to zero
      berg(i,k)=0._r8
      prd = 0._r8

      !condensation loop.

      ! get in-cloud qi and ni after nucleation
      if (icldm(i,k) .gt. 0._r8) then 
         qiic(i,k)=qinew/icldm(i,k)
         niic(i,k)=ninew/icldm(i,k)
      else
         qiic(i,k)=0._r8
         niic(i,k)=0._r8
      endif

      !if T < 0 C then bergeron.
      if (do_cldice .and. (t(i,k).lt.273.15_r8)) then

         !if ice exists
         if (qi(i,k).gt.qsmall) then

            bergtsf = 0._r8 ! bergeron time scale (fraction of timestep)

            qvi = svp_to_qsat(esi(i,k), p(i,k))
            qvl = svp_to_qsat(esl(i,k), p(i,k))

            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
            abi = 1._r8+dqsidt*xxls/cpp

            ! get ice size distribution parameters

            if (qiic(i,k).ge.qsmall) then
               lami(k) = (cons1*ci* &
                    niic(i,k)/qiic(i,k))**(1._r8/di)
               n0i(k) = niic(i,k)*lami(k)

               ! check for slope
               ! adjust vars
               if (lami(k).lt.lammini) then

                  lami(k) = lammini
                  n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
               else if (lami(k).gt.lammaxi) then
                  lami(k) = lammaxi
                  n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
               end if

               epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))

               !if liquid exists  
               if (qc(i,k).gt. qsmall) then 

                  !begin bergeron process
                  !     do bergeron (vapor deposition with RHw=1)
                  !     code to find berg (a rate) goes here

                  ! calculate Bergeron process

                  prd = epsi*(qvl-qvi)/abi

               else
                  prd = 0._r8
               end if

               ! multiply by cloud fraction

               prd = prd*min(icldm(i,k),lcldm(i,k))

               !     transfer of existing cloud liquid to ice

               berg(i,k)=max(0._r8,prd)

            end if  !end liquid exists bergeron

            if (berg(i,k).gt.0._r8) then
               bergtsf=max(0._r8,(qc(i,k)/berg(i,k))/deltat) 

               if(bergtsf.lt.1._r8) berg(i,k) = max(0._r8,qc(i,k)/deltat)

            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (bergtsf.lt.1._r8.or.icldm(i,k).gt.lcldm(i,k)) then

               if (qiic(i,k).ge.qsmall) then

                  ! first case is for case when liquid water is present, but is completely depleted 
                  ! in time step, i.e., bergrsf > 0 but < 1

                  if (qc(i,k).ge.qsmall) then
                     rhin  = (1.0_r8 + relhum(i,k)) / 2._r8
                     if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then
                        prd = epsi*(rhin*qvl-qvi)/abi

                        ! multiply by cloud fraction assuming liquid/ice maximum overlap
                        prd = prd*min(icldm(i,k),lcldm(i,k))

                        ! add to cmei
                        cmei(i,k) = cmei(i,k) + (prd * (1._r8- bergtsf))

                     end if ! rhin 
                  end if ! qc > qsmall

                  ! second case is for pure ice cloud, either no liquid, or icldm > lcldm

                  if (qc(i,k).lt.qsmall.or.icldm(i,k).gt.lcldm(i,k)) then

                     ! note: for case of no liquid, need to set liquid cloud fraction to zero
                     ! store liquid cloud fraction in 'dum'

                     if (qc(i,k).lt.qsmall) then 
                        dum=0._r8 
                     else
                        dum=lcldm(i,k)
                     end if

                     ! set RH to grid-mean value for pure ice cloud
                     rhin = relhum(i,k)

                     if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then

                        prd = epsi*(rhin*qvl-qvi)/abi

                        ! multiply by relevant cloud fraction for pure ice cloud
                        ! assuming maximum overlap of liquid/ice
                        prd = prd*max((icldm(i,k)-dum),0._r8)
                        cmei(i,k) = cmei(i,k) + prd

                     end if ! rhin
                  end if ! qc or icldm > lcldm
               end if ! qiic
            end if ! bergtsf or icldm > lcldm

            !     if deposition, it should not reduce grid mean rhi below 1.0
            if(cmei(i,k) > 0.0_r8 .and. (relhum(i,k)*esl(i,k)/esi(i,k)) > 1._r8 ) &
                 cmei(i,k)=min(cmei(i,k),(q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat)

         end if            !end ice exists loop
         !this ends temperature < 0. loop

         !-------------------------------------------------------------------
      end if  ! 
      !..............................................................

      ! evaporation should not exceed available water

      if ((-berg(i,k)).lt.-qc(i,k)/deltat) berg(i,k) = max(qc(i,k)/deltat,0._r8)

      !sublimation process...
      if (do_cldice .and. ((relhum(i,k)*esl(i,k)/esi(i,k)).lt.1._r8 .and. qiic(i,k).ge.qsmall )) then

         qvi = svp_to_qsat(esi(i,k), p(i,k))
         qvl = svp_to_qsat(esl(i,k), p(i,k))
         dqsidt =  xxls*qvi/(rv*t(i,k)**2)
         abi = 1._r8+dqsidt*xxls/cpp

         ! get ice size distribution parameters

         lami(k) = (cons1*ci* &
              niic(i,k)/qiic(i,k))**(1._r8/di)
         n0i(k) = niic(i,k)*lami(k)

         ! check for slope
         ! adjust vars
         if (lami(k).lt.lammini) then

            lami(k) = lammini
            n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
         else if (lami(k).gt.lammaxi) then
            lami(k) = lammaxi
            n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
         end if

         epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))

         ! modify for ice fraction below
         prd = epsi*(relhum(i,k)*qvl-qvi)/abi * icldm(i,k)
         cmei(i,k)=min(prd,0._r8)

      endif

      ! sublimation should not exceed available ice
      if (cmei(i,k).lt.-qi(i,k)/deltat) cmei(i,k)=-qi(i,k)/deltat

      ! sublimation should not increase grid mean rhi above 1.0 
      if(cmei(i,k) < 0.0_r8 .and. (relhum(i,k)*esl(i,k)/esi(i,k)) < 1._r8 ) &
           cmei(i,k)=min(0._r8,max(cmei(i,k),(q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat))

      ! limit cmei due for roundoff error

      cmei(i,k)=cmei(i,k)*omsm

      ! conditional for ice nucleation 
      if (do_cldice .and. (t(i,k).lt.(tmelt - 5._r8))) then 

         ! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
         ! ice nucleation rate (dum2) has already been calculated and read in (naai)

         dum2i(i,k)=naai(i,k)
      else
         dum2i(i,k)=0._r8
      end if

   end do ! i loop
end do ! k loop


!! initialize sub-step precip flux variables
do i=1,ncol
   !! flux is zero at top interface, so these should stay as 0.
   rflx1(i,1)=0._r8
   sflx1(i,1)=0._r8
   do k=top_lev,pver

      ! initialize normal and sub-step precip flux variables
      rflx1(i,k+1)=0._r8
      sflx1(i,k+1)=0._r8
   end do ! i loop
end do ! k loop
!! initialize final precip flux variables.
do i=1,ncol
   !! flux is zero at top interface, so these should stay as 0.
   rflx(i,1)=0._r8
   sflx(i,1)=0._r8
   do k=top_lev,pver
      ! initialize normal and sub-step precip flux variables
      rflx(i,k+1)=0._r8
      sflx(i,k+1)=0._r8
   end do ! i loop
end do ! k loop

do i=1,ncol
   ltrue(i)=0
   do k=top_lev,pver
      ! skip microphysical calculations if no cloud water

      if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall.or.cmei(i,k).ge.qsmall) ltrue(i)=1
   end do
end do

! assign number of sub-steps to iter
! use 2 sub-steps, following tests described in MG2008
iter = 2

! get sub-step time step
deltat=deltat/real(iter)

! since activation/nucleation processes are fast, need to take into account
! factor mtime = mixing timescale in cloud / model time step
! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk

!        note: mtime for bulk aerosols was set to: mtime=deltat/1200._r8

mtime=1._r8
rate1ord_cw2pr_st(:,:)=0._r8 ! rce 2010/05/01

!!!! skip calculations if no cloud water
do i=1,ncol
   if (ltrue(i).eq.0) then
      tlat(i,1:pver)=0._r8
      qvlat(i,1:pver)=0._r8
      qctend(i,1:pver)=0._r8
      qitend(i,1:pver)=0._r8
      qnitend(i,1:pver)=0._r8
      qrtend(i,1:pver)=0._r8
      nctend(i,1:pver)=0._r8
      nitend(i,1:pver)=0._r8
      nrtend(i,1:pver)=0._r8
      nstend(i,1:pver)=0._r8
      prect(i)=0._r8
      preci(i)=0._r8
      qniic(i,1:pver)=0._r8
      qric(i,1:pver)=0._r8
      nsic(i,1:pver)=0._r8
      nric(i,1:pver)=0._r8
      rainrt(i,1:pver)=0._r8
      goto 300
   end if

   qcsinksum_rate1ord(1:pver)=0._r8 
   qcsum_rate1ord(1:pver)=0._r8 


!!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !.....................................................................................................
   do it=1,iter

      ! initialize sub-step microphysical tendencies

      tlat(i,1:pver)=0._r8
      qvlat(i,1:pver)=0._r8
      qctend(i,1:pver)=0._r8
      qitend(i,1:pver)=0._r8
      qnitend(i,1:pver)=0._r8
      qrtend(i,1:pver)=0._r8
      nctend(i,1:pver)=0._r8
      nitend(i,1:pver)=0._r8
      nrtend(i,1:pver)=0._r8
      nstend(i,1:pver)=0._r8

      ! initialize diagnostic precipitation to zero

      qniic(i,1:pver)=0._r8
      qric(i,1:pver)=0._r8
      nsic(i,1:pver)=0._r8
      nric(i,1:pver)=0._r8

      rainrt(i,1:pver)=0._r8


      ! begin new i,k loop, calculate new cldmax after adjustment to cldm above

      ! initialize vertically-integrated rain and snow tendencies

      qrtot = 0._r8
      nrtot = 0._r8
      qstot = 0._r8
      nstot = 0._r8

      ! initialize precip at surface

      prect(i)=0._r8
      preci(i)=0._r8

      do k=top_lev,pver
      
         qcvar=relvar(i,k)
         cons2=gamma(qcvar+2.47_r8)
         cons3=gamma(qcvar)
         cons9=gamma(qcvar+2._r8)
         cons10=gamma(qcvar+1._r8)
         cons12=gamma(qcvar+1.15_r8) 
         cons15=gamma(qcvar+bc/3._r8)
         cons18=qcvar**2.47_r8
         cons19=qcvar**2
         cons20=qcvar**1.15_r8

         ! set cwml and cwmi to current qc and qi

         cwml(i,k)=qc(i,k)
         cwmi(i,k)=qi(i,k)

         ! initialize precip fallspeeds to zero

         ums(k)=0._r8 
         uns(k)=0._r8 
         umr(k)=0._r8 
         unr(k)=0._r8

         ! calculate precip fraction based on maximum overlap assumption

         ! for sub-columns cldm has already been set to 1 if cloud
         ! water or ice is present, so cldmax will be correctly set below
         ! and nothing extra needs to be done here

         if (k.eq.top_lev) then
            cldmax(i,k)=cldm(i,k)
         else
            ! if rain or snow mix ratio is smaller than
            ! threshold, then set cldmax to cloud fraction at current level

            if (do_clubb_sgs) then
               if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall) then
                  cldmax(i,k)=cldm(i,k)
               else
                  cldmax(i,k)=cldmax(i,k-1)
               end if
            else

               if (qric(i,k-1).ge.qsmall.or.qniic(i,k-1).ge.qsmall) then
                  cldmax(i,k)=max(cldmax(i,k-1),cldm(i,k))
               else
                  cldmax(i,k)=cldm(i,k)
               end if
            endif
         end if

         ! decrease in number concentration due to sublimation/evap
         ! divide by cloud fraction to get in-cloud decrease
         ! don't reduce Nc due to bergeron process

         if (cmei(i,k) < 0._r8 .and. qi(i,k) > qsmall .and. cldm(i,k) > mincld) then
            nsubi(k)=cmei(i,k)/qi(i,k)*ni(i,k)/cldm(i,k)
         else
            nsubi(k)=0._r8
         end if
         nsubc(k)=0._r8


         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%

         if (do_cldice .and. dum2i(i,k).gt.0._r8.and.t(i,k).lt.(tmelt - 5._r8).and. &
              relhum(i,k)*esl(i,k)/esi(i,k).gt. rhmini+0.05_r8) then

            !if NCAI > 0. then set numice = ncai (as before)
            !note: this is gridbox averaged

            nnuccd(k)=(dum2i(i,k)-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
            nnuccd(k)=max(nnuccd(k),0._r8)
            nimax = dum2i(i,k)*icldm(i,k)

            !Calc mass of new particles using new crystal mass...
            !also this will be multiplied by mtime as nnuccd is...

            mnuccd(k) = nnuccd(k) * mi0

            !  add mnuccd to cmei....
            cmei(i,k)= cmei(i,k) + mnuccd(k) * mtime

            !  limit cmei

            qvi = svp_to_qsat(esi(i,k), p(i,k))
            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
            abi = 1._r8+dqsidt*xxls/cpp
            cmei(i,k)=min(cmei(i,k),(q(i,k)-qvi)/abi/deltat)

            ! limit for roundoff error
            cmei(i,k)=cmei(i,k)*omsm

         else
            nnuccd(k)=0._r8
            nimax = 0._r8
            mnuccd(k) = 0._r8
         end if

         !c............................................................................
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
         ! for microphysical process calculations
         ! units are kg/kg for mixing ratio, 1/kg for number conc

         ! limit in-cloud values to 0.005 kg/kg

         qcic(i,k)=min(cwml(i,k)/lcldm(i,k),5.e-3_r8)
         qiic(i,k)=min(cwmi(i,k)/icldm(i,k),5.e-3_r8)
         ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)
         niic(i,k)=max(ni(i,k)/icldm(i,k),0._r8)

         if (qc(i,k) - berg(i,k)*deltat.lt.qsmall) then
            qcic(i,k)=0._r8
            ncic(i,k)=0._r8
            if (qc(i,k)-berg(i,k)*deltat.lt.0._r8) then
               berg(i,k)=qc(i,k)/deltat*omsm
            end if
         end if

         if (do_cldice .and. qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.qsmall) then
            qiic(i,k)=0._r8
            niic(i,k)=0._r8
            if (qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.0._r8) then
               cmei(i,k)=(-qi(i,k)/deltat-berg(i,k))*omsm
            end if
         end if

         ! add to cme output

         cmeout(i,k) = cmeout(i,k)+cmei(i,k)

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! droplet activation
         ! calculate potential for droplet activation if cloud water is present
         ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
         ! number tendency (npccnin) is read in from companion routine

         ! assume aerosols already activated are equal to number of existing droplets for simplicity
         ! multiply by cloud fraction to obtain grid-average tendency

         if (qcic(i,k).ge.qsmall) then   
            npccn(k) = max(0._r8,npccnin(i,k))  
            dum2l(i,k)=(nc(i,k)+npccn(k)*deltat)/lcldm(i,k)
            dum2l(i,k)=max(dum2l(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3  
            ncmax = dum2l(i,k)*lcldm(i,k)
         else
            npccn(k)=0._r8
            dum2l(i,k)=0._r8
            ncmax = 0._r8
         end if

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! get size distribution parameters based on in-cloud cloud water/ice 
         ! these calculations also ensure consistency between number and mixing ratio
         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         !......................................................................
         ! cloud ice

         if (qiic(i,k).ge.qsmall) then

            ! add upper limit to in-cloud number concentration to prevent numerical error
            niic(i,k)=min(niic(i,k),qiic(i,k)*1.e20_r8)

            lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**(1._r8/di)
            n0i(k) = niic(i,k)*lami(k)

            ! check for slope
            ! adjust vars

            if (lami(k).lt.lammini) then

               lami(k) = lammini
               n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
               niic(i,k) = n0i(k)/lami(k)
            else if (lami(k).gt.lammaxi) then
               lami(k) = lammaxi
               n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
               niic(i,k) = n0i(k)/lami(k)
            end if

         else
            lami(k) = 0._r8
            n0i(k) = 0._r8
         end if

         if (qcic(i,k).ge.qsmall) then

            ! add upper limit to in-cloud number concentration to prevent numerical error
            ncic(i,k)=min(ncic(i,k),qcic(i,k)*1.e20_r8)

            ncic(i,k)=max(ncic(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm  

            ! get pgam from fit to observations of martin et al. 1994

            pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
            pgam(k)=1._r8/(pgam(k)**2)-1._r8
            pgam(k)=max(pgam(k),2._r8)
            pgam(k)=min(pgam(k),15._r8)

            ! calculate lamc

            lamc(k) = (pi/6._r8*rhow*ncic(i,k)*gamma(pgam(k)+4._r8)/ &
                 (qcic(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)

            ! lammin, 50 micron diameter max mean size

            lammin = (pgam(k)+1._r8)/50.e-6_r8
            lammax = (pgam(k)+1._r8)/2.e-6_r8

            if (lamc(k).lt.lammin) then
               lamc(k) = lammin
               ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
            else if (lamc(k).gt.lammax) then
               lamc(k) = lammax
               ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
            end if

            ! parameter to calculate droplet freezing

            cdist1(k) = ncic(i,k)/gamma(pgam(k)+1._r8) 

         else
            lamc(k) = 0._r8
            cdist1(k) = 0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! begin micropysical process calculations 
         !.................................................................
         ! autoconversion of cloud liquid water to rain
         ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
         ! minimum qc of 1 x 10^-8 prevents floating point error

         if (qcic(i,k).ge.1.e-8_r8) then

            ! nprc is increase in rain number conc due to autoconversion
            ! nprc1 is decrease in cloud droplet conc due to autoconversion

            ! assume exponential sub-grid distribution of qc, resulting in additional
            ! factor related to qcvar below

            ! hm switch for sub-columns, don't include sub-grid qc
            if (microp_uniform) then

               prc(k) = 1350._r8*qcic(i,k)**2.47_r8* &
                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
               nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
               nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))

            else

               prc(k) = cons2/(cons3*cons18)*1350._r8*qcic(i,k)**2.47_r8* &
                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
               nprc(k) = prc(k)/cons22
               nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))

            end if               ! sub-column switch

         else
            prc(k)=0._r8
            nprc(k)=0._r8
            nprc1(k)=0._r8
         end if

         ! add autoconversion to precip from above to get provisional rain mixing ratio
         ! and number concentration (qric and nric)

         ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)

         dum=0.45_r8
         dum1=0.45_r8

         if (k.eq.top_lev) then
            qric(i,k)=prc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
            nric(i,k)=nprc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
         else
            if (qric(i,k-1).ge.qsmall) then
               dum=umr(k-1)
               dum1=unr(k-1)
            end if

            ! no autoconversion of rain number if rain/snow falling from above
            ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
            ! by the existing rain/snow particles from above

            if (qric(i,k-1).ge.1.e-9_r8.or.qniic(i,k-1).ge.1.e-9_r8) then
               nprc(k)=0._r8
            end if

            qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                 (rho(i,k)*dz(i,k)*((pra(k-1)+prc(k))*lcldm(i,k)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax(i,k))))&
                 /(dum*rho(i,k)*cldmax(i,k))
            nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                 (rho(i,k)*dz(i,k)*(nprc(k)*lcldm(i,k)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax(i,k))))&
                 /(dum1*rho(i,k)*cldmax(i,k))

         end if

         !.......................................................................
         ! Autoconversion of cloud ice to snow
         ! similar to Ferrier (1994)

         if (do_cldice) then
            if (t(i,k).le.273.15_r8.and.qiic(i,k).ge.qsmall) then

               ! note: assumes autoconversion timescale of 180 sec
               
               nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)

               prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
                    (cons23/lami(k)+3._r8*cons24/lami(k)**2+ &
                    6._r8*dcs/lami(k)**3+6._r8/lami(k)**4)*exp(-lami(k)*dcs)
            else
               prci(k)=0._r8
               nprci(k)=0._r8
            end if
         else
            ! Add in the particles that we have already converted to snow, and
            ! don't do any further autoconversion of ice.
            prci(k)  = tnd_qsnow(i, k) / cldm(i,k)
            nprci(k) = tnd_nsnow(i, k) / cldm(i,k)
         end if

         ! add autoconversion to flux from level above to get provisional snow mixing ratio
         ! and number concentration (qniic and nsic)

         dum=(asn(i,k)*cons25)
         dum1=(asn(i,k)*cons25)

         if (k.eq.top_lev) then
            qniic(i,k)=prci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
            nsic(i,k)=nprci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
         else
            if (qniic(i,k-1).ge.qsmall) then
               dum=ums(k-1)
               dum1=uns(k-1)
            end if

            qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
                 (rho(i,k)*dz(i,k)*((prci(k)+prai(k-1)+psacws(k-1)+bergs(k-1))*icldm(i,k)+(prds(k-1)+ &
                 pracs(k-1)+mnuccr(k-1))*cldmax(i,k))))&
                 /(dum*rho(i,k)*cldmax(i,k))

            nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                 (rho(i,k)*dz(i,k)*(nprci(k)*icldm(i,k)+(nsubs(k-1)+nsagg(k-1)+nnuccr(k-1))*cldmax(i,k))))&
                 /(dum1*rho(i,k)*cldmax(i,k))

         end if

         ! if precip mix ratio is zero so should number concentration

         if (qniic(i,k).lt.qsmall) then
            qniic(i,k)=0._r8
            nsic(i,k)=0._r8
         end if

         if (qric(i,k).lt.qsmall) then
            qric(i,k)=0._r8
            nric(i,k)=0._r8
         end if

         ! make sure number concentration is a positive number to avoid 
         ! taking root of negative later

         nric(i,k)=max(nric(i,k),0._r8)
         nsic(i,k)=max(nsic(i,k),0._r8)

         !.......................................................................
         ! get size distribution parameters for precip
         !......................................................................
         ! rain

         if (qric(i,k).ge.qsmall) then
            lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
            n0r(k) = nric(i,k)*lamr(k)

            ! check for slope
            ! adjust vars

            if (lamr(k).lt.lamminr) then

               lamr(k) = lamminr

               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
            else if (lamr(k).gt.lammaxr) then
               lamr(k) = lammaxr
               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
            end if

            ! provisional rain number and mass weighted mean fallspeed (m/s)

            unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
            umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(i,k))

         else
            lamr(k) = 0._r8
            n0r(k) = 0._r8
            umr(k) = 0._r8
            unr(k) = 0._r8
         end if

         !......................................................................
         ! snow

         if (qniic(i,k).ge.qsmall) then
            lams(k) = (cons6*cs*nsic(i,k)/qniic(i,k))**(1._r8/ds)
            n0s(k) = nsic(i,k)*lams(k)

            ! check for slope
            ! adjust vars

            if (lams(k).lt.lammins) then
               lams(k) = lammins
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
               nsic(i,k) = n0s(k)/lams(k)

            else if (lams(k).gt.lammaxs) then
               lams(k) = lammaxs
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
               nsic(i,k) = n0s(k)/lams(k)
            end if

            ! provisional snow number and mass weighted mean fallspeed (m/s)

            ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(i,k))
            uns(k) = min(asn(i,k)*cons7/lams(k)**bs,1.2_r8*rhof(i,k))

         else
            lams(k) = 0._r8
            n0s(k) = 0._r8
            ums(k) = 0._r8
            uns(k) = 0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         ! heterogeneous freezing of cloud water

         if (.not. use_hetfrz_classnuc) then

            if (do_cldice .and. qcic(i,k).ge.qsmall .and. t(i,k).lt.269.15_r8) then

               ! immersion freezing (Bigg, 1953)


               ! subcolumns

               if (microp_uniform) then

                  mnuccc(k) = &
                     pi*pi/36._r8*rhow* &
                     cdist1(k)*gamma(7._r8+pgam(k))* &
                     bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
                     lamc(k)**3/lamc(k)**3

                  nnuccc(k) = &
                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
                     *bimm* &
                     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3

               else

                  mnuccc(k) = cons9/(cons3*cons19)* &
                     pi*pi/36._r8*rhow* &
                     cdist1(k)*gamma(7._r8+pgam(k))* &
                     bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
                     lamc(k)**3/lamc(k)**3

                  nnuccc(k) = cons10/(cons3*qcvar)* &
                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
                     *bimm* &
                     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3
               end if           ! sub-columns


               ! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
               ! dust size and number in 4 bins are read in from companion routine

               tcnt=(270.16_r8-t(i,k))**1.3_r8
               viscosity=1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
               mfp=2.0_r8*viscosity/(p(i,k)  &                   ! Mean free path (m)
                  *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k))))           

               nslip1=1.0_r8+(mfp/rndst(i,k,1))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,1)/mfp))))! Slip correction factor
               nslip2=1.0_r8+(mfp/rndst(i,k,2))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,2)/mfp))))
               nslip3=1.0_r8+(mfp/rndst(i,k,3))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,3)/mfp))))
               nslip4=1.0_r8+(mfp/rndst(i,k,4))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,4)/mfp))))

               ndfaer1=1.381e-23_r8*t(i,k)*nslip1/(6._r8*pi*viscosity*rndst(i,k,1))  ! aerosol diffusivity (m2/s)
               ndfaer2=1.381e-23_r8*t(i,k)*nslip2/(6._r8*pi*viscosity*rndst(i,k,2))
               ndfaer3=1.381e-23_r8*t(i,k)*nslip3/(6._r8*pi*viscosity*rndst(i,k,3))
               ndfaer4=1.381e-23_r8*t(i,k)*nslip4/(6._r8*pi*viscosity*rndst(i,k,4))


               if (microp_uniform) then

                  mnucct(k) = &
                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                     cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4

                  nnucct(k) = (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)

               else

                  mnucct(k) = gamma(qcvar+4._r8/3._r8)/(cons3*qcvar**(4._r8/3._r8))*  &
                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                     cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4

                  nnucct(k) =  gamma(qcvar+1._r8/3._r8)/(cons3*qcvar**(1._r8/3._r8))*  &
                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)

               end if      ! sub-column switch

               ! make sure number of droplets frozen does not exceed available ice nuclei concentration
               ! this prevents 'runaway' droplet freezing

               if (nnuccc(k)*lcldm(i,k).gt.nnuccd(k)) then
                  dum=(nnuccd(k)/(nnuccc(k)*lcldm(i,k)))
                  ! scale mixing ratio of droplet freezing with limit
                  mnuccc(k)=mnuccc(k)*dum
                  nnuccc(k)=nnuccd(k)/lcldm(i,k)
               end if

            else
               mnuccc(k)=0._r8
               nnuccc(k)=0._r8
               mnucct(k)=0._r8
               nnucct(k)=0._r8
            end if

         else
            if (do_cldice .and. qcic(i,k) >= qsmall) then
               con1 = 1._r8/(1.333_r8*pi)**0.333_r8
               r3lx = con1*(rho(i,k)*qcic(i,k)/(rhow*max(ncic(i,k)*rho(i,k), 1.0e6_r8)))**0.333_r8 ! in m
               r3lx = max(4.e-6_r8, r3lx)
               mi0l = 4._r8/3._r8*pi*rhow*r3lx**3_r8
                
               nnuccc(k) = frzimm(i,k)*1.0e6_r8/rho(i,k)
               mnuccc(k) = nnuccc(k)*mi0l 

               nnucct(k) = frzcnt(i,k)*1.0e6_r8/rho(i,k)
               mnucct(k) = nnucct(k)*mi0l 

               nnudep(k) = frzdep(i,k)*1.0e6_r8/rho(i,k)
               mnudep(k) = nnudep(k)*mi0
            else
               nnuccc(k) = 0._r8
               mnuccc(k) = 0._r8

               nnucct(k) = 0._r8
               mnucct(k) = 0._r8

               nnudep(k) = 0._r8
               mnudep(k) = 0._r8
            end if
         endif


         !.......................................................................
         ! snow self-aggregation from passarelli, 1978, used by reisner, 1998
         ! this is hard-wired for bs = 0.4 for now
         ! ignore self-collection of cloud ice

         if (qniic(i,k).ge.qsmall .and. t(i,k).le.273.15_r8) then
            nsagg(k) = -1108._r8*asn(i,k)*Eii* &
                 pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)*rho(i,k)** &
                 ((2._r8+bs)/3._r8)*qniic(i,k)**((2._r8+bs)/3._r8)* &
                 (nsic(i,k)*rho(i,k))**((4._r8-bs)/3._r8)/ &
                 (4._r8*720._r8*rho(i,k))
         else
            nsagg(k)=0._r8
         end if

         !.......................................................................
         ! accretion of cloud droplets onto snow/graupel
         ! here use continuous collection equation with
         ! simple gravitational collection kernel
         ! ignore collisions between droplets/cloud ice
         ! since minimum size ice particle for accretion is 50 - 150 micron

         ! ignore collision of snow with droplets above freezing

         if (qniic(i,k).ge.qsmall .and. t(i,k).le.tmelt .and. &
              qcic(i,k).ge.qsmall) then

            ! put in size dependent collection efficiency
            ! mean diameter of snow is area-weighted, since
            ! accretion is function of crystal geometric area
            ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

            dc0 = (pgam(k)+1._r8)/lamc(k)
            ds0 = 1._r8/lams(k)
            dum = dc0*dc0*uns(k)*rhow/(9._r8*mu(i,k)*ds0)
            eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

            eci = max(eci,0._r8)
            eci = min(eci,1._r8)


            ! no impact of sub-grid distribution of qc since psacws
            ! is linear in qc

            psacws(k) = pi/4._r8*asn(i,k)*qcic(i,k)*rho(i,k)* &
                 n0s(k)*Eci*cons11/ &
                 lams(k)**(bs+3._r8)
            npsacws(k) = pi/4._r8*asn(i,k)*ncic(i,k)*rho(i,k)* &
                 n0s(k)*Eci*cons11/ &
                 lams(k)**(bs+3._r8)
         else
            psacws(k)=0._r8
            npsacws(k)=0._r8
         end if

         ! add secondary ice production due to accretion of droplets by snow 
         ! (Hallet-Mossop process) (from Cotton et al., 1986)

         if (.not. do_cldice) then
            ni_secp   = 0.0_r8
            nsacwi(k) = 0.0_r8
            msacwi(k) = 0.0_r8
         else if((t(i,k).lt.270.16_r8) .and. (t(i,k).ge.268.16_r8)) then
            ni_secp   = 3.5e8_r8*(270.16_r8-t(i,k))/2.0_r8*psacws(k)
            nsacwi(k) = ni_secp
            msacwi(k) = min(ni_secp*mi0,psacws(k))
         else if((t(i,k).lt.268.16_r8) .and. (t(i,k).ge.265.16_r8)) then
            ni_secp   = 3.5e8_r8*(t(i,k)-265.16_r8)/3.0_r8*psacws(k)
            nsacwi(k) = ni_secp
            msacwi(k) = min(ni_secp*mi0,psacws(k))
         else
            ni_secp   = 0.0_r8
            nsacwi(k) = 0.0_r8
            msacwi(k) = 0.0_r8
         endif
         psacws(k) = max(0.0_r8,psacws(k)-ni_secp*mi0)

         !.......................................................................
         ! accretion of rain water by snow
         ! formula from ikawa and saito, 1991, used by reisner et al., 1998

         if (qric(i,k).ge.1.e-8_r8 .and. qniic(i,k).ge.1.e-8_r8 .and. & 
              t(i,k).le.273.15_r8) then

            pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))**2+ &
                 0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(i,k)* &
                 n0r(k)*n0s(k)* &
                 (5._r8/(lamr(k)**6*lams(k))+ &
                 2._r8/(lamr(k)**5*lams(k)**2)+ &
                 0.5_r8/(lamr(k)**4*lams(k)**3)))

            npracs(k) = pi/2._r8*rho(i,k)*ecr*(1.7_r8*(unr(k)-uns(k))**2+ &
                 0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
                 (1._r8/(lamr(k)**3*lams(k))+ &
                 1._r8/(lamr(k)**2*lams(k)**2)+ &
                 1._r8/(lamr(k)*lams(k)**3))

         else
            pracs(k)=0._r8
            npracs(k)=0._r8
         end if

         !.......................................................................
         ! heterogeneous freezing of rain drops
         ! follows from Bigg (1953)

         if (t(i,k).lt.269.15_r8 .and. qric(i,k).ge.qsmall) then

            mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
                 (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamr(k)**3 &
                 /lamr(k)**3

            nnuccr(k) = pi*nric(i,k)*bimm* &
                 (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamr(k)**3
         else
            mnuccr(k)=0._r8
            nnuccr(k)=0._r8
         end if

         !.......................................................................
         ! accretion of cloud liquid water by rain
         ! formula from Khrouditnov and Kogan (2000)
         ! gravitational collection kernel, droplet fall speed neglected

         if (qric(i,k).ge.qsmall .and. qcic(i,k).ge.qsmall) then

            ! include sub-grid distribution of cloud water

            ! add sub-column switch

            if (microp_uniform) then

               pra(k) = 67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
               npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))

            else

               pra(k) = accre_enhan(i,k)*(cons12/(cons3*cons20)*67._r8*(qcic(i,k)*qric(i,k))**1.15_r8)
               npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))

            end if               ! sub-column switch

         else
            pra(k)=0._r8
            npra(k)=0._r8
         end if

         !.......................................................................
         ! Self-collection of rain drops
         ! from Beheng(1994)

         if (qric(i,k).ge.qsmall) then
            nragg(k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
         else
            nragg(k)=0._r8
         end if

         !.......................................................................
         ! Accretion of cloud ice by snow
         ! For this calculation, it is assumed that the Vs >> Vi
         ! and Ds >> Di for continuous collection

         if (do_cldice .and. qniic(i,k).ge.qsmall.and.qiic(i,k).ge.qsmall &
              .and.t(i,k).le.273.15_r8) then

            prai(k) = pi/4._r8*asn(i,k)*qiic(i,k)*rho(i,k)* &
                 n0s(k)*Eii*cons11/ &
                 lams(k)**(bs+3._r8)
            nprai(k) = pi/4._r8*asn(i,k)*niic(i,k)* &
                 rho(i,k)*n0s(k)*Eii*cons11/ &
                 lams(k)**(bs+3._r8)
         else
            prai(k)=0._r8
            nprai(k)=0._r8
         end if

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! calculate evaporation/sublimation of rain and snow
         ! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
         ! in-cloud condensation/deposition of rain and snow is neglected
         ! except for transfer of cloud water to snow through bergeron process

         ! initialize evap/sub tendncies
         pre(k)=0._r8
         prds(k)=0._r8

         ! evaporation of rain
         ! only calculate if there is some precip fraction > cloud fraction

         if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8.or.cldmax(i,k).gt.lcldm(i,k)) then

            ! set temporary cloud fraction to zero if cloud water + ice is very small
            ! this will ensure that evaporation/sublimation of precip occurs over
            ! entire grid cell, since min cloud fraction is specified otherwise
            if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8) then
               dum=0._r8
            else
               dum=lcldm(i,k)
            end if

            ! saturation vapor pressure
            esn=svp_water(t(i,k))
            qsn=svp_to_qsat(esn, p(i,k))

            ! recalculate saturation vapor pressure for liquid and ice
            esl(i,k)=esn
            esi(i,k)=svp_ice(t(i,k))
            ! hm fix, make sure when above freezing that esi=esl, not active yet
            if (t(i,k).gt.tmelt)esi(i,k)=esl(i,k)

            ! calculate q for out-of-cloud region
            qclr=(q(i,k)-dum*qsn)/(1._r8-dum)

            if (qric(i,k).ge.qsmall) then

               qvs=svp_to_qsat(esl(i,k), p(i,k))
               dqsdt = xxlv*qvs/(rv*t(i,k)**2)
               ab = 1._r8+dqsdt*xxlv/cpp
               epsr = 2._r8*pi*n0r(k)*rho(i,k)*Dv(i,k)* &
                    (f1r/(lamr(k)*lamr(k))+ &
                    f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                    sc(i,k)**(1._r8/3._r8)*cons13/ &
                    (lamr(k)**(5._r8/2._r8+br/2._r8)))

               pre(k) = epsr*(qclr-qvs)/ab

               ! only evaporate in out-of-cloud region
               ! and distribute across cldmax
               pre(k)=min(pre(k)*(cldmax(i,k)-dum),0._r8)
               pre(k)=pre(k)/cldmax(i,k)
               am_evp_st(i,k) = max(cldmax(i,k)-dum, 0._r8)
            end if

            ! sublimation of snow
            if (qniic(i,k).ge.qsmall) then
               qvi=svp_to_qsat(esi(i,k), p(i,k))
               dqsidt =  xxls*qvi/(rv*t(i,k)**2)
               abi = 1._r8+dqsidt*xxls/cpp
               epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                    (f1s/(lams(k)*lams(k))+ &
                    f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                    sc(i,k)**(1._r8/3._r8)*cons14/ &
                    (lams(k)**(5._r8/2._r8+bs/2._r8)))
               prds(k) = epss*(qclr-qvi)/abi

               ! only sublimate in out-of-cloud region and distribute over cldmax
               prds(k)=min(prds(k)*(cldmax(i,k)-dum),0._r8)
               prds(k)=prds(k)/cldmax(i,k)
               am_evp_st(i,k) = max(cldmax(i,k)-dum, 0._r8)
            end if

            ! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
            ! get updated RH at end of time step based on cloud water/ice condensation/evap

            qtmp=q(i,k)-(cmei(i,k)+(pre(k)+prds(k))*cldmax(i,k))*deltat
            ttmp=t(i,k)+((pre(k)*cldmax(i,k))*xxlv+ &
                 (cmei(i,k)+prds(k)*cldmax(i,k))*xxls)*deltat/cpp

            !limit range of temperatures!
            ttmp=max(180._r8,min(ttmp,323._r8))

            esn=svp_water(ttmp)  ! use rhw to allow ice supersaturation
            qsn=svp_to_qsat(esn, p(i,k))

            ! modify precip evaporation rate if q > qsat
            if (qtmp.gt.qsn) then
               if (pre(k)+prds(k).lt.-1.e-20_r8) then
                  dum1=pre(k)/(pre(k)+prds(k))
                  ! recalculate q and t after cloud water cond but without precip evap
                  qtmp=q(i,k)-(cmei(i,k))*deltat
                  ttmp=t(i,k)+(cmei(i,k)*xxls)*deltat/cpp
                  esn=svp_water(ttmp) ! use rhw to allow ice supersaturation
                  qsn=svp_to_qsat(esn, p(i,k))
                  dum=(qtmp-qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))
                  dum=min(dum,0._r8)

                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
                  pre(k)=dum*dum1/deltat/cldmax(i,k)

                  ! do separately using RHI for prds....
                  esn=svp_ice(ttmp) ! use rhi to allow ice supersaturation
                  qsn=svp_to_qsat(esn, p(i,k))
                  dum=(qtmp-qsn)/(1._r8 + cons28*qsn/(cpp*rv*ttmp**2))
                  dum=min(dum,0._r8)

                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
                  prds(k)=dum*(1._r8-dum1)/deltat/cldmax(i,k)
               end if
            end if
         end if

         ! bergeron process - evaporation of droplets and deposition onto snow

         if (qniic(i,k).ge.qsmall.and.qcic(i,k).ge.qsmall.and.t(i,k).lt.tmelt) then
            qvi=svp_to_qsat(esi(i,k), p(i,k))
            qvs=svp_to_qsat(esl(i,k), p(i,k))
            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
            abi = 1._r8+dqsidt*xxls/cpp
            epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                 (f1s/(lams(k)*lams(k))+ &
                 f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                 sc(i,k)**(1._r8/3._r8)*cons14/ &
                 (lams(k)**(5._r8/2._r8+bs/2._r8)))
            bergs(k)=epss*(qvs-qvi)/abi
         else
            bergs(k)=0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! conservation to ensure no negative values of cloud water/precipitation
         ! in case microphysical process rates are large

         ! make sure and use end-of-time step values for cloud water, ice, due
         ! condensation/deposition

         ! note: for check on conservation, processes are multiplied by omsm
         ! to prevent problems due to round off error

         ! include mixing timescale  (mtime)

         qce=(qc(i,k) - berg(i,k)*deltat)
         nce=(nc(i,k)+npccn(k)*deltat*mtime)
         qie=(qi(i,k)+(cmei(i,k)+berg(i,k))*deltat)
         nie=(ni(i,k)+nnuccd(k)*deltat*mtime)

         ! conservation of qc

         dum = (prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+ &
              psacws(k)+bergs(k))*lcldm(i,k)*deltat

         if (dum.gt.qce) then
            ratio = qce/deltat/lcldm(i,k)/(prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+psacws(k)+bergs(k))*omsm 

            prc(k) = prc(k)*ratio
            pra(k) = pra(k)*ratio
            mnuccc(k) = mnuccc(k)*ratio
            mnucct(k) = mnucct(k)*ratio  
            msacwi(k) = msacwi(k)*ratio  
            psacws(k) = psacws(k)*ratio
            bergs(k) = bergs(k)*ratio
         end if

         ! conservation of nc

         dum = (nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
              npsacws(k)-nsubc(k))*lcldm(i,k)*deltat

         if (dum.gt.nce) then
            ratio = nce/deltat/((nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+&
                 npsacws(k)-nsubc(k))*lcldm(i,k))*omsm

            nprc1(k) = nprc1(k)*ratio
            npra(k) = npra(k)*ratio
            nnuccc(k) = nnuccc(k)*ratio
            nnucct(k) = nnucct(k)*ratio  
            npsacws(k) = npsacws(k)*ratio
            nsubc(k)=nsubc(k)*ratio
         end if

         ! conservation of qi

         if (do_cldice) then

            frztmp = -mnuccc(k) - mnucct(k) - msacwi(k)
            if (use_hetfrz_classnuc) frztmp = -mnuccc(k)-mnucct(k)-mnudep(k)-msacwi(k)
            dum = ( frztmp*lcldm(i,k) + (prci(k)+prai(k))*icldm(i,k) )*deltat

            if (dum.gt.qie) then

               frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
               if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
               ratio = (qie/deltat + frztmp*lcldm(i,k))/((prci(k)+prai(k))*icldm(i,k))*omsm 
               prci(k) = prci(k)*ratio
               prai(k) = prai(k)*ratio
            end if

            ! conservation of ni
            frztmp = -nnucct(k) - nsacwi(k)
            if (use_hetfrz_classnuc) frztmp = -nnucct(k) - nnuccc(k) - nnudep(k) - nsacwi(k)
            dum = ( frztmp*lcldm(i,k) + (nprci(k)+nprai(k)-nsubi(k))*icldm(i,k) )*deltat

            if (dum.gt.nie) then

               frztmp = nnucct(k) + nsacwi(k)
               if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
               ratio = (nie/deltat + frztmp*lcldm(i,k))/ &  
                     ((nprci(k)+nprai(k)-nsubi(k))*icldm(i,k))*omsm
               nprci(k) = nprci(k)*ratio
               nprai(k) = nprai(k)*ratio
               nsubi(k) = nsubi(k)*ratio
            end if
         end if

         ! for precipitation conservation, use logic that vertical integral 
         ! of tendency from current level to top of model (i.e., qrtot) cannot be negative

         ! conservation of rain mixing rat

         if (((prc(k)+pra(k))*lcldm(i,k)+(-mnuccr(k)+pre(k)-pracs(k))*&
              cldmax(i,k))*dz(i,k)*rho(i,k)+qrtot.lt.0._r8) then

            if (-pre(k)+pracs(k)+mnuccr(k).ge.qsmall) then

               ratio = (qrtot/(dz(i,k)*rho(i,k))+(prc(k)+pra(k))*lcldm(i,k))/&
                    ((-pre(k)+pracs(k)+mnuccr(k))*cldmax(i,k))*omsm 

               pre(k) = pre(k)*ratio
               pracs(k) = pracs(k)*ratio
               mnuccr(k) = mnuccr(k)*ratio
            end if
         end if

         ! conservation of nr
         ! for now neglect evaporation of nr
         nsubr(k)=0._r8

         if ((nprc(k)*lcldm(i,k)+(-nnuccr(k)+nsubr(k)-npracs(k)&
              +nragg(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+nrtot.lt.0._r8) then

            if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).ge.qsmall) then

               ratio = (nrtot/(dz(i,k)*rho(i,k))+nprc(k)*lcldm(i,k))/&
                    ((-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k))*cldmax(i,k))*omsm

               nsubr(k) = nsubr(k)*ratio
               npracs(k) = npracs(k)*ratio
               nnuccr(k) = nnuccr(k)*ratio
               nragg(k) = nragg(k)*ratio
            end if
         end if

         ! conservation of snow mix ratio

         if (((bergs(k)+psacws(k))*lcldm(i,k)+(prai(k)+prci(k))*icldm(i,k)+(pracs(k)+&
              mnuccr(k)+prds(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+qstot.lt.0._r8) then

            if (-prds(k).ge.qsmall) then

               ratio = (qstot/(dz(i,k)*rho(i,k))+(bergs(k)+psacws(k))*lcldm(i,k)+(prai(k)+prci(k))*icldm(i,k)+&
                    (pracs(k)+mnuccr(k))*cldmax(i,k))/(-prds(k)*cldmax(i,k))*omsm

               prds(k) = prds(k)*ratio
            end if
         end if

         ! conservation of ns

         ! calculate loss of number due to sublimation
         ! for now neglect sublimation of ns
         nsubs(k)=0._r8

         if ((nprci(k)*icldm(i,k)+(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax(i,k))*&
              dz(i,k)*rho(i,k)+nstot.lt.0._r8) then

            if (-nsubs(k)-nsagg(k).ge.qsmall) then

               ratio = (nstot/(dz(i,k)*rho(i,k))+nprci(k)*icldm(i,k)+&
                    nnuccr(k)*cldmax(i,k))/((-nsubs(k)-nsagg(k))*cldmax(i,k))*omsm

               nsubs(k) = nsubs(k)*ratio
               nsagg(k) = nsagg(k)*ratio
            end if
         end if

         ! get tendencies due to microphysical conversion processes
         ! note: tendencies are multiplied by appropaiate cloud/precip 
         ! fraction to get grid-scale values
         ! note: cmei is already grid-average values

         qvlat(i,k) = qvlat(i,k)-(pre(k)+prds(k))*cldmax(i,k)-cmei(i,k) 

         tlat(i,k) = tlat(i,k)+((pre(k)*cldmax(i,k)) &
              *xxlv+(prds(k)*cldmax(i,k)+cmei(i,k))*xxls+ &
              ((bergs(k)+psacws(k)+mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(i,k)+(mnuccr(k)+ &
              pracs(k))*cldmax(i,k)+berg(i,k))*xlf)

         qctend(i,k) = qctend(i,k)+ &
              (-pra(k)-prc(k)-mnuccc(k)-mnucct(k)-msacwi(k)- & 
              psacws(k)-bergs(k))*lcldm(i,k)-berg(i,k)

         if (do_cldice) then

            frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
            if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
            qitend(i,k) = qitend(i,k) + frztmp*lcldm(i,k) + &
               (-prci(k)-prai(k))*icldm(i,k) + cmei(i,k) + berg(i,k)

         end if

         qrtend(i,k) = qrtend(i,k)+ &
              (pra(k)+prc(k))*lcldm(i,k)+(pre(k)-pracs(k)- &
              mnuccr(k))*cldmax(i,k)

         qnitend(i,k) = qnitend(i,k)+ &
              (prai(k)+prci(k))*icldm(i,k)+(psacws(k)+bergs(k))*lcldm(i,k)+(prds(k)+ &
              pracs(k)+mnuccr(k))*cldmax(i,k)

         ! add output for cmei (accumulate)
         cmeiout(i,k) = cmeiout(i,k) + cmei(i,k)

         ! assign variables for trop_mozart, these are grid-average
         ! evaporation/sublimation is stored here as positive term

         evapsnow(i,k) = evapsnow(i,k)-prds(k)*cldmax(i,k)
         nevapr(i,k) = nevapr(i,k)-pre(k)*cldmax(i,k)
         nevapr2(i,k) = nevapr2(i,k)-pre(k)*cldmax(i,k)

         ! change to make sure prain is positive: do not remove snow from
         ! prain used for wet deposition
         prain(i,k) = prain(i,k)+(pra(k)+prc(k))*lcldm(i,k)+(-pracs(k)- &
              mnuccr(k))*cldmax(i,k)
         prodsnow(i,k) = prodsnow(i,k)+(prai(k)+prci(k))*icldm(i,k)+(psacws(k)+bergs(k))*lcldm(i,k)+(&
              pracs(k)+mnuccr(k))*cldmax(i,k)

         ! following are used to calculate 1st order conversion rate of cloud water
         !    to rain and snow (1/s), for later use in aerosol wet removal routine
         ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
         !    used to calculate pra, prc, ... in this routine
         ! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
         !                      (no cloud ice or bergeron terms)
         ! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }

         qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k) + (pra(k)+prc(k)+psacws(k))*lcldm(i,k) 
         qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc(i,k) 

         ! microphysics output, note this is grid-averaged
         prao(i,k)=prao(i,k)+pra(k)*lcldm(i,k)
         prco(i,k)=prco(i,k)+prc(k)*lcldm(i,k)
         mnuccco(i,k)=mnuccco(i,k)+mnuccc(k)*lcldm(i,k)
         mnuccto(i,k)=mnuccto(i,k)+mnucct(k)*lcldm(i,k)
         mnuccdo(i,k)=mnuccdo(i,k)+mnuccd(k)*lcldm(i,k)
         msacwio(i,k)=msacwio(i,k)+msacwi(k)*lcldm(i,k)
         psacwso(i,k)=psacwso(i,k)+psacws(k)*lcldm(i,k)
         bergso(i,k)=bergso(i,k)+bergs(k)*lcldm(i,k)
         bergo(i,k)=bergo(i,k)+berg(i,k)
         prcio(i,k)=prcio(i,k)+prci(k)*icldm(i,k)
         praio(i,k)=praio(i,k)+prai(k)*icldm(i,k)
         mnuccro(i,k)=mnuccro(i,k)+mnuccr(k)*cldmax(i,k)
         pracso (i,k)=pracso (i,k)+pracs (k)*cldmax(i,k)

         ! multiply activation/nucleation by mtime to account for fast timescale

         nctend(i,k) = nctend(i,k)+ npccn(k)*mtime+&
              (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k) & 
              -npra(k)-nprc1(k))*lcldm(i,k)      

         if (do_cldice) then

            frztmp = nnucct(k) + nsacwi(k)
            if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
            nitend(i,k) = nitend(i,k) + nnuccd(k)*mtime + & 
                  frztmp*lcldm(i,k) + (nsubi(k)-nprci(k)-nprai(k))*icldm(i,k)

         end if

         nstend(i,k) = nstend(i,k)+(nsubs(k)+ &
              nsagg(k)+nnuccr(k))*cldmax(i,k)+nprci(k)*icldm(i,k)

         nrtend(i,k) = nrtend(i,k)+ &
              nprc(k)*lcldm(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
              +nragg(k))*cldmax(i,k)

         ! make sure that nc and ni at advanced time step do not exceed
         ! maximum (existing N + source terms*dt), which is possible due to
         ! fast nucleation timescale

         if (nctend(i,k).gt.0._r8.and.nc(i,k)+nctend(i,k)*deltat.gt.ncmax) then
            nctend(i,k)=max(0._r8,(ncmax-nc(i,k))/deltat)
         end if

         if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax) then
            nitend(i,k)=max(0._r8,(nimax-ni(i,k))/deltat)
         end if

         ! get final values for precipitation q and N, based on
         ! flux of precip from above, source/sink term, and terminal fallspeed
         ! see eq. 15-16 in MG2008

         ! rain

         if (qric(i,k).ge.qsmall) then
            if (k.eq.top_lev) then
               qric(i,k)=qrtend(i,k)*dz(i,k)/cldmax(i,k)/umr(k)
               nric(i,k)=nrtend(i,k)*dz(i,k)/cldmax(i,k)/unr(k)
            else
               qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                    (rho(i,k)*dz(i,k)*qrtend(i,k)))/(umr(k)*rho(i,k)*cldmax(i,k))
               nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                    (rho(i,k)*dz(i,k)*nrtend(i,k)))/(unr(k)*rho(i,k)*cldmax(i,k))

            end if
         else
            qric(i,k)=0._r8
            nric(i,k)=0._r8
         end if

         ! snow

         if (qniic(i,k).ge.qsmall) then
            if (k.eq.top_lev) then
               qniic(i,k)=qnitend(i,k)*dz(i,k)/cldmax(i,k)/ums(k)
               nsic(i,k)=nstend(i,k)*dz(i,k)/cldmax(i,k)/uns(k)
            else
               qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
                    (rho(i,k)*dz(i,k)*qnitend(i,k)))/(ums(k)*rho(i,k)*cldmax(i,k))
               nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                    (rho(i,k)*dz(i,k)*nstend(i,k)))/(uns(k)*rho(i,k)*cldmax(i,k))
            end if
         else
            qniic(i,k)=0._r8
            nsic(i,k)=0._r8
         end if

         ! calculate precipitation flux at surface
         ! divide by density of water to get units of m/s

         prect(i) = prect(i)+(qrtend(i,k)*dz(i,k)*rho(i,k)+&
              qnitend(i,k)*dz(i,k)*rho(i,k))/rhow
         preci(i) = preci(i)+qnitend(i,k)*dz(i,k)*rho(i,k)/rhow

         ! convert rain rate from m/s to mm/hr

         rainrt(i,k)=qric(i,k)*rho(i,k)*umr(k)/rhow*3600._r8*1000._r8

         ! vertically-integrated precip source/sink terms (note: grid-averaged)

         qrtot = max(qrtot+qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
         qstot = max(qstot+qnitend(i,k)*dz(i,k)*rho(i,k),0._r8)
         nrtot = max(nrtot+nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
         nstot = max(nstot+nstend(i,k)*dz(i,k)*rho(i,k),0._r8)

         ! calculate melting and freezing of precip

         ! melt snow at +2 C

         if (t(i,k)+tlat(i,k)/cpp*deltat > 275.15_r8) then
            if (qstot > 0._r8) then

               ! make sure melting snow doesn't reduce temperature below threshold
               dum = -xlf/cpp*qstot/(dz(i,k)*rho(i,k))
               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.275.15_r8) then
                  dum = (t(i,k)+tlat(i,k)/cpp*deltat-275.15_r8)*cpp/xlf
                  dum = dum/(xlf/cpp*qstot/(dz(i,k)*rho(i,k)))
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qric(i,k)=qric(i,k)+dum*qniic(i,k)
               nric(i,k)=nric(i,k)+dum*nsic(i,k)
               qniic(i,k)=(1._r8-dum)*qniic(i,k)
               nsic(i,k)=(1._r8-dum)*nsic(i,k)
               ! heating tendency 
               tmp=-xlf*dum*qstot/(dz(i,k)*rho(i,k))
               meltsdt(i,k)=meltsdt(i,k) + tmp

               tlat(i,k)=tlat(i,k)+tmp
               qrtot=qrtot+dum*qstot
               nrtot=nrtot+dum*nstot
               qstot=(1._r8-dum)*qstot
               nstot=(1._r8-dum)*nstot
               preci(i)=(1._r8-dum)*preci(i)
            end if
         end if

         ! freeze all rain at -5C for Arctic

         if (t(i,k)+tlat(i,k)/cpp*deltat < (tmelt - 5._r8)) then

            if (qrtot > 0._r8) then

               ! make sure freezing rain doesn't increase temperature above threshold
               dum = xlf/cpp*qrtot/(dz(i,k)*rho(i,k))
               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.(tmelt - 5._r8)) then
                  dum = -(t(i,k)+tlat(i,k)/cpp*deltat-(tmelt-5._r8))*cpp/xlf
                  dum = dum/(xlf/cpp*qrtot/(dz(i,k)*rho(i,k)))
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qniic(i,k)=qniic(i,k)+dum*qric(i,k)
               nsic(i,k)=nsic(i,k)+dum*nric(i,k)
               qric(i,k)=(1._r8-dum)*qric(i,k)
               nric(i,k)=(1._r8-dum)*nric(i,k)
               ! heating tendency 
               tmp = xlf*dum*qrtot/(dz(i,k)*rho(i,k))
               frzrdt(i,k)=frzrdt(i,k) + tmp

               tlat(i,k)=tlat(i,k)+tmp
               qstot=qstot+dum*qrtot
               qrtot=(1._r8-dum)*qrtot
               nstot=nstot+dum*nrtot
               nrtot=(1._r8-dum)*nrtot
               preci(i)=preci(i)+dum*(prect(i)-preci(i))
            end if
         end if

         ! if rain/snow mix ratio is zero so should number concentration

         if (qniic(i,k).lt.qsmall) then
            qniic(i,k)=0._r8
            nsic(i,k)=0._r8
         end if

         if (qric(i,k).lt.qsmall) then
            qric(i,k)=0._r8
            nric(i,k)=0._r8
         end if

         ! make sure number concentration is a positive number to avoid 
         ! taking root of negative

         nric(i,k)=max(nric(i,k),0._r8)
         nsic(i,k)=max(nsic(i,k),0._r8)

         !.......................................................................
         ! get size distribution parameters for fallspeed calculations
         !......................................................................
         ! rain

         if (qric(i,k).ge.qsmall) then
            lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
            n0r(k) = nric(i,k)*lamr(k)

            ! check for slope
            ! change lammax and lammin for rain and snow
            ! adjust vars

            if (lamr(k).lt.lamminr) then

               lamr(k) = lamminr

               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
            else if (lamr(k).gt.lammaxr) then
               lamr(k) = lammaxr
               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
            end if


            ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

            unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
            umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(i,k))

         else
            lamr(k) = 0._r8
            n0r(k) = 0._r8
            umr(k)=0._r8
            unr(k)=0._r8
         end if

         !calculate mean size of combined rain and snow

         if (lamr(k).gt.0._r8) then
            Artmp = n0r(k) * pi / (2._r8 * lamr(k)**3._r8)
         else 
            Artmp = 0._r8
         endif

         if (lamc(k).gt.0._r8) then
            Actmp = cdist1(k) * pi * gamma(pgam(k)+3._r8)/(4._r8 * lamc(k)**2._r8)
         else 
            Actmp = 0._r8
         endif

         if (Actmp.gt.0_r8.or.Artmp.gt.0) then
            rercld(i,k)=rercld(i,k) + 3._r8 *(qric(i,k) + qcic(i,k)) / (4._r8 * rhow * (Actmp + Artmp))
            arcld(i,k)=arcld(i,k)+1._r8
         endif

         !......................................................................
         ! snow

         if (qniic(i,k).ge.qsmall) then
            lams(k) = (cons6*cs*nsic(i,k)/ &
                 qniic(i,k))**(1._r8/ds)
            n0s(k) = nsic(i,k)*lams(k)

            ! check for slope
            ! adjust vars

            if (lams(k).lt.lammins) then
               lams(k) = lammins
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
               nsic(i,k) = n0s(k)/lams(k)

            else if (lams(k).gt.lammaxs) then
               lams(k) = lammaxs
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
               nsic(i,k) = n0s(k)/lams(k)
            end if

            ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)

            ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(i,k))
            uns(k) = min(asn(i,k)*cons7/lams(k)**bs,1.2_r8*rhof(i,k))

         else
            lams(k) = 0._r8
            n0s(k) = 0._r8
            ums(k) = 0._r8
            uns(k) = 0._r8
         end if

         !c........................................................................
         ! sum over sub-step for average process rates

         ! convert rain/snow q and N for output to history, note, 
         ! output is for gridbox average

         qrout(i,k)=qrout(i,k)+qric(i,k)*cldmax(i,k)
         qsout(i,k)=qsout(i,k)+qniic(i,k)*cldmax(i,k)
         nrout(i,k)=nrout(i,k)+nric(i,k)*rho(i,k)*cldmax(i,k)
         nsout(i,k)=nsout(i,k)+nsic(i,k)*rho(i,k)*cldmax(i,k)

         tlat1(i,k)=tlat1(i,k)+tlat(i,k)
         qvlat1(i,k)=qvlat1(i,k)+qvlat(i,k)
         qctend1(i,k)=qctend1(i,k)+qctend(i,k)
         qitend1(i,k)=qitend1(i,k)+qitend(i,k)
         nctend1(i,k)=nctend1(i,k)+nctend(i,k)
         nitend1(i,k)=nitend1(i,k)+nitend(i,k)

         t(i,k)=t(i,k)+tlat(i,k)*deltat/cpp
         q(i,k)=q(i,k)+qvlat(i,k)*deltat
         qc(i,k)=qc(i,k)+qctend(i,k)*deltat
         qi(i,k)=qi(i,k)+qitend(i,k)*deltat
         nc(i,k)=nc(i,k)+nctend(i,k)*deltat
         ni(i,k)=ni(i,k)+nitend(i,k)*deltat

         rainrt1(i,k)=rainrt1(i,k)+rainrt(i,k)

         !divide rain radius over substeps for average
         if (arcld(i,k) .gt. 0._r8) then
            rercld(i,k)=rercld(i,k)/arcld(i,k)
         end if

         !calculate precip fluxes and adding them to summing sub-stepping variables
         !! flux is zero at top interface
         rflx(i,1)=0.0_r8
         sflx(i,1)=0.0_r8

         !! calculating the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
         rflx(i,k+1)=qrout(i,k)*rho(i,k)*umr(k)
         sflx(i,k+1)=qsout(i,k)*rho(i,k)*ums(k)

         !! add to summing sub-stepping variable
         rflx1(i,k+1)=rflx1(i,k+1)+rflx(i,k+1)
         sflx1(i,k+1)=sflx1(i,k+1)+sflx(i,k+1)

         !c........................................................................

      end do ! k loop

      prect1(i)=prect1(i)+prect(i)
      preci1(i)=preci1(i)+preci(i)

   end do ! it loop, sub-step

   do k = top_lev, pver
      rate1ord_cw2pr_st(i,k) = qcsinksum_rate1ord(k)/max(qcsum_rate1ord(k),1.0e-30_r8) 
   end do

300 continue  ! continue if no cloud water
end do ! i loop

! convert dt from sub-step back to full time step
deltat=deltat*real(iter)

!c.............................................................................

do i=1,ncol

   ! skip all calculations if no cloud water
   if (ltrue(i).eq.0) then

      do k=1,top_lev-1
         ! assign zero values for effective radius above 1 mbar
         effc(i,k)=0._r8
         effi(i,k)=0._r8
         effc_fn(i,k)=0._r8
         lamcrad(i,k)=0._r8
         pgamrad(i,k)=0._r8
         deffi(i,k)=0._r8
      end do

      do k=top_lev,pver
         ! assign default values for effective radius
         effc(i,k)=10._r8
         effi(i,k)=25._r8
         effc_fn(i,k)=10._r8
         lamcrad(i,k)=0._r8
         pgamrad(i,k)=0._r8
         deffi(i,k)=0._r8
      end do
      goto 500
   end if

   ! initialize nstep for sedimentation sub-steps
   nstep = 1

   ! divide precip rate by number of sub-steps to get average over time step

   prect(i)=prect1(i)/real(iter)
   preci(i)=preci1(i)/real(iter)

   do k=top_lev,pver

      ! assign variables back to start-of-timestep values before updating after sub-steps 

      t(i,k)=t1(i,k)
      q(i,k)=q1(i,k)
      qc(i,k)=qc1(i,k)
      qi(i,k)=qi1(i,k)
      nc(i,k)=nc1(i,k)
      ni(i,k)=ni1(i,k)

      ! divide microphysical tendencies by number of sub-steps to get average over time step

      tlat(i,k)=tlat1(i,k)/real(iter)
      qvlat(i,k)=qvlat1(i,k)/real(iter)
      qctend(i,k)=qctend1(i,k)/real(iter)
      qitend(i,k)=qitend1(i,k)/real(iter)
      nctend(i,k)=nctend1(i,k)/real(iter)
      nitend(i,k)=nitend1(i,k)/real(iter)

      rainrt(i,k)=rainrt1(i,k)/real(iter)

      ! divide by number of sub-steps to find final values
      rflx(i,k+1)=rflx1(i,k+1)/real(iter)
      sflx(i,k+1)=sflx1(i,k+1)/real(iter)

      ! divide output precip q and N by number of sub-steps to get average over time step

      qrout(i,k)=qrout(i,k)/real(iter)
      qsout(i,k)=qsout(i,k)/real(iter)
      nrout(i,k)=nrout(i,k)/real(iter)
      nsout(i,k)=nsout(i,k)/real(iter)

      ! divide trop_mozart variables by number of sub-steps to get average over time step 

      nevapr(i,k) = nevapr(i,k)/real(iter)
      nevapr2(i,k) = nevapr2(i,k)/real(iter)
      evapsnow(i,k) = evapsnow(i,k)/real(iter)
      prain(i,k) = prain(i,k)/real(iter)
      prodsnow(i,k) = prodsnow(i,k)/real(iter)
      cmeout(i,k) = cmeout(i,k)/real(iter)

      cmeiout(i,k) = cmeiout(i,k)/real(iter)
      meltsdt(i,k) = meltsdt(i,k)/real(iter)
      frzrdt (i,k) = frzrdt (i,k)/real(iter)


      ! microphysics output
      prao(i,k)=prao(i,k)/real(iter)
      prco(i,k)=prco(i,k)/real(iter)
      mnuccco(i,k)=mnuccco(i,k)/real(iter)
      mnuccto(i,k)=mnuccto(i,k)/real(iter)
      msacwio(i,k)=msacwio(i,k)/real(iter)
      psacwso(i,k)=psacwso(i,k)/real(iter)
      bergso(i,k)=bergso(i,k)/real(iter)
      bergo(i,k)=bergo(i,k)/real(iter)
      prcio(i,k)=prcio(i,k)/real(iter)
      praio(i,k)=praio(i,k)/real(iter)

      mnuccro(i,k)=mnuccro(i,k)/real(iter)
      pracso (i,k)=pracso (i,k)/real(iter)

      mnuccdo(i,k)=mnuccdo(i,k)/real(iter)

      ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
      nevapr(i,k) = nevapr(i,k) + evapsnow(i,k)
      prer_evap(i,k) = nevapr2(i,k)
      prain(i,k) = prain(i,k) + prodsnow(i,k)

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! calculate sedimentation for cloud water and ice
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! update in-cloud cloud mixing ratio and number concentration 
      ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
      ! note: these are in-cloud values***, hence we divide by cloud fraction

      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)
      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)/icldm(i,k)
      dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k),0._r8)
      dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)/icldm(i,k),0._r8)

      ! obtain new slope parameter to avoid possible singularity

      if (dumi(i,k).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)

         lami(k) = (cons1*ci* &
              dumni(i,k)/dumi(i,k))**(1._r8/di)
         lami(k)=max(lami(k),lammini)
         lami(k)=min(lami(k),lammaxi)
      else
         lami(k)=0._r8
      end if

      if (dumc(i,k).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)
         ! add lower limit to in-cloud number concentration
         dumnc(i,k)=max(dumnc(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3 
         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         lammax = (pgam(k)+1._r8)/2.e-6_r8
         lamc(k)=max(lamc(k),lammin)
         lamc(k)=min(lamc(k),lammax)
      else
         lamc(k)=0._r8
      end if

      ! calculate number and mass weighted fall velocity for droplets
      ! include effects of sub-grid distribution of cloud water


      if (dumc(i,k).ge.qsmall) then
         unc = acn(i,k)*gamma(1._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+1._r8))
         umc = acn(i,k)*gamma(4._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+4._r8))
         ! fallspeed for output
         vtrmc(i,k)=umc
      else
         umc = 0._r8
         unc = 0._r8
      end if

      ! calculate number and mass weighted fall velocity for cloud ice

      if (dumi(i,k).ge.qsmall) then
         uni =  ain(i,k)*cons16/lami(k)**bi
         umi = ain(i,k)*cons17/(6._r8*lami(k)**bi)
         uni=min(uni,1.2_r8*rhof(i,k))
         umi=min(umi,1.2_r8*rhof(i,k))

         ! fallspeed
         vtrmi(i,k)=umi
      else
         umi = 0._r8
         uni = 0._r8
      end if

      fi(k) = g*rho(i,k)*umi
      fni(k) = g*rho(i,k)*uni
      fc(k) = g*rho(i,k)*umc
      fnc(k) = g*rho(i,k)*unc

      ! calculate number of split time steps to ensure courant stability criteria
      ! for sedimentation calculations

      rgvm = max(fi(k),fc(k),fni(k),fnc(k))
      nstep = max(int(rgvm*deltat/pdel(i,k)+1._r8),nstep)

      ! redefine dummy variables - sedimentation is calculated over grid-scale
      ! quantities to ensure conservation

      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
      dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
      dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)

      if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
      if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8

   end do       !!! vertical loop
   do n = 1,nstep  !! loop over sub-time step to ensure stability

      do k = top_lev,pver
         if (do_cldice) then
            falouti(k) = fi(k)*dumi(i,k)
            faloutni(k) = fni(k)*dumni(i,k)
         else
            falouti(k)  = 0._r8
            faloutni(k) = 0._r8
         end if

         faloutc(k) = fc(k)*dumc(i,k)
         faloutnc(k) = fnc(k)*dumnc(i,k)
      end do

      ! top of model

      k = top_lev
      faltndi = falouti(k)/pdel(i,k)
      faltndni = faloutni(k)/pdel(i,k)
      faltndc = faloutc(k)/pdel(i,k)
      faltndnc = faloutnc(k)/pdel(i,k)

      ! add fallout terms to microphysical tendencies

      qitend(i,k) = qitend(i,k)-faltndi/nstep
      nitend(i,k) = nitend(i,k)-faltndni/nstep
      qctend(i,k) = qctend(i,k)-faltndc/nstep
      nctend(i,k) = nctend(i,k)-faltndnc/nstep

      ! sedimentation tendencies for output
      qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
      qisedten(i,k)=qisedten(i,k)-faltndi/nstep

      dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
      dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
      dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
      dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

      do k = top_lev+1,pver

         ! for cloud liquid and ice, if cloud fraction increases with height
         ! then add flux from above to both vapor and cloud water of current level
         ! this means that flux entering clear portion of cell from above evaporates
         ! instantly

         dum=lcldm(i,k)/lcldm(i,k-1)
         dum=min(dum,1._r8)
         dum1=icldm(i,k)/icldm(i,k-1)
         dum1=min(dum1,1._r8)

         faltndqie=(falouti(k)-falouti(k-1))/pdel(i,k)
         faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(i,k)
         faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(i,k)
         faltndqce=(faloutc(k)-faloutc(k-1))/pdel(i,k)
         faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(i,k)
         faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(i,k)

         ! add fallout terms to eulerian tendencies

         qitend(i,k) = qitend(i,k)-faltndi/nstep
         nitend(i,k) = nitend(i,k)-faltndni/nstep
         qctend(i,k) = qctend(i,k)-faltndc/nstep
         nctend(i,k) = nctend(i,k)-faltndnc/nstep

         ! sedimentation tendencies for output
         qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
         qisedten(i,k)=qisedten(i,k)-faltndi/nstep

         ! add terms to to evap/sub of cloud water

         qvlat(i,k)=qvlat(i,k)-(faltndqie-faltndi)/nstep
         ! for output
         qisevap(i,k)=qisevap(i,k)-(faltndqie-faltndi)/nstep
         qvlat(i,k)=qvlat(i,k)-(faltndqce-faltndc)/nstep
         ! for output
         qcsevap(i,k)=qcsevap(i,k)-(faltndqce-faltndc)/nstep

         tlat(i,k)=tlat(i,k)+(faltndqie-faltndi)*xxls/nstep
         tlat(i,k)=tlat(i,k)+(faltndqce-faltndc)*xxlv/nstep

         dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
         dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
         dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
         dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

         Fni(K)=MAX(Fni(K)/pdel(i,K),Fni(K-1)/pdel(i,K-1))*pdel(i,K)
         FI(K)=MAX(FI(K)/pdel(i,K),FI(K-1)/pdel(i,K-1))*pdel(i,K)
         fnc(k)=max(fnc(k)/pdel(i,k),fnc(k-1)/pdel(i,k-1))*pdel(i,k)
         Fc(K)=MAX(Fc(K)/pdel(i,K),Fc(K-1)/pdel(i,K-1))*pdel(i,K)

      end do   !! k loop

      ! units below are m/s
      ! cloud water/ice sedimentation flux at surface 
      ! is added to precip flux at surface to get total precip (cloud + precip water)
      ! rate

      prect(i) = prect(i)+(faloutc(pver)+falouti(pver))/g/nstep/1000._r8  
      preci(i) = preci(i)+(falouti(pver))/g/nstep/1000._r8

   end do   !! nstep loop

   ! end sedimentation
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   ! get new update for variables that includes sedimentation tendency
   ! note : here dum variables are grid-average, NOT in-cloud

   do k=top_lev,pver

      dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)
      dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)
      dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)
      dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)

      if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
      if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8

      ! calculate instantaneous processes (melting, homogeneous freezing)
      if (do_cldice) then

         if (t(i,k)+tlat(i,k)/cpp*deltat > tmelt) then
            if (dumi(i,k) > 0._r8) then

               ! limit so that melting does not push temperature below freezing
               dum = -dumi(i,k)*xlf/cpp
               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.tmelt) then
                  dum = (t(i,k)+tlat(i,k)/cpp*deltat-tmelt)*cpp/xlf
                  dum = dum/dumi(i,k)*xlf/cpp 
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qctend(i,k)=qctend(i,k)+dum*dumi(i,k)/deltat

               ! for output
               melto(i,k)=dum*dumi(i,k)/deltat

               ! assume melting ice produces droplet
               ! mean volume radius of 8 micron

               nctend(i,k)=nctend(i,k)+3._r8*dum*dumi(i,k)/deltat/ &
                    (4._r8*pi*5.12e-16_r8*rhow)

               qitend(i,k)=((1._r8-dum)*dumi(i,k)-qi(i,k))/deltat
               nitend(i,k)=((1._r8-dum)*dumni(i,k)-ni(i,k))/deltat
               tlat(i,k)=tlat(i,k)-xlf*dum*dumi(i,k)/deltat
            end if
         end if

         ! homogeneously freeze droplets at -40 C

         if (t(i,k)+tlat(i,k)/cpp*deltat < 233.15_r8) then
            if (dumc(i,k) > 0._r8) then

               ! limit so that freezing does not push temperature above threshold
               dum = dumc(i,k)*xlf/cpp
               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.233.15_r8) then
                  dum = -(t(i,k)+tlat(i,k)/cpp*deltat-233.15_r8)*cpp/xlf
                  dum = dum/dumc(i,k)*xlf/cpp
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qitend(i,k)=qitend(i,k)+dum*dumc(i,k)/deltat
               ! for output
               homoo(i,k)=dum*dumc(i,k)/deltat

               ! assume 25 micron mean volume radius of homogeneously frozen droplets
               ! consistent with size of detrained ice in stratiform.F90
               nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*1.563e-14_r8* &
                    500._r8)/deltat
               qctend(i,k)=((1._r8-dum)*dumc(i,k)-qc(i,k))/deltat
               nctend(i,k)=((1._r8-dum)*dumnc(i,k)-nc(i,k))/deltat
               tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)/deltat
            end if
         end if

         ! remove any excess over-saturation, which is possible due to non-linearity when adding 
         ! together all microphysical processes
         ! follow code similar to old CAM scheme

         qtmp=q(i,k)+qvlat(i,k)*deltat
         ttmp=t(i,k)+tlat(i,k)/cpp*deltat

         esn = svp_water(ttmp)  ! use rhw to allow ice supersaturation
         qsn = svp_to_qsat(esn, p(i,k))

         if (qtmp > qsn .and. qsn > 0) then
            ! expression below is approximate since there may be ice deposition
            dum = (qtmp-qsn)/(1._r8+cons27*qsn/(cpp*rv*ttmp**2))/deltat
            ! add to output cme
            cmeout(i,k) = cmeout(i,k)+dum
            ! now add to tendencies, partition between liquid and ice based on temperature
            if (ttmp > 268.15_r8) then
               dum1=0.0_r8
               ! now add to tendencies, partition between liquid and ice based on te
            else if (ttmp < 238.15_r8) then
               dum1=1.0_r8
            else
               dum1=(268.15_r8-ttmp)/30._r8
            end if

            dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                 *qsn/(cpp*rv*ttmp**2))/deltat
            qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
            ! for output
            qcreso(i,k)=dum*(1._r8-dum1)
            qitend(i,k)=qitend(i,k)+dum*dum1
            qireso(i,k)=dum*dum1
            qvlat(i,k)=qvlat(i,k)-dum
            ! for output
            qvres(i,k)=-dum
            tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
         end if
      end if

      !...............................................................................
      ! calculate effective radius for pass to radiation code
      ! if no cloud water, default value is 10 micron for droplets,
      ! 25 micron for cloud ice

      ! update cloud variables after instantaneous processes to get effective radius
      ! variables are in-cloud to calculate size dist parameters

      dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/lcldm(i,k)
      dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/icldm(i,k)
      dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/lcldm(i,k)
      dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/icldm(i,k)

      ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1

      dumc(i,k)=min(dumc(i,k),5.e-3_r8)
      dumi(i,k)=min(dumi(i,k),5.e-3_r8)

      !...................
      ! cloud ice effective radius

      if (dumi(i,k).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)
         lami(k) = (cons1*ci*dumni(i,k)/dumi(i,k))**(1._r8/di)

         if (lami(k).lt.lammini) then
            lami(k) = lammini
            n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*cons1)
            niic(i,k) = n0i(k)/lami(k)
            ! adjust number conc if needed to keep mean size in reasonable range
            if (do_cldice) nitend(i,k)=(niic(i,k)*icldm(i,k)-ni(i,k))/deltat

         else if (lami(k).gt.lammaxi) then
            lami(k) = lammaxi
            n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*cons1)
            niic(i,k) = n0i(k)/lami(k)
            ! adjust number conc if needed to keep mean size in reasonable range
            if (do_cldice) nitend(i,k)=(niic(i,k)*icldm(i,k)-ni(i,k))/deltat
         end if
         effi(i,k) = 1.5_r8/lami(k)*1.e6_r8

      else
         effi(i,k) = 25._r8
      end if

      ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
      ! radius has already been determined from the size distribution.
      if (.not. do_cldice) then
         effi(i,k) = re_ice(i,k) * 1e6_r8      ! m -> um
      end if

      !...................
      ! cloud droplet effective radius

      if (dumc(i,k).ge.qsmall) then

         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! set tendency to ensure minimum droplet concentration
         ! after update by microphysics, except when lambda exceeds bounds on mean drop
         ! size or if there is no cloud water
         if (dumnc(i,k).lt.cdnl/rho(i,k)) then   
            nctend(i,k)=(cdnl/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat   
         end if
         dumnc(i,k)=max(dumnc(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         ! Multiply by omsm to fit within RRTMG's table.
         lammax = (pgam(k)+1._r8)*omsm/2.e-6_r8
         if (lamc(k).lt.lammin) then
            lamc(k) = lammin
            ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
                 gamma(pgam(k)+1._r8)/ &
                 (pi*rhow*gamma(pgam(k)+4._r8))
            ! adjust number conc if needed to keep mean size in reasonable range
            nctend(i,k)=(ncic(i,k)*lcldm(i,k)-nc(i,k))/deltat

         else if (lamc(k).gt.lammax) then
            lamc(k) = lammax
            ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
                 gamma(pgam(k)+1._r8)/ &
                 (pi*rhow*gamma(pgam(k)+4._r8))
            ! adjust number conc if needed to keep mean size in reasonable range
            nctend(i,k)=(ncic(i,k)*lcldm(i,k)-nc(i,k))/deltat
         end if

         effc(i,k) = &
              gamma(pgam(k)+4._r8)/ &
              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
         !assign output fields for shape here
         lamcrad(i,k)=lamc(k)
         pgamrad(i,k)=pgam(k)

      else
         effc(i,k) = 10._r8
         lamcrad(i,k)=0._r8
         pgamrad(i,k)=0._r8
      end if

      ! ice effective diameter for david mitchell's optics
      if (do_cldice) then
         deffi(i,k)=effi(i,k)*rhoi/917._r8*2._r8
      else
         deffi(i,k)=effi(i,k) * 2._r8
      end if


!!! recalculate effective radius for constant number, in order to separate
      ! first and second indirect effects
      ! assume constant number of 10^8 kg-1

      dumnc(i,k)=1.e8_r8

      if (dumc(i,k).ge.qsmall) then
         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         lammax = (pgam(k)+1._r8)/2.e-6_r8
         if (lamc(k).lt.lammin) then
            lamc(k) = lammin
         else if (lamc(k).gt.lammax) then
            lamc(k) = lammax
         end if
         effc_fn(i,k) = &
              gamma(pgam(k)+4._r8)/ &
              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8

      else
         effc_fn(i,k) = 10._r8
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

   end do ! vertical k loop

500 continue

   do k=top_lev,pver
      ! if updated q (after microphysics) is zero, then ensure updated n is also zero

      if (qc(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-nc(i,k)/deltat
      if (do_cldice .and. qi(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-ni(i,k)/deltat
   end do

end do ! i loop

! add snow ouptut
do i = 1,ncol
   do k=top_lev,pver
      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
         dsout(i,k)=3._r8*rhosn/917._r8*(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
      endif
   end do
end do

!calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual
do i = 1,ncol
   do k=top_lev,pver
      !! RAIN
      if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
         reff_rain(i,k)=1.5_r8*(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)*1.e6_r8
      endif
      !! SNOW
      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
         reff_snow(i,k)=1.5_r8*(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)*1.e6_r8
      end if
   end do
end do

! analytic radar reflectivity
! formulas from Matthew Shupe, NOAA/CERES
! *****note: radar reflectivity is local (in-precip average)
! units of mm^6/m^3

do i = 1,ncol
   do k=top_lev,pver
      if (qc(i,k)+qctend(i,k)*deltat.ge.qsmall) then
         dum=((qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)*1000._r8)**2 &
              /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)/1.e6_r8)*lcldm(i,k)/cldmax(i,k)
      else
         dum=0._r8
      end if
      if (qi(i,k)+qitend(i,k)*deltat.ge.qsmall) then
         dum1=((qi(i,k)+qitend(i,k)*deltat)*rho(i,k)/icldm(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(i,k)/cldmax(i,k)
      else 
         dum1=0._r8
      end if

      if (qsout(i,k).ge.qsmall) then
         dum1=dum1+(qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
      end if

      refl(i,k)=dum+dum1

      ! add rain rate, but for 37 GHz formulation instead of 94 GHz
      ! formula approximated from data of Matrasov (2007)
      ! rainrt is the rain rate in mm/hr
      ! reflectivity (dum) is in DBz

      if (rainrt(i,k).ge.0.001_r8) then
         dum=log10(rainrt(i,k)**6._r8)+16._r8

         ! convert from DBz to mm^6/m^3

         dum = 10._r8**(dum/10._r8)
      else
         ! don't include rain rate in R calculation for values less than 0.001 mm/hr
         dum=0._r8
      end if

      ! add to refl

      refl(i,k)=refl(i,k)+dum

      !output reflectivity in Z.
      areflz(i,k)=refl(i,k)

      ! convert back to DBz 

      if (refl(i,k).gt.minrefl) then 
         refl(i,k)=10._r8*log10(refl(i,k))
      else
         refl(i,k)=-9999._r8
      end if

      !set averaging flag
      if (refl(i,k).gt.mindbz) then 
         arefl(i,k)=refl(i,k)
         frefl(i,k)=1.0_r8  
      else
         arefl(i,k)=0._r8
         areflz(i,k)=0._r8
         frefl(i,k)=0._r8
      end if

      ! bound cloudsat reflectivity

      csrfl(i,k)=min(csmax,refl(i,k))

      !set averaging flag
      if (csrfl(i,k).gt.csmin) then 
         acsrfl(i,k)=refl(i,k)
         fcsrfl(i,k)=1.0_r8  
      else
         acsrfl(i,k)=0._r8
         fcsrfl(i,k)=0._r8
      end if

   end do
end do


! averaging for snow and rain number and diameter

qrout2(:,:)=0._r8
qsout2(:,:)=0._r8
nrout2(:,:)=0._r8
nsout2(:,:)=0._r8
drout2(:,:)=0._r8
dsout2(:,:)=0._r8
freqs(:,:)=0._r8
freqr(:,:)=0._r8
do i = 1,ncol
   do k=top_lev,pver
      if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
         qrout2(i,k)=qrout(i,k)
         nrout2(i,k)=nrout(i,k)
         drout2(i,k)=(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)
         freqr(i,k)=1._r8
      endif
      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
         qsout2(i,k)=qsout(i,k)
         nsout2(i,k)=nsout(i,k)
         dsout2(i,k)=(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
         freqs(i,k)=1._r8
      endif
   end do
end do

! output activated liquid and ice (convert from #/kg -> #/m3)
do i = 1,ncol
   do k=top_lev,pver
      ncai(i,k)=dum2i(i,k)*rho(i,k)
      ncal(i,k)=dum2l(i,k)*rho(i,k)
   end do
end do


!redefine fice here....
nfice(:,:)=0._r8
do k=top_lev,pver
   do i=1,ncol
      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
      dumfice=qsout(i,k) + qrout(i,k) + dumc(i,k) + dumi(i,k)  

      if (dumfice.gt.qsmall.and.(qsout(i,k)+dumi(i,k).gt.qsmall)) then
         nfice(i,k)=(qsout(i,k) + dumi(i,k))/dumfice
      endif

      if (nfice(i,k).gt.1._r8) then
         nfice(i,k)=1._r8
      endif

   enddo
enddo


end subroutine micro_mg_tend

#endif

!========================================================================
!UTILITIES
!========================================================================

pure subroutine micro_mg_get_cols(ncol, nlev, top_lev, qcn, qin, &
     mgncol, mgcols)

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.

  integer, intent(in) :: ncol      ! Number of columns with meaningful data
  integer, intent(in) :: nlev      ! Number of levels to use
  integer, intent(in) :: top_lev   ! Top level for microphysics

  real(r8), intent(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)

  integer, intent(out) :: mgncol   ! Number of columns MG will use
  integer, allocatable, intent(out) :: mgcols(:) ! column indices

  integer :: lev_offset  ! top_lev - 1 (defined here for consistency)
  logical :: ltrue(ncol) ! store tests for each column

  integer :: i, ii ! column indices

  if (allocated(mgcols)) deallocate(mgcols)

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know if water is present at any level
  ! in each column.

  ltrue = any(qcn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qin(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)

  ! Scan for true values to get a usable list of indices.

  mgncol = count(ltrue)
  allocate(mgcols(mgncol))
  i = 0
  do ii = 1,ncol
     if (ltrue(ii)) then
        i = i + 1
        mgcols(i) = ii
     end if
  end do

end subroutine micro_mg_get_cols

end module micro_mg1_0
