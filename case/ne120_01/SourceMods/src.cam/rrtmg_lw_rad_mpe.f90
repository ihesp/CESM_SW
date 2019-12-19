!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.6 $
!     created:   $Date: 2008/04/24 16:17:27 $
!

       module rrtmg_lw_rad_mpe

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                              RRTMG_LW                                    *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                   a rapid radiative transfer model                       *
! *                       for the longwave region                            * 
! *             for application to general circulation models                *
! *                                                                          *
! *                                                                          *
! *            Atmospheric and Environmental Research, Inc.                  *
! *                        131 Hartwell Avenue                               *
! *                        Lexington, MA 02421                               *
! *                                                                          *
! *                                                                          *
! *                           Eli J. Mlawer                                  *
! *                        Jennifer S. Delamere                              *
! *                         Michael J. Iacono                                *
! *                         Shepard A. Clough                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                       email:  miacono@aer.com                            *
! *                       email:  emlawer@aer.com                            *
! *                       email:  jdelamer@aer.com                           *
! *                                                                          *
! *        The authors wish to acknowledge the contributions of the          *
! *        following people:  Steven J. Taubman, Karen Cady-Pereira,         *
! *        Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert Bergstrom.  *
! *                                                                          *
! ****************************************************************************

! -------- Modules --------

      use rrtmg_lw_rad, only: inatm
      use shr_kind_mod, only: r8 => shr_kind_r8
      use ppgrid,       only: pcols, begchunk, endchunk

!      use parkind, only : jpim, jprb 
      use rrlw_vsn
      use mcica_subcol_gen_lw, only: mcica_subcol_lw
      use rrtmg_lw_cldprmc, only: cldprmc
! Move call to rrtmg_lw_ini and following use association to 
! GCM initialization area
!      use rrtmg_lw_init, only: rrtmg_lw_ini
      use rrtmg_lw_rtrnmc, only: rtrnmc
      use rrtmg_lw_setcoef, only: setcoef
      use rrtmg_lw_taumol, only: taumol

      implicit none

#include <mpif.h>
! public interfaces/functions/subroutines
      public :: rrtmg_lw_parallel

      type rrtmg_lw_param
          integer :: lchnk                      ! chunk identifier
          integer :: ncol                       ! Number of horizontal columns
          integer :: nlay                       ! Number of model layers
          integer :: icld                    ! Cloud overlap method
                                                            !    0: Clear only
                                                            !    1: Random
                                                            !    2: Maximum/random
                                                            !    3: Maximum
          real(kind=r8), pointer :: play(:,:)            ! Layer pressures (hPa, mb)
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                            !    Dimensions: (ncol,nlay+1)
          real(kind=r8), pointer :: tlay(:,:)            ! Layer temperatures (K)
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: tlev(:,:)            ! Interface temperatures (K)
                                                            !    Dimensions: (ncol,nlay+1)
          real(kind=r8), pointer :: tsfc(:)              ! Surface temperature (K)
                                                            !    Dimensions: (ncol)
          real(kind=r8), pointer :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: cfc11vmr(:,:)        ! CFC11 volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: cfc12vmr(:,:)        ! CFC12 volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: cfc22vmr(:,:)        ! CFC22 volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: ccl4vmr(:,:)         ! CCL4 volume mixing ratio
                                                            !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: emis(:,:)            ! Surface emissivity
                                                            !    Dimensions: (ncol,nbndlw)

          integer :: inflglw                    ! Flag for cloud optical properties
          integer :: iceflglw                   ! Flag for ice particle specification
          integer :: liqflglw                   ! Flag for liquid droplet specification

          real(kind=r8), pointer :: cldfmcl(:,:,:)       ! Cloud fraction
                                           !    Dimensions: (ngptlw,ncol,nlay)
          real(kind=r8), pointer :: ciwpmcl(:,:,:)       ! Cloud ice water path (g/m2)
                                           !    Dimensions: (ngptlw,ncol,nlay)
          real(kind=r8), pointer :: clwpmcl(:,:,:)       ! Cloud liquid water path (g/m2)
                                           !    Dimensions: (ngptlw,ncol,nlay)
          real(kind=r8), pointer :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                           !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                           !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: taucmcl(:,:,:)       ! Cloud optical depth
                                           !    Dimensions: (ngptlw,ncol,nlay)
    !      real(kind=r8), pointer :: ssacmcl(:,:,:)      ! Cloud single scattering albedo
                                         !    Dimensions: (ngptlw,ncol,nlay)
                                         !   for future expansion
                                         !   lw scattering not yet available
    !      real(kind=r8), pointer :: asmcmcl(:,:,:)      ! Cloud asymmetry parameter
                                         !    Dimensions: (ngptlw,ncol,nlay)
                                         !   for future expansion
                                         !   lw scattering not yet available
          real(kind=r8), pointer :: tauaer(:,:,:)        ! aerosol optical depth
                                         !   at mid-point of LW spectral bands
                                         !    Dimensions: (ncol,nlay,nbndlw)
    !      real(kind=r8), pointer :: ssaaer(:,:,:)       ! aerosol single scattering albedo
                                         !    Dimensions: (ncol,nlay,nbndlw)
                                         !   for future expansion 
                                         !   (lw aerosols/scattering not yet available)
    !      real(kind=r8), pointer :: asmaer(:,:,:)       ! aerosol asymmetry parameter
                                         !    Dimensions: (ncol,nlay,nbndlw)
                                         !   for future expansion 
                                         !   (lw aerosols/scattering not yet available)

    ! ----- Output -----

          real(kind=r8), pointer :: uflx(:,:)           ! Total sky longwave upward flux (W/m2)
                                        !    Dimensions: (ncol,nlay+1)
          real(kind=r8), pointer :: dflx(:,:)           ! Total sky longwave downward flux (W/m2)
                                        !    Dimensions: (ncol,nlay+1)
          real(kind=r8), pointer :: hr(:,:)             ! Total sky longwave radiative heating rate (K/d)
                                        !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: uflxc(:,:)          ! Clear sky longwave upward flux (W/m2)
                                        !    Dimensions: (ncol,nlay+1)
          real(kind=r8), pointer :: dflxc(:,:)          ! Clear sky longwave downward flux (W/m2)
                                        !    Dimensions: (ncol,nlay+1)
          real(kind=r8), pointer :: hrc(:,:)            ! Clear sky longwave radiative heating rate (K/d)
                                        !    Dimensions: (ncol,nlay)
          real(kind=r8), pointer :: uflxs(:,:,:)        ! Total sky longwave upward flux spectral (W/m2)
                                        !    Dimensions: (nbndlw,ncol,nlay+1)
          real(kind=r8), pointer :: dflxs(:,:,:)        ! Total sky longwave downward flux spectral (W/m2)
                                        !    Dimensions: (nbndlw,ncol,nlay+1)
    ! ----- JFLFY -----

          integer :: mpi_rank
          integer(8) :: end_type
      end type

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_lw_parallel &
            (lchnk   ,ncol    ,nlay    ,icld    ,                   &
             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,&
             cfc11vmr,cfc12vmr, &
             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
             cldfmcl ,taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, uflxs, dflxs )

      use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol

! ----- Input -----
      integer, intent(in) :: lchnk                      ! chunk identifier
      integer, intent(in) :: ncol                       ! Number of horizontal columns
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(inout) :: icld                    ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum
      real(kind=r8), intent(in), target :: play(:,:)            ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in), target :: tlay(:,:)            ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: tlev(:,:)            ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in), target :: tsfc(:)              ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in), target :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: cfc11vmr(:,:)        ! CFC11 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: cfc12vmr(:,:)        ! CFC12 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: cfc22vmr(:,:)        ! CFC22 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: ccl4vmr(:,:)         ! CCL4 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: emis(:,:)            ! Surface emissivity
                                                        !    Dimensions: (ncol,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in), target :: cldfmcl(:,:,:)       ! Cloud fraction
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in), target :: ciwpmcl(:,:,:)       ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in), target :: clwpmcl(:,:,:)       ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in), target :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in), target :: taucmcl(:,:,:)       ! Cloud optical depth
                                                        !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=r8), intent(in), target :: ssacmcl(:,:,:)      ! Cloud single scattering albedo
                                                        !    Dimensions: (ngptlw,ncol,nlay)
                                                        !   for future expansion
                                                        !   lw scattering not yet available
!      real(kind=r8), intent(in), target :: asmcmcl(:,:,:)      ! Cloud asymmetry parameter
                                                        !    Dimensions: (ngptlw,ncol,nlay)
                                                        !   for future expansion
                                                        !   lw scattering not yet available
      real(kind=r8), intent(in), target :: tauaer(:,:,:)        ! aerosol optical depth
                                                        !   at mid-point of LW spectral bands
                                                        !    Dimensions: (ncol,nlay,nbndlw)
!      real(kind=r8), intent(in), target :: ssaaer(:,:,:)       ! aerosol single scattering albedo
                                                        !    Dimensions: (ncol,nlay,nbndlw)
                                                        !   for future expansion 
                                                        !   (lw aerosols/scattering not yet available)
!      real(kind=r8), intent(in), target :: asmaer(:,:,:)       ! aerosol asymmetry parameter
                                                        !    Dimensions: (ncol,nlay,nbndlw)
                                                        !   for future expansion 
                                                        !   (lw aerosols/scattering not yet available)

! ----- Output -----

      real(kind=r8), intent(out), target :: uflx(:,:)           ! Total sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out), target :: dflx(:,:)           ! Total sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out), target :: hr(:,:)             ! Total sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(out), target :: uflxc(:,:)          ! Clear sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out), target :: dflxc(:,:)          ! Clear sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out), target :: hrc(:,:)            ! Clear sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(out), target :: uflxs(:,:,:)        ! Total sky longwave upward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)
      real(kind=r8), intent(out), target :: dflxs(:,:,:)        ! Total sky longwave downward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)

      type(rrtmg_lw_param) :: param

      integer, external :: slave_setftz, slave_unsetftz
      integer :: mpi_rank, ierr

      call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)

      param%lchnk = lchnk
      param%ncol = ncol
      param%nlay = nlay
      param%icld = icld
      param%play => play
      param%plev => plev
      param%tlay => tlay
      param%tlev => tlev
      param%tsfc => tsfc
      param%h2ovmr => h2ovmr
      param%o3vmr => o3vmr
      param%co2vmr => co2vmr
      param%ch4vmr => ch4vmr
      param%o2vmr => o2vmr
      param%n2ovmr => n2ovmr
      param%cfc11vmr => cfc11vmr
      param%cfc12vmr => cfc12vmr
      param%cfc22vmr => cfc22vmr
      param%ccl4vmr => ccl4vmr
      param%emis => emis
      param%inflglw = inflglw
      param%iceflglw = iceflglw
      param%liqflglw = liqflglw
      param%cldfmcl => cldfmcl
      param%ciwpmcl => ciwpmcl
      param%clwpmcl => clwpmcl
      param%reicmcl => reicmcl
      param%relqmcl => relqmcl
      param%taucmcl => taucmcl
      param%tauaer => tauaer
      param%uflx => uflx
      param%dflx => dflx
      param%hr => hr
      param%uflxc => uflxc
      param%dflxc => dflxc
      param%hrc => hrc
      param%uflxs => uflxs
      param%dflxs => dflxs
      param%mpi_rank = mpi_rank

      call athread_spawn(slave_setftz, 0)
      call athread_join

      call athread_spawn_rrtmg_lw(param)
      call athread_join

      call athread_spawn(slave_unsetftz, 0)
      call athread_join

      end subroutine rrtmg_lw_parallel

      subroutine rrtmg_lw_mpe (param)

! -------- Description --------

! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
! model for application to GCMs, that has been adapted from RRTM_LW for
! improved efficiency.
!
! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
!  area, since this has to be called only once. 
!
! This routine:
!    a) calls INATM to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    b) calls CLDPRMC to set cloud optical depth for McICA based 
!       on input cloud properties 
!    c) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    d) calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands
!    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
!       radiative transfer calculation using McICA, the Monte-Carlo 
!       Independent Column Approximation, to represent sub-grid scale 
!       cloud variability
!    f) passes the necessary fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
!     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM. 
!
!    1) Standard, single forward model calculation (imca = 0)
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!
! This call to RRTMG_LW must be preceeded by a call to the module
!     mcica_subcol_gen_lw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngpt) dimension.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags param%inflglw, param%iceflglw, and param%liqflglw; see text file rrtmg_lw_instructions
!     and subroutine rrtmg_lw_cldprop.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (param%inflglw = 0)
!    2) Input cloud fraction and cloud physical properties (param%inflglw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of param%iceflglw and param%liqflglw
!
! One method of aerosol property input is possible:
!     Aerosol properties can be input in only one way (controlled by input 
!     flag iaer, see text file rrtmg_lw_instructions for further details):
!
!    1) Input aerosol optical depth directly by layer and spectral band (iaer=10);
!       band average optical depth at the mid-point of each spectral band.
!       RRTMG_LW currently treats only aerosol absorption;
!       scattering capability is not presently available. 
!
!
! ------- Modifications -------
!
! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
! set of g-points for application to GCMs.  
!
!-- Original version (derived from RRTM_LW), reduction of g-points, other
!   revisions for use with GCMs.  
!     1999: M. J. Iacono, AER, Inc.
!-- Adapted for use with NCAR/CAM.
!     May 2004: M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability. 
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Conversion to F90 formatting for consistency with rrtmg_sw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays.
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to add longwave aerosol absorption.
!     Apr 2008: M. J. Iacono, AER, Inc.

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi
      use rrlw_wvn, only: ng, ngb, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------
      type(rrtmg_lw_param), intent(inout) :: param

! ----- Local -----

! Control
      integer :: istart                         ! beginning band of calculation
      integer :: iend                           ! ending band of calculation
      integer :: iout                           ! output option flag (inactive)
      integer :: iaer                           ! aerosol option flag
      integer :: iplon                          ! column loop index
      integer :: imca                           ! flag for mcica [0=off, 1=on]
      integer :: ims                            ! value for changing mcica permute seed
      integer :: k                              ! layer loop index
      integer :: ig                             ! g-point loop index

! Atmosphere
      real(kind=r8) :: pavel(param%nlay)              ! layer pressures (mb) 
      real(kind=r8) :: tavel(param%nlay)              ! layer temperatures (K)
      real(kind=r8) :: pz(0:param%nlay)               ! level (interface) pressures (hPa, mb)
      real(kind=r8) :: tz(0:param%nlay)               ! level (interface) temperatures (K)
      real(kind=r8) :: tbound                   ! surface temperature (K)
      real(kind=r8) :: coldry(param%nlay)             ! dry air column density (mol/cm2)
      real(kind=r8) :: wbrodl(param%nlay)             ! broadening gas column density (mol/cm2)
      real(kind=r8) :: wkl(mxmol,param%nlay)          ! molecular amounts (mol/cm-2)
      real(kind=r8) :: wx(maxxsec,param%nlay)         ! cross-section amounts (mol/cm-2)
      real(kind=r8) :: pwvcm                    ! precipitable water vapor (cm)
      real(kind=r8) :: semiss(nbndlw)           ! lw surface emissivity
      real(kind=r8), allocatable :: fracs(:,:)          ! 
      real(kind=r8), allocatable :: taug(:,:)           ! gaseous optical depths
      real(kind=r8), allocatable :: taut(:,:)           ! gaseous + aerosol optical depths

      real(kind=r8) :: taua(param%nlay,nbndlw)        ! aerosol optical depth
!      real(kind=r8) :: ssaa(param%nlay,nbndlw)        ! aerosol single scattering albedo
                                                 !   for future expansion 
                                                 !   (lw aerosols/scattering not yet available)
!      real(kind=r8) :: asma(param%nlay+1,nbndlw)      ! aerosol asymmetry parameter
                                                 !   for future expansion 
                                                 !   (lw aerosols/scattering not yet available)

! Atmosphere - setcoef
      integer :: laytrop                          ! tropopause layer index
      integer :: jp(param%nlay)                         ! lookup table index 
      integer :: jt(param%nlay)                         ! lookup table index 
      integer :: jt1(param%nlay)                        ! lookup table index 
      real(kind=r8) :: planklay(param%nlay,nbndlw)      ! 
      real(kind=r8) :: planklev(0:param%nlay,nbndlw)    ! 
      real(kind=r8) :: plankbnd(nbndlw)           ! 

      real(kind=r8) :: colh2o(param%nlay)               ! column amount (h2o)
      real(kind=r8) :: colco2(param%nlay)               ! column amount (co2)
      real(kind=r8) :: colo3(param%nlay)                ! column amount (o3)
      real(kind=r8) :: coln2o(param%nlay)               ! column amount (n2o)
      real(kind=r8) :: colco(param%nlay)                ! column amount (co)
      real(kind=r8) :: colch4(param%nlay)               ! column amount (ch4)
      real(kind=r8) :: colo2(param%nlay)                ! column amount (o2)
      real(kind=r8) :: colbrd(param%nlay)               ! column amount (broadening gases)

      integer :: indself(param%nlay)
      integer :: indfor(param%nlay)
      real(kind=r8) :: selffac(param%nlay)
      real(kind=r8) :: selffrac(param%nlay)
      real(kind=r8) :: forfac(param%nlay)
      real(kind=r8) :: forfrac(param%nlay)

      integer :: indminor(param%nlay)
      real(kind=r8) :: minorfrac(param%nlay)
      real(kind=r8) :: scaleminor(param%nlay)
      real(kind=r8) :: scaleminorn2(param%nlay)

      real(kind=r8) :: &                          !
                         fac00(param%nlay), fac01(param%nlay), &
                         fac10(param%nlay), fac11(param%nlay) 
      real(kind=r8) :: &                          !
                         rat_h2oco2(param%nlay),rat_h2oco2_1(param%nlay), &
                         rat_h2oo3(param%nlay),rat_h2oo3_1(param%nlay), &
                         rat_h2on2o(param%nlay),rat_h2on2o_1(param%nlay), &
                         rat_h2och4(param%nlay),rat_h2och4_1(param%nlay), &
                         rat_n2oco2(param%nlay),rat_n2oco2_1(param%nlay), &
                         rat_o3co2(param%nlay),rat_o3co2_1(param%nlay)

! Atmosphere/clouds - cldprop
      integer :: ncbands                          ! number of cloud spectral bands
      integer :: inflag                           ! flag for cloud property method
      integer :: iceflag                          ! flag for ice cloud properties
      integer :: liqflag                          ! flag for liquid cloud properties

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=r8), allocatable :: cldfmc(:,:)       ! cloud fraction [mcica]
      real(kind=r8), allocatable :: ciwpmc(:,:)       ! cloud ice water path [mcica]
      real(kind=r8), allocatable :: clwpmc(:,:)       ! cloud liquid water path [mcica]
      real(kind=r8) :: relqmc(param%nlay)             ! liquid particle size (microns)
      real(kind=r8) :: reicmc(param%nlay)             ! ice particle effective radius (microns)
      real(kind=r8) :: dgesmc(param%nlay)             ! ice particle generalized effective size (microns)
      ! real(kind=r8) :: taucmc(ngptlw,param%nlay)      ! cloud optical depth [mcica]
      real(kind=r8), allocatable :: taucmc(:,:)       ! cloud optical depth [mcica]
!      real(kind=r8) :: ssacmc(ngptlw,param%nlay)     ! cloud single scattering albedo [mcica]
                                                !   for future expansion 
                                                !   (lw scattering not yet available)
!      real(kind=r8) :: asmcmc(ngptlw,param%nlay)     ! cloud asymmetry parameter [mcica]
                                                !   for future expansion 
                                                !   (lw scattering not yet available)

! Output
      real(kind=r8) :: totuflux(0:param%nlay)         ! upward longwave flux (w/m2)
      real(kind=r8) :: totdflux(0:param%nlay)         ! downward longwave flux (w/m2)
      real(kind=r8) :: totufluxs(nbndlw,0:param%nlay) ! upward longwave flux spectral (w/m2)
      real(kind=r8) :: totdfluxs(nbndlw,0:param%nlay) ! downward longwave flux spectral (w/m2)
      real(kind=r8) :: fnet(0:param%nlay)             ! net longwave flux (w/m2)
      real(kind=r8) :: htr(0:param%nlay)              ! longwave heating rate (k/day)
      real(kind=r8) :: totuclfl(0:param%nlay)         ! clear sky upward longwave flux (w/m2)
      real(kind=r8) :: totdclfl(0:param%nlay)         ! clear sky downward longwave flux (w/m2)
      real(kind=r8) :: fnetc(0:param%nlay)            ! clear sky net longwave flux (w/m2)
      real(kind=r8) :: htrc(0:param%nlay)             ! clear sky longwave heating rate (k/day)

      allocate(fracs(param%nlay,ngptlw), taug(param%nlay,ngptlw), taut(param%nlay,ngptlw))
      allocate(cldfmc(ngptlw,param%nlay), ciwpmc(ngptlw,param%nlay), clwpmc(ngptlw,param%nlay), taucmc(ngptlw,param%nlay))

! Initializations

      oneminus = 1._r8 - 1.e-6_r8
      pi = 2._r8 * asin(1._r8)
      fluxfac = pi * 2.e4_r8                    ! orig:   fluxfac = pi * 2.d4  
      istart = 1
      iend = 16
      iout = 0
      ims = 1

! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

! *** This version uses McICA (imca = 1) ***

! Set icld to select of clear or cloud calculation and cloud overlap method  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap
! icld = 2, with clouds using maximum/random cloud overlap
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (param%icld.lt.0.or.param%icld.gt.3) param%icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 10, input total aerosol optical depth (tauaer) directly 
      iaer = 10

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 256 to 140 for input absorption coefficient 
! data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_lw_ini

!  This is the main longitude/column loop within RRTMG.
      do iplon = 1, param%ncol

!  Prepare atmospheric profile from GCM for use in RRTMG, and define
!  other input parameters.  

         call inatm (iplon, param%nlay, param%icld, iaer, &
              param%play, param%plev, param%tlay, param%tlev, param%tsfc, param%h2ovmr, &
              param%o3vmr, param%co2vmr, param%ch4vmr, param%o2vmr, param%n2ovmr, param%cfc11vmr, param%cfc12vmr, &
              param%cfc22vmr, param%ccl4vmr, param%emis, param%inflglw, param%iceflglw, param%liqflglw, &
              param%cldfmcl, param%taucmcl, param%ciwpmcl, param%clwpmcl, param%reicmcl, param%relqmcl, param%tauaer, &
              pavel, pz, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfmc, taucmc, ciwpmc, clwpmc, reicmc, dgesmc, relqmc, taua)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprop.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprop.  

         call cldprmc(param%nlay, inflag, iceflag, liqflag, cldfmc, ciwpmc, &
                      clwpmc, reicmc, dgesmc, relqmc, ncbands, taucmc)

! Calculate information needed by the radiative transfer routine
! that is specific to this atmosphere, especially some of the 
! coefficients and indices needed to compute the optical depths
! by interpolating data from stored reference atmospheres. 

         call setcoef(param%nlay, istart, pavel, tavel, tz, tbound, semiss, &
                      coldry, wkl, wbrodl, &
                      laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                      colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                      colbrd, fac00, fac01, fac10, fac11, &
                      rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                      rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                      rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                      selffac, selffrac, indself, forfac, forfrac, indfor, &
                      minorfrac, scaleminor, scaleminorn2, indminor)

!  Calculate the gaseous optical depths and Planck fractions for 
!  each longwave spectral band.

         call taumol(param%nlay, pavel, wx, coldry, &
                     laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                     colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                     colbrd, fac00, fac01, fac10, fac11, &
                     rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                     rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                     rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     minorfrac, scaleminor, scaleminorn2, indminor, &
                     fracs, taug)



! Combine gaseous and aerosol optical depths, if aerosol active
         if (iaer .eq. 0) then
            do k = 1, param%nlay
               do ig = 1, ngptlw 
                  taut(k,ig) = taug(k,ig)
               enddo
            enddo
         elseif (iaer .eq. 10) then
            do k = 1, param%nlay
               do ig = 1, ngptlw 
                  taut(k,ig) = taug(k,ig) + taua(k,ngb(ig))
               enddo
            enddo
         endif

! Call the radiative transfer routine.
! Either routine can be called to do clear sky calculation.  If clouds
! are present, then select routine based on cloud overlap assumption
! to be used.  Clear sky calculation is done simultaneously.
! For McICA, RTRNMC is called for clear and cloudy calculations.

         call rtrnmc(param%nlay, istart, iend, iout, pz, semiss, ncbands, &
                     cldfmc, taucmc, planklay, planklev, plankbnd, &
                     pwvcm, fracs, taut, &
                     totuflux, totdflux, fnet, htr, &
                     totuclfl, totdclfl, fnetc, htrc, totufluxs, totdfluxs )

!  Transfer up and down fluxes and heating rate to output arrays.
!  Vertical indexing goes from bottom to top

         do k = 0, param%nlay
            param%uflx(iplon,k+1) = totuflux(k)
            param%dflx(iplon,k+1) = totdflux(k)
            param%uflxc(iplon,k+1) = totuclfl(k)
            param%dflxc(iplon,k+1) = totdclfl(k)
            param%uflxs(:,iplon,k+1) = totufluxs(:,k)
            param%dflxs(:,iplon,k+1) = totdfluxs(:,k)
         enddo
         do k = 0, param%nlay-1
            param%hr(iplon,k+1) = htr(k)
            param%hrc(iplon,k+1) = htrc(k)
         enddo

      enddo

      end subroutine rrtmg_lw_mpe
      end module rrtmg_lw_rad_mpe

