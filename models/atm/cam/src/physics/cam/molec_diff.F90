module molec_diff

  !------------------------------------------------------------------------------------------------- !
  ! Module to compute molecular diffusivity for various constituents                                 !
  !                                                                                                  !
  ! Public interfaces :                                                                              !
  !                                                                                                  !
  !    init_molec_diff           Initializes time independent coefficients                           !
  !    init_timestep_molec_diff  Time-step initialization for molecular diffusivity                  !
  !    compute_molec_diff        Computes constituent-independent terms for moleculuar diffusivity   !
  !    vd_lu_qdecomp             Computes constituent-dependent terms for moleculuar diffusivity and !
  !                              updates terms in the triadiagonal matrix used for the implicit      !
  !                              solution of the diffusion equation                                  !
  !                                                                                                  !
  !---------------------------Code history---------------------------------------------------------- !
  ! Modularized     :  J. McCaa, September 2004                                                      !
  ! Lastly Arranged :  S. Park,  January.  2010                                                      !
  !                    M. Mills, November  2011
  !------------------------------------------------------------------------------------------------- !

  use perf_mod
  use physconst,    only : mbarv
  use constituents, only : pcnst
  use phys_control, only : waccmx_is             !WACCM-X runtime switch
  use ref_pres,     only : nbot_molec, ntop_molec

  implicit none
  private
  save

  public init_molec_diff
  public init_timestep_molec_diff
  public compute_molec_diff
  public vd_lu_qdecomp

  ! ---------- !
  ! Parameters !
  ! ---------- !

  integer,  parameter   :: r8 = selected_real_kind(12) ! 8 byte real

  real(r8), parameter   :: km_fac = 3.55E-7_r8         ! Molecular viscosity constant [ unit ? ]
  real(r8), parameter   :: pr_num = 1._r8              ! Prandtl number [ no unit ]
  real(r8), parameter   :: pwr    = 2._r8/3._r8        ! Exponentiation factor [ no unit ]
  real(r8), parameter   :: d0     = 1.52E20_r8         ! Diffusion factor [ m-1 s-1 ] molec sqrt(kg/kmol/K) [ unit ? ]
                                                       ! Aerononmy, Part B, Banks and Kockarts (1973), p39
                                                       ! Note text cites 1.52E18 cm-1 ...

  real(r8)              :: rair                        ! Gas constant for dry air
  real(r8)              :: mw_dry                      ! Molecular weight of dry air
  real(r8)              :: n_avog                      ! Avogadro's number [ molec/kmol ]
  real(r8)              :: gravit
  real(r8)              :: cpair
  real(r8)              :: kbtz                        ! Boltzman constant

  real(r8), allocatable :: mw_fac(:)                   ! sqrt(1/M_q + 1/M_d) in constituent diffusivity [  unit ? ]
  real(r8), allocatable :: alphath(:)                  ! Thermal diffusion factor, -0.38 for H, 0 for others

  logical :: waccmx_mode = .false.

contains

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_molec_diff( kind, ncnst, rair_in, mw_dry_in, n_avog_in, gravit_in, &
                              cpair_in, kbtz_in, errstring)

    use constituents,     only : cnst_mw, cnst_get_ind
    use upper_bc,         only : ubc_init
    use physics_buffer,   only : physics_buffer_desc

    integer,  intent(in)  :: kind           ! Kind of reals being passed in
    integer,  intent(in)  :: ncnst          ! Number of constituents
    real(r8), intent(in)  :: rair_in
    real(r8), intent(in)  :: mw_dry_in      ! Molecular weight of dry air
    real(r8), intent(in)  :: n_avog_in      ! Avogadro's number [ molec/kmol ]
    real(r8), intent(in)  :: gravit_in
    real(r8), intent(in)  :: cpair_in
    real(r8), intent(in)  :: kbtz_in        ! Boltzman constant

    character(len=*), intent(out) :: errstring

    ! Local

    integer               :: k              ! Level index
    integer               :: m              ! Constituent index
    integer               :: indx_H         ! Constituent index for H
    integer               :: ierr           ! Allocate error check

    errstring = ' '

    rair       = rair_in
    mw_dry     = mw_dry_in
    n_avog     = n_avog_in
    gravit     = gravit_in
    cpair      = cpair_in
    kbtz       = kbtz_in

    if( kind /= r8 ) then
       errstring = 'inconsistent KIND of reals passed to init_molec_diff'
       return
    end if

    ! Determine whether WACCM-X is on.
    waccmx_mode = waccmx_is('ionosphere') .or. waccmx_is('neutral')

  ! Initialize upper boundary condition variables

    call ubc_init()

  ! Molecular weight factor in constitutent diffusivity
  ! ***** FAKE THIS FOR NOW USING MOLECULAR WEIGHT OF DRY AIR FOR ALL TRACERS ****

    allocate(mw_fac(ncnst))
    do m = 1, ncnst
       mw_fac(m) = d0 * mw_dry * sqrt(1._r8/mw_dry + 1._r8/cnst_mw(m)) / n_avog
    end do

    !--------------------------------------------------------------------------------------------
    ! For WACCM-X, get H data index and initialize thermal diffusion coefficient
    !--------------------------------------------------------------------------------------------
    if ( waccmx_mode ) then

      call cnst_get_ind('H',  indx_H)

      allocate(alphath(ncnst), stat=ierr)
      if ( ierr /= 0 ) then
         errstring = 'allocate failed in init_molec_diff'
         return
      end if
      alphath(:ncnst) = 0._r8
      alphath(indx_H) = -0.38_r8

    endif

  end subroutine init_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_timestep_molec_diff(pbuf2d, state)
    !--------------------------- !
    ! Timestep dependent setting !
    !--------------------------- !
    use upper_bc,     only : ubc_timestep_init
    use physics_types,only: physics_state
    use ppgrid,       only: begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    call ubc_timestep_init( pbuf2d, state)

  end subroutine init_timestep_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  integer function compute_molec_diff( lchnk             ,                                          &
       pcols             , pver                , ncnst     , ncol     , t      , pmid   , pint   ,  &
       zi                , ztodt               , kvm       , kvt      , tint   , rhoi   , tmpi2  ,  &
       kq_scal           , ubc_t               , ubc_mmr   , ubc_flux , dse_top, cc_top ,           &
       cnst_mw_out       , cnst_fixed_ubc_out  , cnst_fixed_ubflx_out , mw_fac_out      ,           &
       ntop_molec_out    , nbot_molec_out      , kvt_returned )

    use upper_bc,        only : ubc_get_vals
    use constituents,    only : cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx
    use physconst,       only : cpairv, rairv, kmvis, kmcnd

    ! --------------------- !
    ! Input-Output Argument !
    ! --------------------- !

    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncnst
    integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
    integer,  intent(in)    :: lchnk                     ! Chunk identifier
    real(r8), intent(in)    :: t(pcols,pver)             ! Temperature input
    real(r8), intent(in)    :: pmid(pcols,pver)          ! Midpoint pressures
    real(r8), intent(in)    :: pint(pcols,pver+1)        ! Interface pressures
    real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface heights
    real(r8), intent(in)    :: ztodt                     ! 2 delta-t

    real(r8), intent(inout) :: kvm(pcols,pver+1)         ! Viscosity ( diffusivity for momentum )
    real(r8), intent(out)   :: kvt(pcols,pver+1)         ! Kinematic molecular conductivity
    real(r8), intent(inout) :: tint(pcols,pver+1)        ! Interface temperature
    real(r8), intent(inout) :: rhoi(pcols,pver+1)        ! Density ( rho ) at interfaces
    real(r8), intent(inout) :: tmpi2(pcols,pver+1)       ! dt*(g*rho)**2/dp at interfaces

    real(r8), intent(out)   :: kq_scal(pcols,pver+1)     ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8), intent(out)   :: ubc_t(pcols)              ! Upper boundary temperature (K)
    real(r8), intent(out)   :: ubc_mmr(pcols,ncnst)      ! Upper boundary mixing ratios [ kg/kg ]
    real(r8), intent(out)   :: ubc_flux(ncnst)           ! Upper boundary flux [ kg/s/m^2 ]
    real(r8), intent(out)   :: cnst_mw_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubc_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubflx_out(ncnst)
    real(r8), intent(out)   :: mw_fac_out(pcols,pver+1,ncnst) ! composition dependent mw_fac on interface level
    real(r8), intent(out)   :: dse_top(pcols)            ! dse on top boundary
    real(r8), intent(out)   :: cc_top(pcols)             ! Lower diagonal at top interface
    integer,  intent(out)   :: ntop_molec_out
    integer,  intent(out)   :: nbot_molec_out
    logical,  intent(out)   :: kvt_returned              ! Whether we actually returned kvt (vs. kvh).

    ! --------------- !
    ! Local variables !
    ! --------------- !

    integer                 :: m                          ! Constituent index
    integer                 :: i                          ! Column index
    integer                 :: k                          ! Level index

    real(r8)                :: mbarvi                     ! mbarv on interface level
    real(r8)                :: km_top(pcols)              ! molecular conductivity at the top

    real(r8)                :: mkvisc                     ! Molecular kinematic viscosity c*tint**(2/3)/rho

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! We don't apply cpairv to kvt if WACCM-X is on.
    kvt_returned = waccmx_mode

  ! Get upper boundary values

    call ubc_get_vals( lchnk, ncol, ntop_molec, pint, zi, ubc_t, ubc_mmr, ubc_flux )

  ! Below are already computed, just need to be copied for output

    cnst_mw_out(:ncnst)          = cnst_mw(:ncnst)
    cnst_fixed_ubc_out(:ncnst)   = cnst_fixed_ubc(:ncnst)
    cnst_fixed_ubflx_out(:ncnst) = cnst_fixed_ubflx(:ncnst)
    ntop_molec_out               = ntop_molec
    nbot_molec_out               = nbot_molec

    ! Zero out constituent ubc's that are not used.
    do m = 1, ncnst
       if (.not. cnst_fixed_ubc(m)) then
          ubc_mmr(:,m) = 0._r8
       end if
    end do

    !
    !  Need variable mw_fac for kvt and constant otherwise
    !
    if ( kvt_returned ) then
      do m = 1, ncnst
        do k = ntop_molec+1, nbot_molec
          do i = 1, ncol
             mbarvi = 0.5_r8 * (mbarv(i,k-1,lchnk)+mbarv(i,k,lchnk))
             mw_fac_out(i,k,m) = d0 * mbarvi * sqrt(1._r8/mbarvi + 1._r8/cnst_mw(m)) / n_avog
          enddo
        enddo
        mw_fac_out(:ncol,ntop_molec,m) = 1.5_r8*mw_fac_out(:ncol,ntop_molec+1,m)-.5_r8*mw_fac_out(:ncol,ntop_molec+2,m)
        do k = nbot_molec+1, pver+1
          mw_fac_out(:ncol,k,m) = mw_fac_out(:ncol,nbot_molec,m)
        enddo
      end do
    else
      do k = 1, pver+1
        do i = 1, ncol
          mw_fac_out(i,k,:ncnst) = mw_fac(:ncnst)
        enddo
      enddo
    endif

  ! Density and related factors for molecular diffusion and ubc.
  ! Always have a fixed upper boundary T if molecular diffusion is active. Why ?
  ! For kvt, set ubc temperature to average of next two lower interface level temperatures

    if ( kvt_returned ) then
      tint(:ncol,ntop_molec) = 1.5_r8*tint(:ncol,ntop_molec+1)-.5_r8*tint(:ncol,ntop_molec+2)
    else
      tint (:ncol,ntop_molec) = ubc_t(:ncol)
    endif

    rhoi (:ncol,ntop_molec) = pint(:ncol,ntop_molec) / ( rairv(:ncol,ntop_molec,lchnk) * tint(:ncol,ntop_molec) )
    tmpi2(:ncol,ntop_molec) = ztodt * ( gravit * rhoi(:ncol,ntop_molec))**2 &
                                    / ( pmid(:ncol,ntop_molec) - pint(:ncol,ntop_molec) )

  ! Compute molecular kinematic viscosity, heat diffusivity and factor for constituent diffusivity
  ! This is a key part of the code.  For WACCM-X, use constituent dependent molecular viscosity and conductivity

    kvt     = 0._r8
    kq_scal = 0._r8
    if ( kvt_returned ) then
      do k = ntop_molec, nbot_molec
         do i = 1, ncol
           mkvisc  = kmvis(i,k,lchnk) / rhoi(i,k)
           kvm(i,k) = kvm(i,k) + mkvisc
           mkvisc  = kmcnd(i,k,lchnk) / rhoi(i,k)
           kvt(i,k) = mkvisc
           kq_scal(i,k) = sqrt(tint(i,k)) / rhoi(i,k)
         end do
      end do
    else
      do k = ntop_molec, nbot_molec
        do i = 1, ncol
          mkvisc   = km_fac * tint(i,k)**pwr / rhoi(i,k)
          kvm(i,k) = kvm(i,k) + mkvisc
          kvt(i,k) = mkvisc * pr_num * cpairv(i,k,lchnk)
          kq_scal(i,k) = sqrt(tint(i,k)) / rhoi(i,k)
        end do
      end do
    endif

  ! Top boundary condition for dry static energy

    dse_top(:ncol) = cpairv(:ncol,ntop_molec,lchnk) * tint(:ncol,ntop_molec) + gravit * zi(:ncol,ntop_molec)

  ! Top value of cc for dry static energy

    if (kvt_returned) then
       do i = 1, ncol
          cc_top(i) = ztodt * gravit**2 * rhoi(i,ntop_molec) * km_fac * &
               ubc_t(i)**pwr / ( pmid(i,1) - pint(i,1) )
       enddo

    else
       cc_top = 0._r8
    end if

    compute_molec_diff = 1
    return
  end function compute_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  function vd_lu_qdecomp( &
       pcols , pver   , ncol       , fixed_ubc  , mw     , &
       kv    , kq_scal, mw_facm    , dpidz_sq   , p      , &
       interface_boundary, molec_boundary, rhoi   ,        &
       tint  , ztodt  , ntop_molec , nbot_molec , nbot   , &
       lchnk , t          , m      , no_molec_decomp)      result(decomp)

    use infnan, only: nan, assignment(=)
    use coords_1d, only: Coords1D
    use linear_1d_operators, only: BoundaryType, TriDiagDecomp
    use vdiff_lu_solver, only: fin_vol_lu_decomp

    !------------------------------------------------------------------------------ !
    ! Add the molecular diffusivity to the turbulent diffusivity for a consitutent. !
    ! Update the superdiagonal (ca(k)), diagonal (cb(k)) and subdiagonal (cc(k))    !
    ! coefficients of the tridiagonal diffusion matrix, also ze and denominator.    !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncol                  ! Number of atmospheric columns

    integer,  intent(in)    :: ntop_molec
    integer,  intent(in)    :: nbot_molec
    integer,  intent(in)    :: nbot

    logical,  intent(in)    :: fixed_ubc             ! Fixed upper boundary condition flag
    real(r8), intent(in)    :: kv(pcols,pver+1)      ! Eddy diffusivity
    real(r8), intent(in)    :: kq_scal(pcols,pver+1) ! Molecular diffusivity ( kq_fac*sqrt(T)*m_d/rho )
    real(r8), intent(in)    :: mw                    ! Molecular weight for this constituent
    real(r8), intent(in)    :: mw_facm(pcols,pver+1) ! composition dependent sqrt(1/M_q + 1/M_d) for this constituent
    real(r8), intent(in)    :: dpidz_sq(ncol,pver+1) ! (g*rho)**2 (square of vertical derivative of pint)
    type(Coords1D), intent(in) :: p                  ! Pressure coordinates
    type(BoundaryType), intent(in) :: interface_boundary ! Boundary on grid edge.
    type(BoundaryType), intent(in) :: molec_boundary ! Boundary at edge of molec_diff region.
    real(r8), intent(in)    :: rhoi(pcols,pver+1)    ! Density at interfaces [ kg/m3 ]
    real(r8), intent(in)    :: tint(pcols,pver+1)    ! Interface temperature [ K ]
    real(r8), intent(in)    :: ztodt                 ! 2 delta-t [ s ]

    integer,  intent(in)    :: lchnk		    ! Chunk number
    real(r8), intent(in)    :: t(pcols,pver)	    ! temperature
    integer,  intent(in)    :: m 		    ! cnst index

    ! Decomposition covering levels without vertical diffusion.
    type(TriDiagDecomp), intent(in) :: no_molec_decomp

    ! LU decomposition information for solver.
    type(TriDiagDecomp) :: decomp

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    ! Level index.
    integer :: k

    ! Molecular diffusivity for constituent.
    real(r8)                :: kmq(ncol,nbot_molec+1)

    ! Term for drift due to molecular separation: (m_i/m - 1) / p
    real(r8)                :: mw_term(ncol,nbot_molec+1)

    ! Diffusion coefficient.
    real(r8)                :: diff_coef(ncol,nbot_molec+1)
    ! Advection velocity.
    real(r8)                :: advect_v(ncol,nbot_molec+1)

    ! 1/mbar * d(mbar)/dp
    real(r8)                :: gradm(ncol,nbot_molec+1)

    ! alphaTh/T * dT/dp, for now alphaTh is non-zero only for H.
    real(r8)                :: gradt(ncol,nbot_molec+1)

    ! mbarv at interface
    real(r8)                :: mbarvi(ncol)

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! --------------------------------------------------------------------- !
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the !
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are !
    ! a combination of ca and cc; they are not required by the solver.      !
    !---------------------------------------------------------------------- !

    call t_startf('vd_lu_qdecomp')

    kmq  = 0._r8
    mw_term = 0._r8
    gradm = 0._r8
    gradt = 0._r8

    ! Compute difference between scale heights of constituent and dry air

    if ( waccmx_mode ) then

       ! Top level first.
       k = ntop_molec
       mbarvi = .75_r8*mbarv(:ncol,k,lchnk)+0.5_r8*mbarv(:ncol,k+1,lchnk) &
            -.25_r8*mbarv(:ncol,k+2,lchnk)
       mw_term(:,k) = (mw/mbarvi - 1._r8) / p%ifc(:,k)
       gradm(:,k) = (mbarv(:ncol,k,lchnk)-mbarvi)/ &
            (p%mid(:,k)-p%ifc(:,k))/ &
            (mbarv(:ncol,k,lchnk)+mbarvi)*2._r8

       if (alphath(m) /= 0._r8) then
          gradt(:,k) = alphath(m)*(t(:ncol,k)-tint(:ncol,k))/ &
               (p%mid(:ncol,k)-p%ifc(:ncol,k))/ &
               (t(:ncol,k)+tint(:ncol,k))*2._r8
       end if

       ! Interior of molecular diffusion region.
       do k = ntop_molec+1, nbot_molec
          mbarvi = 0.5_r8 * (mbarv(:ncol,k-1,lchnk)+mbarv(:ncol,k,lchnk))
          mw_term(:,k) = (mw/mbarvi - 1._r8) / p%ifc(:,k)
          gradm(:,k) = (mbarv(:ncol,k,lchnk)-mbarv(:ncol,k-1,lchnk)) * &
               p%rdst(:,k-1)/mbarvi
       enddo

       if (alphath(m) /= 0._r8) then
          do k = ntop_molec+1, nbot_molec
             gradt(:,k) = alphath(m)*(t(:ncol,k)-t(:ncol,k-1)) &
                  *p%rdst(:,k-1)/tint(:ncol,k)
          end do
       end if

       ! Leave nbot_molec+1 terms as zero, because molecular diffusion is
       ! small at the lower boundary.

    else

       do k = ntop_molec, nbot_molec
          mw_term(:,k) = (mw/mw_dry - 1._r8) / p%ifc(:ncol,k)
       enddo

    endif

    !-------------------- !
    ! Molecular diffusion !
    !-------------------- !

    ! Start with non-molecular portion of diffusion.

    ! Molecular diffusion coefficient.
    do k = ntop_molec, nbot_molec
       kmq(:,k)  = kq_scal(:ncol,k) * mw_facm(:ncol,k)
    end do

    diff_coef = kv(:ncol,:nbot_molec+1) + kmq

    ! "Drift" terms.
    advect_v = kmq*mw_term
    if ( waccmx_mode ) then
       advect_v = advect_v - kmq*gradt - &
            (kv(:ncol,:nbot_molec+1) + kmq)*gradm
    end if

    ! Convert from z to pressure representation.
    diff_coef = dpidz_sq(:,:nbot_molec+1) * diff_coef
    advect_v = dpidz_sq(:,:nbot_molec+1) * advect_v

    if( fixed_ubc ) then
       decomp = fin_vol_lu_decomp(ztodt, p, &
            coef_q_diff=diff_coef, coef_q_adv=advect_v, &
            upper_bndry=interface_boundary, &
            lower_bndry=molec_boundary, &
            graft_decomp=no_molec_decomp)
    else
       decomp = fin_vol_lu_decomp(ztodt, p, &
            coef_q_diff=diff_coef, coef_q_adv=advect_v, &
            lower_bndry=molec_boundary, &
            graft_decomp=no_molec_decomp)
    end if

    call t_stopf('vd_lu_qdecomp')

  end function vd_lu_qdecomp

end module molec_diff
