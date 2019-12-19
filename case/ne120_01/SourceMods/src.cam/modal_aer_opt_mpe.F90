module modal_aer_opt_mpe

! parameterizes aerosol coefficients using chebychev polynomial
! parameterize aerosol radiative properties in terms of
! surface mode wet radius and wet refractive index

! Ghan and Zaveri, JGR 2007.

! uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, pver, pverp
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use ref_pres,          only: top_lev => clim_modal_aero_top_lev
use physconst,         only: rhoh2o, rga, rair
use radconstants,      only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info, rad_cnst_get_aer_mmr, &
                             rad_cnst_get_aer_props, rad_cnst_get_mode_props
use physics_types,     only: physics_state

use physics_buffer, only : pbuf_get_index,physics_buffer_desc, pbuf_get_field
use cam_history,       only: phys_decomp, addfld, add_default, outfld
use cam_history_support, only: fillvalue

use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
use modal_aero_calcsize,    only: modal_aero_calcsize_diag

use modal_aer_opt, only: modal_size_parameters, &
                         xrmin, xrmax, &
                         mod_crefwsw => crefwsw, mod_crefwlw => crefwlw, &
                         dgnumwet_idx, qaerwat_idx

implicit none
#include <mpif.h>
private
save

public :: modal_aero_sw_parallel, modal_aero_lw_parallel, modal_aero_sw_param, modal_aero_lw_param


! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
! in terms of refractive index and wet radius
integer, parameter :: ncoef=5, prefr=7, prefi=10

character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                       '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

type modal_aero_sw_param
   integer :: list_idx
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes, imode
   integer :: nspec
   integer :: mpi_rank

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state), pointer :: state

! --- Summers ---

   real(r8), pointer :: ssavis(:)

   real(r8), pointer :: extinct(:,:)
   real(r8), pointer :: extinctnir(:,:)
   real(r8), pointer :: extinctuv(:,:)
   real(r8), pointer :: absorb(:,:)

   real(r8), pointer :: burden(:)
   real(r8), pointer :: burdendust(:), burdenso4(:), burdenbc(:), &
                        burdenpom(:), burdensoa(:), burdenseasalt(:)

   real(r8), pointer :: aodmode(:)
   real(r8), pointer :: dustaodmode(:)          ! dust aod in aerosol mode

   ! total species AOD
   real(r8), pointer :: dustaod(:), so4aod(:), bcaod(:), &
                        pomaod(:), soaaod(:), seasaltaod(:)

   real(r8), pointer :: aoduv(:)               ! extinction optical depth in uv
   real(r8), pointer :: aoduvst(:)             ! stratospheric extinction optical depth in uv
   real(r8), pointer :: aodnir(:)              ! extinction optical depth in nir
   real(r8), pointer :: aodnirst(:)            ! stratospheric extinction optical depth in nir

! --- Local ---
   integer, pointer :: troplev(:)

   real(r8), pointer :: mass(:,:)        ! layer mass
   real(r8), pointer :: air_density(:,:) ! (kg/m3)

   real(r8), pointer :: radsurf(:,:)    ! aerosol surface mode radius
   real(r8), pointer :: logradsurf(:,:) ! log(aerosol surface mode radius)
   real(r8), pointer :: cheb(:,:,:)

! --- Pointers ---

   real(r8), pointer :: dgnumwet(:,:)     ! number mode wet diameter
   real(r8), pointer :: qaerwat(:,:)      ! aerosol water (g/g)

   real(r8), pointer :: refrtabsw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitabsw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: extpsw(:,:,:,:) ! specific extinction
   real(r8), pointer :: abspsw(:,:,:,:) ! specific absorption
   real(r8), pointer :: asmpsw(:,:,:,:) ! asymmetry factor

   real(r8), pointer :: tauxar(:,:,:) ! layer extinction optical depth
   real(r8), pointer :: wa(:,:,:)     ! layer single-scatter albedo
   real(r8), pointer :: ga(:,:,:)     ! asymmetry factor
   real(r8), pointer :: fa(:,:,:)     ! forward scattered fraction

   integer(8) :: source_mmr_a_ptr
   integer(8) :: idx_mmr_a_ptr
   integer(8) :: idx_props_ptr
   integer(8) :: spectype_ptr
   integer(8) :: state_q_ptr
   integer :: state_q_dim1, state_q_dim2

   integer(8) :: end_type
end type

type modal_aero_lw_param
   integer :: list_idx
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes, imode
   integer :: nspec

   real(8) :: sigma_logr_aer

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state), pointer :: state

! --- Local ---

! --- Pointers ---
   real(r8), pointer :: mass(:,:) ! layer mass

   real(r8), pointer :: dgnumwet(:,:)  ! wet number mode diameter (m)
   real(r8), pointer :: qaerwat(:,:)   ! aerosol water (g/g)

   real(r8), pointer :: refrtablw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitablw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: absplw(:,:,:,:) ! specific absorption

   real(r8), pointer :: tauxar(:,:,:) ! layer extinction optical depth

   integer(8) :: source_mmr_a_ptr
   integer(8) :: idx_mmr_a_ptr
   integer(8) :: idx_props_ptr
   integer(8) :: state_q_ptr
   integer :: state_q_dim1, state_q_dim2

   integer(8) :: end_type
end type

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================

subroutine modal_aero_sw_parallel(list_idx, state, pbuf, nnite, idxnite, &
                         tauxar, wa, ga, fa)

   ! calculates aerosol sw radiative properties
   use rad_constituents, only : modes, ma_list
   use phys_prop, only : physprop, numphysprops
   use tropopause,     only : tropopause_find

   integer,             intent(in) :: list_idx       ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state          ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer,             intent(in) :: nnite          ! number of night columns
   integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns

   real(r8), target, intent(out) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
   real(r8), target, intent(out) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
   real(r8), target, intent(out) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
   real(r8), target, intent(out) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

   ! Local variables
   integer :: m_idx
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec
   integer, target :: troplev(pcols)

   real(r8), target :: mass(pcols,pver)        ! layer mass
   real(r8), target :: air_density(pcols,pver) ! (kg/m3)

   real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
   real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
   real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
   real(r8), pointer :: wetdens_m(:,:,:)  ! 

   real(r8) :: sigma_logr_aer         ! geometric standard deviation of number distribution
   real(r8), target :: radsurf(pcols,pver)    ! aerosol surface mode radius
   real(r8), target :: logradsurf(pcols,pver) ! log(aerosol surface mode radius)
   real(r8), target :: cheb(ncoef,pcols,pver)

   real(r8), pointer :: refrtabsw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitabsw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: extpsw(:,:,:,:) ! specific extinction
   real(r8), pointer :: abspsw(:,:,:,:) ! specific absorption
   real(r8), pointer :: asmpsw(:,:,:,:) ! asymmetry factor

   ! Diagnostics
   real(r8), target :: extinct(pcols,pver)
   real(r8), target :: extinctnir(pcols,pver)
   real(r8), target :: extinctuv(pcols,pver)
   real(r8), target :: absorb(pcols,pver)
   real(r8), target :: aodvis(pcols)               ! extinction optical depth
   real(r8), target :: aodvisst(pcols)             ! stratospheric extinction optical depth
   real(r8), target :: aodabs(pcols)               ! absorption optical depth

   real(r8), target :: aodabsbc(pcols)             ! absorption optical depth of BC

   real(r8), target :: ssavis(pcols)

   real(r8), target :: burden(pcols)
   real(r8), target :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
                       burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols)

   real(r8), target :: aodmode(pcols)
   real(r8), target :: dustaodmode(pcols)          ! dust aod in aerosol mode

   ! total species AOD
   real(r8), target :: dustaod(pcols), so4aod(pcols), bcaod(pcols), &
                       pomaod(pcols), soaaod(pcols), seasaltaod(pcols)

   real(r8), target :: aoduv(pcols)               ! extinction optical depth in uv
   real(r8), target :: aoduvst(pcols)             ! stratospheric extinction optical depth in uv
   real(r8), target :: aodnir(pcols)              ! extinction optical depth in nir
   real(r8), target :: aodnirst(pcols)            ! stratospheric extinction optical depth in nir


   character(len=32) :: outname

   ! debug output
   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf            ! volume fraction of insoluble aerosol
   character(len=*), parameter :: subname = 'modal_aero_sw'

   type(modal_aero_sw_param) :: pm
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8
   wa(:ncol,:,:)     = 0._r8
   ga(:ncol,:,:)     = 0._r8
   fa(:ncol,:,:)     = 0._r8

   ! zero'th layer does not contain aerosol
   tauxar(1:ncol,0,:)  = 0._r8
   wa(1:ncol,0,:)      = 0.925_r8
   ga(1:ncol,0,:)      = 0.850_r8
   fa(1:ncol,0,:)      = 0.7225_r8

   mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
   air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

   ! diagnostics for visible band summed over modes
   extinct(1:ncol,:)     = 0.0_r8
   absorb(1:ncol,:)      = 0.0_r8
   aodvis(1:ncol)        = 0.0_r8
   aodvisst(1:ncol)      = 0.0_r8
   aodabs(1:ncol)        = 0.0_r8
   burdendust(:ncol)     = 0.0_r8
   burdenso4(:ncol)      = 0.0_r8
   burdenpom(:ncol)      = 0.0_r8
   burdensoa(:ncol)      = 0.0_r8
   burdenbc(:ncol)       = 0.0_r8
   burdenseasalt(:ncol)  = 0.0_r8
   ssavis(1:ncol)        = 0.0_r8

   aodabsbc(:ncol)       = 0.0_r8 
   dustaod(:ncol)        = 0.0_r8
   so4aod(:ncol)         = 0.0_r8
   pomaod(:ncol)         = 0.0_r8
   soaaod(:ncol)         = 0.0_r8
   bcaod(:ncol)          = 0.0_r8
   seasaltaod(:ncol)     = 0.0_r8

   ! diags for other bands
   extinctuv(1:ncol,:)   = 0.0_r8
   extinctnir(1:ncol,:)  = 0.0_r8
   aoduv(:ncol)          = 0.0_r8
   aodnir(:ncol)         = 0.0_r8
   aoduvst(:ncol)        = 0.0_r8
   aodnirst(:ncol)       = 0.0_r8
   call tropopause_find(state, troplev)

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx == 0) then
      ! water uptake and wet radius for the climate list has already been calculated
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
   else
      ! If doing a diagnostic calculation then need to calculate the wet radius
      ! and water uptake for the diagnostic modes
      call modal_aero_calcsize_diag(state, pbuf, list_idx, dgnumdry_m)  
      call modal_aero_wateruptake_dr(state, pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m)
   endif

   pm%list_idx = list_idx
   pm%lchnk = lchnk
   pm%ncol = ncol
   pm%nmodes = nmodes

   pm%tauxar => tauxar
   pm%wa => wa
   pm%ga => ga
   pm%fa => fa

   pm%troplev => troplev
   pm%mass => mass
   pm%air_density => air_density

   pm%radsurf => radsurf
   pm%logradsurf => logradsurf
   pm%cheb => cheb

   pm%pbuf => pbuf
   pm%state => state

   pm%ssavis => ssavis
   pm%extinct => extinct
   pm%extinctnir => extinctnir
   pm%extinctuv => extinctuv
   pm%absorb => absorb
   pm%burden => burden
   pm%burdendust => burdendust
   pm%burdenso4 => burdenso4
   pm%burdenbc => burdenbc
   pm%burdenpom => burdenpom
   pm%burdensoa => burdensoa
   pm%burdenseasalt => burdenseasalt
   pm%aodmode => aodmode
   pm%dustaodmode => dustaodmode
   pm%dustaod => dustaod
   pm%so4aod => so4aod
   pm%bcaod => bcaod
   pm%pomaod => pomaod
   pm%soaaod => soaaod
   pm%seasaltaod => seasaltaod
   pm%aoduv => aoduv
   pm%aoduvst => aoduvst
   pm%aodnir => aodnir
   pm%aodnirst => aodnirst

   call MPI_Comm_Rank(MPI_COMM_WORLD, pm%mpi_rank, 0)

   do m = 1, nmodes

      ! diagnostics for visible band for each mode
      burden(:ncol)       = 0._r8
      aodmode(1:ncol)     = 0.0_r8
      dustaodmode(1:ncol) = 0.0_r8

      pm%dgnumwet => dgnumwet_m(:,:,m)
      pm%qaerwat  => qaerwat_m(:,:,m)

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtabsw=refrtabsw , &
         refitabsw=refitabsw, extpsw=extpsw, abspsw=abspsw, asmpsw=asmpsw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      call modal_size_parameters(ncol, sigma_logr_aer, pm%dgnumwet, radsurf, logradsurf, cheb)

      pm%refrtabsw => refrtabsw
      pm%refitabsw => refitabsw
      pm%extpsw => extpsw
      pm%abspsw => abspsw
      pm%asmpsw => asmpsw
      pm%imode = m
      pm%nspec = nspec

      m_idx = ma_list(list_idx)%idx(m)
      pm%source_mmr_a_ptr = loc(modes%comps(m_idx)%source_mmr_a)
      pm%idx_mmr_a_ptr = loc(modes%comps(m_idx)%idx_mmr_a)
      pm%idx_props_ptr = loc(modes%comps(m_idx)%idx_props)
      pm%spectype_ptr = loc(modes%comps(m_idx)%type)
      pm%state_q_ptr = loc(state%q)
      pm%state_q_dim1 = size(state%q, dim=1)
      pm%state_q_dim2 = size(state%q, dim=2)

      call athread_spawn_modal_aero_sw_kern(pm)
      call athread_join

      ! mode diagnostics
      ! The diagnostics are currently only output for the climate list.  Code mods will
      ! be necessary to provide output for the rad_diag lists.
      if (list_idx == 0) then
         do i = 1, nnite
            burden(idxnite(i))  = fillvalue
            aodmode(idxnite(i)) = fillvalue
            dustaodmode(idxnite(i)) = fillvalue
         end do

         write(outname,'(a,i1)') 'BURDEN', m
         call outfld(trim(outname), burden, pcols, lchnk)

         write(outname,'(a,i1)') 'AODMODE', m
         call outfld(trim(outname), aodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'AODDUST', m
         call outfld(trim(outname), dustaodmode, pcols, lchnk)

      end if

   end do ! nmodes

   if (list_idx > 0) then
      deallocate(dgnumdry_m)
      deallocate(dgnumwet_m)
      deallocate(qaerwat_m)
      deallocate(wetdens_m)
   end if

   ! Output visible band diagnostics for quantities summed over the modes
   ! These fields are put out for diagnostic lists as well as the climate list.
   do i = 1, nnite
      extinct(idxnite(i),:) = fillvalue
      absorb(idxnite(i),:)  = fillvalue
      aodvis(idxnite(i))    = fillvalue
      aodabs(idxnite(i))    = fillvalue
      aodvisst(idxnite(i))  = fillvalue
   end do

   call outfld('EXTINCT'//diag(list_idx),  extinct, pcols, lchnk)
   call outfld('ABSORB'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVIS'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODABS'//diag(list_idx),   aodabs,  pcols, lchnk)
   call outfld('AODVISst'//diag(list_idx), aodvisst,pcols, lchnk)

   ! These diagnostics are output only for climate list
   if (list_idx == 0) then
      do i = 1, ncol
         if (aodvis(i) > 1.e-10_r8) then
            ssavis(i) = ssavis(i)/aodvis(i)
         else
            ssavis(i) = 0.925_r8
         endif
      end do

      do i = 1, nnite
         ssavis(idxnite(i))     = fillvalue

         aoduv(idxnite(i))      = fillvalue
         aodnir(idxnite(i))     = fillvalue
         aoduvst(idxnite(i))    = fillvalue
         aodnirst(idxnite(i))   = fillvalue
         extinctuv(idxnite(i),:)  = fillvalue
         extinctnir(idxnite(i),:) = fillvalue

         burdendust(idxnite(i)) = fillvalue
         burdenso4(idxnite(i))  = fillvalue
         burdenpom(idxnite(i))  = fillvalue
         burdensoa(idxnite(i))  = fillvalue
         burdenbc(idxnite(i))   = fillvalue
         burdenseasalt(idxnite(i)) = fillvalue

         aodabsbc(idxnite(i))   = fillvalue

         dustaod(idxnite(i))    = fillvalue
         so4aod(idxnite(i))     = fillvalue
         pomaod(idxnite(i))     = fillvalue
         soaaod(idxnite(i))     = fillvalue
         bcaod(idxnite(i))      = fillvalue
         seasaltaod(idxnite(i)) = fillvalue
       end do

      call outfld('SSAVIS',        ssavis,        pcols, lchnk)

      call outfld('EXTINCTUV',     extinctuv,     pcols, lchnk)
      call outfld('EXTINCTNIR',    extinctnir,    pcols, lchnk)
      call outfld('AODUV',         aoduv,         pcols, lchnk)
      call outfld('AODNIR',        aodnir,        pcols, lchnk)
      call outfld('AODUVst',       aoduvst,       pcols, lchnk)
      call outfld('AODNIRst',      aodnirst,      pcols, lchnk)

      call outfld('BURDENDUST',    burdendust,    pcols, lchnk)
      call outfld('BURDENSO4' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOM' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOA' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBC'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALT', burdenseasalt, pcols, lchnk)

      call outfld('AODABSBC',      aodabsbc,      pcols, lchnk)

      call outfld('AODDUST',       dustaod,       pcols, lchnk)
      call outfld('AODSO4',        so4aod,        pcols, lchnk)
      call outfld('AODPOM',        pomaod,        pcols, lchnk)
      call outfld('AODSOA',        soaaod,        pcols, lchnk)
      call outfld('AODBC',         bcaod,         pcols, lchnk)
      call outfld('AODSS',         seasaltaod,    pcols, lchnk)
   end if

end subroutine modal_aero_sw_parallel

subroutine modal_aero_sw_kern (pm)
   use rad_constituents, only : modes, ma_list
   use phys_prop, only : physprop

   type(modal_aero_sw_param), intent(inout) :: pm

   integer :: list_idx

   ! Local variables
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec
   integer :: troplev(pcols)

   real(r8) :: mass(pcols)        ! layer mass
   real(r8) :: air_density(pcols) ! (kg/m3)

   real(r8),    pointer :: gl_specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specmmr(pcols)         ! species mass mixing ratio
   real(r8)             :: specdens               ! species density (kg/m3)
   complex(r8), pointer :: gl_specrefindex(:)     ! species refractive index
   complex(r8)          :: specrefindex(nswbands) ! species refractive index
   character*32         :: spectype            ! species type
   real(r8)             :: hygro_aer           !

   real(r8) :: radsurf(pcols)    ! aerosol surface mode radius
   real(r8) :: logradsurf(pcols) ! log(aerosol surface mode radius)
   real(r8) :: cheb(ncoef,pcols)

   real(r8)    :: refr(pcols)     ! real part of refractive index
   real(r8)    :: refi(pcols)     ! imaginary part of refractive index
   complex(r8) :: crefin(pcols)   ! complex refractive index

   real(r8) :: vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: watervol(pcols) ! volume concentration of water in each mode (m3/kg)
   real(r8) :: wetvol(pcols)   ! volume concentration of wet mode (m3/kg)

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cext(pcols,ncoef), cabs(pcols,ncoef), casm(pcols,ncoef)
   real(r8) :: pext(pcols)     ! parameterized specific extinction (m2/kg)
   real(r8) :: specpext(pcols) ! specific extinction (m2/kg)
   real(r8) :: dopaer(pcols)   ! aerosol optical depth in layer
   real(r8) :: pabs(pcols)     ! parameterized specific absorption (m2/kg)
   real(r8) :: pasm(pcols)     ! parameterized asymmetry factor
   real(r8) :: palb(pcols)     ! parameterized single scattering albedo

   ! Diagnostics
   real(r8) :: extinct(pcols)
   real(r8) :: extinctnir(pcols)
   real(r8) :: extinctuv(pcols)
   real(r8) :: absorb(pcols)
   real(r8) :: aodvis(pcols)               ! extinction optical depth
   real(r8) :: aodvisst(pcols)             ! stratospheric extinction optical depth
   real(r8) :: aodabs(pcols)               ! absorption optical depth

   real(r8) :: aodabsbc(pcols)             ! absorption optical depth of BC

   real(r8) :: ssavis(pcols)
   real(r8) :: dustvol(pcols)              ! volume concentration of dust in aerosol mode (m3/kg)

   real(r8) :: burden(pcols)
   real(r8) :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
               burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols)

   real(r8) :: aodmode(pcols)
   real(r8) :: dustaodmode(pcols)          ! dust aod in aerosol mode

   real(r8) :: specrefr, specrefi
   real(r8) :: scatdust(pcols), scatso4(pcols), scatbc(pcols), &
               scatpom(pcols), scatsoa(pcols), scatseasalt(pcols)
   real(r8) :: absdust(pcols), absso4(pcols), absbc(pcols), &
               abspom(pcols), abssoa(pcols), absseasalt(pcols)
   real(r8) :: hygrodust(pcols), hygroso4(pcols), hygrobc(pcols), &
               hygropom(pcols), hygrosoa(pcols), hygroseasalt(pcols)

   real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro
   real(r8) :: aodc                        ! aod of component

   ! total species AOD
   real(r8) :: dustaod(pcols), so4aod(pcols), bcaod(pcols), &
               pomaod(pcols), soaaod(pcols), seasaltaod(pcols)

   logical :: savaervis ! true if visible wavelength (0.55 micron)
   logical :: savaernir ! true if near ir wavelength (~0.88 micron)
   logical :: savaeruv  ! true if uv wavelength (~0.35 micron)

   real(r8) :: aoduv(pcols)               ! extinction optical depth in uv
   real(r8) :: aoduvst(pcols)             ! stratospheric extinction optical depth in uv
   real(r8) :: aodnir(pcols)              ! extinction optical depth in nir
   real(r8) :: aodnirst(pcols)            ! stratospheric extinction optical depth in nir

! --- Output ---
   real(r8) :: tauxar(pcols)
   real(r8) :: wa(pcols)
   real(r8) :: ga(pcols)
   real(r8) :: fa(pcols)

! --- Local ---
   complex(r8) :: crefwsw(nswbands)
   real(r8) :: qaerwat(pcols)             ! species mass mixing ratio
   real(r8) :: refrtabsw(prefr)           ! table of real refractive indices for aerosols
   real(r8) :: refitabsw(prefi)           ! table of imag refractive indices for aerosols
   real(r8) :: extpsw(ncoef,prefr,prefi)  ! specific extinction
   real(r8) :: abspsw(ncoef,prefr,prefi)  ! specific absorption
   real(r8) :: asmpsw(ncoef,prefr,prefi)  ! asymmetry factor

! --- RAD_CNST ---
   integer :: m_idx, idx, id
   character(len=1) :: source

   list_idx = pm%list_idx
   ncol = pm%ncol
   nmodes = pm%nmodes
   m = pm%imode
   nspec = pm%nspec

   crefwsw = mod_crefwsw
   troplev = pm%troplev

   ssavis = pm%ssavis
   burden = pm%burden
   burdendust = pm%burdendust
   burdenso4 = pm%burdenso4
   burdenbc = pm%burdenbc
   burdenpom = pm%burdenpom
   burdensoa = pm%burdensoa
   burdenseasalt = pm%burdenseasalt
   aodmode = pm%aodmode
   dustaodmode = pm%dustaodmode
   dustaod = pm%dustaod
   so4aod = pm%so4aod
   bcaod = pm%bcaod
   pomaod = pm%pomaod
   soaaod = pm%soaaod
   seasaltaod = pm%seasaltaod
   aoduv = pm%aoduv
   aoduvst = pm%aoduvst
   aodnir = pm%aodnir
   aodnirst = pm%aodnirst

      do isw = 1, nswbands
         refrtabsw = pm%refrtabsw(:,isw)
         refitabsw = pm%refitabsw(:,isw)
         extpsw = pm%extpsw(:,:,:,isw)
         abspsw = pm%abspsw(:,:,:,isw)
         asmpsw = pm%asmpsw(:,:,:,isw)

         savaervis = (isw .eq. idx_sw_diag)
         savaeruv  = (isw .eq. idx_uv_diag)
         savaernir = (isw .eq. idx_nir_diag)

         do k = top_lev, pver
            mass = pm%mass(:,k)
            air_density = pm%air_density(:,k)
            radsurf = pm%radsurf(:,k)
            logradsurf = pm%logradsurf(:,k)
            cheb = pm%cheb(:,:,k)
            extinct = pm%extinct(:,k)
            extinctnir = pm%extinctnir(:,k)
            extinctuv = pm%extinctuv(:,k)
            absorb = pm%absorb(:,k)
            qaerwat = pm%qaerwat(:,k)


            ! form bulk refractive index
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8
            dustvol(:ncol) = 0._r8

            scatdust(:ncol)     = 0._r8
            absdust(:ncol)      = 0._r8
            hygrodust(:ncol)    = 0._r8
            scatso4(:ncol)      = 0._r8
            absso4(:ncol)       = 0._r8
            hygroso4(:ncol)     = 0._r8
            scatbc(:ncol)       = 0._r8
            absbc(:ncol)        = 0._r8
            hygrobc(:ncol)      = 0._r8
            scatpom(:ncol)      = 0._r8
            abspom(:ncol)       = 0._r8
            hygropom(:ncol)     = 0._r8
            scatsoa(:ncol)      = 0._r8
            abssoa(:ncol)       = 0._r8
            hygrosoa(:ncol)     = 0._r8
            scatseasalt(:ncol)  = 0._r8
            absseasalt(:ncol)   = 0._r8
            hygroseasalt(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec

               m_idx = ma_list(list_idx)%idx(m)

               source = modes%comps(m_idx)%source_mmr_a(l)
               idx    = modes%comps(m_idx)%idx_mmr_a(l)
               if (source /= 'A') call abort
               gl_specmmr => pm%state%q(:,:,idx)

               id = modes%comps(m_idx)%idx_props(l)
               specdens = physprop(id)%density_aer
               gl_specrefindex => physprop(id)%refindex_aer_sw
               spectype = modes%comps(m_idx)%type(l)
               hygro_aer = physprop(id)%hygro_aer

               specmmr = gl_specmmr(:,k)
               specrefindex = gl_specrefindex(:)

               do i = 1, ncol
                  vol(i)      = specmmr(i)/specdens
                  dryvol(i)   = dryvol(i) + vol(i)
                  crefin(i)   = crefin(i) + vol(i)*specrefindex(isw)
               end do

               ! compute some diagnostics for visible band only
               if (savaervis) then

                  specrefr = real(specrefindex(isw))
                  specrefi = aimag(specrefindex(isw))

                  do i = 1, ncol
                     burden(i) = burden(i) + specmmr(i)*mass(k)
                  end do

                  if (trim(spectype) == 'dust') then
                     do i = 1, ncol
                        burdendust(i) = burdendust(i) + specmmr(i)*mass(k)
                        dustvol(i)    = vol(i)
                        scatdust(i)   = vol(i)*specrefr
                        absdust(i)    = -vol(i)*specrefi
                        hygrodust(i)  = vol(i)*hygro_aer
                     end do
                  end if

                  if (trim(spectype) == 'sulfate') then
                     do i = 1, ncol
                        burdenso4(i) = burdenso4(i) + specmmr(i)*mass(k)
                        scatso4(i)   = vol(i)*specrefr
                        absso4(i)    = -vol(i)*specrefi
                        hygroso4(i)  = vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'black-c') then
                     do i = 1, ncol
                        burdenbc(i) = burdenbc(i) + specmmr(i)*mass(k)
                        scatbc(i)   = vol(i)*specrefr
                        absbc(i)    = -vol(i)*specrefi
                        hygrobc(i)  = vol(i)*hygro_aer
                   end do
                  end if
                  if (trim(spectype) == 'p-organic') then
                     do i = 1, ncol
                        burdenpom(i) = burdenpom(i) + specmmr(i)*mass(k)
                        scatpom(i)   = vol(i)*specrefr
                        abspom(i)    = -vol(i)*specrefi
                        hygropom(i)  = vol(i)*hygro_aer
                      end do
                  end if
                  if (trim(spectype) == 's-organic') then
                     do i = 1, ncol
                        burdensoa(i) = burdensoa(i) + specmmr(i)*mass(k)
                        scatsoa(i)   = vol(i)*specrefr
                        abssoa(i)    = -vol(i)*specrefi
                        hygrosoa(i)  = vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'seasalt') then
                     do i = 1, ncol
                        burdenseasalt(i) = burdenseasalt(i) + specmmr(i)*mass(k)
                        scatseasalt(i)   = vol(i)*specrefr
                        absseasalt(i)    = -vol(i)*specrefi
                        hygroseasalt(i)  = vol(i)*hygro_aer
                      end do
                  end if

               end if
            end do ! species loop

            do i = 1, ncol
               watervol(i) = qaerwat(i)/rhoh2o
               wetvol(i) = watervol(i) + dryvol(i)
               if (watervol(i) < 0._r8) then
                  watervol(i) = 0._r8
                  wetvol(i) = dryvol(i)
               end if

               ! volume mixing
               crefin(i) = crefin(i) + watervol(i)*crefwsw(isw)
               crefin(i) = crefin(i)/max(wetvol(i),1.e-60_r8)
               refr(i)   = real(crefin(i))
               refi(i)   = abs(aimag(crefin(i)))
            end do

            ! call t_startf('binterp')

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(extpsw, ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw, refitabsw, &
                         itab, jtab, ttab, utab, cext)
            call binterp(abspsw, ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw, refitabsw, &
                         itab, jtab, ttab, utab, cabs)
            call binterp(asmpsw, ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw, refitabsw, &
                         itab, jtab, ttab, utab, casm)

            ! call t_stopf('binterp')

            ! parameterized optical properties
            do i=1,ncol

               if (logradsurf(i) .le. xrmax) then
                  pext(i) = 0.5_r8*cext(i,1)
                  do nc = 2, ncoef
                     pext(i) = pext(i) + cheb(nc,i)*cext(i,nc)
                  enddo
                  pext(i) = exp(pext(i))
               else
                  pext(i) = 1.5_r8/(radsurf(i)*rhoh2o) ! geometric optics
               endif

               ! convert from m2/kg water to m2/kg aerosol
               specpext(i) = pext(i)
               pext(i) = pext(i)*wetvol(i)*rhoh2o
               pabs(i) = 0.5_r8*cabs(i,1)
               pasm(i) = 0.5_r8*casm(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheb(nc,i)*cabs(i,nc)
                  pasm(i) = pasm(i) + cheb(nc,i)*casm(i,nc)
               enddo
               pabs(i) = pabs(i)*wetvol(i)*rhoh2o
               pabs(i) = max(0._r8,pabs(i))
               pabs(i) = min(pext(i),pabs(i))

               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)

               dopaer(i) = pext(i)*mass(i)
            end do

            if (savaeruv) then
               do i = 1, ncol
                 extinctuv(i) = extinctuv(i) + dopaer(i)*air_density(i)/mass(i)
                 aoduv(i) = aoduv(i) + dopaer(i)
                  if (k.le.troplev(i)) then
                    aoduvst(i) = aoduvst(i) + dopaer(i)
                  end if
               end do
            end if

            if (savaernir) then
               do i = 1, ncol
                  extinctnir(i) = extinctnir(i) + dopaer(i)*air_density(i)/mass(i)
                  aodnir(i) = aodnir(i) + dopaer(i)
                  if (k.le.troplev(i)) then
                    aodnirst(i) = aodnirst(i) + dopaer(i)
                  end if
               end do
            endif

            ! Save aerosol optical depth at longest visible wavelength
            ! sum over layers
            if (savaervis) then
               ! aerosol extinction (/m)
               do i = 1, ncol
                  extinct(i) = extinct(i) + dopaer(i)*air_density(i)/mass(i)
                  absorb(i)  = absorb(i) + pabs(i)*air_density(i)
                  aodvis(i)    = aodvis(i) + dopaer(i)
                  aodabs(i)    = aodabs(i) + pabs(i)*mass(i)
                  aodmode(i)   = aodmode(i) + dopaer(i)
                  ssavis(i)    = ssavis(i) + dopaer(i)*palb(i)
                  if (k.le.troplev(i)) then
                    aodvisst(i) = aodvisst(i) + dopaer(i)
                  end if

                  if (wetvol(i) > 1.e-40_r8) then

                     dustaodmode(i) = dustaodmode(i) + dopaer(i)*dustvol(i)/wetvol(i)

                     ! partition optical depth into contributions from each constituent
                     ! assume contribution is proportional to refractive index X volume

                     scath2o        = watervol(i)*real(crefwsw(isw))
                     absh2o         = -watervol(i)*aimag(crefwsw(isw))
                     sumscat        = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                                      scatdust(i) + scatseasalt(i) + scath2o
                     sumabs         = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                                      absdust(i) + absseasalt(i) + absh2o
                     sumhygro       = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                                      hygrodust(i) + hygroseasalt(i)

                     scatdust(i)    = (scatdust(i) + scath2o*hygrodust(i)/sumhygro)/sumscat
                     absdust(i)     = (absdust(i) + absh2o*hygrodust(i)/sumhygro)/sumabs

                     scatso4(i)     = (scatso4(i) + scath2o*hygroso4(i)/sumhygro)/sumscat
                     absso4(i)      = (absso4(i) + absh2o*hygroso4(i)/sumhygro)/sumabs

                     scatpom(i)     = (scatpom(i) + scath2o*hygropom(i)/sumhygro)/sumscat
                     abspom(i)      = (abspom(i) + absh2o*hygropom(i)/sumhygro)/sumabs

                     scatsoa(i)     = (scatsoa(i) + scath2o*hygrosoa(i)/sumhygro)/sumscat
                     abssoa(i)      = (abssoa(i) + absh2o*hygrosoa(i)/sumhygro)/sumabs

                     scatbc(i)      = (scatbc(i) + scath2o*hygrobc(i)/sumhygro)/sumscat
                     absbc(i)       = (absbc(i) + absh2o*hygrobc(i)/sumhygro)/sumabs

                     scatseasalt(i) = (scatseasalt(i) + scath2o*hygroseasalt(i)/sumhygro)/sumscat
                     absseasalt(i)  = (absseasalt(i) + absh2o*hygroseasalt(i)/sumhygro)/sumabs
                     
                     aodabsbc(i)    = aodabsbc(i) + absbc(i)*dopaer(i)*(1.0_r8-palb(i))

                     aodc           = (absdust(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust(i))*dopaer(i)
                     dustaod(i)     = dustaod(i) + aodc

                     aodc           = (absso4(i)*(1.0_r8 - palb(i)) + palb(i)*scatso4(i))*dopaer(i)
                     so4aod(i)      = so4aod(i) + aodc

                     aodc           = (abspom(i)*(1.0_r8 - palb(i)) + palb(i)*scatpom(i))*dopaer(i)
                     pomaod(i)      = pomaod(i) + aodc

                     aodc           = (abssoa(i)*(1.0_r8 - palb(i)) + palb(i)*scatsoa(i))*dopaer(i)
                     soaaod(i)      = soaaod(i) + aodc

                     aodc           = (absbc(i)*(1.0_r8 - palb(i)) + palb(i)*scatbc(i))*dopaer(i)
                     bcaod(i)       = bcaod(i) + aodc

                     aodc           = (absseasalt(i)*(1.0_r8 - palb(i)) + palb(i)*scatseasalt(i))*dopaer(i)
                     seasaltaod(i)  = seasaltaod(i) + aodc

                  endif

               end do
            endif

            tauxar = pm%tauxar(:,k,isw)
            wa = pm%wa(:,k,isw)
            ga = pm%ga(:,k,isw)
            fa = pm%fa(:,k,isw)
            do i=1,ncol
               tauxar(i) = tauxar(i) + dopaer(i)
               wa(i)     = wa(i)     + dopaer(i)*palb(i)
               ga(i)     = ga(i)     + dopaer(i)*palb(i)*pasm(i)
               fa(i)     = fa(i)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
            end do
            pm%tauxar(:,k,isw) = tauxar
            pm%wa(:,k,isw) = wa
            pm%ga(:,k,isw) = ga
            pm%fa(:,k,isw) = fa

            pm%extinct(:,k) = extinct
            pm%extinctnir(:,k) = extinctnir
            pm%extinctuv(:,k) = extinctuv
            pm%absorb(:,k) = absorb

         end do ! pver

      end do ! sw bands

   pm%ssavis = ssavis
   pm%burden = burden
   pm%burdendust = burdendust
   pm%burdenso4 = burdenso4
   pm%burdenbc = burdenbc
   pm%burdenpom = burdenpom
   pm%burdensoa = burdensoa
   pm%burdenseasalt = burdenseasalt
   pm%aodmode = aodmode
   pm%dustaodmode = dustaodmode
   pm%dustaod = dustaod
   pm%so4aod = so4aod
   pm%bcaod = bcaod
   pm%pomaod = pomaod
   pm%soaaod = soaaod
   pm%seasaltaod = seasaltaod
   pm%aoduv = aoduv
   pm%aoduvst = aoduvst
   pm%aodnir = aodnir
   pm%aodnirst = aodnirst

end subroutine modal_aero_sw_kern

!===============================================================================

subroutine modal_aero_lw_parallel(list_idx, state, pbuf, tauxar)
   use rad_constituents, only : modes, ma_list
   use phys_prop, only : physprop, numphysprops

   ! calculates aerosol lw radiative properties

   integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state    ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), target, intent(out) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

   ! Local variables
   integer :: m_idx
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec

   real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
   real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
   real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
   real(r8), pointer :: wetdens_m(:,:,:)  ! 

   real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution

   real(r8), target :: mass(pcols,pver) ! layer mass

   real(r8), pointer :: refrtablw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitablw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: absplw(:,:,:,:) ! specific absorption

   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf             ! volume fraction of insoluble aerosol

   character(len=*), parameter :: subname = 'modal_aero_lw'

   type(modal_aero_lw_param) :: pm
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8

   ! dry mass in each cell
   mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx == 0) then
      ! water uptake and wet radius for the climate list has already been calculated
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
   else
      ! If doing a diagnostic calculation then need to calculate the wet radius
      ! and water uptake for the diagnostic modes
      call modal_aero_calcsize_diag(state, pbuf, list_idx, dgnumdry_m)  
      call modal_aero_wateruptake_dr(state, pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m)
   endif

   pm%list_idx = list_idx
   pm%lchnk = lchnk
   pm%ncol = ncol
   pm%nmodes = nmodes

   pm%tauxar => tauxar

   pm%mass => mass

   pm%pbuf => pbuf
   pm%state => state

   do m = 1, nmodes

      pm%dgnumwet => dgnumwet_m(:,:,m)
      pm%qaerwat  => qaerwat_m(:,:,m)

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtablw=refrtablw , &
         refitablw=refitablw, absplw=absplw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      ! this is the same calculation that's done in modal_size_parameters, but there
      ! some intermediate results are saved and the chebyshev polynomials are stored
      ! in a array with different index order.  Could be unified.

      pm%refrtablw => refrtablw
      pm%refitablw => refitablw
      pm%absplw => absplw
      pm%imode = m
      pm%nspec = nspec
      pm%sigma_logr_aer = sigma_logr_aer

      m_idx = ma_list(list_idx)%idx(m)
      pm%source_mmr_a_ptr = loc(modes%comps(m_idx)%source_mmr_a)
      pm%idx_mmr_a_ptr = loc(modes%comps(m_idx)%idx_mmr_a)
      pm%idx_props_ptr = loc(modes%comps(m_idx)%idx_props)
      pm%state_q_ptr = loc(state%q)
      pm%state_q_dim1 = size(state%q, dim=1)
      pm%state_q_dim2 = size(state%q, dim=2)

      call athread_spawn_modal_aero_lw_kern(pm)
      call athread_join

   end do ! m = 1, nmodes

   if (list_idx > 0) then
      deallocate(dgnumdry_m)
      deallocate(dgnumwet_m)
      deallocate(qaerwat_m)
      deallocate(wetdens_m)
   end if

end subroutine modal_aero_lw_parallel

subroutine modal_aero_lw_kern (pm)
   use rad_constituents, only : modes, ma_list
   use phys_prop, only : physprop

   type(modal_aero_lw_param), intent(inout) :: pm

   integer :: list_idx
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec

   real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution

   real(r8) :: alnsg_amode             ! log of geometric standard deviation of number distribution
   real(r8) :: xrad(pcols)
   real(r8) :: cheby(ncoef,pcols,pver)  ! chebychef polynomials

   real(r8) :: mass(pcols) ! layer mass

   real(r8),    pointer :: gl_specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specmmr(pcols)         ! species mass mixing ratio
   real(r8)             :: specdens               ! species density (kg/m3)
   complex(r8), pointer :: gl_specrefindex(:)     ! species refractive index
   complex(r8)          :: specrefindex(nlwbands) ! species refractive index

   real(r8) :: vol(pcols)       ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: wetvol(pcols)    ! volume concentration of wet mode (m3/kg)
   real(r8) :: watervol(pcols)  ! volume concentration of water in each mode (m3/kg)
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index
   complex(r8) :: crefin(pcols) ! complex refractive index

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cabs(pcols,ncoef)
   real(r8) :: pabs(pcols)      ! parameterized specific absorption (m2/kg)
   real(r8) :: dopaer(pcols)    ! aerosol optical depth in layer

! --- Local ---
   complex(r8) :: crefwlw(nlwbands)       ! species mass mixing ratio
   real(r8) :: tauxar(pcols)              ! species mass mixing ratio
   real(r8) :: qaerwat(pcols)             ! species mass mixing ratio
   real(r8) :: dgnumwet(pcols)            ! species mass mixing ratio
   real(r8) :: refrtablw(prefr)           ! table of real refractive indices for aerosols
   real(r8) :: refitablw(prefi)           ! table of imag refractive indices for aerosols
   real(r8) :: absplw(ncoef,prefr,prefi)  ! specific absorption

! --- RAD_CNST ---
   integer :: m_idx, idx, id
   character(len=1) :: source

   list_idx = pm%list_idx
   ncol = pm%ncol
   nmodes = pm%nmodes
   m = pm%imode
   nspec = pm%nspec
   sigma_logr_aer = pm%sigma_logr_aer
   crefwlw = mod_crefwlw

      do k = top_lev, pver
         dgnumwet = pm%dgnumwet(:,k)
         do i = 1, ncol
            alnsg_amode = log( sigma_logr_aer )
            ! convert from number diameter to surface area
            xrad(i) = log(0.5_r8*dgnumwet(i)) + 2.0_r8*alnsg_amode*alnsg_amode
            ! normalize size parameter
            xrad(i) = max(xrad(i), xrmin)
            xrad(i) = min(xrad(i), xrmax)
            xrad(i) = (2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
            ! chebyshev polynomials
            cheby(1,i,k) = 1.0_r8
            cheby(2,i,k) = xrad(i)
            do nc = 3, ncoef
               cheby(nc,i,k) = 2.0_r8*xrad(i)*cheby(nc-1,i,k)-cheby(nc-2,i,k)
            end do
         end do
      end do

      do ilw = 1, nlwbands
         refrtablw = pm%refrtablw(:,ilw)
         refitablw = pm%refitablw(:,ilw)
         absplw = pm%absplw(:,:,:,ilw)

         do k = top_lev, pver
            mass = pm%mass(:,k)
            qaerwat = pm%qaerwat(:,k)

            ! form bulk refractive index. Use volume mixing for infrared
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               m_idx = ma_list(list_idx)%idx(m)

               source = modes%comps(m_idx)%source_mmr_a(l)
               idx    = modes%comps(m_idx)%idx_mmr_a(l)
               if (source /= 'A') call abort
               gl_specmmr => pm%state%q(:,:,idx)

               id = modes%comps(m_idx)%idx_props(l)
               specdens = physprop(id)%density_aer
               gl_specrefindex => physprop(id)%refindex_aer_lw

               specmmr = gl_specmmr(:,k)
               specrefindex = gl_specrefindex(:)

               do i = 1, ncol
                  vol(i)    = specmmr(i)/specdens
                  dryvol(i) = dryvol(i) + vol(i)
                  crefin(i) = crefin(i) + vol(i)*specrefindex(ilw)
               end do
            end do

            do i = 1, ncol
               watervol(i) = qaerwat(i)/rhoh2o
               wetvol(i)   = watervol(i) + dryvol(i)
               if (watervol(i) < 0.0_r8) then
                  watervol(i) = 0._r8
                  wetvol(i)   = dryvol(i)
               end if

               crefin(i) = crefin(i) + watervol(i)*crefwlw(ilw)
               if (wetvol(i) > 1.e-40_r8) crefin(i) = crefin(i)/wetvol(i)
               refr(i) = real(crefin(i))
               refi(i) = aimag(crefin(i))
            end do

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(absplw, ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtablw, refitablw, &
                         itab, jtab, ttab, utab, cabs)

            ! parameterized optical properties
            do i = 1, ncol
               pabs(i) = 0.5_r8*cabs(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheby(nc,i,k)*cabs(i,nc)
               end do
               pabs(i)   = pabs(i)*wetvol(i)*rhoh2o
               pabs(i)   = max(0._r8,pabs(i))
               dopaer(i) = pabs(i)*mass(i)
            end do

            tauxar = pm%tauxar(:,k,ilw)
            do i = 1, ncol
               tauxar(i) = tauxar(i) + dopaer(i)
            end do
            pm%tauxar(:,k,ilw) = tauxar

         end do ! k = top_lev, pver

      end do  ! nlwbands
end subroutine modal_aero_lw_kern

!===============================================================================

      subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

!     bilinear interpolation of table
!
      implicit none
      integer im,jm,km,ncol
      real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
      integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic
      real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols), &
             tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)

      if(ix(1).gt.0)go to 30
      if(im.gt.1)then
        do ic=1,ncol
          do i=1,im
            if(x(ic).lt.xtab(i))go to 10
          enddo
   10     ix(ic)=max0(i-1,1)
          ip1=min(ix(ic)+1,im)
          dx=(xtab(ip1)-xtab(ix(ic)))
          if(abs(dx).gt.1.e-20_r8)then
             t(ic)=(x(ic)-xtab(ix(ic)))/dx
          else
             t(ic)=0._r8
          endif
	end do
      else
        ix(:ncol)=1
        t(:ncol)=0._r8
      endif
      if(jm.gt.1)then
        do ic=1,ncol
          do j=1,jm
            if(y(ic).lt.ytab(j))go to 20
          enddo
   20     jy(ic)=max0(j-1,1)
          jp1=min(jy(ic)+1,jm)
          dy=(ytab(jp1)-ytab(jy(ic)))
          if(abs(dy).gt.1.e-20_r8)then
             u(ic)=(y(ic)-ytab(jy(ic)))/dy
          else
            u(ic)=0._r8
          endif
	end do
      else
        jy(:ncol)=1
        u(:ncol)=0._r8
      endif
   30 continue
      do ic=1,ncol
         tu(ic)=t(ic)*u(ic)
         tuc(ic)=t(ic)-tu(ic)
         tcuc(ic)=1._r8-tuc(ic)-u(ic)
         tcu(ic)=u(ic)-tu(ic)
         jp1=min(jy(ic)+1,jm)
         ip1=min(ix(ic)+1,im)
         do k=1,km
            out(ic,k)=tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))   &
               +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1)
	 end do
      enddo
      return
      end subroutine binterp

end module modal_aer_opt_mpe
