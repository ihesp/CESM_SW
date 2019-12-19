module modal_aer_opt_cpe

! parameterizes aerosol coefficients using chebychev polynomial
! parameterize aerosol radiative properties in terms of
! surface mode wet radius and wet refractive index

! Ghan and Zaveri, JGR 2007.

! uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, mod_pver => pver
use spmd_utils,        only: masterproc
use ref_pres,          only: clim_modal_aero_top_lev
use physconst,         only: mod_rhoh2o => rhoh2o
use radconstants,      only: mod_nswbands => nswbands, mod_nlwbands => nlwbands, mod_idx_sw_diag => idx_sw_diag, mod_idx_uv_diag => idx_uv_diag, mod_idx_nir_diag => idx_nir_diag
use rad_constituents,  only: n_diag

use modal_aer_opt, only : modal_size_parameters, &
                          mod_xrmin => xrmin, mod_xrmax => xrmax, &
                          mod_crefwsw => crefwsw, mod_crefwlw => crefwlw, &
                          dgnumwet_idx, qaerwat_idx
use modal_aer_opt_mpe, only : modal_aero_sw_param, modal_aero_lw_param

implicit none
private
save

! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
! in terms of refractive index and wet radius
integer, parameter :: cnst_ncoef=5, cnst_prefr=7, cnst_prefi=10

character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                       '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)
!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================

subroutine modal_aero_sw_kern (gl_pm)
   use phys_prop, only : physprop, numphysprops

   type(modal_aero_sw_param), intent(inout) :: gl_pm

   integer :: list_idx

   ! Local variables
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec
   integer :: troplev(pcols)

   real(r8) :: mass(pcols)        ! layer mass
   real(r8) :: air_density(pcols) ! (kg/m3)

   real(r8)             :: gl_specmmr(1)          ! species mass mixing ratio
   pointer(gl_specmmr_ptr, gl_specmmr)
   real(r8)             :: specmmr(pcols)         ! species mass mixing ratio
   real(r8)             :: specdens               ! species density (kg/m3)
   complex(r8) :: gl_specrefindex(mod_nswbands)     ! species refractive index
   pointer(gl_specrefindex_ptr, gl_specrefindex)
   complex(r8)          :: specrefindex(mod_nswbands) ! species refractive index
   character*32         :: spectype            ! species type
   real(r8)             :: hygro_aer           !

   real(r8) :: radsurf(pcols)    ! aerosol surface mode radius
   real(r8) :: logradsurf(pcols) ! log(aerosol surface mode radius)
   real(r8) :: cheb(cnst_ncoef,pcols)

   real(r8)    :: refr(pcols)     ! real part of refractive index
   real(r8)    :: refi(pcols)     ! imaginary part of refractive index
   complex(r8) :: crefin(pcols)   ! complex refractive index

   real(r8) :: vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: watervol(pcols) ! volume concentration of water in each mode (m3/kg)
   real(r8) :: wetvol(pcols)   ! volume concentration of wet mode (m3/kg)

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cext(pcols,cnst_ncoef), cabs(pcols,cnst_ncoef), casm(pcols,cnst_ncoef)
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
   complex(r8) :: crefwsw(mod_nswbands)
   real(r8) :: pm_qaerwat(pcols,mod_pver)
   pointer(pm_qaerwat_ptr, pm_qaerwat)
   real(r8) :: qaerwat(pcols)                            ! species mass mixing ratio
   real(r8) :: refrtabsw(cnst_prefr)                     ! table of real refractive indices for aerosols
   real(r8) :: refitabsw(cnst_prefi)                     ! table of imag refractive indices for aerosols
   real(r8) :: extpsw(cnst_ncoef,cnst_prefr,cnst_prefi)  ! specific extinction
   real(r8) :: abspsw(cnst_ncoef,cnst_prefr,cnst_prefi)  ! specific absorption
   real(r8) :: asmpsw(cnst_ncoef,cnst_prefr,cnst_prefi)  ! asymmetry factor

! --- Math data ---
   real(8) :: exp_data(528)

! --- RAD_CNST ---
   integer :: m_idx, idx, id
   character(len=1) :: source

   character*1 :: source_mmr_a_buf(gl_pm%nspec)
   integer :: idx_mmr_a_buf(gl_pm%nspec)
   integer :: idx_props_buf(gl_pm%nspec)
   character*32 :: spectype_buf(gl_pm%nspec)

   character*1 :: gl_source_mmr_a(1)
   pointer(gl_source_mmr_a_ptr, gl_source_mmr_a)
   integer :: gl_idx_mmr_a(1), gl_idx_props_a(1)
   pointer(gl_idx_mmr_a_ptr, gl_idx_mmr_a)
   pointer(gl_idx_props_ptr, gl_idx_props_a)
   character*32 :: gl_spectype_a(1)
   pointer(gl_spectype_ptr, gl_spectype_a)

   character*1 :: source_mmr_a(1)
   pointer(source_mmr_a_ptr, source_mmr_a)
   integer :: idx_mmr_a(1), idx_props_a(1)
   pointer(idx_mmr_a_ptr, idx_mmr_a)
   pointer(idx_props_ptr, idx_props_a)
   character*32 :: spectype_a(1)
   pointer(spectype_ptr, spectype_a)

   real(r8) :: density_aer_buf(numphysprops)
   real(r8) :: hygro_aer_buf(numphysprops)
   integer(8) :: refindex_aer_sw_buf(numphysprops)
   real(r8) :: density_aer_a(1), hygro_aer_a(1)
   pointer(density_aer_ptr, density_aer_a)
   pointer(hygro_aer_ptr, hygro_aer_a)
   integer(8) :: refindex_aer_sw_a(1)
   pointer(refindex_aer_sw_ptr, refindex_aer_sw_a)

! --- Slave ---
   integer :: myid, rid, cid
   integer :: cnst
   type(modal_aero_sw_param) :: pm

   integer :: pver
   integer :: top_lev
   real(r8) :: rhoh2o
   integer :: nswbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
   integer :: ncoef, prefr, prefi
   real(r8) :: xrmax

   character*1 :: s_A
   character*4 :: s_dust
   character*7 :: s_sulfate
   character*7 :: s_black_c
   character*9 :: s_p_organic
   character*9 :: s_s_organic
   character*7 :: s_seasalt

interface
   integer function trim_strcmp(a, b)
      character(len=*) :: a, b
   end function
end interface

   call set_exp_data(exp_data)
   ! call set_func_ptr
   call athread_get_id(myid)
   rid = myid / 8
   cid = mod(myid, 8)

   call dmemcpy(pm, gl_pm, loc(gl_pm%end_type) - loc(gl_pm))

   s_A = 'A'
   s_dust = 'dust'
   s_sulfate = 'sulfate'
   s_black_c = 'black-c'
   s_p_organic = 'p-organic'
   s_s_organic = 's-organic'
   s_seasalt = 'seasalt'

   pm_qaerwat_ptr = loc(pm%qaerwat)

   pver = mod_pver
   top_lev = clim_modal_aero_top_lev
   rhoh2o = mod_rhoh2o
   xrmax = mod_xrmax
   nswbands = mod_nswbands
   idx_sw_diag = mod_idx_sw_diag
   idx_uv_diag = mod_idx_uv_diag
   idx_nir_diag = mod_idx_nir_diag

   ncoef = cnst_ncoef
   prefr = cnst_prefr
   prefi = cnst_prefi

   list_idx = pm%list_idx
   ncol = pm%ncol
   nmodes = pm%nmodes
   m = pm%imode
   nspec = pm%nspec

   gl_source_mmr_a_ptr = pm%source_mmr_a_ptr
   gl_idx_mmr_a_ptr = pm%idx_mmr_a_ptr
   gl_idx_props_ptr = pm%idx_props_ptr
   gl_spectype_ptr = pm%spectype_ptr
   source_mmr_a_buf = gl_source_mmr_a(1:nspec)
   call dmemcpy(idx_mmr_a_buf, gl_idx_mmr_a, nspec * kind(idx_mmr_a_buf))
   call dmemcpy(idx_props_buf, gl_idx_props_a, nspec * kind(idx_props_buf))
   call dmemcpy(spectype_buf, gl_spectype_a, nspec * len(spectype_buf))
   source_mmr_a_ptr = loc(source_mmr_a_buf)
   idx_mmr_a_ptr    = loc(idx_mmr_a_buf)
   idx_props_ptr    = loc(idx_props_buf)
   spectype_ptr     = loc(spectype_buf)

   do i = 1, numphysprops
      density_aer_buf(i) = physprop(i)%density_aer
      hygro_aer_buf(i) = physprop(i)%hygro_aer
      refindex_aer_sw_buf(i) = loc(physprop(i)%refindex_aer_sw)
   enddo
   density_aer_ptr = loc(density_aer_buf)
   hygro_aer_ptr = loc(hygro_aer_buf)
   refindex_aer_sw_ptr = loc(refindex_aer_sw_buf)

   call dmemcpy(crefwsw, mod_crefwsw, nswbands * 2*kind(crefwsw))
   call dmemcpy(troplev, pm%troplev, ncol * kind(troplev))

#ifdef DIAGNOSTICS
   call dmemcpy(ssavis, pm%ssavis, ncol * kind(ssavis))
   call dmemcpy(burden, pm%burden, ncol * kind(burden))
   call dmemcpy(burdendust, pm%burdendust, ncol * kind(burdendust))
   call dmemcpy(burdenso4, pm%burdenso4, ncol * kind(burdenso4))
   call dmemcpy(burdenbc, pm%burdenbc, ncol * kind(burdenbc))
   call dmemcpy(burdenpom, pm%burdenpom, ncol * kind(burdenpom))
   call dmemcpy(burdensoa, pm%burdensoa, ncol * kind(burdensoa))
   call dmemcpy(burdenseasalt, pm%burdenseasalt, ncol * kind(burdenseasalt))
   call dmemcpy(aodmode, pm%aodmode, ncol * kind(aodmode))
   call dmemcpy(dustaodmode, pm%dustaodmode, ncol * kind(dustaodmode))
   call dmemcpy(dustaod, pm%dustaod, ncol * kind(dustaod))
   call dmemcpy(so4aod, pm%so4aod, ncol * kind(so4aod))
   call dmemcpy(bcaod, pm%bcaod, ncol * kind(bcaod))
   call dmemcpy(pomaod, pm%pomaod, ncol * kind(pomaod))
   call dmemcpy(soaaod, pm%soaaod, ncol * kind(soaaod))
   call dmemcpy(seasaltaod, pm%seasaltaod, ncol * kind(seasaltaod))
   call dmemcpy(aoduv, pm%aoduv, ncol * kind(aoduv))
   call dmemcpy(aoduvst, pm%aoduvst, ncol * kind(aoduvst))
   call dmemcpy(aodnir, pm%aodnir, ncol * kind(aodnir))
   call dmemcpy(aodnirst, pm%aodnirst, ncol * kind(aodnirst))
#endif

      do isw = cid + 1, nswbands, 8
         call dmemcpy(refrtabsw, pm%refrtabsw(:,isw), prefr * kind(refrtabsw))
         call dmemcpy(refitabsw, pm%refitabsw(:,isw), prefi * kind(refitabsw))
         call dmemcpy(extpsw, pm%extpsw(:,:,:,isw), ncoef*prefr*prefi * kind(extpsw))
         call dmemcpy(abspsw, pm%abspsw(:,:,:,isw), ncoef*prefr*prefi * kind(abspsw))
         call dmemcpy(asmpsw, pm%asmpsw(:,:,:,isw), ncoef*prefr*prefi * kind(asmpsw))

         savaervis = (isw .eq. idx_sw_diag)
         savaeruv  = (isw .eq. idx_uv_diag)
         savaernir = (isw .eq. idx_nir_diag)

         do k = rid + 1, pver, 8
            if (k < top_lev) cycle
            call dmemcpy(mass, pm%mass(:,k), ncol * kind(mass))
            call dmemcpy(air_density, pm%air_density(:,k), ncol * kind(air_density))
            call dmemcpy(radsurf, pm%radsurf(:,k), ncol * kind(radsurf))
            call dmemcpy(logradsurf, pm%logradsurf(:,k), ncol * kind(logradsurf))
            call dmemcpy(cheb, pm%cheb(:,:,k), ncoef*pcols * kind(cheb))
#ifdef DIAGNOSTICS
            call dmemcpy(extinct, pm%extinct(:,k), ncol * kind(extinct))
            call dmemcpy(extinctnir, pm%extinctnir(:,k), ncol * kind(extinctnir))
            call dmemcpy(extinctuv, pm%extinctuv(:,k), ncol * kind(extinctuv))
            call dmemcpy(absorb, pm%absorb(:,k), ncol * kind(absorb))
#endif
            call dmemcpy(qaerwat, pm_qaerwat(:,k), ncol * kind(qaerwat))

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

               source = source_mmr_a(l)
               idx    = idx_mmr_a(l)

               if (trim_strcmp(source, s_A) /= 0) call abort
               gl_specmmr_ptr = pm%state_q_ptr + ((idx-1) * pm%state_q_dim2 + (k-1)) * pm%state_q_dim1 * kind(specmmr)

               id = idx_props_a(l)
               specdens = density_aer_a(id)
               gl_specrefindex_ptr = refindex_aer_sw_a(id)
               cnst = len(spectype)
               call dmemcpy(spectype, spectype_a(l), cnst)
               hygro_aer = hygro_aer_a(id)

               call dmemcpy(specmmr, gl_specmmr, ncol * kind(specmmr))
               call dmemcpy(specrefindex, gl_specrefindex(:), nswbands * 2*kind(specrefindex))
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

                  if (trim_strcmp(spectype, s_dust) == 0) then
                     do i = 1, ncol
                        burdendust(i) = burdendust(i) + specmmr(i)*mass(k)
                        dustvol(i)    = vol(i)
                        scatdust(i)   = vol(i)*specrefr
                        absdust(i)    = -vol(i)*specrefi
                        hygrodust(i)  = vol(i)*hygro_aer
                     end do
                  end if

                  if (trim_strcmp(spectype, s_sulfate) == 0) then
                     do i = 1, ncol
                        burdenso4(i) = burdenso4(i) + specmmr(i)*mass(k)
                        scatso4(i)   = vol(i)*specrefr
                        absso4(i)    = -vol(i)*specrefi
                        hygroso4(i)  = vol(i)*hygro_aer
                     end do
                  end if
                  if (trim_strcmp(spectype, s_black_c) == 0) then
                     do i = 1, ncol
                        burdenbc(i) = burdenbc(i) + specmmr(i)*mass(k)
                        scatbc(i)   = vol(i)*specrefr
                        absbc(i)    = -vol(i)*specrefi
                        hygrobc(i)  = vol(i)*hygro_aer
                   end do
                  end if
                  if (trim_strcmp(spectype, s_p_organic) == 0) then
                     do i = 1, ncol
                        burdenpom(i) = burdenpom(i) + specmmr(i)*mass(k)
                        scatpom(i)   = vol(i)*specrefr
                        abspom(i)    = -vol(i)*specrefi
                        hygropom(i)  = vol(i)*hygro_aer
                      end do
                  end if
                  if (trim_strcmp(spectype, s_s_organic) == 0) then
                     do i = 1, ncol
                        burdensoa(i) = burdensoa(i) + specmmr(i)*mass(k)
                        scatsoa(i)   = vol(i)*specrefr
                        abssoa(i)    = -vol(i)*specrefi
                        hygrosoa(i)  = vol(i)*hygro_aer
                     end do
                  end if
                  if (trim_strcmp(spectype, s_seasalt) == 0) then
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

            call dmemcpy(tauxar, pm%tauxar(:,k,isw), ncol * kind(tauxar))
            call dmemcpy(wa, pm%wa(:,k,isw), ncol * kind(wa))
            call dmemcpy(ga, pm%ga(:,k,isw), ncol * kind(ga))
            call dmemcpy(fa, pm%fa(:,k,isw), ncol * kind(fa))
            do i=1,ncol
               tauxar(i) = tauxar(i) + dopaer(i)
               wa(i)     = wa(i)     + dopaer(i)*palb(i)
               ga(i)     = ga(i)     + dopaer(i)*palb(i)*pasm(i)
               fa(i)     = fa(i)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
            end do
            call dmemcpy(pm%tauxar(:,k,isw), tauxar, ncol * kind(tauxar))
            call dmemcpy(pm%wa(:,k,isw), wa, ncol * kind(wa))
            call dmemcpy(pm%ga(:,k,isw), ga, ncol * kind(ga))
            call dmemcpy(pm%fa(:,k,isw), fa, ncol * kind(fa))

#ifdef DIAGNOSTICS
            call dmemcpy(pm%extinct(:,k), extinct, ncol * kind(extinct))
            call dmemcpy(pm%extinctnir(:,k), extinctnir, ncol * kind(extinctnir))
            call dmemcpy(pm%extinctuv(:,k), extinctuv, ncol * kind(extinctuv))
            call dmemcpy(pm%absorb(:,k), absorb, ncol * kind(absorb))
#endif
         end do ! pver

      end do ! sw bands

#ifdef DIAGNOSTICS
   call dmemcpy(pm%ssavis, ssavis, ncol * kind(ssavis))
   call dmemcpy(pm%burden, burden, ncol * kind(burden))
   call dmemcpy(pm%burdendust, burdendust, ncol * kind(burdendust))
   call dmemcpy(pm%burdenso4, burdenso4, ncol * kind(burdenso4))
   call dmemcpy(pm%burdenbc, burdenbc, ncol * kind(burdenbc))
   call dmemcpy(pm%burdenpom, burdenpom, ncol * kind(burdenpom))
   call dmemcpy(pm%burdensoa, burdensoa, ncol * kind(burdensoa))
   call dmemcpy(pm%burdenseasalt, burdenseasalt, ncol * kind(burdenseasalt))
   call dmemcpy(pm%aodmode, aodmode, ncol * kind(aodmode))
   call dmemcpy(pm%dustaodmode, dustaodmode, ncol * kind(dustaodmode))
   call dmemcpy(pm%dustaod, dustaod, ncol * kind(dustaod))
   call dmemcpy(pm%so4aod, so4aod, ncol * kind(so4aod))
   call dmemcpy(pm%bcaod, bcaod, ncol * kind(bcaod))
   call dmemcpy(pm%pomaod, pomaod, ncol * kind(pomaod))
   call dmemcpy(pm%soaaod, soaaod, ncol * kind(soaaod))
   call dmemcpy(pm%seasaltaod, seasaltaod, ncol * kind(seasaltaod))
   call dmemcpy(pm%aoduv, aoduv, ncol * kind(aoduv))
   call dmemcpy(pm%aoduvst, aoduvst, ncol * kind(aoduvst))
   call dmemcpy(pm%aodnir, aodnir, ncol * kind(aodnir))
   call dmemcpy(pm%aodnirst, aodnirst, ncol * kind(aodnirst))
#endif
   call reset_exp_data

end subroutine modal_aero_sw_kern

subroutine modal_aero_lw_kern (gl_pm)
   use rad_constituents, only : modes, ma_list
   use phys_prop, only : physprop, numphysprops

   type(modal_aero_lw_param), intent(inout) :: gl_pm

   integer :: list_idx
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec

   real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution

   real(r8) :: alnsg_amode             ! log of geometric standard deviation of number distribution
   real(r8) :: xrad(pcols)
   real(r8) :: cheby(cnst_ncoef,pcols,mod_pver)  ! chebychef polynomials

   real(r8) :: mass(pcols) ! layer mass

   real(r8)             :: gl_specmmr(1)          ! species mass mixing ratio
   pointer(gl_specmmr_ptr, gl_specmmr)
   real(r8)             :: specmmr(pcols)         ! species mass mixing ratio
   real(r8)             :: specdens               ! species density (kg/m3)
   complex(r8) :: gl_specrefindex(mod_nlwbands)     ! species refractive index
   pointer(gl_specrefindex_ptr, gl_specrefindex)
   complex(r8)          :: specrefindex(mod_nlwbands) ! species refractive index
   real(r8)             :: hygro_aer           !

   real(r8) :: vol(pcols)       ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: wetvol(pcols)    ! volume concentration of wet mode (m3/kg)
   real(r8) :: watervol(pcols)  ! volume concentration of water in each mode (m3/kg)
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index
   complex(r8) :: crefin(pcols) ! complex refractive index

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cabs(pcols,cnst_ncoef)
   real(r8) :: pabs(pcols)      ! parameterized specific absorption (m2/kg)
   real(r8) :: dopaer(pcols)    ! aerosol optical depth in layer

! --- Local ---
   real(r8) :: pm_dgnumwet(pcols,mod_pver)
   pointer(pm_dgnumwet_ptr, pm_dgnumwet)
   real(r8) :: pm_qaerwat(pcols,mod_pver)
   pointer(pm_qaerwat_ptr, pm_qaerwat)
   complex(r8) :: crefwlw(mod_nlwbands)                  ! species mass mixing ratio
   real(r8) :: tauxar(pcols)                             ! species mass mixing ratio
   real(r8) :: qaerwat(pcols)                            ! species mass mixing ratio
   real(r8) :: dgnumwet(pcols)                           ! species mass mixing ratio
   real(r8) :: refrtablw(cnst_prefr)                     ! table of real refractive indices for aerosols
   real(r8) :: refitablw(cnst_prefi)                     ! table of imag refractive indices for aerosols
   real(r8) :: absplw(cnst_ncoef,cnst_prefr,cnst_prefi)  ! specific absorption

! --- RAD_CNST ---
   integer :: m_idx, idx, id
   character(len=1) :: source

   character*1 :: source_mmr_a_buf(gl_pm%nspec)
   integer :: idx_mmr_a_buf(gl_pm%nspec)
   integer :: idx_props_buf(gl_pm%nspec)

   character*1 :: gl_source_mmr_a(1)
   pointer(gl_source_mmr_a_ptr, gl_source_mmr_a)
   integer :: gl_idx_mmr_a(1), gl_idx_props_a(1)
   pointer(gl_idx_mmr_a_ptr, gl_idx_mmr_a)
   pointer(gl_idx_props_ptr, gl_idx_props_a)

   character*1 :: source_mmr_a(1)
   pointer(source_mmr_a_ptr, source_mmr_a)
   integer :: idx_mmr_a(1), idx_props_a(1)
   pointer(idx_mmr_a_ptr, idx_mmr_a)
   pointer(idx_props_ptr, idx_props_a)

   real(r8) :: density_aer_buf(numphysprops)
   real(r8) :: hygro_aer_buf(numphysprops)
   integer(8) :: refindex_aer_lw_buf(numphysprops)
   real(r8) :: density_aer_a(1), hygro_aer_a(1)
   pointer(density_aer_ptr, density_aer_a)
   pointer(hygro_aer_ptr, hygro_aer_a)
   integer(8) :: refindex_aer_lw_a(1)
   pointer(refindex_aer_lw_ptr, refindex_aer_lw_a)

! --- Slave ---
   integer :: myid, rid, cid
   type(modal_aero_lw_param) :: pm

   integer :: pver
   integer :: top_lev
   real(r8) :: rhoh2o
   integer :: nlwbands
   integer :: ncoef, prefr, prefi
   real(r8) :: xrmax, xrmin

   character*1 :: s_A

interface
   integer function trim_strcmp(a, b)
      character(len=*) :: a, b
   end function
end interface

   ! call set_func_ptr
   call athread_get_id(myid)
   rid = myid / 8
   cid = mod(myid, 8)

   call dmemcpy(pm, gl_pm, loc(gl_pm%end_type) - loc(gl_pm))

   s_A = 'A'

   pm_qaerwat_ptr = loc(pm%qaerwat)
   pm_dgnumwet_ptr = loc(pm%dgnumwet)

   pver = mod_pver
   top_lev = clim_modal_aero_top_lev
   rhoh2o = mod_rhoh2o
   xrmax = mod_xrmax
   xrmin = mod_xrmin
   nlwbands = mod_nlwbands

   ncoef = cnst_ncoef
   prefr = cnst_prefr
   prefi = cnst_prefi

   list_idx = pm%list_idx
   ncol = pm%ncol
   nmodes = pm%nmodes
   m = pm%imode
   nspec = pm%nspec

   sigma_logr_aer = pm%sigma_logr_aer

   gl_source_mmr_a_ptr = pm%source_mmr_a_ptr
   gl_idx_mmr_a_ptr = pm%idx_mmr_a_ptr
   gl_idx_props_ptr = pm%idx_props_ptr
   source_mmr_a_buf = gl_source_mmr_a(1:nspec)
   call dmemcpy(idx_mmr_a_buf, gl_idx_mmr_a, nspec * kind(idx_mmr_a_buf))
   call dmemcpy(idx_props_buf, gl_idx_props_a, nspec * kind(idx_props_buf))
   source_mmr_a_ptr = loc(source_mmr_a_buf)
   idx_mmr_a_ptr    = loc(idx_mmr_a_buf)
   idx_props_ptr    = loc(idx_props_buf)

   do i = 1, numphysprops
      density_aer_buf(i) = physprop(i)%density_aer
      refindex_aer_lw_buf(i) = loc(physprop(i)%refindex_aer_lw)
   enddo
   density_aer_ptr = loc(density_aer_buf)
   refindex_aer_lw_ptr = loc(refindex_aer_lw_buf)

   call dmemcpy(crefwlw, mod_crefwlw, nlwbands * 2*kind(crefwlw))

      do k = rid + 1, pver, 8
         if (k < top_lev) cycle
         call dmemcpy(dgnumwet, pm_dgnumwet(:,k), ncol * kind(dgnumwet))
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

      do ilw = cid + 1, nlwbands, 8
         call dmemcpy(refrtablw, pm%refrtablw(:,ilw), prefr * kind(refrtablw))
         call dmemcpy(refitablw, pm%refitablw(:,ilw), prefi * kind(refitablw))
         call dmemcpy(absplw, pm%absplw(:,:,:,ilw), ncoef*prefr*prefi * kind(absplw))

         do k = rid + 1, pver, 8
            if (k < top_lev) cycle
            call dmemcpy(mass, pm%mass(:,k), ncol * kind(mass))
            call dmemcpy(qaerwat, pm_qaerwat(:,k), ncol * kind(qaerwat))

            ! form bulk refractive index. Use volume mixing for infrared
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               source = source_mmr_a(l)
               idx    = idx_mmr_a(l)

               if (trim_strcmp(source, s_A) /= 0) call abort
               gl_specmmr_ptr = pm%state_q_ptr + ((idx-1) * pm%state_q_dim2 + (k-1)) * pm%state_q_dim1 * kind(specmmr)

               id = idx_props_a(l)
               specdens = density_aer_a(id)
               gl_specrefindex_ptr = refindex_aer_lw_a(id)
               hygro_aer = hygro_aer_a(id)

               call dmemcpy(specmmr, gl_specmmr, ncol * kind(specmmr))
               call dmemcpy(specrefindex, gl_specrefindex(:), nlwbands * 2*kind(specrefindex))

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

            call dmemcpy(tauxar, pm%tauxar(:,k,ilw), ncol * kind(tauxar))
            do i = 1, ncol
               tauxar(i) = tauxar(i) + dopaer(i)
            end do
            call dmemcpy(pm%tauxar(:,k,ilw), tauxar, ncol * kind(tauxar))

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

end module modal_aer_opt_cpe
