
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,92) = 8.00e-14_r8
      rate(:,:,93) = 3.90e-17_r8
      rate(:,:,96) = 4.20e-13_r8
      rate(:,:,97) = 8.50e-2_r8
      rate(:,:,98) = 1.30e-16_r8
      rate(:,:,100) = 1.00e-20_r8
      rate(:,:,101) = 2.58e-04_r8
      rate(:,:,108) = 1.20e-10_r8
      rate(:,:,109) = 2.02e-10_r8
      rate(:,:,110) = 1.204e-10_r8
      rate(:,:,111) = 1.50e-10_r8
      rate(:,:,112) = 9.75e-11_r8
      rate(:,:,113) = 1.50e-11_r8
      rate(:,:,114) = 7.20e-11_r8
      rate(:,:,115) = 1.794e-10_r8
      rate(:,:,116) = 1.628e-10_r8
      rate(:,:,117) = 2.84e-10_r8
      rate(:,:,118) = 1.674e-10_r8
      rate(:,:,119) = 9.60e-11_r8
      rate(:,:,120) = 4.10e-11_r8
      rate(:,:,121) = 1.012e-10_r8
      rate(:,:,122) = 1.20e-10_r8
      rate(:,:,123) = 4.49e-10_r8
      rate(:,:,124) = 2.57e-10_r8
      rate(:,:,125) = 2.14e-11_r8
      rate(:,:,126) = 1.90e-10_r8
      rate(:,:,127) = 1.31e-10_r8
      rate(:,:,128) = 3.50e-11_r8
      rate(:,:,129) = 9.00e-12_r8
      rate(:,:,130) = 1.20e-10_r8
      rate(:,:,131) = 1.50e-10_r8
      rate(:,:,132) = 1.20e-10_r8
      rate(:,:,135) = 7.20e-11_r8
      rate(:,:,136) = 6.90e-12_r8
      rate(:,:,137) = 1.60e-12_r8
      rate(:,:,141) = 1.80e-12_r8
      rate(:,:,144) = 1.80e-12_r8
      rate(:,:,150) = 5.00e-12_r8
      rate(:,:,151) = 7.00e-13_r8
      rate(:,:,152) = 5.00e-11_r8
      rate(:,:,169) = 1.00e-11_r8
      rate(:,:,170) = 2.20e-11_r8
      rate(:,:,171) = 3.50e-12_r8
      rate(:,:,196) = 1.70e-13_r8
      rate(:,:,268) = 9.0e-10_r8
      rate(:,:,269) = 1.0e-10_r8
      rate(:,:,270) = 4.4e-10_r8
      rate(:,:,271) = 4.0e-10_r8
      rate(:,:,272) = 2.0e-10_r8
      rate(:,:,273) = 1.0e-12_r8
      rate(:,:,274) = 6.0e-11_r8
      rate(:,:,275) = 5.0e-16_r8
      rate(:,:,279) = 2.31e-06_r8
      rate(:,:,280) = 2.31e-07_r8
      rate(:,:,281) = 2.31e-07_r8
      rate(:,:,282) = 4.63e-07_r8
      rate(:,:,283) = 4.63e-07_r8
      rate(:,:,284) = 2.31e-07_r8
      rate(:,:,285) = 1.29e-07_r8
      rate(:,:,286) = 1.29e-07_r8
      rate(:,:,287) = 1.29e-07_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,90) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,94) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,95) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,99) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,102) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,103) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,104) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,105) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,106) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,107) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,134) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,138) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:,139) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,140) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,206) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,143) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,145) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,146) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,214) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,242) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,147) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,149) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,153) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,154) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,155) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,156) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,157) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,159) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,178) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,183) = 7.40e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,160) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,215) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,161) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,163) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,189) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,168) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,173) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,175) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,176) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,177) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,179) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,180) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,181) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,182) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,184) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,205) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,213) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,185) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,187) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,212) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,186) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,190) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,191) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,194) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,195) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,197) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,198) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,219) = 1.40e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,199) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,230) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,200) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,201) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,202) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,203) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,204) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,232) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,207) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,208) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,211) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,210) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,216) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,217) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,218) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,268) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,269) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,270) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,271) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,272) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,273) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,274) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,275) = 5.0e-16_r8 * exp_fac(:,:)
      rate(:,:,279) = 2.31e-06_r8 * exp_fac(:,:)
      rate(:,:,280) = 2.31e-07_r8 * exp_fac(:,:)
      rate(:,:,281) = 2.31e-07_r8 * exp_fac(:,:)
      rate(:,:,282) = 4.63e-07_r8 * exp_fac(:,:)
      rate(:,:,283) = 4.63e-07_r8 * exp_fac(:,:)
      rate(:,:,284) = 2.31e-07_r8 * exp_fac(:,:)
      rate(:,:,285) = 1.29e-07_r8 * exp_fac(:,:)
      rate(:,:,286) = 1.29e-07_r8 * exp_fac(:,:)
      rate(:,:,287) = 1.29e-07_r8 * exp_fac(:,:)
      rate(:,:,220) = 1.60e-10_r8 * exp( -260._r8 * itemp(:,:) )
      rate(:,:,221) = 6.00e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,222) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,223) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,224) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,225) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,228) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,239) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,226) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,227) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,229) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,231) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,233) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,234) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,237) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,238) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,240) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,241) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,133), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,142), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,158), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,162), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,164), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,166), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,172), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,188), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,192), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,209), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,236), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)

      rate(:,:kbot,92) = 8.00e-14_r8
      rate(:,:kbot,93) = 3.90e-17_r8
      rate(:,:kbot,98) = 1.30e-16_r8
      rate(:,:kbot,100) = 1.00e-20_r8
      rate(:,:kbot,136) = 6.90e-12_r8
      rate(:,:kbot,150) = 5.00e-12_r8
      rate(:,:kbot,151) = 7.00e-13_r8
      rate(:,:kbot,269) = 1.0e-10_r8
      rate(:,:kbot,270) = 4.4e-10_r8
      rate(:,:kbot,271) = 4.0e-10_r8
      rate(:,:kbot,272) = 2.0e-10_r8
      rate(:,:kbot,273) = 1.0e-12_r8
      rate(:,:kbot,274) = 6.0e-11_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,90) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,94) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,95) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,99) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,102) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,103) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,104) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,134) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,138) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,139) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,140) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,146) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,147) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,153) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,154) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,159) = 3.30e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,160) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,161) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,133) = wrk(:,:)











      end subroutine setrxt_hrates

      end module mo_setrxt
