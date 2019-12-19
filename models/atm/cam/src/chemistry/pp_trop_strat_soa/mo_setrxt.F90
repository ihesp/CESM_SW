
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

      rate(:,:,92) = 1.20e-10_r8
      rate(:,:,93) = 2.02e-10_r8
      rate(:,:,94) = 1.204e-10_r8
      rate(:,:,95) = 1.50e-10_r8
      rate(:,:,96) = 9.75e-11_r8
      rate(:,:,97) = 1.50e-11_r8
      rate(:,:,98) = 7.20e-11_r8
      rate(:,:,99) = 1.794e-10_r8
      rate(:,:,100) = 1.628e-10_r8
      rate(:,:,101) = 2.84e-10_r8
      rate(:,:,102) = 1.674e-10_r8
      rate(:,:,103) = 9.60e-11_r8
      rate(:,:,104) = 4.10e-11_r8
      rate(:,:,105) = 1.012e-10_r8
      rate(:,:,106) = 1.20e-10_r8
      rate(:,:,107) = 4.49e-10_r8
      rate(:,:,108) = 2.57e-10_r8
      rate(:,:,109) = 1.31e-10_r8
      rate(:,:,110) = 3.50e-11_r8
      rate(:,:,111) = 9.00e-12_r8
      rate(:,:,112) = 1.20e-10_r8
      rate(:,:,113) = 1.50e-10_r8
      rate(:,:,114) = 1.20e-10_r8
      rate(:,:,118) = 7.20e-11_r8
      rate(:,:,119) = 6.90e-12_r8
      rate(:,:,120) = 1.60e-12_r8
      rate(:,:,124) = 1.80e-12_r8
      rate(:,:,127) = 1.80e-12_r8
      rate(:,:,151) = 1.00e-11_r8
      rate(:,:,152) = 2.20e-11_r8
      rate(:,:,153) = 3.50e-12_r8
      rate(:,:,178) = 1.70e-13_r8
      rate(:,:,225) = 4.50e-13_r8
      rate(:,:,235) = 1.00e-14_r8
      rate(:,:,238) = 7.00e-13_r8
      rate(:,:,241) = 2.00e-13_r8
      rate(:,:,242) = 6.80e-14_r8
      rate(:,:,251) = 1.00e-12_r8
      rate(:,:,252) = 1.00e-11_r8
      rate(:,:,253) = 1.15e-11_r8
      rate(:,:,256) = 4.00e-14_r8
      rate(:,:,273) = 3.00e-12_r8
      rate(:,:,276) = 6.80e-13_r8
      rate(:,:,277) = 5.40e-11_r8
      rate(:,:,289) = 2.40e-12_r8
      rate(:,:,292) = 1.40e-11_r8
      rate(:,:,295) = 5.00e-12_r8
      rate(:,:,307) = 2.40e-12_r8
      rate(:,:,311) = 1.40e-11_r8
      rate(:,:,313) = 2.40e-12_r8
      rate(:,:,315) = 3.50e-12_r8
      rate(:,:,316) = 4.50e-11_r8
      rate(:,:,323) = 2.40e-12_r8
      rate(:,:,333) = 3.00e-12_r8
      rate(:,:,334) = 1.00e-11_r8
      rate(:,:,338) = 2.3e-11_r8
      rate(:,:,347) = 2.10e-6_r8
      rate(:,:,351) = 7.10e-6_r8
      rate(:,:,357) = 7.10e-6_r8
      rate(:,:,359) = 6.34e-8_r8
      rate(:,:,360) = 6.34e-8_r8
      rate(:,:,361) = 6.34e-8_r8
      rate(:,:,362) = 6.34e-8_r8
      rate(:,:,363) = 6.34e-8_r8
      rate(:,:,364) = 6.34e-8_r8
      rate(:,:,365) = 6.34e-8_r8
      rate(:,:,366) = 6.34e-8_r8
      rate(:,:,367) = 6.34e-8_r8
      rate(:,:,368) = 6.34e-8_r8
      rate(:,:,369) = 6.34e-8_r8
      rate(:,:,370) = 6.34e-8_r8
      rate(:,:,371) = 6.34e-8_r8
      rate(:,:,372) = 6.34e-8_r8
      rate(:,:,373) = 6.34e-8_r8
      rate(:,:,374) = 6.34e-8_r8
      rate(:,:,375) = 6.34e-8_r8
      rate(:,:,376) = 6.34e-8_r8
      rate(:,:,377) = 6.34e-8_r8
      rate(:,:,378) = 6.34e-8_r8
      rate(:,:,379) = 6.34e-8_r8
      rate(:,:,397) = 2.31e-06_r8
      rate(:,:,398) = 2.31e-07_r8
      rate(:,:,399) = 2.31e-07_r8
      rate(:,:,400) = 4.63e-07_r8
      rate(:,:,401) = 4.63e-07_r8
      rate(:,:,402) = 2.31e-07_r8
      rate(:,:,403) = 1.29e-07_r8
      rate(:,:,404) = 1.29e-07_r8
      rate(:,:,405) = 1.29e-07_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,85) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,87) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      rate(:,:,88) = 3.30e-11_r8 * exp( 55._r8 * itemp(:,:) )
      rate(:,:,89) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,90) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,91) = 4.63e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 100._r8 * itemp(:,:) )
      rate(:,:,115) = 7.70e-11_r8 * exp_fac(:,:)
      rate(:,:,136) = 2.10e-11_r8 * exp_fac(:,:)
      rate(:,:,117) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,121) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,233) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,260) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,265) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,278) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,282) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,306) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,319) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,330) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,344) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,122) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,123) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,188) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,126) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,128) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,129) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,196) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,224) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,243) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,263) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,267) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,272) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,284) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,293) = 2.30e-11_r8 * exp_fac(:,:)
      rate(:,:,309) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,321) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,332) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,346) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,130) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,132) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,304) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,134) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      rate(:,:,135) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,137) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,138) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,139) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,141) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,160) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,165) = 7.40e-12_r8 * exp_fac(:,:)
      rate(:,:,246) = 8.10e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,142) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,197) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,143) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,145) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,171) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,150) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,155) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,157) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,158) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,159) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,161) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,162) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,163) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,164) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,166) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,187) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,195) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,167) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,169) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,194) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,168) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,172) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,173) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,176) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,177) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,179) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      rate(:,:,180) = 3.00e-12_r8 * exp( -500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,181) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,208) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,182) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,183) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,184) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,185) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,186) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,210) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,189) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,190) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,193) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,192) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,198) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,199) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,200) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,201) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,202) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,203) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,206) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,217) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,204) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,205) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,207) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,209) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,211) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,212) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,215) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,216) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,218) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,219) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,269) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,220) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,221) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,222) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,223) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,226) = 2.40e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,227) = 2.60e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,228) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,234) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,240) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,261) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,266) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,270) = 8.60e-13_r8 * exp_fac(:,:)
      rate(:,:,283) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,290) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,308) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,314) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,320) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,324) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,331) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,336) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,339) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,345) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,229) = 7.20e-11_r8 * exp( -70._r8 * itemp(:,:) )
      rate(:,:,231) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,236) = 1.60e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,237) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      rate(:,:,239) = 2.60e-12_r8 * exp( 365._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 350._r8 * itemp(:,:) )
      rate(:,:,244) = 4.63e-12_r8 * exp_fac(:,:)
      rate(:,:,337) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,340) = 2.6e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,245) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,258) = 6.50e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,248) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,296) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,249) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,250) = 2.50e-12_r8 * exp_fac(:,:)
      rate(:,:,271) = 7.10e-13_r8 * exp_fac(:,:)
      rate(:,:,297) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,254) = 6.90e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,259) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,262) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,264) = 8.70e-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,274) = 8.40e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,275) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,317) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,279) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,280) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      rate(:,:,281) = 2.30e-12_r8 * exp( -170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,285) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,318) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,286) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,287) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,288) = 1.30e-13_r8 * exp_fac(:,:)
      rate(:,:,294) = 5.30e-12_r8 * exp_fac(:,:)
      rate(:,:,312) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,322) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,291) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,310) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,325) = 5.e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,298) = 4.60e-12_r8 * exp_fac(:,:)
      rate(:,:,299) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,303) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      rate(:,:,305) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      rate(:,:,326) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,327) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,329) = 1.70e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,335) = 2.3e-12_r8 * exp( -193._r8 * itemp(:,:) )
      rate(:,:,341) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,342) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,343) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,353) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,355) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,356) = 1.70e-12_r8 * exp( -710._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,116), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,125), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,133), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,140), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,144), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,146), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,148), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,154), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,170), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,174), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,191), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,214), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.50e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2.0_r8)
      call jpl( rate(1,1,230), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.60e-29_r8 * itemp(:,:)**3.1_r8
      kinf(:,:) = 9.00e-12_r8 * itemp(:,:)**0.85_r8
      call jpl( rate(1,1,232), m, 0.48_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,247), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,257), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,302), m, 0.5_r8, ko, kinf, n )

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


      end subroutine setrxt_hrates

      end module mo_setrxt
