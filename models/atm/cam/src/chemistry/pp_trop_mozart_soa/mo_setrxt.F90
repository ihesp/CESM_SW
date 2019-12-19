
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

      rate(:,:,49) = 1.20e-10_r8
      rate(:,:,51) = 1.20e-10_r8
      rate(:,:,54) = 1.31e-10_r8
      rate(:,:,55) = 3.50e-11_r8
      rate(:,:,56) = 9.00e-12_r8
      rate(:,:,60) = 7.20e-11_r8
      rate(:,:,61) = 6.90e-12_r8
      rate(:,:,62) = 1.60e-12_r8
      rate(:,:,66) = 1.80e-12_r8
      rate(:,:,69) = 1.80e-12_r8
      rate(:,:,88) = 1.00e-11_r8
      rate(:,:,89) = 2.20e-11_r8
      rate(:,:,90) = 3.50e-12_r8
      rate(:,:,107) = 4.50e-13_r8
      rate(:,:,116) = 1.00e-14_r8
      rate(:,:,119) = 7.00e-13_r8
      rate(:,:,122) = 2.00e-13_r8
      rate(:,:,123) = 6.80e-14_r8
      rate(:,:,132) = 1.00e-12_r8
      rate(:,:,133) = 1.00e-11_r8
      rate(:,:,134) = 1.15e-11_r8
      rate(:,:,137) = 4.00e-14_r8
      rate(:,:,154) = 3.00e-12_r8
      rate(:,:,157) = 6.80e-13_r8
      rate(:,:,158) = 5.40e-11_r8
      rate(:,:,170) = 2.40e-12_r8
      rate(:,:,173) = 1.40e-11_r8
      rate(:,:,176) = 5.00e-12_r8
      rate(:,:,188) = 2.40e-12_r8
      rate(:,:,192) = 1.40e-11_r8
      rate(:,:,194) = 2.40e-12_r8
      rate(:,:,196) = 3.50e-12_r8
      rate(:,:,197) = 4.50e-11_r8
      rate(:,:,204) = 2.40e-12_r8
      rate(:,:,214) = 3.00e-12_r8
      rate(:,:,215) = 1.00e-11_r8
      rate(:,:,219) = 2.3e-11_r8
      rate(:,:,228) = 2.10e-6_r8
      rate(:,:,232) = 7.10e-6_r8
      rate(:,:,238) = 7.10e-6_r8
      rate(:,:,240) = 2.31e-06_r8
      rate(:,:,241) = 2.31e-07_r8
      rate(:,:,242) = 2.31e-07_r8
      rate(:,:,243) = 4.63e-07_r8
      rate(:,:,244) = 4.63e-07_r8
      rate(:,:,245) = 2.31e-07_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,46) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,47) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      rate(:,:,48) = 3.30e-11_r8 * exp( 55._r8 * itemp(:,:) )
      rate(:,:,50) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,52) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,53) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,57) = 7.70e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:,59) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,63) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,114) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,141) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,146) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,159) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,163) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,187) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,200) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,211) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,225) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,64) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,65) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:,68) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,70) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,71) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,106) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,124) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,144) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,148) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,153) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,165) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,174) = 2.30e-11_r8 * exp_fac(:,:)
      rate(:,:,190) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,202) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,213) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,227) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,72) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,74) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,185) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,76) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,78) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,127) = 8.10e-12_r8 * exp_fac(:,:)
      rate(:,:,79) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:,80) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:,82) = 1.20e-13_r8 * exp( -2450._r8 * itemp(:,:) )
      rate(:,:,87) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,92) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,94) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,97) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,98) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,99) = 3.40e-11_r8 * exp( -1600._r8 * itemp(:,:) )
      rate(:,:,100) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,101) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,150) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,102) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,103) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,104) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,105) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,108) = 2.40e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,109) = 2.60e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,110) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,115) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,121) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,142) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,147) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,151) = 8.60e-13_r8 * exp_fac(:,:)
      rate(:,:,164) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,171) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,189) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,195) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,201) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,205) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,212) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,217) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,220) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,226) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,112) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,117) = 1.60e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,118) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      rate(:,:,120) = 2.60e-12_r8 * exp( 365._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 350._r8 * itemp(:,:) )
      rate(:,:,125) = 4.63e-12_r8 * exp_fac(:,:)
      rate(:,:,218) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,221) = 2.6e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,126) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,139) = 6.50e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,129) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,177) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,130) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,131) = 2.50e-12_r8 * exp_fac(:,:)
      rate(:,:,152) = 7.10e-13_r8 * exp_fac(:,:)
      rate(:,:,178) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,135) = 6.90e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,140) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,143) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,145) = 8.70e-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,155) = 8.40e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,156) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,198) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,160) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,161) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      rate(:,:,162) = 2.30e-12_r8 * exp( -170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,166) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,199) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,167) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,168) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,169) = 1.30e-13_r8 * exp_fac(:,:)
      rate(:,:,175) = 5.30e-12_r8 * exp_fac(:,:)
      rate(:,:,193) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,203) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,172) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,191) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,206) = 5.e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,179) = 4.60e-12_r8 * exp_fac(:,:)
      rate(:,:,180) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,184) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      rate(:,:,186) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      rate(:,:,207) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,208) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,210) = 1.70e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,216) = 2.3e-12_r8 * exp( -193._r8 * itemp(:,:) )
      rate(:,:,222) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,223) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,224) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,234) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,236) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,237) = 1.70e-12_r8 * exp( -710._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,58), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,67), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,75), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,77), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,81), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,83), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,85), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,91), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,96), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.50e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2.0_r8)
      call jpl( rate(1,1,111), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.60e-29_r8 * itemp(:,:)**3.1_r8
      kinf(:,:) = 9.00e-12_r8 * itemp(:,:)**0.85_r8
      call jpl( rate(1,1,113), m, 0.48_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,128), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,138), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,183), m, 0.5_r8, ko, kinf, n )

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
