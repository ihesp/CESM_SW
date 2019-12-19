




      module mo_lin_matrix

      private
      public :: linmat

      contains

      subroutine linmat01( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)

         mat(518) = -( rxt(3) + rxt(4) + het_rates(1) )

         mat(395) = -( rxt(54) + rxt(55) + rxt(56) + rxt(67) + rxt(68) + rxt(69) &
                 + het_rates(2) )
         mat(272) = rxt(1) + 2.000_r8*rxt(2) + rxt(60) + rxt(61) + rxt(62) &
                      + 2.000_r8*rxt(65) + rxt(72) + rxt(73) + rxt(74) + 2.000_r8*rxt(77)
         mat(513) = rxt(4)
         mat(418) = rxt(6)
         mat(442) = rxt(8)
         mat(46) = rxt(10)
         mat(492) = rxt(12)
         mat(251) = rxt(21)
         mat(313) = rxt(24)
         mat(60) = rxt(25)
         mat(290) = rxt(32)
         mat(217) = rxt(50)
         mat(37) = rxt(51)
         mat(596) = rxt(53)
         mat(332) = rxt(93)

         mat(330) = -( rxt(93) + rxt(97)*y(7) + rxt(98)*y(7) + rxt(100)*y(43) &
                      + rxt(101)*y(44) + rxt(102)*y(45) + rxt(103)*y(46) + rxt(104)*y(47) &
                      + rxt(105)*y(42) + rxt(106)*y(50) + rxt(107)*y(49) + rxt(108)*y(15) &
                      + rxt(109)*y(15) + rxt(110)*y(15) + rxt(111)*y(20) + het_rates(3) )
         mat(270) = rxt(1)
         mat(511) = rxt(3)
         mat(249) = rxt(20)

         mat(269) = -( rxt(1) + rxt(2) + rxt(58) + rxt(60) + rxt(61) + rxt(62) + rxt(65) &
                      + rxt(70) + rxt(72) + rxt(73) + rxt(74) + rxt(77) + het_rates(4) )
         mat(508) = rxt(4)
         mat(489) = rxt(13)
         mat(27) = rxt(88)
         mat(24) = rxt(91) + rxt(92)
         mat(329) = rxt(98)*y(7)

         mat(26) = -( rxt(85) + rxt(88) + rxt(87)*y(51) + het_rates(5) )

         mat(23) = -( rxt(91) + rxt(92) + het_rates(6) )
         mat(502) = rxt(3)
         mat(25) = rxt(85) + rxt(87)*y(51)

         mat(177) = -( rxt(57) + het_rates(8) )
         mat(408) = rxt(6)
         mat(85) = rxt(247)

         mat(419) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(443) = rxt(8) + .500_r8*rxt(205)
         mat(47) = rxt(10)
         mat(493) = rxt(13)
         mat(153) = rxt(256)
         mat(333) = 2.000_r8*rxt(97)*y(7)

         mat(444) = -( rxt(8) + rxt(205) + het_rates(10) )
         mat(48) = rxt(9) + rxt(126)
         mat(134) = rxt(11)
         mat(494) = rxt(12)
         mat(79) = rxt(15) + rxt(135)
         mat(208) = rxt(30)
         mat(98) = rxt(35)

         mat(476) = -( rxt(136)*y(15) + rxt(143)*y(19) + rxt(144)*y(19) + rxt(155)*y(20) &
                      + rxt(208)*y(41) + rxt(209)*y(48) + rxt(210)*y(46) + rxt(211)*y(42) &
                 + het_rates(22) )
         mat(135) = rxt(11)
         mat(80) = rxt(14)
         mat(67) = rxt(16)
         mat(252) = rxt(19)
         mat(116) = 2.000_r8*rxt(22)
         mat(198) = rxt(27)
         mat(142) = rxt(33)
         mat(445) = .500_r8*rxt(205)
         mat(335) = rxt(108)*y(15) + rxt(111)*y(20)

         mat(496) = -( rxt(12) + rxt(13) + rxt(204) + het_rates(11) )
         mat(49) = rxt(9) + rxt(10) + rxt(126)
         mat(81) = rxt(14)
         mat(210) = rxt(29)
         mat(99) = rxt(34)

         mat(132) = -( rxt(11) + het_rates(12) )
         mat(45) = 2.000_r8*rxt(203) + 2.000_r8*rxt(229) + 2.000_r8*rxt(235) &
                      + 2.000_r8*rxt(240)
         mat(485) = rxt(204)
         mat(431) = .500_r8*rxt(205)
         mat(202) = rxt(230) + rxt(236) + rxt(241)
         mat(94) = rxt(231) + rxt(239) + rxt(242)

         mat(75) = -( rxt(14) + rxt(15) + rxt(135) + het_rates(13) )

         mat(44) = -( rxt(9) + rxt(10) + rxt(126) + rxt(203) + rxt(229) + rxt(235) &
                      + rxt(240) + het_rates(14) )

         mat(184) = -( het_rates(16) )
         mat(325) = rxt(108)*y(15)
         mat(461) = rxt(136)*y(15)
         mat(525) = rxt(167)*y(15)

         mat(62) = -( rxt(16) + het_rates(17) )

         mat(230) = -( rxt(17) + rxt(18) + het_rates(18) )
         mat(65) = rxt(16)
         mat(327) = rxt(109)*y(15) + rxt(110)*y(15)

         mat(220) = -( het_rates(21) )
         mat(64) = rxt(16)
         mat(229) = 2.000_r8*rxt(17)
         mat(246) = rxt(19) + 2.000_r8*rxt(21)
         mat(551) = rxt(28)
         mat(326) = rxt(109)*y(15) + rxt(111)*y(20)
         mat(465) = rxt(144)*y(19) + rxt(155)*y(20)
         mat(528) = rxt(162)*y(20)

         mat(356) = -( rxt(206) + het_rates(23) )
         mat(78) = rxt(15) + rxt(135)
         mat(331) = rxt(109)*y(15)
         mat(472) = rxt(143)*y(19) + rxt(208)*y(41) + rxt(211)*y(42)
         mat(534) = rxt(207)*y(41)

         mat(112) = -( rxt(22) + het_rates(24) )
         mat(344) = .500_r8*rxt(206)

         mat(247) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(82) )
         mat(22) = rxt(49)
         mat(467) = rxt(136)*y(15) + rxt(155)*y(20) + rxt(208)*y(41) + rxt(209)*y(48) &
                      + rxt(210)*y(46) + rxt(211)*y(42)

         mat(541) = -( rxt(162)*y(20) + rxt(167)*y(15) + rxt(207)*y(41) + het_rates(27) )
         mat(29) = 2.000_r8*rxt(23)
         mat(319) = rxt(24)
         mat(19) = 2.000_r8*rxt(26)
         mat(199) = rxt(27)
         mat(564) = rxt(28)
         mat(211) = rxt(29)
         mat(31) = rxt(31)
         mat(338) = 3.000_r8*rxt(100)*y(43) + 2.000_r8*rxt(101)*y(44) &
                      + 3.000_r8*rxt(102)*y(45) + rxt(103)*y(46) + 4.000_r8*rxt(104)*y(47)
         mat(479) = rxt(208)*y(41) + 3.000_r8*rxt(209)*y(48) + rxt(210)*y(46)

         mat(28) = -( rxt(23) + het_rates(28) )

         mat(310) = -( rxt(24) + het_rates(29) )
         mat(59) = rxt(25)
         mat(206) = rxt(30)
         mat(18) = 2.000_r8*rxt(178)

         mat(57) = -( rxt(25) + het_rates(30) )

         mat(17) = -( rxt(26) + rxt(178) + het_rates(31) )

         mat(565) = -( rxt(28) + het_rates(32) )
         mat(542) = rxt(162)*y(20) + rxt(167)*y(15) + 2.000_r8*rxt(207)*y(41)

         mat(194) = -( rxt(27) + het_rates(33) )
         mat(203) = rxt(230) + rxt(236) + rxt(241)

         mat(204) = -( rxt(29) + rxt(30) + rxt(230) + rxt(236) + rxt(241) + het_rates(34) &
       )

         mat(30) = -( rxt(31) + het_rates(35) )

         mat(584) = -( het_rates(36) )
         mat(32) = rxt(31)
         mat(298) = rxt(32)
         mat(145) = rxt(33)
         mat(100) = rxt(34)
         mat(340) = rxt(105)*y(42) + rxt(106)*y(50) + rxt(107)*y(49)
         mat(481) = rxt(211)*y(42)

         mat(286) = -( rxt(32) + het_rates(37) )
         mat(96) = rxt(35)

         mat(126) = -( het_rates(38) )

         mat(138) = -( rxt(33) + het_rates(39) )
         mat(95) = rxt(231) + rxt(239) + rxt(242)

         mat(93) = -( rxt(34) + rxt(35) + rxt(231) + rxt(239) + rxt(242) + het_rates(40) &
       )

         mat(103) = -( het_rates(52) )

         mat(146) = -( rxt(256) + het_rates(53) )
         mat(262) = rxt(58) + rxt(70)
         mat(83) = rxt(249)*y(51)

         mat(68) = -( het_rates(54) )
         mat(172) = rxt(57)

         mat(82) = -( rxt(247) + rxt(249)*y(51) + het_rates(55) )
         mat(371) = rxt(54) + rxt(55) + rxt(56) + rxt(67) + rxt(68) + rxt(69)
         mat(258) = rxt(60) + rxt(61) + rxt(62) + rxt(72) + rxt(73) + rxt(74)

         mat(166) = -( het_rates(56) )
         mat(407) = rxt(7)
         mat(84) = rxt(247)
         mat(148) = rxt(256)

         mat(88) = -( het_rates(58) )

         mat(157) = -( het_rates(57) )
         mat(406) = rxt(7)
         mat(381) = rxt(54) + rxt(55) + rxt(56) + rxt(67) + rxt(68) + rxt(69)
         mat(175) = rxt(57)
         mat(263) = rxt(58) + rxt(60) + rxt(61) + rxt(62) + rxt(70) + rxt(72) + rxt(73) &
                      + rxt(74)

         mat(50) = -( rxt(52) + het_rates(59) )

         mat(119) = -( het_rates(60) )
         mat(51) = rxt(52)
         mat(587) = rxt(53)

         mat(605) = -( rxt(53) + het_rates(61) )
         mat(219) = rxt(50)

         mat(214) = -( rxt(50) + het_rates(62) )
         mat(35) = rxt(51)

         mat(34) = -( rxt(51) + het_rates(63) )
         mat(21) = rxt(49)

         mat(20) = -( rxt(49) + het_rates(64) )

         mat(38) = -( het_rates(65) )

         mat(1) = -( het_rates(66) )

         mat(2) = -( het_rates(67) )

         mat(3) = -( het_rates(68) )

         mat(4) = -( het_rates(69) )

         mat(5) = -( het_rates(70) )

         mat(6) = -( het_rates(71) )

         mat(7) = -( het_rates(72) )

         mat(8) = -( het_rates(73) )

         mat(9) = -( het_rates(74) )

         mat(10) = -( het_rates(75) )

         mat(11) = -( het_rates(76) )

         mat(12) = -( het_rates(77) )

         mat(13) = -( het_rates(78) )

         mat(14) = -( het_rates(79) )

         mat(15) = -( het_rates(80) )

         mat(16) = -( het_rates(81) )


      end subroutine linmat01

      subroutine linmat( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)

      call linmat01( mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix
