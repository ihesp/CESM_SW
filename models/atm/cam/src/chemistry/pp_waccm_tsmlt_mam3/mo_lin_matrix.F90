




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

         mat(1016) = -( rxt(3) + rxt(4) + het_rates(1) )

         mat(943) = -( rxt(92) + rxt(93) + rxt(94) + rxt(105) + rxt(106) + rxt(107) &
                 + het_rates(2) )
         mat(904) = rxt(1) + 2.000_r8*rxt(2) + rxt(98) + rxt(99) + rxt(100) &
                      + 2.000_r8*rxt(103) + rxt(110) + rxt(111) + rxt(112) &
                      + 2.000_r8*rxt(115)
         mat(1014) = rxt(4)
         mat(1244) = rxt(6)
         mat(1281) = rxt(8)
         mat(103) = rxt(10)
         mat(1423) = rxt(12)
         mat(1471) = rxt(21)
         mat(1041) = rxt(24)
         mat(137) = rxt(25)
         mat(1189) = rxt(32)
         mat(554) = rxt(88)
         mat(82) = rxt(89)
         mat(808) = rxt(91)
         mat(969) = rxt(131)

         mat(970) = -( rxt(131) + rxt(135)*y(7) + rxt(136)*y(7) + rxt(138)*y(104) &
                      + rxt(139)*y(105) + rxt(140)*y(106) + rxt(141)*y(114) &
                      + rxt(142)*y(115) + rxt(143)*y(107) + rxt(144)*y(112) &
                      + rxt(145)*y(113) + rxt(146)*y(108) + rxt(147)*y(103) &
                      + rxt(148)*y(111) + rxt(149)*y(110) + rxt(150)*y(116) &
                      + rxt(151)*y(117) + rxt(152)*y(118) + rxt(153)*y(119) &
                      + rxt(156)*y(15) + rxt(157)*y(15) + rxt(158)*y(15) + het_rates(3) )
         mat(905) = rxt(1)
         mat(1015) = rxt(3)
         mat(1472) = rxt(20)

         mat(903) = -( rxt(1) + rxt(2) + rxt(96) + rxt(98) + rxt(99) + rxt(100) + rxt(103) &
                      + rxt(108) + rxt(110) + rxt(111) + rxt(112) + rxt(115) &
                 + het_rates(4) )
         mat(1013) = rxt(4)
         mat(1422) = rxt(13)
         mat(54) = rxt(126)
         mat(51) = rxt(129) + rxt(130)
         mat(968) = rxt(136)*y(7)

         mat(53) = -( rxt(123) + rxt(126) + rxt(125)*y(120) + het_rates(5) )

         mat(50) = -( rxt(129) + rxt(130) + het_rates(6) )
         mat(984) = rxt(3)
         mat(52) = rxt(123) + rxt(125)*y(120)

         mat(650) = -( het_rates(21) )
         mat(1490) = rxt(18)
         mat(1465) = rxt(20)
         mat(964) = rxt(158)*y(15)

         mat(602) = -( het_rates(20) )
         mat(1489) = rxt(17) + rxt(18)
         mat(606) = rxt(61)
         mat(636) = 1.340_r8*rxt(67)
         mat(735) = .700_r8*rxt(68)
         mat(661) = rxt(74)
         mat(531) = rxt(76)
         mat(511) = rxt(79)
         mat(256) = .450_r8*rxt(81)
         mat(376) = 2.000_r8*rxt(82)
         mat(145) = rxt(90)
         mat(1137) = rxt(254)*y(102)
         mat(300) = rxt(439)*y(120)

         mat(476) = -( rxt(95) + het_rates(8) )
         mat(1223) = rxt(6)
         mat(299) = rxt(436)

         mat(1252) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(1289) = rxt(8) + .500_r8*rxt(399)
         mat(104) = rxt(10)
         mat(1431) = rxt(13)
         mat(412) = rxt(446)
         mat(977) = 2.000_r8*rxt(135)*y(7)

         mat(1290) = -( rxt(8) + rxt(399) + het_rates(10) )
         mat(105) = rxt(9) + rxt(197)
         mat(1454) = rxt(11)
         mat(1432) = rxt(12)
         mat(218) = rxt(15) + rxt(206)
         mat(565) = rxt(30)
         mat(285) = rxt(36)
         mat(197) = .600_r8*rxt(64) + rxt(311)
         mat(292) = rxt(65) + rxt(357)
         mat(534) = rxt(76)

         mat(1389) = -( rxt(255)*y(102) + rxt(256)*y(109) + rxt(257)*y(107) &
                      + rxt(258)*y(103) + rxt(260)*y(112) + rxt(261)*y(113) &
                      + rxt(262)*y(119) + rxt(263)*y(118) + rxt(266)*y(15) + het_rates(23) &
       )
         mat(1455) = rxt(11)
         mat(219) = rxt(14)
         mat(157) = rxt(16)
         mat(1481) = rxt(19)
         mat(317) = 2.000_r8*rxt(22)
         mat(491) = rxt(27)
         mat(403) = rxt(33)
         mat(265) = rxt(62)
         mat(230) = rxt(63)
         mat(129) = rxt(69)
         mat(43) = rxt(70)
         mat(170) = rxt(71)
         mat(175) = rxt(72)
         mat(132) = rxt(75)
         mat(332) = rxt(83)
         mat(119) = rxt(84)
         mat(165) = rxt(85)
         mat(214) = rxt(86)
         mat(1291) = .500_r8*rxt(399)
         mat(979) = rxt(156)*y(15)

         mat(1434) = -( rxt(12) + rxt(13) + rxt(398) + het_rates(11) )
         mat(106) = rxt(9) + rxt(10) + rxt(197)
         mat(220) = rxt(14)
         mat(567) = rxt(29)
         mat(286) = rxt(35)
         mat(199) = .400_r8*rxt(64)

         mat(1457) = -( rxt(11) + het_rates(12) )
         mat(107) = 2.000_r8*rxt(397) + 2.000_r8*rxt(418) + 2.000_r8*rxt(424) &
                      + 2.000_r8*rxt(429)
         mat(1435) = rxt(398)
         mat(1293) = .500_r8*rxt(399)
         mat(568) = rxt(419) + rxt(425) + rxt(430)
         mat(287) = rxt(420) + rxt(428) + rxt(431)

         mat(215) = -( rxt(14) + rxt(15) + rxt(206) + het_rates(13) )

         mat(102) = -( rxt(9) + rxt(10) + rxt(197) + rxt(397) + rxt(418) + rxt(424) &
                      + rxt(429) + het_rates(14) )

         mat(872) = -( het_rates(16) )
         mat(609) = rxt(61)
         mat(229) = rxt(63)
         mat(196) = .400_r8*rxt(64)
         mat(743) = .300_r8*rxt(68)
         mat(372) = rxt(73)
         mat(967) = rxt(156)*y(15)
         mat(1143) = rxt(213)*y(15)
         mat(435) = rxt(252)*y(15)
         mat(1377) = rxt(266)*y(15)

         mat(154) = -( rxt(16) + het_rates(17) )

         mat(57) = -( het_rates(42) )

         mat(17) = -( het_rates(43) )

         mat(1509) = -( rxt(17) + rxt(18) + het_rates(19) )
         mat(159) = rxt(16)
         mat(267) = rxt(62)
         mat(647) = 1.340_r8*rxt(66)
         mat(177) = rxt(72)
         mat(537) = rxt(76)
         mat(279) = .690_r8*rxt(77)
         mat(621) = rxt(78)
         mat(514) = rxt(79)
         mat(333) = .100_r8*rxt(83)
         mat(183) = rxt(280)
         mat(193) = 2.000_r8*rxt(292)
         mat(983) = rxt(157)*y(15) + rxt(158)*y(15)

         mat(1171) = -( het_rates(22) )
         mat(156) = rxt(16)
         mat(1501) = 2.000_r8*rxt(17)
         mat(1477) = rxt(19) + 2.000_r8*rxt(21)
         mat(830) = rxt(28)
         mat(456) = rxt(34)
         mat(74) = rxt(57)
         mat(975) = rxt(157)*y(15)

         mat(1114) = -( rxt(400) + het_rates(24) )
         mat(217) = rxt(15) + rxt(206)
         mat(610) = rxt(61)
         mat(264) = rxt(62)
         mat(643) = 1.340_r8*rxt(66) + .660_r8*rxt(67)
         mat(128) = rxt(69)
         mat(169) = rxt(71)
         mat(664) = rxt(74)
         mat(533) = rxt(76)
         mat(277) = rxt(77)
         mat(619) = rxt(78)
         mat(512) = 2.000_r8*rxt(79)
         mat(259) = .560_r8*rxt(81)
         mat(377) = 2.000_r8*rxt(82)
         mat(331) = .900_r8*rxt(83)
         mat(213) = rxt(86)
         mat(180) = rxt(280)
         mat(192) = rxt(292)
         mat(973) = rxt(157)*y(15)
         mat(1149) = rxt(254)*y(102) + rxt(259)*y(103)
         mat(1383) = rxt(255)*y(102) + rxt(258)*y(103)

         mat(312) = -( rxt(22) + het_rates(25) )
         mat(1075) = .500_r8*rxt(400)

         mat(1484) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(158) )
         mat(49) = rxt(87)
         mat(1392) = rxt(255)*y(102) + rxt(256)*y(109) + rxt(257)*y(107) + rxt(258)*y(103) &
                      + rxt(262)*y(119) + rxt(266)*y(15)

         mat(1150) = -( rxt(213)*y(15) + rxt(254)*y(102) + rxt(259)*y(103) &
                      + rxt(264)*y(119) + rxt(265)*y(118) + het_rates(28) )
         mat(56) = 2.000_r8*rxt(23)
         mat(1046) = rxt(24)
         mat(22) = 2.000_r8*rxt(26)
         mat(490) = rxt(27)
         mat(829) = rxt(28)
         mat(564) = rxt(29)
         mat(71) = rxt(31)
         mat(68) = rxt(56)
         mat(974) = 2.000_r8*rxt(138)*y(104) + 2.000_r8*rxt(139)*y(105) &
                      + 2.000_r8*rxt(140)*y(106) + 2.000_r8*rxt(141)*y(114) &
                      + rxt(142)*y(115) + rxt(143)*y(107) + rxt(144)*y(112) &
                      + rxt(145)*y(113) + 4.000_r8*rxt(146)*y(108) + rxt(148)*y(111)
         mat(1384) = rxt(255)*y(102) + 3.000_r8*rxt(256)*y(109) + rxt(257)*y(107) &
                      + rxt(260)*y(112) + rxt(261)*y(113)

         mat(55) = -( rxt(23) + het_rates(29) )

         mat(1044) = -( rxt(24) + het_rates(30) )
         mat(138) = rxt(25)
         mat(563) = rxt(30)
         mat(21) = 2.000_r8*rxt(225)

         mat(134) = -( rxt(25) + het_rates(31) )

         mat(20) = -( rxt(26) + rxt(225) + het_rates(32) )


      end subroutine linmat01

      subroutine linmat02( mat, y, rxt, het_rates )
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

         mat(824) = -( rxt(28) + het_rates(33) )
         mat(1141) = rxt(213)*y(15) + 2.000_r8*rxt(254)*y(102) + rxt(259)*y(103) &
                      + rxt(264)*y(119) + rxt(265)*y(118)

         mat(486) = -( rxt(27) + het_rates(34) )
         mat(559) = rxt(419) + rxt(425) + rxt(430)

         mat(560) = -( rxt(29) + rxt(30) + rxt(419) + rxt(425) + rxt(430) + het_rates(35) &
       )

         mat(69) = -( rxt(31) + het_rates(36) )

         mat(839) = -( het_rates(37) )
         mat(70) = rxt(31)
         mat(1187) = rxt(32)
         mat(399) = rxt(33)
         mat(453) = rxt(34)
         mat(282) = rxt(35)
         mat(966) = rxt(147)*y(103) + rxt(148)*y(111) + rxt(149)*y(110) &
                      + 2.000_r8*rxt(150)*y(116) + 2.000_r8*rxt(151)*y(117) &
                      + 3.000_r8*rxt(152)*y(118) + 2.000_r8*rxt(153)*y(119)
         mat(1376) = rxt(258)*y(103) + 2.000_r8*rxt(262)*y(119) + 3.000_r8*rxt(263)*y(118)
         mat(1142) = rxt(259)*y(103) + 2.000_r8*rxt(264)*y(119) + 3.000_r8*rxt(265)*y(118)

         mat(1196) = -( rxt(32) + het_rates(38) )
         mat(284) = rxt(36)

         mat(452) = -( rxt(34) + het_rates(39) )

         mat(397) = -( rxt(33) + het_rates(40) )
         mat(281) = rxt(420) + rxt(428) + rxt(431)

         mat(280) = -( rxt(35) + rxt(36) + rxt(420) + rxt(428) + rxt(431) + het_rates(41) &
       )

         mat(344) = -( het_rates(121) )

         mat(405) = -( rxt(446) + het_rates(122) )
         mat(894) = rxt(96) + rxt(108)
         mat(297) = rxt(439)*y(120)

         mat(201) = -( het_rates(123) )
         mat(471) = rxt(95)

         mat(296) = -( rxt(436) + rxt(439)*y(120) + het_rates(124) )
         mat(923) = rxt(92) + rxt(93) + rxt(94) + rxt(105) + rxt(106) + rxt(107)
         mat(891) = rxt(98) + rxt(99) + rxt(100) + rxt(110) + rxt(111) + rxt(112)

         mat(414) = -( het_rates(125) )
         mat(1219) = rxt(7)
         mat(298) = rxt(436)
         mat(406) = rxt(446)

         mat(222) = -( het_rates(127) )

         mat(425) = -( het_rates(126) )
         mat(1220) = rxt(7)
         mat(930) = rxt(92) + rxt(93) + rxt(94) + rxt(105) + rxt(106) + rxt(107)
         mat(475) = rxt(95)
         mat(896) = rxt(96) + rxt(98) + rxt(99) + rxt(100) + rxt(108) + rxt(110) &
                      + rxt(111) + rxt(112)

         mat(587) = -( het_rates(59) )
         mat(734) = .700_r8*rxt(68)

         mat(494) = -( het_rates(83) )

         mat(442) = -( het_rates(64) )

         mat(607) = -( rxt(61) + het_rates(50) )
         mat(262) = rxt(62)
         mat(127) = rxt(69)
         mat(329) = .400_r8*rxt(83)
         mat(117) = rxt(84)

         mat(319) = -( het_rates(49) )

         mat(260) = -( rxt(62) + het_rates(65) )

         mat(789) = -( het_rates(48) )
         mat(195) = .600_r8*rxt(64) + rxt(311)
         mat(641) = 1.340_r8*rxt(66)
         mat(742) = .300_r8*rxt(68)
         mat(174) = rxt(72)
         mat(371) = rxt(73)
         mat(663) = rxt(74)
         mat(618) = rxt(78)
         mat(187) = rxt(80)
         mat(258) = .130_r8*rxt(81)
         mat(118) = rxt(84)

         mat(227) = -( rxt(63) + het_rates(54) )

         mat(194) = -( rxt(64) + rxt(311) + het_rates(58) )

         mat(150) = -( het_rates(82) )

         mat(84) = -( het_rates(45) )

         mat(233) = -( het_rates(44) )

         mat(23) = -( het_rates(71) )

         mat(288) = -( rxt(65) + rxt(357) + het_rates(81) )

         mat(26) = -( het_rates(70) )

         mat(108) = -( het_rates(73) )

         mat(358) = -( het_rates(84) )

         mat(324) = -( rxt(83) + het_rates(85) )

         mat(184) = -( rxt(80) + het_rates(72) )
         mat(323) = .800_r8*rxt(83)

         mat(335) = -( het_rates(74) )

         mat(115) = -( rxt(84) + het_rates(75) )

         mat(33) = -( het_rates(94) )

         mat(38) = -( het_rates(95) )

         mat(246) = -( het_rates(96) )

         mat(160) = -( rxt(85) + het_rates(97) )

         mat(61) = -( het_rates(98) )

         mat(540) = -( het_rates(100) )

         mat(208) = -( rxt(86) + het_rates(101) )

         mat(254) = -( rxt(81) + het_rates(86) )
         mat(162) = .900_r8*rxt(85)

         mat(375) = -( rxt(82) + het_rates(53) )
         mat(255) = .130_r8*rxt(81)
         mat(163) = .450_r8*rxt(85)

         mat(697) = -( het_rates(88) )

         mat(740) = -( rxt(68) + het_rates(77) )
         mat(276) = .402_r8*rxt(77)
         mat(212) = rxt(86)

         mat(637) = -( rxt(66) + rxt(67) + het_rates(78) )
         mat(273) = .288_r8*rxt(77)
         mat(211) = rxt(86)

         mat(721) = -( het_rates(79) )

         mat(120) = -( het_rates(80) )

         mat(760) = -( het_rates(76) )
         mat(290) = rxt(65) + rxt(357)
         mat(640) = .660_r8*rxt(66)

         mat(462) = -( het_rates(46) )
         mat(186) = rxt(80)

         mat(125) = -( rxt(69) + het_rates(47) )

         mat(303) = -( het_rates(99) )

         mat(29) = -( het_rates(60) )

         mat(517) = -( het_rates(61) )

         mat(166) = -( rxt(71) + het_rates(62) )

         mat(369) = -( rxt(73) + het_rates(63) )
         mat(167) = .820_r8*rxt(71)
         mat(327) = .250_r8*rxt(83)
         mat(209) = .100_r8*rxt(86)

         mat(172) = -( rxt(72) + het_rates(69) )

         mat(268) = -( het_rates(18) )

         mat(75) = -( het_rates(51) )

         mat(510) = -( rxt(79) + het_rates(52) )

         mat(616) = -( rxt(78) + het_rates(66) )

         mat(388) = -( het_rates(55) )

         mat(189) = -( rxt(292) + het_rates(56) )
         mat(42) = rxt(70)

         mat(41) = -( rxt(70) + het_rates(57) )

         mat(139) = -( het_rates(87) )

         mat(625) = -( het_rates(67) )

         mat(662) = -( rxt(74) + het_rates(68) )
         mat(257) = .180_r8*rxt(81)
         mat(164) = .450_r8*rxt(85)

         mat(572) = -( het_rates(89) )

         mat(530) = -( rxt(76) + het_rates(90) )

         mat(677) = -( het_rates(91) )

         mat(130) = -( rxt(75) + het_rates(92) )

         mat(272) = -( rxt(77) + het_rates(93) )

         mat(90) = -( het_rates(128) )

         mat(241) = -( het_rates(129) )

         mat(178) = -( rxt(280) + het_rates(130) )

         mat(44) = -( rxt(55) + het_rates(131) )
         mat(958) = rxt(139)*y(105) + rxt(140)*y(106) + 2.000_r8*rxt(141)*y(114) &
                      + 2.000_r8*rxt(142)*y(115) + rxt(143)*y(107) + rxt(145)*y(113) &
                      + rxt(148)*y(111) + rxt(149)*y(110) + rxt(150)*y(116) &
                      + 2.000_r8*rxt(151)*y(117)
         mat(1302) = rxt(257)*y(107) + rxt(261)*y(113)

         mat(65) = -( rxt(56) + het_rates(132) )
         mat(961) = rxt(138)*y(104) + rxt(140)*y(106) + rxt(144)*y(112)
         mat(1305) = rxt(260)*y(112)

         mat(72) = -( rxt(57) + het_rates(133) )
         mat(432) = rxt(252)*y(15)

         mat(433) = -( rxt(252)*y(15) + het_rates(134) )
         mat(45) = 2.000_r8*rxt(55)
         mat(66) = rxt(56)
         mat(73) = rxt(57)
         mat(962) = rxt(142)*y(115) + rxt(149)*y(110)

         mat(552) = -( rxt(88) + het_rates(135) )
         mat(81) = rxt(89)

         mat(96) = -( het_rates(136) )

         mat(142) = -( rxt(90) + het_rates(137) )

         mat(379) = -( het_rates(138) )
         mat(143) = rxt(90)
         mat(803) = rxt(91)

         mat(805) = -( rxt(91) + het_rates(139) )
         mat(553) = rxt(88)

         mat(80) = -( rxt(89) + het_rates(140) )
         mat(48) = rxt(87)

         mat(47) = -( rxt(87) + het_rates(141) )

         mat(1) = -( het_rates(142) )

         mat(2) = -( het_rates(143) )

         mat(3) = -( het_rates(144) )

         mat(4) = -( het_rates(145) )

         mat(5) = -( het_rates(146) )

         mat(6) = -( het_rates(147) )

         mat(7) = -( het_rates(148) )

         mat(8) = -( het_rates(149) )

         mat(9) = -( het_rates(150) )

         mat(10) = -( het_rates(151) )

         mat(11) = -( het_rates(152) )

         mat(12) = -( het_rates(153) )

         mat(13) = -( het_rates(154) )

         mat(14) = -( het_rates(155) )

         mat(15) = -( het_rates(156) )

         mat(16) = -( het_rates(157) )


      end subroutine linmat02

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
      call linmat02( mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix
