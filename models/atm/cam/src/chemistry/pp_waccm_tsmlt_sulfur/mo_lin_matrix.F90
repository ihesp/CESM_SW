




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

         mat(1198) = -( rxt(3) + rxt(4) + het_rates(1) )

         mat(945) = -( rxt(92) + rxt(93) + rxt(94) + rxt(105) + rxt(106) + rxt(107) &
                 + het_rates(2) )
         mat(906) = rxt(1) + 2.000_r8*rxt(2) + rxt(98) + rxt(99) + rxt(100) &
                      + 2.000_r8*rxt(103) + rxt(110) + rxt(111) + rxt(112) &
                      + 2.000_r8*rxt(115)
         mat(1191) = rxt(4)
         mat(1246) = rxt(6)
         mat(1283) = rxt(8)
         mat(105) = rxt(10)
         mat(1425) = rxt(12)
         mat(1473) = rxt(21)
         mat(1031) = rxt(24)
         mat(139) = rxt(25)
         mat(969) = rxt(32)
         mat(556) = rxt(88)
         mat(84) = rxt(89)
         mat(810) = rxt(91)
         mat(1057) = rxt(131)

         mat(1061) = -( rxt(131) + rxt(135)*y(7) + rxt(136)*y(7) + rxt(138)*y(104) &
                      + rxt(139)*y(105) + rxt(140)*y(106) + rxt(141)*y(114) &
                      + rxt(142)*y(115) + rxt(143)*y(107) + rxt(144)*y(112) &
                      + rxt(145)*y(113) + rxt(146)*y(108) + rxt(147)*y(103) &
                      + rxt(148)*y(111) + rxt(149)*y(110) + rxt(150)*y(116) &
                      + rxt(151)*y(117) + rxt(152)*y(118) + rxt(153)*y(119) &
                      + rxt(156)*y(15) + rxt(157)*y(15) + rxt(158)*y(15) + het_rates(3) )
         mat(910) = rxt(1)
         mat(1195) = rxt(3)
         mat(1477) = rxt(20)

         mat(905) = -( rxt(1) + rxt(2) + rxt(96) + rxt(98) + rxt(99) + rxt(100) + rxt(103) &
                      + rxt(108) + rxt(110) + rxt(111) + rxt(112) + rxt(115) &
                 + het_rates(4) )
         mat(1190) = rxt(4)
         mat(1424) = rxt(13)
         mat(56) = rxt(126)
         mat(53) = rxt(129) + rxt(130)
         mat(1056) = rxt(136)*y(7)

         mat(55) = -( rxt(123) + rxt(126) + rxt(125)*y(120) + het_rates(5) )

         mat(52) = -( rxt(129) + rxt(130) + het_rates(6) )
         mat(1161) = rxt(3)
         mat(54) = rxt(123) + rxt(125)*y(120)

         mat(652) = -( het_rates(21) )
         mat(1492) = rxt(18)
         mat(1467) = rxt(20)
         mat(1052) = rxt(158)*y(15)

         mat(604) = -( het_rates(20) )
         mat(1491) = rxt(17) + rxt(18)
         mat(608) = rxt(61)
         mat(638) = 1.340_r8*rxt(67)
         mat(737) = .700_r8*rxt(68)
         mat(663) = rxt(74)
         mat(533) = rxt(76)
         mat(513) = rxt(79)
         mat(258) = .450_r8*rxt(81)
         mat(378) = 2.000_r8*rxt(82)
         mat(147) = rxt(90)
         mat(996) = rxt(254)*y(102)
         mat(302) = rxt(457)*y(120)

         mat(478) = -( rxt(95) + het_rates(8) )
         mat(1225) = rxt(6)
         mat(301) = rxt(454)

         mat(1254) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(1291) = rxt(8) + .500_r8*rxt(399)
         mat(106) = rxt(10)
         mat(1433) = rxt(13)
         mat(414) = rxt(464)
         mat(1065) = 2.000_r8*rxt(135)*y(7)

         mat(1292) = -( rxt(8) + rxt(399) + het_rates(10) )
         mat(107) = rxt(9) + rxt(197)
         mat(1456) = rxt(11)
         mat(1434) = rxt(12)
         mat(220) = rxt(15) + rxt(206)
         mat(567) = rxt(30)
         mat(287) = rxt(36)
         mat(199) = .600_r8*rxt(64) + rxt(311)
         mat(294) = rxt(65) + rxt(357)
         mat(536) = rxt(76)

         mat(1391) = -( rxt(255)*y(102) + rxt(256)*y(109) + rxt(257)*y(107) &
                      + rxt(258)*y(103) + rxt(260)*y(112) + rxt(261)*y(113) &
                      + rxt(262)*y(119) + rxt(263)*y(118) + rxt(266)*y(15) + het_rates(23) &
       )
         mat(1457) = rxt(11)
         mat(221) = rxt(14)
         mat(159) = rxt(16)
         mat(1483) = rxt(19)
         mat(319) = 2.000_r8*rxt(22)
         mat(493) = rxt(27)
         mat(405) = rxt(33)
         mat(267) = rxt(62)
         mat(232) = rxt(63)
         mat(131) = rxt(69)
         mat(45) = rxt(70)
         mat(172) = rxt(71)
         mat(177) = rxt(72)
         mat(134) = rxt(75)
         mat(334) = rxt(83)
         mat(121) = rxt(84)
         mat(167) = rxt(85)
         mat(216) = rxt(86)
         mat(1293) = .500_r8*rxt(399)
         mat(1067) = rxt(156)*y(15)

         mat(1436) = -( rxt(12) + rxt(13) + rxt(398) + het_rates(11) )
         mat(108) = rxt(9) + rxt(10) + rxt(197)
         mat(222) = rxt(14)
         mat(569) = rxt(29)
         mat(288) = rxt(35)
         mat(201) = .400_r8*rxt(64)

         mat(1459) = -( rxt(11) + het_rates(12) )
         mat(109) = 2.000_r8*rxt(397) + 2.000_r8*rxt(436) + 2.000_r8*rxt(442) &
                      + 2.000_r8*rxt(447)
         mat(1437) = rxt(398)
         mat(1295) = .500_r8*rxt(399)
         mat(570) = rxt(437) + rxt(443) + rxt(448)
         mat(289) = rxt(438) + rxt(446) + rxt(449)

         mat(217) = -( rxt(14) + rxt(15) + rxt(206) + het_rates(13) )

         mat(104) = -( rxt(9) + rxt(10) + rxt(197) + rxt(397) + rxt(436) + rxt(442) &
                      + rxt(447) + het_rates(14) )

         mat(874) = -( het_rates(16) )
         mat(611) = rxt(61)
         mat(231) = rxt(63)
         mat(198) = .400_r8*rxt(64)
         mat(745) = .300_r8*rxt(68)
         mat(374) = rxt(73)
         mat(1055) = rxt(156)*y(15)
         mat(1002) = rxt(213)*y(15)
         mat(437) = rxt(252)*y(15)
         mat(1379) = rxt(266)*y(15)

         mat(156) = -( rxt(16) + het_rates(17) )

         mat(59) = -( het_rates(42) )

         mat(19) = -( het_rates(43) )

         mat(1511) = -( rxt(17) + rxt(18) + het_rates(19) )
         mat(161) = rxt(16)
         mat(269) = rxt(62)
         mat(649) = 1.340_r8*rxt(66)
         mat(179) = rxt(72)
         mat(539) = rxt(76)
         mat(281) = .690_r8*rxt(77)
         mat(623) = rxt(78)
         mat(516) = rxt(79)
         mat(335) = .100_r8*rxt(83)
         mat(185) = rxt(280)
         mat(195) = 2.000_r8*rxt(292)
         mat(1071) = rxt(157)*y(15) + rxt(158)*y(15)

         mat(1152) = -( het_rates(22) )
         mat(158) = rxt(16)
         mat(1503) = 2.000_r8*rxt(17)
         mat(1479) = rxt(19) + 2.000_r8*rxt(21)
         mat(833) = rxt(28)
         mat(458) = rxt(34)
         mat(76) = rxt(57)
         mat(1063) = rxt(157)*y(15)

         mat(1131) = -( rxt(402) + het_rates(24) )
         mat(219) = rxt(15) + rxt(206)
         mat(612) = rxt(61)
         mat(266) = rxt(62)
         mat(644) = 1.340_r8*rxt(66) + .660_r8*rxt(67)
         mat(130) = rxt(69)
         mat(171) = rxt(71)
         mat(666) = rxt(74)
         mat(535) = rxt(76)
         mat(279) = rxt(77)
         mat(621) = rxt(78)
         mat(514) = 2.000_r8*rxt(79)
         mat(261) = .560_r8*rxt(81)
         mat(379) = 2.000_r8*rxt(82)
         mat(333) = .900_r8*rxt(83)
         mat(215) = rxt(86)
         mat(182) = rxt(280)
         mat(194) = rxt(292)
         mat(1062) = rxt(157)*y(15)
         mat(1009) = rxt(254)*y(102) + rxt(259)*y(103)
         mat(1386) = rxt(255)*y(102) + rxt(258)*y(103)

         mat(314) = -( rxt(22) + het_rates(25) )
         mat(1091) = .500_r8*rxt(402)

         mat(1486) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(128) )
         mat(51) = rxt(87)
         mat(1394) = rxt(255)*y(102) + rxt(256)*y(109) + rxt(257)*y(107) + rxt(258)*y(103) &
                      + rxt(262)*y(119) + rxt(266)*y(15)

         mat(1006) = -( rxt(213)*y(15) + rxt(254)*y(102) + rxt(259)*y(103) &
                      + rxt(264)*y(119) + rxt(265)*y(118) + het_rates(28) )
         mat(58) = 2.000_r8*rxt(23)
         mat(1033) = rxt(24)
         mat(23) = 2.000_r8*rxt(26)
         mat(491) = rxt(27)
         mat(830) = rxt(28)
         mat(565) = rxt(29)
         mat(73) = rxt(31)
         mat(69) = rxt(56)
         mat(1059) = 2.000_r8*rxt(138)*y(104) + 2.000_r8*rxt(139)*y(105) &
                      + 2.000_r8*rxt(140)*y(106) + 2.000_r8*rxt(141)*y(114) &
                      + rxt(142)*y(115) + rxt(143)*y(107) + rxt(144)*y(112) &
                      + rxt(145)*y(113) + 4.000_r8*rxt(146)*y(108) + rxt(148)*y(111)
         mat(1383) = rxt(255)*y(102) + 3.000_r8*rxt(256)*y(109) + rxt(257)*y(107) &
                      + rxt(260)*y(112) + rxt(261)*y(113)

         mat(57) = -( rxt(23) + het_rates(29) )

         mat(1034) = -( rxt(24) + het_rates(30) )
         mat(140) = rxt(25)
         mat(566) = rxt(30)
         mat(24) = 2.000_r8*rxt(225)

         mat(136) = -( rxt(25) + het_rates(31) )

         mat(22) = -( rxt(26) + rxt(225) + het_rates(32) )


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

         mat(826) = -( rxt(28) + het_rates(33) )
         mat(1000) = rxt(213)*y(15) + 2.000_r8*rxt(254)*y(102) + rxt(259)*y(103) &
                      + rxt(264)*y(119) + rxt(265)*y(118)

         mat(488) = -( rxt(27) + het_rates(34) )
         mat(561) = rxt(437) + rxt(443) + rxt(448)

         mat(562) = -( rxt(29) + rxt(30) + rxt(437) + rxt(443) + rxt(448) + het_rates(35) &
       )

         mat(71) = -( rxt(31) + het_rates(36) )

         mat(841) = -( het_rates(37) )
         mat(72) = rxt(31)
         mat(967) = rxt(32)
         mat(401) = rxt(33)
         mat(455) = rxt(34)
         mat(284) = rxt(35)
         mat(1054) = rxt(147)*y(103) + rxt(148)*y(111) + rxt(149)*y(110) &
                      + 2.000_r8*rxt(150)*y(116) + 2.000_r8*rxt(151)*y(117) &
                      + 3.000_r8*rxt(152)*y(118) + 2.000_r8*rxt(153)*y(119)
         mat(1378) = rxt(258)*y(103) + 2.000_r8*rxt(262)*y(119) + 3.000_r8*rxt(263)*y(118)
         mat(1001) = rxt(259)*y(103) + 2.000_r8*rxt(264)*y(119) + 3.000_r8*rxt(265)*y(118)

         mat(970) = -( rxt(32) + het_rates(38) )
         mat(286) = rxt(36)

         mat(454) = -( rxt(34) + het_rates(39) )

         mat(399) = -( rxt(33) + het_rates(40) )
         mat(283) = rxt(438) + rxt(446) + rxt(449)

         mat(282) = -( rxt(35) + rxt(36) + rxt(438) + rxt(446) + rxt(449) + het_rates(41) &
       )

         mat(346) = -( het_rates(121) )

         mat(407) = -( rxt(464) + het_rates(122) )
         mat(896) = rxt(96) + rxt(108)
         mat(299) = rxt(457)*y(120)

         mat(203) = -( het_rates(123) )
         mat(473) = rxt(95)

         mat(298) = -( rxt(454) + rxt(457)*y(120) + het_rates(124) )
         mat(925) = rxt(92) + rxt(93) + rxt(94) + rxt(105) + rxt(106) + rxt(107)
         mat(893) = rxt(98) + rxt(99) + rxt(100) + rxt(110) + rxt(111) + rxt(112)

         mat(416) = -( het_rates(125) )
         mat(1221) = rxt(7)
         mat(300) = rxt(454)
         mat(408) = rxt(464)

         mat(224) = -( het_rates(127) )

         mat(427) = -( het_rates(126) )
         mat(1222) = rxt(7)
         mat(932) = rxt(92) + rxt(93) + rxt(94) + rxt(105) + rxt(106) + rxt(107)
         mat(477) = rxt(95)
         mat(898) = rxt(96) + rxt(98) + rxt(99) + rxt(100) + rxt(108) + rxt(110) &
                      + rxt(111) + rxt(112)

         mat(589) = -( het_rates(59) )
         mat(736) = .700_r8*rxt(68)

         mat(496) = -( het_rates(83) )

         mat(444) = -( het_rates(64) )

         mat(609) = -( rxt(61) + het_rates(50) )
         mat(264) = rxt(62)
         mat(129) = rxt(69)
         mat(331) = .400_r8*rxt(83)
         mat(119) = rxt(84)

         mat(321) = -( het_rates(49) )

         mat(262) = -( rxt(62) + het_rates(65) )

         mat(791) = -( het_rates(48) )
         mat(197) = .600_r8*rxt(64) + rxt(311)
         mat(643) = 1.340_r8*rxt(66)
         mat(744) = .300_r8*rxt(68)
         mat(176) = rxt(72)
         mat(373) = rxt(73)
         mat(665) = rxt(74)
         mat(620) = rxt(78)
         mat(189) = rxt(80)
         mat(260) = .130_r8*rxt(81)
         mat(120) = rxt(84)

         mat(229) = -( rxt(63) + het_rates(54) )

         mat(196) = -( rxt(64) + rxt(311) + het_rates(58) )

         mat(152) = -( het_rates(82) )

         mat(86) = -( het_rates(45) )

         mat(235) = -( het_rates(44) )

         mat(25) = -( het_rates(71) )

         mat(290) = -( rxt(65) + rxt(357) + het_rates(81) )

         mat(28) = -( het_rates(70) )

         mat(110) = -( het_rates(73) )

         mat(360) = -( het_rates(84) )

         mat(326) = -( rxt(83) + het_rates(85) )

         mat(186) = -( rxt(80) + het_rates(72) )
         mat(325) = .800_r8*rxt(83)

         mat(337) = -( het_rates(74) )

         mat(117) = -( rxt(84) + het_rates(75) )

         mat(35) = -( het_rates(94) )

         mat(40) = -( het_rates(95) )

         mat(248) = -( het_rates(96) )

         mat(162) = -( rxt(85) + het_rates(97) )

         mat(63) = -( het_rates(98) )

         mat(542) = -( het_rates(100) )

         mat(210) = -( rxt(86) + het_rates(101) )

         mat(256) = -( rxt(81) + het_rates(86) )
         mat(164) = .900_r8*rxt(85)

         mat(377) = -( rxt(82) + het_rates(53) )
         mat(257) = .130_r8*rxt(81)
         mat(165) = .450_r8*rxt(85)

         mat(699) = -( het_rates(88) )

         mat(742) = -( rxt(68) + het_rates(77) )
         mat(278) = .402_r8*rxt(77)
         mat(214) = rxt(86)

         mat(639) = -( rxt(66) + rxt(67) + het_rates(78) )
         mat(275) = .288_r8*rxt(77)
         mat(213) = rxt(86)

         mat(723) = -( het_rates(79) )

         mat(122) = -( het_rates(80) )

         mat(762) = -( het_rates(76) )
         mat(292) = rxt(65) + rxt(357)
         mat(642) = .660_r8*rxt(66)

         mat(464) = -( het_rates(46) )
         mat(188) = rxt(80)

         mat(127) = -( rxt(69) + het_rates(47) )

         mat(305) = -( het_rates(99) )

         mat(31) = -( het_rates(60) )

         mat(519) = -( het_rates(61) )

         mat(168) = -( rxt(71) + het_rates(62) )

         mat(371) = -( rxt(73) + het_rates(63) )
         mat(169) = .820_r8*rxt(71)
         mat(329) = .250_r8*rxt(83)
         mat(211) = .100_r8*rxt(86)

         mat(174) = -( rxt(72) + het_rates(69) )

         mat(270) = -( het_rates(18) )

         mat(77) = -( het_rates(51) )

         mat(512) = -( rxt(79) + het_rates(52) )

         mat(618) = -( rxt(78) + het_rates(66) )

         mat(390) = -( het_rates(55) )

         mat(191) = -( rxt(292) + het_rates(56) )
         mat(44) = rxt(70)

         mat(43) = -( rxt(70) + het_rates(57) )

         mat(141) = -( het_rates(87) )

         mat(627) = -( het_rates(67) )

         mat(664) = -( rxt(74) + het_rates(68) )
         mat(259) = .180_r8*rxt(81)
         mat(166) = .450_r8*rxt(85)

         mat(574) = -( het_rates(89) )

         mat(532) = -( rxt(76) + het_rates(90) )

         mat(679) = -( het_rates(91) )

         mat(132) = -( rxt(75) + het_rates(92) )

         mat(274) = -( rxt(77) + het_rates(93) )

         mat(92) = -( het_rates(129) )

         mat(243) = -( het_rates(130) )

         mat(180) = -( rxt(280) + het_rates(131) )

         mat(46) = -( rxt(55) + het_rates(132) )
         mat(1046) = rxt(139)*y(105) + rxt(140)*y(106) + 2.000_r8*rxt(141)*y(114) &
                      + 2.000_r8*rxt(142)*y(115) + rxt(143)*y(107) + rxt(145)*y(113) &
                      + rxt(148)*y(111) + rxt(149)*y(110) + rxt(150)*y(116) &
                      + 2.000_r8*rxt(151)*y(117)
         mat(1304) = rxt(257)*y(107) + rxt(261)*y(113)

         mat(67) = -( rxt(56) + het_rates(133) )
         mat(1049) = rxt(138)*y(104) + rxt(140)*y(106) + rxt(144)*y(112)
         mat(1307) = rxt(260)*y(112)

         mat(74) = -( rxt(57) + het_rates(134) )
         mat(434) = rxt(252)*y(15)

         mat(435) = -( rxt(252)*y(15) + het_rates(135) )
         mat(47) = 2.000_r8*rxt(55)
         mat(68) = rxt(56)
         mat(75) = rxt(57)
         mat(1050) = rxt(142)*y(115) + rxt(149)*y(110)

         mat(554) = -( rxt(88) + het_rates(141) )
         mat(83) = rxt(89)

         mat(98) = -( het_rates(142) )

         mat(1) = -( rxt(407) + het_rates(143) )

         mat(2) = -( rxt(409) + het_rates(144) )

         mat(3) = -( rxt(410) + het_rates(145) )

         mat(4) = -( rxt(408) + het_rates(140) )

         mat(144) = -( rxt(90) + het_rates(146) )

         mat(381) = -( het_rates(147) )
         mat(145) = rxt(90)
         mat(805) = rxt(91)

         mat(807) = -( rxt(91) + het_rates(148) )
         mat(555) = rxt(88)

         mat(82) = -( rxt(89) + het_rates(149) )
         mat(50) = rxt(87)

         mat(49) = -( rxt(87) + het_rates(150) )

         mat(5) = -( rxt(400) + rxt(403) + het_rates(136) )

         mat(7) = -( rxt(404) + het_rates(137) )
         mat(6) = rxt(400)

         mat(8) = -( rxt(401) + rxt(405) + het_rates(138) )

         mat(10) = -( rxt(406) + het_rates(139) )
         mat(9) = rxt(401)

         mat(11) = -( rxt(411) + het_rates(151) )

         mat(12) = -( rxt(412) + het_rates(152) )

         mat(13) = -( rxt(413) + het_rates(153) )

         mat(14) = -( rxt(414) + het_rates(154) )

         mat(15) = -( rxt(415) + het_rates(155) )

         mat(16) = -( rxt(416) + het_rates(156) )

         mat(17) = -( rxt(417) + het_rates(157) )

         mat(18) = -( rxt(418) + het_rates(158) )


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
