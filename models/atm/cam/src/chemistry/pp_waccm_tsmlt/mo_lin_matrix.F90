




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

         mat(922) = -( rxt(3) + rxt(4) + het_rates(1) )

         mat(960) = -( rxt(89) + rxt(90) + rxt(91) + rxt(102) + rxt(103) + rxt(104) &
                 + het_rates(2) )
         mat(811) = rxt(1) + 2.000_r8*rxt(2) + rxt(95) + rxt(96) + rxt(97) &
                      + 2.000_r8*rxt(100) + rxt(107) + rxt(108) + rxt(109) &
                      + 2.000_r8*rxt(112)
         mat(923) = rxt(4)
         mat(1115) = rxt(6)
         mat(1150) = rxt(8)
         mat(107) = rxt(10)
         mat(1287) = rxt(12)
         mat(824) = rxt(21)
         mat(859) = rxt(24)
         mat(43) = rxt(25)
         mat(788) = rxt(32)
         mat(984) = rxt(128)

         mat(985) = -( rxt(128) + rxt(132)*y(7) + rxt(133)*y(7) + rxt(135)*y(104) &
                      + rxt(136)*y(105) + rxt(137)*y(106) + rxt(138)*y(114) &
                      + rxt(139)*y(115) + rxt(140)*y(107) + rxt(141)*y(112) &
                      + rxt(142)*y(113) + rxt(143)*y(108) + rxt(144)*y(103) &
                      + rxt(145)*y(111) + rxt(146)*y(110) + rxt(147)*y(116) &
                      + rxt(148)*y(117) + rxt(149)*y(118) + rxt(150)*y(119) &
                      + rxt(153)*y(15) + rxt(154)*y(15) + rxt(155)*y(15) + het_rates(3) )
         mat(812) = rxt(1)
         mat(924) = rxt(3)
         mat(825) = rxt(20)

         mat(809) = -( rxt(1) + rxt(2) + rxt(93) + rxt(95) + rxt(96) + rxt(97) + rxt(100) &
                      + rxt(105) + rxt(107) + rxt(108) + rxt(109) + rxt(112) &
                 + het_rates(4) )
         mat(917) = rxt(4)
         mat(1283) = rxt(13)
         mat(62) = rxt(123)
         mat(59) = rxt(126) + rxt(127)
         mat(979) = rxt(133)*y(7)

         mat(61) = -( rxt(120) + rxt(123) + rxt(122)*y(120) + het_rates(5) )

         mat(58) = -( rxt(126) + rxt(127) + het_rates(6) )
         mat(892) = rxt(3)
         mat(60) = rxt(120) + rxt(122)*y(120)

         mat(626) = -( het_rates(21) )
         mat(1303) = rxt(18)
         mat(821) = rxt(20)
         mat(978) = rxt(155)*y(15)

         mat(578) = -( het_rates(20) )
         mat(1302) = rxt(17) + rxt(18)
         mat(582) = rxt(61)
         mat(612) = 1.340_r8*rxt(67)
         mat(711) = .700_r8*rxt(68)
         mat(637) = rxt(74)
         mat(514) = rxt(76)
         mat(494) = rxt(79)
         mat(242) = .450_r8*rxt(81)
         mat(345) = 2.000_r8*rxt(82)
         mat(1008) = rxt(251)*y(102)
         mat(303) = rxt(443)*y(120)

         mat(459) = -( rxt(92) + het_rates(8) )
         mat(1090) = rxt(6)
         mat(302) = rxt(440)

         mat(1120) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(1155) = rxt(8) + .500_r8*rxt(396)
         mat(108) = rxt(10)
         mat(1292) = rxt(13)
         mat(404) = rxt(450)
         mat(989) = 2.000_r8*rxt(132)*y(7)

         mat(1156) = -( rxt(8) + rxt(396) + het_rates(10) )
         mat(109) = rxt(9) + rxt(194)
         mat(842) = rxt(11)
         mat(1293) = rxt(12)
         mat(208) = rxt(15) + rxt(203)
         mat(543) = rxt(30)
         mat(252) = rxt(36)
         mat(200) = .600_r8*rxt(64) + rxt(308)
         mat(226) = rxt(65) + rxt(354)
         mat(516) = rxt(76)

         mat(1253) = -( rxt(252)*y(102) + rxt(253)*y(109) + rxt(254)*y(107) &
                      + rxt(255)*y(103) + rxt(257)*y(112) + rxt(258)*y(113) &
                      + rxt(259)*y(119) + rxt(260)*y(118) + rxt(263)*y(15) + het_rates(23) &
       )
         mat(843) = rxt(11)
         mat(209) = rxt(14)
         mat(171) = rxt(16)
         mat(829) = rxt(19)
         mat(311) = 2.000_r8*rxt(22)
         mat(491) = rxt(27)
         mat(377) = rxt(33)
         mat(283) = rxt(62)
         mat(215) = rxt(63)
         mat(131) = rxt(69)
         mat(54) = rxt(70)
         mat(154) = rxt(71)
         mat(160) = rxt(72)
         mat(136) = rxt(75)
         mat(321) = rxt(83)
         mat(122) = rxt(84)
         mat(149) = rxt(85)
         mat(196) = rxt(86)
         mat(1157) = .500_r8*rxt(396)
         mat(991) = rxt(153)*y(15)

         mat(1295) = -( rxt(12) + rxt(13) + rxt(395) + het_rates(11) )
         mat(110) = rxt(9) + rxt(10) + rxt(194)
         mat(210) = rxt(14)
         mat(545) = rxt(29)
         mat(253) = rxt(35)
         mat(202) = .400_r8*rxt(64)

         mat(836) = -( rxt(11) + het_rates(12) )
         mat(106) = 2.000_r8*rxt(394) + 2.000_r8*rxt(422) + 2.000_r8*rxt(428) &
                      + 2.000_r8*rxt(433)
         mat(1285) = rxt(395)
         mat(1146) = .500_r8*rxt(396)
         mat(538) = rxt(423) + rxt(429) + rxt(434)
         mat(249) = rxt(424) + rxt(432) + rxt(435)

         mat(205) = -( rxt(14) + rxt(15) + rxt(203) + het_rates(13) )

         mat(105) = -( rxt(9) + rxt(10) + rxt(194) + rxt(394) + rxt(422) + rxt(428) &
                      + rxt(433) + het_rates(14) )

         mat(1428) = -( het_rates(16) )
         mat(591) = rxt(61)
         mat(217) = rxt(63)
         mat(204) = .400_r8*rxt(64)
         mat(730) = .300_r8*rxt(68)
         mat(368) = rxt(73)
         mat(995) = rxt(153)*y(15)
         mat(1029) = rxt(210)*y(15)
         mat(396) = rxt(249)*y(15)
         mat(1257) = rxt(263)*y(15)

         mat(168) = -( rxt(16) + het_rates(17) )

         mat(65) = -( het_rates(42) )

         mat(19) = -( het_rates(43) )

         mat(1317) = -( rxt(17) + rxt(18) + het_rates(19) )
         mat(172) = rxt(16)
         mat(284) = rxt(62)
         mat(622) = 1.340_r8*rxt(66)
         mat(161) = rxt(72)
         mat(519) = rxt(76)
         mat(276) = .690_r8*rxt(77)
         mat(596) = rxt(78)
         mat(496) = rxt(79)
         mat(322) = .100_r8*rxt(83)
         mat(166) = rxt(277)
         mat(182) = 2.000_r8*rxt(289)
         mat(993) = rxt(154)*y(15) + rxt(155)*y(15)

         mat(1038) = -( het_rates(22) )
         mat(170) = rxt(16)
         mat(1311) = 2.000_r8*rxt(17)
         mat(827) = rxt(19) + 2.000_r8*rxt(21)
         mat(1062) = rxt(28)
         mat(440) = rxt(34)
         mat(82) = rxt(57)
         mat(987) = rxt(154)*y(15)

         mat(1387) = -( rxt(404) + het_rates(24) )
         mat(211) = rxt(15) + rxt(203)
         mat(590) = rxt(61)
         mat(285) = rxt(62)
         mat(623) = 1.340_r8*rxt(66) + .660_r8*rxt(67)
         mat(132) = rxt(69)
         mat(155) = rxt(71)
         mat(645) = rxt(74)
         mat(520) = rxt(76)
         mat(277) = rxt(77)
         mat(597) = rxt(78)
         mat(497) = 2.000_r8*rxt(79)
         mat(245) = .560_r8*rxt(81)
         mat(347) = 2.000_r8*rxt(82)
         mat(323) = .900_r8*rxt(83)
         mat(197) = rxt(86)
         mat(167) = rxt(277)
         mat(183) = rxt(289)
         mat(994) = rxt(154)*y(15)
         mat(1028) = rxt(251)*y(102) + rxt(256)*y(103)
         mat(1256) = rxt(252)*y(102) + rxt(255)*y(103)

         mat(306) = -( rxt(22) + het_rates(25) )
         mat(1340) = .500_r8*rxt(404)

         mat(822) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(128) )
         mat(1241) = rxt(252)*y(102) + rxt(253)*y(109) + rxt(254)*y(107) + rxt(255)*y(103) &
                      + rxt(259)*y(119) + rxt(263)*y(15)

         mat(1020) = -( rxt(210)*y(15) + rxt(251)*y(102) + rxt(256)*y(103) &
                      + rxt(261)*y(119) + rxt(262)*y(118) + het_rates(28) )
         mat(64) = 2.000_r8*rxt(23)
         mat(861) = rxt(24)
         mat(24) = 2.000_r8*rxt(26)
         mat(489) = rxt(27)
         mat(1061) = rxt(28)
         mat(541) = rxt(29)
         mat(79) = rxt(31)
         mat(76) = rxt(56)
         mat(986) = 2.000_r8*rxt(135)*y(104) + 2.000_r8*rxt(136)*y(105) &
                      + 2.000_r8*rxt(137)*y(106) + 2.000_r8*rxt(138)*y(114) &
                      + rxt(139)*y(115) + rxt(140)*y(107) + rxt(141)*y(112) &
                      + rxt(142)*y(113) + 4.000_r8*rxt(143)*y(108) + rxt(145)*y(111)
         mat(1248) = rxt(252)*y(102) + 3.000_r8*rxt(253)*y(109) + rxt(254)*y(107) &
                      + rxt(257)*y(112) + rxt(258)*y(113)

         mat(63) = -( rxt(23) + het_rates(29) )

         mat(856) = -( rxt(24) + het_rates(30) )
         mat(42) = rxt(25)
         mat(539) = rxt(30)
         mat(23) = 2.000_r8*rxt(222)

         mat(41) = -( rxt(25) + het_rates(31) )

         mat(22) = -( rxt(26) + rxt(222) + het_rates(32) )

         mat(1063) = -( rxt(28) + het_rates(33) )
         mat(1022) = rxt(210)*y(15) + 2.000_r8*rxt(251)*y(102) + rxt(256)*y(103) &
                      + rxt(261)*y(119) + rxt(262)*y(118)

         mat(485) = -( rxt(27) + het_rates(34) )
         mat(535) = rxt(423) + rxt(429) + rxt(434)


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

         mat(536) = -( rxt(29) + rxt(30) + rxt(423) + rxt(429) + rxt(434) + het_rates(35) &
       )

         mat(77) = -( rxt(31) + het_rates(36) )

         mat(878) = -( het_rates(37) )
         mat(78) = rxt(31)
         mat(787) = rxt(32)
         mat(373) = rxt(33)
         mat(437) = rxt(34)
         mat(250) = rxt(35)
         mat(982) = rxt(144)*y(103) + rxt(145)*y(111) + rxt(146)*y(110) &
                      + 2.000_r8*rxt(147)*y(116) + 2.000_r8*rxt(148)*y(117) &
                      + 3.000_r8*rxt(149)*y(118) + 2.000_r8*rxt(150)*y(119)
         mat(1244) = rxt(255)*y(103) + 2.000_r8*rxt(259)*y(119) + 3.000_r8*rxt(260)*y(118)
         mat(1016) = rxt(256)*y(103) + 2.000_r8*rxt(261)*y(119) + 3.000_r8*rxt(262)*y(118)

         mat(782) = -( rxt(32) + het_rates(38) )
         mat(248) = rxt(36)

         mat(435) = -( rxt(34) + het_rates(39) )

         mat(370) = -( rxt(33) + het_rates(40) )
         mat(247) = rxt(424) + rxt(432) + rxt(435)

         mat(246) = -( rxt(35) + rxt(36) + rxt(424) + rxt(432) + rxt(435) + het_rates(41) &
       )

         mat(334) = -( het_rates(121) )

         mat(397) = -( rxt(450) + het_rates(122) )
         mat(803) = rxt(93) + rxt(105)
         mat(300) = rxt(443)*y(120)

         mat(184) = -( het_rates(123) )
         mat(454) = rxt(92)

         mat(299) = -( rxt(440) + rxt(443)*y(120) + het_rates(124) )
         mat(940) = rxt(89) + rxt(90) + rxt(91) + rxt(102) + rxt(103) + rxt(104)
         mat(801) = rxt(95) + rxt(96) + rxt(97) + rxt(107) + rxt(108) + rxt(109)

         mat(406) = -( het_rates(125) )
         mat(1086) = rxt(7)
         mat(301) = rxt(440)
         mat(398) = rxt(450)

         mat(218) = -( het_rates(127) )

         mat(417) = -( het_rates(126) )
         mat(1087) = rxt(7)
         mat(946) = rxt(89) + rxt(90) + rxt(91) + rxt(102) + rxt(103) + rxt(104)
         mat(458) = rxt(92)
         mat(805) = rxt(93) + rxt(95) + rxt(96) + rxt(97) + rxt(105) + rxt(107) + rxt(108) &
                      + rxt(109)

         mat(563) = -( het_rates(59) )
         mat(710) = .700_r8*rxt(68)

         mat(469) = -( het_rates(83) )

         mat(425) = -( het_rates(64) )

         mat(583) = -( rxt(61) + het_rates(50) )
         mat(280) = rxt(62)
         mat(130) = rxt(69)
         mat(319) = .400_r8*rxt(83)
         mat(120) = rxt(84)

         mat(295) = -( het_rates(49) )

         mat(278) = -( rxt(62) + het_rates(65) )

         mat(765) = -( het_rates(48) )
         mat(199) = .600_r8*rxt(64) + rxt(308)
         mat(617) = 1.340_r8*rxt(66)
         mat(718) = .300_r8*rxt(68)
         mat(158) = rxt(72)
         mat(365) = rxt(73)
         mat(639) = rxt(74)
         mat(594) = rxt(78)
         mat(177) = rxt(80)
         mat(244) = .130_r8*rxt(81)
         mat(121) = rxt(84)

         mat(212) = -( rxt(63) + het_rates(54) )

         mat(198) = -( rxt(64) + rxt(308) + het_rates(58) )

         mat(140) = -( het_rates(82) )

         mat(93) = -( het_rates(45) )

         mat(258) = -( het_rates(44) )

         mat(25) = -( het_rates(71) )

         mat(223) = -( rxt(65) + rxt(354) + het_rates(81) )

         mat(28) = -( het_rates(70) )

         mat(111) = -( het_rates(73) )

         mat(352) = -( het_rates(84) )

         mat(314) = -( rxt(83) + het_rates(85) )

         mat(174) = -( rxt(80) + het_rates(72) )
         mat(313) = .800_r8*rxt(83)

         mat(325) = -( het_rates(74) )

         mat(118) = -( rxt(84) + het_rates(75) )

         mat(44) = -( het_rates(94) )

         mat(49) = -( het_rates(95) )

         mat(232) = -( het_rates(96) )

         mat(144) = -( rxt(85) + het_rates(97) )

         mat(69) = -( het_rates(98) )

         mat(523) = -( het_rates(100) )

         mat(191) = -( rxt(86) + het_rates(101) )

         mat(240) = -( rxt(81) + het_rates(86) )
         mat(146) = .900_r8*rxt(85)

         mat(344) = -( rxt(82) + het_rates(53) )
         mat(241) = .130_r8*rxt(81)
         mat(147) = .450_r8*rxt(85)

         mat(673) = -( het_rates(88) )

         mat(716) = -( rxt(68) + het_rates(77) )
         mat(274) = .402_r8*rxt(77)
         mat(195) = rxt(86)

         mat(613) = -( rxt(66) + rxt(67) + het_rates(78) )
         mat(271) = .288_r8*rxt(77)
         mat(194) = rxt(86)

         mat(697) = -( het_rates(79) )

         mat(123) = -( het_rates(80) )

         mat(736) = -( het_rates(76) )
         mat(225) = rxt(65) + rxt(354)
         mat(616) = .660_r8*rxt(66)

         mat(445) = -( het_rates(46) )
         mat(176) = rxt(80)

         mat(128) = -( rxt(69) + het_rates(47) )

         mat(286) = -( het_rates(99) )

         mat(34) = -( het_rates(60) )

         mat(500) = -( het_rates(61) )

         mat(150) = -( rxt(71) + het_rates(62) )

         mat(363) = -( rxt(73) + het_rates(63) )
         mat(151) = .820_r8*rxt(71)
         mat(317) = .250_r8*rxt(83)
         mat(192) = .100_r8*rxt(86)

         mat(156) = -( rxt(72) + het_rates(69) )

         mat(254) = -( het_rates(18) )

         mat(89) = -( het_rates(51) )

         mat(493) = -( rxt(79) + het_rates(52) )

         mat(592) = -( rxt(78) + het_rates(66) )

         mat(380) = -( het_rates(55) )

         mat(179) = -( rxt(289) + het_rates(56) )
         mat(53) = rxt(70)

         mat(52) = -( rxt(70) + het_rates(57) )

         mat(137) = -( het_rates(87) )

         mat(601) = -( het_rates(67) )

         mat(638) = -( rxt(74) + het_rates(68) )
         mat(243) = .180_r8*rxt(81)
         mat(148) = .450_r8*rxt(85)

         mat(548) = -( het_rates(89) )

         mat(513) = -( rxt(76) + het_rates(90) )

         mat(653) = -( het_rates(91) )

         mat(133) = -( rxt(75) + het_rates(92) )

         mat(270) = -( rxt(77) + het_rates(93) )

         mat(99) = -( het_rates(129) )

         mat(266) = -( het_rates(130) )

         mat(162) = -( rxt(277) + het_rates(131) )

         mat(55) = -( rxt(55) + het_rates(132) )
         mat(972) = rxt(136)*y(105) + rxt(137)*y(106) + 2.000_r8*rxt(138)*y(114) &
                      + 2.000_r8*rxt(139)*y(115) + rxt(140)*y(107) + rxt(142)*y(113) &
                      + rxt(145)*y(111) + rxt(146)*y(110) + rxt(147)*y(116) &
                      + 2.000_r8*rxt(148)*y(117)
         mat(1171) = rxt(254)*y(107) + rxt(258)*y(113)

         mat(73) = -( rxt(56) + het_rates(133) )
         mat(975) = rxt(135)*y(104) + rxt(137)*y(106) + rxt(141)*y(112)
         mat(1174) = rxt(257)*y(112)

         mat(80) = -( rxt(57) + het_rates(134) )
         mat(388) = rxt(249)*y(15)

         mat(389) = -( rxt(249)*y(15) + het_rates(135) )
         mat(56) = 2.000_r8*rxt(55)
         mat(74) = rxt(56)
         mat(81) = rxt(57)
         mat(976) = rxt(139)*y(115) + rxt(146)*y(110)

         mat(39) = -( het_rates(141) )

         mat(84) = -( het_rates(142) )

         mat(1) = -( rxt(409) + het_rates(143) )

         mat(31) = -( het_rates(144) )

         mat(2) = -( rxt(411) + het_rates(145) )

         mat(3) = -( rxt(412) + het_rates(146) )

         mat(4) = -( rxt(410) + het_rates(140) )

         mat(5) = -( rxt(397) + rxt(405) + het_rates(136) )

         mat(7) = -( rxt(406) + het_rates(137) )
         mat(6) = rxt(397)

         mat(8) = -( rxt(403) + rxt(407) + het_rates(138) )

         mat(10) = -( rxt(408) + het_rates(139) )
         mat(9) = rxt(403)

         mat(11) = -( rxt(413) + het_rates(147) )

         mat(12) = -( rxt(414) + het_rates(148) )

         mat(13) = -( rxt(415) + het_rates(149) )

         mat(14) = -( rxt(416) + het_rates(150) )

         mat(15) = -( rxt(417) + het_rates(151) )

         mat(16) = -( rxt(418) + het_rates(152) )

         mat(17) = -( rxt(419) + het_rates(153) )

         mat(18) = -( rxt(420) + het_rates(154) )


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
