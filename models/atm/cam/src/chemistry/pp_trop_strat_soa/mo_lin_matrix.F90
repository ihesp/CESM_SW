




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

         mat(833) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(1168) = rxt(84)

         mat(1178) = -( rxt(84) + het_rates(2) )
         mat(843) = rxt(3)
         mat(980) = rxt(5)
         mat(1041) = rxt(6)
         mat(130) = rxt(8)
         mat(885) = rxt(10)
         mat(720) = rxt(19)
         mat(804) = rxt(22)
         mat(81) = rxt(23)
         mat(1062) = rxt(30)
         mat(927) = rxt(87) + rxt(88)
         mat(97) = rxt(135)

         mat(920) = -( rxt(87) + rxt(88) + rxt(90)*y(4) + rxt(91)*y(4) + rxt(93)*y(107) &
                      + rxt(94)*y(108) + rxt(95)*y(109) + rxt(96)*y(117) + rxt(97)*y(118) &
                      + rxt(98)*y(110) + rxt(99)*y(115) + rxt(100)*y(116) &
                      + rxt(101)*y(111) + rxt(102)*y(106) + rxt(103)*y(114) &
                      + rxt(104)*y(113) + rxt(105)*y(119) + rxt(106)*y(120) &
                      + rxt(107)*y(121) + rxt(108)*y(122) + rxt(109)*y(12) &
                      + rxt(110)*y(12) + rxt(111)*y(12) + het_rates(3) )
         mat(836) = rxt(2)
         mat(718) = rxt(18)

         mat(528) = -( het_rates(18) )
         mat(1138) = rxt(16)
         mat(714) = rxt(18)
         mat(911) = rxt(111)*y(12)

         mat(536) = -( het_rates(17) )
         mat(1139) = rxt(15) + rxt(16)
         mat(560) = rxt(56)
         mat(570) = 1.340_r8*rxt(62)
         mat(611) = .700_r8*rxt(63)
         mat(583) = rxt(69)
         mat(477) = rxt(71)
         mat(417) = rxt(74)
         mat(269) = .450_r8*rxt(76)
         mat(354) = 2.000_r8*rxt(77)
         mat(989) = rxt(200)*y(105)

         mat(94) = -( rxt(135) + het_rates(5) )
         mat(929) = rxt(5)

         mat(974) = -( rxt(5) + het_rates(6) )
         mat(1035) = rxt(6) + .500_r8*rxt(350)
         mat(128) = rxt(8)
         mat(879) = rxt(11)
         mat(95) = rxt(135)
         mat(921) = 2.000_r8*rxt(90)*y(4)

         mat(1037) = -( rxt(6) + rxt(350) + het_rates(7) )
         mat(129) = rxt(7) + rxt(147)
         mat(437) = rxt(9)
         mat(881) = rxt(10)
         mat(198) = rxt(13) + rxt(156)
         mat(472) = rxt(28)
         mat(302) = rxt(34)
         mat(222) = .600_r8*rxt(59) + rxt(255)
         mat(293) = rxt(60) + rxt(301)
         mat(480) = rxt(71)

         mat(1274) = -( rxt(201)*y(105) + rxt(202)*y(112) + rxt(203)*y(110) &
                      + rxt(204)*y(106) + rxt(206)*y(115) + rxt(207)*y(116) &
                      + rxt(208)*y(122) + rxt(209)*y(121) + rxt(212)*y(12) + het_rates(20) &
       )
         mat(438) = rxt(9)
         mat(200) = rxt(12)
         mat(212) = rxt(14)
         mat(721) = rxt(17)
         mat(333) = 2.000_r8*rxt(20)
         mat(462) = rxt(25)
         mat(407) = rxt(31)
         mat(288) = rxt(57)
         mat(253) = rxt(58)
         mat(140) = rxt(64)
         mat(78) = rxt(65)
         mat(173) = rxt(66)
         mat(179) = rxt(67)
         mat(144) = rxt(70)
         mat(344) = rxt(78)
         mat(156) = rxt(79)
         mat(206) = rxt(80)
         mat(239) = rxt(81)
         mat(1042) = .500_r8*rxt(350)
         mat(928) = rxt(109)*y(12)

         mat(876) = -( rxt(10) + rxt(11) + rxt(349) + het_rates(8) )
         mat(127) = rxt(7) + rxt(8) + rxt(147)
         mat(197) = rxt(12)
         mat(469) = rxt(27)
         mat(301) = rxt(33)
         mat(221) = .400_r8*rxt(59)

         mat(434) = -( rxt(9) + het_rates(9) )
         mat(126) = 2.000_r8*rxt(348) + 2.000_r8*rxt(380) + 2.000_r8*rxt(386) &
                      + 2.000_r8*rxt(391)
         mat(853) = rxt(349)
         mat(1019) = .500_r8*rxt(350)
         mat(464) = rxt(381) + rxt(387) + rxt(392)
         mat(299) = rxt(382) + rxt(390) + rxt(393)

         mat(195) = -( rxt(12) + rxt(13) + rxt(156) + het_rates(10) )

         mat(125) = -( rxt(7) + rxt(8) + rxt(147) + rxt(348) + rxt(380) + rxt(386) &
                      + rxt(391) + het_rates(11) )

         mat(756) = -( het_rates(13) )
         mat(565) = rxt(56)
         mat(251) = rxt(58)
         mat(220) = .400_r8*rxt(59)
         mat(620) = .300_r8*rxt(63)
         mat(376) = rxt(68)
         mat(914) = rxt(109)*y(12)
         mat(994) = rxt(163)*y(12)
         mat(1260) = rxt(212)*y(12)

         mat(207) = -( rxt(14) + het_rates(14) )

         mat(98) = -( het_rates(39) )

         mat(55) = -( het_rates(40) )

         mat(1152) = -( rxt(15) + rxt(16) + het_rates(16) )
         mat(211) = rxt(14)
         mat(287) = rxt(57)
         mat(580) = 1.340_r8*rxt(61)
         mat(178) = rxt(67)
         mat(482) = rxt(71)
         mat(279) = .690_r8*rxt(72)
         mat(557) = rxt(73)
         mat(419) = rxt(74)
         mat(343) = .100_r8*rxt(78)
         mat(246) = rxt(226)
         mat(118) = 2.000_r8*rxt(236)
         mat(926) = rxt(110)*y(12) + rxt(111)*y(12)

         mat(724) = -( rxt(116) + het_rates(19) )
         mat(209) = rxt(14)
         mat(1141) = 2.000_r8*rxt(15)
         mat(716) = rxt(17) + 2.000_r8*rxt(19)
         mat(894) = rxt(26)
         mat(410) = rxt(32)
         mat(913) = rxt(110)*y(12)

         mat(1131) = -( rxt(358) + het_rates(21) )
         mat(199) = rxt(13) + rxt(156)
         mat(568) = rxt(56)
         mat(286) = rxt(57)
         mat(579) = 1.340_r8*rxt(61) + .660_r8*rxt(62)
         mat(139) = rxt(64)
         mat(172) = rxt(66)
         mat(590) = rxt(69)
         mat(481) = rxt(71)
         mat(278) = rxt(72)
         mat(556) = rxt(73)
         mat(418) = 2.000_r8*rxt(74)
         mat(272) = .560_r8*rxt(76)
         mat(355) = 2.000_r8*rxt(77)
         mat(342) = .900_r8*rxt(78)
         mat(238) = rxt(81)
         mat(729) = rxt(116)
         mat(245) = rxt(226)
         mat(117) = rxt(235) + rxt(236)
         mat(925) = rxt(110)*y(12)
         mat(1005) = rxt(200)*y(105) + rxt(205)*y(106)
         mat(1271) = rxt(201)*y(105) + rxt(204)*y(106)

         mat(327) = -( rxt(20) + het_rates(22) )
         mat(1089) = .500_r8*rxt(358)

         mat(715) = -( rxt(17) + rxt(18) + rxt(19) + het_rates(124) )
         mat(1258) = rxt(201)*y(105) + rxt(202)*y(112) + rxt(203)*y(110) + rxt(204)*y(106) &
                      + rxt(208)*y(122) + rxt(212)*y(12)

         mat(1002) = -( rxt(163)*y(12) + rxt(200)*y(105) + rxt(205)*y(106) &
                      + rxt(210)*y(122) + rxt(211)*y(121) + het_rates(25) )
         mat(89) = 2.000_r8*rxt(21)
         mat(799) = rxt(22)
         mat(60) = 2.000_r8*rxt(24)
         mat(460) = rxt(25)
         mat(902) = rxt(26)
         mat(471) = rxt(27)
         mat(110) = rxt(29)
         mat(922) = 3.000_r8*rxt(93)*y(107) + 2.000_r8*rxt(94)*y(108) &
                      + 3.000_r8*rxt(95)*y(109) + 2.000_r8*rxt(96)*y(117) + rxt(97)*y(118) &
                      + rxt(98)*y(110) + 2.000_r8*rxt(99)*y(115) + rxt(100)*y(116) &
                      + 4.000_r8*rxt(101)*y(111) + rxt(103)*y(114)
         mat(1268) = rxt(201)*y(105) + 3.000_r8*rxt(202)*y(112) + rxt(203)*y(110) &
                      + 2.000_r8*rxt(206)*y(115) + rxt(207)*y(116)

         mat(88) = -( rxt(21) + het_rates(26) )

         mat(793) = -( rxt(22) + het_rates(27) )
         mat(80) = rxt(23)
         mat(468) = rxt(28)
         mat(59) = 2.000_r8*rxt(175)

         mat(79) = -( rxt(23) + het_rates(28) )

         mat(58) = -( rxt(24) + rxt(175) + het_rates(29) )

         mat(899) = -( rxt(26) + het_rates(30) )
         mat(999) = rxt(163)*y(12) + 2.000_r8*rxt(200)*y(105) + rxt(205)*y(106) &
                      + rxt(210)*y(122) + rxt(211)*y(121)

         mat(456) = -( rxt(25) + het_rates(31) )
         mat(465) = rxt(381) + rxt(387) + rxt(392)

         mat(466) = -( rxt(27) + rxt(28) + rxt(381) + rxt(387) + rxt(392) + het_rates(32) &
       )

         mat(108) = -( rxt(29) + het_rates(33) )

         mat(773) = -( het_rates(34) )
         mat(109) = rxt(29)
         mat(1050) = rxt(30)
         mat(402) = rxt(31)
         mat(411) = rxt(32)
         mat(300) = rxt(33)
         mat(915) = rxt(102)*y(106) + rxt(103)*y(114) + rxt(104)*y(113) &
                      + 2.000_r8*rxt(105)*y(119) + 2.000_r8*rxt(106)*y(120) &
                      + 3.000_r8*rxt(107)*y(121) + 2.000_r8*rxt(108)*y(122)
         mat(1261) = rxt(204)*y(106) + 2.000_r8*rxt(208)*y(122) + 3.000_r8*rxt(209)*y(121)
         mat(995) = rxt(205)*y(106) + 2.000_r8*rxt(210)*y(122) + 3.000_r8*rxt(211)*y(121)

         mat(1059) = -( rxt(30) + het_rates(35) )
         mat(303) = rxt(34)

         mat(408) = -( rxt(32) + het_rates(36) )

         mat(400) = -( rxt(31) + het_rates(37) )
         mat(298) = rxt(382) + rxt(390) + rxt(393)


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

         mat(297) = -( rxt(33) + rxt(34) + rxt(382) + rxt(390) + rxt(393) + het_rates(38) &
       )

         mat(500) = -( het_rates(56) )
         mat(610) = .700_r8*rxt(63)

         mat(440) = -( het_rates(80) )

         mat(379) = -( het_rates(61) )

         mat(561) = -( rxt(56) + het_rates(47) )
         mat(284) = rxt(57)
         mat(138) = rxt(64)
         mat(340) = .400_r8*rxt(78)
         mat(154) = rxt(79)

         mat(314) = -( het_rates(46) )

         mat(281) = -( rxt(57) + het_rates(62) )

         mat(703) = -( het_rates(45) )
         mat(219) = .600_r8*rxt(59) + rxt(255)
         mat(575) = 1.340_r8*rxt(61)
         mat(617) = .300_r8*rxt(63)
         mat(176) = rxt(67)
         mat(374) = rxt(68)
         mat(585) = rxt(69)
         mat(555) = rxt(73)
         mat(216) = rxt(75)
         mat(271) = .130_r8*rxt(76)
         mat(155) = rxt(79)

         mat(248) = -( rxt(58) + het_rates(51) )

         mat(218) = -( rxt(59) + rxt(255) + het_rates(55) )

         mat(191) = -( het_rates(79) )

         mat(119) = -( het_rates(42) )

         mat(157) = -( het_rates(41) )

         mat(40) = -( het_rates(68) )

         mat(289) = -( rxt(60) + rxt(301) + het_rates(78) )

         mat(43) = -( het_rates(67) )

         mat(145) = -( het_rates(70) )

         mat(361) = -( het_rates(81) )

         mat(335) = -( rxt(78) + het_rates(82) )

         mat(213) = -( rxt(75) + het_rates(69) )
         mat(334) = .800_r8*rxt(78)

         mat(346) = -( het_rates(71) )

         mat(152) = -( rxt(79) + het_rates(72) )

         mat(68) = -( het_rates(91) )

         mat(73) = -( het_rates(92) )

         mat(259) = -( het_rates(93) )

         mat(201) = -( rxt(80) + het_rates(94) )

         mat(90) = -( het_rates(95) )

         mat(486) = -( het_rates(103) )

         mat(233) = -( rxt(81) + het_rates(104) )

         mat(267) = -( rxt(76) + het_rates(83) )
         mat(203) = .900_r8*rxt(80)

         mat(353) = -( rxt(77) + het_rates(50) )
         mat(268) = .130_r8*rxt(76)
         mat(204) = .450_r8*rxt(80)

         mat(46) = -( het_rates(96) )

         mat(181) = -( het_rates(97) )

         mat(1) = -( het_rates(98) )

         mat(49) = -( het_rates(99) )

         mat(226) = -( het_rates(100) )

         mat(2) = -( het_rates(101) )

         mat(658) = -( het_rates(85) )

         mat(615) = -( rxt(63) + het_rates(74) )
         mat(276) = .402_r8*rxt(72)
         mat(237) = rxt(81)

         mat(571) = -( rxt(61) + rxt(62) + het_rates(75) )
         mat(274) = .288_r8*rxt(72)
         mat(236) = rxt(81)

         mat(636) = -( het_rates(76) )

         mat(131) = -( het_rates(77) )

         mat(676) = -( het_rates(73) )
         mat(291) = rxt(60) + rxt(301)
         mat(574) = .660_r8*rxt(61)

         mat(391) = -( het_rates(43) )
         mat(215) = rxt(75)

         mat(136) = -( rxt(64) + het_rates(44) )

         mat(318) = -( het_rates(102) )

         mat(61) = -( het_rates(57) )

         mat(423) = -( het_rates(58) )

         mat(168) = -( rxt(66) + het_rates(59) )

         mat(372) = -( rxt(68) + het_rates(60) )
         mat(169) = .820_r8*rxt(66)
         mat(338) = .250_r8*rxt(78)
         mat(234) = .100_r8*rxt(81)

         mat(174) = -( rxt(67) + het_rates(66) )

         mat(254) = -( het_rates(15) )

         mat(111) = -( het_rates(48) )

         mat(416) = -( rxt(74) + het_rates(49) )
         mat(116) = rxt(235)

         mat(553) = -( rxt(73) + het_rates(63) )

         mat(307) = -( het_rates(52) )

         mat(115) = -( rxt(235) + rxt(236) + het_rates(53) )
         mat(77) = rxt(65)

         mat(76) = -( rxt(65) + het_rates(54) )

         mat(165) = -( het_rates(84) )

         mat(542) = -( het_rates(64) )

         mat(584) = -( rxt(69) + het_rates(65) )
         mat(270) = .180_r8*rxt(76)
         mat(205) = .450_r8*rxt(80)

         mat(516) = -( het_rates(86) )

         mat(476) = -( rxt(71) + het_rates(87) )

         mat(599) = -( het_rates(88) )

         mat(141) = -( rxt(70) + het_rates(89) )

         mat(273) = -( rxt(72) + het_rates(90) )

         mat(82) = -( het_rates(125) )

         mat(187) = -( het_rates(126) )

         mat(241) = -( rxt(226) + het_rates(127) )

         mat(66) = -( het_rates(142) )

         mat(103) = -( het_rates(143) )

         mat(3) = -( rxt(363) + het_rates(144) )

         mat(52) = -( het_rates(145) )

         mat(4) = -( rxt(369) + het_rates(146) )

         mat(5) = -( rxt(370) + het_rates(147) )

         mat(6) = -( rxt(364) + het_rates(132) )

         mat(7) = -( rxt(365) + het_rates(133) )

         mat(8) = -( rxt(367) + het_rates(134) )

         mat(9) = -( rxt(366) + het_rates(135) )

         mat(10) = -( rxt(368) + het_rates(136) )

         mat(11) = -( het_rates(137) )

         mat(12) = -( het_rates(138) )

         mat(13) = -( het_rates(139) )

         mat(14) = -( het_rates(140) )

         mat(15) = -( het_rates(141) )

         mat(16) = -( rxt(351) + rxt(359) + het_rates(128) )

         mat(18) = -( rxt(360) + het_rates(129) )
         mat(17) = rxt(351)

         mat(19) = -( rxt(357) + rxt(361) + het_rates(130) )

         mat(21) = -( rxt(362) + het_rates(131) )
         mat(20) = rxt(357)

         mat(22) = -( rxt(371) + het_rates(148) )

         mat(23) = -( rxt(372) + het_rates(149) )

         mat(24) = -( rxt(373) + het_rates(150) )

         mat(25) = -( rxt(374) + het_rates(151) )

         mat(26) = -( rxt(375) + het_rates(152) )

         mat(27) = -( rxt(376) + het_rates(153) )

         mat(28) = -( rxt(377) + het_rates(154) )

         mat(29) = -( rxt(378) + het_rates(155) )

         mat(30) = -( rxt(397) + het_rates(156) )

         mat(31) = -( rxt(398) + het_rates(157) )

         mat(32) = -( rxt(399) + het_rates(158) )

         mat(33) = -( het_rates(159) )

         mat(34) = -( rxt(400) + het_rates(160) )

         mat(35) = -( rxt(401) + het_rates(161) )

         mat(36) = -( rxt(402) + het_rates(162) )

         mat(37) = -( rxt(379) + het_rates(165) )

         mat(38) = -( rxt(82) + het_rates(166) )

         mat(39) = -( rxt(83) + het_rates(167) )


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
