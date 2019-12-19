




      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: nlnmat

      contains

      subroutine nlnmat01( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)


!----------------------------------------------
! ... local variables
!----------------------------------------------

!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------


         mat(518) = -(rxt(81)*y(2) + rxt(99)*y(3) + rxt(121)*y(9) + rxt(124)*y(10) &
                      + rxt(146)*y(21) + rxt(151)*y(22) + rxt(158)*y(23) + rxt(161) &
                      *y(27) + rxt(187)*y(36) + rxt(216)*y(60) + rxt(219)*y(61))
         mat(400) = -rxt(81)*y(1)
         mat(337) = -rxt(99)*y(1)
         mat(423) = -rxt(121)*y(1)
         mat(447) = -rxt(124)*y(1)
         mat(226) = -rxt(146)*y(1)
         mat(478) = -rxt(151)*y(1)
         mat(362) = -rxt(158)*y(1)
         mat(540) = -rxt(161)*y(1)
         mat(581) = -rxt(187)*y(1)
         mat(124) = -rxt(216)*y(1)
         mat(601) = -rxt(219)*y(1)

         mat(400) = mat(400) + rxt(80)*y(4)
         mat(276) = rxt(80)*y(2)

         mat(395) = -(rxt(80)*y(4) + rxt(81)*y(1) + 4._r8*rxt(82)*y(2) + rxt(119)*y(9) &
                      + (rxt(122) + rxt(123)) * y(10) + rxt(130)*y(11) + rxt(142) &
                      *y(18) + rxt(150)*y(22) + rxt(157)*y(23) + rxt(160)*y(24) &
                      + rxt(168)*y(29) + rxt(180)*y(32) + rxt(181)*y(33) + rxt(184) &
                      *y(34) + rxt(190)*y(37) + rxt(200)*y(38) + rxt(201)*y(39) &
                      + rxt(202)*y(40) + rxt(212)*y(59) + (rxt(248) + rxt(257) &
                      ) * y(52) + rxt(254)*y(54))
         mat(272) = -rxt(80)*y(2)
         mat(513) = -rxt(81)*y(2)
         mat(418) = -rxt(119)*y(2)
         mat(442) = -(rxt(122) + rxt(123)) * y(2)
         mat(492) = -rxt(130)*y(2)
         mat(235) = -rxt(142)*y(2)
         mat(473) = -rxt(150)*y(2)
         mat(357) = -rxt(157)*y(2)
         mat(115) = -rxt(160)*y(2)
         mat(313) = -rxt(168)*y(2)
         mat(558) = -rxt(180)*y(2)
         mat(197) = -rxt(181)*y(2)
         mat(207) = -rxt(184)*y(2)
         mat(290) = -rxt(190)*y(2)
         mat(129) = -rxt(200)*y(2)
         mat(141) = -rxt(201)*y(2)
         mat(97) = -rxt(202)*y(2)
         mat(54) = -rxt(212)*y(2)
         mat(110) = -(rxt(248) + rxt(257)) * y(2)
         mat(74) = -rxt(254)*y(2)

         mat(332) = (rxt(94)+rxt(95))*y(4)
         mat(272) = mat(272) + (rxt(94)+rxt(95))*y(3) + rxt(116)*y(8) + rxt(253)*y(54) &
                      + rxt(246)*y(55) + rxt(215)*y(60) + rxt(218)*y(61)
         mat(180) = rxt(116)*y(4) + rxt(117)*y(9) + rxt(118)*y(10) + rxt(250)*y(53)
         mat(418) = mat(418) + rxt(117)*y(8)
         mat(442) = mat(442) + rxt(118)*y(8)
         mat(473) = mat(473) + 2.000_r8*rxt(153)*y(22)
         mat(224) = rxt(149)*y(23)
         mat(357) = mat(357) + rxt(149)*y(21)
         mat(152) = rxt(250)*y(8) + 1.150_r8*rxt(259)*y(57)
         mat(74) = mat(74) + rxt(253)*y(4)
         mat(87) = rxt(246)*y(4)
         mat(170) = rxt(258)*y(57)
         mat(162) = 1.150_r8*rxt(259)*y(53) + rxt(258)*y(56)
         mat(122) = rxt(215)*y(4)
         mat(596) = rxt(218)*y(4)

         mat(330) = -((rxt(94) + rxt(95)) * y(4) + rxt(96)*y(82) + rxt(99)*y(1) &
                      + rxt(112)*y(32) + rxt(113)*y(38))
         mat(270) = -(rxt(94) + rxt(95)) * y(3)
         mat(249) = -rxt(96)*y(3)
         mat(511) = -rxt(99)*y(3)
         mat(556) = -rxt(112)*y(3)
         mat(128) = -rxt(113)*y(3)

         mat(270) = mat(270) + rxt(114)*y(58)
         mat(151) = .850_r8*rxt(259)*y(57)
         mat(91) = rxt(114)*y(4)
         mat(161) = .850_r8*rxt(259)*y(53)

         mat(269) = -(rxt(80)*y(2) + rxt(90)*y(6) + rxt(94)*y(3) + rxt(114)*y(58) &
                      + rxt(116)*y(8) + rxt(145)*y(21) + rxt(215)*y(60) + rxt(218) &
                      *y(61) + rxt(246)*y(55) + (rxt(252) + rxt(253)) * y(54) + rxt(255) &
                      *y(52))
         mat(390) = -rxt(80)*y(4)
         mat(24) = -rxt(90)*y(4)
         mat(329) = -rxt(94)*y(4)
         mat(90) = -rxt(114)*y(4)
         mat(178) = -rxt(116)*y(4)
         mat(222) = -rxt(145)*y(4)
         mat(121) = -rxt(215)*y(4)
         mat(591) = -rxt(218)*y(4)
         mat(86) = -rxt(246)*y(4)
         mat(73) = -(rxt(252) + rxt(253)) * y(4)
         mat(108) = -rxt(255)*y(4)

         mat(508) = 2.000_r8*rxt(81)*y(2) + 2.000_r8*rxt(99)*y(3) + rxt(121)*y(9) &
                      + rxt(124)*y(10) + rxt(151)*y(22) + rxt(146)*y(21) &
                      + 2.000_r8*rxt(158)*y(23) + rxt(161)*y(27) + rxt(187)*y(36) &
                      + rxt(216)*y(60) + rxt(219)*y(61)
         mat(390) = mat(390) + 2.000_r8*rxt(81)*y(1) + 2.000_r8*rxt(82)*y(2) + rxt(89) &
                      *y(6) + rxt(122)*y(10) + rxt(150)*y(22) + rxt(130)*y(11) &
                      + rxt(157)*y(23) + rxt(168)*y(29) + rxt(190)*y(37)
         mat(329) = mat(329) + 2.000_r8*rxt(99)*y(1)
         mat(269) = mat(269) + 2.000_r8*rxt(90)*y(6)
         mat(24) = mat(24) + rxt(89)*y(2) + 2.000_r8*rxt(90)*y(4)
         mat(413) = rxt(121)*y(1) + rxt(251)*y(53)
         mat(437) = rxt(124)*y(1) + rxt(122)*y(2)
         mat(468) = rxt(151)*y(1) + rxt(150)*y(2) + rxt(134)*y(13) + rxt(152)*y(23) &
                      + rxt(170)*y(29)
         mat(489) = rxt(130)*y(2) + rxt(132)*y(23)
         mat(77) = rxt(134)*y(22)
         mat(188) = rxt(138)*y(23)
         mat(222) = mat(222) + rxt(146)*y(1) + rxt(148)*y(23)
         mat(352) = 2.000_r8*rxt(158)*y(1) + rxt(157)*y(2) + rxt(152)*y(22) + rxt(132) &
                      *y(11) + rxt(138)*y(16) + rxt(148)*y(21) + 2.000_r8*rxt(159) &
                      *y(23) + rxt(164)*y(27) + rxt(171)*y(29) + rxt(188)*y(36) &
                      + rxt(192)*y(37)
         mat(531) = rxt(161)*y(1) + rxt(164)*y(23)
         mat(308) = rxt(168)*y(2) + rxt(170)*y(22) + rxt(171)*y(23) + ( &
                      + 2.000_r8*rxt(174)+2.000_r8*rxt(175))*y(29) + (rxt(196) &
                       +rxt(197))*y(37)
         mat(571) = rxt(187)*y(1) + rxt(188)*y(23)
         mat(285) = rxt(190)*y(2) + rxt(192)*y(23) + (rxt(196)+rxt(197))*y(29) &
                      + 2.000_r8*rxt(198)*y(37)
         mat(150) = rxt(251)*y(9)
         mat(121) = mat(121) + rxt(216)*y(1)
         mat(591) = mat(591) + rxt(219)*y(1)

         mat(26) = -(rxt(83)*y(2) + rxt(84)*y(4) + rxt(86)*y(1))
         mat(368) = -rxt(83)*y(5)
         mat(256) = -rxt(84)*y(5)
         mat(503) = -rxt(86)*y(5)

         mat(323) = rxt(94)*y(4)
         mat(256) = mat(256) + rxt(94)*y(3)

         mat(23) = -(rxt(89)*y(2) + rxt(90)*y(4))
         mat(367) = -rxt(89)*y(6)
         mat(255) = -rxt(90)*y(6)

         mat(502) = rxt(86)*y(5)
         mat(367) = mat(367) + rxt(83)*y(5)
         mat(255) = mat(255) + rxt(84)*y(5)
         mat(25) = rxt(86)*y(1) + rxt(83)*y(2) + rxt(84)*y(4)

         mat(177) = -(rxt(116)*y(4) + rxt(117)*y(9) + rxt(118)*y(10) + rxt(250)*y(53))
         mat(265) = -rxt(116)*y(8)
         mat(408) = -rxt(117)*y(8)
         mat(433) = -rxt(118)*y(8)
         mat(149) = -rxt(250)*y(8)

         mat(383) = rxt(254)*y(54) + rxt(115)*y(58)
         mat(265) = mat(265) + rxt(252)*y(54)
         mat(107) = 1.100_r8*rxt(260)*y(57)
         mat(72) = rxt(254)*y(2) + rxt(252)*y(4)
         mat(167) = .200_r8*rxt(258)*y(57)
         mat(89) = rxt(115)*y(2)
         mat(159) = 1.100_r8*rxt(260)*y(52) + .200_r8*rxt(258)*y(56)

         mat(419) = -(rxt(117)*y(8) + rxt(119)*y(2) + rxt(120)*y(23) + rxt(121)*y(1) &
                      + rxt(129)*y(11) + rxt(137)*y(16) + rxt(172)*y(29) + rxt(193) &
                      *y(37) + rxt(251)*y(53))
         mat(181) = -rxt(117)*y(9)
         mat(396) = -rxt(119)*y(9)
         mat(358) = -rxt(120)*y(9)
         mat(514) = -rxt(121)*y(9)
         mat(493) = -rxt(129)*y(9)
         mat(190) = -rxt(137)*y(9)
         mat(314) = -rxt(172)*y(9)
         mat(291) = -rxt(193)*y(9)
         mat(153) = -rxt(251)*y(9)

         mat(396) = mat(396) + rxt(122)*y(10)
         mat(273) = rxt(116)*y(8) + rxt(114)*y(58)
         mat(181) = mat(181) + rxt(116)*y(4)
         mat(443) = rxt(122)*y(2) + rxt(220)*y(61)
         mat(92) = rxt(114)*y(4)
         mat(597) = rxt(220)*y(10)

         mat(444) = -(rxt(118)*y(8) + (rxt(122) + rxt(123)) * y(2) + rxt(124)*y(1) &
                      + rxt(125)*y(11) + rxt(127)*y(22) + rxt(133)*y(23) + rxt(173) &
                      *y(29) + rxt(194)*y(37) + rxt(220)*y(61))
         mat(182) = -rxt(118)*y(10)
         mat(397) = -(rxt(122) + rxt(123)) * y(10)
         mat(515) = -rxt(124)*y(10)
         mat(494) = -rxt(125)*y(10)
         mat(475) = -rxt(127)*y(10)
         mat(359) = -rxt(133)*y(10)
         mat(315) = -rxt(173)*y(10)
         mat(292) = -rxt(194)*y(10)
         mat(598) = -rxt(220)*y(10)

         mat(515) = mat(515) + rxt(121)*y(9)
         mat(397) = mat(397) + rxt(119)*y(9) + rxt(130)*y(11)
         mat(420) = rxt(121)*y(1) + rxt(119)*y(2) + 2.000_r8*rxt(129)*y(11) + rxt(137) &
                      *y(16) + rxt(120)*y(23) + rxt(172)*y(29) + rxt(193)*y(37)
         mat(475) = mat(475) + rxt(131)*y(11) + rxt(134)*y(13)
         mat(494) = mat(494) + rxt(130)*y(2) + 2.000_r8*rxt(129)*y(9) + rxt(131)*y(22) &
                      + rxt(132)*y(23)
         mat(79) = rxt(134)*y(22)
         mat(191) = rxt(137)*y(9)
         mat(359) = mat(359) + rxt(120)*y(9) + rxt(132)*y(11)
         mat(315) = mat(315) + rxt(172)*y(9)
         mat(292) = mat(292) + rxt(193)*y(9)

         mat(476) = -(rxt(127)*y(10) + rxt(128)*y(12) + rxt(131)*y(11) + rxt(134) &
                      *y(13) + rxt(139)*y(17) + rxt(141)*y(18) + rxt(150)*y(2) + rxt(151) &
                      *y(1) + rxt(152)*y(23) + (4._r8*rxt(153) + 4._r8*rxt(154) &
                      ) * y(22) + rxt(156)*y(24) + (rxt(169) + rxt(170)) * y(29) &
                      + rxt(179)*y(32) + rxt(183)*y(33) + rxt(185)*y(34) + rxt(191) &
                      *y(37) + rxt(199)*y(38) + rxt(213)*y(59) + rxt(214)*y(60) &
                      + rxt(217)*y(61) + rxt(224)*y(62) + (rxt(226) + rxt(227) &
                      ) * y(65))
         mat(445) = -rxt(127)*y(22)
         mat(135) = -rxt(128)*y(22)
         mat(495) = -rxt(131)*y(22)
         mat(80) = -rxt(134)*y(22)
         mat(67) = -rxt(139)*y(22)
         mat(237) = -rxt(141)*y(22)
         mat(398) = -rxt(150)*y(22)
         mat(516) = -rxt(151)*y(22)
         mat(360) = -rxt(152)*y(22)
         mat(116) = -rxt(156)*y(22)
         mat(316) = -(rxt(169) + rxt(170)) * y(22)
         mat(561) = -rxt(179)*y(22)
         mat(198) = -rxt(183)*y(22)
         mat(209) = -rxt(185)*y(22)
         mat(293) = -rxt(191)*y(22)
         mat(130) = -rxt(199)*y(22)
         mat(55) = -rxt(213)*y(22)
         mat(123) = -rxt(214)*y(22)
         mat(599) = -rxt(217)*y(22)
         mat(218) = -rxt(224)*y(22)
         mat(42) = -(rxt(226) + rxt(227)) * y(22)

         mat(516) = mat(516) + rxt(146)*y(21) + rxt(158)*y(23)
         mat(398) = mat(398) + rxt(142)*y(18) + rxt(157)*y(23) + rxt(160)*y(24) &
                      + rxt(180)*y(32) + rxt(181)*y(33) + rxt(200)*y(38) + rxt(201) &
                      *y(39)
         mat(335) = 2.000_r8*rxt(96)*y(82) + rxt(112)*y(32) + rxt(113)*y(38)
         mat(421) = rxt(120)*y(23)
         mat(495) = mat(495) + rxt(132)*y(23)
         mat(237) = mat(237) + rxt(142)*y(2)
         mat(225) = rxt(146)*y(1) + 2.000_r8*rxt(147)*y(23)
         mat(360) = mat(360) + rxt(158)*y(1) + rxt(157)*y(2) + rxt(120)*y(9) &
                      + rxt(132)*y(11) + 2.000_r8*rxt(147)*y(21) + rxt(165)*y(27)
         mat(116) = mat(116) + rxt(160)*y(2)
         mat(252) = 2.000_r8*rxt(96)*y(3)
         mat(538) = rxt(165)*y(23)
         mat(561) = mat(561) + rxt(180)*y(2) + rxt(112)*y(3)
         mat(198) = mat(198) + rxt(181)*y(2)
         mat(130) = mat(130) + rxt(200)*y(2) + rxt(113)*y(3)
         mat(142) = rxt(201)*y(2)


      end subroutine nlnmat01

      subroutine nlnmat02( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)


!----------------------------------------------
! ... local variables
!----------------------------------------------

!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------


         mat(496) = -(rxt(125)*y(10) + rxt(129)*y(9) + rxt(130)*y(2) + rxt(131)*y(22) &
                      + rxt(132)*y(23) + rxt(140)*y(18) + rxt(228)*y(65))
         mat(446) = -rxt(125)*y(11)
         mat(422) = -rxt(129)*y(11)
         mat(399) = -rxt(130)*y(11)
         mat(477) = -rxt(131)*y(11)
         mat(361) = -rxt(132)*y(11)
         mat(238) = -rxt(140)*y(11)
         mat(43) = -rxt(228)*y(11)

         mat(517) = rxt(124)*y(10)
         mat(399) = mat(399) + rxt(123)*y(10) + rxt(184)*y(34) + rxt(202)*y(40)
         mat(446) = mat(446) + rxt(124)*y(1) + rxt(123)*y(2)
         mat(477) = mat(477) + rxt(128)*y(12) + rxt(185)*y(34)
         mat(136) = rxt(128)*y(22)
         mat(539) = rxt(186)*y(34)
         mat(210) = rxt(184)*y(2) + rxt(185)*y(22) + rxt(186)*y(27)
         mat(99) = rxt(202)*y(2)

         mat(132) = -(rxt(128)*y(22))
         mat(460) = -rxt(128)*y(12)

         mat(431) = rxt(127)*y(22)
         mat(460) = mat(460) + rxt(127)*y(10)
         mat(485) = rxt(140)*y(18) + rxt(228)*y(65)
         mat(228) = rxt(140)*y(11)
         mat(547) = (rxt(232)+rxt(237)+rxt(243))*y(34)
         mat(202) = (rxt(232)+rxt(237)+rxt(243))*y(32)
         mat(39) = rxt(228)*y(11)

         mat(75) = -(rxt(134)*y(22))
         mat(456) = -rxt(134)*y(13)

         mat(429) = rxt(133)*y(23)
         mat(343) = rxt(133)*y(10)


         mat(428) = rxt(125)*y(11)
         mat(484) = rxt(125)*y(10)

         mat(184) = -(rxt(137)*y(9) + rxt(138)*y(23))
         mat(409) = -rxt(137)*y(16)
         mat(347) = -rxt(138)*y(16)

         mat(461) = rxt(139)*y(17)
         mat(63) = rxt(139)*y(22)

         mat(62) = -(rxt(139)*y(22))
         mat(455) = -rxt(139)*y(17)

         mat(183) = rxt(138)*y(23)
         mat(342) = rxt(138)*y(16)

         mat(230) = -(rxt(140)*y(11) + rxt(141)*y(22) + rxt(142)*y(2) + rxt(166)*y(27) &
                      + rxt(189)*y(36))
         mat(487) = -rxt(140)*y(18)
         mat(466) = -rxt(141)*y(18)
         mat(388) = -rxt(142)*y(18)
         mat(529) = -rxt(166)*y(18)
         mat(569) = -rxt(189)*y(18)

         mat(411) = rxt(137)*y(16)
         mat(186) = rxt(137)*y(9)

         mat(220) = -(rxt(145)*y(4) + rxt(146)*y(1) + (rxt(147) + rxt(148) + rxt(149) &
                      ) * y(23))
         mat(267) = -rxt(145)*y(21)
         mat(506) = -rxt(146)*y(21)
         mat(349) = -(rxt(147) + rxt(148) + rxt(149)) * y(21)

         mat(387) = rxt(150)*y(22)
         mat(465) = rxt(150)*y(2) + rxt(141)*y(18) + rxt(213)*y(59) + rxt(214)*y(60) &
                      + rxt(217)*y(61)
         mat(229) = rxt(141)*y(22)
         mat(53) = rxt(213)*y(22)
         mat(120) = rxt(214)*y(22)
         mat(589) = rxt(217)*y(22)

         mat(356) = -(rxt(120)*y(9) + rxt(132)*y(11) + rxt(133)*y(10) + rxt(138)*y(16) &
                      + (rxt(147) + rxt(148) + rxt(149)) * y(21) + rxt(152)*y(22) &
                      + rxt(157)*y(2) + rxt(158)*y(1) + 4._r8*rxt(159)*y(23) + (rxt(164) &
                      + rxt(165)) * y(27) + rxt(171)*y(29) + rxt(188)*y(36) + rxt(192) &
                      *y(37))
         mat(417) = -rxt(120)*y(23)
         mat(491) = -rxt(132)*y(23)
         mat(441) = -rxt(133)*y(23)
         mat(189) = -rxt(138)*y(23)
         mat(223) = -(rxt(147) + rxt(148) + rxt(149)) * y(23)
         mat(472) = -rxt(152)*y(23)
         mat(394) = -rxt(157)*y(23)
         mat(512) = -rxt(158)*y(23)
         mat(534) = -(rxt(164) + rxt(165)) * y(23)
         mat(312) = -rxt(171)*y(23)
         mat(575) = -rxt(188)*y(23)
         mat(289) = -rxt(192)*y(23)

         mat(512) = mat(512) + rxt(151)*y(22)
         mat(394) = mat(394) + rxt(142)*y(18) + rxt(160)*y(24)
         mat(271) = rxt(145)*y(21)
         mat(417) = mat(417) + rxt(137)*y(16)
         mat(472) = mat(472) + rxt(151)*y(1) + rxt(131)*y(11) + rxt(156)*y(24) &
                      + rxt(169)*y(29) + rxt(191)*y(37) + rxt(224)*y(62) &
                      + .500_r8*rxt(226)*y(65)
         mat(491) = mat(491) + rxt(131)*y(22) + rxt(140)*y(18)
         mat(189) = mat(189) + rxt(137)*y(9)
         mat(234) = rxt(142)*y(2) + rxt(140)*y(11) + rxt(166)*y(27) + rxt(189)*y(36)
         mat(223) = mat(223) + rxt(145)*y(4)
         mat(114) = rxt(160)*y(2) + rxt(156)*y(22) + rxt(163)*y(27)
         mat(534) = mat(534) + rxt(166)*y(18) + rxt(163)*y(24)
         mat(312) = mat(312) + rxt(169)*y(22)
         mat(575) = mat(575) + rxt(189)*y(18)
         mat(289) = mat(289) + rxt(191)*y(22)
         mat(216) = rxt(224)*y(22)
         mat(41) = .500_r8*rxt(226)*y(22)

         mat(112) = -(rxt(156)*y(22) + rxt(160)*y(2) + rxt(163)*y(27))
         mat(457) = -rxt(156)*y(24)
         mat(375) = -rxt(160)*y(24)
         mat(524) = -rxt(163)*y(24)

         mat(457) = mat(457) + 2.000_r8*rxt(154)*y(22)
         mat(344) = 2.000_r8*rxt(159)*y(23)

         mat(247) = -(rxt(96)*y(3) + rxt(225)*y(63))
         mat(328) = -rxt(96)*y(82)
         mat(36) = -rxt(225)*y(82)

         mat(467) = 2.000_r8*rxt(153)*y(22) + rxt(128)*y(12) + rxt(134)*y(13) &
                      + rxt(139)*y(17) + rxt(141)*y(18) + rxt(152)*y(23) + rxt(156) &
                      *y(24) + rxt(179)*y(32) + rxt(183)*y(33) + rxt(199)*y(38)
         mat(133) = rxt(128)*y(22)
         mat(76) = rxt(134)*y(22)
         mat(66) = rxt(139)*y(22)
         mat(231) = rxt(141)*y(22)
         mat(221) = rxt(149)*y(23)
         mat(351) = rxt(152)*y(22) + rxt(149)*y(21)
         mat(113) = rxt(156)*y(22)
         mat(552) = rxt(179)*y(22) + (rxt(233)+rxt(238)+rxt(244))*y(33) + (rxt(234) &
                       +rxt(245))*y(39)
         mat(195) = rxt(183)*y(22) + (rxt(233)+rxt(238)+rxt(244))*y(32)
         mat(127) = rxt(199)*y(22)
         mat(139) = (rxt(234)+rxt(245))*y(32)

         mat(541) = -(rxt(161)*y(1) + rxt(163)*y(24) + (rxt(164) + rxt(165)) * y(23) &
                      + rxt(166)*y(18) + rxt(182)*y(33) + rxt(186)*y(34))
         mat(519) = -rxt(161)*y(27)
         mat(117) = -rxt(163)*y(27)
         mat(363) = -(rxt(164) + rxt(165)) * y(27)
         mat(240) = -rxt(166)*y(27)
         mat(199) = -rxt(182)*y(27)
         mat(211) = -rxt(186)*y(27)

         mat(401) = rxt(168)*y(29) + rxt(180)*y(32)
         mat(338) = rxt(112)*y(32)
         mat(424) = rxt(172)*y(29)
         mat(479) = rxt(169)*y(29) + rxt(179)*y(32)
         mat(319) = rxt(168)*y(2) + rxt(172)*y(9) + rxt(169)*y(22) + ( &
                      + 4.000_r8*rxt(174)+2.000_r8*rxt(176))*y(29) + rxt(196)*y(37) &
                      + rxt(221)*y(61)
         mat(564) = rxt(180)*y(2) + rxt(112)*y(3) + rxt(179)*y(22)
         mat(296) = rxt(196)*y(29)
         mat(602) = rxt(221)*y(29)


         mat(523) = rxt(186)*y(34)
         mat(301) = 2.000_r8*rxt(175)*y(29)
         mat(545) = (rxt(233)+rxt(238)+rxt(244))*y(33) + (rxt(232)+rxt(237)+rxt(243)) &
                      *y(34)
         mat(193) = (rxt(233)+rxt(238)+rxt(244))*y(32)
         mat(201) = rxt(186)*y(27) + (rxt(232)+rxt(237)+rxt(243))*y(32)

         mat(310) = -(rxt(168)*y(2) + (rxt(169) + rxt(170)) * y(22) + rxt(171)*y(23) &
                      + rxt(172)*y(9) + rxt(173)*y(10) + (4._r8*rxt(174) + 4._r8*rxt(175) &
                      + 4._r8*rxt(176) + 4._r8*rxt(177)) * y(29) + (rxt(195) + rxt(196) &
                      + rxt(197)) * y(37) + rxt(221)*y(61))
         mat(392) = -rxt(168)*y(29)
         mat(470) = -(rxt(169) + rxt(170)) * y(29)
         mat(354) = -rxt(171)*y(29)
         mat(415) = -rxt(172)*y(29)
         mat(439) = -rxt(173)*y(29)
         mat(287) = -(rxt(195) + rxt(196) + rxt(197)) * y(29)
         mat(593) = -rxt(221)*y(29)

         mat(510) = rxt(161)*y(27)
         mat(392) = mat(392) + rxt(181)*y(33) + rxt(184)*y(34)
         mat(470) = mat(470) + rxt(183)*y(33)
         mat(354) = mat(354) + rxt(165)*y(27)
         mat(532) = rxt(161)*y(1) + rxt(165)*y(23) + rxt(182)*y(33)
         mat(59) = rxt(223)*y(61)
         mat(196) = rxt(181)*y(2) + rxt(183)*y(22) + rxt(182)*y(27)
         mat(206) = rxt(184)*y(2)
         mat(593) = mat(593) + rxt(223)*y(30)

         mat(57) = -(rxt(223)*y(61))
         mat(586) = -rxt(223)*y(30)

         mat(303) = 2.000_r8*rxt(176)*y(29) + rxt(195)*y(37)
         mat(279) = rxt(195)*y(29)


         mat(300) = 2.000_r8*rxt(177)*y(29)

         mat(565) = -(rxt(112)*y(3) + rxt(179)*y(22) + rxt(180)*y(2) + (rxt(232) &
                      + rxt(237) + rxt(243)) * y(34) + (rxt(233) + rxt(238) + rxt(244) &
                      ) * y(33) + (rxt(234) + rxt(245)) * y(39))
         mat(339) = -rxt(112)*y(32)
         mat(480) = -rxt(179)*y(32)
         mat(402) = -rxt(180)*y(32)
         mat(212) = -(rxt(232) + rxt(237) + rxt(243)) * y(32)
         mat(200) = -(rxt(233) + rxt(238) + rxt(244)) * y(32)
         mat(144) = -(rxt(234) + rxt(245)) * y(32)

         mat(480) = mat(480) + rxt(170)*y(29)
         mat(241) = rxt(166)*y(27)
         mat(364) = rxt(164)*y(27)
         mat(118) = rxt(163)*y(27)
         mat(542) = rxt(166)*y(18) + rxt(164)*y(23) + rxt(163)*y(24) + rxt(182)*y(33)
         mat(320) = rxt(170)*y(22)
         mat(200) = mat(200) + rxt(182)*y(27)

         mat(194) = -(rxt(181)*y(2) + rxt(182)*y(27) + rxt(183)*y(22) + (rxt(233) &
                      + rxt(238) + rxt(244)) * y(32))
         mat(384) = -rxt(181)*y(33)
         mat(526) = -rxt(182)*y(33)
         mat(462) = -rxt(183)*y(33)
         mat(549) = -(rxt(233) + rxt(238) + rxt(244)) * y(33)

         mat(462) = mat(462) + rxt(185)*y(34)
         mat(348) = rxt(171)*y(29)
         mat(304) = rxt(171)*y(23)
         mat(203) = rxt(185)*y(22)

         mat(204) = -(rxt(184)*y(2) + rxt(185)*y(22) + rxt(186)*y(27) + (rxt(232) &
                      + rxt(237) + rxt(243)) * y(32))
         mat(385) = -rxt(184)*y(34)
         mat(463) = -rxt(185)*y(34)
         mat(527) = -rxt(186)*y(34)
         mat(550) = -(rxt(232) + rxt(237) + rxt(243)) * y(34)

         mat(434) = rxt(173)*y(29)
         mat(305) = rxt(173)*y(10)


      end subroutine nlnmat02

      subroutine nlnmat03( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)


!----------------------------------------------
! ... local variables
!----------------------------------------------

!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------



         mat(302) = rxt(197)*y(37)
         mat(546) = (rxt(234)+rxt(245))*y(39)
         mat(278) = rxt(197)*y(29)
         mat(137) = (rxt(234)+rxt(245))*y(32)

         mat(584) = -(rxt(187)*y(1) + rxt(188)*y(23) + rxt(189)*y(18))
         mat(521) = -rxt(187)*y(36)
         mat(365) = -rxt(188)*y(36)
         mat(242) = -rxt(189)*y(36)

         mat(403) = rxt(190)*y(37) + rxt(200)*y(38)
         mat(340) = rxt(113)*y(38)
         mat(426) = rxt(193)*y(37)
         mat(481) = rxt(191)*y(37) + rxt(199)*y(38)
         mat(321) = (rxt(195)+rxt(196))*y(37)
         mat(298) = rxt(190)*y(2) + rxt(193)*y(9) + rxt(191)*y(22) + (rxt(195) &
                       +rxt(196))*y(29) + 4.000_r8*rxt(198)*y(37) + rxt(222)*y(61)
         mat(131) = rxt(200)*y(2) + rxt(113)*y(3) + rxt(199)*y(22)
         mat(604) = rxt(222)*y(37)

         mat(286) = -(rxt(190)*y(2) + rxt(191)*y(22) + rxt(192)*y(23) + rxt(193)*y(9) &
                      + rxt(194)*y(10) + (rxt(195) + rxt(196) + rxt(197)) * y(29) &
                      + 4._r8*rxt(198)*y(37) + rxt(222)*y(61))
         mat(391) = -rxt(190)*y(37)
         mat(469) = -rxt(191)*y(37)
         mat(353) = -rxt(192)*y(37)
         mat(414) = -rxt(193)*y(37)
         mat(438) = -rxt(194)*y(37)
         mat(309) = -(rxt(195) + rxt(196) + rxt(197)) * y(37)
         mat(592) = -rxt(222)*y(37)

         mat(509) = rxt(187)*y(36)
         mat(391) = mat(391) + rxt(201)*y(39) + rxt(202)*y(40)
         mat(572) = rxt(187)*y(1)
         mat(140) = rxt(201)*y(2)
         mat(96) = rxt(202)*y(2)

         mat(126) = -(rxt(113)*y(3) + rxt(199)*y(22) + rxt(200)*y(2))
         mat(324) = -rxt(113)*y(38)
         mat(459) = -rxt(199)*y(38)
         mat(377) = -rxt(200)*y(38)

         mat(227) = rxt(189)*y(36)
         mat(345) = rxt(188)*y(36)
         mat(568) = rxt(189)*y(18) + rxt(188)*y(23)

         mat(138) = -(rxt(201)*y(2) + (rxt(234) + rxt(245)) * y(32))
         mat(379) = -rxt(201)*y(39)
         mat(548) = -(rxt(234) + rxt(245)) * y(39)

         mat(346) = rxt(192)*y(37)
         mat(282) = rxt(192)*y(23)

         mat(93) = -(rxt(202)*y(2))
         mat(373) = -rxt(202)*y(40)

         mat(430) = rxt(194)*y(37)
         mat(280) = rxt(194)*y(10)

         mat(103) = -((rxt(248) + rxt(257)) * y(2) + rxt(255)*y(4) + rxt(260)*y(57))
         mat(374) = -(rxt(248) + rxt(257)) * y(52)
         mat(260) = -rxt(255)*y(52)
         mat(155) = -rxt(260)*y(52)

         mat(146) = -(rxt(250)*y(8) + rxt(251)*y(9) + rxt(259)*y(57))
         mat(174) = -rxt(250)*y(53)
         mat(405) = -rxt(251)*y(53)
         mat(156) = -rxt(259)*y(53)

         mat(262) = rxt(255)*y(52) + rxt(252)*y(54) + rxt(246)*y(55)
         mat(104) = rxt(255)*y(4)
         mat(70) = rxt(252)*y(4)
         mat(83) = rxt(246)*y(4)

         mat(68) = -((rxt(252) + rxt(253)) * y(4) + rxt(254)*y(2))
         mat(257) = -(rxt(252) + rxt(253)) * y(54)
         mat(370) = -rxt(254)*y(54)

         mat(82) = -(rxt(246)*y(4))
         mat(258) = -rxt(246)*y(55)

         mat(371) = rxt(257)*y(52) + rxt(254)*y(54)
         mat(101) = rxt(257)*y(2)
         mat(69) = rxt(254)*y(2)

         mat(166) = -(rxt(258)*y(57))
         mat(158) = -rxt(258)*y(56)

         mat(382) = rxt(248)*y(52)
         mat(264) = rxt(253)*y(54)
         mat(176) = rxt(250)*y(53)
         mat(407) = rxt(251)*y(53)
         mat(106) = rxt(248)*y(2)
         mat(148) = rxt(250)*y(8) + rxt(251)*y(9)
         mat(71) = rxt(253)*y(4)

         mat(88) = -(rxt(114)*y(4) + rxt(115)*y(2))
         mat(259) = -rxt(114)*y(58)
         mat(372) = -rxt(115)*y(58)

         mat(372) = mat(372) + rxt(248)*y(52)
         mat(102) = rxt(248)*y(2) + .900_r8*rxt(260)*y(57)
         mat(164) = .800_r8*rxt(258)*y(57)
         mat(154) = .900_r8*rxt(260)*y(52) + .800_r8*rxt(258)*y(56)

         mat(157) = -(rxt(258)*y(56) + rxt(259)*y(53) + rxt(260)*y(52))
         mat(165) = -rxt(258)*y(57)
         mat(147) = -rxt(259)*y(57)
         mat(105) = -rxt(260)*y(57)

         mat(50) = -(rxt(212)*y(2) + rxt(213)*y(22))
         mat(369) = -rxt(212)*y(59)
         mat(454) = -rxt(213)*y(59)

         mat(119) = -(rxt(214)*y(22) + rxt(215)*y(4) + rxt(216)*y(1))
         mat(458) = -rxt(214)*y(60)
         mat(261) = -rxt(215)*y(60)
         mat(504) = -rxt(216)*y(60)

         mat(605) = -(rxt(217)*y(22) + rxt(218)*y(4) + rxt(219)*y(1) + rxt(220)*y(10) &
                      + rxt(221)*y(29) + rxt(222)*y(37) + rxt(223)*y(30))
         mat(482) = -rxt(217)*y(61)
         mat(277) = -rxt(218)*y(61)
         mat(522) = -rxt(219)*y(61)
         mat(451) = -rxt(220)*y(61)
         mat(322) = -rxt(221)*y(61)
         mat(299) = -rxt(222)*y(61)
         mat(61) = -rxt(223)*y(61)

         mat(522) = mat(522) + rxt(216)*y(60)
         mat(404) = rxt(212)*y(59)
         mat(277) = mat(277) + rxt(215)*y(60)
         mat(482) = mat(482) + rxt(214)*y(60)
         mat(56) = rxt(212)*y(2)
         mat(125) = rxt(216)*y(1) + rxt(215)*y(4) + rxt(214)*y(22)

         mat(214) = -(rxt(224)*y(22))
         mat(464) = -rxt(224)*y(62)

         mat(505) = rxt(219)*y(61)
         mat(266) = rxt(218)*y(61)
         mat(435) = rxt(220)*y(61)
         mat(464) = mat(464) + rxt(213)*y(59) + rxt(217)*y(61) + (.500_r8*rxt(226) &
                       +rxt(227))*y(65)
         mat(486) = rxt(228)*y(65)
         mat(306) = rxt(221)*y(61)
         mat(58) = rxt(223)*y(61)
         mat(283) = rxt(222)*y(61)
         mat(52) = rxt(213)*y(22)
         mat(588) = rxt(219)*y(1) + rxt(218)*y(4) + rxt(220)*y(10) + rxt(217)*y(22) &
                      + rxt(221)*y(29) + rxt(223)*y(30) + rxt(222)*y(37)
         mat(40) = (.500_r8*rxt(226)+rxt(227))*y(22) + rxt(228)*y(11)

         mat(34) = -(rxt(225)*y(82))
         mat(244) = -rxt(225)*y(63)

         mat(452) = rxt(224)*y(62)
         mat(213) = rxt(224)*y(22)


         mat(243) = rxt(225)*y(63)
         mat(33) = rxt(225)*y(82)

         mat(38) = -((rxt(226) + rxt(227)) * y(22) + rxt(228)*y(11))
         mat(453) = -(rxt(226) + rxt(227)) * y(65)
         mat(483) = -rxt(228)*y(65)
      end subroutine nlnmat03
      subroutine nlnmat_finit( mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat( 1) = lmat( 1)
         mat( 2) = lmat( 2)
         mat( 3) = lmat( 3)
         mat( 4) = lmat( 4)
         mat( 5) = lmat( 5)
         mat( 6) = lmat( 6)
         mat( 7) = lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = lmat( 12)
         mat( 13) = lmat( 13)
         mat( 14) = lmat( 14)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 17) = lmat( 17)
         mat( 18) = lmat( 18)
         mat( 19) = lmat( 19)
         mat( 20) = lmat( 20)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 23) = mat( 23) + lmat( 23)
         mat( 24) = mat( 24) + lmat( 24)
         mat( 25) = mat( 25) + lmat( 25)
         mat( 26) = mat( 26) + lmat( 26)
         mat( 27) = lmat( 27)
         mat( 28) = lmat( 28)
         mat( 29) = lmat( 29)
         mat( 30) = lmat( 30)
         mat( 31) = lmat( 31)
         mat( 32) = lmat( 32)
         mat( 34) = mat( 34) + lmat( 34)
         mat( 35) = lmat( 35)
         mat( 37) = lmat( 37)
         mat( 38) = mat( 38) + lmat( 38)
         mat( 44) = lmat( 44)
         mat( 45) = lmat( 45)
         mat( 46) = lmat( 46)
         mat( 47) = lmat( 47)
         mat( 48) = lmat( 48)
         mat( 49) = lmat( 49)
         mat( 50) = mat( 50) + lmat( 50)
         mat( 51) = lmat( 51)
         mat( 57) = mat( 57) + lmat( 57)
         mat( 59) = mat( 59) + lmat( 59)
         mat( 60) = lmat( 60)
         mat( 62) = mat( 62) + lmat( 62)
         mat( 64) = lmat( 64)
         mat( 65) = lmat( 65)
         mat( 67) = mat( 67) + lmat( 67)
         mat( 68) = mat( 68) + lmat( 68)
         mat( 75) = mat( 75) + lmat( 75)
         mat( 78) = lmat( 78)
         mat( 79) = mat( 79) + lmat( 79)
         mat( 80) = mat( 80) + lmat( 80)
         mat( 81) = lmat( 81)
         mat( 82) = mat( 82) + lmat( 82)
         mat( 83) = mat( 83) + lmat( 83)
         mat( 84) = lmat( 84)
         mat( 85) = lmat( 85)
         mat( 88) = mat( 88) + lmat( 88)
         mat( 93) = mat( 93) + lmat( 93)
         mat( 94) = lmat( 94)
         mat( 95) = lmat( 95)
         mat( 96) = mat( 96) + lmat( 96)
         mat( 98) = lmat( 98)
         mat( 99) = mat( 99) + lmat( 99)
         mat( 100) = lmat( 100)
         mat( 103) = mat( 103) + lmat( 103)
         mat( 112) = mat( 112) + lmat( 112)
         mat( 116) = mat( 116) + lmat( 116)
         mat( 119) = mat( 119) + lmat( 119)
         mat( 126) = mat( 126) + lmat( 126)
         mat( 132) = mat( 132) + lmat( 132)
         mat( 134) = lmat( 134)
         mat( 135) = mat( 135) + lmat( 135)
         mat( 138) = mat( 138) + lmat( 138)
         mat( 142) = mat( 142) + lmat( 142)
         mat( 145) = lmat( 145)
         mat( 146) = mat( 146) + lmat( 146)
         mat( 148) = mat( 148) + lmat( 148)
         mat( 153) = mat( 153) + lmat( 153)
         mat( 157) = mat( 157) + lmat( 157)
         mat( 166) = mat( 166) + lmat( 166)
         mat( 172) = lmat( 172)
         mat( 175) = lmat( 175)
         mat( 177) = mat( 177) + lmat( 177)
         mat( 184) = mat( 184) + lmat( 184)
         mat( 194) = mat( 194) + lmat( 194)
         mat( 198) = mat( 198) + lmat( 198)
         mat( 199) = mat( 199) + lmat( 199)
         mat( 202) = mat( 202) + lmat( 202)
         mat( 203) = mat( 203) + lmat( 203)
         mat( 204) = mat( 204) + lmat( 204)
         mat( 206) = mat( 206) + lmat( 206)
         mat( 208) = lmat( 208)
         mat( 210) = mat( 210) + lmat( 210)
         mat( 211) = mat( 211) + lmat( 211)
         mat( 214) = mat( 214) + lmat( 214)
         mat( 217) = lmat( 217)
         mat( 219) = lmat( 219)
         mat( 220) = mat( 220) + lmat( 220)
         mat( 229) = mat( 229) + lmat( 229)
         mat( 230) = mat( 230) + lmat( 230)
         mat( 246) = lmat( 246)
         mat( 247) = mat( 247) + lmat( 247)
         mat( 249) = mat( 249) + lmat( 249)
         mat( 251) = lmat( 251)
         mat( 252) = mat( 252) + lmat( 252)
         mat( 258) = mat( 258) + lmat( 258)
         mat( 262) = mat( 262) + lmat( 262)
         mat( 263) = lmat( 263)
         mat( 269) = mat( 269) + lmat( 269)
         mat( 270) = mat( 270) + lmat( 270)
         mat( 272) = mat( 272) + lmat( 272)
         mat( 286) = mat( 286) + lmat( 286)
         mat( 290) = mat( 290) + lmat( 290)
         mat( 298) = mat( 298) + lmat( 298)
         mat( 310) = mat( 310) + lmat( 310)
         mat( 313) = mat( 313) + lmat( 313)
         mat( 319) = mat( 319) + lmat( 319)
         mat( 325) = lmat( 325)
         mat( 326) = lmat( 326)
         mat( 327) = lmat( 327)
         mat( 329) = mat( 329) + lmat( 329)
         mat( 330) = mat( 330) + lmat( 330)
         mat( 331) = lmat( 331)
         mat( 332) = mat( 332) + lmat( 332)
         mat( 333) = lmat( 333)
         mat( 335) = mat( 335) + lmat( 335)
         mat( 338) = mat( 338) + lmat( 338)
         mat( 340) = mat( 340) + lmat( 340)
         mat( 344) = mat( 344) + lmat( 344)
         mat( 356) = mat( 356) + lmat( 356)
         mat( 371) = mat( 371) + lmat( 371)
         mat( 381) = lmat( 381)
         mat( 395) = mat( 395) + lmat( 395)
         mat( 406) = lmat( 406)
         mat( 407) = mat( 407) + lmat( 407)
         mat( 408) = mat( 408) + lmat( 408)
         mat( 418) = mat( 418) + lmat( 418)
         mat( 419) = mat( 419) + lmat( 419)
         mat( 431) = mat( 431) + lmat( 431)
         mat( 442) = mat( 442) + lmat( 442)
         mat( 443) = mat( 443) + lmat( 443)
         mat( 444) = mat( 444) + lmat( 444)
         mat( 445) = mat( 445) + lmat( 445)
         mat( 461) = mat( 461) + lmat( 461)
         mat( 465) = mat( 465) + lmat( 465)
         mat( 467) = mat( 467) + lmat( 467)
         mat( 472) = mat( 472) + lmat( 472)
         mat( 476) = mat( 476) + lmat( 476)
         mat( 479) = mat( 479) + lmat( 479)
         mat( 481) = mat( 481) + lmat( 481)
         mat( 485) = mat( 485) + lmat( 485)
         mat( 489) = mat( 489) + lmat( 489)
         mat( 492) = mat( 492) + lmat( 492)
         mat( 493) = mat( 493) + lmat( 493)
         mat( 494) = mat( 494) + lmat( 494)
         mat( 496) = mat( 496) + lmat( 496)
         mat( 502) = mat( 502) + lmat( 502)
         mat( 508) = mat( 508) + lmat( 508)
         mat( 511) = mat( 511) + lmat( 511)
         mat( 513) = mat( 513) + lmat( 513)
         mat( 518) = mat( 518) + lmat( 518)
         mat( 525) = lmat( 525)
         mat( 528) = lmat( 528)
         mat( 534) = mat( 534) + lmat( 534)
         mat( 541) = mat( 541) + lmat( 541)
         mat( 542) = mat( 542) + lmat( 542)
         mat( 551) = lmat( 551)
         mat( 564) = mat( 564) + lmat( 564)
         mat( 565) = mat( 565) + lmat( 565)
         mat( 584) = mat( 584) + lmat( 584)
         mat( 587) = lmat( 587)
         mat( 596) = mat( 596) + lmat( 596)
         mat( 605) = mat( 605) + lmat( 605)
         mat( 109) = 0._r8
         mat( 111) = 0._r8
         mat( 143) = 0._r8
         mat( 160) = 0._r8
         mat( 163) = 0._r8
         mat( 168) = 0._r8
         mat( 169) = 0._r8
         mat( 171) = 0._r8
         mat( 173) = 0._r8
         mat( 179) = 0._r8
         mat( 185) = 0._r8
         mat( 187) = 0._r8
         mat( 192) = 0._r8
         mat( 205) = 0._r8
         mat( 215) = 0._r8
         mat( 232) = 0._r8
         mat( 233) = 0._r8
         mat( 236) = 0._r8
         mat( 239) = 0._r8
         mat( 245) = 0._r8
         mat( 248) = 0._r8
         mat( 250) = 0._r8
         mat( 253) = 0._r8
         mat( 254) = 0._r8
         mat( 268) = 0._r8
         mat( 274) = 0._r8
         mat( 275) = 0._r8
         mat( 281) = 0._r8
         mat( 284) = 0._r8
         mat( 288) = 0._r8
         mat( 294) = 0._r8
         mat( 295) = 0._r8
         mat( 297) = 0._r8
         mat( 307) = 0._r8
         mat( 311) = 0._r8
         mat( 317) = 0._r8
         mat( 318) = 0._r8
         mat( 334) = 0._r8
         mat( 336) = 0._r8
         mat( 341) = 0._r8
         mat( 350) = 0._r8
         mat( 355) = 0._r8
         mat( 366) = 0._r8
         mat( 376) = 0._r8
         mat( 378) = 0._r8
         mat( 380) = 0._r8
         mat( 386) = 0._r8
         mat( 389) = 0._r8
         mat( 393) = 0._r8
         mat( 410) = 0._r8
         mat( 412) = 0._r8
         mat( 416) = 0._r8
         mat( 425) = 0._r8
         mat( 427) = 0._r8
         mat( 432) = 0._r8
         mat( 436) = 0._r8
         mat( 440) = 0._r8
         mat( 448) = 0._r8
         mat( 449) = 0._r8
         mat( 450) = 0._r8
         mat( 471) = 0._r8
         mat( 474) = 0._r8
         mat( 488) = 0._r8
         mat( 490) = 0._r8
         mat( 497) = 0._r8
         mat( 498) = 0._r8
         mat( 499) = 0._r8
         mat( 500) = 0._r8
         mat( 501) = 0._r8
         mat( 507) = 0._r8
         mat( 520) = 0._r8
         mat( 530) = 0._r8
         mat( 533) = 0._r8
         mat( 535) = 0._r8
         mat( 536) = 0._r8
         mat( 537) = 0._r8
         mat( 543) = 0._r8
         mat( 544) = 0._r8
         mat( 553) = 0._r8
         mat( 554) = 0._r8
         mat( 555) = 0._r8
         mat( 557) = 0._r8
         mat( 559) = 0._r8
         mat( 560) = 0._r8
         mat( 562) = 0._r8
         mat( 563) = 0._r8
         mat( 566) = 0._r8
         mat( 567) = 0._r8
         mat( 570) = 0._r8
         mat( 573) = 0._r8
         mat( 574) = 0._r8
         mat( 576) = 0._r8
         mat( 577) = 0._r8
         mat( 578) = 0._r8
         mat( 579) = 0._r8
         mat( 580) = 0._r8
         mat( 582) = 0._r8
         mat( 583) = 0._r8
         mat( 585) = 0._r8
         mat( 590) = 0._r8
         mat( 594) = 0._r8
         mat( 595) = 0._r8
         mat( 600) = 0._r8
         mat( 603) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 3) = mat( 3) - dti
         mat( 4) = mat( 4) - dti
         mat( 5) = mat( 5) - dti
         mat( 6) = mat( 6) - dti
         mat( 7) = mat( 7) - dti
         mat( 8) = mat( 8) - dti
         mat( 9) = mat( 9) - dti
         mat( 10) = mat( 10) - dti
         mat( 11) = mat( 11) - dti
         mat( 12) = mat( 12) - dti
         mat( 13) = mat( 13) - dti
         mat( 14) = mat( 14) - dti
         mat( 15) = mat( 15) - dti
         mat( 16) = mat( 16) - dti
         mat( 17) = mat( 17) - dti
         mat( 20) = mat( 20) - dti
         mat( 23) = mat( 23) - dti
         mat( 26) = mat( 26) - dti
         mat( 28) = mat( 28) - dti
         mat( 30) = mat( 30) - dti
         mat( 34) = mat( 34) - dti
         mat( 38) = mat( 38) - dti
         mat( 44) = mat( 44) - dti
         mat( 50) = mat( 50) - dti
         mat( 57) = mat( 57) - dti
         mat( 62) = mat( 62) - dti
         mat( 68) = mat( 68) - dti
         mat( 75) = mat( 75) - dti
         mat( 82) = mat( 82) - dti
         mat( 88) = mat( 88) - dti
         mat( 93) = mat( 93) - dti
         mat( 103) = mat( 103) - dti
         mat( 112) = mat( 112) - dti
         mat( 119) = mat( 119) - dti
         mat( 126) = mat( 126) - dti
         mat( 132) = mat( 132) - dti
         mat( 138) = mat( 138) - dti
         mat( 146) = mat( 146) - dti
         mat( 157) = mat( 157) - dti
         mat( 166) = mat( 166) - dti
         mat( 177) = mat( 177) - dti
         mat( 184) = mat( 184) - dti
         mat( 194) = mat( 194) - dti
         mat( 204) = mat( 204) - dti
         mat( 214) = mat( 214) - dti
         mat( 220) = mat( 220) - dti
         mat( 230) = mat( 230) - dti
         mat( 247) = mat( 247) - dti
         mat( 269) = mat( 269) - dti
         mat( 286) = mat( 286) - dti
         mat( 310) = mat( 310) - dti
         mat( 330) = mat( 330) - dti
         mat( 356) = mat( 356) - dti
         mat( 395) = mat( 395) - dti
         mat( 419) = mat( 419) - dti
         mat( 444) = mat( 444) - dti
         mat( 476) = mat( 476) - dti
         mat( 496) = mat( 496) - dti
         mat( 518) = mat( 518) - dti
         mat( 541) = mat( 541) - dti
         mat( 565) = mat( 565) - dti
         mat( 584) = mat( 584) - dti
         mat( 605) = mat( 605) - dti
      end subroutine nlnmat_finit
      subroutine nlnmat( mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
      call nlnmat01( mat, y, rxt )
      call nlnmat02( mat, y, rxt )
      call nlnmat03( mat, y, rxt )
      call nlnmat_finit( mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
