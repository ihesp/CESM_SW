




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


         mat(501) = -(rxt(90)*y(2) + rxt(108)*y(3) + rxt(134)*y(21) + rxt(139)*y(22) &
                      + rxt(147)*y(23) + rxt(160)*y(9) + rxt(163)*y(10) + rxt(175) &
                      *y(27) + rxt(202)*y(36))
         mat(582) = -rxt(90)*y(1)
         mat(524) = -rxt(108)*y(1)
         mat(319) = -rxt(134)*y(1)
         mat(441) = -rxt(139)*y(1)
         mat(390) = -rxt(147)*y(1)
         mat(340) = -rxt(160)*y(1)
         mat(365) = -rxt(163)*y(1)
         mat(606) = -rxt(175)*y(1)
         mat(545) = -rxt(202)*y(1)

         mat(582) = mat(582) + rxt(89)*y(4)
         mat(213) = rxt(89)*y(2)

         mat(585) = -(rxt(89)*y(4) + rxt(90)*y(1) + 4._r8*rxt(91)*y(2) + rxt(138) &
                      *y(22) + rxt(145)*y(20) + rxt(146)*y(23) + rxt(149)*y(24) &
                      + rxt(158)*y(9) + (rxt(161) + rxt(162)) * y(10) + rxt(169)*y(11) &
                      + rxt(182)*y(29) + rxt(195)*y(32) + rxt(196)*y(33) + rxt(199) &
                      *y(34) + rxt(205)*y(37) + rxt(215)*y(38) + rxt(216)*y(39) &
                      + rxt(217)*y(40) + rxt(239)*y(18) + (rxt(266) + rxt(267) &
                      ) * y(64) + rxt(273)*y(66))
         mat(215) = -rxt(89)*y(2)
         mat(504) = -rxt(90)*y(2)
         mat(444) = -rxt(138)*y(2)
         mat(257) = -rxt(145)*y(2)
         mat(393) = -rxt(146)*y(2)
         mat(75) = -rxt(149)*y(2)
         mat(343) = -rxt(158)*y(2)
         mat(368) = -(rxt(161) + rxt(162)) * y(2)
         mat(413) = -rxt(169)*y(2)
         mat(306) = -rxt(182)*y(2)
         mat(486) = -rxt(195)*y(2)
         mat(173) = -rxt(196)*y(2)
         mat(185) = -rxt(199)*y(2)
         mat(280) = -rxt(205)*y(2)
         mat(152) = -rxt(215)*y(2)
         mat(105) = -rxt(216)*y(2)
         mat(69) = -rxt(217)*y(2)
         mat(243) = -rxt(239)*y(2)
         mat(87) = -(rxt(266) + rxt(267)) * y(2)
         mat(52) = -rxt(273)*y(2)

         mat(527) = (rxt(103)+rxt(104))*y(4)
         mat(215) = mat(215) + (rxt(103)+rxt(104))*y(3) + rxt(153)*y(8) + rxt(272) &
                      *y(66) + rxt(264)*y(67)
         mat(166) = rxt(153)*y(4) + rxt(154)*y(9) + rxt(155)*y(10) + rxt(269)*y(65)
         mat(343) = mat(343) + rxt(154)*y(8)
         mat(368) = mat(368) + rxt(155)*y(8)
         mat(444) = mat(444) + 2.000_r8*rxt(141)*y(22)
         mat(321) = rxt(137)*y(23)
         mat(393) = mat(393) + rxt(137)*y(21)
         mat(114) = rxt(269)*y(8) + 1.150_r8*rxt(277)*y(69)
         mat(52) = mat(52) + rxt(272)*y(4)
         mat(97) = rxt(264)*y(4)
         mat(122) = rxt(276)*y(69)
         mat(136) = 1.150_r8*rxt(277)*y(65) + rxt(276)*y(68)

         mat(525) = -((rxt(103) + rxt(104)) * y(4) + rxt(105)*y(71) + rxt(108)*y(1) &
                      + rxt(125)*y(59) + rxt(126)*y(60) + rxt(130)*y(20) + rxt(131) &
                      *y(32) + rxt(132)*y(38))
         mat(214) = -(rxt(103) + rxt(104)) * y(3)
         mat(226) = -rxt(105)*y(3)
         mat(502) = -rxt(108)*y(3)
         mat(14) = -rxt(125)*y(3)
         mat(19) = -rxt(126)*y(3)
         mat(256) = -rxt(130)*y(3)
         mat(484) = -rxt(131)*y(3)
         mat(150) = -rxt(132)*y(3)

         mat(214) = mat(214) + rxt(150)*y(70)
         mat(113) = .850_r8*rxt(277)*y(69)
         mat(57) = rxt(150)*y(4)
         mat(135) = .850_r8*rxt(277)*y(65)

         mat(207) = -(rxt(89)*y(2) + rxt(99)*y(6) + rxt(103)*y(3) + rxt(133)*y(21) &
                      + rxt(150)*y(70) + rxt(153)*y(8) + rxt(264)*y(67) + (rxt(271) &
                      + rxt(272)) * y(66) + rxt(274)*y(64))
         mat(568) = -rxt(89)*y(4)
         mat(5) = -rxt(99)*y(4)
         mat(512) = -rxt(103)*y(4)
         mat(308) = -rxt(133)*y(4)
         mat(55) = -rxt(150)*y(4)
         mat(159) = -rxt(153)*y(4)
         mat(93) = -rxt(264)*y(4)
         mat(51) = -(rxt(271) + rxt(272)) * y(4)
         mat(84) = -rxt(274)*y(4)

         mat(490) = 2.000_r8*rxt(90)*y(2) + 2.000_r8*rxt(108)*y(3) + rxt(160)*y(9) &
                      + rxt(163)*y(10) + rxt(139)*y(22) + rxt(134)*y(21) &
                      + 2.000_r8*rxt(147)*y(23) + rxt(175)*y(27) + rxt(202)*y(36)
         mat(568) = mat(568) + 2.000_r8*rxt(90)*y(1) + 2.000_r8*rxt(91)*y(2) + rxt(98) &
                      *y(6) + rxt(161)*y(10) + rxt(138)*y(22) + rxt(169)*y(11) &
                      + rxt(146)*y(23) + rxt(182)*y(29) + rxt(205)*y(37)
         mat(512) = mat(512) + 2.000_r8*rxt(108)*y(1)
         mat(207) = mat(207) + 2.000_r8*rxt(99)*y(6)
         mat(5) = mat(5) + rxt(98)*y(2) + 2.000_r8*rxt(99)*y(4)
         mat(159) = mat(159) + rxt(157)*y(10)
         mat(327) = rxt(160)*y(1) + rxt(270)*y(65)
         mat(352) = rxt(163)*y(1) + rxt(161)*y(2) + rxt(157)*y(8)
         mat(427) = rxt(139)*y(1) + rxt(138)*y(2) + rxt(173)*y(13) + rxt(140)*y(23) &
                      + rxt(184)*y(29)
         mat(398) = rxt(169)*y(2) + rxt(171)*y(23)
         mat(40) = rxt(173)*y(22)
         mat(447) = rxt(241)*y(23)
         mat(308) = mat(308) + rxt(134)*y(1) + rxt(136)*y(23)
         mat(376) = 2.000_r8*rxt(147)*y(1) + rxt(146)*y(2) + rxt(140)*y(22) + rxt(171) &
                      *y(11) + rxt(241)*y(16) + rxt(136)*y(21) + 2.000_r8*rxt(148) &
                      *y(23) + rxt(178)*y(27) + rxt(185)*y(29) + rxt(203)*y(36) &
                      + rxt(207)*y(37)
         mat(593) = rxt(175)*y(1) + rxt(178)*y(23)
         mat(289) = rxt(182)*y(2) + rxt(184)*y(22) + rxt(185)*y(23) + ( &
                      + 2.000_r8*rxt(189)+2.000_r8*rxt(190))*y(29) + (rxt(211) &
                       +rxt(212))*y(37)
         mat(531) = rxt(202)*y(1) + rxt(203)*y(23)
         mat(264) = rxt(205)*y(2) + rxt(207)*y(23) + (rxt(211)+rxt(212))*y(29) &
                      + 2.000_r8*rxt(213)*y(37)
         mat(111) = rxt(270)*y(9)

         mat(7) = -(rxt(92)*y(2) + rxt(93)*y(4) + rxt(95)*y(1))
         mat(551) = -rxt(92)*y(5)
         mat(198) = -rxt(93)*y(5)
         mat(489) = -rxt(95)*y(5)

         mat(506) = rxt(103)*y(4)
         mat(198) = mat(198) + rxt(103)*y(3)

         mat(4) = -(rxt(98)*y(2) + rxt(99)*y(4))
         mat(550) = -rxt(98)*y(6)
         mat(197) = -rxt(99)*y(6)

         mat(488) = rxt(95)*y(5)
         mat(550) = mat(550) + rxt(92)*y(5)
         mat(197) = mat(197) + rxt(93)*y(5)
         mat(6) = rxt(95)*y(1) + rxt(92)*y(2) + rxt(93)*y(4)

         mat(249) = -(rxt(130)*y(3) + rxt(143)*y(22) + rxt(145)*y(2) + rxt(176)*y(27) &
                      + rxt(219)*y(62))
         mat(515) = -rxt(130)*y(20)
         mat(430) = -rxt(143)*y(20)
         mat(571) = -rxt(145)*y(20)
         mat(596) = -rxt(176)*y(20)
         mat(141) = -rxt(219)*y(20)

         mat(310) = rxt(136)*y(23)
         mat(379) = rxt(136)*y(21)

         mat(58) = -((rxt(235) + rxt(236)) * y(22))
         mat(419) = -(rxt(235) + rxt(236)) * y(19)

         mat(554) = rxt(239)*y(18)
         mat(419) = mat(419) + rxt(238)*y(18)
         mat(396) = rxt(237)*y(18)
         mat(228) = rxt(239)*y(2) + rxt(238)*y(22) + rxt(237)*y(11) + rxt(180)*y(27) &
                      + rxt(204)*y(36)
         mat(588) = rxt(180)*y(18)
         mat(529) = rxt(204)*y(18)

         mat(158) = -(rxt(152)*y(22) + rxt(153)*y(4) + rxt(154)*y(9) + (rxt(155) &
                      + rxt(156) + rxt(157)) * y(10) + rxt(269)*y(65))
         mat(423) = -rxt(152)*y(8)
         mat(206) = -rxt(153)*y(8)
         mat(326) = -rxt(154)*y(8)
         mat(349) = -(rxt(155) + rxt(156) + rxt(157)) * y(8)
         mat(110) = -rxt(269)*y(8)

         mat(564) = rxt(273)*y(66) + rxt(151)*y(70)
         mat(206) = mat(206) + rxt(271)*y(66)
         mat(83) = 1.100_r8*rxt(278)*y(69)
         mat(50) = rxt(273)*y(2) + rxt(271)*y(4)
         mat(118) = .200_r8*rxt(276)*y(69)
         mat(54) = rxt(151)*y(2)
         mat(129) = 1.100_r8*rxt(278)*y(64) + .200_r8*rxt(276)*y(68)

         mat(333) = -(rxt(154)*y(8) + rxt(158)*y(2) + rxt(159)*y(23) + rxt(160)*y(1) &
                      + rxt(168)*y(11) + rxt(187)*y(29) + rxt(208)*y(37) + rxt(240) &
                      *y(16) + rxt(270)*y(65))
         mat(161) = -rxt(154)*y(9)
         mat(575) = -rxt(158)*y(9)
         mat(383) = -rxt(159)*y(9)
         mat(494) = -rxt(160)*y(9)
         mat(403) = -rxt(168)*y(9)
         mat(296) = -rxt(187)*y(9)
         mat(270) = -rxt(208)*y(9)
         mat(453) = -rxt(240)*y(9)
         mat(112) = -rxt(270)*y(9)

         mat(575) = mat(575) + rxt(161)*y(10)
         mat(209) = rxt(153)*y(8) + rxt(150)*y(70)
         mat(161) = mat(161) + rxt(153)*y(4) + 2.000_r8*rxt(156)*y(10) + rxt(152) &
                      *y(22)
         mat(358) = rxt(161)*y(2) + 2.000_r8*rxt(156)*y(8)
         mat(434) = rxt(152)*y(8)
         mat(56) = rxt(150)*y(4)

         mat(359) = -((rxt(155) + rxt(156) + rxt(157)) * y(8) + (rxt(161) + rxt(162) &
                      ) * y(2) + rxt(163)*y(1) + rxt(164)*y(11) + rxt(166)*y(22) &
                      + rxt(172)*y(23) + rxt(188)*y(29) + rxt(209)*y(37))
         mat(162) = -(rxt(155) + rxt(156) + rxt(157)) * y(10)
         mat(576) = -(rxt(161) + rxt(162)) * y(10)
         mat(495) = -rxt(163)*y(10)
         mat(404) = -rxt(164)*y(10)
         mat(435) = -rxt(166)*y(10)
         mat(384) = -rxt(172)*y(10)
         mat(297) = -rxt(188)*y(10)
         mat(271) = -rxt(209)*y(10)

         mat(495) = mat(495) + rxt(160)*y(9)
         mat(576) = mat(576) + rxt(158)*y(9) + rxt(169)*y(11)
         mat(334) = rxt(160)*y(1) + rxt(158)*y(2) + 2.000_r8*rxt(168)*y(11) + rxt(240) &
                      *y(16) + rxt(159)*y(23) + rxt(187)*y(29) + rxt(208)*y(37)
         mat(435) = mat(435) + rxt(170)*y(11) + rxt(173)*y(13)
         mat(404) = mat(404) + rxt(169)*y(2) + 2.000_r8*rxt(168)*y(9) + rxt(170)*y(22) &
                      + rxt(171)*y(23)
         mat(42) = rxt(173)*y(22)
         mat(454) = rxt(240)*y(9)
         mat(384) = mat(384) + rxt(159)*y(9) + rxt(171)*y(11)
         mat(297) = mat(297) + rxt(187)*y(9)
         mat(271) = mat(271) + rxt(208)*y(9)

         mat(438) = -(rxt(138)*y(2) + rxt(139)*y(1) + rxt(140)*y(23) + (4._r8*rxt(141) &
                      + 4._r8*rxt(142)) * y(22) + rxt(143)*y(20) + rxt(144)*y(24) &
                      + rxt(152)*y(8) + rxt(166)*y(10) + rxt(167)*y(12) + rxt(170) &
                      *y(11) + rxt(173)*y(13) + (rxt(183) + rxt(184)) * y(29) + rxt(194) &
                      *y(32) + rxt(198)*y(33) + rxt(200)*y(34) + rxt(206)*y(37) &
                      + rxt(214)*y(38) + (rxt(235) + rxt(236)) * y(19) + rxt(238) &
                      *y(18) + rxt(242)*y(17))
         mat(579) = -rxt(138)*y(22)
         mat(498) = -rxt(139)*y(22)
         mat(387) = -rxt(140)*y(22)
         mat(253) = -rxt(143)*y(22)
         mat(73) = -rxt(144)*y(22)
         mat(164) = -rxt(152)*y(22)
         mat(362) = -rxt(166)*y(22)
         mat(195) = -rxt(167)*y(22)
         mat(407) = -rxt(170)*y(22)
         mat(45) = -rxt(173)*y(22)
         mat(300) = -(rxt(183) + rxt(184)) * y(22)
         mat(480) = -rxt(194)*y(22)
         mat(171) = -rxt(198)*y(22)
         mat(183) = -rxt(200)*y(22)
         mat(274) = -rxt(206)*y(22)
         mat(149) = -rxt(214)*y(22)
         mat(61) = -(rxt(235) + rxt(236)) * y(22)
         mat(238) = -rxt(238)*y(22)
         mat(37) = -rxt(242)*y(22)

         mat(498) = mat(498) + rxt(134)*y(21) + rxt(147)*y(23)
         mat(579) = mat(579) + rxt(145)*y(20) + rxt(239)*y(18) + rxt(146)*y(23) &
                      + rxt(149)*y(24) + rxt(195)*y(32) + rxt(196)*y(33) + rxt(215) &
                      *y(38) + rxt(216)*y(39)
         mat(521) = rxt(130)*y(20) + 2.000_r8*rxt(105)*y(71) + rxt(131)*y(32) &
                      + rxt(132)*y(38)
         mat(253) = mat(253) + rxt(145)*y(2) + rxt(130)*y(3)
         mat(337) = rxt(159)*y(23)
         mat(407) = mat(407) + rxt(171)*y(23)
         mat(238) = mat(238) + rxt(239)*y(2)
         mat(316) = rxt(134)*y(1) + 2.000_r8*rxt(135)*y(23)
         mat(387) = mat(387) + rxt(147)*y(1) + rxt(146)*y(2) + rxt(159)*y(9) &
                      + rxt(171)*y(11) + 2.000_r8*rxt(135)*y(21) + rxt(179)*y(27)
         mat(73) = mat(73) + rxt(149)*y(2)
         mat(224) = 2.000_r8*rxt(105)*y(3) + rxt(218)*y(62)
         mat(603) = rxt(179)*y(23)
         mat(480) = mat(480) + rxt(195)*y(2) + rxt(131)*y(3)
         mat(171) = mat(171) + rxt(196)*y(2)
         mat(149) = mat(149) + rxt(215)*y(2) + rxt(132)*y(3)
         mat(102) = rxt(216)*y(2)
         mat(144) = rxt(218)*y(71)


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


         mat(406) = -(rxt(164)*y(10) + rxt(168)*y(9) + rxt(169)*y(2) + rxt(170)*y(22) &
                      + rxt(171)*y(23) + rxt(237)*y(18))
         mat(361) = -rxt(164)*y(11)
         mat(336) = -rxt(168)*y(11)
         mat(578) = -rxt(169)*y(11)
         mat(437) = -rxt(170)*y(11)
         mat(386) = -rxt(171)*y(11)
         mat(237) = -rxt(237)*y(11)

         mat(497) = rxt(163)*y(10)
         mat(578) = mat(578) + rxt(162)*y(10) + rxt(199)*y(34) + rxt(217)*y(40)
         mat(361) = mat(361) + rxt(163)*y(1) + rxt(162)*y(2)
         mat(437) = mat(437) + rxt(167)*y(12) + rxt(200)*y(34)
         mat(194) = rxt(167)*y(22) + rxt(221)*y(62)
         mat(602) = rxt(201)*y(34)
         mat(182) = rxt(199)*y(2) + rxt(200)*y(22) + rxt(201)*y(27)
         mat(67) = rxt(217)*y(2)
         mat(143) = rxt(221)*y(12)

         mat(189) = -(rxt(167)*y(22) + rxt(221)*y(62))
         mat(426) = -rxt(167)*y(12)
         mat(139) = -rxt(221)*y(12)

         mat(351) = rxt(166)*y(22)
         mat(426) = mat(426) + rxt(166)*y(10)
         mat(397) = rxt(237)*y(18)
         mat(230) = rxt(237)*y(11)
         mat(470) = (rxt(250)+rxt(255)+rxt(261))*y(34)
         mat(178) = (rxt(250)+rxt(255)+rxt(261))*y(32)

         mat(39) = -(rxt(173)*y(22))
         mat(418) = -rxt(173)*y(13)

         mat(346) = rxt(172)*y(23)
         mat(371) = rxt(172)*y(10)


         mat(345) = rxt(164)*y(11)
         mat(395) = rxt(164)*y(10)

         mat(458) = -(rxt(186)*y(29) + rxt(240)*y(9) + rxt(241)*y(23))
         mat(301) = -rxt(186)*y(16)
         mat(338) = -rxt(240)*y(16)
         mat(388) = -rxt(241)*y(16)

         mat(439) = rxt(242)*y(17)
         mat(38) = rxt(242)*y(22)

         mat(33) = -(rxt(242)*y(22))
         mat(417) = -rxt(242)*y(17)

         mat(446) = rxt(241)*y(23)
         mat(370) = rxt(241)*y(16)

         mat(232) = -(rxt(180)*y(27) + rxt(204)*y(36) + rxt(237)*y(11) + rxt(238) &
                      *y(22) + rxt(239)*y(2))
         mat(595) = -rxt(180)*y(18)
         mat(533) = -rxt(204)*y(18)
         mat(400) = -rxt(237)*y(18)
         mat(429) = -rxt(238)*y(18)
         mat(570) = -rxt(239)*y(18)

         mat(328) = rxt(240)*y(16)
         mat(449) = rxt(240)*y(9) + rxt(186)*y(29)
         mat(291) = rxt(186)*y(16)

         mat(311) = -(rxt(133)*y(4) + rxt(134)*y(1) + (rxt(135) + rxt(136) + rxt(137) &
                      ) * y(23))
         mat(208) = -rxt(133)*y(21)
         mat(493) = -rxt(134)*y(21)
         mat(382) = -(rxt(135) + rxt(136) + rxt(137)) * y(21)

         mat(574) = rxt(145)*y(20) + rxt(138)*y(22)
         mat(516) = rxt(130)*y(20)
         mat(250) = rxt(145)*y(2) + rxt(130)*y(3) + rxt(143)*y(22) + rxt(176)*y(27) &
                      + rxt(219)*y(62)
         mat(59) = rxt(235)*y(22)
         mat(160) = rxt(152)*y(22)
         mat(433) = rxt(138)*y(2) + rxt(143)*y(20) + rxt(235)*y(19) + rxt(152)*y(8) &
                      + rxt(238)*y(18)
         mat(234) = rxt(238)*y(22)
         mat(598) = rxt(176)*y(20)
         mat(142) = rxt(219)*y(20)

         mat(385) = -((rxt(135) + rxt(136) + rxt(137)) * y(21) + rxt(140)*y(22) &
                      + rxt(146)*y(2) + rxt(147)*y(1) + 4._r8*rxt(148)*y(23) + rxt(159) &
                      *y(9) + rxt(171)*y(11) + rxt(172)*y(10) + (rxt(178) + rxt(179) &
                      ) * y(27) + rxt(185)*y(29) + rxt(203)*y(36) + rxt(207)*y(37) &
                      + rxt(241)*y(16))
         mat(314) = -(rxt(135) + rxt(136) + rxt(137)) * y(23)
         mat(436) = -rxt(140)*y(23)
         mat(577) = -rxt(146)*y(23)
         mat(496) = -rxt(147)*y(23)
         mat(335) = -rxt(159)*y(23)
         mat(405) = -rxt(171)*y(23)
         mat(360) = -rxt(172)*y(23)
         mat(601) = -(rxt(178) + rxt(179)) * y(23)
         mat(298) = -rxt(185)*y(23)
         mat(540) = -rxt(203)*y(23)
         mat(272) = -rxt(207)*y(23)
         mat(455) = -rxt(241)*y(23)

         mat(496) = mat(496) + rxt(139)*y(22)
         mat(577) = mat(577) + rxt(239)*y(18) + rxt(149)*y(24)
         mat(211) = rxt(133)*y(21)
         mat(60) = rxt(236)*y(22)
         mat(335) = mat(335) + rxt(240)*y(16)
         mat(436) = mat(436) + rxt(139)*y(1) + rxt(236)*y(19) + rxt(170)*y(11) &
                      + rxt(144)*y(24) + rxt(183)*y(29) + rxt(206)*y(37)
         mat(405) = mat(405) + rxt(170)*y(22) + rxt(237)*y(18)
         mat(455) = mat(455) + rxt(240)*y(9) + rxt(186)*y(29)
         mat(236) = rxt(239)*y(2) + rxt(237)*y(11) + rxt(180)*y(27) + rxt(204)*y(36)
         mat(314) = mat(314) + rxt(133)*y(4)
         mat(72) = rxt(149)*y(2) + rxt(144)*y(22) + rxt(177)*y(27)
         mat(601) = mat(601) + rxt(180)*y(18) + rxt(177)*y(24)
         mat(298) = mat(298) + rxt(183)*y(22) + rxt(186)*y(16)
         mat(540) = mat(540) + rxt(204)*y(18)
         mat(272) = mat(272) + rxt(206)*y(22)

         mat(70) = -(rxt(144)*y(22) + rxt(149)*y(2) + rxt(177)*y(27))
         mat(420) = -rxt(144)*y(24)
         mat(556) = -rxt(149)*y(24)
         mat(589) = -rxt(177)*y(24)

         mat(420) = mat(420) + 2.000_r8*rxt(142)*y(22)
         mat(372) = 2.000_r8*rxt(148)*y(23)

         mat(219) = -(rxt(105)*y(3) + rxt(218)*y(62))
         mat(513) = -rxt(105)*y(71)
         mat(140) = -rxt(218)*y(71)

         mat(248) = rxt(143)*y(22)
         mat(428) = rxt(143)*y(20) + 2.000_r8*rxt(141)*y(22) + rxt(167)*y(12) &
                      + rxt(173)*y(13) + rxt(242)*y(17) + rxt(238)*y(18) + rxt(140) &
                      *y(23) + rxt(144)*y(24) + rxt(194)*y(32) + rxt(198)*y(33) &
                      + rxt(214)*y(38)
         mat(190) = rxt(167)*y(22)
         mat(41) = rxt(173)*y(22)
         mat(34) = rxt(242)*y(22)
         mat(231) = rxt(238)*y(22)
         mat(309) = rxt(137)*y(23)
         mat(377) = rxt(140)*y(22) + rxt(137)*y(21)
         mat(71) = rxt(144)*y(22)
         mat(471) = rxt(194)*y(22) + (rxt(251)+rxt(256)+rxt(262))*y(33) + (rxt(252) &
                       +rxt(263))*y(39)
         mat(169) = rxt(198)*y(22) + (rxt(251)+rxt(256)+rxt(262))*y(32)
         mat(147) = rxt(214)*y(22)
         mat(100) = (rxt(252)+rxt(263))*y(32)

         mat(610) = -(rxt(175)*y(1) + rxt(176)*y(20) + rxt(177)*y(24) + (rxt(178) &
                      + rxt(179)) * y(23) + rxt(180)*y(18) + rxt(197)*y(33) + rxt(201) &
                      *y(34))
         mat(505) = -rxt(175)*y(27)
         mat(258) = -rxt(176)*y(27)
         mat(76) = -rxt(177)*y(27)
         mat(394) = -(rxt(178) + rxt(179)) * y(27)
         mat(244) = -rxt(180)*y(27)
         mat(174) = -rxt(197)*y(27)
         mat(186) = -rxt(201)*y(27)

         mat(586) = rxt(182)*y(29) + rxt(195)*y(32)
         mat(528) = rxt(131)*y(32) + rxt(126)*y(60)
         mat(344) = rxt(187)*y(29)
         mat(445) = rxt(183)*y(29) + rxt(194)*y(32)
         mat(464) = rxt(186)*y(29)
         mat(307) = rxt(182)*y(2) + rxt(187)*y(9) + rxt(183)*y(22) + rxt(186)*y(16) + ( &
                      + 4.000_r8*rxt(189)+2.000_r8*rxt(191))*y(29) + rxt(211)*y(37)
         mat(487) = rxt(195)*y(2) + rxt(131)*y(3) + rxt(194)*y(22)
         mat(281) = rxt(211)*y(29)
         mat(20) = rxt(126)*y(3)


         mat(587) = rxt(201)*y(34)
         mat(284) = 2.000_r8*rxt(190)*y(29)
         mat(465) = (rxt(251)+rxt(256)+rxt(262))*y(33) + (rxt(250)+rxt(255)+rxt(261)) &
                      *y(34)
         mat(167) = (rxt(251)+rxt(256)+rxt(262))*y(32)
         mat(175) = rxt(201)*y(27) + (rxt(250)+rxt(255)+rxt(261))*y(32)

         mat(294) = -(rxt(182)*y(2) + (rxt(183) + rxt(184)) * y(22) + rxt(185)*y(23) &
                      + rxt(186)*y(16) + rxt(187)*y(9) + rxt(188)*y(10) + (4._r8*rxt(189) &
                      + 4._r8*rxt(190) + 4._r8*rxt(191) + 4._r8*rxt(192)) * y(29) &
                      + (rxt(210) + rxt(211) + rxt(212)) * y(37))
         mat(573) = -rxt(182)*y(29)
         mat(432) = -(rxt(183) + rxt(184)) * y(29)
         mat(381) = -rxt(185)*y(29)
         mat(451) = -rxt(186)*y(29)
         mat(331) = -rxt(187)*y(29)
         mat(356) = -rxt(188)*y(29)
         mat(268) = -(rxt(210) + rxt(211) + rxt(212)) * y(29)

         mat(492) = rxt(175)*y(27)
         mat(573) = mat(573) + rxt(196)*y(33) + rxt(199)*y(34)
         mat(432) = mat(432) + rxt(198)*y(33)
         mat(381) = mat(381) + rxt(179)*y(27)
         mat(597) = rxt(175)*y(1) + rxt(179)*y(23) + rxt(197)*y(33)
         mat(170) = rxt(196)*y(2) + rxt(198)*y(22) + rxt(197)*y(27)
         mat(180) = rxt(199)*y(2)


         mat(283) = 2.000_r8*rxt(191)*y(29) + rxt(210)*y(37)
         mat(259) = rxt(210)*y(29)


         mat(282) = 2.000_r8*rxt(192)*y(29)

         mat(482) = -(rxt(131)*y(3) + rxt(194)*y(22) + rxt(195)*y(2) + (rxt(250) &
                      + rxt(255) + rxt(261)) * y(34) + (rxt(251) + rxt(256) + rxt(262) &
                      ) * y(33) + (rxt(252) + rxt(263)) * y(39))
         mat(523) = -rxt(131)*y(32)
         mat(440) = -rxt(194)*y(32)
         mat(581) = -rxt(195)*y(32)
         mat(184) = -(rxt(250) + rxt(255) + rxt(261)) * y(32)
         mat(172) = -(rxt(251) + rxt(256) + rxt(262)) * y(32)
         mat(103) = -(rxt(252) + rxt(263)) * y(32)

         mat(255) = rxt(176)*y(27)
         mat(440) = mat(440) + rxt(184)*y(29)
         mat(240) = rxt(180)*y(27)
         mat(389) = rxt(178)*y(27)
         mat(74) = rxt(177)*y(27)
         mat(605) = rxt(176)*y(20) + rxt(180)*y(18) + rxt(178)*y(23) + rxt(177)*y(24) &
                      + rxt(197)*y(33)
         mat(302) = rxt(184)*y(22)
         mat(172) = mat(172) + rxt(197)*y(27)

         mat(168) = -(rxt(196)*y(2) + rxt(197)*y(27) + rxt(198)*y(22) + (rxt(251) &
                      + rxt(256) + rxt(262)) * y(32))
         mat(565) = -rxt(196)*y(33)
         mat(590) = -rxt(197)*y(33)
         mat(424) = -rxt(198)*y(33)
         mat(468) = -(rxt(251) + rxt(256) + rxt(262)) * y(33)

         mat(424) = mat(424) + rxt(200)*y(34)
         mat(375) = rxt(185)*y(29)
         mat(286) = rxt(185)*y(23)
         mat(176) = rxt(200)*y(22)


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


         mat(177) = -(rxt(199)*y(2) + rxt(200)*y(22) + rxt(201)*y(27) + (rxt(250) &
                      + rxt(255) + rxt(261)) * y(32))
         mat(566) = -rxt(199)*y(34)
         mat(425) = -rxt(200)*y(34)
         mat(591) = -rxt(201)*y(34)
         mat(469) = -(rxt(250) + rxt(255) + rxt(261)) * y(34)

         mat(350) = rxt(188)*y(29)
         mat(287) = rxt(188)*y(10)


         mat(285) = rxt(212)*y(37)
         mat(466) = (rxt(252)+rxt(263))*y(39)
         mat(260) = rxt(212)*y(29)
         mat(98) = (rxt(252)+rxt(263))*y(32)

         mat(547) = -(rxt(202)*y(1) + rxt(203)*y(23) + rxt(204)*y(18))
         mat(503) = -rxt(202)*y(36)
         mat(392) = -rxt(203)*y(36)
         mat(242) = -rxt(204)*y(36)

         mat(584) = rxt(205)*y(37) + rxt(215)*y(38)
         mat(526) = rxt(132)*y(38)
         mat(342) = rxt(208)*y(37)
         mat(443) = rxt(206)*y(37) + rxt(214)*y(38)
         mat(305) = (rxt(210)+rxt(211))*y(37)
         mat(279) = rxt(205)*y(2) + rxt(208)*y(9) + rxt(206)*y(22) + (rxt(210) &
                       +rxt(211))*y(29) + 4.000_r8*rxt(213)*y(37)
         mat(151) = rxt(215)*y(2) + rxt(132)*y(3) + rxt(214)*y(22)

         mat(267) = -(rxt(205)*y(2) + rxt(206)*y(22) + rxt(207)*y(23) + rxt(208)*y(9) &
                      + rxt(209)*y(10) + (rxt(210) + rxt(211) + rxt(212)) * y(29) &
                      + 4._r8*rxt(213)*y(37))
         mat(572) = -rxt(205)*y(37)
         mat(431) = -rxt(206)*y(37)
         mat(380) = -rxt(207)*y(37)
         mat(330) = -rxt(208)*y(37)
         mat(355) = -rxt(209)*y(37)
         mat(293) = -(rxt(210) + rxt(211) + rxt(212)) * y(37)

         mat(491) = rxt(202)*y(36)
         mat(572) = mat(572) + rxt(216)*y(39) + rxt(217)*y(40)
         mat(535) = rxt(202)*y(1)
         mat(101) = rxt(216)*y(2)
         mat(65) = rxt(217)*y(2)

         mat(146) = -(rxt(132)*y(3) + rxt(214)*y(22) + rxt(215)*y(2))
         mat(510) = -rxt(132)*y(38)
         mat(422) = -rxt(214)*y(38)
         mat(563) = -rxt(215)*y(38)

         mat(229) = rxt(204)*y(36)
         mat(374) = rxt(203)*y(36)
         mat(530) = rxt(204)*y(18) + rxt(203)*y(23)

         mat(99) = -(rxt(216)*y(2) + (rxt(252) + rxt(263)) * y(32))
         mat(559) = -rxt(216)*y(39)
         mat(467) = -(rxt(252) + rxt(263)) * y(39)

         mat(373) = rxt(207)*y(37)
         mat(262) = rxt(207)*y(23)

         mat(62) = -(rxt(217)*y(2))
         mat(555) = -rxt(217)*y(40)

         mat(347) = rxt(209)*y(37)
         mat(261) = rxt(209)*y(10)

         mat(78) = -((rxt(266) + rxt(267)) * y(2) + rxt(274)*y(4) + rxt(278)*y(69))
         mat(557) = -(rxt(266) + rxt(267)) * y(64)
         mat(201) = -rxt(274)*y(64)
         mat(124) = -rxt(278)*y(64)

         mat(107) = -(rxt(269)*y(8) + rxt(270)*y(9) + rxt(277)*y(69))
         mat(155) = -rxt(269)*y(65)
         mat(323) = -rxt(270)*y(65)
         mat(126) = -rxt(277)*y(65)

         mat(203) = rxt(274)*y(64) + rxt(271)*y(66) + rxt(264)*y(67)
         mat(80) = rxt(274)*y(4)
         mat(48) = rxt(271)*y(4)
         mat(90) = rxt(264)*y(4)

         mat(46) = -((rxt(271) + rxt(272)) * y(4) + rxt(273)*y(2))
         mat(199) = -(rxt(271) + rxt(272)) * y(66)
         mat(552) = -rxt(273)*y(66)

         mat(89) = -(rxt(264)*y(4))
         mat(202) = -rxt(264)*y(67)

         mat(558) = rxt(267)*y(64) + rxt(273)*y(66)
         mat(79) = rxt(267)*y(2)
         mat(47) = rxt(273)*y(2)

         mat(116) = -(rxt(276)*y(69))
         mat(127) = -rxt(276)*y(68)

         mat(561) = rxt(266)*y(64)
         mat(204) = rxt(272)*y(66)
         mat(156) = rxt(269)*y(65)
         mat(324) = rxt(270)*y(65)
         mat(81) = rxt(266)*y(2)
         mat(108) = rxt(269)*y(8) + rxt(270)*y(9)
         mat(49) = rxt(272)*y(4)

         mat(53) = -(rxt(150)*y(4) + rxt(151)*y(2))
         mat(200) = -rxt(150)*y(70)
         mat(553) = -rxt(151)*y(70)

         mat(553) = mat(553) + rxt(266)*y(64)
         mat(77) = rxt(266)*y(2) + .900_r8*rxt(278)*y(69)
         mat(115) = .800_r8*rxt(276)*y(69)
         mat(123) = .900_r8*rxt(278)*y(64) + .800_r8*rxt(276)*y(68)

         mat(128) = -(rxt(276)*y(68) + rxt(277)*y(65) + rxt(278)*y(64))
         mat(117) = -rxt(276)*y(69)
         mat(109) = -rxt(277)*y(69)
         mat(82) = -rxt(278)*y(69)

         mat(12) = -(rxt(125)*y(3))
         mat(507) = -rxt(125)*y(59)

         mat(17) = -(rxt(126)*y(3))
         mat(508) = -rxt(126)*y(60)


         mat(245) = rxt(219)*y(62)
         mat(187) = rxt(221)*y(62)
         mat(216) = rxt(218)*y(62)
         mat(137) = rxt(219)*y(20) + rxt(221)*y(12) + rxt(218)*y(71)

         mat(138) = -(rxt(218)*y(71) + rxt(219)*y(20) + rxt(221)*y(12))
         mat(217) = -rxt(218)*y(62)
         mat(246) = -rxt(219)*y(62)
         mat(188) = -rxt(221)*y(62)

         mat(509) = 2.000_r8*rxt(125)*y(59) + rxt(126)*y(60)
         mat(13) = 2.000_r8*rxt(125)*y(3)
         mat(18) = rxt(126)*y(3)


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
         mat( 4) = mat( 4) + lmat( 4)
         mat( 5) = mat( 5) + lmat( 5)
         mat( 6) = mat( 6) + lmat( 6)
         mat( 7) = mat( 7) + lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = mat( 12) + lmat( 12)
         mat( 13) = mat( 13) + lmat( 13)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 17) = mat( 17) + lmat( 17)
         mat( 18) = mat( 18) + lmat( 18)
         mat( 20) = mat( 20) + lmat( 20)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 23) = lmat( 23)
         mat( 24) = lmat( 24)
         mat( 25) = lmat( 25)
         mat( 26) = lmat( 26)
         mat( 27) = lmat( 27)
         mat( 28) = lmat( 28)
         mat( 29) = lmat( 29)
         mat( 30) = lmat( 30)
         mat( 31) = lmat( 31)
         mat( 32) = lmat( 32)
         mat( 33) = mat( 33) + lmat( 33)
         mat( 35) = lmat( 35)
         mat( 36) = lmat( 36)
         mat( 37) = mat( 37) + lmat( 37)
         mat( 39) = mat( 39) + lmat( 39)
         mat( 42) = mat( 42) + lmat( 42)
         mat( 43) = lmat( 43)
         mat( 44) = lmat( 44)
         mat( 45) = mat( 45) + lmat( 45)
         mat( 46) = mat( 46) + lmat( 46)
         mat( 53) = mat( 53) + lmat( 53)
         mat( 58) = mat( 58) + lmat( 58)
         mat( 62) = mat( 62) + lmat( 62)
         mat( 63) = lmat( 63)
         mat( 64) = lmat( 64)
         mat( 65) = mat( 65) + lmat( 65)
         mat( 66) = lmat( 66)
         mat( 67) = mat( 67) + lmat( 67)
         mat( 68) = lmat( 68)
         mat( 70) = mat( 70) + lmat( 70)
         mat( 73) = mat( 73) + lmat( 73)
         mat( 78) = mat( 78) + lmat( 78)
         mat( 88) = lmat( 88)
         mat( 89) = mat( 89) + lmat( 89)
         mat( 90) = mat( 90) + lmat( 90)
         mat( 91) = lmat( 91)
         mat( 92) = lmat( 92)
         mat( 99) = mat( 99) + lmat( 99)
         mat( 102) = mat( 102) + lmat( 102)
         mat( 104) = lmat( 104)
         mat( 107) = mat( 107) + lmat( 107)
         mat( 108) = mat( 108) + lmat( 108)
         mat( 112) = mat( 112) + lmat( 112)
         mat( 116) = mat( 116) + lmat( 116)
         mat( 128) = mat( 128) + lmat( 128)
         mat( 137) = mat( 137) + lmat( 137)
         mat( 138) = mat( 138) + lmat( 138)
         mat( 145) = lmat( 145)
         mat( 146) = mat( 146) + lmat( 146)
         mat( 148) = lmat( 148)
         mat( 151) = mat( 151) + lmat( 151)
         mat( 153) = lmat( 153)
         mat( 157) = lmat( 157)
         mat( 158) = mat( 158) + lmat( 158)
         mat( 168) = mat( 168) + lmat( 168)
         mat( 171) = mat( 171) + lmat( 171)
         mat( 174) = mat( 174) + lmat( 174)
         mat( 176) = mat( 176) + lmat( 176)
         mat( 177) = mat( 177) + lmat( 177)
         mat( 178) = mat( 178) + lmat( 178)
         mat( 180) = mat( 180) + lmat( 180)
         mat( 181) = lmat( 181)
         mat( 182) = mat( 182) + lmat( 182)
         mat( 186) = mat( 186) + lmat( 186)
         mat( 189) = mat( 189) + lmat( 189)
         mat( 193) = lmat( 193)
         mat( 195) = mat( 195) + lmat( 195)
         mat( 202) = mat( 202) + lmat( 202)
         mat( 203) = mat( 203) + lmat( 203)
         mat( 205) = lmat( 205)
         mat( 207) = mat( 207) + lmat( 207)
         mat( 214) = mat( 214) + lmat( 214)
         mat( 215) = mat( 215) + lmat( 215)
         mat( 219) = mat( 219) + lmat( 219)
         mat( 220) = lmat( 220)
         mat( 221) = lmat( 221)
         mat( 224) = mat( 224) + lmat( 224)
         mat( 226) = mat( 226) + lmat( 226)
         mat( 227) = lmat( 227)
         mat( 228) = mat( 228) + lmat( 228)
         mat( 232) = mat( 232) + lmat( 232)
         mat( 233) = lmat( 233)
         mat( 234) = mat( 234) + lmat( 234)
         mat( 249) = mat( 249) + lmat( 249)
         mat( 267) = mat( 267) + lmat( 267)
         mat( 279) = mat( 279) + lmat( 279)
         mat( 280) = mat( 280) + lmat( 280)
         mat( 294) = mat( 294) + lmat( 294)
         mat( 306) = mat( 306) + lmat( 306)
         mat( 307) = mat( 307) + lmat( 307)
         mat( 311) = mat( 311) + lmat( 311)
         mat( 324) = mat( 324) + lmat( 324)
         mat( 325) = lmat( 325)
         mat( 326) = mat( 326) + lmat( 326)
         mat( 333) = mat( 333) + lmat( 333)
         mat( 343) = mat( 343) + lmat( 343)
         mat( 351) = mat( 351) + lmat( 351)
         mat( 358) = mat( 358) + lmat( 358)
         mat( 359) = mat( 359) + lmat( 359)
         mat( 362) = mat( 362) + lmat( 362)
         mat( 368) = mat( 368) + lmat( 368)
         mat( 372) = mat( 372) + lmat( 372)
         mat( 385) = mat( 385) + lmat( 385)
         mat( 397) = mat( 397) + lmat( 397)
         mat( 398) = mat( 398) + lmat( 398)
         mat( 403) = mat( 403) + lmat( 403)
         mat( 404) = mat( 404) + lmat( 404)
         mat( 406) = mat( 406) + lmat( 406)
         mat( 413) = mat( 413) + lmat( 413)
         mat( 415) = lmat( 415)
         mat( 416) = lmat( 416)
         mat( 428) = mat( 428) + lmat( 428)
         mat( 436) = mat( 436) + lmat( 436)
         mat( 438) = mat( 438) + lmat( 438)
         mat( 439) = mat( 439) + lmat( 439)
         mat( 443) = mat( 443) + lmat( 443)
         mat( 445) = mat( 445) + lmat( 445)
         mat( 458) = mat( 458) + lmat( 458)
         mat( 475) = lmat( 475)
         mat( 482) = mat( 482) + lmat( 482)
         mat( 487) = mat( 487) + lmat( 487)
         mat( 488) = mat( 488) + lmat( 488)
         mat( 490) = mat( 490) + lmat( 490)
         mat( 501) = mat( 501) + lmat( 501)
         mat( 502) = mat( 502) + lmat( 502)
         mat( 504) = mat( 504) + lmat( 504)
         mat( 507) = mat( 507) + lmat( 507)
         mat( 508) = mat( 508) + lmat( 508)
         mat( 509) = mat( 509) + lmat( 509)
         mat( 512) = mat( 512) + lmat( 512)
         mat( 514) = lmat( 514)
         mat( 515) = mat( 515) + lmat( 515)
         mat( 516) = mat( 516) + lmat( 516)
         mat( 517) = lmat( 517)
         mat( 519) = lmat( 519)
         mat( 521) = mat( 521) + lmat( 521)
         mat( 522) = lmat( 522)
         mat( 525) = mat( 525) + lmat( 525)
         mat( 526) = mat( 526) + lmat( 526)
         mat( 527) = mat( 527) + lmat( 527)
         mat( 528) = mat( 528) + lmat( 528)
         mat( 547) = mat( 547) + lmat( 547)
         mat( 558) = mat( 558) + lmat( 558)
         mat( 562) = lmat( 562)
         mat( 585) = mat( 585) + lmat( 585)
         mat( 588) = mat( 588) + lmat( 588)
         mat( 601) = mat( 601) + lmat( 601)
         mat( 604) = lmat( 604)
         mat( 605) = mat( 605) + lmat( 605)
         mat( 608) = lmat( 608)
         mat( 610) = mat( 610) + lmat( 610)
         mat( 85) = 0._r8
         mat( 86) = 0._r8
         mat( 94) = 0._r8
         mat( 95) = 0._r8
         mat( 96) = 0._r8
         mat( 106) = 0._r8
         mat( 119) = 0._r8
         mat( 120) = 0._r8
         mat( 121) = 0._r8
         mat( 125) = 0._r8
         mat( 130) = 0._r8
         mat( 131) = 0._r8
         mat( 132) = 0._r8
         mat( 133) = 0._r8
         mat( 134) = 0._r8
         mat( 154) = 0._r8
         mat( 163) = 0._r8
         mat( 165) = 0._r8
         mat( 179) = 0._r8
         mat( 191) = 0._r8
         mat( 192) = 0._r8
         mat( 196) = 0._r8
         mat( 210) = 0._r8
         mat( 212) = 0._r8
         mat( 218) = 0._r8
         mat( 222) = 0._r8
         mat( 223) = 0._r8
         mat( 225) = 0._r8
         mat( 235) = 0._r8
         mat( 239) = 0._r8
         mat( 241) = 0._r8
         mat( 247) = 0._r8
         mat( 251) = 0._r8
         mat( 252) = 0._r8
         mat( 254) = 0._r8
         mat( 263) = 0._r8
         mat( 265) = 0._r8
         mat( 266) = 0._r8
         mat( 269) = 0._r8
         mat( 273) = 0._r8
         mat( 275) = 0._r8
         mat( 276) = 0._r8
         mat( 277) = 0._r8
         mat( 278) = 0._r8
         mat( 288) = 0._r8
         mat( 290) = 0._r8
         mat( 292) = 0._r8
         mat( 295) = 0._r8
         mat( 299) = 0._r8
         mat( 303) = 0._r8
         mat( 304) = 0._r8
         mat( 312) = 0._r8
         mat( 313) = 0._r8
         mat( 315) = 0._r8
         mat( 317) = 0._r8
         mat( 318) = 0._r8
         mat( 320) = 0._r8
         mat( 322) = 0._r8
         mat( 329) = 0._r8
         mat( 332) = 0._r8
         mat( 339) = 0._r8
         mat( 341) = 0._r8
         mat( 348) = 0._r8
         mat( 353) = 0._r8
         mat( 354) = 0._r8
         mat( 357) = 0._r8
         mat( 363) = 0._r8
         mat( 364) = 0._r8
         mat( 366) = 0._r8
         mat( 367) = 0._r8
         mat( 369) = 0._r8
         mat( 378) = 0._r8
         mat( 391) = 0._r8
         mat( 399) = 0._r8
         mat( 401) = 0._r8
         mat( 402) = 0._r8
         mat( 408) = 0._r8
         mat( 409) = 0._r8
         mat( 410) = 0._r8
         mat( 411) = 0._r8
         mat( 412) = 0._r8
         mat( 414) = 0._r8
         mat( 421) = 0._r8
         mat( 442) = 0._r8
         mat( 448) = 0._r8
         mat( 450) = 0._r8
         mat( 452) = 0._r8
         mat( 456) = 0._r8
         mat( 457) = 0._r8
         mat( 459) = 0._r8
         mat( 460) = 0._r8
         mat( 461) = 0._r8
         mat( 462) = 0._r8
         mat( 463) = 0._r8
         mat( 472) = 0._r8
         mat( 473) = 0._r8
         mat( 474) = 0._r8
         mat( 476) = 0._r8
         mat( 477) = 0._r8
         mat( 478) = 0._r8
         mat( 479) = 0._r8
         mat( 481) = 0._r8
         mat( 483) = 0._r8
         mat( 485) = 0._r8
         mat( 499) = 0._r8
         mat( 500) = 0._r8
         mat( 511) = 0._r8
         mat( 518) = 0._r8
         mat( 520) = 0._r8
         mat( 532) = 0._r8
         mat( 534) = 0._r8
         mat( 536) = 0._r8
         mat( 537) = 0._r8
         mat( 538) = 0._r8
         mat( 539) = 0._r8
         mat( 541) = 0._r8
         mat( 542) = 0._r8
         mat( 543) = 0._r8
         mat( 544) = 0._r8
         mat( 546) = 0._r8
         mat( 548) = 0._r8
         mat( 549) = 0._r8
         mat( 560) = 0._r8
         mat( 567) = 0._r8
         mat( 569) = 0._r8
         mat( 580) = 0._r8
         mat( 583) = 0._r8
         mat( 592) = 0._r8
         mat( 594) = 0._r8
         mat( 599) = 0._r8
         mat( 600) = 0._r8
         mat( 607) = 0._r8
         mat( 609) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 4) = mat( 4) - dti
         mat( 7) = mat( 7) - dti
         mat( 9) = mat( 9) - dti
         mat( 12) = mat( 12) - dti
         mat( 15) = mat( 15) - dti
         mat( 17) = mat( 17) - dti
         mat( 21) = mat( 21) - dti
         mat( 24) = mat( 24) - dti
         mat( 27) = mat( 27) - dti
         mat( 33) = mat( 33) - dti
         mat( 39) = mat( 39) - dti
         mat( 46) = mat( 46) - dti
         mat( 53) = mat( 53) - dti
         mat( 58) = mat( 58) - dti
         mat( 62) = mat( 62) - dti
         mat( 70) = mat( 70) - dti
         mat( 78) = mat( 78) - dti
         mat( 89) = mat( 89) - dti
         mat( 99) = mat( 99) - dti
         mat( 107) = mat( 107) - dti
         mat( 116) = mat( 116) - dti
         mat( 128) = mat( 128) - dti
         mat( 138) = mat( 138) - dti
         mat( 146) = mat( 146) - dti
         mat( 158) = mat( 158) - dti
         mat( 168) = mat( 168) - dti
         mat( 177) = mat( 177) - dti
         mat( 189) = mat( 189) - dti
         mat( 207) = mat( 207) - dti
         mat( 219) = mat( 219) - dti
         mat( 232) = mat( 232) - dti
         mat( 249) = mat( 249) - dti
         mat( 267) = mat( 267) - dti
         mat( 294) = mat( 294) - dti
         mat( 311) = mat( 311) - dti
         mat( 333) = mat( 333) - dti
         mat( 359) = mat( 359) - dti
         mat( 385) = mat( 385) - dti
         mat( 406) = mat( 406) - dti
         mat( 438) = mat( 438) - dti
         mat( 458) = mat( 458) - dti
         mat( 482) = mat( 482) - dti
         mat( 501) = mat( 501) - dti
         mat( 525) = mat( 525) - dti
         mat( 547) = mat( 547) - dti
         mat( 585) = mat( 585) - dti
         mat( 610) = mat( 610) - dti

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
