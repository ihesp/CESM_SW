




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

         mat(649) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(536) = rxt(45)

         mat(533) = -( rxt(45) + het_rates(2) )
         mat(643) = rxt(3)
         mat(668) = rxt(5)
         mat(114) = rxt(7)
         mat(965) = rxt(9)
         mat(439) = rxt(47) + rxt(48)

         mat(438) = -( rxt(47) + rxt(48) + rxt(50) + rxt(52)*y(4) + rxt(53)*y(4) &
                      + rxt(54)*y(16) + rxt(55)*y(16) + rxt(56)*y(16) + het_rates(3) )
         mat(635) = rxt(2)

         mat(617) = -( het_rates(17) )
         mat(605) = rxt(14) + rxt(15)
         mat(448) = rxt(17)
         mat(499) = 1.340_r8*rxt(23)
         mat(550) = .700_r8*rxt(24)
         mat(507) = rxt(30)
         mat(459) = rxt(32)
         mat(370) = rxt(35)
         mat(276) = .450_r8*rxt(37)
         mat(303) = 2.000_r8*rxt(38)

         mat(279) = -( het_rates(11) )
         mat(600) = rxt(15)
         mat(436) = rxt(56)*y(16)

         mat(87) = -( het_rates(116) )

         mat(54) = -( het_rates(117) )

         mat(401) = -( rxt(58) + het_rates(15) )
         mat(144) = rxt(13)
         mat(437) = rxt(55)*y(16)

         mat(943) = -( het_rates(5) )
         mat(678) = rxt(5) + .500_r8*rxt(231)
         mat(116) = rxt(7)
         mat(978) = rxt(10)
         mat(445) = 2.000_r8*rxt(52)*y(4)

         mat(672) = -( rxt(5) + rxt(231) + het_rates(6) )
         mat(115) = rxt(6) + rxt(84)
         mat(234) = rxt(8)
         mat(972) = rxt(9)
         mat(134) = rxt(12) + rxt(93)
         mat(221) = .600_r8*rxt(20) + rxt(136)
         mat(247) = rxt(21) + rxt(182)
         mat(460) = rxt(32)

         mat(979) = -( rxt(9) + rxt(10) + rxt(230) + het_rates(7) )
         mat(117) = rxt(6) + rxt(7) + rxt(84)
         mat(137) = rxt(11)
         mat(225) = .400_r8*rxt(20)

         mat(233) = -( rxt(8) + het_rates(8) )
         mat(113) = 2.000_r8*rxt(229)
         mat(950) = rxt(230)
         mat(662) = .500_r8*rxt(231)

         mat(133) = -( rxt(11) + rxt(12) + rxt(93) + het_rates(9) )

         mat(112) = -( rxt(6) + rxt(7) + rxt(84) + rxt(229) + het_rates(10) )

         mat(759) = -( rxt(94)*y(16) + het_rates(12) )
         mat(235) = rxt(8)
         mat(135) = rxt(11)
         mat(146) = rxt(13)
         mat(102) = 2.000_r8*rxt(16)
         mat(231) = rxt(18)
         mat(208) = rxt(19)
         mat(162) = rxt(25)
         mat(76) = rxt(26)
         mat(126) = rxt(27)
         mat(131) = rxt(28)
         mat(93) = rxt(31)
         mat(314) = rxt(39)
         mat(121) = rxt(40)
         mat(179) = rxt(41)
         mat(242) = rxt(42)
         mat(442) = 2.000_r8*rxt(50) + rxt(54)*y(16)
         mat(673) = .500_r8*rxt(231)

         mat(898) = -( rxt(239) + het_rates(13) )
         mat(136) = rxt(12) + rxt(93)
         mat(453) = rxt(17)
         mat(232) = rxt(18)
         mat(504) = 1.340_r8*rxt(22) + .660_r8*rxt(23)
         mat(163) = rxt(25)
         mat(127) = rxt(27)
         mat(511) = rxt(30)
         mat(462) = rxt(32)
         mat(259) = rxt(33)
         mat(494) = rxt(34)
         mat(372) = 2.000_r8*rxt(35)
         mat(278) = .560_r8*rxt(37)
         mat(305) = 2.000_r8*rxt(38)
         mat(316) = .900_r8*rxt(39)
         mat(243) = rxt(42)
         mat(406) = rxt(58)
         mat(196) = rxt(108)
         mat(111) = rxt(116) + rxt(117)
         mat(444) = rxt(55)*y(16)

         mat(100) = -( rxt(16) + het_rates(14) )
         mat(843) = .500_r8*rxt(239)

         mat(793) = -( het_rates(18) )
         mat(451) = rxt(17)
         mat(209) = rxt(19)
         mat(223) = .400_r8*rxt(20)
         mat(554) = .300_r8*rxt(24)
         mat(328) = rxt(29)
         mat(443) = rxt(54)*y(16)
         mat(760) = rxt(94)*y(16)

         mat(143) = -( rxt(13) + het_rates(19) )

         mat(604) = -( rxt(14) + rxt(15) + het_rates(20) )
         mat(145) = rxt(13)
         mat(230) = rxt(18)
         mat(498) = 1.340_r8*rxt(22)
         mat(130) = rxt(28)
         mat(458) = rxt(32)
         mat(257) = .690_r8*rxt(33)
         mat(491) = rxt(34)
         mat(369) = rxt(35)
         mat(313) = .100_r8*rxt(39)
         mat(193) = rxt(108)
         mat(110) = 2.000_r8*rxt(117)
         mat(440) = rxt(55)*y(16) + rxt(56)*y(16)

         mat(260) = -( het_rates(21) )

         mat(104) = -( het_rates(22) )

         mat(74) = -( rxt(26) + het_rates(28) )

         mat(152) = -( het_rates(23) )

         mat(108) = -( rxt(116) + rxt(117) + het_rates(24) )
         mat(75) = rxt(26)

         mat(295) = -( het_rates(25) )

         mat(203) = -( het_rates(26) )

         mat(368) = -( rxt(35) + het_rates(27) )
         mat(109) = rxt(116)

         mat(57) = -( het_rates(29) )

         mat(359) = -( het_rates(30) )
         mat(200) = rxt(36)

         mat(159) = -( rxt(25) + het_rates(31) )

         mat(447) = -( rxt(17) + het_rates(32) )
         mat(228) = rxt(18)
         mat(161) = rxt(25)
         mat(312) = .400_r8*rxt(39)
         mat(120) = rxt(40)

         mat(818) = -( het_rates(33) )
         mat(224) = .600_r8*rxt(20) + rxt(136)
         mat(502) = 1.340_r8*rxt(22)
         mat(555) = .300_r8*rxt(24)
         mat(132) = rxt(28)
         mat(329) = rxt(29)
         mat(510) = rxt(30)
         mat(493) = rxt(34)
         mat(202) = rxt(36)
         mat(277) = .130_r8*rxt(37)
         mat(122) = rxt(40)

         mat(206) = -( rxt(19) + het_rates(34) )

         mat(422) = -( het_rates(35) )
         mat(543) = .700_r8*rxt(24)

         mat(60) = -( het_rates(36) )

         mat(375) = -( het_rates(37) )

         mat(123) = -( rxt(27) + het_rates(38) )

         mat(331) = -( het_rates(39) )

         mat(226) = -( rxt(18) + het_rates(40) )

         mat(325) = -( rxt(29) + het_rates(41) )
         mat(124) = .820_r8*rxt(27)
         mat(309) = .250_r8*rxt(39)
         mat(238) = .100_r8*rxt(42)

         mat(479) = -( het_rates(42) )

         mat(128) = -( rxt(28) + het_rates(43) )

         mat(42) = -( het_rates(44) )

         mat(164) = -( het_rates(45) )

         mat(45) = -( het_rates(49) )

         mat(344) = -( het_rates(50) )

         mat(307) = -( rxt(39) + het_rates(51) )

         mat(198) = -( rxt(36) + het_rates(46) )
         mat(306) = .800_r8*rxt(39)

         mat(318) = -( het_rates(47) )

         mat(118) = -( rxt(40) + het_rates(48) )

         mat(386) = -( het_rates(52) )

         mat(586) = -( het_rates(53) )

         mat(252) = -( rxt(33) + het_rates(54) )

         mat(547) = -( rxt(24) + het_rates(55) )
         mat(255) = .402_r8*rxt(33)
         mat(241) = rxt(42)

         mat(495) = -( rxt(22) + rxt(23) + het_rates(56) )
         mat(253) = .288_r8*rxt(33)
         mat(240) = rxt(42)

         mat(566) = -( het_rates(57) )

         mat(138) = -( het_rates(58) )

         mat(835) = -( het_rates(59) )
         mat(249) = rxt(21) + rxt(182)
         mat(503) = .660_r8*rxt(22)

         mat(171) = -( het_rates(60) )

         mat(489) = -( rxt(34) + het_rates(61) )

         mat(506) = -( rxt(30) + het_rates(62) )
         mat(275) = .180_r8*rxt(37)
         mat(178) = .450_r8*rxt(41)

         mat(519) = -( het_rates(63) )

         mat(91) = -( rxt(31) + het_rates(64) )

         mat(284) = -( het_rates(65) )

         mat(409) = -( het_rates(66) )

         mat(237) = -( rxt(42) + het_rates(67) )

         mat(66) = -( het_rates(68) )

         mat(71) = -( het_rates(69) )

         mat(265) = -( het_rates(70) )

         mat(174) = -( rxt(41) + het_rates(71) )

         mat(83) = -( het_rates(78) )

         mat(273) = -( rxt(37) + het_rates(79) )
         mat(176) = .900_r8*rxt(41)

         mat(302) = -( rxt(38) + het_rates(80) )
         mat(274) = .130_r8*rxt(37)
         mat(177) = .450_r8*rxt(41)


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

         mat(48) = -( het_rates(72) )

         mat(181) = -( het_rates(73) )

         mat(1) = -( het_rates(74) )

         mat(51) = -( het_rates(75) )

         mat(212) = -( het_rates(76) )

         mat(2) = -( het_rates(77) )

         mat(219) = -( rxt(20) + rxt(136) + het_rates(81) )

         mat(187) = -( het_rates(82) )

         mat(244) = -( rxt(21) + rxt(182) + het_rates(83) )

         mat(466) = -( het_rates(84) )

         mat(456) = -( rxt(32) + het_rates(85) )

         mat(64) = -( het_rates(100) )

         mat(95) = -( het_rates(101) )

         mat(3) = -( het_rates(102) )

         mat(40) = -( het_rates(103) )

         mat(4) = -( het_rates(104) )

         mat(5) = -( het_rates(105) )

         mat(6) = -( het_rates(90) )

         mat(7) = -( het_rates(91) )

         mat(8) = -( het_rates(92) )

         mat(9) = -( het_rates(93) )

         mat(10) = -( het_rates(94) )

         mat(11) = -( het_rates(95) )

         mat(12) = -( het_rates(96) )

         mat(13) = -( het_rates(97) )

         mat(14) = -( het_rates(98) )

         mat(15) = -( het_rates(99) )

         mat(16) = -( rxt(232) + het_rates(86) )

         mat(18) = -( het_rates(87) )
         mat(17) = rxt(232)

         mat(19) = -( rxt(238) + het_rates(88) )

         mat(21) = -( het_rates(89) )
         mat(20) = rxt(238)

         mat(77) = -( het_rates(118) )

         mat(148) = -( het_rates(119) )

         mat(192) = -( rxt(108) + het_rates(120) )

         mat(22) = -( het_rates(106) )

         mat(23) = -( het_rates(107) )

         mat(24) = -( het_rates(108) )

         mat(25) = -( het_rates(109) )

         mat(26) = -( het_rates(110) )

         mat(27) = -( het_rates(111) )

         mat(28) = -( het_rates(112) )

         mat(29) = -( het_rates(113) )

         mat(30) = -( rxt(240) + het_rates(121) )

         mat(31) = -( rxt(241) + het_rates(122) )

         mat(32) = -( rxt(242) + het_rates(123) )

         mat(33) = -( het_rates(124) )

         mat(34) = -( rxt(243) + het_rates(125) )

         mat(35) = -( rxt(244) + het_rates(126) )

         mat(36) = -( rxt(245) + het_rates(127) )

         mat(37) = -( het_rates(128) )

         mat(38) = -( rxt(43) + het_rates(129) )

         mat(39) = -( rxt(44) + het_rates(130) )


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
