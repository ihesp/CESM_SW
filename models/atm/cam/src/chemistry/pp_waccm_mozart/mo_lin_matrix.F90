




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

         mat(501) = -( rxt(3) + rxt(4) + het_rates(1) )

         mat(585) = -( rxt(63) + rxt(64) + rxt(65) + rxt(76) + rxt(77) + rxt(78) &
                 + het_rates(2) )
         mat(215) = rxt(1) + 2.000_r8*rxt(2) + rxt(69) + rxt(70) + rxt(71) &
                      + 2.000_r8*rxt(74) + rxt(81) + rxt(82) + rxt(83) + 2.000_r8*rxt(86)
         mat(504) = rxt(4)
         mat(343) = rxt(6)
         mat(368) = rxt(8)
         mat(32) = rxt(10)
         mat(413) = rxt(12)
         mat(227) = rxt(21)
         mat(306) = rxt(24)
         mat(11) = rxt(25)
         mat(280) = rxt(32)
         mat(527) = rxt(102)

         mat(525) = -( rxt(102) + rxt(106)*y(7) + rxt(107)*y(7) + rxt(109)*y(43) &
                      + rxt(110)*y(44) + rxt(111)*y(45) + rxt(112)*y(53) + rxt(113)*y(54) &
                      + rxt(114)*y(46) + rxt(115)*y(51) + rxt(116)*y(52) + rxt(117)*y(47) &
                      + rxt(118)*y(42) + rxt(119)*y(50) + rxt(120)*y(49) + rxt(121)*y(55) &
                      + rxt(122)*y(56) + rxt(123)*y(57) + rxt(124)*y(58) + rxt(127)*y(15) &
                      + rxt(128)*y(15) + rxt(129)*y(15) + het_rates(3) )
         mat(214) = rxt(1)
         mat(502) = rxt(3)
         mat(226) = rxt(20)

         mat(207) = -( rxt(1) + rxt(2) + rxt(67) + rxt(69) + rxt(70) + rxt(71) + rxt(74) &
                      + rxt(79) + rxt(81) + rxt(82) + rxt(83) + rxt(86) + het_rates(4) )
         mat(490) = rxt(4)
         mat(398) = rxt(13)
         mat(8) = rxt(97)
         mat(5) = rxt(100) + rxt(101)
         mat(512) = rxt(107)*y(7)

         mat(7) = -( rxt(94) + rxt(97) + rxt(96)*y(63) + het_rates(5) )

         mat(4) = -( rxt(100) + rxt(101) + het_rates(6) )
         mat(488) = rxt(3)
         mat(6) = rxt(94) + rxt(96)*y(63)

         mat(249) = -( het_rates(20) )
         mat(233) = rxt(18)
         mat(220) = rxt(20)
         mat(515) = rxt(129)*y(15)

         mat(58) = -( het_rates(19) )
         mat(228) = rxt(17) + rxt(18)
         mat(588) = rxt(222)*y(41)
         mat(88) = rxt(268)*y(63)

         mat(158) = -( rxt(66) + het_rates(8) )
         mat(326) = rxt(6)
         mat(92) = rxt(265)

         mat(333) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(358) = rxt(8) + .500_r8*rxt(245)
         mat(29) = rxt(10)
         mat(403) = rxt(13)
         mat(112) = rxt(275)
         mat(517) = 2.000_r8*rxt(106)*y(7)

         mat(359) = -( rxt(8) + rxt(245) + het_rates(10) )
         mat(30) = rxt(9) + rxt(165)
         mat(193) = rxt(11)
         mat(404) = rxt(12)
         mat(42) = rxt(15) + rxt(174)
         mat(181) = rxt(30)
         mat(66) = rxt(36)

         mat(438) = -( rxt(223)*y(41) + rxt(224)*y(48) + rxt(225)*y(46) + rxt(226)*y(42) &
                      + rxt(228)*y(51) + rxt(229)*y(52) + rxt(230)*y(58) + rxt(231)*y(57) &
                      + rxt(234)*y(15) + het_rates(22) )
         mat(195) = rxt(11)
         mat(45) = rxt(14)
         mat(37) = rxt(16)
         mat(224) = rxt(19)
         mat(73) = 2.000_r8*rxt(22)
         mat(171) = rxt(27)
         mat(102) = rxt(33)
         mat(362) = .500_r8*rxt(245)
         mat(521) = rxt(127)*y(15)

         mat(406) = -( rxt(12) + rxt(13) + rxt(244) + het_rates(11) )
         mat(31) = rxt(9) + rxt(10) + rxt(165)
         mat(44) = rxt(14)
         mat(182) = rxt(29)
         mat(67) = rxt(35)

         mat(189) = -( rxt(11) + het_rates(12) )
         mat(28) = 2.000_r8*rxt(243) + 2.000_r8*rxt(247) + 2.000_r8*rxt(253) &
                      + 2.000_r8*rxt(258)
         mat(397) = rxt(244)
         mat(351) = .500_r8*rxt(245)
         mat(178) = rxt(248) + rxt(254) + rxt(259)
         mat(64) = rxt(249) + rxt(257) + rxt(260)

         mat(39) = -( rxt(14) + rxt(15) + rxt(174) + het_rates(13) )

         mat(27) = -( rxt(9) + rxt(10) + rxt(165) + rxt(243) + rxt(247) + rxt(253) &
                      + rxt(258) + het_rates(14) )

         mat(458) = -( het_rates(16) )
         mat(522) = rxt(127)*y(15)
         mat(604) = rxt(181)*y(15)
         mat(145) = rxt(220)*y(15)
         mat(439) = rxt(234)*y(15)

         mat(33) = -( rxt(16) + het_rates(17) )

         mat(232) = -( rxt(17) + rxt(18) + het_rates(18) )
         mat(35) = rxt(16)
         mat(514) = rxt(128)*y(15) + rxt(129)*y(15)

         mat(311) = -( het_rates(21) )
         mat(36) = rxt(16)
         mat(234) = 2.000_r8*rxt(17)
         mat(221) = rxt(19) + 2.000_r8*rxt(21)
         mat(475) = rxt(28)
         mat(148) = rxt(34)
         mat(26) = rxt(57)
         mat(516) = rxt(128)*y(15)

         mat(385) = -( rxt(246) + het_rates(23) )
         mat(43) = rxt(15) + rxt(174)
         mat(519) = rxt(128)*y(15)
         mat(601) = rxt(222)*y(41) + rxt(227)*y(42)
         mat(436) = rxt(223)*y(41) + rxt(226)*y(42)

         mat(70) = -( rxt(22) + het_rates(24) )
         mat(372) = .500_r8*rxt(246)

         mat(219) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(71) )
         mat(428) = rxt(223)*y(41) + rxt(224)*y(48) + rxt(225)*y(46) + rxt(226)*y(42) &
                      + rxt(230)*y(58) + rxt(234)*y(15)

         mat(610) = -( rxt(181)*y(15) + rxt(222)*y(41) + rxt(227)*y(42) + rxt(232)*y(58) &
                      + rxt(233)*y(57) + het_rates(27) )
         mat(16) = 2.000_r8*rxt(23)
         mat(307) = rxt(24)
         mat(3) = 2.000_r8*rxt(26)
         mat(174) = rxt(27)
         mat(487) = rxt(28)
         mat(186) = rxt(29)
         mat(23) = rxt(31)
         mat(20) = rxt(56)
         mat(528) = 2.000_r8*rxt(109)*y(43) + 2.000_r8*rxt(110)*y(44) &
                      + 2.000_r8*rxt(111)*y(45) + 2.000_r8*rxt(112)*y(53) + rxt(113)*y(54) &
                      + rxt(114)*y(46) + rxt(115)*y(51) + rxt(116)*y(52) &
                      + 4.000_r8*rxt(117)*y(47) + rxt(119)*y(50)
         mat(445) = rxt(223)*y(41) + 3.000_r8*rxt(224)*y(48) + rxt(225)*y(46) &
                      + rxt(228)*y(51) + rxt(229)*y(52)

         mat(15) = -( rxt(23) + het_rates(28) )

         mat(294) = -( rxt(24) + het_rates(29) )
         mat(10) = rxt(25)
         mat(180) = rxt(30)
         mat(2) = 2.000_r8*rxt(193)

         mat(9) = -( rxt(25) + het_rates(30) )

         mat(1) = -( rxt(26) + rxt(193) + het_rates(31) )

         mat(482) = -( rxt(28) + het_rates(32) )
         mat(605) = rxt(181)*y(15) + 2.000_r8*rxt(222)*y(41) + rxt(227)*y(42) &
                      + rxt(232)*y(58) + rxt(233)*y(57)

         mat(168) = -( rxt(27) + het_rates(33) )
         mat(176) = rxt(248) + rxt(254) + rxt(259)

         mat(177) = -( rxt(29) + rxt(30) + rxt(248) + rxt(254) + rxt(259) + het_rates(34) &
       )

         mat(21) = -( rxt(31) + het_rates(35) )

         mat(547) = -( het_rates(36) )
         mat(22) = rxt(31)
         mat(279) = rxt(32)
         mat(104) = rxt(33)
         mat(151) = rxt(34)
         mat(68) = rxt(35)
         mat(526) = rxt(118)*y(42) + rxt(119)*y(50) + rxt(120)*y(49) &
                      + 2.000_r8*rxt(121)*y(55) + 2.000_r8*rxt(122)*y(56) &
                      + 3.000_r8*rxt(123)*y(57) + 2.000_r8*rxt(124)*y(58)
         mat(443) = rxt(226)*y(42) + 2.000_r8*rxt(230)*y(58) + 3.000_r8*rxt(231)*y(57)
         mat(608) = rxt(227)*y(42) + 2.000_r8*rxt(232)*y(58) + 3.000_r8*rxt(233)*y(57)

         mat(267) = -( rxt(32) + het_rates(37) )
         mat(65) = rxt(36)

         mat(146) = -( rxt(34) + het_rates(38) )

         mat(99) = -( rxt(33) + het_rates(39) )
         mat(63) = rxt(249) + rxt(257) + rxt(260)

         mat(62) = -( rxt(35) + rxt(36) + rxt(249) + rxt(257) + rxt(260) + het_rates(40) &
       )

         mat(78) = -( het_rates(64) )

         mat(107) = -( rxt(275) + het_rates(65) )
         mat(203) = rxt(67) + rxt(79)
         mat(90) = rxt(268)*y(63)

         mat(46) = -( het_rates(66) )
         mat(153) = rxt(66)

         mat(89) = -( rxt(265) + rxt(268)*y(63) + het_rates(67) )
         mat(558) = rxt(63) + rxt(64) + rxt(65) + rxt(76) + rxt(77) + rxt(78)
         mat(202) = rxt(69) + rxt(70) + rxt(71) + rxt(81) + rxt(82) + rxt(83)

         mat(116) = -( het_rates(68) )
         mat(324) = rxt(7)
         mat(91) = rxt(265)
         mat(108) = rxt(275)

         mat(53) = -( het_rates(70) )

         mat(128) = -( het_rates(69) )
         mat(325) = rxt(7)
         mat(562) = rxt(63) + rxt(64) + rxt(65) + rxt(76) + rxt(77) + rxt(78)
         mat(157) = rxt(66)
         mat(205) = rxt(67) + rxt(69) + rxt(70) + rxt(71) + rxt(79) + rxt(81) + rxt(82) &
                      + rxt(83)

         mat(12) = -( rxt(55) + het_rates(59) )
         mat(507) = rxt(110)*y(44) + rxt(111)*y(45) + 2.000_r8*rxt(112)*y(53) &
                      + 2.000_r8*rxt(113)*y(54) + rxt(114)*y(46) + rxt(116)*y(52) &
                      + rxt(119)*y(50) + rxt(120)*y(49) + rxt(121)*y(55) &
                      + 2.000_r8*rxt(122)*y(56)
         mat(415) = rxt(225)*y(46) + rxt(229)*y(52)

         mat(17) = -( rxt(56) + het_rates(60) )
         mat(508) = rxt(109)*y(43) + rxt(111)*y(45) + rxt(115)*y(51)
         mat(416) = rxt(228)*y(51)

         mat(24) = -( rxt(57) + het_rates(61) )
         mat(137) = rxt(220)*y(15)

         mat(138) = -( rxt(220)*y(15) + het_rates(62) )
         mat(13) = 2.000_r8*rxt(55)
         mat(18) = rxt(56)
         mat(25) = rxt(57)
         mat(509) = rxt(113)*y(54) + rxt(120)*y(49)


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
