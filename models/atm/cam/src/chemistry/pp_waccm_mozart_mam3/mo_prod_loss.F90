



      module mo_prod_loss

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: exp_prod_loss
      public :: imp_prod_loss

      contains

      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:,:,:)
      real(r8), intent(in) :: rxt(:,:,:)
      real(r8), intent(in) :: het_rates(:,:,:)



!--------------------------------------------------------------------
! ... loss and production for Explicit method
!--------------------------------------------------------------------


         loss(:,:,1) = ((rxt(:,:,108) +rxt(:,:,109) +rxt(:,:,110))* y(:,:,3) &
                  +rxt(:,:,136)* y(:,:,22) +rxt(:,:,167)* y(:,:,27) + rxt(:,:,47) &
                  + rxt(:,:,48) + het_rates(:,:,15))* y(:,:,15)
         prod(:,:,1) = 0._r8
         loss(:,:,2) = ((rxt(:,:,97) +rxt(:,:,98))* y(:,:,3) + rxt(:,:,5) &
                  + het_rates(:,:,7))* y(:,:,7)
         prod(:,:,2) = 0._r8
         loss(:,:,3) = ((rxt(:,:,143) +rxt(:,:,144))* y(:,:,22) + het_rates(:,:,19)) &
                 * y(:,:,19)
         prod(:,:,3) = (rxt(:,:,46) +rxt(:,:,79) +rxt(:,:,249)*y(:,:,55))*y(:,:,51) &
                  +.380_r8*rxt(:,:,48)*y(:,:,15) +rxt(:,:,207)*y(:,:,41)*y(:,:,27)
         loss(:,:,4) = (rxt(:,:,111)* y(:,:,3) +rxt(:,:,155)* y(:,:,22) +rxt(:,:,162) &
                 * y(:,:,27) + het_rates(:,:,20))* y(:,:,20)
         prod(:,:,4) = (1.440_r8*rxt(:,:,48) +rxt(:,:,110)*y(:,:,3))*y(:,:,15)
         loss(:,:,5) = (rxt(:,:,208)* y(:,:,22) +rxt(:,:,207)* y(:,:,27) + rxt(:,:,36) &
                  + het_rates(:,:,41))* y(:,:,41)
         prod(:,:,5) = 0._r8
         loss(:,:,6) = (rxt(:,:,105)* y(:,:,3) +rxt(:,:,211)* y(:,:,22) + rxt(:,:,43) &
                  + het_rates(:,:,42))* y(:,:,42)
         prod(:,:,6) = 0._r8
         loss(:,:,7) = (rxt(:,:,100)* y(:,:,3) + rxt(:,:,39) + het_rates(:,:,43)) &
                 * y(:,:,43)
         prod(:,:,7) = 0._r8
         loss(:,:,8) = (rxt(:,:,101)* y(:,:,3) + rxt(:,:,40) + het_rates(:,:,44)) &
                 * y(:,:,44)
         prod(:,:,8) = 0._r8
         loss(:,:,9) = (rxt(:,:,102)* y(:,:,3) + rxt(:,:,41) + het_rates(:,:,45)) &
                 * y(:,:,45)
         prod(:,:,9) = 0._r8
         loss(:,:,10) = (rxt(:,:,103)* y(:,:,3) +rxt(:,:,210)* y(:,:,22) + rxt(:,:,42) &
                  + het_rates(:,:,46))* y(:,:,46)
         prod(:,:,10) = 0._r8
         loss(:,:,11) = (rxt(:,:,104)* y(:,:,3) + rxt(:,:,37) + het_rates(:,:,47)) &
                 * y(:,:,47)
         prod(:,:,11) = 0._r8
         loss(:,:,12) = (rxt(:,:,209)* y(:,:,22) + rxt(:,:,38) + het_rates(:,:,48)) &
                 * y(:,:,48)
         prod(:,:,12) = 0._r8
         loss(:,:,13) = (rxt(:,:,107)* y(:,:,3) + rxt(:,:,44) + het_rates(:,:,49)) &
                 * y(:,:,49)
         prod(:,:,13) = 0._r8
         loss(:,:,14) = (rxt(:,:,106)* y(:,:,3) + rxt(:,:,45) + het_rates(:,:,50)) &
                 * y(:,:,50)
         prod(:,:,14) = 0._r8
         loss(:,:,15) = (rxt(:,:,249)* y(:,:,55) + rxt(:,:,46) + rxt(:,:,79) &
                  + het_rates(:,:,51))* y(:,:,51)
         prod(:,:,15) = (rxt(:,:,143)*y(:,:,22) +rxt(:,:,144)*y(:,:,22))*y(:,:,19) &
                  +.440_r8*rxt(:,:,48)*y(:,:,15)
         loss(:,:,16) = ( + het_rates(:,:,25))* y(:,:,25)
         prod(:,:,16) = 0._r8
         loss(:,:,17) = ( + het_rates(:,:,26))* y(:,:,26)
         prod(:,:,17) = 0._r8

      end subroutine exp_prod_loss

      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:)
      real(r8), intent(in) :: rxt(:)
      real(r8), intent(in) :: het_rates(:)



!--------------------------------------------------------------------
! ... loss and production for Implicit method
!--------------------------------------------------------------------


         loss(61) = (rxt(81)* y(2) +rxt(99)* y(3) +rxt(121)* y(9) +rxt(124)* y(10) &
                  +rxt(146)* y(21) +rxt(151)* y(22) +rxt(158)* y(23) +rxt(161)* y(27) &
                  +rxt(187)* y(36) +rxt(216)* y(60) +rxt(219)* y(61) + rxt(3) + rxt(4) &
                  + het_rates(1))* y(1)
         prod(61) =rxt(80)*y(4)*y(2)
         loss(56) = (rxt(81)* y(1) + 2._r8*rxt(82)* y(2) +rxt(80)* y(4) +rxt(119) &
                 * y(9) + (rxt(122) +rxt(123))* y(10) +rxt(130)* y(11) +rxt(142) &
                 * y(18) +rxt(150)* y(22) +rxt(157)* y(23) +rxt(160)* y(24) +rxt(168) &
                 * y(29) +rxt(180)* y(32) +rxt(181)* y(33) +rxt(184)* y(34) +rxt(190) &
                 * y(37) +rxt(200)* y(38) +rxt(201)* y(39) +rxt(202)* y(40) &
                  + (rxt(248) +rxt(257))* y(52) +rxt(254)* y(54) +rxt(212)* y(59) &
                  + rxt(54) + rxt(55) + rxt(56) + rxt(67) + rxt(68) + rxt(69) &
                  + het_rates(2))* y(2)
         prod(56) = (rxt(1) +2.000_r8*rxt(2) +rxt(60) +rxt(61) +rxt(62) + &
                 2.000_r8*rxt(65) +rxt(72) +rxt(73) +rxt(74) +2.000_r8*rxt(77) + &
                 rxt(94)*y(3) +rxt(95)*y(3) +rxt(116)*y(8) +rxt(215)*y(60) + &
                 rxt(218)*y(61) +rxt(246)*y(55) +rxt(253)*y(54))*y(4) &
                  + (rxt(117)*y(9) +rxt(118)*y(10) +rxt(250)*y(53))*y(8) &
                  + (rxt(258)*y(56) +1.150_r8*rxt(259)*y(53))*y(57) +rxt(4)*y(1) &
                  +rxt(93)*y(3) +rxt(6)*y(9) +rxt(8)*y(10) +rxt(12)*y(11) +rxt(10) &
                 *y(14) +rxt(149)*y(23)*y(21) +rxt(153)*y(22)*y(22) +rxt(24)*y(29) &
                  +rxt(25)*y(30) +rxt(32)*y(37) +rxt(53)*y(61) +rxt(50)*y(62) +rxt(51) &
                 *y(63) +rxt(21)*y(82)
         loss(54) = (rxt(99)* y(1) + (rxt(94) +rxt(95))* y(4) + (rxt(97) +rxt(98)) &
                 * y(7) + (rxt(108) +rxt(109) +rxt(110))* y(15) +rxt(111)* y(20) &
                  +rxt(112)* y(32) +rxt(113)* y(38) +rxt(105)* y(42) +rxt(100)* y(43) &
                  +rxt(101)* y(44) +rxt(102)* y(45) +rxt(103)* y(46) +rxt(104)* y(47) &
                  +rxt(107)* y(49) +rxt(106)* y(50) +rxt(96)* y(82) + rxt(93) &
                  + het_rates(3))* y(3)
         prod(54) = (rxt(1) +rxt(114)*y(58))*y(4) +rxt(3)*y(1) +.850_r8*rxt(259)*y(57) &
                 *y(53) +rxt(20)*y(82)
         loss(51) = (rxt(80)* y(2) +rxt(94)* y(3) +rxt(90)* y(6) +rxt(116)* y(8) &
                  +rxt(145)* y(21) +rxt(255)* y(52) + (rxt(252) +rxt(253))* y(54) &
                  +rxt(246)* y(55) +rxt(114)* y(58) +rxt(215)* y(60) +rxt(218)* y(61) &
                  + rxt(1) + rxt(2) + rxt(58) + rxt(60) + rxt(61) + rxt(62) + rxt(65) &
                  + rxt(70) + rxt(72) + rxt(73) + rxt(74) + rxt(77) + het_rates(4)) &
                 * y(4)
         prod(51) = (rxt(4) +2.000_r8*rxt(81)*y(2) +2.000_r8*rxt(99)*y(3) + &
                 rxt(121)*y(9) +rxt(124)*y(10) +rxt(146)*y(21) +rxt(151)*y(22) + &
                 2.000_r8*rxt(158)*y(23) +rxt(161)*y(27) +rxt(187)*y(36) + &
                 rxt(216)*y(60) +rxt(219)*y(61))*y(1) + (rxt(132)*y(11) + &
                 rxt(138)*y(16) +rxt(148)*y(21) +rxt(152)*y(22) +rxt(157)*y(2) + &
                 rxt(159)*y(23) +rxt(164)*y(27) +rxt(171)*y(29) +rxt(188)*y(36) + &
                 rxt(192)*y(37))*y(23) + (rxt(82)*y(2) +rxt(89)*y(6) +rxt(122)*y(10) + &
                 rxt(130)*y(11) +rxt(150)*y(22) +rxt(168)*y(29) +rxt(190)*y(37))*y(2) &
                  + (rxt(170)*y(22) +rxt(174)*y(29) +rxt(175)*y(29) +rxt(196)*y(37) + &
                 rxt(197)*y(37))*y(29) + (rxt(91) +rxt(92) +2.000_r8*rxt(90)*y(4)) &
                 *y(6) +rxt(98)*y(7)*y(3) +rxt(88)*y(5) +rxt(251)*y(53)*y(9) +rxt(13) &
                 *y(11) +rxt(134)*y(22)*y(13) +rxt(198)*y(37)*y(37)
         loss(20) = (rxt(86)* y(1) +rxt(83)* y(2) +rxt(84)* y(4) +rxt(87)* y(51) &
                  + rxt(85) + rxt(88) + het_rates(5))* y(5)
         prod(20) =rxt(94)*y(4)*y(3)
         loss(19) = (rxt(89)* y(2) +rxt(90)* y(4) + rxt(91) + rxt(92) + het_rates(6)) &
                 * y(6)
         prod(19) = (rxt(85) +rxt(87)*y(51) +rxt(83)*y(2) +rxt(84)*y(4) +rxt(86)*y(1)) &
                 *y(5) +rxt(3)*y(1)
         loss(43) = (rxt(116)* y(4) +rxt(117)* y(9) +rxt(118)* y(10) +rxt(250)* y(53) &
                  + rxt(57) + het_rates(8))* y(8)
         prod(43) = (rxt(115)*y(58) +rxt(254)*y(54))*y(2) + (.200_r8*rxt(258)*y(56) + &
                 1.100_r8*rxt(260)*y(52))*y(57) +rxt(252)*y(54)*y(4) +rxt(6)*y(9) &
                  +rxt(247)*y(55)
         loss(57) = (rxt(121)* y(1) +rxt(119)* y(2) +rxt(117)* y(8) +rxt(129)* y(11) &
                  +rxt(137)* y(16) +rxt(120)* y(23) +rxt(172)* y(29) +rxt(193)* y(37) &
                  +rxt(251)* y(53) + rxt(6) + rxt(7) + het_rates(9))* y(9)
         prod(57) = (rxt(8) +.500_r8*rxt(205) +rxt(122)*y(2) +rxt(220)*y(61))*y(10) &
                  + (rxt(114)*y(58) +rxt(116)*y(8))*y(4) +2.000_r8*rxt(97)*y(7)*y(3) &
                  +rxt(13)*y(11) +rxt(10)*y(14) +rxt(256)*y(53)
         loss(58) = (rxt(124)* y(1) + (rxt(122) +rxt(123))* y(2) +rxt(118)* y(8) &
                  +rxt(125)* y(11) +rxt(127)* y(22) +rxt(133)* y(23) +rxt(173)* y(29) &
                  +rxt(194)* y(37) +rxt(220)* y(61) + rxt(8) + rxt(205) &
                  + het_rates(10))* y(10)
         prod(58) = (rxt(119)*y(2) +rxt(120)*y(23) +rxt(121)*y(1) + &
                 2.000_r8*rxt(129)*y(11) +rxt(137)*y(16) +rxt(172)*y(29) + &
                 rxt(193)*y(37))*y(9) + (rxt(12) +rxt(130)*y(2) +rxt(131)*y(22) + &
                 rxt(132)*y(23))*y(11) + (rxt(15) +rxt(135) +rxt(134)*y(22))*y(13) &
                  + (rxt(9) +rxt(126))*y(14) +rxt(11)*y(12) +rxt(30)*y(34) +rxt(35) &
                 *y(40)
         loss(59) = (rxt(151)* y(1) +rxt(150)* y(2) +rxt(127)* y(10) +rxt(131)* y(11) &
                  +rxt(128)* y(12) +rxt(134)* y(13) +rxt(136)* y(15) +rxt(139)* y(17) &
                  +rxt(141)* y(18) + (rxt(143) +rxt(144))* y(19) +rxt(155)* y(20) &
                  + 2._r8*(rxt(153) +rxt(154))* y(22) +rxt(152)* y(23) +rxt(156) &
                 * y(24) + (rxt(169) +rxt(170))* y(29) +rxt(179)* y(32) +rxt(183) &
                 * y(33) +rxt(185)* y(34) +rxt(191)* y(37) +rxt(199)* y(38) +rxt(208) &
                 * y(41) +rxt(211)* y(42) +rxt(210)* y(46) +rxt(209)* y(48) +rxt(213) &
                 * y(59) +rxt(214)* y(60) +rxt(217)* y(61) +rxt(224)* y(62) &
                  + (rxt(226) +rxt(227))* y(65) + het_rates(22))* y(22)
         prod(59) = (rxt(142)*y(18) +rxt(157)*y(23) +rxt(160)*y(24) +rxt(180)*y(32) + &
                 rxt(181)*y(33) +rxt(200)*y(38) +rxt(201)*y(39))*y(2) &
                  + (rxt(108)*y(15) +rxt(111)*y(20) +2.000_r8*rxt(96)*y(82) + &
                 rxt(112)*y(32) +rxt(113)*y(38))*y(3) + (rxt(120)*y(9) + &
                 rxt(132)*y(11) +2.000_r8*rxt(147)*y(21) +rxt(158)*y(1) + &
                 rxt(165)*y(27))*y(23) +rxt(146)*y(21)*y(1) +.500_r8*rxt(205)*y(10) &
                  +rxt(11)*y(12) +rxt(14)*y(13) +rxt(16)*y(17) +2.000_r8*rxt(22)*y(24) &
                  +rxt(27)*y(33) +rxt(33)*y(39) +rxt(19)*y(82)
         loss(60) = (rxt(130)* y(2) +rxt(129)* y(9) +rxt(125)* y(10) +rxt(140)* y(18) &
                  +rxt(131)* y(22) +rxt(132)* y(23) +rxt(228)* y(65) + rxt(12) &
                  + rxt(13) + rxt(204) + het_rates(11))* y(11)
         prod(60) = (rxt(29) +rxt(184)*y(2) +rxt(185)*y(22) +rxt(186)*y(27))*y(34) &
                  + (rxt(9) +rxt(10) +rxt(126))*y(14) + (rxt(123)*y(10) + &
                 rxt(202)*y(40))*y(2) +rxt(124)*y(10)*y(1) +rxt(128)*y(22)*y(12) &
                  +rxt(14)*y(13) +rxt(34)*y(40)
         loss(38) = (rxt(128)* y(22) + rxt(11) + het_rates(12))* y(12)
         prod(38) = (rxt(230) +rxt(236) +rxt(241) +rxt(232)*y(32) +rxt(237)*y(32) + &
                 rxt(243)*y(32))*y(34) + (2.000_r8*rxt(203) +2.000_r8*rxt(229) + &
                 2.000_r8*rxt(235) +2.000_r8*rxt(240))*y(14) + (rxt(204) + &
                 rxt(140)*y(18) +rxt(228)*y(65))*y(11) + (rxt(231) +rxt(239) + &
                 rxt(242))*y(40) + (.500_r8*rxt(205) +rxt(127)*y(22))*y(10)
         loss(30) = (rxt(134)* y(22) + rxt(14) + rxt(15) + rxt(135) + het_rates(13)) &
                 * y(13)
         prod(30) =rxt(133)*y(23)*y(10)
         loss(25) = ( + rxt(9) + rxt(10) + rxt(126) + rxt(203) + rxt(229) + rxt(235) &
                  + rxt(240) + het_rates(14))* y(14)
         prod(25) =rxt(125)*y(11)*y(10)
         loss(44) = (rxt(137)* y(9) +rxt(138)* y(23) + het_rates(16))* y(16)
         prod(44) = (rxt(108)*y(3) +rxt(136)*y(22) +rxt(167)*y(27))*y(15) &
                  +rxt(139)*y(22)*y(17)
         loss(28) = (rxt(139)* y(22) + rxt(16) + het_rates(17))* y(17)
         prod(28) =rxt(138)*y(23)*y(16)
         loss(49) = (rxt(142)* y(2) +rxt(140)* y(11) +rxt(141)* y(22) +rxt(166)* y(27) &
                  +rxt(189)* y(36) + rxt(17) + rxt(18) + het_rates(18))* y(18)
         prod(49) = (rxt(109)*y(15) +rxt(110)*y(15))*y(3) +rxt(137)*y(16)*y(9) &
                  +rxt(16)*y(17)
         loss(48) = (rxt(146)* y(1) +rxt(145)* y(4) + (rxt(147) +rxt(148) +rxt(149)) &
                 * y(23) + het_rates(21))* y(21)
         prod(48) = (rxt(144)*y(19) +rxt(155)*y(20) +rxt(141)*y(18) +rxt(150)*y(2) + &
                 rxt(213)*y(59) +rxt(214)*y(60) +rxt(217)*y(61))*y(22) &
                  + (rxt(109)*y(15) +rxt(111)*y(20))*y(3) + (rxt(19) + &
                 2.000_r8*rxt(21))*y(82) +rxt(16)*y(17) +2.000_r8*rxt(17)*y(18) &
                  +rxt(162)*y(27)*y(20) +rxt(28)*y(32)
         loss(55) = (rxt(158)* y(1) +rxt(157)* y(2) +rxt(120)* y(9) +rxt(133)* y(10) &
                  +rxt(132)* y(11) +rxt(138)* y(16) + (rxt(147) +rxt(148) +rxt(149)) &
                 * y(21) +rxt(152)* y(22) + 2._r8*rxt(159)* y(23) + (rxt(164) + &
                 rxt(165))* y(27) +rxt(171)* y(29) +rxt(188)* y(36) +rxt(192)* y(37) &
                  + rxt(206) + het_rates(23))* y(23)
         prod(55) = (rxt(143)*y(19) +rxt(208)*y(41) +rxt(211)*y(42) +rxt(131)*y(11) + &
                 rxt(151)*y(1) +rxt(156)*y(24) +rxt(169)*y(29) +rxt(191)*y(37) + &
                 rxt(224)*y(62) +.500_r8*rxt(226)*y(65))*y(22) + (rxt(140)*y(11) + &
                 rxt(142)*y(2) +rxt(166)*y(27) +rxt(189)*y(36))*y(18) + (rxt(15) + &
                 rxt(135))*y(13) + (rxt(160)*y(2) +rxt(163)*y(27))*y(24) &
                  +rxt(109)*y(15)*y(3) +rxt(145)*y(21)*y(4) +rxt(137)*y(16)*y(9) &
                  +rxt(207)*y(41)*y(27)
         loss(35) = (rxt(160)* y(2) +rxt(156)* y(22) +rxt(163)* y(27) + rxt(22) &
                  + het_rates(24))* y(24)
         prod(35) = (.500_r8*rxt(206) +rxt(159)*y(23))*y(23) +rxt(154)*y(22)*y(22)
         loss(50) = (rxt(96)* y(3) +rxt(225)* y(63) + rxt(19) + rxt(20) + rxt(21) &
                  + het_rates(82))* y(82)
         prod(50) = (rxt(136)*y(15) +rxt(155)*y(20) +rxt(208)*y(41) +rxt(209)*y(48) + &
                 rxt(210)*y(46) +rxt(211)*y(42) +rxt(128)*y(12) +rxt(134)*y(13) + &
                 rxt(139)*y(17) +rxt(141)*y(18) +rxt(152)*y(23) +rxt(153)*y(22) + &
                 rxt(156)*y(24) +rxt(179)*y(32) +rxt(183)*y(33) +rxt(199)*y(38))*y(22) &
                  + (rxt(233)*y(33) +rxt(234)*y(39) +rxt(238)*y(33) +rxt(244)*y(33) + &
                 rxt(245)*y(39))*y(32) +rxt(149)*y(23)*y(21) +rxt(49)*y(64)
         loss(62) = (rxt(161)* y(1) +rxt(167)* y(15) +rxt(166)* y(18) +rxt(162)* y(20) &
                  + (rxt(164) +rxt(165))* y(23) +rxt(163)* y(24) +rxt(182)* y(33) &
                  +rxt(186)* y(34) +rxt(207)* y(41) + het_rates(27))* y(27)
         prod(62) = (rxt(24) +rxt(168)*y(2) +rxt(169)*y(22) +rxt(172)*y(9) + &
                 2.000_r8*rxt(174)*y(29) +rxt(176)*y(29) +rxt(196)*y(37) + &
                 rxt(221)*y(61))*y(29) + (3.000_r8*rxt(100)*y(43) + &
                 2.000_r8*rxt(101)*y(44) +3.000_r8*rxt(102)*y(45) +rxt(103)*y(46) + &
                 4.000_r8*rxt(104)*y(47) +rxt(112)*y(32))*y(3) + (rxt(208)*y(41) + &
                 3.000_r8*rxt(209)*y(48) +rxt(210)*y(46) +rxt(179)*y(32))*y(22) &
                  + (rxt(28) +rxt(180)*y(2))*y(32) +2.000_r8*rxt(23)*y(28) &
                  +2.000_r8*rxt(26)*y(31) +rxt(27)*y(33) +rxt(29)*y(34) +rxt(31)*y(35)
         loss(21) = ( + rxt(23) + het_rates(28))* y(28)
         prod(21) = (rxt(232)*y(34) +rxt(233)*y(33) +rxt(237)*y(34) +rxt(238)*y(33) + &
                 rxt(243)*y(34) +rxt(244)*y(33))*y(32) +rxt(186)*y(34)*y(27) &
                  +rxt(175)*y(29)*y(29)
         loss(53) = (rxt(168)* y(2) +rxt(172)* y(9) +rxt(173)* y(10) + (rxt(169) + &
                 rxt(170))* y(22) +rxt(171)* y(23) + 2._r8*(rxt(174) +rxt(175) + &
                 rxt(176) +rxt(177))* y(29) + (rxt(195) +rxt(196) +rxt(197))* y(37) &
                  +rxt(221)* y(61) + rxt(24) + het_rates(29))* y(29)
         prod(53) = (rxt(161)*y(1) +rxt(165)*y(23) +rxt(182)*y(33))*y(27) &
                  + (rxt(181)*y(33) +rxt(184)*y(34))*y(2) + (rxt(25) +rxt(223)*y(61)) &
                 *y(30) +rxt(183)*y(33)*y(22) +2.000_r8*rxt(178)*y(31) +rxt(30)*y(34)
         loss(27) = (rxt(223)* y(61) + rxt(25) + het_rates(30))* y(30)
         prod(27) = (rxt(176)*y(29) +rxt(195)*y(37))*y(29)
         loss(17) = ( + rxt(26) + rxt(178) + het_rates(31))* y(31)
         prod(17) =rxt(177)*y(29)*y(29)
         loss(63) = (rxt(180)* y(2) +rxt(112)* y(3) +rxt(179)* y(22) + (rxt(233) + &
                 rxt(238) +rxt(244))* y(33) + (rxt(232) +rxt(237) +rxt(243))* y(34) &
                  + (rxt(234) +rxt(245))* y(39) + rxt(28) + het_rates(32))* y(32)
         prod(63) = (rxt(162)*y(20) +rxt(167)*y(15) +2.000_r8*rxt(207)*y(41) + &
                 rxt(163)*y(24) +rxt(164)*y(23) +rxt(166)*y(18) +rxt(182)*y(33))*y(27) &
                  +rxt(170)*y(29)*y(22)
         loss(45) = (rxt(181)* y(2) +rxt(183)* y(22) +rxt(182)* y(27) + (rxt(233) + &
                 rxt(238) +rxt(244))* y(32) + rxt(27) + het_rates(33))* y(33)
         prod(45) = (rxt(230) +rxt(236) +rxt(241) +rxt(185)*y(22))*y(34) &
                  +rxt(171)*y(29)*y(23)
         loss(46) = (rxt(184)* y(2) +rxt(185)* y(22) +rxt(186)* y(27) + (rxt(232) + &
                 rxt(237) +rxt(243))* y(32) + rxt(29) + rxt(30) + rxt(230) + rxt(236) &
                  + rxt(241) + het_rates(34))* y(34)
         prod(46) =rxt(173)*y(29)*y(10)
         loss(22) = ( + rxt(31) + het_rates(35))* y(35)
         prod(22) = (rxt(234)*y(39) +rxt(245)*y(39))*y(32) +rxt(197)*y(37)*y(29)
         loss(64) = (rxt(187)* y(1) +rxt(189)* y(18) +rxt(188)* y(23) + het_rates(36)) &
                 * y(36)
         prod(64) = (rxt(32) +rxt(190)*y(2) +rxt(191)*y(22) +rxt(193)*y(9) + &
                 rxt(195)*y(29) +rxt(196)*y(29) +2.000_r8*rxt(198)*y(37) + &
                 rxt(222)*y(61))*y(37) + (rxt(105)*y(42) +rxt(106)*y(50) + &
                 rxt(107)*y(49) +rxt(113)*y(38))*y(3) + (rxt(211)*y(42) + &
                 rxt(199)*y(38))*y(22) +rxt(200)*y(38)*y(2) +rxt(31)*y(35) +rxt(33) &
                 *y(39) +rxt(34)*y(40)
         loss(52) = (rxt(190)* y(2) +rxt(193)* y(9) +rxt(194)* y(10) +rxt(191)* y(22) &
                  +rxt(192)* y(23) + (rxt(195) +rxt(196) +rxt(197))* y(29) &
                  + 2._r8*rxt(198)* y(37) +rxt(222)* y(61) + rxt(32) + het_rates(37)) &
                 * y(37)
         prod(52) = (rxt(201)*y(39) +rxt(202)*y(40))*y(2) +rxt(187)*y(36)*y(1) &
                  +rxt(35)*y(40)
         loss(37) = (rxt(200)* y(2) +rxt(113)* y(3) +rxt(199)* y(22) + het_rates(38)) &
                 * y(38)
         prod(37) = (rxt(188)*y(23) +rxt(189)*y(18))*y(36)
         loss(39) = (rxt(201)* y(2) + (rxt(234) +rxt(245))* y(32) + rxt(33) &
                  + het_rates(39))* y(39)
         prod(39) = (rxt(231) +rxt(239) +rxt(242))*y(40) +rxt(192)*y(37)*y(23)
         loss(33) = (rxt(202)* y(2) + rxt(34) + rxt(35) + rxt(231) + rxt(239) &
                  + rxt(242) + het_rates(40))* y(40)
         prod(33) =rxt(194)*y(37)*y(10)
         loss(34) = ((rxt(248) +rxt(257))* y(2) +rxt(255)* y(4) +rxt(260)* y(57) &
                  + het_rates(52))* y(52)
         prod(34) = 0._r8
         loss(40) = (rxt(250)* y(8) +rxt(251)* y(9) +rxt(259)* y(57) + rxt(256) &
                  + het_rates(53))* y(53)
         prod(40) = (rxt(58) +rxt(70) +rxt(246)*y(55) +rxt(252)*y(54) +rxt(255)*y(52)) &
                 *y(4) +rxt(249)*y(55)*y(51)
         loss(29) = (rxt(254)* y(2) + (rxt(252) +rxt(253))* y(4) + het_rates(54)) &
                 * y(54)
         prod(29) =rxt(57)*y(8)
         loss(31) = (rxt(246)* y(4) +rxt(249)* y(51) + rxt(247) + het_rates(55)) &
                 * y(55)
         prod(31) = (rxt(54) +rxt(55) +rxt(56) +rxt(67) +rxt(68) +rxt(69) + &
                 rxt(254)*y(54) +rxt(257)*y(52))*y(2) + (rxt(60) +rxt(61) +rxt(62) + &
                 rxt(72) +rxt(73) +rxt(74))*y(4)
         loss(42) = (rxt(258)* y(57) + het_rates(56))* y(56)
         prod(42) = (rxt(256) +rxt(250)*y(8) +rxt(251)*y(9))*y(53) +rxt(248)*y(52) &
                 *y(2) +rxt(253)*y(54)*y(4) +rxt(7)*y(9) +rxt(247)*y(55)
         loss(32) = (rxt(115)* y(2) +rxt(114)* y(4) + het_rates(58))* y(58)
         prod(32) = (rxt(248)*y(2) +.900_r8*rxt(260)*y(57))*y(52) &
                  +.800_r8*rxt(258)*y(57)*y(56)
         loss(41) = (rxt(260)* y(52) +rxt(259)* y(53) +rxt(258)* y(56) &
                  + het_rates(57))* y(57)
         prod(41) = (rxt(58) +rxt(60) +rxt(61) +rxt(62) +rxt(70) +rxt(72) +rxt(73) + &
                 rxt(74))*y(4) + (rxt(54) +rxt(55) +rxt(56) +rxt(67) +rxt(68) + &
                 rxt(69))*y(2) +rxt(57)*y(8) +rxt(7)*y(9)
         loss(26) = (rxt(212)* y(2) +rxt(213)* y(22) + rxt(52) + het_rates(59))* y(59)
         prod(26) = 0._r8
         loss(36) = (rxt(216)* y(1) +rxt(215)* y(4) +rxt(214)* y(22) + het_rates(60)) &
                 * y(60)
         prod(36) =rxt(52)*y(59) +rxt(53)*y(61)
         loss(65) = (rxt(219)* y(1) +rxt(218)* y(4) +rxt(220)* y(10) +rxt(217)* y(22) &
                  +rxt(221)* y(29) +rxt(223)* y(30) +rxt(222)* y(37) + rxt(53) &
                  + het_rates(61))* y(61)
         prod(65) = (rxt(214)*y(22) +rxt(215)*y(4) +rxt(216)*y(1))*y(60) &
                  +rxt(212)*y(59)*y(2) +rxt(50)*y(62)
         loss(47) = (rxt(224)* y(22) + rxt(50) + het_rates(62))* y(62)
         prod(47) = (rxt(217)*y(22) +rxt(218)*y(4) +rxt(219)*y(1) +rxt(220)*y(10) + &
                 rxt(221)*y(29) +rxt(222)*y(37) +rxt(223)*y(30))*y(61) &
                  + (rxt(213)*y(59) +.500_r8*rxt(226)*y(65) +rxt(227)*y(65))*y(22) &
                  +rxt(228)*y(65)*y(11) +rxt(51)*y(63)
         loss(23) = (rxt(225)* y(82) + rxt(51) + het_rates(63))* y(63)
         prod(23) =rxt(224)*y(62)*y(22) +rxt(49)*y(64)
         loss(18) = ( + rxt(49) + het_rates(64))* y(64)
         prod(18) =rxt(225)*y(82)*y(63)
         loss(24) = (rxt(228)* y(11) + (rxt(226) +rxt(227))* y(22) + het_rates(65)) &
                 * y(65)
         prod(24) = 0._r8
         loss(1) = ( + het_rates(66))* y(66)
         prod(1) = 0._r8
         loss(2) = ( + het_rates(67))* y(67)
         prod(2) = 0._r8
         loss(3) = ( + het_rates(68))* y(68)
         prod(3) = 0._r8
         loss(4) = ( + het_rates(69))* y(69)
         prod(4) = 0._r8
         loss(5) = ( + het_rates(70))* y(70)
         prod(5) = 0._r8
         loss(6) = ( + het_rates(71))* y(71)
         prod(6) = 0._r8
         loss(7) = ( + het_rates(72))* y(72)
         prod(7) = 0._r8
         loss(8) = ( + het_rates(73))* y(73)
         prod(8) = 0._r8
         loss(9) = ( + het_rates(74))* y(74)
         prod(9) = 0._r8
         loss(10) = ( + het_rates(75))* y(75)
         prod(10) = 0._r8
         loss(11) = ( + het_rates(76))* y(76)
         prod(11) = 0._r8
         loss(12) = ( + het_rates(77))* y(77)
         prod(12) = 0._r8
         loss(13) = ( + het_rates(78))* y(78)
         prod(13) = 0._r8
         loss(14) = ( + het_rates(79))* y(79)
         prod(14) = 0._r8
         loss(15) = ( + het_rates(80))* y(80)
         prod(15) = 0._r8
         loss(16) = ( + het_rates(81))* y(81)
         prod(16) = 0._r8

      end subroutine imp_prod_loss

      end module mo_prod_loss
