!*************************************************************************
!	> File Name: micro_mg1_0_cpe.F90
!	> Author: Xu Kai
!	> Created Time: 2018年11月28日 星期三 19时08分36秒
! ************************************************************************/

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        subroutine shr_spfn_gamma_nonintrinsic_r8_v2(x, gamma)
        
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !
        ! 7 Feb 2013 -- S. Santos
        ! The following comments are from the original version. Changes have
        ! been made to update syntax and allow inclusion into this module.
        !
        !----------------------------------------------------------------------
        !
        ! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
        !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
        !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
        !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
        !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
        !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
        !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
        !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
        !   MACHINE-DEPENDENT CONSTANTS.
        !
        !
        !*******************************************************************
        !*******************************************************************
        !
        ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
        !
        ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
        ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
        ! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
        !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
        !                  GAMMA(XBIG) = BETA**MAXEXP
        ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
        !          APPROXIMATELY BETA**MAXEXP
        ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
        !          1.0+EPS .GT. 1.0
        ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
        !          1/XMININ IS MACHINE REPRESENTABLE
        !
        !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
        !
        !                            BETA       MAXEXP        XBIG
        !
        ! CRAY-1         (S.P.)        2         8191        966.961
        ! CYBER 180/855
        !   UNDER NOS    (S.P.)        2         1070        177.803
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (S.P.)        2          128        35.040
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (D.P.)        2         1024        171.624
        ! IBM 3033       (D.P.)       16           63        57.574
        ! VAX D-FORMAT   (D.P.)        2          127        34.844
        ! VAX G-FORMAT   (D.P.)        2         1023        171.489
        !
        !                            XINF         EPS        XMININ
        !
        ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
        ! CYBER 180/855
        !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
        ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
        ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
        ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
        !
        !*******************************************************************
        !*******************************************************************
        !
        ! ERROR RETURNS
        !
        !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
        !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
        !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
        !
        !
        !  INTRINSIC FUNCTIONS REQUIRED ARE:
        !
        !     INT, DBLE, EXP, LOG, REAL, SIN
        !
        !
        ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
        !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
        !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
        !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
        !
        !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
        !              SONS, NEW YORK, 1968.
        !
        !  LATEST MODIFICATION: OCTOBER 12, 1989
        !
        !  AUTHORS: W. J. CODY AND L. STOLTZ
        !           APPLIED MATHEMATICS DIVISION
        !           ARGONNE NATIONAL LABORATORY
        !           ARGONNE, IL 60439
        !
        !----------------------------------------------------------------------
        implicit none
        
          real(8), intent(in) :: x
          real(8), intent(out) :: gamma
          real(8) :: fact, res, sum, xden, xnum, y, y1, ysq, z
        
          integer :: i, n
          logical :: negative_odd
        
          ! log(2*pi)/2
          real(8), parameter :: logsqrt2pi = 0.9189385332046727417803297d0


          real(8), parameter :: pi = 3.1415926535897932384626434d0
          real(8), parameter :: xinfr8 = huge(1.d0)
          real(8), parameter :: epsr8 = epsilon(1.d0)
          real(8), parameter :: xminr8 = tiny(1.d0)
          real(8), parameter :: xbig_gamma = 171.624d0
          
        !----------------------------------------------------------------------
        !  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
        !     APPROXIMATION OVER (1,2).
        !----------------------------------------------------------------------
          real(8), parameter :: P(8) = &
               (/-1.71618513886549492533811d+0, 2.47656508055759199108314d+1, &
                 -3.79804256470945635097577d+2, 6.29331155312818442661052d+2, &
                  8.66966202790413211295064d+2,-3.14512729688483675254357d+4, &
                 -3.61444134186911729807069d+4, 6.64561438202405440627855d+4 /)
          real(8), parameter :: Q(8) = &
               (/-3.08402300119738975254353d+1, 3.15350626979604161529144d+2, &
                 -1.01515636749021914166146d+3,-3.10777167157231109440444d+3, &
                  2.25381184209801510330112d+4, 4.75584627752788110767815d+3, &
                 -1.34659959864969306392456d+5,-1.15132259675553483497211d+5 /)
        !----------------------------------------------------------------------
        !  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
        !----------------------------------------------------------------------
          real(8), parameter :: C(7) = &
               (/-1.910444077728d-03,          8.4171387781295d-04, &
                 -5.952379913043012d-04,       7.93650793500350248d-04, &
                 -2.777777777777681622553d-03, 8.333333333333333331554247d-02, &
                  5.7083835261d-03 /)
          
          integer(8) :: sin_agent, log_agent, exp_agent
          real(8) :: tmp1, tmp2
        
          call math_agent_sin_c(sin_agent)
          call math_agent_log_c(log_agent)
          call math_agent_exp_c(exp_agent)

          negative_odd = .false.
          fact = 1.d0
          n = 0
          y = x
          if (y <= 0.d0) then
        !----------------------------------------------------------------------
        !  ARGUMENT IS NEGATIVE
        !----------------------------------------------------------------------
             y = -x
             y1 = aint(y)
             res = y - y1
             if (res /= 0.d0) then
                negative_odd = (y1 /= aint(y1*0.5d0)*2.d0)
                tmp1 = pi*res
                call math_agent_1i1o(sin_agent, tmp1, tmp2)
                fact = -pi / tmp2
                !fact = -pi/sin(pi*res)
                y = y + 1.d0
             else
                gamma = xinfr8
                goto 500
             end if
          end if
        !----------------------------------------------------------------------
        !  ARGUMENT IS POSITIVE
        !----------------------------------------------------------------------
          if (y < epsr8) then
        !----------------------------------------------------------------------
        !  ARGUMENT .LT. EPS
        !----------------------------------------------------------------------
             if (y >= xminr8) then
                res = 1.d0/y
             else
                gamma = xinfr8
                goto 500
             end if
          elseif (y < 12.d0) then
             y1 = y
             if (y < 1.d0) then
        !----------------------------------------------------------------------
        !  0.0 .LT. ARGUMENT .LT. 1.0
        !----------------------------------------------------------------------
                z = y
                y = y + 1.d0
             else
        !----------------------------------------------------------------------
        !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
        !----------------------------------------------------------------------
                n = int(y) - 1
                y = y - real(n, 8)
                z = y - 1.d0
             end if
        !----------------------------------------------------------------------
        !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
        !----------------------------------------------------------------------
             xnum = 0.d0
             xden = 1.d0
             do i=1,8
                xnum = (xnum+P(i))*z
                xden = xden*z + Q(i)
             end do
             res = xnum/xden + 1.d0
             if (y1 < y) then
        !----------------------------------------------------------------------
        !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
        !----------------------------------------------------------------------
                res = res/y1
             elseif (y1 > y) then
        !----------------------------------------------------------------------
        !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
        !----------------------------------------------------------------------
                do i = 1,n
                   res = res*y
                   y = y + 1.d0
                end do
             end if
          else
        !----------------------------------------------------------------------
        !  EVALUATE FOR ARGUMENT .GE. 12.0,
        !----------------------------------------------------------------------
             if (y <= xbig_gamma) then
                ysq = y*y
                sum = C(7)
                do i=1,6
                   sum = sum/ysq + C(i)
                end do
                sum = sum/y - y + logsqrt2pi
                call math_agent_1i1o(log_agent, y, tmp1)
                sum = sum + (y-0.5d0)*tmp1
                !sum = sum + (y-0.5d0)*log(y)
                call math_agent_1i1o(exp_agent, sum, res)
                !res = exp(sum)
             else
                gamma = xinfr8
                goto 500
             end if
          end if
        !----------------------------------------------------------------------
        !  FINAL ADJUSTMENTS AND RETURN
        !----------------------------------------------------------------------
          if (negative_odd)  res = -res
          if (fact /= 1.d0) res = fact/res
          gamma = res
          500 continue
        ! ---------- LAST LINE OF GAMMA ----------
        end subroutine shr_spfn_gamma_nonintrinsic_r8_v2

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        function shr_spfn_gamma_nonintrinsic_r8(X,math_agent,P,Q,C,S) result(gamma)
        
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !
        ! 7 Feb 2013 -- S. Santos
        ! The following comments are from the original version. Changes have
        ! been made to update syntax and allow inclusion into this module.
        !
        !----------------------------------------------------------------------
        !
        ! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
        !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
        !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
        !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
        !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
        !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
        !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
        !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
        !   MACHINE-DEPENDENT CONSTANTS.
        !
        !
        !*******************************************************************
        !*******************************************************************
        !
        ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
        !
        ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
        ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
        ! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
        !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
        !                  GAMMA(XBIG) = BETA**MAXEXP
        ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
        !          APPROXIMATELY BETA**MAXEXP
        ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
        !          1.0+EPS .GT. 1.0
        ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
        !          1/XMININ IS MACHINE REPRESENTABLE
        !
        !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
        !
        !                            BETA       MAXEXP        XBIG
        !
        ! CRAY-1         (S.P.)        2         8191        966.961
        ! CYBER 180/855
        !   UNDER NOS    (S.P.)        2         1070        177.803
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (S.P.)        2          128        35.040
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (D.P.)        2         1024        171.624
        ! IBM 3033       (D.P.)       16           63        57.574
        ! VAX D-FORMAT   (D.P.)        2          127        34.844
        ! VAX G-FORMAT   (D.P.)        2         1023        171.489
        !
        !                            XINF         EPS        XMININ
        !
        ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
        ! CYBER 180/855
        !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
        ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
        ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
        ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
        !
        !*******************************************************************
        !*******************************************************************
        !
        ! ERROR RETURNS
        !
        !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
        !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
        !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
        !
        !
        !  INTRINSIC FUNCTIONS REQUIRED ARE:
        !
        !     INT, DBLE, EXP, LOG, REAL, SIN
        !
        !
        ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
        !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
        !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
        !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
        !
        !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
        !              SONS, NEW YORK, 1968.
        !
        !  LATEST MODIFICATION: OCTOBER 12, 1989
        !
        !  AUTHORS: W. J. CODY AND L. STOLTZ
        !           APPLIED MATHEMATICS DIVISION
        !           ARGONNE NATIONAL LABORATORY
        !           ARGONNE, IL 60439
        !
        !----------------------------------------------------------------------
        implicit none
        
          real(8), intent(in) :: x
          real(8) :: gamma
          real(8) :: fact, res, sum, xden, xnum, y, y1, ysq, z
        
          integer :: i, n
          logical :: negative_odd
        
          !! log(2*pi)/2
          !real(8) :: logsqrt2pi = 0.9189385332046727417803297d0


          !real(8) :: pi = 3.1415926535897932384626434d0
          !real(8) :: xinfr8 = huge(1.d0)
          !real(8) :: epsr8 = epsilon(1.d0)
          !real(8) :: xminr8 = tiny(1.d0)
          !real(8) :: xbig_gamma = 171.624d0
          
        !!----------------------------------------------------------------------
        !!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
        !!     APPROXIMATION OVER (1,2).
        !!----------------------------------------------------------------------
        !  real(8) :: P(8) = &
        !       (/-1.71618513886549492533811d+0, 2.47656508055759199108314d+1, &
        !         -3.79804256470945635097577d+2, 6.29331155312818442661052d+2, &
        !          8.66966202790413211295064d+2,-3.14512729688483675254357d+4, &
        !         -3.61444134186911729807069d+4, 6.64561438202405440627855d+4 /)
        !  real(8) :: Q(8) = &
        !       (/-3.08402300119738975254353d+1, 3.15350626979604161529144d+2, &
        !         -1.01515636749021914166146d+3,-3.10777167157231109440444d+3, &
        !          2.25381184209801510330112d+4, 4.75584627752788110767815d+3, &
        !         -1.34659959864969306392456d+5,-1.15132259675553483497211d+5 /)
        !!----------------------------------------------------------------------
        !!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
        !!----------------------------------------------------------------------
        !  real(8) :: C(7) = &
        !       (/-1.910444077728d-03,          8.4171387781295d-04, &
        !         -5.952379913043012d-04,       7.93650793500350248d-04, &
        !         -2.777777777777681622553d-03, 8.333333333333333331554247d-02, &
        !          5.7083835261d-03 /)

        
          real(8), intent(in) :: P(8), Q(8), C(7), S(6)
          integer(8), intent(in) :: math_agent(3)
          real(8) :: tmp1, tmp2, tmp3, tmp4
          !integer(8) :: sin_agent, log_agent, exp_agent
          !real(8) :: logsqrt2pi
          !real(8) :: pi
          !real(8) :: xinfr8
          !real(8) :: epsr8
          !real(8) :: xminr8
          !real(8) :: xbig_gamma

          !logsqrt2pi    = S(1)   
          !pi            = S(2) 
          !xinfr8        = S(3) 
          !epsr8         = S(4) 
          !xminr8        = S(5) 
          !xbig_gamma    = S(6)

          !sin_agent = math_agent(1)
          !log_agent = math_agent(2)
          !exp_agent = math_agent(3)
        
          !call math_agent_sin_c(sin_agent)
          !call math_agent_log_c(log_agent)
          !call math_agent_exp_c(exp_agent)

          negative_odd = .false.
          fact = 1.d0
          n = 0
          y = x
          if (y <= 0.d0) then
        !----------------------------------------------------------------------
        !  ARGUMENT IS NEGATIVE
        !----------------------------------------------------------------------
             y = -x
             y1 = aint(y)
             res = y - y1
             if (res /= 0.d0) then
                negative_odd = (y1 /= aint(y1*0.5d0)*2.d0)
                tmp1 = S(2)*res
                call math_agent_1i1o(math_agent(1), tmp1, tmp2)
                fact = -S(2) / tmp2
                !fact = -pi/sin(pi*res)
                y = y + 1.d0
             else
                gamma = S(3)
                return
             end if
          end if
        !----------------------------------------------------------------------
        !  ARGUMENT IS POSITIVE
        !----------------------------------------------------------------------
          if (y < S(4)) then
        !----------------------------------------------------------------------
        !  ARGUMENT .LT. EPS
        !----------------------------------------------------------------------
             if (y >= S(5)) then
                res = 1.d0/y
             else
                gamma = S(3)
                return
             end if
          elseif (y < 12.d0) then
             y1 = y
             if (y < 1.d0) then
        !----------------------------------------------------------------------
        !  0.0 .LT. ARGUMENT .LT. 1.0
        !----------------------------------------------------------------------
                z = y
                y = y + 1.d0
             else
        !----------------------------------------------------------------------
        !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
        !----------------------------------------------------------------------
                n = int(y) - 1
                y = y - real(n, 8)
                z = y - 1.d0
             end if
        !----------------------------------------------------------------------
        !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
        !----------------------------------------------------------------------
             xnum = 0.d0
             xden = 1.d0
             do i=1,8
                xnum = (xnum+P(i))*z
                xden = xden*z + Q(i)
             end do
             res = xnum/xden + 1.d0
             if (y1 < y) then
        !----------------------------------------------------------------------
        !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
        !----------------------------------------------------------------------
                res = res/y1
             elseif (y1 > y) then
        !----------------------------------------------------------------------
        !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
        !----------------------------------------------------------------------
                do i = 1,n
                   res = res*y
                   y = y + 1.d0
                end do
             end if
          else
        !----------------------------------------------------------------------
        !  EVALUATE FOR ARGUMENT .GE. 12.0,
        !----------------------------------------------------------------------
             if (y <= S(6)) then
                ysq = y*y
                sum = C(7)
                do i=1,6
                   sum = sum/ysq + C(i)
                end do
                sum = sum/y - y + S(1)
                call math_agent_1i1o(math_agent(2), y, tmp1) ! log10(y)
                call math_agent_1i1o(math_agent(3), 1.d0, tmp2) ! e^1
                call math_agent_1i1o(math_agent(2), tmp2, tmp3) ! log10(e)
                tmp4 = tmp1 / tmp3
                sum = sum + (y-0.5d0)*tmp4
                !sum = sum + (y-0.5d0)*log(y)
                call math_agent_1i1o(math_agent(3), sum, res)
                !res = exp(sum)
             else
                gamma = S(3)
                return
             end if
          end if
        !----------------------------------------------------------------------
        !  FINAL ADJUSTMENTS AND RETURN
        !----------------------------------------------------------------------
          if (negative_odd)  res = -res
          if (fact /= 1.d0) res = fact/res
          gamma = res
        ! ---------- LAST LINE OF GAMMA ----------
        end function shr_spfn_gamma_nonintrinsic_r8

        !subroutine test ( nacon, qrout2, qsout2, drout2, dsout2, pver)
        !implicit none
        !real(8), intent(in) :: rndst(pver, 4), nacon(pver, 4)
        !real(8), intent(out) :: drout2(pver), dsout2(pver)
        !integer, intent(in) :: pver
        !integer :: k
        !do k = 1,pver
        !qrout2(k) = drout2(k) + nacon(k,1)
        !qsout2(k) = drout2(k) + nacon(k,2)
        !drout2(k) = drout2(k) + nacon(k,3)
        !dsout2(k) = drout2(k) + nacon(k,4)
        !enddo
        !end subroutine

        function svp_water(t, lg10, power, tboil) result(es)
          real(8), intent(in) :: t  ! Temperature in Kelvin
          real(8) :: es             ! SVP in Pa

          integer(8), intent(in) :: lg10, power
          real(8), intent(in) :: tboil
          real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5

  !call lwpf_start_test_unit_kernel_a()
          tmp1 = 10.d0
          tmp2 = tboil/t

          call math_agent_1i1o(lg10, tmp2, tmp3)
          tmp4 = 11.344d0*(1.d0-t/tboil)
          call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          tmp4 = (-3.49149d0*(tboil/t-1.d0))
          call math_agent_2i1o(power, tmp1, tmp4, tmp5)
          tmp4 = 1013.246d0
          call math_agent_1i1o(lg10, tmp4, tmp1)

          tmp4 = -7.90298d0*(tboil/t-1.d0)+ &
               5.02808d0*tmp3- &
               1.3816d-7*(tmp2-1.d0)+ &
               8.1328d-3*(tmp5-1.d0)+ &
               tmp1

          tmp1 = 10.d0
          call math_agent_2i1o(power, tmp1, tmp4, tmp2)

          es = tmp2 * 100.d0

    !es = 10.d0**(-7.90298d0*(tboil/t-1.d0)+ &
    !           5.02808d0*log10(tboil/t)- &
    !           1.3816d-7*(10.d0**(11.344d0*(1.d0-t/tboil))-1.d0)+ &
    !           8.1328d-3*(10.d0**(-3.49149d0*(tboil/t-1.d0))-1.d0)+ &
    !           log10(1013.246d0))*100.d0
  !call lwpf_stop_test_unit_kernel_a()

        
        end function svp_water

        function svp_ice(t,lg10,power,h2otrip) result(es)
          real(8), intent(in) :: t  ! Temperature in Kelvin
          real(8) :: es             ! SVP in Pa

          integer(8),intent(in) :: lg10, power
          real(8), intent(in) :: h2otrip
          real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5
        
  !call lwpf_start_test_unit_kernel_a()
          tmp1 = h2otrip/t
          call math_agent_1i1o(lg10, tmp1, tmp2)
          tmp1 = 6.1071d0
          call math_agent_1i1o(lg10, tmp1, tmp3)
          tmp1 = -9.09718d0*(h2otrip/t-1.d0)-3.56654d0* &
               tmp2+0.876793d0*(1.d0-t/h2otrip)+ &
               tmp3
          tmp2 = 10.d0
          call math_agent_2i1o(power, tmp2, tmp1, tmp3)
          es = tmp3*100.d0

          !es = 10.d0**(-9.09718d0*(h2otrip/t-1.d0)-3.56654d0* &
          !     log10(h2otrip/t)+0.876793d0*(1.d0-t/h2otrip)+ &
          !     log10(6.1071d0))*100.d0
        
  !call lwpf_stop_test_unit_kernel_a()
        end function svp_ice

        function svp_to_qsat(es, p, epsilo, omeps) result(qs)
        
          real(8), intent(in) :: es  ! SVP
          real(8), intent(in) :: p   ! Current pressure.
          real(8), intent(in) :: epsilo, omeps
          real(8) :: qs
        
          ! If pressure is less than SVP, set qs to maximum of 1.
          if ( (p - es) <= 0.d0 ) then
             qs = 1.0d0
          else
             qs = epsilo*es / (p - omeps*es)
          end if
        
        end function svp_to_qsat
        subroutine micro_mg1_0_compute ( &
             g, r, rv, cpp, rhow, tmelt, xxlv, xlf, &
             xxls, rhosn, rhoi, ac, bc, as, bs, ai, bi, ar, &   
             br, ci, di, cs, ds, cr, dr, f1s, f2s, Eii, Ecr, &      
             f1r, f2r, DCS, qsmall, bimm, aimm, rhosu, mi0, &
             rin, pi, cons1, cons4, cons5, cons6, cons7, &           
             cons8, cons11, cons13, cons14, cons16, cons17, &       
             cons22, cons23, cons24, cons25, cons27, cons28, &
             lammini, lammaxi, lamminr, lammaxr, &            
             lammins, lammaxs, tmax_fsnow, tmin_fsnow, tt0, csmin, &
             csmax, minrefl, mindbz, rhmini, use_hetfrz_classnuc, &
             microp_uniform, pcols, pver, ncol, top_lev, deltatin,&
             tn, qn, qc, qi, nc,                              &
             ni, p, pdel, cldn, liqcldf,                      &
             relvar, accre_enhan,                             &
             icecldf, rate1ord_cw2pr_st, naai, npccnin,       &
             rndst, nacon, tlat, qvlat, qctend,               &
             qitend, nctend, nitend, effc, effc_fn,           &
             effi, prect, preci, nevapr, evapsnow, am_evp_st, &
             prain, prodsnow, cmeout, deffi, pgamrad,         &
             lamcrad, qsout, dsout, rflx, sflx,               &
             qrout, reff_rain, reff_snow, qcsevap, qisevap,   &
             qvres, cmeiout, vtrmc, vtrmi, qcsedten,          &
             qisedten, prao, prco, mnuccco, mnuccto,          &
             msacwio, psacwso, bergso, bergo, melto,          &
             homoo, qcreso, prcio, praio, qireso,             &
             mnuccro, pracso, meltsdt, frzrdt, mnuccdo,       &
             nrout, nsout, refl, arefl, areflz,               &
             frefl, csrfl, acsrfl, fcsrfl, rercld,            &
             ncai, ncal, qrout2, qsout2, nrout2,              &
             nsout2, drout2, dsout2, freqs, freqr,            &
             nfice, prer_evap, do_cldice,                     &
             tnd_qsnow, tnd_nsnow, re_ice,                    &
             frzimm, frzcnt, frzdep,                          &
             svp_tmelt, svp_h2otrip, svp_tboil,               &
             svp_ttrice, svp_epsilo, myid                     )


        !use wv_sat_methods_cpe, only: &
        !    svp_water => wv_sat_svp_water, &
        !    svp_ice => wv_sat_svp_ice, &
        !    svp_to_qsat => wv_sat_svp_to_qsat, &
        !    wv_sat_methods_init
        !    !shr_spfn_gamma_nonintrinsic_r8

        implicit none

        interface
            real(8) function shr_spfn_gamma_nonintrinsic_r8(x,math_agent,P,Q,C,S)
            real(8) :: x, P(8), Q(8), C(7), S(6)
            integer(8) :: math_agent(3)
            end function
            real(8) function svp_water(t, lg10, power, svp_tboil)
                real(8) :: t, svp_tboil
                integer(8) :: lg10, power
            end function
            real(8) function svp_ice(t, lg10, power, svp_h2otrip)
                real(8) :: t, svp_h2otrip
                integer(8) :: lg10, power
            end function
            real(8) function svp_to_qsat(es, p, svp_epsilo, svp_omeps)
                real(8) :: es, p, svp_epsilo, svp_omeps
            end function
        end interface

        real(8), intent(in) :: g              !gravity
        real(8), intent(in) :: r              !Dry air Gas constant
        real(8), intent(in) :: rv             !water vapor gas contstant
        real(8), intent(in) :: cpp            !specific heat of dry air
        real(8), intent(in) :: rhow           !density of liquid water
        real(8), intent(in) :: tmelt          ! Freezing point of water (K)
        real(8), intent(in) :: xxlv           ! latent heat of vaporization
        real(8), intent(in) :: xlf            !latent heat of freezing
        real(8), intent(in) :: xxls           !latent heat of sublimation
        
        real(8), intent(in) :: rhosn  ! bulk density snow
        real(8), intent(in) :: rhoi   ! bulk density ice
        
        real(8), intent(in) :: ac,bc,as,bs,ai,bi,ar,br  !fall speed parameters 
        real(8), intent(in) :: ci,di    !ice mass-diameter relation parameters
        real(8), intent(in) :: cs,ds    !snow mass-diameter relation parameters
        real(8), intent(in) :: cr,dr    !drop mass-diameter relation parameters
        real(8), intent(in) :: f1s,f2s  !ventilation param for snow
        real(8), intent(in) :: Eii      !collection efficiency aggregation of ice
        real(8), intent(in) :: Ecr      !collection efficiency cloud droplets/rain
        real(8), intent(in) :: f1r,f2r  !ventilation param for rain
        real(8), intent(in) :: DCS      !autoconversion size threshold
        real(8), intent(in) :: qsmall   !min mixing ratio 
        real(8), intent(in) :: bimm,aimm !immersion freezing
        real(8), intent(in) :: rhosu     !typical 850mn air density
        real(8), intent(in) :: mi0       ! new crystal mass
        real(8), intent(in) :: rin       ! radius of contact nuclei
        real(8), intent(in) :: pi       ! pi
        
        ! Additional constants to help speed up code
        
        real(8), intent(in) :: cons1
        real(8), intent(in) :: cons4
        real(8), intent(in) :: cons5
        real(8), intent(in) :: cons6
        real(8), intent(in) :: cons7
        real(8), intent(in) :: cons8
        real(8), intent(in) :: cons11
        real(8), intent(in) :: cons13
        real(8), intent(in) :: cons14
        real(8), intent(in) :: cons16
        real(8), intent(in) :: cons17
        real(8), intent(in) :: cons22
        real(8), intent(in) :: cons23
        real(8), intent(in) :: cons24
        real(8), intent(in) :: cons25
        real(8), intent(in) :: cons27
        real(8), intent(in) :: cons28
        
        real(8), intent(in) :: lammini
        real(8), intent(in) :: lammaxi
        real(8), intent(in) :: lamminr
        real(8), intent(in) :: lammaxr
        real(8), intent(in) :: lammins
        real(8), intent(in) :: lammaxs
        
        ! parameters for snow/rain fraction for convective clouds
        real(8), intent(in) :: tmax_fsnow ! max temperature for transition to convective snow
        real(8), intent(in) :: tmin_fsnow ! min temperature for transition to convective snow
        
        !needed for findsp
        real(8), intent(in) :: tt0       ! Freezing temperature
        
        real(8), intent(in) :: csmin,csmax,minrefl,mindbz
        
        real(8), intent(in) :: rhmini     ! Minimum rh for ice cloud fraction > 0.
        
        logical, intent(in) :: use_hetfrz_classnuc ! option to use heterogeneous freezing
        
        ! input arguments
        logical,  intent(in) :: microp_uniform  ! True = configure uniform for sub-columns  False = use w/o sub-columns (standard)
        integer,  intent(in) :: pcols                ! size of column (first) index
        integer,  intent(in) :: pver                 ! number of layers in columns
        integer,  intent(in) :: ncol                 ! number of columns
        integer,  intent(in) :: top_lev              ! top level microphys is applied
        real(8), intent(in) :: deltatin             ! time step (s)
        real(8), intent(in) :: tn(pver)       ! input temperature (K)
        real(8), intent(in) :: qn(pver)       ! input h20 vapor mixing ratio (kg/kg)
        real(8), intent(in) :: relvar(pver)   ! relative variance of cloud water (-)
        real(8), intent(in) :: accre_enhan(pver) ! optional accretion enhancement factor (-)
        
        ! note: all input cloud variables are grid-averaged
        real(8), intent(inout) :: qc(pver)    ! cloud water mixing ratio (kg/kg)
        real(8), intent(inout) :: qi(pver)    ! cloud ice mixing ratio (kg/kg)
        real(8), intent(inout) :: nc(pver)    ! cloud water number conc (1/kg)
        real(8), intent(inout) :: ni(pver)    ! cloud ice number conc (1/kg)
        real(8), intent(in) :: p(pver)        ! air pressure (pa)
        real(8), intent(in) :: pdel(pver)     ! pressure difference across level (pa)
        real(8), intent(in) :: cldn(pver)     ! cloud fraction
        real(8), intent(in) :: icecldf(pver)  ! ice cloud fraction   
        real(8), intent(in) :: liqcldf(pver)  ! liquid cloud fraction
        
        real(8), intent(out) :: rate1ord_cw2pr_st(pver) ! 1st order rate for direct cw to precip conversion 
        ! used for scavenging
        ! Inputs for aerosol activation
        real(8), intent(in) :: naai(pver)      ! ice nulceation number (from microp_aero_ts) 
        real(8), intent(in) :: npccnin(pver)   ! ccn activated number tendency (from microp_aero_ts)
        real(8), intent(in) :: rndst(pver,4)   ! radius of 4 dust bins for contact freezing (from microp_aero_ts)
        real(8), intent(in) :: nacon(pver,4)   ! number in 4 dust bins for contact freezing  (from microp_aero_ts)
        
        ! Used with CARMA cirrus microphysics
        ! (or similar external microphysics model)
        logical,  intent(in) :: do_cldice             ! Prognosing cldice

        real(8), intent(in) :: tnd_qsnow(pver)       
        real(8), intent(in) :: tnd_nsnow(pver)       
        real(8), intent(in) :: re_ice   (pver)       
        real(8), intent(in) :: frzimm   (pver)       
        real(8), intent(in) :: frzcnt   (pver)       
        real(8), intent(in) :: frzdep   (pver)       
        
        ! output arguments
        
        real(8), intent(out) :: tlat(pver)    ! latent heating rate       (W/kg)
        real(8), intent(out) :: qvlat(pver)   ! microphysical tendency qv (1/s)
        real(8), intent(out) :: qctend(pver)  ! microphysical tendency qc (1/s) 
        real(8), intent(out) :: qitend(pver)  ! microphysical tendency qi (1/s)
        real(8), intent(out) :: nctend(pver)  ! microphysical tendency nc (1/(kg*s))
        real(8), intent(out) :: nitend(pver)  ! microphysical tendency ni (1/(kg*s))
        real(8), intent(out) :: effc(pver)    ! droplet effective radius (micron)
        real(8), intent(out) :: effc_fn(pver) ! droplet effective radius, assuming nc = 1.e8 kg-1
        real(8), intent(out) :: effi(pver)    ! cloud ice effective radius (micron)
        real(8), intent(out) :: prect        ! surface precip rate (m/s)
        real(8), intent(out) :: preci        ! cloud ice/snow precip rate (m/s)
        real(8), intent(out) :: nevapr(pver)  ! evaporation rate of rain + snow
        real(8), intent(out) :: evapsnow(pver)! sublimation rate of snow
        real(8), intent(out) :: am_evp_st(pver)! stratiform evaporation area
        real(8), intent(out) :: prain(pver)   ! production of rain + snow
        real(8), intent(out) :: prodsnow(pver)! production of snow
        real(8), intent(out) :: cmeout(pver)  ! evap/sub of cloud
        real(8), intent(out) :: deffi(pver)   ! ice effective diameter for optics (radiation)
        real(8), intent(out) :: pgamrad(pver) ! ice gamma parameter for optics (radiation)
        real(8), intent(out) :: lamcrad(pver) ! slope of droplet distribution for optics (radiation)
        real(8), intent(out) :: qsout(pver)   ! snow mixing ratio (kg/kg)
        real(8), intent(out) :: dsout(pver)   ! snow diameter (m)
        real(8), intent(out) :: rflx(pver+1)  ! grid-box average rain flux (kg m^-2 s^-1)
        real(8), intent(out) :: sflx(pver+1)  ! grid-box average snow flux (kg m^-2 s^-1)
        real(8), intent(out) :: qrout(pver)     ! grid-box average rain mixing ratio (kg/kg)
        real(8), intent(inout) :: reff_rain(pver) ! rain effective radius (micron)
        real(8), intent(inout) :: reff_snow(pver) ! snow effective radius (micron)
        real(8), intent(out) :: qcsevap(pver) ! cloud water evaporation due to sedimentation
        real(8), intent(out) :: qisevap(pver) ! cloud ice sublimation due to sublimation
        real(8), intent(out) :: qvres(pver) ! residual condensation term to ensure RH < 100%
        real(8), intent(out) :: cmeiout(pver) ! grid-mean cloud ice sub/dep
        real(8), intent(out) :: vtrmc(pver) ! mass-weighted cloud water fallspeed
        real(8), intent(out) :: vtrmi(pver) ! mass-weighted cloud ice fallspeed
        real(8), intent(out) :: qcsedten(pver) ! qc sedimentation tendency
        real(8), intent(out) :: qisedten(pver) ! qi sedimentation tendency
        ! microphysical process rates for output (mixing ratio tendencies)
        real(8), intent(out) :: prao(pver) ! accretion of cloud by rain 
        real(8), intent(out) :: prco(pver) ! autoconversion of cloud to rain
        real(8), intent(out) :: mnuccco(pver) ! mixing rat tend due to immersion freezing
        real(8), intent(out) :: mnuccto(pver) ! mixing ratio tend due to contact freezing
        real(8), intent(out) :: msacwio(pver) ! mixing ratio tend due to H-M splintering
        real(8), intent(out) :: psacwso(pver) ! collection of cloud water by snow
        real(8), intent(out) :: bergso(pver) ! bergeron process on snow
        real(8), intent(out) :: bergo(pver) ! bergeron process on cloud ice
        real(8), intent(out) :: melto(pver) ! melting of cloud ice
        real(8), intent(out) :: homoo(pver) ! homogeneos freezign cloud water
        real(8), intent(out) :: qcreso(pver) ! residual cloud condensation due to removal of excess supersat
        real(8), intent(out) :: prcio(pver) ! autoconversion of cloud ice to snow
        real(8), intent(out) :: praio(pver) ! accretion of cloud ice by snow
        real(8), intent(out) :: qireso(pver) ! residual ice deposition due to removal of excess supersat
        real(8), intent(out) :: mnuccro(pver) ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
        real(8), intent(out) :: pracso (pver) ! mixing ratio tendency due to accretion of rain by snow (1/s)
        real(8), intent(out) :: meltsdt(pver) ! latent heating rate due to melting of snow  (W/kg)
        real(8), intent(out) :: frzrdt (pver) ! latent heating rate due to homogeneous freezing of rain (W/kg)
        real(8), intent(out) :: mnuccdo(pver) ! mass tendency from ice nucleation
        real(8), intent(out) :: nrout(pver) ! rain number concentration (1/m3)
        real(8), intent(out) :: nsout(pver) ! snow number concentration (1/m3)
        real(8), intent(out) :: refl(pver)    ! analytic radar reflectivity        
        real(8), intent(out) :: arefl(pver)  !average reflectivity will zero points outside valid range
        real(8), intent(out) :: areflz(pver)  !average reflectivity in z.
        real(8), intent(out) :: frefl(pver)
        real(8), intent(out) :: csrfl(pver)   !cloudsat reflectivity 
        real(8), intent(out) :: acsrfl(pver)  !cloudsat average
        real(8), intent(out) :: fcsrfl(pver)
        real(8), intent(out) :: rercld(pver) ! effective radius calculation for rain + cloud
        real(8), intent(out) :: ncai(pver) ! output number conc of ice nuclei available (1/m3)
        real(8), intent(out) :: ncal(pver) ! output number conc of CCN (1/m3)
        real(8), intent(out) :: qrout2(pver)
        real(8), intent(out) :: qsout2(pver)
        real(8), intent(out) :: nrout2(pver)
        real(8), intent(out) :: nsout2(pver)
        real(8), intent(out) :: drout2(pver) ! mean rain particle diameter (m)
        real(8), intent(out) :: dsout2(pver) ! mean snow particle diameter (m)
        real(8), intent(out) :: freqs(pver)
        real(8), intent(out) :: freqr(pver)
        real(8), intent(out) :: nfice(pver)
        real(8), intent(out) :: prer_evap(pver)
        
        real(8) :: nevapr2(pver)
        
        !character(128),   intent(out) :: errstring       ! Output status (non-blank for error return)
        
        ! Tendencies calculated by external schemes that can replace MG's native
        ! process tendencies.
        
!        ! Used with CARMA cirrus microphysics
!        ! (or similar external microphysics model)
!        real(8), intent(in), pointer :: tnd_qsnow(:,:) ! snow mass tendency (kg/kg/s)
!        real(8), intent(in), pointer :: tnd_nsnow(:,:) ! snow number tendency (#/kg/s)
!        real(8), intent(in), pointer :: re_ice(:,:)    ! ice effective radius (m)
!        
!        ! From external ice nucleation.
!        real(8), intent(in), pointer :: frzimm(:,:) ! Number tendency due to immersion freezing (1/cm3)
!        real(8), intent(in), pointer :: frzcnt(:,:) ! Number tendency due to contact freezing (1/cm3)
!        real(8), intent(in), pointer :: frzdep(:,:) ! Number tendency due to deposition nucleation (1/cm3)
        
        ! local workspace
        ! all units mks unless otherwise stated
        
        ! Additional constants to help speed up code
        real(8) :: cons2
        real(8) :: cons3
        real(8) :: cons9
        real(8) :: cons10
        real(8) :: cons12
        real(8) :: cons15
        real(8) :: cons18
        real(8) :: cons19
        real(8) :: cons20
        
        ! temporary variables for sub-stepping 
        real(8) :: t1(pver)
        real(8) :: q1(pver)
        real(8) :: qc1(pver)
        real(8) :: qi1(pver)
        real(8) :: nc1(pver)
        real(8) :: ni1(pver)
        real(8) :: tlat1(pver)
        real(8) :: qvlat1(pver)
        real(8) :: qctend1(pver)
        real(8) :: qitend1(pver)
        real(8) :: nctend1(pver)
        real(8) :: nitend1(pver)
        real(8) :: prect1
        real(8) :: preci1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        real(8) :: deltat        ! sub-time step (s)
        real(8) :: omsm    ! number near unity for round-off issues
        real(8) :: dto2    ! dt/2 (s)
        real(8) :: mincld  ! minimum allowed cloud fraction
        real(8) :: q(pver) ! water vapor mixing ratio (kg/kg)
        real(8) :: t(pver) ! temperature (K)
        real(8) :: rho(pver) ! air density (kg m-3)
        real(8) :: dv(pver)  ! diffusivity of water vapor in air
        real(8) :: mu(pver)  ! viscocity of air
        real(8) :: sc(pver)  ! schmidt number
        real(8) :: kap(pver) ! thermal conductivity of air
        real(8) :: rhof(pver) ! air density correction factor for fallspeed
        real(8) :: cldmax(pver) ! precip fraction assuming maximum overlap
        real(8) :: cldm(pver)   ! cloud fraction
        real(8) :: icldm(pver)   ! ice cloud fraction
        real(8) :: lcldm(pver)   ! liq cloud fraction
        real(8) :: icwc    ! in cloud water content (liquid+ice)
        real(8) :: calpha  ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: cbeta ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: cbetah ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: cgamma ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: cgamah ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: rcgama ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: cmec1 ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: cmec2 ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: cmec3 ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: cmec4 ! parameter for cond/evap (Zhang et al. 2003)
        real(8) :: qtmp ! dummy qv 
        real(8) :: dum  ! temporary dummy variable
        
        real(8) :: cme(pver)  ! total (liquid+ice) cond/evap rate of cloud
        
        real(8) :: cmei(pver) ! dep/sublimation rate of cloud ice
        real(8) :: cwml(pver) ! cloud water mixing ratio
        real(8) :: cwmi(pver) ! cloud ice mixing ratio
        real(8) :: nnuccd(pver)   ! ice nucleation rate from deposition/cond.-freezing
        real(8) :: mnuccd(pver)   ! mass tendency from ice nucleation
        real(8) :: qcld              ! total cloud water
        real(8) :: lcldn(pver) ! fractional coverage of new liquid cloud
        real(8) :: lcldo(pver) ! fractional coverage of old liquid cloud
        real(8) :: nctend_mixnuc(pver)
        real(8) :: arg ! argument of erfc
        
        ! for calculation of rate1ord_cw2pr_st
        real(8) :: qcsinksum_rate1ord(pver)   ! sum over iterations of cw to precip sink
        real(8) :: qcsum_rate1ord(pver)    ! sum over iterations of cloud water       
        
        real(8) :: alpha
        
        real(8) :: dum1,dum2   !general dummy variables
        
        real(8) :: npccn(pver)     ! droplet activation rate
        real(8) :: qcic(pver) ! in-cloud cloud liquid mixing ratio
        real(8) :: qiic(pver) ! in-cloud cloud ice mixing ratio
        real(8) :: qniic(pver) ! in-precip snow mixing ratio
        real(8) :: qric(pver) ! in-precip rain mixing ratio
        real(8) :: ncic(pver) ! in-cloud droplet number conc
        real(8) :: niic(pver) ! in-cloud cloud ice number conc
        real(8) :: nsic(pver) ! in-precip snow number conc
        real(8) :: nric(pver) ! in-precip rain number conc
        real(8) :: lami(pver) ! slope of cloud ice size distr
        real(8) :: n0i(pver) ! intercept of cloud ice size distr
        real(8) :: lamc(pver) ! slope of cloud liquid size distr
        real(8) :: n0c(pver) ! intercept of cloud liquid size distr
        real(8) :: lams(pver) ! slope of snow size distr
        real(8) :: n0s(pver) ! intercept of snow size distr
        real(8) :: lamr(pver) ! slope of rain size distr
        real(8) :: n0r(pver) ! intercept of rain size distr
        real(8) :: cdist1(pver) ! size distr parameter to calculate droplet freezing
        ! combined size of precip & cloud drops
        real(8) :: arcld(pver) ! averaging control flag
        real(8) :: Actmp  !area cross section of drops
        real(8) :: Artmp  !area cross section of rain
        
        real(8) :: pgam(pver) ! spectral width parameter of droplet size distr
        real(8) :: lammax  ! maximum allowed slope of size distr
        real(8) :: lammin  ! minimum allowed slope of size distr
        real(8) :: nacnt   ! number conc of contact ice nuclei
        real(8) :: mnuccc(pver) ! mixing ratio tendency due to freezing of cloud water
        real(8) :: nnuccc(pver) ! number conc tendency due to freezing of cloud water
        
        real(8) :: mnucct(pver) ! mixing ratio tendency due to contact freezing of cloud water
        real(8) :: nnucct(pver) ! number conc tendency due to contact freezing of cloud water
        real(8) :: msacwi(pver) ! mixing ratio tendency due to HM ice multiplication
        real(8) :: nsacwi(pver) ! number conc tendency due to HM ice multiplication
        
        real(8) :: prc(pver) ! qc tendency due to autoconversion of cloud droplets
        real(8) :: nprc(pver) ! number conc tendency due to autoconversion of cloud droplets
        real(8) :: nprc1(pver) ! qr tendency due to autoconversion of cloud droplets
        real(8) :: nsagg(pver) ! ns tendency due to self-aggregation of snow
        real(8) :: dc0  ! mean size droplet size distr
        real(8) :: ds0  ! mean size snow size distr (area weighted)
        real(8) :: eci  ! collection efficiency for riming of snow by droplets
        real(8) :: psacws(pver) ! mixing rat tendency due to collection of droplets by snow
        real(8) :: npsacws(pver) ! number conc tendency due to collection of droplets by snow
        real(8) :: uni ! number-weighted cloud ice fallspeed
        real(8) :: umi ! mass-weighted cloud ice fallspeed
        real(8) :: uns(pver) ! number-weighted snow fallspeed
        real(8) :: ums(pver) ! mass-weighted snow fallspeed
        real(8) :: unr(pver) ! number-weighted rain fallspeed
        real(8) :: umr(pver) ! mass-weighted rain fallspeed
        real(8) :: unc ! number-weighted cloud droplet fallspeed
        real(8) :: umc ! mass-weighted cloud droplet fallspeed
        real(8) :: pracs(pver) ! mixing rat tendency due to collection of rain by snow
        real(8) :: npracs(pver) ! number conc tendency due to collection of rain by snow
        real(8) :: mnuccr(pver) ! mixing rat tendency due to freezing of rain
        real(8) :: nnuccr(pver) ! number conc tendency due to freezing of rain
        real(8) :: pra(pver) ! mixing rat tendnency due to accretion of droplets by rain
        real(8) :: npra(pver) ! nc tendnency due to accretion of droplets by rain
        real(8) :: nragg(pver) ! nr tendency due to self-collection of rain
        real(8) :: prci(pver) ! mixing rat tendency due to autoconversion of cloud ice to snow
        real(8) :: nprci(pver) ! number conc tendency due to autoconversion of cloud ice to snow
        real(8) :: prai(pver) ! mixing rat tendency due to accretion of cloud ice by snow
        real(8) :: nprai(pver) ! number conc tendency due to accretion of cloud ice by snow
        real(8) :: qvs ! liquid saturation vapor mixing ratio
        real(8) :: qvi ! ice saturation vapor mixing ratio
        real(8) :: dqsdt ! change of sat vapor mixing ratio with temperature
        real(8) :: dqsidt ! change of ice sat vapor mixing ratio with temperature
        real(8) :: ab ! correction factor for rain evap to account for latent heat
        real(8) :: qclr ! water vapor mixing ratio in clear air
        real(8) :: abi ! correction factor for snow sublimation to account for latent heat
        real(8) :: epss ! 1/ sat relaxation timescale for snow
        real(8) :: epsr ! 1/ sat relaxation timescale for rain
        real(8) :: pre(pver) ! rain mixing rat tendency due to evaporation
        real(8) :: prds(pver) ! snow mixing rat tendency due to sublimation
        real(8) :: qce ! dummy qc for conservation check
        real(8) :: qie ! dummy qi for conservation check
        real(8) :: nce ! dummy nc for conservation check
        real(8) :: nie ! dummy ni for conservation check
        real(8) :: ratio ! parameter for conservation check
        real(8) :: dumc(pver) ! dummy in-cloud qc
        real(8) :: dumnc(pver) ! dummy in-cloud nc
        real(8) :: dumi(pver) ! dummy in-cloud qi
        real(8) :: dumni(pver) ! dummy in-cloud ni
        real(8) :: dums(pver) ! dummy in-cloud snow mixing rat
        real(8) :: dumns(pver) ! dummy in-cloud snow number conc
        real(8) :: dumr(pver) ! dummy in-cloud rain mixing rat
        real(8) :: dumnr(pver) ! dummy in-cloud rain number conc
        ! below are parameters for cloud water and cloud ice sedimentation calculations
        real(8) :: fr(pver)
        real(8) :: fnr(pver)
        real(8) :: fc(pver)
        real(8) :: fnc(pver)
        real(8) :: fi(pver)
        real(8) :: fni(pver)
        real(8) :: fs(pver)
        real(8) :: fns(pver)
        real(8) :: faloutr(pver)
        real(8) :: faloutnr(pver)
        real(8) :: faloutc(pver)
        real(8) :: faloutnc(pver)
        real(8) :: falouti(pver)
        real(8) :: faloutni(pver)
        real(8) :: falouts(pver)
        real(8) :: faloutns(pver)
        real(8) :: faltndr
        real(8) :: faltndnr
        real(8) :: faltndc
        real(8) :: faltndnc
        real(8) :: faltndi
        real(8) :: faltndni
        real(8) :: faltnds
        real(8) :: faltndns
        real(8) :: faltndqie
        real(8) :: faltndqce
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(8) :: relhum(pver) ! relative humidity
        real(8) :: csigma ! parameter for cond/evap of cloud water/ice
        real(8) :: rgvm ! max fallspeed for all species
        real(8) :: arn(pver) ! air density corrected rain fallspeed parameter
        real(8) :: asn(pver) ! air density corrected snow fallspeed parameter
        real(8) :: acn(pver) ! air density corrected cloud droplet fallspeed parameter
        real(8) :: ain(pver) ! air density corrected cloud ice fallspeed parameter
        real(8) :: nsubi(pver) ! evaporation of cloud ice number
        real(8) :: nsubc(pver) ! evaporation of droplet number
        real(8) :: nsubs(pver) ! evaporation of snow number
        real(8) :: nsubr(pver) ! evaporation of rain number
        real(8) :: mtime ! factor to account for droplet activation timescale
        real(8) :: dz(pver) ! height difference across model vertical level
        
        
        !! add precip flux variables for sub-stepping
        real(8) :: rflx1(pver+1)
        real(8) :: sflx1(pver+1)
        
        ! returns from function/subroutine calls
        real(8) :: tsp(pver)      ! saturation temp (K)
        real(8) :: qsp(pver)      ! saturation mixing ratio (kg/kg)
        real(8) :: qsphy(pver)      ! saturation mixing ratio (kg/kg): hybrid rh
        real(8) :: qs            ! liquid-ice weighted sat mixing rat (kg/kg)
        real(8) :: es            ! liquid-ice weighted sat vapor press (pa)
        real(8) :: esl(pver)      ! liquid sat vapor pressure (pa)
        real(8) :: esi(pver)      ! ice sat vapor pressure (pa)
        
        ! sum of source/sink terms for diagnostic precip
        
        real(8) :: qnitend(pver) ! snow mixing ratio source/sink term
        real(8) :: nstend(pver)  ! snow number concentration source/sink term
        real(8) :: qrtend(pver) ! rain mixing ratio source/sink term
        real(8) :: nrtend(pver)  ! rain number concentration source/sink term
        real(8) :: qrtot ! vertically-integrated rain mixing rat source/sink term
        real(8) :: nrtot ! vertically-integrated rain number conc source/sink term
        real(8) :: qstot ! vertically-integrated snow mixing rat source/sink term
        real(8) :: nstot ! vertically-integrated snow number conc source/sink term
        
        ! new terms for Bergeron process
        
        real(8) :: dumnnuc ! provisional ice nucleation rate (for calculating bergeron)
        real(8) :: ninew  ! provisional cloud ice number conc (for calculating bergeron)
        real(8) :: qinew ! provisional cloud ice mixing ratio (for calculating bergeron)
        real(8) :: qvl  ! liquid sat mixing ratio   
        real(8) :: epsi ! 1/ sat relaxation timecale for cloud ice
        real(8) :: prd ! provisional deposition rate of cloud ice at water sat 
        real(8) :: berg(pver) ! mixing rat tendency due to bergeron process for cloud ice
        real(8) :: bergs(pver) ! mixing rat tendency due to bergeron process for snow
        
        !bergeron terms
        real(8) :: bergtsf   !bergeron timescale to remove all liquid
        real(8) :: rhin      !modified RH for vapor deposition
        
        ! diagnostic rain/snow for output to history
        ! values are in-precip (local) !!!!
        
        real(8) :: drout(pver)     ! rain diameter (m)
        
        !averageed rain/snow for history
        real(8) :: dumfice
        
        !ice nucleation, droplet activation
        real(8) :: dum2i(pver) ! number conc of ice nuclei available (1/kg)
        real(8) :: dum2l(pver) ! number conc of CCN (1/kg)
        real(8) :: ncmax
        real(8) :: nimax
        
        real(8) :: qcvar     ! 1/relative variance of sub-grid qc
        
        ! loop array variables
        integer i,k,nstep,n, l
        integer ii,kk, m
        
        ! loop variables for sub-step solution
        integer iter,it,ltrue
        
        ! used in contact freezing via dust particles
        real(8)  tcnt, viscosity, mfp
        real(8)  slip1, slip2, slip3, slip4
        !        real(8)  dfaer1, dfaer2, dfaer3, dfaer4
        !        real(8)  nacon1,nacon2,nacon3,nacon4
        real(8)  ndfaer1, ndfaer2, ndfaer3, ndfaer4
        real(8)  nslip1, nslip2, nslip3, nslip4
        
        ! used in ice effective radius
        real(8)  bbi, cci, ak, iciwc, rvi
        
        ! used in Bergeron processe and water vapor deposition
        real(8)  Tk, deles, Aprpr, Bprpr, Cice, qi0, Crate, qidep
        
        ! mean cloud fraction over the time step
        real(8)  cldmw(pver)
        
        ! used in secondary ice production
        real(8) ni_secp
        
        ! variabels to check for RH after rain evap
        
        real(8) :: esn
        real(8) :: qsn
        real(8) :: ttmp
        
        
        
        real(8) :: rainrt(pver)  ! rain rate for reflectivity calculation
        real(8) :: rainrt1(pver)
        real(8) :: tmp
        
        real(8) dmc,ssmc,dstrn  ! variables for modal scheme.
        
        real(8) :: cdnl    ! cloud droplet number limiter
        
        ! heterogeneous freezing
        real(8) :: mnudep(pver) ! mixing ratio tendency due to deposition of water vapor
        real(8) :: nnudep(pver) ! number conc tendency due to deposition of water vapor
        real(8) :: con1 ! work cnstant
        real(8) :: r3lx ! Mean volume radius (m)
        real(8) :: mi0l
        real(8) :: frztmp
        
        logical  :: do_clubb_sgs

        real(8), intent(in) :: svp_tmelt, svp_h2otrip, svp_tboil, svp_ttrice, svp_epsilo
        real(8) :: svp_omeps
        integer(8) :: lg10, power, exp_agent, sqrt_agent, sin_agent, log_agent
        integer(8) :: math_agent(3)
        real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
        real(8) :: tmpn(32)

        integer, intent(in) :: myid
        integer :: kst

        real(8) :: Pgama(8) 
        real(8) :: Qgama(8) 
        real(8) :: Cgama(7) 
        real(8) :: Sgama(6) 
        
        !real(8) :: Para_water(8) = &
        !     (/11.344d0, -3.49149d0, 1013.246d0, -7.90298d0, 5.02808d0, &
        !     1.3816d-7, 8.1328d-3, 0.d0 /)

  !call lwpf_enter_test_unit()
        real(8) :: dp1_81, dp8_794d_5, dp1_496d_6, dp1_5, dp120_0, dp1_414d3
        real(8) :: dp0_54, dp273_15, dp2_47, dp1_15, dp0_0005714, dp0_2714
        real(8) :: dp_1_79, dp1350_0, dp25d_6, dp9_1, dp269_15, dp270_16
        real(8) :: dp298_0, dp0_85, dp1_8d_5, dp28_96d_3, dp8_314409, dp1_1, dp1_257
        real(8) :: dp1_381d_23, dp268_16, dp265_16, dp180_0, dp323_0, dp233_15, dp3_14
        real(8) :: dp1_563d_14, dp238_15 

        cdnl    = 0.d6
        Pgama(:) = &
             (/-1.71618513886549492533811d+0, 2.47656508055759199108314d+1, &
               -3.79804256470945635097577d+2, 6.29331155312818442661052d+2, &
                8.66966202790413211295064d+2,-3.14512729688483675254357d+4, &
               -3.61444134186911729807069d+4, 6.64561438202405440627855d+4 /)
        Qgama(:) = &
             (/-3.08402300119738975254353d+1, 3.15350626979604161529144d+2, &
               -1.01515636749021914166146d+3,-3.10777167157231109440444d+3, &
                2.25381184209801510330112d+4, 4.75584627752788110767815d+3, &
               -1.34659959864969306392456d+5,-1.15132259675553483497211d+5 /)
        Cgama(:) = &
             (/-1.910444077728d-03,          8.4171387781295d-04, &
               -5.952379913043012d-04,       7.93650793500350248d-04, &
               -2.777777777777681622553d-03, 8.333333333333333331554247d-02, &
                5.7083835261d-03 /)
        Sgama(:) = &
             (/0.9189385332046727417803297d0, 3.1415926535897932384626434d0, &
               huge(1.d0), epsilon(1.d0), tiny(1.d0), 171.624d0 /)
        dp1_81 = 1.81d0
        dp8_794d_5 = 8.794d-5
        dp1_496d_6 = 1.496d-6
        dp1_5 = 1.5d0
        dp120_0 = 120.d0
        dp1_414d3 = 1.414d3
        dp0_54 = 0.54d0
        dp273_15 = 273.15d0
        dp2_47 = 2.47d0
        dp1_15 = 1.15d0
        dp0_0005714 = 0.0005714d0
        dp0_2714 = 0.2714d0
        dp_1_79 = -1.79d0
        dp1350_0 = 1350.d0
        dp25d_6 = 25.d-6
        dp9_1 = 9.1d0
        dp269_15 = 269.15d0
        dp270_16 = 270.16d0
        dp298_0 = 298.d0
        dp0_85 = 0.85d0
        dp1_8d_5 = 1.8d-5
        dp28_96d_3 = 28.96d-3
        dp8_314409 = 8.314409d0
        dp1_1 = 1.1d0
        dp1_257 = 1.257d0
        dp1_381d_23 = 1.381d-23
        dp268_16 = 268.16d0
        dp265_16 = 265.16d0
        dp180_0 = 180.d0
        dp323_0 = 323.d0
        dp233_15 = 233.15d0
        dp3_14 = 3.14d0
        dp1_563d_14 = 1.563d-14
        dp238_15 = 238.15d0
        

        !integer(8) :: ni_address

        !ni_address = loc(ni(:))

        !if ( myid == 0 ) then
        !    print*, "ni_address:", ni_address
        !end if

        !qrout2(:)=0.d0
        !qsout2(:)=0.d0
        !nrout2(:)=0.d0
        !nsout2(:)=0.d0
        !drout2(:)=0.d0
        !dsout2(:)=0.d0
        !pgam(:) = 0.d0
        !pra(:) = 0.d0
        !pre(:) = 0.d0
        !lammin=0.d0
        !lammax=0.d0

        svp_omeps = 1.d0 - svp_epsilo
        call math_agent_log10_c(lg10)
        call math_agent_pow_c(power)
        call math_agent_exp_c(exp_agent)
        call math_agent_sqrt_c(sqrt_agent)
        call math_agent_sin_c(sin_agent)
        call math_agent_log_c(log_agent)
        math_agent(1) = sin_agent
        math_agent(2) = lg10
        !math_agent(2) = log_agent
        math_agent(3) = exp_agent
        !para_water(8) = svp_tboil
        !call wv_sat_methods_init(svp_tmelt, svp_h2otrip, svp_tboil, svp_ttrice, svp_epsilo)

        !do k=1,pver
        !    qrout2(k) = qrout2(k)+nacon(k,1)
        !    qsout2(k) = qsout2(k)+nacon(k,2)
        !    drout2(k) = drout2(k)+nacon(k,3)
        !    dsout2(k) = dsout2(k)+nacon(k,4)
        !enddo
        ! initialize  output fields for number conc qand ice nucleation

        ncai(1:pver)=0.d0 
        ncal(1:pver)=0.d0
        
        !Initialize rain size
        rercld(1:pver)=0.d0
        arcld(1:pver)=0.d0
        
        !initialize radiation output variables
        pgamrad(1:pver)=0.d0 ! liquid gamma parameter for optics (radiation)
        lamcrad(1:pver)=0.d0 ! slope of droplet distribution for optics (radiation)
        deffi  (1:pver)=0.d0 ! slope of droplet distribution for optics (radiation)
        !initialize radiation output variables
        !initialize water vapor tendency term output
        qcsevap(1:pver)=0.d0 
        qisevap(1:pver)=0.d0 
        qvres  (1:pver)=0.d0 
        cmeiout (1:pver)=0.d0
        vtrmc (1:pver)=0.d0
        vtrmi (1:pver)=0.d0
        qcsedten (1:pver)=0.d0
        qisedten (1:pver)=0.d0    
        
        prao(1:pver)=0.d0 
        prco(1:pver)=0.d0 
        mnuccco(1:pver)=0.d0 
        mnuccto(1:pver)=0.d0 
        msacwio(1:pver)=0.d0 
        psacwso(1:pver)=0.d0 
        bergso(1:pver)=0.d0 
        bergo(1:pver)=0.d0 
        melto(1:pver)=0.d0 
        homoo(1:pver)=0.d0 
        qcreso(1:pver)=0.d0 
        prcio(1:pver)=0.d0 
        praio(1:pver)=0.d0 
        qireso(1:pver)=0.d0 
        mnuccro(1:pver)=0.d0 
        pracso (1:pver)=0.d0 
        meltsdt(1:pver)=0.d0
        frzrdt (1:pver)=0.d0
        mnuccdo(1:pver)=0.d0
        
        rflx(:)=0.d0
        sflx(:)=0.d0
        effc(:)=0.d0
        effc_fn(:)=0.d0
        effi(:)=0.d0
        
        ! assign variable deltat for sub-stepping...
        deltat=deltatin
        
        ! parameters for scheme
        
        omsm=0.99999d0
        dto2=0.5d0*deltat
        mincld=0.0001d0
        
        ! initialize multi-level fields
        q(1:pver)=qn(1:pver)
        t(1:pver)=tn(1:pver)

        lami(:) = 0.d0
        n0i(:) = 0.d0
        
        ! initialize time-varying parameters
        
        do k=1,pver
              rho(k)=p(k)/(r*t(k))
              !tmp1 = 1.81d0
              call math_agent_2i1o(power, t(k), dp1_81, tmp2)
              dv(k) = dp8_794d_5*tmp2/p(k)
              !tmp1 = 1.5d0
              call math_agent_2i1o(power, t(k), dp1_5, tmp2)
              mu(k) = dp1_496d_6*tmp2/(t(k)+dp120_0) 
              !dv(k) = 8.794d-5*t(k)**1.81d0/p(k)
              !mu(k) = 1.496d-6*t(k)**1.5d0/(t(k)+120.d0) 
              sc(k) = mu(k)/(rho(k)*dv(k))
              !kap(k) = 1.414d3*1.496d-6*t(k)**1.5d0/(t(k)+120.d0) 
              !kap(k) = 1.414d3*1.496d-6*tmp2/(t(k)+120.d0) 
              kap(k) = dp1_414d3*mu(k) 
        
              ! air density adjustment for fallspeed parameters
              ! includes air density correction factor to the
              ! power of 0.54 following Heymsfield and Bansemer 2007
        
              tmp1 = rhosu/rho(k)
              call math_agent_2i1o(power, tmp1, dp0_54, tmp2)
              rhof(k) = tmp2
              !rhof(k)=(rhosu/rho(k))**0.54d0
        
              arn(k)=ar*rhof(k)
              asn(k)=as*rhof(k)
              acn(k)=ac*rhof(k)
              ain(k)=ai*rhof(k)
        
              ! get dz from dp and hydrostatic approx
              ! keep dz positive (define as layer k-1 - layer k)
        
              dz(k)= pdel(k)/(rho(k)*g)

              !! Check
              !qrout2(k) = rho(k)
              !qsout2(k) = dv(k)
              !nrout2(k) = mu(k)
              !nsout2(k) = sc(k)
              !drout2(k) = kap(k)
              !dsout2(k) = dz(k)
        
        end do
        
        !if ( myid .lt. ncol ) then
        ! initialization
        qc(1:top_lev-1) = 0.d0
        qi(1:top_lev-1) = 0.d0
        nc(1:top_lev-1) = 0.d0
        ni(1:top_lev-1) = 0.d0
        t1(1:pver) = t(1:pver)
        q1(1:pver) = q(1:pver)
        qc1(1:pver) = qc(1:pver)
        qi1(1:pver) = qi(1:pver)
        nc1(1:pver) = nc(1:pver)
        ni1(1:pver) = ni(1:pver)
        
        ! initialize tendencies to zero
        tlat1(1:pver)=0.d0
        qvlat1(1:pver)=0.d0
        qctend1(1:pver)=0.d0
        qitend1(1:pver)=0.d0
        nctend1(1:pver)=0.d0
        nitend1(1:pver)=0.d0
        
        ! initialize precip output
        qrout(1:pver)=0.d0
        qsout(1:pver)=0.d0
        nrout(1:pver)=0.d0
        nsout(1:pver)=0.d0
        dsout(1:pver)=0.d0
        
        drout(1:pver)=0.d0
        
        reff_rain(1:pver)=0.d0
        reff_snow(1:pver)=0.d0
        
        ! initialize variables for trop_mozart
        nevapr(1:pver) = 0.d0
        nevapr2(1:pver) = 0.d0
        evapsnow(1:pver) = 0.d0
        prain(1:pver) = 0.d0
        prodsnow(1:pver) = 0.d0
        cmeout(1:pver) = 0.d0
        
        am_evp_st(1:pver) = 0.d0
        
        ! for refl calc
        rainrt1(1:pver) = 0.d0
        
        ! initialize precip fraction and output tendencies
        cldmax(1:pver)=mincld
        
        !initialize aerosol number
        !        naer2(1:ncol,1:pver,:)=0.d0
        dum2l(1:pver)=0.d0
        dum2i(1:pver)=0.d0
        
        ! initialize avg precip rate
        prect1=0.d0
        preci1=0.d0
        
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !Get humidity and saturation vapor pressures
        
        do k=top_lev,pver
        
        
              ! find wet bulk temperature and saturation value for provisional t and q without
              ! condensation
              
              es = svp_water(t(k), lg10, power, svp_tboil)
              qs = svp_to_qsat(es, p(k), svp_epsilo, svp_omeps)

          !tmp1 = 10.d0
          !tmp2 = svp_tboil/t(k)
          !call math_agent_1i1o(lg10, tmp2, tmp3)
          !tmp4 = 11.344d0*(1.d0-t(k)/svp_tboil)
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !tmp4 = (-3.49149d0*(svp_tboil/t(k)-1.d0))
          !call math_agent_2i1o(power, tmp1, tmp4, tmp5)
          !tmp4 = 1013.246d0
          !call math_agent_1i1o(lg10, tmp4, tmp1)
          !tmp4 = -7.90298d0*(svp_tboil/t(k)-1.d0)+ &
          !     5.02808d0*tmp3- &
          !     1.3816d-7*(tmp2-1.d0)+ &
          !     8.1328d-3*(tmp5-1.d0)+ &
          !     tmp1
          !tmp1 = 10.d0
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !es = tmp2 * 100.d0
          !if ( (p(k) - es) <= 0.d0 ) then
          !   qs = 1.0d0
          !else
          !   qs = svp_epsilo*es / (p(k) - svp_omeps*es)
          !end if

        
              ! Prevents negative values.
              if (qs < 0.0d0) then
                 qs = 1.0d0
                 es = p(k)
              end if
        
              esl(k)=svp_water(t(k), lg10, power, svp_tboil)
              esi(k)=svp_ice(t(k), lg10, power, svp_h2otrip)

          !tmp1 = 10.d0
          !tmp2 = svp_tboil/t(k)
          !call math_agent_1i1o(lg10, tmp2, tmp3)
          !tmp4 = 11.344d0*(1.d0-t(k)/svp_tboil)
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !tmp4 = (-3.49149d0*(svp_tboil/t(k)-1.d0))
          !call math_agent_2i1o(power, tmp1, tmp4, tmp5)
          !tmp4 = 1013.246d0
          !call math_agent_1i1o(lg10, tmp4, tmp1)
          !tmp4 = -7.90298d0*(svp_tboil/t(k)-1.d0)+ &
          !     5.02808d0*tmp3- &
          !     1.3816d-7*(tmp2-1.d0)+ &
          !     8.1328d-3*(tmp5-1.d0)+ &
          !     tmp1
          !tmp1 = 10.d0
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !esl(k) = tmp2 * 100.d0

          !tmp1 = svp_h2otrip/t(k)
          !call math_agent_1i1o(lg10, tmp1, tmp2)
          !tmp1 = 6.1071d0
          !call math_agent_1i1o(lg10, tmp1, tmp3)
          !tmp1 = -9.09718d0*(svp_h2otrip/t(k)-1.d0)-3.56654d0* &
          !     tmp2+0.876793d0*(1.d0-t(k)/svp_h2otrip)+ &
          !     tmp3
          !tmp2 = 10.d0
          !call math_agent_2i1o(power, tmp2, tmp1, tmp3)
          !esi(k) = tmp3*100.d0
        

              ! hm fix, make sure when above freezing that esi=esl, not active yet
              if (t(k).gt.tmelt)esi(k)=esl(k)
        
              relhum(k)=q(k)/qs
        
              ! get cloud fraction, check for minimum
        
              cldm(k)=max(cldn(k),mincld)
              cldmw(k)=max(cldn(k),mincld)
        
              icldm(k)=max(icecldf(k),mincld)
              lcldm(k)=max(liqcldf(k),mincld)

        
              ! subcolumns, set cloud fraction variables to one
              ! if cloud water or ice is present, if not present
              ! set to mincld (mincld used instead of zero, to prevent
              ! possible division by zero errors
        
              if (microp_uniform) then
        
                 cldm(k)=mincld
                 cldmw(k)=mincld
                 icldm(k)=mincld
                 lcldm(k)=mincld
        
                 if (qc(k).ge.qsmall) then
                    lcldm(k)=1.d0           
                    cldm(k)=1.d0
                    cldmw(k)=1.d0
                 end if
        
                 if (qi(k).ge.qsmall) then             
                    cldm(k)=1.d0
                    icldm(k)=1.d0
                 end if
        
              end if               ! sub-columns
        
              ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
        
              nfice(k)=0.d0
              dumfice=qc(k)+qi(k)
              if (dumfice.gt.qsmall .and. qi(k).gt.qsmall) then
                 nfice(k)=qi(k)/dumfice
              endif
        
              if (do_cldice .and. (t(k).lt.tmelt - 5.d0)) then
        
                 ! if aerosols interact with ice set number of activated ice nuclei
                 dum2=naai(k)
        
                 dumnnuc=(dum2-ni(k)/icldm(k))/deltat*icldm(k)
                 dumnnuc=max(dumnnuc,0.d0)
                 ! get provisional ni and qi after nucleation in order to calculate
                 ! Bergeron process below
                 ninew=ni(k)+dumnnuc*deltat
                 qinew=qi(k)+dumnnuc*deltat*mi0
        
                 !T>268
              else
                 ninew=ni(k)
                 qinew=qi(k)
              end if
        
              ! Initialize CME components
        
              cme(k) = 0.d0
              cmei(k)=0.d0
        
        
              !-------------------------------------------------------------------
              !Bergeron process
        
              ! make sure to initialize bergeron process to zero
              berg(k)=0.d0
              prd = 0.d0
        
              !condensation loop.
        
              ! get in-cloud qi and ni after nucleation
              if (icldm(k) .gt. 0.d0) then 
                 qiic(k)=qinew/icldm(k)
                 niic(k)=ninew/icldm(k)
              else
                 qiic(k)=0.d0
                 niic(k)=0.d0
              endif
        
              !! Check
              !qrout2(k) = relhum(k)
              !qsout2(k) = qiic(k)
              !nrout2(k) = niic(k)
              !nsout2(k) = nfice(k)
              !drout2(k) = qinew 
              !dsout2(k) = ninew

              !if T < 0 C then bergeron.
              if (do_cldice .and. (t(k).lt.dp273_15)) then
        
                 !if ice exists
                 if (qi(k).gt.qsmall) then
        
                    bergtsf = 0.d0 ! bergeron time scale (fraction of timestep)
        
                    qvi = svp_to_qsat(esi(k), p(k), svp_epsilo, svp_omeps)
                    qvl = svp_to_qsat(esl(k), p(k), svp_epsilo, svp_omeps)

          !if ( (p(k) - esi(k)) <= 0.d0 ) then
          !   qvi = 1.0d0
          !else
          !   qvi = svp_epsilo*esi(k) / (p(k) - svp_omeps*esi(k))
          !end if

          !if ( (p(k) - esl(k)) <= 0.d0 ) then
          !   qvl = 1.0d0
          !else
          !   qvl = svp_epsilo*esl(k) / (p(k) - svp_omeps*esl(k))
          !end if

      !qrout2(k) = qrout2(k)+esi(k)
      !qsout2(k) = qsout2(k)+esl(k)
      !nrout2(k) = nrout2(k)+qvi
      !nsout2(k) = nsout2(k)+qvl
      !drout2(k) = qinew
      !dsout2(k) = ninew

                    !dqsidt =  xxls*qvi/(rv*t(k)**2)
                    dqsidt =  xxls*qvi/(rv*t(k)*t(k))
                    abi = 1.d0+dqsidt*xxls/cpp
        
                    ! get ice size distribution parameters
        
                    if (qiic(k).ge.qsmall) then
                       tmp1 = 1.d0/di
                       tmp2 = (cons1*ci* &
                            niic(k)/qiic(k))
                       call math_agent_2i1o(power, tmp2, tmp1, lami(k))
                       !lami(k) = (cons1*ci* &
                       !     niic(k)/qiic(k))**(1.d0/di)
                       n0i(k) = niic(k)*lami(k)

        
                       ! check for slope
                       ! adjust vars
                       if (lami(k).lt.lammini) then
        
                          lami(k) = lammini
                          tmp1 = di+1.d0
                          call math_agent_2i1o(power, lami(k), tmp1, tmp2)
                          n0i(k) = tmp2*qiic(k)/(ci*cons1)
                          !n0i(k) = lami(k)**(di+1.d0)*qiic(k)/(ci*cons1)
                       else if (lami(k).gt.lammaxi) then
                          lami(k) = lammaxi
                          tmp1 = di+1.d0
                          call math_agent_2i1o(power, lami(k), tmp1, tmp2)
                          n0i(k) = tmp2*qiic(k)/(ci*cons1)
                          !n0i(k) = lami(k)**(di+1.d0)*qiic(k)/(ci*cons1)
                       end if
        
                       epsi = 2.d0*pi*n0i(k)*rho(k)*Dv(k)/(lami(k)*lami(k))
        
                       !if liquid exists  
                       if (qc(k).gt. qsmall) then 
        
                          !begin bergeron process
                          !     do bergeron (vapor deposition with RHw=1)
                          !     code to find berg (a rate) goes here
        
                          ! calculate Bergeron process
        
                          prd = epsi*(qvl-qvi)/abi
        
                       else
                          prd = 0.d0
                       end if
        
                       ! multiply by cloud fraction
        
                       prd = prd*min(icldm(k),lcldm(k))
        
                       !     transfer of existing cloud liquid to ice
        
                       berg(k)=max(0.d0,prd)
        
                    end if  !end liquid exists bergeron
        
                    if (berg(k).gt.0.d0) then
                       bergtsf=max(0.d0,(qc(k)/berg(k))/deltat) 
        
                       if(bergtsf.lt.1.d0) berg(k) = max(0.d0,qc(k)/deltat)
        
                    endif

      !qrout2(k) = qrout2(k)+bergtsf
      !qsout2(k) = qsout2(k)+qiic(k)
      !nrout2(k) = nrout2(k)+qc(k)
      !nsout2(k) = nsout2(k)+relhum(k)
      !drout2(k) = drout2(k)+icldm(k)
      !dsout2(k) = dsout2(k)+lcldm(k)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
                    if (bergtsf.lt.1.d0.or.icldm(k).gt.lcldm(k)) then
        
                       if (qiic(k).ge.qsmall) then
        
                          ! first case is for case when liquid water is present, but is completely depleted 
                          ! in time step, i.e., bergrsf > 0 but < 1
        
                          if (qc(k).ge.qsmall) then
                             rhin  = (1.0d0 + relhum(k)) / 2.d0
                             if ((rhin*esl(k)/esi(k)) > 1.d0) then
                                prd = epsi*(rhin*qvl-qvi)/abi
        
                                ! multiply by cloud fraction assuming liquid/ice maximum overlap
                                prd = prd*min(icldm(k),lcldm(k))
        
                                ! add to cmei
                                cmei(k) = cmei(k) + (prd * (1.d0- bergtsf))

        
                             end if ! rhin 
                          end if ! qc > qsmall
        
                          ! second case is for pure ice cloud, either no liquid, or icldm > lcldm
        
                          if (qc(k).lt.qsmall.or.icldm(k).gt.lcldm(k)) then
        
                             ! note: for case of no liquid, need to set liquid cloud fraction to zero
                             ! store liquid cloud fraction in 'dum'
        
                             if (qc(k).lt.qsmall) then 
                                dum=0.d0 
                             else
                                dum=lcldm(k)
                             end if
        
                             ! set RH to grid-mean value for pure ice cloud
                             rhin = relhum(k)
        
                             if ((rhin*esl(k)/esi(k)) > 1.d0) then
        
                                prd = epsi*(rhin*qvl-qvi)/abi
        
                                ! multiply by relevant cloud fraction for pure ice cloud
                                ! assuming maximum overlap of liquid/ice
                                prd = prd*max((icldm(k)-dum),0.d0)
                                cmei(k) = cmei(k) + prd
        
                        !qrout2(k) = qrout2(k) + epsi
                        !qsout2(k) = qsout2(k) + rhin*qvl
                        !nrout2(k) = nrout2(k) + qvi
                        !nsout2(k) = nsout2(k) + rhin*qvl-qvi
                        !dsout2(k) = dsout2(k) + abi

                             end if ! rhin
                          end if ! qc or icldm > lcldm
                       end if ! qiic
                    end if ! bergtsf or icldm > lcldm
        
                    !     if deposition, it should not reduce grid mean rhi below 1.0
                    if(cmei(k) > 0.0d0 .and. (relhum(k)*esl(k)/esi(k)) > 1.d0 ) &
                         cmei(k)=min(cmei(k),(q(k)-qs*esi(k)/esl(k))/abi/deltat)
        
                 end if            !end ice exists loop
                 !this ends temperature < 0. loop
        
                 !-------------------------------------------------------------------
              end if  ! 
              !..............................................................
        
              ! evaporation should not exceed available water
        
              if ((-berg(k)).lt.-qc(k)/deltat) berg(k) = max(qc(k)/deltat,0.d0)
        
              !sublimation process...
              if (do_cldice .and. ((relhum(k)*esl(k)/esi(k)).lt.1.d0 .and. qiic(k).ge.qsmall )) then
        
                 qvi = svp_to_qsat(esi(k), p(k), svp_epsilo, svp_omeps)
                 qvl = svp_to_qsat(esl(k), p(k), svp_epsilo, svp_omeps)

          !if ( (p(k) - esi(k)) <= 0.d0 ) then
          !   qvi = 1.0d0
          !else
          !   qvi = svp_epsilo*esi(k) / (p(k) - svp_omeps*esi(k))
          !end if

          !if ( (p(k) - esl(k)) <= 0.d0 ) then
          !   qvl = 1.0d0
          !else
          !   qvl = svp_epsilo*esl(k) / (p(k) - svp_omeps*esl(k))
          !end if

                 !dqsidt =  xxls*qvi/(rv*t(k)**2)
                 dqsidt =  xxls*qvi/(rv*t(k)*t(k))
                 abi = 1.d0+dqsidt*xxls/cpp
        
                 ! get ice size distribution parameters
        
                 tmp1 = 1.d0/di
                 tmp2 = (cons1*ci* &
                 niic(k)/qiic(k))
                 call math_agent_2i1o(power, tmp2, tmp1, lami(k))
                 !lami(k) = (cons1*ci* &
                 !     niic(k)/qiic(k))**(1.d0/di)
                 n0i(k) = niic(k)*lami(k)
        
                 ! check for slope
                 ! adjust vars
                 if (lami(k).lt.lammini) then
        
                    lami(k) = lammini
                    tmp1 = di+1.d0
                    call math_agent_2i1o(power, lami(k), tmp1, tmp2)
                    n0i(k) = tmp2*qiic(k)/(ci*cons1)
                    !n0i(k) = lami(k)**(di+1.d0)*qiic(k)/(ci*cons1)
                 else if (lami(k).gt.lammaxi) then
                    lami(k) = lammaxi
                    tmp1 = di+1.d0
                    call math_agent_2i1o(power, lami(k), tmp1, tmp2)
                    n0i(k) = tmp2*qiic(k)/(ci*cons1)
                    !n0i(k) = lami(k)**(di+1.d0)*qiic(k)/(ci*cons1)
                 end if
        
                 epsi = 2.d0*pi*n0i(k)*rho(k)*Dv(k)/(lami(k)*lami(k))
        
                 ! modify for ice fraction below
                 prd = epsi*(relhum(k)*qvl-qvi)/abi * icldm(k)
                 cmei(k)=min(prd,0.d0)
        
              endif
        
              ! sublimation should not exceed available ice
              if (cmei(k).lt.-qi(k)/deltat) cmei(k)=-qi(k)/deltat
        
              ! sublimation should not increase grid mean rhi above 1.0 
              if(cmei(k) < 0.0d0 .and. (relhum(k)*esl(k)/esi(k)) < 1.d0 ) &
                   cmei(k)=min(0.d0,max(cmei(k),(q(k)-qs*esi(k)/esl(k))/abi/deltat))
        
              ! limit cmei due for roundoff error
        
              cmei(k)=cmei(k)*omsm
        
              ! conditional for ice nucleation 
              if (do_cldice .and. (t(k).lt.(tmelt - 5.d0))) then 
        
                 ! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
                 ! ice nucleation rate (dum2) has already been calculated and read in (naai)
        
                 dum2i(k)=naai(k)
              else
                 dum2i(k)=0.d0
              end if

        end do ! k loop
        !end if
        
        !do k=top_lev,pver
        !      qrout2(k) = qiic(k)
        !      qsout2(k) = niic(k)
        !      nrout2(k) = nfice(k)
        !      nsout2(k) = berg(k)
        !      drout2(k) = cmei(k)
        !      dsout2(k) = dum2i(k)
        !end do
        
        !! initialize sub-step precip flux variables
           !! flux is zero at top interface, so these should stay as 0.
           rflx1(1)=0.d0
           sflx1(1)=0.d0
           do k=top_lev,pver
        
              ! initialize normal and sub-step precip flux variables
              rflx1(k+1)=0.d0
              sflx1(k+1)=0.d0
        end do ! k loop
        !! initialize final precip flux variables.
           !! flux is zero at top interface, so these should stay as 0.
           rflx(1)=0.d0
           sflx(1)=0.d0
           do k=top_lev,pver
              ! initialize normal and sub-step precip flux variables
              rflx(k+1)=0.d0
              sflx(k+1)=0.d0
        end do ! k loop
        
           ltrue=0
           do k=top_lev,pver
              ! skip microphysical calculations if no cloud water
        
              if (qc(k).ge.qsmall.or.qi(k).ge.qsmall.or.cmei(k).ge.qsmall) ltrue=1
           end do

           call ltrue_syn(ltrue, myid, pcols)
        !do k=1,pver
        !      qrout2(k) = qc(k)
        !      qsout2(k) = qi(k)
        !      nrout2(k) = cmei(k)
        !      nsout2(k) = qsmall
        !      drout2(k) = ltrue
        !      dsout2(k) = deltat
        !end do
        
        ! assign number of sub-steps to iter
        ! use 2 sub-steps, following tests described in MG2008
        iter = 2
        
        ! get sub-step time step
        deltat=deltat/real(iter)
        
        ! since activation/nucleation processes are fast, need to take into account
        ! factor mtime = mixing timescale in cloud / model time step
        ! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
        ! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk
        
        !        note: mtime for bulk aerosols was set to: mtime=deltat/1200.d0
        
        mtime=1.d0
        rate1ord_cw2pr_st(:)=0.d0 ! rce 2010/05/01
        
        !!!! skip calculations if no cloud water
           if (ltrue.eq.0) then
              tlat(1:pver)=0.d0
              qvlat(1:pver)=0.d0
              qctend(1:pver)=0.d0
              qitend(1:pver)=0.d0
              qnitend(1:pver)=0.d0
              qrtend(1:pver)=0.d0
              nctend(1:pver)=0.d0
              nitend(1:pver)=0.d0
              nrtend(1:pver)=0.d0
              nstend(1:pver)=0.d0
              prect=0.d0
              preci=0.d0
              qniic(1:pver)=0.d0
              qric(1:pver)=0.d0
              nsic(1:pver)=0.d0
              nric(1:pver)=0.d0
              rainrt(1:pver)=0.d0
              goto 300
           end if
        
           qcsinksum_rate1ord(1:pver)=0.d0 
           qcsum_rate1ord(1:pver)=0.d0 
        
            !do it=1,iter
            !do k=1,pver
            !      qrout2(k) = qrout2(k) + rflx1(k)
            !      qsout2(k) = qsout2(k) + sflx1(k)
            !      nrout2(k) = nrout2(k) + rflx(k)
            !      nsout2(k) = nsout2(k) + sflx(k)
            !      drout2(k) = drout2(k) + qc(k)
            !      dsout2(k) = dsout2(k) + qi(k)
            !end do
            !end do
        
        !!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !.....................................................................................................
           do it=1,iter
        
              ! initialize sub-step microphysical tendencies
        
              tlat(1:pver)=0.d0
              qvlat(1:pver)=0.d0
              qctend(1:pver)=0.d0
              qitend(1:pver)=0.d0
              qnitend(1:pver)=0.d0
              qrtend(1:pver)=0.d0
              nctend(1:pver)=0.d0
              nitend(1:pver)=0.d0
              nrtend(1:pver)=0.d0
              nstend(1:pver)=0.d0
        
              ! initialize diagnostic precipitation to zero
        
              qniic(1:pver)=0.d0
              qric(1:pver)=0.d0
              nsic(1:pver)=0.d0
              nric(1:pver)=0.d0
        
              rainrt(1:pver)=0.d0
        
        
              ! begin new i,k loop, calculate new cldmax after adjustment to cldm above
        
              ! initialize vertically-integrated rain and snow tendencies
        
              qrtot = 0.d0
              nrtot = 0.d0
              qstot = 0.d0
              nstot = 0.d0
        
              ! initialize precip at surface
        
              prect=0.d0
              preci=0.d0

              if ( myid >= 32 ) then
                  call get_k_8(tmpn)
                  qrtot = tmpn(25)
                  nrtot = tmpn(26)
                  qstot = tmpn(27)
                  nstot = tmpn(28)
                  prect = tmpn(29)
                  preci = tmpn(30)
              endif
        
              do k=top_lev,pver

                 !if ( k .eq. top_lev .and. myid >= 32 ) then
                 !    call get_k_1(tmpn, 30, myid)
                 !    qrtot = tmpn(25)
                 !    nrtot = tmpn(26)
                 !    qstot = tmpn(27)
                 !    nstot = tmpn(28)
                 !    prect = tmpn(29)
                 !    preci = tmpn(30)
                 !endif


                 qcvar=relvar(k)
                 cons2=shr_spfn_gamma_nonintrinsic_r8(qcvar+dp2_47,math_agent,Pgama,Qgama,Cgama,Sgama)
                 cons3=shr_spfn_gamma_nonintrinsic_r8(qcvar,math_agent,Pgama,Qgama,Cgama,Sgama)
                 cons9=shr_spfn_gamma_nonintrinsic_r8(qcvar+2.d0,math_agent,Pgama,Qgama,Cgama,Sgama)
                 cons10=shr_spfn_gamma_nonintrinsic_r8(qcvar+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)
                 cons12=shr_spfn_gamma_nonintrinsic_r8(qcvar+dp1_15,math_agent,Pgama,Qgama,Cgama,Sgama) 
                 cons15=shr_spfn_gamma_nonintrinsic_r8(qcvar+bc/3.d0,math_agent,Pgama,Qgama,Cgama,Sgama)
                 tmp1 = dp2_47
                 call math_agent_2i1o(power, qcvar, tmp1, cons18)
                 cons19=qcvar*qcvar
                 tmp1 = dp1_15
                 call math_agent_2i1o(power, qcvar, tmp1, cons20)
                 !cons18=qcvar**2.47d0
                 !cons19=qcvar**2
                 !cons20=qcvar**1.15d0
        
        !if ( it .eq. 1 ) then
        !    qrout2(k) = qrout2(k) + qcvar
        !    qsout2(k) = qsout2(k) + cons2
        !    nrout2(k) = nrout2(k) + cons15
        !    nsout2(k) = nsout2(k) + cons18
        !    drout2(k) = drout2(k) + qc(k)
        !    dsout2(k) = dsout2(k) + qi(k)
        !end if 

                 ! set cwml and cwmi to current qc and qi
        
                 cwml(k)=qc(k)
                 cwmi(k)=qi(k)
        
                 ! initialize precip fallspeeds to zero
        
                 ums(k)=0.d0 
                 uns(k)=0.d0 
                 umr(k)=0.d0 
                 unr(k)=0.d0
        
                 ! calculate precip fraction based on maximum overlap assumption
        
                 ! for sub-columns cldm has already been set to 1 if cloud
                 ! water or ice is present, so cldmax will be correctly set below
                 ! and nothing extra needs to be done here
        
                 !if ( k .eq. 1 ) then
                     !call get_k_1(tmpn, 24, myid)
                 !end if

                 if ( myid < 32 ) then 
                    if (k.eq.top_lev) then
                      cldmax(k)=cldm(k)
                    else
                       ! if rain or snow mix ratio is smaller than
                       ! threshold, then set cldmax to cloud fraction at current level
        
                       if (do_clubb_sgs) then
                          if (qc(k).ge.qsmall.or.qi(k).ge.qsmall) then
                             cldmax(k)=cldm(k)
                          else
                             cldmax(k)=cldmax(k-1)
                          end if
                       else
        
                          if (qric(k-1).ge.qsmall.or.qniic(k-1).ge.qsmall) then
                             cldmax(k)=max(cldmax(k-1),cldm(k))
                          else
                             cldmax(k)=cldm(k)
                          end if
                       endif
                    end if
                else
                    if (k.eq.top_lev) then
                       if (do_clubb_sgs) then
                          if (qc(k).ge.qsmall.or.qi(k).ge.qsmall) then
                             cldmax(k)=cldm(k)
                          else
                             cldmax(k)=tmpn(5)
                          end if
                       else
        
                          if (tmpn(15).ge.qsmall.or.tmpn(3).ge.qsmall) then
                             cldmax(k)=max(tmpn(5),cldm(k))
                          else
                             cldmax(k)=cldm(k)
                          end if
                       endif
                    else
                       ! if rain or snow mix ratio is smaller than
                       ! threshold, then set cldmax to cloud fraction at current level
        
                       if (do_clubb_sgs) then
                          if (qc(k).ge.qsmall.or.qi(k).ge.qsmall) then
                             cldmax(k)=cldm(k)
                          else
                             cldmax(k)=cldmax(k-1)
                          end if
                       else
        
                          if (qric(k-1).ge.qsmall.or.qniic(k-1).ge.qsmall) then
                             cldmax(k)=max(cldmax(k-1),cldm(k))
                          else
                             cldmax(k)=cldm(k)
                          end if
                       endif
                    end if
                endif



         !if ( it .eq. 1 .and. k > 1 ) then
         !    qrout2(k) = qrout2(k) + cldmax(k-1) 
         !    qsout2(k) = qsout2(k) + cldmax(k)
         !    nrout2(k) = nrout2(k) + qniic(k-1)
         !    nsout2(k) = nsout2(k) + qi(k-1)
         !    drout2(k) = drout2(k) + cldm(k)
         !    dsout2(k) = dsout2(k) + qric(k-1)
         !endif
        
                 ! decrease in number concentration due to sublimation/evap
                 ! divide by cloud fraction to get in-cloud decrease
                 ! don't reduce Nc due to bergeron process
        
                 if (cmei(k) < 0.d0 .and. qi(k) > qsmall .and. cldm(k) > mincld) then
                    nsubi(k)=cmei(k)/qi(k)*ni(k)/cldm(k)
                 else
                    nsubi(k)=0.d0
                 end if
                 nsubc(k)=0.d0
        
        
                 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
                 ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
        
                 if (do_cldice .and. dum2i(k).gt.0.d0.and.t(k).lt.(tmelt - 5.d0).and. &
                      relhum(k)*esl(k)/esi(k).gt. rhmini+0.05d0) then
        
                    !if NCAI > 0. then set numice = ncai (as before)
                    !note: this is gridbox averaged
        
                    nnuccd(k)=(dum2i(k)-ni(k)/icldm(k))/deltat*icldm(k)
                    nnuccd(k)=max(nnuccd(k),0.d0)
                    nimax = dum2i(k)*icldm(k)
        
                    !Calc mass of new particles using new crystal mass...
                    !also this will be multiplied by mtime as nnuccd is...
        
                    mnuccd(k) = nnuccd(k) * mi0
        
                    !  add mnuccd to cmei....
                    cmei(k)= cmei(k) + mnuccd(k) * mtime
        
                    !  limit cmei
        
                    qvi = svp_to_qsat(esi(k), p(k), svp_epsilo, svp_omeps)

          !if ( (p(k) - esi(k)) <= 0.d0 ) then
          !   qvi = 1.0d0
          !else
          !   qvi = svp_epsilo*esi(k) / (p(k) - svp_omeps*esi(k))
          !end if

                    !dqsidt =  xxls*qvi/(rv*t(k)**2)
                    dqsidt =  xxls*qvi/(rv*t(k)*t(k))
                    abi = 1.d0+dqsidt*xxls/cpp
                    cmei(k)=min(cmei(k),(q(k)-qvi)/abi/deltat)
        
                    ! limit for roundoff error
                    cmei(k)=cmei(k)*omsm
        
                 else
                    nnuccd(k)=0.d0
                    nimax = 0.d0
                    mnuccd(k) = 0.d0
                 end if
        
                 !c............................................................................
                 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
                 ! for microphysical process calculations
                 ! units are kg/kg for mixing ratio, 1/kg for number conc
        
                 ! limit in-cloud values to 0.005 kg/kg
        
                 qcic(k)=min(cwml(k)/lcldm(k),5.d-3)
                 qiic(k)=min(cwmi(k)/icldm(k),5.d-3)
                 ncic(k)=max(nc(k)/lcldm(k),0.d0)
                 niic(k)=max(ni(k)/icldm(k),0.d0)
        
                 if (qc(k) - berg(k)*deltat.lt.qsmall) then
                    qcic(k)=0.d0
                    ncic(k)=0.d0
                    if (qc(k)-berg(k)*deltat.lt.0.d0) then
                       berg(k)=qc(k)/deltat*omsm
                    end if
                 end if
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+qcic(k)
         !    qsout2(k) = qsout2(k)+qiic(k)
         !    nrout2(k) = nrout2(k)+ncic(k)
         !    nsout2(k) = nsout2(k)+niic(k)
         !    drout2(k) = drout2(k)+qc(k)
         !    dsout2(k) = dsout2(k)+berg(k)
         !endif
                 if (do_cldice .and. qi(k)+(cmei(k)+berg(k))*deltat.lt.qsmall) then
                    qiic(k)=0.d0
                    niic(k)=0.d0
                    if (qi(k)+(cmei(k)+berg(k))*deltat.lt.0.d0) then
                       cmei(k)=(-qi(k)/deltat-berg(k))*omsm
                    end if
                 end if
        
                 ! add to cme output
        
                 cmeout(k) = cmeout(k)+cmei(k)

        !if ( it .eq. 1 ) then
        !    qrout2(k) = qrout2(k) + qcic(k)  
        !    qsout2(k) = qsout2(k) + ncic(k)  
        !    nrout2(k) = nrout2(k) + berg(k)  
        !    nsout2(k) = nsout2(k) + niic(k)  
        !    drout2(k) = drout2(k) + qiic(k)  
        !    dsout2(k) = dsout2(k) + cmeout(k)
        !end if 
        
                 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 ! droplet activation
                 ! calculate potential for droplet activation if cloud water is present
                 ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
                 ! number tendency (npccnin) is read in from companion routine
        
                 ! assume aerosols already activated are equal to number of existing droplets for simplicity
                 ! multiply by cloud fraction to obtain grid-average tendency
        
                 if (qcic(k).ge.qsmall) then   
                    npccn(k) = max(0.d0,npccnin(k))  
                    dum2l(k)=(nc(k)+npccn(k)*deltat)/lcldm(k)
                    dum2l(k)=max(dum2l(k),cdnl/rho(k)) ! sghan minimum in #/cm3  
                    ncmax = dum2l(k)*lcldm(k)
                 else
                    npccn(k)=0.d0
                    dum2l(k)=0.d0
                    ncmax = 0.d0
                 end if
        
                 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 ! get size distribution parameters based on in-cloud cloud water/ice 
                 ! these calculations also ensure consistency between number and mixing ratio
                 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
                 !......................................................................
                 ! cloud ice
        
                 if (qiic(k).ge.qsmall) then
        
                    ! add upper limit to in-cloud number concentration to prevent numerical error
                    niic(k)=min(niic(k),qiic(k)*1.d20)
        
                    tmp1 = 1.d0/di
                    tmp2 = (cons1*ci*niic(k)/qiic(k))
                    call math_agent_2i1o(power, tmp2, tmp1, lami(k))
                    !lami(k) = (cons1*ci*niic(k)/qiic(k))**(1.d0/di)
                    n0i(k) = niic(k)*lami(k)
        
                    ! check for slope
                    ! adjust vars
        
                    if (lami(k).lt.lammini) then
        
                       lami(k) = lammini
                       tmp1 = di+1.d0
                       call math_agent_2i1o(power, lami(k), tmp1, tmp2)
                       n0i(k) = tmp2*qiic(k)/(ci*cons1)
                       !n0i(k) = lami(k)**(di+1.d0)*qiic(k)/(ci*cons1)
                       niic(k) = n0i(k)/lami(k)
                    else if (lami(k).gt.lammaxi) then
                       lami(k) = lammaxi
                       tmp1 = di+1.d0
                       call math_agent_2i1o(power, lami(k), tmp1, tmp2)
                       n0i(k) = tmp2*qiic(k)/(ci*cons1)
                       !n0i(k) = lami(k)**(di+1.d0)*qiic(k)/(ci*cons1)
                       niic(k) = n0i(k)/lami(k)
                    end if
        
                 else
                    lami(k) = 0.d0
                    n0i(k) = 0.d0
                 end if
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+npccn(k)
         !    qsout2(k) = qsout2(k)+ncmax
         !    nrout2(k) = nrout2(k)+dum2l(k)
         !    nsout2(k) = nsout2(k)+lami(k)
         !    drout2(k) = drout2(k)+qcic(k)
         !    dsout2(k) = dsout2(k)+qsmall
         !endif
                 pgam(k) = 0.d0
                 lammin = 0.d0
                 lammax = 0.d0
                 if (qcic(k).ge.qsmall) then
        
                    ! add upper limit to in-cloud number concentration to prevent numerical error
                    ncic(k)=min(ncic(k),qcic(k)*1.d20)
        
                    ncic(k)=max(ncic(k),cdnl/rho(k)) ! sghan minimum in #/cm  
        
                    ! get pgam from fit to observations of martin et al. 1994
        
                    pgam(k)=dp0_0005714*(ncic(k)/1.d6*rho(k))+dp0_2714
                    pgam(k)=1.d0/(pgam(k)*pgam(k))-1.d0
                    pgam(k)=max(pgam(k),2.d0)
                    pgam(k)=min(pgam(k),15.d0)

                    ! calculate lamc
        
                    tmp1 = (pi/6.d0*rhow*ncic(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/ &
                         (qcic(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)))
                    tmp2 = 1.d0/3.d0
                    call math_agent_2i1o(power, tmp1, tmp2, lamc(k))
                    !lamc(k) = (pi/6.d0*rhow*ncic(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0)/ &
                    !     (qcic(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0)))**(1.d0/3.d0)
        
                    ! lammin, 50 micron diameter max mean size
        
                    lammin = (pgam(k)+1.d0)/50.d-6
                    lammax = (pgam(k)+1.d0)/2.d-6
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k) + pgam(k)
         !    qsout2(k) = qsout2(k) + lammin
         !    nrout2(k) = nrout2(k) + lammax
         !    !nsout2(k) = nsout2(k) + qcic(k)
         !    !drout2(k) = drout2(k) + tmp1
         !    !dsout2(k) = dsout2(k) + tmp2
         !endif
        
                    if (lamc(k).lt.lammin) then
                       lamc(k) = lammin
                       ncic(k) = 6.d0*lamc(k)*lamc(k)*lamc(k)*qcic(k)* &
                            shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/(pi*rhow*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama))
                    else if (lamc(k).gt.lammax) then
                       lamc(k) = lammax
                       ncic(k) = 6.d0*lamc(k)*lamc(k)*lamc(k)*qcic(k)* &
                            shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/(pi*rhow*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama))
                    end if
        
                    ! parameter to calculate droplet freezing
        
                    cdist1(k) = ncic(k)/shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama) 
        
                 else
                    lamc(k) = 0.d0
                    cdist1(k) = 0.d0
                 end if

         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+ncic(k)
         !    qsout2(k) = qsout2(k)+pgam(k)
         !    nrout2(k) = nrout2(k)+lamc(k)
         !    nsout2(k) = nsout2(k)+pgam(k)
         !    drout2(k) = drout2(k)+lammin
         !    dsout2(k) = dsout2(k)+lammax
         !endif
        
                 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 ! begin micropysical process calculations 
                 !.................................................................
                 ! autoconversion of cloud liquid water to rain
                 ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
                 ! minimum qc of 1 x 10^-8 prevents floating point error
        
                 if (qcic(k).ge.1.d-8) then
        
                    ! nprc is increase in rain number conc due to autoconversion
                    ! nprc1 is decrease in cloud droplet conc due to autoconversion
        
                    ! assume exponential sub-grid distribution of qc, resulting in additional
                    ! factor related to qcvar below
        
                    ! hm switch for sub-columns, don't include sub-grid qc
                    if (microp_uniform) then
        
                       call math_agent_2i1o(power, qcic(k), dp2_47, tmp1)
                       tmp2 = (ncic(k)/1.d6*rho(k))
                       call math_agent_2i1o(power, tmp2, dp_1_79, tmp3)
                       prc(k) = dp1350_0*tmp1*tmp3
                       !prc(k) = dp1350_0*qcic(k)**dp2_47* &
                       !     (ncic(k)/1.d6*rho(k))**(dp_1_79)
                       !nprc(k) = prc(k)/(4.d0/3.d0*pi*rhow*(dp25d_6)**3)
                       nprc(k) = prc(k)/(4.d0/3.d0*pi*rhow*(dp25d_6)*(dp25d_6)*(dp25d_6))
                       nprc1(k) = prc(k)/(qcic(k)/ncic(k))
        
                    else
        
                       call math_agent_2i1o(power, qcic(k), dp2_47, tmp1)
                       tmp2 = (ncic(k)/1.d6*rho(k))
                       call math_agent_2i1o(power, tmp2, dp_1_79, tmp3)
                       prc(k) = cons2/(cons3*cons18)*dp1350_0*tmp1* &
                            tmp3
                       !prc(k) = cons2/(cons3*cons18)*dp1350_0*qcic(k)**dp2_47* &
                       !     (ncic(k)/1.d6*rho(k))**(dp_1_79)
                       nprc(k) = prc(k)/cons22
                       nprc1(k) = prc(k)/(qcic(k)/ncic(k))
        
                    end if               ! sub-column switch
        
                 else
                    prc(k)=0.d0
                    nprc(k)=0.d0
                    nprc1(k)=0.d0
                 end if
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k) + qcic(k) 
         !    qsout2(k) = qsout2(k) + ncic(k)
         !    nrout2(k) = nrout2(k) + rho(k)
         !    nsout2(k) = nsout2(k) + prc(k)
         !    drout2(k) = drout2(k) + qric(k)
         !    dsout2(k) = dsout2(k) + nric(k)
         !endif
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k) + qcic(k) 
         !    qsout2(k) = qsout2(k) + nprc1(k)
         !    nrout2(k) = nrout2(k) + prc(k)
         !    nsout2(k) = nsout2(k) + rho(k)
         !    drout2(k) = drout2(k) + qric(k)
         !    dsout2(k) = dsout2(k) + nric(k)
         !endif
                 ! add autoconversion to precip from above to get provisional rain mixing ratio
                 ! and number concentration (qric and nric)
        
                 ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)
        
                 dum=0.45d0
                 dum1=0.45d0
        
                 !if ( myid <= 30 ) then
                 !    tmpn(1 ) = umr   (pver)
                 !    tmpn(2 ) = unr   (pver)
                 !    tmpn(3 ) = qniic (pver)
                 !    tmpn(4 ) = rho   (pver)
                 !    tmpn(5 ) = cldmax(pver)
                 !    tmpn(6 ) = pra   (pver)
                 !    tmpn(7 ) = pre   (pver)
                 !    tmpn(8 ) = pracs (pver)
                 !    tmpn(9 ) = mnuccr(pver)
                 !    tmpn(10) = nric  (pver)
                 !    tmpn(11) = nsubr (pver)
                 !    tmpn(12) = npracs(pver)
                 !    tmpn(13) = nnuccr(pver)
                 !    tmpn(14) = nragg (pver)
                 !    tmpn(15) = qric  (pver)
                 !    tmpn(16) = pracs (pver)
                 !end if
                 !if ( it .eq. 1 .and. (k.eq.pver )  ) then
                 !    qrout2(1) = qrout2(1) + tmpn(1) 
                 !    qsout2(1) = qsout2(1) + tmpn(2) 
                 !    nrout2(1) = nrout2(1) + tmpn(3)  
                 !    nsout2(1) = nsout2(1) + tmpn(4) 
                 !    drout2(1) = drout2(1) + tmpn(10)
                 !    dsout2(1) = dsout2(1) + tmpn(15)
                 !endif

         !if ( it .eq. 1 .and. k <=8 .and. k > 1 ) then
         !    qrout2(k) = qrout2(k) + unr(k-1) 
         !    qsout2(k) = qsout2(k) + nric(k-1)
         !    nrout2(k) = nrout2(k) + nsubr(k-1)
         !    nsout2(k) = nsout2(k) + npracs(k-1)
         !    drout2(k) = drout2(k) + nnuccr(k-1)
         !    dsout2(k) = dsout2(k) + nragg(k-1)
         !endif
         !if ( it .eq. 1 .and. k > 1) then
         !    qrout2(k) = qrout2(k) + umr(k-1) 
         !    qsout2(k) = qsout2(k) + qric(k-1)
         !    nrout2(k) = nrout2(k) + pra(k-1)
         !    nsout2(k) = nsout2(k) + pre(k-1)
         !    drout2(k) = drout2(k) + pracs(k-1)
         !    dsout2(k) = dsout2(k) + mnuccr(k-1)
         !endif

                 if ( myid < 32 ) then
                     if (k.eq.top_lev) then
                         qric(k)=prc(k)*lcldm(k)*dz(k)/cldmax(k)/dum
                         nric(k)=nprc(k)*lcldm(k)*dz(k)/cldmax(k)/dum
                     else
                         if (qric(k-1).ge.qsmall) then
                             dum=umr(k-1)
                             dum1=unr(k-1)
                         end if
                         if (qric(k-1).ge.1.d-9.or.qniic(k-1).ge.1.d-9) then
                             nprc(k)=0.d0
                         end if
                         qric(k) = (rho(k-1)*umr(k-1)*qric(k-1)*cldmax(k-1)+ &
                         (rho(k)*dz(k)*((pra(k-1)+prc(k))*lcldm(k)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax(k))))&
                         /(dum*rho(k)*cldmax(k))
                         nric(k) = (rho(k-1)*unr(k-1)*nric(k-1)*cldmax(k-1)+ &
                         (rho(k)*dz(k)*(nprc(k)*lcldm(k)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax(k))))&
                         /(dum1*rho(k)*cldmax(k))
                     endif
                 else
                     if (k.eq.top_lev) then
                         if (tmpn(15).ge.qsmall) then
                             dum=tmpn(1)
                             dum1=tmpn(2)
                         end if
                         if (tmpn(15).ge.1.d-9.or.tmpn(3).ge.1.d-9) then
                             nprc(k)=0.d0
                         end if
                         qric(k) = (tmpn(4)*tmpn(1)*tmpn(15)*tmpn(5)+ &
                         (rho(k)*dz(k)*((tmpn(6)+prc(k))*lcldm(k)+(tmpn(7)-tmpn(16)-tmpn(9))*cldmax(k))))&
                         /(dum*rho(k)*cldmax(k))
                         nric(k) = (tmpn(4)*tmpn(2)*tmpn(10)*tmpn(5)+ &
                         (rho(k)*dz(k)*(nprc(k)*lcldm(k)+(tmpn(11)-tmpn(12)-tmpn(13)+tmpn(14))*cldmax(k))))&
                         /(dum1*rho(k)*cldmax(k))
                     else 
                         if (qric(k-1).ge.qsmall) then
                             dum=umr(k-1)
                             dum1=unr(k-1)
                         end if
                         if (qric(k-1).ge.1.d-9.or.qniic(k-1).ge.1.d-9) then
                             nprc(k)=0.d0
                         end if
                         qric(k) = (rho(k-1)*umr(k-1)*qric(k-1)*cldmax(k-1)+ &
                         (rho(k)*dz(k)*((pra(k-1)+prc(k))*lcldm(k)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax(k))))&
                         /(dum*rho(k)*cldmax(k))
                         nric(k) = (rho(k-1)*unr(k-1)*nric(k-1)*cldmax(k-1)+ &
                         (rho(k)*dz(k)*(nprc(k)*lcldm(k)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax(k))))&
                         /(dum1*rho(k)*cldmax(k))
                     end if
                 endif
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k) + qric(k)
         !endif
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k) + niic(k)
         !    qsout2(k) = qsout2(k) + berg(k)
         !    nrout2(k) = nrout2(k) + qcic(k)
         !    nsout2(k) = nsout2(k) + ncic(k)
         !    drout2(k) = drout2(k) + qric(k)
         !    dsout2(k) = dsout2(k) + nric(k)
         !endif
                 !.......................................................................
                 ! Autoconversion of cloud ice to snow
                 ! similar to Ferrier (1994)
        
                 if (do_cldice) then
                    if (t(k).le.dp273_15.and.qiic(k).ge.qsmall) then
        
                       ! note: assumes autoconversion timescale of 180 sec
                       
                       tmp1 = -lami(k)*dcs
                       call math_agent_1i1o(exp_agent, tmp1, tmp2)
                       nprci(k) = n0i(k)/(lami(k)*dp180_0)*tmp2
                       !nprci(k) = n0i(k)/(lami(k)*dp180_0)*exp(-lami(k)*dcs)
                     
                       prci(k) = pi*rhoi*n0i(k)/(6.d0*dp180_0)* &
                            (cons23/lami(k)+3.d0*cons24/(lami(k)*lami(k))+ &
                            6.d0*dcs/(lami(k)*lami(k)*lami(k))+6.d0/(lami(k)*lami(k)*lami(k)*lami(k)))*tmp2
                       !prci(k) = pi*rhoi*n0i(k)/(6.d0*dp180_0)* &
                       !     (cons23/lami(k)+3.d0*cons24/lami(k)**2+ &
                       !     6.d0*dcs/lami(k)**3+6.d0/lami(k)**4)*exp(-lami(k)*dcs)
                    else
                       prci(k)=0.d0
                       nprci(k)=0.d0
                    end if
                 else
                    ! Add in the particles that we have already converted to snow, and
                    ! don't do any further autoconversion of ice.
                    prci(k)  = tnd_qsnow( k) / cldm(k)
                    nprci(k) = tnd_nsnow( k) / cldm(k)
                 end if
        
                 ! add autoconversion to flux from level above to get provisional snow mixing ratio
                 ! and number concentration (qniic and nsic)
        
                 dum=(asn(k)*cons25)
                 dum1=(asn(k)*cons25)
        
         !if ( it .eq. 1 .and. k > 1 ) then
         !    qrout2(k) = qrout2(k) + uns(k-1)
         !    qsout2(k) = qsout2(k) + nprci(k)
         !    nrout2(k) = nrout2(k) + nsubs(k-1)
         !    nsout2(k) = nsout2(k) + nsagg(k-1)
         !    drout2(k) = drout2(k) + nnuccr(k-1)
         !    dsout2(k) = dsout2(k) + rho(k-1)
         !endif

             if ( myid < 32 ) then
                 if (k.eq.top_lev) then
                    qniic(k)=prci(k)*icldm(k)*dz(k)/cldmax(k)/dum
                    nsic(k)=nprci(k)*icldm(k)*dz(k)/cldmax(k)/dum
                 else
                    if (qniic(k-1).ge.qsmall) then
                       dum=ums(k-1)
                       dum1=uns(k-1)
                    end if
        
                    qniic(k) = (rho(k-1)*ums(k-1)*qniic(k-1)*cldmax(k-1)+ &
                         (rho(k)*dz(k)*((prci(k)+prai(k-1)+psacws(k-1)+bergs(k-1))*icldm(k)+(prds(k-1)+ &
                         pracs(k-1)+mnuccr(k-1))*cldmax(k))))&
                         /(dum*rho(k)*cldmax(k))
        
                    nsic(k) = (rho(k-1)*uns(k-1)*nsic(k-1)*cldmax(k-1)+ &
                         (rho(k)*dz(k)*(nprci(k)*icldm(k)+(nsubs(k-1)+nsagg(k-1)+nnuccr(k-1))*cldmax(k))))&
                         /(dum1*rho(k)*cldmax(k))
        
                 end if
            else 
                 if (k.eq.top_lev) then
                    if (tmpn(3).ge.qsmall) then
                       dum=tmpn(8)
                       dum1=tmpn(17)
                    end if
        
                    qniic(k) = (tmpn(4)*tmpn(8)*tmpn(3)*tmpn(5)+ &
                         (rho(k)*dz(k)*((prci(k)+tmpn(18)+tmpn(19)+tmpn(20))*icldm(k)+(tmpn(21)+ &
                         tmpn(16)+tmpn(9))*cldmax(k))))&
                         /(dum*rho(k)*cldmax(k))
        
                    nsic(k) = (tmpn(4)*tmpn(17)*tmpn(22)*tmpn(5)+ &
                         (rho(k)*dz(k)*(nprci(k)*icldm(k)+(tmpn(23)+tmpn(24)+tmpn(13))*cldmax(k))))&
                         /(dum1*rho(k)*cldmax(k))
                 else
                    if (qniic(k-1).ge.qsmall) then
                       dum=ums(k-1)
                       dum1=uns(k-1)
                    end if
        
                    qniic(k) = (rho(k-1)*ums(k-1)*qniic(k-1)*cldmax(k-1)+ &
                         (rho(k)*dz(k)*((prci(k)+prai(k-1)+psacws(k-1)+bergs(k-1))*icldm(k)+(prds(k-1)+ &
                         pracs(k-1)+mnuccr(k-1))*cldmax(k))))&
                         /(dum*rho(k)*cldmax(k))
        
                    nsic(k) = (rho(k-1)*uns(k-1)*nsic(k-1)*cldmax(k-1)+ &
                         (rho(k)*dz(k)*(nprci(k)*icldm(k)+(nsubs(k-1)+nsagg(k-1)+nnuccr(k-1))*cldmax(k))))&
                         /(dum1*rho(k)*cldmax(k))
        
                 end if
             endif

        
         !!if ( it .eq. 1 .and. k .eq. 1 ) then
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k) + dum
         !    qsout2(k) = qsout2(k) + dum1
         !    nrout2(k) = nrout2(k) + prci(k)
         !    nsout2(k) = nsout2(k) + nprci(k)
         !    drout2(k) = drout2(k) + qniic(k)
         !    dsout2(k) = dsout2(k) + nsic(k)
         !endif
                 ! if precip mix ratio is zero so should number concentration
        
                 if (qniic(k).lt.qsmall) then
                    qniic(k)=0.d0
                    nsic(k)=0.d0
                 end if
        
                 if (qric(k).lt.qsmall) then
                    qric(k)=0.d0
                    nric(k)=0.d0
                 end if
         !if ( it .eq. 1 ) then
         !    qsout2(k) = qsout2(k) + qric(k)
         !endif
        
                 ! make sure number concentration is a positive number to avoid 
                 ! taking root of negative later
        
                 nric(k)=max(nric(k),0.d0)
                 nsic(k)=max(nsic(k),0.d0)
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+prci(k)
         !    qsout2(k) = qsout2(k)+nprci(k)
         !    nrout2(k) = nrout2(k)+qniic(k)
         !    nsout2(k) = nsout2(k)+qric(k)
         !    drout2(k) = drout2(k)+nric(k)
         !    dsout2(k) = dsout2(k)+nsic(k)
         !endif
                 !.......................................................................
                 ! get size distribution parameters for precip
                 !......................................................................
                 ! rain
        
                 if (qric(k).ge.qsmall) then
                    tmp1 = (pi*rhow*nric(k)/qric(k))
                    tmp2 = 1.d0/3.d0
                    call math_agent_2i1o(power, tmp1, tmp2, lamr(k))
                    !lamr(k) = (pi*rhow*nric(k)/qric(k))**(1.d0/3.d0)
                    n0r(k) = nric(k)*lamr(k)
        
                    ! check for slope
                    ! adjust vars
        
                    if (lamr(k).lt.lamminr) then
        
                       lamr(k) = lamminr
        
                       !n0r(k) = lamr(k)**4*qric(k)/(pi*rhow)
                       n0r(k) = lamr(k)*lamr(k)*lamr(k)*lamr(k)*qric(k)/(pi*rhow)
                       nric(k) = n0r(k)/lamr(k)
                    else if (lamr(k).gt.lammaxr) then
                       lamr(k) = lammaxr
                       !n0r(k) = lamr(k)**4*qric(k)/(pi*rhow)
                       n0r(k) = lamr(k)*lamr(k)*lamr(k)*lamr(k)*qric(k)/(pi*rhow)
                       nric(k) = n0r(k)/lamr(k)
                    end if
        
                    ! provisional rain number and mass weighted mean fallspeed (m/s)
        
                    call math_agent_2i1o(power, lamr(k), br, tmp1)
                    unr(k) = min(arn(k)*cons4/tmp1,dp9_1*rhof(k))
                    umr(k) = min(arn(k)*cons5/(6.d0*tmp1),dp9_1*rhof(k))
                    !unr(k) = min(arn(k)*cons4/lamr(k)**br,dp9_1*rhof(k))
                    !umr(k) = min(arn(k)*cons5/(6.d0*lamr(k)**br),dp9_1*rhof(k))
        
                 else
                    lamr(k) = 0.d0
                    n0r(k) = 0.d0
                    umr(k) = 0.d0
                    unr(k) = 0.d0
                 end if
        
                 !......................................................................
                 ! snow
        
                 if (qniic(k).ge.qsmall) then
                    tmp1 = (cons6*cs*nsic(k)/qniic(k))
                    tmp2 = 1.d0/ds
                    call math_agent_2i1o(power, tmp1, tmp2, lams(k))
                    !lams(k) = (cons6*cs*nsic(k)/qniic(k))**(1.d0/ds)
                    n0s(k) = nsic(k)*lams(k)
        
                    ! check for slope
                    ! adjust vars
        
                    if (lams(k).lt.lammins) then
                       lams(k) = lammins
                       tmp1 = ds+1.d0
                       call math_agent_2i1o(power, lams(k), tmp1, tmp2)
                       n0s(k) = tmp2*qniic(k)/(cs*cons6)
                       !n0s(k) = lams(k)**(ds+1.d0)*qniic(k)/(cs*cons6)
                       nsic(k) = n0s(k)/lams(k)
        
                    else if (lams(k).gt.lammaxs) then
                       lams(k) = lammaxs
                       tmp1 = ds+1.d0
                       call math_agent_2i1o(power, lams(k), tmp1, tmp2)
                       n0s(k) = tmp2*qniic(k)/(cs*cons6)
                       !n0s(k) = lams(k)**(ds+1.d0)*qniic(k)/(cs*cons6)
                       nsic(k) = n0s(k)/lams(k)
                    end if
        
                    ! provisional snow number and mass weighted mean fallspeed (m/s)
        
                    call math_agent_2i1o(power, lams(k), bs, tmp1)
                    ums(k) = min(asn(k)*cons8/(6.d0*tmp1),1.2d0*rhof(k))
                    uns(k) = min(asn(k)*cons7/tmp1,1.2d0*rhof(k))
                    !ums(k) = min(asn(k)*cons8/(6.d0*lams(k)**bs),1.2d0*rhof(k))
                    !uns(k) = min(asn(k)*cons7/lams(k)**bs,1.2d0*rhof(k))
        
                 else
                    lams(k) = 0.d0
                    n0s(k) = 0.d0
                    ums(k) = 0.d0
                    uns(k) = 0.d0
                 end if
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+umr(k)
         !    qsout2(k) = qsout2(k)+unr(k)
         !    nrout2(k) = nrout2(k)+lams(k)
         !    nsout2(k) = nsout2(k)+n0s(k)
         !    drout2(k) = drout2(k)+ums(k)
         !    dsout2(k) = dsout2(k)+uns(k)
         !endif
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+cdist1(k)
         !    qsout2(k) = qsout2(k)+pgam(k)
         !    nrout2(k) = nrout2(k)+t(k)
         !    nsout2(k) = nsout2(k)+lamc(k)
         !    !drout2(k) = drout2(k)+use_hetfrz_classnuc
         !    !dsout2(k) = dsout2(k)+do_cldice
         !endif
                 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
                 ! heterogeneous freezing of cloud water
         !nnuccc(k) = 0.d0
         !mnuccc(k) = 0.d0

         !nnucct(k) = 0.d0
         !mnucct(k) = 0.d0

         !nnudep(k) = 0.d0
         !mnudep(k) = 0.d0
        
                 if (.not. use_hetfrz_classnuc) then
        
                    if (do_cldice .and. qcic(k).ge.qsmall .and. t(k).lt.dp269_15) then
        
                       ! immersion freezing (Bigg, 1953)
        
        
                       ! subcolumns
        
                       if (microp_uniform) then
        
                          tmp1 = (aimm*(dp273_15-t(k)))
                          call math_agent_1i1o(exp_agent, tmp1, tmp2)
                          mnuccc(k) = &
                             pi*pi/36.d0*rhow* &
                             cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(7.d0+pgam(k),math_agent,Pgama,Qgama,Cgama,Sgama)* &
                             bimm*(tmp2-1.d0)/ &
                             (lamc(k)*lamc(k)*lamc(k))/(lamc(k)*lamc(k)*lamc(k))
                          !mnuccc(k) = &
                          !   pi*pi/36.d0*rhow* &
                          !   cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(7.d0+pgam(k))* &
                          !   bimm*(exp(aimm*(dp273_15-t(k)))-1.d0)/ &
                          !   lamc(k)**3/lamc(k)**3
        
                          nnuccc(k) = &
                             pi/6.d0*cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama) &
                             *bimm* &
                             (tmp2-1.d0)/(lamc(k)*lamc(k)*lamc(k))
                          !nnuccc(k) = &
                          !   pi/6.d0*cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0) &
                          !   *bimm* &
                          !   (exp(aimm*(dp273_15-t(k)))-1.d0)/lamc(k)**3
        
                       else
        
                          tmp1 = (aimm*(dp273_15-t(k)))
                          call math_agent_1i1o(exp_agent, tmp1, tmp2)
                          mnuccc(k) = cons9/(cons3*cons19)* &
                             pi*pi/36.d0*rhow* &
                             cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(7.d0+pgam(k),math_agent,Pgama,Qgama,Cgama,Sgama)* &
                             bimm*(tmp2-1.d0)/ &
                             (lamc(k)*lamc(k)*lamc(k))/(lamc(k)*lamc(k)*lamc(k))
        
                          nnuccc(k) = cons10/(cons3*qcvar)* &
                             pi/6.d0*cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama) &
                             *bimm* &
                             (tmp2-1.d0)/(lamc(k)*lamc(k)*lamc(k))

                          !mnuccc(k) = 1.d0
                          !nnuccc(k) = 1.d0
                          !mnuccc(k) = cons9/(cons3*cons19)* &
                          !   pi*pi/36.d0*rhow* &
                          !   cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(7.d0+pgam(k))* &
                          !   bimm*(exp(aimm*(dp273_15-t(k)))-1.d0)/ &
                          !   lamc(k)**3/lamc(k)**3
        
                          !nnuccc(k) = cons10/(cons3*qcvar)* &
                          !   pi/6.d0*cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0) &
                          !   *bimm* &
                          !   (exp(aimm*(dp273_15-t(k)))-1.d0)/lamc(k)**3
                       end if           ! sub-columns
        
         !if ( it .eq. 1 ) then
         !    !qrout2(k) = qrout2(k)+nnuccc(k)
         !    !qsout2(k) = qsout2(k)+mnuccc(k)
         !    !nrout2(k) = nrout2(k)+nnucct(k)
         !    !nsout2(k) = nsout2(k)+mnucct(k)
         !    drout2(k) = drout2(k)+nnuccc(k)
         !    dsout2(k) = dsout2(k)+mnuccc(k)
         !endif
        
                       ! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
                       ! dust size and number in 4 bins are read in from companion routine
        
                       tmp1 = (dp270_16-t(k))
                       call math_agent_2i1o(power, tmp1, 1.3d0, tcnt)
                       tmp1 = (t(k)/298.0d0)
                       call math_agent_2i1o(power, tmp1, dp0_85, tmp2)
                       viscosity=dp1_8d_5*tmp2    ! Viscosity (kg/m/s)
                       tmp1 = (8.0d0*dp28_96d_3/(pi*dp8_314409*t(k)))
                       call math_agent_1i1o(sqrt_agent, tmp1, tmp2)
                       mfp=2.0d0*viscosity/(p(k)  &                   ! Mean free path (m)
                          *tmp2)           

                       !tcnt=(dp270_16-t(k))**1.3d0
                       !viscosity=dp1_8d_5*(t(k)/298.0d0)**dp0_85    ! Viscosity (kg/m/s)
                       !mfp=2.0d0*viscosity/(p(k)  &                   ! Mean free path (m)
                       !   *sqrt(8.0d0*dp28_96d_3/(pi*dp8_314409*t(k))))           
        
                       tmp1 = (-(dp1_1*rndst(k,1)/mfp))
                       call math_agent_1i1o(exp_agent, tmp1, tmp2)
                       nslip1=1.0d0+(mfp/rndst(k,1))*(dp1_257+(0.4d0*tmp2))! Slip correction factor
                       tmp1 = (-(dp1_1*rndst(k,2)/mfp))
                       call math_agent_1i1o(exp_agent, tmp1, tmp2)
                       nslip2=1.0d0+(mfp/rndst(k,2))*(dp1_257+(0.4d0*tmp2))
                       tmp1 = (-(dp1_1*rndst(k,3)/mfp))
                       call math_agent_1i1o(exp_agent, tmp1, tmp2)
                       nslip3=1.0d0+(mfp/rndst(k,3))*(dp1_257+(0.4d0*tmp2))
                       tmp1 = (-(dp1_1*rndst(k,4)/mfp))
                       call math_agent_1i1o(exp_agent, tmp1, tmp2)
                       nslip4=1.0d0+(mfp/rndst(k,4))*(dp1_257+(0.4d0*tmp2))
                       !nslip1=1.0d0+(mfp/rndst(k,1))*(dp1_257+(0.4d0*Exp(-(dp1_1*rndst(k,1)/mfp))))! Slip correction factor
                       !nslip2=1.0d0+(mfp/rndst(k,2))*(dp1_257+(0.4d0*Exp(-(dp1_1*rndst(k,2)/mfp))))
                       !nslip3=1.0d0+(mfp/rndst(k,3))*(dp1_257+(0.4d0*Exp(-(dp1_1*rndst(k,3)/mfp))))
                       !nslip4=1.0d0+(mfp/rndst(k,4))*(dp1_257+(0.4d0*Exp(-(dp1_1*rndst(k,4)/mfp))))
        
                       ndfaer1=dp1_381d_23*t(k)*nslip1/(6.d0*pi*viscosity*rndst(k,1))  ! aerosol diffusivity (m2/s)
                       ndfaer2=dp1_381d_23*t(k)*nslip2/(6.d0*pi*viscosity*rndst(k,2))
                       ndfaer3=dp1_381d_23*t(k)*nslip3/(6.d0*pi*viscosity*rndst(k,3))
                       ndfaer4=dp1_381d_23*t(k)*nslip4/(6.d0*pi*viscosity*rndst(k,4))
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+nslip1
         !    qsout2(k) = qsout2(k)+nslip2 
         !    nrout2(k) = nrout2(k)+nslip3 
         !    nsout2(k) = nsout2(k)+ndfaer1 
         !    drout2(k) = drout2(k)+ndfaer2 
         !    dsout2(k) = dsout2(k)+ndfaer3 
         !endif
        
                       if (microp_uniform) then
        
                          mnucct(k) = &
                             (ndfaer1*(nacon(k,1)*tcnt)+ndfaer2*(nacon(k,2)*tcnt)+ &
                             ndfaer3*(nacon(k,3)*tcnt)+ndfaer4*(nacon(k,4)*tcnt))*pi*pi/3.d0*rhow* &
                             cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+5.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/(lamc(k)*lamc(k)*lamc(k)*lamc(k))
                          !mnucct(k) = &
                          !   (ndfaer1*(nacon(k,1)*tcnt)+ndfaer2*(nacon(k,2)*tcnt)+ &
                          !   ndfaer3*(nacon(k,3)*tcnt)+ndfaer4*(nacon(k,4)*tcnt))*pi*pi/3.d0*rhow* &
                          !   cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+5.d0)/lamc(k)**4
        
                          nnucct(k) = (ndfaer1*(nacon(k,1)*tcnt)+ndfaer2*(nacon(k,2)*tcnt)+ &
                             ndfaer3*(nacon(k,3)*tcnt)+ndfaer4*(nacon(k,4)*tcnt))*2.d0*pi*  &
                             cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+2.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/lamc(k)
        
                       else
        
                          tmp1 = 4.d0/3.d0
                          call math_agent_2i1o(power, qcvar, tmp1, tmp2)
                          mnucct(k) = shr_spfn_gamma_nonintrinsic_r8(qcvar+4.d0/3.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/(cons3*tmp2)*  &
                             (ndfaer1*(nacon(k,1)*tcnt)+ndfaer2*(nacon(k,2)*tcnt)+ &
                             ndfaer3*(nacon(k,3)*tcnt)+ndfaer4*(nacon(k,4)*tcnt))*pi*pi/3.d0*rhow* &
                             cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+5.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/(lamc(k)*lamc(k)*lamc(k)*lamc(k))
                          !mnucct(k) = shr_spfn_gamma_nonintrinsic_r8(qcvar+4.d0/3.d0)/(cons3*qcvar**(4.d0/3.d0))*  &
                          !   (ndfaer1*(nacon(k,1)*tcnt)+ndfaer2*(nacon(k,2)*tcnt)+ &
                          !   ndfaer3*(nacon(k,3)*tcnt)+ndfaer4*(nacon(k,4)*tcnt))*pi*pi/3.d0*rhow* &
                          !   cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+5.d0)/lamc(k)**4
                          tmp1 = 1.d0/3.d0
                          call math_agent_2i1o(power, qcvar, tmp1, tmp2)
                          nnucct(k) =  shr_spfn_gamma_nonintrinsic_r8(qcvar+1.d0/3.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/(cons3*tmp2)*  &
                             (ndfaer1*(nacon(k,1)*tcnt)+ndfaer2*(nacon(k,2)*tcnt)+ &
                             ndfaer3*(nacon(k,3)*tcnt)+ndfaer4*(nacon(k,4)*tcnt))*2.d0*pi*  &
                             cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+2.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/lamc(k)
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+nacon(k,1)
         !    qsout2(k) = qsout2(k)+nacon(k,2)
         !    nrout2(k) = nrout2(k)+mnucct(k) 
         !    nsout2(k) = nsout2(k)+nnucct(k) 
         !    drout2(k) = drout2(k)+nacon(k,3)
         !    dsout2(k) = dsout2(k)+nacon(k,4)
         !endif
                          !nnucct(k) =  shr_spfn_gamma_nonintrinsic_r8(qcvar+1.d0/3.d0)/(cons3*qcvar**(1.d0/3.d0))*  &
                          !   (ndfaer1*(nacon(k,1)*tcnt)+ndfaer2*(nacon(k,2)*tcnt)+ &
                          !   ndfaer3*(nacon(k,3)*tcnt)+ndfaer4*(nacon(k,4)*tcnt))*2.d0*pi*  &
                          !   cdist1(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+2.d0)/lamc(k)
        
                       end if      ! sub-column switch
        
                       ! make sure number of droplets frozen does not exceed available ice nuclei concentration
                       ! this prevents 'runaway' droplet freezing
        
                       if (nnuccc(k)*lcldm(k).gt.nnuccd(k)) then
                          dum=(nnuccd(k)/(nnuccc(k)*lcldm(k)))
                          ! scale mixing ratio of droplet freezing with limit
                          mnuccc(k)=mnuccc(k)*dum
                          nnuccc(k)=nnuccd(k)/lcldm(k)
                       end if
        
                    else
                       mnuccc(k)=0.d0
                       nnuccc(k)=0.d0
                       mnucct(k)=0.d0
                       nnucct(k)=0.d0
                    end if
                 else
                    if (do_cldice .and. qcic(k) >= qsmall) then

                       tmp1 = (1.333d0*pi)
                       call math_agent_2i1o(power, tmp1, 0.333d0, tmp2)
                       con1 = 1.d0/tmp2

                       tmp1 = (rho(k)*qcic(k)/(rhow*max(ncic(k)*rho(k), 1.0d6)))
                       call math_agent_2i1o(power, tmp1, 0.333d0, tmp2)
                       r3lx = con1*tmp2

                       !con1 = 1.d0/(1.333d0*pi)**0.333d0
                       !r3lx = con1*(rho(k)*qcic(k)/(rhow*max(ncic(k)*rho(k), 1.0d6)))**0.333d0 ! in m
                       r3lx = max(4.d-6, r3lx)
                       mi0l = 4.d0/3.d0*pi*rhow*r3lx*r3lx*r3lx
                        
                       nnuccc(k) = frzimm(k)*1.0d6/rho(k)
                       mnuccc(k) = nnuccc(k)*mi0l 
        
                       nnucct(k) = frzcnt(k)*1.0d6/rho(k)
                       mnucct(k) = nnucct(k)*mi0l 
        
                       nnudep(k) = frzdep(k)*1.0d6/rho(k)
                       mnudep(k) = nnudep(k)*mi0

         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+mnuccc(k)
         !    qsout2(k) = qsout2(k)+nnuccc(k) 
         !    nrout2(k) = nrout2(k)+mnucct(k) 
         !    nsout2(k) = nsout2(k)+nnucct(k) 
         !    drout2(k) = drout2(k)+mnudep(k) 
         !    dsout2(k) = dsout2(k)+nnudep(k)
         !endif
        
                    else
                       nnuccc(k) = 0.d0
                       mnuccc(k) = 0.d0
        
                       nnucct(k) = 0.d0
                       mnucct(k) = 0.d0
        
                       nnudep(k) = 0.d0
                       mnudep(k) = 0.d0
                    end if
                 endif
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+nnuccc(k)
         !    qsout2(k) = qsout2(k)+mnuccc(k)
         !    nrout2(k) = nrout2(k)+nnucct(k)
         !    nsout2(k) = nsout2(k)+mnucct(k)
         !    drout2(k) = drout2(k)+nnudep(k)
         !    dsout2(k) = dsout2(k)+mnudep(k)
         !endif
                 !.......................................................................
                 ! snow self-aggregation from passarelli, 1978, used by reisner, 1998
                 ! this is hard-wired for bs = 0.4 for now
                 ! ignore self-collection of cloud ice
        
                 if (qniic(k).ge.qsmall .and. t(k).le.dp273_15) then
                    tmp1 = ((1.d0-bs)/3.d0)
                    call math_agent_2i1o(power, pi, tmp1, tmp2)
                    tmp1 = ((-2.d0-bs)/3.d0)
                    call math_agent_2i1o(power, rhosn, tmp1, tmp3)
                    tmp1 = ((2.d0+bs)/3.d0)
                    call math_agent_2i1o(power, rho(k), tmp1, tmp4)
                    call math_agent_2i1o(power, qniic(k), tmp1, tmp5)
                    tmp1 = (nsic(k)*rho(k))
                    tmp6 = ((4.d0-bs)/3.d0)
                    call math_agent_2i1o(power, tmp1, tmp6, tmp7)
                    nsagg(k) = -1108.d0*asn(k)*Eii* &
                         tmp2*tmp3*tmp4*tmp5* &
                         tmp7/ &
                         (4.d0*720.d0*rho(k))
                    !nsagg(k) = -1108.d0*asn(k)*Eii* &
                    !     pi**((1.d0-bs)/3.d0)*rhosn**((-2.d0-bs)/3.d0)*rho(k)** &
                    !     ((2.d0+bs)/3.d0)*qniic(k)**((2.d0+bs)/3.d0)* &
                    !     (nsic(k)*rho(k))**((4.d0-bs)/3.d0)/ &
                    !     (4.d0*720.d0*rho(k))
                 else
                    nsagg(k)=0.d0
                 end if
        
                 !.......................................................................
                 ! accretion of cloud droplets onto snow/graupel
                 ! here use continuous collection equation with
                 ! simple gravitational collection kernel
                 ! ignore collisions between droplets/cloud ice
                 ! since minimum size ice particle for accretion is 50 - 150 micron
        
                 ! ignore collision of snow with droplets above freezing
        
                 if (qniic(k).ge.qsmall .and. t(k).le.tmelt .and. &
                      qcic(k).ge.qsmall) then
        
                    ! put in size dependent collection efficiency
                    ! mean diameter of snow is area-weighted, since
                    ! accretion is function of crystal geometric area
                    ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)
        
                    dc0 = (pgam(k)+1.d0)/lamc(k)
                    ds0 = 1.d0/lams(k)
                    dum = dc0*dc0*uns(k)*rhow/(9.d0*mu(k)*ds0)
                    eci = dum*dum/((dum+0.4d0)*(dum+0.4d0))
        
                    eci = max(eci,0.d0)
                    eci = min(eci,1.d0)
        
        
                    ! no impact of sub-grid distribution of qc since psacws
                    ! is linear in qc
        
                    tmp1 = bs+3.d0
                    call math_agent_2i1o(power, lams(k), tmp1, tmp2)
                    psacws(k) = pi/4.d0*asn(k)*qcic(k)*rho(k)* &
                         n0s(k)*Eci*cons11/ &
                         tmp2
                    npsacws(k) = pi/4.d0*asn(k)*ncic(k)*rho(k)* &
                         n0s(k)*Eci*cons11/ &
                         tmp2
                    !psacws(k) = pi/4.d0*asn(k)*qcic(k)*rho(k)* &
                    !     n0s(k)*Eci*cons11/ &
                    !     lams(k)**(bs+3.d0)
                    !npsacws(k) = pi/4.d0*asn(k)*ncic(k)*rho(k)* &
                    !     n0s(k)*Eci*cons11/ &
                    !     lams(k)**(bs+3.d0)
                 else
                    psacws(k)=0.d0
                    npsacws(k)=0.d0
                 end if
        
                 ! add secondary ice production due to accretion of droplets by snow 
                 ! (Hallet-Mossop process) (from Cotton et al., 1986)
        
                 if (.not. do_cldice) then
                    ni_secp   = 0.0d0
                    nsacwi(k) = 0.0d0
                    msacwi(k) = 0.0d0
                 else if((t(k).lt.dp270_16) .and. (t(k).ge.dp268_16)) then
                    ni_secp   = 3.5d8*(dp270_16-t(k))/2.0d0*psacws(k)
                    nsacwi(k) = ni_secp
                    msacwi(k) = min(ni_secp*mi0,psacws(k))
                 else if((t(k).lt.dp268_16) .and. (t(k).ge.dp265_16)) then
                    ni_secp   = 3.5d8*(t(k)-dp265_16)/3.0d0*psacws(k)
                    nsacwi(k) = ni_secp
                    msacwi(k) = min(ni_secp*mi0,psacws(k))
                 else
                    ni_secp   = 0.0d0
                    nsacwi(k) = 0.0d0
                    msacwi(k) = 0.0d0
                 endif
                 psacws(k) = max(0.0d0,psacws(k)-ni_secp*mi0)
        
                 !.......................................................................
                 ! accretion of rain water by snow
                 ! formula from ikawa and saito, 1991, used by reisner et al., 1998
        
                 if (qric(k).ge.1.d-8 .and. qniic(k).ge.1.d-8 .and. & 
                      t(k).le.dp273_15) then
        
                    tmp1 = ((1.2d0*umr(k)-0.95d0*ums(k))*(1.2d0*umr(k)-0.95d0*ums(k))+ &
                         0.08d0*ums(k)*umr(k))
                    call math_agent_2i1o(power, tmp1, 0.5d0, tmp2)
                    pracs(k) = pi*pi*ecr*(tmp2*rhow*rho(k)* &
                         n0r(k)*n0s(k)* &
                         (5.d0/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k))+ &
                         2.d0/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k)*lams(k))+ &
                         0.5d0/(lamr(k)*lamr(k)*lamr(k)*lamr(k)*lams(k)*lams(k)*lams(k))))
                    !pracs(k) = pi*pi*ecr*(((1.2d0*umr(k)-0.95d0*ums(k))**2+ &
                    !     0.08d0*ums(k)*umr(k))**0.5d0*rhow*rho(k)* &
                    !     n0r(k)*n0s(k)* &
                    !     (5.d0/(lamr(k)**6*lams(k))+ &
                    !     2.d0/(lamr(k)**5*lams(k)**2)+ &
                    !     0.5d0/(lamr(k)**4*lams(k)**3)))
        
                    tmp1 = (1.7d0*(unr(k)-uns(k))*(unr(k)-uns(k))+ &
                         0.3d0*unr(k)*uns(k))
                    call math_agent_2i1o(power, tmp1, 0.5d0, tmp2)
                    npracs(k) = pi/2.d0*rho(k)*ecr*tmp2*n0r(k)*n0s(k)* &
                         (1.d0/(lamr(k)*lamr(k)*lamr(k)*lams(k))+ &
                         1.d0/(lamr(k)*lamr(k)*lams(k)*lams(k))+ &
                         1.d0/(lamr(k)*lams(k)*lams(k)*lams(k)))
                    !npracs(k) = pi/2.d0*rho(k)*ecr*(1.7d0*(unr(k)-uns(k))**2+ &
                    !     0.3d0*unr(k)*uns(k))**0.5d0*n0r(k)*n0s(k)* &
                    !     (1.d0/(lamr(k)**3*lams(k))+ &
                    !     1.d0/(lamr(k)**2*lams(k)**2)+ &
                    !     1.d0/(lamr(k)*lams(k)**3))
        
                 else
                    pracs(k)=0.d0
                    npracs(k)=0.d0
                 end if
        
                 !.......................................................................
                 ! heterogeneous freezing of rain drops
                 ! follows from Bigg (1953)
        
                 if (t(k).lt.dp269_15 .and. qric(k).ge.qsmall) then
        
                    tmp1 = (aimm*(dp273_15-t(k)))
                    call math_agent_1i1o(exp_agent, tmp1, tmp2)
                    mnuccr(k) = 20.d0*pi*pi*rhow*nric(k)*bimm* &
                         (tmp2-1.d0)/(lamr(k)*lamr(k)*lamr(k)) &
                         /(lamr(k)*lamr(k)*lamr(k))
        
                    nnuccr(k) = pi*nric(k)*bimm* &
                         (tmp2-1.d0)/(lamr(k)*lamr(k)*lamr(k))
                    !mnuccr(k) = 20.d0*pi*pi*rhow*nric(k)*bimm* &
                    !     (exp(aimm*(dp273_15-t(k)))-1.d0)/lamr(k)**3 &
                    !     /lamr(k)**3
        
                    !nnuccr(k) = pi*nric(k)*bimm* &
                    !     (exp(aimm*(dp273_15-t(k)))-1.d0)/lamr(k)**3
                 else
                    mnuccr(k)=0.d0
                    nnuccr(k)=0.d0
                 end if

         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+nsagg(k)
         !    !qsout2(k) = qsout2(k)+npsacws(k)
         !    !nrout2(k) = nrout2(k)+pracs(k)
         !    !nsout2(k) = nsout2(k)+npracs(k)
         !    !drout2(k) = drout2(k)+nnuccr(k)
         !    dsout2(k) = dsout2(k)+mnuccr(k)
         !endif
        
                 !.......................................................................
                 ! accretion of cloud liquid water by rain
                 ! formula from Khrouditnov and Kogan (2000)
                 ! gravitational collection kernel, droplet fall speed neglected
        
                 if (qric(k).ge.qsmall .and. qcic(k).ge.qsmall) then
        
                    ! include sub-grid distribution of cloud water
        
                    ! add sub-column switch
        
                    if (microp_uniform) then
        
                       tmp1 = (qcic(k)*qric(k))
                       call math_agent_2i1o(power, tmp1, dp1_15, tmp2)
                       pra(k) = 67.d0*tmp2
                       !pra(k) = 67.d0*(qcic(k)*qric(k))**dp1_15
                       npra(k) = pra(k)/(qcic(k)/ncic(k))
        
                    else
        
                       tmp1 = (qcic(k)*qric(k))
                       call math_agent_2i1o(power, tmp1, dp1_15, tmp2)
                       pra(k) = accre_enhan(k)*(cons12/(cons3*cons20)*67.d0*tmp2)
                       !pra(k) = accre_enhan(k)*(cons12/(cons3*cons20)*67.d0*(qcic(k)*qric(k))**dp1_15)
                       npra(k) = pra(k)/(qcic(k)/ncic(k))
        
                    end if               ! sub-column switch
        
                 else
                    pra(k)=0.d0
                    npra(k)=0.d0
                 end if
        
                 !.......................................................................
                 ! Self-collection of rain drops
                 ! from Beheng(1994)
        
                 if (qric(k).ge.qsmall) then
                    nragg(k) = -8.d0*nric(k)*qric(k)*rho(k)
                 else
                    nragg(k)=0.d0
                 end if
        
                 !.......................................................................
                 ! Accretion of cloud ice by snow
                 ! For this calculation, it is assumed that the Vs >> Vi
                 ! and Ds >> Di for continuous collection
        
                 if (do_cldice .and. qniic(k).ge.qsmall.and.qiic(k).ge.qsmall &
                      .and.t(k).le.dp273_15) then
        
                    tmp1 = bs+3.d0
                    call math_agent_2i1o(power, lams(k), tmp1, tmp2)

                    prai(k) = pi/4.d0*asn(k)*qiic(k)*rho(k)* &
                         n0s(k)*Eii*cons11/ &
                         tmp2
                    nprai(k) = pi/4.d0*asn(k)*niic(k)* &
                         rho(k)*n0s(k)*Eii*cons11/ &
                         tmp2
                    !prai(k) = pi/4.d0*asn(k)*qiic(k)*rho(k)* &
                    !     n0s(k)*Eii*cons11/ &
                    !     lams(k)**(bs+3.d0)
                    !nprai(k) = pi/4.d0*asn(k)*niic(k)* &
                    !     rho(k)*n0s(k)*Eii*cons11/ &
                    !     lams(k)**(bs+3.d0)
                 else
                    prai(k)=0.d0
                    nprai(k)=0.d0
                 end if
        
                 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 ! calculate evaporation/sublimation of rain and snow
                 ! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
                 ! in-cloud condensation/deposition of rain and snow is neglected
                 ! except for transfer of cloud water to snow through bergeron process
        
                 ! initialize evap/sub tendncies
                 pre(k)=0.d0
                 prds(k)=0.d0
        
                 ! evaporation of rain
                 ! only calculate if there is some precip fraction > cloud fraction
        
                 if (qcic(k)+qiic(k).lt.1.d-6.or.cldmax(k).gt.lcldm(k)) then
        
                    ! set temporary cloud fraction to zero if cloud water + ice is very small
                    ! this will ensure that evaporation/sublimation of precip occurs over
                    ! entire grid cell, since min cloud fraction is specified otherwise
                    if (qcic(k)+qiic(k).lt.1.d-6) then
                       dum=0.d0
                    else
                       dum=lcldm(k)
                    end if
        
                    ! saturation vapor pressure
                    esn=svp_water(t(k), lg10, power, svp_tboil)
                    qsn=svp_to_qsat(esn, p(k), svp_epsilo, svp_omeps)

          !tmp1 = 10.d0
          !tmp2 = svp_tboil/t(k)
          !call math_agent_1i1o(lg10, tmp2, tmp3)
          !tmp4 = 11.344d0*(1.d0-t(k)/svp_tboil)
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !tmp4 = (-3.49149d0*(svp_tboil/t(k)-1.d0))
          !call math_agent_2i1o(power, tmp1, tmp4, tmp5)
          !tmp4 = 1013.246d0
          !call math_agent_1i1o(lg10, tmp4, tmp1)
          !tmp4 = -7.90298d0*(svp_tboil/t(k)-1.d0)+ &
          !     5.02808d0*tmp3- &
          !     1.3816d-7*(tmp2-1.d0)+ &
          !     8.1328d-3*(tmp5-1.d0)+ &
          !     tmp1
          !tmp1 = 10.d0
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !esn = tmp2 * 100.d0

          !if ( (p(k) - esn) <= 0.d0 ) then
          !   qsn = 1.0d0
          !else
          !   qsn = svp_epsilo*esn / (p(k) - svp_omeps*esn)
          !end if
        
                    ! recalculate saturation vapor pressure for liquid and ice
                    esl(k)=esn
                    esi(k)=svp_ice(t(k), lg10, power, svp_h2otrip)

          !tmp1 = svp_h2otrip/t(k)
          !call math_agent_1i1o(lg10, tmp1, tmp2)
          !tmp1 = 6.1071d0
          !call math_agent_1i1o(lg10, tmp1, tmp3)
          !tmp1 = -9.09718d0*(svp_h2otrip/t(k)-1.d0)-3.56654d0* &
          !     tmp2+0.876793d0*(1.d0-t(k)/svp_h2otrip)+ &
          !     tmp3
          !tmp2 = 10.d0
          !call math_agent_2i1o(power, tmp2, tmp1, tmp3)
          !esi(k) = tmp3*100.d0

                    ! hm fix, make sure when above freezing that esi=esl, not active yet
                    if (t(k).gt.tmelt)esi(k)=esl(k)
        
                    ! calculate q for out-of-cloud region
                    qclr=(q(k)-dum*qsn)/(1.d0-dum)
        
                    if (qric(k).ge.qsmall) then
        
                       qvs=svp_to_qsat(esl(k), p(k), svp_epsilo, svp_omeps)

          !if ( (p(k) - esl(k)) <= 0.d0 ) then
          !   qvs = 1.0d0
          !else
          !   qvs = svp_epsilo*esl(k) / (p(k) - svp_omeps*esl(k))
          !end if

                       dqsdt = xxlv*qvs/(rv*t(k)*t(k))
                       !dqsdt = xxlv*qvs/(rv*t(k)**2)
                       ab = 1.d0+dqsdt*xxlv/cpp
                       tmp1 = (arn(k)*rho(k)/mu(k))
                       tmp3 = 1.d0/3.d0
                       call math_agent_2i1o(power, tmp1, 0.5d0, tmp2)
                       call math_agent_2i1o(power, sc(k), tmp3, tmp4)
                       tmp5 = (5.d0/2.d0+br/2.d0)
                       call math_agent_2i1o(power, lamr(k), tmp5, tmp6)

                       epsr = 2.d0*pi*n0r(k)*rho(k)*Dv(k)* &
                            (f1r/(lamr(k)*lamr(k))+ &
                            f2r*tmp2* &
                            tmp4*cons13/ &
                            (tmp6))
                       !epsr = 2.d0*pi*n0r(k)*rho(k)*Dv(k)* &
                       !     (f1r/(lamr(k)*lamr(k))+ &
                       !     f2r*(arn(k)*rho(k)/mu(k))**0.5d0* &
                       !     sc(k)**(1.d0/3.d0)*cons13/ &
                       !     (lamr(k)**(5.d0/2.d0+br/2.d0)))
        
                       pre(k) = epsr*(qclr-qvs)/ab
        
                       ! only evaporate in out-of-cloud region
                       ! and distribute across cldmax
                       pre(k)=min(pre(k)*(cldmax(k)-dum),0.d0)
                       pre(k)=pre(k)/cldmax(k)
                       am_evp_st(k) = max(cldmax(k)-dum, 0.d0)

         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+epsr
         !    qsout2(k) = qsout2(k)+qclr
         !    nrout2(k) = nrout2(k)+qvs
         !    nsout2(k) = nsout2(k)+qclr-qvs
         !    drout2(k) = drout2(k)+dum
         !    dsout2(k) = dsout2(k)+pre(k)
         !endif
                    end if
        
                    ! sublimation of snow
                    if (qniic(k).ge.qsmall) then
                       qvi=svp_to_qsat(esi(k), p(k), svp_epsilo, svp_omeps)

          !if ( (p(k) - esi(k)) <= 0.d0 ) then
          !   qvi = 1.0d0
          !else
          !   qvi = svp_epsilo*esi(k) / (p(k) - svp_omeps*esi(k))
          !end if
                       !dqsidt =  xxls*qvi/(rv*t(k)**2)
                       dqsidt =  xxls*qvi/(rv*t(k)*t(k))
                       abi = 1.d0+dqsidt*xxls/cpp
                       tmp1 = (asn(k)*rho(k)/mu(k))
                       tmp3 = 1.d0/3.d0
                       call math_agent_2i1o(power, tmp1, 0.5d0, tmp2)
                       call math_agent_2i1o(power, sc(k), tmp3, tmp4)
                       tmp5 = (5.d0/2.d0+bs/2.d0)
                       call math_agent_2i1o(power, lams(k), tmp5, tmp6)
                       epss = 2.d0*pi*n0s(k)*rho(k)*Dv(k)* &
                            (f1s/(lams(k)*lams(k))+ &
                            f2s*tmp2* &
                            tmp4*cons14/ &
                            (tmp6))
                       !epss = 2.d0*pi*n0s(k)*rho(k)*Dv(k)* &
                       !     (f1s/(lams(k)*lams(k))+ &
                       !     f2s*(asn(k)*rho(k)/mu(k))**0.5d0* &
                       !     sc(k)**(1.d0/3.d0)*cons14/ &
                       !     (lams(k)**(5.d0/2.d0+bs/2.d0)))
                       prds(k) = epss*(qclr-qvi)/abi
        
                       ! only sublimate in out-of-cloud region and distribute over cldmax
                       prds(k)=min(prds(k)*(cldmax(k)-dum),0.d0)
                       prds(k)=prds(k)/cldmax(k)
                       am_evp_st(k) = max(cldmax(k)-dum, 0.d0)
                    end if
        
                    ! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
                    ! get updated RH at end of time step based on cloud water/ice condensation/evap
        
                    qtmp=q(k)-(cmei(k)+(pre(k)+prds(k))*cldmax(k))*deltat
                    ttmp=t(k)+((pre(k)*cldmax(k))*xxlv+ &
                         (cmei(k)+prds(k)*cldmax(k))*xxls)*deltat/cpp
        
                    !limit range of temperatures!
                    ttmp=max(dp180_0,min(ttmp,dp323_0))
        
                    esn=svp_water(ttmp, lg10, power, svp_tboil)  ! use rhw to allow ice supersaturation
                    qsn=svp_to_qsat(esn, p(k), svp_epsilo, svp_omeps)

          !tmp1 = 10.d0
          !tmp2 = svp_tboil/ttmp
          !call math_agent_1i1o(lg10, tmp2, tmp3)
          !tmp4 = 11.344d0*(1.d0-ttmp/svp_tboil)
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !tmp4 = (-3.49149d0*(svp_tboil/ttmp-1.d0))
          !call math_agent_2i1o(power, tmp1, tmp4, tmp5)
          !tmp4 = 1013.246d0
          !call math_agent_1i1o(lg10, tmp4, tmp1)
          !tmp4 = -7.90298d0*(svp_tboil/ttmp-1.d0)+ &
          !     5.02808d0*tmp3- &
          !     1.3816d-7*(tmp2-1.d0)+ &
          !     8.1328d-3*(tmp5-1.d0)+ &
          !     tmp1
          !tmp1 = 10.d0
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !esn = tmp2 * 100.d0

          !if ( (p(k) - esn) <= 0.d0 ) then
          !   qsn = 1.0d0
          !else
          !   qsn = svp_epsilo*esn / (p(k) - svp_omeps*esn)
          !end if
        
                    ! modify precip evaporation rate if q > qsat
                    if (qtmp.gt.qsn) then
                       if (pre(k)+prds(k).lt.-1.d-20) then
                          dum1=pre(k)/(pre(k)+prds(k))
                          ! recalculate q and t after cloud water cond but without precip evap
                          qtmp=q(k)-(cmei(k))*deltat
                          ttmp=t(k)+(cmei(k)*xxls)*deltat/cpp
                          esn=svp_water(ttmp, lg10, power, svp_tboil) ! use rhw to allow ice supersaturation
                          qsn=svp_to_qsat(esn, p(k), svp_epsilo, svp_omeps)

          !tmp1 = 10.d0
          !tmp2 = svp_tboil/ttmp
          !call math_agent_1i1o(lg10, tmp2, tmp3)
          !tmp4 = 11.344d0*(1.d0-ttmp/svp_tboil)
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !tmp4 = (-3.49149d0*(svp_tboil/ttmp-1.d0))
          !call math_agent_2i1o(power, tmp1, tmp4, tmp5)
          !tmp4 = 1013.246d0
          !call math_agent_1i1o(lg10, tmp4, tmp1)
          !tmp4 = -7.90298d0*(svp_tboil/ttmp-1.d0)+ &
          !     5.02808d0*tmp3- &
          !     1.3816d-7*(tmp2-1.d0)+ &
          !     8.1328d-3*(tmp5-1.d0)+ &
          !     tmp1
          !tmp1 = 10.d0
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !esn = tmp2 * 100.d0

          !if ( (p(k) - esn) <= 0.d0 ) then
          !   qsn = 1.0d0
          !else
          !   qsn = svp_epsilo*esn / (p(k) - svp_omeps*esn)
          !end if

                          dum=(qtmp-qsn)/(1.d0 + cons27*qsn/(cpp*rv*ttmp*ttmp))
                          dum=min(dum,0.d0)
        
                          ! modify rates if needed, divide by cldmax to get local (in-precip) value
                          pre(k)=dum*dum1/deltat/cldmax(k)
        
                          ! do separately using RHI for prds....
                          esn=svp_ice(ttmp, lg10, power, svp_h2otrip) ! use rhi to allow ice supersaturation
                          qsn=svp_to_qsat(esn, p(k), svp_epsilo, svp_omeps)

          !tmp1 = svp_h2otrip/ttmp
          !call math_agent_1i1o(lg10, tmp1, tmp2)
          !tmp1 = 6.1071d0
          !call math_agent_1i1o(lg10, tmp1, tmp3)
          !tmp1 = -9.09718d0*(svp_h2otrip/ttmp-1.d0)-3.56654d0* &
          !     tmp2+0.876793d0*(1.d0-ttmp/svp_h2otrip)+ &
          !     tmp3
          !tmp2 = 10.d0
          !call math_agent_2i1o(power, tmp2, tmp1, tmp3)
          !esn = tmp3*100.d0

          !if ( (p(k) - esn) <= 0.d0 ) then
          !   qsn = 1.0d0
          !else
          !   qsn = svp_epsilo*esn / (p(k) - svp_omeps*esn)
          !end if

                          dum=(qtmp-qsn)/(1.d0 + cons28*qsn/(cpp*rv*ttmp*ttmp))
                          dum=min(dum,0.d0)
        
                          ! modify rates if needed, divide by cldmax to get local (in-precip) value
                          prds(k)=dum*(1.d0-dum1)/deltat/cldmax(k)
                       end if
                    end if
                 end if
        
                 ! bergeron process - evaporation of droplets and deposition onto snow
        
                 if (qniic(k).ge.qsmall.and.qcic(k).ge.qsmall.and.t(k).lt.tmelt) then
                    qvi=svp_to_qsat(esi(k), p(k), svp_epsilo, svp_omeps)
                    qvs=svp_to_qsat(esl(k), p(k), svp_epsilo, svp_omeps)

          !if ( (p(k) - esi(k)) <= 0.d0 ) then
          !   qvi = 1.0d0
          !else
          !   qvi = svp_epsilo*esi(k) / (p(k) - svp_omeps*esi(k))
          !end if

          !if ( (p(k) - esl(k)) <= 0.d0 ) then
          !   qvs = 1.0d0
          !else
          !   qvs = svp_epsilo*esl(k) / (p(k) - svp_omeps*esl(k))
          !end if

                    dqsidt =  xxls*qvi/(rv*t(k)*t(k))
                    abi = 1.d0+dqsidt*xxls/cpp

                    tmp1 = (asn(k)*rho(k)/mu(k))
                    tmp3 = 1.d0/3.d0
                    call math_agent_2i1o(power, tmp1, 0.5d0, tmp2)
                    call math_agent_2i1o(power, sc(k), tmp3, tmp4)
                    tmp5 = (5.d0/2.d0+bs/2.d0)
                    call math_agent_2i1o(power, lams(k), tmp5, tmp6)
                    epss = 2.d0*pi*n0s(k)*rho(k)*Dv(k)* &
                         (f1s/(lams(k)*lams(k))+ &
                         f2s*tmp2* &
                         tmp4*cons14/ &
                         (tmp6))
                    !epss = 2.d0*pi*n0s(k)*rho(k)*Dv(k)* &
                    !     (f1s/(lams(k)*lams(k))+ &
                    !     f2s*(asn(k)*rho(k)/mu(k))**0.5d0* &
                    !     sc(k)**(1.d0/3.d0)*cons14/ &
                    !     (lams(k)**(5.d0/2.d0+bs/2.d0)))
                    bergs(k)=epss*(qvs-qvi)/abi
                 else
                    bergs(k)=0.d0
                 end if
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+npra(k)
         !    qsout2(k) = qsout2(k)+nprai(k)
         !    nrout2(k) = nrout2(k)+prai(k)
         !    nsout2(k) = nsout2(k)+pre(k)
         !    drout2(k) = drout2(k)+prds(k)
         !    dsout2(k) = dsout2(k)+bergs(k)
         !endif
                 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 ! conservation to ensure no negative values of cloud water/precipitation
                 ! in case microphysical process rates are large
        
                 ! make sure and use end-of-time step values for cloud water, ice, due
                 ! condensation/deposition
        
                 ! note: for check on conservation, processes are multiplied by omsm
                 ! to prevent problems due to round off error
        
                 ! include mixing timescale  (mtime)
        
                 qce=(qc(k) - berg(k)*deltat)
                 nce=(nc(k)+npccn(k)*deltat*mtime)
                 qie=(qi(k)+(cmei(k)+berg(k))*deltat)
                 nie=(ni(k)+nnuccd(k)*deltat*mtime)
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+qce
         !    qsout2(k) = qsout2(k)+nce
         !    nrout2(k) = nrout2(k)+qie
         !    nsout2(k) = nsout2(k)+nie
         !    !drout2(k) = drout2(k)+prds(k)
         !    !dsout2(k) = dsout2(k)+bergs(k)
         !endif
                 ! conservation of qc
        
                 dum = (prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+ &
                      psacws(k)+bergs(k))*lcldm(k)*deltat

         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+prc(k)
         !    qsout2(k) = qsout2(k)+pra(k)
         !    nrout2(k) = nrout2(k)+mnuccc(k)
         !    nsout2(k) = nsout2(k)+mnucct(k)
         !    drout2(k) = drout2(k)+msacwi(k)
         !    dsout2(k) = dsout2(k)+psacws(k)
         !endif
        
                 if (dum.gt.qce) then
                    ratio = qce/deltat/lcldm(k)/(prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+psacws(k)+bergs(k))*omsm 
        
                    prc(k) = prc(k)*ratio
                    pra(k) = pra(k)*ratio
                    mnuccc(k) = mnuccc(k)*ratio
                    mnucct(k) = mnucct(k)*ratio  
                    msacwi(k) = msacwi(k)*ratio  
                    psacws(k) = psacws(k)*ratio
                    bergs(k) = bergs(k)*ratio
                 end if
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+prc(k)
         !    qsout2(k) = qsout2(k)+pra(k)
         !    nrout2(k) = nrout2(k)+mnuccc(k)
         !    nsout2(k) = nsout2(k)+mnucct(k)
         !    drout2(k) = drout2(k)+msacwi(k)
         !    dsout2(k) = dsout2(k)+psacws(k)
         !endif
                 ! conservation of nc
        
                 dum = (nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
                      npsacws(k)-nsubc(k))*lcldm(k)*deltat
        
                 if (dum.gt.nce) then
                    ratio = nce/deltat/((nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+&
                         npsacws(k)-nsubc(k))*lcldm(k))*omsm
        
                    nprc1(k) = nprc1(k)*ratio
                    npra(k) = npra(k)*ratio
                    nnuccc(k) = nnuccc(k)*ratio
                    nnucct(k) = nnucct(k)*ratio  
                    npsacws(k) = npsacws(k)*ratio
                    nsubc(k)=nsubc(k)*ratio
                 end if
        
                 ! conservation of qi
        
                 if (do_cldice) then
        
                    frztmp = -mnuccc(k) - mnucct(k) - msacwi(k)
                    if (use_hetfrz_classnuc) frztmp = -mnuccc(k)-mnucct(k)-mnudep(k)-msacwi(k)
                    dum = ( frztmp*lcldm(k) + (prci(k)+prai(k))*icldm(k) )*deltat
        
                    if (dum.gt.qie) then
        
                       frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
                       if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
                       ratio = (qie/deltat + frztmp*lcldm(k))/((prci(k)+prai(k))*icldm(k))*omsm 
                       prci(k) = prci(k)*ratio
                       prai(k) = prai(k)*ratio
                    end if
        
                    ! conservation of ni
                    frztmp = -nnucct(k) - nsacwi(k)
                    if (use_hetfrz_classnuc) frztmp = -nnucct(k) - nnuccc(k) - nnudep(k) - nsacwi(k)
                    dum = ( frztmp*lcldm(k) + (nprci(k)+nprai(k)-nsubi(k))*icldm(k) )*deltat
        
                    if (dum.gt.nie) then
        
                       frztmp = nnucct(k) + nsacwi(k)
                       if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
                       ratio = (nie/deltat + frztmp*lcldm(k))/ &  
                             ((nprci(k)+nprai(k)-nsubi(k))*icldm(k))*omsm
                       nprci(k) = nprci(k)*ratio
                       nprai(k) = nprai(k)*ratio
                       nsubi(k) = nsubi(k)*ratio
                    end if
                 end if

         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+ratio
         !    qsout2(k) = qsout2(k)+qnitend(k)
         !    nrout2(k) = nrout2(k)+qrtot
         !    nsout2(k) = nsout2(k)+nrtot
         !    drout2(k) = drout2(k)+qstot
         !    dsout2(k) = dsout2(k)+nstot
         !endif
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+frztmp
         !    qsout2(k) = qsout2(k)+mnucct(k)
         !    nrout2(k) = nrout2(k)+nnucct(k)
         !    nsout2(k) = nsout2(k)+prci(k)
         !    drout2(k) = drout2(k)+prai(k)
         !    dsout2(k) = dsout2(k)+nsubi(k)
         !endif
        
                 ! for precipitation conservation, use logic that vertical integral 
                 ! of tendency from current level to top of model (i.e., qrtot) cannot be negative
        
                 ! conservation of rain mixing rat
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+pre(k)
         !    qsout2(k) = qsout2(k)+pracs(k)
         !    nrout2(k) = nrout2(k)+mnuccr(k)
         !endif
                 if (((prc(k)+pra(k))*lcldm(k)+(-mnuccr(k)+pre(k)-pracs(k))*&
                      cldmax(k))*dz(k)*rho(k)+qrtot.lt.0.d0) then
        
                    if (-pre(k)+pracs(k)+mnuccr(k).ge.qsmall) then
        
                       ratio = (qrtot/(dz(k)*rho(k))+(prc(k)+pra(k))*lcldm(k))/&
                            ((-pre(k)+pracs(k)+mnuccr(k))*cldmax(k))*omsm 
        
                       pre(k) = pre(k)*ratio
                       pracs(k) = pracs(k)*ratio
                       mnuccr(k) = mnuccr(k)*ratio
                    end if
                 end if
         !if ( it .eq. 1 ) then
         !    nsout2(k) = nsout2(k)+ qrtot
         !    drout2(k) = drout2(k)+ pracs(k)
         !    dsout2(k) = dsout2(k)+ mnuccr(k)
         !endif
        
                 ! conservation of nr
                 ! for now neglect evaporation of nr
                 nsubr(k)=0.d0
        
                 if ((nprc(k)*lcldm(k)+(-nnuccr(k)+nsubr(k)-npracs(k)&
                      +nragg(k))*cldmax(k))*dz(k)*rho(k)+nrtot.lt.0.d0) then
        
                    if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).ge.qsmall) then
        
                       ratio = (nrtot/(dz(k)*rho(k))+nprc(k)*lcldm(k))/&
                            ((-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k))*cldmax(k))*omsm
        
                       nsubr(k) = nsubr(k)*ratio
                       npracs(k) = npracs(k)*ratio
                       nnuccr(k) = nnuccr(k)*ratio
                       nragg(k) = nragg(k)*ratio
                    end if
                 end if
        
                 ! conservation of snow mix ratio
        
                 if (((bergs(k)+psacws(k))*lcldm(k)+(prai(k)+prci(k))*icldm(k)+(pracs(k)+&
                      mnuccr(k)+prds(k))*cldmax(k))*dz(k)*rho(k)+qstot.lt.0.d0) then
        
                    if (-prds(k).ge.qsmall) then
        
                       ratio = (qstot/(dz(k)*rho(k))+(bergs(k)+psacws(k))*lcldm(k)+(prai(k)+prci(k))*icldm(k)+&
                            (pracs(k)+mnuccr(k))*cldmax(k))/(-prds(k)*cldmax(k))*omsm
        
                       prds(k) = prds(k)*ratio
                    end if
                 end if
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+ratio
         !    qsout2(k) = qsout2(k)+prds(k)
         !    nrout2(k) = nrout2(k)+qstot
         !    nsout2(k) = nsout2(k)+psacws(k)
         !    drout2(k) = drout2(k)+bergs(k)
         !    dsout2(k) = dsout2(k)+dz(k)
         !endif
        
                 ! conservation of ns
        
                 ! calculate loss of number due to sublimation
                 ! for now neglect sublimation of ns
                 nsubs(k)=0.d0
        
                 if ((nprci(k)*icldm(k)+(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax(k))*&
                      dz(k)*rho(k)+nstot.lt.0.d0) then
        
                    if (-nsubs(k)-nsagg(k).ge.qsmall) then
        
                       ratio = (nstot/(dz(k)*rho(k))+nprci(k)*icldm(k)+&
                            nnuccr(k)*cldmax(k))/((-nsubs(k)-nsagg(k))*cldmax(k))*omsm
        
                       nsubs(k) = nsubs(k)*ratio
                       nsagg(k) = nsagg(k)*ratio
                    end if
                 end if
        
                 ! get tendencies due to microphysical conversion processes
                 ! note: tendencies are multiplied by appropaiate cloud/precip 
                 ! fraction to get grid-scale values
                 ! note: cmei is already grid-average values
        
                 qvlat(k) = qvlat(k)-(pre(k)+prds(k))*cldmax(k)-cmei(k) 
        
                 tlat(k) = tlat(k)+((pre(k)*cldmax(k)) &
                      *xxlv+(prds(k)*cldmax(k)+cmei(k))*xxls+ &
                      ((bergs(k)+psacws(k)+mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(k)+(mnuccr(k)+ &
                      pracs(k))*cldmax(k)+berg(k))*xlf)
        
                 qctend(k) = qctend(k)+ &
                      (-pra(k)-prc(k)-mnuccc(k)-mnucct(k)-msacwi(k)- & 
                      psacws(k)-bergs(k))*lcldm(k)-berg(k)
        
                 if (do_cldice) then
        
                    frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
                    if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
                    qitend(k) = qitend(k) + frztmp*lcldm(k) + &
                       (-prci(k)-prai(k))*icldm(k) + cmei(k) + berg(k)
        
                 end if
        
                 qrtend(k) = qrtend(k)+ &
                      (pra(k)+prc(k))*lcldm(k)+(pre(k)-pracs(k)- &
                      mnuccr(k))*cldmax(k)
        
                 qnitend(k) = qnitend(k)+ &
                      (prai(k)+prci(k))*icldm(k)+(psacws(k)+bergs(k))*lcldm(k)+(prds(k)+ &
                      pracs(k)+mnuccr(k))*cldmax(k)
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+qnitend(k)
         !    qsout2(k) = qsout2(k)+icldm(k)
         !    nrout2(k) = nrout2(k)+lcldm(k)
         !    nsout2(k) = nsout2(k)+mnuccr(k)
         !    drout2(k) = drout2(k)+psacws(k)
         !    dsout2(k) = dsout2(k)+cldmax(k)
         !    !qsout2(k) = qsout2(k)+prai(k)
         !    !nrout2(k) = nrout2(k)+prci(k)
         !    !nsout2(k) = nsout2(k)+bergs(k)
         !    !drout2(k) = drout2(k)+prds(k)
         !    !dsout2(k) = dsout2(k)+pracs(k)
         !endif
        
                 ! add output for cmei (accumulate)
                 cmeiout(k) = cmeiout(k) + cmei(k)
        
                 ! assign variables for trop_mozart, these are grid-average
                 ! evaporation/sublimation is stored here as positive term
        
                 evapsnow(k) = evapsnow(k)-prds(k)*cldmax(k)
                 nevapr(k) = nevapr(k)-pre(k)*cldmax(k)
                 nevapr2(k) = nevapr2(k)-pre(k)*cldmax(k)
        
                 ! change to make sure prain is positive: do not remove snow from
                 ! prain used for wet deposition
                 prain(k) = prain(k)+(pra(k)+prc(k))*lcldm(k)+(-pracs(k)- &
                      mnuccr(k))*cldmax(k)
                 prodsnow(k) = prodsnow(k)+(prai(k)+prci(k))*icldm(k)+(psacws(k)+bergs(k))*lcldm(k)+(&
                      pracs(k)+mnuccr(k))*cldmax(k)
        
                 ! following are used to calculate 1st order conversion rate of cloud water
                 !    to rain and snow (1/s), for later use in aerosol wet removal routine
                 ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
                 !    used to calculate pra, prc, ... in this routine
                 ! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
                 !                      (no cloud ice or bergeron terms)
                 ! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }
        
                 qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k) + (pra(k)+prc(k)+psacws(k))*lcldm(k) 
                 qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc(k) 
        
                 ! microphysics output, note this is grid-averaged
                 prao(k)=prao(k)+pra(k)*lcldm(k)
                 prco(k)=prco(k)+prc(k)*lcldm(k)
                 mnuccco(k)=mnuccco(k)+mnuccc(k)*lcldm(k)
                 mnuccto(k)=mnuccto(k)+mnucct(k)*lcldm(k)
                 mnuccdo(k)=mnuccdo(k)+mnuccd(k)*lcldm(k)
                 msacwio(k)=msacwio(k)+msacwi(k)*lcldm(k)
                 psacwso(k)=psacwso(k)+psacws(k)*lcldm(k)
                 bergso(k)=bergso(k)+bergs(k)*lcldm(k)
                 bergo(k)=bergo(k)+berg(k)
                 prcio(k)=prcio(k)+prci(k)*icldm(k)
                 praio(k)=praio(k)+prai(k)*icldm(k)
                 mnuccro(k)=mnuccro(k)+mnuccr(k)*cldmax(k)
                 pracso (k)=pracso (k)+pracs (k)*cldmax(k)
        
                 ! multiply activation/nucleation by mtime to account for fast timescale
        
                 nctend(k) = nctend(k)+ npccn(k)*mtime+&
                      (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k) & 
                      -npra(k)-nprc1(k))*lcldm(k)      
        
                 if (do_cldice) then
        
                    frztmp = nnucct(k) + nsacwi(k)
                    if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
                    nitend(k) = nitend(k) + nnuccd(k)*mtime + & 
                          frztmp*lcldm(k) + (nsubi(k)-nprci(k)-nprai(k))*icldm(k)
        
                 end if
        
                 nstend(k) = nstend(k)+(nsubs(k)+ &
                      nsagg(k)+nnuccr(k))*cldmax(k)+nprci(k)*icldm(k)
        
                 nrtend(k) = nrtend(k)+ &
                      nprc(k)*lcldm(k)+(nsubr(k)-npracs(k)-nnuccr(k) &
                      +nragg(k))*cldmax(k)

         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+nprc(k)*lcldm(k)
         !    qsout2(k) = qsout2(k)+(nsubr(k)-npracs(k)-nnuccr(k) &
         !     +nragg(k))*cldmax(k)
         !    nrout2(k) = nrout2(k)+nprc(k)*lcldm(k)+(nsubr(k)-npracs(k)-nnuccr(k) &
         !     +nragg(k))*cldmax(k)
         !    nsout2(k) = nsout2(k)+nrtend(k)
         !    !drout2(i,k) = drout2(i,k)+nsubr(k)-npracs(k)-nnuccr(k)
         !    !dsout2(i,k) = dsout2(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
         !    ! +nragg(k))*cldmax(i,k)
         !endif

         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+nprc(k)
         !    qsout2(k) = qsout2(k)+nsubr(k)
         !    nrout2(k) = nrout2(k)+npracs(k)
         !    nsout2(k) = nsout2(k)+nnuccr(k)
         !    drout2(k) = drout2(k)+nragg(k)
         !    dsout2(k) = dsout2(k)+nrtend(k)
         !endif
        
                 ! make sure that nc and ni at advanced time step do not exceed
                 ! maximum (existing N + source terms*dt), which is possible due to
                 ! fast nucleation timescale
        
                 if (nctend(k).gt.0.d0.and.nc(k)+nctend(k)*deltat.gt.ncmax) then
                    nctend(k)=max(0.d0,(ncmax-nc(k))/deltat)
                 end if
        
                 if (do_cldice .and. nitend(k).gt.0.d0.and.ni(k)+nitend(k)*deltat.gt.nimax) then
                    nitend(k)=max(0.d0,(nimax-ni(k))/deltat)
                 end if
        
         !if ( it .eq. 1 ) then
         !    qrout2(k) = qrout2(k)+qrtend(k)
         !    qsout2(k) = qsout2(k)+qnitend(k)
         !    nrout2(k) = nrout2(k)+nstend(k)
         !    nsout2(k) = nsout2(k)+nrtend(k)
         !    drout2(k) = drout2(k)+nctend(k)
         !    dsout2(k) = dsout2(k)+nitend(k)
         !endif
                 ! get final values for precipitation q and N, based on
                 ! flux of precip from above, source/sink term, and terminal fallspeed
                 ! see eq. 15-16 in MG2008
        
                 ! rain
        
                 if (qric(k).ge.qsmall) then
                     if ( myid < 32 ) then
                        if (k.eq.top_lev) then
                           qric(k)=qrtend(k)*dz(k)/cldmax(k)/umr(k)
                           nric(k)=nrtend(k)*dz(k)/cldmax(k)/unr(k)
                        else
                           qric(k) = (rho(k-1)*umr(k-1)*qric(k-1)*cldmax(k-1)+ &
                                (rho(k)*dz(k)*qrtend(k)))/(umr(k)*rho(k)*cldmax(k))
                           nric(k) = (rho(k-1)*unr(k-1)*nric(k-1)*cldmax(k-1)+ &
                                (rho(k)*dz(k)*nrtend(k)))/(unr(k)*rho(k)*cldmax(k))
        
                        end if
                    else 
                        if (k.eq.top_lev) then
                           qric(k) = (tmpn(4)*tmpn(1)*tmpn(15)*tmpn(5)+ &
                                (rho(k)*dz(k)*qrtend(k)))/(umr(k)*rho(k)*cldmax(k))
                           nric(k) = (tmpn(4)*tmpn(2)*tmpn(10)*tmpn(5)+ &
                                (rho(k)*dz(k)*nrtend(k)))/(unr(k)*rho(k)*cldmax(k))
                        else
                           qric(k) = (rho(k-1)*umr(k-1)*qric(k-1)*cldmax(k-1)+ &
                                (rho(k)*dz(k)*qrtend(k)))/(umr(k)*rho(k)*cldmax(k))
                           nric(k) = (rho(k-1)*unr(k-1)*nric(k-1)*cldmax(k-1)+ &
                                (rho(k)*dz(k)*nrtend(k)))/(unr(k)*rho(k)*cldmax(k))
        
                        end if
                    endif
                 else
                    qric(k)=0.d0
                    nric(k)=0.d0
                 end if
        
         !if ( it .eq. 1 ) then
         !    nrout2(k) = nrout2(k) + qric(k)
         !endif
                 ! snow
        
                 if (qniic(k).ge.qsmall) then
                     if ( myid < 32 ) then
                        if (k.eq.top_lev) then
                           qniic(k)=qnitend(k)*dz(k)/cldmax(k)/ums(k)
                           nsic(k)=nstend(k)*dz(k)/cldmax(k)/uns(k)
                        else
                           qniic(k) = (rho(k-1)*ums(k-1)*qniic(k-1)*cldmax(k-1)+ &
                                (rho(k)*dz(k)*qnitend(k)))/(ums(k)*rho(k)*cldmax(k))
                           nsic(k) = (rho(k-1)*uns(k-1)*nsic(k-1)*cldmax(k-1)+ &
                                (rho(k)*dz(k)*nstend(k)))/(uns(k)*rho(k)*cldmax(k))
                        end if
                    else 
                        if (k.eq.top_lev) then
                           qniic(k) = (tmpn(4)*tmpn(8)*tmpn(3)*tmpn(5)+ &
                                (rho(k)*dz(k)*qnitend(k)))/(ums(k)*rho(k)*cldmax(k))
                           nsic(k) = (tmpn(4)*tmpn(17)*tmpn(22)*tmpn(5)+ &
                                (rho(k)*dz(k)*nstend(k)))/(uns(k)*rho(k)*cldmax(k))
                        else
                           qniic(k) = (rho(k-1)*ums(k-1)*qniic(k-1)*cldmax(k-1)+ &
                                (rho(k)*dz(k)*qnitend(k)))/(ums(k)*rho(k)*cldmax(k))
                           nsic(k) = (rho(k-1)*uns(k-1)*nsic(k-1)*cldmax(k-1)+ &
                                (rho(k)*dz(k)*nstend(k)))/(uns(k)*rho(k)*cldmax(k))
                        end if
                    endif
                 else
                    qniic(k)=0.d0
                    nsic(k)=0.d0
                 end if
        !if ( it .eq. 1  ) then
        !    qrout2(k) = qrout2(k)+qniic(k)
        !    qsout2(k) = qsout2(k)+nsic(k)
        !    !nrout2(k) = nrout2(k)+qc(k)
        !    !nsout2(k) = nsout2(k)+qi(k)
        !    !drout2(k) = drout2(k)+nc(k)
        !    !dsout2(k) = dsout2(k)+ni(k)
        !endif
        
                 ! calculate precipitation flux at surface
                 ! divide by density of water to get units of m/s
        
                    !qrout2(k) = qrtend(k)
                    !qsout2(k) = qnitend(k)
                    !nrout2(k) = prect
                    !nsout2(k) = (qrtend(k)*dz(k)*rho(k)+qnitend(k)*dz(k)*rho(k))/rhow
                    !drout2(k) = qrtend(k)*dz(k)*rho(k)
                    !dsout2(k) = qnitend(k)*dz(k)*rho(k)

                 prect = prect+(qrtend(k)*dz(k)*rho(k)+&
                      qnitend(k)*dz(k)*rho(k))/rhow
                 preci = preci+qnitend(k)*dz(k)*rho(k)/rhow
        
                 ! convert rain rate from m/s to mm/hr
        
                 rainrt(k)=qric(k)*rho(k)*umr(k)/rhow*3600.d0*1000.d0
        
                 ! vertically-integrated precip source/sink terms (note: grid-averaged)
        
                 qrtot = max(qrtot+qrtend(k)*dz(k)*rho(k),0.d0)
                 qstot = max(qstot+qnitend(k)*dz(k)*rho(k),0.d0)
                 nrtot = max(nrtot+nrtend(k)*dz(k)*rho(k),0.d0)
                 nstot = max(nstot+nstend(k)*dz(k)*rho(k),0.d0)
        
                 ! calculate melting and freezing of precip
        
                 ! melt snow at +2 C
        
                 if (t(k)+tlat(k)/cpp*deltat > 275.15d0) then
                    if (qstot > 0.d0) then
        
                       ! make sure melting snow doesn't reduce temperature below threshold
                       dum = -xlf/cpp*qstot/(dz(k)*rho(k))
                       if (t(k)+tlat(k)/cpp*deltat+dum.lt.275.15d0) then
                          dum = (t(k)+tlat(k)/cpp*deltat-275.15d0)*cpp/xlf
                          dum = dum/(xlf/cpp*qstot/(dz(k)*rho(k)))
                          dum = max(0.d0,dum)
                          dum = min(1.d0,dum)
                       else
                          dum = 1.d0
                       end if
        
                       qric(k)=qric(k)+dum*qniic(k)
                       nric(k)=nric(k)+dum*nsic(k)
                       qniic(k)=(1.d0-dum)*qniic(k)
                       nsic(k)=(1.d0-dum)*nsic(k)
                       ! heating tendency 
                       tmp=-xlf*dum*qstot/(dz(k)*rho(k))
                       meltsdt(k)=meltsdt(k) + tmp
        
                       tlat(k)=tlat(k)+tmp
                       qrtot=qrtot+dum*qstot
                       nrtot=nrtot+dum*nstot
                       qstot=(1.d0-dum)*qstot
                       nstot=(1.d0-dum)*nstot
                       preci=(1.d0-dum)*preci
                    end if
                 end if
        
         !if ( it .eq. 1 ) then
         !    nsout2(k) = nsout2(k) + qric(k)
         !endif
                 ! freeze all rain at -5C for Arctic
        
                 if (t(k)+tlat(k)/cpp*deltat < (tmelt - 5.d0)) then
        
                    if (qrtot > 0.d0) then
        
                       ! make sure freezing rain doesn't increase temperature above threshold
                       dum = xlf/cpp*qrtot/(dz(k)*rho(k))
                       if (t(k)+tlat(k)/cpp*deltat+dum.gt.(tmelt - 5.d0)) then
                          dum = -(t(k)+tlat(k)/cpp*deltat-(tmelt-5.d0))*cpp/xlf
                          dum = dum/(xlf/cpp*qrtot/(dz(k)*rho(k)))
                          dum = max(0.d0,dum)
                          dum = min(1.d0,dum)
                       else
                          dum = 1.d0
                       end if
        
                       qniic(k)=qniic(k)+dum*qric(k)
                       nsic(k)=nsic(k)+dum*nric(k)
                       qric(k)=(1.d0-dum)*qric(k)
                       nric(k)=(1.d0-dum)*nric(k)
                       ! heating tendency 
                       tmp = xlf*dum*qrtot/(dz(k)*rho(k))
                       frzrdt(k)=frzrdt(k) + tmp
        
                       tlat(k)=tlat(k)+tmp
                       qstot=qstot+dum*qrtot
                       qrtot=(1.d0-dum)*qrtot
                       nstot=nstot+dum*nrtot
                       nrtot=(1.d0-dum)*nrtot
                       preci=preci+dum*(prect-preci)
                    end if
                 end if
        
                 ! if rain/snow mix ratio is zero so should number concentration
        
                 if (qniic(k).lt.qsmall) then
                    qniic(k)=0.d0
                    nsic(k)=0.d0
                 end if
        
                 if (qric(k).lt.qsmall) then
                    qric(k)=0.d0
                    nric(k)=0.d0
                 end if
        
         !if ( it .eq. 1 ) then
         !    drout2(k) = drout2(k) + qric(k)
         !endif
                 ! make sure number concentration is a positive number to avoid 
                 ! taking root of negative
        
                 nric(k)=max(nric(k),0.d0)
                 nsic(k)=max(nsic(k),0.d0)
        
                 !.......................................................................
                 ! get size distribution parameters for fallspeed calculations
                 !......................................................................
                 ! rain
        
                 if (qric(k).ge.qsmall) then
                    tmp1 = (pi*rhow*nric(k)/qric(k))
                    tmp2 = 1.d0/3.d0
                    call math_agent_2i1o(power, tmp1, tmp2, lamr(k))
                    !lamr(k) = (pi*rhow*nric(k)/qric(k))**(1.d0/3.d0)
                    n0r(k) = nric(k)*lamr(k)
        
                    ! check for slope
                    ! change lammax and lammin for rain and snow
                    ! adjust vars
        
                    if (lamr(k).lt.lamminr) then
        
                       lamr(k) = lamminr
        
                       !n0r(k) = lamr(k)**4*qric(k)/(pi*rhow)
                       n0r(k) = lamr(k)*lamr(k)*lamr(k)*lamr(k)*qric(k)/(pi*rhow)
                       nric(k) = n0r(k)/lamr(k)
                    else if (lamr(k).gt.lammaxr) then
                       lamr(k) = lammaxr
                       n0r(k) = lamr(k)*lamr(k)*lamr(k)*lamr(k)*qric(k)/(pi*rhow)
                       !n0r(k) = lamr(k)**4*qric(k)/(pi*rhow)
                       nric(k) = n0r(k)/lamr(k)
                    end if
        
        
                    ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
        
                    call math_agent_2i1o(power, lamr(k), br, tmp1)
                    unr(k) = min(arn(k)*cons4/tmp1,dp9_1*rhof(k))
                    umr(k) = min(arn(k)*cons5/(6.d0*tmp1),dp9_1*rhof(k))
                    !unr(k) = min(arn(k)*cons4/lamr(k)**br,dp9_1*rhof(k))
                    !umr(k) = min(arn(k)*cons5/(6.d0*lamr(k)**br),dp9_1*rhof(k))
        
                 else
                    lamr(k) = 0.d0
                    n0r(k) = 0.d0
                    umr(k)=0.d0
                    unr(k)=0.d0
                 end if
        
                 !calculate mean size of combined rain and snow
        
                 if (lamr(k).gt.0.d0) then
                    !Artmp = n0r(k) * pi / (2.d0 * lamr(k)**3.d0)
                    Artmp = n0r(k) * pi / (2.d0 * lamr(k)*lamr(k)*lamr(k))
                 else 
                    Artmp = 0.d0
                 endif
        
                 if (lamc(k).gt.0.d0) then
                    !Actmp = cdist1(k) * pi * shr_spfn_gamma_nonintrinsic_r8(pgam(k)+3.d0)/(4.d0 * lamc(k)**2.d0)
                    Actmp = cdist1(k) * pi * shr_spfn_gamma_nonintrinsic_r8(pgam(k)+3.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/(4.d0 * lamc(k)*lamc(k))
                 else 
                    Actmp = 0.d0
                 endif
        
                 if (Actmp.gt.0d0.or.Artmp.gt.0) then
                    rercld(k)=rercld(k) + 3.d0 *(qric(k) + qcic(k)) / (4.d0 * rhow * (Actmp + Artmp))
                    arcld(k)=arcld(k)+1.d0
                 endif
        
                 !......................................................................
                 ! snow
        
                 if (qniic(k).ge.qsmall) then
                    tmp1 = (cons6*cs*nsic(k)/ &
                         qniic(k))
                    tmp2 = 1.d0/ds
                    call math_agent_2i1o(power, tmp1, tmp2, lams(k))
                    !lams(k) = (cons6*cs*nsic(k)/ &
                    !     qniic(k))**(1.d0/ds)
                    n0s(k) = nsic(k)*lams(k)
        
                    ! check for slope
                    ! adjust vars
        
                    if (lams(k).lt.lammins) then
                       lams(k) = lammins
                       tmp1 = ds+1.d0
                       call math_agent_2i1o(power, lams(k), tmp1, tmp2)
                       n0s(k) = tmp2*qniic(k)/(cs*cons6)
                       !n0s(k) = lams(k)**(ds+1.d0)*qniic(k)/(cs*cons6)
                       nsic(k) = n0s(k)/lams(k)
        
                    else if (lams(k).gt.lammaxs) then
                       lams(k) = lammaxs
                       tmp1 = ds+1.d0
                       call math_agent_2i1o(power, lams(k), tmp1, tmp2)
                       n0s(k) = tmp2*qniic(k)/(cs*cons6)
                       !n0s(k) = lams(k)**(ds+1.d0)*qniic(k)/(cs*cons6)
                       nsic(k) = n0s(k)/lams(k)
                    end if
        
                    ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
        
                    call math_agent_2i1o(power, lams(k), bs, tmp1)
                    ums(k) = min(asn(k)*cons8/(6.d0*tmp1),1.2d0*rhof(k))
                    uns(k) = min(asn(k)*cons7/tmp1,1.2d0*rhof(k))
                    !ums(k) = min(asn(k)*cons8/(6.d0*lams(k)**bs),1.2d0*rhof(k))
                    !uns(k) = min(asn(k)*cons7/lams(k)**bs,1.2d0*rhof(k))
        
                 else
                    lams(k) = 0.d0
                    n0s(k) = 0.d0
                    ums(k) = 0.d0
                    uns(k) = 0.d0
                 end if

                 !if ( k .eq. pver .and. myid < 32 ) then
                 !    tmpn(1 ) = umr   (pver)
                 !    tmpn(2 ) = unr   (pver)
                 !    tmpn(3 ) = qniic (pver)
                 !    tmpn(4 ) = rho   (pver)
                 !    tmpn(5 ) = cldmax(pver)
                 !    tmpn(6 ) = pra   (pver)
                 !    tmpn(7 ) = pre   (pver)
                 !    tmpn(8 ) = ums   (pver)
                 !    tmpn(9 ) = mnuccr(pver)
                 !    tmpn(10) = nric  (pver)
                 !    tmpn(11) = nsubr (pver)
                 !    tmpn(12) = npracs(pver)
                 !    tmpn(13) = nnuccr(pver)
                 !    tmpn(14) = nragg (pver)
                 !    tmpn(15) = qric  (pver)
                 !    tmpn(16) = pracs (pver)
                 !    tmpn(17) = uns   (pver)
                 !    tmpn(18) = prai  (pver)
                 !    tmpn(19) = psacws(pver)
                 !    tmpn(20) = bergs (pver)
                 !    tmpn(21) = prds  (pver)
                 !    tmpn(22) = nsic  (pver)
                 !    tmpn(23) = nsubs (pver)
                 !    tmpn(24) = nsagg (pver)
                 !    tmpn(25) = qrtot
                 !    tmpn(26) = nrtot
                 !    tmpn(27) = qstot
                 !    tmpn(28) = nstot
                 !    tmpn(29) = prect
                 !    tmpn(30) = preci
                 !    call put_k_1(tmpn, 30, myid)
                 !endif

                 !if ( k .eq. pver ) then
                 !    call put_k_1(tmpn, 24, myid)
                 !end if 
                 !c........................................................................
                 ! sum over sub-step for average process rates
        
                 ! convert rain/snow q and N for output to history, note, 
                 ! output is for gridbox average
        
                 qrout(k)=qrout(k)+qric(k)*cldmax(k)
                 qsout(k)=qsout(k)+qniic(k)*cldmax(k)
                 nrout(k)=nrout(k)+nric(k)*rho(k)*cldmax(k)
                 nsout(k)=nsout(k)+nsic(k)*rho(k)*cldmax(k)
        
                 tlat1(k)=tlat1(k)+tlat(k)
                 qvlat1(k)=qvlat1(k)+qvlat(k)
                 qctend1(k)=qctend1(k)+qctend(k)
                 qitend1(k)=qitend1(k)+qitend(k)
                 nctend1(k)=nctend1(k)+nctend(k)
                 nitend1(k)=nitend1(k)+nitend(k)
        
                 t(k)=t(k)+tlat(k)*deltat/cpp
                 q(k)=q(k)+qvlat(k)*deltat
                 qc(k)=qc(k)+qctend(k)*deltat
                 qi(k)=qi(k)+qitend(k)*deltat
                 nc(k)=nc(k)+nctend(k)*deltat
                 ni(k)=ni(k)+nitend(k)*deltat
        
                 rainrt1(k)=rainrt1(k)+rainrt(k)
        
                 !divide rain radius over substeps for average
                 if (arcld(k) .gt. 0.d0) then
                    rercld(k)=rercld(k)/arcld(k)
                 end if

        !if ( it .eq. 1  ) then
        !    qrout2(k) = qrout2(k)+t(k)
        !    qsout2(k) = qsout2(k)+q(k)
        !    nrout2(k) = nrout2(k)+qc(k)
        !    nsout2(k) = nsout2(k)+qi(k)
        !    drout2(k) = drout2(k)+nc(k)
        !    dsout2(k) = dsout2(k)+ni(k)
        !endif
        
        
                 !calculate precip fluxes and adding them to summing sub-stepping variables
                 !! flux is zero at top interface
                 rflx(1)=0.0d0
                 sflx(1)=0.0d0
        
                 !! calculating the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
                 rflx(k+1)=qrout(k)*rho(k)*umr(k)
                 sflx(k+1)=qsout(k)*rho(k)*ums(k)
        
                 !! add to summing sub-stepping variable
                 rflx1(k+1)=rflx1(k+1)+rflx(k+1)
                 sflx1(k+1)=sflx1(k+1)+sflx(k+1)
        
                 !c........................................................................
        
              
              end do ! k loop

              if ( myid < 32 ) then
                  tmpn(1 ) = umr   (pver)
                  tmpn(2 ) = unr   (pver)
                  tmpn(3 ) = qniic (pver)
                  tmpn(4 ) = rho   (pver)
                  tmpn(5 ) = cldmax(pver)
                  tmpn(6 ) = pra   (pver)
                  tmpn(7 ) = pre   (pver)
                  tmpn(8 ) = ums   (pver)
                  tmpn(9 ) = mnuccr(pver)
                  tmpn(10) = nric  (pver)
                  tmpn(11) = nsubr (pver)
                  tmpn(12) = npracs(pver)
                  tmpn(13) = nnuccr(pver)
                  tmpn(14) = nragg (pver)
                  tmpn(15) = qric  (pver)
                  tmpn(16) = pracs (pver)
                  tmpn(17) = uns   (pver)
                  tmpn(18) = prai  (pver)
                  tmpn(19) = psacws(pver)
                  tmpn(20) = bergs (pver)
                  tmpn(21) = prds  (pver)
                  tmpn(22) = nsic  (pver)
                  tmpn(23) = nsubs (pver)
                  tmpn(24) = nsagg (pver)
                  tmpn(25) = qrtot
                  tmpn(26) = nrtot
                  tmpn(27) = qstot
                  tmpn(28) = nstot
                  tmpn(29) = prect
                  tmpn(30) = preci
                  call put_k_8(tmpn)
              endif
        
              prect1=prect1+prect
              preci1=preci1+preci
        
           end do ! it loop, sub-step
        
           do k = top_lev, pver
              rate1ord_cw2pr_st(k) = qcsinksum_rate1ord(k)/max(qcsum_rate1ord(k),1.0d-30) 
           end do

           call para_syn_dimk(prect1, myid)
           call para_syn_dimk(preci1, myid)
        
        300 continue  ! continue if no cloud water
        
        ! convert dt from sub-step back to full time step
        deltat=deltat*real(iter)
        
        !do k=1,pver
        !      !qrout2(k) = tlat1(k)
        !      !qsout2(k) = qvlat1(k)
        !      !nrout2(k) = qctend1(k)
        !      !nsout2(k) = qitend1(k)
        !      !drout2(k) = nctend1(k)
        !      !dsout2(k) = nitend1(k)
        !        qrout2(k) = prect1
        !        qsout2(k) = preci1
        !        nrout2(k) = sflx1(k+1)
        !        nsout2(k) = rflx1(k+1)
        !        drout2(k) = rainrt1(k)
        !        dsout2(k) = nsout(k)
        !end do
        
        
        !c.............................................................................
        
        
           ! skip all calculations if no cloud water
           if (ltrue.eq.0) then
        
              do k=1,top_lev-1
                 ! assign zero values for effective radius above 1 mbar
                 effc(k)=0.d0
                 effi(k)=0.d0
                 effc_fn(k)=0.d0
                 lamcrad(k)=0.d0
                 pgamrad(k)=0.d0
                 deffi(k)=0.d0
              end do
        
              do k=top_lev,pver
                 ! assign default values for effective radius
                 effc(k)=10.d0
                 effi(k)=25.d0
                 effc_fn(k)=10.d0
                 lamcrad(k)=0.d0
                 pgamrad(k)=0.d0
                 deffi(k)=0.d0
              end do
              goto 500
           end if
        
           ! initialize nstep for sedimentation sub-steps
           nstep = 1
        
           ! divide precip rate by number of sub-steps to get average over time step
        
           prect=prect1/real(iter)
           preci=preci1/real(iter)
        
           do k=top_lev,pver
        
              ! assign variables back to start-of-timestep values before updating after sub-steps 
        
              t(k)=t1(k)
              q(k)=q1(k)
              qc(k)=qc1(k)
              qi(k)=qi1(k)
              nc(k)=nc1(k)
              ni(k)=ni1(k)
        
              ! divide microphysical tendencies by number of sub-steps to get average over time step
        
              tlat(k)=tlat1(k)/real(iter)
              qvlat(k)=qvlat1(k)/real(iter)
              qctend(k)=qctend1(k)/real(iter)
              qitend(k)=qitend1(k)/real(iter)
              nctend(k)=nctend1(k)/real(iter)
              nitend(k)=nitend1(k)/real(iter)
        
              rainrt(k)=rainrt1(k)/real(iter)
        
              ! divide by number of sub-steps to find final values
              rflx(k+1)=rflx1(k+1)/real(iter)
              sflx(k+1)=sflx1(k+1)/real(iter)
        
              ! divide output precip q and N by number of sub-steps to get average over time step
        
              qrout(k)=qrout(k)/real(iter)
              qsout(k)=qsout(k)/real(iter)
              nrout(k)=nrout(k)/real(iter)
              nsout(k)=nsout(k)/real(iter)
        
              ! divide trop_mozart variables by number of sub-steps to get average over time step 
        
              nevapr(k) = nevapr(k)/real(iter)
              nevapr2(k) = nevapr2(k)/real(iter)
              evapsnow(k) = evapsnow(k)/real(iter)
              prain(k) = prain(k)/real(iter)
              prodsnow(k) = prodsnow(k)/real(iter)
              cmeout(k) = cmeout(k)/real(iter)
        
              cmeiout(k) = cmeiout(k)/real(iter)
              meltsdt(k) = meltsdt(k)/real(iter)
              frzrdt (k) = frzrdt (k)/real(iter)
        
        
              ! microphysics output
              prao(k)=prao(k)/real(iter)
              prco(k)=prco(k)/real(iter)
              mnuccco(k)=mnuccco(k)/real(iter)
              mnuccto(k)=mnuccto(k)/real(iter)
              msacwio(k)=msacwio(k)/real(iter)
              psacwso(k)=psacwso(k)/real(iter)
              bergso(k)=bergso(k)/real(iter)
              bergo(k)=bergo(k)/real(iter)
              prcio(k)=prcio(k)/real(iter)
              praio(k)=praio(k)/real(iter)
        
              mnuccro(k)=mnuccro(k)/real(iter)
              pracso (k)=pracso (k)/real(iter)
        
              mnuccdo(k)=mnuccdo(k)/real(iter)
        
              ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
              nevapr(k) = nevapr(k) + evapsnow(k)
              prer_evap(k) = nevapr2(k)
              prain(k) = prain(k) + prodsnow(k)
        
              !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              ! calculate sedimentation for cloud water and ice
              !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
              ! update in-cloud cloud mixing ratio and number concentration 
              ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
              ! note: these are in-cloud values***, hence we divide by cloud fraction
        
              dumc(k) = (qc(k)+qctend(k)*deltat)/lcldm(k)
              dumi(k) = (qi(k)+qitend(k)*deltat)/icldm(k)
              dumnc(k) = max((nc(k)+nctend(k)*deltat)/lcldm(k),0.d0)
              dumni(k) = max((ni(k)+nitend(k)*deltat)/icldm(k),0.d0)
        
              ! obtain new slope parameter to avoid possible singularity
        
              if (dumi(k).ge.qsmall) then
                 ! add upper limit to in-cloud number concentration to prevent numerical error
                 dumni(k)=min(dumni(k),dumi(k)*1.d20)
        
                 tmp1 = (cons1*ci* &
                      dumni(k)/dumi(k))
                 tmp2 = 1.d0/di
                 call math_agent_2i1o(power, tmp1, tmp2, lami(k))

                 !lami(k) = (cons1*ci* &
                 !     dumni(k)/dumi(k))**(1.d0/di)
                 lami(k)=max(lami(k),lammini)
                 lami(k)=min(lami(k),lammaxi)
              else
                 lami(k)=0.d0
              end if
        
              if (dumc(k).ge.qsmall) then
                 ! add upper limit to in-cloud number concentration to prevent numerical error
                 dumnc(k)=min(dumnc(k),dumc(k)*1.d20)
                 ! add lower limit to in-cloud number concentration
                 dumnc(k)=max(dumnc(k),cdnl/rho(k)) ! sghan minimum in #/cm3 
                 pgam(k)=dp0_0005714*(ncic(k)/1.d6*rho(k))+dp0_2714
                 pgam(k)=1.d0/(pgam(k)*pgam(k))-1.d0
                 !pgam(k)=1.d0/(pgam(k)**2)-1.d0
                 pgam(k)=max(pgam(k),2.d0)
                 pgam(k)=min(pgam(k),15.d0)
        
                 tmp1 = (pi/6.d0*rhow*dumnc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/ &
                      (dumc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)))
                 tmp2 = 1.d0/3.d0
                 call math_agent_2i1o(power, tmp1, tmp2, lamc(k))
                 !lamc(k) = (pi/6.d0*rhow*dumnc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0)/ &
                 !     (dumc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0)))**(1.d0/3.d0)
                 lammin = (pgam(k)+1.d0)/50.d-6
                 lammax = (pgam(k)+1.d0)/2.d-6
                 lamc(k)=max(lamc(k),lammin)
                 lamc(k)=min(lamc(k),lammax)
              else
                 lamc(k)=0.d0
              end if
        
              ! calculate number and mass weighted fall velocity for droplets
              ! include effects of sub-grid distribution of cloud water
        
        
              if (dumc(k).ge.qsmall) then
                 call math_agent_2i1o(power, lamc(k), bc, tmp1)
                 unc = acn(k)*shr_spfn_gamma_nonintrinsic_r8(1.d0+bc+pgam(k),math_agent,Pgama,Qgama,Cgama,Sgama)/(tmp1*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama))
                 umc = acn(k)*shr_spfn_gamma_nonintrinsic_r8(4.d0+bc+pgam(k),math_agent,Pgama,Qgama,Cgama,Sgama)/(tmp1*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama))
                 !unc = acn(k)*shr_spfn_gamma_nonintrinsic_r8(1.d0+bc+pgam(k))/(lamc(k)**bc*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0))
                 !umc = acn(k)*shr_spfn_gamma_nonintrinsic_r8(4.d0+bc+pgam(k))/(lamc(k)**bc*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0))
                 ! fallspeed for output
                 vtrmc(k)=umc
              else
                 umc = 0.d0
                 unc = 0.d0
              end if
        
              ! calculate number and mass weighted fall velocity for cloud ice
        
              if (dumi(k).ge.qsmall) then
                 call math_agent_2i1o(power, lami(k), bi, tmp1)
                 uni =  ain(k)*cons16/tmp1
                 umi = ain(k)*cons17/(6.d0*tmp1)
                 !uni =  ain(k)*cons16/lami(k)**bi
                 !umi = ain(k)*cons17/(6.d0*lami(k)**bi)
                 uni=min(uni,1.2d0*rhof(k))
                 umi=min(umi,1.2d0*rhof(k))
        
                 ! fallspeed
                 vtrmi(k)=umi
              else
                 umi = 0.d0
                 uni = 0.d0
              end if
        
              fi(k) = g*rho(k)*umi
              fni(k) = g*rho(k)*uni
              fc(k) = g*rho(k)*umc
              fnc(k) = g*rho(k)*unc
        
              ! calculate number of split time steps to ensure courant stability criteria
              ! for sedimentation calculations
        
              rgvm = max(fi(k),fc(k),fni(k),fnc(k))
              nstep = max(int(rgvm*deltat/pdel(k)+1.d0),nstep)
        
              ! redefine dummy variables - sedimentation is calculated over grid-scale
              ! quantities to ensure conservation
        
              dumc(k) = (qc(k)+qctend(k)*deltat)
              dumi(k) = (qi(k)+qitend(k)*deltat)
              dumnc(k) = max((nc(k)+nctend(k)*deltat),0.d0)
              dumni(k) = max((ni(k)+nitend(k)*deltat),0.d0)
        
              if (dumc(k).lt.qsmall) dumnc(k)=0.d0
              if (dumi(k).lt.qsmall) dumni(k)=0.d0
        
           end do       !!! vertical loop
           
           call para_syn_max_dimk(nstep, myid) 

           !do k=1,pver
           !     qrout2(k) = nstep 
           !     qsout2(k) = nstep 
           !     nrout2(k) = nstep 
           !     nsout2(k) = nstep 
           !     drout2(k) = nstep 
           !     dsout2(k) = nstep 
           !end do

           do n = 1,nstep  !! loop over sub-time step to ensure stability
        
              do k = top_lev,pver
                 if (do_cldice) then
                    falouti(k) = fi(k)*dumi(k)
                    faloutni(k) = fni(k)*dumni(k)
                 else
                    falouti(k)  = 0.d0
                    faloutni(k) = 0.d0
                 end if
        
                 faloutc(k) = fc(k)*dumc(k)
                 faloutnc(k) = fnc(k)*dumnc(k)
              end do
        
          if ( myid < 32 ) then
              ! top of model
        
              k = top_lev
              faltndi = falouti(k)/pdel(k)
              faltndni = faloutni(k)/pdel(k)
              faltndc = faloutc(k)/pdel(k)
              faltndnc = faloutnc(k)/pdel(k)
        
              ! add fallout terms to microphysical tendencies
        
              qitend(k) = qitend(k)-faltndi/nstep
              nitend(k) = nitend(k)-faltndni/nstep
              qctend(k) = qctend(k)-faltndc/nstep
              nctend(k) = nctend(k)-faltndnc/nstep
        
              ! sedimentation tendencies for output
              qcsedten(k)=qcsedten(k)-faltndc/nstep
              qisedten(k)=qisedten(k)-faltndi/nstep
        
              dumi(k) = dumi(k)-faltndi*deltat/nstep
              dumni(k) = dumni(k)-faltndni*deltat/nstep
              dumc(k) = dumc(k)-faltndc*deltat/nstep
              dumnc(k) = dumnc(k)-faltndnc*deltat/nstep

              kst = top_lev+1
          else
              kst = 1
          endif
        
        !if (myid >= 32 ) then 
        !    call get_k_1(tmpn, 1, myid)
        !endif
        !if ( myid < 32 ) then 
        !    call put_k_1(tmpn, 1, myid)
        !endif
              !do k = top_lev+1,pver
              do k = kst,pver
        
                 ! for cloud liquid and ice, if cloud fraction increases with height
                 ! then add flux from above to both vapor and cloud water of current level
                 ! this means that flux entering clear portion of cell from above evaporates
                 ! instantly

                 if ( myid >= 32 .and. k .eq. 1) then
                     call get_k_3(tmpn)
                 else
                     tmpn(1) = lcldm(k-1)
                     tmpn(2) = icldm(k-1)
                     tmpn(3) = falouti(k-1)
                     tmpn(4) = faloutni(k-1)
                     tmpn(5) = faloutc(k-1)
                     tmpn(6) = faloutnc(k-1)
                     tmpn(7) = fni(k-1)
                     tmpn(8) = pdel(k-1)
                     tmpn(9) = fi(k-1)
                     tmpn(10)= fnc(k-1)
                     tmpn(11)= fc(k-1)
                 endif
                     
                 dum=lcldm(k)/tmpn(1)
                 dum=min(dum,1.d0)
                 dum1=icldm(k)/tmpn(2)
                 dum1=min(dum1,1.d0)
        
                 faltndqie=(falouti(k)-tmpn(3))/pdel(k)
                 faltndi=(falouti(k)-dum1*tmpn(3))/pdel(k)
                 faltndni=(faloutni(k)-dum1*tmpn(4))/pdel(k)
                 faltndqce=(faloutc(k)-tmpn(5))/pdel(k)
                 faltndc=(faloutc(k)-dum*tmpn(5))/pdel(k)
                 faltndnc=(faloutnc(k)-dum*tmpn(6))/pdel(k)
        
                 !dum=lcldm(k)/lcldm(k-1)
                 !dum=min(dum,1.d0)
                 !dum1=icldm(k)/icldm(k-1)
                 !dum1=min(dum1,1.d0)
        
                 !faltndqie=(falouti(k)-falouti(k-1))/pdel(k)
                 !faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(k)
                 !faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(k)
                 !faltndqce=(faloutc(k)-faloutc(k-1))/pdel(k)
                 !faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(k)
                 !faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(k)
        
                 ! add fallout terms to eulerian tendencies
        
                 qitend(k) = qitend(k)-faltndi/nstep
                 nitend(k) = nitend(k)-faltndni/nstep
                 qctend(k) = qctend(k)-faltndc/nstep
                 nctend(k) = nctend(k)-faltndnc/nstep
        
                 !if ( nstep .eq. 1 ) then
                 !    qrout2(k) = qrout2(k) + qcsedten(k)
                 !    qsout2(k) = qsout2(k) + faltndc/nstep
                 !    nrout2(k) = nrout2(k) + faltndc
                 !    nsout2(k) = nsout2(k) + nstep
                 !endif
                 ! sedimentation tendencies for output
                 qcsedten(k)=qcsedten(k)-faltndc/nstep
                 qisedten(k)=qisedten(k)-faltndi/nstep
        
                 !if ( nstep .eq. 1 ) then
                 !    drout2(k) = drout2(k) + qcsedten(k)          
                 !    dsout2(k) = dsout2(k) + qisedten(k)
                 !endif
                 ! add terms to to evap/sub of cloud water
        
                 qvlat(k)=qvlat(k)-(faltndqie-faltndi)/nstep
                 ! for output
                 qisevap(k)=qisevap(k)-(faltndqie-faltndi)/nstep
                 qvlat(k)=qvlat(k)-(faltndqce-faltndc)/nstep

                 !if ( n .eq. 2 ) then
                 !    qrout2(k) = qrout2(k) + qcsevap(k) 
                 !    qsout2(k) = qsout2(k) + faltndqce
                 !    nrout2(k) = nrout2(k) + faltndc
                 !    nsout2(k) = nsout2(k) + nstep
                 !    drout2(k) = drout2(k) + (faltndqce-faltndc)/nstep          
                 !    dsout2(k) = dsout2(k) + qcsevap(k)-(faltndqce-faltndc)/nstep
                 !endif

                 ! for output
                 qcsevap(k)=qcsevap(k)-(faltndqce-faltndc)/nstep
        
                 tlat(k)=tlat(k)+(faltndqie-faltndi)*xxls/nstep
                 tlat(k)=tlat(k)+(faltndqce-faltndc)*xxlv/nstep
        
                 dumi(k) = dumi(k)-faltndi*deltat/nstep
                 dumni(k) = dumni(k)-faltndni*deltat/nstep
                 dumc(k) = dumc(k)-faltndc*deltat/nstep
                 dumnc(k) = dumnc(k)-faltndnc*deltat/nstep
        
                 Fni(K)=MAX(Fni(K)/pdel(k),tmpn(7)/tmpn(8))*pdel(k)
                 FI(K)=MAX(FI(K)/pdel(k),tmpn(9)/tmpn(8))*pdel(k)
                 fnc(k)=max(fnc(k)/pdel(k),tmpn(10)/tmpn(8))*pdel(k)
                 Fc(K)=MAX(Fc(K)/pdel(k),tmpn(11)/tmpn(8))*pdel(k)

                 !if ( myid < 32 .and. k .eq. pver ) then
                 !    tmpn(1) = lcldm(pver)
                 !    tmpn(2) = icldm(pver)
                 !    tmpn(3) = falouti(pver)
                 !    tmpn(4) = faloutni(pver)
                 !    tmpn(5) = faloutc(pver)
                 !    tmpn(6) = faloutnc(pver)
                 !    tmpn(7) = fni(pver)
                 !    tmpn(8) = pdel(pver)
                 !    tmpn(9) = fi(pver)
                 !    tmpn(10)= fnc(pver)
                 !    tmpn(11)= fc(pver)
                 !    call put_k_1(tmpn, 11, myid)
                 !endif

                 !Fni(K)=MAX(Fni(K)/pdel(k),Fni(K-1)/pdel(k-1))*pdel(k)
                 !FI(K)=MAX(FI(K)/pdel(k),FI(K-1)/pdel(k-1))*pdel(k)
                 !fnc(k)=max(fnc(k)/pdel(k),fnc(k-1)/pdel(k-1))*pdel(k)
                 !Fc(K)=MAX(Fc(K)/pdel(k),Fc(K-1)/pdel(k-1))*pdel(k)
        
              end do   !! k loop

              if ( myid < 32 ) then
                  tmpn(1) = lcldm(pver)
                  tmpn(2) = icldm(pver)
                  tmpn(3) = falouti(pver)
                  tmpn(4) = faloutni(pver)
                  tmpn(5) = faloutc(pver)
                  tmpn(6) = faloutnc(pver)
                  tmpn(7) = fni(pver)
                  tmpn(8) = pdel(pver)
                  tmpn(9) = fi(pver)
                  tmpn(10)= fnc(pver)
                  tmpn(11)= fc(pver)
                  call put_k_3(tmpn)
              endif


        
              ! units below are m/s
              ! cloud water/ice sedimentation flux at surface 
              ! is added to precip flux at surface to get total precip (cloud + precip water)
              ! rate
        
              prect = prect+(faloutc(pver)+falouti(pver))/g/nstep/1000.d0  
              preci = preci+(falouti(pver))/g/nstep/1000.d0
        
           end do   !! nstep loop
           !do k=1,pver
           !     qrout2(k) = qrout2(k) + qitend(k) 
           !     qsout2(k) = qsout2(k) + nitend(k) 
           !     nrout2(k) = nrout2(k) + qctend(k)
           !     nsout2(k) = nsout2(k) + nctend(k) 
           !     drout2(k) = drout2(k) + qcsedten(k)          
           !     dsout2(k) = dsout2(k) + qisedten(k)
           !end do        
                                
           ! end sedimentation
           !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
           ! get new update for variables that includes sedimentation tendency
           ! note : here dum variables are grid-average, NOT in-cloud
        
           do k=top_lev,pver
        
              dumc(k) = max(qc(k)+qctend(k)*deltat,0.d0)
              dumi(k) = max(qi(k)+qitend(k)*deltat,0.d0)
              dumnc(k) = max(nc(k)+nctend(k)*deltat,0.d0)
              dumni(k) = max(ni(k)+nitend(k)*deltat,0.d0)
        
              if (dumc(k).lt.qsmall) dumnc(k)=0.d0
              if (dumi(k).lt.qsmall) dumni(k)=0.d0
        
              ! calculate instantaneous processes (melting, homogeneous freezing)
              if (do_cldice) then
        
                 if (t(k)+tlat(k)/cpp*deltat > tmelt) then
                    if (dumi(k) > 0.d0) then
        
                       ! limit so that melting does not push temperature below freezing
                       dum = -dumi(k)*xlf/cpp
                       if (t(k)+tlat(k)/cpp*deltat+dum.lt.tmelt) then
                          dum = (t(k)+tlat(k)/cpp*deltat-tmelt)*cpp/xlf
                          dum = dum/dumi(k)*xlf/cpp 
                          dum = max(0.d0,dum)
                          dum = min(1.d0,dum)
                       else
                          dum = 1.d0
                       end if
        
                       qctend(k)=qctend(k)+dum*dumi(k)/deltat
        
                       ! for output
                       melto(k)=dum*dumi(k)/deltat
        
                       ! assume melting ice produces droplet
                       ! mean volume radius of 8 micron
        
                       nctend(k)=nctend(k)+3.d0*dum*dumi(k)/deltat/ &
                            (4.d0*pi*5.12d-16*rhow)
        
                       qitend(k)=((1.d0-dum)*dumi(k)-qi(k))/deltat
                       nitend(k)=((1.d0-dum)*dumni(k)-ni(k))/deltat
                       tlat(k)=tlat(k)-xlf*dum*dumi(k)/deltat
                    end if
                 end if
        
                 ! homogeneously freeze droplets at -40 C
        
                 if (t(k)+tlat(k)/cpp*deltat < dp233_15) then
                    if (dumc(k) > 0.d0) then
        
                       ! limit so that freezing does not push temperature above threshold
                       dum = dumc(k)*xlf/cpp
                       if (t(k)+tlat(k)/cpp*deltat+dum.gt.dp233_15) then
                          dum = -(t(k)+tlat(k)/cpp*deltat-dp233_15)*cpp/xlf
                          dum = dum/dumc(k)*xlf/cpp
                          dum = max(0.d0,dum)
                          dum = min(1.d0,dum)
                       else
                          dum = 1.d0
                       end if
        
                       qitend(k)=qitend(k)+dum*dumc(k)/deltat
                       ! for output
                       homoo(k)=dum*dumc(k)/deltat
        
                       ! assume 25 micron mean volume radius of homogeneously frozen droplets
                       ! consistent with size of detrained ice in stratiform.F90
                       nitend(k)=nitend(k)+dum*3.d0*dumc(k)/(4.d0*dp3_14*dp1_563d_14* &
                            500.d0)/deltat
                       qctend(k)=((1.d0-dum)*dumc(k)-qc(k))/deltat
                       nctend(k)=((1.d0-dum)*dumnc(k)-nc(k))/deltat
                       tlat(k)=tlat(k)+xlf*dum*dumc(k)/deltat
                    end if
                 end if
        
                 ! remove any excess over-saturation, which is possible due to non-linearity when adding 
                 ! together all microphysical processes
                 ! follow code similar to old CAM scheme
        
                 qtmp=q(k)+qvlat(k)*deltat
                 ttmp=t(k)+tlat(k)/cpp*deltat
        
                 esn = svp_water(ttmp, lg10, power, svp_tboil)  ! use rhw to allow ice supersaturation
                 qsn = svp_to_qsat(esn, p(k), svp_epsilo, svp_omeps)

          !tmp1 = 10.d0
          !tmp2 = svp_tboil/ttmp
          !call math_agent_1i1o(lg10, tmp2, tmp3)
          !tmp4 = 11.344d0*(1.d0-ttmp/svp_tboil)
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !tmp4 = (-3.49149d0*(svp_tboil/ttmp-1.d0))
          !call math_agent_2i1o(power, tmp1, tmp4, tmp5)
          !tmp4 = 1013.246d0
          !call math_agent_1i1o(lg10, tmp4, tmp1)
          !tmp4 = -7.90298d0*(svp_tboil/ttmp-1.d0)+ &
          !     5.02808d0*tmp3- &
          !     1.3816d-7*(tmp2-1.d0)+ &
          !     8.1328d-3*(tmp5-1.d0)+ &
          !     tmp1
          !tmp1 = 10.d0
          !call math_agent_2i1o(power, tmp1, tmp4, tmp2)
          !esn = tmp2 * 100.d0

          !if ( (p(k) - esn) <= 0.d0 ) then
          !   qsn = 1.0d0
          !else
          !   qsn = svp_epsilo*esn / (p(k) - svp_omeps*esn)
          !end if
        
                 if (qtmp > qsn .and. qsn > 0) then
                    ! expression below is approximate since there may be ice deposition
                    dum = (qtmp-qsn)/(1.d0+cons27*qsn/(cpp*rv*ttmp*ttmp))/deltat
                    ! add to output cme
                    cmeout(k) = cmeout(k)+dum
                    ! now add to tendencies, partition between liquid and ice based on temperature
                    if (ttmp > 268.15d0) then
                       dum1=0.0d0
                       ! now add to tendencies, partition between liquid and ice based on te
                    else if (ttmp < dp238_15) then
                       dum1=1.0d0
                    else
                       dum1=(268.15d0-ttmp)/30.d0
                    end if
        
                    dum = (qtmp-qsn)/(1.d0+(xxls*dum1+xxlv*(1.d0-dum1))*(xxls*dum1+xxlv*(1.d0-dum1)) &
                         *qsn/(cpp*rv*ttmp*ttmp))/deltat
                    qctend(k)=qctend(k)+dum*(1.d0-dum1)
                    ! for output
                    qcreso(k)=dum*(1.d0-dum1)
                    qitend(k)=qitend(k)+dum*dum1
                    qireso(k)=dum*dum1
                    qvlat(k)=qvlat(k)-dum
                    ! for output
                    qvres(k)=-dum
                    tlat(k)=tlat(k)+dum*(1.d0-dum1)*xxlv+dum*dum1*xxls
                 end if
              end if
        
              !...............................................................................
              ! calculate effective radius for pass to radiation code
              ! if no cloud water, default value is 10 micron for droplets,
              ! 25 micron for cloud ice
        
              ! update cloud variables after instantaneous processes to get effective radius
              ! variables are in-cloud to calculate size dist parameters
        
              dumc(k) = max(qc(k)+qctend(k)*deltat,0.d0)/lcldm(k)
              dumi(k) = max(qi(k)+qitend(k)*deltat,0.d0)/icldm(k)
              dumnc(k) = max(nc(k)+nctend(k)*deltat,0.d0)/lcldm(k)
              dumni(k) = max(ni(k)+nitend(k)*deltat,0.d0)/icldm(k)
        
              ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
        
              dumc(k)=min(dumc(k),5.d-3)
              dumi(k)=min(dumi(k),5.d-3)
        
              !...................
              ! cloud ice effective radius
        
              if (dumi(k).ge.qsmall) then
                 ! add upper limit to in-cloud number concentration to prevent numerical error
                 dumni(k)=min(dumni(k),dumi(k)*1.d20)

                 tmp1 = (cons1*ci*dumni(k)/dumi(k))
                 tmp2 = 1.d0/di
                 call math_agent_2i1o(power, tmp1, tmp2, lami(k))
                 !lami(k) = (cons1*ci*dumni(k)/dumi(k))**(1.d0/di)
        
                 if (lami(k).lt.lammini) then
                    lami(k) = lammini

                    tmp1 = di+1.d0
                    call math_agent_2i1o(power, lami(k), tmp1, tmp2)
                    n0i(k) = tmp2*dumi(k)/(ci*cons1)
                    !n0i(k) = lami(k)**(di+1.d0)*dumi(k)/(ci*cons1)
                    niic(k) = n0i(k)/lami(k)
                    ! adjust number conc if needed to keep mean size in reasonable range
                    if (do_cldice) nitend(k)=(niic(k)*icldm(k)-ni(k))/deltat
        
                 else if (lami(k).gt.lammaxi) then
                    lami(k) = lammaxi
                    tmp1 = di+1.d0
                    call math_agent_2i1o(power, lami(k), tmp1, tmp2)
                    n0i(k) = tmp2*dumi(k)/(ci*cons1)
                    !n0i(k) = lami(k)**(di+1.d0)*dumi(k)/(ci*cons1)
                    niic(k) = n0i(k)/lami(k)
                    ! adjust number conc if needed to keep mean size in reasonable range
                    if (do_cldice) nitend(k)=(niic(k)*icldm(k)-ni(k))/deltat
                 end if
                 effi(k) = 1.5d0/lami(k)*1.d6
        
              else
                 effi(k) = 25.d0
              end if
        
              ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
              ! radius has already been determined from the size distribution.
              if (.not. do_cldice) then
                 effi(k) = re_ice(k) * 1d6      ! m -> um
              end if
        
              !...................
              ! cloud droplet effective radius
        
              if (dumc(k).ge.qsmall) then
        
                 ! add upper limit to in-cloud number concentration to prevent numerical error
                 dumnc(k)=min(dumnc(k),dumc(k)*1.d20)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! set tendency to ensure minimum droplet concentration
                 ! after update by microphysics, except when lambda exceeds bounds on mean drop
                 ! size or if there is no cloud water
                 if (dumnc(k).lt.cdnl/rho(k)) then   
                    nctend(k)=(cdnl/rho(k)*lcldm(k)-nc(k))/deltat   
                 end if
                 dumnc(k)=max(dumnc(k),cdnl/rho(k)) ! sghan minimum in #/cm3 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 pgam(k)=dp0_0005714*(ncic(k)/1.d6*rho(k))+dp0_2714
                 pgam(k)=1.d0/(pgam(k)*pgam(k))-1.d0
                 pgam(k)=max(pgam(k),2.d0)
                 pgam(k)=min(pgam(k),15.d0)
        
                 tmp1 = (pi/6.d0*rhow*dumnc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/ &
                      (dumc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)))
                 tmp2 = 1.d0/3.d0
                 call math_agent_2i1o(power, tmp1, tmp2, lamc(k))
                 !lamc(k) = (pi/6.d0*rhow*dumnc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0)/ &
                 !     (dumc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0)))**(1.d0/3.d0)
                 lammin = (pgam(k)+1.d0)/50.d-6
                 ! Multiply by omsm to fit within RRTMG's table.
                 lammax = (pgam(k)+1.d0)*omsm/2.d-6
                 if (lamc(k).lt.lammin) then
                    lamc(k) = lammin
                    ncic(k) = 6.d0*lamc(k)*lamc(k)*lamc(k)*dumc(k)* &
                         shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/ &
                         (pi*rhow*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama))
                    ! adjust number conc if needed to keep mean size in reasonable range
                    nctend(k)=(ncic(k)*lcldm(k)-nc(k))/deltat
        
                 else if (lamc(k).gt.lammax) then
                    lamc(k) = lammax
                    ncic(k) = 6.d0*lamc(k)*lamc(k)*lamc(k)*dumc(k)* &
                         shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/ &
                         (pi*rhow*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama))
                    ! adjust number conc if needed to keep mean size in reasonable range
                    nctend(k)=(ncic(k)*lcldm(k)-nc(k))/deltat
                 end if
        
                 effc(k) = &
                      shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/ &
                      shr_spfn_gamma_nonintrinsic_r8(pgam(k)+3.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/lamc(k)/2.d0*1.d6
                 !assign output fields for shape here
                 lamcrad(k)=lamc(k)
                 pgamrad(k)=pgam(k)
        
              else
                 effc(k) = 10.d0
                 lamcrad(k)=0.d0
                 pgamrad(k)=0.d0
              end if
        
              ! ice effective diameter for david mitchell's optics
              if (do_cldice) then
                 deffi(k)=effi(k)*rhoi/917.d0*2.d0
              else
                 deffi(k)=effi(k) * 2.d0
              end if
        
        
        !!! recalculate effective radius for constant number, in order to separate
              ! first and second indirect effects
              ! assume constant number of 10^8 kg-1
        
              dumnc(k)=1.d8
        
              if (dumc(k).ge.qsmall) then
                 pgam(k)=dp0_0005714*(ncic(k)/1.d6*rho(k))+dp0_2714
                 pgam(k)=1.d0/(pgam(k)*pgam(k))-1.d0
                 pgam(k)=max(pgam(k),2.d0)
                 pgam(k)=min(pgam(k),15.d0)
        
                 tmp1 = (pi/6.d0*rhow*dumnc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/ &
                      (dumc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0,math_agent,Pgama,Qgama,Cgama,Sgama)))
                 tmp2 = 1.d0/3.d0
                 call math_agent_2i1o(power, tmp1, tmp2, lamc(k))
                 !lamc(k) = (pi/6.d0*rhow*dumnc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0)/ &
                 !     (dumc(k)*shr_spfn_gamma_nonintrinsic_r8(pgam(k)+1.d0)))**(1.d0/3.d0)
                 lammin = (pgam(k)+1.d0)/50.d-6
                 lammax = (pgam(k)+1.d0)/2.d-6
                 if (lamc(k).lt.lammin) then
                    lamc(k) = lammin
                 else if (lamc(k).gt.lammax) then
                    lamc(k) = lammax
                 end if
                 effc_fn(k) = &
                      shr_spfn_gamma_nonintrinsic_r8(pgam(k)+4.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/ &
                      shr_spfn_gamma_nonintrinsic_r8(pgam(k)+3.d0,math_agent,Pgama,Qgama,Cgama,Sgama)/lamc(k)/2.d0*1.d6
        
              else
                 effc_fn(k) = 10.d0
              end if
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!
        
           end do ! vertical k loop
           !do k=1,pver
           !     qrout2(k) = qcreso(k) 
           !     qsout2(k) = qitend(k) 
           !     nrout2(k) = qireso(k)
           !     nsout2(k) = qvlat(k) 
           !     drout2(k) = deffi(k)          
           !     dsout2(k) = effc_fn(k)
           !end do        
        
        500 continue
        
           do k=top_lev,pver
              ! if updated q (after microphysics) is zero, then ensure updated n is also zero
        
              if (qc(k)+qctend(k)*deltat.lt.qsmall) nctend(k)=-nc(k)/deltat
              if (do_cldice .and. qi(k)+qitend(k)*deltat.lt.qsmall) nitend(k)=-ni(k)/deltat
           end do
        
        
        ! add snow ouptut
           do k=top_lev,pver
              if (qsout(k).gt.1.d-7.and.nsout(k).gt.0.d0) then
                 tmp1 = (pi * rhosn * nsout(k)/qsout(k))
                 tmp2 = -1.d0/3.d0
                 call math_agent_2i1o(power, tmp1, tmp2, tmp3)
                 dsout(k)=3.d0*rhosn/917.d0*tmp3
                 !dsout(k)=3.d0*rhosn/917.d0*(pi * rhosn * nsout(k)/qsout(k))**(-1.d0/3.d0)
              endif
           end do
        
        !calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual
           do k=top_lev,pver
              !! RAIN
              if (qrout(k).gt.1.d-7.and.nrout(k).gt.0.d0) then
                 tmp1 = (pi * rhow * nrout(k)/qrout(k))
                 tmp2 = -1.d0/3.d0
                 call math_agent_2i1o(power, tmp1, tmp2, tmp3)
                 reff_rain(k)=1.5d0*tmp3*1.d6
                 !reff_rain(k)=1.5d0*(pi * rhow * nrout(k)/qrout(k))**(-1.d0/3.d0)*1.d6
              endif
              !! SNOW
              if (qsout(k).gt.1.d-7.and.nsout(k).gt.0.d0) then
                 tmp1 = (pi * rhosn * nsout(k)/qsout(k))
                 tmp2 = -1.d0/3.d0
                 call math_agent_2i1o(power, tmp1, tmp2, tmp3)
                 reff_snow(k)=1.5d0*tmp3*1.d6
                 !reff_snow(k)=1.5d0*(pi * rhosn * nsout(k)/qsout(k))**(-1.d0/3.d0)*1.d6
              end if
           end do
        
        ! analytic radar reflectivity
        ! formulas from Matthew Shupe, NOAA/CERES
        ! *****note: radar reflectivity is local (in-precip average)
        ! units of mm^6/m^3
        
           do k=top_lev,pver
              if (qc(k)+qctend(k)*deltat.ge.qsmall) then
                 dum=((qc(k)+qctend(k)*deltat)/lcldm(k)*rho(k)*1000.d0)*((qc(k)+qctend(k)*deltat)/lcldm(k)*rho(k)*1000.d0) &
                      /(0.109d0*(nc(k)+nctend(k)*deltat)/lcldm(k)*rho(k)/1.d6)*lcldm(k)/cldmax(k)
              else
                 dum=0.d0
              end if
              if (qi(k)+qitend(k)*deltat.ge.qsmall) then
                 tmp1 = ((qi(k)+qitend(k)*deltat)*rho(k)/icldm(k)*1000.d0/0.1d0)
                 tmp2 = (1.d0/0.63d0)
                 call math_agent_2i1o(power, tmp1, tmp2, tmp3)
                 dum1=tmp3*icldm(k)/cldmax(k)
                 !dum1=((qi(k)+qitend(k)*deltat)*rho(k)/icldm(k)*1000.d0/0.1d0)**(1.d0/0.63d0)*icldm(k)/cldmax(k)
              else 
                 dum1=0.d0
              end if
        
              if (qsout(k).ge.qsmall) then
                 tmp1 = (qsout(k)*rho(k)*1000.d0/0.1d0)
                 tmp2 = (1.d0/0.63d0)
                 call math_agent_2i1o(power, tmp1, tmp2, tmp3)
                 dum1 = dum1 + tmp3
                 !dum1=dum1+(qsout(k)*rho(k)*1000.d0/0.1d0)**(1.d0/0.63d0)
              end if
        
              refl(k)=dum+dum1
        
              ! add rain rate, but for 37 GHz formulation instead of 94 GHz
              ! formula approximated from data of Matrasov (2007)
              ! rainrt is the rain rate in mm/hr
              ! reflectivity (dum) is in DBz
        
              if (rainrt(k).ge.0.001d0) then
                 call math_agent_2i1o(power, rainrt(k), 6.d0, tmp1)
                 call math_agent_1i1o(lg10, tmp1, tmp2)
                 dum = tmp2 + 16.d0
                 !dum=log10(rainrt(k)**6.d0)+16.d0
        
                 ! convert from DBz to mm^6/m^3
        
                 tmp1 = dum/10.d0
                 call math_agent_2i1o(power, 10.d0, tmp1, dum)
                 !dum = 10.d0**(dum/10.d0)
              else
                 ! don't include rain rate in R calculation for values less than 0.001 mm/hr
                 dum=0.d0
              end if
        
              ! add to refl
        
              refl(k)=refl(k)+dum
        
              !output reflectivity in Z.
              areflz(k)=refl(k)
        
              ! convert back to DBz 
        
              if (refl(k).gt.minrefl) then 
                 refl(k)=10.d0*log10(refl(k))
              else
                 refl(k)=-9999.d0
              end if
        
              !set averaging flag
              if (refl(k).gt.mindbz) then 
                 arefl(k)=refl(k)
                 frefl(k)=1.0d0  
              else
                 arefl(k)=0.d0
                 areflz(k)=0.d0
                 frefl(k)=0.d0
              end if
        
              ! bound cloudsat reflectivity
        
              csrfl(k)=min(csmax,refl(k))
        
              !set averaging flag
              if (csrfl(k).gt.csmin) then 
                 acsrfl(k)=refl(k)
                 fcsrfl(k)=1.0d0  
              else
                 acsrfl(k)=0.d0
                 fcsrfl(k)=0.d0
              end if
        
           end do
        
        
        ! averaging for snow and rain number and diameter
        
        qrout2(:)=0.d0
        qsout2(:)=0.d0
        nrout2(:)=0.d0
        nsout2(:)=0.d0
        drout2(:)=0.d0
        dsout2(:)=0.d0
        freqs(:)=0.d0
        freqr(:)=0.d0
           do k=top_lev,pver
              if (qrout(k).gt.1.d-7.and.nrout(k).gt.0.d0) then
                 qrout2(k)=qrout(k)
                 nrout2(k)=nrout(k)

                 tmp1 = (pi * rhow * nrout(k)/qrout(k))
                 tmp2 = -1.d0/3.d0
                 call math_agent_2i1o(power, tmp1, tmp2, drout2(k))
                 !drout2(k)=(pi * rhow * nrout(k)/qrout(k))**(-1.d0/3.d0)
                 freqr(k)=1.d0
              endif
              if (qsout(k).gt.1.d-7.and.nsout(k).gt.0.d0) then
                 qsout2(k)=qsout(k)
                 nsout2(k)=nsout(k)
                 tmp1 = (pi * rhosn * nsout(k)/qsout(k))
                 tmp2 = -1.d0/3.d0
                 call math_agent_2i1o(power, tmp1, tmp2, dsout2(k))
                 !dsout2(k)=(pi * rhosn * nsout(k)/qsout(k))**(-1.d0/3.d0)
                 freqs(k)=1.d0
              endif
           end do
        
        ! output activated liquid and ice (convert from #/kg -> #/m3)
           do k=top_lev,pver
              ncai(k)=dum2i(k)*rho(k)
              ncal(k)=dum2l(k)*rho(k)
           end do
        
        
        !redefine fice here....
        nfice(:)=0.d0
        do k=top_lev,pver
              dumc(k) = (qc(k)+qctend(k)*deltat)
              dumi(k) = (qi(k)+qitend(k)*deltat)
              dumfice=qsout(k) + qrout(k) + dumc(k) + dumi(k)  
        
              if (dumfice.gt.qsmall.and.(qsout(k)+dumi(k).gt.qsmall)) then
                 nfice(k)=(qsout(k) + dumi(k))/dumfice
              endif
        
              if (nfice(k).gt.1.d0) then
                 nfice(k)=1.d0
              endif
        
        enddo

        end subroutine


