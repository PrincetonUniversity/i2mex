FUNCTION R8_ZEROIN(AX,BX,F,TOL,EPSLON)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, half, one, two
  implicit none
  real(fp) :: R8_ZEROIN
  real(fp), parameter :: THREE = 3.0_fp
  real(fp) :: epslon,zansr
  real(fp) :: AX,BX,F,TOL
  external F
  !
  !      A ZERO OF THE FUNCTION  F(X)  IS COMPUTED IN THE INTERVAL AX,BX .
  !
  !  INPUT..
  !
  !  AX     LEFT ENDPOINT OF INITIAL INTERVAL
  !  BX     RIGHT ENDPOINT OF INITIAL INTERVAL
  !  F      FUNCTION SUBPROGRAM WHICH EVALUATES F(X) FOR ANY X IN
  !         THE INTERVAL  AX,BX
  !  TOL    DESIRED LENGTH OF THE INTERVAL OF UNCERTAINTY OF THE
  !         FINAL RESULT ( .GE. 0.0)
  !
  !
  !  OUTPUT..
  !
  !  R8_ZEROIN ABCISSA APPROXIMATING A ZERO OF  F  IN THE INTERVAL AX,BX
  !
  !
  !      IT IS ASSUMED  THAT   F(AX)   AND   F(BX)   HAVE  OPPOSITE  SIGNS
  !  WITHOUT  A  CHECK.  R8_ZEROIN  RETURNS A ZERO  X  IN THE GIVEN INTERVAL
  !  AX,BX  TO WITHIN A TOLERANCE  4*MACHEPS*ABS(X) + TOL, WHERE MACHEPS
  !  IS THE RELATIVE MACHINE PRECISION.
  !      THIS FUNCTION SUBPROGRAM IS A SLIGHTLY  MODifIED  TRANSLATION  OF
  !  THE ALGOL 60 PROCEDURE  ZERO  GIVEN IN  RICHARD BRENT, ALGORITHMS FOR
  !  MINIMIZATION WITHOUT DERIVATIVES, PRENTICE - HALL, INC. (1973).
  !
  !  dmc 10 Nov 1998 -- added a test to enforces AX.le.R8_ZEROIN.le.BX
  !    or AX.ge.R8_ZEROIN.ge.BX on exit
  !
  !
  real(fp) :: A,B,C,D,E,EPS,FA,FB,FC,TOL1,XM,P,Q,R,S
  !
  !  COMPUTE EPS, THE RELATIVE MACHINE PRECISION
  !
  !      EPS = 1.0
  !   10 EPS = EPS/2.0
  !      TOL1 = 1.0 + EPS
  !      if (TOL1 .GT. 1.0) GO TO 10
  EPS = EPSLON
  !
  ! INITIALIZATION
  !
  A = AX
  B = BX
  FA = F(A)
  FB = F(B)
  !
  ! BEGIN STEP
  !
20 C = A
  FC = FA
  D = B - A
  E = D
30 if (ABS(FC) .GE. ABS(FB)) GO TO 40
  A = B
  B = C
  C = A
  FA = FB
  FB = FC
  FC = FA
  !
  ! CONVERGENCE TEST
  !
40 TOL1 = TWO*EPS*ABS(B) + HALF*TOL
  XM = HALF*(C - B)
  if (ABS(XM) .LE. TOL1) GO TO 90
  if (FB .EQ. ZERO) GO TO 90
  !
  ! IS BISECTION NECESSARY
  !
  if (ABS(E) .LT. TOL1) GO TO 70
  if (ABS(FA) .LE. ABS(FB)) GO TO 70
  !
  ! IS QUADRATIC INTERPOLATION POSSIBLE
  !
  if (A .NE. C) GO TO 50
  !
  ! LINEAR INTERPOLATION
  !
  S = FB/FA
  P = TWO*XM*S
  Q = ONE- S
  GO TO 60
  !
  ! INVERSE QUADRATIC INTERPOLATION
  !
50 Q = FA/FC
  R = FB/FC
  S = FB/FA
  P = S*(TWO*XM*Q*(Q - R) - (B - A)*(R - ONE))
  Q = (Q - ONE)*(R - ONE)*(S - ONE)
  !
  ! ADJUST SIGNS
  !
60 if (P .GT. ZERO) Q = -Q
  P = ABS(P)
  !
  ! IS INTERPOLATION ACCEPTABLE
  !
  if ((TWO*P) .GE. (THREE*XM*Q - ABS(TOL1*Q))) GO TO 70
  if (P .GE. ABS(HALF*E*Q)) GO TO 70
  E = D
  D = P/Q
  GO TO 80
  !
  ! BISECTION
  !
70 D = XM
  E = D
  !
  ! COMPLETE STEP
  !
80 A = B
  FA = FB
  if (ABS(D) .GT. TOL1) B = B + D
  if (ABS(D) .LE. TOL1) B = B + SIGN(TOL1, XM)
  FB = F(B)
  if ((FB*(FC/ABS(FC))) .GT. ZERO) GO TO 20
  GO TO 30
  !
  ! DONE
  !
90 continue
  zansr = B                    ! dmc -- guarantee answer is in range
  zansr = max(zansr,min(ax,bx))
  zansr = min(zansr,max(ax,bx))
  R8_ZEROIN = zansr
  !
  RETURN
END FUNCTION R8_ZEROIN

