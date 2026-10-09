!******************** START FILE ZEROIN.FOR ; GROUP TRKRLIB ************
!C
      REAL FUNCTION ZEROIN(AX,BX,F,TOL,EPSLON)
!
      REAL AX,BX,F,TOL
      EXTERNAL F
!
!      A ZERO OF THE FUNCTION  F(X)  IS COMPUTED IN THE INTERVAL AX,BX .
!
!  INPUT..
!
!  AX     LEFT endPOINT OF INITIAL INTERVAL
!  BX     RIGHT endPOINT OF INITIAL INTERVAL
!  F      FUNCTION SUBPROGRAM WHICH EVALUATES F(X) FOR ANY X IN
!         THE INTERVAL  AX,BX
!  TOL    DESIRED LENGTH OF THE INTERVAL OF UNCERTAINTY OF THE
!         FINAL RESULT ( .GE. 0.0)
!
!
!  OUTPUT..
!
!  ZEROIN ABCISSA APPROXIMATING A ZERO OF  F  IN THE INTERVAL AX,BX
!
!
!      IT IS ASSUMED  THAT   F(AX)   AND   F(BX)   HAVE  OPPOSITE  SIGNS
!  WITHOUT  A  CHECK.  ZEROIN  returnS A ZERO  X  IN THE GIVEN INTERVAL
!  AX,BX  TO WITHIN A TOLERANCE  4*MACHEPS*ABS(X) + TOL, WHERE MACHEPS
!  IS THE RELATIVE MACHINE PRECISION.
!      THIS FUNCTION SUBPROGRAM IS A SLIGHTLY  MODifIED  TRANSLATION  OF
!  THE ALGOL 60 PROCEDURE  ZERO  GIVEN IN  RICHARD BRENT, ALGORITHMS FOR
!  MINIMIZATION WITHOUT DERIVATIVES, PRENTICE - HALL, INC. (1973).
!
!  dmc 10 Nov 1998 -- added a test to enforces AX.le.ZEROIN.le.BX
!    or AX.ge.ZEROIN.ge.BX on exit
!
!
      REAL  A,B,C,D,E,EPS,FA,FB,FC,TOL1,XM,P,Q,R,S
!
!  COMPUTE EPS, THE RELATIVE MACHINE PRECISION
!
!      EPS = 1.0
!   10 EPS = EPS/2.0
!      TOL1 = 1.0 + EPS
!      if (TOL1 .GT. 1.0) goto 10
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
   30 if (ABS(FC) .GE. ABS(FB)) goto 40
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
!
! CONVERGENCE TEST
!
   40 TOL1 = 2.0*EPS*ABS(B) + 0.5*TOL
      XM = .5*(C - B)
      if (ABS(XM) .LE. TOL1) goto 90
      if (FB .EQ. 0.0) goto 90
!
! IS BISECTION NECESSARY
!
      if (ABS(E) .LT. TOL1) goto 70
      if (ABS(FA) .LE. ABS(FB)) goto 70
!
! IS QUADRATIC INTERPOLATION POSSIBLE
!
      if (A .NE. C) goto 50
!
! LINEAR INTERPOLATION
!
      S = FB/FA
      P = 2.0*XM*S
      Q = 1.0 - S
      goto 60
!
! INVERSE QUADRATIC INTERPOLATION
!
   50 Q = FA/FC
      R = FB/FC
      S = FB/FA
      P = S*(2.0*XM*Q*(Q - R) - (B - A)*(R - 1.0))
      Q = (Q - 1.0)*(R - 1.0)*(S - 1.0)
!
! ADJUST SIGNS
!
   60 if (P .GT. 0.0) Q = -Q
      P = ABS(P)
!
! IS INTERPOLATION ACCEPTABLE
!
      if ((2.0*P) .GE. (3.0*XM*Q - ABS(TOL1*Q))) goto 70
      if (P .GE. ABS(0.5*E*Q)) goto 70
      E = D
      D = P/Q
      goto 80
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
      if ((FB*(FC/ABS(FC))) .GT. 0.0) goto 20
      goto 30
!
! doNE
!
   90 continue
      zansr = B                    ! dmc -- guarantee answer is in range
      zansr = max(zansr,min(ax,bx))
      zansr = min(zansr,max(ax,bx))
      ZEROIN = zansr
!
      return
      end
!******************** end FILE ZEROIN.FOR ; GROUP TRKRLIB **************
