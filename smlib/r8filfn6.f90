!******************** START FILE FILFN6.FOR ; GROUP FILTR6 ******************
!
!  FILFN6
!
!  COMPUTE WEIGHTED AVERAGE ABOUT PT X(J0).
!
!  END CONDITIONS:
!     XEND?=0==> FIX END PT. OF REGION BEING SMOOTHED.
!                THE DATA CURVE IS REFLECTED ABOUT THE PT.
!                [X(1),Y(1)] ( OR [X(N),Y(N)])       (E0)
!               SO THAT BY SYMMETRY OF THE WEIGHTING
!                FUNCTION ==> THE END PTS ARE FIXED
!
!      XEND?=1==> DON'T FIX END PT. OF REGION. END PT
!                IS RESET TO WEIGHTED AVERAGE OF NEARBY PTS.
!                THE DATA IS REFLECTED ABOUT THE LINE X=X(1)
!                FOR THE LEFT END, ABOUT THE LINE X=X(N)
!                FOR  THE RIGHT END.                 (E1)
!
!
!  IF XEND DOES NOT EQUAL 0 OR 1 XEND SHOULD BE BETWEEN 0 AND 1;
!   THEN THE 'EXTENSION FORMULA' IS E=(1-XEND)*E0+XEND*E1
!
FUNCTION r8filfn6(X,Y,N,J0,DELTA0,XEND1,XEND2,ICOD)
  !
  !  DMC NOV 1985 - ICOD=1 : NORMAL TRIANGULAR WEIGHTED AVG INTEGRAL
  !
  !    ICOD=2 : BROKEN TWO PIECE TRIANGLE FOR SAWTOOTH CORRELATION
  !             INTEGRAL
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: ZERO, HALF, ONE
  implicit none
 real(fp) :: r8filfn6
 ! Arguments
  integer, intent(in) :: n, j0
  real(fp), dimension(n), intent(in) :: x
  real(fp), dimension(n), intent(in) :: y
  real(fp), dimension(n), intent(in) :: delta0
  real(fp), intent(in) :: xend1, xend2
  integer, intent(in) :: icod
  ! Local variables:
  integer :: j,jr
  real(fp) :: delta,xlim1,xlim2,denom, &
       zsign,wasum,xa,ya,xb,yb,ya0,ya1,xa1,yb0,yb1,xb1

  real(fp), external :: r8wxcint
  real(fp), external :: r8wxfint
  !
  !  CHECK FOR NUMERICALLY ZERO OR NEGATIVE DELTA
  DELTA=min((X(N)-X(1)),DELTA0(J0))
  IF((DELTA.LE.ZERO).OR.((X(J0)+DELTA).EQ.X(J0))) THEN
    IF(ICOD.EQ.1) THEN
      r8filfn6=Y(J0)
    ELSE
      r8filfn6=ZERO
    ENDIF
    RETURN
  ENDIF
  !
  IF(ICOD.EQ.1) THEN
    !  CHECK FOR FIXED ENDPOINTS
    r8filfn6=Y(J0)
    IF((J0.EQ.1).AND.(XEND1.EQ.ZERO)) RETURN
    IF((J0.EQ.N).AND.(XEND2.EQ.ZERO)) RETURN
  ENDIF
  !
  !  EVALUATE DENOMINATOR OF WEIGHTED AVERAGE FORMULA
  !
  XLIM1=X(J0)-DELTA
  XLIM2=X(J0)+DELTA
  !
  DENOM=r8wxcint(XLIM1,XLIM2,ONE,X(J0),DELTA)
  !
  !  EVALUATE NUMERATOR OF WGTED. AVG.
  !
  !  (A) PTS TO THE LEFT OF X(J0)
  !
  ZSIGN=ONE
  WASUM=ZERO
  !
  J=J0
  XA=X(J0)
  YA=Y(J0)
10 J=J-1
  XB=XA
  YB=YA
  IF (J.LE.0) THEN
  !
  !  REFLECT PTS ABOUT X(1),Y(1)
  !
     JR=min(N,(2-J))
     IF((X(J0)-XB).GE.DELTA.or.(J.LE.(1-N))) GO TO 50
     XA=X(1)-(X(JR)-X(1))
     YA0=Y(1)-(Y(JR)-Y(1))
     YA1=Y(JR)
     YA=(ONE-XEND1)*YA0+XEND1*YA1
  ELSE
     IF((X(J0)-XB).GE.DELTA) GO TO 50
     XA=X(J)
     YA=Y(J)
  ENDIF
  !
  !  COMPUTE WEIGHTED AVERAGE INTEGRAL FROM XA1 TO XB
  !
  XA1=max(XA,(X(J0)-DELTA))
  WASUM=WASUM+ZSIGN*r8wxfint(XA1,XB,XA,YA,XB,YB,X(J0),DELTA)
  !
  GO TO 10
!
  ! (B) COMPUTE CONTRIBUTION OF DATA TO RIGHT OF X(J0)
  !
50 CONTINUE
  !
  IF(ICOD.EQ.2) ZSIGN=-ONE
  !
  J=J0-1
  XB=X(J0)
  YB=Y(J0)
60 J=J+1
  XA=XB
  YA=YB
  IF(N.LE.J) THEN
  !
  !  REFLECT ABOUT X(N),Y(N)
  !
     JR=max(1,(N+N-J-1))
     IF((XA-X(J0)).GE.DELTA.or.(J.GE.(2*N-1))) GO TO 100
     XB=X(N)+(X(N)-X(JR))
     YB0=Y(N)+(Y(N)-Y(JR))
     YB1=Y(JR)
     YB=(ONE-XEND2)*YB0+XEND2*YB1
  ELSE
     IF((XA-X(J0)).GE.DELTA) GO TO 100
     XB=X(J+1)
     YB=Y(J+1)
  ENDIF
  !
  !  INTEGRATE FROM XA TO XB1
  !
  XB1=min(XB,(X(J0)+DELTA))
  !
  WASUM=WASUM+ZSIGN*r8wxfint(XA,XB1,XA,YA,XB,YB,X(J0),DELTA)
  !
  GO TO 60
  !
  !
  !  RETURN COMPUTED WEIGHTED AVERAGE
  !
100 CONTINUE
  !
  r8filfn6=WASUM/DENOM
  !
  RETURN
END FUNCTION r8filfn6
