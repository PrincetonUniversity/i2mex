!******************** START FILE FILERS.FOR ; GROUP FILTR6 ******************
!
!  FILERS    SEEK OUT NEXT REGION WHERE SMOOTHED DATA VIOLATES
!            EPSLON ERROR BARS ON ORIGINAL DATA
!
!  JS=LOCATION AT WHICH TO START SEARCH
!
!  IND(6)  INDEX DESCRIPTION OF ERROR  IND(1)=0 ==> NO ERROR
!          FOUND  (OUTPUT)
!          IND(3),IND(4) DELIMIT ACTUAL VIOLATIONS WITHIN A
!          REGION WHERE THE DIVERGENCE (Y-YSM) IS OF THE SAME
!          SIGN   (X,Y (DIM(N)) ARE RAW DATA, YSM IS SMOOTHED DATA)
!          IND(2),IND(5) DELIMIT REGION WHERE (Y-YSM) IS OF THE
!          SAME SIGN
!          IND(1),IND(6) DELIMIT REGION WITHIN DELTA(IND(2))
!            OF X(IND(3)) (ON THE RIGHT); AND DELTA(IND(5))
!            OF X(IND(4)) (ON THE LEFT) ;  THIS IS THE REGION
!            WHERE A CORRECTION WILL BE APPLIED
!
!  EPS(NE) ARRAY CONTAINS ERROR BARS FOR POINTS 1 THRU NE; IF
!          J>NE, EPS(NE) IS USED AS THE ERROR BAR FOR
!          THE JTH DATA POINT.
!
!   EPS2(NE2):  IF A REGION OF SAME SIGN DIVERGENCE (Y-YSM)
!              CONTAINS ND LESS THAN NE2 DATA POINTS, THEN THE
!            ERROR BAR WILL BE MULTIPLIED BY THE FACTOR
!              EPS2(ND) ... THIS ALLOWS A LESS STIFF ERROR
!              BAR CRITERION FOR VIOLATIONS OF SHORT WIDTH
!
SUBROUTINE r8filers(X,Y,YSM,N,DELTA,EPS,NE,EPS2,NE2,IND,JS)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: ZERO, HALF, ONE
  implicit none
  ! Arguments
  integer, intent(in) :: N, NE, NE2, js
  real(fp), dimension(n), intent(in)  :: x, y
  real(fp), dimension(ne), intent(in) :: eps
  real(fp), dimension(1), intent(in)  :: eps2
  real(fp), dimension(n), intent(out) :: ysm, delta
  integer, dimension(6), intent(out) :: ind
  ! Local variables
  integer :: j,jlp,ierrf,isign0,je1,je2,js2,jd2,js1,jd1, &
       j2,jlp2,nse,jst,isign
  real(fp) :: dymax,dy,x0,zemax,zetes
  !
  integer, external :: r8ierfcn
  !
  !  CLEAR ERROR FLAG
  !
  IND(1)=0
  !
  IF(JS.GT.N) RETURN
  !  FIND FIRST EPSLON VIOLATION
  DO J=JS,N
    DYMAX=ZERO
    JLP=J
    IERRF=r8ierfcn(Y,YSM,N,JLP,EPS,NE,DY,ISIGN0)
    IF(IERRF.EQ.0) cycle
    !  POSSIBLE ERROR DETECTED
    DYMAX=max(DYMAX,DY)
    IF(J.LT.N) GO TO 10
    JE1=N
    JE2=N
    JS2=N
    JD2=N
    GO TO 40
    !
10  CONTINUE
    JE1=J
    IF(JE1.GT.1) GO TO 20
    JS1=1
    JD1=1
    !  FIND  JE2=IND(4),  JS2=IND(5)
20  CONTINUE
    JE2=JE1
    DO J2=JE1,N
      JLP2=J2
      IERRF=r8ierfcn(Y,YSM,N,JLP2,EPS,NE,DY,ISIGN)
      DYMAX=max(DYMAX,DY)
      IF(ISIGN.NE.ISIGN0) GO TO 35
      IF(IERRF.EQ.1) JE2=J2
    END DO
    J2=N+1
35  CONTINUE
    JS2=J2-1
    !  FIND JD2 (IND(6))
    X0=X(JE2)
    DO J2=JE2,N
      IF((X(J2)-X0).GT.DELTA(J2)) GO TO 38
    END DO
    J2=N
38  CONTINUE
    JD2=J2
!
    IF(JE1.EQ.1) GO TO 60
!  FIND JS1  (IND(2))
40  CONTINUE
    DO J2=JE1-1,1,-1
      JLP2=J2
      IERRF=r8ierfcn(Y,YSM,N,JLP2,EPS,NE,DY,ISIGN)
      IF(ISIGN.NE.ISIGN0) GO TO 55
    END DO
    J2=0
55  CONTINUE
    JS1=J2+1
!  FIND JD1 (IND(1))
    X0=X(JE1)
    DO J2=JE1,1,-1
      IF((X0-X(J2)).GT.DELTA(J2)) GO TO 58
    END DO
    J2=1
58  CONTINUE
    JD1=J2
    !  DOES EPS2 ARRAY COME TO BEAR ?
60  CONTINUE
    NSE=JS2-JS1+1
    IF(NSE.GT.NE2) GO TO 80
    !  MAYBE
    ZEMAX=EPS(NE)*EPS2(NSE)
    IF(JE1.GE.NE) GO TO 75
    JST=min(JE2,NE)
    ZEMAX=ZERO
    DO J2=JE1,JST
      ZETES=EPS(J2)*EPS2(NSE)
      ZEMAX=max(ZEMAX,ZETES)
    END DO
75  CONTINUE
    IF(ZEMAX.GE.DYMAX) cycle
    !  ERROR DEFINATELY DETECTED
80  CONTINUE
    IND(1)=JD1
    IND(2)=JS1
    IND(3)=JE1
    IND(4)=JE2
    IND(5)=JS2
    IND(6)=JD2
    GO TO 200
    !
    !  END OF LOOPT
    !
  END DO
200 CONTINUE
  RETURN
END SUBROUTINE r8filers
