!******************** START FILE WXCINT.FOR ; GROUP FILTR6 ******************
!
!  WXCINT
!
!  COMPUTE C*INTEGRAL OF WEIGHTING FUNCTION FROM X1 TO X2
!
!  DENOMINATOR IN WEIGHTED AVERAGE FORMULA
!
!  W(T)=1-ABS((XCEN-T)/DELTA)   XCEN-DELTA .LE. T .LE. XCEN+DELTA
!
!  TRIANGULAR SHAPED WEIGHTING FUNCTION CENTERED AT XCEN
!    X1 SHOULD BE .LT. X2 AND X1 AND X2 SHOULD
!    BE WITHIN DELTA OF XCEN
!
FUNCTION r8WXCINT(X1,X2,C,XCEN,DELTA)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: ZERO, HALF, ONE
  implicit none
  REAL(fp) :: r8WXCINT
  ! Arguments
  real(fp), intent(in) :: x1, x2, c
  real(fp), intent(in) :: xcen
  real(fp), intent(in) :: delta
  ! Local variables
  real(fp) :: xa, xb, zdi
  !
  r8wxcint=C
  XA=max(-DELTA,(X1-XCEN))
  XB=min(DELTA,(X2-XCEN))
  IF((XA.GT.DELTA).OR.(XB.LT.-DELTA)) GO TO 1000
  IF(XB.LT.XA) GO TO 1000
  !
  !  NUMERICS CHECK; DELTA <<< X1
  IF(XCEN.EQ.(XCEN+DELTA)) RETURN
  !
  ZDI=ONE/DELTA
  IF((XA.GT.ZERO).OR.(XB.LT.ZERO)) GO TO 100
  !
  r8wxcint=C*ZDI*((XB-XA)-HALF*(XA*XA+XB*XB)*ZDI)
  RETURN
  !
100 CONTINUE
  IF(XA.GT.ZERO) r8wxcint=C*ZDI*((XB-XA)-HALF*(XA+XB)*(XB-XA)*ZDI)
  IF(XB.LT.ZERO) r8wxcint=C*ZDI*((XB-XA)-HALF*(XA+XB)*(XA-XB)*ZDI)
  RETURN
  !
1000 CONTINUE
  write(6,9000)
9000 FORMAT(' SUBROUTINE WXCINT FROM R8FILFN6 FROM R8FILTR6.FOR:'/ &
          ' ILLEGAL INTEGRATION LIMITS')
  call flush(6)
  stop
END FUNCTION R8WXCINT
