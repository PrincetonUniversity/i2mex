!******************** START FILE IERFCN.FOR ; GROUP FILTR6 ******************
!
!  IERFCN   TEST Y(J)-YSM(J)
!   IERFCN=1==> EPSLON VIOLATION
!   IERFCN=0==> NO VIOLATION
!
!   ISIGN=SIGN OF VIOLATION
!   DY=MAGNITUDE OF VIOLATION
!
INTEGER FUNCTION r8ierfcn(Y,YSM,N,JP,EPS,NE,DY,ISIGN)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: ZERO, HALF, ONE
  implicit none
  ! Arguments
  integer, intent(in) :: n, jp, ne
  real(fp), dimension(n), intent(in) :: y
  real(fp), dimension(n), intent(in) :: ysm
  real(fp), dimension(ne), intent(in) :: eps
  real(fp), intent(out) :: dy
  integer, intent(out) :: isign
  ! Local variables
  real(fp) :: zeps
!
  r8ierfcn=0
  ZEPS=EPS(NE)
  IF(JP.LT.NE) ZEPS=EPS(JP)
  DY=Y(JP)-YSM(JP)
  ISIGN=SIGN(ONE,DY)
  DY=ABS(DY)
  IF(DY.LE.ZEPS) RETURN
  r8ierfcn=1
  RETURN
END FUNCTION r8ierfcn
