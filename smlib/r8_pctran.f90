!******************** START FILE R8_PCTRAN.FOR
!
!  R8_PCTRAN   DO PCWISE LINEAR INTERPOLATION TO PROGRAM COORDINATES
!
!  ARGS:
!
!  X  IN    X VALUES TO WHICH TO INTERPOLATE
!  Y  OUT   Y VALUES INTERPOLATED
!  N   IN   DIMENSIONALITY OF X,Y
!
!  XPC  IN  X VALUES FROM WHICH TO INTERPOLATE
!  YPC  IN  Y VALUES FROM WHICH TO INTERPOLATE
!  NPC  IN  # OF X,Y, VALUES
!  MAXPC  IN  DIMENSIONALITY OF XPC,YPC
!
SUBROUTINE R8_PCTRAN(X,Y,N,XPC,YPC,MAXPC,NPC)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: ZERO, HALF, ONE
  implicit none
  ! Arguments
  integer, intent(in) :: n, maxpc, npc
  real(fp), dimension(n), intent(in)    :: x
  real(fp), dimension(n), intent(inout) :: y
  real(fp), dimension(maxpc), intent(in) :: xpc, ypc
  ! Local variables
  real(fp) :: zy
  integer :: i, npm1, j, jp1
  !
  IF(NPC.GT.1) GO TO 10
  IF(NPC.LE.0) ZY=ZERO
  IF(NPC.GT.0) ZY=YPC(1)
  DO I=1,N
    Y(I)=ZY
  END DO
  RETURN
!
10 CONTINUE
!
  DO 15 I=1,N
    IF(X(I).GT.XPC(1)) GO TO 12
    Y(I)=YPC(1)
    GO TO 15
12  CONTINUE
    IF(X(I).LT.XPC(NPC)) GO TO 13
    Y(I)=YPC(NPC)
    GO TO 15
13  CONTINUE
    NPM1=NPC-1
    DO 14 J=1,NPM1
      IF((X(I).GE.XPC(J)).AND.(X(I).LT.XPC(J+1))) GO TO 16
14  CONTINUE
16  CONTINUE
    JP1=J+1
    Y(I)=YPC(J)+(X(I)-XPC(J))*(YPC(JP1)-YPC(J))/(XPC(JP1)-XPC(J))
15 CONTINUE
  RETURN
END SUBROUTINE R8_PCTRAN
 
