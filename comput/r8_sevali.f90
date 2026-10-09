!******************** START FILE SEVALI.FOR ; GROUP TRKRLIB *************
!  MOD DMC SUMMER 1987 JET/GARCHING
!   ROUTINE MAY EVALUATE LINEAR INSTEAD OF SPLINE INTERPOLATION, if
!   COMMON SWITCH ILIN IS SET ***
!....................................................
FUNCTION r8_sevali(N, U, X, Y, B, C, D, dx)
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) :: r8_sevali
  integer :: iput,ilin
  integer :: N
  real(fp) :: U
  real(fp), dimension(N) :: X, Y, B, C, D
  real(fp) :: dx
  !
  !  THIS subroutine EVALUATES THE CUBIC SPLINE FUNCTION
  !
  !    SEVALI = Y(I) + B(I)*(U-X(I)) + C(I)*(U-X(I))**2 + D(I)*(U-X(I))**3
  !
  !    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
  !
  !  if  U .LT. X(1) THEN  I = 1  IS USED.
  !  if  U .GE. X(N) THEN  I = N  IS USED.
  !
  !  INPUT..
  !
  !    N = THE NUMBER OF DATA POINTS
  !    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
  !    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
  !    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE
  !
  !  if  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
  !  BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
  !
  !  DMC GARCHING OCT 1985 KLUGE
  COMMON/ZSEVALI/ IPUT,ILIN
  !  RETURN ZONE INDEX TO CALLER VIA THIS COMMON BLOCK ***
  !
  integer :: I, J, K
  real(fp) :: zslop
  DATA I/1/
  SAVE I

  if ( I .GE. N ) I = 1
  if ( U .LT. X(I) ) GO TO 10
  if ( U .LE. X(I+1) ) GO TO 30
  !
  !  BINARY SEARCH
  !
10 I = 1
  J = N+1
20 K = (I+J)/2
  if ( U .LT. X(K) ) J = K
  if ( U .GE. X(K) ) I = K
  if ( J .GT. I+1 ) GO TO 20
  !
  !  EVALUATE SPLINE
  !
30 DX = U - X(I)
  if(ILIN.EQ.0) then
    r8_sevali = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
  else
    if(I.EQ.N) then
      ZSLOP=(Y(N)-Y(N-1))/(X(N)-X(N-1))
    else
      ZSLOP=(Y(I+1)-Y(I))/(X(I+1)-X(I))
    end if
    r8_sevali=Y(I)+DX*ZSLOP
  end if
  !
  IPUT=I
  !
  return
END FUNCTION r8_sevali
