!******************** START FILE FILTER.FOR ; GROUP FILTR6 ******************
!
!  FILTER
!
!  USE LOCAL AVERAGING FUNCTION OF RADIUS DELTA
!
!  TO SMOOTH DATA CURVE.
!
!-----------------
SUBROUTINE r8filter(X,Y0,Y,N,DELTA,XEND1,XEND2)
  use iso_c_binding, only: fp => c_double
  implicit none
  ! Arguments
  integer, intent(in) :: n
  real(fp), dimension(n), intent(in)  :: x
  real(fp), dimension(n), intent(in)  :: y0
  real(fp), dimension(n), intent(out) :: y
  real(fp), dimension(n), intent(in)  :: delta
  real(fp), intent(in) :: xend1, xend2
  ! Local variables
  integer :: j,jlp
  real(fp), external :: r8filfn6
  !
  DO J=1,N
    JLP=J
    Y(J)=r8filfn6(X,Y0,N,JLP,DELTA,XEND1,XEND2,1)
  END DO
  !
  RETURN
END SUBROUTINE r8filter
