subroutine integ_wts(x,w)

  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: half
  implicit none

  !  return x evaluation points and weights for 
  !  a definite integral over the region [0,1].
  !  Note the x evaluation points returned are not in ascending order.

  real(fp), dimension(10), intent(out) :: x
  real(fp), dimension(10), intent(out) :: w

  !----------------
  real(fp), dimension(5) :: x1
  real(fp), dimension(5) :: w10

  ! gauss-kronrod-patterson quadrature coefficients for use in
  ! quadpack routine qng.  these coefficients were calculated with
  ! 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.

  data x1    (  1) / 0.973906528517171720077964012084452_fp/
  data x1    (  2) / 0.865063366688984510732096688423493_fp/
  data x1    (  3) / 0.679409568299024406234327365114874_fp/
  data x1    (  4) / 0.433395394129247190799265943165784_fp/
  data x1    (  5) / 0.148874338981631210884826001129720_fp/

  data w10   (  1) / 0.066671344308688137593568809893332_fp/
  data w10   (  2) / 0.149451349150580593145776339657697_fp/
  data w10   (  3) / 0.219086362515982043995534934228163_fp/
  data w10   (  4) / 0.269266719309996355091226921569469_fp/
  data w10   (  5) / 0.295524224714752870173892994651338_fp/
  !----------------

  integer :: ii
  real(fp) :: xinc,ww

  do ii=1,5
    xinc = x1(ii)*HALF
    x(ii) = HALF + xinc
    x(5+ii) = HALF - xinc
    ww = w10(ii)*HALF
    w(ii) = ww
    w(5+ii) = ww
  end do

  return
end subroutine integ_wts
