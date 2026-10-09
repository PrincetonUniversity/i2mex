subroutine i2mex_integrate1d(ns, s, f, nsi, si, aifi)

  ! Perform the integration of a function f defined on
  ! the grid s of size ns assuming not-a-knot boundary condirions. 
  ! The result, aifi, is an array whose elements aifi(i) are the 
  ! integral from node 1 to node i, on the target mesh si of size nsi.

  ! A. Pletzer Oct 10 2000

  implicit none

  integer, parameter :: r8=selected_real_kind(12,100)
  integer, intent(in) :: ns
  real(r8), intent(in) :: s(ns)
  real(r8), intent(in) :: f(ns)
  integer, intent(in) :: nsi
  real(r8), intent(in) :: si(nsi)
  real(r8), intent(out) :: aifi(nsi)

  real(r8) wk(ns)
  real(r8) cspl(4, ns)

  cspl = 0.0_r8
  cspl(1, 1:ns) = f(1:ns)

  ! compute the spline coefficients
  call V_SPLINE(0, 0, ns, s, cspl, wk)

  ! integrate & interpolate
  call R8v_ieval(ns, s, nsi, si, cspl, aifi)

end subroutine i2mex_integrate1d

subroutine i2mex_integratePeriodic1d(nt1, t, f, nti1, ti, aifi)

  ! Perform the integration of a periodic function f defined
  ! over the grid t of size nt1: f(nt1)=f(1). The result, aifi, is an 
  ! array whose elements aifi(i) are the integral from node 1 to 
  ! node i, on the target mesh ti of size nti1.

  ! A. Pletzer Oct 10 2000

  implicit none

  integer, parameter :: r8=selected_real_kind(12,100)
  integer, intent(in) :: nt1
  real(r8), intent(in) :: t(nt1)
  real(r8), intent(in) :: f(nt1)
  integer, intent(in) :: nti1
  real(r8), intent(in) :: ti(nti1)
  real(r8), intent(out) :: aifi(nti1)

  real(r8) wk(nt1)
  real(r8) cspl(4, nt1)

  cspl = 0.0_r8
  cspl(1, 1:nt1) = f(1:nt1)

  ! compute the spline coefficients
  call V_SPLINE(-1, 0, nt1, t, cspl, wk)

  ! integrate & interpolate
  call R8v_ieval(nt1, t, nti1, ti, cspl, aifi)

end subroutine i2mex_integratePeriodic1d

subroutine R8v_ieval(N, X, ni, xi, cspl, aifi)
  IMPLICIT NONE
  integer, parameter :: r8=selected_real_kind(12,100)
  INTEGER, intent(in) :: N
  REAL(R8), intent(in) ::  X(N)
  INTEGER, intent(in) :: Ni
  REAL(R8), intent(in) ::  Xi(Ni)
  REAL(R8), intent(in) :: cspl(4,N)
  REAL(R8), intent(out) :: aifi(ni)

  integer i
  !
  !  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE INTEGRAL FROM x(1) TO U
  ! the spline coefficients cspl are those returned by v_spline
  !
  !    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
  !
  !  IF  U .LT. X(1) THEN  I = 1  IS USED.
  !  IF  U .GE. X(N) THEN  I = N  IS USED.
  !
  !  INPUT..
  !
  !    N = THE NUMBER OF DATA POINTS
  !    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
  !    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
  !    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE

  INTEGER nm1
  real(r8), dimension(:), allocatable :: dx, x1, x2

  nm1 = n - 1

  allocate(dx(nm1), x1(nm1), x2(nm1))

  x1 = x(1:n-1)
  x2 = x(2:n)

  do i = 1, ni
     dx = 0.0
     where (x2 < xi(i))
        dx = x2 - x1
     end where
     where ((x1 < xi(i))  .and. (xi(i) <= x2))
        dx = xi(i) - x1
     end where

     aifi(i) = sum( dx*( cspl(1, 1:nm1) &
          & + dx*( cspl(2, 1:nm1)/2.0_r8 &
          & + dx*( cspl(3, 1:nm1)/6.0_r8 &
          & + dx*  cspl(4, 1:nm1)/24.0_r8 ) ) ) )
  enddo

  deallocate(dx, x1, x2)

END subroutine R8v_ieval


subroutine i2mex_integrateAlongFluxSurface(f, nt1, ns, t, psi, aifi, ier)

  ! Perform the integration of a periodic function f defined
  ! on the *ORIGINAL* grid (t_ori, psi_ori). The result, aifi, is a
  ! 2-d array whose elements aifi(i, j) are integrals from node 1 to 
  ! node i, on the target mesh (t, psi) of size nt1*ns.

  ! A. Pletzer April 30 2001

  use i2mex_mod
  implicit none

  integer, parameter :: r8=selected_real_kind(12,100)
  real(r8), intent(in) :: f(i2mex_o%nt1, i2mex_o%ns)
  integer, intent(in) :: nt1
  integer, intent(in) :: ns
  real(r8), intent(in) :: t(nt1)
  real(r8), intent(in) :: psi(ns)
  real(r8), intent(out) :: aifi(nt1, ns)
  integer, intent(out) :: ier ! 0 = OK

  real(r8) :: soff
  real(r8) :: s(ns), zone_integral(nt1-1, ns)
  integer i, j, iok

  ier = 0

  call i2mex_getS(ns, psi, s, iok) 
  call i2mex_error(iok)
 
  soff = s(1)+1.e-6_r8*(s(2) - s(1))
  call eq_flxint_init(0, ns, s, soff, iok)
  call eq_flxint_chinit(0, nt1, t, iok)

  call eq_flxint_arr2(f, &
       &   1, zone_integral, nt1-1, ns, iok)

  ! add contribitions
  do j = 1, ns
     aifi(1, j) = 0.0_r8
     do i = 2, nt1
        aifi(i, j) = sum(zone_integral(1:i-1, j))/i2mex_twopi_r8
     enddo
  enddo


end subroutine i2mex_integrateAlongFluxSurface
