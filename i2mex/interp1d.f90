subroutine i2mex_interp1d(ns, s, f, nsi, si, fi)

  ! Return fi, the interpolation of the original set of ns data
  ! f @ s onto a target grid si of size nsi. Not-a-knot BCs are
  ! applied.

  use ezspline
  use ezspline_obj
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
  
  integer, intent(in) :: ns ! number of input data points
  real(r8), intent(in) :: s(ns) ! original grid
  real(r8), intent(in) :: f(ns) ! the function values at s
  integer, intent(in) :: nsi ! number of output data points
  real(r8), intent(in) :: si(ns) ! the target grid
  real(r8), intent(out) :: fi(ns) ! the interpolated data
  
  type(ezspline1) :: spline_o
  integer ier

  call EZspline_init(spline_o, ns, (/0,0/), ier)
  call EZspline_error(ier)

  spline_o%x1 = s

  call EZspline_setup(spline_o, f, ier)
  call EZspline_error(ier)

  call EZspline_interp(spline_o, nsi, si, fi, ier)  
  call EZspline_error(ier)

  call EZspline_free(spline_o, ier)
  call EZspline_error(ier)

end subroutine i2mex_interp1d
  
subroutine i2mex_interpPeriodic1d(nt1, t, f, nti1, ti, fi)

  ! Return fi, the interpolation of the original set of nt1 data
  ! f @ t onto a target grid ti of size nti1. Periodic BCs are
  ! applied: f(nt1) = f(1).

  use ezspline
  use ezspline_obj
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
  
  integer, intent(in) :: nt1 ! number of input data points
  real(r8), intent(in) :: t(nt1) ! original grid
  real(r8), intent(in) :: f(nt1) ! the function values at t
  integer, intent(in) :: nti1 ! number of output data points
  real(r8), intent(in) :: ti(nti1) ! the target grid
  real(r8), intent(out) :: fi(nti1) ! the interpolated data
  
  type(ezspline1) :: spline_o
  integer ier

  call EZspline_init(spline_o, nt1, (/-1,-1/), ier)
  call EZspline_error(ier)

  spline_o%x1 = t

  call EZspline_setup(spline_o, f, ier)
  call EZspline_error(ier)

  call EZspline_interp(spline_o, nti1, ti, fi, ier)  
  call EZspline_error(ier)

  call EZspline_free(spline_o, ier)
  call EZspline_error(ier)

end subroutine i2mex_interpPeriodic1d
  
