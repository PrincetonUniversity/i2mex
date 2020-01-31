subroutine i2mex_getPhi(ns, psi, res, ier)
  
  ! Return the toroidal flux/(2*pi) in res.
  ! The toroidal flux/(2*pi) is the integral of q d psi

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier

  
  real(i2mex_r8) q(ns)
  integer iok

  ier = 0
  call i2mex_getQ(ns, psi, q, iok)
  call i2mex_integrate1d(ns, psi, q, ns, psi, res)

  if(iok/=0) ier = 110

end subroutine i2mex_getPhi

subroutine i2mex_getPhiP(ns, psi, res, ier)
  
  ! Return d (toroidal flux/(2*pi))/d psi in res
  ! This is just the q profile

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier

  integer iok

  ier = 0
  call i2mex_getQ(ns, psi, res, iok)

  if(iok/=0) ier = 111

end subroutine i2mex_getPhiP

subroutine i2mex_getPhiPP(ns, psi, res, ier)
  
  ! Return d^2 (toroidal flux/(2*pi))/d psi^2 in res
  ! This is just the d q/ d psi

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier

  integer iok

  ier = 0
  call i2mex_getQP(ns, psi, res, iok)

  if(iok/=0) ier = 112

end subroutine i2mex_getPhiPP
