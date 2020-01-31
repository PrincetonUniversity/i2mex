subroutine i2mex_getP(ns, psi, p, ier)

  ! Get the pressure p on the psi mesh (p in [mu0 * Pa]=[Tesla] units)

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: p(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns)
  integer iok

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_p, 0, p, iok)
  if(iok/=0) ier = 13

end subroutine i2mex_getP

subroutine i2mex_getPp(ns, psi, pp, ier)

  ! Get d p /d psi on the psi mesh. psi is the poloidal flux/(2*pi)

  use i2mex_mod
  implicit none
  integer :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8) :: pp(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns), dsdpsi(ns)
  integer iok


  include 'cubinterp.h'

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)
  call i2mex_getDS(ns, psi, dsdpsi, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_p, 1, pp, iok)
  if(iok/=0) ier = 14

  pp(2:ns) = dsdpsi(2:ns) * pp(2:ns)
  pp(1) = FCCCC0(pp(2), pp(3), pp(4), pp(5), &
       &                 s(2),  s(3),  s(4),  s(5), &
       & s(1))

end subroutine i2mex_getPp


subroutine i2mex_getPpp(ns, psi, ppp, ier)

  ! Get d^2 p /d psi^2 on psi mesh

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: ppp(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns), dsdpsi(ns), d2sdpsi2(ns), pp(ns)
  integer iok


  include 'cubinterp.h'

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)
  call i2mex_getDS(ns, psi, dsdpsi, iok)
  call i2mex_error(iok)
  call i2mex_getD2S(ns, psi, d2sdpsi2, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_p, 1, pp, iok)
  if(iok/=0) ier = 18
  call eq_rgetf(ns, s, i2mex_o%id_p, 2, ppp, iok)
  if(iok/=0) ier = 18

  if(psi(1) >  i2mex_o%psi_axis) then
     ppp(1:ns) = d2sdpsi2(1:ns) * pp(1:ns) + dsdpsi(1:ns)**2 * ppp(1:ns)
      
  else
     ppp(2:ns) = d2sdpsi2(2:ns) * pp(2:ns) + dsdpsi(2:ns)**2 * ppp(2:ns)
     if(ns >=5) then
        ppp(1) = FCCCC0(ppp(2), ppp(3), ppp(4), ppp(5), &
       &                 s(2),  s(3),  s(4),  s(5), &
       & s(1))
     else
        ppp(1) = ppp(2)
     endif
  endif

  call i2mex_toAxisFraction(1, ns, (/0._i2mex_r8/), psi, 0.0_i2mex_r8, &
       & i2mex_to_axis_fraction, ppp, iok)
  call i2mex_toEdgeFraction(1, ns, (/0._i2mex_r8/), psi,  &
       & i2mex_to_edge_fraction, ppp, iok)

end subroutine i2mex_getPpp
