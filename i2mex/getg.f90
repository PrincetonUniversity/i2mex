subroutine i2mex_getG(ns, psi, g, ier)

  ! Get g, the covariant toroidal magnetic field component, on the psi mesh

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: g(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns)
  integer iok

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_g, 0, g, iok)
  if(iok/=0) ier = 19

end subroutine i2mex_getG

subroutine i2mex_getGp(ns, psi, gp, ier)

  ! Get d g /d psi on psi mesh

  use i2mex_mod
  implicit none
  integer :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8) :: gp(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns), dsdpsi(ns)
  integer iok


  include 'cubinterp.h'

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)
  call i2mex_getDS(ns, psi, dsdpsi, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_g, 1, gp, iok)
  if(iok/=0) ier = 20

  gp(2:ns) = dsdpsi(2:ns) * gp(2:ns)
  gp(1) = FCCCC0(gp(2), gp(3), gp(4), gp(5), &
       &                 s(2),  s(3),  s(4),  s(5), &
       & s(1))

end subroutine i2mex_getGp

subroutine i2mex_getGpp(ns, psi, gpp, ier)

  ! Get d^2 g /d psi^2 on psi mesh

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: gpp(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns), dsdpsi(ns), d2sdpsi2(ns), gp(ns)
  integer iok


  include 'cubinterp.h'

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)
  call i2mex_getDS(ns, psi, dsdpsi, iok)
  call i2mex_error(iok)
  call i2mex_getD2S(ns, psi, d2sdpsi2, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_g, 1, gp, iok)
  if(iok/=0) ier = 18
  call eq_rgetf(ns, s, i2mex_o%id_g, 2, gpp, iok)
  if(iok/=0) ier = 18

  if(psi(1) >  i2mex_o%psi_axis) then
     gpp(1:ns) = d2sdpsi2(1:ns) * gp(1:ns) + dsdpsi(1:ns)**2 * gpp(1:ns)
      
  else
     gpp(2:ns) = d2sdpsi2(2:ns) * gp(2:ns) + dsdpsi(2:ns)**2 * gpp(2:ns)
     if(ns >=5) then
        gpp(1) = FCCCC0(gpp(2), gpp(3), gpp(4), gpp(5), &
       &                 s(2),  s(3),  s(4),  s(5), &
       & s(1))
     else
        gpp(1) = gpp(2)
     endif
  endif

  call i2mex_toAxisFraction(1, ns, (/0._i2mex_r8/), psi, 0.0_i2mex_r8, &
       & i2mex_to_axis_fraction, gpp, iok)
  call i2mex_toEdgeFraction(1, ns, (/0._i2mex_r8/), psi,  &
       & i2mex_to_edge_fraction, gpp, iok)

end subroutine i2mex_getGpp
