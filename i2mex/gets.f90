subroutine i2mex_getS(ns, psi, s, ier)

  ! s ~ sqrt(psi)
  ! This is the radial coordinate as used internally in i2MEX

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(*)
  real(i2mex_r8), intent(out) :: s(*)
  integer, intent(out) :: ier

  ier=0

!!$  if(psi(1)==0.0_i2mex_r8 .and. psi(ns)==i2mex_o%psi_edge .and. ns>=5) then
     s(1:ns) = sqrt((psi(1:ns) - i2mex_o%psi_axis)/(i2mex_o%psi_edge - i2mex_o%psi_axis))
!!$  else
!!$     ier = 10
!!$  endif
end subroutine i2mex_getS

subroutine i2mex_getDS(ns, psi, ds, ier)

  ! d s/ d psi

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(*)
  real(i2mex_r8), intent(out) :: ds(*)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns)

  include 'cubinterp.h'

  ier=0
  s(1:ns) = sqrt((psi(1:ns) - i2mex_o%psi_axis)/(i2mex_o%psi_edge - i2mex_o%psi_axis))
  if(psi(1)/=0.0_i2mex_r8) then
     ds(1:ns) = 0.5_i2mex_r8 * s(1:ns) / psi(1:ns)
  else
     ds(2:ns) = 0.5_i2mex_r8 * s(2:ns) / psi(2:ns)
     if(ns >= 5) then
        ! cubic extrapolation to axis
        ds(1) = FCCCC0(ds(2), ds(3), ds(4), ds(5), &
             &                 s(2),  s(3),  s(4),  s(5), &
             & s(1))
     else if(ns >= 3) then
        ! linear extrapolation
        ds(1) = ds(2) - (s(2) - s(1))*(ds(3) - ds(2))/(s(3) - s(2))
     else if(ns >= 2) then
        ! const extrapolation
        ds(1) = ds(2)
     else
        ier = 11 
     endif
  endif
        
end subroutine i2mex_getDS

subroutine i2mex_getD2S(ns, psi, d2s, ier)

  ! d^2 s/ d psi^2

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(*)
  real(i2mex_r8), intent(out) :: d2s(*)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns)

  include 'cubinterp.h'

  ier=0
  s(1:ns) = sqrt((psi(1:ns) - i2mex_o%psi_axis)/(i2mex_o%psi_edge - i2mex_o%psi_axis))
  if(psi(1)/=0.0_i2mex_r8) then
     d2s(1:ns) = -0.25_i2mex_r8 * s(1:ns) /psi(1:ns)**2
  else
     d2s(2:ns) = -0.25_i2mex_r8 * s(2:ns) /psi(2:ns)**2
     if(ns >= 5) then
        ! cubic extrapolation to axis
        d2s(1) = FCCCC0(d2s(2), d2s(3), d2s(4), d2s(5), &
             &                   s(2),   s(3),   s(4),   s(5), &
             & s(1))
     else if(ns >= 3) then
        ! linear extrapolation
        d2s(1) = d2s(2) - (s(2) - s(1))*(d2s(3) - d2s(2))/(s(3) - s(2))
     else if(ns >= 2) then
        ! const extrapolation
        d2s(1) = d2s(2) 
     else
        ier = 12
     endif
  endif

end subroutine i2mex_getD2S
