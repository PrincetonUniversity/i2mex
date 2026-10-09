subroutine i2mex_getQ(ns, psi, q, ier)

  ! Get q on psi mesh

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: q(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns)
  integer iok

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_q, 0, q, iok)
  if(iok/=0) ier = 16

end subroutine i2mex_getQ

subroutine i2mex_getQp(ns, psi, qp, ier)

  ! Get d q /d psi on psi mesh

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8) :: qp(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns), dsdpsi(ns)
  integer iok


  include 'cubinterp.h'

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)
  call i2mex_getDS(ns, psi, dsdpsi, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_q, 1, qp, iok)
  if(iok/=0) ier = 17

  if(psi(1) >  i2mex_o%psi_axis) then
      qp(1:ns) = dsdpsi(1:ns) * qp(1:ns)
  else
     qp(2:ns) = dsdpsi(2:ns) * qp(2:ns)
     if(ns >=5) then
     qp(1) = FCCCC0(qp(2), qp(3), qp(4), qp(5), &
       &                 s(2),  s(3),  s(4),  s(5), &
       & s(1))
     else
        qp(1) = qp(2)
     endif
  endif

  call i2mex_toAxisFraction(1, ns, (/0._i2mex_r8/), psi, 0.0_i2mex_r8, &
       & i2mex_to_axis_fraction, qp, ier)
  call i2mex_toEdgeFraction(1, ns, (/0._i2mex_r8/), psi,  &
       & i2mex_to_edge_fraction, qp, ier)

end subroutine i2mex_getQp


subroutine i2mex_getQpp(ns, psi, qpp, ier)

  ! Get d^2 q /d psi^2 on psi mesh

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: qpp(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns), dsdpsi(ns), d2sdpsi2(ns), qp(ns)
  integer iok


  include 'cubinterp.h'

  ier=0
  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)
  call i2mex_getDS(ns, psi, dsdpsi, iok)
  call i2mex_error(iok)
  call i2mex_getD2S(ns, psi, d2sdpsi2, iok)
  call i2mex_error(iok)

  call eq_rgetf(ns, s, i2mex_o%id_q, 1, qp, iok)
  if(iok/=0) ier = 18
  call eq_rgetf(ns, s, i2mex_o%id_q, 2, qpp, iok)
  if(iok/=0) ier = 18

  if(psi(1) >  i2mex_o%psi_axis) then
     qpp(1:ns) = d2sdpsi2(1:ns) * qp(1:ns) + dsdpsi(1:ns)**2 * qpp(1:ns)
      
  else
     qpp(2:ns) = d2sdpsi2(2:ns) * qp(2:ns) + dsdpsi(2:ns)**2 * qpp(2:ns)
     if(ns >=5) then
        qpp(1) = FCCCC0(qpp(2), qpp(3), qpp(4), qpp(5), &
       &                 s(2),  s(3),  s(4),  s(5), &
       & s(1))
     else
        qpp(1) = qpp(2)
     endif
  endif

  call i2mex_toAxisFraction(1, ns, (/0._i2mex_r8/), psi, 0.0_i2mex_r8, &
       & i2mex_to_axis_fraction, qpp, iok)
  call i2mex_toEdgeFraction(1, ns, (/0._i2mex_r8/), psi,  &
       & i2mex_to_edge_fraction, qpp, iok)

end subroutine i2mex_getQpp


subroutine i2mex_getQFromPsiGXZ(nt1, ns, psi, g, x, z, q, ier)

  ! Compute the safety factor profile from g, X and Z and return the result
  ! in q. This routine it typically used when the q profile information needed
  ! to build i2mex is not available from the equilibrium. *NOTE* This routine
  ! does NOT require i2mex to be initialized. Also be warned that the present
  ! implementation suffers from some inaccuracy near axis.

  use i2mex_mod
  use ezspline_obj
  use ezspline
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: g(ns) ! x^2 B.grad phi
  real(i2mex_r8), intent(in) :: psi(ns) ! pol flux/2*pi
  real(i2mex_r8), intent(in) :: x(nt1, ns), z(nt1, ns) ! the grid
  real(i2mex_r8), intent(out) :: q(ns)
  integer, intent(out) :: ier
  

  integer i, iok
  real(r8), dimension(nt1) :: the
  real(r8), dimension(:,:), allocatable :: integrand, xt, xp, zt, zp
  type(ezspline2) :: xspl, zspl

   include 'cubinterp.h'

   ier = 0
   if(ns <= 5) then
      ier = 124
      return
   endif
   if(nt1 <= 5) then
      ier = 125
      return
   endif


   the = i2mex_twopi_r8*(/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)

   allocate(xt(nt1, ns), xp(nt1, ns))
   allocate(zt(nt1, ns), zp(nt1, ns))
   allocate(integrand(nt1, ns))

   call ezspline_init(xspl, nt1, ns, (/-1, -1/), (/0, 0/), iok)
   call ezspline_error(iok)
   xspl%x1 = the
   xspl%x2 = psi
   call ezspline_setup(xspl, x, iok)
   call ezspline_error(iok)
   call EZspline_derivative(xspl, 1, 0, nt1, ns, the, psi, xt, iok)
   call ezspline_error(iok)
   call EZspline_derivative(xspl, 0, 1, nt1, ns, the, psi, xp, iok)
   call ezspline_error(iok)
   call ezspline_free(xspl, iok)
   call ezspline_error(iok)
   
   call ezspline_init(zspl, nt1, ns, (/-1, -1/), (/0, 0/), iok)
   call ezspline_error(iok)
   zspl%x1 = the
   zspl%x2 = psi
   call ezspline_setup(zspl, z, iok)
   call ezspline_error(iok)
   call EZspline_derivative(zspl, 1, 0, nt1, ns, the, psi, zt, iok)
   call ezspline_error(iok)
   call EZspline_derivative(zspl, 0, 1, nt1, ns, the, psi, zp, iok)
   call ezspline_error(iok)
   call ezspline_free(zspl, iok)
   call ezspline_error(iok)

   ! jacobian/x**2
   integrand = abs(xp*zt - xt*zp)/x

   ! extrapolate to axis
    do i=1, nt1
       integrand(i,1) = FCCCC0( &
            &  integrand(i,2), integrand(i,3), integrand(i,4), integrand(i,5), &
            &                 psi(2),  psi(3),  psi(4), psi(5), &
            & psi(1))
    enddo
   
    call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, q, iok) 
    q = q * g/i2mex_twopi_r8

   ! extrapolate to axis
!!$   q(1) = FCCCC0( &
!!$        & q(2), q(3), q(4), q(5), &
!!$        &  psi(2),  psi(3),  psi(4), psi(5), &
!!$        &  psi(1))
   q(1) = 2._r8*q(2) - q(3)

   
   deallocate(xt, xp)
   deallocate(zt, zp)
   deallocate(integrand)
   
 end subroutine i2mex_getQFromPsiGXZ
