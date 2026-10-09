subroutine i2mex_getGlasserGreenJohnsonEFH(ns, psi, e, f, h, ier)

  ! Compute the Glasser, Greene and Johnson coefficients
  ! E, F, and H on the set of poloidal surface psi

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: e(ns), f(ns), h(ns)
  integer, intent(out) :: ier

  integer i, nt1, iok
  real(i2mex_r8) :: zvp(ns), zvpp(ns)
  real(i2mex_r8) :: zb2(ns), zb2og2(ns), zob2g2(ns), &
       & zog2(ns), zob2(ns), zg2ob2(ns)
  real(i2mex_r8) :: qp(ns), pp(ns), g(ns)
  real(i2mex_r8), dimension(:), allocatable :: the
  real(i2mex_r8), dimension(:,:), allocatable :: integrand, jac, x, b2
  real(i2mex_r8), dimension(:,:), allocatable :: gtt, gtp, gpp

  ier = 0

  call i2mex_getOriNt1(nt1, iok)
  call i2mex_error(iok)
  ALLOCATE(the(nt1))
  call i2mex_getOriT(nt1, the, iok)
  call i2mex_error(iok)
  ALLOCATE(integrand(nt1, ns), jac(nt1, ns), x(nt1, ns), b2(nt1, ns))
  ALLOCATE(gtt(nt1, ns), gtp(nt1, ns), gpp(nt1, ns))

  call i2mex_getPp(ns, psi, pp, iok)
  call i2mex_error(iok)
  call i2mex_getQp(ns, psi, qp, iok)
  call i2mex_error(iok)
  call i2mex_getG(ns, psi, g, iok)
  call i2mex_error(iok)
  

  call i2mex_GetJ(nt1, ns, the, psi, jac, iok)
  call i2mex_error(iok)
  call i2mex_getPoloidalIntegral(nt1, ns, the, jac, zvp, iok)
  call i2mex_error(iok)

  call i2mex_GetJp(nt1, ns, the, psi, integrand, iok)
  call i2mex_error(iok)
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, zvpp, iok)
  call i2mex_error(iok)
  zvpp = zvpp/i2mex_twopi_r8

  call i2mex_getMetric(nt1, ns, the, psi, gtt, gtp, gpp, iok)
  call i2mex_error(iok)
  call i2mex_getX(nt1, ns, the, psi, x, iok)
  call i2mex_error(iok)

  b2 = (gpp + spread(g**2, dim=1, ncopies=nt1))/x**2
  
  integrand = jac*b2
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, zb2, iok)
  call i2mex_error(iok)

  integrand = jac*b2/gpp
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, zb2og2, iok)

  integrand = jac/(b2*gpp)
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, zob2g2, iok)

  integrand = jac/gpp
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, zog2, iok)

  integrand = jac/b2
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, zob2, iok)

  integrand = jac*gpp/b2
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, zg2ob2, iok)
  

  e = -pp*zb2og2*(zvpp - g*qp*zvp/zb2)/(i2mex_twopi_r8*qp**2)
  f = pp**2 *( g**2 *(zb2og2*zob2g2 - zog2**2) &
       & + zb2og2*zob2 )/(i2mex_twopi_r8**2 * qp**2)
  h = g*pp*(zog2 - zvp*zb2og2/zb2)/(i2mex_twopi_r8*qp)

  DEALLOCATE(the)
  
  DEALLOCATE(integrand, jac, x, b2)
  DEALLOCATE(gtt, gtp, gpp)


  if(iok /=0) ier = 106
  

end subroutine i2mex_getGlasserGreenJohnsonEFH
