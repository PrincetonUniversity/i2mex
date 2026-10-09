subroutine i2mex_getSpest1(ns, psi, s, ier)
 
  ! Compute
  ! (integral_0^psi d psi q/g)/(integral_0^psimax d psi q/g)
  ! This is the radial coordinate used in PEST1 according to
  ! Manickam.
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: s(ns)
  integer, intent(out) :: ier
 
  real(i2mex_r8) :: f(ns), g(ns)
  integer iok
 
  ier = 0
  call i2mex_getQ(ns, psi, f, iok)
  call i2mex_getG(ns, psi, g, iok)
  f = f /g
  call i2mex_integrate1d(ns, psi, f, ns, psi, s)
  ! normalize
  s = s/s(ns)
 
  if(iok/=0) ier = 146
 
 
end subroutine i2mex_getSpest1
 
subroutine i2mex_getUpsilon(res, ier)
 
  ! Compute (R0/2*pi) * the volume integral dV/X**2 and return result
  ! in 'res'.
 
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8), dimension(:), allocatable :: the, psi, circ
  real(i2mex_r8), dimension(:,:), allocatable :: integrand, x
  real(i2mex_r8) :: R0
 
  integer nt1, nt, ns
  real(i2mex_r8), dimension(:), allocatable :: result
 
  ier = 0
 
  call i2mex_getOriNt1(nt1, iok)
  call i2mex_getOriNs(ns, iok)
  ALLOCATE(the(nt1))
  ALLOCATE(psi(ns), circ(ns))
  ALLOCATE(integrand(nt1, ns), x(nt1, ns))
  call i2mex_getOriT(nt1, the, iok)
  call i2mex_getOriPsi(ns, psi, iok)
 
 
  call i2mex_getJ(nt1, ns, the, psi, integrand, iok)
  call i2mex_getX(nt1, ns, the, psi, x, ier)
  integrand = integrand/x**2
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, circ, iok)
 
  ! integrate over psi
  res = 0.0_i2mex_r8
  do i=1, ns-1
     res = res + &
          & 0.5_i2mex_r8*( circ(i+1) + circ(i  ) )*(psi(i+1) - psi(i  ))
  enddo
 
  call i2mex_getRmajor(R0, iok)
  res = R0 * abs(res)
 
  DEALLOCATE(the)
  DEALLOCATE(psi, circ)
  DEALLOCATE(integrand, x)
 
  if(iok /= 0) ier = 145
 
end subroutine i2mex_getUpsilon
 
subroutine i2mex_getSurface(res, ier)
 
  ! Compute the plasma surface and return result in 'res'.
 
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8), dimension(:), allocatable :: the, psi, circ
  real(i2mex_r8), dimension(:,:), allocatable :: integrand, x
 
  integer nt1, nt, ns
  real(i2mex_r8), dimension(:), allocatable :: result
 
  ier = 0
 
!!$    ! call eq_flxint_chinit(0,nchbdys,chi_bdys,ierr)
!!$    call fluxav_nzones_get(nt, ns)
!!$    allocate(result(ns))
!!$    call eq_flxint('DAREA', 1, result, 1, ns, iok)
!!$    res = sum(result)
!!$    deallocate(result)
 
  call i2mex_getOriNt1(nt1, iok)
  call i2mex_getOriNs(ns, iok)
  ALLOCATE(the(nt1))
  ALLOCATE(psi(ns), circ(ns))
  ALLOCATE(integrand(nt1, ns), x(nt1, ns))
  call i2mex_getOriT(nt1, the, iok)
  call i2mex_getOriPsi(ns, psi, iok)
 
 
  call i2mex_getJ(nt1, ns, the, psi, integrand, iok)
  call i2mex_getX(nt1, ns, the, psi, x, ier)
  integrand = integrand/x
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, circ, iok)
 
  res = 0.0_i2mex_r8
  do i=1, ns-1
     res = res + &
          & 0.5_i2mex_r8*( circ(i+1) + circ(i  ) )*(psi(i+1) - psi(i  ))
  enddo
  res = abs(res)
 
  DEALLOCATE(the)
  DEALLOCATE(psi, circ)
  DEALLOCATE(integrand, x)
 
  if(iok /= 0) ier = 76
 
end subroutine i2mex_getSurface
 
subroutine i2mex_getVolume(res, ier)
 
  ! Compute the plasma volume and return result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer nt, ns, iok
  real(i2mex_r8), dimension(:), allocatable :: v, psi
  real(i2mex_r8), dimension(:), allocatable :: result
 
  ier = 0
 
!!$    ! call eq_flxint_chinit(0,nchbdys,chi_bdys,ierr)
!!$    call fluxav_nzones_get(nt, ns)
!!$    allocate(result(ns))
!!$    call eq_flxint('DVOL', 1, result, 1, ns, iok)
!!$    res = sum(result)
!!$    deallocate(result)
 
  call i2mex_getOriNs(ns, iok)
  call i2mex_error(iok)
 
  ALLOCATE(psi(ns), v(ns))
 
  call i2mex_getOriPsi(ns, psi, iok)
  call i2mex_getV(ns, psi, v, iok)
 
  res = v(ns)
 
  DEALLOCATE(psi, v)
 
  if(iok /= 0) ier = 77
 
end subroutine i2mex_getVolume
 
subroutine i2mex_getV(ns, psi, res, ier)
 
  ! Compute the enclosed volume V(psi) and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier
 
  integer iok
  real(i2mex_r8) :: vp(ns)
 
  ier = 0
 
  call i2mex_getVp(ns, psi, vp, iok)
  call i2mex_error(iok)
  call i2mex_integrate1D(ns, psi, vp, ns, psi, res)
 
  if(iok /= 0) ier = 107
 
end subroutine i2mex_getV
 
 
subroutine i2mex_getVp(ns, psi, res, ier)
 
  ! Compute d/dpsi of enclosed volume V and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier
 
  integer i, iok, nt1
  real(i2mex_r8) :: circ(ns)
  real(i2mex_r8), dimension(:,:), allocatable :: integrand
  real(i2mex_r8), dimension(:), allocatable :: the
 
  ier = 0
 
  call i2mex_getOriNt1(nt1, iok)
  call i2mex_error(iok)
  ALLOCATE(the(nt1))
  ALLOCATE(integrand(nt1, ns))
  the = i2mex_twopi_r8*(/ (real(i-1,i2mex_r8)/real(nt1-1,i2mex_r8), i=1, nt1) /)
 
  call i2mex_getJ(nt1, ns, the, psi, integrand, iok)
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, res, iok)
 
  res = abs(res)*i2mex_twopi_r8
 
  if(iok /= 0) ier = 108
 
  DEALLOCATE(the)
  DEALLOCATE(integrand)
 
end subroutine i2mex_getVp
 
subroutine i2mex_getVolumeAveragedPressure(nt1, ns, the, psi, res, ier)
 
  ! compute <p>, the volume average pressure and return the result
  ! in 'res'. ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8) :: volume, ds
  real(i2mex_r8) :: integrand(nt1, ns), circ(ns), pressure(ns)
 
  ier = 0
  call i2mex_getP(ns, psi, pressure, iok)
  call i2mex_getJ(nt1, ns, the, psi, integrand, iok)
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, circ, iok)
 
  volume = 0.0_i2mex_r8
  res    = 0.0_i2mex_r8
  do i=1, ns-1
     ds = 0.5_i2mex_r8*( circ(i+1) + circ(i  ) )*( psi(i+1) - psi(i  ))
     volume = volume + ds
     res    = res    + ds * 0.5_i2mex_r8*( pressure(i+1) + pressure(i  ) )
  enddo
  res = res / volume
 
  if(iok /= 0) ier = 78
 
end subroutine i2mex_getVolumeAveragedPressure
 
 
subroutine i2mex_getVolumeAveragedBSquare(nt1, ns, the, psi, res, ier)
 
  ! compute <B^2>, the volume average of the square of the magnetic field
  ! and return the result in 'res'. ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8) :: volume, ds1, ds2
  real(i2mex_r8), dimension(:,:), allocatable :: xt,xp,zt,zp,x,jac, integrand
  real(i2mex_r8) :: circ1(ns), circ2(ns), g(ns)
 
  ier = 0
 
  allocate(xt(nt1, ns), xp(nt1, ns), &
       & zt(nt1, ns), zp(nt1, ns), &
       & x(nt1, ns), integrand(nt1, ns), jac(nt1, ns), stat=iok)
 
  call i2mex_getG(ns, psi, g, iok)
  call i2mex_getX(nt1, ns, the, psi, x, iok)
  call i2mex_getGradX(nt1, ns, the, psi, xt, xp, iok)
  call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, iok)
 
  jac = abs(xp*zt - xt*zp)*x
  call i2mex_toAxis(nt1, ns, the, psi, +0.0_i2mex_r8, jac, iok)
 
!!$  call i2mex_getJ(nt1, ns, the, psi, integrand, iok)
!!$ 
!!$  integrand = jac * ( &
!!$       & (xt**2 + zt**2)/(xp*zt-xt*zp)**2 + &
!!$       & spread(g**2, dim=1, ncopies=nt1) &
!!$       & ) / x**2
  integrand = (xt**2 + zt**2)/jac + &
       & jac*spread(g**2, dim=1, ncopies=nt1)/ x**2
  call i2mex_toAxis(nt1, ns, the, psi, +0.0_i2mex_r8, integrand, iok)
 
  call i2mex_getPoloidalIntegral(nt1, ns, the,       jac, circ1, iok)
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, circ2, iok)
 
  volume = 0.0_i2mex_r8
  res    = 0.0_i2mex_r8
  do i=1, ns-1
     ds1 = 0.5_i2mex_r8*( circ1(i+1) + circ1(i  ) )*( psi(i+1) - psi(i  ))
     ds2 = 0.5_i2mex_r8*( circ2(i+1) + circ2(i  ) )*( psi(i+1) - psi(i  ))
     volume = volume + ds1
     res    = res    + ds2
  enddo
  res = res/volume
 
  deallocate(xt, xp, &
       & zt, zp, &
       & x, integrand, jac, stat=iok)
 
  if(iok /= 0) ier = 79
 
 
end subroutine i2mex_getVolumeAveragedBSquare
 
subroutine i2mex_getBeta(nt1, ns, the, psi, res, ier)
 
  ! Compute (the volume averaged) beta and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer iok
  real(i2mex_r8) :: b0
 
  ier = 0
  call i2mex_getVolumeAveragedPressure(nt1, ns, the, psi, res, iok)
  call i2mex_getB0(b0, iok)
  res = 2.0_i2mex_r8 * res/b0**2
 
  if(iok /= 0) ier = 80
 
end subroutine i2mex_getBeta
 
subroutine i2mex_getBetaStar(nt1, ns, the, psi, res, ier)
 
  ! Compute the (volume averaged) 2 sqrt(<p^2>)/B0^2
  ! and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer iok, i
  real(i2mex_r8) :: b0, volume, pressure(ns), circ(ns), ds
  real(i2mex_r8), dimension(:,:), allocatable :: integrand
 
  ier = 0
 
  allocate(integrand(nt1, ns))
 
  call i2mex_getP(ns, psi, pressure, iok)
  call i2mex_getJ(nt1, ns, the, psi, integrand, iok)
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, circ, iok)
 
  volume = 0.0_i2mex_r8
  res    = 0.0_i2mex_r8
  do i=1, ns-1
     ds = 0.5_i2mex_r8*( circ(i+1) + circ(i  ) )*( psi(i+1) - psi(i  ))
     volume = volume + ds
     res    = res    + ds * 0.5_i2mex_r8* &
          & ( pressure(i+1) + pressure(i  ) )**2
  enddo
  call i2mex_getB0(b0, iok)
  res = 2.0_i2mex_r8 * sqrt(res/volume)/b0**2
 
  deallocate(integrand)
 
  if(iok /= 0) ier = 80
 
end subroutine i2mex_getBetaStar
 
subroutine i2mex_getBetaToroidal(nt1, ns, the, psi, res, ier)
 
  ! Compute beta-toroidal and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer iok
  real(i2mex_r8) :: b2_ave, B0
 
  ier = 0
  call i2mex_getVolumeAveragedPressure(nt1, ns, the, psi, res, iok)
  call i2mex_getB0(B0, iok)
  res = 2.0_i2mex_r8 * res / B0**2
 
  if(iok /= 0) ier = 81
 
end subroutine i2mex_getBetaToroidal
 
subroutine i2mex_getBetaPoloidalFreidberg(nt1, ns, the, psi, res, ier)
 
  ! Compute the approximate beta-poroidal and return the result in 'res'. The
  ! formula used is loosely based on Eq(4.32) in Ideal Magnetohydrodynamics
  ! by J P Freidberg, since there appears to be an error in Eq(4.29) therein.
  ! Call i2mex_getBetaPoloidal for a more accurate estimate of beta-poloidal.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer iok
  real(i2mex_r8) :: Ip, B0, area, a
 
  ier = 0
  call i2mex_getVolumeAveragedPressure(nt1, ns, the, psi, res, iok)
  call i2mex_getPlasmaCurrent(nt1, ns, the, psi, Ip, iok)
  call i2mex_getSurface(area, iok)
  call i2mex_getRminor(a, iok)
 
  ! take the poloidal field to be Bp = Ip/(2*pi*a*sqrt(kappa)).
 
  res = 4.0_i2mex_r8 * i2mex_twopi_r8 * area * res / Ip**2
 
  if(iok /= 0) ier = 82
 
end subroutine i2mex_getBetaPoloidalFreidberg
 
subroutine i2mex_getBetaPoloidal(nt1, ns, the, psi, res, ier)
 
  ! Compute beta-toroidal and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer iok, n2 ,i
  real(i2mex_r8) :: circ, Ip, Bp, x(nt1, ns), z(nt1,ns)
  real(i2mex_r8), dimension(i2mex_o%nt1) :: xbound, zbound, t
 
  ier = 0
 
  ! Compute the plasma circumference (integral dl) on the last surface (n2)
 
  circ = 0.0_i2mex_r8
  n2 = i2mex_o%ns
  call i2mex_getOriT(i2mex_o%nt1, t, iok) ! use original grid
  call i2mex_error(iok)
  call i2mex_getX(i2mex_o%nt1, 1, t, psi(ns), xbound, iok)
  call i2mex_error(iok)
  call i2mex_getZ(i2mex_o%nt1, 1, t, psi(ns), zbound, iok)
  call i2mex_error(iok)
  do i = 1, i2mex_o%nt1-1
     circ = circ + sqrt( &
          & ( xbound(i+1) - xbound(i  ) )**2 + &
          & ( zbound(i+1) - zbound(i  ) )**2 &
          &            )
  enddo
 
  ! B-poloidal
 
  call i2mex_getPlasmaCurrent(nt1, ns, the, psi, Ip, iok)
  Bp = Ip/circ
 
  call i2mex_getVolumeAveragedPressure(nt1, ns, the, psi, res, iok)
  res = 2.0_i2mex_r8 * res/ Bp**2
 
  if(iok /= 0) ier = 83
 
end subroutine i2mex_getBetaPoloidal
 
subroutine i2mex_getBetaN(nt1, ns, the, psi, res, ier)
 
  ! Compute the normalized beta (beta_N) and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer iok
  real(i2mex_r8) :: Ip, B0, a
 
  ier = 0
  call i2mex_getBeta(nt1, ns, the, psi, res, iok)
  call i2mex_getB0(B0, iok)
  call i2mex_getPlasmaCurrent(nt1, ns, the, psi, Ip, iok)
  call i2mex_getRminor(a, iok)
 
  ! 100 <beta> a B0/ Ip where Ip is in MA
 
  res = 100.0_i2mex_r8 * res * a * B0 * 0.2_i2mex_r8*i2mex_twopi_r8 / abs(Ip)
 
  if(iok /= 0) ier = 84
 
end subroutine i2mex_getBetaN
 
subroutine i2mex_getFrame(xin, xout, zbot, ztop, ier)
 
  ! Compute the plasma boundary frame; ie the abscissae xin and xout
  ! where the plasma boundary intersects the mid plane, and the ordinates
  ! zbot and ztop where the plasma boundary intersects the vertical plane
  ! passing through the magnetic axis.
 
  use ezspline
  use ezspline_obj
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: xin  ! mid plane intersect inside
  real(i2mex_r8), intent(out) :: xout ! mid plane intersect outside
  real(i2mex_r8), intent(out) :: zbot ! bottom intersect below magnetic axis
  real(i2mex_r8), intent(out) :: ztop ! top intersect above magnetic axis
  integer, intent(out) :: ier
 
  real(i2mex_r8) :: psi(i2mex_o%ns)
  real(i2mex_r8) :: t(i2mex_o%nt1), thetacyl(i2mex_o%nt1), alpha(i2mex_o%nt1)
  real(i2mex_r8) :: xbound(i2mex_o%nt1), zbound(i2mex_o%nt1)
  real(i2mex_r8) :: tcyl0, sgn = 1.0_i2mex_r8, bigdt, dt
  real(i2mex_r8) :: aout, ain, atop, abot, xm, zm
  integer i, iok
 
  type(ezspline1) :: xspl, zspl
 
  ier = 0
  call i2mex_getOriT(i2mex_o%nt1, t, iok) ! use original grid
  call i2mex_error(iok)
  call i2mex_getOriPsi(i2mex_o%ns, psi, iok) ! use original grid
  call i2mex_error(iok)
  call i2mex_getX(i2mex_o%nt1, 1, t, psi(i2mex_o%ns), xbound, iok)
  call i2mex_error(iok)
  call i2mex_getZ(i2mex_o%nt1, 1, t, psi(i2mex_o%ns), zbound, iok)
  call i2mex_error(iok)
  call i2mex_getRmagnetic(xm, iok)
  call i2mex_error(iok)
  call i2mex_getZmagnetic(zm, iok)
  call i2mex_error(iok)
 
  ! Cylindrical angle.
 
  thetacyl = atan2(zbound-zm, xbound-xm)
  thetacyl(i2mex_o%nt1) = thetacyl(1) + i2mex_twopi_r8
 
  ! alpha ~ cylindrical angle but monotonically increasing with xbound
 
  bigdt = i2mex_pi_r8
  if(thetacyl(2) - thetacyl(1) < 0.0_i2mex_r8) sgn = -1._i2mex_r8
 
  alpha(1) = 0.0_i2mex_r8
  do i=2, i2mex_o%nt1 - 1
     dt = abs(thetacyl(i) - thetacyl(i-1))
     if(dt > bigdt) then
!!$          print *,'segment ',i,i-1,' dt = ', dt
        dt = abs(dt - i2mex_twopi_r8)
     endif
!!$       print *,'i=', i, ' dt = ', dt
     alpha(i) = alpha(i-1) + dt
  enddo
  alpha(i2mex_o%nt1) = alpha(1) + i2mex_twopi_r8
 
!!$    ! debug
!!$    do i = 2, i2mex_o%nt1
!!$       if (alpha(i) < alpha(i-1)) then
!!$          print*,'***error** i=',i, alpha(i-1),alpha(i)
!!$          stop 'error'
!!$       endif
!!$    enddo
!!$    ! debug
 
!!$  aout = modulo(                                       thetacyl(1), &
!!$       & i2mex_twopi_r8)
!!$  ain  = modulo(                        i2mex_pi_r8  + thetacyl(1), &
!!$       & i2mex_twopi_r8)
!!$  atop = modulo(             i2mex_pi_r8/2._i2mex_r8 + thetacyl(1), &
!!$       & i2mex_twopi_r8)
!!$  abot = modulo(3._i2mex_r8* i2mex_pi_r8/2._i2mex_r8 + thetacyl(1), &
!!$       & i2mex_twopi_r8)

  aout = 0._i2mex_r8
  ain  = i2mex_pi_r8
  atop  = i2mex_pi_r8/2._i2mex_r8
  abot  = 3._i2mex_r8 * i2mex_pi_r8/2._i2mex_r8
 
  call ezspline_init(xspl, i2mex_o%nt1, (/-1, -1/), iok)
  xspl%x1 = alpha
  call ezspline_setup(xspl, xbound, iok)
  call ezspline_interp(xspl, aout, xout, iok)
  call ezspline_error(iok)
  call ezspline_interp(xspl, ain , xin , iok)
  call ezspline_error(iok)
  call ezspline_free(xspl, iok)
  if(iok/=0) ier = 140
 
  call ezspline_init(zspl, i2mex_o%nt1, (/-1, -1/), iok)
  zspl%x1 = alpha
  call ezspline_setup(zspl, zbound, iok)
  call ezspline_interp(zspl, abot, zbot, iok)
  call ezspline_error(iok)
  call ezspline_interp(zspl, atop, ztop, iok)
  call ezspline_error(iok)
  call ezspline_free(zspl, iok)
  if(iok/=0) ier = 141
 
!!$    print *,'ain, aout, atop, abot=', ain, aout, atop, abot
!!$    print *,'xin, xout, ztop, zbot=', xin, xout, ztop, zbot
 
end subroutine i2mex_getFrame
 
subroutine i2mex_getRminor(res, ier)
 
  ! Compute the minor radius, and return the result in 'res'. The
  ! minor radius is defined as the half distance between Xout and Xin,
  ! where Xin and Xout are the plasma boundary points intersecting the
  ! mid plane.
  !
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  real(i2mex_r8) :: xout, xin, zbot, ztop
 
  ier = 0
  call i2mex_getFrame(xin, xout, zbot, ztop, ier)
  call i2mex_error(ier)
  res = (xout-xin)/2._i2mex_r8
 
end subroutine i2mex_getRminor
 
subroutine i2mex_getRmajor(res, ier)
 
  ! Compute the major radius, and return the result in 'res'. The
  ! major radius is defined as the average of Xout and Xin, where
  ! Xin and Xout are the plasma boundary points intersecting the
  ! mid plane.
  !
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  real(i2mex_r8) :: xout, xin, zbot, ztop
 
  ier = 0
  call i2mex_getFrame(xin, xout, zbot, ztop, ier)
  call i2mex_error(ier)
  res = (xout+xin)/2._i2mex_r8
 
end subroutine i2mex_getRmajor
 
subroutine i2mex_getRmagnetic(res, ier)
 
  ! Compute the magnetic radius position, and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  real(i2mex_r8) :: t(i2mex_o%nt1), xmag(i2mex_o%nt1)
  integer iok
 
  ier = 0
  call i2mex_getOriT(i2mex_o%nt1, t, iok) ! use original grid
  call i2mex_getX(i2mex_o%nt1, 1, t, (/i2mex_o%psi_axis/), xmag, iok)
  call i2mex_error(iok)
 
  res = sum(xmag)/real(i2mex_o%nt1,i2mex_r8)
 
end subroutine i2mex_getRmagnetic
 
subroutine i2mex_getZmagnetic(res, ier)
 
  ! Compute the magnetic Z position, and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  real(i2mex_r8) :: t(i2mex_o%nt1), zmag(i2mex_o%nt1)
  integer iok
 
  ier = 0
  call i2mex_getOriT(i2mex_o%nt1, t, iok) ! use original grid
  call i2mex_getZ(i2mex_o%nt1, 1, t, (/i2mex_o%psi_axis/), zmag, iok)
  call i2mex_error(iok)
 
  res = sum(zmag)/real(i2mex_o%nt1,i2mex_r8)
 
end subroutine i2mex_getZmagnetic
 
subroutine i2mex_getB0(res, ier)
 
  ! Compute the vacuum magnetic field at the major radius position,
  ! and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer iok
  real(i2mex_r8) :: R0
 
 
  ier = 0
  call i2mex_getRmajor(R0, iok)
  call i2mex_getG(1, (/i2mex_o%psi_edge/), res, iok)
  res = res / R0
 
end subroutine i2mex_getB0
 
subroutine i2mex_getJdotBOverBSquare(nt1, ns, the, psi, res, ier)
 
  ! Compute <J.B>/<B^2> and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8) :: g(ns), pp(ns), b2(ns), jave(ns)
  real(i2mex_r8), dimension(:,:), allocatable :: integrand, dum1, dum2
 
  call i2mex_getG(ns, psi, g, iok)
  call i2mex_getPp(ns, psi, pp, iok)
  call i2mex_getGp(ns, psi, res, iok) ! res is g'
  allocate(integrand(nt1, ns))
  allocate(dum1(nt1, ns))
  allocate(dum2(nt1, ns))
  ! |grad psi|^2
  call i2mex_getMetric(nt1, ns, the, psi, dum1, dum2, integrand, ier)
  call i2mex_getX(nt1, ns, the, psi, dum1, ier) ! dum1 is X
  ! B^2
  integrand = (integrand + spread(g*g, dim=1, ncopies=nt1))/dum1**2
  ! Jacobian
  call i2mex_getJ(nt1, ns, the, psi, dum2, ier)
  ! J B^2
  integrand = dum2*integrand
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, b2, ier)
  call i2mex_getPoloidalIntegral(nt1, ns, the, dum2, jave, ier)
  !
  res = - res - g*pp*(jave/b2)
 
  deallocate(integrand)
  deallocate(dum1)
  deallocate(dum2)
 
  if(iok/=0) ier = 109
 
end subroutine i2mex_getJdotBOverBSquare
 
subroutine i2mex_getPolAvrgOneOverRSquare(nt1, ns, the, psi, res, ier)
 
  ! Compute <1/R^2> = integral d theta J /X^2 / integral d theta J 
  ! and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8), dimension(:,:), allocatable :: integrand, jac
  real(i2mex_r8) :: anorm(ns)

  iok = 0
  allocate(integrand(nt1, ns), jac(nt1, ns))

  call i2mex_getX(nt1, ns, the, psi, integrand, ier) ! X
  ! Jacobian
  call i2mex_getJ(nt1, ns, the, psi, jac, ier)
  ! J/ X^2
  integrand = jac / integrand**2
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, res, ier)
  ! normalize
  call i2mex_getPoloidalIntegral(nt1, ns, the, jac, anorm, ier)
  res = res / anorm
  !
  deallocate(integrand)
  deallocate(jac)
 
  if(iok/=0) ier = 159
 
end subroutine i2mex_getPolAvrgOneOverRSquare

subroutine i2mex_getPolAvrgBSquare(nt1, ns, the, psi, res, ier)
 
  ! Compute <B^2> = integral d theta J B^2 / integral d theta J 
  ! and return the result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8) :: g(ns)
  real(i2mex_r8), dimension(:,:), allocatable :: integrand, dum1, dum2
  real(i2mex_r8) :: anorm(ns)

  iok = 0
  allocate(integrand(nt1, ns))
  allocate(dum1(nt1, ns), dum2(nt1, ns))

  call i2mex_getG(ns, psi, g, ier)

  ! |grad psi|^2
  call i2mex_getMetric(nt1, ns, the, psi, dum1, dum2, integrand, ier)
  integrand = integrand + spread(g*g, dim=1, ncopies=nt1)
  ! X
  call i2mex_getX(nt1, ns, the, psi, dum1, ier) ! dum1 = X
  ! Jacobian
  call i2mex_getJ(nt1, ns, the, psi, dum2, ier) ! dum2 = jac
  call i2mex_getPoloidalIntegral(nt1, ns, the, dum2, anorm, ier)
  call i2mex_error(ier)
  ! J/ X^2
  integrand = dum2 * integrand / dum1**2
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, res, ier)
  call i2mex_error(ier)
  res = res / anorm
  !
  deallocate(integrand)
  deallocate(dum1)
  deallocate(dum2)
 
  if(iok/=0) ier = 160
 
end subroutine i2mex_getPolAvrgBSquare

subroutine i2mex_getPlasmaCurrent(nt1, ns, the, psi, res, ier)
 
  ! Compute Ip, the plasma current, and return the result in 'res'.
  ! The current is in [T/m] units; you must divide by 0.4*pi to 
  ! get the current in [MA].
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8) :: ds
 
  real(i2mex_r8), dimension(:,:), allocatable :: integrand, x
  real(i2mex_r8), dimension(:), allocatable :: circ, pp, g, gp
 
!!$    real(i2mex_r8) :: integrand(nt1, ns), x(nt1, ns), &
!!$         & circ(ns), pp(ns), g(ns), gp(ns)
 
  ier = 0
 
  allocate(integrand(nt1, ns), x(nt1, ns))
  allocate(circ(ns), pp(ns), g(ns), gp(ns))
 
  call i2mex_getPP(ns, psi, pp, iok)
  call i2mex_getG(ns, psi, g, iok)
  call i2mex_getGP(ns, psi, gp, iok)
  call i2mex_getJ(nt1, ns, the, psi, integrand, iok)
  call i2mex_getX(nt1, ns, the, psi, x, iok)
  integrand = -integrand*( spread(pp, dim=1, ncopies=nt1) + &
       & spread(g*gp, dim=1, ncopies=nt1)/x**2 )
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, circ, iok)
  res    = 0.0_i2mex_r8
  do i=1, ns-1
     ds = 0.5_i2mex_r8*( circ(i+1) + circ(i  ) )*( psi(i+1) - psi(i  ))
     res    = res    + ds
  enddo
 
  if(iok /= 0) ier = 85
 
  deallocate(integrand, x)
  deallocate(circ, pp, g, gp)
 
end subroutine i2mex_getPlasmaCurrent

subroutine i2mex_getPlasmaCurProf(nt1, ns, the, psi, res, ier)
 
  ! Compute Ip(psi), the enclosed plasma current profile, and return the 
  ! result in 'res'.  The code here is the same as i2mex_getPlasmaCurrent,
  ! except that the entire profile is returned, not just the value
  ! integrated out to the boundary (dmc Jun 12 2003). The current is
  ! in [T/m] units, you must divide the current by 0.4*pi to convert to
  ! [MA].
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8) :: ds
 
  real(i2mex_r8), dimension(:,:), allocatable :: integrand, x
  real(i2mex_r8), dimension(:), allocatable :: circ, pp, g, gp
 
!!$    real(i2mex_r8) :: integrand(nt1, ns), x(nt1, ns), &
!!$         & circ(ns), pp(ns), g(ns), gp(ns)
 
  ier = 0
 
  allocate(integrand(nt1, ns), x(nt1, ns))
  allocate(circ(ns), pp(ns), g(ns), gp(ns))
 
  call i2mex_getPP(ns, psi, pp, iok)
  call i2mex_getG(ns, psi, g, iok)
  call i2mex_getGP(ns, psi, gp, iok)
  call i2mex_getJ(nt1, ns, the, psi, integrand, iok)
  call i2mex_getX(nt1, ns, the, psi, x, iok)
  integrand = -integrand*( spread(pp, dim=1, ncopies=nt1) + &
       & spread(g*gp, dim=1, ncopies=nt1)/x**2 )
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, circ, iok)
  res(1)  = 0.0_i2mex_r8
  do i=1, ns-1
     ds = 0.5_i2mex_r8*( circ(i+1) + circ(i  ) )*( psi(i+1) - psi(i  ))
     res(i+1)    = res(i)    + ds
  enddo
 
  if(iok /= 0) ier = 85
 
  deallocate(integrand, x)
  deallocate(circ, pp, g, gp)
 
end subroutine i2mex_getPlasmaCurProf
 
subroutine i2mex_getLi(nt1, ns, the, psi, res, ier)
 
  ! Compute li, the internal plasma inductance, and return the
  ! result in 'res'.
  ! ier error flag (ok if =0)
 
  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1)
  real(i2mex_r8), intent(in) :: psi(ns)
  real(i2mex_r8), intent(out) :: res
  integer, intent(out) :: ier
 
  integer i, iok
  real(i2mex_r8) :: ds2, circ2(ns), Ip, R0
  real(i2mex_r8), dimension(:,:), allocatable :: x, xt, xp, zt, zp, integrand
 
  ier = 0
 
  allocate(x(nt1, ns), xt(nt1, ns), xp(nt1, ns), &
       & zt(nt1, ns), zp(nt1, ns), &
       & integrand(nt1, ns), stat=iok)
 
  call i2mex_getX(nt1, ns, the, psi, x, iok)
  call i2mex_getGradX(nt1, ns, the, psi, xt, xp, iok)
  call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, iok)
  integrand = -(xp*zt - xt*zp)/(x*(xp**2 + zp**2))
  call i2mex_getPoloidalIntegral(nt1, ns, the, integrand, circ2, iok)
 
  res    = 0.0_i2mex_r8
  do i=1, ns-1
     ds2 = 0.5_i2mex_r8*( circ2(i+1) + circ2(i  ) )*( psi(i+1) - psi(i  ))
     res    = res    + ds2
  enddo
 
  ! res = integral d Volume B-poloidal^2 /(2*pi)
 
  call i2mex_getPlasmaCurrent(nt1, ns, the, psi, ip, iok)
  call i2mex_getRmajor(R0, iok)
 
  res = i2mex_fourpi_r8 * abs(res) /(R0 * ip**2)
 
  deallocate(x, xt, xp, &
       & zt, zp, integrand, stat=iok)
 
  if(iok /= 0) ier = 86
 
end subroutine i2mex_getLi
 
