subroutine i2mex_poloidalRemap(mx, npsi, kb, nt1, ns, the, psi, g, x, z, ier)

  ! Poloidal regridding according to prescribed Jacobian dependence 
  ! J ~ x^mx /(|grad psi|^npsi B^kb). 
  !
  ! This a convenience routine typically used to remap the equilibrium onto 
  ! a grid with a different Jacobian dependence in x, |grad psi| and
  ! the magnetic field B. Some codes will require a specific dependence,
  ! among which 
  ! 
  ! mx npsi kb
  ! 2   0   0  PEST1 coordinates
  ! 1   1   0  Equal Arc
  ! 0   0   0  Hamada
  ! 0   0   2  Boozer 
  !
  ! are popular choices. For accuracy, the optimal choice may be equal-arc
  ! as it prevents excessive bending of the poloidal rays on the outer board.
  !
  ! *NOTE* This routine does not rely on a XPLASMA/I2MEX representation.
  ! It should typically be used BEFORE initializing the i2mex object.

  use ezspline_obj
  use ezspline
  implicit none
  integer, parameter  :: r8 = selected_real_kind(12,100)

  integer, intent(in) :: mx, npsi, kb ! the exponents entering the Jacobian
  integer, intent(in) :: nt1, ns  ! the poloidal/radial grid sizes
  real(r8), intent(in) :: the(nt1), psi(ns) ! the poloidal [rad] and radial [Wb/rad] grids
  real(r8), intent(in) :: g(ns) ! the covariant toroidal B function (Btor*X)
  real(r8), intent(inout) :: x(nt1, ns), z(nt1,ns) ! old and new grid arrays
  integer, intent(out) :: ier ! error flag

  real(r8), parameter :: pi = 3.1415926535897931_r8

  real(r8) :: gradPsi, B, jac, x0, z0
  real(r8), dimension(:), allocatable :: theNew, the_old
  real(r8), dimension(:), allocatable :: psin, s, ds
  real(r8), dimension(:,:), allocatable :: integrand
  real(r8), dimension(:,:,:), allocatable :: dx, dz
  type(ezspline2) :: xspl, zspl
  integer iok, i, j

  include 'cubinterp.h'
  
  ier = 0

  allocate(theNew(nt1), the_old(nt1))
  allocate(psin(ns), s(ns),ds(ns))
  allocate(integrand(nt1, ns))
  allocate(dx(nt1, ns, 2), dz(nt1, ns, 2))

  ! enforce periodicity
  x(nt1,:) = x(1,:)
  z(nt1,:) = z(1,:)

  ! axis 
  x0 = sum(x(:,1))/real(nt1,r8)
  z0 = sum(z(:,1))/real(nt1,r8)

  psin = abs(psi-psi(1))
  s = sqrt( psin/(psin(ns)-psin(1)) )
  ds(2:ns) = 0.5_r8 * s(2:ns) / psin(2:ns)
  ! extrapolate to axis
  ds(1) = FCCCC0(ds(2), ds(3), ds(4), ds(5), &
       &                 s(2),  s(3),  s(4),  s(5), &
       & s(1))

  ! x 
  call ezspline_init(xspl, nt1, ns, (/-1,-1/), (/0,0/), iok)
  call ezspline_error(iok)
  xspl%x1 = the
  xspl%x2 = s
  call ezspline_setup(xspl, x, iok)
  call ezspline_error(iok)
  call EZspline_gradient(xspl, nt1, ns, the, s, dx, iok)
  call ezspline_error(iok)

  if(iok/=0) ier = 132

  ! z
  call ezspline_init(zspl, nt1, ns, (/-1,-1/), (/0,0/), iok)
  call ezspline_error(iok)
  zspl%x1 = the
  zspl%x2 = s
  call ezspline_setup(zspl, z, iok)
  call ezspline_error(iok)
  call EZspline_gradient(zspl, nt1, ns, the, s, dz, iok)
  call ezspline_error(iok)

  if(iok/=0) ier = 133

  ! Jacobian, |gradPsi|, B
  do j=2, ns
     do i=1, nt1
        jac     = abs( x(i,j)*ds(j)*( dx(i,j,1)*dz(i,j,2) - dx(i,j,2)*dz(i,j,1) ) )
        gradPsi = x(i,j)*sqrt(dx(i,j,1)**2 + dz(i,j,1)**2)/jac
        b       = sqrt(gradPsi**2 + g(j)**2)/x(i,j)
        integrand(i,j) = jac* gradPsi**npsi * b**kb/ x(i,j)**mx
     enddo
  enddo
  ! special treatment on axis
  call i2mex_toAxis(nt1, ns, the, psi, 0.0_r8, integrand, ier)

  ! compute new theta on old grid and interpolate to get x and z on
  ! the new grid
  do j = 2, ns
     call i2mex_integratePeriodic1d(nt1, the, integrand(1,j), nt1, the, theNew)
     theNew = 2._r8*pi * theNew/theNew(nt1)

     ! compute old theta for new theta= old theta values
     call i2mex_interp1d(nt1, theNew, the, nt1, the, the_old)

     call EZspline_interp(xspl, nt1, 1, the_old, s(j:j), x(:,j:j), iok)
     call ezspline_error(iok)
     call EZspline_interp(zspl, nt1, 1, the_old, s(j:j), z(:,j:j), iok)
     call ezspline_error(iok)
  enddo

  if(iok/=0) ier = 134

  x(:,1) = x0
  z(:,1) = z0

  call ezspline_free(xspl, iok)
  call ezspline_free(zspl, iok)

  deallocate(theNew, the_old)
  deallocate(psin,s,ds)
  deallocate(integrand)
  deallocate(dx, dz)
  
end subroutine i2mex_poloidalRemap
