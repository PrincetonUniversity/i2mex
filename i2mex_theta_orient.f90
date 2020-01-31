subroutine i2mex_theta_orient(it_orient, ier)

  ! Return in it_orient the poloidal angle orientation: 
  ! either i2mex_clockwise or i2mex_counterclockwise

  use i2mex_mod
  implicit none
  integer, intent(out) :: it_orient
  integer, intent(out) :: ier

  integer ns, nt1, iok
  real(i2mex_r8), allocatable :: psi(:), the(:)
  real(i2mex_r8), allocatable :: x(:,:), z(:,:)
  real(i2mex_r8) det, x0, x1, x2, z0, z1, z2

  ier = 0

  it_orient = i2mex_clockwise

  call i2mex_getOriNs(ns, iok)
  call i2mex_getOriNt1(nt1, iok)

  allocate(psi(ns), the(nt1))
  allocate(x(nt1,ns), z(nt1,ns))

  call i2mex_getOriPsi(ns, psi, iok)
  call i2mex_getOriT(nt1, the, iok)

  call i2mex_getX(nt1, ns, the, psi, x, iok)
  call i2mex_getZ(nt1, ns, the, psi, z, iok)

  x0 = x(1,1)
  z0 = z(1,1)
  x1 = x(1, ns)
  z1 = z(1, ns)
  x2 = x(2, ns)
  z2 = z(2, ns)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)

  if(det < 0._i2mex_r8) then
     it_orient = i2mex_clockwise
  else
     it_orient = i2mex_counterclockwise
  endif
  
  
  

end subroutine i2mex_theta_orient
