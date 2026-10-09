subroutine i2mex_getDelta(nt1, ns, the, psi, res, ier)

  ! Return in res the difference 'delta' between the 
  ! PEST1 straight field line poloidal angle and 'the' 
  ! on the (the, psi) grid.

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
  real(i2mex_r8), intent(out) :: res(nt1, ns)
  integer, intent(out) :: ier

  integer i, nt, j, iok, it
  real(i2mex_r8) g_ori(i2mex_o%ns), q_ori(i2mex_o%ns)
  real(i2mex_r8) psi_ori(i2mex_o%ns), the_ori(i2mex_o%nt1)
  real(i2mex_r8), dimension(:,:), allocatable :: x, xjacob, integrand, result, temp
  real(i2mex_r8) exponent

  ier = 0

  if(i2mex_mx==2 .and.  i2mex_npsi==0 .and. i2mex_kb==0) then
     ! Pest-1 coordinates
     res(1:nt1, 1:ns) = 0.0_i2mex_r8
     return
  endif

  allocate(temp(nt1, ns))

  nt = nt1 - 1

  call i2mex_getOriT(i2mex_o%nt1, the_ori, iok)
  call i2mex_error(iok)
  call i2mex_getOriPsi(i2mex_o%ns, psi_ori, iok)
  call i2mex_error(iok)

  allocate(integrand(i2mex_o%nt1, i2mex_o%ns))
  allocate(x(i2mex_o%nt1, i2mex_o%ns))
  allocate(xjacob(i2mex_o%nt1, i2mex_o%ns))

  call i2mex_getQ(i2mex_o%ns, psi_ori, q_ori, iok)
  call i2mex_error(iok)
  call i2mex_getG(i2mex_o%ns, psi_ori, g_ori, iok)
  call i2mex_error(iok)
  call i2mex_getX(i2mex_o%nt1, i2mex_o%ns, the_ori, psi_ori, x     , iok)
  call i2mex_error(iok)
  call i2mex_getJ(i2mex_o%nt1, i2mex_o%ns, the_ori, psi_ori, xjacob, iok)
  call i2mex_error(iok)

  do i = 1, i2mex_o%ns
     integrand(:, i) = g_ori(i)*xjacob(:,i)/(q_ori(i)* x(:,i)**2)
  enddo

  call i2mex_integrateAlongFluxSurface(integrand, nt1, ns, the, psi, temp, iok)

!!$  call i2mex_toAxisFraction(nt1, ns, the, psi, +0.5_i2mex_r8, &
!!$       & i2mex_to_axis_fraction, temp, iok)
!!$  call i2mex_toEdgeFraction(nt1, ns, the, psi,  &
!!$       & i2mex_to_edge_fraction, temp, iok)

!!$  exponent = 0.5_i2mex_r8
  exponent = 0.0_i2mex_r8
!!$  call i2mex_toAxis(nt1, ns, the, psi, exponent, temp, ier)
  call i2mex_toAxisFraction(nt1, ns, the, psi, +0._i2mex_r8, &
       & i2mex_to_axis_fraction, temp, iok)
  call i2mex_toEdgeFraction(nt1, ns, the, psi,  &
       & i2mex_to_edge_fraction, temp, iok)

  ! subtract angle
  do i = 1, ns
     temp(:,i) = temp(:,i) - the
  enddo

  ! regularize
!!$  temp(1:nt1,1) = sum(temp(1:nt1,1))/real(nt1, i2mex_r8)

  call i2mex_theta_orient(it, iok)

  if(it == i2mex_counterclockwise) then
     print *,'keep theta orientation in getDelta'
     res(1:nt1, 1:ns) = temp(1:nt1, 1:ns)
  else
     print *,'reverse theta orientation in getDelta'
     do i = 1, ns
        do j = 1, nt1
           res(j, i) = - temp(nt1-j+1, i)
        enddo
     enddo
  endif 

  deallocate(temp)

  deallocate(integrand)
  deallocate(x)
  deallocate(xjacob)

  if(iok/=0) ier = 142

end subroutine i2mex_getDelta

subroutine i2mex_getQDeltaT(nt1, ns, the, psi, res, ier)

  ! Return in res the derivative d /dthe of  
  ! 'q*delta' where delta is the difference between 
  ! the PEST1 straight field line poloidal angle and 
  ! 'the' the poloidal angle on the (the, psi) grid.

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
  real(i2mex_r8), intent(out) :: res(nt1, ns)
  integer, intent(out) :: ier

  integer i, nt, j, iok
  real(i2mex_r8) g(ns), q(ns)
  real(i2mex_r8), dimension(:,:), allocatable :: x, xjacob

  ier = 0

  allocate( x(nt1, ns) )
  allocate( xjacob(nt1, ns) )

  call i2mex_getG(ns, psi, g, iok)
  call i2mex_getQ(ns, psi, q, iok)
  call i2mex_getX(nt1, ns, the, psi, x, iok)
  call i2mex_getJ(nt1, ns, the, psi, xjacob, iok)

  res = spread(g, dim=1, ncopies=nt1) * xjacob / (x*x) - &
       & spread(q, dim=1, ncopies=nt1)

  deallocate(x)
  deallocate(xjacob)
  
  if(iok/=0) ier = 163

end subroutine i2mex_getQDeltaT


subroutine i2mex_getQDeltaP(nt1, ns, the, psi, res, ier)

  ! Return in res the derivative d /dpsi of  
  ! 'q*delta' where delta is the difference between 
  ! the PEST1 straight field line poloidal angle and 
  ! 'the' the poloidal angle on the (the, psi) grid.

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
  real(i2mex_r8), intent(out) :: res(nt1, ns)
  integer, intent(out) :: ier

  integer i, nt, j, iok, it
  real(i2mex_r8) g_ori(i2mex_o%ns)
  real(i2mex_r8) g_oriP(i2mex_o%ns), qP(ns)
  real(i2mex_r8) psi_ori(i2mex_o%ns), the_ori(i2mex_o%nt1)
  real(i2mex_r8), dimension(:,:), allocatable :: x, xjacob, integrand, result
  real(i2mex_r8), dimension(:,:), allocatable :: xT, xP, xjacobP
  real(i2mex_r8), dimension(:,:), allocatable :: xsq, temp

  real(i2mex_r8), dimension(ns) :: g, gp
  integer npoints
  real(i2mex_r8) exponent

  ier = 0

  if(i2mex_mx==2 .and.  i2mex_npsi==0 .and. i2mex_kb==0) then
     ! Pest-1 coordinates
     res(1:nt1, 1:ns) = 0.0_i2mex_r8
     return
  endif

  allocate(temp(nt1, ns))

  nt = nt1 - 1

!!$  call i2mex_getOriT(i2mex_o%nt1, the_ori, iok)
!!$  call i2mex_error(iok)
!!$  call i2mex_getOriPsi(i2mex_o%ns, psi_ori, iok)
!!$  call i2mex_error(iok)
!!$
!!$  allocate(integrand(i2mex_o%nt1, i2mex_o%ns))
!!$  allocate(x(i2mex_o%nt1, i2mex_o%ns))
!!$  allocate(xsq(i2mex_o%nt1, i2mex_o%ns))
!!$  allocate(xT(i2mex_o%nt1, i2mex_o%ns))
!!$  allocate(xP(i2mex_o%nt1, i2mex_o%ns))
!!$  allocate(xjacob(i2mex_o%nt1, i2mex_o%ns))
!!$  allocate(xjacobP(i2mex_o%nt1, i2mex_o%ns))
!!$
!!$  call i2mex_getG(i2mex_o%ns, psi_ori, g_ori, iok)
!!$  call i2mex_error(iok)
!!$
!!$  call i2mex_getGP(i2mex_o%ns, psi_ori, g_oriP, iok)
!!$  call i2mex_error(iok)
!!$
!!$  call i2mex_getX(i2mex_o%nt1, i2mex_o%ns, the_ori, psi_ori, x     , iok)
!!$  call i2mex_error(iok)
!!$  call i2mex_getGradX(i2mex_o%nt1, i2mex_o%ns, the_ori, psi_ori, xT, xP, iok)
!!$  call i2mex_error(iok)
!!$
!!$  call i2mex_getJ(i2mex_o%nt1, i2mex_o%ns, the_ori, psi_ori, xjacob, iok)
!!$  call i2mex_error(iok)
!!$  call i2mex_getJP(i2mex_o%nt1, i2mex_o%ns, the_ori, psi_ori, xjacobP, iok)
!!$  call i2mex_error(iok)
!!$
!!$  xsq = x**2
!!$
!!$  do i = 1, i2mex_o%ns
!!$     integrand(:,i) = g_orip(i)*xjacob(:,i)/xsq(:,i) &
!!$          & + g_ori(i)*xjacobp(:,i)/xsq(:,i) &
!!$          & - g_ori(i)*xjacob(:,i)*2._i2mex_r8*xp(:,i)/x(:,i)**3
!!$  enddo
!!$
!!$  exponent = -0.5_i2mex_r8 * i2mex_npsi
!!$  call i2mex_toAxis(i2mex_o%nt1, i2mex_o%ns, the_ori, psi_ori, &
!!$       & exponent, integrand, iok)
!!$  call i2mex_error(iok)
!!$  ! !!$  do i=1, i2mex_o%nt1
!!$  ! !!$     integrand(i,1) = FCCCC0(integrand(i,2), &
!!$  ! !!$          & integrand(i,3), integrand(i,4), integrand(i,5), &
!!$  ! !!$          & psi_ori(2),  psi_ori(3),  psi_ori(4),  psi_ori(5), &
!!$  ! !!$          & psi_ori(1))
!!$  ! !!$  enddo
!!$  integrand(:,1) = sum(integrand(:,1))/real(i2mex_o%nt1, i2mex_r8)
!!$
!!$  call i2mex_integrateAlongFluxSurface(integrand, nt1, ns, the, psi, temp, iok)


  allocate(integrand(nt1, ns))
  allocate(x(nt1, ns))
  allocate(xsq(nt1, ns))
  allocate(xT(nt1, ns))
  allocate(xP(nt1, ns))
  allocate(xjacob(nt1, ns))
  allocate(xjacobP(nt1, ns))
  
  call i2mex_getG(ns, psi, g, iok)
  call i2mex_error(iok)

  call i2mex_getGP(ns, psi, gP, iok)
  call i2mex_error(iok)

  call i2mex_getX(nt1, ns, the, psi, x     , iok)
  call i2mex_error(iok)
  call i2mex_getGradX(nt1, ns, the, psi, xT, xP, iok)
  call i2mex_error(iok)

  call i2mex_getJ(nt1, ns, the, psi, xjacob, iok)
  call i2mex_error(iok)
  call i2mex_getJP(nt1, ns, the, psi, xjacobP, iok)
  call i2mex_error(iok)

  xsq = x**2

  do i=1, ns
     integrand(1:nt1,i) = gp(i)*xjacob(1:nt1,i)/xsq(1:nt1,i) &
          & + g(i)*xjacobp(1:nt1,i)/xsq(1:nt1,i) &
          & - g(i)*xjacob(1:nt1,i)*2._i2mex_r8*xp(1:nt1,i)/x(1:nt1,i)**3
     call i2mex_integratePeriodic1d(nt1, the, integrand(1:nt1,i), &
          & nt1, the, temp(1,i))
  enddo



  call i2mex_getQP(ns, psi, qp, iok)
  call i2mex_error(iok)
  ! subtract angle
  do i = 1, ns
     temp(:,i) = temp(:,i) - qP(i)* the
  enddo

  exponent = -0.5_i2mex_r8  * i2mex_npsi
  call i2mex_toAxis(nt1, ns, the, psi, exponent, temp, ier)
!!$
  call i2mex_toAxisFraction(nt1, ns, the, psi, -0.5_i2mex_r8, &
       & i2mex_to_axis_fraction, temp, iok)
  call i2mex_toEdgeFraction(nt1, ns, the, psi,  &
       & i2mex_to_edge_fraction, temp, iok)
!!$
!!$  ! regularize
!!$  temp(1:nt1,1) = sum(temp(1:nt1,1))/real(nt1, i2mex_r8)

  call i2mex_theta_orient(it, iok)

  if(it == i2mex_counterclockwise) then
     print *,'keep theta orientation in getQDeltaQ'
     res(1:nt1, 1:ns) = temp(1:nt1, 1:ns)
  else
     print *,'reverse theta orientation in getQDeltaQ'
     do i = 1, ns
        do j = 1, nt1
           res(j, i) = - temp(nt1-j+1, i)
        enddo
     enddo
  endif 
  


  deallocate(temp)

  deallocate(integrand)
  deallocate(x)
  deallocate(xsq)
  deallocate(xT)
  deallocate(xP)
  deallocate(xjacob)
  deallocate(xjacobP)

  if(iok/=0) ier = 143

end subroutine i2mex_getQDeltaP
