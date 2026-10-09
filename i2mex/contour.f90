subroutine i2mex_contour(nx, ny, x, y, f, x0, y0, xp, yp, ns, ttot, xs, ys, ier)
  ! History:
  ! 04/22/08 CLF: allocate temp arrays
  !
  ! Given a fct F(x, y) , return the ns values of (xs, ys) where 
  ! f(xs, ys) = f(xp, yp). Ttot is the total time variable after a close contour
  ! (some codes need ttot to compute the rotational transform along the contour).
  ! Method: Hamilton's equation dx/dt=x*dF/dy; dy/dt=-x*dF/fx are integrated
  ! starting from (xp, yp) until the cummulative angle between (x-x0, y-y0) 
  ! and (xnew-x0, ynew-y0) reaches 2*pi, where (xnew, ynew) are new coordinates
  ! along the trajectory.
  ! As a second step, the (xnew, ynew) coordinates are interpolated onto 
  ! the equidistant angle grid to yield the (xs, ys).

  
  use cont_mod
  use imex_ode_mod

  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  integer, intent(in) :: nx, ny
  real(r8), intent(in) :: x(nx), y(ny)
  real(r8), intent(in) :: f(nx, ny)
  real(r8), intent(in) :: x0, y0 ! approximate contour centre
  real(r8), intent(in) :: xp, yp ! starting point 
  integer, intent(in) :: ns
  real(r8), intent(out) :: ttot, xs(ns), ys(ns)
  integer, intent(out) :: ier

  integer, parameter :: neq=2, mf=10
  integer istate, iwork(5)
  real(r8) atol, rtol

  real(r8) :: fs(ns) ! diagnostics

  real(r8), parameter :: twopi=6.28318530717958623_r8, eps=1.e-12_r8
  integer, parameter :: nsteps = 2000000 ! max number of integration steps
                                         ! RGA, increased from 1000000
  integer i, iok, npoints, iok2
  real(r8) :: lenscale, t, dt, tout, angle, a, const, vecprod, scaprod
  real(r8), dimension(ns) :: as
  real(r8), dimension(:, :), allocatable :: xy
  real(r8), dimension(:),    allocatable :: rwork, aa
  real(r8), dimension(:),    allocatable :: dx, dxnew, xx, xx0
  real(r8) :: diff, aver
  real(r8) :: delta(2)

  real(r8), parameter :: ndelta = 1024 ! no of steps (can be increased for more accuracy)
                                       ! RGA 12Oct2004, increased from 512
  type(ezspline1) :: xspl, yspl
  integer::             zsts

  external i2mex_xydot

  ier = 0


  call ezspline_init(fspl, nx, ny, (/0,0/), (/0,0/), iok)
  call ezspline_error(iok)

  fspl%x1 = x
  fspl%x2 = y

  call ezspline_setup(fspl, f, iok)
  call ezspline_error(iok)

  ! imex_ode stuff
  rtol=1.e-12_r8
  atol=1.e-12_r8
  istate = 1

  allocate(dx(neq), dxnew(neq), xx(neq), xx0(neq), aa(nsteps))
  allocate( xy(neq, nsteps),stat=zsts)
  xy = 0._r8
  xx0 = (/x0, y0/)
  xx  = (/xp, yp/)
  i = 1
  npoints= 1
  lenscale = sqrt( (xp-x0)**2 + (yp-y0)**2 )
  
  dt = twopi/(cont_q*real(ndelta,r8))

  angle = 0._r8
  const = 1._r8
  aa = 0._r8
  xy(:,1) = xx
  t = 0.0_r8

  ttot = 0. 
  xs = 0.
  ys = 0.
  
!!$  ! probe direction. We want to integrate counterclockwise.
!!$  dx = xx - xx0
!!$  tout = t + const*dt
!!$  cont_sign = 1._r8 ! integration sign in i2mex_xydot
!!$     call lsode_r8 (i2mex_xydot, neq, xx, t, tout, itol, rtol, atol, itask, &
!!$          & istate, iopt, rwork, lrw, iwork, liw, i2mex_jac, mf)
!!$   dxnew = xx - xx0
!!$   if(dx(1)*dxnew(2)-dx(2)*dxnew(1) > 0._r8) then
!!$      print *,'correct sign'
!!$      cont_sign = -1._r8
!!$   endif
!!$
!!$   t = 0.0_r8
!!$   xx = (/xp, yp/)
  
  allocate(rwork(100 + 21*neq),stat=zsts)

  do while (npoints < nsteps-1 .and. abs(const*dt) > eps)
     dx = xx - xx0
     tout = t + const*dt
     call imex_ode (i2mex_xydot, neq, xx, t, tout, rtol, atol, &
          & istate, rwork, iwork)
     if(istate/=2) then
        print *,'IMEX_ODE error', istate
     endif
     dxnew = xx - xx0
     vecprod = dx(1)*dxnew(2)-dx(2)*dxnew(1)
     scaprod = dx(1)*dxnew(1)+dx(2)*dxnew(2)
     a = atan2(vecprod, scaprod)
     if(a < 0._r8) then
        a = a + twopi
     endif
     if (a>twopi/2.) then
        print *, '?i2mex_contour: contour angle jumped by more then pi probably due to a flat spot'
        ier=1
        return
     end if
     angle = angle + a
     if(abs(angle) > twopi) then
        ! back off
        t = t - const*dt
        const = const/2._r8
        angle = angle - a
        xx = xy(:,npoints)
        ! reset imex_ode state
        istate = 1
     else if(const==1._r8) then
        npoints = npoints + 1 ! foward counter
        aa(npoints) = angle
        xy(:, npoints) = xx
        !t = t + const*dt
     endif
     i = i + 1
  enddo
  deallocate(rwork)

  ttot = t
!!$  ! linear error correction
!!$  delta = xy(:,npoints)-xy(:,1)
!!$  print *,'delta = ', sqrt(delta(1)**2 + delta(2)**2)
!!$  do i = 1, 2
!!$     xy(i,1:npoints) = xy(i,1:npoints) - delta(i)*aa(1:npoints)/twopi
!!$  enddo

  npoints = npoints + 1
  aa(npoints) = angle
  xy(:, npoints) = xx
  

  ! interpolate onto equal angle grid

  call ezspline_init(xspl, npoints, (/-1, -1/), iok)
  call ezspline_error(iok)
  call ezspline_init(yspl, npoints, (/-1, -1/), iok)
  call ezspline_error(iok)

  xspl%x1 = aa(1:npoints)
  yspl%x1 = aa(1:npoints)


  xy(:,npoints) = xy(:,1) ! make periodic

  call ezspline_setup(xspl, xy(1,1:npoints), iok)
  if(iok/=0) then
     ! check if x1 is non-ascendant, which would be symptomatic
     ! of an open contour
     do i = 1, npoints-1
        if(xspl%x1(i) >= xspl%x1(i+1) ) then
           print *,'ERROR: non-ascending abscissa value at i=',i, xspl%x1(i), ' >= ', xspl%x1(i+1), ' (x,y) = ', xy(1,i), xy(2,i)
        endif
     enddo
     print *,'magnetic axis: (x0, y0) = (',x0,',',y0,')'
     print *,'starting point: (xp, yp) = (',xp,',',yp,')'
  endif
  call ezspline_error(iok)
  call ezspline_setup(yspl, xy(2,1:npoints), iok)
  call ezspline_error(iok)

  as = twopi* (/ (real(i-1,r8)/real(ns-1,r8), i=1, ns) /)
  call ezspline_interp(xspl, ns, as, xs, iok)
  call ezspline_error(iok)
  call ezspline_interp(yspl, ns, as, ys, iok)
  call ezspline_error(iok)

  call ezspline_free(xspl, iok)
  call ezspline_error(iok)
  call ezspline_free(yspl, iok)
  call ezspline_error(iok)


  ! check

  call ezspline_interp(fspl, ns, xs, ys, fs, iok)
  call ezspline_error(iok)
  diff = 0.0_r8
  aver = sum(fs)/real(ns, r8)
  diff = sqrt(maxval((fs-aver)**2))
  if (aver/=0.0_r8) then
     if(diff/abs(aver) > 0.001_r8 .and. abs(aver) > 1.e-2_r8) then
        ier = 123
        print *,'i2mex_contour: diff=', diff, ' abs(aver)=', abs(aver)
     endif
  endif

  ! clean up

  call ezspline_free(fspl, iok2)
  call ezspline_error(iok2)
  deallocate(dx, dxnew, xx, xx0, aa)
  deallocate(xy)
  if(iok/=0 .or. iok2/=0) ier = 122
  if(istate/=2 .and. istate/=1) ier = 120
  if(npoints>=nsteps) ier = 121

end subroutine i2mex_contour


subroutine i2mex_xydot(t, y, ydot)

  use cont_mod

  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  integer iok
  real(r8) t, y(2), ydot(2)
  real(r8) df(2)

  call ezspline_gradient(fspl, y(1), y(2), df, iok)
  call ezspline_error(iok)

  ydot(1) = -cont_sign * y(1)*df(2)
  ydot(2) = +cont_sign * y(1)*df(1)

end subroutine i2mex_xydot
