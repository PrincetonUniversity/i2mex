  subroutine i2mex_getDelStar(nt1, ns, the, psi, res, ier)

    ! Return the Grad-Shafranov operator (Del^*) in res

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1)
    real(i2mex_r8), intent(in) :: psi(ns)
    real(i2mex_r8), intent(out) :: res(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8), dimension(ns) :: pp, g, gp
    real(i2mex_r8), dimension(:, :), allocatable :: x
    real(i2mex_r8), dimension(:, :), allocatable :: xt, xp, zt, zp
    real(i2mex_r8), dimension(:, :), allocatable :: xtt, xpt, xpp
    real(i2mex_r8), dimension(:, :), allocatable :: ztt, zpt, zpp
    real(i2mex_r8), dimension(:, :), allocatable :: d
    integer iok, i, j
    real(i2mex_r8), parameter :: SMALL=1.e-6, LARGE=1.e+6


    allocate(x(nt1, ns), xt(nt1, ns), xp(nt1, ns))
    allocate(zt(nt1, ns), zp(nt1, ns), xtt(nt1, ns))
    allocate(xpt(nt1, ns), xpp(nt1, ns), ztt(nt1, ns))
    allocate(zpt(nt1, ns), zpp(nt1, ns), d(nt1, ns))

    ier = 0
    call i2mex_getX(nt1, ns, the, psi, x, iok)
    call i2mex_error(iok)
    call i2mex_getGradX(nt1, ns, the, psi, xt, xp, iok)
    call i2mex_error(iok)
    call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, iok)
    call i2mex_error(iok)
    call i2mex_getXtt(nt1, ns, the, psi, xtt, iok)
    call i2mex_error(iok)
    call i2mex_getXpt(nt1, ns, the, psi, xpt, iok)
    call i2mex_error(iok)
    call i2mex_getXpp(nt1, ns, the, psi, xpp, iok)
    call i2mex_error(iok)
    call i2mex_getZtt(nt1, ns, the, psi, ztt, iok)
    call i2mex_error(iok)
    call i2mex_getZpt(nt1, ns, the, psi, zpt, iok)
    call i2mex_error(iok)
    call i2mex_getZpp(nt1, ns, the, psi, zpp, iok)
    call i2mex_error(iok)    

    d = xp*zt - xt*zp
    call i2mex_toAxis(nt1, ns, the, psi, +0.0_i2mex_r8, d, iok)

    ! laplacian^* psi

    res =          (-(d**2*zt) + d*x*(xpt*xt - xp*xtt + zpt*zt - zp*ztt) + &
     &     x*(-(xp*xt*xtt*zp) + xt**3*zpp - 2*xp*xt**2*zpt - xpp*xt**2*zt - &
     &        xtt*zp**2*zt - xt*zp*zpt*zt + xt*zpp*zt**2 - xp*zpt*zt**2 - &
     &        xpp*zt**3 + xpt*(xt**2*zp + xp*xt*zt + 2*zp*zt**2) + &
     &        xp*(xp*xt + zp*zt)*ztt))/(d**3*x)


    if(iok/=0) ier = 60

    deallocate(x, xt, xp)
    deallocate(zt, zp, xtt)
    deallocate(xpt, xpp, ztt)
    deallocate(zpt, zpp, d)

  end subroutine i2mex_getDelStar

  subroutine i2mex_getXJphi(nt1, ns, the, psi, res, ier)

    ! Return the right hand side (-p'*X^2 - g g') of the 
    ! Grad-Shafranov equation in res

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1)
    real(i2mex_r8), intent(in) :: psi(ns)
    real(i2mex_r8), intent(out) :: res(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8), dimension(ns) :: pp, g, gp
    real(i2mex_r8), dimension(:,:), allocatable :: x, z
    integer iok, i

    ier = 0

    allocate(x(nt1, ns), z(nt1, ns))

    call i2mex_getX(nt1, ns, the, psi, x, iok)
    call i2mex_error(iok)

    call i2mex_getPp(ns, psi, pp, iok)
    call i2mex_error(iok)
    call i2mex_getG(ns, psi, g, iok)
    call i2mex_error(iok)
    call i2mex_getGp(ns, psi, gp, iok)
    call i2mex_error(iok)
    
    do i=1, nt1
       res(i,1:ns) = - pp(1:ns) * x(i,1:ns)**2 -  g(1:ns)*gp(1:ns)
    enddo

    if(iok/=0) ier = 61

    deallocate(x, z)

  end subroutine i2mex_getXJphi

  subroutine i2mex_getGsError(nt1, ns, the, psi, res, ier)

    ! Return the Grad-Shafranov equation error in res. The error 
    ! is normalized to R J_phi on axis.

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1)
    real(i2mex_r8), intent(in) :: psi(ns)
    real(i2mex_r8), intent(out) :: res(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8), dimension(:,:), allocatable :: zdum
    integer iok, i

    ier = 0

    allocate(zdum(nt1, ns))

    call i2mex_getDelStar(nt1, ns, the, psi, res, iok)
    call i2mex_error(iok)
    call i2mex_getXJphi(nt1, ns, the, psi, zdum, iok)
    call i2mex_error(iok)

    if(iok/=0) ier = 62

    res(1:nt1,1) = 0.0_i2mex_r8
    res(1:nt1,2:ns) = real(nt1-1,i2mex_r8)* (res(1:nt1,2:ns) - zdum(1:nt1,2:ns)) &
         & /SUM(zdum(1:nt1-1,1))

    deallocate(zdum)

  end subroutine i2mex_getGsError
