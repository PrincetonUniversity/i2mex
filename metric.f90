  subroutine i2mex_getMetric(nt1, ns, the, psi, gtt, gtp, gpp, ier)
 
    ! Get (|grad the|^2, grad the . grad psi, |grad psi|^2)
    ! on the (the, psi) mesh
    !
    ! gtt = |grad the|^2
    ! gtp = grad the . grad psi
    ! gpp = |grad psi|^2
 
    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: gtt(nt1, ns), gtp(nt1, ns), gpp(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8), dimension(:,:), allocatable :: x, xt, xp
    real(i2mex_r8), dimension(:,:), allocatable :: z, zt, zp
    real(i2mex_r8), dimension(:,:), allocatable :: jac

    real(i2mex_r8), dimension(:,:,:), allocatable :: grad_psi
    real(i2mex_r8), dimension(:,:,:), allocatable :: grad_the

    integer iok

    integer i

    ier=0

    allocate(x(nt1, ns), z(nt1, ns))

    call i2mex_getX(nt1, ns, the, psi, x, ier)
    call i2mex_getZ(nt1, ns, the, psi, z, ier)


    allocate(&
         & xt(nt1, ns), xp(nt1, ns), &
         & zt(nt1, ns), zp(nt1, ns), &
         & jac(nt1, ns))

    call i2mex_getGradX(nt1, ns, the, psi, xt, xp, ier)
    call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, ier)
    call i2mex_getJ(nt1, ns, the, psi, jac, ier)

    if (i2mex_direct==0) then

       gpp = x*x*(xt**2 + zt**2)/jac**2

    else

       ! use direct representation
       allocate(grad_psi(nt1, ns, 2))

       print *,'*** use direct representation for grad psi '

       call eq_gRZ(nt1*ns, x(1,1), z(1,1), 1, i2mex_o%id_psi, nt1*ns, grad_psi, iok)
       if(iok/=0) ier = 155

       gpp = grad_psi(:,:,1)**2 + grad_psi(:,:,2)**2

       deallocate(grad_psi)

    endif

    gtp = -x*x*(zt*zp + xt*xp)/jac**2
    gtt = x*x*(zp*zp + xp*xp)/jac**2

    call i2mex_toAxis(nt1, ns, the, psi, -1.0_i2mex_r8, gtt, iok)
    call i2mex_error(iok)
    gtt(:,1) = abs(gtt(:,1)) ! enforce positive definiteness
    call i2mex_toAxis(nt1, ns, the, psi, +0.0_i2mex_r8, gtp, iok)
    call i2mex_error(iok)
    call i2mex_toAxis(nt1, ns, the, psi, +1.0_i2mex_r8, gpp, iok)
    call i2mex_error(iok)
    gpp(:,1) = abs(gpp(:,1)) ! enforce positive definiteness

    deallocate( &
         & xt, xp, &
         & zt, zp, &
         & jac)



    deallocate(x, z)
 
  end subroutine i2mex_getMetric
 
  subroutine i2mex_getMetricT(nt1, ns, the, psi, gtt, gtp, gpp, ier)
 
    ! Get poloidal derivative (d /d the) of metric quantities
    ! (|grad the|^2, grad the . grad psi, |grad psi|^2)
    ! on the (the, psi) mesh
    !
    ! gtt = (d /d the) |grad the|^2
    ! gtp = (d /d the) grad the . grad psi
    ! gpp = (d /d the) |grad psi|^2
 
    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: gtt(nt1, ns), gtp(nt1, ns), gpp(nt1, ns)
    integer, intent(out) :: ier
 
    real(i2mex_r8), dimension(:,:), allocatable :: x, xt, xp, xtt, xtp
    real(i2mex_r8), dimension(:,:), allocatable :: z, zt, zp, ztt, ztp
    real(i2mex_r8), dimension(:,:), allocatable :: d3
    integer iok
 
    integer i
 
    ier=0
 
    allocate( &
         & x(nt1, ns), xt(nt1, ns), xp(nt1, ns), xtt(nt1, ns), xtp(nt1, ns), &
         & z(nt1, ns), zt(nt1, ns), zp(nt1, ns), ztt(nt1, ns), ztp(nt1, ns), &
         & d3(nt1, ns))
 
    call i2mex_getX(nt1, ns, the, psi, x, ier)
    call i2mex_getGradX(nt1, ns, the, psi, xt, xp, ier)
    call i2mex_getXtt(nt1, ns, the, psi, xtt, ier)
    call i2mex_getXtp(nt1, ns, the, psi, xtp, ier)
 
    call i2mex_getZ(nt1, ns, the, psi, z, ier)
    call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, ier)
    call i2mex_getZtt(nt1, ns, the, psi, ztt, ier)
    call i2mex_getZtp(nt1, ns, the, psi, ztp, ier)
 
    d3 = (-(xt*zp) + xp*zt)**3
    call i2mex_toAxis(nt1, ns, the, psi, +0.0_i2mex_r8, d3, iok)
    call i2mex_error(iok)
    
 
    gtt = (-2*(xt*zp - xp*zt)* &
         &        (xp*xtp + zp*ztp) -  &
         &       2*(xp**2 + zp**2)* &
         &        (-(xtt*zp) + xtp*zt - xt*ztp +  &
         &          xp*ztt))/d3
 
    gtp = (2*(xp*xt + zp*zt)* &
         &        (-(xtt*zp) + xtp*zt - xt*ztp +  &
         &          xp*ztt) + (xt*zp - xp*zt)* &
         &        (xt*xtp + xp*xtt + zt*ztp +  &
         &          zp*ztt))/d3
 
    gpp = (2*(xt**2 + zt**2)* &
         &        (xtt*zp - xtp*zt + xt*ztp -  &
         &          xp*ztt) - 2*(xt*zp - xp*zt)* &
         &        (xt*xtt + zt*ztt))/ &
         &     d3
 
 
    call i2mex_toAxis(nt1, ns, the, psi, -1.0_i2mex_r8, gtt, iok)
    call i2mex_error(iok)
    call i2mex_toAxis(nt1, ns, the, psi, -0.0_i2mex_r8, gtp, iok)
    call i2mex_error(iok)
    call i2mex_toAxis(nt1, ns, the, psi, +1.0_i2mex_r8, gpp, iok)
    call i2mex_error(iok)
 
    deallocate( &
         & x, xt, xp, xtt, xtp, &
         & z, zt, zp, ztt, ztp, &
         & d3)
 
  end subroutine i2mex_getMetricT
 
 
  subroutine i2mex_getMetricP(nt1, ns, the, psi, gtt, gtp, gpp, ier)
 
    ! Get flux derivative (d /d psi) of metric quantities
    ! (|grad the|^2, grad the . grad psi, |grad psi|^2)
    ! on the (the, psi) mesh
    !
    ! gtt = (d /d the) |grad the|^2
    ! gtp = (d /d the) grad the . grad psi
    ! gpp = (d /d the) |grad psi|^2
 
    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: gtt(nt1, ns), gtp(nt1, ns), gpp(nt1, ns)
    integer, intent(out) :: ier
 
    real(i2mex_r8), dimension(:,:), allocatable :: x, xt, xp, xtp, xpp
    real(i2mex_r8), dimension(:,:), allocatable :: z, zt, zp, ztp, zpp
    real(i2mex_r8), dimension(:,:), allocatable :: d3
    integer iok
 
    integer i
 
    ier=0
 
    allocate( &
         & x(nt1, ns), xt(nt1, ns), xp(nt1, ns), xtp(nt1, ns), xpp(nt1, ns), &
         & z(nt1, ns), zt(nt1, ns), zp(nt1, ns), ztp(nt1, ns), zpp(nt1, ns), &
         & d3(nt1, ns))
 
    call i2mex_getX(nt1, ns, the, psi, x, ier)
    call i2mex_getGradX(nt1, ns, the, psi, xt, xp, ier)
    call i2mex_getXtp(nt1, ns, the, psi, xtp, ier)
    call i2mex_getXpp(nt1, ns, the, psi, xpp, ier)
 
    call i2mex_getZ(nt1, ns, the, psi, z, ier)
    call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, ier)
    call i2mex_getZtp(nt1, ns, the, psi, ztp, ier)
    call i2mex_getZpp(nt1, ns, the, psi, zpp, ier)
 
    d3 = (-(xt*zp) + xp*zt)**3
    call i2mex_toAxis(nt1, ns, the, psi, +0.0_i2mex_r8, d3, iok)
    call i2mex_error(iok)
    
    gtt = (2*(xp*xpp + zp*zpp)* &
         &        (-(xt*zp) + xp*zt) -  &
         &       2*(xp**2 + zp**2)* &
         &        (-(xtp*zp) - xt*zpp + xpp*zt +  &
         &          xp*ztp))/d3
    gtp = (2*(xp*xt + zp*zt)* &
         &        (-(xtp*zp) - xt*zpp + xpp*zt +  &
         &          xp*ztp) + (xt*zp - xp*zt)* &
         &        (xpp*xt + xp*xtp + zpp*zt +  &
         &          zp*ztp))/d3
    gpp = (2*(xt**2 + zt**2)* &
         &        (xtp*zp + xt*zpp - xpp*zt -  &
         &          xp*ztp) - 2*(xt*zp - xp*zt)* &
         &        (xt*xtp + zt*ztp))/ &
         &     d3
 
 
    call i2mex_toAxis(nt1, ns, the, psi, -2.0_i2mex_r8, gtt, iok)
    call i2mex_error(iok)
    call i2mex_toAxis(nt1, ns, the, psi, -1.0_i2mex_r8, gtp, iok)
    call i2mex_error(iok)
    call i2mex_toAxis(nt1, ns, the, psi, +0.0_i2mex_r8, gpp, iok)
    call i2mex_error(iok)
 
    deallocate( &
         & x, xt, xp, xtp, xpp, &
         & z, zt, zp, ztp, zpp, &
         & d3)
 
  end subroutine i2mex_getMetricP
