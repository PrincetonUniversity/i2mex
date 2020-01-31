  subroutine i2mex_getX(nt1, ns, the, psi, x, ier)

    ! Get X on the (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: x(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    integer iok

    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    if(iok/=0) ier = 26

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_frhochi(nt1*ns, ss, tt, 1, i2mex_o%id_x, nt1*ns, x, iok)
    if(iok /= 0) ier = 21
    
    !!call i2mex_toAxis(nt1, ns, the, psi,  0.0_i2mex_r8, x, iok)
    !!call i2mex_error(iok)

    deallocate(tt)
    deallocate(ss)
    
  end subroutine i2mex_getX

  subroutine i2mex_getGradX(nt1, ns, the, psi, xt, xp, ier)

    ! Get (d X/d the, d X/ dpsi) on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: xt(nt1, ns), xp(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns), dsdpsi(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    real(i2mex_r8), dimension(:,:,:), allocatable :: dx
    integer iok, i
    
    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    allocate(dx(nt1, ns, 2), stat=iok)
    if(iok/=0) ier = 27

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_grhochi(nt1*ns,ss,tt,1,i2mex_o%id_x,nt1*ns,dx,iok)
    if(iok /= 0) ier = 22
        
    call i2mex_getDS(ns, psi, dsdpsi, iok)
    call i2mex_error(iok)

    xt(1:nt1, 1:ns) = dx(1:nt1, 1:ns,2)
    do i=1, nt1
       xp(i, 1:ns) = dx(i,1:ns,1)*dsdpsi(1:ns)
    enddo

    call i2mex_toAxis(nt1, ns, the, psi, +0.5_i2mex_r8, xt, iok)
    call i2mex_error(iok)
    call i2mex_toAxis(nt1, ns, the, psi, -0.5_i2mex_r8, xp, iok)    
    call i2mex_error(iok)
    
    deallocate(tt)
    deallocate(ss)
    deallocate(dx)
    
  end subroutine i2mex_getGradX

  subroutine i2mex_getXtt(nt1, ns, the, psi, xtt, ier)

    ! Get d^2 X/d the^2 on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: xtt(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    integer iok

    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    if(iok/=0) ier = 28

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_hrhochi(nt1*ns,ss,tt,1,i2mex_o%id_x,                    &
         &                (/0,0,0,0,1,0/),nt1*ns,xtt,iok)
    if(iok /= 0) ier = 23
    
    call i2mex_toAxis(nt1, ns, the, psi,  1.0_i2mex_r8, xtt, iok)
    call i2mex_error(iok)
    
    deallocate(tt)
    deallocate(ss)
    
  end subroutine i2mex_getXtt

  subroutine i2mex_getXtp(nt1, ns, the, psi, xpt, ier)

    ! Get d X/(d psi d the) on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: xpt(nt1, ns)
    integer, intent(out) :: ier


    call i2mex_getXpt(nt1, ns, the, psi, xpt, ier)

  end subroutine i2mex_getXtp

  subroutine i2mex_getXpt(nt1, ns, the, psi, xpt, ier)

    ! Get d X/(d psi d the) on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: xpt(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns), dsdpsi(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    integer iok, i

    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    if(iok/=0) ier = 29

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_hrhochi(nt1*ns,ss,tt,1,i2mex_o%id_x,                    &
         &                (/0,0,0,0,0,1/),nt1*ns,xpt,iok)
    if(iok /= 0) ier = 24

    call i2mex_getDS(ns, psi, dsdpsi, iok)
    call i2mex_error(iok)

    do i=1, nt1
       xpt(i, 1:ns) = xpt(i,1:ns)*dsdpsi(1:ns)
    enddo
    
    call i2mex_toAxis(nt1, ns, the, psi, -0.0_i2mex_r8, xpt, iok)
    call i2mex_error(iok)
    
    deallocate(tt)
    deallocate(ss)
    
  end subroutine i2mex_getXpt


  subroutine i2mex_getXpp(nt1, ns, the, psi, xpp, ier)

    ! Get d^2 X/d psi^2 on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: xpp(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns), dsdpsi(ns), d2sdpsi2(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    real(i2mex_r8), dimension(:,:,:), allocatable :: xs
    integer iok, i

    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    allocate(xs(nt1, ns, 2), stat=iok)
    if(iok/=0) ier = 30

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_hrhochi(nt1*ns,ss,tt,1,i2mex_o%id_x,                    &
         &                (/0,1,0,1,0,0/),nt1*ns,xs,iok)
    if(iok /= 0) ier = 25

    call i2mex_getDS(ns, psi, dsdpsi, iok)
    call i2mex_error(iok)
    call i2mex_getD2S(ns, psi, d2sdpsi2, iok)
    call i2mex_error(iok)
    do i=1, nt1
       xpp(i, 1:ns) = xs(i,1:ns,2)*dsdpsi(1:ns)**2 + &
            & xs(i,1:ns,1)*d2sdpsi2(1:ns)
    enddo
    
    call i2mex_toAxis(nt1, ns, the, psi,  -1.0_i2mex_r8, xpp, iok)
    call i2mex_error(iok)
    
    deallocate(tt)
    deallocate(ss)
    deallocate(xs)
    
  end subroutine i2mex_getXpp
