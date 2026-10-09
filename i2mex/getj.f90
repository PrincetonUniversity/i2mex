  subroutine i2mex_getJ(nt1, ns, the, psi, jac, ier)

    ! Get Jacobian 'jac' on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: jac(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8), dimension(:,:), allocatable :: xt, xp, zt, zp, xa
    real(i2mex_r8) signj
    integer iok

    integer i

    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    real(i2mex_r8), dimension(:), allocatable :: s

    ier=0

    if(i2mex_xzakima > 0) then

       allocate(ss(nt1, ns), tt(nt1, ns))
       allocate(s(ns))
       call i2mex_getS(ns, psi, s, iok)
       ss = spread(  s, dim=1, ncopies=nt1)
       tt = spread(the, dim=2, ncopies=ns )
       call eq_frhochi(nt1*ns,ss,tt,1,i2mex_o%id_jakima,nt1*ns,jac,iok)
       deallocate(ss, tt)
       deallocate(s)

    else
    
       allocate( xt(nt1, ns), xp(nt1, ns), &
            & zt(nt1, ns), zp(nt1, ns), xa(nt1, ns), stat=iok)
       if(iok/=0) ier = 44

       call i2mex_getX(nt1, ns, the, psi, xa, iok)
       call i2mex_getGradX(nt1, ns, the, psi, xt, xp, iok)
       call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, iok)
       signj = 2._i2mex_r8*(i2mex_o%isThetaClockwise-0.5_i2mex_r8)
!!$       jac = -(xp*zt - xt*zp)*xa * signj
       jac = (xp*zt - xt*zp)*xa * signj


!!$  call i2mex_toAxisFraction(nt1, ns, the, psi, 0.0_i2mex_r8, &
!!$       & i2mex_to_axis_fraction, jac, iok)
!!$  call i2mex_toEdgeFraction(nt1, ns, the, psi,  &
!!$       & i2mex_to_edge_fraction, jac, iok)
       call i2mex_toAxis(nt1, ns, the, psi,  0.0_i2mex_r8, jac, iok)    
       call i2mex_error(iok)
       jac(:,1) = abs(jac(:,1)*jac(:,2))/jac(:,2) ! enforce sign

!!$    jac = max(0._i2mex_r8, jac)

       deallocate( xt, xp, &
            & zt, zp, xa )

    endif

    if(iok /= 0) ier = 41

  end subroutine i2mex_getJ

  subroutine i2mex_getJt(nt1, ns, the, psi, jt, ier)

    ! Get d Jacobian /d theta 'jt' on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: jt(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8), dimension(:,:), allocatable :: xt, xp, zt, zp, xa
    real(i2mex_r8), dimension(:,:), allocatable :: xtt, xpt, ztt, zpt
    real(i2mex_r8) signj
    integer iok

    integer i, j

    real(i2mex_r8), dimension(:,:,:), allocatable :: dj
    real(i2mex_r8), dimension(:, :), allocatable :: tt, ss
    real(i2mex_r8), dimension(:), allocatable :: s, dsdp

    ier=0
    
    if (i2mex_xzakima > 0) then

       allocate(dj(nt1, ns, 2))
       allocate(ss(nt1,ns), tt(nt1,ns), s(ns), dsdp(ns))
       call i2mex_getS(ns, psi, s, iok)
       ss = spread(  s, dim=1, ncopies=nt1)
       tt = spread(the, dim=2, ncopies=ns )

       ! use  Akima representation for J to compute dJ/dtheta
       call eq_grhochi(nt1*ns,ss,tt,1,i2mex_o%id_jakima,nt1*ns,dj,iok)
       do j = 1, ns
          do i = 1, nt1
             jt(i,j) = dj(i,j,2)
          enddo
       enddo
       deallocate(dj)
       deallocate(ss, tt, s, dsdp)

    else

       allocate( xt(nt1, ns), xp(nt1, ns), &
            & zt(nt1, ns), zp(nt1, ns), xa(nt1, ns), stat=iok)
       if(iok/=0) ier = 45
       allocate( xtt(nt1, ns), xpt(nt1, ns), ztt(nt1, ns), zpt(nt1, ns), &
            & stat=iok)
       if(iok/=0) ier = 45

       call i2mex_getX(nt1, ns, the, psi, xa, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getGradX(nt1, ns, the, psi, xt, xp, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getXtt(nt1, ns, the, psi, xtt, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getXpt(nt1, ns, the, psi, xpt, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getZtt(nt1, ns, the, psi, ztt, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getZpt(nt1, ns, the, psi, zpt, iok)
       if(iok/=0) call i2mex_error(iok)
       signj = 2._i2mex_r8*(i2mex_o%isThetaClockwise-0.5_i2mex_r8)

!!$       jt = signj*( &
!!$            & -(xp*zt - xt*zp)*xt - (xpt*zt - xtt*zp + xp*ztt - xt*zpt)*xa &
!!$            & )
       jt = signj*( &
            & (xp*zt - xt*zp)*xt + (xpt*zt - xtt*zp + xp*ztt - xt*zpt)*xa &
            & )

       call i2mex_toAxis(nt1, ns, the, psi,  0.0_i2mex_r8, jt, iok)    
       call i2mex_error(iok)

!!$  call i2mex_toAxisFraction(nt1, ns, the, psi, 0.0_i2mex_r8, &
!!$       & i2mex_to_axis_fraction, jt, iok)
!!$  call i2mex_toEdgeFraction(nt1, ns, the, psi,  &
!!$       & i2mex_to_edge_fraction, jt, iok)

       deallocate( xt, xp, &
            & zt, zp, xa )
       deallocate( xtt, xpt, ztt, zpt )

    endif

    if(iok /= 0) ier = 42

  end subroutine i2mex_getJt

  subroutine i2mex_getJp(nt1, ns, the, psi, jp, ier)

    ! Get d Jacobian /d psi 'jp' on the new (the, psi) mesh

    use i2mex_mod

    use ezspline_obj
    use ezspline

    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: jp(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8), dimension(:,:), allocatable :: xt, xp, zt, zp, xa
    real(i2mex_r8), dimension(:,:), allocatable :: xpt, xpp, zpt, zpp
    real(i2mex_r8) signj
    integer iok

    integer i, j

!!$    type(ezspline2_r8) :: jspl
!!$    real(i2mex_r8), dimension(:,:), allocatable :: jac

    real(i2mex_r8), dimension(:,:,:), allocatable :: dj
    real(i2mex_r8), dimension(:, :), allocatable :: tt, ss
    real(i2mex_r8), dimension(:), allocatable :: s, dsdp

    ier=0


!!$    call ezspline_init(jspl, nt1, ns, (/-1, -1/), (/0, 0/), iok)
!!$    call ezspline_error(iok)
!!$    jspl%x1 = the
!!$    jspl%x2 = psi
!!$    jspl%isHermite = 0
!!$    allocate(jac(nt1, ns))
!!$    call i2mex_getJ(nt1, ns, the, psi, jac, iok)
!!$    call i2mex_error(iok)
!!$    call ezspline_setup(jspl, jac, iok)
!!$    deallocate(jac)
!!$    call ezspline_error(iok)
!!$    call ezspline_derivative(jspl, 0, 1, nt1, ns, the, psi, jp, iok)
!!$    call ezspline_error(iok)
!!$    call ezspline_free(jspl, iok)
!!$    call ezspline_error(iok)

    if (i2mex_xzakima > 0) then

       allocate(dj(nt1, ns, 2))
       allocate(ss(nt1,ns), tt(nt1,ns), s(ns), dsdp(ns))
       call i2mex_getS(ns, psi, s, iok)
       ss = spread(  s, dim=1, ncopies=nt1)
       tt = spread(the, dim=2, ncopies=ns )

       ! use  Akima representation for J to compute dJ/dpsi
       call eq_grhochi(nt1*ns,ss,tt,1,i2mex_o%id_jakima,nt1*ns,dj,iok)
       call i2mex_getDS(ns, psi, dsdp, iok)
       do j = 1, ns
          do i = 1, nt1
             jp(i,j) = dj(i,j,1)*dsdp(j)
          enddo
       enddo
       deallocate(dj)
       deallocate(ss, tt, s, dsdp)

    else

       ! compute j-prime from x, z representation. Can lead to some ringing..
    
       allocate( xt(nt1, ns), xp(nt1, ns), &
            & zt(nt1, ns), zp(nt1, ns), xa(nt1, ns), stat=iok )
       if(iok/=0) ier = 46
       allocate( xpt(nt1, ns), xpp(nt1, ns), zpt(nt1, ns), zpp(nt1, ns), &
            & stat=iok )
       if(iok/=0) ier = 46

       call i2mex_getX(nt1, ns, the, psi, xa, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getGradX(nt1, ns, the, psi, xt, xp, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getGradZ(nt1, ns, the, psi, zt, zp, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getXpt(nt1, ns, the, psi, xpt, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getXpp(nt1, ns, the, psi, xpp, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getZpt(nt1, ns, the, psi, zpt, iok)
       if(iok/=0) call i2mex_error(iok)
       call i2mex_getZpp(nt1, ns, the, psi, zpp, iok)
       if(iok/=0) call i2mex_error(iok)
       signj = 2._i2mex_r8*(i2mex_o%isThetaClockwise-0.5_i2mex_r8)

!!$       jp = signj*( &
!!$            & -(xp*zt - xt*zp)*xp - (xpp*zt - xpt*zp + xp*zpt - xt*zpp)*xa &
!!$            & )
       jp = signj*( &
            & (xp*zt - xt*zp)*xp + (xpp*zt - xpt*zp + xp*zpt - xt*zpp)*xa &
            & )

       deallocate( xt, xp, &
            & zt, zp, xa )
       deallocate( xpt, xpp, zpt, zpp )

    endif

!!$    call i2mex_toAxis(nt1, ns, the, psi,  0.0_i2mex_r8, jp, iok)    
!!$    call i2mex_error(iok)

  call i2mex_toAxisFraction(nt1, ns, the, psi, 0.0_i2mex_r8, &
       & i2mex_to_axis_fraction, jp, iok)
  call i2mex_toEdgeFraction(nt1, ns, the, psi,  &
       & i2mex_to_edge_fraction, jp, iok)

    if(iok /= 0) ier = 43

  end subroutine i2mex_getJp
