  subroutine i2mex_getZ(nt1, ns, the, psi, z, ier)

    ! Get Z on the (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: z(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    integer iok

    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    if(iok/=0) ier = 36

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_frhochi(nt1*ns, ss, tt, 1, i2mex_o%id_z, nt1*ns, z, iok)
    if(iok /= 0) ier = 31
    
    !!call i2mex_toAxis(nt1, ns, the, psi,  0.0_i2mex_r8, z, iok)
    !!call i2mex_error(iok)

    deallocate(tt)
    deallocate(ss)
    
  end subroutine i2mex_getZ

  subroutine i2mex_getGradZ(nt1, ns, the, psi, zt, zp, ier)

    ! Get (d Z/d the, d Z/ dpsi) on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: zt(nt1, ns), zp(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns), dsdpsi(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    real(i2mex_r8), dimension(:,:,:), allocatable :: dz
    integer iok, i


    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    allocate(dz(nt1, ns, 2), stat=iok)
    if(iok/=0) ier = 37

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_grhochi(nt1*ns,ss,tt,1,i2mex_o%id_z,nt1*ns,dz,iok)
    if(iok /= 0) ier = 32
        
    call i2mex_getDS(ns, psi, dsdpsi, iok)
    call i2mex_error(iok)

    zt(1:nt1, 1:ns) = dz(1:nt1, 1:ns,2)
    do i=1, nt1
       zp(i, 1:ns) = dz(i,1:ns,1)*dsdpsi(1:ns)
    enddo

    call i2mex_toAxis(nt1, ns, the, psi, +0.5_i2mex_r8, zt, iok)    
    call i2mex_error(iok)
    call i2mex_toAxis(nt1, ns, the, psi, -0.5_i2mex_r8, zp, iok)    
    call i2mex_error(iok)
    
    deallocate(tt)
    deallocate(ss)
    deallocate(dz)
    
  end subroutine i2mex_getGradZ

  subroutine i2mex_getZtt(nt1, ns, the, psi, ztt, ier)

    ! Get d^2 Z/d the^2 on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: ztt(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    integer iok

    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    if(iok/=0) ier = 38

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_hrhochi(nt1*ns,ss,tt,1,i2mex_o%id_z,                    &
         &                (/0,0,0,0,1,0/),nt1*ns,ztt,iok)
    if(iok /= 0) ier = 33
    
    call i2mex_toAxis(nt1, ns, the, psi,  1.0_i2mex_r8, ztt, iok)    
    call i2mex_error(iok)
    
    deallocate(tt)
    deallocate(ss)
    
  end subroutine i2mex_getZtt

  subroutine i2mex_getZtp(nt1, ns, the, psi, zpt, ier)

    ! Get d Z/(d psi d the) on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: zpt(nt1, ns)
    integer, intent(out) :: ier

    call i2mex_getZpt(nt1, ns, the, psi, zpt, ier)

  end subroutine i2mex_getZtp

  subroutine i2mex_getZpt(nt1, ns, the, psi, zpt, ier)

    ! Get d Z/(d psi d the) on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: zpt(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns), dsdpsi(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    integer iok, i

    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    if(iok/=0) ier = 39

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_hrhochi(nt1*ns,ss,tt,1,i2mex_o%id_z,                    &
         &                (/0,0,0,0,0,1/),nt1*ns,zpt,iok)
    if(iok /= 0) ier = 34

    call i2mex_getDS(ns, psi, dsdpsi, iok)
    call i2mex_error(iok)

    do i=1, nt1
       zpt(i, 1:ns) = zpt(i,1:ns)*dsdpsi(1:ns)
    enddo
    
    call i2mex_toAxis(nt1, ns, the, psi, 0.0_i2mex_r8, zpt, iok)    
    call i2mex_error(iok)
   
    deallocate(tt)
    deallocate(ss)
    
  end subroutine i2mex_getZpt


  subroutine i2mex_getZpp(nt1, ns, the, psi, zpp, ier)

    ! Get d^2 Z/d psi^2 on the new (the, psi) mesh

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(out) :: zpp(nt1, ns)
    integer, intent(out) :: ier

    real(i2mex_r8) :: s(ns), dsdpsi(ns), d2sdpsi2(ns)
    real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
    real(i2mex_r8), dimension(:,:,:), allocatable :: zs
    integer iok, i

    ier=0
    call i2mex_getS(ns, psi, s, iok)
    call i2mex_error(iok)

    allocate(tt(nt1, ns), stat=iok)
    allocate(ss(nt1, ns), stat=iok)
    allocate(zs(nt1, ns, 2), stat=iok)
    if(iok/=0) ier = 40

    ss = spread(  s, dim=1, ncopies=nt1)
    tt = spread(the, dim=2, ncopies=ns )

    call eq_hrhochi(nt1*ns,ss,tt,1,i2mex_o%id_z,                    &
         &                (/0,1,0,1,0,0/),nt1*ns,zs,iok)
    if(iok /= 0) ier = 35

    call i2mex_getDS(ns, psi, dsdpsi, iok)
    call i2mex_error(iok)
    call i2mex_getD2S(ns, psi, d2sdpsi2, iok)
    call i2mex_error(iok)
    do i=1, nt1
       zpp(i, 1:ns) = zs(i,1:ns,2)*dsdpsi(1:ns)**2 + &
            & zs(i,1:ns,1)*d2sdpsi2(1:ns)
    enddo
    
    call i2mex_toAxis(nt1, ns, the, psi, -1.0_i2mex_r8, zpp, iok)    
    call i2mex_error(iok)
    
    deallocate(tt)
    deallocate(ss)
    deallocate(zs)
    
  end subroutine i2mex_getZpp
