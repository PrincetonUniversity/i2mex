subroutine i2mex_get_iBoozer(nt1, ns, the, psi, fun, ier)

  ! Get Iboozer on the (the, psi) mesh
  !   I * (Lp/dtheta)**-1 gives contribution to
  !   Bpol; see also i2mex_get_nuBoozer

  !     mod(Bpol) = (dLp/dtheta)**-1 * (g(psi)*d(nu)/dtheta + I)
  !     dLp/dtheta = sqrt((dR/dtheta)**2+(dZ/dtheta)**2)

  !   I is d/dtheta[ id_iboozer profile in xplasma]

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
  real(i2mex_r8), intent(out) :: fun(nt1, ns)   ! Iboozer returned here
  integer, intent(out) :: ier

  !---------------------------
  real(i2mex_r8) :: s(ns)
  real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
  integer iok,ict_dfdth(6)
  integer :: id_iboozer, id_nuboozer
  !---------------------------

  call eqm_gen_boozer(id_iboozer,id_nuboozer,ier)
  if(ier.ne.0) then
     ier = 171
     return
  endif

  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)

  !  set up to get 2d profile

  allocate(tt(nt1, ns), stat=iok)
  allocate(ss(nt1, ns), stat=iok)
  if(iok/=0) ier = 36

  ss = spread(  s, dim=1, ncopies=nt1)
  tt = spread(the, dim=2, ncopies=ns )

  call eq_frhochi(nt1*ns, ss, tt, 1, id_iboozer, nt1*ns, fun, iok)
  if(iok /= 0) then
     ier = 171
     return
  endif
    
  deallocate(tt)
  deallocate(ss)
    
end subroutine i2mex_get_iBoozer

subroutine i2mex_get_nuBoozer(nt1, ns, the, psi, fun, ier)

  ! Get nu(theta,psi) on the (the, psi) mesh
  !   Boozer toroidal coordinate:  zeta = phi + nu

  !   d(nu(theta,psi))/dtheta * (Lp/dtheta)**-1 gives contribution to
  !   Bpol; see also i2mex_get_IBoozer

  !     mod(Bpol) = (dLp/dtheta)**-1 * (g(psi)*d(nu)/dtheta + d(I)/dtheta)
  !     dLp/dtheta = sqrt((dR/dtheta)**2+(dZ/dtheta)**2)

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
  real(i2mex_r8), intent(out) :: fun(nt1, ns)   ! Iboozer returned here
  integer, intent(out) :: ier

  !---------------------------
  real(i2mex_r8) :: s(ns),numin,numax,nu0,nu1
  real(i2mex_r8), dimension(:,:), allocatable :: ss, tt
  integer iok
  integer :: id_iboozer, id_nuboozer
  integer :: is,ith
  real(i2mex_r8), parameter :: zsmall = 1.0d-8
  !---------------------------

  call eqm_gen_boozer(id_iboozer,id_nuboozer,ier)
  if(ier.ne.0) then
     ier = 172
     return
  endif

  call i2mex_getS(ns, psi, s, iok)
  call i2mex_error(iok)

  !  set up to get 2d profile

  allocate(tt(nt1, ns), stat=iok)
  allocate(ss(nt1, ns), stat=iok)
  if(iok/=0) ier = 36

  ss = spread(  s, dim=1, ncopies=nt1)
  tt = spread(the, dim=2, ncopies=ns )

  call eq_frhochi(nt1*ns, ss, tt, 1, id_nuboozer, nt1*ns, fun, iok)
  if(iok /= 0) then
     ier = 171
     return
  endif
    
  deallocate(tt)
  deallocate(ss)

  !  fun is periodic; make it periodic w/endpoints = 0

  do is=1,ns
     nu0=fun(1,is)
     nu1=fun(nt1,is)
     numin=min(nu0,nu1)
     numax=max(nu0,nu1)
     fun(1,is)=0
     fun(nt1,is)=0
     do ith=2,nt1-1
        numin=min(numin,fun(ith,is))
        numax=max(numax,fun(ith,is))
        fun(ith,is)=fun(ith,is)-nu0
     enddo
     
     if(abs(nu1-nu0).gt.zsmall*(numax-numin)) then
        ier = 175
     endif
  enddo

end subroutine i2mex_get_nuBoozer
