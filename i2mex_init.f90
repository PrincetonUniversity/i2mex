subroutine i2mex_init(label, nt1, ns, t, psi, p, g, q, x, z, ier)
 
  ! Initialization and memory allocation. Construct the set of primitive
  ! equilibrium data, given a plasma equilibrium in terms of
  ! the poloidal flux [Wb/Rad], the pressure [Tesla^2], the covariant
  ! toroidal field profile g [T*m], the safety factor profile q and a
  ! set of ns nested flux surfaces (x, z) with nt1 poloidal points.
  !
  ! Under some circumstances, the input data 
  !
  ! ier: Error flag (=0 is ok)
 
  use i2mex_mod
  implicit none
  character*(*), intent(in) :: label
  integer, intent(inout) :: nt1, ns
  real(i2mex_r8), intent(inout) :: t(nt1)
  real(i2mex_r8), intent(inout) :: psi(ns), p(ns), g(ns), q(ns)
  real(i2mex_r8), intent(inout) :: x(nt1, ns), z(nt1, ns)
  integer, intent(out) :: ier
 
  real(i2mex_r8) s(ns), dum1, dum2, smin_eps, x0, z0, x1, z1, x2, z2, det
  real(i2mex_r8) ajac_maxvar
  integer iok
 
  real(i2mex_r8), dimension(nt1) ::   zbcR0, zbcR1,  zbcZ0, zbcZ1

  real(i2mex_r8), dimension(:,:), allocatable :: temp
  integer it, i

  real(i2mex_r8), dimension(:,:,:), allocatable :: dx, dz
  real(i2mex_r8), dimension(:,:), allocatable :: ss, tt, jakima, dsdp

  integer nt1_3, ns_3
  real(i2mex_r8), dimension(:), allocatable :: the_3, psi_3, p_3, g_3, q_3
  real(i2mex_r8), dimension(:,:), allocatable :: x_3, z_3
  
  ier = 0

  ! Decimate the data if necessary
  if(((i2mex_nt1max > 0) .or. (i2mex_nsmax > 0)) .and. &
       & ((ns > i2mex_nsmax) .or. (nt1 > i2mex_nt1max)) ) then

     nt1_3 = min(nt1, max(i2mex_nt1max, 17))
     ns_3  = min(ns,  max(i2mex_nsmax , 21))
     
     
     allocate(the_3(nt1_3), psi_3(ns_3))
     allocate(p_3(ns_3), g_3(ns_3), q_3(ns_3))
     allocate(x_3(nt1_3, ns_3), z_3(nt1_3, ns_3))

     psi_3 = psi(1) + (psi(ns)-psi(1))* (/  ( &
          & real(i-1,i2mex_r8)/real(ns_3-1,i2mex_r8), &
          & i=1,ns_3 ) /)

     the_3 = t(1) + (t(nt1)-t(1))* (/  ( &
          & real(i-1,i2mex_r8)/real(nt1_3-1,i2mex_r8), &
          & i=1,nt1_3 ) /)

     call i2mex_decimate(nt1, ns, t, psi, p, g, q, x, z, &
          &  nt1_3, ns_3, the_3, psi_3, p_3, g_3, q_3, x_3, z_3, &
          &  iok)
     call i2mex_error(iok)

     ! replace
     t(1:nt1_3)  = the_3
     psi(1:ns_3) = psi_3
     p(1:ns_3) = p_3
     q(1:ns_3) = q_3
     g(1:ns_3) = g_3
     x(1:nt1_3, 1:ns_3) = x_3
     z(1:nt1_3, 1:ns_3) = z_3

     ns = ns_3
     nt1 = nt1_3

     deallocate(the_3, psi_3)
     deallocate(p_3, g_3, q_3)
     deallocate(x_3, z_3)

  end if
     
     
 
  i2mex_o%nt1 = nt1
  i2mex_o%ns  = ns
  i2mex_o%psi_edge = psi(ns)
  i2mex_o%psi_axis = psi(1)
 
 
  ! Set up the mesh. Take radial variable ~ sqrt(psi) to avoid singular
  ! Jacobian behaviour near axis
 
  s = sqrt( (psi - psi(1))/(psi(ns) - psi(1)) )
 
  i2mex_o%label = label
  call eqm_select('i2mex-'//label, i2mex_torSym)

  call eqm_rho(s, ns, i2mex_tol, i2mex_o%id_s, iok)
  if(iok /=0) ier = 1
 
  call eqm_chi(t, i2mex_No_autoT, nt1, i2mex_tol, i2mex_o%id_t, iok)
  if(iok /=0) ier = 2
 
 
  call eqm_rhofun(i2mex_cubic, i2mex_o%id_s, 'PSI', psi, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_psi, iok)
  if(iok /=0) ier = 3
 
  call eqm_rhofun(i2mex_cubic, i2mex_o%id_s, 'G', g, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_g, iok)
  if(iok /=0) ier = 5
 
  call eqm_rhofun(i2mex_cubic, i2mex_o%id_s, 'P', p, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_p, iok)
!!$  ! pressure can be discontinuous in slope => use Akima
!!$  call eqm_rhofun(i2mex_akima, i2mex_o%id_s, 'P', p, i2mex_not_a_knot, &
!!$       & dum1, i2mex_not_a_knot, dum2, &
!!$       & i2mex_o%id_p, iok)
  if(iok /=0) then
     ier = 4
  else
     call eqm_mark_pmhd(i2mex_o%id_p,iok)
     if(iok /=0) ier = 4
  endif
 
  call eqm_rhofun(i2mex_cubic, i2mex_o%id_s, 'Q', q, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_q, iok)
  if(iok /=0) ier = 6
 
!!$   call eqm_rzmag(x, z, nt1, ns, 2, i2mex_o%id_x, i2mex_o%id_z, iok)
!!$   if(iok /=0) ier = 7
 
  ! determine if theta is clockwise
  i2mex_o%isThetaClockwise = 0
  x0 = sum(x(1:nt1-1,1))/real(nt1-1)
  z0 = sum(z(1:nt1-1,1))/real(nt1-1)
  x1 = x(1, ns)
  z1 = z(1, ns)
  x2 = x(2, ns)
  z2 = z(2, ns)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
  if(det < 0.0_i2mex_r8) then
     i2mex_o%isThetaClockwise = 1
  endif
 

  ! phase shift 
  if(i2mex_it0 > 0 .and. i2mex_it0 < nt1) then 
     allocate(temp(nt1,ns))
     i = i2mex_it0 !(nt1-1)/4
     print *,'i2mex::WARNING theta origin shifted by ', i, ' nt1=', nt1
     temp(1      :nt1-1-i  ,:)   = x(i+1:nt1-1,:)
     temp(nt1-i  :nt1-1    ,:)   = x(1:i      ,:)
     temp(nt1,:) = temp(1,:)
     x(1:nt1,1:ns) = temp
     temp(1      :nt1-1-i  ,:)   = z(i+1:nt1-1,:)
     temp(nt1-i  :nt1-1    ,:)   = z(1:i      ,:)
     temp(nt1,:) = temp(1,:)
     z(1:nt1,1:ns) = temp
     deallocate(temp)
  endif

  ! load x and z
 
!!$   call eqm_rzmagb(x, z, nt1, ns, 2,  &
!!$     & 0, zbcR0, 0, zbcR1,  &
!!$     & 0, zbcZ0, 0, zbcZ1,  &
!!$     & i2mex_o%id_x, i2mex_o%id_z ,iok)
   call eqm_rzmag(x, z, nt1, ns, 2,  &
     & i2mex_o%id_x, i2mex_o%id_z ,iok)
   if(iok /=0) ier = 7
 
   ! not setting the B direction will prevent the calculation
   ! of Br, Bz, and Bmod.
   if(i2mex_direct > 0 ) call eqm_bset(i2mex_bphi_cw, i2mex_Ip_cw)
 
   smin_eps = s(1) + 1.e-2_i2mex_r8 *(s(2) - s(1))
   call eq_flxint_init(0, ns, s, smin_eps, iok)
   if(iok /=0) ier = 8

   ! Set boundary box coordinates
   i2mex_o%Rleft = minval(x)
   i2mex_o%Rrigh = maxval(x)
   i2mex_o%Zbot = minval(z)
   i2mex_o%Ztop = maxval(z)


   ! the following allows one to build akima splines for the grid 
   ! coordinates, this in order to prevent the ringing problem that
   ! arises for noisy data.

   if(i2mex_xzakima > 0) then

      call eqm_frhochi(i2mex_akima, i2mex_o%id_t, i2mex_o%id_s, 'xakima', x, nt1, &
           & i2mex_not_a_knot, dum1, &
           & i2mex_not_a_knot, dum2, &
           & i2mex_o%id_xakima, iok)
      if(iok /=0) ier = 7

      call eqm_frhochi(i2mex_akima, i2mex_o%id_t, i2mex_o%id_s, 'zakima', z, nt1, &
           & i2mex_not_a_knot, dum1, &
           & i2mex_not_a_knot, dum2, &
           & i2mex_o%id_zakima, iok)
      if(iok /=0) ier = 8

      ! compute the Jacobian and store
      allocate(dx(nt1, ns, 2))
      allocate(dz(nt1, ns, 2))
      allocate(ss(nt1,ns), tt(nt1,ns), jakima(nt1,ns), dsdp(nt1,ns))
      ss = spread(  s, dim=1, ncopies=nt1)
      tt = spread(  t, dim=2, ncopies=ns )
      call eq_grhochi(nt1*ns,ss,tt,1,i2mex_o%id_xakima,nt1*ns,dx,iok)
      call eq_grhochi(nt1*ns,ss,tt,1,i2mex_o%id_zakima,nt1*ns,dz,iok)
      s(1) = 1.e-8_i2mex_r8
      dsdp = spread( 0.5_i2mex_r8 /(s*(psi(ns) - psi(1))), dim=1, ncopies=nt1)
      jakima = abs((dx(:,:,2)*dz(:,:,1) - dx(:,:,1)*dz(:,:,2))*x*dsdp)
      call i2mex_toAxis(nt1, ns, t, psi,  0.0_i2mex_r8, jakima, iok)    
      call eqm_frhochi(i2mex_akima, i2mex_o%id_t, i2mex_o%id_s, &
           & 'jakima', jakima, nt1, &
           & i2mex_not_a_knot, dum1, &
           & i2mex_not_a_knot, dum2, &
           & i2mex_o%id_jakima, iok)
      deallocate(dx)
      deallocate(dz)
      deallocate(ss, tt, jakima, dsdp)

   endif
  
  
  


 
!!$  call eqm_bset(i2mex_bphi_ccw, i2mex_Ip_ccw)
 
 
end subroutine i2mex_init
