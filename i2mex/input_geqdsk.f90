subroutine i2mex_fromGeqdsk(filename, it_orientation, ier)
 
  ! read G-EQDSK file and and perform inverse map
 
   use i2mex_mod
   use geqdsk_mod
   use ezspline_obj
   use ezspline
   use cont_mod
 
 
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
  type(geqdsk) :: geq
  type(ezspline2) :: pspl
 
  real(r8), dimension(:), allocatable :: p, q, g, psi, psi_geq, psis, &
       & the, xgrid, zgrid, tlike, q_est
  real(r8), dimension(:, :), allocatable :: x, z
  real(r8), dimension(:, :), allocatable :: xin, zin
  real(r8) :: pmax, pmin, xmin, xmax, zmin, zmax, surf, xm, zm, xp, zp, &
       & pmag, psep, ttot
  integer iok, i, j, ioutboard(1)
 
  real(r8), dimension(:,:), allocatable :: integrand
  real(r8), dimension(:), allocatable :: qres
  real(r8) x0, z0, x1, z1, x2, z2, det
  real(r8) ztol, delta, dlim
 
  ier = 0
 
  call geq_init(geq, filename, iok)
  call geq_error(iok)
  if(iok/=0) then
     ier = 117
     return
  endif
 
  xmin = geq%RLEFT_M
  xmax = geq%RLEFT_M + geq%RDIM_M
  zmin = geq%ZMID_M - geq%ZDIM_M/2.0_r8
  zmax = geq%ZMID_M + geq%ZDIM_M/2.0_r8
 
 
  allocate(xgrid(geq%nw))
  allocate(zgrid(geq%nh))
 
  xgrid=xmin + (xmax - xmin)*(/ (real(i-1,r8)/real(geq%NW-1,r8), i=1,geq%NW) /)
  zgrid=zmin + (zmax - zmin)*(/ (real(i-1,r8)/real(geq%NH-1,r8), i=1,geq%NH) /)
 
!!$  i2mex_nt1 = 257
!!$  i2mex_ns = 401
 
  allocate(psi(i2mex_ns), psi_geq(geq%NW), p(i2mex_ns), g(i2mex_ns), q(i2mex_ns), the(i2mex_nt1))
  allocate(tlike(i2mex_ns), q_est(i2mex_ns))
  allocate(x(i2mex_nt1, i2mex_ns), z(i2mex_nt1, i2mex_ns))
  allocate(xin(i2mex_nt1, i2mex_ns), zin(i2mex_nt1, i2mex_ns))
  allocate(psis(i2mex_nt1))
 
  xm = geq%RMAXIS_M
  zm = geq%ZMAXIS_M
  ioutboard = maxloc(geq%RBBBS_M) ! most outboard point
  xp = geq%RBBBS_M(ioutboard(1))
  zp = geq%ZBBBS_M(ioutboard(1))
 
  call ezspline_init(pspl, geq%NW, geq%NH, (/0,0/), (/0,0/), iok)
  call ezspline_error(iok)
 
  ! toroidal current sign
  cont_sign = 1.0_r8
  if(geq%SIBRY_Wb__Rad < geq%SIMAG_Wb__Rad) cont_sign = -1.0_r8
 
  ! We like to have psi = 0 on magnetic axis, and increasing
  ! outwards (this will change the current orientation but that's
  ! ok).
  geq%psirz_Wb__Rad = cont_sign*(geq%psirz_Wb__Rad - geq%SIMAG_Wb__Rad)
  geq%SIBRY_Wb__Rad = cont_sign*(geq%SIBRY_Wb__Rad - geq%SIMAG_Wb__Rad)
  geq%SIMAG_Wb__Rad = 0.0_r8
  geq%FFPRIM_T2M2Rad__Wb = cont_sign*geq%FFPRIM_T2M2Rad__Wb
  geq%pprime_NtRad__M2Wb = cont_sign*geq%pprime_NtRad__M2Wb
  geq%CURRENT_A = cont_sign*geq%CURRENT_A
  ! set back
  cont_sign = 1.0_r8
 
  pspl%x1 = xgrid
  pspl%x2 = zgrid
  !pspl%isHermite = 1 ! Akima Hermite
  call ezspline_setup(pspl, geq%PSIRZ_Wb__Rad, iok)
  call ezspline_error(iok)
 
  ! interpolate to get psi on axis
  call ezspline_interp(pspl, xm, zm, pmag, iok)
  call ezspline_error(iok)
  ! psi on separatrix
  call ezspline_interp(pspl, xp, zp, psep, iok)
  call ezspline_error(iok)
 
  do j = 2, i2mex_ns
     ! rough estimate for q (used to determine angle integration step)
     cont_q = geq%qpsi(1) + &
          & real(j,r8)*(geq%qpsi(geq%NW) - geq%qpsi(1))/real(i2mex_ns,r8)
 
     surf = real(j-1, r8)/real(i2mex_ns-1, r8)
     xp = xm + surf * i2mex_LAST_NORM_SURFACE_IS*(geq%RBBBS_M(ioutboard(1)) - xm)
     zp = zm + surf * i2mex_LAST_NORM_SURFACE_IS*(geq%ZBBBS_M(ioutboard(1)) - zm)
     call i2mex_contour(geq%NW, geq%NH, xgrid, zgrid, geq%PSIRZ_Wb__Rad,  &
       & xm, zm, xp, zp, i2mex_nt1, ttot, xin(1, j), zin(1, j), iok)
     call i2mex_error(iok)
 
     ! get psi on contour
     call ezspline_interp(pspl, i2mex_nt1, xin(:, j), zin(:, j), psis, iok)
     call ezspline_error(iok)
 
     ! average psi on contour
     psi(j) = sum(psis)/real(i2mex_nt1, r8)
 
     ! save tot time
     tlike(j) = ttot
 
  enddo
  if(iok/=0) ier = 120
  psi(1) = pmag
 
  xin(1:i2mex_nt1, 1) = xm
  zin(1:i2mex_nt1, 1) = zm
 
  ! normalized to psi on magnetic axis = 0
  psi = psi - psi(1)
  pmax = MAXVAL(psi)
  pmin = MINVAL(psi)
 
  psi_geq = 0._r8 + ((pmax-pmin)/i2mex_LAST_NORM_SURFACE_IS)* &
       & (/ (real(i-1, r8)/real(geq%nw-1,r8), i=1, geq%nw) /)
 
 
  call ezspline_free(pspl, iok)
 
  the = i2mex_twopi_r8* (/ (real(i-1,r8)/real(i2mex_nt1-1,r8), i=1, i2mex_nt1) /)
 
  ! pressure
  call i2mex_interp1d(geq%NW, psi_geq, &
       & geq%PRES_Nt__M2 * 4.e-7_r8 * i2mex_pi_r8, &
       & i2mex_ns, psi, p)
  ! q
  call i2mex_interp1d(geq%NW, psi_geq, &
       & geq%QPSI                                , &
       & i2mex_ns, psi, q)
  ! g
  call i2mex_interp1d(geq%NW, psi_geq, &
       & abs(geq%FPOL_TM)                        , &
       & i2mex_ns, psi, g)
 
  ! estimate for q
  q_est(2:i2mex_ns) = abs(g(2:i2mex_ns))*tlike(2:i2mex_ns)/(i2mex_twopi_r8)
  q_est(1) = q_est(2) + (q_est(3)-q_est(2))*(psi(1)-psi(2))/(psi(3)-psi(2))
 
  q = q_est ! take our own q
 
  ! initialize i2mex
  ! get theta orientation
  x0 = (maxval(xin(:,1))+minval(xin(:,1)))/2._i2mex_r8
  z0 = (maxval(zin(:,1))+minval(zin(:,1)))/2._i2mex_r8
  x1 = xin(1, i2mex_ns)
  z1 = zin(1, i2mex_ns)
  x2 = xin(2, i2mex_ns)
  z2 = zin(2, i2mex_ns)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
 
  if( it_orientation*det < 0.0_i2mex_r8 ) then
     print *,'keep theta orientation'
     do i = 1, i2mex_nt1
        x(i, 1:i2mex_ns) = xin(i, 1:i2mex_ns)
        z(i, 1:i2mex_ns) = zin(i, 1:i2mex_ns)
     enddo
  else
     print *,'reverse theta orientation'
     do i = 1, i2mex_nt1
        x(i, 1:i2mex_ns) = REAL(xin(i2mex_nt1-i+1, 1:i2mex_ns), i2mex_r8)
        z(i, 1:i2mex_ns) = REAL(zin(i2mex_nt1-i+1, 1:i2mex_ns), i2mex_r8)
     enddo
  endif
 
  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & i2mex_nt1, i2mex_ns, the, psi, g, x, z, iok)
  endif
  
  call i2mex_init('GEQDSK_'//filename, &
       i2mex_nt1, i2mex_ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  if(i2mex_direct>0) then
        ! Build Psi(R,Z) representation, if required...
        ztol  = 1.0e-6*xmax
        dlim  = 0.05_r8 * Xmax ! box boundaries
        delta = 0.0_r8         ! smoothing parameter

        i2mex_o%Rleft = xmin
        i2mex_o%Rrigh = xmax
        i2mex_o%Zbot = zmin
        i2mex_o%Ztop = zmax

        call eqm_cbdy(5, &
             & (/xmin, xmax, xmax, xmin, xmin /), &
             & (/zmin, zmin, zmax, zmax, zmin /), &
             & iok)
        if(iok/=0) ier = 154

        call eqm_rzgrid(xgrid, zgrid, -1, -1, geq%NW, geq%NH, ztol, &
             & i2mex_o%id_Rgrid, i2mex_o%id_Zgrid, iok)
        if(iok/=0) then
           ier = 151
           call i2mex_error(ier)
        endif

        ! interpolation order is Akima
        call eqm_rzfunda('PSIRZ', i2mex_o%id_psirz, geq%psirz_Wb__Rad, &
             & geq%NW, geq%NH, i2mex_akima, delta, iok)
        if(iok/=0) ier = 153
  endif
 
  deallocate(xgrid)
  deallocate(zgrid)
  deallocate(psi, psi_geq, p, g, q, the)
  deallocate(tlike, q_est)
  deallocate(x, z)
  deallocate(xin, zin)
  deallocate(psis)
 
  call geq_free(geq, iok)
  if(iok/=0) ier = 150
 
 
end subroutine i2mex_fromGeqdsk
