subroutine i2mex_ToMapdsk(mth1, nosurf, the, psi, filename, ier)
  
  ! save data in MAPDSK.cdf file 

  use ezcdf
  use i2mex_mod
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)

  integer, intent(in) :: mth1 ! no of poloidal rays
  integer, intent(in) :: nosurf ! no of radial grid points
  real(r8), intent(in) :: the(mth1) ! theta grid
  real(r8), intent(in) :: psi(nosurf) ! flux grid
  character*(*), intent(in) :: filename ! data will be saved in 'filename'
  integer, intent(out) :: ier ! returned error flag (0=ok)

  integer mth2, mth, mth1x, i, j, iu, dims(3), iok
  character*8 today
  character*10 now
  character*80 dat
  real(r8) :: xma, zma, x0, b0, ip, ppf, beta, eli, res
  real(r8) :: xmin, xmax, zmin, zmax, r0, z0, xin, xout, zbot, ztop
  real(r8), dimension(:,:), allocatable :: x, xt, xp, xtt, xtp, xpp
  real(r8), dimension(:,:), allocatable :: z, zt, zp, ztt, ztp, zpp
  real(r8), dimension(:,:), allocatable :: jac, jacp
  real(r8), dimension(:,:), allocatable :: gtt, gtp, gpp
  real(r8), dimension(:,:), allocatable :: temp2
  real(r8), dimension(:,:), allocatable :: xsq, xsqdps

  real(r8), allocatable, dimension(:) :: temp1
  real(r8), allocatable, dimension(:) :: tempt

  ! Glasser, Green & Johnson coefficients
  real(r8), dimension(nosurf) :: e, f, h

!!$  real(r8), dimension(nosurf) :: v_of_psi, vp_of_psi
!!$  real(r8), dimension(nosurf) :: g, gp, q, qp, qpp, p, pp
  real(r8), dimension(:), allocatable :: v_of_psi, vp_of_psi
  real(r8), dimension(:), allocatable :: g, gp, q, qp, qpp, p, pp
 
 ier = 0
 
  mth = mth1 - 1
  mth2= mth1 + 1
  mth1x = mth1

  call cdfOpn(iu, filename, 'w')
  if(iu==0) then
     print*,'--MAPDSK-- failed to open ', filename
     ier = 1
     return
  endif
  
  dims=(/80, 1, 1/)
  call cdfDefVar(iu, 'title', dims, 'CHAR')
  call cdfDefVar(iu, 'date', dims, 'CHAR')
  dims=(/132, 1, 1/)
  call cdfDefVar(iu, 'comments', dims, 'CHAR')
  dims=(/1, 1, 1/)
  call cdfDefVar(iu, 'mth', dims, 'INT')
  call cdfDefVar(iu, 'nosurf', dims, 'INT')

  call cdfDefVar(iu, 'remap', dims, 'INT')
  call cdfDefVar(iu, 'mx', dims, 'INT')
  call cdfDefVar(iu, 'npsi', dims, 'INT')
  call cdfDefVar(iu, 'kb', dims, 'INT')

  call cdfDefVar(iu, 'Xmin', dims, 'R8')
  call cdfDefVar(iu, 'Xmax', dims, 'R8')
  call cdfDefVar(iu, 'X0', dims, 'R8')
  call cdfDefVar(iu, 'Xmag', dims, 'R8')

  call cdfDefVar(iu, 'Zmin', dims, 'R8')
  call cdfDefVar(iu, 'Zmax', dims, 'R8')
  call cdfDefVar(iu, 'Z0', dims, 'R8')
  call cdfDefVar(iu, 'Zmag', dims, 'R8')

  call cdfDefVar(iu, 'B0', dims, 'R8')
  call cdfDefVar(iu, 'Ip', dims, 'R8')
  call cdfDefVar(iu, 'Beta', dims, 'R8')
  call cdfDefVar(iu, 'BetaStar', dims, 'R8')
  call cdfDefVar(iu, 'BetaN', dims, 'R8')
  call cdfDefVar(iu, 'li', dims, 'R8')
  call cdfDefVar(iu, 'PPF', dims, 'R8')
  call cdfDefVar(iu, 'Upsiln', dims, 'R8')
  dims=(/nosurf, 1, 1/)
  call cdfDefVar(iu, 'JdotBOverBSquare', dims, 'R8')
  call cdfDefVar(iu, 'ballooning_alpha', dims, 'R8')
  call cdfDefVar(iu, 'local_shear_s', dims, 'R8')
  call cdfDefVar(iu, 'Di', dims, 'R8')
  call cdfDefVar(iu, 'Dr', dims, 'R8')
  call cdfDefVar(iu, 'surface_averaged_radius', dims, 'R8')
  call cdfDefVar(iu, 'PsiBig', dims, 'R8')
  call cdfDefVar(iu, 'V_of_psi', dims, 'R8')
  call cdfDefVar(iu, 'Vp_of_psi', dims, 'R8')
  call cdfDefVar(iu, 'pa', dims, 'R8')
  call cdfDefVar(iu, 'ppa', dims, 'R8')
  call cdfDefVar(iu, 'qa', dims, 'R8')
  call cdfDefVar(iu, 'qpa', dims, 'R8')
  call cdfDefVar(iu, 'qppa', dims, 'R8')
  call cdfDefVar(iu, 'ga', dims, 'R8')
  call cdfDefVar(iu, 'gpa', dims, 'R8')
  call cdfDefVar(iu, 'spest1', dims, 'R8')
!!$  call cdfDefVar(iu, 'dpsi_drho_fa', dims, 'R8')
!!$  call cdfDefVar(iu, 'dpsi2_drho2_fa', dims, 'R8')
  dims=(/mth1x, 1, 1/)
  call cdfDefVar(iu, 'xinf', dims, 'R8')
  call cdfDefVar(iu, 'zinf', dims, 'R8')
  dims=(/mth1x, nosurf, 1/)
  call cdfDefVar(iu, 'xa', dims, 'R8')
  call cdfDefVar(iu, 'za', dims, 'R8')
  call cdfDefVar(iu, 'xpth', dims, 'R8')
  call cdfDefVar(iu, 'zpth', dims, 'R8')
  call cdfDefVar(iu, 'xpsi', dims, 'R8')
  call cdfDefVar(iu, 'zpsi', dims, 'R8')
  call cdfDefVar(iu, 'grpssq', dims, 'R8')
  call cdfDefVar(iu, 'grthsq', dims, 'R8')
  call cdfDefVar(iu, 'grpsth', dims, 'R8')
  call cdfDefVar(iu, 'gptdth', dims, 'R8')
!!$  call cdfDefVar(iu, 'xsqdps', dims, 'R8')
!!$  call cdfDefVar(iu, 'xsqdth', dims, 'R8')
  call cdfDefVar(iu, 'gpsdth', dims, 'R8')
  call cdfDefVar(iu, 'xjacob', dims, 'R8')
  call cdfDefVar(iu, 'xjprym', dims, 'R8')
  call cdfDefVar(iu, 'delta', dims, 'R8')
  call cdfDefVar(iu, 'qdelp', dims, 'R8')
  call cdfDefVar(iu, 'qdelt', dims, 'R8')
  call cdfDefVar(iu, 'gserror', dims, 'R8')

  ALLOCATE( temp1(nosurf) )
  ALLOCATE( tempt(mth1) )
  allocate(v_of_psi(nosurf), vp_of_psi(nosurf))
  allocate(g(nosurf), gp(nosurf), q(nosurf), qp(nosurf), qpp(nosurf), &
       & p(nosurf), pp(nosurf))
  

  
  call cdfPutVar(iu, 'title', i2mex_o%label, iok)
  if(iok/=0) then
     print *,'--MAPDSK--ERROR failed to write title=',i2mex_o%label
     ier = 2
  endif
  call cdfPutVar(iu, 'comments', i2mex_comments, iok)
  if(iok/=0) then
     print *,'--MAPDSK--ERROR failed to write comments=',i2mex_comments
     ier = 2
  endif
  call date_and_time(date=today, time=now)
  dat = today(1:4)//'/'//today(5:6)//'/'// &
       & today(7:8)//' at '//now(1:2)//':'//now(3:4)
  call cdfPutVar(iu, 'date', dat, iok)
  if(iok/=0) then
     print *,'--MAPDSK--ERROR failed to write date=',dat
     ier = 2
  endif

  call cdfPutvar(iu, 'mth', mth, iok)
  call cdfPutvar(iu, 'nosurf', nosurf, iok)

  call cdfPutvar(iu, 'remap', i2mex_remap, iok)
  call cdfPutvar(iu, 'mx', i2mex_mx, iok)
  call cdfPutvar(iu, 'npsi', i2mex_npsi, iok)
  call cdfPutvar(iu, 'kb', i2mex_kb, iok)

  allocate(temp2(mth1x, nosurf))
  allocate(x(mth1x, nosurf), &
       & xt(mth1x, nosurf), xp(mth1x, nosurf), &
       & xtt(mth1x, nosurf), xtp(mth1x, nosurf), xpp(mth1x, nosurf))
  allocate(z(mth1x, nosurf), &
       & zt(mth1x, nosurf), zp(mth1x, nosurf), &
       & ztt(mth1x, nosurf), ztp(mth1x, nosurf), zpp(mth1x, nosurf))


  call i2mex_getX(mth1, nosurf, the, psi, x, ier)
  call i2mex_getGradX(mth1, nosurf, the, psi, xt, xp, ier)
  call i2mex_getXtt(mth1, nosurf, the, psi, xtt, ier)
  call i2mex_getXtp(mth1, nosurf, the, psi, xtp, ier)
  call i2mex_getXpp(mth1, nosurf, the, psi, xpp, ier)
  call i2mex_error(ier)
!!$  x(mth1x,:) = x(2,:)
!!$  xt(mth1x,:) = xt(2,:)
!!$  xp(mth1x,:) = xp(2,:)
!!$  xtt(mth1x,:) = xtt(2,:)
!!$  xtp(mth1x,:) = xtp(2,:)
!!$  xpp(mth1x,:) = xpp(2,:)
  call i2mex_getRmagnetic(xma, ier)
  call i2mex_getZmagnetic(zma, ier)
  call i2mex_getFrame(xin, xout, zbot, ztop, ier)
  call i2mex_error(ier)
  r0 = (xin  + xout)/2._r8
  z0 = (ztop + zbot)/2._r8
!!$  xma = sum(x(1:mth,1))/real(mth, r8)
  xmin= minval(x)
  xmax= maxval(x)
  call cdfPutvar(iu, 'Xmin', xmin, iok)
  call cdfPutvar(iu, 'Xmax', xmax, iok)
!!$  r0 = (xmax+xmin)/2._r8
  
  call cdfPutvar(iu, 'X0', r0, iok)
  call cdfPutvar(iu, 'Xmag', xma, iok)
  call cdfPutvar(iu, 'xa', x, iok)

  call i2mex_getZ(mth1, nosurf, the, psi, z, ier)
  call i2mex_getGradZ(mth1, nosurf, the, psi, zt, zp, ier)
  call i2mex_getZtt(mth1, nosurf, the, psi, ztt, ier)
  call i2mex_getZtp(mth1, nosurf, the, psi, ztp, ier)
  call i2mex_getZpp(mth1, nosurf, the, psi, zpp, ier)
  call i2mex_error(ier)
!!$  z(mth1x,:) = z(2,:)
!!$  zt(mth1x,:) = zt(2,:)
!!$  zp(mth1x,:) = zp(2,:)
!!$  ztt(mth1x,:) = ztt(2,:)
!!$  ztp(mth1x,:) = ztp(2,:)
!!$  zpp(mth1x,:) = zpp(2,:)
!!$  zma = sum(z(1:mth,1))/real(mth, r8)
  zmin= minval(z)
  zmax= maxval(z)
  call cdfPutvar(iu, 'Zmin', zmin, iok)
  call cdfPutvar(iu, 'Zmax', zmax, iok)
  call cdfPutvar(iu, 'Z0', z0, iok)
  call cdfPutvar(iu, 'Zmag', zma, iok)
  call cdfPutvar(iu, 'za', z, iok)

  call i2mex_getB0(res, ier)
  call cdfPutvar(iu, 'B0', res, iok)
  call i2mex_getPlasmaCurrent(mth1, nosurf, the, psi, res, ier)
  call cdfPutvar(iu, 'Ip', res, iok)
  call i2mex_getBeta(mth1, nosurf, the, psi, res, ier)  
  call cdfPutvar(iu, 'Beta', res, iok)
  call i2mex_getBetaStar(mth1, nosurf, the, psi, res, ier)
  call cdfPutvar(iu, 'BetaStar', res, iok)
  call i2mex_getBetaN(mth1, nosurf, the, psi, res, ier)
  call cdfPutvar(iu, 'BetaN', res, iok)
  call i2mex_getLi(mth1, nosurf, the, psi, res, ier)
  call cdfPutvar(iu, 'li', res, iok)
  call i2mex_getP(nosurf, psi, p, iok)
  call i2mex_getVolumeAveragedPressure(mth1, nosurf, the, psi, res, ier)
  if(res /= 0.0_r8) res = p(1)/res
  call cdfPutvar(iu, 'PPF', res, iok)
  call i2mex_getUpsilon(res, iok)
  call cdfPutvar(iu, 'Upsiln', res, iok)
  call i2mex_getJdotBOverBSquare(mth1, nosurf, the, psi, temp1, ier)
  call cdfPutvar(iu, 'JdotBOverBSquare', temp1, iok)
  call i2mex_getV(nosurf, psi, v_of_psi, ier)
  call i2mex_getVp(nosurf, psi, vp_of_psi, ier)
  call i2mex_getQp(nosurf, psi, qp, ier)
  call i2mex_getQpp(nosurf, psi, qpp, ier)
  call i2mex_getQ(nosurf, psi, q, iok)
  temp1 = 2._r8*v_of_psi*temp1/(q*vp_of_psi)
  call cdfPutvar(iu, 'local_shear_s', temp1, iok)

  call i2mex_getGlasserGreenJohnsonEFH(nosurf, psi, e, f, h, ier)
  call i2mex_error(ier)
  temp1 = e + f + h - 0.25_r8
  call cdfPutvar(iu, 'Di', temp1, iok)
  temp1 = e + f + h**2
  call cdfPutvar(iu, 'Dr', temp1, iok)
!!$  call i2mex_getV(nosurf, psi, v_of_psi, ier)
  call i2mex_error(ier)
  temp1 = sqrt(v_of_psi/(2.0_r8*i2mex_pi_r8**2 *r0)) 
  call cdfPutvar(iu, 'surface_averaged_radius', temp1, iok)
  call i2mex_getPp(nosurf, psi, pp, iok)
  temp1 = - temp1 * vp_of_psi * pp/(2.0_r8*i2mex_pi_r8**2)
  call cdfPutvar(iu, 'ballooning_alpha', temp1, iok)
  call cdfPutvar(iu, 'PsiBig', i2mex_twopi_r8*psi, iok)
  call cdfPutvar(iu, 'V_of_psi', v_of_psi, iok)
  call cdfPutvar(iu, 'Vp_of_psi', vp_of_psi, iok)
  call cdfPutvar(iu, 'pa', p, iok)
  call cdfPutvar(iu, 'ppa', pp, iok)
  call cdfPutvar(iu, 'qa', q, iok)
  call cdfPutvar(iu, 'qpa', qp, iok)
  call cdfPutvar(iu, 'qppa', qpp, iok)
  call i2mex_getG(nosurf, psi, g, iok)
  call cdfPutvar(iu, 'ga', g, iok)
  call i2mex_getGp(nosurf, psi, gp, iok)
  call cdfPutvar(iu, 'gpa', gp, iok)
  
  call i2mex_getSpest1(nosurf, psi, temp1, iok)
  call cdfPutvar(iu, 'spest1', temp1, iok)

!!$  call cdfPutvar(iu, 'dpsi_drho_fa', temp1, iok)
!!$  call cdfPutvar(iu, 'dpsi2_drho2_fa', temp1, iok)
  call cdfPutvar(iu, 'xinf', x(:,nosurf), iok)
  call cdfPutvar(iu, 'zinf', z(:,nosurf), iok)
  call cdfPutvar(iu, 'xpth', xt, iok)
  call cdfPutvar(iu, 'zpth', zt, iok)
  call cdfPutvar(iu, 'xpsi', xp, iok)
  call cdfPutvar(iu, 'zpsi', zp, iok)
  allocate(xsq(mth1x, nosurf))
  allocate(xsqdps(mth1x, nosurf))
  xsq = x**2
  xsqdps = 2.0_i2mex_r8* x * xp
!!$  call cdfPutvar(iu, 'xsq', xsq, iok)
!!$  call cdfPutvar(iu, 'xsqdps', xsqdps, iok)
!!$  call cdfPutvar(iu, 'xsqdth', temp2, iok)

  allocate(gtt(mth1x, nosurf), gtp(mth1x, nosurf), gpp(mth1x, nosurf))
  
  ! metric
  ! gtt -> |grad the|^2
  ! gtp -> grad the . grad psi
  ! gpp -> |grad psi|^2
  call i2mex_getMetric(mth1, nosurf, the, psi, gtt, gtp, gpp, ier)
  call cdfPutvar(iu, 'grpssq', gpp, iok)
  call cdfPutvar(iu, 'grthsq', gtt, iok)
  call cdfPutvar(iu, 'grpsth', gtp, iok)

  ! d/d the of metric
  ! gtt -> |grad the|^2
  ! gtp -> grad the . grad psi
  ! gpp -> |grad psi|^2
  call i2mex_getMetricT(mth1, nosurf, the, psi, gtt, gtp, gpp, ier)
  call cdfPutvar(iu, 'gptdth', gtp, iok)
  call cdfPutvar(iu, 'gpsdth', gpp, iok)

  allocate(jac(mth1x, nosurf))
  allocate(jacp(mth1x, nosurf))

  call i2mex_getJ(mth1, nosurf, the, psi, jac, ier)
  call cdfPutvar(iu, 'xjacob', jac, iok)
  call i2mex_getJp(mth1, nosurf, the, psi, jacp, ier)
  call cdfPutvar(iu, 'xjprym', jacp, iok)


  !
  ! compute delta = theta_p - theta
  !
  call i2mex_getDelta(mth1, nosurf, the, psi, temp2, ier)
  call i2mex_error(ier)

!!$  do i=1, nosurf
!!$     tempt(1:mth1) = g(i)*jac(1:mth1,i)/ &
!!$          & (q(i)* xsq(1:mth1,i))
!!$     call i2mex_integratePeriodic1d(mth1, the, tempt, mth1, the, temp2(1,i))
!!$     temp2(1:mth1,i) = temp2(1:mth1,i) - the(1:mth1)
!!$  enddo
!!$
!!$    ! extrapolate to axis
!!$
!!$    do i=1, mth1
!!$       temp2(i,1) = FCCCC0(temp2(i,2), temp2(i,3), temp2(i,4), temp2(i,5), &
!!$            &                 psi(2),  psi(3),  psi(4),  psi(5), &
!!$            & psi(1))
!!$    enddo

  call cdfPutvar(iu, 'delta', temp2, iok)

  !
  ! compute d(q*delta)/d psi
  !
  call i2mex_getQDeltaP(mth1, nosurf, the, psi, temp2, ier)
  call i2mex_error(ier)

!!$  do i=1, nosurf
!!$     tempt(1:mth1) = gp(i)*jac(1:mth1,i)/xsq(1:mth1,i) &
!!$          & + g(i)*jacp(1:mth1,i)/xsq(1:mth1,i) &
!!$          & - g(i)*jac(1:mth1,i)*xsqdps(1:mth1,i)/xsq(1:mth1,i)**2
!!$     call i2mex_integratePeriodic1d(mth1, the, tempt, mth1, the, temp2(1,i))
!!$     temp2(1:mth1,i) = temp2(1:mth1,i) - qp(i)*the(1:mth1)
!!$  enddo
!!$
!!$    ! extrapolate to axis
!!$
!!$    do i=1, mth1
!!$       temp2(i,1) = FCCCC0(temp2(i,2), temp2(i,3), temp2(i,4), temp2(i,5), &
!!$            &                 psi(2),  psi(3),  psi(4),  psi(5), &
!!$            & psi(1))
!!$    enddo
!!$

  call cdfPutvar(iu, 'qdelp', temp2, iok)

  call i2mex_getQDeltaT(mth1, nosurf, the, psi, temp2, ier)
  call i2mex_error(ier)

  call cdfPutvar(iu, 'qdelt', temp2, iok)
  

  call i2mex_getGsError(mth1, nosurf, the, psi, temp2, ier)
  call i2mex_error(ier)
  
  call cdfPutvar(iu, 'gserror', temp2, iok)

  DEALLOCATE(  temp1  )
  DEALLOCATE(  tempt  )
  deallocate(v_of_psi, vp_of_psi, stat=iok)
  if(iok/=0) print*,'failed to deallocate v_of_psi, vp_of_psi '
  deallocate(g, gp, q, qp, qpp, p, pp, stat=iok)
  if(iok/=0) print*,'failed to deallocate g, gp, q, qp, qpp, p, pp '

  deallocate(gtt, gtp, gpp, stat=iok)
  if(iok/=0) print*,'failed to deallocate gtt, gtp, gpp '
  deallocate(xsq, stat=iok)
  if(iok/=0) print*,'failed to deallocate xsq'
  deallocate(xsqdps, stat=iok)
  if(iok/=0) print*,'failed to deallocate xsqdps '
  deallocate(jac, stat=iok)
  if(iok/=0) print*,'failed to deallocate jac '
  deallocate(jacp, stat=iok)
  if(iok/=0) print*,'failed to deallocate jacp'
  deallocate(temp2, stat=iok)
  if(iok/=0) print*,'failed to deallocate temp2 '
  deallocate(x, xt, xp, xtt, xtp, xpp, stat=iok)
  if(iok/=0) print*,'failed to deallocate x, xt, xp, xtt, xtp, xpp'
  deallocate(z, zt, zp, ztt, ztp, zpp, stat=iok)
  if(iok/=0) print*,'failed to deallocate z, zt, zp, ztt, ztp, zpp'

  call cdfCls(iu)

end subroutine i2mex_ToMapdsk
  

  
