subroutine i2mex_ReadEquilibrium(inputFormat, inputFile, it_orientation, ier)
  !
  ! Read Grad Shafranov equibrium data stored in format inputFormat from
  ! file inputFile, and initialize i2mex.
  ! Supported file formats:
  ! inputFormat = -1  inputFile is in CHEASE 'INP1' binary file format
  ! inputFormat = +1  inputFile is in CHEASE netCDF 'inp1.cdf' file format
  ! inputFormat = +2  inputFile is in JSOLVER netCDF 'eqdsk.cdf' file format
  ! inputFormat = +3  inputFile is in EFIT G_EQDSK file format.
  ! inputFormat = +4  inputFile is in EFIT G_EQDSK file format,
  !                   rerunning ESC taking q profile as input
  ! inputFormat = +5  Menard or psipqgRZ format in MKS units
  ! inputFormat = +6  FreeqBe or free BC equilibrium by Belova
  ! inputFormat = +7  read ceo.cdf or MAPDSK.NC 
  ! inputFormat = +8  L. Don Pearlstein's inverse equilibrium code
  !
  ! it_orientation determines the theta orientation:
  ! it_orientation = +1 for clockwise
  ! it_orientation = -1 for counterclockwise.
  !
  use i2mex_mod
  implicit none
  integer, intent(in) :: inputFormat
  character*(*), intent(in) :: inputFile
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier

  ier = 0
  print*,'input file is ', inputFile
  select case(inputFormat)
 
     case(:-2)
        ier = 104
        return
     case(-1)
        call i2mex_fromInp1Binary(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(0)
        ier = 104
        return
     case(1)
        call i2mex_fromInp1CDF(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(2)
        call i2mex_fromEqdskCDF(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(3)
        call i2mex_fromGeqdsk(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(4)
        call i2mex_fromGeqdsk(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(5)
        call i2mex_fromMenard(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(6)
        call i2mex_fromFreeqBe(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(7)
        call i2mex_fromCeo(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(8)
        call i2mex_fromLDPBinary(inputFile, it_orientation, ier)
        call i2mex_error(ier)
     case(9:)
        ier = 104
        return
 
  end select
 
end subroutine i2mex_ReadEquilibrium
 
subroutine i2mex_fromFreeqbe(filename, it_orientation, ier)
 
  ! read free equilibrium format by Elena Belova
 
   use i2mex_mod
   use freeqbe_mod
   use ezspline_obj
   use ezspline
   use cont_mod
 
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
 
  type(freeqbe) :: feb
  type(ezspline2) :: pspl
 
  real(r8), dimension(:), allocatable :: p, q, g, psi, psi_feb, psis, &
       & the, xgrid, zgrid, tlike, q_est
  real(r8), dimension(:, :), allocatable :: x, z
  real(r8), dimension(:, :), allocatable :: xin, zin
  real(r8) :: pmax, pmin, xmin, xmax, zmin, zmax, surf, xm, zm, xp, zp, &
       & pmag, psep, ttot
  integer iok, i, j, ioutboard(1)
 
  real(r8), dimension(:,:), allocatable :: integrand
  real(r8), dimension(:), allocatable :: qres
  real(r8) x0, z0, x1, z1, x2, z2, det, xsep, zsep
 
  integer noutb
  real(r8), dimension(:), allocatable :: xoutb, zoutb, psioutb
 
  ier = 0
 
  call feb_init(feb, filename, iok)
  call feb_error(iok)
  if(iok/=0) then
     ier = 117
     return
  endif
 
  call feb_toMKS(feb, iok)
 
  xmin = minval(feb% r)
  xmax = maxval(feb% r)
  zmin = minval(feb% z)
  zmax = maxval(feb% z)
 
 
  allocate(xgrid(feb% nr))
  allocate(zgrid(feb% nz))
 
  xgrid=xmin + (xmax - xmin)*(/ (real(i-1,r8)/real(feb% nr-1,r8), &
       & i=1,feb% nr) /)
  zgrid=zmin + (zmax - zmin)*(/ (real(i-1,r8)/real(feb% nz-1,r8), &
       & i=1,feb% nz) /)
 
!!$  i2mex_nt1 = 129
!!$  i2mex_ns =  101
 
  allocate(psi(i2mex_ns), psi_feb(feb% npsi), p(i2mex_ns), g(i2mex_ns), q(i2mex_ns), the(i2mex_nt1))
  allocate(tlike(i2mex_ns), q_est(i2mex_ns))
  allocate(x(i2mex_nt1, i2mex_ns), z(i2mex_nt1, i2mex_ns))
  allocate(xin(i2mex_nt1, i2mex_ns), zin(i2mex_nt1, i2mex_ns))
  allocate(psis(i2mex_nt1))
 
  xm = feb% rzaxis(1)
  zm = feb% rzaxis(2)
 
  ! psi - spline setup
 
  call ezspline_init(pspl, feb% nr, feb% nz, (/0,0/), (/0,0/), iok)
  call ezspline_error(iok)
  pspl%x1 = xgrid
  pspl%x2 = zgrid
  !pspl%isHermite = 1 ! Akima Hermite
  call ezspline_setup(pspl, feb% psi2, iok)
  call ezspline_error(iok)
 
  ! Find intersection of separatrix with outer mid plane
  noutb = 101
  allocate(xoutb(noutb), zoutb(noutb), psioutb(noutb))
  xoutb = xm + (xmax-xm)*(/ (real(i-1,r8)/real(noutb, r8), i=1, noutb) /)
  zoutb = (/ (zm, i=1, noutb) /)
  call ezspline_interp(pspl, noutb, xoutb, zoutb, psioutb, iok)
  call ezspline_error(iok)
  ! rough estimate
  ioutboard = minloc((psioutb - feb% psi0b(2))**2)
  ! refine
  call i2mex_interp1d(ioutboard(1)+4, &
       & psioutb(1:ioutboard(1)+4), &
       &   xoutb(1:ioutboard(1)+4), 1, feb% psi0b(2), xsep)
  zsep = zm
  deallocate(xoutb, zoutb, psioutb)
 
  ! toroidal current sign
  cont_sign = 1.0_r8
  if(feb %psi0b(2) < feb% psi0b(1)) cont_sign = -1.0_r8
 
  ! We like to have psi = 0 on magnetic axis, and increasing
  ! outwards (this will change the current orientation but that's
  ! ok).
  feb% psi2     = cont_sign*(feb% psi2     - feb% psi0b(1))
  feb% psi0b(2) = cont_sign*(feb% psi0b(2) - feb% psi0b(1))
  feb% psi0b(1) = 0.0_r8
  ! set back
  cont_sign = 1.0_r8
 
  ! interpolate to get psi on axis
  call ezspline_interp(pspl, xm, zm, pmag, iok)
  call ezspline_error(iok)
  ! psi on separatrix
  call ezspline_interp(pspl, xsep, zsep, psep, iok)
  call ezspline_error(iok)
 
  do j = 2, i2mex_ns
 
     surf = real(j-1, r8)/real(i2mex_ns-1, r8)
     xp = xm + surf * i2mex_LAST_NORM_SURFACE_IS*(xsep - xm)
     zp = zm + surf * i2mex_LAST_NORM_SURFACE_IS*(zsep - zm)
 
     ! rough estimate for q (used to determine angle integration step)
     cont_q = 1.0_r8 + 20._r8* surf**2
 
     call i2mex_contour(feb% nr, feb% nz, xgrid, zgrid, feb% psi2,  &
       & xm, zm, xp, zsep, i2mex_nt1, ttot, xin(1, j), zin(1, j), iok)
     call i2mex_error(iok)
 
     ! get psi on contour
     call ezspline_interp(pspl, i2mex_nt1, xin(:, j), zin(:, j), psis, iok)
     call ezspline_error(iok)
 
     ! average psi on contour
     psi(j) = sum(psis)/real(i2mex_nt1, r8)
 
     ! psi is in Wb/Rad
     !psi(j) = psi(j) / i2mex_twopi_r8
 
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
 
  psi_feb = 0._r8 + ((pmax-pmin)/i2mex_LAST_NORM_SURFACE_IS)* &
       & (/ (real(i-1, r8)/real(feb% npsi-1,r8), i=1, feb% npsi) /)
 
 
  call ezspline_free(pspl, iok)
 
  the = i2mex_twopi_r8* (/ (real(i-1,r8)/real(i2mex_nt1-1,r8), i=1, i2mex_nt1) /)
 
  ! pressure
  call i2mex_interp1d(feb% npsi, psi_feb, &
       & feb% p * 4.e-7_r8 * i2mex_pi_r8, &
       & i2mex_ns, psi, p)
  ! g
  call i2mex_interp1d(feb% npsi, psi_feb, &
       & abs(feb% RB)                        , &
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
 
  call i2mex_init('FREEQBE_'//filename, &
       i2mex_nt1, i2mex_ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  deallocate(xgrid)
  deallocate(zgrid)
  deallocate(psi, psi_feb, p, g, q, the)
  deallocate(tlike, q_est)
  deallocate(x, z)
  deallocate(xin, zin)
  deallocate(psis)
 
  call feb_free(feb, iok)
  if(iok/=0) ier = 149
 
 
end subroutine i2mex_fromFreeqbe
 
subroutine i2mex_fromInp1CDF(filename, it_orientation, ier)
  !
  ! read chease netCDF output file
  !
  use i2mex_mod
  use ezcdf
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
  integer nxx(3)
 
  integer nt1, ns, iok, inp1, isym, i, j, iu
  real(r8) :: pedge, gedge
  real(r8), dimension(:), allocatable :: p, pm, pp, g, gp, q, qp, &
       & the, psi, psim
  real(r8), dimension(:,:), allocatable :: xin, zin, x, z
  real(r8) :: x0, x1, x2, z0, z1, z2, det
 
  ier = 0
 
  call ezcdf_open(iu, filename, 'r', iok)
  IF(iok /= 0) THEN
     print*,'**ERROR** could not open ', filename
     ier = 90
     return
  ENDIF
 
  ! read dimensions
 
  call cdfGetVar(iu, 'nxx', nxx, iok)
  if(iok/=0) then
     ier = 91
     return
  endif
 
  ns   =  nxx(2)
  isym =  nxx(3)
 
  allocate(pm(ns-1), stat=iok)
  allocate(p(ns), stat=iok)
  allocate(pp(ns), stat=iok)
  allocate(q(ns), stat=iok)
  allocate(qp(ns), stat=iok)
  allocate(g(ns), stat=iok)
  allocate(gp(ns), stat=iok)
  allocate(psi(ns), stat=iok)
  allocate(psim(ns-1), stat=iok)
 
  ! read profiles
 
  call cdfGetVar(iu, 'pm', pm, iok) ! zone centre
  call cdfGetVar(iu, 'pp', pp, iok)
  call cdfGetVar(iu, 'q', q, iok)
  call cdfGetVar(iu, 'g', g, iok)
  call cdfGetVar(iu, 'gp', gp, iok)
  call cdfGetVar(iu, 'psival', psi, iok)
  call cdfGetVar(iu, 'psivalm', psim, iok) ! zone centre
  if(iok/=0) then
     ier = 94
     return
  endif
 
  if(isym == 1) then
     nt1 = 2*(nxx(1)-1) + 1
  else
     nt1 = nxx(1) + 1
  endif
 
  allocate(the(nt1))
  allocate(x(nt1, ns))
  allocate(z(nt1, ns))
 
 
  if(isym == 1) then
 
     ! plasma is up
 
     allocate(xin(nxx(1)+3, ns))
     allocate(zin(nxx(1)+3, ns))
     call cdfGetVar(iu, 'x', xin, iok)
     call cdfGetVar(iu, 'z', zin, iok)
     if(iok/=0) then
        ier = 94
        return
     endif
     call cdfCls(iu)
 
     ! get theta orientation
     x0 = (maxval(xin(1:,1))+minval(xin(:,1)))/2._i2mex_r8
     z0 = (maxval(zin(1:,1))+minval(zin(:,1)))/2._i2mex_r8
     x1 = xin(1, ns)
     z1 = zin(1, ns)
     x2 = xin(2, ns)
     z2 = zin(2, ns)
     det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
     if( it_orientation*det < 0.0_i2mex_r8 ) then
        ! keep orientation
        do j = 1, ns
           do i = 1, nxx(1)+1
              x(i,j) = xin(i+2, j)
              z(i,j) = zin(i+2, j)
              x(nt1-i+1,j) = xin(i+2, j)
              z(nt1-i+1,j) =-zin(i+2, j)
           enddo
        enddo
     else
        ! change theta orientation
        do j = 1, ns
           do i = 1, nxx(1)+1
              x(nt1-i+1,j) = xin(i+2, j)
              z(nt1-i+1,j) = zin(i+2, j)
              x(i,j) = xin(i+2, j)
              z(i,j) =-zin(i+2, j)
           enddo
        enddo
     endif
 
     deallocate(xin)
     deallocate(zin)
 
     z(:,1) = z0
 
  else
 
     ! plasma is up-down
 
     allocate(xin(nt1+2, ns), stat=iok)
     allocate(zin(nt1+2, ns), stat=iok)
 
     call cdfGetVar(iu, 'x', xin, iok)
     call cdfGetVar(iu, 'z', zin, iok)
     if(iok/=0) then
        ier = 94
        return
     endif
     call cdfCls(iu)
 
     ! get theta orientation
     x0 = (maxval(xin(:,1))+minval(xin(:,1)))/2._i2mex_r8
     z0 = (maxval(zin(:,1))+minval(zin(:,1)))/2._i2mex_r8
     x1 = xin(1, ns)
     z1 = zin(1, ns)
     x2 = xin(2, ns)
     z2 = zin(2, ns)
     det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
 
     if( it_orientation*det < 0.0_i2mex_r8 ) then
        do i = 1, nt1
           x(i, 1:ns) = xin(i+2, 1:ns)
           z(i, 1:ns) = zin(i+2, 1:ns)
        enddo
     else
        do i = 1, nt1
           x(i, 1:ns) = REAL(xin(nt1+2-i+1, 1:ns), i2mex_r8)
           z(i, 1:ns) = REAL(zin(nt1+2-i+1, 1:ns), i2mex_r8)
        enddo
     endif
 
     deallocate(xin)
     deallocate(zin)
 
  endif
 
  z(nt1,:) = z(1,:)
 
  ! assume uniform mesh in theta
 
  the=i2mex_twopi_r8*(/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)
 
  ! integrate pp & gp to get p, resp. g
 
  pedge = pm(ns-1)  + &
       & (psi(ns)-psim(ns-1))*0.5_i2mex_r8*(pp(ns)+pp(ns-1))! from input file
  call i2mex_integrate1d(ns, psi, pp, ns, psi, p)
  p = p - p(ns) + pedge
  gedge = g(ns)
  call i2mex_integrate1d(ns, psi, gp, ns, psi, g)
  g = g - g(ns) + gedge
 
  ! shift
  psi = abs(psi - psi(1))
 
  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & nt1, ns, the, psi, g, x, z, ier)
  endif
 
  call i2mex_init('inp1_netCDF_input_'//filename, &
       nt1, ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  if ( iok/=0 ) ier = 3
 
  deallocate(p, stat=iok)
  deallocate(pm, stat=iok)
  deallocate(pp, stat=iok)
  deallocate(q, stat=iok)
  deallocate(qp, stat=iok)
  deallocate(g, stat=iok)
  deallocate(gp, stat=iok)
  deallocate(the, stat=iok)
  deallocate(psi, stat=iok)
  deallocate(psim, stat=iok)
  deallocate(x, stat=iok)
  deallocate(z, stat=iok)
 
  if( iok/=0 ) ier = 96
 
 
end subroutine i2mex_fromInp1CDF
 
 
subroutine i2mex_fromEqdskCDF(filename, it_orientation, ier)
  !
  ! read chease JSOLVER output file in netCDF format
  !
  use ezcdf
  use i2mex_mod
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)

  include 'cubinterp.h'

  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier

  integer, dimension(:), allocatable :: nxx
  real(r8), dimension(:), allocatable :: axx

  integer nt1, ns, iok, inp1, isym, i, j, iu, nths, dimlens(3)
  real(r8) :: pedge, gedge, r0
  real(r8), dimension(:), allocatable :: p, pp, g, gp, q, qp, &
       & the, psi
  real(r8), dimension(:,:), allocatable :: xin, zin, x, z
  character*3 :: xtype
  real(r8) :: x0, x1, x2, z0, z1, z2, det
  integer :: nxy(10)
  integer :: iouter(1), index

  ier = 0

  call ezcdf_open(iu, filename, 'r', iok)
  IF(iok /= 0) THEN
     print*,'**ERROR** could not open ', filename
     ier = 34
     return
  ENDIF

  ! read dimensions

  call cdfInqVar(iu, 'nxx', dimlens, xtype, iok)
  ALLOCATE(nxx(dimlens(1)))
  call cdfGetVar(iu, 'nxx', nxx, iok)
  if(iok/=0) then
     ier = 35
     return
  endif

  call cdfGetVar(iu, 'nxy', nxy, iok)
  if(iok/=0) then
     ier = 35
     return
  endif


  call cdfInqVar(iu, 'axx', dimlens, xtype, iok)
  ALLOCATE(axx(dimlens(1)))
  call cdfGetVar(iu, 'axx', axx, iok)
  if(iok/=0) then
     ier = 35
     return
  endif

  ! 5th element of axx appears to contain r0
  r0 = axx(5)

  isym =  nxy(8)

  nt1  =  nxx(1) + 1
  if(isym ==1) nt1 = (nxx(1)-1)*2 + 1

  ns   =  nxx(2)

  call cdfInqVar(iu, 'x', dimlens, xtype, iok)
  nths=dimlens(1)

!!$  DEALLOCATE(nxx)
  DEALLOCATE(axx)


  if(isym == 1) print*,'equilibrium up-down symmetric'

  allocate(p(ns), stat=iok)
  allocate(pp(ns), stat=iok)
  allocate(q(ns), stat=iok)
  allocate(qp(ns), stat=iok)
  allocate(g(ns), stat=iok)
  allocate(gp(ns), stat=iok)
  allocate(the(nt1), stat=iok)
  allocate(psi(ns), stat=iok)
  allocate(x(nt1, ns), stat=iok)
  allocate(z(nt1, ns), stat=iok)

  if( iok/=0 ) then
     print*,'***ERROR** memory allocation '
     ier = 5
     return
  endif

  ! read profiles

  call cdfGetVar(iu, 'p', p, iok)
  call cdfGetVar(iu, 'ppxx', pp, iok)
  call cdfGetVar(iu, 'q', q, iok)
  call cdfGetVar(iu, 'gxx', g, iok)
  ! need to multiply g fcts by r0 (for some reason)
  g = r0* g
  call cdfGetVar(iu, 'gpx', gp, iok)
  gp = r0* gp
  call cdfGetVar(iu, 'psival', psi, iok)
  if(iok/=0) then
     ier = 35
     return
  endif

  ! read grid

  if(isym == 1) then

     ! plasma is up

     allocate(xin(nxx(1)+1, ns), stat=iok)
     if(iok/=0) print *,' failed to allocate xin'
     allocate(zin(nxx(1)+1, ns), stat=iok)
     if(iok/=0) print *,' failed to allocate zin'
     call cdfGetVar(iu, 'x', xin, iok)
     call cdfGetVar(iu, 'z', zin, iok)
     if(iok/=0) then
        ier = 95
        return
     endif

!!$     ! locate outerboard index
!!$     iouter = maxloc(xin(:,ns))

     ! get theta orientation
     x0 = (maxval(xin(1:,1))+minval(xin(:,1)))/2._i2mex_r8
     z0 = (maxval(zin(1:,1))+minval(zin(:,1)))/2._i2mex_r8
     x1 = xin(1, ns)
     z1 = zin(1, ns)
     x2 = xin(2, ns)
     z2 = zin(2, ns)
     det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)

     if( it_orientation*det < 0.0_i2mex_r8 ) then
        ! keep orientation
        do j = 1, ns
           do i = 1, nxx(1)
              x(i,j) = xin(i+1, j)
              z(i,j) = zin(i+1, j)
              x(nt1-i+1,j) = xin(i+1, j)
              z(nt1-i+1,j) =-zin(i+1, j)
           enddo
        enddo
     else
        ! change theta orientation
        do j = 1, ns
           do i = 1, nxx(1)
              x(nt1-i+1,j) = xin(i+1, j)
              z(nt1-i+1,j) =+zin(i+1, j)
              x(i,j) = xin(i+1, j)
              z(i,j) =-zin(i+1, j)
           enddo
        enddo
     endif

     deallocate(xin)
     deallocate(zin)

     z(:,1) = z0
     x(:,1) = x0

  else

     ! plasma is up-down

     allocate(xin(nths, ns), stat=iok)
     allocate(zin(nths, ns), stat=iok)


     call cdfGetVar(iu, 'x', xin, iok)
     call cdfGetVar(iu, 'z', zin, iok)
     if(iok/=0) then
        ier = 35
        return
     endif

     ! locate outerboard index
     iouter = maxloc(xin(:,ns))

     ! reverse theta orientation if necessary

     x0 = sum(xin(1:nt1-1,1))/real(nt1-1)
     z0 = sum(zin(1:nt1-1,1))/real(nt1-1)
     x1 = xin(1, ns)
     z1 = zin(1, ns)
     x2 = xin(2, ns)
     z2 = zin(2, ns)
     det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)

     if( it_orientation*det < 0.0_i2mex_r8 ) then
        print *,'keep theta orientation'
        do i = 1, nt1
           index =  i + iouter(1) - 1
           if(index <   1) index = index + nt1 - 1
           if(index > nt1) index = index - nt1 + 1
           x(i, 1:ns) = xin(index, 1:ns)
           z(i, 1:ns) = zin(index, 1:ns)
        enddo
     else
        print *,'reverse theta orientation'
        do i = 1, nt1
           index = nt1 - i + iouter(1)
           if(index <   1) index = index + nt1 - 1
           if(index > nt1) index = index - nt1 + 1           
           x(i, 1:ns) = xin(index, 1:ns)
           z(i, 1:ns) = zin(index, 1:ns)
        enddo
     endif

     deallocate(xin)
     deallocate(zin)



  endif

  call cdfCls(iu)

  ! assume uniform mesh in theta

  the=i2mex_twopi_r8*(/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)

  ! integrate pp & gp to get p, resp. g

  pedge = p(ns) ! from input file
  call i2mex_integrate1d(ns, psi, pp, ns, psi, p)
  p = p - p(ns) + pedge
  gedge = g(ns)
  call i2mex_integrate1d(ns, psi, gp, ns, psi, g)
  g = g - g(ns) + gedge

  psi = abs(psi - psi(1))

  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & nt1, ns, the, psi, g, x, z, ier)
  endif

  ! Extrapolate to axis. Jsolver may need it...

!!$  do i = i2mex_nqextra, 1, -1
!!$     call i2mex_toAxis(1, ns-i+1, the, psi(i), 0.0_i2mex_r8, q(i), iok)
!!$  enddo
  ! use a few points spread near the axis to extrapolate 

  if(4*i2mex_nqextra < ns .and. i2mex_nqextra > 0) then
     ! to the axis
     do i = 1, i2mex_nqextra
        q(i) = FCCCC0( &
             &   q(  i2mex_nqextra), &
             &   q(2*i2mex_nqextra), &
             &   q(3*i2mex_nqextra), &
             &   q(4*i2mex_nqextra), &
             & psi(  i2mex_nqextra), &
             & psi(2*i2mex_nqextra), &
             & psi(3*i2mex_nqextra), &
             & psi(4*i2mex_nqextra), &
             & psi(i) )
     enddo
  endif

  ! to the edge, last surface is not accurate
     i = ns
     q(i) = FCCCC0( &
          &   q(ns-2), &
          &   q(ns-3), &
          &   q(ns-4), &
          &   q(ns-5), &
          & psi(ns-2), &
          & psi(ns-3), &
          & psi(ns-4), &
          & psi(ns-5), &
          & psi(i) )


 
  call i2mex_init('EQDSK_netCDF_input_'//filename, &
       & nt1, ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  if ( iok/=0 ) ier = 3
 
  deallocate(nxx, stat=iok)
  deallocate(p, stat=iok)
  deallocate(pp, stat=iok)
  deallocate(q, stat=iok)
  deallocate(qp, stat=iok)
  deallocate(g, stat=iok)
  deallocate(gp, stat=iok)
  deallocate(the, stat=iok)
  deallocate(psi, stat=iok)
!!$  deallocate(xin, stat=iok)
!!$  deallocate(zin, stat=iok)
  deallocate(x, stat=iok)
  deallocate(z, stat=iok)
 
  if( iok/=0 ) ier = 4
 
 
end subroutine i2mex_fromEqdskCDF
 
subroutine i2mex_fromInp1Binary(filename, it_orientation, ier)
  !
  ! read chease binary output file
  !
  use i2mex_mod
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
  REAL(r8) AXX(5)
  REAL(r8) AXY
  integer nxx(3)
 
  integer nt1, ns, iok, inp1, isym, i, j
  real(r8) :: pedge, gedge
  real(r8), dimension(:), allocatable :: zdum, p, pm, pp, g, gp, q, qp, &
       & the, psi, psim
  real(r8), dimension(:,:), allocatable :: xin, zin, x, z
  real(r8) :: x0, x1, x2, z0, z1, z2, det
 
  ier = 0
 
  inp1 = 1
  OPEN(inp1, FILE=filename, STATUS='old',                                 &
       &            FORM='unformatted', IOSTAT=iok)
  IF(iok /= 0) THEN
     print*,'**ERROR** could not find file ', filename
     ier = 1
     return
  ENDIF
  REWIND(inp1)
 
  ! read dimensions
 
  read(inp1) (nxx(i),i=1,3)
 
  nt1  =  nxx(1) + 1
  ns   =  nxx(2)
  isym =  nxx(3)
 
  if(isym == 1) then
     print*,'***ERROR** equilibrium appears to be up-down symmetric'
     print*,'this version only supports up-down asymmetric'
     print*,'plasmas'
     ier = 2
     return
  endif
 
  allocate(zdum(ns), stat=iok)
  allocate(p(ns), stat=iok)
  allocate(pm(ns-1), stat=iok)
  allocate(pp(ns), stat=iok)
  allocate(q(ns), stat=iok)
  allocate(qp(ns), stat=iok)
  allocate(g(ns), stat=iok)
  allocate(gp(ns), stat=iok)
  allocate(the(nt1), stat=iok)
  allocate(psi(ns), stat=iok)
  allocate(psim(ns-1), stat=iok)
  allocate(xin(nt1+2, ns), stat=iok)
  allocate(zin(nt1+2, ns), stat=iok)
  allocate(x(nt1, ns), stat=iok)
  allocate(z(nt1, ns), stat=iok)
 
  if( iok/=0 ) then
     print*,'***ERROR** memory allocation '
     ier = 5
     return
  endif
 
 
  ! read profiles
 
  read(inp1) (zdum(i),i=1,5)
  read(inp1) (pm(i),i=1,ns-1) ! pressure at mid-intervals
  p(ns) = 0.0_r8
  pedge = p(ns)
  read(inp1) (pp(i),i=1,ns)
  read(inp1) (q(i),i=1,ns)
  read(inp1) (zdum(i),i=1,ns)
  read(inp1) (g(i),i=1,ns)
  gedge = g(ns)
  read(inp1) (gp(i),i=1,ns)
  read(inp1) (zdum(i),i=1,ns)
  read(inp1) (zdum(i),i=1,ns)
  read(inp1) (psi(i),i=1,ns)
  read(inp1) (psim(i),i=1,ns-1)
  read(inp1) ((xin(i,j),i=1,nt1+2),j=1,ns)
  read(inp1) ((zin(i,j),i=1,nt1+2),j=1,ns)
 
  close(inp1)
 
 
  ! reverse theta orientation if necessary
 
  x0 = sum(xin(1:nt1-1,1))/real(nt1-1)
  z0 = sum(zin(1:nt1-1,1))/real(nt1-1)
  x1 = xin(1, ns)
  z1 = zin(1, ns)
  x2 = xin(2, ns)
  z2 = zin(2, ns)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
  if( it_orientation*det < 0.0_i2mex_r8 ) then
     print *,'keep theta orientation'
     do i = 1, nt1
        x(i, 1:ns) = REAL(xin(i+2, 1:ns), i2mex_r8)
        z(i, 1:ns) = REAL(zin(i+2, 1:ns), i2mex_r8)
     enddo
  else
     print *,'reverse theta orientation'
     do i = 1, nt1
        x(i, 1:ns) = REAL(xin(nt1+2-i+1, 1:ns), i2mex_r8)
        z(i, 1:ns) = REAL(zin(nt1+2-i+1, 1:ns), i2mex_r8)
     enddo
  endif
 
  ! assume uniform mesh in theta
 
  the=i2mex_twopi_r8*(/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)
 
  ! integrate pp & gp to get p, resp. g
 
  pedge = pm(ns-1)  + &
       & (psi(ns)-psim(ns-1))*0.5_i2mex_r8*(pp(ns)+pp(ns-1))! from input file
  call i2mex_integrate1d(ns, psi, pp, ns, psi, p)
  p = p - p(ns) + pedge
  gedge = g(ns)
  call i2mex_integrate1d(ns, psi, gp, ns, psi, g)
  g = g - g(ns) + gedge
 
  psi = abs(psi - psi(1))
 
  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & nt1, ns, the, psi, g, x, z, ier)
  endif
 
  call i2mex_init('INP1_binary_'//filename, &
       & nt1, ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  if ( iok/=0 ) ier = 3
 
  deallocate(zdum, stat=iok)
  deallocate(p, stat=iok)
  deallocate(pm, stat=iok)
  deallocate(pp, stat=iok)
  deallocate(q, stat=iok)
  deallocate(qp, stat=iok)
  deallocate(g, stat=iok)
  deallocate(gp, stat=iok)
  deallocate(the, stat=iok)
  deallocate(psi, stat=iok)
  deallocate(psim, stat=iok)
  deallocate(xin, stat=iok)
  deallocate(zin, stat=iok)
  deallocate(x, stat=iok)
  deallocate(z, stat=iok)
 
  if( iok/=0 ) ier = 4
 
 
end subroutine

subroutine i2mex_fromLDPBinary(filename, it_orientation, ier)
  !
  ! read data from L. Don Pearlstein's inverse equilibrium code.
  !
  use i2mex_mod
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)

  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier

  integer nt1, ns, iok, inp1, i, j
  real(r8) :: pedge, gedge
  real(r8), dimension(:), allocatable :: the, psi, p, g, q
  real(r8), dimension(:,:), allocatable :: x, z, xin, zin
  real(r8) :: x0, x1, x2, z0, z1, z2, det

  ier = 0

  inp1 = 1

  OPEN(inp1, FILE=filename, STATUS='old',                                 &
       &            FORM='unformatted', IOSTAT=iok)
  IF(iok /= 0) THEN
     print*,'**ERROR** could not find file ', filename
     ier = 1
     return
  ENDIF
  REWIND(inp1)

  READ(inp1) ns,nt1

  allocate(psi(ns), stat=iok)
  allocate(p(ns), stat=iok)
  allocate(q(ns), stat=iok)
  allocate(g(ns), stat=iok)
  allocate(the(nt1), stat=iok)
  allocate(x(nt1, ns), stat=iok)
  allocate(z(nt1, ns), stat=iok)
  allocate(xin(nt1, ns), stat=iok)
  allocate(zin(nt1, ns), stat=iok)

  if( iok/=0 ) then
    print*,'***ERROR** memory allocation '
    ier = 5
    return
  endif
!-----------------------------------------------------------------------
!     read binary data.
!-----------------------------------------------------------------------
  READ(inp1)psi
  READ(inp1)g
  READ(inp1)p
  READ(inp1)q
  READ(inp1)xin
  READ(inp1)zin
  close(inp1)
  p = 4.e-7_r8*i2mex_pi_r8*p

  psi = abs(psi - psi(1))
  g = abs(g)

  ! reverse theta orientation if necessary

  x0 = sum(x(1:nt1-1,1))/real(nt1-1)
  z0 = sum(z(1:nt1-1,1))/real(nt1-1)
  x1 = x(1, ns)
  z1 = z(1, ns)
  x2 = x(2, ns)
  z2 = z(2, ns)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)

  if( it_orientation*det < 0.0_i2mex_r8 ) then
     print *,'keep theta orientation'
     do i = 1, nt1
        x(i, 1:ns) = REAL(xin(i, 1:ns), i2mex_r8)
        z(i, 1:ns) = REAL(zin(i, 1:ns), i2mex_r8)
     enddo
  else
     print *,'reverse theta orientation'
     do i = 1, nt1
        x(i, 1:ns) = REAL(xin(nt1-i+1, 1:ns), i2mex_r8)
        z(i, 1:ns) = REAL(zin(nt1-i+1, 1:ns), i2mex_r8)
     enddo
  endif

  ! assume uniform mesh in theta

  the=i2mex_twopi_r8*(/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)

  ! This equilibrium file type has equal angle for a poloidal distribution
  ! it must be remapped
  call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
       & nt1, ns, the, psi, g, x, z, ier)

  call i2mex_init('INP2_binary_'//filename, &
       & nt1, ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)

  call  eqi_jacheck_maxvar_set(0.0001_r8)

  if ( iok/=0 ) ier = 3

  deallocate(psi, stat=iok)
  deallocate(p, stat=iok)
  deallocate(q, stat=iok)
  deallocate(g, stat=iok)
  deallocate(the, stat=iok)
  deallocate(x, stat=iok)
  deallocate(z, stat=iok)
  deallocate(xin, stat=iok)
  deallocate(zin, stat=iok)

  if( iok/=0 ) ier = 4

end subroutine
 
subroutine i2mex_fromMDSPlus(path, time, nt1, ns, it_orientation, ier)
 ! Access equilibrium profile data from MDSPlus tree or local directory.
 
 !
 !    MDS+ path syntax:
 !                 MDS+:<server-name>:<tree-name>(<shot-number>) or
 !                 MDS+:<server-name>:<tree-name>(<tok.yy>,<runid>)
 !
 !    Unix file path syntax:  <path>/<runid> or just <runid>
 !
 !    VMS file path syntax:   <disk>:[<dir>]<runid> or just <runid>
 
 
  use i2mex_mod
  use ezspline
  use ezspline_obj
  implicit none
 
  character*(*), intent(in) :: path ! MDSPlus path (see above)
  real(i2mex_r8), intent(inout) :: time ! time slice in sec, will be corrected
                                      ! if outside shot interval
  integer, intent(in) :: nt1 ! no of poloidal rays + 1
  integer, intent(in) :: ns ! no of radial nodes
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
  integer iok, iwarn, igot, i
  real(i2mex_r8) :: stime, ftime, dtime
 
  real(i2mex_r8), dimension(:,:), allocatable :: x, z
  real(i2mex_r8), dimension(:), allocatable :: tnew, snew, psinew, pnew, gnew, qnew
  real(i2mex_r8), dimension(:,:), allocatable :: xnew, znew
  type(ezspline1) :: f_spl1
  type(ezspline2) :: f_spl2
 
  character*7 string_time
 
  ! *NOTE* default Reals!!!
  real, dimension(:), allocatable :: the, psi, p, g, q, rho
  real, dimension(:,:), allocatable :: xin, zin
  real toroidalFlux, plasmaCurrent
 
  real(i2mex_r8) :: x0, x1, x2, z0, z1, z2, det
 
  ! trread stuff
  CHARACTER*32 label
  CHARACTER*16 units
  INTEGER imulti
  INTEGER istype
  INTEGER IRANK
  INTEGER IDIMS(8)
  CHARACTER*10 ZXABB(8)
 
  integer ns_in
 
  ier = 0
 
  !call trx_msgs(6, 'i2mex')
 
  call trx_connect(trim(path), iok)
  if(iok /= 0) then
     ier = 97
     return
  endif
 
  call trx_tlims(stime,ftime,iok)
  if(stime > time) time = stime
  if(ftime < time) time = ftime
 
 
  dtime = 0.01_i2mex_r8
  call trx_time(time,dtime,iwarn,iok)
  if(iwarn /= 0) then
     ier = 98
  endif
  if(iok /=0 ) then
     ier = 99
     return
  endif
 
  call trx_ready('trx_mhd',iok)
  if(iok /= 0) then
     ier = 103
     return
  endif
 
  ! get the number of radial surfaces
 
  call rplabel('XB',label,units,imulti,istype)
  if(istype.le.0) then
     ier = 100
     return
  endif
!
  call rpdims(istype,irank,idims,zxabb,iok)
  if(iok.ne.0) then
     ier = 101
     return
  endif
 
  ns_in=idims(1)+1
  write(*,*) ' Radial dimension ns_in = ', ns_in
 
  allocate(the(nt1))
  allocate(rho(ns_in), psi(ns_in), p(ns_in), q(ns_in), g(ns_in))
  allocate(x(nt1, ns_in), z(nt1, ns_in))
  allocate(xin(nt1, ns_in), zin(nt1, ns_in))
 
  the = i2mex_twopi_r4* (/ (real(i-1)/real(nt1-1), i=1, nt1) /)
 
  call t1mhdeq(real(time), &
       & real(dtime) ,  &
       & ns_in, nt1, igot,'MKS',the, &
       & xin, zin, rho, psi, p, q, g, toroidalFlux, plasmaCurrent, iok)
  write(*, '(a,F15.3,a)') ' toroidalFlux  (from t1mhdeq) ', toroidalFlux,' [Wb]'
  write(*, '(a,F15.3,a)') ' plasmaCurrent (from t1mhdeq) ', plasmaCurrent,' [A]'
 
!!$  print *,'psi, p, q, g'
!!$  do i=1, ns_in
!!$     print *, psi(i), p(i), q(i), g(i)
!!$  enddo
 
  if(iok /=0) then
     ier = 102
     return
  endif
 
  ! reverse theta orientation if necessary
 
  x0 = sum(xin(1:nt1-1,1))/real(nt1-1)
  z0 = sum(zin(1:nt1-1,1))/real(nt1-1)
  x1 = xin(1, ns_in)
  z1 = zin(1, ns_in)
  x2 = xin(2, ns_in)
  z2 = zin(2, ns_in)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
  if( it_orientation*det < 0.0_i2mex_r8 ) then
     print *,'keep theta orientation'
     do i = 1, nt1
        x(i, 1:ns_in) = REAL(xin(i, 1:ns_in), i2mex_r8)
        z(i, 1:ns_in) = REAL(zin(i, 1:ns_in), i2mex_r8)
     enddo
  else
     print *,'reverse theta orientation'
     do i = 1, nt1
        x(i, 1:ns_in) = REAL(xin(nt1-i+1, 1:ns_in), i2mex_r8)
        z(i, 1:ns_in) = REAL(zin(nt1-i+1, 1:ns_in), i2mex_r8)
     enddo
  endif
 
  ! interpolate on finer grid and apply precision/units conversions
 
  ALLOCATE(snew(ns))
  ALLOCATE(tnew(nt1))
  ALLOCATE(psinew(ns))
  ALLOCATE(pnew(ns))
  ALLOCATE(qnew(ns))
  ALLOCATE(gnew(ns))
  ALLOCATE(xnew(nt1, ns))
  ALLOCATE(znew(nt1, ns))
  tnew = real(the, i2mex_r8)
  snew = (/ (real(i-1,i2mex_r8)/real(ns-1,i2mex_r8), i=1, ns) /)
 
  ! psi
  call EZspline_init(f_spl1, ns_in, (/0,0/), iok)
  call EZspline_setup(f_spl1, real(psi, i2mex_r8), iok)
  call EZspline_interp(f_spl1, ns, snew, psinew, iok)
  call EZspline_free(f_spl1, iok)
 
  ! pressure
  call EZspline_init(f_spl1, ns_in, (/0,0/), iok)
  call EZspline_setup(f_spl1, 4.0e-7_i2mex_r8*i2mex_pi_r8 * real(p, i2mex_r8), iok)
  call EZspline_interp(f_spl1, ns, snew, pnew, iok)
  call EZspline_free(f_spl1, iok)
 
  ! q
  call EZspline_init(f_spl1, ns_in, (/0,0/), iok)
  call EZspline_setup(f_spl1, real(q, i2mex_r8), iok)
  call EZspline_interp(f_spl1, ns, snew, qnew, iok)
  call EZspline_free(f_spl1, iok)
 
  ! g
  call EZspline_init(f_spl1, ns_in, (/0,0/), iok)
  call EZspline_setup(f_spl1, real(g, i2mex_r8), iok)
  call EZspline_interp(f_spl1, ns, snew, gnew, iok)
  call EZspline_free(f_spl1, iok)
 
  ! x
  call EZspline_init(f_spl2, nt1, ns_in, (/-1,-1/), (/0,0/), iok)
  call EZspline_setup(f_spl2, x, iok)
  call EZspline_interp(f_spl2, nt1, ns, tnew, snew, xnew, iok)
  call EZspline_free(f_spl2, iok)
 
  ! z
  call EZspline_init(f_spl2, nt1, ns_in, (/-1,-1/), (/0,0/), iok)
  call EZspline_setup(f_spl2, z, iok)
  call EZspline_interp(f_spl2, nt1, ns, tnew, snew, znew, iok)
  call EZspline_free(f_spl2, iok)
 
 
 
!!$  ! convert to i2mex units, r8 precision and initialize...
!!$  ! psi -> psi/2*pi [Wb/rad]
!!$  ! p -> 4e-7 *p [T^2]
!!$
!!$  call i2mex_init(trim(path), &
!!$       & nt1, ns, &
!!$       & real(the, i2mex_r8), &
!!$       & real(psi, i2mex_r8)/i2mex_twopi_r8, &
!!$       & 4.0e-7_i2mex_r8*i2mex_pi_r8 * real(p, i2mex_r8), &
!!$       & real(g, i2mex_r8), &
!!$       & real(q, i2mex_r8), &
!!$       & x, z, &
!!$       & iok)
!!$  call i2mex_error(iok)
 
 
  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & nt1, ns, tnew, psinew, gnew, xnew, znew, ier)
  endif
 
!!$  print *,'psinew, pnew, qnew, gnew'
!!$  do i=1, ns
!!$     write(*,'(4F12.6)') psinew(i), pnew(i), qnew(i), gnew(i)
!!$  enddo
 
  write(string_time, '(f7.3)') time
  call i2mex_init(trim(path)//'_t='//trim(adjustl(string_time)), &
       & nt1, ns, &
       & tnew, &
       & psinew, &
       & pnew, &
       & gnew, &
       & qnew, &
       & xnew, znew, &
       & iok)
  call i2mex_error(iok)
 
  !  dmc-- I think with TRANSP runs it is better to update the q profile
  !  from Psi(s) with a numerical integration-- done here:

  call eqm_gen_q('q_i2mex',i2mex_o%id_q,iok)
  if(iok.ne.0) iok=6
  call i2mex_error(iok)

  DEALLOCATE(snew)
  DEALLOCATE(tnew)
  DEALLOCATE(psinew)
  DEALLOCATE(pnew)
  DEALLOCATE(qnew)
  DEALLOCATE(gnew)
  DEALLOCATE(xnew)
  DEALLOCATE(znew)
 
  deallocate(the)
  deallocate(rho, psi, p, q, g)
  deallocate(x, z)
  deallocate(xin, zin)
 
end subroutine i2mex_fromMDSPlus
 
subroutine i2mex_fromMenard(filename, it_orientation, ier)
  !
  ! read equilibrium file in psipqgRZ format (Menard)
  !
  use ezcdf
  use i2mex_mod
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
  integer nt1, ns, iok, i, j, iu, dimlens(3)
  real(r8), dimension(:), allocatable :: p, g, q, &
       & the, psi
  real(r8), dimension(:,:), allocatable :: x, z, xin, zin
  character*3 :: xtype
  real(r8) :: x0, x1, x2, z0, z1, z2, det
 
  ier = 0
 
  call ezcdf_open(iu, filename, 'r', iok)
  IF(iok /= 0) THEN
     print*,'**ERROR** could not open ', filename
     ier = 126
     return
  ENDIF
 
  ! read dimensions
 
  call cdfInqVar(iu, 'psi', dimlens, xtype, iok)
  ns = dimlens(1)
  call cdfInqVar(iu, 'R', dimlens, xtype, iok)
  nt1 = dimlens(1)
 
  allocate(p(ns), stat=iok)
  allocate(q(ns), stat=iok)
  allocate(g(ns), stat=iok)
  allocate(the(nt1), stat=iok)
  allocate(psi(ns), stat=iok)
  allocate(x(nt1, ns), stat=iok)
  allocate(z(nt1, ns), stat=iok)
  allocate(xin(nt1, ns), stat=iok)
  allocate(zin(nt1, ns), stat=iok)
 
  if( iok/=0 ) then
     print*,'***ERROR** memory allocation '
     ier = 127
     return
  endif
 
 
  ! read profiles
 
  call cdfGetVar(iu, 'psi', psi, iok)
  psi = abs(psi - psi(1))
  call cdfGetVar(iu, 'p', p, iok)
  p = 2.e-7_r8*i2mex_twopi_r8 * p
  call cdfGetVar(iu, 'q', q, iok)
  call cdfGetVar(iu, 'g', g, iok)
  g = abs(g)
  if(iok/=0) then
     ier = 128
     return
  endif
 
  ! read grid
 
  call cdfGetVar(iu, 'R', xin, iok)
  call cdfGetVar(iu, 'Z', zin, iok)
  if(iok/=0) then
     ier = 129
     return
  endif
 
  call cdfCls(iu)
 
 
  ! reverse theta orientation if necessary
 
  x0 = sum(xin(1:nt1-1,1))/real(nt1-1)
  z0 = sum(zin(1:nt1-1,1))/real(nt1-1)
  x1 = xin(1, ns)
  z1 = zin(1, ns)
  x2 = xin(2, ns)
  z2 = zin(2, ns)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
  if( it_orientation*det < 0.0_i2mex_r8 ) then
     print *,'keep theta orientation'
     do i = 1, nt1
        x(i, 1:ns) = REAL(xin(i, 1:ns), i2mex_r8)
        z(i, 1:ns) = REAL(zin(i, 1:ns), i2mex_r8)
     enddo
  else
     print *,'reverse theta orientation'
     do i = 1, nt1
        x(i, 1:ns) = REAL(xin(nt1-i+1, 1:ns), i2mex_r8)
        z(i, 1:ns) = REAL(zin(nt1-i+1, 1:ns), i2mex_r8)
     enddo
  endif
 
 
  ! assume uniform mesh in theta
 
  the=i2mex_twopi_r8*(/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)
 
  ! initialize
 
  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & nt1, ns, the, psi, g, x, z, ier)
  endif
 
  call i2mex_init('Menard_input_'//filename, &
       & nt1, ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  if ( iok/=0 ) ier = 130
 
  deallocate(p, stat=iok)
  deallocate(q, stat=iok)
  deallocate(g, stat=iok)
  deallocate(the, stat=iok)
  deallocate(psi, stat=iok)
  deallocate(xin, stat=iok)
  deallocate(zin, stat=iok)
  deallocate(x, stat=iok)
  deallocate(z, stat=iok)
 
  if( iok/=0 ) ier = 131
 
 
end subroutine i2mex_fromMenard


subroutine i2mex_fromCEO(filename, it_orientation, ier)
  !
  ! read equilibrium file in ceo format (Pletzer)
  !
  use ezcdf
  use i2mex_mod
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  character*(*), intent(in) :: filename
  integer, intent(in) :: it_orientation
  integer, intent(out) :: ier
 
  integer nt1, ns, iok, i, j, iu, dimlens(3)
  real(r8), dimension(:), allocatable :: p, g, q, &
       & the, psi
  real(r8), dimension(:,:), allocatable :: x, z, xin, zin
  character*3 :: xtype
  real(r8) :: x0, x1, x2, z0, z1, z2, det
  character(20) :: pName, qName, gName, psibigName, rName, zName
 
  ier = 0
 
  call ezcdf_open(iu, filename, 'r', iok)
  IF(iok /= 0) THEN
     print*,'**ERROR** could not open ', filename
     ier = 126
     return
  ENDIF
 

  if( iok/=0 ) then
     print*,'***ERROR** memory allocation '
     ier = 127
     return
  endif
 
  ! look for profile names
  psibigName = 'psibig'
  call cdfInqVar(iu, psibigName, dimlens, xtype, iok)
  if(iok/=0) then
     ! try
     psibigName = 'PsiBig'
     call cdfInqVar(iu, psibigName, dimlens, xtype, iok)
     if(iok/=0) then
        print *,' i2mex::input error cannot find psibig???'
        ier = 164
     endif
  endif

  pName = 'p'
  call cdfInqVar(iu, pName, dimlens, xtype, iok)
  if(iok/=0) then
     ! try
     pName = 'pa'
     call cdfInqVar(iu, pName, dimlens, xtype, iok)
     if(iok/=0) then
        print *,' i2mex::input error cannot find p???'
        ier = 165
     endif
  endif
  ! read in radial dimensions
  ns = dimlens(1)

  qName = 'q'
  call cdfInqVar(iu, qName, dimlens, xtype, iok)
  if(iok/=0) then
     ! try
     qName = 'qa'
     call cdfInqVar(iu, qName, dimlens, xtype, iok)
     if(iok/=0) then
        print *,' i2mex::input error cannot find q???'
        ier = 166
     endif
  endif

  gName = 'g'
  call cdfInqVar(iu, gName, dimlens, xtype, iok)
  if(iok/=0) then
     ! try
     gName = 'ga'
     call cdfInqVar(iu, gName, dimlens, xtype, iok)
     if(iok/=0) then
        print *,' i2mex::input error cannot find g???'
        ier = 167
     endif
  endif

  rName = 'r'
  call cdfInqVar(iu, rName, dimlens, xtype, iok)
  if(iok/=0) then
     ! try
     rName = 'xa'
     call cdfInqVar(iu, rName, dimlens, xtype, iok)
     if(iok/=0) then
        print *,' i2mex::input error cannot find r???'
        ier = 168
     endif
  endif
  ! read in poloidal dimension
  nt1 = dimlens(1)

  zName = 'z'
  call cdfInqVar(iu, zName, dimlens, xtype, iok)
  if(iok/=0) then
     ! try
     zName = 'za'
     call cdfInqVar(iu, zName, dimlens, xtype, iok)
     if(iok/=0) then
        print *,' i2mex::input error cannot find z???'
        ier = 169
     endif
  endif

  ! allocate
  allocate(p(ns), stat=iok)
  allocate(q(ns), stat=iok)
  allocate(g(ns), stat=iok)
  allocate(the(nt1), stat=iok)
  allocate(psi(ns), stat=iok)
  allocate(x(nt1, ns), stat=iok)
  allocate(z(nt1, ns), stat=iok)
  allocate(xin(nt1, ns), stat=iok)
  allocate(zin(nt1, ns), stat=iok)
 
  ! read profiles
 
  call cdfGetVar(iu, psibigName, psi, iok)
  psi = abs(psi - psi(1))
  psi = psi/i2mex_twopi_r8
  call cdfGetVar(iu, pName, p, iok)
  call cdfGetVar(iu, qName, q, iok)
  call cdfGetVar(iu, gName, g, iok)
  g = abs(g)
  if(iok/=0) then
     ier = 128
     return
  endif
 
  ! read grid
 
  call cdfGetVar(iu, rName, xin, iok)
  call cdfGetVar(iu, zName, zin, iok)
  if(iok/=0) then
     ier = 129
     return
  endif
 
  call cdfCls(iu)
 
 
  ! reverse theta orientation if necessary
 
  x0 = sum(xin(1:nt1-1,1))/real(nt1-1)
  z0 = sum(zin(1:nt1-1,1))/real(nt1-1)
  x1 = xin(1, ns)
  z1 = zin(1, ns)
  x2 = xin(2, ns)
  z2 = zin(2, ns)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
  if( it_orientation*det < 0.0_i2mex_r8 ) then
     print *,'keep theta orientation'
     do i = 1, nt1
        x(i, 1:ns) = REAL(xin(i, 1:ns), i2mex_r8)
        z(i, 1:ns) = REAL(zin(i, 1:ns), i2mex_r8)
     enddo
  else
     print *,'reverse theta orientation'
     do i = 1, nt1
        x(i, 1:ns) = REAL(xin(nt1-i+1, 1:ns), i2mex_r8)
        z(i, 1:ns) = REAL(zin(nt1-i+1, 1:ns), i2mex_r8)
     enddo
  endif
 
 
  ! assume uniform mesh in theta
 
  the=i2mex_twopi_r8*(/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)
 
  ! initialize
 
  ! to remap poloidally set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  if(i2mex_remap/=0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & nt1, ns, the, psi, g, x, z, ier)
  endif
 
  call i2mex_init('Ceo_input_'//filename, &
       & nt1, ns, the, psi, p, g, q, x, z, iok)
  call i2mex_error(iok)
 
  if ( iok/=0 ) ier = 170
 
  deallocate(p, stat=iok)
  deallocate(q, stat=iok)
  deallocate(g, stat=iok)
  deallocate(the, stat=iok)
  deallocate(psi, stat=iok)
  deallocate(xin, stat=iok)
  deallocate(zin, stat=iok)
  deallocate(x, stat=iok)
  deallocate(z, stat=iok)
 
  if( iok/=0 ) ier = 170
 
 
end subroutine i2mex_fromCEO
