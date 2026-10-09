subroutine makedens_gneg(nx,x,f,fout,gfac,inside,iedge,imin_patch,itty,istat)

  ! DMC Dec 2010

  ! some RF codes what density gradients near the edge to always be
  ! negative (i.e. always increasing towards the inside).

  ! this routine patches a profile f(1:nx) to enforce this, i.e. if at the
  ! very edge df/dx is forced to have the expected behavior.

  !  *              *
  !    *              *
  !         *   =>      *  *
  !      *  

  ! (the moved point has slightly higher value than the end point)
  ! (if more than one point is moved, the last moved point is replaced with
  ! an average between the preceding point and the first unmoved point).

  ! a derivative value limit is defined:
  !
  !    dlimit = gfac * (f(iedge)-f(inside))/(x(iedge)-x(inside))
  !
  !       with 0.00001 < gfac <= 0.1 expected & enforced
  !
  !    dlimit < 0 is expected; if not the routine returns fout = f
  !    writes a message and sets a status code
  !
  ! the profile fout(1:nx) is f(1:nx) patched so that df/dx in the
  ! edge region is no less than dlimit.  It does this, marching in until
  ! it encounters a point in the original data that has a higher value
  ! than is necessary to satisfy this condition.  This is supposed to
  ! happen before index (imin_patch) is reached.  If not, a status code
  ! is set, fout = f is returned, and a message is printed.

  use iso_c_binding, only: fp => c_double
  implicit none

  !---------------------
  !  arguments

  integer, intent(in) :: nx   ! array size
  real(fp), intent(in) :: x(nx) ! x grid; for j=inside:iedge-1 expect x(j+1)>x(j)
  real(fp), intent(in) :: f(nx) ! input profile

  real(fp), intent(out) :: fout(nx) ! modified profile

  real(fp), intent(in) :: gfac  ! derivative limit factor

  integer, intent(in) :: inside  ! inside reference point
  integer, intent(in) :: iedge   ! edge reference point

  integer, intent(in) :: imin_patch   ! index of innermost modifiable f value

  !  expected:  inside < imin_patch < iedge

  integer, intent(in) :: itty    ! I/O unit for messages
  integer, intent(out) :: istat  ! status code: 0 = OK

  !---------------------
  !  local

  integer :: ix,ix_stop
  real(fp) :: zgfac,dlimit,dx,fnew,wt,fav

  real(fp), parameter :: gmin = 0.00002_fp
  real(fp), parameter :: gmax = 0.2_fp

  !---------------------
  !  error checks

  istat = 0
  fout = f

  !  check gfac value (warning only, if out of range)
  zgfac = max(gmin,min(gmax,gfac))
  if(zgfac.ne.gfac) then
    write(itty,*) ' %makedens_gneg: local value of gfac modified to satisfy constraints:'
    write(itty,*) '  gfac as input: ',gfac,'; gfac as used: ',zgfac
  endif

  !  check reference indices

  if(iedge.le.inside) then
    istat = 1
    write(itty,*) ' ?makedens_gneg: index inputs error, iedge < inside:'
    write(itty,*) '  iedge = ',iedge,' inside = ',inside
    return
  endif

  if(imin_patch.ge.iedge) then
    istat = 1
    write(itty,*) ' ?makedens_gneg: index inputs error, imin_patch >= iedge:'
    write(itty,*) '  iedge = ',iedge,' imin_patch = ',imin_patch
    return
  endif

  if(imin_patch.le.inside) then
    istat = 1
    write(itty,*) ' ?makedens_gneg: index inputs error, imin_patch <= inside:'
    write(itty,*) '  inside = ',inside,' imin_patch = ',imin_patch
    return
  endif

  !  check x order w/in range of indices

  do ix=inside,iedge-1
    if(x(ix+1).le.x(ix)) then
      istat = 2
      exit
    endif
  enddo

  if(istat.gt.0) then
    write(itty,*) ' ?makedens_gneg: x(inside:iedge) grid not strict ascending:'
    write(itty,*) '  ',x(inside:iedge)
    return
  endif

  ! check f(inside) > f(iedge)

  if(f(inside).le.f(iedge)) then
    istat = 3
    write(itty,*) ' ?makedens_gneg:  f(inside) <= f(iedge):'
    write(itty,*) '  f(inside) = ',f(inside),'; f(iedge) = ',f(iedge)
    return
  endif

  !-------------------------------------------------
  !  OK all checks passed
  !-------------------------------------------------

  dlimit = zgfac * (f(iedge)-f(inside))/(x(iedge)-x(inside))

  ix_stop = imin_patch-1

  do ix=iedge-1,imin_patch-1,-1
    dx = x(ix+1)-x(ix)
    fnew = fout(ix+1) - dx*dlimit
    if(fout(ix).lt.fnew) then
      fout(ix)=fnew
    else
      ! OK: found a high enough value that we can stop
      ix_stop = ix+1
      if(ix_stop.lt.iedge-1) then
        wt = (x(ix_stop)-x(ix_stop-1))/(x(ix_stop+1)-x(ix_stop-1))
        fav = fout(ix_stop-1)*(1.0_fp-wt) + fout(ix_stop+1)*wt
        fout(ix_stop) = (fout(ix_stop) + fav)*0.5_fp
      endif
      exit
    endif
  enddo

  if(ix_stop.lt.imin_patch) then
    istat = 99
    write(itty,*) ' ?makedens_gneg: derivative patch failed, patch region too large.'
    fout = f
    return
  endif

  return
end subroutine makedens_gneg
