logical function ltrx_freebdy(idum)

  use trx_module
  implicit NONE

  !  return .TRUE. if free boundary data are available, o.w. .FALSE.

  integer :: idum    ! dummy argument (not currently used)

  !---------------------------

  ltrx_freebdy = (nRfree.gt.0).AND.(nZfree.gt.0)

end function ltrx_freebdy

subroutine trx_freebdy_RZgridSize(inumR,inumZ)

  use trx_module
  implicit NONE

  !  return the {R,Z} grid sizes for free boundary Psi(R,Z) data
  !  If no such data in the current run, return (0,0)

  integer, intent(out) :: inumR  ! R grid size
  integer, intent(out) :: inumZ  ! Z grid size
  !---------------------------

  inumR = nRfree
  inumZ = nZfree

end subroutine trx_freebdy_RZgridSize

subroutine trx_freebdy_RZminmax(zRmin,zRmax,zZmin,zZmax)

  use trx_module
  implicit NONE

  !  return the {R,Z} grid limits for free boundary Psi(R,Z) data
  !  If no such data in the current run, return (0,0,0,0)

  real*8, intent(out) :: zRmin,zRmax  ! R grid limits (m)
  real*8, intent(out) :: zZmin,zZmax  ! Z grid limits (m)
  !---------------------------

  if(nRfree.eq.0) then
     zRmin=0
     zRmax=0
  else
     zRmin=Rgrid_free(1)
     zRmax=Rgrid_free(nRfree)
  endif

  if(nZfree.eq.0) then
     zZmin=0
     zZmax=0
  else
     zZmin=Zgrid_free(1)
     zZmax=Zgrid_free(nZfree)
  endif

end subroutine trx_freebdy_RZminmax

subroutine trx_freebdy_RZgrids(inRmax,inRgot,zRgrid,inZmax,inZgot,zZgrid,ierr)

  use trx_module
  implicit NONE

  !  return the {R,Z} grids for free boundary data.
  !  If no such data in the current run, return zeroes in (inRgot,inZgot);
  !  ierr=0 is returned even if there is no data.

  !  The error return code ierr is set only if the array sizes provided are
  !  to small to receive the grid data.

  integer, intent(in) :: inRmax       ! R grid passed array size
  integer, intent(out) :: inRgot      ! actual number of data points returned
  real*8 :: zRgrid(inRmax)      ! the grid in zRgrid(1:inRgot)

  integer, intent(in) :: inZmax       ! Z grid passed array size
  integer, intent(out) :: inZgot      ! actual number of data points returned
  real*8 :: zZgrid(inZmax)      ! the grid in zZgrid(1:inZgot)

  integer, intent(out) :: ierr        ! completion code (0=OK)
  !---------------------------
  integer :: lunzer
  !---------------------------

  ierr = 0
  inRgot = 0
  inZgot = 0

  if(nRfree.eq.0) return
  if(nZfree.eq.0) return

  if(nRfree.gt.inRmax) then
     write(lunzer(0),*) ' ?trx_freebdy_RZgrids: array size too small:'
     write(lunzer(0),*) '  array size inRmax=',inRmax,'; need at least: ',nRfree
     ierr=1
  endif

  if(nZfree.gt.inZmax) then
     write(lunzer(0),*) ' ?trx_freebdy_RZgrids: array size too small:'
     write(lunzer(0),*) '  array size inZmax=',inZmax,'; need at least: ',nZfree
     ierr=1
  endif

  if(ierr.ne.0) return

  inRgot = nRfree
  inZgot = nZfree
  zRgrid(1:inRgot) = Rgrid_free
  zZgrid(1:inZgot) = Zgrid_free

end subroutine trx_freebdy_RZgrids
