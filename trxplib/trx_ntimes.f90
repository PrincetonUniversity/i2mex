subroutine trx_ntimes(insc,inpr,ierr)

  !   return the number of times in the run

  use trx_module
  implicit NONE

  integer, intent(OUT) :: insc     ! # of times, SCALAR f(t) data
  integer, intent(OUT) :: inpr     ! # of times, PROFILE f(x,t) data
  integer, intent(OUT) :: ierr     ! completion code: 0=OK

  !----------------------------------
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer lunzer
  !----------------------------------

  if(nsurf.eq.0) then
     write(lunzer(0),*) ' ?trx_ntimes: nsurf=0:  use trx_connect first!'
     ierr=1
     return
  endif

  ierr=0
  insc=nsctime
  inpr=nprtime

  return
end subroutine trx_ntimes
