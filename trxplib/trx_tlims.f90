subroutine trx_tlims(stime,ftime,ierr)
 
  use trx_module
  implicit NONE
 
!
!  return the start time and stop time, seconds, of the current run.
!  return ierr.ne.0 if no run is connected.
!
  real*8, intent(out) :: stime      ! start time of run (seconds)
  real*8, intent(out) :: ftime      ! stop time of run (seconds)
!
  integer ierr                      ! completion code, 0=OK
!
!----------------------------------
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
  integer lunzer
!----------------------------------
!
  stime=0.0_R8
  ftime=0.0_R8
!
  if(nsurf.eq.0) then
     write(lunzer(0),*) ' ?trx_tlims: nsurf=0:  use trx_connect first!'
     ierr=1
     return
  endif
!
  ierr=0
  stime=tmin
  ftime=tmax
!
  return
end subroutine trx_tlims
