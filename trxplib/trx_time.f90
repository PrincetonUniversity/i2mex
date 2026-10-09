subroutine trx_time(ztime,zdelta,iwarn,ierr)
  use trx_module
  implicit NONE
!
!   set the time and delta_t for fetching/averaging data
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  real*8, intent(IN) :: ztime      ! time of interest, seconds
  real*8, intent(IN) :: zdelta     ! delta(time), +/- for averaging, seconds
  integer, intent(OUT) :: iwarn    ! warn if ztime-zdelta or ztime+zdelta
                                   ! ...is outside the limit
 
  integer, intent(OUT) :: ierr     ! completion code: 0=OK
!
!-----------------------------
!
  integer lunzer
!
  real ztime_r4,zdelta_r4
!
!-----------------------------
!
  ierr=0
  iwarn=0
!
  ztime_r4=ztime
  zdelta_r4=zdelta
!
  if(nsctime.eq.0) then
     write(lunzer(0),*) ' ?trx_time:  call trx_connect first!'
     ierr=1
     return
  endif
!
!  set module variables...
!
  itset=1
!
  time0=max(tmin,min(tmax,ztime_r4))
  if(time0.ne.ztime_r4) iwarn=1
!
  delta_t=abs(zdelta_r4)
!
  if((time0+delta_t).gt.tmax) iwarn=1
  if((time0-delta_t).lt.tmin) iwarn=1
!
  return
end subroutine trx_time
