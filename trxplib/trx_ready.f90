subroutine trx_ready(zsub,ierr)
!
!  check if we are connected to a run; if not, set error flag
!  and write message
!
!  also check if desired time has been set
!
  use trx_module
  implicit NONE
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  character(*), intent(in) :: zsub    ! calling subroutine name for msgs
  integer, intent(out)     :: ierr    ! completion code, 0 = no error
!
  integer lunzer
!----------------------
!
  if(nsurf.eq.0) then
     write(lunzer(0),*) ' ?',zsub,': nsurf=0:  use trx_connect first!'
     ierr=1
     return
  endif
!
  if(itset.eq.0) then
     write(lunzer(0),*) ' ?',zsub,': itset=0:  use trx_time first!'
     ierr=1
     return
  endif
!
  ierr=0
  return
!
end subroutine trx_ready
