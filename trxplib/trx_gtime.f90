subroutine trx_gtime(ztime,zdelta_t,ierr)
 
  use trx_module
  implicit NONE
 
!
!   return currently set time of interest
!   ierr is set if time of interest has not been set
!
  real*8 ztime             ! time of interest (seconds) (returned)
  real*8 zdelta_t          ! +/- averaging time (seconds) (returned)
  integer ierr             ! completion coee, 0=OK (returned)
!
!-------------------------------
!
  if(itset.eq.0) then
     ztime=0.0
     zdelta_t=0.0
     ierr=1
  else
     ztime=time0
     zdelta_t=delta_t
     ierr=0
  endif
!
  return
end subroutine trx_gtime
