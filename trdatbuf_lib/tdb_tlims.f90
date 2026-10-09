subroutine tdb_tlims(d,zd_tinit,zd_ftime)
  use trdatbuf_obj
  implicit NONE
  !
  !  return time range covered by trdatbuf data
  !
  type (trdatbuf) :: d
  real*8, intent(out) :: zd_tinit  ! start time
  real*8, intent(out) :: zd_ftime  ! stop time

  zd_tinit = d%datbuf(1)
  zd_ftime = d%datbuf(2)

end subroutine tdb_tlims

subroutine tdb_tlim_ok(d,zd_ftime_ok)
  use trdatbuf_obj
  implicit NONE
  !
  !  return time beyond which a code failure should be ignored
  !  (allows experimentalists to set off time after disruption,
  !  without leading to a crash requiring manual intervention 
  !  for run to completion -- FTIME_OK in the TRANSP namelist).
  !
  type (trdatbuf) :: d
  real*8, intent(out) :: zd_ftime_ok  ! stop time

  zd_ftime_ok = d%datbuf(10)

end subroutine tdb_tlim_ok
