subroutine tdb_ripple0(d,icoils,zphi)
  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d

  !  information on TF ripple field -- primary component

  integer, intent(out) :: icoils  ! no. of TF coils (TF ripples)
  real*8, intent(out) :: zphi ! shift of 1st coil rel. phi=0 (radians)

  icoils = d%datbuf(6) + 0.0001d0  ! slightly exceeds integer value
  zphi = d%datbuf(7)

end subroutine tdb_ripple0

subroutine tdb_ripple0_2(d,icoils,zphi)
  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d

  !  information on TF ripple field -- secondary component

  integer, intent(out) :: icoils  ! no. of TF coils (TF ripples)
  real*8, intent(out) :: zphi ! shift of 1st coil rel. phi=0 (radians)

  icoils = d%datbuf(12) + 0.0001d0  ! slightly exceeds integer value
  zphi = d%datbuf(13)

end subroutine tdb_ripple0_2
