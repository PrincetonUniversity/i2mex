subroutine tdb_ccwchk(d,i_bccw,i_jccw)
  !
  ! if available in trdatbuf object, set direction of B_phi and J_phi...
  !
  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  !  ccw means "counter-clockwise looking down from above"
  integer, intent(out) :: i_bccw  ! 1: ccw, -1: cw, 0: unspecified (B_phi)
  integer, intent(out) :: i_jccw  ! 1: ccw, -1: cw, 0: unspecified (J_phi)

  i_bccw=0
  i_jccw=0

  if(abs(d%DATBUF(3)).gt.0.1) then
     if(d%DATBUF(3).gt.0.0) then
        i_bccw=1
     else
        i_bccw=-1
     endif
  endif
  
  if(abs(d%DATBUF(4)).gt.0.1) then
     if(d%DATBUF(4).gt.0.0) then
        i_jccw=1
     else
        i_jccw=-1
     endif
  endif

end subroutine tdb_ccwchk
