subroutine tdb_ebmax(d,zebmax)
  use trdatbuf_obj
  implicit NONE
  !
  !  return the maximum beam energy found in the data
  !
  type (trdatbuf) :: d
  real*8, intent(out) :: zebmax  ! max beam energy

  if(d%nsc.eq.10) then
     zebmax = 0.0d0        ! old trdat files do not have this...
  else
     zebmax = d%datbuf(11)
  endif

end subroutine tdb_ebmax
