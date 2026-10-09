subroutine tdb_get_rzsizes(d,inumR,inumZ)

  !  get sizes of R and Z grids to go with Psi(R,Z) data vs. time

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d

  integer,intent(out) :: inumR,inumZ  ! the grid sizes, or 0 if none found.

  !------------------------------------

  inumR = d%nRpsi
  inumZ = d%nZpsi

end subroutine tdb_get_rzsizes
