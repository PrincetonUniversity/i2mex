subroutine tdb_get_rzgrids(d,inumR,Rgrid,inumZ,Zgrid,ier)

  !  get R and Z grids -- but, passed argument grid sizes must match or
  !  an error code is set and a message is printed.

  !  see sister subroutine "tdb_get_rzsizes"...

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d

  integer,intent(in) :: inumR,inumZ   ! the grid sizes
  real*8 :: Rgrid(inumR),Zgrid(inumZ) ! the grids, set if sizes match
  integer, intent(out) :: ier         ! status code: 0=OK

  !------------------------------------
  integer :: nonlin,lunmsg_tdb
  !------------------------------------

  ier = 0

  nonlin = lunmsg_tdb(0)

  if(inumR.ne.d%nRpsi) then
     ier = ier + 1
     write(nonlin,*) ' ?tdb_get_RZgrids: R grid size mismatch.'
     write(nonlin,*) '  Correct size is: ',d%nRpsi,'; passed size is: ',inumR
  endif

  if(inumZ.ne.d%nZpsi) then
     ier = ier + 1
     write(nonlin,*) ' ?tdb_get_RZgrids: Z grid size mismatch.'
     write(nonlin,*) '  Correct size is: ',d%nZpsi,'; passed size is: ',inumZ
  endif

  if(ier.ne.0) return

  if(inumR.gt.0) then
     Rgrid = d%datbuf(d%lRpsi:d%lRpsi+d%nRpsi-1)
  endif

  if(inumZ.gt.0) then
     Zgrid = d%datbuf(d%lZpsi:d%lZpsi+d%nZpsi-1)
  endif

end subroutine tdb_get_rzgrids
