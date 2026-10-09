integer function lunmsg_tdb(idum)

  ! return i/o unit number for trdatbuf_lib messages

  use tdb_static
  implicit NONE
  integer :: idum

  lunmsg_tdb = lunmsg_trdatbuf

end function lunmsg_tdb

subroutine tdb_lunmsg_set(ilun_new)

  ! set new i/o unit number for trdatbuf_lib messages

  use tdb_static
  implicit NONE
  integer, intent(in) :: ilun_new

  lunmsg_trdatbuf = ilun_new

end subroutine tdb_lunmsg_set
