integer function tdb_nchan_find(d,z2char)

  ! return number of channels for indicated data type, found in trdatbuf buffer

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d                 ! data buffer object
  character*(*), intent(in) :: z2char  ! 2-character code...

  !  a copy of z2char is converted to uppercase before testing
  !  if z2char is uppercase, then,

  !    z2char.eq.'NB' means -- requested information is # of neutral beams
  !
  !    z2char.eq.'EC' means -- requested information is # of ECH antennas
  !       electron cyclotron heating / current drive
  !
  !    z2char.eq.'LH' means -- requested information is # of LH antennas
  !       lower hybrid heating / current drive
  !
  !    z2char.eq.'RF' means -- requested information is # of ICRF antennas
  !       ICRF heating / current drive

  !  if z2char is recognized, the number of beams or antennas is returned;
  !  it could be zero.

  !  if z2char is not recognized, the number -1 is returned and a warning
  !  is written to the tdb log...

  character*3 ztest
  integer :: lunmsg_tdb

  !----------------------

  ztest = z2char(1:min(len(z2char),len(ztest)))
  call uupper(ztest)

  if(ztest.eq.'NB') then
     tdb_nchan_find = d%nbdata

  else if(ztest.eq.'EC') then
     tdb_nchan_find = d%nantech_d

  else if(ztest.eq.'LH') then
     tdb_nchan_find = d%nantlh_d

  else if(ztest.eq.'RF') then
     tdb_nchan_find = d%nantich_d

  else
     tdb_nchan_find = -1
     write(lunmsg_tdb(0),*) &
          ' ?? trdatbuf_lib/tdb_nchan_find -- unrecognized designation: ', &
          trim(ztest)
  endif

end function tdb_nchan_find
