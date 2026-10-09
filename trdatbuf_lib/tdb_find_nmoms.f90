subroutine tdb_find_nmoms(d,imoms,icirc)

  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d

  integer, intent(out) :: imoms,icirc

  ! return imoms = # of moments (0:imoms) in Fourier expansion
  ! return icirc = 1 if circular flux surfaces (minor & major radius only)
  !                are in use (this is very rare).

  imoms = 1
  icirc = 1

  if(d%ldmmx.gt.0) then
     imoms=d%mmax
     icirc=0
  else
     imoms=d%nmomd
     icirc=0
  endif

end subroutine tdb_find_nmoms
