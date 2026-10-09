subroutine tdb_izminmax(d,izmin,izmax)

  !  if d%LDATZIM > 0 & d%NTIME1 > 0:
  !     return izmin and izmax corresponding to min and max Zimp
  !     as stored in d%datbuf(d%LDATZIM:d%LDATZIM+d%NTIME1-1)

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  integer :: izmin,izmax

  !-----------------

  real*8 :: zmin,zmax,ztest
  integer :: ii,indx

  !---------------------------

  izmin=0
  izmax=0

  if(d%ldatzim.eq.0) return
  if(d%ntime1.eq.0) return

  zmin=d%datbuf(d%ldatzim)
  zmax=zmin

  do ii=1,d%ntime1
     indx = d%ldatzim + ii - 1
     zmin=min(zmin,d%datbuf(indx))
     zmax=max(zmax,d%datbuf(indx))
  enddo

  izmin = zmin+0.5
  izmax = zmax
  ztest = izmax

  if(ztest.lt.zmax) izmax=izmax+1

end subroutine tdb_izminmax
