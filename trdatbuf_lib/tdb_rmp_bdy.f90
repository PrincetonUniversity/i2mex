subroutine tdb_rmp_bdy(d,ztime,zr1,zr2)

  !  find midplane boundary intercept radii at the indicated time

  use trdatbuf_obj
  use tdbsub_uts
  implicit NONE

  type(trdatbuf) :: d          ! data buffer object
  real*8, intent(in) :: ztime  ! time (seconds)

  real*8, intent(out) :: zr1   ! inner intercept major radius (cm)
  real*8, intent(out) :: zr2   ! outer intercept major radius (cm)

  !---------------------------------------------
  integer :: it1,iltime,intime,ilr1,ilr2
  real*8 :: zf1
  !---------------------------------------------

  iltime = d%ltime2
  intime = d%ntime2

  call tdbsub_lookup(d%datbuf(iltime:iltime+intime-1),intime,ztime,it1,zf1)

  ilr1 = d%lrmp_bdy1
  ilr2 = d%lrmp_bdy2

  zr1 = (ONE-zf1)*d%datbuf(ilr1+it1-1) + zf1*d%datbuf(ilr1+it1)
  zr2 = (ONE-zf1)*d%datbuf(ilr2+it1-1) + zf1*d%datbuf(ilr2+it1)

end subroutine tdb_rmp_bdy

subroutine tdb_rmp_bdy_minmax(d,zrmin,zrmax)

  !  find the minimum inner boundary intercept radius and the maximum
  !  outer boundary intercept radius over all available times

  use trdatbuf_obj
  implicit NONE

  type(trdatbuf) :: d          ! data buffer object

  real*8, intent(out) :: zrmin ! inner intercept major radius (cm) -- MINIMUM
  real*8, intent(out) :: zrmax ! outer intercept major radius (cm) -- MAXIMUM

  !---------------------------------------------
  integer :: intime,ilr1,ilr2
  !---------------------------------------------

  intime = d%ntime2
  ilr1 = d%lrmp_bdy1
  ilr2 = d%lrmp_bdy2

  zrmin = minval(d%datbuf(ilr1:ilr1+intime-1))
  zrmax = maxval(d%datbuf(ilr2:ilr2+intime-1))

end subroutine tdb_rmp_bdy_minmax
