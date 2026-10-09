subroutine tdbsub_slookup(zt,int,zt1,it1,zfrac1)

  ! call access w/o module: tdbsub_lookup(...)

  use tdbsub_uts
  implicit NONE

  integer, intent(in) :: int
  real*8, intent(in)  :: zt(int)
  real*8, intent(in) :: zt1
  integer, intent(out) :: it1
  real*8, intent(out) :: zfrac1

  call tdbsub_lookup(zt,int,zt1,it1,zfrac1)

end subroutine tdbsub_slookup
