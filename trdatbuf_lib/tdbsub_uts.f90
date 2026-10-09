module tdbsub_uts

  include 'trdatbuf_constants.incl'

  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE = 1.0d0
  real*8, parameter :: HALF = 0.5d0
  real*8, parameter :: EPSINV = 1.0d34

  contains
    subroutine tdbsub_lookup(zt,int,zt1,it1,zfrac1)
      ! in ordered time vector zt(1:int) find zone containing zt1
      ! return index and offset fraction in it1 and zfrac1

      ! **zt ordered & zt1 btw zt(1) and zt(int) assumed w/o checking!!

      implicit NONE
      integer, intent(in) :: int
      real*8, intent(in)  :: zt(int)
      real*8, intent(in) :: zt1
      integer, intent(out) :: it1
      real*8, intent(out) :: zfrac1

      if(zt1.le.zt(1)) then
         it1=1
         zfrac1=0
         return
      else if(zt1.ge.zt(int)) then
         it1=int-1
         zfrac1=1
         return
      endif

      it1 = 1 + (zt1-zt(1))*(int-1)/(zt(int)-zt(1))
      it1 = max(1,min((int-1),it1))
      do
         if(zt1.lt.zt(it1)) then
            it1=it1-1
         else if(zt1.gt.zt(it1+1)) then
            it1=it1+1
         else
            exit
         endif
      enddo
      zfrac1=(zt1-zt(it1))/(zt(it1+1)-zt(it1))

    end subroutine tdbsub_lookup

end module tdbsub_uts
