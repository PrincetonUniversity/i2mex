      subroutine rptime_p(zbuf,ibufsize,iret)
C
C  fetch the timebase for PROFILE functions of time
C
      use datmgr_mod
      use cplotr_mod
C
      real zbuf(ibufsize)               ! buffer into which to copy time vector
C
      integer iret                      ! number of time points copied
C
C  if ibufsize is too small, iret=0 is returned!
C
      iret=0
      if(ibufsize.lt.ntr) then
         call rpbufsiz('rptime_p',ibufsize,ntr)
         return
      endif
C
      do it=1,ntr
         zbuf(it)=time3(it)
      enddo
      iret=ntr
C
      return
      end
 
      subroutine r8_rptime_p(zbufr8,ibufsize,iret)
C
C  fetch timebase -- to R8 array
c
      implicit NONE
c
      integer :: ibufsize
      real*8 :: zbufr8(ibufsize)
      integer :: iret
c
c---------------
c
      real, dimension(:), allocatable :: zbuf
c
c---------------
c
      allocate(zbuf(ibufsize)); zbuf = 0.0
c
      call rptime_p(zbuf,ibufsize,iret)
c
      zbufr8 = zbuf
      deallocate(zbuf)
c
      return
      end
