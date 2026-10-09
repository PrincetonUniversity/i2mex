      subroutine rptime_s(zbuf,ibufsize,iret)
C
      use datmgr_mod
      use cplotr_mod
C
C  fetch the timebase for SCALAR functions of time
C
      real zbuf(ibufsize)               ! buffer into which to copy time vector
C
      integer iret                      ! number of time points copied
C
C  if ibufsize is too small, iret=0 is returned!
C
      iret=0
      if(ibufsize.lt.ntt) then
         call rpbufsiz('rptime_s',ibufsize,ntt)
         return
      endif
C
      do it=1,ntt
         zbuf(it)=time(it)
      enddo
      iret=ntt
C
      return
      end

 
      subroutine r8_rptime_s(zbufr8,ibufsize,iret)
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
      call rptime_s(zbuf,ibufsize,iret)
c
      zbufr8 = zbuf
      deallocate(zbuf)
c
      return
      end
