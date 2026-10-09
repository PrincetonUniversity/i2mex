      subroutine ttintrp(time3,zdat3,ntr,time,zdata,ntt)
C
C  interpolate from profile timebase to scalar timebase
C
      implicit NONE

      integer :: ntr,ntt

      real time3(ntr)                   ! profile timebase (input)
      real zdat3(ntr)                   ! data on profile timebase (input)
C
      real time(ntt)                    ! scalar timebase (input)
      real zdata(ntt)                   ! data interpolated to scalar timebase
C
C------------------------------
      integer :: itr,it
      real :: zfac,zt
C---------------------------------------------------------
C
      itr=1
      zfac=0.0
C
      do it=1,ntt
         zt=time(it)
         if(zt.le.time3(1)) then
            zdata(it)=zdat3(1)
         else if(zt.ge.time3(ntr)) then
            zdata(it)=zdat3(ntr)
         else
C
 10         continue
            if(zt.gt.time3(itr+1)) then
               itr=itr+1
               go to 10
            endif
C
C  zt.le.time3(itr+1) (& .lt. time3(ntr))
C
            if(time3(itr).eq.time3(itr+1)) then
               zfac=0.5                 ! degenerate case
            else
               zfac=(zt-time3(itr))/(time3(itr+1)-time3(itr))
            endif
            zfac=max(0.0,zfac)
C
            zdata(it)=zdat3(itr)+zfac*(zdat3(itr+1)-zdat3(itr))
C
         endif
      enddo
C
      return
      end
