      subroutine vmsieee(lun_warn,zvals,ilim,ufieee,ieee)
c
c     just copy the numbers
c
c  overflow/underflow check provided.
c
      integer lun_warn                  ! LUN for warning messages
      integer ilim
      real zvals(ilim)                  ! data to convert or recover
      real ufieee(ilim)                 ! converted data (ieee format)
      integer ieee                      ! =1:  convert to ieee
c                                       ! =0:  convert from ieee back
      iwarnuf=0
      iwarnof=0
c
      do i=1,ilim
         if(ieee.eq.1) then
C
C  convert to ieee
C  check for over/underflow
C
            if(zvals(i).ne.0.0) then
               if(abs(zvals(i)).lt.1.1754945e-38) then
                  zvals(i)=0.0
                  if((iwarnuf.eq.0).and.(lun_warn.gt.0)) then
                     write(lun_warn,
     >                  '('' %ubwfcmp:  underflow warning.'')')
                     iwarnuf=1
                  endif
               else if(zvals(i).gt.8.4070587e+37) then
                  zvals(i)=8.4070587e+37
                  if((iwarnof.eq.0).and.(lun_warn.gt.0)) then
                     write(lun_warn,
     >                  '('' %ubwfcmp: +overflow warning.'')')
                     iwarnof=1
                  endif
               else if(zvals(i).lt.-8.4070587e+37) then
                  zvals(i)=-8.4070587e+37
                  if((iwarnof.eq.0).and.(lun_warn.gt.0)) then
                     write(lun_warn,
     >                  '('' %ubwfcmp: -overflow warning.'')')
                     iwarnof=1
                  endif
               endif
            endif
            ufieee(i)=zvals(i)
         else
c
c  convert back *from* ieee; no error check
c
            zvals(i)=ufieee(i)
c
         endif
      enddo
c
      return
      end
 
