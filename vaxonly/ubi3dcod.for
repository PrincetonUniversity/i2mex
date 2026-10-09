      subroutine ubi3dcod(ilbuf,istart,ivalue)
C
C  decode a 3 byte sequence as an integer, hi order bits first
C
      implicit NONE
C
      integer ivalue
 
      integer*1 ilbuf(*)
 
      integer istart,ibyte(3)
C
      integer ict,is
C
      ict=0
      do is=istart,istart+2
         ict=ict+1
         if(ilbuf(is).lt.0) then
            ibyte(ict)=256+ilbuf(is)
         else
            ibyte(ict)=ilbuf(is)
         endif
      enddo
C
      ivalue=256*(256*ibyte(1)+ibyte(2))+ibyte(3)
C
      return
      end
