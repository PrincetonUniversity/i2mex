      subroutine ubi3ncod(ivalue,ilbuf,istart)
C
C  encode a non-negative integer (btw 0 and 16M) in 3 bytes,
C   hi order bits stored first.
C
      implicit NONE
C
      integer ivalue
 
      integer*1 ilbuf(*)
 
      integer byte1,byte2,byte3
      integer istart
C
      integer i2to8,i2to16,jvalue
C
      data i2to8/256/
      data i2to16/65536/
C-----------------------------------
C
      jvalue=ivalue
      if(jvalue.lt.0) then
         write(6,'('' ?ubi3ncod:  negative value passed.'')')
         go to 999
      endif
C
      byte1=jvalue/i2to16
      if(byte1.gt.255) then
         write(6,'('' ?ubi3ncod:  passed value too large.'')')
         go to 999
      endif
C
      jvalue=jvalue-i2to16*byte1
C
      byte2=jvalue/i2to8
C
      jvalue=jvalue-i2to8*byte2
C
      byte3=jvalue
C
C-----------
C
      if(byte1.le.127) then
         ilbuf(istart)=byte1
      else
         ilbuf(istart)=byte1-256
      endif
C
      if(byte2.le.127) then
         ilbuf(istart+1)=byte2
      else
         ilbuf(istart+1)=byte2-256
      endif
C
      if(byte3.le.127) then
         ilbuf(istart+2)=byte3
      else
         ilbuf(istart+2)=byte3-256
      endif
C
C debugging:
cx      call ubi3dcod(ilbuf,istart,kvalue)
cx      write(6,'('' ubi3ncod:  '',2(1x,i15))') ivalue,kvalue
C
      go to 1000
C
 999  continue
C
 1000 continue
      return
      end
