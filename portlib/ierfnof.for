C-----------------------------------------------------------------------
C  IERFNOF -- return the integer error code (IOSTAT) for file not found
C    on this machine
C
C    it is "29" only on DEC machines!
C
      integer function IERFNOF(lun)
C
      data iansr/0/
C
      if(iansr.eq.0) then
          open(unit=lun,file='Xyz.GarbageFile',status='old',
     >        iostat= iansr )
          close(unit=lun,err=10)
 10       continue
      endif
C
      ierfnof = iansr
C
      return
      end
