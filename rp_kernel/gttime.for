C******************** START FILE GTTIME.FOR ; GROUP LSIGEN ******************
C-----
C  GET TIME; AVOID NAME CONFLICT "TIME"
C
      SUBROUTINE GTTIME(ZZTIME)
      implicit none

      CHARACTER*8 ZZTIME
      integer :: itime(8)
C
      call clocaltim(itime)
      write(zztime,'(i2.2,a1,i2.2,a1,i2.2)') itime(5),':',
     &     itime(6),':',itime(7)

      !CALL TIME(ZZTIME)
      RETURN
      END
C******************** END FILE GTTIME.FOR ; GROUP LSIGEN ******************
