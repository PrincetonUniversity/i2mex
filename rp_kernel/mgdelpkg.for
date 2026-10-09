      subroutine mgdelpkg(idel)
C
C  delete multigraph package #idel
C  ...taken from subroutine muldef dmc 11 Aug 1999
C
      use cplotr_mod
C
      character*10 zabb
C------------------------------------
C
      zabb=abb(idel)
      write(lunzer(0),1001) abb(idel)
 1001 format(' %mgdelpkg:  deleted package:  ',a)
      call aordr_del(abb,iordrb,NAXMGP,nbal,zabb) ! nbal decremented
C
C  aordr_del fixed abb and iordrb arrays;
C  fixup the associated data...
C
      IF(IDEL.EQ.(NBAL+1)) GO TO 1000
      DO 120 IP=IDEL,NBAL
        IPP1=IP+1
        LABELB(IP)=LABELB(IPP1)
        UNITSB(IP)=UNITSB(IPP1)
        DO 117 J=1,naxmgf
          IFUNB(J,IP)=IFUNB(J,IPP1)
 117    CONTINUE
        IINTB(IP)=IINTB(IPP1)
        INFB(IP)=INFB(IPP1)
 120  CONTINUE
C
 1000 continue
      return
      end
 
