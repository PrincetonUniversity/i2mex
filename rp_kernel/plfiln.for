C-----------------------------------------------------------------
C  PLFILN - RPLOT SUBROUTINE, CONSTRUCT INPUT DATA FILENAME
C
C  PASSED STRING ZZEND IS THE LAST 2 CHARACTERS OF THE MAIN FILENAME
C   FOLLOWED BY "." FOLLOWED BY THE FILENAME EXTENSION; E.G. "MF.PLN"
C
C  RETURNED SRING ZFILN IS THE FILENAME; ZFILN ASSUMED LONG ENUF
C
      SUBROUTINE PLFILN(ZZEND,ZFILN)
C
C
C  DECLARATIONS
C   COMMON BLOCKS FOR PLOTTING---
      use cplotr_mod
C       --------------------
C
      CHARACTER*(*) ZZEND
      CHARACTER*(*) ZFILN
C
      character*40 zdisk
      character*140 zdir
      character*10 zrunid
C
C------------------------------------
C
      luntrm=6
C
      IL=0
C
C  lengths of filename components
C
      ILEN=LEN(ZFILN)
      ILEND=LEN(ZZEND)
C
      if(lrun_x.eq.0) then
         zdisk=fdisk
         ilfdisk=lfdisk
         zdir=fdir
         ilfdir=lfdir
         zrunid=runid
         ilrunid=lrunid
      else
         zdisk=fdisk_x(lrun_x)
         if(zdisk.eq.' ') then
            ilfdisk=0
         else
            ilfdisk=ilnurd(fdisk_x(lrun_x))
         endif
         zdir=fdir_x(lrun_x)
         if(zdir.eq.' ') then
            ilfdir=0
         else
            ilfdir=ilnurd(fdir_x(lrun_x))
         endif
         zrunid=runid_x(lrun_x)
         ilrunid=ilnurd(zrunid)
      endif
C
      if (zdisk(1:4) .eq. 'MDS+') then
cdbg         write(luntrm,9902)
cdbg 9902    format(' ?? PLFILN: skipping disk/dir for MDS+')
         continue
      else
         IF(ILFDISK.GT.0) THEN
            ILP=IL+1
            IL=IL+ILFDISK
            IF(IL.GT.ILEN) GO TO 900
            ZFILN(ILP:IL)=ZDISK(1:ILFDISK)
         ENDIF
C
         IF(ILFDIR.GT.0) THEN
            ILP=IL+1
            IL=IL+ILFDIR
            IF(IL.GT.ILEN) GO TO 900
            ZFILN(ILP:IL)=ZDIR(1:ILFDIR)
         ENDIF
      endif
C
      IF(ILRUNID.GT.0) THEN
        ILP=IL+1
        IL=IL+ILRUNID
        IF(IL.GT.ILEN) GO TO 900
        ZFILN(ILP:IL)=ZRUNID(1:ILRUNID)
      ENDIF
C
      IF(ILEND.GT.0) THEN
        ILP=IL+1
        IL=IL+ILEND
        IF(IL.GE.ILEN) GO TO 900
        ZFILN(ILP:IL)=ZZEND
      ENDIF
C
      IL=IL+1
      DO 10 ILB=IL,ILEN
        ZFILN(ILB:ILB)=' '
 10   CONTINUE
C
      GO TO 1000
C
C  ERROR
C
 900  CONTINUE
      write(luntrm,9901) zzend
 9901 format(' PLFILN:  LENGTH ERROR CONSTRUCTING ',A,' FILENAME.')
C
 1000 CONTINUE
C
      RETURN
      END
