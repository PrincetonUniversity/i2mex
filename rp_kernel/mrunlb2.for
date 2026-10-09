C-----------------------------------------------------------------------
C  MRUNLB2 -- MAKE RUN ID LABEL FOR RPLOT
C
      SUBROUTINE MRUNLB2(ZDIR,ILDIR,ZRUNID,ILRUNID,ZLABEL)
C
      CHARACTER*(*) ZDIR   ! INPUT DIRECTORY NAME
      INTEGER ILDIR        ! INPUT NAME LENGTH
      CHARACTER*(*) ZRUNID  ! INPUT RUN ID
      INTEGER ILRUNID      ! INPUT RUN ID LENGTH
C
      CHARACTER*(*) ZLABEL ! OUTPUT LABEL
C
      CHARACTER*20 ZSUB
C
C
      character*10 zdigs
      data zdigs /'0123456789'/
      integer inum_nodig
C-----------------------------------------------------------------------
C
      ZSUB=' '
C
      ILD=ILDIR
      if(ZDIR(ILD:ILD).eq.'/') ILD=ILD-1
      if (ild-ilrunid .gt. 0) then
       if (zdir(ild-ilrunid+1:ild) .eq. ZRUNID(1:ILRUNID)) then
         ild=ild-ILRUNID-1
         if (zdir(ild-ilrunid+1:ild) .eq. ZRUNID(1:ILRUNID)) then
            ild=ild-ILRUNID-1
         endif
       endif
      endif
      inum_nodig=0  ! count number of not-digit-characters
      DO IC=ILD,1,-1
        IF((ZDIR(IC:IC).EQ.'/').and.(inum_nodig.gt.0)) GO TO 10
        if(index(zdigs,zdir(ic:ic)).le.0) inum_nodig=inum_nodig+1
      ENDDO
      ZSUB='RUN '
      GO TO 90
C
 10   CONTINUE
      ZSUB=ZDIR(IC+1:ILD)
      do ic=1,len(zsub)
         if(zsub(ic:ic).eq.'/') zsub(ic:ic)='.'
      enddo
C
C  NEW RUN LABEL
 90   CONTINUE
      IBLNK=ILRUNID+1
      IHD=len_trim(zsub)
      WRITE(ZLABEL,5020) ZSUB(1:IHD),ZRUNID(1:ILRUNID)
 5020 FORMAT(A,1X,A)
C
      RETURN
      END
