C******************** START FILE DMDLOC.FOR ; GROUP PLDMGR ******************
C-------------------------------------------------------------
C  DMDLOC
C  LOCATE BLOCK OF DATA IN DATA AREA
C    BY NAME (5 CHARS)
C
      SUBROUTINE DMDLOC(ZNAME,IND,ISIZE,IPT)
C
      use datmgr_mod
      CHARACTER*(*) ZNAME
C
      character*50 zname2
C
C-------------------------------------
C
      zname2=zname
      call uupper(zname2)
C
 9    CONTINUE
      JL=1
 10   CONTINUE
      JL=LNEXT(JL)
      IF(DMGLBL(JL).EQ.'%FINI') GO TO 100
      IF(DMGLBL(JL).NE.ZNAME2 .and. DMGLBL(JL).NE.ZNAME) GO TO 10
C  NAME FOUND
      CALL DATREF(JL,IPT,ISIZE)
      IND=JL
      RETURN
C  NAME NOT FOUND
 100  CONTINUE
      if(zname2.eq.'X') then
         zname2='%XXXC'
         go to 9                        ! try to recover
      else if(zname2.eq.'XB') then
         zname2='%XXXB'
         go to 9                        ! try to recover
      endif
      ISIZE=0
      IPT=0
      IND=0
      RETURN
      END
C******************** END FILE DMDLOC.FOR ; GROUP PLDMGR ******************
