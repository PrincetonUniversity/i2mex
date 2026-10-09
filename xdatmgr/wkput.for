C---------------------------------------------------------------
C  WKPUT - COPY NONCONTIGUOUS DATA TO CONTIGUOUS WORKSPACE
C
      SUBROUTINE WKPUT(ILOC,INC,ZWK,INUM)
C
      use datmgr_mod
C
      REAL ZWK(INUM)
C
      IL=ILOC-INC
C
      DO 10 I=1,INUM
        IL=IL+INC
        ZWK(I)=DATBUF(IL)
 10   CONTINUE
C
      RETURN
      END
