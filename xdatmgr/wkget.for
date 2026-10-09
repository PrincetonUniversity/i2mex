C---------------------------------------------------------------
C  WKGET - COPY CONTIGUOUS WORKSPACE TO NONCONTIGUOUS COMMON DATA
C
      SUBROUTINE WKGET(ILOC,INC,ZWK,INUM)
C
      use datmgr_mod
C
      REAL ZWK(INUM)
C
      IL=ILOC-INC
C
      DO 10 I=1,INUM
        IL=IL+INC
        DATBUF(IL)=ZWK(I)
 10   CONTINUE
C
      RETURN
      END
