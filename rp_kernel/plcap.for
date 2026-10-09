C******************** START FILE PLCAP.FOR ; GROUP PLOTR1 ******************
C---------------------------------------------
C  CAPITALIZE RPLOT ENTITY ABBREVIATIONS
C
      SUBROUTINE PLCAP(ABREV,N)
C
      CHARACTER*(*) ABREV(*)
C
      integer iwarn
C
      IF(N.EQ.0) RETURN
      DO 10 I=1,N
        CALL uupper(ABREV(I))
        call idchek0(abrev(i),iwarn,1)
cxx        write(6,*) abrev(i),iwarn
 10   CONTINUE
      RETURN
      END
C******************** END FILE PLCAP.FOR ; GROUP PLOTR1 ******************
