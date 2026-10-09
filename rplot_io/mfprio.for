C******************** START FILE MFPRIO.FOR ; GROUP PLOTR3 ******************
C---------------
C  JACK PRIORITY OF INDICATED DATA AND SAVE PTR INDEX
C
      SUBROUTINE MFPRIO(ILIS,INLIS,ICT,IND)
C
      use datmgr_mod
      INTEGER ILIS(INLIS)
C
      IF(IND.EQ.0) RETURN
      CALL PLPRIO(6,IND)
      ICT=ICT+1
      ILIS(ICT)=IND
      RETURN
      END
C******************** END FILE MFPRIO.FOR ; GROUP PLOTR3 ******************
