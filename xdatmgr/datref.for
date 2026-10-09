C******************** START FILE DATREF.FOR ; GROUP DATMGR ******************
C-------------------------------------------------------------
C  DATREF
C
C  REFERENCE DATA ENTRY J-- UPDATE ITS ACCESS AND RETURN POINTERS
C  TO DATA TO MAIN PROGRAM
C
      SUBROUTINE DATREF(J,ILOC,ISIZ)
C
      use datmgr_mod
C
C  UPDATE ACCESS COUNT
C
      call dmg_macc_incr
      LACC(J)=MACC
C
      IF(LOCD(J).EQ.0) THEN
        WRITE(LUNDMO,1001)
 1001 FORMAT(' ?DATMGR -- REFERENCE TO UNDEFINED DATA')
        call bad_exit
      ENDIF
      ILOC=LOCD(J)
      ISIZ=NWDS(J)
C
      RETURN
      END
C******************** END FILE DATREF.FOR ; GROUP DATMGR ******************
