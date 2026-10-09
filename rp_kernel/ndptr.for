C******************** START FILE NDPTR.FOR ; GROUP PLDMGR ******************
C---------------------------------------------------------------
C  NDPTR
C
C  RETURN PTR TO PROFILE FCN AT A PARTICULAR TIME (IT)
C
C  JF-- FCN ID CODE
C  INX-- NUMBER OF PTS IN PROFILE
C  IT-- AT IT'TH TIME
C
C  DATA MUST BE IN CORE
C
      INTEGER FUNCTION NDPTR(JF,INX,IT)
C
      use datmgr_mod
      use cplotr_mod
C
      NDPTR=0
C
C  LOCATE DATA
      CALL DMDLOC(ABR(JF),IND1,ISIZ1,IPTR)
      IF(IND1.EQ.0) RETURN
      NDPTR=IPTR+(IT-1)*INX
      RETURN
      END
C******************** END FILE NDPTR.FOR ; GROUP PLDMGR ******************
