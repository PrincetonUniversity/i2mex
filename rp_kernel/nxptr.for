C******************** START FILE NXPTR.FOR ; GROUP PLDMGR ******************
C---------------------------------------------------------------
C  NXPTR
C
C  RETURN POINTER TO TIME VARYING X AXIS'S VARIATION AT A PARTICULAR
C  TIME
C
C  JX-- X AXIS ID CODE
C  IT-- AT IT'TH TIME
C
C  DATA MUST BE IN CORE
C
      INTEGER FUNCTION NXPTR(JX,IT)
C
      use datmgr_mod
      use cplotr_mod
C
      CHARACTER*21 ZNAM
C
      INX=NZONEX(JX)
C
      NXPTR=0
C
      IF(JX.EQ.1) ZNAM='%XAZC'
      IF(JX.EQ.2) ZNAM='%XAZB'
      IF(JX.GT.2 .or. .not.NLTRANSP) ZNAM=ABR(NFX(JX))
C  LOCATE DATA
      CALL DMDLOC(ZNAM,IND1,ISIZ1,IPTR)
      IF(IND1.EQ.0) RETURN
      NXPTR=IPTR+(IT-1)*INX
      RETURN
      END
C******************** END FILE NXPTR.FOR ; GROUP PLDMGR ******************
