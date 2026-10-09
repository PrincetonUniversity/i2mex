C-----------------------------------------------------------------
C  ISCAL0 - SCALE DEFAULT LABELING SUBUTILITY USED BY SCALIT (RPLOT)
C
      SUBROUTINE ISCAL0(IXID,IAX,ISCAL)
C
      use cplotr_mod
C
      IAX=NAXISC(IXID)
C
      IF((NSCALC(1,IXID).EQ.1).AND.(NSCALC(2,IXID).EQ.1)) THEN
        ISCAL=1
      ELSE
        ISCAL=2
      ENDIF
C
      RETURN
      END
