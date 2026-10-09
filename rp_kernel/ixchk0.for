C******************** START FILE IXCHK0.FOR ; GROUP IXCALC ******************
C----------------------------
C  IXCHK0
C  CHECK AGAINS F=0 SINGULARITY
C
      SUBROUTINE IXCHK0(F,SMALL)
C
      IF(ABS(F).LT.SMALL) THEN
        IF(F.LT.0.0) F=-SMALL
        IF(F.GE.0.0) F=+SMALL
      ENDIF
      RETURN
      END
C******************** END FILE IXCHK0.FOR ; GROUP IXCALC ******************
