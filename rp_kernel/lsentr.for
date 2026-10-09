C******************** START FILE LSENTR.FOR ; GROUP LSIGEN ******************
C--------------------------------------------
C  (2)  LSENTR
C     FINALIZE ENTRY IF GRAPHS UNDER ITS LABEL WERE DRAWN
C
      SUBROUTINE LSENTR
C
      use cplotr_mod
C
C  WERE GRAPHS DRAWN ???
C
      IDREW=NPAGEG-LSPAGI(NLSENT)+1
      IF(IDREW.GT.0) GO TO 50
C  NO GRAPHS DRAWN; DELETE ENTRY
      NLSENT=NLSENT-1
      RETURN
C  GRAPHS DRAWN; RECORD NUMBER
 50   CONTINUE
      LSPAGL(NLSENT)=IDREW
      RETURN
      	END
C******************** END FILE LSENTR.FOR ; GROUP LSIGEN ******************
