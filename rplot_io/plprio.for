C---------------------------------------------------------------------
C PLPRIO  RESTORE STANDARD PRIORITY OF RPLOT DATA ITEMS
C
C  EXEMPT USER DEFINED FCNS (PRIO 7) & program defined workspaces (prio 10).
C
      SUBROUTINE PLPRIO(IPRIO,IND)
C
      use datmgr_mod
      use cplotr_mod
C
      IF(IND.GT.0) THEN
         IF((MPRIO(IND).NE.7).and.(MPRIO(IND).NE.10)) then
            MPRIO(IND)=IPRIO
         ENDIF
      ENDIF
C
      RETURN
      END
