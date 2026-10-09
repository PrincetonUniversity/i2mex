C---------------------------------------------------------------------
C  PLMMRLS -- REDUCE PRIORITY (RELEASE) OF MOMENTS DATA READ IN
C
      SUBROUTINE PLMMRLS(IDIM)
C
C  ARGUMENT IDIM = 1 -- PARTIAL SET, =2 -- FULL SET.  SEE PLMMRD
C   SUBROUTINE
C
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
C---------------------------------------------------------------------
C
      IF(IMMGEO.EQ.0) THEN
C
C  LOWER PRIORITY ON SYMMETRIC MOMENTS DATA
C
        CALL PLPRIO(5,INDR0)
C
        DO 30 IM=1,IMOM
          CALL PLPRIO(5,INDMM(1,IM))
          IF(IDIM.GT.1) CALL PLPRIO(5,INDMM(2,IM))
 30     CONTINUE
C
      ELSE
C
C  LOWER PRIORITY ON ASYMMETRIC MOMENTS DATA
C
        CALL PLPRIO(5,INDYMP)
        CALL PLPRIO(5,INDRMP)
        IF(IDIM.GT.1) THEN
          CALL PLPRIO(5,INDAMM(1,0))
          CALL PLPRIO(5,INDAMM(3,0))
          DO IM=1,IMOM
            DO IS=1,4
              CALL PLPRIO(5,INDAMM(IS,IM))
            ENDDO
          ENDDO
        ENDIF
C
      ENDIF
C
      RETURN
      END
