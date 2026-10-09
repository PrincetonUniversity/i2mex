C******************** START FILE TKBRY2.FOR ; GROUP TKBLOAT ************
C-- END FILE #TKBRAD# ************************************
CFORMFEEDC-- START FILE #TKBRY2# ************************************
C------------------------------------------------------------------
      SUBROUTINE TKBRY2RP(ISURF,ZRM1,ZYM1,ZRM2,ZYM2,ZR0,ZY0,INIT,N)
C
C  GET DATA ON THE TWO SURFACES INSIDE SURFACE ISURF.
C  UNLESS INIT=0 ASSUME 2ND SURFACE GETS VALUES FROM OLD 1ST SURFACE
C  ARRAYS
C
      REAL ZRM1(N),ZYM1(N),ZRM2(N),ZYM2(N)
C
      ISM1=ISURF-1
      ISM2=ISM1-1
C
      IF(INIT.EQ.0) THEN
        CALL TKBRYSRP(ISM2,ZRM2,ZYM2,ZR0M2,ZY0M2,N)
      ELSE
        DO 10 ITH=1,N
          ZRM2(ITH)=ZRM1(ITH)
          ZYM2(ITH)=ZYM1(ITH)
 10     CONTINUE
      ENDIF
C
      CALL TKBRYSRP(ISM1,ZRM1,ZYM1,ZR0,ZY0,N)
C
      RETURN
      END
C******************** END FILE TKBRY2.FOR ; GROUP TKBLOAT **************
