C-----------------------------------------------------------------------
C  PLTIMI -- FIND PROFILE TIME INTERPOLATION FACTORS
C   **** profile database ****
C
      SUBROUTINE PLTIMI(ZTIME,IT1,IT2,Z1,Z2,IEXP)
C
      use datmgr_mod
C
      use cplotr_mod
C
C  INPUT ZTIME = TIME TO INTERPOLATE TO
C
C  OUTPUT INTERPOLATION FACTORS S.T.
C    INTERPOLATION = Z1*DATA(IT1)+Z2*DATA(IT2)
C
C  PROVIDE FOR FLAT EXTRAPOLATION WITH FLAG IEXP PASSED BACK
C   IEXP=0 INDICATES NO EXTRAPOLATION
C
C  CHECK ENDPTS
      IEXP=0
      ZFIX=ZTIME
C
      IF(ZFIX.GE.TIME3(NTR)) THEN
         ZFIX=TIME3(NTR)
         IT1=NTR
         IT2=NTR
         Z1=1.0
         Z2=0.0
         IF(ZTIME.GT.TIME3(NTR)) IEXP=2
      ELSE IF(ZFIX.LE.TIME3(1)) THEN
         ZFIX=TIME3(1)
         IT1=1
         IT2=1
         Z1=1.0
         Z2=0.0
         IF(ZTIME.LT.TIME3(1)) IEXP=1
      ELSE
C  INTERPOLATE BETWEEN 2 INTERIOR TIME PTS
         INTM1=NTR-1
         DO 50 IT=1,INTM1
            ITP1=IT+1
            IF((TIME3(IT).LE.ZFIX).AND.(TIME3(ITP1).GE.ZFIX))
     >         GO TO 60
 50      CONTINUE
 60      CONTINUE
         ZDTI=1./(TIME3(ITP1)-TIME3(IT))
         IT1=IT
         IT2=ITP1
         Z1=(TIME3(IT2)-ZFIX)*ZDTI
         Z2=(ZFIX-TIME3(IT1))*ZDTI
      ENDIF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  PLTIMI_SC -- FIND scalar TIME INTERPOLATION FACTORS
C   **** scalar database ****
C
      SUBROUTINE PLTIMI_SC(ZTIME,IT1,IT2,Z1,Z2,IEXP)
C
      use datmgr_mod
C
      use cplotr_mod
C
C  INPUT ZTIME = TIME TO INTERPOLATE TO
C
C  OUTPUT INTERPOLATION FACTORS S.T.
C    INTERPOLATION = Z1*DATA(IT1)+Z2*DATA(IT2)
C
C  PROVIDE FOR FLAT EXTRAPOLATION WITH FLAG IEXP PASSED BACK
C   IEXP=0 INDICATES NO EXTRAPOLATION
C
C  CHECK ENDPTS
      IEXP=0
      ZFIX=ZTIME
C
      IF(ZFIX.GE.TIME(NTT)) THEN
         ZFIX=TIME(NTT)
         IT1=NTT
         IT2=NTT
         Z1=1.0
         Z2=0.0
         IF(ZTIME.GT.TIME(NTT)) IEXP=2
      ELSE IF(ZFIX.LE.TIME(1)) THEN
         ZFIX=TIME(1)
         IT1=1
         IT2=1
         Z1=1.0
         Z2=0.0
         IF(ZTIME.LT.TIME(1)) IEXP=1
      ELSE
C  INTERPOLATE BETWEEN 2 INTERIOR TIME PTS
         INTM1=NTT-1
         DO 50 IT=1,INTM1
            ITP1=IT+1
            IF((TIME(IT).LE.ZFIX).AND.(TIME(ITP1).GE.ZFIX))
     >         GO TO 60
 50      CONTINUE
 60      CONTINUE
         ZDTI=1./(TIME(ITP1)-TIME(IT))
         IT1=IT
         IT2=ITP1
         Z1=(TIME(IT2)-ZFIX)*ZDTI
         Z2=(ZFIX-TIME(IT1))*ZDTI
      ENDIF
C
      RETURN
      END
