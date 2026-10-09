! 30Nov1999 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
C******************** START FILE TKBNDR.FOR ; GROUP TKBLOAT ******************
C-- END FILE #TKBMRY# ************************************
CFORMFEEDC-- START FILE #TKBNDR# ************************************
C--------------------------------------------------------------
C  RETURN NUMERIC NORMAL DERIVATIVE
C
      SUBROUTINE r8tkbndr(ITH,ZRM1,ZYM1,ZRM2,ZYM2,
     >                    ZDERV,ZDERVN,ZANGN,ZANGS,I1,I2)
C
C  INPUT - ITH = PT NUMBER
C    ARRAYS ZRM1,ZYM1 - DESCRIBE NEAREST CLOSED SURFACE
C    ARRAYS ZRM2,ZYM2 - DESCRIBE NEXT SURFACE IN
C
C  OUTPUT - ZDERVN - "NORMAL DERIVATIVE";  NORM DOT (DL VECTOR ALONG
C   LINE OF CONST. THETA GOING FROM #2 SURFACE TO #1 SURFACE, I.E.
C     ( (ZRM1(ITH)-ZRM2(ITH)), (ZYM1(ITH)-ZYM2(ITH)) )
C
C  ZANGN - ANGLE IN (R,Y) SPACE OF DL VECTOR
C  ZANGS(2) - CONSTRAINT ANGLES DUE TO SURFACE ELEMENT ON #1
C    ABOUT POINT ITH
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      include 'r8tkbparm.inc'
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER i1,i2,ith,ithm,i,idir
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 zderv,zdervn,zangn,zdr1,zdy1,r8fpolar,zdr2,zdy2,zangsp,
     > zangsn,z,zd2
C============
      REAL*8 ZRM1(NTK),ZYM1(NTK),ZRM2(NTK),ZYM2(NTK),ZANGS(2)
C  STATEMENT FUNCTION
      ITHM(I)=MOD((I+NTK-1),NTK)+1
C
      ZDR1=ZRM1(ITH)-ZRM2(ITH)
      ZDY1=ZYM1(ITH)-ZYM2(ITH)
C
      ZANGN=r8fpolar(ZDR1,ZDY1)
C
      I1=ITHM(ITH+1)
      I2=ITHM(ITH-1)
C
      ZDR2=ZRM1(I2)-ZRM1(I1)
      ZDY2=ZYM1(I2)-ZYM1(I1)
C
      ZANGSP=r8fpolar(ZDR2,ZDY2)
      IF(ZANGSP.LT.ZANGN) THEN
        IDIR=1
      ELSE
        IDIR=-1
      ENDIF
C
 10   CONTINUE
      ZANGSN=ZANGSP+IDIR*3.141593E0_R8
      IF((ZANGN-ZANGSP)*(ZANGN-ZANGSN).GT.0.0E0_R8) THEN
        ZANGSP=ZANGSN
        GO TO 10
      ELSE                    !Avoid a CIVIC bug
        GO TO 20		!Avoid a CIVIC bug
      ENDIF
20    CONTINUE		!Avoid a CIVIC bug
C
C  CONSTRAINT ANGLES FOR NEXT STEP OF THETA LINE - AVOID SINGULARITY
      ZANGS(1)=min(ZANGSN,ZANGSP)
      ZANGS(2)=max(ZANGSN,ZANGSP)
C
      Z=ZDR1*ZDY2-ZDY1*ZDR2
      ZD2=SQRT(ZDR2**2+ZDY2**2)
C
      ZDERVN=ABS(Z)/ZD2
C
      ZDERV=SQRT(ZDR1**2+ZDY1**2)
C
      RETURN
      END
C******************** END FILE TKBNDR.FOR ; GROUP TKBLOAT ******************
