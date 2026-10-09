! 30Nov1999 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
C******************** START FILE TKBRAD.FOR ; GROUP TKBLOAT ******************
C-- END FILE #TKBNDR# ************************************
CFORMFEEDC-- START FILE #TKBRAD# ************************************
C---------------------------------------------------------
C  SUBROUTINE TKBRAD
C
      SUBROUTINE r8tkbrad(ZR,ZY,ZR0,ZY0,ZDISTA,ZRAD)
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      include 'r8tkbparm.inc'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ith,ithm1
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 zy0,zrad,zr0,zdist,zr1,zr2,zy1,zy2,zrad0,zradp,zl,zdl,zradn
C============
      REAL*8 ZR(NTK),ZY(NTK),ZDISTA(NTK)
C
C  RETURN AVG DISTANCE TO SURFACE PTS (ZR,ZY) FROM CTR PT (ZR0,ZY0)
C
      ZDIST(ZR1,ZR2,ZY1,ZY2)=SQRT((ZR2-ZR1)**2+(ZY2-ZY1)**2)
C
      ZRAD0=ZDIST(ZR0,ZR(1),ZY0,ZY(1))
      ZDISTA(1)=ZRAD0
C
      ZRADP=ZRAD0
C
      ZL=0.0E0_R8
      ZRAD=0.0E0_R8
C
      DO 10 ITH=2,NTK
        ITHM1=ITH-1
        ZDL=ZDIST(ZR(ITH),ZR(ITHM1),ZY(ITH),ZY(ITHM1))
        ZL=ZL+ZDL
        ZRADN=ZDIST(ZR0,ZR(ITH),ZY0,ZY(ITH))
        ZDISTA(ITH)=ZRADN
        ZRAD=ZRAD+ZDL*0.5E0_R8*(ZRADN+ZRADP)
C
        ZRADP=ZRADN
 10   CONTINUE
C
C COMPLETE LOOP INTEGRAL
C
      ZDL=ZDIST(ZR(NTK),ZR(1),ZY(NTK),ZY(1))
      ZL=ZL+ZDL
      ZRAD=ZRAD+ZDL*0.5E0_R8*(ZRADP+ZRAD0)
C
      ZRAD=ZRAD/ZL
C
      RETURN
      END
C******************** END FILE TKBRAD.FOR ; GROUP TKBLOAT ******************
