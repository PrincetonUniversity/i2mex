! 30Nov1999 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
C******************** START FILE TKBRY2.FOR ; GROUP TKBLOAT ******************
C-- END FILE #TKBRAD# ************************************
CFORMFEEDC-- START FILE #TKBRY2# ************************************
C------------------------------------------------------------------
      SUBROUTINE r8tkbry2(zrmc,zymc,locmax,nmom,
     >                    ZRM1,ZYM1,ZRM2,ZYM2,ZR0,ZY0,INIT)
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nmom,init,locmax,ism1,ism2,ith
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 zr0,zy0,zr0m2,zy0m2
C============
      REAL*8 zrmc(3,0:locmax,2)
      REAL*8 zymc(3,0:locmax,2)
C
C  mod DMC 29 Dec 1996:
C  passed args added, to avoid later need to refer to COMMON...
C
C  GET DATA ON THE TWO SURFACES INSIDE SURFACE target surface.
C  UNLESS INIT=0 ASSUME 2ND SURFACE GETS VALUES FROM OLD 1ST SURFACE
C  ARRAYS
C
      include 'r8tkbparm.inc'
      REAL*8 ZRM1(NTK),ZYM1(NTK),ZRM2(NTK),ZYM2(NTK)
C
      ISM1=2
      ISM2=1
C
      IF(INIT.EQ.0) THEN
        CALL r8tkbrys(ISM2,zrmc,zymc,locmax,nmom,
     >                ZRM2,ZYM2,ZR0M2,ZY0M2)
      ELSE
        DO 10 ITH=1,NTK
          ZRM2(ITH)=ZRM1(ITH)
          ZYM2(ITH)=ZYM1(ITH)
 10     CONTINUE
      ENDIF
C
      CALL r8tkbrys(ISM1,zrmc,zymc,locmax,nmom,
     >              ZRM1,ZYM1,ZR0,ZY0)
C
      RETURN
      END
C******************** END FILE TKBRY2.FOR ; GROUP TKBLOAT ******************
