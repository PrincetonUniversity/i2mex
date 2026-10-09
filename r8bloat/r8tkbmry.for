! 30Nov1999 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
C******************** START FILE TKBMRY.FOR ; GROUP TKBLOAT ******************
 
 
C-- END FILE #TKBMMC# ************************************
CFORMFEEDC-- START FILE #TKBMRY# ************************************
C-----------------------------------------------------------------
C  DMC - UTILITY ROUTINE FOR MOMENT SURFACES EXTRAPOLATION.
C
C  RETURN THE R,Y PT. FOR GIVEN SURFACE, PASSED THETA VALUE
C
      SUBROUTINE r8tkbmry(nsurf,ISURF,ZTHETA,zrmc,zymc,locmax,nmom,
     >                    ZROUT,ZYOUT)
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER locmax,nmom,isurf,locmax2,imn
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 ztheta,zrout,zyout,zzrmcc,zzrmcs,zzymcc,zzymcs
C============
      integer nsurf
      REAL*8 zrmc(nsurf,0:locmax,2)
      REAL*8 zymc(nsurf,0:locmax,2)
C
      REAL*8 snthtk(locmax),csthtk(locmax)
C
C--------------------------------------------------
C
C
      CALL r8sincos ( ZTHETA, NMOM, SNTHTK, CSTHTK )
C
C  UP-DOWN ASYMMETRIC CASE (B.BALET JAN 94) :
C
      ZROUT = ZRMC(ISURF,0,1)
      ZYOUT = ZYMC(ISURF,0,1)
      DO 10610 IMN = 1, NMOM
         ZZRMCC = ZRMC(ISURF,IMN,1)
         ZZRMCS = ZRMC(ISURF,IMN,2)
         ZZYMCC = ZYMC(ISURF,IMN,1)
         ZZYMCS = ZYMC(ISURF,IMN,2)
         ZROUT = ZROUT + ZZRMCC * CSTHTK(IMN) + ZZRMCS * SNTHTK(IMN)
         ZYOUT = ZYOUT + ZZYMCC * CSTHTK(IMN) + ZZYMCS * SNTHTK(IMN)
10610 CONTINUE
C
      RETURN
      END
C******************** END FILE TKBMRY.FOR ; GROUP TKBLOAT ******************
