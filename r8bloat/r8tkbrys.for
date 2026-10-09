! 30Nov1999 fgtok
! 30Nov1999 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
C******************** START FILE TKBRYS.FOR ; GROUP TKBLOAT ******************
C-- END FILE #TKBRY2# ************************************
CFORMFEEDC-- START FILE #TKBRYS# ************************************
C-----------------------------------------------------------------
C  DMC - UTILITY ROUTINE FOR MOMENT SURFACES EXTRAPOLATION.
C
C  mod dmc 29 Dec 1996:  add passed args, to avoid necessity of later
C  reference to COMMON...
C
C  GIVEN A SURFACE AT THE SPECIFIED MOMENTS INDEX LOCATION THIS
C  ROUTINE RETURNS THE (R,Y) AS A FCN OF THETA, NTK PTS. FROM 0 TO
C  NTK-1/NTK * TWOPI; AND GIVES A ROUGH ESTIMATE OF THE CENTROID
C  (R0,Y0) DEFINED FROM AN APPROXIMATE PT TO PT LINE INTEGRATION.
C
      SUBROUTINE r8tkbrys(ISURF,zrmc,zymc,locmax,nmom,
     >                    ZR,ZY,ZR0,ZY0)
C
C============
      use iso_c_binding, only: fp => c_double 
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER locmax,nmom,isurf,ith,ithm1
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 zr0,zy0,twopi,zl,zth,zdl
C============
      REAL*8 zrmc(3,0:locmax,2)
      REAL*8 zymc(3,0:locmax,2)
C
      include 'r8tkbparm.inc'
      REAL*8 ZR(NTK),ZY(NTK)
C
      DATA TWOPI/6.283185E0_R8/
C
C------------------------------------
C
      CALL r8tkbmry(3,ISURF,0.0E0_R8,zrmc,zymc,locmax,nmom,ZR(1),ZY(1))
C
      ZR0=0.0E0_R8
      ZY0=0.0E0_R8
      ZL=0.0E0_R8
C
      DO 100 ITH=2,NTK
        ITHM1=ITH-1
C
        ZTH=TWOPI*REAL(ITH-1,kind=fp)/REAL(NTK,kind=fp)
        CALL r8tkbmry(3,ISURF,ZTH,zrmc,zymc,locmax,nmom,
     >                ZR(ITH),ZY(ITH))
C
        ZDL=SQRT((ZR(ITH)-ZR(ITHM1))**2+(ZY(ITH)-ZY(ITHM1))**2)
C
        ZL=ZL+ZDL
        ZR0=ZR0+ZDL*0.5E0_R8*(ZR(ITH)+ZR(ITHM1))
        ZY0=ZY0+ZDL*0.5E0_R8*(ZY(ITH)+ZY(ITHM1))
C
 100  CONTINUE
C
C  COMPLETE THE LOOP
C
        ZDL=SQRT((ZR(1)-ZR(NTK))**2+(ZY(1)-ZY(NTK))**2)
C
        ZL=ZL+ZDL
        ZR0=ZR0+ZDL*0.5E0_R8*(ZR(1)+ZR(NTK))
        ZY0=ZY0+ZDL*0.5E0_R8*(ZY(1)+ZY(NTK))
C
        ZR0=ZR0/ZL
        ZY0=ZY0/ZL
C
      RETURN
      END
C******************** END FILE TKBRYS.FOR ; GROUP TKBLOAT ******************
