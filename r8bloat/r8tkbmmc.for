! 30Nov1999 fgtok
! 30Nov1999 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
C******************** START FILE TKBMMC.FOR ; GROUP TKBLOAT ******************
C-- END FILE #TKBCON# ************************************
CFORMFEEDC-- START FILE #TKBMMC# ************************************
C-------------------------------------------------------------
C
      SUBROUTINE r8tkbmmc(NLSYM,zrmc,zymc,locmax,nmom,ZR,ZY)
C
C  local moments coeffs, nmom moments active:
C============
      use iso_c_binding, only: fp => c_double 
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nmom,locmax,itheta,j,ierr
C============
      REAL*8 zrmc(3,0:locmax,2),zymc(3,0:locmax,2)
C
C  updown symmetry flag:
      logical nlsym
C
C  dmc 27 Dec 1996:  remove TRCOM! and add calling arguments
C
C  B.BALET  FEB 1994 : UPDATE TO DEAL WITH UP-DOWN ASYMMETRIC CASE
C
C  FROM R,Y CONTOUR DATA DEFINE THE MOMENTS ARRAYS FOR 1 SURFACE
C
      include 'r8tkbparm.inc'
C
      REAL*8 ZR(NTK),ZY(NTK)  ! input contour
C
      REAL*8 ZS1D(NTK)
      REAL*8 ZC1D(NTK)
C
C---------------------------------
C
      ITHETA=NTK
C
C  IMSL FFT
C
      CALL r8fftsc(ZR, ITHETA, ZS1D, ZC1D, ierr)
C
C  NORMALIZE
C
      zrmc(3,0,1) = ZC1D(1)/(2.0E0_R8*ITHETA)
C
      DO 20310 J = 2, NMOM+1
C
         zrmc(3,J-1,1) = ZC1D(J)/REAL(ITHETA,kind=fp)
         if(.not.nlsym) zrmc(3,J-1,2) = ZS1D(J)/REAL(ITHETA,kind=fp)
C
20310 CONTINUE
C
C  IMSL FFT
C
      CALL r8fftsc(ZY, ITHETA, ZS1D, ZC1D, ierr)
C
C  NORMALIZE
C
      if(.not.nlsym) zymc(3,0,1)=ZC1D(1)/(2.0E0_R8*ITHETA)
C
      DO 20410 J = 2, NMOM+1
C
         if(.not.nlsym) ZYMC(3,J-1,1) = ZC1D(J)/REAL(ITHETA,kind=fp)
         ZYMC(3,J-1,2) = ZS1D(J)/REAL(ITHETA,kind=fp)
C
20410 CONTINUE
C
      RETURN
      END
C******************** END FILE TKBMMC.FOR ; GROUP TKBLOAT ******************
