!******************** START FILE TKBRAD.FOR ; GROUP TKBLOAT ******************
!-- end FILE #TKBNDR# ************************************
!FORMFEEDC-- START FILE #TKBRAD# ************************************
!---------------------------------------------------------
!  subroutine TKBRAD
!
      subroutine TKBRAD(ZR,ZY,ZR0,ZY0,ZDISTA,ZRAD)
use iso_c_binding, only: fp => c_double
!
      include 'tkbparm.inc'
      REAL ZR(NTK),ZY(NTK),ZDISTA(NTK)
!
!  return AVG DISTANCE TO SURFACE PTS (ZR,ZY) FROM CTR PT (ZR0,ZY0)
!
      ZDIST(ZR1,ZR2,ZY1,ZY2)=SQRT((ZR2-ZR1)**2+(ZY2-ZY1)**2)
!
      ZRA_fp=ZDIST(ZR0,ZR(1),ZY0,ZY(1))
      ZDISTA(1)=ZRA_fp
!
      ZRADP=ZRA_fp
!
      ZL=0.0
      ZRAD=0.0
!
      do 10 ITH=2,NTK
        ITHM1=ITH-1
        ZDL=ZDIST(ZR(ITH),ZR(ITHM1),ZY(ITH),ZY(ITHM1))
        ZL=ZL+ZDL
        ZRADN=ZDIST(ZR0,ZR(ITH),ZY0,ZY(ITH))
        ZDISTA(ITH)=ZRADN
        ZRAD=ZRAD+ZDL*0.5*(ZRADN+ZRADP)
!
        ZRADP=ZRADN
 10   continue
!
! COMPLETE LOOP INTEGRAL
!
      ZDL=ZDIST(ZR(NTK),ZR(1),ZY(NTK),ZY(1))
      ZL=ZL+ZDL
      ZRAD=ZRAD+ZDL*0.5*(ZRADP+ZRA_fp)
!
      ZRAD=ZRAD/ZL
!
      return
      end
!******************** end FILE TKBRAD.FOR ; GROUP TKBLOAT ******************
