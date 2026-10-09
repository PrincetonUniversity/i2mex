!******************** START FILE TKBNDR.FOR ; GROUP TKBLOAT ******************
!-- end FILE #TKBMRY# ************************************
!FORMFEEDC-- START FILE #TKBNDR# ************************************
!--------------------------------------------------------------
!  return NUMERIC NORMAL DERIVATIVE
!
      subroutine TKBNDR(ITH,ZRM1,ZYM1,ZRM2,ZYM2, &
                          ZDERV,ZDERVN,ZANGN,ZANGS,I1,I2)
use iso_c_binding, only: fp => c_double
!
!  INPUT - ITH = PT NUMBER
!    ARRAYS ZRM1,ZYM1 - DESCRIBE NEAREST CLOSED SURFACE
!    ARRAYS ZRM2,ZYM2 - DESCRIBE NEXT SURFACE IN
!
!  OUTPUT - ZDERVN - "NORMAL DERIVATIVE";  NORM doT (DL VECTOR ALONG
!   LINE OF CONST. THETA GOING FROM #2 SURFACE TO #1 SURFACE, I.E.
!     ( (ZRM1(ITH)-ZRM2(ITH)), (ZYM1(ITH)-ZYM2(ITH)) )
!
!  ZANGN - ANGLE IN (R,Y) SPACE OF DL VECTOR
!  ZANGS(2) - CONSTRAINT ANGLES DUE TO SURFACE ELEMENT ON #1
!    ABOUT POINT ITH
!
      include 'tkbparm.inc'
!
      REAL ZRM1(NTK),ZYM1(NTK),ZRM2(NTK),ZYM2(NTK),ZANGS(2)
!  STATEMENT FUNCTION
      ITHM(I)=MOD((I+NTK-1),NTK)+1
!
      ZDR1=ZRM1(ITH)-ZRM2(ITH)
      ZDY1=ZYM1(ITH)-ZYM2(ITH)
!
      ZANGN=FPOLAR(ZDR1,ZDY1)
!
      I1=ITHM(ITH+1)
      I2=ITHM(ITH-1)
!
      ZDR2=ZRM1(I2)-ZRM1(I1)
      ZDY2=ZYM1(I2)-ZYM1(I1)
!
      ZANGSP=FPOLAR(ZDR2,ZDY2)
      if(ZANGSP.LT.ZANGN) THEN
        IDIR=1
      else
        IDIR=-1
      end if
!
 10   continue
      ZANGSN=ZANGSP+IDIR*3.141593
      if((ZANGN-ZANGSP)*(ZANGN-ZANGSN).GT.0.0) THEN
        ZANGSP=ZANGSN
        goto 10
      else                    !Avoid a CIVIC bug
        goto 20		!Avoid a CIVIC bug
      end if
20    continue		!Avoid a CIVIC bug
!
!  CONSTRAINT ANGLES FOR NEXT STEP OF THETA LINE - AVOID SINGULARITY
      ZANGS(1)=AMIN1(ZANGSN,ZANGSP)
      ZANGS(2)=AMAX1(ZANGSN,ZANGSP)
!
      Z=ZDR1*ZDY2-ZDY1*ZDR2
      ZD2=SQRT(ZDR2**2+ZDY2**2)
!
      ZDERVN=ABS(Z)/ZD2
!
      ZDERV=SQRT(ZDR1**2+ZDY1**2)
!
      return
      end
!******************** end FILE TKBNDR.FOR ; GROUP TKBLOAT ******************
