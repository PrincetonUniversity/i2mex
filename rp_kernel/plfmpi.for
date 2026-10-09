C======================================================================
C  PLFMPI  ROOT FUNCTION FOR NAG INVERSE MAP / ROOT FINDER
C   DECLARED EXTERNAL IN PLFRYG AND PASSED THRU NAG C05NBE SUBROUTINE
C
      SUBROUTINE PLFMPI(INDIM,ZX,ZF,IFLAG)
C
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
      REAL ZX(INDIM),ZF(INDIM)
C
      REAL ZJACO(2,2)
C
C---------------------------------------
C
      ZXI=ZX(1)
      ZTH=ZX(2)
C
      CALL PLJACO(ZXI,ZTH,0,ZROUT,ZYOUT,ZJACO,ISTAT)
C
      ZF(1)=ZROUT-ZRTARG
      ZF(2)=ZYOUT-ZYTARG
C
      RETURN
      END
