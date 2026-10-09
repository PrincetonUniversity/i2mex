C******************** START FILE PLAMMCV.FOR ; GROUP PLOTMM ******************
C-----------------------------------------------------------------------
C  PLAMMCV
C
C  CALCULATE Asymmetric MOMENTS CONTOURS
 
 
      SUBROUTINE PLAMMCV(ZRCMOM,ZRSMOM, ZYCMOM,ZYSMOM,JMOM,JXL,
     >			   ZRC,ZYC,Jnc,Jnx)
 
C       04/12/94  tbt Made from PLMMCV.

      use cplotr_mod
      use plfmpa_mod
 
      integer :: jmom,jxl,jnc,jnx

      REAL ZRCMOM(0:JMOM),ZYCMOM(0:JMOM)
      REAL ZRSMOM(0:JMOM),ZYSMOM(0:JMOM)
      REAL ZRC(Jnc,Jnx),ZYC(Jnc,Jnx)
C
C-------------------------------
      real zcos(0:jmom),zsin(0:jmom)
C
      real*8 :: zrc8(jnc),zyc8(jnc)
      real*8 :: ztwopi
      integer :: ith,im
C-------------------------------
C
      data ztwopi/6.2831853071795862D+00/
C
      if(JNC.ne.NCOSTABL) go to 100
C
C  use precomputed sin/cos tables
C
 
      DO ITH=1,Jnc
 
         ZRC8(ITH)=0.0
         ZYC8(ITH)=0.0
 
         DO IM=0,JMOM
 
            ZRC8(ITH)=ZRC8(ITH)+ZRCMOM(IM) * ZCosTabl(Im, Ith)
     1           +ZRSMOM(IM) * ZSinTabl(Im, Ith)
            ZYC8(ITH)=ZYC8(ITH)+ZYCMOM(IM) * ZCosTabl(Im, Ith)
     1           +ZYSMOM(IM) * ZSinTabl(Im, Ith)
 
         ENDDO
      ENDDO
C
      go to 1000
C-------------------------------
C
 100  continue
C
C  precomputed table not available; use sins and cosines
C
      zcos(0)=1
      zsin(0)=0
      DO ITH=1,Jnc
 
         ZRC8(ITH)=0.0
         ZYC8(ITH)=0.0
 
         ZTH=Ztwopi*(ITH-1)/FLOAT(Jnc-1)
         call sincos(zth,jmom,zsin(1),zcos(1))
 
         DO IM=0,JMOM
 
            ZRC8(ITH)=ZRC8(ITH)+ZRCMOM(IM)*ZCOS(IM)
     1           +ZRSMOM(IM)*ZSIN(IM)
            ZYC8(ITH)=ZYC8(ITH)+ZYCMOM(IM)*ZCOS(IM)
     1           +ZYSMOM(IM)*ZSIN(IM)
 
         ENDDO
      ENDDO

 1000 continue

      zrc(1:jnc,jxl)=zrc8
      zyc(1:jnc,jxl)=zyc8

      RETURN
      END
C******************** END FILE PLAMMCV.FOR ; GROUP PLOTMM ******************
