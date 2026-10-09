C******************** START FILE PLMMCV.FOR ; GROUP PLOTMM ******************
C-----------------------------------------------------------------------
C  PLMMCV
C
C  CALCULATE MOMENTS CONTOURS
 
 
      SUBROUTINE PLMMCV(ZRMOM0,ZRMOM,ZYMOM,Jmom,Jxl,
     >			ZRC,ZYC,JNC,JNX)
 
 
C  Mod:
C  04/13/94 tbt Added common PLFMPA - ZCosTabl & ZSintabl

      use cplotr_mod
      use plfmpa_mod

      integer :: jmom,jxl,jnc,jnx
 
      REAL ZRMOM(Jmom),ZYMOM(Jmom)
      REAL ZRC(JNC,JNX),ZYC(JNC,JNX)
 
C--------------------------------------------------------------
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
      DO ITH=1,JNC
 
         ZRC8(ITH)=ZRMOM0
         ZYC8(ITH)=0.0
 
         DO IM=1,Jmom
 
            ZRC8(ITH)=ZRC8(ITH)+ZRMOM(IM)*ZCosTabl(Im, Ith)
            ZYC8(ITH)=ZYC8(ITH)+ZYMOM(IM)*ZSinTabl(Im, Ith)

         ENDDO
      ENDDO
C     
      go to 1000
C
C-------------------------
C  precomputed tables not available
C
 100  continue
      zcos(0)=1
      zsin(0)=0
      DO ITH=1,JNC
 
         ZRC8(ITH)=ZRMOM0
         ZYC8(ITH)=0.0
 
         ZTH=ztwopi*(ITH-1)/FLOAT(JNC-1)
         call sincos(zth,jmom,zsin(1),zcos(1))
	
         DO IM=1,Jmom
 
            ZRC8(ITH)=ZRC8(ITH)+ZRMOM(IM)*ZCOS(IM)
            ZYC8(ITH)=ZYC8(ITH)+ZYMOM(IM)*ZSIN(IM)

         ENDDO
      ENDDO

 1000 continue

      zrc(1:jnc,jxl)=zrc8
      zyc(1:jnc,jxl)=zyc8
 
      RETURN
C
      END
C******************** END FILE PLMMCV.FOR ; GROUP PLOTMM ******************
