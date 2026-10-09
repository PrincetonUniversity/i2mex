C******************** START FILE PLMGEO.FOR ; GROUP PLOTMM ******************
C-----------------------------------------------------------------------
C  PLMGEO
C
C  CALCULATE GEOMETRIC QUANTITIES ASSOCIATED WITH MOMENTS SURFACE
 
 
      SUBROUTINE PLMGEO(ZRMOM0,ZRMOM,ZYMOM,JMOM,
     >     ZRADIUS,ZA,ZELONG,ZDELTA,ZINDNT)
 
      use cplotr_mod
      use plfmpa_mod
 
C  INPUT ZRMOM,ZYMOM,ZRMOM0 - THE MOMENTS DESCRIBING A SURFACE
C
      REAL ZRMOM(JMOM),ZYMOM(JMOM)
C  OUTPUT
C    ZRADIUS - *MIDPLANE* RADIUS
C    ZA  - *MIDPLANE* HALFWIDTH
C
C    ZELONG - SURFACE VERTICAL ELONGATION
C    ZDELTA - SURFACE TRIANGULARITY
C    ZINDNT - SURFACE INDENTATION
C
C  NUMBER OF SURFACE EVALUATIONS...
CCC	PARAMETER (IEVAL=60)

      real zra(3), zya(3)
      Real ZRR, ZYY
C ---------------------------------------------------------------
 
C  MIDPLANE INTERCEPTS...
 
      ZR1=ZRMOM0
      ZR2=ZRMOM0
      DO 10 IM=1,JMOM
        ZR2=ZR2+ZRMOM(IM)
        ISIGN=1-2*(IM-2*(IM/2))  ! +1 IF IM EVEN, -1 IF IM ODD
        ZR1=ZR1+ISIGN*ZRMOM(IM)
 10   CONTINUE
      ZRGEO=0.5*(ZR1+ZR2)
C
      ZRADIUS=ZRGEO
      ZA=0.5*(ZR2-ZR1)
C
      ZRMAX=0.0
      ZRMIN=ZR2
      ZZMAX=0.0
      ZRZMAX=0.0
 
      Ieval = NaxMmp
      ict=1
      DO 20 JTH=-1,IEVAL
         if(jth.lt.1) then
            ith=ieval+jth
         else
            ith=jth
         endif
 
         ZRL=ZRMOM0
         ZYL=0.0
 
         DO 30 IM=1,JMOM
 
            ZRL=ZRL+ZRMOM(IM)* ZCosTabl(Im, Ith)
            ZYL=ZYL+ZYMOM(IM)* ZSinTabl(Im, Ith)
 
 30      CONTINUE
 
         ict=ict+1
         if(ict.le.3) then
            zra(ict)=zrl
            zya(ict)=zyl
         else
            zra(1)=zra(2)
            zra(2)=zra(3)
            zra(3)=zrl
            zya(1)=zya(2)
            zya(2)=zya(3)
            zya(3)=zyl
         endif
         if(jth.lt.1) go to 20
 
         ZRMIN=AMIN1(ZRMIN,ZRL)
         ZRMAX=AMAX1(ZRMAX,ZRL)
         IF(ZYL.GT.ZZMAX) THEN
            ZZMAX =ZYL
            ZRZMAX=ZRL
         ENDIF
C
C  improve estimates by parabolic extrapolation
C
         zzR0=zra(2)
         zzAR=(zra(3)-zra(1))*0.5
         zzBR=(zra(3)+zra(1))*0.5 - zzR0
 
         if((zra(2).ge.max(zra(1),zra(3))).and.
     >      (zra(2).gt.1.000001*min(zra(1),zra(3)))) then
            zx=-zzAR/(2.0*zzBR)
            zrx=zzR0+zx*(zzAr+ zx*zzBR)
            ZRMAX=max(ZRMAX,zrx)
         endif
 
         if((zra(2).le.min(zra(1),zra(3))).and.
     >      (zra(2).lt.0.999999*max(zra(1),zra(3)))) then
            zx=-zzAR/(2.0*zzBR)
            zrx=zzR0+zx*(zzAr+ zx*zzBR)
            ZRMIN=min(ZRMIN,zrx)
         endif
 
         zzY0=zya(2)
         zzAY=(zya(3)-zya(1))*0.5
         zzBY=(zya(3)+zya(1))*0.5 - zzY0
 
         if((zya(2).ge.max(zya(1),zya(3))).and.
     >      (zya(2).gt.1.000001*min(zya(1),zya(3)))) then
            zx=-zzAY/(2.0*zzBY)
            zyx=zzY0+zx*(zzAy+ zx*zzBY)
            if(zyx.gt.zzmax) then
               zzmax=zyx
               zrzmax=zzR0+zx*(zzAr+ zx*zzBR)
            endif
         endif
 
 20   CONTINUE
 
      ZZA=0.5*(ZRMAX-ZRMIN)
      ZELONG=ZZMAX/ZZA
      ZINDNT=(ZR1-ZRMIN)/(ZRMAX-ZRMIN)
      ZDELTA=(ZRZMAX-ZRGEO)/(0.5*(ZRMAX-ZRMIN))
 
      RETURN
      END
C******************** END FILE PLMGEO.FOR ; GROUP PLOTMM ******************
 
 
 
 
 
