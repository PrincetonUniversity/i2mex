C******************** START FILE PLAMGEO.FOR ; GROUP PLOTMM ******************
C-----------------------------------------------------------------------
C  PLAMGEO
C
C  CALCULATE GEOMETRIC QUANTITIES ASSOCIATED WITH MOMENTS SURFACE
 
 
 
      SUBROUTINE PLAMGEO(ZRCMOM, ZRSMOM, ZYCMOM, ZYSMOM, JMOM,
     >     ZR0p,ZY0p,ZA,ZELONG,ZDELTA,ZINDNT)
 
      use cplotr_mod
      use plfmpa_mod

C     04/12/94 TBT Created from PLMGEO and GetMPa.for
 
C  INPUT ZRCMOM,ZRSMOM,ZYCMOM, ZYSMOM- THE MOMENTS DESCRIBING A SURFACE
C         R cos, R sin, Y cos,  Y sin
C
      REAL ZRCMOM(0:JMOM), ZYCMOM(0:JMOM)
      Real ZRSmom(0:Jmom), ZYSmom(0:Jmom)
 
C  more input (dmc jun 94)
C    ZR0p,ZY0p - *MIDPLANE* RADIUS and elevation
C    ZA  - *MIDPLANE* HALFWIDTH
C
C  OUTPUT
C    ZELONG - SURFACE VERTICAL ELONGATION
C    ZDELTA - SURFACE TRIANGULARITY
C    ZINDNT - SURFACE INDENTATION
C
 
 
C  NUMBER OF SURFACE EVALUATIONS...

      real zra(3),zya(3)
      Real ZRR, ZYY
 
C	========================================================
C
C  MIDPLANE INTERCEPTS...
 
      ZRMAX=-1.E10
      ZRMIN= 1.E10
      ZZMAX=-1.E10
      ZZMIN= 1.E10
      ZRZMAX=0.0
      ZRZMIN=0.0
 
      Ieval = NaxMmp
      ict=1
      DO 20 JTH=-1,IEVAL
         if(jth.lt.1) then
            ith=ieval+jth
         else
            ith=jth
         endif
 
         ZRL=0.0
         ZYL=0.0
 
         DO 30 IM=0,JMOM
 
            ZRL=ZRL+ZRCMOM(IM) * ZcosTabl(Im,Ith)
     1           +ZRSmom(IM) * ZsinTabl(Im,Ith)
            ZYL=ZYL+ZYCmom(IM) * ZcosTabl(Im,Ith)
     1           +ZYSMOM(IM) * ZsinTabl(Im,Ith)
 
 30      CONTINUE
 
         ict=ict+1
         if(ict.le.3) then
            zra(ict)=zrl
            zya(ict)=zyl-ZY0P
         else
            zra(1)=zra(2)
            zra(2)=zra(3)
            zra(3)=zrl
            zya(1)=zya(2)
            zya(2)=zya(3)
            zya(3)=zyl-ZY0P
         endif
         if(jth.lt.1) go to 20
 
         ZRMIN=AMIN1(ZRMIN,ZRL)
         ZRMAX=AMAX1(ZRMAX,ZRL)
         IF(ZYL-ZY0p.GT.ZZMAX) THEN
            ZZMAX =ZYL-ZY0p
            ZRZMAX=ZRL
         ENDIF
         IF(ZYL-ZY0p.LT.ZZMIN) THEN
            ZZMIN =ZYL-ZY0p
            ZRZMIN=ZRL
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
 
         if((zya(2).le.min(zya(1),zya(3))).and.
     >      (zya(2).lt.0.999999*max(zya(1),zya(3)))) then
            zx=-zzAY/(2.0*zzBY)
            zyx=zzY0+zx*(zzAy+ zx*zzBY)
            if(zyx.lt.zzmin) then
               zzmin=zyx
               zrzmin=zzR0+zx*(zzAr+ zx*zzBR)
            endif
         endif
 
 
 20   CONTINUE
 
      ZELONG=0.5*(ZZMAX-ZZMIN)/ZA
      ZINDNT=((ZR0p-ZA)-ZRMIN)/(ZRMAX-ZRMIN)
      ZDELTA=(0.5*(ZRZMAX+ZRZMIN)-ZR0p)/(0.5*(ZRMAX-ZRMIN))
C
      RETURN
      END
C******************** END FILE PLAMGEO.FOR ; GROUP PLOTMM ******************
