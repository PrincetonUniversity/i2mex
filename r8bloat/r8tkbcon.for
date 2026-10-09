C******************** START FILE TKBCON.FOR ; GROUP TKBLOAT ******************
C
C-- END FILE #QUADFUNS# ************************************
CFORMFEEDC-- START FILE #TKBCON# ************************************
C----------------------------------------------------------
C  TKBCON
C
C  dmc 27 Dec 1996 -- removing TRCOM so this routine can be used
C  outside TRANSP!  adding passed arguments.
C
C  B.BALET  FEB 1994 : ROUTINE UPDATED TO SUPPORT UP-DOWN ASYMMETRIC CASE
C
C  CONSTRUCT THE NEXT CONSTRUCTED MOMENTS SURFACE
C
      SUBROUTINE r8tkbcon(nlsym,xiblo,isurfx,zrmc,zymc,locmax,nmom,
     >      itflag,ziter)
C
C  nlsym -- .TRUE. if equilibrium is up-down symmetric
C  xiblo(...) -- extrapolated grid
C  isurfx -- index to last known grid pt. + 1
C  zrmc(...),zymc(...) -- local moments sets
C  locmax -- local moments max (array dimension)
C  nmom -- number of active moments
C  itflag -- =0 if this is the first extrapolated surface, gt.0 otherwise
C  ziter -- numerical control
C
C============
      use iso_c_binding, only: fp => c_double 
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER locmax,nmom,itflag,isurfx,ith,i1,i2
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 ziter,zr0sv,zy0sv,xpi,twopi,zdxi,zr0,zy0,zrad,zdra,zachk,
     > zrmid,zymid,zthanl,zderv,zdervn,zangn,zth,zth0,r8fpolar,zr0l,
     > zrnl,ztharc,zdang1,zdang2,zdang,zrxb,zdrat,zrata,zfac,zl1,zl0,zl,
     > zdl12,zdl22,zfac2,zang
C============
      logical nlsym
      REAL*8 xiblo(isurfx)
      REAL*8 zrmc(3,0:locmax,2)
      REAL*8 zymc(3,0:locmax,2)
C
C  dmc 12 Aug 1994 -- added argument ZITER
C   try for an extrapolation scheme that doesn't bend the theta lines
C   so violently -- ZITER=0.0 gives this new scheme, ZITER=1.0 gives
C   the old scheme, intermediate values give intermediate solutions.
C   ZITER must be in range 0.0 to 1.0 inclusive
C
C   caller (BLOAT) should loop over values of ZITER going from 0.0 to
C   1.0 until a successful (non-singular) extrapolation is achieved.
C
      include 'r8tkbparm.inc'
C
      REAL*8 ZR(NTK),ZY(NTK)
      REAL*8 ZRM1(NTK),ZYM1(NTK),ZRM2(NTK),ZYM2(NTK),ZDIST(NTK)
      REAL*8 ZANGS(2)
C
      SAVE ZR,ZY,ZRM1,ZYM1,ZRM2,ZYM2,ZDIST,ZANGS  !RTM 22 NOV 88
      SAVE ZR0SV,ZY0SV
C
C========================================================================
C
      data xpi/3.141593E0_R8/
      data twopi/6.283185E0_R8/
C
C  NORMALIZED SPATIAL ELEMENT
C
      ZDXI=XIBLO(ISURFX)-XIBLO(ISURFX-1)
C
C  GET DATA ON NEAREST ALREADY CONSTRUCTED SURFACES
C
      CALL r8tkbry2(zrmc,zymc,locmax,nmom,
     >           ZRM1,ZYM1,ZRM2,ZYM2,ZR0,ZY0,ITFLAG)
C
      if(itflag.eq.0) then
         zr0sv=zr0
         zy0sv=zy0
      endif
C
        CALL r8tkbrad(ZRM1,ZYM1,ZR0,ZY0,ZDIST,ZRAD)
C
      ZDRA=ZDXI*ZRAD/XIBLO(ISURFX-1)
C
C  CONSTRUCT R,Y PTS FOR NEW SURFACE
C
      ZACHK=0.04E0_R8*XPI
C
        IF(.not.NLSYM) THEN            !UP-DOWN ASYMMETRIC CASE
C
        ZRMID=ZR0     ! use line integral from TKBRY2
        ZYMID=ZY0
C
        ELSE                            !UP-DOWN SYMMETRIC CASE
C
          ZRMID=0.5E0_R8*(ZRM1(1)+ZRM1(NTK2))
          ZYMID=0.5E0_R8*(ZYM1(1)+ZYM1(NTK2))
C
        ENDIF
C
C  DMC
C  THIS IS AN EMPIRICAL DIFFERENTIAL EXPANSION OF THE FLUX COORD SYSTEM
C
      ZTHANL=TWOPI/NTK
C
      DO 100 ITH=1,NTK
C
C  I1,I2 GIVE NEIGHBOURING PT PTRS
C
        CALL r8tkbndr(ITH,ZRM1,ZYM1,ZRM2,ZYM2,
     >       ZDERV,ZDERVN,ZANGN,ZANGS,I1,I2)
        ZTH=TWOPI*(ITH-1)/REAL(NTK,kind=fp)
C
        ZTH0=r8fpolar((ZRM1(ITH)-ZRMID),(ZYM1(ITH)-ZYMID))
C
        ZR0L=SQRT((ZRM1(ITH)-ZRMID)**2+(ZYM1(ITH)-ZYMID)**2)
        ZRNL=0.5E0_R8* (
     >        SQRT((ZRM1(I1)-ZRMID)**2+(ZYM1(I1)-ZYMID)**2) +
     >        SQRT((ZRM1(I2)-ZRMID)**2+(ZYM1(I2)-ZYMID)**2) )
C
        ZTHARC=r8fpolar((ZRM1(I1)-ZRMID),(ZYM1(I1)-ZYMID)) -
     >           r8fpolar((ZRM1(I2)-ZRMID),(ZYM1(I2)-ZYMID))
        ZTHARC=ABS(mod((ZTHARC+XPI+TWOPI),TWOPI)-XPI)
C
        ZDANG1=mod((ZTH-ZANGN+XPI+TWOPI),TWOPI)-XPI
        ZDANG2=mod((ZTH-ZTH0 +XPI+TWOPI),TWOPI)-XPI
C
        ZDANG=(ZDANG1+ZDANG2*min(1.0E0_R8,(ZTHANL/ZTHARC)))
C
        ZRXB=1.0E0_R8+ZTHANL/ZTHARC*(ZRNL/ZR0L-1.0E0_R8)
C
C  this number gives the ratio:  avg spacing / local spacing
C   where spacing is dist. btw. flux surfaces.  limit range of factor
C    needed to be fixed for heavily squeezed NSTX cases, dmc Apr 1995
C
        ZDRAT=max(0.4E0_R8,min(2.5E0_R8,(ZDRA/ZDERVN)))
        ZRATA=ZRXB*SQRT(ZRAD/ZDIST(ITH))*ZDRAT**1.5E0_R8
        ZFAC=(1.0E0_R8+ZDXI*(ZRATA-1.0E0_R8))**2
        ZL1=ZDERV*ZFAC
C
        ZDANG=ZDANG*ZFAC *ZITER  ! ziter factor, dmc 12 Aug
C
        zl0=zderv*ZDRAT
        zl= ziter*zl1 + (1.0E0_R8-ziter)*zl0  ! dmc 12 Aug
 
C  MOD DMC - IF PTS GETTING CLOSE TOGETHER, ACCELERATE ANGLE
C  NORMALIZATION (APR 1986)
C
        ZDL12=(ZRM1(I1)-ZRM1(ITH))**2+(ZYM1(I1)-ZYM1(ITH))**2
        ZDL22=(ZRM1(I2)-ZRM1(ITH))**2+(ZYM1(I2)-ZYM1(ITH))**2
C
        ZFAC2=0.5E0_R8*(ZDL12**2+ZDL22**2)/(ZDL12*ZDL22)
C
        ZDANG=ZDANG*ZFAC2
C
        ZANG=min((ZANGS(2)-ZACHK),max((ZANGS(1)+ZACHK),
     >                 (ZANGN+ZDANG)))
C
        ZR(ITH)=ZRM1(ITH)+ZL*COS(ZANG)
          ZY(ITH)=ZYM1(ITH)+ZL*SIN(ZANG)
 100  CONTINUE
C
C  EVALUATE MOMENTS FROM R,Y SEQUENCE
C
cdbg      call r8_grafc3(zr,zy,NTK,zrm1,zym1,NTK,zrm2,zym2,NTK,'cm','cm',
cdbg     >    'TKBCON','NTK PT EXPANSION CONTOURS','DEBUG PLOT')
C
      CALL r8tkbmmc(nlsym,zrmc,zymc,locmax,nmom,ZR,ZY)
C
      RETURN
      END
C******************** END FILE TKBCON.FOR ; GROUP TKBLOAT ******************
