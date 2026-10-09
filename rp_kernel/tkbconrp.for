C******************** START FILE TKBCON.FOR ; GROUP TKBLOAT ************
C
C  TKBCONRP
C
 
      SUBROUTINE TKBCONRP(ISURF,ISURF0,ZITER)
 
C
C  ISURF0 - 1ST CONSTRUCTED (EXTRAPOLATED) SURFACE INDEX
C  ISURF - CURRENT INDEX OF SURFACE TO BE CONSTRUCTED
C
C  ZITER -- expansion method parameter
C
C
C  CONSTRUCT THE NEXT CONSTRUCTED MOMENTS SURFACE
 
 
C  TBT Jan 1994 -- CLeaned up Zangp
C
C  DMC DEC 1988 -- ADAPTED FOR USE AS AN RPLOT SUBROUTINE
C  BASED ON SUBROUTINE TKBCON FROM TRANSP
C
C  THIS ROUTINE CALLS TRANSP ROUTINES (W/O TRANSP COMMON)
C  ALSO TWO RPLOT VERSIONS OF TRANSP ROUTINES, TKBMMCRP.FOR AND
C  TKBMRYRP.FOR, HAVE BEEN CREATED, TO USE RPLOT INSTEAD OF
C  TRANSP COMMON BLOCKS
C
C  dmc Aug 1995 -- for asymmetric cases: added ZITER parameter and
C    copied over techniques learned in TRANSP.  Comments copied:
C
C  dmc 12 Aug 1994 -- added argument ZITER
C   try for an extrapolation scheme that doesn't bend the theta lines
C   so violently -- ZITER=0.0 gives this new scheme, ZITER=1.0 gives
C   the old scheme, intermediate values give intermediate solutions.
C   ZITER must be in range 0.0 to 1.0 inclusive
C
C   caller should loop over values of ZITER going from 0.0 to
C   1.0 until a successful (non-singular) extrapolation is achieved.
C
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
      Integer Itemp    ! tbt
 
      DATA TWOPI/6.283185/
      DATA XPI/3.141593/
C
      SAVE ZR0SV,ZY0SV
C
C=======================================================================
 
      ZANGP = 0.0   ! Initialize for compiler    ! tbt 1/94
 
C
C  NORMALIZED SPATIAL ELEMENT
C
      ZDXI=1.0/FLOAT(INX)
C
C  GET DATA ON NEAREST ALREADY CONSTRUCTED SURFACES
C
      INIT=ISURF-ISURF0
      CALL TKBRY2RP(ISURF,ZRM1,ZYM1,ZRM2,ZYM2,ZR0,ZY0,INIT,ICNTPS)
C
      if(init.eq.0) then
         zr0sv=zr0
         zy0sv=zy0
      endif
C
      CALL TKBRAD(ZRM1,ZYM1,ZR0,ZY0,ZDIST,ZRAD)
C
      ZXI=(ISURF-1)*ZDXI
      ZDRA=ZDXI*ZRAD/ZXI
C
C  CONSTRUCT R,Y PTS FOR NEW SURFACE
C
      ZACHK=0.04*XPI
C
C  gen for asymmetric geometry:
      ZRMID=ZR0
      ZYMID=ZY0
C
C  DMC
C  THIS IS AN EMPIRICAL DIFFERENTIAL EXPANSION OF THE FLUX COORD SYSTEM
C
      ZTHANL=TWOPI/ICNTPS
C
      DO 100 ITH=1,ICNTPS
C
C  I1,I2 GIVE NEIGHBOURING PT PTRS
C
        CALL TKBNDR(ITH,ZRM1,ZYM1,ZRM2,ZYM2,
     >       ZDERV,ZDERVN,ZANGN,ZANGS,I1,I2)
        ZTH=TWOPI*(ITH-1)/FLOAT(ICNTPS)
C
        ZTH0=FPOLAR((ZRM1(ITH)-ZRMID),(ZYM1(ITH)-ZYMID))
C
        ZR0L=SQRT((ZRM1(ITH)-ZRMID)**2+(ZYM1(ITH)-ZYMID)**2)
        ZRNL=0.5* (
     >        SQRT((ZRM1(I1)-ZRMID)**2+(ZYM1(I1)-ZYMID)**2) +
     >        SQRT((ZRM1(I2)-ZRMID)**2+(ZYM1(I2)-ZYMID)**2) )
C
        ZTHARC=FPOLAR((ZRM1(I1)-ZRMID),(ZYM1(I1)-ZYMID)) -
     >           FPOLAR((ZRM1(I2)-ZRMID),(ZYM1(I2)-ZYMID))
        ZTHARC=ABS(AMOD((ZTHARC+XPI+TWOPI),TWOPI)-XPI)
C
        ZDANG1=AMOD((ZTH-ZANGN+XPI+TWOPI),TWOPI)-XPI
        ZDANG2=AMOD((ZTH-ZTH0 +XPI+TWOPI),TWOPI)-XPI
C
        ZDANG=(ZDANG1+ZDANG2*AMIN1(1.0,(ZTHANL/ZTHARC)))
C
        ZDANG=(((1.0+ZDANG*ZDXI)**2)-1.0)*ABS(ZDANG)
C
        ZRXB=1.0+ZTHANL/ZTHARC*(ZRNL/ZR0L-1.0)
C
C  this number gives the ratio:  avg spacing / local spacing
C   where spacing is dist. btw. flux surfaces.  limit range of factor
C    needed to be fixed for heavily squeezed NSTX cases, dmc Apr 1995
C
        ZDRAT=max(0.4,min(2.5,(ZDRA/ZDERVN)))
        ZRATA=ZRXB*SQRT(ZRAD/ZDIST(ITH))*ZDRAT**1.5
        ZFAC=(1.0+ZDXI*(ZRATA-1.0))**2
        ZL1=ZDERV*ZFAC
C
        ZDANG=ZDANG*ZFAC *ZITER  ! ziter factor, dmc 12 Aug
C
        zl0=zderv*ZDRAT
        zl= ziter*zl1 + (1.0-ziter)*zl0  ! dmc 12 Aug
C
C  MOD DMC - IF PTS GETTING CLOSE TOGETHER, ACCELERATE ANGLE
C  NORMALIZATION (APR 1986)
C
        ZDL12=(ZRM1(I1)-ZRM1(ITH))**2+(ZYM1(I1)-ZYM1(ITH))**2
        ZDL22=(ZRM1(I2)-ZRM1(ITH))**2+(ZYM1(I2)-ZYM1(ITH))**2
C
        ZFAC2=0.5*(ZDL12**2+ZDL22**2)/(ZDL12*ZDL22)
C
        ZDANG=ZDANG*ZFAC2
C
        ZANG=AMIN1((ZANGS(2)-ZACHK),AMAX1((ZANGS(1)+ZACHK),
     >                 (ZANGN+ZDANG)))
C
        ZR(ITH)=ZRM1(ITH)+ZL*COS(ZANG)
        ZY(ITH)=ZYM1(ITH)+ZL*SIN(ZANG)
C
C  DMC FOR RPLOT ONLY -- CHECK IF ZANG IS MONOTONIC WITH THETA
C  INDEX;  IF SO SET FLAG
C
        IF(ITH.EQ.1) THEN
          ZANGPP=ZANG
          IFLAG=0
        ELSE IF(ITH.EQ.2) THEN
          IF(ZANG.GT.ZANGP) THEN
            IDIR=1
          ELSE IF(ZANG.LT.ZANGP) THEN
            IDIR=-1
          ELSE
            IDIR=0
          ENDIF
          ZANGPP=ZANGPP+TWOPI*IDIR
        ELSE
          IF(IDIR*(ZANG-ZANGP).LE.0.0) IFLAG=1
        ENDIF
        ZANGP=ZANG
 100  CONTINUE
C
C  COMPLETE RPLOT ZANG MONOTONICITY CHECK
      IF(IDIR*(ZANGPP-ZANGP).LE.0.0) IFLAG=1
C
C  MONOTONIC IF IFLAG=0
      IF(IFLAG.EQ.0) THEN
        IF(ICRVSF.EQ.0) ICRVSF=ISURF
      ENDIF
C
C  EVALUATE MOMENTS FROM R,Y SEQUENCE
C
CX	CALL GRAFX1(ZR,ZY,ICNTPS,'CM','CM','TKBCONRP',
CX     >            'ICNTPS PT EXPANSION CONTOUR','DEBUG PLOT')
C
      CALL TKBMMCRP(ISURF)
C
      RETURN
      END
C******************** END FILE TKBCONRP.FOR ; GROUP TKBLOAT ************
