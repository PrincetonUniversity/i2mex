C******************** START FILE SMOOF1.FOR ; GROUP SMOOF2 ******************
C=============================================================
C  SMOOF1
C
C  SMOOTH A 1D FUNCTION GIVEN A CONST. WEIGHTED AVG WIDTH
C  AND A CONSTANT 'RELATIVE' EPSILON
C
C
      SUBROUTINE SMOOF1(DELTA,EPSR,EPSTYPE,X,N)
C
      use datmgr_mod
C
C argument input:
C
C  DELTA-- WGHTED AVG WIDTH FOR FILTR PROGRAM
C
C  EPSR -- ERROR BAR (cf EPSTYPE), e.g. if EPSTYPE = '%' then
C          RELATIVE ERROR BAR IN PERCENT; E.G. EPSR=20.0 ==>
C          SMOOTHED OUTPUT SHOULD NOT DIFFER FROM UNSMOOTHED
C          INPUT BY MORE THAN 20 %; see EPSTYPE...
C  EPSTYPE -- '%' if EPSR is percentage
C             'R' if EPSR is fractional (specifies max relative change)
C             otherwise EPSR is absolute (specifies max absolute change)
C             if EPSR is positive, or else abs(EPSR) is fractional if EPSR
C             is negative.  If EPSR=0.0, then there is no limit on change
C             caused by smoothing.
C
C  N    -- # OF PTS IN FUNCTION BEING SMOOTHED
C  X(N)    FUNCTION INDEPENDANT VARIABLE
C
C
C  INPUT THRU COMMON:
C   SMWORK( --- , 1)  INPUT FUNCTION TO BE SMOOTHED
C  OUTPUT THRU COMMON:
C   SMWORK( --- , 2)  OUTPUT SMOOTHED VERSION OF FUNCTION
C
C  WORKSPACE IN COMMON:
C   SMWORK( ---, 3)  ARRAY OF DELTA VALUES FOR FILTR6
C   SMWORK  ( ---, 4)  ARRAY OF EPSILON VALUES (ABSOLUTE) FOR FILTR6
C
      character*1 EPSTYPE
C
      DIMENSION X(N)
C
C  A VERY SMALL NUMBER:
C
C  if DELTA=0, just copy data
      if(DELTA.eq.0.0) then
         do i=1,n
            smwork(i,2)=smwork(i,1)
         enddo
         RETURN
      endif
C
C  check EPSR & EPSTYPE arguments
C
      if(EPSR.eq.0.0) then
         ifrac=0
         ZEPSA=ZLARGE
      else if(EPSTYPE.eq.'%') then
         ifrac=1
         ZEFRAC=EPSR/100.
      else if(EPSTYPE.eq.'R') then
         ifrac=1
         ZEFRAC=EPSR
      else
         if(EPSR.lt.0.0) then
            ifrac=1
            ZEFRAC=EPSR
         else
            ifrac=0
            ZEPSA=EPSR
         endif
      endif
C
      ZEMIN=ZLARGE
      ZEMAX=0.0
C
      ZDELTA=ABS(DELTA)
C  DBL INVERT?
      INVRT=0
      IF(ZDELTA.NE.DELTA) INVRT=1
      ISIGN=0
      ISIGNC=0
      DO 820 I=1,N
         SMWORK(I,3)=ZDELTA
         if(ifrac.eq.0) then
            ZEPS=ZEPSA
         else
            ZEPS=ABS(ZEFRAC*SMWORK(I,1))
         endif
         IF(ZEPS.GT.0.0) ZEMIN=AMIN1(ZEPS,ZEMIN)
         ZEMAX=AMAX1(ZEPS,ZEMAX)
         SMWORK(I,4)=ZEPS
C  INVERSION 1?
         IF(INVRT.EQ.1) THEN
            IF(ABS(SMWORK(I,1)).GT.ZSMALL) THEN
               SMWORK(I,1)=1/SMWORK(I,1)
            ELSE
               IF(SMWORK(I,1).LT.0.0) THEN
                  SMWORK(I,1)=-ZLARGE
               ELSE
                  SMWORK(I,1)=ZLARGE
               ENDIF
            ENDIF
         ENDIF
         IF((ISIGN.EQ.0).AND.(SMWORK(I,1).GT.0.0)) ISIGN=1
         IF((ISIGN.EQ.0).AND.(SMWORK(I,1).LT.0.0)) ISIGN=-1
         IF((ISIGN.EQ.1).AND.(SMWORK(I,1).LT.0.0)) ISIGNC=1
         IF((ISIGN.EQ.-1).AND.(SMWORK(I,1).GT.0.0)) ISIGNC=1
 820  CONTINUE
C
      DO 825 I=1,N
         if(ifrac.eq.1) then
            IF(ISIGNC.EQ.1) SMWORK(I,4)=max(0.1*ZEMAX,SMWORK(I,4))
         endif
         IF(SMWORK(I,4).LE.0.0) SMWORK(I,4)=ZEMIN
 825  CONTINUE
C
C SMOOTH
C
      CALL FILTR6(X,SMWORK(1,1),SMWORK(1,2),N,
     >         SMWORK(1,4),N,1.0,0,SMWORK(1,3),0,0.0,0,0.0,
     >         1.0,1.0,ZDUM,NDUM1,NDUM2)
C
C  INVERSION 2?
      IF(INVRT.EQ.1) THEN
         DO 830 I=1,N
            IF(ABS(SMWORK(I,2)).GT.ZSMALL) THEN
               SMWORK(I,2)=1/SMWORK(I,2)
            ELSE
               IF(SMWORK(I,2).LT.0.0) THEN
                  SMWORK(I,2)=-ZLARGE
               ELSE
                  SMWORK(I,2)=ZLARGE
               ENDIF
            ENDIF
 830     CONTINUE
      ENDIF
C
      RETURN
      END
C******************** END FILE SMOOF1.FOR ; GROUP SMOOF2 ******************
