C                                             PLCFIXUP.FOR IN GROUP RPLOT_SUB
 
      SUBROUTINE PLCFIXUP (IDATA, NPTS)
 
C	SUBROUTINE PLCFIXUP FIXES UP ILLEGAL FLOATING POINT NUMBERS GENERATED
C                           BY ARITHMETIC EXCEPTIONS WITH 0.
 
C                  PLCFIXUP IS CALLED BY PLCHNDLR.FOR            TBT   7/90
 
 
      INTEGER IDATA(NPTS)               ! ACTUALLY REAL ARRAY.
      INTEGER ILLOP                     ! ILLEGAL FLT PT NUMBER
#ifdef __F90
      data illop / Z"00008000" /
#else
      data illop / '00008000'X /
#endif
 
C	----------------------------------------------------------
 
      DO 100 I=1,NPTS
         IF (IDATA(I) .EQ. ILLOP)  IDATA(I) = 0
 100  CONTINUE                          ! END DO I
 
      RETURN
      END
