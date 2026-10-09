C                                              PLCHNDLR.FOR IN GROUP RPLOT_SUB
 
      Subroutine PLC_Smry (IDATA, IPTS, ZINPUT, IWARN)
 
 
C	----------------------------------------------  ENTRY POINT PLC_Smry
C	.PRINT OUT THE Summary OF # OF ERRORS
 
C	    INPUT ARGS:  IDATA(IPTS) - ARRAY TO CHECK FOR ILLEGAL VALUES
C			 ZINPUT - CHARACTER STRING CONTAINING INPUT LINE.
C          OUTPUT ARGS:  IDATA       - ARRAY WITH 0 INSTEAD OF ILLEGAL VALUES.
C                        IWARN - =1 if there were any errors; 0 otherwise.
C
C	INTEGER FUNCTION PLCHNDLR IS IN FILE SOURCE:PLCHNDLR.FOR  TBT  7/90
C                                 HANDLES ARITHMETIC EXCEPTIONS FOR RPLOT
C                                 CALCULATOR. STORES 0.0 IS ANY RESULTANT THAT
C                                 HAS AN ARITHMETIC EXCEPTION.
C         Entry Point PLCerror     is called by routines PLCdiv,mult,add,sub
C                                   to simulate error handling & record keeping.
C	  ENTRY POINT PLC_HAND_INI IS SUPPLIED TO INITIALIZE # OF ERRORS TO 0.
C	  ENTRY POINT PLC_Smry IS SUPPLIED TO PRINT OUT #'S ERROR MESSAGES.
C	  SUBROUTINE  PLCFIXUP IS CALLED TO ZERO OUT ILLEGAL VALUES PUT INTO THE
C                              RESULTANT AFTER AN ARITHMETIC EXCEPTION OCCURS.
 
C            	         PLCHNDLR IS ESTABLISHED IN PLCFXT.FOR - THE CALCULATOR.

      use datmgr_mod
      use rpcalc_mod

      LOGICAL  LXOUT
      LOGICAL ILERR2
      Integer Itype
      INTEGER IDATA(IPTS)
      CHARACTER*(*) ZINPUT
      CHARACTER*1   HYPHEN
      CHARACTER*1   CARROT
 
      CHARACTER*32  CARTHEX(IAEXS)
 
      data hyphen/'-'/
      data carrot/'^'/
      data carthex /
     1	          'FLOATING PT OVERFLOWS',
     2		  'FLOATING DIVIDES BY ZERO     ',
     3		  'FLOATING POINT UNDERFLOWS    ',
     4		  'UNDEFINED EXPONENTIATIONS    ',
     5		  'LOGS OF ZERO OR NEG VALUE    ',
     6		  'SQ ROOTS OF NEGATIVE VALUE   ',
     7		  'SIGNIFICANCE LOSSES IN MATH LIB',
     8		  'FLOATING OVERFLOWS IN MATH LIB',
     9            'FLOATING UNDERFLOWS IN MATH LIB',
     A		  'INVALID ARGS IN MATH LIB     '   /
C
C-------------------------------------------------------------
C
      iwarn = 0
C
      ILAST = LEN(ZINPUT)
      LUNTRM = lunzer(0)                ! not LNOUT(); avoid UREAD dependence
      ILERR2 = .FALSE.
	
      IF (INDEX1 .NE. 0)  THEN
C		IF (LXOUT(0)) THEN
         WRITE (LUNTRM,8201) CARTHEX(INDEX1)
 8201    FORMAT(/ ' THE FIRST ERRORS WERE ', A32)
C		END IF  ! LXOUT
      END IF                            ! INDEX1
 
      DO 10 I=1,IAEXS
         IF (INUMERR(I) .NE. 0)  THEN
            ILERR2 = .TRUE.
C	      IF(LXOUT(0)) THEN
 
            WRITE(LUNTRM, 8001) CARTHEX(I), INUMERR(I)
 8001       FORMAT(/' ?PLC_Smry: # OF ', A32, ' =', I8/
     1         '     PARTIAL RESULT SET TO ZERO AT BAD POINTS.')
C	      END IF  ! LXOUT
            INUMERR(I) = 0.
         END IF                         ! INUMERR
 10   CONTINUE                          ! END DO I
 
      IF (ILERROR) THEN
C	    .WRITE OUT WHERE IN THE INPUT LINE THE FIRST ERROR OCCURRED.
C---------------------------
 
C     IF (LXOUT(0)) THEN
         IF ( ILAST .LE. 70)  THEN
            WRITE(LUNTRM,8126) ZINPUT
 8126       FORMAT(/' ?PLC_Smry: POSITION OF',
     1         ' FIRST ERROR IN INPUT LINE:'
     1         /1X, A)
            WRITE(LUNTRM,8127) (HYPHEN, I=3,IERRPOS), CARROT
 8127       FORMAT(1X, 79A1)
         ELSE                           ! INPUT IS MORE THAN ONE LINE
            WRITE(LUNTRM,8126) ZINPUT(1:70)
            IF (IERRPOS .LE. 72)  THEN
C		    .ERROR IS IN FIRST PART OF INPUT.
               WRITE(LUNTRM,8127) (HYPHEN, I=3,IERRPOS),
     1            CARROT
               WRITE(LUNTRM,8128) ZINPUT(71:ILAST)
            ELSE
C		    .ERROR IS IN SECOND PART OF INPUT.
               WRITE(LUNTRM,8128) ZINPUT(71:ILAST)
               WRITE(LUNTRM,8127) (HYPHEN, I=73,IERRPOS),
     1            CARROT
 8128          FORMAT(1X, A)
            END IF                      ! IER
 
         END IF                         ! ILAST
         WRITE (LUNTRM,8129)
 8129    FORMAT(/' ?PLC_Smry: THE RPLOT ACCUMULATOR "$" ',
     1      'HAS UNDEFINED VALUES')
C     END IF  ! LXOUT
C---------------------------
 
      END IF                            !  ILERROR
 
      IF (ILERR2)  THEN            ! ERROR HAS OCCURRED. PAUSE FOR READING MSGS
         CALL PLCFIXUP (IDATA, IPTS)    ! ZERO OUT ILLEGAL VALUES.
         iwarn = 1
      END IF                            ! ILERR2
 
      RETURN
      END
