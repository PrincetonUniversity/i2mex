C                                                          PLCFNXCT.FOR
c  this function takes a list of function tokens in reverse Polish
c  notation (rpn) from parse (SEE PLCPARSE) and treats them like
c  commands to a simple stack hand calculator. It returns values for
c  each point on the X axis for one time.
c  This subroutine is called by PLCFXT.
c
      SUBROUTINE PLCfnxct(                                   ! TBT 8/90
     1                    result, rpnlst, conlst,
     2                    Iarglst, prec, ARGLST, INX,
     3		          IT, ISTAT, ZINPUT,
     4                    error, IERR, LSHIFT)
C
c	BEWARE: BEWARE : IMPLICIT INTEGER(A-Z)
C------------------------------------------------
c
c  the calling parameters are the same as those used in PLCPARSE,
c  with the exception of arglst, which is a real array of values
c  corresponding to the dummy arguments supplied to parse.
c  Also prec is the argument which holds the arg given to disp. No longer used.*
C
c  Arguments:
c	Input:
c	    RPNLIST  - I*4 array - reverse Polish notation tokens.
c	    INPOSLST - I*4 array - the position of the token in the input line
c				   Used in error message.
c	    CONLST   - R*8 array - list of constants from PLCPARSE.
c	    IARGLIST - I*4 array - indirect pointer to argument in ARGLST.
c           PREC     - Dummy argument - not used anymore. From original FNXCT.
c	    ARGLST   - R*4 array(# of points in x direction,# of dummy arguments
c           INX      - I*4       - # of points in the X direction to calculate.
c	    IT       - I*4       - ITth time. Print error messages for IT=1.
c	    ZINPUT   - C*128     - input string to parser (equation). Used to
c		                   print place of error if any.
c	    LFNDFUNC - Logical   - T if a "new" type function was found during
c                                  the parse step ("new" = GRAD, VOLINT, etc..)
c	Output:
c	    RESULT   - R*8 array - Computed results at each of the INX points
c				    on the X axis.
c           ISTAT    - I*4       - Type of the X-axis after calulation.
c	    ERROR    - Logical   - T if an error (non arithmetic exception -
c                                   which are handled seperately in the handler
c				     (PLCHNDLR)) has occurred.
c	    IERR     - I*4       - Error code if ERROR = .TRUE.
c
c          Lshift    - warning logical, set .true. if zone/bdy shift was
c                      forced...
c
c
C	RECENT CHANGES:
C            6/02/95  tbt  Added FLXDIFF (92)
C            6/01/95  TBT  Added ZONE0   (91)- see plcparse.
C            5/31/95  tbt  Added DIV
C	     1/28/94  TBT  Moved 1001 to end of subroutine to avoid warning.
C	     7/27/93  TBT  Added calls to PLCmult,sub,add,div,power...
C            4/02/93  TBT  Added call to Plctype for "-".
C	    11/06/92  TBT  Handler of arithmetic exceptions was messing
C                          up indices in arrays. Putting in subroutines
C                          PlcAdd, PlcSub, PlcMult & PlcDiv.    Tbt
C	     9/30/92  TBT  Added DFDX.
C            9/09/92  TBT  Looking at error handling on / by zero.
C            6/11/92  TBT  Changed TBTLocal to have R*8 values first.
C	     2/28/91  TBT  Added TIMINT.
C	    10/05/90  TBT  Added LFNDFUNC to call list.
C	    10/01/90  TBT  Added code for functions # 75-85.
C			   Added checks for TYPE .ne. 1 = error.
C			   Added calls to PLZCENTR.
C	    09/27/90  TBT  Added code for functions 73-77 (ZONEB, ZONEC, GRAD,
C			   SINT AND VINT). Defined TEMPSTCK.
C                          Added IT to the argument list.
C			   Added call to PLCXCALC.
C			   Added ISTAT to argument list.
C			   Added ITYPE.
C		           Added LERR.
C	    09/26/90  TBT  Made call from PLCFXT once per time. Changed STACK(I)
C                          to STACK(NR0,I). Added all the DO loops over INX.
C			   Put in check of variable type.
C			   Put in INCLUDE CPLOTR to get NR0.
C	    09/24/90  TBT  ADDED MAXSTACK TO KEEP STACK FROM BECOMING TOO LARGE.
c	    06/01/90  TBT  Original version, FNXCT.for gotten from Jane Murphy.
c                          FNXCT is used in the VAX calculator, CALC.
c	-----------------------------------------------------------------------
 
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod

      implicit integer (a-z)
 
      integer rpnlst(*)
CC      INTEGER INPOSLST(*)     ! INPUT - POSITION IN INPUT LINE OF TOKEN.
C  INPOSLST now in COMMON
 
      CHARACTER*(*) ZINPUT    ! 128 Characters long.
      double precision conlst(*), result(NR0)
      REAL  arglst(INX,*)
      INTEGER IARGLST(*)
      logical error
C	LOGICAL ILPRINT		! True to print error message in PLCZCNTR the
                                !    first time we change ZBNDY to ZCENTR.
CC      Logical LFNDFUNC        ! Defined in PLCPARSE as T if a new function
			        !    has been found in the parse string.
C
C  INPOSLST and LFNDFUNC now in RPCALC COMMON
C
      LOGICAL LERR, LSHIFT
 
      PARAMETER (MAXSTACK=25) ! Maximum # of values kept in the stack.
 
      real*8, allocatable :: stack(:,:)       ! NR0 = # of points in X axis
      INTEGER,allocatable :: ITYPE(:)         ! Data type of each STACK level.
      REAL,allocatable :: TEMPSTCK(:,:)       ! Use as a temporary stack.
      Real,allocatable ::  dXaxis(:)          ! delta X
      real*8,allocatable :: wkarry(:)

      Double Precision Sixm1, Stemp
      real value, eps
 
      real zdtl                               ! for time int (dmc 9/96)
 
      data eps/1.0e-5/
 
      DOUBLE PRECISION  FINTEGRAL,  FNM1, FN
 
      Double precision temp0 (NR0)            ! 9/9/92 TBT temporary.
      Common /DUMBO/   Temp0
      Common /DumboI/  IXX
c
c  dmc -- COMMON rearranged to assure alignment of all variables
c    26 May 1992
c
      COMMON /TBTLOCAL/
     1                    FINTEGRAL(NR0),         !    TBT  2/28/91
     2                    FNM1(NR0),
     3			  TOKEN, PSPNT, MSPNT,    ! GET VALUES STORED NOT
     4                    SPNT,  RPNT             !   OPTIMIZED OUT.   TBT 9/90
c
 
      integer Nout      ! Dummy output unit
      Data    Nout  /0/
c     ----------------------------------------------------------------------
c     move arrays from automatic stack allocation to the heap
c
      allocate(stack(NR0,MAXSTACK ), ITYPE(MAXSTACK))
      allocate(TEMPSTCK(NR0,2),dXaxis(NR0),wkarry(NR0))
c
c  initialize pointers and flags
c
 
      lshift = .false.
c
1     error = .false.
      IERR = 0
      rpnt = 1
      spnt = 0
c
c  extract and dispatch the next function token
c
 1000 token = rpnlst( rpnt )
      IF (SPNT .GE. MAXSTACK)  GO TO 9250
      PSPNT = SPNT + 1   ! TBT
      MSPNT = SPNT - 1   ! TBT
      ICURPOS = INPOSLST(RPNT)     ! ESTABLISH CURRENT POSITION IN INPUT LINE
				     !  SO PLCHNDLR (HANDLER) KNOWS WHERE WE ARE
			             !  IN CASE OF AN ERROR.
      rpnt = rpnt + 1
      if ( token .eq. 0 ) goto 1001
      if ( token .lt. 0 .or. token .gt. 93 ) goto 900
                                   ! Note 93 = nfn+50 in PLCPARSE  <======================
      goto (
     1 900, 900,   3,   4,   5,   6,   7,   8,   9,  10,
     1  11,  12,  13,  14,  15, 900, 900, 900, 900, 900,
     1  21,  21,  21,  21,  21,  21,  21,  21,  21,  21,
     1  21,  21,  21,  21,  21,  21,  21,  21,  21,  21,
     1  41,  41,  41,  41,  41,  41,  41,  41,  41,  41,
     1  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,
     1  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
     1  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,
     1     81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
     1     91,  92,  93,                                   900), token
c
c  here with a bad token.  flag it and scream
c
900   call ZERMSG('?PLCFNXCT: Unknown token')
      IERR = ICURPOS
      go to 9999
 
 9250 CALL ZERMSG(' ?PLCFNXCT: Too many nested operands')
      IERR = ICURPOS
      GO TO 9999
 
 9260 CONTINUE
      CALL ZERMSG(' ?PLCFNXCT: X-axis conflict between operands')
      IERR = ICURPOS
      GO TO 9999
 
 9270 CALL ZERMSG(' ?PLCFNXCT: Subexpression is not on zone center')
      IERR = ICURPOS
      GO TO 9999
 
 
 9280 CALL ZERMSG(' ?PLCFNXCT: Subexpression is not on zone boundary')
      IERR = ICURPOS
      GO TO 9999
 
 9290 CALL ZERMSG(
     1   ' ?PLCFNXCT: Subexpression is not on zone center or boundary.')
      IERR = ICURPOS
      GO TO 9999
 
 
 9890 CALL ZERMSG(
     a        ' ?PLCFNXCT: Subexpression is not on Major radius x axis')
      IERR = ICURPOS
      GO TO 9999
 
 9891 CALL ZERMSG(
     A      ' ?PLCFNXCT: Major radius X-axis not found')
      IERR = ICURPOS
      GO TO 9999
 
 9892 CALL ZERMSG(
     A  ' ?PLCFNXCT: Subexpression is not a profile, cannot difference')
      IERR = ICURPOS
      GO TO 9999
 
9999  error = .true.
      goto 8888                                            ! <--------- RETURN
 
 
C	------------------------------------------------------------
 
c
c  ---------------
c  token execution
c  ---------------
c  the following sections of code processes the function tokens
c
c unary -, +, -, *, /, \, ^
c
 3    DO 30 IX=1,INX
         STACK(IX,SPNT) = - STACK(IX,SPNT)
 30   CONTINUE                          ! IX
      goto 1000
c
 4    CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1   TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2   ZINPUT, IT, ICURPOS, LERR, LSHIFT)
      IF (LERR) GO TO 9260
      DO 40 IX=1,INX
Ctbt	    STACK(IX,MSPNT) = STACK(IX,MSPNT) + STACK(IX,SPNT)
          Call PlcAdd (STACK(IX,MSPNT), STACK(IX,SPNT))
 40   CONTINUE    ! IX
      spnt = MSPNT
      goto 1000
c
c
 5    CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 50 IX=1,INX
Ctbt	    STACK(IX,MSPNT) = STACK(IX,MSPNT) - STACK(IX,SPNT)
          Call PlcSub(STACK(IX,MSPNT), STACK(IX,SPNT))
 50   CONTINUE    ! IX
      spnt = MSPNT
      goto 1000
c
 6    CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 6000 IX=1,INX
Ctbt	    STACK(IX,MSPNT) = STACK(IX,MSPNT) * STACK(IX,SPNT)
          Call PlcMult(STACK(IX,MSPNT), STACK(IX,SPNT))
 6000 CONTINUE    ! IX
      spnt = MSPNT
      goto 1000
c
 7    CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 7000 IX=1,INX
Ctbt	    Temp0(IX) = STACK(IX,MSPNT) / STACK(IX,SPNT)
CDEBUG 	    Write (Nout,9156) ix, temp0(ix), stack(ix,mspnt),
CDEBUG     1                     stack(ix,spnt),
CDEBUG     1                     mspnt, spnt
 9156       FORMAT(1x, i3, 3d15.5, 2i8)
Ctbt	    STACK(IX,MSPNT) = Temp0(IX)
          Call PlcDiv(STACK(IX,MSPNT), STACK(IX,SPNT))
 7000 CONTINUE   ! IX
      spnt = MSPNT
      goto 1000
c
 8    CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 800 IX=1,INX
          STACK(IX,MSPNT) = mod( STACK(IX,MSPNT), STACK(IX,SPNT) )
 800  CONTINUE   ! IX
      spnt = MSPNT
      goto 1000
c
 9    CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 810 IX=1,INX
C	    if (mod( STACK(IX,SPNT), 1.d0) .ne. 0.d0) then
C		STACK(IX,MSPNT) = STACK(IX,MSPNT) ** STACK(IX,SPNT)
C	    else
C		STACK(IX,MSPNT) = STACK(IX,MSPNT) ** int(STACK(IX,SPNT))
C	    end if
 
          Call PlcPower (STACK(IX,MSPNT), STACK(IX,SPNT))
 
 810  CONTINUE		! IX
      spnt = MSPNT
      goto 1000
c
c  now for the relational operators (tests are relative unless the test
c  constant is zero.)
c EQ
 10   relate = 0
      CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 100 IX=1,INX
      if (STACK( IX,spnt) .ne. 0.) then
       value = abs(STACK( IX,MSPNT)-STACK( IX,spnt))/
     >       abs(STACK( IX,spnt))
       if ( value.le.eps ) relate = 1
      else
       if (STACK( IX,MSPNT) .eq. STACK( IX,spnt)) relate = 1
      end if
      STACK(IX,MSPNT) = relate
 100  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c LT
 11   relate = 0
      CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 110 IX=1,INX
      if (STACK( IX,spnt) .ne. 0.) then
       value = abs(STACK( IX,MSPNT)-STACK( IX,spnt))/
     >       abs(STACK( IX,spnt))
       if ( ( value.gt.eps ) .and.
     *	 ( STACK(IX,MSPNT) .lt. STACK(IX,SPNT) ) ) relate = 1
      else
       if (STACK( IX,MSPNT) .lt. STACK( IX,spnt)) relate = 1
      end if
      STACK(IX,MSPNT) = relate
 110  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c NE
 12   relate = 0
      CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 120 IX=1,INX
      if (STACK( IX,spnt) .ne. 0.) then
       value = abs(STACK( IX,MSPNT)-STACK( IX,spnt))/
     >       abs(STACK( IX,spnt))
       if ( value.gt.eps ) relate = 1
      else
       if (STACK( IX,MSPNT) .ne. STACK( IX,spnt)) relate = 1
      end if
      STACK(IX,MSPNT) = relate
 120  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c LE
 13   relate = 0
      CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 130 IX=1,INX
      if (STACK( IX,spnt) .ne. 0.) then
       value = abs(STACK( IX,MSPNT)-STACK( IX,spnt))/
     >       abs(STACK( IX,spnt))
       if ( (value.le.eps) .or.
     *	 (STACK(IX,MSPNT) .le. STACK(IX,SPNT)) ) relate = 1
      else
       if (STACK( IX,MSPNT) .le. STACK( IX,spnt)) relate = 1
      end if
      STACK(IX,MSPNT) = relate
 130  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c GT
 14   relate = 0
      CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 140 IX=1,INX
      if (STACK( IX,spnt) .ne. 0.) then
       value = abs(STACK( IX,MSPNT)-STACK( IX,spnt))/
     >       abs(STACK( IX,spnt))
       if ( ( value.gt.eps ) .and.
     *	 ( STACK(IX,MSPNT) .gt. STACK(IX,SPNT) ) ) relate = 1
      else
       if (STACK( IX,MSPNT) .gt. STACK( IX,spnt)) relate = 1
      end if
      STACK(IX,MSPNT) = relate
 140  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c GE
 15   relate = 0
      CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 150 IX=1,INX
      if (STACK( IX,spnt) .ne. 0.) then
       value = abs(STACK( IX,MSPNT)-STACK( IX,spnt))/
     >        abs(STACK( IX,spnt))
       if ( (value.le.eps) .or.
     * 	 (STACK(IX,MSPNT) .ge. STACK(IX,SPNT)) ) relate = 1
      else
       if (STACK( IX,MSPNT) .ge. STACK( IX,spnt)) relate = 1
      end if
      STACK(IX,MSPNT) = relate
 150  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c
c  here to load constants
c
 21   offset = token - 20
      spnt = PSPNT
      ITYPE(SPNT) = 0.      ! => constant.
      DO 210 IX=1,INX
          STACK(IX,SPNT) = conlst( offset )
 210  CONTINUE  ! IX
      goto 1000
 
c
c  here to load user arguments in place of dummy variables
c
 41   spnt = PSPNT
        IOP = IARGLST( RPNLST(RPNT))
      ITYPE(SPNT) = IKOPND(IOP)
 
 
C       .I think this check is not needed because of check after #85.
C	IF (LFNDFUNC .AND.  ITYPE(SPNT) .GT. 2)  GO TO 9260 ! If LFNDFUNC then
                                        ! type must be 1 or 2
      DO 410 IX=1,INX
          STACK(IX,SPNT) = ARGLST(IX,IOP)
 410  CONTINUE  ! IX
      rpnt = rpnt + 1
      goto 1000
 
c
c  here to perform parse defined functions
c
 51   DO 510 IX=1,INX
          STACK(IX,SPNT) = abs( STACK(IX,SPNT) )
 510  CONTINUE  ! IX
      goto 1000
c
 52   DO 520 IX=1,INX
          STACK(IX,SPNT) = aint( STACK(IX,SPNT) )
 520  CONTINUE  ! IX
      goto 1000
 
c
 53   DO 530 IX=1,INX
C	    STACK(IX,SPNT) = sqrt( STACK(IX,SPNT) )
          Call          PlcSqrt( STACK(IX,SPNT) )
 530  CONTINUE  ! IX
      goto 1000
 
c
 54   DO 540 IX=1,INX
C	    STACK(IX,SPNT) = exp( STACK(IX,SPNT) )
          Call          PLCexp( STACK(IX,SPNT) )
 540  CONTINUE  ! IX
      goto 1000
c
 55   DO 550 IX=1,INX
C	    STACK(IX,SPNT) = log( STACK(IX,SPNT) )
          Call          PLClog( STACK(IX,SPNT) )
 550  CONTINUE  ! IX
      goto 1000
c
 56   DO 560 IX=1,INX
C	    STACK(IX,SPNT) = log10( STACK(IX,SPNT) )
          CALL          PLClog10( STACK(IX,SPNT) )
 560  CONTINUE
      goto 1000
c
 57   DO 570 IX=1,INX
          STACK(IX,SPNT) = cos( STACK(IX,SPNT) )
 570  CONTINUE
      goto 1000
c
 58   DO 580 IX=1,INX
          STACK(IX,SPNT) = sin( STACK(IX,SPNT) )
 580  CONTINUE  ! IX
      goto 1000
c
 59   DO 590 IX=1,INX
          STACK(IX,SPNT) = tan( STACK(IX,SPNT) )
 590  CONTINUE  ! IX
      goto 1000
c
 60   spnt = PSPNT
      ITYPE(SPNT) = 0.      ! => constant.
      DO 600 IX=1,INX
          STACK(IX,SPNT) = atan2( 0.d0, -1.d0)
 600  CONTINUE  ! IX
      goto 1000
 
c
 61   CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 610 IX=1,INX
          STACK(IX,MSPNT) = min( STACK(IX,MSPNT), STACK(IX,SPNT) )
 610  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c
 62   CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 620 IX=1,INX
          STACK(IX,MSPNT) = max( STACK(IX,MSPNT), STACK(IX,SPNT) )
 620  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c
 63   spnt = PSPNT
      ITYPE(SPNT) = 0.      ! => constant.
      DO 630 IX=1,INX
          STACK(IX,SPNT) = exp(1.0d0)
 630  CONTINUE  ! IX
      goto 1000
 
c
 64   DO 640 IX=1,INX
C	    STACK(IX,SPNT) = acos( STACK(IX,SPNT) )
          Call          PLCacos( STACK(IX,SPNT) )
 640  CONTINUE  ! IX
      goto 1000
 
c
 65   DO 650 IX=1,INX
C	    STACK(IX,SPNT) = asin( STACK(IX,SPNT) )
          Call          PLCasin( STACK(IX,SPNT) )
 650  CONTINUE  ! IX
      goto 1000
 
c
 66   DO 660 IX=1,INX
          STACK(IX,SPNT) = atan( STACK(IX,SPNT) )
 660  CONTINUE  ! IX
      goto 1000
 
c
 67   spnt = PSPNT
      ITYPE(SPNT) = 0.      ! => constant.
      DO 670 IX=1,INX
          STACK(IX,SPNT) = 180.d0/atan2(0.d0,-1.d0)
 670  CONTINUE  ! IX
      goto 1000
 
c
 68   spnt = PSPNT
      ITYPE(SPNT) = 0.      ! => constant.
      DO 680 IX=1,INX
          STACK(IX,SPNT) = atan2(0.d0,-1.d0)/180.d0
 680  CONTINUE  ! IX
      goto 1000
 
c
 69   CALL PLCTYPE( STACK(1,MSPNT), STACK(1,SPNT),
     1                TEMPSTCK,  INX, ITYPE(MSPNT),  ! Check X-axis type.
     2                ZINPUT, IT, ICURPOS, LERR, LSHIFT)
          IF (LERR) GO TO 9260
      DO 690 IX=1,INX
C	    STACK(IX,MSPNT) = atan2( STACK(IX,MSPNT), STACK(IX,SPNT) )
          Call           PLCatan2( STACK(IX,MSPNT), STACK(IX,SPNT) )
 690  CONTINUE  ! IX
      spnt = MSPNT
      goto 1000
 
c
 70   CALL ZERMSG(' ?PLCFNXCT: PREC operator undefined in RPLOT')
      ierr = icurpos
      GO TO 9999
C	prec = int( STACK(IX,SPNT) )
C	STACK(IX,SPNT) = arglst(27)
C	goto 1000
 
c
 71     DO 710 IX=1,INX
C	  STACK(IX,SPNT) = cosh( STACK(IX,SPNT) )
        Call          PLCcosh( STACK(IX,SPNT) )
 710  CONTINUE  ! IX
      goto 1000
 
c
 72     DO 720 IX=1,INX
C	  STACK(IX,SPNT) = 1./cosh( STACK(IX,SPNT) )
        Call            PLCOcosh( STACK(IX,SPNT) )
 720  CONTINUE  ! IX
      goto 1000
c
 
C	.Change to zone boundaries.                          ZONEB
 73   continue
      if (.not.nltransp) then
         call zermsg(' ?invalid operator for non-TRANSP data.')
         go to 9270
      endif
      IF( ITYPE(SPNT) .EQ. 2) GO TO 1000 ! Aready on zone boundaries.
      IF( ITYPE(SPNT) .GT. 1) GO TO 9270  ! Not on zone centers
			                    ! ITYPE = 0 & -1 (constants&scalars)
					    !     will be changed to zone bndry
 
      DO 730 IX=1,INX
          TEMPSTCK(IX,1) = STACK(IX,SPNT)  ! dbl to single
 730  CONTINUE
 
      CALL XINTZB( TEMPSTCK(1,1), TEMPSTCK(1,2), INX)
 
      DO 731 IX=1,INX
          STACK(IX,SPNT) = TEMPSTCK(IX,2)
 731  CONTINUE
 
      ITYPE(SPNT) = 2  ! Zone boundaries.
      GO TO 1000
 
 
C	.Change to zone centers                              ZONEC
 77   continue
      if (.not.nltransp) then
         call zermsg(' ?invalid operator for non-TRANSP data.')
         go to 9280
      endif
      IF( ITYPE(SPNT) .EQ. 1) GO TO 1000 ! Aready on zone centers.
      IF( ITYPE(SPNT) .GT. 2) GO TO 9280  ! Not on zone boundaries
	                                    ! ITYPE = 0 & -1 (const & scalrs)
					    !   will be changed to zone centers.
 
      DO 770 IX=1,INX
          TEMPSTCK(IX,1) = STACK(IX,SPNT)  ! dbl to single
 770  CONTINUE
 
      CALL XINTZC( TEMPSTCK(1,1), TEMPSTCK(1,2), INX)
 
      DO 771 IX=1,INX
          STACK(IX,SPNT) = TEMPSTCK(IX,2)
 771  CONTINUE
 
      ITYPE(SPNT) = 1    ! Zone centers
      GO TO 1000
 
 
c	New functions from RPLOT integro/calculator
C
 76   CONTINUE		    ! GRAD	iint = -1
 75   CONTINUE		    ! LOGDERIV	iint = -2
 74   CONTINUE		    ! SCALEN	iint = -3
 78   CONTINUE                    ! VOLINT	iint =  1
 79   CONTINUE		    ! FLXINT	iint =  2
 80   CONTINUE		    ! ARINT	iint =  3
 81   CONTINUE		    ! LINAVG	iint =  4
 82   CONTINUE		    ! VOLAVG	iint =  5
 83   CONTINUE		    ! RMSVAVG	iint =  6
 84   CONTINUE		    ! DILINAVG	iint =  7
 85   CONTINUE		    ! DIVOLAVG	iint =  8
 
      if (.not.nltransp) then
         call zermsg(' ?invalid operator for non-TRANSP data.')
         go to 9290
      endif
      IF (ITYPE(SPNT) .GT. 2) GO TO 9290 ! Must be zone center or boundary.
      IF (ITYPE(SPNT) .EQ. 2) then      ! Must be on zone centers-so convert
                                        ! Let ITYPE = 0&-1 (Const&scalar)
                                        !     pass through.
         lshift=.true.
         CALL PLCZCNTR( STACK(1,SPNT), TEMPSTCK, INX,
     2      ITYPE(SPNT),   ILPRINT)
 
      endif
 
      IINT = TOKEN - 77           ! iint = -3to-1 and 1-8
      CALL PLCXCALC( STACK(1,SPNT), TEMPSTCK, INX, IINT, IT)
      ITYPE(SPNT) = 2     ! Now on zone boundary
      GO TO 1000
 
c
 86     DO 860 IX=1,INX
C	  STACK(IX,SPNT) = sinh( STACK(IX,SPNT) )
        Call          PLCsinh( STACK(IX,SPNT) )
 860  CONTINUE  ! IX
      goto 1000
c
c
 87     DO 870 IX=1,INX
        STACK(IX,SPNT) = tanh( STACK(IX,SPNT) )
 870  CONTINUE  ! IX
      goto 1000
c
C                                                     TIMINT (Time integration)
 88     CONTINUE
C
C  dmc 3 Sep 96:  get correct results whether integrating scalar or
C  profile data
C
      if(istat.lt.0) then
         intl=ntt
      else
         intl=ntr
      endif
C
      IF (IT .EQ. 1)  THEN
C	    .Time 0.
          DO 880 IX=1,INX
      	FNM1(IX)       = STACK(IX,SPNT)
      	STACK(IX,SPNT) = 0.
      	FINTEGRAL(IX)  = 0.
 880        CONTINUE
 
      ELSE IF (IT .EQ. intl)  THEN
C	    .Last time.
          if(istat.lt.0) then
             zdtl=time(ntt)-time(ntt-1)
          else
             zdtl=time3(ntr)-time3(ntr-1)
          endif
          DO 883 IX=1,INX
      	STACK(IX,SPNT) = .5*(STACK(IX,SPNT)+FNM1(IX))
     1                             *zdtl
     2                             + FINTEGRAL(IX)
 883      CONTINUE
 
      ELSE
C	    .Time between T0 and last time.
          if(istat.lt.0) then
             zdtl=time(it)-time(it-1)
          else
             zdtl=time3(it)-time3(it-1)
          endif
          DO 885 IX=1,INX
      	FN = STACK(IX,SPNT)
      	STACK(IX,SPNT) = .5*(FN + FNM1(IX))
     1                             *zdtl
     2                             + FINTEGRAL(IX)
      	FNM1(IX) = FN                      ! Save off fct for next time
      	FINTEGRAL(IX) = STACK(IX,SPNT)     ! Save off INT for next time.
 885      CONTINUE  ! IX
      END IF
 
      goto 1000
c
 
C			                                        dF(x)/dx
   89 Continue
		                                  ! dF(x)/dx only for
						  !    x axis = major radius
      Call PLCGetX( It, dXaxis, inx, ItypMr, Ierr) ! Get delta X axis values
	
      if (.not.nltransp) then
         call zermsg(' ?invalid operator for non-TRANSP data.')
         go to 9890
      endif
      If (Itype(spnt) .ne. ITypMr)  go to 9890 !  Need a check of type.
      If (ierr .ne. 0)  go to 9891
      Temp0(1)   = (Stack(2,Spnt)  -Stack(1,Spnt))    /dXaxis(2)
      Temp0(Inx) = (Stack(Inx,Spnt)-Stack(Inx-1,Spnt))/dXaxis(Inx)
 
      InxM = Inx-1
      Do 890 Ix=2,InxM
          Temp0(Ix) = .5 *(
     1                (Stack(Ix,Spnt)-Stack(Ix-1,Spnt)) / dXaxis(Ix) +
     3		      (Stack(Ix+1,Spnt)-Stack(Ix,Spnt)) / dXaxis(Ix+1)  )
  890   Continue   ! Ix
 
      Do 891 Ix=1,Inx
          Stack(Ix,Spnt) = Temp0(Ix)
  891 Continue
 
      go to 1000
 
 
C			                 div(F) = (FLXDIFF(surf*F)/dvol)
 90   Continue
	
C       This is all handled in PLCPARSE by loading in FLXDIFF & Dvol in input.
      CALL ZERMSG(' ?PLCFNXCT: DIV Token found -error')
      ITYPE(SPNT) = 1		! Now on zone center
      call abortt
CCC	go to 1000
 
 
C	.Change to zone centers (assume 0.0 at center of plasma)        ZONE0
 91   continue
      if (.not.nltransp) then
         call zermsg(' ?invalid operator for non-TRANSP data.')
         go to 9280
      endif
      IF( ITYPE(SPNT) .EQ. 1) GO TO 1000 ! Aready on zone centers.
      IF( ITYPE(SPNT) .GT. 2) GO TO 9280  ! Not on zone boundaries
	                                    ! ITYPE = 0 & -1 (const & scalrs)
					    !   will be changed to zone centers.
      DO 970 IX=1,INX
          TEMPSTCK(IX,1) = STACK(IX,SPNT)  ! dbl to single
 970  CONTINUE
 
      CALL XINTZ0( TEMPSTCK(1,1), TEMPSTCK(1,2), INX)
 
      DO 971 IX=1,INX
          STACK(IX,SPNT) = TEMPSTCK(IX,2)
 971  CONTINUE
 
      ITYPE(SPNT) = 1    ! Zone centers
      GO TO 1000
 
 
                                                !              FLXDIFF
 92   Continue
C       . Flux differences
C	.Change to zone boundary (assume 0.0 at center of plasma)
      if (.not.nltransp) then
         call zermsg(' ?invalid operator for non-TRANSP data.')
         go to 9280
      endif
      IF( ITYPE(SPNT) .EQ. 2) GO TO 927    ! Aready on zone boundary
      IF( ITYPE(SPNT) .GT. 2) GO TO 9280   ! Not on zone boundaries
	                                     ! ITYPE = 0 & -1 (const & scalrs)
					     !   will be changed to zone boundary
      DO 924 IX=1,INX
          TEMPSTCK(IX,1) = STACK(IX,SPNT)  ! dbl to single
 924   CONTINUE
 
      CALL XINTZB( TEMPSTCK(1,1), TEMPSTCK(1,2), INX)  ! Put on zone boundary
 
      DO 926 IX=1,INX
          STACK(IX,SPNT) = TEMPSTCK(IX,2)
 926  CONTINUE
 927  Continue
 
      Sixm1 =         Stack(1,Spnt)
      Stack(1,Spnt) = Stack(1,Spnt) - 0.0     ! 0.0 at center of plasma
 
      DO 928 IX=2,INX
         Stemp =          Stack(Ix,Spnt)   ! save for next difference
         Stack(Ix,Spnt) = Stack(Ix,Spnt) - Sixm1
         Sixm1 = Stemp
 928  CONTINUE		! IX
 
      ITYPE(SPNT) = 1    ! Zone centers output
      goto 1000
 
 93   Continue
C  generic finite difference operator D(...)
      if(itype(spnt).le.0) go to 9892
      call PlcDiff(Stack(1,Spnt),WkArry,INX)
      go to 1000
C	-------------------------------------------------------
 
c
c  here at end of rpn function list.  value is on top of stack.
c
 1001 IF (SPNT .NE. 1) THEN
         CALL ZERMSG( ' ?PLCFNXCT: Algorithm error! SPNT .NE. 1')
         ERROR = .TRUE.
         call abortt
         goto 8888
      END IF   ! SPNT
 
      ! spnt = 1                      ! ITYPE(i) = -1  => scalar
      IF ( ITYPE(SPNT) .GT. 0)        ! ITYPE(i) =  0  => constant
     1             ISTAT = ITYPE(SPNT)  ! Keep the ISTAT specified in PLCKIN if
                                        !   this is a calculation of constants.
					!   or sclalars.
 
      DO 10000 IX=1,INX
          result(IX) = STACK(IX,SPNT)
10000  CONTINUE                 ! IX

c
c all exits
c        
 8888  continue
       deallocate(stack,itype,tempstck,dxaxis,wkarry)

      END
