c                                                               PLCPARSE.FOR
 
 
C	  LOGICAL FUNCTION PLCPARSE ( text,
 
        Subroutine PLCPARSE ( OLDtext,             ! Changed 1/24/91
     1 		   vartab,
     1 		   rpnlst,
     1 		   conlst,
     1 		   ierror )
c
CC     A	 INPOSLST,
CC     B         LFNDFUNC,   ! Returned T if a new function (#23-NFN) is found
C
C  the above former arguments are now set via RPCALC COMMON
C
c               PPPP     AAA    RRRR     SSSS   EEEEE
c               P   P   A   A   R   R   S       E
c               PPPP    AAAAA   RRRR     SSS    EEE
c               P       A   A   R   R       S   E
c               P       A   A   R   R   SSSS    EEEEE
c
c
c    RUN TIME FUNCTION PARSING AND EXECUTING ROUTINES FOR FORTRAN
C         - called by PLCFXT
c    ------------------------------------------------------------
c
c  Dave W. Smith		+++++++++++++++++++++++++++++++++++
c  Seaver Computer Center       +++  (c) 1979 to public domain  +++
c  Claremont, CA   91711	+++++++++++++++++++++++++++++++++++
c
c  Dave W. Smith
c  Lawrence Livermore Laboratory
c  P.O. Box 808  L-440
c  Livermore, CA  94550
c
c  version	  what
c
c    1.0	Initial version number for the sake of keeping the
c		various releases straight.
c
c    1.1	Allow "@" by itself to be a valid identifier.
c		This used by CALC to hold the result of the last
c	        invocation of FNXCT.
c
c    1.3	Convert arithmetic to double precision.  Introduce
c		function DISP to modify display mode.  (Done by Charles
c		Karney at PPL.)
c
c    1.4	Reconvert arithmetic to single precision for
c		data base usage (Mike Bell, PPL)
c
c    1.5	Convert appropriate variables to character type for
c		usage on the VAX (Mike Bell, PPL).
c
c    1.6	Add overflow handler, fix -2**-1 bug,
c		convert to G_FLOATING, make DISP work again
c
C    1.7        Code has been modified to accept $ as a character, the same
C               as a letter (A-Z).   TBT  8/09/90
C
C    1.8        ARRAY INPOSLST AND INCLUDE RPCALC PUT IN TO BE ABLE TO PRINT
C                OUT WHERE IN THE INPUT LINE AN ARITHMETIC EXCEPTION OCCURS.
C                THIS WILL BE USED IN THE HANDLER, PLCHNDLR.  TBT  9/6/90
C
C    1.9        TOOK OUT CALL TO HANDL   TBT  9/10/90
C
C    2.0        HANDLES -3**2 SINCE THIS IS REALLY -(3**2) NOT (-3)**2.
C		CHANGED HIERARCHY OF UNARY "-" FROM 5 TO 4.
C               CHANGED HIERARCHY OF "**" FROM 4 TO 5.
C               WE NOW FORCE UNARY MINUS, E.G. 3**-2, UNTO THE STACK to
C 		BE ABLE TO HANDLE -----PI.
C		CHANGE ORDER OF ** HANDLING, FROM LEFT-TO-RIGHT TO RIGHT-TO-LEFT
C                  i.e. A**B**C IS NOW = A**(B**C).  SEE TOKEN=9.
C
C    3.1        Put in NFN=38, Added SINH and TANH.   TBT 10/03/90
C    3.2        Put in argument LFNDFUNC.             TBT 10/05/90
C    3.5        Put in TIMINT, NFN=39                 TBT  2/28/91
C    3.6        Put in dFdX  , NFN=40                 TBT  9/29/92
C    3.7        Put in check for input # too big for  TBT  7/27/93
C               machine accuracy. (See 9080)
C    3.8        Put in div   , NFN=41                 tbt  5/31/95
C    3.9        Put in ZONE0   NFN=42                 tbt  6/01/95
C    4.0        Put in FLXDIFF NFN=43                 tbt  6/02/95
C
CCC	options/g_floating    TBT
      use rpcalc_mod

      implicit integer (a - y)
      implicit double precision (z)
c
      character*(*) OLDtext, vartab(*)
      Character*170    text     ! new to hold DIV(F) expansion
      integer rpnlst(*)
      double precision conlst(*)
c
c  description of arguments ...
c
c  OLDtext:	character string of undefined length containing expression
c		to be parsed, terminated by a semicolon.  converted to
c		upper case in processing.
c
c  vartab:	table of user's dummy variables.  format of this table can
c		be figured out by looking at the data statemtent for fntab,
c		below.  note that the characters in this table *must* be
c		upper case.
c
c  rpnlst:	in which parse returns rpn tokens for function.  later
c		used by fnxct.  this array is also of undefined length, but
c		note that the number of rpn tokens for a function cannot
c		exceed the number of symbols in the function.
C
C  INPOSLST:    ARRAY TO RETURN THE POSITION IN THE INPUT LINE OF EACH TOKEN.
c
c  conlst:	array (of type real) in which parse stores any literal
c		numbers found the in function.  size of array is undefined.
c
C  LFNDFUNC:    Logical returned T is one of the new (#23-NFT) functions is
C     	        found during the parse. Else false.
c
c  ierror:	set to zero if successful or points to offending
c		character if a parsing error occurs.
c
c  function result: true if no overflow, false otherwise
c
c  program limits ....
c
c       1.      20 constants max.
c       2.      binary functions (min,max) cannot be nested more than 10
c	       deep.
c
c  any of the above limits can be changed with a bit of programming.  but
c  please be careful.
c
c
c  required (supplied) subroutines ...
c
c	0.	PLCPARSE	routine to parse a mathematical function
c				and generate rpn token code.
c	1.	FASCII		converts the fortran representation of a
c				character to its ascii equivalent.
c				(ICHAR is used here)
c	2.	CHR_TOUPPER	converts a lower case character string
c				to upper case.
c	3.	PLCFNDSM	symbol table searching routine.
c	4.      PLCFNXCT	routine to execute an rpn token list
c
C	NOTE: PLC..... ROUTINES ARE NOW IN RPLOT_SUB LIBRARY. TBT
 
c  PLCPARSE and PLCFNXCT are the two routines visible directly to the user,
c  although nothing prevents him/her from using the other routines.
c
c  of these, ICHAR and UPCASE are *very* machine dependent, and must be
c  modified to port PARSE, et al.
c
c
c  -----------------------
c  local data declarations
c  -----------------------
c
c  the array string is used to build a character string to
c  be used for matching user defined symbols and internally
c  defined function names.
c
      character string*23, chr*1, nxtchr*1, chr_toupper*1
c
c  the following arrays are used as stacks for the rpn algorithm
C
c
      PARAMETER ( MAXSTACK=25 ) ! MAXIMUM DEPTH OF STACK
      integer stack(MAXSTACK), shier(MAXSTACK)
      INTEGER IPOSTACK(MAXSTACK)      ! POSITION IN INPUT LINE OF THIS TOKEN.
c
c  the following arrays are used as stacks to fake binary functions
c  (e.g.  min(x,y) ).  comtok( cpnt ) holds the token to be substituted
c  for the next comma.  comlvl( cpnt ) holds a level indicator to
c  check the validity of the comma.
c
      integer comtok(10), comlvl(10)
c
c  these flags are used when evaluating numeric constants
c
      logical valid, negate, unary
 
 
c
c  the following arrays are used to hold the table of internally
c  defined functions and the number of arguments required.
c
      parameter (nfn = 44)              ! number of supplied functions
		              ! NOTE: if NFN changes, change NFN+50 in PLCFNXCT.
C			      ! ALSO: Changes should be made to the help file
C			      !       by changing RPCALC.blk and RPCALDAT.for.
 
c
      character*8 fntab(nfn)
      integer nargs(nfn)
c
c  now we set up the function table. *note* that the order in which
c  a function appears in the table is critical, and must match the
c  decoding scheme in 'PLCFNXCT'.  Also see test for FNIDX at label 2044 below.
c  the rpn token for a function is 50 + (index of function in fntab).
C
C  dmc -- note also, in RPLOT context, RPCALC COMMON and rpcaldat.for
C  BLOCK DATA need to be updated to be consistent with this table...
c
      data fntab(1), nargs(1)/'ABS', 1/
      data fntab(2), nargs(2)/'INT', 1/
      data fntab(3), nargs(3)/'SQRT', 1/
      data fntab(4), nargs(4)/'EXP', 1/
      data fntab(5), nargs(5)/'LOG', 1/
      data fntab(6), nargs(6)/'LOG10', 1/
      data fntab(7), nargs(7)/'COS', 1/
      data fntab(8), nargs(8)/'SIN', 1/
      data fntab(9), nargs(9)/'TAN', 1/
      data fntab(10), nargs(10)/'PI', 0/
      data fntab(11), nargs(11)/'MIN', 2/
      data fntab(12), nargs(12)/'MAX', 2/
      data fntab(13), nargs(13)/'EE', 0/
      data fntab(14), nargs(14)/'ACOS', 1/
      data fntab(15), nargs(15)/'ASIN', 1/
      data fntab(16), nargs(16)/'ATAN', 1/
      data fntab(17), nargs(17)/'R2D', 0/
      data fntab(18), nargs(18)/'D2R', 0/
      data fntab(19), nargs(19)/'ATAN2', 2/
      data fntab(20), nargs(20)/'DISP', 1/
      data fntab(21), nargs(21)/'COSH', 1/
      data fntab(22), nargs(22)/'SECH', 1/
      DATA FNTAB(23), NARGS(23)/'ZONEB',1/
      DATA FNTAB(24), NARGS(24)/'SCALEN', 1/
      DATA FNTAB(25), NARGS(25)/'LOGDERIV',1/
      DATA FNTAB(26), NARGS(26)/'GRAD', 1/
      DATA FNTAB(27), NARGS(27)/'ZONEC',1/
      DATA FNTAB(28), NARGS(28)/'VOLINT', 1/
      DATA FNTAB(29), NARGS(29)/'FLXINT', 1/
      DATA FNTAB(30), NARGS(30)/'ARINT', 1/
      DATA FNTAB(31), NARGS(31)/'LINAVG', 1/
      DATA FNTAB(32), NARGS(32)/'VOLAVG', 1/
      DATA FNTAB(33), NARGS(33)/'RMSVAVG', 1/
      DATA FNTAB(34), NARGS(34)/'DILINAVG', 1/
      DATA FNTAB(35), NARGS(35)/'DIVOLAVG', 1/
      data fntab(36), nargs(36)/'SINH', 1/
      data fntab(37), nargs(37)/'TANH', 1/
      data fntab(38), nargs(38)/'TIMINT', 1/
      Data Fntab(39), Nargs(39)/'DFDX  ', 1/
      Data Fntab(40), Nargs(40)/'DIV   ', 1/
      Data Fntab(41), Nargs(41)/'ZONE0', 1/
      Data Fntab(42), Nargs(42)/'FLXDIFF', 1/
      Data Fntab(43), Nargs(43)/'D    ', 1/
      data fntab(NFN),nargs(NFN)/';', -1/
c
C
C
c  --------------
c  initialization
c  --------------
CCC	external handl
CCC	call lib$establish(handl)
CCC	PLCparse=.true.
c
c  convert the expression to upper case (not needed)
c
c	call upcase(text)
c
c  here we initialize all pointers and flags
c
c  we start at the beginning of the expression and assume no errors
c
      text = OLDtext     ! tbt
      txtpnt = 1
      ierror = 0
c
c  now we set flag to look first for a unary operator
c
      unary = .true.
	
c
c  now we reset the pointers used in the rpn algorith
c
      cpnt = 1
      rpnt = 1
      spnt = 0
c
c  reset the binary function stack pointer
c
      fpnt = 0
c
c  zero the parenthesis counts
c
      lparen = 0
      rparen = 0
c
c  Initialize function finder.
 
      LFNDFUNC = .false.
 
      call rpcaldat_exec
c  ---------------------------------------------------------
c  here to examine the next character of the input function.
c  ---------------------------------------------------------
c
c  note that the function is not scanned character by character,
c  but that each section of code advances the text pointer to
c  the character after the token which that section extracts.
c
 
 1000  	chr = chr_toupper(text(txtpnt:txtpnt))
 
c
c  if the character is a semicolon, we've reached the end of the line.
c
      if ( chr .eq. ';' ) goto 5000
c
c  otherwise, kick the pointer and continue
c
      txtpnt = txtpnt + 1
c
c we will also find the next character useful
c
      nxtchr = chr_toupper(text(txtpnt:txtpnt))
c
c  certain cases can be dealt with immediately.  for example, we
c  don't really care about spaces between tokens.
c
      if ( chr .eq. ' ' ) goto 1000
 
c
c  we are looking for one of two cases of tokens; those which we
c  classify as unary (e.g. unary plus & minus, constants, and
c  variables), and those which we classify as binary (e.g. most
c  arithmetic and relational operators).  there are strict rules
c  for the context in which these operators may appear.  for example,
c  it is forbidden for two binary operators to appear together.
c  since we know the context in which the last token was seen, we
c  can use this information to check the legality of the next token.
c
      if ( .not. unary ) goto 3000
c
c  -------------------------------------------------------
c  here to begin checking for unary operators and operands
c  -------------------------------------------------------
c
c  first we check for a left parenthesis
c
 2000 if ( chr .ne. '(' ) goto 2010
        lparen = lparen + 1
        token = 1
        hier = 0
        goto 4000
c
c  now we check for unary plus and minus.  unary plus we ignore.
c
 2010 if ( chr .ne. '+' ) goto 2020
        goto 1000
c
 2020 if ( chr .ne. '-' ) goto 2030
        token = 3
        hier = 4
        goto 4000
 
 
c  any unary operator found now must be followed by a binary
c  operator.  set the flag now to dispatch the next token.
c
 2030 unary = .false.
	
c  check to see if the token is a numeric constant.  if so,
c  extract the number and store it in the caller's constant
c  list.  the rpn token is 20 + (offset of constant).
c
      if ( .not.
     1      ( chr .ge. '0' .and. chr .le. '9' .or. chr .eq. '.' )
     1    ) goto 2040
c
        znum = 0.0d0
        zpoint = 0.1d0
        valid = .false.
        if ( chr .eq. '.' ) goto 2034
        znum = dble( ichar( chr ) - ichar( '0' ) )
        valid = .true.
c
 2032   chr = chr_toupper(text(txtpnt:txtpnt))
        if ( chr .eq. '.' ) goto 2033
        if ( chr .eq. 'E' .or. chr .eq. 'D' ) goto 2035
        if ( chr .lt. '0' .or. chr .gt. '9' ) goto 2039
        IF ( Log10(Abs(Znum)) .GT. ExpMax-1 ) GoTo 9080     ! Tbt
          znum = znum*10.d0 + dble( ichar( chr ) - ichar( '0' ) )
          valid = .true.
          txtpnt = txtpnt + 1
          goto 2032
c
 2033   txtpnt = txtpnt + 1
 2034   chr = chr_toupper(text(txtpnt:txtpnt))
        if ( chr .eq. 'E' .or. chr .eq. 'D' ) goto 2035
        if ( chr .lt. '0' .or. chr .gt. '9' ) goto 2039
          znum = znum + dble( ichar(chr) - ichar('0') ) * zpoint
          zpoint = zpoint * 0.1d0
          valid = .true.
          txtpnt = txtpnt + 1
          goto 2034
c
 2035   txtpnt = txtpnt + 1
        zexp = 0.d0
        negate = .false.
        chr = chr_toupper(text(txtpnt:txtpnt))
        if (.not. (( chr .eq. '-' ) .or. ( chr .eq. '+' ))) goto 2036
          if ( chr .eq. '-' ) negate = .true.
          txtpnt = txtpnt + 1
          chr = chr_toupper(text(txtpnt:txtpnt))
 2036   if (( chr .lt. '0' ) .or. ( chr .gt. '9' )) goto 9002
 2037   If (.Not. Negate  .And. Zexp .Ge. ExpMax-1) Goto 9080    ! Tbt
        zexp = zexp * 10.d0 + dble( ichar(chr) - ichar('0') )
        txtpnt = txtpnt + 1
        chr = chr_toupper(text(txtpnt:txtpnt))
        if (( chr .ge. '0' ) .and. ( chr .le. '9' )) goto 2037
        if ( negate ) zexp = -zexp
        IF ( Zexp .Gt. 0  .And.
     1         Log10(Abs(Znum))+Zexp .GT. ExpMax-1 ) GoTo 9080     ! Tbt
        znum = znum * 10.d0 ** int( zexp )
c
 2039   if ( .not. valid ) goto 9004
        if ( cpnt .gt. 20 ) goto 9006
        conlst( cpnt ) = znum
        token = 20 + cpnt
        hier = -1
        cpnt = cpnt + 1
        goto 4000
c
c  here we check to see if the token is a string.  if so, we
c  assemble into array 'string', and check to see if it's a
c  user defined symbol (dummy argument).  if not, we'll check
c  to see if it's an internally defined symbol (function).
c
2040   IF ( CHR .EQ. '$')  GO TO 20400
         if ( chr .lt. '@' .or. chr .gt. 'Z' ) goto 2050
20400   strpnt = 1
        string(strpnt:strpnt) = chr
        if ( chr .eq. '@' ) goto 2042
c
2041    chr = chr_toupper(text(txtpnt:txtpnt))
        if ( .not. (
     1        ( CHR .EQ. '$')  .OR.
     1        ( chr .ge. 'A' .and. chr .le. 'Z' ) .or.
     1        ( chr .ge. '0' .and. chr .le. '9' ) .or. chr.eq.'_'
     1      )) goto 2042
        strpnt = strpnt + 1
          if ( strpnt .gt. len(string) ) goto 9014
          string(strpnt:strpnt) = chr
          txtpnt = txtpnt + 1
          goto 2041
c
 2042   strlen = strpnt
c
c  see if the string matches one of the callers dummy variables.
c  if so, the rpn token is 41.
c
        call PLCfndsm ( string(1:strlen), vartab, varidx )
        if ( varidx .le. 0 ) goto 2044
          token = 41
          hier = -1
          goto 4000
c
c  the string isn't a user defined symbol, so we'll check to see if
c  it's an internally defined function.  if so, check to see how many
c  arguments it takes, and set up to look for the arguments.
c  the rpn token is 50 + (offset of function in fntab)
c
 2044   call PLCfndsm ( string(1:strlen), fntab, fnidx )
        if ( fnidx .le. 0 ) goto 9010
 
C       --------------------------------
C       Handle DIV(F) which is defined as FLXDIFF(SURF*(F))/DVOL  tbt 5/95
C       Handle by inserting right on stack and switch operator!!
        IF (fnidx .EQ. 40) Then
           fnidx = 42    ! FLXDIFF operator
             Call PLCGRAD(text, txtpnt)     ! Put FLXDIFF on stack and insert rest
          ENDIF                             ! (SURF*( , )/DVOL into text
C       --------------------------------
 
 
        IF ( (FNIDX .GE. 23  .AND.    ! Function is GRAD, VOLINT, ZONEC...
     1          FNIDX .LE. 35) .OR. FNIDX .eq. 41)  LFNDFUNC = .TRUE.
C                                      ZONE0
 
        token = 50 + fnidx
        hier = 6
        fnarg = nargs(fnidx)
        if ( fnarg .ne. 0 ) goto 2045
          goto 4000
 2045   unary = .true.
        if ( chr .ne. '(' ) goto 9020
        if ( fnarg .ne. 1 ) goto 2046
          goto 4000
 2046   if ( fnarg.ne.2) call errmsg_exit(
     1         '?PLCPARSE: Internal error in FNTAB')
        fpnt = fpnt + 1
        if ( fpnt .gt. 10 ) goto 9018
        comlvl( fpnt ) = lparen - rparen + 1
        comtok( fpnt ) = token
        goto 1000
c
c  if we've fallen through to here, we didn't find a unary operator.
c  time to throw in the towel.
c
 2050 goto 9030
c
c  --------------------------------------------
c  here to begin checking for binary operators.
c  --------------------------------------------
c
c  first we check for a right parenthesis
c  check to see that we don't have two many of them, then check to
c  see that we haven't closed a binary function too soon
c
 3000   if ( chr .ne. ')' ) goto 3010
        rparen = rparen + 1
        if ( rparen .gt. lparen ) goto 9040
        if ( fpnt .gt. 0 ) then
              if ( (lparen-rparen+1) .eq. comlvl(fpnt) ) goto 9065
        endif
        token = 2
        hier = 0
        goto 4000
c
c  a unary operator must follow any of the following binary
c  operators.  set the flag now to dispatch the next token.
c
 3010 unary = .true.
c
c  now we check for the traditional arithmetic binary operators.
c
      if ( chr .ne. '+' ) goto 3020
        token = 4
        hier = 2
        goto 4000
c
 3020 if ( chr .ne. '-' ) goto 3030
        token = 5
        hier = 2
        goto 4000
c
 3030 if ( chr .ne. '*' ) goto 3050
        if ( nxtchr .ne. '*' ) goto 3040
          txtpnt = txtpnt + 1
          token = 9
          hier = 5
          goto 4000
 3040   token = 6
        hier = 3
        goto 4000
c
 3050 if ( chr .ne. '/' ) goto 3060
        token = 7
        hier = 3
        goto 4000
c
 3060 continue
      if ( chr .ne. '\\' ) goto 3070
        token = 8
        hier = 3
        goto 4000
c
 3070 if ( chr .ne. '^' ) goto 3080
        token = 9
        hier = 5
        goto 4000
c
c  now we check for the relational operators
c
 3080 if ( chr .ne. '=' ) goto 3085
        token = 10
        hier = 1
        goto 4000
c
 3085 if (chr .ne. '!' ) go to 3090
        if ( nxtchr .ne. '=' ) goto 3090
          txtpnt = txtpnt + 1
          token = 12
          hier = 1
          goto 4000
c
 3090 if ( chr .ne. '<' ) goto 3130
 3110   if ( nxtchr .ne. '=' ) goto 3120
          txtpnt = txtpnt + 1
          token = 13
          hier = 1
          goto 4000
 3120   token = 11
        hier = 1
        goto 4000
c
 3130 if ( chr .ne. '>' ) goto 3160
 3140   if ( nxtchr .ne. '=' ) goto 3150
          txtpnt = txtpnt + 1
          token = 15
          hier = 1
          goto 4000
 3150   token = 14
        hier = 1
        goto 4000
c
c  a comma can be a binary function.  really...
c  it gets replaced with a token for an internally defined binary
c  function.
c
 3160 if ( chr .ne. ',' ) goto 3170
        if ( fpnt .eq. 0 ) goto 9050
        if ( ( lparen - rparen ) .ne. comlvl( fpnt ) ) goto 9050
          token = comtok( fpnt )
          fpnt = fpnt - 1
          hier = 1
          goto 4000
c
c  if we've fallen through to here, we didn't manage to find a
c  binary operator.  time to give up and go home.
c
 3170 goto 9070
c
c  -----------------------------------------------
c  here to apply dijkstra's algorithm to the token
c  -----------------------------------------------
c
c  hier .lt. 0 implies that the token is an operand.  move it
c  directly to the rpn list and branch back for the next token.
c  token 41 (load value from user) is a special case.  it needs
c  to have the index of the value in the user's symbol table
c  loaded next.
c
 4000 continue
      if ( hier .ge. 0 ) goto 4010
        rpnlst( rpnt ) = token
        INPOSLST(RPNT) = TXTPNT-1
        rpnt = rpnt + 1
        if ( token .ne. 41 ) goto 1000
          rpnlst( rpnt ) = varidx
          INPOSLST(RPNT) = TXTPNT-1
          rpnt = rpnt + 1
          goto 1000
c
c  parenthesis get special treatment.  if a left parenthesis,
c  push it onto the operator stack.  if a right parenthesis, pop
c  operators off of the operator stack and into the rpnlst until
c  a left parenthesis is found, then discard both parenthesis.
c
 4010 if ( token .ne. 1 ) goto 4020
        spnt = spnt + 1
        IF (SPNT .GT. MAXSTACK)  GO TO 9071
        stack( spnt ) = token
        IPOSTACK(SPNT) = TXTPNT-1
        shier( spnt ) = hier   !heir -- fixed typo, dmc 21 May 92
        goto 1000
c
 4020 if ( token .ne. 2 ) goto 4030
c
c  here with a right parenthesis.  if we don't find a left
c  parenthesis on the stack, go yell and scream
c
 4021 if ( spnt .eq. 0 ) goto 9040
        if ( stack( spnt ) .eq. 1 ) goto 4022
        rpnlst( rpnt ) = stack( spnt )
        INPOSLST(RPNT) = IPOSTACK(SPNT)
        rpnt = rpnt + 1
        spnt = spnt - 1
        goto 4021
 4022 spnt = spnt - 1
      goto 1000
c
c  here with an operator.  until its hierarchical value is
c  greater than or equal to that of the operator on the top
c  of the stack, or the stack is empty, we will pop operators
c  off of the stack and move them to the rpn list.
c
 4030 if ( spnt .eq. 0 ) goto 4040
      if ( hier .gt. shier( spnt ) ) go to 4040
      IF (TOKEN .EQ. 3)  GO TO 4040                 !FORCE ONTO STACK "-"
                                   ! OPERAND MUST BE IN RPNLST BEFORE "-"
                                   ! THIS ALLOWS MULTIPLE "-", E.G. "$-----2*3"
      IF ( HIER .EQ. SHIER( SPNT ) .AND. TOKEN .EQ. 9) GOTO 4040
                                             ! ** IS RIGHT TO LEFT OPERATOR
                                             ! i.e. A**B**C = A**(B**C)
        rpnlst( rpnt ) = stack( spnt )
        INPOSLST(RPNT) = IPOSTACK(SPNT)
        rpnt = rpnt + 1
        spnt = spnt - 1
        goto 4030
c
c  here to push the token onto the operator stack
c
 4040 spnt = spnt + 1
      IF (SPNT .GT. MAXSTACK)  GO TO 9071
      stack( spnt ) = token
      IPOSTACK(SPNT) = TXTPNT-1
      shier( spnt ) = hier
      goto 1000
c
c  -----------------------------------------------------
c  here at end of text to terminate dijkstra's algorithm
c  -----------------------------------------------------
c
c  here when the end of the input text has been reached.  first check
c  to see that we ended on the right foot (i.e. looking for a binary
c  operator), then check to see that the parenthesis were correctly
c  balanced.  if everything is still o.k. we move whatever remains on
c  the operator stack to the rpn list, mark end of list, and return.
c
 5000 if ( unary ) goto 9030
      if ( lparen .ne. rparen ) goto 9040
c
 5010 if ( spnt .eq. 0 ) goto 5020
        rpnlst( rpnt ) = stack( spnt )
        INPOSLST(RPNT) = IPOSTACK(SPNT)
        rpnt = rpnt + 1
        spnt = spnt - 1
        goto 5010
c
 5020 if ( rpnt .eq. 1 ) goto 9000
      rpnlst( rpnt ) = 0
      INPOSLST(RPNT) = 1
c
c       ------
      return
c       ------
c
c  --------------
c  error routines
c  --------------
c
c  somewhere along the line this turkey blew it.
c
c  error messages are arranged roughly in the order in which
c  they may occur in the program.
c
 9000 call ZERMSG('?PLCPARSE: Null function')
      goto 9999
c
 9002 call ZERMSG(
     1       '?PLCPARSE: Invalid scientific notation in constant')
      goto 9999
c
 9004 call ZERMSG(
     1        '?PLCPARSE: Poorly formed numeric constant')
      goto 9999
c
 9006 call ZERMSG('?PLCPARSE: Too many constants')
      goto 9999
c
 9010 call ZERMSG('?PLCPARSE: Undefined identifier: '
     1             //string(1:strlen))
      goto 9999
c
 9014 call ZERMSG('?PLCPARSE: Identifier too long: '
     1              // string)
      goto 9999
c
 9018 call ZERMSG(
     1       '?PLCPARSE: Binary functions nested too deeply')
      goto 9999
c
 9020 call ZERMSG(
     1       '?PLCPARSE: Function requires parenthesis: '
     1         // string(1:strlen))
      goto 9999
c
 9030 call ZERMSG('?PLCPARSE: Expected unary operator')
      goto 9999
c
 9040 call ZERMSG('?PLCPARSE: Unbalanced parenthesis')
      goto 9999
c
 9050 call ZERMSG('?PLCPARSE: Misplaced comma')
      goto 9999
c
 9065 call ZERMSG(
     1       '?PLCPARSE: Incorrect number of arguments for function')
      goto 9999
c
 9070 call ZERMSG('?PLCPARSE: Expected binary operator')
      goto 9999
c
 9071   CALL ZERMSG('?PLCPARSE: Stack limits exceeded ')
      GO TO  9999
 
 9080   Call Zermsg('?PLCParse: Input number too large ')
      Go to 9999       ! Tbt
 
 
c  common exit routine for error messages
c
 9999 ierror = txtpnt
      rpnlst( 1 ) = 0
c
c       ------
      return
c       ------
      end
C------------------------
