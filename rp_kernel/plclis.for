C-------------------------------------------------------------------
C  PLCLIS  - LIST AVAILABLE OPERATORS, RPLOT CALCULATOR
C
      SUBROUTINE PLCLIS
 
C	Last changed:
C            5/31/95  tbt Added DIV operator
C	    10/02/90  TBT Added text for arithmetic exception.
C
      use rpcalc_mod
C
C-----------------------------------------------------------------
 
      LUNTRM=LUNZER(0)
C
 1000   FORMAT(5X,A10,5X,A, 10x, I10)
C
      WRITE(LUNTRM,1001)
 1001   FORMAT(/'%PLCLIS:' /
     >' MONADIC OPERATORS with USAGE FORMAT:  operator operand')
C
      DO I=1,IMONAD
        IF(IFMT(I).EQ.1) THEN
          WRITE(LUNTRM,1000) ZOP(I),ZOLBL(I)
        ENDIF
      ENDDO
 
C
      WRITE(LUNTRM,1002)
 1002   FORMAT(/
     >' MONADIC OPERATORS with USAGE FORMAT:  operator(operand)'/
     >'   Arithmetic exceptions are handled internally and      '/
     >'   reported at the end of the entire calculation.        ')
C
      DO I=1,18
        IF(IFMT(I).EQ.2) THEN
          WRITE(LUNTRM,1000) ZOP(I),ZOLBL(I)
        ENDIF
      ENDDO
 
      DO I=19,IMONAD
        IF(IFMT(I).EQ.2) THEN
          WRITE(LUNTRM,1000) ZOP(I),ZOLBL(I)
        ENDIF
      ENDDO
C
      WRITE(LUNTRM,1003)
 1003   FORMAT(///
     >' DYADIC OPERATORS with USAGE FORMAT:  operand operator operand',
     A  '    Precedence')
C
      DO I=1,IDYAD
        IF(IFMT2(I).EQ.1) THEN
          WRITE(LUNTRM,1000) ZOP2(I),ZOLBL2(I),IPREC2(I)
        ENDIF
      ENDDO
 
      WRITE(LUNTRM,1100)
 1100 FORMAT(/ ' Equal dyadic operators are parsed left to right ',
     1           'except for "**".'/
     2           ' Examples: A+B*C/D is A+((B*C)/D) but A*B**C**D is ',
     3           'A*(B**(C**D)).'/
     4           ' Relational operations result in 1. for TRUE or 0. ',
     5           'for FALSE.'/
     6           ' The result of (PI=PI=1.) is 1. (TRUE) and'/
     7           ' the result of (1.=PI=PI) is 0. (FALSE) since'/
     8           ' the equation is parsed left to right, i.e. '/
     9           ' (PI=PI=1.) is (PI=PI)=1. which is (1.=1.) which',
     A           ' is 1. (TRUE).'/)
 
      WRITE(LUNTRM,1004)
 1004   FORMAT(//
     >' DYADIC OPERATORS with ',
     >  'USAGE FORMAT:  operator( operand, operand )')
C
      DO I=1,IDYAD
        IF(IFMT2(I).EQ.2) THEN
          WRITE(LUNTRM,1000) ZOP2(I),ZOLBL2(I)
        ENDIF
      ENDDO
 
      WRITE(LUNTRM,1101)
 1101 FORMAT(/9x,' For Arc Tangent:  RESULT = ATAN2(ARG1,ARG2)'/
     1   9x, '    If ARG1>0 then RESULT>0'/
     2   9x, '    If ARG1=0 then if      ARG2>0 then RESULT =  0'/
     3   9x, '                   else if ARG2<0 then RESULT = PI'/
     4   9x, '    If ARG1<0 then RESULT<0'/
     5   9x, '    If ARG2=0 then ABS(RESULT) = PI/2.'/
     6   9x, '    If ARG1=0 and ARG2=0 then ERROR'/
     7   9x, '    In all cases:  -PI <= RESULT <= PI.'//
     8   9x, ' For min & max:  a preparse allows more than 2 args,'/
     9   9x, '    i.e. min(a,b,c,d) becomes min(a,min(b,min(c,d)))'/
     1   //)
C
      RETURN
      END
