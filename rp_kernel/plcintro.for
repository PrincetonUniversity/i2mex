      subroutine plcintro(zchar)
C
C  type out introductory blurb on the RPLOT calculator...
C
      character*1 zchar
C
      LUNT=LUNZER(0)
C
      WRITE(LUNT,1001)
 1001 FORMAT(/
     >' RPLOT Calculator (vsn 4.0 dmc/tbt 22 Sept 1999)'//
     >' Use the calculator by entering expressions and commands'/
     >' at the calculator prompt.  An EXPRESSION consists of'/
     >' a combination of operators and operands written in a'/
     >' traditional algebraic/functional syntax, such as:'/
     >'     0.5*ZONEB(NE)*(ZONEB(TE)**1.5)*GRAD(TE).'/
     >' The list of available operators is given below.  Operands'/
     >' may be any of the following:  numeric constants, "PI",'/
     >' "EE" (=exp(1.)), "R2D" (=180./PI), "D2R" (=PI/180.),'/
     >' RPLOT scalar functions of time, including "TIME",'/
     >' RPLOT profile functions of time, or, subexpressions.'/
     >' Scalar and profile functions may be user defined (see the'/
     >' section on COMMANDS).  The output of the expression evaluation'/
     >' is stored as a "temporary" RPLOT profile or scalar function,'/
     >' called the RPLOT ACCUMULATOR, which, once defined'/
     >' may be referenced as an operand by the symbol "$".'/
     >' Expression operands must be consistent as to profile'/
     >' type (i.e. "f*g" is not valid if "f" is a function of'/
     >' minor radius and "g" a function of major radius); the result'/
     >' of an expression involving functions of (x,t) is itself a'/
     >' function of (x,t).  "Zone-centered" and "Boundary-centered"'/
     >' profile flux functions can be mixed.  Constants and scalar'/
     >' functions can always be included in profile expressions.'/
     >' The expression evaluator works by processing profiles point'/
     >' by point in a loop over time.  Arithmetic errors are trapped'/
     >' with warnings.')
C
      WRITE(LUNT,1002) zchar,zchar,zchar,zchar,zchar,zchar
 1002 format(/
     >' Calculator COMMANDS enable access to functionality not'/
     >' available through the expression evaluator, such as:'/
     >'   * creation, naming and labeling of user defined functions'/
     >'   * time integration and differentiation'/
     >'   * mappings, e.g. flux zone <--> midplane major radius'/
     >'   * smoothing and time averaging'/
     >'   * interpolation.'//
     >' ----> all commands start with the "',a1,'" character.'//
     >' A complete list of commands is given below.  Typical command'/
     >' syntax looks like:'/
     >5x,a1,'SAVE(P_E,"electron pressure","J/cm3",1.5*1.601e-19*NE*TE)'/
     >' (this command saves the results of the expression'/
     >' "1.5*1.601e-19*NE*TE" as a user defined fuction with labels;'/
     >' the arguments are given in positional syntax).  Arguments'/
     >' can also be specified in a keyword based syntax, as in:'/
     >5x,a1,'TIME_AVG(DELTA_T=0.1, EXPR= NE*TE)'/
     >' which sets the RPLOT ACCUMULATOR to NE*TE with each time point'/
     >' time averaged over the range +/-DELTA_T.'/
     >' Implied SAVE commands can also be generated with the syntax:'/
     >'     TMP1=1.5*1.601e-19*NE*TE'/
     >'     P_E,"electron pressure","J/cm3" = ',a1,'TIME_AVG(0.1,TMP1)'/
     >' which creates user defined functions TMP1 and P_E.  Most'/
     >' command arguments have default values which are automatically'/
     >' assumed if arguments are omitted.  For example:'/
     >5x,a1,'SMOOTH(0.1,0.05,,,NE*TE)'/
     >' in this case the 3rd and 4th arguments are defaulted.  Keyword'/
     >' syntax can be used to avoid having to place the commas, as in:'/
     >5x,a1,'SMOOTH(DELTA_T=0.1,DELTA_X=0.05,NE*TE)'/)
C
      return
      end
