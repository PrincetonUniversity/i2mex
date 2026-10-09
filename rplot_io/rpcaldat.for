C-----------------------------------------------     RPCALDAT.for in RPLOT_SUB
      subroutine rpcaldat_exec
C
C	Last changed:
C            8/07/99  dmc added LUNRPC
C            6/02/95  tbt Added FLXDIFF - Bdy to Bdy  differences in fluxes.
C            6/01/95  tbt Added ZONE0 - interpolate to zone centers with 0.0 at center
C            5/31/95  tbt Added DIV
C	     9/29/92  TBT Added dFdX.
C            3/05/91  TBT Added TIMINT.
C   	    10/17/90  TBT Added ZONEB and ZONEC
C	    10/04/90  TBT Added IPREC2.
C	    10/02/90  TBT Added 11 ZOP functions. GRAD, VOLINT,...
C
C  LOAD RPLOT CALCULATOR COMMON BLOCK
C
      use rpcalc_mod
C
      ZOP = (/
     A           '=         ','-         ',
     1           '+         ',
     2           'SQRT      ','LOG       ',
     3           'LOG10     ','EXP       ',
     4           'SIN       ','COS       ',
     5           'TAN       ','ASIN      ',
     6           'ACOS      ','ATAN      ',
     7           'SINH      ','COSH      ',
     8           'TANH      ','SECH      ',
     8           'INT       ','ABS       ',
     9           'GRAD      ','LOGDERIV  ',
     A           'SCALEN    ','VOLINT    ',
     1           'FLXINT    ','ARINT     ',
     2           'LINAVG    ','VOLAVG    ',
     3           'RMSVAVG   ','DILINAVG  ',
     4           'DIVOLAVG  ','ZONEB     ',
     5           'ZONEC     ','TIMINT    ',
     6           'DFDX      ','DIV       ',
     7           'ZONE0     ','FLXDIFF   ',
     8           'D         '/)
 
      ZOLBL = (/
     >'assignment operator           ','unary minus                   ',
     >'unary plus                    ',
     >'square root                   ','natural (base e) logarithm    ',
     >'base 10 logarithm             ','exponential e**x              ',
     >'sin (argument in radians)     ','cos (argument in radians)     ',
     >'tan - overflow caution        ','arc sine                      ',
     >'arc cosine                    ','arc tangent                   ',
     >'hyperbolic sine               ','hyperbolic cosine             ',
     >'hyperbolic tangent            ','1. / hyperbolic cosine        ',
     >'integer value                 ','absolute value                ',
     >'gradient                      ','logarithmic derivative        ',
     >'scale length                  ','volume integral               ',
     >'ptcl/power flux integral      ','area integral                 ',
     >'line average                  ','volume average                ',
     >'RMS volume average            ','double inverse line average   ',
     >'double inverse volume average ','zone boundary (interpolate to)',
     A'zone centered (interpolate to)','time integration (T0=T(1))    ',
     B'dF(x)/dx (for RMAJR x axis)   ','div(F)=drav*grad(surf*F)/dvol ',
     C'zone centered (with 0. at 0.) ','Bdy-to-Bdy difference in fluxs',
     C'finite difference operator    '/)
 
      IFMT(1:3) = 1
      IFMT(4:38)= 2 ! Increase when new function is added <===========
C
      ZOP2 = (/
     >           '**        ','*         ',
     >           '/         ','+         ',
     >           '-         ','==        ',
     >           '!=        ','<         ',
     >           '<=        ','>         ',
     >           '>=        ','MIN       ',
     >           'MAX       ','ATAN2     '
     >              /)
      ZOLBL2 = (/
     >'exponentiation                ','multiplication                ',
     >'division                      ','addition                      ',
     >'subtraction                   ','relational equal              ',
     >'relational not equal          ','relational less than          ',
     >'relational less than or equal ','relational greater than       ',
     >'relational greater than or eq ','pt by pt MIN operator         ',
     >'pt by pt MAX operator         ','arc tan2                      '
     >   /)
 
      IFMT2 = (/ 1,1,1,1,1,1,1,1,1,1,1,2,2,2 /)
      IPREC2 = (/ 5,3,3,2,2,1,1,1,1,1,1,0,0,0 /) ! Precedence of operators.
C
      lunrpc = 6                    ! default LUN, RP calculator messages
C-----------------------------------------------------------------------
C  dmc 9 Aug 1999 --
C  these blocks define rplot calculator "%-commands"
C  put cmds in caps -- parsing is case-insensitive
C
C  modification notes (1) increase parameter ncmdrpp in inclshare/RPCALC
C                     (2) modify data statements
C
      rppcmds(1) = 'SAVE'
      rppcmds(2) = 'TIME_TRACE'
      rppcmds(3) = 'XINTERP'
      rppcmds(4) = 'DELETE'
      rppcmds(5) = 'TIME_DERIV'
      rppcmds(6) = 'RMJMAP'
      rppcmds(7) = 'XMAP'
      rppcmds(8) = 'SMOOTH'
      rppcmds(9) = 'TIME_AVG'
      rppcmds(10) = 'TIME_INT'
      rppcmds(11) = 'LABEL'
      rppcmds(12) = 'MINPRO'
      rppcmds(13) = 'MAXPRO'
      rppcmds(14) = 'XLOCATE'
      rppcmds(15) = 'RFETCH'
      rppcmds(16) = 'MG_CREATE'
      rppcmds(17) = 'MG_DELETE'
      rppcmds(18) = 'MG_ADDFUN'
      rppcmds(19) = 'MG_DELFUN'
      rppcmds(20) = 'GS2FETCH'
C
      lrepacc = (/
     >     .false.,.true.,.true.,.false.,.true.,.true.,.true.,.true.,
     >     .true.,.true.,.false.,.true.,.true.,.true.,.true.,.false.,
     >     .false.,.false.,.false.,.false./)
C
      icmdord = (/
     >   1,4,11,16,17,18,19,2,3,14,12,13,5,10,8,9,7,6,15,20/)
C
C  SAVE description
      rppdescr(1,1) = 'SAVE calculator result as a labeled function.'
      rppdescr(2,1) = ' '
      rppdescr(3,1) = ' '
C  TIME_TRACE description
      rppdescr(1,2) =
     >'for a given X axis index value (integer btw 1 and Nx)'
      rppdescr(2,2) =
     >'select the corresponding time trace.  OK to set INDEX to'
      rppdescr(3,2) =
     >'I<no. btw 0.0 and 1.0>; e.g I0.0 ->x(1), I1.0 -> x(Nx).'
C  XINTERP description
      rppdescr(1,3) =
     >'INTERPolate profile data to specified X value, yielding a'
      rppdescr(2,3) =
     >'scalar function of time.  X can itself be constant or a scalar'
      rppdescr(3,3) =
     >'fcn of time.  EXCEPTION = value to use if X out of bounds.'
C  DELETE description
      rppdescr(1,4) =
     >'DELETE user defined function & recover space.'
      rppdescr(2,4) =
     >'specify FCN as "*" to delete ALL user defined functions.'
      rppdescr(3,4) = ' '
C  TIME_DERIV description
      rppdescr(1,5) =
     >'evaluate TIME DERIVative of data'
      rppdescr(2,5) = ' '
      rppdescr(3,5) = ' '
C  RMJMAP description
      rppdescr(1,6) = 
     >'MAP profile vs flux coordinate to MaJor Radius.'
      rppdescr(2,6) = ' '
      rppdescr(3,6) = ' '
C  XMAP description
      rppdescr(1,7) =
     >'MAP profile vs. major radius to flux coordinate (X).'
      rppdescr(2,7) =
     >'specify TARGET as "CTR" or "BDY" for zone-centered or bdy-'
      rppdescr(3,7) =
     >'centered result.  Map from SIDE "IN", "OUT" or (avg) "INOUT"'
C  SMOOTH description
      rppdescr(1,8) =
     >'SMOOTH data, triangular weighting +/- DELTA_T and +/- DELTA_X'
      rppdescr(2,8) =
     >'set EPS_* to limit smoothing, EPS_X and EPS_T are additive,'
      rppdescr(3,8) =
     >'append "%" to value for EPS_* in per-cent, "R" for relative.'
C  TIME_AVG description
      rppdescr(1,9) =
     >'TIME AVeraGe data +/- DELTA_T (secs) e.g for simple smoothing'
      rppdescr(2,9) =
     >'for double inverse time averaging give DELTA_T.lt.0'
      rppdescr(3,9) = ' '
C  TIME_INT description
      rppdescr(1,10) =
     >'TIME INTegrate from specified time T0 (secs)'
      rppdescr(2,10) = ' '
      rppdescr(3,10) = ' '
C  LABEL description
      rppdescr(1,11) =
     >'Change LABELs of an existing function'
      rppdescr(2,11) = ' '
      rppdescr(3,11) = ' '
C  MINPRO description
      rppdescr(1,12) =
     >'Extract as a scalar function the MINimum PROfile value'
      rppdescr(2,12) =
     >'at each time.'
      rppdescr(3,12) = ' '
C  MINPRO description
      rppdescr(1,13) =
     >'Extract as a scalar function the MAXimum PROfile value'
      rppdescr(2,13) =
     >'at each time.'
      rppdescr(3,13) = ' '
C  XLOCATE description
      rppdescr(1,14) =
     >'In the profile X axis find the X value corresponding to'
      rppdescr(2,14) =
     >'the indicated test value (MIN, MAX, or a numeric constant)'
      rppdescr(3,14) =
     >'NOT_UNIQUE and NOT_FOUND specify exception return values.'
C  RFETCH description
      rppdescr(1,15) =
     >'Fetch a scalar or profile function from another runid'
      rppdescr(2,15) =
     >'and use it to create a named user defined function in'
      rppdescr(3,15) =
     >'the current session.'
C  MG_CREATE description
      rppdescr(1,16) =
     >'Create a multigraph association, providing a label and at least'
      rppdescr(2,16) =
     >'one member function.  Precede member id with a minus sign ("-")'
      rppdescr(3,16) =
     >'to make its additive inverse a member of the multigraph.'
C  MG_DELETE description
      rppdescr(1,17) =
     >'Delete a multigraph association.'
      rppdescr(2,17) = ' '
      rppdescr(3,17) = ' '
C  MG_ADDFUN description
      rppdescr(1,18) =
     >'Add one or more member functions to a multigraph.  Precede the'
      rppdescr(2,18) =
     >'member id with a minus sign ("-") to make its additive inverse'
      rppdescr(3,18) =
     >'a member of the multigraph.'
C  MG_DELFUN description
      rppdescr(1,19) =
     >'Delete one or more member functions from a multigraph.'
      rppdescr(2,19) = ' '
      rppdescr(3,19) = ' '
C  GS2FETCH description
      rppdescr(1,20) =
     >'Read GS2 tabulated results and create three profile functions'
      rppdescr(2,20) =
     >'User specifies path to data, x axis = "X" or "R", and a prefix'
      rppdescr(3,20) =
     >'for output <prefix>AKY, <prefix>OMEGA, <prefix>GAMMA profiles'
C
C---------------------------------------
      rppkeys = ' ' ! base default values
C  SAVE arguments
      ncmdargs(1) = 4
      rppkeys(1,1) = 'FCN_ID'
      rppkeys(2,1) = 'LABEL'
      rppkeys(3,1) = 'UNITS'
      rppkeys(4,1) = 'EXPR'
C  TIME_TRACE arguments
      ncmdargs(2) = 2
      rppkeys(1,2) = 'INDEX'
      rppkeys(2,2) = 'EXPR'
C  XINTERP arguments
      ncmdargs(3) = 3
      rppkeys(1,3) = 'X'
      rppkeys(2,3) = 'EXCEPTION'
      rppkeys(3,3) = 'EXPR'
C  DELETE arguments
      ncmdargs(4) = 1
      rppkeys(1,4) = 'FCN_ID'
C  TIME_DERIV arguments
      ncmdargs(5) = 1
      rppkeys(1,5) = 'EXPR'
C  RMJMAP arguments
      ncmdargs(6) = 1
      rppkeys(1,6) = 'EXPR'
C  XMAP arguments
      ncmdargs(7) = 3
      rppkeys(1,7) = 'TARGET'
      rppkeys(2,7) = 'SIDE'
      rppkeys(3,7) = 'EXPR'
C  SMOOTH arguments
      ncmdargs(8) = 5
      rppkeys(1,8) = 'DELTA_T'
      rppkeys(2,8) = 'DELTA_X'
      rppkeys(3,8) = 'EPS_T'
      rppkeys(4,8) = 'EPS_X'
      rppkeys(5,8) = 'EXPR'
C  TIME_AVG arguments
      ncmdargs(9) = 2
      rppkeys(1,9) = 'DELTA_T'
      rppkeys(2,9) = 'EXPR'
C  TIME_INT arguments
      ncmdargs(10) = 2
      rppkeys(1,10) = 'TO'
      rppkeys(2,10) = 'EXPR'
C  LABEL arguments
      ncmdargs(11) = 3
      rppkeys(1,11) = 'ITEM_ID'
      rppkeys(2,11) = 'LABEL'
      rppkeys(3,11) = 'UNITS'
C  MINPRO arguments
      ncmdargs(12) = 1
      rppkeys(1,12) = 'EXPR'
C  MAXPRO arguments
      ncmdargs(13) = 1
      rppkeys(1,13) = 'EXPR'
C  XLOCATE arguments
      ncmdargs(14) = 4
      rppkeys(1,14) = 'TEST_VALUE'
      rppkeys(2,14) = 'NOT_UNIQUE'
      rppkeys(3,14) = 'NOT_FOUND'
      rppkeys(4,14) = 'EXPR'
C  RFETCH arguments
      ncmdargs(15) = 4
      rppkeys(1,15) = 'PATH'
      rppkeys(2,15) = 'RUN_ID'
      rppkeys(3,15) = 'FCN_ID'
      rppkeys(4,15) = 'LOCAL_ID'
C  MG_CREATE arguments
      ncmdargs(16) = 8
      rppkeys(1,16) = 'PKG_ID'
      rppkeys(2,16) = 'LABEL'
      rppkeys(3,16) = 'FCN1'
      rppkeys(4,16) = 'FCN2'
      rppkeys(5,16) = 'FCN3'
      rppkeys(6,16) = 'FCN4'
      rppkeys(7,16) = 'FCN5'
      rppkeys(8,16) = 'FCN6'
C  MG_DELETE arguments
      ncmdargs(17) = 1
      rppkeys(1,17) = 'PKG_ID'
C  MG_ADDFUN arguments
      ncmdargs(18) = 8
      rppkeys(1,18) = 'PKG_ID'
      rppkeys(2,18) = 'FCN1'
      rppkeys(3,18) = 'FCN2'
      rppkeys(4,18) = 'FCN3'
      rppkeys(5,18) = 'FCN4'
      rppkeys(6,18) = 'FCN5'
      rppkeys(7,18) = 'FCN6'
      rppkeys(8,18) = 'FCN7'
C  MG_DELFUN arguments
      ncmdargs(19) = 8
      rppkeys(1,19) = 'PKG_ID'
      rppkeys(2,19) = 'FCN1'
      rppkeys(3,19) = 'FCN2'
      rppkeys(4,19) = 'FCN3'
      rppkeys(5,19) = 'FCN4'
      rppkeys(6,19) = 'FCN5'
      rppkeys(7,19) = 'FCN6'
      rppkeys(8,19) = 'FCN7'
C  GS2FETCH arguments
      ncmdargs(20) = 3
      rppkeys(1,20) = 'PATH'
      rppkeys(2,20) = 'XID'
      rppkeys(3,20) = 'PREFIX'
C
C  a blank default means no default.  expression type rule also given here.
      rppadfs = ' '  ! base default
C  SAVE argument defaults
      exprtype(1) = 'any'
      rppadfs(2,1) = 'Unknown'
      rppadfs(3,1) = 'Unknown'
      rppadfs(4,1) = '$'
      lrbefor(1) = .false.
      beforop(1) = ' '
      lrafter(1) = .false.
C  TIME_TRACE argument defaults
      exprtype(2) = 'any profile'
      rppadfs(2,2) = '$'
      lrbefor(2) = .false.
      beforop(2) = ' '
      lrafter(2) = .true.
C  XINTERP argument defaults
      exprtype(3) = 'any profile'
      rppadfs(2,3) = '-99.0'
      rppadfs(3,3) = '$'
      lrbefor(3) = .false.
      beforop(3) = ' '
      lrafter(3) = .true.
C  DELETE argument defaults
      exprtype(4) = 'N.A.'
      lrbefor(4) = .false.
      beforop(4) = ' '
      lrafter(4) = .false.
C  TIME_DERIV argument defaults
      exprtype(5) = 'any'
      rppadfs(1,5) = '$'
      lrbefor(5) = .false.
      beforop(5) = ' '
      lrafter(5) = .true.
C  RMJMAP argument defaults
      exprtype(6) = 'profile vs. X'
      rppadfs(1,6) = '$'
      lrbefor(6) = .false.
      beforop(6) = ' '
      lrafter(6) = .true.
C  XMAP argument defaults
      exprtype(7) = 'profile vs. Rmajor'
      rppadfs(1,7) = 'CTR'
      rppadfs(2,7) = 'INOUT'
      rppadfs(3,7) = '$'     
      lrbefor(7) = .false.
      beforop(7) = ' '
      lrafter(7) = .true.
C  SMOOTH argument defaults
      exprtype(8) = 'any'
      rppadfs(2,8) = '0.0'
      rppadfs(3,8) = '0.0'
      rppadfs(4,8) = '0.0'
      rppadfs(5,8) = '$'
      lrbefor(8) = .false.
      beforop(8) = ' '
      lrafter(8) = .true.
C  TIME_AVG argument defaults
      exprtype(9) = 'any'
      rppadfs(2,9) = '$'
      lrbefor(9) = .false.
      beforop(9) = ' '
      lrafter(9) = .true.
C  TIME_INT argument defaults
      exprtype(10) = 'any'
      rppadfs(1,10) = '0.0'
      rppadfs(2,10) = '$'
      lrbefor(10) = .true.
      beforop(10) = 'TIMINT'
      lrafter(10) = .true.
C  LABEL argument defaults
      exprtype(11) = 'N.A.'
      rppadfs(2,11) = '%unchanged'
      rppadfs(3,11) = '%unchanged'
      lrbefor(11) = .false.
      beforop(11) = ' '
      lrafter(11) = .false.
C  MINPRO argument defaults
      exprtype(12) = 'any profile'
      rppadfs(1,12) = '$'
      lrbefor(12) = .false.
      beforop(12) = ' '
      lrafter(12) = .true.
C  MAXPRO argument defaults
      exprtype(13) = 'any profile'
      rppadfs(1,13) = '$'
      lrbefor(13) = .false.
      beforop(13) = ' '
      lrafter(13) = .true.
C  XLOCATE argument defaults
      exprtype(14) = 'any profile'
      rppadfs(2,14) = '-222.0'
      rppadfs(3,14) = '-111.0'
      rppadfs(4,14) = '$'
      lrbefor(14) = .false.
      beforop(14) = ' '
      lrafter(14) = .true.
C  RFETCH argument defaults
      exprtype(15) = 'N.A.'
      rppadfs(1,15) = '.'
      lrbefor(15) = .false.
      beforop(15) = ' '
      lrafter(15) = .true.
C  MG_CREATE argument defaults
      exprtype(16) = 'N.A.'
      rppadfs(4,16) = '%empty'
      rppadfs(5,16) = '%empty'
      rppadfs(6,16) = '%empty'
      rppadfs(7,16) = '%empty'
      rppadfs(8,16) = '%empty'
      lrbefor(16) = .false.
      beforop(16) = ' '
      lrafter(16) = .false.
C  MG_DELETE argument defaults
      exprtype(17) = 'N.A.'
      lrbefor(17) = .false.
      beforop(17) = ' '
      lrafter(17) = .false.
C  MG_ADDFUN argument defaults
      exprtype(18) = 'N.A.'
      rppadfs(2,18) = '%empty'
      rppadfs(3,18) = '%empty'
      rppadfs(4,18) = '%empty'
      rppadfs(5,18) = '%empty'
      rppadfs(6,18) = '%empty'
      rppadfs(7,18) = '%empty'
      rppadfs(8,18) = '%empty'
      lrbefor(18) = .false.
      beforop(18) = ' '
      lrafter(18) = .false.
C  MG_DELFUN argument defaults
      exprtype(19) = 'N.A.'
      rppadfs(2,19) = '%empty'
      rppadfs(3,19) = '%empty'
      rppadfs(4,19) = '%empty'
      rppadfs(5,19) = '%empty'
      rppadfs(6,19) = '%empty'
      rppadfs(7,19) = '%empty'
      rppadfs(8,19) = '%empty'
      lrbefor(19) = .false.
      beforop(19) = ' '
      lrafter(19) = .false.
C  MG_DELFUN argument defaults
      exprtype(20) = 'N.A.'
      lrbefor(20) = .false.
      beforop(20) = ' '
      lrafter(20) = .false.
C------------------------------------------------------------------
C
      gdchar = '%'
C
      END
