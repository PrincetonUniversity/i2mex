!----------------------------------------------   RPCALC.blk    -  INCLUDE
!
! RPLOT CALCULATOR COMMON -- SEE ALSO            BLOCK DATA, RPCALDAT.FOR
!
!   *** dmc 8/08/99 ***
!   general comment:  when adding operators, there are additional tables
!   in rplot_sub/plcparse.for (fntab/nfn) that need to be modified in synch!
!   also in rplot_sub/plcfnxct.for one must add the appropriate jumps in
!   the large computed goto, and the corresponding target code to evaluate
!   the new operators!
!
!
!	LAST CHANGE:
!
!               8/04/99  dmc  breaking up plcfxt subroutine
!                             moving some plcfxt arrays into COMMON
!
!               6/02/95  tbt  Added FLXDIFF- Changed IMONAD to 37
!               6/01/95  tbt  Added ZONEC0 - changed IMONAD to 36
!               5/31/95  tbt  Added DIV -    changed IMONAD to 35
!               7/26/93  TbT  Added RPcalc4 (Inumerr, Index1)
!		9/29/92	 TBT  Changed IMONAD to 34 for dFdX.
!	        3/05/90  TBT  Changed IMONAD to 33 for TIMINT.
!              10/17/90  TBT  Changed IMONAD to 32 for ZONEB & ZONEC.
!              10/04/90  TBT  Added IPREC2.
!	       10/03/90  TBT  Changed IMONAD to 30 for SINH, TANH and +.
!	       10/02/90  TBT  Changed IMONAD from 16 to 27 for GRAD, VOLINT...
!		9/06/90  TBT  ADDED COMMON RPCALC3 TO PASS TO USE IN PLCHNDLR
!
!
!  MONADIC OPERATOR DESCRIPTORS -- SEE BLOCK DATA
!============
module rpcalc_mod
  implicit none

  INTEGER imonad,idyad,imaxop,nptacc,nscalar,iaexs,maxrpn, &
       maxcon,maxabb,ncmdrpp,nargrpp
  !============
  PARAMETER (IMONAD=38)
  CHARACTER*10 ZOP(IMONAD)
  CHARACTER*30 ZOLBL(IMONAD)
  INTEGER IFMT(IMONAD)
  !
  !  DYADIC OPERATOR DESCRIPTORS -- SEE BLOCK DATA
  PARAMETER (IDYAD=14)
  CHARACTER*10 ZOP2(IDYAD)
  CHARACTER*30 ZOLBL2(IDYAD)
  INTEGER IFMT2(IDYAD)
  INTEGER IPREC2(IDYAD)
  !
  !  MAX NUMBER OF OPERATORS, MAX NUMBER OF OPERANDS
  PARAMETER (IMAXOP=25)     !   tbt 8/90
  !
  CHARACTER*65 ZOPND(IMAXOP)  ! OPERANDS
  !
  INTEGER INOPO   ! NUMBER OF OPERATORS
  INTEGER IKOPO(IMAXOP)  ! CLASS OF OPERATORS (MONADIC/DYADIC)
  INTEGER ILOPO(IMAXOP)  ! NUMBER OF OPERATOR W/IN CLASS
  !
  INTEGER INOPND  ! NUMBER OF OPERANDS
  INTEGER IKOPND(IMAXOP) ! CLASS OF OPERANDS (CONSTANT, F(T), ..)
  INTEGER IAOPND(IMAXOP) ! ADDRESS OF OPERAND (IF A FUNCTION)
  REAL ZVOPND(IMAXOP)    ! VALUE OF OPERAND (IF A REAL CONSTANT)
  !
  LOGICAL  ILERROR       ! HAS AT LEAST ONE ARITHMETIC EXCEPTION OCCURRED
  !			       !     SINCE THE LAST INITIALIZATION.
  INTEGER  ICURPOS       ! CURRENT POSITION IN INPUT LINE.
  INTEGER  IERRPOS       ! POSITION IN THE INPUT LINE WHERE FIRST ERROR OCCURRED.

  !
  Parameter (Iaexs = 10)
  integer Index1, INumerr(IAEXS)
  real    ExpMax
  real	ExpMin
 
  !---------------------------------------------------------------------
  !  arrays previously local to plcfxt
  !
  PARAMETER  (MAXRPN=200)        ! MAXIMUM # OF TOKENS IN THE RPN LIST
  PARAMETER  (MAXCON=50)         ! MAXIMUM # OF CONSTANTS IN CONSTANT LIST
  INTEGER	   IRPNLIST(MAXRPN)  ! REVERSE POLISH NOTATION LIST  TBT 8/90
  !				       ! RETURNED FROM PARSE.
  INTEGER      INPOSLST(MAXRPN)  ! LIST OF POSITIONS IN INPUT LINE OF EACH
                                       !      TOKEN.
  DOUBLE PRECISION ZCONLIST(MAXCON)
  !                                       ! LIST OF CONSTANTS RETURNED FROM PARSE
  !
  PARAMETER  (MAXABB=10000)      ! max # of known symbol abbreviations
  CHARACTER*21 ALLABB(MAXABB)
  !
  logical lfndfunc
  !
  integer iposition(imaxop)      ! position in rpnlist of each operand
  integer ipos_str(imaxop)       ! position in input expr of each operand
  integer nops                   ! number of operands
  !
  integer ntypmr                 ! "Major Radius" axis id code
  logical Lrmajm                 ! "Major Radius" axis existence flag
  !
  integer lunrpc                 ! LUN for rplot calculator messages
  !
  !---------------------------------------------------------------------
  !  new dmc 9 Aug 1999 -- define rplot calculator "%"-command definition
  !  and parsing variables
  !
  parameter (ncmdrpp=20)         ! no. of commands currently defined.
  parameter (nargrpp=8)          ! max no. of args
  !
  character*512 cmdbuf           ! command buffer
  integer lcmdbuf                ! length of command buffer
  !
  integer kcmd                   ! selected command
  logical keyfmt                 ! keyword argument format flag
  integer klencmd                ! command str length (@pos=1 in buffer)
  integer kposarg(nargrpp)       ! arg position in buffer
  integer klenarg(nargrpp)       ! arg length in buffer
  !
  !  the following should be initialized in BLOCK DATA rplot_sub/rpcaldat.for
  integer icmdord(ncmdrpp)       ! **ORDERING** for help message (plcmds)
  logical lrepacc(ncmdrpp)       ! .TRUE. if command replaces accumulator
  character*10 rppcmds(ncmdrpp)  ! the command strings
  character*65 rppdescr(3,ncmdrpp) ! brief description of commands
  integer ncmdargs(ncmdrpp)      ! number of arguments
  character*10 rppkeys(nargrpp,ncmdrpp) ! argument keywords
  character*10 rppadfs(nargrpp,ncmdrpp) ! argument defaults
  character*20 exprtype(ncmdrpp) ! expression type constraints
  character*10 beforop(ncmdrpp)  ! monadic operator, if lrbefor(..)=.TRUE.
  !
  !  command characteristics
  !
  logical lrbefor(ncmdrpp)       ! .T:  modify expr with monadic operator
  logical lrafter(ncmdrpp)       ! .T:  exec after expr evaluation
  !
  character*1 gdchar             ! "guard character" for calc. commands
  !                                      cf BLOCK DATA rpcaldat.for, default: %

end module rpcalc_mod
