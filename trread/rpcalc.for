      subroutine rpcal0(zexpr,iwarn,ier)
c
c  calculator call, no results returned...
c
      character*(*) zexpr               ! calculator expression (input)
c
      integer iwarn                     ! arithmetic warning, 0 = OK
      integer ier                       ! completion code, 0 = OK
c
c-------------------------
c
      real zbuf(1)
      integer ibufsize,iret,istype
c
      ibufsize=0
      call rpcalc(zexpr,zbuf,ibufsize,iret,istype,iwarn,ier)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rpcalc(zexpr,zbuf,ibufsize,iret,istype,iwarn,ier)
C
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod
C
C  fetch an rplot calculator result.
C  data only, use rptime_p or rptime_s for timebase.
C  use rpdims for dimensioning information.
C
      character*(*) zexpr               ! calculator expression (input)
C
C  zexpr can also contain a command...
C
      real zbuf(*)                      ! buffer into which to write result
      integer ibufsize                  ! presumed size of buffer
C
C  *** if ibufsize = 0 *** then return no data, but do carry out the
C  rplot calculator calculation and do set error and warning flags
C  and do set istype
C
      integer iret                      ! number of data points written
      integer istype                    ! data type code of result
C
C  istype=-1 -- scalar f(t)
C  istype= 1 -- f(x,t), x is TRANSP zone ctrs if this is TRANSP data
C  istype= 2 -- f(x,t), x is TRANSP zone bdys if this is TRANSP data
C  ..etc..
C  istype = 0 indicates no expression was evaluated ( & iret=0 also).
C
      integer iwarn                     ! arithmetic warning, 0 = OK
      integer ier                       ! completion code, 0 = OK
C
C  ier=1 means:  error in the calculator, messages were written.
C  ier=2 means:  insufficient buffer space
C
C  iwarn.ne.0 means:  arithmetic errors were trapped during evaluation.
C
C  iret=0 if an error occurs.
C
C-------------------------------------------------------------
C
      integer idims(8)
C
C
      logical lshift,ieval,icmdexec,ioutpt
C
      character*10 zxabb(8)
      character*512 zinput,zcmd2,zxpres
C
C  this variable keeps track of status of RPLOT accumulator, i.e.
C  what type of profile is stored there.
C
      save istat
C
      data ict/0/
C-------------------------------------------------------------
C
C  initially assume no output
C
      iret=0
      istype=0
      ioutpt=.false.
C
      zinput=zexpr
      call trcaps(zinput)
C
C  a cobble to force trprofil etc to be linked.
C  for some reason, on SUNs, external stmts are insufficient.
C
      ilenx=len_trim(zinput)
      ilevn=2*(ilenx/2)
      if(ilenx.eq.ilevn) then
         ict=ict+1
      else
         ict=ict-1
      endif
      if(ict.eq.2147483647) then
C  this will never happen...
         call rpcalnk0(zbuf)
      endif
C
C  preparse:  look for commands (loop back to 100 because sometimes a 2nd
C  command is implied by an assignment; e.g.
C    tmp3 = %time_avg(0.1,ne*te)
C  results in a TIME_AVG command and a SAVE command to create TMP3.
C
 100  continue
      ier=0
      ilenx=len_trim(zinput)
C
      call plcprep(gdchar,zinput,zcmd2,zxpres,ieval,icmdexec,ier)
      if(ier.ne.0) go to 900
C
C F(T) DATA CHECK
      CALL DMGFOTX(2,IPT,IER)
      IF(IER.NE.0) RETURN
C
C  get workspaces...
      CALL DMDLOC('%WRK1',IND1,ISIZ1,IWRK1)
      CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)
C
      iwarn=0
C
      if(.not.ieval) go to 105
C
C  parse expression
C
      call plcparsr(zinput, ier)
      if(ier.ne.0) then
         ier=1
         go to 900
      endif
C
C  check/fetch input data
C
      call plckin(istat, Nscalar, ier)
      if(ier.ne.0) then
         ier=1
         go to 900
      endif
C
C  evaluate expression
C
      call plceval(istat,iwrk2,zinput,ilenx,ier,iwarn,lshift)
      if(ier.ne.0) then
         ier=1
         go to 900
      endif
C
      call rpdims(istat,irank,idims,zxabb,ier)
      if(ier.ne.0) then
         ier=3
         go to 900
      endif
C
      ioutpt=.true.
C
 105  continue
      if(icmdexec) call plcexec(istat,iwrk1,iwrk2,ipt,ier)
      if(ier.ne.0) go to 900
C
C  follow up command/expression
      if(zcmd2.ne.' ') then
         zinput=zcmd2
         go to 100
      endif
C
      if(.not.ioutpt) go to 900
C
      isize=1
      do ir=1,irank
         isize=isize*idims(ir)
      enddo
      if((ibufsize.gt.0).and.(ibufsize.lt.isize)) then
         call rpbufsiz('rpcalc',ibufsize,isize)
         ier=2
         go to 900
      endif
C
      IPF=iwrk2
C
      do i=1,isize
         if(ibufsize.gt.0) zbuf(i)=datbuf(ipf+i-1)
         if(istat.eq.-1) datbuf(Nptacc+i-1)=datbuf(ipf+i-1)
      enddo
      if(ibufsize.gt.0) iret=isize
      istype=istat
C
 900  continue
      return
      end
C------------------------------------
      subroutine rpcalnk0(zbuf)
      real zbuf(*)
C
      character*64 zlbl
      character*32 zuns
C
C  this routine should never be called.  It is here to assure that
C  trprofil, trscalar, and trfunid are linked...
C
      ier=-99                           ! if MDS+ tree is opened, leave open
      call trfunid(' ',' ','12345A00','FOO',itype,ier)
C
      call trscalar(' ',' ','12345A00','FOO',1,zlbl,zuns,
     >   itimes,zbuf,zbuf(2),ier)
C
      call trprofil(' ',' ','12345A00','FOO',1,1,zlbl,zuns,
     >   itype,inx,itimes,zbuf,zbuf(2),ier)
C
      return
      end
 
