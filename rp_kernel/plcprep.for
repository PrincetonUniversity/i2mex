      subroutine plcprep(zpre,zinput,zcmd2,zxpres,ieval,icmdexec,ier)
      use rpcalc_mod
C
C  dmc Aug 1999:  RPLOT calculator pre-parse
C
C    0.  ignore comments (starting with unquoted # character)
C    1.  detect commands
C    2.  convert syntax "name=expr" into
C        zpre//'SAVE(name,"unknown","unknown",expr) & detect as a command
C    3.  convert max(a,b,c) -> max(a,max(b,c))
C                min(a,b,c) -> min(a,min(b,c)) ... so dyadic parser works
C                note that a,b,c can be subexpressions of arb. complexity.
C
C  arguments...
C
      character*(*) zpre                ! input command prefix, NO BLANKS
C
C  if the first characters of zinput match zpre exactly, interpret input
C  line as a command
C
      character*(*) zinput              ! input/output...
C  on input zinput contains the user input expression or command
C  on output zinput contains any remaining expression to be evaluated
C     by the calculator parser/evaluator
C
      character*(*) zxpres              ! input/output
C
C  on input the expression that led to the current accumulator contents.
C  ou output:  the expression or cmd-expression combination to be evaluated
C  with numerical result -- i.e. what will be stored in the accumulator
C  after completion of all steps.
C
C  (usually an expression, but could also be TDERIV(expression)...)
C
      character*(*) zcmd2               ! output
C
C  a residual command, with deferred execution...
C
      logical ieval                     ! output .TRUE. if expression needs
C                                         ...evaluation
C
      logical icmdexec                  ! output .TRUE. to exec command
C                                         ...after expression evaluation
C
      integer ier                       ! output status code, 0=OK
C
C  examples:  given zpre = '%'
C
C  on input:  zinput='ne*te'
C  on output:  zinput unchanged
C              zxpres = 'ne*te'
C              zcmd2 = ' '
C              ieval = T
C              icmdexec = F
C              ier = 0
C
C  on input:  zinput= 'tmp = ne*te'
C  on output:  zinput = 'ne*te'
C              zxpres = 'ne*te'
C              zcmd2 = ' '
C              ieval = T
C              icmdexec = T    & command is SAVE(tmp,"TMP",,$)
C              ier = 0
C
C  on input:  zinput= '%delete(tmp)'
C  on output:  zinput= unchanged
C              zxpres unchanged
C              zcmd2 = ' '
C              ieval = F
C              icmdexec = T    & command is to delete TMP user fcn.
C              ier = 0
C
C  on input:  zinput= 'tmp = %time_deriv(ne*te)'
C  on output:  zinput = 'ne*te'
C              zxpres = '%time_deriv(ne*te)'
C              zcmd2 = '%SAVE(tmp,,,$)'   -- deferred evaluation of SAVE (tmp)
C              ieval = T
C              icmdexec = T    & evaluate %tderiv
C              ier = 0
C
C-----------------------------------------------------------------------
C
      character*512 zinbuf,zwkstr
C
C
      character*10 zlhs
C
      character*1 zquot
C
      logical idebug
      data idebug/.false./
C-----------------------------------------------------------------------
C
      ier=0
      zcmd2=' '
      ieval=.false.
      icmdexec=.false.
C
      lunt=lunzer(0)
      ilpre=len_trim(zpre)
      if(ilpre.le.0) then
         call zermsg(' ?plcprep:  call error, zero length prefix.')
         ier=1
         go to 1000
      endif
C
C  look at input line
C
      ilzin=len_trim(zinput)
      if(ilzin.le.0) go to 1000         ! no content in zinput
C
C  drop anything after the "#" comment character if any
C
      zquot=' '
      do ic=1,ilzin
         if(zquot.ne.' ') then
            if(zinput(ic:ic).eq.zquot) zquot=' '
         else if((zinput(ic:ic).eq.'"').or.(zinput(ic:ic).eq.'''')) then
            zquot=zinput(ic:ic)
         else if(zinput(ic:ic).eq.'#') then
            zinput(ic:)=' '
            go to 5
         endif
      enddo
C
 5    continue
      ilzin=len_trim(zinput)
      if(ilzin.le.0) go to 1000         ! no content in zinput
C
C  left shift zinput if necessary
C
      do ic=1,ilzin
         if((zinput(ic:ic).ne.' ').and.(zinput(ic:ic).ne.char(9))) then
            go to 10
         endif
      enddo
      go to 1000                        ! no content in zinput
C
 10   continue
      if(ic.gt.1) then
         ishift=ic-1
         ilzin=ilzin-ishift
         do is=1,ilzin
            zinput(is:is)=zinput(is+ishift:is+ishift)
         enddo
         zinput(ilzin+1:)=' '
      endif
C
C  check for command syntax
C
      iassign=0
      if(zpre(1:ilpre).eq.zinput(1:ilpre)) go to 100
C
C  check for A=B syntax
C  ... look for = sign
C
      ieqs=0
      icommas=0
      zquot=' '
      do ic=1,ilzin
         if(zquot.eq.' ') then
            if(zinput(ic:ic).eq.'"') then
               zquot='"'
            else if(zinput(ic:ic).eq.'''') then
               zquot=''''
            else if(zinput(ic:ic).eq.',') then
               icommas=icommas+1
            else if(zinput(ic:ic).eq.'=') then
               ieqs=ic
               go to 20
            endif
         else
            if(zinput(ic:ic).eq.zquot) then
               zquot=' '                ! closing quote was detected
            endif
         endif
      enddo
C
      if(ieqs.eq.0) then
         ieval=.true.
         zxpres=zinput
         go to 500                      ! zinput contains just an expression
      endif
C
C  check validity of LHS...
C
 20   continue
      if(ieqs.eq.1) go to 29
C
      ilhs=len_trim(zinput(1:ieqs-1))
      icomm1=index(zinput(1:ilhs),',')-1
      if(icomm1.le.0) icomm1=ilhs

      call idchek_setadj(zinput(1:icomm1),iadj)

      call idchek(zinput(1:icomm1),idnum,iclass,iuser,iadj)
      ier=0
      if((iclass.eq.1).and.(iuser.eq.0)) ier=1 ! can't assign:  file scalar
      if((iclass.eq.2).and.(iuser.eq.0)) ier=2 ! can't assign:  file profile
      if(iclass.eq.3) ier=3             ! can't assign:  multigraph name
      if(iclass.gt.3) ier=4
      if(ier.ne.0) then
         if (ier.eq.1) call zermsg(' ?plcprep:  run defined scalar.')
         if (ier.eq.2) call zermsg(' ?plcprep:  run defined profile.')
         if (ier.eq.3) call zermsg(' ?plcprep:  multigraph.')
         call zermsg(' ?plcprep:  cannot assign to '//zinput(1:icomm1)//
     >      ':  name in use.')
         go to 1000
      endif
C
C  OK looks valid; look for nonblank RHS
C
      do ic=ieqs+1,ilzin
         if((zinput(ic:ic).ne.' ').and.(zinput(ic:ic).ne.char(9))) then
            go to 25
         endif
      enddo
C
      ier=1
      call zermsg(' ?plcprep, RHS is empty:  "'//zinput(1:ieqs)//' "')
      go to 1000
C
 25   continue
      icmd2=0
      if(zinput(ic:ic+ilpre-1).eq.zpre(1:ilpre)) icmd2=1
      icrhs=ic
      ilrhs=ilzin-icrhs+1
C
      zinbuf=zpre(1:ilpre)//
     >   'SAVE('//zinput(1:ilhs)//','
      ilzbuf=len_trim(zinbuf)
C
C  accept labeling strings if they are present; otherwise insert
C  commas to get default labeling
C
      if(icommas.eq.0) then
         zinbuf(ilzbuf+1:ilzbuf+4+ilhs)='"'//zinput(1:ilhs)//'",,'
         ilzbuf=ilzbuf+4+ilhs
      else if(icommas.eq.1) then
         zinbuf(ilzbuf+1:ilzbuf+1)=','
         ilzbuf=ilzbuf+1
      endif
C
      if(icmd2.eq.0) then
         if((ilzbuf+ilrhs).gt.len(zinput)) then
            write(lunt,9001) zinput(1:min(80,ilzin))
 9001       format(' ?plcprep:  expression too long in ASSIGN stmt:'/
     >         1x,a)
            ier=1
            go to 1000
         endif
C  perform assignment A=expr as SAVE(A,...,expr) command
         zcmd2=' '
         zinbuf(ilzbuf+1:)=zinput(ic:ilzin)//')'
         zinput=zinbuf
         ilzin=len_trim(zinput)
      else
C  perform assignment as (evaluation) command followed by SAVE
         iassign=1
         zinbuf(ilzbuf+1:ilzbuf+2)='$)'
         zcmd2=zinbuf                   ! deferred SAVE command
         zinbuf=zinput(icrhs:ilzin)     ! RHS is command to evaluate next
         zinput=zinbuf
         ilzin=len_trim(zinput)
      endif
      go to 100                         ! parse as a command
C
 29   continue
      call zermsg(' ?plcprep:  item left of = sign invalid:  '//
     >   zinput(1:ieqs))
      ier=1
      go to 1000
C
C---------------------------------------
C
C  OK parse a command line
C
 100  continue
      istart=ilpre+1
      zinbuf=zinput(istart:ilzin)
      call rppcpars(zpre(1:ilpre),zinbuf,zinput,zxpres,ieval,icmdexec,
     >   ier)
      if(icmdexec.and.(iassign.eq.1)) then
         if(.not.lrepacc(kcmd)) then
            call zermsg(
     >' ?plcprep:  cannot assign, command on RHS has no numeric output')
            ier=1
         endif
      endif
      if(ieval) then
         go to 500
      else
         go to 1000
      endif
C
C---------------------------------------
C
C  OK preparse an expression
C
 500  continue
      if(ier.eq.0) then
         call uupper(zinput)
         call mmxpnd(zinput,zwkstr)     ! expand min/max
      endif
C
C  nothing implemented yet...
C
C
C  exit...
C
 1000 continue
      if(idebug) then
         ilz=max(1,len_trim(zinput))
         write(6,*) 'plcprep:  zinput:  ',zinput(1:ilz)
         ilz=max(1,len_trim(zxpres))
         write(6,*) 'plcprep:  zxpres:  ',zxpres(1:ilz)
         ilz=max(1,len_trim(zcmd2))
         write(6,*) 'plcprep:  zcmd2:  ',zcmd2(1:ilz)
         write(6,*) '...ieval = ',ieval
         write(6,*) '...icmdexec = ',icmdexec
      endif
C
      return
      end
