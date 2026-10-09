      subroutine rppcpars(zpre,cmdstr,expr,descr,ieval,icmdexec,ier)
      use rpcalc_mod
C
C  parse an RPLOT calculator "auxilliary command" string
C
      character*(*) zpre                ! command guard string (input)
C
C  zpre is prepended to command to build descr string output
C
      character*(*) cmdstr              ! command string (input)
C
C  cmdstr starts with actual command (not guard prefix)
C
      character*(*) expr                ! expression to evaluate (output)
      character*(*) descr               ! expression description
C
C  either expr -- if command does not cause a computation which alters
C  the calculator accumulator, or zpre//cmdstr if the command does cause a
C  computation.
C
      logical ieval                     ! output .TRUE. for expr evaluation
      logical icmdexec                  ! output .TRUE. unless error occurs
      integer ier                       ! status code (output, 0=OK)
C
C  return ier=1 if the command is not recognized; =2 if there is an
C    error with the command arguments.
C
C  error messages will be written on lunrpc which defaults to 6 (stdout)
C  but can be modified by a call to plc_msgs.
C
      character*10 ztest
C
      character*1 zquot
C
      integer iargsta(100),iargend(100),iargeqs(100),inarg
C
      logical idebug
      data idebug/.false./
C
C----------------------
C
C  1.  find command verb -- stuff before the first "(".
C
      ieval=.false.
      icmdexec=.false.
      ier=0
      ilen=len_trim(cmdstr)           ! non-blank length of cmdstr
C
      cmdbuf=cmdstr(1:ilen)
      lcmdbuf=ilen
C
      do ic=1,ilen
         if(cmdstr(ic:ic).eq.'(') go to 10
      enddo
C
      ier=1
      write(lunrpc,9001) cmdstr(1:ilen)
 9001 format(' ?rppcpars -- command:  opening parentheses not found.'/
     >   1x,a)
      return
C
C----------------------
C  command found -- is it valid?
C
 10   continue
      iargst=ic
      icmd=ic-1
      call uupper(cmdbuf(1:icmd))
      do i=1,ncmdrpp
         if(cmdbuf(1:icmd).eq.rppcmds(i)) go to 20
      enddo
C
      ier=1
      write(lunrpc,9002) cmdbuf(1:icmd),cmdbuf(1:lcmdbuf)
 9002 format(' ?rppcpars -- unrecognized command:  ',a/
     >   '  command was:  ',a)
      return
C
C-------------------------------------
C  valid command.  Now scan arguments.
C
 20   continue
      kcmd=i
      klencmd=icmd
C
      keyfmt=.false.
C
C  parentheses level
C
      ilevel=1
C
C  scan argument region
      inarg=1
      iargsta(inarg)=0
      iargend(inarg)=0
      iargeqs(inarg)=0
C
      ic=iargst
      zquot=' '
 30   continue
      ic=ic+1
      if (ic.gt.ilen) then
         write(lunrpc,9003) cmdbuf(1:klencmd)
 9003    format(' ?rppcpars -- closing parentheses not found:'/1x,a)
         ier=2
         return
      endif
C
C  check next character for parentheses, etc.
C
      if(zquot.eq.' ') then
C
C  not inside a quote string
C
C  skip blank or tab
C
         if(cmdbuf(ic:ic).eq.' ') go to 30
         if(cmdbuf(ic:ic).eq.char(9)) go to 30
C
C  check for start of new argument, or end of command
C
         if(ilevel.eq.1) then
            if(cmdbuf(ic:ic).eq.',') then
               inarg=inarg+1
               iargsta(inarg)=0
               iargend(inarg)=0
               iargeqs(inarg)=0
               go to 30
            endif
            if(cmdbuf(ic:ic).eq.')') then
               if(ic.lt.ilen) then
                  write(lunrpc,9005) cmdbuf(ic+1:ilen),cmdbuf(1:lcmdbuf)
 9005             format(' %rppcpars -- extra characters ignored:  ',a/
     >               '  command was:  ',a)
                  lcmdbuf=ic
               endif
               go to 100                ! done scanning
            endif
         endif
C
C  ok we have nonblank that is somehow part of the argument
C
         if(iargsta(inarg).eq.0) iargsta(inarg)=ic
         iargend(inarg)=ic
C
C  check for nested parentheses
C
         if(cmdbuf(ic:ic).eq.'(') then
            ilevel=ilevel+1
            go to 30                    ! keep scanning
         endif
         if(cmdbuf(ic:ic).eq.')') then
            ilevel=ilevel-1
            go to 30                    ! keep scanning
         endif
C
C  check for quote string start
C
         if((cmdbuf(ic:ic).eq.'"').or.(cmdbuf(ic:ic).eq.'''')) then
            zquot=cmdbuf(ic:ic)
            go to 30                    ! scan for end of quote string
         endif
C
C  check for equals sign
C
         if(ilevel.eq.1) then
            if(cmdbuf(ic:ic).eq.'=') then
               iargeqs(inarg)=ic
               keyfmt=.true.
               go to 30
            endif
         endif
C
C  keep scanning
C
         go to 30
C
      else
C  looking for end of quote string
         iargend(inarg)=ic
         if(cmdbuf(ic:ic).eq.zquot) then
            zquot=' '
            go to 30                    ! resume scan
         endif
      endif
      go to 30
C
C-----------------------------------------------------------------
C  end of command has been found, and initial information on the
C  location of arguments in the command line.
C
 100  continue
C
      do i=1,nargrpp
         kposarg(i)=0
         klenarg(i)=0
      enddo
C
      if(inarg.gt.ncmdargs(kcmd)) then
         write(lunrpc,9019) ncmdargs(kcmd),cmdbuf(1:lcmdbuf)
 9019    format(' ?rppcpars:  too many arguments, ',i2,' expected.'/
     >      '  command was:  ',a)
         ier=2
         return
      endif
C
      if(.not.keyfmt) then
C
C  positional argument syntax
C
         do i=1,inarg
            if(iargsta(i).gt.0) then
               kposarg(i)=iargsta(i)
               klenarg(i)=iargend(i)-iargsta(i)+1
            endif
         enddo
C
         ier=0
C
      else
C
C  keyword argument syntax
C
         do i=1,inarg
            if(iargsta(i).gt.0) then
               if(iargeqs(i).eq.0) then
                  ier=2
                  write(lunrpc,9009) cmdbuf(1:lcmdbuf)
 9009             format(
     >' ?rppcpars:  cannot mix keyword and positional argument syntax:'/
     >1x,a)
                  return
               endif
               if(iargsta(i).eq.iargeqs(i)) then
                  ier=2
                  write(lunrpc,9010) cmdbuf(1:lcmdbuf)
 9010             format(
     >' ?rppcpars:  missing keyword before "=" sign:'/1x,a)
                  return
               endif
C
C  find the end of the keyword
C
               do ik=iargeqs(i)-1,iargsta(i),-1
                  if((cmdbuf(ik:ik).ne.' ').and.
     >               (cmdbuf(ik:ik).ne.char(9))) go to 120
               enddo
               ik=iargsta(i)
 120           continue
               ztest=cmdbuf(iargsta(i):ik)
               call uupper(ztest)
               imatch=0
               do j=1,ncmdargs(kcmd)
                  if(ztest.eq.rppkeys(j,kcmd)) imatch=j
               enddo
               if(imatch.eq.0) then
                  write(lunrpc,9012) cmdbuf(iargsta(i):ik),
     >               cmdbuf(1:lcmdbuf)
 9012             format(' ?rppcpars:  unrecognized keyword:  ',a/
     >               '  command was:  ',a)
                  ier=2
                  return
               endif
C
               jarg=imatch
               if(kposarg(jarg).gt.0) then
                  write(lunrpc,9013) ztest,cmdbuf(1:lcmdbuf)
 9013             format(
     >               ' ?rppcpars:  keyword appears more than once:  ',a/
     >               '  command was:  ',a)
                  ier=2
                  return
               endif
C
C  OK find start of the argument value
C
               if(iargend(i).gt.iargeqs(i)) then
                  do ik=iargeqs(i)+1,iargend(i)
                     if((cmdbuf(ik:ik).ne.' ').and.
     >                  (cmdbuf(ik:ik).ne.char(9))) go to 140
                  enddo
                  ik=iargend(i)
 140              continue
                  kposarg(jarg)=ik
                  klenarg(jarg)=iargend(i)-ik+1
               endif
            endif                       ! iargsta(i).gt.0 test
         enddo                          ! i = 1 to inarg loop
C
         ier = 0
C
      endif                             ! keyword format test
C
C  strip quotes off arguments, if necessary...
C
      do j=1,ncmdargs(kcmd)
         ia1=kposarg(j)
         if(ia1.gt.0) then
            ia2=ia1+klenarg(j)-1
C  enclosing quotes test
            if((cmdbuf(ia1:ia1).eq.'"').or.
     >         (cmdbuf(ia1:ia1).eq.'''')) then
               if(cmdbuf(ia2:ia2).eq.cmdbuf(ia1:ia1)) then
                  i1=ia1+1
                  i2=ia2-1
                  if(i2.lt.i1) then
                     kposarg(j)=0
                     klenarg(j)=0
                  else
                     kposarg(j)=i1
                     klenarg(j)=i2-i1+1
                  endif
C
C  end enclosing quotes conditional code
C
               endif
            endif
         endif
C
      enddo                             ! arg loop
C
C  check if any arguments have been defaulted for which user input
C  is required.
C
      do j=1,ncmdargs(kcmd)
         if(kposarg(j).eq.0) then
            if(rppadfs(j,kcmd).eq.' ') then
               ier=2
               write(lunrpc,9020) j,rppkeys(j,kcmd),cmdbuf(1:lcmdbuf)
 9020          format(
     >' ?rppcpars -- argument value required for argument #',i2,
     >' keyword = ',a/'  command was:  ',a)
            endif
         endif
      enddo
      if(ier.eq.2) return
C
      if(ier.eq.0) then
C
C  fetch expression argument now...
C
         do j=1,ncmdargs(kcmd)
            if(rppkeys(j,kcmd).eq.'EXPR') then
               if(kposarg(j).ne.0) then
                  ik1=kposarg(j)
                  ik2=ik1+klenarg(j)-1
                  if(lrbefor(kcmd)) then
                     ilop=len_trim(beforop(kcmd))
                     ieval=.true.
                     expr=
     >                  beforop(kcmd)(1:ilop)//'('//cmdbuf(ik1:ik2)//')'
                  else if(cmdbuf(ik1:ik2).ne.'$') then
                     ieval=.true.
                     expr=cmdbuf(ik1:ik2)
                  endif
               endif
            endif
         enddo
C
         icmdexec=.true.
         if(lrafter(kcmd)) then
            descr=zpre//cmdbuf(1:ilen)
         else
            if(ieval) descr=expr
         endif
C
      endif                             ! ier=0
C
      if(idebug) then
         write(6,7700) cmdbuf(1:lcmdbuf)
 7700    format(' command buffer:  ',a)
         write(6,7701) kcmd,rppcmds(kcmd),cmdbuf(1:klencmd)
 7701    format(' rpp debug:  cmd# = ',i2,' cmd = ',a,1x,a)
         if(keyfmt) then
            write(6,*) ' ... argument syntax:  keyword'
         else
            write(6,*) ' ... argument syntax:  positional'
         endif
         do j=1,ncmdargs(kcmd)
            if(kposarg(j).eq.0) then
               write(6,7703) rppkeys(j,kcmd),rppadfs(j,kcmd)
 7703          format(' arg. ',a,' defaulted to:  ',a)
            else
               ia1=kposarg(j)
               ia2=kposarg(j)+klenarg(j)-1
               write(6,7704) rppkeys(j,kcmd),cmdbuf(ia1:ia2)
 7704          format(' arg. ',a,' = ',a)
            endif
         enddo
      endif
C
      return
      end
