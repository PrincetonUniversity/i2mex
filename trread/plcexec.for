      subroutine plcexec(istat,iwrk1,iwrk2,ipt,ier)
C
C  execute a command, as parsed by a prior "rppcpars" call...
C
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod

      implicit NONE

      integer istat                     ! input:  calculator accumulator status
      integer iwrk1                     ! workspace address
      integer iwrk2                      ! address of accumulator data
      integer ipt                       ! input:  f(t) storage area
      integer ier                       ! output:  status code, 0=OK
C
      character*21 zid,zidsave
      character*10 zrunid
      character*22 zidarg
      character*32 zuns
      character*64 zlbl
C
      character*32 zstr
C
      character*128 zbuf,zpath
      character*5 ztest5
C
      character*1 typ_epst,typ_epsx,zchar
C
c
c  dmc: declare all locals, IMPLICIT NONE added 1/2011
c
      integer :: indloc,isizloc,iptloc

      integer :: lunt,lunzer,ipx,ind1,ind2,iadj,indx,iopt,ifr
      integer :: idnum,iclass,iuser,iwarn,itarg,imap,iflag
      integer :: inx,ix,it,ia,ift,ifind_ordr,ipf,indt
      integer :: ifb,imin,ier2,ier3,ilf,ilz,il,ipkg,if,itypb
      integer :: isign,istz,ifin,iscalar,iprofile,ifun
      integer :: i,j,inum,ipos,isiz1,isiz2,ila,isiznew

      real :: zintrp,zx,zconst,zxcept,zdelt,zdelx,zepst,zepsx
      real :: zdelta_t,ztestval,z0,z2,zshift
C----------------------------------
C
      ier=66
      lunt=lunzer(0)
C
C WORKSPACES
C
CXX -- values are passed --      call setup_workspaces
C
C X axis
C
      ipx=0
      if(istat.gt.0) then
         if(nlxfot(istat)) then
            call dmgxot(istat,ind1,ind2)
            ipx=locd(ind1)
            if(istat.eq.2) ipx=locd(ind2)
         endif
      endif
C
C error check
C
      if((kcmd.lt.1).or.(kcmd.gt.ncmdrpp)) then
         ier=88
         call zermsg(
     >      ' ??plcexec -- command parse (rppcpars) not done.')
         go to 1000
      endif
C
      call zermsg(' %plcexec:  command execution starting:  '//
     >   rppcmds(kcmd))
C
C----------------------------------------------------------------------
      if(rppcmds(kcmd).eq.'SAVE') then
C
C  *** SAVE COMMAND ***
C  abbrev:
         call plcgarg(1,zid)
         call idchek_setadj(zid,iadj)
         call idchek(zid,idnum,iclass,iuser,iadj)
         iwarn=0
         ier=0
         if((iclass.eq.1).and.(iuser.eq.0)) ier=1 ! can't assign:  file scalar
         if((iclass.eq.2).and.(iuser.eq.0)) ier=2 ! can't assign:  file profile
         if(iclass.eq.3) ier=3          ! can't assign:  multigraph name
         if(iclass.gt.3) ier=4
         if(ier.eq.4) then
            call zermsg(' ?plcexec: SAVE failed, bad identifier:  '//
     >         zid)
         else if(ier.ne.0) then
            call zermsg(' ?plcexec: SAVE failed, identifier in use:  '//
     >         zid)
         else
            isiznew = -1
            isizloc = -2
            if((istat.ne.0).and.(iclass.gt.0)) then
C  user defined function of same name exists:  delete it.
               call zermsg(' %plcexec:  replacing user function:  '//
     >            zid)
               if(istat.gt.0) then
                  !  profile: if size is same will reuse space
                  isiznew = nzonex(istat)*ntr
                  call dmdloc(zid,indloc,isizloc,iptloc)
               endif
               if(isiznew.ne.isizloc) then
                  !  size not the same (or not a profile)
                  call plcdelfn(zid,ier)
                  call setup_workspaces ! update might be needed
               endif
               iwarn=1
            endif
C
C  label:
            call plcgarg(2,zlbl)
C
C  units:
            call plcgarg(3,zuns)
C
            if(istat.eq.0) then
               ier=1
               call zermsg(
     >' ?plcexec:  SAVE failed, calculator accumulator data not ready.')
            else if(istat.lt.0) then
C  save scalar function
               ier=0
               call plftmk4(time,datbuf(iwrk2),ntt,zlbl,zuns,zid)
            else
C  save profile function
               ier=0
               if(iadj.eq.0) ier=-88
               if(isiznew.ne.isizloc) then
                  ! need to save function in new slot
                  call plsfsave(zid,zlbl,zuns,istat,iwrk2,ier)
               else
                  ! update old slot (isiznew=isizloc means replace at same loc)
                  labelr(idnum)=zlbl
                  unitsr(idnum)=zuns
                  do i=1,isizloc
                     datbuf(iptloc+i-1)=datbuf(iwrk2+i-1)
                  enddo
                  ier=0
               endif
            endif
            if((ier.eq.0).and.(iwarn.eq.1)) then
               call zermsg(' %plcexec:  successfully redefined:  '//zid)
            endif
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'TIME_TRACE') then
C
C  *** TIME_TRACE COMMAND ***
C  extract f(t) at fixed x index from profile
C
         ier=0
         if(istat.le.0) then
            call zermsg(
     >         ' ?plcexec:  TIME_TRACE:  no profile data found.')
            ier=1
         else
C  index argument
            inx=nzonex(istat)
            call plcgarg(1,zstr)
            call smargtr(zstr,zintrp,'I',' ',zchar,ier)
            if(zchar.eq.'I') then
               zx=1.0+float(inx)*zintrp
               ix=zx
            else
               ix=zintrp+0.5
            endif
            ix=max(1,min(inx,ix))
C  extract trace
            do it=1,ntr
               ia=iwrk2+(it-1)*inx+(ix-1)
               datbuf(iwrk1+it-1)=datbuf(ia)
            enddo
C  interpolate from profile timebase to scalar timebase
            call ttintrp(time3,datbuf(iwrk1),ntr,
     >         time,datbuf(iwrk2),ntt)
C  set status:  scalar function of time
            istat=-1
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'XINTERP') then
C
C  *** XINTERP COMMAND ***
C  interpolate f(x,t) -> f(t) at x(t)
C
         ier=0
         if(istat.le.0) then
            call zermsg(
     >         ' ?plcexec:  XINTERP:  no profile data found.')
            ier=1
         else
            call dmgfotx_ww(2,ipt,ier,iwrk1,iwrk2)
         endif
         if(ier.eq.0) then
            call plcgarg(1,zstr)
            zid=zstr(1:len(zid))
            ift=ifind_ordr(abt,iordrt,nft,zid)
            if(ift.eq.0) then
               read(zstr,'(g20.0)',iostat=ier) zconst
               if(ier.ne.0) ier=2
               IPF=0
            else
               IPF=IPT+(IFT-1)*NTT      ! address of X scalar fcn data
            endif
         endif
         if(ier.gt.1) then
            call zermsg(
     >         ' %plcexec:  XINTERP:  "X" argument decode error: '//zid)
         else
            call plcgarg(2,zstr)
            read(zstr,'(g20.0)',iostat=ier) zxcept
            if(ier.ne.0) ier=3
         endif
         if(ier.gt.2) then
            call zermsg(
     >       ' %plcexec:  XINTERP:  "EXCEPTION" argument decode error: '
     >         //zid)
         else
C  OK do the interpolation...
            call plcxntrp(iwrk1,iwrk2,istat,ipf,zconst,zxcept)
            call ttintrp(time3,datbuf(iwrk1),ntr,
     >         time,datbuf(iwrk2),ntt)
            istat=-1
         endif
C  index argument
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'DELETE') then
C
C  *** DELETE COMMAND ***
C  remove user defined function(s)
C
         call plcgarg(1,zid)
         if(zid.eq.'*') then
            ier=0
            call plcdstar               ! delete all user functions
         else
            ier=0
            call plcdelfn(zid,ier) ! delete named user function
            call setup_workspaces ! update might be needed
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'TIME_DERIV') then
C
C  *** TIME DERIVATIVE ***
C
         ier=0
         if(istat.eq.0) then
            ier=1
            call zermsg(
     >         ' ?plcexec:  TIME_DERIV failed, '//
     >         'calculator accumulator data not ready.')
         else
            call plcdfdt(istat,iwrk1,iwrk2,ipt)
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'RMJMAP') then
C
C  *** RMJMAP COMMAND:  flux zone/surface -> RMAJM (midplane major radius)
C
         ier=0
         if(istat.eq.0) then
            ier=1
            call zermsg(
     >         ' ?plcexec:  TIME_DERIV failed, '//
     >         'calculator accumulator data not ready.')
         else
            call plcmjr(istat,ier)
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'XMAP') then
C
C  *** XMAP COMMAND:  RMAJM (midlane Rmajor) -> flux zone/surface
C
         ier=0
         call plcgarg(1,zstr)
         call uupper(zstr)
         if(zstr.eq.'CTR') then
            itarg=1
         else if(zstr.eq.'BDY') then
            itarg=2
         else
            itarg=0
         endif
         if(itarg.eq.0) then
            call zermsg(' ?plcexec: XMAP failed, bad argument:  '//zstr)
            call zermsg(
     >         '  (x target) argument should be "CTR" or "BDY".')
            ier=1
         endif
C
         call plcgarg(2,zstr)
         call uupper(zstr)
         if(zstr.eq.'OUT') then
            imap=1
         else if(zstr.eq.'IN') then
            imap=2
         else if(zstr.eq.'AVG') then
            imap=3
         else
            imap=0
         endif
         if(imap.eq.0) then
            call zermsg(' ?plcexec: XMAP failed, bad argument:  '//zstr)
            call zermsg(
     >         '  (method) argument should be "IN", "OUT", or "AVG".')
            ier=2
         endif
C
         if(ier.eq.0) then
            call plcbnd(istat,itarg,imap,ier)
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'SMOOTH') then
C
C  *** SMOOTH COMMAND ***
C
         ier=0
         if(istat.eq.0) then
            call zermsg(
     >         ' ?plcexec:  SMOOTH aborted, no accumulator data.')
            ier=1
         else
C  delta(t) argument
            call plcgarg(1,zstr)
            call smargtr(zstr,zdelt,'I',' ',zchar,ier)
            if(ier.ne.0) go to 20
            if(zchar.eq.'I') then
               indt=1
            else
               indt=0
            endif
C  delta(x) argument
            call plcgarg(2,zstr)
            call smargtr(zstr,zdelx,'I',' ',zchar,ier)
            if(ier.ne.0) go to 20
            if(zchar.eq.'I') then
               indx=1
            else
               indx=0
            endif
C  eps(t) argument
            call plcgarg(3,zstr)
            call smargtr(zstr,zepst,'%RA','A',typ_epst,ier)
            zepst=abs(zepst)
            if(ier.ne.0) go to 20
C  eps(x) argument
            call plcgarg(4,zstr)
            call smargtr(zstr,zepsx,'%RA','A',typ_epsx,ier)
            zepsx=abs(zepsx)
            if(ier.ne.0) go to 20
C
            call smfunc(indt,zdelt,indx,zdelx,
     >         typ_epst,zepst,typ_epsx,zepsx,istat,iwrk2,ier)
C
 20         continue
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'TIME_AVG') then
C
C  *** TIME_AVG COMMAND ***
C
         ier=0
         if(istat.eq.0) then
            call zermsg(
     >         ' ?plcexec:  TIME_AVG aborted, no accumulator data.')
            ier=1
         else
            call plcgarg(1,zstr)
            read(zstr,'(g20.0)',err=9) zdelta_t
            go to 10
 9          continue
            call zermsg(
     >         ' %plcexec:  TIME_AVG delta_t argument invalid:  '//zstr)
            ier=2
 10         continue
            if(ier.eq.0) then
               iopt=1
               if(zdelta_t.lt.0.0) then
                  iopt=2
                  zdelta_t=-zdelta_t
               endif
               if(istat.lt.0) then
                  inx=-1
               else
                  inx=nzonex(istat)
               endif
               if(abs(zdelta_t).lt.1.0e-10) then
                  write(lunt,*) ' %plcexec:  TIME_AVG:  delta(t) = ',
     >               zdelta_t,' too small; reset.'
                  zdelta_t=1.0e-10
               endif
               if(iopt.eq.2) write(lunt,*) ' (double inverse rule) '
               call smtima(iwrk2,inx,iwrk1,zdelta_t,iopt)
            endif
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'TIME_INT') then
C
C  *** TIME INTEGRAL ***
C
         ier=0
         call plcgarg(1,zstr)           ! fetch T0
         call plc_timi(istat,zstr,ier)
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'LABEL') then
C
C  ** LABEL COMMAND ***
C
         ier=0
C
         call plcgarg(1,zid)
         call plcgarg(2,zlbl)
         call plcgarg(3,zuns)
C
         if((zlbl.eq.'%unchanged').and.(zuns.eq.'%unchanged')) then
            call zermsg(' %plcexec warning:  LABEL:  no changes given.')
         else
C
            iflag=0
            ift=ifind_ordr(abt,iordrt,nft,zid)
            if(ift.gt.0) then
               iflag=1
               if(zlbl.ne.'%unchanged') labelt(ift)=zlbl
               if(zuns.ne.'%unchanged') unitst(ift)=zuns
            else
               ifr=ifind_ordr(abr,iordrr,nfxt,zid)
               if(ifr.gt.0) then
                  iflag=2
                  if(zlbl.ne.'%unchanged') labelr(ifr)=zlbl
                  if(zuns.ne.'%unchanged') unitsr(ifr)=zuns
               else
                  ifb=ifind_ordr(abb,iordrb,nbal,zid)
                  if(ifb.gt.0) then
                     iflag=3
                     if(zlbl.ne.'%unchanged') labelb(ifb)=zlbl
                     if(zuns.ne.'%unchanged') then
                        call zermsg(' %plcexec warning:  LABEL:  '//
     >                     'multigraph units label not changed.')
                        call zermsg('  ...multigraph units label'//
     >                     ' is taken from member functions.')
                     endif
                  endif
               endif
            endif
C
            if(iflag.eq.0) then
               call zermsg(
     >            ' %plcexec warning:  LABEL:  invalid id:  '//zid)
            endif
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'MINPRO') then
         ier=0
         if(istat.le.0) then
            call zermsg(
     >         ' ?plcexec:  MINPRO:  no profile data found.')
            ier=1
         else
            call plcxtrac(istat,iwrk1,iwrk2,ipx,1,0)
            call ttintrp(time3,datbuf(iwrk1),ntr,
     >         time,datbuf(iwrk2),ntt)
            istat=-1                    ! mark as scalar function
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'MAXPRO') then
         ier=0
         if(istat.le.0) then
            call zermsg(
     >         ' ?plcexec:  MAXPRO:  no profile data found.')
            ier=1
         else
            call plcxtrac(istat,iwrk1,iwrk2,ipx,0,0)
            call ttintrp(time3,datbuf(iwrk1),ntr,
     >         time,datbuf(iwrk2),ntt)
            istat=-1                    ! mark as scalar function
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'XLOCATE') then
         ier=0
         if(istat.le.0) then
            call zermsg(
     >         ' ?plcexec:  MINPRO:  no profile data found.')
            ier=1
         else
            call plcgarg(1,zstr)
            if(zstr.eq.'MIN') then
               imin=1
            else if(zstr.eq.'MAX') then
               imin=0
            else
               imin=-1
               read(zstr,'(g20.0)',iostat=ier) ztestval
               if(ier.ne.0) then
                  call zermsg(
     >               ' ?plcexec:  TEST_VALUE decode failure:  '//zstr)
               else
                  call plcgarg(2,zstr)
                  read(zstr,'(g20.0)',iostat=ier2) z2
                  if(ier2.ne.0) then
                     call zermsg(
     >                  ' ?plcexec:  NOT_UNIQUE fallback value decode'//
     >                  ' failed:  '//zstr)
                  endif
                  call plcgarg(3,zstr)
                  read(zstr,'(g20.0)',iostat=ier3) z0
                  if(ier3.ne.0) then
                     call zermsg(
     >                  ' ?plcexec:  NOT_FOUND fallback value decode'//
     >                  ' failed:  '//zstr)
                  endif
                  ier=max(ier,ier2,ier3)
               endif
            endif
            if(ier.eq.0) then
               if(imin.ge.0) then
                  call plcxtrac(istat,iwrk1,iwrk2,ipx,imin,1)
               else
                  call plcxloc(istat,iwrk1,iwrk2,ipx,ztestval,z2,z0)
               endif
               call ttintrp(time3,datbuf(iwrk1),ntr,
     >            time,datbuf(iwrk2),ntt)
               istat=-1                 ! mark as scalar function
            endif
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'RFETCH') then
         ier=0
C
C  path
         call plcgarg(1,zbuf)
         ztest5=zbuf(1:5)
         call uupper(ztest5)
         if(ztest5.eq.'MDS+:') then
            zpath=zbuf
         else if((zbuf.eq.'cwd').or.(zbuf.eq.' ')) then
            zpath=' '
         else if(zbuf.eq.'.') then
C
C  use current run's path if argument = "."
C
            if(nlmds) then
               ilf=len_trim(fdisk)
               ila=index(fdisk,'@')
               zpath='MDS+:'//fdisk(6:ila-1)//':'//fdisk(ila+1:ilf)
               if(lfdir.gt.0) then
                  ilz=len_trim(zpath)
                  zpath(ilz+1:)=':'//fdir(1:lfdir)
               endif
            else
               zpath=' '
               if(lfdir.gt.0) then
                  zpath=fdir(1:lfdir)
               endif
            endif
         else if((zbuf(1:1).eq.'.').and.(zbuf(2:2).ne.'.').and.
     >         (zbuf(2:2).ne.'/')) then
C
C  use RPLOT convention; input of...
C    .tok.yy
C  maps to
C    $RESULDIR/tok.yy
C
            ilz=len_trim(zbuf)
            call ufilnam('RESULTDIR',zbuf(2:ilz),zpath)
         else
            call ufilnam(zbuf,'@#$%^&',zpath)
            il=index(zpath,'@#$%^&')
            zpath(il:)=' '
         endif
C
C  runid
         call plcgarg(2,zrunid)
C
C  function id within run "zrunid"
         call plcgarg(3,zid)
C
C  function id to use in current session
         call plcgarg(4,zidsave)
         call idchek(zidsave,idnum,iclass,iuser,1)
         ier=0
         if((iclass.eq.1).and.(iuser.eq.0)) ier=1 ! can't assign:  file scalar
         if((iclass.eq.2).and.(iuser.eq.0)) ier=2 ! can't assign:  file profile
         if(iclass.eq.3) ier=3          ! can't assign:  multigraph name
         if(iclass.gt.3) ier=4
         if(ier.eq.4) then
            call zermsg(
     >         ' ?plcexec:  illegal character in identifier:  '//
     >         zidsave)
         else if(ier.ne.0) then
            call zermsg(
     >         ' ?plcexec:  RFETCH failed, identifier in use:  '//
     >         zidsave)
            ier=1
         else
C  remove user defined function to make room...
            if(iclass.gt.0) then
               call zermsg(' %plcexec:  replacing user function:  '//
     >         zidsave)
               call plcdelfn(zidsave,ier)
               call setup_workspaces  ! update might be needed
            endif
         endif
C
C  OK go fetch the data
C
         if(ier.eq.0) then
            zshift=0.0
            call plcfget(zpath,zrunid,zid,zidsave,zshift,istat,ier)
         endif
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'GS2FETCH') then
         ier=0
         call plcgarg(1,zpath)
         call plcgarg(2,zid)
         call plcgarg(3,zidsave)
         call trgs2fetch(zpath,zid,zidsave,ier)
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'MG_CREATE') then
         ier=0
         call plcgarg(1,zid)
         ipkg=ifind_ordr(abb,iordrb,nbal,zid)
         if(ipkg.ne.0) then
            call zermsg(
     >         ' ?plcexec:  MG_CREATE:  package already exits:  '//
     >         zid)
            ier=1
         else
C  create package
            zidsave=zid
            call plcgarg(2,zlbl)        ! label
            call plcgarg(3,zidarg)         ! 1st fcn (required)
            isign=1
            istz=1
            ifin=len_trim(zidarg)
C  check sign
            if(zidarg(1:1).eq.'+') istz=istz+1
            if(zidarg(1:1).eq.'-') then
               istz=istz+1
               isign=-1
            endif
            zid=zidarg(istz:ifin)
C  scalar or profile
            iscalar=ifind_ordr(abt,iordrt,nft,zid)
            if(iscalar.eq.0) then
               iprofile=ifind_ordr(abr,iordrr,nfxt,zid)
               if(iprofile.eq.0) then
                  ier=1
                  call zermsg(' ?plcexec:  unrecognized function id: '//
     >               zid)
               else
C  create profile multigraph
                  ifun=iprofile
                  nbal=nbal+1
                  iintb(nbal)=0
                  unitsb(nbal)=unitsr(ifun)
               endif
            else
C  create scalar multigraph
               ifun=iscalar
               nbal=nbal+1
               iintb(nbal)=1
               unitsb(nbal)=unitst(ifun)
            endif
            if(ier.eq.0) then
C  labeling
               abb(nbal)=zidsave
               labelb(nbal)=zlbl
C  insert first function
               infb(nbal)=1
               ifunb(infb(nbal),nbal)=isign*ifun
C  loop to add additional fcns
               do if=4,ncmdargs(kcmd)
                  call mgtryadd(nbal,if)
               enddo
               ila=len_trim(abb(nbal))
               write(lunt,8801)
     >            rppcmds(kcmd),abb(nbal)(1:ila),infb(nbal)
C  maintain index ordering
               call aordr_add(abb,iordrb,nbal)
            endif
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'MG_DELETE') then
         ier=0
         call plcgarg(1,zid)
         ipkg=ifind_ordr(abb,iordrb,nbal,zid)
         if(ipkg.eq.0) then
            call zermsg(' %plcexec:  MG_DELETE:  no such package:  '//
     >         zid)
         else
            call mgdelpkg(ipkg)
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'MG_ADDFUN') then
         ier=0
         call plcgarg(1,zid)
         ipkg=ifind_ordr(abb,iordrb,nbal,zid)
         if(ipkg.eq.0) then
            call zermsg(' ?plcexec:  MG_ADDFUN:  no such package:  '//
     >         zid)
            ier=1
         else
            ila=len_trim(abb(ipkg))
            do if=2,ncmdargs(kcmd)
               call mgtryadd(ipkg,if)
            enddo
            write(lunt,8801) rppcmds(kcmd),abb(ipkg)(1:ila),infb(ipkg)
         endif
C
C----------------------------------------------------------------------
      else if(rppcmds(kcmd).eq.'MG_DELFUN') then
         ier=0
         call plcgarg(1,zid)
         ipkg=ifind_ordr(abb,iordrb,nbal,zid)
         if(ipkg.eq.0) then
            call zermsg(' ?plcexec:  MG_DELFUN:  no such package:  '//
     >         zid)
            ier=1
         else
            ila=len_trim(abb(ipkg))
            itypb=iintb(ipkg)           !0:profile, 1:scalar
            do i=2,7
               call plcgarg(i,zid)
               ifun=0
               if(zid.ne.'%empty') then
                  if(itypb.eq.0) then
                     ifun=ifind_ordr(abr,iordrr,nfxt,zid)
                     if(ifun.eq.0) then
                        call zermsg(' %plcexec:  MG_DELFUN:  not '//
     >                     'a profile function id:  '//zid)
                     endif
                  else
                     ifun=ifind_ordr(abt,iordrt,nft,zid)
                     call zermsg(' %plcexec:  MG_DELFUN:  not '//
     >                  'a scalar function id:  '//zid)
                  endif
               endif
               if(ifun.gt.0) then
                  ipos=0
                  inum=infb(ipkg)
                  do j=1,inum
                     if(abs(ifunb(j,ipkg)).eq.ifun) then
                        ipos=j
                     endif
                  enddo
                  if(ipos.gt.0) then
                     call mgdelfun(ipkg,ifun,ipos)
                  else
                     call zermsg(' %plcexec:  MG_DELFUN:  not a '//
     >                  'member of '//abb(ipkg)(1:ila)//':  '//zid)
                  endif
               endif
            enddo
            write(lunt,8801) rppcmds(kcmd),abb(ipkg)(1:ila),infb(ipkg)
         endif
C
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
      endif                             ! command string test
C
C----------------------------------------------------------------------
 8801 format(
     >   ' %plcexec:  ',a,':  package ',a,' now contains ',
     >   i2,' member functions.')
C----------------------------------------------------------------------
C
      if(ier.eq.66) then
         call zermsg(
     >      ' ??plcexec code error, no branch for kcmd:  '//
     >      rppcmds(kcmd))
         call abortt
      endif
C
C----------------------------------------------------------------------
C
      if(ier.eq.77) then
         call zermsg(' ??plcexec command not yet implemented:  '//
     >      rppcmds(kcmd))
         go to 1000
      endif
C
C----------------------------------------------------------------------
C
 1000 continue
      return

      CONTAINS

        subroutine setup_workspaces
          ! contained routine...

          CALL DMDLOC('%WRK1',IND1,ISIZ1,IWRK1)
          CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)

        end subroutine setup_workspaces

      end
