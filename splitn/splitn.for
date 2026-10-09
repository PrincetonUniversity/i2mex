C******************** START FILE SPLITN.FOR ; GROUP SPLITN *************
C
      INTEGER FUNCTION SPLITN(zpath,zfile,ICRUN)
C
      use splitn_module
C
      implicit NONE
C
C  READS A NAMELIST DATA FILE AND SPLITS THE VARIABLES
C
C  modified WIE 09 Apr 1997 -- add in EFITIN namelist --
C    pass it thru to <runid>TR.EFITIN .. cf ADDNL
C
C  modified DMC 16 Nov 1994 -- directory path argument zpath added
C    this argument specifies the source directory for the tr.dat
C    file; if blank, the current working directory is assumed.
C    if non-blank, the argument should be left justified in the
C    passed character string with only trailing blanks allowed.
C    Full punctuation appropriate to the current OS should be
C    supplied (e.g. include trailing "/" if this is unix).
C
C  ALSO DROPS COMMENTS
C  FIRST DRAFT BY M. THOMPSON, EXPANDED BY D. MC CUNE
C
C  redone dmc Dec 2002 -- partitioning of namelist driven from
C  "nxlist.summary" file...
C
C
C  ***return value***  SPLITN=1 -- successful!
C                             anything else: error
C
C--------------------------------------------------------------------
C
      character*(*) zpath  ! optional path (directory) for namelist
      character*(*) zfile  ! optional filename for namelist
      CHARACTER*(*) ICRUN  ! namelist run-id
C
C--------------------------------------------------------------------
C
      CHARACTER FNAMI*128  ! input filename (generic TRANSP namelist)
      CHARACTER FNAMO*20   ! output filename (fortran readable namelist)
C
C---------------------------------------------------
      integer i,inb,inc,imds,ios
C
C  some for debugging...
C
      integer j,k,il,inum,ia,ilinp,ilinq,iadr
C
C---------------------------------------------------
C
      imds = 1
C
      inc=index(ICRUN,' ')-1
      if(inc.lt.0) inc=len(ICRUN)
C
      if(zpath.eq.' ') then
         if(zfile.eq.' ') then
C try MDSplus extracted file first
            imds = 2
            FNAMI=ICRUN(1:INC)//'TR.DAT_MDS'
         else
            FNAMI=zfile
         endif
      else
         inb=index(zpath,' ')-1
         if(inb.lt.0) inb=len(zpath)
         if(zfile.eq.' ') then
            FNAMI=zpath(1:inb)//ICRUN(1:INC)//'TR.DAT'
         else
            FNAMI=zpath(1:inb)//zfile
         endif
      endif
C
      FNAMO=ICRUN(1:INC)//'TR.ZDA'
C
      i=1
      SPLITN=0
      do while (i .le. imds .and. splitn .eq. 0)
C
C  SPLITN namelist read: quiet operation:
         call read_nl(fnami,ios, quiet=.TRUE.)
C
         IF(IOS.EQ.0) THEN
            SPLITN=1
         else if (imds .eq. 2 .and. i .eq. 1) then
C try regular namelist
            FNAMI=ICRUN(1:INC)//'TR.DAT'
         else
            RETURN              ! TR.DAT NAMELIST FILE NOT OPENED.
         ENDIF
         i=i+1
      end do
C
C  make sure no old TR.ZDA files are left around.
C
      open(unit=lun,file=fnamo,status='old',iostat=ios)
      if(ios.eq.0) then
                                        ! delete the file
         close(unit=lun)
         call fdelete(fnamo,ios)
      endif
C
C  write fortran readable namelists (TR.ZDA file).
C
      if(splitn.eq.1) then
         call splitn_f77(fnamo,ios)
         if(ios.ne.0) splitn=0
      endif
C
#ifdef __DEBUG
      do i=1,nvars
         j=var_order(i)
         write(6,*) '--------------------------------'
         write(6,*) ' name: ',varlist(j)%name
         write(6,*) ' type: ',varlist(j)%type,
     >      ' ... chsize = ',varlist(j)%chsize
         write(6,*) ' belonging to: ',varlist(j)%naml
         write(6,*) ' rank: ',varlist(j)%rank
         inum=1
         do k=1,varlist(j)%rank
            write(6,*) '    dim. ',k,':  (',varlist(j)%dims(1,k),':',
     >         varlist(j)%dims(2,k),')'
            inum=inum*(varlist(j)%dims(2,k)-varlist(j)%dims(1,k)+1)
         enddo
         write(6,*) '    size: ',inum
         if(varlist(j)%short_dflt.ne.' ') then
            write(6,*) ' short default string: ',
     >         varlist(j)%short_dflt
         endif
         k=varlist(j)%long_dflt_addr
         if(k.gt.0) then
            il=len_trim(long_dflts(k))
            write(6,*) ' long default string: ',long_dflts(k)(1:il)
         endif
         write(6,*) ' ======='
         ia=varlist(j)%nlinadr
         ilinp=0
         do k=ia,ia+inum-1
            if(ilines(k).ne.0) then
               if(ilines(k).ne.ilinp) then
                  ilinp=ilines(k)
                  ilinq=ordl(ilinp)
                  il=len_trim(textnl(ilinq))
                  write(6,*) '   file: ',textnl(ilinq)(1:il)
               endif
            endif
         enddo
         write(6,*) ' ======='
         write(6,*) ' values:'
         do k=1,inum
            iadr=varlist(j)%addr + k - 1
            if(varlist(j)%type(1:1).eq.'C') then
               il=len_trim(chbuf(iadr))
               write(6,*) '    ',chbuf(iadr)(1:il)
            else if(varlist(j)%type.eq.'R') then
               write(6,*) '    ',rbuf(iadr)
            else if(varlist(j)%type.eq.'I') then
               write(6,*) '    ',intbuf(iadr)
            else if(varlist(j)%type.eq.'L') then
               write(6,*) '    ',logbuf(iadr)
            else if(varlist(j)%type.eq.'D') then
               write(6,*) '    ',dbuf(iadr)
            endif
         enddo
      enddo
#endif
C     
      return
      end
