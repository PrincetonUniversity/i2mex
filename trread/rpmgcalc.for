      subroutine rpmgcalc(mgnames,mgsigns,nnames,prog,nprog,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C
C  given a list of functions (& signs), execute "prog(1...nprog)"
C  on each function in succession, returning the results in
C  zmgbuf(*,*):  zmgbuf(i,j) is the i'th word returned for the j'th
C  function in the list of functions
C
C  "prog" contains RPLOT calculator expressions and commands + the
C  special symbol "@"; at each occurrence of "@" the current function
C  name is substituted into the calculator command or expression.
C
C  If there is no "prog" (i.e. nprog=0 or all lines blank), just
C  read and return the data without any processing
C
C  input:
      character*(*) mgnames(*)          ! list of function names
      integer mgsigns(*)                ! associated signs (+/- 1)
      integer nnames                    ! number of names in list
C
      character*(*) prog(*)             ! calculator program
      integer nprog                     ! number of steps in program
C
      integer ibufsize                  ! maximum size of each result
C
C  output:
C
      real zmgbuf(ibufsize,*)           ! buffer array -- space for each result
C
C  zmgbuf(1...ret,k) is results of "prog..." applied to mgsigns(k)*mgnames(k)
C
      integer iret                      ! number of values returned per result
      integer istype                    ! results data -- subtype code
C               -1 = scalar f(t)
C               +N = profile f(x,t) of type N; N identifies which x axis.
      integer iwarn                     ! flag if arithmetic errors were caught
C               =0 means no arithmetic errors
      integer ier                       ! completion code, 0 = normal
C
C  iret,istype,iwarn, and ier are related to the same named quantities
C  in the "rpcalc" subroutine.  The values iwarn and ier are the maximum
C  occurring in rpcalc calls in the loop over functions.  iret and istype
C  must have the final values after executing the final calculator command
C  prog(nprog) -- or ier is set.
C
C  ** for short programs one can use
C       rpmg0cal  --  no program -- just read the data
C       rpmg1cal  --  one line program
C       rpmg2cal  --  two line program
C       rpmg3cal  --  three line program
C     which may be convenient for those wishing to avoid setting
C     up the prog(...) array.
C
C----------------------------------------------------------------------
C
      logical ireadonly
C
      character*64 zlbl
      character*32 zuns,zfxpr
      character*512 stmt
C
C
C------------------------------------
C
      iret=0
      istype=0
      iwarn=0
      ier=0
C
      if(nnames.eq.0) then
         call zermsg('%rpmgcalc:  empty multigraph list.')
         return
      endif
C
      ireadonly=.true.
      if(nprog.gt.0) then
         do i=1,nprog
            if(prog(i).ne.' ') ireadonly=.false.
         enddo
      endif
C
C  check function names & labels
C
      lunz=lunzer(0)
      iflag=0
      do i=1,nnames
         call rplabel(mgnames(i),zlbl,zuns,imulti,istyp2)
         if((imulti.eq.1).or.(zlbl.eq.'error')) then
            call zermsg(
     >         ' ?rpmgcalc:  invalid multigraph member name:  '//
     >         mgnames(i))
            ier=ier+1
         else
            if(istype.eq.0) then
               istype=istyp2
            else
               if(istype.ne.istyp2) then
                  iflag=iflag+1         ! type inconsistency
               endif
            endif
         endif
      enddo
      if(iflag.gt.0) then
         call zermsg(
     >      ' ?rpmgcalc:  multigraph member type inconsistencies.')
         ier=ier+1
         write(lunz,'(''  member name     type code'')')
         do i=1,nnames
            call rplabel(mgnames(i),zlbl,zuns,imulti,istyp2)
            write(lunz,'(3x,a,10x,i3)') mgnames(i),istyp2
         enddo
      endif
      if(ier.gt.0) return
C
      if(ireadonly) then
C
C  just read the data
C
         do i=1,nnames
            call rprofile(mgnames(i),zmgbuf(1,i),ibufsize,iret2,ier2)
            if(ier2.eq.0) then
               if(mgsigns(i).eq.-1) then
C  apply sign factor
                  do j = 1,iret2
                     zmgbuf(j,i)=-zmgbuf(j,i)
                  enddo
               endif
            else
C  note error
               ier=max(ier,ier2)
            endif
            if(iret.eq.0) then
               iret=iret2
            else if(iret.ne.iret2) then
C  check that all returned data lenghts are equal
               call zermsg(
     >            ' ??rpmgcalc:  member data length inconsistency.')
               ier=ier+1
            endif
         enddo
         go to 1000
      else
C
C  execute program on each named entity
C
         iret=0
         istype=0
         do i=1,nnames
            if(mgsigns(i).eq.-1) then
C  form negation subexpression
               iln=len_trim(mgnames(i))
               zfxpr='(-1*'//mgnames(i)(1:iln)//')'
            else
C  sign = +1:  just copy the name
               zfxpr=mgnames(i)
            endif
            ilf=len_trim(zfxpr)
C
C  process each program step for this name...
C
            do j=1,nprog
               if(prog(j).ne.' ') then
C  check for "@"-name substitution
                  ilp=len_trim(prog(j))
                  ils=0
                  icp=1
                  do ic=1,ilp
                     if(prog(j)(ic:ic).eq.'@') then
                        if(ic.gt.icp) then
                           stmt(ils+1:)=prog(j)(icp:ic-1) ! before name
                           icp=ic+1
                           ils=len_trim(stmt)
                        endif
                        stmt(ils+1:)=zfxpr(1:ilf) ! the name (expression)
                        ils=len_trim(stmt)
                     endif
                  enddo
                  if(icp.le.ilp) stmt(ils+1:)=prog(j)(icp:ilp) ! after name
C
C  compute...
C
                  call rpcalc(stmt,zmgbuf(1,i),ibufsize,
     >               iret2,istyp2,iwarn2,ier2)
C
                  ier=max(ier,ier2)
                  iwarn=max(iwarn,iwarn2)
                  if(ier.gt.0) go to 1000
               endif
            enddo
C  all stmts executed, no error so far
            if(iret.eq.0) then
               iret=iret2
            else if(iret.ne.iret2) then
C  check that all returned data lenghts are equal
               call zermsg(
     >            ' ??rpmgcalc:  result data length inconsistency.')
               ier=ier+1
            endif
            if(istype.eq.0) then
               istype=istyp2
            else if(istype.ne.istyp2) then
C  check result data type
               call zermsg(
     >            ' ??rpmgcalc:  result data type inconsistency.')
               ier=ier+1
            endif
            if(ier.gt.0) go to 1000
         enddo
      endif                             ! readonly test
C-------------------------
 1000 continue
      return
      end
C----------------------
C  call rpmgcalc with no program (just fetch the raw profiles)
C
      subroutine rpmg0cal(mgnames,mgsigns,nnames,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C  input:
      character*(*) mgnames(*)          ! list of function names
      integer mgsigns(*)                ! associated signs (+/- 1)
      integer nnames                    ! number of names in list
C
      integer ibufsize                  ! maximum size of each result
C
C  output:
C
      real zmgbuf(ibufsize,*)           ! buffer array -- space for each result
C
C  zmgbuf(1...ret,k) is results of "prog..." applied to mgsigns(k)*mgnames(k)
C
      integer iret                      ! number of values returned per result
      integer istype                    ! results data -- subtype code
C               -1 = scalar f(t)
C               +N = profile f(x,t) of type N; N identifies which x axis.
      integer iwarn                     ! flag if arithmetic errors were caught
C               =0 means no arithmetic errors
      integer ier                       ! completion code, 0 = normal
      character*1  dummy(1)
C
      call rpmgcalc(mgnames,mgsigns,nnames,dummy,0,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C
      return
      end
C----------------------
C  call rpmgcalc with a one line program...
C
      subroutine rpmg1cal(mgnames,mgsigns,nnames,stmt,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C  input:
      character*(*) mgnames(*)          ! list of function names
      integer mgsigns(*)                ! associated signs (+/- 1)
      integer nnames                    ! number of names in list
C
      character*(*) stmt                ! one line program
C
      integer ibufsize                  ! maximum size of each result
C
C  output:
C
      real zmgbuf(ibufsize,*)           ! buffer array -- space for each result
C
C  zmgbuf(1...ret,k) is results of "prog..." applied to mgsigns(k)*mgnames(k)
C
      integer iret                      ! number of values returned per result
      integer istype                    ! results data -- subtype code
C               -1 = scalar f(t)
C               +N = profile f(x,t) of type N; N identifies which x axis.
      integer iwarn                     ! flag if arithmetic errors were caught
C               =0 means no arithmetic errors
      integer ier                       ! completion code, 0 = normal
C
C
      character*512 prog(1)
C
      prog(1)=stmt
      call rpmgcalc(mgnames,mgsigns,nnames,prog,1,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C
      return
      end
C----------------------
C  call rpmgcalc with a two line program...
C
      subroutine rpmg2cal(mgnames,mgsigns,nnames,stmt1,stmt2,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C  input:
      character*(*) mgnames(*)          ! list of function names
      integer mgsigns(*)                ! associated signs (+/- 1)
      integer nnames                    ! number of names in list
C
      character*(*) stmt1,stmt2         ! two line program
C
      integer ibufsize                  ! maximum size of each result
C
C  output:
C
      real zmgbuf(ibufsize,*)           ! buffer array -- space for each result
C
C  zmgbuf(1...ret,k) is results of "prog..." applied to mgsigns(k)*mgnames(k)
C
      integer iret                      ! number of values returned per result
      integer istype                    ! results data -- subtype code
C               -1 = scalar f(t)
C               +N = profile f(x,t) of type N; N identifies which x axis.
      integer iwarn                     ! flag if arithmetic errors were caught
C               =0 means no arithmetic errors
      integer ier                       ! completion code, 0 = normal
C
C
      character*512 prog(2)
C
      prog(1)=stmt1
      prog(2)=stmt2
      call rpmgcalc(mgnames,mgsigns,nnames,prog,2,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C
      return
      end
C----------------------
C  call rpmgcalc with a three line program...
C
      subroutine rpmg3cal(mgnames,mgsigns,nnames,
     >   stmt1,stmt2,stmt3,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C  input:
      character*(*) mgnames(*)          ! list of function names
      integer mgsigns(*)                ! associated signs (+/- 1)
      integer nnames                    ! number of names in list
C
      character*(*) stmt1,stmt2,stmt3   ! three line program
C
      integer ibufsize                  ! maximum size of each result
C
C  output:
C
      real zmgbuf(ibufsize,*)           ! buffer array -- space for each result
C
C  zmgbuf(1...ret,k) is results of "prog..." applied to mgsigns(k)*mgnames(k)
C
      integer iret                      ! number of values returned per result
      integer istype                    ! results data -- subtype code
C               -1 = scalar f(t)
C               +N = profile f(x,t) of type N; N identifies which x axis.
      integer iwarn                     ! flag if arithmetic errors were caught
C               =0 means no arithmetic errors
      integer ier                       ! completion code, 0 = normal
C
C
      character*512 prog(3)
C
      prog(1)=stmt1
      prog(2)=stmt2
      prog(3)=stmt3
      call rpmgcalc(mgnames,mgsigns,nnames,prog,3,zmgbuf,
     >   ibufsize,iret,istype,iwarn,ier)
C
      return
      end
