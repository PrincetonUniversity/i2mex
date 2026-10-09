      subroutine tgetprof(luout,idrun,abbrev,maxtimes,maxdata,
     >                    label,units,
     >                    itype,inx,ntimes,times,data,ierr)
C
      use datmgr_mod
      use cplotr_mod
      implicit none
C
C  dmc 19 Nov 1997
C  transfer profile data to calling argument arrays
C  ...support routine for "trprofil".
C
C  the run labels and timebases are already in memory; the named
C  function data may have to be read in.
C
C  input arguments:
C
      integer luout                     ! lun for messages
      integer idrun                     ! run COMMON data index
      character*(*) abbrev              ! function id
      integer maxtimes                  ! max no. of time pts. (array dim)
      integer maxdata                   ! max no. of data pts. (array dim)
C
C  output arguments:
      character*(*) label               ! function label
      character*(*) units               ! function units
      integer ntimes                    ! number of times (actual)
      integer itype                     ! x axis "type" code
      integer inx                       ! no. of pts in x axis
      real times(maxtimes)              ! the timebase
      real data(maxdata)                ! the data
C
      integer ierr                      ! completion code, 0 = normal
C
C clf: declare after implicit none
      integer inumf, inumt, if, ifcn, ineed, ind, ipt
      integer idata, it,    iadr0, ix, iadr
C
      character*50 zlbl
C
      ierr=0
      ntimes=0
      label=' '
      units=' '
      itype=0
      inx=0
C
      zlbl = rlbl(idrun)(1:lrlbl(idrun))//'!'//abbrev
C
      inumf=nfxt_x(idrun)
      inumt=ntr_x(idrun)
      if(inumt.gt.maxtimes) then
         ierr=2
         write(luout,9901) inumt,maxtimes
 9901    format(/' ?trprofil -- not enough room in passed arrays'/
     >      '  passed array dimension (maxtimes)  = ',i6/
     >      '  number of time points in data file = ',i6)
         go to 1000
      endif
C
      do if=1,inumf
         if(abbrev.eq.abr_x(if,idrun)) go to 100
      enddo
C
C  not found in list
C
      ierr=3
      write(luout,9902) abbrev,runid_x(idrun)
 9902 format(/' %trprofil -- profile id "',a,'" not found in ',a)
      go to 1000
C
C  found the function id
C
 100  continue
      ifcn=if
      itype=itypr_x(ifcn,idrun)
      inx=nzonex_x(itype,idrun)
C
      ineed=inx*inumt
      if(ineed.gt.maxdata) then
         ierr=4
         write(luout,9903) ineed,maxdata,inx,inumt
 9903    format(/' ?trprofil -- not enough room in passed arrays'/
     >      '  passed array dimension (maxdata)   = ',i6/
     >      '  number of words needed             = ',i6/
     >      '  (',i3,' x axis pts * ',i5,' time pts)')
         go to 1000
      endif
C
C  OK -- no error
C
C  get the data
C
      lrun_x=idrun
      call dmgfxt(ifcn,ind)
      if (ntr .lt. 0) then
         ierr=1
         return
      endif
      lrun_x=0
C
      ipt=locd(ind)
C
      label=labelr_x(ifcn,idrun)
      units=unitsr_x(ifcn,idrun)
      ntimes=inumt
C
      idata=0
      do it=1,ntimes
         times(it)=time3_x(it,idrun)
         iadr0=ipt+(it-1)*inx
         do ix=1,inx
            iadr=iadr0+ix-1
            idata=idata+1
            data(idata)=datbuf(iadr)
         enddo
      enddo
C
 1000 continue
      return
      end

C
C =========================================================
C
      subroutine tgetprof_connect(luout,idrun,abbrev,
     >                    label,units,ifcn,
     >                    itype,inx,ntimes,ierr)
C
      use datmgr_mod
      use cplotr_mod
C
C  rga 30 Jan 2009
C  dmc 19 Nov 1997
C  transfer profile label and dimension data to calling arguments
C  ...support routine for "trprofil_connect".
C
C  the run labels and timebases are already in memory; the named
C  function data may have to be read in.
C
C  input arguments:
C
      integer luout                     ! lun for messages
      integer idrun                     ! run COMMON data index
      character*(*) abbrev              ! function id
C
C  output arguments:
      character*(*) label               ! function label
      character*(*) units               ! function units
      integer ifcn                      ! function index
      integer itype                     ! x axis "type" code
      integer inx                       ! no. of pts in x axis
      integer ntimes                    ! number of times (actual)
C
      integer ierr                      ! completion code, 0 = normal
C
C-----------------------------------------------
C
      character*50 zlbl
C
C-----------------------------------------------
C
      label=' '
      units=' '
      ifcn=0
      itype=0
      inx=0
      ntimes=0
      ierr=0
C
      zlbl = rlbl(idrun)(1:lrlbl(idrun))//'!'//abbrev
C
      inumf=nfxt_x(idrun)
      inumt=ntr_x(idrun)
C
      do if=1,inumf
         if(abbrev.eq.abr_x(if,idrun)) go to 100
      enddo
C
C  not found in list
C
      ierr=3
      write(luout,9902) abbrev,runid_x(idrun)
 9902 format(/' %trprofil_connect -- profile id "',a,
     &     '" not found in ',a)
      go to 1000
C
C  found the function id
C
 100  continue
      ifcn=if
      itype=itypr_x(ifcn,idrun)
      inx=nzonex_x(itype,idrun)
C
      lrun_x=0
C
      label=labelr_x(ifcn,idrun)
      units=unitsr_x(ifcn,idrun)
      ntimes=inumt
C
 1000 continue
      return
      end
C
C =========================================================
C
      subroutine tgetprof_fetch(luout,idrun,abbrev,ifcn,
     >                    maxtimes,maxdata,times,data,ierr)
C
      use datmgr_mod
      use cplotr_mod
C
C  rga 30 Jan 2009
C  dmc 19 Nov 1997
C  transfer profile data to calling argument arrays
C  ...support routine for "trprofil_fetch".
C
C  the run labels and timebases are already in memory; the named
C  function data may have to be read in.
C
C  input arguments:
C
      integer luout                     ! lun for messages
      integer idrun                     ! run COMMON data index
      character*(*) abbrev              ! function id
      integer ifcn                      ! function index
      integer maxtimes                  ! max no. of time pts. (array dim)
      integer maxdata                   ! max no. of data pts. (array dim)
C
C  output arguments:
      real times(maxtimes)              ! the timebase
      real data(maxdata)                ! the data
C
      integer ierr                      ! completion code, 0 = normal
C
C-----------------------------------------------
C
      character*50 zlbl
C
C-----------------------------------------------
C
      ierr=0
C
      zlbl = rlbl(idrun)(1:lrlbl(idrun))//'!'//abbrev
C
      inumf=nfxt_x(idrun)
      inumt=ntr_x(idrun)
      if(inumt.gt.maxtimes) then
         ierr=2
         write(luout,9901) inumt,maxtimes
 9901    format(/' ?trprofil_fetch -- not enough room in passed arrays'/
     >      '  passed array dimension (maxtimes)  = ',i8/
     >      '  number of time points in data file = ',i8)
         go to 1000
      endif
C
      if (ifcn>=1 .and. ifcn<=inumf) then
         if(abbrev.eq.abr_x(ifcn,idrun)) go to 100
      end if
C
C  not found in list
C
      ierr=3
      write(luout,9902) abbrev,runid_x(idrun),ifcn
 9902 format(/' %trprofil_fetch -- profile id "',a,'" not found in ',
     >     a,' at expected index',i6)
      go to 1000
C
C  found the function id
C
 100  continue
      itype=itypr_x(ifcn,idrun)
      inx=nzonex_x(itype,idrun)
C
      ineed=inx*inumt
      if(ineed.gt.maxdata) then
         ierr=4
         write(luout,9903) maxdata,ineed,inx,inumt
 9903    format(/' ?trprofil_fetch -- not enough room in passed arrays'/
     >      '  passed array dimension (maxdata)   = ',i9/
     >      '  number of words needed             = ',i9/
     >      '  (',i3,' x axis pts * ',i5,' time pts)')
         go to 1000
      endif
C
C  OK -- no error
C
C  get the data
C
      lrun_x=idrun
      call dmgfxt(ifcn,ind)
      lrun_x=0
C
      ipt=locd(ind)
C
      idata=0
      do it=1,inumt
         times(it)=time3_x(it,idrun)
         iadr0=ipt+(it-1)*inx
         do ix=1,inx
            iadr=iadr0+ix-1
            idata=idata+1
            data(idata)=datbuf(iadr)
         enddo
      enddo
C
 1000 continue
      return
      end

