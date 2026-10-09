      subroutine tgetscal(luout,idrun,abbrev,maxtimes,label,units,
     >                    ntimes,times,sdata,ierr)
C
      use datmgr_mod
      use cplotr_mod
C
C  dmc 18 Nov 1997
C  transfer scalar data (already in memory) to calling argument arrays
C  ...support routine for "trscalar".
C
C  input arguments:
C
      integer luout                     ! lun for messages
      integer idrun                     ! run COMMON data index
      character*(*) abbrev              ! function id
      integer maxtimes                  ! max no. of time pts. (array dim)
C
C  output arguments:
      character*(*) label               ! function label
      character*(*) units               ! function units
      integer ntimes                    ! number of times (actual)
      real times(maxtimes)              ! the timebase
      real sdata(maxtimes)              ! the data
C
      integer ierr                      ! completion code, 0 = normal
C
      character*50 zlbl
C
      ierr=0
      ntimes=0
      label=' '
      units=' '
C
      zlbl = rlbl(idrun)(1:lrlbl(idrun))//'!%F(T)'
C
      call dmdloc(zlbl,ind,isize,ipt)
      if(ipt.eq.0) then
C  if we get here there is a bug -- see subroutine tconnect...
         ierr=1
         write(luout,*) '?trscalar -- f(t) data never read: ',
     >      runid_x(idrun)
         go to 1000
      endif
C
      inumf=nft_x(idrun)
      inumt=ntt_x(idrun)
      if(inumt.gt.maxtimes) then
         ierr=2
         write(luout,9901) maxtimes,inumt
 9901    format(/' ?trscalar -- not enough room in passed arrays'/
     >      '  passed array dimension (maxtimes)  = ',i6/
     >      '  number of time points in data file = ',i6)
         go to 1000
      endif
C
      do if=1,inumf
         if(abbrev.eq.abt_x(if,idrun)) go to 100
      enddo
C
C  not found in list
C
      ierr=3
      write(luout,9902) abbrev,runid_x(idrun)
 9902 format(/' %trscalar -- scalar id "',a,'" not found in ',a)
      go to 1000
C
C  got it
C
 100  continue
      ifcn=if
      label=labelt_x(ifcn,idrun)
      units=unitst_x(ifcn,idrun)
      ntimes=inumt
C
      iadr0=ipt+ntimes*(ifcn-1)
      do it=1,ntimes
         iadr=iadr0+it-1
         sdata(it)=datbuf(iadr)
         times(it)=time_x(it,idrun)
      enddo
C
 1000 continue
      return
      end
C
C =============================================================
C
      subroutine tgetscal_connect(luout,idrun,abbrev,
     >                    label,units,ifcn,
     >                    ntimes,ierr)
C
      use datmgr_mod
      use cplotr_mod
C
C  rga 30 Jan 2009
C  dmc 18 Nov 1997
C  transfer scalar label and dimension data to calling argument arrays
C  ...support routine for "trscalar_connect".
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
      ierr=0
      ifcn=0
      ntimes=0
      label=' '
      units=' '
C
      zlbl = rlbl(idrun)(1:lrlbl(idrun))//'!%F(T)'
C
      call dmdloc(zlbl,ind,isize,ipt)
      if(ipt.eq.0) then
C  if we get here there is a bug -- see subroutine tconnect...
         ierr=1
         write(luout,*) '?trscalar_connect -- f(t) data never read: ',
     >      runid_x(idrun)
         go to 1000
      endif
C
      inumf=nft_x(idrun)
      inumt=ntt_x(idrun)
C
      do if=1,inumf
         if(abbrev.eq.abt_x(if,idrun)) go to 100
      enddo
C
C  not found in list
C
      ierr=3
      write(luout,9902) abbrev,runid_x(idrun)
 9902 format(/' %trscalar_connect -- scalar id "',a,
     &     '" not found in ',a)
      go to 1000
C
C  got it
C
 100  continue
      ifcn=if
      label=labelt_x(ifcn,idrun)
      units=unitst_x(ifcn,idrun)
      ntimes=inumt
C
 1000 continue
      return
      end
C
C ====================================================
C
      subroutine tgetscal_fetch(luout,idrun,abbrev,ifcn,
     >                    maxtimes,times,sdata,ierr)
C
      use datmgr_mod
      use cplotr_mod
C
C  rga 30 Jan 2009
C  dmc 18 Nov 1997
C  transfer scalar data to calling argument arrays
C  ...support routine for "trscalar_fetch".
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
C
C  output arguments:
      real times(maxtimes)              ! the timebase
      real sdata(maxtimes)              ! the data
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
      zlbl = rlbl(idrun)(1:lrlbl(idrun))//'!%F(T)'
C
      call dmdloc(zlbl,ind,isize,ipt)
      if(ipt.eq.0) then
C  if we get here there is a bug -- see subroutine tconnect...
         ierr=1
         write(luout,*) '?trscalar -- f(t) data never read: ',
     >      runid_x(idrun)
         go to 1000
      endif
C
      inumf=nft_x(idrun)
      inumt=ntt_x(idrun)
      if(inumt.gt.maxtimes) then
         ierr=2
         write(luout,9901) maxtimes,inumt
 9901    format(/' ?trscalar_fetch -- not enough room in passed arrays'/
     >      '  passed array dimension (maxtimes)  = ',i6/
     >      '  number of time points in data file = ',i6)
         go to 1000
      endif
C
      if (ifcn>=1 .and. ifcn<=inumf) then
         if(abbrev.eq.abt_x(ifcn,idrun)) go to 100
      end if

C
C  not found in list
C
      ierr=3
      write(luout,9902) abbrev,runid_x(idrun),ifcn
 9902 format(/' %trscalar_fetch -- scalar id "',a,'" not found in ',
     >     a, ' at expected index',i6)
      go to 1000
C
C  got it
C
 100  continue
C
      iadr0=ipt+inumt*(ifcn-1)
      do it=1,inumt
         iadr=iadr0+it-1
         sdata(it)=datbuf(iadr)
         times(it)=time_x(it,idrun)
      enddo
C
 1000 continue
      return
      end
