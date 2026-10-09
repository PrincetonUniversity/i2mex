      logical function tdb_xyprof_present(d,zid)
 
      ! return TRUE if data for id (ZID) exits; otherwise FALSE
      ! no error or warning if ZID is not a known name
 
      use trdatbuf_obj
 
      implicit NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
      type (trdatbuf) :: d
 
      character*(*), intent(in) :: zid  ! (3 letter) profile ID, e.g. "PSI"
      !---------------------------------------
      character*10 zidi
      !---------------------------------------
 
      zidi = zid
      call uupper(zidi)
 
      tdb_xyprof_present = .FALSE.
 
C=>TRDATGEN+LOGIC
!
!    ******************************************
!    * TRDATGEN GENERATED CODE -- DO NOT EDIT *
!    ******************************************
!
!    code generation ends at C=>TRDATGEN- line
!
!
      if(zidi.eq.'FD0') then
         tdb_xyprof_present = (d%LFFD0.gt.0)
      endif
 
      if(zidi.eq.'FDB') then
         tdb_xyprof_present = (d%LFFDB.gt.0)
      endif
 
      if(zidi.eq.'FDP') then
         tdb_xyprof_present = (d%LFFDP.gt.0)
      endif
 
      if(zidi.eq.'FDQ') then
         tdb_xyprof_present = (d%LFFDQ.gt.0)
      endif
 
      if(zidi.eq.'FDR') then
         tdb_xyprof_present = (d%LFFDR.gt.0)
      endif
 
      if(zidi.eq.'FDS') then
         tdb_xyprof_present = (d%LFFDS.gt.0)
      endif
 
      if(zidi.eq.'PSI') then
         tdb_xyprof_present = (d%LFPSI.gt.0)
      endif
 
C
C=>TRDATGEN-
 
      return
      end
 
C===========================================================================
      subroutine tdb_xysizes(d,zid,ztyp,inx,iny,ierr)
 
      !  find grid sizes in trdatbuf object (d) for profile (zid) of type
      !    ztyp='RZ' for f(R,Z); ztyp='Ex' for f(E,x).
 
      use trdatbuf_obj
 
      implicit NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
      type (trdatbuf) :: d
 
      character*(*), intent(in) :: zid  ! (3 letter) profile ID, e.g. "PSI"
      character*(*), intent(in) :: ztyp ! (2 letter) profile type:
      !   "RZ" for f(R,Z); "Ex" for f(E,x); case blind comparison is used
 
      integer, intent(out) :: inx   ! size of first grid e.g. R if type f(R,Z)
      integer, intent(out) :: iny   ! size of 2nd grid e.g. Z if type f(R,Z)
 
      integer, intent(out) :: ierr  ! completion code (0=OK)
 
      !---------------------------------------
      ! local information:
 
      character*10 :: zidi
      character*3 :: ztypi
 
      integer :: nonlin,lunmsg_tdb,imatch
 
      !---------------------------------------
 
      nonlin = lunmsg_tdb(0)
 
      zidi = zid
      call uupper(zidi)
 
      ztypi = ztyp
      call uupper(ztypi)
 
      imatch = 0
 
c  set initial values for output variables...
 
      ierr=0
      inx=0
      iny=0
C
C=>TRDATGEN+NUM
!
!    ******************************************
!    * TRDATGEN GENERATED CODE -- DO NOT EDIT *
!    ******************************************
!
!    code generation ends at C=>TRDATGEN- line
!
!
C  FD0
      if(zidi.eq.'FD0') then
         imatch=1
         if(d%LFFD0.le.0) then
            write(nonlin,*) ' %tdb_xysizes: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xysizes subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         inx = d%NEFD0
         iny = d%NXFD0
      endif
C  FDB
      if(zidi.eq.'FDB') then
         imatch=1
         if(d%LFFDB.le.0) then
            write(nonlin,*) ' %tdb_xysizes: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xysizes subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         inx = d%NEFDB
         iny = d%NXFDB
      endif
C  FDP
      if(zidi.eq.'FDP') then
         imatch=1
         if(d%LFFDP.le.0) then
            write(nonlin,*) ' %tdb_xysizes: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xysizes subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         inx = d%NEFDP
         iny = d%NXFDP
      endif
C  FDQ
      if(zidi.eq.'FDQ') then
         imatch=1
         if(d%LFFDQ.le.0) then
            write(nonlin,*) ' %tdb_xysizes: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xysizes subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         inx = d%NEFDQ
         iny = d%NXFDQ
      endif
C  FDR
      if(zidi.eq.'FDR') then
         imatch=1
         if(d%LFFDR.le.0) then
            write(nonlin,*) ' %tdb_xysizes: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xysizes subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         inx = d%NEFDR
         iny = d%NXFDR
      endif
C  FDS
      if(zidi.eq.'FDS') then
         imatch=1
         if(d%LFFDS.le.0) then
            write(nonlin,*) ' %tdb_xysizes: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xysizes subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         inx = d%NEFDS
         iny = d%NXFDS
      endif
C  PSI
      if(zidi.eq.'PSI') then
         imatch=1
         if(d%LFPSI.le.0) then
            write(nonlin,*) ' %tdb_xysizes: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'RZ') then
            write(nonlin,*) 'error in tdb_xysizes subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: RZ'
            ierr=1
            return
         endif
         inx = d%NRPSI
         iny = d%NZPSI
      endif
C
C=>TRDATGEN-
C
      if(imatch.eq.0) then
         write(nonlin,*) 'error in tdb_xysizes subroutine:'
         write(nonlin,*) ' ID '//trim(zidi)//' not recognized.'
         write(nonlin,*) ' ID of f(t,x,y) profile was expected.'
         ierr=2
      endif
 
      return
      end
 
C===========================================================================
      subroutine tdb_xygrids(d,zid,ztyp,x,inxi,inxgot,y,inyi,inygot,
     >     ierr)
 
      !  retrieve grids in trdatbuf object (d) for profile (zid) of type
      !    ztyp='RZ' for f(R,Z); ztyp='Ex' for f(E,x).
 
      use trdatbuf_obj
 
      implicit NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
      type (trdatbuf) :: d
 
      character*(*), intent(in) :: zid  ! (3 letter) profile ID, e.g. "PSI"
      character*(*), intent(in) :: ztyp ! (2 letter) profile type:
      !   "RZ" for f(R,Z); "Ex" for f(E,x); case blind comparison is used
 
      integer, intent(in) :: inxi   ! size of X array provided
      real*8, intent(out) :: x(inxi) ! X array
      integer, intent(out) :: inxgot ! number of values returned
 
      integer, intent(in) :: inyi   ! size of Y array provided
      real*8, intent(out) :: y(inyi) ! Y array
      integer, intent(out) :: inygot ! number of values returned
 
      integer, intent(out) :: ierr  ! completion code (0=OK)
 
      !---------------------------------------
      ! local information:
 
      character*10 :: zidi
      character*3 :: ztypi
 
      integer :: nonlin,lunmsg_tdb,imatch
      integer :: iloc1,iloc2
 
      !---------------------------------------
 
      nonlin = lunmsg_tdb(0)
 
      zidi = zid
      call uupper(zidi)
 
      ztypi = ztyp
      call uupper(ztypi)
 
      imatch = 0
 
c  set initial values for output variables...
 
      ierr=0
      inxgot=0
      inygot=0
c
      x=0.0d0
      y=0.0d0
C
C=>TRDATGEN+ARRAY
!
!    ******************************************
!    * TRDATGEN GENERATED CODE -- DO NOT EDIT *
!    ******************************************
!
!    code generation ends at C=>TRDATGEN- line
!
!
C  FD0
      if(zidi.eq.'FD0') then
         imatch=1
         if(d%LFFD0.le.0) then
            write(nonlin,*) ' %tdb_xygrids: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         if(inxi.lt.d%NEFD0) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' x array size too small for: NEFD0'
            write(nonlin,*) ' got: ',inxi,' need: ',d%NEFD0
            ierr=1
            return
         endif
         if(inyi.lt.d%NXFD0) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' y array size too small for: NXFD0'
            write(nonlin,*) ' got: ',inyi,' need: ',d%NXFD0
            ierr=1
            return
         endif
         iloc1 = d%LEFD0
         inxgot = d%NEFD0
         iloc2 = iloc1 + inxgot - 1
         x(1:inxgot) = d%datbuf(iloc1:iloc2)
         iloc1 = d%LXFD0
         inygot = d%NXFD0
         iloc2 = iloc1 + inygot - 1
         y(1:inygot) = d%datbuf(iloc1:iloc2)
      endif
C  FDB
      if(zidi.eq.'FDB') then
         imatch=1
         if(d%LFFDB.le.0) then
            write(nonlin,*) ' %tdb_xygrids: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         if(inxi.lt.d%NEFDB) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' x array size too small for: NEFDB'
            write(nonlin,*) ' got: ',inxi,' need: ',d%NEFDB
            ierr=1
            return
         endif
         if(inyi.lt.d%NXFDB) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' y array size too small for: NXFDB'
            write(nonlin,*) ' got: ',inyi,' need: ',d%NXFDB
            ierr=1
            return
         endif
         iloc1 = d%LEFDB
         inxgot = d%NEFDB
         iloc2 = iloc1 + inxgot - 1
         x(1:inxgot) = d%datbuf(iloc1:iloc2)
         iloc1 = d%LXFDB
         inygot = d%NXFDB
         iloc2 = iloc1 + inygot - 1
         y(1:inygot) = d%datbuf(iloc1:iloc2)
      endif
C  FDP
      if(zidi.eq.'FDP') then
         imatch=1
         if(d%LFFDP.le.0) then
            write(nonlin,*) ' %tdb_xygrids: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         if(inxi.lt.d%NEFDP) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' x array size too small for: NEFDP'
            write(nonlin,*) ' got: ',inxi,' need: ',d%NEFDP
            ierr=1
            return
         endif
         if(inyi.lt.d%NXFDP) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' y array size too small for: NXFDP'
            write(nonlin,*) ' got: ',inyi,' need: ',d%NXFDP
            ierr=1
            return
         endif
         iloc1 = d%LEFDP
         inxgot = d%NEFDP
         iloc2 = iloc1 + inxgot - 1
         x(1:inxgot) = d%datbuf(iloc1:iloc2)
         iloc1 = d%LXFDP
         inygot = d%NXFDP
         iloc2 = iloc1 + inygot - 1
         y(1:inygot) = d%datbuf(iloc1:iloc2)
      endif
C  FDQ
      if(zidi.eq.'FDQ') then
         imatch=1
         if(d%LFFDQ.le.0) then
            write(nonlin,*) ' %tdb_xygrids: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         if(inxi.lt.d%NEFDQ) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' x array size too small for: NEFDQ'
            write(nonlin,*) ' got: ',inxi,' need: ',d%NEFDQ
            ierr=1
            return
         endif
         if(inyi.lt.d%NXFDQ) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' y array size too small for: NXFDQ'
            write(nonlin,*) ' got: ',inyi,' need: ',d%NXFDQ
            ierr=1
            return
         endif
         iloc1 = d%LEFDQ
         inxgot = d%NEFDQ
         iloc2 = iloc1 + inxgot - 1
         x(1:inxgot) = d%datbuf(iloc1:iloc2)
         iloc1 = d%LXFDQ
         inygot = d%NXFDQ
         iloc2 = iloc1 + inygot - 1
         y(1:inygot) = d%datbuf(iloc1:iloc2)
      endif
C  FDR
      if(zidi.eq.'FDR') then
         imatch=1
         if(d%LFFDR.le.0) then
            write(nonlin,*) ' %tdb_xygrids: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         if(inxi.lt.d%NEFDR) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' x array size too small for: NEFDR'
            write(nonlin,*) ' got: ',inxi,' need: ',d%NEFDR
            ierr=1
            return
         endif
         if(inyi.lt.d%NXFDR) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' y array size too small for: NXFDR'
            write(nonlin,*) ' got: ',inyi,' need: ',d%NXFDR
            ierr=1
            return
         endif
         iloc1 = d%LEFDR
         inxgot = d%NEFDR
         iloc2 = iloc1 + inxgot - 1
         x(1:inxgot) = d%datbuf(iloc1:iloc2)
         iloc1 = d%LXFDR
         inygot = d%NXFDR
         iloc2 = iloc1 + inygot - 1
         y(1:inygot) = d%datbuf(iloc1:iloc2)
      endif
C  FDS
      if(zidi.eq.'FDS') then
         imatch=1
         if(d%LFFDS.le.0) then
            write(nonlin,*) ' %tdb_xygrids: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'EX') then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: Ex'
            ierr=1
            return
         endif
         if(inxi.lt.d%NEFDS) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' x array size too small for: NEFDS'
            write(nonlin,*) ' got: ',inxi,' need: ',d%NEFDS
            ierr=1
            return
         endif
         if(inyi.lt.d%NXFDS) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' y array size too small for: NXFDS'
            write(nonlin,*) ' got: ',inyi,' need: ',d%NXFDS
            ierr=1
            return
         endif
         iloc1 = d%LEFDS
         inxgot = d%NEFDS
         iloc2 = iloc1 + inxgot - 1
         x(1:inxgot) = d%datbuf(iloc1:iloc2)
         iloc1 = d%LXFDS
         inygot = d%NXFDS
         iloc2 = iloc1 + inygot - 1
         y(1:inygot) = d%datbuf(iloc1:iloc2)
      endif
C  PSI
      if(zidi.eq.'PSI') then
         imatch=1
         if(d%LFPSI.le.0) then
            write(nonlin,*) ' %tdb_xygrids: no '//trim(zidi)//' data.'
            return
         endif
         if(ztypi.ne.'RZ') then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' value of ztyp argument: ',ztyp
            write(nonlin,*) ' for ',zidi,' profile expected: RZ'
            ierr=1
            return
         endif
         if(inxi.lt.d%NRPSI) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' x array size too small for: NRPSI'
            write(nonlin,*) ' got: ',inxi,' need: ',d%NRPSI
            ierr=1
            return
         endif
         if(inyi.lt.d%NZPSI) then
            write(nonlin,*) 'error in tdb_xygrids subroutine:'
            write(nonlin,*) ' y array size too small for: NZPSI'
            write(nonlin,*) ' got: ',inyi,' need: ',d%NZPSI
            ierr=1
            return
         endif
         iloc1 = d%LRPSI
         inxgot = d%NRPSI
         iloc2 = iloc1 + inxgot - 1
         x(1:inxgot) = d%datbuf(iloc1:iloc2)
         iloc1 = d%LZPSI
         inygot = d%NZPSI
         iloc2 = iloc1 + inygot - 1
         y(1:inygot) = d%datbuf(iloc1:iloc2)
      endif
C
C=>TRDATGEN-
C
      if(imatch.eq.0) then
         write(nonlin,*) 'error in tdb_xygrids subroutine:'
         write(nonlin,*) ' ID '//trim(zidi)//' not recognized.'
         write(nonlin,*) ' ID of f(t,x,y) profile was expected.'
         ierr=2
      endif
 
      return
      end
 
