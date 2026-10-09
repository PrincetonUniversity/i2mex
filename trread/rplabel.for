      subroutine rplabel(abbrev,zlabel,zunits,imulti,istype)
C
C  given an RPLOT name (abbrev) return the label (zlabel) and units (zunits)
C
      use cplotr_mod

      character*(*) abbrev              ! input name (C*10)
      character*(*) zlabel              ! output label (C*32)
      character*(*) zunits              ! output units label (C*16)
C
      integer imulti                    ! output =1 if a multigraph
      integer istype                    ! output data subtype code
C
C         istype = -1   -- scalar f(t) function or multigraph
C         istype = 1,2,... -- various f(x,t) profile functions / multigraphs
C
C  if the name is not recognized,
C     zlabel = zunits = 'error' are returned, and
C     imulti = istype = -1000   are returned.
C
C------------------------------------
C
      character*10 zabbr
C------------------------------------
C
      zlabel='error'
      zunits='error'
      imulti=-1000
      istype=-1000
C
      zabbr=abbrev
      call trcaps(zabbr)
C
      ilbl=min(len(zlabel),len(labelt(1)))
      ilbu=min(len(zunits),len(unitst(1)))
C
C  try for a scalar
C
      iadr=ifind_ordr(abt,iordrt,nft,zabbr)
      if(iadr.gt.0) then
         zlabel=labelt(iadr)(1:ilbl)
         zunits=unitst(iadr)(1:ilbu)
         imulti=0
         istype=-1
         go to 100
      endif
C
C  try for a profile
C
      iadr=ifind_ordr(abr,iordrr,nfxt,zabbr)
      if(iadr.gt.0) then
         zlabel=labelr(iadr)(1:ilbl)
         zunits=unitsr(iadr)(1:ilbu)
         imulti=0
         istype=itypr(iadr)
         go to 100
      endif
C
C  try for a multigraph
C
      iadr=ifind_ordr(abb,iordrb,nbal,zabbr)
      if(iadr.gt.0) then
         zlabel=labelb(iadr)(1:ilbl)
         zunits=unitsb(iadr)(1:ilbu)
         imulti=1
         if(iintb(iadr).eq.1) then
            istype=-1
         else
            if1=iabs(ifunb(1,iadr))
            istype=itypr(if1)
         endif
         go to 100
      endif
C
C  exit
C
 100  continue
C
      if(istype.eq.-1000) then
         call zermsg(' ?rplabel:  unrecognized item name:  '//abbrev)
      endif
C
      return
      end

C-----------------------------------------------
C  new DMC June 2007: rpexist routines...

C-----------------------------------------------
      subroutine rpexist_scalar(abbrev,exist)

      !  check for existence of a scalar function

      use cplotr_mod

      implicit NONE

      character*(*), intent(in) :: abbrev ! f(t) scalar name, case insensitive
      logical, intent(out) :: exist       ! .TRUE. if it is in the data

C------------------------------------
C
      character*10 zabbr
C
      integer :: iadr,ifind_ordr
C------------------------------------
C
      zabbr=abbrev
      call trcaps(zabbr)

      iadr=ifind_ordr(abt,iordrt,nft,zabbr)

      exist = (iadr.gt.0)

      end

C-----------------------------------------------
      subroutine rpexist_profile(abbrev,exist)

      !  check for existence of a profile function

      use cplotr_mod

      implicit NONE

      character*(*), intent(in) :: abbrev ! f(x,t) profile name, case insensitive
      logical, intent(out) :: exist       ! .TRUE. if it is in the data

C------------------------------------
      character*10 zabbr
C
      integer :: iadr,ifind_ordr
C------------------------------------
C
      zabbr=abbrev
      call trcaps(zabbr)

      iadr=ifind_ordr(abr,iordrr,nfxt,zabbr)

      exist = (iadr.gt.0)

      end

C-----------------------------------------------
      subroutine rpexist_multi(abbrev,exist,istype)

      !  check for existence of a multigraph set
      !  if it exists return its type
      !    if size (#members) also desired, use rpsize_multi(...) (below)

      use cplotr_mod

      implicit NONE

      character*(*), intent(in) :: abbrev ! f(x,t) profile name, case insensitive
      logical, intent(out) :: exist       ! .TRUE. if it is in the data

      integer, intent(out) :: istype      ! type of multigraph, if exist=.TRUE.
      !  same interpretation as in subroutine rplabel

C------------------------------------
      character*10 zabbr
C
      integer :: iadr,ifind_ordr,if1
C
C------------------------------------
C
      zabbr=abbrev
      call trcaps(zabbr)

      iadr=ifind_ordr(abb,iordrb,nbal,zabbr)
      if(iadr.gt.0) then
         exist = .TRUE.
         if(iintb(iadr).eq.1) then
            istype=-1
         else
            if1=iabs(ifunb(1,iadr))
            istype=itypr(if1)
         endif
      else
         exist = .FALSE.
         istype = -1000
      endif

      end

C-----------------------------------------------
      subroutine rpsize_multi(abbrev,exist,istype,isize)

      !  check for existence of a multigraph set
      !  if it exists return its type and size

      use cplotr_mod
      implicit NONE

      character*(*), intent(in) :: abbrev ! f(x,t) profile name, case insensitive
      logical, intent(out) :: exist       ! .TRUE. if it is in the data

      integer, intent(out) :: istype      ! type of multigraph, if exist=.TRUE.
      !  same interpretation as in subroutine rplabel

      integer, intent(out) :: isize       ! #of members (0 if exist=.FALSE.)
C------------------------------------
      character*10 zabbr
C
      integer :: iadr,ifind_ordr,if1
C
C------------------------------------
C
      zabbr=abbrev
      call trcaps(zabbr)

      iadr=ifind_ordr(abb,iordrb,nbal,zabbr)
      if(iadr.gt.0) then
         exist = .TRUE.
         if(iintb(iadr).eq.1) then
            istype=-1
         else
            if1=iabs(ifunb(1,iadr))
            istype=itypr(if1)
         endif
         isize = infb(iadr)
      else
         exist = .FALSE.
         istype = -1000
         isize = 0
      endif

      end
