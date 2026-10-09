module datmgr_mod

  implicit NONE
  SAVE
  PUBLIC

  !*************** START FILE DATMGR.BLK ; GROUP DATMGR *************
  !--------------------------------------------------------------
  !  DATMGR
  !
  !  dmc -- f90 free form compatibility -- 27 Feb 2000
  !    use "!" for all comments; & @73 / & @5 for continuations ...!
  !
  !  dmc -- converted to f90 module -- 4 Nov 2009 -- data buffer is now
  !  explandable, to adapt to data storage needs encountered.
  !
  !  BLOCK FOR MANAGING A DATA BUFFER for RPLOT -- traditional TRANSP output.
  !
  !  TIME VECTOR now included here
  !
  !  BUFFER WILL ALLOW INSERTION OF VARIABLE SIZE BLOCKS OF DATA
  !  SPACE ALLOCATED ON LAST IN LAST OUT BASIS
  !
  !  BUFFER SIZE
  !============
  ! idecl:  explicitize implicit INTEGER declarations:

  INTEGER :: macc,lavail

  logical :: no_delete = .FALSE.

  integer :: lundmo=6

  integer :: nwsmin=16384 ! min "workspace" size w/in DATBUF(...)

  !============

  integer, parameter :: ndbsiz_min =   50000000 ! -> 200MB
  integer, parameter :: ndbsiz_max = 2001001000 ! -> 8GB

  !============
  ! These values refer to "priorities" used by rplot library routines
  ! to identify persistent uses of datbuf memory.  Persistent items (protected
  ! against replacement) accumulate up from the low end in DATBUF memory, and
  ! accumulate down from the high end...

  integer :: hi_end_prio = 7  ! "high end" storage priority
  integer :: lo_end_prio = 10 ! "low end" storage priority

  integer :: ndent=0

  integer :: NDBSIZ=0

  !  MAXIMUM NUMBER OF DISTINCT ENTRIES

  integer, parameter :: MAXENT=8192

  !
  !  DATBUF IS THE DATA BUFFER (FULL SIZE REAL WORDS)

  real, dimension(:), allocatable :: datbuf

  !  FOR THE JTH ENTRY IN THE BUFFER,
  !   LOCD(J)=START ADDRESS OF DATA IN BUFFER
  !   NWDS(J)=SIZE OF DATA ENTRY
  !   LACC(J)=LAST ACCESS CODE (UPDATED EACH TIME DATA IS ACCESSED)
  !   LPREV(J)=PTR TO PREVIOUS DATA ENTRY
  !   LNEXT(J)=PTR TO NEXT DATA ENTRY
  !   MPRIO(J)="PRIORITY" OF ENTRY J-- LOW-PRIORITY ITEMS CANNOT DISPLACE
  !     HIGH PRIORITY ITEMS
  !
  !   NDENT=NUMBER OF ENTRIES
  !   MACC=CURRENT ACCESS CODE
  !
  !  NOTE ENTRIES 1 AND NDENT ARE DUMMIES WITH 0 SIZE
  !    TO START SUCCESSOR/PREDECESSOR CHAINS; INITIALLY NDENT=2
  !
  !  LAVAIL= PTR TO FIRST FREE ENTRY SLOT IN LOCD
  !
  !  LUNDMO= OUTPUT CHANNEL FOR ERROR MESSAGES

  INTEGER :: LOCD(MAXENT),NWDS(MAXENT),LACC(MAXENT),LPREV(MAXENT),     &
       LNEXT(MAXENT),MPRIO(MAXENT)

  !  char string label for each buffer entry

  CHARACTER*50 DMGLBL(MAXENT)

  !-----------------------------------
  ! timebase data

  INTEGER :: NTIME=0  ! (current allocation)
  INTEGER :: NSM=0

  INTEGER :: NTCORR   ! for poplot-- error recovery counter
  logical, dimension(:), allocatable :: LTWRIT

  real, dimension(:), allocatable :: time  ! scalar data time vector
  real, dimension(:), allocatable :: time3 ! f(x,t) profile time vector
  real, dimension(:), allocatable :: xtime,work1,work2

  real, dimension(:,:), allocatable :: time_x,time3_x
  real, dimension(:,:), allocatable :: smwork,smwork2

  !-------------------------------------------------------------------

CONTAINS

  subroutine dmg_tinit

    !  (re)initialize time data

    logical :: re_init

    !-------------------------------

    re_init = (ntime.gt.0)

    call dmg_texpand(ntime)  ! if ntime==0 this allocates arrays for first time

    NTCORR=0

    if(re_init) then
       LTWRIT(1:ntime) = .TRUE.
       time(1:ntime) = 0.0
       time3(1:ntime) = 0.0
       xtime(1:ntime) = 0.0
       work1(1:ntime) = 0.0
       work2(1:ntime) = 0.0

       ! leaving time_x & time3_x undisturbed... might contain data from
       ! 2ndary runs that is best kept available.

    endif

  end subroutine dmg_tinit

  subroutine dmg_texpand(inumt)

    ! expand time vector and related arrays

    integer, intent(in) :: inumt

    ! if inumt=0: expand time vector to max(inumt_min,2*ntime)
    ! if inumt>0: expand time vector to inumt, if inumt>ntime

    !---------------------------------
    integer, parameter :: inumt_min = 2048
    integer :: inumt_new,inumt_old,inxtra,inr0

    logical, dimension(:), allocatable :: logtmp
    real, dimension(:), allocatable :: r4tmp
    real, dimension(:,:), allocatable :: r4tmp2d
    !---------------------------------

    inumt_old = ntime

    if(inumt.gt.0) then
       if(inumt.le.ntime) return
       inumt_new = max(inumt_min,inumt)

    else
       inumt_new = max(inumt_min,2*ntime)
    endif

    call iget_naxxvr(inxtra,inr0)

    ntime=inumt_new
    if(inumt_old.gt.0) then
       allocate(logtmp(inumt_old),r4tmp(inumt_old),r4tmp2d(inumt_old,inxtra))
    endif

    ! ltwrit

    if(inumt_old.gt.0) then
       logtmp = ltwrit
       deallocate(ltwrit)
    endif

    allocate(ltwrit(ntime))
    if(inumt_old.gt.0) then
       ltwrit(1:inumt_old)=logtmp
    endif
    ltwrit(inumt_old+1:ntime) = .TRUE.

    if(allocated(logtmp)) deallocate(logtmp)
    !====================

    ! time

    if(inumt_old.gt.0) then
       r4tmp = time
       deallocate(time)
    endif

    allocate(time(ntime))
    if(inumt_old.gt.0) then
       time(1:inumt_old)=r4tmp
    endif
    time(inumt_old+1:ntime)=0.0

    ! time3

    if(inumt_old.gt.0) then
       r4tmp = time3
       deallocate(time3)
    endif

    allocate(time3(ntime))
    if(inumt_old.gt.0) then
       time3(1:inumt_old)=r4tmp
    endif
    time3(inumt_old+1:ntime)=0.0

    ! xtime

    if(inumt_old.gt.0) then
       r4tmp = xtime
       deallocate(xtime)
    endif

    allocate(xtime(ntime))
    if(inumt_old.gt.0) then
       xtime(1:inumt_old)=r4tmp
    endif
    xtime(inumt_old+1:ntime)=0.0

    ! work1

    if(inumt_old.gt.0) then
       r4tmp = work1
       deallocate(work1)
    endif

    allocate(work1(ntime))
    if(inumt_old.gt.0) then
       work1(1:inumt_old)=r4tmp
    endif
    work1(inumt_old+1:ntime)=0.0

    ! work2

    if(inumt_old.gt.0) then
       r4tmp = work2
       deallocate(work2)
    endif

    allocate(work2(ntime))
    if(inumt_old.gt.0) then
       work2(1:inumt_old)=r4tmp
    endif
    work2(inumt_old+1:ntime)=0.0

    if(allocated(r4tmp)) deallocate(r4tmp)
    !====================

    ! time_x

    if(inumt_old.gt.0) then
       r4tmp2d = time_x
       deallocate(time_x)
    endif

    allocate(time_x(ntime,inxtra))
    if(inumt_old.gt.0) then
       time_x(1:inumt_old,1:inxtra)=r4tmp2d(1:inumt_old,1:inxtra)
    endif
    time_x(inumt_old+1:ntime,1:inxtra)=0.0

    ! time3_x

    if(inumt_old.gt.0) then
       r4tmp2d = time3_x
       deallocate(time3_x)
    endif

    allocate(time3_x(ntime,inxtra))
    if(inumt_old.gt.0) then
       time3_x(1:inumt_old,1:inxtra)=r4tmp2d(1:inumt_old,1:inxtra)
    endif
    time3_x(inumt_old+1:ntime,1:inxtra)=0.0

    if(allocated(r4tmp2d)) deallocate(r4tmp2d)

    ! also re-allocate smoothing workspaces (but do not preserve old contents)

    if(allocated(smwork)) then
       deallocate(smwork,smwork2)
    endif

    nsm=max(ntime,inr0)
    allocate(smwork(nsm,4),smwork2(2,nsm))

  end subroutine dmg_texpand

  subroutine dmgini

    ! init call -- no size specified

    call dmgini_sized(0)

  end subroutine dmgini

  subroutine dmgini_sized(isize)

    ! init call -- size specified

    integer, intent(in) :: isize

    !-------------------------
    integer :: j
    !-------------------------

    if((.NOT.ALLOCATED(datbuf)).or.(isize.gt.ndbsiz)) then
       write(6,*) ' dmg_datbuf_expand call from dmgini_sized: isize=',isize
       call dmg_datbuf_expand(isize)
    endif

    MACC=0
    DO J=1,MAXENT
       LOCD(J)=0
       NWDS(J)=0
       LACC(J)=0
       LPREV(J)=0
       LNEXT(J)=0
       MPRIO(J)=0
    ENDDO

    NDENT=2
    LOCD(1)=1
    LOCD(2)=NDBSIZ+1
    LNEXT(1)=2
    LPREV(2)=1

    ! INITIAL NAMES FOR DATA BLOCKS
    DMGLBL(1)='%INIT'
    DMGLBL(2)='%FINI'

    LAVAIL=3

  end subroutine dmgini_sized

  subroutine dmg_datbuf_expand(isize)

    ! allocate/expand DATBUF(...)

    integer, intent(in) :: isize

    !--------------------

    integer :: isizu,isizo,ialloc,istat,ihi,idiff,j
    real, dimension(:), allocatable :: tmp_datbuf

    !--------------------
    !  compute new DATBUF size

    isizu = ndbsiz_min
    isizu = max(isizu,min(2*ndbsiz,ndbsiz+2*ndbsiz_min))
    isizu = max(isizu,isize)
    isizu = min(isizu,ndbsiz_max)

    isizo=0
    ialloc=0
    ihi=isizu+1

    if(allocated(datbuf)) then
       isizo = size(datbuf)
       if(isizo.ge.isizu) then
          write(lundmo,*) ' old size = ',isizo
          write(lundmo,*) ' new size = ',isizu
          write(lundmo,*) ' max size = ',ndbsiz_max
          call errmsg_exit(' ?datmgr_mod(dmg_datbuf_expand): DATBUF expansion not available at ndbsiz_max')
       else
          write(lundmo,*) ' '
          write(lundmo,*) ' %datmgr_mod: expanding DATBUF(...):'
          write(lundmo,*) '  old size = ',isizo
          write(lundmo,*) '  new size = ',isizu
          write(lundmo,*) ' '
       endif

       allocate(tmp_datbuf(isizo),stat=istat)
       if(istat.ne.0) then
          call errmsg_exit(' ?datmgr_mod(dmg_datbuf_expand): tmp_datbuf ALLOCATE error!')
       endif
       ialloc=1
       tmp_datbuf = datbuf
       deallocate(datbuf)

       ! find high memory storage

       do j=1,ndent
          if(mprio(j).eq.hi_end_prio) then
             ihi=min(ihi,locd(j))
          endif
       enddo
    endif

    idiff = isizu-isizo

    ndbsiz = isizu
    allocate(datbuf(isizu),stat=istat)
    if(istat.ne.0) then
       call errmsg_exit(' ?datmgr_mod(dmg_datbuf_expand): datbuf ALLOCATE error!')
    endif

    if(ialloc.eq.1) then
       if(ihi.gt.isizu) then
          datbuf(1:isizo) = tmp_datbuf
          datbuf(isizo+1:ndbsiz)=0.0
       else
          datbuf(1:ihi-1) = tmp_datbuf(1:ihi-1)
          datbuf(ihi+idiff:ndbsiz) = tmp_datbuf(ihi:isizo)
          datbuf(ihi:ihi+idiff-1) = 0.0
          do j=1,ndent
             if(locd(j).ge.ihi) locd(j)=locd(j)+idiff
          enddo
       endif
       deallocate(tmp_datbuf)
    else
       datbuf(1:ndbsiz)=0.0
    endif

    LOCD(2)=NDBSIZ+1   ! terminator

  end subroutine dmg_datbuf_expand

end module datmgr_mod
