subroutine mds_cache_twr(ixflag,iasize,intt,intr,tt,tr,ier)

  ! write timebase data to cache

  use cplotr_mod
  implicit NONE

  integer, intent(in) :: ixflag   ! primary or 2ndary run indicator
  integer, intent(in) :: iasize   ! timebase array size

  integer, intent(in) :: intt     ! number of f(t) time pts
  integer, intent(in) :: intr     ! number of f(x,t) time pts

  real :: tt(iasize)              ! f(t) timebase
  real :: tr(iasize)              ! f(x,t) timebase

  integer, intent(out) :: ier     ! completion code (0=OK)

  !--------------------------------------------------------
  integer :: lunzer
  character*200 :: zfiln
  !--------------------------------------------------------
  !  write the timebase sizes first...

  call mds_cache_fname(ixflag,'numtimes.DAT',zfiln,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?mds_cache_twr: mds_cache_fname error.'
     return
  endif

  open(unit=lun_tf,file=zfiln,status='unknown',iostat=ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?mds_cache_twr: open failure: '//trim(zfiln)
     return
  endif

  write(lun_tf,'(2(2x,i8))') intt,intr
  close(unit=lun_tf)

  !------------
  !  OK now write the timebase data

  call mds_cache_write(ixflag,'time1d',-2,tt,intt,ier)
  if(ier.eq.0) then
     call mds_cache_write(ixflag,'time2d',-3,tr,intr,ier)
  endif

end subroutine mds_cache_twr

subroutine mds_cache_trd_size(ixflag,intt,intr,ier)

  ! read timebase sizes from cache

  use cplotr_mod
  implicit NONE

  integer, intent(in) :: ixflag   ! primary or 2ndary run indicator

  integer, intent(out) :: intt    ! number of f(t) time pts
  integer, intent(out) :: intr    ! number of f(x,t) time pts

  integer, intent(out) :: ier     ! completion code (0=OK)
  !  ier.ne.0 usually indicates a cache miss

  !--------------------------------------------------------
  integer :: lunzer
  character*200 :: zfiln
  !--------------------------------------------------------
  !  read the timebase sizes first...

  intt=0
  intr=0

  call mds_cache_fname(ixflag,'numtimes.DAT',zfiln,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?mds_cache_trd: mds_cache_fname error.'
     return
  endif

  open(unit=lun_tf,file=zfiln,status='old',iostat=ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' %mds_cache_trd: no timebase size file.'
     return
  endif

  read(lun_tf,'(2(2x,i8))') intt,intr
  close(unit=lun_tf)

end subroutine mds_cache_trd_size

subroutine mds_cache_trd(ixflag,iasize,intt,intr,tt,tr,ier)

  ! read timebase data from cache

  use cplotr_mod
  implicit NONE

  integer, intent(in) :: ixflag   ! primary or 2ndary run indicator
  integer, intent(in) :: iasize   ! timebase array size -- if 0 just read sizes

  integer, intent(out) :: intt    ! number of f(t) time pts
  integer, intent(out) :: intr    ! number of f(x,t) time pts

  real :: tt(iasize)              ! f(t) timebase
  real :: tr(iasize)              ! f(x,t) timebase

  integer, intent(out) :: ier     ! completion code (0=OK)
  !  ier.ne.0 usually indicates a cache miss

  !--------------------------------------------------------
  integer :: lunzer
  character*200 :: zfiln
  !--------------------------------------------------------
  !  read the timebase sizes first...

  intt=0
  intr=0

  call mds_cache_fname(ixflag,'numtimes.DAT',zfiln,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?mds_cache_trd: mds_cache_fname error.'
     return
  endif

  open(unit=lun_tf,file=zfiln,status='old',iostat=ier)
  if(ier.ne.0) then
     ! return silently, if iasize=0
     if(iasize.gt.0) then
        write(lunzer(0),*) ' %mds_cache_trd: no timebase size file.'
     endif
     return
  endif

  read(lun_tf,'(2(2x,i8))') intt,intr
  close(unit=lun_tf)

  if(iasize.eq.0) return

  !------------
  !  OK now read the timebase data

  call mds_cache_read(ixflag,'time1d',-2,tt,intt,ier)
  if(ier.eq.0) then
     call mds_cache_read(ixflag,'time2d',-3,tr,intr,ier)
  endif

end subroutine mds_cache_trd
