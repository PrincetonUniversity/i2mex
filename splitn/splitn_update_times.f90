subroutine splitn_update_ntimes(inum)

  ! return the number of update times in the current namelist

  use splitn_module
  implicit NONE

  integer, intent(out) :: inum

  !--------------------------

  if((nlines.le.0).or.(.not.have_database)) then
     write(6,*) ' ?splitn_update_times: no namelist has been read.'
     inum = 0
  else
     inum = nupdate
  endif

end subroutine splitn_update_ntimes

subroutine splitn_update_times(inum,ztimes,ierr)

  ! return the times of update blocks in namelist

  use splitn_module
  implicit NONE

  integer, intent(in) :: inum
  real*8, intent(out) :: ztimes(inum)
  integer, intent(out) :: ierr  ! status return code, 0=OK

  !----------------------------
  ! return ierr=1 & ztimes=0.0 if there are no update blocks.
  ! return ierr=2 & ztimes=0.0 if there are update blocks but nupdate > inum
  ! otherwise return ierr=0 and the times in ztimes(1:nupdate), and zeroes 
  !   in ztimes(j), j>nupdate.
  !----------------------------

  ztimes = 0.0d0
  ierr = 0

  if(nupdate.eq.0) then
     write(6,*) &
          ' %splitn_update_times: there are no update blocks in the namelist.'
     ierr = 1
     return
  endif

  if(inum.lt.nupdate) then
     write(6,*) &
          ' ?splitn_update_times: array dimension for update times too small:'
     write(6,*) '  received: inum=',inum,'; actual #updates = ',nupdate
     ierr = 2
     return
  endif

  ztimes(1:nupdate) = tup(1:nupdate)

end subroutine splitn_update_times

subroutine splitn_update_nxtime(prev_time,next_time)

  ! return the first update time after the passed time (prev_time)

  use splitn_module
  implicit NONE

  real*8, intent(in) :: prev_time
  real*8, intent(out) :: next_time

  !----------------------------------------
  integer :: ii
  !----------------------------------------

  next_time = 1.0d34

  do ii=nupdate,1,-1
     if(tup(ii).gt.prev_time) next_time=tup(ii)
  enddo

end subroutine splitn_update_nxtime

subroutine splitn_update_blkid(ztime,blkstr)

  ! return update block ID: "[0]" for initial block prior to any update, 
  !   "[1]" for 1st block
  !   "[2]" for 2nd block
  !      ..etc..  according to the time provided

  implicit NONE

  real*8, intent(in) :: ztime
  character*(*), intent(out) :: blkstr

  !------------------------
  integer :: inum
  real*8, dimension(:), allocatable :: ztimes
  integer :: ii,iblk,ierr
  !------------------------

  blkstr="[0]"

  call splitn_update_ntimes(inum)

  allocate(ztimes(inum)); ztimes=0.0d0

  call splitn_update_times(inum,ztimes,ierr)
  if(ierr.ne.0) return

  if(ztime.lt.ztimes(1)) then
     deallocate(ztimes)
     return
  endif

  do ii=1,inum
     if(ztimes(ii).gt.ztime) exit
     iblk=ii
  enddo

  if(iblk.lt.10) then
     write(blkstr,'("[",i1,"]")') iblk
  else if(iblk.lt.100) then
     write(blkstr,'("[",i2,"]")') iblk
  else if(iblk.lt.1000) then
     write(blkstr,'("[",i3,"]")') iblk
  else if(iblk.lt.10000) then
     write(blkstr,'("[",i4,"]")') iblk
  else if(iblk.lt.100000) then
     write(blkstr,'("[",i5,"]")') iblk
  else
     write(blkstr,'("[",i6,"]")') iblk
  endif

  deallocate(ztimes)

end subroutine splitn_update_blkid


