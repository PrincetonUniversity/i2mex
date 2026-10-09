subroutine mds_cache_putdev(zfiln,zdev)

  ! store device name in cache file; ignore errors...

  use cplotr_mod
  implicit none

  character*(*), intent(in) :: zfiln   ! filename
  character*(*), intent(in) :: zdev    ! device, i.e. tokamak label

  !----------------------------
  integer :: istat,ilenz
  !----------------------------

  ilenz = len(trim(zdev))
  if(ilenz.eq.0) return

  open(unit=lun_tf,file=zfiln,status='unknown',iostat=istat)
  if(istat.eq.0) then
     write(lun_tf,'(1x,A)') trim(zdev)
     close(unit=lun_tf)
  endif

end subroutine mds_cache_putdev

subroutine mds_cache_getdev(zfiln,zdev,ilt)

  ! retrieve device name from cache file; set ilt=0 on error

  use cplotr_mod
  implicit none

  character*(*), intent(in) :: zfiln   ! filename
  character*(*), intent(out) :: zdev   ! device, i.e. tokamak label
  integer, intent(out) :: ilt          ! non-blank length of zdev; 0 on error

  !----------------------------
  integer :: istat
  !----------------------------

  zdev = ' '
  ilt = 0

  open(unit=lun_tf,file=zfiln,status='old',iostat=istat)
  if(istat.eq.0) then
     read(lun_tf,'(1x,A)') zdev
     close(unit=lun_tf)
     ilt = len(trim(zdev))
  endif

end subroutine mds_cache_getdev
