subroutine trx_wr_stimes(spath,ierr)

  ! write file of sawtooth times -- 1st record = #of sawteeth
  ! subsequently, one time per record

  ! use data recorded in module -- report error if none is found

  use trx_module
  implicit NONE

  character*(*), intent(in) :: spath
  integer, intent(out) :: ierr

  !------------------------------------------
  integer :: it,inum,ilun
  integer :: lunzer
  !------------------------------------------

  ierr=0
  if(.not.allocated(kevent)) then
     write(lunzer(0),*) &
          ' ?trx_wr_stimes: no sawtooth data; connect to run first.'
     ierr=1
     return
  endif

  inum=0
  do it=1,nsctime
     if(kevent(it).eq.2) inum=inum+1
  enddo

  call find_io_unit(ilun)

  open(unit=ilun,file=trim(spath),status='unknown',iostat=ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) &
          ' ?trx_wr_stimes: open failure, filename was: '//trim(spath)
     return
  endif

  write(ilun, &
       '(1x,i5,"  ! #sawteeth (1x,i5) // times (1x,1pe13.6) one per line")') &
       inum

  do it=1,nsctime
     if(kevent(it).eq.2) then
        write(ilun,'(1x,1pe13.6)') time_saw(it)
     endif
  enddo

  close(unit=ilun)

end subroutine trx_wr_stimes
