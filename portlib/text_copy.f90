subroutine text_copy(path1,path2,ierr)

  ! copy a text file (path1) to a new text file (path2).
  !   done line-by-line with fortran I/O, the output file could have
  !   have different trailing white space than the original
  !   **line width limit: 200 characters**

  implicit NONE

  character*(*), intent(in) :: path1  ! file to copy *from*
  character*(*), intent(in) :: path2  ! file to copy *to*

  integer, intent(out) :: ierr  ! completion status, 0=OK
  !  ierr = 1 means: read-access to (path1) failed
  !  ierr = 2 means: write-access to (path2) failed

  !--------------------
  integer :: iluns(2),istat
  character*200 :: line
  !--------------------

  ierr = 0

  call find_io_unit_list(2,iluns)

  open(unit=iluns(1),file=path1,status='old',action='read',iostat=istat)
  if(istat.ne.0) then
     ierr=1
     return
  endif

  open(unit=iluns(2),file=path2,status='unknown',iostat=istat)
  if(istat.ne.0) then
     ierr=2
     close(unit=iluns(1))
     return
  endif

  do
     line=' '
     read(iluns(1),'(A)',end=99) line
     write(iluns(2),'(A)') trim(line)
  enddo

99 continue
  close(unit=iluns(1))
  close(unit=iluns(2))
  return

end subroutine text_copy
