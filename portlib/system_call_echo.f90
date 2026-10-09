integer function system_call_echo(cmd)
  implicit none
  character(len=*), intent(in) :: cmd   ! input shell command

  ! issue a system command and capture its output with ">" to a
  ! temporary file.  Print the contents of the file to standard output
  ! and then delete the file.

  character(len=40) :: tmpfile
  character(len=20) :: pidstr
  character(len=120) :: lbuf
  integer :: ilenp,ilun,istat

  !--------------------------------
  ! form temporary filename

  call sget_pid_str(pidstr,ilenp)

  tmpfile = pidstr(1:ilenp)//'_sys.err'

  ! execute

  call execute_command_line(trim(cmd)//' > '//trim(tmpfile),exitstat=istat)

  ! return system call status
  system_call_echo = istat

  ! try to open file 

  write(6,*) ' '
  write(6,*) ' %system_call_echo: '//trim(cmd)
  write(6,*) ' %system_call_echo: status: ',istat
  write(6,*) ' %system_call_echo: stdout:'
  write(6,*) ' '

  call find_io_unit(ilun)

  open(unit=ilun,file=trim(tmpfile),status='old',iostat=istat)
  if(istat.eq.0) then
    ! echo file contents
    do
      read(ilun,'(A)',iostat=istat) lbuf
      if(istat.ne.0) exit
      write(6,'(1x,A)') trim(lbuf)
    end do
    ! close and delete file
    close(unit=ilun,status='delete')
  endif

end function system_call_echo
