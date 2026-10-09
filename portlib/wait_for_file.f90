subroutine wait_for_file(luntmp,filename,luntty,idelete,istry,iswarn,isfail, &
     istat)

  ! look for a file "filename" every "istry" seconds.
  ! if not found, issue warning every "iswarn" seconds.
  ! time out after "isfail" seconds.
  !
  ! if idelete=1, delete the file when it is found.
  ! if idelete=0, the file is left open on exit.
  !
  ! return istat=0 if file was found; istat=1 if timed out.

  ! defaults:  istry = 0 => try once per second
  !            iswarn = 0 => warn once every 60 seconds
  !            isfail = 0 => timeout after 300 seconds

  implicit NONE

  !----------------
  !  arguments
  integer, intent(in) :: luntmp           ! Fortran LUN for opening "filename".
  character*(*), intent(in) :: filename   ! name of file sought
  integer, intent(in) :: luntty           ! Fortran LUN for messages
  integer, intent(in) :: idelete          ! =1 to delete file once found
  integer, intent(in) :: istry            ! try every "istry" seconds
  integer, intent(in) :: iswarn           ! warn every "iswarn" seconds
  integer, intent(in) :: isfail           ! time out after "isfail" seconds

  integer, intent(out) :: istat           ! completion status, 0=OK, 1=time out

  !----------------
  !  local
  character*160 file_list(1)
  integer, parameter :: inlist = 1
  integer :: igot

  !----------------

  file_list(1)=filename

  call wait_for_file_list(luntmp,inlist,file_list, &
       luntty,idelete,istry,iswarn,isfail,igot,istat)

end subroutine wait_for_file

subroutine wait_for_file_list(luntmp,inlist,file_list, &
     luntty,idelete,istry,iswarn,isfail, &
     igot,istat)

  ! look for any of a list of files every "istry" seconds.
  ! if not found, issue warning every "iswarn" seconds.
  ! time out after "isfail" seconds.
  !
  ! if idelete=1, delete the file when it is found.
  ! if idelete=0, the file is left open on exit.
  !
  ! return igot=<index to file that was opened> or 0 if timed out
  ! return istat=0 if file was found; istat=1 if timed out.

  ! defaults:  istry = 0 => try once per second
  !            iswarn = 0 => warn once every 60 seconds
  !            isfail = 0 => timeout after 300 seconds

  implicit NONE

  !----------------
  !  arguments
  integer, intent(in) :: luntmp           ! Fortran LUN for opening "filename".
  integer, intent(in) :: inlist           ! #of files in list
  character*(*), intent(in) :: file_list(inlist)   ! names of files sought
  integer, intent(in) :: luntty           ! Fortran LUN for messages
  integer, intent(in) :: idelete          ! =1 to delete file once found
  integer, intent(in) :: istry            ! try every "istry" seconds
  integer, intent(in) :: iswarn           ! warn every "iswarn" seconds
  integer, intent(in) :: isfail           ! time out after "isfail" seconds

  integer, intent(out) :: igot            ! # of file that was opened, or 0
  integer, intent(out) :: istat           ! completion status, 0=OK, 1=time out

  !----------------
  !  local

  integer :: ilist,itry,iwarn,ifail,ictot,icwarn,io
  !----------------

  itry=max(0,istry); if(itry.eq.0) itry=1
  iwarn=max(0,iswarn); if(iwarn.eq.0) iwarn=60
  ifail=max(0,isfail); if(ifail.eq.0) ifail=300

  igot=0
  istat=0
  close(unit=luntmp,iostat=io)

  ictot=-1
  icwarn=-1
  do
     ictot=ictot+1
     icwarn=icwarn+1

     do ilist=1,inlist
        open(unit=luntmp,file=trim(file_list(ilist)),status='old',iostat=io)
        if(io.eq.0) exit
     enddo
     if(io.eq.0) exit

     close(unit=luntmp,iostat=io)
  
     if(ictot.eq.ifail) then
        write(luntty,*) ' %wait_for_file_list: after ',ifail,' seconds:'
        write(luntty,*) '  TIMED OUT waiting for file: exiting!'
        istat=1
        igot=0
        io=999
        exit
     endif

     if((icwarn.eq.iwarn).or.(ictot.le.3)) then
        if(icwarn.eq.iwarn) icwarn=0
        if(ictot.eq.0) then
           write(luntty,*) ' %wait_for_file_list: looking for any of the following files:'
           do ilist=1,inlist
              write(luntty,*) '    ',ilist,'.  ',trim(file_list(ilist))
           enddo
        else
           write(luntty,*) '  still waiting after ',ictot,' seconds...'
        endif
     endif

     call sleep(itry) ! wait before retry

  enddo

  if((istat.eq.0).and.(io.eq.0)) then
     igot=ilist
     if(idelete.eq.1) then
        close(unit=luntmp,status='delete')
     endif
  endif

end subroutine wait_for_file_list
