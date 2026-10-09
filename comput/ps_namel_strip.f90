subroutine ps_namel_strip(ifile,ps_tmpfile,luntmp,iout,ierr)

  ! strip comments from commented namelist
  !   commented namelist: ifile (to be read on luni)
  !   stripped namelist: ps_tmpfile (to be written on luntmp)

  ! The following are comments:  (a) strings to the right of the
  ! exclamation point "!" comment character (not imbedded in quotes);
  ! (b) lines that are all blank after "!" comments are stripped.

  !--------------
  ! arguments:

  character(len=*), intent(in) :: ifile ! file containing unstripped namelist.

  character(len=*), intent(in) :: ps_tmpfile ! file to contain stripped namelist.

  integer, intent(out) :: luntmp     ! I/O unit number OUTPUT of open
  !  file containing fortran-readable namelist

  integer, intent(in) :: iout        ! I/O unit for error messages
  integer, intent(out) :: ierr       ! status code (0=OK)

  !--------------
  ! local:

  character(len=150) :: buf
  character(len=200) :: my_tmpfile
  character(len=1), parameter :: squot = "'"
  character(len=1), parameter :: dquot = '"'
  character(len=1), parameter :: cchar = "!"
  character(len=1) :: cquot

  integer :: ic,icmt,ilnb,ilz
  integer :: luni

  !--------------
  !  MOD DMC: if input ps_tmpfile is blank, generate a temporary filename here
  !--------------

  if(ps_tmpfile.eq.' ') then
    call tmpfile('ps_namel_strip_',my_tmpfile,ilz)
  else
    my_tmpfile = ps_tmpfile
  endif

  luntmp=-1

  call find_io_unit(luni)
  open(unit=luni,file=ifile,status='old',iostat=ierr)
  if(ierr.ne.0) then
    write(iout,*) ' ?ps_namel_strip: namelist input file open failure.'
    write(iout,*) '  filename was: '//trim(ifile)
    return
  endif

  call find_io_unit(luntmp)
  open(unit=luntmp,file=trim(my_tmpfile),status='unknown',position='rewind',iostat=ierr)
  if(ierr.ne.0) then
    write(iout,*) ' ?ps_namel_strip: stripped namelist output file open failure.'
    write(iout,*) '  filename was: '//trim(my_tmpfile)
    return
  endif

  cquot=' '

  do 
    read(luni,'(A)',end=100) buf
    ilnb = len(trim(buf))
    if(ilnb.eq.0) cycle  ! all blank

    icmt=0
    do ic=1,ilnb
      if(cquot.ne.' ') then
        ! looking for closing quote
        if(buf(ic:ic).eq.cquot) then
          cquot=' '
        endif
        ! replace NULLs or other unprintables inside quotes
        if((ichar(buf(ic:ic)).lt.ichar(' ')).or. &
             (ichar(buf(ic:ic)).eq.127)) buf(ic:ic)=' '
      else
        ! looking for "!"
        if(buf(ic:ic).eq.cchar) then
          icmt=ic
          exit
        else if(buf(ic:ic).eq.squot) then
          cquot=squot
        else if(buf(ic:ic).eq.dquot) then
          cquot=dquot
        endif
      endif
    enddo

    if(icmt.gt.0) then
      buf(icmt:ilnb)=' '
      ilnb = len(trim(buf))
      if(ilnb.eq.0) cycle  ! all blank after removal of comment
    endif

    if(cquot.ne.' ') then
      ! insert closing quote
      ic=min(len(buf),ilnb+1)
      buf(ic:ic)=cquot
    endif
    write(luntmp,'(A)') buf(1:ilnb)
  enddo

100 continue

  close(unit=luni)

  rewind(unit=luntmp)   ! position this file to be read by caller...

  return
end subroutine ps_namel_strip
