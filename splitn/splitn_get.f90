! MODIFICATION -- DMC March 2009 -- account for namelist update blocks
!   if the zname argument is appended with an update index reference of the
!   form "[<num>]", this means, fetch the namelist value from the <num>'th
!   update block.  If this construct is "[0]" or omitted, fetch the initial
!   namelist value (as before).

subroutine splitn_iget(zname,isize,array,ierr)

  !  get an INTEGER namelist object
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  integer, intent(out) :: array(isize)  ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = 0

  call splitn_bdecode1('splitn_iget',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  iblock=max(0,min(nupdate,iblock))

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_iget: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_iget: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'I') then
     write(6,*) '?splitn_iget: not an INTEGER object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_iget: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
  else
     ii=nint + (iblock-1)*nint_st + varlist(jj)%addr_st
  endif

  ierr = 0
  array = intbuf(ii:ii+isize-1)

end subroutine splitn_iget

subroutine splitn_lget(zname,isize,array,ierr)

  !  get a LOGICAL namelist object
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  logical, intent(out) :: array(isize)  ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = .FALSE.

  call splitn_bdecode1('splitn_lget',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  iblock=max(0,min(nupdate,iblock))

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_lget: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_lget: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'L') then
     write(6,*) '?splitn_lget: not a LOGICAL object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_lget: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
  else
     ii=nlog + (iblock-1)*nlog_st + varlist(jj)%addr_st
  endif

  ierr = 0
  array = logbuf(ii:ii+isize-1)

end subroutine splitn_lget

subroutine splitn_rget(zname,isize,array,ierr)

  !  get an REAL namelist object
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  real, intent(out) :: array(isize)     ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = 0

  call splitn_bdecode1('splitn_rget',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  iblock=max(0,min(nupdate,iblock))

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_rget: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_rget: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'R') then
     write(6,*) '?splitn_rget: not a REAL object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_rget: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
  else
     ii=nreal + (iblock-1)*nreal_st + varlist(jj)%addr_st
  endif

  ierr = 0
  array = rbuf(ii:ii+isize-1)

end subroutine splitn_rget

subroutine splitn_dget(zname,isize,array,ierr)

  !  get an REAL*8 namelist object
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  real*8, intent(out) :: array(isize)   ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = 0

  call splitn_bdecode1('splitn_dget',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  iblock=max(0,min(nupdate,iblock))

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_dget: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_dget: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'D') then
     write(6,*) '?splitn_dget: not a REAL*8 object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_dget: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
  else
     ii=nr8 + (iblock-1)*nr8_st + varlist(jj)%addr_st
  endif

  ierr = 0
  array = dbuf(ii:ii+isize-1)

end subroutine splitn_dget

subroutine splitn_cget(zname,iwid,isize,array,ierr)

  !  get an CHARACTER*n namelist object
  !  isize *must* match the actual size of the object.
  !  iwid must be .ge. the actual width of object elements.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: iwid           ! width of object array elements
  integer, intent(in) :: isize          ! size of object (# of elements)
  character*(iwid), intent(out) :: array(isize)   ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = ' '

  call splitn_bdecode1('splitn_cget',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  iblock=max(0,min(nupdate,iblock))

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_cget: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_cget: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type(1:1).ne.'C') then
     write(6,*) '?splitn_cget: not a CHARACTER*nn object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_cget: array size mismatch: expected=',isize, &
          ' found=',ii
     return
  endif

  if(iwid.lt.varlist(jj)%chsize) then
     write(6,*) '?splitn_cget: passed character element size too small: ', &
          'have=',iwid,' need=',varlist(jj)%chsize
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
  else
     ii=nchv + (iblock-1)*nchv_st + varlist(jj)%addr_st
  endif

  ierr = 0
  do i=1,isize
     array(i) = chbuf(ii+i-1)
  enddo

end subroutine splitn_cget
