subroutine splitn_igetd(zname,isize,array,ierr)

  !  get an INTEGER namelist object *default values*
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object ("[<num>]" ignored).
  integer, intent(in) :: isize          ! size of object
  integer, intent(out) :: array(isize)  ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = 0

  call splitn_bdecode1('splitn_igetd',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_igetd: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_igetd: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'I') then
     write(6,*) '?splitn_igetd: not an INTEGER object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_igetd: size mismatch: expected=',isize,' found=',ii
     return
  endif

  ii=varlist(jj)%addr

  ierr = 0
  array = intbuf_d(ii:ii+isize-1)

end subroutine splitn_igetd

subroutine splitn_lgetd(zname,isize,array,ierr)

  !  get a LOGICAL namelist object *default values*
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object ("[<num>]" ignored).
  integer, intent(in) :: isize          ! size of object
  logical, intent(out) :: array(isize)  ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = .FALSE.

  call splitn_bdecode1('splitn_lgetd',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_lgetd: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_lgetd: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'L') then
     write(6,*) '?splitn_lgetd: not a LOGICAL object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_lgetd: size mismatch: expected=',isize,' found=',ii
     return
  endif

  ii=varlist(jj)%addr

  ierr = 0
  array = logbuf_d(ii:ii+isize-1)

end subroutine splitn_lgetd

subroutine splitn_rgetd(zname,isize,array,ierr)

  !  get an REAL namelist object *default values*
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object ("[<num>]" ignored).
  integer, intent(in) :: isize          ! size of object
  real, intent(out) :: array(isize)     ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = 0

  call splitn_bdecode1('splitn_rgetd',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_rgetd: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_rgetd: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'R') then
     write(6,*) '?splitn_rgetd: not a REAL object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_rgetd: size mismatch: expected=',isize,' found=',ii
     return
  endif

  ii=varlist(jj)%addr

  ierr = 0
  array = rbuf_d(ii:ii+isize-1)

end subroutine splitn_rgetd

subroutine splitn_dgetd(zname,isize,array,ierr)

  !  get an REAL*8 namelist object *default values*
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object ("[<num>]" ignored).
  integer, intent(in) :: isize          ! size of object
  real*8, intent(out) :: array(isize)   ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = 0

  call splitn_bdecode1('splitn_dgetd',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_dgetd: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_dgetd: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'D') then
     write(6,*) '?splitn_dgetd: not a REAL*8 object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_dgetd: size mismatch: expected=',isize,' found=',ii
     return
  endif

  ii=varlist(jj)%addr

  ierr = 0
  array = dbuf_d(ii:ii+isize-1)

end subroutine splitn_dgetd

subroutine splitn_cgetd(zname,iwid,isize,array,ierr)

  !  get an CHARACTER*n namelist object  *default values*
  !  isize *must* match the actual size of the object.
  !  iwid must be .ge. the actual width of object elements.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object ("[<num>]" ignored).
  integer, intent(in) :: iwid           ! width of object array elements
  integer, intent(in) :: isize          ! size of object (# of elements)
  character*(iwid), intent(out) :: array(isize)   ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = ' '

  call splitn_bdecode1('splitn_cgetd',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_cgetd: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_cgetd: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type(1:1).ne.'C') then
     write(6,*) '?splitn_cgetd: not a CHARACTER*nn object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_cgetd: array size mismatch: expected=',isize, &
          ' found=',ii
     return
  endif

  if(iwid.lt.varlist(jj)%chsize) then
     write(6,*) '?splitn_cgetd: passed character element size too small: ', &
          'have=',iwid,' need=',varlist(jj)%chsize
     return
  endif

  ii=varlist(jj)%addr

  ierr = 0
  do i=1,isize
     array(i) = chbuf_d(ii+i-1)
  enddo

end subroutine splitn_cgetd
