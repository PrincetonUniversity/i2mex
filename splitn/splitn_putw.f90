subroutine splitn_iputw(zname,isize,array,ierr)

  !  put an INTEGER namelist object in the write buffer;
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  integer, intent(in) :: array(isize)   ! object data (copied to write buffer)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_iputw',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_iputw: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_iputw: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'I') then
     write(6,*) '?splitn_iputw: not an INTEGER object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_iputw: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))

  if(iblock.eq.0) then
     ii = varlist(jj)%addr
  else
     ii = nint + (iblock-1)*nint_st + varlist(jj)%addr_st
  endif

  ierr = 0
  intbuf_w(ii:ii+isize-1) = array

end subroutine splitn_iputw

subroutine splitn_lputw(zname,isize,array,ierr)

  !  put a LOGICAL namelist object in the write buffer;
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  logical, intent(in) :: array(isize)   ! object data (copied to write buffer)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_lputw',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_lputw: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_lputw: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'L') then
     write(6,*) '?splitn_lputw: not a LOGICAL object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_lputw: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))

  if(iblock.eq.0) then
     ii = varlist(jj)%addr
  else
     ii = nlog + (iblock-1)*nlog_st + varlist(jj)%addr_st
  endif

  ierr = 0
  logbuf_w(ii:ii+isize-1) = array

end subroutine splitn_lputw

subroutine splitn_rputw(zname,isize,array,ierr)

  !  put an REAL namelist object in the write buffer;
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  real, intent(in) :: array(isize)      ! object data (copied to write buffer)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_rputw',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_rputw: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_rputw: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'R') then
     write(6,*) '?splitn_rputw: not a REAL object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_rputw: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))

  if(iblock.eq.0) then
     ii = varlist(jj)%addr
  else
     ii = nreal + (iblock-1)*nreal_st + varlist(jj)%addr_st
  endif

  ierr = 0
  rbuf_w(ii:ii+isize-1) = array

end subroutine splitn_rputw

subroutine splitn_dputw(zname,isize,array,ierr)

  !  put an REAL*8 namelist object in the write buffer;
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  real*8, intent(in) :: array(isize)    ! object data (copied to write buffer)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_dputw',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_notilda('splitn_dputw',znam32,ierr)
  if(ierr.ne.0) return

  ierr = 1  ! assume for now...

  if(.not.have_database) then
     write(6,*) '?splitn_dputw: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_dputw: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'D') then
     write(6,*) '?splitn_dputw: not a REAL*8 object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_dputw: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))

  if(iblock.eq.0) then
     ii = varlist(jj)%addr
  else
     ii = nr8 + (iblock-1)*nr8_st + varlist(jj)%addr_st
  endif

  if(isize.eq.1) then
     if((varlist(jj)%name.eq.'TINIT').and.(array(1).ge.tup_w(1))) then
        write(6,*) '?splitn_dputw: new TINIT value .ge. 1st update time:', &
             array(1),tup_w(1)
        return
     else
        tinit_w=array(1)
     endif
  endif

  ierr = 0
  dbuf_w(ii:ii+isize-1) = array

end subroutine splitn_dputw

subroutine splitn_cputw(zname,iwid,isize,array,ierr)

  !  put a CHARACTER*n namelist object in the write buffer;
  !  isize *must* match the actual size of the object.
  !  iwid must be .le. the actual width of object elements.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: iwid           ! width of object array elements
  integer, intent(in) :: isize          ! size of object (# of elements)
  character*(iwid), intent(in) :: array(isize)   ! object data (copied)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch,iblock,iwuse
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_cputw',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_cputw: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_cputw: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type(1:1).ne.'C') then
     write(6,*) '?splitn_cputw: not a CHARACTER*nn object: ',znam32
     return
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_cputw: array size mismatch: expected=',isize, &
          ' found=',ii
     return
  endif

  iwuse=iwid
  if(iwid.gt.varlist(jj)%chsize) then
     iwuse=1
     do i=1,isize
        iwuse=max(iwuse,len_trim(array(i)))
     enddo
     if(iwuse.gt.varlist(jj)%chsize) then
        write(6,*) &
             '?splitn_cputw: passed character element size too large: ', &
             ' passed=',iwid,' allowed=',varlist(jj)%chsize
        write(6,*) ' max non-blank width = ',iwuse
        write(6,*) ' trailing blanks ignored in size comparison.'
        return
     endif
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))

  if(iblock.eq.0) then
     ii = varlist(jj)%addr
  else
     ii = nchv + (iblock-1)*nchv_st + varlist(jj)%addr_st
  endif

  ierr = 0
  do i=1,isize
     chbuf_w(ii+i-1)=array(i)(1:iwuse)
  enddo

end subroutine splitn_cputw
