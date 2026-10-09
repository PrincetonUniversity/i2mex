!
!  The *_*getf(...) subroutines only return value elements explicitly set
!  in the namelist file.  Elements not explicitly set are left unmodified.
!
!  dmc Jan 2005: added splitn_getf_flag, see below...
!
!  MODIFIED dmc Mar 2009: update blocks: some namelist quantities are now
!    updatable, i.e., can be modified via read of an update namelist, which
!    is trackable in splitn.  Append the update block number in syntax
!    "[<num>]" to have the routine check for updates of the quantity in
!    block #<num>.  Where updates are specified, element values of the
!    passed array are set.  If this suffix is omitted or "[0]" the initial 
!    namelist is referred to, as before.  If a quantity is not updatable 
!    but "[<num>]" is specified with <num> .gt. 0, then, no quantities are set.
!
subroutine splitn_igetf(zname,isize,array,ierr)

  !  get an INTEGER namelist object
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  integer, intent(inout) :: array(isize)  ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ia,ivar,imatch,iblock,ilin
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_igetf',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_igetf: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_igetf: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'I') then
     write(6,*) '?splitn_igetf: not an INTEGER object: ',znam32
     return
  endif

  iblock=max(0,min(nupdate,iblock))
  if((iblock.gt.0).AND.varlist(jj)%steerable.le.1) then
     return  ! no changes possible
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_igetf: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
     ilin=varlist(jj)%nlinadr
  else
     ii = nint + (iblock-1)*nint_st + varlist(jj)%addr_st
     ilin = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + &
          varlist(jj)%nlinadr_st
  endif

  ierr = 0
  do ia=1,isize
     if(ilines(ilin+ia-1).ne.0) then
        array(ia) = intbuf(ii+ia-1)
     endif
  enddo

end subroutine splitn_igetf

subroutine splitn_lgetf(zname,isize,array,ierr)

  !  get a LOGICAL namelist object
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  logical, intent(inout) :: array(isize)  ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ia,ivar,imatch,iblock,ilin
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_lgetf',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_lgetf: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_lgetf: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'L') then
     write(6,*) '?splitn_lgetf: not a LOGICAL object: ',znam32
     return
  endif

  iblock=max(0,min(nupdate,iblock))
  if((iblock.gt.0).AND.varlist(jj)%steerable.le.1) then
     return  ! no changes possible
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_lgetf: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
     ilin=varlist(jj)%nlinadr
  else
     ii = nlog + (iblock-1)*nlog_st + varlist(jj)%addr_st
     ilin = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + &
          varlist(jj)%nlinadr_st
  endif

  ierr = 0
  do ia=1,isize
     if(ilines(ilin+ia-1).ne.0) then
        array(ia) = logbuf(ii+ia-1)
     endif
  enddo

end subroutine splitn_lgetf

subroutine splitn_rgetf(zname,isize,array,ierr)

  !  get an REAL namelist object
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  real, intent(inout) :: array(isize)     ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ia,ivar,imatch,iblock,ilin
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_rgetf',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_rgetf: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_rgetf: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'R') then
     write(6,*) '?splitn_rgetf: not a REAL object: ',znam32
     return
  endif

  iblock=max(0,min(nupdate,iblock))
  if((iblock.gt.0).AND.varlist(jj)%steerable.le.1) then
     return  ! no changes possible
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_rgetf: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
     ilin=varlist(jj)%nlinadr
  else
     ii = nreal + (iblock-1)*nreal_st + varlist(jj)%addr_st
     ilin = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + &
          varlist(jj)%nlinadr_st
  endif

  ierr = 0
  do ia=1,isize
     if(ilines(ilin+ia-1).ne.0) then
        array(ia) = rbuf(ii+ia-1)
     endif
  enddo

end subroutine splitn_rgetf

subroutine splitn_dgetf(zname,isize,array,ierr)

  !  get an REAL*8 namelist object
  !  isize *must* match the actual size of the object.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  real*8, intent(inout) :: array(isize)   ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ia,ivar,imatch,iblock,ilin
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_dgetf',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_dgetf: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_dgetf: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type.ne.'D') then
     write(6,*) '?splitn_dgetf: not a REAL*8 object: ',znam32
     return
  endif

  iblock=max(0,min(nupdate,iblock))
  if((iblock.gt.0).AND.varlist(jj)%steerable.le.1) then
     return  ! no changes possible
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_dgetf: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
     ilin=varlist(jj)%nlinadr
  else
     ii = nr8 + (iblock-1)*nr8_st + varlist(jj)%addr_st
     ilin = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + &
          varlist(jj)%nlinadr_st
  endif

  ierr = 0
  do ia=1,isize
     if(ilines(ilin+ia-1).ne.0) then
        array(ia) = dbuf(ii+ia-1)
     endif
  enddo

end subroutine splitn_dgetf

subroutine splitn_cgetf(zname,iwid,isize,array,ierr)

  !  get an CHARACTER*n namelist object
  !  isize *must* match the actual size of the object.
  !  iwid must be .ge. the actual width of object elements.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: iwid           ! width of object array elements
  integer, intent(in) :: isize          ! size of object (# of elements)
  character*(iwid), intent(inout) :: array(isize)   ! object data (returned)
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ia,ivar,imatch,iblock,ilin
  character*32 znam32
  !---------------------------

  call splitn_bdecode1('splitn_cgetf',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_cgetf: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_cgetf: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  if(varlist(jj)%type(1:1).ne.'C') then
     write(6,*) '?splitn_cgetf: not a CHARACTER*nn object: ',znam32
     return
  endif

  iblock=max(0,min(nupdate,iblock))
  if((iblock.gt.0).AND.varlist(jj)%steerable.le.1) then
     return  ! no changes possible
  endif

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_cgetf: array size mismatch: expected=',isize, &
          ' found=',ii
     return
  endif

  if(iwid.lt.varlist(jj)%chsize) then
     write(6,*) '?splitn_cgetf: passed character element size too small: ', &
          'have=',iwid,' need=',varlist(jj)%chsize
     return
  endif

  if(iblock.eq.0) then
     ii=varlist(jj)%addr
     ilin=varlist(jj)%nlinadr
  else
     ii = nchv + (iblock-1)*nchv_st + varlist(jj)%addr_st
     ilin = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + &
          varlist(jj)%nlinadr_st
  endif

  ierr = 0
  do ia=1,isize
     if(ilines(ilin+ia-1).ne.0) then
        array(ia) = chbuf(ii+ia-1)
     endif
  enddo

end subroutine splitn_cgetf

subroutine splitn_getf_flag(zname,isize,array,ierr)

  !  set each element of a logical array according as a quantity
  !  is defined in the file or from the defaults.  The quantity can be
  !  of any data type.
  ! 
  !  on output, array(j)=.TRUE. means the corresponding namelist element
  !  is defined in the namelist file; .FALSE. means the definition comes
  !  from the namelist default settings.  If the quantity is explicitly
  !  in the namelist file, .TRUE. is returned for corresponding array
  !  elements, even if the explictly assigned value matches the default
  !  value.
  !
  !  the information pertains to the initial namelist settings; update
  !  blocks are not taken into account.  See splitn_getu_flag(...)

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  logical, intent(out) :: array(isize)  ! file definition flags
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ia,ivar,imatch,iblock
  character*32 znam32
  !---------------------------

  array = .FALSE.

  call splitn_bdecode1('splitn_getf_flag',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_getf_flag: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_getf_flag: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_getf_flag: size mismatch: expected=',isize,' found=',ii
     return
  endif

  ierr = 0
  do ia=1,isize
     if(ilines(varlist(jj)%nlinadr+ia-1).ne.0) then
        array(ia) = .TRUE.
     endif
  enddo

end subroutine splitn_getf_flag

subroutine splitn_getu_flag(zname,isize,array,ierr)

  !  set each element of a logical array according as a quantity
  !  is specified in a namelist file update blocks or in any of a set
  !  of blocks.  The quantity can be of any data type.
  ! 
  !  on output, array(j)=.TRUE. means the corresponding namelist element
  !  is set (even if the value is not actually changed); .FALSE. means it 
  !  does not appear in the specified update block(s).

  !  The "zname" argument controls which update blocks to check:

  !  if the bare name is specified, or if it is appended with "[*]" or "[:]"
  !  or "[0]", all update blocks are checked.  If "[2]" is appended, only 
  !  update block #2 is checked; if "[2:4]" is appended, 2,3, and 4 are 
  !  checked.  (see subroutine splitn_bdecode2).

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(in) :: isize          ! size of object
  logical, intent(out) :: array(isize)  ! file definition flags
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ia,ivar,imatch,iblock1,iblock2,iblock,ilin
  character*32 znam32
  !---------------------------

  array = .FALSE.

  call splitn_bdecode2('splitn_getf_flag',zname,znam32,iblock1,iblock2,ierr)
  if(ierr.ne.0) return

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?splitn_getf_flag: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_getf_flag: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  if(ii.ne.isize) then
     write(6,*) '?splitn_getf_flag: size mismatch: expected=',isize,' found=',ii
     return
  endif

  if(nupdate.eq.0) then
     return  ! there are no updates in the file.
  endif

  if(varlist(jj)%steerable.le.1) then
     return  ! the quantity is not updatable: so, no updates.
  endif

  if((iblock1.eq.0).and.(iblock2.eq.0)) then
     iblock1=1
     iblock2=nupdate
  endif

  if(iblock1.eq.0) iblock1 = 1

  ierr = 0

  do iblock = iblock1,iblock2
     ilin = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + &
          varlist(jj)%nlinadr_st
     do ia=1,isize
        if(ilines(ilin+ia-1).ne.0) then
           array(ia) = .TRUE.
        endif
     enddo
  enddo

end subroutine splitn_getu_flag
