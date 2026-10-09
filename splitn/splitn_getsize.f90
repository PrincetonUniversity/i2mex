subroutine splitn_getsize(zname,isize,ierr)

  !  get the size (number of array elements or 1 if scalar) of a namelist 
  !  object.  Information on multiple dimension not returned; total size only.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(out) :: isize         ! size of object
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch
  character*32 znam32
  !---------------------------

  ierr = 1

  znam32=zname
  call uupper(znam32)

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

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo

  isize=ii
  ierr = 0

end subroutine splitn_getsize

subroutine splitn_getdims(zname,irank,irank_max,dims,isize,ierr)

  !  get the size (number of array elements or 1 if scalar) of a namelist 
  !  object along with array dimension information

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer, intent(out) :: irank         ! rank of object (0=scalar, 1=vector..)
  integer, intent(in) :: irank_max      ! max rank (dim of dim array)
  integer, intent(out) :: dims(2,irank_max) ! dimensioning info
  integer, intent(out) :: isize         ! total size of object
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,ii,jj,ivar,imatch
  character*32 znam32
  !---------------------------

  ierr = 1

  znam32=zname
  call uupper(znam32)

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

  irank = varlist(jj)%rank
  dims = 0
  isize = 0

  if(irank.gt.irank_max) then
     write(6,*) '?splitn_getdims: max rank exceeded: ',irank_max,irank
  else
     ii=1
     do i=1,irank
        dims(1:2,i)=varlist(jj)%dims(1:2,i)
        ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
     enddo
  endif

  isize=ii
  ierr = 0

end subroutine splitn_getdims

subroutine splitn_get_type(zname,ztype,ichsize,ierr)

  !  get the datatype (real/integer/logical/double/character) of a namelist 
  !  object.  If type Character*nnn return nnn also.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  character*(*), intent(out) :: ztype   ! (CHARACTER*5) datatype code:
  !  R/I/L/D/C*nnn
  integer, intent(out) :: ichsize       ! if C*n -- number of characters n
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,jj,ivar,imatch
  character*32 znam32
  !---------------------------

  ierr = 1

  znam32=zname
  call uupper(znam32)

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

  ztype = varlist(jj)%type
  ichsize = 0
  if(ztype(1:2).eq.'C*') ichsize = varlist(jj)%chsize
  ierr = 0

end subroutine splitn_get_type

subroutine splitn_get_dstr(zname,dstr,ierr)

  !  get the string which defines the object's default value.
  !  if the passed string is too short a truncated string is returned
  !  without error or warning indicated.  SPLITN default strings are
  !  never more than 512 characters long.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  character*(*), intent(out) :: dstr    ! default specifier
  integer, intent(out) :: ierr          ! completion code (0=OK)

  !---------------------------
  integer i,jj,ivar,imatch,ilen,jlen,klen,iadr
  character*32 znam32
  !---------------------------

  ierr = 1

  znam32=zname
  call uupper(znam32)

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

  dstr = ' '
  iadr=varlist(jj)%long_dflt_addr
  if(iadr.eq.0) then
     ilen=len_trim(varlist(jj)%short_dflt)
     jlen=len(dstr)
     klen=min(ilen,jlen)
     dstr(1:klen)=varlist(jj)%short_dflt(1:klen)
  else
     ilen=len_trim(long_dflts(iadr))
     jlen=len(dstr)
     klen=min(ilen,jlen)
     dstr(1:klen)=long_dflts(iadr)(1:klen)
  endif
  ierr = 0

end subroutine splitn_get_dstr

subroutine splitn_get_trustat(zname,iupdate,ierr)

  ! return iupdate=1 if quantity (zname) is TRANSP-namelist-updatable.
  ! otherwise return iupdate=0.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer,intent(out) :: iupdate        ! =0: not updatable;
                                   ! =1: updatable via TRANSP update namelist

  integer,intent(out) :: ierr           ! status return code 0=OK

  !-----------------------------------
  character*32 :: znam32
  integer :: indbr1,ivar,imatch,jj
  !-----------------------------------

  ierr = 0

  indbr1 = index(zname,'[')
  if(indbr1.gt.1) then
     znam32 = zname(1:indbr1-1)
  else
     znam32 = zname
  endif
  call uupper(znam32)
 
  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     ierr=1
     write(6,*) '?splitn_get_trustat: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)
 
  iupdate = 0

  if(varlist(jj)%steerable.eq.2) then
     iupdate = 1
  endif

end subroutine splitn_get_trustat

subroutine splitn_get_udetails(zname,iupdatable, &
     idefault,inum_update,inmax,iupdate_vector,ierr)

  !  get default/update status information on item:
  !    (a) whether it is updatable (trdat only or trdat & TRANSP);
  !    (b) whether its (initial) value is defaulted
  !    (c) whether value change updates exist after initial setting
  !    (d) vector indicating which blocks contain updates

  !  for arrays-- a non-default initial value for any array element
  !    means, the whole array is considered to have a non-default value;
  !    similarly, if any array element is changed on update, the whole
  !    array is considered to have been changed.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname    ! name of object
  integer,intent(out) :: iupdatable     ! =0: not updatable; =1: trdat only
                                        ! >1: updatable in TRANSP

  integer,intent(out) :: idefault       ! =-1: no text, default initial value
                                        ! =0: text present, still the default
                                        ! =1: non-default value
                                        !     (this implies text is present)

  integer,intent(out) :: inum_update    ! number of update blocks w/changes

  integer,intent(in) :: inmax           ! array size for iupdate_vector

  integer,intent(out) :: iupdate_vector(inmax) ! update status vector:
                                        ! iupdate_vector(j) for update block #j
                                        !   -1 means: no text, no value change
                                        !    0 means: text but no value change
                                        !    1 means: value change
                                        !      (this implies text is present).

  integer,intent(out) :: ierr           ! status return code 0=OK

  !---------------------------
  integer i,ii,jj,ivar,imatch,ilen,isize,ichsize,itext,iup,iadl0,iadl1,iadl2,ic
  character*32 znam32
  character*40 zcombo
  character*5 :: ztype
  character*1 :: dtype
  character*6 :: zi6
  character*8 :: zupid

  real, dimension(:), allocatable :: r0,r1
  real*8, dimension(:), allocatable :: d0,d1
  integer, dimension(:), allocatable :: i0,i1
  logical, dimension(:), allocatable :: L0,L1
  integer, parameter :: iwid=128
  character*(iwid), dimension(:), allocatable :: ch0,ch1

  !---------------------------
  ! default outputs:

  ierr = 1

  iupdatable = 0
  idefault = -1
  inum_update = 0
  iupdate_vector = -1

  !---------------------------
  ! check name
  znam32=zname
  call uupper(znam32)

  if(.not.have_database) then
     write(6,*) '?splitn_get_udetails: no namelist has been read.'
     return
  endif

  call iorder(znam32,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?splitn_get_udetails: unrecognized name: ',znam32
     return
  endif

  jj = var_order(ivar)
  
  !---------------------------
  !  Name is OK

  iupdatable = varlist(jj)%steerable

  if((iupdatable.ge.2).and.(nupdate.gt.inmax)) then
     write(6,*) '?splitn_get_udetails: update status vector size too small:'
     write(6,*) ' received: ',inmax,'; need: ',nupdate,'.'
     return
  endif

  !---------------------------
  !  Error checks OK

  ierr=0

  !  get size for internal use

  ii=1
  do i=1,varlist(jj)%rank
     ii=ii*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
  enddo
  isize=ii

  !  look for text.  If no text, then, default values apply and there
  !  are no updates and this is known without further checking...

  itext = 0
  do ii=1,isize
     if(ilines(varlist(jj)%nlinadr+ii-1).ne.0) then
        itext = 1
        idefault = 0   ! text detected
        exit
     endif
  enddo

  iadl0 = nreal + nr8 + nint + nlog + nchv
  iadl1 = nreal_st + nr8_st + nint_st + nlog_st + nchv_st

  if((itext.eq.0).and.(iupdatable.ge.2)) then
     ! look for update text...
     do iup=1,nupdate
        do ii=1,isize
           iadl2 = iadl0 + (iup-1)*iadl1 + varlist(jj)%nlinadr_st + ii - 1
           if(ilines(iadl2).ne.0) then
              itext=1
              iupdate_vector(iup)=0  ! text detected
              exit
           endif
        enddo
     enddo
  endif

  if(itext.eq.0) then
     ! the quantity is never referenced in the namelist text; therefore,
     ! its default value applies throughout -- we are done here.
     return
  endif

  ztype = varlist(jj)%type
  ichsize = 0
  if(ztype(1:2).eq.'C*') ichsize = varlist(jj)%chsize

  !-----------------------------
  ! OK, text exists: scan through values to determine changes

  dtype = ztype(1:1)

  if(dtype.eq.'R') then
     allocate(r0(isize),r1(isize))
  else if(dtype.eq.'I') then
     allocate(i0(isize),i1(isize))
  else if(dtype.eq.'D') then
     allocate(d0(isize),d1(isize))
  else if(dtype.eq.'L') then
     allocate(L0(isize),L1(isize))
  else if(dtype.eq.'C') then
     allocate(ch0(isize),ch1(isize))
  endif

  do
     ! check defaults...
     if(dtype.eq.'R') then
        call splitn_rgetd(zname,isize,r0,ierr); if(ierr.ne.0) exit
        r1=r0
        call splitn_rget(zname,isize,r1,ierr); if(ierr.ne.0) exit
        do ii=1,isize
           if(r0(ii).ne.r1(ii)) then
              idefault=1  ! value change detected
              exit
           endif
        enddo

     else if(dtype.eq.'I') then
        call splitn_igetd(zname,isize,i0,ierr); if(ierr.ne.0) exit
        i1=i0
        call splitn_iget(zname,isize,i1,ierr); if(ierr.ne.0) exit
        do ii=1,isize
           if(i0(ii).ne.i1(ii)) then
              idefault=1  ! value change detected
              exit
           endif
        enddo

     else if(dtype.eq.'D') then
        call splitn_dgetd(zname,isize,d0,ierr); if(ierr.ne.0) exit
        d1=d0
        call splitn_dget(zname,isize,d1,ierr); if(ierr.ne.0) exit
        do ii=1,isize
           if(d0(ii).ne.d1(ii)) then
              idefault=1  ! value change detected
              exit
           endif
        enddo

     else if(dtype.eq.'L') then
        call splitn_Lgetd(zname,isize,L0,ierr); if(ierr.ne.0) exit
        L1=L0
        call splitn_Lget(zname,isize,L1,ierr); if(ierr.ne.0) exit
        do ii=1,isize
           if((L0(ii).AND.(.not.L1(ii))) .OR. ((.not.L0(ii)).AND.L1(ii))) then
              idefault=1  ! value change detected
              exit
           endif
        enddo

     else if(dtype.eq.'C') then
        call splitn_cgetd(zname,iwid,isize,ch0,ierr); if(ierr.ne.0) exit
        ch1=ch0
        call splitn_cget(zname,iwid,isize,ch1,ierr); if(ierr.ne.0) exit
        do ii=1,isize
           if(ch0(ii).ne.ch1(ii)) then
              idefault=1  ! value change detected
              exit
           endif
        enddo

     endif

     ! check update blocks...
     if(iupdatable.le.1) exit

     do iup=1,nupdate
        
        ! update reference string
        zi6=' '
        write(zi6,'(i6)') iup
        ic=1
        do ii=5,1,-1
           if(zi6(ii:ii).eq.' ') then
              ic=ii+1
              exit
           endif
        enddo

        zupid = '['//zi6(ic:6)//']'
        zcombo = trim(zname)//trim(zupid)

        if(dtype.eq.'R') then

           r0=r1
           call splitn_rget(zcombo,isize,r1,ierr); if(ierr.ne.0) exit
           do ii=1,isize
              if(r0(ii).ne.r1(ii)) then
                 inum_update = inum_update + 1
                 iupdate_vector(iup)=1  ! value change detected
                 exit
              endif
           enddo

        else if(dtype.eq.'I') then

           i0=i1
           call splitn_iget(zcombo,isize,i1,ierr); if(ierr.ne.0) exit
           do ii=1,isize
              if(i0(ii).ne.i1(ii)) then
                 inum_update = inum_update + 1
                 iupdate_vector(iup)=1  ! value change detected
                 exit
              endif
           enddo

        else if(dtype.eq.'D') then

           d0=d1
           call splitn_dget(zcombo,isize,d1,ierr); if(ierr.ne.0) exit
           do ii=1,isize
              if(d0(ii).ne.d1(ii)) then
                 inum_update = inum_update + 1
                 iupdate_vector(iup)=1  ! value change detected
                 exit
              endif
           enddo

        else if(dtype.eq.'L') then

           L0=L1
           call splitn_lget(zcombo,isize,L1,ierr); if(ierr.ne.0) exit
           do ii=1,isize
              if((L0(ii).AND.(.not.L1(ii))) .OR. &
                   ((.not.L0(ii)).AND.L1(ii))) then
                 inum_update = inum_update + 1
                 iupdate_vector(iup)=1  ! value change detected
                 exit
              endif
           enddo

        else if(dtype.eq.'C') then

           ch0=ch1
           call splitn_cget(zcombo,isize,ch1,ierr); if(ierr.ne.0) exit
           do ii=1,isize
              if(ch0(ii).ne.ch1(ii)) then
                 inum_update = inum_update + 1
                 iupdate_vector(iup)=1  ! value change detected
                 exit
              endif
           enddo

        endif
     enddo

     exit
  enddo

  if(dtype.eq.'R') then
     deallocate(r0,r1)
  else if(dtype.eq.'I') then
     deallocate(i0,i1)
  else if(dtype.eq.'D') then
     deallocate(d0,d1)
  else if(dtype.eq.'L') then
     deallocate(L0,L1)
  else if(dtype.eq.'C') then
     deallocate(ch0,ch1)
  endif

end subroutine splitn_get_udetails

subroutine splitn_update_report

  ! report updated quantities -- #of updates

  use splitn_module
  implicit NONE

  integer :: ii,jj,kk,iany

  integer :: inmax,idefault,iupdatable,inum_update,ierr,isum
  integer, dimension(:), allocatable :: iupdate_vector
  !---------------------------------------

  inmax=max(1,nupdate)
  allocate(iupdate_vector(inmax))

  iany = 0
  do ii=1,nvars
     jj=var_order(ii)
     if(varlist(jj)%name.eq.'~UPDATE_TIME') cycle

     call splitn_get_udetails(varlist(jj)%name, iupdatable, &
          idefault, inum_update, inmax, iupdate_vector, ierr)

     isum=0
     do kk=1,inmax
        if(iupdate_vector(kk).gt.-1) isum=isum+1
     enddo

     if(isum.gt.0) then
        iany = iany + 1
        if(iany.eq.1) then
           write(6,*) ' '
           write(6,*) ' %splitn_update_report: update namelists report:'
           write(6,*) '  (some references are not updates-- i.e. if value matches prior setting).'
           write(6,*) ' '
        endif
        write(6,*) varlist(jj)%name,' ...has ',isum,' reference(s) and ',inum_update,' update(s).'
     endif
  enddo

  if(iany.gt.0) write(6,*) ' '

end subroutine splitn_update_report
