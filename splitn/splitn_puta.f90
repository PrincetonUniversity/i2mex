subroutine splitn_iput_ar(zname,ivalue,irank,ind0,inum,ierr)

  ! assign new value to elements of an INTEGER array

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of array variable
  integer, intent(in) :: ivalue       ! new value to assign
  integer, intent(in) :: irank        ! rank of array (must match)
  integer, intent(in) :: ind0(irank)  ! start index for assignment
  integer, intent(in) :: inum         ! number of elements to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ioff,isize,iblock,ilinadr,iadst,ilinst
  character*32 znam32,istr
  !--------------------------------------

  call splitn_bdecode1('splitn_iput_ar',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_iput_ar',znam32,jj,ierr)
  if(ierr.ne.0) return

  if(varlist(jj)%type.ne.'I') then
     write(6,*) '?splitn_iput_ar: not INTEGER: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.eq.0) then
     write(6,*) '?splitn_iput_ar: use splitn_iput_sc to assign SCALAR.'
     ierr=1
  endif

  if(varlist(jj)%rank.ne.irank) then
     write(6,*) '?splitn_iput_ar: rank mismatch: ',trim(znam32)
     write(6,*) ' object rank = ',varlist(jj)%rank,'; passed rank = ',irank
     ierr=1
  endif

  if(inum.le.0) then
     write(6,*) '?splitn_iput_ar: inum=',inum
     write(6,*) ' number of elements to copy (inum) must be a positive number.'
     ierr=1
  endif

  if(ierr.ne.0) return

  do ii=1,irank
     if((ind0(ii).lt.varlist(jj)%dims(1,ii)).or. &
          (ind0(ii).gt.varlist(jj)%dims(2,ii)) ) then
        ierr=1
        write(6,*) '?splitn_iput_ar: subscript bounds on dimension #',ii
        write(6,*) ' passed value = ',ind0(ii),' not in range ', &
             varlist(jj)%dims(1,ii),':',varlist(jj)%dims(2,ii)
     endif
  enddo
  if(ierr.ne.0) return

  ioff=0
  isize=1
  do ii=irank,1,-1
     ioff = ioff + (ind0(ii)-varlist(jj)%dims(1,ii))*isize
     isize = isize*(varlist(jj)%dims(2,ii)-varlist(jj)%dims(1,ii)+1)
  enddo

  if(ioff+inum.gt.isize) then
     write(6,*) '?splitn_iput_ar: ind0 = (',ind0,'); inum=',inum
     write(6,*) ' write reference beyond end of array.'
     ierr=1
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock
 
  ilinst = varlist(jj)%nlinadr_st
  iadst  = varlist(jj)%addr_st
  if(iblock.eq.0) then
     ii = varlist(jj)%addr + ioff
     ilinadr = varlist(jj)%nlinadr
  else
     ii = nint + (iblock-1)*nint_st + iadst + ioff
     ilinadr = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
  endif

  intbuf(ii:ii+inum-1)=ivalue
  if(varlist(jj)%steerable.eq.2) then
     call iset_ucop(intbuf,ii,inum,iadst+ioff,ilinst+ioff,kupdate,nupdate_max)
  endif

  istr=' '
  write(istr,'(I16)') ivalue
  call splitn_put_ar_str(jj,isize,ioff,inum,znam32,istr,ilinadr)

end subroutine splitn_iput_ar

subroutine splitn_lput_ar(zname,lvalue,irank,ind0,inum,ierr)

  ! assign new value to elements of an LOGICAL array

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of array variable
  logical, intent(in) :: lvalue       ! new value to assign
  integer, intent(in) :: irank        ! rank of array (must match)
  integer, intent(in) :: ind0(irank)  ! start index for assignment
  integer, intent(in) :: inum         ! number of elements to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ioff,isize,iblock,ilinadr,iadst,ilinst
  character*32 znam32
  character*2 lstr
  !--------------------------------------

  call splitn_bdecode1('splitn_lput_ar',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_lput_ar',znam32,jj,ierr)
  if(ierr.ne.0) return

  if(varlist(jj)%type.ne.'L') then
     write(6,*) '?splitn_lput_ar: not LOGICAL: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.eq.0) then
     write(6,*) '?splitn_lput_ar: use splitn_lput_sc to assign SCALAR.'
     ierr=1
  endif

  if(varlist(jj)%rank.ne.irank) then
     write(6,*) '?splitn_lput_ar: rank mismatch: ',trim(znam32)
     write(6,*) ' object rank = ',varlist(jj)%rank,'; passed rank = ',irank
     ierr=1
  endif

  if(inum.le.0) then
     write(6,*) '?splitn_lput_ar: inum=',inum
     write(6,*) ' number of elements to copy (inum) must be a positive number.'
     ierr=1
  endif

  if(ierr.ne.0) return

  do ii=1,irank
     if((ind0(ii).lt.varlist(jj)%dims(1,ii)).or. &
          (ind0(ii).gt.varlist(jj)%dims(2,ii)) ) then
        ierr=1
        write(6,*) '?splitn_lput_ar: subscript bounds on dimension #',ii
        write(6,*) ' passed value = ',ind0(ii),' not in range ', &
             varlist(jj)%dims(1,ii),':',varlist(jj)%dims(2,ii)
     endif
  enddo
  if(ierr.ne.0) return

  ioff=0
  isize=1
  do ii=irank,1,-1
     ioff = ioff + (ind0(ii)-varlist(jj)%dims(1,ii))*isize
     isize = isize*(varlist(jj)%dims(2,ii)-varlist(jj)%dims(1,ii)+1)
  enddo

  if(ioff+inum.gt.isize) then
     write(6,*) '?splitn_lput_ar: ind0 = (',ind0,'); inum=',inum
     write(6,*) ' write reference beyond end of array.'
     ierr=1
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock
 
  ilinst = varlist(jj)%nlinadr_st
  iadst  = varlist(jj)%addr_st
  if(iblock.eq.0) then
     ii = varlist(jj)%addr + ioff
     ilinadr = varlist(jj)%nlinadr
  else
     ii = nlog + (iblock-1)*nlog_st + iadst + ioff
     ilinadr = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
  endif

  logbuf(ii:ii+inum-1)=lvalue
  if(varlist(jj)%steerable.eq.2) then
     call lset_ucop(logbuf,ii,inum,iadst+ioff,ilinst+ioff,kupdate,nupdate_max)
  endif

  if(lvalue) then
     lstr='.T'
  else
     lstr='.F'
  endif
  call splitn_put_ar_str(jj,isize,ioff,inum,znam32,lstr,ilinadr)

end subroutine splitn_lput_ar

subroutine splitn_rput_ar(zname,rvalue,irank,ind0,inum,ierr)

  ! assign new REAL value to elements of a REAL or REAL*8 array

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of array variable
  real, intent(in) :: rvalue          ! new value to assign
  integer, intent(in) :: irank        ! rank of array (must match)
  integer, intent(in) :: ind0(irank)  ! start index for assignment
  integer, intent(in) :: inum         ! number of elements to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ioff,isize,iblock,ilinadr,iadst,ilinst
  character*32 znam32
  character*20 rstr
  character*1 echar
  !--------------------------------------

  call splitn_bdecode1('splitn_rput_ar',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_rput_ar',znam32,jj,ierr)
  if(ierr.ne.0) return

  if((varlist(jj)%type.ne.'R').and.(varlist(jj)%type.ne.'D')) then
     write(6,*) '?splitn_rput_ar: not FLOATING POINT: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.eq.0) then
     write(6,*) '?splitn_rput_ar: use splitn_rput_sc to assign SCALAR.'
     ierr=1
  endif

  if(varlist(jj)%rank.ne.irank) then
     write(6,*) '?splitn_rput_ar: rank mismatch: ',trim(znam32)
     write(6,*) ' object rank = ',varlist(jj)%rank,'; passed rank = ',irank
     ierr=1
  endif

  if(inum.le.0) then
     write(6,*) '?splitn_rput_ar: inum=',inum
     write(6,*) ' number of elements to copy (inum) must be a positive number.'
     ierr=1
  endif

  if(ierr.ne.0) return

  do ii=1,irank
     if((ind0(ii).lt.varlist(jj)%dims(1,ii)).or. &
          (ind0(ii).gt.varlist(jj)%dims(2,ii)) ) then
        ierr=1
        write(6,*) '?splitn_rput_ar: subscript bounds on dimension #',ii
        write(6,*) ' passed value = ',ind0(ii),' not in range ', &
             varlist(jj)%dims(1,ii),':',varlist(jj)%dims(2,ii)
     endif
  enddo
  if(ierr.ne.0) return

  ioff=0
  isize=1
  do ii=irank,1,-1
     ioff = ioff + (ind0(ii)-varlist(jj)%dims(1,ii))*isize
     isize = isize*(varlist(jj)%dims(2,ii)-varlist(jj)%dims(1,ii)+1)
  enddo

  if(ioff+inum.gt.isize) then
     write(6,*) '?splitn_rput_ar: ind0 = (',ind0,'); inum=',inum
     write(6,*) ' write reference beyond end of array.'
     ierr=1
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock
 
  ilinst = varlist(jj)%nlinadr_st
  iadst  = varlist(jj)%addr_st
  if(iblock.eq.0) then
     ii = varlist(jj)%addr + ioff
     ilinadr = varlist(jj)%nlinadr
  else
     if(varlist(jj)%type.eq.'R') then
        ii = nreal + (iblock-1)*nreal_st + iadst + ioff
     else
        ii = nr8 + (iblock-1)*nr8_st + iadst + ioff
     endif
     ilinadr = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
  endif

  if(varlist(jj)%type.eq.'R') then
     rbuf(ii:ii+inum-1)=rvalue
     if(varlist(jj)%steerable.eq.2) then
        call rset_ucop(rbuf,ii,inum,iadst+ioff,ilinst+ioff, &
             kupdate,nupdate_max)
     endif
     echar='e'
  else
     dbuf(ii:ii+inum-1)=rvalue
     if(varlist(jj)%steerable.eq.2) then
        call dset_ucop(dbuf,ii,inum,iadst+ioff,ilinst+ioff, &
             kupdate,nupdate_max)
     endif
     echar='d'
  endif

  rstr=' '
  write(rstr,'(1pe13.6)') rvalue
  call splitn_fput_clean(rstr,echar)
  call splitn_put_ar_str(jj,isize,ioff,inum,znam32,rstr,ilinadr)

end subroutine splitn_rput_ar

subroutine splitn_dput_ar(zname,dvalue,irank,ind0,inum,ierr)

  ! assign new REAL*8 value to elements of a REAL or REAL*8 array

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of array variable
  real*8, intent(in) :: dvalue        ! new value to assign
  integer, intent(in) :: irank        ! rank of array (must match)
  integer, intent(in) :: ind0(irank)  ! start index for assignment
  integer, intent(in) :: inum         ! number of elements to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ioff,isize,iblock,ilinadr,iadst,ilinst
  character*32 znam32
  character*20 dstr
  character*1 echar
  !--------------------------------------

  call splitn_bdecode1('splitn_dput_ar',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_dput_ar',znam32,jj,ierr)
  if(ierr.ne.0) return

  if((varlist(jj)%type.ne.'R').and.(varlist(jj)%type.ne.'D')) then
     write(6,*) '?splitn_dput_ar: not FLOATING POINT: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.eq.0) then
     write(6,*) '?splitn_dput_ar: use splitn_dput_sc to assign SCALAR.'
     ierr=1
  endif

  if(varlist(jj)%rank.ne.irank) then
     write(6,*) '?splitn_dput_ar: rank mismatch: ',trim(znam32)
     write(6,*) ' object rank = ',varlist(jj)%rank,'; passed rank = ',irank
     ierr=1
  endif

  if(inum.le.0) then
     write(6,*) '?splitn_dput_ar: inum=',inum
     write(6,*) ' number of elements to copy (inum) must be a positive number.'
     ierr=1
  endif

  if(ierr.ne.0) return

  do ii=1,irank
     if((ind0(ii).lt.varlist(jj)%dims(1,ii)).or. &
          (ind0(ii).gt.varlist(jj)%dims(2,ii)) ) then
        ierr=1
        write(6,*) '?splitn_dput_ar: subscript bounds on dimension #',ii
        write(6,*) ' passed value = ',ind0(ii),' not in range ', &
             varlist(jj)%dims(1,ii),':',varlist(jj)%dims(2,ii)
     endif
  enddo
  if(ierr.ne.0) return

  ioff=0
  isize=1
  do ii=irank,1,-1
     ioff = ioff + (ind0(ii)-varlist(jj)%dims(1,ii))*isize
     isize = isize*(varlist(jj)%dims(2,ii)-varlist(jj)%dims(1,ii)+1)
  enddo

  if(ioff+inum.gt.isize) then
     write(6,*) '?splitn_dput_ar: ind0 = (',ind0,'); inum=',inum
     write(6,*) ' write reference beyond end of array.'
     ierr=1
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock
 
  ilinst = varlist(jj)%nlinadr_st
  iadst  = varlist(jj)%addr_st
  if(iblock.eq.0) then
     ii = varlist(jj)%addr + ioff
     ilinadr = varlist(jj)%nlinadr
  else
     if(varlist(jj)%type.eq.'R') then
        ii = nreal + (iblock-1)*nreal_st + iadst + ioff
     else
        ii = nr8 + (iblock-1)*nr8_st + iadst + ioff
     endif
     ilinadr = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
  endif

  if(varlist(jj)%type.eq.'R') then
     rbuf(ii:ii+inum-1)=dvalue
     if(varlist(jj)%steerable.eq.2) then
        call rset_ucop(rbuf,ii,inum,iadst+ioff,ilinst+ioff, &
             kupdate,nupdate_max)
     endif
     echar='e'
  else
     dbuf(ii:ii+inum-1)=dvalue
     if(varlist(jj)%steerable.eq.2) then
        call dset_ucop(dbuf,ii,inum,iadst+ioff,ilinst+ioff, &
             kupdate,nupdate_max)
     endif
     echar='d'
  endif

  dstr=' '
  write(dstr,'(1pd19.12)') dvalue
  call splitn_fput_clean(dstr,echar)
  call splitn_put_ar_str(jj,isize,ioff,inum,znam32,dstr,ilinadr)

end subroutine splitn_dput_ar

subroutine splitn_chput_ar(zname,chvalue,irank,ind0,inum,ierr)

  ! assign new value to elements of an CHARACTER*n array

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of array variable
  character*(*), intent(in) :: chvalue ! new value to assign
  integer, intent(in) :: irank        ! rank of array (must match)
  integer, intent(in) :: ind0(irank)  ! start index for assignment
  integer, intent(in) :: inum         ! number of elements to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ioff,isize,iblock,ilinadr,iadst,ilinst
  integer :: ilen,ic,it
  character*32 znam32
  character*1 zdelim
  !--------------------------------------

  call splitn_bdecode1('splitn_chput_ar',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_chput_ar',znam32,jj,ierr)
  if(ierr.ne.0) return

  if(varlist(jj)%type(1:1).ne.'C') then
     write(6,*) '?splitn_chput_ar: not CHARACTER: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.eq.0) then
     write(6,*) '?splitn_chput_ar: use splitn_chput_sc to assign SCALAR.'
     ierr=1
  endif

  if(varlist(jj)%rank.ne.irank) then
     write(6,*) '?splitn_chput_ar: rank mismatch: ',trim(znam32)
     write(6,*) ' object rank = ',varlist(jj)%rank,'; passed rank = ',irank
     ierr=1
  endif

  if(len_trim(chvalue).gt.varlist(jj)%chsize) then
     write(6,*) '?splitn_chput_ar: string value too long: "',trim(chvalue),'"'
     write(6,*) ' ',trim(znam32),' is CHARACTER*',varlist(jj)%chsize
     ierr=1
  endif

  if(inum.le.0) then
     write(6,*) '?splitn_chput_ar: inum=',inum
     write(6,*) ' number of elements to copy (inum) must be a positive number.'
     ierr=1
  endif

  if(ierr.ne.0) return

  do ii=1,irank
     if((ind0(ii).lt.varlist(jj)%dims(1,ii)).or. &
          (ind0(ii).gt.varlist(jj)%dims(2,ii)) ) then
        ierr=1
        write(6,*) '?splitn_chput_ar: subscript bounds on dimension #',ii
        write(6,*) ' passed value = ',ind0(ii),' not in range ', &
             varlist(jj)%dims(1,ii),':',varlist(jj)%dims(2,ii)
     endif
  enddo
  if(ierr.ne.0) return

  ioff=0
  isize=1
  do ii=irank,1,-1
     ioff = ioff + (ind0(ii)-varlist(jj)%dims(1,ii))*isize
     isize = isize*(varlist(jj)%dims(2,ii)-varlist(jj)%dims(1,ii)+1)
  enddo

  if(ioff+inum.gt.isize) then
     write(6,*) '?splitn_chput_ar: ind0 = (',ind0,'); inum=',inum
     write(6,*) ' write reference beyond end of array.'
     ierr=1
     return
  endif

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock
 
  ilinst = varlist(jj)%nlinadr_st
  iadst  = varlist(jj)%addr_st
  if(iblock.eq.0) then
     ii = varlist(jj)%addr + ioff
     ilinadr = varlist(jj)%nlinadr
  else
     ii = nchv + (iblock-1)*nchv_st + iadst + ioff
     ilinadr = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
  endif

  chbuf(ii:ii+inum-1)=chvalue
  if(varlist(jj)%steerable.eq.2) then
     call chset_ucop(chbuf,ii,inum,iadst+ioff,ilinst+ioff, &
          kupdate,nupdate_max)
  endif

  zdelim = "'"
  if(index(chvalue,zdelim).gt.0) then
     if(index(chvalue,'"').eq.0) zdelim = '"'
  endif

  chval=' '
  chval(1:1)=zdelim
  it=1

  do ic=1,max(1,len_trim(chvalue))
     it=it+1
     chval(it:it)=chvalue(ic:ic)
     if(chvalue(ic:ic).eq.zdelim) then
        it=it+1
        chval(it:it)=chvalue(ic:ic)
     endif
  enddo

  it=it+1
  chval(it:it)=zdelim

  call splitn_put_ar_str(jj,isize,ioff,inum,znam32,chval,ilinadr)

end subroutine splitn_chput_ar

subroutine splitn_put_ar_str(jj,isize,ioff,inum,zname,zvalue,ilinadr)

  ! insert into namelist file text arrays
  !   text to reflect change in namelist data
  !   a minimal edit of existing text is the goal.

  use splitn_module
  implicit NONE

  integer, intent(in) :: jj               ! namelist array variable index
  integer, intent(in) :: isize            ! size of array
  integer, intent(in) :: ioff             ! offset to first assigned element
  integer, intent(in) :: inum             ! number of elements to assign
  character*(*), intent(in) :: zname      ! variable name
  character*(*), intent(in) :: zvalue     ! assigned value
  integer, intent(in) :: ilinadr          ! ptr to text line address

  !  value could have both leading and trailing blanks; these need
  !  to be removed

  integer iv1,iv2,ic,ilen
  integer in1,in2,ir1,ir2,ii1,ii2,i,j,ia
  integer iv1a,iv2a,ir1a,ir2a,ilenc,iv1r,iv2r
  integer knum_new,knum_old
  character*10 zrepeat
  character*30 zindices
  character*140 ztrailer,zcmt
  integer itrailer,itrailv1
  integer ilina,inumpq,imatch,isign
  integer, dimension(:), allocatable :: ilinps,joffs,joffx
  integer idiff,icmt,irank,joff,jnum,koff,knum,koffx,joffz
  integer :: ilinqs

  logical aorder,istraddle

  !-----------------------------------

  allocate(ilinps(2*isize),joffs(2*isize),joffx(2*isize))

  in1 = 1   ! name field assumed to be left justified
  in2 = len_trim(zname)

  iv2 = len_trim(zvalue)

  do ic=1,len(zvalue)
     if(zvalue(ic:ic).ne.' ') then
        iv1=ic
        exit
     endif
  enddo

  ! non-blank part is zvalue(iv1:iv2).

  ilen=iv2-iv1+1

  !-----------
  irank = varlist(jj)%rank  ! dimensionality of array object
  ilina = ilinadr

  inumpq = 0  ! no. of file lines which define elements in this array
  do i=1,isize
     if(ilines(ilina+i-1).ne.0) then
        imatch=0
        do j=1,inumpq
           if(ilines(ilina+i-1).eq.ilinps(j)) then
              imatch=j 
              exit
           endif
        enddo
        if(imatch.eq.0) then
           inumpq=inumpq+1
           joffs(inumpq)=i-1
           ilinps(inumpq)=ilines(ilina+i-1)
           joffx(inumpq)=i-1
        else
           joffx(imatch)=i-1
        endif
     endif
  enddo
  aorder=.TRUE.
  do j=2,inumpq
     if(joffs(j).le.joffx(j-1)) then
        aorder=.FALSE.
        write(6,*) ' %splitn_put_ar_str warning: for namelist array ', &
             trim(zname)
        write(6,*) '  duplicate value assignments cannot be separated;'
        write(6,*) '  new values written to end of file.'
        write(6,*) '  user should edit file to clean up duplications.'
        exit
     endif
  enddo
     
  if((inumpq.eq.0).or.(.not.aorder)) then
     !  new item -- not previously referenced in namelist file text

     ilina = ilinadr + ioff
     call get_newline
     call add_new_line(ioff,inum,newline)  ! put to end of block
     
  else
     !  item previously referenced

     joff=ioff
     koff=ioff+inum  ! beyond end of assignments
     jnum=inum

     if(joff.lt.joffs(1)) then
        !  some elements of index lower than the lowest previously
        !  assigned in the file text ... are assigned here.

        ilina = ilinadr + joff
        
        !  insert in line before the line which previously assigned
        !  the lowest indexed value.
        call add_new_line(joff,min(koff-joff,joffs(1)-joff),ilinps(1))

        joff=joffs(1)   ! start of values still to be assigned
        jnum=koff-joff  ! number of such values (0 if less than 0)

     endif

     !  now loop over previous assignments; if block being reassigned
     !  intersects previous assignment, edit it; if not, insert a new line.

     do j=1,inumpq
        ilinqs = ordl(ilinps(j))
        if(jnum.le.0) exit

        isign=1
        if(lenl(ilinqs).lt.0) isign=-1

        if(j.eq.inumpq) then
           koffx=koff
        else
           koffx=joffs(j+1)
        endif
        if(joff.lt.koffx) then
           ! something needs to be written in vicinity of line (ilinps(j)).

           knum=min(koff-joff,koffx-joff)  ! #of elements to be assigned here
           joffz=joff+knum-1

           if(joff.le.joffx(j)) then
              ! need to modify line (ilinps(j)) stored at (ilinqs)

              lcmt=cmtfld(ilinqs)  ! save off comment (if any)
              if(lcmt.gt.0) then
                 zcmt=textnl(ilinqs)(lcmt:)
                 textnl(ilinqs)(lcmt:)=' '
                 ilenc=len_trim(zcmt)
                 icmt=lcmt
              else
                 zcmt=' '
                 ilenc=0
                 icmt=0
              endif

              ilina = ilinadr

              ! tentative range of value text to be replaced:
              if(irrange(1,ilina+joff).gt.0) then
                 iv1r=irrange(1,ilina+joff)
              else
                 iv1r=ivrange(1,ilina+joff)
              endif
              val2=valfld(2,ilinqs)
              iv2r=val2
              if(joffz.lt.joffx(j)) then
                 if(irrange(1,ilina+joff+1).gt.0) then
                    iv2r=irrange(1,ilina+joff+1)-1
                 else
                    iv2r=ivrange(1,ilina+joff+1)-1
                 endif
              endif

              itrailer=0
              itrailv1=0
              ztrailer=' '
              istraddle=.FALSE.
              if((joffz.lt.joffx(j)).and.(joffs(j).lt.joff)) then
                 if((krepeat(ilina+joff).gt.0).and. &
                      (irrange(1,ilina+joff).eq.irrange(1,ilina+joffz))) then
                    istraddle = .TRUE.
                    !  this is set if the SAME repeat count block covers
                    !  both the first and last elements of the new assignment.
                 endif
              endif

              if(joffz.lt.joffx(j)) then
                 ! at the end of the assignment line are some values that
                 ! are to remain unchanged; find them...
                 knum_new=0
                 ir1a=irrange(1,ilina+joffz)
                 ir2a=irrange(2,ilina+joffz)
                 if((krepeat(ilina+joffz).gt.1).and. &
                      (irrange(1,ilina+joffz).eq.irrange(1,ilina+joffz+1))) &
                      then
                    knum_old=krepeat(ilina+joffz)
                    knum_new=1
                    do
                       if(irrange(1,ilina+joffz).ne. &
                            irrange(1,ilina+joffz+knum_new+1)) exit
                       knum_new=knum_new+1
                    enddo
                 endif
                 if(knum_new.eq.1) then
                    ! lose repeat count
                    itrailv1=2
                    ztrailer= &
                         ','//textnl(ilinqs)(ivrange(1,ilina+joffz):val2)
                    if(.not.istraddle) then
                       textnl(ilinqs)(irrange(1,ilina+joffz):val2)=' '
                    else
                       textnl(ilinqs)(ivrange(2,ilina+joffz)+1:val2)=' '
                    endif
                    krepeat(ilina+joffz+1)=0
                    irrange(1:2,ilina+joffz+1)=0  ! ivrange not known yet
                 else if(knum_new.gt.1) then
                    ! reduce repeat count
                    call gen_zrepeat(knum_new)
                    itrailv1=2+(ir2-ir1)+1
                    ztrailer=','//zrepeat(ir1:ir2)// &
                         textnl(ilinqs)(ivrange(1,ilina+joffz):val2)
                    if(.not.istraddle) then
                       textnl(ilinqs)(irrange(1,ilina+joffz):val2)=' '
                    else
                       textnl(ilinqs)(ivrange(2,ilina+joffz)+1:val2)=' '
                    endif
                    krepeat(ilina+joffz+1:ilina+joffz+knum_new)=knum_new
                    irrange(1,ilina+joffz+1:ilina+joffz+knum_new)= &
                         irrange(2,ilina+joffz+1:ilina+joffz+knum_new) - &
                         (ir2-ir1)   ! final ivrange/irrange not known yet
                 else
                    ! no intersecting repeat count
                    ii1=iv2r+1
                    ztrailer=','//textnl(ilinqs)(ii1:val2)
                    itrailv1=2+ivrange(1,ilina+joffz+1)-ii1
                    textnl(ilinqs)(ii1:val2)=' '
                 endif
                 itrailer=len_trim(ztrailer)
              endif

              if(joff.gt.joffs(j)) then
                 ! at the start of block of values assigned by this
                 ! line are some (at least one) which should not
                 ! be changed.  However, the last sub-block could
                 ! have a repeat count that needs to be modified.
                 knum_new=0
                 if((krepeat(ilina+joff).gt.1).and. &
                      (irrange(1,ilina+joff).eq.irrange(1,ilina+joff-1))) then
                    knum_old=krepeat(ilina+joff)
                    knum_new=1
                    do
                       if(irrange(1,ilina+joff).ne. &
                            irrange(1,ilina+joff-knum_new-1)) exit
                       knum_new=knum_new+1
                    enddo
                 endif
                 ir1a=irrange(1,ilina+joff)
                 ir2a=irrange(2,ilina+joff)
                 iv1a=ivrange(1,ilina+joff)
                 iv2a=ivrange(2,ilina+joff)
                 if(knum_new.eq.1) then
                    ! lose repeat count
                    aline(iv1a:iv2a)=textnl(ilinqs)(iv1a:iv2a)
                    textnl(ilinqs)(ir1a:iv2a)=' '
                    ir2a=ir1a+(iv2a-iv1a)+1
                    textnl(ilinqs)(ir1a:ir2a)=aline(iv1a:iv2a)//','
                    iv1r=ir2a+1  ! replace after...
                    krepeat(ilina+joff-1)=0
                    irrange(1:2,ilina+joff-1)=0
                    ivrange(1,ilina+joff-1)=ir1a
                    ivrange(2,ilina+joff-1)=ir2a-1
                 else if(knum_new.gt.1) then
                    ! reduce repeat count
                    call gen_zrepeat(knum_new)
                    ir1=ir2-(ir2a-ir1a)
                    textnl(ilinqs)(ir1a:ir2a)=zrepeat(ir1:ir2)
                    krepeat(ilina+joff-knum_new:ilina+joff-1)=knum_new
                    textnl(ilinqs)(iv2a+1:iv2a+1)=','
                    iv1r=iv2a+2  ! replace after...
                 endif
              endif

              !  replacement section...

              val1=iv1r
              textnl(ilinqs)(val1:)=' '
              if(knum.gt.1) then
                 call gen_zrepeat(knum)
                 val2=val1+(ir2-ir1)+1+(iv2-iv1)
                 iv1a=val1+(ir2-ir1)+1
                 iv2a=val2
              else
                 val2=val1+iv2-iv1
                 iv1a=val1
                 iv2a=val2
              endif
              call insert_val(ilinps(j),ilinqs,knum,ilina+joff)
              len1=len_trim(textnl(ilinqs))

              ! re-insert trailer if necessary

              if(itrailer.gt.0) then
                 textnl(ilinqs)(len1+1:len1+itrailer)=ztrailer(1:itrailer)
                 idiff=len1+itrailv1-ivrange(1,ilina+joffz+1)
                 do ia=ilina+joffz+1,ilina+joffx(j)
                    if(irrange(1,ia).gt.0) then
                       irrange(1:2,ia)=irrange(1:2,ia)+idiff
                    endif
                    ivrange(1:2,ia)=ivrange(1:2,ia)+idiff
                 enddo
                 len1=len_trim(textnl(ilinqs))
              endif

              valfld(2,ilinqs)=len1

              ! re-insert comment -- decide current or next line.

              cmtfld(ilinqs)=0
              if(ilenc.gt.0) then
                 lcmt = max(lcmt,len1+2)
                 if(lcmt+len_trim(zcmt)-1.lt.80) then
                    textnl(ilinqs)(lcmt:)=zcmt
                    cmtfld(ilinqs)=lcmt
                 else
                    aline=' '
                    aline(icmt:icmt+ilenc-1)=zcmt(1:ilenc)
                    call clear_fields
                    lcmt=icmt
                    newline = ilinps(j)+1
                    call add_nl_line(aline)
                 endif
              endif

              lenl(ilinqs)=isign*max(1,len_trim(textnl(ilinqs)))

           else
              ! new values being assigned are not already assigned in this line
              ! insert new line just after (ilinps(j))
              
              ilina = ilinadr + joff
              call add_new_line(joff,knum,ilinps(j)+1)

           endif

           joff=koffx   ! start of values still to be assigned
           jnum=koff-joff  ! number of such values (0 if less than 0)

        endif
     enddo
              
  endif

  deallocate(ilinps,joffs,joffx)

  contains
     subroutine add_new_line(ioff,inum,ipos)

       integer, intent(in) :: ioff ! offset rel. start address of array
       integer, intent(in) :: inum ! no. of values (repeat count)
       integer, intent(in) :: ipos ! where to put line in file (0 for end)

       integer ip,inump,ishift

       !  write  <name>(<indices>) = <repeat>*<value>
       !  line into file text memory

       newline = ipos

       call clear_fields
       if(inum.eq.1) then
          if(ioff.eq.0) then
             aline = ' '//zname(in1:in2)//' = '//zvalue(iv1:iv2)
             nam1=2
             nam2=2+in2-in1
          else
             call gen_zindices(ioff)
             aline = ' '//zname(in1:in2)//zindices(ii1:ii2)// &
                  ' = '//zvalue(iv1:iv2)
             nam1=2
             nam2=2+in2-in1+1+ii2-ii1
          endif
          val1=nam2+4
          val2=val1+iv2-iv1
          iv1a=val1
          iv2a=val2
       else
          call gen_zrepeat(inum)
          if(ioff.eq.0) then
             aline = ' '//zname(in1:in2)// &
                  ' = '//zrepeat(ir1:ir2)//zvalue(iv1:iv2)
             nam1=2
             nam2=2+in2-in1
          else
             call gen_zindices(ioff)
             aline = ' '//zname(in1:in2)//zindices(ii1:ii2)// &
                  ' = '//zrepeat(ir1:ir2)//zvalue(iv1:iv2)
             nam1=2
             nam2=2+in2-in1+1+ii2-ii1
          endif
          val1=nam2+4
          val2=val1+ir2-ir1+1+iv2-iv1
          ir2a=val1+ir2-ir1
          ir1a=val1
          iv2a=val2
          iv1a=val1+ir2-ir1+1
       endif

       call add_nl_line(aline)

       inump=min(nlines,ipos)

       ! maintain sort order consistency in line pointer arrays

       do ip=1,inumpq
          if(ilinps(ip).ge.inump) ilinps(ip)=ilinps(ip)+1
       enddo

       ilines(ilina:ilina+inum-1)=inump

       ivrange(1,ilina:ilina+inum-1)=iv1a
       ivrange(2,ilina:ilina+inum-1)=iv2a
       if(inum.gt.1) then
          krepeat(ilina:ilina+inum-1)=inum
          irrange(1,ilina:ilina+inum-1)=ir1a
          irrange(2,ilina:ilina+inum-1)=ir2a
       else
          krepeat(ilina)=0
          irrange(1:2,ilina)=0
       endif

     end subroutine add_new_line

     subroutine gen_zrepeat(inum)
       integer, intent(in) :: inum
       integer k
       ! generate repeat count substring & set non-blank limits ir1:ir2

       write(zrepeat,'(I9,"*")') inum
       ir2=10
       do k=1,ir2
          if(zrepeat(k:k).ne.' ') then
             ir1=k
             exit
          endif
       enddo
     end subroutine gen_zrepeat

     subroutine gen_zindices(ioff)
       integer,intent(in) :: ioff   ! offset from array base address

       character*10 zinda(irank)
       integer kblock,ind,iwk,irat,ii,ic,iz

       !  generate array subscript string and set non-blank indices ii1:ii2

       iwk=ioff
       do ii=irank,1,-1
          kblock=(varlist(jj)%dims(2,ii)-varlist(jj)%dims(1,ii))+1
          irat=iwk/kblock
          ind = (iwk - irat*kblock) + varlist(jj)%dims(1,ii)
          write(zinda(ii),'(I10)') ind
          iwk=irat
       enddo

       zindices='('
       iz=1

       do ii=1,irank
          do ic=1,10
             if(zinda(ii)(ic:ic).ne.' ') then
                iz=iz+1
                zindices(iz:iz)=zinda(ii)(ic:ic)
             endif
          enddo
          iz=iz+1
          if(ii.eq.irank) then
             zindices(iz:iz)=')'
          else
             zindices(iz:iz)=','
          endif
       enddo

       ii1=1
       ii2=iz
     end subroutine gen_zindices

     subroutine insert_val(ilinp,ilinq,knum,ilina)
       integer,intent(in) :: ilinp   ! line no. in file
       integer,intent(in) :: ilinq   ! text array address
       integer,intent(in) :: knum    ! no. of values (repeat count)
       integer,intent(in) :: ilina   ! element address in ilines, etc.

       if(knum.gt.1) then
          textnl(ilinq)(val1:val2)= zrepeat(ir1:ir2)//zvalue(iv1:iv2)
          krepeat(ilina:ilina+knum-1)=knum
          irrange(1,ilina:ilina+knum-1)=val1
          irrange(2,ilina:ilina+knum-1)=val1+(ir2-ir1)
       else
          textnl(ilinq)(val1:val2)= zvalue(iv1:iv2)
          krepeat(ilina)=0
          irrange(1,ilina)=0
          irrange(2,ilina)=0
       endif
       ilines(ilina:ilina+knum-1)=ilinp
       ivrange(1,ilina:ilina+knum-1)=iv1a
       ivrange(2,ilina:ilina+knum-1)=iv2a
       
     end subroutine insert_val

end subroutine splitn_put_ar_str
