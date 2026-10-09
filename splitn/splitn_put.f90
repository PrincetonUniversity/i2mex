! MOD DMC Apr 2009 -- update namelist variables can also be modified.
!  append construct "[<num>]" to the argument (zname) to specify which
!  namelist update block to change.  If the construct is omitted or
!  matches "[0]", then, edit the main namelist (as before).  If "[1]"
!  edit the instance in the 1st update block, "[2]" to edit the instance
!  in the 2nd update block, etc.

!  for non-updatable quantities the appended construct can still be present
!  but it is ignored; the main namelist value (the only one that exists) is
!  modified.


subroutine splitn_iput_sc(zname,ivalue,ierr)

  ! assign a new value to an INTEGER scalar

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of scalar variable
  integer, intent(in) :: ivalue       ! new value to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ilina,iblock,ilinst,iadst
  character*32 znam32,istr
  !--------------------------------------

  call splitn_bdecode1('splitn_iput_sc',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_iput_sc',znam32,jj,ierr)
  if(ierr.ne.0) return

  if(varlist(jj)%type.ne.'I') then
     write(6,*) '?splitn_iput_sc: not INTEGER: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.ne.0) then
     write(6,*) '?splitn_iput_sc: not a SCALAR: ',trim(znam32)
     ierr=1
  endif
  if(ierr.ne.0) return

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock

  ilinst = varlist(jj)%nlinadr_st
  iadst  = varlist(jj)%addr_st
  if(iblock.eq.0) then
     ii = varlist(jj)%addr
     ilina = varlist(jj)%nlinadr
  else
     ii = nint + (iblock-1)*nint_st + iadst
     ilina = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
  endif

  intbuf(ii)=ivalue
  if(varlist(jj)%steerable.eq.2) then
     call iset_ucop(intbuf,ii,1,iadst,ilinst,kupdate,nupdate_max)
  endif

  istr=' '
  write(istr,'(I16)') ivalue
  call splitn_put_sc_str(ilina,znam32,istr)

end subroutine splitn_iput_sc

subroutine splitn_lput_sc(zname,lvalue,ierr)

  ! assign a new value to an LOGICAL scalar

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of scalar variable
  logical, intent(in) :: lvalue       ! new value to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ilina,iblock,ilinst,iadst
  character*32 znam32
  character*2 lstr
  !--------------------------------------

  call splitn_bdecode1('splitn_lput_sc',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_lput_sc',znam32,jj,ierr)
  if(ierr.ne.0) return

  if(varlist(jj)%type.ne.'L') then
     write(6,*) '?splitn_lput_sc: not LOGICAL: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.ne.0) then
     write(6,*) '?splitn_iput_sc: not a SCALAR: ',trim(znam32)
     ierr=1
  endif
  if(ierr.ne.0) return

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock

  ilinst = varlist(jj)%nlinadr_st
  iadst  = varlist(jj)%addr_st
  if(iblock.eq.0) then
     ii = varlist(jj)%addr
     ilina = varlist(jj)%nlinadr
  else
     ii = nlog + (iblock-1)*nlog_st + iadst
     ilina = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
  endif

  logbuf(ii)=lvalue
  if(varlist(jj)%steerable.eq.2) then
     call lset_ucop(logbuf,ii,1,iadst,ilinst,kupdate,nupdate_max)
  endif

  if(lvalue) then
     lstr='.T'
  else
     lstr='.F'
  endif

  call splitn_put_sc_str(ilina,znam32,lstr)

end subroutine splitn_lput_sc

subroutine splitn_rput_sc(zname,rvalue,ierr)

  ! assign a new REAL value to a REAL or REAL*8 scalar

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of scalar variable
  real, intent(in) :: rvalue          ! new value to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ilina,iblock,ilinst,iadst
  character*32 znam32
  character*20 rstr
  character*1 echar
  real*8 :: dtest
  !--------------------------------------

  call splitn_bdecode1('splitn_rput_sc',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_rput_sc',znam32,jj,ierr)
  if(ierr.ne.0) return

  call splitn_put_notilda('splitn_rput_sc',znam32,ierr)
  if(ierr.ne.0) return

  if((varlist(jj)%type.ne.'R').and.(varlist(jj)%type.ne.'D')) then
     write(6,*) '?splitn_rput_sc: not FLOATING POINT: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.ne.0) then
     write(6,*) '?splitn_iput_sc: not a SCALAR: ',trim(znam32)
     ierr=1
  endif
  if(ierr.ne.0) return

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock

  ilinst = varlist(jj)%nlinadr_st
  iadst  = varlist(jj)%addr_st
  if(varlist(jj)%type.eq.'R') then
     if(iblock.eq.0) then
        ii = varlist(jj)%addr
        ilina = varlist(jj)%nlinadr
     else
        ii = nreal + (iblock-1)*nreal_st + iadst
        ilina = nreal+nint+nlog+nr8+nchv + &
             (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
     endif
     rbuf(ii)=rvalue
     echar='e'

     if(varlist(jj)%steerable.eq.2) then
        call rset_ucop(rbuf,ii,1,iadst,ilinst,kupdate,nupdate_max)
     endif
  else
     if(iblock.eq.0) then
        ii = varlist(jj)%addr
        ilina = varlist(jj)%nlinadr
     else
        ii = nr8 + (iblock-1)*nr8_st + iadst
        ilina = nreal+nint+nlog+nr8+nchv + &
             (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
     endif
     dtest=rvalue
     if(varlist(jj)%name.eq.'TINIT') then
        if(dtest.ge.tup(1)) then
           write(6,*) '?splitn_rput_sc: new TINIT value >= update time: ', &
                dtest,tup(1)
           ierr=2
        else
           dbuf(ii)=dtest
           tinit=dtest
        endif
     else
        dbuf(ii)=dtest
     endif
     echar='d'

     if(varlist(jj)%steerable.eq.2) then
        call dset_ucop(dbuf,ii,1,iadst,ilinst,kupdate,nupdate_max)
     endif
  endif
  if(ierr.ne.0) return

  rstr=' '
  write(rstr,'(1pe13.6)') rvalue
  call splitn_fput_clean(rstr,echar)
  call splitn_put_sc_str(ilina,znam32,rstr)

end subroutine splitn_rput_sc

subroutine splitn_dput_sc(zname,dvalue,ierr)

  ! assign a new REAL*8 value to a REAL or REAL*8 scalar

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of scalar variable
  real*8, intent(in) :: dvalue        ! new value to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ilina,iblock,ilinst,iadst
  character*32 znam32
  character*20 dstr
  character*1 echar
  !--------------------------------------

  call splitn_bdecode1('splitn_dput_sc',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_dput_sc',znam32,jj,ierr)
  if(ierr.ne.0) return

  call splitn_put_notilda('splitn_dput_sc',znam32,ierr)
  if(ierr.ne.0) return

  if((varlist(jj)%type.ne.'R').and.(varlist(jj)%type.ne.'D')) then
     write(6,*) '?splitn_dput_sc: not FLOATING POINT: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.ne.0) then
     write(6,*) '?splitn_iput_sc: not a SCALAR: ',trim(znam32)
     ierr=1
  endif
  if(ierr.ne.0) return

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock

  iadst  = varlist(jj)%addr_st
  ilinst = varlist(jj)%nlinadr_st
  if(varlist(jj)%type.eq.'R') then
     if(iblock.eq.0) then
        ii = varlist(jj)%addr
        ilina = varlist(jj)%nlinadr
     else
        ii = nreal + (iblock-1)*nreal_st + iadst
        ilina = nreal+nint+nlog+nr8+nchv + &
             (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
     endif
     rbuf(ii)=dvalue
     echar='e'

     if(varlist(jj)%steerable.eq.2) then
        call rset_ucop(rbuf,ii,1,iadst,ilinst,kupdate,nupdate_max)
     endif
  else
     if(iblock.eq.0) then
        ii = varlist(jj)%addr
        ilina = varlist(jj)%nlinadr
     else
        ii = nr8 + (iblock-1)*nr8_st + iadst
        ilina = nreal+nint+nlog+nr8+nchv + &
             (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
     endif
     if(varlist(jj)%name.eq.'TINIT') then
        if(dvalue.ge.tup(1)) then
           write(6,*) '?splitn_rput_sc: new TINIT value >= update time: ', &
                dvalue,tup(1)
           ierr=2
        else
           dbuf(ii)=dvalue
           tinit=dvalue
        endif
     else
        dbuf(ii)=dvalue
     endif
     echar='d'

     if(varlist(jj)%steerable.eq.2) then
        call dset_ucop(dbuf,ii,1,iadst,ilinst,kupdate,nupdate_max)
     endif
  endif
  if(ierr.ne.0) return

  dstr=' '
  write(dstr,'(1pd19.12)') dvalue
  call splitn_fput_clean(dstr,echar)
  call splitn_put_sc_str(ilina,znam32,dstr)

end subroutine splitn_dput_sc

subroutine splitn_chput_sc(zname,chvalue,ierr)

  ! assign a new value to an CHARACTER*n scalar

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zname  ! name of scalar variable
  character*(*), intent(in) :: chvalue ! new value to assign

  integer, intent(out) :: ierr        ! completion code, 0=OK

  !--------------------------------------
  integer :: jj,ii,ilen,ic,it
  integer :: ilina,iblock,ilinst,iadst
  character*32 znam32
  character*1 zdelim
  !--------------------------------------

  call splitn_bdecode1('splitn_chput_sc',zname,znam32,iblock,ierr)
  if(ierr.ne.0) return

  call splitn_put_chk('splitn_chput_sc',znam32,jj,ierr)
  if(ierr.ne.0) return

  if(varlist(jj)%type(1:1).ne.'C') then
     write(6,*) '?splitn_chput_sc: not CHARACTER: ',trim(znam32)
     ierr=1
  endif

  if(varlist(jj)%rank.ne.0) then
     write(6,*) '?splitn_chput_sc: not a SCALAR: ',trim(znam32)
     ierr=1
  endif

  if(len_trim(chvalue).gt.varlist(jj)%chsize) then
     write(6,*) '?splitn_chput_sc: string value too long: "',trim(chvalue),'"'
     write(6,*) ' ',trim(znam32),' is CHARACTER*',varlist(jj)%chsize
     ierr=1
  endif
  if(ierr.ne.0) return

  if(varlist(jj)%steerable.ne.2) iblock=0
  iblock=max(0,min(nupdate,iblock))
  kupdate = iblock

  iadst  = varlist(jj)%addr_st
  ilinst = varlist(jj)%nlinadr_st
  if(iblock.eq.0) then
     ii = varlist(jj)%addr
     ilina = varlist(jj)%nlinadr
  else
     ii = nchv + (iblock-1)*nchv_st + iadst
     ilina = nreal+nint+nlog+nr8+nchv + &
          (iblock-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + ilinst
  endif

  chbuf(ii)=chvalue
  if(varlist(jj)%steerable.eq.2) then
     call chset_ucop(chbuf,ii,1,iadst,ilinst,kupdate,nupdate_max)
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

  call splitn_put_sc_str(ilina,znam32,chval)

end subroutine splitn_chput_sc

subroutine splitn_put_chk(subname,zname,jj,ierr)

  ! verify existence of a named item; return its index if found.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: subname   ! name of caller (for error msg)
  character*(*), intent(in) :: zname     ! name of namelist variable 
  integer, intent(out) :: jj             ! index to variable descriptor
  integer, intent(out) :: ierr           ! completion code (0=OK)

  !------------------------------
  integer ivar,imatch
  !------------------------------

  ierr = 1

  if(.not.have_database) then
     write(6,*) '?'//trim(subname)//': no namelist has been read.'
     return
  endif

  ierr = 0

  call iorder(zname,ivar,imatch)
  if(imatch.eq.0) then
     write(6,*) '?'//trim(subname)//': unrecognized name: ',trim(zname)
     ierr = 1
  endif

  if(.not.edit_enabled) then
     write(6,*) '?'//trim(subname)//': namelist edit not enabled.'
     ierr = 1
  endif

  if(ierr.ne.0) then
     jj = 0
  else
     jj = var_order(ivar)
     call splitn_put_edit_mark(0)
  endif

end subroutine splitn_put_chk

subroutine splitn_put_edit_mark(iforce)

  use splitn_module
  implicit NONE

  ! if this is the first change:
  ! insert edit marker @ EOF

  integer, intent(in) :: iforce  ! if =1, force a mark; if =0, chk edit_started

  integer :: icur

  if((.not.edit_started).or.(iforce.eq.1)) then
     edit_started=.TRUE.
     !  add a line marking the edit

     call clear_fields
     lcmt=3
     aline = '  ! '//trim(edit_program)//' ('
     icur = len_trim(aline) + 1
     call getlog(aline(icur:len(aline)))
     icur = len_trim(aline) + 1
     aline(icur:icur)=','
     icur = icur+2
     call fdate(aline(icur:len(aline)))
     icur = len_trim(aline) + 1
     aline(icur:icur)=')'

     newline = nlines + 1
     call add_nl_line(aline)  ! add to end of file text array
  endif

end subroutine splitn_put_edit_mark

subroutine splitn_put_notilda(subname,zname,ierr)

  character*(*), intent(in) :: subname  ! name of calling subroutine
  character*(*), intent(in) :: zname    ! variable name
  integer, intent(out) :: ierr          ! status code returned, 0=OK

  ierr=0
  if(zname(1:1).eq.'~') then
     write(6,*) ' ?'//trim(subname)//': update of '//trim(zname)
     write(6,*) '  cannot be done with this call.'
     ierr = 1
  endif

end subroutine splitn_put_notilda

subroutine splitn_fput_clean(dstr,echar)

  ! shorten an encoded floating point string
  ! n.nnnnEmm -> nn.nnn where appropriate; drop trailing 0's, etc.

  !  examples:  1.200000E+01 -> 12.0
  !             7.000000E-01 ->  0.7
  !             6.300000E+23 ->  6.3e23

  implicit NONE
  character*(*), intent(inout) :: dstr  ! floating pt no. in string form
  character*1, intent(in) :: echar      ! character to use for exponent symbol

  !--------------------
  character*4 cexp
  integer iexp,ic,ilen,ie1,ie2,iechar,ine,istop,imove
  !--------------------
  ! first find and decode the exponent; substitute "echar"

  ilen=len(dstr)
  ie1=0
  ie2=0
  do ic=ilen,1,-1
     if(ie2.eq.0) then
        if(dstr(ic:ic).ne.' ') then
           ie2=ic
        endif
     else
        if((dstr(ic:ic).eq.'E').or.(dstr(ic:ic).eq.'e').or. &
             (dstr(ic:ic).eq.'D').or.(dstr(ic:ic).eq.'d')) then
           iechar=ic
           ie1=ic+1
           exit
        endif
     endif
  enddo

  cexp=' '
  ine=ie2-ie1+1
  cexp(4-ine+1:4)=dstr(ie1:ie2)
  read(cexp,'(I4)') iexp  ! encode exponent
  dstr(iechar:)=' '       ! delete exponent from string
  !debug    write(6,*) ' iexp = ',iexp,' dstr = ',dstr

  istop=5
  if(iexp.eq.1) istop=6
  if(iexp.eq.2) istop=7
  if(iexp.eq.-1) istop=4

  !  now strip off any trailing zeroes

  imove=iechar
  do ic=iechar-1,istop,-1
     if(dstr(ic:ic).eq.'0') then
        imove=ic
     else
        exit
     endif
  enddo

  dstr(imove:iechar-1)=' '
  if(iexp.eq.-1) then
     do ic=imove-1,4,-1
        dstr(ic+1:ic+1)=dstr(ic:ic)
     enddo
     dstr(4:4)=dstr(2:2)
     dstr(2:2)='0'
  else if(iexp.eq.0) then
     continue
  else if(iexp.eq.1) then
     dstr(3:3)=dstr(4:4)
     dstr(4:4)='.'
  else if(iexp.eq.2) then
     dstr(3:3)=dstr(4:4)
     dstr(4:4)=dstr(5:5)
     dstr(5:5)='.'
  else
     dstr(imove:imove)=echar
     ic=imove+1
     if(iexp.gt.0) then
        if(iexp.lt.10) then
           write(dstr(ic:ic),'(I1)') iexp
        else if(iexp.lt.100) then
           write(dstr(ic:ic+1),'(I2)') iexp
        else
           write(dstr(ic:ic+2),'(I3)') iexp
        endif
     else
        if(iexp.gt.-10) then
           write(dstr(ic:ic+1),'(I2)') iexp
        else if(iexp.gt.-100) then
           write(dstr(ic:ic+2),'(I3)') iexp
        else
           write(dstr(ic:ic+3),'(I4)') iexp
        endif
     endif
  endif

end subroutine splitn_fput_clean

subroutine splitn_put_sc_str(ilina,zname,zvalue)

  ! insert into namelist file text arrays

  use splitn_module
  implicit NONE

  integer, intent(in) :: ilina            ! ptr to line information
  character*(*), intent(in) :: zname      ! variable name
  character*(*), intent(in) :: zvalue     ! assigned value

  !  value could have both leading and trailing blanks; these need
  !  to be removed

  integer iv1,iv2,ic,ilen
  integer in1,in2
  integer ilinp,ilinq
  integer idiff,icmt,isign,ishift

  !-----------------------------------

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

  ilinp = ilines(ilina)
  if(ilinp.eq.0) then
     !  new item -- not previously referenced in namelist file text
     isign=1
     call clear_fields
     aline = ' '//zname(in1:in2)//' = '//zvalue(iv1:iv2)

     call get_newline

     lcmt = 0

     nam1 = 2
     nam2 = nam1 + (in2-in1)

     val1 = nam2 + 4
     val2 = val1 + (iv2-iv1)
     
     call add_nl_line(aline)

     ilines(ilina) = newline
     irrange(:,ilina) = 0
     krepeat(ilina) = 0
     ivrange(1,ilina) = val1
     ivrange(2,ilina) = val2

  else
     ilinq = ordl(ilinp)
     isign=1
     if(lenl(ilinq).lt.0) isign=-1
     lcmt=cmtfld(ilinq)
     val1=valfld(1,ilinq)
     val2=valfld(2,ilinq)
     textnl(ilinq)(val1:val2)=' '
     val2=val1+iv2-iv1
     valfld(2,ilinq)=val2

     ivrange(1,ilina)=val1
     ivrange(2,ilina)=val2

     len1=len_trim(textnl(ilinq))
     if((lcmt.eq.0).or.(lcmt.gt.val2+1)) then
        textnl(ilinq)(val1:val2)=zvalue(iv1:iv2)
     else
        idiff = val2+2 - lcmt
        if(len1+idiff.lt.80) then
           !  move comment over to the right
           aline(lcmt+idiff:len1+idiff)=textnl(ilinq)(lcmt:len1)
           textnl(ilinq)(lcmt:)=' '
           textnl(ilinq)(val1:val2)=zvalue(iv1:iv2)
           textnl(ilinq)(lcmt+idiff:len1+idiff)=aline(lcmt+idiff:len1+idiff)
           cmtfld(ilinq)=lcmt+idiff
        else
           !  put comment on next line
           newline = ilinp+1
           aline=' '
           aline(lcmt:len1)=textnl(ilinq)(lcmt:len1)
           textnl(ilinq)(lcmt:)=' '
           textnl(ilinq)(val1:val2)=zvalue(iv1:iv2)
           cmtfld(ilinq)=0
           icmt=lcmt
           call clear_fields
           lcmt=icmt
           call add_nl_line(aline)
        endif
     endif
     lenl(ilinq)=isign*max(1,len_trim(textnl(ilinq)))

  endif
end subroutine splitn_put_sc_str
