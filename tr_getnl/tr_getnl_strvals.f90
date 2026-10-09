subroutine tr_getnl_strvals(zname,svalues,maxval,nvalues,ierr)
!
!  find values (as character strings) associated with namelist
!  vector item named by "zname".
!
!  all namelist items are considered to be 1d vectors; scalars are
!  1d vectors of length 1.
!
!  syntax "zname = value"  means assignment to zname(1)
!  syntax "zname = val1,val2" means zname(1)=val1, zname(2)=val2
!  syntax "zname(2) = val1,val2" means zname(2)=val1, zname(3)=val2
!  syntax "zname(3) = 2*val1,val2" means zname(3)=val1, zname(4)=val1, and,
!                                        zname(5)=val2.
!
!  this code is oriented towards parsing of TRANSP namelist (TR.DAT) files.
!  this is *not* a general purpose namelist parser
!
!  the TRANSP namelist TR.DAT files contains comments (starting with "!")
!  the files also contain no more than 1 left-hand-side,
!     and 1 or more assigned values per line
!
  use tr_getnl
  implicit NONE
!
! input:
!
  character*(*) zname       ! name of item wanted
  integer maxval            ! size of value buffer
!
! output:
!
  character*(*) svalues(maxval) ! buffer of values for named item
  integer nvalues           ! largest element of "zname" assigned
!
  integer ierr              ! completion code (0=OK)
!
!   ierr=1 -- namelist not ready (never was read)
!   ierr=2 -- number of values exceeds maxval
!
!   note: nvalues=0 means no values corresponding to "zname" were
!   found.  This is *not* considered an error.
!
!-------------------------------------------------
  integer i,ilen,ic,ier1
!-------------------------------------------------
!
  if((nltext_nlines.eq.0).or.(nltext_status.ne.0)) then
     ierr=1
     return
  else
     ierr=0
  endif
!
! clear values
!
  nvalues=0
  do i=1,maxval
     svalues(i)=' '
  enddo
!
  do i=1,nltext_nlines
     ilen=nltext_lens(i)
     do ic=1,ilen
        if(nltext(i)(ic:ic).eq.'!') exit
        if(nltext(i)(ic:ic).eq.'=') then
           call tr_getnl_1line(i,ic,zname,svalues,maxval,nvalues,ier1)
           ierr=max(ierr,ier1)
        endif
     enddo
  enddo
!
  return
end subroutine tr_getnl_strvals
!--------------------------------------------------------
subroutine tr_getnl_1line(iline,ieqs,zname,svalues,maxval,nvalues,ierr)
!
!  parse one line of the namelist, where ieqs gives location of = sign
!  if LHS matches zname, then determine max index assigned to
!     this could involve parsing the element specification on the LHS
!     and counting the number of values on the RHS
!
  use tr_getnl
  implicit NONE
!
! input:
!
  integer iline             ! line in namelist being scanned
  integer ieqs              ! location of "=" sign.
  character*(*) zname       ! name of item wanted
  integer maxval            ! size of value buffer
!
! output:
!
  character*(*) svalues(maxval) ! buffer of values for named item
  integer nvalues           ! largest element of "zname" assigned
!
  integer ierr              ! completion code (0=OK)
!
!   ierr=1 -- namelist not ready (never was read)
!   ierr=2 -- number of values exceeds maxval
!
!   note: nvalues=0 means no values corresponding to "zname" were
!   found.  This is *not* considered an error.
!
!-------------------------------------------------
  integer ilen,ic,icc,ilparen,irparen
  integer ie1,ie2,ivals,ierck,iblank,icln,irepeat
  integer iv,iv1,iv2
  integer lt,lunzer
!
  character(32) :: ztest1,ztest2
!
!-------------------------------------------------
!
  ierr=0
!
  lt=lunzer(0)
!
  ilen=nltext_lens(iline)
!
  ztest1=zname
  call uupper(ztest1)
!
! get name on LHS (bound by "(" or "=")
!
  ztest2=' '
  icc=0
  ilparen=0
  do ic=1,ieqs-1
     if(nltext(iline)(ic:ic).eq.'(') then
        ilparen=ic
        exit
     endif
     if(nltext(iline)(ic:ic).ne.' ') then
        icc=icc+1
        ztest2(icc:icc)=nltext(iline)(ic:ic)
     endif
  enddo
!
  call uupper(ztest2)
!
! see if name matches...
!
  if(ztest1.ne.ztest2) return        ! not the right item
!
! OK... have a line that assigns to the desired named item
!
  ie1=1
  ie2=0
!
! see if parentheses indicates elements or range of elements in LHS
!
  if(ilparen.gt.0) then
     irparen=0
     do ic=ilparen+1,ieqs-1
        if(nltext(iline)(ic:ic).eq.')') then
           irparen=ic
           exit
        endif
     enddo
     if(irparen.eq.0) then
!
!  oops...
!
        write(lt,*) ' ?tr_getnl_strvalues:  syntax error in line:'
        write(lt,*) ' "',nltext(iline)(1:ilen),'"'
        write(lt,*) ' parentheses left of equal sign.'
        ierr=1
        return
     endif
     if(irparen.lt.(ieqs-1)) then
        if(nltext(iline)(irparen+1:ieqs-1).ne.' ') then
!
!  oops...
!
           write(lt,*) ' ?tr_getnl_strvalues:  syntax error in line:'
           write(lt,*) ' "',nltext(iline)(1:ilen),'"'
           write(lt,*) ' extra chars after closing paren, left of equal sign.'
           ierr=1
           return
        endif
     endif
     iblank=0
     if(irparen.eq.ilparen+1) then
        iblank=1
     else if(nltext(iline)(ilparen+1:irparen-1).eq.' ') then
        iblank=1
     endif
     if(iblank.eq.1) then
!
!  oops...
!
        write(lt,*) ' ?tr_getnl_strvalues:  syntax error in line:'
        write(lt,*) ' "',nltext(iline)(1:ilen),'"'
        write(lt,*) ' nothing between parentheses left of equal sign.'
        ierr=1
        return
     endif
!
!  OK, parse #s btw parens
!
     icln=index(nltext(iline)(ilparen:irparen),':')
     if(icln.eq.0) then
        ztest1=nltext(iline)(ilparen+1:irparen-1)
        read(ztest1,'(I10)',iostat=ierck) ie1
        ierr=max(ierr,ierck)
     else
        ztest1=nltext(iline)(ilparen+1:max((ilparen+1),(icln-1)))
        read(ztest1,'(I10)',iostat=ierck) ie1
        ierr=max(ierr,ierck)
        ztest2=nltext(iline)(icln+1:max((icln+1),(irparen-1)))
        read(ztest2,'(I10)',iostat=ierck) ie2
        ierr=max(ierr,ierck)
        if(ie2.lt.ie1) then
           write(lt,*) ' %2nd index less than first.'
           ierr=1
        endif
        if(ie1.lt.1) then
           write(lt,*) ' %1st index not positive.'
           ierr=1
        endif
        if(ie1.gt.maxval) then
           write(lt,*) ' %1st index exceeds maxval = ',maxval
           ierr=1
        endif
        if(ie2.gt.maxval) then
           write(lt,*) ' %2nd index exceeds maxval = ',maxval
           ierr=1
        endif
     endif
!
     if(ierr.ne.0) then
!
!  oops...
!
        write(lt,*) ' ?tr_getnl_strvalues:  syntax error in line:'
        write(lt,*) ' "',nltext(iline)(1:ilen),'"'
        write(lt,*) ' error reading "',zname,'" element indexing'
        ierr=1
        return
     endif
  endif
!
!  OK...
!    now, ie1=1 & ie2 = 0 if just the name appeared left of = sign
!         ie1=n & ie2 = 0 if "name(n)" appeared left of = sign
!         ie1=n & ie2 = m if "name(n:m)" appeared left of = sign
!
!    now scan values (treated as char strings) right of = sign
!    find the total number of values; if ie2.ne.0 it should match m-n+1.
!    values are delimited by commas
!    values can be repeated by N* syntax
!    blank values not allowed.
!
  ivals=0
  ic=ieqs
  do while(ic.lt.ilen)
     ic=ic+1
     call tr_getnl_nexval(iline,ic,ilen,irepeat,ztest1,ierr)
     if(ierr.ne.0) return       ! blank value
     iv1=ie1+ivals
     iv2=ie1+ivals+irepeat-1
     if(iv2.gt.maxval) then
!
!  oops...
!
        write(lt,*) ' ?tr_getnl_strvalues:  syntax error in line:'
        write(lt,*) ' "',nltext(iline)(1:ilen),'"'
        write(lt,*) ' too many values right of equals sign.'
        write(lt,*) ' maxval = ',maxval
        ierr=1
        return
     endif
     do iv=iv1,iv2
        svalues(iv)=ztest1
     enddo
     nvalues=iv2
     ivals=ivals+irepeat
  enddo
!
  return
end subroutine tr_getnl_1line
!--------------------------------------------------------
subroutine tr_getnl_nexval(iline,ic,ilen,irepeat,zval,ierr)
!
!  find next value for namelist item.  Must be non-blank
!  start scan at char. posn. "ic".
!
  use tr_getnl
  implicit NONE
!
! input:
!
  integer iline            ! index to namelist text line
!
! input/output
!
  integer ic               ! line scanning index
!
! on input:  start of scan
! on output:  loc. of delimitting comma or last non-blank char. +1
!
! input:
!
  integer ilen             ! loc. last non-blank char.
!
! output:
!
  integer irepeat          ! repeat count on value (default: 1)
  character*(*) zval       ! the value (string)
!
  integer ierr             ! completion code, 0=OK
!
! this routine writes an error message and sets ierr if no non-blank
! value is found.
!
!---------------------------------------------
  integer icc,ilval,iq0,imul
  integer lt,lunzer
!
  character*1 zquot
  character*20 zmul
!---------------------------------------------
!
  ierr=0
  irepeat=0
  iq0=0
  zval=' '
  zquot=' '
!
  ilval=len(zval)
!
  icc=0
  ic=ic-1
10 continue
  ic=ic+1
  if(ic.gt.ilen) go to 100
  if(zquot.eq.' ') then
     if(nltext(iline)(ic:ic).eq.',') go to 100
     if(nltext(iline)(ic:ic).eq.'!') then
        ic=ilen+1
        go to 100
     endif
     if(nltext(iline)(ic:ic).eq.'''') zquot=''''
     if(nltext(iline)(ic:ic).eq.'"') zquot='"'
     if(zquot.ne.' ') then
        if(iq0.eq.0) iq0=ic
     endif
  else
     if(nltext(iline)(ic:ic).eq.zquot) zquot=' '
  endif
  if((nltext(iline)(ic:ic).ne.' ').or.(zquot.ne.' ')) then
     icc=icc+1
     if(icc.le.ilval) zval(icc:icc)=nltext(iline)(ic:ic)
  endif
  go to 10
100 continue
  if(icc.eq.0) then
     ierr=1
     lt=lunzer(0)
     write(lt,*) ' ?tr_getnl_strvals:  syntax error in line:'
     write(lt,*) ' "',nltext(iline)(1:ilen),'"'
     write(lt,*) ' blank where value was expected.'
     return
  endif
  if(zquot.ne.' ') then
     ierr=1
     lt=lunzer(0)
     write(lt,*) ' ?tr_getnl_strvals:  syntax error in line:'
     write(lt,*) ' "',nltext(iline)(1:ilen),'"'
     write(lt,*) ' unterminated quote string.'
     return
  endif
  if(icc.gt.ilval) then
     ierr=1
     lt=lunzer(0)
     write(lt,*) ' ?tr_getnl_strvals:  syntax error in line:'
     write(lt,*) ' "',nltext(iline)(1:ilen),'"'
     write(lt,*) ' value field too long.'
     return
  endif
!
! OK -- check for repeat count
!
  if(iq0.eq.0) iq0=icc
  imul=index(zval(1:iq0),'*')
  if(imul.eq.0) then
     irepeat=1
  else if(imul.eq.1) then
     ierr=1
     lt=lunzer(0)
     write(lt,*) ' ?tr_getnl_strvals:  syntax error in line:'
     write(lt,*) ' "',nltext(iline)(1:ilen),'"'
     write(lt,*) ' missing repeat count.'
     return
  else if(imul.eq.icc) then
     ierr=1
     lt=lunzer(0)
     write(lt,*) ' ?tr_getnl_strvals:  syntax error in line:'
     write(lt,*) ' "',nltext(iline)(1:ilen),'"'
     write(lt,*) ' blank value after repeat count.'
     return
  else
     zmul=zval(1:imul-1)
     read(zmul,'(I10)',iostat=ierr) irepeat
     if(ierr.ne.0) then
        ierr=1
        lt=lunzer(0)
        write(lt,*) ' ?tr_getnl_strvals:  syntax error in line:'
        write(lt,*) ' "',nltext(iline)(1:ilen),'"'
        write(lt,*) ' repeat count read error.'
        return
     endif
     zval=zval(imul+1:icc)
  endif
  ierr=0
  return
end subroutine tr_getnl_nexval
