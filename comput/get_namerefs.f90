subroutine get_namerefs(zname,iref)

  use namerefs
  implicit none

  character(len=*), intent(in) :: zname  ! input name
  integer, intent(out) :: iref           ! output reference count

  !  given name, return reference count.
  !  add to list if necessary...

  !--------------------------------------
  integer :: i,i1,i2,ihalf,indx,jlen
  !--------------------------------------

  wk1=zname
  call namerefs_cln(wk1,jlen)
  wk3=wk1

  if(nsize.eq.0) then

    call expand
    nsize=1
    clist(nsize)=wk1
    chord(nsize)=nsize
    nrefs(nsize)=1
    iref=nrefs(nsize)

  else
    wk2=clist(chord(nsize))
    if(iupper.eq.1) then 
      call uupper(wk1)  ! case blind compare...
      call uupper(wk2)
    endif

    ! now: wk1 = input string converted to UC if needed
    !      wk2 = last string in list converted to UC if needed
    !      wk3 = input string, cleaned but not converted to UC

    if(wk1.gt.wk2) then

      !  beyond last element in list
      if(nsize.eq.nmax) call expand
      nsize=nsize+1
      clist(nsize)=wk3
      chord(nsize)=nsize
      nrefs(nsize)=1
      iref=nrefs(nsize)

    else
      !  binary search to find matching element or element before which
      !  to insert a new name

      i1=1
      i2=nsize
      do
        ihalf = (i1+i2)/2
        indx=chord(ihalf)
        wk2=clist(indx)
        if(iupper.eq.1) call uupper(wk2)
        if(i1.eq.i2) exit

        if(wk1.le.wk2) then
          i2=ihalf
        else
          i1=ihalf+1
        endif
      enddo

      if(wk1.eq.wk2) then
        nrefs(indx)=nrefs(indx)+1  ! match: increment reference count
        iref=nrefs(indx)
      else
        if(nsize.eq.nmax) call expand
        nsize=nsize+1              ! add new entry to list
        clist(nsize)=wk3
        nrefs(nsize)=1
        iref=nrefs(nsize)
        do i=nsize,i1+1,-1         ! keep ordering
          chord(i)=chord(i-1)
        enddo
        chord(i1)=nsize
      endif

    endif
  endif

  return
end subroutine get_namerefs

subroutine get_namerefs_str(instr,iref,outstr)

  !  given input string <name> 
  !  return             <name>_#
  !    where # = reference count returned by get_namerefs

  use namerefs
  implicit none

  character(len=*), intent(in) :: instr   ! input string
  integer, intent(out) :: iref            ! reference count
  character(len=*), intent(out) :: outstr ! output string

  !---------------------------------------------
  integer :: ilen,idigs,jlen
  character(len=10) :: zdigs
  !---------------------------------------------

  call get_namerefs(instr,iref)

  outstr = ' '
  ilen=len(outstr)

  wk3=instr
  call namerefs_cln(wk3,jlen)
  !  namerefs'wk3 contains the string, "cleaned up..."; jlen=non-blank
  !  length after left justification...

  if(iref.le.9) then
     idigs=1; write(zdigs(1:idigs),'(I1)') iref
  else if(iref.le.99) then
     idigs=2; write(zdigs(1:idigs),'(I2)') iref
  else if(iref.le.999) then
     idigs=3; write(zdigs(1:idigs),'(I3)') iref
  else if(iref.le.9999) then
     idigs=4; write(zdigs(1:idigs),'(I4)') iref
  else if(iref.le.99999) then
     idigs=5; write(zdigs(1:idigs),'(I5)') iref
  else if(iref.le.999999) then
     idigs=6; write(zdigs(1:idigs),'(I6)') iref
  else if(iref.le.9999999) then
     idigs=7; write(zdigs(1:idigs),'(I7)') iref
  else if(iref.le.99999999) then
     idigs=8; write(zdigs(1:idigs),'(I8)') iref
  else if(iref.le.999999999) then
     idigs=9; write(zdigs(1:idigs),'(I9)') iref
  else
     idigs=10; write(zdigs(1:idigs),'(I10)') iref
  endif

  jlen=min(jlen,ilen-idigs-1)

  if(jlen.lt.0) then
     outstr='?'
  else if(jlen.eq.0) then
     outstr='_'//zdigs(1:idigs)
  else
     outstr=wk3(1:jlen)//'_'//zdigs(1:idigs)
  endif

  return
end subroutine get_namerefs_str
