subroutine namerefs_upper

  use namerefs
  implicit none

  iupper=1  ! case blind compare of names (set to uppercase during compare)

  return
end subroutine namerefs_upper

subroutine namerefs_lower

  use namerefs
  implicit none

  iupper=0  ! case sensitive compare of names

  return
end subroutine namerefs_lower

subroutine namerefs_free

  use namerefs
  implicit none

  if(nmax.gt.0) then
    nmax=0
    nsize=0
    deallocate(clist,chord,nrefs)
  endif

  return
end subroutine namerefs_free

subroutine namerefs_cln(zstr,jlen)

  !  remove non-printable characters from zstr; left justify zstr

  implicit none
  character(len=*), intent(inout) :: zstr  ! string to be "cleaned"
  integer, intent(out) :: jlen          ! non-blank length

  integer i,ilen,iblank,i0,i1,i2

  ilen=len(zstr)
  iblank=ichar(' ')
  i0=0
  i1=0

  !  convert unprintables to blank; find limits of non-blank characters...
  do i=1,ilen

    if(ichar(zstr(i:i)).lt.iblank) then
      zstr(i:i)=' '
    endif

    if(zstr(i:i).ne.' ') then
      if(i0.eq.0) i0=i
      i1=i
    endif

  enddo

  !  length
  jlen=0
  if(i0.gt.0) then
    jlen=i1-i0+1
  endif

  !  left justification
  if(i0.gt.1) then
    do i=1,jlen
      i2=i0+i-1
      zstr(i:i)=zstr(i2:i2)
    enddo
    zstr(jlen+1:i1)=' '
  endif

  !  within the specified length, replace blanks with underscores
  do i=1,jlen
    if(zstr(i:i).eq.' ') zstr(i:i)='_'
  enddo

  return
end subroutine namerefs_cln
