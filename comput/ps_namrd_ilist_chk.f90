subroutine ps_namrd_ilist_chk(listname,nmax,list,nfound,iout,ierr)

  ! determine length of list by number of non-blank elements in the list
  ! once a blank element is detected, all subsequent elements must also
  ! be blank, or, an error flag is set.

  implicit none

  !-----------------
  ! arguments:

  character(len=*), intent(in) :: listname   ! name of list (for error message)
  integer, intent(in) :: nmax                ! size of list array
  character(len=*), intent(in) :: list(nmax) ! the list itself

  integer, intent(out) :: nfound          ! address of last non-blank element

  integer, intent(in) :: iout             ! I/O unit for error messages
  integer, intent(out) :: ierr            ! exit status (0=OK)

  !------------------
  ! local:

  integer :: ii,iblank
  !------------------

  nfound=0
  ierr=0

  iblank=0
  do ii=1,nmax
     if(list(ii).eq.' ') then
        if(iblank.eq.0) iblank=ii
     else
        if(iblank.gt.0) then
           write(iout,*) ' ?ps_namrd_ilist_chk: in list "'//trim(listname)//'":'
           write(iout,*) '  non-blank element #',ii, &
                ' follows blank element #',iblank
           ierr=1
           exit
        endif
     endif
  enddo

  if(ierr.eq.0) then
     nfound = iblank - 1
  endif

end subroutine ps_namrd_ilist_chk
