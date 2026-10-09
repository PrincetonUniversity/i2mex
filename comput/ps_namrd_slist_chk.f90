subroutine ps_namrd_slist_chk(listname,nmax,iAMU,iZatom,iZion,nfound,iout,ierr)

  ! determine length of species list by number elements in the list
  ! with positive atomic weight (iAMU).

  ! perform sanity check on Z values.  If Zatom is defined but not Zion,
  ! set Zion=Zatom; if Zion is defined but not Zatom, set Zatom=Zion.

  use iso_c_binding, only: fp => c_double
  implicit none

  !-----------------
  ! arguments:

  character(len=*), intent(in) :: listname   ! name of list (for error message)
  integer, intent(in) :: nmax             ! size of list arrays

  integer, intent(in) :: iAMU(nmax)       ! AMU of species

  ! these can be corrected, hence intent(inout):
  integer, intent(inout) :: iZatom(nmax)     ! Atomic number of species
  integer, intent(inout) :: iZion(nmax)      ! Ionic charge number of species

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
    if(iAMU(ii).le.0) then
      if(iblank.eq.0) iblank=ii
    else
      if(iblank.gt.0) then
        write(iout,*) ' ?ps_namrd_slist_chk: in list "'//trim(listname)//'":'
        write(iout,*) '  non-empty species element #',ii, &
             ' follows empty element #',iblank
        ierr=1
        exit
      endif
    endif
  enddo

  if(ierr.eq.0) then
    nfound = iblank - 1
  else
    return
  endif

  do ii=1,nfound
    if(iZatom(ii).le.0) iZatom(ii)=iZion(ii)
    if(iZion(ii).le.0) iZion(ii)=iZatom(ii)

    if(iZatom(ii).le.0) then
      write(iout,*) ' ?ps_namrd_slist_chk: in list "'//trim(listname)//'":'
      write(iout,*) '  species element #',ii,' has non-zero mass ', &
           ' but atomic number is not positive:'
      ierr=1
      exit
    endif

    if((iZion(ii).lt.1).or.(iZion(ii).gt.iZatom(ii))) then
      write(iout,*) ' ?ps_namrd_slist_chk: in list "'//trim(listname)//'":'
      write(iout,*) '  species element #',ii,' atomic number: ',iZatom(ii)
      write(iout,*) '  species ionic charge: ',iZion(ii), &
           ' should be in range between 1 and the atomic number, but:'
      ierr=1
      exit
    endif

    if((iAMU(ii).lt.iZatom(ii)).or.(iAMU(ii).gt.3*iZatom(ii))) then
      write(iout,*) ' ?ps_namrd_slist_chk: in list "'//trim(listname)//'":'
      write(iout,*) '  species element  #',ii,' fails sanity check:'
      ierr=1
      exit
    endif
  enddo

  if(ierr.ne.0) then
    write(iout,*) '  iAMU=',iAMU(ii),' iZatom=',iZatom(ii),' iZion=',iZion(ii)
  endif

  return
end subroutine ps_namrd_slist_chk
