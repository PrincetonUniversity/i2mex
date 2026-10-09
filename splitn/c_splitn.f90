
! ----------- I/O ------------
subroutine c_splitn_read(cfilename, ios)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cfilename(*)    ! C name
  integer      :: ios             ! error flag

  character*256 :: filename   ! fortran name

  call cstring(filename,cfilename,'2F')
  call splitn_read(filename,ios)
end subroutine c_splitn_read

subroutine c_splitn_write(cfilename, ios)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cfilename(*)    ! C name
  integer      :: ios             ! error flag

  character*256 :: filename   ! fortran name

  call cstring(filename,cfilename,'2F')
  call splitn_write(filename,ios)
end subroutine c_splitn_write

subroutine c_splitn_merge_write(cfnam_in, cfnam_bck, cfnam_out, ctag, ios)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cfnam_in(*)     ! C name of 2nd input
  character(kind=c_char) :: cfnam_bck(*)    ! C name of backup file
  character(kind=c_char) :: cfnam_out(*)    ! C name of output file
  character(kind=c_char) :: ctag(*)         ! C name of tag
  integer      :: ios             ! error flag

  character*256 :: fnam_in   ! fortran name
  character*256 :: fnam_bck  ! fortran name
  character*256 :: fnam_out  ! fortran name
  character*20  :: tag       ! fortran name

  call cstring(fnam_in,cfnam_in,'2F')
  call cstring(fnam_bck,cfnam_bck,'2F')
  call cstring(fnam_out,cfnam_out,'2F')
  call cstring(tag,ctag,'2F')
  call splitn_merge_write(fnam_in,fnam_bck,fnam_out,tag,ios)
end subroutine c_splitn_merge_write

subroutine c_splitn_edit_enable(cprogram_name)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cprogram_name(*)    ! C name
 
  character*32 :: program_name  ! fortran name
  call cstring(program_name, cprogram_name, '2F')
  call splitn_edit_enable(program_name)
end subroutine c_splitn_edit_enable

! ----------------- info ----------------
subroutine c_splitn_get_info(cname,itype,irank,nd,idims,isize,ichsize,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! input, C name
  integer      :: itype       ! output, 0->int, 1->logical, 2->real, 3->double, 4->character, -1->?
  integer      :: irank       ! output, rank of variable
  integer      :: nd          ! input,  dimension of idims
  integer      :: idims(2,nd) ! output, dimensions (fortran ordering)
  integer      :: isize       ! output, total number of array elements
  integer      :: ichsize     ! output, number of characters for character type
  integer      :: ierr        ! output, nonzero on error

  character*64 :: zname  ! fortran name
  character*8  :: ztype  ! character type
  character*1  :: zt     ! first character of type
  
  call cstring(zname,cname,'2F')
  itype=-1
  irank=0
  ichsize=0
  idims=0
  
  call splitn_get_type(zname,ztype,ichsize,ierr)
  zt = ztype(1:1)
  if (ierr==0) then
    if (zt=='I') then
      itype=0
    else if (zt=='L') then
      itype=1
    else if (zt=='R') then
      itype=2
    else if (zt=='D') then
      itype=3
    else if (zt=='C') then
      itype=4
    else 
      itype=-1
      ierr=1
      return
    end if
    
    call splitn_getdims(zname,irank,nd,idims,isize,ierr)
  end if
end subroutine c_splitn_get_info

! ---------------- get defaults --------------
subroutine c_splitn_igetd(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  integer      :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_igetd(zname,isize,array,ierr)
end subroutine c_splitn_igetd

!
! boolean get which returns integer 0,1 depending on .false.,.true. of data
!
subroutine c_splitn_bgetd(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  integer      :: array(*)    ! data, 0->false, 1->true
  integer      :: ierr        ! nonzero on error

  integer      :: i
  character*64 :: zname          ! fortran name
  logical      :: zarray(isize)  ! fortran logical

  call cstring(zname,cname,'2F')
  call splitn_lgetd(zname,isize,zarray,ierr)
  if (ierr==0) then
     do i=1,isize
        if (zarray(i)) then
           array(i)=1
        else
           array(i)=0
        end if
     end do
  end if
end subroutine c_splitn_bgetd

subroutine c_splitn_lgetd(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  logical      :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_lgetd(zname,isize,array,ierr)
end subroutine c_splitn_lgetd


subroutine c_splitn_rgetd(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  real         :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_rgetd(zname,isize,array,ierr)
end subroutine c_splitn_rgetd


subroutine c_splitn_dgetd(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  real*8       :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_dgetd(zname,isize,array,ierr)
end subroutine c_splitn_dgetd


subroutine c_splitn_cgetd(cname,iwid,isize,carray,ierr)
  use iso_c_binding, only: c_char, c_null_char
  implicit none

  character(kind=c_char) :: cname(*)    ! C name
  integer      :: iwid        ! C width, 1 more then fortran width
  integer      :: isize       ! size of object
  character(kind=c_char) :: carray(*)   ! data
  integer      :: ierr        ! nonzero on error

  integer :: i
  character*64 :: zname               ! fortran name
  character*(iwid-1) :: array(isize)  ! fortran data

  call cstring(zname,cname,'2F')

  call splitn_cgetd(zname,iwid-1,isize,array,ierr)
  if (ierr==0) then
     do i=1,isize
        call cstring(array(i),carray((i-1)*iwid+1),'2C')
        carray(i*iwid)=c_null_char   ! 0 terminated
     end do
  end if
end subroutine c_splitn_cgetd

!
! --------------------- getf -------------------
!
!  The *_*getf(...) subroutines only return value elements explicitly set
!  in the namelist file.  Elements in the array argument which are not 
!  explicitly set are left unmodified.
!
subroutine c_splitn_igetf(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  integer      :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_igetf(zname,isize,array,ierr)
end subroutine c_splitn_igetf


subroutine c_splitn_lgetf(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  logical      :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_lgetf(zname,isize,array,ierr)
end subroutine c_splitn_lgetf


subroutine c_splitn_rgetf(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  real         :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_rgetf(zname,isize,array,ierr)
end subroutine c_splitn_rgetf


subroutine c_splitn_dgetf(cname,isize,array,ierr)
  use iso_c_binding, only: c_char, c_null_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  real*8       :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_dgetf(zname,isize,array,ierr)
end subroutine c_splitn_dgetf


subroutine c_splitn_cgetf(cname,iwid,isize,carray,ierr)
  use iso_c_binding, only: c_char, c_null_char
  implicit none

  character(kind=c_char) :: cname(*)    ! C name
  integer      :: iwid        ! C width, 1 more then fortran width
  integer      :: isize       ! size of object
  character(kind=c_char) :: carray(*)   ! data
  integer      :: ierr        ! nonzero on error

  integer :: i
  character*64 :: zname             ! fortran name
  character*(iwid-1) :: array(isize)  ! fortran data
  character*(iwid-1) :: ctest       ! test string for being set by function call

  do i=1, iwid-1
     ctest(i:i) = "&"
  end do

  call cstring(zname,cname,'2F')
  do i=1,isize
     array(i) = ctest            ! mark string as unset
  end do

  call splitn_cgetf(zname,iwid-1,isize,array,ierr)
  if (ierr==0) then
     do i=1,isize
        if (array(i) /= ctest) then
           call cstring(array(i),carray((i-1)*iwid+1),'2C')  ! only change if different
           carray(i*iwid)=c_null_char   ! 0 terminated
        end if
     end do
  end if
end subroutine c_splitn_cgetf

!
! --------------------- get -------------------
!
!  The *_*get(...) subroutines return all of the values in the namelist buffer.
!
subroutine c_splitn_iget(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  integer      :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_iget(zname,isize,array,ierr)
end subroutine c_splitn_iget

!
! boolean get which returns integer 0,1 depending on .false.,.true. of data
!
subroutine c_splitn_bget(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  integer      :: array(*)    ! data, 0->false, 1->true
  integer      :: ierr        ! nonzero on error

  integer      :: i
  character*64 :: zname          ! fortran name
  logical      :: zarray(isize)  ! fortran logical

  call cstring(zname,cname,'2F')
  call splitn_lget(zname,isize,zarray,ierr)
  if (ierr==0) then
     do i=1,isize
        if (zarray(i)) then
           array(i)=1
        else
           array(i)=0
        end if
     end do
  end if
end subroutine c_splitn_bget


subroutine c_splitn_rget(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  real         :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_rget(zname,isize,array,ierr)
end subroutine c_splitn_rget


subroutine c_splitn_dget(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  real*8       :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_dget(zname,isize,array,ierr)
end subroutine c_splitn_dget

!
!  get an CHARACTER*n namelist object
!  isize *must* match the actual size of the object.
!  iwid must be .ge. the actual width of object elements.
!
subroutine c_splitn_cget(cname,iwid,isize,carray,ierr)
  use iso_c_binding, only: c_char, c_null_char
  implicit none

  character(kind=c_char) :: cname(*)    ! C name
  integer      :: iwid        ! C width, 1 more then fortran width
  integer      :: isize       ! size of object
  character(kind=c_char) :: carray(*)   ! data of size at least iwid*isize
  integer      :: ierr        ! nonzero on error

  integer :: i
  character*64 :: zname             ! fortran name
  character*(iwid-1) :: array(isize)  ! fortran data

  call cstring(zname,cname,'2F')
  call splitn_cget(zname,iwid-1,isize,array,ierr)
  if (ierr==0) then
     do i=1,isize
        call cstring(array(i),carray((i-1)*iwid+1),'2C')  ! only change if different
        carray(i*iwid)=c_null_char   ! 0 terminated
     end do
  end if
end subroutine c_splitn_cget


! --------------------- putw -------------------
subroutine c_splitn_iputw(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  integer      :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_iputw(zname,isize,array,ierr)
end subroutine c_splitn_iputw


subroutine c_splitn_lputw(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  logical      :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_lputw(zname,isize,array,ierr)
end subroutine c_splitn_lputw

!
! boolean put which transfers .false,.true. depending on 0,nonzero of integer argument
!
subroutine c_splitn_bputw(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  integer      :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  integer      :: i
  character*64 :: zname         ! fortran name
  logical      :: zarray(isize) ! fortran array

  call cstring(zname,cname,'2F')
  zarray = (/ (array(i)/=0, i=1,isize) /)
  call splitn_lputw(zname,isize,zarray,ierr)
end subroutine c_splitn_bputw


subroutine c_splitn_rputw(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  real         :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_rputw(zname,isize,array,ierr)
end subroutine c_splitn_rputw


subroutine c_splitn_dputw(cname,isize,array,ierr)
  use iso_c_binding, only: c_char
  implicit none
  character(kind=c_char) :: cname(*)    ! C name
  integer      :: isize       ! size of object
  real*8       :: array(*)    ! data
  integer      :: ierr        ! nonzero on error

  character*64 :: zname  ! fortran name
  call cstring(zname,cname,'2F')
  call splitn_dputw(zname,isize,array,ierr)
end subroutine c_splitn_dputw


subroutine c_splitn_cputw(cname,iwid,isize,carray,ierr)
  use iso_c_binding, only: c_char
  implicit none

  character(kind=c_char) :: cname(*)    ! C name
  integer      :: iwid        ! C width, 1 more then fortran width
  integer      :: isize       ! size of object
  character(kind=c_char) :: carray(*)   ! data
  integer      :: ierr        ! nonzero on error

  integer :: i
  character*64       :: zname         ! fortran name
  character*(iwid-1) :: array(isize)  ! fortran data

  call cstring(zname,cname,'2F')
  do i=1,isize
     call cstring(array(i),carray((i-1)*iwid+1),'2F')
  end do

  call splitn_cputw(zname,iwid-1,isize,array,ierr)
end subroutine c_splitn_cputw
