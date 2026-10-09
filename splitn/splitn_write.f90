subroutine splitn_write(fname,ios)

  !  WRITE a TRANSP namelist to file...

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: fname     ! filename
  integer, intent(out) :: ios            ! completion status code, 0=OK

  !--------------------------------------

  call write_nl(fname,ios)

end subroutine splitn_write

subroutine splitn_write_with_removals(fname,nd,dlist,ios)

  !  WRITE a TRANSP namelist to file, commenting out items dlist(1:nd)
  !  namelist is read back for test purposes and to restore self consistency
  !  of memory representation.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: fname     ! filename
  integer, intent(in) :: nd              ! size of list of deletion names
  character*(*), intent(in) :: dlist(nd) ! names of items to delete
  integer, intent(out) :: ios            ! completion status code, 0=OK

  !--------------------------------------
  character*32 :: dtest1,dtest2
  integer :: i,j,k,imatch,icn1,icn2
  character*150 :: locbuf
  !--------------------------------------

  do i=1,nlines
     j=ordl(i)
     if(eqsfld(j).gt.0) then

        !  name field present in text line; get name
        icn1 = namfld(1,j)
        icn2 = namfld(2,j)
        dtest1 = textnl(j)(icn1:icn2)
        call uupper(dtest1)  ! uppercase

        !  is name in the deletion list?
        imatch = 0
        do k=1,nd
           dtest2 = dlist(k)
           call uupper(dtest2)  ! uppercase
           if(dtest1.eq.dtest2) then
              imatch = k
              exit
           endif
        enddo

        if(imatch.gt.0) then
           !  comment it out
           locbuf = ' !!(ptransp) '//trim(textnl(j))
           textnl(j) = locbuf
           lenl(j) = len(trim(locbuf))
        endif
     endif
  enddo

  ! write...

  call write_nl(fname,ios)

  if(ios.ne.0) then
     write(6,*) ' ?splitn_write_with_removals: error during write... '

  else
     ! read back...

     write(6,*) ' =>splitn_write_with_removals: write OK, readback... '
     call read_nl(fname,ios)
     if(ios.ne.0) then
        write(6,*) ' ?splitn_write_with_removals: error on read back... '
     else
        write(6,*) ' =>splitn_write_with_removals: readback OK... '
     endif
  endif

end subroutine splitn_write_with_removals  

subroutine splitn_write_alt(fname,ios)

  !  WRITE a TRANSP namelist into the splitn_module ...
  !  USING alternate information 

  !  *** this is for debugging only ***

  !  generally the output of this should match the output of a
  !  splitn_write call, but some style changes (handling of comma-
  !  delimiters for example) could occur

  !  specifically, it is known that splitn_write and splitn_write_alt
  !  output won't match if any of the following are true:
  !    a) there are tabs (ctrl^I's) in the input namelist file (tabs
  !       are allowed, but not tracked in the data arrays used here).
  !    b) there are quantities whose values are assigned more than once
  !       in the namelist file (warnings are also generated), even if the
  !       duplicate assignment does not change the value.  A blank field
  !       will appear in the test output.
  !    c) there are instances where the last character of a value
  !       (not the last value on the line) and subsequent delimitting
  !       comma are separated by whitespace (legal but not tracked).
  !    d) thare are instances where the last value on an input line is
  !       followed by a comma (which is legal but not required, and not
  !       saved in the data arrays used by splitn_write_alt).

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: fname     ! filename
  integer, intent(out) :: ios            ! completion status code, 0=OK

  !--------------------------------------
  integer ilf,i,j,k,istop,ilinp,ilinq
  character*150, dimension(:), allocatable :: ztext
  !--------------------------------------

  open(unit=lun,file=fname,status='old',iostat=ios)
  if(ios.eq.0) then
     ! file exists: rename it
     close(unit=lun)
     ilf=len_trim(fname)
     call frename(fname(1:ilf),fname(1:ilf)//'~',ios)
  endif

  open(unit=lun,file=fname,status='new',iostat=ios)
  if(ios.ne.0) return

  !------------> allocate temporary array; set to blank

  allocate(ztext(nlines)); ztext=' '

  !------------> put LHS and equal signs in place

  do i=1,nlines
     j=ordl(i)
     if(eqsfld(j).gt.0) then
        k=eqsfld(j)
        ztext(i)(k:k)='='
        nam1=namfld(1,j)
        nam2=namfld(2,j)
        ztext(i)(nam1:nam2)=textnl(j)(nam1:nam2)
     endif
  enddo

  !------------> insert values with repeat counts & trailing commas

  i=0
  istop=maxlines(0)
  do
     i=i+1
     if(i.gt.istop) exit

     ilinp=ilines(i)
     if(ilinp.eq.0) cycle

     ilinq=ordl(ilinp)

     if(krepeat(i).gt.0) then
        val1=irrange(1,i)
        val2=ivrange(2,i)
        ztext(ilinp)(val1:val2+1)=textnl(ilinq)(val1:irrange(2,i))// &
             textnl(ilinq)(ivrange(1,i):val2)//','
        i=i+krepeat(i)-1
     else
        val1=ivrange(1,i)
        val2=ivrange(2,i)
        ztext(ilinp)(val1:val2+1)=textnl(ilinq)(val1:val2)//','
     endif
  enddo

  !------------> erase end-of-line trailing commas

  do i=1,nlines
     j = max(1,len_trim(ztext(i)))
     if(ztext(i)(j:j).eq.',') ztext(i)(j:j)=' '
  enddo

  !------------> put comments in place

  do i=1,nlines
     j=ordl(i)
     if(cmtfld(j).gt.0) then
        ztext(i)(cmtfld(j):) = textnl(j)(cmtfld(j):)
     endif
  enddo

  !------------> write

  do i=1,nlines
     j = max(1,len_trim(ztext(i)))
     write(lun,'(A)') ztext(i)(1:j)
  enddo

  close(unit=lun)

end subroutine splitn_write_alt
