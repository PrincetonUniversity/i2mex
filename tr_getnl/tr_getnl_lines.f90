!.........................................................................
!         SUBROUTINE TR_GETNL_LINES
!.........................................................................
 
subroutine tr_getnl_lines(ilun,ierr)
!
  use tr_getnl
  implicit NONE
!
  integer ilun                ! lun for file read (input)
  integer ierr                ! completion code:  0=OK
!
  integer iline,lt,lunzer
  character*120 zline
!------------------------------
!
  lt=lunzer(0)
 
10 continue
 
  read(ilun,'(A)',end=99) zline
  nltext_nlines=nltext_nlines+1
  go to 10
 
99 continue
 
  if(nltext_nlines.eq.0) then
     ierr=1
     nltext_nlines=0
     nltext_status=1
     write(lt,*) ' ?? trx_get_nltext:  namelist file has no lines.'
     go to 199
  endif
 
  rewind(ilun)
  if(allocated(nltext)) deallocate(nltext)
  if(allocated(nltext_lens)) deallocate(nltext_lens)
  allocate(nltext(nltext_nlines))
  allocate(nltext_lens(nltext_nlines))
 
  if(allocated(mltext)) deallocate(mltext)
  if(allocated(mltext_lens)) deallocate(mltext_lens)
  allocate(mltext(nltext_nlines+50))
  allocate(mltext_lens(nltext_nlines+50))
  mltext_nlines=nltext_nlines
 
  do iline=1,nltext_nlines
     read(ilun,'(A)') nltext(iline)
     call tr_getnl_cleanup(nltext(iline))
     nltext_lens(iline)=len_trim(nltext(iline))
     mltext(iline)=nltext(iline)
  enddo
 
199 continue
 
  close(unit=ilun)
 
  return
end subroutine tr_getnl_lines
 
subroutine tr_getnl_cleanup(zline)
  ! detab lines
 
  implicit NONE
 
  character*(*) zline
 
  integer i,il
  character*1 ztab,zquot,c
 
  !--------------------------------------
 
  ztab=char(9)
  zquot=' '
 
  il=len_trim(zline)
  do i=1,il
     c=zline(i:i)
     if(zquot.ne.' ') then
        if(c.eq.zquot) zquot=' '  ! closing quote
     else
        if(c.eq.'"') then
           zquot = c              ! opening double quote
        else if(c.eq."'") then
           zquot = c              ! opening single quote
        else if(c.eq.'!') then
           exit                   ! terminating comment
        else if(c.eq.ztab) then
           zline(i:i)=' '         ! detab lines (not cmts or quoted strings)
        endif
     endif
  enddo
 
end subroutine tr_getnl_cleanup
