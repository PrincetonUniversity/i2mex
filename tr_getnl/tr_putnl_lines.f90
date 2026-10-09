!.....................................................................
!         SUBROUTINE TR_PUTNL_LINES
!......................................................................
 
subroutine tr_putnl_lines(ilun,ierr)
!
  use tr_getnl
  implicit NONE
!
  integer ilun                ! lun for file write (input argument)
  integer ierr                ! completion code:  0=OK
!
  integer lt,lunzer
!------------------------------
 
  lt = lunzer(0)

  call tr_putnl_write(ilun,ierr)

  if(ierr.ne.0) then
     write(lt,*) ' ?tr_putnl_lines: write error.'
     if(mltext_nlines.eq.0) then
        write(lt,*) '  module does not contain a namelist.'
     endif
  endif
 
  if(allocated(nltext)) deallocate(nltext)
  if(allocated(nltext_lens)) deallocate(nltext_lens)
  if(allocated(mltext)) deallocate(mltext)
  if(allocated(mltext_lens)) deallocate(mltext_lens)
 
  close(unit=ilun)
 
  return
end subroutine tr_putnl_lines

subroutine tr_putnl_write(ilun,ierr)
!
  use tr_getnl
  implicit NONE
!
  integer ilun                ! lun for file write (input argument)
  integer ierr                ! completion code:  0=OK
!
  integer iline,ilen
!------------------------------
  
  if(mltext_nlines.eq.0) then
     ierr = 1

  else
     ierr = 0

     do iline=1,mltext_nlines
        ilen=len_trim(mltext(iline))
        write(ilun,'(A)') mltext(iline)(1:ilen)
        mltext_lens(iline)=len_trim(mltext(iline))
     enddo
  endif

end subroutine tr_putnl_write
