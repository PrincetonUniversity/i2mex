!.........................................................
!            SUBROUTINE TR_PUTNL_FTEXT
!.........................................................
 
subroutine tr_putnl_ftext(zpath,zrunid,iluni,ierr)
!
  use tr_getnl
  implicit NONE
!
  character*(*) zpath         ! path to file (blank if cwd), input
  character*(*) zrunid        ! runid, input
  integer iluni               ! abs(iluni) = lun for file read (input)
  integer ierr                ! completion code:  0=OK
!
! if iluni.lt.0 the final tr_putnl_lines call is skipped
!
!    Write the TRANSP namelist -- filepath and runid provided
!
  integer lt,ilz,ilr,ilzf
  character*196 zfile,zfile2
!
  integer lunzer,ilun
!
!------------------------------------------
!
  ilun=abs(iluni)

  lt=lunzer(0)

  ilr=len_trim(zrunid)
  if(zpath.eq.' ') then
     zfile=zrunid(1:ilr)//'TR.DAT'
     zfile2=zfile
     call ulower(zfile2)
  else
     ilz=len_trim(zpath)
     if((zpath(ilz:ilz).ne.'/').and.(zpath(ilz:ilz).ne.']')) then
        zfile=zpath(1:ilz)//'/'//zrunid(1:ilr)//'TR.DAT'
     else
        zfile=zpath(1:ilz)//zrunid(1:ilr)//'TR.DAT'
     endif
     ilzf=len_trim(zfile)
     zfile2=zfile
     call ulower(zfile2(ilz+1:ilzf))
  endif

  write(*,'(a)')'File opened for writing: '
  write(*,'(a)') zfile

  open(unit=ilun,file=zfile,status='replace',iostat=ierr)
 
  if(iluni.gt.0) call tr_putnl_lines(ilun,ierr)
 
  return
end subroutine tr_putnl_ftext
