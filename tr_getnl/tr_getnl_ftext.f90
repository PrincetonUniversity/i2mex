!..............................................................
!                  SUBROUTINE TR_GETNL_FTEXT
!..............................................................
 
subroutine tr_getnl_ftext(zpath,zrunid,iluni,ierr)
!
  use tr_getnl
  implicit NONE
!
  character*(*) zpath         ! path to file (blank if cwd), input
  character*(*) zrunid        ! runid, input
  integer iluni               ! abs(iluni) = lun for file read (input)
  integer ierr                ! completion code:  0=OK
!
! if iluni.lt.0 the final tr_getnl_lines call is skipped
!
!    read the TRANSP namelist -- filepath and runid provided
!
  integer lt,ilz,ilr,ilzf
  character*196 zfile,zfile2
!
  integer lunzer,ilun
!
!------------------------------------------
!
  ilun=abs(iluni)
!
  lt=lunzer(0)
!
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
!
  open(unit=ilun,file=zfile,status='old',iostat=ierr)
  if(ierr.ne.0) open(unit=ilun,file=zfile2,status='old',iostat=ierr)
  if(ierr.ne.0) then
     write(lt,*) ' ?? trx_get_nltext:  namelist file open failure.'
     write(lt,*) '    filename was:  ',zfile
     ierr=1
     nltext_nlines=0
     nltext_status=1
     return
  endif
!
!  OK
!
  nltext_status=0
  nltext_nlines=0
!
  if(iluni.gt.0) call tr_getnl_lines(ilun,ierr)
  return
end subroutine tr_getnl_ftext
