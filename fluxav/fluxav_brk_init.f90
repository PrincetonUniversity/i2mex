subroutine fluxav_brk_init(control)
 
  use fluxav
  implicit NONE
 
!
!  initialize break lists to contain 2 elements each:  min & max limiting
!  points of the corresponding coordinate dimension
!
  character*(*) control
!
  integer icheck
!
!---------------------------------------------------
  icheck=0
!
  if((control.eq.'all').or.(control.eq.'rho')) then
     if (n_rhomax.eq.0) call fluxav_balloc('rho',4)
     n_rhobrk=0
     icheck=icheck+1
  endif
!
  if((control.eq.'all').or.(control.eq.'th')) then
     if (n_thmax.eq.0)  call fluxav_balloc('th',4)
     n_thbrk=0
     icheck=icheck+1
  endif
!
  if((control.eq.'all').or.(control.eq.'phi')) then
     if (n_phimax.eq.0) call fluxav_balloc('phi',4)
     n_phibrk=0
     icheck=icheck+1
  endif
!
  if(icheck.eq.0) then
     write(lunerr,*) ' ?fluxav_brk_init:  no action taken:'
     write(lunerr,*) '  invalid control string:  "',control,'"'
  endif
!
  return
!
end subroutine fluxav_brk_init
