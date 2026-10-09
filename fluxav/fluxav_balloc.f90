subroutine fluxav_balloc(control,ineed)
 
  use fluxav
  implicit NONE
 
!   expand the indicated break list
!
  character*(*) control        ! specifies which list
!
!-----------------------------------------------------
  integer, parameter :: isize0 = 32
  integer :: icheck,iprev,isize
  integer, intent(in) :: ineed
!-----------------------------------------------------
!
  icheck=0
  if(control.eq.'rho') then
     if(n_rhomax.eq.0) then
        isize=max(isize0,ineed)
        allocate(rhobrk(isize))
        rhobrk(1:isize)=0.0_R8
        n_rhomax=isize
     else
        iprev = n_rhomax
        n_rhomax=max(ineed,iprev*2)
!
        allocate(rhobrk_tmp(n_rhomax))
        rhobrk_tmp(1:iprev)=rhobrk(1:iprev)
        rhobrk_tmp(iprev+1:n_rhomax)=0.0_R8
        deallocate(rhobrk)
        rhobrk => rhobrk_tmp
     endif
     icheck=icheck+1
!
  else if(control.eq.'th') then
     if(n_thmax.eq.0) then
        isize=max(isize0,ineed)
        allocate(thbrk(isize))
        thbrk(1:isize)=0.0_R8
        n_thmax=isize
     else
        iprev = n_thmax
        n_thmax=max(ineed,iprev*2)
!
        allocate(thbrk_tmp(n_thmax))
        thbrk_tmp(1:iprev)=thbrk(1:iprev)
        thbrk_tmp(iprev+1:n_thmax)=0.0_R8
        deallocate(thbrk)
        thbrk => thbrk_tmp
     endif
     icheck=icheck+1
!
  else if(control.eq.'phi') then
     if(n_phimax.eq.0) then
        isize=max(isize0,ineed)
        allocate(phibrk(isize))
        phibrk(1:isize)=0.0_R8
        n_phimax=isize
     else
        iprev = n_phimax
        n_phimax=max(ineed,iprev*2)
!
        allocate(phibrk_tmp(n_phimax))
        phibrk_tmp(1:iprev)=phibrk(1:iprev)
        phibrk_tmp(iprev+1:n_phimax)=0.0_R8
        deallocate(phibrk)
        phibrk => phibrk_tmp
     endif
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
end subroutine fluxav_balloc
