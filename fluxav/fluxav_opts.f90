subroutine fluxav_lunerr_set(ilun)
 
!  set fortran logical unit no. for error messages
 
  use fluxav
  implicit NONE
  integer ilun
 
  lunerr = ilun
 
  return
 
end subroutine fluxav_lunerr_set
 
 
subroutine fluxav_lunerr_get(ilun)
 
!  get fortran logical unit no. for error messages
 
  use fluxav
  implicit NONE
  integer ilun
 
  ilun = lunerr
 
  return
 
end subroutine fluxav_lunerr_get
 
 
 
subroutine fluxav_thlims_set(thmin,thmax,ierr)
 
!  set upper & lower limits of poloidal angle coordinate
!  these should differ by something very close to twopi.
 
  use fluxav
  implicit NONE
  real*8 thmin,thmax
  integer ierr
 
  ierr = 0
  if(abs((thmax-thmin)-twopi).gt.1.0e-4) then
     write(lunerr,*) ' ?fluxav_thlims_set:  thmax-thmin .ne. twopi:'
     write(lunerr,*) '  thmin = ',thmin,'  thmax = ',thmax
     write(lunerr,*) '   diff = ',thmax-thmin,'  twopi = ',twopi
     write(lunerr,*) '  no action taken.'
     ierr=1
     return
  endif
 
  theta_min = thmin
  theta_max = thmax
 
  return
 
end subroutine fluxav_thlims_set
 
 
subroutine fluxav_phlims_set(phmin,phmax,ierr)
 
!  set upper & lower limits of toroidal angle coordinate
!  these should differ by something very close to twopi.
 
  use fluxav
  implicit NONE
  real*8 phmin,phmax
  integer ierr
 
  ierr = 0
  if(abs((phmax-phmin)-twopi).gt.1.0e-4) then
     write(lunerr,*) ' ?fluxav_phlims_set:  phmax-phmin .ne. twopi:'
     write(lunerr,*) '  phmin = ',phmin,'  phmax = ',phmax
     write(lunerr,*) '   diff = ',phmax-phmin,'  twopi = ',twopi
     write(lunerr,*) '  no action taken.'
     ierr=1
     return
  endif
 
  phi_min = phmin
  phi_max = phmax
 
  return
 
end subroutine fluxav_phlims_set
