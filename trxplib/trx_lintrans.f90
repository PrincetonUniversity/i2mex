subroutine trx_lintrans(id_orig,factor,offset,new_name,id_new,ierr)
!
!  form a new function of the form factor*[old function]+offset
!    i.e. a linear transformation of an old function
!
!  at present only items of form f(rho) are supported, and the
!  interpolation order must be at least 0 (piecewise linear).
!
!  these restrictions can be moved later, if there is a reason to do so...
!  dmc 1 May 2001
!
  implicit NONE
 
  integer, intent(in) :: id_orig        ! xplasma id of original function
  real*8, intent(in) :: factor,offset   ! defining the linear transformation
  character*(*), intent(in) :: new_name ! desired name of new function
 
  integer, intent(out) :: id_new        ! xplasma id of resulting new function
  integer, intent(out) :: ierr          ! completion code, 0=OK
 
!-----------------------------
 
  integer irank,id_rho,id_chi,id_phi,iorder
  integer nrho,idum,ierdum
  integer ibc1,ibc2
 
  integer lunzer
  character*20 zname
!
  real*8 bcs(2),rhobdys(2)
  real*8, dimension(:), allocatable :: rho,f
!-----------------------------
!
!  look up old function...
!
  call eq_xid_mag(id_orig,irank,id_rho,id_chi,id_phi,iorder,ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) '?trx_lintrans:  invalid id_orig:  ',id_orig
     return
  else
     if((irank.ne.1).or.(iorder.lt.0).or.(id_rho.eq.0)) then
        call eq_get_fname(id_orig,zname)
        write(lunzer(0),*) '?trx_lintrans:  cannot transform:  ',zname
        if((irank.ne.1).or.(id_rho.eq.0)) &
             write(lunzer(0),*) ' function is not of the form f(rho).'
        if(iorder.lt.0) &
             write(lunzer(0),*) ' invalid interpolation order: ',iorder
        ierr=1
        return
     endif
  endif
!
!  get old function's rho grid
!
  call eq_ngrid(id_rho,nrho)
!
  allocate(rho(nrho),f(nrho))
  call eq_grid(id_rho,rho,nrho,idum,ierdum)
  rhobdys(1)=rho(1)
  rhobdys(2)=rho(nrho)
!
!  evaluate function on it's old grid
!
  call eq_rgetf(nrho,rho,id_orig,0,f,ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) '?trx_lintrans: eq_rgetf error (1).'
     go to 999
  endif
!
!  get BCs
!
  if(iorder.ge.1) then
     call eq_rgetf(2,rhobdys,id_orig,1,bcs,ierr)
     if(ierr.ne.0) then
        write(lunzer(0),*) '?trx_lintrans: eq_rgetf error (2).'
        go to 999
     endif
     ibc1=1
     ibc2=1
  else
     ibc1=0
     ibc2=0
     bcs=0
  endif
!
!  transform f
!
  f= factor*f + offset
!
!  transform BCs
!
  bcs= factor*bcs
!
!  set up new function...
!
  call eqm_rhofun(iorder,id_rho,new_name,f, &
       ibc1,bcs(1),ibc2,bcs(2), &
       id_new,ierr)
!
!  all done
!
999 continue
  deallocate(f,rho)
  return
end subroutine trx_lintrans
