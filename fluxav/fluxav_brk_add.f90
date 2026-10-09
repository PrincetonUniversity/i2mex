subroutine fluxav_brk_add(control,brklist,inum,ierr)
 
  use fluxav
  implicit NONE
 
!
!   add internal break point(s) for numeric integration, in the
!   given coordinate direction.
!
!   these would be typically at the spline nodes, but, additional points
!   can be added as needed.
!
!   the idea is to break up integrands to subregions with smoothly varying
!   behaviour, for accuracy of the Gauss-Legendre quadrature method.
!
!   input:
!
  character*(*) control          ! coordinate select ("rho", "theta", "phi")
!
  integer inum                   ! number of break points in list
  real*8 brklist(inum)           ! the break points in *ascending* order
!
!   output:
!
  integer ierr                   ! completion code, 0=OK
!
!-----------------------------------------------------------
!
  integer i,inew
  real*8 brknew(inum+max(n_rhobrk,n_thbrk,n_phibrk))
!
!-----------------------------------------------------------
!
  ierr=0
  do i=1,inum
     if(i.gt.1) then
        if(brklist(i-1).ge.brklist(i)) then
           ierr=1
           write(lunerr,*) ' ?fluxav_brk_add:  ',control,' break list'
           write(lunerr,*) '  out of order:'
           write(lunerr,*) '    brklist(',i-1,')=',brklist(i-1)
           write(lunerr,*) '    brklist(',i,')=',brklist(i)
        endif
     endif
     if(control.eq.'th') then
        if((brklist(i).lt.theta_min).or.(brklist(i).gt.theta_max)) then
           ierr=1
           write(lunerr,*) ' ?fluxav_brk_add:  ',control,' break list item'
           write(lunerr,*) '  out of range:  ',theta_min,' to ',theta_max
           write(lunerr,*) '    brklist(',i,')=',brklist(i)
        endif
     endif
     if(control.eq.'phi') then
        if((brklist(i).lt.phi_min).or.(brklist(i).gt.phi_max)) then
           ierr=1
           write(lunerr,*) ' ?fluxav_brk_add:  ',control,' break list item'
           write(lunerr,*) '  out of range:  ',phi_min,' to ',phi_max
           write(lunerr,*) '    brklist(',i,')=',brklist(i)
        endif
     endif
  enddo
!
  if(ierr.ne.0) return
!
!  OK -- add the items
!
  if(control.eq.'rho') then
     call fluxav_brk_add0(rhobrk,n_rhobrk,brklist,inum,brknew,inew)
     if(inew.gt.n_rhomax) call fluxav_balloc('rho',inew)
     rhobrk(1:inew)=brknew(1:inew)
     n_rhobrk=inew
!
  else if(control.eq.'th') then
     call fluxav_brk_add0(thbrk,n_thbrk,brklist,inum,brknew,inew)
     if(inew.gt.n_thmax) call fluxav_balloc('th',inew)
     thbrk(1:inew)=brknew(1:inew)
     n_thbrk=inew
!
  else if(control.eq.'phi') then
     call fluxav_brk_add0(phibrk,n_phibrk,brklist,inum,brknew,inew)
     if(inew.gt.n_phimax) call fluxav_balloc('phi',inew)
     phibrk(1:inew)=brknew(1:inew)
     n_phibrk=inew
!
  else
     ierr=1
     write(lunerr,*) ' ?fluxav_brk_add:  invalid control string:  ',control
  endif
!
  return
!
end subroutine fluxav_brk_add
