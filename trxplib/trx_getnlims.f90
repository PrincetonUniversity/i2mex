subroutine trx_getnlims(ncirlm,nlinlm,ierr)
!
!  read TRANSP namelist to get NLINLM (# of line limiters)
!                   and to get NCIRLM (# of circle limiters)
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer ncirlm,nlinlm   ! TRANSP limiter counts
  integer ierr            ! completion code:  0 = OK
!
!  *caution* ncirlm=0 and/or nlinlm=0 is possible
!----------------------------------------------------
!
  ierr=0
  call tr_getnl_intvec('NLINLM',nlinlm,1,istat)
  if(istat.eq.0) then
     nlinlm=0
  else if(istat.eq.1) then
     if(nlinlm.lt.0) ierr=ierr+1
  else
     ierr=ierr+1
  endif
!
  call tr_getnl_intvec('NCIRLM',ncirlm,1,istat)
  if(istat.eq.0) then
     ncirlm=0
  else if(istat.eq.1) then
     if(ncirlm.lt.0) ierr=ierr+1
  else
     ierr=ierr+1
  endif
!
  return
  end
