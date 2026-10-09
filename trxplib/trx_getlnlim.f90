subroutine trx_getlnlim(nlinlm,alnlimR,alnlimZ,alnlimt,ierr)
!
!  read TRANSP line limiter information from namelist
!
  implicit NONE
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer, intent(in) :: nlinlm        ! expected number of line limiters
!
! output:
!
  real*8 alnlimR(nlinlm),alnlimZ(nlinlm) ! (R,Z) of point on line limiter
  real*8 alnlimt(nlinlm)               ! angle of orientation of line
!
  integer, intent(out) :: ierr         ! completion code, 0=OK
!
!---------------------------------
!
  integer istat
  integer lt,lunzer
!
!---------------------------------
!
  ierr = 0
!
  lt=lunzer(0)
  call tr_getnl_r8vec('ALNLMR',alnlimR,nlinlm,istat)
  if(istat.ne.nlinlm) then
     write(lt,*) ' %trx_getlnlim:  tr_getnl_r8vec(ALNLMR) returned istat = ', &
          istat
     write(lt,*) '  expected:  istat=nlinlm=',nlinlm
     if(istat.lt.nlinlm) then
        write(lt,*) ' %too few values; 0.0 assumed for remainder.'
        alnlimR(istat+1:nlinlm)=0.0_R8
     else
        write(lt,*) ' ?too many values.'
        ierr=1
        return
     endif
  endif
 
!
  call tr_getnl_r8vec('ALNLMY',alnlimZ,nlinlm,istat)
  if(istat.ne.nlinlm) then
     write(lt,*) ' %trx_getlnlim:  tr_getnl_r8vec(ALNLMZ) returned istat = ', &
          istat
     write(lt,*) '  expected:  istat=nlinlm=',nlinlm
     if(istat.lt.nlinlm) then
        write(lt,*) ' %too few values; 0.0 assumed for remainder.'
        alnlimZ(istat+1:nlinlm)=0.0_R8
     else
        write(lt,*) ' ?too many values.'
        ierr=1
        return
     endif
  endif
 
!
  call tr_getnl_r8vec('ALNLMT',alnlimt,nlinlm,istat)
  if(istat.ne.nlinlm) then
     write(lt,*) ' %trx_getlnlim:  tr_getnl_r8vec(ALNLMT) returned istat = ', &
          istat
     write(lt,*) '  expected:  istat=nlinlm=',nlinlm
     if(istat.lt.nlinlm) then
        write(lt,*) ' %too few values; 0.0 assumed for remainder.'
        alnlimt(istat+1:nlinlm)=0.0_R8
     else
        write(lt,*) ' ?too many values.'
        ierr=1
        return
     endif
  endif
 
!
! convert to meters; standardize orientation angle
!
  alnlimR=0.01_R8*alnlimR
  alnlimZ=0.01_R8*alnlimZ
  alnlimt=180.0_R8-alnlimt
!
  return
  end
