subroutine trx_getcrlim(ncirlm,crlimR,crlimZ,crlimrad,ierr)
!
!  read TRANSP circle limiter information from namelist
!
  implicit NONE
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer, intent(in) :: ncirlm        ! expected number of circular limiters
!
! output:
!
  real*8 crlimR(ncirlm),crlimZ(ncirlm) ! (R,Z) of centers of circles
  real*8 crlimrad(ncirlm)              ! radii of circular limiters
!
  integer, intent(out) :: ierr         ! completion code, 0=OK
!
!---------------------------------
!
  integer istat
  integer lt,lunzer
!
!---------------------------------
  ierr = 0
!
  lt=lunzer(0)
  call tr_getnl_r8vec('CRLMR1',crlimR,ncirlm,istat)
  if(istat.ne.ncirlm) then
     write(lt,*) ' %trx_getcrlim:  tr_getnl_r8vec(CRLMR1) returned istat = ', &
          istat
     write(lt,*) '  expected:  istat=ncirlm=',ncirlm
     if(istat.lt.ncirlm) then
        write(lt,*) ' %too few values; 0.0 assumed for remainder.'
        crlimR(istat+1:ncirlm)=0.0_R8
     else
        write(lt,*) ' ?too many values.'
        ierr=1
        return
     endif
  endif
 
!
  call tr_getnl_r8vec('CRLMY1',crlimZ,ncirlm,istat)
  if(istat.ne.ncirlm) then
     write(lt,*) ' %trx_getcrlim:  tr_getnl_r8vec(CRLMY1) returned istat = ', &
          istat
     write(lt,*) '  expected:  istat=ncirlm=',ncirlm
     if(istat.lt.ncirlm) then
        write(lt,*) ' %too few values; 0.0 assumed for remainder.'
        crlimZ(istat+1:ncirlm)=0.0_R8
     else
        write(lt,*) ' ?too many values.'
        ierr=1
        return
     endif
  endif
 
!
  call tr_getnl_r8vec('CRLMRD',crlimrad,ncirlm,istat)
  if(istat.ne.ncirlm) then
     write(lt,*) ' %trx_getcrlim:  tr_getnl_r8vec(CRLMRD) returned istat = ', &
          istat
     write(lt,*) '  expected:  istat=ncirlm=',ncirlm
     if(istat.lt.ncirlm) then
        write(lt,*) ' %too few values; 0.0 assumed for remainder.'
        crlimrad(istat+1:ncirlm)=0.0_R8
     else
        write(lt,*) ' ?too many values.'
        ierr=1
        return
     endif
  endif
!
! convert to meters
!
  crlimR=0.01_R8*crlimR
  crlimZ=0.01_R8*crlimZ
  crlimrad=0.01_R8*crlimrad
!
  return
  end
