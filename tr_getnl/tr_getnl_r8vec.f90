subroutine tr_getnl_r8vec(zname,r8vec,maxval,istat)
!
  use tr_getnl
  implicit NONE
!
!  get a vector of REALs corresponding to namelist name "zname"
!
  character*(*), intent(in) :: zname       ! namelist name
  integer, intent(in) :: maxval            ! dimension of results vector
  real*8, dimension(maxval), intent(out) :: r8vec  ! results vector
  integer, intent(out) :: istat
!
!  on output:
!
! normal returns:
!  istat = 0 -- no values found
!  istat = N, 1.le.N.le.maxval -- N values found, stored in r8vec(1:N)
!
! error returns:
!  istat = -1 -- error (e.g. namelist was never read)
!  istat = -N, N.gt.maxval -- N values were found, r8vec(1:maxval) set
!          to the first maxval of them, there are too many values to
!          return them all
!
!-------------------------------------
!
  integer i,ierr,iersum,ils
!
  character(32), dimension(maxval) :: svalues
  integer nvalues
!
  integer lunzer,lt
!-------------------------------------
!
  call tr_getnl_strvals(zname,svalues,maxval,nvalues,ierr)
!
  if(ierr.ne.0) then
     istat=-1
     return
  endif
!
  if(nvalues.eq.0) then
     istat=0
     return
  endif
!
  iersum=0
  do i=1,min(nvalues,maxval)
     read(svalues(i),'(G32.0)',iostat=ierr) r8vec(i)
     if(ierr.ne.0) then
        lt=lunzer(0)
        write(lt,*) '%tr_getnl_r8vec: REAL decode error:'
        ils=len_trim(zname)
        write(lt,*) ' namelist item:  ',zname(1:ils),'(',i,')=',svalues(i)
        iersum=iersum+ierr
     endif
  enddo
!
  if(iersum.gt.0) then
     istat=-1
     return
  endif
!
  if(nvalues.gt.maxval) then
     istat=-nvalues
  else
     istat=nvalues
  endif
!
  return
  end
