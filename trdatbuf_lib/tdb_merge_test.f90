subroutine tdb_merge_test(id,ivala,ivalb,inum,itest)

  ! steering merge test primitive for two trdatbuf objects
  !
  ! return TRUE if item is steerable (value allowed to change)
  !  or if item is not steerable but all values are identical.

  character*(*), intent(in) :: id  ! identifier
  integer, intent(in) :: inum      ! array sizes
  integer, intent(in) :: ivala(inum),ivalb(inum)  ! actual values (dataset a&b)
  logical, intent(out) :: itest    ! result of current test

  !--------------
  logical tdb_steerxcep
  integer :: i,ilev,imin,imax
  !--------------

  itest = tdb_steerxcep(id,ilev)
  if(itest) then
     if(ilev.eq.1) then
        continue  ! change in value is allowed (unconditional)
     else
        ! not allowed to have one element zero and the other non-zero, but,
        ! other changes in value are allowed
        do i=1,inum
           imin=min(abs(ivala(i)),abs(ivalb(i)))
           imax=max(abs(ivala(i)),abs(ivalb(i)))
           if((imin.le.0).and.(imax.ne.0)) then
              itest=.FALSE.
              exit
           endif
        enddo
     endif
  else
     ! change in value not allowed: test
     do i=1,inum
        itest = ivala(i).eq.ivalb(i)
        if(.not.itest) exit
     enddo
  endif

end subroutine tdb_merge_test

subroutine tdb_merge_test1(id,ivala,ivalb,itest)

  ! steering merge test primitive for two trdatbuf objects
  !
  ! return TRUE if item is steerable (value allowed to change)
  !  or if item is not steerable but all values are identical.

  character*(*), intent(in) :: id  ! identifier
  integer, intent(in) :: ivala,ivalb  ! actual values (dataset a&b)
  logical, intent(out) :: itest    ! result of current test

  !--------------
  logical tdb_steerxcep
  integer :: i,ilev,imin,imax
  !--------------

  itest = tdb_steerxcep(id,ilev)
  if(itest) then
     if(ilev.ne.1) then ! change in value is allowed if eq 1 (unconditional)
        ! not allowed to have one element zero and the other non-zero, but,
        ! other changes in value are allowed
       imin=min(abs(ivala),abs(ivalb))
       imax=max(abs(ivala),abs(ivalb))
       if((imin.le.0).and.(imax.ne.0)) then
         itest=.FALSE.
       endif
     endif
  else
     ! change in value not allowed: test
    itest = ivala.eq.ivalb
  endif

end subroutine tdb_merge_test1
