!  DMC -- get event table data via trdatbuf object

!  See trdatusub/set_sevent.f90 for code where this data is defined
!  See use in datcor/trgdat_exec.for from which these routines are called

subroutine tdb_sevtbl_adr0(d,iadr0)
  use trdatbuf_obj
  implicit NONE

  !  return start location of SEV table
  type (trdatbuf) :: d
  integer, intent(out) :: iadr0

  iadr0 = d%lsevent
end subroutine tdb_sevtbl_adr0

subroutine tdb_sevtbl_ntot(d,ntot)
  use trdatbuf_obj
  implicit NONE

  !  return start location of SEV table
  type (trdatbuf) :: d
  integer, intent(out) :: ntot

  if(d%lsevent.eq.0) then
     ntot=0
  else
     ntot = 1 + d%nsevtbl + d%nsevtbl*d%nsevper*d%nsevwds
  endif

end subroutine tdb_sevtbl_ntot

subroutine tdb_sevtbl_stats(d,nwds,ntabls,ntims)
  use trdatbuf_obj
  implicit NONE

  !  return start location of SEV table
  type (trdatbuf) :: d
  integer, intent(out) :: nwds,ntabls,ntims

  if(d%lsevent.eq.0) then
     nwds=0
     ntabls=0
     ntims=0
  else
     nwds = d%nsevwds
     ntims= d%nsevper
     ntabls = d%nsevtbl
  endif

end subroutine tdb_sevtbl_stats

subroutine tdb_sevtbl_get1(data,ibase, nwds,ntabls,ntims, itabl,itim, dentry)

  ! return a single table entry (nwds words)

  implicit NONE

  real*8, intent(in) :: data(*)
  integer, intent(in) :: ibase   ! location of first word of SEV data
  !  (in d%datbuf this is iadr0 returned by tdb_sevtbl_adr0)

  integer, intent(in) :: nwds,ntabls,ntims  ! SEV data stats
  !  (as returned by tdb_sevtbl_stats)

  integer, intent(in) :: itabl  ! table number (btw 1 and ntabls inclusive)
  integer, intent(in) :: itim   ! time index (btw 1 and ntims inclusive)

  real*8, intent(out) :: dentry(nwds)  ! the table entry 

  ! if any of the stats are zero return all zeroes

  !---------------------------
  integer :: iadr0,iadr1,iadr2
  !---------------------------

  if(nwds.le.0) return
  if(min(ntabls,ntims).le.0) then
     dentry=0.0d0
     return
  endif

  iadr0 = ibase + ntabls + 1  ! address of first table entry

  iadr1 = iadr0 + ((itabl-1)*ntims + (itim-1))*nwds
  iadr2 = iadr1 + nwds - 1

  dentry = data(iadr1:iadr2)

end subroutine tdb_sevtbl_get1
