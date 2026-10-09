      subroutine trread_dd(zdisk,zdir,zrunid,ilmds,ierr)

      ! return FDISK, FDIR, and MDS flag setting for most recently opened
      ! TRANSP run

      use cplotr_mod

      implicit NONE

      character*(*), intent(out) :: zdisk  ! run path data "FDISK" in CPLOTR
      character*(*), intent(out) :: zdir   ! run path data "FDIR" in CPLOTR
      character*(*), intent(out) :: zrunid ! run path data "RUNID" in CPLOTR
      logical, intent(out) :: ilmds        ! NLMDS in CPLOTR
      integer, intent(out) :: ierr         ! status code 0=OK

      ierr = 0

      ilmds = nlmds
      zdisk = fdisk
      zdir = fdir
      zrunid = runid

      return
      end

