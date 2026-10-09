subroutine trxplib_clear_mdescr

  ! clear list of machine description files

  use trxplib_ps_options
  implicit NONE

  !---------------------------

  naux_mdescr = 0
  aux_mdescr = ' '

end subroutine trxplib_clear_mdescr

subroutine trxplib_add_mdescr(mdescr_file,ierr)

  ! add to list of machine description files.
  ! re-initialize ps_aux & ps -- if error occurs, restore list to prior state

  !   *** side effects: ps & ps_aux modified in plasma_state_mod ***
  !   For use with trxpl and trxpl2ps only...

  use trxplib_ps_options
  implicit NONE

  character*(*), intent(in) :: mdescr_file
  integer, intent(out) :: ierr

  !-------------------------------

  if(naux_mdescr.eq.max_naux) then
     write(0,*) ' ?? trxplib_ms_mdescr:  too many machine description files. '
     ierr=99
     return
  endif

  naux_mdescr = naux_mdescr + 1
  aux_mdescr(naux_mdescr) = mdescr_file

  call trx_init_state(ierr)

  if(ierr.ne.0) then
     ierr = 1
     aux_mdescr(naux_mdescr) = ' '
     naux_mdescr = naux_mdescr - 1
  endif

end subroutine trxplib_add_mdescr

subroutine trxplib_getnum_mdescr(inum)

  ! return number of recorded machine description files

  use trxplib_ps_options
  implicit NONE

  integer, intent(out) :: inum
  !---------------------------

  inum =  naux_mdescr

end subroutine trxplib_getnum_mdescr

subroutine trxplib_getfull_mdescr(indx,mdescr,ierr)

  ! return the indx'th mdescr file, full file path.
  ! ierr=2 and mdescr is returned blank if indx is out of range
  ! ierr=1 and mdescr is returned blank if no file is found in the
  !   (quasi-hard coded) search path

  use trxplib_ps_options
  implicit NONE

  integer, intent(in) :: indx
  character*(*), intent(out) :: mdescr
  integer,intent(out) :: ierr

  !---------------------------
  !  local...

  integer, parameter :: nsrch = 4
  character*30 :: dpath(nsrch) = (/ "LOCAL/tables/mdescr_nbi       ", &
                                    "CODESYSDIR/tables/mdescr_nbi  ", &
                                    "TRANSP_LOCATION/mdescr        ", &
                                    "                              "/)

  integer :: isrch,ilun,istat,ilen
  character*150 :: locdescr
  !---------------------------

  mdescr = ' '
  if((indx.le.0).OR.(indx.gt.naux_mdescr)) then
     ierr=1
     return
  endif

  locdescr = aux_mdescr(indx)

  call find_io_unit(ilun)

  do isrch=1,nsrch
     ilen=max(1,len(trim(dpath(isrch))))
     call ufilnam(dpath(isrch)(1:ilen),trim(locdescr),mdescr)
     open(unit=ilun,file=mdescr,status='old',action='read',iostat=istat)
     if(istat.eq.0) exit ! success
  enddo

  if(istat.eq.0) then
     close(unit=ilun)
     return  ! found file that exists
  else
     mdescr=' '
     ierr=1
  endif

end subroutine trxplib_getfull_mdescr

subroutine trxplib_getname_mdescr(indx,mdescr)

  ! return the indx'th mdescr file, full file path.
  ! return blank if indx out of range


  use trxplib_ps_options
  implicit NONE

  integer, intent(in) :: indx
  character*(*), intent(out) :: mdescr
  
  !--------------------

  mdescr = ' '
  if(indx.le.0) return
  if(indx.gt.naux_mdescr) return

  mdescr = aux_mdescr(indx)

end subroutine trxplib_getname_mdescr
