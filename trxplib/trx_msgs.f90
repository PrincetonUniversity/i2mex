subroutine trx_msgs(lun,zfile)
!
!  redirect warning and error messages to the indicated lun & file.
!  if the file on lun has already been opened, pass zfile = ' '
!
  implicit NONE
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer lun           ! fortran lun for messages (input)
  character*(*) zfile   ! file for messages, or blank (input)
!
!-------------------------------
  integer ierr
!
  call eqm_msgs(lun,zfile,ierr)
  if(ierr.ne.0) then
     write(lun,*) ' ?trx_msgs:  eqm_msgs returned ierr = ',ierr
  endif
!
  call plc_msgs(lun,' ')
!
  return
end subroutine trx_msgs
