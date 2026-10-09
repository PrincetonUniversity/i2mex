subroutine errset_mpi(myid,istat)

  ! set an error status: istat>0 expected & enforced
  ! A CALL TO THIS ROUTINE KILLS YOUR RUN... EVEN IF mpi_errstat_enable=.FALSE.

  use mpi_env_mod
  implicit NONE

  integer, intent(in) :: myid   ! process ID: 0,1,2,...<Nprocs>-1
  ! NOTE: pass myid=-1 if the process ID is unknown.

  integer, intent(in) :: istat  ! error code (user defined), .gt.0 expected...

  !-----------------------
  integer :: jstat,id,mpierr
  !-----------------------

  jstat = max(1,istat)

  id=myid
  if(id.lt.0) then
#ifdef __MPI
     call MPI_COMM_RANK(mpi_env_comm,id,mpierr)
#else
     id=0
#endif
  endif

#ifdef __MPI
  write(0,*) ' ERRSET_MPI call by process myid=',id,'; istat = ',istat
  if(istat.ne.jstat) write(0,*) ' ...istat reset to: ',jstat
  write(0,*) ' **** myid=',id,' MPI_ABORT CALL **** '
  call MPI_ABORT(mpi_env_comm,jstat,mpierr)
  call exit(1)
#else
  write(0,*) ' ERRSET called, istat=',istat,' program exiting... '
  if(istat.ne.jstat) write(0,*) ' ...istat reset to: ',jstat
  call exit(1)
#endif

end subroutine errset_mpi
