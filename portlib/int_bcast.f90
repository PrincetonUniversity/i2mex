subroutine int_bcast(ival,icomm)

  ! broadcast a single integer -- follow with MPI_BARRIER, to be sure
  ! (in a paranoid way) that all processes have received the data.
  ! sender: CPU0; receiver: all others

#ifdef __MPI
  use mpi
#endif
  implicit NONE
  integer :: ival   ! broadcast value (input on cpu0, output on others)
  integer, intent(in) :: icomm   ! MPI communicator

#ifdef __MPI

  integer :: ierr

  call MPI_BCAST(ival,1, MPI_INTEGER, 0, icomm, ierr)
  if(ierr.ne.0) then
     call errmsg_exit(' ?int_bcast: MPI_BCAST error detected!')
  endif

  call MPI_BARRIER(icomm, ierr)
  if(ierr.ne.0) then
     call errmsg_exit(' ?int_bcast: MPI_BARRIER error detected!')
  endif

#endif

end subroutine int_bcast

subroutine intvec_bcast(ivals,numvals,icomm)

  ! broadcast an integer vector -- follow with MPI_BARRIER, to be sure
  ! (in a paranoid way) that all processes have received the data.
  ! sender: CPU0; receiver: all others

#ifdef __MPI
  use mpi
#endif
  implicit NONE
  integer, intent(in) :: numvals ! number of integer values to send/receive
  integer :: ivals(numvals) ! broadcast value (input on cpu0, output on others)
  integer, intent(in) :: icomm   ! MPI communicator

#ifdef __MPI

  integer :: ierr

  call MPI_BCAST(ivals,numvals, MPI_INTEGER, 0, icomm, ierr)
  if(ierr.ne.0) then
     call errmsg_exit(' ?int_bcast: MPI_BCAST error detected!')
  endif

  call MPI_BARRIER(icomm, ierr)
  if(ierr.ne.0) then
     call errmsg_exit(' ?int_bcast: MPI_BARRIER error detected!')
  endif

#endif

end subroutine intvec_bcast
