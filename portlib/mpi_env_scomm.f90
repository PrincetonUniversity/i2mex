subroutine mpi_env_comm_get(icomm)

  !  return the mpi_env_mod MPI communicator variable
  !  NOTE: this is MPI_COMM_WORLD, always, with the current version
  !        of the code; there is no provision for changing this.
  !        (See updated source/doc/mpi_transp.doc and the new
  !        mpi_proc_data module for multiple-communicator support).

  use mpi_env_mod
  implicit NONE

  integer, intent(out) :: icomm

  icomm = mpi_env_comm

end subroutine mpi_env_comm_get
