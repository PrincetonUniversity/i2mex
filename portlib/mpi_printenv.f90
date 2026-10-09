subroutine mpi_printenv0

  use mpi_env_mod
  implicit NONE

  !-----------------------------
  !  find current processor ID; write on unit 0
  !-----------------------------

  call mpi_printenv(-1,0)  ! unknown processor ID; write on unit 0

end subroutine mpi_printenv0

subroutine mpi_printenv(iproc_id,ilun)

  use mpi_env_mod
  implicit NONE

  ! display module contents as formatted ASCII written to unit (ilun)
  ! output is for processor #(iproc_id)

  !-------------------------------------
  ! passed:
  integer, intent(in) :: iproc_id  ! (if unknown give a negative value)
  integer, intent(in) :: ilun      ! Fortran I/O unit number for output

  !-------------------------------------
  ! local:
  integer :: i, ilenn,ilenv

  integer :: myid,mpierr
  !-------------------------------------

  myid = iproc_id
  mpierr=0
  if(myid.lt.0) then
#ifdef __MPI
     call MPI_COMM_RANK(mpi_env_comm,myid,mpierr)
#else
     myid=0
#endif
  endif

  write(ilun,'(A)') ' -------------------------------------------- '
  write(ilun,*) 'process_id: ',myid,' mpi_printenv:'

  do i=1,n_active_names

     ilenn = len(trim(mpi_env_names(i)))
     ilenv = len(trim(mpi_env_vals(i)))

     if(ilenn + ilenv.gt.72) then
        write(ilun,'(2x,a," ="/4x,a)') trim(mpi_env_names(i)), &
             trim(mpi_env_vals(i))
     else
        write(ilun,'(2x,a," = ",a)') trim(mpi_env_names(i)), &
             trim(mpi_env_vals(i))
     endif

  enddo

  write(ilun,'(A)') ' end mpi_printenv'
  write(ilun,'(A)') ' ------------------------ '
  write(ilun,'(A)') ' '

end subroutine mpi_printenv
