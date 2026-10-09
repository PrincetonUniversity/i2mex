subroutine cmpi_sset_env(cname_in,cvalue_in,ierr)

  ! C callable mpi_sset_env

  use iso_c_binding, only: c_char
  use mpi_env_mod
  implicit NONE


  character(kind=c_char) cname_in(*), cvalue_in(*)
  integer, intent(out) :: ierr

  !-------------------------------

  character*(maxlen_name) :: name_in
  character*(maxlen_val) :: value

  !-------------------------------

  call cstring(name_in,cname_in,'2F')  ! cname_in -> name_in
  call cstring(value,cvalue_in,'2F')   ! cvalue_in -> value

  call mpi_sset_env(name_in,value,ierr)

end subroutine cmpi_sset_env

subroutine mpi_sset_env(name_in,value,ierr)

  use mpi_env_mod
  implicit NONE

  ! set environment variable module and in OS

  character*(*), intent(in) :: name_in   ! name
  character*(*), intent(in) :: value     ! value
  integer, intent(out) :: ierr

  !--------------------------------------------------------
  integer :: iertmp,ilen,iadr
  logical :: inew
  !--------------------------------------------------------

  call get_loc_myid
  if(mdebug) write(0,*) ' proc #',loc_myid,' entering mpi_sset_env: ', &
       'name_in,value: ',name_in,' ',value

  call legal_name(name_in,ilen,ierr)
  if(ierr.ne.0) return

  call legal_val(value,ilen,ierr)
  if(ierr.ne.0) return

  call lookup_alloc(name_in,iadr,inew)

  if(.not.inew) then
     if(value(1:ilen).ne.mpi_env_vals(iadr)) then
        inew=.TRUE.
     else
        return
     endif
  endif

  mpi_env_vals(iadr) = value(1:ilen)

  call sset_env(name_in,value(1:ilen),iertmp)

  if(mdebug) write(0,*) ' proc #',loc_myid,' exiting mpi_sset_env: ', &
       'ierr: ',ierr

end subroutine mpi_sset_env
