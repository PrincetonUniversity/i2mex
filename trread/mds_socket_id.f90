subroutine mds_socket_id(isock_id)

  !  Return the MDS+ socket ID if most recently opened TRANSP run 
  !  was opened using MDS+; if not, return zero.

  use cplotr_mod

  if(nlmds) then

     isock_id = ismds

  else

     isock_id = 0

  endif

end subroutine mds_socket_id
