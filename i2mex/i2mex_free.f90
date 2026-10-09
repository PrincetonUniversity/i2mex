subroutine i2mex_free(ier)
  
  ! deallocate buffer and clean up

  use i2mex_mod
  implicit none
  integer, intent(out) :: ier

  ier = 0
  call eqm_free

!!$  i2mex_o%nt1 = 0
!!$  i2mex_o%ns = 0
!!$
  i2mex_o%id_s = 0
  i2mex_o%id_t = 0
  i2mex_o%id_psi = 0
  i2mex_o%id_p = 0
  i2mex_o%id_g = 0
  i2mex_o%id_q = 0
  i2mex_o%id_x = 0
  i2mex_o%id_z = 0
!!$
!!$  i2mex_o%psi_edge = 0.0_i2mex_r8
!!$  i2mex_o%psi_axis = 0.0_i2mex_r8

end subroutine i2mex_free
