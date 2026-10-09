module cont_mod

  use ezspline_obj
  use ezspline

  SAVE
  
  integer, parameter :: cont_r8 = selected_real_kind(12,100)
  real(cont_r8) :: cont_sign
  real(cont_r8) :: cont_q ! estimate of the q (safety factor)
  type(ezspline2) :: fspl

end module cont_mod
