module trx_bxtr_options

  ! mini-options module for trx_bxtr

  ! extrap_method = 1 for faster, ad hoc method (default: recommended)
  ! extrap_method = 2 for del.B=0, del x B = 0 "while smooth" extrapolator

  integer :: extrap_method = 1  

  ! edge smoothing parameter (meters) -- default: none
  ! note: xplasma will not permit edge_smooth to be more than a small
  !       fraction of the plasma width.

  real*8 :: edge_smooth = 0.0d0

end module trx_bxtr_options
