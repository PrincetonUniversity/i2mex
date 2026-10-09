subroutine i2mex_decimate( &
     & nt1_1, ns_1, t_1, psi_1, p_1, g_1, q_1, x_1, z_1, &
     & nt1_2, ns_2, t_2, psi_2, p_2, g_2, q_2, x_2, z_2, &
     & ier)

  ! Decimate equilibrium data. Too much resolution can cause 
  ! ringing.

  use ezspline_obj
  use ezspline

  implicit none


  integer, parameter :: r8 = selected_real_kind(12,100)

  integer, intent(in) :: nt1_1, ns_1 ! input grid sizes
  real(r8), intent(in) :: t_1(nt1_1), psi_1(ns_1) ! input grid
  real(r8), intent(in), dimension(ns_1) :: p_1, g_1, q_1 !input profiles
  real(r8), intent(in), dimension(nt1_1, ns_1) :: x_1, z_1 ! input mesh

  integer, intent(in) :: nt1_2, ns_2 ! output grid sizes
  real(r8), intent(in) :: t_2(nt1_2), psi_2(ns_2) ! output grid
  real(r8), intent(out), dimension(ns_2) :: p_2, g_2, q_2 ! output profiles
  real(r8), intent(out), dimension(nt1_2, ns_2) :: x_2, z_2 ! output mesh

  integer, intent(out) :: ier

  type(ezspline2) :: spl
  integer iok


  ier = 0

  call i2mex_interp1d(ns_1, psi_1, p_1, ns_2, psi_2, p_2)
  call i2mex_interp1d(ns_1, psi_1, g_1, ns_2, psi_2, g_2)
  call i2mex_interp1d(ns_1, psi_1, q_1, ns_2, psi_2, q_2)

  ! 2-d 

  ! x
  call ezspline_init(spl, nt1_1, ns_1, (/-1,-1/), (/0,0/), iok)
  call ezspline_error(iok)

  spl%x1 = t_1
  spl%x2 = psi_1

  call ezspline_setup(spl, x_1, iok)
  call ezspline_error(iok)

  call ezspline_interp(spl, nt1_2, ns_2, t_2, psi_2, x_2, iok)
  call ezspline_error(iok)

  call ezspline_free(spl, iok)
  call ezspline_error(iok)

  ! z
  call ezspline_init(spl, nt1_1, ns_1, (/-1,-1/), (/0,0/), iok)
  call ezspline_error(iok)

  spl%x1 = t_1
  spl%x2 = psi_1

  call ezspline_setup(spl, z_1, iok)
  call ezspline_error(iok)

  call ezspline_interp(spl, nt1_2, ns_2, t_2, psi_2, z_2, iok)
  call ezspline_error(iok)

  call ezspline_free(spl, iok)
  call ezspline_error(iok)

  if(iok /= 0) ier = 156


end subroutine i2mex_decimate
