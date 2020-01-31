  subroutine i2mex_scaleQ(psii, qi, ier)

    ! Apply Bateman's scaling to get q=qi at surface psii

    use i2mex_mod
    implicit none
    real(i2mex_r8), intent(in) :: psii, qi
    integer, intent(out) :: ier

    real(i2mex_r8) const, g_old, psi_pos, q_old, dum1, dum2
    real(i2mex_r8) s(i2mex_o%ns), psi(i2mex_o%ns), g(i2mex_o%ns), q(i2mex_o%ns)
    integer iok, idum
    
    ier=0
    if(qi==0.0_i2mex_r8) return

    psi_pos=psii
    if(psii < i2mex_o%psi_axis) psi_pos=i2mex_o%psi_axis
    if(psii > i2mex_o%psi_edge) psi_pos=i2mex_o%psi_edge

    call i2mex_getG(1, psi_pos, g_old, ier)
    call i2mex_error(ier)
    call i2mex_getQ(1, psi_pos, q_old, ier)
    call i2mex_error(ier)

    const = g_old**2 * (qi**2/q_old**2 - 1.0_i2mex_r8)

!!$    call eq_grid(i2mex_o%id_psi, psi, i2mex_o%ns, idum, iok)
    call i2mex_getOriPsi(i2mex_o%ns, psi, iok)
    if(iok/=0) ier = 49
    call i2mex_getS(i2mex_o%ns, psi, s, ier)
    call i2mex_error(ier)

    call i2mex_getG(i2mex_o%ns, psi, g, ier)
    call i2mex_getQ(i2mex_o%ns, psi, q, ier)
    q = sqrt(g**2 + const) * q / g
    g = sqrt(g**2 + const)
    !
    ! recompute spline coefficients for q and g profiles
    ! 
    ! force the spline coefficients to be recomputed (99)

  call eqm_rhofun(99, i2mex_o%id_s, 'G', g, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_g, iok)
  if(iok /=0) ier = 47

  call eqm_rhofun(99, i2mex_o%id_s, 'Q', q, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_q, iok)
  if(iok /=0) ier = 48

  end subroutine i2mex_scaleQ
