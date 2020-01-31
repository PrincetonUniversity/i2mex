  subroutine i2mex_scaleG(psii, gi, ier)

    ! Rescale the equilibrium by setting g(psii) = gi. 
    !
    ! *NOTE* This is a rescaling in the magnetic field strength, which 
    ! changes p and psi but leaves q invariant. Since psi is an independent
    ! variable in i2mex, you must make sure that subsequent interpolation
    ! calls use the *new* psi mesh. Therefore, it is in general a good idea
    ! to follow a call to i2mex_scaleG by i2mex_getOriPsi to reset the 
    ! radial grid.
    ! 

    use i2mex_mod
    implicit none
    real(i2mex_r8), intent(in) :: psii, gi
    integer, intent(out) :: ier

    real(i2mex_r8) const, g_old, psi_pos, dum1, dum2
    real(i2mex_r8) psi(i2mex_o%ns), p(i2mex_o%ns), g(i2mex_o%ns)
    integer iok
    
    ier=0
    if(gi <= 0.0_i2mex_r8) return

    psi_pos=psii
    if(psii < i2mex_o%psi_axis) psi_pos=i2mex_o%psi_axis
    if(psii > i2mex_o%psi_edge) psi_pos=i2mex_o%psi_edge

    call i2mex_getG(1, psi_pos, g_old, ier)
    call i2mex_error(iok)

   call i2mex_getOriPsi(i2mex_o%ns, psi, iok)
    call i2mex_error(iok)
    call i2mex_getP(i2mex_o%ns, psi, p, iok)
    call i2mex_error(iok)
    call i2mex_getG(i2mex_o%ns, psi, g, iok)
    call i2mex_error(iok)

    const = gi/g_old

    ! rescale in B

    psi = const    * psi
    i2mex_o%psi_axis = psi(1)
    i2mex_o%psi_edge = psi(i2mex_o%ns)

    p   = const**2 * p
    g   = const    * g

    if(iok/=0) ier = 136
    
    !
    ! recompute spline coefficients for q and g profiles
    ! 
    ! force the spline coefficients to be recomputed (99)

  call eqm_rhofun(99, i2mex_o%id_s, 'PSI', psi, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_psi, iok)
  if(iok /=0) ier = 137

  call eqm_rhofun(99, i2mex_o%id_s, 'G', g, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_g, iok)
  if(iok /=0) ier = 139

  call eqm_rhofun(99, i2mex_o%id_s, 'P', p, i2mex_not_a_knot, &
       & dum1, i2mex_not_a_knot, dum2, &
       & i2mex_o%id_p, iok)
  if(iok /=0) then
     ier = 138
  else
     call eqm_mark_pmhd(i2mex_o%id_p, iok)
     if(iok /=0) ier = 138
  endif

end subroutine i2mex_scaleG
