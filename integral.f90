  subroutine i2mex_getPoloidalIntegral(nt1, ns, the, integrand, res, ier)

    ! Compute the poloidal integral of 'integrand' and 
    ! return the results in 'res'. 

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1)
    real(i2mex_r8), intent(in) :: integrand(nt1, ns)
    real(i2mex_r8), intent(out) :: res(ns)
    integer, intent(out) :: ier

    integer i, j
    real(i2mex_r8) :: array(1)

    ier = 0

    do j = 1, ns
       call i2mex_integratePeriodic1d(nt1, the, integrand(1, j), &
            & 1, the(nt1), res(j))
    enddo

!!$    ier = 0
!!$    do j=1, ns
!!$       res(j) = 0.0_i2mex_r8
!!$       do i=1, nt1-1
!!$          res(j) = res(j) + &
!!$               & 0.5_i2mex_r8*( integrand(i+1,j)+integrand(i  ,j) ) * &
!!$               & (the(i+1)-the(i  ))
!!$       enddo
!!$    enddo

  end subroutine i2mex_getPoloidalIntegral
