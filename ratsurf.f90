subroutine i2mex_findRationalSurface(mpols, ntors, psimin, psimax, psis, ier)
  
  ! Search for rational surface psis where q=mpols/ntors in interval 
  ! [psimin, psimax]. 
  ! *NOTE* Only the leftmost rational surface will be found.

  use i2mex_mod
  implicit none

  integer, intent(in) :: mpols ! resonant poloidal mode
  integer, intent(in) :: ntors ! resonant toroidal mode
  real(i2mex_r8), intent(inout) :: psimin ! lower bound
  real(i2mex_r8), intent(inout) :: psimax ! upper bound
  real(i2mex_r8), intent(out) :: psis ! singular surface
  integer, intent(out) :: ier 

  real(i2mex_r8) :: qs

  integer ns, i, iok
  real(i2mex_r8) :: psi2(2), qp2(2), qpp2(2), a, b, c, d, dq, psiqq_a, psiqq_b
  real(i2mex_r8) :: psi3(3), q3(3), qp3(3), tmp
  real(i2mex_r8), dimension(:), allocatable :: psi, q


  ier = 0

  call i2mex_getOriNs(ns, iok)
  call i2mex_error(iok)
  ALLOCATE(psi(ns), q(ns))
  call i2mex_getOriPsi(ns, psi, iok)
  call i2mex_error(iok)
  call i2mex_getQ(ns, psi, q, iok)
  call i2mex_error(iok)

  qs = real(mpols, i2mex_r8)/real(ntors, i2mex_r8)


  psis = psimin
     
  loop: do i = 1, ns-1

     if (psi(i) < psimin) cycle

     if( (qs-q(i))*(q(i+1)-qs) > 0.0_i2mex_r8 ) then
        
        ! rational surface between psi(i) and psi(i+1).

        psimin = psi(i)
        psimax = psi(i+1)

        ! psis is now determined by cubic interpolation (q, psi) 
        ! within this interval using num recipes' formula
        ! on p. 106.

        psi2 = (/ psi(i), psi(i+1) /)
        call i2mex_getQp(2, psi2, qp2, iok)
        call i2mex_getQpp(2, psi2, qpp2, iok)

        psiqq_a = - qpp2(1)/qp2(1)**3
        psiqq_b = - qpp2(2)/qp2(2)**3

        dq = q(i+1) - q(i)
        a = (q(i+1) - qs)/dq
        b = 1.0_i2mex_r8 - a
        c = a*(a**2-1.0_i2mex_r8)*dq**2/6.0_i2mex_r8
        d = b*(b**2-1.0_i2mex_r8)*dq**2/6.0_i2mex_r8
        
        psis = a*psi(i) + b*psi(i+1) + c*psiqq_a + d*psiqq_b

        ! now refine using Newton's scheme

        psi3 = (/psi(i) , psis, psi(i+1)/)
        call i2mex_getQ(3, psi3, q3, iok)
        call i2mex_getQp(3, psi3, qp3, iok)
        if(qp3(2) /= 0.0_i2mex_r8) then
           psis = psis + (qs-q3(2))/qp3(2)
        endif
        
        exit loop

     endif
  enddo loop

     
  DEALLOCATE(psi, q)

  if(psis == psimin) ier = 105

end subroutine i2mex_findRationalSurface

  
  
  
  
  
