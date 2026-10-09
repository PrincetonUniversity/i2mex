      subroutine mkquin2(sub,x,nx,th,nth,fquin)
use iso_c_binding, only: fp => c_double
!
!  create a data set for interpolation, from evaluation of function and
!  derivatives & 2nd derivatives
!
      external sub                      ! passed subroutine (x,th,output)
      real x(nx)                        ! x coordinate array
      real th(nth)                      ! th coordinate array
!
      real fquin(0:5,nx,nth)            ! function data / spline coeff array
!
!
!  sub's interface:  subroutine sub(xi,th,ra(6))
!    ra(1) -- value of fcn f
!    ra(2) -- df/dx
!    ra(3) -- df/dth
!    ra(4) -- d2f/dx2
!    ra(5) -- d2f/dth2
!    ra(6) -- d2f/dx.dth
!
!----------------------------
!
      do ix=1,nx
         do ith=1,nth
            call sub(x(ix),th(ith),fquin(0,ix,ith))
         end do
      end do
!
      return
      end
