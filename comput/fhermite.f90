      real function fhermite(x,zf0,zs0,zf1,zs1)
!
!  C1 Hermite interpolation:  Using cubic basis functions for
!      point and slope matching.
!
!      Interpolation is onto interval [0,1].
!
!      D. McCune 31 Oct 1996
!
!  basis functions:
!     0 at x=0, 1 at x=1, slope 0 at both ends:
!        f1(x)=3*x**2 - 2*x**3
!     1 at x=0, 0 at x=1, slope 0 at both ends:
!        f1(x)=3*xbar**2 - 2*xbar**3
!
      f(x)=x*x*(3.0-2.0*x)
!
!     0 at endpts, slope 1 at x=0, slope 0 at x=1:
!        f1(x)=x - 2*x**2 + x**3
!     0 at endpts, slope 0 at x=0, slope 1 at x=1:
!        f1(x)=xbar - 2*xbar**2 + xbar**3
!
      s(x)=x*(1.0-x*(2.0-x))
!
!-----------------------------------------
!
      if((x.lt.0.0).or.(x.gt.1.0)) then
         write(6,1001) x
 1001    format(' ?fhermite:  x=',1pe14.7,' outside interval [0,1]')
         fhermite=0.0
         return
      end if
!
      xbar=1.0-x
!
      zans=zf0*f(xbar) + zf1*f(x) + zs0*s(x) - zs1*s(xbar)
!
      fhermite=zans
!
      return
      end
