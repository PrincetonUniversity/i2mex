      real function qhermitd(xparam,f0,fx0,fxx0,f1,fx1,fxx1,h,hi, &
                             dfdx,d2fdx2)
!
!  evaluate a quintic Hermite interpolation; compute derivatives as well
!  as function value.
!
!  this routine takes a function and its 1st and 2nd derivatives at two
!  successive points on the grid.
!
!  note that the argument x is a **normalized parameter**, x=0 as at the
!  left end of the interval with fcn value f0, x=1 is at the right end of
!  the interval with fcn value f1.  That is, in the caller, there will
!  be some code like
!
!      xparam = (xwant-x(i))/(x(i+1)-x(i))  ! x(i).le.xwant.le.x(i+1)
!
!      h = x(i+1)-x(i)
!      hi = 1/h
!
!      ans=qhermitd(x,f(i),fd(i),fdd(i),f(i+1),fd(i+1),fdd(i+1),h,hi,
!     >       dfdx,d2fdx2)
!
!  the grid spacing h,hi are necessary due to the implicit transformation
!  of variables from [x(i),x(i+1)] to [0,1]; hi=1/h is provided for speed,
!  to avoid extra floating divides on repeated evaluations in the same
!  width intervals.
!
!  dfdx is in "units of" f/h; d2fdx2 in f/(h*h).
!
!  input:
!
      real xparam                       ! value (in [0,1]) to evaluate at
!
      real f0,fx0,fxx0                  ! fcn value, df/dx, d2fdx2 @ xparam=0
      real f1,fx1,fxx1                  ! fcn value, df/dx, d2fdx2 @ xparam=1
!
      real h,hi                         ! actual and inverse grid spacing
!
!  output:
!
!     (function value)                  ! interpolated value of f
      real dfdx                         ! interpolated value of df/dx
      real d2fdx2                       ! interpolated value of d2f/dx2
!
!----------------------------
!
      x=xparam
      x2=x*x
      x3=x2*x
!
      xb=1.0-x
      xb2=xb*xb
      xb3=xb*xb2
!
!  quintic interpolation is done by summing basis polynomials that have
!  nice properties, e.g. alpha(0)=alpha'(0)=alpha''(0)=0,
!                        alpha(1)=1, alpha'(1)=alpha''(1)=0
!
!  there are six such polynomials, each one having unit value or
!  derivative at 0 or 1, and 0 for all other values & derivatives
!
!  alpha(1)=1   alphabar(0)=1
!  beta'(1)=1   betabar'(0)=1
!  gamma''(1)=1 gammabar''(0)=1
!
      alpha=x3*(10.0-x*(15.0-6.0*x))    ! 6*x**5 -15*x**4 +10*x**3
      alphabar=1.0-alpha                ! alphabar = 1 -alpha  is easy to show
!
      dalpha=x2*(30.0-x*(60.0-30.0*x))  ! 30*x**4 -60*x**3 +30*x**2
      dalphabar=-dalpha
!
      ddalpha=x*(60.0-x*(180.0-120.0*x)) ! 120*x**3 -180*x**2 +60*x
      ddalphabar=-ddalpha
!
!
      beta=x3*(-4.0+x*(7.0-3.0*x))      ! -3*x**5 +7*x**4 -4*x**3
      betabar=xb3*(4.0-xb*(7.0-3.0*xb)) ! betabar(x) = -beta(xb)  easy to show
!
      dbeta=x2*(-12.0+x*(28.0-15.0*x))  ! -15*x**4 +28*x**3 -12*x**2
      dbetabar=xb2*(-12.0+xb*(28.0-15.0*xb)) ! betabar'(x) = beta'(xb)
!
      ddbeta=x*(-24.0+x*(84.0-60.0*x))  ! -60*x**3 +84*x**2 -24*x
      ddbetabar=xb*(24.0-xb*(84.0-60.0*xb))
!
!
      gamma=x3*(0.5-x*(1.0-0.5*x))      ! (1/2)*x**5 -x**4 +(1/2)*x**3
      gammabar=xb3*(0.5-xb*(1.0-0.5*xb)) ! gammabar(x) = gamma(xb) easy to show
!
      dgamma=x2*(1.5-x*(4.0-2.5*x))     ! (5/2)*x**4 -4*x**3 +(3/2)*x**2
      dgammabar=xb2*(-1.5+xb*(4.0-2.5*xb)) ! gammabar'(x) = -gamma'(xb)
!
      ddgamma=x*(3.0-x*(12.0-10.0*x))   ! 10*x**3 -12*x**2 + 3*x
      ddgammabar=xb*(3.0-xb*(12.0-10.0*xb))
!
!
!---------------------------
!
!  evaluate the function
!
      sum = alphabar*f0 + alpha*f1 &
         +h*(betabar*fx0 + beta*fx1) &
         +h*h*(gammabar*fxx0 + gamma*fxx1)
!
      qhermitd = sum
!
!
!  evaluate the first derivative
!
      dfdx = hi*(dalphabar*f0 +dalpha*f1) &
         +(dbetabar*fx0 +dbeta*fx1) &
         +h*(dgammabar*fxx0 +dgamma*fxx1)
!
!  evaluate the second derivative
!
      d2fdx2 = hi*hi*(ddalphabar*f0 +ddalpha*f1) &
         +hi*(ddbetabar*fx0 +ddbeta*fx1) &
         +(ddgammabar*fxx0 +ddgamma*fxx1)
!
      return
      end
