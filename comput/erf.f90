      function erf(realarg)
!
!  single precision error function
!     simply making use of intrinsic double precision function
!
!      intrinsic derf
      real erf
      real realarg
!
      double precision darg,dans
!
      darg=realarg
      dans=derf(darg)
!
      erf=dans
      return
      end
