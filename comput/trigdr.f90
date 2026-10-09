!---------------------------------------
!  these are single precision REAL trig routines
!  sindr,cosdr,tandr -- arguments in degrees
!  asindr,acosdr,atandr -- inverse routines -- result in degrees
!
!  these are intended as (single precision only) substitutes for the
!  non-standard f77 intrinsic functions sind,cosd,tand,asind,acosd,atand
!  the latter not being available on all machines.
!
!  dmc 23 June 1999
!
!---------------------------------------
      real function sindr(zdeg)
      real zdeg
!
!  sine, arg in degrees, single precision
!
      data zpi/3.14159265/

      zrad=zdeg*zpi/180.0
      sindr=sin(zrad)

      return
      end
!---------------------------------------
      real function cosdr(zdeg)
      real zdeg
!
!  cosine, arg in degrees, single precision
!
      data zpi/3.14159265/

      zrad=zdeg*zpi/180.0
      cosdr=cos(zrad)

      return
      end
!---------------------------------------
      real function tandr(zdeg)
      real zdeg
!
!  tangent, arg in degrees, single precision
!
      data zpi/3.14159265/

      zrad=zdeg*zpi/180.0
      tandr=tan(zrad)

      return
      end
!---------------------------------------
      real function asindr(zarg)
!
!  arcsine, returned in degrees, single precision
!
      real zarg

      data zpi/3.14159265/

      zrad = asin(zarg)
      zdeg = zrad*180.0/zpi

      asindr = zdeg
      return
      end
!---------------------------------------
      real function acosdr(zarg)
!
!  arc-cosine, returned in degrees, single precision
!
      real zarg

      data zpi/3.14159265/

      zrad = acos(zarg)
      zdeg = zrad*180.0/zpi

      acosdr = zdeg
      return
      end
!---------------------------------------
      real function atandr(zarg)
!
!  arctangent, returned in degrees, single precision
!
      real zarg

      data zpi/3.14159265/

      zrad = atan(zarg)
      zdeg = zrad*180.0/zpi

      atandr = zdeg
      return
      end
