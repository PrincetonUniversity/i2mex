!-----------------------------------------------
!  utility subroutine for blims -- find min/max singularity
!  use parabolic fit from 3 data pts
!
subroutine r8_blims_mmx(zth0,zdth,zbm,zb0,zbp,zbsing,zthsing)
  !
  ! input:
  !  zth0 -- theta of ctrmost data pt (zb0)
  !  zdth -- delta(theta) from zth0 to zbm (-) and zbp (+) data points
  !  zbm,zb0,zbp -- sequence of three B field data pts
  !
  !  expected:  either zbm.ge.zb0.and.zbp.ge.zb0
  !               or   zbm.le.zb0.and.zbp.le.zb0
  !
  ! output:
  !  zbsing -- min (or max) B based on parabolic fit
  !  zthsing -- theta location (normalized to range (-pi,pi)) of the
  !            singularity.
  !
  ! fit:
  !
  !  th from zth0-zdth to zth0+zdth
  !
  !  f(th) = zb0 + (zbp-zbm)/(2*zdth) * (th-zth0)
  !              + (zbp+zbm-2*zb0)/(2*zdth**2) * (th-zth0)**2
  !
  !  f'(th) = (zbp-zbm)/(2*zdth) + (zbp+zbm-2*zb0)/(zdth**2) * (th-zth0)
  !
  !  solve f'(th)=0 to find zthsing; evaluate f(zthsing).
  !
  !-------------------
  !  1st check for degenerate case
  !
  !============
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: xpi => pi, twopi
  implicit none
  real(fp) :: zdth,zbm,zb0,zbp,zbsing,zthsing,zth0,zanum
  real(fp) :: za,zb,zdelth
  !
  zanum=(zbp+zbm-2*zb0)
  if(zanum.eq.0.0D0) then
    zbsing=zb0
    zthsing=zth0
    go to 10
  end if
  !
  za=zanum/(2*zdth**2)
  zb=(zbp-zbm)/(2*zdth)
  !
  zdelth = -zdth*zb/(2*za)
  zthsing= zth0 + zdelth
  !
  zbsing = zb0 + zb*zdelth + za*zdelth**2
  !
  !  standardize zthsing range
  !
10 continue
  if(zthsing.lt.-xpi) then
    zthsing=min(xpi,zthsing+twopi)
    go to 10
  else if(zthsing.gt.xpi) then
    zthsing=max(-xpi,zthsing-twopi)
    go to 10
  end if
  !
  return
end subroutine r8_blims_mmx
