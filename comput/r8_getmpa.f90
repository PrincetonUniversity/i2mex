subroutine r8_getmpa(zrmc,zymc,mimom,imom,zrin,zrout,zymp)
  !
  !  DMC 18 March 1994:
  !  given the updown asymmetric moments of a surface, define an
  !  approximate midplane to the surface.  Return the inner and
  !  outer R intercepts of this midplane with the surface, and
  !  return the height of the surface above y=0.
  !
  !  rev dmc 25 Mar 1994:  use
  !
  !============
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, twopi
  implicit none
  integer :: imom,mimom,imm,ith
  real(fp) :: zth,zr,zy
  real(fp), dimension(0:mimom,2) :: zrmc  ! R moments
  real(fp), dimension(0:mimom,2) :: zymc  ! Y moments
  !
  !  mimom gives the moments 1st array dimension; imom gives the number
  !   of *active* moments; imom.le.mimom ...
  !
  real(fp) :: zrin            ! inner intercept R (output)
  real(fp) :: zrout           ! outer intercept R (output)
  real(fp) :: zymp            ! midplane vertical displacement (output)
  !
  integer, parameter :: inth = 400
  real(fp), dimension(inth) :: zrcon, zycon
  real(fp), dimension(0:imom) :: zcos,zsin
  !---------------------------------------------------------------------
  !
  zcos(0)=1
  zsin(0)=0
  do ith=1,inth-1
    zth=(twopi*(ith-1))/(inth-1)
    zr=zero
    zy=zero
    call r8sincos(zth,imom,zsin(1:imom),zcos(1:imom))
    do imm=0,imom
      zr=zr+zrmc(imm,1)*zcos(imm)+zrmc(imm,2)*zsin(imm)
      zy=zy+zymc(imm,1)*zcos(imm)+zymc(imm,2)*zsin(imm)
    end do
    zrcon(ith)=zr
    zycon(ith)=zy
  end do
  !
  zrcon(inth)=zrcon(1)
  zycon(inth)=zycon(1)
  !
  call r8_getmpa_ry(zrcon,zycon,inth,zrin,zrout,zymp)
  !
  return
end subroutine r8_getmpa
!---------------------------------------------------------------------
!
subroutine r8_getmpa_ry(zrcon,zycon,inth,zrin,zrout,zymp)
  !
  !  find midplane elevation and R intercepts based on centroid, starting
  !  from a closed contour of (R,Y) pairs
  !
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero
  implicit none

  integer, intent(in) :: inth
  real(fp), intent(in) :: zrcon(inth),zycon(inth)  ! closed contour

  real(fp), intent(out) :: zrin,zrout    ! midplane intercepts (approx.)
  real(fp), intent(out) :: zymp          ! midplane elevation (centroid)
  !
  !  local:
  !
  real(fp) :: zrmp,zr,zy,zrp,zyp,zytest,zrans
  integer :: ith
  !
  !  find centroid
  !
  call r8_plcentr(zrcon,zycon,inth,zrmp,zymp)
  !
  !  find approximate midplane intercept locations
  !
  zr=zrcon(1)
  zy=zycon(1)
  !
  do ith=2,inth
    !
    zrp=zr
    zyp=zy
    zr=zrcon(ith)
    zy=zycon(ith)
    !
    zytest=(zy-zymp)*(zymp-zyp)
    if(zytest.ge.zero) then
      zrans=(zr*(zymp-zyp)+zrp*(zy-zymp))/(zy-zyp)
      if(zrans.gt.zrmp) then
        zrout=zrans
      else
        zrin=zrans
      end if
    end if
  end do
  !
  return
end subroutine r8_getmpa_ry
