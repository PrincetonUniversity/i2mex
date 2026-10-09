!
! --------- alphabeta_to_thetaphi----------
! This is the formula that should be used for ITER. NOTE: called in expert files
! NOTE: consider deleting this
!
subroutine alphabeta_to_thetaphi(zalpha, zbeta, ztheta, zphi)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: rad2deg
  implicit none

  real(fp), intent(in)  :: zalpha, zbeta  ! xpoldeg,xtordeg as input to namelist
  real(fp), intent(out) :: ztheta, zphi   ! xpoldeg,xtordeg as input to floatinbeam(2),floatinbeam(1) 

  ztheta = rad2deg * asin(cos(zbeta/rad2deg)*sin(zalpha/rad2deg))
  zphi   = rad2deg * (-atan(tan(zbeta/rad2deg)/cos(zalpha/rad2deg)))

  return
end subroutine alphabeta_to_thetaphi

!
! --------- thetaphi_to_alphabeta ----------
! This does the inverse of the conversion of xpoldeg,xtordeg as done for ITER
! So it allows constructing the alpha,beta needed given the desired theta, phi.
! NOTE: consider deleting this
!
subroutine thetaphi_to_alphabeta(ztheta, zphi, zalpha, zbeta)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: rad2deg, zero, one, two
  implicit none

  real(fp), intent(in)  :: ztheta, zphi   ! xpoldeg,xtordeg as input to floatinbeam(2),floatinbeam(1) 
  real(fp), intent(out) :: zalpha, zbeta  ! xpoldeg,xtordeg as input to namelist

  zalpha = rad2deg * acos(cos(ztheta/rad2deg)/sqrt(one+(tan(zphi/rad2deg)*sin(ztheta/rad2deg))**two))
  if (ztheta<zero) zalpha=-zalpha

  zbeta  = rad2deg * (-atan(cos(zalpha/rad2deg)*tan(zphi/rad2deg)))

  return
end subroutine thetaphi_to_alphabeta

!
! ---------- psangle_to_thetaphi ---------
! go from plasma state aiming angles to theta and phi desired by torbeam
!
subroutine psangle_to_thetaphi(ztheta_ec, zphi_ec, ztheta, zphi)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only : pi
  implicit none

  real(fp), intent(in)  :: ztheta_ec, zphi_ec  ! plasma state aiming angles
  real(fp), intent(out) :: ztheta, zphi        ! xpoldeg,xtordeg as input to floatinbeam(2),floatinbeam(1)

  real(fp), parameter :: eps = 1.0e-10_fp
  real(fp)::alpha,beta

  ztheta = ztheta_ec - 90._fp
  zphi   = -zphi_ec  + 180._fp
  if (zphi>270._fp) zphi=zphi-360._fp

end subroutine psangle_to_thetaphi

!
! ---------- thetaphi_to_psangle ----------
! go from torbeam theta and phi to plasma state aiming angles
subroutine thetaphi_to_psangle(ztheta, zphi, ztheta_ec, zphi_ec)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only : pi
  implicit none

  real(fp), intent(in)  :: ztheta, zphi        ! xpoldeg,xtordeg as input to floatinbeam(2),floatinbeam(1)
  real(fp), intent(out) :: ztheta_ec, zphi_ec  ! plasma state aiming angles
  real(fp), parameter :: eps = 1.0e-10_fp
  real(fp)::alpha,beta

  ztheta_ec = ztheta + 90._fp
  zphi_ec   = -zphi  +180._fp
  if (zphi_ec>270._fp) zphi_ec=zphi_ec-360._fp

  return
end subroutine thetaphi_to_psangle
