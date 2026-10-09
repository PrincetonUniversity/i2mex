subroutine rf_egridtr(eminev,zdemin,egrmax,maxen,enerev)
  !
  !  real(fp) interface
  !
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, one
  implicit none
  real(fp), intent(in) :: eminev      ! min energy (eV) (usually zero)
  real(fp), intent(in) :: zdemin      ! energy step at eminev
  real(fp), intent(in) :: egrmax      ! max energy (eV)
  integer, intent(in) :: maxen      ! no. of grid zones
  real(fp), intent(out) :: enerev(0:maxen) ! energy zone *bdys*
  !
  !  note 0:maxen indexing of enerev!
  !
  !------------------------------------------
  !  for TRANSP (dmc Feb 96) control the max energy of the non-uniform
  !  grid -- mod dmc Jan 1997)
  !
  !  the energy grid is linearly spaced in steps of (zdemin) until the
  !  first step of a logarithmic spacing from the top of the linearly
  !  spaced region is larger than the linear step (zdemin) would be.
  !
  !------------------------------------------
  !
  real(fp) :: zEbase,zErat,deove,zlebase,zlemax,zle
  integer i,i0,ileft
  !
  !------------------------------------------
  !
  enerev(0)=max(ZERO,eminev)
  if(enerev(0).eq.ZERO) then
    enerev(1)=zdemin
    i0=1
  else
    i0=0
  end if

  !  linear part-- until log spacing yields bigger step...

  do
    zEbase=enerev(i0)
    zErat=egrmax/zEbase
    ileft=maxen-i0
    deove=exp(log(zErat)/ileft)
    if(zEbase*(deove-ONE).lt.zdemin) then
      enerev(i0+1)=enerev(i0)+zdemin
      i0=i0+1
    else
      exit
    end if
  end do

  zlebase = log(zEbase)
  zlemax = log(egrmax)
  do i=i0+1,maxen
    zle=(zlebase*(maxen-i)+zlemax*(i-i0))/(maxen-i0)
    enerev(i)=exp(zle)
  end do
!
  return
end subroutine rf_egridtr

!-----------------------------
subroutine rf_egridtr_real(eminev,zdemin,egrmax,maxen,enerev)
  !
  !  REAL interface
  !
  use iso_c_binding, only: dp => c_double, sp => c_float
  implicit none
  real(sp), intent(in) :: eminev        ! min energy (eV) (usually zero)
  real(sp), intent(in) :: zdemin        ! energy step at eminev
  real(sp), intent(in) :: egrmax        ! max energy (eV)
  integer, intent(in) :: maxen      ! no. of grid zones
  real(sp), intent(out) :: enerev(0:maxen) ! energy zone *bdys*
  !
  real(dp) :: zeminev,zegrmax,zzdemin
  real(dp) :: zenerev(0:maxen)
  !
  zeminev=eminev
  zzdemin=zdemin
  zegrmax=egrmax
  !
  call rf_egridtr(zeminev,zzdemin,zegrmax,maxen,zenerev)
  !
  enerev=zenerev
  !
  return
end subroutine rf_egridtr_real
