subroutine rfx_antcen(nant,ra,za,rcen,zcen,hh)

  ! compute antenna center and half-height from array info
  ! hh = (zmax-zmin)/2; rcen = r value of location in ra(:),za(:) closest to
  !                            center of structure

  ! mod DMC May 2008: hh is now half the distance btw end pts; 
  !   (rcen,zcen) is the center -- allow off-midplane antenna.

  use iso_c_binding, only: fp => c_double
  implicit none

  !---------
  integer, intent(in) :: nant
  real(fp), intent(in) :: ra(nant),za(nant)
  real(fp), intent(out) :: rcen,zcen
  real(fp), intent(out) :: hh

  !---------
  real(fp) :: zcmin,zctest,zrr,zzz
  integer :: k,icmin

  !---------

  zcmin=(ra(nant)-ra(1))**2 + (za(nant)-za(1))**2
  hh = 0.5_fp * sqrt(zcmin)

  icmin=1

  do k=2,nant
    zctest = (ra(1)+ra(nant)-2.0_fp*ra(k))**2 + (za(1)+za(nant)-2.0_fp*za(k))**2
    if(zctest.lt.zcmin) then
      zcmin=zctest
      icmin=k
    endif
  enddo

  rcen=ra(icmin)
  zcen=za(icmin)

  return
end subroutine rfx_antcen
