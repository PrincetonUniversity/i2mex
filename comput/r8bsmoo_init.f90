subroutine r8bsmoo_init(imj,icentr,iedge,zsmoo,zxi)
  use iso_c_binding, only: fp => c_double
  use r8bsmoo_mod
  implicit none

  integer, intent(in) :: imj   ! grid array size
  integer, intent(in) :: icentr,iedge   ! active grid limits
  real(fp), intent(in) :: zsmoo    ! smoothing parameter
  real(fp), intent(in) :: zxi(imj,2)        ! the grid (TRCOM style)

  !  ***** initialize the module *****

  mj=imj
  lcentr=icentr
  ledge=iedge

  if(allocated(xi)) deallocate(xi)
  allocate(xi(mj,2))

  xi = zxi

  dxbsmoo = zsmoo

  lcp1 = lcentr+1
  lep1 = ledge+1

  nzones = ledge-lcentr + 1

  return
end subroutine r8bsmoo_init
 
subroutine r8bsmoo_get(zsmoo)
  use iso_c_binding, only: fp => c_double
  use r8bsmoo_mod

  implicit none

  real(fp), intent(out) :: zsmoo

  zsmoo = dxbsmoo

  return
end subroutine r8bsmoo_get
