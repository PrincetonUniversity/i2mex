subroutine r8_trirad(p1,p2,p3,p0,rad,ier)
  !
  !  given 3 verticies of a triangle (in a plane), return the center point
  !  and radius of a circle which passes through these points.
  !
  !  if the points are colinear, return ier=1 -- no circle with finite
  !  radius exists
  !
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero
  implicit none
  !
  real(fp), dimension(2), intent(in) :: p1,p2,p3    ! the verticies
  real(fp), dimension(2), intent(out) :: p0         ! center of the circle
  real(fp), intent(out) :: rad                      ! radius of circle
  !
  integer, intent(out) :: ier    ! completion code: 0=OK, 1=no circle exists
  !
  !-----------------------------------------------------
  !
  real(fp) :: delx21,delx31,dely21,dely31,avx21,avx31,avy21,avy31
  real(fp) :: pcross
  !
  !-----------------------------------------------------
  !
  p0 = zero
  rad = zero
  !
  delx21 = p2(1) - p1(1)
  delx31 = p3(1) - p1(1)
  dely21 = p2(2) - p1(2)
  dely31 = p3(2) - p1(2)
  !
  pcross=(delx21*dely31-dely21*delx31)
  if(pcross.eq.zero) then
    ier=1      ! colinear points
    return
  else
    ier=0
  endif
  !
  avx21 = (p2(1)+p1(1))/2
  avx31 = (p3(1)+p1(1))/2
  avy21 = (p2(2)+p1(2))/2
  avy31 = (p3(2)+p1(2))/2
  !
  p0(1) = (delx21*dely31*avx21 + dely21*dely31*avy21 &
          -delx31*dely21*avx31 - dely31*dely21*avy31)/pcross
  !
  p0(2) =-(delx21*delx31*avx21 + dely21*delx31*avy21 &
          -delx31*delx21*avx31 - dely31*delx21*avy31)/pcross
  !
  rad = sqrt((p1(1)-p0(1))**2+(p1(2)-p0(2))**2)
  !
  !debug  write(6,*) 'p0= ',p0
  !debug  write(6,*) 'rad=',rad
  !debug  write(6,*) 'rad=',sqrt((p2(1)-p0(1))**2+(p2(2)-p0(2))**2)
  !debug  write(6,*) 'rad=',sqrt((p3(1)-p0(1))**2+(p3(2)-p0(2))**2)
  !
  return
end subroutine r8_trirad
