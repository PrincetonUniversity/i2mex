function get_j_from_jdotb(rbt,r2iav,b2av,riav,jdotb)
  !
  !  in axisymmetric toroidal system...
  !  convert from <j dot B> to <Jt>a
  !
  !     toroidal current = sum [ dArea * <Jt>a ]
  !
  !     < > means differential volume average over flux surface
  !     < >a means differential area average over flux surface
  !
  !  <Jt>a = <Jt/R>/<1/R> = (R*Bt)*<1/R**2>/(<B**2><1/R>) * <j dot B>
  !
  !  this routine is 0d -- does just one point location.
  !
  !  user's choice of units (but be consistent)
 
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) :: get_j_from_jdotb
  real(fp), intent(in) ::  rbt           ! R*Bt, e.g. Tesla*meter
  real(fp), intent(in) ::  r2iav         ! <1/R**2>, e.g. meter**-2
  real(fp), intent(in) ::  b2av          ! <B**2>, e.g. Tesla**2
  real(fp), intent(in) ::  riav          ! <1/R>, e.g. meter**-1
 
  real(fp), intent(in) ::  jdotb         ! <j dot B>, e.g. Tesla * Amps/meter**2
 
  !  output (function value):  <Jt>a (same units as j in jdotb)
  !------------------------------------------------------------
 
  get_j_from_jdotb = jdotb * (rbt*r2iav)/(b2av*riav)
 
  !------------------------------------------------------------
 
end function get_j_from_jdotb
 
function get_jdotb_from_j(rbt,r2iav,b2av,riav,jta)

  !
  !  in axisymmetric toroidal system...
  !  convert from <j dot B> to <Jt>a
  !
  !     toroidal current = sum [ dArea * <Jt>a ]
  !
  !     < > means differential volume average over flux surface
  !     < >a means differential area average over flux surface
  !
  !  <Jt>a = <Jt/R>/<1/R> = (R*Bt)*<1/R**2>/(<B**2><1/R>) * <j dot B>
  !
  !  this routine is 0d -- does just one point location.
  !
  !  user's choice of units (but be consistent)

  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) :: get_jdotb_from_j
  real(fp), intent(in) ::  rbt           ! R*Bt, e.g. Tesla*meter
  real(fp), intent(in) ::  r2iav         ! <1/R**2>, e.g. meter**-2
  real(fp), intent(in) ::  b2av          ! <B**2>, e.g. Tesla**2
  real(fp), intent(in) ::  riav          ! <1/R>, e.g. meter**-1

  real(fp), intent(in) ::  jta           ! <Jt>a, e.g. Amps/meter**2

  !  output (function value):  <J dot B>  ! units derived from above
  !------------------------------------------------------------

  get_jdotb_from_j = jta * (b2av*riav)/(rbt*r2iav)

  !------------------------------------------------------------

end function get_jdotb_from_j
