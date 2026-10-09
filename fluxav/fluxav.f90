module fluxav
!
  implicit NONE
  SAVE
!
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!  dmc 8 Sep 2000
!  support module for flux surface & flux zone averaging and integration
!  in context of a toroidal system
!
!  dmc 15 Aug 2006 -- module simplified; reduced role in xplasma2
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!
  integer :: lunerr  = 6        ! unit for error messages
!
!-------------------------------------------------------------------
! theta and phi coordinate limits (should differ by 2pi)
!
  real*8, parameter :: twopi = 6.2831853071795862_R8
 
  real*8 :: theta_min = 0.0_R8       ! minimum poloidal angle
  real*8 :: theta_max = twopi        ! max poloidal angle
 
  real*8 :: phi_min =   0.0_R8       ! minimum toroidal angle
  real*8 :: phi_max =   twopi        ! max toroidal angle
!
!-------------------------------------------------------------------
! break-point lists, in each dimension...
!
  real*8, dimension(:), pointer :: rhobrk,rhobrk_tmp
  real*8, dimension(:), pointer :: thbrk,thbrk_tmp
  real*8, dimension(:), pointer :: phibrk,phibrk_tmp
!
! no. of break-points used, each dimension...
!
  integer :: n_rhobrk = 0
  integer :: n_thbrk = 0
  integer :: n_phibrk = 0
!
! size of (expandable) allocated list array, each dimension...
!
  integer :: n_rhomax = 0
  integer :: n_thmax = 0
  integer :: n_phimax = 0
!
!-------------------------------------------------------------------
end module fluxav
