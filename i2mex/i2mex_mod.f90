module i2mex_mod
  !
  ! Standard toroidal equilibrium data interface.
  !
  ! OBJECT MEMBERS
  !
  ! pletzer@pppl.gov Thu Aug 17 10:36:31 EDT 2000
  !
  implicit none
 
  integer, parameter :: i2mex_clockwise = +1, i2mex_counterclockwise = -1
 
  integer, parameter  :: i2mex_r8 = selected_real_kind(12,100)
  integer, parameter  :: i2mex_r4 = selected_real_kind(6,37)
  real(i2mex_r8), parameter :: i2mex_pi_r8  = 3.1415926535897931_i2mex_r8
  real(i2mex_r8), parameter :: i2mex_twopi_r8  = 6.2831853071795865_i2mex_r8
  real(i2mex_r8), parameter :: i2mex_fourpi_r8 = 12.56637061435917_i2mex_r8
  real(i2mex_r4), parameter :: i2mex_pi_r4  = 3.1415926535897931_i2mex_r4
  real(i2mex_r4), parameter :: i2mex_twopi_r4  = 6.2831853071795865_i2mex_r4
  real(i2mex_r4), parameter :: i2mex_fourpi_r4 = 12.56637061435917_i2mex_r4
  integer, parameter :: i2mex_torSym=1, i2mex_No_autoT = 0
  integer, parameter :: i2mex_cubic=2, i2mex_akima=1, i2mex_linear=0
  integer, parameter :: i2mex_not_a_knot=0
  integer, parameter :: i2mex_bphi_ccw = +1, i2mex_bphi_cw = -1
  integer, parameter :: i2mex_Ip_ccw = +1, i2mex_Ip_cw = -1
  real(i2mex_r8), parameter :: i2mex_tol = 1.0e-4_i2mex_r8

  ! extrapolation to axis
  real(i2mex_r8) :: i2mex_to_axis_fraction = 0.02_i2mex_r8
  
  ! extrapolation to edge
  real(i2mex_r8) :: i2mex_to_edge_fraction = 0.01_i2mex_r8


  ! for direct equilibrium only
  integer :: i2mex_ns = 401
  integer :: i2mex_nt1 = 129
  real(i2mex_r8) :: i2mex_LAST_NORM_SURFACE_IS = 0.999_i2mex_r8
 
  ! The following can be used to remap the grid poloidal according
  ! to the dependence of the Jacobian ~ x^mx /(|grad psi|^npsi B^kb).
  ! remap=0 means no remapping performed.
  !
  integer :: i2mex_remap=0
  integer :: i2mex_mx=1, i2mex_npsi=1, i2mex_kb=0
 
  ! set this to 1 if direct representation BR, BZ, BPhi = BR, BZ, BPhi(R, Z)
  ! should be allowed
  integer :: i2mex_direct = 0
 
  ! no of points used to exrapolate q to axis (getting obsolete, use i2mex_to_axis_fraction instead)
  integer :: i2mex_nqextra = 0
 
  character*132 :: i2mex_comments=''

  ! set this to a value /=0 to shift the theta=0 ray
  integer :: i2mex_it0 = 0

  ! set this to 1 if J-prime is to be computed using Akima
  integer :: i2mex_xzakima = 0

  ! prevents ringing to occur for equilibria that are too fine
  ! by setting any of these different from zero.
  integer :: i2mex_nt1max = 0 ! 257
  integer :: i2mex_nsmax  = 0 !401
 
  type i2mex_obj
 
     ! Primitive toroidal plasma repesentation object as spline
     ! coefficients stored in an Xplasma buffer array. Only the
     ! indices into this buffer are stored for the radial
     ! position (id_s), the poloidal angle (id_t), the
     ! covariant toroidal B field (id_g), the safety factor
     ! (id_q) and the (X, Z)-grid (id_x and id_z).
 
     character(80) :: label
     integer :: nt1, ns
     real(i2mex_r8) :: psi_axis, psi_edge
     ! box boundaries for direct representation     
     real(i2mex_r8) :: rleft, rrigh, zbot, ztop 
     integer :: id_s, id_t, id_psi, id_p, id_g, id_q, id_x, id_z
     integer :: id_xakima, id_zakima, id_jakima
     integer :: id_Rgrid, id_Zgrid, id_psirz
     integer :: isThetaClockwise
 
  end type i2mex_obj
 
  type(i2mex_obj) :: i2mex_o
 
end module i2mex_mod
 
 
 
