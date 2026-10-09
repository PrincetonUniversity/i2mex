module trx_module
 
!  module "COMMON" for trx library -- TRANSP data access and interpolation
!
!-----------------------
!  mod DMC Feb 2008 -- use trdatbuf_module
!  create ability to access TRDAT (TRANSP input) datasets
!-----------------------

  use trdatbuf_module
  implicit NONE
  SAVE

!  New Feb 2008:

  Type (trdatbuf) :: d
  logical :: d_init_flag = .FALSE.
  logical :: d_data_avail = .FALSE.

  logical :: symZ = .TRUE.   ! .true. to symmetrize the grid limits in Z.  This causes problems
                             ! when using a free boundary grid which is not centered in Z.
!-------------------------------
!  original trx_module variables follow...
!
!  basic statistics
!
  integer :: nsctime         ! number of points in scalar timebase
  integer :: nprtime         ! number of points in profile timebase
  integer :: nxmax           ! max no. of points, any x axis
  integer :: nmax            ! max no. of points,
  integer :: nsurf           ! no. of flux surfaces, including mag. axis
  integer :: nrmajm          ! no. of major radius points (RMAJM)
  integer :: nrmjsym         ! no. of major radius points (RMJSYM)
!
  integer :: ksym            ! =0: updown symmetric; =1: updown asymmetric
  integer :: nmoms           ! no. of moments, TRANSP equilibrium.
!
!  profile type identification
!
  integer :: itype_x         ! type code:  zone centered profiles
  integer :: itype_xb        ! type code:  bdy centered profiles
  integer :: itype_Rmajm     ! type code:  vs. RMAJM
  integer :: itype_Rmjsym    ! type code:  vs. RMJSYM
!
!  free boundary information
!
  integer :: nRfree,nZfree
  integer :: ifound_Psi0=0
  real*8, dimension(:), allocatable :: Rgrid_free,Zgrid_free
  real*8, dimension(:,:), allocatable :: PsiRZ_free
  real*8 :: raxis_mhd,zaxis_mhd
  real*8 :: Psi0_mhd = 0.0d0
!
!  grid identification -- (rho,chi), (R,Z)
!
  integer :: id_rho          ! flux surface grid
  integer :: id_rhozc        ! flux zone centered grid augmented w/ axis & bdy
  integer :: id_chi          ! poloidal angle grid
  integer :: id_R            ! R grid
  integer :: id_Z            ! Z grid
  integer :: id_Rmajm        ! TRANSP RMAJM grid
  integer :: id_Rmjsym       ! TRANSP RMJSYM grid
!
!  (the following is used to enable reversing direction of increasing
!   theta given by TRANSP-- for software testing purposes)
!
  integer :: th_reverse = 0  ! set =1 to reverse theta order
!
  integer :: isign_rzbc = 1  ! set =-1 to get finite difference BC applied
                             ! to R(theta,rho) and Z(theta,rho) profiles.
!
!  extrapolation factors
!
  real*8 afac0,afaclin       ! zone 1 ctr -> axis
  real*8 efac0,efaclin       ! zone N ctr -> edge
  real*8 afac0b,afaclinb     ! surface 1 -> axis
!
!  controls for spline boundary conditions
!
  integer :: ibc_axis        ! axial boundary conditions
  real*8  :: zbc_axis        ! axial BC data item
  integer :: ibc_edge        ! edge boundary condition
  real*8  :: zbc_edge        ! edge BC data item
!
!  extrapolation over-ride controls
!
  integer :: ival_axis       ! axial value over-ride flag
  real*8  :: zval_axis       ! axial over-ride value
  integer :: ival_edge       ! edge value over-ride flag
  real*8  :: zval_edge       ! edge over-ride value
!
!  extrapolation limit controls --
!    don't let extrapolation change sign or get too close to zero
!
  integer :: isign_check     ! check that extrapolation doesn't change sign
  real*8  :: zsign_lim       ! zsign_lim*(edge value) = nearest approach to 0.0
!
!  edge value means non-extrapolated value nearest to the edge.
!  -------------
!  TRANSP run id
!
  character(40) :: run_label
!
!  tokamak id
!
  character(4) :: tdev
!
!  time information
!
  real :: tmin               ! minimum time in run
  real :: tmax               ! maximum time in run
!
  real, dimension(:), allocatable :: time_sc,time_pr,time_saw
  integer, dimension(:), allocatable :: kevent ! (indexed on scalar time grid)
  !  kevent(j) = 1 means: time_sc(j) at start of sawtooth;
  !  kevent(j) = 2 means: time_sc(j) at end of sawtooth; event at time_saw(j).
  !  kevent(j) = 0 means: between sawteeth, or, unknown

  integer :: itset           ! flag that slice time has been set (=1 if so)
  real :: time0              ! time of interest
  real :: delta_t            ! +/- delta_t, time to average over
!
! units information
!
  integer, parameter :: maxid=4096
!
  integer xaxis_id(maxid)          ! xplasma x-axis id for 1d objects
  character(10) prof_units(maxid)  ! units label, vs. xplasma id#
  integer units_status(maxid)      ! status code
!   -1 = undefined
!    0 = OK, MKS conversion succeeded
!    1 = MKS conversion failed, original TRANSP units & label are used.
  integer prof_iord(maxid)         ! interp. order -- e.g. 1 for Hermite
!
!--------------------------------------
! TRANSP data file access information

  logical :: mds_arch_flag  ! =.TRUE. -- using MDS+ TRANSP data
  character*120 file_path   ! path to TRANSP files (if mds_sock_id = 0)

end module trx_module
