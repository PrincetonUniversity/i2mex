module trxplib_ps_options

  implicit NONE
  SAVE

  ! module to hold data options for accessing Plasma State data from
  ! TRANSP archives

  !---------------------------------
  ! path information

  logical :: lmds = .TRUE.  ! .TRUE. to use MDSplus server
  !                         ! ...or .FALSE. to use files

  character*200 :: rpath = ' '   ! path to run:
  !  if(lmds) syntax is <server>:<treename>(<run-ident>)
  !  if(.NOT.lmds) syntax is <directory-path>/<runid> or just <runid>

  character*200 :: opath = ' '   ! directory path for output file data
  !  (blank for current working directory

  character*100 :: ps_file = ' ' ! Plasma State filename

  character*10 :: ps_prefix= ' ' ! prefix for time series PS output

  !---------------------------------
  ! time slice selection information
  real*8 :: tselect = 0.0d0      ! time of interest

  integer :: saw_hint = 0        ! -1: pre-sawtooth; +1: post; 0: no hint
  real*8 :: delta_t = 0.0d0      ! +/- averaging time
  !  time average is over range [tselect-deltat:tselect+deltat] but without
  !  crossing sawtooth event boundary

  !---------------------------------
  ! equilibrium data representation & other options

  integer :: n_theta = 151       ! #poloidal angle pts in flux surfaces

  !  the following will be reset to match data, if the TRANSP run accessed
  !  has free boundary information

  integer :: nR = 101            ! #pts, R axis, (R,Z) overlay grid
  integer :: nZ = 101            ! #pts, Z axis, (R,Z) overlay grid

  integer :: Bccw_hint = 0       ! +1: force B_phi ccw viewed from above;
  !                              ! -1: force B_phi cw viewed from above
  !    leave at zero to use the orientation data saved in the archived run data

  integer :: Jccw_hint = 0       ! +1: force J_phi ccw viewed from above;
  !                              ! -1: force J_phi cw viewed from above
  !    leave at zero to use the orientation data saved in the archived run data

  !---------------------------------
  ! Plasma State options

  logical :: lheavy=.TRUE.   ! .TRUE. to output "heavy-weight" states (incl 2d 
                      ! MHD equilibrium data; .FALSE. to output "light-weight"
                      ! states

  !---------------------------------
  ! DMC Nov 2010
  ! auxilliary machine description files (max 100)
  ! contents are made available in "ps_aux" during trx_gen_state_geq execution
  ! the machine description files are read in trx_init_state_obj --
  ! see trxplib/trx_gen_state.f90

  integer, parameter :: max_naux = 100
  integer :: naux_mdescr = 0
  character*150 aux_mdescr(max_naux)
  !---------------------------------

  ! internal status

  logical :: ready = .FALSE.

  !---------------------------------
  ! internal labeling data

  character*48 :: geqdsk_lbl
  character*20 :: runid

  !-------------------------------------------------------------------------
CONTAINS

  subroutine reset
    ! restore default settings

    lmds = .TRUE.
    rpath = ' '
    opath = ' '
    ps_file = ' '

    tselect = 0.0d0
    saw_hint = 0
    delta_t = 0.0d0

    n_theta = 151
    nR = 101
    nZ = 101
    Bccw_hint = 0
    Jccw_hint = 0

    lheavy = .TRUE.

    ready = .FALSE.

  end subroutine reset

end module trxplib_ps_options
