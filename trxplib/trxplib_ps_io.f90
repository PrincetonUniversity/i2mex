!  ** trxplib_ps_io ** 
!  Control extraction of Plasma States from TRANSP archives
!
!  For many applications the following call sequence should serve:
!
!  (first call to connect to run -- be sure to check ierr return value):
!      call trxplib_ps_open(rdtype,rdpath,outpath,oprfix,ierr)
!         character*(*), intent(in) :: rdtype  !  "MDS+" or "FILE"
!         character*(*), intent(in) :: rdpath  ! path to run data
!         character*(*), intent(in) :: outpath ! path to directory for output
!         character*(*), intent(in) :: oprfix  ! PS file prefix (blank OK)
!         integer, intent(out) :: ierr       ! completion status, 0=OK
!
!            (this will write protimes.dat, sawtimes.dat, and a "test" state
!             including machine description, in trim(outpath)... ).
!
!  (multiple calls for each desired time_of_interest...)
!      call trxplib_ps_write(time_of_interest,ps_filename,ierr)
!         real*8, intent(in) :: time_of_interest    ! desired time
!         character*(*), intent(in) :: ps_filename  ! desired output filename
!         integer, intent(out) :: ierr       ! completion status, 0=OK
!
!    --- OR ---
!
!      call trxplib_ps_get(time_of_interest,<state_object>,ierr)
!         this returns a filled-in state object (no file is written)
!
!      there is no time averaging if these calls are used.
!
!  (final call to disconnect from run)
!      call trxplib_ps_close
!
!---------------------------------------------------
!  control routine, call any time:
!     trxplib_mhdeq_opt(...)  ! specify if output states should include
!                             ! MHD equilibrium flux surface geometry data
!
!     (see code, below, for explanation of calling arguments).
!---------------------------------------------------
!  additional routines, to be called PRIOR to trxplib_ps_open
!
!     trxplib_ps_reset  ! set defaults; if connected to run, disconnect now
!     trxplib_ps_nTheta(...)  ! set n_theta (generally default value is OK)
!     trxplib_ps_nRZ(...)  ! set nR and nZ (defaults OK: will be overridden
!                          ! in cases where TRANSP free boundary data exists).
!     trxplib_ccw_options(...)  ! CCW control options
!
!  Except for trxplib_ps_reset these routines have no effect on currently
!  opened run, if called after trxplib_ps_open(...)
!
!     (see code, below, for explanation of calling arguments).
!---------------------------------------------------
!  additional write routine:
!
!     trxplib_ps_write_opt(...)
!
!  works like trxplib_ps_write but with control over time averaging
!
!       --- OR ---
!
!     trxplib_ps_get_opt(...)
!
!  works like trxplib_ps_get but with control over time averaging
!
!     (see code, below, for explanation of calling arguments).
!---------------------------------------------------
!  The software may write messages-- to stdout (fortran unit 6) by default.
!  It shares message I/O control with "rplot" (TRANSP data access) libraries.
!  To redirect such messages, "rplot" library calls are used:
!
!     call plc_msgs(<lun>,filename) ...to open <filename> on integer unit <lun>
!                                   and write messages there, or,
!
!     call plc_msgs(<lun>,' ') ...to just write on integeru unit <lun>;
!                              presumably the caller has opened a file.


!--------------------------------------------------------------------------
subroutine trxplib_ps_open(rdtype,rdpath,outpath,oprfix,ierr)

  ! open connection to run...

  use plasma_state_mod
  use trxplib_ps_options
  implicit NONE

  character*(*), intent(in) :: rdtype  !  "MDS+" or "FILE"
  character*(*), intent(in) :: rdpath  ! path to run data
  character*(*), intent(in) :: outpath ! path to directory for output
  character*(*), intent(in) :: oprfix  ! PS time series output file prefix

  integer, intent(out) :: ierr         ! completion status, 0=OK

  !-------------------
  !  local:

  type (plasma_state) :: wkps

  integer :: nonlin,lunzer,ilun,ilen,istat,ierloc
  character*10 :: ztest,pidstr
  character*200 :: tfile

  !-------------------------------------------------
  ! executable code:

  ierr = 0

  nonlin = lunzer(0)  ! LUN for messages

  call find_io_unit(ilun)  ! LUN for output files
  call sget_pid_str(pidstr,ilen)  ! PID of current process

  !-------------------
  ! mark run as NOT open
  !-------------------

  ready = .FALSE.

  !-------------------
  ! test argument: rdtype
  !-------------------

  ztest = rdtype
  call uupper(ztest)

  if(ztest.eq.'MDS+') then
     lmds = .TRUE.
  else if(ztest.eq.'FILE') then
     lmds = .FALSE.
  else

     ierr=1
     write(nonlin,*) &
          ' ?trxplib_ps_open: 1st argument (rdtype) invalid or blank: '// &
          trim(rdtype)
     return
  endif

  !-------------------
  ! test argument: opath
  !   (have to be able to write a file here)
  !-------------------

  tfile = trim(outpath)//'/'//pidstr(1:ilen)//'.tmp'
  open(unit=ilun,file=tfile,status='unknown',iostat=istat)

  if(istat.ne.0) then
     ierr=1
     write(nonlin,*) &
          ' ?trxplib_ps_open: file open test unsuccessful in output directory:'
     write(nonlin,*) '   '//trim(outpath)
     return
  endif

  write(unit=ilun,fmt='(A)',iostat=istat) ' this is a test '

  if(istat.ne.0) then
     ierr=1
     write(nonlin,*) &
          ' ?trxplib_ps_open: file write test failed in output directory:'
     write(nonlin,*) '   '//trim(outpath)
  endif

  close(unit=ilun,status='delete')
  if(ierr.ne.0) return

  !-------------------
  !  OK...
  !-------------------

  rpath = rdpath
  opath = outpath

  call ps_init_user_state(wkps, "wkps_trxplib", ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_open: work state initialization failure.'
     return
  endif

  ps_prefix = oprfix

  call trxplib_ps_connect(wkps,ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_open: connect failure: '//trim(rdpath)
     return
  endif

  call ps_free_user_state(wkps, ierloc)
  if(ierloc.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_open: error deallocating work state (ignored).'
  endif

end subroutine trxplib_ps_open

!--------------------------------------------------------------------------
subroutine trxplib_ps_write(time,filename,ierr)

  ! write a Plasma State file -- no time averaging, no sawtooth hint

  use plasma_state_mod
  use trxplib_ps_options
  implicit NONE

  real*8, intent(in) :: time    ! time at which to extract data
  character*(*), intent(in) :: filename  ! root name of PS file
                                ! NOTE: ".cdf" appended, NetCDF file

  integer, intent(out) :: ierr  ! exit status code, 0=OK

  !---------------------------
  ! local:

  integer :: nonlin,lunzer
  integer :: ierloc
  type (plasma_state) :: wkps
  !---------------------------

  ierr = 0

  nonlin = lunzer(0)  ! LUN for messages

  saw_hint = 0
  delta_t = 0.0d0

  tselect = time
  ps_file = filename//'.cdf'

  call ps_init_user_state(wkps, "wkps_trxplib", ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write: work state initialization failure.'
     return
  endif

  call trxplib_ps_xplasma_ini(ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write: trxplib_ps_xplasma_ini failure.'
     return
  endif

  call trxplib_ps_write1(wkps, .FALSE., ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write: trxplib_ps_write1 failure.'
     return
  endif

  call ps_free_user_state(wkps, ierloc)
  if(ierloc.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write: error deallocating work state (ignored).'
  endif

end subroutine trxplib_ps_write

!--------------------------------------------------------------------------
subroutine trxplib_ps_get(time,ss,ierr)

  ! write a Plasma State file -- no time averaging, no sawtooth hint

  use plasma_state_mod
  use trxplib_ps_options
  implicit NONE

  real*8, intent(in) :: time    ! time at which to extract data
  type (plasma_state) :: ss     ! plasma state object, retuned with data
  !  CAUTION: any old data previously in (ss) is replaced.

  integer, intent(out) :: ierr  ! exit status code, 0=OK

  !---------------------------
  ! local:
  integer :: nonlin,lunzer
  !---------------------------

  ierr = 0

  nonlin = lunzer(0)  ! LUN for messages

  saw_hint = 0
  delta_t = 0.0d0

  tselect = time
  ps_file = 'NONE'

  call trxplib_ps_xplasma_ini(ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_get: trxplib_ps_xplasma_ini failure.'
     return
  endif

  call trxplib_ps_write1(ss, .FALSE., ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_get: trxplib_ps_write1 failure.'
     return
  endif

end subroutine trxplib_ps_get

!--------------------------------------------------------------------------
subroutine trxplib_ps_write_opt(time,filename,nhint_saw,delt_avg,ierr)

  ! write a Plasma State file -- no time averaging, no sawtooth hint

  use plasma_state_mod
  use trxplib_ps_options
  implicit NONE

  real*8, intent(in) :: time    ! time at which to extract data
  character*(*), intent(in) :: filename  ! root name of PS file
                                ! NOTE: ".cdf" appended, NetCDF file

  integer, intent(in) :: nhint_saw  ! sawtooth hint:
  !                             ! =1: look for post-sawtooth time(s)
  !                             ! =-1: look for pre-sawtooth time(s)
  !                             ! =0: no hint, use exact time provided

  real*8, intent(in) :: delt_avg  ! +/- averaging time for extracted data

  integer, intent(out) :: ierr  ! exit status code, 0=OK

  !---------------------------
  ! local:

  integer :: nonlin,lunzer
  integer :: ierloc
  type (plasma_state) :: wkps
  !---------------------------

  ierr = 0

  nonlin = lunzer(0)  ! LUN for messages

  saw_hint = nhint_saw
  saw_hint = max(-1, min(1, saw_hint))

  delta_t = abs(delt_avg)

  tselect = time
  ps_file = filename//'.cdf'

  call ps_init_user_state(wkps, "wkps_trxplib", ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write_opt: work state initialization failure.'
     return
  endif

  call trxplib_ps_xplasma_ini(ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write_opt: trxplib_ps_xplasma_ini failure.'
     return
  endif

  call trxplib_ps_write1(wkps, .FALSE., ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write_opt: trxplib_ps_write1 failure.'
     return
  endif

  call ps_free_user_state(wkps, ierloc)
  if(ierloc.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write_opt: error deallocating work state (ignored).'
  endif

end subroutine trxplib_ps_write_opt

!--------------------------------------------------------------------------
subroutine trxplib_ps_get_opt(time,ss,nhint_saw,delt_avg,ierr)

  ! write a Plasma State file -- no time averaging, no sawtooth hint

  use plasma_state_mod
  use trxplib_ps_options
  implicit NONE

  real*8, intent(in) :: time    ! time at which to extract data
  type (plasma_state) :: ss     ! plasma state object, retuned with data
  !  CAUTION: any old data previously in (ss) is replaced.

  integer, intent(in) :: nhint_saw  ! sawtooth hint:
  !                             ! =1: look for post-sawtooth time(s)
  !                             ! =-1: look for pre-sawtooth time(s)
  !                             ! =0: no hint, use exact time provided

  real*8, intent(in) :: delt_avg  ! +/- averaging time for extracted data

  integer, intent(out) :: ierr  ! exit status code, 0=OK

  !---------------------------
  ! local:
  integer :: nonlin,lunzer
  !---------------------------

  ierr = 0

  nonlin = lunzer(0)  ! LUN for messages

  saw_hint = nhint_saw
  saw_hint = max(-1, min(1, saw_hint))

  delta_t = abs(delt_avg)

  tselect = time
  ps_file = 'NONE'

  call trxplib_ps_xplasma_ini(ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_get_opt: trxplib_ps_xplasma_ini failure.'
     return
  endif

  call trxplib_ps_write1(ss, .FALSE., ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_get_opt: trxplib_ps_write1 failure.'
     return
  endif

end subroutine trxplib_ps_get_opt

!--------------------------------------------------------------------------
subroutine trxplib_ps_close

  ! disconnect from run (no arguments)

  use trxplib_ps_options
  implicit NONE

  !---------------------------

  ready = .FALSE.

end subroutine trxplib_ps_close

!--------------------------------------------------------------------------
subroutine trxplib_mhdeq_opt(option)

  ! specify TRUE for "heavy" state with MHD flux surfaces include
  ! otherwise specify FALSE

  ! if this routine is never called, TRUE is in effect

  use trxplib_ps_options
  implicit NONE

  logical, intent(in) :: option

  !---------------------------

  lheavy = option

end subroutine trxplib_mhdeq_opt

!--------------------------------------------------------------------------
subroutine trxplib_ps_reset

  ! reset I/O module defaults (no arguments)
  ! (this disconnects from any currently open run)

  use trxplib_ps_options
  implicit NONE

  !---------------------------

  call reset

end subroutine trxplib_ps_reset

!--------------------------------------------------------------------------
subroutine trxplib_ccw_options(Bccw,Jccw)

  ! control CCW options (CCW stands for "counter-clockwise")

  ! toroidal field:
  !   Bccw=1 means, force B_phi to point CCW in tokamak as viewed from above
  !   Bccw=-1 means force B_phi to point CW
  !   Bccw=0 means, use the orientation found in the TRANSP run data

  ! toroidal current:
  !   Jccw=1 means, force J_phi to point CCW in tokamak as viewed from above
  !   Jccw=-1 means force J_phi to point CW
  !   Jccw=0 means, use the orientation found in the TRANSP run data

  use trxplib_ps_options
  implicit NONE

  integer, intent(in) :: Bccw  ! B_phi orientation control
  integer, intent(in) :: Jccw  ! J_phi orientation control

  !---------------------------

  Bccw_hint = max(-1,min(1, Bccw))

  Jccw_hint = max(-1,min(1, Jccw))

end subroutine trxplib_ccw_options

!--------------------------------------------------------------------------
subroutine trxplib_ps_ntheta(nth)

  ! specify number of poloidal angle points, in MHD equilibrium flux surfaces
  ! numerical representation

  use trxplib_ps_options
  implicit NONE

  integer, intent(in) :: nth

  !---------------------------
  !  local:
  integer :: nonlin,lunzer
  !---------------------------

  nonlin = lunzer(0)

  n_theta = max(65,nth)
  
  if(n_theta.ne.nth) then
     write(nonlin,*) ' %trxplib_ps_ntheta: minimum value (65) imposed.'
     write(nonlin,*) '  input value was overridden: ',nth
  endif

end subroutine trxplib_ps_ntheta

!--------------------------------------------------------------------------
subroutine trxplib_ps_nRZ(nRi,nZi)

  ! specify number of (R,Z) cartesian grid points for
  ! numerical representation of free boundary MHD equilibria

  use trxplib_ps_options
  implicit NONE

  integer, intent(in) :: nRi,nZi   ! the R & Z grid sizes, respectively

  !---------------------------
  !  local:
  integer :: nonlin,lunzer
  !---------------------------

  nonlin = lunzer(0)

  nR = max(65,nRi)
  nZ = max(65,nZi)
  
  if((nR.ne.nRi).OR.(nZ.ne.nZi)) then
     write(nonlin,*) ' %trxplib_ps_nRZ: minimum value (65) imposed.'
     write(nonlin,*) '  input value(s) overridden: ',nR,nZ
  endif

end subroutine trxplib_ps_nRZ
