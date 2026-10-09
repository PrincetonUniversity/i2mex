subroutine trxplib_ps_connect(ss, ierr)

  ! initialize connection to TRANSP run & write test data
  ! use the working state (ss) as passed

  ! this routine should only be called internally within trxplib

  use plasma_state_mod
  use trxplib_ps_options
  implicit NONE

  !-----------------------------------
  ! arguments:

  type (plasma_state) :: ss     ! working state
  integer, intent(out) :: ierr  ! status code 0=normal

  !-----------------------------------
  ! local:

  character*250 :: tpath,spath

  integer :: lunzer,nonlin

  logical :: lheavy_save

  real*8 :: zt1,zt2
  real*8 :: Rmin,Rmax,Zmin,Zmax

  !-----------------------------------
  ! load all needed run data in memory...
  !   (rplot memory now uses dynamic expansion capability as an option)

  call datmgr_no_delete

  !-----------------------------------
  ! construct tpath

  ierr = 0
  nonlin = lunzer(0)

  lheavy_save = lheavy
  lheavy = .TRUE.

  if(lmds) then
     tpath = "MDS+:"//trim(rpath)
  else
     tpath = "FILE:"//trim(rpath)
  endif

  do 
     !  connect to run data...
     call trx_connect(tpath, ierr)
     if(ierr.ne.0) then
        write(nonlin,*) ' ?trxplib_ps_connect: trx_connect error.'
        exit
     endif

     !  discover time range...
     call trx_tlims(zt1,zt2, ierr)
     if(ierr.ne.0) exit

     !  write sawtooth times...
     spath = trim(opath)//'/sawtimes.dat'
     call trx_wr_stimes(spath, ierr)
     if(ierr.ne.0) exit

     !  write profile times...
     spath = trim(opath)//'/protimes.dat'
     call trx_wr_protimes(spath, ierr)
     if(ierr.ne.0) exit

     !  write test data: "heavy" state at 1st available time

     saw_hint = 0
     delta_t = 0.0d0

     tselect = zt1

     !  this routine also sets internal labels: geqdsk_lbl, runid
     call trxplib_ps_xplasma_ini(ierr)
     if(ierr.ne.0) then
        write(nonlin,*) ' ?trxplib_ps_connect: trxplib_ps_xplasma_ini error.'
        exit
     endif

     ps_file = trim(runid)//'_init.cdf'

     call trxplib_ps_write1(ss, .TRUE., ierr)
     if(ierr.ne.0) then
        write(nonlin,*) ' ?trxplib_ps_connect: _init state write error.'
        exit
     endif

     exit
  enddo

  lheavy = lheavy_save

end subroutine trxplib_ps_connect
