subroutine trxplib_ps_write1(ss, ini_flag, ierr)

  ! write a state file  
  ! use the working state (ss) as passed

  ! this routine should only be called internally within trxplib

  use plasma_state_mod
  use trxplib_ps_options
  implicit NONE

  !-----------------------------------
  ! arguments:

  type (plasma_state) :: ss     ! working state
  logical, intent(in) :: ini_flag  ! .TRUE. if called from trxplib_ps_connect
  integer, intent(out) :: ierr  ! status code 0=normal

  !-----------------------------------
  ! local:

  logical :: out_flag
  integer :: lunzer,nonlin,ilun,iws
  character*250 :: spath
  !-----------------------------------

  ierr = 0

  nonlin = lunzer(0)

  if(lheavy) then
     ilun=1
  else
     ilun=0
  endif

  call trx_init_state_obj(ss,ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?trxplib_ps_write1: trx_init_state_obj error.'
     return
  endif

  out_flag = (ps_file.ne."NONE")

  if(ini_flag) then
     iws=2
  else
     if(out_flag) then
        iws=1
     else
        iws=0
        ilun=0
     endif
  endif

  spath = trim(opath)//'/'//trim(ps_file)

  if(ilun.eq.1) call find_io_unit(ilun)

  call trx_gen_state_obj(ss, ilun, geqdsk_lbl, spath, iws, ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' trxplib_ps_write1: trx_gen_state_obj error detected.'
  endif

end subroutine trxplib_ps_write1
