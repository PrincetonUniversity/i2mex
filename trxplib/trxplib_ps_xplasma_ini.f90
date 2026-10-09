subroutine trxplib_ps_xplasma_ini(ierr)

  ! init TRANSP MHDeq xplasma representation

  use trxplib_ps_options
  implicit NONE

  !------------------------------
  ! passed:

  integer, intent(out) :: ierr   ! exit code, 0=OK

  !------------------------------
  ! local:

  integer :: iwarn,ilun_wall,irzflag,ilparen,ilblank,lunzer,ifind
  integer :: ilt,ilz
  character*40 :: errstr,zrun_label
  character*10 :: zdate
  character*12 :: ztime
  character*4 :: ztok

  logical :: ltrx_freebdy ! test for free bdy data

  real*8 :: Psi0,Rmin,Rmax,Zmin,Zmax

  !------------------------------

  zdate=' '
  call c9date(zdate)
  call tget_rlbl(ztok,zrun_label)

  ilparen=index(zrun_label,'(')

  if(ilparen.gt.0) then
     runid=zrun_label(1:max(1,ilparen-1))
  else
     ilblank=index(zrun_label,' ')
     runid=zrun_label(ilblank+1:len(trim(zrun_label)))
  endif

  if(saw_hint.eq.-1) then
     call trx_chk_saw_r8(tselect,1)
  else if(saw_hint.eq.1) then
     call trx_chk_saw_r8(tselect,2)
  else
     call trx_chk_saw_r8(tselect,0)
  endif

  do
     errstr = "trx_time"
     call trx_time(tselect, delta_t, iwarn, ierr)
     if(ierr.ne.0) exit

     errstr = "trx_init_xplasma"
     call trx_init_xplasma(ierr)
     if(ierr.ne.0) exit

     call trx_set_threverse(0)

     !  get the core equilibrium

     errstr = "trx_mhd"
     call trx_mhd(n_theta, ierr)
     if(ierr.ne.0) exit

     !  get the limiter; determine (R,Z) range if necessary

     call find_io_unit(ilun_wall)
     if(ltrx_freebdy(0)) then

        errstr = "trx_wall_freebdy"
        call trx_wall_freebdy(ilun_wall,ierr)

     else

        Rmin = 0.0d0
        Rmax = 0.0d0

        Zmin = 0.0d0
        Zmax = 0.0d0

        errstr = "trx_wall_RZ"
        call trx_wall_RZ(ilun_wall,Rmin,Rmax,nR,Zmin,Zmax,nZ,ierr)

     endif
     if(ierr.ne.0) exit

     !  set up field extrapolation; specify signs...
     irzflag=1
     errstr = "trx_bxtr"
     call trx_bxtr(Bccw_hint,Jccw_hint,irzflag,ierr)
     if(ierr.ne.0) exit

     call trx_psi0(ifind,Psi0)
     if(ifind.eq.0) then
        geqdsk_lbl='TRXPL '//zdate//'*'
     else
        ! trxplib has Psi0 data only for "fb" free boundary cases
        geqdsk_lbl='TRXfb '//zdate//'*'
     endif

     ilt=len_trim(geqdsk_lbl)+2
     geqdsk_lbl(ilt:)=ztok

     ilt=len_trim(geqdsk_lbl)+2
     geqdsk_lbl(ilt:)=runid

     ilt=len_trim(geqdsk_lbl)+2
     geqdsk_lbl(ilt:)='t='

     ilt=len_trim(geqdsk_lbl)
     if(delta_t.gt.0) geqdsk_lbl(ilt:ilt)='~'
     ilz=min(12,(len(geqdsk_lbl)-ilt))

     write(ztime,'(f12.5)') tselect
     geqdsk_lbl(ilt+1:ilt+ilz)=ztime(1:ilz)

     exit
  enddo

  if(ierr.ne.0) then
     write(lunzer(0),*) ' ?trxplib_ps_xplasma_ini: error in subroutine: '//trim(errstr)
  endif

end subroutine trxplib_ps_xplasma_ini
