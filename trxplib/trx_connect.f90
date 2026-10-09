!
!  connect to a TRANSP run.
!
!---------------------------------------------------------------
subroutine trx_connect(path,ierr)
!
!  ** fortran interface **
!
  use trx_module
  implicit NONE
!
!  connect to a TRANSP run at specified "path:"
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  character(*), intent(IN) :: path   ! path / runid
  integer, intent(OUT) :: ierr       ! completion code, 0=OK
!
!    MDS+ path syntax:
!                 MDS+:<server-name>:<tree-name>(<shot-number>) or
!                 MDS+:<server-name>:<tree-name>(<tok.yy>,<runid>)
!
!    file path syntax:  <path>/<runid> or just <runid>
!
!--------------------------
!  local:
  character(64) :: zlbl
  character(32) :: zuns
  character(10) :: zxabb(8)
!
  integer imulti,istype,idims(10),irank,ilen
!
  integer lunzer
!
  integer, parameter :: maxa1=5001
  real*8 :: array1(maxa1)
  real*8 :: ztime
  real*8, parameter :: ZERO = 0.0d0
!
  logical :: iexist
!
!--------------------------
! ...initialize, then,
! ...use trread routine; set stats on run size...
!
  call trx_init
!
  call kconnect(path,run_label,nsctime,nprtime,nxmax,nmax, ierr)
  if(ierr.ne.0) then
     call trx_init
     return
  endif

  mds_arch_flag = (path(1:4).eq.'MDS+')

  if(.NOT.mds_arch_flag) then
     call ufilnam(path(6:),' ',file_path)
     !  remove trailing "/": path string ends in runid, not subdirectory
     ilen = len(trim(file_path))
     if(file_path(ilen:ilen).eq.'/') file_path(ilen:ilen)=' '
  else
     file_path = ' '
  endif

!
! get better labels...
!
  call tget_rlbl(tdev,run_label)
!
! OK... set tmin & tmax
!
  call trx_connect_tlims(ierr)
  if(ierr.ne.0) return
!
! OK... set no. of flux surfaces...
!
!
  call rplabel('XB',zlbl,zuns,imulti,istype)
  if(istype.le.0) then
     write(lunzer(0),*) ' ?trx_connect:  no "XB" profile in run:  ',run_label
     ierr=1
     return
  else
     itype_xb=istype
  endif
!
  call rpdims(istype,irank,idims,zxabb,ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) ' ?trx_connect:  "rpdims" error in run:  ',run_label
     ierr=1
     return
  endif
!
  nsurf=idims(1)+1
  call eqm_select('trx_connect',1)  ! clear xplasma
!
! also:  get X axis type id
!
  call rplabel('X',zlbl,zuns,imulti,istype)
  if(istype.le.0) then
     write(lunzer(0),*) ' ?trx_connect:  no "X" profile in run:  ',run_label
     ierr=1
     return
  else
     itype_x=istype
  endif
!
!-----------------------------------------
! check for free boundary data

  write(lunzer(0),*) ' %trx_connect: probing for free boundary data.'

  ztime = 0.5d0*(tmin+tmax)

  call rpexist_profile('RGRID',iexist)
  if(iexist) then
     call r8_t1profil('RGRID',zlbl,zuns,ztime,ZERO, &
          istype,array1,maxa1,nRfree,ierr)
     allocate(Rgrid_free(nRfree))
     Rgrid_free = array1(1:nRfree)*0.01_R8  ! -> m

     call r8_t1profil('ZGRID',zlbl,zuns,ztime,ZERO, &
          istype,array1,maxa1,nZfree,ierr)
     if(ierr.ne.0) then
        write(lunzer(0),*) &
             ' ?trx_connect: unexpected ZGRID read error after RGRID OK.'
        call trx_init
        ierr=1
        return
     endif

     allocate(Zgrid_free(nZfree))
     Zgrid_free = array1(1:nZfree)*0.01_R8  ! -> m
     
     allocate(PsiRZ_free(nRfree,nZfree))

     call r8_t1profil('PSIRZ',zlbl,zuns,ztime,ZERO, &
          istype,PsiRZ_free,nRfree*nZfree,irank,ierr)

     if ((PsiRZ_free(1,1).eq.PsiRZ_free(nRfree,1)).and.&
          (PsiRZ_free(1,1).eq.PsiRZ_free(1,nZfree)).and.&
          (PsiRZ_free(1,1).eq.PsiRZ_free(nRfree,nZfree))) then

        ierr = 0
        write(lunzer(0),*) ' (no real free boundary data).'

        deallocate(Rgrid_free, Zgrid_free, PsiRZ_free)
        nRfree = 0;  nZfree = 0
     else
        PsiRZ_free = ZERO

        ifound_psi0=0
        call rpexist_scalar('PSI0_TR',iexist); if(iexist) ifound_psi0=1
     endif
  else
     ierr=0
     write(lunzer(0),*) ' (no free boundary data).'
  endif

  return
end subroutine trx_connect
!--------------------------------------------
subroutine trx_connect_tlims(ierr)
!
!  fetch time limits
!
  use trx_module
  implicit NONE
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer ierr
!
! temporary array for sawtooth times
!
  integer iret
  integer lunzer
!
  logical :: iexist
  integer :: it,jt
!------------------------------
!
  ierr=0
!
  if(allocated(time_sc)) then
     deallocate(time_sc,time_pr,time_saw,kevent)
  endif

  allocate(time_sc(nsctime),time_saw(nsctime),time_pr(nprtime))
  allocate(kevent(nsctime)); kevent = 0

  call rptime_s(time_sc,nsctime,iret)
  if(iret.ne.nsctime) then
     write(lunzer(0),*) ' ?trx_connect:  internal inconsistency!'
     write(lunzer(0),*) '  nsctime = ',nsctime,' ...but, iret = ',iret
     ierr=1
     return
  endif
  tmin=time_sc(1)
  tmax=time_sc(nsctime)
!
  call rptime_p(time_pr,nprtime,iret)
  if(iret.ne.nprtime) then
     write(lunzer(0),*) ' ?trx_connect:  internal inconsistency!'
     write(lunzer(0),*) '  nprtime = ',nprtime,' ...but, iret = ',iret
     ierr=1
     return
  endif
  tmin=min(tmin,time_pr(1))
  tmax=max(tmax,time_pr(nprtime))
!
!  see if sawtooth times are available
  call rpexist_scalar('tlastsaw',iexist)
  jt = 0
  if(iexist) then
     call rpscalar('tlastsaw',time_saw,nsctime,iret,ierr)
     do it=2,nsctime
        if(time_saw(it).gt.time_saw(it-1)) then
           kevent(it)=2
           kevent(it-1)=1
           jt = jt + 1
        endif
     enddo
  endif
  if(jt.eq.0) then
     write(lunzer(0),*) ' %trx_connect_tlims: no sawteeth found.'
  else
     write(lunzer(0),*) ' %trx_connect_tlims: found ',jt,' sawteeth.'
  endif

  return
end subroutine trx_connect_tlims
!----------------------------------------------------------------------------
subroutine trx_init_xplasma(ierr)
!
!  read rho-axis (TRANSP XB data)
!
  use trx_module
  implicit NONE
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer, intent(OUT) :: ierr       ! completion code, 0=OK
!
  real zbuf(nsurf)  ! standard profiles
  real zbufr(nxmax) ! profiles vs. major radius
!
  real*8 zrho(nsurf+1),zdavg,zsq1,zsq2,zrho1,zrho2,zconv,ztime_r8
  real*8 zbufr8(nxmax)
!
  character(64) zlbl
  character(32) zuns
!
  integer istype,imax,igot,ieven,i,iwarn
  integer lunzer
!
!----------------------------
! re-initialize xplasma & trxplib units data
!
  units_status = -1
  xaxis_id = -1
  prof_units = ' '
  prof_iord = -1
!
! ...will read 2d axisymmetric transp mhd eq.
  call eqm_select('transp(trxplib) '//trim(run_label),1)
  ztime_r8 = time0
  call eqm_time(ztime_r8)
!
! ---> create grid of flux surfaces, augment TRANSP data w/ mag. axis
!
  imax=nsurf-1
  call t1profil('XB',zlbl,zuns,time0,delta_t, &
       & istype,zbuf,imax,igot,ierr)
  if(ierr.ne.0) return
!
  zrho(1)=0.0_R8
  zrho(2:nsurf-1)=zbuf(1:imax-1)   ! copy & convert to REAL*8
  zrho(nsurf)=1.0_R8               ! but make sure =1 at bdy
!
  call eqm_rho(zrho,nsurf,1.0e-4_R8,id_rho,ierr)
!
!  ---> create complementary grid of flux zone ctrs augmented with mag. axis
!  and plasma bdy
!
  call t1profil('X',zlbl,zuns,time0,delta_t, &
       & istype,zbuf,imax,igot,ierr)
  if(ierr.ne.0) return
!
  zrho(1)=0.0_R8
  zrho(2:nsurf)=zbuf(1:imax)   ! copy & convert to REAL*8
  zrho(nsurf+1)=1.0_R8
!
  call eqm_uaxis('rho_zc',id_rho,0,zrho,nsurf+1,1.0e-6_R8,id_rhozc,ierr)
!
! ---> extrapolation factors
!
  ieven=1
  zdavg=(zrho(nsurf)-zrho(1))/imax
  do i=1,imax
     if(abs(zdavg-(zrho(i+1)-zrho(i))).gt.(1.0E-4_R8*imax*zdavg)) ieven=0
  enddo
  if(ieven.eq.1) then
     afac0b=1.0E0_R8/3.0E0_R8
     afaclinb=1.0E0_R8
     afac0=1.0E0_R8/8.0E0_R8
     afaclin=0.5E0_R8
     efac0=afac0
     efaclin=afaclin
  else
     zsq1=(zrho(2)-zrho(1))**2
     zsq2=(zrho(3)-zrho(1))**2
     afac0b=zsq1/(zsq2-zsq1)
     afaclinb=(zrho(2)-zrho(1))/(zrho(3)-zrho(2))
!
     zrho1=0.5_R8*(zrho(1)+zrho(2))
     zrho2=0.5_R8*(zrho(2)+zrho(3))
     zsq1=(zrho1-zrho(1))**2
     zsq2=(zrho2-zrho(1))**2
     afac0=zsq1/(zsq2-zsq1)
     afaclin=(zrho1-zrho(1))/(zrho2-zrho1)
!
     zrho1=0.5_R8*(zrho(nsurf-1)+zrho(nsurf))
     zrho2=0.5_R8*(zrho(nsurf-2)+zrho(nsurf-1))
     zsq1=(zrho(nsurf)-zrho1)**2
     zsq2=(zrho(nsurf)-zrho2)**2
     efac0=zsq1/(zsq2-zsq1)
     efaclin=(zrho(nsurf)-zrho1)/(zrho1-zrho2)
  endif
!
! ---> major radius grids (optional): standard grid, TRDAT grid
!
  call t1profil('RMAJM',zlbl,zuns,time0,delta_t, &
       & istype,zbufR,nxmax,igot,ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) &
          & 'trx_init_xplasma warning: RMAJM (major radius) grid: not found.'
     ierr=0  ! OK
     return
  endif
  call trx_mks_conv('RMAJM',zconv,zuns,iwarn)
!
  zbufR8(1:igot)=zconv*zbufR(1:igot)
  call eqm_uaxis('Rmajm',0,0,zbufR8,igot, &
       & 1.0e-6_R8*zbufR8(igot),id_rmajm,ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) &
          & 'trx_init_xplasma warning: eqm_uaxis error on RMAJM.'
     ierr=0
     return
  endif
  itype_rmajm=istype
  nrmajm=igot
!
! TRDAT grid
!
  call t1profil('RMJSYM',zlbl,zuns,time0,delta_t, &
       & istype,zbufR,nxmax,igot,ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) &
          & 'trx_init_xplasma warning: RMJSYM (major radius) grid: not found.'
     ierr=0  ! OK
     return
  endif
  call trx_mks_conv('RMJSYM',zconv,zuns,iwarn)
!
  zbufR8(1:igot)=zconv*zbufR(1:igot)
  call eqm_uaxis('Rmjsym',0,0,zbufR8,igot, &
       & 1.0e-6_R8*zbufR8(igot),id_rmjsym,ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) &
          & 'trx_init_xplasma warning: eqm_uaxis error on RMJSYM.'
     ierr=0
     return
  endif
  itype_rmjsym=istype
  nrmjsym=igot
 
  return
end subroutine trx_init_xplasma
