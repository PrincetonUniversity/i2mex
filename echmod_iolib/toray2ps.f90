subroutine toray2ps(filenam,powech,curfac,iant,ss,iout,ierr)

  ! extract data from one TORAY run, load Plasma State output profiles
  ! and map TORAY output from netCDF file to the Plasma State.

  use plasma_state_mod
  use ezcdf

  implicit NONE

  character*(*), intent(in) :: filenam  ! filename, "toray.nc" used if blank

  real*8, intent(in) :: powech   ! ECH injected power, this antenna/run
  real*8, intent(in) :: curfac   ! anomaly factor for current, =1.0 usually

  integer, intent(in) :: iant    ! antenna index

  type (plasma_state) :: ss      ! Plasma State object

  integer, intent(in) :: iout    ! I/O unit for error messages

  integer, intent(out) :: ierr   ! completion status, 0=OK

  !---------------------
  ! local:

  character*100 zfilenam

  integer :: ncid,dimlens(3),iz,inz,inz_exp
  character*5 xtype

  real*8, dimension(:), allocatable :: pei_data,curi_data,xdata,xtest,tmp_prof
  real*8 :: xtol,zcp,zpp,zincr

  !---------------------

  ierr = 0

  if(filenam.eq.' ') then
     zfilenam = "toray.nc"
  else
     zfilenam = filenam
  endif

  !---------------------
  ! error check...

  if(ss%nrho_ecrf.eq.0) then
     write(iout,*) ' ?toray2ps: ECRF grid has not been set up.'
     ierr=1
  endif

  if(ss%necrf_src.eq.0) then
     write(iout,*) ' ?toray2ps: no ECRF sources in Plasma State (ss).'
     ierr=1
  else if((iant.lt.1).or.(iant.gt.ss%necrf_src)) then
     write(iout,*) ' ?toray2ps: Antenna index: ',iant
     write(iout,*) '  not in expected range: 1 to necrf_src = ',ss%necrf_src
     ierr=1
  endif

  if(.not.allocated(ss%peech_src)) then
     write(iout,*) ' ?toray2ps: Plasma State ECH profiles not allocated.'
     ierr = 1
  endif

  if(ierr.ne.0) return

  ss%peech_src(:,iant) = 0.0d0
  ss%curech_src(:,iant) = 0.0d0

  !---------------------
  inz_exp = ss%nrho_ecrf - 1   ! expected number of bins

  !---------------------
  ! open the file...

  call cdf_open(ncid,trim(zfilenam),'r',ierr)
  if(ierr.ne.0) then
     write(iout,*) ' ?toray2ps: NetCDF open failure: '//trim(zfilenam)
     return
  endif

  !---------------------
  ! read the profile size...

  call cdf_inquire(ncid,'xmrho',dimlens,xtype)
  inz = dimlens(1)

  !---------------------
  ! allocate buffers...

  allocate(pei_data(inz),curi_data(inz),xdata(inz),xtest(inz+1))

  !---------------------
  ! read data & close file...

  call cdf_read(ncid,'xmrho',xdata)
  call cdf_read(ncid,'tpowde',pei_data)
  call cdf_read(ncid,'tidept',curi_data)

  call cdf_close(ncid)

  !---------------------
  ! grid  for re-zoning
  xtest(1)=0.0d0
  xtest(2:inz) = 0.5d0*(xdata(1:inz-1)+xdata(2:inz))
  xtest(inz+1)=1.0d0

  !---------------------
  ! OK: load the data: W/bin, A/bin
  allocate(tmp_prof(inz))
  zpp=0.0d0
  zcp=0.0d0
  tmp_prof=0.0d0

  do iz=1,inz
     zincr = pei_data(iz) - zpp
     tmp_prof(iz)= zincr*powech
     zpp = pei_data(iz)
  enddo
  call ps_user_rezone1(xtest, ss%rho_ecrf, tmp_prof,&
       ss%peech_src(1:inz_exp,iant), ierr, state=ss)

  tmp_prof=0.0d0
  do iz=1,inz
     zincr = curi_data(iz) - zcp
     tmp_prof(iz) = zincr*powech*curfac
     zcp = curi_data(iz)
  enddo
  call ps_user_rezone1(xtest, ss%rho_ecrf, tmp_prof,&
       ss%curech_src(1:inz_exp,iant), ierr, state=ss)
  !
  !---------------------
  !  all done.


  deallocate(pei_data,curi_data,xdata,xtest,tmp_prof)

end subroutine toray2ps
