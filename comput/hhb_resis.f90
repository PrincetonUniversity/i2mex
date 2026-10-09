subroutine hhb_resis(zneCM3i,zteEVi,zzeffi, &
     deltai,lpolCM,philimWB,rho,dVdrhoCM3,xiotb,vacRBphiTCM, &
     eta_spitzer,eta_nc,iwarn)

  ! 0d CGS implementation of resistivity calculation
  ! This is HHB resistivity as used in TRANSP, magcor/resis.for

  ! Reference:S. P. Hirshman, R. J. Hawryluk, B. Birge,
  ! Nucl. Fusion 17, 611 (1977). ("HHB")).

  ! All information passed through arguments; no modules or COMMON
  ! REAL*8 precision

  ! implemented in TRANSP/NTCC comput library
  ! reference to real(fp) function cloge_r8_fcn, also in comput library

  !----------------------------
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, one, twopi
  implicit none

  real(fp), parameter :: xmu = 1.25664E-6_fp

  !----------------------------
  ! arguments:

  real(fp) :: zneCM3i  ! electron density, cm**-3, iwarn=1 set if .le. 1 ptcl/cm3
  real(fp) :: zteEVi   ! electron temperature, eV, iwarn=1 set if .le. 1 eV
  real(fp) :: zZeffi   ! Zeff; iwarn=1 set if .le.1

  ! if iwarn=1 is set, a non-zero resistivity is still returned, but it is
  ! likely to be rubbish

  ! for the following, iwarn=2 is set and zero resistivity is returned, if
  ! any of these arguments are zero or negative

  real(fp) :: deltai   ! inverse aspect ratio r/R
  real(fp) :: LpolCM   ! poloidal path length, cm
  real(fp) :: philimWb ! enclosed toroidal flux at plasma boundary, Wb
  real(fp) :: rho      ! flux surface label: sqrt(PhiWB/PhilimWB) where
  ! phiWB is the enclosed toroidal flux at the rho surface
  real(fp) :: dVdrhoCM3    ! dVol/drho at the flux surface, cm**3
  real(fp) :: xiotb    ! iota(bar) = 1/q = dPsi/dPhi
  real(fp) :: vacRBphiTCM  ! toroidal field (R*B_phi) at vacuum boundary, T*cm

  ! output...

  real(fp), intent(out) :: eta_spitzer  ! Spitzer resistivity, Ohm*cm
  real(fp), intent(out) :: eta_nc   ! HHB neoclassical resistivity, Ohm*cm

  integer, intent(out) :: iwarn   ! warning flag, 0=normal

  !----------------------------
  !  local variables (several names taken from TRANSP's magcor/resis.for)

  real(fp) :: zneCM3,zteEV,zZeff

  real(fp) :: cloge_r8_fcn,cloge

  real(fp) :: zgamm,zconst

  real(fp) :: zbpav,zcuror

  real(fp) :: zvstae,zcr,zxi,zd1m,zft,zfmn,zfnc,zte32i

  integer :: jwarn

  !----------------------------

  iwarn = 0
  eta_spitzer = ZERO
  eta_nc = ZERO

  ! check plasma parameter inputs

  if(znecm3i.lt.ONE) then
    zneCM3=ONE
    iwarn=1
  else
    zneCM3=znecm3i
  endif

  if(zteevi.lt.ONE) then
    zteEV=ONE
    iwarn=1
  else
    zteEV=zteevi
  endif

  if(zzeffi.lt.ONE) then
    zZeff=ONE
    iwarn=1
  else
    zZeff=zzeffi
  endif

  ! check field and geometry inputs

  if(deltai.le.ZERO) iwarn=2
  if(LpolCM.le.ZERO) iwarn=2
  if(PhilimWb.le.ZERO) iwarn=2
  if(rho.le.ZERO) iwarn=2

  if(dVdrhoCM3.le.ZERO) iwarn=2
  if(xiotb.le.ZERO) iwarn=2
  if(vacRBphiTCM.le.ZERO) iwarn=2

  if(iwarn.eq.2) return

  !----------------------------------
  ! OK-- get electron Coulomb log

  cloge = cloge_r8_fcn(zneCM3,zteEV,zZeff,jwarn)

  !----------------------------------
  ! Spitzer resistivity

  !  for traditional (HHB) form of Spitzer resistivity
  ZGAMM=ZZEFF*0.581_fp*(2.67_fp+ZZEFF)/(1.13_fp+ZZEFF)
  zconst=5.22E-3_fp

  zte32i = ONE/(zteEV*sqrt(zteEV))

  eta_spitzer = zconst*cloge*ZGAMM*zte32i

  !----------------------------------
  ! NC resistivity 

  ! use Bpol volume average
  ! Bpol = (1/R)*grad(Psi) = (1/R)*grad(rho)*dPsi/drho
  !    dPsi/drho = xiotb*rho*PhilimWb/pi
  ! <Bpol> = [2pi*int(dl*R/grad(rho)*Bpol]/dVdrho
  !        =  2*Lpol*xiotb*rho*PhilimWb/dVdrho

  zBPAV = PhilimWB*2.E4_fp*rho*xiotb*LpolCM/dVdrhoCM3

  ! this formulation uses
  ! mu0*Ip = 2pi*r*<Bpol>; zcuror = Ip/r = 2pi*<Bpol>/mu0; convert to Amps/cm

  zcuror = TWOPI*zBPAV/(100.0_fp*XMU)

  ! code transcribed from magcor/resis.for

  ZVSTAE=3.46E-9_fp*zneCM3*cloge*vacRBPHITCM/(zcuror*zteEV*zteEV*sqrt(deltai))
  ZCR=0.56_fp*(3.0_fp-ZZEFF)/((3.0_fp+ZZEFF)*ZZEFF)
  ZXI=0.58_fp+0.2_fp*ZZEFF
  ZD1M=ONE-deltai
  ZFT=ONE-ZD1M*ZD1M/(SQRT(ONE-deltai*deltai)*(ONE+1.46_fp*SQRT(deltai)))
  ZFMN=ZFT/(ONE+ZXI*ZVSTAE)
  ZFNC=ONE/((ONE-ZFMN)*(ONE-ZCR*ZFMN))

  eta_nc = ZFNC*eta_spitzer

  return
end subroutine hhb_resis
