module physconst_mod

  use iso_c_binding, only: fp => c_double

  implicit none

  !-----------------------------------------
  real(fp), parameter :: zero = 0.0_fp
  real(fp), parameter :: half = 0.5_fp
  real(fp), parameter :: one  = 1.0_fp
  real(fp), parameter :: two  = 2.0_fp
  !-----------------------------------------

  real(fp), parameter :: pi   = atan2(0.0_fp,-1.0_fp)
  real(fp), parameter :: twopi= PI+PI
  real(fp), parameter :: pio2 = half*pi
  real(fp), parameter :: pio4 = half*pio2
  real(fp), parameter :: sqpi = sqrt(pi)
  real(fp), parameter :: kb_K = 1.38064852e-23_fp    !  J/K
  real(fp), parameter :: kb_eV = 1.60217653e-19_fp   !  J/eV
  real(fp), parameter :: kb_keV = 1.60217653e-16_fp  !  J/keV
  real(fp), parameter :: qe_eV  = 1.60217653e-19_fp  !  qe (eV) 
  real(fp), parameter :: e = 2.7182818285_fp
  real(fp), parameter :: me_Kg  = 9.10938215e-31_fp  ! electron mass (kg)
  real(fp), parameter :: mp_Kg  = 1.672621637e-27_fp ! proton mass (kg)
  real(fp), parameter :: mp_g   = 1.672621637e-24_fp ! proton mass (g)
  real(fp), parameter :: amu_g  = 1.660539066e-24_fp ! AMU_g 
  !-----------------------------------------
  ! Physics constants, based on CODATA 2006
  real(fp), parameter :: ZEL  = 4.803206799125e-10_fp !ELECTRON CHARGE (STATCOULOMBS)
  real(fp), parameter :: AEE  = 1.602176487e-19_fp ! elementary charge
  real(fp), parameter :: AME  = 9.10938215e-31_fp  ! electron mass (kg)
  real(fp), parameter :: AMP  = 1.672621637e-27_fp ! proton mass (kg)
  real(fp), parameter :: VC   = 2.99792458e8_fp    ! speed of light
  real(fp), parameter :: RMU0 = (4.0e-7_fp)*PI     ! permeability
  real(fp), parameter :: EPS0 = ONE/(VC*VC*RMU0)   ! permittivity
  real(fp), parameter :: usdp = rmu0
  real(fp), parameter :: aee_amp = 9.578833391e+7  ! electron_charge/proton_mass (C*kg^-1)
  real(fp), parameter :: evptocm_sec = 1.384112291e+6_fp !sqrt(2*kb_keV*10^7/mp_g) 
                         !cm/sec for 1eV proton, note J=10^7 erg
  !Define conversion factors
  real(fp), parameter :: zcmtom = 1.0e-2_fp  ! Centimeter to meter
  real(fp), parameter :: zmtocm = 1.0e+2_fp  ! meters to centimeters
  real(fp), parameter :: zcm2tom2 = 1.0e-4_fp  ! cm^2 to m^2
  real(fp), parameter :: zm2tocm2 = 1.0e+4_fp  ! m^2 to cm^2
  real(fp), parameter :: zcm3tom3 = 1.0e-6_fp  ! cm^3 to m^3 | 1/m^3 to 1/cm^3
  real(fp), parameter :: zm3tocm3 = 1.0e+6_fp  ! m^3 to cm^3 | 1/cm^3 to 1/m^3
  real(fp), parameter :: zwtomw = 1.0e-6_fp  ! Watts to Mega-watts
  real(fp), parameter :: zatoma = 1.0e-6_fp  ! Amperes to Mega-amperes
  real(fp), parameter :: zev2kev = 1.0e-3_fp ! eV to keV
  real(fp), parameter :: zkev2ev = 1.0e+3_fp ! keV to eV
  real(fp), parameter :: cgs2mks = 1.0E-7_fp ! CGS TO MKS
  real(fp), parameter :: t2gauss = 1.0E+4_fp ! Tesla to Gauss
  real(fp), parameter :: gauss2t = 1.0E-4_fp ! Gauss to Tesla

  real(fp), parameter :: rad2deg = 180.0_fp/pi
  real(fp), parameter :: deg2rad = pi/180.0_fp


  real(fp), parameter :: epslon = 1.0e-34_fp    ! small number
  real(fp), parameter :: epsinv = 1.0e+34_fp    ! large number

!contains

  ! potential area for conversion utils etc

end module physconst_mod
