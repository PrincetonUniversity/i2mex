subroutine i2mex_error(ier)
 
  ! Print error message
 
  implicit none
  integer, intent(in) :: ier
 
  select case(ier)
 
  case(1)
     print*,'--i2mex ERROR-- occurred while setting up S'
  case(2)
     print*,'--i2mex ERROR-- occurred while setting up T'
  case(3)
     print*,'--i2mex ERROR-- occurred while setting up PSI'
  case(4)
     print*,'--i2mex ERROR-- occurred while setting up P'
  case(5)
     print*,'--i2mex ERROR-- occurred while setting up G'
  case(6)
     print*,'--i2mex ERROR-- occurred while setting up Q'
  case(7)
     print*,'--i2mex ERROR-- occurred while setting up X'
  case(8)
     print*,'--i2mex ERROR-- occurred while setting up Z'
 
  case(9)
     print*,'--i2mex ERROR-- occurred while while attempting to free the buffer'
 
  case(10)
     print*,'--i2mex ERROR-- occurred computing S (i2mex_getS) from PSI: PSI(NS)=PSI(1)'
  case(11)
     print*,'--i2mex ERROR-- occurred computing DS (i2mex_getDS) from PSI: PSI(NS)=PSI(1)'
  case(12)
     print*,'--i2mex ERROR-- occurred computing D2S (i2mex_getD2S) from PSI: PSI(NS)=PSI(1)'
 
  case(13)
     print*,'--i2mex ERROR-- occurred while computing P (i2mex_getP)'
  case(14)
     print*,'--i2mex ERROR-- occurred while computing PP (i2mex_getPP)'
 
  case(16)
     print*,'--i2mex ERROR-- occurred while computing Q (i2mex_getQ)'
  case(17)
     print*,'--i2mex ERROR-- occurred while computing QP (i2mex_getQP)'
  case(18)
     print*,'--i2mex ERROR-- occurred while computing QPP (i2mex_getQPP)'
 
  case(19)
     print*,'--i2mex ERROR-- occurred while computing G (i2mex_getG)'
  case(20)
     print*,'--i2mex ERROR-- occurred while computing GP (i2mex_getGP)'
 
  case(21)
     print*,'--i2mex ERROR-- occurred while computing X (i2mex_getX)'
  case(22)
     print*,'--i2mex ERROR-- occurred while computing grad X (i2mex_getGradX)'
  case(23)
     print*,'--i2mex ERROR-- occurred while computing d^2 X/d t^2 (i2mex_getXtt)'
  case(24)
     print*,'--i2mex ERROR-- occurred while computing d^2 X/d t d p (i2mex_getXpt)'
  case(25)
     print*,'--i2mex ERROR-- occurred while computing d^2 X/d p^2 (i2mex_getXpp)'
 
  case(26)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getX'
  case(27)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getGradX'
  case(28)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getXtt'
  case(29)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getXpt'
  case(30)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getXpp'
 
 
  case(31)
     print*,'--i2mex ERROR-- occurred while computing Z (i2mex_getZ)'
  case(32)
     print*,'--i2mex ERROR-- occurred while computing grad Z (i2mex_getGradZ)'
  case(33)
     print*,'--i2mex ERROR-- occurred while computing d^2 Z/d t^2 (i2mex_getZtt)'
  case(34)
     print*,'--i2mex ERROR-- occurred while computing d^2 Z/d t d p (i2mex_getZpt)'
  case(35)
     print*,'--i2mex ERROR-- occurred while computing d^2 Z/d p^2 (i2mex_getZpp)'
 
  case(36)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getZ'
  case(37)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getGradZ'
  case(38)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getZtt'
  case(39)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getZpt'
  case(40)
     print*,'--i2mex ERROR-- lacking memory in i2mex_getZpp'
 
 
  case(41)
     print*,'--i2mex ERROR-- occurred while computing Jacobian (i2mex_getJ)'
  case(42)
     print*,'--i2mex ERROR-- occurred while computing d Jacobian/d t (i2mex_getJt)'
  case(43)
     print*,'--i2mex ERROR-- occurred while computing d Jacobian/d p (i2mex_getJp)'
  case(44)
     print*,'--i2mex ERROR-- lacked memory in i2mex_getJ'
  case(45)
     print*,'--i2mex ERROR-- lacked memory in i2mex_getJt'
  case(46)
     print*,'--i2mex ERROR-- lacked memory in i2mex_getJp'
 
  case(47)
     print*,'--i2mex ERROR--  occurred while setting up G in i2mex_scaleQ'
  case(48)
     print*,'--i2mex ERROR--  occurred while setting up Q in i2mex_scaleQ'
  case(49)
     print*,'--i2mex ERROR--  occurred while re-computing spline coefficients for PSI in i2mex_scaleQ'
 
  case(50)
     print*,'--i2mex ERROR--  occurred while saving object in i2mex_save'
  case(51)
     print*,'--i2mex ERROR--  occurred while restoring object in i2mex_load'
  case(52)
     print*,'--i2mex ERROR--  occurred while requesting id_s in i2mex_load'
  case(53)
     print*,'--i2mex ERROR--  occurred while requesting id_t in i2mex_load'
  case(54)
     print*,'--i2mex ERROR--  occurred while requesting id_psi in i2mex_load'
  case(55)
     print*,'--i2mex ERROR--  occurred while requesting id_p in i2mex_load'
  case(56)
     print*,'--i2mex ERROR--  occurred while requesting id_g in i2mex_load'
  case(57)
     print*,'--i2mex ERROR--  occurred while requesting id_q in i2mex_load'
  case(58)
     print*,'--i2mex ERROR--  occurred while requesting id_x in i2mex_load'
  case(59)
     print*,'--i2mex ERROR--  occurred while requesting id_z in i2mex_load'
 
  case(60)
     print*,'--i2mex ERROR--  occurred in i2mex_getDelStar'
  case(61)
     print*,'--i2mex ERROR--  occurred in i2mex_getXJparallel'
  case(62)
     print*,'--i2mex ERROR--  occurred in i2mex_getGsError'
 
  case(63)
     print*,'--i2mex ERROR--  1st integrand dimension too small in i2mex_getPoloidalIntegral'
  case(64)
     print*,'--i2mex ERROR--  2nd integrand dimension too small in i2mex_getPoloidalIntegral'
 
  case(65)
     print*,'--i2mex ERROR--  occurred in i2mex_getOriT'
  case(66)
     print*,'--i2mex ERROR--  occurred in i2mex_getOriPsi'
  case(67)
     print*,'--i2mex ERROR--  occurred in i2mex_getOriP'
  case(68)
     print*,'--i2mex ERROR--  occurred in i2mex_getOriG'
  case(69)
     print*,'--i2mex ERROR--  occurred in i2mex_getOriQ'
 
  case(70)
     print*,'--i2mex ERROR--  incompatible NT1 in i2mex_getOriX'
  case(71)
     print*,'--i2mex ERROR--  incompatible NS in i2mex_getOriX'
  case(72)
     print*,'--i2mex ERROR--  error in i2mex_getOriX'
  case(73)
     print*,'--i2mex ERROR--  incompatible NT1 in i2mex_getOriZ'
  case(74)
     print*,'--i2mex ERROR--  incompatible NS in i2mex_getOriZ'
  case(75)
     print*,'--i2mex ERROR--  error in i2mex_getOriZ'
 
  case(76)
     print*,'--i2mex ERROR--  error in i2mex_getSurface'
  case(77)
     print*,'--i2mex ERROR--  error in i2mex_getVolume'
  case(78)
     print*,'--i2mex ERROR--  error in i2mex_getVolumeAveragedPressure'
  case(79)
     print*,'--i2mex ERROR--  error in i2mex_getVolumeAveragedBSquare'
  case(80)
     print*,'--i2mex ERROR--  error in i2mex_getBeta'
  case(81)
     print*,'--i2mex ERROR--  error in i2mex_getBetaToroidal'
  case(82)
     print*,'--i2mex ERROR--  error in i2mex_getBetaPoloidalFreidberg'
  case(83)
     print*,'--i2mex ERROR--  error in i2mex_getBetaPoloidal'
  case(84)
     print*,'--i2mex ERROR--  error in i2mex_getBetaN'
  case(85)
     print*,'--i2mex ERROR--  error in i2mex_getPlasmaCurrent'
  case(86)
     print*,'--i2mex ERROR--  error in i2mex_getLi'
 
 
  case(87)
     print*,'--i2mex ERROR--  error in i2mex_'
  case(88)
     print*,'--i2mex ERROR--  error in i2mex_get'
  case(89)
     print*,'--i2mex ERROR--  error in i2mex_get'
 
  case(90)
     print*,'--i2mex ERROR--  could not find equilibrium input file '
  case(91)
     print*,'--i2mex ERROR--  an error occurred while calling cdfGetVar '
  case(92)
     print*,'--i2mex ERROR--  plasma representation appears to cover only half plane '
  case(93)
     print*,'--i2mex ERROR--  lacking memory in i2mex_fromInp1CDF'
  case(94)
     print*,'--i2mex ERROR--  failed to read data in i2mex_fromInp1CDF'
  case(95)
     print*,'--i2mex ERROR--  failed to call i2mex_init in i2mex_fromInp1CDF'
  case(96)
     print*,'--i2mex ERROR--  failed free memory in i2mex_fromInp1CDF'
 
  case(97)
     print*,'--i2mex ERROR--  could not connect to MDSPlus path in i2mex_fromMDSPlus'
   case(98)
     print*,'--i2mex ERROR--  a minor problem occurred when setting the time in i2mex_fromMDSPlus'
     print*,'Most probably,the specified time slice does not fit in the time interval'
 
   case(99)
     print*,'--i2mex ERROR--  a serious problem occurred when setting the time in i2mex_fromMDSPlus'
   case(100)
     print*,'--i2mex ERROR--  no XB profile found in i2mex_fromMDSPlus'
   case(101)
     print*,'--i2mex ERROR--  RPDIMS error occurred in i2mex_fromMDSPlus'
   case(102)
     print*,'--i2mex ERROR--  Failed to access MHD profiles in i2mex_fromMDSPlus'
   case(103)
     print*,'--i2mex ERROR--  MHD profiles not ready in i2mex_fromMDSPlus'
 
   case(104)
     print*,'--i2mex ERROR--  unsupported format in i2mex_ReadEquilibrium'
 
   case(105)
     print*,'--i2mex WARNING--  no rational surface was found '
 
   case(106)
     print*,'--i2mex ERROR-- occurred in i2mex_getGlasserGreenJohnsonEFH '
 
   case(107)
     print*,'--i2mex ERROR-- occurred in i2mex_getV '
   case(108)
     print*,'--i2mex ERROR-- occurred in i2mex_getVp '
 
   case(109)
     print*,'--i2mex ERROR-- occurred in i2mex_getJdotBOverBSquare '
 
   case(110)
     print*,'--i2mex ERROR-- occurred in i2mex_getPhi '
   case(111)
     print*,'--i2mex ERROR-- occurred in i2mex_getPhiP '
   case(112)
     print*,'--i2mex ERROR-- occurred in i2mex_getPhiPP '
 
   case(113)
     print*,'--i2mex ERROR-- while attempting to access the old geometry and profiles in i2mex_refineESCPAndQ'
   case(114)
     print*,'--i2mex ERROR-- while interpolating profiles onto ESC grid in i2mex_refineESCPAndQ'
   case(115)
     print*,'--i2mex WARNING-- while running ESC in p & j// mode in i2mex_refineESCPAndQ'
   case(116)
     print*,'--i2mex WARNING-- while running ESC in p & q mode in i2mex_refineESCPAndQ'
 
   case(117)
     print*,'--i2mex ERROR-- could not find GEQDSK file in i2mex_fromGeqdsk'
   case(118)
     print*,'--i2mex WARNING-- failed to compute an accurate equilibrium in i2mex_fromGeqdsk'
   case(119)
     print*,'--i2mex WARNING-- Psi not monotonic increasing in i2mex_fromGeqdsk.'
     print*,'                Most probably due to high triangularity and/or'
     print*,'                proximity to an X point. Radial mesh will be '
     print*,'                truncated at max(Psi)!'
 
  case(120)
     print*,'--i2mex ERROR-- attempt to integrate along a closed contour failed in i2mex_contour:'
     print*,'                LSODE error flag raised'
  case(121)
     print*,'--i2mex ERROR-- attempt to integrate along a closed contour failed in i2mex_contour: '
     print*,'                Max no of integration steps reached. '
     print*,'                May need to increase the nsteps in contour.f90 '
     print*,'                This could also be symptomatic of an open contour.'
  case(122)
     print*,'--i2mex ERROR-- attempt to integrate along a closed contour failed in i2mex_contour: '
     print*,'                Spline interpolation error. '
  case(123)
     print*,'--i2mex WARNING-- inaccurate closed contour computation in i2mex_contour: '
 
  case(124)
     print*,'--i2mex ERROR-- Need at least 5 radial points in i2mex_getQFromPsiGXZ '
  case(125)
     print*,'--i2mex ERROR-- Need at least 5 poloidal points in i2mex_getQFromPsiGXZ '
 
  case(126)
     print*,'--i2mex ERROR-- Could not open file in i2mex_fromMenard'
  case(127)
    print*,'--i2mex ERROR-- Memory allocation error in i2mex_fromMenard'
  case(128)
    print*,'--i2mex ERROR-- Failed to read 1-D profiles in i2mex_fromMenard'
  case(129)
    print*,'--i2mex ERROR-- Failed to read 2-D profiles in i2mex_fromMenard'
  case(130)
    print*,'--i2mex ERROR-- Failed to initialize i2mex in i2mex_fromMenard'
  case(131)
    print*,'--i2mex ERROR-- Unkown error in i2mex_fromMenard'
 
  case(132)
    print*,'--i2mex ERROR-- Interpolation error for X in i2mex_poloidalRemap'
  case(133)
    print*,'--i2mex ERROR-- Interpolation error for Z in i2mex_poloidalRemap'
  case(134)
    print*,'--i2mex ERROR-- Remapping error in i2mex_poloidalRemap'
  case(135)
    print*,'--i2mex ERROR-- Unknown error in i2mex_poloidalRemap'
 
  case(136)
    print*,'--i2mex ERROR-- Unknown error in i2mex_scaleG'
  case(137)
    print*,'--i2mex ERROR-- Failed to recompute spline coefficients for PSI in i2mex_scaleG'
  case(138)
    print*,'--i2mex ERROR-- Failed to recompute spline coefficients for P in i2mex_scaleG'
  case(139)
    print*,'--i2mex ERROR-- Failed to recompute spline coefficients for G in i2mex_scaleG'
 
  case(140)
     print*,'--i2mex ERROR-- Occurred in i2mex_getFrame while attempting to interpolate xbound'
  case(141)
     print*,'--i2mex ERROR-- Occurred in i2mex_getFrame while attempting to interpolate zbound'
 
  case(142)
     print*,'--i2mex ERROR-- Occurred in i2mex_getDelta'
  case(143)
     print*,'--i2mex ERROR-- Occurred in i2mex_getQDeltaP'
 
  case(144)
     print*,'--i2mex ERROR-- Occurred in i2mex_toAxis: Need >=5 points to extrpolate to axis'
 
  case(145)
     print*,'--i2mex ERROR--  error in i2mex_getUpsilon'
 
  case(146)
     print*,'--i2mex ERROR--  error in i2mex_getSpest1'
 
  case(147)
     print*,'--i2mex ERROR--  error in i2mex_i2mex_fromFreeqBe: failed to read file'
 
  case(148)
     print*,'--i2mex ERROR-- could not find FREEQBE file in i2mex_fromFreeqBe'
  case(149)
     print*,'--i2mex ERROR-- Some error occurred in i2mex_fromFreeqbe'
 
  case(150)
     print*,'--i2mex ERROR-- Some error occurred in i2mex_fromGeqdsk'

  case(151)
     print*,'--i2mex ERROR-- error occurred after calling xplasma::eqm_rzgrid in i2mex_fromGeqdsk'
  case(152)
     print*,'--i2mex ERROR-- error occurred after calling xplasma::eqm_dbdy in i2mex_fromGeqdsk'
  case(153)
     print*,'--i2mex ERROR-- error occurred after calling xplasma::eqm_rzfunda in i2mex_fromGeqdsk'
  case(154)
     print*,'--i2mex ERROR-- error occurred after calling xplasma::eqm_cbdy in i2mex_fromGeqdsk'

  case(155)
     print*,'--i2mex ERROR-- fraction too high in i2mex_toAxisFraction'

  case(156)
     print*,'--i2mex ERROR-- occurred in i2mex_decimate'
  case(157)
     print*,'--i2mex WARNING--  to few points in i2mex_toEdgeFraction to extrapolate'
  case(158)
     print*,'--i2mex ERROR--  fraction too high in i2mex_toEdgeFraction'
  case(159)
     print*,'--i2mex ERROR--  error occurred after calling i2mex_getPoloidalAverageOneOverRSquare'
  case(160)
     print*,'--i2mex ERROR--  error occurred after calling i2mex_getPoloidalAverageBSquare'
 
  case(161)
     print*,'--i2mex ERROR-- too few points in i2mex_toAxisFraction'

  case(162)
     print*,'--i2mex ERROR-- failed to compute |grad psi|^2 from direct representation in i2mex_getMetric'

  case(163)
     print*,'--i2mex ERROR-- failed to compute d(q delta)/d the in i2mex_getQDeltaT'

  case(164)
     print*,'--i2mex ERROR-- Cannot find psibig or any variation of it in i2mex_input'
  case(165)
     print*,'--i2mex ERROR-- Cannot find p or any variation of it in i2mex_input'
  case(166)
     print*,'--i2mex ERROR-- Cannot find q or any variation of it in i2mex_input'
  case(167)
     print*,'--i2mex ERROR-- Cannot find g or any variation of it in i2mex_input'
  case(168)
     print*,'--i2mex ERROR-- Cannot find r or any variation of it in i2mex_input'
  case(169)
     print*,'--i2mex ERROR-- Cannot find z or any variation of it in i2mex_input'
  case(170)
    print*,'--i2mex ERROR-- error in i2mex_fromCEO'

  case(171)
    print*,'--i2mex ERROR-- error in i2mex_get_iboozer'

  case(172)
    print*,'--i2mex ERROR-- error in i2mex_get_nuboozer'

  case(173)
    print*,'--i2mex ERROR-- error in i2mex_get_iboozer: theta(orig) not found.'

  case(174)
    print*,'--i2mex ERROR-- error in i2mex_get_iboozer: not monotonic in theta'

  case(175)
    print*,'--i2mex ERROR-- error in i2mex_get_nuboozer: not periodic in theta'

  case(176)
    print*,'--i2mex ERROR-- ESC is not supported in this version'

  end select
 
end subroutine i2mex_error
 
