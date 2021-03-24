program test

  write(*,'(//a//)') '**TEST 1**'
  call test1

  write(*,'(//a//)') '**TEST 2**'
  call test2

end program test

!-----------------------------------------------------------------------------

subroutine test1

  use i2mex_mod
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  integer nt1, ns
  real(r8), dimension(:), allocatable ::  t,  psi
  real(r8), dimension(:), allocatable ::  e, f, h, dr, di, v, vp
 

  integer i, j, ier, inputFormat
  character*130 inputFile
  real(r8) time

  integer, parameter :: nti1=129, nsi=101
  real(r8) :: psimin, psimax, psis, qs, psia, psib
  real(r8) ti(nti1), psii(nsi), qi(nsi), gi(nsi), presi(nsi), phii(nsi)
  real(r8), dimension(nti1, nsi) :: xi, zi
  real(r8), dimension(nti1, nsi) :: xti, zti, xpi, zpi
  real(r8), dimension(nti1, nsi) ::  xtti, ztti, xpti, zpti, xppi, zppi
  real(r8), dimension(nti1, nsi) ::  jac, jt, jp
  real(r8), dimension(nti1, nsi) ::  gserror
  integer m, n, tic, toc, tics_per_second, mth

  real(r8) surface, res, rmajor, rma
  
  write(*,*)' Enter the input equilibrium file format'
  write(*,*)' -1  for CHEASE INP1 format.'
  write(*,*)'  0  for TRANSP data in UFILE format (requires MDSPlus connection)'
  write(*,*)' +1  for CHEASE inp1.cdf format'
  write(*,*)' +2  for JSOLVER eqdsk.cdf format'
  write(*,*)' +3  for EFIT G-EQDSK format'
  write(*,*)' +4  for EFIT G-EQDSK format, rerunning equilibrium through ESC-Q'
  write(*,*)' +5  for Menard''s psipqgRZ netCDF format'
  write(*,*)' +6  for Belova''s freeqbe netCDF format'
  read(5,*) inputFormat
  if (inputFormat<-1 .or. inputFormat>6) inputFormat=1
  if (inputFormat==0) then
     write(*,*) ' Enter the path or MDSPlus node (requires MDSPlus connection)'
     write(*,*) ' Example:'
     write(*,*) 'MDS+:TRANSPGRID.PPPL.GOV:TRANSP_NSTX(NSTX.05,116313P03)'
     read(5,'(a80)') inputFile
     write(*,*) ' Enter the time slize in sec (data will be interpolated)'
     read(5,*) time
     call i2mex_fromMDSPlus(inputFile, &
          & time, nti1, nsi, i2mex_counterclockwise, ier)
  else
     write(*,*) ' Enter input file name '
     if(inputFormat == 1) then
        write(*,*) ' Example: inp1_d.cdf'
     endif
     if(inputFormat == 3 .or. inputFormat==4) then
        write(*,*) ' This can be an MDSPlus path as in '
        write(*,*) ' "MDS+/REDUCE:BIRCH.PPPL.GOV:8501:EFIT01(103984;t=0.23)"'
        write(*,*) ' for instance'
     endif
     if(inputFormat==-1) then
        write(*,*) ' Depending on your platform, you have the choice between:'
        write(*,*) ' INP1_ieee-le for binary in little endian format (Alphas, ix86...), or'
        write(*,*) ' INP1_ieee-be for binary in big endian format (Sun, ...). '
     endif
     read(5,*) inputFile
     call i2mex_ReadEquilibrium(inputFormat, inputFile, i2mex_counterclockwise, ier)

  endif

  call i2mex_error(ier)
  if(ier/=0) then
     write(*,*)' An error occurred while attempting to load in equilibrium data.'
     write(*,*)' Please check that the filename or MDSPlus path are correct.'
     stop
  endif


!!$  call i2mex_refineESCPAndQ(4, ier)
!!$  call i2mex_error(ier)

!!$  call i2mex_save('drive.i2mex', ier)
!!$  call i2mex_error(ier)
!!$
!!$  call i2mex_free(ier)
!!$
!!$  call i2mex_load('drive.i2mex', ier)
!!$  call i2mex_error(ier)

  call i2mex_getOriNt1(nt1, ier)
  call i2mex_error(ier)
  call i2mex_getOriNs(ns, ier)
  call i2mex_error(ier)

  allocate(t(nt1), psi(ns))
  allocate(e(ns), f(ns), h(ns), dr(ns), di(ns), v(ns))

  call i2mex_getOriPsi(ns, psi, ier)
  call i2mex_error(ier)
  call i2mex_getOriT(nt1, t, ier)
  call i2mex_error(ier)

  call i2mex_getV(ns, psi, v, ier)
  call i2mex_error(ier)

  call i2mex_getGlasserGreenJohnsonEFH(ns, psi, e, f, h, ier)
  call i2mex_error(ier)
  dr = e + f + h**2
  di = e + f + h - 0.25_i2mex_r8
  write(*,*) &
       & '      PSI            E             F             H            DR            DI'
  do i = 1, ns
     write(*,'(6e14.6)') psi(i), e(i), f(i), h(i), dr(i), di(i)
  enddo

  ti = i2mex_twopi_r8* (/ (real(i-1, r8)/real(nti1-1, r8), i=1, nti1) /)
  psimin = minval(psi)
  psimax = maxval(psi)
  psii = psimin + (psimax-psimin)* &
       & (/ (real(i-1, r8)/real(nsi-1, r8), i=1, nsi) /)

  write(*,'(/a)')' B-field rescaling g-edge=R0'
  call i2mex_getRmajor(rmajor, ier)
  call i2mex_error(ier)

  call i2mex_getRmagnetic(rma, ier)
  call i2mex_error(ier)
  
  write(*,'(a,f17.6,a,f17.6)') ' Major radius R0=',rmajor,' Magnetic axis Rm=', rma
  

  ! rescale B
  call i2mex_scaleG(psimax, rmajor, ier)
  call i2mex_error(ier)

  ! now we need to recompute psi (psi scales with B)
  call i2mex_getOriPsi(ns, psi, ier)
  call i2mex_error(ier)
  psimin = minval(psi)
  psimax = maxval(psi)
  psii = psimin + (psimax-psimin)* &
       & (/ (real(i-1, r8)/real(nsi-1, r8), i=1, nsi) /)


  
  psia = 0.0_i2mex_r8
  psib = psimax
  call i2mex_findRationalSurface(3, 2, psia, psib, psis, ier)
  call i2mex_error(ier)
  write(*,'(a,f17.10)')' Rational surface (3,2) @ psi=',psis
!!$  call i2mex_getQ(1, psis, qs, ier)
!!$  call i2mex_error(ier)
!!$  write(*,'(a,f17.10)')' q(psis)                     =', qs
  

  call i2mex_getX(nti1, nsi, ti, psii, xi, ier)
  call i2mex_error(ier)
  call i2mex_getGradX(nti1, nsi, ti, psii, xti, xpi, ier)
  call i2mex_error(ier)
  call i2mex_getXtt(nti1, nsi, ti, psii, xtti, ier)
  call i2mex_error(ier)
  call i2mex_getXpt(nti1, nsi, ti, psii, xpti, ier)
  call i2mex_error(ier)
  call i2mex_getXpp(nti1, nsi, ti, psii, xppi, ier)
  call i2mex_error(ier)

  call i2mex_getZ(nti1, nsi, ti, psii, zi, ier)
  call i2mex_error(ier)
  call i2mex_getGradZ(nti1, nsi, ti, psii, zti, zpi, ier)
  call i2mex_error(ier)
  call i2mex_getZtt(nti1, nsi, ti, psii, ztti, ier)
  call i2mex_error(ier)
  call i2mex_getZpt(nti1, nsi, ti, psii, zpti, ier)
  call i2mex_error(ier)
  call i2mex_getZpp(nti1, nsi, ti, psii, zppi, ier)
  call i2mex_error(ier)

  call i2mex_getJ(nti1, nsi, ti, psii, jac, ier)
  call i2mex_error(ier)
  call i2mex_getJt(nti1, nsi, ti, psii, jt, ier)
  call i2mex_error(ier)
  call i2mex_getJp(nti1, nsi, ti, psii, jp, ier)
  call i2mex_error(ier)

  i = nti1/2
  j = nsi/3

  write(*,'(/a,f15.10,a,f15.10)')' xmax=', maxval(xi)
  write(*,'(a,f15.10,a,f15.10/)')' xmin=', minval(xi)

  write(*,'(a,f15.10,a,f15.10)') &
       & ' xti =',xti(i,j),' approx->', (xi(i+1,j)-xi(i-1,j))/(ti(i+1)-ti(i-1))
  write(*,'(a,f15.10,a,f15.10)')&
       & ' xpi =',xpi(i,j),' approx->', (xi(i,j+1)-xi(i,j-1))/(psii(j+1)-psii(j-1))
  write(*,'(a,f15.10,a,f15.10)')' xtti=',xtti(i,j), &
       & ' approx->', 4._r8*(xi(i+1,j)-2._r8*xi(i,j)+xi(i-1,j))/(ti(i+1)-ti(i-1))**2
  write(*,'(a,f15.10,a,f15.10,2f15.10)')' xpti=',xpti(i,j), &
       & ' approx->', (xi(i+1,j+1)-xi(i-1,j+1)-xi(i+1,j-1)+xi(i-1,j-1))/ &
       & ((ti(i+1)-ti(i-1))*(psii(j+1)-psii(j-1))), &
       & (xpi(i+1,j)-xpi(i-1,j))/(ti(i+1)-ti(i-1)), &
       & (xti(i,j+1)-xti(i,j-1))/(psii(j+1)-psii(j-1))
  write(*,'(a,f15.10,a,f15.10)') ' xppi=',xppi(i,j), &
       & ' approx->', 4._r8*(xi(i,j+1)-2._r8*xi(i,j)+xi(i,j-1))/(psii(j+1)-psii(j-1))**2
  
  write(*,'(/a,f15.10,a,f15.10)')' zmax=', maxval(zi)
  write(*,'(a,f15.10,a,f15.10/)')' zmin=', minval(zi)

  write(*,'(a,f15.10,a,f15.10)') &
       & ' zti =',zti(i,j),' approx->', (zi(i+1,j)-zi(i-1,j))/(ti(i+1)-ti(i-1))
  write(*,'(a,f15.10,a,f15.10)')&
       & ' zpi =',zpi(i,j),' approx->', (zi(i,j+1)-zi(i,j-1))/(psii(j+1)-psii(j-1))
  write(*,'(a,f15.10,a,f15.10)')' ztti=',ztti(i,j), &
       & ' approx->', 4._r8*(zi(i+1,j)-2._r8*zi(i,j)+zi(i-1,j))/(ti(i+1)-ti(i-1))**2
  write(*,'(a,f15.10,a,f15.10,2f15.10)')' zpti=',zpti(i,j), &
       & ' approx->', (zi(i+1,j+1)-zi(i-1,j+1)-zi(i+1,j-1)+zi(i-1,j-1))/ &
       & ((ti(i+1)-ti(i-1))*(psii(j+1)-psii(j-1))), &
       & (zpi(i+1,j)-zpi(i-1,j))/(ti(i+1)-ti(i-1)), &
       & (zti(i,j+1)-zti(i,j-1))/(psii(j+1)-psii(j-1))
  write(*,'(a,f15.10,a,f15.10)') ' zppi=',zppi(i,j), &
       & ' approx->', 4._r8*(zi(i,j+1)-2._r8*zi(i,j)+zi(i,j-1))/(psii(j+1)-psii(j-1))**2

  write(*,'(/a,f15.10)')' Jacobian=', jac(i,j)
  write(*,'(a,f15.10,a,f15.10)') ' jti =', jt(i,j), &
       & ' approx->', (jac(i+1,j)-jac(i-1,j))/(ti(i+1)-ti(i-1))
  write(*,'(a,f15.10,a,f15.10)') ' jpi =', jp(i,j), &
       & ' approx->', (jac(i,j+1)-jac(i,j-1))/(psii(j+1)-psii(j-1))

  call i2mex_getGsError(nti1, nsi, ti, psii, gserror, ier)
  call i2mex_error(ier)
  write(*,'(/a,e15.7)')' Average relative Grad-Shafranov error: ', &
       & sum(sum(gserror, dim=1))/real(nti1*nsi, r8)

  call i2mex_getG(nsi, psii, gi, ier)
  call i2mex_error(ier)
  call i2mex_getQ(nsi, psii, qi, ier)
  call i2mex_error(ier)
  write(*,'(/a,f15.10,a,f15.10)')' g-axis=',gi(1),' g-edge=',gi(nsi)
  write(*,'(a,f15.10,a,f15.10)')' q-axis=',qi(1),' q-edge=',qi(nsi)


  ! just checking that we still have the expected g vlaues
  call i2mex_getG(nsi, psii, gi, ier)
  call i2mex_error(ier)
  write(*,'(/a,f15.10,a,f15.10)')' g-axis=',gi(1),' g-edge=',gi(nsi)  
  
!!$  write(*,'(/a)')' Bateman rescaling of q profile'
!!$
!!$  call i2mex_scaleQ(psimin, 1.1_r8, ier)
!!$  call i2mex_error(ier)
!!$
!!$  call i2mex_getG(nsi, psii, gi, ier)
!!$  call i2mex_error(ier)
!!$  call i2mex_getQ(nsi, psii, qi, ier)
!!$  call i2mex_error(ier)
!!$  write(*,'(/a,f15.10,a,f15.10)')' g-axis=',gi(1),' g-edge=',gi(nsi)
!!$  write(*,'(a,f15.10,a,f15.10)')' q-axis=',qi(1),' q-edge=',qi(nsi)
!!$
!!$  call i2mex_getGsError(nti1, nsi, ti, psii, gserror, ier)
!!$  call i2mex_error(ier)
!!$  write(*,'(/a,e15.7)')' Average relative Grad-Shafranov error: ', &
!!$       & sum(sum(gserror, dim=1))/real(nti1*nsi, r8)

  call i2mex_getPhi(nsi, psii, phii, ier)
  call i2mex_error(ier)
  call i2mex_getP(nsi, psii, presi, ier)
  call i2mex_error(ier)
  call i2mex_getQ(nsi, psii, qi, ier)
  call i2mex_error(ier)
  call i2mex_getG(nsi, psii, gi, ier)
  call i2mex_error(ier)
  write(*,*)'      psi            phi             p              q              g'
  do i=1, nsi
     write(*,'(5e15.7)') psii(i), phii(i), presi(i), qi(i), gi(i)
  enddo

  ! global parameters
  
!!$  call system_clock(tic, tics_per_second)
  call i2mex_getSurface(res, ier)
!!$  call system_clock(toc, tics_per_second)
!!$  print*,'time to compute surface ', real(toc-tic)/real(tics_per_second),' sec'
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' Surface               =', res
  call i2mex_getVolume(res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' Volume                =', res
  call i2mex_getVolumeAveragedPressure(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' VolumeAveragedPressure=',  res
  call i2mex_getVolumeAveragedBSquare(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' VolumeAveragedBSquare =',  res
  call i2mex_getBeta(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' Beta                  =',  res
  call i2mex_getBetaStar(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' BetaStar              =',  res
  call i2mex_getRminor(res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' Rminor                =',  res
  call i2mex_getRmajor(res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' Rmajor                =',  res
  call i2mex_getB0(res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' B0                    =',  res
  call i2mex_getPlasmaCurrent(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' PlasmaCurrent         =',  res
  call i2mex_getBetaToroidal(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' BetaToroidal          =',  res
  call i2mex_getBetaPoloidalFreidberg(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' BetaPoloidalFreidberg =',  res
  call i2mex_getBetaPoloidal(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' BetaPoloidal          =',  res
  call i2mex_getBetaN(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' BetaN                 =',  res
  call i2mex_getLi(nti1, nsi, ti, psii, res, ier)
  call i2mex_error(ier)
  write(*,'(a,f15.10)')' Li                    =',  res

  call i2mex_ToMapdsk(nti1, nsi, ti, psii, 'mapdsk.cdf', ier)

  ! plot results

  call i2mex_2plotmtv(nti1, nsi, ti, psii, 'drive.mtv', ier)
  call i2mex_error(ier)

  ! plot error

  call i2mex_2Dx(nti1, nsi, ti, psii, gserror, 'gserror.dx', ier)
  call i2mex_error(ier)

  ! clean up
  call i2mex_free(ier)

  deallocate(t, psi)
  deallocate(e, f, h, dr, di, v)


  if(ier==0) then
     print *,'Done! Test 1 was successful'
     print *,'The results can now be visualized using a variety of tools.'
     print *,'If you have "plotmtv" in your path then type "plotmtv drive.mtv"'
     print *,'If you have matlab installed with the mexCDF toolbox, type '
     print *,'   mapdsk at the matlab prompt, this will display a variety'
     print *,'   of profile and metric quantities.'
     print *,'If you have dx installed, run the gserror.net visual program'
     print *,'   to inspect the normalized Grad-Shafranov error.'
  else
     print *,'Error ier=', ier,' occurred during test run'
  endif


end subroutine test1
  

!----------------------------------------------------------------------------

subroutine test2

  ! Typical interface to M3D

  use i2mex_mod
  implicit none

  integer, parameter :: r8 = selected_real_kind(12,100)
  integer nt1, ns
  real(r8), dimension(:), allocatable ::  t, psi, q, g, p
  real(r8), dimension(:,:), allocatable :: x, z, j_phi, b_phi
  integer i, j, ier

  ! Read from CHEASE file inp1.cdf
  ! (set NIDEAL=3 in CHEASE namelist)

  call i2mex_fromInp1CDF('inp1sym.cdf', i2mex_clockwise, ier)
  call i2mex_error(ier)

  ! Get the original grid sizes to
  ! allocate space for the data

  call i2mex_getOriNt1(nt1, ier)
  call i2mex_error(ier)
  call i2mex_getOriNs(ns, ier)
  call i2mex_error(ier)

  write(*,'(a,i5,a,i5)') ' Grid sizes from CHEASE nt1*ns=',nt1,'*',ns

  ! allocate memory for the data

  allocate(t(nt1), psi(ns), q(ns), p(ns), g(ns))
  allocate(x(nt1, ns), z(nt1, ns), j_phi(nt1, ns), b_phi(nt1, ns))

  ! get the original grids

  call i2mex_getOriPsi(ns, psi, ier)
  call i2mex_error(ier)
  call i2mex_getOriT(nt1, t, ier)
  call i2mex_error(ier)

  ! Get and compute the data

  call i2mex_getQ(ns, psi, q, ier)
  call i2mex_error(ier)

  call i2mex_getG(ns, psi, g, ier)
  call i2mex_error(ier)

  call i2mex_getP(ns, psi, p, ier)
  call i2mex_error(ier)

  call i2mex_getX(nt1, ns, t, psi, x, ier)
  call i2mex_error(ier)

  call i2mex_getZ(nt1, ns, t, psi, z, ier)
  call i2mex_error(ier)

  call i2mex_getXJphi(nt1, ns, t, psi, j_phi, ier)
  call i2mex_error(ier)
  j_phi = j_phi/x

  b_phi = spread(g, dim=1, ncopies=nt1)/x


  ! Clean-up

  call i2mex_free(ier)

  if(ier==0) then
     print *,'Test 2 was successful'
  else
     print *,'Error ier=', ier,' occurred during test run'
  endif

  deallocate(t, psi, q, p, g)
  deallocate(x, z, j_phi, b_phi)

end subroutine test2
  
