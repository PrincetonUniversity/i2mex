program mex2eqs
 
  use ezcdf
  use i2mex_mod
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
 
 
  integer :: nt1=129, ns=101
 
  real(r8) :: psi_min, psi_max, xmax, xmin, zmax, zmin, z0, rerror
  real(r8) :: eps, area, vol, elong, beta, betap_, betat, betan, ima, li
  real(r8) :: r0, b0old, b0new
  real(r8), dimension(:), allocatable ::  psii
  real(r8), dimension(:), allocatable ::  t,  psi, p, pp, g, gp, q, phi, jpar
  real(r8), dimension(:), allocatable :: bsquare, oneOverRsquare
  real(r8), dimension(:,:), allocatable ::  xa, za, xj, i2d, nu2d
  real(r8), dimension(:,:), allocatable ::  gtt, gtp, gpp
  real(r8), dimension(:,:), allocatable ::  gserror
  real(r8), dimension(:,:), allocatable ::  temp
 
  logical :: idebug
  integer :: is
  real(r8), dimension(:,:), allocatable :: ipgq
  real(r8), dimension(:,:), allocatable :: jb2
  real(r8), dimension(:), allocatable :: wk


  integer i, j, ier, inputFormat, nsi
  character*130 inputFile
  character*12 filename
  real(r8) time
  integer m, n, k, tic, toc, tics_per_second, dims(3), iu
 
  character*8 today
  character*10 now
  character*80 lon_date
 
  real(r8) :: xma_, zma_
 
 
  integer :: nths, nsf
 
  ! uread stuff
  integer lnout, idcod
  real(r8) r8dcod
  logical lxout,lyes,lchar,ireduce,ilmds
  character*150 zfile

  ! x plasma stuff
  integer :: NR, NZ, iauto, id_R, id_Z
  real(r8) :: ztol, dlim, deltaW, deltaS, deltaE, deltaN, edge_smoothing
  real(r8), dimension(:), allocatable :: Rgrid, Zgrid, Psigrid
  real(r8), allocatable :: Rcoord(:,:), Zcoord(:,:), psiRZ(:,:),BRBZBPhi(:,:,:)
  integer, allocatable :: nregion(:)
  external eqm_brz_adhoc, eqm_bpsi, eqm_ball

  call uinita(0)
  nths = 133
  nsf = 194

  ! to support Boozer coords, tolerate much more fluctuation in det(J) than
  ! is the default for i2mex/xplasma...

  call eqm_select('mex2eqs (preliminary)',1)
  call eqi_jacheck_maxvar_Set(0.0000001_r8)
 
  if(lxout(0)) then
     write(lnout(0),*) ' '
     write(lnout(0),*) ' MEX2EQS:  retrieve and map equilibrium data for '
     write(lnout(0),*) ' postprocessing by EQS, the front end to ORBIT'
     write(lnout(0),*) ' '
     write(lnout(0),*) ' Important: EQS is statically linked so in order to '
     write(lnout(0),*) ' proceed you need to know the NTHS and NSF parameters'
     write(lnout(0),*) ' that have been set in the EQS commonblock o.cln.'
     write(lnout(0),*) ' NTHS = max no of poloidal rays + 5 (133)'
     write(lnout(0),*) ' NSF  = max no of radial nodes + 1 (194)'
     write(lnout(0),*) ' '
 
  endif
      
  CALL UREAD('Enter NTHS$')
  NTHS = IDCOD(IER)
  CALL UREAD('Enter NSF$')
  NSF  = IDCOD(IER)
 
  if (nths<21) nths = 133
  if (nsf<21) nsf = 194
 
  write(lnout(0),*) ' NTHS * NSF = ', NTHS,' * ', NSF
 
  ! set the grid sizes
  ns = NSF -1
  nt1 = NTHS - 4
 
  write(lnout(0),*) ' No of poloidal sections = ', nt1-1
  write(lnout(0),*) ' No of radial nodes      = ', ns
 
 
  ! set jacobian for Boozer coordinates
  i2mex_remap=1
  i2mex_mx=1; i2mex_npsi=1; i2mex_kb=0 ! Equal-Arc
  i2mex_direct = 0  ! =1 to allow direct representation to co-exist

  write(lnout(0),*)' Choose your theta variable by specifying the Jacobian dependence'
  write(lnout(0),*)' J ~ X**m/(|grad psi|**n |B|**k). Popular choices are:'
  write(lnout(0),*)' PEST     : m=2, n=0, k=0'
  write(lnout(0),*)' Boozer   : m=0, n=0, k=2'
  write(lnout(0),*)' Equal-arc: m=1, n=1, k=0'
  call ureadi('Enter m$', i2mex_mx)
  call ureadi('Enter n$', i2mex_npsi)
  call ureadi('Enter k$', i2mex_kb)

  idebug = .FALSE.

  do
     write(lnout(0),*) ' ** enter "C" to clear debug plot flag;'
     write(lnout(0),*) ' ** enter "D" to set debug plot flag: ',idebug
     write(lnout(0),*) '    (for plot of I+g*q vs. det(J)*B**2 at end).'
     write(lnout(0),*) ' '
     write(lnout(0),*)' The equilibrium can be read from various sources and formats'
     write(lnout(0),*)' -1  for CHEASE INP1 format.'
     write(lnout(0),*)'  0  for TRANSP data in UFILE format (requires MDSPlus connection)'
     write(lnout(0),*)' +1  for CHEASE inp1.cdf format'
     write(lnout(0),*)' +2  for JSOLVER eqdsk.cdf format'
     write(lnout(0),*)' +3  for EFIT G-EQDSK format'
     write(lnout(0),*)' +4  for EFIT G-EQDSK format, using G-EQDSK vacuum psi'
     write(lnout(0),*)' +5  for Menard''s psipqgRZ netCDF format'
     write(lnout(0),*)' +6  for Belova''s freeqbe netCDF format'
 
     CALL UREAD('Enter number$')
     if(lchar('C')) then
        idebug=.FALSE.
        write(lnout(0),*) ' debug flag CLEARED.'
        cycle
     else if(lchar('D')) then
        idebug=.TRUE.
        write(lnout(0),*) ' debug flag SET.'
        cycle
     else
        exit
     endif
  enddo
  inputFormat = IDCOD(IER)
 
  if (inputFormat<-1 .or. inputFormat>5) inputFormat=1
  if (inputFormat==0) then
     write(lnout(0),*) ' You have chosen to read in TRANSP results (UFILES), these'
     write(lnout(0),*) ' can be stored locally in which case you have to provide'
     write(lnout(0),*) ' a <RUNID> string (e.g. "11115P07", "<my_dir>/11115P07)", or'
     write(lnout(0),*) ' in the MDSPlus tree in which case you need to pass a string'
     write(lnout(0),*) ' of the format MDS+:TRANSPGRID.PPPL.GOV:TRANSP(NSTX.00,11115P07)'
     write(lnout(0),*) ' for example.'
     CALL UREADLL(' Enter "filename" or MDSPlus path$', inputFile)
     write(lnout(0),*) 'Input file entered is ', inputFile
     CALL UREAD('Enter the time slice in sec$')
     time = R8DCOD(IER)
 
     call i2mex_fromMDSPlus(inputFile, &
          & time, nt1, ns, i2mex_counterclockwise, ier)
 
  else
     if(inputFormat==-1) then
        write(lnout(0),*) ' Depending on your platform, you have the choice between:'
        write(lnout(0),*) ' INP1_ieee-le for binary in little endian format (Alphas, ix86...), or'
        write(lnout(0),*) ' INP1_ieee-be for binary in big endian format (Sun, ...). '
     endif
     if(inputFormat==3 .or. inputFormat==4) then
        write(lnout(0),*) ' You have chosen to read in a GEQDSK file, which can '
        write(lnout(0),*) ' either can be stored locally (<path>/<filename>) or '
        write(lnout(0),*) ' remotely in the MDSPlus tree. An example of MDSPlus '
        write(lnout(0),*) ' path is:'
        write(lnout(0),*) ' MDS+/REDUCE:SKYLARK.PPPL.GOV:8501:EFIT06(116313;t=0.9)'
     endif
     IF (InputFormat==4) i2mex_direct = 1  ! direct representation psi(R,z)
     CALL UREADLL('Enter file name or MDSPlus path$',inputFile)
 
     write(lnout(0),*) ' InputFormat = ', InputFormat, ' inputFile: ', inputFile
     call i2mex_ReadEquilibrium(inputFormat, inputFile, i2mex_clockwise, ier)

  endif
 
  call i2mex_error(ier)
  if(ier /=0) then
     write(lnout(0),*) ' An error occurred while attempting to load in the'
     write(lnout(0),*) ' equilibrium data. Please check that the file name'
     write(lnout(0),*) ' or MDSPlus path is correct.'
     stop
  endif
 
  call i2mex_getOriNs(nsi, ier)
 
  ALLOCATE(t(nt1))
  ALLOCATE(psii(nsi))
  ALLOCATE(psi(ns), p(ns), pp(ns), q(ns), g(ns), gp(ns), phi(ns), jpar(ns))
  ALLOCATE(xa(nt1, ns), za(nt1, ns), xj(nt1, ns))
  ALLOCATE(nu2d(nt1,ns), i2d(nt1,ns))  ! Boozer terms
  ALLOCATE(gtt(nt1, ns), gtp(nt1, ns), gpp(nt1, ns))
  ALLOCATE(gserror(nt1, ns))
  ALLOCATE(temp(nths, nsf))
  ALLOCATE(bsquare(ns), oneOverRSquare(ns))
 
  t = i2mex_twopi_r8* (/ (real(i-1, r8)/real(nt1-1, r8), i=1, nt1) /)
  call i2mex_getOriPsi(nsi, psii, ier)
  psi_min = psii(1); psi_max = psii(nsi)
  psi = psi_min + (psi_max-psi_min)* &
       & (/ (real(i-1, r8)/real(ns-1, r8), i=1, ns) /)
 
 
  call i2mex_getP(ns, psi, p, ier)
  call i2mex_error(ier)
  call i2mex_getPp(ns, psi, pp, ier)
  call i2mex_error(ier)
  call i2mex_getQ(ns, psi, q, ier)
  call i2mex_error(ier)
  call i2mex_integrate1d(ns, psi, q, ns, psi, phi, ier)
  call i2mex_error(ier)
  call i2mex_getG(ns, psi, g, ier)
  call i2mex_error(ier)
  call i2mex_getGp(ns, psi, gp, ier)
  call i2mex_error(ier)
  call i2mex_getX(nt1, ns, t, psi, xa, ier)
  call i2mex_error(ier)
  call i2mex_getZ(nt1, ns, t, psi, za, ier)
  call i2mex_error(ier)
  call i2mex_getJ(nt1, ns, t, psi, xj, ier)

  call i2mex_getPolAvrgBSquare(nt1, ns, t, psi, bsquare, ier)
  call i2mex_error(ier)
  call i2mex_getPolAvrgOneOverRSquare(nt1, ns, t, psi, oneOverRSquare, ier)
  call i2mex_error(ier)
  jpar = - (gp*bsquare/g + pp)/oneOverRSquare
 
  call i2mex_getMetric(nt1, ns, t, psi, gtt, gtp, gpp, ier)
  call i2mex_error(ier)
 
  ! check GS error
  call i2mex_getGsError(nt1, ns, t, psi, gserror, ier)
  call i2mex_error(ier)
  call i2mex_2Dx(nt1, ns, t, psi, gserror, 'gserror.dx', ier)
  call i2mex_error(ier)
  rerror = sum(sum(gserror, dim=1))/(real(nt1,r8)*real(ns,r8))
  write(lnout(0),*)' '
  write(lnout(0),'(a,e10.2)')' Average Grad-Shafranov error: ', rerror
  if(abs(rerror)>0.01_r8) write(lnout(0),*) ' --Warning--: large rel GS error > 1%!'
  if(abs(rerror)>1._r8  ) then
     write(lnout(0),*) ' --WARNING--: LARGE rel GS error > 100%!'
     write(lnout(0),*) '   (GS error measure lacks rotation term)'
     write(lnout(0),*) ' Continue at your own risk.'
  endif
 
 
  xmax = maxval(maxval(xa,dim=1))
  xmin = minval(minval(xa,dim=1))
  zmax = maxval(maxval(za,dim=1))
  zmin = minval(minval(za,dim=1))
  eps = (xmax-xmin)/(xmax+xmin)
  r0 = (xmax+xmin)/2.0_r8
  elong = (zmax-zmin)/(xmax-xmin)
  call i2mex_getSurface(area, ier)
  call i2mex_getVolume(vol, ier)
  call i2mex_getBeta(nt1, ns, t, psi, beta, ier)
  call i2mex_getBetaPoloidal(nt1, ns, t, psi, betap_, ier)
  call i2mex_getBetaToroidal(nt1, ns, t, psi, betat, ier)
  call i2mex_getBetaN(nt1, ns, t, psi, betan, ier)
  call i2mex_getPlasmaCurrent(nt1, ns, t, psi, ima, ier)
  ima = ima/(0.4_r8*i2mex_pi_r8)
  call i2mex_getLi(nt1, ns, t, psi, li, ier)
  write(lnout(0),*)' '
  write(lnout(0),*)'  a/R0    Area   Volume   Elong   Beta'
  write(lnout(0),'(5f8.3)') eps, area, vol, elong, beta
  write(lnout(0),*)' '
  write(lnout(0),*)' Beta-p  Beta-t  Beta-N    I-MA     li'
  write(lnout(0),'(5f8.3)') betap_, betat, betan, ima, li
  write(lnout(0),*)' '
 
  write(lnout(0),*)' '

  write(lnout(0),*)' Dimensions of X, Z grid.'
  CALL UREAD('Enter NX$')
  NR = IDCOD(IER)
  CALL UREAD('Enter NZ$')
  NZ = IDCOD(IER)
  
 iauto = 0 ! 0 = user supplied
 ztol  = 1.0e-6*xmax
 dlim  = 100._r8 * Xmax

 ! Set up (R,Z) grid in order to compute 2-d functions of (R, Z)
 allocate(Rgrid(NR), Zgrid(NZ))

 write(lnout(0),*)' The size of the computational box encompasses the plasma boundary'
 write(lnout(0),*)' plus some buffer region. Choose the buffer size (dW, dE, dS, dN) '
 write(lnout(0),*)' such that '
 write(lnout(0),*)' XW = (1-dW)*xmin (west boundary of the computational box)'
 write(lnout(0),*)' XE = (1+dE)*xmax (east)'
 write(lnout(0),*)' ZS = (1+dS)*zmin (south)'
 write(lnout(0),*)' ZN = (1+dN)*zmax (north)'
 write(lnout(0),*)' where xmin/xmax/zmin/zmax are min/max values of the plasma boundary.'

 CALL UREAD('Enter dW$')
 deltaW = min(0.99_r8, max(0._r8, R8DCOD(IER)))
 CALL UREAD('Enter dE$')
 deltaE = max(0.01_r8, R8DCOD(IER))
 CALL UREAD('Enter dS$')
 deltaS = max(0.01_r8, R8DCOD(IER))
 CALL UREAD('Enter dN$')
 deltaN = max(0.01_r8, R8DCOD(IER))
 
 ! Build direct representation 

 if((InputFormat==3 .OR. InputFormat==4) .AND. i2mex_direct>0) then

    ! -> Use psi(r,z) data from EFIT file. In this case the psi(r,z)
    ! data should have already been built in i2mex input routine, 
    ! provided i2mex_direct has been set to 1.

    
    ! take original grid
!!$    Rgrid = i2mex_o%Rleft + &
!!$         & (i2mex_o%Rrigh-i2mex_o%Rleft)* &
!!$         & (/ (real(i-1,r8)/real(NR-1,r8), i=1, NR) /)
!!$    Zgrid = i2mex_o%Zbot  + &
!!$         & (i2mex_o% Ztop-i2mex_o% Zbot)* &
!!$         & (/ (real(i-1,r8)/real(NZ-1,r8), i=1, NZ) /)

    Xmin = max( (1._r8-deltaW)* Xmin, i2mex_o%Rleft)
    Xmax = min( (1._r8+deltaE)* Xmax, i2mex_o%Rrigh)
    Z0 = (Zmax+Zmin)/2._r8
    Zmin = max( Z0 + (Zmin-Z0)*(1._r8+deltaS), i2mex_o%Zbot )
    Zmax = min( Z0 + (Zmax-Z0)*(1._r8+deltaN), i2mex_o%Ztop )
    Rgrid = Xmin + (Xmax-Xmin)*(/ (real(i-1,r8)/real(NR-1,r8), i=1, NR) /)
    Zgrid = Zmin + (Zmax-Zmin)*(/ (real(i-1,r8)/real(NZ-1,r8), i=1, NZ) /)

    ! now compute br, bz, bphi
    edge_smoothing = 0.05_r8
     call eqm_brz(eqm_ball, edge_smoothing, ier)

 else

    ! -> Assume that psi(r,z) data is not available. So must compute psi(r,z)
    ! from flux data using xplasma routines.
 
    ! set up grid FIRST

    ! grid constrained by plasma size + halo
    Xmin = (1._r8-deltaW)* Xmin
    Xmax = (1._r8+deltaE)* Xmax
    Zmin = min( (1._r8+deltaS)* Zmin, (1._r8-deltaS)* Zmin)
    Zmax = max( (1._r8+deltaN)* Zmax, (1._r8-deltaN)* Zmax)
    Rgrid = Xmin + (Xmax-Xmin)*(/ (real(i-1,r8)/real(NR-1,r8), i=1, NR) /)
    Zgrid = Zmin + (Zmax-Zmin)*(/ (real(i-1,r8)/real(NZ-1,r8), i=1, NZ) /)
    call eqm_bset(i2mex_bphi_cw, i2mex_Ip_cw)
    
    call eqm_cbdy(5,(/Xmin, Xmax, Xmax, Xmin, Xmin/), &
         &          (/Zmin, Zmin, Zmax, Zmax, Zmin/), &
         &          ier)
    if(ier/=0) print*,'**Error after eqm_cbdy', ier

    call eqm_rzgrid(Rgrid,Zgrid,iauto,iauto,NR,NZ,ztol, &
     &                id_R,id_Z,ier)
    if(ier/=0) print*,'**Error after eqm_rzgrid', ier
      

    edge_smoothing = 0.05
    call eqm_brz(eqm_brz_adhoc, edge_smoothing, ier)
 endif
 if(ier/=0) print*,'**Error after eqm_brz', ier

 call i2mex_get_Iboozer(nt1, ns, t, psi, i2d, ier)
 call i2mex_error(ier)
 call i2mex_get_nuboozer(nt1, ns, t, psi, nu2d, ier)
 call i2mex_error(ier)

 call eqdbg_plot

 ! Get the B field = [BR, BZ, Bphi]
 allocate(Rcoord(NR, NZ), Zcoord(NR, NZ), nregion(NR*NZ), &
      & BRBZBPhi(3, NR, NZ), psiRZ(NR, NZ))
 Rcoord = spread(Rgrid, dim=2, ncopies=NZ)
 Zcoord = spread(Zgrid, dim=1, ncopies=NR)
 CALL eqrz_bget(NR*NZ, &
      & Rcoord(1,1), & 
      & Zcoord(1,1), & 
      & spread(0.0_r8, dim=1, ncopies=NR*NZ), & ! Phi
      & ztol,BRBZBPhi,nregion,ier)
 if(ier/=0) print*,'**Error after eqrz_bget', ier

 ! inverse sign!!
 BRBZBPhi = -BRBZBPhi

 ! get Psi on R Z grid
 if (i2mex_direct.eq.1) then
    call eq_fRZ(NR*NZ, Rcoord(1,1), Zcoord(1,1), 1, i2mex_o%id_psirz, &
         NR*NZ, psiRZ(1,1), ier)
 else
    call eq_fRZ(NR*NZ, Rcoord(1,1), Zcoord(1,1), 1, i2mex_o%id_psi, &
         NR*NZ, psiRZ(1,1), ier)
 endif
 if(ier/=0) print*,'**Error after eq_fRZ'

   
 
  xma_ = sum(xa(1:nt1-1,1))/real(nt1-1,r8)
  zma_ = sum(za(1:nt1-1,1))/real(nt1-1,r8)
  b0old = g(1)/xma_
  write(lnout(0),'(a,f10.4)') &
       & ' The toroidal vacuum magnetic field on magnetic axis Bm = ', b0old
       CALL UREAD('Enter new value for Bm$')
       b0new = R8DCOD(IER)
  if(b0new > 0._r8) then
     write(lnout(0),*) ' ** rescale: b0new, b0old = ',b0new,b0old
     ! rescale
     p = b0new**2 * p / b0old**2
     pp = b0new * pp / b0old
     g = b0new* g / b0old
     i2d = b0new* i2d / b0old
     gp = gp
     psi = b0new* psi / b0old
     phi = b0new* phi / b0old
     jpar = b0new* jpar/ b0old
     gpp = b0new**2 * gpp / b0old**2
     xj = b0old * abs(xj) / b0new
     BRBZBPhi = b0new * BRBZBPhi / b0old
     psiRZ = b0new * psiRZ / b0old
  endif
 
  ! save in map01.cdf file
 
  filename = 'map01.cdf'
  write(lnout(0),*)' '
  write(lnout(0),*)' Data will now be saved in file ', filename
 
  call cdfOpn(iu, filename, 'w')
  if(iu==0) then
     print*,'--Oh oh...-- failed to open ', filename
     ier = 1
  endif
 
  dims=(/1, 1, 1/)
  call cdfDefVar(iu, 'mth', dims, 'INT')
  call cdfDefVar(iu, 'nosurf', dims, 'INT')
  call cdfDefVar(iu, 'mx', dims, 'INT')
  call cdfDefVar(iu, 'npsi', dims, 'INT')
  call cdfDefVar(iu, 'kb', dims, 'INT')
  call cdfDefVar(iu, 'nr', dims, 'INT')
  call cdfDefVar(iu, 'nz', dims, 'INT')
  call cdfDefVar(iu, 'time', dims, 'R8')
!!$  call cdfDefVar(iu, 'r', dims, 'R8')
  call cdfDefVar(iu, 'xma', dims, 'R8')
  call cdfDefVar(iu, 'zma', dims, 'R8')
  call cdfDefVar(iu, 'Bma', dims, 'R8')
  call cdfDefVar(iu, 'Bm_experimental', dims, 'R8')
  dims=(/80, 1, 1/)
  call cdfDefVar(iu, 'title', dims, 'CHAR')
  call cdfDefVar(iu, 'date', dims, 'CHAR')
 
  dims=(/ns, 1, 1/)
  call cdfDefVar(iu, 'psibar', dims, 'R8')
  call cdfDefVar(iu, 'phibar', dims, 'R8')
!!$  call cdfDefVar(iu, 'psival', dims, 'R8')
  call cdfDefVar(iu, 'p', dims, 'R8')
  call cdfDefVar(iu, 'g', dims, 'R8')
  call cdfDefVar(iu, 'q', dims, 'R8')
  dims=(/nths, nsf, 1/)
  call cdfDefVar(iu, 'x', dims, 'R8')
  call cdfDefVar(iu, 'z', dims, 'R8')
  call cdfDefVar(iu, 'xjacob', dims, 'R8')
  call cdfDefVar(iu, 'grpssq', dims, 'R8')

  call cdfDefVar(iu, 'I2d', dims, 'R8')   ! I(theta,psi)
  call cdfDefVar(iu, 'nu2d', dims, 'R8')  ! nu(theta,psi) = zeta - phi

  call cdfDefVar(iu, 'xcoord', (/nr,nz,1/), 'R8')
  call cdfDefVar(iu, 'zcoord', (/nr,nz,1/), 'R8')
  call cdfDefVar(iu, 'psixz', (/nr,nz,1/), 'R8')
  call cdfDefVar(iu, 'Bx', (/nr,nz,1/), 'R8')
  call cdfDefVar(iu, 'Bz', (/nr,nz,1/), 'R8')
  call cdfDefVar(iu, 'Bphi', (/nr,nz,1/), 'R8')
  
 
  call cdfPutVar(iu, 'title', i2mex_o%label, ier)
  call cdfPutVar(iu, 'time', time, ier)
  call date_and_time(date=today, time=now)
  lon_date = today(1:4)//'/'//today(5:6)//'/'// &
       & today(7:8)//' at '//now(1:2)//':'//now(3:4)
  call cdfPutVar(iu, 'date', lon_date, ier)
 
  call cdfPutVar(iu, 'mth', nt1-1, ier)
  call cdfPutVar(iu, 'nosurf', ns, ier)
  call cdfPutVar(iu, 'mx', i2mex_mx, ier)
  call cdfPutVar(iu, 'npsi', i2mex_npsi, ier)
  call cdfPutVar(iu, 'kb', i2mex_kb, ier)
  call cdfPutVar(iu, 'nr', NR, ier)
  call cdfPutVar(iu, 'nz', NZ, ier)
!!$  call cdfPutVar(iu, 'r', 1.0_r8, ier) ! we don't rescale so r=1?
  call cdfPutVar(iu, 'xma', xma_, ier)
  call cdfPutVar(iu, 'zma', zma_, ier)
  call cdfPutVar(iu, 'Bma', g(1)/xma_, ier)
  ! old value, before rescaling 
  call cdfPutVar(iu, 'Bm_experimental', b0old, ier)
  ! save both: pol flux in Wb/rad=psibar and in Wb=psival to avoid confusion
  call cdfPutVar(iu, 'psibar', psi, ier)       ! in [Wb/rd]
  call cdfPutVar(iu, 'phibar', phi, ier)       ! in [Wb/rd]
!!$  call cdfPutVar(iu, 'psival', psi*i2mex_twopi_r8, ier) ! in [Wb]
 
  call cdfPutVar(iu, 'p', p, ier)
  call cdfPutVar(iu, 'g', g, ier)
  call cdfPutVar(iu, 'q', q, ier)
  temp = 0.0_r8
  temp(1:nt1, 1:ns) = xa
  call cdfPutVar(iu, 'x', temp, ier)
  temp(1:nt1, 1:ns) = za
  call cdfPutVar(iu, 'z', temp, ier)
  temp(1:nt1, 1:ns) = xj
  call cdfPutVar(iu, 'xjacob', temp, ier)
  temp(1:nt1, 1:ns) = gpp
  call cdfPutVar(iu, 'grpssq', temp, ier)

  temp(1:nt1, 1:ns) = i2d
  call cdfPutVar(iu, 'I2d', temp, ier)
  temp(1:nt1, 1:ns) = nu2d
  call cdfPutVar(iu, 'nu2d', temp, ier)

  if(idebug) then
     allocate(wk(2*max(nt1,ns)),ipgq(nt1,ns),jb2(nt1,ns)); wk=0
     do is=1,ns
        ipgq(1:nt1,is)=g(is)*q(is)+i2d(1:nt1,is)
        jb2(1:nt1,is)=abs(xj(1:nt1,is))*(g(is)*g(is)+gpp(1:nt1,is))/ &
             (xa(1:nt1,is)**2)
     enddo
     write(lnout(0),*) ' ...at plasma bdy:'
     write(lnout(0),'("   theta=",1pe12.5,"; I,g,q=",3(1pe12.5,1x))') &
          t(1),i2d(1,ns),g(ns),q(ns)
     call r8_grf3f2(t,ipgq,jb2,psi,wk,nt1,ns,nt1, &
          'theta','I+g*q','det(J)*B**2','Psi', &
          'rad','T*m','Wb', &
          'mex2eqs debug plot','mex2eqs',0)
     deallocate(wk,ipgq,jb2)
  endif

  call cdfPutVar(iu, 'xcoord', Rcoord, ier)
  call cdfPutVar(iu, 'zcoord', Zcoord, ier)
  call cdfPutVar(iu, 'psixz', psiRZ, ier)
  call cdfPutVar(iu, 'Bx', BrBzBphi(1,:,:), ier)
  call cdfPutVar(iu, 'Bz', BrBzBphi(2,:,:), ier)
  call cdfPutVar(iu, 'Bphi', BrBzBphi(3,:,:), ier)

 
  call cdfCls(iu)

  ! save file for esc 
  call cdfOpn(iu,'esc.nc', 'w')
  call cdfDefVar(iu, 'incr', (/1,1,1/), 'INT')
  call cdfDefVar(iu, 'xbin', (/nt1,1,1/), 'R8')
  call cdfDefVar(iu, 'zbin', (/nt1,1,1/), 'R8')
  call cdfDefVar(iu, 'sin', (/ns,1,1/), 'R8')
  call cdfDefVar(iu, 'pin', (/ns,1,1/), 'R8')
  call cdfDefVar(iu, 'jin', (/ns,1,1/), 'R8')

  call cdfPutVar(iu, 'incr', 0) ! 0=<j.B>/R<B.gradPhi>
  call cdfPutVar(iu, 'xbin', xa(1:nt1,ns))
  call cdfPutVar(iu, 'zbin', za(1:nt1,ns))
  call cdfPutVar(iu, 'sin', sqrt(phi/maxval(phi)))
  call cdfPutVar(iu, 'pin',  p)
  call cdfPutVar(iu, 'jin',  jpar/r0)
  call cdfCls(iu)
  
  
 
  DEALLOCATE(t)
  DEALLOCATE(psii)
  DEALLOCATE(psi, p, pp, q, g, gp, phi, jpar)
  DEALLOCATE(bsquare, oneOverRSquare)
  DEALLOCATE(xa, za, xj, i2d, nu2d)
  DEALLOCATE(gtt, gtp, gpp)
  DEALLOCATE(gserror)
  DEALLOCATE(temp)
  deallocate(Rgrid, Zgrid)
  deallocate(Rcoord, Zcoord, nregion, &
      & BRBZBPhi, psiRZ)
 
  call i2mex_free(ier)
 
  print *, 'All done.'
 
end program
