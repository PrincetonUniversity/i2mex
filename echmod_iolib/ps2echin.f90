subroutine ps2echin(ss,iant,idamp,inray,iout,filenam,ierr)

  ! from a Plasma State (ss) create an echin file for TORAY, for ECH
  ! antenna #iant.

  ! any messages are written on unit iout; status code returned in ierr.

  use plasma_state_mod
  implicit NONE

  type (plasma_state) :: ss

  integer, intent(in) :: iant  ! ECH antenna number
  integer, intent(in) :: idamp ! damping model number
  integer, intent(in) :: inray ! number of rays to be traced, this antenna

  integer, intent(in) :: iout  ! fortarn I/O unit number for messages
  character*(*), intent(in) :: filenam  ! filename.  If blank, use "echin"

  integer, intent(out) :: ierr ! return status code, 0=OK

  !--------------------------------------
  integer :: inprof,inm1,io,ios,inbfld,inmin,istat,inside

  ! profiles interpolated to TORAY grid:
  real*8 :: znorm,zbzxr,zr0,zb0,zrmins,ztflux
  real*8 :: zpi
  real*8, dimension(:), allocatable :: zne_tor,zte_tor,zef_tor,zplflx_tor,zwk
  real*8, dimension(:), allocatable :: znewk

  real*8, parameter :: zgfac = 0.10d0  ! for patch routine

  character*100 zfilenam
  !--------------------------------------

  call find_io_unit(io)

  if(filenam.eq.' ') then
     zfilenam='echin'
  else
     zfilenam=filenam
  endif

  ierr = 0

  zpi = 3.1415926535897931D+00

  !------------------------------
  ! interpolate input profiles...

  inprof = ss%nrho_ecrf
  inm1 = inprof-1

  allocate(zne_tor(inprof),zte_tor(inprof),zef_tor(inprof),zplflx_tor(inprof))
  allocate(zwk(inm1))

  ! interpolate the Psi(rho) 1d spline...
  call ps_intrp_1d_vec(ss%rho_ecrf, ss%id_psipol, zplflx_tor, ierr, &
       state=ss)
  if(ierr.ne.0) then
     write(iout,*) ' ?ps2echin: Psi(rho) spline interpolation error.'
     return
  endif

  ! convert to normalized poloidal flux profile

  znorm = 1.0d0/zplflx_tor(inprof)
  zplflx_tor(1:inprof) = zplflx_tor(1:inprof)*znorm

  ! smoothed rezone of plasma parameter profiles...
  !   rezone to zone ctr then shift to bdys... CONTAINED subroutine...

  call rezone_intrp(ss%id_ns(0),zne_tor)  ! density...
  call rezone_intrp(ss%id_Ts(0),zTe_tor)  ! temperature: KeV ok
  call rezone_intrp(ss%id_Zeff,zef_tor)   ! Zeff

  zne_tor = 1.0d-6*zne_tor  ! -> cm**-3

  !------------------------------
  ! open the file

  open(unit=io,file=zfilenam,status='unknown',iostat=ios)
  if(ios.ne.0) then
     write(iout,*) ' ?ps2echin: "'//trim(zfilenam)//'" file open failure.'
     ierr=1
     return
  endif

  !------------------------------
  ! scalars...

  zbzxr = ss%g_eq(ss%nrho_eq)  ! vacuum R*B_phi (always positive), T*m
  zr0 = ss%r_axis              ! axis location, m
  zb0 = zbzxr/zr0              ! vacuum |B_phi|, T

  ztflux = ss%phit(ss%nrho_eq) ! enclosed toroidal flux, Wb
  zrmins = 100.0d0*sqrt(ztflux/(zpi*zb0)) ! nominal radius, cm

  zr0 = 100.0d0*zr0   ! ->cm
  zb0 = 1.0d4*zb0     ! ->Gauss

  inbfld = 3

  write(io, 1001) ss%t0
  write(io, 1000) idamp,inprof,inray,inbfld

  !RGA, give toray the same sign of B0 as in the psiin file.  Correct
  !     driven current sign in ps2trcom_ech()

  if (ss%kccw_Bphi<0) zb0=-zb0

  write(io, 1001) ss%freq_ec(iant), ss%ec_Omode_Fraction(iant), &
       100.0d0*ss%R_ec_launch(iant), 100.0d0*ss%Z_ec_launch(iant), &
       ss%ec_theta_aim(iant), ss%ec_phi_aim(iant), &
       ss%ec_half_power_angle(iant), ss%ec_beam_elongation(iant), &
       zr0,zb0,zrmins

  ! DMC Dec 2010: patch density to prevent dne/dx > 0 in edge

  if(allocated(znewk)) deallocate(znewk)
  allocate(znewk(inprof))
  znewk = zne_tor

  call find_inmin(ss%rho_ecrf,znewk,inside,inmin)

  call makedens_gneg(inprof,ss%rho_ecrf,znewk,zne_tor,zgfac, &
       inside,inprof,inmin, iout, istat)


  if(istat.ne.0) then
     call errmsg_exit(' ?ps2echin: error in density patch, makedens_gneg!')
  endif

  write(io, 1001) zplflx_tor(1:inprof)
  write(io, 1001) zef_tor(1:inprof)
  write(io, 1001) zne_tor(1:inprof)
  write(io, 1001) zte_tor(1:inprof)

1000 format (20i4)
1001 format (5e16.9)

  close(unit=io)

CONTAINS
  subroutine rezone_intrp(id,zans)
    integer, intent(in) :: id
    real*8, dimension(:) :: zans

    ! smoothed rezone to zone ctrs in "zwk"

    !------------
    real*8 :: zlim1,zlim2
    integer :: ix
    !------------

    call ps_rho_rezone(ss%rho_ecrf, id, zwk, ierr, &
         state=ss, zonesmoo=.TRUE.)

    ! interpolation to interior boundaries

    do ix=2,inprof-1
       zans(ix)=0.5d0*(zwk(ix-1)+zwk(ix))
    enddo

    ! constrained parabolic extrapolation to boundaries
    ! note: all profiles are positive quantities: {n,T,Zeff}

    zlim1=0.5d0*zwk(1)
    zlim2=1.5d0*zwk(1)

    zans(1)=(9*zwk(1)-zwk(2))/8
    zans(1)=max(zlim1,min(zlim2,zans(1)))

    zlim1=0.5d0*zwk(inm1)
    zlim2=1.5d0*zwk(inm1)

    zans(inprof)=(9*zwk(inm1)-zwk(inm1-1))/8
    zans(inprof)=max(zlim1,min(zlim2,zans(inprof)))

  end subroutine rezone_intrp

  subroutine find_inmin(zx,zne,ix1,inmin)

    ! find first point that is both inside x=0.9 and inboard of a point
    ! which is different from the edge value by more than 1/5 * (nmax-nmin)

    real*8, intent(in), dimension(:) :: zx,zne
    integer, intent(out) :: inmin

    integer :: ix,ix1,istart
    real*8 :: znmin,znmax,zdeln

    logical :: idiff_found,il_gneg

    do ix1=inprof-1,2,-1
       if(zx(ix1).le.0.9d0) exit
    enddo
    il_gneg=((zne(ix1)-zne(2)).lt.0d0) !Rough estimation of gradient sign: 
                                      !.true. negative gradient in 0<x<0.9

    if(il_gneg) then
       istart=2
       ix1=1
    else
       istart=ix1
    endif
    znmin=zne(ix1)
    znmax=znmin
    do ix=istart,inprof
       znmin=min(znmin,zne(ix))
       znmax=max(znmax,zne(ix))
    enddo
    
    zdeln = (znmax-znmin)/5
    
    idiff_found = .FALSE.
    
    if(il_gneg) then
       do ix=inprof-1,istart,-1
          if(abs(zne(ix)-zne(inprof)).gt.zdeln) idiff_found = .TRUE.
          if(idiff_found.AND.(zx(ix).le.0.9d0)) then
             inmin=ix
             exit
          endif
       enddo
    else
       do ix=inprof-1,istart,-1
          if(abs(zne(ix)-zne(inprof)).gt.zdeln) idiff_found = .TRUE.
          if(idiff_found) then
             inmin=ix
             exit
          endif
       enddo
    endif

  end subroutine find_inmin

end subroutine ps2echin

subroutine echin_ck_dampmod(nant,ndamp_in,ndamp_use,iout)

  implicit NONE

  ! check damping model selection

  integer, intent(in) :: nant   ! number of antennas
  integer, intent(in) :: ndamp_in(nant)   ! input damping model selection
  integer, intent(out) :: ndamp_use(nant) ! output selections to be used
  integer, intent(in) :: iout   ! Fortran I/O unit for messages
  
  ! these rules are based on the presumption that the same damping model
  ! is desirable for all antennas, in most simulations.

  ! Rules:

  !   (a) find the first antenna index iant such that ndamp_in(iant).gt.0;
  !       if no such index exists set idamp=0 otherwise idamp=ndamp_in(iant).
  !   (b) for all indices jant, the selected damping model is:
  !           idamp if ndamp_in(jant).gt.0
  !           abs(ndamp_in(jant)) if ndamp_in(jant).le.0

  !----------------
  integer :: iant,jant,idamp,iwarn
  !----------------

  ndamp_use = 0

  idamp = 0
  iwarn = 0

  do iant=1,nant
     if(ndamp_in(iant).gt.0) then
        idamp=ndamp_in(iant)
        exit
     endif
  enddo

  do jant=1,nant
     if(ndamp_in(jant).gt.0) then
        ndamp_use(jant)=idamp
     else
        ndamp_use(jant)=abs(ndamp_in(jant))
     endif
     if(ndamp_use(jant).ne.ndamp_in(jant)) iwarn = iwarn + 1
  enddo

  if(iwarn.gt.0) then
     write(iout,*) ' %echin_ck_dampmod: damping models adjusted: '
     write(iout,*) '  antenna#  input     use'
     do jant=1,nant
        write(iout,'(3(4x,i6))') jant,ndamp_in(jant),ndamp_use(jant)
     enddo
  endif

end subroutine echin_ck_dampmod
