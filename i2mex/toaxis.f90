subroutine i2mex_toAxis(nt1, ns, the, psi, power, f, ier)

  ! Cubic extrapolation to axis for 1st, psi(1) node. 
  ! *NOTE* Need to have at least 5 nodes

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
  real(i2mex_r8), intent(in) :: power ! f ~ psi**power as psi-> 0
  real(i2mex_r8), intent(inout) :: f(nt1, ns)
  integer, intent(out) :: ier

  real(i2mex_r8), parameter :: SMALL=1.e-6, LARGE=1.e+6
  integer i

  include 'cubinterp.h'

  ier = 0
  if(ns <5) return ! don't try to extrapolate if too few points are available ier = 144

  if(psi(1) == 0._i2mex_r8) then

     if(power > 0._i2mex_r8) then

        ! goes to zero

        f(1:nt1, 1) = SMALL

     else if(power < 0._i2mex_r8) then

        ! goes to infinity

        f(1:nt1, 1) = LARGE

     else

        ! finite

        do i=1, nt1
           f(i,1) = FCCCC0( &
                & f(i,2), &
                & f(i,3), &
                & f(i,4), &
                & f(i,5), &
                &                 psi(2),  psi(3),  psi(4),  psi(5), &
                & psi(1))
        enddo

     endif

  else

     ! apply scaling

     do i=1, nt1

        f(i,1) = psi(1)**power * FCCCC0( &
             & f(i,2)/psi(2)**power, &
             & f(i,3)/psi(3)**power, &
             & f(i,4)/psi(4)**power, &
             & f(i,5)/psi(5)**power, &
             &                 psi(2),  psi(3),  psi(4),  psi(5), &
             & psi(1))
     enddo

  endif

end subroutine i2mex_toAxis

subroutine i2mex_toAxisFraction(nt1, ns, the, psi, power, fraction, f, ier)

  ! Cubic extrapolation to axis for first fraction of points
  ! fraction should be < 0.5

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
  real(i2mex_r8), intent(in) :: power ! f ~ psi**power as psi-> 0
  real(i2mex_r8), intent(in) :: fraction ! fraction of grid where fct is to 
  ! be extrapolated
  real(i2mex_r8), intent(inout) :: f(nt1, ns)
  integer, intent(out) :: ier

  real(i2mex_r8), parameter :: SMALL=1.e-6, LARGE=1.e+6
  integer i, j
  integer :: npts = 1

  ier = 0
  if(fraction == 0.0_i2mex_r8) return


  if(fraction >=0.5_i2mex_r8) ier = 155

  npts = int(fraction * real(ns, i2mex_r8))

  if(npts+4 > ns) then
     ier = 161
     return
  endif

  if(power == 0.0_i2mex_r8) then

     do j = 1, npts
        do i = 1, nt1
           f(i,j) = f(i,npts+1)
        enddo
     enddo

  else

     do j = 2, npts
        do i = 1, nt1
           f(i,j) = psi(j)**power *f(i,npts+1)/psi(npts+1)**power
        enddo
     enddo

     if(psi(1) == 0._i2mex_r8) then
        if(power > 0._i2mex_r8) then
           ! goes to zero
           f(1:nt1, 1) = SMALL
        else if(power < 0._i2mex_r8) then
           ! goes to infinity
           f(1:nt1, 1) = LARGE
        endif
     endif

  endif



end subroutine i2mex_toAxisFraction

subroutine i2mex_toEdgeFraction(nt1, ns, the, psi, fraction, f, ier)

  ! Cubic extrapolation to edge for last fraction of points
  ! fraction should be < 0.5

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
  real(i2mex_r8), intent(in) :: fraction ! fraction of grid where fct is to 
  ! be extrapolated
  real(i2mex_r8), intent(inout) :: f(nt1, ns)
  integer, intent(out) :: ier

  integer i, j
  integer :: npts = 1

  include 'cubinterp.h'

  ier = 0
  if(fraction == 0.0_i2mex_r8) return

  if(fraction >=0.5_i2mex_r8) ier = 158

  npts = int(fraction * real(ns, i2mex_r8))

  if(npts+4 > ns) then
     ier = 157
     return
  endif

  do j = 0, npts-1

     do i = 1, nt1

        f(i,ns-j) = FCCCC0( &
             & f(i,ns-npts-1), &
             & f(i,ns-npts-3), &
             & f(i,ns-npts-5), &
             & f(i,ns-npts-7), &
             & psi(ns-npts-1), psi(ns-npts-3), psi(ns-npts-5), psi(ns-npts-7), &
             & psi(ns-j))
     enddo

  enddo


end subroutine i2mex_toEdgeFraction
