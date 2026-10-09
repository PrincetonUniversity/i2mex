subroutine r8filtr90(x,y1,y,n,epslon,delta,iddf,ddv,iblf,blv,xenda,xendb)
  !
  !  a simplified interface to filtr6 -- fortran-90
  !
  !  all of the arguments correspond to filtr6.for arguments (see the
  !  extensive comments in filtr6.for).
  !
  !  the following simplifications are present:
  !    1.  epslon is a scalar; epslon=0.0 => epslon feature not used;
  !        there is no "eps2" feature.
  !    2.  delta is a scalar
  !    3.  the output quantities sm,ndrop,nblrs are dropped.
  !
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: ZERO, HALF, ONE
  implicit none
  ! Arguments
  integer, intent(in) :: n
  real(fp), dimension(n), intent(in)  :: x  ! indep. coordinate (input), x must be strict ascending
  real(fp), dimension(n), intent(in)  :: y1 ! data to be smoothed (input)
  real(fp), dimension(n), intent(out) :: y  ! the data as smoothed (output)
  real(fp), intent(in) :: epslon            ! limit on abs(y(i)-y1(i)), epslon.le.0  =>  epslon = infinity
  real(fp), intent(in) :: delta             ! smoothing parameter
  !  on output each y(j) is a triangular weighted average of the piecewise
  !  linear interpolant of y1 over the range x(j)-delta to x(j)+delta.
  !
  integer,  intent(in) :: iddf              ! data drop flag (=1 to enable)
  real(fp), intent(in) :: ddv               ! drop criterion
  integer,  intent(in) :: iblf              ! baseline flag (=1 to enable)
  real(fp), intent(in) :: blv               ! baseline value
  real(fp), intent(in) :: xenda             ! LHS boundary condition
  real(fp), intent(in) :: xendb             ! RHS boundary condition
  !
  !  xend[a,b]=0.0 => endpoint is fixed, e.g. y(1)=y1(1) on output
  !  xend[a,b]=1.0 => dy/dx --> 0 as x--> the boundary
  !  intermediate values => intermediate action at end point.
  !
  ! Local variables
  real(fp) :: sm
  integer :: ndrop,nblrs
  !
  real(fp) :: epslim
  real(fp), parameter :: r8large = 1.0e36_fp
  real(fp), dimension(:), allocatable :: epsarr, delarr
  !
  integer :: stat_alloc
  !
  ! Test for NaNs in input data
  if (any(y1(:)/=y1(:))) then
     write(6,*) ' !filtr90: warning: NaNs detected in data, cannot smooth!'
     where (y1(:)==y1(:))
        y = y1
     elsewhere
        y = 0.0
     endwhere
     return
  endif
  !
  allocate(epsarr(n),stat=stat_alloc)
  if(stat_alloc.ne.0) then
    write(6,*) ' ?filtr90:  epsarr array allocation failed!'
    call flush(6)
    stop
  endif
  !
  allocate(delarr(n),stat=stat_alloc)
  if(stat_alloc.ne.0) then
    write(6,*) ' ?filtr90:  delarr array allocation failed!'
    call flush(6)
    stop
  endif
  !
  epslim=r8large
  if(epslon.le.ZERO) then
    epsarr=epslim
  else
    epsarr=epslon
  endif
  !
  delarr=delta
  !
  call r8filtr6(x,y1,y,n,epsarr,n,epsarr,0,delarr,iddf,ddv, &
       iblf,blv,xenda,xendb,sm,ndrop,nblrs)
  !
  deallocate(epsarr)
  deallocate(delarr)
  !
  return
end subroutine r8filtr90
