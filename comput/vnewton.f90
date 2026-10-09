subroutine vnewton(xmin,xmax,ffpsub,nvec,xvec,fvec,tol,nauxi,auxi, &
     nauxo,auxo,ier)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero
  implicit none

  ! compute a set of solutions { f(xvec(j))=fvec(j) }

  real(fp), intent(in) :: xmin,xmax  ! (xmin,xmax) as 2-vector, xmax .ne. xmin

  external ffpsub  ! f,f' function evaluation routine (see description below)

  integer, intent(in) :: nvec      ! solution vector size
  real(fp), intent(inout) :: xvec(nvec)! initial/final solution estimates(in/out)
  real(fp), intent(in) :: fvec(nvec) ! function values sought

  real(fp), intent(in) :: tol        ! completion criterion:
  !   change in x < tol*(xmax-xmin)

  !  for ffpsub:

  integer :: nauxi                 ! # of auxilliary inputs per vector element
  real(fp), intent(in) :: auxi(nauxi,nvec)

  integer :: nauxo                 ! # of auxilliary outputs per vector element
  real(fp), intent(in) :: auxo(nauxo,nvec)

  integer, intent(out) :: ier      ! completion code (0=OK)

  !-----------------
  !  user provides:
  !  subroutine <fppsub>(nvec,xvec,feval,dfdx,nauxi,auxi,nauxo,auxo)
  !    (nvec,xvec,nauxi,auxi,nauxo,auxo) as declared above).
  !    real(fp) feval(nvec)  ! f evaluation at each pt (to compare with fvec).
  !    real(fp) dfdx(nvec)   ! df/dx at each evaluation point.
  !-----------------

  real(fp) :: xlo(nvec),xhi(nvec),fxlo(nvec),fxhi(nvec)
  real(fp) :: feval(nvec),dfdx(nvec)

  real(fp) :: dxmaxa,df,dx(nvec),ftest

  integer, parameter :: itmax=50

  integer :: i,iter,isplit
  logical ihave_lims

  !-------------------------------------------------------
  !  rudimentary error checks

  ier=0

  if(xmax.le.xmin) then
    ier=1
    return
  endif

  do i=1,nvec
    if((xvec(i).lt.xmin).or.(xvec(i).gt.xmax)) ier=2
  enddo
  if(ier.gt.0) return

  xlo=xmin
  xhi=xmax

  iter=0

  do
    ! evaluate vector of f, df/dx

    call ffpsub(nvec,xvec,feval,dfdx,nauxi,auxi,nauxo,auxo)

    ! compute Newton steps & maximum step; recover from out of range steps

    dxmaxa=ZERO

    do i=1,nvec
      isplit=0
      df=fvec(i)-feval(i)
      if(df.eq.ZERO) cycle
      if((xmax-xmin)*abs(dfdx(i)).lt.abs(df)) then
        isplit=1
      else
        dx(i)=df/dfdx(i)
        if(xvec(i)+dx(i).gt.xhi(i)) then
          isplit=1
        else if(xvec(i)-dx(i).lt.xlo(i)) then
          isplit=1
        else
          dxmaxa=max(dxmaxa,abs(dx(i)))
        endif
      endif

      if(isplit.eq.1) then

        ! Newton step out of range, fall back to binary search

        if(.not.ihave_lims) then
          ihave_lims=.TRUE.
          call ffpsub(nvec,xlo,fxlo,dfdx,nauxi,auxi,nauxo,auxo)
          call ffpsub(nvec,xhi,fxhi,dfdx,nauxi,auxi,nauxo,auxo)
        endif
        ftest=(fvec(i)-fxlo(i))*df
        if(ftest.eq.ZERO) then
          dx(i)=xlo(i)-xvec(i)  ! match at lo endpoint
          cycle
        else if(ftest.lt.ZERO) then

          ! target lies between last evaluation point and xlo
          xhi(i)=xvec(i)

        else

          ftest=df*(fvec(i)-fxhi(i))
          if(ftest.eq.ZERO) then
            dx(i)=xhi(i)-xvec(i)  ! match at hi endpoint
            cycle
          else

            ! target lies between last evaluation point and xhi
            xlo(i)=xvec(i)

          endif
        endif
        dx(i)=(xlo(i)+xhi(i))/2 - xvec(i)
        dxmaxa=max(dxmaxa,abs(dx(i)))
      endif

    enddo

    ! OK ready for next step...

    if(dxmaxa.gt.tol*(xmax-xmin)) then

      iter=iter+1
      if(iter.gt.itmax) then

        ier=3  ! convergence failure
        return

      endif

      xvec = xvec + dx

    else

      exit  ! convergence OK

    endif

  enddo  ! end of Newton convergence loop

  return
end subroutine vnewton
 
subroutine vnewton_error(lun,ier)

  implicit none

  !  echo error message (immediate exit, no message, iff ier=0)

  integer, intent(in) :: lun   ! Fortran LUN where to write message
  integer, intent(in) :: ier   ! error code

  if(ier.eq.0) then
    return

  else if(ier.eq.1) then
    write(lun,*) &
         ' ?vnewton: xmax .le. xmin'

  else if(ier.eq.2) then
    write(lun,*) &
         ' ?vnewton: one or more xvec(...) values outside range [xmin,xmax]'

  else if(ier.eq.3) then
    write(lun,*) &
         ' ?vnewton: for one or more xvec(...) max iteration count exceeded.'

  else
    write(lun,*) ' ?vnewton:  error (unspecified).'
  endif

  write(lun,*) '  (error code ',ier,')'

  return
end subroutine vnewton_error
