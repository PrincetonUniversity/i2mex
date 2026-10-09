subroutine zriddery (ivec, iok0, x1, x2, eps, eta, subr, x, ifail, &
     ivecd,zinput,ninput,zoutput,noutput)
  !
  !  *** see also zridderx.for ***
  !  this is a copy of the zridderx code -- use this instead of zridderx,
  !  if the external routine is itself likely to generate a zridderx call
  !  (avoid seg faults due to inadvertant and unsupported recursion).
  !
  !  **vectorized dmc Apr 2000**
  !  real(fp) root finder -- extended to allow extra arguments to
  !  subr:  zinput(ninput) & zoutput(noutput) [dmc 4 Nov 1999]
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  !
  integer :: ivec                    ! vector size
  integer :: ifail               ! error code
  real(fp), dimension(ivec) :: x1
  real(fp), dimension(ivec) :: x2
  real(fp) :: eps, eta
  real(fp), dimension(ivec) :: x
  logical,  dimension(ivec) :: iok0  ! mask -- .true. to skip element
  external subr
  integer :: ivecd
  integer :: ninput,noutput
  real(fp), dimension(ivecd,ninput)  :: zinput
  real(fp), dimension(ivecd,noutput) :: zoutput
  !
  !     returns an approximate zero for the function (cf subr) given an
  !     initial interval (x1,x2) which brackets a root.
  !     The approximate zero is determined so that x is within the accuracy
  !     of eps or |func(x)| < eta.
  !     The algorithm is based on Ridders' method
  !     see Numerical Recipes 2nd ed. p351
  !
  !  external subroutine header:
  !     subroutine subr(ivec,iok,x,ansr,ivecd,zinput,ninput,zoutput,noutput)
  !        integer ivec
  !        logical iok(ivec)  ! if iok(i) TRUE then SKIP element i of vector
  !        real(fp) x(ivec)     ! input vector {x(i)}
  !        real(fp) ansr(ivec)  ! output vector {f(x(i))}
  !        integer ivecd      ! auxilliary information vector size (.ge.ivec)
  !        integer ninput,noutput
  !        real(fp) zinput(ivecd,ninput),zoutput(ivecd,noutput)
  !
  ! ----------------------------------------------------------------------
  !
  real(fp), dimension(:), allocatable :: xl,xh,fl,fh,xm,fm,xnew,fnew,s
  logical, dimension(:), allocatable :: iok
  !
  integer nc
  !
  integer, parameter :: itermax = 100
  real(fp), parameter :: rbig = 1.024d99
  integer, parameter :: stderr = 6
  integer :: iter,imsg,i,isave
  !
  character(len=80), dimension(3) :: errmsg
  errmsg(1) = 'ZRIDDER: root not bracketed'
  errmsg(2) = 'ZRIDDER: anomalous behavior of func=0.(?)'
  errmsg(3) = 'ZRIDDER: exceeded maximum (100) iterations'
  !
  !---------------------------------------------------------
  allocate(iok(ivec))
  allocate(xl(ivec),xh(ivec),fl(ivec),fh(ivec))
  allocate(xm(ivec),fm(ivec),xnew(ivec),fnew(ivec),s(ivec))
  !
  do i=1,ivec
    iok(i)=iok0(i)
    xl(i) = x1(i)
    xh(i) = x2(i)
  end do
  !
  !  eval at endpts of intervals; check for acceptance of endpts...
  !
  call subr(ivec,iok,xl,fl,ivecd,zinput,ninput,zoutput,noutput)
  do i=1,ivec
    if(.not.iok(i)) then
      if (ABS (fl(i)) .le. eta) then
        iok(i)=.TRUE.
        x(i)     = xl(i)
      end if
    end if
  end do
  !
  call subr(ivec,iok,xh,fh,ivecd,zinput,ninput,zoutput,noutput)
  !
  nc=0
  imsg=0
  do i=1,ivec
    if(.not.iok(i)) then
      if (ABS (fh(i)) .le. eta) then
        iok(i)=.TRUE.
        x(i)     = xh(i)
      else
        x(i) = rbig
        nc=nc+1
        if((fl(i).gt.0.d0).and.(fh(i).gt.0.d0).or. &
             (fl(i).lt.0.d0).and.(fh(i).lt.0.d0)) then
          imsg=1
          isave=i
        end if
      end if
    end if
  end do
  if(imsg.eq.1) go to 9999
  if(nc.eq.0) then
    ifail=0
    go to 10000                    ! all done
  end if
  !
  do iter=1,itermax
    do i=1,ivec
      if(.not.iok(i)) then
        xm(i)=(xl(i)+xh(i))/2.d0
      end if
    end do
    call subr(ivec,iok,xm,fm,ivecd,zinput,ninput,zoutput,noutput)
    do i=1,ivec
      if(.not.iok(i)) then
        s(i)  = SQRT(max(0.d0,(fm(i)**2-fl(i)*fh(i))))
        if (s(i) .eq. 0.d0) then
          isave=i
          imsg = 2
          x(i) = -rbig
        else
          xnew(i)=xm(i)+(xm(i)-xl(i))*fm(i)/s(i) * &
               SIGN (1.d0,fl(i)-fh(i))
        end if
      end if
    end do
    call subr(ivec,iok,xnew,fnew,ivecd,zinput,ninput,zoutput,noutput)
    nc=0
    do i=1,ivec
      if(.not.iok(i)) then
        if (ABS (x(i)-xnew(i)) .le. eps .or. &
             ABS (fnew(i)) .le. eta) then
          iok(i)=.TRUE.
        end if
        x(i) = xnew(i)
      end if
      if(.not.iok(i)) then
        if (SIGN (fm(i),fnew(i)) .ne. fm(i)) then
          xl(i) = xm(i)
          fl(i) = fm(i)
          xh(i) = xnew(i)
          fh(i) = fnew(i)
        else if (SIGN (fl(i),fnew(i)) .ne. fl(i)) then
          xh(i) = xnew(i)
          fh(i) = fnew(i)
        else if (SIGN (fh(i),fnew(i)) .ne. fh(i)) then
          xl(i) = xnew(i)
          fl(i) = fnew(i)
        else
          call errmsg_exit('subroutine ZRIDDERY: programming error')
        end if
        if (ABS (xl(i)-xh(i)) .le. eps) then
          iok(i)=.TRUE.
        else
          nc=nc+1
        end if
      end if
    end do
    if(nc.eq.0) then
      ifail=0
      go to 10000
    end if
  end do
  !
  !  get here if for some point, convergence failed
  !
  imsg = 3
  do i=1,ivec
    if(.not.iok(i)) isave=i
  end do
  !
9999 continue
  if((imsg.gt.0).and.(imsg.le.3)) then
    write (stderr, '(a)')  errmsg(imsg)
    write (stderr, '(a)') 'i xl(i), fl(i), xh(i), fh(i):'
    write (stderr,   *  )  isave, xl(isave), fl(isave), &
         xh(isave), fh(isave)
  else
    call errmsg_exit('subroutine ZRIDDERY: unspecified error')
  end if
  !
  ifail = imsg
  !
10000 continue
  deallocate(xl,xh,fl,fh,xm,fm,xnew,fnew,s,iok)
  return
end subroutine zriddery
