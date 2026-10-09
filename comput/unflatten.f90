subroutine unflatten(darray,nd)

  !--------------------------

  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero
  implicit none

  integer, intent(in) :: nd   ! size of darray
  real(fp), intent(inout) :: darray(nd)   ! data to be "unflattened"

  !--------------------------

  ! without affenting minval(darray) or maxval(darray), introduce
  ! slight tilts to any flat spots, to prevent any occurrence of
  ! darray(j)-darray(j-1) = 0 precisely.

  ! this can be useful as protection when communicating a profile
  ! from a "REAL" code to a "REAL(FP)" code which cannot tolerate
  ! flat spots or "zero gradient regions".  In this case, when
  ! there are differences between two successive REAL*4 numbers, there
  ! is lots of room for introducing a small tilt in REAL(FP) space.

  ! example of effect:

  ! 1.0000000 2.0000000 2.0000000 2.0000000 1.0000000
  !   becomes
  ! 1.0000000 1.9999999 2.0000000 1.9999999 1.0000000

  ! 1.0000000 2.0000000 2.0000000 2.0000000 3.0000000
  !   becomes
  ! 1.0000000 1.9999999 2.0000000 2.0000001 3.0000000

  !--------------------------
  ! local variables:

  real(fp) :: global_min,global_max,global_diff

  real(fp) :: anchor,end_right,end_left,local_diff,incr,incr1,incr2

  integer :: i,lflat,nflat,lanchor

  real(fp), parameter :: ZSMALL = 1.0e-5_fp
  real(fp), parameter :: ZTINY = 1.0e-9_fp

  !--------------------------
  ! executable code:

  global_max=maxval(darray)
  global_min=minval(darray)
  global_diff = global_max - global_min

  if(global_diff.eq.ZERO) return

  !--------------------------
  ! check left bdy region

  nflat=1
  i=1
  do
    i=i+1
    if(darray(i).ne.darray(i-1)) exit
    nflat=nflat+1
  enddo

  if(nflat.gt.1) then

    local_diff=darray(i)-darray(i-1)
    incr=ZSMALL*local_diff
    incr=min(ZTINY*global_diff,max(-ZTINY*global_diff,incr))/ &
         (nflat-1)

    do i=2,nflat
      darray(i)=darray(i-1)+incr
    enddo
  endif

  !--------------------------
  ! check right bdy region

  nflat=1
  i=nd
  do
    i=i-1
    if(darray(i).ne.darray(i+1)) exit
    nflat=nflat+1
  enddo

  if(nflat.gt.1) then

    local_diff=darray(i)-darray(i+1)
    incr=ZSMALL*local_diff
    incr=min(ZTINY*global_diff,max(-ZTINY*global_diff,incr))/ &
         (nflat-1)

    do i=nd-1,nd-nflat+1,-1
      darray(i)=darray(i+1)+incr
    enddo
  endif

  !--------------------------
  ! check for interior flat regions

  lflat=1
  do
    lflat=lflat+1
    if(lflat.ge.nd) exit

    if(darray(lflat).ne.darray(lflat+1)) cycle

    !----------------------------------
    ! start of flat region detected...
    !----------------------------------

    nflat=1
    i = lflat
    do
      i=i+1
      if(darray(i).ne.darray(i-1)) exit
      nflat=nflat+1
    enddo

    ! darray(lflat:lflat+nflat-1) is "flat"; nflat.ge.2 to get here...
    ! find "middle" of flat spot

    lanchor = lflat + nflat/2
    anchor=darray(lanchor)

    incr1=anchor-darray(lflat-1)      ! guaranteed non-zero
    incr2=darray(lflat+nflat)-anchor  ! ditto

    if(incr1*incr2.gt.zero) then

      ! incr1 & incr2 of same sign:
      ! flat point is an "inflection"

      if(incr1.gt.ZERO) then
        incr=min(incr1,incr2)  ! both positive, take smaller modulus
      else
        incr=max(incr1,incr2)  ! both negative, take smaller modulus
      endif
      incr=ZSMALL*incr
      incr=min(ZTINY*global_diff,max(-ZTINY*global_diff,incr))/ &
           (nflat/2)

      !  add a tilt to the flat region, in the direction set by the
      !  neighbouring points

      do i=lanchor-1,lflat,-1
        darray(i)=darray(i+1)-incr
      enddo

      do i=lanchor+1,lflat+nflat-1,1
        darray(i)=darray(i-1)+incr
      enddo

    else

      ! flat spot is local maximum (if incr1.gt.0) or minimum

      if(incr1.gt.ZERO) then
        incr=min(incr1,-incr2)         ! local maximum (incr > 0)
      else
        incr=max(incr1,-incr2)         ! local minimum (incr < 0)
      endif
      incr=ZSMALL*incr
      incr=min(ZTINY*global_diff,max(-ZTINY*global_diff,incr))/(nflat/2)

      !  from the anchor point tilt down/up in both directions according
      !  as this is a local maximum/minimum

      do i=lanchor-1,lflat,-1
        darray(i)=darray(i+1)-incr
      enddo

      do i=lanchor+1,lflat+nflat-1,1
        darray(i)=darray(i-1)-incr
      enddo

    endif

    lflat = lflat + nflat - 1

  enddo

  return
end subroutine unflatten
