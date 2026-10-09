subroutine chk_interpolate(zprof,xi,nx)
  !
  !  patch "holes" (zero values) in profile,
  !  then smooth.
  !
  !  done for MC radial profiles that are formed by summing of the
  !  profile itself and its normalization
  !
  !    zprof(j) = [MC sum]f*dw
  !    zwsum(j) = [MC sum]dw
  !
  !  and at normalization time
  !
  !    if (zwsum(j).gt.0.0) then
  !         zprof(j)=zprof(j)/zwsum(j)
  !    else
  !         zprof(j)=0.0
  !    endif
  !
  !  which can leave "holes" in the profile in case of lousy statistics.
  !
  !  *** SO ***
  !  fill these holes with linear extrapolation & linear interpolation
  !
  !-----------------------------
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  integer :: nx
  real(fp), dimension(nx) :: zprof,xi
  !
  !
  real(fp) :: zx1,zf1,zx2,zf2,zx,zf
  !
  integer :: j,jnz,jj,isrch,jsave1,jsave2
  !
  !-----------------------------
  !
  !  1.  fill in from left
  !
  do j=1,nx
    if(zprof(j).ne.0.0_fp) go to 10
  enddo
  !
  go to 1000                        ! all zero:  exit now...
  !
10 continue
  jnz=j
  do jj=1,jnz-1
    zprof(jj)=zprof(jnz)
  enddo
  !
  !  2.  fill in from right
  !
  do j=nx,1,-1
    if(zprof(j).ne.0.0_fp) go to 20
  enddo
  !
20 continue
  jnz=j
  do jj=nx,jnz+1,-1
    zprof(jj)=zprof(jnz)
  enddo
  !
  !  3.  fill in remaining gaps with linear interpolation
  !
  isrch=1
  do j=1+1,nx-1
    if(isrch.eq.1) then
      !  find start of zeroed section
      if(zprof(j).eq.0.0_fp) then
        jsave1=j-1
        zx1=xi(jsave1)
        zf1=zprof(jsave1)
        isrch=2
      endif
    else if(isrch.eq.2) then
      !  find end of zeroed section
      if(zprof(j).ne.0.0_fp) then
        jsave2=j
        zx2=xi(jsave2)
        zf2=zprof(jsave2)
        isrch=1
        !  apply patch
        do jj=jsave1+1,jsave2-1
          zx=xi(jj)
          zf=zf1+(zf2-zf1)*(zx-zx1)/(zx2-zx1)
          zprof(jj)=zf
          zprof(jj)=max(0.0_fp,zprof(jj))
        enddo
      endif
    endif                          ! isrch mode
  enddo                             ! j loop
  !
1000 continue
  return
end subroutine chk_interpolate
