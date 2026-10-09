subroutine r8bsmooh(zprof,znorm)
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp), dimension(*) :: zprof, znorm

  call r8bsmoop(zprof,znorm)   ! patch holes

  call r8bsmoo(zprof,znorm)    ! smooth result
 
  return
end subroutine r8bsmooh

!------------------------------------------------------------
subroutine r8bsmoop(zprof,znorm)
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
  !    end if
  !
  !  which can leave "holes" in the profile in case of lousy statistics.
  !
  !  *** SO ***
  !  fill these holes with flat extrapolation & linear interpolation
  !  before smoothing with "r8bsmoo"
  !
  !-----------------------------
  !
  use iso_c_binding, only: fp => c_double
  use r8bsmoo_mod
  implicit none
  !
  real(fp) :: zprof(MJ)
  real(fp) :: znorm(MJ)
  !
  real(fp) :: z1(mj),zx1,zf1,zx2,zf2,zx,zf
  !
  integer :: j,jnz,jj,isrch,jsave1,jsave2
  !
  !-----------------------------
  !
  !  0.  normalize
  !
  do j=lcentr,ledge
    zprof(j)=zprof(j)/znorm(j)
    z1(j)=1.0_fp
  end do
  !
  !  1.  fill in from left
  !
  do j=lcentr,ledge
    if(zprof(j).ne.0.0_fp) go to 10
  end do
  !
  go to 1000                        ! all zero:  exit now...
  !
10 continue
  jnz=j
  do jj=lcentr,jnz-1
    zprof(jj)=zprof(jnz)
  end do
  !
  !  2.  fill in from right
  !
  do j=ledge,lcentr,-1
    if(zprof(j).ne.0.0_fp) go to 20
  end do
  !
20 continue
  jnz=j
  do jj=ledge,jnz+1,-1
    zprof(jj)=zprof(jnz)
  end do
  !
  !  3.  fill in remaining gaps with linear interpolation
  !
  isrch=1
  do j=lcentr+1,ledge-1
    if(isrch.eq.1) then
      !  find start of zeroed section
      if(zprof(j).eq.0.0_fp) then
        jsave1=j-1
        zx1=xi(jsave1,2)
        zf1=zprof(jsave1)
        isrch=2
      end if
    else if(isrch.eq.2) then
      !  find end of zeroed section
      if(zprof(j).ne.0.0_fp) then
        jsave2=j
        zx2=xi(jsave2,2)
        zf2=zprof(jsave2)
        isrch=1
        !  apply patch
        do jj=jsave1+1,jsave2-1
          zx=xi(jj,2)
          zf=zf1+(zf2-zf1)*(zx-zx1)/(zx2-zx1)
          zprof(jj)=zf
        end do
      end if
    end if                          ! isrch mode
  end do                            ! j loop
  !
  !  un-normalize
  !
  do j=lcentr,ledge
    zprof(j)=zprof(j)*znorm(j)
  end do
  !
1000 continue
  return
end subroutine r8bsmoop
