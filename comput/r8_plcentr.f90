subroutine r8_plcentr(R,Z,N,RCENTR,ZCENTR)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, half
  implicit none

  ! 
  !     USE GREEN'S THEOREM TO CALCULATE CENTROIDS
  !     DICK WIELAND
  !
  !     RGA: Jun2010, switch to centroid of inscribed polygon
  !
  real(fp), dimension(*) :: R, Z
  real(fp) :: RCENTR, ZCENTR
  integer :: N

  real(fp) :: AREA,LINT,MR,MZ
  integer :: I

  !     CALCULATE AREA ENCLOSED
  AREA = zero
  MR   = zero
  MZ   = zero
  do I = 1,N-1
    LINT = R(I)*Z(I+1)-R(I+1)*Z(I)
    AREA = AREA + LINT
    MR   = MR + (R(I)+R(I+1))*LINT
    MZ   = MZ + (Z(I)+Z(I+1))*LINT
  end do
  AREA = half*AREA

  if (abs(AREA)>1.e-10_fp) then
    MR = MR/(6.0_fp*AREA)
    MZ = MZ/(6.0_fp*AREA)
  else
    MR=zero
    MZ=zero
    do i = 1, N
      MR = MR+R(I)
      MZ = MZ+Z(I)
    end do
    MR = MR/max(1,N)
    MZ = MZ/max(1,N)
  end if

  RCENTR = MR
  ZCENTR = MZ

  return
end subroutine r8_plcentr

