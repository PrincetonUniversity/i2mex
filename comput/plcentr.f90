      subroutine PLCENTR (R,Z,N,RCENTR,ZCENTR)
use iso_c_binding, only: fp => c_double

      implicit none

!     USE GREEN'S THEOREM TO CALCULATE CENTROIDS
!     DICK WIELAND
!
!     RGA: Jun2010, switch to centroid of inscribed polygon
!
      REAL R(*),Z(*),RCENTR,ZCENTR
      integer :: N

      REAL AREA,LINT,MR,MZ
      integer :: I

!     CALCULATE AREA ENCLOSED
      AREA = 0.
      MR   = 0.
      MZ   = 0.
      do I = 1,N-1
          LINT = R(I)*Z(I+1)-R(I+1)*Z(I)
          AREA = AREA + LINT
          MR   = MR + (R(I)+R(I+1))*LINT
          MZ   = MZ + (Z(I)+Z(I+1))*LINT
      end do
      AREA = AREA/2.

      if (abs(AREA)>1.e-10) then
         MR = MR/(6.*AREA)
         MZ = MZ/(6.*AREA)
      else
         MR=0.
         MZ=0.
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
      end

