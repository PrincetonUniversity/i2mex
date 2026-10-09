      subroutine getmpa(zrmc,zymc,mimom,imom,zrin,zrout,zymp)
use iso_c_binding, only: fp => c_double
!
!  DMC 18 March 1994:
!  given the updown asymmetric moments of a surface, define an
!  approximate midplane to the surface.  Return the inner and
!  outer R intercepts of this midplane with the surface, and
!  return the height of the surface above y=0.
!
!  rev dmc 25 Mar 1994:  use
!
!============
      implicit none
      integer :: imom,mimom,inth,imm,ith
!============
      REAL zth,zr,zy
!============
      REAL zrmc(0:mimom,2)  ! R moments
      REAL zymc(0:mimom,2)  ! Y moments
!
!  mimom gives the moments 1st array dimension; imom gives the number
!   of *active* moments; imom.le.mimom ...
!
      REAL zrin            ! inner intercept R (output)
      REAL zrout           ! outer intercept R (output)
      REAL zymp            ! midplane vertical displacement (output)
!
      parameter (inth=400)
      REAL zrcon(inth),zycon(inth)
      REAL zcos(0:imom),zsin(0:imom)
!---------------------------------------------------------------------
!
      zcos(0)=1
      zsin(0)=0
      do ith=1,inth-1
        zth=(6.2831853071795862*(ith-1))/(inth-1)
        zr=0.0
        zy=0.0
        call sincos(zth,imom,zsin(1:imom),zcos(1:imom))
        do imm=0,imom
          zr=zr+zrmc(imm,1)*zcos(imm)+zrmc(imm,2)*zsin(imm)
          zy=zy+zymc(imm,1)*zcos(imm)+zymc(imm,2)*zsin(imm)
        end do
        zrcon(ith)=zr
        zycon(ith)=zy
      end do
!
      zrcon(inth)=zrcon(1)
      zycon(inth)=zycon(1)
!
      call getmpa_ry(zrcon,zycon,inth,zrin,zrout,zymp)
!
      return
      end
!---------------------------------------------------------------------
!
      subroutine getmpa_ry(zrcon,zycon,inth,zrin,zrout,zymp)
use iso_c_binding, only: fp => c_double
!
!  find midplane elevation and R intercepts based on centroid, starting
!  from a closed contour of (R,Y) pairs
!
      implicit none

      integer, intent(in) :: inth
      real, intent(in) :: zrcon(inth),zycon(inth)  ! closed contour

      real, intent(out) :: zrin,zrout    ! midplane intercepts (approx.)
      real, intent(out) :: zymp          ! midplane elevation (centroid)
!
!  local:
!
      real :: zrmp,zr,zy,zrp,zyp,zytest,zrans
      integer :: ith
!
!  find centroid
!
      call plcentr(zrcon,zycon,inth,zrmp,zymp)
!
!  find approximate midplane intercept locations
!
      zr=zrcon(1)
      zy=zycon(1)
!
      do ith=2,inth
!
        zrp=zr
        zyp=zy
        zr=zrcon(ith)
        zy=zycon(ith)
!
        zytest=(zy-zymp)*(zymp-zyp)
        if(zytest.ge.0.0) then
          zrans=(zr*(zymp-zyp)+zrp*(zy-zymp))/(zy-zyp)
          if(zrans.gt.zrmp) then
            zrout=zrans
          else
            zrin=zrans
          end if
        end if
      end do
!
      return
      end
