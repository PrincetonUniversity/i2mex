!-----------------------------------------------------------------------
!  PLFMPA - Common Block
!
!
!  Mod tbt 13 Apr 1994 -- Added ZCosTabl and ZSinTabl to be able to
!                         calculate Sin and Cos functions once.
!  mod dmc 31 Mar 1994 -- support updown asymmetric moments sets
!  from TRANSP
!
!  COMMON BLOCK FOR PLF ROUTINES WITHIN RPLOT SUBROUTINE LIBRARY
!  USED FOR GENERATING MAP OF DATA ON 2D (R,Y) GRID, IN ROUTINES
!  CALLED FROM LIBRARY ROUTINE PLFMRY
!
!  COMMON BLOCK SUPPORTS COMPUTATION OF FAST MOMENTS INVERSE MAP TO
!  DO THIS PROBLEM.  THE PLASMA EQUILIBRIUM SURFACES ARE "EXPANDED"
!  TO EXTEND THE DEFINED SPATIAL RANGE OF THE MAP (CF SUBROUTINE
!  TKBCONRP ET AL)
!
!  MUST BE DECLARED AFTER 'CPLOTR' COMMON TO GET VALUE OF PARAMETER
!  NAXMOM
!
module plfmpa_mod
  use cplotr_mod, only: naxmom, naxmmp
  implicit none

  integer, parameter :: ICNTPS=128   ! NO. OF PTS IN SURFACE CONTOURS
  integer, parameter :: IGRIDL=1001  ! MAX NO. OF PTS IN R,Y GRIDS
  integer, parameter :: IMAXSF=401   ! MAX NO. OF SURFACES

  integer :: ICUR, &
       IMMGEO, &
       INDR0,ILR0,INDMM(2,NAXMOM),ILCMM(2,NAXMOM),IMOM, &
       INDAMM(4,0:NAXMOM),ILCAMM(4,0:NAXMOM), &
       INDRMP,ILRMP,INDYMP,ILYMP, &
       IT1,IT2,INX, &
       INUMR,INUMY, &
       INUMSF,ICRVSF, &
       nCosTabl, &
       INXP1

  real :: ZR0,ZY0, &
       ZR(ICNTPS),ZY(ICNTPS),ZCURV(ICNTPS),ZDIST(ICNTPS),ZANGS(2), &
       ZRM1(ICNTPS),ZYM1(ICNTPS), &
       ZRM2(ICNTPS),ZYM2(ICNTPS), &
       ZXTRAP,ZTFIX, &
       Z1,Z2, &
       ZRGRID(IGRIDL),ZYGRID(IGRIDL), &
       ZRMM0A(IMAXSF),ZRMOMA(NAXMOM,IMAXSF),ZYMOMA(NAXMOM,IMAXSF), &
       ZRMC(0:NAXMOM,2),ZYMC(0:NAXMOM,2), &
       ZRMCA(0:NAXMOM,2,IMAXSF),ZYMCA(0:NAXMOM,2,IMAXSF), &
       SNTHTK(NAXMOM),CSTHTK(NAXMOM), &
       ZRTARG,ZYTARG, &
       ZCosTabl(0:NaxMom,NaxMMP), ZSinTabl(0:NaxMom,NaxMMP), &
       ZFMOMSR(IMAXSF,4,2,0:NAXMOM), &
       ZFMOMSY(IMAXSF,4,2,0:NAXMOM),ZFMOMSX(IMAXSF)

!
!  ICUR -- INDEX TO CURRENT EXPANDED SURFACE
!
!  ZR0,ZY0  -- REFERENCE POINT IN "MIDDLE OF" REGION BOUNDED BY
!    CURRENT SURFACE
!
!  ZR(..),ZY(..)  THE SURFACE CONTOUR (EVALUATED MOMENTS EXPANSION)
!   ZCURV(..)     CURVATURE OF SURFACE AT EACH EVALUATED POINT
!
!   ZDIST(..)  METRIC DATA USED BY TKBCON I.E. TKBCONRP, RPLOT VERSION
!   ZANGS(..)  ANGLE DATA USED BY TKBCON
!
!  ZRM1(..),ZYM1(..) -- PRECEDING SURFACE CONTOUR (TKBCON)
!  ZRM2(..),ZYM2(..) -- PRE-PRECEDING SURFACE CONTOUR
!
!  ZFMOMSR/Y -- table of functions xi**M*[Mth moments] and spline coeffs
!  ZFMOMSX -- shifted XI axis for ZFMOMS
!
!-------------------
!
!  ILR0 -- FCN ID 0'TH R MOMENT
!   INDR0 -- R0 fcn number for memory management
!  ILCMM -- FCN ID'S HIGHER R AND Y MOMENTS
!   INDMM -- moments fcn numbers for memory management
!  IMOM -- NUMBER OF HIGHER MOMENTS
!
!-------------------
!
!  New March 1994 -- to support updown asymmetry
!
!  IMMGEO -- =0 for symmetric geometry, =1 for asymmetric geometry
!
!  if IMMGEO=1, then ILR0,INDR0,ILCMM,INDMM are not used!
!  if IMMGEO=0, then ILCAMM,INDAMM,ILRMP,ILYMP,INDRMP,INDYMP not used!
!
!  ILCAMM -- fcn IDs for all moments functions
!  INDAMM -- memory management fcn no.s for all moments functions
!    1st index:  1 for R cos moments
!                2 for R sin moments
!                3 for Y cos moments
!                4 for Y sin moments
!
!  ILRMP,INDRMP -- ditto for R(midplane)
!  ILYMP,INDYMP -- ditto for Y(midplane)
!
!-------------------
!
!  IT1,IT2,Z1,Z2 -- TIME INTERPOLATION DATA
!
!    INTERPOLATION = Z1*DATBUF(IADR(IT1))+Z2*DATBUF(IADR(IT2))
!      WHERE IADR(I) YIELDS THE DESIRED DATA ADDRESS FOR EACH TIME
!      INDEX
!
!  INX = NUMBER OF SURFACES
!  INXP1 = INX+1
!
!-------------------
!
!  ZRGRID,ZYGRID  -- R,Y GRID TO MAP
!   INUMR,INUMY --   NUMBERS IN USE
!
!  ZRTARG,ZYTARG  -- CURRENT TARGET FOR INVERSE MAP
!
!-------------------
! for symmetric geometry only:
!  ZRMM0A  -- R0 MOMENT FOR EACH SURFACE INCLUDING EXTRAPOLATED SURFACES
!  ZRMOMA  -- R HIGHER MOMENTS FOR EACH SURFACE
!  ZYMOMA  -- Y HIGHER MOMENTS FOR EACH SURFACE
! for updown asymmetric geometry:
!  ZRMC    -- R moments
!  ZYMC    -- Y moments
! for all geometries:
!  INUMSF  -- NUMBER OF SURFACES (INCLUDING EXTRAPOLATED SURFACES)
!  ICRVSF  -- INDEX TO FIRST SURFACE WHICH IS CONCAVE EVERYWHERE
!
!  SNTHTK,CSTHTK -- SIN/COS LADDER FOR EVALUATING MOMENTS
!
!------------------
! April 1994, TBT
!
!        ZCosTabl(0:NaxMom,NaxMMP), ZSinTabl(0:NanMom,NaxMMP)
!  Cos and Sin tables so terms are not calculated every time.
!   NCosTabl -- size of precomputed ZCosTabl,ZSinTabl
!
!-----------------------------------------------------------------------

end module plfmpa_mod
