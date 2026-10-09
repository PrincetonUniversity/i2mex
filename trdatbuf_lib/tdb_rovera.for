C----------------------------------------------
C  COMPUTE R/A VECTOR FOR DATA INTERPOLATION
C
      SUBROUTINE TDB_ROVERA(t,INRI,iarat,ZX,drshaf,izones)
      use trdatbuf_aux
      IMPLICIT NONE
      type (profget) :: t
C
C  dmc 4 May 1993 -- for INRI <= 4 this is a midplane spatial coordinate
C  dmc June 2005 -- if IARAT <> 0 this is a normalized aspect ratio
C
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER inri,iz,j,icen,iarat
!============
      integer :: izones
      REAL*8 ZX(izones),drshaf(izones)
      real*8 :: zr1,zr2,zr1p,zr2p,za1,za2,zra
C
      icen=izones
      if(inri.ge.5) then
C
C  flux coordinate
C
         drshaf = 0.0_R8
         if(t%ibdy) then
            zx=t%xibdys
         else
            do j=1,izones-1
               zx(j)=0.5_R8*(t%xibdys(j)+t%xibdys(j+1))
            enddo
            zx(izones)=t%xibdys(izones)+0.5_R8*
     >           (t%xibdys(izones)-t%xibdys(izones-1))
         endif
      else if(iarat.ge.0) then
C
C  minor or major radius coordinate
C
         za2=t%rmajmp(t%nrmaj)-t%rmajmp(1)
         zra=0.5_R8*(t%rmajmp(t%nrmaj)+t%rmajmp(1))
         if(t%ibdy) then
            zx(1)=0.0_R8
            drshaf(1)=t%rmajmp(icen)-zra
            do j=1,izones-1
               zr1=t%rmajmp(icen-j)
               zr2=t%rmajmp(icen+j)
               zx(j+1)=(zr2-zr1)/za2
               drshaf(j+1)=0.5_R8*(zr2+zr1)-zra
            enddo
         else
            zr1=t%rmajmp(icen)
            zr2=zr1
            do j=1,izones-1
               zr1p=zr1
               zr2p=zr2
               zr1=t%rmajmp(icen-j)
               zr2=t%rmajmp(icen+j)
               zx(j)=0.5_R8*(zr2+zr2p-(zr1+zr1p))/za2
               drshaf(j)=0.25_R8*(zr2+zr2p+zr1+zr1p) - zra
            enddo
            zx(izones)=2.0_r8-zx(izones-1)
            drshaf(izones) = 0.0_R8
         endif
      else
C
C  normalized aspect ratio -- used for final mapping of presymmetrized data
C
         drshaf = 0.0_R8
C
C  bdy aspect ratio for normalization:
C
         za2 = (t%rmajmp(t%nrmaj)-t%rmajmp(1))/
     >        (t%rmajmp(t%nrmaj)+t%rmajmp(1))
C
         if(t%ibdy) then
            zx(1)=0.0_R8
            do j=1,izones-1
               if(iarat.gt.0) then
                  zr1=t%rmajmp(icen-j)
                  zr2=t%rmajmp(icen+j)
               else
                  zr1=1.0_R8/t%bmidp(icen-j)
                  zr2=1.0_R8/t%bmidp(icen+j)
               endif
               za1=(zr2-zr1)/(zr2+zr1)
               zx(j+1)=za1/za2
            enddo
         else
            if(iarat.gt.0) then
               zr1=t%rmajmp(icen)
            else
               zr1=1.0_R8/t%bmidp(icen)
            endif
            zr2=zr1
            do j=1,izones-1
               zr1p=zr1
               zr2p=zr2
               if(iarat.gt.0) then
                  zr1=t%rmajmp(icen-j)
                  zr2=t%rmajmp(icen+j)
               else
                  zr1=1.0_R8/t%bmidp(icen-j)
                  zr2=1.0_R8/t%bmidp(icen+j)
               endif
               ! avg aspect ratio half way btw two surfaces
               za1=((zr2+zr2p)-(zr1+zr1p))/((zr2+zr2p)+(zr1+zr1p))
               zx(j)=za1/za2
            enddo
            zx(izones)=2.0_r8-zx(izones-1) ! 1/2 zone width beyond...
         endif
      endif
      RETURN
      END
