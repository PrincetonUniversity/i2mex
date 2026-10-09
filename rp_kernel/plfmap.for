C-----------------------------------------------------------------------
C  PLFMAP  GET MOMENTS AT CURRENT TIME; GET BOUNDARY SURFACE
C
      SUBROUTINE PLFMAP
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
C-----------------------------------------------------------------------
C
C  geometry generalization -- dmc 16 Jun 1994
C
      if(immgeo.eq.1) go to 100
C
      DO 50 IX=1,INX
C
        I1=NDPTR(ILR0,INX,IT1)+IX-1
        I2=NDPTR(ILR0,INX,IT2)+IX-1
        ZRMM0A(IX)=DATBUF(I1)*Z1+DATBUF(I2)*Z2
C
        DO 70 IM=1,IMOM
C  R MOMENTS
          I1=NDPTR(ILCMM(1,IM),INX,IT1)+IX-1
          I2=NDPTR(ILCMM(1,IM),INX,IT2)+IX-1
          ZRMOMA(IM,IX)=DATBUF(I1)*Z1+DATBUF(I2)*Z2
C  Y MOMENTS
          I1=NDPTR(ILCMM(2,IM),INX,IT1)+IX-1
          I2=NDPTR(ILCMM(2,IM),INX,IT2)+IX-1
          ZYMOMA(IM,IX)=DATBUF(I1)*Z1+DATBUF(I2)*Z2
 70     CONTINUE
C
 50   CONTINUE
C
C  DEFINE THE CONTOUR
      CALL PLMMCV(ZRMM0A(INX),ZRMOMA(1,INX),ZYMOMA(1,INX),IMOM,
     >    1,ZR,ZY,ICNTPS,1)
C
      go to 500
C-----------------------------------------------------------
C  asymmetric geometry
C
 100  continue
C
      do ix=1,inx
         do im=0,imom
            i1=ndptr(ilcamm(1,im),inx,it1)+ix-1
            i2=ndptr(ilcamm(1,im),inx,it2)+ix-1
            zrmca(im,1,ix)=datbuf(i1)*z1+datbuf(i2)*z2  ! R Cos
            i1=ndptr(ilcamm(2,im),inx,it1)+ix-1
            i2=ndptr(ilcamm(2,im),inx,it2)+ix-1
            zrmca(im,2,ix)=datbuf(i1)*z1+datbuf(i2)*z2  ! R Sin
            i1=ndptr(ilcamm(3,im),inx,it1)+ix-1
            i2=ndptr(ilcamm(3,im),inx,it2)+ix-1
            zymca(im,1,ix)=datbuf(i1)*z1+datbuf(i2)*z2  ! Y Cos
            i1=ndptr(ilcamm(4,im),inx,it1)+ix-1
            i2=ndptr(ilcamm(4,im),inx,it2)+ix-1
            zymca(im,2,ix)=datbuf(i1)*z1+datbuf(i2)*z2  ! Y Sin
         enddo
      enddo
C
      call plammcv(
     1   zrmca(0,1,inx),zrmca(0,2,inx),zymca(0,1,inx),zymca(0,2,inx),
     2   imom,1,zr,zy,icntps,1)
C
C-----------------------------------------------------------
 500  continue
      RETURN
      END
