C-------------------------------------------------------------
C  TKBCON SUBROUTINE - DMC - VERSION MODIFIED FOR RPLOT USE
C   DECEMBER 1988
C
C   TKBMMC.FOR - TRANSP VERSION
C   TKBMMCRP.FOR - RPLOT VERSION
C     DIFFERENT COMMONS
C
C   GIVEN (R,Y) CLOSED CONTOUR USE FFT TO GET FIT MOMENTS
C   THE NUMBER OF MOMENTS IS ASSUMED NOT TO INCREASE; HIGH ORDER
C   MOMENTS WILL GET DROPPED BUT SHOULD NOT BE NEEDED
C
      SUBROUTINE TKBMMCRP(ISURF)
C
C  FROM R,Y CONTOUR DATA DEFINE THE MOMENTS ARRAYS FOR 1 SURFACE
C
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
      REAL ZS1D(ICNTPS)
      REAL ZC1D(ICNTPS)
C
C------------------------------------------------------------------
C
      ITHETA=ICNTPS
C
C  FFT
C
      call r4fftsc(zr, itheta, zs1d, zc1d, ier)
      if(ier.ne.0) call errmsg_exit(' ?tkbmmcrp: FFT error (1)!')
C
C  NORMALIZE
C
      if(immgeo.eq.0) then
        ZRMM0A(ISURF) = ZC1D(1)/(2.0*ITHETA)
C
        DO 20300 J = 2, IMOM+1
C
          ZRMOMA(J-1,ISURF) = ZC1D(J)/FLOAT(ITHETA)
C
20300   CONTINUE
C
      else
        zrmca(0,1,isurf)=zc1d(1)/(2.0*itheta)
        zrmca(0,2,isurf)=0.0
C
        do j = 2, imom+1
           jm1=j-1
           zrmca(jm1,1,isurf)=zc1d(j)/float(itheta)
           zrmca(jm1,2,isurf)=zs1d(j)/float(itheta)
        enddo
C
      endif
C
C  FFT
C
      call r4fftsc(zy, itheta, zs1d, zc1d, ier)
      if(ier.ne.0) call errmsg_exit(' ?tkbmmcrp: FFT error (1)!')
C
C  NORMALIZE
C
C
      if(immgeo.eq.0) then
        DO 20400 J = 2, IMOM+1
C
          ZYMOMA(J-1,ISURF) = ZS1D(J)/FLOAT(ITHETA)
C
20400   CONTINUE
C
      else
        zymca(0,1,isurf)=zc1d(1)/(2.0*itheta)
        do j = 2, imom+1
           jm1=j-1
           zymca(jm1,1,isurf)=zc1d(j)/float(itheta)
           zymca(jm1,2,isurf)=zs1d(j)/float(itheta)
        enddo
      endif
C
      RETURN
      END
C******************** END FILE TKBMMC.FOR ; GROUP TKBLOAT **************
