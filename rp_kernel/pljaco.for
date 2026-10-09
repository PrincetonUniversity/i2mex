C=====================================================================
C  PLJACO  compute (R,Y) and the 2 x 2 Jacobian  d(R,Y)/d(xi,theta)
C          given (xi,theta) as input
C
C  use the bloated surfaces set up in the PLFMPA COMMON.  This gives
C  a geometry that extends well beyond the plasma.
C
C  this is similar to the TRANSP MOMRY subroutine.
C
C  mod dmc Sept 1995 -- use splines for zxiin < 1
C
      subroutine PLJACO(zxiin,zthin,ijac,zrout,zyout,zjaco,istat)

      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod

      integer ijac,istat
      real zxii,zthi,zrout,zyout,zjaco(2,2)
C
C input
C  ZXIIN  --  XI coordinate (sqrt(trflux)/sqrt(trflux-lim))
C  ZTHIN  --  THETA coordinate (going poloidally around the surface
C  IJAC   --  1 if the Jacobian is wanted, 0 otherwise
C
C output
C  ZROUT --  R(xi,theta)
C  ZYOUT --  Y(xi,theta)
C  ZJACO(2,2) -- d(R,Y)/d(xi,theta)
C    ZJACO(1,1) = dR/dxi         ZJACO(1,2) = dR/dtheta
C    ZJACO(2,1) = dY/dxi         ZJACO(2,2) = dY/dtheta
C
C  ISTAT -- =0 -- inside plasma, ZXI.le.1.0
C           =1 -- outside plasma but inside bloated surfaces
C           =2 -- beyond the bloated surfaces (no Jacobian computed)
C
C--------------------------------------------------------------------
C
C  for the Jacobian:
C
C  numerical d/dxi of moments
C
        real zdrmdxi(2,0:NAXMOM)  ! (cos:sin,moment index)
        real zdymdxi(2,0:NAXMOM)
C
C  the interpolated moments themselves
C
        real zrmma(2,0:NAXMOM)
        real zymma(2,0:NAXMOM)
C
        logical jactest
C
      data ztwopi/6.2831853072/
C
C--------------------------------------------------------------------
C
        ZXI=ZXIIN
        ZTH=ZTHIN
C
      IF(ZXI.LT.0.0) THEN
        ZXI=-ZXI
        ZTH=ZTH+3.141593
      ENDIF
C
      IF(ZTH.GE.ZTWOPI) THEN
         ZTH=ZTH-ZTWOPI
      ELSE IF(ZTH.LT.0.0) THEN
         ZTH=ZTH+ZTWOPI
      ENDIF
C
      ZDXI=1.0/FLOAT(INX)
      IND=ZXI/ZDXI
      INDP=IND+1
C
      IF(INDP.GT.INUMSF) THEN
        IND=INUMSF
        INDP=INUMSF
        ZI=0.0
        ZIP=1.0
        ZFACX=ZXI/(INUMSF*ZDXI)
          ISTAT=2
      ELSE IF(IND.EQ.0) THEN
        IND=1
        INDP=1
        ZI=0.0
        ZIP=1.0
        ZFACX=1.0
          ISTAT=0
      ELSE
        ZXI1=IND*ZDXI
        ZIP=(ZXI-ZXI1)/ZDXI
        ZI=1.0-ZIP
        ZFACX=1.0
          ISTAT=0
          if(ZXI.ge.1.0) ISTAT=1
      ENDIF
C
        jactest = (ijac.eq.1).and.(istat.lt.2)
C
      ZXII=IND*ZDXI
      ZXIP=INDP*ZDXI
C
      CALL SINCOS(ZTH,IMOM,SNTHTK,CSTHTK)
C
      IF(IMMGEO.EQ.0) THEN
C
C  0th moment, symmetric case
C
         ZY0SUM=0.0
           if(istat.eq.0) then
              call pljaco_Rspl(zxi,1.0e0,0.0e0,1,0,ZR0SUM,ZR0DOX)
           else
              ZR0SUM=ZI*ZRMM0A(IND)+ZIP*ZRMM0A(INDP)
           endif
           if(jactest) then
              zrmma(2,0)=0.0
              zymma(2,0)=0.0
              zymma(1,0)=0.0
              zrmma(1,0)=ZR0SUM
              zdrmdxi(2,0)=0.0
              zdymdxi(2,0)=0.0
              zdymdxi(1,0)=0.0
              if(istat.eq.0) then
                 zdrmdxi(1,0)=ZR0DOX
              else
                 if(indp.eq.1) then
                    zdrmdxi(1,0)=0.0
                 else
                    zdrmdxi(1,0)=(zrmm0a(indp)-zrmm0a(ind))/zdxi
                 endif
              endif
           endif
      ELSE
C
C   0th moment, asymmetric case
C
           if(istat.eq.0) then
              call pljaco_Rspl(zxi,1.0e0,0.0e0,1,0,ZR0SUM,ZR0DOX)
              call pljaco_Yspl(zxi,1.0e0,0.0e0,1,0,ZY0SUM,ZY0DOX)
           else
              ZR0SUM=ZI*ZRMCA(0,1,IND)+ZIP*ZRMCA(0,1,INDP)
              ZY0SUM=ZI*ZYMCA(0,1,IND)+ZIP*ZYMCA(0,1,INDP)
           endif
           if(jactest) then
              zrmma(2,0)=0.0
              zymma(2,0)=0.0
              zymma(1,0)=ZY0SUM
              zrmma(1,0)=ZR0SUM
              zdrmdxi(2,0)=0.0
              zdymdxi(2,0)=0.0
              if(istat.eq.0) then
                 zdrmdxi(1,0)=ZR0DOX
                 zdymdxi(1,0)=ZY0DOX
              else
                 if(indp.eq.1) then
                    zdrmdxi(1,0)=0.0
                    zdymdxi(1,0)=0.0
                 else
                    zdrmdxi(1,0)=(zrmca(0,1,indp)-zrmca(0,1,ind))/zdxi
                    zdymdxi(1,0)=(zymca(0,1,indp)-zymca(0,1,ind))/zdxi
                 endif
              endif
           endif
      ENDIF
C
C   higher moments
C
      DO IM=1,IMOM
         ZXIIX=ZXII**IM        ! bdy values
         ZXIPX=ZXIP**IM        ! bdy values
           if(IM.EQ.1) then
              ZXIXM1=1.0
           else
              ZXIXM1=ZXI**(IM-1)  ! xi**(IM-1)
           endif
         ZXIX=ZXIXM1*ZXI        ! xi**(IM)
C
           if(istat.eq.0) then
              call pljaco_rspl(zxi,zxix,zxixm1,
     >             1,im,zrmma(1,im),zdrmdxi(1,im))
              call pljaco_rspl(zxi,zxix,zxixm1,
     >             2,im,zrmma(2,im),zdrmdxi(2,im))
              call pljaco_yspl(zxi,zxix,zxixm1,
     >             1,im,zymma(1,im),zdymdxi(1,im))
              call pljaco_yspl(zxi,zxix,zxixm1,
     >             2,im,zymma(2,im),zdymdxi(2,im))
              ZR0SUM=ZR0SUM+
     >          (zrmma(1,im)*CSTHTK(IM)+zrmma(2,im)*SNTHTK(IM))
              ZY0SUM=ZY0SUM+
     >          (zymma(1,im)*CSTHTK(IM)+zymma(2,im)*SNTHTK(IM))
         else IF(IMMGEO.EQ.0) THEN
C  symmetric case
           IF(INDP.GT.1) THEN
C  zones bound away from axis
      	ZRPAR=
     >             (ZI*ZRMOMA(IM,IND)/ZXIIX+ZIP*ZRMOMA(IM,INDP)/ZXIPX)
                ZRMM=ZXIX*ZRPAR
      	ZYPAR=
     >             (ZI*ZYMOMA(IM,IND)/ZXIIX+ZIP*ZYMOMA(IM,INDP)/ZXIPX)
                ZYMM=ZXIX*ZYPAR
                if(jactest) then
                   zrmma(2,IM)=0.0
                   zymma(1,IM)=0.0
                   zrmma(1,IM)=ZRMM
                   zymma(2,IM)=ZYMM
                   zdrmdxi(2,IM)=0.0
                   zdymdxi(1,IM)=0.0
                   zdrmdxi(1,IM)=IM*ZXIXM1*ZRPAR +
     >                 ZXIX/ZDXI*
     >                   (ZRMOMA(IM,INDP)/ZXIPX-ZRMOMA(IM,IND)/ZXIIX)
                   zdymdxi(2,IM)=IM*ZXIXM1*ZYPAR +
     >                 ZXIX/ZDXI*
     >                   (ZYMOMA(IM,INDP)/ZXIPX-ZYMOMA(IM,IND)/ZXIIX)
                endif
           ELSE
C  innermost zone
      	ZFAC0=ZXIX/ZXIPX
      	ZRMM=ZRMOMA(IM,INDP)*ZFAC0
      	ZYMM=ZYMOMA(IM,INDP)*ZFAC0
                if(jactest) then
                   zrmma(2,IM)=0.0
                   zymma(1,IM)=0.0
                   zrmma(1,IM)=ZRMM
                   zymma(2,IM)=ZYMM
                   zdrmdxi(2,IM)=0.0
                   zdymdxi(1,IM)=0.0
                   zdrmdxi(1,IM)=IM*ZXIXM1*ZRMOMA(IM,INDP)/ZXIPX
                   zdymdxi(2,IM)=IM*ZXIXM1*ZYMOMA(IM,INDP)/ZXIPX
                endif
           ENDIF
           ZR0SUM=ZR0SUM+ZFACX*ZRMM*CSTHTK(IM)
           ZY0SUM=ZY0SUM+ZFACX*ZYMM*SNTHTK(IM)
         ELSE
C  Asymmetric case
           IF(INDP.GT.1) THEN
C  Zones bound away from axis
      	ZRMCPAR=
     >             (ZI*ZRMCA(IM,1,IND)/ZXIIX+ZIP*ZRMCA(IM,1,INDP)/ZXIPX)
                ZRMCL=ZXIX*ZRMCPAR
      	ZRMSPAR=
     >             (ZI*ZRMCA(IM,2,IND)/ZXIIX+ZIP*ZRMCA(IM,2,INDP)/ZXIPX)
                ZRMSL=ZXIX*ZRMSPAR
      	ZYMCPAR=
     >             (ZI*ZYMCA(IM,1,IND)/ZXIIX+ZIP*ZYMCA(IM,1,INDP)/ZXIPX)
                ZYMCL=ZXIX*ZYMCPAR
      	ZYMSPAR=
     >             (ZI*ZYMCA(IM,2,IND)/ZXIIX+ZIP*ZYMCA(IM,2,INDP)/ZXIPX)
                ZYMSL=ZXIX*ZYMSPAR
                if(jactest) then
                   zrmma(2,IM)=ZRMSL
                   zymma(1,IM)=ZYMCL
                   zrmma(1,IM)=ZRMCL
                   zymma(2,IM)=ZYMSL
                   zdrmdxi(1,IM)=IM*ZXIXM1*ZRMCPAR +
     >                 ZXIX/ZDXI*
     >                   (ZRMCA(IM,1,INDP)/ZXIPX-ZRMCA(IM,1,IND)/ZXIIX)
                   zdrmdxi(2,IM)=IM*ZXIXM1*ZRMSPAR +
     >                 ZXIX/ZDXI*
     >                   (ZRMCA(IM,2,INDP)/ZXIPX-ZRMCA(IM,2,IND)/ZXIIX)
                   zdymdxi(1,IM)=IM*ZXIXM1*ZYMCPAR +
     >                 ZXIX/ZDXI*
     >                   (ZYMCA(IM,1,INDP)/ZXIPX-ZYMCA(IM,1,IND)/ZXIIX)
                   zdymdxi(2,IM)=IM*ZXIXM1*ZYMSPAR +
     >                 ZXIX/ZDXI*
     >                   (ZYMCA(IM,2,INDP)/ZXIPX-ZYMCA(IM,2,IND)/ZXIIX)
                endif
           ELSE
C  innermost zone
      	ZFAC0=ZXIX/ZXIPX
      	ZRMCL=ZRMCA(IM,1,INDP)*ZFAC0
      	ZRMSL=ZRMCA(IM,2,INDP)*ZFAC0
      	ZYMCL=ZYMCA(IM,1,INDP)*ZFAC0
      	ZYMSL=ZYMCA(IM,2,INDP)*ZFAC0
                if(jactest) then
                   zrmma(2,IM)=ZRMSL
                   zymma(1,IM)=ZYMCL
                   zrmma(1,IM)=ZRMCL
                   zymma(2,IM)=ZYMSL
                   zdrmdxi(1,IM)=IM*ZXIXM1*ZRMCA(IM,1,INDP)/ZXIPX
                   zdrmdxi(2,IM)=IM*ZXIXM1*ZRMCA(IM,2,INDP)/ZXIPX
                   zdymdxi(1,IM)=IM*ZXIXM1*ZYMCA(IM,1,INDP)/ZXIPX
                   zdymdxi(2,IM)=IM*ZXIXM1*ZYMCA(IM,2,INDP)/ZXIPX
                endif
           ENDIF
           ZR0SUM=ZR0SUM+ZFACX*(ZRMCL*CSTHTK(IM)+ZRMSL*SNTHTK(IM))
           ZY0SUM=ZY0SUM+ZFACX*(ZYMCL*CSTHTK(IM)+ZYMSL*SNTHTK(IM))
         ENDIF
        enddo
C
C  the (R,Y) position, output
C
        ZROUT=ZR0SUM
        ZYOUT=ZY0SUM
C
C  clear the Jacobian matrix
C
        ZJACO(1,1)=0.0
        ZJACO(2,1)=0.0
        ZJACO(1,2)=0.0
        ZJACO(2,2)=0.0
C
        if(.not.jactest) return
C
C  Evaluate the Jacobian:
C    for the d/dxi part sum the numerical d/dxi's evaluated above
C    use analytic formula for d/dtheta part
C
        ZJACO(1,1)=zdrmdxi(1,0)
        ZJACO(2,1)=zdymdxi(1,0)
C
        do im=1,imom
           ZJACO(1,1)=ZJACO(1,1) +
     >         (zdrmdxi(1,im)*csthtk(im)+zdrmdxi(2,im)*snthtk(im))
           ZJACO(2,1)=ZJACO(2,1) +
     >         (zdymdxi(1,im)*csthtk(im)+zdymdxi(2,im)*snthtk(im))
C
           ZJACO(1,2)=ZJACO(1,2) +
     >         im*(-zrmma(1,im)*snthtk(im)+zrmma(2,im)*csthtk(im))
           ZJACO(2,2)=ZJACO(2,2) +
     >         im*(-zymma(1,im)*snthtk(im)+zymma(2,im)*csthtk(im))
C
        enddo
C
        RETURN
        END
c-----------------------------------
      subroutine pljaco_rspl(zxi,zxim,zxixm1,ics,im,zm,zdmdx)
c
c  spline-interpolation of a moment:
c   a) spline interpolate the function x**-m*[m moment]
c   b) multiply in x**m
c
c  input:  zxi  -- x
c          zxim -- x**m  (evaluated once in caller for efficiency)
c          zxixm1 -- x**(m-1) (also evaluated by caller)
c          ics  -- sin or cos moment
c          im   -- moment number
c
c output:  zm   -- the moment value
c          zdmdx-- d/dx of the moment
c-----------------------------------------
C
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
c-----------------------------------------
c  use SEVALI COMMON
C
C  DMC GARCHING OCT 1985 KLUGE
      COMMON/ZSEVALI/ IPUT,ILIN
c
c----------------------------------------
c
        if(im.gt.0) then
           if((immgeo.eq.0).and.(ics.eq.2)) then
              zm=0.0
              zdmdx=0.0
              return
           endif
        endif
c
        ilin=0
c
        zmxm=sevali(inxp1,zxi,zfmomsx,zfmomsR(1,1,ics,im),
     >    zfmomsR(1,2,ics,im),zfmomsR(1,3,ics,im),zfmomsR(1,4,ics,im),
     >     dx)
c
        zm=zxim*zmxm
c
c  the derivative
c
        zdm=(3.0*zfmomsR(iput,4,ics,im)*dx +
     >       2.0*zfmomsR(iput,3,ics,im))*dx +
     >       zfmomsR(iput,2,ics,im)
c
        zdmdx = im*zxixm1*zmxm + zxim*zdm
c
        return
        end
c-----------------------------------
      subroutine pljaco_yspl(zxi,zxim,zxixm1,ics,im,zm,zdmdx)
c
c  spline-interpolation of a moment:
c   a) spline interpolate the function x**-m*[m moment]
c   b) multiply in x**m
c
c  input:  zxi  -- x
c          zxim -- x**m  (evaluated once in caller for efficiency)
c          zxixm1 -- x**(m-1) (also evaluated by caller)
c          ics  -- sin or cos moment
c          im   -- moment number
c
c output:  zm   -- the moment value
c          zdmdx-- d/dx of the moment
c-----------------------------------------
C
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
c-----------------------------------------
c  use SEVALI COMMON
C
C  DMC GARCHING OCT 1985 KLUGE
      COMMON/ZSEVALI/ IPUT,ILIN
c
c----------------------------------------
c
        if(im.gt.0) then
           if((immgeo.eq.0).and.(ics.eq.1)) then
              zm=0.0
              zdmdx=0.0
              return
           endif
        endif
c
        ilin=0
c
        zmxm=sevali(inxp1,zxi,zfmomsx,zfmomsY(1,1,ics,im),
     >    zfmomsY(1,2,ics,im),zfmomsY(1,3,ics,im),zfmomsY(1,4,ics,im),
     >    dx)
c
        zm=zxim*zmxm
c
c  the derivative
c
        zdm=(3.0*zfmomsY(iput,4,ics,im)*dx +
     >       2.0*zfmomsY(iput,3,ics,im))*dx +
     >       zfmomsY(iput,2,ics,im)
c
        zdmdx = im*zxixm1*zmxm + zxim*zdm
c
        return
        end
