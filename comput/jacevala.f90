      subroutine jacevala(DRMDROC,DYMDROC,DRMDROS,DYMDROS, &
                          RMSPLC,YMSPLC,RMSPLS,YMSPLS, &
                          ZSNTHTK,ZCSTHTK,NMOM,ZJAC)
use iso_c_binding, only: fp => c_double
!
!  evaluate 2x2 jacobian using moments and derivatives
!    ** asymmetric formula **
!
!  evaluate moments position:  asymmetric eq. moments formula
!
      real rmsplc(0:nmom),ymsplc(0:nmom),rmspls(0:nmom),ymspls(0:nmom)
      real drmdroc(0:nmom),dymdroc(0:nmom)
      real drmdros(0:nmom),dymdros(0:nmom)
      real zcsthtk(nmom),zsnthtk(nmom)
!
      real zjac(2,2)
!
!--------------------------------
!
      ZJAC(1,1)= DRMDROC(0)
      ZJAC(1,2)= 0.0
!
      ZJAC(2,1)= DYMDROC(0)
      ZJAC(2,2)= 0.0
!
      do 151 IMUL=1,NMOM
!
         ZC = ZCSTHTK(IMUL)
         ZS = ZSNTHTK(IMUL)
!
         ZJAC(1,1) = ZJAC(1,1) + DRMDROC(IMUL)*ZC &
            + DRMDROS(IMUL)*ZS
         ZJAC(1,2) = ZJAC(1,2) - RMSPLC(IMUL)*IMUL*ZS &
            + RMSPLS(IMUL)*IMUL*ZC
!
         ZJAC(2,1) = ZJAC(2,1) + DYMDROC(IMUL)*ZC &
            + DYMDROS(IMUL)*ZS
         ZJAC(2,2) = ZJAC(2,2) - YMSPLC(IMUL)*IMUL*ZS &
            + YMSPLS(IMUL)*IMUL*ZC
!
 151  continue
!
      return
      end

