      subroutine ryevala(rmsplc,ymsplc,rmspls,ymspls, &
                         zsnthtk,zcsthtk,nmom,zr99,zy99)
use iso_c_binding, only: fp => c_double
!
!  evaluate moments position:  asymmetric eq. moments formula
!
      real rmsplc(0:nmom),ymsplc(0:nmom),rmspls(0:nmom),ymspls(0:nmom)
      real zcsthtk(nmom),zsnthtk(nmom)
!
!--------------------------------
!
      ZR99 = RMSPLC(0)
      ZY99 = YMSPLC(0)
      do 51 IMUL=1,NMOM
         ZR99 = ZR99 + RMSPLC(IMUL)*ZCSTHTK(IMUL) &
            + RMSPLS(IMUL)*ZSNTHTK(IMUL)
         ZY99 = ZY99 + YMSPLC(IMUL)*ZCSTHTK(IMUL) &
            + YMSPLS(IMUL)*ZSNTHTK(IMUL)
 51   continue
!
      return
      end
