      subroutine ryeval(r0spl,rmspl,ymspl,zsnthtk,zcsthtk,nmom, &
                         zr99,zy99)
use iso_c_binding, only: fp => c_double
!
!  evaluate moments position:  symmetric eq. moments formula
!
      real r0spl,rmspl(nmom),ymspl(nmom)
      real zcsthtk(nmom),zsnthtk(nmom)
!
!--------------------------------
!
      ZR99=R0SPL
      ZY99 = 0.
      do 50 IMUL=1,NMOM
         ZR99=ZR99+RMSPL(IMUL)*ZCSTHTK(IMUL)
         ZY99=ZY99+YMSPL(IMUL)*ZSNTHTK(IMUL)
 50   continue
!
      return
      end
