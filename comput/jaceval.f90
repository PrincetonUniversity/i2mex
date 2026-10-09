      subroutine jaceval(DR0DRO,DRMDRO,DYMDRO,RMSPL,YMSPL, &
                          ZSNTHTK,ZCSTHTK,NMOM,ZJAC)
use iso_c_binding, only: fp => c_double
!
!  evaluate 2x2 jacobian using moments and derivatives
!    ** symmetric formula **
!
      real rmspl(nmom),ymspl(nmom)
      real dr0dro,drmdro(nmom),dymdro(nmom)
      real zcsthtk(nmom),zsnthtk(nmom)
!
      real zjac(2,2)
!
!--------------------------------
!
      ZJAC(1,1)=DR0DRO
      ZJAC(1,2)=0.0
!
      ZJAC(2,1)=0.0
      ZJAC(2,2)=0.0
!
      do 150 IMUL=1,NMOM
!
         ZJAC(1,1)=ZJAC(1,1)+DRMDRO(IMUL)*ZCSTHTK(IMUL)
         ZJAC(1,2)=ZJAC(1,2)-RMSPL(IMUL)*IMUL*ZSNTHTK(IMUL)
!
         ZJAC(2,1)=ZJAC(2,1)+DYMDRO(IMUL)*ZSNTHTK(IMUL)
         ZJAC(2,2)=ZJAC(2,2)+YMSPL(IMUL)*IMUL*ZCSTHTK(IMUL)
!
 150  continue
!
      return
      end

