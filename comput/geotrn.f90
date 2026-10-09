      subroutine GEOTRN(MODE, XIN, EIN, DIN, EOUT, doUT)
use iso_c_binding, only: fp => c_double
!***********************************************************************
!*****GEOTRN TRANSFORMS ELONGATION AND TRIANGULARITY FROM A MOMENTS TO *
!*****A GEOMETRICAL REPRESENTATION AND VICE VERSA.                     *
!*****REFERENCES:                                                      *
!*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *
!*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*
!*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *
!*****CALCULATED parameterS:                                           *
!*****EOUT-OUTPUT ELONGATION.                                          *
!*****doUT-OUTPUT TRIANGULARITY.                                       *
!*****INPUT parameterS:                                                *
!*****MODE-DESIGNATES DIRECTION OF TRANSFORMATION.                     *
!*****    =1 GEOMETRICAL => MOMENTS.                                   *
!*****    =2 MOMENTS => GEOMETRICAL.                                   *
!*****XIN-INPUT REDUCED MINOR RADIUS                                   *
!*****EIN-INPUT ELONGATION.                                            *
!*****DIN-INPUT TRIANGULARITY.                                         *
!***********************************************************************
      if (XIN.LE.0.0) goto 30
      if (MODE.EQ.2) goto 20
!*****GEOMETRICAL ==> MOMENTS REPRESENTATION.
      doUT = DIN/4.0
      CTC = 0.0
      do 10 I=1,10
           CTC = 4.0*doUT/(SQRT(XIN**2+32.0*DOUT**2)+XIN)
           doUT = XIN*DIN/(4.0-6.0*CTC**2)
   10 continue
      EOUT = XIN*EIN/(SQRT(1.0-CTC*CTC)*(XIN+2.0*doUT*CTC))
      return
!*****MOMENTS ==> GEOMETRICAL REPRESENTATION.
   20 CTC = 4.0*DIN/(SQRT(XIN**2+32.0*DIN**2)+XIN)
      STC = SQRT(1.0-CTC*CTC)
      S2TC = 2.0*STC*CTC
      C2TC = 2.0*CTC*CTC - 1.0
      EOUT = EIN*(XIN*STC+DIN*S2TC)/XIN
      doUT = DIN*(1.0-3.0*C2TC)/XIN
      return
   30 EOUT = EIN
      doUT = DIN
      return
      end
