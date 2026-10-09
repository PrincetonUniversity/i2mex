C-----------------------------------------------------------------------
C  TDB_SYMXMAP -- MAP ECE FREQUENCIES TO NORMALIZED RADIUS FOR PRESYMMETRIZ-
C    ATION OF ECE DATA
C
C  modularization -- DMC May 2005
C
C  R vector defined from frequency vector of ECE data, using vacuum map,
C  for presymmetrization purposes.  The data will be symmetrized to a
C  coordinate of relative midplane aspect ratio, but in the context of
C  ECE data this becomes ratio ((1/B1)-(1/B2))/((1/B1)+(1/B2)) and will
C  be interpreted at such when the presymmetrized data is used-- i.e.
C  by this means presymmetrization can be done early, in trdat, before
C  internal adjustments to the field strength are known...
C
C  at end the nominal R is converted to the usual x = (R-Rmaj)/Rmin
C
C  (DMC June 2005)
C
      SUBROUTINE TDB_SYMXMAP(d,zr1,zr2,zrbz,ILX,INX,IWORK)
C
      use trdatbuf_obj
      IMPLICIT NONE

      type (trdatbuf) :: d
      real*8, intent(in) :: zrbz ! external (R*Bz) used for B estimate here.
      real*8, intent(in) :: zr1,zr2 ! midplane bdy intercept locations
      integer, intent(in) :: ilx ! address of frequencies
      integer, intent(in) :: inx ! no. of frequencies
      integer, intent(in) :: iwork ! where to put corresponding midplane coord

      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
C  PASSED INPUT
C  ------------
C    ILX   -- ADDRESS OF ECE FREQUENCIES
C    INX   -- NUMBER OF ECE FREQUENCIES
C    IWORK -- ADDRESS WHERE CORRESPONDING RADII MAY BE STORED
C
C  OUTPUT TO WORKBUF(IWORK ... IWORK+INX-1)
C
C  note DMC July 1992 -- to support calls from TRDAT debug test code,
C    allow the possibility that the time TG1 and TG2 profile quantities
C    below have not been computed, in which case the arrays BIMDP1 and
C    BMIDP2 are assumed to contain the B field data, set up by the
C    TRDAT test code.
C
C-----------------------------------------------------------------------
C  local...
      integer :: ix
      real*8 :: zfreq,zbfreq,zrfreq,zzrmaj,zzrmin,zx
C
C-----------------------------------------------------------------------
C
C  MAP THE FREQUENCIES to vacuum field based "major radii"
C
      zzrmaj=(zr2+zr1)/2
      zzrmin=(zr2-zr1)/2
C
      DO 300 IX=1,INX           ! scan frequencies (low->high)
         ZFREQ=d%DATBUF(ILX+IX-1)
         ZBFREQ=ZFREQ/(28.0_R8*d%NHECFT)
         ZRFREQ=ZRBZ/ZBFREQ
C
C  SO THE CORRESPONDING X VALUE IS...
C
         ZX=(ZRFREQ-ZZRMAJ)/ZZRMIN
C
C  SAVE IN REVERSED ORDER
C
         d%WORKBUF(IWORK+INX-IX)=ZX ! (low->high; R ~ 1/f roughly)
C
 300  CONTINUE
C
      RETURN
      END
