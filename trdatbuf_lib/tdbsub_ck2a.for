C-----------------------------------------------------------------------
C  TDBSUB_CK2A -- CHECK THE POINTERS AND CONTROLS ON A TRDAT PROFILE
C
C  DMC MAY 1991 -- FOR NEW TRDAT/TRANSP INTERFACE
C
      SUBROUTINE TDBSUB_CK2A(d,ZDATID,IFPTR,IXPTR,IXLEN,INRI,INSY,IER)
C
C  PASSED INPUT:
C
C    ZDATID -- DATA PROFILE ID
C    IFPTR  -- PTR TO DATA -- NONZERO IF DATA EXISTS
C    IXPTR  -- PTR TO X AXIS IF DATA EXISTS
C    IXLEN  -- NO. PTS IN X AXIS IF DATA EXISTS
C    INRI   -- X AXIS TYPE CODE IF DATA EXISTS
C    INSY   -- DATA SYMMETRIZATION CODE IF DATA EXISTS
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      use trdatbuf_obj
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ixptr,ixlen,inri,insy,ier,ifptr,ia,iabs,ix,iwarn
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zxmin,zxmax,zxatest
!============
      type (trdatbuf) :: d
      character*(*) ZDATID
C
      integer :: nonlin,lunmsg_tdb
C
C-----------------------------------------------------------------------
C
      IER=0
C
      nonlin=lunmsg_tdb(0)
C
      IF(IFPTR.LE.0) GO TO 1000
C
      IF((IXPTR.LE.0).OR.(IXLEN.LE.0)) THEN
         IER=1
         WRITE(NONLIN,9901) ZDATID,IXPTR,IXLEN
 9901    FORMAT(
     >        ' ?TDBSUB_CK2A -- ERROR IN TRDAT? -- "',A,
     >        '" X AXIS DATA INVALID'/
     >        '   ADDRESS = ',I10,'  LENGTH = ',I10)
      ENDIF
C
      IA=IABS(INRI)
      IF((IA.EQ.0).OR.(IA.GT.8)) THEN
         IER=1
 
         WRITE(NONLIN,9902) ZDATID,INRI
 9902    FORMAT(
     >        ' ?TDBSUB_CK2A:  NAMELIST CONTROL NRI',A,
     >        ' = ',I10,' ** INVALID **')
      ENDIF
C
      IF((IA.EQ.3).AND.((INSY.LE.0).OR.(INSY.GT.4))) THEN
         IER=1
 
         WRITE(NONLIN,9903) ZDATID,INRI,ZDATID,INSY
 9903    FORMAT(
     >        ' ?TDBSUB_CK2A:  NAMELIST ERROR:  ',
     >        'SYMMETRIZATION CONTROL NOT SET'/
     >        '  NRI',A,' = ',I2,' ... BUT NSY',A,' = ',I10)
      ENDIF
C
      IF((IA.NE.3).AND.(INSY.NE.0)) THEN
         IER=1
 
         WRITE(NONLIN,9904) ZDATID,INRI,ZDATID,INSY
 9904    FORMAT(
     >        ' ?TDBSUB_CK2A:  NAMELIST ERROR:  ',
     >        'SYMMETRIZATION CONTROL *IS* SET'/
     >        '  NRI',A,' = ',I2,' (not 2-sided) ... BUT NSY',A,
     >        ' = ',I10)
      ENDIF
C
C  if IER=0 to this point, check X axis values and ordering also...
C
      if(ier.eq.0) then
         zxmin=1.0d34
         zxmax=-zxmin
         iwarn=0
         do ix=1,ixlen
            if(ix.gt.1) then
               if(d%DATBUF(ixptr+ix-1).le.d%DATBUF(ixptr+ix-2)) then
                  iwarn=iwarn+1
                  if(iwarn.eq.1) then
                     write(nonlin,*) ' ?TDBSUB_CK2A: "',trim(zdatid),
     >                    '" x axis not in strict ascending order:'
                     write(nonlin,
     >                    '(a,i4,a,i4,a,1pe12.5,a,i4,a,1pe12.5,a)')
     >                    '  nx=',ixlen,'  x(',ix-1,')=',
     >                    d%DATBUF(ixptr+ix-2),'  x(',ix,')=',
     >                    d%DATBUF(ixptr+ix-1),'.'
                  else if(iwarn.eq.2) then
                     write(nonlin,*) ' [...multiple violations...]'
                  endif
                  ier=1
               endif
            endif
            zxmin=min(zxmin,d%DATBUF(ixptr+ix-1))
            zxmax=max(zxmin,d%DATBUF(ixptr+ix-1))
         enddo
C
         zxatest = max(abs(zxmin),abs(zxmax))
         if(zxatest.gt.5.0_R8) then
            if((inri.ge.5).or.(inri.lt.0)) then
               ier=1
 
               write(nonlin,9910) zdatid,inri,zxatest
            endif
         endif
 9910    format(' ?TDBSUB_CK2A:  NRI',A,' = ',I2,
     >        ' specifies a normalized'/
     >        '  X axis (abs(x(i)) btw 0 and 1 + epsilon), but:'/
     >        '  max(abs(x(i))) = ',1pe11.4,' in data, i=1 to nx.')
C
      endif
C
 1000 CONTINUE
      RETURN
      END
