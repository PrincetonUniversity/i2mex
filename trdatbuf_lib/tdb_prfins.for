C----------------------------------------------------------------
C  TDB_PRFINS
C
C  INTERPOLATE SYMMETRIZED DATA TO DESIRED X COORDINATES
C   CONVERT <X> FCN TO A "SHAFRANOV SHIFT" PROFILE
C   INPUT (PASSED):
C     INUMR--  NUMBER OF PTS IN SYMMETRIZED UNINTERPOLATED DATA PROFILES
C     ZXF(INXF)-- TARGET X LOCATIONS AT WHICH TO INTERPOLATE PROFILES
C         (COMMON):
C     WORKBUF(LBX(2)..): DATA (UNINTERPOLATED)
C     WORKBUF(LBX(3)..): CORRESPONDING X LOCATIONS
C     WORKBUF(LBX(5)..): AVG <X> FCN (SEE COMMENTS IN CALLER, PRESYM)
C
C   OUTPUT (COMMON):
C     WORKBUF(LBX(7)..): DATA INTERPOLATED TO PTS ZXF(1..INXF)
C     WORKBUF(LBX(8)..): INTERPOLATED "EXPERIMENTAL FLUX SURFACE SHIFT"
C                       AT EACH PT ZXF(IX)
C
C---------
C
      SUBROUTINE TDB_PRFINS(d,zr1,zr2,IMAP,INUMR,ZXF,INXF)
C
      use trdatbuf_obj
      IMPLICIT NONE

      type (trdatbuf) :: d
      real*8, intent(in) :: zr1,zr2 ! bdy midplane intercepts

      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER inxf,inumr,inm1,i,im1,ix,ix1,ix2,iy1,iy2,is1,is2,ixm1
      integer :: inumleft,imap
      integer :: lunmsg_tdb
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zxf,zx,zlfac,rminrd
!============
      DIMENSION ZXF(INXF)
C
C------------------------------------------------------------
      real*8 :: za_bdy          ! aspect ratio at boundary
      real*8 :: zr0             ! major radius defined from boundary @midplane
      real*8 :: zrmin           ! minor radius defined from boundary @midplane
C
      real*8 :: za1,za2,zx1,zx2 ! aspect ratios; normalized aspect ratios
      real*8 :: zrmj_tmp,zrmn_tmp ! local major and minor radii
C
C------------------------------------------------------------
C
      zr0=(zr1+zr2)/2
      zrmin=(zr2-zr1)/2
      za_bdy=(zr2-zr1)/(zr2+zr1)
C
C   INPUT DATA VALUES MONOTONICALLY **INCREASING**
C   INPUT X PTS MONOTONICALLY **DECREASING**
C
      INM1=INUMR-1
      DO 10 I=1,INM1
         IM1=I-1
         IF(d%WORKBUF(d%LBX(2)+I).LE.d%WORKBUF(d%LBX(2)+IM1)) GO TO 9100
         IF(d%WORKBUF(d%LBX(3)+I).GE.d%WORKBUF(d%LBX(3)+IM1)) GO TO 9100
 10   CONTINUE
C
C  OK-- INSERT DATA IN STANDARD FORM
C
      inumleft=inumr
      DO 50 IX=1,INXF
         ZX=ZXF(IX)
C
         DO 40 I=inumleft,2,-1
C  ADRESSES  X, DATA(X), SHIFT(DATA(X))
            IX1=d%LBX(3)+I-1
            IX2=IX1-1
            IY1=d%LBX(2)+I-1
            IY2=IY1-1
            IS1=d%LBX(5)+I-1
            IS2=IS1-1
C
            zrmn_tmp=zrmin*d%workbuf(ix1)
            if(imap.eq.1) then
               zx1=zrmn_tmp/zrmin
            else
               zrmj_tmp=zrmin*d%workbuf(is1)+zr0
               za1=zrmn_tmp/zrmj_tmp
               zx1=za1/za_bdy
            endif
C
            zrmn_tmp=zrmin*d%workbuf(ix2)
            if(imap.eq.1) then
               zx2=zrmn_tmp/zrmin
            else
               zrmj_tmp=zrmin*d%workbuf(is2)+zr0
               za2=zrmn_tmp/zrmj_tmp
               zx2=za2/za_bdy
            endif
C
            IF((zx1.LE.ZX).AND.(zx2.GE.ZX))  GO TO 45
 40      CONTINUE
C     PT OUT OF RANGE
         if(ix.eq.inxf) then
            if(abs(zx2-zx).lt.1.0d-5) then
               i=2
               go to 45
            endif
         endif
         write(lunmsg_tdb(0),*) 
     >        '? SR TDB_PRFINS<PROFE2<PROFEX X PT OUT OF RANGE:'
         write(lunmsg_tdb(0),*) ' x-target ',inxf,zxf(1:inxf)
         GO TO 9500
C--/OK
 45      CONTINUE
C
         inumleft=i
         ZLFAC=(ZX-ZX1)/(ZX2-ZX1)
C
         IXM1=IX-1
         d%WORKBUF(d%LBX(7)+IXM1)=
     >        d%WORKBUF(IY1)+(d%WORKBUF(IY2)-d%WORKBUF(IY1))*ZLFAC
C  STORE DEDUCED SHAFRANOV SHIFT
         d%WORKBUF(d%LBX(8)+IXM1)=zrmin*
     >        (d%WORKBUF(IS1)+(d%WORKBUF(IS2)-d%WORKBUF(IS1))*ZLFAC)
 50   CONTINUE
      RETURN
C------
C  ERRORS
C
 9100 CONTINUE
      write(lunmsg_tdb(0),*) 
     >  '? SR TDB_PRFINS<PROFE2<PROFEX MONOTONICITY VIOLATION'
 9500 CONTINUE
      write(lunmsg_tdb(0),*) 'WORKBUF:F',inumr,
     >     d%WORKBUF(d%LBX(2):d%LBX(2)+inumr-1)
      write(lunmsg_tdb(0),*) 'WORKBUF:X',inumr,
     >     d%WORKBUF(d%LBX(3):d%LBX(2)+inumr-1)
      call bad_exit   
C
      END
