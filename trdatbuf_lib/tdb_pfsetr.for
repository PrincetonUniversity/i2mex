C-----------------------------------------------------------------------
C  TDB_PFSETR
C
C  SET UP (OPTIONALLY SHIFTED) INTERPOLATION TARGET X VECTOR
C  FOR ELECTRON TEMPERATURE OR DENSITY DATA -- or other profile data
C
C  dmc May 2005 -- modular version: no trcom...
C
C  dmc Feb 1997 -- support INRI=6 for data vs. normalized poloidal flux
C  dmc Feb 2002 -- support INRI=7 for data vs. sqrt(normalized poloidal flux)
C  dmc Feb 2002 -- support INRI=8 for data vs. normalized toroidal flux
C        (INRI=5 supports data vs. sqrt(normalized toroidal flux), as always)
C
C  dmc revised Oct 1995 -- to tolerate interior B(R) non-monotonicity
C    singularity
C
C  NEW VERSION DMC OCT 1983 -- X COORDINATE FOR GENERALIZED GEOMETRY
C
C  REVISED DMC MAR 1992 -- FOR USE WITH PROFLI
C    ARGUMENTS VARIABLES ARE FOR THE MOST PART DESCRIBED IN PROFLI.FOR
C
C    INPUTS:  INRI,INSY,IHECFT(THE ECE HARMONIC IF >0),IBDY
C    OUTPUTS:  ILOOP,INOUT
C
      SUBROUTINE TDB_PFSETR(d,t,INRI,INSY,IHECFT,IBDY,ILOOP,INOUT)
C
      use trdatbuf_obj
      use trdatbuf_aux
      IMPLICIT NONE

      type (trdatbuf) :: d
      type (profget) :: t

      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER insy,ihecft,ibdy,iloop,inri,ishift,izb,inria,il,jz,j
      INTEGER jp1,iadr,izp1,inm1,jr,istart,nrmaj,nzones,iarat
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zdr0,zsign,zshaf,zpf,zfp0,zfp,zdr,zr,zfbdy,zfn,zftol
      real*8 :: rminrd
!============
      INTEGER INOUT(2),lunmsg_tdb,nonlin
C
C  REV DMC FEB 1987
C   IF IBDY.EQ.1 SETUP FOR DATA INTERPOLATION TO ZONE BDY'S NOT ZONE
C   CTR'S
C
C  REV DMC JUNE 1984
C
C  IF IHECFT.GT.0 THEN:
C------------------------------
C  INTERPRET STORED DATA AS ECE DATA VS. FREQUENCY
C
C  IHECFT IS THE HARMONIC
C
C  COMMON VECTOR RMAJMP CONTAINS A SERIES OF MIDPLANE RADII COVERING
C  THE ENTIRE PLASMA; COMMON VECTOR FBX WHEN MULTIPLIED WITH THE
C  EXTERNAL FIELD YIELDS THE LOCAL FIELD STRENGTH AT THESE MAJOR RADIAL
C  LOCATIONS
C
C   THIS FEATURE MAY BE USED ONLY IF LEVGEO.GT.1 AND AT LEAST A
C  SHAFRANOV SHIFTED CIRCLES EQUILIBRIUM IS CALCULATED, INCLUDING A
C  PARA/DIAMAGNETIC EFFECT ON THE TOROIDAL FIELD
C
C   SUBROUTINE FBCALM WAS CALLED TO CALCULATE THE RATIO FBX OF |B|
C  THE TOTAL FIELD TO THE EXTERNAL FIELD (BZXR/R) AT EACH
C  MAJOR RADIUS POINT ALONG THE MIDPLANE AS STORED IN ARRAY RMAJMP
C
C   THE TOTAL FIELD STRENGTH INCLUDES THE EFFECT OF THE POLOIDAL
C  FIELD INDUCED BY TOROIDAL PLASMA CURRENT, AS WELL AS THE PARA/DIA
C  MAGNETIC EFFECT DUE TO POLOIDAL PLASMA CURRENT.
C
C  dmc Oct 18 1995:  in presence of strong enough pressure local Beta,
C   the B(R) gradient can reverse, causing a gap in Te coverage.  Now
C   fill this gap with interpolation instead of aborting.  Save the
C   gap indices here
C
C  dmc Nov 8, 1998:  now use pfsetsng to patch the ECE map -- old
C   ztdb_pfsetr COMMON and patching code was removed.
C
C----------------------------
      REAL*8, dimension(:), allocatable :: zwk1,zwk2,zx,drshaf
C----------------------------
C
C  PASSED INPUT ARGUMENTS:
C   INRI-- TYPE OF X VECTOR OF DATA TO BE INTERPOLATED (FROM COMMON
C     NRITER (TEMPERATURE) OR NRINER (DENSITY))
C     =1 OR 4:  DATA COVERS RANGE 0 TO 1 OF X = (R-R0)/A
C         IF ECE DATA - COVERS OUTSIDE HALF OF PLASMA MIDPLANE
C         =4:  data was originally vs. "minor radius".
C     =2     :  DATA COVERS RANGE -1 TO 0  OF X
C         IF ECE DATA - COVERS INSIDE HALF OF PLASMA MIDPLANE
C     =3     :  DATA COVERS RANGE -1 TO +1 (I.E. SYMMETRIZABLE
C               2-SIDED DATA)
C     =5     :  DATA HAS ALREADY BEEN MAPPED TO "XI" FLUX SURFACES
C     =6     :  data is mapped on surfaces labeled by normalized
C               poloidal magnetic flux
C     =7     :  data is mapped on surfaces labeled by normalized
C               SQRT(poloidal magnetic flux)
C     =8     :  data is mapped on surfaces labed by normalized
C               toroidal flux (no SQRT).
C
C  PASSED OUTPUT ARGUMENTS:
C   ILOOP= LOOP LIMIT =1 IF DATA IS SINGLE SIDED; 2 IF 2-SIDED
C     IJ=1 TO ILOOP ...
C   INOUT(IJ) FLAG FOR IJ'TH TARGET X VECTOR
C     =1 -- COVERS RANGE 0 TO 1 AND IS MONOTONIC INCREASING
C     =2 -- COVERS RANGE -1 TO 0 AND IS MONOTONIC INCREASING X:
C           INDEX REVERSAL NEEDED TO GET INCREASING MINOR RADIUS INDEX
C  ECE DATA - =1: COVERS INSIDE HALF OF MIDPLANE, MONOTONIC INCREASING
C             =2:  COVERS OUTSIDE HALF, IS MONOTONIC INCREASING IN
C             FREQUENCY SPACE - REVERSE INDICES TO GET MONOTONIC
C             INCREASING MINOR RADIUS ZONE INDEX
C
C  COMMON OUTPUT:
C   WORKBUF(LBX(IJ)...)  IS THE IJ'TH TARGET X VECTOR
C
C  MONOTONIC INCREASING TARGETS ARE REQUIRED FOR INTERPOLATION ROUTINE
C  INT2D USED IN PROFE1
C
C----------------------------------------------------------------------
C
C  FLUX SURFACE SHIFT EFFECT IS INCLUDED UNLESS
C   (1) THE DATA IS ABEL INVERTED, STORED VS. NORMALIZED MINOR RADIUS
C   (2) THE DATA HAS BEEN PRESYMMETRIZED AND MAPPED TO NORMALIZED
C         MINOR RADIUS
C
C  ISHIFT IS NOT USED IN CASE OF ECE DATA MAPPING
C
      nonlin=lunmsg_tdb(0)
C
      nrmaj=t%nrmaj
      nzones = t%nzones
      izp1 = nzones + 1
C
      iarat=0
      IF(IHECFT.LE.0) THEN
         ISHIFT=1
         IF(INRI.GE.4) ISHIFT=0
         IF((INRI.EQ.3).AND.(INSY.EQ.1)) then
                                !  data was presymmetrized...
            ISHIFT=0
            iarat=1             ! use normalized aspect ratio
            if(ihecft.lt.0) then
                                !  presymmetrized data was ECE originally;
               iarat=-1         ! use (1/B) version of aspect ratio
            endif
         endif
      ENDIF
C
      IF((INRI.EQ.3).AND.(INSY.EQ.2 .or. INSY.EQ.3 .or. INSY.EQ.4)) THEN
        ILOOP=2
        INOUT(1)=1
        INOUT(2)=2
      ELSE
        ILOOP=1
        IF(INRI.EQ.2) THEN
          INOUT(1)=2
        ELSE
          INOUT(1)=1
        ENDIF
        IF(IHECFT.GT.0) INOUT(1)=3-INOUT(1)
      ENDIF
C
C  NORMALIZED POSITION (X) SPACE MAP:
      allocate(zx(izp1),drshaf(izp1))
      CALL TDB_ROVERA(t,INRI,iarat,ZX,drshaf,izp1)
      if(ishift.eq.1) then
         rminrd=0.5_R8*(t%rmajmp(nrmaj)-t%rmajmp(1))
      endif
C
      IZB=2-IBDY  ! ZONE/BDY SHIFT INDEX
C
      inria=abs(inri)
      IF(IHECFT.LE.0) THEN
C------------------------------------------------
        if((inria.lt.6).or.(inria.gt.7)) then
           ZDR0=0.0_R8          ! may need to adjust for free bdy at some pt.
           DO 100 IL=1,ILOOP
C
              DO 99 JZ=1,izp1
C IN/OUT MAPPING
                 IF(INOUT(IL).EQ.1) THEN
                    ZSIGN=1.0_R8
                    J=JZ
                 ELSE
                    ZSIGN=-1.0_R8
                    J=izp1+1-JZ
                 ENDIF
C NORMALIZED SHIFT
                 IF(ISHIFT.EQ.1) THEN
C  MHD EQUILIBRIUM
                    zshaf = (zdr0 + drshaf(j))/rminrd
                 ELSE
C  NO SHIFT
                    ZSHAF=0.0_R8
                 ENDIF
C ADDRESS OF TARGET X
                 IADR=d%LBX(IL)+JZ-1
C TARGET X
                 d%WORKBUF(IADR)=ZSIGN*ZX(J)+ZSHAF
                 if(inria.eq.8) d%WORKBUF(iadr)=
     >              d%WORKBUF(iadr)*d%WORKBUF(iadr)
C LOOP END
 99           CONTINUE
 100       CONTINUE
        else
C-------------------------------------------
C  poloidal flux mapping (inria=6,7)
C
           if(iloop.ne.1) then
              write(nonlin,*)
     1           ' ?tdb_pfsetr:  inria=6 or 7 & iloop.ne.1 unexpected.'
              iloop=0
              go to 999
           endif
           do j=1,izp1
              iadr=d%LBX(iloop)+j-1
              if(ibdy.eq.1) then
C     pol. flux @ bdy
                 zpf=t%plflxg(j)
              else
C     pol. flux @ zone ctr
                 jp1=j+1
                 if(j.le.nzones) then
                    zpf=0.5_R8*(t%plflxg(j)+t%plflxg(jp1))
                 else
C  extrapolate
                    zpf=t%plflxg(izp1)+
     >                   0.5_R8*(t%plflxg(izp1)-t%plflxg(nzones))
                 endif
              endif
C  normalize
              zpf=max(0.0_R8,zpf)
              zpf=zpf/t%plflxg(izp1)
C  store
              d%WORKBUF(iadr)=zpf
              if(inria.eq.7) d%WORKBUF(iadr)=sqrt(zpf)
           enddo
        endif
C-------------------------------------------
C  FREQUENCY SPACE MAP
      ELSE
C
        IZP1=NZONES+1
C
        ZFP0=28.0_R8*IHECFT*t%BMIDP(IZP1)
C
        allocate(zwk1(izp1),zwk2(izp1))
C
        DO 200 IL=1,ILOOP
C
C  ZFBDY IS A FREQUENCY PAST THE EDGE OF THE PLASMA;
C  ZFP IS THE FREQUENCY AT THE PREVIOUS ZONE BDY
C
C  IN THIS DOUBLE LOOP A MONOTONICALLY INCREASING SEQUENCE OF
C  FREQUENCIES ARE DEFINED WHICH CORRESPOND TO TRANSP ZONE CTRS.
C  HOWEVER THE ORDERING WILL BE REVERSED (I.E. MOVING FROM EDGE TO
C  CENTER OF PLASMA) IF INOUT(IL)=2 AND THE FREQUENCIES CORRESPOND
C  TO THE OUTER HALF OF THE PLASMA MIDPLANE INTERSECTION
C
           ZFP=ZFP0
C
           IF(INOUT(IL).EQ.1) THEN
              ZFBDY=28.0_R8*IHECFT*
     >             (t%bmidp(1)-0.5_R8*(t%bmidp(2)-t%bmidp(1)))
              d%WORKBUF(d%LBX(IL)+nzones)=ZFBDY
           ELSE
              INM1=NRMAJ-1
              ZFBDY=28.0_R8*IHECFT*
     >             (t%bmidp(nrmaj)+
     >             0.5_R8*(t%bmidp(nrmaj)-t%bmidp(inm1)))
              d%WORKBUF(d%LBX(IL))=ZFBDY
           ENDIF
C
           DO 199 J=1,nzones
C
              IF(INOUT(IL).EQ.1) THEN
                 IADR=d%LBX(IL)+J-1
                 JR=izp1-J
              ELSE
                 IADR=d%LBX(IL)+izp1-J
                 JR=IZP1+J
              ENDIF
C
              ZFN=28.0_R8*IHECFT*t%bmidp(JR)
              d%WORKBUF(IADR)=0.5_R8*(ZFN+ZFP)
C
              ZFP=ZFN
C
 199       CONTINUE
C
           zftol=1.000001_R8
           istart=0
C
           call tdb_pfsetsng(d%WORKBUF(d%LBX(il):d%LBX(il)+izp1-1),zx,
     1          izp1,zwk1,zwk2,'tdb_pfsetr',nonlin)
C
C  the pfsetsng routine now enforces monotonicity of the ECE
C  freq->R map, so the non-monotonicity patch loop previously
C  coded is no longer needed.  It turned out to be a continuity
C  hazard with respect to Te (and hence Zeff in certain operational
C  modes).
C
 200    CONTINUE
        deallocate(zwk1,zwk2)
C------------------------------------------------
      ENDIF
C
 999  continue
      deallocate(zx,drshaf)
      RETURN
      END
C-----------------------------------------------------------------------
C
C  this routine removes flat spots, from an expected monotonic function,
C  preventing any slope less than a fixed amount below the average
C  slope from happening anywhere.
C  dmc 8 Nov 1998
C
      subroutine tdb_pfsetsng(zfdata,zx,inum,zwk1,zwk2,runid,nonlin)
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nonlin,inum,itmax,iter,ict,i,ip
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zfac,zmark,zsmin,zs,zdf
!============
      REAL*8 zfdata(inum)		! frequencies
      REAL*8 zx(inum)			! indep. coordinate
      REAL*8 zwk1(inum),zwk2(inum)	! work arrays
      character*(*) runid		! runid label
c
      data zfac/0.1_R8/
      data zmark/-1.0E35_R8/
      data itmax/1000/
c
c------------------------------------------
c
      zsmin=zfac*(zfdata(inum)-zfdata(1))/(zx(inum)-zx(1))
c
      iter=0
 10   continue
      iter=iter+1
      if(iter.gt.itmax) then
         write(nonlin,9901) itmax
 9901    format(
     >        ' ?tdb_pfsetr:  more than ',i5,' iterations in pfsetsng')
         call bad_exit
      endif
      ict=0
      zwk1(inum)=zfdata(inum)
      do i=1,inum-1
         zwk1(i)=zfdata(i)
         zwk2(i)=zmark
         ip=i+1
         zs=(zfdata(ip)-zfdata(i))/(zx(ip)-zx(i))
         if(zs.lt.zsmin) then
            ict=ict+1
            zwk2(i)=0.5_R8*(zfdata(i)+zfdata(ip))
         endif
      enddo
      if(ict.eq.0) go to 100
      do i=1,inum-1
         if(zwk2(i).ne.zmark) then
            ip=i+1
            zs=(zfdata(ip)-zfdata(i))/(zx(ip)-zx(i))
            zdf=0.5_R8*(1.0001_R8*zsmin-zs)*(zx(ip)-zx(i))
            zdf=max(zdf,0.00001_R8*max(abs(zwk1(i)),abs(zwk1(ip))))
            zwk1(i)=zwk1(i)-zdf
            zwk1(ip)=zwk1(ip)+zdf
         endif
      enddo
c
      do i=1,inum
         zfdata(i)=zwk1(i)
      enddo
      go to 10
c
 100  continue
      return
      end
