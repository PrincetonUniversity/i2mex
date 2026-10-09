C-----------------------------------------------------------------------
C
      SUBROUTINE TDB_PROFLI(d, t, INRI,INSY,
     >        ILX,INX,ILF,
     >        ILXSY,INXSY,ILFSY,ILSSY,
     >        ierr)
C
      use trdatbuf_obj
      use trdatbuf_aux
      use trdatbuf_iface, only: tdb_xilmp, tdb_xisymp
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
C  dmc -- May 2005 -- modular implementation -- no trcom
C         f90 derived types are employed...
C
      type (trdatbuf) :: d
      type (profget) :: t
C
C  DMC -- **REDESIGNED INTERFACE** -- 17 MAR 1992
C   THE CODE IS BEING RESTRUCTURED FOR USE WITH THE OUTPUT OF THE NEW
C   VERSION TRDAT PROGRAM.  TRDAT OUTPUTS AND TRANSP COMMON NOW CONTAIN
C   THE SAME THINGS FOR EACH PROFILE INPUT DATA TYPE.  TDB_PROFLI PROCESSES
C   THE INTERPOLATION OF PROFILE DATA USING ALL OF THESE FORMS OF INFO
C   FOR ALL TYPES OF INPUT PROFILE DATA
C
C   INPUTS:  CONTROL OPTIONS, TIME TO INTERPOLATE TO, POINTERS TO
C     PROFILE DATA IN DATBUF
C
C   OUTPUTS:  THE interpolated DATA 
C
C   optional debug outputs, if t%idebug = T
C     DATA ASSYMETRY profile IF SYMMETRIZED HERE.
C     THE SHIFT PROFILE (IF SLICE AND STACK PRESYMMETRIZATION WAS USED)
C     AND MAPPING CHECKS VS. MAJOR RADIUS -- the original data vs. major
C     radius; the symmetrized data vs. major radius, for comparison.
C
C  DMC - 1986 - NEW ROUTINE - INTERPOLATE PROFILE DATA
C   COPIED FROM OLD CODE, PROFE1.FOR, JAN 1986.
C   GENERALIZED ROUTINE ALLOWS SAME TECHNIQUES TO BE USED IN
C   INTERPOLATING DATA OTHER THAN ELECTRON TEMPERATURE AND DENSITY.
C
C  SUBROUTINES OF TDB_PROFLI - CF PFSETR.FOR (DEFINE INTERPOLATION TARGET
C  DATA) AND INT2D.FOR (EXECUTE INTERPOLATION ON 2D ARRAY)
C
C  ARGUMENTS (INPUT)
C
C    INRI - TYPE OF NORMALIZED RADIAL VECTOR USED TO STORE THE 2D
C     DATA ARRAY (WHICH IS TO BE INTERPOLATED HERE)
C      =1:  MIDPLANE RADIUS - OUTSIDE HALF OF PLASMA
C          SHAFRANOV CORRECTION MAY BE EMPLOYED
C      =2:  MIDPLANE RADIUS - INSIDE HALF OF PLASMA
C      =3:  MIDPLANE RADIUS, 2 SIDED DATA, BOTH SIDES OF PLASMA
C      =4:  ABEL INVERTED DATA MINOR RADIUS; NO SHAFRANOV CORRECTION
C    DATA MAY ALSO BE TEMPERATURE VS ECE FREQUENCY; SEE IECE SWITCH...
C    NOTE:  NEGATIVE INRI VALUES MAY BE ENCOUNTERED IN TRDAT -- THESE
C      INDICATE AN INPUT UFILE WITH SPATIAL COORDINATE ALREADY NORM-
C      ALIZED; THE ABSOLUTE VALUE IS USED IN THIS ROUTINE
C      =5:  (new dmc 4 May 1993)  *** DATA HAS ALREADY BEEN MAPPED ***
C          TO normalized sqrt(tor.flux) "XI" coordinate.
C      =6:  data is mapped to normalized poloidal flux (no SQRT)
C      =7:  data mapped to normalized SQRT(pol.flux)
C      =8:  data mapped to normalized toroidal flux (no SQRT)
C
C    INSY - SYMMETRIZATION OPTION
C     VALID FOR INRI=3 ONLY:
C      =1:  SLICE AND STACK; DATA HAS BEEN PRESYMMETRIZED.
C      =2:  IN/OUT AVERAGE, PERFORMED HERE.
C      =3:  weighted IN/OUT AVERAGE, PERFORMED HERE ALSO
C      =4:  "R weighted in/out average".  Use the TRANSP MHD 
C           geometry to determine the locations of the two intercepts of each 
C           flux surface with the midplane. The weighted average of the profile 
C           at the inner and outer flux surfaces is given by
C               Navg = (Nin*Rin+Nout*Rout)/(Rin+Rout)
C           where Rin,Rout are the major radii at the inner and outer flux surface
C           midplane intercepts.  Same Notes applies as for NSYxxx=2.  This option
C           is mainly provided as a compatibility feature, since old versions of
C           TRANSP applied this weight when symmetrizing electron density (NER)
C           data and other density profiles, under NSYxxx=2.
C
C
C    in "t"
C     t%iece is .TRUE. IFF data is Te vs. ECE frequency 
C
C    t%ibdy = T TO INTERPOLATE DATA DIRECTLY TO ZONE BDY'S THEN FILL IN
C      ZONE CTRS  O.W. USE OLD ALGORITHM: INTERPOLATE DIRECTLY TO CTR'S
C      AND FILL IN BDY'S -- USE THIS FOR CHI(E) AND OTHER BDY - ORIENTED
C      DATA
C    +++> t%IBDY .AND. t%IECE NOT COMPATIBLE <+++
C
C    ILX - POINTER TO 2D STORED DATA X (RADIUS) PTS
C    INX - NUMBER OF X PTS IN 2D STORED DATA
C    ILF - POINTER TO 2D STORED DATA
C
C  THESE POINTERS POINT TO OUTPUT OF SLICE AND STACK PRESYMMETRIZER:
C    ILXSY - POINTER TO 2D STORED DATA X PTS -- PRESYMMETRIZED
C    INXSY - NUMBER OF X PTS IN 2D STORED DATA -- PRESYMMETRIZED
C    ILFSY - POINTER TO 2D STORED DATA -- PRESYMMETRIZED
C    ILSSY - POINTER TO 2D STORED DATA SHIFT PROFILE -- PRESYMMETRIZED
C
C
C
C    t%TIME - ** TIME ** TO INTERPOLATE DATA
C
C
C
C --------------
C  OUTPUT DATA
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER insy,ibdy,ilx,inx,ilf,ilxsy,inxsy,ilfsy
      INTEGER ilssy,inri,ierr,inria,iabs,j,ihece,iloop,iz,i,ilxuse
      INTEGER inxuse,ilfuse,iltuse,intime,il,izones,iixtrp,isgn
      INTEGER i2,iadr,icenr,i1,ioff,jp1,iorig,jm1
      integer :: nonlin,nout,lunmsg_tdb
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ztime,ztol,ztl1,ztl2,zslp,zf2,zf1
!============
C
C  WORKSPACE
C
      real*8, dimension(:), allocatable :: zdatw
C
      REAL*8 ZW1,ZW2,zr1,zr2,zr1p,zr2p
C
      INTEGER INOUT(2)
C
C  PFSETR ROUTINE ALSO USED PARTS OF COMMON ARRAY WORKBUF FOR WORKSPACE
C
      LOGICAL ILDBUG  ! .TRUE. FOR EXTRA DEBUG OUTPUT
      DATA ILDBUG/.FALSE./
C
C----------------------
C
C  EXECUTABLE CODE
C
C    PART 1.  COMPUTE ZDATA AND ZDATW -- PROFILES MAPPED TO TRANSP
C    FLUX SURFACES, AND AN OPTIONAL ASSYMETRY FACTOR
C
C    DEFINE THE INTERPOLATION TARGET POSITIONS
C    IN THE COORDINATES IN WHICH THE 2D DATA IS STORED...
C
      ztime = t%time

      IZONES=t%NZONES+1
      allocate(zdatw(izones)); zdatw = 0.0_R8

      if(t%ibdy) then
         ibdy = 1
      else
         ibdy = 0
      endif

      IERR=0

      if(t%idebug) then
         call tdb_xilmp(t%nzones,t%xibdys,t%xilmp)
         call tdb_xisymp(t%nzones,t%xibdys,t%rmajmp,t%xirsym,t%rmjsym)
      endif

      nonlin = lunmsg_tdb(0)
      nout = nonlin

      if (ilf .eq. 0 .and. ilfsy .eq. 0) then
         write(nonlin,*)
     >        'DATA POINTERS PASSED TO TDB_PROFLI ARE ZERO'
         ierr = 1
         goto 40100
      end if
C
      INRIA=IABS(INRI)
C
      IF((INRIA.LE.0).OR.(INRIA.GT.8)) THEN
         WRITE(NONLIN,9000) INRI
 9000    FORMAT(' ?TDB_PROFLI:  ILLEGAL INPUT DATA RADIAL TYPE INRI=',
     >        I12)
         IERR=1
      ENDIF
C
      IF((t%iece).AND.(IBDY.EQ.1)) THEN
        WRITE(NONLIN,9001)
 9001   FORMAT(
     >       ' ?TDB_PROFLI ERROR:  BOTH FLAGS IECE AND IBDY ARE ACTIVE')
        IERR=1
      ENDIF
C
      IF((INRIA.EQ.3).AND.((INSY.LT.1).OR.(INSY.GT.4))) THEN
        WRITE(NONLIN,9002) INSY
 9002   FORMAT(' ?TDB_PROFLI:  ILLEGAL SYMMETRIZATION CODE INSY=',I12)
        IERR=1
      ENDIF
C
      IF((INRIA.NE.3).AND.(INSY.NE.0)) THEN
        WRITE(NONLIN,9003) INSY
 9003   FORMAT(' ?TDB_PROFLI:  NON-ZERO SYMMETRIZATION CODE INSY=',I12)
        IERR=1
      ENDIF
C
      IF(IERR.NE.0) GO TO 40100
C
      IF(t%iece) THEN
        IHECE=d%NHECFT
      ELSE
        IHECE=0
      ENDIF
C
      IF(IHECE.GT.0) THEN
C  IF THE ECE DATA WAS MAPPED IN PRESYMMETRIZATION, DO NOT TRY
C  TO MAP IT AGAIN...!
        IF((INRIA.EQ.3).AND.(INSY.EQ.1)) IHECE=-1
      ENDIF
C
      CALL TDB_PFSETR(d,t,INRIA,INSY,IHECE,IBDY,ILOOP,INOUT)
      if(iloop.eq.0) then
         ierr=777
         go to 40100
      endif
C
C  CLEAR THE OUTPUT ARRAYS
C
      t%data_zc = 0.0_R8
      t%data_zb = 0.0_R8
C
      if(t%idebug) then
         t%asym_zc = 0.0_R8
         t%asym_zb = 0.0_R8
         t%shift = 0.0_R8
         t%datrsym = 0.0_R8
         t%datusym = 0.0_R8
      endif
C
C  ASSUME LINEAR INTERPOLATION IS OK
C
C  CHOOSE DATA SOURCE
      IF((INRIA.EQ.3).AND.(INSY.EQ.1)) THEN
C  PRESYMMETRIZED
        iorig=0
        ILXUSE=ILXSY
        INXUSE=INXSY
        ILFUSE=ILFSY
        ILTUSE=d%LTIME2
        INTIME=d%NTIME2
C  ERROR CHECK ON TIME RANGE:
        ZTL1=d%DATBUF(d%LTIME2)
        ZTL2=d%DATBUF(d%LTIME2+d%NTIME2-1)

        ztol = 1.0e-12_R8*max(1.0_R8,(ZTL2-ZTL1))

        if((ZTIME.LT.ZTL1).and.(ZTIME.ge.(ZTL1-ztol))) then
           ztime=ZTL1
        endif

        if((ZTIME.GT.ZTL2).and.(ZTIME.le.(ZTL2+ztol))) then
           ztime=ZTL2
        endif

        IF((ZTIME.LT.ZTL1).OR.(ZTIME.GT.ZTL2)) THEN
          WRITE(NONLIN,9901) ZTL1,ZTL2,ZTIME
          ierr=999
          go to 40100
        ENDIF

 9901   FORMAT(/
     >' ?TDB_PROFLI -- CALL TO INTERPOLATE PRESYMMETRIZED DATA IS OUT'/
     >'  OF BOUNDS IN TIME.  PRESYMMETRIZED DATA RANGE IS'/
     >'  ',1PE14.7,' TO ',1PE14.7,' SECONDS; INTERPOLATION IS'/
     >'  REQUESTED AT ',1PE14.7,' SECONDS.'/)

      ELSE
C  ORIGINAL
        iorig=1
        ILXUSE=ILX
        INXUSE=INX
        ILFUSE=ILF
        ILTUSE=d%LTIME2
        INTIME=d%NTIME2
      ENDIF
C
C  LOOP OVER INNER/OUTER HALF; OR ONE HALF ONLY -- SEE PFSETR ROUTINE
C
      DO IL=1,ILOOP
C
         ZSLP=0.0_R8
         IIXTRP=1
C
         if(iorig.eq.1) then
            CALL INT2D(d%DATBUF(ILTUSE),INTIME,
     1           d%DATBUF(ILXUSE),INXUSE,
     2           d%DATBUF(ILFUSE),INTIME,INXUSE,
     3           ZTIME,
     4           d%WORKBUF(d%LBX(IL):d%LBX(IL)+IZONES-1), IZONES,
     4           d%WORKBUF(d%LBX(3):d%LBX(3)+IZONES-1),
     5           ZSLP,
     6           IIXTRP, ILDBUG, IERR, NOUT)
         else
            CALL INT2D(d%DATBUF(ILTUSE),INTIME,
     1           d%DATBUF(ILXUSE),INXUSE,
     2           d%DATBUF(ILFUSE),INTIME,INXUSE,
     3           ZTIME,
     4           d%WORKBUF(d%LBX(IL):d%LBX(IL)+IZONES-1), IZONES,
     4           d%WORKBUF(d%LBX(3):d%LBX(3)+IZONES-1),
     5           ZSLP,
     6           IIXTRP, ILDBUG, IERR, NOUT)
         endif
C
         IF(IERR.NE.0) GO TO 40100
C
C		--------------------------------
C		1.6	STORE NEW DATA VALUES
C		--------------------------------
C
C  SET UP FOR EXTRACTION OF ASYMMETRY TERM.
C  SIGN=+1 FOR FIRST SIDE (IL=1), -1 FOR SECOND SIDE (IL=2)
C
         ISGN=3-2*IL
         ZF2=.5_R8*(ILOOP-1)
C
         ZF1=1.0_R8/ILOOP
C
         DO I = 1,izones
C  CHECK FOR INDEX ORDER REVERSAL BTW. MINOR AND MAJOR RADIAL COORD
            I2=I-1
            IF(INOUT(IL).EQ.2) I2=izones-I
            IADR=d%LBX(3)+I2
            if(ibdy.eq.0) then
               t%data_zc(i)=t%data_zc(i) + ZF1*d%WORKBUF(IADR)
            else
               t%data_zb(i)=t%data_zb(i) + ZF1*d%WORKBUF(IADR)
            endif
            if(t%idebug) then
               if((i.gt.1).or.(ibdy.eq.0)) then
                  zdatw(i)=zdatw(i)+ ISGN*ZF2*d%WORKBUF(IADR)
               endif
            endif
         enddo
C
C  END OF INNER/OUTER DATA-HALF LOOP
C
      enddo
C
C  IF DATA IS 2-SIDED ADJUST ZDATA TO BE THE VOLUME-WEIGHTED FLUX
C  SURFACE AVERAGE - ASSUMING A COS(POL. ANGLE) VARIATION IN DENSITY
C  ALONG THE FLUX SURFACE OR WHEN INSY=3, WEIGHT BASED ON FLUX
C  SURFACE SPACING.
C
      IF(INRIA.EQ.3) THEN
         IF (INSY.EQ.3) THEN
            DO J=1,izones
               if((ibdy.eq.0).or.(j.gt.1)) then
                  call profli_nsy3(J,ZW1,ZW2)
                  if(ibdy.eq.0) then
                     t%data_zc(j)=t%data_zc(j) + 
     >                    zdatw(j)*(ZW2-ZW1)/(ZW2+ZW1)
                  END IF
               endif
            END DO
         ELSE IF (INSY.EQ.4) THEN
            zr1=t%rmajmp(izones)
            zr2=zr1
            DO J=1,izones
               if((ibdy.eq.0).or.(j.gt.1)) then
                  if(ibdy.eq.1) then
                     zr1=t%rmajmp(izones-j+1)
                     zr2=t%rmajmp(izones+j-1)
                     t%data_zb(j)=t%data_zb(j)+
     >                    0.5_R8*zdatw(j)*(zr2-zr1)/(zr2+zr1)
                  else
                     zr1p=zr2
                     zr2p=zr2
                     zr1=t%rmajmp(izones-j+1)
                     zr2=t%rmajmp(izones+j-1)
                     t%data_zc(j)=t%data_zc(j)+
     >                    0.5_R8*zdatw(j)*(zr2+zr2p-(zr1+zr1p))/
     >                    (zr2+zr2p+zr1+zr1p)
                  endif
               endif
            enddo
         END IF
      ENDIF
C
C  FILL IN VALUES NOT INTERPOLATED FROM DATA
C   I.E. BDY'S IF CTR'S WERE GOTTEN FROM DATA
C   OR CTR'S IF BDY'S WERE GOTTEN FROM DATA
C
      if(t%idebug) then
         if(ibdy.eq.0) then
            t%asym_zc = zdatw
         else
            t%asym_zb = zdatw
         endif
      endif

      IF(IBDY.EQ.0) THEN
         t%data_zb(1)=t%data_zc(1)
         if(t%idebug) t%asym_zb(1)=0.0_R8
         do j=2,izones
            t%data_zb(j)=0.5_R8*(t%data_zc(j-1)+t%data_zc(j))
            if(t%idebug) then
               t%asym_zb(j)=0.5_R8*(t%asym_zc(j-1)+t%asym_zc(j))
            endif
         enddo
      else
         t%data_zc(izones)=t%data_zb(izones)
         if(t%idebug) t%asym_zc(izones)=t%asym_zb(izones)
         DO J=2,izones
            t%data_zc(j-1)=0.5_R8*(t%data_zb(j-1)+t%data_zb(j))
            if(t%idebug) then
               t%asym_zc(j-1)=0.5_R8*(t%asym_zb(j-1)+t%asym_zb(j))
            endif
         enddo
C
      ENDIF
C
C-----------------------------------------------------------------------
C  PART 2.  COMPUTE THE INFERRED SHIFT -- ACTUALLY, INTERPOLATE THE
C  OUTPUT OF THE PRESYMMETRIZER
C
      IF((INRIA.EQ.3).AND.(INSY.EQ.1).and.(t%idebug)) THEN
C
        ZSLP=0.0_R8
        CALL INT2D(d%DATBUF(d%LTIME2),d%NTIME2,
     1			d%DATBUF(ILXSY),INXSY,
     2			d%DATBUF(ILSSY),d%NTIME2,INXSY,
     3			ZTIME,
     4			t%xibdys, IZONES,
     4			t%shift,
     5			ZSLP,
     6			IIXTRP, ILDBUG, IERR, NOUT)
C
        IF(IERR.NE.0) GO TO 40100
C
      ENDIF
C
C-----------------------------------------------------------------------
C  PART 3.  COMPUTE THE DATA FOR THE OUTPUT SYMMETRIZATION / MAPPING
C  MULTIGRAPH  ("TECOM" STYLE MULTIGRAPH FOR RPLOT)
C
C  THE ORIGINAL DATA IS SIMPLY MAPPED TO THE TWO SIDED OUTPUT ARRAY
C
      if(t%idebug) then

         ICENR=2*t%nzones + 3
         t%datrsym(icenr) = t%data_zb(1)

         t%datrsym(1:2) = t%data_zc(izones)
         t%datrsym(t%nrsym-1:t%nrsym) = t%data_zc(izones)
C
         IOFF=0
         DO J=2,izones
            JM1=J-1
            IOFF=IOFF+1
            t%datrsym(ICENR+IOFF)=t%data_zc(JM1)
            t%datrsym(ICENR-IOFF)=t%data_zc(JM1)
            IOFF=IOFF+1
            t%datrsym(ICENR+IOFF)=t%data_zb(J)
            t%datrsym(ICENR-IOFF)=t%data_zb(J)
         enddo
      endif
C
C  THE UNMAPPED DATA WILL BE INTERPOLATED DIRECTLY FROM THE ORIGINAL
C  DATA -- WITH A POSSIBLE ECE MAP
C
      if(t%idebug) then
         if(t%iece) IHECE = d%nhecft ! even if data was presymmetrized...
         CALL TDB_UNMAP(d,t,INRIA,IHECE,ZTIME,ILX,INX,ILF,ierr)
      endif
C
40100 continue
      deallocate(zdatw)
      RETURN
 
      contains
 
      !
      ! --------------- profli_nsy3 --------------
      ! compute the weighting functions for INSY=3 option
      ! w1=Rin*dRin, w2=Rout*dRout
      !
      subroutine profli_nsy3(J,W1,W2)
      implicit none
 
      integer, intent(in) :: j     ! flux zone/boundary
      real*8, intent(out) :: w1,w2 ! inside and outside weight factors
      integer :: jrout, jop, jom   ! JR index for outer boundary and +-
      integer :: jrin,  jip, jim   ! JR index for inner boundary and +-
      real*8  :: zrout, zrin       ! major radius at outer and inner
      real*8  :: zdout, zdin       ! width of major radius at outer and inner
      integer :: jmax              ! max JR index
 
      jmax = 2*t%nzones+1
 
      jrout = (t%nzones+1)+(j-1)
      jop   = min(jmax,jrout+1)
      jom   = max(1,jrout-1)
 
      jrin  = (t%nzones+1)-(j-1)
      jip   = min(jmax,jrin+1)
      jim   = max(1,jrin-1)
 
      if (ibdy.eq.1) then
         ! boundary
         zrout = t%rmajmp(jrout)
         zrin  = t%rmajmp(jrin)
         zdout = (t%rmajmp(jop)-t%rmajmp(jom))/(jop-jom)
         zdin  = (t%rmajmp(jip)-t%rmajmp(jim))/(jip-jim)
      else
         zrout = 0.5D0*(t%rmajmp(jrout)+t%rmajmp(jop))
         zrin  = 0.5D0*(t%rmajmp(jrin)+t%rmajmp(jim))
         zdout = (t%rmajmp(jop)-t%rmajmp(jop-1))
         zdin  = (t%rmajmp(jim+1)-t%rmajmp(jim))
      end if
 
      w1 = zrin*zdin
      w2 = zrout*zdout
 
      end subroutine profli_nsy3
 
      END
