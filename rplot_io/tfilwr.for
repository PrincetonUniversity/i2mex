C******************** START FILE TFILWR.FOR ; GROUP TFILIO *************
C=============================================================
C  MODIF 28 JUL 81  D. MC CUNE
C   A SECOND PAIR OF ROUTINES "TFILW2" AND "TFILR2" HAVE BEEN
C   ADDED TO SUPPORT PASSING ADDED INFORMATION NEEDED FOR
C   MULTIPLE X-AXES OF VARYING LENGTHS.
C
C  NEW INFO-- RECORD POSITION OFFSET DATA FOR FCNS OF TIME + 1
C  ADDL COORDINATE.  LABELING AND SIZE DATA FOR VARIOUS X-AXES,
C  AND X-AXIS ABBREVIATIONS WHICH INDICATE WHICH PHYSICAL VARIABLE
C  IN THE MF FILE DEFINES THE X-AXIS PHYSICAL VALUES AND UNITS.
C    THE OLD ROUTINES ARE MAINTAINED FOR BACKWARDS COMPATIBILITY
C  TO OUTPUT OF TRANSP RUNS PRE 9-1-81... D. MC CUNE
C
C-------
C  TFILWR(LUN)   WRITE TRANSP PLOTTING TF FILE (CONTAINING
C		 LABELING AND DIMENSIONING INFORMATION)
C		 ON LOGICAL UNIT LUN
C
      SUBROUTINE TFILWR(LUN,ier)
C
C  COMMON BLOCKS  (PLOTTING SYSTEM) ---
      use cplotr_mod
      implicit none
C
      CHARACTER*3 FMT
C
      integer ifx(naxxvr)
      integer ier, ios,  ioff,  i,  ityp, inrec
      integer IXR, IDUM, IDUM2, IR, IB,   J,    LUN
C----------------------------------------------------
C  DMC JUNE 1988--
C    CLEANUP.  SUPPORT NEW TF FILE FORMAT WITH LONGER LABELS
C    CONSOLIDATE OLD TFILWR/TFILW2 ROUTINES AND USE F77 IF/THEN/ELSE
C
C----------------------------------------------------
C
C  CREATE AND OPEN FILE
      IF(FILNC.EQ.'N') THEN
C  unix:  check for old copy of file and delete it if it exists
        print *, '?tfilwr: FILNC=="N" is not allowed on unix'
        call bad_exit
      ELSE IF(FILNC.EQ.'Y') THEN
C  unix:  check for old copy of file and delete it if it exists
        open(unit=lun,file=tfiln,status='old',iostat=ios)
        if(ios.eq.0) then
          close(unit=lun,status='delete')
        else
          close(unit=lun)
        endif
        open(unit=lun,file=tfiln,status='new',iostat=ier)
        if(ier.ne.0) return
      ENDIF
C
      FMT='F09'  ! 2009 DMC, replacing 'F88' format
c  mainly: changes in tfilr2/tfilw2 subroutines
c
c  dmc -- reconstruct NFR, NROFFF,NROFFX if necessary
c    (for xfrevert TF.PLN file regen from NetCDF or MDSplus data)
c
      ifx=nfx
      if(nlxvar) then
         if(nfr.eq.0) then
            if(nzones.eq.0) nzones=nzonex(1)
            ifx=0
            ioff=2                      ! 1st "record" -> time
            do i=1,nfxt
               nrofff(i)=ioff
               ityp=itypr(i)
               ifx(ityp)=ifx(ityp)+1
               inrec=1+(nzonex(ityp)-1)/nzones
               ioff=ioff+inrec
               if(abr(i).eq.xndabb(ityp)) then
                  nroffx(ityp)=nrofff(i)
                  nrecx(ityp)=inrec
               endif
            enddo
            nfr=nrofff(nfxt)-1
         endif
      endif
 
C
C   CALL NEW VERSION "TFILW2"
C    IF FLEXIBLE X-AXIS FEATURES ARE IN USE
C
      IF(NLXVAR) THEN
        IXR=-NXR
      ELSE
        IXR=NXR
      ENDIF
C
C  WRITE INTEGER DESCRIPTORS
      CALL TFILHD(LUN,NRUN,RUNID,NSHOT,NZONES,NFT,NFR,IXR,1,FMT)
C
      IF(NLXVAR) THEN
C
C  FLEXIBLE GEOMETRY & X AXES...
C
C  WRITE ACTUAL # OF FCNS OF TIME + 1 ADDL COORDINATE (NFR IS NOW
C  THE TOTAL # OF MF FILE RECORDS CONTAINING SUCH DATA WHICH IS
C  NO LONGER NECESSARILY THE SAME AS THE # OF FCNS SINCE FCNS OF
C  DIFFERENT LENGTH ARE NOW SUPPORTED)
C
        WRITE(LUN,2000) NFXT
 2000   FORMAT(I5)
        WRITE(LUN,2002) RMAJOR,RMINOR
      ELSE
C  WRITE PLASMA TORUS MAJOR AND MINOR RADIUS (OLD FIXED GEO RUNS)
        WRITE(LUN,2002) RMAJOR,RMINOR
C  WRITE ZONE CENTER AND ZONE BDY RADIAL ARRAYS
        DO 10 I=1,NZONES
          WRITE(LUN,2002) RZON(I),RBOUN(I)
 10     CONTINUE
      ENDIF
 2002 FORMAT(2E10.4)
C
      IDUM=0
      IDUM2=0
C
C  WRITE LABELS, UNITS AND ABBREVIATIONS OF FUNCTIONS OF TIME
      DO 20 I=1,NFT
        CALL TFILW2(LUN,LABELT(I),UNITST(I),ABT(I)(1:10),
     >                IDUM,IDUM2,1)
 20   CONTINUE
C
C  WRITE LABELS, UNITS AND ZONE/BDY CODE, AND ABBREVIATIONS
C  OF FUNCTIONS OF TIME AND RADIUS
      IF(NLXVAR) THEN
        IR=NFXT
      ELSE
        IR=NFR
      ENDIF
C
      DO 30 I=1,IR
        CALL TFILW2(LUN,LABELR(I),UNITSR(I),ABR(I)(1:10),
     >                ITYPR(I),NROFFF(I),2)
 30   CONTINUE
C
C  MULTIGRAPH PACKAGES
      WRITE(LUN,2005) NBAL
 2005 FORMAT(1x,I5)
C
      IF(NBAL.EQ.0) GO TO 100

c  DMC Jan 2009 -- support MG sets of size .gt. 15 -- write out multiple
c  lines of member IDs now.  New limit is 60 -- for now.
c  RGA Jun 2011 -- how about 126 and expand if then block the lazy way

      DO 50 IB=1,NBAL
        CALL TFILW2(LUN,LABELB(IB),UNITSB(IB),
     >      ABB(IB),IINTB(IB),INFB(IB),3)
        if(infb(ib).le.15) then
           WRITE(LUN,2008) (IFUNB(J,IB),J=1,INFB(IB))
        else if(infb(ib).le.30) then
           WRITE(LUN,2018) (IFUNB(J,IB),J=1,15)
           WRITE(LUN,2008) (IFUNB(J,IB),J=16,INFB(IB))
        else if(infb(ib).le.45) then
           WRITE(LUN,2018) (IFUNB(J,IB),J=1,15)
           WRITE(LUN,2018) (IFUNB(J,IB),J=16,30)
           WRITE(LUN,2008) (IFUNB(J,IB),J=31,INFB(IB))
        else if(infb(ib).le.60) then
           WRITE(LUN,2018) (IFUNB(J,IB),J=1,15)
           WRITE(LUN,2018) (IFUNB(J,IB),J=16,30)
           WRITE(LUN,2018) (IFUNB(J,IB),J=31,45)
           WRITE(LUN,2008) (IFUNB(J,IB),J=46,INFB(IB))
        else if(infb(ib).le.75) then
           WRITE(LUN,2018) (IFUNB(J,IB),J=1,15)
           WRITE(LUN,2018) (IFUNB(J,IB),J=16,30)
           WRITE(LUN,2018) (IFUNB(J,IB),J=31,45)
           WRITE(LUN,2018) (IFUNB(J,IB),J=46,60)
           WRITE(LUN,2008) (IFUNB(J,IB),J=61,INFB(IB))
        else if(infb(ib).le.90) then
           WRITE(LUN,2018) (IFUNB(J,IB),J=1,15)
           WRITE(LUN,2018) (IFUNB(J,IB),J=16,30)
           WRITE(LUN,2018) (IFUNB(J,IB),J=31,45)
           WRITE(LUN,2018) (IFUNB(J,IB),J=46,60)
           WRITE(LUN,2018) (IFUNB(J,IB),J=61,75)
           WRITE(LUN,2008) (IFUNB(J,IB),J=76,INFB(IB))
        else if(infb(ib).le.105) then
           WRITE(LUN,2018) (IFUNB(J,IB),J=1,15)
           WRITE(LUN,2018) (IFUNB(J,IB),J=16,30)
           WRITE(LUN,2018) (IFUNB(J,IB),J=31,45)
           WRITE(LUN,2018) (IFUNB(J,IB),J=46,60)
           WRITE(LUN,2018) (IFUNB(J,IB),J=61,75)
           WRITE(LUN,2018) (IFUNB(J,IB),J=76,90)
           WRITE(LUN,2008) (IFUNB(J,IB),J=91,INFB(IB))
        else if(infb(ib).le.120) then
           WRITE(LUN,2018) (IFUNB(J,IB),J=1,15)
           WRITE(LUN,2018) (IFUNB(J,IB),J=16,30)
           WRITE(LUN,2018) (IFUNB(J,IB),J=31,45)
           WRITE(LUN,2018) (IFUNB(J,IB),J=46,60)
           WRITE(LUN,2018) (IFUNB(J,IB),J=61,75)
           WRITE(LUN,2018) (IFUNB(J,IB),J=76,90)
           WRITE(LUN,2018) (IFUNB(J,IB),J=91,105)
           WRITE(LUN,2008) (IFUNB(J,IB),J=106,INFB(IB))
        else                                                      ! see CPLOTR NAXMGF   
           WRITE(LUN,2018) (IFUNB(J,IB),J=1,15)
           WRITE(LUN,2018) (IFUNB(J,IB),J=16,30)
           WRITE(LUN,2018) (IFUNB(J,IB),J=31,45)
           WRITE(LUN,2018) (IFUNB(J,IB),J=46,60)
           WRITE(LUN,2018) (IFUNB(J,IB),J=61,75)
           WRITE(LUN,2018) (IFUNB(J,IB),J=76,90)
           WRITE(LUN,2018) (IFUNB(J,IB),J=91,105)
           WRITE(LUN,2018) (IFUNB(J,IB),J=106,120)
           WRITE(LUN,2008) (IFUNB(J,IB),J=121,min(126,INFB(IB)))  ! same as CPLOTR NAXMGF
           if(infb(ib).gt.126) then
              write(6,*) ' ** TFILWR warning: infb(ib).gt.126: '  ! same as CPLOTR NAXMGF
     &             ,infb(ib)
              write(6,*) ' ** multigraph truncated!'
           endif
        endif

cxx 2008   FORMAT(15I4)

c  "continuation line" implemented -- DMC Jan 2009 -- MG limit incr to 60
 2008   format('x',15i5)                ! modified DMC 30 Apr 2002
 2018   format('c',15i5)                ! DMC Jan 2009 -- #fcns > 15 possible

 50   CONTINUE
C
 100  CONTINUE
C
      IF(NLXVAR) THEN
C
C  DMC "F09" format adjustments here...
C
C  WRITE OUT X-AXIS SPECIFICATIONS
C  (QUANTITY NXR ALREADY WRITTEN OUT, FIRST LINE)
        DO 110 I=1,NXR
        WRITE(LUN,2007) ifx(I),NROFFX(I),NRECX(I),NZONEX(I),
     >	       XNDABB(I),XLAB(I)
 110    CONTINUE
 2007   FORMAT('x',I6,1X,I6,1X,I6,1X,I8,1X,A/6x,A)
      ELSE
C
C  LABELS AND VALUES FOR INDEPENDANT VARIABLES OTHER THAN TIME
        DO 170 I=1,NXR
          WRITE(LUN,2009) XLAB(I),XLABU(I)
 170    CONTINUE
        DO 180 I=1,NXR
          WRITE(LUN,2010) (XARRY(J,I),J=1,NZONES)
 180    CONTINUE
 2009   FORMAT(1x,A/1x,A)
 2010   FORMAT(5(1PE10.3))
      ENDIF
C
      CLOSE(UNIT=LUN,STATUS='KEEP')
C
      RETURN
     	END
C******************** END FILE TFILWR.FOR ; GROUP TFILIO ***************
