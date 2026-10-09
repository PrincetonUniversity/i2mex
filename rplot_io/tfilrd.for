C******************** START FILE TFILRD.FOR ; GROUP TFILIO *************
C=============================================================
C  TFILRD(LUN,IER)   READ TRANSP PLOTTING TF FILE (CONTAINING
C		 LABELING AND DIMENSIONING INFORMATION)
C		 ON LOGICAL UNIT LUN
C
C  dmc 14 Nov 1997 -- important change:  ixflag argument added.
C    meaning:  ixflag=0, function as before
C              ixflag.gt.0, read in label set for an auxilliary run
C
C    plan is to use the same TF.PLN read routine, whether reading
C      the labels for the "main" run (cf rplot) or a "secondary"
C      run (cf trprofil/trscalar, or rplot auxilliary run data input
C      options).
C
C    the auxilliary run labels now have space for them in CPLOTR, in
C    arrays ABR_X, ABT_X, etc., instead of ABR, ABT, etc.
C
C----------------------------------------------------------------
C
      SUBROUTINE TFILRD(LUN,IER,ixflag_in)
C
C  COMMON BLOCKS---
      use cplotr_mod
C
C  local stuff
C
      integer ixflag_in, ixflag
      logical my_skipid
      CHARACTER*3 FMT
      CHARACTER*10 ZRUNID
      logical ilxvar,ilxsave
C
      character*8 zint8
      character*10 zabr
      character*32 zunits
      character*64 zlabel
      integer idum(NAXMGF)
      real zdum(100)
      integer k
      character*132 zbuf
C
C----------------------------------------------------------------
C
C  CPLOTR initialization -- BLOCK DATA -- freeshare/cplset.for
C
C----------------------------------------------------------------
C
C  DMC JUNE 1988  -- CLEANUP -- REORGANIZE
C
C   MOD TO SUPPORT LONGER LABELS FOR RPLOT DATA FUNCTIONS
C
C----------------------------------------------------------------
C
C   OPEN FILE
C
      ixflag=abs(ixflag_in)
      if (ixflag_in .lt. 0) then
         my_skipid = .true.
      else
         my_skipid = .false.
      endif
      if(ier.lt.0) then
         iquiet=1
      else
         iquiet=0
      endif
C
      IF(FILNC.EQ.'N') THEN
        OPEN(UNIT=LUN,IOSTAT=IOS
     >      ,ACCESS='SEQUENTIAL',STATUS='OLD')
        IER=IOS
      ELSE IF(FILNC.EQ.'Y') THEN
        CALL TFILOP(LUN,TFILN,IER)
      ENDIF
C
      IF(IER.NE.0) GO TO 900
C
C  READ INTEGER DESCRIPTORS
CX	(old) CALL TFILHD(LUN,NRUN,RUNID,NSHOT,NZONES,NFT,NFR,NXR,0,FMT)
C
      CALL TFILHD(LUN,IRUN,ZRUNID,ISHOT,IZONES,IFT,IFR,IXR,0,FMT)
C
C  CHECK FOR OLD/NEW FORMAT BY SIGN OF VARIABLE "NXR"
C  GET NFXT= NO. OF PROFILE FCNS; FIX SIGN OF NXR
C
      IF(IXR.LT.0) THEN
         ILXVAR=.TRUE.
         IXR=-IXR
         READ(LUN,2000) IFXT
 2000    FORMAT(I5)
      ELSE
C  OLD FORMAT (FIXED GEOMETRY, FIXED PROFILE FCN SIZE)
         ILXVAR=.FALSE.
         IFXT=IFR
      ENDIF
C
      if(ixflag.eq.0) then
         NLXVAR=ILXVAR
         NXR=IXR
         NFXT=IFXT
         NFR=IFR
         NFT=IFT
         NZONES=IZONES
         NSHOT=ISHOT
         if (.not. my_skipid) RUNID=ZRUNID
         NRUN=IRUN
      else 
         ilxsave=NLXVAR
         NLXVAR=ILXVAR
         NXR_X(ixflag)=IXR
         NFXT_X(ixflag)=IFXT
         NFT_X(ixflag)=IFT
         if (.not. my_skipid) RUNID_X(ixflag)=ZRUNID
         NZONES_X(ixflag)=IZONES
         NFR_X(ixflag)=IFR
C 10/12/00 CAL: Handle "old style" PLN files
C               if called from trprofil: ixflag is always > 0
         if (nfr .eq. 0) NFR=IFR
      endif
C
C  READ PLASMA TORUS MAJOR AND MINOR RADIUS
      READ(LUN,2002) ZMAJOR,ZMINOR
 2002 FORMAT(2E10.4)
      if(ixflag.eq.0) then
         RMAJOR=ZMAJOR
         RMINOR=ZMINOR
      endif
C
C  OLD FORMAT...
C  READ ZONE CENTER AND ZONE BDY RADIAL ARRAYS
      IF(.NOT.ILXVAR) THEN
         DO 10 I=1,IZONES
            READ(LUN,2002) ZRZON,ZRBOUN
            if(ixflag.eq.0) then
               RZON(I)=ZRZON
               RBOUN(I)=ZRBOUN
            endif
 10      CONTINUE
C  PLOTR ABBREVIATIONS
         if(ixflag.eq.0) then
            XNDABB(1)='RZON '
            XNDABB(2)='RBOUN'
         endif
      ENDIF
C
C  READ LABELS, UNITS AND ABBREVIATIONS OF FUNCTIONS OF TIME
      DO 20 I=1,IFT
         CALL TFILR2(LUN,zlabel,zunits,zabr,
     1     IDUM1,IDUM2,1,FMT)
         if(ixflag.eq.0) then
            LABELT(I)=zlabel
            UNITST(I)=zunits
            ABT(I)=zabr
         else
            LABELT_X(I,ixflag)=zlabel
            UNITST_X(I,ixflag)=zunits
            ABT_X(I,ixflag)=zabr
         endif
 20   CONTINUE
C
C  READ LABELS, UNITS AND ZONE/BDY CODE, AND ABBREVIATIONS
C  OF FUNCTIONS OF TIME AND RADIUS
      DO 30 I=1,IFXT
         CALL TFILR2(LUN,zlabel,zunits,zabr,
     1     iitype,irofff,2,FMT)
         if(ixflag.eq.0) then
            LABELR(I)=zlabel
            UNITSR(I)=zunits
            ABR(I)=zabr
            ITYPR(I)=iitype
            NROFFF(I)=irofff
         else
            LABELR_X(I,ixflag)=zlabel
            UNITSR_X(I,ixflag)=zunits
            ABR_X(I,ixflag)=zabr
            ITYPR_X(I,ixflag)=iitype
            NROFFF_X(I,ixflag)=irofff
         endif
 30   CONTINUE
C
C  MULTIGRAPH PACKAGE LABELS -- only saved for "main" run
C
      read(lun,'(A)') zint8
      zint8=adjustr(zint8)
      read(zint8,'(I8)') ibal

      if(ixflag.eq.0) then
         NBAL=IBAL
      endif
C
      IF(IBAL.EQ.0) GO TO 100
C
      DO 50 IB=1,IBAL
         CALL TFILR2(LUN,zlabel,zunits,zabr,
     1     iintl,iinfb,3,FMT)
         call tfilrd_2008(lun,idum,iinfb)   ! bugfix dmc 30 Apr 2002
         if(ixflag.eq.0) then
            LABELB(IB)=zlabel
            UNITSB(IB)=zunits
            ABB(IB)=zabr
            IINTB(IB)=iintl
            INFB(IB)=iinfb
            do j=1,iinfb
               IFUNB(J,IB)=idum(j)
            enddo
C  BUG FIX D. MC CUNE 12 OCT 1982
            IF(IINTB(IB).EQ.1) UNITSB(IB)=UNITST(IABS(IFUNB(1,IB)))
         endif
 50   CONTINUE
C
 100  CONTINUE
C
      IF(.NOT.ILXVAR) THEN
C   SET NEW-STYLE COMMON VARIABLES FOR FORWARD COMPATIBILITY
         DO 90 I=1,IFXT
            if(ixflag.eq.0) then
               NROFFF(I)=I+1
            else
               NROFFF_X(I,ixflag)=I+1
            endif
 90      CONTINUE
C
         INXR=IXR
         IF(INXR.EQ.0) INXR=2
C
         DO 95 I=1,INXR
            if(ixflag.eq.0) then
               NZONEX(I)=IZONES
               NRECX(I)=1
               NLXFOT(I)=.FALSE.
            else
               NZONEX_X(I,ixflag)=IZONES
               NRECX_X(I,ixflag)=1
            endif
 95      CONTINUE
C
C  READ LABELS AND VALUES FOR INDEPENDANT VARIABLES OTHER THAN
C  TIME.  NXR=# OF SUCH VARIABLES; IF NXR=0 THEN WE ARE READING
C  OLD-STYLE FILE WHICH HAD FIXED # OF NON-TIME INDEPENDANT VAR-
C  IABLES,
C  IABLES EQUAL TO 2.  FIX UP COMMON TO WORK FOR OLD STYLE
C  FILES TOO.
C
         IF(IXR.NE.0) GO TO 130
         IXR=2
         if(ixflag.eq.0) then
            NXR=IXR
            XLAB(1)='RADIUS    '
            XLABU(1)='CM        '
            XLAB(2)='RADIAL BDY'
            XLABU(2)='CM        '
C
            DO 120 I=1,IZONES
               XARRY(I,1)=RZON(I)
               XARRY(I,2)=RBOUN(I)
 120        CONTINUE
         endif
C
         GO TO 190
C
 130     CONTINUE
C  LABELS AND VALUES FOR INDEPENDANT VARIABLES OTHER THAN TIME
         DO 170 I=1,IXR

            IF(FMT.EQ.'   ') THEN
               zlabel=' '
               zunits=' '
               read(LUN,2009) zlabel(1:20),zunits(1:10)
            ELSE IF(FMT.EQ.'F88') THEN
               zlabel=' '
               zunits=' '
               READ(LUN,2009) zlabel(1:32),zunits(1:16)
            else if(FMT.eq.'F09') then
               read(lun,'(1x,a/1x,a)') zlabel,zunits
            ENDIF

            if(ixflag.eq.0) then
               XLAB(I)=zlabel
               XLABU(I)=zunits
               CALL TFILND(XLABU(I))
            endif
 170     CONTINUE
C
         DO 180 I=1,IXR
            READ(LUN,2010) (zdum(j),j=1,izones)
            if(ixflag.eq.0) then
               do j=1,izones
                  XARRY(J,I)=zdum(j)
               enddo
            endif
 180     CONTINUE
C
 2009    FORMAT(A,A)
 2010    FORMAT(5(1PE10.3))
C
 190     CONTINUE
      ELSE
C--------------
C  VARIABLE X AXIS FORMAT
         DO 110 I=1,IXR
            IF(FMT.EQ.'   ') THEN
               zlabel=' '
               zabr=' '
               read(lun,2020) ifxi,iroffxi,irecxi,izonexi,
     1              zlabel(1:10),zabr(1:5)
            ELSE IF(FMT.EQ.'F88') THEN
               zlabel=' '
               read(lun,'(a)') zbuf
               k=index(zbuf,'x')
               if (k<=0 .or. k>4) then
                  read(zbuf,2020) ifxi,iroffxi,irecxi,izonexi, ! original F88 format
     1                 zlabel(1:32),zabr
               else
                  read(zbuf(k+1:),2021) ifxi,iroffxi,irecxi,izonexi,
     1                 zlabel(1:32),zabr                                    ! RGA, modified Jan2009 -- expect 'x' in first column
               end if
            ELSE IF(FMT.eq.'F09') THEN
               zlabel=' '
               read(lun,'(a)') zbuf
               k=index(zbuf,'x')
               read(zbuf(k+1:),2022) ifxi,iroffxi,irecxi,izonexi,zabr
               read(lun,'(6x,a)') zlabel
            ENDIF
            if(ixflag.eq.0) then
               NFX(I)=ifxi
               NROFFX(I)=iroffxi
               NRECX(I)=irecxi
               NZONEX(I)=izonexi
               XLAB(I)=zlabel
               XNDABB(I)=zabr
            else
               NZONEX_X(I,ixflag)=izonexi
               NRECX_X(I,ixflag)=irecxi
               xndabb_x(i,ixflag)=zabr
            endif
 110     CONTINUE
 2020    FORMAT(I4,1X,I4,1X,I3,1X,I4,1X,A,1X,A)
 2021    FORMAT(I6,1X,I6,1X,I6,1X,I8,1X,A,1X,A)   ! allow NZONEX>10000 (PSIRZ) and bump up the others
 2022    FORMAT(I6,1X,I6,1X,I6,1X,I8,1X,A) ! allow NZONEX>10000 (PSIRZ) and bump up the others
      ENDIF
C
C-----------------------------
C  EXIT
C
      CLOSE(UNIT=LUN,Status='KEEP')
C
      if(ixflag.ne.0) NLXVAR=ilxsave
C
      RETURN
C
C  ERRORS
C
 900  CONTINUE
      if(iquiet.eq.0) then
         WRITE(6,9001,iostat=ier) TFILN
 9001    FORMAT('  COULD NOT OPEN FILE:  ',A)
      endif
      IER=1
      RETURN
      END
C******************** END FILE TFILRD.FOR ; GROUP TFILIO ***************
      subroutine tfilrd_2008(lun,idum,iinfb)
      implicit NONE
 
      integer lun,iinfb
      integer idum(iinfb)
c
c  dmc 30 Apr 2008 -- modified this read to allow longer integers.
c  if the first non-blank character in the ascii file is "x", then
c  the longer integers are required; otherwise not.
c
c  dmc Jan 2009 -- expanded max multigraph list size to 60; multiple
c  ascii lines now needed to read list of IDs

      integer i,j,iswitch,istart,i1,i2
      character*80 buf
c
c------------------------------
c
      idum=0

      buf=' '
      i1 = -14

      do
         read(lun,'(A)') buf
c
         do i=1,len(buf)
            if(buf(i:i).ne.' ') then
               if(buf(i:i).eq.'x') then
                  i1=i1+15
                  i2=min(iinfb,i1+14)
                  iswitch=1
                  istart=i+1
               else if(buf(i:i).eq.'c') then
                  i1=i1+15
                  i2=min(iinfb,i1+14)
                  iswitch=2  ! should be another line to read...
                  istart=i+1
               else
                  iswitch=0
                  istart=1
               endif
               exit
            endif
         enddo
 
         if(iswitch.eq.0) then
            read(buf(istart:),2008) (idum(J),J=1,iinfb)
 2008       FORMAT(15I4)        ! old way
            exit
         else
            read(buf(istart:),2009) (idum(J),J=i1,i2)
 2009       format(15i5)        ! new way
            if(iswitch.eq.1) exit
         endif
      enddo
 
      return
      end
