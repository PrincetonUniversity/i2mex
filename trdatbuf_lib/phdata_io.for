C-----------------------------------------------------------------------
C  PHDA00 -- INITIALIZE READ/WRITE OPERATION ON PH.DAT FILE
C    ROUTINE IS CALLED FROM GENERATED SUBROUTINE PHDATA
C
      SUBROUTINE PHDA00(KLUN)
C
C  KLUN -- L.U.N. OF OPEN PH.DAT FILE
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
C---------------------------------------------------------------
C  COMMON BLOCKS SHARED BTW PHDAXX SUBROUTINES
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER klun
!============
      INTEGER IFLAG,IBUF(6)
      CHARACTER*16 ZNAMSV,ZNAM16
      CHARACTER*80 ZBUFSV,ZBUF
      COMMON/PHDAIO1/ ZNAMSV,ZNAM16,ZBUFSV,ZBUF
      COMMON/PHDAIO2/ IFLAG,IBUF
C
C---------------------------------------------------------------
C
      IFLAG=0
      ZNAMSV='$BEGIN  '
C
      RETURN
      END
C-----------------------------------------------------------------------
C  PHDAWR -- WRITE OPERATION ON PH.DAT FILE
C    ROUTINE IS CALLED FROM GENERATED SUBROUTINE PHDATA
C
      SUBROUTINE PHDAWR(KLUN,ZNAME,IVAL,INVAL)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inval,klun,iwrote,iwp1,ilft
!============
      CHARACTER*(*) ZNAME
      INTEGER IVAL(INVAL)
C
C  KLUN -- L.U.N. OF OPEN PH.DAT FILE
C  ZNAME -- CHARACTER*(*) NAME OF ITEM TO WRITE
C  IVAL(INVAL) -- DIMENSION AND VALUE OF ITEM (INTEGER ARRAY) TO WRITE
C
C---------------------------------------------------------------
C  COMMON BLOCKS SHARED BTW PHDAXX SUBROUTINES
      INTEGER IFLAG,IBUF(6)
      CHARACTER*16 ZNAMSV,ZNAM16
      CHARACTER*80 ZBUFSV,ZBUF
      COMMON/PHDAIO1/ ZNAMSV,ZNAM16,ZBUFSV,ZBUF
      COMMON/PHDAIO2/ IFLAG,IBUF
C
C---------------------------------------------------------------
C
      ZNAM16=ZNAME
      IF(ZNAM16(1:8).EQ.'ZZZ-END-') THEN
C  WRITE END OF HEADER MARKER
        WRITE(KLUN,'(A)') ZNAME
      ELSE
C  WRITE DATA
        IWROTE=0
 10     CONTINUE
        IF(IWROTE.LT.INVAL) THEN
          IF(IWROTE.EQ.0) THEN
C  WRITE 1ST RECORD INCLUDING ITEM NAME
            ZBUF=ZNAM16//':'
            CALL PHDAW1(KLUN,2,IVAL,INVAL)
            IWROTE=4
          ELSE
C  WRITE ADDL RECORDS TO INCLUDE ALL DATA, IF NECESSARY
            ZBUF=' '
            IWP1=IWROTE+1
            ILFT=INVAL-IWROTE
            CALL PHDAW1(KLUN,0,IVAL(IWP1),ILFT)
            IWROTE=IWROTE+6
          ENDIF
          GO TO 10
        ENDIF
      ENDIF
      ZNAMSV=ZNAM16
C
      RETURN
      END
C-----------------------------------------------------------------------
C  PHDAW1 -- WRITE DATA TO PH.DAT FILE
C
      SUBROUTINE PHDAW1(KLUN,IOFF,IVAL,INVAL)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ioff,inval,klun,ifld,ic,imax,inum,i,icx,il,ilnb
!============
      INTEGER IVAL(INVAL)
C
C  KLUN -- L.U.N. OF OPEN PH.DAT FILE
C  IOFF --  FIELD OFFSET (E.G. 1 MEANS SKIP FIRST FIELD)
C  IVAL(INVAL) -- DATA TO WRITE (INTEGERS)
C
C  THIS ROUTINE WRITES ONE LINE OF DATA, WHICH IS AT MOST 6
C  ITEMS ENCODED I13
C
C---------------------------------------------------------------
C  COMMON BLOCKS SHARED BTW PHDAXX SUBROUTINES
      INTEGER IFLAG,IBUF(6)
      CHARACTER*16 ZNAMSV,ZNAM16
      CHARACTER*80 ZBUFSV,ZBUF
      COMMON/PHDAIO1/ ZNAMSV,ZNAM16,ZBUFSV,ZBUF
      COMMON/PHDAIO2/ IFLAG,IBUF
C
C---------------------------------------------------------------
C
      IFLD=13
C
      IC=IOFF*IFLD+1
C
      IMAX=6-IOFF
C
      INUM=min(INVAL,IMAX)
C
      IC=IC-IFLD
      DO 10 I=1,INUM
        IC=IC+IFLD
        ICX=IC+IFLD-1
        WRITE(ZBUF(IC:ICX),'(I13)') IVAL(I)
 10   CONTINUE
C
      IL=LEN(ZBUF)
      DO 20 IC=IL,1,-1
        IF(ZBUF(IC:IC).NE.' ') GO TO 30
 20   CONTINUE
      IC=1
 30   ILNB=IC
C
      WRITE(KLUN,'(A)') ZBUF(1:ILNB)
C
      RETURN
      END
C-----------------------------------------------------------------------
C  PHDAWR1 -- WRITE OPERATION ON PH.DAT FILE
C    ROUTINE IS CALLED FROM GENERATED SUBROUTINE PHDATA
C
      SUBROUTINE PHDAWR1(KLUN,ZNAME,IVAL)
      IMPLICIT NONE
      INTEGER klun
      CHARACTER*(*) ZNAME
      INTEGER IVAL
      INTEGER IVAL2(1)
      IVAL2(1) = IVAL
      CALL PHDAWR(KLUN,ZNAME,IVAL2,1)
      RETURN
      END
C-----------------------------------------------------------------------
C  PHDARD -- READ OPERATION ON PH.DAT FILE
C    ROUTINE IS CALLED FROM GENERATED SUBROUTINE PHDATA
C
      SUBROUTINE PHDARD(KLUN,ZNAME,IVAL,INVAL,IER)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inval,ier,klun,i,iread,irp1,ilft
!============
      CHARACTER*(*) ZNAME
      INTEGER IVAL(INVAL)
C
C  KLUN -- L.U.N. OF OPEN PH.DAT FILE
C  ZNAME -- CHARACTER*(*) NAME OF ITEM TO READ
C  IVAL(INVAL) -- DIMENSION AND VALUE OF ITEM (INTEGER ARRAY) TO READ
C
C  IER -- RETURN CODE =0 IF OK OR CAUTIONS ONLY
C         =1 IF UNEXPECTED E-O-F
C         =2 IF OTHER ERROR
C
C---------------------------------------------------------------
C  COMMON BLOCKS SHARED BTW PHDAXX SUBROUTINES
      INTEGER IFLAG,IBUF(6)
      CHARACTER*16 ZNAMSV,ZNAM16
      CHARACTER*80 ZBUFSV,ZBUF
      COMMON/PHDAIO1/ ZNAMSV,ZNAM16,ZBUFSV,ZBUF
      COMMON/PHDAIO2/ IFLAG,IBUF
C
      integer :: lunmsg_tdb
C
C---------------------------------------------------------------
C  DEFAULT OUTPUT -- CLEAR VALUES TO ZERO
C
      DO 5 I=1,INVAL
        IVAL(I)=0
 5    CONTINUE
C
      ZNAM16=ZNAME
C
 10   CONTINUE
C  GET NEXT NAME FROM FILE...
      CALL PHDAR1(KLUN,1,IER)
      IF(IER.NE.0) GO TO 1000
C
      IF(ZNAMSV.EQ.ZNAM16) THEN
C  GO READ THE DATA
        GO TO 100
C
      ELSE IF(ZNAMSV.LT.ZNAM16) THEN
        WRITE(LUNMSG_TDB(0),1001) ZNAMSV
 1001 FORMAT(
     >' %PHDARD-- IGNORED UNRECOGNIZED NAME "',A,'" IN PH.DAT FILE')
        IFLAG=0
        GO TO 10
C
      ELSE
        WRITE(LUNMSG_TDB(0),1002) ZNAM16
 1002 FORMAT(
     >' %PHDARD-- DATA FOR "',A,'" NOT FOUND IN PH.DAT FILE')
        IER = 2
        GO TO 1000
      ENDIF
C
C  HAVE THE RIGHT NAME -- PROCESS INPUT DATA
C
 100  CONTINUE
      IF(ZNAM16(1:8).EQ.'ZZZ-END-') GO TO 1000
C
      IREAD=0
 110  CONTINUE
      IF(IREAD.LT.INVAL) THEN
        IF(IREAD.EQ.0) THEN
C  READ DATA ON SAME LINE WITH NAME
          CALL PHDAR2(KLUN,2,IVAL,INVAL,IER)
          IF(IER.NE.0) GO TO 1000
          IREAD=4
          IFLAG=0
        ELSE
C  READ ADDL RECORDS TO GET ALL DATA, IF NECESSARY
          IRP1=IREAD+1
          ILFT=INVAL-IREAD
          CALL PHDAR1(KLUN,0,IER)
          IF(IER.NE.0) GO TO 110
          IF(IFLAG.EQ.2) THEN
C  ONLY PROCESS IF THIS IS A DATA RECORD
            CALL PHDAR2(KLUN,0,IVAL(IRP1),ILFT,IER)
            IF(IER.NE.0) GO TO 1000
            IREAD=IREAD+6
            IFLAG=0
          ELSE
            WRITE(LUNMSG_TDB(0),1003) ZNAM16
 1003 FORMAT(
     >' %PHDARD -- MISSING PH.DAT DATA RECORD DETECTED, NAME="',A,'"')
            GO TO 1000
          ENDIF
        ENDIF
        GO TO 110
      ENDIF
C
C  EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  PHDARD1 -- READ OPERATION ON PH.DAT FILE
C    ROUTINE IS CALLED FROM GENERATED SUBROUTINE PHDATA
C
      SUBROUTINE PHDARD1(KLUN,ZNAME,IVAL,IER)
      IMPLICIT NONE
      INTEGER ier,klun
      CHARACTER*(*) ZNAME
      INTEGER IVAL
      INTEGER IVAL2(1)
      CALL PHDARD(KLUN,ZNAME,IVAL2,1,IER)
      IVAL = IVAL2(1)
      RETURN
      END
C-----------------------------------------------------------------------
C  PHDAR1 -- READ RECORD FROM PH.DAT FILE
C
      SUBROUTINE PHDAR1(KLUN,ICALL,IER)
C
C  KLUN -- L.U.N. OF OPEN PH.DAT FILE
C  ICALL -- =1:  LOOK FOR A RECORD WITH A NAME HEADER
C           =0:  JUST READ A RECORD
C
C  IER -- RETURN CODE =0 IF OK
C         =1 IF UNEXPECTED E-O-F
C         =2 IF OTHER ERROR
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
C---------------------------------------------------------------
C  COMMON BLOCKS SHARED BTW PHDAXX SUBROUTINES
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER icall,ier,klun
!============
      INTEGER IFLAG,IBUF(6)
      CHARACTER*16 ZNAMSV,ZNAM16
      CHARACTER*80 ZBUFSV,ZBUF
      COMMON/PHDAIO1/ ZNAMSV,ZNAM16,ZBUFSV,ZBUF
      COMMON/PHDAIO2/ IFLAG,IBUF
C
C---------------------------------------------------------------
C  ON ENTRY, IFLAG=1 MEANS I HAVE ALREADY READ A RECORD WITH A NAME
C  FIELD; IFLAG=2 MEANS I HAVE ALREADY READ A RECORD WITH DATA ONLY
C
 10   CONTINUE
      IF((ICALL.EQ.1).AND.(IFLAG.NE.1)) IFLAG=0
C
C  NOW IFLAG.GT.0 MEANS WE ALREADY HAVE READ THE DATA
C
      IER=0
      IF(IFLAG.GT.0) RETURN
C
      READ(KLUN,'(A)',END=202,ERR=203) ZBUF
C
      IF((ZBUF(1:1).GE.'A').AND.(ZBUF(1:1).LE.'Z')) THEN
        IFLAG=1
        ZNAMSV=ZBUF(1:16)
      ELSE
        IFLAG=2
      ENDIF
C
      GO TO 10
C
 202  CONTINUE
      IER=1
      RETURN
C
 203  CONTINUE
      IER=2
      RETURN
C
      END
C
C-----------------------------------------------------------------------
C  PHDAR2 -- READ/INTERPRET DATA FROM PH.DAT FILE
C
      SUBROUTINE PHDAR2(KLUN,IOFF,IVAL,INVAL,IER)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ioff,inval,ier,klun,ifld,ic,imax,inum,i,icx
!============
      INTEGER IVAL(INVAL)
C
C  KLUN -- L.U.N. OF OPEN PH.DAT FILE
C  IOFF -- FIELD OFFSET
C  IVAL(INVAL) -- DIMENSION AND VALUE OF ITEM (INTEGER ARRAY) TO READ
C
C  IER -- RETURN CODE =0 IF OK OR CAUTIONS ONLY
C         =1 IF UNEXPECTED E-O-F
C         =2 IF OTHER ERROR
C
C  READ INTEGER DATA BY DECODING 13 CHARACTER FIELDS
C  SEE ROUTINE PHDAW1 WHICH WROTE THE DATA BEING READ HERE
C
C---------------------------------------------------------------
C  COMMON BLOCKS SHARED BTW PHDAXX SUBROUTINES
      INTEGER IFLAG,IBUF(6)
      CHARACTER*16 ZNAMSV,ZNAM16
      CHARACTER*80 ZBUFSV,ZBUF
      COMMON/PHDAIO1/ ZNAMSV,ZNAM16,ZBUFSV,ZBUF
      COMMON/PHDAIO2/ IFLAG,IBUF
C
C---------------------------------------------------------------
      IER=0
C
      IFLD=13
C
      IC=IOFF*IFLD+1
C
      IMAX=6-IOFF
C
      INUM=min(INVAL,IMAX)
C
      IC=IC-IFLD
      DO 10 I=1,INUM
        IC=IC+IFLD
        ICX=IC+IFLD-1
        READ(ZBUF(IC:ICX),'(I13)',ERR=203) IVAL(I)
 10   CONTINUE
C
      RETURN
C
 203  CONTINUE
      IER=2
      RETURN
C
      END
C-----------------------------------------------------------------------
C******************** END FILE OLYMPS.FOR ; GROUP OLYMPS ***************
C*DECK Z2
! 19jan2003 fgtok -s r8_precision.sub all.sub "r8con.csh conversion"
