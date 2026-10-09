C	**************** START FILE  TRSCALAR.FOR  ;  GROUP  TRREAD_LIB ******
 
 
      SUBROUTINE TRSCALAR( DISK,
     1			     DIR,     RUNIDIN, ZABBREV,MAXTIMES,
     2			     LABEL,   UNITS,
     3                       NTIMES,  TIMES,   SDATA,  IERR     )
 
      use tconnect_mod
 
C	SUBROUTINE TRSCALAR READS IN A SPECIFIED SCALAR FUNCTION AND ASSOCIATED
C			    TIMES, LABEL AND UNITS FROM A TRANSP RUN.
C
C  ** NOTE NOTE NOTE **
C
C     for performance reasons it is desirable to suppress the trailing
C     MDS+ tree close, if another read in the same tree is to follow
C     immediately.  To get this behaviour, set IERR=-99 on *just before
C     calling* this routine.  If there is NO ERROR (IERR=0 on exit), the
C     tree is left open for subsequent calls.  If there is an error, the
C     tree is closed.  On the last of a series of calls be sure IERR is
C     not -99 so that the trailing MDSplus close does occur.
C
C	ARGUMENTS:
C		   NOTE 1: UNIT LUN (91) IS USED TO READ IN TRANSP FILES.
C		   NOTE 2: UNIT LUOUT (6) IS USED TO OUTPUT ERROR MESSAGES.
C
C	LAST CHANGED:
C
C           11/18/97  DMC  ** rewrite ** for NetCDF support and performance
C                          use the RPLOT buffering mechanism
C                          the user interface is unchanged.
C
C           11/04/96  TBT  HAPPY 21ST BIRTHDAY JENNIFER
C                          Corrected Format 2011 for 8/13 Change
c            8/13/96  dmc  RUNID increased to 10 chars; ZTEMP2 -> 14 chars
C            6/12/95  tbt  NOTE: On a unix system, DISK is ignored
C	     6/02/93  TbT  Added call to TrCaps and arg Zabbrev.
C            2/02/92  tbt  Added check for FMT.eq.' ' to read in ru 8225.
C	     1/22/91  tbt  Changed test for MAXTIMES. Added label 97&98.
C	     3/11/92  TBT  Add capability to read Old files such as
C                          RUNDATA:[transp.pdx.82]9712.
C
C	   INPUT:  DISK     - C*64 DISK NAME WHERE DATA IS LOCATED.
C			      IF DISK=' ', ASSUME NAME IS LOGICAL DISK "RUNDATA
C                             If Unix system, ignore DISK.
C
C		   DIR      - C*64 FILE DIRECTORY, E.G. 'TRANSP.TFTR.86
C			      NO '[' OR ']' ALLOWED IN DIRECTORY NAME.
C		   RUNIDIN  - C*(*) INPUT RUNID
C		   RUNID    - C*8 VARIABLE DEFINING THE RUNID, E.G. '4058  '.
C			      DIR & RUNID WILL BE USED TO MAKE TRANSP FILE NAMES
C	                       E.G. 'RUNDATA:[TRANSP.TFTR.86]4058MF.PLN    '
C			            'RUNDATA:[TRANSP.TFTR.86]4058TF.PLN    '
C			            'RUNDATA:[TRANSP.TFTR.86]4058NF.PLN    '
C	           ZABBREV  - C*(*max 10) ABBREVIATION OF THE SCALAR TO READ IN
C		   MAXTIMES - INTEGER DIMENSION OF ARRAY TIMES.
C
C	   RETURNED:
C		   LABEL    - LABEL FOR SCALAR (was C*32 now can be longer)
C		   UNITS    - UNITS OF SCALAR (was C*16 now can be longer)
C                  NTIMES   - INTEGER # OF TIMES READ IN FROM TRANSP FILE.
C		   TIMES    - REAL ARRAY OF DIMENSION MAXTIMES INTO WHICH
C		              TIMES ARE READ FROM DISK.
C		   SDATA    - REAL ARRAY OF DIMENSION MAXTIMES INTO WHICH NTIMES
C			      SCALAR VALUES ARE READ.
C		   IERR     - INTEGER ERROR MESSAGE. RETURNED = 0 IF OK.
C			      0 = OK
C			      2 = ERROR READING NF FILE.
C			      3 = UNEXPECTED EOF READING NF FILE.
C			      5 = TRIED TO READ OLD FORMAT TF FILE.
C			     15 = SCALAR ABBREVIATION NOT FOUND
C			     16 = COULDN'T OPEN NF FILE
C			     17 = COULDN'T OPEN TF FILE
C                            21 = TRSCALAR: - FMT .NE. blank or F88
C			    102 = # OF SCALARS > DIMENSION OF BUFFER.
C		                  NO DATA READ IN.
C			    101 = # OF TIMES IN  FILE > MAXTIMES.
C				  MAXTIMES TIMES WERE READ IN.
 
      IMPLICIT NONE
C
      external TFILRD		! link the label reader
      external TRFUNID          ! link this one too
C
      CHARACTER*(*) DIR
      CHARACTER*(*)  DISK
      character*150   diskl     ! or MDS+ path
      CHARACTER*(*) RUNIDIN
      CHARACTER*10   RUNID
      Character*(*) Zabbrev
      CHARACTER*10   ABBREV
      CHARACTER*(*)  LABEL
      CHARACTER*(*)  UNITS
      INTEGER	      MAXTIMES
      INTEGER       NTIMES
      REAL	      TIMES(MAXTIMES)
      REAL          SDATA(MAXTIMES)
      INTEGER       IERR
 
      integer idrun
      integer i,ierr0
 
      INTEGER  LUOUT,lunzer
 
      integer ildsk,ildir,ilrun
 
C	-----------------------------------------------------------------------
	
      IERR0 = IERR
C
      LUOUT = lunzer(0)
C
      Abbrev = Zabbrev
      Call TrCaps(Abbrev)     ! Make all capitals
 
      IERR  = 0     ! INITIALIZE
      LABEL = ' '
      UNITS = ' '
      NTIMES= 0
      RUNID = RUNIDIN   ! SET C*(*) TO C*10
C
      DISKL=DISK
C
      DO I=1,MAXTIMES
          TIMES(I) = 0
          SDATA(I) = 0.
      END DO  ! I
 
      ildsk=max(1,len_trim(diskl))
      ildir=max(1,len_trim(dir))
      ilrun=len_trim(runid)
 
      call tconnect(luout, diskl, dir, runid, idrun, ierr)
      if(ierr.ne.0) then
         write(luout,9901) diskl(1:ildsk),dir(1:ildir),runid(1:ilrun)
 9901    format(' %trscalar:  failed to connect to run:'/
     1     '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
      else
C
         call tgetscal(luout, idrun,abbrev,maxtimes,label,units,
     1      ntimes,times,sdata,ierr)
         if(ierr.ne.0) then
            write(luout,9902) abbrev,diskl(1:ildsk),dir(1:ildir),
     1           runid(1:ilrun)
 9902       format(' %trscalar:  failed to read named scalar:  ',a/
     1         '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
         endif
      endif
C
C  If MDSplus: close secondary run    ! CAL 1/28/00
      if (disk(1:4) .eq. 'MDS+') then
         if((ierr.ne.0).or.(ierr0.ne.-99)) then
            write(luout,*) ' %close MDSplus tree for 2nd run'
            call tconnect_close(ierr)
            if(ierr.ne.0) then
               write(luout,9903) abbrev,diskl(1:ildsk),dir(1:ildir),
     1              runid(1:ilrun)
 9903          format(' %trscalar:  failed to close MDS Tree ',a/
     1            '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
            endif
         endif
      endif
C
      return
      end

C
C ======================================================================
C
      SUBROUTINE TRSCALAR_CONNECT( DISK,
     1			     DIR,     RUNIDIN, ZABBREV,
     2			     IDRUN, LABEL, UNITS, IFCN,
     3                       NTIMES,  IMDS, IERR     )
 
      use tconnect_mod
 
C	SUBROUTINE TRSCALAR_CONNECT CONNECTS AND READS IN THE DIMENSIONS OF 
C                           A SPECIFIED SCALAR FUNCTION OF TIME AS WELL AS THE ASSOCIATED
C			    LABEL AND UNITS FROM A TRANSP RUN.  IT IS EXPECTED
C                           THAT TRSCALAR_FETCH WILL BE CALLED AFTER THIS
C                           FUNCTION IF IERR==0 OR IERR==-99.
C
C  ** NOTE NOTE NOTE **
C
C     Set IERR=0 or -99 on entering this subroutine.  If there is no error
C     than IERR will not change and the tree will remain open for the 
C     subsequent call to TRSCALAR_FETCH.  TRSCALAR_FETCH will then close the
C     tree if IERR=0 or leave it open if IERR=-99.
C
C
C	ARGUMENTS:
C		   NOTE 1: UNIT LUN (91) IS USED TO READ IN TRANSP FILES.
C		   NOTE 2: UNIT LUOUT (6) IS USED TO OUTPUT ERROR MESSAGES.
C
C	LAST CHANGED:
C
C
C	   INPUT:  DISK     - C*64 DISK NAME WHERE DATA IS LOCATED.
C			      IF DISK=' ', ASSUME NAME IS LOGICAL DISK "RUNDATA
C                             If Unix system, ignore DISK.
C
C		   DIR      - C*64 FILE DIRECTORY, E.G. 'TRANSP.TFTR.86
C			      NO '[' OR ']' ALLOWED IN DIRECTORY NAME.
C		   RUNIDIN  - C*(*) INPUT RUNID
C			      DIR & RUNID WILL BE USED TO MAKE TRANSP FILE NAMES
C	                       E.G. 'RUNDATA:[TRANSP.TFTR.86]4058MF.PLN    '
C			            'RUNDATA:[TRANSP.TFTR.86]4058TF.PLN    '
C			            'RUNDATA:[TRANSP.TFTR.86]4058NF.PLN    '
C	           ZABBREV  - C*(*max 10) ABBREVIATION OF THE SCALAR TO READ IN
C
C	   RETURNED:
C                  IDRUN    - INTEGER Index of run in CPLOTR
C		   LABEL    - LABEL FOR SCALAR (was C*32 now can be longer)
C		   UNITS    - UNITS OF SCALAR (was C*16 now can be longer)
C                  IFCN     - INTEGER index of function
C                  NTIMES   - INTEGER # OF TIMES READ IN FROM TRANSP FILE.
C                  IMDS     - INTEGER NONZERO IF THIS IS AN MDSPLUS CONNECTION
C		   IERR     - INTEGER ERROR MESSAGE. RETURNED = 0 IF OK.
C			      0 = OK
C			      2 = ERROR READING NF FILE.
C			      3 = UNEXPECTED EOF READING NF FILE.
C			      5 = TRIED TO READ OLD FORMAT TF FILE.
C			     15 = SCALAR ABBREVIATION NOT FOUND
C			     16 = COULDN'T OPEN NF FILE
C			     17 = COULDN'T OPEN TF FILE
C                            21 = TRSCALAR: - FMT .NE. blank or F88
C			    102 = # OF SCALARS > DIMENSION OF BUFFER.
C		                  NO DATA READ IN.
C			    101 = # OF TIMES IN  FILE > MAXTIMES.
C				  MAXTIMES TIMES WERE READ IN.
 
      IMPLICIT NONE
C
      external TFILRD		! link the label reader
      external TRFUNID          ! link this one too
C
      CHARACTER*(*)  DIR
      CHARACTER*(*)  DISK
      character*150  diskl     ! or MDS+ path
      CHARACTER*(*)  RUNIDIN
      CHARACTER*10   RUNID
      Character*(*)  Zabbrev
      CHARACTER*10   ABBREV
      CHARACTER*(*)  LABEL
      CHARACTER*(*)  UNITS

      INTEGER       IDRUN
      INTEGER       IFCN
      INTEGER       NTIMES
      INTEGER       IMDS
      INTEGER       IERR
 
      integer i,ierr0,ierrc
 
      INTEGER  LUOUT,lunzer
 
      integer ildsk,ildir,ilrun
 
C	-----------------------------------------------------------------------
	
      IERR0 = IERR
C
      LUOUT = lunzer(0)
C
      Abbrev = Zabbrev
      Call TrCaps(Abbrev)     ! Make all capitals
 
      IERR  = 0     ! INITIALIZE
      LABEL = ' '
      UNITS = ' '
      IFCN  = 0
      NTIMES= 0
      RUNID = RUNIDIN   ! SET C*(*) TO C*10
C
      DISKL=DISK
C
      if (disk(1:4) .eq. 'MDS+') then
         IMDS=1
      else
         IMDS=0
      end if

      ildsk=max(1,len_trim(diskl))
      ildir=max(1,len_trim(dir))
      ilrun=len_trim(runid)
 
      call tconnect(luout, diskl, dir, runid, idrun, ierr)

      if(ierr.ne.0) then
         idrun=-1
         write(luout,9901) diskl(1:ildsk),dir(1:ildir),runid(1:ilrun)
 9901    format(' %trscalar_connect:  failed to connect to run:'/
     1     '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
      else
C
         call tgetscal_connect(luout, idrun,abbrev,label,units,
     1      ifcn,ntimes,ierr)

         if(ierr.ne.0) then
            idrun=-1
            write(luout,9902) abbrev,diskl(1:ildsk),dir(1:ildir),
     1           runid(1:ilrun)
 9902       format(' %trscalar_connect:  failed to '
     1           //'read named scalar:  ',a/
     2           '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
         endif
      endif
C
C  If MDSplus: close secondary run    ! CAL 1/28/00
      if (imds/=0) then
         if((ierr.ne.0)) then
            write(luout,*) ' %close MDSplus tree for 2nd run'
            call tconnect_close(ierrc)
            if(ierrc.ne.0) then
               write(luout,9903) abbrev,diskl(1:ildsk),dir(1:ildir),
     1              runid(1:ilrun)
 9903          format(' %trscalar_connect:  '
     1              //'failed to close MDS Tree ',a/
     2              '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
            endif
         endif
      endif

      if (ierr0==-99 .and. ierr==0) ierr=-99   ! restore for trscalar_fetch
C
      return
      end
C
C =======================================================================
C
      SUBROUTINE TRSCALAR_FETCH_CSTRING( IDRUN, ZCABBREV, IFCN, IMDS,
     1                           MAXTIMES, TIMES, SDATA, IERR  )
      use iso_c_binding, only: c_char
      use tconnect_mod
      implicit none

      character(kind=c_char) ::  zcabbrev(*)

      INTEGER       IDRUN
      INTEGER       IFCN
      INTEGER       IMDS
      INTEGER	    MAXTIMES
      INTEGER       IERR

      REAL	    TIMES(MAXTIMES)
      REAL          SDATA(MAXTIMES)

      CHARACTER*32 :: ZABBREV
      integer :: ilabb

      call cstring(ZABBREV,ZCABBREV, '2F')

      ZABBREV = ADJUSTL(ZABBREV)
      ilabb = len_trim(ZABBREV)
      if(ilabb.eq.0) then
         ZABBREV(1:1) = " "
         ilabb = 1
      endif

      CALL TRSCALAR_FETCH( IDRUN, ZABBREV(1:ilabb), IFCN, IMDS,
     1                           MAXTIMES, TIMES, SDATA, IERR  )
      END


      SUBROUTINE TRSCALAR_FETCH( IDRUN, ZABBREV, IFCN, IMDS,
     1                           MAXTIMES, TIMES, SDATA, IERR  )
 
      use tconnect_mod
 
C	SUBROUTINE TRSCALAR_FETCH READS SHOULD BE CALLED AFTER TRSCALAR_CONNECT TO
C                           READ IN THE DATA OF A SPECIFIED FUNCTION OF TIME AND
C			    ANOTHER COORDINATE AS WELL AS THE ASSOCIATED
C			    TIMES.
C
C  ** NOTE NOTE NOTE **
C
C       If there was no error in TRSCALAR_CONNECT, then IERR=0 or -99 entering
C       this subroutine.  The tree will be closed if IERR=0 or if there was an error
C       otherwise the tree will remain open.  With IERR=-99, TRSCALAR_CONNECT/TRSCALAR_FETCH 
C       should eventually be called with IERR=0 to close the tree.
C
C	ARGUMENTS:
C		   NOTE 1: UNIT LUN (91) IS USED TO READ IN TRANSP FILES.
C		   NOTE 2: UNIT LUOUT (6) IS USED TO OUTPUT ERROR MESSAGES.
C
C	LAST CHANGED:
C
C	   INPUT:  IDRUN    - INTEGER Index to run in CPLOTR from TRPROFIL_CONNECT
C	           ZABBREV  - C*Max 10 ABBREVIATION OF THE PROFILE TO READ IN.
C                  IFCN     - INTEGER INDEX OF FUNCTION
C                  IMDS     - INTEGER NONZERO IF THIS IS AN MDSPLUS CONNECTION
C		   MAXTIMES - INTEGER DIMENSION OF ARRAY TIMES.
C
C	   RETURNED:
C		   TIMES    - REAL ARRAY OF DIMENSION MAXTIMES INTO WHICH
C		              TIMES ARE READ FROM DISK.
C		   SDATA    - REAL ARRAY OF DIMENSION MAXTIMES INTO WHICH NTIMES
C			      SCALAR VALUES ARE READ.
C		   IERR     - INTEGER ERROR MESSAGE. RETURNED = 0 IF OK.
C			      0 = OK
C			      2 = ERROR READING NF FILE.
C			      3 = UNEXPECTED EOF READING NF FILE.
C			      5 = TRIED TO READ OLD FORMAT TF FILE.
C			     15 = SCALAR ABBREVIATION NOT FOUND
C			     16 = COULDN'T OPEN NF FILE
C			     17 = COULDN'T OPEN TF FILE
C                            21 = TRSCALAR: - FMT .NE. blank or F88
C			    102 = # OF SCALARS > DIMENSION OF BUFFER.
C		                  NO DATA READ IN.
C			    101 = # OF TIMES IN  FILE > MAXTIMES.
C				  MAXTIMES TIMES WERE READ IN.
 
      IMPLICIT NONE
C
C
      Character*(*) Zabbrev
      CHARACTER*10  ABBREV

      INTEGER       IDRUN
      INTEGER       IFCN
      INTEGER       IMDS
      INTEGER	    MAXTIMES
      INTEGER       IERR

      REAL	    TIMES(MAXTIMES)
      REAL          SDATA(MAXTIMES)
 
      integer i,ierr0,ierrc
 
      INTEGER  LUOUT,lunzer
 
C     -----------------------------------------------------------------------
	
      IERR0 = IERR
C
      LUOUT = lunzer(0)
C
      Abbrev = Zabbrev
      Call TrCaps(Abbrev)     ! Make all capitals
 
      IERR  = 0     ! INITIALIZE
      TIMES = 0.
      SDATA = 0.
 
      if (idrun<0) then
         write(luout,'(a)') ' %trscalar_fetch: idrun<0, '
     &        //'trscalar_connect must have failed'
         ierr=1 ; return
      end if

      call tgetscal_fetch(luout,idrun,abbrev,ifcn,maxtimes,
     1     times,sdata,ierr)

      if(ierr.ne.0) then
         write(luout,9902) abbrev
 9902    format(' %trscalar_fetch:  unexpectedly failed to read '
     &        //'named scalar:  ',a)
      endif
C
C  If MDSplus: close secondary run    ! CAL 1/28/00
      if (imds/=0) then
         if((ierr.ne.0).or.(ierr0.ne.-99)) then
            write(luout,*) ' %close MDS tree for 2nd run'
            call tconnect_close(ierrc)
            if(ierrc.ne.0) then
               if (ierr==0) ierr=1
               write(luout,9903) 
 9903          format(' %trscalar_fetch:  failed to close MDS Tree ')
            endif
         endif
      endif
C
      return
      end

