C	****************  START FILE TRPROFIL.FOR  ;  GROUP TRREAD_LIB  *******
 
 
      SUBROUTINE TRPROFIL( DISK,
     1			     DIR,    RUNIDIN, ZABBREV,MAXTIMES, MAXDATA,
     2			     LABEL,  UNITS,   ITYPE,  INX,
     3                       NTIMES, TIMES,   DATA,   IERR     )
 
      use tconnect_mod
 
C	SUBROUTINE TRPROFIL READS IN A SPECIFIED FUNCTION OF TIME AND
C			    ANOTHER COORDINATE AS WELL AS THE ASSOCIATED
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
C	UPDATES:
C
C           02/15/00  CAL  If MDSplus, close tree here, since lrun_x
C                          is set to 0 in tconnect and idrun is local.
C           11/18/97  DMC  ** rewrite ** for NetCDF support and performance
C                          use the RPLOT buffering mechanism
C                          the user interface is unchanged.
C
c                 ..8/13/96  dmc  RUNID now 10 chars
C                   6/12/95  tbt  DISK is ignore if unix system.
C                   6/02/93  TbT  Added call to TrCaps so abbrev can be in
C                                 lower case. Changed arg to ZAbbrev.
C	                          Took out Write statements - check for ierr.
C		    5/15/92  TBT  Changed call to INITPL2 from INITPL so as
C				  not to call GFLIB routines for Mikklesen.
C		    3/11/92  TBT  Add capability to read old FORMAT via
C     			          TRGETTM & MFBLKI
C		    1/18/91  TBT  PUT INERROR RETURN FOR NTIMES>MAXTIMES
C				  & TOOK OUT OF CLOSE (DISP="SAVE")
C		    4/10/90  TBT  ADDED CALL TO TFILHD.
C	            4/09/90  TBT  PUT IN CHECK FOR INX = 0
C		    3/27/90  TBT  CHANGED RUNID TO RUNIDIN, C*6->C*8.
C
C	ARGUMENTS:
C		   NOTE 1: UNIT LUN (77) IS USED TO READ IN THE TRANSP FILES.
C		   NOTE 2: UNIT LUOUT (6) IS USED TO OUTPUT ERROR MESSAGES.
C
C
C	   INPUT:  DISK     - C*64 DISK WHERE DATA RESIDES.
C                             IF = ' ', ASSUME LOGICAL DISK "RUNDATA"
C                             If Unix system - DISK is ignored.
C		   DIR      - C*64 FILE DIRECTORY, E.G. 'TRANSP.TFTR.86'
C			      NO '[' OR ']' IS ALLOWED IN DIRECTORY.
C		   RUNIDIN  - C*10 VARIABLE DEFINING THE RUNID, E.G. '4058  '.
C			      DIR & RUNID WILL BE USED TO MAKE TRANSP FILE NAMES
C	                       E.G. 'RUNDATA:[TRANSP.TFTR.86]4058MF.PLN    '
C			            'RUNDATA:[TRANSP.TFTR.86]4058TF.PLN    '
C			            'RUNDATA:[TRANSP.TFTR.86]4058NF.PLN    '
C	           ZABBREV  - C*Max 10 ABBREVIATION OF THE PROFILE TO READ IN.
C		   MAXTIMES - INTEGER DIMENSION OF ARRAY TIMES.
C		   MAXDATA  - INTEGER DIMENSION OF ARRAY DATA
C
C	   RETURNED:
C		   LABEL    - LABEL FOR PROFILE (was C*32 can now be longer)
C		   UNITS    - UNITS OF PROFILE (was C*16 can now be longer)
C		   ITYPE    - INTEGER TYPE OF PROFILE.
C		   INX      - # OF POINTS IN THE X DIRECTION
C                  NTIMES   - INTEGER # OF TIMES READ FROM DISK.
C		   TIMES    - REAL ARRAY OF DIMENSION MAXTIMES INTO WHICH
C		              TIMES ARE READ FROM DISK.
C		   DATA     - REAL ARRAY OF DIMENSION MAXDATA INTO WHICH
C			      NTIMES * INX PROFILE
C			      VALUES ARE READ.(( DATA(J,I) J=1,INX), I=1,NTIMES)
C		   IERR     - INTEGER ERROR MESSAGE. RETURNED = 0 IF OK.
C			      0 = OK
C			      2 = ERROR READING MF FILE.
C			      3 = UNEXPECTED EOF READING MF FILE.
C			      5 = TRIED TO READ OLD FORMAT TF FILE.
C			     15 = PROFILE ABBREVIATION NOT FOUND
C			     16 = COULDN'T OPEN MF FILE
C			     17 = COULDN'T OPEN TF FILE
C                            26 = RUNID IN FILE IS DIFFERENT THAN ASKED FOR
C			    101 = # OF TIMES IN  FILE > MAXTIMES.
C				  MAXTIMES TIMES WERE READ IN.
C			    102 = # OF PROFILES > DIMENSION OF BUFFER.
C		                  NO DATA READ IN.
c			    103 = TIMES RETURNED ARE NOT MONOTONICALLY
C				  INCREASING. DATA IS ALSO RETURNED.
C			    104 = # OF TIMES * # IN X DIRECTIONS > MAXDATA
C			    105 = # OF POINTS IN X DIRECTION = 0
C
C	SEE SUBROUTINE TFILRD FOR ORIGINAL VERSION.....
C
C
      IMPLICIT NONE
C
      external TFILRD		! link the label reader
      external TRFUNID          ! link this one too
C
      CHARACTER*(*)  DISK
      character*150   diskl     ! or MDS+ path
      CHARACTER*(*)  DIR
      CHARACTER*10   RUNID
      CHARACTER*(*) RUNIDIN
        Character*(*) ZAbbrev
      CHARACTER*10   ABBREV
      CHARACTER*(*)  LABEL
      CHARACTER*(*)  UNITS
	
      INTEGER       IERR
      INTEGER	      ITYPE
      INTEGER       INX
      INTEGER	      MAXTIMES
      INTEGER       MAXDATA
      INTEGER       NTIMES
	
      REAL	      TIMES(MAXTIMES)
      REAL          DATA (MAXDATA)
	
      INTEGER  LUOUT,lunzer
 
      integer i
      integer idrun,ierr0
 
C	-----------------------------------------------------------------------
 
      IERR0 = ierr
C
      LUOUT = lunzer(0)
C
      Abbrev = Zabbrev
      Call TrCaps(Abbrev)
	
      RUNID = RUNIDIN              ! MAKE CHARACTER*8
      IERR  = 0     ! INITIALIZE
C
      DISKL=DISK
C
      LABEL = ' '
      UNITS = ' '
      ITYPE = 0
      INX   = 0
      NTIMES= 0
 
      DO I=1,MAXTIMES
          TIMES(I) = 0
      END DO   ! I
 
      DO I=1,MAXDATA
          DATA(I) = 0.
      END DO   ! I
 
      call tconnect(luout, diskl, dir, runid, idrun, ierr)
      if(ierr.ne.0) then
         write(luout,9901) diskl,dir,runid
 9901    format(' %trprofil:  failed to connect to run:'/
     1     '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
      else
C
         call tgetprof(luout, idrun,abbrev,maxtimes,maxdata,label,units,
     1      itype,inx,ntimes,times,data,ierr)
C
         if(ierr.ne.0) then
            write(luout,9902) abbrev,diskl,dir,runid
 9902       format(' %trprofil:  failed to read profile:  ',a/
     1         '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
         endif
      endif
C
C  If MDSplus: close secondary run    ! CAL 1/28/00
      if (disk(1:4) .eq. 'MDS+') then
         if((ierr.ne.0).or.(ierr0.ne.-99)) then
            write(luout,*) ' %close MDS tree for 2nd run'
            call tconnect_close(ierr)
            if(ierr.ne.0) then
               write(luout,9903) abbrev,diskl,dir,runid
 9903          format(' %trprofil:  failed to close MDS Tree ',a/
     1            '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
            endif
         endif
      endif
C
      return
      end

C
C ===============================================================
C
 
      SUBROUTINE TRPROFIL_CONNECT( ZDISK,
     1			     ZDIR, ZRUNID, ZABBREV,
     2			     IDRUN, ZLABEL,ZUNITS, IFCN, ITYPE,  INX, 
     3                       NTIMES, IMDS, IERR     )
 
      use iso_c_binding, only: c_char
      use tconnect_mod
 
C	SUBROUTINE TRPROFIL_CONNECT   CONNECTS AND READS IN THE DIMENSIONS OF
C                           A SPECIFIED FUNCTION OF TIME AND
C			    ANOTHER COORDINATE AS WELL AS THE ASSOCIATED
C			    LABEL AND UNITS FROM A TRANSP RUN.  IT IS EXPECTED
C                           THAT TRPROFIL_FETCH WILL BE CALLED AFTER THIS
C                           FUNCTION IF IERR==0 OR IERR==-99.
C
C  ** NOTE NOTE NOTE **
C
C     Set IERR=0 or -99 on entering this subroutine.  If there is no error
C     than IERR will not change and the tree will remain open for the 
C     subsequent call to TRPROFIL_FETCH.  TRPROFIL_FETCH will then close the
C     tree if IERR=0 or leave it open if IERR=-99.
C
C	ARGUMENTS:
C		   NOTE 1: UNIT LUN (77) IS USED TO READ IN THE TRANSP FILES.
C		   NOTE 2: UNIT LUOUT (6) IS USED TO OUTPUT ERROR MESSAGES.
C
C
C	   INPUT:  DISK     - C*64 DISK WHERE DATA RESIDES.
C                             IF = ' ', ASSUME LOGICAL DISK "RUNDATA"
C                             If Unix system - DISK is ignored.
C		   DIR      - C*64 FILE DIRECTORY, E.G. 'TRANSP.TFTR.86'
C			      NO '[' OR ']' IS ALLOWED IN DIRECTORY.
C		   RUNID    - C*10 VARIABLE DEFINING THE RUNID, E.G. '4058  '.
C			      DIR & RUNID WILL BE USED TO MAKE TRANSP FILE NAMES
C	                       E.G. 'RUNDATA:[TRANSP.TFTR.86]4058MF.PLN    '
C			            'RUNDATA:[TRANSP.TFTR.86]4058TF.PLN    '
C			            'RUNDATA:[TRANSP.TFTR.86]4058NF.PLN    '
C	            ABBREV  - C*Max 10 ABBREVIATION OF THE PROFILE TO READ IN.
C
C	   RETURNED:
C                  IDRUN    - INTEGER Index of run in CPLOTR
C		   LABEL    - LABEL FOR PROFILE (was C*32 can now be longer)
C		   UNITS    - UNITS OF PROFILE (was C*16 can now be longer)
C                  IFCN     - INTEGER index of function
C		   ITYPE    - INTEGER TYPE OF PROFILE.
C		   INX      - INTEGER # OF POINTS IN THE X DIRECTION
C                  NTIMES   - INTEGER # OF TIMES READ FROM DISK.
C                  IMDS     - INTEGER NONZERO IF THIS IS AN MDSPLUS CONNECTION
C		   IERR     - INTEGER ERROR MESSAGE. RETURNED = 0 IF OK.
C		       0 or -99 = OK
C			      2 = ERROR READING MF FILE.
C			      3 = UNEXPECTED EOF READING MF FILE.
C			      5 = TRIED TO READ OLD FORMAT TF FILE.
C			     15 = PROFILE ABBREVIATION NOT FOUND
C			     16 = COULDN'T OPEN MF FILE
C			     17 = COULDN'T OPEN TF FILE
C                            26 = RUNID IN FILE IS DIFFERENT THAN ASKED FOR
C			    101 = # OF TIMES IN  FILE > MAXTIMES.
C				  MAXTIMES TIMES WERE READ IN.
C			    102 = # OF PROFILES > DIMENSION OF BUFFER.
C		                  NO DATA READ IN.
c			    103 = TIMES RETURNED ARE NOT MONOTONICALLY
C				  INCREASING. DATA IS ALSO RETURNED.
C			    104 = # OF TIMES * # IN X DIRECTIONS > MAXDATA
C			    105 = # OF POINTS IN X DIRECTION = 0
C
C	SEE SUBROUTINE TFILRD FOR ORIGINAL VERSION.....
C
C
      IMPLICIT NONE
C
      character(kind=c_char)  :: ZDISK(*)
      character*150 :: diskl     ! or MDS+ path
      character(kind=c_char)  :: ZDIR(*)
      CHARACTER*150 :: DIR
      character(kind=c_char)  :: ZRUNID(*)
      CHARACTER*10  :: RUNID
      character(kind=c_char)  :: ZAbbrev(*)
      CHARACTER*10  :: ABBREV
      character(kind=c_char)  :: ZLABEL(*)
      CHARACTER*128 :: LABEL
      character(kind=c_char)  :: ZUNITS(*)
      CHARACTER*64  :: UNITS
      
      INTEGER       IERR
      INTEGER       IDRUN
      INTEGER       IFCN
      INTEGER	    ITYPE
      INTEGER       INX
      INTEGER       NTIMES
      INTEGER       IMDS
		
      INTEGER  LUOUT,lunzer
 
      integer i
      integer ierr0,ierrc
 
C	-----------------------------------------------------------------------
 
      IERR0 = ierr
C
      LUOUT = lunzer(0)
C
      call cstring(DISKL,  ZDISK,   '2F')
      call cstring(DIR,    ZDIR,    '2F')
      call cstring(RUNID,  ZRUNID,  '2F')
      call cstring(ABBREV, ZABBREV, '2F')
      
      Call TrCaps(Abbrev)
	
c$$$      print *, 'DISKL  = ',trim(DISKL)
c$$$      print *, 'DIR    = ',trim(DIR)
c$$$      print *, 'RUNID  = ',trim(RUNID)
c$$$      print *, 'ABBREV = ',trim(ABBREV)

      IERR  = 0     ! INITIALIZE
C
      if (diskl(1:4) .eq. 'MDS+') then
         IMDS=1
      else
         IMDS=0
      end if

      LABEL = ' '
      UNITS = ' '
      IFCN  = 0
      ITYPE = 0
      INX   = 0
      NTIMES= 0
  
      call tconnect(luout, diskl, dir, runid, idrun, ierr)
      if(ierr.ne.0) then
         idrun=-1
         write(luout,9901) diskl,dir,runid
 9901    format(' %trprofil_connect:  failed to connect to run:'/
     1     '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
      else
C
         call tgetprof_connect(luout,idrun,abbrev,label,units,
     1      ifcn,itype,inx,ntimes,ierr)
C
         if(ierr.ne.0) then
            idrun=-1
            write(luout,9902) abbrev,diskl,dir,runid
 9902       format(' %trprofil:  failed to read named profile:  ',a/
     1           '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
         else
            call fstring2c(LABEL, ZLABEL)
            call fstring2c(UNITS, ZUNITS)    
         endif
      endif
C
C  If MDSplus: close secondary run    ! CAL 1/28/00
      if (imds/=0) then
         if((ierr.ne.0)) then
            write(luout,*) ' %close MDS tree for 2nd run'
            call tconnect_close(ierrc)
            if(ierrc.ne.0) then
               write(luout,9903) abbrev,diskl,dir,runid
 9903          format(' %trprofil_fetch:  failed to close MDS Tree ',a/
     1            '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
            endif
         endif
      endif

      if (ierr0==-99 .and. ierr==0) ierr=-99   ! restore for trprofil_fetch
C
c$$$      print *, 'IDRUN = ', IDRUN
c$$$      print *, 'LABEL = ', trim(LABEL)
c$$$      print *, 'UNITS = ', trim(UNITS)
c$$$      print *, 'IFCN  = ', IFCN
c$$$      print *, 'ITYPE = ', ITYPE
c$$$      print *, 'INX   = ', INX
c$$$      print *, 'NTIMES = ', NTIMES
c$$$      print *, 'IMDS   = ', IMDS
c$$$      print *, 'IERR   = ', IERR

      return
      
      contains

C
C     move a fortran string to a C string.  It is assumed the cstring is null terminated.  The fortran string
C     will only copy up to the null of the C string otherwise it will add a new null at the end of the copy.
C      
      subroutine fstring2c(fstr, cstr)
      use iso_c_binding, only: c_char, c_null_char
      character*(*) :: fstr
      character(kind=c_char)  :: cstr(*)

      integer :: i       ! temp
      integer :: nfstr   ! trimmed length of input string

      nfstr = len_trim(fstr)

      do i=1, nfstr
         if (cstr(i) .eq. c_null_char) return   ! no more room in C string
         cstr(i) = fstr(i:i)
      end do
      cstr(nfstr+1) = c_null_char  ! whole fortran string was copied so null terminate
      end subroutine fstring2c
      
      end
C
C
C ==============================================================
C
      SUBROUTINE TRPROFIL_FETCH_CSTRING( IDRUN, ZCABBREV, IFCN, IMDS,
     3                           MAXTIMES, MAXDATA, TIMES, DATA, IERR )
 
      use iso_c_binding, only: c_char
      use tconnect_mod
      implicit none

      character(kind=c_char) ::  zcabbrev(*)
	
      INTEGER       IDRUN
      INTEGER       IFCN
      INTEGER       IMDS
      INTEGER	    MAXTIMES
      INTEGER       MAXDATA
      INTEGER       IERR
	
      REAL	    TIMES(MAXTIMES)
      REAL          DATA (MAXDATA)

      CHARACTER*32 :: ZABBREV
      integer :: ilabb

      call cstring(ZABBREV,ZCABBREV, '2F')

      ZABBREV = ADJUSTL(ZABBREV)
      ilabb = len_trim(ZABBREV)
      if(ilabb.eq.0) then
         ZABBREV(1:1) = " "
         ilabb = 1
      endif

      CALL TRPROFIL_FETCH( IDRUN, ZABBREV(1:ilabb), IFCN, IMDS,
     3     MAXTIMES, MAXDATA, TIMES, DATA, IERR )
      END

      SUBROUTINE TRPROFIL_FETCH( IDRUN, ZABBREV, IFCN, IMDS,
     3                           MAXTIMES, MAXDATA, TIMES, DATA, IERR )
 
      use tconnect_mod
 
C	SUBROUTINE TRPROFIL_FETCH  SHOULD BE CALLED AFTER TRPROFIL_CONNECT TO
C                           READ IN THE DATA OF A SPECIFIED FUNCTION OF TIME AND
C			    ANOTHER COORDINATE AS WELL AS THE ASSOCIATED
C			    TIMES.
C
C  ** NOTE NOTE NOTE **
C
C       If there was no error in TRPROFIL_CONNECT, then IERR=0 or -99 entering
C       this subroutine.  The tree will be closed if IERR=0 or if there was an error
C       otherwise the tree will remain open.  With IERR=-99, TRPROFIL_CONNECT/TRPROFIL_FETCH 
C       should eventually be called with IERR=0 to close the tree.
C
C	ARGUMENTS:
C		   NOTE 1: UNIT LUN (77) IS USED TO READ IN THE TRANSP FILES.
C		   NOTE 2: UNIT LUOUT (6) IS USED TO OUTPUT ERROR MESSAGES.
C
C
C	   INPUT:  IDRUN    - INTEGER Index to run in CPLOTR from TRPROFIL_CONNECT
C	           ZABBREV  - C*Max 10 ABBREVIATION OF THE PROFILE TO READ IN.
C                  IFCN     - INTEGER INDEX OF FUNCTION
C                  IMDS     - INTEGER NONZERO IF THIS IS AN MDSPLUS CONNECTION
C		   MAXTIMES - INTEGER DIMENSION OF ARRAY TIMES.
C		   MAXDATA  - INTEGER DIMENSION OF ARRAY DATA
C
C	   RETURNED:
C		   TIMES    - REAL ARRAY OF DIMENSION MAXTIMES INTO WHICH
C		              TIMES ARE READ FROM DISK.
C		   DATA     - REAL ARRAY OF DIMENSION MAXDATA INTO WHICH
C			      NTIMES * INX PROFILE
C			      VALUES ARE READ.(( DATA(J,I) J=1,INX), I=1,NTIMES)
C		   IERR     - INTEGER ERROR MESSAGE. RETURNED = 0 IF OK.
C			      0 = OK
C			      2 = ERROR READING MF FILE.
C			      3 = UNEXPECTED EOF READING MF FILE.
C			      5 = TRIED TO READ OLD FORMAT TF FILE.
C			     15 = PROFILE ABBREVIATION NOT FOUND
C			     16 = COULDN'T OPEN MF FILE
C			     17 = COULDN'T OPEN TF FILE
C                            26 = RUNID IN FILE IS DIFFERENT THAN ASKED FOR
C			    101 = # OF TIMES IN  FILE > MAXTIMES.
C				  MAXTIMES TIMES WERE READ IN.
C			    102 = # OF PROFILES > DIMENSION OF BUFFER.
C		                  NO DATA READ IN.
c			    103 = TIMES RETURNED ARE NOT MONOTONICALLY
C				  INCREASING. DATA IS ALSO RETURNED.
C			    104 = # OF TIMES * # IN X DIRECTIONS > MAXDATA
C			    105 = # OF POINTS IN X DIRECTION = 0
C
C
      IMPLICIT NONE
C
C
      Character*(*) ZAbbrev
      CHARACTER*10  ABBREV
	
      INTEGER       IDRUN
      INTEGER       IFCN
      INTEGER       IMDS
      INTEGER	    MAXTIMES
      INTEGER       MAXDATA
      INTEGER       IERR
	
      REAL	    TIMES(MAXTIMES)
      REAL          DATA (MAXDATA)
	
      INTEGER  LUOUT,lunzer
 
      integer i
      integer ierr0,ierrc
 
C	-----------------------------------------------------------------------
 
      IERR0 = ierr
C
      LUOUT = lunzer(0)
C
      Abbrev = Zabbrev
      Call TrCaps(Abbrev)
	
      IERR  = 0     ! INITIALIZE
      TIMES = 0.
      DATA  = 0.

      if (idrun<0) then
         write(luout,'(a)') ' %trprofil_fetch: idrun<0, '
     &        //'trprofil_connect must have failed'
         ierr=1 ; return
      end if

      call tgetprof_fetch(luout,idrun,abbrev,ifcn,maxtimes,maxdata,
     1     times,data,ierr)
C
      if(ierr.ne.0) then
         write(luout,9902) abbrev
 9902    format(' %trprofil_fetch:  unexpectedly failed to read '
     &        //'named profile:  ',a)
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
 9903          format(' %trprofil_fetch:  failed to close MDS Tree ')
            endif
         endif
      endif
C
      return
      end
