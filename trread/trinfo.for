C   *************** START FILE  TRINFO.FOR  ;  GROUP  TRREAD_LIB ******
C
C
      SUBROUTINE TRINFO( DISK,
     1			     DIR,     RUNIDIN, ZABBREV, ITYPE,
     2                       label, units, irank, idims, zxnames, IERR)
 
      use tconnect_mod
 
C
C	SUBROUTINE TRINFO reads in the labeling information on a run
C            and then returns a type code for a particular data item
C            (scalar or profile function) in that run.
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
C	   INPUT:  DISK     - C*64 DISK NAME WHERE DATA IS LOCATED.
C			      IF DISK=' ', ASSUME NAME IS LOGICAL DISK "RUNDATA
C                             If Unix system, ignore DISK.
C
C		   DIR      - C*64 FILE DIRECTORY, E.G. 'TRANSP.TFTR.86
C			      NO '[' OR ']' ALLOWED IN DIRECTORY NAME.
C		   RUNIDIN  - C*(*) INPUT RUNID
C
C           ZABBREV  - C*(*max 10) ABBREVIATION OF THE SCALAR TO READ IN.
C		   MAXTIMES - INTEGER DIMENSION OF ARRAY TIMES.
C
C	   RETURNED:
C                  ITYPE    - 0 if there is an error or fcn not found
C                            -1 if a scalar function of time
C                            +N if a profile with X axis id #N
C
C		   LABEL    - LABEL FOR SCALAR (was C*32 can now be longer)
C		   UNITS    - UNITS OF SCALAR (was C*16 can now be longer)
C
C                  IRANK  --- dimensionality:  1 for f(t), 2 for f(x,t), etc.
C                  IDIMS  --- size(s) of x dimension(s) & time dimension
C                             idims(irank)=#times  #x's at idims(1:irank-1)
C                  ZXNAMES -- name(s) of x dimension(s), zxnames(1:irank-1)
C
C		   IERR     - INTEGER ERROR MESSAGE. RETURNED = 0 IF OK.
C			      0 = OK
C			      2 = ERROR READING NF FILE.
C			      3 = UNEXPECTED EOF READING NF FILE.
C			      5 = TRIED TO READ OLD FORMAT TF FILE.
C			     15 = SCALAR ABBREVIATION NOT FOUND
C			     16 = COULDN'T OPEN NF FILE
C			     17 = COULDN'T OPEN TF FILE
C                            21 = TRSCALAR: - FMT .NE. blank or F88
 
      IMPLICIT NONE
C
      external TFILRD		! link the label reader
C
      CHARACTER*(*) DIR
      CHARACTER*(*)  DISK
      character*150   diskl     ! or MDS+ path
      CHARACTER*(*) RUNIDIN
      CHARACTER*10   RUNID
      Character*(*) Zabbrev
      CHARACTER*10   ABBREV
C
C  output...
C
      integer       ITYPE
C
      character*(*) label
      character*(*) units
C
      integer irank
      integer idims(*)
      character(*) zxnames(*)           ! generally expect character*10
C
      INTEGER       IERR
C
      integer idrun,iflag
      integer i,ierr0
 
      INTEGER  LUN
      INTEGER  LUOUT,lunzer
      DATA LUN  /91/            ! UNIT NUMBER FOR INPUT NF & TF FILES.
 
C	-----------------------------------------------------------------------
	
      IERR0 = IERR
      LUOUT = lunzer(0)
C
      Abbrev = Zabbrev
      Call TrCaps(Abbrev)     ! Make all capitals
 
      IERR  = 0     ! INITIALIZE
      RUNID = RUNIDIN   ! SET C*(*) TO C*10
C
      DISKL=DISK
C
      iflag=0
      call tconnect(-luout, diskl, dir, runid, idrun, ierr)
      if(idrun.eq.0) then
         iflag=1                        ! had to open connection
         call tconnect(luout, diskl, dir, runid, idrun, ierr)
      endif
C
      if(ierr.ne.0) then
         write(luout,9901) diskl,dir,runid
 9901    format(' %trinfo:  failed to connect to run:'/
     1     '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
      else
C
         call trinfo0(luout, idrun,abbrev, itype,
     >      label, units, irank, idims, zxnames, IERR)
      endif
C
C
C  if we connected...
C  If MDSplus: close secondary run    ! CAL 1/28/00
C
      if(iflag.eq.1) then
         if (disk(1:4) .eq. 'MDS+') then
            if((ierr.ne.0).or.(ierr0.ne.-99)) then
               write(luout,*) ' %close MDSplus tree for 2nd run'
               call tconnect_close(ierr)
               if(ierr.ne.0) then
                  write(luout,9903) abbrev,diskl,dir,runid
 9903             format(' %trinfo:  failed to close MDS Tree ',a/
     1               '  disk:  ',a/'  dir:   ',a/'  runid: ',a)
               endif
            endif
         endif
      endif
C
      return
      end
C----------------------------------------------
      subroutine trinfo0(luout, idrun, abbrev, itype,
     >   label, units, irank, idims, zxnames, IERR)
C
C  lookup a function inside an auxilliary run
C
      use cplotr_mod
C
C  input:
      integer luout                     ! i/o unit
      integer idrun                     ! run index (set by tconnect)
      character*(*) abbrev              ! id of function to look for
C
C  output:
      integer itype                     ! type of function, if found
C
      character*(*) label
      character*(*) units
C
      integer irank
      integer idims(*)
      character(*) zxnames(*)           ! generally expect character*10
C
      integer ier                       ! =0 if found, =15 if not found
C
C  itype=0 if ier.ne.0 on exit
C  itype=-1 indicates a scalar
C  itype=+N indicates a type N profile (vs. x axis id #N).
C
      itype=0
C
C  check scalar functions
C
      ifuns=nft_x(idrun)
      do if=1,ifuns
         if(abbrev.eq.abt_x(if,idrun)) then
            label=labelt_x(if,idrun)
            units=unitst_x(if,idrun)
            irank=1
            idims(1)=ntt_x(idrun)
            itype=-1
            go to 1000
         endif
      enddo
C
      ifuns=nfxt_x(idrun)
      do if=1,ifuns
         if(abbrev.eq.abr_x(if,idrun)) then
            label=labelr_x(if,idrun)
            units=unitsr_x(if,idrun)
            irank=2
            idims(2)=ntr_x(idrun)
            itype=itypr_x(if,idrun)
            idims(1)=nzonex_x(itype,idrun)
            zxnames(1)=xndabb_x(itype,idrun)
            go to 1000
         endif
      enddo
C
      write(luout,
     >   '('' %trinfo:  no function named "'',a,''" in run.'')') abbrev
      ier=1
C
 1000 continue
      return
      end
