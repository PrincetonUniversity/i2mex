!*************** START FILE CPLOTR.BLK ; GROUP CPLOTR ******************
!==============================================================
!  dmc -- NLTRANSP logical introduced (orig. by Rob Andre).  This
!    switch is used to support non-TRANSP databases being examined
!    by RPLOT.  TRANSP databases assumes the 1st two X axes for
!    profile functions are 1d "zone centers" X and "zone boundaries"
!    XB respectively; general databases do not have this assumption.
!    NLTRANSP = .TRUE. is set in rp_kernel/pconnect.for, iff the X/XB
!    assumption is validated in the run database.  cmt: Aug 2003.
!
!  dmc 6 Nov 2009 -- cleanup some unused quantities: NRLIM, NTLIM, NBLIM,
!     MAXFOT, MAXFXT, MAXMGP, MAXXVR, MAXXPT, MAXTPT (note that PARAMETERs
!     NAXFOT etc. are used instead); move NTIME-dimensioned arrays to 
!     datmgr_mod for dynamic allocatability...
!
!  dmc 12 Dec 2009 -- LABEL EXPANSION: 
!     description labels: 32chars -> 64chars  (LABELR, etc...)
!     units labels: 16chars -> 32chars        (UNITSR, etc...)
!
!  dmc -- f90 free form compatibility -- 27 Feb 2000
!    use "!" for all comments; & @73 / & @5 for continuations ...!
!
!  12/17/99 CAL -- Added definitions for MDSplus
!
!  13 Nov 1997 -- DMC -- Adding arrays to support storing in memory
!            the contents of several runs simultaneously.
!            Old arrays are defined as before; new arrays take
!            the old array names and append "_X", for storage of
!            information from "extra" runs.  So there is one "main"
!            run and up to several "extra" runs.
!
!  09/24/97 CAL -- Added definitions for netCDF
!
!  REV DMC MARCH 1990 -- EXPANDING RUN ID LABEL RUNID, RUNLB2
!
!--------------------------------------------------------------
!  REV D. MC CUNE  24 JULY 1981
!     STARTED ADDING VARIABLES TO GENERALIZE NON-TIME
!  COORDINATE OF PLOTTING VARIABLES TO HAVE NOT ONLY VARIABLE
!  LABELS/VALUES BUT VARIABLE NUMBER OF ZONES AS WELL.  ALSO
!  AN ASSOCIATION OF X-AXIS INDEX WITH A TIME-VARYING X-AXIS
!  VARIABLE FROM THE PLOT FILE WILL BE SUPPORTED
!  THE MEANING OF SOME OLD COMMON VARIABLES SHALL CHANGE.
!
!  PLOTR    D. MC CUNE    OCT 12 1979   PPPL
!        PLOTTING ROUTINE FOR TRANSP GRAPHICS OUTPUT
!           *VAX-11  VERSION #01**
!
!    *******    COMMON BLOCKS    ******
!
!    THIS FILE CONTAINS COMMON BLOCKS FOR TRANSP GRAPHICS
!     OUTPUT SYSTEM,
!
!    PROGRAMS:  RPLOT -- INTERACTIVE/BATCH TRANSP DATA ACCESS + PLOTTING
!              POPLOT -- POST-PROCESS ASCII TRANSP OUTPUT TO BINARY
!                        FORMAT READABLE BY RPLOT
!              PLABEL -- GENERATE LABEL/MAP FILE FOR TRANSP OUTPUT
!      plotgen python script -- GENERATE CODE USED BY PLABEL AND TRANSP TO
!                        CREATE LABELED TRANSP OUTPUT
!
!    FILES:
!     ####TF.PLN    LABELING AND FILE ORGANIZATION INFORMATION
!                   (SEQUENTIAL ACCESS ASCII)
!     ####NF.PLN    SCALAR FUNCTIONS OF TIME OUTPUT BY TRANSP
!                    (SEQUENTIAL ACCESS BINARY, CONTAINING TIMES
!                   AND A VARIABLE NUMBER OF SCALAR FUNCTIONS)
!     ####MF.PLN    PROFILE FUNCTIONS OF TIME AND ADDITIONAL COORDINATE
!                   OUTPUT BY TRANSP, REFORMATTED BY POPLOT.  (DIRECT
!                   ACCESS BINARY).
!     ####.CDF      netCDF file combining all 3 .PLN files
!
! [D. McCune 12 Dec 2009 -- this file is 30 years old!!]
! [We plan a NetCDF conversion s.t. PLN files and POPLOT post-processing
!  may no longer be needed...]
!
!=======================================================================
!  ALPHABETIZED DMC NOV 1988
!  INDEX OF COMMON VARIABLES:
!
!   ABB    10 CHARACTER ABBREVIATIONS FOR MULTIGRAPH PACKAGES
!   ABR    10 CHARACTER ABBREVIATIONS FOR PROFILE FUNCTIONS OF TIME
!          AND 1 ADDITIONAL COORDINATE
!   ABT    10 CHARACTER ABBREVIATIONS FOR SCALAR FCNS OF TIME
!
!   ADELT     IF KATFLG=1, AND THIS IS A MULTIPLOT VS. X (NOT TIME),
!             PLOT PROFILES AVERAGED OVER +/- ADELT SECONDS FROM THE
!             USER-SPECIFIED PLOTTING TIME
!
!  CDFILN  netCDF  FILE NAME
!------------------
! FOR VERY OLD TRANSP RUNS WITH TIME INDEPENDENT GEOMETRY, THE FOLLOWING
! 3 QUANTITIES ARE COMPUTED FROM RADIUS DATA READ IN FROM THE TF.PLN
! FILE...
!
!   DAREA  CROSS-SECTIONAL AREA OF ZONE  CM**2
!   DRAV   WIDTH OF ZONE  CM
!   DVOL   VOLUME OF ZONE  CM**3
!
! IN NEWER RUNS, SUCH GEOMETRIC INFORMATION IS COMPUTED IN TRANSP AND
! WRITTEN TO OUTPUT AS RPLOT DATA FUNCTIONS WITH NAMES SUCH AS 'DVOL'
!
! SEVERAL OTHER QUANTITIES FALL INTO THIS CATEGORY:  FORMERLY ASSUMED
! TIME INDEPENDENT, NOW CONSIDERED TIME DEPENDENT.  SEE NLTGEO LOGICAL
! SWITCH.
!------------------
!
!   FDIR   DIRECTORY WHERE RPLOT INPUT DATA (.PLN FILES) ARE LOCATED
!          OR BLANK IF THESE ARE IN THE CURRENT DEFAULT DIRECTORY
!          SEE ALSO LFDIR
!   FDISK  DISK WHERE RPLOT INPUT DATA FILES (.PLN FILES) ARE LOCATED
!          OR BLANK IF THESE ARE ON THE CURRENT DEFAULT DEVICE
!          SEE ALSO LFDISK
!
!   FILNC  FLAG WHETHER OR NOT FILENAMES ARE KNOWN-- POPLOT WORKS
!          BY EXTERNALLY ASSIGNING FILES TO LOGICAL UNIT NUMBERS
!
!   IDABS  PLOT SCALING CONTROLS -- ABSCISSAE
!   IDATE  CURRENT DATE (FORTRAN FORMAT)
!   IDCDF  netCDF File Id
!   IDORD  PLOT SCALING CONTROLS -- ORDINATE
!
!   IFUNB   LIST OF FUNCTION NUMBERS ASSOCIATED WITH EACH
!           MULTI-GRAPH.
!   IINTB   IINTB(IB)=0 MEANS MULTIGRAPH IB IS A PROFILE MULTIGRAPH
!                    =1 MEANS MULTIGRAPH IB IS A SCALAR MULTIGRAPH
!   INFB    INFB(J)= # OF CURVES ASSOCIATED WITH MULTI-GRAPH
!           PACKAGE NO. J
!
!   ISMDS   MDSplus Server socket
!
!   ITITLE  GENERIC PLOT LABEL
!
!   ITYPR   FCNS OF TIME AND RADIUS FLAG (ARRAY) INDICATING
!           WHETHER FCN IS DEFINED AT TRANSP ZONES OR AT
!           TRANSP BOUNDARIES
!           ... OR OTHER NON-TEMPORAL COORDINATE
!
!   KATFLG--  =1 ==> CALCULATE TIME AVERAGE FOR MULTIPLOTS
!                    VS. RADIUS; AVERAGE SEVERAL PROFILES
!                    TO GENERATE THE PROFILES GRAPHED
!             =0 ==> STANDARD OPERATION (JUST FETCH DATA AT
!                    REQUESTED TIME
!
!   LABELR  64 CHARACTER LABEL, FCNS OF TIME AND ADDL COORD. (ARRAY)
!   LABELT   DITTO, SCALAR FCNS OF TIME
!   LABELB  64 CHARACTER LABEL
!           FOR MULTI-GRAPH PACKAGES  (PACKAGES ENABLE THE
!           PLOTTING OF SEVERAL CURVES ON ONE SET OF AXES)
!
! THE FOLLOWING 3 ARE DEFINED IF NLTGEO=.TRUE., I.E. IF GEOMETRY
! IS TIME DEPENDENT:
!   LDAREA  FUNCTION NUMBER OF DATA FUNCTION CONTAINING D(AREA) VS
!           TIME AND FLUX ZONE INDEX
!   LDRAV   FUNCTION NUMBER OF DATA FUNCTION CONTAINING D(R)
!   LDRAVFAC  function no. of data containing <dr>*<1/dr> correction
!   LDVOL   FUNCTION NUMBER OF DATA FUNCTION CONTAINING D(VOLUME)
!   LELDOT  scalar function number, "ELDOT" data (for time derivatives)
!
!   LX      ptr to "xi" (zone ctrs)
!   LXB     ptr to "xi" (zone bdys)
!
!   LFDIR   (SEE FDIR) LENGTH OF DIRECTORY NAME WHERE RPLOT INPUT FILES
!           ARE TO BE FOUND; 0 FOR CURRENT DEFAULT DIRECTORY
!   LFDISK  (SEE FDISK) LENGTH OF DISK NAME WHERE RPLOT INPUT FILES
!           ARE TO BE FOUND; 0 FOR CURRENT DEFAULT DISK
!   LRUNID  LENGTH OF RUN ID (NO. OF CHARACTERS)
!
!   (SEE ALSO NLSENT= NUMBER OF INDEXED PLOTS SO FAR)...
!   (THIS IS TO SUPPORT THE RPLOT PLOT INDEXING FEATURE:  ALL RPLOT
!   PLOTS HAVE PAGE NUMBERS.  AT THE END OF AN RPLOT SESSION AN INDEX
!   CORRELATING PLOT LABELS WITH PAGE NUMBERS MAY BE GENERATED)
!
!    LSINDX(NLSIZ)-- ALPHABETIC-ORDERING ARRAY OF LABELS OF PLOTTED DATA
!           FOR GENERATING INDEX OF RPLOT GENERATED PLOTS FOR DISPLAY
!    LSLABL(NLSIZ)-- abbrev:label OF PLOTS MADE SO FAR (=NLSENT)
!    LSUNTS(NLSIZ)-- UNITS OF PLOTS MADE SO FAR (=NLSENT)
!    LSPAGI(NLSIZ)-- FIRST PAGE LABELED LSLABL(J) IS # LSPAGI(J)
!    LSPAGL(NLSIZ)-- PAGES LSPAGI(J) -- LSPAGI(J)+LSPAGL(J)-1 HAVE
!            LABEL LSLABL(J)
!
!    LSURF (IF NLTGEO IS .TRUE.) INDEX TO DATA FUNCTION CONTAINING
!           FLUX SURFACES VS FLUX SURFACE INDEX AND TIME
!
!    (moved to datmgr_mod):
!    LTWRIT  (POPLOT PROGRAM ONLY) INTERNAL COMMUNICATIONS FOR
!           DESELECTING OUT-OF-SEQUENCE TIME RECORDS
!           SEE ALSO NTCORR=NO. OF OUT-OF-SEQUENCE RECORDS TO DROP
!
!   MDSPULSE   MDSplus Pulse Number
!
!   MFILN   'MF' FILE NAME
!
!   NBAL NUMBER OF MULTI-GRAPH PACKAGES AVAILABLE
!
!   NFILN   'NF' FILE NAME
!
!   NFR  *USE IS HISTORICAL* NUMBER OF PROFILE FUNCTION "RECORDS"
!   NFT  NUMBER OF SCALAR FCNS OF TIME
!   NFT0  "ORIGINAL" NUMBER OF SCALAR FCNS OF TIME
!    DIFFERENCE BETWEEN NFT AND NFT0 INDICATES PRESENCE OF USER DEFINED
!    FUNCTIONS
!   NFTX  NUMBER USED IN GATHERING USER DEFINED SCALAR FUNCTIONS:
!    THESE ARE DEFINED OUTSIDE THE RPLOT SCALAR FUNCTION PLOTTING
!    CODE; EACH TIME THE SCALAR PLOTTING CODE IS ENTERED, ANY NEW
!    USER DEFINED FUNCTIONS ARE "GATHERED" AND NFTX IS RESET.  USER
!    DEFINED FUNCTIONS ARE GOTTEN OUTSIDE OF THE SCALAR CODE EITHER
!    BY READING A DATA UFILE OR BY EXTRACTING A SLICE FROM A PROFILE
!    FUNCTION, OR, WITH THE RPLOT CALCULATOR
!
!   NFX ... (PLTRGN) THE NUMBER OF FUNCTIONS DEFINED FOR EACH X
!    AXIS DEFINITION
!       ... (RPLOT) ID NO. OF DATA FUNCTION CONTAINING THE X
!    COORDINATE
!
!   NFXT ... TOTAL # OF FCNS OF TIME + 1 ADDL. COORD.
!   NFXT0 ... ORIGINAL TOTAL NFXT IN RUN RESULTS DATABASE.  IF NFXT
!    .GT. NFXT0, USER DEFINED PROFILE FUNCTIONS ARE PRESENT.  THESE
!    ARE READ FROM UFILES WRITTEN BY RPLOT E.G. FROM OTHER RELATED
!    TRANSP RUN DATABASES, OR, COMPUTED WITH THE RPLOT CALCULATOR
!
!------
!    NLCDF = .TRUE. if netCDF File
!------
!    NLMDS = .TRUE. if MDSplus
!------
!    NLSENT-- # OF DISTINCT-LABELED GRAPH GROUPS, SAVED SO THAT
!      AN INDEX OF PLOTS MADE CAN OPTIONALLY BE GENERATED AT THE END
!      OF THE RPLOT SESSION.  SEE LSPAGI, ET AL.
!
!-------
!   NLTGEO
!
!  SUPPORT TIME VARYING GEOMETRY
!  NLTGEO = .TRUE. TO INDICATE TIME VARYING GEOMETRY
!  DVOL,DAREA,DRAV,SURF-- DV,DA,DL, SURF AREAS TO USE IN GEOMETRIC
!    OPERATIONS IF NLTGEO=.FALSE.; IF NLTGEO=.TRUE. CORRESPONDING
!    QUANTITIES ARE READ FROM DATA FILE
!  IF NLTGEO=.TRUE. LDVOL,...,LSURF GIVE LOCATION IN MF FILE
!
!--------
!    NLTRANSP -- .TRUE. if this is a transp data set.  This is to support
!                other datasets such as forcebal.
!
!---
!  NLXFOT(J)=.TRUE. INDICATES THAT THE JTH NONTEMPORAL COORDINATE
!    IS INDEED A FCN OF TIME -- AND MUST BE READ FROM THE MF FILE
!     NLXFOT(J)=.FALSE. MEANS X AXIS IS CONSTANT, STORED IN XARRY
!  NLXFTD(J)=.TRUE. MEANS THE ORIGINAL (DEFAULT) X AXIS IS A FCN OF
!    TIME AND MUST BE READ FROM THE MF FILE
!           =.FALSE. MEANS THE ORIGINAL IS CONSTANT, STORED IN XARRY
!  NLXMON(J)=.TRUE. VERIFIES THAT JTH NONTEMPORAL COORDINATE IS
!    A MONOTONICALLY INCREASING FUNCTION AT ALL TIMES
!
!  THE X AXIS MAY BE INTERACTIVELY REDEFINED BY THE USER TO BE
!    ANY FCN DEFINED WITH RESPECT TO THE ORIGINAL (DEFAULT) X AXIS.
!
!  SEE XNDABB,XF,XFLAB,XFLABU ...
!
!  NLXVAR FLAGS POSSIBLE TIME VARYING GEOMETRY (CF FIXGEO ROUTINE)
!---
!  NOWNRS IS THE NUMBER OF PRIVELAGED "OWNERS" OF TRANSP RUNS
!  WHO MAY SAVE MODIFIED TF.PLN FILES
!
!    NPAGEG-- CURRENT PAGE # FOR PLOT OUTPUT -- REFERRED TO IN
!      RPLOT GENERATED DISPLAY INDEX
!
!     INTEGER NRECX(1...NXR)  # OF MF FILE RECORDS PER FUNCTION
!      DEFINED WITH RESPECT TO SPECIFIED X-AXIS ... SEE NZONEX
!      AND NFR ... PERTAINS TO "OLD" MF FILE FORMAT--
! ... DMC 1988 ... THE TRANSP PLOT DATA POSTPROCESSING PROGRAM
! CONVERTS TRANSP ASCII OUTPUT TO AN INTERMEDIATE FILE FORMAT
! A BINARY FILE WITH FIXED LENGTH RECORDS STORING ALL PROFILE DATA
! AT TIME T1, THEN TIME T2, ETC.  THIS IS THE **OLD** MF FILE FORMAT
! STILL IN THE ARCHIVES FOR OLDER TRANSP RUNS.  THEN ON A SECOND PASS
! THIS DATA IS CONVERTED TO A REFORMATTED BINARY FILE, THE **NEW** MF
! FILE FORMAT, WITH LARGER FIXED LENGTH RECORDS, ORGANIZED CONTAINING
! THE FILE RECORD TIMES T1, T2, ... TN, THEN THE COMPLETE TIME-X
! VARIATION OF PROFILE FUNCTION F1, THEN F2, ... UP TO NFXT0'TH F
! PROFILE FUNCTION IN THE FILE, THE RUN RESULT DATABASE.
!  RPLOT WILL WORK FOR BOTH OLD AND NEW FORMAT MF FILES BUT THE
!  NEW FORMAT YIELDS FASTER ACCESS TIMES.  SEE ALSO VARIABLES NFR
!  AND NZONEX WHICH ARE USED ESPECIALLY IN ACCESSING THE OLD FORMAT
!  MF FILE (STILL USED TODAY IN THE PLOT POSTPROCESSING PROGRAM)
!
!    NR0 --- MAX NO. OF PTS IN NON-TEMPORAL INDEPENDENT COORDINATE
!---
!    NROFFF -- ADDRESSING IN OLD MF FILE FORMAT
!  TO ACCOMODATE VARIABLE-SIZE X-AXIS, THE OFFSET OF RECORDS
!  FOR A GIVEN FCN OF (X,T) IS NO LONGER GIVEN BY THE FCN NUMBER.
!  SO THE FOLLOWING ARRAYS ARE PROVIDED TO SPECIFY THE OFFSET
!   NROFFF(J)-- OFFSET FOR FIRST RECORD FOR FCN J
!     THE JTH FCN'S VARIATION WITH RESPECT TO ITS INDEPENDANT
!    NONTEMPORAL COORDINATE, AT THE KTH TIME, IS LOCATED IN
!    THE MF FILE, STARTING AT RANDOM ACCESS RECORD
!    IREC1=(K-1)*(NRF+1) + NROFFF(J)
!    AND CONTINUING TO
!    IREC2=IREC1+NRECX(ITYPR(J)) -1
!  THE NUMBER OF RECORDS BEING CHARACTERISTIC OF THE NUMBER OF
!  PTS IN THE NONTEMPORAL COORDINATE
!    NROFFX -- SIMILAR ADDRESSING FOR X AXIS DATA
!---
!   NRUN    TRANSP RUN NUMBER ** SUPERSEDED BY 6 CHAR RUNID **
!
!   NSCALE  NO. OF SCALING DEFAULT RULES SAVED
!   NSCALC  SCALING DEFAULTS SWITCH ARRAY
!
!   NSHOT   TRANSP RUN DATA SHOT NUMBER
!
!      MAX NO. OF TIME PTS. -- moved: see datmgr_mod "NTIME", dynamic now
!
!   NTR  NUMBER OF TIME PTS, FCNS OF RADIUS AND TIME
!   NTT  NUMBER OF TIME PTS, FCNS OF TIME
!
!   NXR  NUMBER OF KNOWN NON-TEMPORAL INDEPENDENT COORDINATES
!
!   NZONES  NUMBER OF RADIAL ZONES
!        *7-24-81 NUMBER OF WORDS IN MF RANDOM ACCESS FILE
!        RECORDS (RECORDSIZE) FOR OLD MF FILE FORMAT
!   NZONEX
!     INTEGER NZONEX(1...NXR) # OF PTS IN FUNCTION FOR SPEC.
!      X AXIS.  NOTE NRECX=(NZONEX-1)/NFR  + 1
!
!   OWNER
!   OWNERS
!   OWNUIC  SYMBOLS USED TO DETERMINE IF CURRENT USER OF RPLOT
!     IS "PRIVELEGED" TO HAVE WRITE ACCESS TO THE TF.PLN FILE
!     E.G. TO PERMANENTLY DEFINE A NEW MULTIGRAPH ASSOCIATION OR
!     LABEL CHANGE
!
!  PFLAB, PFLABU USED TO SET THE CORRECT PLOT LABEL -- E.G. WITH UNITS
!    TRANSFORM DUE TO INTEGRAL OPERATORS -- FOR THE UPCOMING PLOT
!
!  PLTABB -- ABBREVIATION (RPLOT DATA/MULTIGRAPH ID) FOR UPCOMING PLOT
!
!----------------
!   SOME STUFF BEING PHASED OUT DUE TO INTRODUCTION OF TIME-VARYING
!    GEOMETRY.  SEE FIXGEO SUBROUTINE
!
!   RBOUN   TRANSP RADIAL ZONE BOUNDARIES (ARRAY)
!      *7-24-81  NOT USED
!   RMINOR  PLASMA MINOR RADIUS OF TRANSP RUN
!   RMAJOR  PLASMA MAJOR RADIUS
!   RZON    TRANSP RADIAL ZONE CENTERS (ARRAY)
!      *7-24-81  NOT USED
!----------------
!   RSCALC  PLOT SCALING DEFAULTS INFORMATION VECTOR
!
!   RUNID   MAX 10 CHAR RUN ID OF RUN BEING EXAMINED WITH RPLOT
!   RUNLB2  EXTENDED LABEL "TOK.YY NNNN" FOR PLOTS
!
!   SSELEC  SUBSTRING SELECTOR OPTIONALLY USED WHEN GENERATING LISTS
!     OF RPLOT NAMES
!
!   SURF    FLUX SURFACE SURFACE AREAS *TIME INDEPENDENT GEOMETRY RUNS*
!           SEE NLTGEO, AND DVOL, ETC.
!
!   SXMIN,SXMAX  PLOT SCALE DEFAULT DATA
!   SYMIN,SYMAX  MORE PLOT SCALE DEFAULT DATA
!
!   TFILN   'TF' FILE NAME
!
!  moved to datmgr_mod:
!   TIME     ARRAY FOR TIME VALUES-- FCNS OF TIME
!   TIME3    ARRAY FOR TIME VALUES-- FCNS OF TIME AND ADDL COORDINATE
!
!   TIMLAB   ASCII LABEL FOR TIME
!   TIMUNS   UNITS FOR TIME
!
!   UNITSB   PHYSICAL UNITS LABEL 32 CHARACTERS FOR MULTIGRAPHS
!   UNITSR   PHYSICAL UNITS LABEL 32 CHARS FOR PROFILE FUNCTIONS
!   UNITST   PHYSICAL UNITS LABEL 32 CHARS FOR SCALAR FUNCTIONS
!
!    XARRY(--,--)  VALUES OF INDEPENDANT VARIABLES FOR
!                  PLOTTING-- IF TIME INDEPENDENT
!
!  ARRAY XF AND LABEL ARRAY XFLAB, XFLABU ARE USED TO SPECIFY THE
!    (CURRENTLY IN USE) X AXIS TO PLOTTING ROUTINES,
!    XFABB IS THE ABBREVIATION DISPLAYED AS A PLOT LABEL
!    XFPABB IS ANOTHER ABBREVIATION PLOT LABEL
!
!    XLAB(--)  64 CHAR LABEL FOR EACH INDEPENDANT VARIABLE
!    XLABU(--) 32 CHAR PHYSICAL UNITS LABEL FOR EACH INDEP. COORD.
!
!    XLABD,XLABDU  ARE UTILITY LABELS FOR PASSING TO PLOT ROUTINES
!
!     CHARACTER*5 XNDABB(1...NXR) -- 10 CHARACTER ABBREVIATION
!       GIVING ASSOCIATION OF X AXIS TO FUNCTION IN MF FILE
!       WHICH DEFINES THE X AXIS IN PHYS. UNITS AS A FCN OF TIME
!
!  THE CHARACTER ARRAY XNDABB STORES THE ABBREVIATION ID OF THE
!    DEFAULT X AXIS
!
!-------------------------------------------------------------------
!
!  CPLOTR LIMITATIONS HAVE BEEN QUANTIZED IN BLOCK TRPLIM
!   VARIABLES MAXFOT ETC.  SEE BLOCK-DATA SUBROUTINE CPLSET *****
!     THESE MAXIMA ARE DETERMINED BY THE SIZES OF ARRAYS
!   DECLARED IN THIS FILE ****
!
!
!-------------------------------------------------------------------
!  DECLARATIONS
!  INDEX BUFFER SIZE
module cplotr_mod
  implicit none

  integer, PARAMETER :: NLSIZ = 4096
  !  PLOTTING BUFFER SIZE--
  integer, PARAMETER :: NR0=3*257*257 ! MAX NO. NON-TEMPORAL AXIS POINTS
                                    ! RGA, account for PSIRZ(R,Z) output

  !  NUMBER OF FCNS AND PACKAGES LIMITS
  integer, PARAMETER :: NAXFOT=4096   ! MAX NO. OF SCALAR FUNCTIONS
  integer, PARAMETER :: NAXFXT=4096   ! MAX NO. OF PROFILE FUNCTIONS
  integer, PARAMETER :: NAXMGP=2000   ! MAX NO. OF MULTIGRAPH GROUPINGS

  integer, PARAMETER :: NAXMGF=126    ! max no. items per multigraph*
  ! *display limit may be lower.  Note: some routines w/o CPLOTR
  ! declare this dimension with matching constant.  Check these
  ! libraries:  tfplot rplot_io -- the ascii I/O routines tfilrd.for,
  ! tfils7.f00 and tfilwr.for will need work, if this PARAMETER value expands.
  ! Also search on 'parameter.*naxmgf','CPLOTR.*NAXMGF' and look in plotgen.py
  !   cd $CODESYSDIR/source ; find . -maxdepth 2 -exec grep -H -i 'CPLOTR.*NAXMGF' {} \;

  integer, PARAMETER :: NAXXVR=50     ! MAX NO. OF NON-TIME COORDS DEFINED
  integer, PARAMETER :: NAXXVT=NAXXVR+2    ! MAX NO. OF COORDS INCL. TIME
  !  MAX NO OF MOMENTS AND NUMBER OF CONTOUR PTS, MOMENTS PLOTTING
  integer, PARAMETER :: NAXMOM=64
  integer, PARAMETER :: NAXMMP=301
  !  NO OF FIXED PTS FOR 3D PLOTS REQUIRING INTERPOLATED FIXED X AXIS
  integer, PARAMETER :: NAXFIX=81
  !---------------

  INTEGER nbal,nrun,nshot,nzones,nft0,nft,nftx,       &
       nfr,ntt,ntr,nownrs,lrunid,lfdisk,lfdir,lun_tf,lun_mf,lun_nf,     &
       nscale,katflg,nxr,nfxt,nfxt0,npageg,nlsent,                      &
       ldvol,ldarea,ldrav,lsurf,ldravfac,leldot,lx,lxb

  integer, dimension(NAXXVT)   :: IDORD
  integer, dimension(NAXXVT)   :: IDABS
  integer, dimension(2,NAXXVT) :: NSCALC  ! SCALING DEFAULTS SWITCH ARRAY
  integer, dimension(NAXXVT)   :: NAXISC

  REAL :: rminor,rmajor,adelt
  real, dimension(NAXXVT) :: sxmin, sxmax
  real, dimension(NAXXVT) :: symin, symax
  real, dimension(4,NAXXVT) :: rscalc
  !
  CHARACTER*64 XLAB(NAXXVR),XLABD(NAXXVR),XFLAB,PFLAB
  !
  CHARACTER*32 XLABU(NAXXVR),XLABDU(NAXXVR),XFLABU,PFLABU
  !
  CHARACTER*20 OWNERS(10),OWNUIC(2,5,10)
  !
  CHARACTER*10 XNDABB(NAXXVR),XFABB,XFPABB,PLTABB
  !
  CHARACTER*10 RUNID
  CHARACTER*9  IDATE
  CHARACTER*60 ITITLE
  !
  CHARACTER*64 TIMLAB
  CHARACTER*32 TIMUNS
  CHARACTER*140 TFILN,MFILN,NFILN,FDISK,FDIR,CDFILN,TMP_FILN
  CHARACTER*30 RUNLB2
  CHARACTER*20 OWNER
  CHARACTER*1 FILNC
  !
  character*140 curstring       ! composite string for "main" run
  !
  CHARACTER*64 LABELR(NAXFXT),LABELT(NAXFOT),LABELB(NAXMGP)
  CHARACTER*32 UNITSR(NAXFXT),UNITST(NAXFOT),UNITSB(NAXMGP)
  CHARACTER*21 ABR(NAXFXT),ABT(NAXFOT)
  real ttagr(naxfxt),ttagt(naxfot)	
  CHARACTER*10 ABB(NAXMGP)
  !
  CHARACTER*76 LSLABL(NLSIZ)  ! size = size(abbrev) + 2 + size(label)
  CHARACTER*32 LSUNTS(NLSIZ)
  !  FCN/MG NAME LISTING SELECTOR SUBSTRING
  CHARACTER*10 SSELEC,SSELEC_SAVE
  !
  !  new dmc 5 Aug 1999 -- sorting arrays
  !
  integer iordrt(naxfot),iordrr(naxfxt),iordrb(naxmgp)
  !
  INTEGER KSELEC      ! list output option control ...
  !
  REAL RZON(NR0),RBOUN(NR0)
  REAL XARRY(NR0,NAXXVR),XF(NR0),XPATCH(NAXXVR)
  REAL DVOL(NR0),DAREA(NR0),DRAV(NR0),SURF(NR0)
  !
  !! see datmgr_mod      REAL TIME(NTIME),TIME3(NTIME),WORK1(NTIME),WORK2(NTIME)
  !
  INTEGER IDCDF, ISMDS, MDSPULSE
  INTEGER NFX(NAXXVR),NRECX(NAXXVR),NZONEX(NAXXVR),NROFFX(NAXXVR),NROFFF(NAXFXT)
  !
  INTEGER ITYPR(NAXFXT),IINTB(NAXMGP),INFB(NAXMGP),IFUNB(NAXMGF,NAXMGP)
  !
  INTEGER LSPAGI(NLSIZ),LSPAGL(NLSIZ),LSINDX(NLSIZ)
  !
  LOGICAL NLCDF, NLMDS
  LOGICAL NLTGEO, NLTRANSP
  LOGICAL NLXFOT(NAXXVR),NLXMON(NAXXVR),NLXFTD(NAXXVR)
  LOGICAL NLXVAR
  logical nlxtrap0
  !
  !  new dmc 10 March 2000 -- units conversion factor sorting arrays
  !
  character*32 units_orig(naxfxt)
  character*32 units_mks(naxfxt)
  real conv_units(naxfxt)     ! *conv_units(..):  units_orig -> units_mks
  integer iordru(naxfxt),nu_mks
  !
  !---------------------------------
  !  dmc -- data for "extra" runs
  !
  integer, parameter :: naxxtra = 10
  !
  integer nrun_x                ! number of extra runs info stored
  integer krun_x                ! index of last extra run read in
  integer lrun_x          ! index of extra run currently in use
  !
  !  abbreviations, units, labels
  !
  character*10 abr_x(naxfxt,naxxtra),abt_x(naxfot,naxxtra)
  character*32 unitsr_x(naxfxt,naxxtra),unitst_x(naxfot,naxxtra)
  character*64 labelr_x(naxfxt,naxxtra),labelt_x(naxfot,naxxtra)
  !
  !  data disk,dir,runid
  !
  character*140 fdir_x(naxxtra),fdisk_x(naxxtra)
  character*10 runid_x(naxxtra)
  !
  !  mdsplus info
  !
  logical nlmds_x(naxxtra)
  integer ismds_x(naxxtra), mdspulse_x(naxxtra)
  !
  !  cdf info
  !
  logical nlcdf_x(naxxtra)
  integer idcdf_x(naxxtra)
  !
  !  timebase info
  !
  character*64 timlab_x(naxxtra)
  character*32 timuns_x(naxxtra)
  !! see datmgr_mod...      real time_x(ntime,naxxtra),time3_x(ntime,naxxtra)
  !
  !  runid labeling
  !
  character*40 rlbl(naxxtra)  ! labels
  integer lrlbl(naxxtra)      ! lengths of labels
  !
  !  # of time pts, scalar fcns, profile fcns
  !
  integer ntt_x(naxxtra),ntr_x(naxxtra)
  !
  !  more #s of things
  !
  integer nft_x(naxxtra)            ! no.s of f(t) fcns
  integer nfxt_x(naxxtra)           ! no.s of f(x,t) fcns
  integer nfx_x(naxxvr,naxxtra)     ! x axis fcn ptr
  integer nxr_x(naxxtra)            ! no.s of x-axis classes
  integer nzonex_x(naxxvr,naxxtra)  ! no.s of zones / x-axis
  integer itypr_x(naxfxt,naxxtra)   ! x-axis of each f(x,t)
  !
  !  for ancient style MF.PLN files
  !
  integer nfr_x(naxxtra)            ! no.'s of file records
  integer nrofff_x(naxfxt,naxxtra)  ! record offsets
  integer nrecx_x(naxxvr,naxxtra)   ! record counts
  integer nzones_x(naxxtra)         ! record sizes
  character*10 xndabb_x(naxxvr,naxxtra)  ! x axis ids
  !
  !  old NetCDF files contain 1 time axis, newer ones contain two time
  !  axes... keep track of what we have open currently...
  !
  integer ncdft
  integer ncdft_x(naxxtra)
  !
  !  mds+ data cache support items
  !
  logical mds_cache               ! TRUE if mds_cache option is active
  integer nbcache_lim             ! max binary items in cache dir.
  logical mds_cache_only          ! TRUE if mds_cache is the ONLY data
  !
  character*100 mds_cache_root    ! mds+ cache system, root directory
  character*200 mds_cache_dir     ! mds+ cache directory, current run
  character*200 mds_cache_dir_x(naxxtra) ! extra runs' cache directories
  !
  integer nbcache_act             ! actual #binary items in cache
  integer nbcache_act_x(naxxtra)  ! ditto, 2ndary runs
  !
  integer nbcache_list(naxfxt)    ! list of binary items in cache
  integer nbcache_list_x(naxfxt,naxxtra) ! ditto, 2ndary runs
  !
  integer nbcache_lun             ! lun for cache related i/o
  !
  !----------------------------------------------------------------
  !  COMMON stub for imbedding inside executing TRANSP run
  !
  real transp_tinit       ! transp start time
  integer transp_state(2) ! transp data state
  !       element(1) refers to scalar data
  !       element(2) refers to profile data
  !         meaning of values:
  !               0 = no data
  !               1 = one time point
  !               2 = two time points
  logical transp_imbed    ! .TRUE. if running imbedded
  !
  !----------------------------------------------------------------
  !  list of deferred-closed runs
  !
  integer, parameter :: mds_nclist_max = 100
  integer mds_nclist, mds_shot_clist(mds_nclist_max)
  character*40 mds_tree_clist(mds_nclist_max)
  character*150 mds_save_server

end module cplotr_mod
