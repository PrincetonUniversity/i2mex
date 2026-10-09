C---------------------------------------------------------------------
C  INIWRK - COMPUTE SIZE AND INITIATE RPLOT WORKSPACES IN MEMORY
C
      SUBROUTINE INIWRK(ISMIN)
C
C  DMC JUNE 1989 -- PASSED ARGUMENT ISMIN SPECIFIES MINIMUM WORKSPACE
C    SIZE
C
C  DMC 1 SEPT 1987 -- BROKEN OUT OF INIRUN.FOR
C
C  CALLED AFTER SCALAR DATA IS READ IN.  FIRST THREE ITEMS IN
C  DATA MEMORY (DATBUF) SHOULD BE:
C
C  "%F(T)"  SCALAR DATA
C  "%WRK1"  1ST WORKSPACE
C  "%WRK2"  2ND WORKSPACE
C
C  ON EXIT FROM THIS ROUTINE.
C
C  THESE REGIONS OF DATBUF ARE ALLOCATED WITH "PRIORITY 10" TO
C  PROHIBIT DISPLACEMENT BY DATA READ IN RANDOMLY FROM PROFILE DATA-
C  BASE (WHICH I/O TAKES PLACE AT "PRIORITY 5"; LEAST RECENTLY
C  REFERENCED ITEMS AT EQUAL OR LOWER PRIORITY MAY BE OVERWRITTEN)
C  SEE DATMGR
C
C  Note DMC 27 Dec 2009: delete lock failure prevented, when called from
C  dmgfotx, by prior copying of prio=5 datasets "out of the way", to allow
C  for expansion of %F(T) and workspace prio=10 section at low end of DATBUF.

      use datmgr_mod
      use cplotr_mod
      implicit NONE

C---------------------
C
      LOGICAL MOMRUN
C
      integer :: ismin,iptt,jsizt,jtim,igeo,idum0,ix
      integer :: isizt,isizxt,ixmax1,ixmax2,isize,isize1,isize2
      integer :: iadr1,iadr2,iloc,il1,il2,ilocn,jloc
C
      integer :: lunzer
C
C---------------------
C
      integer,  save :: ismini = 0
C
C---------------------
C
      IF(ISMIN.GT.0) ISMINI=ISMIN
      ISMINI=MAX(16384,ISMINI)
C
C---------------------
      CALL DMDLOC('%F(T)',JTIM,JSIZT,IPTT)
      IF(IPTT.EQ.0) THEN
         write(lunzer(0),*)
     >        '?RPLOT/INIWRK - CODE ERROR, SCALAR DATA NOT IN MEMORY!'
         call bad_exit
      ENDIF
C
      MPRIO(JTIM)=10
C
      ISIZT=16*NTT
      isizxt=nzonex(1)*ntr
C
      IXMAX1=2+NAXFIX
      IXMAX2=2+2*NAXFIX
      DO 20 IX=1,NXR
         IXMAX1=MAX0(IXMAX1,NZONEX(IX))
         IXMAX2=MAX0(IXMAX2,NZONEX(IX))
 20   CONTINUE
C
      ISIZE1=MAX0(ISIZT,(MAX0(16,NTR)*IXMAX1))
      ISIZE2=MAX0(ISIZT,(MAX0(16,NTR)*IXMAX2))
C
C  IF RUN HAS MOMENTS EQUILIBRIUM DATA, WORKSPACES ARE USED
C  TO GENERATE EQUILIBRIUM SURFACE CONTOURS; ALLOCATE ENUF SPACE
C
      IF(MOMRUN(IDUM0,IGEO)) THEN
         ISIZE1=MAX0(ISIZE1,(NAXMMP*MAX0(NTR,IXMAX1)))
         ISIZE2=MAX0(ISIZE2,(NAXMMP*MAX0(NTR,IXMAX1)))
      ENDIF
C
C-------------------------
      ISIZE=MAX(ISIZE1,ISIZE2,ISMINI)
C-------------------------
C
C  ADDRESS RANGE OF TWO WORKSPACES
      IADR1=IPTT+JSIZT
      IADR2=IADR1+2*(ISIZE+isizxt)-1
C
C  DELETE ANYTHING IN THE WAY!
      ILOC=1
 50   CONTINUE
      ILOC=LNEXT(ILOC)
 55   CONTINUE
      IF(DMGLBL(ILOC).EQ.'%FINI') GO TO 100
      IL1=LOCD(ILOC)
      IL2=IL1+NWDS(ILOC)-1
      IF((IL1.GT.IADR2).OR.(IL2.LT.IADR1)) GO TO 50
      ILOCN=LNEXT(ILOC)
      CALL DMIDEL(ILOC)                 ! delete lock => failure !
      ILOC=ILOCN
      GO TO 55
C
 100  CONTINUE
C
C  WORKSPACES ALLOCATED WITH PRIORITY 10-- THE HIGHEST
C   ALL OTHER DATA ALLOCATED WITH PRIORITY 5-- AND MAY BE SWAPPED IN
C   AND OUT
C  ARGUMENTS GERRYMANDERED TO FORCE ALLOCATION OF SPACE CONTIGUOUSLY
C  AFTER %F(T)
C
C  FIRST WORKSPACE
      CALL DMGBSF(-ISIZE,JLOC,10)
C  RESERVE THE ALLOCATED SPACE FOR LATER USE
      NWDS(JLOC)=ISIZE
      DMGLBL(JLOC)='%WRK1'
C  SECOND WORKSPACE--
      CALL DMGBSF(-ISIZE,JLOC,10)
      NWDS(JLOC)=ISIZE
      DMGLBL(JLOC)='%WRK2'
C
C  x axes workspaces -- TRANSP runs only...
C  (mark as empty initially)
C
      if(NLTRANSP) then
         call dmgbsf(-isizxt,jloc,10)
         nwds(jloc)=isizxt
         dmglbl(jloc)='%X_XC'           ! becomes %XXXC when filled
C
         call dmgbsf(-isizxt,jloc,10)
         nwds(jloc)=isizxt
         dmglbl(jloc)='%X_XB'           ! becomes %XXXB when filled
C
         call dmgbsf(-isizxt,jloc,10)
         nwds(jloc)=isizxt
         dmglbl(jloc)='%X_ZC'           ! becomes %XAZC when filled
C
         call dmgbsf(-isizxt,jloc,10)
         nwds(jloc)=isizxt
         dmglbl(jloc)='%X_ZB'           ! becomes %XAZB when filled
      endif
C
      RETURN
      END
