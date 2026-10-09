C--------------------------------------------------------------
C  PLOTR-- DATA MGR INTERFACE --  **DMGFOTX**
C
C  DMGFOT-- READ SCALAR FCNS OF TIME INTO DATA AREA
C  DMGFOTX-- extend scalar fcns of time block with user defined fcns
C  DMGXOT-- READ TIME-VARYING X AXIS INTO DATA AREA
C  DMGGEO-- READ TIME-VARYING GEOMETRY INFO INTO DATA AREA
C  DMGFXT-- READ FCN OF TIME + ADDL COORDINATE INTO DATA AREA
C  DMDLOC-- LOCATE DATA IN DATA AREA
C  DMPRIN-- PRINT OUT CONTENTS OF DATA AREA
C
C
      SUBROUTINE DMGFOTX(ICALL,IPT,IER)
C
C  CREATE EXTENDED SPACE FOR
C  USER DEFINED FUNCTIONS
C
C        ICALL =2 CALL FROM TPLT2D -- COLLECT USER DEFINED FCNS;
C    EXTEND ALLOCATION OF SCALAR F(T) DATA BLOCK IF NECESSARY
C	       =3 CALL FROM READSCAL - SAME AS 2.  TBT 10/6/89
C
C  ICALL = 1 & ICALL = 4 are illegal!  use dmgfot
C
C  IER OUTPUT  = COMPLETION CODE, 0 DENOTES SUCCESS
C
      use datmgr_mod
      use cplotr_mod
      implicit NONE
C
C-----------------------------
C
      integer, intent(in) :: icall
      integer, intent(out) :: ipt
      integer, intent(out) :: ier
C
C-----------------------------
C
      integer :: ifxtnd,isize,jtim,ift,iadr,iptt,idum2,indt,it
      integer :: iptw,iwrk1,iwrk2,idiff,itest
C
      integer :: lunzer
C
      CHARACTER*5 ZABB
C
      logical idlock
C
C----------------------------------------
C  dmc 26 Dec 2009 -- support "no_delete" memory management mode;
C    if prio=10 worspace area is expanded, copy "in the way" profile
C    datasets up and out of the way; expand buffer if necessary
C
C  dmc 17 Nov 1997 -- support concept of primary or secondary run
C    primary run (lrun_x = 0 in COMMON) -- RPLOT's current runid
C    secondary run (lrun_x .gt. 0 in COMMON) -- for TRPROFIL or TRSCALAR
C
C  REDONE DMC SEPT 1987
C   SCALAR F(T) DATA IS MADE PERMANENTLY RESIDENT IN DATMGR
C
      IER=0
C
      IFXTND=16			! ALLOCATION EXTENSION, NO. USER FCNS
C
      IF(ICALL.EQ.1 .or. icall .eq. 4) THEN
         write(lunzer(0),*) ' ?dmgfotx call error (ICALL) - use dmgfot!'
         call bad_exit
      ELSE
         if(lrun_x.ne.0) then
            write(lunzer(0),*)
     1           ' ?dmgfotx: CODE ERROR, lrun_x .ne. 0'
            call bad_exit
         endif
         CALL DMDLOC('%F(T)',JTIM,ISIZE,IPT)

         DO IFT=1,NFTX
            IADR=IPT+(NFT+IFT-1)*NTT
            WRITE(ZABB,'(''%T'',I3.3)') IFT
            CALL DMDLOC(ZABB,INDT,IDUM2,IPTT)
            DO IT=1,NTT
               DATBUF(IADR+IT-1)=DATBUF(IPTT+IT-1)
            ENDDO
            idlock=no_delete
            no_delete=.FALSE.
            CALL DMIDEL(INDT)   ! delete temporary copy
            no_delete=idlock
            dmglbl(indt) = '(deleted)'
         ENDDO
C
         NFT=NFT+NFTX
         itest = NTT*(NFT + IFXTND/2)
         IF(itest.GT.ISIZE) THEN
            WRITE(lunzer(0),2001)
 2001       FORMAT(/' % DMGFOTX - EXPANDING SCALAR DATA MEMORY AREA')
C
C  this will fail in case of a delete lock...
C
            idiff = itest + NTT*IFXTND - isize
            call dmg_wkxpand(idiff,ier)
            if(ier.ne.0) then
               call errmsg_exit(
     >              ' ?dmgfotx: unexpected dmg_wkxpand error!')
            endif
C
            NWDS(JTIM)=itest + NTT*IFXTND
C
         ENDIF
C
      ENDIF
C
      NFTX=0
      RETURN
C
      END
C----------------------------------------
      subroutine dmgfotx_ww(icall,ipt,ier,iwrk1,iwrk2)
C
C  supplement DMGFOTX call with lookup of (possibly modified) workspace
C  addresses
C
      use datmgr_mod
      implicit NONE
C
C-----------------------------
C  passed:
C
      integer, intent(in) :: icall  ! as in DMGFOTX call
      integer, intent(out) :: ipt   ! as in DMGFOTX call
      integer, intent(out) :: ier   ! as in DMGFOTX call

      integer, intent(out) :: iwrk1,iwrk2
C
C-----------------------------
C  local:
C
      integer :: ind1,ind2,isiz1,isiz2
C-----------------------------
C
      call dmgfotx(icall,ipt,ier)
C
C WORKSPACES
C
      CALL DMDLOC('%WRK1',IND1,ISIZ1,IWRK1)
      CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)
C
      return
      end


