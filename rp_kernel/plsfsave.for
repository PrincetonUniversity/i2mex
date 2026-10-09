      subroutine plsfsave(zinput,zlbl,zuns,jtyp,iadr,ier)

      use datmgr_mod
      use cplotr_mod

C  save calculator buffer (datbuf(iadr...)) as labeled function
C
      character*(*) zinput              ! fcn id
      character*(*) zlbl                ! 32 char label
      character*(*) zuns                ! 16 char units
      integer jtyp                      ! X axis type
      integer iadr                      ! where data is currently stored
C
      integer ier                       ! output completion code, 0=OK
C
C-----------------------------------------------
C
      if(ier.eq.-88) then
         ier=0                          ! skip abbrev. check
      else
         ier=0
         CALL PLABCK(ZINPUT,IER)        ! Check is name exists: Ier=1
         IF(IER.EQ.1) GO TO 1000
      endif
C
      NFXT=NFXT+1
      IND=NFXT
C
      ITYPR(IND)=JTYP
C
      ABR(IND)=ZINPUT
C
      LABELR(IND)=ZLBL
      UNITSR(IND)=ZUNS
      call aordr_add(abr,iordrr,nfxt)
C
C  CREATE PERMANENT SPACE FOR USER'S DATA -- IDENTIFIED AS SUCH BY
C  MPRIO(..)=7 IN DATMGR COMMON
C
C  save at end of datbuf...
 
      call plsfadd(jtyp,iadr1,zinput)   ! get permanent storage address
C
C  copy data
C
      ISIZ=NTR*NZONEX(JTYP)
C
      DO 200 IL=1,ISIZ
         DATBUF(IADR1+IL-1)=DATBUF(IADR+IL-1)
 200  CONTINUE
C
 1000 continue
      return
      end
