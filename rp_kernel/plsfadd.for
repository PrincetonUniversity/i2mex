      subroutine plsfadd(jtyp,iadr1,zlbl)

      use datmgr_mod
      use cplotr_mod

C  create space for new user function of type "jtyp"
C  the space is created at iadr1
C
      integer jtyp                      ! fcn type (in)
      integer iadr1                     ! fcn data address (out)
      character*(*) zlbl                ! data id label (in)
C
C---------------------------------
C
 15   CONTINUE
C
      ILDEND=NDBSIZ+1
      J7=0
      JHIGH = 0
C
      DO 25 J=1,NDENT
        IF(MPRIO(J).EQ.7) THEN
          IF(LOCD(J).LT.ILDEND) THEN
            ILDEND=LOCD(J)
            J7=LPREV(J)
            JHIGH = J
          ENDIF
        ENDIF
 25   CONTINUE
C
      ISIZ=NTR*NZONEX(JTYP)
      INX = NZONEX(JTYP)
C
      IADR1=ILDEND-ISIZ
      IADR2=ILDEND-1
C
C  CLEAR OUT ANY LOWER PRIORITY STUFF THAT MIGHT BE IN THE WAY
c   (if not possible, get more space...)
c
      ILOC=1
 50   CONTINUE
      ILOC=LNEXT(ILOC)
 55   CONTINUE
      IF(DMGLBL(ILOC).EQ.'%FINI') GO TO 100
      IL1=LOCD(ILOC)
      IL2=IL1+NWDS(ILOC)-1
      IF((IL1.GT.IADR2).OR.(IL2.LT.IADR1)) GO TO 50
      ILOCN=LNEXT(ILOC)
      IF((MPRIO(ILOC).GT.5).or.NO_DELETE) THEN
         write(lunzer(0),*) 
     >        ' %plsfadd: expand buffer to make room for user fcn.'
         call dmg_datbuf_expand(0)
         go to 15  ! rescan...
      ELSE
         ! delete to make room...
         CALL DMIDEL(ILOC)
      ENDIF
      ILOC=ILOCN
      GO TO 55
C
 100  CONTINUE
C
C  IF THIS IS FIRST ONE, INSERT AT END
C
      IF(J7.EQ.0) THEN
          J7=LPREV(ILOC)
      ELSE
          J7 = LPREV(JHIGH)
      END IF  ! J7
C
      CALL DMINEW(J7,IND2,7)
      LOCD(IND2)=IADR1
      NWDS(IND2)=ISIZ
      DMGLBL(IND2)=ZLBL
C
      return
      end
