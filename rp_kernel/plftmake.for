      subroutine plftmake(lunt,zabr,zlbl,zuns,iordr,zamax,
     >   ztime,zdata,inpts)

      use datmgr_mod
      use cplotr_mod

C
C  create the scalar function in memory
C  called from plftmk*.for in rplot_sub
C
C  input
      integer lunt                      ! i/o channel for messages
      character*(*) zabr                ! function abbreviation
      character*(*) zlbl                ! function label
      character*(*) zuns                ! function units
      integer iordr                     ! timebase not strict ascending =>
C                                       ! iordr.gt.0
      real zamax                        ! max(abs(zdata(i),i=1 to inpts))
      real ztime(inpts)                 ! the timebase
      real zdata(inpts)                 ! the data
C
      CHARACTER*5 ZABBL
C
      REAL*8, PARAMETER :: ZERO = 0.0d0
      REAL ZANSI(60),ztmin,ztmax
      integer :: ix0 = 0   ! define for xinter
C
      EXTERNAL XIDENT
C
C--------------------------------
C
      NFTX=NFTX+1
      IND=NFT+NFTX
      ABT(IND)=zabr
      call aordr_add(abt,iordrt,ind)
C
      LABELT(IND)=zlbl
      UNITST(IND)=zuns
C
      WRITE(ZABBL,'(''%T'',I3.3)') NFTX
      CALL RP_DMGALO(NTT,JLOC,6)
C
      ztmin = minval(ztime(1:inpts))
      ztmax = maxval(ztime(1:inpts))
C
      IPT=LOCD(JLOC)
      DO IT=1,NTT
        if(nlxtrap0.and.
     >        ((time(it).lt.ztmin).or.(time(it).gt.ztmax))) then
           zans = ZERO
        else IF(IORDR.EQ.0) THEN
          CALL XINTER(XIDENT,TIME(IT),ZTIME,INPTS,
     >        IX0,IX0P1,ZXI,ZXIC,IEX)
          ZANS=ZDATA(IX0)*ZXIC+ZDATA(IX0P1)*ZXI
        ELSE
C  GENERALIZED INTERPOLATION FOR NONMONOTONIC T AXIS
C  MULTIPLE ANSWERS POSSIBLE
          IANS=0
          ZTIMI=TIME(IT)
          IF(ZTIMI.LT.ZTIME(1)) THEN
            IANS=IANS+1
            ZANSI(IANS)=ZDATA(1)
          ENDIF
          IF(ZTIMI.GE.ZTIME(INPTS)) THEN
            IANS=IANS+1
            ZANSI(IANS)=ZDATA(INPTS)
          ENDIF
          DO IT2=2,INPTS
            IT2M1=IT2-1
            IF( ((ZTIME(IT2M1).LE.ZTIMI).AND.(ZTIMI.LT.ZTIME(IT2)))
     >         .OR.((ZTIME(IT2M1).GE.ZTIMI).AND.(ZTIMI.GT.ZTIME(IT2))))
     >         THEN
               IF(ZTIME(IT2).EQ.ZTIME(IT2M1)) THEN
                  IANS=IANS+1
                  ZANSI(IANS)=ZDATA(IT2M1)
                  IANS=IANS+1
                  ZANSI(IANS)=ZDATA(IT2)
               ELSE
                  IANS=IANS+1
                  ZF=(ZTIMI-ZTIME(IT2M1))/(ZTIME(IT2)-ZTIME(IT2M1))
                  ZANSI(IANS)=(1.0-ZF)*ZDATA(IT2M1)+ZF*ZDATA(IT2)
               ENDIF
            ENDIF
          ENDDO
C  AVERAGE THE RESULTS AND CHECK THE VARIANCE
          ZANS=ZANSI(1)
          ZANSMN=ZANS
          ZANSMX=ZANS
          DO II=2,IANS
            ZANS=ZANS+ZANSI(II)
            ZANSMN=AMIN1(ZANSMN,ZANSI(II))
            ZANSMX=AMAX1(ZANSMX,ZANSI(II))
          ENDDO
          ZANS=ZANS/IANS
C
          IF(ABS(ZANSMN-ZANSMX).GT.0.0001*ZAMAX) THEN
            WRITE(LUNT,9009) ZTIMI,ZANS,ZANSMN,ZANSMX
 9009 FORMAT(' %PLFTMAKE:  T=',1PE11.4,' NON-UNIQUE INTERPOLATION,'/
     >'  AVG=',1PE11.4,' IS USED; MIN=',1PE11.4,' MAX=',1PE11.4)
          ENDIF
        ENDIF
C
C  STORE THE ANSWER .........
        ITI=IPT+IT-1
        DATBUF(ITI)=ZANS
      ENDDO
C
      NWDS(JLOC)=NTT
      DMGLBL(JLOC)=ZABBL
C
      return
      end
