C******************** START FILE SMTIMA.FOR ; GROUP SMOOF2 ******************
C------------------------------------------------------------
C  SMTIMA
C
C  TIME AVERAGE PROFILE DATA
C
C   TIME3(NTR) DATA PTS IN TIME
C
C   INX X PTS
C
C  DATA AT DATBUF(IPF)
C  WORKSPACE AT DATBUF(IPF2)
C
C  TIME AVG AT EACH TIME PT. +/- ZTDELT
C   IF IOPT=1 DO STRAIGHT TIME AVG; IF IOPT=2 DO DOUBLE INV. TIME AVG
C
      SUBROUTINE SMTIMA(IPF,INX,IPF2,ZTDELT,IOPT)
C
      use datmgr_mod
      use cplotr_mod
C
      if(inx.le.0) then
         call smtima2(ipf,1,ipf2,ztdelt,iopt,time,ntt)
      else
         call smtima2(ipf,inx,ipf2,ztdelt,iopt,time3,ntr)
      endif
      return
      end
C-----------------------------------------------------------------------
      subroutine smtima2(IPF,INX,IPF2,ZTDELT,IOPT,ztime,inumt)
C
      use datmgr_mod
C
      real ztime(inumt)
C
C  WORKSPACE CONTAINS MANY*INUMT DATA PTS (SEE SUBROUTINE INIRUN,
C  IN PLOTR1.FOR)
C
C  COPY FCN TO BE AVERAGED INTO WORKSPACE
      INTOT=INUMT*INX
      CALL copyr4(DATBUF(IPF),DATBUF(IPF2),INTOT)
C  USE FOR WEIGHTING FACTORS OF TIME AVG
C
      call smtima0(ztime,datbuf(ipf2),inumt,datbuf(ipf),inumt,
     >   inx,1,ztdelt,iopt)
C
      return
      end
C----------------------------------------------------------------------
      subroutine smtima0(ztime,zdata_in,inumt,zdata_out,intout,
     >   inx,itstart,ztdelt,iopt)
C
      use datmgr_mod
      use cplotr_mod
C
      real ztime(inumt)                 ! timebase
      real zdata_in(inx,inumt)          ! data input, to be averaged
      real zdata_out(inx,intout)        ! data output, time averaged
C
      integer itstart                   ! first time index for output
      real ztdelt                       ! averaging window +/-
      integer iopt                      ! regular or =2 for dbl inv averaging
C
C-------------------------------------------------
C
      real, dimension(:), allocatable :: zdts
C
C-------------------------------------------------
C
      allocate(zdts(ntime))
      zdts=0.0
C
      ZTMIN=ZTIME(1)-0.5*(ZTIME(2)-ZTIME(1))
      ZTMAX=ZTIME(INUMT)+0.5*(ZTIME(INUMT)-ZTIME(INUMT-1))
      ZT0=ZTMIN
      IT0=1
C
      zdtmin=1.0e-8*(ztime(inumt)-ztime(1))
      if(zdtmin.eq.0.0) zdtmin=1.0e-8
C
C  LOOP OVER TIME
      itend = itstart+intout-1
      DO 200 JT=itstart,itend
         itout=jt-itstart+1
C  INTERVAL
         ZT1=AMAX1(ZTMIN,(ZTIME(JT)-ZTDELT))
         ZT2=AMIN1(ZTMAX,(ZTIME(JT)+ZTDELT))
C  LOOP IN TIME TO LOCATE PTS WHICH CONTRIBUTE TO AVG AT TIME JT
         ZWSUM=0.0
         ZTP=ZT0
C
         IT1=0
         IT2=0
         DO 100 JT2=IT0,INUMT
C
            ZT=ZTMAX
            JT2P1=JT2+1
            IF(JT2.LT.INUMT) ZT=0.5*(ZTIME(JT2)+ZTIME(JT2P1))
C  CHECK FOR INTERSECTION
            IF(ZT.LE.ZT1) GO TO 90
C  FLAG INTERSECTION
            IF(IT1.EQ.0) THEN
               IT1=JT2
               ZT0=ZTP
            ENDIF
C  EXTENT OF INTERSECTION
            ZDT=max(zdtmin,(AMIN1(ZT,ZT2)-AMAX1(ZTP,ZT1)))
C  STORE
            ZDTS(JT2)=ZDT
            ZWSUM=ZWSUM+ZDT
C  CHECK FOR END OF INTERSECTION
            IF(ZT.GE.ZT2) GO TO 110
C  NEXT INTERVAL
 90         CONTINUE
            ZTP=ZT
C
 100     CONTINUE
         JT2=INUMT
 110     CONTINUE
         IT2=JT2
         IT0=IT1
C  NOW CALCULATE AVERAGE FOR EACH X AT CURRENT TIME
         DO 150 IX=1,INX
            ZASUM=0.0
            DO 140 IT=IT1,IT2
               ZDAT=zdata_in(ix,it)
C  CHECK FOR DBL INV AVERAGING
               IF(IOPT.EQ.2) THEN
                  if(zdat.eq.0.0) zdat=1.e-35
                  IF(ABS(ZDAT).LT.1.E-35) ZDAT=sign(1.E-35,zdat)
                  ZDAT=1.0/ZDAT
               ENDIF
               ZASUM=ZASUM+ZDAT*ZDTS(IT)
 140        CONTINUE
C
            ZAVG=ZASUM/ZWSUM
C  CHECK FOR DBL INVERT BACK
            IF(IOPT.EQ.2) THEN
               IF(ABS(ZAVG).LT.1.E-35) ZAVG=1.E-35
               ZAVG=1.0/ZAVG
            ENDIF
            zdata_out(ix,itout)=ZAVG
 150     CONTINUE
 200  CONTINUE
      RETURN
      END
C******************** END FILE SMTIMA.FOR ; GROUP SMOOF2 ******************
