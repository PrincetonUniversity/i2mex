      subroutine plc_timi(istat,zt0,ier)

      use datmgr_mod
      use cplotr_mod

C  time integrate quantity in RPLOT calculator accumulator
C
C  input
C     istat = type of quantity; -1: scalar, >0: profile
c     zt0   = lower time limit of integration, in character form
C
C  output
C     ier   = completion code, 0 = OK
C
      integer istat
      character*(*) zt0
      integer ier
C
      real zbuf0(nr0)
C-----------------------------------------
C
      ier=0
      if(istat.eq.0) then
         call zermsg(' %TIME_INT error:  calculator data not ready.')
         ier=ier+1
      endif
C
      read(zt0,'(g20.0)',err=9) time0
      go to 10
C
 9    continue
      call zermsg(' %TIME_INT error on decode of TIME0:  '//zt0)
      ier=ier+1
C
 10   continue
      if(ier.gt.0) return
C
      if(istat.lt.0) then
         ztim0=time(1)
         ztimn=time(ntt)
         intl=ntt
         inx=1
      else if(istat.gt.0) then
         ztim0=time3(1)
         ztimn=time3(ntr)
         intl=ntr
         inx=nzonex(istat)
      endif
C
      IF (TIME0 .LE. ztim0)  THEN
         TIME0 = ztim0
         INDEXT= 1
         TRATIO= 0.0
      ELSE IF (TIME0 .GE. ztimn)  THEN
         TIME0 = ztimn
         INDEXT= intl-1
         TRATIO= 1.0
      ELSE                              !  t(1) <= t0 <= t(NTR)
         DO 15 IT=2,intl
            if(istat.lt.0) then
               ztimel=time(it)
               ztimelp=time(it-1)
            else
               ztimel=time3(it)
               ztimelp=time3(it-1)
            endif
            IF (TIME0 .LT. ztimel)  THEN
               INDEXT = IT-1
               GO TO 16
            END IF                      ! time0
 15      CONTINUE
         CALL ZERMSG(' ?PLCFXT: Loop 15 completed - error')
         go to 10
 16      CONTINUE
	
         TRATIO = (TIME0-ztimelp) / (ztimel-ztimelp)
 
      END IF                            ! TIME0
C
      write(lunzer(0),1001) time0
 1001 format(' %plc_timi:  time integrating from TIME0 = ',1pe15.8)
C
C  find accumulator workspace
C
      CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)
C
C  no. of x pts
C
      if(istat.lt.0) then
         inx=1
      else
         inx=nzonex(istat)
      endif
C
C	.Find value of $ at T0. ($ IN WRK2)
      DO 400 IX=1,INX
         ID1 = (IWRK2-1) + (INDEXT-1)*INX + IX
         ID2 = (IWRK2-1) +  INDEXT   *INX + IX
 
         zbuf0(IX) = (1-TRATIO)* DATBUF(ID1) + TRATIO * DATBUF(ID2)
 400  CONTINUE
 
C	.Set integral = 0 at T0.
      DO 430 IT=1,intl
         DO 420 IX=1,INX
            IADL = (IWRK2-1) + (IT-1)*INX + IX
 
            DATBUF(IADL) = DATBUF(IADL) - zbuf0(IX)
 420     CONTINUE
 430  CONTINUE
 
C
      return
      end
