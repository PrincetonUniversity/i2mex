C******************** START FILE GETTM3.FOR ; GROUP PLOTR3 ******************
C===========================================================
C  GETTM3  GET TIME VECTOR FOR 3D PROFILE FCNS OF TIME AND
C          RADIUS
 
 
 
      SUBROUTINE GETTM3(IER,ITEXP)
 
c
c  mod dmc May 1994 -- added argument ITEXP = no. of time pts expected.
c   if ITEXP=0, the expected number is not known
c
c   I use this to avoid some error checking code that may not be reliable
c   on all machines.  this mechanism improves the reliability of poplt2.
c
 
C  Mod TbT Jan 1994 -- Deleted Zdummy - Put in 2 Ztime's.
C  MOD DMC NOV 1988 -- AS POPLOT SUBROUTINE -- ABILITY TO IGNORE
C  SOME OUT OF SEQUENCE TIME POINTS
 
C  COMMON BLOCKS---

      use datmgr_mod
      use cplotr_mod
      use mfblok_mod
C
      logical itestmf
C
C----------------------
C
      if(itexp.gt.0) call dmg_texpand(itexp)
C
      NT0=NTIME
C
      IER=0
C
      if(lrun_x.eq.0) then
         itestmf=MFBLKI
         imflun=MFLUNI
      else
         itestmf=MFBLKI_X(lrun_x)
         imflun=MFLUNI_X(lrun_x)
      endif
C
      IF(itestmf) THEN
        write(lunzer(0),*)
     >      ' ?GETTM3, NOT READY FOR MFBLKI=.TRUE. FORMAT'
        call bad_exit
      ENDIF
C
C  READ SEQUENCE OF TIMES ON FIRST PASS FOR GIVEN RUN #
      ICT=0
      ICTOT=0
      ICUSE=0
      IREC=-NFR
C
      ZTIMEP=-1.0e10
C
 10   CONTINUE
      ICT=ICT+1
      ICTOT=ICTOT+1
      IREC=IREC+NFR+1
      if(itexp.gt.0) then
         if(ictot.gt.itexp) go to 90
      endif

      if(ict.gt.ntime) call dmg_texpand(0)

      READ(IMFLUN,REC=IREC,IOSTAT=IERR,ERR=90) ZTIME
      IF(IERR.NE.0) GO TO 90
      IF(LTWRIT(ICT).or.(lrun_x.ne.0)) THEN
        ICUSE=ICUSE+1
        if(lrun_x.eq.0) then
           TIME3(ICUSE)=ZTIME
           ZTIMEP=ZTIME
        else
           TIME3_X(icuse,lrun_x)=ZTIME
           ZTIMEP=ZTIME
        endif
      ELSE
         ICTOT=ICTOT-1
      ENDIF
      IF(ICUSE.GE.NT0) GO TO 50
      GO TO 10
C
 50   CONTINUE
      ICUSE=NT0
      WRITE(lunzer(0),9001) NT0,ZTIMEP
 9001 FORMAT(
     >  ' %GETTM3:  PROGRAM ARRAY CAPACITY OF ',I4,' TIME POINTS'/
     >  '   EXCEEDED, TIMES AFTER T=',1PE10.3,' SECONDS IGNORED')
      IER=1
      GO TO 95
C
 90   CONTINUE
C  CHECK FOR ERROR
      IF((ITEXP.EQ.0).OR.(ICUSE.NE.ITEXP)) THEN
        READ(IMFLUN,REC=25000000,IOSTAT=IERR2,ERR=91) Ztime ! = dummy
 91     CONTINUE
        IF(IERR.EQ.IERR2) GO TO 95
        WRITE(lunzer(0),8001) IERR
 8001   FORMAT(
     >' %GETTM3:  "MF" FILE ERROR, IOSTAT CODE = ',I5/
     >  '   EXECUTION CONTINUING, BUT DATA MAY BE LOST OR DAMAGED')
      ENDIF
C
 95   CONTINUE
C  CHECK FOR FILE ERROR-- INCOMPLETE TRANSFER OF LAST TIME RECORD
      IREC=IREC-1
      READ(IMFLUN,REC=IREC,ERR=99) Ztime   ! = dummy
      GO TO 100
C  ERROR IN LAST RECORD
 99   CONTINUE
      ICUSE=ICUSE-1
      WRITE(lunzer(0),8003) ICUSE
 8003 FORMAT(
     >' ?GETTM3:  IMCOMPLETE TIME RECORD IN "MF" FILE IGNORED--'/
     >  '   DATA LOST.  ',I5,' TIME RECORDS RETAINED')
      IER=3
C
C  TIME VECTOR IN HAND; PROCEED TO FETCH DESIRED FUNCTIONS
C
 100  CONTINUE
      if(lrun_x.eq.0) then
         NTR=ICUSE
      else
         NTR_X(lrun_x)=ICUSE
      endif
      RETURN
      	END
C******************** END FILE GETTM3.FOR ; GROUP PLOTR3 ******************
