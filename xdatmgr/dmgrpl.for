C-------------------------------------------------------------
C  DMGRPL
C
C  REPLACE EXISTING ENTRY OR ENTRIES WITH NEW ENTRY
 
 
      SUBROUTINE DMGRPL(ISIZE,JLOC,IPRIO)
 
 
C Update:
C  1/20/94 TBT Added Include RPLOTR and check for user function deletion.
C 09/28/93 TBT Switched printout variables Isize & priority.
 
 
      use datmgr_mod
 
C
C  IPRIO -- PRIORITY OF ALLOCATION
C  ISIZE -- SIZE OF NEW ENTRY (MUST BE .GT. 0)
C  JLOC  -- OUTPUT PTR TO DESCRIPTOR OF NEW ENTRY
C
C  A LARGE NUMBER...
      DATA ZLARGE/1.E35/
C-------------------------------------------
 
C  	.UPDATE ACCESS CODE
      call dmg_macc_incr

      do
         JST=1
         JREPL1=0
         JREPL2=0
         ZAGEM=0.
 
C  	.OUTER LOOP-- SCAN POSSIBLE STARTING POINTS FOR NEW SLOT
 10      CONTINUE
         if(no_delete) go to 200  ! no replacements -- expand buffer
         ILOC0=LOCD(JST)+NWDS(JST)
         IF((NDBSIZ-ILOC0+1).LT.ISIZE) GO TO 200
 
C  	.INNER LOOP-- FIND MINIMUM AGE IN GROUP OF ENTRIES TO BE REPLACED
         ZAGMIN=ZLARGE
         J2=LNEXT(JST)
 20      CONTINUE
 
C  	.CHECK FOR SUFFICIENT PRIORITY TO OVERWRITE ENTRY
         IF(IPRIO.LT.MPRIO(J2)) THEN
            ZAGMIN=-1.0
            GO TO 50
         ENDIF
 
C  	.PRIORITY OK.  GET STATS ON THIS SLOT
         CALL DMISTA(J2,ILOC2,ISIZ2,INPREV,INNEXT)
 
C  	.AGE FACTOR -- ENHANCED BY PRIORITY DIFFERENCE
         ZAGE=FLOAT(MACC-LACC(J2))*(1+(IPRIO-MPRIO(J2)))
 
C  	.INCREASE BY 1/(SPACE "WASTE") FACTOR
         ZWASTE=FLOAT(ISIZ2)/FLOAT(ISIZ2-MIN0(INPREV,INNEXT))
         ZAGE=ZAGE*ZWASTE
         ZAGMIN=AMIN1(ZAGMIN,ZAGE)
 
C  	.SEE IF ELIMINATION OF THIS ENTRY WOULD CREATE ENOUGH SPACE FOR
C  	.NEW ENTRY
         ISIZA=ILOC2+ISIZ2-ILOC0
         IF(ISIZA.GE.ISIZE) GO TO 50
 
C  	.CHECK NEXT ENTRY IF NEED MORE SPACE
         J2=LNEXT(J2)
         if(j2.eq.0) go to 200
         GO TO 20
 
 50      CONTINUE
C  	.IF THIS SET IS OLDER THAN PREVIOUSLY ENCOUNTERED SETS, SAVE IT
C  	.AS ELIGIBLE FOR REPLACEMENT
         IF(ZAGMIN.LT.ZAGEM) GO TO 90
C
         JREPL1=JST
         JREPL2=J2
         ZAGEM=ZAGMIN
C  	.END OF OUTER LOOP
 90      CONTINUE
         JST=LNEXT(JST)
         GO TO 10
 
C  	.LOOP EXIT
 
 200     CONTINUE
C       .CHECK FOR FAILURE TO ALLOCATE
         IF(JREPL1.EQ.0) THEN
            write(6,*) ' dmg_datbuf_expand in dmgrpl.for, isize=',isize
            call dmg_datbuf_expand(0)
            isiza=abs(isize)
            call dmgbsf(isiza,jloc,iprio)
            if(jloc.gt.0) return

            cycle  ! still no slot found

         ELSE
            exit
         ENDIF
      enddo
C     
C  	.DELETE ENTRIES TO BE SUPERSEDED
C
      J=JREPL1
 210  CONTINUE
      J=LNEXT(J)
 
      CALL DMIDEL(J)
      IF(J.EQ.JREPL2) GO TO 250
      GO TO 210
C
 250  CONTINUE
C  	.INSERT NEW ENTRY
      CALL DMINEW(JREPL1,JLOC,IPRIO)
 
      RETURN
      END
