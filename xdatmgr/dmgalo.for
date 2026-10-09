C******************** START FILE DMGALO.FOR ; GROUP DATMGR ******************
C-----------------------------------------------------------------
C  DMGALO
C
C  ALLOCATE SPACE FOR NEW DATA ENTRY
C
      SUBROUTINE DMGALO(ISIZE,JLOC,IPRIO)
C
C  ISIZE-- INPUT.  IF ISIZE.GT.0 A SLOT IS ALLOCATED IN BUFFER ON
C    BEST AVAILABLE FIT BASIS (IF AN UNUSED SLOT IS AVAILABLE), OR,
C    IF AN OLD ENTRY MUST BE OVERWRITTEN A SUITABLE ENTRY IS FOUND
C    (SEE SUBROUTINE DMGRPL)
C         IF ISIZE.EQ.0 THE LARGEST AVAILABLE FREE SLOT IS ALLOCATED
C  JLOC-- OUTPUT.  POINTER TO DESCRIPTORS OF ALLOCATED SPACE
C  IPRIO-- INPUT.  PRIORITY CODE FOR THIS ALLOCATION
C    CANNOT DISPLACE ITEMS OF HIGHER PRIORITY
C
      use datmgr_mod
C
      itest = 10*isize
      itest2 = itest
      if(itest.gt.(4*ndbsiz_min)) then
         itest = 3*isize + ndbsiz_min/2
         itest2 = 3*isize + ndbsiz_min
      endif

      if(itest.gt.NDBSIZ) then
         write(6,*) ' dmgalo dmg_datbuf_expand isize,itest=',isize,itest
         call dmg_datbuf_expand(itest2)
      endif

      do
         JLOC=0
C  CHECK FOR BEST FIT-- FREE SLOT
         ISIZA=IABS(ISIZE)
         CALL DMGBSF(ISIZA,JLOC,IPRIO)
C  FOUND A FREE SLOT?
         IF(JLOC.GT.0) RETURN
C  NO FREE SLOT-- REPLACE EXISTING ENTRIES IN BUFFER
         CALL DMGRPL(ISIZE,JLOC,IPRIO)
         IF(JLOC.EQ.0) THEN
            write(6,*) ' dmg_datbuf_expand (jloc=0 after dmgrpl).'
            call dmg_datbuf_expand(0)
            cycle
         ELSE
            exit
         ENDIF
      enddo

      RETURN
      END
C******************** END FILE DMGALO.FOR ; GROUP DATMGR ******************
