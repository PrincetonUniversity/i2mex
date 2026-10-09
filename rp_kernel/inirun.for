C******************** START FILE INIRUN.FOR ; GROUP PLOTR1 *************
C--------------------------------------------------------------
C  INITIALIZE TO LOOK AT ONE RUN
C
C  mod CAL oct 1997 -- NetCDF interface
C  Mod TBT Jan 1994 -- Cleanup IxMax
C  mod dmc jul 1992 -- unix operability
C
      SUBROUTINE INIRUN(IER)
C
C  DMC JUNE 1989 -- NEW PASSED ARGUMENT ISMIN GIVES THE MINIMUM MEMORY
C  WORKSPACE SIZE IN DATMGR, SET IN RPLOT MAIN.
C
C  dmc 1999 -- ISMIN argument is now NWSMIN, in DATMGR COMMON
C   *** ISMIN argument removed, now is NWSMIN in COMMON
C
C  11/10/99 CAL -- check for char = '0'
C
      use datmgr_mod
C
      use cplotr_mod
C
      CHARACTER(len(fdisk)) :: ZDISK
      CHARACTER(len(fdir))  :: ZDIR
      integer      iwait
C
      logical itransp
C
C
C---------------------------------------------------------
C
C 02/25/00 CAL: use ier=-777 as "wait for restore" flag
      if (ier .eq. -777) then
         iwait= ier
      else
         iwait=0
      endif
      IWARN=0
C
      itransp = transp_imbed .and. (lrun_x.eq.0)
C
      IBLNK=INDEX(RUNID,' ')
      IF(IBLNK.EQ.0) IBLNK=LEN(RUNID)+1
      ILNB=IBLNK-1
C
      LRUNID=ILNB
C
C  RUN I.D. LABEL
      if(runlb2.eq.' ' .or. ichar(runlb2(1:1)) .le. 32) then
         IF(LFDIR.EQ.0) THEN
            CALL SHOWDEFL(ZDISK,ILDSK,ZDIR,ILD)
         ELSE
            ILD=LFDIR
            ZDIR=FDIR
         ENDIF
C
C  NEW RUN LABEL
         CALL MRUNLB2(ZDIR,ILD,RUNID,LRUNID,RUNLB2)
      endif
C
      CALL C9DATE(IDATE)
C
      ITITLE=' '
      WRITE(ITITLE,5005) IDATE
 5005 FORMAT('RPLOT GENERATED PLOT ',A9)
C
C
C  CONSTRUCT FILENAMES
      CALL PLFILN('.CDF',CDFILN)
      CALL PLFILN('TF.PLN',TFILN)
      CALL PLFILN('MF.PLN',MFILN)
      CALL PLFILN('NF.PLN',NFILN)
C
C  CLEAR INDEX-GENERATION COUNTER AND PAGE # COUNTER
      NPAGEG=0
      NLSENT=0
C
C  CLEAR DATA AREA
      CALL DMGINI
      if(itransp) no_delete=.TRUE.
C
C  clear time shift tags
      ttagr = 0.0
      ttagt = 0.0
      nlxtrap0=.TRUE.  ! set flag for extrapolation-to-zero in time
C
C  first check access permission
C
      if(.not.itransp) then
         call vprotec(runlb2,ier)
         if(ier.ne.0) then
            ier=1
            go to 1000
         endif
      endif
C
C  READ LABELS, NFT,NFR, ETC.
C
      ier=iwait
      call pconnect(CDFILN,TFILN,MFILN,NFILN,ier,iwarn)
      if(ier.ne.0) go to 990
C
C  DETERMINE SIZE OF WORK SPACE AND ALLOCATE
C
      CALL INIWRK(NWSMIN)
C
      if(itransp) then
C
C  create DATBUF slots, all zero data, for each function
C  TRANSP will provide data later.
C
         do if=1,nfxt
            call dmgfxt(if,ind)
         enddo
      endif
C
C%	CALL DMPRIN
C
C  SET UP GEOMETRY FACTORS FOR INTEGRATIONS, AVERAGES, ETC.
      CALL FIXGEO
C
C  compute initial sort order on names
C  ...names can be added later.
C
      call aordr(iordrt,abt,nft)
      call aordr(iordrr,abr,nfxt)
      call aordr(iordrb,abb,nbal)
C
      GO TO 1000
C------------------------------------------
C  ERRORS...
C
 990  CONTINUE
C
C------------------------------------------
 1000 CONTINUE
C
      RETURN
      END
C******************** END FILE INIRUN.FOR ; GROUP PLOTR1 ***************
