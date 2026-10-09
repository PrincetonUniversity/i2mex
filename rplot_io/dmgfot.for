C******************** START FILE DMGFOT.FOR ; GROUP PLDMGR ***********
C--------------------------------------------------------------
C  PLOTR-- DATA MGR INTERFACE
C
C  DMGFOT-- READ SCALAR FCNS OF TIME INTO DATA AREA
C  DMGXOT-- READ TIME-VARYING X AXIS INTO DATA AREA
C  DMGGEO-- READ TIME-VARYING GEOMETRY INFO INTO DATA AREA
C  DMGFXT-- READ FCN OF TIME + ADDL COORDINATE INTO DATA AREA
C  DMDLOC-- LOCATE DATA IN DATA AREA
C  DMPRIN-- PRINT OUT CONTENTS OF DATA AREA
C
C  dmc 18 Nov 1997 -- support use in context of trprofil/trscalar as
C  well as rplot.
C
      SUBROUTINE DMGFOT(ICALL,IPT,IER)
C
C  READ SCALAR FCNS INTO DATA AREA AND CREATE EXTENDED SPACE FOR
C  USER DEFINED FUNCTIONS
C
C  ICALL INPUT =1 CALL FROM INIRUN -- INITIAL ALLOCATION AND READ
C    FILE
C              =4 Same as 1 for netCDF             cal 09/24/97
C
C  ICALL=2 and ICALL=3 are illegal -- use subroutine dmgfotx instead!
C
C  IPT OUTPUT  = PTR TO SCALAR DATA BLOCK
C
C  IER OUTPUT  = COMPLETION CODE, 0 DENOTES SUCCESS
C
      use datmgr_mod
      use cplotr_mod
C
      character*50 zdlbl
      logical      lmds, itransp
C----------------------------------------
C  dmc 17 Nov 1997 -- support concept of primary or secondary run
C    primary run (lrun_x = 0 in COMMON) -- RPLOT's current runid
C    secondary run (lrun_x .gt. 0 in COMMON) -- for TRPROFIL or TRSCALAR
C
C  REDONE DMC SEPT 1987
C   SCALAR F(T) DATA IS MADE PERMANENTLY RESIDENT IN DATMGR
C
      IER=0
C
      itransp = transp_imbed .and. (lrun_x.eq.0)
C
      if(lrun_x.eq.0) then
         IFXTND=16			! ALLOCATION EXTENSION, NO. USER FCNS
         zdlbl='%F(T)'
         lmds=nlmds
      else
         lmds=nlmds_x(lrun_x)
         IFXTND=0
         zdlbl=RLBL(LRUN_X)(1:LRLBL(LRUN_X))//'!%F(T)'
         call dmdloc(zdlbl,indo,isizo,ipto)
         if(indo.gt.0) then
            ipt=ipto			! already in the buffer
            go to 900
         endif
      endif
C
      IF(ICALL.EQ.1 .or. icall .eq. 4) THEN
         isize=0
C NETcdf
         if(icall.eq.4) then
            if(lrun_x.eq.0) then
               iidcdf=idcdf
               ineed=ntt*nft
            else
               iidcdf=idcdf_x(lrun_x)
               ineed=ntt_x(lrun_x)*nft_x(lrun_x)
            endif
            isize=ineed
C MDSplus
         else if (lmds) then
            if(lrun_x.eq.0) then
               isize=ntt*nft
            else
               isize=ntt_x(lrun_x)*nft_x(lrun_x)
            endif
            ineed=isize
C TRANSP
         else if (itransp) then
            ntt=2
            isize=ntt*nft
            time(1)=transp_tinit-1
            time(2)=transp_tinit
            ineed=isize
         else
C  have to read NF.PLN file to determine space needed
            ineed=0             ! determine by reading file
            if(lrun_x.ne.0) then
               call nfread(lun_nf,0,ineed,ier)
               if(ier.ne.0) go to 900
               isize=ineed
            endif
         endif
         if(lrun_x.eq.0) then
            iarg=-1
            iprio=10
            call dmgbsf(iarg,JTIM,iprio) ! f(t) data for run
         else
            iprio=5                     ! normal priority
            call rp_dmgalo(ineed,JTIM,iprio)
            if(ier.ne.0) go to 900
         endif
C  ERROR CHECKING
         JNEX=LNEXT(JTIM)
         if(lrun_x.eq.0) then
            IF((LOCD(JTIM).NE.1).OR.(LOCD(JNEX).LE.NDBSIZ)) THEN
               write(lunzer(0),*)
     1	    '?RPLOT/DMGFOT - DATA BUFFER INITIALIZATION ERROR!'
               IER=1
               GO TO 900
            ENDIF
         endif
C  READ FCNS INTO DATA AREA
C
         IPT1=LOCD(JTIM)
C
C  If netCDF :
         if (icall .eq. 4) then
            call cdfnfrd(iidcdf,ipt1,isize,ier)
            IF(IER.NE.0) GO TO 900
C  If MDSplus:
         else if (lmds) then
            call mdsnfrd(ipt1,isize,ier)
            IF(IER.NE.0) GO TO 900
C  If .PLN:
         else if (.not.itransp) then
            if(lrun_x.eq.0) then
               isize=0		! will be determined on read
C                       ...and data needs to be inverted after nfread
            else
               isize=ineed      ! nfread will order data
            endif

            if(2*isize.gt.NDBSIZ) then
               write(6,*) ' expand before NFREAD, isize=',isize
               call dmg_datbuf_expand(3*isize)
            endif
            CALL NFREAD(lun_nf,IPT1,ISIZE,IER)
            IF(IER.NE.0) GO TO 900

            if(2*isize.gt.NDBSIZ) then
               write(6,*) ' expand after NFREAD, isize=',isize
               call dmg_datbuf_expand(3*isize)
            endif

C     ERROR CHECK-- initializing space for new primary run
            if(lrun_x.eq.0) then
               IPT2=ISIZE+1
C       COPY DATA INTO SECOND AREA, BUT INVERT S.T. DATA IS MADE
C       TIME-CONTIGUOUS
               DO 100 I=1,NFT
                  DO J=1,NTT
                     IP1=IPT1+(J-1)*NFT+I-1
                     IP2=IPT2+(I-1)*NTT+J-1
                     DATBUF(IP2)=DATBUF(IP1)
                  enddo
 100           CONTINUE
C     COPY DATA BACK INTO FIRST AREA, IN CORRECT ORDER
               CALL copyr4(DATBUF(IPT2),DATBUF(IPT1),ISIZE)
            endif
C
         endif
C       SET SIZE RESERVING ROOM FOR SOME USER FCNS
         NWDS(JTIM)=ISIZE + NTT*IFXTND
         DMGLBL(JTIM)=zdlbl

         if((locd(JTIM)+NWDS(JTIM)).gt.NDBSIZ) then
            write(6,*) ' dmg_datbuf_expand: ',locd(jtim),nwds(jtim)
            call dmg_datbuf_expand(0)
         endif
         IPT=IPT1
C
C       STORE IN NFT0 THE ORIG. NO. OF FCNS IN THE FILE
         if(lrun_x.eq.0) then
            NFT0=NFT
         endif
C
C  icall for packing user fcns, not data read:
C
      ELSE
         write(lunzer(0),*)
     >        ' ?dmgfot -- ICALL error -- call dmgfotx instead!'
         call bad_exit
      ENDIF
C
 900  CONTINUE
C
      if(lrun_x.eq.0) then
         NFTX=0
      endif
C
      RETURN
C
      END
C******************** END FILE DMGFOT.FOR ; GROUP PLDMGR *************
