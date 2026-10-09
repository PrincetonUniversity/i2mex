      subroutine pconnect(zcdfile,ztfile,zmfile,znfile,
     1                   ier,iwarn)
C
C  dmc 17 Nov 1997 -- connect to a run
C
C  this can be used by rplot_sub/inirun.for to connect to a "primary"
C  run for RPLOT use, or, by trprofil/trscalar, to connect to a "secondary"
C  run via the callable data interface
C
C  if a "primary" run then CPLOTR COMMON variable LRUN_X = 0 must be set
C  if a "secondary" run the LRUN_X .gt. 0 must be set
C
C  the code was extracted from inirun.for 17 Nov 1997.
C
C  History:
C  01/06/00 CAL:  support MDSplus
C  02/25/00 CAL:  kludge for MDSplus wait-restore
C
C  ier flag:  on entry
C        ier=-777 -- wait for restore of offline data
C          otherwise set ier=2 if a (slow) restore request is issued
C  ier flag:  on exit,
C        ier=0 -- normal completion
C        ier=1 -- error opening or reading data
C        ier=2 -- data is offline, being restored
C-----------------------------------------------------------------------
C
      use datmgr_mod
      use cplotr_mod
      use mfblok_mod
C
      CHARACTER*24 ZRUNID  ! lengthened dmc 13 Aug 1996
C
C  arguments
C
C input:
C
      character*(*) zcdfile,ztfile,zmfile,znfile ! filenames
C                       .CDF  TF.PLN MF.PLN NF.PLN
C output:
C
      integer ier               ! completion code, 0= normal
      integer iwarn             ! warning flag, 0= no warning
C
C------------------------------
C  local stuff
C
      integer iidcdf
      logical inlcdf, inlmds, itransp
C
      logical imfblk, lwait, lcache
C
C--------------------------------------------------------------------
C
C     call cplset_exec          ! link the block data...
C     call rpcaldat_exec        ! link the block data...
C
C 02/25/00 CAL: use ier= -777 as "wait" flag
      if (ier .eq. -777) then
         lwait = .true.
      else
         lwait = .false.
      endif
      luntrm=lunzer(0)
      ier=0
C
      itransp = transp_imbed .and. (lrun_x.eq.0)
C
C  iwarn is not cleared -- value may have been set prior to pconnect call
C
      ZRUNID=RUNID
C
      if(lrun_x.eq.0) then
         iidcdf=idcdf
         inlcdf=nlcdf
         inlmds=nlmds
         ilf=len_trim(zcdfile)-4
         curstring=zcdfile(1:max(1,ilf)) ! stash ID, current run
c        ! curstring updated in mds_connect, if this is an MDS+ run
      else
         iidcdf=idcdf_x(lrun_x)
         inlcdf=nlcdf_x(lrun_x)
         inlmds=nlmds_x(lrun_x)
      endif
C
C Check if open netCDF file: close
C
      if (iidcdf .gt. 0) then
         call cdfilcl(iidcdf, ier)
         ier = 0
         inlcdf = .false.
         iidcdf = 0
      end if
C
C if inside TRANSP go straight to TF.PLN file ***
C
      if(itransp) then
         ier=1                          ! force TF read
         go to 10
      endif
C
C Check if MDSplus
C=================
      if (inlmds) then
         write(luntrm,*) ' ...connecting to MDSplus'
         lcache=.TRUE.
         call mds_connect(lwait,lcache,ier)
C Tree successfully opened
         if (ier .eq. 0) then
            write(luntrm,*) ' ...reading MDSplus header data...'
            call mdshrd(luntrm,ier,lrun_x)
            if (ier .eq. 0) then
               go to 17
            else
               if (lrun_x.eq.0) then
                  nlmds = .false.
               else
                  nlmds_x(lrun_x) = .false.
               endif
               ier=1
               return
            endif
C  Run is being restored
         else if (ier .eq. 2) then
            write(luntrm,*) ' %PCONNECT: '//zrunid//
     >         ' is beeing restored'
         else if (ier .eq. -1) then
            ier=1
            write(luntrm,*) ' ?PCONNECT: '//zrunid//' does not exist'
         else if (ier .eq. 4) then
            ier=1
            write(luntrm,*) ' ?PCONNECT: no privilege accessing '//
     >            zrunid
         else
            ier=1
            write(luntrm,*) ' ?PCONNECT: failed reading MDSplus'
         endif
         if (lrun_x.eq.0) then
            nlmds = .false.
         else
            nlmds_x(lrun_x) = .false.
         endif
         return
      endif
C
C Try to open netCDF file first
C
      call cdfilop(ZCDFILE, iidcdf, ier)
      if(ier.ne.0) then
C (unix only):  try lowercase filename
         write(luntrm,*) '(retry folding filename to lowercase)'
         ilcdf=ilnurd(zcdfile)
         do ic=ilcdf,1,-1
            if(zcdfile(ic:ic).eq.'/') then
               ick=ic+1
               go to 11
            endif
         enddo
         ick=1
 11      continue
         call ulower(zcdfile(ick:ilcdf))
         call cdfilop(zcdfile,iidcdf,ier)
      endif
      if (ier .eq. 0) then
         inlcdf = .true.
         call cdfhrd(iidcdf,ier,LRUN_X)
         if (ier .eq. 0) then
            go to 17
         else
            ier=1
            call cdfilcl(iidcdf, i)
            inlcdf = .false.
         end if
      else
         ier=1
      end if
C
C If error or skipcdf, or private TF.PLN, try .PLN
C
 10   continue
C
      if (ier .ne. 0) then
         iidcdf = 0
         write(luntrm,*) ' ...reading TF.PLN header data...'
         tfiln=ztfile
         ier=0
         CALL TFILRD(lun_tf,IER,LRUN_X)
 
         IF(IER.EQ.0) GO TO 17
         ier=1
 
         IBL=INDEX(ZTFILE,' ')-1
         IF(IBL.LE.0) IBL=LEN(ZTFILE)
         WRITE(LUNTRM,2000) ZTFILE(1:IBL)
 2000    FORMAT(/' ? PLOT LABEL FILE "',A,'" NOT FOUND -- '/
     >      '   MDSplus/NetCDF access also unsuccessful.')
C
         GO TO 990
      end if
C
C  labels have been read successfully *****
C
 17   CONTINUE
C
      if(lrun_x.eq.0) then
         idcdf=iidcdf
         nlcdf=inlcdf
      else
         idcdf_x(lrun_x)=iidcdf
         nlcdf_x(lrun_x)=inlcdf
      endif
C
C  RUNID warning (primary runs only)
C
      IF(ZRUNID.EQ.RUNID) GO TO 16
      if(LRUN_X.eq.0) then
         if (inlcdf) then
            ilnb=index(zcdfile,' ')-1
            WRITE(LUNTRM,9916) ZCDFILE(1:ilnb),RUNID
         else
            ilnb=index(znfile,' ')-1
            WRITE(LUNTRM,9916) ZNFILE(1:ilnb),RUNID
         end if
 9916    FORMAT(' % RUN ID IN "',A,'" IS ',A)
         IWARN=1
      endif
      RUNID=ZRUNID              ! restore RUNID
 
 16   CONTINUE
C
C  DMC AUG 1985 - MAKESURE ALL ABBREVIATIONS ARE CAPITALIZED
C
      if(LRUN_X.eq.0) then
         CALL PLCAP(XNDABB,NXR)
         CALL PLCAP(ABR,NFXT)
         CALL PLCAP(ABT,NFT)
         CALL PLCAP(ABB,NBAL)
 
         ! decide whether this is a transp data set
         if (nxr<2) then
            nltransp = .false.
         else
            nltransp = ((xndabb(1)=='X' .and. xndabb(2)=='XB')
     >           .or. (xndabb(1)=='RZON' .and. xndabb(2)=='RBOUN'))
     >           .and. (nzonex(1)==nzonex(2))
         end if
      else
         CALL PLCAP(abr_x(1,lrun_x),nfxt_x(lrun_x))
         CALL PLCAP(abt_x(1,lrun_x),nft_x(lrun_x))
      endif
C
C  OCT 1987 -- ADDING USER DEFINED FCNS; BE SURE TO SAVE NO. OF
C  FCNS IN ORIGINAL DATA FILE F(X,T) (primary run only)
C
      if(LRUN_X.eq.0) then
         NFXT0=NFXT
      endif
C
      if (.not.inlcdf .and. .not. inlmds .and. .not. itransp) then
C  OPEN PLOTTING DATA FILES
C
C  GET RECORD SIZE IN BYTES; IN OPEN STATEMENT USE RECSIZE IN LONGWORDS
C
         if(lrun_x.eq.0) then
            call mfinq(ZMFILE,NZONES,IBLKSZ,IRECSZ,MFBLKI)
            imfblk=MFBLKI
         else
            call mfinq(ZMFILE,NZONES_X(LRUN_X),IBLKSZ,
     >         IRECSZ,MFBLKI_X(LRUN_X))
            imfblk=MFBLKI_X(LRUN_X)
         endif
C
C  OPEN NF FILE
C
         close(unit=lun_nf,err=911)
 911     continue
         CALL GENOPEN(lun_nf,ZNFILE,'OLD','BINARY',0,IOS)
C
         IF(IOS.NE.0) GO TO 91
C
         if(LRUN_X.eq.0) then
            MFLUNI=lun_mf
            imflun=MFLUNI
         else
            imflun=MFLUNI_X(LRUN_X)
         endif
C
         close(unit=imflun,err=921)
 921     continue
         CALL GENOPEN(imflun,ZMFILE,'OLD','DIRECT',IRECSZ,IOS)
C
         IF(IOS.NE.0) GO TO 92
C
C  OK
         GO TO 95
C     ERRORS
 91      CONTINUE
         WRITE(LUNTRM,9001) ZNFILE
 9001    FORMAT(/' ? FAILED TO OPEN F(T) DATA FILE - ',A)
         IER=1
         GO TO 990
C
 92      CONTINUE
         CLOSE(UNIT=lun_nf)
         WRITE(LUNTRM,9002) ZMFILE
 9002    FORMAT(/' ? FAILED TO OPEN F(X,T) DATA FILE - ',A)
         IER=1
         GO TO 990
C
 95      CONTINUE
C
C  IF USING NEW BLOCKED MF FILE STRUCTURE READ AND CHECK HEADER
C
         IF(IMFBLK) THEN
            if(lrun_x.eq.0) then
               call mfhdrd(imflun,NFXT,IER)
               if(ier.ne.0) ier=1
            else
               call mfhdrd(imflun,NFXT_X(LRUN_X),IER)
               if(ier.ne.0) ier=1
            endif
            if(ier.ne.0) then
               CLOSE(UNIT=lun_nf)
               CLOSE(UNIT=imflun)
               go to 990
            endif
         endif
C
C  READ MF TIME VECTOR AND SET UP WORKSPACE
         IER=0
         IF(IMFBLK) THEN
            if(lrun_x.eq.0) then
               isizb=1+(mfhdr(2)-1)/iblksz
               call dmg_texpand(isizb*iblksz)
               call mftimrd(imflun,mfhdr,iblksz,time3,NTIME,ntr)
            else
               isizb=1+(mfhdr_x(2,lrun_x)-1)/iblksz
               call dmg_texpand(isizb*iblksz)
               call mftimrd(imflun,mfhdr_x(1,lrun_x),iblksz,
     >            time3_x(1,lrun_x),NTIME,ntr_x(lrun_x))
            endif
         ELSE
            CALL GETTM3(IWARN,0)
         ENDIF
         if(lrun_x.eq.0) then
            IF(NTR.EQ.0) IER=1
         else
            if(ntr_x(lrun_x).eq.0) IER=1
         endif
      end if
C  End .PLN files
C-----------------
C  default profile timebase -- when imbedded in TRANSP
C
      if(itransp) then
         ntr=2
         time3(1)=transp_tinit-1
         time3(2)=transp_tinit
      endif
C
C  READ SCALAR DATA
      if (inlcdf) then
         call dmgfot(4,ipt,ier)
         if(ier.ne.0) ier=1
      else
         CALL DMGFOT(1,IPT,IER)
         if(ier.ne.0) ier=1
      end if
C
      IF(IER.NE.0) THEN
         WRITE(LUNTRM,9014)
 9014    FORMAT(' ?PCONNECT -- I/O ERROR ON scalar f(t) data')
         IER=1
         CLOSE(UNIT=lun_nf)
         CLOSE(UNIT=imflun)
         if (inlcdf) then
            call cdfilcl(iidcdf,ier2)
            inlcdf = .false.
         end if
         GO TO 990
      ENDIF
C
C  CHECK TIME AXES
C
      if(lrun_x.eq.0) then
         CALL TIMCK1(TIME3,NTR)
      else
         CALL TIMCK1(TIME3_X(1,lrun_x),ntr_x(lrun_x))
      endif
      go to 1000
C
C  an error occurred...
C
 990  continue
      ier=max(1,ier)
      return
C
 1000 continue
C
      if(lrun_x.eq.0) then
         idcdf=iidcdf
         nlcdf=inlcdf
      else
         idcdf_x(lrun_x)=iidcdf
         nlcdf_x(lrun_x)=inlcdf
      endif
C
      return
      end
 
 
 
 
 
 
 
 
 
