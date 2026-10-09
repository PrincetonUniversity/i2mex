      subroutine tconnect(luout, zdisk, zdir, zrunid, idrun, ierr)
C
      use tconnect_mod
      use cplotr_mod
C
C  dmc 18 Nov 1997
C  (trread library routine -- support trscalar & trprofil)
C  connect to a run.  This means:  read the labels, read the
C  timebases, read the scalar functions, open the profile functions.
C
C  Both NetCDF and traditional TRANSP output data formats are supported.
C  MDS+ support added early 2000...
C
      integer luout                     ! lun for messages
C  if(luout.lt.0) check if run is known w/out connecting...
      character*(*) zdisk               ! disk
      character*(*) zdir                ! directory (path)
      character*(*) zrunid              ! TRANSP runid
C
      integer idrun                     ! run id code (returned)
      integer ierr                      ! error code (returned, 0 = no error)
C
C----------------------------
C  local stuff
C
      logical lwait,lcache,ixcache
      integer ilunt
C
      character*40 zrlbl
C
      character*36 zstr
C
      character*128 zdiri
      character*48 zdiski
C
      CHARACTER*64 MDS_SERVER
      CHARACTER*24 MDS_TREE
      INTEGER      MDS_SHOT
      integer socket
      integer stat
      integer ilz1, ilz2
C
      logical :: mds_cache_exclusive
      external :: mds_cache_exclusive
C
      data zstr/'abcdefghijklmnopqrstuvwxyz0123456789'/
C------------------------------------------------------
C  initializations...
C
      ilunt=abs(luout)
C
      call initcf
      call dmgini_chk
C
C  make label string; get non-blank string lengths...
C
      ilz=len(zrlbl)
C
      idrun=0
      ierr=0
      iwarn=0
C
      call tmklbl(zdisk,ildisk,zdir,ildir,zrunid,ilrunid,zrlbl,ilrlbl)
C
      zdiski=zdisk(1:ildisk)
      if(zdiski(1:4).eq.'MDS+') then
         ixcache = mds_cache_exclusive(0)
         if(ixcache) then
            write(ilunt,*) ' ?Tconnect -- secondary run access unsafe;'
            write(ilunt,*) '  MDS+ RPLOT_CACHE_ONLY appears to be set.'
            ierr=1000
            return
         endif
      endif
C
C  standardize punctuation...
C
      zdiri=zdir(1:ildir)
      if(zdir(1:ildir).ne.' ' .and. zdisk(1:4) .ne. 'MDS+') then
         if(zdir(1:1).eq.'$') then
            ilz1 = index(zdir,'/')-1
            if (ilz1 .le. 0) ilz1=ildir
            call get_environment_variable(zdir(2:ilz1),zdiri)
            ilz2 = index(zdiri,' ')
            if(zdir(ildir:ildir).ne.'/') then
               zdiri(ilz2:)=zdir(ilz1+1:ildir)//'/'
            else 
               zdiri(ilz2:)=zdir(ilz1+1:ildir)
            endif
         else if(zdir(ildir:ildir).ne.'/') then
            zdiri=zdir(1:ildir)//'/'
         endif
      endif
C
      icount=0
C
C  check for prior tconnect to this run...
C
      iskip=0
      if (zdiski(1:4) .eq. 'MDS+') then
         if(tc_mds_open.eq.1) then
            iskip=1
            if(save_disk.ne.zdisk) iskip=0
            if(save_dir.ne.zdir) iskip=0
            if(save_runid.ne.zrunid) iskip=0
            if(zdiski(6:20) .eq. 'TRANSPGRID.PPPL'  .or.    
     >         zdiski(6:21) .eq. 'TRANSPGRID1.PPPL' .or.    
     >         zdiski(6:21) .eq. 'TRANSPGRID2.PPPL')
     >         iskip=0
         endif
      endif
C
 5    continue
C
C  look to see if run label info is already stored...
C
      do ir=1,nrun_x
         if(ilrlbl.eq.lrlbl(ir)) then
            if(zrlbl.eq.rlbl(ir)) then
               iok=1                    ! verify the match...
               icount=icount+1
               if(zdiski.ne.fdisk_x(ir)) iok=0
               if(zdiri.ne.fdir_x(ir)) iok=0
               if(iok.eq.1) then
                  idrun=ir
                  go to 100             ! match
               else if(icount.eq.1) then
                  ipc=min(ilz,ilrlbl+2)
                  ilrlbl=ipc
                  zrlbl(ipc-1:ipc)='_'//zstr(icount:icount)
                  go to 5               ! restart search
               else if(icount.gt.1) then
                  zrlbl(ipc:ipc)=zstr(icount:icount)
                  go to 5               ! restart search
              endif
            endif
         endif
      enddo
C
C  no match -- read labels, etc...
C
      if(luout.lt.0) return
C
      inrp=nrun_x                       ! save for error recovery
      ikrp=krun_x
C
      if(nrun_x.lt.naxxtra) then
         nrun_x=nrun_x+1
         krun_x=krun_x+1
      else
         krun_x=krun_x+1
         if(krun_x.gt.naxxtra) krun_x=1
      endif
C
      lrun_x=krun_x
      lrlbl(lrun_x)=ilrlbl
      rlbl(lrun_x)=zrlbl
C
      fdisk_x(lrun_x)=zdiski
C
      fdir_x(lrun_x)=zdiri
C
      runid_x(lrun_x)=zrunid(1:ilrunid)
 
      if (zdiski(1:4) .eq. 'MDS+') then
         nlmds_x(lrun_x) = .true.
      endif
C
      if(iskip.eq.0) then
         call zconnect(ierr,iwarn)
         if(ierr.eq.0) then
            if(zdiski(1:4) .eq. 'MDS+') then
               call mark_mds_open
            endif
         endif
      else
         ierr=0
         iwarn=0
      endif
C
      if(ierr.ne.0) then
         write(ilunt,'('' %tconnect -- error opening run '',a)') zrunid
         nrun_x=inrp
         krun_x=ikrp
         rlbl(lrun_x)='?invalid'
         lrlbl(lrun_x)=ilnurd(rlbl(lrun_x))
         go to 1000
      else if(iwarn.ne.0) then
         write(ilunt,'('' %tconnect -- warning flag, run '',a)') zrunid
      endif
C
      idrun=lrun_x
C
      go to 1000
C
C  match -- but f(t) data may need to be reread.
C
 100  continue
      if(luout.lt.0) return
C
      lrun_x=idrun
C
      if(iskip.eq.0) then
         if(nlmds_x(lrun_x)) then
C
C  reconnect to MDS+ run
C
C  connect to the correct server...
C
            lwait=.false.
            lcache=.false.
            call mds_connect(lwait,lcache,ierr)
            if(ierr.ne.0) go to 1000
            call mark_mds_open
C
C  OK
C
         endif
C
         if(nlcdf_x(lrun_x)) then
            call dmgfot(4,ipt,ierr)     ! NetCDF read (if needed)
         else
            call dmgfot(1,ipt,ierr)     ! TF.PLN read (if needed)
         endif
      endif
C
      go to 1000
C
 1000 continue
      lrun_x=0
      return
 
      contains
        subroutine zconnect(ierr,iwarn)
C
           character*140 ztfile,znfile,zmfile,zcdfile
 
           CALL PLFILN('.CDF',ZCDFILE)
           CALL PLFILN('TF.PLN',ZTFILE)
           CALL PLFILN('MF.PLN',ZMFILE)
           CALL PLFILN('NF.PLN',ZNFILE)

           call pconnect(zcdfile,ztfile,zmfile,znfile,ierr,iwarn)
 
           return
 
        end subroutine zconnect
 
        subroutine mark_mds_open

           write(luout,*) ' %connected to 2nd runid via MDSplus'
           tc_mds_open=1
           tc_idrun=lrun_x
           save_disk=zdisk
           save_dir=zdir
           save_runid=zrunid

        end subroutine mark_mds_open

      end
