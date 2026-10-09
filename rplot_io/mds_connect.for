      subroutine mds_connect(lwait,lcache,ier)
C
C Connect to MDSplus Server and Open Tree
C
C 02/25/00 CAL
C
      use cplotr_mod
      implicit none
C
      logical       lwait   ! if true:  wait for restore
C                           ! if false: restore in no-wait
      logical       lcache  ! if true:  re-init cache variables
                            ! if false: do not re-init cache variables
 
      integer       ier
 
 
      integer        mds_value, mds_open, mds_close
      integer        mdsfconnect, MdsSetSocket
C Local
      CHARACTER*7  TOKYR
      CHARACTER*64 MDS_SERVER
      CHARACTER*24 MDS_TREE
      INTEGER      MDS_SHOT
      character*10 rid
      character(len=200) :: zcmd
 
      integer      stat, il, il0, mds_stat, retlen, istat
      integer      ilr, str_length
      integer      luntrm, dscr, socket, save_socket
      logical      connect, disconnect, close_all
      integer      lunzer, ilc, i, ilt, iermds, ils
C----------------------------------------------------------------------
      luntrm=lunzer(0)
      ier=0
      connect= .false.
      disconnect= .false.
      close_all = .false.
C
C
C Parse server / tree
C=====================
C
      il0=index(fdisk,'@')-1
C Primary run
      if(lrun_x.eq.0) then
         il=il0
         mds_server=fdisk(6:il)
         mds_tree=fdisk(il+2:)
         if (lfdir .gt. 2) then  ! if fdir = TOK.YY
           if (fdir(lfdir-2:lfdir-2) .eq. "." .and.
     &        fdisk(il-8:il) .eq. ".PPPL.GOV" ) then
             if ((lfdisk .lt. il+8) .and. 
     &          (fdisk(il+2:il+7) .eq. "TRANSP")) then
               mds_tree=fdisk(il+2:il+7)//'_'// fdir(1:lfdir-3)
             endif  
           endif
         endif
         curstring='MDS+:'//fdisk(6:il)//';'//fdisk(il+2:lfdisk)//'('
         ilc=str_length(curstring)
         if(fdisk(5:5).eq.'+') then
            curstring(ilc+1:)=runid(1:lrunid)//')'
         else
            curstring(ilc+1:)=fdir(1:lfdir)//','//runid(1:lrunid)//')'
         endif
         write(luntrm,*) 'mds_connect: ',trim(curstring)
C Scondary runs
      else
         il=index(fdisk_x(lrun_x),'@')-1
         mds_server=fdisk_x(lrun_x)(6:il)
         mds_tree=fdisk_x(lrun_x)(il+2:)
      endif
      write(luntrm,*) ' MDS+ tree & server: ',trim(mds_tree)//' '//
     >     trim(mds_server)

C============================================
C -- Primary Run
      if(lrun_x.eq.0) then
         mdspulse=0
C -- is mds_shot number supplied (e.g. CMOD) ?
         if (fdisk(5:5) .eq. '+') then
            read(runid,'(i10)',iostat=istat) mds_shot
            if(istat.ne.0) then
               write(luntrm,*) '?Mds_connect: not an integer: ',runid
               write(luntrm,*) ' Need '//
     >              'syntax (<tok>,<runid>) or (<tok.yy>,<runid>)'//
     >              ' inside the parentheses.'
               ier=1
               return
            else
               write(luntrm,*) '%Mds_connect: mds_shot =',mds_shot
            endif
         else
            mds_shot=0
            rid = runid
            ilr=lrunid
            tokyr = fdir(1:lfdir)
            call mkbeastid(runid,mds_shot,rid)
            if (mds_shot .eq. 0) then
               write(luntrm,*) 
     >              '?mds_connect: -E-  error translating'//runid
               ier=1
               return
            endif
         endif
      else
C -- Secondary Runs
         mdspulse_x(lrun_x)=0
C -- is mds_shot number supplied (e.g. CMOD) ?
         if (fdisk_x(lrun_x)(5:5) .eq. '+') then
            read(runid_x(lrun_x),'(i10)',iostat=istat) mds_shot
            if(istat.ne.0) then
               write(luntrm,*) '?Mds_connect: not an integer: ',
     >              runid_x(lrun_x)
               write(luntrm,*) ' Need '//
     >              'syntax (<tok>,<runid>) or (<tok.yy>,<runid>)'//
     >              ' inside the parentheses.'
               ier=1
               return
            else
               write(luntrm,*) '%Mds_connect: mds_shot =',mds_shot
            endif
         else
            mds_shot=0
            rid = runid_x(lrun_x)
            il=index(fdir_x(lrun_x),' ')-1
            if (il.le.0) il=len(fdir_x(lrun_x))
            ilr=index(runid_x(lrun_x),' ')-1
            if (ilr.le.0) il=len(runid_x(lrun_x))
            tokyr = fdir_x(lrun_x)(1:il)
            call mkbeastid(runid_x(lrun_x),mds_shot,rid)
            if (mds_shot .eq. 0) then
               write(luntrm,*) 
     >              '?mds_connect: -E-  error translating'//runid
               ier=1
               return
            endif
         endif
      endif
C

      if (lrun_x.eq.0) then

         if(.NOT.mds_cache_only) then
C have previous connection?
            if(mds_save_server.ne.' ') then
C different server?
               if(mds_server.ne.mds_save_server) then
C Don't think we need to close  close_all=.TRUE.
                  if(mds_save_server.ne.'LOCAL') then
                     socket = MdsSetSocket(-1)
                  endif
               endif
            else
C No: first time
               if(mds_server.ne.'LOCAL') connect=.true.
            endif
            mds_save_server=mds_server
         endif
C Scondary runs
      else

         if(.NOT.mds_cache_only) then
c  different server ?
            if(mds_server.ne.fdisk(6:il0)) then
               if(mds_server.eq.'LOCAL') then
                  socket = MdsSetSocket(-1)
                  socket=0
               else
                  socket=0
                  do i=1,naxxtra
                     if (ismds_x(i).gt.0 .and. i .ne. lrun_x ) then
                        il=index(fdisk_x(i),'@')-1
                        if (fdisk_x(i)(6:il).eq.mds_server) then 
                           socket = MdsSetSocket(-1)
                           socket = MdsSetSocket(ismds_x(i))
                           exit
                        endif
                     endif
                  enddo
                  if (socket .eq. 0) connect = .true.
               endif
C     same server
            else
               socket = ismds   ! same server...
               connect=.false.
            endif
            ismds_x(lrun_x) = socket
         endif
      endif
C
C  close trees on server, if necessary.
      if (close_all) then
         do i=1,mds_nclist
            stat = mds_close(mds_tree_clist(i),mds_shot_clist(i))
            if(mod(stat,2).ne.1) then
               call mdserr(luntrm,'?mds_connect(close): ',stat)
            endif
            mds_tree_clist(i)=' '
            mds_shot_clist(i)=0
         enddo
         mds_nclist=0
      endif
C

C  remote server?
      if(connect) then
         save_socket = MdsSetSocket(-1)
         il=index(mds_server,' ')
         mds_server(il:il)=CHAR(0)
         write(luntrm,*) 'Mds_connect: Connecting to server:  ',
     1        mds_server(1:il-1)
         socket = mdsfconnect(mds_server)
C if error?
         if (socket.lt.0) then
            write(luntrm,*) '?mdsfconnect:  error setting server:  ',
     >           mds_server(1:il)
            ier = 1
            if(lrun_x.eq.0) then
               mds_save_server=' '
               ismds=0
            else
               ismds_x(lrun_x)=0
            endif
            return
C Connected:
         else
            stat = MdsSetSocket(socket)
            if (lrun_x .eq. 0) then
               ismds=socket
            else
               ismds_x(lrun_x)=socket
            endif
         endif
      endif

C  open tree
      if(.NOT.mds_cache_only) then
         il=index(mds_tree,' ')-1
         if(il.le.0) il=len(mds_tree)
C
         if (lrun_x .eq. 0) then
            mdspulse=mds_shot
         else
            mdspulse_x(lrun_x)=mds_shot
         endif
         write(luntrm,*) 'Mds_connect: mds_open:  ',
     1        mds_tree(1:il),' ',mds_shot
         stat=mds_open(mds_tree(1:il),mds_shot)
         if(mod(stat,2).ne.1) then
            call mdserr(luntrm,'?mds_connect: mds_open:  ',stat)
            ier = 1
            return
         endif
      endif
C
C----------------------------------------------------------------------
      if(mds_cache.and.lcache) then
c
c  init cache; make sure cache directory exists
c
         if(lrun_x.eq.0) then
            call mds_cache_init(mds_server,mds_tree,mds_shot,
     >         tokyr,rid,mds_cache_dir,nbcache_act,nbcache_list,naxfxt)
           zcmd = 'mkdir -p '//trim(mds_cache_root)//trim(mds_cache_dir)
         else
            call mds_cache_init(mds_server,mds_tree,mds_shot,
     >         tokyr,rid,mds_cache_dir_x(lrun_x),
     >         nbcache_act_x(lrun_x),nbcache_list_x(1,lrun_x),naxfxt)
            zcmd = 'mkdir -p'//trim(mds_cache_root)//
     >                         trim(mds_cache_dir_x(lrun_x))
         endif
         call execute_command_line(zcmd,exitstat=iermds)
c
         if(iermds.ne.0) then
            write(luntrm,*) ' %mds_connect:  cache directory error'
            write(luntrm,*) '  mds_cache_root:  ',mds_cache_root
            if(lrun_x.eq.0) then
               write(luntrm,*) '  subdirectory:  ',mds_cache_dir
            else
               write(luntrm,*) '  subdirectory:  ',
     >            mds_cache_dir_x(lrun_x)
            endif
            write(luntrm,*) '  MDS+ caching *disabled*'
            mds_cache=.FALSE.
         endif
c
      endif
C
      return
      end
C--------------------------------------------------------------
