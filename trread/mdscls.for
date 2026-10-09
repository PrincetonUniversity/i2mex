      subroutine mdscls(idrun,ier)
C
C 01/28/00 CAL : close a MDSplus pulse
C
      use cplotr_mod
      integer       idrun      ! index into secondary run(s)
      integer       ier        ! returned 0 = success
C                              !          1 = error
#ifndef __NOMDSPLUS
C
C Local
      character*64  mds_server
      character*24  mds_tree,mds_tree_orig
      integer       mds_shot,mds_shot_orig,stat
      integer       il, luntrm, iflag
C
 
      luntrm=lunzer(0)
 
C Parse tree: main run
      il=index(fdisk,'@')-1
      mds_tree_orig=fdisk(il+2:)
      mds_shot_orig=mdspulse
C
C this run:
      if(idrun.eq.0) then
         mds_tree=mds_tree_orig
         mds_shot=mds_shot_orig
      else
         il=index(fdisk_x(idrun),'@')-1
         mds_tree=fdisk_x(idrun)(il+2:)
         mds_shot=mdspulse_x(idrun)
      endif
C close pulse -- or, if on main run server, just re-open main run
C  and defer close (for performance reasons...
      il=index(mds_tree,' ')-1
      if(il.le.0) ilt=len(mds_tree)
      iflag=0
      if((idrun.gt.0).and.(ismds.gt.0)) then
C same server
         if(ismds_x(idrun).eq.ismds) then
            iflag=1
            call mds_close_defer(mds_tree(1:il),mds_shot,stat)
            print *,' ...defer close'
            stat=mds_open(trim(mds_tree_orig),mds_shot_orig)
C different server: Set context of primary
         else if(ismds_x(idrun).gt. 0) then
            call mds_close_defer(mds_tree(1:il),mds_shot,stat)
            print *,' ...defer close'
            stat = MdsSetSocket(ismds)
            if(stat.eq.0) then
               write(luntrm,*) '?mdscls -- primary run server ',
     >              'reconnect failed!'
            else
               stat=mds_open(trim(mds_tree_orig),mds_shot_orig)
               iflag=1
            endif
         endif
      endif
      if(iflag.eq.0) then
         !  different server -- do the close now.
         stat=mds_close(mds_tree(1:il),mds_shot)
      endif
      if (mod(stat,2).ne.1) then
         call mdserr(luntrm,'%mdscls:  ',stat)
         ierr = 1
      else
         ierr = 0
      endif
C
C  02/23/04  CAL: Should never get here
C-----------------------
C  if necessary:  disconnect from 2ndary run server, reconnect to
C  primary run server.  If this happens, above logic also assured that
C  the tree was really closed.
C
      if(iflag.eq.0 .and. idrun.gt.0) then
         if(ismds.gt.0) then
            if((ismds_x(idrun).gt.0).and.
     >         (ismds_x(idrun).ne.ismds)) then

               ! sts = MdsDisconnect()
               mds_server="!"
               mds_server(2:2)=char(0)
               sts = MdsCacheConnect(mds_server)  ! disconnects all cached servers

               il=index(fdisk,'@')-1
               mds_server=fdisk(6:il)
               il=index(mds_server,' ')
               mds_server(il:il)=char(0)
               sts = MdsCacheConnect(mds_server)
C               ismds = MdsSetSocket(sts)
               if(sts.le.0) then
                  write(luntrm,*) '?mdscls -- primary run server ',
     >               'reconnect failed!'
                  ierr=1
               endif
            endif
         endif
C
      endif
C
#endif
      return
      end
C----------------------------------------------------------------------
      subroutine mds_close_defer(mds_tree,mds_shot,istat)
C
      use cplotr_mod
c
c  add tree to deferred-close list
c
      character*(*), intent(in) :: mds_tree ! tree name: deferred close
      integer, intent(in) :: mds_shot   ! tree id: deferred close.
      integer, intent(out) :: istat     ! status code.
#ifndef __NOMDSPLUS
c
c-------------------------
c
      integer i,imatch,luntrm,mds_close
c
c-------------------------
c  do not add if it is already in the list.
c  if the list is full, close the tree now.
c
      istat=0
c
      luntrm=lunzer(0)
c
      imatch=0
      do i=1,mds_nclist
         if((mds_tree_clist(i).eq.mds_tree) .and.
     >      (mds_shot_clist(i).eq.mds_shot)) then
            imatch=i
            exit
         endif
      enddo
c
c  exit now if tree already on list...
      if(imatch.gt.0) return
c
c  check list...
c
      if(mds_nclist.eq.mds_nclist_max) then
c
c  list is full
c
         write(luntrm,*)
     >      ' %mds_close_defer: MDSplus deferred close list full.'
         istat=mds_close(trim(mds_tree),mds_shot)
         if (mod(istat,2).ne.1) then
            call mdserr(luntrm,'%mds_close_defer:  ',istat)
         endif
      else
c
c  list not full: add to list
c
         mds_nclist=mds_nclist+1
         mds_tree_clist(mds_nclist) = mds_tree
         mds_shot_clist(mds_nclist) = mds_shot
c
      endif
c
#endif
      return
      end
