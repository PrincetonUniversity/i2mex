#define BUFSIZ 300

      subroutine gmkdir(root_in,rpath_in,ier)
c
c  make a multi-level subdirectory starting from root
c
c  example:  call gmkdir('HOME','.cache/sub1/sub2/sub3',ier)
c    => $HOME/.cache/sub1/sub2/sub3 directory is created (if it does
c       not already exist).  Each higher level directory is also created
c       as needed.
c
      implicit NONE
c
c  input:
c
      character*(*) root_in             ! root -- start point (must exist)
      character*(*) rpath_in            ! subdirectory chain
c
c  output:
      integer ier                       ! error code, 0 = OK
c
c-------------------------------------------------------------
c
      character*BUFSIZ start
      character*BUFSIZ next
c
      character*250 root,rpath
c
c  DMC -- procedure to check directories /u & /p can be unreliable
c  and is not needed.  If ' ' is passed for root, start check a 
c  little further down the path than "/".
c
      integer, parameter :: iroot_param_srch = 6
      integer :: ii,istat
c
      character*512 zcmd
c-------------------------------------------------------------
c
      if(root_in.ne.' ') then
         root = root_in
         rpath = rpath_in
      else
         if(rpath_in(1:1).eq.'/') then
            do ii=iroot_param_srch,1,-1
               if(rpath_in(ii:ii).eq.'/') then
                  root=rpath_in(1:ii)
                  rpath=rpath_in(ii+1:)
                  exit
               endif
            enddo
         else if(rpath_in(1:1).eq.'~') then
            if(rpath_in(2:2).eq.'/') then
               root='HOME'
               rpath=rpath_in(3:)
            else
               ! (I don't think this can work: ~<username>)
               root=' '
               rpath = rpath_in
            endif
         else
            call getcwd(root)
            rpath = rpath_in    ! path relative to $cwd
         endif
      endif

c!      write(0,*) ' ** root_in = "',trim(root_in),'"'
c!      write(0,*) ' ** rpath_in = "',trim(rpath_in),'"'
c!      write(0,*) ' ** root = "',trim(root),'"'
c!      write(0,*) ' ** rpath = "',trim(rpath),'"'

      ier=0
      call gmkdir_init(root,start,rpath,next)

c!      write(0,*) ' gmkdir_init: '
c!      write(0,*) '   root = "',trim(root),'"'
c!      write(0,*) '   start = "',trim(start),'"'
c!      write(0,*) '   rpath = "',trim(rpath),'"'
c!      write(0,*) '   next = "',trim(next),'"'

 10   continue
      if((ier.eq.0).and.(next.ne.' ')) then
         call gmkdir_next(start,next,ier)

c!!         write(0,*) ' gmkdir_next: '
c!!         write(0,*) '   start = "',trim(start),'"'
c!!         write(0,*) '   next = "',trim(next),'"'
c!!         write(0,*) '   ier = ',ier

         go to 10
      endif
c
      return
      end
c
c-------------------------------------------------------------
c
      subroutine gmkdir_init(root,start,rpath,next)
c
      implicit NONE
c
      character*(*) root                ! root (as passed)
      character*(*) start               ! starting directory (standardized)
c
      character*(*) rpath               ! relative path (as passed)
      character*(*) next                ! starting rel. path (standardized)
c-------------
      integer :: iclen,ic,ic1,ic2,istat,iplen,io
c
      character*BUFSIZ ztest
      character*15 cpid
      character*30 tmpfile
c---------------------------------
c
      start=' '
      next=' '
c
c  root directory...
c
      call sget_pid_str(cpid,iplen)
      tmpfile = 'gmkdir_'//cpid(1:iplen)//'.tmp'
c
      if(root.ne.' ') call ufilnam(root,' ',start)
c
c  relative path from there...
c
      iclen=len_trim(rpath)
c
c  does requested directory already exist?
c
      call find_io_unit(io)
c
      ztest = trim(start)//'/'//trim(rpath)//'/'//trim(tmpfile)
      open(unit=io,file=trim(ztest),status='unknown',iostat=istat)
      if(istat.eq.0) then
                                !  open successful; directory must exist.
         close(unit=io,status='delete',iostat=istat)
         start = trim(start)//'/'//trim(rpath)
         next = ' '
         return
      endif
c
c  work up through rpath parents...
c
      ic=iclen
      do
         ic=ic-1
         if(ic.le.0) exit
         if(rpath(ic:ic).ne.'/') cycle

         ztest = trim(start)//'/'//rpath(1:ic)//'/'//trim(tmpfile)
         open(unit=io,file=trim(ztest),status='unknown',iostat=istat)
         if(istat.eq.0) then
                                !  open successful; directory must exist.
            close(unit=io,status='delete',iostat=istat)
            start = trim(start)//'/'//rpath(1:ic)
            next = rpath(ic+1:)
            return
         endif
      enddo
c
c  try start directory
c
      ztest = trim(start)//'/'//trim(tmpfile)
      open(unit=io,file=trim(ztest),status='unknown',iostat=istat)
      if(istat.eq.0) then
                                !  open successful; directory must exist.
         close(unit=io,status='delete',iostat=istat)
         next=rpath
         return
      endif
c
c  this code reached if directory path (start) is not fully in place.
c  (this would be unusual).
c
      ic=0
 10   continue
      ic=ic+1
      if(ic.gt.iclen) go to 1000        ! reached the end
      if((ichar(rpath(ic:ic)).ne.0).and.(rpath(ic:ic).ne.' ').and.
     >   (rpath(ic:ic).ne.'~').and.(rpath(ic:ic).ne.'/')) then
         ic1=ic
      else
         go to 10
      endif
c
      ic=iclen+1
 20   continue
      ic=ic-1
      if(rpath(ic:ic).ne.'/') then
         ic2=ic
      else
         go to 20
      endif
c
      next=rpath(ic1:ic2)
      go to 1000
c-------------------
 1000 continue
      return
      end
c----------------------------------------------------------------------
c
      subroutine gmkdir_next(start,next,ier)
      use iso_c_binding, only: c_char
c
c check the next directory
c
      implicit NONE
c
      character*BUFSIZ start            ! current root
      character*BUFSIZ next             ! current relative path
c
      integer ier
c
c-----------------
c
      integer cmkdir
      integer :: ils,ilx,inext,ild,iln,istat
c
      character*1 zdelim,zterm
      character*BUFSIZ dirfile,newdir
c
      character*512 zcmd
      character(kind=c_char), dimension(:), allocatable :: cpath
c-----------------
c
      zdelim='/'
c
      ils=len_trim(start)
      ilx=len_trim(next)
c
      inext=index(next,zdelim)-1
      if(inext.le.0) then
         inext=ilx
      endif
c
      if(start.ne.' ') then
         dirfile=start(1:ils)//next(1:inext)
      else
c
c  blank root
c
         dirfile=next(1:inext)
      endif
      ild=len_trim(dirfile)
      newdir=dirfile(1:ild)//'/'
      ild=ild+1
      iln=len_trim(newdir)
c-------------------------------------------
c make the shell command
c   ...mod DMC: call cmkdir(<path>)
      if(allocated(cpath)) deallocate(cpath)
      allocate(cpath(ild+1))
      call cstring(newdir(1:ild),cpath,'2C')
      istat = cmkdir(cpath)
c
      ier=istat
c--------------------------------------------
c  new root and rel. path
c
      start=newdir
      if(inext.lt.ilx) then
         next=next(inext+2:ilx)
      else
         next=' '
      endif
c
      return
      end
 
