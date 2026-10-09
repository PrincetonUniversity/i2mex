      subroutine mds_cache_list(ixflag,zfiln,ildir,cache_list,inum,imax,
     >     ier)
c
      use cplotr_mod
c
c  read cache list file; or, if not found:  create empty cache list file.
c
c  input:
      integer ixflag                    ! main run / 2ndary run index
      character*(*) zfiln               ! cache list filename
      integer ildir                     ! length of path part of filename
c
      integer imax                      ! max list size
c
c  output:
      integer inum                      ! current actual list size
      character*(*) cache_list(imax)    ! the list...
      integer ier                       ! completion code, 0=OK
c
c  list items:  -1 => scalar function data
c               +n => profile function #n
c
c---------------------------
c
      integer str_length,idescr_cstring
c
      character*20 zname
c
      character*40 zdate
      character*40 ztest
c---------------------------
c
c  validation algorithm -- version string
c  change this to invalidate all file caches, e.g. because of caching
c  format change or bugfixes
c
      character*4, parameter :: cc_vsn = 'v1.0'
      character*4 cc_test
c
c  dmc: additional validation: non-zero timebases
c
      real :: tdum1(1),tdum2(1)
      integer :: intt,intr
c
c---------------------------
c  get MDS+ run date (for cache validation)
c
      zdate=' '
      ztest=' '
      zdate='???'
c
      if(.NOT.mds_cache_only) then
         istat = Mds_Value('DATE_RUN',idescr_cstring(zdate),isize)
         if (mod(istat,2).ne.1) then
            call mdserr(lun, ':DATE_RUN', istat)
            ier=1
            write(lunzer(0),*)
     >      ' ?mds_cache_list: DATE_RUN not valid; will try DATE_BEGAN!'
            istat = Mds_Value('DATE_BEGAN',idescr_cstring(zdate),isize)
            if (mod(istat,2).ne.1) then
               call mdserr(lun, ':DATE_BEGAN', istat)
               ier=1
               write(lunzer(0),*)
     >         ' ?mds_cache_list:  cache validation data access failed!'
               return
            endif
         endif
      endif
c
c  look for timebase sizes (used in validation procedure)
c  (if error occurs, the timebase sizes will be zero, and the cache
c  will not be validated).
c
      call mds_cache_trd(ixflag,0,intt,intr,tdum1,tdum2,ier)
c
c  try first to open for read...
c
      ier=0
      inum=0
      ilf=str_length(zfiln)
c
      if(min(intt,intr).eq.0) then
         write(lunzer(0),*) ' %mds_cache_list: CACHE data invalid'
         write(lunzer(0),*) '  (zero length time vector detected)'
         ios=1
      else
         open(unit=lun_tf,file=zfiln(1:ilf),status='OLD',
     >        iostat=ios)
      endif
c
 5    continue
c
      if(ios.ne.0) then
c
c  open for read failed, so, start new cache directory...
c
         call fclean_dir(zfiln(1:ildir),ier)
c
         if(mds_cache_only) then
            write(lunzer(0),*) ' ?mds_cache_list: CACHE data not found.'
            write(lunzer(0),*) '  RPLOT_CACHE_ONLY = TRUE -> error.'
            ier=1
            return
         endif
c
         open(unit=lun_tf,file=zfiln(1:ilf),status='NEW',
     >      iostat=ios)
c
         if(ios.ne.0) then
            write(lunzer(0),*)
     >         ' ?mds_cache_list:  item list file create failed.'
            ier=1
         else
            ilzd=str_length(zdate)
            write(lun_tf,'(1x,''%validation:         '',a4)') cc_vsn
            write(lun_tf,'(1x,a)') zdate(1:ilzd)
            close(unit=lun_tf)
            ier=0
         endif
         return
      endif
c
c  read list
c
      icount=0
 10   continue
      icount=icount+1
      read(lun_tf,'(1x,a,t23,a4)',end=99) zname,cc_test
      if(icount.eq.1) then
c
c  cache entry validation
c
         ios=0
         if(zname(1:1).eq.'%') then
            read(lun_tf,'(1x,a)') ztest

            ! skip run date test if mds_cache_only is set, but, still
            ! check the cache software version ID

            if(.NOT.mds_cache_only) then
               if(ztest.ne.zdate) then
                  ios=1
                  write(lunzer(0),*)
     >                 ' %MDS+ cache invalid (rundate change)'
               endif
            endif

            if(cc_test.ne.cc_vsn) then
               ios=1
               write(lunzer(0),*)
     >            ' %MDS+ cache invalid (code version update: ',
     >            cc_vsn,')'
            endif
         else
            ios=1
               write(lunzer(0),*)
     >            ' %MDS+ cache invalid (list file syntax error)'
         endif
         if(ios.ne.0) then
            ! write accessible cache: action is to flush on validation failure:
            write(lunzer(0),*)
     >           ' %MDS+ cache flushed:  validation failure.'
            close(unit=lun_tf)
            go to 5
         else
            go to 10                    ! OK:  keep reading
         endif
      endif
c
      if(zname(1:1).eq.'!') go to 10
c
      if(inum.ge.imax) go to 10
c
      inum=inum+1
      cache_list(inum)=zname
c
      go to 10
c
 99   continue
      close(unit=lun_tf)
c
      return
      end
c-----------------------------
      subroutine mds_cache_list_convert(ixflag,cache_list,ilc,
     >   inum,ilist,ier)
c
      use cplotr_mod
c
c  convert list of names to numeric form
c
      implicit NONE
c
c  input
c
      integer ixflag                    ! main run / 2ndary run index
      integer ilc                       ! number of names in list
      character*(*) cache_list(ilc)     ! items in cache list
c
c  output
c
      integer inum                      ! list length (copied)
      integer ilist(ilc)                ! numeric list
      integer ier                       ! error code
c
c--------------------------------
c
      integer iordr_tmp(naxfxt)
c
      character*20 zname
c
      integer :: i,indx
      integer :: ifind_ordr
c--------------------------------
c
      ier=0
      if(ilc.eq.0) then
         inum=0
         return
      endif
c
      if(ixflag.ne.0) then
         call aordr(iordr_tmp,abr_x(1,ixflag),nfxt_x(ixflag))
      else
         call aordr(iordr_tmp,abr,nfxt)
      endif
c
      inum=0
      do i=1,ilc
         zname=cache_list(i)
         if(zname.eq.'RUN_SCALARS') then
c
            inum=inum+1
            ilist(inum)=-1
c
         else if(zname.eq.'time1d') then
c
            inum=inum+1
            ilist(inum)=-2
c
         else if(zname.eq.'time2d') then
c
            inum=inum+1
            ilist(inum)=-3
c
         else
c
            pltabb=zname
            if(ixflag.eq.0) then
               indx=ifind_ordr(abr,iordr_tmp,nfxt,pltabb)
            else
               indx=ifind_ordr(abr_x(1,ixflag),iordr_tmp,nfxt_x(ixflag),
     >            pltabb)
            endif
c
            inum=inum+1
            ilist(inum)=indx
c
         endif
      enddo
c
      return
      end
c-----------------------------
      subroutine mds_cache_list_precon(cache_list,ilc,
     >   inum,ilist,ier)
c
      use cplotr_mod
c
      implicit NONE
c
c  convert list of names to numeric form
c
c  input
c
      integer ilc                       ! number of names in list
      character*(*) cache_list(ilc)     ! items in cache list
c
c  output
c
      integer inum                      ! list length (copied)
      integer ilist(ilc)                ! numeric list
      integer ier                       ! error code
c
c--------------------------------
c
      integer iordr_tmp(naxfxt)
c
      character*20 zname
c
      integer :: i
c--------------------------------
c
      ier=0
      if(ilc.eq.0) then
         inum=0
         return
      endif
c
      inum=0
      do i=1,ilc
         zname=cache_list(i)
         if(zname.eq.'RUN_SCALARS') then
c
            inum=inum+1
            ilist(inum)=-1
c
         else if(zname.eq.'time1d') then
c
            inum=inum+1
            ilist(inum)=-2
c
         else if(zname.eq.'time2d') then
c
            inum=inum+1
            ilist(inum)=-3
c
         else
c
            inum=inum+1
c
         endif
      enddo
c
      return
      end
