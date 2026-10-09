      subroutine mds_cache_read(ixflag,zname,ifnum,zdata,isize,imiss)
c
      use cplotr_mod
c
c  MDS+ cache read attempt
c
c  input:
c    ixflag -- primary / 2ndary run index code
c    zname -- name of item / binary file
c    ifnum -- id number of item
c                -1 for f(t) dataset
c                -2 for scalar f(t) timebase
c                -3 for profile f(x,t) timebase
c
c  output:
c    zdata(1:isize) -- the data to read
c    imiss -- =1 on exit if read fails, =0 if successful
c             =2: item size too large for cache
c
      integer ixflag,ifnum,isize,imiss
      character*(*) zname
      real zdata(isize)
c
c-----------------
c
      character*200 zfiln
c
      integer str_length
c
      integer :: ilim_pbin,ilim_size,lunzer
c
      character*3 :: encod
c-----------------
c
      call mds_cache_lims(ilim_pbin,ilim_size)
      if(isize.gt.ilim_size) then
         imiss=2
         write(lunzer(0),*) 
     >        ' ?mds_cache_write: single item size too large: ',isize
         return
      endif
c
c  first:  is item in the list?
c
      imiss=1
      if(ixflag.eq.0) then
         do i=1,nbcache_act
            if(ifnum.eq.nbcache_list(i)) imiss=0
         enddo
      else
         do i=1,nbcache_act_x(ixflag)
            if(ifnum.eq.nbcache_list_x(i,ixflag)) imiss=0
         enddo
      endif
c
      if(imiss.eq.1) return
c
c  item is on list; data file should be available
c
      ilz=str_length(zname)
      call mds_cache_fname(ixflag,zname(1:ilz)//'.DAT',zfiln,ier)
      if(ier.ne.0) then
         imiss=1
         return
      endif
c
c  try to read the data
c
      ilf=str_length(zfiln)
c
      if(isize.le.ilim_pbin) then
         call pbinrd(lunzer(0),lun_tf,zfiln(1:ilf),zdata,isize,igot,
     >        ier)
         if(ier.ne.0) then
            imiss=1
            return
         endif
         if(igot.ne.isize) then
            write(lunzer(0),*)
     >           ' ?mds_cache_read:  unexpected size, item:  ',
     >           zname(1:ilz)
            imiss=1
            return
         endif

      else
         ia1 = 1-ilim_pbin
         ict = -1
         do 
            ia1 = ia1 + ilim_pbin
            if(ia1.gt.isize) exit

            ia2 = min(isize,(ia1+ilim_pbin-1))
            isiz2 = ia2-ia1+1

            ict = ict + 1
            write(encod,'("p",i2.2)') ict

            call pbinrd(lunzer(0),lun_tf,zfiln(1:ilf)//encod,
     >           zdata(ia1),isiz2,igot,ier)
            if(ier.ne.0) exit
            if(igot.ne.isiz2) then
               ier=1
               write(lunzer(0),*)
     >              ' ?mds_cache_read:  unexpected size, item:  ',
     >              zname(1:ilz),' page: '//encod
               exit
            endif
         enddo
         if(ier.ne.0) imiss=1
      endif
c
c  ok
c
      return
      end
