      subroutine mds_cache_write(ixflag,zname,ifnum,zdata,isize,ier)
c
      use cplotr_mod
c
c  MDS+ cache write attempt
c
c  input:
c    ixflag -- primary / 2ndary run index code
c    zname -- name of item / binary file
c    zdata(1:isize) -- the data to write
c    ifnum -- id number of item
c
c  output:
c    ier -- completion code, 0=normal, successful
c
      integer ixflag,ifnum,isize
      character*(*) zname
      real zdata(isize)
c
c-----------------
c
      character*200 zfiln
c
      integer str_length
c
      integer :: ilim_pbin,ilim_size,isiz2
      character*3 :: encod
c-----------------
c
c  first:  is item in the list?
c
      ier=0
c
      call mds_cache_lims(ilim_pbin,ilim_size)
      if(isize.gt.ilim_size) then
         ier=1
         write(lunzer(0),*) 
     >        ' ?mds_cache_write: single item size too large: ',isize
         return
      endif
c
      item=0
      if(ixflag.eq.0) then
         inum=nbcache_act
         do i=1,inum
            if(ifnum.eq.nbcache_list(i)) item=i
         enddo
      else
         inum=nbcache_act_x(ixflag)
         do i=1,inum
            if(ifnum.eq.nbcache_list_x(i,ixflag)) item=i
         enddo
      endif
c
c  is the cache full?
c
      if(item.eq.0) then
         if(inum.ge.nbcache_lim) return
      endif
c
c  cache not full or item should be in cache, so, write...
c
      ilz=str_length(zname)
      call mds_cache_fname(ixflag,zname(1:ilz)//'.DAT',zfiln,ier)
      if(ier.ne.0) then
         ier=1
         return
      endif
c
c  try to write the data
c
      ilf=str_length(zfiln)
c
      if(isize.le.ilim_pbin) then
         call pbinwr(lunzer(0),lun_tf,zfiln(1:ilf),zdata,isize,ier)
         if(ier.ne.0) then
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

            call pbinwr(lunzer(0),lun_tf,zfiln(1:ilf)//encod,
     >           zdata(ia1),isiz2,ier)
            if(ier.ne.0) exit
         enddo
         if(ier.ne.0) return
      endif
c
c  ok:  update internal data & index file if necessary...
c
      if(item.eq.0) then
         if(ixflag.eq.0) then
            nbcache_act=nbcache_act+1
            nbcache_list(nbcache_act)=ifnum
         else
            nbcache_act_x(ixflag)=nbcache_act_x(ixflag)+1
            nbcache_list_x(nbcache_act_x(ixflag),ixflag)=ifnum
         endif
c
         call mds_cache_fname(ixflag,'ITEMS.LIST',zfiln,ier)
         if(ier.ne.0) return
c
         ilf=str_length(zfiln)
         open(unit=lun_tf,file=zfiln(1:ilf),status='OLD',
#if __F90
     >      position='APPEND',
#else
     >      access='APPEND',
#endif
     >      iostat=ier)
         if(ier.eq.0) then
            write(lun_tf,'(1x,a,t22,i9)') zname(1:ilz),isize
            close(unit=lun_tf)
         endif
c
      endif
c
      return
      end

      subroutine mds_cache_lims(ilim1,ilim100)

c  return single item cache size limits (determined by pbinwr I/O routine)

      integer, intent(out) :: ilim1,ilim100

      ilim1 = 256*256*256 - 1
      ilim100 = 100*ilim1

      return
      end
