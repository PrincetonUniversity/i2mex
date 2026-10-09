      subroutine mds_cache_init(mds_server,mds_tree,mds_shot,
     >   tokyr,runid,cache_dir,ncache,ncache_list,listsize)
c
c  init cache data for new run
c
c  input:
c
      character*(*) mds_server          ! MDS+ server
      character*(*) mds_tree            ! MDS+ tree name
      integer mds_shot                  ! MDS+ shot, or, if zero...
      character*(*) tokyr               ! tok.yy directory id
      character*(*) runid               ! runid
c
c  output:
c
      character*(*) cache_dir           ! cache directory
      integer ncache                    ! cache binary items count (set to 0)
      integer ncache_list(listsize)     ! cache binary items list (cleared)
c
c----------------------------------------------
c
      character*20 zpre
      character*10 zshot
      character*1 zdelim
c
      integer str_length
c
c----------------------------------------------
c
c  subdirectory delimiter character
c
      zdelim='/'
      zpre='.rplot_cache'
c
      ilp=str_length(zpre)
      ils=str_length(mds_server)
      ilt=str_length(mds_tree)
c
      if(mds_shot.eq.0) then
         ily=str_length(tokyr)
         ilr=str_length(runid)
c
         cache_dir=zpre(1:ilp)//zdelim//
     >      mds_server(1:ils)//zdelim//mds_tree(1:ilt)//
     >      zdelim//tokyr(1:ily)//zdelim//runid(1:ilr)
c
      else
         write(zshot,'(i10)') mds_shot
         ilz=str_length(zshot)
         do i=1,ilz
            if(zshot(i:i).ne.' ') then
               inb=i
               exit
            endif
         enddo
c
         cache_dir=zpre(1:ilp)//zdelim//
     >      mds_server(1:ils)//zdelim//mds_tree(1:ilt)//
     >      zdelim//zshot(inb:ilz)
c
      endif
c
c  because of VMS syntax, change "." to "_" in server name part
c
      do ic=ilp+2,ilp+1+ils
         if(cache_dir(ic:ic).eq.'.') cache_dir(ic:ic)='_'
      enddo
c
c  clear cache list variables
c
      ncache=0
      do i=1,listsize
         ncache_list(i)=0
      enddo
c
      return
      end
