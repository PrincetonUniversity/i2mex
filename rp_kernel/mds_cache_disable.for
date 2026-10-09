      subroutine mds_cache_disable
c
c  disable MDS+ cache
c
      use cplotr_mod
c
      write(6,*) '[mds_cache_disable:  MDS+ cache disabled.]'
      mds_cache=.false.
      mds_cache_only=.false.
c
      return
      end
