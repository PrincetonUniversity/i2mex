      logical function mds_cache_active(idum)
c
c  return TRUE if the MDS_CACHE flag is set, FALSE otherwise
c
      use cplotr_mod
c
      mds_cache_active = mds_cache
c
      return
      end
c
c-------------------------------------------------------------------
c
      logical function mds_cache_exclusive(idum)
c
c  return TRUE if all MDS data to be read through CACHE -- i.e. no
c  cache validation; cache miss is an error.
c
      use cplotr_mod
c
      mds_cache_exclusive = mds_cache .AND. mds_cache_only
c
      return
      end
