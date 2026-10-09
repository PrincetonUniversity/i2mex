      subroutine mdsmfrd(ifcn, ind, retn)
C
C Read "ifcn" Profile Function into COMMON buffer
C
C 01/06/00  CAL
 
      use datmgr_mod
      use cplotr_mod
 
C Input
      integer ifcn
      integer ind
C Return
      integer retn
 
C Processing
      integer       ipt, isize
      integer       Mds_Value
      character*60  expr
      integer       status, size
      integer       idims(2), dsc
      character*21  name
C----------------------------------------------------------------------
 
      ipt = locd(ind)
C
      if(lrun_x.eq.0) then
         isize=nzonex(itypr(ifcn)) * ntr
         expr='.TRANSP_OUT:' // abr(ifcn)
         idims(1)=nzonex(itypr(ifcn))
         idims(2)=ntr
         name=abr(ifcn)
      else
         itype=itypr_x(ifcn,lrun_x)
         izonex=nzonex_x(itype,lrun_x)
         isize=izonex*ntr_x(lrun_x)
         idims(1)=izonex
         idims(2)=ntr_x(lrun_x)
         expr='.TRANSP_OUT:' // abr_x(ifcn,lrun_x)
         name=abr_x(ifcn,lrun_x)
      endif
C
      nwds(ind) = isize
c
c  try the cache...
c
      imiss=1
      if(mds_cache) then
         call mds_cache_read(lrun_x,name,ifcn,datbuf(ipt),isize,imiss)
      endif
c
      if(imiss.gt.0) then
         if(mds_cache_only) then
            write(lunzer(0),*) ' RPLOT_CACHE_ONLY = TRUE -> error.'
            write(lunzer(0),*) ' failed to read f(t) data in cache: ',
     >           trim(name)
            retn=imiss
            return
         endif
c
c  cache miss, read over network with MDS+
c
         write(lunzer(0),*)
     >      ' %mdsmfrd:  cache miss, reverting to MDS+: ',name
         dsc = idescr_floatarr(datbuf(ipt),idims,2)
         status = Mds_Value(expr,dsc,size)
         if (mod(status,2).ne.1) then
            luntrm=lunzer(0)
            call mdserr(luntrm, expr, status)
            if (status .ne. 265388258) then
               retn = 1
               return
            else
               write(luntrm,*) expr,' empty:  taken as ZERO'
               datbuf(ipt:ipt+isize-1)=0.0
            endif
         endif
         if ( size .le. 1 ) then
            write(lunzer(0),*)
     >      ' ?mdsmfrd:  got no data from mdsplus - error'
            write(lunzer(0),*)
     >      '            check if $ARCDIR is accessible'
            retn=1
            return
         endif
C
C  write cache entry
C
         if(mds_cache.and.(imiss.eq.1)) then
            call mds_cache_write(lrun_x,name,ifcn,datbuf(ipt),isize,ier)
            if(ier.ne.0) 
     >           write(luntrm,*) ' %cache record error -- ignorable.'
         endif
C
      endif
C
      return
      end
 
 
 
