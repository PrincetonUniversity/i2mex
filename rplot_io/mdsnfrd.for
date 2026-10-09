      subroutine mdsnfrd(ipt, isize, retn)
C
C Read Scalar Functions into COMMON buffer
C
C 01/06/00  CAL
 
      use datmgr_mod
      use cplotr_mod
 
C Input
      integer ipt           ! pointer into datbuf
C Return
      integer isize         ! size of all Scalar Functions
      integer retn
 
C Processing
      integer       i, ipt1, ier
      integer       Mds_Value
      character*60  expr
      integer       status, size
      integer       idims(1), dsc
      integer       imiss
C----------------------------------------------------------------------
C
C Read all Scalar Functions into datbuf
C--------------------------------------
 
 
cdbg      type *,'Read Scalar Functions',nft,ntt
 
      luntrm=lunzer(0)
      ier=0
 
      if(lrun_x.eq.0) then
         isize=ntt*nft
         inft=nft
         intt=ntt
      else
         isize=ntt_x(lrun_x)*nft_x(lrun_x)
         inft=nft_x(lrun_x)
         intt=ntt_x(lrun_x)
      endif
c
      imiss=1
      if(mds_cache) then
         call mds_cache_read(lrun_x,'RUN_SCALARS',-1,datbuf(ipt),isize,
     >      imiss)
      endif
c
      if(imiss.eq.1) then
         if(mds_cache_only) then
            write(lunzer(0),*) ' RPLOT_CACHE_ONLY = TRUE -> error.'
            write(lunzer(0),*) ' failed to read f(t) data in cache.'
            retn=imiss
            return
         endif
c
c  read over network via MDS+
c
         write(lunzer(0),*) ' ...reading TRANSP f(t) scalars with MDS+'
         write(lunzer(0),*) '    (this may take a while)'
         ipt1 = ipt
         idims(1)=intt
         do i = 1, inft
            dsc = idescr_floatarr(datbuf(ipt1),idims,1)
            if(lrun_x.eq.0) then
               expr = '.TRANSP_OUT:'//abt(i)
            else
               expr = '.TRANSP_OUT:'//abt_x(i,lrun_x)
            endif
            status = Mds_Value(expr,dsc,size)
            if (mod(status,2).ne.1) then
               call mdserr(luntrm, expr, status)
C Don't stop if TREE$-E-NODATA, no data available for this node
               if (status .ne. 265388258) then
                  retn = 1
                  return
               else
                  write(luntrm,*) expr,' empty:  taken as ZERO'
                  datbuf(ipt1:ipt1+intt-1)=0.0
                  ipt1 = ipt1 + intt
                  ier=ier+1
               endif
            else
               if ( size .le. 1 ) then
                  write(lunzer(0),*)
     >                 '?mdsnfrd:  got no data from mdsplus - error'
                  write(lunzer(0),*)
     >                 '           check if $ARCDIR is accessible'
                  retn=1
                  return
               endif
               ipt1 = ipt1 + intt
            endif
         end do
         if (ier .eq. 0) then
            write(lunzer(0),*)
     >      ' ...TRANSP f(t) scalars read completed.'
         else
            write(lunzer(0),900) ier
 900        format(' ...TRANSP f(t) scalars read completed with',
     >      i4,' warnings')
         endif
c
c  write cache entry
c
         if(mds_cache) then
            call mds_cache_write(lrun_x,'RUN_SCALARS',-1,
     >         datbuf(ipt),isize,ier)
         endif
c
      endif
c
      retn = 0
      return
c
      end
