      subroutine mdstrd(ixflag,ier)
c
c  read the run's timebases.  if ixflag=0:  primary run
c  if ixflag.ne.0:  2ndary run
c
      use datmgr_mod
      use cplotr_mod
c
      implicit NONE
c
      integer ixflag,ier
c
      integer       Mds_Value
      integer       idescr_floatarr,idescr_long
      integer       lunzer
c
      integer idims(1), dsc, status,imiss,itest,inumt,inumr,iret
      character*120 expr
c
c----------------------------------------------
c
      ier=0
c
      if (ixflag .eq. 0) then

         imiss=1
         if(mds_cache) then
            call mds_cache_trd_size(ixflag,ntt,ntr,imiss)
         endif

         if(imiss.ne.0) then
            ! get timebase sizes from MDS+

            dsc = idescr_long(inumt)
            expr = 'size(.TRANSP_OUT:TIME1D)'
            status = Mds_Value(expr,dsc,NTT) ! read scalar timebase
            if (mod(status,2).ne.1) then
               expr = 'size(.OUTPUTS.ONE_D:TIME1D)'
               status = Mds_Value(expr,dsc,NTT) ! try again
               if (mod(status,2).ne.1) then
                  ier=1
                  call mdserr(lunzer(0), expr, status)
               endif
            endif

            dsc = idescr_long(inumr)
            expr = 'size(.TRANSP_OUT:TIME2D)'
            status = Mds_Value(expr,dsc,NTT) ! read scalar timebase
            if (mod(status,2).ne.1) then
               expr = 'size(.OUTPUTS.ONE_D:TIME2D)'
               status = Mds_Value(expr,dsc,NTT) ! try again
               if (mod(status,2).ne.1) then
                  ier=1
                  call mdserr(lunzer(0), expr, status)
               endif
            endif

            if(ier.ne.0) then
               write(lunzer(0),*) ' ?mdsTrd: timebase size read error.'
               ntt=0
               ntr=0
               return
            else
               ntt=inumt
               ntr=inumr
            endif
         endif

         itest=max(ntt,ntr)
         if(itest.gt.ntime) call dmg_texpand(itest)

         if(mds_cache) then
            call mds_cache_trd(ixflag,ntime,ntt,ntr,time,time3,imiss)
            if(imiss.eq.0) return

            if(mds_cache_only) then
               ier=1
               write(lunzer(0),*) ' ?mdstrd: data times not in cache.'
               write(lunzer(0),*) '  RPLOT_CACHE_ONLY = TRUE -> error.'
               return
            endif
         endif
            
         idims(1)= NTIME
         dsc = idescr_floatarr(time,idims,1)
         expr = '.TRANSP_OUT:TIME1D'
         status = Mds_Value(expr,dsc,NTT) ! read scalar timebase
         if (mod(status,2).ne.1) then
            expr = '.OUTPUTS.ONE_D:TIME1D'
            status = Mds_Value(expr,dsc,NTT) ! try again
            if (mod(status,2).ne.1) then
               ier=1
               call mdserr(lunzer(0), expr, status)
            endif
         endif
         print*,'NTT =',NTT,time(1),' - ',time(NTT)
         dsc = idescr_floatarr(time3,idims,1)
         expr = '.TRANSP_OUT:TIME2D'
         status = Mds_Value(expr,dsc,NTR) ! read profiles timebase
         if (mod(status,2).ne.1) then
            expr = '.OUTPUTS.ONE_D:TIME2D'
            status = Mds_Value(expr,dsc,NTR) ! try again
            if (mod(status,2).ne.1) then
               call mdserr(lunzer(0), expr, status)
               call copyr4(time,time3,NTT)
            endif
         endif
         print*,'NTR =',NTR,time3(1),' - ',time3(NTR)
C
         if(mds_cache) then
            call mds_cache_twr(ixflag,ntime,ntt,ntr,time,time3,ier)
         endif
C
      else
C
         itest=max(ntt_x(ixflag),ntr_x(ixflag))
         if(itest.gt.ntime) call dmg_texpand(itest)

         if(mds_cache) then
            call mds_cache_trd(ixflag,ntime,ntt_x(ixflag),ntr_x(ixflag),
     >           time_x(1,ixflag),time3_x(1,ixflag),imiss)
            if(imiss.eq.0) return

            if(mds_cache_only) then
               ier=1
               write(lunzer(0),*) ' ?mdstrd: data times not in cache.'
               write(lunzer(0),*) '  RPLOT_CACHE_ONLY = TRUE -> error.'
               return
            endif
         endif
            
         idims(1)= NTIME
         dsc = idescr_floatarr(time_x(1,ixflag),idims,1)
         expr = '.TRANSP_OUT:TIME1D'
         status = Mds_Value(expr,dsc,NTT_X(ixflag)) ! read scalars timebase
         if (mod(status,2).ne.1) then
            expr = '.OUTPUTS.ONE_D:TIME1D'
            status = Mds_Value(expr,dsc,NTT_X(ixflag)) ! try again
            if (mod(status,2).ne.1) then
               ier=1
               call mdserr(lunzer(0), expr, status)
            endif
         endif
         print*,'NTT_X(',ixflag,') =',NTT_X(ixflag),time_x(1,ixflag),
     >        ' - ', time_x(NTT_X(ixflag),ixflag)
         dsc = idescr_floatarr(time3_x(1,ixflag),idims,1)
         expr = '.TRANSP_OUT:TIME2D'
         status = Mds_Value(expr,dsc,NTR_X(ixflag)) ! read profiles timebase
         if (mod(status,2).ne.1) then
            expr = '.OUTPUTS.ONE_D:TIME2D'
            status = Mds_Value(expr,dsc,NTR_X(ixflag)) ! try again
            if (mod(status,2).ne.1) then
               call mdserr(lunzer(0), expr, status)
               call copyr4(time_x(1,ixflag),time3_x(1,ixflag),
     >              NTT_X(ixflag))
            endif
         endif
         print*,'NTR_X(',ixflag,') =',NTR_X(ixflag),time3_x(1,ixflag),
     >        ' - ',time3_x(NTR_X(ixflag),ixflag)
C
         if(mds_cache) then
            call mds_cache_twr(ixflag,ntime,ntt_x(ixflag),ntr_x(ixflag),
     >           time_x(1,ixflag),time3_x(1,ixflag),ier)
         endif
C
      endif
c
      return
      end
