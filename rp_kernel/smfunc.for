      subroutine smfunc(indt,zdelt,indx,zdelx,
     >   typ_epst,zepst,typ_epsx,zepsx,
     >   istat,ipf,ier)

      use datmgr_mod
      use cplotr_mod

C  non-interactive smoothing routine-- implement calculator %SMOOTH
C  command.
C
C  inputs:
C
      integer indt                      ! "time indicial" flag (if =1)
      real zdelt                        ! time delta
      integer indx                      ! "x axis indicial" flag (if =1)
      real zdelx                        ! x delta
C
      character*1 typ_epst              ! epsilon type (R, A, or %)
      real zepst                        ! epsilon -- max change on time smooth
      character*1 typ_epsx              ! epsilon type (R, A, or %)
      real zepsx                        ! epsilon -- max change on x smooth
C
      integer istat                     ! type of profile being smoothed
C
      integer ipf                       ! location of data to be smoothed
C
C  output:
C
      integer ier                       ! completion code, 0 = normal
C
C     DATBUF(ipf...) are modified in COMMON
C
C  more comments on arguments:
C  "indicial spaces" range from 0 to 1 with equal spacing for all N points
C  in the corresponding dimension
C
C smoothing widths
C  zdelt & zdelx -- positive:  standard smooth, negative:  dbl inv. smooth
C
C |f-f~| limits
C  zepst & zepsx -- positive:  absolute limit,  negative:  fractional limit
C
      real, dimension(:), allocatable :: ztime
      real zxarr(nr0)
C
C----------------------------------------------------
C
      lunt=lunzer(0)
      ier=0
C
      allocate(ztime(ntime)); ztime=0.0
C
C  compare parameters in SMWORK_BLK with CPLOTR parameters
C  NB dmc SMWORK_BLK arrays moved to datmgr_mod -- DMC Nov. 2009
C
      imaxp=max(NR0,NTIME)
      if(imaxp.gt.NSM) then
         ier=ier+1
         write(lunt,9001) nsm,imaxp
 9001    format(
     >      ' ?smfunc -- error in NSM value: NSM=',i6,' IMAXP=',i6)
      endif
C
C  check eps type parameters
C
      if(index('%RA',typ_epst).eq.0) then
         ier=ier+1
         call zermsg(' ?smfunc:  unexpected eps(t) type:  '//typ_epst)
      endif
C
      if(index('%RA',typ_epsx).eq.0) then
         ier=ier+1
         call zermsg(' ?smfunc:  unexpected eps(x) type:  '//typ_epsx)
      endif
C
      if(istat.eq.0) then
         ier=ier+1
         call zermsg(' ?smfunc:  smoother input data not ready.')
      endif
C
      if(ier.gt.0) return
C
C-----------------------------
C
      iwarn=0
C
      iflag=0
      if(istat.lt.0) then
         inumt=ntt
         inx=1
         zxarr(1)=0.0
         ipx=0
      else
         inumt=ntr
         inx=nzonex(istat)
         zxarr(1)=0.0
         ipx=0
         if(inx.gt.1) then
            if((indx.eq.0).and.(nlxfot(istat))) then
               CALL DMGXOT(ISTAT,IND1,IND2)
               IPX=LOCD(IND1)           ! ptr for time varying x axis...
               IF(ISTAT.EQ.2) IPX=LOCD(IND2)
            else
               ipx=0                    ! time invariant x axis...
               do ix=1,inx
                  if(indx.eq.1) then
                     zxarr(ix)=float(ix-1)/(inx-1)
                  else
                     zxarr(ix)=xarry(ix,istat)
                  endif
               enddo
            endif                       ! x axis fixed or varying
         endif                          ! no. of x pts .gt. 1
      endif                             ! scalar or profile
C
      do it=1,inumt
         if(indt.eq.1) then
            ztime(it)=float(it-1)/float(ntt-1)
         else
            if(istat.lt.0) then
               ztime(it)=time(it)
            else
               ztime(it)=time3(it)
            endif
            if(it.gt.1) then
               if(ztime(it).le.ztime(it-1)) iflag=iflag+1
            endif
         endif
      enddo
C
      if(iflag.gt.0) then
         call zermsg(' ?smfunc:  time axis not monotonic increasing.')
         call zermsg('  for smoothing delta(t) requires "I" option.')
         ier=99
         return
      endif
C
C----------------------
C
C  1.  smooth vs. x
C
      if((inx.gt.1).and.(zdelx.ne.0.0)) then
C
         do it=1,ntr
            if((indx.eq.0).and.(nlxfot(istat))) then
               ixl=ipx+(it-1)*inx
               call copyr4(datbuf(ixl),zxarr,inx)
            endif
C  monotonicity check
            iflag=0
            do ix=2,inx
               if(zxarr(ix).le.zxarr(ix-1)) iflag=iflag+1
            enddo
            if(iflag.gt.0) then
               iwarn=iwarn+1
               if(iwarn.le.3) then
                  write(lunt,9005) ztime(it)
 9005             format(
     >               ' %smfunc:  x axis not monotonic at t=',1pe11.4/
     >               '           smoothing vs. x skipped at that time.')
               else if(iwarn.eq.4) then
                  call zermsg(' %smfunc: max warning count reached.')
               endif
            else
C  OK smooth
               ixf=ipf+(it-1)*inx - 1
               do ix=1,inx
                  smwork(ix,1)=datbuf(ixf+ix)
               enddo
C
               call smoof1(zdelx,zepsx,typ_epsx,zxarr,inx)
C
               do ix=1,inx
                  datbuf(ixf+ix)=smwork(ix,2)
               enddo
            endif
         enddo
C
      endif
C
C  2.  smooth vs. t
C
      if(zdelt.ne.0.0) then
         do ix=1,inx
C
            do it=1,inumt
               smwork(it,1)=datbuf(ipf+(it-1)*inx+(ix-1))
            enddo
C
            call smoof1(zdelt,zepst,typ_epst,ztime,inumt)
C
            do it=1,inumt
               datbuf(ipf+(it-1)*inx+(ix-1))=smwork(it,2)
            enddo
C
         enddo
      endif
C
      return
      end
