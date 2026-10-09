      subroutine rptimav(zmgbuf,ibufsize,nfcns,istype,ttarg,delta_t,
     >   outbuf,ioutsize,ier)
C
C  interpolate or time average data -- removing time dependence and
C  reducing data to a single time point.
C
      use datmgr_mod
      use cplotr_mod
C
      integer ibufsize,nfcns,istype,ioutsize
      real zmgbuf(ibufsize,nfcns)
      real outbuf(ioutsize,nfcns)
      real ttarg,delta_t
C
C  where:  zmgbuf = the data to be interpolated or averaged
C                   from a prior RPMGCALC or RPMG*CAL call, or an
C                   individual function or calculator result (set nfcns=1).
C
C          ibufsize = 1st array dimension of zmgbuf, (max) size per item.
C
C          nfcns = number of functions to time average
C
C          istype = function subtype code (will determine actual item size)
C
C          ttarg = time to interpolate to or average around
C          delta_t = average +/- delta_t around ttarg
C              if delta_t is NEGATIVE then 1/[time average of (1/data)]
C              is computed.  ("double inverse averaging").
C
C output:  outbuf = array into which to write results of interpolation
C          ioutsize = (max) size per result
C
C          ier = completion code, 0 = normal
C
C
C               RPTIMAV declares
C                   real zmgbuf(ibufsize,nfcns)
C                   real outbuf(ioutsize,nfcns)
C
C  local:
C
      integer idims(10)                 ! 2 would be enough dmc 9/1999
C
      real zbuf1(nr0),zbuf2(nr0)
C
      character*10 zxabb(8)
C
C---------------
C
      lunt=lunzer(0)
      call rpdims(istype,irank,idims,zxabb,ier)
      if(ier.ne.0) return
C
      itot=1
      do i=1,irank
         itotp=itot                     ! itotp = size per time point
         itot=itot*idims(i)
      enddo
C
      if(itot.gt.ibufsize) then
         write(lunt,9001) istype,itot,ibufsize
 9001    format(' ?RPTIMAV:  type ',i3,' functions need ',i6,' words,'/
     >          '            but only ',i6,' were provided.')
         ier=1
         return
      endif
C
      if(itotp.gt.ioutsize) then
         write(lunt,9002) istype,itotp,ioutsize
 9002    format(' ?RPTIMAV:  type ',i3,' functions need ',i5,' words,'/
     >          '       per time point but only ',i5,' were provided.')
         ier=1
         return
      endif
C
      if(irank.eq.1) then
         ztmin=time(1)
         ztmax=time(ntt)
      else
         ztmin=time3(1)
         ztmax=time3(ntr)
      endif
C
      if((ttarg.lt.ztmin).or.(ttarg.gt.ztmax)) then
         write(lunt,9003) ttarg,ztmin,ztmax
 9003    format(' ?RPTIMAV:  target time = ',1pe12.5,' seconds'/
     >      '   not in range of data:  ',1pe12.5,' to ',1pe12.5,
     >      ' seconds.')
         ier=2
         return
      endif
C
      if(irank.eq.1) then
C  scalar function -- find time bin
         if(ntt.eq.1) then
            it1=1
            it2=1
            zfact=0.0
         else
            do it=1,ntt-1
               if((time(it).le.ttarg).and.(ttarg.le.time(it+1))) then
                  it1=it
                  it2=it+1
                  if(time(it).eq.time(it+1)) then
                     zfact=0.5
                  else
                     zfact=(ttarg-time(it))/(time(it+1)-time(it))
                  endif
                  go to 50
               endif
            enddo
         endif
      else
C  profile function -- find time bin
         if(ntr.eq.1) then
            it1=1
            it2=1
            zfact=0.0
         else
            do it=1,ntr-1
               if((time3(it).le.ttarg).and.(ttarg.le.time3(it+1))) then
                  it1=it
                  it2=it+1
                  if(time3(it).eq.time3(it+1)) then
                     zfact=0.5
                  else
                     zfact=(ttarg-time3(it))/(time3(it+1)-time3(it))
                  endif
                  go to 50
               endif
            enddo
         endif
      endif
C
 50   continue                          ! loop breakout
C
      if(delta_t.eq.0.0) then
         iopt=0
      else if(delta_t.gt.0.0) then
         iopt=1
         zdelta=delta_t
      else if(delta_t.lt.0.0) then
         iopt=2
         zdelta=abs(delta_t)
      endif
C
      do if=1,nfcns
         if(irank.eq.1) then
            itotp=1
            if(iopt.eq.0) then
               zbuf1(1)=zmgbuf(it1,if)
               zbuf2(1)=zmgbuf(it2,if)
            else
               call smtima0(time,zmgbuf(1,if),ntt,zbuf1,1,
     >            itotp,it1,zdelta,iopt)
               call smtima0(time,zmgbuf(1,if),ntt,zbuf2,1,
     >            itotp,it2,zdelta,iopt)
            endif
         else
            if(iopt.eq.0) then
               do ix=1,itotp
                  iadr=(it1-1)*itotp + ix
                  zbuf1(ix)=zmgbuf(iadr,if)
                  iadr=(it2-1)*itotp + ix
                  zbuf2(ix)=zmgbuf(iadr,if)
               enddo
            else
               call smtima0(time3,zmgbuf(1,if),ntr,zbuf1,1,
     >            itotp,it1,zdelta,iopt)
               call smtima0(time3,zmgbuf(1,if),ntr,zbuf2,1,
     >            itotp,it2,zdelta,iopt)
            endif
         endif
         do ix=1,itotp
            outbuf(ix,if)=(1.0-zfact)*zbuf1(ix)+zfact*zbuf2(ix)
         enddo
      enddo
C
      return
      end
