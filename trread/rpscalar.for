      subroutine rpscalar(zname,zbuf,ibufsize,iret,ier)
C
C  fetch a scalar function of time -- data only, use rptime_s for timebase.
C     use rplabel for labels.
C
      use datmgr_mod
      use cplotr_mod

      character*(*) zname
      real zbuf(ibufsize)               ! buffer into which to copy data
C
      integer iret                      ! number of data points copied
      integer ier                       ! completion code, 0 = OK
C
C  ier=1 means:  invalid name
C  ier=2 means:  insufficient buffer space
C  ier=3 means:  rplot internal error (seek help)
C
C  iret=0 if an error occurs.
C-------------------------------------------------------------
C
      character*10 zabbr
C
C-------------------------------------------------------------
C
      iret=0
      if(ibufsize.lt.ntt) then
         call rpbufsiz('rpscalar',ibufsize,ntt)
         ier=2
         return
      endif
C
      zabbr=zname
      call trcaps(zabbr)
      iadr=ifind_ordr(abt,iordrt,nft,zabbr)
C
      if(iadr.eq.0) then
         call zermsg(
     >      ' ?rpscalar:  not a scalar function name:  '//zname)
         ier=1
         return
      endif
C
      call dmgfotx(2,ipt,ier)
      if(ier.ne.0) then
         ier=3
         return
      endif
C
      ier=0
C
      ipf=ipt+(iadr-1)*ntt
C
      do it=1,ntt
         zbuf(it)=datbuf(ipf+it-1)
      enddo
      iret=ntt
C
      return
      end
