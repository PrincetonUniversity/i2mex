      subroutine rprofile(zname,zbuf,ibufsize,iret,ier)
C
C  fetch a profile function of time + addl coordinate--
C  data only, use rptime_p for timebase.  use rplabel for labels,
C  rpdims for dimensioning information.
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
      character*10 zabbr,zxabb(8)
C
      character*64 zlbl
      character*32 zuns
C
      integer idims(8)
C
C-------------------------------------------------------------
C
      iret=0
C
      zabbr=zname
      call trcaps(zabbr)
      iadr=ifind_ordr(abr,iordrr,nfxt,zabbr)
C
      if(iadr.eq.0) then
         call zermsg(
     >      ' ?rprofile:  not a profile function name:  '//zname)
         ier=1
         return
      endif
C
      call rplabel(zabbr,zlbl,zuns,imulti,istype)
      call rpdims(istype,irank,idims,zxabb,ier)
      if(ier.ne.0) then
         ier=3
         return
      endif
C
      isize=1
      do ir=1,irank
         isize=isize*idims(ir)
      enddo
      if(ibufsize.lt.isize) then
         ier=2
         call rpbufsiz('rprofile',ibufsize,isize)
         return
      endif
C
      CALL DMGFXT(IADR,IND)
      IPF=LOCD(IND)
C
      ier=0
      do i=1,isize
         zbuf(i)=datbuf(ipf+i-1)
      enddo
      iret=isize
C
      return
      end
