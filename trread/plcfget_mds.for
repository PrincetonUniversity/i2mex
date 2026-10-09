      subroutine plcfget_mds(zpathi,zdisk,zdir,ier)
c
c  create MDS+ syntax wanted by trprofil & trscalar, from
c  %RFETCH MDS+ syntax:
c
c    zpath = MDS+:<server>:<treename>:<tok.yy> -->
c
c     zdisk = MDS+-<server>@<treename>
c     zdir = <tok.yy>
c
c    zpath = MDS+:<server>:<treename> -->
c
c     zdisk = MDS++<server>@<treename>
c     zdir = ' '
c
c  also support variant syntax forms:
c
c    MDS+:<server>:<treename>(<shot#>)  -->
c
c     zdisk = MDS++<server>@<treename>
c     zdir = ' '
c
c    MDS+:<server>:<treename>(<tok.yy>,<runid>)  -->
c
c     zdisk = MDS+-<server>@<treename>
c     zdir = <tok.yy>
c
c
      character*(*) zpathi              ! MDS+ path info (input)
c
      character*(*) zdisk               ! trprofil/trscalar disk arg (output)
      character*(*) zdir                ! trprofil/trscalar dir arg (output)
c
      integer ier                       ! completion code: 0=OK
c
c------------------------------
c  it is already known that the string starts with "MDS+:".
c
      character*140 zpath
      integer icolon(4),ilen,idigs,itest,ios
      character*6 ztestnn
c
c------------------------------
c
      iparen=index(zpathi,'(')
      if(iparen.eq.0) then
         zpath=zpathi
      else
         icomma=index(zpathi,',')
         if(icomma.gt.0) then
            zpath='MDS+-'//zpathi(6:iparen-1)//':'//
     >         zpathi(iparen+1:icomma-1)
         else
            zpath='MDS++'//zpathi(6:iparen-1)
         endif
      endif
c
      ilen=len_trim(zpath)
      if(zpath(ilen:ilen).eq.':') ilen=ilen-1 ! ignore trailing ":"
c
      icolon(1)=5
      icc=1
c
      do i=6,ilen
         if(zpath(i:i).eq.':') then
            icc=icc+1
            if(icc.gt.4) go to 900
            icolon(icc)=i
         endif
      enddo
c
      if(icc.lt.2) go to 900
c
      if(icc.gt.2) then
c
c  might have a port number...
c
         idigs=icolon(3)-icolon(2)-1
         if((idigs.gt.0).and.(idigs.le.6)) then
            ztestnn='000000'
            ztestnn(6-idigs+1:6)=zpath(icolon(2)+1:icolon(3)-1)
            read(ztestnn,'(i6)',iostat=ios) itest
            if(ios.eq.0) then
c
c  yes, a port number (pure integer between colons); absorb this as
c  a part of the server name
c
               do i=3,icc
                  icolon(i-1)=icolon(i)
               enddo
               icc=icc-1
c
            endif   ! ios
         endif   ! idigs
      endif   ! icc
c
      if(icc.eq.2) then
         zdir=' '
         icolon(3)=ilen+1
         zdisk(1:5)='MDS++'
      else
         zdir=zpath(icolon(3)+1:ilen)
         call uupper(zdir)
         zdisk(1:5)='MDS+-'
      endif
c
      zdisk(6:)=zpath(icolon(1)+1:icolon(2)-1)//'@'//
     >   zpath(icolon(2)+1:icolon(3)-1)
      call uupper(zdisk)
c
      ier=0
      return
c
c---------------------
c  parse error
c
 900  continue
      lunt=lunzer(0)
      write(lunt,*) ' ?plcfget (RFETCH) -- invalid MDS+ path string:'
      write(lunt,*) '  ',zpath
c
      ier=1
      return
      end
