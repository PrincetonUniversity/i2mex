      subroutine rconnect(zrunstr,ier)
C
      use cplotr_mod
C
C  connect to a TRANSP run
C
      character*(*) zrunstr
C
C    zrunstr -- runid:  <path>/<runid> or just <runid>
C
C  *** modified *** dmc 29 Feb 2000:  MDS+ support.  the above syntax
C  can be used, or the above prepended with FILE:, as in
C
C    FILE:<unix-file-spec>
C
C    MDS+:<server-name>:<tree-name>(<shot-number>) or
C    MDS+:<server-name>:<tree-name>(<tok.yy>,<runid>)
C
      integer ier                       ! error code returned; 0 = OK
C
C  assume "initcpl" has been called.
C
      character*140 zbuf
      character*5 ztest
      character*6 zext(4)
      character*6  zhold
C
      character*60 mds_server
      character*30 mds_tree
      character*30 mds_id1
      character*30 mds_id2
      character*20 mds_tstr
      character*1 zchar
C
C
C--------------------------
      zext(1) = ".CDF  "
      zext(2) = "TF.PLN"
      zext(3) = "MF.PLN"
      zext(4) = "NF.PLN"
C
C     call cplset_exec          ! link the block data...
C     call rpcaldat_exec        ! link the block data...
C
C  verify RPLOT system initialization
C
      call initcpl
C
      fdisk=' '
      lfdisk=0
      fdir=' '
      lfdir=0
C
      nlmds=.FALSE.
C
      ztest=zrunstr(1:5)
      call uupper(ztest)
      if(ztest(1:5).eq.'FILE:') then
         zbuf=zrunstr(6:)
      else if(ztest(1:5).eq.'MDS+:') then
         zbuf=zrunstr(6:)
         call trmds_parse(zbuf,mds_server,mds_tree,mds_id1,mds_id2,
     >      mds_tstr,ierr)
         if(ierr.ne.0) go to 999
c
         if(mds_id2.eq.' ') then
            zchar='+'
         else
            zchar='-'
         endif
c
         ils=len_trim(mds_server)
         fdisk='MDS+'//zchar//mds_server(1:ils)//'@'//mds_tree
         lfdisk=len_trim(fdisk)
c
         if(mds_id2.eq.' ') then
            runid=mds_id1
         else
            runid=mds_id2
            fdir=mds_id1
            lfdir=len_trim(fdir)
         endif
c
         nlmds=.TRUE.
c
      else
         zbuf=zrunstr
      endif
C
      if(.not.nlmds) then
         ilen=len_trim(zbuf)
c
         ! remove extensions
         do ix = 1, size(zext)
            ix_len = len_trim(zext(ix))  ! length of test string
            i = max(1,ilen+1-ix_len)     ! start of extension in buffer
            zhold = zbuf(i:ilen)         ! copy to holding string
            call uupper(zhold)           ! make upper case
            if (zhold == zext(ix)) then
               zbuf(i:) = ' '            ! spaces filled to end
               ilen = i-1
            end if
         end do
c
         do i=ilen,1,-1
            if(zbuf(i:i).eq.'/') go to 15
         enddo
         ifi=1
         go to 20
C
 15      continue
         ifi=i+1
         call ufilnam(zbuf(1:i),' ',fdir)
         lfdir=len_trim(fdir)
C
 20      continue
         runid=zbuf(ifi:)
         lrunid=len_trim(runid)
      endif
C
cx      write(6,*) ' lfdir,fdir =   ',lfdir,' ',fdir(1:max(1,lfdir))
cx      write(6,*) ' lfdisk,fdisk = ',lfdisk,' ',fdisk(1:max(1,lfdisk))
cx      write(6,*) ' runid =        ',runid
C
      ier=0
      runlb2=' '
 25   call inirun(ier)
      if(ier.eq.2) then
         call arc_chk(runid,iwait)
         if(iwait.eq.1) then
            NLMDS = .TRUE.
            ier=-777
            go to 25
         endif
      else if(ier.ne.0) then
         return
      endif
C
C  OK...
C  setup calculator
C
      call rpctini
C
      return
C
C  MDS+ syntax error
C
 999  continue
C
      ier=1
      il=len_trim(zrunstr)
      write(lunzer(0),9991) zrunstr(1:il)
 9991 format(' ?rconnect:  MDS+ file spec error, was expecting:'/
     >   '   "MDS+:<server-name>:<tree-name>(<shot-number>)" or'/
     >   '   "MDS+:<server-name>:<tree-name>(<tok.yy>,<runid>)"; got:'/
     >   6x,a)
C
      return
      end
