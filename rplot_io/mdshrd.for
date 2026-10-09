      subroutine mdshrd(lun,ier,ixflag)
C
c  lun -- i/o unit for messages
c  ier -- completion code (returned)
c  ixflag -- 1st run or "secondary" run flag/index
c
C Get Function and Multigraph names, units and labels
C equivalent of reading TF.PLN file or cdfhrd
C
C 01/06/00 CAL
c
c  dmc 23 Mar 2000 -- adding TF.PLN based caching, for performance
c
C
C-----------------------------------------------------------------------
      use cplotr_mod
      implicit none
 
      integer       lun, ier, ixflag

C Functions 
      integer       Mds_Value
      integer       idescr_cstring, idescr_long,  idescr_longarr
      integer       idescr_cstringarr
      integer       ifind_ordr

      character*120 expr
 
      integer       status, size,  year
      integer       idims(1), dsc
 
      character*10  zstr
      character*12  zdev, zdev2, nodes(15)
 
      integer       i, j, k, lzdev, ilt, ily, iltf
      integer       nzall, nz, nzl, nzu, nzx
      character*7   tokyr
      character*24  mds_tree
      integer       ift, ifxt
      integer       str_length
      integer       ildir, ilc, iok, jj, il, IXR, ictf, idum
C
      character*10, allocatable :: zx(:), cont(:)
      character*64, allocatable :: zlabs(:)
      character*10, allocatable :: zabs(:)
      character*32, allocatable :: zunits(:)
      integer,      allocatable :: iordt(:), iordr(:), iordx(:)
C
      character*1 ztmp_char
C
      character*20 cache_name_list(naxfxt)
      logical      ztf_contents
      integer      ier_cache, iwr_cache, ier_wr 
C
C------------------------------------------------------------------
C
      allocate(zlabs(NAXFXT+NAXFOT))
      allocate(zabs(NAXFXT+NAXFOT))
      allocate(zunits(NAXFXT+NAXFOT))
      allocate(iordt(NAXFOT))
      allocate(iordr(NAXFXT))
      allocate(iordx(NAXXVR))
      allocate(zx(NAXFXT), cont(NAXMGF))
C
      ztf_contents=.false.
C make runlabel
      if (ixflag.eq.0) then
         ily=index(fdisk,'@')+1
         mds_tree=fdisk(ily:)
         ily=lfdisk-ily+1
         tokyr=fdir
         ilt=lfdir
         if (tokyr .eq. ' ') then
            ilt=0
            if(mds_cache) then
               call mds_cache_fname(ixflag,'DEVICE_NAME.DAT',tfiln,
     >              ier_cache)
               if(ier_cache.ne.0) then
                  if(mds_cache_only) then
                     ier=1
                     write(lun,*) ' ?mds_cache_fname(DEVICE) error '//
     >                    '(irrecoverable: mds_cache_only=T).'
                     goto 777
                  else
                     mds_cache=.FALSE.
                     write(lun,*) ' %mds_cache_fname(DEVICE) error,'//
     >                    ' cache disabled.'
                  endif
               else
                  call mds_cache_getdev(tfiln,zdev,ilt)
                  if(mds_cache_only) then
                     if(ilt.eq.0) then
                        write(lun,*) ' ?mdshrd: mds_cache_only=T '//
     >                       'and mds_cache_getdev failed.'
                        ier=1
                        goto 777
                     endif
                  endif
               endif
            endif
            if(ilt.eq.0) then
               status = Mds_Value('machine()',idescr_cstring(zdev2),
     >              size)
               if (mod(status,2).ne.1) then
                  call mdserr(lun, 'machine()', status)
                  zdev2='????'
               endif
               status = Mds_Value(':DEVICE',idescr_cstring(zdev),size)
               if (mod(status,2).ne.1) then
                  call mdserr(lun, ':DEVICE', status)
                  zdev=zdev2
               else
                  if(zdev.ne.zdev2) then
                     write(lun,*)
     >                    ' %mdshrd: mds_value("machine()",...) = "',
     >                    zdev2,'"'
                     write(lun,*)
     >                    ' %mdshrd: mds_value(":device",...) = "',
     >                    zdev,'"'
                     if (zdev2(1:1) .eq. ' ') then
                        write(lun,*) ' --> using ":device".'
                     else
                        write(lun,*) ' --> using "machine()" value.'
                        zdev=zdev2
                     endif
                  endif
               endif
               if(mds_cache) then
                  call mds_cache_putdev(tfiln,zdev)
               endif
            endif
            lzdev = index(zdev,' ')-1
            RUNLB2=zdev(1:lzdev)//'.'//mds_tree(1:ily)//
     >          ' '//runid(1:LRUNID)//' (MDS+)'
         else
            RUNLB2=tokyr(1:ilt)//' '// RUNID(1:lrunid)//' (MDS+)'
         endif
      endif

c--------------------------------------------------
      nzall = (NAXFXT+NAXFOT)
      ier = 0
      if(ixflag.eq.0) then
         nlxvar=.TRUE.
      endif
c
c  cache file
c
      ictf=0
      tmp_filn=tfiln
      ztmp_char=filnc
      iwr_cache=0

      if(mds_cache) then
c
         filnc='Y'
c
c  binary objects cache directory file -- & cache validation
c
         call mds_cache_fname(ixflag,'ITEMS.LIST',tfiln,ier_cache)
         if(ier_cache.eq.0) then
c
c  get list of cached names -- ilc = #, can be zero
c
            ildir=str_length(tfiln)-10  ! path to "ITEMS.LIST": tfiln(1:ildir)
            call mds_cache_list(ixflag,tfiln,ildir,cache_name_list,ilc,
     >         naxfxt,ier_cache)
c
c  mark cached time files, if available
c
            if(ier_cache.eq.0) then
               if(ixflag.eq.0) then
                  call mds_cache_list_precon(cache_name_list,ilc,
     >                 nbcache_act,nbcache_list,ier_cache)
               else
                  call mds_cache_list_precon(cache_name_list,ilc,
     >                 nbcache_act_x(ixflag),nbcache_list_x(1,ixflag),
     >                 ier_cache)
               endif
            endif
         endif
c
         if(ier_cache.ne.0) then
            if(mds_cache_only) then
               write(lun,*) 
     >              ' RPLOT_CACHE_ONLY = TRUE -> cannot recover.'
               ier=1
               goto 777
            endif
            write(lun,*)
     >           ' %mdshrd: cache disabled, software error!'
            mds_cache=.false.
         endif
c
c  label cache file
c
         call mds_cache_fname(ixflag,'TF.PLN',tfiln,ier_cache)
         ictf=1
         iwr_cache=0
         iok=0
         if(ier_cache.eq.0) then
            write(lun,*) ' %attempting MDS+ label cache read'
            iwr_cache=-1                ! suppress error messages
            call tfilrd(lun_tf,iwr_cache,lrun_x)
            if(iwr_cache.ne.0) then
c
c  will try subsequent cache write -- but cannot for ixflag.ne.0
c
               if(mds_cache_only) then 
                  write(lun,*) 
     >                 ' RPLOT_CACHE_ONLY = TRUE -> cannot recover.'
                  ier=1
                  goto 777
               endif

               write(lun,*) ' %cache miss, reverting to MDS+'
               if(ixflag.ne.0) iwr_cache=0
c
            else
c
c  cache read successful
c  OK have all labels now:  read timebase
c
               write(lun,*) ' %OK-- reading labels from cache'
               call mdstrd(ixflag,ier)
               iok=1                    ! signal that labels are read already
c
            endif
         endif
c
c  exit if labels read from cache
c
         if(iok.eq.1) go to 2000  ! have labels from cache...
C
      endif
C
      if(mds_cache_only) then
         ier=1
         write(lun,*) ' RPLOT_CACHE_ONLY = TRUE -> cannot recover.'
         go to 2000
      endif
C
C-----------------------------------------------------------------
C Try to extract TF.PLN from MDS node .CONTENTS
C
      if(ictf.eq.0) then
         ! if ictf=1 use cache filename for CONTENTS
         call tmpfile('TFPLN_',tfiln,iltf)
      else
         iltf=len(trim(tfiln))
      endif

      call mds_getfile(lun,'CONTENTS',tfiln(1:iltf),ier)

      if (ier .eq. -1) then
         write(lun,*)
     >        ' %mdshrd: node CONTENTS does not exist'
         if(ictf.eq.1) then
            ictf=0
            call tmpfile('TFPLN_',tfiln,iltf)
         endif
      else if (ier .ne. 0) then
         write(lun,*)
     >        ' %mdshrd: error reading CONTENTS'
         if(ictf.eq.1) then
            ictf=0
            call fdelete(tfiln,idum)
            call tmpfile('TFPLN_',tfiln,iltf)
         endif
      else
         write(lun,*)
     >        ' %mdshrd: TF.PLN written from CONTENTS'
         if(ictf.eq.0) ztf_contents=.true.  ! set to delete tmp file
C     read TF.PLN extracted from CONTENTS node
         if(ictf.eq.1) then
            write(lun,*) ' %mdshrd:  read TF.PLN'
         else
            write(lun,*) ' %mdshrd:  read ',tfiln(1:iltf)
         endif
         ier=0
C    dont over-write runid_x in CPLOTR common: "-lrun_x" 
         call tfilrd(lun_tf,ier,-lrun_x)
         if (ier .eq. 0) then
            call mdstrd(ixflag,ier)
            if (ier .eq. 0) then
               write(lun,*) ' %mdshrd: got labels from CONTENTS'
               iwr_cache=1     ! write TF cache file
               go to 400       ! have labels from TF.PLN
            else
               write(lun,*) ' %mdshrd:  mdstrd returned error'
            endif
         else
            write(lun,*) ' %mdshrd:  tfilrd returned error'
         endif
      endif
c
c
c----------------------------------
C
C Read from MDSplus
C
      if(ixflag.eq.0) then
C Get shot and run
         status = Mds_Value(':SOURCE_SHOT',idescr_long(NSHOT),size)
         if (mod(status,2).ne.1) then
            call mdserr(lun, ':SOURCE_SHOT', status)
            nshot=99999
         endif
 
         ily=index(fdisk,'@')+1
         mds_tree=fdisk(ily:)
         ily=lfdisk-ily+1
         tokyr=fdir
         ilt=lfdir
      else
         ily=index(fdisk_x(ixflag),'@')+1
         mds_tree=fdisk_x(ixflag)(ily:)
         ily=index(mds_tree,' ')-1
         tokyr=fdir_x(ixflag)
         ilt=index(tokyr,' ')-1
         if (ilt.ge.0) ilt=len(tokyr)
      endif
C
C Get Time arrays
      call mdstrd(ixflag,ier)
      if(ier.ne.0) go to 2000
C
C One-dimensional Functions
C--------------------------
      idims(1) = NAXFOT
      expr = 'GETNCI(".OUTPUTS.ONE_D:*","NODE_NAME")'
 
      if (ixflag.eq.0) then
         dsc = idescr_cstringarr(ABT(1),idims,1)
         status = Mds_Value(expr,dsc,NFT)
cxx      print *,'MDSHRD: NFT =',NFT
      else
         dsc = idescr_cstringarr(ABT_X(1,ixflag),idims,1)
         status = Mds_Value(expr,dsc,NFT_X(ixflag))
cxx      print *,'MDSHRD: NFT =',NFT,' NFT_X =',NFT_X(ixflag)
      endif
      if (mod(status,2).ne.1) then
         call mdserr(lun, expr, status)
      endif
 
 
 
C Two-dimensional Functions
      idims(1) = NAXFXT
      expr = 'GETNCI(".OUTPUTS.TWO_D:*","NODE_NAME")'
      if (ixflag.eq.0) then
         dsc = idescr_cstringarr(ABR(1),idims,1)
         status = Mds_Value(expr,dsc,NFXT)
cxx      print *,'MDSHRD: NFXT =',NFXT
      else
         dsc = idescr_cstringarr(ABR_X(1,ixflag),idims,1)
         status = Mds_Value(expr,dsc,NFXT_X(ixflag))
cxx      print *,'MDSHRD: NFXT =',NFXT,' NFXT_X', NFXT_X(ixflag)
      endif
 
 
 
C Get all ABBREVIATIONS
      idims(1)=nzall
      expr = 'GETNCI(".TRANSP_OUT:*","NODE_NAME")'
      dsc = idescr_cstringarr(zabs(1),idims,1)
      status = Mds_Value(expr,dsc,nz)
      if (mod(status,2).ne.1) then
         call mdserr(lun, expr, status)
      endif
cxx      print *,'MDSHRD: all ABBR =',nz
 
C Get all LABELS
      expr = 'GETNCI(".TRANSP_OUT:*:RPLABEL","RECORD")'
      dsc = idescr_cstringarr(zlabs(1),idims,1)
      status = Mds_Value(expr,dsc,nzl)
      if (mod(status,2).ne.1) then
         call mdserr(lun, expr, status)
      endif
cxx      print *,'MDSHRD: all LABELS =',nzl
C
      if (nz .ne. nzl .or.
     >    nz .ne. (nft+nfxt)) then
         call mdshrd_slow(lun,ier,ixflag)
         if (ixflag .eq. 0) then
            call aordr(iordr,ABR,NFXT)
            call aordr(iordt,ABT,NFT)
         else
            call aordr(iordr,ABR_X(1,ixflag),NFXT_X(ixflag))
            call aordr(iordt,ABT_X(1,ixflag),NFT_X(ixflag))
         endif
         goto 1000
      endif
 
C Get all UNITS
      expr = 'GETNCI(".TRANSP_OUT:*:UNITS","RECORD")'
      dsc = idescr_cstringarr(zunits(1),idims,1)
      status = Mds_Value(expr,dsc,nzu)
      if (mod(status,2).ne.1) then
         call mdserr(lun, expr, status)
      endif
cxx      print *,'MDSHRD: all UNITS =',nzu
 
      if (nzu .ne. nzl) then
         call mdshrd_slow(lun,ier,ixflag)
         if (ixflag .eq. 0) then
            call aordr(iordr,ABR,NFXT)
            call aordr(iordt,ABT,NFT)
         else
            call aordr(iordr,ABR_X(1,ixflag),NFXT_X(ixflag))
            call aordr(iordt,ABT_X(1,ixflag),NFT_X(ixflag))
         endif
         goto 1000
      endif
C----------
cxx      print *,' ... rearranging labels, etc'
      if (ixflag .eq. 0) then
         call aordr(iordr,ABR,NFXT)
         call aordr(iordt,ABT,NFT)
         do j=1,nz
            i=ifind_ordr(ABT,iordt,NFT,zabs(j))
            if (i.gt.0) then
               LABELT(i) = zlabs(j)
               UNITST(i) = zunits(j)
               go to 500
            endif
            i=ifind_ordr(ABR,iordr,NFXT,zabs(j))
            if (i.gt.0) then
               LABELR(i) = zlabs(j)
               UNITSR(i) = zunits(j)
               go to 500
            endif
            write(lun,*) 'no match for ',zabs(j)
 500        continue
         end do
      else
         call aordr(iordr,ABR_X(1,ixflag),NFXT_X(ixflag))
         call aordr(iordt,ABT_X(1,ixflag),NFT_X(ixflag))
         do j=1,nz
            i=ifind_ordr(ABT_X(1,ixflag),iordt,NFT_X(ixflag),
     >        zabs(j))
            if (i.gt.0) then
               LABELT_X(i,ixflag) = zlabs(j)
               UNITST_X(i,ixflag) = zunits(j)
               go to 502
            endif
            i=ifind_ordr(ABR_X(1,ixflag),iordr,NFXT_X(ixflag),
     >        zabs(j))
            if (i.gt.0) then
               LABELR_X(i,ixflag) = zlabs(j)
               UNITSR_X(i,ixflag) = zunits(j)
               go to 502
            endif
 502        continue
         end do
      endif
 
C      if (ixflag .eq. 0) then
C         do i=1,20
C            write (6,*) ABT(i),' : ',UNITST(i),':',LABELT(i)
C         enddo
C      else
C         do i=1,20
C            write (6,*) ABT_X(i,ixflag),' : ',UNITST_X(i,ixflag),
C     >            ':',LABELT_X(i,ixflag)
C         enddo
C      endif
C
C Get all XAXIS
      idims(1)=NAXFXT
      expr = 'GETNCI(".TRANSP_OUT:*:XAXIS","RECORD")'
      dsc = idescr_cstringarr(zx(1),idims,1)
      status = Mds_Value(expr,dsc,nzx)
      if (mod(status,2).ne.1) then
         call mdserr(lun, expr, status)
      endif
cxx      print *,'MDSHRD: all XAXIS =',nzx
 
C
C Profile Functions of Time and X
C---------------------------------
 
      IXR=2
      if (ixflag .eq. 0) then
         IFXT=NFXT
      else
         IFXT=NFXT_X(ixflag)
      endif
cxx      print *,'MDSHRD: NFXT =',IFXT
cxx      print *,'MDSHRD: fill XAXIS'
      do i=1,IFXT
C check if X-axis
         if (ixflag .eq. 0) then
            zstr=ABR(i)(1:10)
         else
            zstr = ABR_X(i,ixflag)
         endif
         if (zstr .eq. zx(i)) then
            if (zstr .eq. 'X') then
               jj=1
            else if  (zstr .eq. 'XB') then
               jj=2
            else
               IXR=IXR+1
               jj=IXR
            endif
            if (ixflag .eq. 0) then
               XNDABB(jj) = zx(i)
               XLAB(jj)=LABELR(i)
               NFX(jj)=i
               dsc = idescr_long(NZONEX(jj))
            else
               XNDABB_X(jj,ixflag) = zx(i)
               NFX_X(jj,ixflag)=i
               dsc = idescr_long(NZONEX_X(jj,ixflag))
            endif
            il = index(zstr  ,' ')-1
            if (il .le. 0) il=10
            expr='size(.TRANSP_OUT:'// zstr(1:il) //',0)'
            status = Mds_Value(expr,dsc,size)
            if (mod(status,2).ne.1) then
               call mdserr(lun, expr, status)
            endif
         endif
      end do
 
cxx      if (ixflag .eq. 0) then
cxx         print *,'MDSHRD: NXR =',IXR
cxx         do i=1, IXR
cxx            print *, XNDABB(i),NFX(i),' ',XLAB(i), NZONEX(i)
cxx         end do
cxx      else
cxx         print *,'MDSHRD: NXR_X(',ixflag,') =',IXR
cxx         do i=1, IXR
cxx            print *, XNDABB(i),NFX_X(i,ixflag),' ',XLAB(i),
cxx     >             NZONEX_X(i,ixflag)
cxx         end do
cxx      endif
 
C Find index into XNDABB
      if(ixflag.eq.0) then
         call aordr(iordx,XNDABB,IXR)
      else
         call aordr(iordx,XNDABB_X(1,ixflag),ixr)
      endif
      do i=1,IFXT
         if(ixflag.eq.0) then
            j=ifind_ordr(XNDABB,iordx,IXR,zx(i))
         else
            j=ifind_ordr(XNDABB_X(1,ixflag),iordx,IXR,zx(i))
         endif
         if (j.gt.0) then
            if (ixflag .eq. 0) then
               ITYPR(i)=j
            else
               ITYPR_X(i,ixflag)=j
            endif
            go to 200
         endif
         write(lun,*) zx(i),' could not find XNDABB'
 200     continue
      end do
 
      if (ixflag .eq. 0) then
         NXR = IXR
      else
         NXR_X(ixflag)=IXR
      endif
 
C---------------------------------------------------------------
 1000 continue
C if secondary run, don't read Multigraphs
      if (ixflag .gt. 0) go to 400
C
C Check for multigraphs
C======================
cxx      print *,'MDSHRD: Reading Multigraphs'
      idims(1) = 15
      dsc = idescr_cstringarr(nodes(1),idims,1)
      expr = 'GETNCI(".*","NODE_NAME")'
      status = Mds_Value(expr,dsc,size)
      if (mod(status,2).ne.1) then
         call mdserr(lun, expr, status)
      else
         do i=1,size
            if (nodes(i)(1:11) .eq. 'MULTIGRAPHS') then
               go to 300
            endif
         end do
         NBAL = 0
         write(lun,*) 'MDSHRD: no Multigraphs'
         go to 400
      endif
 
 300  continue
C Multigraphs
C-----------
      idims(1) = NAXMGP
      dsc = idescr_cstringarr(ABB(1),idims,1)
      expr = 'GETNCI(".MULTIGRAPHS.*","NODE_NAME")'
      status = Mds_Value(expr,dsc,NBAL)
      if (mod(status,2).ne.1) then
         call mdserr(lun, expr, status)
      endif
C
cxx      print *,'MDSHRD: NBAL =',NBAL
cxx      print *,'MDSHRD: Reading Labels of Multigraphs'
C Label
      dsc = idescr_cstringarr(LABELB(1),idims,1)
      expr = 'GETNCI(".MULTIGRAPHS.*:LABEL","RECORD")'
      status = Mds_Value(expr,dsc,size)
      if (mod(status,2).ne.1) then
         call mdserr(lun, expr, status)
      endif
      if (size .ne. NBAL) then
         write(lun,*) '?MDSHRD: number of Labels =',size,
     >      ' ,NBAL=',NBAL
      endif
 
      do i=1,NBAL
         il = index(ABB(i),' ')-1
         if (il .le. 0) il=10
C CONSIGN
         idims(1)=naxmgf
         dsc = idescr_longarr(ifunb(1,i),idims,1)
         expr='.MULTIGRAPHS.'//ABB(i)(1:il)//':CONSIGN'
         status = Mds_Value(expr,dsc,INFB(i))
         if (mod(status,2).ne.1) then
            call mdserr(lun, expr, status)
         endif
 
C Contents
         dsc = idescr_cstringarr(cont(1),idims,1)
         expr='.MULTIGRAPHS.'//ABB(i)(1:il)//':CONTENT'
         status = Mds_Value(expr,dsc,INFB(i))
         if (mod(status,2).ne.1) then
            call mdserr(lun, expr, status)
         endif
C Find index into 1D or 2D Functions
         do 100 j=1,INFB(i)
            k=ifind_ordr(ABT,iordt,NFT,cont(j))
            if (k.gt.0) then
               go to 110
            endif
            k=ifind_ordr(ABR,iordr,NFXT,cont(j))
            if (k.gt.0) then
               go to 120
            endif
            write(lun,*) cont(j),' does not match anything'
            go to 100
C 1D
 110        IFUNB(j,i) = ifunb(j,i) * k
            IINTB(i) = 1
            UNITSB(i) = UNITST(k)
            go to 100
C 2D Profile
 120        IFUNB(j,i) = ifunb(j,i) * k
            IINTB(i) = 0
            UNITSB(i) = UNITSR(k)
            go to 100
 100     end do
      end do
cxx     do i=1,20
cxx        write (6,*) ABB(i),' : ',UNITSB(i),':',LABELB(i),INFB(i),
cxx    >    IINTB(i)
cxx     enddo
C
cxx      print *,'MDSHRD: done'
C End of MDS reads
C------------------------------------------------------------------
C
 400  continue
C Check if need to delete temporary TF.PLN file
      if (ztf_contents) call fdelete(tfiln(1:iltf),ier)
C Check if need to write TF Cache file
      if((ier.eq.0).and.(ixflag.eq.0)) then
         if(iwr_cache.ne.0) then
            call mds_cache_fname(ixflag,'TF.PLN',tfiln,ier_cache)
            if(ier_cache.eq.0) then
               call tfilwr(lun_tf,ier_wr)
               if(ier_wr.ne.0) then
                  write(lun,*) ' ?mdshrd:  cache write error!'
               endif
            endif
         endif
      endif
C
 2000 continue
c
      tfiln=tmp_filn
      filnc=ztmp_char
c
      if(mds_cache) then
         if(ixflag.eq.0) then
            call mds_cache_list_convert(ixflag,cache_name_list,ilc,
     >         nbcache_act,nbcache_list,ier_cache)
         else
            call mds_cache_list_convert(ixflag,cache_name_list,ilc,
     >         nbcache_act_x(ixflag),nbcache_list_x(1,ixflag),
     >         ier_cache)
         endif
      endif
c
 777  deallocate(zlabs,zabs,zunits,iordt,iordr,iordx,zx,cont)
c
      return
      end
C---------------------------------------------------------------------------
      subroutine mdshrd_slow(lun,ier,ixflag)
C
C If some subnodes are missing e.g. old CMOD trees,
C need to read individual nodes
      use cplotr_mod
      implicit none
 
      integer       lun, ier, ixflag
 
      integer       ift, ifxt
      character*120 expr
      integer       status, size
      character*10  zstr, zab, zx(NAXFXT)
      integer       idims(1), dsc, zdsc
      integer       iordx(NAXXVR)
 
      integer       i, ii, il, j, jj, ixr
      integer       ifind_ordr
C Functions
      integer       Mds_Value      
      integer       idescr_cstring, idescr_long
C-----------------------------------------------------------

      
      print *,'MDSHRD: Tree contains empty nodes!'
      print *,'        will read Labels of individual nodes'

C One-dimensional Functions
C--------------------------
 
      if (ixflag.eq.0) then
         ift=NFT
         ifxt=NFXT
      else
         ift=NFT_X(ixflag)
         ifxt=NFXT_X(ixflag)
      endif
cxx      print *,'MDSHRD: NFT =',IFT
cxx      print *,'MDSHRD: Reading Labels of f(t) individually'
      i=0
      do ii=1,IFT
         i=i+1
         if (ixflag.eq.0) then
            dsc = idescr_cstring(LABELT(i))
            zab=ABT(i)(1:10)
         else
            dsc = idescr_cstring(LABELT_X(i,ixflag))
            zab=ABT_X(i,ixflag)
         endif
         il=index(zab,' ')-1
         if (il .le. 0) il=10
         expr='.TRANSP_OUT:'//zab(1:il)//':RPLABEL'
         status = Mds_Value(expr,dsc,size)
         if (mod(status,2).ne.1) then
            if (status .ne. 265388144)
     >           call mdserr(lun, expr, status)
C -- remove node if there are no UNITS
C -- for CMOD, the node is in MODEL Tree, but not in run
            IFT=IFT-1
            if (ixflag.eq.0) then
               do j=i,IFT
                  ABT(j)=ABT(j+1)
               end do
            else
               do j=i,IFT
                  ABT_X(j,ixflag)=ABT_X(j+1,ixflag)
               end do
            endif
            i=i-1
         else
            if (ixflag.eq.0) then
               dsc = idescr_cstring(UNITST(i))
            else
               dsc = idescr_cstring(UNITST_X(i,ixflag))
            endif
            expr='.TRANSP_OUT:'//zab(1:il)//':UNITS'
            status = Mds_Value(expr,dsc,size)
            if (mod(status,2).ne.1) then
               if (status .ne. 265388144)
     >            call mdserr(lun, expr, status)
C -- remove node if there are no UNITS
C -- for CMOD, the node is in MODEL Tree, but not in run
               IFT=IFT-1
               if (ixflag.eq.0) then
                  do j=i,IFT
                     ABT(j)=ABT(j+1)
                  end do
               else
                  do j=i,IFT
                     ABT_X(j,ixflag)=ABT_X(j+1,ixflag)
                  end do
               endif
               i=i-1
            endif
         endif
      end do
cxx      print *,'MDSHRD: NFT =',IFT
 
      if (ixflag.eq.0) then
         NFT=ift
cxx         do i=1,20
cxx            write (6,*) ABT(i),' : ',UNITST(i),':',LABELT(i)
cxx         enddo
      else
         NFT_X(ixflag)=ift
      endif
 
C
C Profile Functions of Time and X
C---------------------------------
cxx      print *,
cxx     > 'MDSHRD: Reading Labels of Time and addl. Profile Functions'
      idims(1) = NAXFXT
 
      i=0
      IXR=2
cxx      print *,'MDSHRD: NFXT =',NFXT
      do ii=1,IFXT
         i=i+1
         if (ixflag.eq.0) then
            dsc = idescr_cstring(LABELR(i))
            zab=ABR(i)(1:10)
         else
            dsc = idescr_cstring(LABELR_X(i,ixflag))
            zab=ABR_X(i,ixflag)
         endif
         il=index(zab,' ')-1
         if (il .le. 0) il=10
         expr='.TRANSP_OUT:'//zab(1:il)//':RPLABEL'
         status = Mds_Value(expr,dsc,size)
         if (mod(status,2).ne.1) then
            if (status .ne. 265388144)
     >           call mdserr(lun, expr, status)
C     -- remove node if there are no RPLABEL
C     -- for CMOD, the node is in MODEL Tree, but not in run
            IFXT=IFXT-1
            if (ixflag.eq.0) then
               do j=i,IFXT
                  ABR(j)=ABR(j+1)
               end do
            else
               do j=i,IFXT
                  ABR_X(j,ixflag)=ABR_X(j+1,ixflag)
               end do
            endif
            i=i-1
         else
            if(ixflag.eq.0) then
               dsc = idescr_cstring(UNITSR(i))
            else
               dsc = idescr_cstring(UNITSR_X(i,ixflag))
            endif
            expr='.TRANSP_OUT:'//zab(1:il)//':UNITS'
            status = Mds_Value(expr,dsc,size)
            if (mod(status,2).ne.1) then
               if (status .ne. 265388144)
     >          call mdserr(lun, expr, status)
C -- remove node if there are no UNITS
C -- for CMOD, the node is in MODEL Tree, but not in run
               IFXT=IFXT-1
               if (ixflag.eq.0) then
                  do j=i,IFXT
                     ABR(j)=ABR(j+1)
                  end do
               else
                  do j=i,IFXT
                     ABR_X(j,ixflag)=ABR_X(j+1,ixflag)
                  end do
               endif
               i=i-1
            else
               zdsc = idescr_cstring(zstr)
               expr='.TRANSP_OUT:'//zab(1:il)//':XAXIS'
               status = Mds_Value(expr,zdsc,size)
               if (mod(status,2).ne.1) then
                  call mdserr(lun, expr, status)
               else
                  zx(i)=zstr
C     check if X-axis
                  if (zab .eq. zx(i)) then
                     if (zab .eq. 'X') then
                        jj=1
                     else if  (zab .eq. 'XB') then
                        jj=2
                     else
                        IXR=IXR+1
                        jj=IXR
                     endif
                     if (ixflag .eq. 0) then
                        XNDABB(jj) = zx(i)
                        XLAB(jj)=LABELR(I)
                        NFX(jj)=i
                        ITYPR(i)=jj
                        dsc = idescr_long(NZONEX(jj))
                     else
                        xndabb_x(jj,ixflag) = zx(i)
                        NFX_X(jj,ixflag)=i
                        ITYPR_X(i,ixflag)=jj
                        dsc = idescr_long(NZONEX_X(jj,ixflag))
                     endif
                     expr='size(.TRANSP_OUT:'// zab(1:il) //',0)'
                     status = Mds_Value(expr,dsc,size)
                     if (mod(status,2).ne.1) then
                        call mdserr(lun, expr, status)
                     endif
                  endif
               endif
            endif
         endif
      end do
      if (ixflag.eq.0) then
         NFXT=ifxt
         NXR=ixr
      else
         NFXT_X(ixflag)=ifxt
         NXR_X(ixflag)=ixr
      endif
 
cxx      print *,'MDSHRD: NXR =',IXR
cxx      do i=1, IXR
cxx         print *, XNDABB(i),NFX(i),' ',XLAB(i), NZONEX(i)
cxx      end do
 
C Find index into XNDABB
      if(ixflag.eq.0) then
         call aordr(iordx,XNDABB,IXR)
      else
         call aordr(iordx,XNDABB_X(1,ixflag),ixr)
      endif
      do i=1,IFXT
         if(ixflag.eq.0) then
            j=ifind_ordr(XNDABB,iordx,IXR,zx(i))
         else
            j=ifind_ordr(XNDABB_X(1,ixflag),iordx,IXR,zx(i))
         endif
         if (j.gt.0) then
            if (ixflag .eq. 0) then
               ITYPR(i)=j
            else
               ITYPR_X(i,ixflag)=j
            endif
            go to 200
         endif
         write(lun,*) zx(i),' could not find XNDABB'
 200     continue
      end do
 
cxx      print *,' '
cxx      print *,'MDSHRD: NFXT =',IFXT
 
C      do i=1,20
C         write (6,*) ABR(i),' : ',UNITSR(i),':',LABELR(i),
C     >    ': ',XNDABB(itypr(i))
C      enddo
 
      return
      end
