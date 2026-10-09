      subroutine cdfhrd(file_id, retn, ixflag)
C Read netCDF Header
C
C modification 14 Nov 1997 -- dmc
C dmc:  ixflag argument added, for support of multiple runs open
C              at the same time;
C              the CPLOTR COMMON has been expanded, with capacity to
C              hold label information for multiple runs.
C              ixflag = 0 -- store info in the usual arrays
C              ixflag.gt.0 -- store info in "extra run" arrays
C
C dmc 14 Nov 1997 -- bugfix, ntt.eq.ntr was assumed, incorrectly.
C
C 09/22/97  CAL
 

      use datmgr_mod
      use cplotr_mod
 
      implicit NONE
      include "netcdf.inc"

C Input
      integer  file_id
      integer  ixflag
C Return
      integer  retn
C
C Processing
      integer  ndims, t_id, t3_id
      integer  d_ids(2)
      integer  nvars, ngatts, unlimid

      integer :: luntrm,lunzer
      integer :: inum,intmax,inta,ixr,izonex,itypex,infbc
      integer :: ishot,ifxt,ift,ibal,izones
 
      integer*2   dummy(naxmgf)  ! To convert int*2 into int*4
C
      integer*1 ibyte        ! To convert byte into int*4
C
      integer  i, j, k
      real     atreal(2)
C
      character*10 zrunid
      character*64 zlabel
      character*32 zunits
      character*10 zabr
      character*10 zndabb(naxxvr)
C
      character*40 zstrng
C--------------------------------------------------------------------------
C
      luntrm=lunzer(0)
      write(luntrm,*) ' ... reading NetCDF header data ...'
C
C Inquire about file
      retn = nf_inq(file_id, ndims, nvars, ngatts, unlimid)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9002) NF_STRERROR(retn)
 9002    format(' ? CDFHRD: nf_inq',A)
         return
      end if
cdbg      write(6,*) nvars,' variables,',ndims,' dimensions', ngatts,' att'
 
      if(ixflag.eq.0) then
         nlxvar=.TRUE.
      endif
 
C Get Global Attributes
C======================
C      do i=1,ngatts
C         retn = nf_inq_attname(file_id,NF_GLOBAL,i,atname)
C         write(6,*) atname
C         retn = nf_inq_att(file_id,NF_GLOBAL,atname,xtyp,j)
C         write(6,*) xtyp,j
C      end do
 
 
C  shot number
      retn = nf_get_att_int(file_id,NF_GLOBAL,'shot',ishot)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9001) NF_STRERROR(retn)
 9001    format(' ? CDFHRD: nf_get_att_int '/,'     ',A)
         return
      else if(ixflag.eq.0) then
         nshot=ishot
      end if
 
C  no. of fcns f(x,t)
      retn = nf_get_att_int(file_id,NF_GLOBAL,'NFXT',ifxt)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9001) NF_STRERROR(retn)
         return
      else if(ixflag.eq.0) then
         nfxt=ifxt
      else
         nfxt_x(ixflag)=ifxt
      end if
 
C  no. of fcns f(t)
      retn = nf_get_att_int(file_id,NF_GLOBAL,'NFT',ift)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9001) NF_STRERROR(retn)
         return
      else if(ixflag.eq.0) then
         nft=ift
      else
         nft_x(ixflag)=ift
      end if
 
C  no. of mg packages
      retn = nf_get_att_int(file_id,NF_GLOBAL,'NBAL',ibal)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9001) NF_STRERROR(retn)
         return
      else if(ixflag.eq.0) then
         nbal=ibal
      end if
 
C  std. no. of zones
      retn = nf_get_att_int(file_id,NF_GLOBAL,'NZONES',izones)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9001) NF_STRERROR(retn)
         return
      else if(ixflag.eq.0) then
         nzones=izones
      else
         nzones_x(ixflag)=izones
      end if
 
      retn = nf_get_att_text(file_id,NF_GLOBAL,'Runid',zrunid)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9003) NF_STRERROR(retn)
 9003    format(' % CDFHRD: nf_get_att_text -- '/,
     >        '     ',A)
CCC         return         ! warning only; keep going
      end if
C
      retn = nf_get_att_real(file_id,NF_GLOBAL,'R',atreal)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9013) NF_STRERROR(retn)
 9013    format(' ? CDFHRD: nf_get_att_real '/,'     ',A)
         return
      else if(ixflag.eq.0) then
         rmajor = atreal(1)
         rminor = atreal(2)
      end if
 
cdbg      write(6,*) 'shot',nshot,' ',runid,rmajor,rminor
 
C Read Profile Functions
C=======================
C
C Get Dimensions
C----------------
C Time
      intmax = 0
      t_id = 1
      t3_id = 1
      retn = nf_inq_dimlen(file_id,t_id,inum)
      if (retn .ne. NF_NOERR) then
         write(luntrm,9007) NF_STRERROR(retn)
 9007    format(' ? CDFHRD: nf_inq_dimlen -- Time'/,
     >        '     ',a)
         return
      end if
C
      intmax=max(intmax,inum)
C
C  assume only one time axis for now...
      inta=1
      if(ixflag.eq.0) then
         ntt = inum
         ntr = ntt
      else
         ntt_x(ixflag) = inum
         ntr_x(ixflag) = inum
      endif
C
      call get_vatt(file_id, 1, 'TIME', zlabel, zunits,
     >     luntrm)
      if(ixflag.eq.0) then
         timlab=zlabel
         timuns=zunits
      else
         timlab_x(ixflag)=zlabel
         timuns_x(ixflag)=zunits
      endif
C
C  check for 2nd time axis
C
      retn = nf_inq_dim(file_id,2,zabr,inum)
      if(zabr.eq.'TIME3') then
         inta=2
         t3_id=2
         if(ixflag.eq.0) then
            ntr=inum
         else
            ntr_x(ixflag)=inum
         endif
         intmax=max(intmax,inum)
      endif
C
C  check "xdatmgr" time vector sizes
C
      call dmg_texpand(intmax)
C
C  store in COMMON:  whether 1 or 2 time dimensions were read:
C
      if(ixflag.eq.0) then
         ncdft=inta
      else
         ncdft_x(ixflag)=inta
      endif
C
C X axes
      ixr = ndims - inta
      if(ixflag.eq.0) then
         nxr=ixr
      else
         nxr_x(ixflag)=ixr
      endif
 
      do i=1,ixr
         j = i+inta
         retn = nf_inq_dim(file_id,j,zabr,izonex)
         if (retn .ne. NF_NOERR) then
            write(luntrm,9006) j, NF_STRERROR(retn)
 9006       format(' ? CDFHRD: nf_inq_dim',i2/,
     >           '     ',a)
            return
         else if(ixflag.eq.0) then
            xndabb(i)=zabr
            nzonex(i)=izonex
         else
            nzonex_x(i,ixflag)=izonex
            xndabb_x(i,ixflag)=zabr
         end if
         zndabb(i)=zabr                 ! locally used
      end do
 
C
C
C Get Functions
C----------------
cdbg      write(6,*) ' '
cdbg      write(6,*) 'Profile Functions',nfxt
      j = inta
      do i = 1,ifxt
         j = j+1
         retn = nf_inq_varname(file_id, j, zabr)
         if (retn .ne. NF_NOERR) then
            write(luntrm,9004) j, NF_STRERROR(retn)
 9004       format(' ? CDFHRD: nf_inq_varname ',i3/,
     >           '     ',a)
            return
         end if
 
         retn = nf_inq_vardimid(file_id, j, d_ids)
         itypex = d_ids(1)-inta
         if (retn .ne. NF_NOERR) then
            write(luntrm,9005) j, NF_STRERROR(retn)
 9005       format(' ? CDFHRD: nf_inq_vardimid ',i3/,
     >           '     ',a)
            return
         else if(ixflag.eq.0) then
            itypr(i) = itypex           !  First dimension is Time
         else
            itypr_x(i,ixflag) = itypex    !  First dimension is Time
         end if
 
C Get Attributes
         call get_vatt(file_id, j, zabr, zlabel,  zunits,
     >                 luntrm)
         if(ixflag.eq.0) then
            abr(i)=zabr
            labelr(i)=zlabel
            unitsr(i)=zunits
         else
            abr_x(i,ixflag)=zabr
            labelr_x(i,ixflag)=zlabel
            unitsr_x(i,ixflag)=zunits
         end if
 
C Is Variable a Dimension ?
         if (zabr .eq. zndabb(itypex)) then
            if(ixflag.eq.0) then
               nfx(itypex) = i
               xlab(itypex) = labelr(i)
            else
               nfx_x(itypex,ixflag) = i
            endif
         end if
 
C         write(6,101) labelr(i),unitsr(i),abr(i),itypex,i
C 101     format(3a,2i5)
 
      end do
 
cdbg      write(6,*) ' '
cdbg      write(6,*) 'X Dimensions'
 
C      do i = 1, nxr
C         write(6,102) nroffx(i), nzonex(i),xlab(i), xndabb(i)
C 102      format(2i5,1x,2a)
C       end do
 
C
C Read Scalar Functions
C======================
cdbg      write(6,*) ' '
cdbg      write(6,*) 'Scalar Functions',nft
 
      j = ifxt+inta
      do i=1,ift
         j=j+1
C     Get name
         retn = nf_inq_varname(file_id, j, zabr)
         if (retn .ne. NF_NOERR) then
            write(luntrm,9011) j, NF_STRERROR(retn)
 9011       format(' ? CDFHRD - nf_inq_varname ',i3/,
     >           '     ',a)
            return
         end if
         call get_vatt(file_id, j, zabr, zlabel, zunits, luntrm)
         if(ixflag.eq.0) then
            abt(i)=zabr
            labelt(i)=zlabel
            unitst(i)=zunits
         else
            abt_x(i,ixflag)=zabr
            labelt_x(i,ixflag)=zlabel
            unitst_x(i,ixflag)=zunits
         end if
C         write(6,103)labelt(i),unitst(i),abt(i)
 103     format(1x,3a)
      end do
C
C Read Multi Graphs
C==================
cdbg      write(6,*) ' '
cdbg      write(6,*) 'Multi Graphs',nbal
 
      do i = 1,ibal
         j = j+1
C     Get name
         retn = nf_inq_varname(file_id, j, zabr)
         if (retn .ne. NF_NOERR) then
            write(luntrm,9004) j, NF_STRERROR(retn)
            return
         end if
C Get Attributes
         call get_vatt(file_id, j, zabr, zlabel, zunits, luntrm)
         if(ixflag.eq.0) then
            abb(i)=zabr
            labelb(i)=zlabel
            unitsb(i)=zunits
         endif
         retn = nf_inq_attlen(file_id, j, 'Fct_Ids', infbc)
         if (retn .ne. NF_NOERR) then
            write(luntrm,9008) abb(i), NF_STRERROR(retn)
 9008       format(' ? CDFHRD: nf_inq_attlen - ',a/
     >           '     ',a)
            return
         else if(ixflag.eq.0) then
            infb(i)=infbc
         end if
C
         retn = nf_get_att_int2(file_id, j,'Fct_Ids', dummy)
         if (retn .ne. NF_NOERR) then
            write(luntrm,9009) abb(i), NF_STRERROR(retn)
 9009       format(' ? CDFHRD: nf_get_att - Fct_Ids - ',a/
     >           '     ',a)
            return
         end if
         if(ixflag.eq.0) then
            do k=1,infbc
               ifunb(k,i) = dummy(k)
            end do
         endif
 
C        Read into byte before converting into integer
         retn = nf_get_var_int1(file_id, j, ibyte)
         if (retn .ne. NF_NOERR) then
            write(luntrm,9010) abb(i), NF_STRERROR(retn)
 9010       format(' ? CDFHRD: nf_get_var_int1 - ',a/
     >           '     ',a)
            return
         end if
         if(ixflag.eq.0) then
            iintb(i) = ibyte
         endif 
C         write(6,104)labelb(i),unitsb(i),iintb(i),infb(i),abb(i)
C         write(6,105)(ifunb(k,i),k=1,infb(i))
C 104     format(1x,2a,2i5,1x,a)
C 105     format(1x,20i4)
      end do
 
C Read Time Dimension
C====================
      if(allocated(time)) then
         write(6,*) ' cdfhrd: size(time) = ',size(time)
      else
         write(6,*) ' ?? cdfhrd: time UNALLOCATED!'
      endif
      if(allocated(time3)) then
         write(6,*) ' cdfhrd: size(time3) = ',size(time3)
      else
         write(6,*) ' ?? cdfhrd: time3 UNALLOCATED!'
      endif
      if(ixflag.eq.0) then
         retn = nf_get_var_real(file_id, t_id, time)
         if( retn .eq. NF_NOERR ) then
            if(t_id.ne.t3_id) then
               retn = nf_get_var_real(file_id, t3_id, time3)
            else
               call copyr4(time,time3,ntt)
            endif
         endif
      else
         retn = nf_get_var_real(file_id, t_id, time_x(1,ixflag))
         if( retn .eq. NF_NOERR ) then
            if(t_id.ne.t3_id) then
               retn = nf_get_var_real(file_id, t3_id, time3_x(1,ixflag))
            else
               call copyr4(time_x(1,ixflag),time3_x(1,ixflag),
     >            ntt_x(ixflag))
            endif
         endif
      endif
      if (retn .ne. NF_NOERR) then
         write(luntrm,9014) NF_STRERROR(retn)
 9014    format(' % CDFHRD: nf_get_vara_real - Time'/,
     >        '     ',a)
         ntr = 0
         return
      end if
 
      return
      end
C--------------------------------------------------------------------
      subroutine get_vatt(fid,vid,name,label,units,luntrm)
C
C
      include "netcdf.inc"
 
C Input
      integer fid     ! File Id
      integer vid     ! Variable Id
      character*(*)  name  ! for error message
      integer luntrm  ! LUN for error message
C Returns
      character*(*)  label
      character*(*)  units
 
      integer retn,retn2
C
C  dmc  --  try to handle NetCDF flakiness problem
C    sometimes NetCDF can only store the 1st four letters of
C    the attribute name, it seems.  Or maybe a DEC UNIX bug?
C
      character*20 zstrng,zstrng0
C----------------------------------------------------------------------
C
C     Get Units
C
      call get_vattck(fid,vid,name,units,'units',1,luntrm,ier)
      if(ier.ne.0) return
C
      call get_vattck(fid,vid,name,label,'long_name',2,luntrm,ier)
      if(ier.ne.0) return
C
      return
      end
C
C--------------------------------------------------------------------
      subroutine get_vattck(fid,vid,name,label,lblname,lblid,luntrm,ier)
C
C  get a single attribute with error check / recovery
C
      include "netcdf.inc"
C
C Input
      integer fid     ! File Id
      integer vid     ! Variable Id
      character*(*)  name  ! for error message
      character*(*) lblname  ! label desired:  "units" or "long_name"...
      integer lblid   ! backup path:  attribute # for desired label
      integer luntrm  ! LUN for error message
C Returns
      character*(*)  label              ! label string value returned
      integer ier                       ! completion code, 0=normal
C
      integer ilenl,ilena
C
      integer retn,retn2
C
C  dmc  --  try to handle NetCDF flakiness problem
C    sometimes NetCDF can only store the 1st four letters of
C    the attribute name, it seems.  Or maybe a DEC UNIX bug?
C
      character*20 zstrng
C----------------------------------------------------------------------
C
      itry=0
      ilenl=len(label)
      zstrng=lblname
 10   continue
C
      ier=0
      itry=itry+1
      retn = nf_inq_attlen(fid,vid,zstrng,ilena)
      if(retn.eq.NF_NOERR) then
         if(ilena.lt.ilenl) then
            label(ilena+1:ilenl)=' '
         endif
         ilena=min(ilena,ilenl)
         retn = nf_get_att_text(fid,vid,zstrng,label(1:ilena))
      endif
      if (retn .ne. NF_NOERR) then
         if(itry.eq.1) then
            retn2=nf_inq_attname(fid,vid,lblid,zstrng)
            if(retn2 .eq. NF_NOERR) then
               if(zstrng(1:4).eq.lblname(1:4)) go to 10  ! retry
            else
               zstrng='att-error'
            endif
         endif
         write(luntrm,9005) name, NF_STRERROR(retn)
 9005    format(' ? CDFHRD: nf_get_att_text ',a/,
     >      '     ',a)
         ier=1
      end if
C
      if((itry.gt.1).and.(ier.eq.0)) then
         ilna=index(name,' ')-1
         if(ilna.le.0) ilna=len(name)
         ilnn=index(lblname,' ')-1
         if(ilnn.le.0) ilnn=len(lblname)
         ilnb=index(zstrng,' ')-1
         if(ilnb.le.0) ilnb=len(zstrng)
         write(luntrm,9901) name(1:ilna),lblname(1:ilnn),zstrng(1:ilnb)
 9901    format(
     >      ' %cdfhrd warning:  ',a,' "',a,'" attribute name was:  "',
     >      a,'".')
      endif
C
      return
      end
