subroutine trgs2fetch(zpath,zxid,zprefix,ier)
  use datmgr_mod
  use cplotr_mod
  use rpcalc_mod

  implicit NONE
 
  !  read tabulated GS2 data in format established by R. Budny Sep 2001
  !  (example:  see source/trread/gs2.dat).
  !
  !  an MDS+ access option is planned but has not been coded yet.
 
  character*(*) zpath                   ! path to GS2 tabulated data
  character*(*) zxid                    ! x axis id: 1st letter "X" or "R"
  character*(*) zprefix                 ! prefix for output profile names
  integer ier
 
  external XIDENT
 
  !---------------------------------------------------
  !  local...
 
  integer ilz,ilp,ist,ic,ibrk,ilpre
  integer ind2,isiz2,iwrk2,ifcnx,icur,idum
  integer lunt,lunzer
  integer itype,inx,istat,iopen,iclass
 
  character*1 zxtest
  character*10 zxab,zfab
  character*150 zfile,zline
  character*20 zbuf
 
  integer ipass,inumx,imaxx,inumt,istat2,istat3,ii,inz
 
  real, dimension(:), allocatable :: gs2time
  real, dimension(:), allocatable :: gs2x1,gs2r1,gs2ak1,gs2omeg1,gs2gamm1
  integer nt_gs2
 
  integer, dimension(:), allocatable :: nx_gs2
  real, dimension(:,:), allocatable :: gs2xa,gs2aky,gs2omega,gs2gamma
 
  integer ind,ipt,it,ifcn,ifind_ordr
  integer ix,ia0,ia1
 
  real zti,ztic
  integer it0,it0p1,iextrap
 
  integer nflds
  character*30 cfields(6)
  real zfields(6)
 
  !---------------------------------------------------
 
  ier=0
  !
  ! WORKSPACE
 
  CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)
 
  !---------------------------------------------------
  !  ...basic checks...
 
  lunt=lunzer(0)
 
  zxtest=zxid(1:1)
  call uupper(zxtest)
  if((zxtest.ne.'X').and.(zxtest.ne.'R')) then
     write(lunt,*) ' ?trGS2fetch: XID argument not "R" or "X" as expected.'
     ier=ier+1
  else if(zxtest.eq.'X') then
     itype=1
     inx=nzonex(itype)
     zxab='X'
  else
     itype=ntypmr
     if(itype.le.0) then
        write(lunt,*) ' ?trGS2fetch:  XID="R" but no "R" X axis available.'
        ier=ier+1
     else
        inx=nzonex(itype)
     endif
     zxab='RMAJM'
  endif
 
  ilpre=len_trim(zprefix)
  if(ilpre.gt.5) then
     write(lunt,*) ' ?trGS2fetch: PREFIX argument max length = 5 characters.'
     write(lunt,*) '  value given was "',zprefix(1:ilz),'".'
     ier=ier+1
  else
     call idchek0(zprefix,iclass,1)
     if(iclass.ne.0) then
        write(lunt,*) ' ?trGS2fetch: error with PREFIX argument.'
        ier=ier+1
     endif
  endif
 
  !---------------------------------------------------
  ! ....... get filename, open file
 
  ilp=len_trim(zpath)
  ibrk=0
  do ic=ilp,1,-1
     if(zpath(ic:ic).eq.'/') then
        ibrk=ic
        exit
     endif
  enddo
  if(ibrk.eq.0) then
     if(zpath(1:1).eq.'$') then
        call get_environment_variable(zpath(2:ilp),zfile)
     else
        zfile=zpath
     endif
  else
     ist=1
     if(zpath(1:1).eq.'$') ist=2
     call ufilnam(zpath(ist:ibrk),zpath(ibrk+1:ilp),zfile)
  endif
 
  iopen=0
  open(unit=lun_nf,file=zfile,status='OLD',iostat=istat)
  if(istat.ne.0) then
     ilz=len_trim(zfile)
     write(lunt,*) ' ?trGS2fetch:  open failure: ',zfile(1:ilz)
     if(zpath.ne.zfile) then
        write(lunt,*) '   PATH argument was: ',zpath(1:ilp)
     endif
     ier=ier+1
  else
     iopen=1
  endif
 
  !---------------------------------------------------
 
  !  write(6,*) ' ier = ',ier
  !  write(6,*) ' zpath = "',zpath(1:max(1,len_trim(zpath))),'"'
  !  write(6,*) ' zxid = "',zxid(1:max(1,len_trim(zxid))),'"'
  !  write(6,*) ' zprefix = "',zprefix(1:max(1,len_trim(zprefix))),'"'
 
  if(ier.ne.0) go to 1000
 
  !---------------------------------------------------
  !  ........scan file for times & max no. of spatial pts at any time pt
  !
  !  first assume 100 timepoints; if that doesn't work allow for 1000 time pts
  !  then 10000 etc...
 
  ipass=0
  imaxx=0
  nt_gs2=10
  do
     ipass=ipass+1
     if(ipass.gt.1) rewind lun_nf
 
     nt_gs2=nt_gs2*10
     allocate(gs2time(nt_gs2),stat=istat)
     if(istat.ne.0) then
        write(lunt,*) &
             ' ?trGS2fetch:  gs2time(nt_gs2) allocate failure, nt_gs2=',nt_gs2
        ier=1
        go to 1000
     endif
 
     inumt=0
     inumx=0
     istat=0
     istat2=0
     istat3=0
     imaxx=0
     do
        read(lun_nf,'(A)',iostat=istat) zline
        if(istat.ne.0) then
           imaxx=max(imaxx,inumx)
           exit
        endif
 
        ii=index(zline,'#')
        if(ii.gt.0) zline(ii:len(zline))=' '  ! ignore comment fields
 
        if(index(zline,'.').eq.0) cycle   ! ignore lines without decimal pts
 
        ii=index(zline,'time=')
        if(ii.gt.0) then
 
           ! read next time point
 
           imaxx=max(imaxx,inumx)
           inumx=0
           inumt=inumt+1
           if(inumt.gt.nt_gs2) exit
           ilz=len_trim(zline)
           inz=ilz-(ii+5)+1
           zbuf=' '
           zbuf(20-inz+1:20)=zline(ii+5:ilz)
           read(zbuf,'(G20.0)',iostat=istat2) gs2time(inumt)
           if(istat2.ne.0) exit
           if(inumt.gt.1) then
              if(gs2time(inumt).le.gs2time(inumt-1)) then
                 istat3=1
                 exit
              endif
           endif
        else
           inumx=inumx+1  ! count line with decimal data, no "time="
        endif
     enddo
 
     if(inumt.gt.nt_gs2) then
        deallocate(gs2time)
        cycle
     endif
 
     if(istat.ne.0) exit
     if(istat2.ne.0) exit
     if(istat3.ne.0) exit
 
  enddo
 
  if(istat2.gt.0) then
     write(lunt,*) ' ?trGS2fetch: floating point decode error in line: '
     write(lunt,*) ' ',zline(1:ilz)
     ier=90
  endif
 
  if(inumt.eq.0) then
     write(lunt,*) ' ?trGS2fetch:  no time pts in file'
     ier=91
  else if(istat3.gt.0) then
     write(lunt,*) ' ?trGS2fetch: time not in ascending order: ',gs2time(inumt)
     ier=92
  endif
 
  if(imaxx.eq.0) then
     write(lunt,*) ' ?trGS2fetch: max no. of x pts is zero in file.'
     ier=93
  endif
 
  if(ier.gt.0) go to 1000
 
  !---------------------------------
  !  OK... can read data now
 
  rewind lun_nf
  nt_gs2=inumt
 
  allocate(nx_gs2(nt_gs2))
  allocate(gs2x1(imaxx),gs2r1(imaxx), &
       gs2ak1(imaxx),gs2omeg1(imaxx),gs2gamm1(imaxx))
  allocate(gs2xa(inx,nt_gs2))
  allocate(gs2aky(inx,nt_gs2),gs2omega(inx,nt_gs2),gs2gamma(inx,nt_gs2))
 
  !---------------------------------
  !  form standard time evolving x axis grid for data
  !  interpolate from run database's time evolving grid to GS2 times
 
  ifcnx=ifind_ordr(abr,iordrr,nfxt,zxab)
  call dmgfxt(ifcnx,ind)
  ipt=locd(ind)
 
  do it=1,nt_gs2
 
     call xinter(XIDENT,gs2time(it),TIME3,NTR, &
          it0,it0p1,zti,ztic,iextrap)
 
     do ix=1,inx
 
        ia0=ipt+(it0-1)*inx+ix-1
        ia1=ipt+(it0p1-1)*inx+ix-1
 
        gs2xa(ix,it)=ztic*datbuf(ia0)+zti*datbuf(ia1)
 
     enddo
  enddo
 
  !-----------------------------------
  ! now scan the file for data
 
  it=0
 
  do
     read(lun_nf,'(A)',iostat=istat) zline
     if(istat.ne.0) then
        if(it.gt.0) then
           if(nx_gs2(it).gt.0) call insert
        endif
        exit
     endif
 
     ii=index(zline,'#')
     if(ii.gt.0) zline(ii:len(zline))=' '  ! ignore comment fields
 
     if(index(zline,'.').eq.0) cycle   ! ignore lines without decimal pts
 
     ii=index(zline,'time=')
     if(ii.gt.0) then
        if(it.gt.0) then
           if(nx_gs2(it).gt.0) call insert
        endif
 
        it=it+1           ! next time point
        nx_gs2(it)=0
 
        gs2x1=0
        gs2r1=0
        gs2ak1=0
        gs2omeg1=0
        gs2gamm1=0
 
     else if(it.gt.0) then
 
        call fields(nflds,cfields,zfields,istat2)
        if(istat2.ne.0) then
           write(lunt,*) ' ?trGS2fetch: decode error, field ',istat2
           write(lunt,*) ' ',zline(1:len_trim(zline))
           go to 1000
        endif
 
        if(nflds.eq.6) then
 
           ! full data line detected
 
           nx_gs2(it)=nx_gs2(it)+1
 
           ix=nx_gs2(it)
           gs2x1(ix)=zfields(2)
           gs2r1(ix)=zfields(3)
           gs2ak1(ix)=zfields(4)
           gs2omeg1(ix)=zfields(5)
           gs2gamm1(ix)=zfields(6)
 
        endif
 
     endif
  enddo
 
  !  it is possible that some times have few or even zero spatial points
  !  warn about less than 2 spatial points; remove entirely if the number
  !  of time points is zero
 
  icur=0
  do it=1,nt_gs2
     if(nx_gs2(it).lt.2) then
        write(lunt,*) ' %trGS2fetch warning:  GS2 data at time ',gs2time(it), &
             ' has less than 2 spatial locations.'
     endif
     if(nx_gs2(it).gt.0) then
        icur=icur+1
        if(icur.lt.it) then
           gs2time(icur)=gs2time(it)
           gs2xa(1:inx,icur)=gs2xa(1:inx,it)
           gs2aky(1:inx,icur)=gs2aky(1:inx,it)
           gs2omega(1:inx,icur)=gs2omega(1:inx,it)
           gs2gamma(1:inx,icur)=gs2gamma(1:inx,it)
        endif
     endif
  enddo
  nt_gs2=icur
  if(icur.eq.0) then
     write(lunt,*) ' ?trGS2fetch error: no valid data found.'
     ier=99
     go to 1000
  endif
 
  !  OK -- now we can copy to new TRANSP variables...
 
  zfab=zprefix(1:ilpre)//'AKY'
  call uupper(zfab)
  ifcn=ifind_ordr(abr,iordrr,nfxt,zfab)
  if(ifcn.gt.0) then
     write(lunt,*) ' %trGS2fetch warning: ',zfab,' replaced.'
     call plcdelfn(zfab,idum)
     CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)  ! check for update
  endif
  call trintrp(inx,nt_gs2,gs2time,gs2aky,gs2xa,ifcnx,inx,datbuf(iwrk2))
  call plsfsave(zfab,'('//zprefix(1:ilpre)//') GS2 AKY',' ',itype,iwrk2,idum)
 
  zfab=zprefix(1:ilpre)//'OMEGA'
  call uupper(zfab)
  ifcn=ifind_ordr(abr,iordrr,nfxt,zfab)
  if(ifcn.gt.0) then
     write(lunt,*) ' %trGS2fetch warning: ',zfab,' replaced.'
     call plcdelfn(zfab,idum)
     CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)  ! check for update
  endif
  call trintrp(inx,nt_gs2,gs2time,gs2omega,gs2xa,ifcnx,inx,datbuf(iwrk2))
  call plsfsave(zfab,'('//zprefix(1:ilpre)//') GS2 OMEGA',' ',itype,iwrk2,idum)
 
  zfab=zprefix(1:ilpre)//'GAMMA'
  call uupper(zfab)
  ifcn=ifind_ordr(abr,iordrr,nfxt,zfab)
  if(ifcn.gt.0) then
     write(lunt,*) ' %trGS2fetch warning: ',zfab,' replaced.'
     call plcdelfn(zfab,idum)
     CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)  ! check for update
  endif
  call trintrp(inx,nt_gs2,gs2time,gs2gamma,gs2xa,ifcnx,inx,datbuf(iwrk2))
  call plsfsave(zfab,'('//zprefix(1:ilpre)//') GS2 GAMMA',' ',itype,iwrk2,idum)
 
1000 continue
 
  if(iopen.eq.1) close(unit=lun_nf)
 
  deallocate(gs2time,stat=istat)
  deallocate(nx_gs2,stat=istat)
  deallocate(gs2x1,gs2r1,gs2ak1,gs2omeg1,gs2gamm1,stat=istat)
  deallocate(gs2xa,gs2aky,gs2omega,gs2gamma,stat=istat)
 
  return
 
  contains
 
    subroutine fields(nflds,cfields,zfields,istat2)
 
      ! parse out the fields from an input line...
 
      integer, intent(out) :: nflds             ! no. of fields found
      character*(*), intent(out) :: cfields(*)  ! character values
      real, intent(out) ::  zfields(*)          ! decoded REAL values
      integer, intent(out) ::  istat2           ! decode error flag
 
      integer i,ic,ic1,ic2,inc,ilz,imode,idcod
      character*20 zbuf
      character*1 ztab
 
      !--------------------------
 
      nflds=0
      imode=0         ! 0: looking for non-blank; 1: looking for blank
 
      do i=1,6
         cfields(i)=' '
         zfields(i)=0
      enddo
 
      istat2=0
 
      ztab=char(9)
 
      ilz=len_trim(zline)+1  ! include a trailing blank
      !  write(6,*) zline(1:ilz)
 
      do ic=1,ilz
 
         if(imode.eq.0) then
            if((zline(ic:ic).ne.' ').and.(zline(ic:ic).ne.ztab)) then
               ic1=ic
               imode=1
            endif
 
         else
            if((zline(ic:ic).eq.' ').or.(zline(ic:ic).eq.ztab)) then
               ic2=ic-1
               imode=0
 
               ! have a field:  decode it (max 6 fields)
 
               nflds=nflds+1
               if(nflds.gt.6) then
                  istat2=7
                  return
               endif
 
               cfields(nflds)=zline(ic1:ic2)
               inc=ic2-ic1+1
               zbuf=' '
               zbuf(20-inc+1:20)=zline(ic1:ic2)
               !  write(6,*) ' "',zbuf,'"'
 
               if(nflds.eq.1) then
 
                  ! integer decode
                  read(zbuf,'(I20)',iostat=istat2) idcod
                  if(istat2.ne.0) then
                     istat2=nflds
                     return
                  endif
                  zfields(nflds)=idcod
 
               else
 
                  ! floating decode
                  read(zbuf,'(G20.0)',iostat=istat2) zfields(nflds)
                  if(istat2.ne.0) then
                     istat2=nflds
                     return
                  endif
 
               endif
            endif   ! blank test
         endif   ! imode test
      enddo
    end subroutine fields
 
    subroutine insert
 
      !  interpolate the GS2 profile at one time
      !  to the TRANSP grid (previously interpolated to that time)
 
      integer ix0,ix0p1,ixtrap
      real zxi,zxic
 
      if(nx_gs2(it).eq.1) then
         gs2aky(1:inx,it)=gs2ak1(1)
         gs2omega(1:inx,it)=gs2omeg1(1)
         gs2gamma(1:inx,it)=gs2gamm1(1)
      else
         do ix=1,inx
 
            ! piecewise linear interpolation
 
            if(zxtest.eq.'R') then
               call xinter(XIDENT,gs2xa(ix,it),gs2r1,nx_gs2(it), &
                    ix0,ix0p1,zxi,zxic,ixtrap)
            else
               call xinter(XIDENT,gs2xa(ix,it),gs2x1,nx_gs2(it), &
                    ix0,ix0p1,zxi,zxic,ixtrap)
            endif
 
            gs2aky(ix,it)=zxic*gs2ak1(ix0)+zxi*gs2ak1(ix0p1)
            gs2omega(ix,it)=zxic*gs2omeg1(ix0)+zxi*gs2omeg1(ix0p1)
            gs2gamma(ix,it)=zxic*gs2gamm1(ix0)+zxi*gs2gamm1(ix0p1)
 
         enddo
      endif
    end subroutine insert
 
end subroutine trgs2fetch
