      subroutine plcfget(zpath,zrunid,zid,zidsave,zshift,istat,ier)
C
C  read a function (could be profile or could be scalar) from a 2nd runid,
C  then, save under a new name
C
      use tconnect_mod
      use datmgr_mod
      use cplotr_mod
C
C  input:
      character*(*) zpath               ! path to 2nd run
      character*(*) zrunid              ! id of 2nd run
      character*(*) zid                 ! fcn inside the 2nd run
      character*(*) zidsave             ! new name to use in this session
      real :: zshift                    ! time shift to apply to 2nd run's time
C                      ! zshift=0.0 usually
C  output:
      integer istat                     ! type of fcn read (scalar/profile vs.)
      integer ier                       ! completion code, 0 = OK
C
C-------------------------------------------------------------------
C  local:
C
      character*40 zdisk
      character*100 zdir
C
C
      character*64 zlbl
      character*32 zuns
      character*10 zxname
C
      real, dimension(:), allocatable :: ztime
      real, dimension(:,:), allocatable :: zfdat,zxdat
C
C-------------------------------------------------------------------
C
      ier0=ier
      ier=0
      istat=0
C
C  get logical i/o unit for messages
C
      lunt=lunzer(0)
C
C  make disk/dir args for trprofil/trscalar
C
      call plcfget_dd(zpath,zdisk,zdir,ier)
      if(ier.ne.0) return
C
C WORKSPACES
      CALL DMDLOC('%WRK1',IND1,ISIZ1,IWRK1)
      CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)
C
C  try reading scalar first
C
      call zermsg(' %plcfget:  runid='//zrunid//'  fcn id='//zid)
C
C  determine type of function
C
      ier=-99                           ! if opened, leave tree open
      call trfunid(zdisk,zdir,zrunid,zid,itype,ier)
      if(ier.ne.0) return
C Check if pppl_transp_open changed server
      if (zpath(1:4) .eq. 'MDS+') then
         il=len_trim(fdisk_x(tc_idrun))
         if (zdisk .ne. fdisk_x(tc_idrun)(1:il)) 
     >       zdisk=fdisk_x(tc_idrun)(1:il)
      endif

C
      if(itype.gt.0) go to 100
C
C  scalar read
C
      allocate(ztime(NTIME))
      call trscalar(zdisk,zdir,zrunid,zid,
     >   NTIME,zlbl,zuns,ntimes,ztime,datbuf(iwrk1),ier)

      ztime = ztime  + zshift
C
      if(ier.ne.0) then
C  error in trscalar, other than scalar "zid" not found
         call zermsg(' %plcfget:  read failed.')
         deallocate(ztime)
         return
      endif
C
C  OK we got a scalar
C
      call plftmk4(ztime,datbuf(iwrk1),ntimes,zlbl,zuns,zidsave)
      ICALL = 3
      CALL DMGFOTX_ww(ICALL, IPT, IERR,iwrk1,iwrk2)
                                ! TAKE SCALARS FROM TEMPORARY AREA
                                ! AND PUT IN F(T) AREA IN MEMORY.
      IF (IERR .NE. 0) THEN
         CALL ZERMSG(' ?plcfget: ERROR RETURNED FROM DMGFOTX')
         ier=1
      ENDIF
C
C  load accumulator
C
      istat=-1
      ifcn=ifind_ordr(abt,iordrt,nft,zidsave)
      do it=1,ntt
         datbuf(iwrk2+it-1)=datbuf(ipt+ntt*(ifcn-1)+it-1)
      enddo
      ttagt(ifcn)=zshift
C
      deallocate(ztime)
      return
C
C------------------------------
C  try for a profile...
C
 100  continue
C
C  get x axis information
C
      call trinfg2(zdisk,zdir,zrunid,zid,zlbl,zuns,
     >   inxnew,intimes,zxname,ixref,itype,inx,ier)
      if(ier.ne.0) go to 1000
C
C  read data
C
      allocate(ztime(intimes))
      allocate(zfdat(inxnew,intimes),zxdat(inxnew,intimes))
      call trprofx(zdisk,zdir,zrunid,zid,zxname,
     >   inxnew,intimes,ztime,zfdat,zxdat,ier)
      if(ier.ne.0) go to 1000

      ztime = ztime  + zshift
C
C  OK time interpolate the data
C
      call trintrp(inxnew,intimes,ztime,zfdat,zxdat,ixref,inx,
     >   datbuf(iwrk2))
C
      istat=itype
C
C  save under the given new name
C
      ier=ier0                          ! can suppress name check...
      call plsfsave(zidsave,zlbl,zuns,itype,iwrk2,ier)
      if(ier.eq.0) then
         ifcn=ifind_ordr(abr,iordrr,nfxt,zidsave)
         ttagr(ifcn)=zshift
      endif
C
 1000 continue
C
      deallocate(zfdat,stat=idum)
      deallocate(zxdat,stat=idum)
      deallocate(ztime,stat=idum)
C
      return
      end
