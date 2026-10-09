subroutine trx_gen_state(lun_geq,geqdsk_lbl,tokdev,runid, &
     inzones,spath,iws,itrx,ier)

  ! perform trx_gen_state_obj1 call on "ps" object in plasma_state_mod
  ! for further description of calling arguments see trx_gen_state_obj1

  use plasma_state_mod
  use trx_module
  implicit NONE

  !-----------------------------------------

  integer, intent(in) :: lun_geq     ! lun for G-EQDSK (0 to suppress)
  character*(*), intent(in) :: geqdsk_lbl ! label for G-EQDSK

  character*(*), intent(in) :: tokdev,runid

  integer, intent(in) :: inzones     ! #radial zones

  character*(*), intent(in) :: spath ! output file path
  integer, intent(in) :: iws,itrx   ! data I/O options (see trx_gen_state_obj1)

  integer, intent(out) :: ier

  !-----------------------------------------

  call trx_gen_state_obj1(ps, &
       lun_geq,geqdsk_lbl,tokdev,runid, &
       inzones,spath,iws,itrx,ier)

end subroutine trx_gen_state

subroutine trx_gen_state_obj(ss,lun_geq,geqdsk_lbl,spath,iws,ier)

  ! trx_gen_state call but with state object passed, and, otherwise,
  ! a simplified calling argument list

  ! perform trx_gen_state_obj1 call on passed object
  ! for further description of calling arguments see trx_gen_state_obj1

  use plasma_state_mod
  use trx_module
  implicit NONE

  !-----------------------------------------

  type (plasma_state) :: ss           ! the state object to be filled...

  integer, intent(in) :: lun_geq     ! lun for G-EQDSK (0 to suppress)
  character*(*), intent(in) :: geqdsk_lbl ! label for G-EQDSK

  character*(*), intent(in) :: spath ! output file path
  integer, intent(in) :: iws         ! data output options

  integer, intent(out) :: ier

  !-----------------------------------------
  !  local:

  integer :: inzones,itrx,idx,inrho
  integer :: lunzer,nonlin

  integer :: istyle,irparen,ilparen,ilblank

  character*4 :: tokdev
  character*40 :: zrun_label
  character*20 :: runid
  
  !-----------------------------------------

  ier = 0
  nonlin = lunzer(0)

  ! to fill in missing arguments for trx_gen_state_obj1...

  !------------------------------------
  !  (1) always use available data
  itrx=1

  !------------------------------------
  !  (2) use TRANSP specified grid size for radial grid; this is in xplasma now

  call eq_ganum('__RHO',idx)
  if(idx.le.0) then
     ier=1
     write(nonlin,*) ' ?trx_gen_state_obj: __RHO not found in xplasma. '
     return
  endif

  call eq_ngrid(idx,inrho)
  inzones = inrho - 1

  !------------------------------------
  !  (3) get labeling...

  call tget_rlbl(tokdev,zrun_label)

  ilparen=index(zrun_label,'(')

  if(ilparen.gt.0) then
     istyle=1
     irparen=index(zrun_label,')')
     if(irparen.le.0) then
        write(nonlin,*) ' %zrun_label:  no closing parentheses: ',zrun_label
        irparen=max(ilparen+2,len(trim(zrun_label)))
     endif
     runid=zrun_label(1:max(1,ilparen-1))
  else
     istyle=2
     ilblank=index(zrun_label,' ')
     runid=zrun_label(ilblank+1:len(trim(zrun_label)))
  endif

  !------------------------------------

  call trx_gen_state_obj1(ss, &
       lun_geq,geqdsk_lbl,tokdev,runid, &
       inzones,spath,iws,itrx,ier)

end subroutine trx_gen_state_obj

subroutine trx_gen_state_obj1(ss,lun_geq,geqdsk_lbl,tokdev,runid, &
     inzones,spath,iws,itrx,ier)

  use plasma_state_mod
  use trx_module
  implicit NONE

  !  generate PLASMA STATE file (ala SWIM Fusion Simulation Project)
  !  dmc Nov. 2006.  An ASCII G-EQDSK file is also generated.

  !=============================================================
  !  Summary of logical controls:

  !    lun_geq > 0 -- write GEQ on this fortran LUN
  !                   Plasma state contains Psi(R,Z) and core equilibrium
  !    lun_geq = 0 -- no GEQ output; state to contain 1d profiles only!

  !    inzones -- number of radial zones

  !    iws = 0 -- output of GEQ file or no file
  !        = 1 -- output of state file with or without GEQ file
  !        = 2 -- output of state file with or without GEQ file, and,
  !               machine description and shot configuration namelists

  !    itrx= 1 -- use all available input data
  !    itrx= 2 -- use TRANSP "rplot" output data only

  !=============================================================

  type (plasma_state) :: ss           ! the state object to be filled...

  integer, intent(in) :: lun_geq      ! lun for ASCII G-EQDSK file
  !  set =0 to suppress G-EQDSK output -- in this case write a skinny
  !  state with no equilibrium

  character*(*), intent(in) :: geqdsk_lbl   ! label for G-EQDSK file
  !  (ignored, if lun_geq=0)

  character*(*), intent(in) :: tokdev,runid ! TOKAMAK and RUNID labels

  integer, intent(in) :: inzones      ! number of zones in state profiles

  character*(*), intent(in) :: spath  ! name (or path) of state file to write
  ! NOTE: filename extension ignored.  I.e. if spath = "/a/b/c/foo.cdf" then
  !       the path information "/a/b/c/" and the filename root "foo" are used
  !       but the extension ".cdf" is ignored.  Standard filename extensions
  !       are used for all output files.
  ! If filename is "%ECH" or "/a/b/c/%ECH", extra files, helpful to TORAY
  ! GENRAY runs, are written in the indicated path or $cwd if no path is given.

  integer, intent(in) :: iws          ! =1 to write state file; O.W. don't.
  !  note: G-EQDSK gets written regardless of value of iws, unless lun_geq=0
  !  Note ALSO:  set iws=2 to write machine description and shot configuration
  !  namelists as well as the state file & G-EQDSK

  !  data source options:
  integer, intent(in) :: itrx         ! =1: use all available data
  !         =2: use TRANSP output data only; do not use namelist or
  !             trdatbuf data; retain data from old state.

  integer, intent(out) :: ier         ! status code on exit (0=OK)

  !---------------------------------------------
  character*200 fpath
  character*40 froot,gname
  integer :: ics,ii,ic,ilen,icdot,id_pmhd,id_q,id1,lunzer,jtrx
  integer :: iflag_ech,iflag_lh
  character*32 zunits,ztest
  real*8 :: zcur
  !---------------------------------------------
  !  for G-EQDSK filename from state filename:  replace .ext with .geq, or,
  !  add the extension .geq

  iflag_ech=0
  iflag_lh=0
  jtrx=itrx

  fpath = spath
  call ptrim(fpath,froot,ics)
  gname = ' '

  ztest = froot
  call uupper(ztest)
  if(ztest.eq.'%ECH') then
     write(lunzer(0),*) ' %trx_gen_state: %ECH detected... '
     iflag_ech=1
     jtrx=1
     froot='cur_state.cdf'
  endif

  ilen=len(trim(froot))

  icdot=0
  do ic=ilen,1,-1
     if(froot(ic:ic).eq.'.') then
        icdot=ic
        exit
     endif
  enddo

  if(icdot.eq.0) then
     gname = trim(froot) // '.geq'
  else
     gname = froot(1:icdot) // 'geq'
     froot(icdot:) = ' '
  endif

  ! get additional data needed for GEQ file:

  call trx_scal('PCUR',zunits,zcur,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state: PCUR not found.'
     return
  endif

  call eq_gfnum('pmhd',id_pmhd)
  call eq_gfnum('q',id_q)

  !-------------------------
  ! write GEQ file...

  if(fpath.eq.' ') then
     if(lun_geq.gt.0) then
        write(lunzer(0),*) ' GEQ FILE:   ',trim(froot)//'.geq'
     endif
     if(iws.gt.0) then
        write(lunzer(0),*) ' STATE FILE: ',trim(froot)//'.cdf'
     endif
     if(iws.gt.1) then
        write(lunzer(0),*) ' MACHINE DESCR: ',trim(froot)//'.mdescr'
        write(lunzer(0),*) ' SHOT CONFIG:   ',trim(froot)//'.sconfig'
     endif
     if(lun_geq.gt.0) then
        open(unit=lun_geq,file=trim(froot)//'.geq',status='unknown',iostat=ier)
     endif
  else
     if(lun_geq.gt.0) then
        write(lunzer(0),*) ' GEQ FILE:   ',trim(fpath)//'/'//trim(froot)//'.geq'
     endif
     if(iws.gt.0) then
        write(lunzer(0),*) ' STATE FILE: ',trim(fpath)//'/'//trim(froot)//'.cdf'
     endif
     if(iws.gt.1) then
        write(lunzer(0),*) ' MACHINE DESCR: ',trim(fpath)//'/'//trim(froot)//'.mdescr'
        write(lunzer(0),*) ' SHOT CONFIG:   ',trim(fpath)//'/'//trim(froot)//'.sconfig'
     endif
     if(lun_geq.gt.0) then
        open(unit=lun_geq,file=trim(fpath)//'/'//trim(froot)//'.geq', &
             status='unknown',iostat=ier)
     endif
  endif

  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state: GEQ file open failure.'
     return
  endif

  if(lun_geq.gt.0) then
     call eq_geqdsk_ipq(lun_geq,geqdsk_lbl,zcur,id_pmhd,id_q,ier)
     if(ier.eq.0) then
        close(unit=lun_geq)
     else
        close(unit=lun_geq,status='DELETE')
        write(lunzer(0),*) ' ?trx_gen_state: GEQ file write failure.'
        return
     endif
  endif

  !-------------------------
  !  init and (optionally) write the state...

  call trx_gen_state_geq(ss, &
       (lun_geq.gt.0), &
       inzones,fpath,froot,tokdev,runid,iws,jtrx,iflag_ech,iflag_lh, &
       ier)

  return

  contains

    !  find tail of file path...

    subroutine ptrim(path,rootname,ics)

      character*(*) path,rootname
      integer :: ics

      ! remove prepended path from filename, if any...

      integer :: ilen,ic,is

      ilen=len(trim(path))
      rootname = " "

      ics=0
      do ic=ilen,1,-1
         if(path(ic:ic).eq.'/') then
            ics=ic+1
            exit
         endif
      enddo

      if(ics.eq.0) then
         ! empty path; input is filename only
         rootname = path
         path = ""
      else

         is=0
         do ic=ics,ilen
            is=is+1
            rootname(is:is)=path(ic:ic)
         enddo

         ! strip trailer from path, including last "/"
         path(ics-1:) = " "
      endif

    end subroutine ptrim

end subroutine trx_gen_state_obj1

!=============================================================

subroutine trx_gen_state_geq(ss, igeq, &
     inzones,fpath,froot,tokdev,runid,iws,itrx,iflag_ech,iflag_lh, &
     ier)

  !  generate PLASMA STATE file (ala SWIM Fusion Simulation Project)
  !  dmc Nov. 2006; G-EQDSK file exists and is given in "gpath" argument.

  use plasma_state_mod
  use trx_module
  use trdatbuf_module

  implicit NONE

  !------------------------
  ! note on MKS conversion: trxplib routines do this already...
  ! so further conversion of RPLOT/TRANSP time dependent data not needed
  ! but for namelist data it is not automatic so an MKS conversion is
  ! needed in some cases.
  !------------------------

  type (plasma_state) :: ss           ! the state object to be filled...

  logical, intent(in) :: igeq         ! .TRUE. to read GEQDSK file
  !  if .TRUE., the output state contains Psi(R,Z) and a core equilibrium
  !  {R,Z}(rho,theta) with equal arc theta.  If .FALSE. a "skinny" state is
  !  written containing 1d profiles only.

  integer, intent(in) :: inzones      ! desired #zones in radial state profiles

  character*(*), intent(in) :: fpath  ! path to directory where to write files
  !  (if blank write in current working directory)

  character*(*), intent(in) :: froot  ! filename root for output files
  !  full output filenames are of form trim(froot)//'.cdf' ...or...
  !                        trim(fpath)//'/'//trim(froot)//'.cdf'.

  character*(*), intent(in) :: tokdev, runid  ! TOKAMAK and RUNID labels

  integer, intent(in) :: iws          ! =1 to write state file; O.W. don't.
  !  if iws.ne.1, state is initialized but not written.
  !  state filename is reset to spath, regardless of value of (iws).
  !  NOTE: if iws=2 also write machine description and shot configuration
  !        namelist files

  !  data source options:
  integer, intent(in) :: itrx         ! =1: use all available data
  !         =2: use TRANSP output data only; do not use namelist or
  !             trdatbuf data.

  !  extra outputs options:
  integer, intent(in) :: iflag_ECH    ! =1 for extra output for ECH models
  integer, intent(in) :: iflag_LH     ! =1 for extra output for LH models

  integer, intent(out) :: ier         ! status code on exit (0=OK)
  !------------------------

  type (pwrget) :: zpwr

  integer n_species,n_thi,n_thx,n_bi,n_rfi,n_fusi
  integer n_fast,n_therm,ii,inum,inumtot,intrace
  integer :: lunzer
  integer :: i,j,ic,ib,ibc,ib0
  integer :: ith,jth,jig,idum,ix,inx,inth,id1,ia,iz,id_g,iertmp
  integer :: jrf,itok,izimp1,izimp2,izatom1,izatom2,iaimp1,iaimp2,iaux,id_q
  integer :: jnbi,inbi,infi,jfus,izmax_th,isize,isign
  integer :: izth,iath,in0,jsc0,jig0,jx

  logical :: iflag,isw_test
  integer :: iwarn,inb_trdat,irf_trdat,iec_trdat,ilh_trdat
  integer :: ilun_tf   ! tmp file lun
  character*50 zftmp0  ! tmp filename root
  character*200 zftmp,zftmp1  ! tmp filenames
  integer :: iltmp     ! tmp filename length

  character*200 fullpath,fullpath_geq,fullpath_mdescr,fullpath_sconfig,tmpath
  integer :: icdf

  integer, dimension(:), allocatable :: itype,iZc
  character*20, dimension(:), allocatable :: slbl
  real*8 :: Zi,zvsur,ziz,zeffinc,zf1,zf2,zlim1,zlim2,zfrac,ztest,zconv
  real*8 :: zdenom,zavg,zdatum
  real*8 :: zamax,zzmax,zzmin
  real*8 :: ztime0,zdelta_t,zitem,znloc
  integer, parameter :: max_zmbuf=1000
  real*8 :: zmbuf(max_zmbuf)
  integer :: imap0(max_zmbuf)

  integer :: igot,jgot,izloc,ig

  real*8, dimension(:), allocatable :: znxsum,zqsave
  real*8, dimension(:), allocatable :: Zcharg,Amass,zgrid,zgridc,zprof,zprofc
  logical, dimension(:), allocatable :: izprof
  real*8, dimension(:,:), allocatable :: zcharga
  real*8, dimension(:), allocatable :: zprofc1,zprofc2,zprofc3,zprofc4,Zia
  real*8, dimension(:), allocatable :: zvol,zarea,zzimp,zaimp,abeama,xzbeama
  real*8, dimension(:), allocatable :: zrminor
  real*8, dimension(:), allocatable :: gprofc,omegac,zvpllc,sgas,srcy
  real*8, dimension(:), allocatable :: zzbuf,ffulla,fhalfa,zspectrum
  real*8, dimension(:,:), allocatable :: kvfrac
  real*8 :: kvminm(3),sgrcy_sum,zsc
  real*8 :: zemmx,zemmx_std(2),zemmx_neg(2)
  real*8 :: t0recyc,t0gasfl,om0recyc,om0gasfl,zRedge
  real*8, dimension(:), allocatable :: e0in
  real*8, dimension(:,:), allocatable :: zn0

  character*1, dimension(:), allocatable :: cgas_abbrev

  integer :: ibminm(3),iemmx(2)
  integer, dimension(:,:), allocatable :: ibi
  integer, dimension(:), allocatable :: idns,idts,id_eprps,id_eplls
  integer, dimension(:), allocatable :: ianbi,iznbi,iwk
  integer, dimension(:), allocatable :: ianfusi,iznfusi
  
  logical, dimension(:), allocatable :: ilco
  integer, dimension(:), allocatable :: ibmap  ! th species -> beam species map
  integer, dimension(:), allocatable :: ifmap  ! th species -> fusn species map
  integer :: ishap,indx
  logical :: iomega
  
  integer :: ibeamx,ibeami,ifusx,ifusi,ival(1)
  integer :: imj1,imj2,ix1,ix2  ! ion specie category indices
  integer :: ivtr1,ivtr2        ! ion species w/VTOR data
  integer :: ivpl1,ivpl2        ! ion species w/VPOL data

  integer :: istat,imatch
  integer,dimension(:), allocatable :: inb_sublist

  integer :: nicrf
  logical :: nltoray,nlgen_ech,nltorbeam,ec_aim_data

  integer :: levgeo,nsomod,nmdifb(1)

  character*3, dimension(:), allocatable :: thsuffix  ! therm. species suffix
  character*1, dimension(:), allocatable :: thsuff_1  ! therm. species 1 letter
  character*1, dimension(:), allocatable :: bsuffix   ! beam species suffix
  character*1, dimension(:), allocatable :: fsuffix   ! fusion species suffix
  character*1 :: bsuff0,fsuff0,thsuff0   ! 1st non-blank suffix of each type
  character*10 profname

  integer, dimension(:), allocatable :: ntrace       ! trace element beam info
  integer, dimension(:), allocatable :: imap,imapi

  real*8, dimension(:), allocatable :: ftrace,ftraci ! trace fractions

  !  support trace element beams:
  !    in TRANSP indexing these count as separate beams; in Plasma State
  !    these are not counted separately but treated as optional attributes.
  !
  !  for ntrace(ib).gt.0, ib0=abs(ntrace(ib)) is the main TRANSP beam index
  !  ib0=imap(ib) -> ib0 is zero if TRANSP index ib is a trace beam, otherwise
  !      ib0 is the main beam number in the plasma state; ib=1:ss%nbeam+ntrace
  !  ibc=imapi(ib) -> TRANSP index to ib'th non-trace beam; ib=1:ss%nbeam

  integer :: iorder = 1   ! using Hermite for experimental data profiles...

  character*32 zunits,zname
  character*64 zlabel

  logical :: erase_state,merge_lists
  logical :: splitn_avail,splitn_try
  logical :: trdatbuf_avail,trdatbuf_try
  logical :: iexist,jexist,rsn_exist,rs2_exist
  logical :: all_neutrals,reco_neutrals
  logical :: ptransp_flag

  logical :: iflg_pwr,iflg_vlt,iflg_frq
  logical :: iflg_full,iflg_half

  real*8 :: zRmin,zRmax,zYmin,zYmax
  real*8 :: zpcur,zpcurc,zpcur_diff
  real*8 :: znebar

  real*8, dimension(:), allocatable :: zRlim,zZlim
  integer :: imax_lim,inum_lim

  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE = 1.0d0
  real*8, parameter :: TWOPI = 6.283185307179586d0
  real*8, parameter :: ZSMALL= 1.0d-7

  real*8, dimension(1) :: ztemp
  integer, dimension(1) :: itemp
  logical, dimension(1) :: ltemp

  !----------------------------------------------------
  splitn_try = .TRUE.
  trdatbuf_try = .TRUE.
  erase_state = .TRUE.
  merge_lists = .TRUE.

  if((itrx.eq.2).AND.(max(iflag_lh,iflag_ech).le.0)) then
     ! test_client mode, only use TRANSP output record as input here
     splitn_try = .FALSE.
     trdatbuf_try = .FALSE.
     erase_state = .FALSE.
     merge_lists = .FALSE.
  endif

  ptransp_flag = .FALSE.   ! assume, for now...

  if(fpath.eq.' ') then
     fullpath = trim(froot)//'.cdf'
     fullpath_geq = trim(froot)//'.geq'
     fullpath_mdescr = trim(froot)//'.mdescr'
     fullpath_sconfig = trim(froot)//'.sconfig'
  else
     fullpath = trim(fpath)//'/'//trim(froot)//'.cdf'
     fullpath_geq = trim(fpath)//'/'//trim(froot)//'.geq'
     fullpath_mdescr = trim(fpath)//'/'//trim(froot)//'.mdescr'
     fullpath_sconfig = trim(fpath)//'/'//trim(froot)//'.sconfig'
  endif

  !----------------------------------------------------
  write(lunzer(0),*) ' '
  write(lunzer(0),*) ' %trx_gen_state: building state from TRANSP data.'
  if(iflag_lh.gt.0) write(lunzer(0),*) ' %trx_gen_state: iflag_LH is set.'
  if(iflag_ech.gt.0) write(lunzer(0),*) ' %trx_gen_state: iflag_ECH is set.'
  write(lunzer(0),*) ' '
  !----------------------------------------------------
  !  grid sizes (TRANSP: equilibrium and plasma grids the same)

  call ps_init_tag  ! be sure to init state module

  if(erase_state) then
     call trx_init_state_obj(ss,ier)  ! be sure state is empty...
     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' %trx_gen_state: trx_init_state error code: ',ier
        go to 1000
     endif
  endif

  ss%tokamak_id = tokdev
  ss%runid = runid

  ss%geometry = ' '
  write(ss%geometry,1001) trim(tokdev),trim(run_label),time0,delta_t
1001 format("TRANSP(",a,",'",a,"'; time0=",1pe13.6,"s; delta_t=",1pe13.6,"s)")

  ztime0 = time0     ! real*8
  zdelta_t = delta_t ! real*8

  inx=inzones+1
  ss%nrho=inx
  ss%nrho_eq=inx
  ss%nrho_rad=inx
  ss%nrho_gas=inx

  if(igeq) then
     call eq_ngrid(id_chi,inth)
     ss%nth_eq=max((2*inzones+1),inth)

     !  EQDSK file (path not included)

     ss%eqdsk_file = trim(froot)//'.geq'

  else

     ss%eqdsk_file = 'NONE'

  endif

  !  times: set t0 = time0-delta_t; t1 = time0+delta_t
  write(lunzer(0),*) ' ---------------------------- '
  write(lunzer(0),'(2x,"Extraction time: ",1pe11.4," +/- ",1pe11.4," secs")') &
       time0,delta_t
  write(lunzer(0),*) ' ---------------------------- '

  ss%t0 = time0-delta_t
  ss%t1 = time0+delta_t

  !----------------------------------------------------
  !  get species counts

  call rd_nspecies(n_species,n_thi,n_thx,n_bi,n_rfi,n_fusi)
  !----------------------------------------------------
  !  get all plasma species

  allocate(itype(n_species),iZc(n_species),ibmap(n_species),ifmap(n_species))
  allocate(slbl(n_species))
  allocate(Zcharg(n_species),Amass(n_species))
  allocate(zcharga(inx-1,n_species),izprof(n_species)); izprof=.FALSE.
  allocate(Zia(inx-1))

  allocate(thsuffix(n_species),bsuffix(n_species),fsuffix(n_species))
  thsuffix='   '; bsuffix=' '; fsuffix=' '
  allocate(thsuff_1(n_species))
  thsuff_1=' '
  thsuff0 =' '; bsuff0=' '; fsuff0=' '

  call trx_spec_lbl(n_species,idum,slbl,itype,Zcharg,Amass,iZc)

  ! trx_gen_state code assumes an ordering of species types: verify...
  call check_itype_order(ier)
  if(ier.ne.0) go to 1000

  ! NOTE: the code below assumes that non-impurity thermal ions come
  ! before impurity thermal ions, in the species list returned.
  ! This behaviour verified in trx_spec_lbl, dmc Dec 2007...

  !----------------------------------------------------
  !  see if TRANSP "tokamakium" is present; if so, increment #therm by 1
  !   this allows amalgum impurity with non-integer Z to be split into
  !   two species with integer Zs while preserving Zeff and quasineutrality

  itok=0
  do i=1,n_species
     if(itype(i).eq.ps_tokamakium) then
        itok=1
     endif
  enddo

  n_therm = n_thi + n_thx + itok  ! bulk plasma and impurity thermal ions
  n_fast = n_bi + n_rfi + n_fusi  ! beam ions, RF minority ions, fusion ions

  imj1=1
  imj2=n_thi
  ix1=imj2+1
  ix2=imj2+n_thx+itok

  !  copy total thermal and non-thermal species counts into state & allocate

  ss%nspec_beam = n_bi

  ss%nspec_rfmin = n_rfi
  ss%kdens_rfmin = "data" ! RF minority density profiles (if any) used directly

  ss%nspec_fusion = n_fusi

  ss%nspec_th = n_therm

  !----------------------------------------------------
  ! set some code_info labels here (no other options in TRANSP at present:
  ! DMC Apr 2010).

  if(n_bi.gt.0) then
     ss%nbi_code_info = 'TRANSP:NUBEAM'
     ss%nbi_data_info = 'TRANSP:input_data'
  endif

  if(n_fusi.gt.0) then
     ss%fus_code_info = 'TRANSP:NUBEAM'
  endif

  !----------------------------------------------------
  !  initial allocate of state arrays

  call ps_alloc_plasma_state(ier, state=ss)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state_geq: ps_alloc_plasma_state status: ',&
          ier
     go to 1000
  endif

  !----------------------------------------------------
  !  insert radial grid for plasma profiles.
  !  equilibrium grids are handled below by the ps_update_equilibrium call

  allocate(zgrid(inx),zgridc(inx-1),zprof(inx),zprofc(inx-1))
  allocate(zprofc1(inx-1),zprofc2(inx-1),zprofc3(inx-1),zprofc4(inx-1))
  allocate(zvol(inx),zarea(inx),zzimp(inx-1),zaimp(inx-1))
  allocate(zrminor(inx))

  do ix=1,inx
     zgrid(ix)=(ix-1)*ONE/(inx-1)
     if(ix.gt.1) zgridc(ix-1)=(zgrid(ix)+zgrid(ix-1))/2
  enddo
  ss%rho = zgrid
  ss%rho_rad = zgrid
  ss%rho_gas = zgrid

  ! test if TRANSP data includes the full description of neutral gas sources
  call rpexist_profile('BALN0_SRC',all_neutrals)
  call rpexist_profile('SERECO',reco_neutrals)

  if(n_bi.gt.0) then
     call mk_rho_nbi
     allocate(ianbi(n_bi),iznbi(n_bi))
  endif
  if(n_rfi.gt.0) call mk_rho_icrf
  if(n_fusi.gt.0) then
     allocate(ianfusi(n_fusi),iznfusi(n_fusi))
     call mk_rho_fus
  endif
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state_geq: ps_alloc_plasma_state status: ',&
          ier
     go to 1000
  endif

  !  <A> and <Z> profiles if tokamakium present in TRANSP run

  if(itok.eq.1) then
     call get1('AIMPJ')
     if(ier.ne.0) then
        write(lunzer(0),*) ' %trx_gen_state:  OK, fall back to AIMP scalar...'
        call trx_scal('AIMP',zunits,zdatum,ier)
        if(ier.ne.0) go to 1000
        zprofc = zdatum
     endif

     zaimp = zprofc

     call get1('XZIMPJ')
     if(ier.ne.0) then
        write(lunzer(0),*) ' %trx_gen_state:  OK, fall back to XZIMP scalar...'
        call trx_scal('XZIMP',zunits,zdatum,ier)
        if(ier.ne.0) go to 1000
        zprofc = zdatum
     endif

     zzimp = zprofc

     zamax=ZERO
     zzmax=ZERO
     zzmin=10000.0d0

     do ix=1,inx-1
        zamax=max(zamax,zaimp(ix))
        zzmax=max(zzmax,zzimp(ix))
        zzmin=min(zzmin,zzimp(ix))
     enddo

     izimp1=zzmin
     izimp2=zzmax+0.99d0
     if(izimp2.eq.izimp1) izimp2=izimp1+1

     iaimp1 = izimp1*(zamax/zzmax) + 0.5d0
     iaimp2 = izimp2*(zamax/zzmax) + 0.5d0

     izatom1 = max(izimp1,(iaimp1/2))
     izatom2 = max(izimp2,(iaimp2/2))

     !  zzimp will be constrained to this range to prevent a zero density...

     zlim1 = izimp1 + 1.0d-6
     zlim2 = izimp2 - 1.0d-6

  endif

  allocate(zvpllc(inx-1),gprofc(inx-1))
  allocate(zqsave(inx))

  !  q profile
  call eq_gfnum('q',id_q)
  call eq_rgetf(inx,zgrid,id_q,0,zprof,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?? trx_gen_state: no q profile ??'
     ss%iota = ZERO
     zqsave = ZERO
  else
     ss%iota = ONE/zprof(1:inx)
     zqsave(1:inx) = zprof(1:inx)
  endif

  !  g profile

  call eq_gfnum('G',id_g)
  call eq_rgetf(inx-1,zgridc,id_g,0,gprofc,ier)
  if(ier.ne.0) go to 1000

  !  velocity data

  ss%vtor_omp = ZERO
  ss%vpol_omp = ZERO

  ss%vtor_inmp = ZERO
  ss%vpol_inmp = ZERO

  ivtr1=-1
  ivtr2=-1
  ivpl1=-1
  ivpl2=-1

  call rpexist_profile('OMEGA',iexist)  ! check for toroidal angular velocity

  if(.not.iexist) then
     write(lunzer(0),*) &
          ' %trx_gen_state: no angular velocity found; assuming ZERO rotation.'
     zvpllc=0
     ss%omegat = 0
     iomega=.FALSE.

  else

     call get1('OMEGA')
     ss%omegat = zprofc
     zvpllc=zprofc  ! this is the angular velocity; it will be converted...
     iomega=.TRUE.

     call get1('GBR2')  ! new runs will have computed <B*R**2>...
     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' ?trx_gen_state: <B*R**2> not found.'
        write(lunzer(0),*) &
             '  use approximation: <v_pll> = omega/<1/R>.'
        call get1('GRI')
        if(ier.ne.0) then
           write(lunzer(0),*) &
                ' ?trx_gen_state: <1/R> not found; setting <v_pll> = 0.0'
           zvpllc = 0
        else
           zvpllc = zvpllc/zprofc
        endif
     else
        !  <vpll> = <(B/B_phi)*v_phi> = (1/g)*<R*B*v_phi> = (omega/g)*<B*R**2>
        zvpllc = zvpllc * zprofc/gprofc
     endif

  endif

  ss%rho_bdy_omegat = ONE
  call set_bdy(ss%omegat,ss%omegat_bdy)

  !  electrostatic potential

  call rpexist_profile('VRPOT',iexist)
  if(iexist) then
     call get1b('VRPOT')
     ss%epot = zprof*0.001d0    ! "volts" -> KeV
  endif

  !  resistivity from the TRANSP run

  call get1b('ETA_USE')
  ss%eta_parallel = zprof    ! ohm*cm -> ohm*m

  !-----------------------------------
  !  set species label

  jth = 0  ! for thermal species list
  jig = 0  ! non-impurity thermal species list
  jnbi= 0  ! for beam ions
  jfus= 0  ! fusion product ions
  jrf = 0  ! for RF minority ions
  izmax_th = 0

  ibmap = 0
  ifmap = 0
  do i=1,n_species
     if(itype(i).eq.ps_beam_ion) then
        iz=Zcharg(i)+0.1d0  ! integer charge
        ia=Amass(i)+0.1d0   ! integer AMU
        do j=1,n_species
           if(itype(j).eq.ps_therm_ion) then
              izth=Zcharg(j)+0.1d0  ! integer charge
              iath=Amass(j)+0.1d0   ! integer AMU

              if((iath.eq.ia).and.(izth.eq.iz)) then
                 ! species #j = thermalized beam species #i
                 ibmap(j)=i
              endif
           endif
        enddo
     endif
     if(itype(i).eq.ps_fusion_ion) then
        iz=Zcharg(i)+0.1d0  ! integer charge
        ia=Amass(i)+0.1d0   ! integer AMU
        do j=1,n_species
           if(itype(j).eq.ps_therm_ion) then
              izth=Zcharg(j)+0.1d0  ! integer charge
              iath=Amass(j)+0.1d0   ! integer AMU

              if((iath.eq.ia).and.(izth.eq.iz)) then
                 ! species #j = thermalized fusion species #i
                 ifmap(j)=i
              endif
           endif
        enddo
     endif
  enddo

  do i=1,n_species
     if(itype(i).eq.ps_electron) then
        !  electron
        call ps_species_convert(-1,0,0, ss%qatom_s(0), ss%q_s(0), ss%m_s(0), &
             ier)

     else if((itype(i).eq.ps_therm_ion).or.(itype(i).eq.ps_impurity)) then
        !  thermal ion or impurity
        if(itype(i).eq.ps_therm_ion) then
           jig=jig+1  ! non-impurity
           izmax_th=max(izmax_th,iZc(i))
        endif
        jth=jth+1
        iz=Zcharg(i)+0.1d0  ! integer charge
        ia=Amass(i)+0.5d0   ! integer AMU
        call ps_species_convert(iZc(i),iz,ia, &
             ss%qatom_s(jth), ss%q_s(jth), ss%m_s(jth), ier)

        if(itype(i).eq.ps_therm_ion) then
           call set_suffix3(iz,ia,thsuffix(i))
           call set_suffix1(iz,ia,thsuff_1(i),' ',thsuff0)
        endif

     else if(itype(i).eq.ps_tokamakium) then
        !  TOKAMAKIUM -- impurity amalgum
        !    gets split into two species each with integer Z values
        jth=jth+1
        call ps_species_convert(izatom1,izimp1,iaimp1, &
             ss%qatom_s(jth), ss%q_s(jth), ss%m_s(jth), ier)

        jth=jth+1
        call ps_species_convert(izatom2,izimp2,iaimp2, &
             ss%qatom_s(jth), ss%q_s(jth), ss%m_s(jth), ier)

     else
        !  fast ion 
        iz=Zcharg(i)+0.1d0  ! integer charge
        ia=Amass(i)+0.5d0   ! integer AMU

        if(itype(i).eq.ps_rf_minority) then
           jrf=jrf+1
           call ps_species_convert(iZc(i),iz,ia, &
                ss%qatom_rfmin(jrf), ss%q_rfmin(jrf), ss%m_rfmin(jrf), ier)

        else if(itype(i).eq.ps_beam_ion) then
           jnbi=jnbi+1
           iznbi(jnbi)=iz
           ianbi(jnbi)=ia
           if(iz.eq.0) then
              iz=iZc(i)
              izprof(i)=.TRUE.
              iznbi(jnbi)=iz
           endif
           call ps_species_convert(iZc(i),iz,ia, &
                ss%qatom_snbi(jnbi), ss%q_snbi(jnbi), ss%m_snbi(jnbi), ier)
           call set_suffix1(iz,ia,bsuffix(i),'B',bsuff0)

        else if(itype(i).eq.ps_fusion_ion) then
           jfus=jfus+1
           iznfusi(jfus)=iz
           ianfusi(jfus)=ia
           call ps_species_convert(iZc(i),iz,ia, &
                ss%qatom_sfus(jfus), ss%q_sfus(jfus), ss%m_sfus(jfus), ier)
           call set_suffix1(iz,ia,fsuffix(i),'F',fsuff0)

        endif

     endif

     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' ?trx_gen_state_geq: error during species definition loop.'
        go to 1000
     endif
  enddo

  !  set species labels in state

  ss%Z0max = izmax_th
  call ps_label_species(ier, state=ss)
  if(ier.ne.0) then
     write(lunzer(0),*) &
          ' ?trx_gen_state_geq: ps_label_species status: ', ier
     go to 1000
  endif

  if(merge_lists) then

     !  merge species lists

     call ps_merge_species_lists(ier, state=ss)
     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' ?trx_gen_state_geq: ps_merge_species_lists status: ', ier
        go to 1000
     endif

     !  neutral gas species

     call ps_neutral_species(ier, state=ss)
     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' ?trx_gen_state_geq: ps_neutral_species status: ', ier
        go to 1000
     endif

     jig0=jig
     if(jig0.ne.ss%nspec_gas) then
        write(lunzer(0),*) &
             ' ?trx_gen_state_geq: number of thermal species does not match'
        write(lunzer(0),*) &
             '  number of neutral gas species as expected: ',jig0,ss%nspec_gas,'.'
        ier=1
        go to 1000
     endif

     !  make list of neutral sources per TRANSP style:
     !    two for each gas specie, 1 for "recycling" and 1 for "gas flow"

     call mk_ngsc0
     if(ier.ne.0) go to 1000

  endif

  if(igeq) then
     !----------------------------------------------------
     !  acquire equilibrium within the plasma state (read GEQ file)

     call ps_update_equilibrium(ier, g_filepath=fullpath_geq, state=ss)
     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' ?trx_gen_state_geq: ps_update_equilibrium status: ',ier
        go to 1000
     endif

     !--------------------
     !  supplement equilibrium with computation of flux surface averages and
     !  NCLASS moments

     ss%nrho_eq_geo = ss%nrho_eq
     CALL ps_alloc_plasma_state(ier, state=ss)
     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' ?trx_gen_state_geq: ps_alloc_plasma_state status: ',ier
        go to 1000
     endif

     ss%rho_eq_geo = ss%rho_eq  ! use same grid
     call ps_mhdeq_derive('Everything',ier, state=ss)
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state_geq: ps_mhdeq_derive status: ',&
             ier
        go to 1000
     endif

     !----------------------------------------------------
     !  fetch volumes and areas for local use

     ier=0
     call ps_intrp_1d(ss%rho,ss%id_vol,zvol,iertmp, state=ss)
     ier=max(ier,iertmp)

     call ps_intrp_1d(ss%rho,ss%id_area,zarea,iertmp, state=ss)
     ier=max(ier,iertmp)

     if(ier.ne.0) then
        write(lunzer(0),*) ' ?ps_intrp_1d (volumes and areas) failed, ier=',ier
        go to 1000
     endif

  else
     !  use trxpl old_xplasma volumes and areas

     ier=0
     call eq_volume(inx,ss%rho,0,zvol,iertmp)
     ier=max(ier,iertmp)
     call eq_area(inx,ss%rho,0,zarea,iertmp)
     ier=max(ier,iertmp)

     if(ier.ne.0) then
        write(lunzer(0),*) ' ?xplasma eq (volumes and areas) failed, ier=',ier
        go to 1000
     endif
     
     !  only this 1d MHD eq data is set:

     ss%rho_eq = ss%rho
     ss%vol = zvol
     ss%area = zarea

     call eq_bchk_sign( ss%kccw_Bphi, ss%kccw_Jphi )

     imax_lim=151
     allocate(zRlim(imax_lim),zZlim(imax_lim))

     call eq_limcon(imax_lim,inum_lim,zRlim,zZlim,zsmall,ier)
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?eq_limcon from trx_gen_state_geq: ier=',ier
        go to 1000
     endif

     ss%num_rzlim = inum_lim
     call ps_alloc_plasma_state(ier, state=ss)
     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' ?ps_alloc_plasma_state (for limiter): ier=',ier
        go to 1000
     endif

     ss%rlim = zRlim(1:inum_lim)
     ss%zlim = zZlim(1:inum_lim)

     deallocate(zRlim,zZlim)

     call ps_state_memory_update(ier, state=ss)
     if(ier.ne.0) then
        write(lunzer(0),*) &
             ' ?ps_state_memory_update (volumes and areas): ier=',ier
        go to 1000
     endif

  endif

  !----------------------------------------------------
  !  can now fetch: bootstrap current, etc.

  !  Ohmic heating

  
  call get1p('POH',1)
  ss%pohme = zprofc

  !  Ion electron coupling

  call get1p('QIE',1)
  ss%qie = -zprofc           ! PS sign convention = -[TRANSP convention]

  !  Bootstrap current

  call get1p('CURBS',2)
  ss%curr_bootstrap = zprofc

  !  Ohmic current

  call get1p('CUROH',2)
  ss%curr_ohmic = zprofc
  !
  !----------------------------------------------------
  !  gather charge profile data as needed
  do i=1,n_species
     if(izprof(i)) then
        if(itype(i).ne.ps_beam_ion) then
           write(lunzer(0),*) ' ?trx_gen_state: charge profile for non-beam-specie '
           ier=99
           go to 1000
        endif

        call get1('AVGZ_'//bsuffix(i))
        zcharga(1:inx-1,i)=zprofc

     else
        zcharga(1:inx-1,i)=ZERO
     endif
  enddo

  !----------------------------------------------------
  !  profiles...
  !  Hermites of all densities

  allocate(idns(n_species),idts(n_species))
  allocate(id_eprps(n_species),id_eplls(n_species))

  call trx_spec_prof(n_species,idum,iorder,'N',idns,ier)
  if(ier.ne.0) go to 1000

  !  Hermites of all temperatures (only used for thermal species)

  call trx_spec_prof(n_species,idum,iorder,'T',idts,ier)
  if(ier.ne.0) go to 1000

  !  Hermites of all flux surf avg <Eperp> profiles (only for fast species)

  call trx_spec_prof(n_species,idum,iorder,'<Eperp>',id_eprps,ier)
  if(ier.ne.0) go to 1000

  !  Hermites of all flux surf avg <Epll> profiles (only for fast species)

  call trx_spec_prof(n_species,idum,iorder,'<Epll>',id_eplls,ier)
  if(ier.ne.0) go to 1000

  jth = 0  ! for thermal species list
  jnbi= 0  ! for beam ions
  jfus= 0  ! fusion product ions
  jrf = 0  ! for RF minority ions
 
  ss%zeff = 0
  ss%zeff_fi = 0
  ss%zeff_th = 0

  ss%fi_depletion = 0

  do i=1,n_species
     if(itype(i).eq.ps_electron) then
        !  electron
        call eq_rgetf(inx-1,zgridc,idns(i),0,zprofc,ier)
        if(ier.ne.0) exit
        ss%ns(1:inx-1,ps_elec_index)=zprofc  ! ne
        call eq_rgetf(inx-1,zgridc,idts(i),0,zprofc,ier)
        if(ier.ne.0) exit
        ss%Ts(1:inx-1,ps_elec_index)=zprofc  ! Te

        ss%v_pars(1:inx-1,ps_elec_index)=zvpllc   ! <vpll> the same for all species for now

        ss%rho_bdy_Te = ONE
        call set_bdy(ss%Ts(:,ps_elec_index),ss%Te_bdy)

     else if(itype(i).le.ps_impurity) then
        !  thermal specie
        jth=jth+1
        call eq_rgetf(inx-1,zgridc,idns(i),0,zprofc,ier)
        if(ier.ne.0) exit
        ss%ns(1:inx-1,jth)=zprofc  ! ni
        call eq_rgetf(inx-1,zgridc,idts(i),0,zprofc,ier)
        if(ier.ne.0) exit
        ss%Ts(1:inx-1,jth)=zprofc  ! Ti

        ss%v_pars(1:inx-1,jth)=zvpllc ! <vpll> the same for all species for now

        Zi = Zcharg(i)
        ss%zeff(1:inx-1)=ss%zeff(1:inx-1) + ss%ns(1:inx-1,jth)*Zi*Zi
        ss%zeff_th(1:inx-1)=ss%zeff_th(1:inx-1) + ss%ns(1:inx-1,jth)*Zi*Zi

     else if(itype(i).eq.ps_tokamakium) then
        !  TRANSP tokamakium
        jth=jth+1  ! will increment again...
        call eq_rgetf(inx-1,zgridc,idns(i),0,zprofc,ier)
        if(ier.ne.0) exit
        do ix=1,inx-1
           ziz=max(zlim1,min(zlim2,zzimp(ix)))
           call split_imp(ziz,ziz*ziz,izimp1,izimp2,zf1,zf2,ier)
           if(ier.ne.0) then
              write(lunzer(0),*) ' ?trx_gen_state: split_imp error code: ',ier
              exit
           endif
           ss%ns(ix,jth)  = zf1*zprofc(ix)
           ss%ns(ix,jth+1)= zf2*zprofc(ix)
           zeffinc = zprofc(ix)*(zf1*izimp1*izimp1 + zf2*izimp2*izimp2)
           ss%zeff(ix)=ss%zeff(ix)+zeffinc
           ss%zeff_th(ix)=ss%zeff_th(ix)+zeffinc
        enddo
        if(ier.ne.0) exit

        call eq_rgetf(inx-1,zgridc,idts(i),0,zprofc,ier)
        if(ier.ne.0) exit
        ss%Ts(1:inx-1,jth)=zprofc  ! Ti
        ss%Ts(1:inx-1,jth+1)=zprofc  ! Ti

        ss%v_pars(1:inx-1,jth)=zvpllc ! <vpll> the same for all species for now
        ss%v_pars(1:inx-1,jth+1)=zvpllc ! <vpll> the same for all species

        jth=jth+1
     else
        !  fast specie
        call eq_rgetf(inx-1,zgridc,idns(i),0,zprofc,ier)
        if(ier.ne.0) exit

        if(itype(i).eq.ps_rf_minority) then
           jrf=jrf+1
           ss%nmini(1:inx-1,jrf)=zprofc

        else if(itype(i).eq.ps_beam_ion) then
           jnbi=jnbi+1
           ss%nbeami(1:inx-1,jnbi)=zprofc

        else if(itype(i).eq.ps_fusion_ion) then
           jfus=jfus+1
           ss%nfusi(1:inx-1,jfus)=zprofc
        endif

        if(iZprof(i)) then
           Zia = Zcharga(:,i) ! (impurity beam species <Z>)
        else
           Zia = Zcharg(i)
        endif

        ss%zeff(1:inx-1)=ss%zeff(1:inx-1) + zprofc*Zia*Zia
        ss%zeff_fi(1:inx-1)=ss%zeff_fi(1:inx-1) + zprofc*Zia*Zia
        ss%fi_depletion(1:inx-1)=ss%fi_depletion(1:inx-1) + zprofc*Zia

        if(iZprof(i)) then
           zdenom = sum(zprofc*(zvol(2:inx)-zvol(1:inx-1)))
           if(zdenom.gt.ZERO) then
              Zavg = sum(zprofc*Zia*(zvol(2:inx)-zvol(1:inx-1)))/zdenom
              ss%q_snbi(jnbi) = Zavg*ps_xe
           endif
        endif

        call eq_rgetf(inx-1,zgridc,id_eprps(i),0,zprofc,ier)
        if(ier.ne.0) exit

        if(itype(i).eq.ps_rf_minority) then
           ss%eperp_mini(1:inx-1,jrf)=zprofc    ! <Eperp>RF
        else if(itype(i).eq.ps_beam_ion) then
           ss%eperp_beami(1:inx-1,jnbi)=zprofc  ! <Eperp>NBI
        else if(itype(i).eq.ps_fusion_ion) then
           ss%eperp_fusi(1:inx-1,jfus)=zprofc   ! <Eperp>FUSI
        endif

        call eq_rgetf(inx-1,zgridc,id_eplls(i),0,zprofc,ier)
        if(ier.ne.0) exit

        if(itype(i).eq.ps_rf_minority) then
           ss%epll_mini(1:inx-1,jrf)=zprofc    ! <Epll>RF
        else if(itype(i).eq.ps_beam_ion) then
           ss%epll_beami(1:inx-1,jnbi)=zprofc  ! <Epll>NBI
        else if(itype(i).eq.ps_fusion_ion) then
           ss%epll_fusi(1:inx-1,jfus)=zprofc   ! <Epll>FUSI
        endif

        if(itype(i).eq.ps_rf_minority) then
           zfrac = ss%nmini(1,jrf)/ss%ns(1,ps_elec_index)  ! ni/ne
           ss%fracmin(jrf) = zfrac
           do ix=2,inx-1
              zfrac = ss%nmini(ix,jrf)/ss%ns(ix,ps_elec_index)  ! ni/ne
              ztest = abs(zfrac-ss%fracmin(jrf))/max(zfrac,ss%fracmin(jrf))
              if(ztest.gt.1.0d-4) then
                 write(lunzer(0),*) &
                      ' %trx_gen_state: minority density fraction has radial variation.'
                 ss%fracmin(jrf) = ZERO
                 exit
              endif
           enddo
        endif
     endif
  enddo

  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state: profile lookup error.'
     go to 1000
  endif

  !  Zeff = sum[ni*Zi*Zi]/ne
  ss%zeff = ss%zeff/ss%ns(1:inx-1,ps_elec_index)
  ss%zeff_th = ss%zeff_th/ss%ns(1:inx-1,ps_elec_index)
  ss%zeff_fi = ss%zeff_fi/ss%ns(1:inx-1,ps_elec_index)
  ss%fi_depletion = ss%fi_depletion/ss%ns(1:inx-1,ps_elec_index)

  !  velocity data

  if(iomega) then
     call get_vtor_vpol
  endif

  !  compute single averaged Ti and total thermal ion density
  call ps_ti_fetch(ier, state=ss, ni=ss%ni, Ti=ss%Ti)

  !  set remaining boundary data
  
  ss%rho_bdy_Ti = ONE
  call set_bdy(ss%Ti,ss%Ti_bdy)

  do i=0,ss%nspec_th
     ss%rho_bdy_ns(i) = ONE
     call set_bdy(ss%ns(:,i),ss%ns_bdy(i))
  enddo


  !---------------------------------------------------
  ! line averaged electron density
  !
  if(igeq) then
    do ix = 1, inx
     zRminor(ix) = 0.5d0*(ss%R_midp_out(ix)-ss%R_midp_in(ix))
    end do
  else
    do ix = 1, inx
      call eq_glimRZ(ss%rho(ix),zRmin,zRmax,zYmin,zYmax,ier); if(ier.ne.0) go to 1000
      zRminor(ix) = 0.5d0*(zRmax-zRmin)
    end do
  endif

  znebar = 0.0d0
  do ix=1,inx-1
    ! Line averaged density
    znebar = znebar + ss%ns(ix,0)*abs(zrminor(ix+1)-zrminor(ix))
  end do
  ss%nebar = znebar/(zrminor(inx)+1.0d-6)


  !  radiated power...
  call get1p('PRAD',1); if(ier.ne.0) go to 1000
  ss%prad = zprofc

  !  attempt to determine source -- initially assume there is input data
  ss%rad_data_info='TRANSP:input_data (bolometer)'

  zprofc1=0.0; zprofc2=0.0; zprofc3=0.0; zprofc=0.0
  call get1p('PRAD0',1)
  if(ier.eq.0) then
     zprofc1=zprofc
     call get1p('PRAD_ADJ',1)
     if(ier.eq.0) zprofc2=zprofc
     call get1p('PRADC',1)
     if(ier.eq.0) zprofc3=zprofc
     if(matchNZ(ss%prad,zprofc1,ZERO).or.matchNZ(ss%prad,zprofc2,ZERO)) then
        write(lunzer(0),*) ' %trx_gen_state: PRAD matches bolometer data.'
     else if(matchNZ(ss%prad,zprofc3,ZERO)) then
        write(lunzer(0),*) ' %trx_gen_state: PRAD matches theoretical model.'
        ss%rad_data_info='TRANSP:theoretical_estimate'
     else
        write(lunzer(0),*) ' %trx_gen_state: PRAD match not found.'
        ss%rad_data_info='TRANSP:unknown_source'
     endif
  else
     write(lunzer(0),*) ' %trx_gen_state: old run, assume PRAD from bolometer.'
  endif
  ier=0

  ! theory contributions
  call get1p('PRAD_CY',1)
  if(ier.ne.0) then
     write(lunzer(0),*) ' %trx_gen_state: only total PRAD is available.'
     ier = 0
  else
     ss%prad_cy = zprofc
     call get1p('PRAD_LI',1); if(ier.ne.0) go to 1000
     ss%prad_li = zprofc
     call get1p('PRAD_BR',1); if(ier.ne.0) go to 1000
     ss%prad_br = zprofc
  endif

  !----------------------------------------------------
  !  transport profiles -- powers, W/zone, ang. momentum Nt*m/zone
  !                        particles, (#/sec)/zone
  !  sign convention: positive denotes loss

  call get1p('PCNDE',1)
  ss%pe_trans = zprofc
  call get1p('PCNVE',1)
  ss%pe_trans = ss%pe_trans + zprofc

  call get1p('PCOND',1)
  ss%pi_trans = zprofc
  call get1p('PCONV',1)
  ss%pi_trans = ss%pi_trans + zprofc

  call rpexist_profile('MVISC',iexist)
  if(iexist) then
     call get1p('MVISC',1)
     ss%tq_trans = zprofc
     call get1p('MCONV',1)
     ss%tq_trans = ss%tq_trans + zprofc
  else
     ss%tq_trans = ZERO
  endif

  allocate(znxsum(inx-1)); znxsum = ZERO

  izth=0
  iath=0
  jth=0
  do i=1,n_species
     ith = -1
     profname=' '
     if(itype(i).eq.ps_electron) then
        ith=0
        profname='DIVFE'

     else if(itype(i).eq.ps_therm_ion) then
        jth=jth + 1
        ith=jth

        izth=Zcharg(i)+0.1d0  ! integer charge
        iath=Amass(i)+0.1d0   ! integer AMU

        if(izth.eq.1) then
           if(iath.eq.1) then
              profname='DIVHT'
           else if(iath.eq.2) then
              profname='DIVFD'
           else if(iath.eq.3) then
              profname='DIVFT'
           endif
        else if(izth.eq.2) then
           if(iath.eq.3) then
              profname='DIVHE3T'
           else if(iath.eq.4) then
              profname='DIVHE4T'
           endif
        else if(izth.eq.3) then
           profname='DIVLITHT'
        endif

     else if((itype(i).eq.ps_impurity).or.(itype(i).eq.ps_tokamakium)) then
        jth=jth + 1
        znxsum = znxsum + ss%ns(:,jth)
     endif

     if(ith.ne.-1) then
        if(profname.eq.' ') then
           write(lunzer(0),*) ' ?trx_gen_state: transport profile name error:'
           write(lunzer(0),*) '  species Z, A = ',izth,iath
        else
           call get1p(profname,1)
           ss%sn_trans(:,ith) = zprofc
        endif
     endif
  enddo

  ! divide the impurity transport by density

  call get1p('DFIMP',1)

  jth=0
  do i=1,n_species
     if(itype(i).eq.ps_therm_ion) then
        jth=jth+1
     else if((itype(i).eq.ps_impurity).or.(itype(i).eq.ps_tokamakium)) then
        jth=jth+1
        do j=1,inx-1
           if(znxsum(j).gt.ZERO) then
              znloc=ss%ns(j,jth)
              ss%sn_trans(j,jth) = zprofc(j)*znloc/znxsum(j)
           else
              ss%sn_trans(j,jth) = ZERO
           endif
        enddo
     endif
  enddo

  deallocate(znxsum)

  !----------------------------------------------------
  !  gas flow and recycling sources

  !  try to fetch T0recyc, etc. (available in newer runs); if not available
  !  set to zero.

  if(igeq) then
     zRedge = 0.5d0*(ss%R_midp_in(inx)+ss%R_midp_out(inx))
  else
     call eq_glimRZ(ONE,zRmin,zRmax,zYmin,zYmax,ier); if(ier.ne.0) go to 1000
     zRedge = 0.5d0*(zRmin+zRmax)
  endif

  call trx_scal('T0RECYC',zunits,t0recyc,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' %trx_gen_state: edge data {T0recyc, T0gasfl, OM0recyc, OM0gasfl}'
     write(lunzer(0),*) '  not found (old TRANSP run): zeroed for now; namelist will be examined.'
     t0recyc = ZERO
     t0gasfl = ZERO
     om0recyc = ZERO
     om0gasfl = ZERO
  else
     write(lunzer(0),'(1x,a,1pe12.5,a)') &
          ' %trx_gen_state: edge data available, T0recyc = ',t0recyc,' KeV.'
     call trx_scal('T0GASFL',zunits,t0gasfl,ier)
     write(lunzer(0),'(1x,a,1pe12.5,a)') &
          '  T0gasfl = ',T0gasfl,' KeV.'
     call trx_scal('OM0RECYC',zunits,om0recyc,ier)
     write(lunzer(0),'(1x,a,1pe12.5,a)') &
          '  OM0recyc = ',OM0recyc,' rad/sec.'
     call trx_scal('OM0GASFL',zunits,om0gasfl,ier)
     write(lunzer(0),'(1x,a,1pe12.5,a)') &
          '  OM0gasfl = ',OM0gasfl,' rad/sec.'
  endif

  if(ss%nspec_gas.ge.jig0) then
     allocate(sgas(jig0),srcy(jig0))
     call trx_spec_edge(jig0,idum,sgas,srcy,ier)

     in0=ss%nspec_gas

     sgrcy_sum=ZERO
     do jth=1,in0
        ss%sc0(jth)=srcy(jth)
        ss%e0_av(jth)=1.5d0*T0recyc
        ss%vphi0_av(jth)=om0recyc*zRedge

        ss%sc0(jth+in0)=sgas(jth)
        ss%e0_av(jth+in0)=1.5d0*T0gasfl
        ss%vphi0_av(jth+in0)=om0gasfl*zRedge

        sgrcy_sum = sgrcy_sum + sgas(jth) + srcy(jth)
     enddo
     if(sgrcy_sum.lt.ONE) sgrcy_sum=ONE

     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state: trx_spec_edge error.'
        go to 1000
     endif
  else
     write(lunzer(0),*) ' %trx_gen_state: ss%nspec_gas = ',ss%nspec_gas,'; #thermal=',jig0
     write(lunzer(0),*) '  gas flow and recycling not set (OK).'
  endif

  !----------------------------------------------------
  !  surface voltage & loop voltage profile

  call trx_scal('VSURC',zunits,zvsur,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state: VSURC not found.'
     go to 1000
  endif
  ss%vsur = zvsur

  call trx_scal('VSUR',zunits,zvsur,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state: VSURC not found.'
     go to 1000
  endif

  if(abs(zvsur-ss%vsur).LT.0.001d0) then
     ss%vsur_data_info = 'TRANSP:Surface_voltage:matches_input_data'
  else
     ss%vsur_data_info = 'TRANSP:Surface_voltage:predicted'
  endif

  call get1b('V')
  ss%v_loop = zprof

  !  related...

  ss%tf_data_info = 'TRANSP:vacuum_toroidal_field_is_input'

  !----------------------------------------------------
  !  Now, look for data for items to be read from the TRANSP namelist...

  call trread_lun(ilun_tf)
  call tr_getnl_ready(iflag)
  if(.not.iflag) then
     call tr_getnl_text(ilun_tf,iwarn)
     if(iwarn.eq.0) iflag=.TRUE.
  endif

  if(iflag) then
     zftmp0 = trim(run_label)//'TR.DAT'
     do ic=1,len(trim(zftmp0))
        if(zftmp0(ic:ic).eq.' ') zftmp0(ic:ic)='_'
     enddo

     call tmpfile_d(zftmp0,zftmp,iltmp)

     open(unit=ilun_tf,file=zftmp,status='UNKNOWN',iostat=iwarn)
     if(iwarn.ne.0) iflag=.FALSE.
  endif

  if(iflag) then
     call tr_putnl_write(ilun_tf,iwarn)
     close(unit=ilun_tf)
     if(iwarn.ne.0) iflag=.FALSE.
  endif

  if(splitn_try) then
     if(iflag) then
        !  OK: namelist in tmp file: load SPLITN
        call splitn_read(zftmp,iwarn)
        if(iwarn.ne.0) then
           iflag=.FALSE.
        else
           ! cleanup
           open(unit=ilun_tf,file=zftmp,status='OLD',iostat=iwarn)
           if(iwarn.eq.0) then
              close(unit=ilun_tf,status='delete')
           endif
        endif
     endif

     if(iflag) then
        write(lunzer(0),*) ' %trx_gen_state: TRANSP namelist available; splitn loaded.'
        splitn_avail=.TRUE.
        call splitn_iget('LPredictive_Mode',1,ival,ier)
        if(ier.ne.0) then
           write(lunzer(0),*) ' ?trx_gen_state: LPredictive_Mode access error.'
           go to 1000
        endif

        if((ival(1).eq.0).or.(ival(1).eq.10)) then
           ptransp_flag=.FALSE.
        else
           ptransp_flag=.TRUE.
           write(lunzer(0),*) ' ------------------------------- '
           write(lunzer(0),*) ' *** trx_gen_state: PTRANSP mode '
           write(lunzer(0),*) ' *** trx_gen_state: metadata extraction not complete! '
           write(lunzer(0),*) ' ------------------------------- '
        endif

     else
        write(lunzer(0),*) ' %trx_gen_state: no access to TRANSP namelist data.'
        write(lunzer(0),*) '  namelist quantities omitted.'
        splitn_avail=.FALSE.
     endif
  else
     write(lunzer(0),*) ' %trx_gen_state: caller requests no access to splitn.'
     splitn_avail=.FALSE.
  endif

  if(max(iflag_LH,iflag_ECH).gt.0) then
     if(.not.splitn_avail) then
        write(lunzer(0),*) ' ?trx_gen_state: splitn needed, not accessible.'
        write(lunzer(0),*) '  iflag_ECH or iflag_LH is set.'
        ier=1
        go to 1000
     endif
  endif

  !----------------------------------
  !  metadata... coding for PTRANSP not finished...

  ss%ts_is_input = 0
  ss%ns_is_input = 0
  ss%vtor_is_input = 0
  ss%vpol_is_input = 0

  if(ptransp_flag) then
     ss%ts_data_info = '(PTRANSP)'
     ss%ns_data_info = '(PTRANSP)'
     ss%zeff_data_info = '(PTRANSP)'
     ss%vtor_data_info = '(PTRANSP)'
     ss%vpol_data_info = '(PTRANSP)'
  else
     ! electron temperature...
     call rpexist_profile('TEPRO',iexist)
     if(.not.iexist) then
        ss%ts_is_input(ps_elec_index) = 1 
        ss%ts_data_info = 'TRANSP:{Te_is_input;'  ! more to be added...
     else
        zprofc1 = ZERO; zprofc2 = ZERO
        call get1('TEPRO')
        if(ier.eq.0) zprofc1 = zprofc
        call get1('TE')
        if(ier.eq.0) zprofc2 = zprofc
        ier=0
        if(matchNZ(zprofc1,zprofc2,ZERO)) then
           ss%ts_is_input(ps_elec_index) = 1 
           ss%ts_data_info = 'TRANSP:{Te_is_input;'  ! more to be added...
        else
           ss%ts_data_info = 'TRANSP:{Te_is_predicted;'  ! more to be added...
        endif
     endif

     ! ion temperature...
     call rpexist_profile('TIPRO',iexist)
     if(.not.iexist) then
        ss%ts_data_info = trim(ss%ts_data_info)//'Ti_is_predicted}'
     else
        zprofc1 = ZERO; zprofc2 = ZERO; zprofc3 = ZERO; zprofc4 = ZERO
        call get1('TIPRO')
        if(ier.eq.0) zprofc1 = zprofc
        call get1('TIAV')
        if(ier.eq.0) zprofc2 = zprofc
        call get1('TMJ')
        if(ier.eq.0) zprofc3 = zprofc
        call get1('TX')
        if(ier.eq.0) zprofc4 = zprofc
        ier=0

        if(matchNZ(zprofc1,zprofc2,ZERO)) then
           ss%ts_data_info = trim(ss%ts_data_info)//'Ti(average)_is_input}'
           ss%ts_is_input(imj1:ix2) = 1 
        else if(matchNZ(zprofc1,zprofc3,ZERO)) then
           ss%ts_data_info = trim(ss%ts_data_info)//'Ti(majority)_is_input}'
           ss%ts_is_input(imj1:imj2) = 1 
        else if(matchNZ(zprofc1,zprofc4,ZERO)) then
           ss%ts_data_info = trim(ss%ts_data_info)//'Ti(impurity)_is_input}'
           ss%ts_is_input(ix1:ix2) = 1 
        else
           ss%ts_data_info = trim(ss%ts_data_info)//'Ti_is_predicted}'
        endif
     endif

     !  electron density (always input in traditional TRANSP)
     ss%ns_data_info = 'TRANSP:{ne_is_input'  ! more to be added...
     ss%ns_is_input(0) = 1

     call set_zeffdens_info
     if(ier.ne.0) go to 1000

     call set_vprof_info
     if(ier.ne.0) go to 1000
  endif

  !----------------
  !  metadata: plasma current

  call trx_scal('PCUR',zunits,zpcur,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state: PCUR not found.'
     go to 1000
  endif

  call trx_scal('PCURC',zunits,zpcurc,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state: PCUR not found.'
     go to 1000
  endif

  zpcur_diff = abs(zpcur-zpcurc)/max(ONE,max(abs(zpcur),abs(zpcurc)))

  if(.not.splitn_avail) then
     if(zpcur_diff.lt.1.0d-4) then
        ss%cur_data_info= 'TRANSP:total_current_matches;no_profile_information'
     else
        ss%cur_data_info= 'TRANSP:total_current:no_match'
     endif
  else
     call cur_info(zpcur_diff)
     if(ier.ne.0) go to 1000
  endif
  
  !----------------------------------------------------
  !  Now, look for data for items to be read from the TRANSP input data
  !  (time dependent diagnostic shot data -- if available...)

  if(trdatbuf_try) then

     call trx_trdatbuf_connect(iwarn)

     if(d_data_avail) then
        write(lunzer(0),*) ' %trx_gen_state: TRANSP input data available; trdatbuf loaded.'
        trdatbuf_avail=.TRUE.
     else
        write(lunzer(0),*) ' %trx_gen_state: no access to TRANSP time dependent input data.'
        write(lunzer(0),*) '   trdatbuf quantities omitted.'
        trdatbuf_avail=.FALSE.
     endif

  else
     write(lunzer(0),*) ' %trx_gen_state: caller requests no access to trdatbuf.'
     trdatbuf_avail=.FALSE.
  endif

  !----------------------------------------------------
  !  mod DMC: instead of using namelist, set stop time from last 
  !  output time point.  Use as default for start time as well but
  !  override with namelist TINIT value if available
  !----------------------------------------------------

  ss%tinit = tmin
  ss%tfinal = tmax
  
  if(splitn_avail) then
     !  misc. namelist items...

     ztemp(1) = ss%tinit
     call splitn_dget('TINIT',1,ztemp,ier)
     ss%tinit = ztemp(1)
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state_geq: TINIT fetch failed.'
        go to 1000
     endif

     call splitn_iget('NSHOT',1,itemp,ier)
     ss%shot_number = itemp(1)
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state_geq: NSHOT fetch failed.'
        go to 1000
     endif

     ! data on neutrals, for older runs...

     if(t0recyc.eq.ZERO) then
        ztemp(1) = t0recyc
        call splitn_dget('TIEDGE',1,ztemp(1),ier)
        t0recyc = ztemp(1)
        if(ier.ne.0) then
           write(lunzer(0),*) ' ?trx_gen_state_geq: TIEDGE fetch failed.'
           go to 1000
        endif
        t0recyc = 0.001d0*t0recyc  ! -> KeV
        if(t0recyc.lt.1.0d-3) t0recyc = 1.0d-3
        write(lunzer(0),'(1x,a,1pe12.5,a)') &
             ' %trx_gen_state: T0recyc = ',t0recyc,' KeV set from TIEDGE in namelist.'

        in0=ss%nspec_gas
        call splitn_getsize('E0IN',isize,ier)
        if(ier.ne.0) then
           write(lunzer(0),*) ' ?trx_gen_state_geq: E0in size fetch failed.'
           go to 1000
        endif
        allocate(e0in(isize))
        call splitn_dget('E0IN',isize,e0in,ier)
        if(ier.ne.0) then
           write(lunzer(0),*) ' ?trx_gen_state_geq: E0in fetch failed.'
           go to 1000
        endif

        t0gasfl=ZERO
        do i=1,in0
           t0gasfl=t0gasfl + 0.001d0*e0in(3*i)
        enddo
        t0gasfl=t0gasfl/in0
        if(t0gasfl.lt.1.0d-3) t0gasfl=1.0d-3
        write(lunzer(0),'(1x,a,1pe12.5,a)') &
             ' %trx_gen_state: T0gasfl = ',t0gasfl,' KeV set from E0IN in namelist.'
        do jth=1,in0
           ss%e0_av(jth)=1.5d0*T0recyc
           ss%e0_av(jth+in0)=1.5d0*T0gasfl
        enddo

        deallocate(e0in)
     endif
  endif

  !----------------------------------------------------
  !  the following items are added only if BOTH the namelist and the exp. data
  !  are available

  iflag = splitn_avail.AND.trdatbuf_avail

  if(iflag) then
     !  some metadata labels...

     ss%eq_data_info='TRANSP:bdy_input_data'
     ss%eq_code_info='TRANSP:legacy'

     call splitn_iget('LEVGEO',1,itemp,ier)
     levgeo = itemp(1)
     if(ier.ne.0) go to 1000

     if(levgeo.eq.8) then
        ss%eq_code_info='TRANSP:input_data(scrunch2)'
        ss%eq_data_info='TRANSP:input_data(entire_equilibrium)'
     else if(levgeo.eq.11) then
        ss%eq_code_info='TRANSP:TEQ'
     else if(levgeo.eq.12) then
        ss%eq_code_info='TRANSP:ISOLVER'
     else if(levgeo.gt.12) then
        ss%eq_code_info='TRANSP:unknown'
     endif
     write(ss%eq_code_info(40:),'("LEVGEO=",I2)') levgeo

     call splitn_iget('NSOMOD',1,itemp,ier)
     nsomod = itemp(1)
     if(ier.ne.0) go to 1000

     ss%gas_code_info='TRANSP:unknown'
     if(nsomod.eq.1) then
        ss%gas_code_info='TRANSP:FRANTIC'
     endif

     ss%rad_code_info = 'TRANSP:legacy'

     call splitn_iget('nmdifb',1,nmdifb,ier)
     if(ier.ne.0) go to 1000

     if(nmdifb(1).gt.0) then
        ss%anom_code_info = 'TRANSP:unknown'
        if(nmdifb(1).le.3) then
           ss%anom_code_info = 'TRANSP:1d_input_data:{D,v}[x,t]'
        else if(nmdifb(1).eq.4) then
           ss%anom_code_info = 'TRANSP:input_data:{multi-D[E,x,t]}'
        endif
     else  if(nmdifb(1).eq.-3) then
        ss%anom_code_info = 'TRANSP:1d_input_data:{v}[x,t]'
    endif

  endif

  intrace = 0

  do ibc=1,max_zmbuf
     imap0(ibc)=ibc
  enddo

  if(iflag) then

     !  #beams, #ICRF antennas, #ECRF antennas, #LH antennas:
     !  NOTE: something other than NLTORAY may need to be tested for EC
     !        (eventually -- but then retain NLTORAY test for older runs).

     iwarn=0
     call getnum_heater('NB','NBEAM','NLBEAM',ss%nbeam,inb_trdat,iwarn)
     call getnum_heater('RF','NICHA','NLICRF',ss%nicrf_src,irf_trdat,iwarn)
     call getnum_heater('LH','NANTLH','NLLH',ss%nlhrf_src,ilh_trdat,iwarn)

     call getnum_heater('EC','NANTECH','NLTORAY',ss%necrf_src,iec_trdat,iwarn)
     if(ss%necrf_src.eq.0) then
        ! try this also...
        call getnum_heater('EC','NANTECH','NLGEN_ECH',ss%necrf_src, &
             iec_trdat,iwarn)
     endif
     if(ss%necrf_src.eq.0) then
        ! and this ...
        call getnum_heater('EC','NANTECH','NLTORBEAM',ss%necrf_src, &
             iec_trdat,iwarn)
     endif

     if(iwarn.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state: inconsistency in heating sources.'
        write(lunzer(0),*) '  state contents truncated.'

        iflag = .FALSE.
        ss%nbeam  =0
        ss%nicrf_src =0
        ss%necrf_src =0
        ss%nlhrf_src =0
     endif
  endif

  if(.not.allocated(imap)) then
     allocate(imap(size(imap0)))
     allocate(imapi(size(imap0)))
     imap=imap0
     imapi=imap0
  endif

  if(iflag) then
        
     !  count types of auxilliary heating
     !  also, set some PS "metadata" labels

     iaux = 0

     !  set rho grid sizes to match TRANSP

     if(ss%nbeam.gt.0) then
        call mk_rho_nbi
        iaux = iaux + 1
     endif

     if(ss%nicrf_src.gt.0) then
        call mk_rho_icrf
        iaux = iaux + 1

        ss%ic_code_info = 'TRANSP:unknown'
        call splitn_iget('NICRF',1,itemp,ier)
        nicrf = itemp(1)
        if(ier.ne.0) go to 1000
        
        if(nicrf.eq.6) then
           ier=1
           go to 1000
        else if(nicrf.gt.6) then
           ss%ic_code_info = 'TRANSP:TORICv5+'
        endif

        ss%ic_data_info = 'TRANSP:input_data'
     endif

     if(ss%necrf_src.gt.0) then
        call mk_rho_ecrf
        iaux = iaux + 1

        ss%ec_code_info = 'TRANSP:unknown'

        call splitn_lget('NLTORAY',1,ltemp,ier)
        nltoray = ltemp(1)
        if(ier.ne.0) go to 1000

        if(nltoray) then
           ss%ec_code_info = 'TRANSP:TORAY(GA)'
        endif

        call splitn_lget('NLGEN_ECH',1,ltemp,ier)
        nlgen_ech = ltemp(1)
        if(ier.ne.0) go to 1000

        if(nlgen_ech) then
           ss%ec_code_info = 'TRANSP:GENRAY(ech)'
        endif

        call splitn_lget('NLTORBEAM',1,ltemp,ier)
        nltorbeam = ltemp(1)
        if(ier.ne.0) go to 1000

        if(nltorbeam) then
           ss%ec_code_info = 'TRANSP:TORBEAM(ech)'
        endif

        ss%ec_data_info = 'TRANSP:input_data'
     endif

     if(ss%nlhrf_src.gt.0) then
        call mk_rho_lhrf
        iaux = iaux + 1

        ss%lh_code_info = 'TRANSP:LSC'
        ss%lh_data_info = 'TRANSP:input_data'
     endif

     if(iaux.gt.0) then
        !  incremental allocation for auxilliary heating profiles...
        !    (possibly done in mk_rho* calls also): ps_alloc_plasma_state

        if(ier.eq.0) call ps_alloc_plasma_state(ier, state=ss)
        if(ier.ne.0) then
           write(lunzer(0),*) ' ?trx_gen_state_geq: incremental state allocation failure.'
           write(lunzer(0),*) ' ?trx_gen_state_geq: ps_alloc_plasma_state status: ',&
                ier
           go to 1000
        endif

        ! get names of injectors & antennas...
        ! also set their output grids (just use the TRANSP grid).
        
        ier=0
        call splitn_blabel('bnames', 'B_', ss%nbeam, ss%nbi_src_name, &
             imapi, iertmp)
        ier=max(ier,iertmp)
        call splitn_blabel('icnames', 'IC', ss%nicrf_src, ss%icrf_src_name, &
             imap0, iertmp)
        ier=max(ier,iertmp)
        call splitn_blabel('ecnames', 'EC', ss%necrf_src, ss%ecrf_src_name, &
             imap0, iertmp)
        ier=max(ier,iertmp)
        call splitn_blabel('lhnames', 'LH', ss%nlhrf_src, ss%lhrf_src_name, &
             imap0, iertmp)
        ier=max(ier,iertmp)

        if(ier.ne.0) go to 1000

        !  set rho grids to match TRANSP
        !  fill in state scalar data

        call tdb_pwrget_init(zpwr)
        zpwr%ztime1 = ss%t0
        zpwr%ztime2 = ss%t1

        if(ss%nbeam.gt.0) then
           
           !  neutral beam data: species of each beam

           inum=ss%nbeam
           inumtot = inum + intrace

           allocate(abeama(inumtot),xzbeama(inumtot))
           call splitn_dvec('abeama',inumtot,abeama,iertmp)
           ier=max(ier,iertmp)
           call splitn_dvec('xzbeama',inumtot,xzbeama,iertmp)
           ier=max(ier,iertmp)

           if(ier.ne.0) go to 1000

           do ib=1,inumtot
              ia=abeama(ib)+0.5d0
              iz=xzbeama(ib)+0.5d0
              if(abs(ntrace(ib)).eq.0) then
                 ib0=imap(ib)
                 do i=1,n_bi
                    if((ia.eq.ianbi(i)).and.(iz.eq.iznbi(i))) then
                       if(iz.eq.1) then
                          if(ia.eq.1) then
                             ss%nbion(ib0)='H'
                          else if(ia.eq.2) then
                             ss%nbion(ib0)='D'
                          else if(ia.eq.3) then
                             ss%nbion(ib0)='T'
                          endif
                       else if(iz.eq.2) then
                          if(ia.eq.3) then
                             ss%nbion(ib0)='HE3'
                          else if(ia.eq.4) then
                             ss%nbion(ib0)='HE4'
                          endif
                       else if(iz.eq.10) then
                          ss%nbion(ib0)='Ne'   ! Neon
                       else if(iz.eq.18) then
                          ss%nbion(ib0)='Ar'   ! Argon
                       else if(iz.eq.36) then
                          ss%nbion(ib0)='Kr'   ! Krypton
                       else if(iz.eq.54) then
                          ss%nbion(ib0)='Xe'   ! Xenon
                       endif
                       exit
                    endif
                 enddo
              else
                 ib0=abs(ntrace(ib))
                 ib0=imap(ib0)
                 do i=1,n_bi
                    if((ia.eq.ianbi(i)).and.(iz.eq.iznbi(i))) then
                       if(iz.eq.1) then
                          if(ia.eq.1) then
                             ss%nbion_trace(ib0)='H'
                          else if(ia.eq.2) then
                             ss%nbion_trace(ib0)='D'
                          else if(ia.eq.3) then
                             ss%nbion_trace(ib0)='T'
                          endif
                       else if(iz.eq.2) then
                          if(ia.eq.3) then
                             ss%nbion_trace(ib0)='HE3'
                          else if(ia.eq.4) then
                             ss%nbion_trace(ib0)='HE4'
                          endif
                       else if(iz.eq.10) then
                          ss%nbion_trace(ib0)='Ne'   ! Neon
                       else if(iz.eq.18) then
                          ss%nbion_trace(ib0)='Ar'   ! Argon
                       else if(iz.eq.36) then
                          ss%nbion_trace(ib0)='Kr'   ! Krypton
                       else if(iz.eq.54) then
                          ss%nbion_trace(ib0)='Xe'   ! Xenon
                       endif
                       exit
                    endif
                 enddo
              endif
           enddo

           !  look for trdat data; if not available will use rplot or
           !  namelist data...

           iflg_pwr  = tdb_logchk_nbi(d,'PWR',idum)
           iflg_vlt  = tdb_logchk_nbi(d,'VLT',idum)
           iflg_full = tdb_logchk_nbi(d,'FUL',idum)
           iflg_half = tdb_logchk_nbi(d,'HLF',idum)

           if(iflg_vlt.or.iflg_full.or.iflg_half) then
              if(.not.iflg_pwr) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: trdatbuf beam power: missing.'
                 write(lunzer(0),*) &
                      '  but voltage and/or energy fractions vs. t exist.'
                 ier=1
                 go to 1000
              endif
           endif

           ! TRDAT data is for the main beams only, no trace beams;
           ! so the TRDAT indexing of beams matches that of the Plasma State
           !   (ib=imap(ib) not needed).

           if(iflg_pwr) then
              !  beam power & voltage.  need to pre-fetch on/off times...
              call get_onoff('NB',inum,ier)  ! nbi on off times
              if(ier.ne.0) go to 1000

              zpwr%nbeam=inum
           
              zpwr%item='PWR'     ! beam powers
              zpwr%pweight=.FALSE.

              call tdb_pwrdata_avg(d,zpwr,ier)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: tdb_pwrdata_avg "PWR" error!'
                 go to 1000
              endif
              ss%power_nbi(1:inum) = zpwr%zparam(1:inum)
           else
              ! get powers from RPLOT
              do ii=1,inumtot
                 ib=imap(ii)
                 if(ib.gt.0) then
                    zname=' '
                    write(zname,'("PINJ",i2.2)') ii
                    call trx_scal(trim(zname),zunits,zitem,ier)
                    if(ier.ne.0) then
                       write(lunzer(0),*) &
                            ' ?trx_gen_state: could not read in rplot: ', &
                            trim(zname)
                       go to 1000
                    endif
                    ss%power_nbi(ib)=zitem
                 endif
              enddo
           endif

           ! correct for trace beam powers if necessary
           do ib=1,inum
              if(ftraci(ib).gt.ZERO) then
                 ss%power_nbi_trace(ib) = ftraci(ib)*ss%power_nbi(ib)
                 ss%power_nbi(ib) = (ONE-ftraci(ib))*ss%power_nbi(ib)
              endif
           enddo

           if(iflg_vlt) then
              zpwr%item='VLT'     ! beam voltages
              zpwr%pweight=.TRUE. ! power weighted, all NBI data share timebase

              call tdb_pwrdata_avg(d,zpwr,ier)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: tdb_pwrdata_avg "VLT" error!'
                 go to 1000
              endif
              ss%kvolt_nbi(1:inum) = zpwr%zparam(1:inum)*0.001d0  ! -> keV
           else
              ! get powers from RPLOT
              do ii=1,inumtot
                 ib=imap(ii)
                 if(ib.gt.0) then
                    zname=' '
                    write(zname,'("EINJ",i2.2,"_E1")') ii
                    call trx_scal(trim(zname),zunits,zitem,ier)
                    if(ier.ne.0) then
                       write(lunzer(0),*) &
                            ' ?trx_gen_state: could not read in rplot: ', &
                            trim(zname)
                       go to 1000
                    endif
                    ss%kvolt_nbi(ib)=zitem    ! in KeV already
                 endif
              enddo
           endif

           ! put energy fractions in local arrays for the moment...
           allocate(ffulla(inum),fhalfa(inum))

           if(iflg_full) then
              zpwr%item='FUL'     ! full energy fraction
              zpwr%pweight=.TRUE. ! power weighted, all NBI data share timebase

              call tdb_pwrdata_avg(d,zpwr,ier)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: tdb_pwrdata_avg "FUL" error!'
                 go to 1000
              endif
              ffulla(1:inum) = zpwr%zparam(1:inum)
           else
              zconv=ONE
              call nbget_alg1('FFULL','FFULLA',ffulla)
              if(ier.ne.0) go to 1000
           endif

           if(iflg_half) then
              zpwr%item='HLF'     ! half energy fraction
              zpwr%pweight=.TRUE. ! power weighted, all NBI data share timebase

              call tdb_pwrdata_avg(d,zpwr,ier)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: tdb_pwrdata_avg "HLF" error!'
                 go to 1000
              endif
              fhalfa(1:inum) = zpwr%zparam(1:inum)
           else
              zconv=ONE
              call nbget_alg1('FHALF','FHALFA',fhalfa)
              if(ier.ne.0) go to 1000
           endif

           !  radial grid for profile outputs

           ss%rho_nbi = ss%rho

           !  beam elements of machine description

           zconv=0.01d0  ! cm -> m

           !  aperture dimensions (if nbap_shape="Circle" only ap_halfwidth
           !  applies)

           call nbget_alg2('XZEDGE','XZPEDGA',ss%ap_halfheight)
           if(ier.ne.0) go to 1000

           call nbget_alg2('REDGE','RAPEDGA',ss%ap_halfwidth)
           if(ier.ne.0) go to 1000

           call nbget_alg0('XZPEDG2',ss%ap2_halfheight)
           call nbget_alg0('RAPEDG2',ss%ap2_halfwidth)

           call nbget_alg0('XRAPOFFA',ss%ap_horiz_offset)
           call nbget_alg0('XZAPOFFA',ss%ap_vert_offset)

           call nbget_alg0('XRAPOFF2',ss%ap2_horiz_offset)
           call nbget_alg0('XZAPOFF2',ss%ap2_vert_offset)

           !  ion source dimensions

           call nbget_alg1('BMWIDR','BMWIDRA',ss%b_halfwidth)
           if(ier.ne.0) go to 1000

           call nbget_alg1('BMWIDZ','BMWIDZA',ss%b_halfheight)
           if(ier.ne.0) go to 1000

           !  focal lengths (both used if ion source is rectangular)_

           call nbget_alg1('FOCLR','FOCLRA',ss%b_hfocal_length)
           if(ier.ne.0) go to 1000

           call nbget_alg1('FOCLZ','FOCLZA',ss%b_vfocal_length)
           if(ier.ne.0) go to 1000

           !  beam divergences (both used if ion source is rectangular)_

           zconv=360.0d0*sqrt(2.0d0)/twopi ! radians -> (1/e half angle, degrees)

           call nbget_alg1('DIVR','DIVRA',ss%b_hdivergence)
           if(ier.ne.0) go to 1000

           call nbget_alg1('DIVZ','DIVZA',ss%b_vdivergence)
           if(ier.ne.0) go to 1000

           !  distances, source to aperture & source to tangency radius

           zconv=0.01d0  ! cm -> m

           call nbget_alg2('XLNBAP','XLBAPA',ss%lbscap)
           if(ier.ne.0) go to 1000

           call nbget_alg0('XLBAPA2',ss%lbscap2)

           call nbget_alg2('XLNBTN','XLBTNA',ss%lbsctan)
           if(ier.ne.0) go to 1000

           !  signed tangency radius:
           !    positive => injection to impart momentum in CCW direction
           !      when torus is viewed from above; o.w. negative

           call nbget_alg2('RTCEN','RTCENA',ss%srtcen)
           if(ier.ne.0) go to 1000
           
           !  set ss%srtcen sign:
           !        ss%kccw_Jphi = 1 (Jphi flows CCW viewed from above)
           !            NLCO(ib)=.TRUE. => ss%srtcen(ib) > 0
           !            NLCO(ib)=.FALSE.=> ss%srtcen(ib) < 0
           !        ss%kccw_Jphi = -1 (Jphi flows CW viewed from above)
           !            NLCO(ib)=.TRUE. => ss%srtcen(ib) < 0
           !            NLCO(ib)=.FALSE.=> ss%srtcen(ib) > 0

           call splitn_getsize('NLCO',isize,ier)
           if(ier.ne.0) then
              write(lunzer(0),*) ' ?trx_gen_state_geq: NLCO size fetch failed.'
              go to 1000
           endif
           allocate(ilco(isize))
           call splitn_lget('NLCO',isize,ilco,ier)
           if(ier.ne.0) then
              write(lunzer(0),*) ' ?trx_gen_state: failed to read NLCO namelist data.'
              deallocate(ilco)
              go to 1000
           endif

           do ib=1,inumtot
              ib0=imap(ib)
              if(ib0.gt.0) then
                 isign=ss%kccw_Jphi
                 if(.not.ilco(ib)) isign = -isign

                 ss%srtcen(ib0) = isign * ss%srtcen(ib0)
              endif
           enddo

           deallocate(ilco)

           ! toroidal position of beam source

           zconv = ONE   ! degrees in TRANSP namelist & in plasma state

           call nbget_alg0('XBZETA',ss%Phibsc)
           allocate(iwk(isize))

           ! shape of aperture(s)

           call splitn_iget('NBAPSHA',isize,iwk,ier)
           if(ier.ne.0) then
              write(lunzer(0),*) ' ?trx_gen_state: failed to read NBAPSHA namelist data.'
              deallocate(iwk)
              go to 1000
           endif

           do ib=1,inumtot
              ib0=imap(ib)
              if(ib0.gt.0) then
                 if(iwk(ib).eq.1) then
                    ss%nbap_shape(ib0)='Rectangle'
                 else
                    ss%nbap_shape(ib0)='Circle'
                 endif
              endif
           enddo

           call splitn_iget('NBAPSH2',isize,iwk,ier)
           if(ier.ne.0) then
              write(lunzer(0),*) ' ?trx_gen_state: failed to read NBAPSHA namelist data.'
              deallocate(iwk)
              go to 1000
           endif

           do ib=1,inum
              ib0=imap(ib)
              if(ib0.gt.0) then
                 if(iwk(ib).le.0) then
                    ss%nbap2_shape(ib0)='None'
                 else
                    if(iwk(ib).eq.1) then
                       ss%nbap2_shape(ib0)='Rectangle'
                    else
                       ss%nbap2_shape(ib0)='Circle'
                    endif
                 endif
              endif
           enddo

           ! shape of beam source

           call splitn_iget('NBSHAPA',isize,iwk,ier)
           if(ier.ne.0) then
              write(lunzer(0),*) ' ?trx_gen_state: failed to read NBSHAPA namelist data.'
              deallocate(iwk)
              go to 1000
           endif

           call splitn_iget('NBSHAP',1,itemp,ier)
           ishap = itemp(1)
           if(ier.ne.0) then
              write(lunzer(0),*) ' ?trx_gen_state: failed to read NBSHAP namelist data.'
              deallocate(iwk)
              go to 1000
           endif

           do ib=1,inum
              ib0=imap(ib)
              if(ib0.gt.0) then
                 if(iwk(ib).eq.0) iwk(ib)=ishap
                 if(iwk(ib).eq.1) then
                    ss%nbshape(ib0)='Rectangle'
                 else
                    ss%nbshape(ib0)='Circle'
                 endif
              endif
           enddo

           ! height of aperture & beam source

           zconv=0.01d0  ! cm -> m

           call nbget_alg0('XYBAPA',ss%zbap)
           if(ier.ne.0) go to 1000

           call nbget_alg0('XYBSCA',ss%zbsc)
           if(ier.ne.0) go to 1000

           !------------------------
           !  compare beamlist to reference description, if available
           if(aux%nbeam .gt. 0) then
              allocate(inb_sublist(inum))
              call trx_nb_sublist(ss,aux,inum,inb_sublist, &
                   lunzer(0),imatch,istat)
              if(istat.ne.0) then
                 write(lunzer(0),*) &
                      ' %trx_gen_state warning, beam sublist order not found.'
              endif
           endif

           !------------------------
           !  energy fractions -- one for each distinct voltage & species
           !  first set beam types: Helium or Standard or Negative_ion
           !  Standard type has full:half:third beam energy fractions
           !  i.e. fractions of the beam current that carry full voltage,
           !  1/2 voltage, and 1/3 voltage beam particles...

           allocate(kvfrac(inum+1,3),ibi(inum,3))
           kvfrac = ZERO
           ibi = 0

           zemmx_std(1)=200.0d0
           zemmx_std(2)=10.0d0

           zemmx_neg(1)=4000.0d0
           zemmx_neg(2)=20.0d0

           do ib=1,inum
              if((ss%nbion(ib).eq.'HE3').or.(ss%nbion(ib).eq.'HE4')) then
                 ss%beam_type(ib)='Helium'
              else if((ss%nbion(ib).eq.'H').or.(ss%nbion(ib).eq.'D').or. &
                   (ss%nbion(ib).eq.'T')) then
                 if(ffulla(ib).ge.(0.9999d0)) then
                    ss%beam_type(ib)='Negative_ion'  ! high energy beams
                    !  with 100% full energy fraction
                 else
                    ss%beam_type(ib)='Standard'
                    !  standard beam with mixed energy fractions...
                    if(ss%nbion(ib).eq.'H') indx=1
                    if(ss%nbion(ib).eq.'D') indx=2
                    if(ss%nbion(ib).eq.'T') indx=3
                    do ic=1,inum
                       if(kvfrac(ic,indx).eq.ZERO) then
                          kvfrac(ic,indx)=ss%kvolt_nbi(ib)
                          ibi(ic,indx)=ib
                          exit
                       else if(kvfrac(ic,indx).eq.ss%kvolt_nbi(ib)) then
                          exit
                       endif
                    enddo
                 endif
              else
                 ss%beam_type(ib)='Impurity'
              endif

              if(ss%power_nbi(ib).gt.ZERO) then
                 iemmx(1)=ss%kvolt_nbi(ib)/10  ! 10 KeV increments...
                 iemmx(1)=max(2,iemmx(1))
                 iemmx(2)=iemmx(1)+1
                 iemmx(1)=iemmx(1)/2
                 iemmx=iemmx*10  ! lower & upper limits on beam energy
                 ! in 10 KeV increments

                 if(ss%beam_type(ib).eq.'Negative_ion') then
                    zemmx=iemmx(1)
                    zemmx_neg(1)=min(zemmx_neg(1),zemmx)
                    zemmx=iemmx(2)
                    zemmx_neg(2)=max(zemmx_neg(2),zemmx)
                 else
                    zemmx=iemmx(1)
                    zemmx_std(1)=min(zemmx_std(1),zemmx)
                    zemmx=iemmx(2)
                    zemmx_std(2)=max(zemmx_std(2),zemmx)
                 endif
              endif

           enddo

           ! MOD DMC Nov 2010: option to substitute energy fraction
           ! tables from pre-read machine description, in "aux" state.

           if((aux%n_D_Einj_Standard.gt.0).OR. &
                (aux%n_H_Einj_Standard.gt.0).OR. &
                (aux%n_T_Einj_Standard.gt.0)) then

              ss%n_H_Einj_Standard = aux%n_H_Einj_Standard
              ss%n_D_Einj_Standard = aux%n_D_Einj_Standard
              ss%n_T_Einj_Standard = aux%n_T_Einj_Standard
 
              call ps_alloc_plasma_state(ier, state=ss)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state_geq: ps_alloc_plasma_state (nbi): ',&
                      ier
                 go to 1000
              endif
             
              if(ss%n_H_Einj_Standard.gt.0) then
                 ss%H_Einj_standard = aux%H_Einj_standard
                 ss%H_FFull_standard = aux%H_FFull_standard
                 ss%H_FHalf_standard = aux%H_FHalf_standard
              endif
             
              if(ss%n_D_Einj_Standard.gt.0) then
                 ss%D_Einj_standard = aux%D_Einj_standard
                 ss%D_FFull_standard = aux%D_FFull_standard
                 ss%D_FHalf_standard = aux%D_FHalf_standard
              endif
             
              if(ss%n_T_Einj_Standard.gt.0) then
                 ss%T_Einj_standard = aux%T_Einj_standard
                 ss%T_FFull_standard = aux%T_FFull_standard
                 ss%T_FHalf_standard = aux%T_FHalf_standard
              endif

           else
              ! no "aux" state data, so...
              ! form sorted list of injection energies for each species
              ! & corresponding fractions; first get list sizes

              do indx=1,3
                 kvminm(indx)=1.0d30
                 ibminm(indx)=0
                 do ic=1,inum+1
                    if(kvfrac(ic,indx).eq.ZERO) exit
                    if(kvfrac(ic,indx).lt.kvminm(indx)) then
                       kvminm(indx)=kvfrac(ic,indx)
                       ibminm(indx)=ibi(ic,indx)
                    endif
                 enddo
                 if(indx.eq.1) ss%n_H_Einj_Standard = ic - 1
                 if(indx.eq.2) ss%n_D_Einj_Standard = ic - 1
                 if(indx.eq.3) ss%n_T_Einj_Standard = ic - 1
              enddo

              call ps_alloc_plasma_state(ier, state=ss)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state_geq: ps_alloc_plasma_state status(nbi): ',&
                      ier
                 go to 1000
              endif

              if(ss%n_H_Einj_Standard.gt.0) then
                 call fracsort(kvminm(1),ibminm(1),kvfrac(:,1),ibi(:,1), &
                      ss%H_Einj_standard, &
                      ss%H_FFull_standard,ss%H_Fhalf_standard)
              endif

              if(ss%n_D_Einj_Standard.gt.0) then
                 call fracsort(kvminm(2),ibminm(2),kvfrac(:,2),ibi(:,2), &
                      ss%D_Einj_standard, &
                      ss%D_FFull_standard,ss%D_Fhalf_standard)
              endif

              if(ss%n_T_Einj_Standard.gt.0) then
                 call fracsort(kvminm(3),ibminm(3),kvfrac(:,3),ibi(:,3), &
                      ss%T_Einj_standard, &
                      ss%T_FFull_standard,ss%T_Fhalf_standard)
              endif
           endif

           do ib=1,inum
              if(ss%beam_type(ib).eq.'Negative_ion') then
                 ss%einj_min(ib)=zemmx_neg(1)
                 ss%einj_max(ib)=zemmx_neg(2)
              else
                 ss%einj_min(ib)=zemmx_std(1)
                 ss%einj_max(ib)=zemmx_std(2)
              endif
           enddo

        endif

        if(ss%nicrf_src.gt.0) then

           !  ICRF: some antenna data & other scalar information

           inum=ss%nicrf_src

           iflg_pwr  = tdb_logchk_special(d,'RFP',idum)
           iflg_frq  = tdb_logchk_special(d,'RFF',idum)

           if(iflg_pwr) then
              !  get RF antenna powers and frequencies; get on/off times first.
              call get_onoff('RF',inum,ier)
              if(ier.ne.0) go to 1000

              zpwr%nbeam=inum
           
              zpwr%item='RFP'     ! RF antenna powers
              zpwr%pweight=.FALSE.

              call tdb_pwrdata_avg(d,zpwr,ier)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: tdb_pwrdata_avg "RFP" error!'
                 go to 1000
              endif
              ss%power_ic(1:inum) = zpwr%zparam(1:inum)

           else
              !  get ICRF powers from RPLOT data
              do ii=1,inum
                 zname=' '
                 if(ii.lt.10) then
                    write(zname,'("PICHA",i1)') ii
                 else
                    write(zname,'("PICHA",i2)') ii
                 endif
                 call trx_scal(trim(zname),zunits,zitem,ier)
                 if(ier.ne.0) then
                    write(lunzer(0),*) &
                         ' ?trx_gen_state: could not read in rplot: ', &
                         trim(zname)
                    go to 1000
                 endif
                 ss%power_ic(ii)=zitem
              enddo
           endif

           if(iflg_frq.and.iflg_pwr) then
              zpwr%item='RFF'       ! RF frequencies
              zpwr%pweight=.FALSE.  ! In general the RF frequencies and powers
              !  are on different timebases, so, no power weighting option.

              call tdb_pwrdata_avg(d,zpwr,ier)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: tdb_pwrdata_avg "RFF" error!'
                 go to 1000
              endif
              ss%freq_ic(1:inum) = zpwr%zparam(1:inum) ! Hz

           else
              !  get ICRF frequencies from RPLOT data
              do ii=1,inum
                 zname=' '
                 if(ii.lt.10) then
                    write(zname,'("FREQA",i1)') ii
                 else
                    write(zname,'("FREQA",i2)') ii
                 endif
                 call trx_scal(trim(zname),zunits,zitem,ier)
                 if(ier.ne.0) then
                    write(lunzer(0),*) &
                         ' ?trx_gen_state: could not read in rplot: ', &
                         trim(zname)
                    go to 1000
                 endif
                 ss%freq_ic(ii)=zitem
              enddo
           endif

           ss%rho_icrf = ss%rho

           call antgeo
           if(ier.ne.0) go to 1000

        endif

        if(ss%necrf_src.gt.0) then

           !  ECRF: some antenna data & other scalar information

           inum=ss%necrf_src

           iflg_pwr  = tdb_logchk_special(d,'ECP',idum)

           if(iflg_pwr) then
              !  get RF antenna powers and frequencies; get on/off times first.
              call get_onoff('EC',inum,ier)
              if(ier.ne.0) go to 1000

              zpwr%nbeam=inum
           
              zpwr%item='ECP'     ! RF antenna powers
              zpwr%pweight=.FALSE.

              call tdb_pwrdata_avg(d,zpwr,ier)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: tdb_pwrdata_avg "RFP" error!'
                 go to 1000
              endif
              ss%power_ec(1:inum) = zpwr%zparam(1:inum)

           else
              !  get ECH powers from RPLOT data
              do ii=1,inum
                 zname=' '
                 if(ii.lt.10) then
                    write(zname,'("PECIN",i1)') ii
                 else
                    write(zname,'("PECIN",i2)') ii
                 endif
                 call trx_scal(trim(zname),zunits,zitem,ier)
                 if(ier.ne.0) then
                    write(lunzer(0),*) &
                         ' ?trx_gen_state: could not read in rplot: ', &
                         trim(zname)
                    go to 1000
                 endif
                 ss%power_ec(ii)=zitem
              enddo
           endif

           ! probe for time dependent EC beam aiming data...
           call rpexist_scalar("THETECH1",ec_aim_data)
           if(ec_aim_data) then
              do ii=1,inum
                 zname=' '
                 if(ii.lt.10) then
                    write(zname,'("THETECH",i1)') ii
                 else
                    write(zname,'("THETECH",i2)') ii
                 endif
                 call trx_scal(trim(zname),zunits,zitem,ier)
                 ss%EC_theta_aim(ii) = zitem       ! degrees

                 zname=' '
                 if(ii.lt.10) then
                    write(zname,'("PHAIECH",i1)') ii
                 else
                    write(zname,'("PHAIECH",i2)') ii
                 endif
                 call trx_scal(trim(zname),zunits,zitem,ier)
                 ss%EC_phi_aim(ii) = zitem         ! degrees
              enddo
           endif

           ztest = ZERO
           ! (no frequency vs. t in trdat data but there might be someday)

           if(ztest.eq.ZERO) then
              ! frequency data appears to be missing from the trdat data;
              ! fall back to namelist

              allocate(zzbuf(inum))
              call splitn_dvec('FREQECH',inum,zzbuf,ier)
              ss%freq_ec(1:inum) = zzbuf(1:inum)         ! Hz
              deallocate(zzbuf)
              if(ier.ne.0) go to 1000
           endif

           !  additional EC related quantities
           allocate(zzbuf(inum))

           call splitn_dvec('XECECH',inum,zzbuf,ier)
           ss%R_EC_launch(1:inum) = zzbuf(1:inum)*1.0d-2 ! cm -> M

           call splitn_dvec('ZECECH',inum,zzbuf,ier)
           ss%Z_EC_launch(1:inum) = zzbuf(1:inum)*1.0d-2 ! cm -> M

           ss%Phi_EC_launch(1:inum) = 0.0 ! (degrees) (no such data in TRANSP)

           call splitn_dvec('BHALFECH',inum,zzbuf,ier)
           ss%EC_Half_Power_Angle(1:inum) = zzbuf(1:inum) ! degrees

           ! Note: unsure of "beam elongation" definition in TORAY_GA
           ! as used by TRANSP: print warning if value other than 1 is
           ! detected

           call splitn_dvec('BSRATECH',inum,zzbuf,ier)
           if(((ONE-minval(zzbuf)).gt.1.0d-4).or. &
                ((maxval(zzbuf)-ONE).gt.1.0d-4)) then
              write(lunzer(0),*) &
                   ' Caution: TRANSP namelist BSRATECH deviates from 1.00:'
              write(lunzer(0),*) ' minval: ',minval(zzbuf), &
                   '  maxval: ',maxval(zzbuf)
              write(lunzer(0),*) &
                   ' Compatibility of TRANSP/GA_TORAY definition with '// &
                   ' Plasma State: should be checked.'
           endif
           ss%EC_Beam_Elongation(1:inum) = zzbuf(1:inum)  ! dimensionless

           if(.NOT.ec_aim_data) then
              !  not found in TRANSP output, reverting to namelist...
              call splitn_dvec('THETECH',inum,zzbuf,ier)
              ss%EC_theta_aim(1:inum) = zzbuf(1:inum)       ! degrees

              call splitn_dvec('PHAIECH',inum,zzbuf,ier)
              ss%EC_phi_aim(1:inum) = zzbuf(1:inum)         ! degrees
           endif

           call splitn_dvec('RFMODECH',inum,zzbuf,ier)
           ss%EC_Omode_fraction(1:inum) = zzbuf(1:inum)   ! dimensionless

           deallocate(zzbuf)

           call mk_rho_ecrf
        endif

        if(ss%nlhrf_src.gt.0) then

           !  LHRF: some antenna data & other scalar information

           inum=ss%nlhrf_src

           iflg_pwr  = tdb_logchk_special(d,'ECP',idum)

           if(inum.gt.1) then
              write(lunzer(0),*) &
                   ' **CAUTION, it is not clear that TRANSP properly handles'
              write(lunzer(0),*) &
                   '   cases with multiple LH antennas: trx_gen_state.f90'
           endif

           if(iflg_pwr) then
              !  get RF antenna powers and frequencies; get on/off times first.
              call get_onoff('LH',inum,ier)
              if(ier.ne.0) go to 1000

              zpwr%nbeam=inum
           
              zpwr%item='LHP'     ! RF antenna powers
              zpwr%pweight=.FALSE.

              call tdb_pwrdata_avg(d,zpwr,ier)
              if(ier.ne.0) then
                 write(lunzer(0),*) &
                      ' ?trx_gen_state: tdb_pwrdata_avg "RFP" error!'
                 go to 1000
              endif
              ss%power_lh(1:inum) = zpwr%zparam(1:inum)

           else
              !  get LH powers from RPLOT data
              do ii=1,inum
                 zname=' '
                 if(ii.lt.10) then
                    write(zname,'("PLHANT",i1)') ii
                 else
                    write(zname,'("PLHANT",i2)') ii
                 endif
                 call trx_scal(trim(zname),zunits,zitem,ier)
                 if(ier.ne.0) then
                    write(lunzer(0),*) &
                         ' ?trx_gen_state: could not read in rplot: ', &
                         trim(zname)
                    if(inum.eq.1) then
                       write(lunzer(0),*) ' #antennas = 1: recovery attempt:'
                       call trx_scal('PLH',zunits,zitem,ier)
                       if(ier.ne.0) then
                          write(lunzer(0),*) ' ?PLH also not available.'
                          go to 1000
                       else
                          write(lunzer(0),*) '  OK: PLH available.'
                       endif
                    else
                       go to 1000
                    endif
                 endif
                 ss%power_lh(ii)=zitem
              enddo
           endif

           ztest = ZERO
           !  Note: LHRF frequency vs. time not yet in TRANSP trdatbuf...

           if(ztest.eq.ZERO) then
              ! frequency data appears to be missing from the trdat data;
              ! fall back to namelist

              allocate(zzbuf(1))
              call splitn_dvec('FGHZLH',inum,zzbuf,ier)
              ! (here deal with code change: FGHZLH was scalar, is now vector)
              do ii=2,inum
                 if(zzbuf(ii).le.ZERO) zzbuf(ii)=zzbuf(ii-1)
              enddo
              ss%freq_lh(1:inum) = zzbuf(1)*1.0d9         ! -> Hz
              deallocate(zzbuf)
              if(ier.ne.0) go to 1000
           endif

           call mk_rho_lhrf
        endif

     endif

  endif

  !----------------------------------------------------
  !  ICRF power & spectrum absorbed

  if(ss%nicrf_src.gt.0) then
     call rpexist_profile('XGRID_NPHI',iexist)
     if(iexist) then
        call r8_t1profil('XGRID_NPHI',zlabel,zunits,ztime0,zdelta_t, &
             idum,zmbuf,max_zmbuf,igot,ier)  ! load zmbuf(1:igot)
        allocate(zspectrum(igot))
     endif

     call get_coupled_spectrum

     if(iexist) deallocate(zspectrum)
  endif  ! ICRF present

  !----------------------------------------------------
  !  grab heating and current drive profiles recorded in TRANSP run
  !  this will work even if the namelist (splitn) and input data (trdatbuf)
  !  could not be read...

  call rpexist_profile('QICHE',iexist)
  if(iexist) then
     call mk_rho_icrf
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state: allocation for rho_icrf failed.'
        go to 1000
     endif
     !  RF heating profiles detected: electron heating
     call get1p('QICHE',1)
     ss%picrf_totals(1:inx-1,ps_elec_index) = zprofc

     call get1p('QICHMC',1)  ! add "mode conversion" to electron heating
     ss%picrf_totals(1:inx-1,ps_elec_index) = ss%picrf_totals(1:inx-1,ps_elec_index) + zprofc

     call get1p('QICHI',1)   ! direct ion heating
     ss%picth = zprofc

     if(n_rfi.gt.0) then
        call get1p('QMINE',1)   ! electron heating via minority
        ss%pmine = zprofc

        call get1p('QMINI',1)   ! ion heating via minority
        ss%pmini = zprofc
     endif
  endif

  call rpexist_profile('HHCUR',iexist)
  if(iexist) then
     call mk_rho_icrf
     call get1p('HHCUR',2)
     ss%curich = ss%curich + zprofc
  endif

  call rpexist_profile('ICCUR_P',iexist)
  if(iexist) then
     call mk_rho_icrf
     call get1p('ICCUR_P',2)
     ss%curich = ss%curich + zprofc
  endif

  call rpexist_profile('PBI',iexist)
  if(iexist) then
     call mk_rho_nbi
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state: allocation for rho_nbi failed.'
        go to 1000
     endif

     !  beam heating profiles detected
     call get1p('PBE',1)  ! electron heating
     ss%pbe = zprofc

     call get1p('PBI',1)  ! ion heating
     ss%pbi = zprofc

     !  thermalization power...

     call rpexist_profile('PBTHA',jexist)
     if(jexist) then
        call get1p('PBTHA',1) ! rotation sce friction contribution
        ss%pbth = zprofc
     else
        ss%pbth = 0
     endif

     !  add in contributions from each beam species
     do i=1,n_species
        if(bsuffix(i).ne.' ') then
           call get1p('PBTH_'//bsuffix(i),1)
           ss%pbth = ss%pbth + zprofc
        endif
     enddo

     !  halo sources
     if(all_neutrals) then
        call rpexist_profile('P0HALO',jexist)
        if(jexist) then
           call get1p('P0HALO',1)
           ss%pb0_halo = zprofc
           call get1p('PIHALO',1)
           ss%psc_halo = zprofc - ss%pb0_halo
           call get1p('PCXHALO',1)
           ss%pcx_halo = -zprofc  ! source to neutrals = sink to ions
        endif

        call rpexist_profile('TQ0HALO',jexist)
        if(jexist) then
           call get1p('TQ0HALO',1)
           ss%tqb0_halo = zprofc
           call get1p('TQIHALO',1)
           ss%tqsc_halo = zprofc - ss%tqb0_halo
           call get1p('TQCXHALO',1)
           ss%tqcx_halo = -zprofc  ! source to neutrals = sink to ions
        endif
     else
        call rpexist_profile('PBCX',jexist)
        if(jexist) then
           call get1p('PBCX',1)
           ss%pb0_halo = zprofc
           call get1p('PSC_HALO',1)
           ss%psc_halo = zprofc
           call get1p('PCX_HALO',1)
           ss%pcx_halo = zprofc
        endif

        call rpexist_profile('TQBCX',jexist)
        if(jexist) then
           call get1p('TQBCX',1)
           ss%tqb0_halo = zprofc
           call get1p('TQSC_HALO',1)
           ss%tqsc_halo = zprofc
           call get1p('TQCX_HALO',1)
           ss%tqcx_halo = zprofc
        endif
     endif

     !  current drive

     call get1p('CURB',2)
     ss%curbeam = zprofc
  endif

  call rpexist_profile('CURB_P',iexist)
  if(iexist) then
     call mk_rho_nbi
     call get1p('CURB_P',2)
     ss%curbeam = ss%curbeam + zprofc
  endif

  call rpexist_profile('TQBI',iexist)
  if(iexist) then
     call mk_rho_nbi
     !  beam torque profiles detected
     call get1p('TQBE',1)  ! collisional torque to electrons
     ss%tqbe = zprofc

     call get1p('TQBI',1)  ! collisional torque to ions
     ss%tqbi = zprofc

     call get1p('TQBTH',1)   ! thermalization momentum source
     ss%tqbth = zprofc

     call get1p('TQJXB',1)   ! JxB torque
     ss%tqbjxb = zprofc

     call rpexist_profile('TQABDF',iexist)
     if(iexist) then
        write(lunzer(0),*) ' %trx_gen_state: add TQABDF to TQJXB.'
        call get1p('TQABDF',1)   ! JxB torque
        ss%tqbjxb = ss%tqbjxb + zprofc
     endif

  endif

  call rpexist_profile('PFI',iexist)
  if(iexist) then
     call mk_rho_fus
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state: allocation for rho_fus failed.'
        go to 1000
     endif
     !  fusion product heating profiles detected
     call get1p('PFE',1)  ! electron heating
     ss%pfuse = zprofc

     call get1p('PFI',1)  ! ion heating
     ss%pfusi = zprofc

     !  thermalization power...
     ss%pfusth = 0

     !  add in contributions from each beam species
     do i=1,n_species
        if(fsuffix(i).ne.' ') then
           call get1p('PFTH_'//fsuffix(i),1)
           ss%pfusth = ss%pfusth + zprofc
        endif
     enddo

     !  current drive

     call get1p('CURFI',2)
     ss%curfusn = zprofc
  endif

  ! recombination neutral source profiles: heat, momentum
  if(reco_neutrals) then
     call get1p('P0RECO',1)
     ss%p0_reco = zprofc
     call get1p('PIRECO',1)
     ss%psc_reco = zprofc - ss%p0_reco
     call get1p('PCXRECO',1)
     ss%pcx_reco = -zprofc  ! source to neutrals <--> sink to ions

     call get1p('TQ0RECO',1)
     ss%tq0_reco = zprofc
     call get1p('TQIRECO',1)
     ss%tqsc_reco = zprofc - ss%tq0_reco
     call get1p('TQCXRECO',1)
     ss%tqcx_reco = -zprofc  ! source to neutrals <--> sink to ions
  endif

  call rpexist_profile('PEECH',iexist)
  if(iexist.AND.(iflag_ech.le.0)) then
     call mk_rho_ecrf
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state: allocation for rho_ecrf failed.'
        go to 1000
     endif
     ! ECH/ECCD
     call get1p('PEECH',1)  ! electron heating
     ss%peech = zprofc

     call get1p('ECCUR',2)  ! current drive
     ss%curech = zprofc

     do ia=1,ss%necrf_src
        zname=' '
        if(ia.le.9) then
           write(zname,'("PEECH",i1)') ia
        else if(ia.le.99) then
           write(zname,'("PEECH",i2)') ia
        else
           write(zname,'("PEECH",i3)') ia
        endif
        call rpexist_profile(zname,jexist)
        if(jexist) then
           call get1p(zname,1)
           ss%peech_src(:,ia) = zprofc
        endif

        zname=' '
        if(ia.le.9) then
           write(zname,'("ECCUR",i1)') ia
        else if(ia.le.99) then
           write(zname,'("ECCUR",i2)') ia
        else
           write(zname,'("ECCUR",i3)') ia
        endif
        call rpexist_profile(zname,jexist)
        if(jexist) then
           call get1p(zname,2)
           ss%curech_src(:,ia) = zprofc
        endif
     enddo
  endif

  call rpexist_profile('ECCUR_P',iexist)
  if(iexist.AND.(iflag_ech.le.0)) then
     call mk_rho_ecrf
     call get1p('ECCUR_P',2)  ! current drive
     ss%curech = ss%curech + zprofc
  endif

  call rpexist_profile('PELH',iexist)
  if(iexist.AND.(iflag_lh.le.0)) then
     call mk_rho_lhrf
     if(ier.ne.0) then
        write(lunzer(0),*) ' ?trx_gen_state: allocation for rho_lhrf failed.'
        go to 1000
     endif
     ! LHH/LHCD
     call get1p('PELH',1)  ! electron heating
     ss%pelh = zprofc

     call get1p('PILH',1)  ! ion heating
     ss%pilh = zprofc

     call get1p('LHCUR',2)  ! current drive
     ss%curlh = zprofc

     ! NOTE: as of June 2010, no per-antenna LH data in TRANSP yet.

  endif

  call rpexist_profile('LHCUR_P',iexist)
  if(iexist.AND.(iflag_lh.le.0)) then
     call mk_rho_lhrf
     call get1p('LHCUR_P',2)  ! current drive
     ss%curlh = ss%curlh + zprofc
  endif

  !---------------------
  !  profiles indexed by species

  allocate(zn0(inx-1,ss%nspec_gas)); zn0 = ZERO  ! accumulate total n0(:,gas)

  jth = 0
  jnbi = 0
  jfus = 0
  jrf = 0

  if(ss%nrho_nbi.gt.0) then
     ss%rate_sinb0x = 0
     ss%rate_sinb0xs = 0
     ss%rate_sinb0i = 0
  endif

  if(ss%nrho_fus.gt.0) then
     ss%rate_sinf0x = 0
     ss%rate_sinf0xs = 0
     ss%rate_sinf0i = 0
  endif

  if(ss%nspec_tha.eq.0) then
     write(lunzer(0),*) ' %trx_gen_state: SA_TH species list not available.'
     write(lunzer(0),*) '  species indexed profiles not saved (OK).'
  else
     do i=1,n_species
        iz=Zcharg(i)+0.1d0  ! integer charge
        if(itype(i).eq.ps_electron) then

           if(ss%nrho_nbi.gt.0) then

              ss%sbedep = ZERO
              ss%sbehalo = ZERO
              ss%sbsce(1:inx-1,ps_elec_index) = ZERO
              if(ss%nrho_fus.gt.0) then
                 ss%sfsce(1:inx-1,ps_elec_index) = ZERO
              endif

              if(all_neutrals) then
                 call rpexist_profile('SEHALO',iexist)
                 if(iexist) then
                    call get1p('SEHALO',1)
                    ss%sbehalo(1:inx-1) = zprofc
                    ss%sbsce(1:inx-1,ps_elec_index) = zprofc
                 endif
              else
                 call rpexist_profile('SCEV',iexist)
                 if(iexist) then
                    call get1p('SCEV',1)
                    ss%sbehalo(1:inx-1) = zprofc
                    ss%sbsce(1:inx-1,ps_elec_index) = zprofc
                 endif
              endif

              call rpexist_profile('SBE',iexist)
              if(iexist) then
                 call get1p('SBE',1)
                 ss%sbedep(1:inx-1) = zprofc
                 ss%sbsce(1:inx-1,ps_elec_index) = &
                      ss%sbsce(1:inx-1,ps_elec_index) + zprofc
              endif
           else if(ss%nrho_fus.gt.0) then
              ! fusion products present, but, not beam ions...
              ss%sfsce(1:inx-1,ps_elec_index) = ZERO
           endif

           if(reco_neutrals) then
              ! recombination neutrals: net electron source...
              call get1p('S0RECO',1)
              ss%s0reco_e(1:inx-1) = -zprofc
              call get1p('SERECO',1)
              ss%s0reco_e(1:inx-1) = ss%s0reco_e(1:inx-1) + zprofc
           endif

        else if(itype(i).eq.ps_impurity) then
           jth = jth + 1

        else if(itype(i).eq.ps_tokamakium) then
           jth = jth + 2

        else if(itype(i).eq.ps_therm_ion) then
        
           jth = jth + 1
           izloc=iZc(i)
           if(ss%nrho_nbi.gt.0) then

              if(.not.all_neutrals) then
                 zname = 'SV'//trim(thsuffix(i))
                 call get1p(zname,1)
                 ss%sbsce(1:inx-1,jth) = zprofc
              endif

              zname = 'SBCX'//trim(thsuffix(i))
              call get1p(zname,1)
              ss%sb0halo(1:inx-1,jth) = zprofc

              if(all_neutrals) then
                 zname = 'SIHALO_'//trim(thsuffix(i))
                 call get1p(zname,1)
                 ss%sb0halo_recap(1:inx-1,jth) = zprofc
              endif

              inbi=ibmap(i)
              infi=ifmap(i)
              if(inbi.gt.0) then
                 zname = 'SBTH_'//trim(bsuffix(inbi))
                 call get1p(zname,1)
              else if(infi.gt.0) then
                 zname = 'SFTH_'//trim(fsuffix(infi))
                 call get1p(zname,1)
              else
                 zprofc=0.0
              endif

              if(.not.all_neutrals) then
                 ss%sb0halo_recap(1:inx-1,jth) = &
                      ss%sbsce(1:inx-1,jth) - zprofc + ss%sb0halo(1:inx-1,jth)
              else
                 if(infi.gt.0) then
                    ss%sfsce(1:inx-1,jth) = zprofc
                    ss%sbsce(1:inx-1,jth) = ss%sb0halo_recap(1:inx-1,jth) - &
                         ss%sb0halo(1:inx-1,jth)
                 else
                    if(ss%nrho_fus.gt.0) then
                       ss%sfsce(1:inx-1,jth) = ZERO
                    endif
                    ss%sbsce(1:inx-1,jth) = ss%sb0halo_recap(1:inx-1,jth) - &
                         ss%sb0halo(1:inx-1,jth) + zprofc
                 endif
              endif

              if(all_neutrals) then
                 zname = 'N0BH_'//trim(thsuffix(i))
                 call get1(zname)
                 zn0(1:inx-1,jth) = zn0(1:inx-1,jth) + zprofc
                 ss%n0_halo(1:inx-1,jth) = zprofc

                 zname = 'T0BH_'//trim(thsuffix(i))
                 call get1(zname)
                 ss%T0_halo(1:inx-1,jth) = zprofc

                 if(iomega) then
                    zname = 'OM0BH_'//trim(thsuffix(i))
                    call get1(zname)
                    ss%omeg0_halo(1:inx-1,jth) = zprofc
                 else
                    ss%omeg0_halo(1:inx-1,jth) = ZERO
                 endif
              else
                 ! old run, not all neutrals data saved;
                 ! assume "volume" neutrals are beam-dominated...
                 zname = 'DN0V'//trim(thsuffix(i))
                 call get1(zname)
                 zn0(1:inx-1,jth) = zn0(1:inx-1,jth) + zprofc
                 ss%n0_halo(1:inx-1,jth) = zprofc

                 zname = 'T0V'//trim(thsuffix(i))
                 call get1(zname)
                 ss%T0_halo(1:inx-1,jth) = zprofc

                 if(iomega) then
                    zname = 'OM0V'//trim(thsuffix(i))
                    call get1(zname)
                    ss%omeg0_halo(1:inx-1,jth) = zprofc
                 else
                    ss%omeg0_halo(1:inx-1,jth) = ZERO
                 endif
              endif

           else if(ss%nrho_fus.gt.0) then
              ! fusion products present, but, not beam ions...
              infi=ifmap(i)
              if(infi.gt.0) then
                 zname = 'SFTH_'//trim(fsuffix(infi))
                 call get1p(zname,1)
              else
                 zprofc= ZERO
              endif
              ss%sfsce(1:inx-1,jth) = zprofc

           endif

           if(ss%nrho_gas.gt.0) then

              if(all_neutrals) then
                 in0=ss%nspec_gas
                 do ig=1,in0
                    ! here, loop over source gasses; current gas is
                    !   affected specie

                    ! recycling:
                    zsc = max(ONE,srcy(ig))

                    zname='SIRC_'//cgas_abbrev(ig)//'_'//thsuff_1(i)
                    call get1p(zname,1)
                    ss%sprof0(1:inx-1,jth,ig) = zprofc/zsc

                    zname='N0RC_'//cgas_abbrev(ig)//'_'//thsuff_1(i)
                    call get1(zname)
                    ss%n0norm(1:inx-1,jth,ig) = zprofc/zsc
                    zn0(1:inx-1,jth) = zn0(1:inx-1,jth) + zprofc

                    zname='T0RC_'//cgas_abbrev(ig)//'_'//thsuff_1(i)
                    call get1(zname)
                    ss%T0sc0(1:inx-1,jth,ig) = zprofc

                    zname='OM0RC_'//cgas_abbrev(ig)//'_'//thsuff_1(i)
                    call get1(zname)
                    ss%omeg0sc0(1:inx-1,jth,ig) = zprofc

                    ! gas flow:
                    zsc = max(ONE,sgas(ig))

                    zname='SIGF_'//cgas_abbrev(ig)//'_'//thsuff_1(i)
                    call get1p(zname,1)
                    ss%sprof0(1:inx-1,jth,ig+in0) = zprofc

                    zname='N0GF_'//cgas_abbrev(ig)//'_'//thsuff_1(i)
                    call get1(zname)
                    ss%n0norm(1:inx-1,jth,ig+in0) = zprofc/zsc
                    zn0(1:inx-1,jth) = zn0(1:inx-1,jth) + zprofc

                    zname='T0GF_'//cgas_abbrev(ig)//'_'//thsuff_1(i)
                    call get1(zname)
                    ss%T0sc0(1:inx-1,jth,ig+in0) = zprofc

                    zname='OM0GF_'//cgas_abbrev(ig)//'_'//thsuff_1(i)
                    call get1(zname)
                    ss%omeg0sc0(1:inx-1,jth,ig+in0) = zprofc
                 enddo

                 ! here the current gas is the source gas

                 ig=jth

                 ! recycling:

                 zsc = max(ONE,srcy(ig))

                 zname='SERC_'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%sprof0e(1:inx-1,ig) = zprofc/zsc

                 zname='PIRC_'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%qioniz(1:inx-1,ig) = zprofc/zsc

                 zname='TQIRC_'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%tqioniz(1:inx-1,ig) = zprofc/zsc

                 zname='CFPCX_RC'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%qqcx(1:inx-1,ig) = zprofc/zsc

                 zname='CFTCX_RC'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%tqqcx(1:inx-1,ig) = zprofc/zsc

                 zname='T0CX_RC'//cgas_abbrev(ig)
                 call get1(zname)
                 ss%t0cx(1:inx-1,ig) = zprofc

                 zname='OM0CX_RC'//cgas_abbrev(ig)
                 call get1(zname)
                 ss%omeg0cx(1:inx-1,ig) = zprofc

                 ! gas flow:

                 zsc = max(ONE,sgas(ig))

                 zname='SEGF_'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%sprof0e(1:inx-1,ig+in0) = zprofc/zsc

                 zname='PIGF_'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%qioniz(1:inx-1,ig+in0) = zprofc/zsc

                 zname='TQIGF_'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%tqioniz(1:inx-1,ig+in0) = zprofc/zsc

                 zname='CFPCX_GF'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%qqcx(1:inx-1,ig+in0) = zprofc/zsc

                 zname='CFTCX_GF'//cgas_abbrev(ig)
                 call get1p(zname,1)
                 ss%tqqcx(1:inx-1,ig+in0) = zprofc/zsc

                 zname='T0CX_GF'//cgas_abbrev(ig)
                 call get1(zname)
                 ss%t0cx(1:inx-1,ig+in0) = zprofc

                 zname='OM0CX_GF'//cgas_abbrev(ig)
                 call get1(zname)
                 ss%omeg0cx(1:inx-1,ig+in0) = zprofc

              else
                 ! only partial data available on neutral sources...

                 ! for the edge sources-- the same normalized neutral densities
                 ! and temperatures and angular velocities assigned for all
                 ! sources -- more detailed data not available.  This is enough
                 ! however to do a good job on beam CX loss, in case the state
                 ! data is fed to nubeam_comp_exec...

                 zname = 'DN0W'//trim(thsuffix(i))
                 call get1(zname)
                 zn0(1:inx-1,jth) = zn0(1:inx-1,jth) + zprofc
                 do jsc0=1,ss%ngsc0
                    ss%n0norm(1:inx-1,jth,jsc0) = zprofc/sgrcy_sum
                 enddo

                 zname = 'T0W'//trim(thsuffix(i))
                 call get1(zname)
                 do jsc0=1,ss%ngsc0
                    ss%T0sc0(1:inx-1,jth,jsc0) = zprofc
                 enddo

                 if(iomega) then
                    zname = 'OM0W'//trim(thsuffix(i))
                    call get1(zname)
                    do jsc0=1,ss%ngsc0
                       ss%omeg0sc0(1:inx-1,jth,jsc0) = zprofc
                    enddo
                 else
                    do jsc0=1,ss%ngsc0
                       ss%omeg0sc0(1:inx-1,jth,jsc0) = 0
                    enddo
                 endif

                 ! the above code assumes that the neutral gas species indices
                 ! match the low end of the thermal gas species list indices
                 ! i.e. ordering is:
                 !         <low Z thermal speces w/ matching neutral species>
                 !    then <impurities>
              endif

              if(reco_neutrals) then
                 zname = 'S0V0'//trim(thsuffix(i))
                 call get1p(zname,1)
                 ss%s0reco(1:inx-1,jth) = zprofc

                 zname = 'SIV0'//trim(thsuffix(i))
                 call get1p(zname,1)
                 ss%s0reco_recap(1:inx-1,jth) = zprofc

                 zname = 'N0V0_'//trim(thsuffix(i))
                 call get1(zname)
                 ss%n0_reco(1:inx-1,jth) = zprofc
                 zn0(1:inx-1,jth) = zn0(1:inx-1,jth) + zprofc

                 zname = 'T0V0_'//trim(thsuffix(i))
                 call get1(zname)
                 ss%T0_reco(1:inx-1,jth) = zprofc

                 zname = 'OM0V0_'//trim(thsuffix(i))
                 call get1(zname)
                 ss%omeg0_reco(1:inx-1,jth) = zprofc
              endif
           endif

           if((n_bi.gt.0).or.(n_fusi.gt.0)) then

              if(n_bi.gt.0) then
                 call rpexist_profile( &
                      'RSNBI_'//thsuff_1(i)//'_'//bsuff0,rsn_exist)
              else
                 call rpexist_profile( &
                      'RSNFI_'//thsuff_1(i)//'_'//fsuff0,rsn_exist)
              endif
                 
              if(rsn_exist) then
                 ibeami=0
                 ifusi=0

                 do ii=1,n_species
                    if(bsuffix(ii).ne.' ') then
                       ibeami=ibeami+1

                       profname='RSNBI_'//thsuff_1(i)//'_'//bsuffix(ii)

                       call rpexist_profile(profname,rs2_exist)
                       if(.NOT.rs2_exist) cycle  ! no such profiles: impurity beams
                       
                       call get1(profname)
                       ss%rate_sinb0i(:,jth) = ss%rate_sinb0i(:,jth) + zprofc

                       profname='RSNBX_'//thsuff_1(i)//'_'//bsuffix(ii)
                       call get1(profname)
                       ss%rate_sinb0xs(:,jth,ibeami) = zprofc
                       ss%rate_sinb0x(:,jth) = ss%rate_sinb0x(:,jth) + zprofc
                    else if(fsuffix(ii).ne.' ') then
                       ifusi=ifusi+1

                       profname='RSNFI_'//thsuff_1(i)//'_'//fsuffix(ii)
                       call get1(profname)
                       ss%rate_sinf0i(:,jth) = ss%rate_sinf0i(:,jth) + zprofc

                       profname='RSNFX_'//thsuff_1(i)//'_'//fsuffix(ii)
                       call get1(profname)
                       ss%rate_sinf0xs(:,jth,ifusi) = zprofc
                       ss%rate_sinf0x(:,jth) = ss%rate_sinf0x(:,jth) + zprofc
                    endif
                 enddo

              else
                 ! rough reconstruction-- only total sink rates are available
                 !   as sink/[neutral density] for each neutral species

                 ! choose beam or fusion product specie to which to attribute
                 ! the (fast ion neutralizing) charge exchange and impact
                 ! ionization sinks

                 call sinb_select(ibeamx,ifusx,ibeami,ifusi)

                 zname = 'SB0X'//trim(thsuffix(i))
                 call get1(zname)

                 if(ibeamx.gt.0) then
                    do jx=1,inx-1
                       if(zn0(jx,jth).gt.ZERO) then
                          ss%rate_sinb0x(jx,jth) = zprofc(jx)/zn0(jx,jth)
                       else
                          ss%rate_sinb0x(jx,jth) = ZERO
                       endif
                       ss%rate_sinb0xs(jx,jth,ibeamx) = ss%rate_sinb0x(jx,jth)
                    enddo
                 else if(ifusx.gt.0) then
                    do jx=1,inx-1
                       if(zn0(jx,jth).gt.ZERO) then
                          ss%rate_sinf0x(jx,jth) = zprofc(jx)/zn0(jx,jth)
                       else
                          ss%rate_sinf0x(jx,jth) = ZERO
                       endif
                       ss%rate_sinf0xs(jx,jth,ifusx) = ss%rate_sinf0x(jx,jth)
                    enddo
                 endif

                 zname = 'SB0I'//trim(thsuffix(i))
                 call get1(zname)

                 if(ibeami.gt.0) then
                    do jx=1,inx-1
                       if(zn0(jx,jth).gt.ZERO) then
                          ss%rate_sinb0i(jx,jth) = zprofc(jx)/zn0(jx,jth)
                       else
                          ss%rate_sinb0i(jx,jth) = ZERO
                       endif
                    enddo
                 else
                    do jx=1,inx-1
                       if(zn0(jx,jth).gt.ZERO) then
                          ss%rate_sinf0i(jx,jth) = zprofc(jx)/zn0(jx,jth)
                       else
                          ss%rate_sinf0i(jx,jth) = ZERO
                       endif
                    enddo
                 endif
              endif   ! rsn_exist

           endif   ! fast ion specie existence

        else if(itype(i).eq.ps_rf_minority) then
           jrf = jrf + 1

        else if(itype(i).eq.ps_beam_ion) then

           jnbi = jnbi + 1
           zname = 'SBTH_'//trim(bsuffix(i))
           call get1p(zname,1)
           ss%sbtherm(1:inx-1,jnbi) = zprofc

        else if(itype(i).eq.ps_fusion_ion) then

           jfus = jfus + 1
           zname = 'SFTH_'//trim(fsuffix(i))
           call get1p(zname,1)
           ss%sftherm(1:inx-1,jfus) = zprofc

        endif
     enddo
     write(lunzer(0),*) ' %trx_gen_state -- source averaged data used for'
     write(lunzer(0),*) '  recycling and gas flow source neutral densities.'
  endif

  !----------------------------------------------------
  !  Anomalous fast ion transport data...

  if(iflag) then
     call mk_anom(ier)
     if(ier.ne.0) go to 1000
  endif

  !----------------------------------------------------
  !  OK, store the state to file -- or, update without storing

  if(iws.gt.0) then
     call ps_store_plasma_state(ier,filename=fullpath, state=ss)
     if((iws.gt.1).and.(ier.eq.0)) then
        call ps_mdescr_write(fullpath_mdescr,ier, state=ss)
        if(ier.eq.0) then
           call ps_sconfig_write(fullpath_sconfig,ier, state=ss)
        endif
     endif
  else
     call ps_state_memory_update(ier, state=ss)
  endif
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state_geq: ps_store_plasma_state status: ',&
          ier
     go to 1000
  endif

1000 continue
  if(ier.ne.0) then
     write(lunzer(0),*) ' ** error detected in trx_gen_state.f90 ** '
  endif

  if(allocated(zgrid)) deallocate(zgrid,zgridc,zprof,zprofc,zvol,zarea,zaimp,zzimp)
  if(allocated(zprofc1)) deallocate(zprofc1,zprofc2,zprofc3,zprofc4)
  if(allocated(idns)) deallocate(idns,idts,id_eprps,id_eplls)
  if(allocated(itype)) deallocate(itype,iZc,Zcharg,Amass,slbl,ibmap,ifmap)
  if(allocated(gprofc)) deallocate(gprofc,zvpllc)
  if(allocated(abeama)) deallocate(abeama,xzbeama)
  if(allocated(bsuffix)) deallocate(bsuffix,fsuffix,thsuffix)
  if(allocated(ffulla)) deallocate(ffulla,fhalfa)
  if(allocated(kvfrac)) deallocate(kvfrac,ibi)
  if(allocated(cgas_abbrev)) deallocate(cgas_abbrev)
  if(allocated(srcy)) deallocate(sgas,srcy)
  if(allocated(zqsave)) deallocate(zqsave)
  if(allocated(zcharga)) deallocate(zcharga,izprof,zia)
  if(allocated(imap)) deallocate(imap,imapi)
  if(allocated(zrminor)) deallocate(zrminor)

  if(iflag_ech.gt.0) call trx_gen_ech_files(ss, fpath,runid,ier)
  if(iflag_lh.gt.0) call trx_gen_lh_files(ss, fpath,runid,ier)

  contains

    subroutine check_itype_order(ier)

      integer, intent(out) :: ier

      !  check "itype" (species type) array ordering...

      integer :: ii,itype_prev

      itype_prev = -9999

      !-------------------
      ier = 0

      do ii=1,size(itype)
         if(itype_prev.eq.-9999) then
            if(itype(ii).ne.ps_electron) then
               ier = ier + 1
               write(lunzer(0),*) ' %trx_gen_state(check_itype_order): '// &
                    'electrons not at index 1.'
            endif
         else if(itype_prev.eq.ps_electron) then
            if(itype(ii).ne.ps_therm_ion) then
               ier = ier + 1
               write(lunzer(0),*) ' %trx_gen_state(check_itype_order): '// &
                    'thermal ions expected after electrons.'
            endif
         else if(itype_prev.eq.ps_therm_ion) then
            if(itype(ii).eq.ps_electron) then
               ier = ier + 1
               write(lunzer(0),*) ' %trx_gen_state(check_itype_order): '// &
                    'electrons out of order'
            endif
         else if((itype_prev.eq.ps_impurity).or. &
              (itype_prev.eq.ps_tokamakium)) then
            if(itype(ii).eq.ps_electron) then
               ier = ier + 1
               write(lunzer(0),*) ' %trx_gen_state(check_itype_order): '// &
                    'electrons out of order'
            endif
            if(itype(ii).eq.ps_therm_ion) then
               ier = ier + 1
               write(lunzer(0),*) ' %trx_gen_state(check_itype_order): '// &
                    'ions out of order'
            endif
         else
            if(itype(ii).eq.ps_electron) then
               ier = ier + 1
               write(lunzer(0),*) ' %trx_gen_state(check_itype_order): '// &
                    'electrons out of order'
            endif
            if(itype(ii).eq.ps_therm_ion) then
               ier = ier + 1
               write(lunzer(0),*) ' %trx_gen_state(check_itype_order): '// &
                    'ions out of order'
            endif
            if((itype(ii).eq.ps_impurity).or. &
              (itype(ii).eq.ps_tokamakium)) then
               ier = ier + 1
               write(lunzer(0),*) ' %trx_gen_state(check_itype_order): '// &
                    'impurities out of order'
            endif
         endif

         itype_prev = itype(ii)
      enddo

    end subroutine check_itype_order

    logical function matchNZ(prof1,prof2,zero_in)

      ! "non-zero-match" -- compare two profiles; to get match
      ! both must be non-zero and max diff < 1e-4*max(abs. value)
      !   mod -- compare to "zero_in"; a value other than zero could be
      !          passed

      real*8, dimension(:), intent(in) :: prof1,prof2
      real*8, intent(in) :: zero_in

      !---------------------
      real*8 :: maxa1,maxa2,mtol,maxdiff
      !---------------------

      matchNZ = .FALSE.
      if(size(prof1).ne.size(prof2)) return

      maxa1 = maxval(abs(prof1))
      maxa2 = maxval(abs(prof2))

      if((maxa1.eq.ZERO_IN).or.(maxa2.eq.ZERO_IN)) return

      mtol = 1.0d-4*max(maxa1,maxa2)
      if(abs(maxa2-maxa1).gt.mtol) return

      maxdiff = maxval(abs(prof1-prof2))

      if(maxdiff.gt.mtol) return

      matchNZ=.TRUE.

    end function matchNZ

    subroutine mk_rho_nbi

      ier = 0

      if(ss%nrho_nbi.eq.0) then
         ss%nrho_nbi=inx
         call ps_alloc_plasma_state(iertmp, state=ss)
         if(iertmp.ne.0) then
            ier=iertmp
         else
            ss%rho_nbi = zgrid
         endif
      endif

    end subroutine mk_rho_nbi

    subroutine mk_rho_icrf

      ier = 0

      if(ss%nrho_icrf.eq.0) then
         ss%nrho_icrf=inx
         call ps_alloc_plasma_state(iertmp, state=ss)
         if(iertmp.ne.0) then
            ier=iertmp
         else
            ss%rho_icrf = zgrid
         endif
      endif

    end subroutine mk_rho_icrf

    subroutine mk_rho_fus

      ier = 0

      if(ss%nrho_fus.eq.0) then
         ss%nrho_fus=inx
         call ps_alloc_plasma_state(iertmp, state=ss)
         if(iertmp.ne.0) then
            ier=iertmp
         else
            ss%rho_fus = zgrid
         endif
      endif

    end subroutine mk_rho_fus

    subroutine mk_rho_ecrf

      ! Note: only make the grid if iflag_ech.le.0 -- in which case the
      ! ECH/ECCD heating and current drive profiles are read out of the
      ! TRANSP output.  If iflag_ech.gt.0, leave these empty and available
      ! for recomputation by standalone TORAY or GENRAY

      ier = 0

      if(ss%nrho_ecrf.eq.0) then
         if(iflag_ech.le.0) then
            ss%nrho_ecrf=inx
            call ps_alloc_plasma_state(iertmp, state=ss)
            if(iertmp.ne.0) then
               ier=iertmp
            else
               ss%rho_ecrf = zgrid
            endif
         endif
      endif

    end subroutine mk_rho_ecrf

    subroutine mk_rho_lhrf

      ! Note: only make the grid if iflag_lh.le.0 -- in which case the
      ! LH/LHCD heating and current drive profiles are read out of the
      ! TRANSP output.  If iflag_lh.gt.0, leave these empty and available
      ! for recomputation by standalone GENRAY or other code(s).

      ier = 0

      if(ss%nrho_lhrf.eq.0) then
         if(iflag_lh.le.0) then
            ss%nrho_lhrf=inx
            call ps_alloc_plasma_state(iertmp, state=ss)
            if(iertmp.ne.0) then
               ier=iertmp
            else
               ss%rho_lhrf = zgrid
            endif
         endif
      endif

    end subroutine mk_rho_lhrf

    subroutine cur_info(zpcur_diff)

      ! plasma current profile: metadata; check the namelist

      real*8, intent(in) :: zpcur_diff   ! rel. difference btw input data
      ! and simulation, total plasma current

      integer :: isiz_qmod
      integer, dimension(:), allocatable :: nqmoda
      real*8, dimension(:), allocatable :: tqmoda
      real*8 :: ztime
      integer :: ii,indx
      logical :: ilmdif(1)

      !----------------
      ztime = (ss%t0 + ss%t1)/2

      call splitn_getsize('TQMODA',isiz_qmod,ier)
      if(ier.ne.0) return

      allocate(nqmoda(isiz_qmod),tqmoda(isiz_qmod))

      call splitn_iget('NQMODA',isiz_qmod,nqmoda,ier)
      if(ier.ne.0) return

      call splitn_dget('TQMODA',isiz_qmod,tqmoda,ier)
      if(ier.ne.0) return

      if(tqmoda(1).gt.0.9d34) then
         ! use time-invariant control
         call splitn_lget('NLMDIF',1,ilmdif,ier)
         if(ier.ne.0) return

      else  
         if(ztime.lt.tqmoda(1)) then
            indx=1
         else
            indx=isiz_qmod
            do ii=isiz_qmod-1,1,-1
               if(ztime.ge.tqmoda(ii)) then
                  indx=ii+1
                  exit
               endif
            enddo
         endif

         ilmdif(1) = nqmoda(indx).eq.1
      endif

      if(ilmdif(1)) then
         if(zpcur_diff.lt.1.0d-4) then
            ss%cur_data_info = 'TRANSP:total_current_matched:profile_predicted'
         else
            ss%cur_data_info = 'TRANSP:profile_predicted'
         endif
      else
         if(zpcur_diff.lt.1.0d-4) then
            ss%cur_data_info = 'TRANSP:total_current_matched:profile_from_input_data'
         else
            ss%cur_data_info = 'TRANSP:profile_from_input_data'
         endif
      endif

    end subroutine cur_info
               
    subroutine set_zeffdens_info

      !----------------------------------
      ! set Zeff and density metadata...
      !----------------------------------

      integer :: ig,ngmax,inz,iA,iZ,is,iA2,iZ2,indx
      real*8, dimension(:), allocatable :: aplasm,backz,tzefmod
      integer, dimension(:), allocatable :: ndefine,nzefmod
      character*3 nami
      real*8 :: ztime

      LOGICAL :: NLZEFM,NLZFIN,NLZFI2,NLZVBR,NLZVB2,NLZFXI,NLZSIM,NLZEFA

      !----------------------------------
      ! First check for ndefine(ig)=2 thermal ion species
      !   (rare but not unheard of)

      call splitn_getsize('NDEFINE',ngmax,ier)
      if(ier.ne.0) return

      allocate(ndefine(ngmax),aplasm(ngmax),backz(ngmax))

      call splitn_iget('NDEFINE',ngmax,ndefine,ier)
      if(ier.ne.0) return

      call splitn_dget('APLASM',ngmax,aplasm,ier)
      if(ier.ne.0) return

      call splitn_dget('BACKZ',ngmax,backz,ier)
      if(ier.ne.0) return

      call splitn_iget('NGMAX',1,itemp(1),ier)  ! (perhaps) restrict to #used
      ngmax = itemp(1)
      if(ier.ne.0) return

      do ig=1,ngmax
         if(ndefine(ig).eq.2) then
            iA = aplasm(ig) + 0.1d0
            iZ = backZ(ig) + 0.1d0

            nami='?'
            if(iZ.eq.1) then
               if(iA.eq.1) then
                  nami='H'
               else if(iA.eq.2) then
                  nami='D'
               else if(iA.eq.3) then
                  nami='T'
               endif
            else if(iZ.eq.2) then
               if(iA.eq.3) then
                  nami='He3'
               else if(iA.eq.4) then
                  nami='He4'
               endif
            else if(iZ.eq.3) then
               nami='Li'
            endif

            ss%ns_data_info = trim(ss%ns_data_info)//';n'//trim(nami)//'_is_input'
            do is=1,ss%nspec_th
               iz2=ss%q_s(is)/ps_xe + 0.1d0
               ia2=ss%m_s(is)/ps_mp + 0.1d0

               if((iz2.eq.iZ).AND.(ia2.eq.iA)) then
                  ss%ns_is_input(is) = 1
                  exit
               endif
            enddo
         endif
      enddo

      deallocate(ndefine,aplasm,backz)
                  
      !----------------------------------
      ! OK now look at Zeff options

      !----------------
      ztime = (ss%t0 + ss%t1)/2

      call splitn_getsize('TZEFMOD',inz,ier)
      if(ier.ne.0) return

      allocate(nzefmod(inz),tzefmod(inz))

      call splitn_iget('NZEFMOD',inz,nzefmod,ier)
      if(ier.ne.0) return

      call splitn_dget('TZEFMOD',inz,tzefmod,ier)
      if(ier.ne.0) return

      if(tzefmod(1).gt.0.9d34) then
         ! use time-invariant controls

         call splitn_lget('NLZEFM',1,ltemp(1),ier)
         NLZEFM = ltemp(1)
         if(ier.ne.0) return
         call splitn_lget('NLZFIN',1,ltemp(1),ier)
         NLZFIN = ltemp(1)
         if(ier.ne.0) return
         call splitn_lget('NLZFI2',1,ltemp(1),ier)
         NLZFI2 = ltemp(1)
         if(ier.ne.0) return
         call splitn_lget('NLZVBR',1,ltemp(1),ier)
         NLZVBR = ltemp(1)
         if(ier.ne.0) return
         call splitn_lget('NLZVB2',1,ltemp(1),ier)
         NLZVB2 = ltemp(1)
         if(ier.ne.0) return
         call splitn_lget('NLZFXI',1,ltemp(1),ier)
         NLZFXI = ltemp(1)
         if(ier.ne.0) return
         call splitn_lget('NLZSIM',1,ltemp(1),ier)
         NLZSIM = ltemp(1)
         if(ier.ne.0) return
         call splitn_lget('NLZEFA',1,ltemp(1),ier)
         NLZEFA = ltemp(1)
         if(ier.ne.0) return

      else
         if(ztime.lt.tzefmod(1)) then
            indx=1
         else
            indx=inz
            do ii=inz-1,1,-1
               if(ztime.ge.tzefmod(ii)) then
                  indx=ii+1
                  exit
               endif
            enddo
         endif

         call splitn_zeff_switches(nzefmod(indx),1, &
              NLZEFM,NLZFIN,NLZFI2,NLZVBR,NLZVB2,NLZFXI,NLZSIM,NLZEFA)
      endif

      deallocate(tzefmod,nzefmod)

      if((nlzfxi.OR.nlzsim).AND.(.NOT.nlzvbr).AND.(.NOT.nlzvb2)) then
         ss%ns_is_input(ix1:ix2)=1
         ss%ns_data_info = trim(ss%ns_data_info)//';nImpurity_is_input}'
      else
         ss%ns_data_info = trim(ss%ns_data_info)//'}'
      endif

      !  Zeff info

      if(NLZEFM) then
         ss%zeff_data_info='TRANSP:Zeff_from_resistivity_Vsur_match'
      else if(NLZFIN) then
         if(nlzefa) then
            ss%zeff_data_info='TRANSP:Zeff_from_scalar_input;profile_is_shaped'
         else
            ss%zeff_data_info='TRANSP:flat_Zeff_from_scalar_input'
         endif
      else if(NLZFI2) then
         ss%zeff_data_info='TRANSP:Zeff_from_profile_input'
      else if(NLZVBR) then
         if(nlzefa) then
            ss%zeff_data_info='TRANSP:Zeff_from_VB_chord;profile_is_shaped'
         else if(nlzfxi.or.nlzsim) then
            ss%zeff_data_info='TRANSP:Zeff_from_VB_chord;shape_from_impurity_density'
         else
            ss%zeff_data_info='TRANSP:flat_Zeff_from_VB_chord'
         endif
      else if(NLZVB2) then
         ss%zeff_data_info='TRANSP:Zeff_from_VB_profile'
      else if(nlzfxi) then
         ss%zeff_data_info='TRANSP:Zeff_from_impurity_density'
      else if(nlzsim) then
         ss%zeff_data_info='TRANSP:Zeff_from_impurity_densities'
      else
         ss%zeff_data_info='TRANSP:Zeff_from_unknown_source'
      endif

    end subroutine set_zeffdens_info

    subroutine set_vprof_info

      ! set metadata information for velocity profiles

      logical :: NLVWNC(1)
      integer :: int,it
      real*8 :: ztime

      !-----------------

      if(.not.iomega) then
         ss%vtor_data_info='TRANSP:no_data'
         ss%vpol_data_info='TRANSP:no_data'
         return
      endif

      !  toroidal angular velocity information exists...

      ztime = (ss%t0 + ss%t1)/2

      call splitn_lget('NLVWNC',1,NLVWNC,ier)    ! NC analysis flag
      if(ier.ne.0) return

      ! apply labeling

! fmp - probably need to adapt to pt-solver
!     if(nvphmod.gt.0) then     
!        ss%vtor_data_info='TRANSP:predicted:angular_momentum_equation'
!        ss%vpol_data_info='TRANSP:input_profile_not_used'

!     else if(nlvwnc(1)) then
      if(nlvwnc(1)) then
         if(ivtr1.eq.ix1) then
            ss%vtor_is_input(ix1:ix2)=1
            ss%vtor_data_info='TRANSP:impurity_toroidal_rotation_data'
         else if(ivtr1.eq.-1) then
            ss%vtor_data_info='TRANSP:generic_toroidal_angular_velocity_data'
         else
            ss%vtor_data_info='TRANSP:specie_toroidal_rotation_data'
            ss%vtor_is_input(ivtr1:ivtr2)=1
         endif

         if(ivpl1.eq.ix1) then
            ss%vpol_is_input(ix1:ix2)=1
            ss%vpol_data_info='TRANSP:impurity_poloidal_rotation_data'
         else if(ivpl1.eq.-1) then
            ss%vpol_data_info='TRANSP:input_profile_not_used'
         else
            ss%vpol_data_info='TRANSP:specie_poloidal_rotation_data'
            ss%vpol_is_input(ivpl1:ivpl2)=1
         endif
      else
         ss%vtor_data_info='TRANSP:generic_toroidal_angular_velocity_data'
         ss%vpol_data_info='TRANSP:input_profile_not_used'
      endif

    end subroutine set_vprof_info

    subroutine get_vtor_vpol

      ! get toroidal and poloidal velocities on outer half plane
      ! record presence of input data profiles (if any)
      ! RPLOT calculator expression (rpcal0 call) used to extract
      ! outer half plane variation from data vs. major radius

      ! the code is dependent on the TRANSP output names being set as
      ! they are in outcor/plotgen.i

      character*30, dimension(:), allocatable :: vtor_names
      character*30, dimension(:), allocatable :: vpol_names
      character*30 :: tmpname,testname
      character*200 :: calc_expr
      integer, dimension(:), allocatable :: isigns

      character*64 :: mglbl,mguns

      integer :: istype,isize_vtor,isize_vpol,ii,jj,ispec_nc,ispec_data
      integer :: iwarn1,iwarn2

      logical :: iexist_nc

      !---------------------

      call rpsize_multi('VTORMP',iexist,istype,isize_vtor)
      if(isize_vtor.eq.0) return

      call rpsize_multi('VPOLMP',iexist,istype,isize_vpol)
      if(isize_vpol.eq.0) return

      allocate(vtor_names(isize_vtor))
      allocate(vpol_names(isize_vpol))
      allocate(isigns(max(isize_vtor,isize_vpol)))

      call rpmulti('VTORMP',istype,mglbl,mguns,isize_vtor,isigns,vtor_names, &
           ier)
      if(ier.ne.0) return

      call rpmulti('VPOLMP',istype,mglbl,mguns,isize_vpol,isigns,vpol_names, &
           ier)
      if(ier.ne.0) return

      ! loop over VTORMP names; for those that match known thermal species
      ! compute the outer midplane profile and store.  Also, look for input
      ! data profiles.

      ! VTORMP...

      ispec_data=-2

      do jj = 1,isize_vtor
         if(vtor_names(jj).eq.'VTOR_AVG') cycle
         
         imatch=0
         ispec_nc=-2

         if(vtor_names(jj).eq.'VTORE_NC') then
            imatch=1
            ispec_nc=0
         else if(vtor_names(jj).eq.'VTORX_NC') then
            imatch=1
            ispec_nc=-1
         else if(vtor_names(jj).eq.'VTORE') then
            imatch=1
            if(ispec_data.gt.-2) ispec_data=-3
            if(ispec_data.eq.-2) ispec_data=0
         else if(vtor_names(jj).eq.'VTORX') then
            imatch=1
            if(ispec_data.gt.-2) ispec_data=-3
            if(ispec_data.eq.-2) ispec_data=-1
         else
            do ii=imj1,imj2
               if(thsuff_1(ii+1).eq.'L') then
                  if(vtor_names(jj).eq.'VTORLINC') then
                     imatch=1
                     ispec_nc=ii
                  else if(vtor_names(jj).eq.'VTORLI') then
                     imatch=1
                     if(ispec_data.gt.-2) ispec_data=-3
                     if(ispec_data.eq.-2) ispec_data=ii
                  endif
               else
                  ! +1 in index because electrons are @1 in thsuff_1(...)
                  !    order match of non-impurity thermal species is assured,
                  !    see check_itype_order & species extraction loops...
                  tmpname='VTOR'//thsuff_1(ii+1)
                  if(vtor_names(jj).eq.trim(tmpname)//'_NC') then
                     imatch=1
                     ispec_nc=ii
                  else if(vtor_names(jj).eq.tmpname) then
                     imatch=1
                     if(ispec_data.gt.-2) ispec_data=-3
                     if(ispec_data.eq.-2) ispec_data=ii
                  endif
               endif
            enddo
         endif

         if(imatch.eq.0) then
            ! presumed generic impurity element name e.g. VTORC12
            if(ispec_data.gt.-2) ispec_data=-3
            if(ispec_data.eq.-2) ispec_data=ii
         endif

         if(ispec_nc.gt.-2) then
            ! v_phi, outer half midplane
            tmpname= 'X'//trim(vtor_names(jj))
            call rpexist_profile(tmpname,iexist_nc)

            if(iexist_nc) then
               iwarn1=0
               iwarn2=0
            else
               calc_expr = trim(tmpname)//',"Mapped NC velocity",cm/sec = '// &
                    '%XMAP("CTR","OUT",'//trim(vtor_names(jj))//')'
               call rpcal0(calc_expr,iwarn1,iwarn2)
            endif

            if(max(iwarn1,iwarn2).gt.0) then
               write(lunzer(0),*) ' ---------------------- '
               write(lunzer(0),*) ' outer half plan mapping error:'
               write(lunzer(0),*) ' '//trim(calc_expr)
               write(lunzer(0),*) ' (velocity profile data ignored).'
               write(lunzer(0),*) ' ---------------------- '
            else
               call get1(tmpname)
               if(ispec_nc.eq.-1) then
                  do ii=ix1,ix2
                     ss%vtor_omp(:,ii) = zprofc
                  enddo
               else
                  ss%vtor_omp(:,ispec_nc) = zprofc
               endif
            endif

            ! v_phi, inner half midplane
            tmpname= 'Y'//trim(vtor_names(jj))
            call rpexist_profile(tmpname,iexist_nc)

            if(iexist_nc) then
               iwarn1=0
               iwarn2=0
            else
               calc_expr = trim(tmpname)//',"Mapped NC velocity",cm/sec = '// &
                    '%XMAP("CTR","IN",'//trim(vtor_names(jj))//')'
               call rpcal0(calc_expr,iwarn1,iwarn2)
            endif

            if(max(iwarn1,iwarn2).gt.0) then
               write(lunzer(0),*) ' ---------------------- '
               write(lunzer(0),*) ' inner half plan mapping error:'
               write(lunzer(0),*) ' '//trim(calc_expr)
               write(lunzer(0),*) ' (velocity profile data ignored).'
               write(lunzer(0),*) ' ---------------------- '
            else
               call get1(tmpname)
               if(ispec_nc.eq.-1) then
                  do ii=ix1,ix2
                     ss%vtor_inmp(:,ii) = zprofc
                  enddo
               else
                  ss%vtor_inmp(:,ispec_nc) = zprofc
               endif
            endif
         endif

      enddo

      if(ispec_data.eq.-3) then
         write(lunzer(0),*) ' ---------------------- '
         write(lunzer(0),*) ' VTORMP: input data ambiguity due to '
         write(lunzer(0),*) '   multiple input profiles; input ID skipped.'
         write(lunzer(0),*) ' ---------------------- '
      else if(ispec_data.eq.-1) then
         ivtr1=ix1
         ivtr2=ix2
      else if(ispec_data.gt.-2) then
         ivtr1=ispec_data
         ivtr2=ispec_data
      endif

      ! VPOLMP...

      ispec_data=-2

      do jj = 1,isize_vpol
         if(vpol_names(jj).eq.'VPOL_AVG') cycle
         
         imatch=0
         ispec_nc=-2

         if(vpol_names(jj).eq.'VPOLE_NC') then
            imatch=1
            ispec_nc=0
         else if(vpol_names(jj).eq.'VPOLX_NC') then
            imatch=1
            ispec_nc=-1
         else if(vpol_names(jj).eq.'VPOLE') then
            imatch=1
            if(ispec_data.gt.-2) ispec_data=-3
            if(ispec_data.eq.-2) ispec_data=0
         else if(vpol_names(jj).eq.'VPOLX') then
            imatch=1
            if(ispec_data.gt.-2) ispec_data=-3
            if(ispec_data.eq.-2) ispec_data=-1
         else
            do ii=imj1,imj2
               if(thsuff_1(ii+1).eq.'L') then
                  if(vpol_names(jj).eq.'VPOLLINC') then
                     imatch=1
                     ispec_nc=ii
                  else if(vpol_names(jj).eq.'VPOLLI') then
                     imatch=1
                     if(ispec_data.gt.-2) ispec_data=-3
                     if(ispec_data.eq.-2) ispec_data=ii
                  endif
               else
                  ! +1 in index because electrons are @1 in thsuff_1(...)
                  !    order match of non-impurity thermal species is assured,
                  !    see check_itype_order & species extraction loops...
                  tmpname='VPOL'//thsuff_1(ii+1)
                  if(vpol_names(jj).eq.trim(tmpname)//'_NC') then
                     imatch=1
                     ispec_nc=ii
                  else if(vpol_names(jj).eq.tmpname) then
                     imatch=1
                     if(ispec_data.gt.-2) ispec_data=-3
                     if(ispec_data.eq.-2) ispec_data=ii
                  endif
               endif
            enddo
         endif

         if(imatch.eq.0) then
            ! presumed generic impurity element name e.g. VPOLC12
            if(ispec_data.gt.-2) ispec_data=-3
            if(ispec_data.eq.-2) ispec_data=ii
         endif

         if(ispec_nc.gt.-2) then
            ! outer half midplane
            tmpname= 'X'//trim(vpol_names(jj))
            call rpexist_profile(tmpname,iexist_nc)

            if(iexist_nc) then
               iwarn1=0
               iwarn2=0
            else
               calc_expr = trim(tmpname)//',"Mapped NC velocity",cm/sec = '// &
                    '%XMAP("CTR","OUT",'//trim(vpol_names(jj))//')'
               call rpcal0(calc_expr,iwarn1,iwarn2)
            endif

            if(max(iwarn1,iwarn2).gt.0) then
               write(lunzer(0),*) ' ---------------------- '
               write(lunzer(0),*) ' outer half plan mapping error:'
               write(lunzer(0),*) ' '//trim(calc_expr)
               write(lunzer(0),*) ' (velocity profile data ignored).'
               write(lunzer(0),*) ' ---------------------- '
            else
               call get1(tmpname)
               if(ispec_nc.eq.-1) then
                  do ii=ix1,ix2
                     ss%vpol_omp(:,ii) = zprofc
                  enddo
               else
                  ss%vpol_omp(:,ispec_nc) = zprofc
               endif
            endif

            ! inner half midplane
            tmpname= 'Y'//trim(vpol_names(jj))
            call rpexist_profile(tmpname,iexist_nc)

            if(iexist_nc) then
               iwarn1=0
               iwarn2=0
            else
               calc_expr = trim(tmpname)//',"Mapped NC velocity",cm/sec = '// &
                    '%XMAP("CTR","IN",'//trim(vpol_names(jj))//')'
               call rpcal0(calc_expr,iwarn1,iwarn2)
            endif

            if(max(iwarn1,iwarn2).gt.0) then
               write(lunzer(0),*) ' ---------------------- '
               write(lunzer(0),*) ' inner half plan mapping error:'
               write(lunzer(0),*) ' '//trim(calc_expr)
               write(lunzer(0),*) ' (velocity profile data ignored).'
               write(lunzer(0),*) ' ---------------------- '
            else
               call get1(tmpname)
               if(ispec_nc.eq.-1) then
                  do ii=ix1,ix2
                     ss%vpol_inmp(:,ii) = zprofc
                  enddo
               else
                  ss%vpol_inmp(:,ispec_nc) = zprofc
               endif
            endif
         endif

      enddo

      if(ispec_data.eq.-3) then
         write(lunzer(0),*) ' ---------------------- '
         write(lunzer(0),*) ' VPOLMP: input data ambiguity due to '
         write(lunzer(0),*) '   multiple input profiles; input ID skipped.'
         write(lunzer(0),*) ' ---------------------- '
      else if(ispec_data.eq.-1) then
         ivpl1=ix1
         ivpl2=ix2
      else if(ispec_data.gt.-2) then
         ivpl1=ispec_data
         ivpl2=ispec_data
      endif

    end subroutine get_vtor_vpol

    subroutine set_bdy(prof,bdy)

      ! set boundary value from profile.
      ! this is a half zone width linear extrapolation, except modified
      ! as necessary to prevent sign change in profile

      real*8, dimension(:), intent(in) :: prof
      real*8, intent(out) :: bdy

      !-------------------
      integer :: isize,ism1
      real*8 :: zprod,zxtrap,zlim
      !-------------------

      isize = size(prof)
      ism1 = isize-1

      ! half zone width linear extrapolation:
      zxtrap = prof(isize) + 0.5d0*(prof(isize)-prof(ism1))

      ! sign check
      zprod = prof(isize)*prof(ism1)

      if(zprod.le.0.0) then
         bdy = zxtrap   ! just use (sign already changes or all zero)

      else
         ! last 2 pts are of same sign, non-zero;
         ! constrain extrapolation to prevent sign change

         zlim = 0.1d0*prof(isize)
         if(zlim.lt.ZERO) then
            bdy = min(zxtrap,zlim)  ! prevent extrap going positive
         else
            bdy = max(zxtrap,zlim)  ! prevent extrap going negative
         endif

      endif

    end subroutine set_bdy

    subroutine getnum_heater(chanID,namlnum,namlsw,inum,inum_trdat,iwarn)

      ! get number of heating elements (beams or RF antennas)

      character*(*), intent(in) :: chanID  ! trdat channel ID
      character*(*), intent(in) :: namlnum ! number in namelist
      character*(*), intent(in) :: namlsw  ! switch in namelist or " "

      integer, intent(out) :: inum         ! number of items
      integer, intent(out) :: inum_trdat   ! number of trdat channels

      integer, intent(inout) :: iwarn      ! warning flag

      ! (a) if namlsw is present and .FALSE., inum=0; if absent switchval TRUE
      ! (b) if switchval is TRUE:
      !     1) get number from namelist -- this will be the number of items
      !     2) get number in trdat
      !     3) if trdat number > 0, it must match namelist number; set warning
      !        flag if this condition is not satisfied.

      ! if chanID = "NB", trace beam support is used

      !--------------------------
      ! local:
      logical :: switchval(1)
      integer :: ierloc,isizb,inb,ib0
      real*8 :: trace_min = 1.0d-5
      integer :: itemp(1)
      !--------------------------

      inum=0
      inum_trdat = 0

      ! (a) ...

      switchval(1) = .TRUE.
      if(namlsw.ne.' ') then
         call splitn_lget(namlsw,1,switchval,ierloc)
         if(ierloc.ne.0) then
            write(lunzer(0),*) ' ?trx_gen_state: splitn read error on: '// &
                 trim(namlsw)
            iwarn=1
            switchval(1) = .FALSE.
         endif
      endif

      if(.not.switchval(1)) return

      ! (b) ...

      itemp(1) = inum
      call splitn_iget(namlnum,1,itemp,ierloc)
      inum = itemp(1)
      if(ierloc.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state: splitn read error on: '// &
              trim(namlnum)
         inum=0
         iwarn=1
         return
      endif

      inum_trdat = tdb_nchan_find(d,chanID)

      !-------------------------------
      !  neutral beams only: inum_trdat < inum could occur due to presence
      !  of trace beams (this is rare but needs to be checked)

      if(chanID.eq.'NB') then
         call splitn_getsize('NTRACE',isizb,ierloc)
         if(ierloc.ne.0) then
            write(lunzer(0),*) ' ?trx_gen_state: splitn read error on: '// &
                 trim(namlsw)
            iwarn=1
            allocate(ntrace(inum)); ntrace=0
            allocate(ftrace(inum)); ftrace=ZERO
            allocate(ftraci(inum)); ftraci=ZERO
            allocate(imap(inum))
            allocate(imapi(inum))
            do ibc=1,inum
               imap(ibc)=ibc
               imapi(ibc)=ibc
            enddo
         else
            allocate(ntrace(isizb),ftrace(isizb))
            call splitn_iget('NTRACE',isizb,ntrace,ierloc)
            if(ierloc.ne.0) then
               write(lunzer(0),*) ' ?trx_gen_state: splitn NTRACE read error!'
               iwarn=1
               inum=0
               inum_trdat=0
               return
            endif
            call splitn_dget('FTRACE',isizb,ftrace,ierloc)
            if(ierloc.ne.0) then
               write(lunzer(0),*) ' ?trx_gen_state: splitn FTRACE read error!'
               iwarn=1
               inum=0
               inum_trdat=0
               return
            endif

            !  count number of main beams and trace beams
            allocate(imap(isizb),imapi(isizb)); imap=0; imapi=0
            allocate(ftraci(isizb)); ftraci=ZERO

            intrace=0
            inb=0

            do ibc=1,isizb
               if(ibc.gt.inum) exit
               if(abs(ntrace(ibc)).ne.0) then
                  intrace=intrace+1
               else
                  inb=inb+1
                  imap(ibc)=inb
                  imapi(inb)=ibc
               endif
            enddo

            if(intrace.gt.0) then
               write(lunzer(0),*) &
                    ' *** trx_gen_state: ',intrace,' trace beams! '
               do ibc=1,isizb
                  if(ibc.gt.inum) exit
                  ib0=abs(ntrace(ibc))
                  if(ib0.ne.0) then
                     if(imap(ib0).eq.0) then
                        write(lunzer(0),*) ' ?trx_gen_state: trace beam error.'
                        iwarn=1
                        inum=0
                        inum_trdat=0
                        return
                     endif
                     ib0=imap(ib0)
                     ftraci(ib0)=max(trace_min,ftrace(ibc))
                  endif
               enddo

               inum=inum-intrace  ! corrected # of main beams

            endif
         endif ! NTRACE size check
      endif ! "NB" channel check

      if(inum_trdat.gt.0) then
         if(inum.ne.inum_trdat) then
            write(lunzer(0),*) ' ?trx_gen_state: apparent inconsistency btw TRDAT and namelist:'
            write(lunzer(0),*) '  '//trim(namlnum)//' = ',inum
            write(lunzer(0),*) '  trdat '//trim(chanID)//' = ',inum_trdat
            iwarn=1
         else if(inum.gt.max_zmbuf) then
            write(lunzer(0),*) ' ?trx_gen_state: too many heating sources: ',inum
            inum=max_zmbuf
            iwarn=1
         endif
      endif

    end subroutine getnum_heater

    subroutine splitn_blabel(tname, prefix, inum, cvals, imapi, ier)

      ! get namelist labels for aux. heating devices

      character*(*) :: tname       ! TRANSP namelist name
      character*(*) :: prefix      ! default prefix for naming
      integer, intent(in) :: inum  ! number of item names to be returned
      character*(*) :: cvals(*)    ! labels returned
      integer, dimension(:) :: imapi  ! index mapping

      integer, intent(out) :: ier  ! status code returned

      !-------------------
      integer :: ii,irank,idims(2,10)
      character*32, dimension(:), allocatable :: chbuf
      !-------------------

      ier = 0
      if(inum.eq.0) return

      idims=0
      call splitn_getdims(tname,irank,10,idims,isize,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: error finding "', &
              tname,'" in namelist.'
         return
      endif
      if(isize.lt.inum) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: "',tname, &
              '" namelist array size: ',isize,'; needed: ',inum
         ier=99
         return
      endif

      allocate(chbuf(isize))
      call splitn_cget(tname,len(chbuf(1)),isize,chbuf,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: error fetching "', &
              tname,'" values from namelist.'
         return
      endif
      
      do i=1,inum
         ii=imapi(i)
         if(chbuf(ii).ne.' ') then
            cvals(i)=chbuf(ii)
         else
            if(i.le.999) then
               write(cvals(i),'(A,I3.3)') trim(prefix),i
            else if(i.le.9999) then
               write(cvals(i),'(A,I4.4)') trim(prefix),i
            else
               write(cvals(i),'(A,I5.5)') trim(prefix),i
            endif
         endif
      enddo

      deallocate(chbuf)
    end subroutine splitn_blabel

    subroutine splitn_dvec(tname, inum, vals, ier)

      ! get 1d vector -- first inum values of namelist 1d vector

      character*(*) :: tname       ! TRANSP namelist name
      integer, intent(in) :: inum  ! number of item names to be returned
      real*8 :: vals(*)            ! values returned
      integer, intent(out) :: ier  ! status code returned

      !-------------------
      integer :: irank,idims(2,10)
      real*8, dimension(:), allocatable :: dbuf
      !-------------------

      ier = 0
      if(inum.eq.0) return

      idims=0
      call splitn_getdims(tname,irank,10,idims,isize,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: error finding "', &
              tname,'" in namelist.'
         return
      endif
      if(isize.lt.inum) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: "',tname, &
              '" namelist array size: ',isize,'; needed: ',inum
         ier=99
         return
      endif

      allocate(dbuf(isize))
      call splitn_dget(tname,isize,dbuf,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: error fetching "', &
              tname,'" values from namelist.'
         return
      endif
      
      vals(1:inum) = dbuf(1:inum)

      deallocate(dbuf)
    end subroutine splitn_dvec

    subroutine get_onoff(zid,inum,ier)
      character*(*), intent(in) :: zid    ! id of items for which on/off times
      !  are sought -- trdatbuf code -- e.g. "NB" for neutral beams.
      integer, intent(in) :: inum         ! number of items
      integer, intent(out) :: ier         ! status return code, 0=OK

      character*32 zname1,zname2
      real*8 :: ptfix,dttfix,buf(1)

      zname1 = 'p'//trim(zid)//'tfix'
      zname2 = 'dt'//trim(zid)//'tfix'

      call splitn_dget(zname1,1,buf,ier); if(ier.ne.0) return
      ptfix = buf(1)

      call splitn_dget(zname2,1,buf,ier); if(ier.ne.0) return
      dttfix = buf(1)

      call tdb_onoff_times(d,zid,ptfix,dttfix, &
           zpwr%tbon(1:inum),zpwr%tboff(1:inum),ier)
      if(ier.ne.0) return

    end subroutine get_onoff

    subroutine get1(pname)

      character*(*), intent(in) :: pname  ! name of zone oriented profile to 
      ! fetch -- note, MKS conversion is done in trx_prof

      call eq_gfnum(pname,id1)
      if(id1.eq.0) then
         call trx_prof(pname,zunits,iorder,id1,ier)
         if(ier.ne.0) id1=0
      endif

      if(id1.eq.0) then
         ier=1
         write(lunzer(0),*) ' ?trx_gen_state(get1): failed to find: ', &
              trim(pname)
         zprofc=0

      else
         call eq_rgetf(inx-1,zgridc,id1,0,zprofc,ier)
         if(ier.ne.0) zprofc=0
      endif

    end subroutine get1

    subroutine get1b(pname)

      character*(*), intent(in) :: pname  ! name of bdy oriented profile to 
      ! fetch -- note, MKS conversion is done in trx_prof

      call eq_gfnum(pname,id1)
      if(id1.eq.0) then
         call trx_prof(pname,zunits,iorder,id1,ier)
         if(ier.ne.0) id1=0
      endif

      if(id1.eq.0) then
         ier=1
         write(lunzer(0),*) ' ?trx_gen_state(get1): failed to find: ', &
              trim(pname)
         zprof=0

      else
         call eq_rgetf(inx,zgrid,id1,0,zprof,ier)
         if(ier.ne.0) zprof=0
      endif

    end subroutine get1b

    subroutine get1p(pname,inorm)

      character*(*), intent(in) :: pname  ! name of profile to fetch
      integer :: inorm   ! =1 -- multiply *dV; =2 -- multiply *dA

      character*32 zunits

      call eq_gfnum(pname,id1)
      if(id1.eq.0) then
         call trx_prof(pname,zunits,iorder,id1,ier)
         if(ier.ne.0) id1=0
      endif

      if(id1.eq.0) then
         ier=1
         write(lunzer(0),*) ' %trx_gen_state(get1): failed to find: ', &
              trim(pname)
         zprofc=0

      else
         call eq_rgetf(inx-1,zgridc,id1,0,zprofc,ier)
         if(ier.ne.0) zprofc=0
      endif

      if(inorm.eq.1) then
         zprofc = zprofc*(zvol(2:inx)-zvol(1:inx-1))
      else if(inorm.eq.2) then
         zprofc = zprofc*(zarea(2:inx)-zarea(1:inx-1))
      else
         ier=99
         write(lunzer(0),*) ' ?trx_gen_state_geq(get1p) unexpected inorm: ',inorm
      endif

    end subroutine get1p

    subroutine set_suffix3(iz,ia,lbl3)
      integer, intent(in) :: iz,ia   ! Z & A
      character*3, intent(out) :: lbl3  ! "H","D","T","HE3","HE4", or blank

      lbl3 = ' '

      if(iz.eq.1) then
         if(ia.eq.1) then
            lbl3='H'
         else if(ia.eq.2) then
            lbl3='D'
         else if(ia.eq.3) then
            lbl3='T'
         endif
      else if(iz.eq.2) then
         if(ia.eq.3) then
            lbl3='HE3'
         else if(ia.eq.4) then
            lbl3='HE4'
         endif
      else if(iz.eq.3) then
         lbl3='LI'
      endif

    end subroutine set_suffix3

    subroutine set_suffix1(iz,ia,lbl1,cfb,lbl0)
      integer, intent(in) :: iz,ia   ! Z & A
      character*1, intent(out) :: lbl1  ! "H","D","T","3","4", or blank
      character*1, intent(in) :: cfb
      character*1, intent(inout) :: lbl0  ! set if non-blank

      !  supplement -- DMC May 2010 -- recognize impurity beam species
      !   N:Neon, A:Argon, K:Krypton, X:Xenon
      !    Note, these do not set lbl0

      lbl1 = ' '

      if(iz.eq.1) then
         if(ia.eq.1) then
            if(cfb.eq.'F') then
               lbl1='P'  ! "proton"
            else
               lbl1='H'  ! "H" nucleus
            endif
         else if(ia.eq.2) then
            lbl1='D'
         else if(ia.eq.3) then
            lbl1='T'
         endif
      else if(iz.eq.2) then
         if(ia.eq.3) then
            lbl1='3'
         else if(ia.eq.4) then
            lbl1='4'
         endif
      else if(iz.eq.3) then
         lbl1='L'
      endif

      if(lbl1.ne.' ') then
         if(lbl0.eq.' ') lbl0=lbl1
      endif

      if(lbl1.eq.' ') then
         if(iz.eq.10) then
            lbl1='N'
         else if(iz.eq.18) then
            lbl1='A'
         else if(iz.eq.36) then
            lbl1='K'
         else if(iz.eq.54) then
            lbl1='X'
         endif
      endif

    end subroutine set_suffix1

    subroutine nbget_alg0(vecnam,psvec)

      ! get array of machine beam description data -- no scalar fallback
      ! variable in TRANSP namelist definition.

      character*(*), intent(in) :: vecnam ! TRANSP namelist vector name
      real*8, intent(out), dimension(:) :: psvec  ! plasma state data (output)

      !-------------------
      integer :: irank,idims(2,10),ii
      real*8, dimension(:), allocatable :: dbuf
      !-------------------

      ier = 0
      if(inum.eq.0) return

      idims=0
      call splitn_getdims(vecnam,irank,10,idims,isize,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) &
              ' ?trx_gen_state(nbget_alg0): namelist read failure: '// &
              'size of: '//trim(vecnam)
         return
      else if(irank.ne.1) then
         ier=99
         write(lunzer(0),*) &
              ' ?trx_gen_state(nbget_alg0): namelist vec unexpected rank: ', &
              irank,' for: '//trim(vecnam)
         return
      endif

      allocate(dbuf(isize))

      do
         call splitn_dget(vecnam,isize,dbuf,ier)
         if(ier.ne.0) then
            write(lunzer(0),*) &
                 ' ?trx_gen_state(nbget_alg1): namelist read failure: '// &
                 trim(vecnam)
            exit
         endif

         do ib=1,inum+intrace
            ib0=imap(ib)
            if(ib0.gt.0) then
               psvec(ib0)=zconv*dbuf(ib)
            endif
         enddo

         exit
      enddo

      deallocate(dbuf)

    end subroutine nbget_alg0

    subroutine nbget_alg1(scnam,vecnam,psvec)

      ! get beam machine description data from splitn (TRANSP namelist)
      ! "algorithm 1" rule:  zero-valued elements of "vecnam" are replaced
      ! with the "scnam" value.

      ! zconv converts the data to MKS units for the plasma state

      character*(*), intent(in) :: scnam  ! TRANSP namelist scalar name
      character*(*), intent(in) :: vecnam ! TRANSP namelist vector name

      real*8, intent(out), dimension(:) :: psvec  ! plasma state data (output)

      !-------------------
      integer :: irank,idims(2,10),ii
      real*8, dimension(:), allocatable :: dbuf
      real*8 :: scval(1)
      !-------------------

      ier = 0
      if(inum.eq.0) return

      idims=0
      call splitn_getdims(vecnam,irank,10,idims,isize,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) &
              ' ?trx_gen_state(nbget_alg1): namelist read failure: '// &
              'size of: '//trim(vecnam)
         return
      else if(irank.ne.1) then
         ier=99
         write(lunzer(0),*) &
              ' ?trx_gen_state(nbget_alg1): namelist vec unexpected rank: ', &
              irank,' for: '//trim(vecnam)
         return
      endif

      allocate(dbuf(isize))

      do
         call splitn_dget(vecnam,isize,dbuf,ier)
         if(ier.ne.0) then
            write(lunzer(0),*) &
                 ' ?trx_gen_state(nbget_alg1): namelist read failure: '// &
                 trim(vecnam)
            exit
         endif

         call splitn_dget(scnam,1,scval,ier)
         if(ier.ne.0) then
            write(lunzer(0),*) &
                 ' ?trx_gen_state(nbget_alg1): namelist read failure: '// &
                 trim(scnam)
            exit
         endif
      
         do ii=1,inum+intrace
            ib0=imap(ii)
            if(ib0.gt.0) then
               if(dbuf(ii).eq.ZERO) then
                  psvec(ib0)=scval(1)*zconv
               else
                  psvec(ib0)=dbuf(ii)*zconv
               endif
            endif
         enddo

         exit
      enddo

      deallocate(dbuf)

    end subroutine nbget_alg1

    subroutine nbget_alg2(scnam,vecnam,psvec)

      ! get beam machine description data from splitn (TRANSP namelist)
      ! "algorithm 2" rule:  if 1st element of "vecnam" is 0, apply
      ! "scnam" value to all array elements

      ! zconv converts the data to MKS units for the plasma state

      character*(*), intent(in) :: scnam  ! TRANSP namelist scalar name
      character*(*), intent(in) :: vecnam ! TRANSP namelist vector name

      real*8, intent(out), dimension(:) :: psvec  ! plasma state data (output)

      !-------------------
      integer :: irank,idims(2,10)
      real*8, dimension(:), allocatable :: dbuf
      !-------------------

      ier = 0
      if(inum.eq.0) return

      idims=0
      call splitn_getdims(vecnam,irank,10,idims,isize,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) &
              ' ?trx_gen_state(nbget_alg2): namelist read failure: '// &
              'size of: '//trim(vecnam)
         return
      else if(irank.ne.1) then
         ier=99
         write(lunzer(0),*) &
              ' ?trx_gen_state(nbget_alg2): namelist vec unexpected rank: ', &
              irank,' for: '//trim(vecnam)
         return
      endif

      allocate(dbuf(isize))

      do
         call splitn_dget(vecnam,isize,dbuf,ier)
         if(ier.ne.0) then
            write(lunzer(0),*) &
                 ' ?trx_gen_state(nbget_alg2): namelist read failure: '// &
                 trim(vecnam)
            exit
         endif

         if(dbuf(1).eq.ZERO) then
            call splitn_dget(scnam,1,dbuf(1),ier)
            if(ier.ne.0) then
               write(lunzer(0),*) &
                    ' ?trx_gen_state(nbget_alg2): namelist read failure: '// &
                    trim(scnam)
               exit
            else
               dbuf(2:inum+intrace)=dbuf(1)
            endif
         endif

         do ib=1,inum+intrace
            ib0=imap(ib)
            if(ib0.gt.0) then
               psvec(ib0) = zconv*dbuf(ib)
            endif
         enddo
         exit
      enddo

      deallocate(dbuf)

    end subroutine nbget_alg2

    subroutine fracsort(zemin,ibmin,zelist,iblist,Evec,ffvec,fhvec)

      ! create table of beam fractions vs. beam energy, based on available
      ! data...

      real*8, intent(in) :: zemin   ! minimum energy
      integer, intent(in) :: ibmin  ! beam index w/ minimum energy
      real*8, dimension(:), intent(in) :: zelist  ! list of energies
      integer, dimension(:), intent(in) :: iblist ! list of indices

      real*8, dimension(:), intent(out) :: Evec   ! sorted energies out
      real*8, dimension(:), intent(out) :: ffvec  ! full energy fractions
      real*8, dimension(:), intent(out) :: fhvec  ! half energy fractions

      !----------------
      integer :: ii,ict,ibnext
      real*8 :: zenext
      !----------------

      Evec(1)=zemin
      ffvec(1)=ffulla(ibmin)
      fhvec(1)=fhalfa(ibmin)

      ict = 1

      do
         zenext=Evec(ict)
         do ii=1,size(zelist)
            if(zelist(ii).eq.ZERO) exit
            if(zelist(ii).le.Evec(ict)) cycle

            if(zenext.eq.Evec(ict)) then
               zenext=zelist(ii)
               ibnext=iblist(ii)

            else if(zelist(ii).lt.zenext) then
               zenext=zelist(ii)
               ibnext=iblist(ii)
            endif
         enddo
         if(zenext.eq.Evec(ict)) exit

         ict = ict + 1
         Evec(ict)=zenext
         ffvec(ict)=ffulla(ibnext)
         fhvec(ict)=fhalfa(ibnext)
      enddo

    end subroutine fracsort

    !==============================================================
    subroutine antgeo

      !---------------------------------------------------------
      ! get antenna geometry & n_phi information from TRANSP namelist
      ! the code here should match the handling of defaulted data
      ! items in the same way as TRANSP's trcore/datckich.for ...

      ! at time of call: inum = # of RF antennas...

      ! replicate TRANSP namelist arrays related to ICRF antenna geometry
      ! and n_phi spectrum...
      !---------------------------------------------------------

      integer, dimension(:), allocatable :: num_nphi,msym_nphi
      integer, dimension(:,:), allocatable :: nnphi
      real*8, dimension(:,:), allocatable :: wnphi

      real*8, dimension(:), allocatable :: sepicha,widicha
      real*8, dimension(:,:), allocatable :: phicha

      real*8, dimension(:), allocatable :: thicha,rmjicha,rmnicha

      integer :: ngeoant(1)
      real*8, dimension(:), allocatable :: rgeoant,ygeoant
      real*8 :: antrcen_0,antycen_0

      integer, dimension(:), allocatable :: ngeoant_a
      real*8, dimension(:,:), allocatable :: rgeoant_a,ygeoant_a
      real*8, dimension(:), allocatable :: antrcen_a,antycen_a

      integer :: idim1,idim2,iant,iminn,imaxn,icenn
      real*8 :: zdum,zsum,zdist,zrelph,zdegrad,zth0,zth1,zth
      !----------------------
      
      call splitn_iget('NGEOANT',1,ngeoant,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state: not in namelist: NGEOANT.'
         return
      endif

      call splitn_get1dim('RGEOANT',0,idim1,ier)
      if(ier.ne.0) return

      allocate(rgeoant(idim1),ygeoant(idim1))

      call splitn_dget1('RGEOANT',rgeoant,ier)
      if(ier.ne.0) return
      call splitn_dget1('YGEOANT',ygeoant,ier)
      if(ier.ne.0) return

      antrcen_0 = 0
      antycen_0 = 0

      !-----------

      call splitn_get1dim('THICHA',0,idim1,ier)
      if(ier.ne.0) return

      allocate(thicha(idim1),rmjicha(idim1),rmnicha(idim1))

      call splitn_dget1('THICHA',thicha,ier)
      if(ier.ne.0) return
      call splitn_dget1('RMJICHA',rmjicha,ier)
      if(ier.ne.0) return
      call splitn_dget1('RMNICHA',rmnicha,ier)
      if(ier.ne.0) return

      !-----------

      call splitn_get1dim('SEPICHA',inum,idim1,ier)
      if(ier.ne.0) return

      allocate(sepicha(idim1),widicha(idim1))

      call splitn_dget1('SEPICHA',sepicha,ier)
      if(ier.ne.0) return
      call splitn_dget1('WIDICHA',widicha,ier)
      if(ier.ne.0) return

      call splitn_get2dim('PHICHA',2,inum,idim1,idim2,ier)
      if(ier.ne.0) return

      allocate(phicha(idim1,idim2))

      call splitn_dget2('PHICHA',phicha,ier)
      if(ier.ne.0) return

      !-----------

      call splitn_get2dim('RGEOANT_A',15,inum,idim1,idim2,ier)
      if(ier.ne.0) return

      allocate(rgeoant_a(idim1,idim2),ygeoant_a(idim1,idim2))
      allocate(ngeoant_a(idim2)); ngeoant_a=0
      allocate(antrcen_a(idim2)); antrcen_a=0
      allocate(antycen_a(idim2)); antycen_a=0

      call splitn_dget2('RGEOANT_A',rgeoant_a,ier)
      if(ier.ne.0) return
      call splitn_dget2('YGEOANT_A',ygeoant_a,ier)
      if(ier.ne.0) return

      !-----------

      call splitn_get2dim('NNPHI',1,inum,idim1,idim2,ier)
      if(ier.ne.0) return

      allocate(num_nphi(idim2),msym_nphi(idim2))
      allocate(nnphi(idim1,idim2),wnphi(idim1,idim2))

      call splitn_iget1('NUM_NPHI',num_nphi,ier)
      if(ier.ne.0) return
      call splitn_iget1('MSYM_NPHI',msym_nphi,ier)
      if(ier.ne.0) return

      call splitn_iget2('NNPHI',nnphi,ier)
      if(ier.ne.0) return
      call splitn_dget2('WNPHI',wnphi,ier)
      if(ier.ne.0) return

      !----------------------------------------------------
      ! data acquired.  process...

      !--------------------------------
      ! per antenna geometry data:

      iminn=size(rgeoant_a,1) + 1
      imaxn=0

      do iant=1,inum
         j=0
         do
            if(rgeoant_a(j+1,iant).le.0.0_rspec) exit
            j=j+1
            if(j.eq.size(rgeoant_a,1)) exit   ! maximum index reached
         enddo

         ngeoant_a(iant)=j
         if(j.gt.0) call rfx_antcen(j,rgeoant_a(1:j,iant),ygeoant_a(1:j,iant),&
              antRcen_a(j),antYcen_a(j),zdum)
         iminn=min(iminn,j)
         imaxn=max(imaxn,j)
      enddo

      !--------------------------------
      ! general antenna geometry data: used if per antenna data is defaulted...

      if(ngeoant(1).eq.1) then
         ngeoant(1)=15
         zdegrad=TWOPI/360.0_rspec
         zth0=-zdegrad*thicha(1)/2.0_rspec
         zth1=-zth0
         do j=1,ngeoant(1)
            zth=zth0 + (j-1)*(zth1-zth0)/(ngeoant(1)-1)
            if(rmnicha(1).ge.0.0_rspec) then
               ! r > 0: R0 +r*cos(theta)
               rgeoant(j) = Rmjicha(1) + cos(zth)*rmnicha(1) ! Low Field Side
            else
               ! r < 0: vertical flat antenna
               rgeoant(j) = Rmjicha(1) + rmnicha(1)   ! High Field Side
            endif
            ygeoant(j) =   sin(zth)*rmnicha(1)
         enddo
         antrcen_0 = Rmjicha(1) + Rmnicha(1)
         antycen_0 = 0.0_rspec
      else if(ngeoant(1).gt.1) then
         call rfx_antcen(ngeoant(1),rgeoant(1:ngeoant(1)),ygeoant(1:ngeoant(1)), &
              antrcen_0,antycen_0,zdum)
      else
         if(iminn.eq.0) then
            write(lunzer(0),*) ' ?trx_gen_state: missing RF antenna data:'
            write(lunzer(0),*) '  If not all individual antenna geometries'
            write(lunzer(0),*) '  are specified, {rgeoant,ygeoant} must be.'
            ier=100
            return
         endif
      endif

      !--------------------------------
      ! copy general antenna data where needed

      do iant=1,inum
         if(ngeoant_a(iant).eq.0) then
            antrcen_a(iant)=antrcen_0
            antycen_a(iant)=antycen_0
            ngeoant_a(iant)=ngeoant(1)
            rgeoant_a(1:ngeoant(1),iant) = rgeoant(1:ngeoant(1))
            ygeoant_a(1:ngeoant(1),iant) = ygeoant(1:ngeoant(1))
         endif
      enddo

      !--------------------------------
      ! now deal with n_phi.  Note: in the TRANSP namelist, typically
      ! only the positive n_phi values are given and +/- symmetry is assumed.
      ! ...In the plasma state, the full spectrum will be stored.

      do iant=1,inum
         if(num_nphi(iant).eq.0) then
            num_nphi(iant)=1
            msym_nphi(iant)=1
            wnphi(1,iant)=1.0_rspec

            zrelph=phicha(2,iant)-phicha(1,iant)
            if(abs(zrelph) .gt. 90.0_rspec)then
               if(sepicha(iant).le.0.0_rspec) then
                  write(lunzer(0),*) ' ?trx_gen_state: missing RF antenna data.'
                  write(lunzer(0),*) '  iant=',iant,' num_nphi(iant)=0 but also'
                  write(lunzer(0),*) '  sepicha(iant)=',sepicha(iant)
                  ier=101
                  return
               endif
               nnphi(1,iant) = antRcen_a(iant)*TWOPI/(2*sepicha(iant)) + &
                    0.5_rspec
            else
               zdist=sqrt(widicha(iant)**2+sepicha(iant)**2)
               if(zdist.eq.0.0_rspec) then
                  write(lunzer(0),*) ' ?trx_gen_state: missing RF antenna data.'
                  write(lunzer(0),*) '  iant=',iant,' num_nphi(iant)=0 but also'
                  write(lunzer(0),*) '  sepicha(iant)=widicha(iant)=0.0'
                  ier=102
                  return
               else
                  nnphi(1,i)=antRcen_a(iant)*TWOPI/4.0_rspec/zdist + 0.5_rspec
               endif
            endif
         endif
      enddo

      !----------------------------------------------------
      ! OK, copy data into plasma state; de-symmetrize n_phi spectra

      ss%nrz_antgeo(1:inum) = ngeoant_a(1:inum)

      do iant=1,inum
         if(msym_nphi(iant).eq.1) then
            ss%num_nphi(iant) = 2*num_nphi(iant)
         else
            ss%num_nphi(iant) = num_nphi(iant)
         endif
      enddo

      call ps_alloc_plasma_state(ier, state=ss)
      if(ier.ne.0) return

      do iant=1,inum
         j=ngeoant_a(iant)
         ss%R_antgeo(1:j,iant) = 0.01_rspec*Rgeoant_a(1:j,iant)
         ss%Z_antgeo(1:j,iant) = 0.01_rspec*Ygeoant_a(1:j,iant)

         icenn=num_nphi(iant)
         if(msym_nphi(iant).eq.1) then
            do j=1,num_nphi(iant)
               ss%nphi(icenn+1-j,iant) = -nnphi(j,iant)
               ss%nphi(icenn+j,iant) = nnphi(j,iant)
               ss%wt_nphi(icenn+1-j,iant) = 0.5_rspec*wnphi(j,iant)
               ss%wt_nphi(icenn+j,iant) = 0.5_rspec*wnphi(j,iant)
            enddo
         else
            j=num_nphi(iant)
            ss%nphi(1:j,iant) = nnphi(1:j,iant)
            ss%wt_nphi(1:j,iant) = wnphi(1:j,iant)
         endif

         zsum=0.0_rspec
         do j=1,ss%num_nphi(iant)
            if(ss%wt_nphi(j,iant).le.0.0_rspec) then
               write(lunzer(0),*) ' ?trx_gen_state: weight <= 0 for ICRF n_phi'
               write(lunzer(0),*) '  iant=',iant,' #=',j,' n_phi=', &
                    ss%nphi(j,iant)
               ier=ier+1
            else
               zsum = zsum + ss%wt_nphi(j,iant)
            endif
            if(ier.ne.0) then
               ier=103
               exit
            endif
         enddo
         if(ier.ne.0) return

         j=ss%num_nphi(iant)
         ss%wt_nphi(1:j,iant) = ss%wt_nphi(1:j,iant)/zsum

      enddo

    end subroutine antgeo

    subroutine mk_ngsc0

      ! set up neutral gas sources
      !    two for each gas specie, 1 for "recycling" and 1 for "gas flow"

      character*3 ZLi  ! available Li gas label
      character*2 ZLtest
      integer :: iu

      ss%ngsc0 = 2*ss%nspec_gas
      call ps_alloc_plasma_state(ier, state=ss)
      if(ier.ne.0) then
         write(lunzer(0),*) &
              ' ?trx_gen_state_geq: ps_alloc_plasma_state(mk_ngsc0) status: ',&
              ier
         return
      endif

      in0=ss%nspec_gas
      allocate(cgas_abbrev(in0))

      ZLi='Li'
      do jig=1,size(ss%sgas_name)
         iu=index(ss%sgas_name(jig),'_')
         if(iu.le.0) cycle
         if(iu.gt.4) cycle
         ZLtest = ss%sgas_name(jig)(1:2)
         call uupper(ZLtest)
         if(ZLtest.eq.'LI') then
            ZLi =ss%sgas_name(jig)(1:iu-1)
         endif
      enddo

      do jig=1,in0
         iz=ss%qatom_sgas(jig)/ps_xe + 0.1d0
         ia=ss%m_sgas(jig)/ps_mp + 0.1d0

         ss%is_recycling(jig) = 1
         ss%is_recycling(jig+in0) = 0

         if(iz.eq.1) then
            if(ia.eq.1) then
               ss%gas_atom(jig)='H'
               ss%gas_atom(jig+in0)='H'
               ss%gs_name(jig)='H0rcy'
               ss%gs_name(jig+in0)='H0gf'
               cgas_abbrev(jig)='H'
            else if(ia.eq.2) then
               ss%gas_atom(jig)='D'
               ss%gas_atom(jig+in0)='D'
               ss%gs_name(jig)='D0rcy'
               ss%gs_name(jig+in0)='D0gf'
               cgas_abbrev(jig)='D'
            else if(ia.eq.3) then
               ss%gas_atom(jig)='T'
               ss%gas_atom(jig+in0)='T'
               ss%gs_name(jig)='T0rcy'
               ss%gs_name(jig+in0)='T0gf'
               cgas_abbrev(jig)='T'
            else
               ier=1
            endif
         else if(iz.eq.2) then
            if(ia.eq.3) then
               ss%gas_atom(jig)='HE3'
               ss%gas_atom(jig+in0)='HE3'
               ss%gs_name(jig)='HE3_0rcy'
               ss%gs_name(jig+in0)='HE3_0gf'
               cgas_abbrev(jig)='3'
            else if(ia.eq.4) then
               ss%gas_atom(jig)='HE4'
               ss%gas_atom(jig+in0)='HE4'
               ss%gs_name(jig)='HE4_0rcy'
               ss%gs_name(jig+in0)='HE4_0gf'
               cgas_abbrev(jig)='4'
            else
               ier=1
            endif
         else if(iz.eq.3) then
            ss%gas_atom(jig)=ZLi
            ss%gas_atom(jig+in0)=ZLi
            ss%gs_name(jig)=trim(ZLi)//'_rcy'
            ss%gs_name(jig+in0)=trim(ZLi)//'_gf'
            cgas_abbrev(jig)='L'
         else
            ier=1
         endif

         if(ier.ne.0) then
            write(lunzer(0),*) ' ?trx_gen_state_geq: ia=',ia,' iz=',iz
            write(lunzer(0),*) '  A & Z not of a supported neutral specie.'
            exit
         endif
      enddo
      if(ier.ne.0) return

      call ps_gsc0_species_map(ier, state=ss)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state: ps_gsc0_species_map error!'
         return
      endif

    end subroutine mk_ngsc0

    subroutine sinb_select(ibeamx,ifusx,ibeami,ifusi)

      ! only partial data available on neutral sinks-- namely, the total 
      ! amount for each neutral species, summed over all fast species.
      ! So, choose first eligible fast species and attribute to it alone
      ! the total rate.  The index of the chosen species, for two reactions,
      ! are output.

      ! neutralizing charge exchange sink: Z0 >= Zi
      integer, intent(out) :: ibeamx  ! chosen beam species, or zero
      integer, intent(out) :: ifusx   ! chosen fusion product species, or zero

      ! impact ionization + non-neutralizing CX sink
      integer, intent(out) :: ibeami  ! chosen beam species, or zero
      integer, intent(out) :: ifusi   ! chosen fusion product species, or zero

      ! main routine variables already set on entry to this contained routine:
      !   iz -- atomic number of neutral species
      !   n_bi -- number of beam ions
      !   n_fusi -- number of fusion ions
      !   at least one of {n_bi,n_fusi} is greater than zero.

      integer :: izf,iif

      write(lunzer(0),*) ' %trx_gen_state: reconstructing neutral sinks from partial data (OK)'

      ibeamx=0
      ifusx=0

      ibeami=0
      ifusi=0

      do iif=1,n_bi
         izf = ss%qatom_snbi(iif)/ps_xe + 0.1d0
         if(izf.le.iz) then
            ibeamx=iif
            exit
         endif
      enddo

      if(ibeamx.eq.0) then
         do iif=1,n_fusi
            izf = ss%qatom_sfus(iif)/ps_xe + 0.1d0
            if(izf.le.iz) then
               ifusx=iif
               exit
            endif
         enddo
      endif

      if(n_bi.gt.0) then
         ibeami=1
      else
         ifusi=1
      endif
      
    end subroutine sinb_select
      
    subroutine splitn_get1dim(znam,imin,idim,ier)
      !  get the dimension of a 1d array:
      !    report error of named item is not 1d or if idim < imin

      character*(*), intent(in) :: znam
      integer, intent(in) :: imin
      integer, intent(out) :: idim
      integer, intent(out) :: ier

      !-------------------------------
      integer :: irank,idims(2,10),isize
      !-------------------------------

      idim=0
      ier=0

      call splitn_getdims(znam,irank,10,idims,isize,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: error finding "', &
              trim(znam),'" in namelist.'
         return
      endif

      if(irank.ne.1) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: namelist item not 1d: ', &
              trim(znam)
         ier=1
         return
      endif

      idim=isize
      if(idim.lt.imin) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: "',trim(znam), &
              '" 1d item dimension too small:'
         write(lunzer(0),*) '  needed: ',imin,'  got: ',idim
         ier=2
         return
      endif

    end subroutine splitn_get1dim

    subroutine splitn_get2dim(znam,imin1,imin2,idim1,idim2,ier)
      !  get the dimension of a 1d array:
      !    report error of named item is not 1d or if idim < imin

      character*(*), intent(in) :: znam
      integer, intent(in) :: imin1,imin2
      integer, intent(out) :: idim1,idim2
      integer, intent(out) :: ier

      !-------------------------------
      integer :: irank,idims(2,10),isize
      !-------------------------------

      idim1=0
      idim2=0
      ier=0

      call splitn_getdims(znam,irank,10,idims,isize,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: error finding "', &
              trim(znam),'" in namelist.'
         return
      endif

      if(irank.ne.2) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: namelist item not 2d: ', &
              trim(znam)
         ier=1
         return
      endif

      idim1=idims(2,1)-idims(1,1)+1
      idim2=idims(2,2)-idims(1,2)+1

      if((idim1.lt.imin1).or.(idim2.lt.imin2)) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: "',trim(znam), &
              '" 2d item dimension too small:'
         write(lunzer(0),*) '  1st dim needed: ',imin1,'  got: ',idim1
         write(lunzer(0),*) '  2nd dim needed: ',imin2,'  got: ',idim2
         ier=2
         return
      endif

    end subroutine splitn_get2dim

    subroutine splitn_dget1(znam,zarray,ier)
      ! get real*8 splitn item

      character*(*), intent(in) :: znam
      real*8, dimension(:), intent(out) :: zarray
      integer, intent(out) :: ier

      !------------------------

      ier=0
      zarray = 0.0

      call splitn_dget(znam,size(zarray),zarray,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: could not get: ', &
              trim(znam)
      endif
    end subroutine splitn_dget1

    subroutine splitn_iget1(znam,iarray,ier)
      ! get real*8 splitn item

      character*(*), intent(in) :: znam
      integer, dimension(:), intent(out) :: iarray
      integer, intent(out) :: ier

      !------------------------

      ier=0
      iarray = 0.0

      call splitn_iget(znam,size(iarray),iarray,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: could not get: ', &
              trim(znam)
      endif
    end subroutine splitn_iget1

    subroutine splitn_dget2(znam,zarray,ier)
      ! get real*8 splitn item

      character*(*), intent(in) :: znam
      real*8, dimension(:,:), intent(out) :: zarray
      integer, intent(out) :: ier

      !------------------------
      real*8, dimension(:), allocatable :: zbuf
      integer :: isize,idim1,idim2,indx,ii,jj
      !------------------------

      ier=0
      zarray = 0.0

      idim1=size(zarray,1)
      idim2=size(zarray,2)

      isize = idim1*idim2
      allocate(zbuf(isize))

      call splitn_dget(znam,isize,zbuf,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: could not get: ', &
              trim(znam)
      endif

      indx=0
      do jj=1,idim2
         do ii=1,idim1
            indx = indx + 1
            zarray(ii,jj) = zbuf(indx)
         enddo
      enddo

      deallocate(zbuf)

    end subroutine splitn_dget2

    subroutine splitn_iget2(znam,iarray,ier)
      ! get real*8 splitn item

      character*(*), intent(in) :: znam
      integer, dimension(:,:), intent(out) :: iarray
      integer, intent(out) :: ier

      !------------------------
      integer, dimension(:), allocatable :: ibuf
      integer :: isize,idim1,idim2,indx,ii,jj
      !------------------------

      ier=0
      iarray = 0.0

      idim1=size(iarray,1)
      idim2=size(iarray,2)

      isize = idim1*idim2
      allocate(ibuf(isize))

      call splitn_iget(znam,isize,ibuf,ier)
      if(ier.ne.0) then
         write(lunzer(0),*) ' ?trx_gen_state_geq: could not get: ', &
              trim(znam)
      endif

      indx=0
      do jj=1,idim2
         do ii=1,idim1
            indx = indx + 1
            iarray(ii,jj) = ibuf(indx)
         enddo
      enddo

      deallocate(ibuf)

    end subroutine splitn_iget2

    subroutine get_coupled_spectrum

      ! find coupled spectrum data if available (caller set iexist logical)
      integer :: inum,iant,ierr
      integer :: iph,inph

      inum = ss%nicrf_src
      do iant=1,inum
         zname=' '
         if(iant.lt.10) then
            write(zname,'("PICHA",i1)') iant
         else
            write(zname,'("PICHA",i2)') iant
         endif
         call trx_scal(zname,zunits,ss%picrf_abs(iant),ier)

         if(allocated(ss%wt_nphi_abs)) then
            ! use vacuum spectrum data for now
            ss%wt_nphi_abs(:,iant) = ss%wt_nphi(:,iant)
        
            ! look for coupled spectrum data
           
            if(iexist) then
               if(iant.lt.10) then
                  write(zname,'("CPLSPEC",i1)') iant
               else
                  write(zname,'("CPLSPEC",i2)') iant
               endif

               call r8_t1profil(zname,zlabel,zunits,ztime0,zdelta_t, &
                    idum,zspectrum,igot,jgot,ierr)

               if(ierr.eq.0) then
                  do iph=1,ss%num_nphi(iant)
                     inph=ss%nphi(iph,iant)
                     call find_nphi_val(inph,zmbuf(1:igot),zspectrum, &
                          ss%wt_nphi_abs(iph,iant))
                  enddo
               endif

            endif  ! xgrid_nphi exists
         endif  ! wt_nphi allocated
      enddo  ! antenna loop

    end subroutine get_coupled_spectrum

    subroutine find_nphi_val(inph,xnphi,wspec,wgt1)

      integer, intent(in) :: inph   ! integer Nphi value
      real*8, dimension(:), intent(in) :: xnphi  ! real*8 "Nphi" grid
      real*8, dimension(:), intent(in) :: wspec  ! spectrum over xnphi grid

      real*8, intent(out) :: wgt1   ! spectrum value @inph (nearest value)

      !--------------------------
      integer :: ii,isize
      real*8 ztest_diff,zmin_diff,znph
      !--------------------------

      isize=size(xnphi)

      wgt1 = ZERO
      zmin_diff = 1.0d20

      znph = inph

      do ii=1,isize
         ztest_diff = abs(znph-xnphi(ii))
         if(ztest_diff.lt.zmin_diff) then
            zmin_diff = ztest_diff
            wgt1 = wspec(ii)
         endif
      enddo

    end subroutine find_nphi_val

    subroutine mk_anom(ierr)

      ! gather fast ion anomalous transport data; fill ANOM component in PS
      ! this routine is only called if both namelist & trdatbuf are available

      integer, intent(out) :: ierr

      !----------------------------------
      real*8 :: tmul_h_beam_ions(1)
      real*8 :: tmul_d_beam_ions(1)
      real*8 :: tmul_t_beam_ions(1)
      real*8 :: tmul_he3_beam_ions(1)
      real*8 :: tmul_he4_beam_ions(1)

      real*8 :: tmul_h_fusion_ions(1)
      real*8 :: tmul_t_fusion_ions(1)
      real*8 :: tmul_he3_fusion_ions(1)
      real*8 :: tmul_he4_fusion_ions(1)

      real*8, dimension(:), allocatable :: zmul_fi,zmul_bi

      logical :: nlqlim0(1)
      real*8 :: danom_qlim0(1)
      real*8 :: x_qlim0
      real*8 :: qlim0(1)
      real*8 :: dq_plim0(1)

      integer :: ndifbe,ndifbe2,nkdifb(1),nrho2(1),inume2,inx2,isize,ie,ix,ii
      integer :: idume,idumx,imdifb(1)

      logical :: tdb_xyprof_present
      logical :: have_d1d,have_d2d,have_e2d
      logical :: have_fd0,have_fdb,have_fdp,have_fdq,have_fdr,have_fds

      real*8, dimension(:), allocatable :: zdif_qlim0
      real*8, dimension(:), allocatable :: zedifb,zfdifb
      real*8, dimension(:), allocatable :: zx2d,ze2d
      real*8, dimension(:,:), allocatable :: zd2d
      !----------------------------------
      ierr=0

      !  gather information...

      !  availability of 1d transport data...
      call rpexist_profile('DIFB',have_d1d)

      !  multipliers for 2d data...
      call splitn_dget('tmul_h_beam_ions',1,tmul_h_beam_ions,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('tmul_d_beam_ions',1,tmul_d_beam_ions,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('tmul_t_beam_ions',1,tmul_t_beam_ions,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('tmul_he3_beam_ions',1,tmul_he3_beam_ions,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('tmul_he4_beam_ions',1,tmul_he4_beam_ions,ierr)
      if(ierr.ne.0) go to 999

      call splitn_dget('tmul_h_fusion_ions',1,tmul_h_fusion_ions,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('tmul_t_fusion_ions',1,tmul_t_fusion_ions,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('tmul_he3_fusion_ions',1,tmul_he3_fusion_ions,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('tmul_he4_fusion_ions',1,tmul_he4_fusion_ions,ierr)
      if(ierr.ne.0) go to 999

      if(n_bi.gt.0) then
         allocate(zmul_bi(n_bi))
         zmul_bi=0.0
         do ii=1,n_bi
            if(iznbi(ii).eq.1) then
               if(ianbi(ii).eq.1) then
                  zmul_bi(ii)=tmul_h_beam_ions(1)
               else if(ianbi(ii).eq.2) then
                  zmul_bi(ii)=tmul_d_beam_ions(1)
               else if(ianbi(ii).eq.3) then
                  zmul_bi(ii)=tmul_t_beam_ions(1)
               endif
            else if(iznbi(ii).eq.2) then
               if(ianbi(ii).eq.3) then
                  zmul_bi(ii)=tmul_he3_beam_ions(1)
               else if(ianbi(ii).eq.4) then
                  zmul_bi(ii)=tmul_he4_beam_ions(1)
               endif
            endif
         enddo
      endif

      if(n_fusi.gt.0) then
         allocate(zmul_fi(n_fusi))
         zmul_fi=0.0
         do ii=1,n_fusi
            if(iznfusi(ii).eq.1) then
               if(ianfusi(ii).eq.1) then
                  zmul_fi(ii)=tmul_h_fusion_ions(1)
               else if(ianfusi(ii).eq.3) then
                  zmul_fi(ii)=tmul_t_fusion_ions(1)
               endif
            else if(iznfusi(ii).eq.2) then
               if(ianfusi(ii).eq.3) then
                  zmul_fi(ii)=tmul_he3_fusion_ions(1)
               else if(ianfusi(ii).eq.4) then
                  zmul_fi(ii)=tmul_he4_fusion_ions(1)
               endif
            endif
         enddo
      endif

      !  application control for 1d data...
      call splitn_iget('nkdifb',1,nkdifb,ierr)
      if(ierr.ne.0) go to 999

      !  possible QLIM0 diffusivity enhancement...

      allocate(zdif_qlim0(inx)); zdif_qlim0=ZERO

      danom_qlim0(1)=0.0
      dq_plim0(1)=0.0
      x_qlim0=0.0

      call splitn_lget('nlqlim0',1,nlqlim0,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('danom_qlim0',1,danom_qlim0,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('dq_plim0',1,dq_plim0,ierr)
      if(ierr.ne.0) go to 999
      call splitn_dget('qlim0',1,qlim0,ierr)
      if(ierr.ne.0) go to 999

      if(nlqlim0(1)) then
         call trx_scal('X_QLIM0',zunits,x_qlim0,ierr)
         if(ierr.ne.0) then
            write(lunzer(0),*) ' ?trx_gen_state(mk_anom): X_QLIM0 not found.'
            go to 999
         endif

         danom_qlim0(1)=1.0d-4*danom_qlim0(1)  ! -> m^2/s

         if(x_qlim0.gt.ZERO) then
            do ix=1,inx
               if(ss%rho(ix).le.x_qlim0) then
                  zdif_qlim0(ix)=danom_qlim0(1)
               else if(zqsave(ix).ge.qlim0(1)+dq_plim0(1)) then
                  zdif_qlim0(ix)=danom_qlim0(1)
               else
                  exit  ! end of enhanced diffusivity...
               endif
            enddo
         endif

      endif

      !  energy scaling of 1d data (optional)

      call splitn_getsize('edifbe',isize,ierr)
      if(ierr.ne.0) go to 999

      allocate(zedifb(isize+1),zfdifb(isize+1)); zedifb=0.0; zfdifb=0.0
      call splitn_dget('edifbe',isize,zedifb(2:),ierr)
      if(ierr.ne.0) go to 999

      ndifbe=0
      do ie=2,isize+1
         if(zedifb(ie).ge.1.0d30) then
            exit
         endif
         ndifbe=ndifbe+1
      enddo
      if(ndifbe.gt.0) then
         ndifbe=ndifbe+1
         call splitn_dget('fdifbe',isize,zfdifb(2:),ierr)
         zfdifb(1)=zfdifb(2)
      endif

      ! 2d Diffusivity vs. (x,E) -- beam ions and fusion products only, for now
      !   nmdifb=4 must be set to activate this feature...

      nrho2(1)=0
      if((n_bi+n_fusi).gt.0) then
         call splitn_iget('nmdifb',1,imdifb,ierr)
         if(ierr.ne.0) go to 999
         if(imdifb(1).eq.4) then
            ! Db(x,E) indicated...
            call splitn_iget('nzone_nb',1,nrho2,ierr)
            if(ierr.ne.0) go to 999
            if(nrho2(1).eq.1) then
               nrho2(1)=0
            else
               nrho2(1)=nrho2(1) + 1  ! (want zone bdys here)
            endif
         endif
      endif

      inume2=0
      have_d2d=.FALSE.

      if((nrho2(1).gt.0).AND.tdb_xyprof_present(d,'FD0')) then
         have_fd0=.TRUE.
         call tdb_xysizes(d,'FD0','Ex',inume2,inx2,ierr)
         if(ierr.ne.0) go to 999
         call ckd2d(have_d2d,inume2,inx2,ndifbe2,ierr)
         if(ierr.ne.0) go to 999
      else
         have_fd0=.FALSE.
      endif

      if((nrho2(1).gt.0).AND.tdb_xyprof_present(d,'FDB')) then
         have_fdb=.TRUE.
         call tdb_xysizes(d,'FDB','Ex',inume2,inx2,ierr)
         if(ierr.ne.0) go to 999
         call ckd2d(have_d2d,inume2,inx2,ndifbe2,ierr)
         if(ierr.ne.0) go to 999
      else
         have_fdb=.FALSE.
      endif

      if((nrho2(1).gt.0).AND.tdb_xyprof_present(d,'FDP')) then
         have_fdp=.TRUE.
         call tdb_xysizes(d,'FDP','Ex',inume2,inx2,ierr)
         if(ierr.ne.0) go to 999
         call ckd2d(have_d2d,inume2,inx2,ndifbe2,ierr)
         if(ierr.ne.0) go to 999
      else
         have_fdp=.FALSE.
      endif

      if((nrho2(1).gt.0).AND.tdb_xyprof_present(d,'FDQ')) then
         have_fdq=.TRUE.
         call tdb_xysizes(d,'FDQ','Ex',inume2,inx2,ierr)
         if(ierr.ne.0) go to 999
         call ckd2d(have_d2d,inume2,inx2,ndifbe2,ierr)
         if(ierr.ne.0) go to 999
      else
         have_fdq=.FALSE.
      endif

      if((nrho2(1).gt.0).AND.tdb_xyprof_present(d,'FDR')) then
         have_fdr=.TRUE.
         call tdb_xysizes(d,'FDR','Ex',inume2,inx2,ierr)
         if(ierr.ne.0) go to 999
         call ckd2d(have_d2d,inume2,inx2,ndifbe2,ierr)
         if(ierr.ne.0) go to 999
      else
         have_fdr=.FALSE.
      endif

      if((nrho2(1).gt.0).AND.tdb_xyprof_present(d,'FDS')) then
         have_fds=.TRUE.
         call tdb_xysizes(d,'FDS','Ex',inume2,inx2,ierr)
         if(ierr.ne.0) go to 999
         call ckd2d(have_d2d,inume2,inx2,ndifbe2,ierr)
         if(ierr.ne.0) go to 999
      else
         have_fds=.FALSE.
      endif

      ! exit now, if there is no ANOM data

      if(.NOT.(have_d1d.or.have_d2d)) go to 1000

      ! OK -- ANOM Plasma State array dimensions...

      if(have_d1d) then
         ss%nrho_anom = inx
         ss%nefi_anom = ndifbe
      else
         ss%nrho_anom = 0
         ss%nefi_anom = 0
      endif

      if(have_d2d) then
         ss%nrho_anom2=nrho2(1)
         ss%nefi_anom2=inume2
      else
         ss%nrho_anom2=0
         ss%nefi_anom2=0
      endif

      call ps_alloc_plasma_state(ierr, state=ss)
      if(ierr.ne.0) go to 999

      !------------------------------------
      ! OK load the data...

      ! 1d diffusivity

      if(have_d1d) then

         ss%rho_anom = ss%rho  ! same grid as main plasma

         call get1b('DIFB')

         if((nkdifb(1).eq.1).or.(nkdifb(1).eq.3)) then
            ss%difb_nbi = zprof
            ss%difb_rfmi = zprof
         else
            ss%difb_nbi = ZERO
            ss%difb_rfmi = ZERO
         endif

         if((nkdifb(1).eq.2).or.(nkdifb(1).eq.3)) then
            ss%difb_fusi = zprof
         else
            ss%difb_fusi = ZERO
         endif

         ss%difb_nbi = max(zdif_qlim0, ss%difb_nbi)
         ss%difb_rfmi = max(zdif_qlim0, ss%difb_rfmi)
         ss%difb_fusi = max(zdif_qlim0, ss%difb_fusi)

         ! 1d radial velocity
         call get1b('VELB')

         if((nkdifb(1).eq.1).or.(nkdifb(1).eq.3)) then
            ss%velb_nbi = zprof
            ss%velb_rfmi = zprof
         else
            ss%velb_nbi = ZERO
            ss%velb_rfmi = ZERO
         endif

         if((nkdifb(1).eq.2).or.(nkdifb(1).eq.3)) then
            ss%velb_fusi = zprof
         else
            ss%velb_fusi = ZERO
         endif

         ! energy dependent scaling of 1d profiles...
         if(ss%nefi_anom.gt.0) then
            ss%E_anom = 0.001d0*zedifb(1:ndifbe)
            ss%anom_evar = zfdifb(1:ndifbe)
         endif
      endif

      !------------------------------------
      ! 2d data

      if(have_d2d) then
         have_e2d = .FALSE.

         ! the radial grid...
         ss%rho_anom2(1)=ZERO
         ss%rho_anom2(nrho2(1))=ONE
         do ix=2,nrho2(1)-1
            ss%rho_anom2(ix) = ((nrho2(1)-ix)*ZERO + (ix-1)*ONE)/(nrho2(1)-1)
         enddo

         if(have_fd0) then
            call tdb_xysizes(d,'FD0','Ex',inume2,inx2,ierr)
            if(ierr.ne.0) go to 999

            allocate(ze2d(inume2),zx2d(inx2))

            call tdb_xygrids(d,'FD0','Ex', &
                 ze2d,inume2,idume, zx2d,inx2,idumx, ierr)
            if(ierr.ne.0) go to 999

            ze2d=0.001d0*ze2d ! -> keV

            if(.not.have_e2d) then
               have_e2d=.TRUE.
               ss%E_anom2=ze2d
            else
               if(maxval(abs(ss%E_anom2 - ze2d)).gt. &
                    1.0d-8*maxval(abs(ze2d))) then
                  ierr=1
                  write(lunzer(0),*) ' ?trx_gen_state(mk_anom): '
                  write(lunzer(0),*) &
                       '  {FD0,FDB,FDP,FDQ,FDR,FDS} E-grids mismatch.'
               endif
            endif

            allocate(zd2d(inume2,inx2))
            call tdb_xyprof(d,'FD0',ss%t0,ss%t1,inume2,inx2,zd2d,ierr)
            if(ierr.ne.0) go to 999

            do ii=1,n_bi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_bi(ii), &
                    ss%bidiff_dtrap(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            do ii=1,n_fusi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_fi(ii), &
                    ss%fidiff_dtrap(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            deallocate(ze2d,zx2d,zd2d)
         endif

         if(have_fdb) then
            call tdb_xysizes(d,'FDB','Ex',inume2,inx2,ierr)
            if(ierr.ne.0) go to 999

            allocate(ze2d(inume2),zx2d(inx2))

            call tdb_xygrids(d,'FDB','Ex', &
                 ze2d,inume2,idume, zx2d,inx2,idumx, ierr)
            if(ierr.ne.0) go to 999

            ze2d=0.001d0*ze2d ! -> keV

            if(.not.have_e2d) then
               have_e2d=.TRUE.
               ss%E_anom2=ze2d
            else
               if(maxval(abs(ss%E_anom2 - ze2d)).gt. &
                    1.0d-8*maxval(abs(ze2d))) then
                  ierr=1
                  write(lunzer(0),*) ' ?trx_gen_state(mk_anom): '
                  write(lunzer(0),*) &
                       '  {FD0,FDB,FDP,FDQ,FDR,FDS} E-grids mismatch.'
               endif
            endif

            allocate(zd2d(inume2,inx2))
            call tdb_xyprof(d,'FDB',ss%t0,ss%t1,inume2,inx2,zd2d,ierr)
            if(ierr.ne.0) go to 999

            do ii=1,n_bi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_bi(ii), &
                    ss%bidiff_btrap(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            do ii=1,n_fusi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_fi(ii), &
                    ss%fidiff_btrap(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            deallocate(ze2d,zx2d,zd2d)
         endif

         if(have_fdp) then
            call tdb_xysizes(d,'FDP','Ex',inume2,inx2,ierr)
            if(ierr.ne.0) go to 999

            allocate(ze2d(inume2),zx2d(inx2))

            call tdb_xygrids(d,'FDP','Ex', &
                 ze2d,inume2,idume, zx2d,inx2,idumx, ierr)
            if(ierr.ne.0) go to 999

            ze2d=0.001d0*ze2d ! -> keV

            if(.not.have_e2d) then
               have_e2d=.TRUE.
               ss%E_anom2=ze2d
            else
               if(maxval(abs(ss%E_anom2 - ze2d)).gt. &
                    1.0d-8*maxval(abs(ze2d))) then
                  ierr=1
                  write(lunzer(0),*) ' ?trx_gen_state(mk_anom): '
                  write(lunzer(0),*) &
                       '  {FD0,FDB,FDP,FDQ,FDR,FDS} E-grids mismatch.'
               endif
            endif

            allocate(zd2d(inume2,inx2))
            call tdb_xyprof(d,'FDP',ss%t0,ss%t1,inume2,inx2,zd2d,ierr)
            if(ierr.ne.0) go to 999

            do ii=1,n_bi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_bi(ii), &
                    ss%bidiff_bpass_co(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            do ii=1,n_fusi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_fi(ii), &
                    ss%fidiff_bpass_co(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            deallocate(ze2d,zx2d,zd2d)
         endif

         if(have_fdq) then
            call tdb_xysizes(d,'FDQ','Ex',inume2,inx2,ierr)
            if(ierr.ne.0) go to 999

            allocate(ze2d(inume2),zx2d(inx2))

            call tdb_xygrids(d,'FDQ','Ex', &
                 ze2d,inume2,idume, zx2d,inx2,idumx, ierr)
            if(ierr.ne.0) go to 999

            ze2d=0.001d0*ze2d ! -> keV

            if(.not.have_e2d) then
               have_e2d=.TRUE.
               ss%E_anom2=ze2d
            else
               if(maxval(abs(ss%E_anom2 - ze2d)).gt. &
                    1.0d-8*maxval(abs(ze2d))) then
                  ierr=1
                  write(lunzer(0),*) ' ?trx_gen_state(mk_anom): '
                  write(lunzer(0),*) &
                       '  {FD0,FDB,FDP,FDQ,FDR,FDS} E-grids mismatch.'
               endif
            endif

            allocate(zd2d(inume2,inx2))
            call tdb_xyprof(d,'FDQ',ss%t0,ss%t1,inume2,inx2,zd2d,ierr)
            if(ierr.ne.0) go to 999

            do ii=1,n_bi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_bi(ii), &
                    ss%bidiff_dpass_co(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            do ii=1,n_fusi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_fi(ii), &
                    ss%fidiff_dpass_co(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            deallocate(ze2d,zx2d,zd2d)
         endif

         if(have_fdr) then
            call tdb_xysizes(d,'FDR','Ex',inume2,inx2,ierr)
            if(ierr.ne.0) go to 999

            allocate(ze2d(inume2),zx2d(inx2))

            call tdb_xygrids(d,'FDR','Ex', &
                 ze2d,inume2,idume, zx2d,inx2,idumx, ierr)
            if(ierr.ne.0) go to 999

            ze2d=0.001d0*ze2d ! -> keV

            if(.not.have_e2d) then
               have_e2d=.TRUE.
               ss%E_anom2=ze2d
            else
               if(maxval(abs(ss%E_anom2 - ze2d)).gt. &
                    1.0d-8*maxval(abs(ze2d))) then
                  ierr=1
                  write(lunzer(0),*) ' ?trx_gen_state(mk_anom): '
                  write(lunzer(0),*) &
                       '  {FD0,FDB,FDP,FDQ,FDR,FDS} E-grids mismatch.'
               endif
            endif

            allocate(zd2d(inume2,inx2))
            call tdb_xyprof(d,'FDR',ss%t0,ss%t1,inume2,inx2,zd2d,ierr)
            if(ierr.ne.0) go to 999

            do ii=1,n_bi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_bi(ii), &
                    ss%bidiff_bpass_ctr(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            do ii=1,n_fusi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_fi(ii), &
                    ss%fidiff_bpass_ctr(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            deallocate(ze2d,zx2d,zd2d)

         else
            ! copy fdp -> fdr

            do ii=1,n_bi
               ss%bidiff_bpass_ctr(:,:,ii) = ss%bidiff_bpass_co(:,:,ii)
            enddo

            do ii=1,n_fusi
               ss%fidiff_bpass_ctr(:,:,ii) = ss%fidiff_bpass_co(:,:,ii)
            enddo

         endif

         if(have_fds) then
            call tdb_xysizes(d,'FDS','Ex',inume2,inx2,ierr)
            if(ierr.ne.0) go to 999

            allocate(ze2d(inume2),zx2d(inx2))

            call tdb_xygrids(d,'FDS','Ex', &
                 ze2d,inume2,idume, zx2d,inx2,idumx, ierr)
            if(ierr.ne.0) go to 999

            ze2d=0.001d0*ze2d ! -> keV

            if(.not.have_e2d) then
               have_e2d=.TRUE.
               ss%E_anom2=ze2d
            else
               if(maxval(abs(ss%E_anom2 - ze2d)).gt. &
                    1.0d-8*maxval(abs(ze2d))) then
                  ierr=1
                  write(lunzer(0),*) ' ?trx_gen_state(mk_anom): '
                  write(lunzer(0),*) &
                       '  {FD0,FDB,FDP,FDQ,FDR,FDS} E-grids mismatch.'
               endif
            endif

            allocate(zd2d(inume2,inx2))
            call tdb_xyprof(d,'FDS',ss%t0,ss%t1,inume2,inx2,zd2d,ierr)
            if(ierr.ne.0) go to 999

            do ii=1,n_bi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_bi(ii), &
                    ss%bidiff_dpass_ctr(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            do ii=1,n_fusi
               call ps_user_1dintrp(ss%rho_anom2,zx2d,zd2d*zmul_fi(ii), &
                    ss%fidiff_dpass_ctr(:,:,ii),ierr, &
                    iswap=.TRUE.,iswapa=.TRUE.)
               if(ierr.ne.0) go to 999
            enddo

            deallocate(ze2d,zx2d,zd2d)

         else
            ! copy fdq -> fds

            do ii=1,n_bi
               ss%bidiff_dpass_ctr(:,:,ii) = ss%bidiff_dpass_co(:,:,ii)
            enddo

            do ii=1,n_fusi
               ss%fidiff_dpass_ctr(:,:,ii) = ss%fidiff_dpass_co(:,:,ii)
            enddo

         endif

      endif
      ierr=0
      go to 1000

999   continue
      write(lunzer(0),*) &
           ' ?trx_gen_state(mk_anom): error acquiring fast ion transport data.'
      ierr=1

1000  continue
      if(allocated(zdif_qlim0)) deallocate(zdif_qlim0)
      if(allocated(zedifb)) deallocate(zedifb,zfdifb)
      if(allocated(ze2d)) deallocate(ze2d,zx2d)
      if(allocated(zd2d)) deallocate(zd2d)

      return

    end subroutine mk_anom

    subroutine ckd2d(have_d2d,inume2,inx2,ndifbe2,ierr)
      ! set grid sizes, or, check that grid sizes match

      logical, intent(inout) :: have_d2d  ! .FALSE. if this is first one
      integer, intent(in) :: inume2,inx2  ! grid sizes (in)
      integer, intent(inout) :: ndifbe2   ! energy grid size (saved)
      integer, intent(out) :: ierr        ! exit status, 0=OK

      ierr=0

      if(.not.have_d2d) then
         have_d2d=.TRUE.
         ndifbe2=inume2
      else
         if(ndifbe2.ne.inume2) then
            write(lunzer(0),*) &
                 ' ?trx_gen_state(mk_anom): energy grid size mismatch: ', &
                 ndifbe2,inume2
            ierr=ierr+1
         endif
      endif

      if(inume2.lt.2) then
         write(lunzer(0),*) &
                 ' ?trx_gen_state(mk_anom): energy grid size too small: ', &
                 ndifbe2,inume2
         ierr=ierr + 1
      endif

      if(inx2.lt.2) then
         write(lunzer(0),*) &
                 ' ?trx_gen_state(mk_anom): x grid size too small: ', &
                 inx2
         ierr=ierr + 1
      endif

    end subroutine ckd2d

end subroutine trx_gen_state_geq

subroutine trx_init_state(ier)

  !  create an empty PLASMA STATE (ala SWIM Fusion Simulation Project)
  !  only a label is inserted. -- use "ps" instance from plasma_state_mod

  use plasma_state_mod
  use trx_module
  implicit NONE

  !------------------------  
  integer, intent(out) :: ier
  !------------------------  

  call trx_init_state_obj(ps,ier)

end subroutine trx_init_state

subroutine trx_init_state_obj(ss,ier)

  !  create an empty PLASMA STATE (ala SWIM Fusion Simulation Project)
  !  only a label is inserted. -- use "ps" instance from plasma_state_mod

  !  idea is to have trx_gen_state processing start with an empty state -- 
  !  no trace of possible prior contents

  !  mod DMC -- look for trxplib specified machine description files
  !       if found load them into aux

  use plasma_state_mod
  use trx_module
  implicit NONE

  !------------------------  
  type (plasma_state) :: ss
  integer, intent(out) :: ier
  integer :: lunzer

  integer :: inum_mdescr,idescr
  character*150 :: mdescr
  !------------------------  

  call ps_init_tag
  call plasma_state_erase(ss,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state_geq: plasma_state_erase status: ',ier
     return
  endif
  
  call plasma_state_erase(aux,ier)
  if(ier.ne.0) then
     write(lunzer(0),*) ' ?trx_gen_state_geq: "aux" erasure status: ',ier
     write(lunzer(0),*) '  error setting up for read of machine description.'
     return
  endif

  call trxplib_getnum_mdescr(inum_mdescr)
  if(inum_mdescr.gt.0) then
     
     do idescr=1,inum_mdescr
        call trxplib_getfull_mdescr(idescr,mdescr,ier)
        if(ier.ne.0) then
           call trxplib_getname_mdescr(idescr,mdescr)
           write(lunzer(0),*) ' ?trx_gen_state_geq: could not OPEN machine description file: '
           write(lunzer(0),*) '  '//trim(mdescr)
           exit
        endif

        call ps_mdescr_read(mdescr,ier, state=aux)
        if(ier.ne.0) then
           call trxplib_getname_mdescr(idescr,mdescr)
           write(lunzer(0),*) ' ?trx_gen_state_geq: could not READ machine description file: '
           write(lunzer(0),*) '  '//trim(mdescr)
           exit
        endif
     enddo
     if(ier.ne.0) return

  endif

  !----------------------------------------------------
  !  OK: label the new state...

  ss%global_label = trim(run_label)//' (trxplib)'
  ss%plasma_code_info = 'TRANSP'

end subroutine trx_init_state_obj
