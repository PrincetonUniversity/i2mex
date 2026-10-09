subroutine trx_gen_ech_files(ss, fpath,runid,ierr)

  ! (called from trx_gen_state) -- generate files for driving ECH codes
  ! when Plasma State derived data is available

  ! The plasma state (PS) has been loaded!
  ! SPLITN access to TRANSP namelist variables is also available.

  ! maintenance note: this routine generates TORAY input files so the 
  ! code has much in common with heatlib/toraymain.f90

  use plasma_state_mod
  use trx_module

  implicit NONE

  type (plasma_state) :: ss           ! state object to use...

  character*(*), intent(in) :: fpath  ! directory in which to write files
  character*(*), intent(in) :: runid  ! TRANSP runid
  integer, intent(out) :: ierr   ! completion status, 0=OK

  !  local...
  !-----------------------------------------------------
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  !-----------------------------------------------------
  real*8, parameter :: zminpow = 0.001_R8

  logical :: dirsw
  character*512 :: zdsave
  character*120 :: zlabel
  integer :: inuma  ! number of antennas or sources
  integer :: nonlin,lunzer,iant,isize,iray,iset,ig,jk,idamp,ilun
  integer :: istat,icta

  character*6 :: ichant
  character*6, dimension(:), allocatable :: ichaa ! ascii encoded antenna #s

  real*8 :: zfmu0,zrfmod,zx00,zz00,zthet,zphai,zbhalf,zbsratio

  character*1 :: bslash = '\\'  ! unix needs escaped backslash

  !-----------------------------------------------------
  !  TORAY/ECH namelist variables

  integer, dimension(:), allocatable :: ndampech,idampech,nrayech
  logical, dimension(:), allocatable :: iactive
  logical :: iswitch

  !-----------------------------------------------------
  !
  !         GEQDSK file related variables and arrays
  !
 
  character*80 z_geqroot,z_psiin
  
  integer inR,inZ,inB,inprof

  real*8 :: zgafsep   ! "GAFSEP" for TORAY: keep small, 1e-5/#zones
  ! (see comments in toray.f90)

  !-----------------------

  ierr = 0

  nonlin = lunzer(0)

  write(nonlin,*) ' '
  write(nonlin,*) ' %trx_gen_ech_files: look for TORAY and GENRAY inputs...'

  call loc_splitn_lget1('nltoray',iswitch)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?NLTORAY access error (ignored)'
     ierr=0
  else
     if(iswitch) then
        write(nonlin,*) ' ...the TRANSP run used TORAY'
     else
        write(nonlin,*) ' ...the TRANSP run did not use TORAY: TORAY controls may not be set correctly.'
     endif
  endif

  call loc_splitn_lget1('nlgen_ech',iswitch)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?NLGEN_ECH access error (ignored)'
     ierr=0
  else
     if(iswitch) then
        write(nonlin,*) ' ...the TRANSP run used GENRAY for ECH'
     else
        write(nonlin,*) ' ...the TRANSP run did not use GENRAY for ECH;'
        write(nonlin,*) '    default GENRAY template file will be sought.'
     endif
  endif

  call loc_splitn_lget1('nltorbeam',iswitch)
  if(ierr.ne.0) then
     write(nonlin,*) ' ?NLTORBEAM access error (ignored)'
     ierr=0
  else
     if(iswitch) then
        write(nonlin,*) ' ...the TRANSP run used TORBEAM for ECH'
     else
        write(nonlin,*) ' ...the TRANSP run did not use TORBEAM for ECH'
     endif
  endif

  call find_io_unit(ilun)

  ! # antennas; form ascii encoded list of antenna numbers

  inuma=ss%necrf_src
  if(inuma.le.0) then
     write(nonlin,*) ' %trx_gen_ech_files: no ECH sources found.'
     return
  endif

  allocate(ichaa(inuma),iactive(inuma))

  icta=0
  do iant=1,inuma
     write(ichant,'(i6)') iant
     ichaa(iant)=trim(adjustl(ichant))
     if(ss%power_ec(iant).gt.zminpow) then
        icta=icta+1
        write(nonlin,*) ' trx_gen_ech_files: antenna #',iant,' has power.'
        iactive(iant)=.TRUE.
     else
        iactive(iant)=.FALSE.
     endif
  enddo

  if(icta.le.0) then
     write(nonlin,*) ' %trx_gen_ech_files: no ECH antennas have power.'
     return
  endif

  !----------------------------------------------------
  if((fpath.ne.' ').AND.(fpath.ne.'.')) then
     ! cd to output directory
     dirsw=.TRUE.
     call getcwd(zdsave)
     call sset_cwd(fpath,ierr)
     if(ierr.ne.0) return
  else
     dirsw=.FALSE.
  endif

  do
     !----------------------------------------------------
     ! TORAY profile size

     call loc_splitn_iget1('NPROFTOR',inprof); if(ierr.ne.0) exit

     ! put rho_ecrf grid into ps

     ss%nrho_ecrf = inprof
     call ps_alloc_plasma_state(ierr, state=ss)
     if(ierr.ne.0) then
        write(nonlin,*) ' ?trx_gen_ech_files: rho_ecrf allocate error.'
        exit
     endif

     ss%rho_ecrf(1)=0.0_R8
     ss%rho_ecrf(inprof)=1.0_R8
     do jk=2,inprof-1
        ss%rho_ecrf(jk)=(jk-1)*1.0_R8/(inprof-1)
     enddo

     ! equilibrium data...

     zgafsep = 1.0d-5/inprof   ! "GAFSEP" for TORAY: keep small, 1e-5/#zones

     z_geqroot='eqdskin'
     z_psiin='psiin'

     inR=129
     inZ=129
     inB=129

     !  write G-eqdsk:

     write(nonlin,*) '   ...file: '//trim(z_geqroot)
     call ps_geqdsk_write_nRZ(ss,z_geqroot, &
          'trxpl:trx_gen_ech_files, TRANSP id: '//trim(runid)//' ', &
          inR,inZ,inB,nonlin,ierr)
     if(ierr.ne.0) then
        write(nonlin,*) '?trx_gen_ech_files: ps_geqdsk_write_nRZ error.'
        ierr=1
        exit
     endif

     !  read G-eqdsk; write psiin:

     write(nonlin,*) '   ...file: '//trim(z_psiin)
     call ps_psiin_write(ss,nonlin,z_geqroot,z_psiin,inprof,zgafsep, &
          ierr)
     if(ierr.ne.0) then
        write(nonlin,*) '?trx_gen_ech_files: ps_psiin_write error.'
        ierr=1
        exit
     endif

     !  write toray.in & gafit.in

     write(nonlin,*) ' '
     write(nonlin,*) ' ...acquire TRANSP namelist variables; write files: toray.in and gafit.in'

     zlabel = 'ECH and ECCD (trxpl:trx_gen_ech_files)'
     ierr=0
     call trx_toray_in_search(nonlin,ierr)
     if (ierr.ne.0) then
        write(nonlin,*) ' ?trx_gen_ech_files: ', &
             ' failed to copy toray.in'
        call errmsg_exit(' ?trx_gen_ech_files: no template toray.in')
     end if
     
     ! Copy GAFIT template namelist into toray.in
     call trx_gafit_in_search(nonlin,ierr)
     if (ierr.ne.0) then
        write(nonlin,*) ' ?trx_gen_ech_files: ', &
             ' copy_gafitnml: no template found'
        call errmsg_exit(' ?trx_gen_ech_files:  no template gafit.in')
     end if
      
     ! compute number of rays implied by mray(...) settings; this will
     ! override TRANSP namelist nrayech(...) control

     call splitn_getsize('NDAMPECH',isize,ierr)
     if(ierr.ne.0) then
        write(nonlin,*) ' %trx_gen_ech_files: namelist array size error.'
        ierr=1
        exit
     endif

     allocate(ndampech(isize),idampech(isize),nrayech(isize))

     call splitn_iget('NDAMPECH',isize,ndampech,ierr)
     if(ierr.ne.0) then
        write(nonlin,*) ' %trx_gen_ech_files: namelist array read error.'
        ierr=1
        exit
     endif

     call splitn_iget('NRAYECH',isize,nrayech,ierr)
     if(ierr.ne.0) then
        write(nonlin,*) ' %trx_gen_ech_files: namelist array read error.'
        ierr=1
        exit
     endif

     call echin_ck_dampmod(ss%necrf_src,ndampech,idampech,nonlin)

     ! loop to prepare antenna-specific input files

     do iant=1,ss%necrf_src

        if (iactive(iant)) then

           idamp=idampech(iant)
           zfmu0= ss%freq_ec(iant)
           zrfmod= ss%EC_Omode_fraction(iant) 
           zx00= 100.0_R8*ss%R_EC_launch(iant)
           zz00= 100.0_R8*ss%Z_EC_launch(iant)
           zthet=ss%EC_theta_aim(iant)
           zphai=ss%EC_phi_aim(iant)
           zbhalf=ss%EC_Half_Power_Angle(iant)
           zbsratio=ss%EC_Beam_Elongation(iant)
           
           write(nonlin,*) ' ... write: echin_ant'//trim(ichaa(iant))
           call ps2echin(ss,iant,idamp,nrayech(iant), nonlin, &
                'echin_ant'//ichaa(iant), ierr)
        else
           ! if the power is zero, do not leave echin file hanging around
           call fdelete('echin_ant'//ichaa(iant),istat)
        endif
     enddo
     ! end of pseudo-loop to write ECH model input files
     exit
  enddo

  !----------------------------------------------------
  !  write script to drive TORAY

  open(unit=ilun,file='toray_job.csh',status='unknown')
  write(nonlin,*) ' ...write toray_job.csh -- supply TORAY executable as argument...'

  write(ilun,1001)
1001 format('#! /bin/csh -f'/ &
          '#'/'# toray driver script '/ &
          '# $1 is the command to run GA TORAY'/'#'/ &
          '  set toray_command = $1')

  do iant=1,ss%necrf_src
     if(iactive(iant)) then
        write(ilun,1002) trim(ichaa(iant)),trim(ichaa(iant))
1002    format(/'  cp echin_ant',a,' echin'/'  $toray_command'/ &
             '  mv toray.nc toray.nc_ant',a)
     endif
  enddo

  close(unit=ilun)

  !----------------------------------------------------
  !  GENRAY...

  do
     ! look for genray.in template file

     write(nonlin,*) ' ...look for GENRAY template namelist... '
     call trx_genray_in_search(nonlin,ierr)
     if(ierr.ne.0) then
        ! this is kind of normal: only some tokamaks have GENRAY template
        ! data available...
        write(nonlin,*) ' %trx_gen_ech_files: GENRAY not available.'
        ierr=0
        exit
     endif

     exit
  enddo

  !----------------------------------------------------
  !  write script to drive GENRAY

  open(unit=ilun,file='genray_job.csh',status='unknown')
  write(nonlin,*) ' ...write genray_job.csh -- supply GENRAY executable as argument...'

  write(ilun,2001)
2001 format('#! /bin/csh -f'/ &
          '#'/'# genray driver script '/ &
          '# $1 is the command to run GENRAY'/'#'/ &
          '  set genray_command = $1'/ &
          '#'/'  prepare_genray_input "init" "EC" "1" "./genray.in" ', &
            '"disabled" "disabled" "yes" >& prep0.log '/ &
          '  if ( $status ) then'/ &
          '    echo "prepare_genray_input init error (see prep0.log)."'/ &
          '    exit 1'/ &
          '  endif'/'#')

  icta=0
  do iant=1,ss%necrf_src
     if(iactive(iant)) then
        icta=icta+1
        write(ilun,*) '######################'
        write(ilun,*) '#  antenna # '//trim(ichaa(iant))
        write(ilun,*) '######################'
        write(ilun,*) ' echo "antenna # '//trim(ichaa(iant))//'"'
        write(ilun,2002) trim(ichaa(iant)),trim(ichaa(iant)), &
             trim(ichaa(iant)),bslash,bslash, &
             trim(ichaa(iant)),trim(ichaa(iant)), &
             trim(ichaa(iant)),trim(ichaa(iant)),trim(ichaa(iant)), &
             trim(ichaa(iant))
2002    format(/'  cp genray_ech.template genray.in'/ &
             '  prepare_genray_input "step" "EC" "',a, &
                '" "./genray.in" "disabled" "disabled" "yes" >& prep',a,'.log'/ &
             '  if ( $status ) then'/ &
             '    echo "prepare_genray_input step error. (see prep',a,'.log)"'/ &
             '    exit 1'/ &
             '  endif'/ &
             '  #...dmc Oct 2010: merge_namelist step added.'/ &
             '  mv genray.in genray.fortran'/ &
             '  merge_namelist  template=genray_ech.template ',a/ &
             '     fortran=genray.fortran  output=genray.in ',a/ &
             '       >& merge',a,'.log'/ &
             '  if ( $status ) then'/ &
             '    echo "merge_namelist error. (see merge',a,'.log)"'/ &
             '    exit 1'/ &
             '  endif'/ &
             '  cp genray.in genray.in_ant',a// &
             '  $genray_command >& genray',a,'.log'/ &
             '  if ( $status ) then'/ &
             '    echo "genray step error (see genray',a,'.log)"'/ &
             '    exit 1'/ &
             '  endif'/ &
             '  cp genray.nc genray.nc_ant',a)

        if(icta.eq.1) then
           ichant = '-'//trim(ichaa(iant))
        else
           ichant = ichaa(iant)
        endif

        write(ilun,2003) trim(ichant)
2003    format(/'  process_genray_output "EC" "',a,'"'/ &
             '  if ( $status ) then'/ &
             '    echo "process_genray_output error."'/ &
             '    exit 1'/ &
             '  endif')
     endif
  enddo

  close(unit=ilun)

  !----------------------------------------------------
  if(dirsw) then
     ! cd back to original directory
     call sset_cwd(zdsave,ierr)
  endif

CONTAINS

  subroutine loc_splitn_iget1(zname,ival)
    character*(*), intent(in) :: zname
    integer, intent(out) :: ival
    integer, dimension(1) :: itemp

    itemp(1) = 0
    call splitn_iget(zname,1,itemp,ierr)
    ival = itemp(1)
    if(ierr.ne.0) then
       write(nonlin,*) ' %trx_gen_ech_files: namelist integer value error: ', &
            trim(zname)
       ierr=1
    endif
  end subroutine loc_splitn_iget1

  subroutine loc_splitn_dget1(zname,dval)
    character*(*), intent(in) :: zname
    real*8, intent(out) :: dval

    dval = 0.0d0
    call splitn_dget(zname,1,dval,ierr)
    if(ierr.ne.0) then
       write(nonlin,*) ' %trx_gen_ech_files: namelist real*8 value error: ', &
            trim(zname)
       ierr=1
    endif
  end subroutine loc_splitn_dget1

  subroutine loc_splitn_lget1(zname,lval)
    character*(*), intent(in) :: zname
    logical, intent(out) :: lval

    lval = .FALSE.
    call splitn_lget(zname,1,lval,ierr)
    if(ierr.ne.0) then
       write(nonlin,*) ' %trx_gen_ech_files: namelist logical value error: ', &
            trim(zname)
       ierr=1
    endif
  end subroutine loc_splitn_lget1

end subroutine trx_gen_ech_files
