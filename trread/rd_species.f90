subroutine rd_species(maxn,n_species,zlbla,abray,itype,ifast,zz,aa,izc)
!
!   find all plasma species -- return information in passed arrays
!   NOTE ordering on output lists:
!     electrons
!     then non-impurity thermal ions
!     them impurity thermal ions
!     then non-thermal ions
!
  implicit NONE
!
  integer, intent(in) :: maxn    ! max no. of species (array sizes)
  integer, intent(out) :: n_species  ! actual no. of species (-1 if error)
!
  character*20, intent(out) :: zlbla(maxn)      ! descr. label for each specie
!
  character*10, intent(out) :: abray(4,maxn)    ! specie profile names
!
!   for specie j:
!     abray(1,j) = name of density profile-- for all species (n/cm**3)
!     abray(2,j) = name of temperature or "average energy" profile.
!     abray(3,j) = name of perpendicular energy density profile or
!                  average perp energy per particle or " ".
!     abray(4,j) = name of parallel energy density profile or
!                  average pll energy per particle or " ".
!        (generally abray(3:4,j) will be set for fast species only)
!
  integer, intent(out) :: itype(maxn)           ! species type code
!
!   itype(j)=-1  -- electrons
!
!   itype(j)=+1  -- non-impurity thermal specie, usually H or He isotope
!                   can be Li
!   itype(j)=+2  -- impurity thermal specie:  Z and A are known, constant
!   itype(j)=+3  -- impurity thermal specie:  Z and A are functions of time
!
!   itype(j)=+4  -- beam ion specie
!   itype(j)=+5  -- rf tail ion specie
!   itype(j)=+6  -- fusion product ion specie
!
  integer, intent(out) :: ifast(maxn)           ! thermal/fast flag
!
!   ifast(j)=0   -- thermal electron or ion specie,
!                   abray(2,j) contains temperature, eV, abray(3:4,j)=" ".
!   ifast(j)=1   -- fast specie, abray(2,j) is avg energy / ptcl <E>, eV
!                   abray(3:4,j) contain perp, pll energy density, Jles/cm**3
!   ifast(j)=2   -- fast specie, abray(2,j) is "temperature" = 2/3 <E>, eV
!                   abray(3:4,j) contain perp, pll <E>/ptcl, eV
!
!        (other fast specie configurations may be added later)
!
  real, intent(out) :: zz(maxn)                 ! charge of specie
  real, intent(out) :: aa(maxn)                 ! mass of specie (amu)
  integer, intent(out) :: izc(maxn)             ! atomic number of species
                                                ! (-1 for electrons)
!
!---------------------------------
!
  integer n_thi,n_thx,n_bi,n_rfi,n_fusi,inum,iprev
!
  integer i,izz,iaa,iu
!
!---------------------------------
!
  abray=' '
  itype=0
  ifast=0
  zz=0
  aa=0
  izc=0
!
  call rd_nspecies(n_species,n_thi,n_thx,n_bi,n_rfi,n_fusi)
!
  if(n_species.lt.0) return
  if(n_species.gt.maxn) then
     call zermsg(' ?rd_species:  passed array dimension too small!')
     n_species=-1
     return
  endif
!
  zlbla(1)='electrons'
  abray(1,1)='NE'
  abray(2,1)='TE'
  itype(1)=-1
  ifast(1)=0
  zz(1)=-1
  aa(1)=1.0/1836.2
  izc(1) = -1
  iprev=1
!
!  thermal non impurities
!
  inum=n_thi
  call rd_thspec(inum,abray(1:4,iprev+1:iprev+inum), &
       zz(iprev+1:iprev+inum),aa(iprev+1:iprev+inum),&
       izc(iprev+1:iprev+inum),n_thi)
  itype(iprev+1:iprev+inum)=1
  ifast(iprev+1:iprev+inum)=0
  do i=1,inum
     zlbla(iprev+i)=abray(1,iprev+i)(2:)
  enddo
  iprev=iprev + inum
!
!  impurities
!
  inum=n_thx
  if(n_thx.gt.0) then
     call rd_thxspec(inum,abray(1:4,iprev+1:iprev+inum), &
          zz(iprev+1:iprev+inum),aa(iprev+1:iprev+inum), &
          izc(iprev+1:iprev+inum),n_thx)
     if((inum.eq.1).and.(zz(iprev+inum).le.0.001)) then
        itype(iprev+inum)=3
     else
        itype(iprev+1:iprev+inum)=2
     endif
     ifast(iprev+1:iprev+inum)=0
     do i=1,inum
        if(abray(1,iprev+i)(1:5).eq.'NIMP_') then
           zlbla(iprev+i)=abray(1,iprev+i)(6:10)
           iu=len_trim(zlbla(iprev+i))
           zlbla(iprev+i)(iu+1:)='_impurity'
        else
           zlbla(iprev+i)='model_impurity'
        endif
     enddo
     iprev=iprev + inum
  endif
!
!  beam ions
!
  inum=n_bi
  if(n_bi.gt.0) then
     call rd_bmspec(inum,abray(1:4,iprev+1:iprev+inum), &
          zz(iprev+1:iprev+inum),aa(iprev+1:iprev+inum),&
          izc(iprev+1:iprev+inum),n_bi)
     itype(iprev+1:iprev+inum)=4
     ifast(iprev+1:iprev+inum)=1
     do i=1,inum
        iu=index(abray(1,iprev+i),'_')
        zlbla(iprev+i)=abray(1,iprev+i)(iu+1:)
        if(zlbla(iprev+i).eq.'3') zlbla(iprev+i)='HE3'
        if(zlbla(iprev+i).eq.'4') zlbla(iprev+i)='HE4'
        if(zlbla(iprev+i).eq.'P') zlbla(iprev+i)='H'
        iu=len_trim(zlbla(iprev+i))
        zlbla(iprev+i)(iu+1:)='_beam_ion'
     enddo
     iprev=iprev + inum
  endif
!
!  rf tail ions
!
  inum=n_rfi
  if(n_rfi.gt.0) then
     call rd_rfspec(inum,abray(1:4,iprev+1:iprev+inum), &
          zz(iprev+1:iprev+inum),aa(iprev+1:iprev+inum),&
          izc(iprev+1:iprev+inum),n_rfi)
     itype(iprev+1:iprev+inum)=5
     ifast(iprev+1:iprev+inum)=2
     do i=1,inum
        iu=index(abray(1,iprev+i),'_')
        zlbla(iprev+i)=abray(1,iprev+i)(iu+1:)
        if(zlbla(iprev+i).eq.'3') zlbla(iprev+i)='HE3'
        if(zlbla(iprev+i).eq.'4') zlbla(iprev+i)='HE4'
        if(zlbla(iprev+i).eq.'P') zlbla(iprev+i)='H'
        iu=len_trim(zlbla(iprev+i))
        zlbla(iprev+i)(iu+1:)='_rf_minority'
     enddo
     iprev=iprev + inum
  endif
!
!  fusion products
!
  inum=n_fusi
  if(n_fusi.gt.0) then
     call rd_fuspec(inum,abray(1:4,iprev+1:iprev+inum), &
          zz(iprev+1:iprev+inum),aa(iprev+1:iprev+inum),&
          izc(iprev+1:iprev+inum),n_fusi)
     itype(iprev+1:iprev+inum)=6
     ifast(iprev+1:iprev+inum)=1
     do i=1,inum
        iu=index(abray(1,iprev+i),'_')
        zlbla(iprev+i)=abray(1,iprev+i)(iu+1:)
        if(zlbla(iprev+i).eq.'3') zlbla(iprev+i)='HE3'
        if(zlbla(iprev+i).eq.'4') zlbla(iprev+i)='HE4'
        if(zlbla(iprev+i).eq.'P') zlbla(iprev+i)='H'
        iu=len_trim(zlbla(iprev+i))
        zlbla(iprev+i)(iu+1:)='_fusion_product'
     enddo
     iprev=iprev + inum
  endif
!
end subroutine rd_species
!------------------------------------------------------------
subroutine rd_thspec(maxn,abray,zz,aa,izc,ngot)
!
!  return list of profile names & Z & A of thermal species
!  non impurity
!
  implicit NONE
!
! input:
  integer maxn                 ! max no. of species
!
! output:
  character*10 abray(4,maxn)   ! TRANSP profile names
  real zz(maxn)                ! Z of species
  real aa(maxn)                ! A of species
  integer izc(maxn)            ! atomic number of species
  integer ngot                 ! number of species found (-1 if error).
!
!------------------------------
  character*10 zabr,ztabr
  character*32 zunits
  character*64 zlabel
  integer imulti,istype
!
  integer, parameter :: max_species=126  ! matches CPLOTR NAXMGF parameter now
!
  character*10 mgfuns(max(maxn,max_species))
  integer isigns(max(maxn,max_species))
  integer ias(max(maxn,max_species))
  integer izs(max(maxn,max_species))
!
  integer infuns
  integer ierr
!
  integer i,j
!
!------------------------------
!
!  This code is specific to TRANSP and may need to be updated
!  if TRANSP changes, e.g. to have different ion temperatures for
!  each thermal specie.
!
  abray=' '
  zz=0
  aa=0
  izc=0
  ngot=-1
!
!  find TI to use for non-impurity thermal species
!
  ztabr='TMJ'
  call rplabel(ztabr,zlabel,zunits,imulti,istype)
  if((imulti.ne.0).or.(istype.le.0)) then
     ztabr='TI'
     call rplabel(ztabr,zlabel,zunits,imulti,istype)
     if((imulti.ne.0).or.(istype.le.0)) then
        call zermsg('?rd_thspec:  neither "TI" nor "TMJ" exist.')
        return         ! need TI to be defined...
     endif
  endif
 
  zabr='PDENS'
  call rplabel(zabr,zlabel,zunits,imulti,istype)
  if((imulti.le.0).or.(istype.le.0)) return
!
!  OK: get contents of PDENS multigraph
!
  call rpmulti(zabr,istype,zlabel,zunits,infuns,isigns,mgfuns,ierr)
  if(ierr.ne.0) return
!
!  derive the thermal species list from the PDENS contents
!
  call rdi_ckpdens(infuns,mgfuns,ias,izs,ierr)
  if(ierr.ne.0) then
     call zermsg(' ?rdi_ckpdens:  unknown Z and A for: '//mgfuns(ierr))
     return
  endif
!
  j=0
  do i=1,infuns
     if(izs(i).gt.0) then
        j=j+1
        if(j.gt.maxn) then
           call zermsg( &
                ' ?rd_thspec: array dimension too small for species list.')
           return
        endif
        abray(1,j)=mgfuns(i)
        abray(2,j)=ztabr
        zz(j)=izs(i)
        aa(j)=ias(i)
        izc(j)=izs(i)
     endif
  enddo
  ngot=j
!
  return
!------------------------------------------------------------
end subroutine rd_thspec
!------------------------------------------------------------
subroutine rd_thxspec(maxn,abray,zz,aa,izc,ngot)
!
!  return list of profile names & Z & A of thermal species
!  impurity
!    return zz=aa=0 if this is a TRANSP-style single model
!    impurity with zz and aa functions of time
!
  use periodic_table_mod
!
  implicit NONE
!
! input:
  integer maxn                 ! max no. of species
!
! output:
  character*10 abray(4,maxn)   ! TRANSP profile names
  real zz(maxn)                ! Z of species
  real aa(maxn)                ! A of species
  integer izc(maxn)            ! atomic number of species
  integer ngot                 ! number of species found
!
!------------------------------------------------------------
!
  character*32 zunits
  character*64 zlabel
  integer imulti,istype
!
!
  character*10 ztabr,zabr
  character*10, dimension(:), allocatable :: zimplist
  integer isize,jsize,ierr,izimp,iznum,i,j,idum
  real aimp
!------------------------
!
!  This code is specific to TRANSP and may need to be updated
!  if TRANSP changes, e.g. to have different ion temperatures for
!  each thermal specie.
!
  abray=' '
  zz=0
  aa=0
  izc=0
  ngot=-1
!
!  find TI to use for non-impurity thermal species
!
  ztabr='TX'
  call rplabel(ztabr,zlabel,zunits,imulti,istype)
  if((imulti.ne.0).or.(istype.le.0)) then
     ztabr='TI'
     call rplabel(ztabr,zlabel,zunits,imulti,istype)
     if((imulti.ne.0).or.(istype.le.0)) then
        call zermsg('?rd_thspec:  neither "TI" nor "TX" exist.')
        return         ! need TI to be defined...
     endif
  endif
!
  call rpnlist('profile:NIMP_',1,isize)
!
! check for multiple impurities
!
  if(isize.gt.0) then
     allocate(zimplist(isize))
     call rplist('profile:NIMP_',1,zimplist,isize,jsize,ierr)
     if(ierr.ne.0) then
        call zermsg( &
            'TRANSP impurity profile information access error.')
        deallocate(zimplist)
        return
     endif
  endif
!
  jsize=0
  do i=1,isize
     call uupper(zimplist(i))
     if((zimplist(i)(1:5).eq.'NIMP_').and. &
          (index(zimplist(i)(6:10),'_').gt.0)) then
!  label has the TRANSP form for a specific impurity ion
        call rplabel(zimplist(i),zlabel,zunits,imulti,istype)
        j=index(zlabel,' ')-1
!  impuirty code has form like "C+6", "Ni+28", etc.
        call inv_periodic_table(zlabel(1:j),.TRUE., &
             iznum,aimp,izimp)
        jsize=jsize+1
        if(jsize.gt.maxn) then
           call zermsg( &
                ' ?rd_thspec: array dimension too small for species list.')
           deallocate(zimplist)
           return
        endif
        abray(1,jsize)=zimplist(i)
        abray(2,jsize)=ztabr
        aa(jsize)=aimp
        zz(jsize)=izimp
        izc(jsize)=iznum
     endif
  enddo
!
! no impurities or one "model impurity (tokamakium)"
!
  if(jsize.eq.0) then
     ngot=0
     zabr='NIMP'
     call rplabel(zabr,zlabel,zunits,imulti,istype)
     if((imulti.ne.0).or.(istype.le.0)) then
        continue
     else
        jsize=1
        abray(1,1)=zabr
        abray(2,1)=ztabr
                     ! single model impurity; zz & aa can vary in time
     endif
  endif
!
  ngot=jsize
  deallocate(zimplist,stat=idum)
  return
end subroutine rd_thxspec
!------------------------------------------------------------
subroutine rd_bmspec(maxn,abray,zz,aa,izc,ngot)
!
!  return list of profile names & Z & A of beam species
!
  implicit NONE
!
! input:
  integer maxn                 ! max no. of species
!
! output:
  character*10 abray(4,maxn)   ! TRANSP profile names
  real zz(maxn)                ! Z of species
  real aa(maxn)                ! A of species
  integer izc(maxn)            ! atomic number of species
  integer ngot                 ! number of species found
!
!---------------------
  call rd_fspec(maxn,abray,zz,aa,izc,ngot, &
       'BDENS_','EBEAM_','UBPRP_','UBPAR_')
  return
!------------------------------------------------------------
end subroutine rd_bmspec
!------------------------------------------------------------
subroutine rd_rfspec(maxn,abray,zz,aa,izc,ngot)
!
!  return list of profile names & Z & A of rf tail ion species
!
  implicit NONE
!
! input:
  integer maxn                 ! max no. of species
!
! output:
  character*10 abray(4,maxn)   ! TRANSP profile names
  real zz(maxn)                ! Z of species
  real aa(maxn)                ! A of species
  integer izc(maxn)            ! atomic number of species
  integer ngot                 ! number of species found
!
!---------------------
  call rd_fspec(maxn,abray,zz,aa,izc,ngot, &
       'NMINI_','TMINI_','EMINPER_','EMINPAR_')
  return
!------------------------------------------------------------
end subroutine rd_rfspec
!------------------------------------------------------------
subroutine rd_fuspec(maxn,abray,zz,aa,izc,ngot)
!
!  return list of profile names & Z & A of fusion product ion species
!
  implicit NONE
!
! input:
  integer maxn                 ! max no. of species
!
! output:
  character*10 abray(4,maxn)   ! TRANSP profile names
  real zz(maxn)                ! Z of species
  real aa(maxn)                ! A of species
  integer izc(maxn)            ! atomic number of species
  integer ngot                 ! number of species found
!
!---------------------
  call rd_fspec(maxn,abray,zz,aa,izc,ngot, &
       'FDENS_','EFUSN_','UFPRP_','UFPAR_')
  return
!------------------------------------------------------------
end subroutine rd_fuspec
!============================================================
subroutine rd_fspec(maxn,abray,zz,aa,izc,ngot,nroot,eroot,pproot,plroot)
!
!  return list of profile names & Z & A of fast ion species
!
  implicit NONE
!
! input:
  integer maxn                 ! max no. of species
!
! output:
  character*10 abray(4,maxn)   ! TRANSP profile names
  real zz(maxn)                ! Z of species
  real aa(maxn)                ! A of species
  integer izc(maxn)            ! atomic number of species
  integer ngot                 ! number of species found
!
! input:
  character*(*) nroot          ! root of density profile names
  character*(*) eroot          ! root of <E> or T profile names
  character*(*) pproot         ! root of perp energy profile names
  character*(*) plroot         ! root of pll energy profile names
!
!---------------------
  integer iln,ile,ilpp,ilpl
  integer isize,isizex,isize2,idum,ierr,i,j
!
  character*10, dimension(:), allocatable :: zflist
  character*5, dimension(:), allocatable :: zcod
  character*5 ztest
!---------------------
!
  abray=' '
  zz=0
  aa=0
  izc=0
  ngot=-1
!
  iln=len_trim(nroot)
  ile=len_trim(eroot)
  ilpp=len_trim(pproot)
  ilpl=len_trim(plroot)
!
  call rpnlist('profile:'//nroot(1:iln),1,isize)
  if(isize.eq.0) then
     ngot=0
     return
  endif
!
! found some...
!
  isizex=max(50,isize)
  allocate(zflist(isizex),zcod(isize))
!
  call rplist('profile:'//nroot(1:iln),1,zflist,isize,idum,ierr)
!
  j=0
  do i=1,isize
     zcod(i)=zflist(i)(iln+1:)
     if((zcod(i).eq.'H').or.(zcod(i).eq.'D').or.(zcod(i).eq.'T').or. &
          (zcod(i).eq.'P').or.(zcod(i).eq.'3').or.(zcod(i).eq.'4').or. &
          (zcod(i).eq.'N').or.(zcod(i).eq.'A').or.(zcod(i).eq.'K').or. &
          (zcod(i).eq.'X')) then
        j=j+1
        abray(1,j)=zflist(i)
        zcod(j)=zcod(i)
        if((zcod(i).eq.'H').or.(zcod(i).eq.'P')) then
           zz(j)=1
           aa(j)=1
           izc(j)=1
        else if(zcod(i).eq.'D') then
           zz(j)=1
           aa(j)=2
           izc(j)=1
        else if(zcod(i).eq.'T') then
           zz(j)=1
           aa(j)=3
           izc(j)=1
        else if(zcod(i).eq.'3') then
           zz(j)=2
           aa(j)=3
           izc(j)=2
        else if(zcod(i).eq.'4') then
           zz(j)=2
           aa(j)=4
           izc(j)=2
        else if(zcod(i).eq.'N') then
           zz(j)=0  ! avg <Z> to be fetched
           aa(j)=20.02  ! Neon, 20.02 x proton mass
           izc(j)=10    ! Neon
        else if(zcod(i).eq.'A') then
           zz(j)=0  ! avg <Z> to be fetched
           aa(j)=39.63  ! Argon, 39.63 x proton mass
           izc(j)=18    ! Argon
        else if(zcod(i).eq.'K') then
           zz(j)=0  ! avg <Z> to be fetched
           aa(j)=83.14  ! Krypton, 83.14 x proton mass
           izc(j)=36    ! Krypton
        else if(zcod(i).eq.'X') then
           zz(j)=0  ! avg <Z> to be fetched
           aa(j)=130.27 ! Xenon, 130.27 x proton mass
           izc(j)=54    ! Xenon
        endif
     else
        call zermsg( &
             ' %rd_species: apparent fast ion density not recognized:'// &
             zflist(i))
     endif
  enddo
  ngot=j
!
  if(ngot.gt.0) then
!
!  get <E> or "T" profile names...
     call rplist('profile:'//eroot(1:ile),1,zflist,isizex,isize2,ierr)
     if(ierr.ne.0) then
        call zermsg('?rd_fspec:  unexpected code error!')
        deallocate(zflist,zcod)
        return
     endif
     do i=1,isize2
        ztest=zflist(i)(ile+1:)
        do j=1,ngot
           if(ztest.eq.zcod(j)) then
              abray(2,j)=zflist(i)
           endif
        enddo
     enddo
!
!  get <Eperp> or Uperp profile names...
     call rplist('profile:'//pproot(1:ilpp),1,zflist,isizex,isize2,ierr)
     if(ierr.ne.0) then
        call zermsg('?rd_fspec:  unexpected code error!')
        deallocate(zflist,zcod)
        return
     endif
     do i=1,isize2
        ztest=zflist(i)(ilpp+1:)
        do j=1,ngot
           if(ztest.eq.zcod(j)) then
              abray(3,j)=zflist(i)
           endif
        enddo
     enddo
!
!  get <Epll> or Upll profile names...
     call rplist('profile:'//plroot(1:ilpl),1,zflist,isizex,isize2,ierr)
     if(ierr.ne.0) then
        call zermsg('?rd_fspec:  unexpected code error!')
        deallocate(zflist,zcod)
        return
     endif
     do i=1,isize2
        ztest=zflist(i)(ilpl+1:)
        do j=1,ngot
           if(ztest.eq.zcod(j)) then
              abray(4,j)=zflist(i)
           endif
        enddo
     enddo
  endif

  deallocate(zflist,zcod)

  return
end subroutine rd_fspec

!------------------------------------------------------------
subroutine rd_th_scedg(nmax,agas,arcy,zz,aa,izc,ngot)

  !  return list of edge sources & Z & A of thermal species
  !  non impurity

  !  NOTE: these are "effective sources"; the numerical values correspond
  !  to N/sec ionized inside the plasma core; this number is somewhat smaller
  !  than the number of neutrals entering across the plasma boundary, the
  !  difference being due to escaping charge exchange neutrals

  implicit NONE

  ! input:
  integer nmax                 ! max no. of species

  ! output:
  character*10 agas(nmax)      ! TRANSP gasflow sources (rplot names)
  character*10 arcy(nmax)      ! TRANSP recycling sources (rplot names)
  real zz(nmax)                ! Z of species
  real aa(nmax)                ! A of species
  integer izc(nmax)            ! atomic number of species
  integer ngot                 ! number of species found (-1 if error).

  !------------------------------

  character*10 :: abray(4,nmax)
  character*5 :: slbl
  character*1 :: suffix1
  integer :: i,igot

  !------------------------------

  call rd_thspec(nmax,abray,zz,aa,izc,ngot)
  if(ngot.lt.0) return

  igot=ngot
  do i=1,igot
     slbl=abray(1,i)(2:6)
     call uupper(slbl)
     if(slbl.eq.'H') then
        agas(i)='GASH'
        arcy(i)='RCYH'
     else if(slbl.eq.'D') then
        agas(i)='GASD'
        arcy(i)='RCYD'
     else if(slbl.eq.'T') then
        agas(i)='GAST'
        arcy(i)='RCYT'
     else if(slbl.eq.'HE3') then
        agas(i)='GAS3'
        arcy(i)='RCY3'
     else if(slbl.eq.'HE4') then
        agas(i)='GAS4'
        arcy(i)='RCY4'
     else if(slbl.eq.'LITH') then
        agas(i)='GASL'
        arcy(i)='RCYL'
     else
        call zermsg( &
             ' ?rd_th_scedg: unrecognized thermal specie density ID: '// &
             abray(1,i))
        ngot=-1
     endif
  enddo

end subroutine rd_th_scedg
