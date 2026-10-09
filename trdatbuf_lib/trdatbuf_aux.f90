module trdatbuf_aux

  !  define auxilliary types to facilitate access to trdatbuf data...

  integer, parameter :: tdbaux_maxnb = 200

  type pwrget

     !  access to power vs. time datasets (e.g. beam or RF powers)

     character*3 :: item    ! NBI data item to fetch/average
     real*8 :: ztime1       ! start of time averaging window
     real*8 :: ztime2       ! end of time averaging window
     logical :: pweight     ! .TRUE. for power weighted average
     integer :: nbeam       ! no. of beams -- expected to match trdatbuf nbdata
     real*8 :: tbon(tdbaux_maxnb)  ! on times for beams (1:nbeam)
     real*8 :: tboff(tdbaux_maxnb) ! off times for beams (1:nbeam)
     integer :: nguard      ! if not equal to 123456789 assume not initialized.

     real*8 :: zparam(tdbaux_maxnb) ! return vector...

  end type pwrget

  type pelget
     real*8 :: tpel         ! pellet injection time (seconds)
     real*8 :: tpele        ! pellet injector end time (seconds)
     real*8 :: freqpel      ! pellet injection frequency (Hz)
     real*8 :: apel         ! A of pellet main component (H isotope)
     real*8 :: apel2        ! A of 2nd pellet component (H isotope)
     real*8 :: fpel2        ! fraction of 2nd component (usually 0.0)
     integer :: nmpela      ! pellet treatment option request (to TRANSP)
                            ! 0: density discontinuity, no ablation calculation
                            ! 1: density discontinuity, ablation calculated
                            ! 2: no density continuity, ablation calculated
     real*8 :: plrsta       ! major radius of injection into VV (cm)
     real*8 :: plysta       ! elevation above midplane of injection (cm)
     real*8 :: plphia       ! toroidal angle of injection (degrees)
                            ! 0 = dead perp, inward; + -> with plasma current
     real*8 :: plthea       ! poloidal angle of injection (degrees)
                            ! 0 = horizonta, + -> pointed upward
     real*8 :: pelvel       ! pellet velocity (cm/sec)
     real*8 :: pelrad       ! pellet radius (cm)
     integer :: kpellet     ! pellet ablation model
                            
  end type pelget

  type profget

     !  access to profile datasets

     !----------------------------------------
     !  input information...

     character*3 :: item      ! profile item desired...  e.g. "NER"

     integer :: iselect       ! profile selection w/in item (3d sets only)
                              !  (nb dmc June 2005 -- only an issue for
                              !  multi-impurity data) -- default value =0 OK

     integer :: nzones        ! number of zones in target grid

     !  to map data vs. normalized toroidal flux or norm. sqrt. tor. flux:

     real*8, dimension(:), pointer :: xibdys  ! (nzones+1) x @ zone bdys
     !      x is sqrt(Phi/Philim) ref enclosed toroidal flux at zone surfaces
     !      xibdys(1)=0; xibdys(nzones+1)=1 expected.

     !  to map data vs. normalized poloidal flux or norm. sqrt. pol. flux:

     real*8, dimension(:), pointer :: plflxg  ! (nzones+1) enclosed poloidal
     !      flux "psi" (Wb/2pi); Bpol = (1/R)*grad(psi)
     !      at zone bdys; plflxg(1)=0 expected.

     !  to map data vs. major or minor radius coordinate:

     real*8, dimension(:), pointer :: rmajmp  ! (2*nzones+1) major radii at
     !      flux zone boundary / midplane intercepts (cm).  RMAJMP(nzones+1)
     !      gives the magnetic axis; RMAJMP(1)= plasma boundary at midplane,
     !      small major radius side; RMAJMP(2*nzones+1) = plasma boundary at
     !      midplane, large major radius side.

     !  to map data vs. ECE frequency:

     real*8, dimension(:), pointer :: bmidp   ! (2*nzones+1) total B field (T)
     !      vs. major radius (rmajmp) along the midplane

     !  for all data items:

     real*8 :: time           ! time at which data is required (seconds)

     real*8 :: delta_t        ! average from time-delta_t to time+delta_t
                              ! default: delta_t = ZERO

     logical :: ibdy          ! T to interpolate to boundaries & fill in zones
                              ! F to interpolate to zone centers & fill in bdys
     !----------------------------------------
     !  internal:

     logical :: idebug        ! T to request extra output pertaining to
                              ! mapping or symmetrization of data
                              ! default: F

     logical :: iece          ! T for data vs. EC resonance frequency
     ! (trdatbuf object knows the harmonic frequency)

     integer :: nguard        ! used to verify object initialization
     integer :: nrmaj         ! 2*nzones + 1 -- rmajmp dimension
     integer :: nrsym         ! 4*nzones + 5 -- rmjsym dimension

     !  the following can be treated as output if desired by caller...
     !  these are only defined if idebug=T

     real*8, dimension(:), pointer :: xilmp
                              ! signed normalized sqrt tor flux on rmajmp grid

     !----------------------------------------
     !  output profiles: data & asymmetry profiles in data units...

     real*8, dimension(:), pointer :: data_zc ! (nzones+1) profile at zone
     !  centers, including a dummy zone, 1/2 zone width beyond plasma boundary

     real*8, dimension(:), pointer :: data_zb ! (nzones+1) profile at zone
     !  boundaries

     !----------------------------------------
     !  additional outputs created only if idebug=T...
     !
     !  for data that is input 2-sided vs. major radius and symmetrized
     !  by in/out averaging (nsyxxx=2 option):

     real*8, dimension(:), pointer :: asym_zc ! (nzones+1) profile asymmetry
     !  at zone centers (averaged value - measured value on large major radius
     !  side)

     real*8, dimension(:), pointer :: asym_zb ! (nzones+1) profile asymmetry
     !  at zone boundaries

     !  for data that is input 2-sided vs. major radius and symmetrized by
     !  the "slice and stack" method:

     real*8, dimension(:), pointer :: shift ! (nzones+1) inferred shift of
     !  flux surfaces relative to plasma boundary (cm)

     !  data mapped to fine major radius grid (in data units)

     real*8, dimension(:), pointer :: rmjsym
                              ! fine (nrsym) radial grid for symmetrization
                              ! comparison outputs
     real*8, dimension(:), pointer :: xirsym
                              ! signed normalized sqrt tor flux on rmjsym grid
                              ! range slightly greater than [-1,-1].

     real*8, dimension(:), pointer :: datrsym ! symmetrized data vs. rmjsym
     
     real*8, dimension(:), pointer :: datusym ! unsymmetrized data va. rmjsym

     real*8 :: ecegap   ! maximal B(R) monotonicity gap (cm)
     !  in high beta tokamak shots dB/dR can change sign due to diamagnetic
     !  effects, screening part of the plasma from view of an ECE temperature
     !  diagnostic; this number estimates the maximum such gap occurring
     !  it is calculated only if t%iece = .TRUE. and t%idebug = .TRUE.

  end type profget

  !----------------------------------------------

  contains

    subroutine tdb_pwrget_init(znbi)

      type(pwrget) :: znbi

      znbi%item = '???'
      znbi%ztime1 = 0
      znbi%ztime2 = -1
      znbi%pweight = .FALSE.
      znbi%nbeam = 0
      znbi%nguard = 123456789
      znbi%tbon(1) = 0
      znbi%tboff(1) = -1

    end subroutine tdb_pwrget_init

    subroutine tdb_profget_init(t,izones,idebug)

      type(profget) :: t
      integer, intent(in) :: izones
      logical, intent(in) :: idebug

      t%nguard = 123456789         ! indicate initialization of item
      t%item = '???'
      t%iselect = 0

      t%nzones = izones
      call tdb_nzones(t%nzones,t%nrmaj,t%nrsym)

      t%idebug = idebug

      allocate(t%xibdys(t%nzones+1)); t%xibdys = 0
      allocate(t%plflxg(t%nzones+1)); t%plflxg = 0
      allocate(t%rmajmp(t%nrmaj)); t%rmajmp = 0
      allocate(t%bmidp(t%nrmaj)); t%bmidp = 0

      if(t%idebug) then
         allocate(t%xilmp(t%nrmaj)); t%xilmp = 0
         allocate(t%xirsym(t%nrsym)); t%xirsym = 0
         allocate(t%rmjsym(t%nrsym)); t%rmjsym = 0
      else
         nullify(t%xilmp)
         nullify(t%xirsym)
         nullify(t%rmjsym)
      endif

      t%time = 0
      t%delta_t = 0
      t%ibdy = .FALSE.
      t%iece = .FALSE.

      !  space for output...

      t%ecegap = 0
      allocate(t%data_zc(t%nzones+1)); t%data_zc = 0
      allocate(t%data_zb(t%nzones+1)); t%data_zb = 0

      if(t%idebug) then
         allocate(t%asym_zc(t%nzones+1)); t%asym_zc = 0
         allocate(t%asym_zb(t%nzones+1)); t%asym_zb = 0
         allocate(t%shift(t%nzones+1)); t%shift = 0
         allocate(t%datrsym(t%nrsym)); t%datrsym = 0
         allocate(t%datusym(t%nrsym)); t%datusym = 0
      else
         nullify(t%asym_zc)
         nullify(t%asym_zb)
         nullify(t%shift)
         nullify(t%datrsym)
         nullify(t%datusym)
      endif
    end subroutine tdb_profget_init

    subroutine tdb_profget_free(t)

      type(profget) :: t

      deallocate(t%xibdys,t%plflxg,t%rmajmp,t%bmidp)
      deallocate(t%data_zc,t%data_zb)
      if(t%idebug) then
         deallocate(t%xilmp,t%xirsym,t%rmjsym)
         deallocate(t%asym_zc,t%asym_zb,t%shift,t%datrsym,t%datusym)
      endif

      t%nguard = 0   ! object de-initialized...

    end subroutine tdb_profget_free

end module trdatbuf_aux
