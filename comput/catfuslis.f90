!-----------------------------------------------------------------------
!  CATFUSLIS -- PROVIDE DATA ON FUSION BY FUSION PRODUCT
!
!  modified (dmc Feb 2002 -- catfuslis.for -- now based on arguments
!  instead of COMMON).  Coordinate any changes with comput/catfus*.*
!
!  FUSION REACTION INDEXING IS AS DESCRIBED IN COMMENTS OF subroutine
!  FUSION (***** see nubeam/fusion.for comments *****).
!  see also comput/catfus*.*
!
!  Caution:  Feb 2002 (dmc):  although catfus has a notion of a "D"
!  RF minority, it is not supported here.
!
!  THIS ROUTINE RECEIVES AS ARGUMENT THE INDEX TO A FAST ION SPECIES
!  USED TO MODEL FUSION PRODUCT FAST IONS.  WHICH SPECIES INDICATES
!  THE FUSION PRODUCT AND THUS ONE OR MORE FUSION REACTIONS INVOLVED.
!
!  THIS ROUTINE GIVES AN ESTIMATE OF THE INITIAL ENERGY OF THE FUSION
!  PRODUCT ** DMC SEPT 1990 ** THIS SHOULD BE CONSISTENT WITH THE MAX
!  ENERGY SET IN S.R. DATCHK FOR THE FUSION PRODUCT SPECIES' COMPUTED
!  DISTRIBUTION FUNCTION.
!
!  ENERGY DATA IS FROM THE NRL PLASMA FORMULARY ... DMC SEPT 1990
!
subroutine CATFUSLIS_r8(IDENT,ILREACT,ZEREACT,idim, &
     nlftrit,nlfhe3,nlrfdeut,nlrftrit,nlrfhe3,ierr)
  !
  use iso_c_binding, only: fp => c_double
  use catfus_mod
  implicit none
  !
  !  input:
  character(len=*) :: ident            ! fusion product ID
  !             "P" "T" "He3" "He4" "Alpha"
  !
  integer :: idim                      ! dimension of output arrays
  !
  logical :: nlftrit                   ! TRUE if fusion product T present
  logical :: nlfhe3                    ! TRUE if fusion product He3 present
  !
  logical :: nlrfdeut                  ! TRUE if rf minority D present
  logical :: nlrftrit                  ! TRUE if rf minority T present
  logical :: nlrfhe3                   ! TRUE if rf minority He3 present
  !
  !  OUTPUT:
  !
  !  SET .TRUE. FOR REACTIONS THAT PRODUCE A DESIRED FUSION PRODUCT:
  logical, dimension(idim) :: ILREACT
  !
  !  INITIAL ENERGY (EV) OF FUSION PRODUCT SO PRODUCED:
  real(fp), dimension(idim) :: ZEREACT
  !
  integer ierr                      ! completion code: 0=OK
  !
  !-------------------------------------------------------------
  character(len=6) :: test
  !-------------------------------------------------------------
  !
  !  INPUT:
  !
  !    IDENT = FUSION PRODUCT INDEX (MUST BE NON-ZERO)
  !
  test=ident
  call uupper(test)                 ! uppercase conversion
  !
  ierr=0
  ilreact=.FALSE.
  zereact=0.0_fp
  !
  if(idim.lt.nreact) then
    write(6,*) ' ?CatFusLis:  idim = ',idim,' too small.'
    write(6,*) '  idim .ge. nreact = ',nreact,' was expected.'
    ierr=1
    return
  end if
  !
  if(test.eq.'T') then
    !
    !  TRITIUM COMES FROM THE D+D REACTION (PROTON BRANCH)
    !
    ILREACT(3)=.TRUE.
    ZEREACT(3)=1.01e6_fp
    !
  else if (test.eq.'P') then
    !
    !  Proton COMES FROM THE D+D REACTION (PROTON BRANCH)
    !
    ilreact(3) = .true.
    zereact(3) = 3.02e6_fp
    !
  else if (test.eq.'P14') then
    !  D+HE3 AND HE3+D (He3 thermal)

    ilreact(2) = .true.
    zereact(2) = 14.7e6_fp
    ilreact(8) = .true.
    zereact(8) = 14.7e6_fp

    if(NLFHE3) then
      !  fusion product He3 -- He3+D burnup reactions
      ILREACT(10)=.TRUE.
      ZEREACT(10)=14.7e6_fp          ! proton energy
    end if
    !
    if(NLRFHE3) then
      !  RF minority He3
      !  He3 + D target
      ILREACT(12)=.TRUE.
      ZEREACT(12)=14.7e6_fp          ! proton energy
    end if
    !
  else if(test.eq.'HE3') then
    !
    !  HELIUM-3 COMES FROM THE D+D REACTION (NEUTRON BRANCH)
    !
    ILREACT(4)=.TRUE.
    ZEREACT(4)=0.82e6_fp
    !
  else if((test.eq.'HE4').or.(test.eq.'ALPHA')) then
    !
    !  HELIUM-4 COMES FROM D+T, D+HE3, AND T+T REACTIONS
    !
    !  D+T AND T+D --
    ILREACT(1)=.TRUE.
    ZEREACT(1)=3.5e6_fp
    ILREACT(7)=.TRUE.
    ZEREACT(7)=3.5e6_fp
    !
    !  D+HE3 AND HE3+D
    ILREACT(2)=.TRUE.
    ZEREACT(2)=3.6e6_fp
    ILREACT(8)=.TRUE.
    ZEREACT(8)=3.6e6_fp
    !
    !  T+T (ACTUALLY I doN'T KNOW WHAT ENERGY TO GIVE ...)
    !
    ILREACT(5)=.TRUE.
    ZEREACT(5)=3.75e6_fp
    !
    if(NLFTRIT) then
      !  fusion product T -- T+D burnup reactions
      ILREACT(9)=.TRUE.
      ZEREACT(9)=3.5e6_fp
    end if
    !
    if(NLFHE3) then
      !  fusion product He3 -- He3+D burnup reactions
      ILREACT(10)=.TRUE.
      ZEREACT(10)=3.6e6_fp
    end if
    !
    if(NLRFTRIT) then
      !  RF minority T -- T+D burnup reactions
      ILREACT(11)=.TRUE.
      ZEREACT(11)=3.5e6_fp
    end if
    !
    if(NLRFHE3) then
      !  RF minority He3 -- He3+D burnup reactions
      ILREACT(12)=.TRUE.
      ZEREACT(12)=3.6e6_fp
    end if
    !
    if(NLRFDEUT) then
      !  RF minority D -- D+T burnup reactions
      ILREACT(13)=.TRUE.
      ZEREACT(13)=3.5e6_fp
      !  RF minority D -- D+He3 burnup reactions
      ILREACT(14)=.TRUE.
      ZEREACT(14)=3.6e6_fp
    end if
  else
    ierr=2
    write(6,*) ' ?CatFusLis:  ident="',ident,'" not recognized.'
  end if
  !
  return
end subroutine CATFUSLIS_r8

!--------------------------------------------
!  REAL interface
!
subroutine catfuslis(IDENT,ILREACT,ZEREACT,idim, &
     nlftrit,nlfhe3,nlrfdeut,nlrftrit,nlrfhe3,ierr)
  !
  use iso_c_binding, only: fp => c_double, sp => c_float
  implicit none
  !
  !  input:
  character(len=*) :: ident            ! fusion product ID
  !             "P" "T" "He3" "He4" "Alpha"
  !
  integer :: idim                      ! dimension of output arrays
  !
  logical :: nlftrit                   ! TRUE if fusion product T present
  logical :: nlfhe3                    ! TRUE if fusion product He3 present
  !
  logical :: nlrfdeut                  ! TRUE if rf minority D present
  logical :: nlrftrit                  ! TRUE if rf minority T present
  logical :: nlrfhe3                   ! TRUE if rf minority He3 present
  !
  !  OUTPUT:
  !
  !  SET .TRUE. FOR REACTIONS THAT PRODUCE A DESIRED FUSION PRODUCT:
  logical, dimension(idim) :: ILREACT
  !
  !  INITIAL ENERGY (EV) OF FUSION PRODUCT SO PRODUCED:
  real(sp), dimension(idim) :: ZEREACT
  !
  integer :: ierr                      ! completion code: 0=OK
  !
  !-------------------------------------------------
  real(fp), dimension(idim) :: zereact_r8
  !-------------------------------------------------
  !
  zereact_r8=0.0_fp
  !
  call catfuslis_r8(IDENT,ILREACT,ZEREACT_R8,idim, &
       nlftrit,nlfhe3,nlrfdeut,nlrftrit,nlrfhe3,ierr)

  zereact = zereact_r8
  !
  return
end subroutine catfuslis
