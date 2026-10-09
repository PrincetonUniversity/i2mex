module periodic_table_mod
  use iso_c_binding, only: fp => c_double
  implicit none

  ! ****** MODIFICATION ******
  ! DMC July 2008 -- "Atomic Mass" outputs in "AMU" are now returned in units of
  ! Proton Mass -- the numerical value is slightly smaller, (1/1.0079) of the AMU
  ! values used by chemists, which are the numbers stored in the periodic table
  ! below.  AMU(Hydrogen)=1.00000 and all other AMUs are scaled accordingly!
  !
  ! RGA June 2013 -- Table represents isotopic average AMU in (u).  Separate subroutine
  ! atom_amu_to_proton() to convert to proton mass while accounting for charge state.
  ! **************************
  !
  ! --- a structure which describes an element ---
  !
  logical :: ignore_user_AMU = .true.      ! change how to handle a user's AMU request.  Should be consistent 
                                           ! with plasma_state_mod::ignore_user_AMU
                                           !   .true. -> old behaviour, ignore user's AMU when the element does not
                                           !             have known isotopes (Z>3)
                                           !   .false.-> allow user AMU for Z>3

  logical :: ignore_electron_AMU = .true.  ! .true. to adjust neutral atom AMU by electron weight

  logical :: ps_params = .true.            ! .true. to use plasma_state compatible parameters, .false. for more accuracy

  type, private :: ELM
    character(len=2) :: symbol   ! symbol for the element
    integer          :: number   ! atomic number of the element
    real(fp)         :: weight   ! isotopic average atomic weight of the element
    integer          :: isoindex ! index of first isotope in isotope list or 0 if none
  end type ELM

  !
  ! --- a structure which describes an isotope ---
  !
  type, private :: ISO
    character(len=2)  :: symbol  ! symbol for the isotope or blank for element symbol
    character(len=13) :: name    ! name for the isotope or blank for element name
    integer           :: number  ! atomic number of the isotope
    real(fp)          :: weight  ! atomic weight of the isotope
  end type ISO

  !
  ! --- parameters ---
  ! http://physics.nist.gov/cuu/Constants/Table/allascii.txt
  !
  real(fp), parameter, private :: electronCharge = 1.602176565e-19_fp    ! (C), also J/eV
  real(fp), parameter, private :: proton_amu     = 1.007276466812_fp     ! (u)
  real(fp), parameter, private :: electron_amu   = 0.000548579909_fp     ! (u)
  real(fp), parameter, private :: neutron_amu    = 1.00866491600_fp      ! (u)
  real(fp), parameter, private :: speedLight     = 2.99792458e8_fp       ! (m/s)
  real(fp), parameter, private :: protonMass     = 1.672621777e-27_fp    ! (kg)

  integer, parameter,private  :: max_elements = 109 ! number of elements

  !
  ! weights are isotopic average in AMU for neutral atom
  !
  type (elm), dimension(max_elements), parameter, private :: elex = (/ &
       elm('Ac', 89,227.0000_fp,0) , elm('Ag', 47,107.8680_fp,0) , elm('Al', 13, 26.9815_fp,0) ,  &
       elm('Am', 95,243.0000_fp,0) , elm('Ar', 18, 39.9480_fp,0) , elm('As', 33, 74.9216_fp,0) ,  &
       elm('At', 85,210.0000_fp,0) , elm('Au', 79,196.9665_fp,0) , elm('B ',  5, 10.8100_fp,0) ,  &
       elm('Ba', 56,137.3400_fp,0) , elm('Be',  4,  9.0122_fp,0) , elm('Bh',107,262.0000_fp,0) ,  &
       elm('Bi', 83,208.9804_fp,0) , elm('Bk', 97,247.0000_fp,0) , elm('Br', 35, 79.9040_fp,0) ,  &
       elm('C ',  6, 12.0110_fp,0) , elm('Ca', 20, 40.0800_fp,0) , elm('Cd', 48,112.4000_fp,0) ,  &
       elm('Ce', 58,140.1200_fp,0) , elm('Cf', 98,251.0000_fp,0) , elm('Cl', 17, 35.4530_fp,0) ,  &
       elm('Cm', 96,247.0000_fp,0) , elm('Co', 27, 58.9332_fp,0) , elm('Cr', 24, 51.9960_fp,0) ,  &
       elm('Cs', 55,132.9054_fp,0) , elm('Cu', 29, 63.5460_fp,0) , elm('Db',105,262.0000_fp,0) ,  &
       elm('Dy', 66,162.5000_fp,0) , elm('Er', 68,167.2600_fp,0) , elm('Es', 99,252.0000_fp,0) ,  &
       elm('Eu', 63,151.9600_fp,0) , elm('F ',  9, 18.9984_fp,0) , elm('Fe', 26, 55.8470_fp,0) ,  &
       elm('Fm',100,257.0000_fp,0) , elm('Fr', 87,223.0000_fp,0) , elm('Ga', 31, 69.7200_fp,0) ,  &
       elm('Gd', 64,157.2500_fp,0) , elm('Ge', 32, 72.5900_fp,0) , elm('H ',  1,  1.0079_fp,2) ,  &
       elm('He',  2,  4.0026_fp,5) , elm('Hf', 72,178.4900_fp,0) , elm('Hg', 80,200.5900_fp,0) ,  &
       elm('Ho', 67,164.9304_fp,0) , elm('Hs',108,265.0000_fp,0) , elm('I ', 53,126.9045_fp,0) ,  &
       elm('In', 49,114.8200_fp,0) , elm('Ir', 77,192.2200_fp,0) , elm('K ', 19, 39.0980_fp,0) ,  &
       elm('Kr', 36, 83.8000_fp,0) , elm('La', 57,138.9055_fp,0) , elm('Li',  3,  6.9410_fp,7) ,  &
       elm('Lr',103,262.0000_fp,0) , elm('Lu', 71,174.9700_fp,0) , elm('Md',101,258.0000_fp,0) ,  &
       elm('Mg', 12, 24.3050_fp,0) , elm('Mn', 25, 54.9380_fp,0) , elm('Mo', 42, 95.9400_fp,0) ,  &
       elm('Mt',109,265.0000_fp,0) , elm('N ',  7, 14.0067_fp,0) , elm('Na', 11, 22.9898_fp,0) ,  &
       elm('Nb', 41, 92.9064_fp,0) , elm('Nd', 60,144.2400_fp,0) , elm('Ne', 10, 20.1790_fp,0) ,  &
       elm('Ni', 28, 58.7000_fp,0) , elm('No',102,259.0000_fp,0) , elm('Np', 93,237.0482_fp,0) ,  &
       elm('O ',  8, 15.9994_fp,0) , elm('Os', 76,190.2000_fp,0) , elm('P ', 15, 30.9738_fp,0) ,  &
       elm('Pa', 91,231.0359_fp,0) , elm('Pb', 82,207.2000_fp,0) , elm('Pd', 46,106.4000_fp,0) ,  &
       elm('Pm', 61,145.0000_fp,0) , elm('Po', 84,209.0000_fp,0) , elm('Pr', 59,140.9077_fp,0) ,  &
       elm('Pt', 78,195.0900_fp,0) , elm('Pu', 94,244.0000_fp,0) , elm('Ra', 88,226.0254_fp,0) ,  &
       elm('Rb', 37, 85.4678_fp,0) , elm('Re', 75,186.2070_fp,0) , elm('Rf',104,261.0000_fp,0) ,  &
       elm('Rh', 45,102.9055_fp,0) , elm('Rn', 86,222.0000_fp,0) , elm('Ru', 44,101.0700_fp,0) ,  &
       elm('S ', 16, 32.0600_fp,0) , elm('Sb', 51,121.7500_fp,0) , elm('Sc', 21, 44.9559_fp,0) ,  &
       elm('Se', 34, 78.9600_fp,0) , elm('Sg',106,263.0000_fp,0) , elm('Si', 14, 28.0860_fp,0) ,  &
       elm('Sm', 62,150.4000_fp,0) , elm('Sn', 50,118.6900_fp,0) , elm('Sr', 38, 87.6200_fp,0) ,  &
       elm('Ta', 73,180.9479_fp,0) , elm('Tb', 65,158.9254_fp,0) , elm('Tc', 43, 97.0000_fp,0) ,  &
       elm('Te', 52,127.6000_fp,0) , elm('Th', 90,232.0381_fp,0) , elm('Ti', 22, 47.9000_fp,0) ,  &
       elm('Tl', 81,204.3700_fp,0) , elm('Tm', 69,168.9342_fp,0) , elm('U ', 92,238.0290_fp,0) ,  &
       elm('V ', 23, 50.9414_fp,0) , elm('W ', 74,183.5000_fp,0) , elm('Xe', 54,131.3000_fp,0) ,  &
       elm('Y ', 39, 88.9059_fp,0) , elm('Yb', 70,173.0400_fp,0) , elm('Zn', 30, 65.3800_fp,0) ,  &
       elm('Zr', 40, 91.2200_fp,0) &
       /)

  character(len=26), parameter, private :: lower  = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter, private :: upper  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(len=10), parameter, private :: number = '0123456789'

  !
  ! known isotopes of select elements and electron
  ! weights are in AMU for neutral atom
  ! https://en.wikipedia.org/wiki/Isotopes_of_hydrogen
  !
  integer, parameter, private :: max_isotopes = 8

  type (iso), dimension(max_isotopes), parameter, private :: isox = (/ &
       iso('e', 'Electron', -1, 0.000548579909_fp), &   ! 1
       iso('H', 'Hydrogen',  1, 1.00782503207_fp),  &   ! 2
       iso('D', 'Deuterium', 1, 2.0141017778_fp),   &   ! 3
       iso('T', 'Tritium',   1, 3.0160492777_fp),   &   ! 4
       iso('',  '',          2, 3.016029_fp),       &   ! 5  He3
       iso('',  '',          2, 4.002602_fp),       &   ! 6  He4
       iso('',  '',          3, 6.015122_fp),       &   ! 7  Li6
       iso('',  '',          3, 7.016004_fp)        &   ! 8  Li7
       /)

  ! 
  ! element* arrays indexed by atomic number
  !
  integer, dimension(max_elements), private :: element_index = (/ &  ! index into elex array by atomic number
       39,  40,  &
       51,  11,   9,  16,  59,  67,  32,  63,  &
       60,  55,   3,  90,  69,  85,  21,   5,  &
       48,  17,  87,  99, 103,  24,  56,  33,  &
       23,  64,  26, 108,  36,  38,   6,  88,  15,  49,  &
       79,  93, 106, 109,  61,  57,  96,  84,  &
       82,  72,   2,  18,  46,  92,  86,  97,  45, 105,  &
       25,  10,  50,  19,  75,  62,  73,  91,  &
       31,  37,  95,  28,  43,  29, 101, 107,  &
       53,  41,  94, 104,  80,  68,  47,  76,   8,  42, &
       100,  71,  13,  74,   7,  83,  &
       35,  78,   1,  98,  70, 102,  66,  77,  &
       4,  22,  14,  20,  30,  34,  54,  65,  &
       52,  81,  27,  89,  12,  44,  58 &
       /)

  character(len=2), dimension(max_elements), parameter, private :: element = (/ &
       'H ', 'He', &
       'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
       'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
       'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', &
       'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
       'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', &
       'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
       'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', &
       'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
       'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
       'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
       'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', &
       'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', &
       'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' &
       /)

  character(len=13), dimension(max_elements), parameter, private :: element_name = (/ &
       'Hydrogen     ', 'Helium       ', &
       'Lithium      ', 'Beryllium    ', 'Boron        ', 'Carbon       ', &
       'Nitrogen     ', &
       'Oxygen       ', 'Fluorine     ', 'Neon         ', &
       'Sodium       ', 'Magnesium    ', 'Aluminum     ', 'Silicon      ', &
       'Phosphorus   ', 'Sulfur       ', 'Chlorine     ', 'Argon        ', &
       'Potassium    ', 'Calcium      ', 'Scandium     ', 'Titanium     ', &
       'Vanadium     ', 'Chromium     ', 'Manganese    ', 'Iron         ', &
       'Cobalt       ', 'Nickel       ', 'Copper       ', 'Zinc         ', &
       'Gallium      ', 'Germanium    ', 'Arsenic      ', 'Selenium     ', &
       'Bromine      ', 'Krypton      ', &
       'Rubidium     ', 'Strontium    ', 'Yttrium      ', 'Zirconium    ', &
       'Niobium      ', 'Molybdenum   ', 'Technetium   ', 'Ruthenium    ', &
       'Rhodium      ', 'Palladium    ', 'Silver       ', 'Cadmium      ', &
       'Indium       ', 'Tin          ', 'Antimony     ', 'Tellurium    ', &
       'Iodine       ', 'Xenon        ', &
       'Cesium       ', 'Barium       ', 'Lanthanum    ', 'Cerium       ', &
       'Praseodymium ', 'Neodymium    ', 'Promethium   ', 'Samarium     ', &
       'Europium     ', 'Gadolinium   ', 'Terbium      ', 'Dysprosium   ', &
       'Holmium      ', 'Erbium       ', 'Thulium      ', 'Ytterbium    ', &
       'Lutetium     ', 'Hafnium      ', 'Tantalum     ', 'Tungsten     ', &
       'Rhenium      ', 'Osmium       ', 'Iridium      ', 'Platinum     ', &
       'Gold         ', 'Mercury      ', &
       'Thallium     ', 'Lead         ', 'Bismuth      ', 'Polonium     ', &
       'Astatine     ', 'Radon        ', &
       'Francium     ', 'Radium       ', 'Actinium     ', 'Thorium      ', &
       'Protactinium ', 'Uranium      ', 'Neptunium    ', 'Plutonium    ', &
       'Americium    ', 'Curium       ', 'Berkelium    ', 'Californium  ', &
       'Einsteinium  ', 'Fermium      ', 'Mendeleev    ', 'Nobelium     ', &
       'Lawrencium   ', 'Rutherfordium', 'Dubnium      ', 'Seaborgium   ', &
       'Bohrium      ', 'Hassium      ', 'Meitnerium   ' &
       /)

  interface to_periodic_table
    module procedure to_periodic_table_r4, to_periodic_table_r8
  end interface to_periodic_table
  interface to_pseudo_periodic_table
    module procedure to_pseudo_periodic_table_r4, to_pseudo_periodic_table_r8
  end interface to_pseudo_periodic_table
  interface   name_periodic_table
    module procedure name_periodic_table_r4 ,name_periodic_table_r8
  end interface name_periodic_table
  interface name_pseudo_periodic_table
    module procedure name_pseudo_periodic_table_r4, name_pseudo_periodic_table_r8
  end interface name_pseudo_periodic_table
  interface inv_periodic_table
    module procedure inv_periodic_table_r4,inv_periodic_table_r8 
  end interface inv_periodic_table
  interface inv_pseudo_periodic_table
    module procedure inv_pseudo_periodic_table_r4,inv_pseudo_periodic_table_r8 
  end interface inv_pseudo_periodic_table
  interface standard_amu  
    module procedure  standard_amu_r4, standard_amu_r8
  end interface standard_amu
  interface tr_species_convert  
    module procedure  tr_species_convert_int, tr_species_convert_r8, tr_species_convert_r4
  end interface tr_species_convert

  !private elindex, use_amu
contains
  !
  !
  ! -------------------------- to_periodic_table ---------------------------
  !
  ! this function returns the periodic table symbol for some element,
  ! if there is an error there will be a ? in the first position of
  ! the returned string followed by a short message
  !
  ! examples:    to_periodic_table(6,0,-1,1)       ->  'C           '
  !              to_periodic_table(6,12.0107,-1,1) ->  'C12         '
  !              to_periodic_table(6,0,4,1)        ->  'C+4         '
  !              to_periodic_table(6,12.0107,4,1)  ->  'C12+4       '
  !              to_periodic_table(6,12.0107,14,1) ->  '?Charge>Z   '
  !
  !   z = atomic number -- must be greater than 0
  !   a = atomic weight or 0
  !   charge = charge on element or -1
  !
  function to_periodic_table_r4(z, a, charge, isymbol) result(out)
    implicit none

    integer, intent(in) :: z       ! the atomic number of the element
    real,    intent(in) :: a       ! the atomic weight of the element
    integer, intent(in) :: charge  ! the charge on the element -- ignore if <0
    integer, intent(in) :: isymbol ! symbol to use to separate element name and charge
    ! 0 -> _  this is an underscore
    ! 1 -> +
    character(len=12)       :: out     ! the symbol returned
    character(len=12)       :: swork   ! working character string
    character(len=1)        :: csymbol ! charge symbol
    real(fp) :: aa
    aa=a
    out=to_periodic_table_r8(z, aa, charge, isymbol)
  end function to_periodic_table_r4

  function to_periodic_table_r8(z, a, charge, isymbol) result(out)
    implicit none

    integer, intent(in) :: z       ! the atomic number of the element
    real(fp),    intent(in) :: a       ! the atomic weight of the element
    integer, intent(in) :: charge  ! the charge on the element -- ignore if <0
    integer, intent(in) :: isymbol ! symbol to use to separate element name and charge
    ! 0 -> _  this is an underscore
    ! 1 -> +
    character(len=12)       :: out     ! the symbol returned
    character(len=12)       :: swork   ! working character string
    character(len=1)        :: csymbol ! charge symbol

    if (z<=0) then
      out = '? Z<=0      '
      return
    end if

    if (z>max_elements) then
      out = '?Z too large'
      return
    end if

    out = element(z)   ! form the symbol
    if (a>0) then
      write(swork, '(I6)') nint(a)
      out = trim(out) // trim(adjustl(swork))
    endif


    ! add in charge state
    if (charge>=0) then
      if (charge>z) then
        out = '?Charge>Z   '
        return
      end if
      write(swork, '(I6)') charge
      if (isymbol > 0) then
        csymbol = '+'
      else
        csymbol = '_'   ! this is an underscore
      end if
      out = trim(out) // csymbol // trim(adjustl(swork))
    end if
  end function to_periodic_table_r8

  !
  ! ------------------------- inv_periodic_table -------------------
  !
  ! this subroutine maps a character string returned by the function
  ! to_periodic_table to the atomic number, atomic weight and charge
  ! or zero if the data is unavailable.
  !
  ! examples:   inv_periodic_table('C',    .true.,  z, a, c) -> z=6,  a=12.0110, c=-1
  !             inv_periodic_table('C',    .false., z, a, c) -> z=6,  a=0.,      c=-1
  !             inv_periodic_table('Mg',   .true.,  z, a, c) -> z=12, a=24.3050, c=-1
  !             inv_periodic_table('C12',  .true.,  z, a, c) -> z=6,  a=12.,     c=-1
  !             inv_periodic_table('C+4',  .true.,  z, a, c) -> z=6,  a=12.0110, c=4
  !             inv_periodic_table('C14+4',.true.,  z, a, c) -> z=6,  a=14,      c=4
  !             inv_periodic_table('C14_4',.true.,  z, a, c) -> z=6,  a=14,      c=4
  !             inv_periodic_table('Q+7',  .true.,  z, a, c) -> z=0,  a=0,       c=-1
  !             inv_periodic_table('',     .true.,  z, a, c) -> z=0,  a=0,       c=-1
  !
  subroutine inv_periodic_table_r4(symbol, default, z, a, charge)
    character(len=*), intent(in) :: symbol   ! string containing element (e.g. 'C12+4 ')
    ! must be left justified
    logical, intent(in)       :: default  ! .true. to return an atomic weight (isotopic average) even
    ! if data is not contained in the string

    integer, intent(out) :: z             ! atomic number of the element
    real,    intent(out) :: a             ! atomic weight of the element in amu
    integer, intent(out) :: charge        ! charge state of the element

    character(len=2) :: as     ! atomic symbol in string
    character(len=10) :: zread ! integer read buffer

    integer :: il         ! length of input string without trailing whitespace
    integer :: is         ! length of atomic symbol part in string
    integer :: j1,j2,j    ! temporary indices into element array

    il = len_trim(symbol)     ! length of string
    if (il == 0) goto 10      ! huh?

    ! --- get 2 character symbol id ---
    is = min(il,2)                                  ! length of symbol part
    as = symbol(1:is)                               ! atomic symbol
    if (index(upper,symbol(1:1)) == 0) goto 10      ! first letter is not upper case
    if (il<2 .or. index(lower,symbol(is:is)) == 0) then  ! lowercase second character?
      is=1           ! symbol has only one character
      as(2:2) = ' '  ! make sure padded with space
    end if

    ! --- search for element ---
    j1 = 1
    j2 = max_elements

    do while ((j2-j1)>1)            ! stop when j2=j1+1
      j = (j1+j2)/2                ! binary section search
      if (as<elex(j)%symbol) then  ! alphabetic ordering
        j2=j
      else
        j1=j
        if (as == elex(j1)%symbol) exit  ! exact match
      end if
    end do
    if (as == elex(j1)%symbol) then        ! matches j1,j2 or nothing
      j = j1
    else if (as == elex(j2)%symbol) then
      j = j2
    else
      goto 10     ! not found in list of elements
    end if

    ! --- found element, set defaults ---
    z = elex(j)%number
    charge = -1
    if (default) then
      a = elex(j)%weight   ! isotopic average
    else
      a = 0.
    end if

    ! --- get atomic weight ---
    j1 = is+1        ! first character of atomic weight
    j2 = j1          ! first character after atomic weight
    do while(j2<=il)
      if (index(number,symbol(j2:j2)) == 0) exit  ! character is not a number
      j2 = j2+1
    end do
    if (j1<j2) then
      zread=' '
      zread(10-(j2-j1)+1:10)=symbol(j1:j2-1)
      read(zread,'(i10)') j               ! convert string to atomic weight
      a = j
    end if

    ! --- get charge ---
    if (j2<il) then
      if (symbol(j2:j2) == '+' .or. symbol(j2:j2) == '_') then
        j1 = j2+1  ! first character of charge
        j2 = j1    ! first character after charge
        do while(j2<=il)
          if (index(number,symbol(j2:j2)) == 0) exit ! not a number
          j2 = j2+1
        end do
        if (j1<j2) then
          zread=' '
          zread(10-(j2-j1)+1:10)=symbol(j1:j2-1)
          read(zread,'(i10)') charge
        end if
      end if
    end if

    return  ! ok

    ! ----- error exit ------
10  continue
    z = 0
    a = 0.
    charge = -1

  end subroutine inv_periodic_table_r4

  subroutine inv_periodic_table_r8(symbol, default, z, a, charge)
    character(len=*), intent(in) :: symbol   ! string containing element (e.g. 'C12+4 ')
    ! must be left justified
    logical, intent(in)       :: default  ! .true. to return an atomic weight (isotopic average) even
    ! if data is not contained in the string

    integer, intent(out) :: z             ! atomic number of the element
    real(fp),    intent(out) :: a             ! atomic weight of the element
    integer, intent(out) :: charge        ! charge state of the element

    character(len=2) :: as      ! atomic symbol in string
    character(len=10) :: zread  ! integer read buffer

    integer :: il         ! length of input string without trailing whitespace
    integer :: is         ! length of atomic symbol part in string
    integer :: j1,j2,j    ! temporary indices into element array

    il = len_trim(symbol)     ! length of string
    if (il == 0) goto 10      ! huh?

    ! --- get 2 character symbol id ---
    is = min(il,2)                                  ! length of symbol part
    as = symbol(1:is)                               ! atomic symbol
    if (index(upper,symbol(1:1)) == 0) goto 10      ! first letter is not upper case
    if (il<2 .or. index(lower,symbol(is:is)) == 0) then  ! lowercase second character?
      is=1           ! symbol has only one character
      as(2:2) = ' '  ! make sure padded with space
    end if

    ! --- search for element ---
    j1 = 1
    j2 = max_elements
    do while ((j2-j1)>1)            ! stop when j2=j1+1
      j = (j1+j2)/2                ! binary section search
      if (as<elex(j)%symbol) then  ! alphabetic ordering
        j2=j
      else
        j1=j
        if (as == elex(j1)%symbol) exit  ! exact match
      end if
    end do
    if (as == elex(j1)%symbol) then        ! matches j1,j2 or nothing
      j = j1
    else if (as == elex(j2)%symbol) then
      j = j2
    else
      goto 10     ! not found in list of elements
    end if

    ! --- found element, set defaults ---
    z = elex(j)%number
    charge = -1
    if (default) then
      a = elex(j)%weight    ! isotopic average weight
    else
      a = 0._fp
    end if

    ! --- get atomic weight ---
    j1 = is+1        ! first character of atomic weight
    j2 = j1          ! first character after atomic weight
    do while(j2<=il)
      if (index(number,symbol(j2:j2)) == 0) exit  ! character is not a number
      j2 = j2+1
    end do
    if (j1<j2) then
      zread=' '
      zread(10-(j2-j1)+1:10)=symbol(j1:j2-1)
      read(zread,'(i10)') j               ! convert string to atomic weight
      a = j
    end if

    ! --- get charge ---
    if (j2<il) then
      if (symbol(j2:j2) == '+' .or. symbol(j2:j2) == '_') then
        j1 = j2+1  ! first character of charge
        j2 = j1    ! first character after charge
        do while(j2<=il)
          if (index(number,symbol(j2:j2)) == 0) exit ! not a number
          j2 = j2+1
        end do
        if (j1<j2) then
          zread=' '
          zread(10-(j2-j1)+1:10)=symbol(j1:j2-1)
          read(zread,'(i10)') charge
        end if
      end if
    end if

    return  ! ok

    ! ----- error exit ------
10  continue
    z = 0
    a = 0._fp
    charge = -1

  end subroutine inv_periodic_table_r8

  !
  !
  ! -------------------------- to_pseudo_periodic_table ---------------------------
  ! return the periodic table symbol except the isotopes of hydrogen
  ! are returned as H,D,T and an electron (z<0) is returned as e if ipseudo/=0
  !
  function to_pseudo_periodic_table_r4(ipseudo, z, a, charge, isymbol) result(out)
    implicit none

    integer, intent(in) :: ipseudo ! 0 -> standard periodic table
    ! 1  -> return H,D,T for hydrogen isotopes
    integer, intent(in) :: z       ! the atomic number of the element
    real,    intent(in) :: a       ! the atomic weight of the element
    integer, intent(in) :: charge  ! the charge on the element
    integer, intent(in) :: isymbol ! symbol to use to separate element name and charge
    ! 0 -> _  this is an underscore
    ! 1 -> +
    character(len=12)        :: cout    ! the symbol returned from to_periodic_table
    character(len=12)        :: out     ! the symbol returned
    character(len=12)        :: swork   ! working character string
    real(fp) :: aa
    aa=a
    out=to_pseudo_periodic_table_r8(ipseudo, z, aa, charge, isymbol)
  end function to_pseudo_periodic_table_r4

  function to_pseudo_periodic_table_r8(ipseudo, z, a, charge, isymbol) result(out)
    implicit none

    integer, intent(in) :: ipseudo ! 0 -> standard periodic table
    ! 1  -> return H,D,T for hydrogen isotopes
    integer, intent(in) :: z       ! the atomic number of the element
    real(fp),    intent(in) :: a       ! the atomic weight of the element
    integer, intent(in) :: charge  ! the charge on the element
    integer, intent(in) :: isymbol ! symbol to use to separate element name and charge
    ! 0 -> _  this is an underscore
    ! 1 -> +
    character(len=12)        :: cout    ! the symbol returned from to_periodic_table
    character(len=12)        :: out     ! the symbol returned
    character(len=12)        :: swork   ! working character string

    if (ipseudo/=0 .and. z<0) then
      out = "e"
    else
      cout = to_periodic_table(z, a, charge, isymbol)
      if (ipseudo/=0.and.cout(1:1) == 'H') then
        if (cout(2:2) == '1') then
          out = 'H'//cout(3:12)//' '
        else if (cout(2:2) == '2') then
          out = 'D'//cout(3:12)//' '
        else if (cout(2:2) == '3') then
          out = 'T'//cout(3:12)//' '
        else
          out = cout
        end if
      else
        out = cout
      end if
    end if
  end function to_pseudo_periodic_table_r8


  !
  ! ---------------------- inv_pseudo_periodic_table --------------------
  ! return the inv_periodic_table data where the hydrogen isotopes
  ! are given by H,D,T and an electron is given by e if ipseudo/=0
  !
  subroutine inv_pseudo_periodic_table_r4(ipseudo, symbol, default, z, a, charge)
    implicit none

    integer, intent(in) :: ipseudo        ! 0 -> standard periodic table
    ! 1 -> use H,D,T for hydrogen isotopes, e for electron
    character(len=*), intent(in) :: symbol   ! string containing element (e.g. 'C12+4 ')
    ! must be left justified
    logical, intent(in)       :: default  ! .true. to return an atomic weight even
    ! if data is not contained in the string

    integer, intent(out) :: z             ! atomic number of the element
    real,    intent(out) :: a             ! atomic weight of the element in AMU
    integer, intent(out) :: charge        ! charge state of the element
    integer              :: k             ! isotope number of hydrogen

    character(len=12)           :: out,lsymbol     ! construct a new symbol string
    character(len=3), parameter :: first  = 'HDT'  ! allowed hydrogen first character
    character(len=3), parameter :: second = ' +_'  ! allowed hydrogen second character

    out = symbol              ! copy, length is at least 12
    if (ipseudo/=0) then
      lsymbol = adjustl(out)          ! left justify
      k = index(first,lsymbol(1:1))   ! look for hydrogen

      if (k/=0 .and. index(second,lsymbol(2:2))/=0) then
        if (k==1) then
          out = 'H1'//lsymbol(2:11)
        else if (k==2) then
          out = 'H2'//lsymbol(2:11)
        else
          out = 'H3'//lsymbol(2:11)
        end if
      else if (lsymbol(1:2)=='e ') then
        z = -1
        a = electron_amu
        charge = -1
        return
      end if
    end if
    call inv_periodic_table(out, default, z, a, charge)

  end subroutine inv_pseudo_periodic_table_r4

  subroutine inv_pseudo_periodic_table_r8(ipseudo, symbol, default, z, a, charge)
    implicit none

    integer, intent(in) :: ipseudo        ! 0 -> standard periodic table
    ! 1 -> use H,D,T for hydrogen isotopes, e for electron
    character(len=*), intent(in) :: symbol   ! string containing element (e.g. 'C12+4 ')
    ! must be left justified
    logical, intent(in)       :: default  ! .true. to return an atomic weight even
    ! if data is not contained in the string

    integer, intent(out) :: z             ! atomic number of the element
    real(fp),    intent(out) :: a             ! atomic weight of the element in AMU
    integer, intent(out) :: charge        ! charge state of the element
    integer              :: k             ! isotope number of hydrogen

    character(len=12)           :: out,lsymbol     ! construct a new symbol string
    character(len=3), parameter :: first  = 'HDT'  ! allowed hydrogen first character
    character(len=3), parameter :: second = ' +_'  ! allowed hydrogen second character

    out = symbol              ! copy, length is at least 12
    if (ipseudo/=0) then
      lsymbol = adjustl(out)          ! left justify
      k = index(first,lsymbol(1:1))   ! look for hydrogen

      if (k/=0 .and. index(second,lsymbol(2:2))/=0) then
        if (k==1) then
          out = 'H1'//lsymbol(2:11)
        else if (k==2) then
          out = 'H2'//lsymbol(2:11)
        else
          out = 'H3'//lsymbol(2:11)
        end if
      else if (lsymbol(1:2)=='e ') then
        z = -1
        a = electron_amu
        charge = -1
        return
      end if
    end if
    call inv_periodic_table(out, default, z, a, charge)

  end subroutine inv_pseudo_periodic_table_r8

  !
  ! ----------------------- name_periodic_table --------------
  ! return the name of the element with atomic number 'z' or
  ! a '?' if it is out of bounds -- add atomic weight to end
  ! if a>0
  ! example:
  !    name_periodic_table(20,0)     -> Calcium
  !    name_periodic_table(20,40.08) -> Calcium40
  !
  function name_periodic_table_r4(z,a) result(out)
    integer, intent(in) :: z       ! the atomic number of the element
    real,    intent(in) :: a       ! the atomic weight of the element
    character(len=16)   :: out     ! name of the element
    integer             :: i       ! length of trimmed name
    character(len=3)    :: as      ! a as a string
    real(fp) :: aa
    aa=a
    out= name_periodic_table_r8(z,aa)
  end function name_periodic_table_r4

  function name_periodic_table_r8(z,a) result(out)
    integer, intent(in)  :: z       ! the atomic number of the element
    real(fp), intent(in) :: a       ! the atomic weight of the element
    character(len=16)    :: out     ! name of the element
    integer              :: i       ! length of trimmed name
    character(len=3)     :: as      ! a as a string

    if (z>0 .and. z<=max_elements) then
      out(1:13) = element_name(z)
      out(14:16) = '   '          ! len_trim fails without this
      if (a>0 .and. a<1000) then
        i = len_trim(out)
        write(as,'(I3)') int(nint(a))
        out(i+1:i+3) = adjustl(as)
      end if
    else
      out = "?"
    end if
  end function name_periodic_table_r8


  !
  ! ----------------------- name_pseudo_periodic_table --------------
  ! return the name of the element with atomic number 'z' or
  ! a '?' if it is out of bounds but return Deuterium, Tritium
  ! as appropriate for hydrogen
  ! example:
  !    name_periodic_table(1,20,0)     -> Calcium
  !    name_periodic_table(1,20,40.08) -> Calcium40
  !    name_periodic_table(1,1,1)      -> Hydrogen
  !    name_periodic_table(1,1,2)      -> Deuterium
  !    name_periodic_table(1,1,3)      -> Tritium
  !    name_periodic_table(0,1,3)      -> Hydrogen3
  function name_pseudo_periodic_table_r4(ipseudo,z,a) result(out)
    integer, intent(in) :: ipseudo ! 0 -> standard periodic table
    ! 1 -> use H,D,T for hydrogen isotopes
    integer, intent(in) :: z       ! the atomic number of the element
    real,    intent(in) :: a       ! the atomic weight of the element
    character(len=16)   :: out     ! name of the element
    integer             :: ia      ! a as an integer
    real(fp) :: aa
    aa=a
    out=name_pseudo_periodic_table_r8(ipseudo,z,aa)
    !
  end function name_pseudo_periodic_table_r4

  function name_pseudo_periodic_table_r8(ipseudo,z,a) result(out)
    integer, intent(in) :: ipseudo ! 0 -> standard periodic table
    ! 1 -> use H,D,T for hydrogen isotopes
    integer,  intent(in) :: z       ! the atomic number of the element
    real(fp), intent(in) :: a       ! the atomic weight of the element
    character(len=16)    :: out     ! name of the element
    integer              :: ia      ! a as an integer

    if (z==1 .and. ipseudo>0) then
      ia = nint(a)
      if (ia==1) then
        out = 'Hydrogen'
      else if (ia==2) then
        out = 'Deuterium'
      else
        out = 'Tritium'
      end if
    else
      out = name_periodic_table_r8(z,a)
    end if
  end function name_pseudo_periodic_table_r8

  !
  ! -------------------- elindex -----------------
  ! private function returns the index into elex by atomic number or 0 if not found
  !
  function elindex(z) result(ii)
    integer, intent(in) :: z    ! atomic number
    integer             :: ii   ! index into elex

    integer :: j,kk  ! temp

    ii = 0
    if (z<1 .or. z>max_elements) return     ! z out of range

    kk = element_index(z)
    if (kk>=1 .and. kk<=max_elements) then
      if (elex(kk)%number==z) then         ! already found
        ii = kk
        return
      end if
    end if

    ! search
    do j=1, max_elements
      if (elex(j)%number==z) then
        ii = j
        element_index(z) = j
        return
      end if
    end do
  end function elindex

  !
  ! ------------------------- standard_amu -------------------
  !
  ! this subroutine maps a charge number (z) to the corresponding
  !  atom's standard mass in amu, A; A=0.0 is returned if z is out
  !  of range.
  !
  ! examples:   call standard_amu(6,A)  --> A=12.011  (carbon)
  !             call standard_amu(26,A) --> A=55.847  (iron)
  !
  subroutine standard_amu_r4(z,a) 
    integer, intent(in) :: z     ! atomic number (in)
    real, intent(out) :: a       ! standard atomic mass, amu (out)

    ! -----------------

    real(fp) :: a8

    ! -----------------
    call standard_amu_r8(z,a8)
    a = real(a8)
  end subroutine standard_amu_r4

  subroutine standard_amu_r8(z,a) 
    integer, intent(in) :: z     ! atomic number (in)
    real(fp), intent(out) :: a     ! standard atomic mass, amu (out)

    integer :: kk 

    kk = elindex(z)
    if (kk<1 .or. kk>max_elements) then
      a = 0._fp
    else
      a = elex(kk)%weight
    end if
  end subroutine standard_amu_r8


  !
  ! -------------------- atom_amu_to_proton -------------------
  ! return the mass of an atom in protons -- account for electrons as if it makes a difference
  ! the only thing missing is the binding energy which would need to be added to the result
  !
  function atom_amu_to_proton(z,ich,a) result(p)
    integer, intent(in) :: z    ! atomic number
    integer, intent(in) :: ich  ! charge state of atom
    real(fp),intent(in) :: a    ! atomic weight in amu of neutral atom

    real(fp) :: p               ! mass in protons

    if (ignore_electron_AMU) then
      p = max(0._fp, a/proton_amu)
    else
      p = max(0._fp, (a-max(0,min(ich,z))*electron_amu)/proton_amu)
    end if

    if (z==1 .and. ich==1) then
      p = p + ((13.6_fp*electronCharge)/speedLight**2)/protonMass  ! add hydrogen binding energy
      if (abs(p-1._fp)<2.e-3_fp) then
        p=1.0_fp   ! so hydrogen comes out exactly as 1, would be 1.000074 using H=1.0079 (isotope avg)
        !          or the more physical 0.9999999854548374 using the more accurate H1= 1.00782503207
        !          and finally get      0.9999999999474358 when the binding energy is added
      end if
    end if
  end function atom_amu_to_proton


  !
  ! --------------------- use_amu -------------------
  ! convert the user requested amu to the one which should be used,
  ! set argument a=0. to get the istopic average amu for the element,
  ! return 0. on an error
  !
  function use_amu(z,a) result(u)
    integer,  intent(in) :: z    ! atomic number
    real(fp), intent(in) :: a    ! atomic weight in amu as requested

    real(fp) :: u    ! atomic weight to use

    integer :: kk   ! index into elex array
    integer :: ii   ! index into isox array of first isotope of element
    integer :: i    ! isox index

    real(fp) :: adiff      ! difference between user amu and isotope amu
    real(fp) :: adiff_min  ! minimum adiff encountered so far

    u = 0._fp    ! default in case of error

    if (a<(z*proton_amu-.1_fp)) then
      call standard_amu_r8(z,u)    ! user amu is too small, use isotopic average
    else
      kk = elindex(z)
      if (kk>=1 .and. kk<=max_elements) then
        ii = elex(kk)%isoindex
        if (ii>0) then
          adiff_min = 1.e20_fp
          do i = ii, max_isotopes
            if (isox(i)%number/=z) exit

            adiff = abs(a-isox(i)%weight)

            if (adiff<adiff_min) then
              u = isox(i)%weight         ! best match so far between a and known isotopes
              adiff_min = adiff
            end if
          end do
        else
          if (ignore_user_AMU) then
            call standard_amu_r8(z,u)     ! ignore user AMU, use isotopic average
          else
            u = max(z*proton_amu,a)       ! element does not have known isotopes so use user supplied amu
          end if
        end if
      end if
    end if
  end function use_amu


  SUBROUTINE tr_species_convert_int(iZatom, iZcharge, iAMU, qatom, qcharge, mass, &
       ierr)

    !--------------------------------
    ! D. McCune Dec. 2007
    !    Convert "integer" species specification
    !      {atomic number, ionic charge number, atomic weight in AMU}
    !    to "floating point" representation
    !      {atomic charge (C), ionic charge (C), atomic mass in kg}
    !
    !    If iZatom=-1 provide output quantites for electrons.
    !
    ! ...to be used during initialization of plasma state system at the
    ! start of a simulation...
    !
    !--------------------------------

    integer, intent(in) :: iZatom    ! atomic number: 1 for H, 2 for He, ...
    !  special: -1 for electrons

    integer, intent(in) :: iZcharge  ! ionic charge btw 1 and iZatom
    !  ignored for electrons

    integer, intent(in) :: iAMU      ! atomic weight in AMU, nearest integer

    real(fp), intent(out) :: qatom    ! charge of fully stripped ion
    real(fp), intent(out) :: qcharge  ! actual charge
    real(fp), intent(out) :: mass     ! mass in kg

    real(fp),parameter::  ps_epslon  = 1.0e-34_fp       ! small number
    real(fp),parameter::  ps_epsinv  = 1.0e+34_fp       ! large number
    integer, intent(out) :: ierr     ! status code, 0=OK

    !---------------------------
    integer  :: iq,iqc
    real(fp) :: a, ps_me, ps_mp, ps_xe
    !---------------------------

    ierr = 0

    qatom   = 0.0_fp
    qcharge = 0.0_fp
    mass    = 0.0_fp

    if (ps_params) then
      ps_me  = 9.1094e-31_fp        ! electron mass, kg
      ps_mp  = 1.6726e-27_fp        ! proton mass, kg
      ps_xe  = 1.6022e-19_fp        ! Joule/ptcl-eV (also, electron charge)
    else
      ps_me  = protonMass*(electron_amu/proton_amu)
      ps_mp  = protonMass
      ps_xe  = electronCharge
    end if

    if(iZatom.eq.-1) then
      qatom   = -ps_xe
      qcharge = -ps_xe
      mass    =  ps_me

    else if((iZatom.eq.0).or.(iZatom.lt.-1)) then
      ierr = 1
      ! ' ?tr_species_convert: atomic number iZatom out of range low:', iZatom

    else
      a = use_amu(iZatom, real(iAmu,fp))

      if (a<0.5_fp) then 
        ierr = 2
        !' ?tr_species_convert: atomic number iZatom out of range high:', iZatom

      else
        iq  = iZatom
        iqc = min(iq,iZcharge)

        qatom   = iq*ps_xe
        qcharge = iqc*ps_xe
        mass    = atom_amu_to_proton(iq,iqc,a)*ps_mp
      endif

    endif
  END SUBROUTINE tr_species_convert_int

  SUBROUTINE tr_species_convert_r4(iZatom, iZcharge, zAMU, qatom, qcharge, mass, &
       ierr)

    !--------------------------------
    ! D. McCune Dec. 2007
    !    Convert "integer" species specification
    !      {atomic number, ionic charge number, atomic weight in AMU}
    !    to "floating point" representation
    !      {atomic charge (C), ionic charge (C), atomic mass in kg}
    !
    !    If iZatom=-1 provide output quantites for electrons.
    !
    ! ...to be used during initialization of plasma state system at the
    ! start of a simulation...
    !
    !--------------------------------

    integer, intent(in) :: iZatom    ! atomic number: 1 for H, 2 for He, ...
    !  special: -1 for electrons

    integer, intent(in) :: iZcharge  ! ionic charge btw 1 and iZatom
    !  ignored for electrons

    REAL, intent(in) :: zAMU      ! atomic weight in AMU, nearest integer

    real(fp), intent(out) :: qatom    ! charge of fully stripped ion
    real(fp), intent(out) :: qcharge  ! actual charge
    real(fp), intent(out) :: mass     ! mass in kg

    real(fp),parameter::  ps_epslon  = 1.0e-34_fp       ! small number
    real(fp),parameter::  ps_epsinv  = 1.0e+34_fp       ! large number
    integer, intent(out) :: ierr     ! status code, 0=OK

    !---------------------------
    integer :: iq,iqc
    REAL(fp) :: a, ps_me, ps_mp, ps_xe
    !---------------------------

    ierr = 0

    qatom   = 0.0_fp
    qcharge = 0.0_fp
    mass    = 0.0_fp

    if (ps_params) then
      ps_me  = 9.1094e-31_fp        ! electron mass, kg
      ps_mp  = 1.6726e-27_fp        ! proton mass, kg
      ps_xe  = 1.6022e-19_fp        ! Joule/ptcl-eV (also, electron charge)
    else
      ps_me  = protonMass*(electron_amu/proton_amu)
      ps_mp  = protonMass
      ps_xe  = electronCharge
    end if

    if(iZatom.eq.-1) then
      qatom   = -ps_xe
      qcharge = -ps_xe
      mass    =  ps_me

    else if((iZatom.eq.0).or.(iZatom.lt.-1)) then
      ierr = 1
      ! ' ?tr_species_convert: atomic number iZatom out of range low:', iZatom

    else
      a = use_amu(iZatom, real(zAmu,fp))

      if (a<0.5_fp) then 
        ierr = 2
        !' ?tr_species_convert: atomic number iZatom out of range high:', iZatom

      else
        iq  = iZatom
        iqc = min(iq,iZcharge)

        qatom   = iq*ps_xe
        qcharge = iqc*ps_xe
        mass    = atom_amu_to_proton(iq,iqc,a)*ps_mp
      endif

    endif
  END SUBROUTINE tr_species_convert_r4

  SUBROUTINE tr_species_convert_r8(iZatom, iZcharge, zAMU, qatom, qcharge, mass, &
       ierr)

    !--------------------------------
    ! D. McCune Dec. 2007
    !    Convert "integer" species specification
    !      {atomic number, ionic charge number, atomic weight in AMU}
    !    to "floating point" representation
    !      {atomic charge (C), ionic charge (C), atomic mass in kg}
    !
    !    If iZatom=-1 provide output quantites for electrons.
    !
    ! ...to be used during initialization of plasma state system at the
    ! start of a simulation...
    !
    !--------------------------------

    integer, intent(in) :: iZatom    ! atomic number: 1 for H, 2 for He, ...
    !  special: -1 for electrons

    integer, intent(in) :: iZcharge  ! ionic charge btw 1 and iZatom
    !  ignored for electrons

    real(fp), intent(in) :: zAMU      ! atomic weight in AMU, nearest integer

    real(fp), intent(out) :: qatom    ! charge of fully stripped ion
    real(fp), intent(out) :: qcharge  ! actual charge
    real(fp), intent(out) :: mass     ! mass in kg

    real(fp), parameter::  ps_epslon  = 1.0e-34_fp       ! small number
    real(fp), parameter::  ps_epsinv  = 1.0e+34_fp       ! large number
    integer, intent(out) :: ierr     ! status code, 0=OK

    !---------------------------
    integer  :: iq,iqc
    real(fp) :: a, ps_me, ps_mp, ps_xe
    !---------------------------

    ierr = 0

    qatom   = 0.0_fp
    qcharge = 0.0_fp
    mass    = 0.0_fp

    if (ps_params) then
      ps_me  = 9.1094e-31_fp        ! electron mass, kg
      ps_mp  = 1.6726e-27_fp        ! proton mass, kg
      ps_xe  = 1.6022e-19_fp        ! Joule/ptcl-eV (also, electron charge)
    else
      ps_me  = protonMass*(electron_amu/proton_amu)
      ps_mp  = protonMass
      ps_xe  = electronCharge
    end if

    if(iZatom.eq.-1) then
      qatom   = -ps_xe
      qcharge = -ps_xe
      mass    =  ps_me

    else if((iZatom.eq.0).or.(iZatom.lt.-1)) then
      ierr = 1
      ! ' ?tr_species_convert: atomic number iZatom out of range low:', iZatom

    else
      a = use_amu(iZatom, zAMU)

      if (a<0.5_fp) then 
        ierr = 2
        !' ?tr_species_convert: atomic number iZatom out of range high:', iZatom

      else
        iq  = iZatom
        iqc = min(iq,iZcharge)

        qatom   = iq*ps_xe
        qcharge = iqc*ps_xe
        mass    = atom_amu_to_proton(iq,iqc,a)*ps_mp

        !if (.false. .and. iq==1 .and. a>1.5_fp .and. a<2.5_fp) then   ! deuterium
        !   print *, 'was = ', mass/1.6726E-27_fp
        !   mass = 1.99831530905844*1.6726E-27_fp        ! old periodic_table
        !   mass = 2.0*1.6726E-27_fp                     ! pre periodic_table
        !   mass = 2.014102*1.6726E-27_fp                ! deuterium neutral atom
        !   mass = 1.9990077004379938_fp*1.6726E-27_fp   ! newer periodic_table
        !   print *, 'now = ', mass/1.6726E-27_fp
        !end if
      endif

    endif
  END SUBROUTINE tr_species_convert_r8
end module periodic_table_mod
 
