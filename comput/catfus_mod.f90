module catfus_mod
  ! Fusion reaction indexing data
  ! the categorization comes originally from TRANSP at PPPL
  ! (http://w3.pppl.gov/transp)
 
  integer, parameter :: nreact=16    ! no. of indexed reactions
 
  ! reagents isotope species:
 
  character(len=3), parameter, dimension(nreact) :: r1 = &
       (/ 'D  ','D  ','D  ','D  ','T  ','T  ','T  ','HE3', &
       'T  ','HE3','T  ','HE3','D  ','D  ','D  ','D  '   /)
 
  character(len=3), parameter, dimension(nreact) :: r2 = &
       (/ 'T  ','HE3','D:P','D:N','T  ','HE3','D  ','D  ', &
       'D  ','D  ','D  ','D  ','T  ','HE3','D:P','D:N'   /)
 
  ! reagent categories:
  ! 0->thermal, 1->beam, 2->fusion product, 3->RF tail
 
  ! valid reagent1 categories, by reaction index  (lower indicies 
  ! allow beam or thermal;
  ! higher indices are special purpose, generally speaking).
 
  integer, parameter, dimension(nreact) :: typmin = &
       (/  0,  0,    0,    0,    0,  0,    1,  1,   &
       2,  2,    3,  3,    3,  3,    3,  3    /)
 
  integer, parameter, dimension(nreact) :: typmax = &
       (/  1,  1,    1,    1,    1,  1,    1,  1,   &
       2,  2,    3,  3,    3,  3,    3,  3    /)
 
  ! maximum valid reagent2 categories (1 means beam or thermal; 
  ! 0 means thermal only).
 
  integer, parameter, dimension(nreact) :: typ2max = &
       (/  1,  1,    1,    1,    1,  1,    0,  0,   &
       1,  1,    1,  1,    1,  1,    1,  1    /)
 
end module catfus_mod
