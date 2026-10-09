subroutine rd_nspecies(n_species,n_thi,n_thx,n_bi,n_rfi,n_fusi)
!
!  return the number of species in the (currently open) TRANSP run
!  return n_species=0 if no run is currently open.
!  return n_species=-1 if some other error occurred (rare).
!
  implicit NONE
!
  integer, intent(out) :: n_species  ! total no. of species including electrons
  integer, intent(out) :: n_thi      ! no. of "non-impurity" therm. ion species
  integer, intent(out) :: n_thx      ! no. of "impurity" thermal ion species
  integer, intent(out) :: n_bi       ! no. of "beam" ion species
  integer, intent(out) :: n_rfi      ! no. of "RF tail" ion species
  integer, intent(out) :: n_fusi     ! no. of "fusion product" ion species.
!
!-------------------------
!
  integer, parameter :: maxabr = 200
  character*10 abrdum(4,maxabr)
  real zz(maxabr)
  real amu(maxabr)
  integer izc(maxabr)
!
  character*10 abr
  character*64 zlabel
  character*32 zunits
!
  integer imulti,istype,inum
!
!----------------------------------
!
  n_species=0
  n_thi=0
  n_thx=0
  n_bi=0
  n_rfi=0
  n_fusi=0
!
  abr='NE'
  call rplabel(abr,zlabel,zunits,imulti,istype)
  if((imulti.ne.0).or.(istype.le.0)) then
     call zermsg('?rd_nspecies: no "NE" profile found.')
     return         ! NE profile not found...
  endif
!
  abr='TE'
  call rplabel(abr,zlabel,zunits,imulti,istype)
  if((imulti.ne.0).or.(istype.le.0)) then
     call zermsg('?rd_nspecies: no "TE" profile found.')
     return         ! NE profile not found...
  endif
!
!  there are electrons...
!
  n_species=1
!
!  get list of non-impurity thermal species
!
  call rd_thspec(maxabr,abrdum,zz,amu,izc,n_thi)
  if(n_thi.le.0) then                  ! expect at least one...
     n_species=-1
     return
  else
     n_species=n_species+n_thi
  endif
!
!  get list of impurity thermal species
!
  call rd_thxspec(maxabr,abrdum,zz,amu,izc,n_thx)
  if(n_thx.lt.0) then                  ! zero or more...
     n_species=-1
     return
  else
     n_species=n_species+n_thx
  endif
!
!  get list of beam ion species
!
  call rd_bmspec(maxabr,abrdum,zz,amu,izc,n_bi)
  if(n_bi.lt.0) then                  ! zero or more...
     n_species=-1
     return
  else
     n_species=n_species+n_bi
  endif
!
!  get list of RF tail ion species
!
  call rd_rfspec(maxabr,abrdum,zz,amu,izc,n_rfi)
  if(n_rfi.lt.0) then                  ! zero or more...
     n_species=-1
     return
  else
     n_species=n_species+n_rfi
  endif
!
!  get list of fusion product ion species
!
  call rd_fuspec(maxabr,abrdum,zz,amu,izc,n_fusi)
  if(n_fusi.lt.0) then                  ! zero or more...
     n_species=-1
     return
  else
     n_species=n_species+n_fusi
  endif
!
  return
end subroutine rd_nspecies
