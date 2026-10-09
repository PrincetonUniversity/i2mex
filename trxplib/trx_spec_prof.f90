subroutine trx_spec_prof(maxn,n_species,iorder,zident,id_array,ierr)
!
!  set up the indicated profiles:
!    zident.eq.'N'  --  densities, for all plasma species
!    zident.eq.'T'  --  temperatures, for all plasma species
!    zident.eq.'<EPERP>' -- average perpendicular energy / particle, all...
!    zident.eq.'<EPLL>' -- average parallel energy / particle, all...
!
!  it is allowed for some or all of the profiles to be already defined.
!
  implicit NONE
!
  integer, intent(in) :: maxn           ! max no. of species expected
  integer, intent(out) :: n_species     ! actual no. of species found
  integer, intent(in) :: iorder         ! fit order (e.g. 1 for Hermite)
  character*(*), intent(in) :: zident   ! "N" or "T" or "<EPERP>" or "<EPLL>"
!    (case insensitive test is used)
!
  integer, intent(out) :: id_array(maxn)  ! ids of indicated profiles
  integer, intent(out) :: ierr          ! completion code, 0=OK
!
!----------------------------------
!
  character*20, dimension(:), allocatable :: slbl
  character*10, dimension(:,:), allocatable :: abray
  integer, dimension(:), allocatable :: ifast,itype,izc
  real, dimension(:), allocatable :: aar4,zzr4
  integer :: i,id=0,idum,iwarn,lunzer
  real*8 factor,offset
  character*10 ztest
  character*16 zunits
  character*20 zname
  logical :: iexist
!
!----------------------------------
!  check number of species available...
!
  id_array=0
  call trx_nspec(n_species)
  if(n_species.eq.0) then
     ierr=1
     return
  endif
  if(n_species.gt.maxn) then
     write(lunzer(0),*) '?trx_spec_prof:  maxn=',maxn,' but: n_species=', &
          n_species
     n_species=0
     ierr=1
     return
  endif
!
!  allocate arrays, get species information, with trread call...
!
  allocate(slbl(n_species))
  allocate(abray(4,n_species))
  allocate(ifast(n_species),itype(n_species),izc(n_species))
  allocate(aar4(n_species),zzr4(n_species))
!
  call rd_species(n_species,idum,slbl,abray,itype,ifast,zzr4,aar4,izc)
!
  ztest=zident
  call uupper(ztest)       ! uppercase conversion (just in case)
!
  do i=1,n_species
     if(ztest.eq.'N') then
!
!  get densities
!
        call eq_gfnum(abray(1,i),id)
        if(id.ne.0) then
           continue        ! accept the old definition
        else
           call trx_prof(abray(1,i),zunits,iorder,id,ierr)
           if(ierr.ne.0) then
              write(lunzer(0),*) '?trx_spec_prof:  trx_prof error on: ', &
                   abray(1,i)
              exit
           endif
        endif
 
     else if(ztest.eq.'T') then
!
!  get temperatures.  for ifast(i).eq.1, convert from <E>
!
        call eq_gfnum(abray(2,i),id)
        if(id.ne.0) then
           continue        ! accept the old definition
        else
           call trx_prof(abray(2,i),zunits,iorder,id,ierr)
           if(ierr.ne.0) then
              write(lunzer(0),*) '?trx_spec_prof:  trx_prof error on: ', &
                   abray(2,i)
              exit
           endif
        endif
!
!  T = 2/3<E>
!
        if(ifast(i).eq.1) then
           idum=id
           zname='T23_'//abray(2,i)
           call eq_gfnum(zname,id)
           if(id.eq.0) then
              factor=2.0d0/3.0d0
              offset=0.0d0
              call trx_lintrans(idum,factor,offset,zname,id,ierr)
              if(ierr.ne.0) then
                 write(lunzer(0),*) &
                      '?trx_spec_prof:  trx_lintrans error on: ',zname
                 exit
              endif
           endif
        endif
!
     else if(ztest.eq.'<EPERP>') then
!
!  <Eperp>
!
        if(ifast(i).eq.0) then
!  thermal population:  <Eperp>=T
           call eq_gfnum(abray(2,i),id)
           if(id.ne.0) then
              continue        ! accept the old definition
           else
              call trx_prof(abray(2,i),zunits,iorder,id,ierr)
              if(ierr.ne.0) then
                 write(lunzer(0),*) &
                      '?trx_spec_prof:  trx_prof error on: ',abray(2,i)
                 exit
              endif
           endif
        else if(ifast(i).eq.1) then
!  beam or fusion product
!  transp/rplot calculator expression => <Eperp>, EV
           zname='E_'//abray(3,i)
           call eq_gfnum(zname,id)
           if(id.eq.0) then
              call rpexist_profile(zname,iexist)
              if(.NOT.iexist) then
                 call trx_calc(zname,'Eperp','EV', &
                      abray(3,i)//'*'//abray(1,i)// &
                      '/(1.602e-19*max(1,'//abray(1,i)//'**2))', &
                      iwarn,ierr)
                 if(ierr.ne.0) then
                    write(lunzer(0),*) &
                         '?trx_spec_prof:  trx_calc error on: ',abray(3,i)
                    exit
                 endif
              endif
              call trx_prof(zname,zunits,iorder,id,ierr)
              if(ierr.ne.0) then
                 write(lunzer(0),*) &
                      '?trx_spec_prof:  trx_prof error on: '//zname
                 exit
              endif
           endif
        else if(ifast(i).eq.2) then
!  RF tail:  just use the function
           call eq_gfnum(abray(3,i),id)
           if(id.ne.0) then
              continue        ! accept the old definition
           else
              call trx_prof(abray(3,i),zunits,iorder,id,ierr)
              if(ierr.ne.0) then
                 write(lunzer(0),*) &
                      '?trx_spec_prof:  trx_prof error on: ',abray(3,i)
                 exit
              endif
           endif
        endif
     else if(ztest.eq.'<EPLL>') then
!
!  <Epll>
!
        if(ifast(i).eq.0) then
!  thermal population:  <Epll>=0.5*T
           zname='Epll_'//abray(2,i)
           call eq_gfnum(zname,id)
           if(id.ne.0) then
              continue        ! accept the old definition
           else
              call eq_gfnum(abray(2,i),id)
              if(id.eq.0) then
                 call trx_prof(abray(2,i),zunits,iorder,id,ierr)
                 if(ierr.ne.0) then
                    write(lunzer(0),*) &
                         '?trx_spec_prof: trx_prof error on: Epll_'//abray(2,i)
                    exit
                 endif
              endif
              idum=id
              factor=1.0d0/2.0d0
              offset=0.0d0
              call trx_lintrans(idum,factor,offset,zname,id,ierr)
              if(ierr.ne.0) then
                 write(lunzer(0),*) &
                      '?trx_spec_prof:  trx_lintrans error on: ',zname
                 exit
              endif
           endif
        else if(ifast(i).eq.1) then
!  beam or fusion product
!  transp/rplot calculator expression => <Epll>, EV
           zname='E_'//abray(4,i)
           call eq_gfnum(zname,id)
           if(id.eq.0) then
              call rpexist_profile(zname,iexist)
              if(.NOT.iexist) then
                 call trx_calc(zname,'Epll','EV', &
                      abray(4,i)//'*'//abray(1,i)// &
                      '/(1.602e-19*max(1,'//abray(1,i)//'**2))', &
                      iwarn,ierr)
                 if(ierr.ne.0) then
                    write(lunzer(0),*) &
                         '?trx_spec_prof:  trx_calc error on: ',abray(4,i)
                    exit
                 endif
              endif
              call trx_prof(zname,zunits,iorder,id,ierr)
              if(ierr.ne.0) then
                 write(lunzer(0),*) &
                      '?trx_spec_prof:  trx_prof error on: ',zname
                 exit
              endif
           endif
        else if(ifast(i).eq.2) then
!  RF tail:  just use the function
           call eq_gfnum(abray(4,i),id)
           if(id.ne.0) then
              continue        ! accept the old definition
           else
              call trx_prof(abray(4,i),zunits,iorder,id,ierr)
              if(ierr.ne.0) then
                 write(lunzer(0),*) &
                      '?trx_spec_prof:  trx_prof error on: ',abray(4,i)
                 exit
              endif
           endif
        endif
     else
        ierr=1
        write(lunzer(0),*) '?trx_spec_prof: ident unknown: ',zident
        exit
     endif
!
     id_array(i)=id
!
  enddo
!
  deallocate(slbl,abray)
  deallocate(ifast,itype,izc)
  deallocate(aar4,zzr4)
!
  return
!
end subroutine trx_spec_prof
