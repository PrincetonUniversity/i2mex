subroutine trx_gen_lh_files(ss, fpath,runid,ierr)

  ! (called from trx_gen_state) -- generate files for driving LH codes
  ! when Plasma State derived data is available

  ! The plasma state (SS) has been loaded!
  ! SPLITN access to TRANSP namelist variables is also available.

  use plasma_state_mod
  use trx_module

  implicit NONE

  type (plasma_state) :: ss           ! state object...

  character*(*), intent(in) :: fpath  ! directory in which to write files
  character*(*), intent(in) :: runid  ! TRANSP runid
  integer, intent(out) :: ierr   ! completion status, 0=OK

  !-----------------------
  logical :: dirsw
  character*512 :: zdsave

  integer :: inuma  ! number of antennas or sources
  integer :: lunzer,iant

  character*6 :: ichant
  character*6, dimension(:), allocatable :: ichaa ! ascii encoded antenna #s
  !-----------------------

  ierr = 0

  ! # antennas; form ascii encoded list of antenna numbers

  inuma=ss%nlhrf_src
  if(inuma.le.0) then
     write(lunzer(0),*) ' %trx_gen_lh_files: no LH sources found.'
     return
  endif

  allocate(ichaa(inuma))

  do iant=1,inuma
     write(ichant,'(i6)') iant
     ichaa(iant)=trim(adjustl(ichant))
  enddo

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

  !----------------------------------------------------
  if(dirsw) then
     ! cd back to original directory
     call sset_cwd(zdsave,ierr)
  endif

end subroutine trx_gen_lh_files
