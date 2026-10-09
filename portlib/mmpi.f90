
   module mmpi
!. mmpi                              basic MPI quantities
!
!
!  Author                            jim.conboy@jet.uk
!  Version                           1.00,  28Nov2011
!  Modifications
!
!  1.00  28Nov2011                   Jim.Conboy@ccfe.ac.uk
!                                    Written
!----_^---------------------------------=========================================|
!
!
   use logmod
#ifdef __MPI
   use mpi
#endif
   implicit  none
   private
   public  :: mmpi_init         &
             ,mmpi_end          &
             ,cmpitag
!
   public  :: mpi_id            &
             ,mpi_np            &
             ,mpi_cw            &
             ,iam_0, not_0      &
             ,mpi_run           &
             ,not_mpi
!
!..GLOBAL
!  variable                          :: name                ! description!
   integer  :: mpi_id    = 0     & ! process id
              ,mpi_np    = 0     & ! # processes     
              ,mpi_CW    = 0       ! mpi_comm_world
!
   logical  :: iam_0   = .true.  & ! process 0
              ,not_0   = .false. & ! not ..
              ,mpi_run = .false. & ! MPI active
              ,not_mpi = .true.    !     inactive
!
   character(len=8)  :: cmpi_id             ! id as char string
   contains

  subroutine mmpi_init( istat )
  !.. mmpi_init                        initialise mpi & mmpi
  !                                    codesys/source/portlib/mmpi.mmpi_init.
  !
  !
  !  Author                            Jim.Conboy@ccfe.ac.uk
  !  Version                           1.00,  28Nov2011
  !  Modifications
  !
  !  1.00  28Nov2011                   Jim.Conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..Arguments
  !
  integer, intent(out)  :: istat          !  return status, 0=OK
  !
  !..local
  integer  :: i              ! loop index 
  !
  character(len=*),parameter  :: cSr='portlib/mmpi.mmpi_init'
  !----------------------------------::================!.........................|
  !
  istat = 0
  if( mpi_np /= 0 )  return
  !
#ifdef __MPI
  call MPI_INIT( istat )
  if( istat /= 0 )  then
!     call errorlog( cSr, ' error from MPI_INIT ' )
     call bad_exit()
  endif
  !
  mpi_CW  =  MPI_COMM_WORLD
  call MPI_COMM_RANK( mpi_CW, mpi_id, istat )
  call MPI_COMM_SIZE( mpi_CW, mpi_np, i     )
  istat = max( istat, i )
  if( istat /= 0 )  then
!     call warnlog( cSr, istat,' from MPI_COMM_Rank or SIZE  ' )
  endif
  !
  iam_0 = mpi_id == 0
  not_0 = mpi_id /= 0
  mpi_run = .true.
  not_mpi = .false.
  write(cmpi_id,'(i8)') mpi_id
  cmpi_id = ' '//adjustl( cmpi_id)
  !
#else
!  call warnlog( cSr//' called but not MPI version ' )
  mpi_id = 0
  mpi_np = 1    ! it had better be..
  mpi_run = .false.  
  not_mpi = .true.
  !
#endif
  write(6,'(2a,i4,a,i4)')  &
    cSr,' -I- MPI ID ', mpi_id,' # proc =', mpi_np
  !
  end subroutine mmpi_init

  subroutine mmpi_end( istat, labort )
  !.. mmpi_end                         Terminate mpi run
  !                                    codesys/source/portlib/mmpi.mmpi_end.
  !
  !
  !  Author                            Jim.Conboy@ccfe.ac.uk
  !  Version                           1.00,  12Dec2011
  !  Modifications
  !
  !  1.00  12Dec2011                   Jim.Conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..Arguments
  !
  integer, intent(out)  :: istat          !  return status, 0=OK
  logical, intent(in), optional  :: labort         ! call mpi_abort, if .T
  !
  !..local
  integer  :: i              ! loop index 
  !
  character(len=*),parameter  :: cSr='portlib/mmpi.mmpi_end'
  !----------------------------------::================!.........................|
  !
#ifdef __MPI
  if( present( labort ))  then
     if( labort )  then
!        call warnlog( cSr, cmpi_id//' calling MPI_ABORT' )
        call MPI_ABORT( mpi_CW, istat, istat )
  return
  endif
  endif
  !
!  call infolog( cSr, cmpi_id//' calling MPI_FINALIZE ' )
  write(6,*) cSr, ' -I- ID',cmpi_id//' calling MPI_FINALIZE '
  call MPI_BARRIER( mpi_CW, istat )
  call MPI_FINALIZE( istat )
  !
#endif
  end subroutine mmpi_end


      function cMPItag( id  )
  !.. cMPItag                             Return character string ( mpi_nnnn ) for MPI 
  !                                       process id
  !                                       codesys/source/portlib/mmpi.cMPItagcMPItag.
  !
  !
  !     Author                            Jim.Conboy@ccfe.ac.uk
  !     Version                           1.00,  01Dec2011
  !     Modifications
  !
  !     1.01  21May2012                   Jim.Conboy@ccfe.ac.uk
  !                                       Change Serial tag ( to agree with csh scripts )
  !     1.00  01Dec2011                   Jim.Conboy@ccfe.ac.uk
  !                                       Written
  !----_^---------------------------------=========================================|
  !
      implicit  none
  !
  !..GLOBAL
  !--   variable                       :: name           ! description
  !
  !..Arguments
  !
      character(len=8)  :: cMPItag         ! 
      integer, intent(in)  :: id 
  !
      character(len=*),parameter  :: cFn='portlib/mmpi.cMPItag'
  !----_^------------------------------::================!.........................|
  !
#ifdef __MPI
  if( abs(id) < 1000 )  then
    write(cMPItag,'(a,i4.4)')   'mpi_',abs(id)
  else
    write(cMPItag,'(a,i6.6)')   'm_',abs(id)
  endif
  !
#else
  cMPItag = 'ser-0'
#endif
   end function cMPItag

   end module mmpi
