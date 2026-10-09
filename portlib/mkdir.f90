  subroutine mkdir( cdir, istat )
  !.. mkdir                            Create a directory ( tree); definition may
  !                                    contain environment variables
  !                                    codesys/source/portlib/mkdir.f90
  !
  !                                    fpopen is used to expand the directory name using 
  !                                    echo, so reproducing the (Bourne) shell syntax 
  !                                    
  !  Author                            jim.conboy@ccfe.ac.uk
  !  Version                           1.00,  20May2011
  !  Modifications
  !
  !  1.00  20May2011                   jim.conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  character(len=*), intent(inout)    :: cdir           ! the directory name
  integer, intent(out)               :: istat          ! return status, 0=OK
  !
  !..local
  logical                            :: lexist
  character(len=1)                   :: char
  character(len=256)                 :: cdir_ex &
                                       ,cbuf
  !
  character(len=*),parameter         :: cSr='portlib/mkdir'
  !----------------------------------::================!.........................|
  !
  !
  call fpopen( 'echo '//cdir//char(0), cdir_ex )             ! expand string
  print *,cSr,' -I- ',trim(cdir_ex )
  !
  call fpopen( 'mkdir -vp '//trim(cdir_ex)//char(0), cbuf )  ! & create directory
  !
  !   Test for success - a test that the file is in fact a directory
  !                      would be better
  !
  inquire( file=trim(cdir_ex), exist=lexist )                ! .. did we ?
  if( lexist )  then
    istat = 0
    cdir = cdir_ex
  else
    print *,cSr,' -E- could not create ',trim(cdir)
    istat = 1
  endif
  !
  return
  end subroutine mkdir
