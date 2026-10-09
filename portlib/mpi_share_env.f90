subroutine mpi_share_env(myid_in,ierr)

  use mpi_env_mod
  implicit NONE

  !------------------------------------------------
  ! share environment information btw processes in an MPI job
  !   (non MPI version: do almost nothing)
  !------------------------------------------------

  integer, intent(in) :: myid_in   ! processor id number: 0 for "head node"
  !  if caller does not know, it should set myid_in = -1

  integer, intent(out) :: ierr     ! status code returned: 0 for OK

  !------------------------------------------------
  !  debugging option added, DMC May 22, 2009
  !  extra output of (myid.eq.0) process sees MPI_SHARE_ENV_DEBUG
  !  environment variable-- no effect except in MPI runs.
  !
  !  the communicator in this module is MPI_COMM_WORLD (no option).
  !    mpi_env_comm = MPI_COMM_WORLD
  !------------------------------------------------

#ifdef __MPI
  integer :: isizes(4)
  character, dimension(:), allocatable :: chbuf

  integer :: idbg_lines = 0
  integer, parameter :: max_dbg_lines = 100
  integer :: ilen_test,ilen_use,iolun
  character*80 dbg_lines(max_dbg_lines)
  character*40 :: ftest

  logical :: idebug
#endif
  integer :: idum
  character*30 :: dbg_val

  integer i,j,k,ilen,inumch,icwd
  integer :: myid,mpierr

  integer, save :: icalls = 0

  !------------------------------------------------
  ! clear error code (this is the only action if no __MPI):

  icalls = icalls + 1

  call get_loc_myid
  if(mdebug) write(0,*) ' proc #',loc_myid,' entering mpi_share_env, ', &
       'myid_in = ',myid_in

  ierr = 0

  myid = myid_in
  if(myid.lt.0) myid=loc_myid

  !  serial, or MPI all ranks:
  call mpi_sget_env('PPPLCLUSTER_HOME',dbg_val,idum)
  if(dbg_val.ne.' ') then
     call mpi_set_filesys('PPPLCLUSTER')
  endif

#ifdef __MPI

  if(myid.eq.0) then

     !--------------------------------------------------------
     ! do these first: it can change the isizes(1:3) counts...
     call mpi_sget_env('MPI_SHARE_ENV_DEBUG',dbg_val,idum)
     if(dbg_val.eq.'TRUE') then
        isizes(4) = 1
     else
        isizes(4) = 0
     endif

     !--------------------------------------------------------

     isizes(1)=n_alloc_names
     isizes(2)=n_active_names

     isizes(3)=0
     do i=1,n_active_names
        isizes(3) = isizes(3) + 1 + len(trim(mpi_env_names(i)))
        isizes(3) = isizes(3) + 1 + len(trim(mpi_env_vals(i)))
     enddo
  endif

  ! blocking integer MPI_BCAST:
  call intvec_bcast(isizes, 4, mpi_env_comm)

  idebug = ( isizes(4).eq.1 )
  if(idebug) then
     dbg_lines = ' '
     idbg_lines = 1
     write(dbg_lines(idbg_lines),'(1x,a,i8,a,i8)') &
          ' mpi_share_env debug output: myid = ', &
          myid,'; mpi_env_comm = ',mpi_env_comm
     idbg_lines = 2
     write(dbg_lines(idbg_lines),'(1x,a,4(1x,i5))') &
          ' 1st MPI_BCAST isizes = ',isizes
  endif

  inumch = isizes(3)

  allocate(chbuf(inumch))

  if(myid.gt.0) then
     n_alloc_names=isizes(1)
     n_active_names=isizes(2)

     if(allocated(mpi_env_names)) deallocate(mpi_env_names,mpi_env_vals)
     allocate(mpi_env_names(n_alloc_names),mpi_env_vals(n_alloc_names))
     if(mdebug) write(0,*) ' proc# ',loc_myid, &
          ' n_alloc_names,n_active_names,inumch:',&
          n_alloc_names,n_active_names,inumch

     do i=1,n_alloc_names
        mpi_env_names(i)=' '
        mpi_env_vals(i)=' '
     enddo
  endif

  if(myid.eq.0) then
     j=0
     do i=1,n_active_names
        call trans(mpi_env_names(i))
        call trans(mpi_env_vals(i))
     enddo
  endif

  call MPI_BCAST(chbuf, inumch, MPI_CHARACTER, 0, mpi_env_comm, ierr)
  if(ierr.ne.0) then
     call err('MPI_BCAST(chbuf)')
     return
  endif

  if(idebug) then
     idbg_lines = idbg_lines + 1
     write(dbg_lines(idbg_lines),'(1x,a,1x,10a1)') 'start of chbuf:', &
          chbuf(1:10)
  endif

  if(myid.gt.0) then
     j=0
     do i=1,n_active_names
        call rcv(mpi_env_names(i))
        call rcv(mpi_env_vals(i))
        if(mdebug) write(0,*) myid,loc_myid,' '//trim(mpi_env_names(i))// &
             ' '//trim(mpi_env_vals(i))
     enddo
  endif

  if(idebug) then
     if(mdebug) write(0,*) ' idebug=T in proc #',loc_myid
     do i=1,n_active_names
        if(idbg_lines.lt.max_dbg_lines) then
           idbg_lines = idbg_lines + 1
           ilen_test = 6 + len(trim(mpi_env_names(i))) + len(trim(mpi_env_vals(i)))
           if(ilen_test.le.72) then
              write(dbg_lines(idbg_lines),'(1x,a)') &
                   trim(mpi_env_names(i))//' = "'//trim(mpi_env_vals(i))//'"'
           else
              ilen_test=ilen_test-72  ! excess
              ilen_use=len(trim(mpi_env_vals(i))) - ilen_test - 3
              write(dbg_lines(idbg_lines),'(1x,a)') &
                   trim(mpi_env_names(i))//' = "'// &
                   mpi_env_vals(i)(1:ilen_use)//'..."'
           endif
        endif
     enddo
  endif
#else
  ! check mpi caller is linked with mpi_portlib 
  if (loc_myid .eq. 0 .and. myid_in .gt. 0) then
     write(0,*) ' '
     write(0,*) &
          ' >>>> portlib: non-MPI version of mpi_share_env called by processor # ',&
          loc_myid,' with myid_in =', myid_in
     write(0,*) ' >>> please use mpi_portlib library (not serial portlib) in MPI programs.'
     write(0,*) ' '
     call bad_exit
  endif
#endif

  !06/24/2011 CLF: 
  ! this section will "always" create a root directory, 
  ! even if the root directory is not writable by a user
  ! and the error is Read-only file system 
  ! 
  ! serial version also...

  if(mdebug) write(0,*) ' call mpi_env_comm_add(mpi_env_comm,...) '
  call mpi_env_comm_add(mpi_env_comm,ierr)

  icwd=0
  do i=1,n_active_names
     if(mpi_env_names(i).eq.'_cwd') icwd=i
  enddo

  if(icwd.gt.0) then
     ! go to directory specified as "cwd"
     ! make sure it exists

     call mpi_mkdir(mpi_env_comm,mpi_env_vals(icwd),ierr)
     if(ierr.ne.0) then
        call err('mpi_mkdir(cwd='//trim(mpi_env_vals(icwd))//')')
        return
     endif

     call sset_cwd(mpi_env_vals(icwd),ierr)
     if(ierr.ne.0) then
        call err('sset_cwd(cwd='//trim(mpi_env_vals(icwd))//')')
        return
     else
        if(icalls.eq.1) then
           write(0,*) ' (mpi_share_env) process myid=',myid,' cwd: ', &
                trim(mpi_env_vals(icwd))
        endif
     endif

  endif
#ifdef __MPI
  if(idebug) then
     call find_io_unit(iolun)
     ftest=' '
     write(ftest,'("mpi_share_env_debug_",i6.6,".txt")') myid
     open(unit=iolun,file=ftest,status='new',iostat=ierr)
     if(ierr.ne.0) then
        call err('mpi_share_env_debug file open')
        return
     endif
     do i=1,idbg_lines
        write(iolun,'(1x,a)') trim(dbg_lines(i))
     enddo
     close(unit=iolun)
  endif
  call MPI_BARRIER(mpi_env_comm,ierr)
#endif

  if(mdebug) write(0,*) ' proc #',loc_myid,' exiting mpi_share_env ierr=',ierr

CONTAINS

#ifdef __MPI

  subroutine trans(str)
    character*(*), intent(in) :: str

    ! load non-blank part of str into buffer
    ! side effect: j is modified... chbuf becomes concatenated set of
    !   non-blank strings terminated by nulls

    ilen=len(trim(str))
    do k=1,ilen
       j=j+1
       chbuf(j) = str(k:k)
    enddo
    j=j+1
    chbuf(j)=char(0)
  end subroutine trans

  subroutine rcv(str)
    character*(*), intent(inout) :: str

    ! fill in non-blank parts of str from buffer; done when null is found
    ! side effect: j is modified...

    integer :: icv

    k=0
    do 
       j=j+1
       icv = ichar(chbuf(j))
       if(icv.ne.0) then
          k=k+1
          str(k:k) = chbuf(j)
       else
          exit
       endif
    enddo

  end subroutine rcv

#endif

  subroutine err(errstr)

    character*(*), intent(in) :: errstr

    write(0,*) ' ?? mpi_share_env: error on process myid = ',myid
    write(0,*) '    message: '//trim(errstr)

  end subroutine err

end subroutine mpi_share_env
