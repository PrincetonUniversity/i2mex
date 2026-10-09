 subroutine mpi_mkgroup(comm, names, start, n, ranks, n_tot, new_group, group_comm, ierr)
!
! Create an mpi communicator group
! 01/07/2012 CLF

#ifdef __MPI
 use mpi
 use mpi_proc_data
 use logmod
 implicit none

! Input
 integer, intent(in)               :: comm  ! communicator
 character(len=*), intent(in), dimension(*)  :: names    ! print name for new_group
 integer, intent(in), dimension(*) :: start ! start of ranks in array ranks  
 integer, intent(in), dimension(*) :: n     ! number of elements in new group  
 integer, intent(in), dimension(*) :: ranks ! ranks of processes in orig_group to appear in new_group
 integer, intent(in)               :: n_tot ! number of sections in group  
! Output
 integer, intent(out)       :: new_group         ! new group handle 
 integer, intent(out)       :: group_comm        ! new group communicator handle
 integer, intent(out)       :: ierr

 integer        :: ierr111=111

! private
 integer        ::  orig_group 	     ! group (handle)
 integer        ::  myrank           ! rank of original process
 integer        ::  size             ! size of original group
 integer        ::  i_incl, n_incl, i, ii
 integer        ::  ranges(4,3)
 logical        ::  member           ! true if in group

 myrank=-1
 member=.false.

! Extract the original group handle
  call MPI_COMM_GROUP(comm, orig_group, ierr)
  call ckerr('MPI_COMM_GROUP',ierr)

  call MPI_COMM_RANK(comm, myrank, ierr)
  call ckerr('MPI_COMM_RANK',ierr)
  call MPI_COMM_SIZE(comm, size, ierr)

! add processes
  ii=0
  check_if_in: do i=1,n_tot
     ii=ii+1
     if (myrank .ge. ranks(start(i)) .and. &
         myrank .le. ranks(start(i)+n(i)-1)) then
         member=.true.
         exit check_if_in
     endif
   enddo check_if_in
   i_incl = ii
   n_incl = start(ii)+n(ii)-1
   
   if (member) then
      call MPI_GROUP_INCL(orig_group, n(i_incl), &
                          ranks(start(i_incl):n_incl),new_group, ierr)
      call ckerr('MPI_GROUP_INCL',ierr)
   else
      if (n_tot .gt. 1) then
         ii=0
         do i=1,n_tot
            if (i .ne.  i_incl) then
               ii=ii+1
               ranges(1,ii)=start(i)-1
               ranges(2,ii)=start(i)+n(i)-2
               ranges(3,ii)=1
            endif
         enddo
         call MPI_GROUP_RANGE_INCL(orig_group, n_tot-1, &
                          ranges,new_group, ierr)
      else
         call MPI_GROUP_INCL(orig_group, size-n(i_incl), &
                          ranks(start(1)+n(1):size),new_group, ierr)
      endif
   endif

! Create Communicator
! Must be called by all processes in "comm"
   
   call MPI_COMM_CREATE(comm, new_group, group_comm, ierr)
   call ckerr('MPI_COMM_CREATE',ierr)

! Give it a name
!   call MPI_COMM_SET_NAME(group_comm, name, ierr)
!   call ckerr('MPI_COMM_SET_NAME',ierr)

! add to storage
   if (member) then
      call mpi_env_comm_add(group_comm, ierr, trim(names(i_incl)))
   else
      call mpi_env_comm_add(MPI_COMM_NULL, ierr, trim(names(i_incl)))
   endif
   if (ierr .ne. 0) then
      call errorLog('mpi_mkgroup: mpi_env_comm_add failed. Error=',ierr)
   endif


CONTAINS

  subroutine ckerr(str, ierr)
    character*(*), intent(in) :: str
    integer, intent(in)       :: ierr     ! error if .ne.0

    if(ierr.eq.0) return
! if error abort
     call errorLog('mpi_mkgroup: '//trim(str)//' ierr=',ierr)
     call errorLog('failed to create group '//names(i_incl)// &
            ' issues MPI_ABORT(...) I am rank',myrank)
    call MPI_ABORT(comm,ierr,ierr111)
  end subroutine ckerr
#endif

 end subroutine mpi_mkgroup
