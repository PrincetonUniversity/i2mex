module mpi_proc_data
!. mpi_proc_data                      module to store process rank and node data on a 
!                   short list of MPI  communicators  (the serial code has an especially short list),
!                   & handle inter-process file transmission
!
!
!  Author                            D. McCune @ pppl.gov
!  Version                           1.01,  05Dec2011
!  Modifications
!
!  1.02  10Feb2012                   C.Ludescher-Furth @ pppl.gov
!                                    add comm-name to mpi_env_comm_add
!                                    interface mpi_proc_info to get comm from name
!                                    interface find_comm to get index from name
!
!  1.01  05Dec2011                   Jim.Conboy@ccfe.ac.uk
!                                    mpi_env_comm_add, mpi_bcast_list -
!                                    Add more comments, label major loops
!                                    NO functional changes
!
!  1.00  05May2011                   D. McCune @ pppl.gov
!                                    Committed, svn rev 14909
!
!----_^---------------------------------=========================================|
!
!
#ifdef __MPI
  use mpi
#endif

  implicit NONE
  PRIVATE
  SAVE
!
!..GLOBAL
!
  PUBLIC :: mpi_proc_info,mpi_env_comm_add,mpi_print_info
  PUBLIC :: mpi_proc_info_comm, mpi_proc_info_name
  PUBLIC :: mpi_set_filesys
  PUBLIC :: mpi_mkdir,mpi_file_bcast,mpi_file_transfer
  PUBLIC :: mpi_file_bcast_single,mpi_file_bcast_list
  PUBLIC :: nodepath_translator
  PUBLIC :: check_mpi_io_node, mpi_parallel_file_system

  private :: comm_arg        !  Translates comm ID to 0 ifndef _MPI
  private :: find_comm

  INTERFACE mpi_file_bcast
     module procedure mpi_file_bcast_single,mpi_file_bcast_list
  END INTERFACE

  INTERFACE mpi_proc_info
     module procedure mpi_proc_info_comm, mpi_proc_info_name
  END INTERFACE

  INTERFACE find_comm
     module procedure find_index_comm, find_index_name
  END INTERFACE


#ifdef __MPI
!  integer, parameter :: max_comm = 250
  integer, parameter :: max_comm = 20   ! CLF: dont think we need 250
#else
  integer, parameter :: max_comm = 1
#endif

  logical, parameter :: loc_debug = .FALSE.
! logical, parameter :: loc_debug = .TRUE.

  integer, parameter :: ncnode = 140       ! max characters in node name

  type :: comdat
     ! saved data for one communicator, this process
     integer :: comm        ! communicator 
                            ! if not in group: MPI_COMM_NULL (=2)
     integer :: group       ! group (handle)
     integer :: num_procs   ! #processors (size)
     integer :: myproc      ! #this process (rank in this group)
                            ! if not in group: MPI_UNDEFINED
     integer :: num_nodes   ! #nodes
     integer :: mynode      ! #this node

     character(len=18)     :: comm_name    ! print name of communicator
     character(len=ncnode) :: mynode_name  ! name this node

     integer :: mynode_numprocs    ! #processors on this node (for comm)
     integer :: maxnode_numprocs   ! max #processors on any node (for comm)

     integer, dimension(:),   allocatable :: mynode_proclist    ! list of ranks
     integer, dimension(:),   allocatable :: allnode_numprocs   ! #s on each node
     integer, dimension(:,:), allocatable :: allnode_proclists  ! all ranks
                                                                !    allnode_proclists(j,k) = rank of j'th process on k'th node
                                                                !      =-1 if j > allnode_numprocs(k)

     integer, dimension(:), allocatable :: node_map  ! rank->node map
                                                     !    rank r is on node allnode_names(node_map(r))

     character(len=ncnode), dimension(:), allocatable :: allnode_names ! all nodenames

     character(len=255) :: save_path
  end type comdat

  !=====================================================================
  !  **** module variables ****
  !
  !  for all known communicators...

  integer :: ncomm = 0

  type (comdat) :: comm_data(max_comm)

  !  file system information...

  logical      :: pll_file_sys_flag = .TRUE.    ! (default value)
  character*20 :: file_sys_ident    = 'UNKNOWN' ! (default value)

  !  if(pll_file_sys_flag) -- rank(0) of each communicator is "I/O process".
  !  if(.NOT.pll_file_sys_flag) -- 1st rank on each node for each communicator
  !                                serves as "I/O process".
  !
  !  file_sys_ident:
  !  "PPPLCLUSTER" file system -- pll_file_sys_flag = .FALSE.
  !-------------------------------------
  !  (for module local use only)...

  integer :: indx = 0   ! index into comm_data (private use)

#ifdef __MPI
  integer :: istatus(MPI_STATUS_SIZE)  ! (local for MPI_RECV calls)
#endif

  character*512 :: cwd    ! current working directory (this rank) ( GLOBAL !! )
  character*512 :: cwd_0  ! current working directory (rank 0)

  character(len=*), parameter :: cMd = 'portlib/mpi_proc_data'
  !=====================================================================

CONTAINS

  !---------------------------------------------------
  !  PUBLIC...

! 02/07/2012 CLF: How about sanity check for group?
!  subroutine mpi_env_comm_add(comm, group, ierr)
  subroutine mpi_env_comm_add(comm, ierr, name, group_in)

    ! add a communicator to the list
    ! (if already there, do nothing)
    ! 02/15/2012 CLF: add by group-name if process is not in group
    !                 in this case comm=MPI_COMM_NULL
    ! error returns:
    !   ierr=1 means the communicator data setup failed
    !   ierr=2 means the maximum number of communicators has been exceeded.
    !   CAUTION: for collective operation errors MPI_ABORT may be issued...

    ! MPI version: 
    ! it is assumed that the caller has already done MPI_INIT(...).

    ! serial version:
    ! single dummy communicator data setup; hostname still fetched.

    integer, intent(in)  :: comm
    integer, intent(out) :: ierr
    character(len=*), intent(in),optional   :: name       ! group name 
    integer, intent(in),optional            :: group_in   ! group handle

    !----------------------------------
    integer :: icomm,ii,jj,kk,imatch,ierloc,iam,ilen

    integer :: ilen_max         ! length of node name array
    integer :: isize_all        ! = num_procs, total # mpi processes ( in comm world )
    integer :: ierr_my, iabort         
    integer :: ifrom, ito, itag, ibase
    integer :: ichname(140)
    integer :: nact, imax

    character(len=140), dimension(:), allocatable :: allnames    ! node names, subscripted by process mpi ID
    character(len=140), dimension(:), allocatable :: actnames    ! node names - unique list ( 1 entry each node )

    integer, dimension(:),   allocatable :: iactcount
    integer, dimension(:,:), allocatable :: iactnames   ! character representation of actnames strings
    !
    character(len=*), parameter :: cSr = cMd//'.mpi_env_comm_add'
    !-------------------------------------------------------------------------------------------

    ierr = 0
    icomm = comm_arg(comm)                   ! mpi_proc_data.comm_arg

    if (present(name)) then
       call find_comm(name,imatch)
    else   
       call find_comm(icomm,imatch)
    endif

    if(imatch.gt.0) return
    
    ! no match, add to list
    
    if(ncomm.eq.max_comm)    then
       ierr=2
       if(loc_debug) write(0,*) ' ?mpi_proc_data: max_comm exceeded. '
       return
    endif

    ncomm = ncomm + 1
    
    comm_data(ncomm)%comm = icomm
    comm_data(ncomm)%mynode_name = ' '
    
    indx = ncomm

    
#ifdef __MPI

   !Check for group name
    if (present(name)) then
       comm_data(ncomm)%comm_name=name
    else if (icomm .eq. MPI_COMM_WORLD) then
       comm_data(ncomm)%comm_name='MPI_COMM_WORLD'
    else
       comm_data(ncomm)%comm_name=' '
    endif
    if (icomm .eq. MPI_COMM_NULL) then
       comm_data(ncomm)%myproc=-1
       return
    endif
      
    l_get_proc_names: do          !  one trip (only) loop; uses   exit <loop> for error jump
       iabort = 0

       call MPI_COMM_GROUP(icomm, comm_data(ncomm)%group, ierloc)
       call ckerr('MPI_COMM_GROUP',ierloc,ierr)

       call MPI_COMM_RANK(icomm, comm_data(ncomm)%myproc, ierloc)
       call ckerr('MPI_COMM_RANK',ierloc,ierr)

       call MPI_COMM_SIZE(icomm, comm_data(ncomm)%num_procs, ierloc)
       call ckerr('MPI_COMM_SIZE',ierloc,ierr)

       if(ierr.gt.0) exit l_get_proc_names

       isize_all = comm_data(ncomm)%num_procs
       iam = comm_data(ncomm)%myproc

       !  this processor (node name)...
       comm_data(ncomm)%mynode_name = ' '
       call MPI_GET_PROCESSOR_NAME(comm_data(ncomm)%mynode_name, ilen, ierloc)
       call ckerr('MPI_GET_PROCESSOR_NAME',ierloc,ierr)

       if(ierr.gt.0) exit l_get_proc_names

       if(loc_debug.AND.(iam.eq.0)) then
          write(0,*) ' %mpi_proc_data(debug): icomm,isize_all = ', &
               icomm,isize_all
       endif

       !  length of name -- 
       ilen = len(trim(comm_data(ncomm)%mynode_name))

       !  length of name string --                                 [ Must be ncnode ?? jec
!      ilen_max = len(comm_data(ncomm)%mynode_name)                [ yes CLF
       ilen_max = ncnode

       ibase = 0

       iabort = 1  ! subsequent collective failures => MPI_ABORT call...

       !--------------------------
       !  gather all node names on proc 0
       !--------------------------

       i_get_allnames:  if(iam.eq.0) then        !  V  iam 0
          if(allocated(allnames)) deallocate(allnames,actnames,iactcount)
          allocate(allnames(0:isize_all-1),actnames(isize_all))
          allocate(iactcount(isize_all))

          allnames(iam) = comm_data(ncomm)%mynode_name

          l_rcv_names: do ifrom=1,isize_all-1

             itag = ibase + ifrom
             call MPI_RECV(ichname, ilen_max, MPI_INTEGER, ifrom, &
                itag, icomm, istatus, ierloc)
             call ckerr('MPI_RECV',ierloc,ierr)
             if(ierr.gt.0)   exit

             allnames(ifrom)=' '
             do ii=1,ilen_max
                if(ichname(ii).eq.0)   exit
                allnames(ifrom)(ii:ii) = char(ichname(ii))
             enddo
          enddo l_rcv_names

       else                                      !  V   not 0 ( i_get_allnames 

          ichname = 0
          do ii=1,ilen
             ichname(ii) = ichar(comm_data(ncomm)%mynode_name(ii:ii))
          enddo

          itag = ibase + iam
          call MPI_SEND(ichname, ilen_max, MPI_INTEGER, 0, &
               itag, icomm, ierloc)
          call ckerr('MPI_SEND',ierloc,ierr)
          if(ierr.gt.0) exit

       endif  i_get_allnames

       if(ierr.gt.0) exit l_get_proc_names

       !--------------------------
       !  OK, build node table (on iam==0) and broadcast
       !--------------------------

       i_nodetbl: if(iam.eq.0) then        !  V  iam 0
          nact = 1
          imatch = 1
          actnames(nact) = allnames(iam)
          iactcount(imatch) = 1

          l_match_all: do ii=1,isize_all-1
             if(allnames(ii).eq.allnames(ii-1)) then
                ! same node as previous process - skip search
                iactcount(imatch)=iactcount(imatch) + 1
                cycle l_match_all
             endif

             imatch=0
             l_match_1: do jj=1,nact-1
                if(allnames(ii).eq.actnames(jj)) then
                   imatch=jj
                   iactcount(imatch)=iactcount(imatch) + 1     ! increment proc count, this node
                   exit l_match_1
                endif
             enddo l_match_1

             if(imatch.gt.0) cycle l_match_all

             ! we have a new node ( host ) name.. add to actnames
             nact = nact + 1
             actnames(nact)=allnames(ii)
             imatch = nact
             iactcount(imatch) = 1                             ! initialise proc count, this node

          enddo l_match_all

          ! maximum process count over all nodes:

          imax = maxval(iactcount(1:nact))
          
          ! build integer representation
          
          if(allocated(iactnames)) deallocate(iactnames)
          allocate(iactnames(ilen_max,nact))
          do ii=1,nact
             do jj=1,ilen_max
                if(actnames(ii)(jj:jj).eq.' ') then
                   iactnames(jj:ilen_max,ii)=0
                   exit
                else
                   iactnames(jj,ii) = ichar(actnames(ii)(jj:jj))
                endif
             enddo
          enddo
          if(loc_debug ) &
               write(0,*) '%mpi_proc_data(debug): ' &
               ,'MPI_BCAST #nodes and maxval(proc/node):',nact,imax
       endif i_nodetbl                                                          ! ^ iam 0


       !  broadcast #nodes...                                                   ! V   ALL

       call MPI_BCAST(nact, 1, MPI_INTEGER, 0, icomm, ierloc)
       call ckerr('MPI_BCAST(nact)',ierloc,ierr)
       if(ierr.gt.0)   exit l_get_proc_names

       !  broadcast max #process on a single node
       
       call MPI_BCAST(imax, 1, MPI_INTEGER, 0, icomm, ierloc)
       call ckerr('MPI_BCAST(imax)',ierloc,ierr)
       if(ierr.gt.0)   exit l_get_proc_names
       
       !  all ranks record these

       comm_data(ncomm)%num_nodes = nact
       comm_data(ncomm)%maxnode_numprocs = imax

       !--------------------------
       !  all ranks can now allocate...
       !--------------------------
       
       allocate( comm_data(ncomm)%allnode_numprocs(0:nact-1) )
       allocate( comm_data(ncomm)%allnode_names(0:nact-1) )
       allocate( comm_data(ncomm)%mynode_proclist(imax) )
       allocate( comm_data(ncomm)%allnode_proclists(imax,0:nact-1) )
       comm_data(ncomm)%save_path = ' '
       allocate( comm_data(ncomm)%node_map(0:isize_all-1) )

       i_alloc: if(iam.gt.0) then        ! V  not 0
          if(allocated(iactnames)) deallocate(iactnames)
          allocate(iactnames(ilen_max,nact))
          if(allocated(actnames)) deallocate(actnames)
          allocate(actnames(nact))
          if(allocated(iactcount)) deallocate(iactcount)
          allocate(iactcount(nact))
       endif i_alloc                                                            !  ^ not 0
       
       !  broadcast all node proclist lengths                                   !  V ALL

       if(loc_debug.AND.(iam.eq.0)) then
          write(0,*) '%mpi_proc_data(debug): broadcast proclist lengths...'
       endif

       call MPI_BCAST(iactcount, nact, MPI_INTEGER, 0, icomm, ierloc)
       call ckerr('MPI_BCAST(iactcount)',ierloc,ierr)
       if(ierr.gt.0)   exit l_get_proc_names

       !  broadcast all node names as integer

       if(loc_debug.AND.(iam.eq.0)) then
          write(0,*) '%mpi_proc_data(debug): broadcast node names as int...'
       endif

       call MPI_BCAST(iactnames, nact*ilen_max, MPI_INTEGER, 0, icomm, ierloc)
       call ckerr('MPI_BCAST(iactnames)',ierloc,ierr)
       if(ierr.gt.0) exit l_get_proc_names

       !  iam.gt.0: convert int to CHARACTER strings                            !  ^ ALL

       if(iam.gt.0) then        !  V not 0
          do ii=1,nact
             actnames(ii) = ' '
             do jj=1,ilen_max
                if(iactnames(jj,ii).eq.0) exit
                actnames(ii)(jj:jj)=char(iactnames(jj,ii))
             enddo
          enddo
        endif                                                                    !  ^ not 0

       !---------------------------
       !  all ranks
       !  can now set list of node names and identify mynode
       !  can also set list sizes
       !---------------------------                                             !  V ALL

       do ii=1,nact
          comm_data(ncomm)%allnode_names(ii-1) = actnames(ii)
          comm_data(ncomm)%allnode_numprocs(ii-1) = iactcount(ii)
          if(actnames(ii).eq.comm_data(ncomm)%mynode_name) then
             comm_data(ncomm)%mynode = ii-1
             comm_data(ncomm)%mynode_numprocs = iactcount(ii)
          endif
       enddo
       
       !---------------------------
       !  rank 0 now fills in the process lists for each node; which will
       !  be broadcast
       !---------------------------
       
       !  iactcount gets recomputed...
       
       i_recount: if(iam.eq.0) then        ! V iam 0
          comm_data(ncomm)%allnode_proclists = -1

          imatch=1
          iactcount(1:nact)=0
          do ii=0,isize_all-1

             if(allnames(ii).eq.actnames(imatch)) then
                continue
             else
                do kk=1,nact
                   if(allnames(ii).eq.actnames(kk)) then
                      imatch=kk
                      exit
                   endif
                enddo
             endif
             
             iactcount(imatch)=iactcount(imatch)+1
             jj=iactcount(imatch)
             comm_data(ncomm)%allnode_proclists(jj,imatch-1) = ii
          enddo
       endif i_recount                                                          ! ^ iam 0

       ! broadcast process lists                                                ! V ALL
       
       if(loc_debug.AND.(iam.eq.0)) then
          write(0,*) '%mpi_proc_data(debug): broadcast node proclists...'
       endif
       
       call MPI_BCAST( comm_data(ncomm)%allnode_proclists, nact*imax, &
            MPI_INTEGER, 0, icomm, ierloc)

       call ckerr('MPI_BCAST(allnode_proclists)',ierloc,ierr)
       if(ierr.gt.0) exit l_get_proc_names

       ii = comm_data(ncomm)%mynode
       
       comm_data(ncomm)%mynode_proclist = comm_data(ncomm)%allnode_proclists(:,ii)

       exit
    enddo  l_get_proc_names  !====================================================================

    if((ierr.ne.0).and.(iabort.gt.0))    then
       write(0,*) ' ?mpi_proc_data: ierr=',ierr,' for communicator: ',icomm
       write(0,*) '  process #',comm_data(ncomm)%myproc, &
            ' issues MPI_ABORT(...)'
       call MPI_ABORT(icomm,ierr,ierloc)
    endif

    !  MPI_ALLREDUCE'd error code...
    
    call collective_errchk('mpi_env_comm_add',icomm,ierr)
    if(ierr.gt.0) then
       ierr=1
       if(loc_debug.and.(iam.eq.0)) &
            write(0,*) ' ?mpi_proc_data: setup failed for comm = ',comm
    else
       if(loc_debug)   &
            write(0,*) ' %mpi_proc_data: MPI comm_data(',ncomm,')'&
            ,' setup OK on ', iam
    endif
    
    if(ierr.eq.0)    then
       !  build node map
       
       do ii=0,nact-1
          do jj=1,comm_data(ncomm)%allnode_numprocs(ii)
             kk = comm_data(ncomm)%allnode_proclists(jj,ii)
             comm_data(ncomm)%node_map(kk) = ii
          enddo
       enddo
    endif
    
#else
    comm_data(ncomm)%num_procs = 1
    comm_data(ncomm)%num_nodes = 1
    comm_data(ncomm)%myproc = 0
    comm_data(ncomm)%mynode = 0

    call hostnm( comm_data(ncomm)%mynode_name )

    comm_data(ncomm)%mynode_numprocs = 1
    comm_data(ncomm)%maxnode_numprocs = 1

    allocate( comm_data(ncomm)%mynode_proclist(1) )
    allocate( comm_data(ncomm)%allnode_numprocs(0:0) )
    allocate( comm_data(ncomm)%allnode_proclists(1,0:0) )
    allocate( comm_data(ncomm)%allnode_names(0:0) )
    allocate( comm_data(ncomm)%node_map(0:0) )

    comm_data(ncomm)%mynode_proclist(1) = 0
    comm_data(ncomm)%allnode_numprocs(0) = 1
    comm_data(ncomm)%allnode_proclists(1,0) = 0
    comm_data(ncomm)%allnode_names(0) = comm_data(ncomm)%mynode_name
    comm_data(ncomm)%node_map(0) = 0

    if(loc_debug) write(0,*) ' %mpi_proc_data: serial comm_data(1) setup OK.'
#endif

  end subroutine mpi_env_comm_add

  !---------------------------------------------------
  !  PUBLIC...

  subroutine mpi_print_info(comm,proc_id,ierr)

    ! print available information for communicator (comm)
    ! on process (proc_id); issue MPI_BARRIER call so that 
    ! output from a single processor only is output without
    ! concurrent outputs from other processors

    integer, intent(in)  :: comm     ! communicator
    integer, intent(in)  :: proc_id  ! process number
    integer, intent(out) :: ierr     ! completion status, 0=normal

    !  error returns:
    !    ierr=1 -- comm unknown
    !    ierr=2 -- proc_id unknown for (comm)

    !-------------------------------------------------
    integer :: icomm,imax,iam,ierloc,inum,ii,jj,kk
    logical :: okmap
    !-------------------------------------------------

    ierr = 0
    icomm = comm_arg(comm)
    call find_comm(icomm,indx)
    
    iam = -1
    
    if(indx.eq.0) then
       ierr=1

    else
       iam = comm_data(indx)%myproc
       imax = comm_data(indx)%num_procs - 1
       if((proc_id.lt.0).OR.(proc_id.gt.imax)) then
          ierr=2
       endif
    endif

    if((ierr.eq.0).AND.(iam.eq.proc_id)) then

       write(0,*) ' '
       write(0,*) '-------------> node/proc info for proc_id = ',proc_id
       write(0,*) '  communicator:  ',icomm
       write(0,*) '  group:         ',comm_data(indx)%comm_name
       write(0,*) '  internal indx: ',indx
       write(0,*) '    myproc:      ',comm_data(indx)%myproc
       write(0,*) '    mynode:      ',comm_data(indx)%mynode
       write(0,*) '    num_procs:   ',comm_data(indx)%num_procs
       write(0,*) '    num_nodes:   ',comm_data(indx)%num_nodes
       write(0,*) '    node name:   ',comm_data(indx)%mynode_name
       inum = comm_data(indx)%mynode_numprocs
       write(0,*) '    #procs on node:  ',inum
       write(0,*) '    max #procs/node: ',comm_data(indx)%maxnode_numprocs
       write(0,*) '    mynode proc#s:   ', &
            comm_data(indx)%mynode_proclist(1:inum)
       okmap = .TRUE.
       do ii=0,comm_data(indx)%num_nodes-1
          inum = comm_data(indx)%allnode_numprocs(ii)
          write(0,*) '    '//trim(comm_data(indx)%allnode_names(ii))// &
               ':  ',comm_data(indx)%allnode_proclists(1:inum,ii)
          do jj=1,inum
             kk=comm_data(indx)%allnode_proclists(jj,ii)
             if(comm_data(indx)%node_map(kk).ne.(ii)) okmap=.FALSE.
          enddo
       enddo
       write(0,*) '    node_map check OK = ',okmap
       write(0,*) ' '

   endif

#ifdef __MPI
    if(loc_debug) write(0,*) 'mpi_print_info: I am ',iam,' call MPI_BARRIER'
    call MPI_BARRIER(icomm,ierloc)
#endif

  end subroutine mpi_print_info

  !---------------------------------------------------
  !  PUBLIC...

  subroutine mpi_proc_info_comm(comm,ierr, &
       myproc,num_procs,mynode,num_nodes, &
       maxnode_numprocs, &
       node_number, &
       node_name,node_numprocs,node_proclist, &
       parallel_file_system, &
       file_system_ident )

    !  public interface to retrieve module data for a single communicator: comm

    integer, intent(in)  :: comm   ! communicator
    integer, intent(out) :: ierr   ! completion status (0=OK)

    !  error returns:
    !    ierr=1 -- comm unknown
    !    ierr=2 -- an proclist array size is too small
    !    ierr=3 -- node_number value out of range

    integer, intent(out), optional :: myproc     ! process # (rank)
                                                 !   range of values (0:num_procs-1)
    integer, intent(out), optional :: num_procs  ! #processors in (comm)

    integer, intent(out), optional :: mynode     ! node # 
                                                 !   range of values (0:num_nodes-1)
    integer, intent(out), optional :: num_nodes  ! #nodes in (comm)

    integer, intent(out), optional :: maxnode_numprocs  ! max #procs/node
                                                        !  this is the maximum #processors on any one node in (comm)

    !-----------
    !  set this to select node of interest.  If omitted, mynode is assumed!

    integer, intent(in), optional :: node_number

    !  NOTE: valid range is 0:num_nodes-1
    !-----------

    character*(*), intent(out), optional :: node_name

    integer, intent(out), optional :: node_numprocs

    integer, intent(out), dimension(:), optional :: node_proclist
    !  where dimension exceeds list size: extra array element values set to -1

    !-------------
    ! optionally return internal setting of parallel_file_system flag
    !-------------

    logical, intent(out), optional :: parallel_file_system

    character*(*), intent(out), optional :: file_system_ident

    !---------------------------------------
    integer :: isize1,isize1a,icomm,inode
    !---------------------------------------

    ierr=0
    icomm = comm_arg(comm)
    call find_comm(icomm,indx)

    if(indx.eq.0)    then
       ierr=1
       return
   endif

    if(present(myproc)) myproc = comm_data(indx)%myproc
    if(present(num_procs)) num_procs = comm_data(indx)%num_procs
    if(present(mynode)) mynode = comm_data(indx)%mynode
    if(present(num_nodes)) num_nodes = comm_data(indx)%num_nodes

    if(present(maxnode_numprocs)) maxnode_numprocs = comm_data(indx)%maxnode_numprocs

    if(present(parallel_file_system)) parallel_file_system = pll_file_sys_flag

    if(present(file_system_ident)) file_system_ident = file_sys_ident

    inode = comm_data(indx)%mynode
    if(present(node_number)) then
       inode = node_number
       if((inode.lt.0).OR.(inode.ge.comm_data(indx)%num_nodes)) then
          ierr=3
          if(present(node_name)) node_name = ' '
          if(present(node_numprocs)) node_numprocs = 0
          if(present(node_proclist)) node_proclist = -1
          return
       endif
    endif

    if(present(node_name)) node_name = comm_data(indx)%allnode_names(inode)

    if(present(node_numprocs)) node_numprocs = comm_data(indx)%allnode_numprocs(inode)

    if(present(node_proclist)) then
       node_proclist = -1
       isize1 = comm_data(indx)%allnode_numprocs(inode)
       isize1a = size(node_proclist)
       if(isize1a.lt.isize1) then
          ierr=2
       else
          node_proclist(1:isize1) = comm_data(indx)%allnode_proclists(1:isize1,inode)
       endif
   endif

  end subroutine mpi_proc_info_comm

  !------------------------------------------------
  !  (public routine)

  subroutine mpi_proc_info_name(name,status,icomm,irank,num_procs)

    !  find communicator and rank for comm-name
#ifdef __MPI
    use mpi
#endif

    character(len=*), intent(in)   :: name       ! group name 
    logical, intent(out)           :: status     ! status of group = -1 if not created
    integer, intent(out)           :: icomm      ! group communicator
    integer, intent(out), optional :: irank      ! rank of process in group
    integer, intent(out), optional :: num_procs  ! # of processes in group

    !------------------
    integer :: ii
    !------------------

    status= .false.
    icomm = -1
    if(present(irank)) irank = -1
    if(present(num_procs)) num_procs = -1

#ifdef __MPI
    
    ii=0
    call find_comm(trim(name),ii)
    if (ii .gt. 0) then
       status = .true.
       icomm = comm_data(ii)%comm       
       if (comm_data(ii)%comm .ne. MPI_COMM_NULL) then
          if(present(irank)) irank = comm_data(ii)%myproc
          if(present(num_procs)) num_procs = comm_data(ii)%num_procs
       endif
    endif
#endif    
    
  end subroutine mpi_proc_info_name



  !---------------------------------------------------
  !  PUBLIC...

  subroutine mpi_set_filesys(filesys_id)

    ! set file system identification

    character*(*), intent(in) :: filesys_id

    !---------------------------------
    logical :: assume_pll
    !---------------------------------

    file_sys_ident = filesys_id

    assume_pll = .TRUE.
    if(filesys_id.eq.'PPPLCLUSTER') assume_pll = .FALSE.

    call mpi_parallel_file_system(assume_pll)

  end subroutine mpi_set_filesys

  !---------------------------------------------------
  !  private:
  subroutine mpi_parallel_file_system(iswitch)

    ! set/clear parallel file system flag

    logical, intent(in) :: iswitch

    pll_file_sys_flag = iswitch

  end subroutine mpi_parallel_file_system

  !---------------------------------------------------
  !  PUBLIC...

  logical function check_mpi_io_node(comm,ierr)

    ! return .TRUE. if for (comm) this process is an I/O node
    !   if(pll_file_sys_flag) this is rank 0 only
    !   if(.NOT.pll_file_sys_flag) it is 1st process on each node list

    integer, intent(in) :: comm    ! communicator
    integer, intent(out) :: ierr   ! status, 0=normal; 1 if (comm) unknown

    !----------------------------
    integer :: icomm,iam
    logical :: iaction
    !----------------------------

    check_mpi_io_node = .FALSE.

    ierr = 0
    icomm = comm_arg(comm)
    call find_comm(icomm,indx)
    if(indx.eq.0) then
       ierr=1
       return
    endif

    iam = comm_data(indx)%myproc
    iaction = (iam.eq.0)  ! root process always acts...

    if(.NOT.pll_file_sys_flag) then
       !  only act on 1st process on each node...
       iaction = iam.eq.comm_data(indx)%mynode_proclist(1)
    endif

    check_mpi_io_node = iaction

  end function check_mpi_io_node

  !---------------------------------------------------
  !  PUBLIC...

  subroutine nodepath_translator(path_in,path_out,comm,irank)

    ! for given rank on given communicator, give a file path translation
    !
    ! for a system with pll_file_sys_flag=T, path_in = path_out ALWAYS

    character*(*), intent(in) :: path_in   ! input path
    character*(*), intent(out) :: path_out ! translated path

    integer, intent(in) :: comm    ! communicator
    integer, intent(in) :: irank   ! rank of interest

    !--------------------------------------
    integer :: icomm,indx,inum,inode_indx,idm1
    character*140 :: node_name
    !--------------------------------------
    path_out = path_in

    if(pll_file_sys_flag) return

    icomm = comm_arg(comm)
    call find_comm(icomm,indx)
    if(indx.eq.0) return

    inum = comm_data(indx)%num_procs
    if((irank.lt.0).or.(irank.ge.inum)) then
       write(0,*) ' %nodepath_translator: rank=',irank, ' out of range 0 to ',inum-1
       return
    endif

    inode_indx = comm_data(indx)%node_map(irank)
    node_name = comm_data(indx)%allnode_names(inode_indx)

    if(file_sys_ident.eq.'PPPLCLUSTER') then
       if(path_in(1:7).eq.'/local/') then
          idm1 = index(node_name,'.') - 1
          if(idm1.le.0) idm1 = len(trim(node_name))
          
          path_out = '/l/'//node_name(1:idm1)//'/'//path_in(8:)
          
       endif
    endif
   
  end subroutine nodepath_translator

  !---------------------------------------------------
  !  PUBLIC...

  subroutine mpi_mkdir(comm,path,ierr,clean)

    ! create a file path for processes in (comm)
    !   if(pll_file_sys_flag) only rank 0 does this
    !   if(.NOT.pll_file_sys_flag) the 1st process on each node does it

    integer, intent(in) :: comm    ! communicator
    character*(*), intent(in) :: path

    integer, intent(out) :: ierr   ! completion status (0=normal)

    !---------------------------------------------------------------------
    !  CAREFUL: could be dangerous...

    logical, intent(in), OPTIONAL :: clean  ! remove existing files in path

    !  (this removes files in the path directory, if it already exists;
    !  hidden files and subdirectories are not touched)
    !---------------------------------------------------------------------

    !------------------------------
    integer :: icomm,iam,ierloc
    logical :: iclean,iaction
    !------------------------------
    character(len=128) :: zcmd

    ierr=0
    iclean = .FALSE.
    if(present(clean)) iclean = clean
    
    iaction = check_mpi_io_node(comm,ierr)
    
    icomm = comm_data(indx)%comm

    !  if processed this path on preceding call for (comm), no action...
    if(.NOT.iclean) then
       if(path.eq.comm_data(indx)%save_path) then
          iaction=.FALSE.
       else
          comm_data(indx)%save_path = path
       endif
    endif

    if(iaction) then
       zcmd = 'mkdir -p '//trim(path)
       call execute_command_line(zcmd,exitstat=ierr)
       if(ierr.eq.0) then
          if(iclean) then
             call fclean_dir(path,ierr)
          endif
       endif
    endif

    ! make sure directory creation is complete on all nodes before continuing
    ! use collective
    
    call collective_errchk('mpi_mkdir',icomm,ierr)
    
  end subroutine mpi_mkdir
  !---------------------------------------------------
  !  PUBLIC...

  subroutine mpi_file_bcast_single(comm,filenam,ierr,newname)
    
    ! broadcast single file from [node 0] to remaining nodes
    ! assume only [node 0] knows the names.
    ! no action if(pll_file_sys_flag)
    
    integer,       intent(in)  :: comm    ! communicator
    character*(*), intent(in)  :: filenam
    integer,       intent(out) :: ierr    ! return status, 0=normal
                                          !  1 means unknown communicator; 2 means a file transfer failed although
                                          !  the code may not survive such a failure...

    character*(*), intent(in), optional :: newname  ! set this if copied
                                                    ! file should have a different name than the original file; if omitted
                                                    ! the original name is kept.

    !----------------------------------------------
    character*200 :: flist(1),nname(1)
    !----------------------------------------------

    flist(1) = filenam

    if(present(newname))    then
       nname(1)=newname
       call mpi_file_bcast_list(comm,flist,ierr,newname=nname)

    else
       call mpi_file_bcast_list(comm,flist,ierr)
    endif

  end subroutine mpi_file_bcast_single
    
  subroutine mpi_file_bcast_list( comm, filelist, ierr, newname)

    ! Broadcast files from [node 0] to remaining nodes, optionally renaming output files
    ! [ portlib/portlib_parcopy performs the actual send/receive data ]
    !
    ! Assume only [node 0] knows the names.  No action if(pll_file_sys_flag)

    integer, intent(in) :: comm  ! communicator

    character(len=*), intent(in), dimension(:) :: filelist  ! (input) file names

    integer, intent(out) :: ierr       ! return status, 0=normal
                                       ! 1 means unknown communicator; 
                                       ! 2 means a file transfer failed although
                                       !   the code may not survive such a failure...
                                                                 
    character(len=*), intent(in), optional, dimension(:) :: newname  ! new/output file names df to same as input
                                                                     ! set this optional list if each copied file should 
                                                                     ! have a different name than the original file; 
                                                                     ! if omitted original name is used for all copies.
                                                                     ! MAY HAVE LESS ENTRIES than filelist ; 
                                                                     ! >> renamed files must preceed non-renamed files in filelist <<
                                                                 

    integer :: icomm,innodes,ii,jj 

    integer :: innew          ! number of entries in newname array
    integer :: ierloc,iam       
    integer :: inum_files     ! Number of file transfers requested
    integer :: ilenc            
    integer :: inum_actual    ! number which must be sent
    integer :: ina            !  "  "   + 1
    integer :: iact           !  0 - none ; 1 - local copy ; 2 - MPI transfer

    logical :: iaction      
      
    integer, dimension(:), allocatable :: ifrom,ito

    integer, parameter :: ncfnm = 256  ! # chars/file name

    character(len=ncfnm),dimension(:),  allocatable :: from_list  ! (local) copy of filelist
    character(len=ncfnm),dimension(:),  allocatable :: to_list    ! filelist OR newname
    character(len=ncfnm),dimension(:,:),allocatable :: ffrom      ! copy of from_list, for bcast/recieve
    character(len=ncfnm),dimension(:,:),allocatable :: fto        ! copy of   to_list, for bcast/recieve
    !-------------------------------------

    ierr=0

    icomm = comm_arg(comm)
    call find_comm(icomm,indx)

    if(indx.eq.0) then
       ierr=1
       return
    endif

 
    call mpi_cwd_check(indx,ierr)
    
    if(ierr.ne.0) then
       write(0,*) '??mpi_file_bcast_list: ierr =',ierr,indx
       write(0,*) '           return'
       return
    endif
    
    iam=comm_data(indx)%myproc
    if(pll_file_sys_flag) then
       innodes = 1
    else
       innodes=comm_data(indx)%num_nodes
    endif

    allocate(ifrom(innodes))
    allocate(ito(innodes))

    if(present(newname)) then
       innew=size(newname)
    else
       innew=0
    endif

    i_filenm: if(iam.eq.0) then    !  V  iam 0
       inum_files = size(filelist)

       inum_actual = 0
       allocate(from_list(inum_files),to_list(inum_files))
       
       l_filenm: do ii=1,inum_files
          if(filelist(ii).eq.' ')  then
             write(0,*) ' ?? mpi_file_bcast_list: blank filename received.'
             write(0,*) '    "filelist" list index: ',ii
             ierr = 2
             exit l_filenm
          endif
          if(ii.le.innew) then
             if(newname(ii).eq.' ') then
                write(0,*) ' ?? mpi_file_bcast_list: blank new name received.'
                write(0,*) '    "newname" list index: ',ii
                ierr = 2
                exit l_filenm
             endif
          endif
          
          ina = inum_actual+1
          from_list(ina) = filelist(ii)
          if(ii.le.innew) then
             to_list(ina) = newname(ii)
          else
             to_list(ina) = filelist(ii)
          endif

          !  what to do with this pair ??
          call mpi_copy_status(from_list(ina),to_list(ina),iact)

          if(iact.eq.0) then
             continue                                          !    no action
          else if(iact.eq.1) then
             ! rank 0 only
             call fcopy(from_list(ina),to_list(ina),ierr)      !    local copy only
             if(ierr.ne.0)   exit l_filenm
          else if(iact.eq.2)    then
             !  MPI distribution needed
             inum_actual = inum_actual + 1                     !     MPI copy needed - increment counter
          endif

       enddo l_filenm
    else    !  ^ iam 0
       !    set counters                                                        !  V not 0
       inum_files = 0
       inum_actual = 0
    endif i_filenm 

    call collective_errchk('mpi_file_bcast_list',icomm,ierr)
    if(ierr.ne.0) return

    l_actual: do ii=1,innodes
       !         broadcast actual # of transfers ( those with  iact = 2 )
       jj=ii-1
       ifrom(ii) = 0                                                            !   all transfers are from  root
       if(pll_file_sys_flag) then
          ito(ii) = 0
       else
          ito(ii) = comm_data(indx)%allnode_proclists(1,jj)
       endif
       
       if(ito(ii).eq.ifrom(ii)) cycle l_actual

#ifdef __MPI
       if(iam.eq.ifrom(ii)) then                      !  V  iam 0
          call MPI_SEND(inum_actual,1,MPI_INTEGER,ito(ii), &
               2001,icomm,ierloc)
          call ckerr('MPI_SEND(inum_actual)',ierloc,ierr)
       else if(iam.eq.ito(ii))    then                ! V not 0
          call MPI_RECV(inum_actual,1,MPI_INTEGER,ifrom(ii), &
               2001,icomm, istatus, ierloc)
          call ckerr('MPI_RECV(inum_actual)',ierloc,ierr)
       endif
#endif
    enddo l_actual

    call collective_errchk('mpi_file_bcast_list',icomm,ierr)                    !  V ALL
    if(ierr.ne.0) return

    i_actual_nm: if(inum_actual.gt.0) then

       ! exchange file names
       ! only procs actively involved in MPI file copying enter this section

       allocate(ffrom(inum_actual,innodes),fto(inum_actual,innodes))
       ilenc = len(ffrom(1,1))

       l_copy_node_nm: do ii=1,innodes
          l_copy_file_nm: do jj=1,inum_actual
             
             ! send/receive file names
             
             if(iam.eq.ifrom(ii)) then
                ffrom(jj,ii) = from_list(jj)
                fto(jj,ii) = to_list(jj)
             endif
             
             if(ito(ii).eq.ifrom(ii)) cycle l_copy_file_nm
             
#ifdef __MPI
             if(iam.eq.ifrom(ii)) then
                call MPI_SEND(ffrom(jj,ii),ilenc,MPI_CHARACTER,ito(ii), &
                     2002,icomm,ierloc)
                call ckerr('MPI_SEND(ffrom)',ierloc,ierr)
                call MPI_SEND(fto(jj,ii),ilenc,MPI_CHARACTER,ito(ii), &
                     2003,icomm,ierloc)
                call ckerr('MPI_SEND(fto)',ierloc,ierr)
             else if(iam.eq.ito(ii)) then
                call MPI_RECV(ffrom(jj,ii),ilenc,MPI_CHARACTER,ifrom(ii), &
                     2002,icomm, istatus, ierloc)
                call ckerr('MPI_RECV(ffrom)',ierloc,ierr)
                call MPI_RECV(fto(jj,ii),ilenc,MPI_CHARACTER,ifrom(ii), &
                     2003,icomm, istatus, ierloc)
                call ckerr('MPI_RECV(fto)',ierloc,ierr)
                !  else
                !  other nodes - do nowt
             endif
#endif
          enddo l_copy_file_nm
       enddo l_copy_node_nm
    endif i_actual_nm

    call collective_errchk('mpi_file_bcast_list',icomm,ierr)
    if(ierr.ne.0) return

    i_actual_cp: if(inum_actual.gt.0) then

       ! perform copy -
       ! only procs actively involved in file copying enter this section
       
       l_node_cp: do ii=1,innodes
          l_file_cp: do jj=1,inum_actual

             if(ifrom(ii).eq.ito(ii)) then
                if(ffrom(jj,ii).eq.fto(jj,ii)) cycle l_file_cp
             endif
             if (loc_debug)    then
                write(0,*) 'mpi_file_bcast_list: iam=',iam, &
                     trim(ffrom(jj,ii)),trim(fto(jj,ii))
             endif
             call parcopy(icomm,iam,ifrom(ii),ito(ii), &                        ! portlib/parcopy 
                  trim(ffrom(jj,ii)),trim(fto(jj,ii)),ierloc)                   !    -> portlib/portlib_parcopy
             if(ierloc.ne.0) ierr=2
             
          enddo l_file_cp
       enddo l_node_cp

    endif i_actual_cp

    ! make sure file transfer is complete on all nodes before continuing
    
    call collective_errchk('mpi_file_bcast',icomm,ierr)
    
  end subroutine mpi_file_bcast_list
  
  !---------------------------------------------------
  !  PUBLIC...
  
  subroutine mpi_file_transfer(ierr)
    
    integer, intent(out) :: ierr
    
    write(0,*) ' ?? mpi_file_transfer: not implemented, do not use! '
    ierr = 1
    
  end subroutine mpi_file_transfer
  
  !------------------------------------------------
  !  (private routine, deal with serial vs. MPI)
  
  integer function comm_arg(comm)
    
    integer, intent(in) :: comm
    
#ifdef __MPI
    comm_arg = comm
#else
    comm_arg = 0   ! serial code: ignore
#endif
  end function comm_arg
  !------------------------------------------------
  !  (private)

  subroutine mpi_cwd_check(indx,ierr)
    
    ! verify that all ranks are in the same $cwd
    ! if not: report error
    
    integer, intent(in) :: indx  ! index to communicator
    integer, intent(out) :: ierr ! completion status (0=normal)
    
    !-----------------------------------
    integer :: iam,icomm,ierloc
    
    !-----------------------------------
    
    iam=comm_data(indx)%myproc
    icomm=comm_data(indx)%comm

    ierr = 0

    call getcwd(cwd)
    if(iam.eq.0) cwd_0 = cwd

#ifdef __MPI
    if(loc_debug) write(0,*) 'mpi_cwd_check(',iam,') call MPI_BCAST',trim(cwd_0)
    call MPI_BCAST(cwd_0, len(cwd_0), MPI_CHARACTER, 0, icomm, ierloc)
    if((iam.ne.0).AND.(ierloc.ne.0)) cwd_0='?MPI_BCAST_error'
    ! 07/22/2011 CLF: synch
    if(loc_debug) write(0,*) ' mpi_cwd_check(',iam,') call MPI_BARRIER: ',trim(cwd_0)
    call MPI_BARRIER(icomm, ierr)
    if(ierr.ne.0)    then
       write(0,*) '?mpi_cwd_check(',iam,') MPI_BARRIER error detected!'
    endif

    if(cwd.ne.cwd_0)    then
       ierr=1 
       write(0,*) ' ?mpi_cwd_check(',iam,') cwd mismatch: ',trim(cwd)
       write(0,*) ' ?                                      ',trim(cwd_0)
    endif

    call collective_errchk('mpi_cwd_check',icomm,ierr)
#endif

    if((iam.eq.0).AND.(ierr.ne.0)) then

       write(0,*) ' ================================ '
       write(0,*) ' mpi_cwd_check rank(0) cwd: '//trim(cwd)
       write(0,*) ' ================================ '

    endif
    
  end subroutine mpi_cwd_check

  subroutine mpi_copy_status(cfrom,cto,iact)

    ! file transfer status assessment

    character*(*), intent(in) :: cfrom  ! "from" path
    character*(*), intent(in) :: cto    ! "to" path
    integer, intent(out)      :: iact   ! assessment

    !  iact=0 means, no transfer required
    !  iact=1 means, only local transfer required
    !  iact=2 means, MPI transfer required
    
    iact=0
    if(cfrom.ne.cto) iact=1
    
    if(pll_file_sys_flag) then

       ! MPI transfer never required
       ! local transfer (root proc) if cfrom.ne.cto

       CONTINUE

    else

       if(file_sys_ident.eq.'PPPLCLUSTER') then
          if(cto(1:7).eq.'/local/') then
             iact=2
          else
             if((cto(1:1).ne.'/').AND.(cto(1:1).ne.'~')) then
                if(cwd(1:7).eq.'/local/') then
                   iact=2
                endif
             endif
          endif
       endif
       
    endif
    
  end subroutine mpi_copy_status

  !------------------------------------------------
  !  (private routine, error handling)

  subroutine ckerr(thiscall, ierloc, ierr)

    character*(*), intent(in) :: thiscall   ! label point of error
    integer, intent(in)       :: ierloc     ! error if ierloc.ne.0
    integer, intent(inout)    :: ierr       ! increment if ierloc.ne.0

    if(ierloc.eq.0) return

    if(loc_debug) write(0,*) ' ?mpi_proc_data: error in '//trim(thiscall)
    ierr = ierr + 1
    
  end subroutine ckerr
  !------------------------------------------------
  !  (private routine)

  subroutine find_index_comm(icomm,imatch)

    !  find index to data structure for communicator (icomm)

    integer, intent(in) :: icomm
    integer, intent(out) :: imatch

    !------------------
    integer :: ii
    !------------------

    imatch = 0
    do ii=1,ncomm
       if(comm_data(ii)%comm.eq.icomm) then
          imatch = ii
          exit
       endif
    enddo
    
  end subroutine find_index_comm

  subroutine find_index_name(name,imatch)

    !  find index to data structure for communicator (icomm)

    character(len=*), intent(in)  :: name       ! group name 
    integer, intent(out)          :: imatch     ! index into comm_data
    
    !------------------
    integer :: ii
    !------------------

    imatch = 0
    do ii=1,ncomm
       if(trim(comm_data(ii)%comm_name).eq.trim(name)) then
          imatch = ii
          exit
       endif
    enddo
    
  end subroutine find_index_name


  !------------------------------------------------
  !  (private routine)
  
  subroutine collective_errchk(subr,icomm,ierr)

    character*(*), intent(in) :: subr  ! calling subroutine
    integer, intent(in)    :: icomm    ! verified communicator
    integer, intent(inout) :: ierr     ! this proc error code on input, collective
    !                                  ! error code, on output.

    !---------------------
    integer :: ierloc,iertmp,iam
    !---------------------

    iam = comm_data(indx)%myproc

    !  prevent MPI_ALLREDUCE starting while other procs still
    !  have MPI_SEND/RECV calls pending

    if(loc_debug) write(0,*) ' I, ',iam,' ',comm_data(indx)%comm_name, &
                  ' am in collective_errchk: ierr = ',ierr,' for '//trim(subr)

#ifdef __MPI
    iertmp = ierr
    CALL MPI_ALLREDUCE(iertmp,ierr,1,MPI_INTEGER,MPI_MAX,icomm,ierloc)
    if (ierloc.gt.0) then
       write(0,*) 'collective_errchk:',iam, &
            ' MPI_ALLREDUCE retuned error =',ierloc
    endif
#endif

  end subroutine collective_errchk
  
  !------------------------------------------------
end module mpi_proc_data
