!
! Defines more verbose and flexible system like calls.  See porttest/execsystem_test.f90 
! for some examples.
!

module execsystem
  use logmod, only: openlog, infolog, errorlog, setloglevel, &
       enterlog, exitlog

  implicit none

  integer, save :: EXEC_METHOD = 3     ! method for making the jsystem and jexec calls.  Probably best to set this at the
                                       ! beginning of the program and leave it. 
                                       !   0 -> make call within the process using heap arguments
                                       !   1 -> use static memory for child exec() arguments -- a custom
                                       !        environment is not available and there will be no verbose messages from
                                       !        the child.  Also, the full path to the executable is required.  For jsystem(),
                                       !        the string '/bin/' will prepended to SHELL_EXEC.
                                       !   2 -> make call to a shell server.  The shell server must have already been
                                       !        started.  A custom environment is not available.  The shell server working
                                       !        directory will be set to the cwd of the calling process.
                                       !   3 -> call Fortran execute_command_line

  integer, save :: EXEC_VERBOSE = 0    ! set >0 for verbose messages to stdout, <0 to stderr (not buffered?)

  character(len=32), save :: SHELL_EXEC   = "sh"   ! default shell, room for some path
  character(len=8),  save :: SHELL_LAUNCH = "-c"   ! argument for launching inline shell command

  interface jsystem
     module procedure jsystem_env, jsystem_one
  end interface

  interface jexec
     module procedure jexec_env, jexec_one
  end interface

  private jsystem_env, jsystem_one, jexec_env, jexec_one, EXEC_RUNID, EXEC_IS_SERVER

  character(len=16), save :: EXEC_RUNID = "                "   ! filled with runid after connection to shell server

  logical, save :: EXEC_IS_SERVER = .false.   ! .true. if this process is connected to a shell server

  integer, parameter, private :: EXEC_MAX_STR = 513  ! maximum length of a string with 0 terminator
contains
  !
  ! --------- jconnect_server -------
  ! Connect to a running shell server.  Will attempt to connect to the fifos
  !   <runid>_command.fifo  -> send shell,exec commands as strings to this fifo
  !   <runid>_result.fifo   -> retrieve integer result from this fifo as a string
  !
  ! Only one process can be connected at a time.  A test message will be sent during the open
  ! to check the connection.
  !
  ! For an MPI process, only the ranks in the communicator which have been identified as I/O
  ! ranks for each of the nodes by mpi_proc_data will try to connect to a shell server.  So there
  ! should be one shell server running on each of the nodes.
  !
  ! The shell server should be started with,
  !    shell_server  fifo <verbose> <output> <runid>_command.fifo <runid>_result.fifo
  !        <verbose> = controls verbosity of server to <output>
  !                      -1 verbose output to stderr
  !                       0 no verbose
  !                       1 verbose output to stdout
  !        <output>  = name of file where stdout, stderr for system commands is sent
  !
  !        <runid>_command.fifo = name of fifo for commands.  This is created and read by the 
  !                               shell server
  !        <runid>_result.fifo  = name of fifo used to return the single integer status from the 
  !                               system command.  This is deleted by this subroutine and created
  !                               by the shell server just before returning the test results.
  !
  subroutine jconnect_server(runid, ier, comm, barrier)
    use iso_c_binding, only: c_char
    use mpi_proc_data, only: check_mpi_io_node

#ifdef __MPI
    use mpi, only: mpi_comm_world
#endif

    character(len=*), intent(in)  :: runid    ! runid used to identify fifos for the server, empty to close
    integer,          intent(out) :: ier      ! nonzero on error in this process (not all reduced)

    integer, optional, intent(in) :: comm     ! communicator used as basis for opening a connection.  Only the
                                              ! unique ranks on each node associated with this communicator
                                              ! will open a connection.  This communicator must have been added to
                                              ! mpi_proc_data with mpi_proc_data.mpi_env_comm_add().  If not present,
                                              ! MPI_COMM_WORLD will be used.
    logical, optional, intent(in) :: barrier  ! for an mpi call, an mpi_barrier will be called before exiting
                                              ! this subroutine if barrier=.true.

    character(len=32) :: com   ! name of command fifo
    character(len=32) :: res   ! name of result fifo
    character(len=64) :: msg   ! message to use when testing the connection

    character(kind=c_char) :: ccom(33), cres(33), cmsg(65)  ! c versions
    
    integer, external :: c_connect_shell_server ! c code

    logical :: try_connect   ! .true. if will try a connection 
    integer :: use_comm      ! communicator to use
    integer :: iber          ! barrier error

    ier = 0

    call enterLog("jconnect_server")

    try_connect = .true.

    EXEC_IS_SERVER = .false.

#ifdef __MPI
    if (present(comm)) then
       use_comm = comm
    else
       use_comm = MPI_COMM_WORLD   
    end if

    try_connect = check_mpi_io_node(use_comm,ier)  ! is this an I/O node?

    if (ier/=0) then
       call errorLog("jconnect_server: unknown communicator =",use_comm)
       ier=1
       try_connect = .false.
    end if
#endif

    if (try_connect) then
       if (len_trim(runid)/=0) then
          if (len_trim(EXEC_RUNID)/=0) then
             call errorLog('jconnect_server: A connection to a shell server has already been made')
             ier=1 ; goto 900
          else
             EXEC_RUNID = adjustl(runid)
             call uupper(EXEC_RUNID)
             
             com = trim(EXEC_RUNID)//"_command.fifo"
             res = trim(EXEC_RUNID)//"_result.fifo"
             msg = "opening server for TRANSP run "//trim(EXEC_RUNID)
             
             call cstring(trim(com),ccom,'2C')
             call cstring(trim(res),cres,'2C')
             call cstring(trim(msg),cmsg,'2C')
             
             call infoLog("jconnect_server: making shell server connection to "//trim(com)//", "//trim(res))

             ier = c_connect_shell_server(0, ccom, cres, cmsg)    ! open the connection

             EXEC_IS_SERVER = ier==0
          end if
       else
          com = " " ; res = " " ; msg = " "
          call cstring(trim(com),ccom,'2C')
          call cstring(trim(res),cres,'2C')
          call cstring(trim(msg),cmsg,'2C')
          
          call infoLog("jconnect_server: closing shell server connection")

          ier = c_connect_shell_server(-1, ccom, cres, cmsg)   ! close the open connection
          
          EXEC_RUNID = " "
       end if
    end if

900 continue

#ifdef __MPI
    if (present(barrier)) then
       if (barrier) then
          call MPI_BARRIER(use_comm, iber)
          if (iber/=0) then
             call errorLog("jconnect_server: barrier error")
             ier=1
          end if
       end if
    end if
#endif

    call exitLog("jconnect_server")
  end subroutine jconnect_server

  !
  ! --------- jsystem ----------
  !
  ! Make a system call with an alternate environment.
  !
  ! For an MPI environment, the communicator (default MPI_COMM_WORLD) must have been added to mpi_proc_data
  ! and only I/O ranks will execute the system call.
  !
  function jsystem_env(str, ne, env, comm, barrier, allcpus, method) result(ier)
    use mpi_proc_data, only:

#ifdef __MPI
    use mpi, only:
#endif

    character(len=*), intent(in) :: str       ! shell command to execute
    integer,          intent(in) :: ne        ! number of environment variables or -1 to use the default
    character(len=*), intent(in) :: env(max(0,ne)) ! list of new NAME=VALUE environment to replace current environment

    integer, optional, intent(in) :: comm     ! communicator used as basis for opening a connection.  Only the
                                              ! unique ranks on each node associated with this communicator
                                              ! will open a connection.  This communicator must have been added to
                                              ! mpi_proc_data with mpi_proc_data.mpi_env_comm_add().  If not present,
                                              ! MPI_COMM_WORLD will be used.
    logical, optional, intent(in) :: barrier  ! for an mpi call, an mpi_barrier will be called before exiting
                                              ! this subroutine if barrier=.true.
    logical, optional, intent(in) :: allcpus  ! .true. if command should be executed for all the cpus
                                              ! .false. if only the io-node cpus should execute the command (default)
    integer, optional, intent(in) :: method   ! if present, use this in place of EXEC_METHOD -- only valid if method=0,1

    character(len=48) :: myshell    ! shell as used here

    integer :: ier             ! returned status code of shell command
    integer :: asmethod        ! exec method to use

    character(len=EXEC_MAX_STR) :: arg(3)  ! exec arguments

    ier = 0

    call enterLog("jsystem_env")

    if (len_trim(adjustl(str))>(EXEC_MAX_STR-1)) then
       call errorLog("jsystem_env: shell command string is too long")
       ier=123 ; goto 900
    end if

    asmethod = EXEC_METHOD
    if (present(method)) then
       if (method==0 .or. method==1) asmethod=method
    end if

    if (asmethod/=1) then
       myshell = adjustl(SHELL_EXEC)           ! using execvp(), path is searched
    else
       myshell = "/bin/"//adjustl(SHELL_EXEC)  ! using execv(), full path to executable is required
    end if

    myshell = adjustl(myshell)
    arg(1)  = adjustl(SHELL_EXEC)
    arg(2)  = adjustl(SHELL_LAUNCH)
    arg(3)  = adjustl(str)


    call infoLog("jsystem_env: executing system call, "//trim(adjustl(str)))

    ier = jexec_env(myshell, 3, arg, ne, env, comm, barrier, allcpus, method, noelog=.true.)

900 continue

    call exitLog("jsystem_env")
  end function jsystem_env

  
  !
  ! --------- jsystem_one ----------
  !
  ! Make a system call with the default environment
  !
  ! For an MPI environment, the communicator (default MPI_COMM_WORLD) must have been added to mpi_proc_data
  ! and only I/O ranks will execute the system call.
  !
  function jsystem_one(str, comm, barrier, allcpus, method) result(ier)
    character(len=*), intent(in) :: str            ! shell command to execute

    integer, optional, intent(in) :: comm     ! communicator used as basis for opening a connection.  Only the
                                              ! unique ranks on each node associated with this communicator
                                              ! will open a connection.  This communicator must have been added to
                                              ! mpi_proc_data with mpi_proc_data.mpi_env_comm_add().  If not present,
                                              ! MPI_COMM_WORLD will be used.
    logical, optional, intent(in) :: barrier  ! for an mpi call, an mpi_barrier will be called before exiting
                                              ! this subroutine if barrier=.true.
    logical, optional, intent(in) :: allcpus  ! .true. if command should be executed for all the cpus
                                              ! .false. if only the io-node cpus should execute the command (default)
    integer, optional, intent(in) :: method   ! if present, use this in place of EXEC_METHOD -- only valid if method=0,1

    integer :: ier    ! returned status code of shell command

    character(len=8) :: env(0)  ! no environment

    ier = jsystem_env(str,-1,env,comm=comm, barrier=barrier, allcpus=allcpus, method=method)
  end function jsystem_one
  

  !
  ! -------------------- jexec_env --------------------
  ! Execute a command on the path in a forked process.
  ! Generally, the first argument should be the name of the command being executed.
  !
  ! For an MPI environment, the communicator (default MPI_COMM_WORLD) must have been added to mpi_proc_data
  ! and only I/O ranks will execute the system call.
  !
  function jexec_env(exe, na, arg, ne, env, comm, barrier, allcpus, method, noelog) result(ier)
    use iso_c_binding, only: c_char, c_null_char
    use mpi_proc_data, only: check_mpi_io_node, mpi_proc_info_comm

#ifdef __MPI
    use mpi, only: mpi_comm_world, mpi_status_size, mpi_integer, mpi_max

    integer, parameter :: STATUS_SIZE = MPI_STATUS_SIZE
#else
    integer, parameter :: STATUS_SIZE = 1
#endif

    character(len=*), intent(in) :: exe            ! command to execute, will be searched on path
    integer,          intent(in) :: na             ! number of arguments
    character(len=*), intent(in) :: arg(max(0,na)) ! list of arguments
    integer,          intent(in) :: ne             ! number of environment variables or -1 to use the default
    character(len=*), intent(in) :: env(max(0,ne)) ! list of new NAME=VALUE environment to replace current environment

    integer, optional, intent(in) :: comm     ! communicator used as basis for opening a connection.  Only the
                                              ! unique ranks on each node associated with this communicator
                                              ! will open a connection.  This communicator must have been added to
                                              ! mpi_proc_data with mpi_proc_data.mpi_env_comm_add().  If not present,
                                              ! MPI_COMM_WORLD will be used.
    logical, optional, intent(in) :: barrier  ! for an mpi call, an mpi_barrier will be called before exiting
                                              ! this subroutine if barrier=.true.
    logical, optional, intent(in) :: allcpus  ! .true. if command should be executed for all the cpus
                                              ! .false. if only the io-node cpus should execute the command (default)
    integer, optional, intent(in) :: method   ! if present, use this in place of EXEC_METHOD -- only valid if method=0,1
    logical, optional, intent(in) :: noelog   ! if present and .false. do not execute enterLog,exitLog,warnLog

    integer :: ier    ! returned status code of shell command

    integer :: ks,ka,ke        ! length of strings

    integer, external :: c_execsystem      ! c code

    character(kind=c_char), allocatable :: exes(:)   ! c version of exec name
    character(kind=c_char), allocatable :: args(:)   ! arguments for exec
    character(kind=c_char), allocatable :: envs(:)   ! enviornment variables

    logical :: try_exec      ! .true. if will try an exec call 
    logical :: try_multi     ! .true. if forwarding system calls to the shell server cpu
    logical :: ielog         ! .true. for enter,exit log messages

    integer :: use_comm      ! communicator to use
    integer :: iber          ! barrier error
    integer :: em            ! bounded exec method
    integer :: asmethod      ! exec method to use
    integer :: mynode        ! node number
    integer :: num_procs     ! number of processors in comm
    integer :: key           ! controls ordering of cpus in splitcomm
    integer :: splitcomm     ! new communicator among cpus on each node
    integer :: mysplit       ! cpu number on splitcomm
    integer :: cer           ! nonzero if splitcomm not ordered correctly
    integer :: aer           ! reduced cer
    integer :: numsplits     ! number of cpus on splitcomm
    integer :: maxna, maxne  ! maximum number of na,ne across all processors
    integer :: maxstr        ! number of strings in imash
    integer :: maxint        ! number of integers needed in imash
    integer :: i,k,ic,j      ! loop variables

    integer :: ixna,ixne,ixcd,ixexe,ixarg,ixenv  ! indices into imash

    integer, allocatable :: imash(:)   ! mash up of arguments and cwd to transfer by mpi    
    integer, allocatable :: iall(:,:)  ! mash up of arguments from each cpu

    character(len=EXEC_MAX_STR) :: mycwd  ! cwd of this cpu

    integer :: fromna, fromne  ! na,ne from other cpu

    character(len=EXEC_MAX_STR) :: fromcwd, fromexe   ! cwd, exe from other cpu

    character(len=EXEC_MAX_STR),allocatable :: fromarg(:), fromenv(:)  ! arg and env from other cpu

    integer :: status(STATUS_SIZE)  ! status from mpi call

    ier = 0
    
    ielog = .true.
    if (present(noelog)) then
       ielog = .not. noelog
    end if

    if (ielog) call enterLog("jexec_env")

    try_exec  = .true.
    try_multi = .false.
    asmethod  = EXEC_METHOD
    mycwd     = ' '            ! considered set if the first character is not a space

    if (present(method)) then
       if (method==0 .or. method==1) asmethod=method
    end if

#ifdef __MPI
    if (present(comm)) then
       use_comm = comm
    else
       use_comm = MPI_COMM_WORLD   
    end if

    try_exec = check_mpi_io_node(use_comm,ier)  ! is this an I/O node?

    if (ier/=0) then
       call errorLog("jexec_env: unknown communicator =",use_comm)
       ier=1
       try_exec = .false.
    else
       if (present(allcpus)) then
          if (asmethod==2 .and. allcpus) then
             ! -- shell server and using all cpus --
             try_exec  = .false.
             try_multi = .true.
          else
             try_exec = try_exec .or. allcpus
          end if
       end if
    end if
#endif

    !
    ! execution of the system call from this cpu
    !
    if (try_exec) then
      if (asmethod == 3) then
        !write(6,*) "System call with Fortran execute_command_line"
        !write(6,*) trim(arg(na)) ; call flush(6)
        call execute_command_line(arg(na),exitstat=ier)
      else
        ks = len_trim(adjustl(exe))

        ka = 0
        do i=1, na
          ka = ka + len_trim(adjustl(arg(i)))+1
        end do
       
        ke = 0
        if (ne>=0) then
          do i=1, ne
            ke = ke + len_trim(adjustl(env(i)))+1
          end do
        end if
       
        allocate(exes(ks+1),args(ka+1),envs(ke+1))  ! room for null termination
       
        exes = c_null_char
        args = c_null_char
        envs = c_null_char
       
        call cstring(trim(adjustl(exe)), exes, '2C')
       
        ka = 1
        do i=1, na
          call cstring(trim(adjustl(arg(i))), args(ka:),'2C')
          ka = ka+len_trim(adjustl(arg(i)))+1
        end do
       
        if (ne>=0) then
          ke = 1
          do i=1, ne
            call cstring(trim(adjustl(env(i))), envs(ke:),'2C')
            ke = ke+len_trim(adjustl(env(i)))+1
          end do
        end if

        call infoLog("jexec_env: executing exec call, "//trim(exe))

        em = max(0, min(2, asmethod))

        ier = c_execsystem(em, EXEC_VERBOSE, exes, na, args, ne, envs)
       
        deallocate(exes, args, envs)
      end if
    else if (.not. try_multi) then
      call infoLog("jexec_env: not an I/O cpu -- no exec()")
    end if

#ifdef __MPI
    !
    ! forward system calls to the cpu on the node connected to the shell server
    !
    if (try_multi) then
       if (asmethod/=2) then
          call errorLog("jexec_env: what am I doing here? -- asmethod=",asmethod)  
          ier=1 ; goto 900
       end if

       !
       ! create a communicator for each node with the shell server cpu at rank 0
       !
       call mpi_proc_info_comm(use_comm, ier, mynode=mynode, num_procs=num_procs)
       if (ier/=0) then
          call errorLog("jexec_env: mpi_proc_info_comm error on use_comm")
          ier=1 ; goto 900
       end if

       if (EXEC_IS_SERVER) then
          key = 0    ! connected to shell server, want this as rank 0
       else
          key = 1   
       end if

       call mpi_comm_split(use_comm, mynode, key, splitcomm, ier)
       if (ier/=0) then
          call errorLog("jexec_env: mpi_comm_split error")
          ier=1 ; goto 900
       end if

       !
       ! have splitcomm, exit through 800 
       !

       ! -- check how ranks are set up --
       call mpi_comm_rank(splitcomm, mysplit, ier)
       if (ier/=0) then
          call errorLog("jexec_env: mpi_comm_rank error on splitcomm")
          ier=1 ; goto 800
       end if

       if (mysplit==0 .and. .not. EXEC_IS_SERVER) then
          call errorLog("jexec_env: cpu is rank 0 in splitcomm but not connected to a shell server")
          cer=1
       else
          cer = 0
       end if
       
       call mpi_allreduce(cer, aer, 1, MPI_INTEGER, MPI_MAX, use_comm, ier)   ! collect error codes from all cpus
       if (ier/=0) aer=1
       if (aer/=0) then
          call errorLog("jexec_env: mpi_allreduce error or there is a rank 0 cpu not connected to a shell server")
          ier=1 ; goto 800
       end if
       
       ! -- mash up input data into imash --
       call mpi_comm_size(splitcomm, numsplits, ier)
       if (ier/=0) then
          call errorLog("jexec_env: mpi_comm_size error on splitcomm")
          ier=1 ; goto 800
       end if

       call getcwd(mycwd)    ! fetch working directory
       mycwd = adjustl(mycwd)  ! so first character is not a space to mark as set

       if (ier==0) call mpi_allreduce(na, maxna, 1, MPI_INTEGER, MPI_MAX, splitcomm, ier)
       if (ier==0) call mpi_allreduce(ne, maxne, 1, MPI_INTEGER, MPI_MAX, splitcomm, ier) 
       if (ier/=0 .or. maxna<0) then
          call errorLog("jexec_env: mpi_allreduce error on maxna or maxne or maxna<0")
          ier=1 ; goto 800
       end if

       maxne = max(0,maxne)                 ! maxne<0 possible if no environment variables

       maxstr = maxna + maxne + 2           ! number of strings, args, envs, exe, cwd
       maxint = 2 + EXEC_MAX_STR*maxstr     ! number of integers needed to store arguments -- waste of bytes

       allocate(imash(maxint))
       imash=0                              ! so all strings are zero terminated

       ixna  = 1                            ! set up pointers into imash, actual pointers would have been more clever
       ixne  = 2
       ixcd  = 3
       ixexe = ixcd  + EXEC_MAX_STR
       ixarg = ixexe + EXEC_MAX_STR
       ixenv = ixarg + maxna*EXEC_MAX_STR

       imash(ixna) = na
       imash(ixne) = ne

       if (len_trim(mycwd)>(EXEC_MAX_STR-1)) then                           ! need space for 0 terminator
          call errorLog("jexec_env: working directory string is too long")
          ier=1 ; goto 800
       end if
          
       ic = ixcd
       do i = 1, len_trim(mycwd)
          imash(ic) = ichar(mycwd(i:i)) ; ic=ic+1
       end do
       
       if (len_trim(exe)>(EXEC_MAX_STR-1)) then
          call errorLog("jexec_env: exe string is too long")
          ier=1 ; goto 800
       end if
          
       ic = ixexe
       do i = 1, len_trim(exe)
          imash(ic) = ichar(exe(i:i)) ; ic=ic+1
       end do

       do k=1, na
          if (len_trim(arg(k))>(EXEC_MAX_STR-1)) then
             call errorLog("jexec_env: argument string is too long")
             ier=1 ; goto 800
          end if
          
          ic = ixarg + (k-1)*EXEC_MAX_STR
          do i = 1, len_trim(arg(k))
             imash(ic) = ichar(arg(k)(i:i)) ; ic=ic+1
          end do
       end do

       do k=1, ne
          if (len_trim(env(k))>(EXEC_MAX_STR-1)) then
             call errorLog("jexec_env: environment string is too long")
             ier=1 ; goto 800
          end if
          
          ic = ixenv + (k-1)*EXEC_MAX_STR
          do i = 1, len_trim(env(k))
             imash(ic) = ichar(env(k)(i:i)) ; ic=ic+1
          end do
       end do

       ! -- send imash to rank 0 --
       allocate(iall(maxint,numsplits))

       call infoLog("jexec_env: sending command to rank 0 of split communicator")

       call mpi_gather(imash, maxint, mpi_integer, iall, maxint ,mpi_integer, 0, splitcomm, ier)
       if (ier/=0) then
          call errorLog("jexec_env: mpi_gather error")
          ier=1 ; goto 800
       end if

       cer = 0   ! system call error flag

       if (mysplit==0) then
          !
          ! -- execute system commands in iall --
          !
          allocate(fromarg(maxna), fromenv(maxne))

          call infoLog("jexec_env: starting system calls on behalf of siblings")
          do j=1, numsplits
             !
             ! exits through 700 with ier set
             !
             call infoLog("jexec_env: exec for cpu ",j-1)

             imash = iall(:,j)

             fromna = imash(ixna)
             fromne = imash(ixne)

             if (fromna>maxna .or. fromne>maxne .or. fromna<0) then
                call errorLog("jexec_env: bad na,ne from sibling ranks")   ! assert
                ier=1 ; goto 700
             end if
                
             call tostr(imash(ixcd),  fromcwd) ; if (ier/=0) goto 700
             call tostr(imash(ixexe), fromexe) ; if (ier/=0) goto 700

             call infoLog("jexec_env: cd to '"//trim(adjustl(fromcwd))//"'")

             call sset_cwd(trim(adjustl(fromcwd)),ier)

             if (ier/=0) then
                call errorLog("jexec_env: error changing directory")  
                call errorLog(trim(adjustl(fromcwd))) 
                ier=1 ; goto 700
             end if
                
             do k=1, fromna
                call tostr(imash(ixarg+(k-1)*EXEC_MAX_STR), fromarg(k)) ; if (ier/=0) goto 700
             end do

             do k=1, fromne
                call tostr(imash(ixenv+(k-1)*EXEC_MAX_STR), fromenv(k)) ; if (ier/=0) goto 700
             end do

             ks = len_trim(adjustl(fromexe))

             ka = 0
             do i=1, fromna
                ka = ka + len_trim(adjustl(fromarg(i)))+1
             end do
       
             ke = 0
             if (fromne>=0) then
                do i=1, fromne
                   ke = ke + len_trim(adjustl(fromenv(i)))+1
                end do
             end if
             
             allocate(exes(ks+1),args(ka+1),envs(ke+1))  ! room for null termination
             
             exes = c_null_char
             args = c_null_char
             envs = c_null_char
             
             call cstring(trim(adjustl(fromexe)), exes, '2C')
             
             ka = 1
             do i=1, fromna
                call cstring(trim(adjustl(fromarg(i))), args(ka:),'2C')
                ka = ka+len_trim(adjustl(fromarg(i)))+1
             end do
             
             if (fromne>=0) then
                ke = 1
                do i=1, fromne
                   call cstring(trim(adjustl(fromenv(i))), envs(ke:),'2C')
                   ke = ke+len_trim(adjustl(fromenv(i)))+1
                end do
             end if
             
             call infoLog("jexec_env: executing exec call, "//trim(fromexe))
             
             aer = c_execsystem(asmethod, EXEC_VERBOSE, exes, fromna, args, fromne, envs)
                   
             call infoLog("jexec_env: status of system call =", aer)

             if (j==1) then
                cer = aer   ! record rank 0 result
             else
                call infoLog("jexec_env: sending status to cpu ",j-1)
                call mpi_send(aer, 1, mpi_integer, j-1, 0, splitcomm, ier)
                if (ier/=0) then
                   call errorLog("jexec_env: mpi_send error using splitcomm to cpu ",j-1)
                   ier=0 ; cer=123    ! mark as a rank 0 error
                end if
             end if

             deallocate(exes, args, envs)
             cycle  ! ok

             ! -- loop error exit --
700          continue
             if (j==1) then
                call errorLog("jexec_env: eror setting up system call on rank 0")
                cer = ier
                ier = 0
             else
                call errorLog("jexec_env: eror setting up sibling system call, sending status to cpu ",j-1)
                call mpi_send(ier, 1, mpi_integer, j-1, 0, splitcomm, aer)  ! send loop error to sibling
                if (aer/=0) then
                   call errorLog("jexec_env: mpi_send error after error using splitcomm to cpu ",j-1)
                   ier=0 ; cer=123    ! mark as a rank 0 error
                end if
             end if
          end do
          call infoLog("jexec_env: done with calls on behalf of siblings")
       else
          !
          ! -- siblings just get back the status code --
          !
          call infoLog("jexec_env: waiting for status from rank 0")
          call mpi_recv(cer, 1, mpi_integer, 0, 0, splitcomm, status, ier)
          if (ier/=0) then
             call errorLog("jexec_env: mpi_recv error using splitcomm")
             ier=1 ; goto 800
          end if
          call infoLog("jexec_env: received status =",aer)
       end if

       ier = cer  ! made it through ok, set subroutine error flag to system call error
          
       ! -- done with splitcomm --
800    continue
       if (mysplit==0 .and. mycwd(1:1)/=' ') then
          call sset_cwd(trim(mycwd),aer)  ! restore cwd after being changed to sibling cwd
          if (aer/=0) then
             call errorLog("jexec_env: error restoring rank 0 cwd")
             ier=1 
          end if
       end if

       if (allocated(fromarg)) deallocate(fromarg,fromenv)
       if (allocated(iall))    deallocate(iall)
       if (allocated(imash))   deallocate(imash)

       call mpi_comm_free(splitcomm, aer)
       if (aer/=0) then
          call errorLog("jexec_env: mpi_comm_free error")
          ier=1 
       end if
    end if
#endif

900 continue
    
#ifdef __MPI
    if (present(barrier)) then
       if (barrier) then
          call MPI_BARRIER(use_comm, iber)
          if (iber/=0) then
             call errorLog("jexec_env: barrier error")
             ier=1
          end if
       end if
    end if
#endif
    
    if (ielog) call exitLog("jexec_env")

  contains
    !
    ! --------- tostr ---------
    ! extracts the string in the integer array
    !
    subroutine tostr(iin, str)
      integer,       intent(in)  :: iin(*)   ! part of imash, max length EXEC_MAX_STR including terminating 0
      character(len=*), intent(out) :: str      ! place string here

      integer :: m   ! loop variable

      ier = 0
      str = ' '

      do m=1, EXEC_MAX_STR
         if (iin(m)==0) return
         str(m:m) = char(iin(m))
      end do

      call errorLog("jexec_env.tostr: did not find terminating 0 in integer representation of string")
      ier=1
    end subroutine tostr
  end function jexec_env
  
  !
  ! --------- jexec_one ----------
  !
  ! Make an exec call with the default environment.
  ! Generally, the first argument should be the name of the command being executed.
  !
  ! For an MPI environment, the communicator (default MPI_COMM_WORLD) must have been added to mpi_proc_data
  ! and only I/O ranks will execute the system call.
  !
  function jexec_one(exe, na, arg, comm, barrier, allcpus, method) result(ier)
    character(len=*), intent(in) :: exe            ! command to execute, will be searched on path
    integer,          intent(in) :: na             ! number of arguments
    character(len=*), intent(in) :: arg(na)        ! list of arguments

    integer, optional, intent(in) :: comm     ! communicator used as basis for opening a connection.  Only the
                                              ! unique ranks on each node associated with this communicator
                                              ! will open a connection.  This communicator must have been added to
                                              ! mpi_proc_data with mpi_proc_data.mpi_env_comm_add().  If not present,
                                              ! MPI_COMM_WORLD will be used.
    logical, optional, intent(in) :: barrier  ! for an mpi call, an mpi_barrier will be called before exiting
                                              ! this subroutine if barrier=.true.
    logical, optional, intent(in) :: allcpus  ! .true. if command should be executed for all the cpus
                                              ! .false. if only the io-node cpus should execute the command (default)
    integer, optional, intent(in) :: method   ! if present, use this in place of EXEC_METHOD -- only valid if method=0,1

    integer :: ier    ! returned status code of command

    character(len=8) :: env(0)  ! no environment

    ier = jexec_env(exe,na,arg,-1,env,comm=comm,barrier=barrier,allcpus=allcpus,method=method)
  end function jexec_one
  
end module execsystem
