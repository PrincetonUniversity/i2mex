
!
! A set of utilities for writing to a log file
!
! example of usage:
!
!    integer :: i = 47
!    real*8  :: d = 3.14159e17
!
!    call openLog("mylog.log",icpu=7,numcpu=16)     ! if not opened, logging goes to portlib_<pid>.log
!
!    call setLogLevel(INFO_LEVEL)  ! default log level = WARN_LEVEL to show warnings and errors only
!
!    call infoLog("just an info message")
!    call infoLog("just an info message with i = ",i)
!    call infoLog("just an info message with d = ",d)
!
!    call warnLog("just a warn message")
!    call sleep(1)
!    call warnLog("just a warn message with i = ",i)
!    call warnLog("just a warn message with d = ",d)
!
!    call errorLog("just an error message")
!    call errorLog("just an error message with i = ",i)
!    call sleep(1)
!    call errorLog("just an error message with d = ",d)
!   
!    call warnLog("enable buffering and no time stamps")
!
!    if (logLevel()>=INFO_LEVEL) then
!      call setLogBuffered(.true.)     ! do not flush after each log message
!      call setLogTime(.false.)        ! do not time stamp log messages
!
!      call infoLog("some more detailed messages")
!
!      write(logline,'(a,i2,a,e12.4)') 'vec(',i,') = ',d
!      call infoLog(logline)
!
!      call setLogTime(.true.)
!      call setLogBuffered(.false.)     ! flush buffers
!    end if
!   
!    call warnLog("time stamps enabled")
!
! with output in mylog_07.log,
!    
!    23Jul2009 12:32:17.26 [07] -openLog: PID=17767
!    23Jul2009 12:32:17.26 [07] %just an info message
!    23Jul2009 12:32:17.26 [07] %just an info message with i =47
!    23Jul2009 12:32:17.26 [07] %just an info message with d =3.141593e+17
!    23Jul2009 12:32:17.26 [07] !just a warn message
!    23Jul2009 12:32:18.26 [07] !just a warn message with i =47
!    23Jul2009 12:32:18.26 [07] !just a warn message with d =3.141593e+17
!    23Jul2009 12:32:18.26 [07] ?just an error message
!    23Jul2009 12:32:18.26 [07] ?just an error message with i =47
!    23Jul2009 12:32:19.26 [07] ?just an error message with d =3.141593e+17
!    23Jul2009 12:32:19.26 [07] !enable buffering and disable time stamps
!    [07] %some more detailed messages
!    [07] %vec(47) =   0.3142E+18
!    23Jul2009 12:32:19.26 [07] !time stamps enabled
!
! Mods
!    22Dec2009    Jim.Conboy@cffe.ac.uk
!                 cstring(trim(),  => cstr_2C(trim(), 
!                 [ Former invalid under debug compilation ]
!
!    16Dec2009    Jim.Conboy@cffe.ac.uk
!                 Added errLogOpen to log file open errors
!----------------------------------------------------------------------------
module logmod
#ifdef __MPI
  use mpi
#endif
  integer, parameter, public :: INFO_LEVEL  = 0   ! log level for info messages
  integer, parameter, public :: WARN_LEVEL  = 1   ! log level for warning messages
  integer, parameter, public :: ERROR_LEVEL = 2   ! log level for error messages
  integer, parameter, public :: NOMSG_LEVEL = 3   ! log level for no messages

  character*132, save :: logline   ! a utility string

  ! --- internal ---
  interface infoLog
     module procedure infoLogMsg, infoLogInt, infoLogFloat, &
                      infoLogMsg2,infoLogInt2,infoLogFloat2, &
                      infoLogLogical
  end interface

  interface warnLog
     module procedure warnLogMsg, warnLogInt, warnLogFloat
  end interface

  interface errorLog
     module procedure errorLogMsg, errorLogInt, errorLogFloat
  end interface

  interface setLogLevel
     module procedure setLogLevelInt, setLogLevelString
  end interface
  
  private
  public logline, infoLog, warnLog, errorLog
  public openLog, closeLog, getLogName, setLogLevel, setLogBuffered, logLevel, setLogTime
  public errLogOpen, errLogAlloc
  public enterLog, exitLog, isOpenLog

contains
  !
  ! ------ openLog ------
  ! Open a log file.  The status flag is similar to fortran's open
  ! statment.  A log file will be automatically opened to portlib_<processID>.log
  ! if a message is written but a log file is not currently open.
  !
  !   status = 'ONCE'    -> do nothing if already open otherwise like 'REPLACE' (default)
  !          = 'NEW'     -> new file created, error if file already exists
  !          = 'OLD'     -> reopen an existing file and append new messages
  !          = 'REPLACE' -> if old file exists, delete and open a new file
  !          = 'RENAME' ->  if old file exists, rename old file as <file.name>~  and open a new file
  !
  ! If there is an error and iostat is not present, a message will be written
  ! to stderr and the error will be ignored.
  !
  ! If ICPU>=0 is given, then the file name will be modified to accommodate the cpu index
  ! and the cpu index will be included in each message which is written.  As an example
  !   FILE='foo.log', ICPU=7, NUMCPU=16  will create the log file 'foo_07.log'
  !
  ! After opening a new log file, a header line will be printed with the flag '-' and
  ! it will contain the process id, such as
  !    23Jul2009 13:43:02.09 -openLog: PID=18475
  !
  ! The LOG_LEVEL environment variable will be used to set the logging level in openLog
  ! unless getenv=.false. is present,
  !    LOG_LEVEL=0   -> info
  !    LOG_LEVEL=1   -> warning
  !    LOG_LEVEL=2   -> error
  !    LOG_LEVEL=3   -> no log messages except the header line
  !
  subroutine openLog(file,status,iostat,ic_myid,comm_world,icpu, numcpu,getenv)
    use iso_c_binding, only: c_char
    character*(*),          intent(in)  :: file    ! name of the log file
    character*(*),optional, intent(in)  :: status  ! open status flag
    integer,      optional, intent(out) :: iostat  ! 0->ok, nonzero->error
    integer,      optional, intent(out) :: ic_myid ! return myid# incorporated into the log file name
    integer,      optional, intent(in)  :: comm_world ! communicator
    integer,      optional, intent(in)  :: icpu    ! index of cpu or -1 if no cpus -- this is indexed from 0
    integer,      optional, intent(in)  :: numcpu  ! maximum number of cpus
    logical,      optional, intent(in)  :: getenv  ! if .true., then fetch the environment variable LOG_LEVEL
                                                   ! and use it to set the initial logging level (default .true.)

    character*10 :: lstatus  ! lower case status

    character(kind=c_char) :: name(len(file)+1) ! C file name
    character(kind=c_char) :: st(2)             ! single character status for C

    integer :: ier      ! error flag
    integer :: ig       ! nonzero to set log level from LOG_LEVEL environment variable
    integer :: ic, mc   ! cpu index and maximum number of cpus
    integer :: mpierr
    integer :: myid 

    integer, SAVE :: numprocs = 0
    integer, SAVE :: icpu_get    ! index of cpu or -1 if no cpus -- this is indexed from 0
    integer, SAVE :: numcpu_get  ! maximum number of cpus when icpu>=0
   
    if (present(comm_world)) then !ignore icpu and numcpu
       if(numprocs.eq.0) then
#ifdef __MPI
          call MPI_COMM_RANK(comm_world,myid,mpierr)
          call MPI_COMM_SIZE(comm_world,numprocs,mpierr)
#else
          numprocs=1
          myid=0
#endif

          numcpu_get= numprocs
          icpu_get=myid

       endif

       ic=icpu_get
       
       mc=numcpu_get

    else
       ic = -1
       if (present(icpu)) ic=icpu
       
       mc = -1
       if (present(numcpu)) mc=numcpu
       
    endif
    ier=0

    lstatus = 'once'
    if (present(status)) lstatus = status

    ig=1
    if (present(getenv)) then
       if (.not. getenv) ig=0
    end if

    call ulower(lstatus)
    if (lstatus=='once') then
       call cstr_2C('1', st, '2C') 
    else if (lstatus=='new') then
       call cstr_2C('n', st, '2C')
    else if (lstatus=='old') then
       call cstr_2C('a', st, '2C')
    else if (lstatus=='replace') then
       call cstr_2C('r', st, '2C')
    else if (lstatus=='rename') then
       call cstr_2C('m', st, '2C')
    else
       print '(a)',"?openLog: unknown status option = "//trim(lstatus)
       ier=1 ; goto 100
    end if

    call cstr_2C(trim(file),name,'2C')
       
    call c_openlog(name(1), st(1), ier,comm_world ,ic, mc, ig)

100 continue
    if (present(iostat)) iostat=ier
    if (present(ic_myid)) ic_myid=ic
  end subroutine openLog

  !
  ! ----- closeLog -----
  ! Close an open log file.
  !
  ! If there is an error and iostat is not present, a message will be written
  ! to stderr and the error will be ignored.
  !
  subroutine closeLog(iostat)
    integer,      optional, intent(out) :: iostat  ! 0->ok, nonzero->error
    
    integer :: ier   ! error flag

    ier=0

    call c_closelog(ier)
    if (present(iostat)) iostat=ier
  end subroutine closeLog

  !
  ! ------- isOpenLog --------
  ! Return true if the log file is open
  !
  function isOpenLog() result(isopen)
    logical :: isopen   ! true if log file is open

    integer, external :: c_isopenlog

    isopen = c_isopenlog()/=0
  end function isOpenLog

  !
  ! ------ getLogName -----
  ! Return the name of the open log file.
  !
  subroutine getLogName(file, ier) 
    use iso_c_binding, only: c_char
    character*(*), intent(out) :: file   ! returned file name
    integer,       intent(out) :: ier    ! nonzero on error
                                         !   -1 -> no open log file
                                         !   -2 -> file name was too small

    character(kind=c_char) :: name(len(file)+1) ! C file name

    file = " "

    call c_getlogname(size(name), name, ier)
    if (ier<0) return
    ier=0
       
    call cstring(file, name, '2F')
    file = adjustl(trim(file))
  end subroutine getLogName

  !
  ! ------ setLogLevelString ------
  ! Set the level of logging by name.  The levels should be used as,
  !
  ! level --  "INFO","WARN","ERR","NOMSG" or "0", "1", "2", "3"
  !
  !   0 -> info level, detailed information not to be used for production runs
  !   1 -> warn level (default), typically used to mark entrance and exit
  !        of major subroutines or other significant messages
  !   2 -> error level, used to write an error before aborting the execution
  !   3 -> no log messages
  !
  subroutine setLogLevelString(level)
    character*(*),optional, intent(in) :: level  ! level as a string

    character*10 :: l_level  ! lower case level
    integer      :: ilevel   ! 0->info, 1->warn, 2->error, 3->no message

    ilevel = 1

    if (present(level)) then
       l_level = adjustl(level)
       call ulower(l_level)

       if (len_trim(l_level)>2) then
          if (l_level(1:3)=='inf') ilevel=0
          if (l_level(1:3)=='war') ilevel=1
          if (l_level(1:3)=='err') ilevel=2
          if (l_level(1:3)=='nom') ilevel=3
       else if (len_trim(l_level)>0) then
          if (l_level(1:1)=='0') ilevel=0
          if (l_level(1:1)=='1') ilevel=1
          if (l_level(1:1)=='2') ilevel=2
          if (l_level(1:1)=='3') ilevel=3
       end if
    end if
    call c_setloglevel(ilevel)
  end subroutine setLogLevelString

  !
  ! ------ setLogLevelInt ------
  ! Set the level of logging by integer.  The levels should be used as,
  !
  !   0 -> info level, detailed information not to be used for production runs
  !   1 -> warn level (default), typically used to mark entrance and exit
  !        of major subroutines or other significant messages
  !   2 -> error level, used to write an error before aborting the execution
  !   3 -> no log messages
  !
  subroutine setLogLevelInt(klevel)
    integer,intent(in) :: klevel  ! 0->info, 1->warn, 2->error, 3->no message

    integer :: ilevel 
    
    ilevel = min(3,max(0,klevel))

    call c_setloglevel(ilevel)
  end subroutine setLogLevelInt

  !
  ! ------- logLevel ------
  ! Return the logging level
  !
  integer function logLevel()
    integer, external :: c_loglevel

    logLevel = c_loglevel()
  end function logLevel

  !
  ! ------- setLogBuffered -----
  ! Change whether the log file is being buffered.  The default is unbuffered.
  !
  subroutine setLogBuffered(ibuff)
    logical, intent(in) :: ibuff   ! .true. for buffered output, .false. for unbuffered
    integer :: ib                  ! integer logical

    ib=0
    if (ibuff) ib=1 

    call c_setlogbuffered(ib)
  end subroutine setLogBuffered

  !
  ! ------- setLogTime -----
  ! Change whether the log file time stamps the messages.  Default is true.
  !
  subroutine setLogTime(itime)
    logical, intent(in) :: itime   ! .true. for time stamped output
    integer :: it                  ! integer logical

    it=0
    if (itime) it=1 

    call c_setlogtime(it)
  end subroutine setLogTime

  !
  ! ----- enterLog ------
  ! Write out a warn level log message of the form
  ! ... \Enter\ msg
  !
  ! This is used for tagging entry to a subroutine or code block and should
  ! be paired with an exitLog(msg)
  !
  subroutine enterLog(msg)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg   ! subroutine or code block name

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (WARN_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_enterlog(cmsg)
  end subroutine enterLog

  !
  ! ----- exitLog ------
  ! Write out a warn level log message of the form
  ! ... /Exit/ msg
  !
  ! This is used for tagging entry to a subroutine or code block and should
  ! be paired with an enterLog(msg)
  !
  subroutine exitLog(msg)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg   ! subroutine or code block name

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (WARN_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_exitlog(cmsg)
  end subroutine exitLog

  !
  ! ----- infoLog ------
  ! Write out an info level log message optionally with an integer (format I9)
  ! or double (format ES24.16) value.
  !
  ! INFO_LEVEL messages are tagged with the flag '%'
  !
  subroutine infoLogMsg2(cSr, cmsg)
    character*(*), intent(in) :: cSr    ! log message
    character*(*), intent(in) :: cmsg   ! log message
  !
    call infoLogMsg( cSr//' -I- '//cmsg )
    return
  end subroutine infoLogMsg2

  subroutine infoLogMsg(msg)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg   ! log message

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (INFO_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_infolog(cmsg)
  end subroutine infoLogMsg

  subroutine infoLogInt2(cSr, cmsg, ivalue )
    character*(*), intent(in) :: cSr      ! log message
    character*(*), intent(in) :: cmsg     ! log message
    integer,       intent(in) :: ivalue   ! integer value
  !
    call infoLogInt( cSr//' -I- '//cmsg, ivalue )
    return
  end subroutine infoLogInt2

  subroutine infoLogInt(msg, ivalue)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg     ! log message
    integer,       intent(in) :: ivalue  ! integer value

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (INFO_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_infologi(cmsg,ivalue)
  end subroutine infoLogInt

  subroutine infoLogLogical(msg, lvalue)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg     ! log message
    logical,       intent(in) :: lvalue  ! logical value
    integer                   :: ivalue

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (INFO_LEVEL<c_loglevel()) return

    ivalue = 0
    if (lvalue) ivalue = 1

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_infologi(cmsg,ivalue)
  end subroutine infoLogLogical

  subroutine infoLogFloat2(cSr, cmsg, dvalue )
    character*(*), intent(in) :: cSr      ! log message
    character*(*), intent(in) :: cmsg     ! log message
    real*8,        intent(in) :: dvalue   ! double value
  !
    call infoLogFloat( cSr//' -I- '//cmsg, dvalue )
    return
  end subroutine infoLogFloat2

  subroutine infoLogFloat(msg, dvalue)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg     ! log message
    real*8,        intent(in) :: dvalue  ! double value

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (INFO_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_infologd(cmsg,dvalue)
  end subroutine infoLogFloat

  !
  ! ----- warnLog ------
  ! Write out a warn level log message optionally with an integer (format I9)
  ! or double (format ES24.16) value.
  !
  ! WARN_LEVEL messages are tagged with the flag '!'
  !
  subroutine warnLogMsg(msg)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg   ! log message

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (WARN_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_warnlog(cmsg)
  end subroutine warnLogMsg

  subroutine warnLogInt(msg, ivalue)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg     ! log message
    integer,       intent(in) :: ivalue  ! integer value

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (WARN_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_warnlogi(cmsg,ivalue)
  end subroutine warnLogInt

  subroutine warnLogFloat(msg, dvalue)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg     ! log message
    real*8,        intent(in) :: dvalue  ! double value

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (WARN_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_warnlogd(cmsg,dvalue)
  end subroutine warnLogFloat

  !
  ! ----- errorLog ------
  ! Write out an error level log message optionally with an integer (format I9)
  ! or double (format ES24.16) value.
  !
  ! ERROR_LEVEL messages are tagged with the flag '?'
  !
  subroutine errorLogMsg(msg)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg   ! log message

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (ERROR_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_errorlog(cmsg)
  end subroutine errorLogMsg

  subroutine errorLogInt(msg, ivalue)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg     ! log message
    integer,       intent(in) :: ivalue  ! integer value

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (ERROR_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_errorlogi(cmsg,ivalue)
  end subroutine errorLogInt

  subroutine errorLogFloat(msg, dvalue)
    use iso_c_binding, only: c_char
    character*(*), intent(in) :: msg     ! log message
    real*8,        intent(in) :: dvalue  ! double value

    integer, external :: c_loglevel       ! get current log level
    character(kind=c_char)      :: cmsg(len(msg)+1) ! C message

    if (ERROR_LEVEL<c_loglevel()) return

    call cstr_2C(trim(msg),cmsg,'2C')
    call c_errorlogd(cmsg,dvalue)
  end subroutine errorLogFloat

  subroutine errLogOpen( cr, iostat, cFile, ilun, Status, lFatal )
  !--   Log a file open error 
    implicit  none
    character(len=*),intent(in)        :: cr      ! calling routine
    integer, intent(in)                :: iostat  ! open return ; do nothing if == 0
    character(len=*),intent(in)        :: cFile   ! tried to open this file
    integer, intent(in)                :: ilun    !   ..  on unit
    character(len=*), &
       intent(in), optional            :: Status  !  of the file
    logical, intent(in), optional      :: lFatal  !  call bad_exit if this is set..
  !
    integer                            :: ifatal  !  {1:2} for notfatal|fatal}
    character(len=4),dimension(3)      :: cEF = (/' -E-',' -F-',' -I-'/)
    character(len=256)                 :: cbuf   &
                                         ,cmsg
  !
    if( iostat == 0 )                       then
       write(cbuf,'(5a,i4)') &
          cr,cEF(3),' opened ',trim(cFile),' on ',ilun
       return
    endif
  !
    ifatal = 1
    if( present(lFatal) )  then
       if( lFatal ) ifatal = 2
    endif
    write(cbuf,'(2a,1x,i6,3a,i4)') &
      cr,cEF(ifatal),iostat,' opening ',trim(cFile),' on ',ilun 
    if( present(status))  cbuf = cbuf//' status '//Status
  !
    call cstr_2C(trim(cbuf),cmsg,'2C')
    call c_errorlog(cmsg)
  !
    if( ifatal == 2 .and. iostat /= 0 )    call bad_exit()
  !
  end subroutine errLogOpen

  subroutine errLogAlloc( cr, istat, cArray, isize, lFatal )
  !--  Log an array allocation error
    implicit  none
    character(len=*),intent(in)        :: cr      ! calling routine
    integer, intent(in)                :: istat   ! allocate return ; do nothing if == 0
    character(len=*),intent(in)        :: cArray  ! tried to allocate this array
    integer, intent(in)                :: isize   !   ..  # words
    logical, intent(in), optional      :: lFatal  !  call bad_exit if this is set..
  !
    integer                            :: ifatal  !  {1:2} for notfatal|fatal}
    character(len=4),dimension(3)      :: cEF = (/' -E-',' -F-',' -I-'/)
    character(len=256)                 :: cbuf   &
                                         ,cmsg
  !
    if( istat == 0 )                     return
  !
    ifatal = 1
    if( present(lFatal) )  then
       if( lFatal ) ifatal = 2
    endif
    write(cbuf,'(2a,1x,i8,3a,i10,a)') &
      cr,cEF(ifatal),istat,' allocating ',trim(cArray),' (',isize,')'
  !
    call cstr_2C(trim(cbuf),cmsg,'2C')
    call c_errorlog(cmsg)
  !
    if( ifatal == 2 .and. istat /= 0 )    call bad_exit()
  !
  end subroutine errLogAlloc
end module logmod
