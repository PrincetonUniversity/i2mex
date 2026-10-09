module splitn_module

!     Modifications
!
!     1.23  05Apr2009 -- dmccune@pppl.gov -- splitn has been extensively
!                modified to support updatable namelist elements, i.e.
!                a portion of the TRANSP (trcom only) namelist can now
!                be read as an update, part way through a TRANSP run.
!                At the same time, trcom names maximum length was increased
!                from 16 to 32 characters; name declarations in splitn
!                lengthened accordingly.
!
!    -1.10  29Aug2007                   jim.conboy@jet.uk
!                                       Rolled back 1.10 ; see
!                                       msplitn_par.f90, splitn_io_list.f90
!
!     1.10  22Aug2007                   jim.conboy@jet.uk
!                                       Add splitn_io_list
!                                       Add nc_pre, nc_ext ( string length ) parameters
!                                       #ifdef __JET nc_ext -> 72
!----_^---------------------------------=========================================|
  implicit NONE
!
  save

  ! module to contain the representation of the TRANSP namelist with
  ! default settings, and...
  !
  ! the text of an actual TRANSP namelist (TR.DAT file).

  ! FORTRAN LUN

  integer :: lun = 82      ! TRANSP namelist file i/o LUN

  ! as with the old
  ! f77 splitn, this module and the routines which use it can convert
  ! a TRANSP namelist to the TR.ZDA form needed by f77 TRANSP COMMON-
  ! based namelist reader codes -- but will also support "editing 
  ! operations", and will contain a copy of the contents of the TRANSP
  ! namelist default values and actual values as encoded in the TR.DAT file.

  ! the text of an actual TRANSP namelist (TR.DAT file):

  integer :: curmax = 2000  ! current max no. of lines (expandable)
  integer :: nlines         ! no. of lines in TR.DAT file
  integer :: newline        ! where the next line goes...
  logical :: addblank       ! flag to insert a blank line

  character*150 aline,bline ! a 1-line buffer & a copy
  character*150, dimension(:), allocatable :: textnl  ! TR.DAT file text lines
  integer, dimension(:), allocatable :: lenl  ! abs() = length of each line
  integer, dimension(:), allocatable :: ordl  ! ordering of lines
  integer, dimension(:,:), allocatable :: namfld  ! location of name field
  integer, dimension(:), allocatable :: eqsfld    ! location of "=" sign
  integer, dimension(:,:), allocatable :: valfld  ! location of name field
  integer, dimension(:), allocatable :: cmtfld    ! location of cmt field

  ! this code would write the file out in its proper order:
  !  do i=1,nlines
  !    jl=ordl(i)
  !    if(lenl(jl).ne.0) then
  !      write(lun,'(A)') textnl(jl)(1:abs(lenl(jl)))
  !    endif
  !  enddo
  ! sign(lenl(...)) flags "special" namelist segments e.g. EFITIN
  ! if lenl(jl).eq.0, the line was marked as deleted.

  !-------------------------------------------------
  !  information to enable namelist editing via splitn
  !  cf splitn_put.f90 & splitn_edit_enable.f90

  character*30 :: edit_program = " "
  logical :: edit_enabled = .FALSE.    ! namelist edit lock
  logical :: edit_started = .FALSE.    ! flag if edit(s) have occurred

  !-------------------------------------------------
  ! namelist database (TRANSP & TRDAT combined)
  ! list of known namelist variables, their types, rank, dimensionality, etc.
  ! this is the nxlist.summary namelist, written by TRANSP code generators;
  ! therefore, its syntactic correctness is tested and assured.

  character*512 :: nlfile              ! name of namelist database file
  logical :: have_database = .FALSE.   ! initialization flag

  character*8 cur_naml      ! current namelist (in namelist database file)

  ! the following pertain to output (fortran readable) namelists:

  integer, parameter :: max_namls = 400 ! max no. of distinct namelists
  integer :: nnamls                     ! actual number
  character*8 all_namls(max_namls)      ! the list of namelist names...

  ! number of variables & updatable variables in each namelist
  integer, dimension(:), allocatable :: num_naml_vars
  integer, dimension(:), allocatable :: num_naml_vars_st

  ! index to variables & updatable variables in each namelist
  integer, dimension(:,:), allocatable :: indx_naml_vars
  integer, dimension(:,:), allocatable :: indx_naml_vars_st

  ! formerly in msplitn_par:
  ! now the same for both JET and PPPL versions:

  integer, parameter            :: nc_tri=3    &  !   # char, trigraph
                                  ,nc_pre=16   &  !           pre variable
                                  ,nc_ext=72      !           ext variable

  ! "maxrank" gives the maximum dimensionality or rank of any namelist
  ! array item.  ***Caution*** all of splitn should be checked if this
  ! ever wants to be increased!

  integer, parameter :: maxrank = 4

  type :: var

     character*32 :: name ! name of quantity
     character*5 :: type  ! type: R/I/L/D or C*n C*nn C*nnn
                          ! for real, integer, logical, real*8, character...
     character*8 :: naml  ! namelist to which quantity belongs

     integer :: steerable ! 0 means: not changeable on restart; 1 means:
                          ! changeable on TRDAT rerun; 2 means: changeable
                          ! on TRANSP restart.

     integer :: chsize    ! (char data only): size n of C*n strings.

     integer :: rank      ! dimensionality or rank {0,1,2,3, or 4}
     integer :: dims(2,maxrank) ! dims(1:2,1:rank) -- size of each dimensions
                          ! lower_limit:upper_limit

     character*20 short_dflt   ! default value string (if short)
     integer :: long_dflt_addr ! ptr to default value string (if long)

     integer :: addr      ! ptr to data (as per type)
     integer :: nlinadr   ! ptrs to TR.DAT file lines containing var. elements

     integer :: addr_st   ! update (steering block) offset data address
     integer :: nlinadr_st  ! ptrs to TR.DAT file lines containing update data

  end type var

  !--------------------------------------

  integer, parameter :: len_short=20  ! len(short_dflt) elements

  integer                                :: nvars      ! size of variables list
  type (var), dimension(:), allocatable  :: varlist
  integer, dimension(:), allocatable     :: var_order  ! alphabetic ordering

  ! list of deleted names

  integer :: ndnames  ! number of known deleted names
  character*32, dimension(:), allocatable :: dnames
  ! (names are put on this list to have splitn automatically comment them
  ! out instead of reporting an error, when encountering an undefined name).

  ! data buffers:  current namelist, namelist defaults (_d), write buffer (_w),
  !                previous namelist (_p)

  !--------------------------------------
  ! MOD dmc Mar 2009: namelists are now updatable; update sections can
  ! be added.  The number of update sections can vary from namelist to 
  ! namelist.  The default number of updates is zero.  The n<xxx>_st
  ! variables give the data sizes per data type per update section.

  ! #updates in data records:
  integer :: nupdate,nupdate_w,nupdate_p
  integer, parameter :: nupdate_d = 0

  ! to ease data mgmt, assume an interim maximum number of updates-- 
  ! this is expandable, but allows all data buffer sizes to be kept the
  ! same: for logbuf, nlog + nupdate_max*nlog_st; 
  !       for intbuf, nint + nupdate_max*nint_st, etc.

  integer :: nupdate_max = 2

  integer :: kupdate ! current update block (used by add_nl_line)

  ! updates occur in time sequences; the update times are stored here:
  real*8, dimension(:), allocatable :: tup,tup_d,tup_w,tup_p

  ! updates should occur after TINIT:
  real*8 :: tinit,tinit_d,tinit_w,tinit_p

  real*8 :: tlarge = 654.3d210  ! a very large time...
  ! tlarge is used to indicate absence of actual update blocks in namelists.

  integer :: nxblock  ! index of next update block 

  ! update block locations
  integer, dimension(:), allocatable :: linup

  ! provide upper limit on namelist file size (#lines)
  integer, parameter :: ilarge = 1000000000

  !--------------------------------------

  integer nlog,nlog_st
  logical, dimension(:), allocatable :: logbuf,logbuf_d,logbuf_w,logbuf_p

  integer nint,nint_st
  integer, dimension(:), allocatable :: intbuf,intbuf_d,intbuf_w,intbuf_p

  integer nreal,nreal_st
  real, dimension(:), allocatable :: rbuf,rbuf_d,rbuf_w,rbuf_p

  integer nr8,nr8_st
  real*8, dimension(:), allocatable :: dbuf,dbuf_d,dbuf_w,dbuf_p

  ! the length of a chbuf element is the presumed maximum length of any
  ! CHARACTER namelist variable (should match ch_maxlen):

  integer nchv,nchv_st
  character*128, dimension(:), allocatable :: chbuf,chbuf_d,chbuf_w,chbuf_p, &
       chtmp
  integer, parameter :: ch_maxlen=512
  character*128 :: chval

  !--------------------------------------
  ! indices of TR.DAT file lines which define var data elements
  integer, dimension(:), allocatable :: ilines
  !   note: this is the file position, not the storage position;
  !   text for element j at ordl(ilines(j)).

  ! this array allows for deferred assignmetns to ilines:
  integer :: iline_nstack
  integer, dimension(:,:), allocatable :: iline_stack

  ! and the range of characters within the line giving the element values.
  integer, dimension(:,:), allocatable :: ivrange

  ! associated repeat counts (0 if none)
  integer, dimension(:), allocatable :: krepeat

  ! range of characters specifying repeat count "nnn*"
  integer, dimension(:,:), allocatable :: irrange

  ! long (>len_short char) default value strings
  integer nlong
  character*512, dimension(:), allocatable :: long_dflts

  logical :: previous   ! T if there is a "previous" namelist

  !-----------------------------------------------------------------
  ! parsing-related information

  logical efitin_flag   ! .TRUE. iff reading EFITIN segment of namelist file.
  integer :: nam1,nam2  ! name field limits
  integer :: val1,val2  ! value field limits
  integer :: lcmt,len1  ! comment field start; length of line

  logical merge_flag    ! .TRUE. if namelist file merge operation was done;
  !  if merge_flag is set, only subsequent allowed action is to WRITE...

  integer, parameter :: maxwarn = 100
  character*32 dblasg1(maxwarn)  ! items doubly assigned (no value change)
  character*32 dblasg2(maxwarn)  ! items doubly assigned (value does change)
  character*32 lackdp(maxwarn)   ! floating point item w/o decimal pt
  integer :: kdblasg1(maxwarn)   ! block w/dblasg1 warning
  integer :: kdblasg2(maxwarn)   ! block w/dblasg2 warning
  integer :: klackdp(maxwarn)    ! block w/lackdp warning
  integer :: ndblasg1,ndblasg2,nlackdp  ! sizes of above lists.

  contains
    subroutine init_nl(clear_warn)
      logical, intent(in) :: clear_warn   ! .true. to clear the warning lists

      !  initialize memory-- to receive text image of TRANSP file

      if(allocated(textnl)) then
         deallocate(textnl,lenl,ordl,cmtfld,namfld,eqsfld,valfld)
      endif

      nlines=0

      !  initially assume no "update blocks" in namelist.
      kupdate = 0
      nupdate = 0
      nxblock = 1

      if(allocated(linup)) deallocate(linup)
      allocate(linup(nupdate_max)); linup=ilarge

      if(allocated(tup)) deallocate(tup)
      allocate(tup(nupdate_max)); tup = tlarge

      allocate(textnl(curmax)); textnl=' '
      allocate(lenl(curmax)); lenl=0
      allocate(ordl(curmax)); ordl=0
      allocate(cmtfld(curmax)); cmtfld=0
      allocate(eqsfld(curmax)); eqsfld=0
      allocate(namfld(2,curmax)); namfld=0
      allocate(valfld(2,curmax)); valfld=0

      !  clear warning lists
      if (clear_warn) then
         ndblasg1=0
         ndblasg2=0
         nlackdp=0
         kdblasg1=0
         kdblasg2=0
         klackdp=0
         dblasg1=' '
         dblasg2=' '
         lackdp=' '
      end if

      efitin_flag=.FALSE.
      merge_flag =.FALSE.

    end subroutine init_nl

    subroutine read_nl(filename,ios, quiet)

      ! read a TR.DAT file

      character*(*),intent(in) :: filename  ! path to namelist (TR.DAT) file
      integer, intent(out) :: ios           ! i/o status on open; 0=OK

      logical, intent(in), optional :: quiet

      !----------------------
      integer ieof,iblock
      logical :: iquiet
      !----------------------
      ! first make sure that namelist database
      ! containing namelist variable definitions and default values
      ! has been read in; re-initialize stored namelist values to 
      ! their defaults.

      call read_db(ios)
      if(ios.ne.0) return

      !----------------------
      ! open namelist file

      iquiet=.FALSE.
      if(present(quiet)) iquiet=quiet

      open(unit=lun,file=filename,status='old',iostat=ios)
      if(ios.ne.0) then
         if(.not.iquiet) then
            write(6,*) ' splitn_module(read_nl): open failure, iostat=',ios
            write(6,*) '   filename: ',trim(filename)
         endif
         return
      endif

      !----------------------
      ! read file...

      call init_nl(.true.)

      efitin_flag=.FALSE.

      do
         ! since (A) reads anything, assume iostat is set only for EOF

         read(lun,'(A)',iostat=ieof) aline
         if(ieof.ne.0) exit
        
         iblock=nxblock

         newline= nlines+1

         call parse_line(aline,.TRUE.,ios)
         if(ios.ne.0) exit  ! this is an actual error

         call add_nl_line(aline)

         if(iblock.lt.nxblock) then
            ! an update block was detected -- mark line position
            linup(iblock)=nlines  ! =newline
         endif

      enddo

      close(unit=lun)

      if(efitin_flag) then
         write(6,*) ' %splitn(read_nl): EOF during read of EFITIN namelist', &
              ' (corrected)'
         call clear_fields
         lcmt=1
         newline = nlines+1
         call add_nl_line('/')
         efitin_flag = .FALSE.
      endif

      if(ios.eq.0) then
         call splitn_check_misc(ios)
      endif

      !  save copy of namelist data in write buffers

      logbuf_w=logbuf
      intbuf_w=intbuf
      rbuf_w=rbuf
      dbuf_w=dbuf
      chbuf_w=chbuf

      tup_w=tup
      nupdate_w=nupdate
      tinit_w=tinit

    end subroutine read_nl

    subroutine splitn_check_misc(ierr)

      integer, intent(out) :: ierr

      ! check for double assignments, etc.; can be considered an error
      ! if environment variable is set.

      !--------------------------------------
      integer :: iwarn
      character*20 :: loc_testval
      !--------------------------------------

      ierr = 0

      if(ndblasg1.gt.0) call splitn_printwarn( &
           'namelist elements assigned the same value more than once:', &
           ndblasg1,len(dblasg1(1)),dblasg1,kdblasg1)
      if(ndblasg2.gt.0) call splitn_printwarn( &
           'namelist elements assigned more than once WITH change of value:', &
           ndblasg2,len(dblasg2(1)),dblasg2,kdblasg2)
      if(nlackdp.gt.0) call splitn_printwarn( &
           'namelist element value field(s): decimal point(s) inserted:', &
           nlackdp,len(lackdp(1)),lackdp,klackdp)

      !  DMC: option to make change of value an ERROR instead of a WARNING:
      !  (added Feb 2009):

      call mpi_sget_env('SPLITN_NO_REASSIGN',loc_testval,iwarn)
      call uupper(loc_testval)
      if((loc_testval.eq.'TRUE').AND.(ndblasg2.gt.0)) then
         ierr=999
         write(6,*) ' '
         write(6,*) ' -> Environment variable SPLITN_NO_REASSIGN = "TRUE",'
         write(6,*) ' therefore element value reassignment is deemed an error.'
         write(6,*) ' '
      endif

    end subroutine splitn_check_misc

    subroutine merge_read_nl(filename,tagname,ios)

      ! read a namelist fragment-- append to end; remove earlier copy
      ! of fragment as determined by tag and/or definition of variables

      character*(*),intent(in) :: filename  ! path to namelist (TR.DAT) file
      character*(*),intent(in) :: tagname   ! fragment tag name
      integer, intent(out) :: ios           ! i/o status on open; 0=OK

      !----------------------
      integer ieof,insave,inames
      character*18 fulltag
      character*32, dimension(:), allocatable :: znames
      character*80 nbuf
      character*9 zdate
      character*40 zuser
      integer :: iparen,ilen,i,j,k,ifind
      !----------------------

      if(merge_flag) then
         write(6,*) ' ?splitn_module: namelist already merged; next operation must be a write.'
         ios=1
         return
      endif

      if((nlines.le.0).or.(.not.have_database)) then
         write(6,*) ' ?splitn_module: normal namelist read must be done first'
         write(6,*) '  before attempting merge with 2nd namelist fragment.'
         ios=1
         return
      endif

      kupdate = 0  ! deal with main block only (no update blocks)

      !----------------------
      fulltag = '!$'//tagname
      call ulower(fulltag)
      !----------------------
      ! open namelist file

      open(unit=lun,file=filename,status='old',iostat=ios)
      if(ios.ne.0) then
         write(6,*) ' splitn_module(merge_read_nl): open failure, iostat=',ios
         write(6,*) '   filename: ',trim(filename)
         return
      endif

      !----------------------
      ! read file...

      insave = nlines
      inames = 0

      allocate(znames(nreal+nint+nlog+nr8+nchv))

      call c9date(zdate)
      call getlog(zuser)

      call clear_fields
      lcmt=3

      call get_newline
      call add_nl_line('  !------------------------------------------------'//&
           ' '//trim(fulltag))

      bline = '  ! '//trim(fulltag)//' inserted on '//zdate//' by user '//zuser
      call get_newline
      call add_nl_line(bline)

      do
         ! since (A) reads anything, assume iostat is set only for EOF

         read(lun,'(A)',iostat=ieof) aline
         if(ieof.ne.0) exit

         call get_newline

         call addtag

         call parse_line(aline,.FALSE.,ios)
         if(ios.ne.0) exit  ! this is an actual error

         if(nam1.gt.0) call addname

         call add_nl_line(aline)

      enddo

      close(unit=lun)

      if(efitin_flag) then
         write(6,*) ' %splitn(read_nl): EOF during read of EFITIN namelist', &
              ' (corrected)'
         call clear_fields
         lcmt=1
         call get_newline
         call add_nl_line('/')
         efitin_flag = .FALSE.
      endif

      !  mark pre-merge lines for deletion if either tag is present or if RHS
      !  variable is on list...

      do i=1,insave
         if(chktag(textnl(i))) then
            lenl(i)=0
         else
            if(namfld(1,i).eq.0) cycle
            ilen=namfld(2,i)-namfld(1,i)+1
            nbuf=textnl(i)(namfld(1,i):namfld(2,i))
            iparen=index(nbuf,'(')
            if(iparen.gt.0) then
               nbuf(iparen:ilen)=' '
            endif
            ilen=len(trim(nbuf))
            call ulower(nbuf(1:ilen))

            ifind=0
            do j=1,inames
               if(nbuf(1:ilen).eq.znames(j)) then
                  ifind=j
                  exit
               endif
            enddo

            if(ifind.gt.0) then
               lenl(i)=0
            endif
         endif
      enddo

      deallocate(znames)
      merge_flag = .TRUE.

      contains
        subroutine addtag

          !  add tag comment to inserted line...

          integer :: ilen,ilenf,insrt

          character*1 :: squot = "'"
          character*1 :: dquot = '"'
          character*1 :: cexcl = '!'
          character*1 :: zquot

          integer :: ic,icmt

          !----------------------------
          ! first look to see if comment is already present.

          if(chktag(aline)) return

          ilen=len(trim(aline))
          ilenf=len(fulltag)

          icmt = 0
          zquot = ' '
          do ic=1,ilen
             if(zquot.eq.' ') then
                if(aline(ic:ic).eq.cexcl) then
                   icmt=ic
                   exit
                endif
                if(aline(ic:ic).eq.squot) zquot=squot
                if(aline(ic:ic).eq.dquot) zquot=dquot
             else
                if(aline(ic:ic).eq.zquot) zquot=' '
             endif
          enddo

          if(icmt.gt.0) then
             !  comment field present -- append to end
             insrt = min((len(aline)-ilenf+1),ilen+2)
          else
             !  no coment field found
             insrt=max(65,min((len(aline)-ilenf+1),ilen+2))
          endif

          aline(insrt:insrt+ilenf-1)=fulltag

        end subroutine addtag

        logical function chktag(zline)

          character*(*), intent(in) :: zline

          bline=zline
          call ulower(bline)
          chktag = index(bline,fulltag).gt.0

        end function chktag

        subroutine addname

          nbuf = aline(nam1:nam2)
          iparen = index(nbuf,'(')
          if(iparen.gt.0) then
             nbuf(iparen:len(nbuf))=' '
          endif

          ilen=len(trim(nbuf))

          call ulower(nbuf(1:ilen))

          ifind=0
          do i=1,inames
             if(nbuf(1:ilen).eq.znames(i)) then
                ifind=i
                exit
             endif
          enddo

          if(ifind.eq.0) then
             inames=inames+1
             znames(inames)=nbuf(1:ilen)
          endif
        end subroutine addname

    end subroutine merge_read_nl

    subroutine get_newline
      ! get new line location before next update block

      addblank = .FALSE.

      if(kupdate.lt.nupdate) then
         newline = linup(kupdate+1)
         if(isblank(newline-1)) then
            newline=newline-1
         else
            addblank=.TRUE.
         endif
      else
         newline = nlines + 1
      endif

    end subroutine get_newline

    !-------------------------------------------------------------------
    subroutine parse_line(aline,ublock,ios)

      ! parse a single namelist input line
      ! <name>[(<indices>,...)] = [<values>,...] [! <comment>]
      ! LHS                       RHS             comment
      !
      ! only one LHS per line.  RHS can have one or more values but
      ! the number should be appropriate as per the storage associated
      ! with the variable specified on the LHS.  The value field can be
      ! null, in which case no assignment is made and the default value
      ! remains in effect.  commas "," separate individual values, and
      ! repeat counts can be used, e.g.
      !
      ! generally, aline is not modified.  EXCEPTION: decimal points will
      ! be inserted, with warning, in non-null floating point value fields
      ! that lack them.  EXCEPTION2: non-printable characters other than
      ! "tab" removed!
      !
      ! dmc Jan 2003 ...

      character*(*), intent(inout) :: aline  ! line from namelist file
      logical, intent(in) :: ublock          ! .TRUE. to allow update block
      integer, intent(out) :: ios            ! parse status, 0=OK

      !-------------------
      integer i,j,k,icmt,ic1,ic2,iwarncc
      integer ic1a,ic2a,ieqs,ierr,ivals,iblock,iadr
      integer imatch,jj,indx(maxrank),ipass
      character*150 wk
      character*1 zquot
      character*26 zalphab(2)
      real*8 :: uptime

      data zalphab/ &
           'abcdefghijklmnopqrstuvwxyz', &
           'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      data iwarncc/-3/

      !-------------------

      ipass=0

5     continue  ! re-enter here, if line is commented out on the fly

      ios=0

      ipass = ipass + 1
      if(ipass.gt.2) then
         ios = 99
         write(6,*) ' ?splitn_module(parse_line) ipass>2: internal error.'
         return
      endif

      if(merge_flag) then
         write(6,*) ' ?splitn_module: line parse unavailable after merge.'
         ios=1
         return
      endif

      call clear_fields

      !  get non-blank length of "aline" -- remove any non-printable characters

      len1=0
      do i=len(aline),1,-1
         j=ichar(aline(i:i))
         if(j.ne.9) then
            if((j.lt.32).or.(j.eq.127)) then
               if(iwarncc.lt.0) then
                  iwarncc=iwarncc+1
                  write(6,*) &
                       ' %splitn_module:parse: non-printable character: (', &
                       j,') replaced with blank.'
               endif
               aline(i:i)=' '
               j=32
            endif
         endif
         if(len1.eq.0) then
            if((j.ne.9).and.(j.ne.32)) len1=i
         endif
      enddo

      wk=aline       ! local, writable copy of line.

      zquot=' '      ! use for quote balance checking
      icmt=0         ! location of "!", start of comment
      ic1=0          ! start of non-blank non-"!"-comment section
      ic1a=0         ! last non-blank before "=" sign
      ic2a=0         ! first non-blank after "=" sign
      ic2=0          ! end of non-blank non-"!"-comment section
      ieqs=0         ! first "=" sign in non-comment non-quoted section
      
      i=0

      ! scan the line:  look for comment fields and non-comment section.
      ! Do quote balance checking.  Convert non-comment non-quoted sections
      ! to uppercase.  Detab.

      do
         i=i+1
         if(i.gt.len1) exit

         if(zquot.ne.' ') then
            ! inside quote string: look for matching, closing quote ONLY
            if(wk(i:i).eq.zquot) then
               zquot=' '
               if(ieqs.ne.0) ic2=i
            endif
            cycle

         else if((wk(i:i).eq."'").or.(wk(i:i).eq.'"')) then
            ! start of quote string
            zquot=wk(i:i)
            if((ieqs.ne.0).and.(ic2a.eq.0)) ic2a=i
            cycle

         endif

         ! outside quote string, before "!" comment:

         if(wk(i:i).eq.char(9)) wk(i:i)=' '  ! detab
         if(wk(i:i).eq.' ') cycle

         ! non-blank

         if((wk(i:i).eq.'!').or.(wk(i:i).eq.'$').or. &
              (wk(i:i).eq.'&').or.(wk(i:i).eq.'/')) then
            icmt=i           ! comment field found; exit
            exit
         endif

         ! bounds of non-blank non-comment section

         if(ic1.eq.0) ic1=i
         ic2=i

         if(wk(i:i).eq.'=') then
            if(ieqs.eq.0) ieqs=i
            cycle
         endif

         if(ieqs.eq.0) ic1a=i
         if((ieqs.ne.0).and.(ic2a.eq.0)) ic2a=i

         ! uppercase conversion

         j=index(zalphab(1),wk(i:i))
         if(j.gt.0) wk(i:i)=zalphab(2)(j:j)
      enddo

      ! check for unbalanced quotes:
      if(zquot.ne.' ') then
         write(6,*) '?splitn_module(parse_line): unbalanced quotes detected:'
         write(6,*) ' ',aline(1:len1)
         ios=1
         return
      endif

      lcmt=icmt

      ! OK: check for &<namelist-name> &END or similar construct

      if(icmt.gt.0) then
         if(wk(icmt:icmt).eq.'/') then
            if (ieqs>0 .and. ieqs<icmt) then
               write(6,*) '?splitn_module(parse_line): end of namelist character "/" not allowed with variable definition:'
               write(6,*) ' ',aline(1:len1)
               ios=1
            else
               efitin_flag=.FALSE.  ! end-of-namelist indicator
            end if
            return
         endif

         if((wk(icmt:icmt).eq.'$').or.(wk(icmt:icmt).eq.'&')) then
            if(wk(icmt+1:icmt+4).eq.'END ') then
               efitin_flag=.FALSE.  ! end-of-namelist indicator
            else if(wk(icmt+1:icmt+7).eq.'EFITIN ') then
               if(kupdate.gt.0) then
                  write(6,*) '?splitn_module(parse_line): EFITIN namelist:'
                  write(6,*) ' cannot be positioned after start of update blocks.'
                  ios=1
               else
                  efitin_flag=.TRUE.
               endif
            endif
            return
         endif
      endif

      if(ic1.eq.0) return  ! line is comment-only

      ! it is more than just a comment, so...
      ! OK: must have name = value combination or name = <NULL> combination.

      if(.not.efitin_flag) then
         if((ieqs.eq.0).or.(ic1.ge.ieqs)) then
            write(6,*) &
                 '?splitn_module(parse_line): expected <LHS> = <RHS> form.'
            if(ieqs.eq.0) then
               write(6,*) ' "=" sign not found.'
            else if(ic1.ge.ieqs) then
               write(6,*) ' LHS (left hand side) not found.'
            endif
            write(6,*) ' ',aline(1:len1)
            ios=1
            return
         endif

         nam1=ic1
         nam2=ic1a
         val1=ic2a
         val2=ic2
         if(ic2.le.ieqs) then
            val1=ieqs+1
            val2=ieqs+1   ! null RHS (blank)
            wk(val1:val2)=' '
         else if((ic2a.eq.ic2).and.(wk(ic2:ic2).eq.',')) then
            val1=ieqs+1
            val2=ieqs+1   ! null RHS (",")
            wk(val1:val2)=' '
         endif

         ! parse LHS -- check that name is known & indexing is OK

         call splitn_lhs(wk(nam1:nam2),imatch,indx,maxrank,ierr)
         if(ierr.ne.0) then
            write(6,*) ' ?splitn_module (parse_line): LHS parse error:'
            write(6,*) ' ',aline(1:len1)
            ios=1
            return
         endif

         if(imatch.eq.0) then
            !  ierr=0 && imatch=0 means detected a known deleted namelist name 
            bline = aline
            aline = '  !(splitn_deleted)  '//trim(bline)
            go to 5  ! reprocess line as commented out
         endif

         jj=var_order(imatch)

         ! check value fields of all floating point types for decimal pts.

         if((varlist(jj)%type.eq.'R').or.(varlist(jj)%type.eq.'D')) then
            call splitn_chkdp(aline,wk)
         endif

         if((varlist(jj)%steerable.eq.3).and.(.NOT.ublock)) then
            write(6,*) ' ?splitn_module (parse_line): update block not supported in merge list:'
            write(6,*) ' ',aline(1:len1)
            ios=1
            return
         endif

         ! parse RHS

         if(varlist(jj)%rank.eq.0) then

            if(varlist(jj)%steerable.eq.3) then
               !  may want to make room for more blocks, before RHS parse...
               if(nxblock.ge.nupdate_max) call spxb_expand
            endif

            ! scalar

            call splitn_rhs0(aline(1:len1),wk(val1:val2),ivals,jj)
            if(ivals.ne.1) then
               write(6,*) ' ?splitn_module (parse line): error in line:'
               if(ivals.gt.1) write(6,*) '  too many values for scalar:'
               write(6,*) ' ',aline(1:len1)
               ios=1
               return
            endif

            if(varlist(jj)%steerable.eq.3) then
               iblock=nxblock-1
               if(iblock.eq.0) then
                  iadr = varlist(jj)%addr
               else
                  iadr = nr8 + (iblock-1)*nr8_st + varlist(jj)%addr_st
               endif

               uptime = dbuf(iadr)
               if(iblock.gt.0) then
                  if(uptime.le.tup(iblock)) then
                     write(6,*) ' ?splitn_module (parse line): update block time not in ascending order:'
                     write(6,*) ' ',aline(1:len1)
                     write(6,*) '  prior update time: ',tup(iblock)
                     ios=1
                     return
                  endif
               else
                  if(uptime.le.tinit) then
                     write(6,*) ' ?splitn_module (parse line): first update block time <= TINIT:'
                     write(6,*) ' ',aline(1:len1)
                     write(6,*) '  TINIT in namelist: ',tinit
                     ios=1
                     return
                  endif
               endif

               write(6,*) ' %splitn_module: update block detected, t=',uptime

               tup(nxblock) = uptime
               nupdate = nxblock
               kupdate = nupdate
               nxblock = nxblock + 1
            endif

            if(varlist(jj)%name.eq.'TINIT') then
               iadr = varlist(jj)%addr
               tinit = dbuf(iadr)
            endif
                     
         else

            call splitn_rhs(aline(1:len1),wk(val1:val2),jj,indx,maxrank,ios)

         endif

      else

         ! inside EFITIN namelist...

         if(ic1.gt.0) then
            if(wk(ic1:ic1).eq.'~') then
               write(6,*) ' ?splitn_module (parse_line): error in line:'
               write(6,*) '  "~" variable in EFITIN namelist.'
               write(6,*) ' ',aline(1:len1)
               ios=1
            endif
         endif
         
      endif  ! efitin_flag

    end subroutine parse_line

    !-------------------------------------------------------------------
    subroutine spxb_expand
      ! increase nupdate_max; re-allocate and expand arrays

      integer, parameter :: incr=10  ! add 10 slots per call
      integer :: inumpre             ! previous number
      integer :: isizepre,isizenew   ! for buffers: old & new sizes

      integer :: ii,jj

      real, dimension(:), allocatable :: r4tmp
      real*8, dimension(:), allocatable :: r8tmp
      integer, dimension(:), allocatable :: itmp
      logical, dimension(:), allocatable :: ltmp
      character*128, dimension(:), allocatable :: chtmp
      integer, dimension(:,:), allocatable :: itmp2

      !------------------------------

      inumpre = nupdate_max
      nupdate_max = nupdate_max + incr

      !------------------------------
      !  update times & update marker line counts

      allocate(r8tmp(inumpre))
      allocate(itmp(inumpre))

      r8tmp = tup
      deallocate(tup); allocate(tup(nupdate_max))
      tup(1:inumpre)=r8tmp
      tup(inumpre+1:)=tlarge

      r8tmp = tup_p
      deallocate(tup_p); allocate(tup_p(nupdate_max))
      tup_p(1:inumpre)=r8tmp
      tup_p(inumpre+1:)=tlarge

      r8tmp = tup_d
      deallocate(tup_d); allocate(tup_d(nupdate_max))
      tup_d(1:inumpre)=r8tmp
      tup_d(inumpre+1:)=tlarge

      r8tmp = tup_w
      deallocate(tup_w); allocate(tup_w(nupdate_max))
      tup_w(1:inumpre)=r8tmp
      tup_w(inumpre+1:)=tlarge

      itmp = linup
      deallocate(linup); allocate(linup(nupdate_max))
      linup(1:inumpre)=itmp
      linup(inumpre+1:)=ilarge

      deallocate(r8tmp,itmp)

      !------------------------------
      !  update data buffers...

      if(nlog_st.gt.0) then
         ! logical...

         isizepre = nlog + nlog_st*inumpre
         isizenew = nlog + nlog_st*nupdate_max

         allocate(ltmp(isizepre))

         ltmp=logbuf
         deallocate(logbuf); allocate(logbuf(isizenew))
         logbuf(1:isizepre)=ltmp
         call lset_ucop(logbuf,1,nlog_st,1,0,inumpre,nupdate_max)

         ltmp=logbuf_p
         deallocate(logbuf_p); allocate(logbuf_p(isizenew))
         logbuf_p(1:isizepre)=ltmp
         call lset_ucop(logbuf_p,1,nlog_st,1,0,inumpre,nupdate_max)

         ltmp=logbuf_d
         deallocate(logbuf_d); allocate(logbuf_d(isizenew))
         logbuf_d(1:isizepre)=ltmp
         call lset_ucop(logbuf_d,1,nlog_st,1,0,inumpre,nupdate_max)

         ltmp=logbuf_w
         deallocate(logbuf_w); allocate(logbuf_w(isizenew))
         logbuf_w(1:isizepre)=ltmp
         call lset_ucop(logbuf_w,1,nlog_st,1,0,inumpre,nupdate_max)

         deallocate(ltmp)
      endif

      if(nint_st.gt.0) then
         ! integer...

         isizepre = nint + nint_st*inumpre
         isizenew = nint + nint_st*nupdate_max

         allocate(itmp(isizepre))

         itmp=intbuf
         deallocate(intbuf); allocate(intbuf(isizenew))
         intbuf(1:isizepre)=itmp
         call iset_ucop(intbuf,1,nint_st,1,0,inumpre,nupdate_max)

         itmp=intbuf_p
         deallocate(intbuf_p); allocate(intbuf_p(isizenew))
         intbuf_p(1:isizepre)=itmp
         call iset_ucop(intbuf_p,1,nint_st,1,0,inumpre,nupdate_max)

         itmp=intbuf_d
         deallocate(intbuf_d); allocate(intbuf_d(isizenew))
         intbuf_d(1:isizepre)=itmp
         call iset_ucop(intbuf_d,1,nint_st,1,0,inumpre,nupdate_max)

         itmp=intbuf_w
         deallocate(intbuf_w); allocate(intbuf_w(isizenew))
         intbuf_w(1:isizepre)=itmp
         call iset_ucop(intbuf_w,1,nint_st,1,0,inumpre,nupdate_max)

         deallocate(itmp)
      endif

      if(nreal_st.gt.0) then
         ! real*4...

         isizepre = nreal + nreal_st*inumpre
         isizenew = nreal + nreal_st*nupdate_max

         allocate(r4tmp(isizepre))

         r4tmp=rbuf
         deallocate(rbuf); allocate(rbuf(isizenew))
         rbuf(1:isizepre)=r4tmp
         call rset_ucop(rbuf,1,nreal_st,1,0,inumpre,nupdate_max)

         r4tmp=rbuf_p
         deallocate(rbuf_p); allocate(rbuf_p(isizenew))
         rbuf_p(1:isizepre)=r4tmp
         call rset_ucop(rbuf_p,1,nreal_st,1,0,inumpre,nupdate_max)

         r4tmp=rbuf_d
         deallocate(rbuf_d); allocate(rbuf_d(isizenew))
         rbuf_d(1:isizepre)=r4tmp
         call rset_ucop(rbuf_d,1,nreal_st,1,0,inumpre,nupdate_max)

         r4tmp=rbuf_w
         deallocate(rbuf_w); allocate(rbuf_w(isizenew))
         rbuf_w(1:isizepre)=r4tmp
         call rset_ucop(rbuf_w,1,nreal_st,1,0,inumpre,nupdate_max)

         deallocate(r4tmp)
      endif

      if(nr8_st.gt.0) then
         ! real*8...

         isizepre = nr8 + nr8_st*inumpre
         isizenew = nr8 + nr8_st*nupdate_max

         allocate(r8tmp(isizepre))

         r8tmp=dbuf
         deallocate(dbuf); allocate(dbuf(isizenew))
         dbuf(1:isizepre)=r8tmp
         call dset_ucop(dbuf,1,nr8_st,1,0,inumpre,nupdate_max)

         r8tmp=dbuf_p
         deallocate(dbuf_p); allocate(dbuf_p(isizenew))
         dbuf_p(1:isizepre)=r8tmp
         call dset_ucop(dbuf_p,1,nr8_st,1,0,inumpre,nupdate_max)

         r8tmp=dbuf_d
         deallocate(dbuf_d); allocate(dbuf_d(isizenew))
         dbuf_d(1:isizepre)=r8tmp
         call dset_ucop(dbuf_d,1,nr8_st,1,0,inumpre,nupdate_max)

         r8tmp=dbuf_w
         deallocate(dbuf_w); allocate(dbuf_w(isizenew))
         dbuf_w(1:isizepre)=r8tmp
         call dset_ucop(dbuf_w,1,nr8_st,1,0,inumpre,nupdate_max)

         deallocate(r8tmp)
      endif

      if(nchv_st.gt.0) then
         ! character*128...

         isizepre = nchv + nchv_st*inumpre
         isizenew = nchv + nchv_st*nupdate_max

         allocate(chtmp(isizepre))

         chtmp=chbuf
         deallocate(chbuf); allocate(chbuf(isizenew))
         chbuf(1:isizepre)=chtmp
         call chset_ucop(chbuf,1,nchv_st,1,0,inumpre,nupdate_max)

         chtmp=chbuf_p
         deallocate(chbuf_p); allocate(chbuf_p(isizenew))
         chbuf_p(1:isizepre)=chtmp
         call chset_ucop(chbuf_p,1,nchv_st,1,0,inumpre,nupdate_max)

         chtmp=chbuf_d
         deallocate(chbuf_d); allocate(chbuf_d(isizenew))
         chbuf_d(1:isizepre)=chtmp
         call chset_ucop(chbuf_d,1,nchv_st,1,0,inumpre,nupdate_max)

         chtmp=chbuf_w
         deallocate(chbuf_w); allocate(chbuf_w(isizenew))
         chbuf_w(1:isizepre)=chtmp
         call chset_ucop(chbuf_w,1,nchv_st,1,0,inumpre,nupdate_max)

         deallocate(chtmp)
      endif

      !------------------------------
      !  line tracking arrays...

      isizepre = nreal+nint+nlog+nr8+nchv + &
           inumpre*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st)

      isizenew = nreal+nint+nlog+nr8+nchv + &
           nupdate_max*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st)

      allocate(itmp(isizepre),itmp2(2,isizepre))

      deallocate(iline_stack); allocate(iline_stack(2,isizenew))

      itmp = ilines
      deallocate(ilines); allocate(ilines(isizenew))
      ilines(1:isizepre) = itmp
      ilines(isizepre+1:isizenew) = 0 

      itmp = krepeat
      deallocate(krepeat); allocate(krepeat(isizenew))
      krepeat(1:isizepre) = itmp
      krepeat(isizepre+1:isizenew) = 0 

      itmp2 = ivrange
      deallocate(ivrange); allocate(ivrange(2,isizenew))
      do jj=1,isizenew
         do ii=1,2
            if(jj.le.isizepre) then
               ivrange(ii,jj)=itmp2(ii,jj)
            else
               ivrange(ii,jj)=0
            endif
         enddo
      enddo

      itmp2 = irrange
      deallocate(irrange); allocate(irrange(2,isizenew))
      do jj=1,isizenew
         do ii=1,2
            if(jj.le.isizepre) then
               irrange(ii,jj)=itmp2(ii,jj)
            else
               irrange(ii,jj)=0
            endif
         enddo
      enddo

      deallocate(itmp,itmp2)

    end subroutine spxb_expand
    !-------------------------------------------------------------------
    subroutine clear_fields

      nam1=0; nam2=0 ! name field 
      val1=0; val2=0 ! value field
      lcmt=0         ! start of comment field

      iline_nstack=0 ! clear deferred ilines assignments

    end subroutine clear_fields

    !-------------------------------------------------------------------
    subroutine write_nl(filename,ios, quiet)

      ! write a namelist file

      character*(*),intent(in) :: filename  ! path to namelist (TR.DAT) file
      integer, intent(out) :: ios           ! i/o status on open; 0=OK

      logical, intent(in), optional :: quiet

      integer i,j,ilf
      logical :: iquiet

      iquiet = .FALSE.
      if(present(quiet)) iquiet=quiet

      open(unit=lun,file=filename,status='old',iostat=ios)
      if(ios.eq.0) then
         ! file exists: rename it
         close(unit=lun)
         ilf=len_trim(filename)
         call frename(filename(1:ilf),filename(1:ilf)//'~',ios)
      endif

      open(unit=lun,file=filename,status='new',iostat=ios)
      if(ios.ne.0) then
         if(.not.iquiet) then
            write(6,*) ' splitn_module(write_nl): open failure, iostat=',ios
            write(6,*) '   filename: ',trim(filename)
         endif
         return
      endif

      do i=1,nlines
         j=ordl(i)
         if(lenl(j).ne.0) then
            write(lun,'(A)') textnl(j)(1:abs(lenl(j)))
         endif
      enddo

      close(unit=lun)

    end subroutine write_nl

    logical function isblank(ilini)
      ! return TRUE if indicated line is all blank

      integer, intent(in) :: ilini
      integer :: indx

      if(ilini.eq.0) then
         isblank=.FALSE.
         return
      endif

      indx=ordl(ilini)

      isblank = (namfld(1,indx).eq.0).AND.(eqsfld(indx).eq.0).AND. &
           (valfld(1,indx).eq.0).AND.(cmtfld(indx).eq.0)

    end function isblank

    subroutine add_nl_line(zline)

      ! add a line to the memory image of the file
      ! increase the available memory, if necessary

      ! insert line at location (newline)

      ! NOTE: this routine expects not just the "zline" text is input
      !   but also that nam1:nam2 and val1:val2 and possibly lcmt are set...

      character*(*), intent(in) :: zline   ! line to add

      ! ipos = 0 means: add line at end of current update block (EOB)
      ! ipos = -1 means: add line at end of file (EOF)
      ! ipos > 0 menas: insert line at indicated location; line currently
      !                 at location and subsequent lines are shifted down one.

      integer, dimension(:), allocatable :: ilenl,iordl,icmtfld,ieqsfld
      integer, dimension(:,:), allocatable :: inamfld,ivalfld
      character*150, dimension(:), allocatable :: tmptext
      integer klines,i,j,iposi,iposj,ishift
      integer :: ii,ilin

      !------------------------------------

      iposi=newline
      ishift=1
      if(addblank) ishift=2

      if(nlines+ishift-1.ge.curmax) then
         klines=nlines
         !  will expand memory -- make temporary copies --
         allocate(iordl(klines),ilenl(klines),ieqsfld(klines),tmptext(klines))
         allocate(inamfld(2,klines),ivalfld(2,klines),icmtfld(klines))
         do i=1,klines
            j=ordl(i)    ! copies put in sorted order
            iordl(i)=i
            ilenl(i)=lenl(j)
            tmptext(i)=textnl(j)
            inamfld(1:2,i)=namfld(1:2,j)
            ieqsfld(i)=eqsfld(j)
            ivalfld(1:2,i)=valfld(1:2,j)
            icmtfld(i)=cmtfld(j)
         enddo

         curmax=2*curmax           ! double the allocated memory

         call init_nl(.false.)     ! fetch & initialize

         nlines=klines             ! copy back from temp arrays
         ordl(1:klines)=iordl
         lenl(1:klines)=ilenl
         textnl(1:klines)=tmptext
         namfld(1:2,1:klines)=inamfld
         eqsfld(1:klines)=ieqsfld
         valfld(1:2,1:klines)=ivalfld
         cmtfld(1:klines)=icmtfld

         deallocate(iordl,ilenl,tmptext,ieqsfld)  ! free temp arrays
         deallocate(inamfld,ivalfld,icmtfld)

      endif

      nlines = nlines+1
      textnl(nlines)=zline
      lenl(nlines)=max(1,len_trim(zline))
      if(efitin_flag) lenl(nlines)=-lenl(nlines)
      namfld(1,nlines)=nam1
      namfld(2,nlines)=nam2
      if(namfld(1,nlines).gt.0) then
         eqsfld(nlines)=index(zline,'=')
      else
         eqsfld(nlines)=0
      endif
      valfld(1,nlines)=val1
      valfld(2,nlines)=val2
      cmtfld(nlines)=lcmt

      if(addblank) then
         nlines = nlines + 1
         textnl(nlines)=' '
         lenl(nlines)=1
         if(efitin_flag) lenl(nlines)=-lenl(nlines)
         namfld(1:2,nlines)=0
         eqsfld(nlines)=0
         valfld(1:2,nlines)=0
         cmtfld(nlines)=0
      endif

      if((iposi.ge.nlines).or.(iposi.le.0)) then
         ordl(nlines)=nlines
         if(addblank) ordl(nlines-1)=nlines-1
      else
         do i=nlines,iposi+ishift,-1
            ordl(i)=ordl(i-ishift)
         enddo
         if(addblank) then
            ordl(iposi)=nlines-1
            ordl(iposi+1)=nlines
         else
            ordl(iposi)=nlines
         endif

         do i=1,maxlines(0)
            if(ilines(i).ge.iposi) ilines(i)=ilines(i)+ishift
         enddo

         !  check block bdys also
         do i=1,nupdate_max
            if((linup(i).ge.iposi).and.(linup(i).lt.ilarge)) then
               linup(i)=linup(i)+ishift
            endif
         enddo

      endif

      ! deferred ilines assignment now...
      do i=1,iline_nstack
         ii=iline_stack(1,i)
         ilin=iline_stack(2,i)
         ilines(ii)=ilin
      enddo
      iline_nstack = 0

      addblank = .FALSE.

    end subroutine add_nl_line

    integer function maxlines(idum)

      integer, intent(in) :: idum

      maxlines = nreal+nint+nlog+nr8+nchv
      maxlines = maxlines + &
           nupdate_max*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st)

    end function maxlines

    subroutine read_db(ios)

      ! open file and read namelist database;
      ! acquire default values & initialize actual values for real namelist
      ! file (about to be read).  If database has already been read, then,
      ! just re-initialize the actual values.

      ! namelist file is expected in "standard place", as per the
      ! code, below

      ! mod DMC -- "~UPDATE_TIME" is hardwired in

      integer, intent(out) :: ios           ! i/o status on open; 0=OK

      ! ---------------------------
      character*32 zname         ! var. name, this line
      character*5 ztype          ! var. type, this line
      integer :: ksteer          ! steerability flag
      integer iirank,iidims(2,4) ! rank (0 for scalar, 1 for vector, ...)
                                 ! & dimensions (match dim element of var type)
      integer ichsize            ! size (n) of C*n strings 
      character*512 zdflt        ! default value string

      character*84 iline         ! line buffer (actual file line length <= 80)
      integer i,j,ii,inaml
      integer isize,ivecsz

      integer ivars,ilong,ilog,iint,ichv,ir8,ireal,ilinp,ilf
      integer :: jlog,jint,jchv,jr8,jreal
      integer :: ilog_st,iint_st,ichv_st,ir8_st,ireal_st,ilinp_st
      integer :: ilen_naml,imaxlen_naml,ios2
      integer :: ibefore,imatch

      ! ---------------------------
      ios = 0
      if(have_database) then

         ! not the first namelist read... store old namelist in "_p" buffers
         logbuf_p = logbuf
         intbuf_p = intbuf
         rbuf_p = rbuf
         dbuf_p = dbuf
         chbuf_p = chbuf

         tup_p=tup
         nupdate_p=nupdate
         tinit_p=tinit

         previous = .TRUE.

         ! re-initialize to default state (before real namelist read)
         logbuf = logbuf_d
         intbuf = intbuf_d
         rbuf = rbuf_d
         dbuf = dbuf_d
         chbuf = chbuf_d

         tup=tup_d
         nupdate=nupdate_d  ! zero
         kupdate=0
         nxblock=1

         tinit=tinit_d

         ilines = 0
         ivrange = 0
         krepeat = 0
         irrange = 0

         return

      endif

      ! ---------------------------
      ! find namelist database file

      nlfile='nxlist.summary'
      ilf=len_trim(nlfile)

      open(unit=lun,file=nlfile(1:ilf),status='old', &
           action='READ',iostat=ios)

      if(ios.ne.0) then
         call ufilnam('CONFIGDIR','nxlist.summary',nlfile)
         ilf=len_trim(nlfile)

         open(unit=lun,file=nlfile(1:ilf),status='old', &
              action='READ',iostat=ios)
      endif

      if(ios.ne.0) then
         call ufilnam('NTCCHOME','bin/nxlist.summary',nlfile)
         ilf=len_trim(nlfile)
         open(unit=lun,file=nlfile(1:ilf), &
              status='old',iostat=ios)
      endif
      if(ios.ne.0) then
         call ufilnam('TRANSP_LOCATION','nxlist.summary',nlfile)
         ilf=len_trim(nlfile)
         open(unit=lun,file=nlfile(1:ilf), &
              status='old',iostat=ios)
      endif
      
      if(ios.ne.0) then
         write(6,*) ' ?splitn_module(read_db): cannot find "nxlist.summary".'
         write(6,*) '  (TRANSP namelist definition database not found).'
         return
      endif

      have_database = .TRUE.

      ! the database file is read in two passes.  In the first pass, 
      ! the sizes of lists and data structures are determined; then,
      ! after the necessary arrays have been allocated, the file is
      ! reread with the data being stored in the allocated structures.

      nvars=0
      nlong=0

      nlog=0
      nint=0
      nchv=0
      nr8=0
      nreal=0

      nlog_st=0
      nint_st=0
      nchv_st=0
      nr8_st=0
      nreal_st=0

      nnamls=0

      nupdate=0; nupdate_w=0; nupdate_p=0  ! parameter defines nupdate_d = 0 
      kupdate=0
      nxblock=1

      ilen_naml=0
      imaxlen_naml=0

      ! allocate and define initial update time vectors (dflt: no updates)
      allocate(tup(nupdate_max),tup_d(nupdate_max), &
           tup_p(nupdate_max),tup_w(nupdate_max))

      !-------------------
      !  add space for splitn-defined variables

      nvars = nvars + 1
      nr8 = nr8 + 1
      nr8_st = nr8_st + 1

      !-------------------
      !  find the requisite space

      do 
         read(lun,'(A)',iostat=ios) iline
         if(ios.ne.0) exit

         j=0
         do i=1,len(iline)
            if((iline(i:i).ne.' ').and.(iline.ne.char(9))) then
               j=i
               exit  ! non-whitespace found
            endif
         enddo
         if(j.eq.0) cycle   ! ignore all-blank line

         if(iline(j:j).eq.'#') cycle   ! ignore comment line

         if(iline(j:j).eq.'*') then
            cur_naml=iline(j+1:j+len(cur_naml))  ! set current namelist 
            cycle
         endif

         nvars=nvars+1
         call read1_db(iline,j,zname,ksteer,ztype,ichsize,iirank,iidims,zdflt)

         isize=igsize(iirank,iidims)
         if(ztype.eq.'R') then
            nreal=nreal+isize
            if(ksteer.eq.2) nreal_st=nreal_st+isize
         else if(ztype.eq.'I') then
            nint=nint+isize
            if(ksteer.eq.2) nint_st=nint_st+isize
         else if(ztype.eq.'L') then
            nlog=nlog+isize
            if(ksteer.eq.2) nlog_st=nlog_st+isize
         else if(ztype.eq.'D') then
            nr8=nr8+isize
            if(ksteer.eq.2) nr8_st=nr8_st+isize
         else if(ztype(1:1).eq.'C') then
            nchv=nchv+isize
            if(ksteer.eq.2) nchv_st=nchv_st+isize
         else
            call errmsg_exit( &
                 '?splitn_module: namelist database: uknown type: '//ztype)
         endif
         if(len_trim(zdflt).gt.len_short) then
            nlong=nlong+1
         endif
      enddo

      nreal=max(1,nreal)
      nr8=max(1,nr8)
      nlog=max(1,nlog)
      nint=max(1,nint)
      nchv=max(1,nchv)

      jreal=nreal + nupdate_max*nreal_st
      jr8=nr8 + nupdate_max*nr8_st
      jlog=nlog + nupdate_max*nlog_st
      jint=nint + nupdate_max*nint_st
      jchv=nchv + nupdate_max*nchv_st

      ! OK: have sizes, now allocate and initialize storage arrays

      allocate(intbuf(jint),logbuf(jlog),rbuf(jreal),dbuf(jr8))
      intbuf=0
      logbuf=.FALSE.
      rbuf=0
      dbuf=0
      allocate(chbuf(jchv),long_dflts(nlong)); chbuf=' '; long_dflts=' '
      allocate(varlist(nvars),var_order(nvars)); var_order=0
      allocate(ilines(jreal+jint+jlog+jr8+jchv)); ilines=0
      allocate(ivrange(2,jreal+jint+jlog+jr8+jchv)); ivrange=0
      allocate(krepeat(jreal+jint+jlog+jr8+jchv)); krepeat=0
      allocate(irrange(2,jreal+jint+jlog+jr8+jchv)); irrange=0
      allocate(dnames(nvars)); dnames=' '; ndnames=0

      allocate(iline_stack(2,jreal+jint+jlog+jr8+jchv))

      do i=1,nvars
         varlist(i)%name=' '
         varlist(i)%type=' '
         varlist(i)%naml=' '
         varlist(i)%steerable=0
         varlist(i)%chsize=0
         varlist(i)%rank=0
         varlist(i)%dims=0
         varlist(i)%short_dflt=' '
         varlist(i)%long_dflt_addr=0
         varlist(i)%addr=0
         varlist(i)%nlinadr=0
         varlist(i)%addr_st=0
         varlist(i)%nlinadr_st=0
      enddo

      ivars=0
      ilog=1   ! for data of various types...
      iint=1
      ichv=1
      ireal=1
      ir8=1

      ilong=1  ! for long default strings

      ilinp=1  ! for namelist file line ptrs

      ilog_st=1   ! for data of various types w/in update blocks
      iint_st=1
      ichv_st=1
      ireal_st=1
      ir8_st=1

      ilinp_st=1

      !-------------------
      ! insert splitn-defined variable

      ivars = ivars + 1
      zname = '~UPDATE_TIME'
      ksteer=3    ! special for (splitn)
      ztype='D'
      ichsize=0
      iirank=0
      iidims=0
      zdflt = '654.3d210'

      varlist(ivars)%name=zname
      call insert_var(ivars)   ! maintain alphabetic order
      varlist(ivars)%steerable = ksteer
      varlist(ivars)%type=ztype
      varlist(ivars)%chsize=ichsize
      varlist(ivars)%naml='(splitn)'
      varlist(ivars)%rank=iirank
      varlist(ivars)%dims=iidims
      varlist(ivars)%short_dflt=zdflt(1:len_short)

      isize=1
      ivecsz=1
      varlist(ivars)%addr=ir8
      varlist(ivars)%nlinadr=ilinp
      varlist(ivars)%addr_st=ir8_st
      varlist(ivars)%nlinadr_st=ilinp_st

      call dset_dflt(dbuf(ir8:ir8+isize-1),isize,ivecsz,ivars)
      call dset_ucop(dbuf,ir8,isize,ir8_st,0,0,nupdate_max)

      tlarge = dbuf(ir8)  ! make sure this matches decode

      ir8 = ir8 + isize
      ir8_st = ir8_st + isize
      ilinp = ilinp + isize
      ilinp_st = ilinp_st + isize

      tup=tlarge; tup_d=tlarge; tup_w=tlarge; tup_p=tlarge

      !-------------------
      ! now rewind & reread file; capturing information

      rewind lun

      do 
         read(lun,'(A)',iostat=ios) iline
         if(ios.ne.0) exit

         j=0
         do i=1,len(iline)
            if((iline(i:i).ne.' ').and.(iline.ne.char(9))) then
               j=i
               exit  ! non-whitespace found
            endif
         enddo
         if(j.eq.0) cycle   ! ignore all-blank line

         if(iline(j:j).eq.'#') cycle   ! ignore comment line

         if(iline(j:j).eq.'*') then
            ! start new namelist; record length of prior list
            imaxlen_naml = max(imaxlen_naml,ilen_naml)
            ilen_naml = 0

            cur_naml=iline(j+1:j+len(cur_naml))  ! set current namelist 
            call addlist
            cycle
         endif

         ivars=ivars+1
         call read1_db(iline,j,zname,ksteer,ztype,ichsize,iirank,iidims,zdflt)

         varlist(ivars)%name=zname
         call insert_var(ivars)   ! maintain alphabetic order
         varlist(ivars)%steerable=ksteer
         varlist(ivars)%type=ztype
         varlist(ivars)%chsize=ichsize
         ilen_naml = ilen_naml + 1
         varlist(ivars)%naml=cur_naml
         varlist(ivars)%rank=iirank
         varlist(ivars)%dims=iidims
         if(len_trim(zdflt).gt.len_short) then
            long_dflts(ilong)=zdflt
            varlist(ivars)%long_dflt_addr=ilong
            ilong=ilong+1
         else
            varlist(ivars)%short_dflt=zdflt(1:len_short)
         endif
            
         isize=igsize(iirank,iidims)

         varlist(ivars)%nlinadr=ilinp
         ilinp=ilinp+isize
         if(ksteer.eq.2) then
            varlist(ivars)%nlinadr_st = ilinp_st
            ilinp_st=ilinp_st+isize
         endif

         ivecsz=iidims(2,1)-iidims(1,1)+1
         if(ztype.eq.'R') then
            varlist(ivars)%addr=ireal
            call rset_dflt(rbuf(ireal:ireal+isize-1),isize,ivecsz,ivars)
            if(ksteer.eq.2) then
               varlist(ivars)%addr_st = ireal_st
               call rset_ucop(rbuf,ireal,isize,ireal_st,0,0,nupdate_max)
               ireal_st=ireal_st+isize
            endif
            ireal=ireal+isize
         else if(ztype.eq.'I') then
            varlist(ivars)%addr=iint
            call iset_dflt(intbuf(iint:iint+isize-1),isize,ivecsz,ivars)
            if(ksteer.eq.2) then
               varlist(ivars)%addr_st = iint_st
               call iset_ucop(intbuf,iint,isize,iint_st,0,0,nupdate_max)
               iint_st=iint_st+isize
            endif
            iint=iint+isize
         else if(ztype.eq.'L') then
            varlist(ivars)%addr=ilog
            call lset_dflt(logbuf(ilog:ilog+isize-1),isize,ivecsz,ivars)
            if(ksteer.eq.2) then
               varlist(ivars)%addr_st = ilog_st
               call lset_ucop(logbuf,ilog,isize,ilog_st,0,0,nupdate_max)
               ilog_st=ilog_st+isize
            endif
            ilog=ilog+isize
         else if(ztype.eq.'D') then
            varlist(ivars)%addr=ir8
            call dset_dflt(dbuf(ir8:ir8+isize-1),isize,ivecsz,ivars)
            if(ksteer.eq.2) then
               varlist(ivars)%addr_st = ir8_st
               call dset_ucop(dbuf,ir8,isize,ir8_st,0,0,nupdate_max)
               ir8_st=ir8_st+isize
            endif

            ! capture default TINIT
            if(varlist(ivars)%name.eq.'TINIT') then
               tinit=dbuf(ir8)
               tinit_d=tinit
               tinit_w=tinit
               tinit_p=tinit
            endif

            ir8=ir8+isize
         else if(ztype(1:1).eq.'C') then
            varlist(ivars)%addr=ichv
            allocate(chtmp(ivecsz)); chtmp=' '
            call chset_dflt(chbuf(ichv:ichv+isize-1),chtmp, &
                 isize,ivecsz,ivars)
            deallocate(chtmp)
            if(ksteer.eq.2) then
               varlist(ivars)%addr_st = ichv_st
               call chset_ucop(chbuf,ichv,isize,ichv_st,0,0,nupdate_max)
               ichv_st=ichv_st+isize
            endif
            ichv=ichv+isize
         endif

      enddo

      ios=0
      close(unit=lun)

      ! allocate and save copies of default values;
      ! the originals intbuf, etc., will be modified when a real namelist
      ! is read in.

      allocate(logbuf_d(jlog),intbuf_d(jint))
      allocate(rbuf_d(jreal),dbuf_d(jr8))
      allocate(chbuf_d(jchv))

      logbuf_d=logbuf
      intbuf_d=intbuf
      rbuf_d=rbuf
      dbuf_d=dbuf
      chbuf_d=chbuf

      ! allocate write buffers -- to support program driven namelist
      ! editing operations

      allocate(logbuf_w(jlog),intbuf_w(jint))
      allocate(rbuf_w(jreal),dbuf_w(jr8))
      allocate(chbuf_w(jchv))

      logbuf_w=.FALSE.
      intbuf_w=0
      rbuf_w=0
      dbuf_w=0
      chbuf_w=' '

      ! allocate "previous" buffers -- to support comparison-of-namelist
      ! operations

      allocate(logbuf_p(jlog),intbuf_p(jint))
      allocate(rbuf_p(jreal),dbuf_p(jr8))
      allocate(chbuf_p(jchv))

      logbuf_p=.FALSE.
      intbuf_p=0
      rbuf_p=0
      dbuf_p=0
      chbuf_p=' '

      previous = .FALSE.

      ! compile data on namelist-by-namelist basis

      allocate(num_naml_vars(nnamls),num_naml_vars_st(nnamls))
      allocate(indx_naml_vars(imaxlen_naml,nnamls))
      allocate(indx_naml_vars_st(imaxlen_naml,nnamls))

      num_naml_vars = 0
      num_naml_vars_st = 0

      indx_naml_vars = 0
      indx_naml_vars_st = 0

      do i=1,nvars
         j=var_order(i)
         if(varlist(j)%steerable.eq.3) cycle  ! skip ~UPDATE_TIME variable

         cur_naml=varlist(j)%naml
         inaml = 0
         do ii=1,nnamls
            if(cur_naml.eq.all_namls(ii)) then
               inaml=ii
               exit
            endif
         enddo

         if(inaml.eq.0) then
            write(6,*) ' ?splitn -- failed to locate namelist association for:'
            write(6,*) '  ',varlist(j)%name
         endif

         ii = num_naml_vars(inaml) + 1
         indx_naml_vars(ii,inaml) = j
         num_naml_vars(inaml) = ii

         if(varlist(j)%steerable .eq. 2) then
            ii = num_naml_vars_st(inaml) + 1
            indx_naml_vars_st(ii,inaml) = j
            num_naml_vars_st(inaml) = ii
         endif
      enddo

    end subroutine read_db

    integer function igsize(iirank,iidims)
      
      ! calculate size of object

      integer,intent(in) :: iirank   ! rank of object
      integer,intent(in) :: iidims(2,*) ! dimensioning info

      integer idim

      !-----------------

      igsize=1
      if(iirank.eq.0) return

      do idim=1,iirank
         igsize=igsize*(iidims(2,idim)-iidims(1,idim)+1)
      enddo
    end function igsize

    !--------
    subroutine insert_var(ivars)

      integer, intent(in) :: ivars

      ! maintain alphabetic ordering of namelist variables

      integer :: ibefore,imatch,isave,ii

      !--------------------
      isave=nvars    ! save total no. of vars
      nvars=ivars-1  ! this is the no. currently known

      call iorder(varlist(ivars)%name,ibefore,imatch)  ! position this name
      if(imatch.eq.1) then
         call errmsg_exit( &
              ' ?splitn_module: namelist database: duplicate name: '// &
              varlist(ivars)%name)
      else
         do ii=nvars,ibefore,-1
            var_order(ii+1)=var_order(ii)
         enddo
         var_order(ibefore)=ivars
      endif

      nvars=isave   ! restore total no.

    end subroutine insert_var

    subroutine iorder(zname,ibefore,imatch)

      character*(*), intent(in) :: zname  ! name to find
      integer, intent(out) :: ibefore     ! position found
      integer, intent(out) :: imatch      ! =1 if exact match

      character*32 znam32
      integer jj,imin,imax,itry

      !---------------------------------------------

      imatch=0
      ibefore=0

      znam32=zname

      if(nvars.eq.0) then
         ibefore=1
      else 
         jj=var_order(nvars)
         if(znam32.gt.varlist(jj)%name) ibefore=nvars+1
      endif

      if(ibefore.eq.0) then

         ! binary search

         imin=1
         imax=nvars

         do
            if(imin.eq.imax) then
               ibefore=imin
               jj=var_order(imin)
               if(znam32.eq.varlist(jj)%name) imatch=1
               exit
            endif

            itry=(imin+imax)/2

            jj=var_order(itry)
            if(znam32.le.varlist(jj)%name) then
               imax=itry
            else
               imin=itry+1
            endif
         enddo

      endif
    end subroutine iorder

    !--------
    ! routines to set namelist default values in data buffers

    subroutine set_zdflt(ivar,zdflt,ildflt)

      ! get default string for specified namelist variable (whether
      ! long or short form)

      integer, intent(in) :: ivar           ! the namelist item index
      character*(*), intent(out) :: zdflt   ! the fetched string
      integer, intent(out) :: ildflt        ! non-blank length of string

      !---------------------

      if(varlist(ivar)%long_dflt_addr.gt.0) then
         zdflt=long_dflts(varlist(ivar)%long_dflt_addr)
         ildflt=len_trim(zdflt)
      else
         ildflt=len_trim(varlist(ivar)%short_dflt)
         zdflt(1:ildflt)=varlist(ivar)%short_dflt(1:ildflt)
         zdflt(ildflt+1:ildflt+1)=' '
      endif

    end subroutine set_zdflt


    subroutine rset_dflt(ra,ina,iblk,ivar)

      integer, intent(in) :: ina         ! size of array to reset
      real, intent(out) :: ra(ina)       ! the array
      integer, intent(in) :: iblk        ! set iblk array elements at a time
      integer, intent(in) :: ivar        ! variable owning these defaults

      !-------------------------------

      character*20 ibuf
      real blk(iblk),ztmp

      character*512 zdflt
      integer ildflt

      integer i,ios,inum,icomma,icprev

      !-------------------------------

      ! set default string

      call set_zdflt(ivar,zdflt,ildflt)

      ! get value(s) for one block

      icomma=0
      call icomma_get('R',icprev,icomma,zdflt(1:ildflt))
      if(icomma.gt.ildflt) then
         ! just one value to decode
         ibuf=' '
         ibuf(20-ildflt+1:20)=zdflt(1:ildflt)
         call uupper(ibuf)
         if(ibuf(18:20).EQ.'_R8') then
            read(ibuf(1:17),'(G17.0)',iostat=ios) ztmp
         else
            read(ibuf,'(G20.0)',iostat=ios) ztmp
         endif
         if(ios.eq.0) blk(1:iblk)=ztmp
      else
         do i=1,iblk
            inum=icomma-1-icprev
            ibuf=' '
            ibuf(20-inum+1:20)=zdflt(icprev+1:icomma-1)
            call uupper(ibuf)
            if(ibuf(18:20).EQ.'_R8') then
               read(ibuf(1:17),'(G17.0)',iostat=ios) ztmp
            else
               read(ibuf,'(G20.0)',iostat=ios) ztmp
            endif
            if(ios.ne.0) exit
            blk(i)=ztmp
            if(i.lt.iblk) then
               ! char position of next comma-delimiter
               call icomma_get('R',icprev,icomma,zdflt(1:ildflt))
            endif
         enddo
      endif
      if(ios.ne.0) then
         call errmsg_exit( &
              ' ?splitn_module: unexpected decode error: '//zdflt(1:ildflt))
      endif

      ! apply block to whole array

      do i=1,ina,iblk
         ra(i:i+iblk-1)=blk(1:iblk)
      enddo

    end subroutine rset_dflt

    subroutine dset_dflt(ra,ina,iblk,ivar)

      integer, intent(in) :: ina         ! size of array to reset
      real*8, intent(out) :: ra(ina)     ! the array
      integer, intent(in) :: iblk        ! set iblk array elements at a time
      integer, intent(in) :: ivar        ! variable owning these defaults

      !-------------------------------

      character*20 ibuf
      real*8 blk(iblk),ztmp

      character*512 zdflt
      integer ildflt

      integer i,ios,inum,icomma,icprev

      !-------------------------------

      ! set default string

      call set_zdflt(ivar,zdflt,ildflt)

      ! get value(s) for one block

      icomma=0
      call icomma_get('R',icprev,icomma,zdflt(1:ildflt))
      if(icomma.gt.ildflt) then
         ! just one value to decode
         ibuf=' '
         ibuf(20-ildflt+1:20)=zdflt(1:ildflt)
         call uupper(ibuf)
         if(ibuf(18:20).EQ.'_R8') then
            read(ibuf(1:17),'(G17.0)',iostat=ios) ztmp
         else
            read(ibuf,'(G20.0)',iostat=ios) ztmp
         endif
         if(ios.eq.0) blk(1:iblk)=ztmp
      else
         do i=1,iblk
            inum=icomma-1-icprev
            ibuf=' '
            ibuf(20-inum+1:20)=zdflt(icprev+1:icomma-1)
            call uupper(ibuf)
            if(ibuf(18:20).EQ.'_R8') then
               read(ibuf(1:17),'(G17.0)',iostat=ios) ztmp
            else
               read(ibuf,'(G20.0)',iostat=ios) ztmp
            endif
            if(ios.ne.0) exit
            blk(i)=ztmp
            if(i.lt.iblk) then
               ! char position of next comma-delimiter
               call icomma_get('R',icprev,icomma,zdflt(1:ildflt))
            endif
         enddo
      endif
      if(ios.ne.0) then
         call errmsg_exit( &
              ' ?splitn_module: unexpected decode error: '//zdflt(1:ildflt))
      endif

      ! apply block to whole array

      do i=1,ina,iblk
         ra(i:i+iblk-1)=blk(1:iblk)
      enddo

    end subroutine dset_dflt

    subroutine iset_dflt(ra,ina,iblk,ivar)

      integer, intent(in) :: ina         ! size of array to reset
      integer, intent(out) :: ra(ina)    ! the array
      integer, intent(in) :: iblk        ! set iblk array elements at a time
      integer, intent(in) :: ivar        ! variable owning these defaults

      !-------------------------------

      character*20 ibuf
      integer blk(iblk),ztmp

      character*512 zdflt
      integer ildflt

      integer i,ios,inum,icomma,icprev

      !-------------------------------

      ! set default string

      call set_zdflt(ivar,zdflt,ildflt)

      ! get value(s) for one block

      icomma=0
      call icomma_get('R',icprev,icomma,zdflt(1:ildflt))
      if(icomma.gt.ildflt) then
         ! just one value to decode
         ibuf=' '
         ibuf(20-ildflt+1:20)=zdflt(1:ildflt)
         read(ibuf,'(I20)',iostat=ios) ztmp
         if(ios.eq.0) blk(1:iblk)=ztmp
      else
         do i=1,iblk
            inum=icomma-1-icprev
            ibuf=' '
            ibuf(20-inum+1:20)=zdflt(icprev+1:icomma-1)
            read(ibuf,'(I20)',iostat=ios) ztmp
            if(ios.ne.0) exit
            blk(i)=ztmp
            if(i.lt.iblk) then
               ! char position of next comma-delimiter
               call icomma_get('R',icprev,icomma,zdflt(1:ildflt))
            endif
         enddo
      endif
      if(ios.ne.0) then
         call errmsg_exit( &
              ' ?splitn_module: unexpected decode error: '//zdflt(1:ildflt))
      endif

      ! apply block to whole array

      do i=1,ina,iblk
         ra(i:i+iblk-1)=blk(1:iblk)
      enddo

    end subroutine iset_dflt

    subroutine lset_dflt(ra,ina,iblk,ivar)

      integer, intent(in) :: ina         ! size of array to reset
      logical, intent(out) :: ra(ina)    ! the array
      integer, intent(in) :: iblk        ! set iblk array elements at a time
      integer, intent(in) :: ivar        ! variable owning these defaults

      !-------------------------------

      character*20 ibuf
      logical blk(iblk),ztmp

      character*512 zdflt
      integer ildflt

      integer i,ios,inum,icomma,icprev

      !-------------------------------

      ! set default string

      call set_zdflt(ivar,zdflt,ildflt)

      ! get value(s) for one block

      icomma=0
      call icomma_get('R',icprev,icomma,zdflt(1:ildflt))
      if(icomma.gt.ildflt) then
         ! just one value to decode
         ibuf=' '
         ibuf(20-ildflt+1:20)=zdflt(1:ildflt)
         read(ibuf,'(L20)',iostat=ios) ztmp
         if(ios.eq.0) blk(1:iblk)=ztmp
      else
         do i=1,iblk
            inum=icomma-1-icprev
            ibuf=' '
            ibuf(20-inum+1:20)=zdflt(icprev+1:icomma-1)
            read(ibuf,'(L20)',iostat=ios) ztmp
            if(ios.ne.0) exit
            blk(i)=ztmp
            if(i.lt.iblk) then
               ! char position of next comma-delimiter
               call icomma_get('R',icprev,icomma,zdflt(1:ildflt))
            endif
         enddo
      endif
      if(ios.ne.0) then
         call errmsg_exit( &
              ' ?splitn_module: unexpected decode error: '//zdflt(1:ildflt))
      endif

      ! apply block to whole array

      do i=1,ina,iblk
         ra(i:i+iblk-1)=blk(1:iblk)
      enddo

    end subroutine lset_dflt

    subroutine chset_dflt(cha,chblk,ina,iblk,ivar)

      integer, intent(in) :: ina         ! size of array to reset
      character*(*), intent(inout) :: cha(ina)       ! the array
      integer, intent(in) :: iblk        ! set iblk array elements at a time
      character*(*), intent(inout) :: chblk(iblk)    ! a sub-block
      integer, intent(in) :: ivar        ! naml. variable owning these defaults

      !-------------------------------

      character*512 zdflt
      integer ildflt

      integer i,j,ios,inum,icomma,icprev

      !-------------------------------

      ! set default string

      call set_zdflt(ivar,zdflt,ildflt)

      ! get value(s) for one block

      icomma=0
      call icomma_get('C',icprev,icomma,zdflt(1:ildflt))
      if(icomma.gt.ildflt) then
         ! just one value to decode
         ! do not include quotes
         chblk(1:iblk)(1:ildflt-2)=zdflt(2:ildflt-1)
      else
         do i=1,iblk
            inum=icomma-icprev-3
            chblk(i)(1:inum)=zdflt(icprev+2:icomma-2)
            if(i.lt.iblk) then
               ! char position of next comma-delimiter
               call icomma_get('C',icprev,icomma,zdflt(1:ildflt))
            endif
         enddo
      endif

      ! apply block to whole array

      do i=1,ina,iblk
         cha(i:i+iblk-1)=chblk(1:iblk)
      enddo

    end subroutine chset_dflt

    subroutine rset_ucop(buf,isource,inum,itarg,ilin0,inum0,inmax)

      ! propagate values into update blocks

      real, dimension(:) :: buf
      integer :: isource,inum  ! loc & number to copy
      integer :: itarg         ! loc of 1st copy
      integer :: ilin0         ! if .gt.0: check ilines 
      ! (stop propagation of update copy if text line sets a value)

      integer :: inum0         ! .gt.0: copy from update block instead of main
      integer :: inmax         ! index of last copy

      !--------------
      ! local:
      integer :: ii,inc,ism1,itarg_loc,itarg_lin,iinc_lin
      logical, dimension(:), allocatable :: icopy
      logical :: icontin
      !--------------

      allocate(icopy(inum)); icopy=.TRUE.
      icontin=.TRUE.

      if(inum0.eq.0) then
         ism1 = isource-1
      else
         ism1 = nreal + (inum0-1)*nreal_st + itarg - 1
      endif

      itarg_loc = nreal + (inum0-1)*nreal_st + itarg - 1  ! incremented...

      if(ilin0.eq.0) then
         itarg_lin=0
      else
         iinc_lin = nreal_st+nint_st+nlog_st+nr8_st+nchv_st
         itarg_lin = nreal+nint+nlog+nr8+nchv + (inum0-1)*iinc_lin + ilin0 - 1
      endif

      do ii=inum0+1,inmax
         if(ilin0.gt.0) then
            itarg_lin = itarg_lin + iinc_lin
            icontin=.FALSE.
            do inc=1,inum
               if(ilines(itarg_lin+inc).gt.0) icopy(inc)=.FALSE.
               icontin=icontin.or.icopy(inc)
            enddo
         endif
         if(.not.icontin) exit

         itarg_loc = itarg_loc+nreal_st
         do inc=1,inum
            if(icopy(inc)) then
               buf(itarg_loc+inc) = buf(ism1+inc)
            endif
         enddo
      enddo

    end subroutine rset_ucop

    subroutine dset_ucop(buf,isource,inum,itarg,ilin0,inum0,inmax)

      ! propagate values into update blocks

      real*8, dimension(:) :: buf
      integer :: isource,inum  ! loc & number to copy
      integer :: itarg         ! loc of 1st copy
      integer :: ilin0         ! if .gt.0: check ilines 
      ! (stop propagation of update copy if text line sets a value)

      integer :: inum0         ! .gt.0: copy from update block instead of main
      integer :: inmax         ! update block index of last copy

      !--------------
      ! local:
      integer :: ii,inc,ism1,itarg_loc,itarg_lin,iinc_lin
      logical, dimension(:), allocatable :: icopy
      logical :: icontin
      !--------------

      allocate(icopy(inum)); icopy=.TRUE.
      icontin=.TRUE.

      if(inum0.eq.0) then
         ism1 = isource-1
      else
         ism1 = nr8 + (inum0-1)*nr8_st + itarg - 1
      endif

      itarg_loc = nr8 + (inum0-1)*nr8_st + itarg - 1  ! incremented...

      if(ilin0.eq.0) then
         itarg_lin=0
      else
         iinc_lin = nreal_st+nint_st+nlog_st+nr8_st+nchv_st
         itarg_lin = nreal+nint+nlog+nr8+nchv + (inum0-1)*iinc_lin + ilin0 - 1
      endif

      do ii=inum0+1,inmax
         if(ilin0.gt.0) then
            itarg_lin = itarg_lin + iinc_lin
            icontin=.FALSE.
            do inc=1,inum
               if(ilines(itarg_lin+inc).gt.0) icopy(inc)=.FALSE.
               icontin=icontin.or.icopy(inc)
            enddo
         endif
         if(.not.icontin) exit

         itarg_loc = itarg_loc+nr8_st
         do inc=1,inum
            if(icopy(inc)) then
               buf(itarg_loc+inc) = buf(ism1+inc)
            endif
         enddo
     enddo

    end subroutine dset_ucop

    subroutine iset_ucop(buf,isource,inum,itarg,ilin0,inum0,inmax)

      ! propagate values into update blocks

      integer, dimension(:) :: buf
      integer :: isource,inum  ! loc & number to copy
      integer :: itarg         ! loc of 1st copy
      integer :: ilin0         ! if .gt.0: check ilines 
      ! (stop propagation of update copy if text line sets a value)

      integer :: inum0         ! .gt.0: copy from update block instead of main
      integer :: inmax         ! index of last copy

      !--------------
      ! local:
      integer :: ii,inc,ism1,itarg_loc,itarg_lin,iinc_lin
      logical, dimension(:), allocatable :: icopy
      logical :: icontin
      !--------------

      allocate(icopy(inum)); icopy=.TRUE.
      icontin=.TRUE.

      if(inum0.eq.0) then
         ism1 = isource-1
      else
         ism1 = nint + (inum0-1)*nint_st + itarg - 1
      endif

      itarg_loc = nint + (inum0-1)*nint_st + itarg - 1  ! incremented...

      if(ilin0.eq.0) then
         itarg_lin=0
      else
         iinc_lin = nreal_st+nint_st+nlog_st+nr8_st+nchv_st
         itarg_lin = nreal+nint+nlog+nr8+nchv + (inum0-1)*iinc_lin + ilin0 - 1
      endif

      do ii=inum0+1,inmax
         if(ilin0.gt.0) then
            itarg_lin = itarg_lin + iinc_lin
            icontin=.FALSE.
            do inc=1,inum
               if(ilines(itarg_lin+inc).gt.0) icopy(inc)=.FALSE.
               icontin=icontin.or.icopy(inc)
            enddo
         endif
         if(.not.icontin) exit

         itarg_loc = itarg_loc+nint_st
         do inc=1,inum
            if(icopy(inc)) then
               buf(itarg_loc+inc) = buf(ism1+inc)
            endif
         enddo
      enddo

    end subroutine iset_ucop

    subroutine lset_ucop(buf,isource,inum,itarg,ilin0,inum0,inmax)

      ! propagate values into update blocks

      logical, dimension(:) :: buf
      integer :: isource,inum  ! loc & number to copy
      integer :: itarg         ! loc of 1st copy
      integer :: ilin0         ! if .gt.0: check ilines 
      ! (stop propagation of update copy if text line sets a value)

      integer :: inum0         ! .gt.0: copy from update block instead of main
      integer :: inmax         ! index of last copy

      !--------------
      ! local:
      integer :: ii,inc,ism1,itarg_loc,itarg_lin,iinc_lin
      logical, dimension(:), allocatable :: icopy
      logical :: icontin
      !--------------

      allocate(icopy(inum)); icopy=.TRUE.
      icontin=.TRUE.

      if(inum0.eq.0) then
         ism1 = isource-1
      else
         ism1 = nlog + (inum0-1)*nlog_st + itarg - 1
      endif

      itarg_loc = nlog + (inum0-1)*nlog_st + itarg - 1  ! incremented...

      if(ilin0.eq.0) then
         itarg_lin=0
      else
         iinc_lin = nreal_st+nint_st+nlog_st+nr8_st+nchv_st
         itarg_lin = nreal+nint+nlog+nr8+nchv + (inum0-1)*iinc_lin + ilin0 - 1
      endif

      do ii=inum0+1,inmax
         if(ilin0.gt.0) then
            itarg_lin = itarg_lin + iinc_lin
            icontin=.FALSE.
            do inc=1,inum
               if(ilines(itarg_lin+inc).gt.0) icopy(inc)=.FALSE.
               icontin=icontin.or.icopy(inc)
            enddo
         endif
         if(.not.icontin) exit

         itarg_loc = itarg_loc+nlog_st
         do inc=1,inum
            if(icopy(inc)) then
               buf(itarg_loc+inc) = buf(ism1+inc)
            endif
         enddo
      enddo

    end subroutine lset_ucop

    subroutine chset_ucop(buf,isource,inum,itarg,ilin0,inum0,inmax)

      ! propagate values into update blocks

      character*128, dimension(:) :: buf
      integer :: isource,inum  ! loc & number to copy
      integer :: itarg         ! loc of 1st copy
      integer :: ilin0         ! if .gt.0: check ilines 
      ! (stop propagation of update copy if text line sets a value)

      integer :: inum0         ! .gt.0: copy from update block instead of main
      integer :: inmax         ! index of last copy

      !--------------
      ! local:
      integer :: ii,inc,ism1,itarg_loc,itarg_lin,iinc_lin
      logical, dimension(:), allocatable :: icopy
      logical :: icontin
      !--------------

      allocate(icopy(inum)); icopy=.TRUE.
      icontin=.TRUE.

      if(inum0.eq.0) then
         ism1 = isource-1
      else
         ism1 = nchv + (inum0-1)*nchv_st + itarg - 1
      endif

      itarg_loc = nchv + (inum0-1)*nchv_st + itarg - 1  ! incremented...

      if(ilin0.eq.0) then
         itarg_lin=0
      else
         iinc_lin = nreal_st+nint_st+nlog_st+nr8_st+nchv_st
         itarg_lin = nreal+nint+nlog+nr8+nchv + (inum0-1)*iinc_lin + ilin0 - 1
      endif

      do ii=inum0+1,inmax
         if(ilin0.gt.0) then
            itarg_lin = itarg_lin + iinc_lin
            icontin=.FALSE.
            do inc=1,inum
               if(ilines(itarg_lin+inc).gt.0) icopy(inc)=.FALSE.
               icontin=icontin.or.icopy(inc)
            enddo
         endif
         if(.not.icontin) exit

         itarg_loc = itarg_loc+nchv_st
         do inc=1,inum
            if(icopy(inc)) then
               buf(itarg_loc+inc) = buf(ism1+inc)
            endif
         enddo
      enddo

    end subroutine chset_ucop

    subroutine icomma_get(ztype,icprev,icomma,dstr)

      character*1, intent(in) :: ztype   ! "C" -- character, o.w. non-character
      integer, intent(out) :: icprev     ! previous comma location
      integer, intent(inout) :: icomma   ! next comma location (updated)
      character*(*) dstr                 ! default string

      ! find comma separating fields;
      ! if type is "C" also check for enclosing quotes & ignore commas
      ! inside quoted strings.

      character*1 qchar
      integer i,ild

      ild=len(dstr)

      icprev=icomma
      if(ztype.ne.'C') then

         ! non-character string defaults -- look for comma directly.

         icomma=index(dstr(icprev+1:),',')
         if(icomma.le.0) then
            icomma=ild+1
         else
            icomma=icprev+icomma
         endif
         return

      endif

      ! deal with case of character string defaults

      i=icprev+1
      if((dstr(i:i).ne."'").and.(dstr(i:i).ne.'"')) then
         call errmsg_exit( &
              ' ?splitn_module: expected quote character at start or after'// &
              ' comma: '//dstr)
      endif

      qchar=dstr(i:i)
      do
         i=i+1
         if(i.eq.ild) then
            if(dstr(i:i).ne.qchar) then
               call errmsg_exit( &
                    ' ?splitn_module: namelist database: no closing quote: '//&
                    dstr)
            else
               icomma=ild+1
               exit
            endif
         endif
         if(dstr(i:i).eq.qchar) then
            if(dstr(i+1:i+1).eq.',') then
               icomma=i+1
               exit
            else
               call errmsg_exit( &
                    ' ?splitn_module: expected comma after closing quote: '//&
                    dstr)
            endif
         endif
      enddo
    end subroutine icomma_get

    subroutine read1_db(iline,istart,zname, &
         ksteer,ztype,ichsize,iirank,iidims,zdflt)

      ! parse a line from the database; read continuation lines to
      ! complete parsing if necessary.

      ! error checking is minimal: the database is supposed to have
      ! been written by code generators that have been well debugged.

      character*(*), intent(inout) :: iline   ! input line
      integer, intent(in) :: istart           ! 1st non-blank char in iline

      character*(*), intent(out) :: zname     ! name of namelist variable
      integer, intent(out) :: ksteer          ! T if item is steerable
      character*(*), intent(out) :: ztype     ! type of namelist variable
      integer, intent(out) :: ichsize         ! size of C*nn vars
                                              ! for non-char data: ichsize=0
      integer, intent(out) :: iirank          ! rank (dimensionality) of var.
      integer, intent(out) :: iidims(2,4)     ! actual dimensions
      character*(*), intent(out) :: zdflt     ! namelist var. default values

      !------------------------------------------------
      integer :: ifld   ! field ctr
      integer :: ist    ! start of field
      integer :: ilast_fld  ! last field count (default string)
      integer i,inc,ilen
      integer idim,ic1,ic2,il,icolon,ilz
      character*10 ibuf
      !------------------------------------------------

      iirank=0
      iidims=0
      ichsize=0

      ist=istart
      ifld=0

      ilen=len(iline)

      do
         ! loop over fields

         i=ist

         do
            ! loop over characters: find whitespace
            i=i+1
            if(iline(i:i).eq.' ') exit
            if(iline(i:i).eq.char(9)) exit
            if(i.eq.ilen) then
               call errmsg_exit(' ?splitn_module: missing fields: '//iline)
            endif
         enddo

         ifld = ifld+1
         if(ifld.eq.1) then
            ! name field
            zname=iline(ist:i-1)
            call uupper(zname)

         else if(ifld.eq.2) then
            ! type field (R/I/L/D/C*n)
            ! steerability flag added (dmc Aug 2005)

            ksteer=0
            if(iline(ist:ist).eq.'~') then
               if(cur_naml(1:3).eq.'TRN') then
                  ksteer=2
               else
                  ksteer=1
               endif
               ist=ist+1
            endif

            if(i.gt.(ist+5)) then
               call errmsg_exit(' ?splitn_module: error in type field: '// &
                    iline)
            endif
            ztype=iline(ist:i-1)
            if(ztype(1:2).eq.'C*') then
               ! get size n of C*n variable
               if(i.eq.(ist+3)) then
                  read(iline(ist+2:ist+2),'(I1)') inc
               else if(i.eq.(ist+4)) then
                  read(iline(ist+2:ist+3),'(I2)') inc
               else if(i.eq.(ist+5)) then
                  read(iline(ist+2:ist+4),'(I3)') inc
               endif
               if((inc.le.0).or.(inc.gt.ch_maxlen)) then
                  call errmsg_exit( &
                       ' ?splitn_module: C*n: n<0 or n>ch_maxlen: '//iline)
               endif
               ichsize = inc
            endif

         else if(ifld.eq.3) then
            if(i.gt.ist+1) then
               call errmsg_exit(' ?splitn_module: error in rank field: '// &
                    iline)
            else
               read(iline(i-1:i-1),'(I1)') iirank
               if(iirank.gt.4) THEN
                  call errmsg_exit(' ?splitn_module: exceeded max rank=4: '// &
                       iline)
               else
                  ilast_fld = 4+iirank
               endif
            endif
                  
         else if(ifld.eq.ilast_fld) then
            if(iline(ist:ist).ne.'/') then
               call errmsg_exit( &
                    ' ?splitn_module: cannot find default value string: '// &
                    iline)
            endif

            exit       ! parse default value string in separate loop

         else
            idim = ifld-3   ! which dimension this is...

            ! expect field of form "(n1:n2)"

            if(iline(ist:ist).ne.'(') call errmsg_exit( &
                 ' ?splitn_module: expected "(" at start of dimension field: '//iline)
            if(iline(i-1:i-1).ne.')') call errmsg_exit( &
                 ' ?splitn_module: expected ")" at end of dimension field: '//iline)
            icolon=index(iline(ist:i-1),':')
            if(icolon.le.0) call errmsg_exit( &
                 ' ?splitn_module: expected ":" in dimension field: '//iline)
            icolon=icolon+ist-1
            ic1=ist+1
            ic2=max(ic1,icolon-1)
            il=ic2-ic1+1
            ibuf=' '
            ibuf(10-il+1:10)=iline(ic1:ic2)
            read(ibuf,'(I10)') iidims(1,idim)

            ic1=icolon+1
            ic2=max(ic1,i-2)
            il=ic2-ic1+1
            ibuf=' '
            ibuf(10-il+1:10)=iline(ic1:ic2)
            read(ibuf,'(I10)') iidims(2,idim)

         endif

         ! find start of next field

         do 
            i=i+1
            if((iline(i:i).ne.' ').and.(iline(i:i).ne.char(9))) exit
            if(i.eq.ilen) then
               call errmsg_exit(' ?splitn_module: missing fields: '//iline)
            endif
         enddo

         ist=i       ! non-blank start of next field

      enddo

      ! get default value string

      ilen = len_trim(iline)
      ist = ist+1

      ! first non-blank

      do i=ist,ilen
         if((iline(i:i).ne.' ').and.(iline(i:i).ne.char(9))) exit
      enddo
      if(iline(i:i).eq.'/') call errmsg_exit( &
           ' ?splitn_module: empty default value string: '//iline)

      ilz=1
      ist=i

      ! check for continuation line; check for end; read if necessary

      do
         if(iline(ilen:ilen).eq.'/') then
            ! NO continuation
            zdflt(ilz:)=iline(ist:ilen-1)
            exit

         else if(iline(ilen:ilen).eq.char(92)) then
            ! backslash: continuation
            zdflt(ilz:)=iline(ist:ilen-1)

         else 
            call errmsg_exit( &
                 ' ?splitn_module: default value string syntax error: '// &
                 iline)
         endif

         ilz=len_trim(zdflt) + 1
         read(lun,'(A)') iline

         ilen=len_trim(iline)
         do i=1,ilen
            if((iline(i:i).ne.' ').and.(iline(i:i).ne.char(9))) exit
         enddo
         if(i.eq.ilen) then
            if(iline(i:i).eq.'/') exit  ! terminator found.
         endif

         ist=i  ! and loop back to insert text from this line
      enddo

    end subroutine read1_db

    subroutine addlist

      !  add namelist "cur_naml" to list of namelists, with the following
      !  caveats which are related to backwards compatibility with the
      !  old "splitn" software:
      !
      !   (a) ignore "GENSIS"
      !   (b) put "TRNCHR" at beginning of list
      !   (c) put "TRDATA" at end of list
      !   (d) put TRN001--TRNnnn after "TRNCHR", in order found
      !   (e) put everything else before "TRDATA", in order found.

      integer i,ipos,itest,ios
      integer :: innn=0
      save innn

      !--------------------------------------------
      if(nnamls.eq.0) innn=0

      ! ...ignore GENSIS
      if(cur_naml.eq.'GENSIS') return

      ! ...ignore duplicates, although these are not expected.
      do i=1,nnamls
         if(cur_naml.eq.all_namls(i)) return
      enddo

      ! ...see if form of namelist name is TRNnnn, nnn = 3 digit integer
      itest=-1
      if(cur_naml(1:3).eq.'TRN') then
         read(cur_naml(4:6),'(I3)',iostat=ios) itest
         if(ios.ne.0) itest=-1
      endif

      ! ...find position in list for this namelist
      if(cur_naml.eq.'TRNCHR') then
         ipos=1
         innn=innn+1
      else if(cur_naml.eq.'TRDATA') then
         ipos=nnamls+1
      else if(itest.gt.-1) then
         ! TRNnnn namelist
         innn=innn+1
         ipos=innn
      else
         ! ...none of the above... insert before TRDATA
         if(nnamls.eq.0) then
            ipos=1
         else if(all_namls(nnamls).eq.'TRDATA') then
            ipos=nnamls
         else
            ipos=nnamls+1
         endif
      endif

      ! insert namelist name now

      do i=nnamls,ipos,-1
         all_namls(i+1)=all_namls(i)
      enddo
      all_namls(ipos)=cur_naml

      nnamls=nnamls + 1

    end subroutine addlist

    logical function is_deleted(zname)

      character*(*), intent(in) :: zname

      !  return TRUE if passed name is a known deleted name

      !--------------------------------
      character*32 :: znam32
      integer :: inam
      !--------------------------------

      znam32 = zname
      call uupper(znam32)

      is_deleted = .FALSE.

      do inam=1,ndnames
         if(znam32.eq.dnames(inam)) then
            is_deleted = .TRUE.
            exit
         endif
      enddo

    end function is_deleted

end module splitn_module
