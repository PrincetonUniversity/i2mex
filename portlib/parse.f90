  subroutine parse( cbuf, istat )
  !.. parse                            Resolve variable references in a string
  !                                    codesys/source/portlib/parse.
  !
  !
  !  Author                            jim.conboy@ccfe.ac.uk
  !  Version                           1.00,  07Dec2011
  !  Modifications
  !
  !  1.00  07Dec2011                   jim.conboy@ccfe.ac.uk
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
  character(len=*),intent(inout)     :: cbuf           ! inp
  integer, intent(out)               :: istat          ! return status, 0=OK
  !
  !..local
  integer                            :: ixi1, ixi2   & ! index to input string
                                       ,ixo          & !          output
                                       ,ixv1, ixv2   & ! 1st, last char of the variable name
                                       ,li, lo         ! lengths of strings
  !
#ifdef _DEBUG
  logical                            :: ldebug = .true.
#else
  logical                            :: ldebug = .false.
#endif
  character(len=2048)                :: cout           ! Output string
  character(len=64)                  :: cVnms          ! Characters allowed in variable names
  !
  character(len=1)                   :: cV           & ! Variable identifier
                                       ,cV1, cV2       !    "     delimiter
  !
  character(len=3), parameter        :: cVdf = '${}'
  !
  character(len=1), parameter        :: cUs  = '_'   &
                                       ,cBl  = ' '   &
                                       ,cQ   = '"'
  !
  character(len=*), parameter        :: cU = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'  &
                                       ,cL = 'abcdefghijklmnopqrstuvwxyz' 
  character(len=*),parameter         :: cSr='portlib/parse'
  !----------------------------------::================!.........................|
  !
  do ixi1=1,len(cout), 8  ; cout(ixi1:)  = '        ' ; enddo
  do ixi1=1,len(cVnms),8  ; cVnms(ixi1:) = '        ' ; enddo
  !
  cVnms = cU//cL//cUs
  cV = cVdf(1:1) ; cV1 = cVdf(2:2) ; cV2 = cVdf(3:3)
  !
  istat = 0
  ixi1 = 0 ; ixi2 = 1 ; ixo = 1
  l_var: do
     i_found: &
       if( ixi2 >= len_trim(cbuf ) )                 then
       !    end of buffer
                                                     exit l_var
       elseif( cbuf(ixi2:ixi2) == cV )               then
       !    Special case - string commences with variable
       ixi1 = ixi2
                                                        else
       ! look for next variable
       ixi1 = index(  cbuf(ixi2:), cV ) + ixi2 - 1
  !
     if( ixi1 == ixi2-1 )                       then
  !     No more variables ..
        ixi1 = len_trim( cbuf(ixi2:) ) + ixi2 
  !d    print *,'End "',cbuf(ixi2:ixi1),cQ
        cout(ixo:) = cbuf( ixi2:ixi1 )
                                                exit l_var
                                                else
  !  collect preceeding text
           cout(ixo:)  = cbuf(ixi2:ixi1-1)
           ixo = ixo + ixi1-ixi2
                                                endif
                                                        endif i_found
  !  variable delimiter ($) in cbuf(ixi1)
  !d print *,ixi1, cbuf(ixi1:ixi1+3)
     i_delim: if( cbuf(ixi1+1:ixi1+1) == cV1 )          then
     !  form is ${VARIABLE}
        ixv1 = ixi1 + 2
        ixv2 = index(  cbuf(ixi1:), cV2 ) + ixi1 -2
        if( ldebug) print *,ixv1,ixv2,cbuf(ixv1:ixv2)
        if( ixv2 < ixv1 )                       then
           istat = 2                                                    ! istat = 2 ; missing delimiter
                                                exit l_var
                                                else
           ixi2 = ixv2 + 2
                                                endif
        
                                                        else
     !     end delimiter is any invalid name character
           ixv1 = ixi1 +1
           ixv2 = ixv1  &
                + verify( trim(cbuf(ixv1:)), cVnms ) - 2
           if( ixv2 == ixv1-2 )                 then
     !          special case - variable name terminates string; ixi2 points to last char
                ixv2 = ixv1 + len_trim(cbuf(ixv1:)) -1
                ixi2 = ixv2
                                                else
     !          point to start of next substring
                ixi2 = ixv2+1
                                                endif
                                                        endif i_delim
     !
     !  cbuf(ixv1:ixv2) is the variable name
           if( ldebug ) print *,ixi1,ixv1,ixv2,ixi2, cQ,cbuf(ixv1:ixv2),cQ
           if( ixv2 < ixv1 )                            then
              istat = 3                                                 ! istat = 3 ; zero length name
                                                        exit l_var
                                                        endif
     !====     Name Translation here    ====
#ifdef __MPI
           call mpi_sget_env(  cbuf(ixv1:ixv2), cout(ixo:), istat )
#else
           call get_environment_variable( cbuf(ixv1:ixv2), cout(ixo:) )
#endif
           li = len_trim(cout(ixo:))
           if( li == 0 )                        then
     !         no string returned..
               cout(ixo:) = cbuf(ixi1:ixv2)
               istat = 1                                                ! istat = 1 ; name not defined
                                                endif
           ixo = ixo + li
           if( ldebug ) print *,ixi1,ixi2,ixo,cQ,cout(:ixo-1),cQ
     enddo l_var
     !
     ixo = len_trim(cout)
     if( istat > 0 )   &
       write(6,'(4x,2a,i3)')      cSr,' -W- istat = ',istat

     if(ldebug .or.  istat > 1 )                        then
        cout(ixo+ixi1:) = '1'
        cout(ixo+ixi2:) = '2'
        if( ixv1 > 0 )   cout(ixo+ixv1:) = 'v-'
        if( ixv2 > 0 )   cout(ixo+ixv2-1:) = '-v'
        write(6,100)    cQ,trim(cbuf),cQ
        write(6,100)    cQ,cout(ixo+1:ixo+len_trim(cbuf)),cQ
                                                        endif
!
     cbuf = cout(:ixo)
!
100 format(4x,3a)
  end subroutine parse
#ifdef __TEST
  program parse_test
  !.. parse_test                       quick test of parse
  !                                    codesys/source/portlib/parse_test.
  !
  !
  !  Author                            <author>@<place>
  !  Version                           1.00,  08Dec2011
  !  Modifications
  !
  !  1.00  08Dec2011                   <author>@<place>
  !                                    Written
  !
  !  Usage
  !     >cpp --traditional -D__TEST -D_DEBUG =1/parse.f90 ./parse.f90 
  !     >ifort -g parse.f90 -l:$CS/lib/portlib.a -o parse
  !     >./parse
  !  OR
  !     >cpp --traditional -D__TEST -D_DEBUG -D__MPI =1/parse.f90 ./parse.f90 
  !     > mpif90 -g parse.f90 -I$CS/obj/mpi_portlib -I $CS/mod -l:$CS/lib/mpi_portlib.a -o mpi_parse
  !----_^---------------------------------=========================================|
  !
#ifdef __MPI
  use mmpi
  use mpi_env_mod
#endif
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..local
  integer                                   :: istat
  character(len=1)                          :: cQ = '"'
  character(len=256)                        :: cbuf
  !
  integer                            :: iex            ! loop index 
  integer, parameter                 :: nex  = 8     & ! # of examples
                                       ,nexm = 1       !      mpi examples
#ifdef __MPI
  logical                            :: lmpi = .true.
#else
  logical                            :: lmpi = .false.
#endif
  !
  character(len=255),dimension(nex)  :: cex = &
                                       (/'No/variables/here           ' &
                                        ,'/tmp/$USER/Transp           ' &
                                        ,'/tmp/${USER}_files          ' &
                                        ,'${USER}_on_$HOST            ' &
                                        ,'$USER$HOST                  ' &
                                        ,'Undefined-$NO_SUCH_VARIABLE!' &
                                        ,'MissingRB-${HOME/directory  ' &
                                        ,'$_cwd/Toric                 ' &
                                                                       /)
  !
  character(len=*),parameter         :: cSr='portlib/parse_test'
  !----------------------------------::================!.........................|
  !
#ifdef __MPI
  call mmpi_init(istat)
  !
  call env_mod_init( istat )
  call mpi_proc_env( istat )
  !
  call mpi_env_prt( 6 )
  !
#endif
  l_ex: do iex = 1,nex
     call parse_one( cex(iex), istat )
      write(6,*) '----------------------------------------------------------------'
     if( .not. lmpi .and. iex+nexm == nex )   exit l_ex
  enddo l_ex
  !
#ifdef __MPI
  call mmpi_end( istat )
#endif
  !
  contains

  subroutine parse_one( cbuf, istat )
  !.. parse_one                        runs one example
  !                                    codesys/source/portlib/parse_test.parse_one.
  !
  !
  !  Author                            <author>@<place>
  !  Version                           1.00,  08Dec2011
  !  Modifications
  !
  !  1.00  08Dec2011                   <author>@<place>
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
  character(len=*), intent(inout)    :: cbuf           ! string to parse
  integer, intent(out)               :: istat          ! return status, 0=OK
  !
  !..local
  integer                            :: il             ! length of string
  !
  character(len=*), parameter        :: ckey = '123456789012345678901234567890'
  character(len=*), parameter        :: cSr='portlib/parse_test.parse_one'
  !----------------------------------::================!.........................|
  !
  il = len_trim(cbuf)
  write(6,100) cQ,cbuf(:il),cQ
#ifdef _DEBUG
  write(6,100) cQ,ckey(:il),cQ
#endif
  call parse( cbuf, istat )
  write(6,100) cQ,trim(cbuf),cQ,istat
  !
100 format(10x,3a)
  end subroutine parse_one

  end program parse_test
#endif
  
