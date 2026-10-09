    subroutine splitn_io_list(  ixspl, cspl, nspl, ncspl, istat )
!.. splitn_io_list                      Return list of variables defining transp inputs
!
!
!     Author                            jim.conboy@jet.uk
!     Version                           1.11,  29Aug2007
!     Modifications
!
!     1.11  29Aug2007                   jim.conboy@jet.uk
!                                       Extracted from splitn_module
!
!     1.10  31Jul2007                   jim.conboy@jet.uk
!                                       ReWritten
!----_^---------------------------------=========================================|

      use splitn_module                                ! namelist interface ( iorder, splitn_cget
!
      implicit  none
!
!..GLOBAL
!--   variable                       :: name           ! description
!
!..Arguments
!
      integer, intent(inout)         :: ixspl          ! index to next entry
                                                       ! 0  => Start
                                                       ! -1 => end
      integer, intent(in)            :: nspl, &        ! dimension of cspl
                                        ncspl          ! length of (element of ) cspl

      character(len=ncspl),                        &
       dimension(1:nspl),intent(out) :: cspl           ! Array of
                                                       ! Ctri//Cpre//Cext

      integer, intent(out)           :: istat          ! return status, 0=OK
!
!..local
      integer, save                  :: iobg,  &       ! 1st io index
                                        ioen           ! Last io index
!
      integer                           i,     &
                                        im1,   &       ! 1st match flag
                                        im2,   &       ! 2nd match flag
                                        iv,    &       ! variable index
                                        ivo,   &       ! ordered variable ix
                                        idt,   &       ! chdatum index
                                        isz,   &       ! char array size
                                        ios,   &       ! io status
                                        ixcspl         ! index into cspl
!
      character(len=64)              :: cbuf
      character(len=nc_pre)         & 
                       ,dimension(8) :: cdatpre        ! Returned pre array
      character(len=nc_ext)         &
                       ,dimension(8) :: cdatext
!
      character(len=16)              :: cbl  = '                '
      character(len=12)               :: cix  = ' 12345678   '
      character(len=*),parameter     :: zcr='splitn_io_list'
!     character(len=*),parameter     :: cStop=cM//cr
!----_^------------------------------::================!.........................|
!
! nspl  =  size( cspl, 1 )                             ! only if in module
!
  if( ixspl == 0 )                                      then
     
!  find start & end of list

     call iorder( 'PREAAA', iobg, im1 )
     call iorder( 'PREZZZ', ioen, im2 )
                                                        endif

!d  write(6,*) iobg, ioen, im1, im2
 

     ixcspl  =  0
     istat   =  1
l_spec:  do                              ! loop over potential file specifiers
!
     idt = mod(ixspl, 8)
     if( idt == 0 )                      then
         ixspl = ixspl + 8
     else
         idt = idt + 1
     endif
     iv  = iobg+ixspl/8
!d   print *,ixcspl, ixspl/8
     if( iv >= ioen )                    then
!                          All done
        do i=ixcspl+1,nspl
           cspl(i)(1:4) = '   '
        enddo
        istat  =  0 ; ixcspl = -1
                                         exit l_spec
                                         endif
!
     ivo = var_order(iv)
     if( varlist(ivo)%rank > 0 )         then
       isz = varlist(ivo)%dims(2,1)
       if( idt == 0 ) idt = 1
       if( idt > isz)                 then
          idt = 0
                                      cycle l_spec
                                      endif
                                         else
       isz = 1 
       idt = 0
                                         endif
     ixspl = ixspl + idt
!
     if(     ixcspl == 0                     &
        .or. idt < 2           )                                then
        cdatpre = cbl
        cdatext = repeat(cbl(:8),nc_ext/8)
        call splitn_cget(varlist(ivo)%name,  &
                      len(cdatpre(1)), isz, cdatpre, ios)
        if( ios /= 0 ) &
          print *,zcr,': ',iv,ivo, varlist(ivo)%name

        call splitn_cget('EXT'//varlist(ivo)%name(4:6),  &
                      len(cdatext(1)), isz, cdatext, ios)
        if( ios /= 0 ) &
          print *,zcr,': ',iv,ivo, 'EXT'//varlist(ivo)%name(4:6) &
                          ,varlist(ivo)%name, cdatpre
                                                                endif
     if( cdatpre(max(1,idt))(:8) .eq. 'ZZZZZZZZ' )              then
        ixspl = 8*(ixspl/8)
                                                                cycle l_spec  ! skip if defaulted
                                                                endif

!
     ixcspl  =  ixcspl + 1
!d*  print *,ixcspl,ixspl,idt,ivo,varlist(ivo)%name
     cspl(ixcspl)(              1:nc_tri+1)        = varlist(ivo)%name(4:6)//cix(idt+1:idt+1)  ! trigraph
     cspl(ixcspl)(       nc_tri+2:nc_tri+1+nc_pre) = cdatpre(max(1,idt))                       ! pre
     cspl(ixcspl)(nc_tri+nc_pre+2:)                = cdatext(max(1,idt))                       ! ext
!
     if( ixcspl >= nspl)                 then
                                         exit l_spec
                                         endif
!
!d   print *,trim(cspl(ixcspl))
!
  enddo  l_spec

  end subroutine splitn_io_list

