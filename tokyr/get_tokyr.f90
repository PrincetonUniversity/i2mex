subroutine get_tokyr(lunfile,lunmsg,ztok,nshot,zyear)

  ! given tokamak id and shot number, return shot-year, if possible.
  ! if unable to do so return zyear=' '.

  ! History:
  !
  ! CLF July 2012 -- try <tok>_YEAR_FILE first
  !                              to accomodate G. Tardini
  ! DMC July 2009 -- added explicito convertor for D3D!

  ! for NSTX, TFTR and JET, explicit convertor routines are maintained.
  ! for other tokamaks:  an environment variable
  !    <tok>_YEAR_FILE
  ! can be defined, pointing to a file containing an ordered list
  ! of records <SHOT#> <YEAR> where each record gives the first SHOT#
  ! in the given year for the given <tok>.  The "#" character can be
  ! used to denote a comment line or comments after a valid record.
  ! blank lines are ignored.  A valid file would look like (ignore
  ! leading "  ! ":

  ! # shot-year data for XYZ tokamak
  !
  !   1      1998
  !   10001  1999   # first shot of 1999
  !   12345  2000   # first shot of 2000
  !   23456  2001   # first shot of 2001
  ! <EOF>

  ! (the above file in some generally accessible place and pointed to
  !  by the XYZ_YEAR_FILE environment variable).
  !
  ! obviously, such a file would need to be updated about once a year
  ! for an operating tokamak... shot numbers can be up to 6 digits long.

  ! **modification DMC June 2006** If shot year is not found by above method,
  !   a backup year may be defined from the PRETR_YY environment variable.

  implicit none

  integer, intent(in) :: lunfile      ! FORTRAN LUN for shot year data file
  ! (file may be opened on this LUN).

  integer, intent(in) :: lunmsg       ! FORTRAN LUN for messages (or 0 to
  ! suppress messages)

  character(len=*), intent(in) :: ztok   ! tok id, e.g "JET" or "NSTX".
  integer, intent(in) :: nshot        ! shot number

  character(len=*), intent(out) :: zyear ! shot year (output)

  ! tok ids are general 3 or 4 characters long, all CAPITALS

  ! if len(zyear).lt.4) the last two digits of the shot year are
  !   placed right justified in zyear (w/ leading blank if need be)
  ! otherwise the full four digit shot year is placed right justified
  !   in zyear (w/ leading blanks as needed).

  ! ============================================

  integer :: iltok,ios,ils,ily,ilz,ilen,ilen_3
  character(len=20) :: envar  ! environment variable
  character(len=100) :: yfile ! filename

  character(len=80) :: zline
  integer :: iishot,i1,i2,i3,i4,iycmod,ilyy
  character(len=4) :: zycmod
  character(len=10) :: zzyear,zzyearp
  character(len=20) :: zzshot

  ! for some tokamaks a shot-year routine is explicitly maintained...

1001 format(' %get_tokyr: subroutine ',a,': tokamak ',a, &
       ' shot ',i6,' => shot year: ',a)
1002 format(' %get_tokyr: using ',a,': tokamak ',a, &
       ' shot ',i6,' => shot year: ',a)

  ilen=len(zyear)
  ilen_3 = max(1,ilen-3)    ! RGA, out of range errors


  ! look in a file indicated by environment variable <tok>_YEAR_FILE
  ! if it exists.

  zyear=' '

  yfile=' '
  iltok=len_trim(ztok)
  envar=ztok(1:iltok)//'_YEAR_FILE'
  call get_environment_variable(trim(envar), yfile)

  if(yfile.ne.' ') then
    ! try to open file

     open(unit=lunfile,file=yfile,status='OLD',iostat=ios,action='READ')
     if(ios.ne.0) then
        close(unit=lunfile,iostat=ios)
        if(lunmsg.gt.0) then
           write(lunmsg,*) ' %get_tokyr: ',envar,' = ',yfile
           write(lunmsg,*) '  but file could not be opened.'
           write(lunmsg,*) '  => trying alternative'
        endif
     else
        
        ! read file; expect order list of records:
        ! <shot #> <year>
        ! meaning <shot #> is the first shot of year <year>.
        
        ilz=len(zzshot)
        zzyear='????'
        do
           read(lunfile,'(A)',iostat=ios) zline
           if(ios.ne.0) exit  ! assuming EOF
           
           call parslin(zline,i1,i2,i3,i4) ! find 2 fields
           if(i1.eq.0) cycle
           if(zline(i1:i1).eq.'#') cycle
           if(i3.eq.0) then
              if(lunmsg.gt.0) then
                 write(lunmsg,*) ' %get_tokyr: in file: ',yfile
                 write(lunmsg,*) '  missing field in line: ',zline
                 cycle
              endif
           endif
           
           ils=i2-i1+1
           zzshot=' '
           zzshot(ilz-ils+1:ilz)=zline(i1:i2)
           read(zzshot,'(i20)',iostat=ios) iishot
           if(ios.ne.0) then
              if(lunmsg.gt.0) then
                 write(lunmsg,*) ' %get_tokyr: in file: ',yfile
                 write(lunmsg,*) '  integer decode error in line: ',zline
                 cycle
              endif
           endif
           
           if(zzyear.eq.'????') zzyear=zline(i3:i4)
           zzyearp=zzyear
           zzyear=zline(i3:i4)
           
           if(iishot.gt.nshot) then
              zzyear=zzyearp
              exit
           endif
        enddo
        
        close(unit=lunfile)   ! done with file.
        
        if(zzyear.eq.'????') then
           if(lunmsg.gt.0) then
              write(lunmsg,*) ' %get_tokyr: found file: ',yfile
              write(lunmsg,*) '  but could not determine shot year.'
           endif
           
        else
           
           ! year extracted OK
           
           ily=len_trim(zzyear)
           ilen=len(zyear)
           if(ilen.le.3) then
              zyear(ilen-1:ilen)=zzyear(ily-1:ily)
           else
              zyear(ilen-3:ilen)=zzyear(ily-3:ily)
           endif
           if(lunmsg.gt.0) then
              ilz=len_trim(envar)
              write(lunmsg,1002) envar(1:ilz),ztok,nshot,zyear(ilen-3:ilen)
           endif
        endif
     endif
  endif
  
  
  if(zyear.eq.' ') then
     if(ztok.eq.'JET') then
        call jetyr(nshot,zyear)
        if(lunmsg.gt.0) write(lunmsg,1001) 'jetyr.for',ztok,nshot, &
             zyear(ilen_3:ilen)
     else if(ztok.eq.'MAST') then
        call mastyr(nshot,zyear)
        if(lunmsg.gt.0) write(lunmsg,1001) 'mastyr.for',ztok,nshot, &
             zyear(ilen_3:ilen)
    else if(ztok.eq.'TFTR') then
       call tftryr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'tftryr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'NSTX') then
       call nstxyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'nstxyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'NSTU') then
       call nstuyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'nstuyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'TCV') then
       call tcvyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'tcvyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'D3D') then
       call d3dyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'd3dyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'DEMO') then
       call demoyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'demoyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'AUGD') then
       call augdyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'augdyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'EAST') then
       call eastyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'eastyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'KSTR') then
       call kstryr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'kstryr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'HL2M') then
       call hl2myr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'hl2myr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'ITER') then
       call iteryr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'iteryr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'ARC') then
       call arcyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'arcyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'ARIS') then
       call arisyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'arisyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'SSSP') then
       call ssspyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'ssspyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'LTX') then
       call ltxyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'ltxyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'GLMT') then
       call glmtyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'glmtyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'STEP') then
       call stepyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'stepyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'SPRC') then
       call sprcyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'sprcyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'ST40') then
       call st40yr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'st40yr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'ST80') then
       call st80yr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'st80yr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'STE1') then
       call ste1yr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'ste1yr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'PI3') then
       call pi3yr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'pi3yr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'SMRT') then
       call smrtyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'smrtyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'SHPD') then
       call shpdyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'shpdyr.f90',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'TCAB') then
       call tcabyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'tcabyr.f90',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'WEST') then
       call westyr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'westyr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if(ztok.eq.'STAR') then
       call staryr(nshot,zyear)
       if(lunmsg.gt.0) write(lunmsg,1001) 'staryr.for',ztok,nshot, &
            zyear(ilen_3:ilen)
    else if((ztok.eq.'CMOD').and.(nshot.gt.9999999)) then
       !  assume CMOD shot encodes year...
       iycmod=nshot/10000000
       iycmod=iycmod+1900
       write(zycmod,'(i4.4)') iycmod
       zyear=' '
       i4=len(zyear)
       if(i4.lt.4) then
          zyear(i4-1:i4)=zycmod(3:4)
       else
          zyear(i4-3:i4)=zycmod
       endif
    else 
       if(zyear.eq.' ') then
          
          !  DMC June 2006 -- if YY not found yet, try environment variable...
          
          ilen=min(4,len(zyear))
          zline=' '
          call get_environment_variable("PRETR_YY", yfile)
          if(zline.ne.' ') then
             ilyy=len(trim(zline))
             if(ilyy.le.ilen) then
                zyear(ilen-ilyy+1:ilen)=zline(1:ilyy)
             endif
          endif
       endif
    endif
    
 endif
 
 if(zyear.eq.' ' .and. yfile.eq.' ') then
    
    ! environment variable UNDEFINED
    
    if(lunmsg.gt.0) then
       write(lunmsg,*) ' %get_tokyr: not defined: ',envar
       write(lunmsg,*) '  => unable to set shot year for ',ztok,' ',nshot
    endif
 endif
 
contains
  
  subroutine parslin(zline,i1,i2,i3,i4)
    
    ! find 1st and 2nd field in line (if possible)
    
    character(len=*), intent(in) :: zline
    
    integer, intent(out) :: i1,i2  ! zline(i1:i2) -- first field
    integer, intent(out) :: i3,i4  ! zline(i3:i4) -- second field
    
    integer i,ilen
    
    i1=0; i2=0; i3=0; i4=0
    if(zline.eq.' ') return
    
    ilen=len_trim(zline)
    do i=1,ilen
       if((zline(i:i).ne.' ').and.(zline(i:i).ne.char(9))) then
          i1=i
          exit
       endif
    enddo
    if(i1.eq.0) return
    
    do i=i1,ilen
       if((zline(i:i).eq.' ').or.(zline(i:i).eq.char(9))) then
          i2=i-1
          exit
       endif
    enddo
    if(i2.eq.0) then
       i2=ilen
       return
    endif
    
    do i=i2+1,ilen
       if((zline(i:i).ne.' ').and.(zline(i:i).ne.char(9))) then
          i3=i
          exit
       endif
    enddo
    
    i4=ilen
    do i=i3,ilen
       if((zline(i:i).eq.' ').or.(zline(i:i).eq.char(9))) then
          i4=i-1
          exit
       endif
    enddo
    
  end subroutine parslin
  
end subroutine get_tokyr
