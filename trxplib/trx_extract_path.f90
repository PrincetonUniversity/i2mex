!
! RGA, parse the path string for info
!
! apath has the format for unix file  <path>/<runid>
!                      for mds        MDS+:<server-name>:<tree-name>(<shot-number>)
!                                     MDS+:<server-name>:<tree-name>(<tok.yy>,<runid>)
!
! returns as much information as possible
!
subroutine trx_extract_path(apath, ishot, aserver, atree, ashot, adir, adevice, arunid, ayear)
  implicit none

  character*(*), intent(in)  :: apath    ! parse this path
  integer,       intent(in)  :: ishot    ! shot number from namelist if available, used to guess runid

  character*(*), intent(out) :: aserver  ! mdsplus server
  character*(*), intent(out) :: atree    ! mdsplus server
  integer,       intent(out) :: ashot    ! mdsplus shot number or 0
  character*(*), intent(out) :: adir     ! directory name
  character*(*), intent(out) :: adevice  ! tokamak
  character*(*), intent(out) :: arunid   ! best guess at runid
  integer,       intent(out) :: ayear    ! year if known otherwise -1
  
  character*26,parameter :: alpha = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'  

  logical      :: ifixtree             ! .true. if a tree called transp should be called transp_<tok>
  character*5  :: a5
  character*16 :: rin                  ! temporary
  character*20 :: ccshot,tokyy         ! shot number and <tok>.<yy>
  integer      :: k,n,m
  integer      :: i1,i2,ic1,ic2,ix,iy  ! temp indices into apath
  
  aserver="" ; atree=""   ; ashot=0
  adevice="" ; arunid=""  ; ayear=-1
  adir=""
  
  i1=0 ; i2=0 ; ic1=0 ; ic2=0
  tokyy    = ""
  ifixtree = .false.
  
  ! look for MDS+:<stuff>:<stuff>(<stuff>)
  if (len(apath)>5) then
     a5 = apath(1:5)
     call ulower(a5)
     if (a5 == 'mds+:') then
        i1 = index(apath,':')
        if (i1>0 .and. len(apath)>i1) then
           i2 = index(apath(i1+1:),':')+i1
           if (i2>i1 .and. len(apath)>i2) then
              ic1 = index(apath(i2+1:),'(')+i2
              if (ic1>i2 .and. len(apath)>ic1+1) then
                 ic2 = index(apath(ic1+1:),')')+ic1
              end if
           end if
        end if
     end if
  end if
  
  if (i1>0 .and. i2>i1 .and. ic1>i2 .and. ic2>ic1) then
     !
     ! MDS+:<stuff>:<stuff>(<stuff>) format
     !
     aserver = adjustl(apath(i1+1:i2-1))
     atree   = adjustl(apath(i2+1:ic1-1))
     call ulower(atree)
     call ulower(aserver)

     ix = index(apath(ic1+1:ic2-1),',')+ic1  ! distinguish between (<shot-number>) and (<tok.yy>,<runid>)
     if (ix>ic1) then
        ! <tok.yy>,<runid>
        if (ix>ic1+1) tokyy  = adjustl(apath(ic1+1:ix-1))
        if (ic2>ix+1) arunid = adjustl(apath(ix+1:ic2-1))
        ifixtree = .true.
     else
        ! <shot-number>
        rin = adjustl(apath(ic1+1:ic2-1))
        read(rin,'(i10)',err=20) ashot
        goto 10
20      ashot=0
10      continue
        if (ashot>0 .and. ishot>0) then     ! guess runid from mds_shot and transp shot number
           k = ashot - ishot*10000          ! for 37065Z47,  2647 = 370652647 - 370650000
           n = k/100                        ! 'Z'
           m = mod(k,100)                   ! 47
           if (k>=100 .and. k<=2699 .and. n>=1 .and. n<=26 .and. m>=0) then   ! A00 to Z99
              write(ccshot,'(i9)') ishot
              write(a5,'(i2)') m
              ccshot = adjustl(ccshot)
              a5     = adjustl(a5)
              if (m>9) then
                 arunid = trim(ccshot)//alpha(n:n)//a5(1:2)
              else
                 arunid = trim(ccshot)//alpha(n:n)//'0'//a5(1:1)
              end if
           end if
        end if

        ! look for device if transpgrid server and tree transp_<tok>
        if (len_trim(aserver)>=10 .and. len_trim(atree)>7) then
           if (aserver(1:10)=="transpgrid" .and. atree(1:7)=="transp_") then
              adevice = atree(8:)
           end if
        end if
     end if
  else if (i1==0) then
     !
     ! unix path, might be in form
     ! <stuff>/<tok>.<yy>/<runid>
     !
     i1 = index(apath,'/',.true.)
     if (i1>0) then
        if (i1<len_trim(apath)) then
           arunid = trim(adjustl(apath(i1+1:)))
           call uupper(arunid)
           ix = index(arunid,'.CDF',.true.)     ! might be 37065Z45.CDF
           if (ix>1) arunid = arunid(:ix-1)
        end if
        i1 = max(1,i1-1)            ! in case of '/'
        adir = adjustl(apath(1:i1))
        i2 = index(adir,'/',.true.)
        if (i2>1 .and. i2<len_trim(adir)) tokyy = trim(adir(i2+1:))  ! might be tokyy
     end if
  end if
  
  !
  ! extract from <tok>.<yy>
  !
  if (len_trim(tokyy)>0) then
     iy = index(tokyy,'.')
     if (iy>1 .and. iy<len_trim(tokyy)) then
        adevice = trim(adjustl(tokyy(:iy-1)))
        rin = adjustl(tokyy(iy+1:))
        read(rin,'(i9)',err=40) ayear
        goto 30
40      ayear=-1
        adevice = ""
30      continue
        if (ayear>=0 .and. ayear<1900) then
           if (ayear<70) then
              ayear=2000+ayear
           else
              ayear=1900+ayear
           end if
        end if
        if (ayear<=1900 .or. ayear>2100) then
           ayear=-1 ; adevice=""
        end if
     end if
  end if

  call uupper(adevice)
  call uupper(arunid) 

  if (ifixtree .and. atree=="transp" .and. len_trim(adevice)>0) then  ! convert to standard mdsplus tree format
     atree = "transp_"//trim(adevice)
     call ulower(atree)
     if (len_trim(arunid)>0) then
        call mkbeastid(arunid,ashot,ccshot)
     end if
  end if

  !
  ! if can't find runid use ishot
  !
  if (len_trim(arunid)==0 .and. len(arunid)>=9) then
     write(arunid,'(i9)') ishot  ! punt
  end if
end subroutine trx_extract_path
