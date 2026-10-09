subroutine run_mg(mgname,rpitem,numruns,runPath,runAlias,inum_cur,iwarn,ier)
 
  ! create a multigraph of a single item (rpitem, e.g. "TE" or "PCUR")
  ! read from several different runs.
 
  ! create the member functions for the multigraph, forming labels of the
  ! form <rpitem>$<runAlias(i)> for i = 1 to numruns
 
  ! inum_cur.ne.0 indicates that this runPath corresponds to the "main"
  ! run to which the session is connected (via prior call to rconnect or
  ! kconnect), hence we just copy the data item and do not have to read
  ! from a separate run database (no trprofil or trscalar call needed)

  use datmgr_mod
  use cplotr_mod
 
  implicit NONE
 
  character*(*), intent(in) :: mgname   ! name of multigraph (mg)
  character*(*), intent(in) :: rpitem   ! item out of which to build mg
 
  ! the multigraph label and phys. units will be taken from rpitem in
  ! the current run (which must exist)
 
  integer, intent(in) :: numruns        ! no. of runs (.le. 15)
  character*(*), intent(in) :: runPath(numruns)  ! run paths (rconnect syntax)
  ! mod DMC June 2006: runPath(...) can have a suffix of the form:
  !    |<time-shift>, with <time-shift> encoded 1pe14.7, to specify a time
  !    shift for the corresponding run input data.

  character*(*), intent(in) ::  runAlias(numruns)  ! run aliases (max 10 chars)
 
  ! for description of runPath syntax see rconnect subroutine
  ! note each runAlias must be unique
 
  integer, intent(in) :: inum_cur       ! if non-zero -- index to current run
 
  integer, intent(out) :: iwarn(numruns)  ! warn on data access failures
  integer, intent(out) :: ier           ! completion code, 0=OK
 
  ! ier.ne.0 means an error with the arguments or that *all* data access
  ! attempts for *all* runs failed...
 
  !    ier=1  mg <mgname> already exists
  !    ier=2  item <rpitem> does *not* exist in current run
  !    ier=3  numruns.le.0  .or.  numruns.gt.15
  !    ier=4  inum_cur.lt.0  .or.  inum_cur.gt.numruns
  !    ier=5  alias error:
  !       runAlias(i).eq.runAlias(j) for some i.ne.j
  !    ier=6  **all** data access attempts failed
  !    ier=7  (... some other error ...)
 
  !   iwarn(j)=0 is the normal return
  !   iwarn(j)=1 -- data access failure on runPath(j)
  !   iwarn(j)=2 -- <rpitem>$<runAlias(j)> already exists -- reused.
 
  !---------------------------------------------
  !  SAVED local...
  integer, parameter :: max_history = 2000

  integer, save :: num_history = 0
  character*10, save :: alias_history(max_history)
  real, save :: tshift_history(max_history)

  !---------------------------------------------
  !  local...
 
  integer ifind_ordr,ifcn,ind,indr,indt,lunt,lunzer
  integer i,j,idup,il
 
  integer ilr,ili
 
  character*21 newids(max(1,min(15,numruns)))
  character*10 runids(max(1,min(15,numruns)))
  character*130 runpaths(max(1,min(15,numruns))),runpath1
  character*14 zdcod
  logical iok(max(1,min(15,numruns)))
  real zttags(max(1,min(15,numruns)))
  integer :: indbr
 
  integer ipt,iacc,itype,itype2
  integer icount
 
  character*64 zlabel
  !---------------------------------------------
 
  iwarn=0
  ier=0
 
  lunt=lunzer(0)
 
  !---------------------------------------------
  ! ---> basic error checks...
 
  ind=ifind_ordr(abb,iordrb,nbal,mgname)
  if(ind.gt.0) then
     ier=1
     write(lunt,*) '?run_mg: multigraph already defined: ',mgname
     write(lunt,*) ' calculator command available:  %MG_DELETE(<mg_name>)'
  endif
 
  indr=ifind_ordr(abr,iordrr,nfxt,mgname)
  indt=ifind_ordr(abt,iordrt,nft,mgname)
  if(max(indr,indt).gt.0) then
     ier=1
     write(lunt,*) &
          '?run_mg: multigraph name "',mgname,'" cannot be used.'
     write(lunt,*) &
          ' name is already used for a function id.'
  endif
 
  indr=ifind_ordr(abr,iordrr,nfxt,rpitem)
  indt=ifind_ordr(abt,iordrt,nft,rpitem)
  if((max(indr,indt).eq.0).or.(indr.gt.nfxt0).or.(indt.gt.nft0)) then
     ier=2
     write(lunt,*) &
          '?run_mg: function item not in current run database: ',rpitem
  endif
 
  if((numruns.le.0).or.(numruns.gt.15)) then
     ier=3
     write(lunt,*) '?run_mg:  invalid multigraph run count: ',numruns
  endif
 
  if((inum_cur.lt.0).or.(inum_cur.gt.15)) then
     ier=4
     write(lunt,*) '?run_mg:  current run identifier out of range: ',inum_cur
  endif
 
  if(ier.ne.0) return
 
  idup=0
  do i=1,numruns
     do j=1,numruns
        if(i.ne.j) then
           if(runAlias(i).eq.runAlias(j)) idup=i
        endif
     enddo
  enddo
 
  if(idup.gt.0) then
     write(lunt,*) '?run_mg:  run alias "', &
          runAlias(idup),'" occurs more than once.'
     ier=5
  endif
 
  if(ier.ne.0) return
 
  !---------------------------------------------
  ! --> get/check new function ids; separate path and runid information
 
  ilr=len_trim(rpitem)
 
  indr=ifind_ordr(abr,iordrr,nfxt,rpitem)
  indt=ifind_ordr(abt,iordrt,nft,rpitem)
 
  do i=1,numruns
 
     if(i.ne.inum_cur) then
        runpath1 = runpath(i)
        zttags(i) = 0.0
        indbr = index(runpath1,'|')
        if(indbr.gt.0) then
           zdcod = runpath1(indbr+1:)
           read(zdcod,'(1pe14.7)') zttags(i)
           runpath1(indbr:)=' '
        endif
        call run_path_sep(runPath1,runpaths(i),runids(i),ier)
        if(ier.ne.0) then
           ier=7
           il=len_trim(runPath(i))
           write(lunt,*) ' ?run_mg:  run path syntax error: ',runPath(i)(1:il)
           return
        endif
     else
        runpaths(i)=' '
        runids(i)=runid
        zttags(i) = 0.0
     endif

     ! following call deletes old copies of data, if time shift has changed...
     call history_check(runAlias(i),zttags(i),ier)
     if(ier.ne.0) return
 
     newids(i)=rpitem(1:ilr)//'$'//runAlias(i)
 
     ind=ifind_ordr(abr,iordrr,nfxt,newids(i))
     if(ind.gt.0) iwarn(i)=2
     ind=ifind_ordr(abt,iordrt,nft,newids(i))
     if(ind.gt.0) iwarn(i)=2
     if(iwarn(i).eq.2) then
        write(lunt,*) &
             '%run_mg:  function will be reused:  ',newids(i)
     endif

  enddo
 
  !---------------------------------------------
  ! --> if copying data from current run, do it now
 
  iok = .FALSE.
  icount=0
 
  if(indt.gt.0) zlabel = labelt(indt)
  if(indr.gt.0) zlabel = labelr(indr)
 
  if(inum_cur.gt.0) then
 
     iok(inum_cur)=.TRUE.
     icount=1
 
     if(indt.gt.0) then
 
        call dmgfotx(2,ipt,ier)
        if(ier.ne.0) call errmsg_exit(' ?? run_mg: internal error [2].')
 
        if (iwarn(inum_cur).eq.0) then
           iacc=ipt+(indt-1)*ntt
 
           call plftmk4(time,datbuf(iacc),ntt,labelt(indt),unitst(indt), &
                newids(inum_cur))
        endif
 
        itype=-1
 
     else if(indr.gt.0) then
 
        call dmgfxt(indr,ind)
        ipt=locd(ind)
        itype=itypr(indr)
 
        if(iwarn(inum_cur).eq.0) then
           ier=-88  ! suppress id check; accept longer name w/"$"
           call plsfsave(newids(inum_cur),labelr(indr),unitsr(indr),itype, &
                ipt,ier)
           if(ier.ne.0) call errmsg_exit(' ?? run_mg: internal error [3].')
        endif
 
     else if(iwarn(inum_cur).eq.0) then
        call errmsg_exit(' ?? run_mg: internal error [4].')
     endif
 
  endif
 
  !---------------------------------------------
  ! --> OK, retrieve data from other runs...
 
  do i=1,numruns
 
     if(i.eq.inum_cur) cycle
     if(iwarn(i).eq.2) then
        iok(i)=.TRUE.
        icount=icount+1
        cycle
     endif
 
     ier=-88  ! suppress id check; accept longer name w/"$"
     call plcfget(runpaths(i),runids(i),rpitem,newids(i),zttags(i),itype2,ier)
     if(ier.ne.0) then
        write(lunt,*) &
             ' %run_mg:  could not read "',rpitem,'" from ',runAlias(i)
        iwarn(i)=1
     else if(itype2.ne.itype) then
        write(lunt,*) &
             ' %run_mg:  "',rpitem,'" in ',runAlias(i),' type inconsistency.'
        iwarn(i)=1
     else
        iok(i)=.TRUE.
        icount=icount+1
     endif
 
  enddo
 
  if(icount.eq.0) then
     write(lunt,*) ' ?run_mg:  all reads failed, "',mgname,'" not created.'
     ier=6
     return
  endif
 
  !---------------------------------------------
  ! --> OK, build the multigraphs via calculator interface
 
  call run_mg_make(mgname,zlabel,numruns,newids(1:numruns),iok(1:numruns), &
       lunt,ier)

  return
 
CONTAINS
 
  subroutine run_path_sep(rp,rproot,rptail,ier)
      
    !  extract runid from full path specification
    !  if MDS+, just copy rp -> rproot
    !  if file, copy directory path -> rproot
 
    character*(*), intent(in) :: rp
    character*(*), intent(out) :: rproot  ! path
    character*(*), intent(out) :: rptail  ! runid
    integer, intent(out) :: ier  ! completion code: 0=OK
 
    character*4 ztest
 
    integer ilparen,irparen,icomma,ic,il,ibrk
    character*1 c
 
    !----------------------------
 
    ier=0
 
    ztest=rp(1:4)
    call uupper(ztest)
    if(ztest.eq.'MDS+') then
 
       !  MDS+ path spec
 
       rproot = rp
       ilparen=index(rp,'(')
       irparen=index(rp,')')
       icomma=index(rp,',')
       if(icomma.eq.0) then
          rptail=rp(ilparen+1:irparen-1)
       else
          rptail=rp(icomma+1:irparen-1)
       endif
 
    else
 
       !  file spec
 
       il=len_trim(rp)
       ibrk=0
       do ic=il,1,-1
          c=rp(ic:ic)
          if((c.eq.'/').or.(c.eq.']').or.(c.eq.'>')) then
             ibrk=ic
             exit
          endif
       enddo
 
       rproot=' '
       if(ibrk.gt.0) rproot=rp(1:ibrk)
       rptail=rp(ibrk+1:il)
 
    endif
 
  end subroutine run_path_sep

  subroutine history_check(alias,tshift,ier)

    ! check history of alias; detect time shift change
    ! if time shift change is detected, old copies of data are deleted!

    character*(*), intent(in) :: alias
    real, intent(in) :: tshift
    integer, intent(out) :: ier

    !-------------------------------
    integer :: ii,imatch
    !-------------------------------

    ier = 0

    imatch = 0
    do ii=1,num_history
       if(alias.eq.alias_history(ii)) then
          imatch = ii
          exit
       endif
    enddo

    if(imatch.eq.0) then
       ! first time for this alias; add to list if possible
       if(num_history.eq.max_history) then
          write(lunt,*) ' ?run_mg: maximum alias history count exceeded!'
          write(lunt,*) '  Need to increase parameter "max_history" & rebuild.'
          ier=1
       else
          num_history = num_history + 1
          alias_history(num_history) = alias
          tshift_history(num_history) = tshift
       endif

       return
    endif

    if(tshift.ne.tshift_history(imatch)) then
       write(lunt,*) &
            ' %run_mg: time shift value changes for run Alias: '//alias
       write(lunt,*) '  (old time shift: ',tshift_history(imatch),')'
       write(lunt,*) '  (new time shift: ',tshift,')'
       write(lunt,*) '  ...removing old data with incorrect shift.'
       tshift_history(imatch)=tshift
       call cleanup(alias)
    endif

  end subroutine history_check

  subroutine cleanup(alias)
    character*(*), intent(in) :: alias

    character*20 atest
    integer :: ii,idollr,idum

    ii=nft+1
    do 
       ii=ii-1
       if(ii.le.0) exit

       idollr = index(abt(ii),'$')
       if(idollr.le.0) cycle

       atest = abt(ii)(idollr+1:)
       if(atest.eq.alias) then
          call plcdelfn(abt(ii),idum)
       endif
    enddo

    ii=nfxt+1
    do
       ii=ii-1
       if(ii.le.0) exit

       idollr = index(abr(ii),'$')
       if(idollr.le.0) cycle

       atest = abr(ii)(idollr+1:)
       if(atest.eq.alias) then
          call plcdelfn(abr(ii),idum)
       endif
    enddo

  end subroutine cleanup

end subroutine run_mg

subroutine run_mg_make(mgname,label,numruns,newids,iok,lunt,ier)

  !  create multigraph from passed list of data items

  implicit NONE

  character*(*), intent(in) :: mgname   ! name of MG
  character*(*), intent(in) :: label    ! label for MG (units from data items)

  integer, intent(in) :: numruns        ! number of data items (=#runs)
  character*(*), intent(in) :: newids(numruns)   ! data item names
  logical, intent(in) :: iok(numruns)   ! TRUE to include item

  integer, intent(in) :: lunt           ! LUN for messages
  integer, intent(out) :: ier           ! status code on exit, 0=OK

  !------------------------
  character*500 calc_line
  integer :: ilm,ilc,il,i
  integer icount,ict,iclim,iprev,iwarn2
  !------------------------

  ilm=len_trim(mgname)
  calc_line='%MG_CREATE('//mgname(1:ilm)//',"'//trim(label)//'"'
 
  ilc=len_trim(calc_line)
 
  ict=0

  icount=0
  do i=1,numruns
     if(iok(i)) icount=icount + 1
  enddo
 
  iclim=5
  do i=1,numruns
 
     if(iok(i)) then
        calc_line(ilc+1:)=','//newids(i)
        ilc=len_trim(calc_line)
        ict=ict+1
        iprev=i
        if(ict.eq.iclim) exit
     endif
 
  enddo
 
  calc_line(ilc+1:ilc+1)=')'
 
  call rpcal0(calc_line,iwarn2,ier)
  if(ier.ne.0) then
     write(lunt,*) ' ??run_mg:  multigraph "',mgname,'" creation failure.'
     ier=7
     return
  endif
 
  do while(ict.lt.icount)
 
     iclim=min(icount,ict+5)
 
     calc_line='%MG_ADDFUN('//mgname(1:ilm)
 
     ilc=len_trim(calc_line)
 
     do i=iprev+1,numruns
 
        if(iok(i)) then
           calc_line(ilc+1:)=','//newids(i)
           ilc=len_trim(calc_line)
           ict=ict+1
           iprev=i
           if(ict.eq.iclim) exit
        endif
 
     enddo
 
     calc_line(ilc+1:ilc+1)=')'
 
     call rpcal0(calc_line,iwarn2,ier)
     if(ier.ne.0) then
        write(lunt,*) ' ??run_mg:  multigraph "',mgname,'" completion failure.'
        ier=7
        return
     endif
 
  enddo
 
end subroutine run_mg_make
