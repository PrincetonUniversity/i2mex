subroutine splitn_update_check(t_prev_update,t_try,filename,istat,ierr)

  ! check for time to write update namelist.  If it is time, write it.

  use splitn_module
  implicit NONE

  !--------------------------------------

  real*8, intent(inout) :: t_prev_update   ! time of last update
  !  (this is reset to the namelist update time if a namelist is written).

  real*8, intent(in) :: t_try              ! time to check to:
  !  if the update time T satisfies t_prev_update < T <= t_try, write update

  character*(*), intent(in) :: filename    ! name of (path to) file to write

  integer, intent(out) :: istat            ! output status report:
  !  =0: no update written
  !  =1: update written, t_prev_update updated.
  !  =2: update written, t_prev_update updated; additional update time(s)
  !      exist that are <= t_try.

  integer, intent(out) :: ierr             ! error status output (0=OK)

  !--------------------------------------
  integer :: ii,jj,kk,iupdate,ilun
  integer :: inaml,ivar,ilnum,jnaml,ifound
  integer :: ilin_min,ilin_max,ilen,inb
  integer :: inum_base,inum_st,i1,i2
  character*1 :: zquot,char1

  integer, dimension(:), allocatable :: inum_naml
  integer, dimension(:,:), allocatable :: ilins_naml
  character*32 :: znam
  !--------------------------------------

  istat=0
  ierr=0

  do ii=1,nupdate
     if((t_prev_update.lt.tup(ii)).AND.(tup(ii).le.t_try)) then
        istat = 1
        iupdate = ii
        if(ii.lt.nupdate) then
           if(tup(ii+1).le.t_try) then
              istat = 2
           endif
        endif

        exit
     endif
  enddo

  if(istat.gt.0) then

     call find_io_unit(ilun)
     open(unit=ilun,file=filename,status='unknown',iostat=ierr)

     if(ierr.ne.0) then
        write(6,*) ' ?splitn_update_check: could not open: ',trim(filename)
     else
        t_prev_update = tup(iupdate)
        write(6,*) ' %splitn_update_check: writing update: ',trim(filename), &
             ' @ t=',t_prev_update

        ! find range of update lines

        ilin_min=curmax + 1
        ilin_max=0

        inum_base = nreal+nint+nlog+nr8+nchv
        inum_st   = nreal_st+nint_st+nlog_st+nr8_st+nchv_st

        i1 = inum_base + (iupdate-1)*inum_st + 1
        i2 = inum_base + iupdate*inum_st

        do ii=i1,i2
           if(ilines(ii).ne.0) then
              ilin_min=min(ilin_min,ilines(ii))
              ilin_max=max(ilin_max,ilines(ii))
           endif
        enddo

        ! sort lines by namelist

        ilnum=ilin_max-ilin_min+1
        allocate(ilins_naml(ilnum,nnamls),inum_naml(nnamls))
        inum_naml = 0

        jnaml=0
        do ii=ilin_min,ilin_max
           jj=ordl(ii)
           if(lenl(jj).eq.0) cycle

           ! carry over lines without LHS...

           if(namfld(1,jj).eq.0) then
              if(jnaml.gt.0) then
                 kk=inum_naml(jnaml)+1
                 ilins_naml(kk,jnaml)=jj
                 inum_naml(jnaml)=kk
              endif
              cycle
           endif

           znam = textnl(jj)(namfld(1,jj):namfld(2,jj))
           if(znam(1:1).eq.'~') cycle  !skip ~UPDATE_TIME

           !  exclude possible element index reference -- want name only
           kk=index(znam,'(')
           if(kk.gt.0) znam(kk:)=' '

           call uupper(znam)

           ! find namelist to which name belongs; write warning if not found.
           ! (it should always be found).

           do inaml=1,nnamls
              ifound=0
              do kk=1,num_naml_vars_st(inaml)
                 ivar=indx_naml_vars_st(kk,inaml)
                 if(varlist(ivar)%name.eq.znam) then
                    ifound=ivar
                    exit
                 endif
              enddo
              if(ifound.gt.0) then
                 jnaml=inaml
                 kk=inum_naml(jnaml)+1
                 ilins_naml(kk,jnaml)=jj
                 inum_naml(jnaml)=kk
                 exit
              endif
           enddo

           if(ifound.eq.0) then
              write(6,*) ' ?splitn_update_check: sort to namelist failed:'
              write(6,*) '  ',trim(textnl(jj))
           endif
        enddo

        ! write each namelist...
        do inaml=1,nnamls
           if(inum_naml(inaml).eq.0) cycle
           !  OK, not empty...
           cur_naml = all_namls(inaml)
           write(ilun,*) ' '
           write(ilun,*) '&'//trim(cur_naml)

           !  screen out comments & blank lines
           do ii=1,inum_naml(inaml)
              jj=ilins_naml(ii,inaml)
              ilen=abs(lenl(jj))
              zquot=' '
              inb=0
              do kk=1,ilen
                 char1 = textnl(jj)(kk:kk)
                 if((char1.ne.' ').and.(char1.ne.char(9))) inb=kk
                 if(zquot.ne.' ') then
                    if(char1.eq.zquot) zquot=' '
                 else
                    if(char1.eq.'!') then
                       ilen=kk-1
                       if(textnl(jj)(1:ilen).eq.' ') ilen=0
                       exit
                    endif
                 endif
              enddo
              if(inb.eq.0) ilen=0
              if(ilen.eq.0) cycle
              write(ilun,'(A)') trim(textnl(jj)(1:ilen))
           enddo
           write(ilun,*) '/'
        enddo

        close(unit=ilun)

        if(istat.eq.2) write(6,*) ' %splitn: updates pending...'
     endif
  endif

end subroutine splitn_update_check
