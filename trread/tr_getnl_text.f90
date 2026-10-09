subroutine tr_getnl_text(ilun,ierr)
!
  use tr_getnl
  implicit NONE
!
  integer ilun                ! lun for file read (input)
  integer ierr                ! completion code:  0=OK
!
!    read the TRANSP namelist
!
  character*196 zpath,zfile,zfile2
  character*10 zrunid
!
  integer ilr,ilz,iline,ier_cache,imiss
!
! MDSplus Fortran interface
!
  integer Mds_Value,idescr_long,idescr_cstringarr,idescr_cstring,isize
!
! cache status routine
!
  logical mds_cache_active,mds_cache_exclusive
  external mds_cache_active,mds_cache_exclusive
!
! for MDS+ data access...
!
  integer status
  integer idims(2)
!
  integer lunzer,lt
  integer ialloc
!---------------------------------
!
  lt=lunzer(0)
!
  ierr=0
  if(nltext_status.ge.0) then
     ierr=nltext_status
     if(ierr.eq.99) then
        write(lt,*) ' ?tr_getnl_text: code error, tr_getnl_clear never called.'
     endif
     return
  endif
!
  call tgetpath(zpath,zrunid)
!
  if(zpath.eq.'MDS+') then
!
!  MDS+ access ...
!
!  try for cached file first...
!
     imiss=0
     if(mds_cache_active(0)) then
        call mds_cache_fname(0,'NAMELIST_TR.DAT',zfile,ier_cache)
        if(ier_cache.ne.0) then
           if(mds_cache_exclusive(0)) then
              write(lt,*) ' ?trx_get_nltext: namelist not in cache;'
              write(lt,*) '  RPLOT_CACHE_ONLY = TRUE -> no MDS+ recovery.'
              ierr=1
              go to 10
           else
              write(lt,*) ' %trx_get_nltext:  cache error, trying MDS+ ...'
           endif
        else
           open(unit=ilun,file=zfile,status='old',iostat=ier_cache)
           if(ier_cache.ne.0) then
              if(mds_cache_exclusive(0)) then
                 write(lt,*) ' ?trx_get_nltext: namelist not in cache;'
                 write(lt,*) '  RPLOT_CACHE_ONLY = TRUE -> no MDS+ recovery.'
                 ierr=1
                 go to 10
              else
                 imiss=1      ! normal cache miss
              endif
           else
!  use the cached namelist file -- read as ordinary text file.
              nltext_nlines=0
              nltext_status=0
              go to 10
           endif
        endif
     endif
!
!  MDS+ read.  if imiss=1, also write cache file
!-----------------------------------
!
     status = Mds_Value('size(NAME_LIST)',idescr_long(nltext_nlines),isize)
     if (mod(status,2).ne.1) then
        call mdserr(lt, 'Mds_Value(''size(NAME_LIST)'',..)', status)
        ierr=1
        nltext_nlines=0
        nltext_status=1
        return
     else if(nltext_nlines.le.0) then
        write(lt,*) ' ?Mds_Value returned non-positive namelist length:  ', &
             nltext_nlines
        ierr=1
        nltext_nlines=0
        nltext_status=1
        return
     else
!  OK
        if(allocated(nltext)) deallocate(nltext)
        if(allocated(nltext_lens)) deallocate(nltext_lens)
        allocate(nltext(nltext_nlines), STAT=ialloc)
        if(ialloc  /=  0) then
           call errmsg_exit('?tr_getnl_text: ERROR allocating nltext')
        endif
        allocate(nltext_lens(nltext_nlines), STAT=ialloc ); nltext_lens=0
        if(ialloc  /=  0) then
           print*,'nltext_nlines =',nltext_nlines
           call errmsg_exit('?tr_getnl_text: ERROR allocating nltext_lens')
        endif
 
        if(allocated(mltext)) deallocate(mltext)
        if(allocated(mltext_lens)) deallocate(mltext_lens)
        allocate(mltext(nltext_nlines+50))
        allocate(mltext_lens(nltext_nlines+50)); mltext_lens=0
!
        idims(1)=nltext_nlines
        idims(2)=0
!
        status = Mds_Value('NAME_LIST',idescr_cstringarr(NLTEXT,idims,1),isize)
        if (mod(status,2).ne.1) then
           call mdserr(lt, 'Mds_Value(''size(NAME_LIST)'',..)', status)
           ierr=1
           nltext_nlines=0
           nltext_status=1
           return
        else
           do iline=1,nltext_nlines
              call tr_getnl_cleanup(nltext(iline))
              nltext_lens(iline)=len_trim(nltext(iline))
              mltext(iline)=nltext(iline)
           enddo
        endif
     endif
!
!  OK write cache file
!
     if(imiss.eq.1) then
!
        open(unit=ilun,file=zfile,status='replace', &
             iostat=ier_cache)
!
        if(ier_cache.eq.0) then
           do iline=1,nltext_nlines
              ilz=nltext_lens(iline)
              if(ilz.eq.0) ilz=len_trim(nltext(iline))
              write(ilun,'(A)') nltext(iline)(1:ilz)
           enddo
           close(unit=ilun)
        else
           write(lt,*) ' %trx_get_nltext:  open failure, no cache write.'
        endif
!
     endif
!
     nltext_status=0        ! read was successful
     mltext_nlines=nltext_nlines
!
     return
!
!-----------------------------------------
!  non-MDS+
  else
     call tr_getnl_ftext(zpath,zrunid,-ilun,ierr)
  endif
!
!  read text file (normal namelist or cache file)
!  get length of file
 
10 continue
 
  if(ierr.eq.0) then
     call tr_getnl_lines(ilun,ierr)
  endif
 
  return
 
end

