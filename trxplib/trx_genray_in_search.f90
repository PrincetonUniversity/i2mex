subroutine trx_genray_in_search(nonlin,ierr)

  ! look for genray template input file; write copy as "genray_ech.template"
  !
  ! if connected to run via MDSplus, try run cache data and then mds_get_inf
  ! if connected via files look for <runid>_genray_ech.in
  !
  ! The above will not succeed if <runid> did not use GENRAY for ECH.  In
  ! this case, copy the file from $LOCAL/tables/nml or $TRANSP_LOCATION...

  use trx_module

  implicit NONE

  integer, intent(in) :: nonlin ! I/O unit for messages
  integer, intent(out) :: ierr  ! completion status, 0=OK

  !-------------------------------
  character*40 :: ztarget
  character*150 :: genray_ech_template
  character*140 :: fdisk,fdir,mds_server,mds_tree
  character*200 :: arch_file,zfile
  character*20 :: runid,tfile
  logical :: ilmds

  integer :: ilun,kmdsplus,istat,imiss,iat,ilen,ier_cache,idum

  ! cache status routines

  logical mds_cache_active,mds_cache_exclusive
  external mds_cache_active,mds_cache_exclusive
  !-------------------------------

  ierr=0

  call find_io_unit(ilun)

  tfile = "genray_ech.template"

  call trread_dd(fdisk,fdir,runid,ilmds,ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' %trx_genray_in_search: run path lookup error.'
     return
  endif

  if(ilmds) then
     write(nonlin,*) ' ...look for GENRAY template in MDS+ run data... '
     imiss=0
     if(mds_cache_active(0)) then
        call mds_cache_fname(0,'GENRAY_ECH.DAT',zfile,ier_cache)
        if(ier_cache.ne.0) then
           if(mds_cache_exclusive(0)) then
              write(nonlin,*) ' ?trx_genray_in_search: namelist not in cache;'
              write(nonlin,*) '  RPLOT_CACHE_ONLY = TRUE -> no MDS+ recovery.'
              ierr=1
              return
           else
              write(nonlin,*) ' %trx_genray_in_search: cache miss, use MDS+.'
              imiss=1  ! normal cache miss
           endif
        else
           write(nonlin,*) ' '
           write(nonlin,*) ' (...attempting copy from MDS+ cache...) '
           call text_copy(trim(zfile),trim(tfile),istat)
           if(istat.ne.0) then
              if(mds_cache_exclusive(0)) then
                 write(nonlin,*) &
                      ' ?trx_genray_in_search: namelist not in cache;'
                 write(nonlin,*) &
                      '  RPLOT_CACHE_ONLY = TRUE -> no MDS+ recovery.'
                 ierr=1
                 return
              else
                 write(nonlin,*) &
                      ' %trx_genray_in_search: OK: cache miss, use MDS+.'
                 imiss=1  ! normal cache miss
              endif
           else
              write(nonlin,*) &
                   ' %trx_genray_in_search: MDS+ archived GENRAY template file was copied.'
              return
           endif
        endif
     endif

     write(nonlin,*) ' ...contact MDS+ server... '
     kmdsplus=2
     iat = index(fdisk,'@')
     ilen = len(trim(fdisk))
     if(iat.le.6) then
        ierr=1
        write(nonlin,*) ' ?trx_genray_in_search: MDS+ path data parse error.'
        write(nonlin,*) '   fdisk: ',trim(fdisk)
        write(nonlin,*) '    fdir: ',trim(fdir)
        return
     endif

     mds_server = fdisk(6:iat-1)
     mds_tree   = fdisk(iat+1:ilen)

  else
     kmdsplus=0
     mds_server=' '
     mds_tree=' '

     write(nonlin,*) ' ...look for GENRAY template namelist in run archive files... '
     arch_file = trim(fdir)//'/'//trim(runid)//'_genray_ech.in'
     open(unit=ilun,file=arch_file,status='old',action='READ',iostat=istat)
     if(istat.ne.0) then
        write(nonlin,*) ' %trx_genray_in_search: GENRAY template not found in archive files (OK).'
        write(nonlin,*) '  A '//trim(tdev)//' standard template will be sought.'
     else
        close(unit=ilun)
        call text_copy(trim(arch_file),trim(tfile),istat)
        if(istat.ne.0) then
           write(nonlin,*) ' ?trx_genray_in_search: archive file copy error.'
        else
           write(nonlin,*) ' %trx_genray_in_search: archived GENRAY template file was copied.'
           return
        endif
     endif
  endif
  
  call splitn_cget('GENRAY_ECH_TEMPLATE',len(genray_ech_template),1, &
       genray_ech_template,ierr)
  if(ierr.ne.0) then
     write(nonlin,*) ' %trx_genray_in_search: could not fetch GENRAY_ECH_TEMPLATE namelist variable.'
     return
  endif

  if(genray_ech_template.ne.' ') then
     write(nonlin,*) ' ...TRANSP namelist points at non-default GENRAY template ID: '//trim(genray_ech_template)
  endif

  call copy_gennml_sub(tdev,runid,genray_ech_template, &
       kmdsplus,mds_server,mds_tree, &
       tfile, &
       nonlin,ierr)

  ! insert in CACHE if indicated...

  if(ierr.eq.0) then
     if((kmdsplus.gt.0).AND.(imiss.eq.1)) then
        write(nonlin,*) ' GENRAY template from MDS+ server ... saved in local disk cache... '
        call text_copy(trim(tfile),trim(zfile),idum)
     endif
  endif

end subroutine trx_genray_in_search
