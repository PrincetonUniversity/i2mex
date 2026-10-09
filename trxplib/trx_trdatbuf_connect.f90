subroutine trx_trdatbuf_connect(iwarn)
 
  ! connect to TRDATBUF (TRANSP input data) dataset, if available

  ! mod DMC May 2009 -- add PH.CDF data to MDS+ disk cache, if flags are set

  use trx_module
  implicit NONE

  integer, intent(out) :: iwarn

  !------------------------------------
  integer :: lunzer,icdf,imiss,ier_cache
  character*200 zfile

  logical mds_cache_active,mds_cache_exclusive
  external mds_cache_active,mds_cache_exclusive

  !------------------------------------

  iwarn=0
  if(d_data_avail) return

  if(.not.d_init_flag) then
     call trdatbuf_init(d)
     d_init_flag = .TRUE.
  else
     call trdatbuf_reInit(d)
  endif

  if(mds_arch_flag) then

     iwarn=1
     imiss=0
     if(mds_cache_active(0)) then
        ! try to read from MDS+ cache
        call mds_cache_fname(0,'PH.CDF',zfile,ier_cache)
        if(ier_cache.ne.0) then
           if(mds_cache_exclusive(0)) then
              write(lunzer(0),*) &
                   ' ?trx_trdatbuf_connect: cache software error.'
              write(lunzer(0),*) &
                   '  RPLOT_CACHE_ONLY = TRUE -> no trdatbuf data.'
              iwarn=1
              return
           else
              write(lunzer(0),*) &
                   ' %trx_trdatbuf_connect: cache software error; revert to MDS+'
           endif
        else
           ! filename OK; open then read...
           call trdatbuf_open(zfile,"R",icdf)
           if(icdf.eq.0) then
              if(mds_cache_exclusive(0)) then
                 write(lunzer(0),*) &
                      ' ?trx_trdatbuf_connect: trdatbuf cache file not found.'
                 write(lunzer(0),*) &
                      '  RPLOT_CACHE_ONLY = TRUE -> no trdatbuf data.'
                 iwarn=1
                 return
              else
                 write(lunzer(0),*) &
                      ' %trx_trdatbuf_connect: cache miss, revert to MDS+'
                 imiss=1
              endif
           else
              ! file opened OK; read...
              call trdatbuf_read(icdf,d,iwarn)
              call trdatbuf_close(icdf)
              if(iwarn.ne.0) then
                 if(mds_cache_exclusive(0)) then
                    write(lunzer(0),*) &
                         ' ?trx_trdatbuf_connect: cache file read error.'
                    write(lunzer(0),*) &
                         '  RPLOT_CACHE_ONLY = TRUE -> no trdatbuf data.'
                    iwarn=1
                    return
                 else
                    write(lunzer(0),*) &
                         ' %trx_trdatbuf_connect: cache read error, revert to MDS+'
                    imiss=1
                 endif
              else
                 d_data_avail = .TRUE.  ! cache data read OK; iwarn=0 now
              endif
           endif
        endif

     endif

     if(iwarn.ne.0) then
        ! try read from MDS+ tree

        call rd_trdatbuf_mds(d,iwarn)
        if(iwarn.ne.0) then
           write(lunzer(0),*) ' %trx_trdatbuf_connect: no experimental data in MDS+ tree'
        else
           d_data_avail = .TRUE.
           if(imiss.eq.1) then
              ! cache is active; attempt cache file write
              call trdatbuf_open(zfile,"W",icdf)
              if(icdf.eq.0) then
                 ier_cache=1
              else
                 call trdatbuf_write(icdf,d,ier_cache)
                 call trdatbuf_close(icdf)
              endif
              if(ier_cache.ne.0) then
                 write(lunzer(0),*) ' %trx_datbuf_connect: cache write failed.'
              endif
           endif
        endif
     endif

  else
     
     ! try read from file

     zfile = trim(file_path)//'PH.CDF'

     call trdatbuf_open(zfile,"r",icdf)
     call trdatbuf_read(icdf,d,iwarn)
     call trdatbuf_close(icdf)

     if(iwarn.ne.0) then
        write(lunzer(0),*) ' %trx_trdatbuf_connect: experimental data file not found.'
     else
        d_data_avail = .TRUE.
     endif

  endif

end subroutine trx_trdatbuf_connect
