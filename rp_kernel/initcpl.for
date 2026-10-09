      subroutine initcpl
C
C  init CPLOTR COMMON for RPLOT or RPLOT CALCULATOR
C
      use datmgr_mod            ! timebased moved here...
      use cplotr_mod
C
      character*10 ztest,ztest2
C
      data icall/0/
C-------------------------------------------------
C
C  only needs to be called once...
C
      if(icall.ne.0) return
      icall=1
C
C  these calls just link in block data objects...
C
      call cplset_exec
      call rpcaldat_exec
C
C  COMMON CONSTANTS -- array sizes etc.
      call initcf
C
      TIMLAB='TIME'
      TIMUNS='SECONDS'
C
C  clear the run label
C
      RUNLB2=' '
C
C  DEFAULT, NO CORRECTIONS TO TIMEBASE -- init timebase arrays
C
      call dmg_tinit  ! sets NTCORR & LTWRIT(...); assures allocation of
                      ! time vectors...
C
C  INITIAL DEFAULT-- ALL FCN NAMES WILL BE LISTED
      SSELEC='*'
      SSELEC_SAVE = SSELEC
      KSELEC=1   ! .ascii output only, now...
C
C  INITIAL DEFAULT = GET PLOT DATA FROM RMS DEFAULT DIRECTORY
C
      LFDISK=0
      LFDIR=0
C
      call initnws
C
C---------------------
C  TRANSP imbed stuff
C
      transp_imbed=.FALSE.
      transp_state=0
      transp_tinit=0
C
      NLTRANSP=.TRUE.
C
C---------------------
C  MDS+ cache options
C
      mds_cache=.false.
      mds_cache_only=.false.
      mds_cache_root='HOME'
C
      call get_environment_variable('RPLOT_CACHE',tmp_filn)
      if(tmp_filn.ne.' ') then
c
c  cache flag is set
c
         mds_cache=.true.
c
c  set cache-only flag, if additional environment variable is TRUE
c
         call get_environment_variable('RPLOT_CACHE_ONLY',ztest)
         call uupper(ztest)
         mds_cache_only = (ztest.eq.'TRUE')
c
c  get cache root directory
c
         ztest=tmp_filn
         call uupper(ztest)
         write(lunzer(0),*) ' %initcpl:  MDS_CACHE enabled.'
         if((ztest.ne.'TRUE').and.(ztest.ne.'HOME').and.
     >      (ztest.ne.'~')) then
            mds_cache_root=tmp_filn
            ils=len_trim(mds_cache_root)
            write(lunzer(0),*) 'MDS_CACHE root directory:  ',
     >         mds_cache_root(1:ils)
         endif
c
c  get cache binary items limit
c
         nbcache_lim=500
         call get_environment_variable('RPLOT_CACHE_MAXB',ztest)
         if(ztest.ne.' ') then
            ilz=len_trim(ztest)
            ztest2=' '
            ztest2(10-ilz+1:10)=ztest(1:ilz)
            read(ztest2,'(i10)',iostat=istat) ilim
            if(istat.ne.0) then
               write(lunzer(0),*) ' ?initcpl:  RPLOT_CACHE_MAXB = "',
     >            ztest(1:ilz),'" integer decode failed.'
            else
               nbcache_lim=max(10,min(naxfxt,ilim))
               write(lunzer(0),*) ' %initcpl:  cache item limit set',
     >              ' RPLOT_CACHE_MAXB = ',nbcache_lim
            endif
         endif
      endif
C
C---------------------
C  MKS units conversion (dmc March 00)
C
      if(nu_mks.eq.0) then
         call u_mks_add('                ','                ',1.0)
         call u_mks_add('EV              ','KeV             ',1.0e-3)
         call u_mks_add('EV/CM           ','KeV/m           ',1.0e-1)
         call u_mks_add('EV/SEC          ','KeV/s           ',1.0e-3)
         call u_mks_add('N/CM3           ','#/m^3           ',1.0e6)

         call u_mks_add('N/CM**3         ','#/m^3           ',1.0e6)
         call u_mks_add('N/CM4           ','#/m^4           ',1.0e8)
         call u_mks_add('N/CM**4         ','#/m^4           ',1.0e8)
         call u_mks_add('N               ','#               ',1.0)
         call u_mks_add('WATTS/CM3       ','W/m^3           ',1.0e6)  ! 10

         call u_mks_add('W/CM**3         ','W/m^3           ',1.0e6)
         call u_mks_add('WATTS/CM2       ','W/m^2           ',1.0e4)
         call u_mks_add('WATTS           ','W               ',1.0)
         call u_mks_add('N/CM3/SEC       ','#/m^3/sec       ',1.0e6)
         call u_mks_add('N/CM2/SEC       ','#/m^2/sec       ',1.0e4)

         call u_mks_add('N/SEC           ','#/sec           ',1.0)
         call u_mks_add('JLES/CM3        ','J/m^3           ',1.0e6)
         call u_mks_add('JLES/CM4        ','J/m^4           ',1.0e8)
         call u_mks_add('JOULES          ','J               ',1.0)
         call u_mks_add('GRAMS/CM3       ','kg/m^3          ',1.0e3)  ! 20

         call u_mks_add('GRAMS           ','kg              ',1.0e-3)
         call u_mks_add('G/CM3/SEC       ','kg/m^3/s        ',1.0e3)
         call u_mks_add('Nt-M/CM3        ','Nt/m^2          ',1.0e6)
         call u_mks_add('Nt-M/CM2        ','Nt/m            ',1.0e4)
         call u_mks_add('Nt-M            ','Nt*m            ',1.0)

         call u_mks_add('NtM-S/CM3       ','Nt*s/m^2        ',1.0e6)
         call u_mks_add('Nt-M-SEC        ','Nt*m*s          ',1.0)
         call u_mks_add('NtMS2/CM3       ','Nt*s^2/m^2      ',1.0e6)
         call u_mks_add('NT-M-SEC2       ','Nt*m*s^2        ',1.0)
         call u_mks_add('NEWTONS/CM      ','Nt/m            ',1.0e2)  ! 30

         call u_mks_add('AMPS/CM3        ','A/m^3           ',1.0e6)
         call u_mks_add('AMPS/CM2        ','A/m^2           ',1.0e4)
         call u_mks_add('AMPS/CM         ','A/m             ',1.0e2)
         call u_mks_add('AMPS            ','A               ',1.0)
         call u_mks_add('TESLA*CM        ','T*m             ',1.0e-2)

         call u_mks_add('TESLA*CM2       ','T*m^2           ',1.0e-4)
         call u_mks_add('TESLA/CM        ','T/m             ',1.0e2)
         call u_mks_add('TESLA           ','T               ',1.0)
         call u_mks_add('TESLA/SEC       ','T/s             ',1.0)
         call u_mks_add('WEBERS          ','Wb              ',1.0)    ! 40

         call u_mks_add('WEBERS/RAD      ','Wb/rad          ',1.0)
         call u_mks_add('CM              ','m               ',1.0e-2)
         call u_mks_add('CM/SEC          ','m/sec           ',1.0e-2)
         call u_mks_add('/CM             ','/m              ',1.0e2)
         call u_mks_add('CM**-1          ','/m              ',1.0e2)

         call u_mks_add('CM**-2          ','/m^2            ',1.0e4)
         call u_mks_add('CM**-3          ','/m^3            ',1.0e6)
         call u_mks_add('CM**-4          ','/m^4            ',1.0e8)
         call u_mks_add('CM**2           ','m^2             ',1.0e-4)
         call u_mks_add('CM**3           ','m^3             ',1.0e-6) ! 50

         call u_mks_add('CM**2/SEC       ','m^2/s           ',1.0e-4)
         call u_mks_add('PASCALS         ','Pa              ',1.0)
         call u_mks_add('PASCALS/CM      ','Pa/m            ',1.0e2)
         call u_mks_add('PASCLS/CM       ','Pa/m            ',1.0e2)
         call u_mks_add('PASCALS/SEC     ','Pa/s            ',1.0)

         call u_mks_add('PAS/SEC         ','Pa/s            ',1.0)
         call u_mks_add('VOLTS           ','V               ',1.0)
         call u_mks_add('VOLTS/CM        ','V/m             ',1.0e2)
         call u_mks_add('SEC**-1         ','/s              ',1.0)
         call u_mks_add('1/SEC           ','/s              ',1.0)    ! 60

         call u_mks_add('RAD/SEC         ','rad/s           ',1.0)
         call u_mks_add('SEC**-2         ','/s^2            ',1.0)
         call u_mks_add('/SEC/CM         ','/s/m            ',1.0e2)
         call u_mks_add('1/SEC/CM        ','/s/m            ',1.0e2)
         call u_mks_add('SECONDS         ','s               ',1.0)

         call u_mks_add('SECS            ','s               ',1.0)
         call u_mks_add('GHz             ','Hz              ',1.0e9)
         call u_mks_add('MHz             ','Hz              ',1.0e6)
         call u_mks_add('RADIANS         ','rad             ',1.0)
         call u_mks_add('AMP*TESLA/CM2   ','A*T/m^2         ',1.0e4)  ! 70

         call u_mks_add('VOLT*TESLA/CM   ','V*T/m           ',1.0e2)
         call u_mks_add('OHMS            ','ohm             ',1.0)
         call u_mks_add('OHM*CM          ','ohm*m           ',1.0e-2)
         call u_mks_add('M**-2           ','m^2             ',1.0)
         call u_mks_add('WATTS/CM3/EV    ','W/m^3/KeV       ',1.0e9)

c
c  items in "degrees" to be left as such...
c
         call u_mks_add('Nt-M/CM3/(RAD/S)','Nt/m^2/(rad/s)  ',1.0e6)
         call u_mks_add('0=NORMAL        ','0=NORMAL        ',1.0)
         call u_mks_add('HOURS           ','s               ',3600.0)
         call u_mks_add('ARB.UNITS       ','Arb. Units      ',1.0)
         call u_mks_add('DEGREES         ','deg             ',1.0)
c
      endif
C
c------------------------------
c  MDSplus stuff...
c
      mds_save_server=' '
      mds_nclist=0
      mds_tree_clist=' '
      mds_shot_clist=0
c
c------------------------------
      return
      end
