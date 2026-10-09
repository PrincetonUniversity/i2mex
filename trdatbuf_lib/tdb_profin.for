      subroutine tdb_profin(d, t, ierr)

      !  fetch profile data from trdatbuf buffer (d) according to
      !  specifications in target structure (t).

      !  time averaging t%delta_t > 0 supported...

      use trdatbuf_obj
      use trdatbuf_aux
      use trdatbuf_iface, only: tdb_range_ntimes2,tdb_range_times2,
     >   tdb_rmp_bdy

      implicit NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

      type (trdatbuf) :: d
      type (profget) :: t
      integer, intent(out) :: ierr  ! completion code (0=OK)

      !---------------------------------------
      ! local information: for time averaging...

      integer :: int_avg        ! no. of time pts
      real*8 :: zt1,zt2         ! averaging time range
      real*8, dimension(:), allocatable :: zt_avg,zw_avg ! times & weights

      type(profget) :: t_use,t_sum
      integer :: it
      integer :: nonlin,lunmsg_tdb
      !---------------------------------------

      if(t%delta_t.le.0.0_R8) then
         call tdb_profin1(d, t, ierr) ! no time averaging
      else
         ! *** time averaging needed ***

         nonlin = lunmsg_tdb(0)

         call tdb_profget_init(t_use,t%nzones,t%idebug)

         t_use%item = t%item
         t_use%xibdys = t%xibdys
         t_use%plflxg = t%plflxg
         t_use%rmajmp = t%rmajmp
         t_use%bmidp = t%bmidp

         ! t_use%time will be set below...

         t_use%delta_t = 0
         t_use%ibdy = t%ibdy
         call tdb_profget_init(t_sum,t%nzones,t%idebug) ! sums are zeroed now.

         zt1 = t%time - t%delta_t
         zt2 = t%time + t%delta_t

                                ! fetch timebase for averaging...

         call tdb_range_ntimes2(d,zt1,zt2,int_avg) ! #time pts in database

         int_avg = int_avg+2    ! add in interval end points zt1,zt2

         allocate(zt_avg(int_avg),zw_avg(int_avg))

         zt_avg(1) = zt1
         zt_avg(int_avg) = zt2
         call tdb_range_times2(d,zt1,zt2,zt_avg(2:int_avg-1),ierr)

         if(ierr.ne.0) then
            write(nonlin,*)
     >           ' ?? tdb_profin with time averaging: unexpected error!'
         else

                                ! compute weighting factors proportional to dt

            zw_avg(1)=(zt_avg(2)-zt_avg(1))/2/(zt2-zt1)
            zw_avg(int_avg)=
     >           (zt_avg(int_avg)-zt_avg(int_avg-1))/2/(zt2-zt1)
            do it=2,int_avg-1
               zw_avg(it)=(zt_avg(it+1)-zt_avg(it-1))/2/(zt2-zt1)
            enddo

                                ! compute sum...

            do it=1,int_avg
               t_use%time = zt_avg(it)

               call tdb_profin1(d, t_use, ierr)
               if(ierr.ne.0) exit

               t_sum%data_zc = t_sum%data_zc + zw_avg(it)*t_use%data_zc
               t_sum%data_zb = t_sum%data_zb + zw_avg(it)*t_use%data_zb

               if(t%idebug) then
                  t_sum%ecegap = t_sum%ecegap+zw_avg(it)*t_use%ecegap
                  t_sum%asym_zc = t_sum%asym_zc+zw_avg(it)*t_use%asym_zc
                  t_sum%asym_zb = t_sum%asym_zb+zw_avg(it)*t_use%asym_zb
                  t_sum%shift = t_sum%shift+zw_avg(it)*t_use%shift
                  t_sum%datrsym = t_sum%datrsym+zw_avg(it)*t_use%datrsym
                  t_sum%datusym = t_sum%datusym+zw_avg(it)*t_use%datusym
                  t_sum%rmjsym = t_sum%rmjsym+zw_avg(it)*t_use%rmjsym
               endif
            enddo
                                ! all done; transfer output

            t%data_zc = t_sum%data_zc
            t%data_zb = t_sum%data_zb

            if(t%idebug) then
               t%ecegap = t_sum%ecegap
               t%asym_zc = t_sum%asym_zc
               t%asym_zb = t_sum%asym_zb
               t%shift = t_sum%shift
               t%datrsym = t_sum%datrsym
               t%datusym = t_sum%datusym
               t%xirsym = t_use%xirsym ! does not vary in time
               t%rmjsym = t_sum%rmjsym
            endif
         endif

         ! release storage

         call tdb_profget_free(t_use)
         call tdb_profget_free(t_sum)

         deallocate(zt_avg,zw_avg)

      endif

      end
!-------------------------------------------------------------------------
      subroutine tdb_profin1(d, t, ierr)

      !  fetch profile data from trdatbuf buffer (d) according to
      !  specifications in target structure (t).

      !  *** no time averaging *** t%delta_t = ZERO expected

      use trdatbuf_obj
      use trdatbuf_aux
      implicit NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

      type (trdatbuf) :: d
      type (profget) :: t
      integer, intent(out) :: ierr  ! completion code (0=OK)

      !  ierr is set if the profile id is unrecognized or if it is recognized
      !  but there is no such data...

      integer :: i,j,iadr,i2zp1
      character*3 :: zid

      integer :: nonlin,lunmsg_tdb

      real*8 zrmj_old,zrmj_new,zr1,zr2,zr1p,zr2p,zdelr,zf

      integer,save :: iwarn=5

      !----------------------------------

      ierr=1

      i=t%iselect
      zid = t%item
      call uupper(zid)
C
      nonlin = lunmsg_tdb(0)
C
C  tdb_profin1 check:
C
      if(t%delta_t.gt.0.0_R8) then
         write(nonlin,*) ' ? tdb_profin1 -- t%delta_t = ',t%delta_t
         write(nonlin,*)
     >        '   this routine does not support time averaging.'
         return
      endif
C
C  check that initialization was done...
C
      if(d%lbx(1).eq.0) then
         ierr=100
         write(nonlin,*)
     >        ' ? tdb_profin: error exit -- "d%lbx" unitialized'
         write(nonlin,*)
     >        ' ? call to "tdb_symini" is needed.'
         return
      endif

C
C  check normalization of R grid -- adjust, as may be necessary when
C  a free boundary code seeks data from a fixed boundary TRANSP dataset...
C
      if(t%nzones.le.0) then
         ierr=200
         write(nonlin,*)
     >        ' ? tdb_profin: t%nzones <= ZERO'
         return
      endif
C
      i2zp1=2*t%nzones + 1
      if((t%rmajmp(i2zp1)-t%rmajmp(1)).le.0.0_R8) then
         ierr=300
         write(nonlin,*)
     >        ' ? tdb_profin -- no midplane radii?  t%rmajmp = '
         write(nonlin,*) t%rmajmp(1:i2zp1)
         return
      endif
      if(min(t%rmajmp(i2zp1),t%rmajmp(1)).le.0.0_R8) then
         ierr=300
         write(nonlin,*)
     >        ' ? tdb_profin -- midplane radii <= 0?  t%rmajmp = '
         write(nonlin,*) t%rmajmp(1:i2zp1)
         return
      endif
C
      call tdb_rmp_bdy(d,t%time,zr1,zr2)
C
      zdelr=max(abs(zr1-t%rmajmp(1)),abs(zr2-t%rmajmp(i2zp1)))/
     >     max(zr1,zr2)
      if(zdelr.gt.0.01_R8) then
         if(iwarn.gt.0) then
            iwarn=iwarn-1
            write(nonlin,*) ' %tdb_profin:  t%rmajmp and input data',
     >           ' major radius grids not matched:'
            write(nonlin,*) ' R1(t%rmajmp, data): ',t%rmajmp(1),zr1
            write(nonlin,*) ' R2(t%rmajmp, data): ',t%rmajmp(i2zp1),zr2
            write(nonlin,*) '   (fixup applied).'
         endif
      endif
C
      zr1p=t%rmajmp(1)
      zr2p=t%rmajmp(i2zp1)
      do i=1,i2zp1
         zf = (t%rmajmp(i)-zr1p)/(zr2p-zr1p)
         zrmj_old=t%rmajmp(i)
         zrmj_new=zr1*(1.0_R8-zf)+zr2*zf
         t%rmajmp(i)=zrmj_new
         t%bmidp(i)=zrmj_old*t%bmidp(i)/zrmj_new
      enddo
      i=t%iselect    ! needed for SIM
C
C  if Te(R) data is requested -- provide automatic switching between TER and
C  ECF
C
      if(zid.eq.'TER') then
         if(d%lfter.eq.0) then
            if(d%lfecf.ne.0) then
               zid='ECF'
            endif
         endif
      else if(zid.eq.'ECF') then
         if(d%lfecf.eq.0) then
            if(d%lfter.ne.0) then
               zid='TER'
            endif
         endif
      endif
C
C  Similar equivalence for toroidal rotation data VP2 or OMG
C  (toroidal velocity or angular velocity input signals are all converted
C  to angular velocity in trdat).
C
      if(zid.eq.'OMG') then
         if(d%lfomg.eq.0) then
            if(d%lfvp2.ne.0) then
               zid='VP2'
            endif
         endif
      else if(zid.eq.'VP2') then
         if(d%lfvp2.eq.0) then
            if(d%lfomg.ne.0) then
               zid='OMG'
            endif
         endif
      endif
C
C=>TRDATGEN+
!
!    ******************************************
!    * TRDATGEN GENERATED CODE -- DO NOT EDIT *
!    ******************************************
!
!    code generation ends at C=>TRDATGEN- line
!
!
C
C------------------
      IF(ZID.EQ.'BOL') THEN
        t%IECE=.FALSE.
        IF(D%LFBOL.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXBOL.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIBOL,D%NSYBOL,
     >        D%LXBOL,D%NXBOL,D%LFBOL,
     >        D%LXSYBOL,D%NXSYBOL,D%LFSYBOL,D%LSSYBOL,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'BPA') THEN
        t%IECE=.FALSE.
        IF(D%LFBPA.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXBPA.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIBPA,D%NSYBPA,
     >        D%LXBPA,D%NXBPA,D%LFBPA,
     >        D%LXSYBPA,D%NXSYBPA,D%LFSYBPA,D%LSSYBPA,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'BPB') THEN
        t%IECE=.FALSE.
        IF(D%LFBPB.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXBPB.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIBPB,D%NSYBPB,
     >        D%LXBPB,D%NXBPB,D%LFBPB,
     >        D%LXSYBPB,D%NXSYBPB,D%LFSYBPB,D%LSSYBPB,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'D2F') THEN
        t%IECE=.FALSE.
        IF(D%LFD2F.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXD2F.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRID2F,D%NSYD2F,
     >        D%LXD2F,D%NXD2F,D%LFD2F,
     >        D%LXSYD2F,D%NXSYD2F,D%LFSYD2F,D%LSSYD2F,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'DE2') THEN
        t%IECE=.FALSE.
        IF(D%LFDE2.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXDE2.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIDE2,D%NSYDE2,
     >        D%LXDE2,D%NXDE2,D%LFDE2,
     >        D%LXSYDE2,D%NXSYDE2,D%LFSYDE2,D%LSSYDE2,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'DF3') THEN
        t%IECE=.FALSE.
        IF(D%LFDF3.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXDF3.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIDF3,D%NSYDF3,
     >        D%LXDF3,D%NXDF3,D%LFDF3,
     >        D%LXSYDF3,D%NXSYDF3,D%LFSYDF3,D%LSSYDF3,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'DF4') THEN
        t%IECE=.FALSE.
        IF(D%LFDF4.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXDF4.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIDF4,D%NSYDF4,
     >        D%LXDF4,D%NXDF4,D%LFDF4,
     >        D%LXSYDF4,D%NXSYDF4,D%LFSYDF4,D%LSSYDF4,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'DF6') THEN
        t%IECE=.FALSE.
        IF(D%LFDF6.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXDF6.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIDF6,D%NSYDF6,
     >        D%LXDF6,D%NXDF6,D%LFDF6,
     >        D%LXSYDF6,D%NXSYDF6,D%LFSYDF6,D%LSSYDF6,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'DFD') THEN
        t%IECE=.FALSE.
        IF(D%LFDFD.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXDFD.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIDFD,D%NSYDFD,
     >        D%LXDFD,D%NXDFD,D%LFDFD,
     >        D%LXSYDFD,D%NXSYDFD,D%LFSYDFD,D%LSSYDFD,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'DFH') THEN
        t%IECE=.FALSE.
        IF(D%LFDFH.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXDFH.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIDFH,D%NSYDFH,
     >        D%LXDFH,D%NXDFH,D%LFDFH,
     >        D%LXSYDFH,D%NXSYDFH,D%LFSYDFH,D%LSSYDFH,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'DFT') THEN
        t%IECE=.FALSE.
        IF(D%LFDFT.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXDFT.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIDFT,D%NSYDFT,
     >        D%LXDFT,D%NXDFT,D%LFDFT,
     >        D%LXSYDFT,D%NXSYDFT,D%LFSYDFT,D%LSSYDFT,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'ECF') THEN
        t%IECE=.TRUE.
        IF(D%LFECF.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXECF.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIECF,D%NSYECF,
     >        D%LXECF,D%NXECF,D%LFECF,
     >        D%LXSYECF,D%NXSYECF,D%LFSYECF,D%LSSYECF,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'GRB') THEN
        t%IECE=.FALSE.
        IF(D%LFGRB.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXGRB.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIGRB,D%NSYGRB,
     >        D%LXGRB,D%NXGRB,D%LFGRB,
     >        D%LXSYGRB,D%NXSYGRB,D%LFSYGRB,D%LSSYGRB,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'LF3') THEN
        t%IECE=.FALSE.
        IF(D%LFLF3.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXLF3.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRILF3,D%NSYLF3,
     >        D%LXLF3,D%NXLF3,D%LFLF3,
     >        D%LXSYLF3,D%NXSYLF3,D%LFSYLF3,D%LSSYLF3,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'LF4') THEN
        t%IECE=.FALSE.
        IF(D%LFLF4.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXLF4.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRILF4,D%NSYLF4,
     >        D%LXLF4,D%NXLF4,D%LFLF4,
     >        D%LXSYLF4,D%NXSYLF4,D%LFSYLF4,D%LSSYLF4,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'LF6') THEN
        t%IECE=.FALSE.
        IF(D%LFLF6.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXLF6.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRILF6,D%NSYLF6,
     >        D%LXLF6,D%NXLF6,D%LFLF6,
     >        D%LXSYLF6,D%NXSYLF6,D%LFSYLF6,D%LSSYLF6,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'LFD') THEN
        t%IECE=.FALSE.
        IF(D%LFLFD.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXLFD.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRILFD,D%NSYLFD,
     >        D%LXLFD,D%NXLFD,D%LFLFD,
     >        D%LXSYLFD,D%NXSYLFD,D%LFSYLFD,D%LSSYLFD,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'LFH') THEN
        t%IECE=.FALSE.
        IF(D%LFLFH.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXLFH.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRILFH,D%NSYLFH,
     >        D%LXLFH,D%NXLFH,D%LFLFH,
     >        D%LXSYLFH,D%NXSYLFH,D%LFSYLFH,D%LSSYLFH,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'LFT') THEN
        t%IECE=.FALSE.
        IF(D%LFLFT.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXLFT.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRILFT,D%NSYLFT,
     >        D%LXLFT,D%NXLFT,D%LFLFT,
     >        D%LXSYLFT,D%NXSYLFT,D%LFSYLFT,D%LSSYLFT,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NER') THEN
        t%IECE=.FALSE.
        IF(D%LFNER.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNER.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINER,D%NSYNER,
     >        D%LXNER,D%NXNER,D%LFNER,
     >        D%LXSYNER,D%NXSYNER,D%LFSYNER,D%LSSYNER,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NI3') THEN
        t%IECE=.FALSE.
        IF(D%LFNI3.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNI3.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINI3,D%NSYNI3,
     >        D%LXNI3,D%NXNI3,D%LFNI3,
     >        D%LXSYNI3,D%NXSYNI3,D%LFSYNI3,D%LSSYNI3,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NI4') THEN
        t%IECE=.FALSE.
        IF(D%LFNI4.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNI4.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINI4,D%NSYNI4,
     >        D%LXNI4,D%NXNI4,D%LFNI4,
     >        D%LXSYNI4,D%NXSYNI4,D%LFSYNI4,D%LSSYNI4,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NI6') THEN
        t%IECE=.FALSE.
        IF(D%LFNI6.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNI6.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINI6,D%NSYNI6,
     >        D%LXNI6,D%NXNI6,D%LFNI6,
     >        D%LXSYNI6,D%NXSYNI6,D%LFSYNI6,D%LSSYNI6,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NID') THEN
        t%IECE=.FALSE.
        IF(D%LFNID.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNID.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINID,D%NSYNID,
     >        D%LXNID,D%NXNID,D%LFNID,
     >        D%LXSYNID,D%NXSYNID,D%LFSYNID,D%LSSYNID,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NIH') THEN
        t%IECE=.FALSE.
        IF(D%LFNIH.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNIH.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINIH,D%NSYNIH,
     >        D%LXNIH,D%NXNIH,D%LFNIH,
     >        D%LXSYNIH,D%NXSYNIH,D%LFSYNIH,D%LSSYNIH,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NIM') THEN
        t%IECE=.FALSE.
        IF(D%LFNIM.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNIM.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINIM,D%NSYNIM,
     >        D%LXNIM,D%NXNIM,D%LFNIM,
     >        D%LXSYNIM,D%NXSYNIM,D%LFSYNIM,D%LSSYNIM,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NIT') THEN
        t%IECE=.FALSE.
        IF(D%LFNIT.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNIT.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINIT,D%NSYNIT,
     >        D%LXNIT,D%NXNIT,D%LFNIT,
     >        D%LXSYNIT,D%NXSYNIT,D%LFSYNIT,D%LSSYNIT,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'NMR') THEN
        t%IECE=.FALSE.
        IF(D%LFNMR.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXNMR.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRINMR,D%NSYNMR,
     >        D%LXNMR,D%NXNMR,D%LFNMR,
     >        D%LXSYNMR,D%NXSYNMR,D%LFSYNMR,D%LSSYNMR,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'OMG') THEN
        t%IECE=.FALSE.
        IF(D%LFOMG.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXOMG.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIOMG,D%NSYOMG,
     >        D%LXOMG,D%NXOMG,D%LFOMG,
     >        D%LXSYOMG,D%NXSYOMG,D%LFSYOMG,D%LSSYOMG,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'PRS') THEN
        t%IECE=.FALSE.
        IF(D%LFPRS.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXPRS.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIPRS,D%NSYPRS,
     >        D%LXPRS,D%NXPRS,D%LFPRS,
     >        D%LXSYPRS,D%NXSYPRS,D%LFSYPRS,D%LSSYPRS,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'QPR') THEN
        t%IECE=.FALSE.
        IF(D%LFQPR.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXQPR.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIQPR,D%NSYQPR,
     >        D%LXQPR,D%NXQPR,D%LFQPR,
     >        D%LXSYQPR,D%NXSYQPR,D%LFSYQPR,D%LSSYQPR,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'SBI') THEN
        t%IECE=.FALSE.
        IF(D%LFSBI.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXSBI.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRISBI,D%NSYSBI,
     >        D%LXSBI,D%NXSBI,D%LFSBI,
     >        D%LXSYSBI,D%NXSYSBI,D%LFSYSBI,D%LSSYSBI,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'TER') THEN
        t%IECE=.FALSE.
        IF(D%LFTER.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXTER.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRITER,D%NSYTER,
     >        D%LXTER,D%NXTER,D%LFTER,
     >        D%LXSYTER,D%NXSYTER,D%LFSYTER,D%LSSYTER,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'TI2') THEN
        t%IECE=.FALSE.
        IF(D%LFTI2.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXTI2.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRITI2,D%NSYTI2,
     >        D%LXTI2,D%NXTI2,D%LFTI2,
     >        D%LXSYTI2,D%NXSYTI2,D%LFSYTI2,D%LSSYTI2,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'TI3') THEN
        t%IECE=.FALSE.
        IF(D%LFTI3.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXTI3.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRITI3,D%NSYTI3,
     >        D%LXTI3,D%NXTI3,D%LFTI3,
     >        D%LXSYTI3,D%NXSYTI3,D%LFSYTI3,D%LSSYTI3,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'TQI') THEN
        t%IECE=.FALSE.
        IF(D%LFTQI.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXTQI.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRITQI,D%NSYTQI,
     >        D%LXTQI,D%NXTQI,D%LFTQI,
     >        D%LXSYTQI,D%NXSYTQI,D%LFSYTQI,D%LSSYTQI,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'V2F') THEN
        t%IECE=.FALSE.
        IF(D%LFV2F.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXV2F.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIV2F,D%NSYV2F,
     >        D%LXV2F,D%NXV2F,D%LFV2F,
     >        D%LXSYV2F,D%NXSYV2F,D%LFSYV2F,D%LSSYV2F,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VB2') THEN
        t%IECE=.FALSE.
        IF(D%LFVB2.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVB2.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVB2,D%NSYVB2,
     >        D%LXVB2,D%NXVB2,D%LFVB2,
     >        D%LXSYVB2,D%NXSYVB2,D%LFSYVB2,D%LSSYVB2,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VC3') THEN
        t%IECE=.FALSE.
        IF(D%LFVC3.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVC3.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVC3,D%NSYVC3,
     >        D%LXVC3,D%NXVC3,D%LFVC3,
     >        D%LXSYVC3,D%NXSYVC3,D%LFSYVC3,D%LSSYVC3,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VC4') THEN
        t%IECE=.FALSE.
        IF(D%LFVC4.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVC4.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVC4,D%NSYVC4,
     >        D%LXVC4,D%NXVC4,D%LFVC4,
     >        D%LXSYVC4,D%NXSYVC4,D%LFSYVC4,D%LSSYVC4,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VC6') THEN
        t%IECE=.FALSE.
        IF(D%LFVC6.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVC6.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVC6,D%NSYVC6,
     >        D%LXVC6,D%NXVC6,D%LFVC6,
     >        D%LXSYVC6,D%NXSYVC6,D%LFSYVC6,D%LSSYVC6,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VCD') THEN
        t%IECE=.FALSE.
        IF(D%LFVCD.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVCD.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVCD,D%NSYVCD,
     >        D%LXVCD,D%NXVCD,D%LFVCD,
     >        D%LXSYVCD,D%NXSYVCD,D%LFSYVCD,D%LSSYVCD,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VCH') THEN
        t%IECE=.FALSE.
        IF(D%LFVCH.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVCH.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVCH,D%NSYVCH,
     >        D%LXVCH,D%NXVCH,D%LFVCH,
     >        D%LXSYVCH,D%NXSYVCH,D%LFSYVCH,D%LSSYVCH,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VCT') THEN
        t%IECE=.FALSE.
        IF(D%LFVCT.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVCT.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVCT,D%NSYVCT,
     >        D%LXVCT,D%NXVCT,D%LFVCT,
     >        D%LXSYVCT,D%NXSYVCT,D%LFSYVCT,D%LSSYVCT,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VEE') THEN
        t%IECE=.FALSE.
        IF(D%LFVEE.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVEE.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVEE,D%NSYVEE,
     >        D%LXVEE,D%NXVEE,D%LFVEE,
     >        D%LXSYVEE,D%NXSYVEE,D%LFSYVEE,D%LSSYVEE,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VIE') THEN
        t%IECE=.FALSE.
        IF(D%LFVIE.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVIE.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVIE,D%NSYVIE,
     >        D%LXVIE,D%NXVIE,D%LFVIE,
     >        D%LXSYVIE,D%NXSYVIE,D%LFSYVIE,D%LSSYVIE,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VMO') THEN
        t%IECE=.FALSE.
        IF(D%LFVMO.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVMO.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVMO,D%NSYVMO,
     >        D%LXVMO,D%NXVMO,D%LFVMO,
     >        D%LXSYVMO,D%NXSYVMO,D%LFSYVMO,D%LSSYVMO,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VP2') THEN
        t%IECE=.FALSE.
        IF(D%LFVP2.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVP2.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVP2,D%NSYVP2,
     >        D%LXVP2,D%NXVP2,D%LFVP2,
     >        D%LXSYVP2,D%NXSYVP2,D%LFSYVP2,D%LSSYVP2,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VPO') THEN
        t%IECE=.FALSE.
        IF(D%LFVPO.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVPO.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVPO,D%NSYVPO,
     >        D%LXVPO,D%NXVPO,D%LFVPO,
     >        D%LXSYVPO,D%NXSYVPO,D%LFSYVPO,D%LSSYVPO,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VPR') THEN
        t%IECE=.FALSE.
        IF(D%LFVPR.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVPR.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVPR,D%NSYVPR,
     >        D%LXVPR,D%NXVPR,D%LFVPR,
     >        D%LXSYVPR,D%NXSYVPR,D%LFSYVPR,D%LSSYVPR,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'VTR') THEN
        t%IECE=.FALSE.
        IF(D%LFVTR.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXVTR.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIVTR,D%NSYVTR,
     >        D%LXVTR,D%NXVTR,D%LFVTR,
     >        D%LXSYVTR,D%NXSYVTR,D%LFSYVTR,D%LSSYVTR,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'ZF2') THEN
        t%IECE=.FALSE.
        IF(D%LFZF2.le.0) then    ! no data
          write(nonlin,*) ' ?? tdb_profli: lf'//zid//' = 0, no data.'
          return
        endif
        IF(D%NXZF2.le.1) then
          write(nonlin,*) ' -> tdb_profli: nx'//zid//' <=1, no profile.'
          ierr=0
          return
        endif
        CALL TDB_PROFLI(d, t, D%NRIZF2,D%NSYZF2,
     >        D%LXZF2,D%NXZF2,D%LFZF2,
     >        D%LXSYZF2,D%NXSYZF2,D%LFSYZF2,D%LSSYZF2,
     >        ierr)
      ENDIF
C
C------------------
      IF(ZID.EQ.'SIM') THEN
        t%IECE=.FALSE.
        IF(D%LFSIM(I).le.0) return  ! no data
        J = D%NISSIM(I)
        CALL TDB_PROFLI(d, t, D%NRISIM(J),D%NSYSIM(J),
     >        D%LXSIM(J),D%NXSIM(J),D%LFSIM(I),
     >        D%LXSYSIM(J),D%NXSYSIM(J),D%LFSYSIM(I),D%LSSYSIM(I),
     >        ierr)
      ENDIF
C
C=>TRDATGEN-
C
      return
      end
