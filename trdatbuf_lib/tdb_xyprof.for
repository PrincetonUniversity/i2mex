      subroutine tdb_xyprof(d,zid,ztime1,ztime2,inx,iny,zdata,ierr)
 
      !  time interpolate f(t,x,y) profile to zdata(x,y) -- no space interp.
      !  data array sizes (inx,iny) must match.
 
      !  if ztime1.ne.ztime2, time average between the two times.
 
      use trdatbuf_obj
 
      implicit NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
      type (trdatbuf) :: d
 
      character*(*), intent(in) :: zid  ! data identification e.g. "PSI"
      real*8, intent(in) :: ztime1,ztime2  ! target time(s) for interpolation
      ! or averaging
 
      integer, intent(in) :: inx   ! size of first grid e.g. R if type f(R,Z)
      integer, intent(in) :: iny   ! size of 2nd grid e.g. Z if type f(R,Z)
 
      real*8, intent(out) :: zdata(inx,iny)  ! time interpolated data
 
      integer, intent(out) :: ierr  ! completion code (0=OK)
 
      real*8 :: tdbsub_i1  ! time integration utility routine
      !---------------------------------------
      ! local information:
 
      character*10 :: zidi
 
      integer :: nonlin,lunmsg_tdb,imatch
      integer :: inxd,inyd
      integer :: intd,it1,ia1,ia2,ix,iy
      real*8 :: zfract,zdat1,zdat2,zt1,zt2
 
      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter ::  ONE = 1.0d0
      !---------------------------------------
 
      nonlin = lunmsg_tdb(0)
 
      zidi = zid
      call uupper(zidi)
 
      imatch = 0
 
c  set initial values for output variables...
 
      ierr=0
      zdata = ZERO
 
      zt1=max(ztime1,ztime2)
      zt2=min(ztime1,ztime2)
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
C  FD0
      if(zidi.eq.'FD0') then
         imatch=1
         if(d%LFFD0.le.0) then
            write(nonlin,*) ' %tdb_xyprof: no '//trim(zidi)//' data.'
            return
         endif
         intd = d%NTFD0
         inxd = d%NEFD0
         inyd = d%NXFD0
         if((inx.ne.inxd).or.(iny.ne.inyd)) then
            write(nonlin,*) 'error in tdb_xyprof call: FD0'
            write(nonlin,*) ' array sizes mismatch:'
            write(nonlin,*) ' data nx = ',inxd,'; passed nx=',inx
            write(nonlin,*) ' data ny = ',inyd,'; passed ny=',iny
            ierr=1
            return
         endif
         if(zt2.eq.zt1) then
            call tdbsub_slookup(d%datbuf(d%LTFD0),intd,zt1,it1,zfract)
         else
            it1=1
         endif
         do iy = 1,iny
            do ix = 1,inx
               ia1 = ((iy-1)*inx + (ix-1))*intd + (it1-1)
               ia2 = ia1 + 1
               if(zt2.eq.zt1) then
                  zdat1 = d%datbuf(d%LFFD0+ia1)
                  zdat2 = d%datbuf(d%LFFD0+ia2)
                  zdata(ix,iy) = zdat1*(ONE-zfract) + zfract*zdat2
               else
                  zdata(ix,iy) = tdbsub_i1(zt1,zt2,
     >                 d%datbuf(d%LTFD0),intd,d%datbuf(d%LFFD0+ia1))/
     >                 (zt2-zt1)
               endif
            enddo
         enddo
      endif
C  FDB
      if(zidi.eq.'FDB') then
         imatch=1
         if(d%LFFDB.le.0) then
            write(nonlin,*) ' %tdb_xyprof: no '//trim(zidi)//' data.'
            return
         endif
         intd = d%NTFDB
         inxd = d%NEFDB
         inyd = d%NXFDB
         if((inx.ne.inxd).or.(iny.ne.inyd)) then
            write(nonlin,*) 'error in tdb_xyprof call: FDB'
            write(nonlin,*) ' array sizes mismatch:'
            write(nonlin,*) ' data nx = ',inxd,'; passed nx=',inx
            write(nonlin,*) ' data ny = ',inyd,'; passed ny=',iny
            ierr=1
            return
         endif
         if(zt2.eq.zt1) then
            call tdbsub_slookup(d%datbuf(d%LTFDB),intd,zt1,it1,zfract)
         else
            it1=1
         endif
         do iy = 1,iny
            do ix = 1,inx
               ia1 = ((iy-1)*inx + (ix-1))*intd + (it1-1)
               ia2 = ia1 + 1
               if(zt2.eq.zt1) then
                  zdat1 = d%datbuf(d%LFFDB+ia1)
                  zdat2 = d%datbuf(d%LFFDB+ia2)
                  zdata(ix,iy) = zdat1*(ONE-zfract) + zfract*zdat2
               else
                  zdata(ix,iy) = tdbsub_i1(zt1,zt2,
     >                 d%datbuf(d%LTFDB),intd,d%datbuf(d%LFFDB+ia1))/
     >                 (zt2-zt1)
               endif
            enddo
         enddo
      endif
C  FDP
      if(zidi.eq.'FDP') then
         imatch=1
         if(d%LFFDP.le.0) then
            write(nonlin,*) ' %tdb_xyprof: no '//trim(zidi)//' data.'
            return
         endif
         intd = d%NTFDP
         inxd = d%NEFDP
         inyd = d%NXFDP
         if((inx.ne.inxd).or.(iny.ne.inyd)) then
            write(nonlin,*) 'error in tdb_xyprof call: FDP'
            write(nonlin,*) ' array sizes mismatch:'
            write(nonlin,*) ' data nx = ',inxd,'; passed nx=',inx
            write(nonlin,*) ' data ny = ',inyd,'; passed ny=',iny
            ierr=1
            return
         endif
         if(zt2.eq.zt1) then
            call tdbsub_slookup(d%datbuf(d%LTFDP),intd,zt1,it1,zfract)
         else
            it1=1
         endif
         do iy = 1,iny
            do ix = 1,inx
               ia1 = ((iy-1)*inx + (ix-1))*intd + (it1-1)
               ia2 = ia1 + 1
               if(zt2.eq.zt1) then
                  zdat1 = d%datbuf(d%LFFDP+ia1)
                  zdat2 = d%datbuf(d%LFFDP+ia2)
                  zdata(ix,iy) = zdat1*(ONE-zfract) + zfract*zdat2
               else
                  zdata(ix,iy) = tdbsub_i1(zt1,zt2,
     >                 d%datbuf(d%LTFDP),intd,d%datbuf(d%LFFDP+ia1))/
     >                 (zt2-zt1)
               endif
            enddo
         enddo
      endif
C  FDQ
      if(zidi.eq.'FDQ') then
         imatch=1
         if(d%LFFDQ.le.0) then
            write(nonlin,*) ' %tdb_xyprof: no '//trim(zidi)//' data.'
            return
         endif
         intd = d%NTFDQ
         inxd = d%NEFDQ
         inyd = d%NXFDQ
         if((inx.ne.inxd).or.(iny.ne.inyd)) then
            write(nonlin,*) 'error in tdb_xyprof call: FDQ'
            write(nonlin,*) ' array sizes mismatch:'
            write(nonlin,*) ' data nx = ',inxd,'; passed nx=',inx
            write(nonlin,*) ' data ny = ',inyd,'; passed ny=',iny
            ierr=1
            return
         endif
         if(zt2.eq.zt1) then
            call tdbsub_slookup(d%datbuf(d%LTFDQ),intd,zt1,it1,zfract)
         else
            it1=1
         endif
         do iy = 1,iny
            do ix = 1,inx
               ia1 = ((iy-1)*inx + (ix-1))*intd + (it1-1)
               ia2 = ia1 + 1
               if(zt2.eq.zt1) then
                  zdat1 = d%datbuf(d%LFFDQ+ia1)
                  zdat2 = d%datbuf(d%LFFDQ+ia2)
                  zdata(ix,iy) = zdat1*(ONE-zfract) + zfract*zdat2
               else
                  zdata(ix,iy) = tdbsub_i1(zt1,zt2,
     >                 d%datbuf(d%LTFDQ),intd,d%datbuf(d%LFFDQ+ia1))/
     >                 (zt2-zt1)
               endif
            enddo
         enddo
      endif
C  FDR
      if(zidi.eq.'FDR') then
         imatch=1
         if(d%LFFDR.le.0) then
            write(nonlin,*) ' %tdb_xyprof: no '//trim(zidi)//' data.'
            return
         endif
         intd = d%NTFDR
         inxd = d%NEFDR
         inyd = d%NXFDR
         if((inx.ne.inxd).or.(iny.ne.inyd)) then
            write(nonlin,*) 'error in tdb_xyprof call: FDR'
            write(nonlin,*) ' array sizes mismatch:'
            write(nonlin,*) ' data nx = ',inxd,'; passed nx=',inx
            write(nonlin,*) ' data ny = ',inyd,'; passed ny=',iny
            ierr=1
            return
         endif
         if(zt2.eq.zt1) then
            call tdbsub_slookup(d%datbuf(d%LTFDR),intd,zt1,it1,zfract)
         else
            it1=1
         endif
         do iy = 1,iny
            do ix = 1,inx
               ia1 = ((iy-1)*inx + (ix-1))*intd + (it1-1)
               ia2 = ia1 + 1
               if(zt2.eq.zt1) then
                  zdat1 = d%datbuf(d%LFFDR+ia1)
                  zdat2 = d%datbuf(d%LFFDR+ia2)
                  zdata(ix,iy) = zdat1*(ONE-zfract) + zfract*zdat2
               else
                  zdata(ix,iy) = tdbsub_i1(zt1,zt2,
     >                 d%datbuf(d%LTFDR),intd,d%datbuf(d%LFFDR+ia1))/
     >                 (zt2-zt1)
               endif
            enddo
         enddo
      endif
C  FDS
      if(zidi.eq.'FDS') then
         imatch=1
         if(d%LFFDS.le.0) then
            write(nonlin,*) ' %tdb_xyprof: no '//trim(zidi)//' data.'
            return
         endif
         intd = d%NTFDS
         inxd = d%NEFDS
         inyd = d%NXFDS
         if((inx.ne.inxd).or.(iny.ne.inyd)) then
            write(nonlin,*) 'error in tdb_xyprof call: FDS'
            write(nonlin,*) ' array sizes mismatch:'
            write(nonlin,*) ' data nx = ',inxd,'; passed nx=',inx
            write(nonlin,*) ' data ny = ',inyd,'; passed ny=',iny
            ierr=1
            return
         endif
         if(zt2.eq.zt1) then
            call tdbsub_slookup(d%datbuf(d%LTFDS),intd,zt1,it1,zfract)
         else
            it1=1
         endif
         do iy = 1,iny
            do ix = 1,inx
               ia1 = ((iy-1)*inx + (ix-1))*intd + (it1-1)
               ia2 = ia1 + 1
               if(zt2.eq.zt1) then
                  zdat1 = d%datbuf(d%LFFDS+ia1)
                  zdat2 = d%datbuf(d%LFFDS+ia2)
                  zdata(ix,iy) = zdat1*(ONE-zfract) + zfract*zdat2
               else
                  zdata(ix,iy) = tdbsub_i1(zt1,zt2,
     >                 d%datbuf(d%LTFDS),intd,d%datbuf(d%LFFDS+ia1))/
     >                 (zt2-zt1)
               endif
            enddo
         enddo
      endif
C  PSI
      if(zidi.eq.'PSI') then
         imatch=1
         if(d%LFPSI.le.0) then
            write(nonlin,*) ' %tdb_xyprof: no '//trim(zidi)//' data.'
            return
         endif
         intd = d%NTPSI
         inxd = d%NRPSI
         inyd = d%NZPSI
         if((inx.ne.inxd).or.(iny.ne.inyd)) then
            write(nonlin,*) 'error in tdb_xyprof call: PSI'
            write(nonlin,*) ' array sizes mismatch:'
            write(nonlin,*) ' data nx = ',inxd,'; passed nx=',inx
            write(nonlin,*) ' data ny = ',inyd,'; passed ny=',iny
            ierr=1
            return
         endif
         if(zt2.eq.zt1) then
            call tdbsub_slookup(d%datbuf(d%LTPSI),intd,zt1,it1,zfract)
         else
            it1=1
         endif
         do iy = 1,iny
            do ix = 1,inx
               ia1 = ((iy-1)*inx + (ix-1))*intd + (it1-1)
               ia2 = ia1 + 1
               if(zt2.eq.zt1) then
                  zdat1 = d%datbuf(d%LFPSI+ia1)
                  zdat2 = d%datbuf(d%LFPSI+ia2)
                  zdata(ix,iy) = zdat1*(ONE-zfract) + zfract*zdat2
               else
                  zdata(ix,iy) = tdbsub_i1(zt1,zt2,
     >                 d%datbuf(d%LTPSI),intd,d%datbuf(d%LFPSI+ia1))/
     >                 (zt2-zt1)
               endif
            enddo
         enddo
      endif
C
C=>TRDATGEN-
C
      if(imatch.eq.0) then
         write(nonlin,*) 'error in tdb_xyprof subroutine:'
         write(nonlin,*) ' ID '//trim(zidi)//' not recognized.'
         write(nonlin,*) ' ID of f(t,x,y) profile was expected.'
         ierr=2
      endif
 
      return
      end
