      subroutine rp_geteq(ztime,zdelta,inxi,ntheta,theta,zRbuf,zZbuf,
     >   ierr)

      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
c
c  get the equilibrium R(theta,x),z(theta,x), at a particular time
c  ... set zdelta.gt.0 for time average +/- zdelta (secs) around ztime (secs)
c
c  snippets of code extracted from rplot_sub/plotmm -- dmc 27 Feb 2000
c
c  input:
      real ztime                        ! time at which to fetch equilibrium
      real zdelta                       ! +/- time to average over
c
c  input:
      integer inxi                      ! array dimension, must match inx
      integer ntheta                    ! array dimension, .le. NaxMMP
      real theta(ntheta)                ! theta values at which to evaluate
c
c  output:
      real zRbuf(ntheta,inxi)           ! R(theta,x), cm
      real zZbuf(ntheta,inxi)           ! Z(theta,x), cm
c
      integer ierr                      ! completion code, 0=OK
c
c----------------------------------------------
c  local declarations
c
      logical momrun
c
c----------------------------------------------
c
      ierr=0
c
      IF( (.not.NLTRANSP) .or. (.NOT.MOMRUN(INUM,IMMGEO)) ) THEN
         CALL ZERMSG(
     >      '? RP_GETEQ -- MOMENTS EQUILIBRIUM DATA NOT AVAILABLE!')
         IERR=1
         RETURN
      ENDIF
c
      jtyp=itypr(inum)                  ! x axis type
      inx=nzonex(jtyp)                  ! no. of zones
c
      if(inx.ne.inxi) then
         write(lunzer(0),*) ' ?rp_geteq:  array size mismatch:'
         write(lunzer(0),*)
     >      '  passed inxi=',inxi,' but COMMON inx=',inx
         ierr=1
         return
      endif
c
      if(ntheta.gt.naxmmp) then
         write(lunzer(0),*) ' ?rp_geteq:  COMMON array size exceeded:'
         write(lunzer(0),*) '  need ntheta = ',ntheta,' have NAXMMP = ',
     >      naxmmp
         ierr=1
         return
      endif
c
      IF(NLXFOT(JTYP)) THEN
C
         CALL DMGXOT(JTYP,INDX1,INDX2)
C
         CALL PLPRIO(6,INDX1)
         CALL PLPRIO(6,INDX2)
C
      ENDIF
C
C  READ THE MOMENTS DATA.  WANT BOTH R AND Y DIMENSIONS (IDIM=2)
C  WANT STORAGE PRIORITY TO REMAIN BOOSTED (IBOOST=1)
C
C  NO. OF MOMENTS (IMOM) AND ADDRESSES ARE WRITTEN TO PLFMPA COMMON
C
      IDIM=2
      IBOOST=1
      CALL PLMMRD(IDIM,IBOOST)
c
c  get sin-cos table
c
      call plmcalc_tbl(NaxMom,imom,ntheta,theta,zCosTabl,zSinTabl)
      ncostabl=ntheta
c
      if(immgeo.eq.1) then
c
c  get asymmetric equilibrium
c
         call pltammt(jtyp,ztime,zdelta,zRbuf,zZbuf,ntheta,inxi)
      else
c
c  get symmetric equilibrium
c
         call pltmmt(jtyp,ztime,zdelta,zRbuf,zZbuf,ntheta,inxi)
      endif
c
c  all done
c
C  EXIT
C
C  RESTORE NORMAL STORAGE PRIORITY TO MOMENTS AND X AXIS DATA
C
      CALL PLPRIO(5,INDX1)
      CALL PLPRIO(5,INDX2)
C
      IDIM=2
      CALL PLMMRLS(IDIM)
C
      return
      end
c------------------------------------------------
      subroutine rp_eq_symflag(isym)
      use cplotr_mod
      use plfmpa_mod
c
c  for use after rp_geteq... returns the symmetry code
c
      integer isym                      ! output symmetry flag, 0=symmetric
c
c  on output, isym=0 means the equilibrium is up-down symmetric;
c             isym=1 means the equilibrium is up-down asymmetric
c
      isym=immgeo
c
      return
      end
 
c------------------------------------------------
      subroutine rp_eq_nmoms(imoms)
      use cplotr_mod
      use plfmpa_mod
c
c  for use after rp_geteq... returns the no. of moments in the run
c
      integer imoms                     ! output no. of moments
c
c  on output, isym=0 means the equilibrium is up-down symmetric;
c             isym=1 means the equilibrium is up-down asymmetric
c
      imoms=imom
c
      return
      end
 
