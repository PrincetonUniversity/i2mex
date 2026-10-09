      subroutine kconnect(path,rlabel,ntscal,ntprof,nxmax,ndmax,ier)
c
c  connect to a TRANSP run; return a generic label for the run,
c  and max dimensions for timebases, x axes and data items.
c
      character*(*) path                ! (IN) path to run, path/runid
c     ...cf subroutine RCONNECT, rplot_sub.
c
      character*(*) rlabel              ! (OUT) generic run label, ~C*60
c
      integer ntscal                    ! (OUT) #pts in scalar timebase
      integer ntprof                    ! (OUT) #pts in profile timebase
      integer nxmax                     ! (OUT) max #pts in any x axis
      integer ndmax                     ! (OUT) max size of any one item
c
c  usually ndmax = nxmax*ntprof
c
      integer ier                       ! (OUT) completion code:  0=OK
c
c-------------------------------------------------------
c
      integer istart
c
      data istart/0/
c
c-------------------------------------------------------
c
c  initialize COMMON if necessary
c
      if(istart.eq.0) then
         istart=1
         call initcpl
      endif
c
c  connect to a run
c
      ntscal=0
      ntprof=0
      nxmax=0
      ndmax=0
c
      call rconnect(path,ier)
c
      if(ier.ne.0) return
c
c-------------------------------------------------------
c  clear namelist
c
      call tr_getnl_clear
c
c-------------------------------------------------------
c
c  get generic label
c
      call grunlb2(rlabel)
c
c-------------------------------------------------------
c
c  get timebase xbase and data item sizes
c
      call rpstats(ntscal,ntprof,nxmax,ndmax)
c
c-------------------------------------------------------
c
      return
      end
