      subroutine rpcalcd(zexpr, rlabel,
     >   zdata, zxdata, ndmax, ztdata, ntmax,
     >   nxgot, ntgot, istype, iwarn, ier)
c
c  call the rplot calculator routine rpcalc.
c  return a copy of the data computed, the corresponding x axis data
c  if applicable, and the corresponding timebase data.
c
c  input:
      character*(*) zexpr               ! calculator expression (cf rpcalc).
      character*(*) rlabel              ! run label (for error message)
c
c  storage buffers filled on output:
      integer ndmax                     ! max size of data item (INPUT)
      integer ntmax                     ! max size of timebase (INPUT)
c
      real zdata(ndmax)                 ! buffer for calculator results data
      real zxdata(ndmax)                ! buffer for x axis data
      real ztdata(ntmax)                ! buffer for timebase
c
c  scalars set on output:
c
      integer nxgot                     ! no. of x pts, 1 for scalar f(t) data
      integer ntgot                     ! no. of time pts in timebase
c
c  ntgot words are written in ztdata;
c  nxgot*ntgot words are written in zdata
c  if (istype.ge.1) nxgot*ntgot words are written in zxdata
c  ...the x variation is stored contiguously, i.e. the fortran 2d
c  array declaration would be array(nxgot,ntgot)
c
      integer istype                    ! data item type code
c
c  istype=-1 -- scalar f(t)
C  istype= 1 -- f(x,t), x is TRANSP zone ctrs if this is TRANSP data
C  istype= 2 -- f(x,t), x is TRANSP zone bdys if this is TRANSP data
C  ..etc..
C  istype = 0 indicates an error
C
      integer iwarn                     ! arithmetic warning, 0 = OK
      integer ier                       ! completion code, 0 = OK
c
c--------------------------------------------------------------------
c
      integer lunzer
c
      integer irank,idims(10)
      character*10 zxabb(8)             ! higher rank objects, eventually...
c
c--------------------------------------------------------------------
c
c  OK... call the calculator
c
      ilunz=lunzer(0)
c
      call rpcalc(zexpr,zdata,ndmax,igot,istype,iwarn,ier)
c
c  message on warning or error
c
      if(max(ier,iwarn).ne.0) then
         ile=len_trim(zexpr)
         ilb=len_trim(rlabel)
         write(ilunz,1001) rlabel(1:ilb),zexpr(1:ile),iwarn,ier
 1001    format(/
     >      '%rpcalcd:  runid:  ',a/
     >      ' expression:  ',a/
     >      ' iwarn=',i6,'      ier=',i6)
      endif
      if(ier.ne.0) then
         nxgot=0
         ntgot=0
         return
      endif
c
      if(istype.eq.-1) then
c
c  if a scalar function result:  have data already, just get timebase
c
         nxgot=1
         ntgot=igot
         call rptime_s(ztdata,ntmax,igot)
         return
      else
c
c  if a profile function result:  have data already, get profile
c  timebase and x axis function
c
         call rpdims(istype,irank,idims,zxabb,ier)
         if(irank.gt.2) then
            write(ilunz,*)
     >         ' ??rpcalcd:  irank.gt.2 object not supported!'
            nxgot=0
            ntgot=0
            ier=99
         else
            nxgot=idims(1)
            ntgot=idims(2)
         endif
c  check time dimension
         if(ntgot.gt.ntmax) then
            write(ilunz,*) ' ??rpcalcd:  ntgot.gt.ntmax, ntgot=',ntgot,
     >         ' ntmax=',ntmax
            nxgot=0
            ntgot=0
            ier=88
         else
c  get timebase
            call rptime_p(ztdata,ntmax,igot)
c  get x axis
            call rprofile(zxabb,zxdata,ndmax,iret, ier)
         endif
c
      endif
c
      return
      end
