! 30Nov1999 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
      subroutine r8tksbcon(xiblo,rmb0,mtbl,rmb2,ymb2,mimom,nmom,
     >                   jzb,jzb0,ziter)
C
C  control routine for Fourier Moments surface extrapolation,
C  up-down symmetry assumed.  Arguments have same meanings as in caller
C  bloats.for, plus:
C
C   jzb -- extrapolated surface to compute now
C   jzb0 -- 1st extrapolated surface
C   ziter -- numerical control
C
C  DMC mod-- constrain # of moments in extrapolation (bugfix for stability)
C            use new local variable "inmom"
C
C  passed:
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER mimom,nmom,jzb,jzb0,mtbl,locmax,ics,im,ii,isurf0,iter
      INTEGER inmom
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 ziter
C============
      REAL*8 xiblo(mtbl),rmb0(mtbl)
      REAL*8 rmb2(mtbl,mimom),ymb2(mtbl,mimom)
C
C  local:  storage for 3 surfaces worth of moments...
      parameter (locmax=64)
      REAL*8 zrmc(3,0:locmax,2)
      REAL*8 zymc(3,0:locmax,2)
C
      logical nlsym
      data nlsym/.true./
C
C--------------------------------------------
C
      if(mimom.gt.locmax)
     >   call errmsg_exit('tksbcon:  locmax too small!')
C
C  copy in
C
      do ics=1,2
         do im=0,nmom
            do ii=1,3
               zrmc(ii,im,ics)=0.0E0_R8
               zymc(ii,im,ics)=0.0E0_R8
            enddo
         enddo
      enddo
C
      isurf0=jzb-3
      do ii=1,2
         zrmc(ii,0,1)=rmb0(isurf0+ii)
         inmom=nmom
         if(isurf0+ii.ge.jzb0) inmom=min(16,inmom)
         if(isurf0+ii.ge.jzb0+2) inmom=min(8,inmom)
         do im=1,inmom
            zrmc(ii,im,1)=rmb2(isurf0+ii,im)
            zymc(ii,im,2)=ymb2(isurf0+ii,im)
         enddo
      enddo
C
      iter=jzb-jzb0
      call r8tkbcon(nlsym,xiblo,jzb,zrmc,zymc,locmax,inmom,iter,ziter)
C
C  copy back
C
      rmb0(jzb)=zrmc(3,0,1)
      do im=1,nmom
         rmb2(jzb,im)=zrmc(3,im,1)
         ymb2(jzb,im)=zymc(3,im,2)
      enddo
C
      return
      end
 
 
