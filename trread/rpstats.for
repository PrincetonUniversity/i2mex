      subroutine rpstats(isctime,iprtime,ixmax,ifmax)
C
C  return some basic size statistics
C
      use cplotr_mod

      integer isctime                   ! no. of scalar timepoints
      integer iprtime                   ! no. of profile timepoints
      integer ixmax                     ! largest number of points, any x axis
      integer ifmax                     ! largest size of any function
C
      isctime=ntt
      iprtime=ntr
C
      ifmax=isctime                     ! size of scalar functions
      ixmax=0
C
      do iax=1,nxr
         ixmax=max(ixmax,nzonex(iax))
         ifmax=max(ifmax,ntr*nzonex(iax))
      enddo
C
      return
      end
