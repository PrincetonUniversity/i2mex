      subroutine tget_rlbl(ztok,zlbl)
      use cplotr_mod
c
c  return label information to identify the run currently being accessed.
c
c  ztok -- experiment id ("TFTR","D3D",...)
c
c  zlbl -- (PPPL legacy style):  "tok.yy runid"
c          (MDS+ / CMOD style):  "tok transp##(shot#)"
c
c
      character*(*) ztok                ! exper. id
      character*(*) zlbl                ! run label
c
      character*196 zpath
      character*10 zrunid
c
      call tgetpath(zpath,zrunid)
c
      if(zpath.ne.'MDS+') then
c
c  traditional
c
         zlbl=runlb2
         indx=index(runlb2,'.')-1
         if(indx.le.0) indx=index(runlb2,' ')-1
         ztok=runlb2(1:indx)
c
      else
c
c  MDS+
c
         if(fdir.eq.' ') then
C
C  (update DMC Apr 2009; replacing obsolete code).
C  CMOD style <tree>(<shot>) -- use runlb2, excluding "(MDS+)" suffix.
C
            indx=index(runlb2,'(')-1
            if(indx.gt.0) then
               zlbl = runlb2(1:indx)
            else
               zlbl = runlb2
            endif

            indx=index(zlbl,'.')-1
            if(indx.gt.0) then
               ztok=zlbl(1:indx)
            else
               ztok='????'
            endif
C
         else
C
C  PPPL legacy style <tok.yy> <runid>
C
            ilr=len_trim(runid)
            zlbl=fdir(1:lfdir)//' '//runid(1:ilr)
            indx=index(fdir,'.')-1
            if(indx.le.0) then
               ztok=fdir(1:4)
            else
               ztok=fdir(1:min(4,indx))
            endif
         endif
c
      endif
c
      return
      end

c-----------------------------------------------------------------
      subroutine tget_runid(zrunid)

      use cplotr_mod
c
c  return RUNID string from CPLOTR COMMON
c
      character*(*) :: zrunid

      zrunid = RUNID

      return
      end
