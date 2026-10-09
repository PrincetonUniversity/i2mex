      subroutine mds_cache_fname(ixflag,ztail,zfiln,ier)
c
      use cplotr_mod
c
c  form a full filename, given the tail name (i.e. -- prepend the
c  path)
c
c  input:
      integer ixflag                    ! primary / 2ndary run index
      character*(*) ztail               ! tail filename
c
c  output:
      character*(*) zfiln               ! full filename
      integer ier                       ! error code, 0=OK
c
c--------------------
c
      character*200 zroot,zpath
      character(len=200) :: zcmd
c
      integer str_length
c---------------------
c
      if(ixflag.eq.0) then
         call gmkdir_init(mds_cache_root,zroot,mds_cache_dir,zpath)
      else
         call gmkdir_init(mds_cache_root,zroot,mds_cache_dir_x(ixflag),
     >       zpath)
      endif
c
      ilroot=str_length(zroot)
c
      ilpath=str_length(zpath)
c
      iltail=str_length(ztail)
c
      if(ilroot+ilpath+iltail+3.gt.len(zfiln)) then
         write(lunzer(0),*)
     >      ' ?mds_cache_fname:  character buffer exceeded'
         write(lunzer(0),*)
     >      '  root:  ',zroot(1:ilroot)
         write(lunzer(0),*)
     >      '  path:  ',zpath(1:ilpath)
         ier=1
         return
      else
         ier=0
      endif
c
      zfiln=zroot(1:ilroot)//zpath(1:ilpath)//'/'//ztail(1:iltail)
      return
      end
