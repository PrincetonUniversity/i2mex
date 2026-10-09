C Open netCDF file
C
C 09/24/97 CAL
C
      subroutine cdfilop(fname, file_id, ierr)
 
      include "netcdf.inc"
 
C Input
      character*(*) fname
C Output
      integer file_id
      integer ierr
C
C Local
      character*80 zstr
C
C--------------
C
        luntrm=lunzer(0)
C
        ierr = nf_open(fname,NF_NOWRITE,file_id)
        if (ierr .eq. -7777) return     ! no NetCDF library
        if (ierr .ne. NF_NOERR) then
           zstr=NF_STRERROR(IERR)
           if((index(zstr,'No such file or directory').eq.0).and.
     >        (index(zstr,'no such file or directory').eq.0)) then
              write(luntrm,9001) zstr
 9001         format(' cdfilop error: '/,'     ',A)
           endif
           if(ierr.eq.0) ierr=1
        else
           ierr=0
        end if
C
        return
        end
 
