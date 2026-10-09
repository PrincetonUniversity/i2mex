C Close netCDF file
C
C 09/24/97 CAL
C
      subroutine cdfilcl(file_id, ierr)
 
      use cplotr_mod
      include "netcdf.inc"
 
C Input
      integer file_id
C Output
      integer ierr
C
      luntrm=lunzer(0)
C
      ierr = nf_close(file_id)
      if (ierr .ne. NF_NOERR) then
         write(luntrm,9001) file_id,NF_STRERROR(IERR)
 9001    format(' ? CDFILCL: nf_close ',i7/,'     ',A)
         if(ierr.eq.0) ierr=1
      else
         ierr=0
      end if
C
      file_id = 0
C
      return
      end
