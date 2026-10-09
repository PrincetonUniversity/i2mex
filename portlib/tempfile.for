C--------------------------------------------------------------------
C  TEMPFILE -- create a temporary filename
C
      subroutine tempfile(zprefix,zfile,ilz)
C
      character*(*) zprefix             ! filename prefix (no path) (input)
      character*(*) zfile               ! full filename (output)
      integer ilz                       ! nonblank length of filename (output)
C
      call tmpfile_d(zprefix,zfile,ilz)
C
      return
      end
