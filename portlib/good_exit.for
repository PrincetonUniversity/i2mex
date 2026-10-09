      subroutine good_exit
C
C  normal exit from an f77 program
C
      call err_end                      ! turn off SUN fpe error handling msg
C
      stop  ! normal exit, stop, status = 0

      end
