      subroutine bad_exit
C
#ifdef __MPI
      write(6,*) 
     >' %bad_exit:  generic f77 error exit call (errset_mpi status 999)'
#else
      write(6,*) 
     >' %bad_exit:  generic f77 error exit call (errset status 999)'
#endif
C
      call err_end
      flush(6)
C
      call errset_mpi(-1,999)
C
      STOP                              ! this line not reached.
      END
