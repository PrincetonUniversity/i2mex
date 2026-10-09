      subroutine parcopy(npcopy_comm,myid,read_id,write_id,inputfile, &
           outputfile,ierr)

      use iso_c_binding, only: c_char
      implicit none

      integer :: npcopy_comm  ! MPI communicator
      integer :: myid         ! current process ID
      integer :: read_id      ! ID of process that should READ file
      integer :: write_id     ! ID of process that should WRITE file

      character*(*),intent(in) :: inputfile,outputfile  ! filenames

      integer, intent(out) :: ierr

      !---------------
      integer :: portlib_parcopy
      logical :: idebug = .FALSE.
      !---------------


      character(kind=c_char) cinput(1+len(inputfile))
      character(kind=c_char) coutput(1+len(outputfile))

      if(idebug) then
         if(myid.eq.read_id) then 
            write(0,*) ' I, ',myid,' am reader: '
            write(0,*) '     inputfile: '//trim(inputfile)
            write(0,*) '    outputfile: '//trim(outputfile)
         endif
         if(myid.eq.write_id) then 
            write(0,*) ' I, ',myid,' am writer: '
            write(0,*) '     inputfile: '//trim(inputfile)
            write(0,*) '    outputfile: '//trim(outputfile)
         endif
      endif

      call cstring(inputfile(1:len_trim(inputfile)),cinput,'2C')
      call cstring(outputfile(1:len_trim(outputfile)),coutput,'2C')

      ierr = portlib_parcopy(npcopy_comm,myid,read_id,write_id,cinput,coutput)

      return 
!     stop
      end
!

