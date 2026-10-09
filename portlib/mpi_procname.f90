      subroutine mpi_procname(myid,numprocs,filenam)

#ifdef __MPI
      use mpi
#endif
      implicit NONE

! local communicator

#ifdef __MPI
     integer :: nbi_comm = MPI_COMM_WORLD
#endif

     !-------------------------------------
     character*(*), intent(in) :: filenam
     character*180 hostnam,hostnam1,proc_name
     integer :: istat,myid,numprocs,ierr,jwrite
     integer :: ii,jj,mywid,maxwid,iresult
     integer :: lunout = 77
     integer, dimension(:), allocatable :: iwrite

     integer, dimension(:,:), allocatable :: intbuf,intbufr
     !-------------------------------------
     ! executable code

     open(unit=lunout,file=trim(adjustl(filenam)),status='replace',iostat=istat)

     if (allocated(iwrite))deallocate(iwrite)
     allocate(iwrite(0:numprocs-1))
     iwrite(:)=0

#ifdef __MPI

     call MPI_GET_PROCESSOR_NAME(proc_name,iresult,ierr)
     call ckerr('MPI_GET_PROCESSOR_NAME')
!    write(0,'(1x,i8,1x,a)') myid, trim(proc_name)

     mywid = len(trim(proc_name))

     call MPI_ALLREDUCE(mywid,maxwid,1, &
       MPI_INTEGER,MPI_MAX,nbi_comm,ierr)
     call ckerr('MPI_ALLREDUCE(mywid)')

     allocate(intbuf(maxwid,0:numprocs-1)); intbuf = 0
     allocate(intbufr(maxwid,0:numprocs-1))

     do ii=1,mywid
       intbuf(ii,myid) = ichar(proc_name(ii:ii))
     enddo

     call MPI_ALLREDUCE(intbuf,intbufr,numprocs*maxwid, &
       MPI_INTEGER,MPI_SUM,nbi_comm,ierr)
     call ckerr('MPI_REDUCE(intbuf)')

!    if(myid.eq.0) then
      hostnam1=' '
      do jj=0,numprocs-1
        iwrite(jj)=0
        hostnam=' '
        do ii=1,maxwid
           if(intbufr(ii,jj).eq.0) exit
           hostnam(ii:ii) = char(intbufr(ii,jj))
        enddo
        
        if (trim(hostnam).ne.trim(hostnam1))iwrite(jj)=jj
        if (jj.eq.0) iwrite(jj)=0
        hostnam1=hostnam

        if (jj.eq.myid) then
          jwrite=iwrite(jj)
        else
          jwrite=0
        endif
!       if (myid.eq.0) then
!        write(0,'(1x,i8,1x,a,1x,i8)') jj, trim(hostnam), iwrite(jj)
!       endif
      enddo
!    endif

      write(lunout,'(i8)') iwrite(myid) 

#else
     call hostnm(hostnam)
    
!    write(lunout,'(1x,i8,1x,a,1x,i8)') myid, trim(hostnam),iwrite(myid) 
     write(lunout,'(i8)') iwrite(myid) 
#endif

      close(unit=lunout)
 
      deallocate(iwrite)

CONTAINS

      subroutine ckerr(callnam)

        character*(*) :: callnam

        if(ierr.ne.0) then
         write(0,*) ' error ierr=',ierr,' in call "',trim(callnam), &
            '" on proc #',myid
         call bad_Exit
        endif

       end subroutine ckerr


      end subroutine mpi_procname
