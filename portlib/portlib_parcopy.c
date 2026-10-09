/* Adapted from pcp from SUT, in MPI */ 
#ifdef __MPI
#include "mpi.h" 
#endif
#include <stdio.h> 
#include <sys/types.h> 
#include <sys/stat.h> 
#include <fcntl.h> 
#include "f77name.h"

#define BUFSIZE    256*1024 
#define CMDSIZE    80 

/*  
    This is an MPI program that copies a file (any format) 
    'cinput' accessible by proc of rank (rd_rank) into a file 'coutput' on the
    file space seen by (wr_rank). Note that if the file names coutput 
    and cinput are different, only proc (rd_rank) has a file named cinput 
    while (wr_rank) has a file named coutput. The contents, however, are 
    the same.

    DMC May 2011: only using MPI_SEND & MPI_RECV here.  If myrank is neither
    rd_rank nor wr_rank the routine immediately returns; there is no collective
    synchronization imposed here.

    npcopy_comm: the communicator passed from the calling program
    myrank: rank of the calling process
    rd_rank: rank of process that reads (sends) the file (named cinput)
    wr_rank: rank of process that writes (receives) the file (named coutput)

    cinput: input file name to be copied
    coutput: output file name; coutput can be the same as cinput
*/
int F77NAME(portlib_parcopy)(npcopy_comm,myrank,rd_rank,wr_rank,cinput,coutput)
     int *npcopy_comm;
     int *myrank;
     int *rd_rank;
     int *wr_rank;

     const char *cinput;
     const char *coutput;
{ 
  int mystatus, done, numread; 

  int infd, outfd, iwrite, istat;
  char buf[BUFSIZE]; 
#ifdef __MPI
  MPI_Comm new_comm;
  MPI_Status status;
#endif
    
  /*  ranks exit exit immediately if not involved in transfer  */

  if ( *myrank != *rd_rank && *myrank != *wr_rank )   return( 0 );

#ifdef __MPI
  new_comm=MPI_Comm_f2c(*npcopy_comm);
#endif

  /*   Set iwrite to 1 if want diagnostic output */

  iwrite=0;

  /*   Open the input file */ 
  /*   If input file open fails, exit */

  if ( *myrank == *rd_rank ) {  
    if ( (infd = open( cinput, O_RDONLY ) ) == -1 ) {
      fprintf( stderr, "input file %s does not exist on %d\n", cinput, *rd_rank );
      mystatus = -1;
    }
    else {
      mystatus = 0;
      if (iwrite==1) fprintf( stderr, "myrank: %d - input file - %s - exists\n", *myrank, cinput );
    }
    /* file open handshake */
#ifdef __MPI
    if ( *rd_rank != *wr_rank ) {
      istat = MPI_Send( &mystatus, 1, MPI_INT, *wr_rank, 1000, new_comm );
      istat = MPI_Recv( &mystatus, 1, MPI_INT, *wr_rank, 1001, new_comm, &status );
    }
#endif
  }
  if ( *myrank == *wr_rank ) {
#ifdef __MPI
    if ( *rd_rank != *wr_rank ) {
      istat = MPI_Recv( &mystatus, 1, MPI_INT, *rd_rank, 1000, new_comm, &status );
    }
#endif
    if ( mystatus == 0 ) {
      if ( (outfd = open( coutput, O_CREAT|O_WRONLY|O_TRUNC, 0644 ) ) == -1 ) {
	fprintf( stderr, "cannot open output file %s on %d\n", coutput, *rd_rank );
	mystatus = -1;
      }
      else {
	mystatus = 0;
	  if (iwrite==1) fprintf( stderr, "myrank: %d - output file - %s - open\n", *myrank, coutput );
      }
    }
#ifdef __MPI
    if ( *rd_rank != *wr_rank ) {
      istat = MPI_Send( &mystatus, 1, MPI_INT, *rd_rank, 1001, new_comm );
    }
#endif
  }

  if ( mystatus < 0 ) return ( -1 );

  /*  OK files opened successfully  */

  /*At this point, *rd_rank has file open for input; *wr_rank has file open
    for output*/

  if (iwrite==1) fprintf( stderr, "myrank %d, reached read/write loop.\n", *myrank);

  mystatus = 0;
  done = 0;

  while ( !done ) {
    if ( *myrank == *rd_rank ) {
      numread = read( infd, buf, BUFSIZE );
#ifdef __MPI
      if ( *rd_rank != *wr_rank ) {
	istat = MPI_Send( &numread, 1, MPI_INT, *wr_rank, 1002, new_comm );
      }
#endif
      if ( numread > 0 ) {
#ifdef __MPI
	if ( *rd_rank != *wr_rank ) {
	  istat = MPI_Send( buf, numread, MPI_BYTE, *wr_rank, 1003, new_comm );
	}
#endif
      }
      else {
	close( infd );
	done = 1;
      }
    }
    if ( *myrank == *wr_rank ) {
#ifdef __MPI
      if ( *rd_rank != *wr_rank ) {
	istat = MPI_Recv( &numread, 1, MPI_INT, *rd_rank, 1002, new_comm, &status );
      }
#endif
      if ( numread > 0 ) {
#ifdef __MPI
	if ( *rd_rank != *wr_rank ) {
	  istat = MPI_Recv( buf, numread, MPI_BYTE, *rd_rank, 1003, new_comm, &status );
	}
#endif
	write( outfd, buf, numread );
      }
      else {
	close( outfd );
	done = 1;
      }
    }
  }

  if (iwrite==1) fprintf( stderr, "myrank %d, completed read/write loop.\n", *myrank);
  return 0;

} 

int F77NAME(portlib_parcopy_)(comm,myrank,rd_rank,wr_rank,cinput,coutput)
     int *comm;
     int *myrank;
     int *rd_rank;
     int *wr_rank;

     const char *cinput;
     const char *coutput;
{ 
  F77NAME(portlib_parcopy)(comm,myrank,rd_rank,wr_rank,cinput,coutput);
}
