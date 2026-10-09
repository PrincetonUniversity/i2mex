/* clocaltim.c:
 
Fortran interface to c localtime library routine
 
   (c) D. R. Ernst <dernst@pppl.gov>  7/29/95, 5/10/99
   Permission granted only for royalty-free not-for profit distribution.
 
   modified from original c24date.c by D. McCune
 
   integer itime(8)
   call clocaltim(itime)
 
   output
 
     itime(1) = 4 digit year             ***> C code:  itime[0]
     itime(2) = month (1=Jan to 12=Dec)
     itime(3) = day of month (1 to 31)
     itime(4) = 0
     itime(5) = hours since midnight (0 to 23)
     itime(6) = minutes after hour (0 to 59)
     itime(7) = seconds after minute (0 to 59); 60 or 61 if leap seconds
     itime(8) = 0
 
  the indexing is compatible with the itime array output of the fortran-90
  date and time routine, BUT:  elements 4 and 8 are not set.
 
*/
 
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include "f77name.h"
 
void F77NAME(clocaltim)(itime)
int itime[];
{
#define LEN 24
 
    time_t TimeNow;
    struct tm *TimeNowStruct;
 
#ifdef __DEBUG_DATE
    printf("%s\n", "entry" );
#endif
 
    time ( &TimeNow );
 
#ifdef __DEBUG_DATE
    printf("Time is now %s\n", ctime( &TimeNow));
    printf("%s\n", "Calling localtime...");
#endif
 
    TimeNowStruct = localtime ( &TimeNow );
 
#ifdef __DEBUG_DATE
    printf("%s\n", "return from localtime...");
    printf( "asctime gives %s" , asctime ( localtime ( &TimeNow ) ) ) ;
    printf( "TimeNowStruct gives %s" , asctime ( TimeNowStruct) );
#endif
 
    itime[0] = 1900 + TimeNowStruct->tm_year;
    itime[1] = 1 + TimeNowStruct->tm_mon;
    itime[2] = TimeNowStruct->tm_mday;
    itime[3] = 0;
    itime[4] = TimeNowStruct->tm_hour;
    itime[5] = TimeNowStruct->tm_min;
    itime[6] = TimeNowStruct->tm_sec;
    itime[7] = 0;
 
#ifdef __DEBUG_DATE
    printf(" year:  %d\n", itime[0]);
    printf(" month: %d\n", itime[1]);
    printf(" day:   %d\n", itime[2]);
    printf(" hour:  %d\n", itime[4]);
    printf(" min:   %d\n", itime[5]);
    printf(" sec:   %d\n", itime[6]);
#endif
}
