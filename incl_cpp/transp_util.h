#ifndef TRANSP_UTIL_H
#define TRANSP_UTIL_H
//
// -------------------- typedefs ----------------------------
//

typedef int  ip_32 ;
typedef int  lp_32 ;    // logical*4, LOGICAL required to have same size as INTEGER
typedef char cp_char ;

typedef unsigned char Byte;

const lp_32 lp_32_TRUE  = -1 ; // logical*4=.TRUE. -- not trustworthy, see portlib/get_fortran_const()
const lp_32 lp_32_FALSE = 0 ;  // logical*4=.FALSE.

typedef double fp_64;
typedef float  fp_32;

#endif
