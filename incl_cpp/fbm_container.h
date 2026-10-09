#ifndef FBM_CONTAINER
#define FBM_CONTAINER_H

#include <string>
#include "f77name.h"

//
// these subroutines define a function for returning the fast ion distribution function
// which can be changed during program execution.
//


//
// function type for returning the fast ion distribution function.
//  jflr      = 0 if point is at gc, 1 if at flr
//  beam_zone = beam zone index of point
//  point     = x,y,z position in plasma either at gc or flr depending on jflr
//  xksid     = v_parallel/v relative to Ip
//  engy      = energy in keV
//  iznum     = atomic number of fast ion species
//  aweight   = atomic weight of fast ion species
//  dist      = returned distribution value
//  ier       = returned nonzero on error
//
typedef void (*fbm_container)(int* jflr, int* beam_zone, double* point, double* xksid, 
			      double* engy, int* iznum, double* aweight, 
			      double* dist, int* ier) ;

//
// set the default function to be called by fbm_container_call
//
void fbm_container_register(fbm_container foo) ;

//
// calls the function stored with fbm_container_register to return the fast ion distribution function
//  jflr      = 0 if point is at gc, !=0 if at flr
//  beam_zone = beam zone index
//  point     = x,y,z position in plasma either at gc or flr depending on nlfbmflr
//  xksid     = v_parallel/v relative to Ip
//  engy      = energy in keV
//  iznum     = atomic number of fast ion species
//  aweight   = atomic weight of fast ion species
//  dist      = returned distribution value
//  ier       = returned nonzero on error
//
extern "C" void F77NAME(fbm_container_call)(int* jflr, int* beam_zone, double* point, double* xksid, 
					    double* engy, int* iznum, double* aweight, 
					    double* dist, int* ier) ;
#endif
