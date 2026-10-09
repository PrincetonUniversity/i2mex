// DO NOT MODIFY, GENERATED CODE!!!


#ifndef ncbox_H
#define ncbox_H

#include "cvarcontainer.h"
#include "f77name.h"

//
// create ncbox.h,ncbox.cpp,ncbox.f90 by executing python varcontainer.py -n ncbox
//
// see varcontainer.h for documentation on the following methods represented in C 
//
namespace ncbox_space {
  //
  // return the unique instance of the ncbox object pointer
  //
  extern "C" void F77NAME(ncbox_cvc_pointer_fetch)(cVarContainer** object) ;

  //
  // return the unique instance of the ncbox object pointer for C++
  //
  extern "C" VarContainer* ncbox_object() ;
} ;

//
// a convenience function to call from totalview's evaluate window
// calls cvc_local() and will print out the variables contents if iout!=0.
//
extern "C" void c_ncbox_local(const char* cname, int iout) ;

//
// a convenience function for showing the container's variables from totalview
//
extern "C" void c_ncbox_show() ;

#endif
