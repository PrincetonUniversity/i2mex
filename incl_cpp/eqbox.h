// DO NOT MODIFY, GENERATED CODE!!!


#ifndef eqbox_H
#define eqbox_H

#include "cvarcontainer.h"
#include "f77name.h"

//
// create eqbox.h,eqbox.cpp,eqbox.f90 by executing python varcontainer.py -n eqbox
//
// see varcontainer.h for documentation on the following methods represented in C 
//
namespace eqbox_space {
  //
  // return the unique instance of the eqbox object pointer
  //
  extern "C" void F77NAME(eqbox_cvc_pointer_fetch)(cVarContainer** object) ;

  //
  // return the unique instance of the eqbox object pointer for C++
  //
  extern "C" VarContainer* eqbox_object() ;
} ;

//
// a convenience function to call from totalview's evaluate window
// calls cvc_local() and will print out the variables contents if iout!=0.
//
extern "C" void c_eqbox_local(const char* cname, int iout) ;

//
// a convenience function for showing the container's variables from totalview
//
extern "C" void c_eqbox_show() ;

#endif
