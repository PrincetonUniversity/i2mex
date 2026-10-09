#ifndef TRANSPERROR_H
#define TRANSPERROR_H

#include <exception>
#include <string> 

//
// ----------------------- TranspError ----------------------
// exception class for fatal errors
//
class TranspError : public std::exception {
public:
  TranspError(std::string info) throw() ;
  TranspError(const TranspError&) throw() ;
  TranspError& operator=(const TranspError&) throw() ;
  virtual ~TranspError() throw() ;
  
  virtual const char*        what()  const throw() ;
  virtual const std::string  whats() const throw() ;
protected:
  std::string description ;   // string for use in error message
} ;

//
// Executing this function will write out the string to cerr
// then call bad_exit()
//
void ExitTranspError(std::string info) ;

extern "C"  void F77NAME(bad_exit)() ;  // from portlib

//
// define exceptions as  ... TRANSP_THROW("?preact::module_name: Informative text") ;
// so that exception handling can be turned on or off by defining/undefining THROW_BAD
// preprocessor symbol.
//
// #define THROW_BAD -> bad_exit() will be called instead of throwing an exception
//
// #undef  THROW_BAD -> a TranspError() exception will be thrown which can be caught
//
#ifdef __GNUC__
#if __GNUC__ == 2
#define THROW_BAD  // gcc 2.96 sometimes hangs when thrown an exception
#endif
#endif

#ifdef THROW_BAD
#define TRANSP_THROW(s) ExitTranspError(s)
#else
#define TRANSP_THROW(s) throw TranspError(s)
#endif

//
// for wrapping C++ code called by fortran code which does not catch C++
// exceptions we write out the error message and let the fortran crash.
//   try {
//        ... code which could crash
//   } FCATCH() ;
//
#define FCATCH() catch(std::exception& e) { cerr<<e.what()<<endl ; F77NAME(bad_exit)() ;}

#endif

