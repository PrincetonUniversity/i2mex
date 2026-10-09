#ifndef GATEKEEPER_H
#define GATEKEEPER_H

#include <vector>
#include <list>
#include <string>
#include <utility>
#include <iostream>

using namespace std ;

//
// ------------------------ GateKeeperError -----------------------
// class for GateKeeper errors -- these are always bad
//
class GateKeeperError : public exception {
 public:
  GateKeeperError(string container, string info) throw() ;
  GateKeeperError(const GateKeeperError&) throw() ;
  GateKeeperError& operator=(const GateKeeperError&) throw() ;
  virtual ~GateKeeperError() throw() ;
  
  virtual const char*   what()  const throw() ;
  virtual const string  whats() const throw() ;
 protected:
  std::string description ;   // string for use in error message
} ;
 
const size_t GATEKEEPER_MAX_RANK=7 ;   // do not allow ranks larger than this in the GateKeepers
  
//
// ------------------------- GateKeeper -------------------------
// class for holding onto pairs of (string,fwrapPtr) objects which can be added,removed
// and quickly indexed.  You will get an error if you use an index which does not point
// to a filled slot.
//  
template<class FP> class GateKeeper {
 public:
  GateKeeper(string name, size_t n) ;              // construct with a default number of entries and a name for error messages
  ~GateKeeper() ;

  FP getPoint(size_t i) const ;                    // retrieve the FP at the index or an exception if empty
  FP getPoint(size_t i, string checkname) const ;  // retrieve the FP at the index or an exception if empty
                                                   // or if the name does not agree with the second argument
  string getName(size_t i) const ;                 // retrieve the name at the index or an exception if empty

  size_t add(string name, FP point) ;              // add a new pair and return the slot index
  void   remove(size_t i) ;                        // remove the pair at this index, exception if empty
  void   remove(size_t i, string checkname) ;      // remove the pair at this index, exception if empty or the name 
                                                   // does not agree with the argument
  void   clear() ;                                 // remove all the pairs

  size_t size()     const ;           // number of used slots
  size_t capacity() const ;           // number of open and used slots, grows dynamically
  void   check()    const ;           // consistency check

  void   write(ostream&) const ;      // write to the stream

  list<size_t> indices() const ;      // list of all the indices in use

 private:
  typedef pair<string,FP> value ;   // what we are holding

  GateKeeper(const GateKeeper&) ;
  GateKeeper& operator=(const GateKeeper&) ;

  void resize(size_t nsize) ;       // resize the containers
  
  // --- data ---
  string cname ;                    // name for error messages
  size_t nstart ;                   // initial size

  vector<value*> keys ;             // hold pointers to values
  list<size_t>   ghosts ;           // holds open indices
} ;

template<class FP> ostream& operator<<(ostream& s, const GateKeeper<FP>& f) {
  f.write(s) ;
  return s ;
} ;


#include "gatekeeper_temp.h"

#endif
