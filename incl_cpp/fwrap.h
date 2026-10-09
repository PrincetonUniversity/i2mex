#ifndef FWRAP_H
#define FWRAP_H

#include <string> 
#include <iostream> 
#include <exception>
#include <map>
#include <list>
#include "f77name.h"
#include "transp_util.h"
#include "transperror.h"

using namespace std ;

//
// ---------------- dimsize ------------------
// helper class for defining different sizes for the logical
// structure of a dimension and the physical storage size.  Use
//   dimsize(n)
// for constructing a dimsize with logical=physical dimensions
// otherwise use
//   dimsize(logical_size, physical_size)
// The fwrap constructor will check that physical_size>=logical_size
//
class dimsize {
public:
  explicit dimsize(size_t a, size_t b=0) : logical(a), physical(b) { if (physical==0) physical=logical ; } ;
  size_t logical ;
  size_t physical ;
} ;

// -------- attributes ---------
//
// generic attribute value storage.  Given one of these, the caller is expected to
// test the type and cast.
//
class fwrap_att {
 public:
  virtual ~fwrap_att() { } ;
  enum att_type { none=0, itype=1, dtype=2, stype=3 } ;  // integer,double,string

  typedef pair<att_type, size_t> attInfo ; // attribute info

  static string att_name(att_type) ;       // convert type to string

  virtual att_type get_type() const=0 ;  // return the type of this attribute
  virtual size_t   get_size() const=0 ;  // return the number of elements in this attribute
  virtual attInfo  get_info() const=0 ;  // get type and size

  virtual fwrap_att* clone() const=0 ;     // duplicate this attribute returning a new object

  virtual ostream& put(ostream& s) const=0 ; // write value to a stream

  static string despace(string) ;          // return the input string without any spaces 
  static string despace0(string) ;         // return the input string without any spaces and throw
                                           // an fwrap_att_error if the result is the empty string

  // ----- error ----
  class fwrap_att_error : public std::exception {
  public:
    fwrap_att_error(string msg) throw() { errmsg=msg ; } ;
    fwrap_att_error(const fwrap_att_error& q) throw() : std::exception(q), errmsg(q.errmsg) {} ;
    fwrap_att_error& operator=(const fwrap_att_error& q) throw();
    virtual ~fwrap_att_error() throw() ;
    virtual const char* what() const throw()  ;

    static fwrap_att_error notScalar(string stype, size_t psize) ;               // return an error when expecting a scalar
    static fwrap_att_error badIndex(string stype, size_t psize, size_t ix) ;     // return an error for an out of range index
    static fwrap_att_error tooSmall(string stype, size_t psize, size_t asize) ;  // return an error when argument array is too small
  private:
    string errmsg ;
  } ;
} ;

inline ostream& operator<<(ostream& s, const fwrap_att& r) {
  return r.put(s) ;
}

ostream& operator<<(ostream& s, fwrap_att::att_type& t) ;

class VarWorkerWriter ;
 
//
// an integer attribute
//
class fwrap_int_att : public fwrap_att {
 public:
  virtual ~fwrap_int_att() ;

  fwrap_int_att(int i=0) ;
  fwrap_int_att(const int* iptr, size_t isize) ;

  virtual att_type   get_type() const ;
  virtual size_t     get_size() const { return psize ; } ;
  virtual attInfo    get_info() const ;  
  
  virtual fwrap_att* clone()    const ;

  int  get_int() const ;  // get the scalar int, error if not scalar
  void set_int(int i) ;   // set as a scalar

  void set_int(const int* iptr, size_t isize) ;     // set as a vector
  size_t fill_int(int* iptr, size_t isize) const ;  // fill the argument with the data, return number copied

  int  getx_int(size_t ix) const ;
  void setx_int(size_t ix, int i) ;

  virtual ostream& put(ostream& s) const ; 
 private:
  size_t psize ;  // size of data
  int*   ptr ;    // data is always allocated

  friend class VarWorkerWriter ;
} ;

//
// a double attribute
//
class fwrap_double_att : public fwrap_att {
 public:
  virtual ~fwrap_double_att() ;

  fwrap_double_att(double i=0.) ;
  fwrap_double_att(const double* iptr, size_t isize) ;

  virtual att_type   get_type() const ;
  virtual size_t     get_size() const { return psize ; } ;
  virtual attInfo    get_info() const ;  
  
  virtual fwrap_att* clone()    const ;

  double get_double() const ;     // get the scalar double, error if not scalar
  void   set_double(double i) ;   // set as a scalar

  void   set_double(const double* iptr, size_t isize) ;   // set as a vector
  size_t fill_double(double* iptr, size_t isize) const ;  // fill the argument with the data, return number copied

  double  getx_double(size_t ix) const ;
  void    setx_double(size_t ix, double i) ;

  virtual ostream& put(ostream& s) const ; 
 private:
  size_t  psize ;  // size of data
  double* ptr ;    // data is always allocated

  friend class VarWorkerWriter ;
} ;

//
// a string attribute
//
class fwrap_string_att : public fwrap_att {
 public:
  virtual ~fwrap_string_att() ;

  fwrap_string_att(string i="") ;
  fwrap_string_att(const string* iptr, size_t isize) ;

  virtual att_type   get_type() const ;
  virtual size_t     get_size() const { return psize ; } ;
  virtual attInfo    get_info() const ;  
  
  size_t  get_maxLength() const ;  // maximum length of all the strings
  
  virtual fwrap_att* clone()    const ;

  string get_string() const ;     // get the scalar string, error if not scalar
  void   set_string(string i) ;   // set as a scalar

  void   set_string(const string* iptr, size_t isize) ;   // set as a vector
  size_t fill_string(string* iptr, size_t isize) const ;  // fill the argument with the data, return number copied

  string  getx_string(size_t ix) const ;
  void    setx_string(size_t ix, string i) ;

  virtual ostream& put(ostream& s) const ; 
 private:
  size_t  psize ;  // size of data
  string* ptr ;    // data is always allocated

  friend class VarWorkerWriter ;
} ;

//
// storage of attributes
// the attribute names will have all the spaces removed and if the result is the empty
// string then an exception will be thrown (fwrap_att::fwrap_att_error)
//
class fwrap_map_att {
 public:
  fwrap_map_att() ;                                // create an empty attribute map
  ~fwrap_map_att() ;                               // destroy the stored attributes
  fwrap_map_att(const fwrap_map_att&) ;            // copy an existing map
  fwrap_map_att& operator=(const fwrap_map_att&) ; // remove all the existing attributes and copy in the ones from the argument

  size_t size() const { return atts.size() ; } ;   // number of attributes
  list<string> keys() const ;                      // return a list of the keys
  void clear() ;                                   // remove all attributes

  fwrap_att*       get(string name) ;              // return an attribute with this name or null if there is no such attribute
  const fwrap_att* get(string name) const ;        // return an attribute with this name or null if there is no such attribute

  // --- scalar methods ---
  void  set(string name, int value) ;              // set this attribute as an integer attribute
  void  set(string name, double value) ;           // set this attribute as a double attribute
  void  set(string name, string value) ;           // set this attribute as a string attribute

  int    get_int(string name) const ;              // return the int    attribute or throw an exception
  double get_double(string name) const ;           // return the double attribute or throw an exception
  string get_string(string name) const ;           // return the string attribute or throw an exception

  int    get_int(string name, int def) const ;        // return the int    attribute or use the default
  double get_double(string name, double def) const ;  // return the double attribute or use the default
  string get_string(string name, string def) const ;  // return the string attribute or use the default

  // --- vector methods ---
  void  set(string name, const int*    value, size_t vsize) ; // set this attribute as a vector integer attribute, data is copied
  void  set(string name, const double* value, size_t vsize) ; // set this attribute as a vector double attribute, data is copied
  void  set(string name, const string* value, size_t vsize) ; // set this attribute as a vector string attribute, data is copied

  void setx(string name, size_t ix, int value) ;      // change the value at an index of an existing vector integer attribute
  void setx(string name, size_t ix, double value) ;   // change the value at an index of an existing vector double attribute
  void setx(string name, size_t ix, string value) ;   // change the value at an index of an existing vector string attribute

  int     getx_int(string name,    size_t ix) const ; // return the value at an index of an existing vector integer attribute
  double  getx_double(string name, size_t ix) const ; // return the value at an index of an existing vector double attribute
  string  getx_string(string name, size_t ix) const ; // return the value at an index of an existing vector string attribute
  
  //
  // in the following the argument size must be at least as large as the source size
  // the size_t returned is the actual size of data which was copied
  //
  size_t fill(string name, int*    value, size_t vsize) const ; // fill the argument with the integer attribute value
  size_t fill(string name, double* value, size_t vsize) const ; // fill the argument with the double attribute value
  size_t fill(string name, string* value, size_t vsize) const ; // fill the argument with the string attribute value

  size_t fill(string name, int*    value, size_t vsize, const int*    def, size_t dsize) const ; // fill the argument with the integer or default attribute value
  size_t fill(string name, double* value, size_t vsize, const double* def, size_t dsize) const ; // fill the argument with the double or default attribute value
  size_t fill(string name, string* value, size_t vsize, const string* def, size_t dsize) const ; // fill the argument with the string or default attribute value

  // --- other ---
  void  remove(string name) ;                      // remove the named attribute

  ostream& write(ostream&, string lead) const ;    // write attribute on each line with "lead" printed first

  // ----- error ----
  class fwrap_map_error : public std::exception {
  public:
    fwrap_map_error(string msg) throw() { errmsg=msg ; } ;
    fwrap_map_error(const fwrap_map_error& q) throw() : std::exception(q), errmsg(q.errmsg) {} ;
    fwrap_map_error& operator=(const fwrap_map_error& q) throw();
    virtual ~fwrap_map_error() throw() ;
    virtual const char* what() const throw()  ;

    static fwrap_map_error notScalar(string name, string stype, size_t psize) ;               // return an error when expecting a scalar
    static fwrap_map_error badIndex(string name, string stype, size_t psize, size_t ix) ;     // return an error for an out of range index
    static fwrap_map_error tooSmall(string name, string stype, size_t psize, size_t asize) ;  // return an error when argument array is too small
  private:
    string errmsg ;
  } ;

 private:
  map<string,fwrap_att*> atts ;

  friend class VarWorkerWriter ;
} ;

ostream& operator<<(ostream& s, const fwrap_map_att&) ;        // simple write

//
// ------------------ fwrap ---------------------------
// Wraps and unwraps a pointer to an array of floats,ints,... .  The advantage
// over valarray<> is that the constructor can take a pointer of
// floats which could be supplied by a fortran subroutine.  Access
// is relative to 0 using fortran ordering.  A pointer supplied in
// a constructor is NOT freed on destruction (i.e. it is borrowed).
//
// Objects of this class have value semantics.  The fwrapPtr class
// wraps fwrap objects to provide reference semantics.
//
// Update:
// operator=() assignment will cause the destruction (if owned) of 
// the current memory and allocation of new memory to hold
// the copied data.  The assign() method is available for simply
// copying into existing memory if the shapes are the same which
// was the original operator=() behaviour.
//
// rank 0 and N fwraps are now allowed and size 0 dimensions are allowed
//

#ifndef NO_GOT_VALARRAY

#include <valarray>
#include <vector>
typedef std::valarray<size_t> Shape ;  // list of dimension sizes -- like fortran
inline void set_vector_int(Shape& s, size_t j) { s=j ; } ;

typedef std::vector<Byte> Bytes ;

#else

#include <vector>
typedef std::vector<size_t> Shape ;
inline void set_vector_int(Shape& s, size_t j) { for(size_t i=0 ; i<s.size() ; i++) s[i]=j ; } ;

typedef std::vector<Byte> Bytes ;

#endif

template <typename FP> class fwrapPtr ;

//
// --------------------------- fwrap ---------------------------
//
template <typename FP> class fwrap {
public:
  // ------- constructors ---------
  fwrap(FP* fp=0) ;                        // create a 0d array (scalar) using FP* as the storage 
                                           // if nonnull

  fwrap(dimsize d, FP* fp=0) ;             // create a 1d array using FP* as the storage if nonnull
  fwrap(size_t  d, FP* fp=0) ;             // create a 1d array using FP* as the storage if nonnull
                                           // with logical_size=physical_size

  fwrap(dimsize d, dimsize e, FP* fp=0) ;  // create a 2d array using FP* as the storage if nonnull
  fwrap(size_t  d, size_t  e, FP* fp=0) ;  // create a 2d array using FP* as the storage if nonnull
                                           // with logical_size=physical_size

  fwrap(const Shape&, FP* fp=0) ;          // create an Nd array using FP* as the storage if nonnull
                                           // with logical = physical
  fwrap(const Shape& logical, const Shape& physical, FP* fp=0) ; // create an Nd array using FP* as 
                                                                 // the storage if nonnull, 
                                                                 // logical.size()==physical.size()
                                                                 // is required

  fwrap(const fwrap<FP>&) ;                // copy constructor, copies attributes also

  ~fwrap() ;                               // destructor -- free memory only if internally created

  // ---- assignment ---
  void assign(const fwrap<FP>&) ;        // only succeeds if arrays agree in actual size.  The whole
                                         // storage is copied but the shape is not changed and need not agree with
                                         // the argument shape  -- the storage location (pointer) is not modified by this call.
                                         // The attributes are not changed by this call.

  void assign(FP) ;                      // copies the scalar FP to the whole storage area

  // -- probe --
  bool   ismine()     const {return myptr ; } ; // true when this fwrap owns the pointer and will delete it
  size_t rank()       const {return prank ; } ; // number of dimensions

  const Shape& shape()  const { return dim ;   } ;  // logical shape of array
  const Shape& pshape() const { return store ; } ;  // physical shape of array
 

  size_t size(size_t) const ;                   // logical size of each dimension -- checked
  size_t size()       const ;                   // total logical size

  size_t pstore(size_t) const ;                 // physical storage size of each dimension -- checked
  size_t pstore()       const ;                 // total physical storage size

  FP&       operator[](const Shape&) ;          // Shape must have a size==rank -- checked
  const FP& operator[](const Shape&) const ;

  // -- dangerous methods --
  // unchecked access -- caller must know rank and limits and use
  // the appropriate number and size of arguments with these methods.
  //
  size_t xsize(size_t i) const { return dim[i] ; } ;        // logical size of each dimension
  size_t xlast_size()    const { return dim[prank-1] ; } ;  // logical size of last dimension

  FP& operator()()                    { return *p ; } ;                   // rank should be 0
  FP& operator()(size_t i)            { return *(p+i) ; } ;               // rank should be 1
  FP& operator()(size_t i, size_t j)  { return *(p+(store[0]*j+i)) ; } ;  // rank should be 2
  FP& operator()(size_t i, size_t j, size_t k)  { 
    return *(p+((k*store[1]+j)*store[0]+i)) ; 
  } ;  // rank should be 3
  FP& operator()(size_t i, size_t j, size_t k, size_t m)  { 
    return *(p+(((m*store[2]+k)*store[1]+j)*store[0]+i)) ; 
  } ;  // rank should be 4
  FP& operator()(size_t i, size_t j, size_t k, size_t m, size_t n)  { 
    return *(p+((((n*store[3]+m)*store[2]+k)*store[1]+j)*store[0]+i)) ; 
  } ;  // rank should be 5
  FP& operator()(size_t i, size_t j, size_t k, size_t m, size_t n, size_t o)  { 
    return *(p+(((((o*store[4]+n)*store[3]+m)*store[2]+k)*store[1]+j)*store[0]+i)) ; 
  } ;  // rank should be 6

  const FP& operator()()                   const { return *p ; } ;                   
  const FP& operator()(size_t i)           const { return *(p+i) ; } ;              // const versions
  const FP& operator()(size_t i, size_t j) const { return *(p+(store[0]*j+i)) ; } ;
  const FP& operator()(size_t i, size_t j, size_t k) const  { 
    return *(p+((k*store[1]+j)*store[0]+i)) ; 
  } ;  // rank should be 3
  const FP& operator()(size_t i, size_t j, size_t k, size_t m) const  { 
    return *(p+(((m*store[2]+k)*store[1]+j)*store[0]+i)) ; 
  } ;  // rank should be 4
  const FP& operator()(size_t i, size_t j, size_t k, size_t m, size_t n) const { 
    return *(p+((((n*store[3]+m)*store[2]+k)*store[1]+j)*store[0]+i)) ; 
  } ;  // rank should be 5
  const FP& operator()(size_t i, size_t j, size_t k, size_t m, size_t n, size_t o) const { 
    return *(p+(((((o*store[4]+n)*store[3]+m)*store[2]+k)*store[1]+j)*store[0]+i)) ; 
  } ;  // rank should be 6

  FP*       fptr()       { return p ; } ; // for use in fortran calls -- do not delete
  const FP* fptr() const { return p ; } ;

  fwrap<FP> last_slice(size_t i) ;      // returns an fwrap built using a pointer to the slice in the
                                        // last dimension with index i ( i<size(rank-1) ).  This can be
                                        // used for working on just the continguous data in the last slice but ...
                                        // WARNING: If the 'this' fwrap is deleted, then the pointer in the 
                                        //          last_slice() fwrap may no be longer valid which could produce
                                        //          a segmentation fault (if your lucky).
 
  // --- attribute storage ---
  // an attribute is an fwrap_att associated with a string and stored in an fwrap_map_att
  //
  size_t att_size() const ;             // return the number of attributes
  fwrap_map_att*       att() ;          // return the object for storing attributes
  const fwrap_map_att* att() const ;    // return the object for storing attributes

private:
  void set_size(size_t i, size_t n) ;   // change the logical size of dimension i to n,
                                        // an exception will be thrown if i>=rank or n<1 
                                        // or n>pstore(i)

  fwrap& operator=(const fwrap<FP>&) ;   // assignment -- will cause allocation of new memory if
                                         // the existing memory is not owned or the total size ends up changing
                                         // avoid using this, I haven't decided whether to keep it

  size_t total_check() const ; // used during construction -- checks dimensions and
                               // returns total physical size

  // -- data --
  bool   myptr ;    // true when the pointer p was created internally and
                    // should be freed on destruction
  size_t prank ;    // number of dimensions

  Shape dim ;       // list of dimensions (logical size of array)
  Shape store ;     // size of storage in each dimension -- store[i]>=dim[i]

  FP* p ;           // float storage

  unsigned count ;  // for fwrapPtr reference counting use only

  mutable fwrap_map_att* atts ;   // attributes for this array if not null

  friend class fwrapPtr<FP> ;
} ;


template<typename FP> ostream& operator<<(ostream&, const fwrap<FP>&) ; // for debugging

//
// copy the logical elements from f to g respecting the physical shape.  If the logical
// size of f is larger then the physical size in g for any dimension and needall=true, an exception will
// be thrown.
//
template<typename FP, typename GP> void fwrapCopy(const fwrap<FP>& f, fwrap<GP>& g, bool needall=false) ;

//
// same as fwrapCopy but selects what values to store in the target based on the
// truth or falseness of the source
//
template<typename FP, typename GP> void fwrapCopy_bool(const fwrap<FP>& f, fwrap<GP>& g, GP truevalue, GP falsevalue, 
						       bool needall=false) ;


template<typename FP> Shape fwrapCharShape(const fwrap<FP>& f) ;   // return the effective logical shape interpreted as string arrays

template<typename FP> size_t fwrapCharLength(const fwrap<FP>& f) ; // return the effective character length.  Only makes sense
                                                                   // for fwrap<char> -- all others should be nonzero
template<typename FP> size_t fwrapCharMaxLength(const fwrap<FP>& f) ; // maximum length of string.  Only makes sense for
                                                                      // fwrap<string> and fwrap<char>
template<typename FP, typename GP> void fwrapCharCopy(const fwrap<FP>& f, fwrap<GP>& g, const Shape& s) ; // copy one element at the index
                                                                   // as interpreted as characters

//
// Like the default fwrapCopy however the first dimension of an fwrap<char> is interpreted like a fortran character string.
// If the source string length is greater than the target it will be truncated independent of needall.
// If the source string length is less than the target, the remainder of the target will be filled with spaces.
//
template<typename FP, typename GP> void fwrapCopy_char(const fwrap<FP>& f, fwrap<GP>& g, bool needall=false) ;

//
// ----------- fwrapPtr exceptions --------------
// exception for when there is a null fwrapPtr
//
class fwrapPtrNull : public std::exception {
public:
  fwrapPtrNull() throw() {} ;
  fwrapPtrNull(const fwrapPtrNull& q) throw() : std::exception(q) {} ;
  fwrapPtrNull& operator=(const fwrapPtrNull& q) throw();
  virtual ~fwrapPtrNull() throw() ;
  virtual const char* what() const throw()  ;
private:
  static const char* errmsg ;
} ;

//
// --------------------- fwrapPtr -----------------------------
// This provides a reference counted wrapper around an fwrap.  This can not be
// used to wrap an existing fwrap object, rather, a new fwrap object is created
// and stored within the fwrapPtr object.
//
template <typename FP> class fwrapPtr {
public:
  // constructors
  static fwrapPtr build(FP* fp=0) ;                        // create a 0d array (scalar) using FP* as the storage 
                                                           // if nonnull

  static fwrapPtr build(dimsize d, FP* fp=0) ;             // create a 1d array using FP* as the storage if nonnull
  static fwrapPtr build(size_t  d, FP* fp=0) ;             // create a 1d array using FP* as the storage if nonnull
                                                           // with logical_size=physical_size

  static fwrapPtr build(dimsize d, dimsize e, FP* fp=0) ;  // create a 2d array using FP* as the storage if nonnull
  static fwrapPtr build(size_t  d, size_t  e, FP* fp=0) ;  // create a 2d array using FP* as the storage if nonnull
                                                           // with logical_size=physical_size

  static fwrapPtr build(const Shape&, FP* fp=0) ;          // create an Nd array using FP* as the storage if nonnull
                                                            // with logical = physical
  static fwrapPtr build(const Shape& logical, const Shape& physical, FP* fp=0) ; // create an Nd array using FP* as 
                                                                          // the storage if nonnull, 
                                                                          // logical.size()==physical.size()
                                                                          // is required
  static fwrapPtr build(const fwrap<FP>& p) ;              // copy the argument to a new array

  // pointer stuff
  fwrap<FP>* operator->()             { if (q==0) throw fwrapPtrNull() ; return q  ; } ;
  fwrap<FP>& operator*()              { if (q==0) throw fwrapPtrNull() ; return *q ; } ;
  const fwrap<FP>* operator->() const { if (q==0) throw fwrapPtrNull() ; return q  ; } ;
  const fwrap<FP>& operator*()  const { if (q==0) throw fwrapPtrNull() ; return *q ; } ;
  fwrap<FP>* pointer()                { if (q==0) throw fwrapPtrNull() ; return q  ; } ;

  // object management
  bool isnull() const { return q==0 ; } ;
  unsigned count() const { return (q!=0) ? q->count : 0 ; } ;
  
  fwrapPtr<FP>() : q(0) {} ;
  fwrapPtr<FP>(const fwrapPtr<FP>& p) : q(p.q) { if (q!=0) ++q->count ; } ;
  fwrapPtr<FP>& operator=(const fwrapPtr<FP>& p) {
    if (p.q!=0) ++p.q->count ;
    if (q!=0 && --q->count==0) delete q ;
    q=p.q ;
    return *this ;
  } ;
  ~fwrapPtr<FP>() { if (q!=0 && --q->count==0) delete q ; } ;
private:
  fwrapPtr(fwrap<FP>* p) : q(p) { ++q->count ; } ; // null pointer not allowed

  fwrap<FP>* q ; 
} ;

template<typename FP> ostream& operator<<(ostream&, const fwrapPtr<FP>&) ; // for debugging

//
// ---------- insame_shape -----------
// used in React to compare two Shapes for equality
//
inline bool insame_shape(const Shape& a, const Shape& b) {
  if (a.size() == 1) return (b.size()==1 && a[0]==b[0]) ;  // quick return for 1d
  if (a.size() != b.size()) return false ;
  for (size_t i=0 ; i<a.size() ; i++) {
    if (a[i]!=b[i]) return false ;
  } ;
  return true ;
} ;

//
// ---------- insame_last_shape -----------
// used in React to compare the shape of the first
// argument to the shape of the second ignoring the
// last dimension (as if using b.last_slice())
inline bool insame_last_shape(const Shape& a, const Shape& b) {
  if (b.size() == 0) return false ;                        // guard against 0d b
  if (a.size() == 1) return (b.size()==2 && a[0]==b[0]) ;  // quick return for 1d
  if (a.size() != b.size()-1) return false ;
  for (size_t i=0 ; i<a.size() ; i++) {
    if (a[i]!=b[i]) return false ;
  } ;
  return true ;
} ;


//
// --------- typedefs for fwrap ----------
// first set are fortran types, the second are C++ only types
//
typedef fwrap<lp_32>   alp_32 ;       // array of logical*4
typedef fwrap<ip_32>   aip_32 ;       // array of ints
typedef fwrap<fp_64>   afp_64 ;       // array of doubles
typedef fwrap<fp_32>   afp_32 ;       // array of singles
typedef fwrap<cp_char> acp_char ;     // array of chars, first dimension is the length

typedef fwrap<bool>    alp_bool ;     // array of C++ bool
typedef fwrap<Bytes>   abp_bytes ;    // array of byte arrays
typedef fwrap<string>  acp_string ;   // array of C++ strings

//
// --------- typedefs for fwrapPtr ----------
// first set are fortran types, the second are C++ only types
//
typedef fwrapPtr<lp_32>   alpPtr_32 ;      // array of logical*4
typedef fwrapPtr<ip_32>   aipPtr_32 ;      // array of ints
typedef fwrapPtr<fp_64>   afpPtr_64 ;      // array of doubles
typedef fwrapPtr<fp_32>   afpPtr_32 ;      // array of singles
typedef fwrapPtr<cp_char> acpPtr_char ;    // array of chars, first dimension is the length

typedef fwrapPtr<bool>    alpPtr_bool ;    // array of C++ bool
typedef fwrapPtr<Bytes>   abpPtr_bytes ;   // array of byte arrays
typedef fwrapPtr<string>  acpPtr_string ;  // array of C++ strings

typedef fwrap<fp_64> ap_64 ;   // array of doubles -- deprecated name
typedef fwrap<fp_32> ap_32 ;   // array of singles -- deprecated name


//
// ---------------- opaque C access -----------------
// these are raw unchecked operations, 
// the indexing has fortran order starting from 0 for cw_* functions and
// from 1 for fw_* functions.
//
extern "C" int   F77NAME(fw_ip32rank)(aip_32** self) ;                                  // return the rank
extern "C" void  F77NAME(fw_ip32shape)(aip_32** self, int* ishape) ;                    // return the logical shape in the array ishape
extern "C" ip_32 F77NAME(cw_ip32get0)(aip_32** self) ;                                  // get a rank 0 value
extern "C" ip_32 F77NAME(cw_ip32get1)(aip_32** self, int*) ;                            // get a rank 1 value
extern "C" ip_32 F77NAME(cw_ip32get2)(aip_32** self, int*, int*) ;                      // get a rank 2 value
extern "C" ip_32 F77NAME(cw_ip32get3)(aip_32** self, int*, int*, int*) ;                // get a rank 3 value
extern "C" ip_32 F77NAME(cw_ip32get4)(aip_32** self, int*, int*, int*, int*) ;          // get a rank 4 value
extern "C" void  F77NAME(cw_ip32set0)(aip_32** self, ip_32*) ;                          // set a rank 0 value
extern "C" void  F77NAME(cw_ip32set1)(aip_32** self, int*, ip_32*) ;                    // set a rank 1 value
extern "C" void  F77NAME(cw_ip32set2)(aip_32** self, int*, int*, ip_32*) ;              // set a rank 2 value
extern "C" void  F77NAME(cw_ip32set3)(aip_32** self, int*, int*, int*, ip_32*) ;        // set a rank 3 value
extern "C" void  F77NAME(cw_ip32set4)(aip_32** self, int*, int*, int*, int*, ip_32*) ;  // set a rank 4 value
extern "C" ip_32 F77NAME(fw_ip32get0)(aip_32** self) ;                                  // get a rank 0 value
extern "C" ip_32 F77NAME(fw_ip32get1)(aip_32** self, int*) ;                            // get a rank 1 value
extern "C" ip_32 F77NAME(fw_ip32get2)(aip_32** self, int*, int*) ;                      // get a rank 2 value
extern "C" ip_32 F77NAME(fw_ip32get3)(aip_32** self, int*, int*, int*) ;                // get a rank 3 value
extern "C" ip_32 F77NAME(fw_ip32get4)(aip_32** self, int*, int*, int*, int*) ;          // get a rank 4 value
extern "C" void  F77NAME(fw_ip32set0)(aip_32** self, ip_32*) ;                          // set a rank 0 value
extern "C" void  F77NAME(fw_ip32set1)(aip_32** self, int*, ip_32*) ;                    // set a rank 1 value
extern "C" void  F77NAME(fw_ip32set2)(aip_32** self, int*, int*, ip_32*) ;              // set a rank 2 value
extern "C" void  F77NAME(fw_ip32set3)(aip_32** self, int*, int*, int*, ip_32*) ;        // set a rank 3 value
extern "C" void  F77NAME(fw_ip32set4)(aip_32** self, int*, int*, int*, int*, ip_32*) ;  // set a rank 4 value

extern "C" int   F77NAME(fw_fp32rank)(afp_32** self) ;                                  // return the rank
extern "C" void  F77NAME(fw_fp32shape)(afp_32** self, int* ishape) ;                    // return the logical shape in the array ishape
extern "C" fp_32 F77NAME(cw_fp32get0)(afp_32** self) ;                                  // get a rank 0 value
extern "C" fp_32 F77NAME(cw_fp32get1)(afp_32** self, int*) ;                            // get a rank 1 value
extern "C" fp_32 F77NAME(cw_fp32get2)(afp_32** self, int*, int*) ;                      // get a rank 2 value
extern "C" fp_32 F77NAME(cw_fp32get3)(afp_32** self, int*, int*, int*) ;                // get a rank 3 value
extern "C" fp_32 F77NAME(cw_fp32get4)(afp_32** self, int*, int*, int*, int*) ;          // get a rank 4 value
extern "C" fp_32 F77NAME(cw_fp32get5)(afp_32** self, int*, int*, int*, int*,int*) ;     // get a rank 5 value
extern "C" void  F77NAME(cw_fp32set0)(afp_32** self, fp_32*) ;                          // set a rank 0 value
extern "C" void  F77NAME(cw_fp32set1)(afp_32** self, int*, fp_32*) ;                    // set a rank 1 value
extern "C" void  F77NAME(cw_fp32set2)(afp_32** self, int*, int*, fp_32*) ;              // set a rank 2 value
extern "C" void  F77NAME(cw_fp32set3)(afp_32** self, int*, int*, int*, fp_32*) ;        // set a rank 3 value
extern "C" void  F77NAME(cw_fp32set4)(afp_32** self, int*, int*, int*, int*, fp_32*) ;  // set a rank 4 value
extern "C" void  F77NAME(cw_fp32set5)(afp_32** self, int*, int*, int*, int*, int*, fp_32*) ;  // set a rank 5 value
extern "C" fp_32 F77NAME(fw_fp32get0)(afp_32** self) ;                                  // get a rank 0 value
extern "C" fp_32 F77NAME(fw_fp32get1)(afp_32** self, int*) ;                            // get a rank 1 value
extern "C" fp_32 F77NAME(fw_fp32get2)(afp_32** self, int*, int*) ;                      // get a rank 2 value
extern "C" fp_32 F77NAME(fw_fp32get3)(afp_32** self, int*, int*, int*) ;                // get a rank 3 value
extern "C" fp_32 F77NAME(fw_fp32get4)(afp_32** self, int*, int*, int*, int*) ;          // get a rank 4 value
extern "C" fp_32 F77NAME(fw_fp32get5)(afp_32** self, int*, int*, int*, int*, int*) ;    // get a rank 5 value
extern "C" void  F77NAME(fw_fp32set0)(afp_32** self, fp_32*) ;                          // set a rank 0 value
extern "C" void  F77NAME(fw_fp32set1)(afp_32** self, int*, fp_32*) ;                    // set a rank 1 value
extern "C" void  F77NAME(fw_fp32set2)(afp_32** self, int*, int*, fp_32*) ;              // set a rank 2 value
extern "C" void  F77NAME(fw_fp32set3)(afp_32** self, int*, int*, int*, fp_32*) ;        // set a rank 3 value
extern "C" void  F77NAME(fw_fp32set4)(afp_32** self, int*, int*, int*, int*, fp_32*) ;  // set a rank 4 value
extern "C" void  F77NAME(fw_fp32set5)(afp_32** self, int*, int*, int*, int*, int*, fp_32*) ;  // set a rank 5 value

extern "C" int   F77NAME(fw_fp64rank)(afp_64** self) ;                                  // return the rank
extern "C" void  F77NAME(fw_fp64shape)(afp_64** self, int* ishape) ;                    // return the logical shape in the array ishape
extern "C" fp_64 F77NAME(cw_fp64get0)(afp_64** self) ;                                  // get a rank 0 value
extern "C" fp_64 F77NAME(cw_fp64get1)(afp_64** self, int*) ;                            // get a rank 1 value
extern "C" fp_64 F77NAME(cw_fp64get2)(afp_64** self, int*, int*) ;                      // get a rank 2 value
extern "C" fp_64 F77NAME(cw_fp64get3)(afp_64** self, int*, int*, int*) ;                // get a rank 3 value
extern "C" fp_64 F77NAME(cw_fp64get4)(afp_64** self, int*, int*, int*, int*) ;          // get a rank 4 value
extern "C" fp_64 F77NAME(cw_fp64get5)(afp_64** self, int*, int*, int*, int*, int*) ;    // get a rank 5 value
extern "C" void  F77NAME(cw_fp64set0)(afp_64** self, fp_64*) ;                          // set a rank 0 value
extern "C" void  F77NAME(cw_fp64set1)(afp_64** self, int*, fp_64*) ;                    // set a rank 1 value
extern "C" void  F77NAME(cw_fp64set2)(afp_64** self, int*, int*, fp_64*) ;              // set a rank 2 value
extern "C" void  F77NAME(cw_fp64set3)(afp_64** self, int*, int*, int*, fp_64*) ;        // set a rank 3 value
extern "C" void  F77NAME(cw_fp64set4)(afp_64** self, int*, int*, int*, int*, fp_64*) ;  // set a rank 4 value
extern "C" void  F77NAME(cw_fp64set5)(afp_64** self, int*, int*, int*, int*, int*, fp_64*) ;  // set a rank 5 value
extern "C" fp_64 F77NAME(fw_fp64get0)(afp_64** self) ;                                  // get a rank 0 value
extern "C" fp_64 F77NAME(fw_fp64get1)(afp_64** self, int*) ;                            // get a rank 1 value
extern "C" fp_64 F77NAME(fw_fp64get2)(afp_64** self, int*, int*) ;                      // get a rank 2 value
extern "C" fp_64 F77NAME(fw_fp64get3)(afp_64** self, int*, int*, int*) ;                // get a rank 3 value
extern "C" fp_64 F77NAME(fw_fp64get4)(afp_64** self, int*, int*, int*, int*) ;          // get a rank 4 value
extern "C" fp_64 F77NAME(fw_fp64get5)(afp_64** self, int*, int*, int*, int*, int*) ;    // get a rank 5 value
extern "C" void  F77NAME(fw_fp64set0)(afp_64** self, fp_64*) ;                          // set a rank 0 value
extern "C" void  F77NAME(fw_fp64set1)(afp_64** self, int*, fp_64*) ;                    // set a rank 1 value
extern "C" void  F77NAME(fw_fp64set2)(afp_64** self, int*, int*, fp_64*) ;              // set a rank 2 value
extern "C" void  F77NAME(fw_fp64set3)(afp_64** self, int*, int*, int*, fp_64*) ;        // set a rank 3 value
extern "C" void  F77NAME(fw_fp64set4)(afp_64** self, int*, int*, int*, int*, fp_64*) ;  // set a rank 4 value
extern "C" void  F77NAME(fw_fp64set5)(afp_64** self, int*, int*, int*, int*, int*, fp_64*) ;  // set a rank 5 value

extern "C" int   F77NAME(fw_lp32rank)(alp_32** self) ;                                  // return the rank
extern "C" void  F77NAME(fw_lp32shape)(alp_32** self, int* ishape) ;                    // return the logical shape in the array ishape
extern "C" lp_32 F77NAME(cw_lp32get0)(alp_32** self) ;                                  // get a rank 0 value
extern "C" lp_32 F77NAME(cw_lp32get1)(alp_32** self, int*) ;                            // get a rank 1 value
extern "C" lp_32 F77NAME(cw_lp32get2)(alp_32** self, int*, int*) ;                      // get a rank 2 value
extern "C" lp_32 F77NAME(cw_lp32get3)(alp_32** self, int*, int*, int*) ;                // get a rank 3 value
extern "C" lp_32 F77NAME(cw_lp32get4)(alp_32** self, int*, int*, int*, int*) ;          // get a rank 4 value
extern "C" void  F77NAME(cw_lp32set0)(alp_32** self, lp_32*) ;                          // set a rank 0 value
extern "C" void  F77NAME(cw_lp32set1)(alp_32** self, int*, lp_32*) ;                    // set a rank 1 value
extern "C" void  F77NAME(cw_lp32set2)(alp_32** self, int*, int*, lp_32*) ;              // set a rank 2 value
extern "C" void  F77NAME(cw_lp32set3)(alp_32** self, int*, int*, int*, lp_32*) ;        // set a rank 3 value
extern "C" void  F77NAME(cw_lp32set4)(alp_32** self, int*, int*, int*, int*, lp_32*) ;  // set a rank 4 value
extern "C" lp_32 F77NAME(fw_lp32get0)(alp_32** self) ;                                  // get a rank 0 value
extern "C" lp_32 F77NAME(fw_lp32get1)(alp_32** self, int*) ;                            // get a rank 1 value
extern "C" lp_32 F77NAME(fw_lp32get2)(alp_32** self, int*, int*) ;                      // get a rank 2 value
extern "C" lp_32 F77NAME(fw_lp32get3)(alp_32** self, int*, int*, int*) ;                // get a rank 3 value
extern "C" lp_32 F77NAME(fw_lp32get4)(alp_32** self, int*, int*, int*, int*) ;          // get a rank 4 value
extern "C" void  F77NAME(fw_lp32set0)(alp_32** self, lp_32*) ;                          // set a rank 0 value
extern "C" void  F77NAME(fw_lp32set1)(alp_32** self, int*, lp_32*) ;                    // set a rank 1 value
extern "C" void  F77NAME(fw_lp32set2)(alp_32** self, int*, int*, lp_32*) ;              // set a rank 2 value
extern "C" void  F77NAME(fw_lp32set3)(alp_32** self, int*, int*, int*, lp_32*) ;        // set a rank 3 value
extern "C" void  F77NAME(fw_lp32set4)(alp_32** self, int*, int*, int*, int*, lp_32*) ;  // set a rank 4 value

//
// these transfers the characters in the first dimension as a string
//
void fwrapCharCharCopy(size_t lsrc, const cp_char* src, size_t ldest, cp_char* dest) ;  // string like copy from src to destination

extern "C" int   F77NAME(fw_charrank)(acp_char** self) ;                                  // return the rank includes length in first dimension
extern "C" void  F77NAME(fw_charshape)(acp_char** self, int* ishape) ;                    // return the logical shape in the array ishape includes length
extern "C" void  F77NAME(cw_charget0)(acp_char** self, const int* ilen, cp_char* c) ;                          // get a rank 0 value
extern "C" void  F77NAME(cw_charget1)(acp_char** self, int*, const int* ilen, cp_char* c) ;                    // get a rank 1 value
extern "C" void  F77NAME(cw_charget2)(acp_char** self, int*, int*, const int* ilen, cp_char* c) ;              // get a rank 2 value
extern "C" void  F77NAME(cw_charget3)(acp_char** self, int*, int*, int*, const int* ilen, cp_char* c) ;        // get a rank 3 value
extern "C" void  F77NAME(cw_charget4)(acp_char** self, int*, int*, int*, int*, const int* ilen, cp_char* c) ;  // get a rank 4 value
extern "C" void  F77NAME(cw_charset0)(acp_char** self, const int* ilen, cp_char* c) ;                          // set a rank 0 value
extern "C" void  F77NAME(cw_charset1)(acp_char** self, int*, const int* ilen, cp_char* c) ;                    // set a rank 1 value
extern "C" void  F77NAME(cw_charset2)(acp_char** self, int*, int*, const int* ilen, cp_char* c) ;              // set a rank 2 value
extern "C" void  F77NAME(cw_charset3)(acp_char** self, int*, int*, int*, const int* ilen, cp_char* c) ;        // set a rank 3 value
extern "C" void  F77NAME(cw_charset4)(acp_char** self, int*, int*, int*, int*, const int* ilen, cp_char* c) ;  // set a rank 4 value
extern "C" void  F77NAME(fw_charget0)(acp_char** self, const int* ilen, cp_char* c) ;                          // get a rank 0 value
extern "C" void  F77NAME(fw_charget1)(acp_char** self, int*, const int* ilen, cp_char* c) ;                    // get a rank 1 value
extern "C" void  F77NAME(fw_charget2)(acp_char** self, int*, int*, const int* ilen, cp_char* c) ;              // get a rank 2 value
extern "C" void  F77NAME(fw_charget3)(acp_char** self, int*, int*, int*, const int* ilen, cp_char* c) ;        // get a rank 3 value
extern "C" void  F77NAME(fw_charget4)(acp_char** self, int*, int*, int*, int*, const int* ilen, cp_char* c) ;  // get a rank 4 value
extern "C" void  F77NAME(fw_charset0)(acp_char** self, const int* ilen, cp_char* c) ;                          // set a rank 0 value
extern "C" void  F77NAME(fw_charset1)(acp_char** self, int*, const int* ilen, cp_char* c) ;                    // set a rank 1 value
extern "C" void  F77NAME(fw_charset2)(acp_char** self, int*, int*, const int* ilen, cp_char* c) ;              // set a rank 2 value
extern "C" void  F77NAME(fw_charset3)(acp_char** self, int*, int*, int*, const int* ilen, cp_char* c) ;        // set a rank 3 value
extern "C" void  F77NAME(fw_charset4)(acp_char** self, int*, int*, int*, int*, const int* ilen, cp_char* c) ;  // set a rank 4 value

#include "fwrap_temp.h"        // stupid compiler
//-----------------------------------------------------------------------------

#endif

