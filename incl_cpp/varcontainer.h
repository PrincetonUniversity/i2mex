#ifndef VARCONTAINER_H
#define VARCONTAINER_H

#include <string>
#include <map>
#include <list>
#include "transp_util.h"
#include "fwrap.h"
#include "gatekeeper.h"

//
// When the c++ compiler meets the new 2003 standard, unordered_map<> will be
// selected for the mapping.  Otherwise the deprecated hash_map<> will be used if using g++
// and finally use map<> if all else fails.  To use unordered_map<> with g++ before the standard
// is implemented, compile with -std=c++11 or '-std=gnu++0x' or '-std=c++0x'.  In TRANSP, this can be done
// by placing the flags in the GXXFLAGS envrionment variable.
//
#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
#define HAVE_UNOMAP
#endif

#ifdef HAVE_UNOMAP

// -- preferred standard hash map --
//#pragma message "VarContainer using unordered_map"
#include <unordered_map>
#define VARCONTAINER_UNOMAP

#else

#if defined(__GNUC__) && __GNUC__ >3

// -- original gcc hash map --
#pragma message "VarContainer using ext/hash_map"
#include <ext/hash_map>
#define VARCONTAINER_HASHMAP

using namespace __gnu_cxx;

// the hash function for looking up names
struct hash_string {
  size_t operator() (const string& x) const { return hash<const char*>()(x.c_str()) ; } 
} ;

#else

// -- if nothing else available --
#pragma message "VarContainer using regular old map"
#include <map>
#define VARCONTAINER_MAP

#endif
#endif

using namespace std ;

typedef GateKeeper<aipPtr_32>    gkInt     ;   // integer*4
typedef GateKeeper<afpPtr_32>    gkFloat   ;   // real*4
typedef GateKeeper<afpPtr_64>    gkDouble  ;   // real*8
typedef GateKeeper<alpPtr_32>    gkLogical ;   // logical*4
typedef GateKeeper<acpPtr_char>  gkChar ;      // character*(*) 
typedef GateKeeper<abpPtr_bytes> gkBytes ;     // valarray<Byte>

//
// ------------------ VarWriter, VarReader ----------
// Interfaces for objects which can send/receive a stream of bytes
// for use in serializing a VarContainer.  The VarContainer::write() and VarContainer::read() methods
// will call begin() first and either done() or error() last.
//
class VarWriter {
public:
    virtual ~VarWriter() { } ;
    
    virtual void begin() = 0;  // called when starting the serialization of a VarContainer
    
    virtual void send(const Byte* p, size_t nbytes) = 0 ;  // write bytes to stream, might throw an exception
                                                           // this will be called between begin() and done()/error()
    
    virtual void done() = 0 ;  // called when successfully finished with serialization, no more bytes will be written.
                               // This should not cause an exception
                               
    virtual void error(string msg) = 0 ;  // will be called once by the VarContainer::send() method
                                          // if an exception is thrown or will be thrown. No more bytes will 
                                          // be written after this is called. This should not cause an exception.
} ;

class VarReader {
public:
    virtual ~VarReader() { } ;
    
    virtual void begin() = 0;  // called when starting the deserialization of a VarContainer
                                                                                        
    virtual void receive(Byte* p, size_t nbytes) = 0 ;  // read nbytes into the pointer, might throw an exception                                                       

    virtual void done() = 0 ;  // called when finished with serialization, no more bytes will be read.
                               // This should not cause an exception  

    virtual void error(string msg) = 0 ;  // called once by the VarContainer::read() method if
                                          // an error or exception is thrown.  No more bytes will be read.
                                          // This should not cause an exception
} ;


//
// ----------------------- VarContainer --------------------
// a place for storing and looking up variables.  The name must be unique
// among all the variable types.  The scalar methods only work with rank 0 variables
// and transfers data as values.  A scalar add() creates a new scalar fwrapPtr for the 
// scalar storage.  The array methods work with scalar and array data
// and transfers or copies the array references as an fwrap or fwrapPtr,
//
//   get<type>(string name, fwrapPtr def)          -> retrieves the fwrapPtr or returns the default if name not found,
//                                                    throws an exception if the name is not found and the 
//                                                    default fwrapPtr is null.  The stored type must be the same
//                                                    as <type>.
//   get(string name, fwrap&   out, fwrapPtr def)
//   get(string name, fwrapPtr out, fwrapPtr def)  -> copies the stored fwrap to the out argument or copies def
//                                                    if name not found, throws an exception if the def 
//                                                    fwrapPtr is null.  A stored integer can be copied to
//                                                    an integer or floating point "out".  A stored floating point
//                                                    can be copied to a floating point "out".
//   add(string name, fwrapPtr fptr)               -> stores this fwrapPtr in the container, throws an exception if
//                                                    name already exists
//   set(string name, const fwrap& sptr)
//   set(string name, fwrapPtr sptr)  -> copies the data from sptr to the stored fptr using fwrapCopy(*sptr,*fptr),
//                                       expect an excpetion if the name is not found or if the ranks are 
//                                       different or if any of sptr dimensions are larger then fptr's dimensions.
//                                       An integer argument can set a stored integer or floating point.  A 
//                                       floating point can set a stored floating point.
//
//   bool replace(string name, int ioption, fwrapPtr fptr) -> an add() will fail if the variable already exists and
//                                       a set will fail if the name does not exist.  This method combines the 
//                                       add()/remove()/set()/get() for the common operation of setting a variable
//                                       independent of whether it exists or not.  The 'ioption' argument guides
//                                       the operation.  The function returns true if the fwrapPtr was added to the
//                                       container.  The temp_att() attribute is set as in the add method.
//                                       A replace of fwrapPtr<bool> and fwrapPtr<string> will always require a copy.
//                                         ioption =  0 -> equivalent to remove() & add(), always returns true
//                                                 =  1 -> if the variable exists with the correct type,rank,shape
//                                                         then a copy of the data is performed otherwise the pointer
//                                                         is added to the container possibly replacing an existing pointer.
//                                                 =  2 -> same as 1 but the type can be a copyable_type (float<->double)
//                                                 =  3 -> same as 2 but the contained variable need only agree in copyable_type and
//                                                         rank and have a large enough shape in each dimension to hold
//                                                         all of the argument's data
//
//  bool replace(string name, int ioption, <type> d)       -> same as replace but use a copy of the value d
//  bool replace(string name, int ioption, const fwrap& f) -> same as replace but use a copy of f's data and attributes
//
//  remove(string name) -> delete the variable called name.  This must be done before a new
//                         variable can be added with the same name
//
// The only way to change the storage location of a variable is to remove it first and then re-add it.
//
// --- attributes ---
// There can be global attributes and variable attributes.  These are stored using the fwrap_map_att class 
// inside the fwrap class for variable attributes and inside the container for global attributes.  If the fwrap changes 
// for a variable then the attributes will also change.  However, when data is copied into an fwrap (using fwrapCopy for
// the get,set functions described above), only the numerical data is copied.  So the attribute functions listed 
// below should be used to get access to the attributes. Global attributes are referenced with an empty variable 
// name.  Currently string, double and integer attributes are supported.
//   att_size(string name)  -> return the number of attributes for a variable or the number of 
//                             global attributes when the name = ''
//   att(string name)       -> returns a reference to the fwrap_map_att class for manipulating the attributes.
//                             this reference should be considered temporary since it depends on the existence
//                             of the underlying fwrap.
//
// --- netcdf ---
// The netcdf functions provide a method for saving and loading the container data in a netcdf file.
//
// netcdf can associate one dimensional variables with the axis of other variables by defining the dimension
// name to have the same name as the axis variable name.  This association can be made in the VarContainer
// by defining the variable attributes , ":0", ":1", ... to be the name of the container variables to be used for
// each of its axis or by storing the axis names in the ":" string array attribute.  For a container variable to be able 
// to act as an axis, the rank must be one and the length must match the dimension length.  
// For a container variable without a ":0",":1",... or ":" attribute, a dimension in the
// netcdf file with the name "dim_<n>" will be created where <n> is the length of the axis.
//
// The variable attribute "ctype" will be created to describe the type of the container variable which was
// the source of the data ("int","float","double","logical").  This is used as a clue for the type of variable to
// be created when the data is loaded from a netcdf file.
//
// The fwrap structure orders the dimensions as in fortran while the netcdf file orders the dimensions as in C.
//
//   ncsave(filename,option) -> save the container to a netcdf file in one shot.  The file will be created from scratch
//                              each time this method is called.  Option 0 and 2 will overwrite an existing file.  Option
//                              2 writes to an intermediate file then copies to the destination file.  Option 1 will
//                              throw an exception if the file already exists.
//
//   ncload(filename,option) -> load data from the netcdf file into the container
//                              option=0  clears the container before reading the data
//                              option=1  replaces global attributes and variables from cdf file
//                              option=2  only replace the data of existing variables with matching names in the cdf file,
//                                        the dimensions and attributes are not changed
//                              option=3  similar to 1 except if key already exists, copy the data only
//                              option=4  clear container and load global attributes only
//                              option=5  similar to 3, do not replace global attributes and copy the data if the key already exists
//
// ------------ update -------------
// Same (or similar) options as ncload but for loading a container from another container.
//
//   update(container,option) -> load data from another VarContainer
//                              option=0  clears the container before reading the data
//                              option=1  replaces global attributes and variables from cdf file
//                              option=2  only replace the data of existing variables with matching names in the cdf file,
//                                        the dimensions and attributes are not changed
//                              option=3  similar to 1 except if key already exists, copy the data only
//                              option=4  clear container and load global attributes only
//                              option=5  similar to 3, do not replace global attributes and copy the data if the key already exists
//
class VarContainer {
public:
  //
  // constructor, give this container a name.
  //
  VarContainer(string name) ;

  //
  // constructor with the initial size for each type, give this container a name.
  //
  VarContainer(string name, size_t n_start_size) ;

  ~VarContainer() ;

  //
  // the variable types
  //
  enum varType {NONE=0, I32=1, F32=2, F64=3, L32=4, CHAR=5, BYTES=6} ;

  //
  // return the name
  //
  string name() const { return cname ; } ;

  //
  // return true if this variable is in the container
  //
  bool has_key(string name) const ;
  
  //
  // delete the variable -- no error if not found
  //
  void remove(string name) ;

  //
  // remove all variables from the container and all global attributes
  //
  void clear() ;

  //
  // return a list of the names of a given type or all the names
  //
  list<string> names(varType t=NONE) const ;

  //
  // return the number of variables of a given type
  //
  int size(varType t=NONE) const ;

  //
  // return the type, rank and shape of the variable or throw an exception if
  // the variable does not exist
  //
  varType      getType(string name) const ;
  size_t       getRank(string name) const ;
  const Shape& getShape(string name) const ;

  //
  // return all the info in one shot in the arguments -- type,rank,shape are 
  // only valid if has_key() returns true.  The shape argument must have at
  // least as many elements as the rank of the variable.
  //
  void getInfo(string name, bool& has_key, varType& type, size_t& rank, Shape& shape) const ;

  //
  // add a new variable filled with all zeros, false, spaces or empty (for Bytes)
  //
  void add0(string name, varType type, const Shape& shape) ;

  //
  // integer variable access
  //
  ip_32 getScalarInt32(string name) const ;
  ip_32 getDefScalarInt32(string name, ip_32 def=0) const ;
  void  set(string name, ip_32 d) ;
  void  add(string name, ip_32 d) ;

  aipPtr_32 getInt32(string name, aipPtr_32 def=aipPtr_32()) const ;
  void      get(string name, aip_32&   out, aipPtr_32 def=aipPtr_32()) const ;
  void      get(string name, aipPtr_32 out, aipPtr_32 def=aipPtr_32()) const ;
  void      set(string name, const aip_32& d) ;
  void      set(string name, const aipPtr_32 d) ;
  void      add(string name, aipPtr_32 d) ;
  bool      replace(string name, int ioption, aipPtr_32 d) ;

  inline bool replace(string name, int ioption, ip_32 d) {
    aipPtr_32 p = aipPtr_32::build() ;
    (*p)() = d ;
    return replace(name, ioption, p) ;
  } ;

  inline bool replace(string name, int ioption, const aip_32& d) {
    aipPtr_32 p = aipPtr_32::build(d) ;
    return replace(name, ioption, p) ;
  }

  //
  // floating point variable access
  //
  fp_32 getScalarFloat32(string name) const ;
  fp_32 getDefScalarFloat32(string name, fp_32 def=0) const ;
  void  set(string name, fp_32 d) ;
  void  add(string name, fp_32 d) ;

  afpPtr_32 getFloat32(string name, afpPtr_32 def=afpPtr_32()) const ;
  void      get(string name, afp_32&   out, afpPtr_32 def=afpPtr_32()) const ;
  void      get(string name, afpPtr_32 out, afpPtr_32 def=afpPtr_32()) const ;
  void      set(string name, const afp_32& d) ;
  void      set(string name, const afpPtr_32 d) ;
  void      add(string name, afpPtr_32 d) ;
  bool      replace(string name, int ioption, afpPtr_32 d) ;

  inline bool replace(string name, int ioption, fp_32 d) {
    afpPtr_32 p = afpPtr_32::build() ;
    (*p)() = d ;
    return replace(name, ioption, p) ;
  } ;

  inline bool replace(string name, int ioption, const afp_32& d) {
    afpPtr_32 p = afpPtr_32::build(d) ;
    return replace(name, ioption, p) ;
  }


  fp_64 getScalarFloat64(string name) const ;
  fp_64 getDefScalarFloat64(string name, fp_64 def=0) const ;
  void  set(string name, fp_64 d) ;
  void  add(string name, fp_64 d) ;

  afpPtr_64 getFloat64(string name, afpPtr_64 def=afpPtr_64()) const ;
  void      get(string name, afp_64&   out, afpPtr_64 def=afpPtr_64()) const ;
  void      get(string name, afpPtr_64 out, afpPtr_64 def=afpPtr_64()) const ;
  void      set(string name, const afp_64& d) ;
  void      set(string name, const afpPtr_64 d) ;
  void      add(string name, afpPtr_64 d) ;
  bool      replace(string name, int ioption, afpPtr_64 d) ;

  inline bool replace(string name, int ioption, fp_64 d) {
    afpPtr_64 p = afpPtr_64::build() ;
    (*p)() = d ;
    return replace(name, ioption, p) ;
  } ;

  inline bool replace(string name, int ioption, const afp_64& d) {
    afpPtr_64 p = afpPtr_64::build(d) ;
    return replace(name, ioption, p) ;
  }

  //
  // logical variable access, these methods only work on names defined as Logical.
  // The bool methods set the stored L32 value to lp_32_TRUE or lp_32_FALSE according to
  // the bool value.  The stored L32 is converted to bool with static_cast<bool>().
  //
  lp_32 getScalarLogical32(string name) const ;
  lp_32 getDefScalarLogical32(string name, lp_32 def=0) const ;
  void  setLogical(string name, lp_32 d) ;
  void  addLogical(string name, lp_32 d) ;

  bool  getScalarLogicalBool(string name) const ;
  bool  getDefScalarLogicalBool(string name, bool def=false) const ;
  void  setLogical(string name, bool d) ;
  void  addLogical(string name, bool d) ;

  alpPtr_32 getLogical32(string name, alpPtr_32 def=alpPtr_32()) const ;
  void      getLogical(string name, alp_32&   out, alpPtr_32 def=alpPtr_32()) const ;
  void      getLogical(string name, alpPtr_32 out, alpPtr_32 def=alpPtr_32()) const ;
  void      setLogical(string name, const alp_32& d) ;
  void      setLogical(string name, const alpPtr_32 d) ;
  void      addLogical(string name, alpPtr_32 d) ;
  bool      replaceLogical(string name, int ioption, alpPtr_32 d) ;

  void      getLogical(string name, alp_bool&   out, alpPtr_bool def=alpPtr_bool()) const ;
  void      getLogical(string name, alpPtr_bool out, alpPtr_bool def=alpPtr_bool()) const ;
  void      setLogical(string name, const alp_bool& d) ;
  void      setLogical(string name, const alpPtr_bool d) ;
  bool      replaceLogical(string name, int ioption, alpPtr_bool d) ;  // always copies d data

  // a way for fortran code to set what it thinks true and false is
  inline void  setLogicalTrue(lp_32 truevalue)   { logical_true=truevalue ; } ;
  inline void  setLogicalFalse(lp_32 falsevalue) { logical_false=falsevalue ; } ;

  inline lp_32 getLogicalTrue()  const { return logical_true  ; } ;
  inline lp_32 getLogicalFalse() const { return logical_false ; } ;

  inline bool replaceLogical(string name, int ioption, lp_32 d) {
    alpPtr_32 p = alpPtr_32::build() ;
    (*p)() = d ;
    return replaceLogical(name, ioption, p) ;
  } ;

  inline bool replaceLogical(string name, int ioption, bool d) {
    alpPtr_32 p = alpPtr_32::build() ;
    (*p)() = d ? logical_true : logical_false ;
    return replaceLogical(name, ioption, p) ;
  } ;

  inline bool replaceLogical(string name, int ioption, const alp_32& d) {
    alpPtr_32 p = alpPtr_32::build(d) ;
    return replaceLogical(name, ioption, p) ;
  }
  
  inline bool replaceLogical(string name, int ioption, const alp_bool& d) {
    alpPtr_bool p = alpPtr_bool::build(d) ;
    return replaceLogical(name, ioption, p) ;
  }
  
  // characters, the first dimension represents the fortran character*(len) length.
  // An N dimensional fwrap<string> corresponds to an N+1 dimensional fwrap<char>.  
  // The fwrapCopy_char procedure is used to perform the copies between character variables.
  cp_char getScalarChar(string name) const ;
  cp_char getDefScalarChar(string name, cp_char def=' ') const ;
  void    set(string name, cp_char d) ;
  void    add(string name, cp_char d) ;

  string getScalarString(string name) const ;
  string getDefScalarString(string name, string def="") const ;
  void   set(string name, string d) ;
  void   add(string name, string d) ;

  acpPtr_char getChar(string name, acpPtr_char def=acpPtr_char()) const ;
  void      get(string name, acp_char&   out, acpPtr_char def=acpPtr_char()) const ;
  void      get(string name, acpPtr_char out, acpPtr_char def=acpPtr_char()) const ;
  void      set(string name, const acp_char& d) ;
  void      set(string name, const acpPtr_char d) ;
  void      add(string name, acpPtr_char d) ;
  bool      replace(string name, int ioption, acpPtr_char d) ;

  void      get(string name, acp_string&   out, acpPtr_string def=acpPtr_string()) const ;
  void      get(string name, acpPtr_string out, acpPtr_string def=acpPtr_string()) const ;
  void      set(string name, const acp_string& d) ;
  void      set(string name, const acpPtr_string d) ;
  bool      replace(string name, int ioption, acpPtr_string d) ;

  inline bool replace(string name, int ioption, cp_char d) {
    acpPtr_char p = acpPtr_char::build() ;
    (*p)() = d ;
    return replace(name, ioption, p) ;
  } ;

  inline bool replace(string name, int ioption, const acp_char& d) {
    acpPtr_char p = acpPtr_char::build(d) ;
    return replace(name, ioption, p) ;
  }

  inline bool replace(string name, int ioption, string d) {
    acpPtr_string p = acpPtr_string::build() ;
    (*p)() = d ;
    return replace(name, ioption, p) ;
  } ;

  inline bool replace(string name, int ioption, const acp_string& d) {
    acpPtr_string p = acpPtr_string::build(d) ;
    return replace(name, ioption, p) ;
  }

  // Bytes, each element is a variable size vector<Byte>.  One application of this is to
  // store the contents of a file so file methods are also provided.  This ofcourse can not
  // be based on fortran.  Note the following treat Bytes as values (like the other wrappers)
  //  so involve copies of the data.
  Bytes getScalarBytes(string name) const ;
  Bytes getDefScalarBytes(string name, Bytes def=Bytes()) const ;
  void  set(string name, Bytes d) ;
  void  add(string name, Bytes d) ;
  
  void getScalarBytesToFile(string name,   string filename) ;  // copy the scalar bytes to a file overwriting an existing file
  void setScalarBytesFromFile(string name, string filename) ;  // move the contents of the file into a scalar Bytes
  void addScalarBytesFromFile(string name, string filename) ;  // add the file to the container as a scalar Bytes
  bool replaceScalarBytesFromFile(string name, int ioption, string filename) ; // replace a scalar Bytes with the file data
  
  size_t sizeScalarBytes(string name) ;  // return the number of bytes in the scalar Bytes

  abpPtr_bytes getBytes(string name, abpPtr_bytes def=abpPtr_bytes()) const ;
  void      get(string name, abp_bytes&   out, abpPtr_bytes def=abpPtr_bytes()) const ;
  void      get(string name, abpPtr_bytes out, abpPtr_bytes def=abpPtr_bytes()) const ;
  void      set(string name, const abp_bytes& d) ;
  void      set(string name, const abpPtr_bytes d) ;
  void      add(string name, abpPtr_bytes d) ;
  bool      replace(string name, int ioption, abpPtr_bytes d) ;

  inline bool replace(string name, int ioption, Bytes d) {
    abpPtr_bytes p = abpPtr_bytes::build() ;
    (*p)() = d ;
    return replace(name, ioption, p) ;
  } ;

  inline bool replace(string name, int ioption, const abp_bytes& d) {
    abpPtr_bytes p = abpPtr_bytes::build(d) ;
    return replace(name, ioption, p) ;
  }

    
  //
  // ----------- attributes ------------
  // the reference returned by these functions should be considered temporary
  // because they will become invalid if the underlying fwrap ends up being deleted.
  //
  // name = variable name or the empty string for global attributes
  //
  // after calling att(string), the fwrap_map_att* is stored in an internal variable which
  // can be retrieved with temp_att().  This pointer is also set to the attributes of any newly added variable with 
  // an add() or replace() method but the pointer is cleared by any method which could change the 
  // container (e.g. set() methods).
  //
  fwrap_map_att*       att(string vname) ;          // return the object for storing attributes or 0 if name is bad
  const fwrap_map_att* att(string vname) const ;    // return the object for storing attributes or 0 if name is bad

  fwrap_map_att*       att() ;                      // return the object for storing the global attributes
  const fwrap_map_att* att() const ;                // return the object for storing the global attributes

  fwrap_map_att*       temp_att() { return patt ; } ;                   // after an att(string) call, the temporary attribute
  const fwrap_map_att* temp_att(string vname) const { return patt ; } ; // is set to this attribute pointer

  //
  // --------- share ---------
  // copy pointers from the argument box to this container using the same name.
  // The variables end up sharing the same memory address so changing the value in one
  // container changes it in the other.
  void share(const VarContainer& box, const list<string>& names, int option=0) ;
                               // box    = source of wrapped variable pointers
                               // names  = names of variables in box which will be put in this container
                               // option = 0  -> replace variable if the name is already in use 
                               //        = 1  -> do not replace existing variables
                               //        = 2  -> throw an exception if there is an existing variable                                                     

  //
  // --------- netcdf ---------
  // 
  void ncsave(string filename, int option=0, list<string> yes=list<string>(), list<string> no=list<string>()) const ;  
                                                      // save container to a netcdf file
                                                      //   option = 0 -> overwrite an existing file
                                                      //   option = 1 -> do not overwrite an existing file
                                                      //   option = 2 -> save to <filename>_temp then do a system move 
                                                      //                 to <filename> overwriting the existing file
                                                      //   yes = if not empty, for a variable to be included in the netcdf file, it 
                                                      //         must contain at least one of the attributes in this list as a scalar integer 
                                                      //         attribute with a nonzero value
                                                      //   no  = for a variable to be included in the netcdf file, it must not
                                                      //         contain a variable in this list as a scalar integer attribute with a 
                                                      //         nonzero value

  // ------ ncload & update & receive matrix ------
  //     clear     load global atts    duplicate_vars  duplicate_atts  add_new_vars
  //  0    X             X                  R&A              X             X           R&A->remove existing variable and add
  //  1                  X                  R&A              X             X             C->copy to existing variable
  //  2                                      C                                           X->data/attributes changed to source
  //  3                  X                   C               X             X             
  //  4    X             X 
  //  5                                      C               X             X  
  //  6                  X                   C                             X  
  //
  void ncload(string filename, int option=0) ;        // load container variables from a netcdf file
                                                      // option = 0 -> clear the container then load
                                                      // option = 1 -> for duplicates, replace existing variables
                                                      // option = 2 -> only copy data into existing variables, don't change attributes
                                                      // option = 3 -> copy duplicates, everything else is replaced including global attributes
                                                      // option = 4 -> clear the container then load global attributes only
                                                      // option = 5 -> like 3 but don't load global attributes

  // ------- update -------
  void update(const VarContainer& box, int option=0) ; // load container from another container
                                                      // option = 0 -> clear the container then load
                                                      // option = 1 -> for duplicates, replace variables and change attributes
                                                      // option = 2 -> for duplicates, copy data to existing container
                                                      // option = 3 -> copy duplicates, everything else is replaced
                                                      // option = 4 -> clear the container then load global attributes only
                                                      // option = 5 -> like 3 but don't load global attributes

  // --------- Serialization -------
  void send(VarWriter* writer, list<string> names=list<string>(),
                               list<string> yes=list<string>(), 
                               list<string> no=list<string>() ) const ;
                                                      //   names = if not empty, only write variable/attribute data
                                                      //           for the listed names. Global attributes are always written.
                                                      //           if names contains only "GLOBAL_ATTS", no variable names are written.
                                                      //   yes =   if not empty, for a variable to be included in the write, it 
                                                      //           must contain at least one of the attributes in this list as a scalar integer 
                                                      //           attribute with a nonzero value.
                                                      //   no  =   for a variable to be included in the write, it must not
                                                      //           contain a variable in this list as a scalar integer attribute with a 
                                                      //           nonzero value.
  
  void receive(VarReader* reader, int option, list<string> keep=list<string>(), 
                                              list<string> protect=list<string>())  ;
                                                      // option = 0 -> clear the container then load (ignores keep,protect)
                                                      // option = 1 -> for duplicates, replace existing variables (but see keep,protect)
                                                      // option = 2 -> only copy data into existing variables (but see protect), don't change attributes
                                                      // option = 3 -> copy duplicates, everything else is replaced (but see keep,protect) including global attributes
                                                      // option = 4 -> clear the container (ignores keep,protect) then load global attributes only
                                                      // option = 5 -> like 3 but don't load global attributes
                                                      // option = 6 -> like 3 but don't change existing variable attributes
                                                      // keep   = if not empty, a variable which contains a scalar integer attribute
                                                      //          with one of these names with nonzero value will not be replaced but 
                                                      //          will still allow the data to be set, attributes will not be changed
                                                      // protect= if not empty, a variable which contains a scalar integer attribute
                                                      //          with one of these names with nonzero value will not allow the
                                                      //          data or attributes to be changed.  Appropriate for variables considered "constant"
  class ContainerError : public exception {
  public:
    ContainerError(string container, string info) throw() ;
    ContainerError(const ContainerError&) throw() ;
    ContainerError& operator=(const ContainerError&) throw() ;
    virtual ~ContainerError() throw() ;
  
    virtual const char*   what()  const throw() ;
    virtual const string  whats() const throw() ;
  protected:
    std::string description ;   // string for use in error message
  } ;

  class KeyNotFound : public ContainerError {
  public:
    // vtype=-1 for logical, vtype=-2 for character, vtype=-3 for bytes, vtype=1 for numerical
    KeyNotFound(string container, string varname, int vtype=0) throw() ;
    KeyNotFound(string container, string varname, varType need, varType got) throw() ;
    KeyNotFound(const KeyNotFound&) throw() ;
    KeyNotFound& operator=(const KeyNotFound&) throw() ;
    virtual ~KeyNotFound() throw() ;
  } ;

  // basically a ContainerError which distinguishes the ioness of the error.
  class IOError : public ContainerError {
  public:
    IOError(string container, string info) throw() ;
    IOError(const IOError&) throw() ;
    IOError& operator=(const IOError&) throw() ;
    virtual ~IOError() throw() ;
  } ;

  // ----------- i/o -----------
  friend ostream& operator<<(ostream& s, const VarContainer& v) ;

  static string typeWrite(varType t) ;   // return a type specific string
private:
  // --- static ---
  static const char*  dnames[] ;                // names of attribute dimensions ":0",":1",...,":6"
  static const size_t N_START_SIZE ;            // beginning size of each of the gatekeepers for default constructor
  static const size_t N_START_MULT ;            // beginning size of the keymaster is N_START_MULT*N_START_SIZE
  static const int    STORE_AXIS_COLON_ARRAY ;  // nonzero to store axis in a ":" array attribute, 0 to use original ":n" attributes

  // --- helper ---
  string dimname(int i) const ;                      // return the default name for an unnamed dimension
  void ncerror(int serror, string extra="") const ;  // throw an exception if the argument indicates a netcdf error
  void ncaddatt(int ncidp, int varid, const fwrap_map_att* a, string ctype) const ;  // add attributes to a variable
                                                                                     // ctype is an extra string attribute if 
                                                                                     // nonempty

  static bool ncinc(const fwrap_map_att* a, const list<string>& yes, const list<string>& no) ;  // return true if the variable
                                                                   // should be included in the ncsave() netcdf file based on the attributes

  static int nckeep(const fwrap_map_att* a, const list<string>& keep, const list<string>& protect) ; 
                            // return 0 if no attribute of a is in keep or protect
                            //        1 if an attribute of a is in keep but not protect
                            //        2 if an attribute of a is in protect
  fwrap_map_att ncgetatt(int ncidp, int varid, int ncatt) const ;  // get attributes from netcdf file, handle 
     
  // -- IntStore --
  // int storage space which automatically grows then cleans up at destruction  
  class IntStore {
  public:
    IntStore(size_t initsize) : dsize(initsize), data(0) { data = new int[initsize] ; } ;
    ~IntStore() { delete[] data ; } ;
    
    void need(size_t nsize) ;        // insure there is this much space, will copy existing data when resizing
    int* ptr() const { return data ; };
  private:
    IntStore(const IntStore&) ;
    IntStore& operator=(const IntStore&) ;  // not allowed
    
    size_t dsize ;   // number allocated
    int*   data ;    // int[dsize]
  } ;

  // -- ByteStore --
  // byte storage space which automatically grows then cleans up at destruction  
  class ByteStore {
  public:
    ByteStore(size_t initsize) : dsize(initsize), data(0) { data = new Byte[initsize] ; } ;
    ~ByteStore() { delete[] data ; } ;
    
    void need(size_t nsize) ;        // insure there is this much space, will copy existing data when resizing
    Byte* ptr() const { return data ; };
  private:
    ByteStore(const ByteStore&) ;
    ByteStore& operator=(const ByteStore&) ;  // not allowed
    
    size_t dsize ;   // number allocated
    Byte*  data ;    // Byte[dsize]
  } ;

  // -- DoubleStore --
  // double storage space which automatically grows then cleans up at destruction
  class DoubleStore {
  public:
    DoubleStore(size_t initsize) : dsize(initsize), data(0) { data = new double[initsize] ; } ;
    ~DoubleStore() { delete[] data ; } ;
    
    void need(size_t nsize) ;        // insure there is this much space, will copy existing data when resizing
    double* ptr() const { return data ; } ;
  private:
    DoubleStore(const DoubleStore&) ;
    DoubleStore& operator=(const DoubleStore&) ;  // not allowed
    
    size_t   dsize ;   // number allocated
    double*  data ;    // double[dsize]
  } ;
  
  // -- StringStore --
  // string storage space which automatically grows then cleans up at destruction
  // This class can accept an array of strings and convert to a single char* with
  // the strings deliminated by zeros.
  class StringStore {
  public:
    StringStore(size_t initsize) ;        // initialize with the default number of strings
    ~StringStore() ;
    
    void need(size_t nsize) ;        // insure that this many strings can be stored -- does NOT copy
    void cneed(size_t nsize) ;       // insure a char array of this size (including terminating 0) can be stored
                                     // does NOT copy
    
    size_t stringtochar(size_t isize) ;  // convert the isize number of strings to the char array and return the 
                                         // char size including 0, zeros in the strings are set to spaces
    size_t chartostring(size_t csize) ;  // convert the char array of size csize including 0 or not to an array of 
                                         // strings and return the number of strings.  A 0 argument results in a 0 output.
    
    string* ptr()   const { return sdata ; } ; // get the pointer to the array of strings
    char*   cptr()  const { return cdata ; } ; // get the pointer to the array of chars
  private:
    StringStore(const StringStore&) ;
    StringStore& operator=(const StringStore&) ;  // not allowed
    
    static const char CHAR_TERM ;   // character used to terminate a string in an array (1)
    static const char CHAR_LAST ;   // character used to terminate final string in an array (0)
    static const char CHAR_EMPTY ;  // character used to indicate the string is really empty (2)
    static const char CHAR_SUB ;    // character used in place of CHAR_TERM,CHAR_EMPTY in strings (space)

    size_t  ssize ;    // number of strings allocated
    string* sdata ;    // string array
    size_t  csize ;    // number of chars allocated
    char*   cdata ;    // char array
  } ;

  // isolate type dependent netcdf variable access
  static int nc_get_var_ip32(int ncidp, int varid, ip_32* p) ;   // get for ip_32 type
  static int nc_get_var_fp32(int ncidp, int varid, fp_32* p) ;   // get for fp_32 type
  static int nc_get_var_fp64(int ncidp, int varid, fp_64* p) ;   // get for fp_64 type
  static int nc_get_var_lp32(int ncidp, int varid, lp_32* p) ;   // get for lp_32 type
  static int nc_get_var_char(int ncidp, int varid, cp_char* p) ; // get for lp_32 type

  static int nc_put_var_ip32(int ncidp, int varid, const ip_32* p) ;   // put for ip_32 type
  static int nc_put_var_fp32(int ncidp, int varid, const fp_32* p) ;   // put for fp_32 type
  static int nc_put_var_fp64(int ncidp, int varid, const fp_64* p) ;   // put for fp_64 type
  static int nc_put_var_lp32(int ncidp, int varid, const lp_32* p) ;   // put for lp_32 type
  static int nc_put_var_char(int ncidp, int varid, const cp_char* p) ; // put for cp_char type
  
  static int nc_put_var_bytes(int ncidp, int svarid, int dvarid, size_t np, const Bytes* p) ;      // put for Bytes type
  static int nc_get_var_bytes(int ncidp, int svarid, int dvarid, size_t np, size_t nb, Bytes* p) ; // get for Bytes type

  // illegal methods
  VarContainer(const VarContainer&) ;
  VarContainer& operator=(const VarContainer&) ;

  // types
  typedef pair<varType,size_t> iPair ;                     // container type, index into GateKeeper

#ifdef VARCONTAINER_UNOMAP
  typedef unordered_map<string,iPair> KeyMaster ;          // map from string to location of fwrapPtr in the GateKeepers
#else
#ifdef VARCONTAINER_HASHMAP
  typedef hash_map<string,iPair,hash_string> KeyMaster ; 
#else
  typedef map<string,iPair> KeyMaster ;
#endif
#endif
  
  void replaceInfo(string name, varType type, size_t k, size_t& rank, Shape& shape) const ; // info after hash search
  void replaceRemove(KeyMaster::iterator& i, string name, varType t, size_t k) ;            // replace after hash search

  // --- data ---
  string    cname ;      // name of the container

  gkInt     igate ;      // hold integer*4
  gkFloat   fgate ;      // hold real*4
  gkDouble  dgate ;      // host real*8
  gkLogical lgate ;      // hold logical*4
  gkChar    cgate ;      // holds character*(*)
  gkBytes   bgate ;      // holds valarray<Byte>

  KeyMaster zool ;       // contains all of the keys for finding the fwrapPtr

  lp_32  logical_true ;  // hold fortran TRUE
  lp_32  logical_false ; // hold fortran FALSE

  fwrap_map_att matt ;   // global attributes

  mutable fwrap_map_att* patt ;  // pointer for C interface access to attributes, this is cleared by any method
                                 // which could change the container
} ;

ostream& operator<<(ostream& s, const VarContainer& v) ;


#endif
