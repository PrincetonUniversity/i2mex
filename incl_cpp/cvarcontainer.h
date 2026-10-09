
#include "varcontainer.h"

//
// C callable interface to VarContainer
//

//
// ======================== cVarContainer ======================
//
// This is the class a C interface will work with.  It contains storage for the list<string>
// results which C can not handle.
//
class cVarContainer {
 public:
  cVarContainer(string name) ;

  cVarContainer(string name, int nsize) ; // construct with a default size

  ~cVarContainer() ;

  VarContainer*       object()       { return container ; } ;
  const VarContainer* object() const { return container ; } ;

  // --- data ---
  int                           names_size ;   // number of stored names
  list<string>                  names ;        // temporary for holding variable names
  list<string>::const_iterator  names_citer ;  // temporary for iterating over variable names

  int                           att_names_size ;   // number of stored attribute names
  list<string>                  att_names ;        // temporary for holding attribute names
  list<string>::const_iterator  att_names_citer ;  // temporary for iterating over attribute names

  aipPtr_32 local_i32 ;        // for cvc_local() method
  afpPtr_32 local_f32 ;
  afpPtr_64 local_f64 ;
  alpPtr_32 local_l32 ;
  
  acpPtr_char  local_char ;
  abpPtr_bytes local_bytes ;

 private:
  VarContainer* container ;    // only instance, read only

  cVarContainer(const cVarContainer& a ) ;               // no copying allowed
  cVarContainer& operator=(const cVarContainer& a ) ;
} ;

//
// a singleton class, lazily constructs the VarContainer on first access with cobject()
//
class cVarContainerSingle {
 public:
  cVarContainerSingle(string name) ;
  cVarContainerSingle(string name, int nsize) ;

  ~cVarContainerSingle() ;
  
  cVarContainer* cobject() ;
 private:
  string         sname ;
  int            ssize ;
  cVarContainer* single ;
} ;

//
// ------------------ object stuff --------------------
//
// build a new instance
//   name = name of cVarContainer to build
//   self = returned new instance
//
extern "C" void  F77NAME(cvc_new)(const char* name, cVarContainer** self, int* ier)  ;

//
// build a new instance specifying the hash set size
//   name = name of cVarContainer to build
//   nsize = number of quantities expected in the container
//   self = returned new instance
//
extern "C" void  F77NAME(cvc_new_size)(const char* name, int* nsize, cVarContainer** self, int* ier)  ;

//
// delete the instance
//   self = instance
//
extern "C" void F77NAME(cvc_delete)(cVarContainer** self, int* ier) ;

//
// copy the instance pointer -- note data is not copied.  This exists to aid
// fortran which knows nothing about the size of a cVarContainer*
//   self = instance
//   you  = empty
//
extern "C" void F77NAME(cvc_copy)(cVarContainer** self, cVarContainer** you) ;

//
// hold onto a pointer for me to avoid fortran f90 copying
//   value = fortran pointer
//
extern "C" void F77NAME(cvc_hold)(void* value) ;

//
// return nonzero if there is a held value
//   self  = instance
//
extern "C" int F77NAME(cvc_isheld)() ;

//
// ----------- management ----------
//
// print out to stdout
//   self = instance
//   ier  = nonzero on error
//
extern "C" void F77NAME(cvc_show)(cVarContainer** self, int* ier) ;

//
// Return the number of variables by type
//   self  = instance
//   type  = type code (0 for all types)
//   isize = returned size
//   ier   = nonzero if there was an error
extern "C" void F77NAME(cvc_size)(cVarContainer** self, int* type, int* isize, int* ier) ;

//
// Return the number of variables by type and store their names internally.  Retrieve
// the names with cvc_names_fetch()
//   self  = instance
//   type  = type code (0 for all types)
//   isize = number of names
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_names)(cVarContainer** self, int* type, int* isize, int* ier) ;

//
// After a cvc_names() call, this will retrieve each name in succession like an iterator.
//   self  = instance
//   idone = nonzero when finished
//   clen  = size of cname variable
//   cname = return the name here
//   ier   = nonzero on error   
extern "C" void F77NAME(cvc_names_fetch)(cVarContainer** self, int* idone, int* clen, char* cname, int* ier) ;

//
// return nonzero if this name has a key
//   self   = instance
//   name   = look for this name
//   result = nonzero when container has this key
//   ier    = nonzero on error
extern "C" void F77NAME(cvc_has_key)(cVarContainer** self, const char* name, int* result, int* ier) ;

//
// remove this name from the container, no error if the name is not in the container
//   self   = instance
//   name   = name to be removed
//   ier    = nonzero on error
extern "C" void F77NAME(cvc_remove)(cVarContainer** self, const char* name, int* ier) ;

//
// remove everything from the container
//   self  = instance
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_clear)(cVarContainer** self, int* ier) ;

//
// return the type of a variable,  
//    0 for none
//    1 for Integer*4
//    2 for Real*4
//    3 for Real*8
//    4 for Logical*4
//  self = instance
//  name = name of variable
//  type = returned type
//  ier  = nonzero on error
extern "C" void F77NAME(cvc_gettype)(cVarContainer** self, const char* name, int* type, int* ier) ;

//
// return the rank of a variable
//   self = instance
//   name = name of variable
//   rank = returned rank
//   ier  = nonzero on error
extern "C" void F77NAME(cvc_getrank)(cVarContainer** self, const char* name, int* rank, int* ier) ;

//
// return the shape of a variable
//   self  = instance
//   name  = name of variable
//   isize = size of array being passed in
//   shape = filled with the shape of the variable
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_getshape)(cVarContainer** self, const char* name, const int* isize, int* shape, int* ier) ;

// return general info about this variable
//   self  = instance
//   name  = name of variable
//   has_key = nonzero when variable is in the container
//   type    = type of the variable
//   rank    = rank of variable
//   isize   = size of shape variable
//   shape   = filled with shape of variable
//   ier     = nonzero on error
extern "C" void F77NAME(cvc_getinfo)(cVarContainer** self, const char* name, int* has_key, int* type, int* rank,
				     const int* isize, int* shape, int* ier) ;

//
// Fills the local_* instance variables with the variable data -- for debugging.
//   self  = instance
//   name  = name of the variable
//   iout  = nonzero to write out the variable to stdout
//   ier   = nonzero on error   
extern "C" void F77NAME(cvc_local)(cVarContainer** self, const char* name, int* iout, int* ier);

//
// ---------------- attributes ---------------
//
// integer,real*8 and character*(*) attributes can be associated with each variable
// and globally with the container.  Some attributes are special
//  "ctype"     used for netcdf operations, do not use this as a general attribute
// for a variable attribute
//  ":0" -> name of the variable for the first axis
//  ":1" -> name of the variable for the second axis
//  ... and so on ...
// The ":<n>" attributes provide an association of an axis of a variable with another
// rank 1 variable with length equal to the length of the axis.  These attributes are translated
// to dimension entries in the netcdf file.
//

//
// Set up for editing the attributes of a variable or the global attributes.
// Call this before any other attribute access subroutine.  The attributes of a variable
// are also setup for editing right after the variable is added to the container.  The internal
// attribute pointer will be cleared when the container is modified (e.g. a set method).
//   self = instance
//   name = name of variable
//   ier  = nonzero on error 
// 
extern "C" void F77NAME(cvc_att)(cVarContainer** self, const char* name, int* ier) ;
  
// internal
fwrap_map_att* cvc_test_att_temp(cVarContainer** self, int* ier) ;

//
// Copy the attributes from the source variable to the destination variable overwriting all of
// the destination attributes.
//   self = instance
//   name_src  = name of variable which is the source of the attributes
//   name_dest = name of variable for destination of the attributes
//   ier       = nonzero if source or destination variables do not exist
//
extern "C" void F77NAME(cvc_att_copy)(cVarContainer** self, const char* name_src, const char* name_dest, int* ier) ;

//
// get number of attributes
//   self = instance
//   inum = returned number of instance
//   ier  = nonzero on error
extern "C" void F77NAME(cvc_att_size)(cVarContainer** self, int* inum, int* ier) ;

//
// remove all the attributes from the current variable
//   self = instance
//   ier  = nonzero on error
extern "C" void F77NAME(cvc_att_clear)(cVarContainer** self, int* ier) ;

//
// return the type of an attribute
//   self  = instance
//   name  = name of the attribute
//   itype = type of attribute
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_type)(cVarContainer** self, const char* name, int* itype, int* ier) ;

//
// return the number of elements in an attribute
//   self  = instance
//   name  = name of the attribute
//   isize = number of values for this attribute
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_vsize)(cVarContainer** self, const char* name, int* isize, int* ier) ;

//
// return info about an attribute
//   self      = instance
//   name      = name of the attribute
//   itype     = type of attribute
//   isize     = number of values for this attribute
//   maxLength = for a string attribute this is the size of the longest string
//   ptr       = return the fwrap_att* -- this pointer is temporary and should be used immediately then discarded
//   ier       = nonzero on error
extern "C" void F77NAME(cvc_att_info)(cVarContainer** self, const char* name, int* itype, int* isize, int* maxLength,
                                      fwrap_att** ptr, int* ier) ;

//
// return the value of an integer scalar attribute
//   self  = instance
//   name  = name of integer attribute
//   ival  = integer value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_int)(cVarContainer** self, const char* name, int* ival, int* ier)  ;

//
// return the value of an integer vector attribute -- cvc_att_info should be called before this to get the pointer and size
//   ptr   = pointer to attribute
//   ilen  = number of elements in ival argument
//   ival  = integer value array
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_vint)(fwrap_att** ptr, int* ilen, int* ival, int* ier)  ;


//
// return the value of a double attribute
//   self  = instance
//   name  = name of double attribute
//   dval  = double value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_double)(cVarContainer** self, const char* name, double* dval, int* ier)  ;

//
// return the value of a double vector attribute -- cvc_att_info should be called before this to get the pointer and size
//   ptr   = pointer to attribute
//   dlen  = number of elements in dval argument
//   dval  = double value array
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_vdouble)(fwrap_att** ptr, int* dlen, double* dval, int* ier)  ;


//
// return the value of a double attribute
//   self  = instance
//   name  = name of double attribute
//   dval  = double value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_fp64)(cVarContainer** self, const char* name, fp_64* dval, int* ier) ;

//
// return the value of a double vector attribute -- cvc_att_info should be called before this to get the pointer and size
//   ptr   = pointer to attribute
//   dlen  = number of elements in dval argument
//   dval  = double value array
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_vfp64)(fwrap_att** ptr, int* dlen, fp_64* dval, int* ier)  ;


//
// return the length of a string attribute
//   self  = instance
//   name  = name of attribute
//   ilen  = returned length
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_slength)(cVarContainer** self, const char* name, int* ilen, int* ier);

//
// return a string attribute
//   self  = instance
//   name  = name of attribute
//   clen  = length of passed char array
//   cname = returned string
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_string)(cVarContainer** self, const char* name, int* clen, char* cname, int* ier);

//
// return one value of a string vector attribute -- cvc_att_info should be called before this to get the pointer and size
//   ptr   = pointer to attribute
//   ix    = index of string in the attribute
//   clen  = length of char argument cval
//   cval  = string will be stored here, a trailing 0 will be added if there is room
//   slen  = actual string length stored in cval
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_get_vxstring)(fwrap_att** ptr, int* ix, int* clen, char* cval, int* slen, int* ier)  ;


//
// return the value of an integer attribute or the default
//   self  = instance
//   name  = name of integer attribute
//   ival  = integer value
//   ier   = nonzero on error
//   def   = integer default
extern "C" void F77NAME(cvc_att_get_intdef)(cVarContainer** self, const char* name, int* ival, int* ier, int* def)  ;

//
// return the value of an integer vector attribute or the default -- cvc_att_info should be called before this to get the pointer and size
//   ptr   = pointer to attribute
//   ilen  = number of elements in ival argument
//   ival  = integer value array
//   ier   = nonzero on error
//   dlen  = length of default array
//   dval  = default integer array
//   iout  = number of items written to ival
extern "C" void F77NAME(cvc_att_get_vintdef)(fwrap_att** ptr, int* ilen, int* ival, int* ier, int* dlen, int* dval, int* iout)  ;


//
// return the value of a double attribute or the default
//   self  = instance
//   name  = name of double attribute
//   dval  = double value
//   ier   = nonzero on error
//   def   = double default
extern "C" void F77NAME(cvc_att_get_doubledef)(cVarContainer** self, const char* name, double* dval, int* ier, double* def)  ;

//
// return the value of an double vector attribute or the default -- cvc_att_info should be called before this to get the pointer and size
//   ptr   = pointer to attribute
//   ilen  = number of elements in ival argument
//   ival  = double value array
//   ier   = nonzero on error
//   dlen  = length of default array
//   dval  = default double array
//   iout  = number of items written to ival
extern "C" void F77NAME(cvc_att_get_vdoubledef)(fwrap_att** ptr, int* ilen, double* ival, int* ier, int* dlen, double* dval, int* iout)  ;

//
// return the value of a double attribute or the default
//   self  = instance
//   name  = name of double attribute
//   dval  = double value
//   ier   = nonzero on error
//   def   = default
extern "C" void F77NAME(cvc_att_get_fp64def)(cVarContainer** self, const char* name, fp_64* dval, int* ier, fp_64* def) ;

//
// return the value of an double vector attribute or the default -- cvc_att_info should be called before this to get the pointer and size
//   ptr   = pointer to attribute
//   ilen  = number of elements in ival argument
//   ival  = double value array
//   ier   = nonzero on error
//   dlen  = length of default array
//   dval  = default double array
//   iout  = number of items written to ival
extern "C" void F77NAME(cvc_att_get_vfp64def)(fwrap_att** ptr, int* ilen, fp_64* ival, int* ier, int* dlen, fp_64* dval, int* iout)  ;


//
// return the length of a string attribute  or the default
//   self  = instance
//   name  = name of attribute
//   ilen  = returned length
//   ier   = nonzero on error
//   def   = default
extern "C" void F77NAME(cvc_att_get_slengthdef)(cVarContainer** self, const char* name, int* ilen, int* ier, const char* def);

//
// return a string attribute or the default
//   self  = instance
//   name  = name of attribute
//   clen  = length of passed char array
//   cname = returned string
//   ier   = nonzero on error
//   def   = default
extern "C" void F77NAME(cvc_att_get_stringdef)(cVarContainer** self, const char* name, int* clen, char* cname, int* ier, 
					    const char* def);

// no vector string default since this seems to be rather language specific

//
// remove an attribute
//   self  = instance
//   name  = name of attribute
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_remove)(cVarContainer** self, const char* name, int* ier);


//
// set the value of an integer scalar attribute
//   self  = instance
//   name  = name of integer attribute
//   ival  = integer value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_set_int)(cVarContainer** self, const char* name, int* ival, int* ier) ;

//
// set the value of an integer vector attribute
//   self  = instance
//   name  = name of integer attribute
//   ilen  = number of elements in ival
//   ival  = integer value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_set_vint)(cVarContainer** self, const char* name, int* ilen, int* ival, int* ier) ;


//
// set the value of a double scalar attribute
//   self  = instance
//   name  = name of integer attribute
//   dval  = double value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_set_double)(cVarContainer** self, const char* name, double* dval, int* ier)  ;

//
// set the value of a double vector attribute
//   self  = instance
//   name  = name of double attribute
//   ilen  = number of elements in ival
//   ival  = double value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_set_vdouble)(cVarContainer** self, const char* name, int* ilen, double* ival, int* ier) ;


//
// set the value of a double attribute
//   self  = instance
//   name  = name of integer attribute
//   dval  = double value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_set_fp64)(cVarContainer** self, const char* name, fp_64* dval, int* ier)  ;

//
// set the value of a fp_64 vector attribute
//   self  = instance
//   name  = name of double attribute
//   ilen  = number of elements in ival
//   ival  = double value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_set_vfp64)(cVarContainer** self, const char* name, int* ilen, fp_64* ival, int* ier) ;


//
// set the value of a string attribute
//   self  = instance
//   name  = name of integer attribute
//   cval  = string value
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_set_string)(cVarContainer** self, const char* name, const char* cval, int* ier);

//
// create a string attribute with a fixed number of strings and return the pointer
//   self      = instance
//   name      = name of the attribute
//   isize     = number of strings in attribute
//   ptr       = return the fwrap_att* -- this pointer is temporary and should be used immediately then discarded
//   ier       = nonzero on error
extern "C" void F77NAME(cvc_att_create_vstring)(cVarContainer** self, const char* name, int* isize, fwrap_att** ptr, int* ier) ;

//
// set one value of a string vector attribute -- cvc_att_create_vstring should be called before this to get the pointer to the new attribute
//   ptr   = pointer to attribute
//   ix    = index of string in the attribute
//   clen  = length of char argument cval
//   cval  = string to be stored without the terminating 0
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_set_vxstring)(fwrap_att** ptr, int* ix, int* clen, const char* cval, int* ier)  ;


//
// setup for reading the names of the attributes for the current variable with cvc_att_name_fetch()
//   self  = instance
//   isize = returned number of attribute names
//   ier   = nonzero on error
extern "C" void F77NAME(cvc_att_names)(cVarContainer** self, int* isize, int* ier) ;

//
// After a cvc_att_names() call, this will retrieve each name in succession like an iterator.
//   self  = instance
//   idone = nonzero when finished
//   clen  = size of cname variable
//   cname = return the name here
//   ier   = nonzero on error   
extern "C" void F77NAME(cvc_att_name_fetch)(cVarContainer** self, int* idone, int* clen, char* cname, int* ier);

//
// ---------------- netcdf -----------------
//
// Save all the container variables and attributes in a netcdf file.  The variable ":<n>" attributes
// are translated to dimension entries.  Variables without a named axis will have dimensions called dim_<m>
// where <m> is the length of the axis.  The variables stored can be restricted by specifying scalar nonzero integer
// attributes which the variable must have (yes attributes) and must not have (no attributes).
//
//     self     = instance
//     filename = name of netcdf file to create
//     option   = 0 -> overwrite existing file
//                1 -> do not overwrite an existing file
//                2 -> save to a temporary file then copy to filename
//     ier      = nonzero on error
//
//     nyes  = number of yes attributes
//     kyes  = array of indices into syes where the ith attribute is stored
//     syes  = null terminated attribute strings indexed by syes+kyes[i] for i=0, nyes
//
//     nno   = number of no attributes
//     kno   = array of indices into sno where the ith attribute is stored
//     sno   = null terminated attribute strings indexed by sno+kno[i] for i=0, nno
 
extern "C" void F77NAME(cvc_ncsave)(cVarContainer** self, const char* filename, int* option, int* ier,
				    const int* nyes, const int* kyes, const char* syes,
				    const int* nno,  const int* kno,  const char* sno) ;

// Load the container from the netcdf file.  Multidimensional attributes will be ignored.  Dimensions
// pointing to valid variables will be translated to ":<n>" attributes.
//
//     self     = instance
//     filename = name of an existing netcdf file
//     option   = 0 -> clear the container before loading the file
//                1 -> clear global attributes and replace variables and attributes contained in the file
//                     but do not modify variables not referenced in the file
//                2 -> do not change dimension or attributes, only set the container data based on the
//                     data in the netcdf file
//                3 -> similar to 1 except existing variables have their data copied (like 2) instead of replaced
//                4 -> clear the container and read in the global attributes only
extern "C" void F77NAME(cvc_ncload)(cVarContainer** self, const char* filename, int* option, int* ier) ;


// ======================= Variables ========================
// ------------------- scalar only ------------------
// Scalar container operations act like values -- the address of the
// variable is not stored in the container.
//
//   <type> = int32=1,float32=2,float64=3,logical32=4,char=5,bytes=6
//
// 
// getscalar<type>(self,name,value,ier)
//   self   = instance
//   name   = name of scalar
//   value  = returned value
//   ier    = return nonzero on error
//
// getdefscalar<type>(self,name,value,def,ier)
//   self   = instance
//   name   = name of scalar
//   value  = returned value
//   def    = same type as <type>, will be used without error
//            if a variable called name does not exist
//   ier    = return nonzero on error
//
// setscalar<type>(self,name,value,ier)
//   self   = instance
//   name   = name of scalar
//   value  = of type <type>, set contained variable to this value
//   ier    = return nonzero on error
//
// addscalar<type>(self,name,value,ier)
//   self   = instance
//   name   = name of scalar
//   value  = of type <type>, add a new scalar contained variable and set it to this value,
//            the argument address is not used for the storage
//   ier    = return nonzero on error
//
//
//
// ------------------- scalar or arrays -------------------
//
// Retrieve data from the container by copying.
//
// get<type>(self,name,rank,shape,array,ier)
//   self   = instance
//   name   = name of variable 
//   rank   = rank of array argument
//   shape  = shape of array argument
//   array  = storage to recieve the data, must have the same rank and extents
//            equal to or larger then the contained variable
//   ier    = nonzero on error
//
// getdef<type>(self,name,rank,shape,array,drank,dshape,darray,ier)
//   self    = instance
//   name    = name of variable 
//   rank    = rank of array argument
//   shape   = shape of array argument
//   array   = storage to recieve the data, must have the same rank and extents
//            equal to or larger then the contained variable
//   drank   = rank of array argument
//   dshape  = shape of array argument
//   darray  = optional array input, this will be used without error if the variable
//             does not exist
//   ier     = nonzero on error
//
//
// Copy the argument to the contained variable.
//
// set<type>(self, name, rank, shape, array, ier)
//   self    = instance
//   name    = name of variable 
//   rank    = rank of array argument
//   shape   = shape of array argument
//   array   = data to send to the contained variable, must have the same rank and
//             extents equal or less then the contained variable
//   ier     = nonzero on error
// 
//
// Creates a new named variable.  There must not already be a variable with that
// name in the container. If inew=0, the argument array is used for the container
// storage for this variable otherwise new storage is allocated and the argument array
// is copied into it.  Attribute functions can be used for the variable right after it has
// been successfully added to the container.
//
// NOTE: because of fortran 90 copying of arguments, the function cvc_hold(array)
//       must be called immediately before this add() function if inew=0.  The purpose is 
//       to temporarily store the pointer to the array in a f77/C way.  The add() method will 
//       then store the pointer to the array in the container along with the type,rank and shape info.
//
// add<type>(self,name, rank, shape, array, inew, ier)
//   self    = instance
//   name    = name of variable 
//   rank    = rank of array argument
//   shape   = shape of array argument
//   array   = the address of this variable will be stored in the container unless
//             inew is zero
//   inew    = nonzero implies a new array of the same shape will be created for storage 
//             zero implies this array address will be used for storage (remember to
//             call cvc_hold(array) first when inew=nonzero)
//   ier     = nonzero on error
//  
// Creates a new named variable filled with zeros or false.
//
// add0<type>(self, name, rank, shape, ier)
//   self    = instance
//   name    = name of variable 
//   rank    = rank of array argument
//   shape   = shape of array argument
//   ier     = nonzero on error
//
//
// An add() will fail if the variable already exists and  a set will fail if the name does not exist.  
// This method combines the add()/remove()/set()/get() for the common operation of setting a variable
// independent of whether it exists or not.  The 'ioption' argument guides
// the operation.  The function returns true if the fwrapPtr was added to the
// container.  The temp_att() attribute is set as in the add method.
//      ioption =  0 -> equivalent to remove() & add(), always returns true
//              =  1 -> if the variable exists with the correct type,rank,shape
//                      then a copy of the data is performed otherwise the pointer
//                      is added to the container possibly replacing an existing pointer.
//              =  2 -> same as 1 but the type can be a copyable_type (float<->double)
//              =  3 -> same as 2 but the contained variable need only agree in copyable_type and
//                      rank and have a large enough shape in each dimension to hold
//                      all of the argument's data
//
// replace<type>(self, name,rank, shape, array, inew, ioption, istatus, ier) 
//   self    = instance
//   name    = name of variable 
//   rank    = rank of array argument
//   shape   = shape of array argument
//   array   = the address of this variable will be stored in the container unless
//             inew is zero
//   inew    = nonzero implies a new array of the same shape will be created for storage 
//             zero implies this array address will be used for storage (remember to
//             call cvc_hold(array) first when inew=nonzero)
//   ioption = see above
//   istatus = nonzero when the array was added, zero on copy (but see ioption==4)
//   ier     = nonzero on error
//  
//
//
// Fetches the address of the fwrap<type> which contains the data for this variable.  The opaque C
// functions in fwrap.h can be used to index into the data directly.
//
// array<type>(self, name, you, ier)
//   self    = instance as a cVarContainer**
//   name    = name of variable 
//   you     = fwrap* is store in this variable
//   ier     = nonzero on error

//
// ----------------- Integer*4 ---------------
//
extern "C" void F77NAME(cvc_getscalarint32)(cVarContainer** self, const char* name, ip_32* value, int* ier) ;

extern "C" void F77NAME(cvc_getdefscalarint32)(cVarContainer** self, const char* name, ip_32* value, const ip_32* def, int* ier) ;

extern "C" void F77NAME(cvc_setscalarint32)(cVarContainer** self, const char* name, const ip_32* value, int* ier);

extern "C" void F77NAME(cvc_addscalarint32)(cVarContainer** self, const char* name, const ip_32* value, int* ier);

extern "C" void F77NAME(cvc_getint32)(cVarContainer** self, const char* name, const int* rank, const int* shape, ip_32* value, 
				      int* ier);

extern "C" void F77NAME(cvc_getdefint32)(cVarContainer** self, const char* name, const int* rank, const int* shape, ip_32* value, 
					 const int* drank, const int* dshape, ip_32* dvalue, int* ier) ;

extern "C" void F77NAME(cvc_setint32)(cVarContainer** self, const char* name, const int* rank, const int* shape, 
				      ip_32* value, int* ier) ;

extern "C" void F77NAME(cvc_addint32)(cVarContainer** self, const char* name, const int* rank, const int* shape, ip_32* avalue, 
				      int* inew, int* ier);

extern "C" void F77NAME(cvc_arrayint32)(cVarContainer** self, const char* name, aip_32** you, int* ier) ;

extern "C" void F77NAME(cvc_add0int32)(cVarContainer** self, const char* name, const int* rank, const int* shape, int* ier) ;
extern "C" void F77NAME(cvc_replaceint32)(cVarContainer** self, const char* name, const int* rank, const int* shape, ip_32* avalue, 
					  int* inew, int* ioption, int* istatus, int* ier);
//
// ----------------- Real*4 ---------------
//
extern "C" void F77NAME(cvc_getscalarfloat32)(cVarContainer** self, const char* name, fp_32* value, int* ier) ;

extern "C" void F77NAME(cvc_getdefscalarfloat32)(cVarContainer** self, const char* name, fp_32* value, const fp_32* 
						 def, int* ier);

extern "C" void F77NAME(cvc_setscalarfloat32)(cVarContainer** self, const char* name, const fp_32* value, int* ier);

extern "C" void F77NAME(cvc_addscalarfloat32)(cVarContainer** self, const char* name, const fp_32* value, int* ier) ;

extern "C" void F77NAME(cvc_getfloat32)(cVarContainer** self, const char* name, const int* rank, const int* shape, fp_32* value, 
					int* ier);

extern "C" void F77NAME(cvc_getdeffloat32)(cVarContainer** self, const char* name, const int* rank, const int* shape, 
					   fp_32* value, const int* drank, const int* dshape, fp_32* dvalue, int* ier);

extern "C" void F77NAME(cvc_setfloat32)(cVarContainer** self, const char* name, const int* rank, const int* shape, fp_32* value, 
					int* ier);

extern "C" void F77NAME(cvc_addfloat32)(cVarContainer** self, const char* name, const int* rank, const int* shape, fp_32* avalue, 
					int* inew, int* ier) ;

extern "C" void F77NAME(cvc_arrayfloat32)(cVarContainer** self, const char* name, afp_32** you, int* ier) ;

extern "C" void F77NAME(cvc_add0float32)(cVarContainer** self, const char* name, const int* rank, const int* shape, int* ier) ;

extern "C" void F77NAME(cvc_replacefloat32)(cVarContainer** self, const char* name, const int* rank, const int* shape, fp_32* avalue, 
					    int* inew, int* ioption, int* istatus, int* ier);
//
// ----------------- Real*8 ---------------
//
extern "C" void F77NAME(cvc_getscalarfloat64)(cVarContainer** self, const char* name, fp_64* value, int* ier);

extern "C" void F77NAME(cvc_getdefscalarfloat64)(cVarContainer** self, const char* name, fp_64* value, const fp_64* def, 
						 int* ier);

extern "C" void F77NAME(cvc_setscalarfloat64)(cVarContainer** self, const char* name, const fp_64* value, int* ier);

extern "C" void F77NAME(cvc_addscalarfloat64)(cVarContainer** self, const char* name, const fp_64* value, int* ier);

extern "C" void F77NAME(cvc_getfloat64)(cVarContainer** self, const char* name, const int* rank, const int* shape, fp_64* value, 
					int* ier);

extern "C" void F77NAME(cvc_getdeffloat64)(cVarContainer** self, const char* name, const int* rank, const int* shape, fp_64* value, 
					   const int* drank, const int* dshape, fp_64* dvalue, int* ier) ;

extern "C" void F77NAME(cvc_setfloat64)(cVarContainer** self, const char* name, const int* rank, const int* shape, 
					fp_64* value, int* ier);

extern "C" void F77NAME(cvc_addfloat64)(cVarContainer** self, const char* name, const int* rank, const int* shape, fp_64* avalue, 
					int* inew, int* ier);

extern "C" void F77NAME(cvc_arrayfloat64)(cVarContainer** self, const char* name, afp_64** you, int* ier) ;

extern "C" void F77NAME(cvc_add0float64)(cVarContainer** self, const char* name, const int* rank, const int* shape, int* ier) ;

extern "C" void F77NAME(cvc_replacefloat64)(cVarContainer** self, const char* name, const int* rank, const int* shape, fp_64* avalue, 
					    int* inew, int* ioption, int* istatus, int* ier);
//
// ----------------- Logical*4 ---------------
//
extern "C" void F77NAME(cvc_settruefalse)(cVarContainer** self, lp_32* truevalue, lp_32* falsevalue, int* ier);

extern "C" void F77NAME(cvc_getscalarlogical32)(cVarContainer** self, const char* name, lp_32* value, int* ier);

extern "C" void F77NAME(cvc_getdefscalarlogical32)(cVarContainer** self, const char* name, lp_32* value, const lp_32* def, 
						   int* ier);

extern "C" void F77NAME(cvc_setscalarlogical32)(cVarContainer** self, const char* name, const lp_32* value, int* ier);

extern "C" void F77NAME(cvc_addscalarlogical32)(cVarContainer** self, const char* name, const lp_32* value, int* ier) ;

extern "C" void F77NAME(cvc_getlogical32)(cVarContainer** self, const char* name, const int* rank, const int* shape, 
					  lp_32* value, int* ier);

extern "C" void F77NAME(cvc_getdeflogical32)(cVarContainer** self, const char* name, const int* rank, const int* shape, 
					     lp_32* value, const int* drank, const int* dshape, lp_32* dvalue, int* ier);

extern "C" void F77NAME(cvc_setlogical32)(cVarContainer** self, const char* name, const int* rank, const int* shape, 
					  lp_32* value, int* ier);

extern "C" void F77NAME(cvc_addlogical32)(cVarContainer** self, const char* name, const int* rank, const int* shape, 
					  lp_32* avalue, int* inew, int* ier);

extern "C" void F77NAME(cvc_arraylogical32)(cVarContainer** self, const char* name, alp_32** you, int* ier) ;

extern "C" void F77NAME(cvc_add0logical32)(cVarContainer** self, const char* name, const int* rank, const int* shape, int* ier) ;

extern "C" void F77NAME(cvc_replacelogical32)(cVarContainer** self, const char* name, const int* rank, const int* shape, lp_32* avalue, 
					      int* inew, int* ioption, int* istatus, int* ier);

//
// ----------------- Character ---------------
// the cvc_*char subroutines assume the length of the strings are in the first dimension of the shape
// and is included in the rank

extern "C" void F77NAME(cvc_getscalarstring)(cVarContainer** self, const char* name, const int* len, cp_char* value, int* ier) ;

extern "C" void F77NAME(cvc_getdefscalarstring)(cVarContainer** self, const char* name, const int* len, cp_char* value, const int* dlen, const cp_char* def, int* ier) ;

extern "C" void F77NAME(cvc_setscalarstring)(cVarContainer** self, const char* name, const int* len, const cp_char* value, int* ier);

extern "C" void F77NAME(cvc_addscalarstring)(cVarContainer** self, const char* name, const int* len, const cp_char* value, int* ier);

extern "C" void F77NAME(cvc_getchar)(cVarContainer** self, const char* name, const int* rank, const int* shape, cp_char* value, 
				      int* ier);

extern "C" void F77NAME(cvc_getdefchar)(cVarContainer** self, const char* name, const int* rank, const int* shape, cp_char* value, 
					 const int* drank, const int* dshape, cp_char* dvalue, int* ier) ;

extern "C" void F77NAME(cvc_setchar)(cVarContainer** self, const char* name, const int* rank, const int* shape, 
				      cp_char* value, int* ier) ;

extern "C" void F77NAME(cvc_addchar)(cVarContainer** self, const char* name, const int* rank, const int* shape, cp_char* avalue, 
				      int* inew, int* ier);

extern "C" void F77NAME(cvc_arraychar)(cVarContainer** self, const char* name, acp_char** you, int* ier) ;

extern "C" void F77NAME(cvc_add0char)(cVarContainer** self, const char* name, const int* rank, const int* shape, int* ier) ;
extern "C" void F77NAME(cvc_replacechar)(cVarContainer** self, const char* name, const int* rank, const int* shape, cp_char* avalue, 
					  int* inew, int* ioption, int* istatus, int* ier);

//
// ----------------- Bytes ---------------
// Bytes is vector<Byte> which acts like a value in the VarContainer.  All operations here
// are restricted to a scalar Bytes and will involve copying.  The olen argument is the number
// of Byte(s) copied from the scalar Bytes.  It is an error if the len argument is smaller than
// the number of Byte(s) stored in the scalar Bytes.

extern "C" void F77NAME(cvc_sizescalarbytes)(cVarContainer** self, const char* name, int* olen, int* ier) ;  // size of scalar Bytes

extern "C" void F77NAME(cvc_getscalarbytes)(cVarContainer** self, const char* name, const int* len, Byte* value, int* olen, int* ier) ;

extern "C" void F77NAME(cvc_getdefscalarbytes)(cVarContainer** self, const char* name, const int* len, Byte* value, const int* dlen, const Byte* def, int* olen, int* ier) ;

extern "C" void F77NAME(cvc_setscalarbytes)(cVarContainer** self, const char* name, const int* len, const Byte* value, int* ier);

extern "C" void F77NAME(cvc_addscalarbytes)(cVarContainer** self, const char* name, const int* len, const Byte* value, int* ier);

extern "C" void F77NAME(cvc_replacescalarbytes)(cVarContainer** self, const char* name, const int* len, const Byte* value, 
                          const int* ioption, int* istatus, int* ier);

extern "C" void F77NAME(cvc_getscalarbytestofile)(cVarContainer** self, const char* name, const char* file, int* ier) ;

extern "C" void F77NAME(cvc_setscalarbytesfromfile)(cVarContainer** self, const char* name, const char* file, int* ier);

extern "C" void F77NAME(cvc_addscalarbytesfromfile)(cVarContainer** self, const char* name, const char* file, int* ier);

extern "C" void F77NAME(cvc_replacescalarbytesfromfile)(cVarContainer** self, const char* name, const char* file, 
                          const int* ioption, int* istatus, int* ier);

extern "C" void F77NAME(cvc_add0bytes)(cVarContainer** self, const char* name, int* ier) ;
