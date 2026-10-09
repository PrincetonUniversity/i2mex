#ifndef VARCONTAINER_IO_H
#define VARCONTAINER_IO_H

#include "varcontainer.h"   

//
// --- VarTbuf<T> ---
// acts as an expanding buffer which will be freed on destruction
template<class T> class VarTbuf {
public:
    VarTbuf() : buf(0), nbuf(0) {} ;
    ~VarTbuf() ;
      
    void need(size_t n) ;  // expand buffer if necessary for this many elements
      
    T*     buf ;   // the buffer
    size_t nbuf ;  // number of elements currently allocated
} ;
        
//
// ------------------ VarWorkerWriter -------------------
// Helper for writing out actual data from a VarContainer.  If you want
// to handle endian issues do it here.
// This is used internally by the VarContainer. The writer->done()/error() will be
// called if necessary on destruction.
//
class VarWorkerWriter {
public:
    VarWorkerWriter(string container, VarWriter* writer) ;  // construct with name of container and writer
    ~VarWorkerWriter() ;
    
    //
    // VarWriter methods
    //
    void begin() ;        // writes out VERSION<n> after calling writer->begin()
    
    inline void send(const Byte* p, size_t nbytes) { writer->send(p,nbytes) ; }
    
    void done() ;             // writes out DONE and then calls writer->done()
    
    void error(string msg) ;  // calls writer->error(msg)
    
    //
    // these are more structured writes of data which includes the size of
    // of the object being sent
    //
    template<typename T> void write(const T*, size_t nt) ;

    void write(const string*, size_t nstring) ;
    void write(const Shape*,  size_t nshape) ;
    void write(const Bytes*,  size_t nbytes) ;
    
    void write(const fwrap_map_att* atts) ;
    
    //
    // write out the data and attributes of the pointer
    //
    template<typename FP> void writeFwrap(const fwrapPtr<FP>&) ;

    void writeFwrapBytes(const fwrapPtr<Bytes>&) ;  
    void writeFwrapLogical(const fwrapPtr<lp_32>&) ;  // write as sequence of 0 or 1 bytes
    
    static const int VERSION  ;  // version of data format
private:
    VarWorkerWriter(const VarWorkerWriter&) ;
    VarWorkerWriter& operator=(const VarWorkerWriter&) ;
    
    // --- data --
    string     container ;   // name of the container being operated on
    VarWriter* writer ;
    bool       wasbegun ;    // true if begin() was called but done()/error() has not been called
    
    VarTbuf<size_t> bsizet ;
    VarTbuf<Byte>   bbyte ;      
} ;


//
// ------------------ VarWorkerReader -------------------
// Helper for reading data for a VarContainer.  If you want
// to handle endian issues do it here.  This must be compatible with VarWorkerWriter.
// This is used internally by the VarContainer.  The reader->done()/error() will be
// called if necessary on destruction.
//
class VarWorkerReader {
public:
    VarWorkerReader(string container, VarReader* reader) ;  // construct with name of container and reader
    ~VarWorkerReader() ;
       
    //
    // VarReader methods
    //
    void begin() ;  // will check against VarWorkerWriter::VERSION
    
    inline void receive(Byte* p, size_t nbytes) { reader->receive(p,nbytes) ; } ;  
                                          
    void done() ; 

    void error(string msg) ; 

    //
    // read the next set of data as nt elements of type T.  The pointer should
    // already have been allocated. This is NOT the inverse of 
    // template<typename T> void VarWorkerWriter::write(T*, size_t)
    // which first sends the number of elements before sending the data.
    //
    template<typename T> void read(T* p, size_t nt) ;

    //
    // these read<type> methods expect to first read a single size_t containing the 
    // number of elements and then will read the elements into an internal
    // temporary buffer which can be fetched with ptr<type>.  Returns the number of
    // elements read.
    //
    size_t  readSizet() ;  
    size_t* ptrSizet() { return bsizet.buf ; } ;
    
    size_t  readByte() ;  
    Byte*   ptrByte() { return bbyte.buf ; } ;
    
    size_t  readInt() ;  
    int*    ptrInt() { return bint.buf ; } ;
    
    size_t  readDouble() ;  
    double* ptrDouble() { return bdouble.buf ; } ;
    
    size_t  readString() ;  
    string* ptrString() { return bstring.buf ; } ;
    
    size_t  readShape() ;  
    Shape*  ptrShape() { return bshape.buf ; } ;
    
    size_t  readBytes() ;  
    Bytes*  ptrBytes() { return bbytes.buf ; } ;
    
    void read(fwrap_map_att* atts) ;  // read in the attributes and set them in the argument map

    //
    // read the data and attributes of the pointer
    //
    template<typename FP> fwrapPtr<FP> readFwrap() ;   
    
    fwrapPtr<Bytes> readFwrapBytes() ;  
    fwrapPtr<lp_32> readFwrapLogical(lp_32 lp_true, lp_32 lp_false) ;  // set to values considered true and false
 
    void expect(size_t x, size_t t) ;  // throw an exception if the sizes are not the same
                                       // used when expect x elements but t elements were sent   
private:
    VarWorkerReader(const VarWorkerReader&) ;
    VarWorkerReader& operator=(const VarWorkerReader&) ;    
    
     
    void toobig(size_t g, size_t b, string msg) ;  // throw an exception if g>b with the given message,
                                                   // used to help prevent allocating ridiculously large data
                                                    
    // --- data ---
    string     container ;   // name of the container being operated on
    VarReader* reader ;
    bool       wasbegun ;    // true if begin() was called but done()/error() was not called
    
    VarTbuf<size_t> bsizet ;    
    VarTbuf<int>    bint ;    
    VarTbuf<double> bdouble ;    
    VarTbuf<string> bstring ;    
    VarTbuf<Shape>  bshape ;    
    VarTbuf<Bytes>  bbytes ;  
    
    VarTbuf<char>   bchar ;  // for strings
    VarTbuf<Byte>   bbyte ;  // for Bytes
} ;

//
// ============================ VarReader's and VarWriter's =============================
//
// --------------- VarFileWriter ------------------
// VarWriter for writing to a file.  File will be opened with begin() and closed with done() or error()
// or destruction.
//
class VarFileWriter : public VarWriter {
public:
    VarFileWriter(string fname) ;  // name of the file
    virtual ~VarFileWriter() ;     // will close open file if necessary
    
    virtual void begin() ;         // opens the file
    
    virtual void send(const Byte* p, size_t nbytes) ; 
    
    virtual void done()  ;   // will close the open file
    
    virtual void error(string msg) ;   // will write out message to cerr and then close the open file  
private:
    string    fname ;   // name of the file 
    ofstream* fs ;      // filled when opened  
} ;

//
// --------------- VarFileReader ------------------
// VarReader for reading from a file.  File will be opened with begin() and closed with done() or error()
// or destruction.
//
class VarFileReader : public VarReader {
public:
    VarFileReader(string fname) ;  // name of the file
    virtual ~VarFileReader() ;     // will close open file if necessary
    
    virtual void begin() ;         // opens the file
    
    virtual void receive(Byte* p, size_t nbytes) ; 
    
    virtual void done()  ;   // will close the open file
    
    virtual void error(string msg) ;   // will write out message to cerr and then close the open file  
private:
    string    fname ;   // name of the file 
    ifstream* fs ;      // filled when opened  
} ;

#endif
