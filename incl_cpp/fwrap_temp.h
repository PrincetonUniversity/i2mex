
//
// template definitions for fwrap.h
//
#include <iomanip> 

using namespace std ;

// ------------------ fwrap -------------------
template<typename FP> size_t fwrap<FP>::total_check() const {
  size_t total = 1 ;        // accumulate total size
  for (size_t i=0 ; i<prank ; i++) {
    if (dim[i]>store[i]) TRANSP_THROW("? fwrap::fwrap(): Logical storage size exceeds physical size") ;
    total *= store[i] ;
  } ;
  return total ;
} ;

template<typename FP> fwrap<FP>::fwrap(FP* fp) : 
  myptr(fp==0), prank(0), dim(), store(), p(fp), count(0), atts(0) {
  size_t total = 1 ;   // accumulate total size
  if (myptr) p = new FP[total] ;
} ;

template<typename FP> fwrap<FP>::fwrap(dimsize d1, FP* fp) : 
  myptr(fp==0), prank(1), dim(1), store(1), p(fp), count(0), atts(0) {
  dim[0]   = d1.logical ;
  store[0] = d1.physical ;
  size_t total = total_check() ;   // accumulate total size
  if (myptr) p = new FP[total] ;
} ;

template<typename FP> fwrap<FP>::fwrap(size_t d1, FP* fp) : 
  myptr(fp==0), prank(1), dim(1), store(1), p(fp), count(0), atts(0) {
  dim[0]   = d1 ;
  store[0] = d1 ;
  if (myptr) {
    size_t total = total_check() ;   // accumulate total size
    p = new FP[total] ;
  } ;
} ;

template<typename FP> fwrap<FP>::fwrap(dimsize d1, dimsize d2, FP* fp) : 
  myptr(fp==0), prank(2), dim(2), store(2), p(fp), count(0), atts(0) {
  dim[0]   = d1.logical ;
  dim[1]   = d2.logical ;
  store[0] = d1.physical ;
  store[1] = d2.physical ;
  size_t total = total_check() ;   // accumulate total size
  if (myptr) p = new FP[total] ;
} ;

template<typename FP> fwrap<FP>::fwrap(size_t d1, size_t d2, FP* fp) : 
  myptr(fp==0), prank(2), dim(2), store(2), p(fp), count(0), atts(0) {
  dim[0]   = d1 ;
  dim[1]   = d2 ;
  store[0] = d1 ;
  store[1] = d2 ;
  if (myptr) {
    size_t total = total_check() ;   // accumulate total size
    p = new FP[total] ;
  } ;
} ;

template<typename FP> fwrap<FP>::fwrap(const Shape& s, FP* fp) : 
  myptr(fp==0), prank(s.size()), dim(s), store(s), p(fp), count(0), atts(0) {
  size_t total = total_check() ;   // accumulate total size
  if (myptr) p = new FP[total] ;
} ;

template<typename FP> fwrap<FP>::fwrap(const Shape& logical, const Shape& physical, FP* fp) : 
  myptr(fp==0), prank(logical.size()), dim(logical), store(physical), p(fp), count(0), atts(0) {
  if (logical.size()!=physical.size()) TRANSP_THROW("?fwrap::fwrap(): Mismatch in size of logical and physical shapes") ;
  size_t total = total_check() ;   // accumulate total size
  if (myptr) p = new FP[total] ;
} ;

template<typename FP> fwrap<FP>::fwrap(const fwrap<FP>& f) : myptr(true), prank(f.prank),
							     dim(f.dim), store(f.store), p(0), count(0), atts(0) {
  size_t total = total_check() ;
  p = new FP[total] ;
  assign(f) ;
  if (f.atts!=0) {
    atts = new fwrap_map_att() ;
    atts->operator=(*f.atts) ;
  } ;
} ;

template<typename FP> fwrap<FP>& fwrap<FP>::operator=(const fwrap<FP>& f) {
  if (&f != this) {
    size_t total = f.pstore() ;    // size of new memory
    if ((!myptr) || (pstore() != total)) {
      // need to (maybe) free old memory and get new memory
      if (myptr) delete[] p ;
      p = new FP[total] ;
      myptr = true ;
    } ;
    prank = f.prank ;
    dim   = f.dim ;
    store = f.store ;
    FP* ip = p ;    // iterate from beginning
    FP* jp = f.p ;
    for (size_t i=0 ; i<total ; i++) *(ip++)=*(jp++) ;  // copy
    if (f.atts==0) {
      if (atts!=0) atts->clear() ;
    }
    else {
      if (atts==0) atts = new fwrap_map_att() ;
      atts->operator=(*f.atts) ;
    } ;
  } ;
  return *this ;
} ;

template<typename FP> void fwrap<FP>::assign(const fwrap<FP>& f) {
  //
  // note: the physical storage location and therefore myptr does not
  // change in this assignment.
  //
  if (&f != this) {
    size_t total = pstore() ;
    if (total != f.pstore()) TRANSP_THROW("?fwrap::assign(): Mismatch in size of fwrap storage in assignment") ;
    FP* ip = p ;    // iterate from beginning
    FP* jp = f.p ;
    for (size_t i=0 ; i<total ; i++) *(ip++)=*(jp++) ;  // copy
  } ;
} ;

template<typename FP> void fwrap<FP>::assign(FP f) {
  FP* ip = p ;    // iterate from beginning
  size_t total = pstore() ;
  for (size_t i=0 ; i<total ; i++) *(ip++)=f ;  // copy
} ;   

template<typename FP> fwrap<FP>::~fwrap() {
  //cout << "kaboom" << this << endl ;
  if (myptr) delete[] p ;   // only delete internally generated memory
  if (atts) {
    atts->clear() ;
    delete atts ;
  } ;
} ; 

template<typename FP>  size_t fwrap<FP>::size(size_t i) const {
  if (i>=prank) TRANSP_THROW("?fwrap::size(): Argument>=rank of array") ;
  return dim[i] ;
} ;

template<typename FP>  size_t fwrap<FP>::size() const {
  size_t total = 1 ;
  for (size_t i=0 ; i<prank ; i++) total *= dim[i] ;
  return total ;
} ;

template<typename FP>  size_t fwrap<FP>::pstore(size_t i) const {
  if (i>=prank) TRANSP_THROW("?fwrap::pstore(): Argument>=rank of array") ;
  return store[i] ;
} ;

template<typename FP>  size_t fwrap<FP>::pstore() const {
  size_t total = 1 ;
  for (size_t i=0 ; i<prank ; i++) total *= store[i] ;
  return total ;
} ;

template<typename FP> FP& fwrap<FP>::operator[](const Shape& ix) {
  if (ix.size() != prank) TRANSP_THROW("?fwrap::operator[]: Argument has the wrong rank") ;
  bool   ok = true ;      // index error flag
  size_t k  = 0 ;         // accumulated linear index
  switch (prank) {
  case 0:
    k = 0 ;               // faster 0d
    break ;
  case 1:
    ok = ix[0]<dim[0] ;   // faster 1d
    k  = ix[0] ;
    break ;
  case 2:
    ok = (ix[0]<dim[0]) && (ix[1]<dim[1]) ;  // faster 2d
    k  = ix[0] + store[0]*ix[1] ;
    break ;
  default:
    size_t i  = prank-1 ;   // dimension
    while(ok) {             // arbitrary d
      ok = ix[i]<dim[i] ;
      k += ix[i] ;
      if (i==0) break ;
      k *= store[--i] ;
    } ;
  } ;
  if (!ok) TRANSP_THROW("?fwrap::operator[]: Index exceeded logical dimension") ;
  return p[k] ;
} ;

template<typename FP> const FP& fwrap<FP>::operator[](const Shape& ix) const {
  if (ix.size() != prank) TRANSP_THROW("?fwrap::operator[]: Argument has the wrong rank") ;
  bool   ok = true ;      // index error flag
  size_t k  = 0 ;         // accumulated linear index
  switch (prank) {
  case 0:
    k = 0 ;               // faster 0d
    break ;
  case 1:
    ok = ix[0]<dim[0] ;   // faster 1d
    k  = ix[0] ;
    break ;
  case 2:
    ok = (ix[0]<dim[0]) && (ix[1]<dim[1]) ;  // faster 2d
    k  = ix[0] + store[0]*ix[1] ;
    break ;
  default:
    size_t i  = prank-1 ;   // dimension
    while(ok) {             // arbitrary d
      ok = ix[i]<dim[i] ;
      k += ix[i] ;
      if (i==0) break ;
      k *= store[--i] ;
    } ;
  } ;
  if (!ok) TRANSP_THROW("?fwrap::operator[]: Index exceeded logical dimension") ;
  return p[k] ;
} ;

template<typename FP>  void fwrap<FP>::set_size(size_t i, size_t n) {
  if (i>=prank)   TRANSP_THROW("?fwrap::set_size(): Argument>=rank of array") ;
  if (n>store[i]) TRANSP_THROW("?fwrap::set_size(): New logical dimension is out of bounds") ;
  dim[i] = n ;
} ;

template<typename FP> fwrap<FP> fwrap<FP>::last_slice(size_t i) {
  if (prank<1) TRANSP_THROW("?fwrap::last_slice(): Can not slice a rank 0 array") ;

  size_t nrank = prank-1 ;  // rank of new fwrap
  if (i>=dim[nrank]) TRANSP_THROW("?fwrap::last_slice(): Slice index is out of bounds") ;

  Shape logic(nrank), physic(nrank) ;  // new logical and physical storage
  size_t k=i ;                         // start of last slice = i*store[prank-2]*store[prank-3]...*store[0]
  for (size_t j=0 ; j<nrank ; j++) {
    k*=store[j] ;
    logic[j]  = dim[j] ;
    physic[j] = store[j] ;
  } ;

  return fwrap(logic,physic,p+k) ;
} ;

template<typename FP> size_t fwrap<FP>::att_size() const {
  if (atts==0) return 0 ;
  return atts->size() ;
} ;     

template<typename FP> fwrap_map_att* fwrap<FP>::att() {
  if (atts==0) atts = new fwrap_map_att() ;
  return atts ;
} ;     

template<typename FP> const fwrap_map_att* fwrap<FP>::att() const {
  if (atts==0) atts = new fwrap_map_att() ;
  return atts ;
} ;     


template<typename FP> ostream& operator<<(ostream& s, const fwrap<FP>& f) {
  // -- write rank and indices ---
  size_t rank = f.rank() ;
  s << "rank: " << rank << endl ;
  for (size_t i=0 ; i<rank ; i++) s << "dim[" << setw(2) << i << "]:   " << setw(6) << f.size(i) 
				    << "  " << setw(6) << f.pstore(i) << endl ;
  if (rank==0) {
    s << f() << endl ;
  }
  else {
    Shape loop(rank) ;     // initialize to zero -- represents loop index
#ifndef NO_GOT_VALARRAY
    loop = 0 ;
#else
    set_vector_int(loop,0) ;
#endif
    size_t k= f.size()>0 ? 0 : rank  ;   // dimension being incremented
    while(k<rank) {                      // terminate when try to increment nonexistent k=rank dimension
      // -- loop body --
      for (size_t i=0 ; i<rank ; i++) s << setw(2) << loop[i] << "  " ; // indices
      s << "  " << f[loop] << endl ;                                    // value
      
      // -- loop management --
      // increment loop index -- when hit the dimension limit
      // for dimension k, that dimension is set to zero and
      // the next dimension needs incrementing.  Store f.size(*)
      // in a Shape for time critical loops.
      //
      k = 0 ;  // start out incrementing smallest dimension
      while (k<rank && (++loop[k] >= f.size(k))) loop[k++]=0 ;
    } ;
  } ;
  
  if (f.att_size()>0) f.att()->write(s, string(" & ")) ;
  return s ;
} ;

template<> ostream& operator<<(ostream& s, const fwrap<Bytes>& f) ;
template<> ostream& operator<<(ostream& s, const fwrap<cp_char>& f) ;

template<typename FP> ostream& operator<<(ostream& s, const fwrapPtr<FP>& f) {
  if (f.isnull()) {
    s << "null fwrapPtr" ;
  }
  else {
    s << *f ;
  } ;
  return s ;
} ;

template<typename FP, typename GP> void fwrapCopy(const fwrap<FP>& f, fwrap<GP>& g, bool needall) {
  size_t rank = f.rank() ;
  if (g.rank()!=rank) TRANSP_THROW("?fwrapCopy: argument ranks are different") ;

  size_t gsize = g.size() ;
  size_t fsize = f.size() ;

  if (gsize==0 || fsize==0) {
    if (needall && fsize!=0) TRANSP_THROW("?fwrapCopy: need all data but target fwrap has zero size") ;
    return ;
  } ;

  const FP* fptr = f.fptr() ;
  GP*       gptr = g.fptr() ;
  if (rank==0) {
    *gptr = static_cast<GP>(*fptr) ;
    return ;
  } ;

  bool  oquick        = true ;       // true when a quick copy is possible
  const Shape& fshape = f.shape()  ; // logical  shape of f
  const Shape& pshape = f.pshape() ; // physical shape of f
  const Shape& gshape = g.pshape() ; // physical shape of g
  Shape tshape(rank) ;               // only loop over these extents

  size_t kmax=1 ;                    // total physical size of f

  size_t fstride[rank] ;
  size_t gstride[rank] ;

  for (size_t i=0 ; i<rank ; i++) {
    size_t fs = fshape[i] ;
    size_t ps = pshape[i] ;
    size_t gs = gshape[i] ;
    if (needall && fs>gs) TRANSP_THROW("?fwrapCopy: need all data but source fwrap has a larger dimension than the target") ;
    oquick = oquick && (fs==ps)&&(fs==gs) ;
    kmax *= ps ;
    tshape[i] = (fs>gs) ? gs : fs ;  // minimum
    if (tshape[i]==0) return ;       // should have been already caught
    if (i==0) {
      fstride[0]=1 ;
      gstride[0]=1 ;
    }
    else {
      size_t ip=i-1 ;
      fstride[i]=fstride[ip]*pshape[ip] ;
      gstride[i]=gstride[ip]*gshape[ip] ;
    } ;
  } ;

  if (oquick) {
    // the logical and physical shape of f is the same as the physical shape of g
    // so the arrays are identical in memory layout
    //cout << "doing a quick copy with rank " << rank << " and " << kmax << " = " << f.pstore() << endl; 
    for (size_t k=0 ; k<kmax ; k++) *gptr++ = static_cast<GP>(*fptr++) ;
    return ;
  } ;
  //cout << "no quickness" << endl ;

  Shape loop(rank) ;
  set_vector_int(loop,0) ;

  size_t jf=0 ;  // index into f
  size_t jg=0 ;  // index into g
  while(1) {
    gptr[jg] = static_cast<GP>(fptr[jf]) ;
    
    size_t k=0 ;
    while(k<rank) {
      size_t ii = loop[k] ;
      if (++loop[k]>=tshape[k]) {
	jf -= fstride[k]*ii ;
	jg -= gstride[k]*ii ;
	loop[k]=0 ;
	k++ ;
      }
      else {
	jf += fstride[k] ;
	jg += gstride[k] ;
	break ;
      } ;
    } ;
    if (k>=rank) break ;
  } ;
} ;


template<typename FP, typename GP> void fwrapCopy_bool(const fwrap<FP>& f, fwrap<GP>& g, GP truevalue, GP falsevalue, 
						       bool needall) {
  size_t rank = f.rank() ;
  if (g.rank()!=rank) TRANSP_THROW("?fwrapCopy: argument ranks are different") ;

  size_t gsize = g.size() ;
  size_t fsize = f.size() ;

  if (gsize==0 || fsize==0) {
    if (needall && fsize!=0) TRANSP_THROW("?fwrapCopy: need all data but target fwrap has zero size") ;
    return ;
  } ;

  const FP* fptr = f.fptr() ;
  GP* gptr = g.fptr() ;
  if (rank==0) {
    *gptr = (*fptr) ? truevalue : falsevalue ;
    return ;
  } ;

  const Shape& fshape = f.shape()  ;
  const Shape& pshape = f.pshape() ;
  const Shape& gshape = g.pshape() ;
  Shape tshape(rank) ;              // only loop over these extents

  size_t fstride[rank] ;
  size_t gstride[rank] ;

  for (size_t i=0 ; i<rank ; i++) {
    if (needall && fshape[i]>gshape[i]) TRANSP_THROW("?fwrapCopy: need all data but source fwrap has a larger dimension than the target") ;
    tshape[i] = (fshape[i]>gshape[i]) ? gshape[i] : fshape[i] ;  // minimum
    if (tshape[i]==0) return ;
    if (i==0) {
      fstride[0]=1 ;
      gstride[0]=1 ;
    }
    else {
      size_t ip=i-1 ;
      fstride[i]=fstride[ip]*pshape[ip] ;
      gstride[i]=gstride[ip]*gshape[ip] ;
    } ;
  } ;

  Shape loop(rank) ;
  set_vector_int(loop,0) ;

  size_t jf=0 ;  // index into f
  size_t jg=0 ;  // index into g

  while(1) {
    gptr[jg] = (fptr[jf]) ? truevalue : falsevalue ;
    
    size_t k=0 ;
    while(k<rank) {
      size_t ii = loop[k] ;
      if (++loop[k]>=tshape[k]) {
	jf -= fstride[k]*ii ;
	jg -= gstride[k]*ii ;
	loop[k]=0 ;
	k++ ;
      }
      else {
	jf += fstride[k] ;
	jg += gstride[k] ;
	break ;
      } ;
    } ;
    if (k>=rank) break ;
  } ;
} ;

template<typename FP> Shape fwrapCharShape(const fwrap<FP>& f) {
  return Shape(f.shape()) ;
}

template<> Shape fwrapCharShape(const fwrap<cp_char>& f);

template<typename FP> size_t fwrapCharLength(const fwrap<FP>& f) {  // only relevant for FP=char which is specialized
  return 1 ;
}

template<> size_t fwrapCharLength(const fwrap<cp_char>& f) ;

template<typename FP, typename GP> void fwrapCharCopy(const fwrap<FP>& f, fwrap<GP>& g, const Shape& s) {
  g[s] = f[s] ;
}

template<typename FP> size_t fwrapCharMaxLength(const fwrap<FP>& f) {
  return 1 ;
} ;

template<> void fwrapCharCopy(const fwrap<cp_char>& f, fwrap<cp_char>& g, const Shape& s) ;
template<> void fwrapCharCopy(const fwrap<cp_char>& f, fwrap<string>&  g, const Shape& s) ;
template<> void fwrapCharCopy(const fwrap<string>&  f, fwrap<cp_char>& g, const Shape& s) ;

template<> size_t fwrapCharMaxLength(const fwrap<cp_char>& f) ;
template<> size_t fwrapCharMaxLength(const fwrap<string>& f) ;

template<typename FP, typename GP> void fwrapCopy_char(const fwrap<FP>& f, fwrap<GP>& g, bool needall) {
  Shape fshape = fwrapCharShape(f) ;
  Shape gshape = fwrapCharShape(g) ;

  size_t rank = fshape.size() ;

  if (gshape.size()!=rank) TRANSP_THROW("?fwrapCopy_char: argument character ranks are different") ;

  size_t flen = fwrapCharLength(f) ;  // length of character element if fixed otherwise nonzero
  size_t glen = fwrapCharLength(g) ; 

  if (glen==0 || flen==0) return ;    // nothing to do

  size_t fsize = 1 ;    // number of character string elements
  size_t gsize = 1 ; 
  for (size_t i=0 ; i<rank ; i++) {
    fsize *= fshape[i] ;
    gsize *= gshape[i] ;
  }

  if (gsize==0 || fsize==0) {
    if (needall && fsize!=0) TRANSP_THROW("?fwrapCopy_char: need all data but target fwrap has zero size") ;
    return ;
  } ;

  Shape tshape(rank) ;        // only loop over these extents

  for (size_t i=0 ; i<rank ; i++) {
    size_t fs = fshape[i] ;
    size_t gs = gshape[i] ;

    if (needall && fs>gs) TRANSP_THROW("?fwrapCopy_char: need all data but source fwrap has a larger dimension than the target") ;

    tshape[i] = (fs>gs) ? gs : fs ;  // minimum
    if (tshape[i]==0) return ;       // already tested, loop below will fail if this condition is not met
  } ;

  Shape loop(rank) ;         // running indices to fill
  set_vector_int(loop,0) ;

  while(1) {
    fwrapCharCopy(f, g, loop) ;
    
    size_t k=0 ;
    while(k<rank) {
      if (++loop[k]>=tshape[k]) {
	loop[k]=0 ;
	k++ ;
      }
      else {
	break ;
      } ;
    } ;
    if (k>=rank) break ;
  } ;
}



template<typename FP> fwrapPtr<FP> fwrapPtr<FP>::build(FP* fp) {
  return fwrapPtr(new fwrap<FP>(fp)) ;
} ;

template<typename FP> fwrapPtr<FP> fwrapPtr<FP>::build(dimsize d, FP* fp){
  return fwrapPtr(new fwrap<FP>(d,fp)) ;
} ;

template<typename FP> fwrapPtr<FP> fwrapPtr<FP>::build(size_t  d, FP* fp){
  return fwrapPtr(new fwrap<FP>(d,fp)) ;
} ;

template<typename FP> fwrapPtr<FP> fwrapPtr<FP>::build(dimsize d, dimsize e, FP* fp){
  return fwrapPtr(new fwrap<FP>(d,e,fp)) ;
} ;

template<typename FP> fwrapPtr<FP> fwrapPtr<FP>::build(size_t  d, size_t  e, FP* fp){
  return fwrapPtr(new fwrap<FP>(d,e,fp)) ;
} ;

template<typename FP> fwrapPtr<FP> fwrapPtr<FP>::build(const Shape& s, FP* fp) {
  return fwrapPtr(new fwrap<FP>(s,fp)) ;
} ;

template<typename FP> fwrapPtr<FP> fwrapPtr<FP>::build(const Shape& logical, const Shape& physical, FP* fp){
  return fwrapPtr(new fwrap<FP>(logical,physical,fp)) ;
} ;

template<typename FP> fwrapPtr<FP> fwrapPtr<FP>::build(const fwrap<FP>& p){
  return fwrapPtr(new fwrap<FP>(p)) ;
} ;

