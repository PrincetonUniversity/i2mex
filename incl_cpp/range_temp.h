
//
// template definitions for range.h
//
#include <iomanip> 


template<typename FP> bool Range::check(const fwrap<FP>& x, fwrap<FP>& f) const {
  FP my_min = minimum ;        // save on repetitious conversions
  FP my_max = maximum ;

  size_t frank = f.rank() ;    // rank of f
  size_t xrank = x.rank() ;    // must be 1
  if (xrank != 1) TRANSP_THROW("?preact::Range::check(): Axis must be a rank 1 fwrap object") ;
  if (frank < 1)  TRANSP_THROW("?preact::Range::check(): Result must be at least a rank 1 fwrap object") ;

  size_t nx = x.size(0) ;      // size of data to operate on
  if (nx != f.size(frank-1)) TRANSP_THROW("?preact::Range::check(): Axis dimension and last function dimension must be the same size") ;

  unsigned int under = 0 ;  // count number of out of bounds requiring action
  unsigned int over  = 0 ;

  for (size_t i=0 ; i<nx ; i++) {
    if (x(i) < my_min) {
      if (minimum_policy == IGNORE_CLEAR || minimum_policy == WARN_CLEAR)
	if (frank==1) { f(i)=clear_minimum ; } else { f.last_slice(i).assign(clear_minimum) ; } ;
      if (minimum_policy != IGNORE_CLEAR && minimum_policy != IGNORE) under++ ;
    }
    else if (x(i) > my_max) {
      if (maximum_policy == IGNORE_CLEAR || maximum_policy == WARN_CLEAR)
	if (frank==1) { f(i)=clear_maximum ; } else { f.last_slice(i).assign(clear_maximum) ; } ;
      if (maximum_policy != IGNORE_CLEAR && maximum_policy != IGNORE) over++ ;
    } ;
  } ;

  if (under>0) {
    cerr << "?preact::Range::check(): " << under << " out-of-range low value(s) in '"
         << label << "' evaluation." << endl;
    cerr << "   axis minimum is: " << minimum << endl;
  } ;

  if (over>0) {
    cerr << "?preact::Range::check(): " << over << " out-of-range high value(s) in '"
         << label << "' evaluation." << endl;
    cerr << "   axis maximum is: " << maximum << endl;
  } ;

  if ((under>0 && minimum_policy == HALT) ||
      (over >0 && maximum_policy == HALT)   ) TRANSP_THROW("?preact::Range::check(): Axis bounds exceeded with HALT policy -- " + label) ;

  return (over>0 || under>0) ;
} ;


template<typename FP> bool Range::check(const FP& x, fwrap<FP>& f) const {
  unsigned int under = 0 ;  // count number of out of bounds requiring action
  unsigned int over  = 0 ;

  if (x < minimum) {
    if (minimum_policy == IGNORE_CLEAR || minimum_policy == WARN_CLEAR)
      f.assign(clear_minimum) ;              // clear function values
    if (minimum_policy != IGNORE_CLEAR && minimum_policy != IGNORE) under++ ;
  }
  else if (x > maximum) {
    if (maximum_policy == IGNORE_CLEAR || maximum_policy == WARN_CLEAR)
      f.assign(clear_maximum) ;              // clear function values
    if (maximum_policy != IGNORE_CLEAR && maximum_policy != IGNORE) over++ ;
  } ;

  if (under>0) {
    cerr << "?preact::Range::check(): " << under << " out-of-range low value(s) in '"
	 << label << "' evaluation." << endl;
    cerr << "   axis minimum is: " << minimum << endl;
  } ;

  if (over>0) {
    cerr << "?preact::Range::check(): " << over << " out-of-range high value(s) in '"
	 << label << "' evaluation." << endl;
    cerr << "   axis maximum is: " << maximum << endl;
  } ;

  if ((under>0 && minimum_policy == HALT) ||
      (over >0 && maximum_policy == HALT)   ) TRANSP_THROW("?preact::Range::check(): Axis bounds exceeded with HALT policy -- " + label) ;
  return (over>0 || under>0) ;
} ;

template<typename FP> bool Range::check(const FP& x, FP& f) const {
  unsigned int under = 0 ;  // count number of out of bounds requiring action
  unsigned int over  = 0 ;

  if (x < minimum) {
    if (minimum_policy == IGNORE_CLEAR || minimum_policy == WARN_CLEAR)
      f=clear_minimum ;              // clear function values
    if (minimum_policy != IGNORE_CLEAR && minimum_policy != IGNORE) under++ ;
  }
  else if (x > maximum) {
    if (maximum_policy == IGNORE_CLEAR || maximum_policy == WARN_CLEAR)
      f=clear_maximum ;              // clear function values
    if (maximum_policy != IGNORE_CLEAR && maximum_policy != IGNORE) over++ ;
  } ;

  if (under>0) {
    cerr << "?preact::Range::check(): " << under << " out-of-range low value(s) in '"
	 << label << "' evaluation." << endl;
    cerr << "   axis minimum is: " << minimum << endl;
  } ;

  if (over>0) {
    cerr << "?preact::Range::check(): " << over << " out-of-range high value(s) in '"
	 << label << "' evaluation." << endl;
    cerr << "   axis maximum is: " << maximum << endl;
  } ;

  if ((under>0 && minimum_policy == HALT) ||
      (over >0 && maximum_policy == HALT)   ) TRANSP_THROW("?preact::Range::check(): Axis bounds exceeded with HALT policy -- " + label) ;
  return (over>0 || under>0) ;
} ;
