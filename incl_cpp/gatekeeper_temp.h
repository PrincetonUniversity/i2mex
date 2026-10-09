
//
// template definitions for gatekeeper.h
//
template<class FP> GateKeeper<FP>::GateKeeper(string name, size_t n) : cname(name), nstart(n), keys(n), ghosts() {
  if (n==0) throw GateKeeperError(cname, "not allowed to construct an empty GateKeeper") ;
  for (size_t i=0 ; i<n ; i++) {
    keys[i]=0 ;
    ghosts.push_front(i) ;
  } ;
  check() ;
} ;

template<class FP> GateKeeper<FP>::~GateKeeper() {
  for (size_t i=0 ; i<keys.size() ; i++) {
    //if (keys[i]!=0) cout << "GateKeeper<FP>::~GateKeeper() deleting " << keys[i]->first << endl ;  // debug
    delete keys[i] ;
  } ;
} ;

template<class FP>  FP GateKeeper<FP>::getPoint(size_t i) const {
  if (i>=keys.size()) throw GateKeeperError(cname,"?getPoint: attempt to index past end of array") ;
  value* p = keys[i] ;
  if (p==0) throw GateKeeperError(cname, "?getPoint: indexed into empty slot of the GateKeeper") ;
  return p->second ;
};
  
template<class FP>  FP GateKeeper<FP>::getPoint(size_t i, string checkname) const { 
  if (i>=keys.size()) throw GateKeeperError(cname, "?getPoint: attempt to index past end of array") ;
  value* p = keys[i] ;
  if (p==0) throw GateKeeperError(cname, "?getPoint: indexed into empty slot of the GateKeeper") ;
  if (checkname != p->first) throw GateKeeperError(cname, "?getPoint: string check failed when indexing into the GateKeeper, want '"
						  +checkname+"' but have '"+p->first+"'") ;
  return p->second ;
} ;
                                               
template<class FP>  string GateKeeper<FP>::getName(size_t i) const {
  if (i>=keys.size()) throw GateKeeperError(cname, "?getName: attempt to index past end of array") ;
  value* p = keys[i] ;
  if (p==0) throw GateKeeperError(cname, "?getName: indexed into empty slot of the GateKeeper") ;
  return p->first ;
} ;     

template<class FP> size_t GateKeeper<FP>::add(string name, FP point) { 
  if (point->rank()>GATEKEEPER_MAX_RANK) throw GateKeeperError(cname,"?add: attempt to add an array with rank>GATEKEEPR_MAX_RANK") ;
  if (ghosts.empty()) resize(4*keys.size()) ;
  size_t i = ghosts.back() ;
  if (i>=keys.size()) throw GateKeeperError(cname,"?add: logic error, integer slot number exceeds key vector size") ;
  if (keys[i]!=0) throw GateKeeperError(cname, "?add: logic error, slot unexpectedly contains data with string"+keys[i]->first) ;
  ghosts.pop_back() ;
  keys[i] = new value(name,point) ;
  return i ;
} ;    

template<class FP>  void GateKeeper<FP>::remove(size_t i) {
  if (i>=keys.size()) throw GateKeeperError(cname, "?remove: attempt to index past end of array") ;
  value* p = keys[i] ;
  if (p==0) throw GateKeeperError(cname, "?remove: indexed into empty slot of the GateKeeper") ;
  delete p ;
  keys[i]=0 ;
  ghosts.push_back(i) ;   // return index so it can be reused
} ;
         
template<class FP>  void GateKeeper<FP>::remove(size_t i, string checkname) {
  if (i>=keys.size()) throw GateKeeperError(cname, "?remove: attempt to index past end of array") ;
  value* p = keys[i] ;
  if (p==0) throw GateKeeperError(cname, "?remove: indexed into empty slot of the GateKeeper") ;
  if (checkname != p->first) throw GateKeeperError(cname, "?remove: string check failed when indexing into the GateKeeper, want '"
						  +checkname+"' but have '"+p->first+"'") ;
  delete p ;
  keys[i]=0 ;
  ghosts.push_back(i) ;   // return index so it can be reused
} ;
         
template<class FP>  void GateKeeper<FP>::clear() {
  ghosts.clear() ;
  for (size_t i=0 ; i<keys.size() ; i++) {
    delete keys[i] ;
    keys[i]=0 ;
    ghosts.push_front(i) ;
  } ;
} ;                  

template<class FP>  size_t GateKeeper<FP>::size()     const { 
  size_t k = keys.size() ;
  size_t g = ghosts.size() ;
  if (k<g) throw GateKeeperError(cname, "?size(): logic error, the ghosts has too many indices") ;
  return (k-g) ;
} ;       

template<class FP>  size_t GateKeeper<FP>::capacity() const {
  return keys.size() ;
} ;       

template<class FP> void GateKeeper<FP>::check() const {
  size_t n = keys.size() ;
  for (list<size_t>::const_iterator q = ghosts.begin() ; q!=ghosts.end() ; ++q) {
    size_t i = *q ;
    if (i>=n) throw GateKeeperError(cname, "?check(): ghosts contains an index which is too large") ;
    if (keys[i]!=0) throw GateKeeperError(cname, "?check(): failed check, index in ghosts points to nonempty slot") ;
  } ;
  size_t nc=0 ;     // number of nonempty slots in keys
  for (size_t i=0 ; i<keys.size() ; i++) {
    if (keys[i]!=0) nc++ ;
  } ;
  if ((ghosts.size()+nc) != keys.size()) throw GateKeeperError(cname, "?check(): there is a slot containing zero which should be occupied") ;
} ;       

template<class FP>  void GateKeeper<FP>::resize(size_t nsize) {
  check() ;
  size_t n = keys.size() ;
  if (nsize<=n) return ;   // do nothing
  keys.resize(nsize,0) ;   // increase with empty slots
  for (size_t i=n ; i<keys.size() ; i++) ghosts.push_front(i) ;
  check() ;
} ;

template<class FP> void GateKeeper<FP>::write(ostream& s) const {
  check() ;
  s << "----------------- GateKeeper ---------------" << endl ;
  s << " Name:      " << cname << endl ;
  s << " Size:      " << size() << endl ;
  s << " Available: " << (capacity()-size()) << endl ; 
  s << endl ;
  s << " Keys: -->" << endl ;
  for (size_t i=0 ; i<keys.size() ; i++) {
    value* p = keys[i] ;
    if (p) s << " " << p->first << endl ;
  } ;
  s << "--------------- End GateKeeper -------------" << endl ;
} ;

template<class FP> list<size_t> GateKeeper<FP>::indices() const {
  list<size_t> ii ;
  for (size_t i=0 ; i<keys.size() ; i++) {
    if (keys[i]!=0) ii.push_back(i) ;
  } ;
  return ii ;
} ;

