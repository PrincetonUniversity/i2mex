//
// template definitions for impreaction.h
//
template<typename FP> void ImpReaction::coronal_ptr(size_t z1, const FP ioniz[], const FP recom[], FP coron[]) {
  size_t z = z1-1 ; // atomic number

  // -- estimate largest populated charge state --
  size_t jc=0 ;   // charge state with largest population
  while(jc<z && (recom[jc]<ioniz[jc])) jc++ ;

  // -- get coronal equilibrium --
  // second pass if we did not estimate jc well
  for (size_t ipass=1 ; ipass<3 ; ipass++) {
    for (size_t j=0 ; j<z1 ; j++) coron[j]=0. ;
    coron[jc] = 1.0 ;    // start at maximum
    
    // -- go up in charge states --
    for (size_t jx=jc+1 ; jx<z1 ; jx++) {
      if (recom[jx]==0.) TRANSP_THROW("?ImpReaction:coronal(): Zero recombination rate (huh?)") ;
      coron[jx] = coron[jx-1]*(ioniz[jx-1]/recom[jx]) ;
      if (coron[jx] < 1.e-10*coron[jc]) break ;                  // avoid underflow
    } ;
    
    // -- go down in charge states --
    for (size_t jx=jc ; jx>0 ; jx--) {
      if (ioniz[jx-1]==0.) TRANSP_THROW("?ImpReaction:coronal(): Zero ionization rate (huh?)") ;
      coron[jx-1] = coron[jx]*(recom[jx]/ioniz[jx-1]) ;
      if (coron[jx] < 1.e-10*coron[jc]) break ;                 // avoid underflow
    } ;
    
    if (ipass==2) break ;
    // --- look for most populated state ---
    FP     zmax = 0. ;  // maximum population
    size_t jmax = 0 ;   // charge state with maximum
    for (size_t jx=0 ; jx<z1 ; jx++) {
      if (coron[jx]>zmax) {
	zmax = coron[jx] ;
	jmax = jx ;
      } ;
    } ;
    //cout << "ipass: " << ipass << "  jc: " << jc << endl ;
    if (jc == jmax) break ;   // already at maximum
    jc = jmax ;               // try it again
  } ;

  // --- normalize ---
  FP sum = 0. ;
  for (size_t jx=0 ; jx<z1 ; jx++) sum+=coron[jx] ;
  for (size_t jx=0 ; jx<z1 ; jx++) coron[jx]/=sum ;
} ;

template<typename FP> void ImpReaction::coronal(const fwrap<FP>& ioniz, const fwrap<FP>& recom, fwrap<FP>& coron) {
  if (ioniz.rank()==1) {
    // ------------ single mode ----------------
    // --- check ---
    if (ioniz.rank()!=1 || recom.rank()!=1 || coron.rank()!=1)
      TRANSP_THROW("?ImpReaction:coronal(): Arguments must be rank 1") ;
  
    size_t z1 = ioniz.xsize(0) ;   // Z+1
    if (z1==0 || recom.xsize(0)!=z1 || coron.xsize(0)!=z1)
      TRANSP_THROW("?ImpReaction:coronal(): Arguments have inconsistent size") ;

    coronal_ptr(z1, &ioniz(0), &recom(0), &coron(0)) ;
  }
  else {
    // ------------ vector mode ----------------
    // --- check ---
    if (ioniz.rank()!=2 || recom.rank()!=2 || coron.rank()!=2)
      TRANSP_THROW("?ImpReaction:coronal(): Arguments must be rank 2") ;
  
    size_t z1 = ioniz.xsize(0) ;   // Z+1
    if (z1==0 || recom.xsize(0)!=z1 || coron.xsize(0)!=z1)
      TRANSP_THROW("?ImpReaction:coronal(): Arguments have inconsistent size") ;
  
    size_t n = ioniz.xsize(1) ;   // number of evaluations
    if (recom.xsize(1)!=n || coron.xsize(1)!=n)
      TRANSP_THROW("?ImpReaction:coronal(): Arguments have inconsistent size") ;

    for (size_t j=0 ; j<n ; j++) {
      coronal_ptr(z1, &ioniz(0,j), &recom(0,j), &coron(0,j)) ;
    } ;
  } ;
} ;

template<typename FP> void ImpReaction::coronal_rad(const fwrap<FP>& coron, const fwrap<FP>& rad, fwrap<FP>& total) {
  if (coron.rank()==1) {
    // -- single mode --
    // --- check ---
    if (coron.rank()!=1 || rad.rank()!=1 || total.rank()!=0)
      TRANSP_THROW("?ImpReaction:coronal_rad(): Arguments must be rank 1, result rank 0") ;
    
    size_t z1 = coron.xsize(0) ;   // Z+1
    if (rad.xsize(0)!=z1) TRANSP_THROW("?ImpReaction:coronal_rad(): Arguments have inconsistent size") ;

    FP result = 0. ;   // final radiation
    FP sum    = 0. ;   // normalize coronal equilibrium
    for (size_t jx=0 ; jx<z1 ; jx++) {
      result += coron(jx)*rad(jx) ;
      sum    += coron(jx) ;
    } ;
    total() = (sum>0.) ? (result/sum) : result ;
  }
  else {
    // -- vector mode --
    // --- check ---
    if (coron.rank()!=2 || rad.rank()!=2 || total.rank()!=1)
      TRANSP_THROW("?ImpReaction:coronal_rad(): Arguments must be rank 1, result rank 0") ;
    
    size_t z1 = coron.xsize(0) ;   // Z+1
    if (rad.xsize(0)!=z1) TRANSP_THROW("?ImpReaction:coronal_rad(): Arguments have inconsistent size") ;

    size_t n = coron.xsize(1) ;  
    if (rad.xsize(1)!=n || total.xsize(0)!=n) 
      TRANSP_THROW("?ImpReaction:coronal_rad(): Arguments have inconsistent size") ;

    for (size_t j=0 ; j<n ; j++) {
      FP result = 0. ;   // final radiation
      FP sum    = 0. ;   // normalize coronal equilibrium
      for (size_t jx=0 ; jx<z1 ; jx++) {
	result += coron(jx,j)*rad(jx,j) ;
	sum    += coron(jx,j) ;
      } ;
      total(j) = (sum>0.) ? (result/sum) : result ;
    } ;
  } ;
} ;


