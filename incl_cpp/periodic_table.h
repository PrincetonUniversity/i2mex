#ifndef PERIODIC_TABLE_H
#define PERIODIC_TABLE_H

#include "f77name.h"
#include "transp_util.h"
#include <cstring>
#include <string>
#include <iostream>

//
// ------------ periodic table ------------
//
// declaration of fortran code -- see cperiodic_table.f90 for documentation
//
extern "C" void F77NAME(cto_periodic_table)(const int& ipseudo, const int& z, const fp_64& a, 
					    const int& charge,  const int& isymbol, char* cout) ;

extern "C" void F77NAME(cinv_periodic_table)(const int& ipseudo, const char* symbol, const int& idefault, 
					     int& z, fp_64& a, int& charge) ;

extern "C" void F77NAME(cname_periodic_table)(const int& ipseudo, const int& z, const fp_64& a, char* cout) ;

namespace TRANSP {
  using namespace std;
  // 
  // equivalent (hopefully) to fortran len_trim on C strings
  //
  inline size_t clen_trim(const char* c) {
    for(size_t j=strlen(c) ; j>0 ; j--) if (c[j-1]>' ') return j  ;
    return 0 ;
  } ;

  //
  // convenient c++ forms using strings
  //
  //    int   ipseudo   ! 0 -> standard periodic table
  //                    ! 1 -> use H,D,T for hydrogen isotopes
  //    int   idefault  ! nonzero to return an atomic weight even
  //                    ! if data is not contained in the string
  //    int   isymbol   ! symbol to use to separate element name and charge
  //                    ! 0 -> _  this is an underscore
  //                    ! 1 -> +
  //    int   z         ! atomic number of the element
  //    fp_** a         ! atomic weight of the element
  //    int   charge    ! charge state of the element
  //    string symbol   ! periodic table symbol on input e.g. 'C12+3'
  //

  // --------- float ----------
  inline string sto_periodic_table(const int& ipseudo, const int& z, const fp_32& a, 
				   const int& charge,  const int& isymbol) {
    fp_64 b = a ;
    char cc[13] ;
    F77NAME(cto_periodic_table)(ipseudo, z, b, charge, isymbol, cc) ;
    return string(cc,clen_trim(cc)) ;
  } ;

  inline void sinv_periodic_table(const int& ipseudo, string symbol, const int& idefault, 
				  int& z, fp_32& a, int& charge) {
    fp_64 b = 0. ;
    F77NAME(cinv_periodic_table)(ipseudo, symbol.c_str(), idefault, z, b, charge) ;
    a = b;
  } ;

  inline string sname_periodic_table(const int& ipseudo, const int& z, const fp_32& a) {
    char cc[17] ;
    fp_64 b = a ;
    F77NAME(cname_periodic_table)(ipseudo, z, b, cc) ;
    return string(cc,clen_trim(cc)) ;
  } ;

  // ----------- double ----------
  inline string sto_periodic_table(const int& ipseudo, const int& z, const fp_64& a, 
				   const int& charge,  const int& isymbol) {
    char cc[13] ;
    F77NAME(cto_periodic_table)(ipseudo, z, a, charge, isymbol, cc) ;
    return string(cc,clen_trim(cc)) ;
  } ;

  inline void sinv_periodic_table(const int& ipseudo, string symbol, const int& idefault, 
				  int& z, fp_64& a, int& charge) {
    F77NAME(cinv_periodic_table)(ipseudo, symbol.c_str(), idefault, z, a, charge) ;
  } ;

  inline string sname_periodic_table(const int& ipseudo, const int& z, const fp_64& a) {
    char cc[17] ;
    F77NAME(cname_periodic_table)(ipseudo, z, a, cc) ;
    return string(cc,clen_trim(cc)) ;
  } ;

} ;

#endif
