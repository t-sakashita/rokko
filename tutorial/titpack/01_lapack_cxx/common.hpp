/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// C++ version of TITPACK Ver.2 by H. Nishimori

#ifndef TITPACK_COMMON_HPP
#define TITPACK_COMMON_HPP

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::matrix<double> matrix_type;

//
// configurations with the specified sz
//

int sz(int n, double szval, std::vector<int>& list1, std::vector<std::vector<int> >& list2);
// return value @  dimension of the matrix
// n            @  lattice size
// szval        @  total sz
// list1(i)     #  i-th spin configuration
// list2        #  inverse list of list1 expressed by the
//                 2-dim search method of M. Ogata and H.Q. Lin.

//
// data check of pairs of sites
//

void datack(std::vector<int> const& ipair, int n);

//
// eigenvalues by the bisection method
//

void bisec(std::vector<double> const& alpha, std::vector<double> const& beta, int ndim,
           std::vector<double>& E, int ne, double eps,
           int& m, int& nsplit, std::vector<int>& iblock, std::vector<int>& isplit,
           std::vector<double>& w);
// alpha  @ diagonal element
// beta   @ subdiagonal element
// ndim   @ matrix dimension
// E      # eigenvalues
// ne     @ number of eigenvalues to calculate
// eps    @ limit of error

//
// eigenvector of a tridiagonal matrix by inverse iteration for the large/medium routines
// 

void vec12(std::vector<double> const& alpha, std::vector<double> const& beta, int ndim,
           std::vector<double> const& E, matrix_type& z,
           std::vector<int>& iblock, std::vector<int>& isplit, std::vector<double>& w);
// E(4)       @  4 lowest eigenvalues
// ndim       @  matrix dimension
// nvec       @  number of vectors to calculate

//
// xx correlation function
//

void xcorr(int n, std::vector<int> const& npair, matrix_type const& x, int xindex,
           std::vector<double>& sxx, std::vector<int>& list1,
           std::vector<std::vector<int> >& list2);
// n           @ lattice size
// npair       @ pair of sites (k,l) <Sx(k)Sx(l)>
// x           @ eigenvetor
// sxx         # xx correlation function
// list1,list2 @ spin configurations generated in 'sz'

//
// ************* zz correlation function **************
//

void zcorr(int n, std::vector<int> const& npair, matrix_type const& x, int xindex,
           std::vector<double>& szz, std::vector<int>& list1);
// n           @ lattice size
// npair       @ pair of sites (k,l) <Sz(k)Sz(l)>
// x           @ eigenvetor
// szz         # zz correlation function
// list1,list2 @ spin configurations generated in 'sz'

//
// Orthogonalization of the eigenvectors
//

void orthg(matrix_type& ev, std::vector<double>& norm, int& idgn, int numvec);
// ideclr  @  declared array size in the main program
// ev      @# vectors to be orthogonalized / orthogonalized vectors
// norm(j) #  norm of the j-th vector returned
// idgn    #  degree of degenearcy
// numvec  @  number of vectors to be checked

//
// configurations with the specified sz
//

void szdy(int n, int idim, double szval, std::vector<int>& list1,
          std::vector<std::vector<int> >& list2);
// n          @  lattice size
// idim       @  dimension of the matrix
// szval      @  total sz
// list1(i)   #  i-th spin configuration
// list2      #  inverse list of list1 expressed by the
//               2-dim search method of M. Ogata and H.Q. Lin.
// ==============================================================
//      This routine is equivalent to sz but is faster than sz.
//      This routine has been developed by Daijiro Yoshioka,
//      University of Tokyo.  The copyright of szdy belongs to him.
//                                              1993/5/10
// ==============================================================

#endif
