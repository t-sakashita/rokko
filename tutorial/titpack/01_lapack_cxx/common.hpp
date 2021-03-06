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
#include <tuple>

using matrix_type = boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>;
using i_matrix_type = boost::numeric::ublas::matrix<int, boost::numeric::ublas::column_major>;

//
// configurations with the specified sz
//
// return value @  dimension of the matrix
// n            @  lattice size
// szval        @  total sz
// list1(i)     #  i-th spin configuration
// list2        #  inverse list of list1 expressed by the
//                 2-dim search method of M. Ogata and H.Q. Lin.

int sz(int n, double szval, std::vector<int>& list1, std::vector<std::pair<int, int>>& list2);

//
// data check of pairs of sites
//

void datack(std::vector<int> const& ipair, int n);

//
// eigenvalues by the bisection method
//
// return value # m and nsplit
// alpha  @ diagonal element
// beta   @ subdiagonal element
// ndim   @ matrix dimension
// E      # eigenvalues
// ne     @ number of eigenvalues to calculate
// eps    @ limit of error

std::tuple<int, int>
bisec(std::vector<double> const& alpha, std::vector<double> const& beta, int ndim,
      std::vector<double>& E, int ne, double eps, std::vector<int>& iblock,
      std::vector<int>& isplit, double *w);

//
// eigenvector of a tridiagonal matrix by inverse iteration for the large/medium routines
// 
// E(4)       @  4 lowest eigenvalues
// ndim       @  matrix dimension
// nvec       @  number of vectors to calculate

void vec12(std::vector<double> const& alpha, std::vector<double> const& beta, int ndim,
           std::vector<double> const& E, int nvec, matrix_type& z,
           std::vector<int>& iblock, std::vector<int>& isplit, double *w);

//
// xx correlation function
//
// n           @ lattice size
// npair       @ pair of sites (k,l) <Sx(k)Sx(l)>
// x           @ eigenvetor
// sxx         # xx correlation function
// list1,list2 @ spin configurations generated in 'sz'

void xcorr(int n, std::vector<int> const& npair, const double *x,
           std::vector<double>& sxx, std::vector<int>& list1,
           std::vector<std::pair<int, int>>& list2);

void xcorr(int n, std::vector<int> const& npair, std::vector<double> const& x,
           std::vector<double>& sxx, std::vector<int>& list1,
           std::vector<std::pair<int, int>>& list2);

void xcorr(int n, std::vector<int> const& npair, matrix_type const& x, int xindex,
           std::vector<double>& sxx, std::vector<int>& list1,
           std::vector<std::pair<int, int>>& list2);

//
// ************* zz correlation function **************
//
// n           @ lattice size
// npair       @ pair of sites (k,l) <Sz(k)Sz(l)>
// x           @ eigenvetor
// szz         # zz correlation function
// list1,list2 @ spin configurations generated in 'sz'

void zcorr(int n, std::vector<int> const& npair, const double *x,
           std::vector<double>& szz, std::vector<int>& list1);

void zcorr(int n, std::vector<int> const& npair, matrix_type const& x, int xindex,
           std::vector<double>& szz, std::vector<int>& list1);

void zcorr(int n, std::vector<int> const& npair, std::vector<double> const& x,
           std::vector<double>& szz, std::vector<int>& list1);

//
// Orthogonalization of the eigenvectors
//
// return value #  degree of degenearcy
// ev      @# vectors to be orthogonalized / orthogonalized vectors
// norm(j) #  norm of the j-th vector returned
// numvec  @  number of vectors to be checked

int orthg(matrix_type& ev, std::vector<double>& norm, int numvec);

//
// configurations with the specified sz
//
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

void szdy(int n, int idim, double szval, std::vector<int>& list1,
          std::vector<std::pair<int, int>>& list2);

#endif
