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

#ifndef TITPACK_SMALL_HPP
#define TITPACK_SMALL_HPP

#include "common.hpp"
#include <vector>

//
// matrix elements
//

void elm3(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
          std::vector<double> const& zrtio, matrix_type& elemnt,
          std::vector<int> const& list1, std::vector<std::pair<int, int> > const& list2);
// n           @ lattice size
// ipair       @ pairs of sites connected by bonds
// bondwt      @ exchange interaction of each bond Jxy
// zrtio       @ ratio of Jz to Jxy
// elemnt      # matrix elementsle
// list1,list2 @ spin configurations generated in 'sz'

//
// eigenvalues of a small matrix
//

void diag(matrix_type& elemnt, std::vector<double>& E, matrix_type& v, int nvec);
// elemnt      @ matrix elements
// E           # eigenvalues
// v           # eigenvector
// ne          @ number of eigenvalues to calculate
// nvec        @ number of eigenvectors to calculate

//
// check of the eigenvector and eigenvalue
//

double check3(matrix_type const& elemnt, matrix_type const& x, int xindex);
// elemnt    @ nonzero elements
// x         @ eigenvector to be checked
// xindex
// return value: Hexpec <x*H*x>

#endif
