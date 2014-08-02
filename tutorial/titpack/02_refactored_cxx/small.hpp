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
#include "hamiltonian.hpp"
#include <vector>

//
// matrix elements
//

void elm3(hamiltonian const& hop, matrix_type& elemnt);

//
// eigenvalues of a small matrix
//
// elemnt      @ matrix elements
// E           # eigenvalues
// v           # eigenvector
// ne          @ number of eigenvalues to calculate
// nvec        @ number of eigenvectors to calculate

void diag(matrix_type& elemnt, std::vector<double>& E, matrix_type& v, int nvec);
//
// check of the eigenvector and eigenvalue
//
// elemnt    @ nonzero elements
// x         @ eigenvector to be checked
// xindex
// return value: Hexpec <x*H*x>

double check3(matrix_type const& elemnt, matrix_type const& x, int xindex);

#endif
