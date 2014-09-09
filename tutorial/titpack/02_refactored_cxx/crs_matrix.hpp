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

#ifndef TITPACK_CRS_MATRIX_HPP
#define TITPACK_CRS_MATRIX_HPP

#include "common.hpp"
#include "hamiltonian.hpp"

//
// crs_matrix
//

class crs_matrix {
public:
  crs_matrix(hamiltonian const& hop);
  int dimension() const { return elemnt_.size1(); }
  double multiply(const double *v1, double *v0) const;
private:
  matrix_type elemnt_;
  i_matrix_type loc_;
};

#endif
