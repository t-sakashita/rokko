/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAGMA_CORE_HPP
#define ROKKO_MAGMA_CORE_HPP

#include <rokko/magma/diagonalize.hpp>
#include <iostream>

namespace rokko {
namespace magma {

struct dsyev {};
struct dsyevd {};
struct dsyevx {};

template<typename ROUTINE>
class solver {
public:
  void initialize(int& argc, char**& argv) {}

  void finalize() {}

  void diagonalize(localized_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                   localized_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in);

  void diagonalize(localized_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                   localized_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in);
};

template<>
inline void solver<rokko::magma::dsyev>::diagonalize(localized_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                                                     localized_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in) {
  rokko::magma::diagonalize(mat, eigvals, eigvecs, timer_in);
}

template<>
inline void solver<rokko::magma::dsyev>::diagonalize(localized_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                                                     localized_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in) {
  rokko::magma::diagonalize(mat, eigvals, eigvecs, timer_in);
}

} // namespace magma
} // namespace rokko


#endif // ROKKO_MAGMA_CORE_HPP

