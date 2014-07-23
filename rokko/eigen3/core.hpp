/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Synge Todo <wistaria@comp-phys.org>,
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN3_CORE_HPP
#define ROKKO_EIGEN3_CORE_HPP

#include <rokko/eigen3/diagonalize.hpp>
#include <iostream>

namespace rokko {
namespace eigen3 {

class solver {
public:
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  template<typename MATRIX_MAJOR, typename TIMER>
  void diagonalize(localized_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                   localized_matrix<MATRIX_MAJOR>& eigvecs, TIMER& timer_in) {
    rokko::eigen3::diagonalize(mat, eigvals, eigvecs, timer_in);
  }
};

} // namespace eigen3
} // namespace rokko


#endif // ROKKO_EIGEN3_CORE_HPP
