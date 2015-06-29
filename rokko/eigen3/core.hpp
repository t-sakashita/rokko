/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN3_CORE_HPP
#define ROKKO_EIGEN3_CORE_HPP

#include <rokko/eigen3/diagonalize.hpp>

namespace rokko {
namespace eigen3 {

class solver {
public:
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  // eigenvalues/eigenvectors
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine, localized_matrix<T, MATRIX_MAJOR>& mat,
		   VEC& eigvals,
		   localized_matrix<T, MATRIX_MAJOR>& eigvecs,
		   rokko::parameters const& params, timer& timer) {
    if ((routine == "") || (routine == "qr")) {
      rokko::eigen3::diagonalize(mat, eigvals, eigvecs, params, timer);
    } else {
      std::cerr << "error: " << routine << " is not eigen3 routine" << std::endl;
      throw;
    }
  }
  // only eigenvalues
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine, localized_matrix<T, MATRIX_MAJOR>& mat,
		   VEC& eigvals,
		   rokko::parameters const& params, timer& timer) {
    if ((routine == "") || (routine == "qr")) {
      rokko::eigen3::diagonalize(mat, eigvals, params, timer);
    } else {
      std::cerr << "error: " << routine << " is not eigen3 routine" << std::endl;
      throw;
    }
  }
};

} // namespace eigen3
} // namespace rokko

#endif // ROKKO_EIGEN3_CORE_HPP
