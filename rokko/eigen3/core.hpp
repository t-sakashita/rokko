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
  // -------------------------standard eigenvalue problem-----------------------------
  // eigenvalues/eigenvectors
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(localized_matrix<T, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 localized_matrix<T, MATRIX_MAJOR>& eigvecs,
			 rokko::parameters const& params) {
    std::string routine = "";
    if(params.defined("routine")) {
      routine = params.get_string("routine");
    }
    if ((routine == "") || (routine == "qr")) {
      rokko::eigen3::diagonalize(mat, eigvals, eigvecs, params);
    } else {
      std::cerr << "error: " << routine << " is not eigen3 routine" << std::endl;
      throw;
    }
  }
  // only eigenvalues
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(localized_matrix<T, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 rokko::parameters const& params) {
    std::string routine = "";
    if(params.defined("routine")) {
      routine = params.get_string("routine");
    }
    if ((routine == "") || (routine == "qr")) {
	rokko::eigen3::diagonalize(mat, eigvals, params);
    } else {
      std::cerr << "error: " << routine << " is not eigen3 routine" << std::endl;
      throw;
    }
  }
  // -------------------------generalized eigenvalue problem-----------------------------
  // eigenvalues/eigenvectors
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(localized_matrix<T, MATRIX_MAJOR>& mata, localized_matrix<T, MATRIX_MAJOR>& matb,
			 VEC& eigvals,
			 localized_matrix<T, MATRIX_MAJOR>& eigvecs,
			 rokko::parameters const& params) {  
    std::cerr << "error: " << "eigen3 does not have routine for generalized eigenvalue problem" << std::endl;
    throw;
  }
  // only eigenvalues
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(localized_matrix<T, MATRIX_MAJOR>& mata, localized_matrix<T, MATRIX_MAJOR>& matb,
			 VEC& eigvals,
			 rokko::parameters const& params) {
    std::cerr << "error: " << "eigen3 does not have routine for generalized eigenvalue problem" << std::endl;
    throw;
  }
};

} // namespace eigen3
} // namespace rokko

#endif // ROKKO_EIGEN3_CORE_HPP
