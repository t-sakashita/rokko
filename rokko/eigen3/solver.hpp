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

#pragma once

#include <rokko/eigen3/diagonalize.hpp>

namespace rokko {
namespace eigen3 {

class solver {
public:
  void initialize(int& /* argc */, char**& /* argv */) {}
  void finalize() {}
  // -------------------------standard eigenvalue problem-----------------------------
  // eigenvalues/eigenvectors
  template<typename T, int MATRIX_MAJOR, typename VEC>
  parameters diagonalize(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			 rokko::parameters const& params) {
    std::string routine = params.defined("routine") ? params.get_string("routine") : "";

    if ((routine == "") || (routine == "qr")) {
      return rokko::eigen3::diagonalize(mat, eigvals, eigvecs, params);
    } else {
      throw std::invalid_argument("eigen3::diagonalize() : " + routine + " is invalid routine name");
    }
  }
  // only eigenvalues
  template<typename T, int MATRIX_MAJOR, typename VEC>
  parameters diagonalize(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 rokko::parameters const& params) {
    std::string routine = params.defined("routine") ? params.get_string("routine") : "";

    if ((routine == "") || (routine == "qr")) {
      return rokko::eigen3::diagonalize(mat, eigvals, params);
    } else {
      throw std::invalid_argument("eigen3::diagonalize() : " + routine + " is invalid routine name");
    }
  }
};

} // namespace eigen3
} // namespace rokko
