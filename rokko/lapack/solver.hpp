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

#ifndef ROKKO_LAPACK_SOLVER_HPP
#define ROKKO_LAPACK_SOLVER_HPP

#include <rokko/parameters.hpp>
#include <rokko/lapack/diagonalize_syev.hpp>
#include <rokko/lapack/diagonalize_syevd.hpp>
#include <rokko/lapack/diagonalize_syevr.hpp>
#include <rokko/lapack/diagonalize_syevx.hpp>
#include <rokko/lapack/diagonalize_bisection.hpp>
#include <rokko/lapack/diagonalize_qr.hpp>
#include <rokko/lapack/diagonalize_sygv.hpp>
#include <rokko/lapack/diagonalize_sygvd.hpp>
#include <rokko/lapack/diagonalize_sygvx.hpp>
#include <rokko/lapack/diagonalize_bisection_dsygvx.hpp>

namespace rokko {
namespace lapack {

class solver {
public:
  void initialize(int& /* argc */, char**& /* argv */) {}
  void finalize() {}

  // standard eigenvalue problem, only eigenvalues, eigenvalue:VEC&
  template<typename T, int MATRIX_MAJOR, typename VEC>
  parameters diagonalize(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 parameters const& params);
  // standard eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:VEC&
  template<typename T, int MATRIX_MAJOR, typename VEC>
  parameters diagonalize(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat,
			 VEC& eigvals, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			 parameters const& params);

private:
  //std::string routine_;  // record routine_;
};

// -------------------------standard eigenvalue problem-----------------------------
// standard eigenvalue problem, only eigenvalues
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat,
                   VEC& eigvals,
			       parameters const& params) {
  const std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if (routine=="syev") {
    return rokko::lapack::diagonalize_syev(mat, eigvals, params);
  } else if ((routine=="syevr") || (routine=="mr3")) {
    return rokko::lapack::diagonalize_syevr(mat, eigvals, params);
  } else if ((routine=="syevd") || (routine=="dc")) {
    return rokko::lapack::diagonalize_syevd(mat, eigvals, params);
  } else if (routine=="syevx") {
    return rokko::lapack::diagonalize_syevx(mat, eigvals, params);
  } else if (routine=="bisection") {
    return rokko::lapack::diagonalize_bisection(mat, eigvals, params);
  } else if (routine=="qr") {
    return rokko::lapack::diagonalize_qr(mat, eigvals, params);
  } else if (routine=="") {
    if (lapack::is_interval(params)) {
      return rokko::lapack::diagonalize_syevr(mat, eigvals, params);
    } else {
      return rokko::lapack::diagonalize_syev(mat, eigvals, params);
    }
  } else {
    throw std::invalid_argument("lapack::diagonalize() : " + routine + " is not lapack routine");
  }
}

// standard eigenvalue problem, eigenvalues/eigenvectors
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat,
                               VEC& eigvals, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  const std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if (routine=="syev") {
    return rokko::lapack::diagonalize_syev(mat, eigvals, eigvecs, params);
  } else if ((routine=="syevr") || (routine=="mr3")) {
    return rokko::lapack::diagonalize_syevr(mat, eigvals, eigvecs, params);
  } else if ((routine=="syevd") || (routine=="dc")) {
    return rokko::lapack::diagonalize_syevd(mat, eigvals, eigvecs, params);
  } else if (routine=="syevx") {
    return rokko::lapack::diagonalize_syevx(mat, eigvals, eigvecs, params);
  } else if (routine=="bisection") {
    return rokko::lapack::diagonalize_bisection(mat, eigvals, eigvecs, params);
  } else if (routine=="qr") {
    return rokko::lapack::diagonalize_qr(mat, eigvals, eigvecs, params);
  } else if (routine=="") {
    if (is_interval(params)) {
      return rokko::lapack::diagonalize_syevr(mat, eigvals, eigvecs, params);
    } else {
      return rokko::lapack::diagonalize_syev(mat, eigvals, eigvecs, params);
    }
  } else {
    throw std::invalid_argument("lapack::diagonalize() : " + routine + " is not lapack routine");
  }
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_SOLVER_HPP
