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

#ifndef ROKKO_LAPACK_CORE_HPP
#define ROKKO_LAPACK_CORE_HPP

#include <rokko/lapack/diagonalize_dsyev.hpp>
#include <rokko/lapack/diagonalize_dsyevd.hpp>
#include <rokko/lapack/diagonalize_dsyevr.hpp>
#include <rokko/lapack/diagonalize_dsyevx.hpp>
#include <rokko/lapack/diagonalize_bisection.hpp>

namespace rokko {
namespace lapack {

class solver {
public:
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  template<typename MATRIX_MAJOR>
  void diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
		   double* eigvals,
                   rokko::parameters const& params, timer& timer);
  template<typename MATRIX_MAJOR>  
  void diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
		   double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		   rokko::parameters const& params, timer& timer);
  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
		   VEC& eigvals,
                   rokko::parameters const& params, timer& timer);
  template<typename MATRIX_MAJOR, typename VEC>  
  void diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
		   VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		   rokko::parameters const& params, timer& timer);
private:
  //std::string routine_;  // record routine_;
};

// only eigenvalues
template<typename MATRIX_MAJOR>
void solver::diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
			 double* eigvals,
			 rokko::parameters const& params, timer& timer) {
  if ((routine=="dsyev") || (routine=="qr")) {
    rokko::lapack::diagonalize_dsyev(mat, eigvals, params, timer);
  } else if ((routine=="dsyevr") || (routine=="mr3")) {
    rokko::lapack::diagonalize_dsyevr(mat, eigvals, params, timer);
  } else if ((routine=="dsyevd") || (routine=="dc")) {
    rokko::lapack::diagonalize_dsyevd(mat, eigvals, params, timer);
  } else if (routine=="dsyevx") {
    rokko::lapack::diagonalize_dsyevx(mat, eigvals, params, timer);
  } else if (routine=="bisection") {
    rokko::lapack::diagonalize_bisection(mat, eigvals, params, timer);
  } else if (routine=="") {
    if (params.defined("upper_limit") || params.defined("lower_limit")) {
      rokko::lapack::diagonalize_dsyevr(mat, eigvals, params, timer);
    }
    else {
      rokko::lapack::diagonalize_dsyev(mat, eigvals, params, timer);
    }
  } else {
    std::cerr << "error: " << routine << " is not lapack routine" << std::endl;
    throw;
  }
}

template<typename MATRIX_MAJOR, typename VEC>
void solver::diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 rokko::parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  timer.stop(timer_id::diagonalize_initialize);
  return solver::diagonalize(routine, mat, &eigvals[0], params, timer);
}

// eigenvalues/eigenvectors
template<typename MATRIX_MAJOR>
void solver::diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
			 double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 rokko::parameters const& params, timer& timer) {
  if ((routine=="dsyev") || (routine=="qr")) {
    rokko::lapack::diagonalize_dsyev(mat, eigvals, eigvecs, params, timer);
  } else if ((routine=="dsyevr") || (routine=="mr3")) {
    rokko::lapack::diagonalize_dsyevr(mat, eigvals, eigvecs, params, timer);
  } else if ((routine=="dsyevd") || (routine=="dc")) {
    rokko::lapack::diagonalize_dsyevd(mat, eigvals, eigvecs, params, timer);
  } else if (routine=="dsyevx") {
    rokko::lapack::diagonalize_dsyevx(mat, eigvals, eigvecs, params, timer);
  } else if (routine=="bisection") {
    rokko::lapack::diagonalize_bisection(mat, eigvals, eigvecs, params, timer);
  } else if (routine=="") {
    if (params.defined("upper_limit") || params.defined("lower_limit")) {
      rokko::lapack::diagonalize_dsyevr(mat, eigvals, eigvecs, params, timer);
    }
    else {
      rokko::lapack::diagonalize_dsyev(mat, eigvals, eigvecs, params, timer);
    }
  } else {
    std::cerr << "error: " << routine << " is not lapack routine" << std::endl;
    throw;
  }
}

template<typename MATRIX_MAJOR, typename VEC>
void solver::diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 rokko::parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  timer.stop(timer_id::diagonalize_initialize);
  return solver::diagonalize(routine, mat, &eigvals[0], eigvecs, params, timer);
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_CORE_HPP
