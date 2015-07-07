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

#include <rokko/parameters.hpp>
#include <rokko/lapack/diagonalize_dsyev.hpp>
#include <rokko/lapack/diagonalize_dsyevd.hpp>
#include <rokko/lapack/diagonalize_dsyevr.hpp>
#include <rokko/lapack/diagonalize_dsyevx.hpp>
#include <rokko/lapack/diagonalize_bisection.hpp>
#include <rokko/lapack/diagonalize_dsygv.hpp>
#include <rokko/lapack/diagonalize_dsygvd.hpp>
#include <rokko/lapack/diagonalize_dsygvx.hpp>
#include <rokko/lapack/diagonalize_bisection_dsygvx.hpp>

namespace rokko {
namespace lapack {

class solver {
public:
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  // standard eigenvalue problem, only eigenvalues, eigenvalue:double*
  template<typename MATRIX_MAJOR>
  void diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
		   double* eigvals,
                   parameters const& params, timer& timer);
  // standard eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:double*
  template<typename MATRIX_MAJOR>  
  void diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
		   double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		   parameters const& params, timer& timer);
  // standard eigenvalue problem, only eigenvalues, eigenvalue:VEC&
  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
		   VEC& eigvals,
                   parameters const& params, timer& timer);
  // standard eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:VEC&
  template<typename MATRIX_MAJOR, typename VEC>  
  void diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
		   VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		   parameters const& params, timer& timer);
  // generalized eigenvalue problem, only eigenvalues, eigenvalue:double*
  template<typename MATRIX_MAJOR>
  void diagonalize(std::string const& routine,
		   localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
		   double* eigvals,
                   parameters const& params, timer& timer);
  // generalized eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:double*
  template<typename MATRIX_MAJOR>  
  void diagonalize(std::string const& routine,
		   localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
		   double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		   parameters const& params, timer& timer);
  // generalized eigenvalue problem, only eigenvalues, eigenvalue:VEC&
  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine,
		   localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
		   VEC& eigvals,
                   parameters const& params, timer& timer);
  // generalized eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:VEC&
  template<typename MATRIX_MAJOR, typename VEC>  
  void diagonalize(std::string const& routine,
		   localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
		   VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		   parameters const& params, timer& timer);
private:
  //std::string routine_;  // record routine_;
};

// -------------------------standard eigenvalue problem-----------------------------

// standard eigenvalue problem, only eigenvalues
template<typename MATRIX_MAJOR>
void solver::diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
			 double* eigvals,
			 parameters const& params, timer& timer) {
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
    if (lapack::is_interval(params)) {
      rokko::lapack::diagonalize_dsyevr(mat, eigvals, params, timer);
    } else {
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
			 parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  timer.stop(timer_id::diagonalize_initialize);
  return solver::diagonalize(routine, mat, &eigvals[0], params, timer);
}

// standard eigenvalue problem, eigenvalues/eigenvectors
template<typename MATRIX_MAJOR>
void solver::diagonalize(std::string const& routine, localized_matrix<double, MATRIX_MAJOR>& mat,
			 double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params, timer& timer) {
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
    if (is_interval(params)) {
      rokko::lapack::diagonalize_dsyevr(mat, eigvals, eigvecs, params, timer);
    } else {
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
			 parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  timer.stop(timer_id::diagonalize_initialize);
  return solver::diagonalize(routine, mat, &eigvals[0], eigvecs, params, timer);
}

// -------------------------generalized eigenvalue problem-----------------------------
// generalized eigenvalue problem, only eigenvalues
template<typename MATRIX_MAJOR>
void solver::diagonalize(std::string const& routine,
			 localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			 double* eigvals,
			 parameters const& params, timer& timer) {
  if ((routine=="dsygv") || (routine=="qr")) {
    rokko::lapack::diagonalize_dsygv(mata, matb, eigvals, params, timer);
  } else if ((routine=="dsygvd") || (routine=="dc")) {
    rokko::lapack::diagonalize_dsygvd(mata, matb, eigvals, params, timer);
  } else if (routine=="dsygvx") {
    rokko::lapack::diagonalize_dsygvx(mata, matb, eigvals, params, timer);
  } else if (routine=="bisection") {
      rokko::lapack::diagonalize_bisection(mata, matb, eigvals, params, timer);
  } else if (routine=="") {
    if (lapack::is_interval(params)) {
      rokko::lapack::diagonalize_dsygvx(mata, matb, eigvals, params, timer);
    } else {
      rokko::lapack::diagonalize_dsygv(mata, matb, eigvals, params, timer);
    }
  } else {
    std::cerr << "error: " << routine << " is not lapack routine" << std::endl;
    throw;
  }
}

template<typename MATRIX_MAJOR, typename VEC>
void solver::diagonalize(std::string const& routine,
			 localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			 VEC& eigvals,
			 parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int dim = mata.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  timer.stop(timer_id::diagonalize_initialize);
  return solver::diagonalize(routine, mata, matb, &eigvals[0], params, timer);
}

// generalized eigenvalue problem, eigenvalues/eigenvectors
template<typename MATRIX_MAJOR>
void solver::diagonalize(std::string const& routine,
			 localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			 double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params, timer& timer) {
  if ((routine=="dsygv") || (routine=="qr")) {
    rokko::lapack::diagonalize_dsygv(mata, matb, eigvals, eigvecs, params, timer);
  } else if ((routine=="dsygvd") || (routine=="dc")) {
    rokko::lapack::diagonalize_dsygvd(mata, matb, eigvals, eigvecs, params, timer);
  } else if (routine=="dsygvx") {
    rokko::lapack::diagonalize_dsygvx(mata, matb, eigvals, eigvecs, params, timer);
  } else if (routine=="bisection") {
      rokko::lapack::diagonalize_bisection(mata, matb, eigvals, eigvecs, params, timer);
  } else if (routine=="") {
    if (is_interval(params)) {
      rokko::lapack::diagonalize_dsygvx(mata, matb, eigvals, eigvecs, params, timer);
    } else {
      rokko::lapack::diagonalize_dsygv(mata, matb, eigvals, eigvecs, params, timer);
    }
  } else {
    std::cerr << "error: " << routine << " is not lapack routine" << std::endl;
    throw;
  }
}

template<typename MATRIX_MAJOR, typename VEC>
void solver::diagonalize(std::string const& routine,
			 localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			 VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int dim = mata.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  timer.stop(timer_id::diagonalize_initialize);
  return solver::diagonalize(routine, mata, matb, &eigvals[0], eigvecs, params, timer);
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_CORE_HPP
