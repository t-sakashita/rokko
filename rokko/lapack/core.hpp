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
  parameters diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat,
			 double* eigvals,
			 parameters const& params);
  // standard eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:double*
  template<typename MATRIX_MAJOR>  
  parameters diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat,
			 double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params);
  // standard eigenvalue problem, only eigenvalues, eigenvalue:VEC&
  template<typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 parameters const& params);
  // standard eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:VEC&
  template<typename MATRIX_MAJOR, typename VEC>  
  parameters diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params);
  // generalized eigenvalue problem, only eigenvalues, eigenvalue:double*
  template<typename MATRIX_MAJOR>
  parameters diagonalize(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			 double* eigvals,
			 parameters const& params);
  // generalized eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:double*
  template<typename MATRIX_MAJOR>  
  parameters diagonalize(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			 double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params);
  // generalized eigenvalue problem, only eigenvalues, eigenvalue:VEC&
  template<typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			 VEC& eigvals,
			 parameters const& params);
  // generalized eigenvalue problem, eigenvalues/eigenvectors, eigenvalue:VEC&
  template<typename MATRIX_MAJOR, typename VEC>  
  parameters diagonalize(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			 VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params);
private:
  //std::string routine_;  // record routine_;
};

// -------------------------standard eigenvalue problem-----------------------------
// standard eigenvalue problem, only eigenvalues
template<typename MATRIX_MAJOR>
parameters solver::diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat,
			       double* eigvals,
			       parameters const& params) {
  parameters params_out;
  std::string routine = "";
  if(params.defined("routine")) {
    routine = params.get_string("routine");
  }
  if ((routine=="dsyev") || (routine=="qr")) {
    return rokko::lapack::diagonalize_dsyev(mat, eigvals, params);
  } else if ((routine=="dsyevr") || (routine=="mr3")) {
    rokko::lapack::diagonalize_dsyevr(mat, eigvals, params);
  } else if ((routine=="dsyevd") || (routine=="dc")) {
    return rokko::lapack::diagonalize_dsyevd(mat, eigvals, params);
  } else if (routine=="dsyevx") {
    return rokko::lapack::diagonalize_dsyevx(mat, eigvals, params);
  } else if (routine=="bisection") {
    return rokko::lapack::diagonalize_bisection(mat, eigvals, params);
  } else if (routine=="") {
    if (lapack::is_interval(params)) {
      return rokko::lapack::diagonalize_dsyevr(mat, eigvals, params);
    } else {
      return rokko::lapack::diagonalize_dsyev(mat, eigvals, params);
    }
  } else {
    std::cerr << "error: " << routine << " is not lapack routine" << std::endl;
    throw;
  }
  return params_out;
}

template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  return solver::diagonalize(mat, &eigvals[0], params);
}

// standard eigenvalue problem, eigenvalues/eigenvectors
template<typename MATRIX_MAJOR>
parameters solver::diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat,
			       double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  parameters params_out;
  std::string routine = "";
  if(params.defined("routine")) {
    routine = params.get_string("routine");
  }
  if ((routine=="dsyev") || (routine=="qr")) {
    return rokko::lapack::diagonalize_dsyev(mat, eigvals, eigvecs, params);
  } else if ((routine=="dsyevr") || (routine=="mr3")) {
    return rokko::lapack::diagonalize_dsyevr(mat, eigvals, eigvecs, params);
  } else if ((routine=="dsyevd") || (routine=="dc")) {
    return rokko::lapack::diagonalize_dsyevd(mat, eigvals, eigvecs, params);
  } else if (routine=="dsyevx") {
    return rokko::lapack::diagonalize_dsyevx(mat, eigvals, eigvecs, params);
  } else if (routine=="bisection") {
    return rokko::lapack::diagonalize_bisection(mat, eigvals, eigvecs, params);
  } else if (routine=="") {
    if (is_interval(params)) {
      return rokko::lapack::diagonalize_dsyevr(mat, eigvals, eigvecs, params);
    } else {
      return rokko::lapack::diagonalize_dsyev(mat, eigvals, eigvecs, params);
    }
  } else {
    std::cerr << "error: " << routine << " is not lapack routine" << std::endl;
    throw;
  }
  return params_out;
}

template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  return solver::diagonalize(mat, &eigvals[0], eigvecs, params);
}


// -------------------------generalized eigenvalue problem-----------------------------
// generalized eigenvalue problem, only eigenvalues
template<typename MATRIX_MAJOR>
parameters solver::diagonalize(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			       double* eigvals,
			       parameters const& params) {
  parameters params_out;
  std::string routine = "";
  if(params.defined("routine")) {
    routine = params.get_string("routine");
  }
  if ((routine=="dsygv") || (routine=="qr")) {
    return rokko::lapack::diagonalize_dsygv(mata, matb, eigvals, params);
  } else if ((routine=="dsygvd") || (routine=="dc")) {
    return rokko::lapack::diagonalize_dsygvd(mata, matb, eigvals, params);
  } else if (routine=="dsygvx") {
    return rokko::lapack::diagonalize_dsygvx(mata, matb, eigvals, params);
  } else if (routine=="") {
    if (lapack::is_interval(params)) {
      return rokko::lapack::diagonalize_dsygvx(mata, matb, eigvals, params);
    } else {
      return rokko::lapack::diagonalize_dsygv(mata, matb, eigvals, params);
    }
  } else {
    std::cerr << "error: " << routine << " is not lapack routine" << std::endl;
    throw;
  }
  return params_out;
}

template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			       VEC& eigvals,
			       parameters const& params) {
  int dim = mata.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  return solver::diagonalize(mata, matb, &eigvals[0], params);
}

// generalized eigenvalue problem, eigenvalues/eigenvectors
template<typename MATRIX_MAJOR>
parameters solver::diagonalize(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			       double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  parameters params_out;
  std::string routine = "";
  if(params.defined("routine")) {
    routine = params.get_string("routine");
  }
  if ((routine=="dsygv") || (routine=="qr")) {
    return rokko::lapack::diagonalize_dsygv(mata, matb, eigvals, eigvecs, params);
  } else if ((routine=="dsygvd") || (routine=="dc")) {
    return rokko::lapack::diagonalize_dsygvd(mata, matb, eigvals, eigvecs, params);
  } else if (routine=="dsygvx") {
    return rokko::lapack::diagonalize_dsygvx(mata, matb, eigvals, eigvecs, params);
  } else if (routine=="") {
    if (is_interval(params)) {
      return rokko::lapack::diagonalize_dsygvx(mata, matb, eigvals, eigvecs, params);
    } else {
      return rokko::lapack::diagonalize_dsygv(mata, matb, eigvals, eigvecs, params);
    }
  } else {
    std::cerr << "error: " << routine << " is not lapack routine" << std::endl;
    throw;
  }
  return params_out;
}

template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			       VEC& eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  int dim = mata.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  return solver::diagonalize(mata, matb, &eigvals[0], eigvecs, params);
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_CORE_HPP
