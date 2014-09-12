/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_PARALLEL_SPARSE_SOLVER_HPP
#define ROKKO_PARALLEL_SPARSE_SOLVER_HPP

#include <rokko/factory.hpp>
#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>

#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {

namespace detail {
    
class ps_solver_base {
public:
  virtual ~ps_solver_base() {}
  virtual void initialize(int& argc, char**& argv) = 0;
  virtual void finalize() = 0;
  virtual void diagonalize(rokko::distributed_crs_matrix& mat,
			   int num_evals, int block_size, int max_iters, double tol, timer& timer) = 0;
  virtual void diagonalize(rokko::distributed_crs_matrix& mat,
			   int num_evals, int block_size, int max_iters, double tol) = 0;
  virtual void diagonalize(rokko::distributed_mfree& mat,
			   int num_evals, int block_size, int max_iters, double tol, timer& timer) = 0;
  virtual void diagonalize(rokko::distributed_mfree& mat,
			   int num_evals, int block_size, int max_iters, double tol) = 0;
  virtual double eigenvalue(int i) const = 0;
  virtual void eigenvector(int i, std::vector<double>& vec) const = 0;
  virtual int num_conv() const = 0;
  virtual rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim, int col_dim) = 0;
};

template<typename SOLVER>
class ps_solver_wrapper : public ps_solver_base {
  typedef SOLVER solver_type;
public:
  ps_solver_wrapper() : solver_impl_() {}
  virtual ~ps_solver_wrapper() {}
  void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
  void finalize() { solver_impl_.finalize(); }
  void diagonalize(rokko::distributed_crs_matrix& mat,
		   int num_evals, int block_size, int max_iters, double tol, timer& timer) {
    solver_impl_.diagonalize(mat, num_evals, block_size, max_iters, tol, timer);
  }
  void diagonalize(rokko::distributed_crs_matrix& mat,
		   int num_evals, int block_size, int max_iters, double tol) {
    timer_dumb timer;
    solver_impl_.diagonalize(mat, num_evals, block_size, max_iters, tol, timer);
  }
  void diagonalize(rokko::distributed_mfree& mat,
		   int num_evals, int block_size, int max_iters, double tol, timer& timer) {
    solver_impl_.diagonalize(mat, num_evals, block_size, max_iters, tol, timer);
  }
  void diagonalize(rokko::distributed_mfree& mat,
		   int num_evals, int block_size, int max_iters, double tol) {
    timer_dumb timer;
    solver_impl_.diagonalize(mat, num_evals, block_size, max_iters, tol, timer);
  }
  double eigenvalue(int i) const {
    return solver_impl_.eigenvalue(i);
  }

  void eigenvector(int i, std::vector<double>& vec) const {
    solver_impl_.eigenvector(i, vec);
  }

  int num_conv() const {
    return solver_impl_.num_conv();
  }

  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim, int col_dim) {
    return solver_impl_.create_distributed_crs_matrix(row_dim, col_dim);
  }

private:
  solver_type solver_impl_;
};

typedef factory<ps_solver_base> ps_solver_factory;
  
} // end namespace detail
  
class parallel_sparse_solver {
public:
  parallel_sparse_solver(std::string const& solver_name) {
    solver_impl_ = detail::ps_solver_factory::instance()->make_product(solver_name);
  }
  parallel_sparse_solver() {
    solver_impl_ = detail::ps_solver_factory::instance()->make_product();
  }
  void initialize(int& argc, char**& argv) { solver_impl_->initialize(argc, argv); }
  void finalize() { solver_impl_->finalize(); }

  void diagonalize(rokko::distributed_crs_matrix& mat,
		   int num_evals, int block_size, int max_iters, double tol, timer& timer) {
    solver_impl_->diagonalize(mat, num_evals, block_size, max_iters, tol, timer);
  }

  void diagonalize(rokko::distributed_crs_matrix& mat,
		   int num_evals, int block_size, int max_iters, double tol) {
    solver_impl_->diagonalize(mat, num_evals, block_size, max_iters, tol);
  }

  void diagonalize(rokko::distributed_mfree& mat,
		   int num_evals, int block_size, int max_iters, double tol, timer& timer) {
    solver_impl_->diagonalize(mat, num_evals, block_size, max_iters, tol, timer);
  }

  void diagonalize(rokko::distributed_mfree& mat,
		   int num_evals, int block_size, int max_iters, double tol) {
    solver_impl_->diagonalize(mat, num_evals, block_size, max_iters, tol);
  }

  double eigenvalue(int i) const {
    return solver_impl_->eigenvalue(i);
  }

  void eigenvector(int i, std::vector<double>& vec) const {
    solver_impl_->eigenvector(i, vec);
  }

  int num_conv() const {
    return solver_impl_->num_conv();
  }

  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim, int col_dim) {
    return solver_impl_->create_distributed_crs_matrix(row_dim, col_dim);
  }

  static std::vector<std::string> solvers() {
    return detail::ps_solver_factory::product_names();
  }
  static std::string default_solver() {
    return detail::ps_solver_factory::default_product_name();
  }
private:
  detail::ps_solver_factory::product_pointer_type solver_impl_;
};

} // end namespace rokko

#define ROKKO_REGISTER_PARALLEL_SPARSE_SOLVER(solver, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  typedef rokko::factory<rokko::detail::ps_solver_base> factory; \
  typedef rokko::detail::ps_solver_wrapper<solver> product; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

#endif // ROKKO_PARALLEL_SPARSE_SOLVER_HPP
