/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014  Rokko Developers https://github.com/t-sakashita/rokko
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
#include <rokko/distributed_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/parameters.hpp>

namespace rokko {

namespace detail {
    
class ps_ev_base {
public:
  virtual ~ps_ev_base() {}
  virtual void initialize(int& argc, char**& argv) = 0;
  virtual void finalize() = 0;
  virtual parameters diagonalize(rokko::distributed_crs_matrix& mat, rokko::parameters const& params) = 0;
  virtual parameters diagonalize(rokko::distributed_mfree& mat, rokko::parameters const& params) = 0;
  virtual double eigenvalue(int k) const = 0;
  virtual void eigenvector(int k, std::vector<double>& vec) const = 0;
  virtual void eigenvector(int k, double* vec) const = 0;
  virtual void eigenvector(int k, distributed_vector<double>& vec) const = 0;
  virtual int get_num_conv() const = 0;
  virtual rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
										    int col_dim) = 0;
  virtual rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
										  int col_dim, int num_entries_per_row) = 0;
  
};

template<typename SOLVER>
class ps_ev_wrapper : public ps_ev_base {
  using solver_type = SOLVER;
public:
  ps_ev_wrapper() : solver_impl_() {}
  virtual ~ps_ev_wrapper() {}
  void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
  void finalize() { solver_impl_.finalize(); }
  parameters diagonalize(rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    return solver_impl_.diagonalize(mat, params);
  }
  parameters diagonalize(rokko::distributed_mfree& mat, rokko::parameters const& params) {
    return solver_impl_.diagonalize(mat, params);
  }
  double eigenvalue(int k) const { return solver_impl_.eigenvalue(k); }
  void eigenvector(int k, std::vector<double>& vec) const { solver_impl_.eigenvector(k, vec); }
  void eigenvector(int k, double* vec) const { solver_impl_.eigenvector(k, vec); }
  void eigenvector(int k, distributed_vector<double>& vec) const { solver_impl_.eigenvector(k, vec); }
  int get_num_conv() const { return solver_impl_.get_num_conv(); }
  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim) {
    return solver_impl_.create_distributed_crs_matrix(row_dim, col_dim);
  }
  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim, int num_entries_per_row) {
    return solver_impl_.create_distributed_crs_matrix(row_dim, col_dim, num_entries_per_row);
  }
private:
  solver_type solver_impl_;
};

using ps_solver_factory = factory<ps_ev_base>;
  
} // end namespace detail
  
class parallel_sparse_ev {
public:
  void construct(std::string const& solver_name) {
    solver_name_ = solver_name;
    solver_impl_ = detail::ps_solver_factory::instance()->make_product(solver_name);
  }
  parallel_sparse_ev(std::string const& solver_name) {
    this->construct(solver_name);
  }
  parallel_sparse_ev() {
    this->construct(this->default_solver());
  }
  void initialize(int& argc, char**& argv) {
    solver_impl_->initialize(argc, argv);
  }
  void finalize() {
    solver_impl_->finalize();
  }
  template<typename MAT>
  parameters diagonalize(MAT& mat, rokko::parameters const& params) {
    return solver_impl_->diagonalize(mat, params);
  }
  double eigenvalue(int k) const { return solver_impl_->eigenvalue(k); }
  void eigenvector(int k, std::vector<double>& vec) const { solver_impl_->eigenvector(k, vec); }
  void eigenvector(int k, double* vec) const { solver_impl_->eigenvector(k, vec); }
  void eigenvector(int k, distributed_vector<double>& vec) const { solver_impl_->eigenvector(k, vec); }
  int get_num_conv() const { return solver_impl_->get_num_conv(); }
  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim) {
    return solver_impl_->create_distributed_crs_matrix(row_dim, col_dim);
  }
  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim, int num_entries_per_row) {
    return solver_impl_->create_distributed_crs_matrix(row_dim, col_dim, num_entries_per_row);
  }
  std::string get_solver_name() { return solver_name_; }
  static std::vector<std::string> solvers() { return detail::ps_solver_factory::product_names(); }
  static std::string default_solver() { return detail::ps_solver_factory::default_product_name(); }
private:
  detail::ps_solver_factory::product_pointer_type solver_impl_;
  std::string solver_name_;
};

} // end namespace rokko

#define ROKKO_REGISTER_PARALLEL_SPARSE_SOLVER(solver, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  using factory = rokko::factory<rokko::detail::ps_ev_base>;  \
  using product = rokko::detail::ps_ev_wrapper<solver>; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

#endif // ROKKO_PARALLEL_SPARSE_SOLVER_HPP
