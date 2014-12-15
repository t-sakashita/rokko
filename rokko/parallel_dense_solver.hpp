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

#ifndef ROKKO_PARALLEL_DENSE_SOLVER_HPP
#define ROKKO_PARALLEL_DENSE_SOLVER_HPP

#include <rokko/factory.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {

namespace detail {
    
class pd_solver_base {
public:
  virtual ~pd_solver_base() {}
  virtual bool is_available_grid_major(grid_row_major_t const& grid_major) = 0;
  virtual bool is_available_grid_major(grid_col_major_t const& grid_major) = 0;
  virtual void initialize(int& argc, char**& argv) = 0;
  virtual void finalize() = 0;
  virtual void diagonalize(distributed_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                           distributed_matrix<matrix_row_major>& eigvecs) = 0;
  virtual void diagonalize(distributed_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                           distributed_matrix<matrix_row_major>& eigvecs, timer& timer_in) = 0;
  virtual void diagonalize(distributed_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                           distributed_matrix<matrix_col_major>& eigvecs) = 0;
  virtual void diagonalize(distributed_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                           distributed_matrix<matrix_col_major>& eigvecs, timer& timer_in) = 0;
  virtual void optimized_matrix_size(distributed_matrix<matrix_row_major>& mat) = 0;
  virtual void optimized_matrix_size(distributed_matrix<matrix_col_major>& mat) = 0;
};
  
template<typename SOLVER>
class pd_solver_wrapper : public pd_solver_base {
  typedef SOLVER solver_type;
public:
  pd_solver_wrapper() : solver_impl_() {}
  virtual ~pd_solver_wrapper() {}
  bool is_available_grid_major(grid_row_major_t const& grid_major) {
    return solver_impl_.is_available_grid_major(grid_major);
  }
  bool is_available_grid_major(grid_col_major_t const& grid_major) {
    return solver_impl_.is_available_grid_major(grid_major);
  }
  void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
  void finalize() { solver_impl_.finalize(); }
  void diagonalize(distributed_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                   distributed_matrix<matrix_row_major>& eigvecs) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, *global_timer::instance());
  }
  void diagonalize(distributed_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                   distributed_matrix<matrix_row_major>& eigvecs, timer& timer_in) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  void diagonalize(distributed_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                   distributed_matrix<matrix_col_major>& eigvecs) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, *global_timer::instance());
  }
  void diagonalize(distributed_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                   distributed_matrix<matrix_col_major>& eigvecs, timer& timer_in) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  void optimized_matrix_size(distributed_matrix<matrix_row_major>& mat) {
    solver_impl_.optimized_matrix_size(mat);
  }
  void optimized_matrix_size(distributed_matrix<matrix_col_major>& mat) {
    solver_impl_.optimized_matrix_size(mat);
  }
private:
  solver_type solver_impl_;
};

typedef factory<pd_solver_base> pd_solver_factory;
  
} // end namespace detail
  
class parallel_dense_solver {
public:
  parallel_dense_solver(std::string const& solver_name) {
    solver_impl_ = detail::pd_solver_factory::instance()->make_product(solver_name);
  }
  parallel_dense_solver() {
    solver_impl_ = detail::pd_solver_factory::instance()->make_product();
  }
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) {
    return solver_impl_->is_available_grid_major(grid_major);
  }
  void initialize(int& argc, char**& argv) { solver_impl_->initialize(argc, argv); }
  void finalize() { solver_impl_->finalize(); }
  template<typename MATRIX_MAJOR>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
    rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
    solver_impl_->diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  template<typename MATRIX_MAJOR>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
    rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs) {
    solver_impl_->diagonalize(mat, eigvals, eigvecs, *global_timer::instance());
  }
  template <typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) const {
    solver_impl_->optimized_matrix_size(mat);
  }
  static std::vector<std::string> solvers() {
    return detail::pd_solver_factory::product_names();
  }
  static std::string default_solver() {
    return detail::pd_solver_factory::default_product_name();
  }
private:
  detail::pd_solver_factory::product_pointer_type solver_impl_;
};

} // end namespace rokko


#define ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(solver, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  typedef rokko::factory<rokko::detail::pd_solver_base> factory; \
  typedef rokko::detail::pd_solver_wrapper<solver> product; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

#endif // ROKKO_PARALLEL_DENSE_SOLVER_HPP
