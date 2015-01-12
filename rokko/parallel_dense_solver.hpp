/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_PARALLEL_DENSE_SOLVER_HPP
#define ROKKO_PARALLEL_DENSE_SOLVER_HPP

#include <rokko/factory.hpp>
#include <rokko/mapping_bc.hpp>
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
                           distributed_matrix<matrix_row_major>& eigvecs, timer& timer) = 0;
  virtual void diagonalize(distributed_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                           distributed_matrix<matrix_col_major>& eigvecs, timer& timer) = 0;
  virtual void optimized_matrix_size(distributed_matrix<matrix_row_major>& mat) = 0;
  virtual void optimized_matrix_size(distributed_matrix<matrix_col_major>& mat) = 0;
  virtual mapping_bc optimized_mapping() = 0;
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
  void initialize(int& argc, char**& argv) {
    solver_impl_.initialize(argc, argv);
  }
  void finalize() { solver_impl_.finalize(); }
  void diagonalize(distributed_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                   distributed_matrix<matrix_row_major>& eigvecs, timer& timer) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer);
  }
  void diagonalize(distributed_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                   distributed_matrix<matrix_col_major>& eigvecs, timer& timer) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer);
  }
  void optimized_matrix_size(distributed_matrix<matrix_row_major>& mapping_bc) {
    solver_impl_.optimized_matrix_size(map);
  }
  void optimized_matrix_size(distributed_matrix<matrix_col_major>& mapping_bc) {
    solver_impl_.optimized_matrix_size(mapping_bc);
  }
  mapping_bc optimized_matrix_size() {
    return solver_impl_.optimized_matrix_size();
  }
private:
  solver_type solver_impl_;
};

typedef factory<pd_solver_base> pd_solver_factory;
  
} // end namespace detail
  
class parallel_dense_solver {
public:
  void construct(std::string const& solver_name, timer& timer) {
    if (!timer.has(rokko::timer_id::solver_construct))
      timer.registrate(rokko::timer_id::solver_construct, "solver::construct");
    timer.start(rokko::timer_id::solver_construct);
    solver_impl_ = detail::pd_solver_factory::instance()->make_product(solver_name);
    timer.stop(rokko::timer_id::solver_construct);
  }
  parallel_dense_solver(std::string const& solver_name, timer& timer) {
    this->construct(solver_name, timer);
  }
  parallel_dense_solver(std::string const& solver_name) {
    this->construct(solver_name, *global_timer::instance());
  }
  parallel_dense_solver(timer& timer) {
    this->construct(this->default_solver(), timer);
  }
  parallel_dense_solver() {
    this->construct(this->default_solver(), *global_timer::instance());
  }
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) {
    return solver_impl_->is_available_grid_major(grid_major);
  }
  void initialize(int& argc, char**& argv, timer& timer) {
    if (!timer.has(rokko::timer_id::solver_initialize))
      timer.registrate(rokko::timer_id::solver_initialize, "solver::initialize");
    timer.start(rokko::timer_id::solver_initialize);
    solver_impl_->initialize(argc, argv);
    timer.stop(rokko::timer_id::solver_initialize);
  }
  void initialize(int& argc, char**& argv) {
    this->initialize(argc, argv, *global_timer::instance());
  }
  void finalize(timer& timer) {
    if (!timer.has(rokko::timer_id::solver_finalize))
      timer.registrate(rokko::timer_id::solver_finalize, "solver::finalize");
    timer.start(rokko::timer_id::solver_finalize);
    solver_impl_->finalize();
    timer.stop(rokko::timer_id::solver_finalize);
  }
  void finalize() { this->finalize(*global_timer::instance()); }
  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, VEC& eigvals,
    rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs,
    timer& timer = *global_timer::instance()) {
    if (!timer.has(rokko::timer_id::diagonalize_initialize))
      timer.registrate(rokko::timer_id::diagonalize_initialize, "diagonalize::initialize");
    if (!timer.has(rokko::timer_id::diagonalize_diagonalize))
      timer.registrate(rokko::timer_id::diagonalize_diagonalize, "diagonalize::diagonalize");
    if (!timer.has(rokko::timer_id::diagonalize_finalize))
      timer.registrate(rokko::timer_id::diagonalize_finalize, "diagonalize::finalize");
    solver_impl_->diagonalize(mat, eigvals, eigvecs, timer);
  }
  // template <typename MATRIX_MAJOR>
  // void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) const {
  //   solver_impl_->optimized_matrix_size(mat);
  // }
  template <typename MATRIX_MAJOR>
  mapping_bc optimized_mapping() const {
    return solver_impl_->optimized_matrix_size();
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
