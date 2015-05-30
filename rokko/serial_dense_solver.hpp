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

#ifndef ROKKO_SERIAL_DENSE_SOLVER_HPP
#define ROKKO_SERIAL_DENSE_SOLVER_HPP

#include <rokko/factory.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/parameters.hpp>

namespace rokko {

namespace detail {

class sd_solver_base {
public:
  virtual ~sd_solver_base() {}
  virtual void initialize(int& argc, char**& argv) = 0;
  virtual void finalize() = 0;
  virtual void diagonalize(localized_matrix<double, matrix_row_major>& mat,
                           localized_vector<double>& eigvals,
                           localized_matrix<double, matrix_row_major>& eigvecs,
			   rokko::parameters const& params, timer& timer) = 0;
  virtual void diagonalize(localized_matrix<double, matrix_col_major>& mat,
                           localized_vector<double>& eigvals,
                           localized_matrix<double, matrix_col_major>& eigvecs,
			   rokko::parameters const& params, timer& timer) = 0;
  virtual void diagonalize(localized_matrix<double, matrix_row_major>& mat, std::vector<double>& eigvals,
                           localized_matrix<double, matrix_row_major>& eigvecs,
			   rokko::parameters const& params, timer& timer) = 0;
  virtual void diagonalize(localized_matrix<double, matrix_col_major>& mat, std::vector<double>& eigvals,
                           localized_matrix<double, matrix_col_major>& eigvecs,
			   rokko::parameters const& params, timer& timer) = 0;
  // no parameters
  virtual void diagonalize(localized_matrix<double, matrix_row_major>& mat,
                           localized_vector<double>& eigvals,
                           localized_matrix<double, matrix_row_major>& eigvecs,
			   timer& timer) = 0;
  virtual void diagonalize(localized_matrix<double, matrix_col_major>& mat,
                           localized_vector<double>& eigvals,
                           localized_matrix<double, matrix_col_major>& eigvecs,
			   timer& timer) = 0;
  virtual void diagonalize(localized_matrix<double, matrix_row_major>& mat, std::vector<double>& eigvals,
                           localized_matrix<double, matrix_row_major>& eigvecs,
			   timer& timer) = 0;
  virtual void diagonalize(localized_matrix<double, matrix_col_major>& mat, std::vector<double>& eigvals,
                           localized_matrix<double, matrix_col_major>& eigvecs,
			   timer& timer) = 0;
};
  
template<typename SOLVER>
class sd_solver_wrapper : public sd_solver_base {
  typedef SOLVER solver_type;
public:
  sd_solver_wrapper() : solver_impl_() {}
  virtual ~sd_solver_wrapper() {}
  void initialize(int& argc, char**& argv) {
    solver_impl_.initialize(argc, argv);
  }
  void finalize() { solver_impl_.finalize(); }
  void diagonalize(localized_matrix<double, matrix_row_major>& mat, localized_vector<double>& eigvals,
                   localized_matrix<double, matrix_row_major>& eigvecs,
		   rokko::parameters const& params, timer& timer) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, params, timer);
  }
  void diagonalize(localized_matrix<double, matrix_col_major>& mat, localized_vector<double>& eigvals,
                   localized_matrix<double, matrix_col_major>& eigvecs,
		   rokko::parameters const& params, timer& timer) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, paramas, timer);
  }
  void diagonalize(localized_matrix<double, matrix_row_major>& mat, std::vector<double>& eigvals,
                   localized_matrix<double, matrix_row_major>& eigvecs,
		   rokko::parameters const& params, timer& timer) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, params, timer);
  }
  void diagonalize(localized_matrix<double, matrix_col_major>& mat, std::vector<double>& eigvals,
                   localized_matrix<double, matrix_col_major>& eigvecs,
		   rokko::parameters const& params, timer& timer) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, params, timer);
  }
private:
  solver_type solver_impl_;
};
    
typedef factory<sd_solver_base> sd_solver_factory;
  
} // end namespace detail
  
class serial_dense_solver {
public:
  void construct(std::string const& solver_name, timer& timer) {
    if (!timer.has(rokko::timer_id::solver_construct))
      timer.registrate(rokko::timer_id::solver_construct, "solver::construct");
    timer.start(rokko::timer_id::solver_construct);
    solver_impl_ = detail::sd_solver_factory::instance()->make_product(solver_name);
    timer.stop(rokko::timer_id::solver_construct);
  }
  serial_dense_solver(std::string const& solver_name, timer& timer) {
    this->construct(solver_name, timer);
  }
  serial_dense_solver(std::string const& solver_name) {
    this->construct(solver_name, *global_timer::instance());
  }
  serial_dense_solver(timer& timer) {
    this->construct(this->default_solver(), timer);
  }
  serial_dense_solver() {
    this->construct(this->default_solver(), *global_timer::instance());
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
  void diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat, VEC& eigvals,
                   localized_matrix<double, MATRIX_MAJOR>& eigvecs,
                   rokko::parameters const& params = null_parameters, timer& timer = *global_timer::instance()) {
    if (!timer.has(rokko::timer_id::diagonalize_initialize))
      timer.registrate(rokko::timer_id::diagonalize_initialize, "diagonalize::initialize");
    if (!timer.has(rokko::timer_id::diagonalize_diagonalize))
      timer.registrate(rokko::timer_id::diagonalize_diagonalize, "diagonalize::diagonalize");
    if (!timer.has(rokko::timer_id::diagonalize_finalize))
      timer.registrate(rokko::timer_id::diagonalize_finalize, "diagonalize::finalize");
    solver_impl_->diagonalize(mat, eigvals, eigvecs, timer);
  }
  static std::vector<std::string> solvers() {
    return detail::sd_solver_factory::product_names();
  }
  static std::string default_solver() {
    return detail::sd_solver_factory::default_product_name();
  }
private:
  detail::sd_solver_factory::product_pointer_type solver_impl_;
};

} // end namespace rokko

#define ROKKO_REGISTER_SERIAL_DENSE_SOLVER(solver, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  typedef rokko::factory<rokko::detail::sd_solver_base> factory; \
  typedef rokko::detail::sd_solver_wrapper<solver> product; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

#endif // ROKKO_SERIAL_DENSE_SOLVER_HPP
