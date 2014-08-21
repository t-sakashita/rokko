/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
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

namespace rokko {

namespace detail {

class sd_solver_base {
public:
  virtual ~sd_solver_base() {}
  virtual void initialize(int& argc, char**& argv) = 0;
  virtual void finalize() = 0;
  virtual void diagonalize(localized_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                           localized_matrix<matrix_row_major>& eigvecs) = 0;
  virtual void diagonalize(localized_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                           localized_matrix<matrix_row_major>& eigvecs, timer& timer_in) = 0;
  virtual void diagonalize(localized_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                           localized_matrix<matrix_col_major>& eigvecs) = 0;
  virtual void diagonalize(localized_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                           localized_matrix<matrix_col_major>& eigvecs, timer& timer_in) = 0;

  virtual void diagonalize(localized_matrix<matrix_row_major>& mat, std::vector<double>& eigvals,
                           localized_matrix<matrix_row_major>& eigvecs) = 0;
  virtual void diagonalize(localized_matrix<matrix_row_major>& mat, std::vector<double>& eigvals,
                           localized_matrix<matrix_row_major>& eigvecs, timer& timer_in) = 0;
  virtual void diagonalize(localized_matrix<matrix_col_major>& mat, std::vector<double>& eigvals,
                           localized_matrix<matrix_col_major>& eigvecs) = 0;
  virtual void diagonalize(localized_matrix<matrix_col_major>& mat, std::vector<double>& eigvals,
                           localized_matrix<matrix_col_major>& eigvecs, timer& timer_in) = 0;
};
  
template<typename SOLVER>
class sd_solver_wrapper : public sd_solver_base {
  typedef SOLVER solver_type;
public:
  sd_solver_wrapper() : solver_impl_() {}
  virtual ~sd_solver_wrapper() {}
  void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
  void finalize() { solver_impl_.finalize(); }
  void diagonalize(localized_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                   localized_matrix<matrix_row_major>& eigvecs) {
    timer_dumb timer;
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer);
  }
  void diagonalize(localized_matrix<matrix_row_major>& mat, localized_vector& eigvals,
                   localized_matrix<matrix_row_major>& eigvecs, timer& timer_in) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  void diagonalize(localized_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                   localized_matrix<matrix_col_major>& eigvecs) {
    timer_dumb timer;
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer);
  }
  void diagonalize(localized_matrix<matrix_col_major>& mat, localized_vector& eigvals,
                   localized_matrix<matrix_col_major>& eigvecs, timer& timer_in) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
  }

  void diagonalize(localized_matrix<matrix_row_major>& mat, std::vector<double>& eigvals,
                   localized_matrix<matrix_row_major>& eigvecs) {
    timer_dumb timer;
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer);
  }
  void diagonalize(localized_matrix<matrix_row_major>& mat, std::vector<double>& eigvals,
                   localized_matrix<matrix_row_major>& eigvecs, timer& timer_in) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  void diagonalize(localized_matrix<matrix_col_major>& mat, std::vector<double>& eigvals,
                   localized_matrix<matrix_col_major>& eigvecs) {
    timer_dumb timer;
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer);
  }
  void diagonalize(localized_matrix<matrix_col_major>& mat, std::vector<double>& eigvals,
                   localized_matrix<matrix_col_major>& eigvecs, timer& timer_in) {
    solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
  }
private:
  solver_type solver_impl_;
};
    
typedef factory<sd_solver_base> sd_solver_factory;
  
} // end namespace detail
  
class serial_dense_solver {
public:
  serial_dense_solver(std::string const& solver_name) {
    solver_impl_ = detail::sd_solver_factory::instance()->make_product(solver_name);
  }
  void initialize(int& argc, char**& argv) { solver_impl_->initialize(argc, argv); }
  void finalize() { solver_impl_->finalize(); }
  template<typename MATRIX_MAJOR>
  void diagonalize(localized_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                   localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
    solver_impl_->diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  template<typename MATRIX_MAJOR>
  void diagonalize(localized_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                   localized_matrix<MATRIX_MAJOR>& eigvecs) {
    timer_dumb timer_in;
    solver_impl_->diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  template<typename MATRIX_MAJOR>
  void diagonalize(localized_matrix<MATRIX_MAJOR>& mat, std::vector<double>& eigvals,
                   localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
    solver_impl_->diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  template<typename MATRIX_MAJOR>
  void diagonalize(localized_matrix<MATRIX_MAJOR>& mat, std::vector<double>& eigvals,
                   localized_matrix<MATRIX_MAJOR>& eigvecs) {
    timer_dumb timer_in;
    solver_impl_->diagonalize(mat, eigvals, eigvecs, timer_in);
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

#endif // ROKKO_SERIAL_DENSE_SOLVER_HPP

#define ROKKO_REGISTER_SERIAL_DENSE_SOLVER(solver, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  typedef rokko::factory<rokko::detail::sd_solver_base> factory; \
  typedef rokko::detail::sd_solver_wrapper<solver> product; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }
