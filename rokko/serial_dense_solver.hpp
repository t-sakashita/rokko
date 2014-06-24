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

#include <rokko/solver_factory.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <boost/shared_ptr.hpp>

namespace rokko {

class serial_dense_solver {
public:
  serial_dense_solver(std::string const& solver_name) {
    solver_impl_ = solver_factory::instance()->make_serial_dense_solver(solver_name);
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
  static std::vector<std::string> solvers() {
    return solver_factory::serial_dense_solver_names();
  }
  static std::string default_solver() {
    return solver_factory::default_serial_dense_solver_name();
  }
private:
  solver_factory::serial_dense_solver_pointer_type solver_impl_;
};

} // end namespace rokko

#endif // ROKKO_SERIAL_DENSE_SOLVER_HPP
