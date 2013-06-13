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

#ifndef ROKKO_SOLVER_HPP
#define ROKKO_SOLVER_HPP

#include <rokko/solver_factory.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <boost/shared_ptr.hpp>

namespace rokko {

class solver {
public:
  solver(std::string const& solver_name) {
    solver_impl_ = solver_factory::instance()->make_solver(solver_name);
  }
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) { return solver_impl_->is_available_grid_major(grid_major);  }
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
    timer_dumb timer_in;
    solver_impl_->diagonalize(mat, eigvals, eigvecs, timer_in);
  }
  template <typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) const {
    solver_impl_->optimized_matrix_size(mat);
  }
private:
  solver_factory::solver_pointer_type solver_impl_;
};

} // end namespace rokko

#endif // ROKKO_SOLVER_HPP
