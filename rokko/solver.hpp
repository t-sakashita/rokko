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
#include <boost/shared_ptr.hpp>

namespace rokko {

class solver {
public:
  solver(std::string const& solver_name) {
    solver_impl_ = solver_factory::instance()->make_solver(solver_name);
  }
  void initialize(int& argc, char**& argv) { solver_impl_->initialize(argc, argv); }
  void finalize() { solver_impl_->finalize(); }
  template<typename MATRIX_MAJOR>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
    rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs) {
    solver_impl_->diagonalize(mat, eigvals, eigvecs);
  }
  //void optimized_matrix_size(int dim, int nprow, int npcol, int& mb, int& nb, int& lld, int& len_array) {
  //  solver_impl_->optimized_matrix_size(dim, nprow, npcol, mb, nb, lld, len_array);
  //}
  template <typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) const {
    solver_impl_->optimized_matrix_size(mat);
  }
private:
  solver_factory::solver_pointer_type solver_impl_;
};

} // end namespace rokko

#endif // ROKKO_SOLVER_HPP
