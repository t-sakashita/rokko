/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef PYROKKO_SERIAL_DENSE_EV_HPP
#define PYROKKO_SERIAL_DENSE_EV_HPP

#include <rokko/pyrokko_localized_matrix.hpp>
#include <rokko/pyrokko_parameters.hpp>
#include <rokko/serial_dense_ev.hpp>

namespace rokko {

class wrap_serial_dense_ev : public serial_dense_ev {
public:
  wrap_serial_dense_ev(std::string const& solver_name) : serial_dense_ev(solver_name) {}
  
  wrap_serial_dense_ev() {}

  void initialize() {
    int num = 1;
    char** ptr = NULL;
    serial_dense_ev::initialize(num, ptr);
  }
  
  template<typename VEC>
  wrap_parameters diagonalize(wrap_localized_matrix& mat, VEC& eigvals, wrap_localized_matrix& eigvecs,
			 wrap_parameters const& params) {
    assert(mat.is_major_col() == eigvecs.is_major_col());
    if (mat.is_major_col())
      return serial_dense_ev::diagonalize(mat.col_ver(), eigvals, eigvecs.col_ver(), parameters(params));
    else
      return serial_dense_ev::diagonalize(mat.row_ver(), eigvals, eigvecs.row_ver(), parameters(params));
  }

  template<typename VEC>
  wrap_parameters diagonalize(wrap_localized_matrix& mat, VEC& eigvals,
			 wrap_parameters const& params) {
    if (mat.is_major_col())
      return serial_dense_ev::diagonalize(mat.col_ver(), eigvals, parameters(params));
    else
      return serial_dense_ev::diagonalize(mat.row_ver(), eigvals, parameters(params));
  }

};

} // end namespace rokko

#endif // PYROKKO_SERIAL_DENSE_EV_HPP
