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

#ifndef ROKKO_CHECK_ORTHOGONALITY_HPP
#define ROKKO_CHECK_ORTHOGONALITY_HPP

#include <rokko/distributed_matrix.hpp>

namespace rokko {

template<typename T, typename MATRIX_MAJOR>
int check_orthogonality(rokko::distributed_matrix<T, MATRIX_MAJOR>& eigvecs) {

  rokko::distributed_matrix<T, MATRIX_MAJOR> mat2(mat.get_m_global(), mat.get_n_global(), mat.g);
  product(mat, false, mat, true, 1, 0, mat2);
  std::cout << "check_orthogonality";
  mat2.print();

  return 0;
}

} // namespace rokko

#endif // ROKKO_CHECK_ORTHOGONALITY_HPP
